#!/usr/bin/env python

# import ndjson
import json
import logging
import os
import sys

import numpy as np
import pandas as pd
import scipy.stats as stats

# from pandas.io.json import json_normalize
import statsmodels.stats.multitest as mt
import yaml
from gcp import upload_dataframe_to_bq
from google.cloud import bigquery

# Create a handler for Google Cloud Logging.
logging.basicConfig(level=logging.INFO)

with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

# Accessing GCP configuration
project = config["GCP"]["PROJECT"]
bucket_name = config["GCP"]["BUCKET"]
samples_dic = config["SAMPLES"]
all_samples = [item for sublist in samples_dic.values() for item in sublist]

# ASM variables
max_p_value = config["ASM"]["MAX_P_VALUE"]
max_bh_threshold = config["ASM"]["MAX_BH_THRESHOLD"]
min_nb_cpg_same_direction = config["ASM"]["MIN_NB_CPG_SAME_DIRECTION"]
min_nb_consecutive_cpg_same_direction = config["ASM"][
    "MIN_NB_CONSECUTIVE_CPG_SAME_DIRECTION"
]
min_asm_region_effect = config["ASM"]["MIN_ASM_REGION_EFFECT"]

# Retrieve Job-defined env vars
BATCH_TASK_INDEX = int(os.getenv("BATCH_TASK_INDEX", 0))
samples_dataset = os.getenv("SAMPLES_DATASET")

# Retrieve sample
sample = all_samples[BATCH_TASK_INDEX]

# Define the path to the JSON credentials file
credentials_path = "/appuser/.config/gcloud/application_default_credentials.json"

# Check if the JSON credentials file exists
if os.path.exists(credentials_path):
    # Set the GOOGLE_APPLICATION_CREDENTIALS environment variable to the path of the JSON key file
    os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = credentials_path
    # Assuming 'project' is already defined somewhere in your script
    os.environ["GOOGLE_CLOUD_PROJECT"] = project
    samples_dataset = "samples_250bp"

# Initialize the Google Cloud Storage client
bq_client = bigquery.Client(project=project)


def compute_asm_effect_and_wilcoxon_pvalue(row):
    # Obtain the fractional methylation of reads based on their tag ref or alt
    ref_values = [item["ref"] for item in row if item["ref"] is not None]
    alt_values = [item["alt"] for item in row if item["alt"] is not None]
    # Ensure there are enough samples to compute the test
    if len(ref_values) == 0 or len(alt_values) == 0:
        raise ValueError(
            "Ref and/or Alt values are empty after filtering None. Cannot perform statistical test."
        )
    try:
        result = stats.mannwhitneyu(ref_values, alt_values, alternative="two-sided")
        p_value = np.round(result.pvalue, 5)
    # If the ref and alt datasets are equal or one is included in the other one:
    except ValueError:
        p_value = 1
    read_asm_effect = np.round(
        sum(alt_values) / len(alt_values) - sum(ref_values) / len(ref_values), 5
    )
    return {"read_asm_effect": read_asm_effect, "wilcoxon_pvalue": p_value}


def find_max_consecutive_positions(qualifying_cpg_pos, all_cpg_pos):
    qualifying_lookup = {value: i for i, value in enumerate(qualifying_cpg_pos)}
    # To track the maximum length of consecutive elements
    max_consecutive = 0
    current_streak = 0
    # Iterate over all_cpg_pos to count consecutive qualifying elements
    for pos in all_cpg_pos:
        if pos in qualifying_lookup:
            if current_streak == 0:  # Starting a new consecutive sequence
                current_streak = 1
            else:
                # Check if this is truly consecutive in terms of qualifying_cpg_pos
                last_pos = all_cpg_pos.index(pos) - 1
                if last_pos >= 0 and all_cpg_pos[last_pos] in qualifying_lookup:
                    current_streak += 1
                else:
                    current_streak = 1
        else:
            if current_streak > max_consecutive:
                max_consecutive = current_streak
            current_streak = 0
    # Final check to update max_consecutive at the end of the loop
    if current_streak > max_consecutive:
        max_consecutive = current_streak
    return max_consecutive


def count_elements_and_consecutive(cpg_data, pvalue_threshold, read_asm_effect):
    """
    Count the number of elements that have a p-value below a given threshold and
    an effect sign that matches the sign of a specified effect value. Also counts consecutive
    elements meeting these criteria based on position continuity.

    Parameters:
    - cpg_data (list of dicts or numpy structured array): Array of records with 'pos', 'fisher_pvalue', and 'effect'.
    - pvalue_threshold (float): The threshold below which the Fisher p-value must fall.
    - read_asm_effect (float): The specified effect value, whose sign will determine the sign of the effect to filter by.

    Returns:
    - dict: A dictionary with 'total_count' of elements matching criteria and 'consecutive_count' of consecutive elements.
    """
    # Determine the sign of the effect from read_asm_effect
    effect_sign = "positive" if read_asm_effect > 0 else "negative"
    # Filter based on p-value and effect sign
    qualifying_cpg_pos = [
        d["pos"]
        for d in cpg_data
        if d["fisher_pvalue"] < pvalue_threshold
        and (
            (effect_sign == "positive" and d["effect"] > 0)
            or (effect_sign == "negative" and d["effect"] < 0)
        )
    ]
    all_cpg_pos = [d["pos"] for d in cpg_data]
    total_sig_cpgs = len(qualifying_cpg_pos)
    consecutive_sig_cpgs = find_max_consecutive_positions(
        qualifying_cpg_pos, all_cpg_pos
    )
    return {
        "total_sig_cpgs": total_sig_cpgs,
        "consecutive_sig_cpgs": consecutive_sig_cpgs,
    }


def find_asm(
    x,
    wilcoxon_to_use,
    max_p_value,
    min_nb_cpg_same_direction,
    min_nb_consecutive_cpg_same_direction,
    min_asm_region_effect,
):
    if (
        x[wilcoxon_to_use] <= max_p_value
        and abs(x["read_asm_effect"]) >= min_asm_region_effect
        and x["consecutive_sig_cpgs"] >= min_nb_consecutive_cpg_same_direction
        and x["total_sig_cpgs"] >= min_nb_cpg_same_direction
    ):
        return 1
    return 0


def main():

    logging.info("Importing the sample data necessary to compute ASM")

    query = f"SELECT sample, chr, region_inf, region_sup, clustering_index, snp_id, snp_pos, cpg_w_snp, reads_w_snp \
                FROM {project}.{samples_dataset}.regions_and_snps \
                WHERE sample = '{sample}' \
                "

    # Execute Query and store as DF
    df = bq_client.query(query).to_dataframe()

    logging.info(f"Number of rows for sample {sample}: {len(df)}")

    logging.info("Computing Wilcoxon p-value")
    # Apply the function and get a series of tuples
    df[["read_asm_effect", "wilcoxon_pvalue"]] = df["reads_w_snp"].apply(
        lambda x: pd.Series(compute_asm_effect_and_wilcoxon_pvalue(x))
    )

    logging.info(
        "Calculate p-value corrected for multiple testing using Benjamini–Hochberg"
    )

    df["corrected_wilcoxon_pvalue"] = np.round(
        mt.multipletests(
            df["wilcoxon_pvalue"], alpha=max_bh_threshold, method="fdr_bh"
        )[1],
        5,
    )

    logging.info(
        "Calculate the number of consecutive and total CpGs that are significative in the same direction of the read ASM"
    )
    df[["total_sig_cpgs", "consecutive_sig_cpgs"]] = df.apply(
        lambda x: pd.Series(
            count_elements_and_consecutive(
                x["cpg_w_snp"], max_p_value, x["read_asm_effect"]
            )
        ),
        axis=1,
    )

    logging.info("Use all the parameters to find out if ASM is present")
    df["asm"] = df.apply(
        lambda x: find_asm(
            x,
            "corrected_wilcoxon_pvalue",
            max_p_value,
            min_nb_cpg_same_direction,
            min_nb_consecutive_cpg_same_direction,
            min_asm_region_effect,
        ),
        axis=1,
    )

    df["asm_not_corrected"] = df.apply(
        lambda x: find_asm(
            x,
            "wilcoxon_pvalue",
            max_p_value,
            min_nb_cpg_same_direction,
            min_nb_consecutive_cpg_same_direction,
            min_asm_region_effect,
        ),
        axis=1,
    )

    upload_dataframe_to_bq(
        bq_client,
        df[
            [
                "sample",
                "chr",
                "region_inf",
                "region_sup",
                "clustering_index",
                "snp_id",
                "snp_pos",
                "read_asm_effect",
                "wilcoxon_pvalue",
                "corrected_wilcoxon_pvalue",
                "total_sig_cpgs",
                "consecutive_sig_cpgs",
                "asm",
                "asm_not_corrected",
            ]
        ],
        f"{samples_dataset}.asm_flagged",
    )

    logging.info("END OF SCRIPT")


# Start script
if __name__ == "__main__":
    try:
        main()
    except Exception as err:
        message = f"Task #{BATCH_TASK_INDEX} failed: {str(err)}"

        print(json.dumps({"message": message, "severity": "ERROR"}))
        sys.exit(1)  # Retry Job Task by exiting the process
