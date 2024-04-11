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
from google.cloud import bigquery, storage
from pandas import json_normalize

from gcp_utils import (  # export_dataframe_to_gcs_as_json,
    convert_types,
    create_df_from_json_for_index_file,
    upload_dataframe_to_bq,
)

# Create a handler for Google Cloud Logging.
logging.basicConfig(level=logging.INFO)

with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

# Accessing GCP configuration
project_id = config["GCP"]["PROJECT_ID"]
bucket_name = config["GCP"]["BUCKET_NAME"]
dataset_name = config["GCP"]["CLOUDASM_STANDARD_REGIONS_DATASET"]
max_p_value = config["ASM"]["MAX_P_VALUE"]
max_bh_threshold = config["ASM"]["MAX_BH_THRESHOLD"]
min_nb_cpg_same_direction = config["ASM"]["MIN_NB_CPG_SAME_DIRECTION"]
min_nb_consecutive_cpg_same_direction = config["ASM"][
    "MIN_NB_CONSECUTIVE_CPG_SAME_DIRECTION"
]
min_asm_region_effect = config["ASM"]["MIN_ASM_REGION_EFFECT"]

# Retrieve Job-defined env vars
BATCH_TASK_INDEX = int(os.getenv("BATCH_TASK_INDEX", 0))

# Initialize the Google Cloud Storage client
bq_client = bigquery.Client(project=project_id)
storage_client = storage.Client(project=project_id)


def consecutive_cpg(row, p_value):
    if int(row["nb_sig_cpg"]) > 1:
        flat_cpg = json_normalize(row["cpg"])
        # logging.info(f"Flat CpG: {flat_cpg}")
        max_nb_consec = {"positive": 0, "negative": 0}
        current_nb_consec = {"positive": 0, "negative": 0}
        for index, row in flat_cpg.iterrows():
            if index > 0:
                prev_row = flat_cpg.iloc[index - 1]
                for direction, sign in [("positive", 1), ("negative", -1)]:
                    if (
                        prev_row.fisher_pvalue < p_value
                        and row.fisher_pvalue < p_value
                        and np.sign(prev_row.effect) == sign
                        and np.sign(row.effect) == sign
                    ):
                        if current_nb_consec[direction] == 0:
                            current_nb_consec[direction] = 2
                        else:
                            current_nb_consec[direction] += 1
                        max_nb_consec[direction] = max(
                            max_nb_consec[direction], current_nb_consec[direction]
                        )
                    else:
                        current_nb_consec[direction] = 0
        return max_nb_consec
    else:
        return {"positive": 0, "negative": 0}


def wilcoxon_pvalue(row):
    # logging.info(f"row ref {row['ref']}")
    # logging.info(f"row ref {row['alt']}")
    # logging.info(f"JSON nromalized row ref {json_normalize(row['ref'])}")
    try:
        _, pvalue = stats.mannwhitneyu(
            json_normalize(row["ref"]),
            json_normalize(row["alt"]),
            alternative="two-sided",
        )
        return np.round(pvalue[0], 5)
    # If the ref and alt datasets are equal or one is included in the other one:
    except ValueError:
        return 1


def find_asm(
    x,
    max_p_value,
    min_nb_cpg_same_direction,
    min_nb_consecutive_cpg_same_direction,
    min_asm_region_effect,
):
    if (
        x["wilcoxon_corr_p_value"] < max_p_value
        and abs(x["asm_region_effect"]) > min_asm_region_effect
    ):
        positive_condition = (
            x["pos_sig_cpg"] >= min_nb_cpg_same_direction
            and x["nb_consec_pos_sig_asm"] >= min_nb_consecutive_cpg_same_direction
        )
        negative_condition = (
            x["neg_sig_cpg"] >= min_nb_cpg_same_direction
            and x["nb_consec_neg_sig_asm"] >= min_nb_consecutive_cpg_same_direction
        )
        if positive_condition or negative_condition:
            return 1
    return 0


def main():

    logging.info("Importing the JSON file as a dataframe")
    # logging.info(f"bucket name: {bucket_name}")
    # logging.info(f"dataset name: {dataset_name}")
    # logging.info(f"batch task index: {BATCH_TASK_INDEX}")
    df, file_name = create_df_from_json_for_index_file(
        storage_client, bucket_name, dataset_name, BATCH_TASK_INDEX, 1
    )

    logging.info(f"File names: {file_name}")

    logging.info(f"Number of rows in DF: {len(df)}")

    for column in df.columns:
        if column not in [
            "sample",
            "asm_region_effect",
            "snp_id",
            "cpg",
            "ref",
            "alt",
            # "wilcoxon_pvalue",
            # "wilcoxon_corr_pvalue",
        ]:
            df[column] = df[column].astype(int)
    df["asm_region_effect"] = df["asm_region_effect"].astype(float)

    for col in ["sample", "snp_id"]:
        df[col] = df[col].astype(int)

    type_map = {
        "pos": int,
        "ref_cov": int,
        "ref_meth": int,
        "alt_cov": int,
        "alt_meth": int,
        "effect": float,
        "fisher_pvalue": float,
        "methyl_perc": float,
    }

    logging.info("Converting the nested arrays to the right data types")
    for col in ["cpg", "ref", "alt"]:
        df[col] = df[col].apply(lambda x: [convert_types(item, type_map) for item in x])

    logging.info("Computing Wilcoxon p-value")
    df["wilcoxon_pvalue"] = df.apply(wilcoxon_pvalue, axis=1)

    logging.info(
        "Calculate p-value corrected for multiple testing using Benjamini–Hochberg"
    )

    df["wilcoxon_corr_pvalue"] = np.round(
        mt.multipletests(
            df["wilcoxon_pvalue"], alpha=max_bh_threshold, method="fdr_bh"
        )[1],
        5,
    )

    logging.info(
        "Calculate number of significant consecutive CpGs in the same direction."
    )

    results = df.apply(lambda row: consecutive_cpg(row, max_p_value), axis=1)

    df["nb_consec_pos_sig_asm"] = results.apply(lambda x: x["positive"])
    df["nb_consec_neg_sig_asm"] = results.apply(lambda x: x["negative"])

    df["asm"] = df.apply(
        lambda x: find_asm(
            x,
            max_p_value,
            min_nb_cpg_same_direction,
            min_nb_consecutive_cpg_same_direction,
            min_asm_region_effect,
        ),
        axis=1,
    )

    logging.info("Saving table to BigQuery")
    upload_dataframe_to_bq(bq_client, df, f"{dataset_name}.asm_flagged")


# Start script
if __name__ == "__main__":
    try:
        main()
    except Exception as err:
        message = f"Task #{BATCH_TASK_INDEX} failed: {str(err)}"

        print(json.dumps({"message": message, "severity": "ERROR"}))
        sys.exit(1)  # Retry Job Task by exiting the process
