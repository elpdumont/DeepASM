#!/usr/bin/env python

# import ndjson
import json
import logging
import os
import sys

import numpy as np
import scipy.stats as stats

# from pandas.io.json import json_normalize
import statsmodels.stats.multitest as mt
import yaml
from google.cloud import bigquery, storage
from pandas import json_normalize

from gcp_utils import (  # export_dataframe_to_gcs_as_json,
    create_df_from_json_for_index_file,
    upload_dataframe_to_bq,
)

# Initialize the Google Cloud Storage client
bq_client = bigquery.Client()
storage_client = storage.Client()

# Create a handler for Google Cloud Logging.
logging.basicConfig(level=logging.INFO)

with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

# Accessing GCP configuration
bucket_name = config["GCP"]["BUCKET_NAME"]
dataset_name = config["GCP"]["CLOUDASM_STANDARD_REGIONS_DATASET"]
p_value = config["ASM"]["P_VALUE"]
bh_threshold = config["ASM"]["BH_THRESHOLD"]

# Retrieve Job-defined env vars
BATCH_TASK_INDEX = int(os.getenv("BATCH_TASK_INDEX", 0))
# nb_files_per_task = int(os.getenv("NB_FILES_PER_TASK", 0))


# def consecutive_cpg(row, direction):
#     if int(row["nb_sig_cpg"]) > 1:
#         flat_cpg = json_normalize(row["cpg"])
#         max_nb_consec = 0
#         current_nb_consec = 0
#         for index, row in flat_cpg.iterrows():
#             if index > 0:
#                 if (
#                     flat_cpg.iloc[index - 1].fisher_pvalue < p_value
#                     and row.fisher_pvalue < p_value
#                     and np.sign(flat_cpg.iloc[index - 1].effect) == direction
#                     and np.sign(row.effect) == direction
#                 ):
#                     if current_nb_consec == 0:
#                         current_nb_consec = 2
#                     else:
#                         current_nb_consec = current_nb_consec + 1
#                     max_nb_consec = max(max_nb_consec, current_nb_consec)
#                 else:
#                     current_nb_consec = 0
#         return max_nb_consec
#     else:
#         return 0


def consecutive_cpg(row, p_value):
    if int(row["nb_sig_cpg"]) > 1:
        flat_cpg = json_normalize(row["cpg"])
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
    try:
        _, pvalue = stats.mannwhitneyu(
            json_normalize(row["ref"]),
            json_normalize(row["alt"]),
            alternative="two-sided",
        )
        return round(pvalue, 5)
    # If the ref and alt datasets are equal or one is included in the other one:
    except ValueError:
        return 1


def main():

    logging.info("Importing the JSON file as a dataframe")
    logging.info(f"bucket name: {bucket_name}")
    logging.info(f"dataset name: {dataset_name}")
    logging.info(f"batch task index: {BATCH_TASK_INDEX}")
    df, file_name = create_df_from_json_for_index_file(
        storage_client, bucket_name, dataset_name, BATCH_TASK_INDEX, 1
    )

    logging.info(f"Head of dataframe: {df.head()}")

    # # load from file-like objects
    # with open(INPUT_FILE) as f:
    #     data = ndjson.load(f)

    # # Converte the JSON file in dataframe
    # df = json_normalize(data)

    ################################## Calculate p-value of asm_region

    # Function to extract Wilcoxon p-value (5-digit rounding)

    logging.info("Computing Wilcoxon p-value")
    df["wilcoxon_pvalue"] = df.apply(wilcoxon_pvalue, axis=1)

    ################################## Calculate p-value corrected for multiple testing using Benjamini–Hochberg

    df["wilcoxon_corr_pvalue"] = mt.multipletests(
        df["wilcoxon_pvalue"], alpha=bh_threshold, method="fdr_bh"
    )[1]
    df["wilcoxon_corr_pvalue"] = df["wilcoxon_corr_pvalue"].round(5)

    logging.info(
        "Calculate number of significant consecutive CpGs in the same direction."
    )

    results = df.apply(lambda row: consecutive_cpg(row, p_value), axis=1)

    df["nb_consec_pos_sig_asm"] = results.apply(lambda x: x["positive"])
    df["nb_consec_neg_sig_asm"] = results.apply(lambda x: x["negative"])

    logging.info("Saving table to BigQuery")
    upload_dataframe_to_bq(bq_client, df, "wilcoxon")
    # Save to BigQuery table
    # df.to_json(OUTPUT_FILE, orient="records", lines=True)


# Start script
if __name__ == "__main__":
    try:
        main()
    except Exception as err:
        message = f"Task #{BATCH_TASK_INDEX} failed: {str(err)}"

        print(json.dumps({"message": message, "severity": "ERROR"}))
        sys.exit(1)  # Retry Job Task by exiting the process
