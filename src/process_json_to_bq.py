# # Kernel functions
import json
import logging

# File and variable management
import os
import sys

# Python packages for data, stats
from google.cloud import bigquery, storage

# Initialize the Google Cloud Storage client
storage_client = storage.Client()
bq_client = bigquery.Client()

# Initialize logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")


# Retrieve Job-defined env vars
TASK_INDEX = int(os.getenv("CLOUD_RUN_TASK_INDEX", 0))
TASK_ATTEMPT = os.getenv("CLOUD_RUN_TASK_ATTEMPT", 0)

# Retrieve User-defined env vars
BUCKET_NAME = os.getenv("BUCKET_NAME")
DATASET_TYPE = os.getenv("DATASET_TYPE")
BUCKET_FOLDER_PATH = os.getenv("BUCKET_FOLDER_PATH")
BQ_ML_DATASET_NAME = os.getenv("BQ_ML_DATASET_NAME")
BQ_ML_TABLE_NAME = os.getenv("BQ_ML_TABLE_NAME")


# Function to load JSON files into BigQuery
def load_json_to_bigquery(file_uris, bq_dataset_name, bq_table_name):
    """Loads JSON files into BigQuery."""

    table_id = f"{bq_client.project}.{bq_dataset_name}.{bq_table_name}"

    job_config = bigquery.LoadJobConfig()
    job_config.source_format = bigquery.SourceFormat.NEWLINE_DELIMITED_JSON
    job_config.autodetect = True
    job_config.write_disposition = bigquery.WriteDisposition.WRITE_APPEND

    for uri in file_uris:
        load_job = bq_client.load_table_from_uri(uri, table_id, job_config=job_config)
        logging.info(f"Starting job {load_job.job_id} for file {uri}")

        load_job.result()  # Waits for the job to complete.
        logging.info(f"Job finished for file {uri}")

        destination_table = bq_client.get_table(table_id)
        logging.info(f"Loaded {destination_table.num_rows} rows into {table_id}")


# Main process
def main():

    # Delete the table if it exists
    delete_table_if_exists(BQ_ML_DATASET_NAME, BQ_ML_TABLE_NAME)

    # Create the table as empty
    table_id = f"{bq_client.project}.{BQ_ML_DATASET_NAME}.{BQ_ML_TABLE_NAME}"
    table = bigquery.Table(table_id)
    table = bq_client.create_table(table)

    # Process files
    files_to_load = list_files(BUCKET_NAME, BUCKET_FOLDER_PATH, DATASET_TYPE)
    if files_to_load:
        load_json_to_bigquery(files_to_load, BQ_ML_DATASET_NAME, BQ_ML_TABLE_NAME)
    else:
        logging.info("No files found matching the criteria.")


# Start script
if __name__ == "__main__":
    try:
        main()
    except Exception as err:
        message = (
            f"Task #{TASK_INDEX}, " + f"Attempt #{TASK_ATTEMPT} failed: {str(err)}"
        )

        print(json.dumps({"message": message, "severity": "ERROR"}))
        sys.exit(1)  # Retry Job Task by exiting the process
