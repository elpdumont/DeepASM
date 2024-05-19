import json
import logging

# File and variable management
import os
import random

# import re
import sys
import time

# Python packages for data, stats
import pandas as pd

# import pandas_gbq
from google.api_core.exceptions import Forbidden, TooManyRequests
from google.cloud import bigquery

# Create a handler for Google Cloud Logging.
logging.basicConfig(level=logging.INFO)


# def delete_bq_table(bq_dataset_id, bq_table_id):
#     """Deletes the BigQuery table if it exists."""
#     table_id = f"{bq_client.project}.{bq_dataset_id}.{bq_table_id}"
#     try:
#         bq_client.delete_table(table_id)
#         logging.info(f"Deleted table '{table_id}'.")
#     except NotFound:
#         logging.info(f"Table {table_id} not found, no deletion performed.")


# def append_df_to_bq_table(df, bq_dataset_id, bq_table_id):
#     pandas_gbq.to_gbq(
#         df,
#         f"{bq_dataset_id}.{bq_table_id}",
#         if_exists="append",
#         progress_bar=True,
#     )


def convert_types(d, type_map):
    return {
        key: type_map[key](value) if key in type_map else value
        for key, value in d.items()
    }


def upload_blob(storage_client, bucket_name, source_file_name, folder_path):
    """
    Uploads a file to a specified folder within a Google Cloud Storage bucket, keeping the original file name.

    Parameters:
    - bucket_name: str. The name of the bucket to upload to.
    - source_file_name: str. The path to the file to upload.
    - folder_path: str. The folder path within the bucket where the file will be uploaded.

    Returns:
    None
    """

    # Get the bucket object
    bucket = storage_client.bucket(bucket_name)
    # Extract the file name from the source file path
    file_name = os.path.basename(source_file_name)

    # logging.info(f"File name {file_name}")
    # Create the full destination path
    destination_blob_name = (
        os.path.join(folder_path, file_name) if folder_path else file_name
    )

    # logging.info(f"Destination blob name: {destination_blob_name}")
    # Create a new blob and upload the file's content
    blob = bucket.blob(destination_blob_name)

    blob.upload_from_filename(source_file_name)

    logging.info(f"File {source_file_name} uploaded to {destination_blob_name}.")


def download_blob(storage_client, bucket_name, source_blob_name, destination_file_path):
    """
    Downloads a file from a specified folder within a Google Cloud Storage bucket.
    Parameters:
    - storage_client: storage.Client. The Google Cloud Storage client.
    - bucket_name: str. The name of the bucket to download from.
    - source_blob_name: str. The blob name in the bucket to download.
    - destination_file_path: str. The local path where the file will be saved.
    Returns:
    None
    """
    # Get the bucket object
    bucket = storage_client.bucket(bucket_name)
    # Create a new blob object
    blob = bucket.blob(source_blob_name)
    # Download the file to the destination path
    blob.download_to_filename(destination_file_path)
    logging.info(f"File {source_blob_name} downloaded to {destination_file_path}.")


def create_df_from_json_for_index_file(
    storage_client, bucket_name, folder_path, task_index, num_files_to_download
):

    # Define the prefix to search within a specific folder, ensuring it ends with '/'
    prefix = folder_path if folder_path.endswith("/") else f"{folder_path}/"

    # List all blobs in the specified folder
    blobs = storage_client.list_blobs(bucket_name, prefix=prefix)

    # Filter blobs that match the '*.json' pattern and sort them
    filtered_blobs = sorted(
        (blob for blob in blobs if blob.name.endswith(".json")),
        key=lambda blob: blob.name,
    )

    # logging.info(f"Filtered blobs: {filtered_blobs}")

    # Calculate start and end index for files to download
    start_index = task_index * num_files_to_download
    end_index = start_index + num_files_to_download

    # Handle the case for the last task index if we have less than the number of files to download
    if start_index >= len(filtered_blobs):
        logging.info("Task index is out of range.")
        sys.exit(1)
    if end_index > len(filtered_blobs):
        end_index = len(filtered_blobs)

    # Process and download files for the given range
    all_data = []
    file_names = []
    for blob in filtered_blobs[start_index:end_index]:
        # logging.info(f"Processing file: {blob.name}")

        # Download the file as bytes and decode it to a string
        file_contents = blob.download_as_bytes().decode("utf-8")

        # Process each line as a separate JSON object
        processed_data = [json.loads(line) for line in file_contents.splitlines()]
        all_data.extend(processed_data)
        file_names.append(os.path.basename(blob.name))

    # Combine all data into a single DataFrame
    combined_data = pd.DataFrame(all_data)
    return combined_data, file_names


def export_dataframe_to_gcs_as_json(
    storage_client, df, bucket_name, folder_path, index, file_name_base
):
    """
    Export a pandas DataFrame to a Google Cloud Storage bucket as a JSON file.

    Parameters:
    - df: pandas.DataFrame to be exported.
    - bucket_name: Name of the GCS bucket.
    - folder_path: Path within the bucket to save the file (excluding the file name).
    - index: An index to include in the file name.
    - file_name_base: A base string to include in the file name.

    Returns:
    - None
    """

    # Get the bucket
    bucket = storage_client.bucket(bucket_name)

    # Convert the DataFrame to a JSON string
    json_string = df.to_json()

    # Construct the full path including the file name
    # Ensuring folder_path ends with '/'
    if not folder_path.endswith("/"):
        folder_path += "/"
    file_name = f"{folder_path}{file_name_base}_{index}.json"

    # Create a blob (GCS object) for the file
    blob = bucket.blob(file_name)

    if blob.exists():
        logging.info(
            f"File gs://{bucket_name}/{file_name} already exists and will be overwritten."
        )

    # Upload the JSON string
    blob.upload_from_string(json_string, content_type="application/json")

    logging.info(f"DataFrame exported to: gs://{bucket_name}/{file_name}")


def upload_dataframe_to_bq(
    bq_client,
    dataframe,
    table_id,
    schema=None,
    partition_field=None,
    cluster_fields=None,
):
    """
    Upload a dataframe to Google BigQuery, with options for schema definition, partitioning, and clustering.
    Parameters:
    - bq_client: The BigQuery client.
    - dataframe: The Pandas DataFrame to upload.
    - table_id: The BigQuery table ID where the DataFrame will be uploaded.
    - schema: The schema of the table, if not set, schema will be autodetected.
    - partition_field: The field name to be used for partitioning.
    - cluster_fields: A list of field names to be used for clustering.
    Returns:
    - None
    """
    logging.info(f"Uploading dataframe to {table_id}")
    job_config = bigquery.LoadJobConfig()
    if not schema:
        job_config.autodetect = True
    else:
        job_config.schema = schema
    job_config.write_disposition = "WRITE_APPEND"
    # Configure partitioning
    if partition_field:
        start, end, interval = (0, 4000, 1)
        job_config.range_partitioning = bigquery.RangePartitioning(
            field=partition_field,
            range_=bigquery.PartitionRange(start=start, end=end, interval=interval),
        )
    # Configure clustering
    if cluster_fields:
        job_config.clustering_fields = cluster_fields
    for attempt in range(1, 7):  # Retry up to 5 times
        try:
            job = bq_client.load_table_from_dataframe(
                dataframe, table_id, job_config=job_config
            )
            result = job.result()  # Wait for the job to complete
            logging.info(f"Load job result for {table_id}: {result}")
            break  # Success, exit the retry loop
        except TooManyRequests:
            logging.warning("Caught TooManyRequests; applying backoff.")
        except Forbidden as e:
            if "rateLimitExceeded" in str(e):
                logging.warning(
                    "Caught Forbidden with rateLimitExceeded; applying backoff."
                )
            else:
                logging.error("Forbidden error not related to rate limits.")
                raise
        except Exception as e:
            logging.error(f"An unexpected error occurred: {e}")
            raise
        # Handle retry for both TooManyRequests and rateLimitExceeded Forbidden errors
        if attempt < 6:
            base_sleep = 2**attempt  # Exponential backoff formula
            random_sleep = random.uniform(
                0, 4
            )  # Add randomness between 0 and 3 seconds
            sleep_time = base_sleep + random_sleep
            logging.info(f"Rate limit exceeded. Retrying in {sleep_time} seconds.")
            time.sleep(sleep_time)
        else:
            logging.error("Maximum retry attempts reached. Job failed.")
            raise Exception("Maximum retry attempts reached.")


def fetch_chunk_from_bq_as_dataframe_w_hmmvar(
    dataset_id,
    table_id,
    task_index,
    total_tasks,
    project_id,
    hmm_var,
    ml_mode,
    ml_nb_datapoints_for_testing,
):
    client = bigquery.Client(project=project_id)
    # Construct SQL to divide the table into chunks
    query = f"""
    SELECT *
    FROM `{dataset_id}.{table_id}`
    WHERE {hmm_var} IS NOT NULL AND MOD(ABS(FARM_FINGERPRINT(CAST(clustering_index AS STRING))), {total_tasks}) = {task_index}
    """
    if ml_mode == "TESTING":
        logging.info("In testing mode. Adding a limit to the import.")
        query += f"LIMIT {ml_nb_datapoints_for_testing}"

    # Start the query and wait for it to complete, then load it into a DataFrame
    query_job = client.query(query)
    dataframe = query_job.to_dataframe()
    return dataframe
