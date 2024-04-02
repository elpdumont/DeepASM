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
from google.cloud import bigquery, storage

# from google.cloud.exceptions import NotFound

bq_client = bigquery.Client()
storage_client = storage.Client()


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


def upload_blob(bucket_name, source_file_name, folder_path):
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
    # Create the full destination path
    destination_blob_name = (
        os.path.join(folder_path, file_name) if folder_path else file_name
    )
    # Create a new blob and upload the file's content
    blob = bucket.blob(destination_blob_name)

    blob.upload_from_filename(source_file_name)

    print(f"File {source_file_name} uploaded to {destination_blob_name}.")


def create_df_from_json_for_index_file(
    bucket_name, folder_path, task_index, num_files_to_download
):
    # Initialize the GCP Storage client
    storage_client = storage.Client()

    # Define the prefix to search within a specific folder, ensuring it ends with '/'
    prefix = folder_path if folder_path.endswith("/") else f"{folder_path}/"

    # List all blobs in the specified folder
    blobs = storage_client.list_blobs(bucket_name, prefix=prefix)

    # Filter blobs that match the '*.json' pattern and sort them
    filtered_blobs = sorted(
        (blob for blob in blobs if blob.name.endswith(".json")),
        key=lambda blob: blob.name,
    )

    # Calculate start and end index for files to download
    start_index = task_index * num_files_to_download
    end_index = start_index + num_files_to_download

    # Handle the case for the last task index if we have less than the number of files to download
    if start_index >= len(filtered_blobs):
        # logging.info("Task index is out of range.")
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
    df, bucket_name, folder_path, index, file_name_base
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

    # Upload the JSON string
    blob.upload_from_string(json_string, content_type="application/json")

    print(f"DataFrame exported to: gs://{bucket_name}/{file_name}")


def upload_dataframe_to_bq(bq_client, dataframe, table_id, schema=None):
    logging.info(f"Uploading dataframe to {table_id}")
    if not schema:
        job_config = bigquery.LoadJobConfig(autodetect=True)
    else:
        job_config = bigquery.LoadJobConfig(schema=schema)

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
