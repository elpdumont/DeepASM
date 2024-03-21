import json
import logging

# File and variable management
import os

# import re
import sys

# Python packages for data, stats
import pandas as pd
import pandas_gbq
from google.cloud import bigquery, storage
from google.cloud.exceptions import NotFound

bq_client = bigquery.Client()
storage_client = storage.Client()


def delete_bq_table(bq_dataset_id, bq_table_id):
    """Deletes the BigQuery table if it exists."""
    table_id = f"{bq_client.project}.{bq_dataset_id}.{bq_table_id}"
    try:
        bq_client.delete_table(table_id)
        logging.info(f"Deleted table '{table_id}'.")
    except NotFound:
        logging.info(f"Table {table_id} not found, no deletion performed.")


def append_df_to_bq_table(df, bq_dataset_id, bq_table_id):
    pandas_gbq.to_gbq(
        df,
        f"{bq_dataset_id}.{bq_table_id}",
        if_exists="append",
        progress_bar=True,
    )


# def export_df_to_gcs_json(df, data_type, bucket_name, bucket_path, file_name):
#     # Convert DataFrame to JSON format
#     json_data = df.to_json(orient="records", lines=True)

#     file_name = data_type + "_" + file_name
#     logging.info(f"Saving JSON data as {file_name}")
#     with open(file_name, "w") as file:
#         file.write(json_data)

#     # Cloud filename includes the bucket_path
#     cloud_filename = os.path.join(bucket_path, file_name)

#     # Initialize a Google Cloud Storage client
#     storage_client = storage.Client()

#     # Get the bucket
#     bucket = storage_client.bucket(bucket_name)

#     # Create a blob (GCS object) and upload the file
#     blob = bucket.blob(cloud_filename)
#     blob.upload_from_filename(file_name)

#     logging.info(f"File {file_name} uploaded to {cloud_filename}.")


# def list_json_files_in_bucket_folder(bucket_name, folder_path):
#     """Lists all the JSON files in the specified bucket folder."""
#     prefix = folder_path if folder_path.endswith("/") else f"{folder_path}/"
#     files = storage_client.list_blobs(bucket_name, prefix=prefix)
#     pattern = re.compile(rf"{re.escape(prefix)}.*\.json$")
#     return [
#         f"gs://{bucket_name}/{file.name}" for file in files if pattern.match(file.name)
#     ]


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
        logging.info("Task index is out of range.")
        sys.exit(1)
    if end_index > len(filtered_blobs):
        end_index = len(filtered_blobs)

    # Process and download files for the given range
    all_data = []
    file_names = []
    for blob in filtered_blobs[start_index:end_index]:
        logging.info(f"Processing file: {blob.name}")

        # Download the file as bytes and decode it to a string
        file_contents = blob.download_as_bytes().decode("utf-8")

        # Process each line as a separate JSON object
        processed_data = [json.loads(line) for line in file_contents.splitlines()]
        all_data.extend(processed_data)
        file_names.append(os.path.basename(blob.name))

    # Combine all data into a single DataFrame
    combined_data = pd.DataFrame(all_data)
    return combined_data, file_names


# def create_df_from_json_for_index_file(bucket_name, folder_path, task_index):
#     # Initialize the GCP Storage client
#     storage_client = storage.Client()

#     # Define the prefix to search within a specific folder, ensuring it ends with '/'
#     prefix = folder_path if folder_path.endswith("/") else f"{folder_path}/"

#     # List all blobs in the specified folder
#     blobs = storage_client.list_blobs(bucket_name, prefix=prefix)

#     # Filter blobs that match the '*.json' pattern and sort them
#     filtered_blobs = sorted(
#         (blob for blob in blobs if blob.name.endswith(".json")),
#         key=lambda blob: blob.name,
#     )
#     # Ensure task_index is within the range of available files
#     if 0 <= task_index < len(filtered_blobs):
#         # Get the Blob object at the specified index
#         selected_blob = filtered_blobs[task_index]
#     else:
#         logging.info(f"Task index {task_index} is out of range.")
#         sys.exit(1)

#     # Log the processing info
#     logging.info(
#         f"Processing the bucket {bucket_name} with the folder path {folder_path} and the file name {selected_blob.name}"
#     )

#     logging.info("Download the file as bytes and decode it to a string")
#     file_contents = selected_blob.download_as_bytes().decode("utf-8")

#     logging.info("Process each line as a separate JSON object")
#     processed_data = []
#     for line in file_contents.splitlines():
#         data = json.loads(line)  # Parse each line as JSON
#         processed_data.append(data)

#     return pd.DataFrame(processed_data), os.path.basename(selected_blob.name)
