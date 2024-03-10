# # Kernel functions
from sklearn.neighbors import KernelDensity
from numpy import asarray
from numpy import exp
import numpy as np
import pandas as pd
import dask
import dask.dataframe as dd
import os
import random
import json
from google.cloud import storage
import logging
import sys


# Initialize the Google Cloud Storage client
storage_client = storage.Client()

# Initialize logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

# Disable warning log for dask
dask.config.set({'dataframe.query-planning-warning': False})


# Retrieve Job-defined env vars
TASK_INDEX = int(os.getenv("CLOUD_RUN_TASK_INDEX", 0))
TASK_ATTEMPT = os.getenv("CLOUD_RUN_TASK_ATTEMPT", 0)

# Retrieve User-defined env vars
BUCKET_NAME = os.getenv("BUCKET_NAME")
FOLDER_PATH = os.getenv("FOLDER_PATH")


# Define main script
def main(bucket_name, folder_path, max_digits = 12):

    # Define the prefix to search within a specific folder, ensuring it ends with '/'
    prefix = folder_path if folder_path.endswith('/') else f"{folder_path}/"

    # List all blobs in the specified folder
    blobs = storage_client.list_blobs(bucket_name, prefix=prefix)

    # Filter blobs that match the 'raw-*.json' pattern and sort them
    filtered_blobs = sorted(
        (blob for blob in blobs if blob.name.startswith(prefix + "raw-") and blob.name.endswith(".json")),
        key=lambda blob: blob.name
    )

    # Ensure task_index is within the range of available files
    if 0 <= TASK_INDEX < len(filtered_blobs):
        # Get the file name of the blob at the specified index
        blob_name = filtered_blobs[TASK_INDEX]
    else:
        logging.info(f"Task index {TASK_INDEX} is out of range.")
        sys.exit(1)
    # Get the bucket and blob (file) from Google Cloud Storage
    
    logging.info(f"Processing the bucket {bucket_name} with the folder path {folder_path} and the file name {blob_name.name}")

    bucket = storage_client.bucket(bucket_name)

    # Download the file as bytes and decode it to a string
    file_contents = blob_name.download_as_bytes().decode('utf-8')

    # Process each line as a separate JSON object
    processed_data = []
    for line in file_contents.splitlines():
        data = json.loads(line)  # Parse each line as JSON
        processed_data.append(data)  # Replace this with your actual processing logic

    logging.info(processed_data[0])

    return processed_data[0] 

# Start script
if __name__ == "__main__":
    try:
        main(BUCKET_NAME, FOLDER_PATH)
    except Exception as err:
        message = (
            f"Task #{TASK_INDEX}, " + f"Attempt #{TASK_ATTEMPT} failed: {str(err)}"
        )

        print(json.dumps({"message": message, "severity": "ERROR"}))
        sys.exit(1)  # Retry Job Task by exiting the process
