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
TASK_INDEX = os.getenv("CLOUD_RUN_TASK_INDEX", 0)
TASK_ATTEMPT = os.getenv("CLOUD_RUN_TASK_ATTEMPT", 0)

# Retrieve User-defined env vars
BUCKET_NAME = os.getenv("BUCKET_NAME")
FILE_PATH = os.getenv("FILE_PATH")


# Define main script
def main(bucket_name, file_path):
    """Program that simulates work using the sleep method and random failures.

    Args:
        sleep_ms: number of milliseconds to sleep
        fail_rate: rate of simulated errors
    """

    logging.info(f"Fetching file from bucket: {bucket_name}, file path: {file_path}")

    print(f"Fetching file from bucket: {bucket_name}, file path: {file_path}")

    # Get the bucket and blob (file) from Google Cloud Storage
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(file_path)

    # Download the file as bytes and decode it to a string
    file_contents = blob.download_as_bytes().decode('utf-8')

    # Process each line as a separate JSON object
    processed_data = []
    for line in file_contents.splitlines():
        data = json.loads(line)  # Parse each line as JSON
        processed_data.append(data)  # Replace this with your actual processing logic

    logging.info(processed_data[0])

    return processed_data[0] #jsonify(processed_data), 200
   
# Start script
if __name__ == "__main__":
    try:
        main(BUCKET_NAME, FILE_PATH)
    except Exception as err:
        message = (
            f"Task #{TASK_INDEX}, " + f"Attempt #{TASK_ATTEMPT} failed: {str(err)}"
        )

        print(json.dumps({"message": message, "severity": "ERROR"}))
        sys.exit(1)  # Retry Job Task by exiting the process
