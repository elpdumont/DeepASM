# # Kernel functions
from sklearn.neighbors import KernelDensity
from numpy import asarray
from numpy import exp
from flask import Flask, request, jsonify
import numpy as np
import pandas as pd
import dask
import dask.dataframe as dd
import os
import json
from google.cloud import storage
import logging
import argparse



app = Flask(__name__)

# Initialize the Google Cloud Storage client
storage_client = storage.Client()

# Initialize logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

# Disable warning log for dask
dask.config.set({'dataframe.query-planning-warning': False})


@app.route('/process', methods=['GET'])
def process_file():
    logging.info("Processing request")
    # Get the bucket name and file path from the request
    global args
    bucket_name = args.bucket
    file_path = args.file_path

    logging.info(f"Fetching file from bucket: {bucket_name}, file path: {file_path}")

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

    return processed_data[0] #jsonify(processed_data), 200

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some JSON files.')
    parser.add_argument('--bucket', type=str, help='The name of the GCS bucket')
    parser.add_argument('--file_path', type=str, help='The path to the JSON file in the GCS bucket')

    args = parser.parse_args()

    app.run(debug=True, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))
