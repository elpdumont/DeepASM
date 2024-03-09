# # Kernel functions
from sklearn.neighbors import KernelDensity
from numpy import asarray
from numpy import exp
from flask import Flask, request, jsonify
import numpy as np
import pandas as pd
import dask.dataframe as dd
import os
import json
from google.cloud import storage
import logging


app = Flask(__name__)

# Initialize the Google Cloud Storage client
storage_client = storage.Client()

# Initialize logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')


@app.route('/process', methods=['GET'])
def process_file():
    logging.info("Processing request")
    # Get the bucket name and file path from the request
    bucket_name = request.args.get('bucket')
    file_path = request.args.get('file_path')

    logging.info("")

    if bucket_name and file_path:
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
    else:
        error_message = "Bucket name or file path not provided"
        logging.error(error_message)
        return jsonify({"error": "Bucket name or file path not provided"}), 400

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))
