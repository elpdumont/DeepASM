
# import numpy as np
# import pandas as pd
# import dask.dataframe as dd


# # Kernel functions
# from sklearn.neighbors import KernelDensity
# from numpy import asarray
# from numpy import exp


from flask import Flask, request, jsonify
import os
import json

app = Flask(__name__)

@app.route('/process', methods=['GET'])
def process_file():
    # Your specific Python job logic here
    # For example, processing a JSON file:
    filename = request.args.get('filename')
    if filename:
        filepath = os.path.join('path_to_your_files', filename)
        with open(filepath, 'r') as file:
            data = json.load(file)
            # Process your JSON data here
            processed_data = data  # This should be replaced by your actual processing logic

        return jsonify(processed_data), 200
    else:
        return jsonify({"error": "No filename provided"}), 400

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))