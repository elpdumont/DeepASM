# Python packages for data, stats
import ast
import json
import logging
import os
import sys

import numpy as np
import pandas as pd
import yaml
from gcp import upload_blob, upload_dataframe_to_bq
from google.cloud import bigquery, storage
from hmmlearn.hmm import GaussianHMM

# from hmmlearn.vhmm import VariationalGaussianHMM
from joblib import dump
from sklearn.utils import check_random_state

# Initialize the Google Cloud Storage client
storage_client = storage.Client()
bq_client = bigquery.Client()

# Initialize random state
rs = check_random_state(546)

# Create a handler for Google Cloud Logging.
logging.basicConfig(level=logging.INFO)

# Import all other variables from the config file
with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

# Accessing GCP configuration
project = config["GCP"]["PROJECT"]
bucket = config["GCP"]["BUCKET"]

# Samples for splitting the dataset
samples_dic = config["SAMPLES"]
dataset_types = list(samples_dic.keys())

# HMM variables
model_type = config["HMM"]["MODEL_TYPE"]
n_model_loop = config["HMM"]["N_MODEL_LOOP"]
covariance = config["HMM"]["COVARIANCE"]
n_iterations = config["HMM"]["N_ITERATIONS"]
algorithm = config["HMM"]["ALGORITHM"]
hmm_var = config["HMM"]["VAR_NAME"]

# model_dic = {
#     "Gaussian": GaussianHMM,
#     "Variational Gaussian": VariationalGaussianHMM,
# }

# Retrieve Job-defined env vars
ml_dataset = os.getenv("ML_DATASET")
model_path = os.getenv("MODEL_PATH")
short_sha = os.getenv("SHORT_SHA")
n_states = int(os.getenv("BATCH_TASK_INDEX", 0)) + 2
home_directory = os.path.expanduser("~")

# Initialize client
credentials_path = "/appuser/.config/gcloud/application_default_credentials.json"
if os.path.exists(credentials_path):
    # Set the GOOGLE_APPLICATION_CREDENTIALS environment variable to the path of the JSON key file
    os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = credentials_path
    # Assuming 'project' is already defined somewhere in your script
    os.environ["GOOGLE_CLOUD_PROJECT"] = project

# Initialize the Google Cloud Storage client
storage_client = storage.Client()
bq_client = bigquery.Client()


def prepare_data_for_hmm(sequence):
    """
    Prepares a sequence of data for processing with a Hidden Markov Model (HMM).

    This function ensures the input sequence is in a 2D NumPy array format required by HMM processing routines, handling both single-dimensional sequences (interpreting them as a sequence of scalar observations) and two-dimensional sequences (interpreting them as a sequence of vector observations). It also calculates the length of the sequence, which is necessary for some HMM algorithms.

    Parameters:
    - sequence (np.ndarray): The input sequence to be processed. This can be either a 1D array of scalar observations or a 2D array of vector observations, where each row represents a timestep.

    Returns:
    - sequence (np.ndarray): The input sequence reshaped into a 2D NumPy array format, with individual observations along rows.
    - lengths (list of int): A list containing a single integer, which is the length of the input sequence. This is used by HMM algorithms that require the lengths of sequences being processed.

    Raises:
    - ValueError: If the input `sequence` has more than two dimensions, indicating it's not in an acceptable format for HMM processing.
    """
    if sequence.ndim == 1:
        sequence = np.atleast_2d(sequence).T
    elif sequence.ndim > 2:
        raise ValueError(
            "Sequence must be 1D (for single float sequence) or 2D (for sequence of vectors)"
        )
    # Determine the length of the sequence dynamically
    sequence_length = sequence.shape[0]
    # For a single sequence, the lengths list contains just one element: the sequence length
    lengths = [sequence_length]
    return sequence, lengths


def save_HMM_model_to_bucket(directory, model, short_sha, bucket, model_path, n_states):
    file_name = (
        directory + "/hmm_model_" + short_sha + "_" + str(n_states) + "states.joblib"
    )
    dump(model, file_name)
    upload_blob(storage_client, bucket, file_name, model_path)
    return None


def main():

    dataset_for_hmm = "TRAINING"
    quoted_samples = ",".join(
        [f"'{sample}'" for sample in samples_dic[dataset_for_hmm]]
    )
    query = f"SELECT {hmm_var} FROM {project}.{ml_dataset}.features_wo_hmm WHERE sample IN ({quoted_samples}) AND {hmm_var} IS NOT NULL"
    # Execute the query and store in dic
    df = bq_client.query(query).to_dataframe()
    df["cpg_directional_fm"] = df["cpg_directional_fm"].apply(
        lambda x: ast.literal_eval(x.strip('"'))
    )
    logging.info("Creating a unique sequence for training the HMM")
    all_obs = np.concatenate(df["cpg_directional_fm"].tolist())
    logging.info(f"Number of CpGs to be used in training: {len(all_obs)}")
    reshaped_data, lengths = prepare_data_for_hmm(all_obs)
    best_ll = None
    best_model = None
    for i in range(n_model_loop):
        logging.info(f"Number of states: {n_states} and iteration: {i}")
        h = GaussianHMM(
            n_states,
            n_iter=n_iterations,
            covariance_type=covariance,
            algorithm=algorithm,
            random_state=rs,
        )
        h.fit(reshaped_data, lengths)
        score = h.score(reshaped_data)
        if best_ll is None or score > best_ll:
            best_ll = score
            best_model = h
    logging.info("Obtain aic, bic, ll")
    aic, bic, ll = best_model.aic(reshaped_data), best_model.bic(reshaped_data), best_ll
    # Store values in DF
    df = pd.DataFrame(
        {
            "n_states": [n_states],
            "short_sha": [short_sha],
            "n_model_loop": [n_model_loop],
            "aic": [aic],
            "bic": [bic],
            "ll": [ll],
        }
    )
    logging.info("Upload results to BQ")
    upload_dataframe_to_bq(bq_client, df, f"{ml_dataset}.hmm_results")
    logging.info("Save model to bucket")
    save_HMM_model_to_bucket(
        home_directory, best_model, short_sha, bucket, model_path, n_states
    )
    logging.info("END OF SCRIPT")


# Start script
if __name__ == "__main__":
    try:
        main()
    except Exception as err:
        message = f"Number of states #{n_states} failed: {str(err)}"

        print(json.dumps({"message": message, "severity": "ERROR"}))
        sys.exit(1)  # Retry Job Task by exiting the process
