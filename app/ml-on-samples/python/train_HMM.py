# Python packages for data, stats
import ast
import json
import logging

# import multiprocessing
import os
import random
import sys

import numpy as np
import pandas as pd
import yaml
from gcp import upload_blob, upload_dataframe_to_bq
from google.cloud import bigquery, storage
from hmmlearn.hmm import GaussianHMM
from hmmlearn.vhmm import VariationalGaussianHMM

# from hmmlearn.vhmm import VariationalGaussianHMM
from joblib import dump
from sklearn.utils import check_random_state

# Initialize the Google Cloud Storage client
storage_client = storage.Client()
bq_client = bigquery.Client()

# Create a handler for Google Cloud Logging.
logging.basicConfig(level=logging.INFO)

# Required to avoid bugs.
# os.environ["OPENBLAS_NUM_THREADS"] = "16"


# Import all other variables from the config file
with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

# Accessing GCP configuration
project = config["GCP"]["PROJECT"]
bucket = config["GCP"]["BUCKET"]

# Samples for splitting the dataset
samples_dic = config["SAMPLES"]
dataset_types = list(samples_dic.keys())

# Retrieve Job-defined env vars
ml_dataset = os.getenv("ML_DATASET")
model_path = os.getenv("MODEL_PATH")
short_sha = os.getenv("SHORT_SHA")
home_directory = os.path.expanduser("~")
ml_mode = os.getenv("ML_MODE")

# ML/HMM variables
ml_nb_datapoints_for_testing = config["ML"]["NB_DATA_POINTS_TESTING"]
n_states = config["ML"]["HMM"]["N_STATES"]
model_type_str = config["ML"]["HMM"]["MODEL_TYPE"]
covariance = config["ML"]["HMM"]["COVARIANCE"]
algorithm = config["ML"]["HMM"]["ALGORITHM"]
hmm_var = config["ML"]["HMM"]["VAR_NAME"]
n_model_loop = config["ML"][ml_mode]["HMM_N_MODEL_LOOP"]
n_iterations = config["ML"][ml_mode]["N_HMM_ITERATIONS"]
hmm_n_clusters = config["ML"][ml_mode]["HMM_N_CLUSTERS"]

# Initialize random state
base_seed = 546  # Example value, adjust as needed
rs = check_random_state(base_seed)
# andom_seeds = [base_seed + i for i in range(n_model_loop)]

dic_model = {
    "VariationalGaussianHMM": VariationalGaussianHMM,
    "GaussianHMM": GaussianHMM,
}
model_type = dic_model[model_type_str]

# Initialize client
credentials_path = "/appuser/.config/gcloud/application_default_credentials.json"
if os.path.exists(credentials_path):
    # Set the GOOGLE_APPLICATION_CREDENTIALS environment variable to the path of the JSON key file
    os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = credentials_path
    # Assuming 'project' is already defined somewhere in your script
    os.environ["GOOGLE_CLOUD_PROJECT"] = project
    # ml_mode = "TESTING"

# Initialize the Google Cloud Storage client
storage_client = storage.Client()
bq_client = bigquery.Client()


def generate_random_integers(n_clusters):
    # Generate a random seed
    seed = random.randint(0, 10000)
    random.seed(seed)
    # Generate a list of N random integers between 0 and 3601
    random_integers = [random.randint(0, 3601) for _ in range(n_clusters)]
    return random_integers


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


def save_HMM_model_to_bucket(
    directory,
    model_name,
    ml_dataset,
    model,
    short_sha,
    bucket,
    model_path,
    n_states,
    covariance,
    ml_mode,
):
    file_name = (
        directory
        + model_name
        + "_"
        + str(n_states)
        + "states_"
        + str(covariance)
        + "_"
        + short_sha
        + "_"
        + ml_dataset
        + "_"
        + ml_mode
        + ".joblib"
    )
    dump(model, file_name)
    upload_blob(storage_client, bucket, file_name, model_path)
    return None


def fit_hmm(
    model_type,
    n_states,
    n_iter,
    covariance,
    algorithm,
    reshaped_data,
    lengths,
    rs,
):
    # logging.info(
    #     f"Number of states: {n_states} and random seed for generation: {rs_seed}"
    # )
    # rs = check_random_state(rs_seed)
    h = model_type(
        n_components=n_states,
        n_iter=n_iterations,
        covariance_type=covariance,
        algorithm=algorithm,
        random_state=rs,
    )
    h.fit(reshaped_data, lengths)
    score = h.score(reshaped_data)
    return score, h


def main():
    logging.info(f"ML MODE: {ml_mode}")
    logging.info(f"Config file: {config}")
    dataset_for_hmm = "TRAINING"
    quoted_samples = ",".join(
        [f"'{sample}'" for sample in samples_dic[dataset_for_hmm]]
    )
    logging.info(f"Preparing a query with these samples: {quoted_samples}")
    list_rand_int = generate_random_integers(hmm_n_clusters)
    quoted_list_rand = ", ".join(
        map(str, list_rand_int)
    )  # Convert integers to string for query
    query = f"""
        SELECT {hmm_var}
        FROM {project}.{ml_dataset}.features_wo_hmm
        WHERE
            sample IN ({quoted_samples}) AND
            {hmm_var} IS NOT NULL AND
            clustering_index in ({quoted_list_rand})
            """
    logging.info("Importing dataset...")
    df = bq_client.query(query).to_dataframe()
    # logging.info("Randomizing the dataset...")
    # df = df.sample(frac=1, random_state=base_seed).reset_index(drop=True)
    logging.info(f"Number of rows in the DF: {len(df):,}")

    logging.info("Converting the sequence from str to floats...")
    df["cpg_directional_fm"] = df["cpg_directional_fm"].apply(
        lambda x: ast.literal_eval(x.strip('"'))
    )
    logging.info("Creating a unique sequence for training the HMM")
    all_obs = np.concatenate(df["cpg_directional_fm"].tolist())
    nb_cpgs_in_training = len(all_obs)
    logging.info(f"Number of CpGs to be used in training: {nb_cpgs_in_training:,}")
    reshaped_data, lengths = prepare_data_for_hmm(all_obs)

    # pool = multiprocessing.Pool(processes=n_model_loop)
    best_ll = None
    best_model = None
    # for rs_seed in random_seeds:
    for i in range(n_model_loop):
        logging.info(f"Iteration: {i}")
        # results.append(
        # pool.apply_async(
        score, h = fit_hmm(
            model_type,
            n_states,
            n_iterations,
            covariance,
            algorithm,
            reshaped_data,
            lengths,
            rs,
        )
        logging.info(f"Score of iteration {i}: {score:,}")
        if best_ll is None or score > best_ll:
            logging.info("Found new a best score!")
            best_ll = score
            best_model = h

    # aic, bic, ll = best_model.aic(reshaped_data), best_model.bic(reshaped_data), best_ll
    # Store values in DF
    df = pd.DataFrame(
        {
            "model_name": [model_type_str],
            "n_states": [n_states],
            "covariance": [covariance],
            "short_sha": [short_sha],
            "n_model_loop": [n_model_loop],
            "n_cpgs": [nb_cpgs_in_training],
            "n_iterations": [n_iterations],
            # "aic": [aic],
            # "bic": [bic],
            "log_loss": [np.round(best_ll)],
            "ml_mode": [ml_mode],
        }
    )
    logging.info("Upload results to BQ")
    upload_dataframe_to_bq(bq_client, df, f"{ml_dataset}.hmm_results")
    logging.info("Save model to bucket")
    save_HMM_model_to_bucket(
        home_directory,
        model_name,
        ml_dataset,
        best_model,
        short_sha,
        bucket,
        model_path,
        n_states,
        covariance,
        ml_mode,
    )
    logging.info(f"END OF SCRIPT for {n_states} states")


# Start script
if __name__ == "__main__":
    try:
        main()
    except Exception as err:
        message = f"Number of states #{n_states} failed: {str(err)}"

        print(json.dumps({"message": message, "severity": "ERROR"}))
        sys.exit(1)  # Retry Job Task by exiting the process
