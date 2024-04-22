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
from hmmlearn.vhmm import VariationalGaussianHMM
from joblib import dump
from scipy.stats import entropy
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

# Genomic variables
genomic_length = config["GENOMICS"]["GENOMIC_LENGTH"]

# Samples for splitting the dataset
samples_dic = config["SAMPLES"]
dataset_types = list(samples_dic.keys())

# HMM variables
model_type = config["HMM"]["MODEL_TYPE"]
n_states = config["HMM"]["N_STATES"]
n_model_loop = config["HMM"]["N_MODEL_LOOP"]
covariance = config["HMM"]["COVARIANCE"]
n_iterations = config["HMM"]["N_ITERATIONS"]
algorithm = config["HMM"]["ALGORITHM"]
hmm_var = config["HMM"]["VAR_NAME"]

model_dic = {
    "Gaussian": GaussianHMM,
    "Variational Gaussian": VariationalGaussianHMM,
}

# Retrieve Job-defined env vars
ml_dataset = os.getenv("ML_DATASET")
model_path = os.getenv("MODEL_PATH")
short_sha = os.getenv("SHORT_SHA")
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


def predict_hidden_states_for_sequences(model, sequences, log_frequency=100000):
    """
    Predicts the most likely hidden states for each element in a series of sequences using a trained Hidden Markov Model (HMM).
    Each sequence of observable data points is processed individually to determine the sequence's hidden states. These states represent the underlying process assumed by the HMM.
    Parameters:
    - model: The trained HMM model used for prediction.
    - sequences: A list of sequences, where each sequence is an array of numeric observable data points. Each inner array represents a sequence to be processed.
    - log_frequency: Specifies the interval of sequences processed at which progress is logged. Logging occurs every 'log_frequency' sequences.
    Returns:
    - predicted_states: A list where each element is an array of the predicted hidden states for a corresponding input sequence.
    """
    predicted_states = []
    for i, sequence in enumerate(sequences):
        # Convert sequence to a NumPy array if it's not already one
        sequence = np.asarray(sequence)
        # Ensure seq is in the right shape (n_samples, n_features)
        if sequence.ndim == 1:
            sequence = np.atleast_2d(sequence).T
        elif sequence.ndim > 2:
            raise ValueError(
                "Sequence must be 1D (for single float sequence) or 2D (for sequence of vectors)"
            )
        # Predict the most likely hidden states for the sequence
        states = model.predict(sequence)
        predicted_states.append(states)
        # Log progress
        if (i + 1) % log_frequency == 0:
            logging.info(f"Processed sequence number: {i + 1}")
    return predicted_states


def extract_features(hidden_states_sequences):
    """
    Extracts a comprehensive set of features from sequences of hidden states in a Hidden Markov Model (HMM) and provides descriptive names for each feature.
    Parameters:
    - hidden_states_sequences (list of list of int): A list of sequences, where each sequence is a list of integers representing the hidden states visited by the HMM.
    Returns:
    - Tuple: A tuple containing a 2D NumPy array of extracted features for each sequence and a list of descriptive names for each feature.
    """
    features = []
    feature_names = []  # To store names of the features
    # Determine all unique states across sequences for consistent ordering
    unique_states = sorted({state for seq in hidden_states_sequences for state in seq})
    for state in unique_states:
        feature_names.append(f"count_state_{state}")
        feature_names.append(f"proportion_state_{state}")
    for _, state_i in enumerate(unique_states):
        for _, state_j in enumerate(unique_states):
            feature_names.append(f"transition_from_{state_i}_to_{state_j}")
    feature_names.extend(["start_state", "end_state", "nb_state_changes"])
    for state in unique_states:
        feature_names.append(f"mean_duration_state_{state}")
        feature_names.append(f"variance_duration_state_{state}")
    for _, state_i in enumerate(unique_states):
        for _, state_j in enumerate(unique_states):
            feature_names.append(f"transition_probability_from_{state_i}_to_{state_j}")
    feature_names.append("entropy_state_distribution")
    for seq in hidden_states_sequences:
        seq = np.array(seq)
        sequence_features = []
        # Original features: Counts, Proportions, and Transitions
        for state in unique_states:
            count = np.sum(seq == state)
            proportion = count / len(seq)
            sequence_features.extend([count, proportion])
        # Count transitions for later use in calculating probabilities
        transitions_matrix = np.zeros((len(unique_states), len(unique_states)))
        for i, state_i in enumerate(unique_states):
            for j, state_j in enumerate(unique_states):
                transition_count = np.sum((seq[:-1] == state_i) & (seq[1:] == state_j))
                transitions_matrix[i, j] = transition_count
                sequence_features.append(transition_count)
        start_state = seq[0]
        end_state = seq[-1]
        sequence_features.extend([start_state, end_state])
        # New features
        state_changes = np.sum(seq[:-1] != seq[1:])
        sequence_features.append(state_changes)
        # Duration in States, Mean & Variance of Stay Durations
        state_durations = {state: [] for state in unique_states}
        current_state = seq[0]
        current_duration = 1
        for i in range(1, len(seq)):
            if seq[i] == current_state:
                current_duration += 1
            else:
                state_durations[current_state].append(current_duration)
                current_state = seq[i]
                current_duration = 1
        state_durations[current_state].append(current_duration)  # for the last state
        for _, durations in state_durations.items():
            if durations:
                mean_duration = np.mean(durations)
                var_duration = np.var(durations)
            else:
                mean_duration = 0
                var_duration = 0
            sequence_features.extend([mean_duration, var_duration])
        # Transition Probabilities
        for i in range(len(unique_states)):
            for j in range(len(unique_states)):
                total_transitions_from_i = np.sum(transitions_matrix[i, :])
                if total_transitions_from_i > 0:
                    transition_prob = (
                        transitions_matrix[i, j] / total_transitions_from_i
                    )
                else:
                    transition_prob = 0
                sequence_features.append(transition_prob)
        # Entropy of State Distribution
        state_counts = np.array(
            [np.sum(seq == state) for state in unique_states], dtype=float
        )
        state_probs = (
            state_counts / state_counts.sum()
            if state_counts.sum() > 0
            else np.zeros_like(state_counts)
        )
        sequence_entropy = entropy(
            state_probs
        )  # scipy's entropy function calculates from probabilities
        sequence_features.append(sequence_entropy)
        features.append(sequence_features)
    return np.array(features), feature_names


def save_HMM_model_to_bucket(directory, model, short_sha, bucket, model_path):
    file_name = directory + "/hmm_model_" + short_sha + ".joblib"
    dump(model, file_name)
    upload_blob(storage_client, bucket, file_name, model_path)
    return None


def main():

    dic_data = {key: {} for key in dataset_types}

    logging.info(f"Dataset types: {dataset_types}")
    logging.info(f"Samples dictionary: {samples_dic}")

    for dataset_name in dataset_types:
        logging.info(f"Importing dataset {dataset_name}")
        quoted_samples = ",".join(
            [f"'{sample}'" for sample in samples_dic[dataset_name]]
        )
        query = f"SELECT * FROM {project}.{ml_dataset}.features_wo_hmm WHERE sample IN ({quoted_samples}) AND {hmm_var} IS NOT NULL"
        # Execute the query and store in dic
        dic_data[dataset_name]["imported"] = bq_client.query(query).to_dataframe()
        dic_data[dataset_name]["imported"][hmm_var + "_no_string"] = dic_data[
            dataset_name
        ]["imported"][hmm_var].apply(lambda x: ast.literal_eval(x.strip('"')))

    logging.info("Creating a unique sequence for training the HMM")
    training_seq = np.concatenate(
        dic_data["TRAINING"]["imported"][hmm_var + "_no_string"].tolist()
    )

    logging.info("Reshaping training data")
    reshaped_data, lengths = prepare_data_for_hmm(training_seq)

    best_ll = None
    best_model = None
    for i in range(n_model_loop):
        print(f"Iteration: {i}")
        h = model_dic[model_type](
            n_states,
            n_iter=n_iterations,
            covariance_type=covariance,
            algorithm=algorithm,
            random_state=rs,
        )  # tol=1e-4,
        h.fit(reshaped_data, lengths)
        score = h.score(reshaped_data)
        print(f"score: {score}")
        if not best_ll or best_ll < best_ll:
            best_ll = score
            best_model = h

    logging.info("Saving model in the bucket")
    save_HMM_model_to_bucket(home_directory, best_model, short_sha, bucket, model_path)

    for dataset_name in dataset_types:
        df_imported = dic_data[dataset_name]["imported"].copy(deep=True)
        logging.info(f"Computing hidden states for dataset: {dataset_name}")
        dic_data[dataset_name]["hidden_states"] = predict_hidden_states_for_sequences(
            best_model, df_imported[hmm_var + "_no_string"]
        )
        logging.info(
            f"Compiling features based on the hidden state for dataset {dataset_name}"
        )
        features, feature_names = extract_features(
            dic_data[dataset_name]["hidden_states"]
        )
        # Use feature_names directly for column naming in the DataFrame
        hs_features_df = pd.DataFrame(features, columns=feature_names).reset_index(
            drop=True
        )
        # Round float values
        hs_features_df = np.round(hs_features_df.astype(float), 4)
        # Form final dataframe
        df_export = pd.concat(
            [df_imported, hs_features_df],
            axis=1,
        )
        # Drop columns for upload
        df_export.drop(hmm_var + "_no_string", axis=1, inplace=True)
        logging.info("Exportind dataset to BQ and Bucket")
        upload_dataframe_to_bq(bq_client, df_export, f"{ml_dataset}.{dataset_name}")
    logging.info("END OF SCRIPT")


# Start script
if __name__ == "__main__":
    try:
        main()
    except Exception as err:
        message = f"Script failed: {str(err)}"

        print(json.dumps({"message": message, "severity": "ERROR"}))
        sys.exit(1)  # Retry Job Task by exiting the process
