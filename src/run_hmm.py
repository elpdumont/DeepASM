# Python packages for data, stats
import json
import logging
import os
import sys

import numpy as np
import pandas as pd
import yaml
from google.cloud import bigquery, storage
from hmmlearn.hmm import GaussianHMM
from joblib import dump
from scipy.stats import entropy
from sklearn.utils import check_random_state

from gcp_utils import (
    export_dataframe_to_gcs_as_json,
    upload_blob,
    upload_dataframe_to_bq,
)

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
project_id = config["GCP"]["PROJECT_ID"]
bucket_name = config["GCP"]["BUCKET_NAME"]
model_folder = config["GCP"]["MODEL_FOLDER"]

# Genomic variables
genomic_length = config["GENOMICS"]["GENOMIC_LENGTH"]

# Samples for splitting the dataset
samples_dic = config["SAMPLES"]
dataset_types = list(samples_dic.keys())

# HMM variables
n_states = config["HMM"]["N_STATES"]
covariance = config["HMM"]["COVARIANCE"]
n_iterations = config["HMM"]["N_ITERATIONS"]
algorithm = config["HMM"]["ALGORITHM"]
hhm_var = config["HMM"]["VAR_NAME"]

# Retrieve Job-defined env vars
BATCH_TASK_INDEX = int(os.getenv("BATCH_TASK_INDEX", 0))
ml_dataset_id = os.getenv("ML_DATASET_ID")


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
            logging(f"Processed sequence number: {i + 1}")

    return predicted_states


def extract_features(hidden_states_sequences):
    """
    Extracts a comprehensive set of features from sequences of hidden states in a Hidden Markov Model (HMM).

    This function processes each sequence of hidden states to extract various features such as counts, proportions, transitions between states, start and end states, the number of state changes, durations in each state, and transition probabilities. It also calculates the mean and variance of the durations spent in each state and the entropy of the state distribution across each sequence.

    Parameters:
    - hidden_states_sequences (list of list of int): A list of sequences, where each sequence is a list of integers representing the hidden states visited by the HMM.

    Returns:
    - np.array: A 2D NumPy array where each row represents the extracted features for a single sequence. The features include counts and proportions of each state, transitions between states, the start and end state of each sequence, the number of state changes within the sequence, the mean and variance of durations spent in each state, transition probabilities between each pair of states, and the entropy of the state distribution.
    """

    features = []

    # Determine all unique states across sequences for consistent ordering
    unique_states = sorted({state for seq in hidden_states_sequences for state in seq})

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

        for state, durations in state_durations.items():
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

    return np.array(features)


def save_HMM_model_to_bucket(model):
    file_name = "hmm_model.joblib"
    dump(model, file_name)
    upload_blob(bucket_name, model_folder, file_name)


def main():

    dic_data = {}

    for dataset_name in dataset_types:
        logging.info(f"Importing dataset {dataset_name}")
        quoted_samples = ",".join(
            [f"'{sample}'" for sample in samples_dic[dataset_name]]
        )

        query = f"SELECT * FROM {project_id}.{ml_dataset_id}.tabular WHERE sample IN ({quoted_samples}) LIMI 100"

        dic_data[dataset_name]["imported"] = bq_client.query(query).to_dataframe()

    logging.info("Creating a unique sequence for training the HMM")
    training_seq = np.array(
        np.concatenate(dic_data["train"]["imported"][hhm_var].values)
    )

    reshaped_data, lengths = prepare_data_for_hmm(training_seq)

    model = GaussianHMM(
        n_components=n_states,
        covariance_type=covariance,
        n_iter=n_iterations,
        algorithm=algorithm,
        random_state=rs,
    )

    logging.info("Fitting the model")
    model.fit(reshaped_data, lengths)

    logging.info("Saving model in the bucket")
    save_HMM_model_to_bucket(model)

    for dataset_type in dataset_types:
        logging(f"Computing hidden states for dataset: {dataset_type}")
        dic_data[dataset_type]["hidden_states"] = predict_hidden_states_for_sequences(
            model, dic_data[dataset_type]["imported"][hhm_var]
        )

        logging.info(
            f"Compiling features based on the hidden state for dataset {dataset_type}"
        )

        dic_data[dataset_type]["hs_features"] = extract_features(
            dic_data[dataset_type]["hidden_states"]
        )

        # Determine the number of features dynamically from the shape of the hs_features array
        num_features = dic_data[dataset_type]["hs_features"].shape[1]

        # Convert the N x num_features NumPy array to a DataFrame with dynamic column naming
        hs_features_df = pd.DataFrame(
            dic_data[dataset_type]["hs_features"],
            columns=[f"hs_{i+1}" for i in range(num_features)],
        )

        # Create dataframe to export to BigQuery
        dic_data[dataset_type]["to_export"] = pd.concat(
            [dic_data[dataset_type]["imported"], hs_features_df], axis=1
        )

    logging.info("Exporting dataset to BigQuery")
    for dataset_type in dataset_types:

        df = dic_data[dataset_type]["to_export"]
        upload_dataframe_to_bq(bq_client, df, f"{ml_dataset_id}.{dataset_type}")
        export_dataframe_to_gcs_as_json(
            df,
            bucket_name,
            ml_dataset_id,
            BATCH_TASK_INDEX,
            dataset_type,
        )

    logging.info("END OF SCRIPT")


# Start script
if __name__ == "__main__":
    try:
        main()
    except Exception as err:
        message = f"Task #{BATCH_TASK_INDEX} failed: {str(err)}"

        print(json.dumps({"message": message, "severity": "ERROR"}))
        sys.exit(1)  # Retry Job Task by exiting the process
