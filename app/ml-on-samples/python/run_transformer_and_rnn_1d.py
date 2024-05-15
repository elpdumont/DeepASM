#!/usr/bin/env python

# import ndjson
import ast
import json
import logging
import os
import random

# import re
import sys

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

# from pandas.io.json import json_normalize
import yaml
from gcp import upload_blob, upload_dataframe_to_bq
from google.cloud import bigquery, storage

# from joblib import dump
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.utils.class_weight import compute_class_weight
from torch.utils.data import DataLoader, TensorDataset

# import torch.nn.functional as F
# import torch.optim as optim


# Create a handler for Google Cloud Logging.
logging.basicConfig(level=logging.INFO)

with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

# Accessing GCP configuration
project = config["GCP"]["PROJECT"]
bucket = config["GCP"]["BUCKET"]

# samples_dic = config["SAMPLES"]
# all_samples = [item for sublist in samples_dic.values() for item in sublist]

# ML variables
label_var = config["ML"]["LABEL_NAME"]
n_random_search = config["ML"]["N_RANDOM_SEARCH_1D"]
batch_size = config["ML"]["BATCH_SIZE"]
dropout_rate = config["ML"]["DROPOUT_RATE"]
n_epochs = config["ML"]["N_EPOCHS"]
padding_value = config["ML"]["PADDING_VALUE"]

# Obtain sample list
samples_dic = config["SAMPLES"]
dataset_types = list(samples_dic.keys())

# Retrieve Job-defined env vars
BATCH_TASK_INDEX = int(os.getenv("BATCH_TASK_INDEX", 0))
model_path = os.getenv("MODEL_PATH")
ml_dataset = os.getenv("ML_DATASET")
short_sha = os.getenv("SHORT_SHA")
home_directory = os.path.expanduser("~")

# Define the path to the JSON credentials file
credentials_path = "/appuser/.config/gcloud/application_default_credentials.json"

# Check if the JSON credentials file exists
if os.path.exists(credentials_path):
    # Set the GOOGLE_APPLICATION_CREDENTIALS environment variable to the path of the JSON key file
    os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = credentials_path
    # Assuming 'project' is already defined somewhere in your script
    os.environ["GOOGLE_CLOUD_PROJECT"] = project
    samples_dataset = "samples_250bp"
    BATCH_TASK_INDEX = 0
    model_path = "samples_250bp/models"
    ml_dataset = "ml_250bp_3"
    short_sha = "test"

# Initialize the Google Cloud Storage client
bq_client = bigquery.Client(project=project)
storage_client = storage.Client(project=project)

device = torch.device(
    f"cuda:{BATCH_TASK_INDEX}" if torch.cuda.is_available() else "cpu"
)


class TransformerModel(nn.Module):
    def __init__(
        self,
        d_model,
        nhead,
        dim_feedforward,
        max_sequence_length,
        num_layers,
        dropout,
        device,
        **kwargs,
    ):
        self.d_model = d_model
        self.device = device
        # self.max_sequence_length = max_sequence_length
        super(TransformerModel, self).__init__()
        self.encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=dim_feedforward,
            dropout=dropout,  # Adding dropout to the encoder layer
            layer_norm_eps=1e-6,  # Using a smaller epsilon for layer normalization for more precise calculations
        )
        self.transformer_encoder = nn.TransformerEncoder(
            self.encoder_layer, num_layers=num_layers
        )
        self.position_embeddings = nn.Embedding(
            max_sequence_length, d_model
        )  # Prepare position embeddings
        self.fc = nn.Linear(d_model, 1)

    def forward(self, x):
        x = x.to(self.device)
        x = x.unsqueeze(-1).repeat(
            1, 1, self.d_model
        )  # Extend features to match d_model
        seq_length, N = x.shape[1], x.shape[0]
        positions = (
            torch.arange(0, seq_length, device=x.device).unsqueeze(0).repeat(N, 1)
        )
        x += self.position_embeddings(positions)
        x = x.permute(1, 0, 2)  # Reshape x to [seq_length, batch, features]
        x = self.transformer_encoder(x)
        x = self.fc(x[-1, :, :])  # Take the last sequence output
        # return torch.sigmoid(x.view(-1))  # Flatten the output for compatibility with target
        return x.view(-1)

    def predict(self, x, threshold=0.5):
        x = x.to(self.device)
        logits = self.forward(x)
        predictions = torch.sigmoid(
            logits
        )  # Apply sigmoid to convert logits to probabilities
        return (
            predictions > threshold
        ).float()  # Threshold the probabilities to get binary predictions


class RNNModel(nn.Module):
    def __init__(
        self,
        input_size,
        hidden_size,
        output_size,
        max_sequence_length,
        dropout_rate,
        device,
        **kwargs,
    ):
        self.device = device
        super(RNNModel, self).__init__()
        self.to(device)
        self.max_sequence_length = max_sequence_length
        self.rnn = nn.LSTM(
            input_size=input_size,
            hidden_size=hidden_size,
            batch_first=True,
            dropout=dropout_rate,
            device=device,
        )
        self.dropout = nn.Dropout(dropout_rate)
        self.fc = nn.Linear(hidden_size, output_size)

    def forward(self, x):
        x.to(self.device)
        # Check input length
        if x.shape[1] != self.max_sequence_length:
            raise ValueError(
                f"Expected sequence length {self.max_sequence_length}, but got {x.shape[1]}"
            )
        # LSTM output shape: (batch_size, sequence_length, hidden_size)
        x = x.unsqueeze(
            -1
        )  # Increase feature dimension from (batch, seq_len) to (batch, seq_len, features)
        output, _ = self.rnn(x)
        # Applying dropout to the output of the LSTM
        output = self.dropout(
            output[:, -1, :]
        )  # Applying dropout to the last output of the sequence
        # Pass the last output to the fully connected layer
        logits = self.fc(output)
        return logits.view(-1)  # Flatten the output to shape (batch,)

    def predict(self, x, threshold=0.5):
        x = x.to(self.device)
        logits = self.forward(x)
        predictions = torch.sigmoid(
            logits
        )  # Apply sigmoid to convert logits to probabilities
        return (
            predictions > threshold
        ).float()  # Threshold the probabilities to get binary predictions


def create_model(model_number, **kwargs):
    # logging.info(kwargs)
    if model_number == 0:
        return TransformerModel(**kwargs)
    elif model_number == 1:
        return RNNModel(**kwargs)
    else:
        raise ValueError("Unsupported model number")


def evaluate_model(model, dataloader, criterion, device):
    model = model.to(device)
    model.eval()  # Set the model to evaluation mode
    total_loss = 0
    all_preds = []
    all_targets = []
    with torch.no_grad():
        for data, targets in dataloader:
            data, targets = data.to(device), targets.to(device)
            # Forward pass
            outputs = model.forward(data)
            preds = model.predict(data)
            loss = criterion(outputs, targets)
            # Update lists
            all_preds.extend(preds.view(-1).cpu().numpy())
            all_targets.extend(targets.view(-1).cpu().numpy())
            # Update loss
            total_loss += loss.item() * data.size(0)
    # Calculate the average loss
    average_loss = total_loss / len(dataloader.dataset)
    logging.info(f"Average loss: {average_loss}")
    # Generate the classification report
    report = classification_report(all_targets, all_preds, output_dict=True, digits=4)
    confusion = confusion_matrix(all_targets, all_preds)
    # Return the report as a dictionary for further analysis if needed
    return report["0"]["f1-score"] + report["1"]["f1-score"], report, confusion


def train_seq_model(
    model, num_epochs, training_data_dataloader, optimizer, criterion, device
):
    for epoch in range(num_epochs):
        total_loss = 0
        num_batches = 0
        for data, targets in training_data_dataloader:
            data, targets = data.to(device), targets.to(device)
            # Training Transformer
            model.zero_grad()
            outputs = model(data)
            # Compute loss
            loss = criterion(outputs, targets)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
            num_batches += 1
        average_loss = total_loss / num_batches
        if (epoch + 1) % 10 == 0:
            logging.info(f"Epoch {epoch+1}, Average loss: {average_loss:.4f}")
    logging.info("End of model training")
    return None


dic_model = {
    0: {
        "model": "TransformerModel",
        "weight_name": "class_weight",
        "grid": {
            "d_model": [8, 16, 32],  # Dimensionality of the model
            "nhead": [2, 4, 8],  # Number of heads in the multi-head attention models
            "num_layers": [
                2,
                4,
                8,
                16,
                32,
            ],  # Number of sub-encoder-layers in the transformer
            "dim_feedforward": [
                16,
                32,
                64,
                128,
                256,
            ],  # Size of the feedforward model in nn.TransformerEncoder
            # Learning rate for the optimizer
            "dropout": [dropout_rate],
            "learning_rate": [
                0.0001,
                0.001,
                0.01,
                0.1,
            ],
            "num_epochs": [n_epochs],
            "weight_decay": [0.001],  # Number of training epochs
        },
    },
    1: {
        "model": "RNNModel",
        "weight_name": "class_weight",
        "grid": {
            "input_size": [1],
            "hidden_size": [8, 16, 32, 64, 128],
            "output_size": [1],
            "dropout_rate": [dropout_rate],
            "subsample": [0.2, 0.4, 0.6, 0.8, 1.0],
            "learning_rate": [
                0.0001,
                0.001,
                0.01,
                0.1,
            ],
            "num_epochs": [n_epochs],
            "weight_decay": [0.001],
        },
    },
}

datasets = ["TRAINING", "VALIDATION", "TESTING", "TRAINING_AND_VALIDATION"]

# Set a random seed for reproducibility
random_seed = 42
random.seed(random_seed)


def compute_classes(dic_data, device):
    class_weight_dic = {"TRAINING": {}, "TRAINING_AND_VALIDATION": {}}
    for data in ["TRAINING", "TRAINING_AND_VALIDATION"]:
        labels = dic_data[data]["labels"]
        weights = np.round(
            compute_class_weight("balanced", classes=[0, 1], y=labels), 2
        )
        class_weight_dic[data]["class_weight"] = {0: weights[0], 1: weights[1]}
        scale_pos_weight = len(labels[labels == 0]) / len(labels[labels == 1])
        class_weight_dic[data]["scale_pos_weight"] = scale_pos_weight
        class_weight_dic[data]["weight_tensor"] = torch.tensor(
            weights, dtype=torch.float32
        ).to(device)
    return class_weight_dic


def save_1d_model_to_bucket(
    directory, model_name, model, short_sha, bucket, model_path, storage_client
):
    file_name = directory + "/" + model_name + "_" + short_sha + ".pth"
    # Save model locally
    torch.save(model.state_dict(), file_name)
    # Save model in bucket
    upload_blob(storage_client, bucket, file_name, model_path)
    return None


def main():
    logging.info(f"Using device: {device}")
    dic_data = {dataset: {} for dataset in dataset_types}
    max_sequence_length = 0
    for dataset in dataset_types:
        logging.info(f"Processing {dataset} dataset...")
        quoted_samples = ",".join([f"'{sample}'" for sample in samples_dic[dataset]])
        logging.info(f"Importing the samples: {quoted_samples}")
        query = f"""
            SELECT * EXCEPT (region_sup, clustering_index, region_nb_cpg, cpgs_w_padding)
            FROM {project}.{ml_dataset}.features_wo_hmm
            WHERE
                cpg_directional_fm IS NOT NULL AND
                {label_var} IS NOT NULL AND
                sample IN ({quoted_samples})
            LIMIT 10000
            """
        df = bq_client.query(query).to_dataframe()
        dic_data[dataset]["labels"] = df[label_var].astype(int)
        dic_data[dataset]["region_info"] = df[
            ["asm", "sample", "chr", "region_inf", "nb_cpg_found", "nb_reads"]
        ]
        dic_data[dataset]["1d_seq"] = df["cpg_directional_fm"].apply(
            lambda x: ast.literal_eval(x.strip('"'))
        )
        current_max_sequence_length = max(dic_data[dataset]["1d_seq"].apply(len))
        logging.info(
            f"Max sequence length in the dataset: {current_max_sequence_length}"
        )
        if current_max_sequence_length > max_sequence_length:
            max_sequence_length = current_max_sequence_length

    for data in ["labels", "region_info", "1d_seq"]:  # '1d_seq', '2d_seq',
        dic_data["TRAINING_AND_VALIDATION"][data] = pd.concat(
            [dic_data["TRAINING"][data], dic_data["VALIDATION"][data]]
        )

    logging.info("Adding a dummy row of variables to adjust the max sequence length")
    zero_sequence = [0] * max_sequence_length  # This is the sequence to append
    for dataset in dic_data.keys():
        # Append a new sequence of zeros to the '1d_seq' array
        dic_data[dataset]["1d_seq"] = pd.concat(
            [dic_data[dataset]["1d_seq"], pd.Series([zero_sequence])], ignore_index=True
        )
        # Assuming you need to add a corresponding new label
        dic_data[dataset]["labels"] = pd.concat(
            [dic_data[dataset]["labels"], pd.Series([0])], ignore_index=True
        )

    for dataset in datasets:
        logging.info(f"Processing Data Loader for: {dataset}")
        sequences = [list(row) for row in dic_data[dataset]["1d_seq"]]
        labels = list(dic_data[dataset]["labels"])
        # Convert sequences to tensors and pad them
        padded_sequences = torch.nn.utils.rnn.pad_sequence(
            [torch.tensor(s) for s in sequences],
            batch_first=True,
            padding_value=padding_value,
        )
        # Convert labels to a tensor
        labels = torch.tensor(labels, dtype=torch.float32)
        # Create DataLoader for batch processing
        data = TensorDataset(padded_sequences, labels)
        dic_data[dataset]["dataloader"] = DataLoader(
            data, batch_size=batch_size, shuffle=True
        )

    logging.info("Computing class weights")
    class_weight_dic = compute_classes(dic_data, device)

    logging.info("Defining criterion")
    criterion = nn.BCEWithLogitsLoss(
        pos_weight=class_weight_dic["TRAINING"]["weight_tensor"][1]
    )

    model_info = dic_model[BATCH_TASK_INDEX]
    # model = model_params["model"]
    model_name = model_info["model"]
    # model_name = re.search(r"\.(\w+)[^\.]*>$", model_name).group(1)
    model_grid = model_info["grid"]

    random_hyperparameters_list = []
    for _ in range(n_random_search):
        random_model_params = {
            key: random.choice(value) for key, value in model_grid.items()
        }
        random_hyperparameters_list.append(
            {
                **random_model_params,
                **{"max_sequence_length": max_sequence_length, "device": device},
            }
        )

    results = []
    for idx, params in enumerate(random_hyperparameters_list):
        logging.info(f"Dictionary {idx+1}: {params} for model {model_name}")
        model_w_params = create_model(BATCH_TASK_INDEX, **params)
        model_w_params = model_w_params.to(device)
        logging.info("Defining optimizer")
        optimizer = torch.optim.AdamW(
            model_w_params.parameters(),
            lr=params["learning_rate"],
            weight_decay=params["weight_decay"],
        )
        logging.info("Start training...")
        train_seq_model(
            model_w_params,
            params["num_epochs"],
            dic_data["TRAINING"]["dataloader"],
            optimizer,
            criterion,
            device,
        )
        logging.info("Evaluate model on f_1 score")
        sumf1, report, confusion = evaluate_model(
            model_w_params, dic_data["VALIDATION"]["dataloader"], criterion, device
        )
        logging.info(f"Sum of F1s: {sumf1}")
        # save results
        results.append({"parameters": params, "sumf1": sumf1})

    # get the best hyperparameters
    best_hyperparameters = max(results, key=lambda x: x["sumf1"])["parameters"]
    logging.info(f"Best hyperparameters: {best_hyperparameters} for model {model_name}")
    # del best_hyperparameters['sumf1']

    # Retrain model on training and validation
    criterion = nn.BCEWithLogitsLoss(
        pos_weight=class_weight_dic["TRAINING_AND_VALIDATION"]["weight_tensor"][1]
    )

    best_model = create_model(BATCH_TASK_INDEX, **best_hyperparameters)
    best_model = best_model.to(device)
    optimizer = torch.optim.AdamW(
        best_model.parameters(),
        lr=best_hyperparameters["learning_rate"],
        weight_decay=best_hyperparameters["weight_decay"],
    )

    train_seq_model(
        best_model,
        n_epochs,
        dic_data["TRAINING_AND_VALIDATION"]["dataloader"],
        optimizer,
        criterion,
        device,
    )

    # Save model
    save_1d_model_to_bucket(
        home_directory,
        model_name,
        best_model,
        short_sha,
        bucket,
        model_path,
        storage_client,
    )

    sumf1, report, confusion = evaluate_model(
        best_model, dic_data["TESTING"]["dataloader"], criterion, device
    )
    logging.info(f"GENERALIZATION ERROR, {sumf1}, {confusion}, {report}")

    # Create a dictionary to upload to BQ
    confusion_dict = {
        "True Negative": confusion[0, 0],
        "False Positive": confusion[0, 1],
        "False Negative": confusion[1, 0],
        "True Positive": confusion[1, 1],
    }

    # Change keys for bigquery
    report["class_0"], report["class_1"] = report["0"], report["1"]
    del report["0"], report["1"]

    dic_results = {
        **{
            "model": model_name,
            "short_sha": short_sha,
            "n_random_search": n_random_search,
            "sum_f1": sumf1,
            "hyper_parameters": str(best_hyperparameters),
        },
        **confusion_dict,
        **report,
    }

    results_df = pd.DataFrame([dic_results])  # The dictionary is put inside a list
    upload_dataframe_to_bq(bq_client, results_df, f"{ml_dataset}.model_results")


# Start script
if __name__ == "__main__":
    try:
        main()
    except Exception as err:
        message = f"Error: {str(err)}"

        print(json.dumps({"message": message, "severity": "ERROR"}))
        sys.exit(1)  # Retry Job Task by exiting the process
