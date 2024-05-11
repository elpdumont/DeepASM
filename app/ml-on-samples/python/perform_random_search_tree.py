#!/usr/bin/env python

# import ndjson
import json
import logging
import os
import random
import sys

import numpy as np
import pandas as pd

# from pandas.io.json import json_normalize
import yaml
from gcp import upload_dataframe_to_bq
from google.cloud import bigquery
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from xgboost import XGBClassifier

# Create a handler for Google Cloud Logging.
logging.basicConfig(level=logging.INFO)

with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

# Accessing GCP configuration
project = config["GCP"]["PROJECT"]
# samples_dic = config["SAMPLES"]
# all_samples = [item for sublist in samples_dic.values() for item in sublist]

# ML variables
n_random_search = config["ML"]["N_RANDOM_SEARCH"]

# Retrieve Job-defined env vars
BATCH_TASK_INDEX = int(os.getenv("BATCH_TASK_INDEX", 0))
ml_dataset = os.getenv("ML_DATASET")

# Define the path to the JSON credentials file
credentials_path = "/appuser/.config/gcloud/application_default_credentials.json"

# Check if the JSON credentials file exists
if os.path.exists(credentials_path):
    # Set the GOOGLE_APPLICATION_CREDENTIALS environment variable to the path of the JSON key file
    os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = credentials_path
    # Assuming 'project' is already defined somewhere in your script
    os.environ["GOOGLE_CLOUD_PROJECT"] = project
    samples_dataset = "samples_250bp"

# Initialize the Google Cloud Storage client
bq_client = bigquery.Client(project=project)

dic_model = {
    0: {
        "model": RandomForestClassifier,
        "weight_name": "class_weight",
        "grid": {
            "max_depth": list(range(10, 110, 10)) + [None],
            "n_estimators": list(range(10, 160, 10)),
            "min_samples_split": [2, 5, 10],
            "min_samples_leaf": [1, 2, 4],
            "bootstrap": [True, False],
        },
    },
    1: {
        "model": XGBClassifier,
        "weight_name": "scale_pos_weight",
        "grid": {
            "max_depth": list(range(5, 80, 5)) + [None],
            "n_estimators": list(range(10, 160, 10)),
            "learning_rate": [0.001, 0.01, 0.1],
            "min_child_weight": [1, 5, 10],
            "subsample": [0.2, 0.4, 0.6, 0.8, 1.0],
            "colsample_bytree": [0.6, 0.8, 1.0],
            "reg_alpha": [0, 0.001, 0.01, 0.1],
        },
    },
}

datasets = ["TRAINING", "VALIDATION", "TESTING", "TRAINING_AND_VALIDATION"]

# Set a random seed for reproducibility
random_seed = 42
random.seed(random_seed)


def compute_class_weight(dic_data):
    class_weight_dic = {"TRAINING": {}, "TRAINING_AND_VALIDATION": {}}
    for data in ["TRAINING", "TRAINING_AND_VALIDATION"]:
        labels = dic_data[data]["labels"]["asm"]
        weights = np.round(
            compute_class_weight("balanced", classes=[0, 1], y=labels), 2
        )
        class_weight_dic[data]["class_weight"] = weights
        scale_pos_weight = len(labels[labels == 0]) / len(labels[labels == 1])
        class_weight_dic[data]["scale_pos_weight"] = scale_pos_weight
    return class_weight_dic


def evaluate_model_for_trees(dic_data, dataset, model):
    data = dic_data[dataset]["TABULAR"]
    labels = dic_data[dataset]["labels"].squeeze()
    data = np.array(data)
    predictions = model.predict(data)
    confusion = confusion_matrix(labels, predictions)
    report = classification_report(labels, predictions, output_dict=True)
    sum_f1 = np.round(report["0.0"]["f1-score"] + report["1.0"]["f1-score"], 3)
    # report = classification_report(labels, predictions, digits=3)
    return sum_f1, confusion, report


def main():
    dic_data = {dataset: {} for dataset in datasets}

    for dataset in ["TRAINING", "VALIDATION", "TESTING"]:
        logging.info(f"Processing {dataset} dataset...")
        query = f"""
            SELECT * EXCEPT (region_sup, clustering_index, region_nb_cpg, cpg_directional_fm, cpgs_w_padding)
            FROM {project}.{ml_dataset}.{dataset}
            WHERE cpg_directional_fm IS NOT NULL AND asm IS NOT NULL
            LIMIT 10000
            """
        df = bq_client.query(query).to_dataframe()
        dic_data[dataset]["labels"] = df[["asm"]]
        dic_data[dataset]["region_info"] = df[
            ["asm", "sample", "chr", "region_inf", "nb_cpg_found", "nb_reads"]
        ]
        dic_data[dataset]["tabular"] = df.drop(
            columns=[
                "asm",
                "sample",
                "chr",
                "region_inf",
                "nb_cpg_found",
                "nb_reads",
            ]
        )

    for data in ["labels", "region_info", "tabular"]:  # '1d_seq', '2d_seq',
        dic_data["TRAINING_AND_VALIDATION"][data] = pd.concat(
            [dic_data["TRAINING"][data], dic_data["VALIDATION"][data]]
        )

    class_weight_dic = compute_class_weight(dic_data)

    model_params = dic_model[BATCH_TASK_INDEX]
    model = model_params["model"]
    grid = model_params["grid"]
    weight_name = model_params["weight_name"]

    weight_dic = {weight_name: class_weight_dic["TRAINING"][weight_name]}

    random_hyperparameters_list = []
    for _ in range(n_random_search):
        random_params = {key: random.choice(value) for key, value in grid.items()}
        random_hyperparameters_list.append(random_params)

    # Prepare training
    x_train = np.array(dic_data["TRAINING"]["tabular"])
    y_train = dic_data["TRAINING"]["labels"].squeeze()

    results = []
    for idx, params in enumerate(random_hyperparameters_list):
        logging.info(f"Dictionary {idx+1}: {params}")

        model = model(**weight_dic, **params)
        model.fit(x_train, y_train)

        # evaluate model on f_1 score
        sumf1, confusion, report = evaluate_model_for_trees(
            dic_data, "VALIDATION", model
        )
        logging.info(f"Sum of F1s: {sumf1}")

        # save results
        results.append({"parameters": params, "sumf1": sumf1})

    # get the best hyperparameters
    best_hyperparameters = max(results, key=lambda x: x["sumf1"])["parameters"]
    logging.info(best_hyperparameters)
    # del best_hyperparameters['sumf1']

    # Retrain model on training and validation
    x_train = np.array(dic_data["TRAINING_AND_VALIDATION"]["tabular"])
    y_train = dic_data["TRAINING_AND_VALIDATION"]["labels"].squeeze()
    weight_dic = {weight_name: class_weight_dic["TRAINING_AND_VALIDATION"][weight_name]}

    best_model = model(weight_dic, **best_hyperparameters)

    best_model.fit(x_train, y_train)
    sumf1, confusion, report = evaluate_model_for_trees(dic_data, "TESTING", best_model)
    logging.info(f"GENERALIZATION ERROR, {sumf1}, {confusion}, {report}")

    # Create a dictionary to upload to BQ
    confusion_dict = {
        "True Negative": confusion[0, 0],
        "False Positive": confusion[0, 1],
        "False Negative": confusion[1, 0],
        "True Positive": confusion[1, 1],
    }
    dic_results = {
        **{"model": str(model), "sum_f1": sumf1},
        **confusion_dict,
        **report,
        **best_hyperparameters,
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
