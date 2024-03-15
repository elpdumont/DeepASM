#!/bin/bash


#---------------------------------------------------------------
# Prepare 3 datasets
# tabular
# methylation-sequence (nucleotide x 2 variables: Cpg presence, methylation presence)
# methylation-matrix (nucleotide x depth x 2 variables: Cpg presence, methylation presence)

gcloud run jobs deploy process-json \
 --image us-east1-docker.pkg.dev/hmh-em-deepasm/docker-repo/python:23d6622 \
 --args="python,/app/process_raw_json_to_bq.py" \
 --tasks 2 \
 --max-retries 1 \
 --cpu 4 \
 --memory 16Gi \
 --task-timeout 2000 \
 --region us-east1 \
 --project=hmh-em-deepasm

gcloud run jobs execute process-json

#------------------------------------------------------------------
# Load tabular data in BigQuery

gcloud run jobs deploy load-json-to-bq \
 --image us-east1-docker.pkg.dev/hmh-em-deepasm/docker-repo/python:latest \
 --args="python,/app/process_json_to_bq.py" \
 --tasks 1 \
 --set-env-vars BUCKET_NAME="hmh_deepasm" \
 --set-env-vars DATASET_TYPE="tabular" \
 --set-env-vars BUCKET_FOLDER_PATH="ml_datasets/" \
 --set-env-vars BQ_ML_DATASET_NAME="ml" \
 --set-env-vars BQ_ML_TABLE_NAME="tabular" \
 --max-retries 0 \
 --region us-east1 \
 --task-timeout 2000 \
 --project=hmh-em-deepasm

gcloud run jobs execute load-json-to-bq

