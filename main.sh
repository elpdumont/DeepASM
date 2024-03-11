#!/bin/bash


#---------------------------------------------------------------
# Prepare 3 datasets
# tabular
# methylation-sequence (nucleotide x 2 variables: Cpg presence, methylation presence)
# methylation-matrix (nucleotide x depth x 2 variables: Cpg presence, methylation presence)

gcloud run jobs deploy process-json \
 --image us-east1-docker.pkg.dev/hmh-em-deepasm/docker-repo/python:latest \
 --args="python,/app/process_json.py" \
 --tasks 116 \
 --parallelism 10 \
 --set-env-vars BUCKET_NAME="hmh_deepasm" \
 --set-env-vars INPUT_DATA_FOLDER_PATH="bq_tables/250bp_asm_labelled/" \
 --set-env-vars OUTPUT_DATA_FOLDER_PATH="ml_datasets/" \
 --set-env-vars GENOMIC_LENGTH="250" \
 --set-env-vars MIN_CPG_COV="20" \
 --set-env-vars KERNEL_FM_NB_VALUES="10" \
 --set-env-vars KERNEL_FM_BANDWIDTH="0.1" \
 --set-env-vars KERNEL_COV_NB_MAX="200" \
 --set-env-vars KERNEL_COV_NB_STEP="10" \
 --set-env-vars KERNEL_COV_BANDWIDTH="5" \
 --max-retries 3 \
 --cpu 4 \
 --memory 16Gi \
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
 --max-retries 3 \
 --region us-east1 \
 --project=hmh-em-deepasm

gcloud run jobs execute load-json-to-bq

