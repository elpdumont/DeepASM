#!/bin/bash


# Image TAG
IMAGE_TAG="4b10ec2"

# Path to the YAML file
CONFIG_FILE="config/config.yaml"

# Use `yq` to read the values. Ensure `yq` is installed and in your PATH.
PROJECT_ID=$(yq e ".GCP.PROJECT_ID" "${CONFIG_FILE}")
IMAGE=$(yq e ".GCP.IMAGE" "${CONFIG_FILE}")

GENOMIC_LENGTH=$(yq e '.GENOMICS.GENOMIC_LENGTH' "${CONFIG_FILE}")
MIN_CPG_COV=$(yq e '.GENOMICS.MIN_CPG_COV' "${CONFIG_FILE}")

# Use variables to create custom variables
ML_DATASET_ID="ml_${GENOMIC_LENGTH}_${MIN_CPG_COV}"

#---------------------------------------------------------------

# Check if the dataset exists
if bq ls --project_id="${PROJECT_ID}" | grep -w "${ML_DATASET_ID}"; then
    echo "Dataset ${ML_DATASET_ID} already exists in project ${PROJECT_ID}."
else
    # Create the dataset since it does not exist
    bq mk --project_id="${PROJECT_ID}" --dataset "${PROJECT_ID}:${ML_DATASET_ID}"
    echo "Dataset ${ML_DATASET_ID} created in project ${PROJECT_ID}."
fi

# List of table names to be used for ML
TABLE_NAMES=("tabular" "sequence_cpg_fm" "sequence_cpg_cov_and_methyl")



#---------------------------------------------------------------
# Prepare 3 datasets
# tabular
# methylation-sequence (nucleotide x 2 variables: Cpg presence, methylation presence)
# methylation-matrix (nucleotide x depth x 2 variables: Cpg presence, methylation presence)
#!/bin/bash

# Delete the tables if they exist
for TABLE_NAME in "${TABLE_NAMES[@]}"; do
    bq rm -f -t "${PROJECT_ID}:${ML_DATASET_ID}.${TABLE_NAME}"
done


gcloud run jobs deploy process-json \
 --image "${IMAGE}":"${IMAGE_TAG}" \
 --args="python,/app/process_raw_json_to_bq.py" \
 --set-env-vars ML_DATASET_ID="${ML_DATASET_ID}" \
 --tasks 1 \
 --max-retries 0 \
 --cpu 4 \
 --memory 16Gi \
 --task-timeout 2000 \
 --region us-east1 \
 --project="${PROJECT_ID}" \
 --execute-now

#------------------------------------------------------------------
# Load tabular data in BigQuery

# gcloud run jobs deploy load-json-to-bq \
#  --image us-east1-docker.pkg.dev/hmh-em-deepasm/docker-repo/python:latest \
#  --args="python,/app/process_json_to_bq.py" \
#  --tasks 1 \
#  --set-env-vars BUCKET_NAME="hmh_deepasm" \
#  --set-env-vars DATASET_TYPE="tabular" \
#  --set-env-vars BUCKET_FOLDER_PATH="ml_datasets/" \
#  --set-env-vars BQ_ML_DATASET_NAME="ml" \
#  --set-env-vars BQ_ML_TABLE_NAME="tabular" \
#  --max-retries 0 \
#  --region us-east1 \
#  --task-timeout 2000 \
#  --project=hmh-em-deepasm

# gcloud run jobs execute load-json-to-bq

