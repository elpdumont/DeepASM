#!/bin/bash

# you need this package to parse the YAML file
# brew install yq


SHORT_SHA="$(git rev-parse --short HEAD)"
echo "SHORT_SHA: ${SHORT_SHA}"

# Submit the build to Google Cloud Build
gcloud builds submit --config=cloudbuild.yaml . --substitutions=SHORT_SHA="${SHORT_SHA}"

export PROJECT_ID=$(yq e '.GCP.PROJECT_ID' config/config.yaml)
export REGION=$(yq e '.GCP.REGION' config/config.yaml)
export PYTHON_IMAGE=$(yq e '.GCP.IMAGE' config/config.yaml)
export GENOMIC_LENGTH=$(yq e '.GENOMICS.GENOMIC_LENGTH' config/config.yaml)
export BUCKET_NAME=$(yq e '.GCP.BUCKET_NAME' config/config.yaml)
export REFERENCE_GENOME=$(yq e '.GENOMICS.REFERENCE_GENOME' config/config.yaml)

export CLOUDASM_DATASET="${REFERENCE_GENOME}_${GENOMIC_LENGTH}_cloudasm"
export ML_DATASET="${REFERENCE_GENOME}_${GENOMIC_LENGTH}_ml_test2"

echo "CloudASM dataset: ${CLOUDASM_DATASET}"
echo "ML dataset: ${ML_DATASET}"

# Update the jobs file
# Copy the template to a new file that can be safely modified
mkdir jobs
cp jobs_templates/* jobs/

# Replace placeholders with actual values
sed -i '' "s#PYTHON_IMAGE_PLACEHOLDER#${PYTHON_IMAGE}#g" jobs/process_json.json
sed -i '' "s/IMAGE_TAG_PLACEHOLDER/${SHORT_SHA}/g" jobs/process_json.json
sed -i '' "s/ML_DATASET_ID_PLACEHOLDER/${ML_DATASET}/g" jobs/process_json.json
sed -i '' "s/CLOUDASM_DATASET_ID_PLACEHOLDER/${CLOUDASM_DATASET}/g" jobs/process_json.json



# Create respective folders in BigQuery and Cloud Storage if they do not exist

if bq ls --project_id="${PROJECT_ID}" | grep -w "${ML_DATASET}"; then
	echo "Dataset ${ML_DATASET} already exists in project ${PROJECT_ID}."
else
	# Create the dataset since it does not exist
	bq mk --project_id="${PROJECT_ID}" --dataset "${PROJECT_ID}:${ML_DATASET}"
	echo "Dataset ${ML_DATASET} created in project ${PROJECT_ID}."
fi
