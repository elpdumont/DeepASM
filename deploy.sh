#!/bin/bash

# you need this package to parse the YAML file
# brew install yq


SHORT_SHA="$(git rev-parse --short HEAD)"
echo "SHORT_SHA: ${SHORT_SHA}"

# Submit the build to Google Cloud Build
gcloud builds submit --config=cloudbuild.yaml . --substitutions=SHORT_SHA="${SHORT_SHA}"

export PROJECT_ID=$(yq e '.GCP.PROJECT_ID' config/config.yaml)
export PYTHON_IMAGE=$(yq e '.GCP.IMAGE' config/config.yaml)
export GENOMIC_LENGTH=$(yq e '.GENOMICS.GENOMIC_LENGTH' config/config.yaml)
export BUCKET_NAME=$(yq e '.GCP.BUCKET_NAME' config/config.yaml)
export REFERENCE_GENOME=$(yq e '.GENOMICS.REFERENCE_GENOME' config/config.yaml)

export CLOUDASM_DATASET="${REFERENCE_GENOME}_${GENOMIC_LENGTH}_cloudasm"
export ML_DATASET="${REFERENCE_GENOME}_${GENOMIC_LENGTH}_ml"

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

