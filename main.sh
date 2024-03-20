#!/bin/bash


#---------------------------------------------------------------


# List of table names to be used for ML
TABLE_NAMES=("tabular" "sequence_cpg_fm" "sequence_cpg_cov_and_methyl")

#---------------------------------------------------------------
# Export CloudASM output to bucket into JSON shards

gsutil -m rm gs://"${BUCKET_NAME}"/"${CLOUDASM_DATASET}"/*.json

bq extract --destination_format=NEWLINE_DELIMITED_JSON \
  --compression=NONE \
  "${PROJECT_ID}":"${CLOUDASM_DATASET}".hg_19_250_all_samples \
  gs://"${BUCKET_NAME}"/"${CLOUDASM_DATASET}"/cloudasm-*.json



# List, filter, and count JSON files
NUM_JSON_FILES=$(gsutil ls gs://"${BUCKET_NAME}"/"${CLOUDASM_DATASET}"/*.json | wc -l)

# Echo the count for demonstration
echo "Number of JSON files: ${NUM_JSON_FILES}"

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

gcloud batch jobs submit process-json \
	--location "${REGION}" \
	--config jobs/process_json.json

# To delete the job
gcloud batch jobs delete process-json --location "${REGION}"

# gcloud run jobs deploy process-json \
#  --image "${IMAGE}":"${IMAGE_TAG}" \
#  --args="python,/app/process_raw_json_to_bq.py" \
#  --set-env-vars ML_DATASET_ID="${ML_DATASET_ID}" \
#  --tasks 1 \
#  --max-retries 0 \
#  --cpu 8 \
#  --memory 32Gi \
#  --task-timeout 2000 \
#  --region "${REGION}" \
#  --project "${PROJECT_ID}" \
#  --execute-now

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
