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


sed -i '' "s/TASK_COUNT_PLACEHOLDER/${NUM_JSON_FILES}/g" jobs/process_json.json
sed -i '' "s/NB_FILES_PER_TASK_PLACEHOLDER/5/g" jobs/process_json.json

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
	bq rm -f -t "${PROJECT_ID}:${ML_DATASET}.${TABLE_NAME}"
done

JOB_NAME="process-json-${SHORT_SHA}"

gcloud batch jobs submit "${JOB_NAME}" \
	--location "${REGION}" \
	--config jobs/process_json.json