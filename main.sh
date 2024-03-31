#!/bin/bash


#---------------------------------------------------------------


# List of table names to be used for ML
TABLE_NAMES=("tabular" "sequence_cpg_fm" "sequence_cpg_fm_nonzeros" "sequence_cpg_cov_and_methyl" "sequence_cpg_cov_and_methyl_nonzeros")

#---------------------------------------------------------------
# Export CloudASM output to bucket into JSON shards

gsutil -m rm gs://"${BUCKET_NAME}"/"${CLOUDASM_DATASET}"/*.json

bq extract --destination_format=NEWLINE_DELIMITED_JSON \
  --compression=NONE \
  "${PROJECT_ID}":"${CLOUDASM_DATASET}".hg_19_250_all_samples \
  gs://"${BUCKET_NAME}"/"${CLOUDASM_DATASET}"/cloudasm-*.json



# List, filter, and count JSON files
NUM_JSON_FILES=$(gsutil ls gs://"${BUCKET_NAME}"/"${CLOUDASM_DATASET}"/*.json | wc -l)


#sed -i '' "s/TASK_COUNT_PLACEHOLDER/1/g" jobs/process_json.json
echo "Number of JSON files: ${NUM_JSON_FILES}"
sed -i '' "s/NB_FILES_PER_TASK_PLACEHOLDER/5/g" jobs/process_json.json
sed -i '' "s/TASK_COUNT_PLACEHOLDER/220/g" jobs/process_json.json
#sed -i '' "s/TASK_COUNT_PLACEHOLDER/1/g" jobs/process_json.json
# Echo the count for demonstration

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

# Delete the JSON files on the bucket
gsutil -m rm gs://"${BUCKET_NAME}"/"${ML_DATASET}"/*

JOB_NAME="process-json-${SHORT_SHA}"

# "${JOB_NAME}"

gcloud batch jobs submit "${JOB_NAME}"-2 \
	--location "${REGION}" \
	--config jobs/process_json.json



#---- For debugging
# To obtain failed jobs
gcloud batch tasks list --job="${JOB_NAME}" --location "${REGION}" --filter="STATE=FAILED"

gcloud batch tasks describe 69 \
  --location="${REGION}" \
  --job=process23 \
  --task_group=group0


# Count the number of regions excluded from the analysis.
bq show --format=prettyjson ${PROJECT_ID}:${ML_DATASET}.tabular | jq '.numRows'
bq show --format=prettyjson ${PROJECT_ID}:${CLOUDASM_DATASET}.hg_19_250_all_samples | jq '.numRows'



# Calculate the percentage difference
percentageDifference=$(echo "scale=2; (100 * ($numRowsML - $numRowsCloudASM) / (($numRowsML + $numRowsCloudASM) / 2))" | bc -l)

# Absolute value to handle negative values
percentageDifference=$(echo $percentageDifference | tr -d '-')

echo "The percentage difference in the number of rows is ${percentageDifference}%."