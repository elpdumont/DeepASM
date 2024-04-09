#!/bin/bash

# Import environmental variables
source src/import_env_variables.sh

echo "Preparing standard genomic regions using the reference genome"
chmod +x src/prepare_regions_using_ref_genome.sh
src/prepare_regions_using_ref_genome.sh


# Update the jobs file
# Copy the template to a new file that can be safely modified
mkdir jobs
cp jobs_templates/* jobs/

# Replace placeholders with actual values
for file in jobs/process_json.json jobs/run_hmm.json; do
    sed -i '' "s#PYTHON_IMAGE_PLACEHOLDER#${PYTHON_IMAGE}#g" "${file}"
    sed -i '' "s/IMAGE_TAG_PLACEHOLDER/${SHORT_SHA}/g" "${file}"
    sed -i '' "s/ML_DATASET_ID_PLACEHOLDER/${ML_DATASET}/g" "${file}"
done

sed -i '' "s/CLOUDASM_DATASET_ID_PLACEHOLDER/${CLOUDASM_DATASET}/g" jobs/process_json.json



# Create respective folders in BigQuery and Cloud Storage if they do not exist

if bq ls --project_id="${PROJECT_ID}" | grep -w "${ML_DATASET}"; then
	echo "Dataset ${ML_DATASET} already exists in project ${PROJECT_ID}."
else
	# Create the dataset since it does not exist
	bq mk --project_id="${PROJECT_ID}" --dataset "${PROJECT_ID}:${ML_DATASET}"
	echo "Dataset ${ML_DATASET} created in project ${PROJECT_ID}."
fi

# List of table names to be used for ML
TABLE_NAMES=("tabular" "sequence_cpg_cov_and_methyl")


#---------------------------------------------------------------
# Prepare reference genome annotation


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


# Delete the tables if they exist
for TABLE_NAME in "${TABLE_NAMES[@]}"; do
	bq rm -f -t "${PROJECT_ID}:${ML_DATASET}.${TABLE_NAME}"
done

# Delete the JSON files on the bucket
gsutil -m rm gs://"${BUCKET_NAME}"/"${ML_DATASET}"/*.json

JOB_NAME="process-json-${SHORT_SHA}"

# "${JOB_NAME}"

gcloud batch jobs submit "${JOB_NAME}" \
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
# Fetch the number of rows in each table as strings but formatted as numeric values
numRowsML=$(bq show --format=prettyjson ${PROJECT_ID}:${ML_DATASET}.tabular | jq -r '.numRows')
numRowsCloudASM=$(bq show --format=prettyjson ${PROJECT_ID}:${CLOUDASM_DATASET}.hg_19_250_all_samples | jq -r '.numRows')

# The `-r` option with `jq` ensures the output is raw, making it suitable for numeric calculations in `bc`
percentageDifference=$(echo "scale=2; ($numRowsML - $numRowsCloudASM) / (($numRowsML + $numRowsCloudASM) / 2) * 100" | bc -l | tr -d '-')

echo "The percentage difference in the number of rows is ${percentageDifference}%."



#---------------------------------------------------------------
# Run HMM and split the dataset into train, validation, and test.

DATASET_NAMES=("TRAINING" "VALIDATION" "TESTING")

for TABLE_NAME in "${DATASET_NAMES[@]}"; do
	bq rm -f -t "${PROJECT_ID}:${ML_DATASET}.${TABLE_NAME}"
done

JOB_NAME="run-hmm-${SHORT_SHA}"

gcloud batch jobs submit "${JOB_NAME}" \
	--location "${REGION}" \
	--config jobs/run_hmm.json
