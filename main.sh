#!/bin/bash

SHORT_SHA="$(git rev-parse --short HEAD)"
echo "SHORT_SHA: ${SHORT_SHA}"

# Import environmental variables
source scripts/import_env_variables.sh

# Create container images on GCP (for Python and bash)
deploy/deploy.sh



docker run -it \
-v ~/.config/gcloud/application_default_credentials.json:/appuser/.config/gcloud/application_default_credentials.json:ro \
-e GOOGLE_APPLICATION_CREDENTIALS=/appuser/.config/gcloud/application_default_credentials.json \
us-east1-docker.pkg.dev/hmh-em-deepasm/docker-repo/bash:"${SHORT_SHA}" \
/bin/bash

docker run -it \
-v ~/.config/gcloud/application_default_credentials.json:/appuser/.config/gcloud/application_default_credentials.json:ro \
-e GOOGLE_APPLICATION_CREDENTIALS=/appuser/.config/gcloud/application_default_credentials.json \
us-east1-docker.pkg.dev/hmh-em-deepasm/docker-repo/python:"9f17e4a" \
/bin/bash


echo "Preparing standard genomic regions using the reference genome"
app/prepare-ref-genome-bash/prepare_regions_w_cpg_using_refg.sh

echo "Formatting CloudASM output to standard genomic regions"
app/prepare-samples/bash/format_cloudasm_to_standard_regions.sh


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






#---------------------------------------------
# PREPARE FEATURES FOR ML. REQUIRE THE JSON FILES OUT THE CLOUDASM SCRIPT (FORMAT CLOUDASM TO STANDARD REGIONS)

echo "Preparing the features for all the regions (even if they do not have ASM flagged)"
NB_FILES_PER_TASK="50"
NB_FILES_TO_PROCESS=$(gsutil ls gs://"${BUCKET}"/"${DATA_PATH}"/after_cloudasm/all_regions/* | wc -l | awk '{print $1}')
TASK_COUNT=$(( (${NB_FILES_TO_PROCESS} + ${NB_FILES_PER_TASK} - 1) / ${NB_FILES_PER_TASK} ))
sed -i '' "s#NB_FILES_PER_TASK_PH#${NB_FILES_PER_TASK}#g" "batch-jobs/prepare_features_for_ML.json"
sed -i '' "s#TASK_COUNT_PH#${TASK_COUNT}#g" "batch-jobs/prepare_features_for_ML.json"


gcloud batch jobs submit "prepare-features-for-ml-${SHORT_SHA}" \
	--location "${REGION}" \
	--config batch-jobs/prepare_features_for_ML.json


nb_regions=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${ML_DATASET}.features_wo_hmm")
nb_regions_w_data=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${ML_DATASET}.features_wo_hmm WHERE cpg_directional_fm IS NOT NULL ")
nb_regions_w_data_and_asm_flagged=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${ML_DATASET}.features_wo_hmm WHERE cpg_directional_fm IS NOT NULL AND asm IS NOT NULL")
echo "Features were prepared for ${nb_regions} regions. Among these, we could extract features for ${nb_regions_w_data} regions. Among these, we have ${nb_regions_w_data_and_asm_flagged} regions with features and ASM."

echo "Exporting the dataset with features (excluding HMM) to the bucket"
bq extract --destination_format=NEWLINE_DELIMITED_JSON "${PROJECT}:${ML_DATASET}.features_wo_hmm" gs://"${BUCKET}"/"${DATA_PATH}"/features_wo_hmm/features_wo_hmm_*.json




#---------------------------------------------
# FIT TRANSFORMER AND RNN

gcloud batch jobs submit "transformer-rnn-1d-${SHORT_SHA}" \
	--location "${REGION}" \
	--config batch-jobs/run_transformer_and_rnn_1d.json

#---------------------------------------------
# Find the number of states that works best
NB_STATES="10"
sed -i '' "s#TASK_COUNT_PH#${NB_STATES}#g" "batch-jobs/find_nb_states_for_HMM.json"

gcloud batch jobs submit "nb-states-hmm-${SHORT_SHA}-4" \
	--location "${REGION}" \
	--config batch-jobs/find_nb_states_for_HMM.json


#---------------------------------------------
echo "Fitting an HMM model on the training set and infering the states-based features for all datasets"
TOTAL_TASKS="120"
sed -i '' "s#TOTAL_TASK_PH#${TOTAL_TASKS}#g" "batch-jobs/derive_features_from_HMM.json"

gcloud batch jobs submit "compute-hmm-${SHORT_SHA}-2" \
	--location "${REGION}" \
	--config batch-jobs/derive_features_from_HMM.json



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


#-----------------------------------------------------------

# ML: PERFORM RANDOM SEARCH FOR TREE MODELS

gcloud batch jobs submit "tree-search-${SHORT_SHA}" \
	--location "${REGION}" \
	--config batch-jobs/perform_random_search_tree.json


