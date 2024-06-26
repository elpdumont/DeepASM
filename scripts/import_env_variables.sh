#!/bin/bash

# requires yq and readarray package

config_file="config.yaml"

# GCP variables
PROJECT=$(yq e '.GCP.PROJECT' "${config_file}")
export PROJECT

REGION=$(yq e '.GCP.REGION' "${config_file}")
export REGION

BUCKET=$(yq e '.GCP.BUCKET' "${config_file}")
export BUCKET

MODEL_FOLDER=$(yq e '.GCP.MODEL_FOLDER' "${config_file}")
export MODEL_FOLDER

ARTIFACT_REGISTRY_REPO=$(yq e '.GCP.ARTIFACT_REGISTRY_REPO' "${config_file}")
export ARTIFACT_REGISTRY_REPO

BQ_DATASET_EXPIRATION_SEC=$(yq e '.GCP.BQ_DATASET_EXPIRATION_SEC'  "${config_file}")
export BQ_DATASET_EXPIRATION_SEC

CLOUDASM_DATASET=$(yq e '.GCP.CLOUDASM_DATASET'  "${config_file}")
export CLOUDASM_DATASET

# Genomic variables
REFERENCE_GENOME=$(yq e '.GENOMICS.REFERENCE_GENOME' "${config_file}")
export REFERENCE_GENOME

GENOMIC_LENGTH=$(yq e '.GENOMICS.GENOMIC_LENGTH' "${config_file}")
export GENOMIC_LENGTH

ENCODE_BLACKLIST_REGIONS_URL=$(yq e '.GENOMICS.ENCODE_BLACKLIST_REGIONS_URL'  "${config_file}")
export ENCODE_BLACKLIST_REGIONS_URL

NB_NUCLEOTIDES_PER_CLUSTER=$(yq e '.GENOMICS.NB_NUCLEOTIDES_PER_CLUSTER'  "${config_file}")
export NB_NUCLEOTIDES_PER_CLUSTER

MIN_NB_CPG_PER_REGION_IN_REF_GENOME=$(yq e '.GENOMICS.MIN_NB_CPG_PER_REGION_IN_REF_GENOME' "${config_file}")
export MIN_NB_CPG_PER_REGION_IN_REF_GENOME

# ASM variables
MIN_NB_CPG_FOR_ASM=$(yq e '.ASM.MIN_NB_CPG_FOR_ASM' "${config_file}")
export MIN_NB_CPG_FOR_ASM

MIN_NB_CPG_FOR_ASM=$(yq e '.ASM.MIN_NB_CPG_FOR_ASM' "${config_file}")
export MIN_NB_CPG_FOR_ASM

MAX_CPG_COV=$(yq e '.ASM.MAX_CPG_COV' "${config_file}")
export MAX_CPG_COV

MIN_CPG_COV=$(yq e '.ASM.MIN_CPG_COV' "${config_file}")
export MIN_CPG_COV

MAX_P_VALUE=$(yq e '.ASM.MAX_P_VALUE' "${config_file}")
export MAX_P_VALUE

MAX_BH_THRESHOLD=$(yq e '.ASM.MAX_BH_THRESHOLD' "${config_file}")
export MAX_BH_THRESHOLD

MIN_ASM_REGION_EFFECT=$(yq e '.ASM.MIN_ASM_REGION_EFFECT' "${config_file}")
export MIN_ASM_REGION_EFFECT

MIN_NB_CPG_SAME_DIRECTION=$(yq e '.ASM.MIN_NB_CPG_SAME_DIRECTION' "${config_file}")
export MIN_NB_CPG_SAME_DIRECTION

MIN_NB_CONSECUTIVE_CPG_SAME_DIRECTION=$(yq e '.ASM.MIN_NB_CONSECUTIVE_CPG_SAME_DIRECTION' "${config_file}")
export MIN_NB_CONSECUTIVE_CPG_SAME_DIRECTION

# Function to read samples into an array from yq output
read_samples() {
    local samples=()
    while IFS= read -r line; do
        samples+=("${line}")
    done < <(yq e "$1" "${config_file}")
    echo "${samples[@]}"
}

# Initialize arrays for each sample type
training_samples=($(read_samples '.SAMPLES.TRAINING[]'))
validation_samples=($(read_samples '.SAMPLES.VALIDATION[]'))
testing_samples=($(read_samples '.SAMPLES.TESTING[]'))

# Concatenate all samples into a single list
SAMPLE_LIST=("${training_samples[@]}" "${validation_samples[@]}" "${testing_samples[@]}")
export SAMPLE_LIST
echo "List of samples: ${SAMPLE_LIST[*]}"


# Form specific variables used in the scripts
export SAMPLES_DATASET="samples_${GENOMIC_LENGTH}bp"
export MODEL_PATH="${SAMPLES_DATASET}/models"
export DATA_PATH="${SAMPLES_DATASET}/datasets"
export ML_DATASET="ml_${GENOMIC_LENGTH}bp_${SHORT_SHA}"
export REFG_DATASET="${REFERENCE_GENOME}_${GENOMIC_LENGTH}bp_refgenome"
export PYTHON_IMAGE="${REGION}-docker.pkg.dev/${PROJECT}/${ARTIFACT_REGISTRY_REPO}/python:${SHORT_SHA}"
export BASH_IMAGE="${REGION}-docker.pkg.dev/${PROJECT}/${ARTIFACT_REGISTRY_REPO}/bash:${SHORT_SHA}"
export CLOUDASM_TABLES=("context_filtered" "cpg_asm" "cpg_read_genotype")
export NB_SAMPLES=${#SAMPLE_LIST[@]}
export NB_CLOUDASM_TABLES=${#CLOUDASM_TABLES[@]}

nb_tasks_for_cloudasm_tables=$((NB_SAMPLES * NB_CLOUDASM_TABLES))

# Create datasets in BQ if they do not exist.
echo "Creating the datasets with an expiration"
DATASET_LIST=("${SAMPLES_DATASET}" "${REFG_DATASET}" "${ML_DATASET}")
for DATASET in "${DATASET_LIST[@]}"; do
    echo "Processing dataset: ${DATASET}"
    if bq ls --project_id="${PROJECT}" | grep -w "${DATASET}"; then
        echo "Dataset ${DATASET} already exists in project ${PROJECT}."
    else
        # Create the dataset since it does not exist
        bq mk   --dataset \
                --location="${REGION}" \
                --default_table_expiration="${BQ_DATASET_EXPIRATION_SEC}" \
                "${PROJECT}:${DATASET}"

        echo "Dataset ${DATASET} created in project ${PROJECT}."
    fi
done

export LC_NUMERIC="en_US.UTF-8"

function format_number_with_comma() {
    local number=$1
    # Remove existing non-numeric characters (e.g., commas) to ensure the input is a valid integer
    #clean_number=$(echo "${number_string}" | sed 's/[^0-9]//g')
    #local clean_number
    printf "%'d\n" "${number}"
}

function execute_query() {
    local QUERY=$1
    #echo "${QUERY}"
    local OUTPUT=$(bq query --use_legacy_sql=false --format=json "${QUERY}")
    #echo "${OUTPUT}"
    local NUMBER=$(echo "${OUTPUT}" | jq -r '.[0] | to_entries | .[0].value')
    #echo "${NUMBER}"
    # Format the number with commas and echo it
    format_number_with_comma "${NUMBER}"
}

# Copy the template to a new file that can be safely modified
(cd batch-jobs && for file in *.json.tpl; do cp "${file}" "${file%.json.tpl}.json"; done)

echo ${ML_DATASET}

# Replace placeholders with actual values
for file in batch-jobs/*.json; do
    echo "Processing file ${file}"
    sed -i '' "s#REGION_PH#${REGION}#g" "${file}"
    sed -i '' "s#PROJECT_PH#${PROJECT}#g" "${file}"
    sed -i '' "s#PYTHON_IMAGE_PH#${PYTHON_IMAGE}#g" "${file}"
    sed -i '' "s#BASH_IMAGE_PH#${BASH_IMAGE}#g" "${file}"
    sed -i '' "s#ML_DATASET_PH#${ML_DATASET}#g" "${file}"
    sed -i '' "s#CLOUDASM_DATASET_PH#${CLOUDASM_DATASET}#g" "${file}"
    sed -i '' "s#ML_DATASET_PH#${ML_DATASET}#g" "${file}"
    sed -i '' "s#SAMPLES_DATASET_PH#${SAMPLES_DATASET}#g" "${file}"
    sed -i '' "s#nb_tasks_for_cloudasm_tables_ph#${nb_tasks_for_cloudasm_tables}#g" "${file}"
    sed -i '' "s#CLOUDASM_TABLES_PH#${CLOUDASM_TABLES[*]}#g" "${file}"
    sed -i '' "s#SAMPLE_LIST_PH#${SAMPLE_LIST[*]}#g" "${file}"
    sed -i '' "s#NB_SAMPLES_PH#${NB_SAMPLES}#g" "${file}"
    sed -i '' "s#SHORT_SHA_PH#${SHORT_SHA}#g" "${file}"
    sed -i '' "s#MODEL_PATH_PH#${MODEL_PATH}#g" "${file}"
    sed -i '' "s#ML_MODE_PH#${ML_MODE}#g" "${file}"
    sed -i '' "s#HMM_MODEL_PH#${HMM_MODEL}#g" "${file}"
done


(cd deploy && for file in *.yaml.tpl; do cp "${file}" "${file%.yaml.tpl}.yaml"; done)

deploy_file="deploy/cloudbuild.yaml"

sed -i '' "s#REGION_PH#${REGION}#g" "${deploy_file}"
sed -i '' "s#PROJECT_PH#${PROJECT}#g" "${deploy_file}"
sed -i '' "s#ARTIFACT_REGISTRY_REPO_PH#${ARTIFACT_REGISTRY_REPO}#g" "${deploy_file}"
sed -i '' "s#IMAGE_TAG_PH#${SHORT_SHA}#g" "${deploy_file}"


# Make all bash files executable
chmod +x app/prepare-ref-genome-bash/*.sh
chmod +x app/prepare-samples/bash/*.sh