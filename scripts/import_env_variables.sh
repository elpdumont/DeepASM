#!/bin/bash

# requires yq and readarray package

config_file="config/config.yaml"


# GCP variables
export PROJECT_ID=$(yq e '.GCP.PROJECT_ID' config/config.yaml)
export REGION=$(yq e '.GCP.REGION' config/config.yaml)
export PYTHON_IMAGE=$(yq e '.GCP.IMAGE' config/config.yaml)
export BUCKET_NAME=$(yq e '.GCP.BUCKET_NAME' config/config.yaml)
export BQ_DATASET_EXPIRATION_SEC=$(yq e '.GCP.BQ_DATASET_EXPIRATION_SEC'  config/config.yaml)
export CLOUDASM_OUTPUT_DATASET=$(yq e '.GCP.CLOUDASM_OUTPUT_DATASET'  config/config.yaml)
export CLOUDASM_STANDARD_REGIONS_DATASET=$(yq e '.GCP.CLOUDASM_STANDARD_REGIONS_DATASET'  config/config.yaml)


# Genomic variables
export REFERENCE_GENOME=$(yq e '.GENOMICS.REFERENCE_GENOME' config/config.yaml)
export GENOMIC_LENGTH=$(yq e '.GENOMICS.GENOMIC_LENGTH' config/config.yaml)
export ENCODE_BLACKLIST_REGIONS_URL=$(yq e '.GENOMICS.ENCODE_BLACKLIST_REGIONS_URL'  config/config.yaml)
export NB_NUCLEOTIDES_PER_CLUSTER=$(yq e '.GENOMICS.NB_NUCLEOTIDES_PER_CLUSTER'  config/config.yaml)
export MIN_NB_CPG_PER_REGION_IN_REF_GENOME=$(yq e '.GENOMICS.MIN_NB_CPG_PER_REGION_IN_REF_GENOME' config/config.yaml)

# ASM variables
export MIN_NB_CPG_FOR_ASM=$(yq e '.ASM.MIN_NB_CPG_FOR_ASM' config/config.yaml)
export MAX_CPG_COV=$(yq e '.ASM.MAX_CPG_COV' config/config.yaml)
export MIN_NB_CPG_FOR_ASM=$(yq e '.ASM.MIN_NB_CPG_FOR_ASM' config/config.yaml)
export P_VALUE=$(yq e '.ASM.P_VALUE' config/config.yaml)
export BH_THRESHOLD=$(yq e '.ASM.BH_THRESHOLD' config/config.yaml)
export ASM_REGION_EFFECT=$(yq e '.ASM.ASM_REGION_EFFECT' config/config.yaml)
export NB_CPG_SAME_DIRECTION_ASM=$(yq e '.ASM.NB_CPG_SAME_DIRECTION_ASM' config/config.yaml)
export CONSECUTIVE_CPG_WITH_SAME_DIRECTION_ASM=$(yq e '.ASM.CONSECUTIVE_CPG_WITH_SAME_DIRECTION_ASM' config/config.yaml)


# Function to read samples into an array from yq output
read_samples() {
    local samples=()
    while IFS= read -r line; do
        samples+=("$line")
    done < <(yq e "$1" "$config_file")
    echo "${samples[@]}"
}

# Initialize arrays for each sample type
training_samples=($(read_samples '.SAMPLES.TRAINING[]'))
validation_samples=($(read_samples '.SAMPLES.VALIDATION[]'))
testing_samples=($(read_samples '.SAMPLES.TESTING[]'))

# Concatenate all samples into a single list
SAMPLE_LIST=("${training_samples[@]}" "${validation_samples[@]}" "${testing_samples[@]}")



# Form specific folders used in the scripts
export CLOUDASM_DATASET="${REFERENCE_GENOME}_${GENOMIC_LENGTH}_cloudasm"
export ML_DATASET="${REFERENCE_GENOME}_${GENOMIC_LENGTH}_ml"
export REFG_FOLDER="${REFERENCE_GENOME}_${GENOMIC_LENGTH}bp_refgenome"


# Create datasets in BQ if they do not exist.
echo "Creating the datasets with an expiration"
DATASET_LIST=("$CLOUDASM_STANDARD_REGIONS_DATASET" "$REFG_FOLDER")
for DATASET in "${DATASET_LIST[@]}"; do
    echo "Processing dataset: ${DATASET}"
    if bq ls --project_id="${PROJECT_ID}" | grep -w "${DATASET}"; then
        echo "Dataset ${DATASET} already exists in project ${PROJECT_ID}."
    else
        # Create the dataset since it does not exist
        bq mk --project_id="${PROJECT_ID}" \
                --dataset \
                --location="${REGION}" \
                --default_table_expiration="${BQ_DATASET_EXPIRATION_SEC}" \
                "${PROJECT_ID}:${DATASET}"

        echo "Dataset ${DATASET} created in project ${PROJECT_ID}."
    fi
done

export LC_NUMERIC="en_US.UTF-8"

function format_number_with_comma() {
    local number_string=$1
    # Remove existing non-numeric characters (e.g., commas) to ensure the input is a valid integer
    local clean_number=$(echo "$number_string" | sed 's/[^0-9]//g')
    printf "%'d\n" "${clean_number}"
}

function execute_query() {
    local QUERY=$1
    local OUTPUT=$(bq query --use_legacy_sql=false --format=json "${QUERY}")
    #echo "${OUTPUT}"
    local NUMBER=$(echo "${OUTPUT}" | jq -r '.[0] | to_entries | .[0].value')
    #echo "${NUMBER}"
    # Format the number with commas and echo it
    format_number_with_comma "${NUMBER}"
}


# Update the jobs file
# Copy the template to a new file that can be safely modified
mkdir jobs
cp jobs_templates/* jobs/

# Replace placeholders with actual values
for file in jobs/process_json.json jobs/run_hmm.json jobs/calculate_wilcoxon_for_regions.json; do
    echo "Processing file ${file}"
    sed -i '' "s#PYTHON_IMAGE_PLACEHOLDER#${PYTHON_IMAGE}#g" "${file}"
    sed -i '' "s/IMAGE_TAG_PLACEHOLDER/${SHORT_SHA}/g" "${file}"
    sed -i '' "s/ML_DATASET_ID_PLACEHOLDER/${ML_DATASET}/g" "${file}"
    sed -i '' "s/CLOUDASM_DATASET_ID_PLACEHOLDER/${CLOUDASM_DATASET}/g" "${file}"
    sed -i '' "s/CLOUDASM_STANDARD_REGIONS_DATASET_PLACEHOLDER/${CLOUDASM_STANDARD_REGIONS_DATASET}/g" "${file}"
done




# Make all files executable
chmod +x src/*.sh