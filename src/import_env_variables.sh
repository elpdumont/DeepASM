#!/bin/bash

export PROJECT_ID=$(yq e '.GCP.PROJECT_ID' config/config.yaml)
export REGION=$(yq e '.GCP.REGION' config/config.yaml)
export PYTHON_IMAGE=$(yq e '.GCP.IMAGE' config/config.yaml)
export GENOMIC_LENGTH=$(yq e '.GENOMICS.GENOMIC_LENGTH' config/config.yaml)
export BUCKET_NAME=$(yq e '.GCP.BUCKET_NAME' config/config.yaml)
export REFERENCE_GENOME=$(yq e '.GENOMICS.REFERENCE_GENOME' config/config.yaml)
export MIN_NB_CPG_PER_REGION_IN_REF_GENOME=$(yq e '.GENOMICS.MIN_NB_CPG_PER_REGION_IN_REF_GENOME' config/config.yaml)
export NB_NUCLEOTIDES_PER_CLUSTER=$(yq e '.GENOMICS.NB_NUCLEOTIDES_PER_CLUSTER'  config/config.yaml)
export REF_GENOME_DATASET_EXPIRATION_SEC=$(yq e '.GENOMICS.REF_GENOME_DATASET_EXPIRATION_SEC'  config/config.yaml)
export ENCODE_BLACKLIST_REGIONS_URL=$(yq e '.GENOMICS.ENCODE_BLACKLIST_REGIONS_URL'  config/config.yaml)

# Form specific folders used in the scripts
export CLOUDASM_DATASET="${REFERENCE_GENOME}_${GENOMIC_LENGTH}_cloudasm"
export ML_DATASET="${REFERENCE_GENOME}_${GENOMIC_LENGTH}_ml"
export REFG_FOLDER="${REFERENCE_GENOME}_${GENOMIC_LENGTH}bp_refgenome"


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
    local NUMBER=$(echo "${OUTPUT}" | jq -r '.[0].f0_')
    # Format the number with commas and echo it
    format_number_with_comma "${NUMBER}"
}

# Define the format_number_with_comma function
