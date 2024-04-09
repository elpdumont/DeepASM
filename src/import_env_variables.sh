#!/bin/bash

export PROJECT_ID=$(yq e '.GCP.PROJECT_ID' config/config.yaml)
export REGION=$(yq e '.GCP.REGION' config/config.yaml)
export PYTHON_IMAGE=$(yq e '.GCP.IMAGE' config/config.yaml)
export GENOMIC_LENGTH=$(yq e '.GENOMICS.GENOMIC_LENGTH' config/config.yaml)
export BUCKET_NAME=$(yq e '.GCP.BUCKET_NAME' config/config.yaml)
export REFERENCE_GENOME=$(yq e '.GENOMICS.REFERENCE_GENOME' config/config.yaml)
export MIN_NB_CPG_PER_REGION=$(yq e '.GENOMICS.MIN_NB_CPG_PER_REGION' config/config.yaml)
export NB_NUCLEOTIDES_PER_CLUSTER=$(yq e '.GENOMICS.NB_NUCLEOTIDES_PER_CLUSTER'  config/config.yaml)

export REF_GENOME_DATASET_EXPIRATION_SEC=$(yq e '.GENOMICS.REF_GENOME_DATASET_EXPIRATION_SEC'  config/config.yaml)

export CLOUDASM_DATASET="${REFERENCE_GENOME}_${GENOMIC_LENGTH}_cloudasm"
export ML_DATASET="${REFERENCE_GENOME}_${GENOMIC_LENGTH}_ml"