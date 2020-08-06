#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_cpg_read_${GENOMIC_INTERVAL}bp_name \
    --replace=true \
    "
    SELECT *, '${SAMPLE}' AS sample 
    FROM ${DATASET_PRED}.${SAMPLE}_cpg_read_${GENOMIC_INTERVAL}bp
    "