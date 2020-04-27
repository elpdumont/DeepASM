#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_OUT}.asm \
    --append_table \
    "
    SELECT
        '${SAMPLE}' AS sample,
        *
    FROM ${DATASET_IN}.${SAMPLE}_asm_snp
    "