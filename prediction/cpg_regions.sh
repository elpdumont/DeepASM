#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_cpg_regions_${CHR}_${LOWER_B}_${UPPER_B} \
    --replace=true \
    "
    WITH 
        CPG_REGIONS AS (
            SELECT * 
            FROM ${DATASET_PRED}.regions_${GENOMIC_INTERVAL}bp
            WHERE 
                chr_region = '${CHR}'
                AND region_sup <= ${UPPER_B}
                AND region_inf >= ${LOWER_B}
        ),
        CONTEXT AS (
            SELECT *
            FROM ${DATASET_CONTEXT}.${SAMPLE}_context_filtered
            WHERE 
                chr = '${CHR}'
                AND pos <= ${UPPER_B}
                AND pos >= ${LOWER_B}
        )
        SELECT * EXCEPT(chr_region) FROM CPG_REGIONS
        INNER JOIN CONTEXT
        ON pos >= region_inf AND pos < region_sup
    "