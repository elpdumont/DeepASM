#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_OUT}.asm_read_cpg_tf_${CHR} \
    --replace=true \
    "
    WITH 
    ASM AS (
        SELECT * FROM ${DATASET_OUT}.asm_read_cpg_arrays
        WHERE chr = '${CHR}'
    ),
    TF AS (
        SELECT 
            chr AS chr_tf, 
            chr_start, 
            chr_end, 
            name AS tf_name,
            score AS tf_score
        FROM ${DATASET_EPI}.encode_ChiP_V2
        WHERE chr = '${CHR}'
    )
    SELECT * EXCEPT(chr_start, chr_end, chr_tf)
    FROM ASM 
    LEFT JOIN TF
    ON chr_start <= region_sup + ${EPI_REGION}
        AND chr_end >= region_inf - ${EPI_REGION}
        AND chr = chr_tf
    "