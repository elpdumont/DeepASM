#!/bin/bash


bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_OUT}.asm_read_cpg_dnase_${CHR} \
    --replace=true \
    "
    WITH 
    ASM AS (
        SELECT * 
        FROM ${DATASET_OUT}.asm_read_cpg_arrays
        WHERE chr = '${CHR}'
    ),
    DNASE AS (
        SELECT 
            chr AS chr_dnase,
            chr_start,
            chr_end, 
            score AS score_dnase
        FROM ${DATASET_EPI}.dnase
        WHERE chr = '${CHR}'
    )
    SELECT * EXCEPT(chr_start, chr_end, chr_dnase)
    FROM ASM 
    LEFT JOIN DNASE
    ON chr_start <= region_sup + ${EPI_REGION}
    AND chr_end >= region_inf - ${EPI_REGION}
    AND chr = chr_dnase
    "
