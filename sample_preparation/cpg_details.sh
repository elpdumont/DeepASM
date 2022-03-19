#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_genomic_window_${GENOMIC_INTERVAL}bp \
    --replace=true \
    "
    WITH TMP AS (
        SELECT 
            chr,
            region_inf,
            region_sup,
            read_id,
            ROUND(SUM(meth)/SUM(cov),3) AS read_fm,
            COUNT(*) AS nb_cpg,
            ARRAY_AGG (pos-region_inf + 1) AS pos_array,
            ARRAY_AGG (meth) AS meth_array
        FROM ${DATASET_PRED}.${SAMPLE}_cpg_regions_${GENOMIC_INTERVAL}bp_clean
        GROUP BY chr, region_inf, region_sup, read_id
        )
        SELECT 
            chr,
            region_inf,
            region_sup,
            ARRAY_AGG(
                STRUCT(read_id, read_fm, nb_cpg, pos_array, meth_array)
            ) AS window_details
        FROM TMP
        GROUP BY chr, region_inf, region_sup
    "

