#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_regions \
    --replace=true \
    "
    WITH 
        CPG_INFO AS (
            SELECT * 
            FROM ${DATASET_PRED}.${SAMPLE}_cpg_fm
        ),
        READ_INFO AS (
            SELECT 
                chr AS chr_read,
                region_inf AS region_inf_read,
                region_sup AS region_sup_read,
                read
            FROM ${DATASET_PRED}.${SAMPLE}_read_fm
        )
        SELECT 
            chr,
            region_inf,
            region_sup,
            CAST(FLOOR((region_inf + region_sup)/2) AS INT64) AS enrich_ref,
            nb_cpg_found,
            cpg,
            read
        FROM CPG_INFO
        INNER JOIN READ_INFO
        ON 
            region_inf =  region_inf_read 
            AND region_sup = region_sup_read
            AND chr_read = chr
    "
