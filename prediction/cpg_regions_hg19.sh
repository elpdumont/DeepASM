#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_EPI}.cpg_regions_${CHR}_${LOWER_B}_${UPPER_B} \
    --replace=true \
    "
    WITH 
        GENOMIC_REGIONS AS (
            SELECT * 
            FROM ${DATASET_EPI}.regions_${GENOMIC_INTERVAL}bp
            WHERE 
                chr_region = '${CHR}'
                AND region_sup <= ${UPPER_B}
                AND region_inf >= ${LOWER_B}
        ),
        CPG_POS AS (
            SELECT inf
            FROM ${DATASET_EPI}.hg19_CpG_pos
            WHERE 
                chr = '${CHR}'
                AND inf <= ${UPPER_B}
                AND inf >= ${LOWER_B}
        ),
        REGION_CPG AS (
            SELECT * EXCEPT(inf) FROM GENOMIC_REGIONS
            INNER JOIN CPG_POS
            ON inf >= region_inf AND inf < region_sup
        ),
        REGION_CPG_COUNT AS (
            SELECT 
                chr_region, 
                region_inf, 
                region_sup, 
                COUNT(*) AS region_nb_cpg 
            FROM REGION_CPG
            GROUP BY chr_region, region_inf, region_sup
        )
        SELECT * 
        FROM REGION_CPG_COUNT
        WHERE region_nb_cpg >= 3
    "