#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_EPI}.hg19_cpg_regions_${CHR}_${LOWER_B}_${UPPER_B} \
    --replace=true \
    "
    WITH 
        -- Select genomic regions within the specificed chromosome and inf/sup limits
        GENOMIC_REGIONS AS (
            SELECT chr_region, region_inf, region_sup, annotate_ref
            FROM ${DATASET_EPI}.hg19_regions_${GENOMIC_INTERVAL}bp
            WHERE 
                chr_region = '${CHR}'
                AND region_sup <= ${UPPER_B}
                AND region_inf >= ${LOWER_B}
        ),
        CPG_POS AS (
            -- Select the CpG coordinates in the table where there is chr, CpG pos (inf), CpG pos (sup)
            SELECT region_inf
            FROM ${DATASET_EPI}.hg19_CpG_pos
            WHERE 
                chr_region = '${CHR}'
                AND region_inf <= ${UPPER_B}
                AND region_inf >= ${LOWER_B}
        ),
        -- Find all the CpG in each genomic region
        REGION_CPG AS (
            SELECT chr_region, GENOMIC_REGIONS.region_inf, GENOMIC_REGIONS.region_sup, annotate_ref FROM GENOMIC_REGIONS
            INNER JOIN CPG_POS
            ON GENOMIC_REGIONS.region_inf <= CPG_POS.region_inf AND CPG_POS.region_inf <= GENOMIC_REGIONS.region_sup
        ),
        -- Count the number of CpGs in each region
        REGION_CPG_COUNT AS (
            SELECT 
                chr_region AS chr, 
                region_inf, 
                region_sup, 
                annotate_ref,
                COUNT(*) AS region_nb_cpg 
            FROM REGION_CPG
            GROUP BY chr_region, region_inf, region_sup, annotate_ref
        )
        -- Keep the regions where there are at least 3 CpGs in the reference genome
        SELECT * 
        FROM REGION_CPG_COUNT
        WHERE region_nb_cpg >= 3
    "