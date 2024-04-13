#!/bin/bash
bq query \
    --use_legacy_sql=false \
    --destination_table "${FOLDER}".regions_w_cpg \
    --replace=true \
    --clustering_fields=chr \
    --range_partitioning=clustering_index,0,4000,1 \
    "
    WITH 
        -- Select genomic regions within the specificed chromosome and inf/sup limits
        GENOMIC_REGIONS AS (
            SELECT chr, region_inf, region_sup, region_center, clustering_index
            FROM ${FOLDER}.regions
        ),
        CPG_POS AS (
            -- Select the CpG coordinates in the table where there is chr, CpG pos (inf), CpG pos (sup)
            SELECT chr, region_inf, region_sup, clustering_index
            FROM ${FOLDER}.CpG_pos
        ),
        -- Find all the CpG in each genomic region
        REGION_CPG AS (
            SELECT
                GENOMIC_REGIONS.chr,
                GENOMIC_REGIONS.region_inf,
                GENOMIC_REGIONS.region_sup,
                GENOMIC_REGIONS.region_center,
                GENOMIC_REGIONS.clustering_index
            FROM GENOMIC_REGIONS
            INNER JOIN CPG_POS
            ON GENOMIC_REGIONS.clustering_index = CPG_POS.clustering_index AND
               GENOMIC_REGIONS.region_inf <= CPG_POS.region_inf AND
               CPG_POS.region_inf <= GENOMIC_REGIONS.region_sup
        ),
        -- Count the number of CpGs in each region
        REGION_CPG_COUNT AS (
            SELECT 
                chr,
                region_inf, 
                region_sup, 
                region_center,
                COUNT(*) AS region_nb_cpg,
                clustering_index
            FROM REGION_CPG
            GROUP BY chr, region_inf, region_sup, region_center, clustering_index
        )
        -- Keep the regions where there are at least MIN_NB_CPG_PER_REGION CpGs in the reference genome
        SELECT * 
        FROM REGION_CPG_COUNT
        WHERE region_nb_cpg >= ${MIN_NB_CPG_PER_REGION_IN_REF_GENOME}
    "