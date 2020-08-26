#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_cpg_regions_${GENOMIC_INTERVAL}bp_clean \
    --replace=true \
    "
    WITH 
        GROUP_CPG AS (
        -- Creates a list of all CpG sites with their coverage
            SELECT 
                chr, 
                pos, 
                SUM(cov) AS cov,
            FROM ${DATASET_PRED}.${SAMPLE}_cpg_regions_${GENOMIC_INTERVAL}bp
            GROUP BY chr, pos
        ),
        GOOD_CPG AS (
            SELECT
                chr,
                pos
            FROM GROUP_CPG
            WHERE cov >= ${MIN_CPG_COV} AND cov <= ${MAX_CPG_COV}
        )
        -- recreate a long table with CpG and read information
        SELECT
            t1.chr,
            t1.region_inf,
            t1.region_sup,
            t1.region_nb_cpg,
            t1.dnase,
            t1.encode_ChiP_V2,
            t1.tf_motifs,
            t1.pos,
            t1.meth,
            t1.cov,
            t1.read_id
        FROM ${DATASET_PRED}.${SAMPLE}_cpg_regions_${GENOMIC_INTERVAL}bp t1
        INNER JOIN GOOD_CPG t2
        ON t1.chr = t2.chr AND t1.pos = t2.pos
        "