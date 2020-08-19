#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_cpg_fm_${GENOMIC_INTERVAL}bp \
    --replace=true \
    "
    WITH 
        DATASETS_JOINED AS (
            SELECT * 
            FROM ${DATASET_PRED}.${SAMPLE}_cpg_regions_${GENOMIC_INTERVAL}bp
        ),
        CPG_FRAC_METHYL AS (
        -- Creates a list of all CpG sites with their respective fractional
        -- methylation and their CpG region
            SELECT 
                chr, 
                pos, 
                SUM(cov) AS cov,
                ROUND(SUM(meth)/SUM(cov),3) AS fm,
                region_inf,
                region_sup,
                region_nb_cpg, 
                dnase,
                encode_ChiP_V2,
                tf_motifs
            FROM DATASETS_JOINED
            GROUP BY chr, pos, region_inf, region_sup, region_nb_cpg, dnase, encode_ChiP_V2, tf_motifs
        ),
        GROUP_CPG_INFO_BY_REGION AS (
        SELECT
            chr,
            region_inf,
            region_sup,
            region_nb_cpg,
            dnase,
            encode_ChiP_V2,
            tf_motifs,
            COUNT(*) AS nb_cpg_found,
            ARRAY_AGG(
                STRUCT(fm, cov, pos)
                ) AS cpg
        FROM CPG_FRAC_METHYL
        WHERE cov >= 10
        GROUP BY region_inf, region_sup, chr, region_nb_cpg, dnase, encode_ChiP_V2, tf_motifs
        )
        SELECT * FROM GROUP_CPG_INFO_BY_REGION
        WHERE nb_cpg_found >= 3
        "