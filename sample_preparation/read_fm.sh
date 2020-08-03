#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_read_fm_${GENOMIC_INTERVAL}bp \
    --replace=true \
    "
    WITH 
        DATASETS_JOINED AS (
            SELECT * 
            FROM ${DATASET_PRED}.${SAMPLE}_cpg_regions_${GENOMIC_INTERVAL}bp
        ),
        READ_FRAC_METHYL AS (
            SELECT 
                ROUND(SUM(meth)/SUM(cov),3) AS fm,
                chr,
                region_inf,
                region_sup,
                region_nb_cpg, 
                dnase,
                encode_ChiP_V2,
                tf_motifs
            FROM DATASETS_JOINED
            GROUP BY read_id, chr, region_inf, region_sup, region_nb_cpg, dnase, encode_ChiP_V2, tf_motifs
        )
        SELECT
            chr,
            region_inf,
            region_sup,
            region_nb_cpg, 
            dnase,
            encode_ChiP_V2,
            tf_motifs,
            ARRAY_AGG (STRUCT (fm)) AS read
        FROM READ_FRAC_METHYL
        GROUP BY chr, region_inf, region_sup, region_nb_cpg, dnase, encode_ChiP_V2, tf_motifs
    "