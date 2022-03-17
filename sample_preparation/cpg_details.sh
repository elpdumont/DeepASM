#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_cpg_details_${GENOMIC_INTERVAL}bp \
    --replace=true \
    "
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
            STRUCT(pos, meth, read_id)
            ) AS cpg_details

    FROM ${DATASET_PRED}.${SAMPLE}_cpg_regions_${GENOMIC_INTERVAL}bp_clean
    GROUP BY chr, region_inf, region_sup, region_nb_cpg, dnase, encode_ChiP_V2, tf_motifs

        "