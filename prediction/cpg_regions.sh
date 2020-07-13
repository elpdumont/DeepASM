#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_cpg_regions_${CHR}_${LOWER_B}_${UPPER_B} \
    --replace=true \
    "
    WITH 
        CPG_REGIONS AS (
            SELECT 
                chr AS chr_region, 
                region_inf, 
                region_sup, 
                region_nb_cpg,
                dnase,
                encode_ChiP_V2,
                tf_motifs
            FROM ${DATASET_EPI}.hg19_cpg_regions_${GENOMIC_INTERVAL}bp_annotated
            WHERE 
                chr = '${CHR}'
                AND region_sup <= ${UPPER_B}
                AND region_inf >= ${LOWER_B}
        ),
        CONTEXT AS (
            SELECT *
            FROM ${DATASET_CONTEXT}.${SAMPLE}_context_filtered
            WHERE 
                chr = '${CHR}'
                AND pos <= ${UPPER_B}
                AND pos >= ${LOWER_B}
        )
        SELECT 
            chr,
            region_inf,
            region_sup,
            region_nb_cpg,
            dnase,
            encode_ChiP_V2,
            tf_motifs,
            pos,
            meth, 
            cov, 
            read_id
        FROM CPG_REGIONS
        INNER JOIN CONTEXT
        ON pos >= region_inf AND pos < region_sup
    "