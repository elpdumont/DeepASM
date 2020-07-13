#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_cpg_read_${GENOMIC_INTERVAL}bp \
    --replace=true \
    "
    WITH 
        CPG_INFO AS (
            SELECT * 
            FROM ${DATASET_PRED}.${SAMPLE}_cpg_fm_${GENOMIC_INTERVAL}bp
        ),
        READ_INFO AS (
            SELECT 
                chr AS chr_read,
                region_inf AS region_inf_read,
                region_sup AS region_sup_read,
                read
            FROM ${DATASET_PRED}.${SAMPLE}_read_fm_${GENOMIC_INTERVAL}bp
        )
        SELECT 
            chr,
            region_inf,
            region_sup,
            region_nb_cpg,
            nb_cpg_found,
            dnase, 
            encode_ChiP_V2, 
            tf_motifs,
            cpg,
            read
        FROM CPG_INFO
        INNER JOIN READ_INFO
        ON 
            region_inf =  region_inf_read 
            AND region_sup = region_sup_read
            AND chr_read = chr
    "
