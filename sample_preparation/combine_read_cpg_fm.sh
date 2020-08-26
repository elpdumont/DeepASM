#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_cpg_read_${GENOMIC_INTERVAL}bp \
    --replace=true \
    "
    SELECT 
        t1.chr,
        t1.region_inf,
        t1.region_sup,
        t1.region_nb_cpg,
        t1.nb_cpg_found,
        ARRAY_LENGTH(t2.read) AS nb_reads,
        t1.dnase, 
        t1.encode_ChiP_V2, 
        t1.tf_motifs,
        t1.cpg,
        t2.read
    FROM ${DATASET_PRED}.${SAMPLE}_cpg_fm_${GENOMIC_INTERVAL}bp t1
    INNER JOIN ${DATASET_PRED}.${SAMPLE}_read_fm_${GENOMIC_INTERVAL}bp t2
    ON 
        t1.region_inf =  t2.region_inf AND
        t1.region_sup = t2.region_sup AND
        t1.chr = t2.chr AND
        nb_cpg_found >= ${MIN_NB_CPG}
    "
