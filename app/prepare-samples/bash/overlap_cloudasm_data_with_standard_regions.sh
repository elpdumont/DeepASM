#!/bin/bash

TEMP_TABLE="${SAMPLES_DATASET}.${sample}_${table}_temp"
TEMP_TABLE_2="${SAMPLES_DATASET}.${sample}_${table}_temp2"

# Create and populate the temporary table, adding the 'sample' column
bq query \
    --use_legacy_sql=false \
    --destination_table="${TEMP_TABLE}" \
    --replace=true \
    --clustering_fields=sample,chr \
    --range_partitioning=clustering_index,0,4000,1 \
    "
    SELECT
        '${sample}' AS sample,
        CAST(FLOOR((c.absolute_nucleotide_pos + p.pos) / ${NB_NUCLEOTIDES_PER_CLUSTER}) AS INT64) AS clustering_index,
        CAST(p.chr AS INT64) AS chr,
        p.* EXCEPT(chr)
        FROM ${PROJECT}.${CLOUDASM_DATASET}.${sample}_${table} p
        JOIN ${REFG_FOLDER}.chr_length c
        ON SAFE_CAST(p.chr AS INT64) = c.chr
        WHERE REGEXP_CONTAINS(p.chr, r'^\\d+$')
    "

bq query \
    --use_legacy_sql=false \
    --destination_table="${TEMP_TABLE_2}" \
    --replace=true \
    --clustering_fields=sample,chr \
    --range_partitioning=clustering_index,0,4000,1 \
    "
    SELECT c.region_inf, c.region_sup, c.region_nb_cpg, p.*
    FROM ${TEMP_TABLE} p
    INNER JOIN ${REFG_FOLDER}.regions_w_cpg_no_blacklist c
    ON p.chr = c.chr AND p.clustering_index = c.clustering_index AND pos >= region_inf AND pos < region_sup
    "

# Append table
bq cp --append_table \
    "${TEMP_TABLE_2}" \
    "${PROJECT}":"${SAMPLES_DATASET}"."${table}"

# Delete temp table
bq rm -f -t "${TEMP_TABLE}"
bq rm -f -t "${TEMP_TABLE_2}"