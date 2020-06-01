#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_OUT}.asm_read_cpg_motifs_${CHR} \
    --replace=true \
    "
    WITH 
    ASM AS (
        SELECT
            sample,
            snp_id,
            region_inf,
            region_sup,
            chr
        FROM ${DATASET_OUT}.asm_read_cpg_arrays
        WHERE chr = '${CHR}'
    ),
    MOTIF_DB AS (
        SELECT 
            chr AS chr_motif, 
            motif_start, 
            motif_end, 
            motif
        FROM ${DATASET_EPI}.asm_motifs_locations
        WHERE chr = '${CHR}'
    ),
    COMBINED AS (
        SELECT * FROM ASM 
        LEFT JOIN MOTIF_DB
        ON
            motif_end <= region_sup + ${EPI_REGION} AND
            motif_start >= region_inf - ${EPI_REGION} AND
            chr = chr_motif
        )
    SELECT * EXCEPT(motif_start, motif_end, chr_motif) FROM COMBINED
    "
