#!/bin/bash

# Get the fractional methylation and coverage of each read
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_OUT}.reads \
    --append_table \
    "
    SELECT
        '${SAMPLE}' AS sample,
        snp_id,
        snp_pos, 
        chr, 
        alt_reads, 
        ref_reads, 
        ref, 
        alt
    FROM ${DATASET_IN}.${SAMPLE}_asm_region_pvalue
    "

# Previously used
                        # ref_meth + alt_meth as meth,
                        # ROUND(alt_meth/alt_cov-ref_meth/ref_cov,3) AS effect,
                        # fisher_pvalue,
                        # ref_cov,
                        # ref_meth,
                        # alt_cov,
                        # alt_meth