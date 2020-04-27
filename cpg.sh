#!/bin/bash

# Get the fractional methylation and coverage of each CpG
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_OUT}.cpgs \
    --append_table \
    "
    SELECT
        '${SAMPLE}' AS sample,
        snp_id,
        ANY_VALUE(snp_pos) AS snp_pos,
        ANY_VALUE(chr) AS chr,
        ARRAY_AGG(
            STRUCT(
                pos,
                ROUND((alt_meth+ref_meth)/(ref_cov+alt_cov),3) AS frac_methyl,
                ref_cov + alt_cov AS cov
            )
            ORDER BY pos
        ) AS cpg
    FROM ${DATASET_IN}.${SAMPLE}_cpg_asm
    GROUP BY snp_id
    "


# Previously used
                        # ref_meth + alt_meth as meth,
                        # ROUND(alt_meth/alt_cov-ref_meth/ref_cov,3) AS effect,
                        # fisher_pvalue,
                        # ref_cov,
                        # ref_meth,
                        # alt_cov,
                        # alt_meth