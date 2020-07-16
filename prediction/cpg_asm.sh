#!/bin/bash

# Create a table of CpG array per (genomic window, snp id) combination
# We ask that there are the minimum number of CpGs per combination
# (specified as a variable)
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_cpg_asm_${CHR} \
    --replace=true \
    "
    WITH 
        CPG_REGIONS AS (
            SELECT 
                chr AS chr_region, 
                region_inf, 
                region_sup, 
                region_nb_cpg
            FROM ${DATASET_EPI}.hg19_cpg_regions_${GENOMIC_INTERVAL}bp_annotated
            WHERE chr = '${CHR}'
        ),
        CONTEXT_ASM AS (
            SELECT *
            FROM ${DATASET_CONTEXT}.${SAMPLE}_cpg_asm
            WHERE chr = '${CHR}'
        ),
        REGION_CPG_JOINED AS(
            SELECT 
                chr,
                region_inf,
                region_sup,
                region_nb_cpg,
                pos, 
                snp_id, 
                snp_pos,
                ref_cov, 
                ref_meth, 
                alt_cov, 
                alt_meth, 
                fisher_pvalue
            FROM CPG_REGIONS
            INNER JOIN CONTEXT_ASM
            ON pos >= region_inf AND pos <= region_sup
        ),
        GROUP_BY_SNPID AS (
            SELECT
                chr,
                region_inf,
                region_sup,
                region_nb_cpg,
                snp_id,
                snp_pos,
                COUNT(*) AS nb_cpg,
                ARRAY_AGG(
                    STRUCT(
                        pos, 
                        ref_cov, 
                        ref_meth, 
                        alt_cov, 
                        alt_meth, 
                        ROUND(alt_meth/alt_cov-ref_meth/ref_cov,3) AS effect,
                        fisher_pvalue
                        )
                        ORDER BY pos
                    ) AS cpg_asm
            FROM REGION_CPG_JOINED
            GROUP BY chr, snp_id, snp_pos, region_inf, region_sup, region_nb_cpg
        ),
        GROUP_BY_SNPID_MIN_CPG AS (
        SELECT * 
        FROM GROUP_BY_SNPID
        WHERE nb_cpg >= ${CPG_PER_ASM_REGION}
        )
        SELECT 
            chr AS chr_asm_region,
            region_inf,
            region_sup,
            region_nb_cpg,
            snp_id AS snp_id_asm_region, 
            snp_pos,
            nb_cpg,
            (SELECT COUNT(fisher_pvalue) FROM UNNEST(cpg_asm) WHERE fisher_pvalue < ${P_VALUE}) AS nb_sig_cpg, 
            (SELECT COUNT(fisher_pvalue) FROM UNNEST(cpg_asm) WHERE fisher_pvalue < ${P_VALUE} AND SIGN(effect) = 1) AS pos_sig_cpg,
            (SELECT COUNT(fisher_pvalue) FROM UNNEST(cpg_asm) WHERE fisher_pvalue < ${P_VALUE} AND SIGN(effect) = -1) AS neg_sig_cpg, 
            (SELECT MIN(pos) FROM UNNEST(cpg_asm) WHERE fisher_pvalue < ${P_VALUE}) AS asm_region_inf,
            (SELECT MAX(pos) FROM UNNEST(cpg_asm) WHERE fisher_pvalue < ${P_VALUE}) AS asm_region_sup,
            (SELECT MIN(pos) FROM UNNEST(cpg_asm)) AS min_cpg,
            (SELECT MAX(pos) FROM UNNEST(cpg_asm)) AS max_cpg,
            cpg_asm
        FROM GROUP_BY_SNPID_MIN_CPG
        ORDER BY region_inf
    "