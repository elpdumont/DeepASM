#!/bin/bash

# Create a table of CpG array per (genomic window, snp id) combination
# We ask that there are the minimum number of CpGs per combination
# (specified as a variable)
bq query \
    --use_legacy_sql=false \
    --destination_table "${CLOUDASM_STANDARD_REGIONS_DATASET}".cpg_asm_over_regions \
    --replace=true \
    "
    WITH 
        CPG_REGIONS AS (
            SELECT 
                chr,
                region_inf, 
                region_sup,
                region_center,
                region_nb_cpg,
                clustering_index
            FROM ${REFG_FOLDER}.regions_w_cpg_wo_blacklisted_regions
        ),
        CONTEXT_ASM AS (
            SELECT *
            FROM ${CLOUDASM_STANDARD_REGIONS_DATASET}.cpg_asm
            -- Note: we do not require that there is a min coverage because that was done by CloudASM
            WHERE (ref_cov + alt_cov < ${MAX_CPG_COV})
        ),
        REGION_CPG_JOINED AS(
            SELECT 
                sample,
                c.clustering_index,
                c.chr,
                region_inf,
                region_sup,
                region_center,
                region_nb_cpg,
                pos, 
                snp_id, 
                snp_pos,
                ref_cov, 
                ref_meth, 
                alt_cov, 
                alt_meth, 
                fisher_pvalue
            FROM CPG_REGIONS p
            INNER JOIN CONTEXT_ASM c
            ON pos >= region_inf AND pos <= region_sup AND c.chr = p.chr AND c.clustering_index = p.clustering_index
        ),
        GROUP_BY_SNPID AS (
            SELECT
                sample,
                clustering_index,
                chr,
                region_inf,
                region_sup,
                region_center,
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
            GROUP BY sample, clustering_index, chr, snp_id, snp_pos, region_inf, region_sup, region_center, region_nb_cpg
        ),
        GROUP_BY_SNPID_MIN_CPG AS (
        SELECT * 
        FROM GROUP_BY_SNPID
        WHERE nb_cpg >= ${MIN_NB_CPG_FOR_ASM}
        )
        SELECT
            sample,
            clustering_index,
            chr,
            region_inf,
            region_sup,
            region_nb_cpg,
            snp_id,
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
    "