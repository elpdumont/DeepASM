#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_cpg_asm_${GENOMIC_INTERVAL}bp_${CHR} \
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
                dnase,
                encode_ChiP_V2,
                tf_motifs,
                pos, 
                snp_id, 
                ref_cov, 
                ref_meth, 
                alt_cov, 
                alt_meth, 
                fisher_pvalue
            FROM CPG_REGIONS
            INNER JOIN CONTEXT_ASM
            ON pos >= region_inf AND pos < region_sup
        ),
        GROUP_BY_SNPID AS (
            SELECT
                chr,
                region_inf,
                region_sup,
                snp_id,
                region_nb_cpg,
                dnase,
                encode_ChiP_V2,
                tf_motifs,
                COUNT(*) AS nb_cpg_found,
                ARRAY_AGG(
                    STRUCT(
                        pos, 
                        ref_cov, 
                        ref_meth, 
                        alt_cov, 
                        alt_meth, 
                        fisher_pvalue
                        )
                    ) AS cpg_asm
            FROM REGION_CPG_JOINED
            GROUP BY chr, snp_id, region_inf, region_sup, region_nb_cpg, dnase, encode_ChiP_V2, tf_motifs
        )
        SELECT * 
        FROM GROUP_BY_SNPID
        WHERE nb_cpg_found >= 3
    "