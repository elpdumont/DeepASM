#!/bin/bash

# This query overlaps regions and snps (and associated cpgs and reads)

bq query \
    --use_legacy_sql=false \
    --destination_table "${SAMPLES_DATASET}".regions_and_snps \
    --replace=true \
    --clustering_fields=sample,chr \
    --range_partitioning=clustering_index,0,4000,1 \
    "
    WITH CPG_WITH_ASM_TMP AS(
        SELECT *
        FROM ${SAMPLES_DATASET}.cpg_asm
        WHERE (ref_cov + alt_cov < ${MAX_CPG_COV}) AND (ref_cov + alt_cov >= ${MIN_CPG_COV})
    ),
    CPG_WITH_ASM AS (
        SELECT
            sample,
            clustering_index,
            chr,
            snp_id,
            snp_pos,
            region_inf,
            region_sup,
            region_nb_cpg,
            COUNT(*) AS nb_cpg_found_w_snp,
            ARRAY_AGG(
                STRUCT(
                    pos,
                    ref_cov,
                    ref_meth,
                    alt_cov,
                    alt_meth,
                    ROUND(alt_meth/alt_cov-ref_meth/ref_cov,5) AS effect,
                    fisher_pvalue
                    )
                    ORDER BY pos
                ) AS cpg_w_snp
        FROM CPG_WITH_ASM_TMP
        GROUP BY sample, clustering_index, chr, region_inf, region_sup, region_nb_cpg, snp_id, snp_pos
        HAVING nb_cpg_found_w_snp >= ${MIN_NB_CPG_FOR_ASM}
    ),
    READS_WITH_ASM_TMP AS (
        SELECT
            sample,
            chr,
            clustering_index,
            snp_id,
            snp_pos,
            region_inf,
            region_sup,
            read_id,
            allele,
            ROUND(SAFE_DIVIDE(SUM(meth),SUM(cov)),5) AS methyl_perc,
        FROM ${SAMPLES_DATASET}.cpg_read_genotype
        GROUP BY sample, chr, clustering_index, snp_id, snp_pos, read_id, allele, region_inf, region_sup
    ),
    READS_WITH_ASM AS (
        SELECT
            sample,
            chr,
            clustering_index,
            snp_id,
            snp_pos,
            region_inf,
            region_sup,
            COUNT(*) AS nb_reads_found_w_snp,
            ARRAY_AGG(
                STRUCT(
                    CASE WHEN allele = 'REF' THEN methyl_perc END AS ref,
                    CASE WHEN allele = 'ALT' THEN methyl_perc END AS alt
                )
            ) AS reads_w_snp
        FROM READS_WITH_ASM_TMP
        GROUP BY sample, chr, clustering_index, snp_id, snp_pos, region_inf, region_sup
    )
    SELECT
        p.sample,
        p.chr,
        p.clustering_index,
        p.snp_id,
        p.snp_pos,
        p.region_inf,
        p.region_sup,
        p.nb_cpg_found_w_snp,
        c.nb_reads_found_w_snp,
        p.cpg_w_snp,
        c.reads_w_snp
    FROM CPG_WITH_ASM p
    INNER JOIN READS_WITH_ASM c
    ON p.sample = c.sample AND p.chr = c.chr AND p.clustering_index = c. clustering_index AND p.region_inf = c.region_inf AND p.region_sup = c.region_sup AND p.snp_id = c.snp_id AND p.snp_pos = c.snp_pos
    "