#!/bin/bash

echo "Create a table of CpG array per (genomic window, snp id) combination"
# We ask that there are the minimum number of CpGs per combination
# (specified as a variable)
bq query \
    --use_legacy_sql=false \
    --destination_table "${SAMPLES_DATASET}".all_cpgs_w_asm_computed \
    --replace=true \
    --clustering_fields=sample,chr,clustering_index \
    --range_partitioning=clustering_index,0,4000,1 \
    "
    WITH CONTEXT_ASM AS (
        SELECT *
        FROM ${SAMPLES_DATASET}.cpg_asm
        -- Note: we do not require that there is a min coverage because that was done by CloudASM
        WHERE (ref_cov + alt_cov < ${MAX_CPG_COV}) AND (ref_cov + alt_cov >= ${MIN_CPG_COV})
        ),
        GROUP_BY_SNPID AS (
            SELECT
                sample,
                clustering_index,
                chr,
                region_inf,
                region_sup,
                region_nb_cpg,
                snp_id,
                snp_pos,
                COUNT(*) AS nb_cpg_found,
                SUM(IF(fisher_pvalue < ${MAX_P_VALUE}, COUNT(*), 0)) AS nb_sig_cpg
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
                    ) AS cpg_asm
            FROM CONTEXT_ASM
            GROUP BY sample, clustering_index, chr, snp_id, snp_pos, region_inf, region_sup, region_nb_cpg
            HAVING nb_cpg_found >= ${MIN_NB_CPG_FOR_ASM}
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
            nb_cpg_found,
            (SELECT COUNT(fisher_pvalue) FROM UNNEST(cpg_asm) WHERE fisher_pvalue < ${MAX_P_VALUE}) AS nb_sig_cpg,
            (SELECT COUNT(fisher_pvalue) FROM UNNEST(cpg_asm) WHERE fisher_pvalue < ${MAX_P_VALUE} AND SIGN(effect) = 1) AS pos_sig_cpg,
            (SELECT COUNT(fisher_pvalue) FROM UNNEST(cpg_asm) WHERE fisher_pvalue < ${MAX_P_VALUE} AND SIGN(effect) = -1) AS neg_sig_cpg,
            (SELECT MIN(pos) FROM UNNEST(cpg_asm) WHERE fisher_pvalue < ${MAX_P_VALUE}) AS min_pos_sig_cpg,
            (SELECT MAX(pos) FROM UNNEST(cpg_asm) WHERE fisher_pvalue < ${MAX_P_VALUE}) AS max_pos_sig_cpg,
            (SELECT MIN(pos) FROM UNNEST(cpg_asm)) AS min_cpg_pos,
            (SELECT MAX(pos) FROM UNNEST(cpg_asm)) AS max_cpg_pos,
            cpg_asm
        FROM GROUP_BY_SNPID
    "

 
echo "Create a list of CpGs - reads ID relevant to the windows where ASM needs to be computed"
# Meaning that CpGs need to be within 2 significant CpGs (that can have different ASM directions) 
# or the two extreme CpGs (but in that case, there won't be ASM)
bq query \
    --use_legacy_sql=false \
    --destination_table "${SAMPLES_DATASET}".cpg_for_asm_region_effect \
    --replace=true \
    --clustering_fields=sample,chr,clustering_index \
    "
    WITH 
        ALL_CPG AS (
            SELECT *
            FROM ${SAMPLES_DATASET}.cpg_read_genotype
        ),
        -- Import all CpGs that are at least 5x covered on each allele and for which we have a fisher p-value
        WELL_COVERED_CPG AS (
            SELECT *
            FROM ${SAMPLES_DATASET}.all_cpgs_w_asm_computed
        )
        SELECT DISTINCT
            c.sample,
            c.clustering_index,
            c.chr,
            c.region_inf,
            c.region_sup,
            pos,
            meth,
            cov,
            allele,
            read_id,
            c.snp_id,
            c.snp_pos,
            min_cpg_pos,
            max_cpg_pos,
            min_pos_sig_cpg,
            max_pos_sig_cpg,
            nb_cpg_found,
            nb_sig_cpg,
            pos_sig_cpg,
            neg_sig_cpg
        FROM ALL_CPG c
        INNER JOIN WELL_COVERED_CPG p
        ON
            c.sample = p.sample AND
            c.chr = p.chr AND
            c.clustering_index = p.clustering_index AND
            c.snp_id = p.snp_id
        WHERE
            pos >= IF(min_pos_sig_cpg IS NULL, min_cpg_pos, min_pos_sig_cpg) AND
            pos <= IF(max_pos_sig_cpg IS NULL, max_cpg_pos, max_pos_sig_cpg)
        "


bq query \
    --use_legacy_sql=false \
    --destination_table "${SAMPLES_DATASET}".reads_asm \
    --replace=true \
    --clustering_fields=sample,chr,clustering_index \
    "
    WITH 
        QUALIFYING_CPG_WILCOX AS (
            SELECT * FROM ${SAMPLES_DATASET}.cpg_for_asm_region_effect
        ),
        METHYL_PER_READ_WILCOX AS (
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
            FROM QUALIFYING_CPG_WILCOX
            GROUP BY sample, chr, clustering_index, snp_id, snp_pos, read_id, allele, region_inf, region_sup
        ),
        SNP_METHYL_ARRAY_REF_WILCOX AS (
            SELECT
                sample,
                chr,
                clustering_index,
                snp_id,
                snp_pos,
                region_inf,
                region_sup,
                ARRAY_AGG(STRUCT(methyl_perc)) AS ref
            FROM METHYL_PER_READ_WILCOX
            WHERE allele = 'REF'
            GROUP BY sample, chr, clustering_index, snp_id, snp_pos, region_inf, region_sup
        ),
        SNP_METHYL_ARRAY_ALT_WILCOX AS (
            SELECT
                sample,
                chr,
                clustering_index,
                snp_id,
                snp_pos,
                region_inf,
                region_sup,
                ARRAY_AGG(STRUCT(methyl_perc)) AS alt
            FROM METHYL_PER_READ_WILCOX
            WHERE allele = 'ALT'
            GROUP BY sample, chr, clustering_index, snp_id, snp_pos, region_inf, region_sup
        ),
        SNP_METHYL_JOIN_WILCOX AS (
            SELECT
                c.sample,
                c.chr,
                c.clustering_index,
                c.snp_id,
                c.snp_pos,
                c.region_inf,
                c.region_sup,
                ref,
                alt
            FROM SNP_METHYL_ARRAY_REF_WILCOX c
            INNER JOIN SNP_METHYL_ARRAY_ALT_WILCOX p
            ON
               c.sample = p.sample AND
               c.clustering_index = p.clustering_index AND
               c.snp_id = p.snp_id AND
               c.snp_pos = p.snp_pos AND
               c.chr = p.chr AND
               c.region_inf = p.region_inf AND
               c.region_sup = p.region_sup
        )
        SELECT
            sample,
            chr,
            clustering_index,
            snp_id,
            snp_pos,
            region_inf,
            region_sup,
            ROUND(((SELECT AVG(methyl_perc) FROM UNNEST(alt)) - (SELECT AVG(methyl_perc) FROM UNNEST(ref))),5) AS asm_region_effect,
            ARRAY_LENGTH(ref) AS ref_reads, 
            ARRAY_LENGTH(alt) AS alt_reads,
            ref, 
            alt
        FROM SNP_METHYL_JOIN_WILCOX
    "


echo "Merge the 2 tables of CpG fractional methyl and read fractional methyl arrays."
bq query \
    --use_legacy_sql=false \
    --destination_table "${SAMPLES_DATASET}".cpg_reads_asm \
    --replace=true \
    --clustering_fields=sample,chr,clustering_index \
    "
    SELECT
        t1.sample,
        t1.chr,
        t1.clustering_index,
        t1.region_inf,
        t1.region_sup,
        t1.snp_id,
        t1.snp_pos,
        t1.region_nb_cpg,
        t1.nb_cpg,
        t1.nb_sig_cpg,
        t1.pos_sig_cpg,
        t1.neg_sig_cpg,
        t1.min_cpg,
        t1.max_cpg,
        t1.asm_region_inf,
        t1.asm_region_sup,
        t2.asm_region_effect,
        t2.ref_reads,
        t2.alt_reads,
        t1.cpg_asm As cpg,
        t2.ref,
        t2.alt
    FROM ${SAMPLES_DATASET}.cpg_asm_over_regions t1
    JOIN ${SAMPLES_DATASET}.reads_asm t2
    ON 
        t1.sample = t2.sample AND
        t1.clustering_index = t2.clustering_index AND
        t1.chr = t2.chr AND 
        t1.region_inf = t2.region_inf AND 
        t1.region_sup = t2.region_sup AND
        t1.snp_id = t2.snp_id AND
        t1.snp_pos = t2.snp_pos
    "

# # Export file to JSON format in the bucket
# # (nested arrays are not supported in)
# bq extract \
#     --destination_format NEWLINE_DELIMITED_JSON \
#     "${CLOUDASM_STANDARD_REGIONS_DATASET}".cpg_reads_asm \
#     gs://"${BUCKET_NAME}"/"${CLOUDASM_STANDARD_REGIONS_DATASET}"/cpg_reads_asm*.json