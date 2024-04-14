#!/bin/bash

echo "Create a table of CpG array per (genomic window, snp id) combination"
# We ask that there are the minimum number of CpGs per combination
# (specified as a variable)
bq query \
    --use_legacy_sql=false \
    --destination_table "${CLOUDASM_STANDARD_REGIONS_DATASET}".cpg_asm_over_regions \
    --replace=true \
    --clustering_fields=sample,chr,clustering_index \
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
                        ROUND(alt_meth/alt_cov-ref_meth/ref_cov,5) AS effect,
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
            region_center,
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


echo "Create a list of CpGs - reads ID relevant to the windows where ASM needs to be computed"
# Meaning that CpGs need to be within 2 significant CpGs (that can have different ASM directions) 
# or the two extreme CpGs (but in that case, there won't be ASM)
bq query \
    --use_legacy_sql=false \
    --destination_table "${CLOUDASM_STANDARD_REGIONS_DATASET}".cpg_for_asm_region_effect \
    --replace=true \
    --clustering_fields=sample,chr,clustering_index \
    "
    WITH 
        -- Import the list of CpGs with their genotype as a function of read_id and snp_id
        ALL_CPG AS (
            SELECT
                sample,
                clustering_index,
                chr,
                pos,
                meth,
                cov,
                snp_id,
                allele,
                read_id 
            FROM ${CLOUDASM_STANDARD_REGIONS_DATASET}.cpg_read_genotype
        ),
        -- Import all CpGs that are at least 5x covered on each allele and for which we have a fisher p-value
        WELL_COVERED_CPG_ARRAY AS (
            SELECT 
                sample,
                clustering_index,
                chr,
                region_inf,
                region_sup,
                region_center,
                snp_id,
                snp_pos,
                asm_region_inf,
                asm_region_sup,
                min_cpg,
                max_cpg,
                nb_cpg,
                nb_sig_cpg,
                pos_sig_cpg,
                neg_sig_cpg,
                ARRAY(
                    (SELECT pos FROM UNNEST(cpg_asm))
                    ) AS pos_array
            FROM ${CLOUDASM_STANDARD_REGIONS_DATASET}.cpg_asm_over_regions
        ),
        WELL_COVERED_CPG AS (
            SELECT
                sample,
                clustering_index,
                chr,
                region_inf,
                region_sup,
                region_center,
                snp_id,
                snp_pos,
                asm_region_inf,
                asm_region_sup,
                min_cpg,
                max_cpg, 
                nb_cpg,
                nb_sig_cpg,
                pos_sig_cpg,
                neg_sig_cpg,
                pos_well_cov
            FROM WELL_COVERED_CPG_ARRAY, UNNEST(pos_array) AS pos_well_cov
        )
        -- Keep the combination of CpG, read_id, allele for which CpGs are at least 5x covered on each allele
            SELECT DISTINCT
                c.sample,
                c.clustering_index,
                c.chr,
                region_inf,
                region_sup,
                pos,
                meth,
                cov,
                allele,
                read_id,
                c.snp_id,
                snp_pos,
                asm_region_inf,
                asm_region_sup,
                min_cpg,
                max_cpg,
                nb_cpg,
                nb_sig_cpg,
                pos_sig_cpg,
                neg_sig_cpg
            FROM ALL_CPG c
            INNER JOIN WELL_COVERED_CPG p
            ON 
                c.sample = p.sample AND
                c.chr = p.chr AND
                c.clustering_index = p.clustering_index AND
                --c.pos = pos_well_cov AND
                c.snp_id = p.snp_id
            WHERE 
                pos >= IF(asm_region_inf IS NULL, min_cpg, asm_region_inf) AND
                pos <= IF(asm_region_sup IS NULL, max_cpg, asm_region_sup)
        "


bq query \
    --use_legacy_sql=false \
    --destination_table "${CLOUDASM_STANDARD_REGIONS_DATASET}".reads_asm \
    --replace=true \
    --clustering_fields=sample,chr,clustering_index \
    "
    WITH 
        QUALIFYING_CPG_WILCOX AS (
            SELECT * FROM ${CLOUDASM_STANDARD_REGIONS_DATASET}.cpg_for_asm_region_effect
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
    --destination_table "${CLOUDASM_STANDARD_REGIONS_DATASET}".cpg_reads_asm \
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
    FROM ${CLOUDASM_STANDARD_REGIONS_DATASET}.cpg_asm_over_regions t1
    JOIN ${CLOUDASM_STANDARD_REGIONS_DATASET}.reads_asm t2 
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