#!/bin/bash


# Create a list of CpGs - reads ID relevant to the windows where ASM needs to be computed
# Meaning that CpGs need to be within 2 significant CpGs (that can have different ASM directions) 
# or the two extreme CpGs (but in that case, there won't be ASM)
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_cpg_for_asm_region_effect \
    --replace=true \
    "
    WITH 
        -- Import the list of CpGs with their genotype as a function of read_id and snp_id
        ALL_CPG AS (
            SELECT 
                chr AS chr_cpg,
                pos AS pos_cpg,
                meth,
                cov,
                snp_id AS snp_id_cpg,
                allele,
                read_id 
            FROM ${DATASET_CLOUDASM}.${SAMPLE}_cpg_read_genotype
        ),
        -- Import all CpGs that are at least 5x covered on each allele and for which we have a fisher p-value
        WELL_COVERED_CPG_ARRAY AS (
            SELECT 
                chr_asm_region,
                region_inf,
                region_sup,
                snp_id_asm_region,
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
            FROM ${DATASET_PRED}.${SAMPLE}_cpg_asm
        ),
        WELL_COVERED_CPG AS (
            SELECT chr_asm_region,
                region_inf,
                region_sup,
                snp_id_asm_region,
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
                chr_cpg,
                pos_cpg,
                meth,
                cov,
                allele,
                read_id,
                region_inf,
                region_sup,
                snp_id_asm_region AS snp_id,
                snp_pos,
                asm_region_inf,
                asm_region_sup,
                min_cpg,
                max_cpg,
                nb_cpg,
                nb_sig_cpg,
                pos_sig_cpg,
                neg_sig_cpg
            FROM ALL_CPG
            INNER JOIN WELL_COVERED_CPG
            ON 
                chr_cpg = chr_asm_region AND 
                pos_cpg = pos_well_cov AND
                snp_id_cpg = snp_id_asm_region
            WHERE 
                pos_cpg >= IF(asm_region_inf IS NULL, min_cpg, asm_region_inf) AND
                pos_cpg <= IF(asm_region_sup IS NULL, max_cpg, asm_region_sup)
        "


bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_reads_asm \
    --replace=true \
    "
    WITH 
        QUALIFYING_CPG_WILCOX AS (
            SELECT * FROM ${DATASET_PRED}.${SAMPLE}_cpg_for_asm_region_effect
        ),
        METHYL_PER_READ_WILCOX AS (
            SELECT 
                snp_id,
                snp_pos,
                chr_cpg,
                ANY_VALUE(asm_region_inf) AS asm_region_inf,
                ANY_VALUE(asm_region_sup) AS asm_region_sup,
                read_id,
                allele,
                ROUND(SAFE_DIVIDE(SUM(meth),SUM(cov)),5) AS methyl_perc,
                ANY_VALUE(nb_cpg) AS nb_cpg,
                ANY_VALUE(nb_sig_cpg) AS nb_sig_cpg,
                ANY_VALUE(pos_sig_cpg) AS pos_sig_cpg,
                ANY_VALUE(neg_sig_cpg) AS neg_sig_cpg
            FROM QUALIFYING_CPG_WILCOX
            GROUP BY snp_id, snp_pos, read_id, allele, chr_cpg
        ),
        SNP_METHYL_ARRAY_REF_WILCOX AS (
            SELECT
                snp_id,
                ANY_VALUE(snp_pos) AS snp_pos,
                ANY_VALUE(chr_cpg) AS chr,
                ANY_VALUE(asm_region_inf) AS asm_region_inf,
                ANY_VALUE(asm_region_sup) AS asm_region_sup,
                ARRAY_AGG(STRUCT(methyl_perc)) AS ref,
                ANY_VALUE(nb_cpg) AS nb_cpg,
                ANY_VALUE(nb_sig_cpg) AS nb_sig_cpg,
                ANY_VALUE(pos_sig_cpg) AS pos_sig_cpg,
                ANY_VALUE(neg_sig_cpg) AS neg_sig_cpg
            FROM METHYL_PER_READ_WILCOX
            WHERE allele = 'REF'
            GROUP BY snp_id
        ),
        SNP_METHYL_ARRAY_ALT_WILCOX AS (
            SELECT
                snp_id AS snp_id_alt,
                ARRAY_AGG(STRUCT(methyl_perc)) AS alt
            FROM METHYL_PER_READ_WILCOX
            WHERE allele = 'ALT'
            GROUP BY snp_id
        ),
        SNP_METHYL_JOIN_WILCOX AS (
            SELECT * FROM SNP_METHYL_ARRAY_REF_WILCOX
            INNER JOIN SNP_METHYL_ARRAY_ALT_WILCOX
            ON snp_id = snp_id_alt
        ),
        SNP_METHYL_EFFECT AS (
            SELECT
                snp_id AS snp_id_effect,
                ROUND(((SELECT AVG(methyl_perc) FROM UNNEST(alt)) - (SELECT AVG(methyl_perc) FROM UNNEST(ref))),3) AS effect
            FROM SNP_METHYL_JOIN_WILCOX
        )
        SELECT 
            snp_id,
            snp_pos,
            chr,
            asm_region_inf,
            asm_region_sup,
            effect AS asm_region_effect,
            ARRAY_LENGTH(ref) AS ref_reads, 
            ARRAY_LENGTH(alt) AS alt_reads,
            ref, 
            alt,
            nb_cpg,
            nb_sig_cpg,
            pos_sig_cpg,
            neg_sig_cpg
        FROM SNP_METHYL_EFFECT
        INNER JOIN SNP_METHYL_JOIN_WILCOX
        ON snp_id = snp_id_effect
        "


# Delete intermediary files
bq rm -f -t ${DATASET_ID}.${SAMPLE}_snp_cpg_array
bq rm -f -t ${DATASET_ID}.${SAMPLE}_cpg_for_asm_region_effect
bq rm -f -t ${DATASET_ID}.${SAMPLE}_cpg_for_wilcox

# Export file to JSON format in the bucket
# (nested arrays are not supported in)
bq extract \
    --destination_format NEWLINE_DELIMITED_JSON \
    ${DATASET_ID}.${SAMPLE}_snp_for_asm_region \
    gs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_snp_for_asm_region.json