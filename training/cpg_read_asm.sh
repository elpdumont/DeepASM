#!/bin/bash


# Merge the 2 tables of CpG fractional methyl and read fractional methyl arrays.
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_cpg_reads_asm \
    --replace=true \
    "
    SELECT
        t1.chr_asm_region AS chr,
        t1.region_inf AS region_inf,
        t1.region_sup AS region_sup,
        t1.snp_id_asm_region AS snp_id,
        t1.snp_pos AS snp_pos,
        t1.region_nb_cpg AS region_nb_cpg,
        t1.nb_cpg AS nb_cpg,
        t1.nb_sig_cpg AS nb_sig_cpg,
        t1.pos_sig_cpg AS pos_sig_cpg,
        t1.neg_sig_cpg AS neg_sig_cpg,
        t1.min_cpg AS min_cpg,
        t1.max_cpg AS max_cpg,
        t1.asm_region_inf AS asm_region_inf,
        t1.asm_region_sup AS asm_region_sup,
        t2.asm_region_effect AS asm_region_effect,
        t2.ref_reads AS ref_reads,
        t2.alt_reads AS alt_reads,
        t1.cpg_asm As cpg,
        t2.ref,
        t2.alt
    FROM ${DATASET_PRED}.${SAMPLE}_cpg_asm t1
    JOIN ${DATASET_PRED}.${SAMPLE}_reads_asm t2 
    ON 
        t1.chr_asm_region = t2.chr AND 
        t1.region_inf = t2.region_inf AND 
        t1.region_sup = t2.region_sup AND
        t1.snp_id_asm_region = t2.snp_id AND
        t1.snp_pos = t2.snp_pos
    "

# Export file to JSON format in the bucket
# (nested arrays are not supported in)
bq extract \
    --destination_format NEWLINE_DELIMITED_JSON \
    ${DATASET_PRED}.${SAMPLE}_cpg_reads_asm \
    gs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_snp_for_asm_region.json