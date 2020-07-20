#!/bin/bash

# Import the file where the p-value of each cpg
# is calculated.
bq --location=US load \
               --autodetect \
               --replace=true \
               --source_format=NEWLINE_DELIMITED_JSON \
                ${DATASET_ID}.${SAMPLE}_asm_region_pvalue \
               gs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_asm_region_pvalue.json 

# Query to select the SNPs with at least 3 significant CpGs in the same direction
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_asm_snp \
    --replace=true \
    "
    SELECT
        chr,
        region_inf,
        region_sup,
        snp_id,
        snp_pos,
        IF (
            wilcoxon_corr_pvalue < ${P_VALUE} 
            AND (
                (pos_sig_cpg >= ${CPG_SAME_DIRECTION_ASM} AND nb_consec_pos_sig_asm >= ${CONSECUTIVE_CPG} AND asm_region_effect > ${ASM_REGION_EFFECT})
                OR (neg_sig_cpg >= ${CPG_SAME_DIRECTION_ASM} AND nb_consec_neg_sig_asm >= ${CONSECUTIVE_CPG} AND asm_region_effect < -${ASM_REGION_EFFECT})
                ), TRUE, FALSE) AS asm_snp,
        asm_region_inf,
        asm_region_sup,
        ref_reads AS nb_ref_reads,
        alt_reads AS nb_alt_reads,
        asm_region_effect,
        wilcoxon_corr_pvalue,
        region_nb_cpg,
        nb_cpg,
        nb_sig_cpg,
        pos_sig_cpg AS nb_pos_sig_cpg,
        neg_sig_cpg AS nb_neg_sig_cpg,
        nb_consec_pos_sig_asm,
        nb_consec_neg_sig_asm
    FROM ${DATASET_ID}.${SAMPLE}_asm_region_pvalue
    "
