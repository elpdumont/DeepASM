#!/bin/bash


# Query to select the SNPs with at least 3 significant CpGs in the same direction
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_all_info_${GENOMIC_INTERVAL}bp \
    --replace=true \
    "
    SELECT
        t1.sample,
        t1.sample_category,
        t1.chr,
        t1.region_inf,
        t1.region_sup,
        t1.snp_id,
        t1.snp_pos,
        t1.asm_snp,
        t1.wilcoxon_corr_pvalue,
        t1.asm_region_effect,
        t1.region_nb_cpg,
        t1.nb_cpg_found,
        t1.nb_reads,
        t1.dnase,
        t1.encode_ChiP_V2,
        t1.tf_motifs,
        t1.global_cpg_fm,
        t1.tot_nb_cpg,
        t1.tot_nb_reads,
        t1.cpg,
        t1.read,
        t2.window_details
    FROM ${DATASET_PRED}.${SAMPLE}_cpg_read_asm_${GENOMIC_INTERVAL}bp t1
    LEFT JOIN ${DATASET_PRED}.${SAMPLE}_genomic_window_${GENOMIC_INTERVAL}bp t2
    ON t1.chr = t2.chr AND t1.region_inf = t2.region_inf AND t1.region_sup = t2.region_sup
    "