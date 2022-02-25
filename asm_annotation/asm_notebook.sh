#!/bin/bash


# Query to select the SNPs with at least 3 significant CpGs in the same direction
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_cpg_read_asm_${GENOMIC_INTERVAL}bp \
    --replace=true \
    "
    WITH 
    TOTAL_CPG AS (
        SELECT 
            *,
            (
                SELECT SUM(nb_cpg_found) 
                FROM ${DATASET_PRED}.${SAMPLE}_cpg_read_${GENOMIC_INTERVAL}bp
            ) AS tot_nb_cpg,
            (
                SELECT SUM(nb_reads) 
                FROM ${DATASET_PRED}.${SAMPLE}_cpg_read_${GENOMIC_INTERVAL}bp
            ) AS tot_nb_reads,
            (
                SELECT SUM(fm) 
                FROM UNNEST(cpg)
            ) AS tot_region_cpg_fm, 
        FROM ${DATASET_PRED}.${SAMPLE}_cpg_read_${GENOMIC_INTERVAL}bp
    ),
    TOTAL_FRAC_METHYL AS (
        SELECT
            *,
            (
            SELECT SUM(tot_region_cpg_fm) 
            FROM TOTAL_CPG
        ) AS tot_cpg_fm,
        FROM TOTAL_CPG
    )
    SELECT
        '${SAMPLE}' AS sample,
        '${SAMPLE_STATUS}' AS sample_category,
        t1.chr,
        t1.region_inf,
        t1.region_sup,
        t2.snp_id,
        t2.snp_pos,
        t2.asm_snp,
        t1.region_nb_cpg,
        t1.nb_cpg_found,
        t1.nb_reads,
        t1.dnase,
        t1.encode_ChiP_V2,
        t1.tf_motifs,
        ROUND(tot_cpg_fm/tot_nb_cpg,3) AS global_cpg_fm,
        t1.tot_nb_cpg,
        t1.tot_nb_reads,
        t1.cpg,
        t1.read
    FROM TOTAL_FRAC_METHYL t1
    LEFT JOIN ${DATASET_PRED}.${SAMPLE}_asm_snp t2
    ON 
        t1.chr = t2.chr AND
        t1.region_inf = t2.region_inf AND
        t1.region_sup = t2.region_sup
    "
