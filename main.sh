# Import from Big Query

DATASET_ID="DeepASM_gm12878"
SAMPLE="gm12878"

# Starting point: SAMPLE_context_filtered and SAMPLE_vcf_read

# Note: later we will not need VCF. We use the VCF to group reads by SNP for now.


# Extract the number of methylation states per ref and alt for each CpG.

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_snp_cpg_details \
    --replace=true \
    "
    WITH 
    -- SNPs with their respective arrays of CpGs (each coming with their fisher p-value)
        SNP_CPG_ARRAY AS (
            SELECT
                snp_id,
                ANY_VALUE(snp_pos) AS snp_pos,
                ANY_VALUE(chr) AS chr,
                ARRAY_AGG(
                    STRUCT(
                        pos,
                        ROUND((alt_meth+ref_meth)/(ref_cov+alt_cov),3) AS frac_methyl,
                        ref_cov + alt_cov AS cov,
                        ref_meth + alt_meth as meth,
                        ROUND(alt_meth/alt_cov-ref_meth/ref_cov,3) AS effect,
                        fisher_pvalue,
                        ref_cov,
                        ref_meth,
                        alt_cov,
                        alt_meth
                    )
                    ORDER BY pos
                ) AS cpg
            FROM ${DATASET_ID}.${SAMPLE}_cpg_asm
            GROUP BY snp_id
        )
        SELECT * FROM SNP_CPG_ARRAY
    "


# Combine with asm_snp table.

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_asm_cpg_array \
    --replace=true \
    "
    WITH ASM_SNP AS (
        SELECT IF(asm_snp=True, 1, 0) as asm_snp_value, *
        FROM ${DATASET_ID}.${SAMPLE}_asm_snp
        ),
    CPG_DETAILS AS (
        SELECT snp_id AS snp_id_tmp, cpg
        FROM ${DATASET_ID}.${SAMPLE}_snp_cpg_details
    ),
    TOGETHER AS (
        SELECT * EXCEPT(asm_snp)
        FROM ASM_SNP 
        INNER JOIN CPG_DETAILS 
        ON snp_id = snp_id_tmp
    )
    SELECT asm_snp_value AS asm_snp, * EXCEPT (snp_id_tmp, asm_snp_value) FROM TOGETHER 
    "



#-------------------------------
# USE THE METHYLATION DATA ACROSS READS.
#-------------------------------



# Make a table with all snp_ids and their reads (fractional methylation and number of reads)

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_snp_reads \
    --replace=true \
    'SELECT snp_id, snp_pos, chr, alt_reads, ref_reads, ref, alt
    FROM hackensack-tyco.cloudasm_gm12878.gm12878_asm_region_pvalue
    '

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_asm_reads \
    --replace=true \
    "
    WITH 
    SNP_READS AS (
        SELECT * 
        FROM ${DATASET_ID}.${SAMPLE}_snp_reads
    ),
    ASM_RESULTS AS (
        SELECT IF(asm_snp=True, 1, 0) as asm_snp_value, snp_id AS snp_id_results
        FROM ${DATASET_ID}.${SAMPLE}_asm_snp
    ),
    JOINED_ARRAY AS (
        SELECT * 
        FROM SNP_READS INNER JOIN ASM_RESULTS 
        ON snp_id = snp_id_results
    )
    SELECT asm_snp_value AS asm_snp, 
    snp_id, snp_pos, chr, ref_reads, alt_reads, ref, alt, 
    (SELECT ARRAY 
        (SELECT methyl_perc FROM UNNEST(REF) 
                UNION ALL SELECT methyl_perc FROM UNNEST(ALT)
            )
        ) AS read_fm
    FROM JOINED_ARRAY
    "

# Built a table with reads array and cpg arrays

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_asm_reads_cpg_arrays \
    --replace=true \
    "
    WITH
    ASM_READS AS (
        SELECT *
        FROM ${DATASET_ID}.${SAMPLE}_asm_reads
    ),
    ASM_CPG AS (
        SELECT snp_id AS snp_id_cpg, nb_ref_reads + nb_alt_reads AS nb_reads, nb_cpg, 
            (
                (SELECT max(pos) FROM UNNEST(cpg)) - (SELECT min(pos) FROM UNNEST(cpg))
             ) AS region_length,
            cpg
        FROM ${DATASET_ID}.${SAMPLE}_asm_cpg_array
    ),
    ASM_READS_CPG_RAW AS (
        SELECT * FROM ASM_READS 
        INNER JOIN ASM_CPG 
        ON snp_id = snp_id_cpg
    )
    SELECT asm_snp, '${SAMPLE}' AS sample, chr, nb_reads, nb_cpg, region_length, read_fm, 
           (SELECT ARRAY 
                (SELECT frac_methyl FROM UNNEST(cpg))) AS cpg_fm,  
           snp_id
    FROM ASM_READS_CPG_RAW
    "