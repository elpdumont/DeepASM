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


# Unnest for nb_cpg = 6. (40k regions, 1k ASM)

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_extract_flatten \
    --replace=true \
    "
    WITH SIX_CPG AS (
        SELECT * FROM ${DATASET_ID}.${SAMPLE}_asm_cpg_array
        WHERE nb_cpg = 6
    )
    SELECT asm_snp, snp_id, snp_pos, chr, nb_ref_reads, nb_alt_reads, 
            asm_region_effect, wilcoxon_corr_pvalue, nb_cpg, nb_sig_cpg, cpg_index,
            pos AS cpg_pos, frac_methyl AS cpg_frac_methyl, effect AS cpg_effect, 
            fisher_pvalue AS cpg_fisher_pvalue, ref_cov AS cpg_ref_cov, 
            ref_meth AS cpg_ref_meth, alt_cov AS cpg_alt_cov, alt_meth AS cpg_alt_meth
    FROM SIX_CPG, UNNEST (cpg) WITH OFFSET AS cpg_index
    ORDER BY snp_id
    "

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_for_tf \
    --replace=true \
    "
    WITH LONG_TABLE AS (
        SELECT * FROM ${DATASET_ID}.${SAMPLE}_extract_flatten
    ),
    FIRST_CPG AS (
    SELECT asm_snp, snp_id, snp_pos, chr, nb_ref_reads, nb_alt_reads, 
            asm_region_effect, wilcoxon_corr_pvalue, nb_cpg, nb_sig_cpg,
            cpg_pos AS one_cpg_pos, cpg_frac_methyl AS one_cpg_frac_methyl,
            cpg_effect AS one_cpg_effect, 
            cpg_fisher_pvalue AS one_cpg_fisher_pvalue,
            cpg_ref_cov AS one_cpg_ref_cov, cpg_ref_meth AS one_cpg_ref_meth,
            cpg_alt_cov AS one_cpg_alt_cov, cpg_alt_meth AS one_cpg_alt_meth
    FROM LONG_TABLE
    WHERE cpg_index = 0
    ),
    SECOND_CPG AS (
    SELECT snp_id AS snp_id_two, 
            cpg_pos AS two_cpg_pos, cpg_frac_methyl AS two_cpg_frac_methyl, 
            cpg_effect AS two_cpg_effect, 
            cpg_fisher_pvalue AS two_cpg_fisher_pvalue,
            cpg_ref_cov AS two_cpg_ref_cov, cpg_ref_meth AS two_cpg_ref_meth,
            cpg_alt_cov AS two_cpg_alt_cov, cpg_alt_meth AS two_cpg_alt_meth
    FROM LONG_TABLE
    WHERE cpg_index = 1
    ),
    THIRD_CPG AS (
    SELECT snp_id AS snp_id_three, 
            cpg_pos AS three_cpg_pos, cpg_frac_methyl AS three_cpg_frac_methyl,
            cpg_effect AS three_cpg_effect, 
            cpg_fisher_pvalue AS three_cpg_fisher_pvalue,
            cpg_ref_cov AS three_cpg_ref_cov, cpg_ref_meth AS three_cpg_ref_meth,
            cpg_alt_cov AS three_cpg_alt_cov, cpg_alt_meth AS three_cpg_alt_meth
    FROM LONG_TABLE
    WHERE cpg_index = 2
    ),
    FOURTH_CPG AS (
    SELECT snp_id AS snp_id_four, 
            cpg_pos AS four_cpg_pos, cpg_frac_methyl AS four_cpg_frac_methyl,
            cpg_effect AS four_cpg_effect, 
            cpg_fisher_pvalue AS four_cpg_fisher_pvalue,
            cpg_ref_cov AS four_cpg_ref_cov, cpg_ref_meth AS four_cpg_ref_meth,
            cpg_alt_cov AS four_cpg_alt_cov, cpg_alt_meth AS four_cpg_alt_meth
    FROM LONG_TABLE
    WHERE cpg_index = 3
    ),
    FIFTH_CPG AS (
    SELECT snp_id AS snp_id_five, 
            cpg_pos AS five_cpg_pos, cpg_frac_methyl AS five_cpg_frac_methyl,
            cpg_fisher_pvalue AS five_cpg_fisher_pvalue,
            cpg_ref_cov AS five_cpg_ref_cov, cpg_ref_meth AS five_cpg_ref_meth,
            cpg_alt_cov AS five_cpg_alt_cov, cpg_alt_meth AS five_cpg_alt_meth
    FROM LONG_TABLE
    WHERE cpg_index = 4
    ),
    SIXTH_CPG AS (
    SELECT snp_id AS snp_id_six, 
            cpg_pos AS six_cpg_pos, cpg_frac_methyl AS six_cpg_frac_methyl,
            cpg_effect AS six_cpg_effect, 
            cpg_fisher_pvalue AS six_cpg_fisher_pvalue,
            cpg_ref_cov AS six_cpg_ref_cov, cpg_ref_meth AS six_cpg_ref_meth,
            cpg_alt_cov AS six_cpg_alt_cov, cpg_alt_meth AS six_cpg_alt_meth
    FROM LONG_TABLE
    WHERE cpg_index = 5
    ),
    ONE_TWO AS (
        SELECT * FROM FIRST_CPG INNER JOIN SECOND_CPG ON snp_id = snp_id_two
    ),
    ONE_TWO_THREE AS (
        SELECT * FROM ONE_TWO INNER JOIN THIRD_CPG ON snp_id = snp_id_three
    ),
    ONE_TWO_THREE_FOUR AS (
        SELECT * FROM ONE_TWO_THREE INNER JOIN FOURTH_CPG ON snp_id = snp_id_four
    ),
    ONE_TWO_THREE_FOUR_FIVE AS (
        SELECT * FROM ONE_TWO_THREE_FOUR INNER JOIN FIFTH_CPG ON snp_id = snp_id_five
    ),
    ONE_TWO_THREE_FOUR_FIVE_SIX AS (
        SELECT * FROM ONE_TWO_THREE_FOUR_FIVE INNER JOIN SIXTH_CPG ON snp_id = snp_id_six
    )
    SELECT * EXCEPT(snp_id_two, snp_id_three, snp_id_four, snp_id_five) 
    FROM ONE_TWO_THREE_FOUR_FIVE_SIX
    "