#!/bin/bash

# Here, the instructions to annotate the windows with ASM information
# for the ENCODE samples

# Requires the following files from CloudASM:

# SO WE PROBABLY DO NOT NEED CONTEXT_FILTERED

# - context_filtered: all cpgs across all reads (~1B rows per sample with CpG pos, methylation status, and read ID)
# - cpg_asm: cpgs near a SNP with the corresponding SNP, allele methylation, coverage, and fisher p-value. (~2M CpGs per sample)
# - cpg_read_genotype: all cpgs across all reads with corresponding allele (REF OR ALT). ~100M CpGs per sample (50x coverage on average)

script_folder="app/prepare-samples"


#-----------------------------------------------
echo "Assemble the 3 datasets (context_filtered, cpg_asm, cpg_read_genotype) for each sample"

# Append to a new table for each sample
# Loop over all samples
job_name="format-cloudasm-${SHORT_SHA}"
gcloud batch jobs submit "${job_name}" \
	--location "${REGION}" \
	--config batch-jobs/overlap_cloudasm_data_with_standard_regions.json


# Export these 3 datasets to the bucket where long-term storage is less expensive than on BQ
echo "Exporting the 3 tables used for ASM calculation and for forming all regions"
for TABLE in "${CLOUDASM_TABLES[@]}"; do
    echo "Exporting table ${TABLE}"
    bq extract --destination_format=NEWLINE_DELIMITED_JSON "${PROJECT}:${SAMPLES_DATASET}.${TABLE}" gs://"${BUCKET}"/"${DATA_PATH}"/before_cloudasm/"${TABLE}"/"${TABLE}"_*.json

    # USE THE CODE BELOW TO IMPORT THESE TABLES FROM THE BUCKET
    # bq query \
    # --use_legacy_sql=false \
    # --destination_table "${SAMPLES_DATASET}".${TABLE} \
    # --replace=true \
    # --clustering_fields=sample,chr \
    # --range_partitioning=clustering_index,0,4000,1 \
    # "
    # SELECT * FROM ${SAMPLES_DATASET}.${TABLE}_tmp
    # "

done

echo "Keep the CpGs with a min and max coverage (defined in config file) for quality insurance"

bq query \
    --use_legacy_sql=false \
    --destination_table "${SAMPLES_DATASET}".all_cpgs_flagged_w_regions \
    --replace=true \
    --clustering_fields=sample,chr \
    --range_partitioning=clustering_index,0,4000,1 \
    "
    WITH
        GOOD_CPG AS (
        -- Creates a list of all CpG sites with their coverage
            SELECT 
                chr, 
                pos, 
                SUM(cov) AS sum_cov,
                clustering_index,
                sample
            FROM ${SAMPLES_DATASET}.context_filtered
            GROUP BY chr, pos, sample, clustering_index
            HAVING sum_cov >= ${MIN_CPG_COV} AND sum_cov <= ${MAX_CPG_COV}
        )
        -- recreate a long table with CpG and read information
        SELECT p.*
        FROM ${SAMPLES_DATASET}.context_filtered p
        INNER JOIN GOOD_CPG c
        ON p.chr = c.chr AND p.pos = c.pos AND c.clustering_index = p.clustering_index AND c.sample = p.sample
    "

#-----------------------------------------------
echo "For each region, group cpg, reads in a nested structure"

# NEED TO BETTER CALCULATE THE FM OF READS TO AVOID COUNTING THE SAME CPG TWICE
# IF IT IS LINKED TO 2 DIFFERENT SNPS

# FORM "DUPLICATE" REGIONS DEPENDING ON WHICH SNP WE LOOK AT
# ONCE ASM IS EVALUATED, COMBINE THE TWO TABLES.

"${script_folder}"/bash/form_regions_with_arrays.sh

bq extract --destination_format=NEWLINE_DELIMITED_JSON "${PROJECT}:${SAMPLES_DATASET}.regions_w_arrays" gs://"${BUCKET}"/"${DATA_PATH}"/before_cloudasm/regions_w_arrays/regions_w_arrays_*.json


# Construct the query to count rows in the table
nb_regions=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${SAMPLES_DATASET}.regions_w_arrays")
echo "Across ${NB_SAMPLES} samples, we have ${nb_regions} regions"

"${script_folder}"/bash/form_regions_with_snps.sh

bq extract --destination_format=NEWLINE_DELIMITED_JSON "${PROJECT}:${SAMPLES_DATASET}.regions_and_snps" gs://"${BUCKET}"/"${DATA_PATH}"/before_cloudasm/regions_and_snps/regions_and_snps_*.json


nb_regions=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${SAMPLES_DATASET}.regions_and_snps")
nb_regions_unique=$(execute_query "SELECT COUNT(*) AS unique_row_count FROM ( SELECT sample, chr, region_inf, region_sup, clustering_index FROM ${PROJECT}.${SAMPLES_DATASET}.regions_and_snps GROUP BY sample, chr, region_inf, region_sup, clustering_index)")
echo "Across ${NB_SAMPLES} samples, we have ${nb_regions} regions overlapping a SNP, of which ${nb_regions_unique} regions are unique"


#-----------------------------------------------
echo "Evaluating ASM in regions that have a SNP. Correcting Wilcoxon p-value across each sample"

gcloud batch jobs submit "evaluate-asm-${SHORT_SHA}" \
	--location "${REGION}" \
	--config batch-jobs/evaluate_asm_in_regions_w_snp.json


nb_regions=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${SAMPLES_DATASET}.asm_flagged")
nb_regions_w_asm=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${SAMPLES_DATASET}.asm_flagged WHERE asm = 1")
nb_regions_w_asm_not_corrected=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${SAMPLES_DATASET}.asm_flagged WHERE asm_not_corrected = 1")

echo "Evaluated ${nb_regions} regions, of which only ${nb_regions_w_asm} have ASM (Wilcoxon Corrected) and ${nb_regions_w_asm_not_corrected} have ASM (Wilcoxon NOT corrected)"

# echo "Eliminate regions where ASM is ambiguous (2 SNPs, different result), partitioning the data, and adding the results to the main table."


# bq query \
#     --use_legacy_sql=false \
#     --destination_table "${SAMPLES_DATASET}".unique_regions_w_asm_flagged \
#     --replace=true \
#     --clustering_fields=sample,chr \
#     --range_partitioning=clustering_index,0,4000,1 \
#     "
#     WITH ASM AS (
#         SELECT
#             sample,
#             chr,
#             clustering_index,
#             region_inf,
#             region_sup,
#             CASE
#                 WHEN AVG(asm) > 0 THEN 1 --we require finding at least one asm
#                 ELSE 0
#             END AS asm
#         FROM ${SAMPLES_DATASET}.asm_flagged
#         GROUP BY sample, chr, clustering_index, region_inf, region_sup
#     )
#     SELECT c.*, p.asm FROM ${SAMPLES_DATASET}.regions_w_arrays c
#     LEFT JOIN ASM p
#     ON c.sample = p.sample AND c.chr = p.chr AND c.clustering_index = p.clustering_index AND c.region_inf = p.region_inf AND c.region_sup = p.region_sup
#     "

echo "Eliminate regions where ASM is ambiguous (2 SNPs, different result), partitioning the data, and adding the results to the main table."
bq query \
    --use_legacy_sql=false \
    --destination_table "${SAMPLES_DATASET}".unique_regions_w_asm_flagged \
    --replace=true \
    --clustering_fields=sample,chr \
    --range_partitioning=clustering_index,0,4000,1 \
    "
    WITH ASM AS (
        SELECT
            sample,
            chr,
            clustering_index,
            region_inf,
            region_sup,
            ANY_VALUE(wilcoxon_pvalue) AS wilcoxon_pvalue,
            ANY_VALUE(corrected_wilcoxon_pvalue) AS corrected_wilcoxon_pvalue,
            ANY_VALUE(total_sig_cpgs) AS total_sig_cpgs,
            ANY_VALUE(consecutive_sig_cpgs) AS consecutive_sig_cpgs,
            ANY_VALUE(read_asm_effect) AS read_asm_effect,
            ANY_VALUE(asm) AS asm,
            ANY_VALUE(asm_not_corrected) AS asm_not_corrected
        FROM ${SAMPLES_DATASET}.asm_flagged
        GROUP BY sample, chr, clustering_index, region_inf, region_sup
        HAVING COUNT(*) = 1
    )
    SELECT
        c.*,
        p.asm,
        p.asm_not_corrected,
        p.wilcoxon_pvalue,
        p.corrected_wilcoxon_pvalue,
        p.total_sig_cpgs,
        p.consecutive_sig_cpgs,
        p.read_asm_effect
    FROM ${SAMPLES_DATASET}.regions_w_arrays c
    LEFT JOIN ASM p
    ON c.sample = p.sample AND c.chr = p.chr AND c.clustering_index = p.clustering_index AND c.region_inf = p.region_inf AND c.region_sup = p.region_sup
    "

nb_regions=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${SAMPLES_DATASET}.unique_regions_w_asm_flagged")
nb_regions_evaluated=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${SAMPLES_DATASET}.unique_regions_w_asm_flagged WHERE asm IS NOT NULL")
nb_regions_w_asm=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${SAMPLES_DATASET}.unique_regions_w_asm_flagged WHERE asm = 1")
nb_regions_w_asm_not_corrected=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${SAMPLES_DATASET}.unique_regions_w_asm_flagged WHERE asm_not_corrected = 1")

echo "We have ${nb_regions} regions. We evaluated ${nb_regions_evaluated} regions for ASM. Among these, there are ${nb_regions_w_asm} regions where we found ASM. If we do not correct the wilcoxon p-value, then we find ${nb_regions_w_asm_not_corrected} regions with ASM."

echo "Exporting the data to the bucket"
bq extract --destination_format=NEWLINE_DELIMITED_JSON "${PROJECT}:${SAMPLES_DATASET}.unique_regions_w_asm_flagged" gs://"${BUCKET}"/"${DATA_PATH}"/after_cloudasm/all_regions/unique_regions_w_asm_flagged_*.json

nb_files=$(gsutil ls gs://${BUCKET}/${DATA_PATH}/after_cloudasm/all_regions/* | wc -l)
echo "There are ${nb_files} that cover all the regions with an ASM flag when ASM could be evaluated"


# bq query \
#     --use_legacy_sql=false \
#     --destination_table "${SAMPLES_DATASET}".unique_regions_w_asm_0_or_1 \
#     --replace=true \
#     "
#     SELECT * FROM ${SAMPLES_DATASET}.unique_regions_w_asm_flagged
#     WHERE asm IS NOT NULL
#     "

# echo "Exporting the data to the bucket with ASM being 0 or 1"
# bq extract --destination_format=NEWLINE_DELIMITED_JSON "${PROJECT}:${SAMPLES_DATASET}.unique_regions_w_asm_0_or_1" gs://"${BUCKET}"/"${DATA_PATH}"/after_cloudasm/regions_with_known_asm/unique_regions_w_asm_0_or_1_*.json



#--------------------------------------------------------------------------
# Comparison with ASM found by mQTL (only for T-cells)
#--------------------------------------------------------------------------

# Table of ASM found by mQTL in T-cells:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4863666/bin/mmc2.xlsx

# # After uploading the file to a bucket:
# bq --location=US load \
#     --replace=true \
#     --autodetect \
#     --source_format=CSV \
#     --skip_leading_rows 1 \
#     tcells_2020.tcells_mQTL \
#     gs://deepasm/mQTL_tcells.csv


# DATASET_ID="tcells_2020" # "deepasm_june2020"
# TABLE="all_regions_2020_08_26_encode_v3" # "deepasm_june2020.all_encode_with_pred"

# # Table with all the mQTL regions and their intersect with DeepASM
# bq query \
#     --use_legacy_sql=false \
#     --destination_table ${DATASET_ID}.tcells_mQTL_DeepASM \
#     --replace=true \
#     "
#     WITH LONG_TABLE AS(
#         SELECT
#             t2.rank,
#             t2.r_square,
#             t2.probe_coor,
#             t2.strongest_snp,
#             t2.dist_cpg_snp,
#             t2.chr,
#             t1.asm_probability,
#             t1.asm_snp,
#             t1.sample,
#             t1.region_inf,
#             t1.region_sup        
#         FROM ${DATASET_ID}.${TABLE} t1
#         RIGHT JOIN ${DATASET_ID}.tcells_mQTL t2
#         ON 
#             t1.region_inf <= t2.probe_coor AND 
#             t1.region_sup >= t2. probe_coor AND
#             t1.chr = t2.chr
#     )
#     SELECT
#         rank,
#         r_square,
#         probe_coor,
#         strongest_snp,
#         dist_cpg_snp,
#         ANY_VALUE(region_inf) AS region_inf,
#         ANY_VALUE(region_sup) AS region_sup,
#         MAX(asm_probability) AS max_asm_probability,
#         MAX(asm_snp) AS cloudasm_asm,
#     FROM LONG_TABLE
#     GROUP BY rank, probe_coor, strongest_snp, dist_cpg_snp, r_square
#     "


# # Metrics to compare mQTL, DeepASM, and CloudASM


# # Total number of CpGs evaluated by mQTL and found to have ASM at the single CpG level
# # 1,440
# bq query --use_legacy_sql=false \
#     "
#     SELECT COUNT(*)
#     FROM ${DATASET_ID}.tcells_mQTL_DeepASM
#     "

# #-------- CLOUDASM
# # Number of regions evaluated by CloudASM 
# # 526 (36% of all possible regions with ASM)
# bq query --use_legacy_sql=false \
#     "
#     SELECT COUNT(*)
#     FROM ${DATASET_PRED}.tcells_mQTL_DeepASM
#     WHERE cloudasm_asm > -1 
#     "

# # Number of ASM found by CloudASM 
# # 169 (32% of all regions evaluated by CloudASM)
# bq query --use_legacy_sql=false \
#     "
#     SELECT COUNT(*)
#     FROM ${DATASET_PRED}.tcells_mQTL_DeepASM
#     WHERE cloudasm_asm = 1 
#     "

# #-------- DEEPASM
# # Number of regions evaluated by DeepASM 
# # 1,268 (88%)
# bq query --use_legacy_sql=false \
#     "
#     SELECT COUNT(*)
#     FROM ${DATASET_PRED}.tcells_mQTL_DeepASM
#     WHERE max_asm_probability IS NOT NULL
#     "

# # Number of regions evaluated by DeepASM and NOT CloudASM
# # 742 (51%)
# bq query --use_legacy_sql=false \
#     "
#     SELECT COUNT(*)
#     FROM ${DATASET_PRED}.tcells_mQTL_DeepASM
#     WHERE max_asm_probability IS NOT NULL AND cloudasm_asm = -1
#     "

# # Within these 742, regions with ASM found by DeepASM
# # 301 (40%)
# bq query --use_legacy_sql=false \
#     "
#     SELECT COUNT(*)
#     FROM ${DATASET_PRED}.tcells_mQTL_DeepASM
#     WHERE cloudasm_asm = -1 AND max_asm_probability > 0.5
#     "







