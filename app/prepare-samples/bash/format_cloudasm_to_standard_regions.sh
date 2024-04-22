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

# Construct the query to count rows in the table
nb_regions=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${SAMPLES_DATASET}.regions_w_arrays")
echo "Across ${NB_SAMPLES} samples, we have ${nb_regions} regions"

"${script_folder}"/bash/form_regions_with_snps.sh

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

echo "Evaluated ${nb_regions} regions, of which only ${nb_regions_w_asm} have ASM"

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
            CASE
                WHEN AVG(asm) > 0 THEN 1 --we require finding at least one asm
                ELSE 0
            END AS asm
        FROM ${SAMPLES_DATASET}.asm_flagged
        GROUP BY sample, chr, clustering_index, region_inf, region_sup
    )
    SELECT c.*, p.asm FROM ${SAMPLES_DATASET}.regions_w_arrays c
    LEFT JOIN ASM p
    ON c.sample = p.sample AND c.chr = p.chr AND c.clustering_index = p.clustering_index AND c.region_inf = p.region_inf AND c.region_sup = p.region_sup
    "

nb_regions=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${SAMPLES_DATASET}.unique_regions_w_asm_flagged")
nb_regions_evaluated=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${SAMPLES_DATASET}.unique_regions_w_asm_flagged WHERE asm IS NOT NULL")
nb_regions_w_asm=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${SAMPLES_DATASET}.unique_regions_w_asm_flagged WHERE asm = 1")
echo "We have ${nb_regions} regions. We evaluated ${nb_regions_evaluated} regions for ASM. Among these, there are ${nb_regions_w_asm} regions where we found ASM"

echo "Exporting the data to the bucket"
bq extract --destination_format=NEWLINE_DELIMITED_JSON "${PROJECT}:${SAMPLES_DATASET}.unique_regions_w_asm_flagged" gs://"${BUCKET}"/"${DATA_PATH}"/after_cloudasm/all_regions/unique_regions_w_asm_flagged_*.json

nb_files=$(gsutil ls gs://${BUCKET}/${DATA_PATH}/all_regions/* | wc -l)
echo "There are ${nb_files} that cover all the regions with an ASM flag when ASM could be evaluated"


bq query \
    --use_legacy_sql=false \
    --destination_table "${SAMPLES_DATASET}".unique_regions_w_asm_0_or_1 \
    --replace=true \
    "
    SELECT * FROM ${SAMPLES_DATASET}.unique_regions_w_asm_flagged
    WHERE asm IS NOT NULL
    "

echo "Exporting the data to the bucket with ASM being 0 or 1"
bq extract --destination_format=NEWLINE_DELIMITED_JSON "${PROJECT}:${SAMPLES_DATASET}.unique_regions_w_asm_0_or_1" gs://"${BUCKET}"/"${DATA_PATH}"/after_cloudasm/regions_with_known_asm/unique_regions_w_asm_0_or_1_*.json

#---------------------------------------------
echo "Preparing the features for all the regions (even if they do not have ASM flagged)"
NB_FILES_PER_TASK="50"
NB_FILES_TO_PROCESS=$(gsutil ls gs://"${BUCKET}"/"${DATA_PATH}"/all_regions/* | wc -l | awk '{print $1}')
TASK_COUNT=$(( (${NB_FILES_TO_PROCESS} + ${NB_FILES_PER_TASK} - 1) / ${NB_FILES_PER_TASK} ))
sed -i '' "s#NB_FILES_PER_TASK_PH#${NB_FILES_PER_TASK}#g" "batch-jobs/prepare_features_for_ML.json"
sed -i '' "s#TASK_COUNT_PH#${TASK_COUNT}#g" "batch-jobs/prepare_features_for_ML.json"


gcloud batch jobs submit "prepare-features-for-ml-${SHORT_SHA}" \
	--location "${REGION}" \
	--config batch-jobs/prepare_features_for_ML.json

nb_regions=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${ML_DATASET}.features_wo_hmm")
nb_regions_w_data=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${ML_DATASET}.features_wo_hmm WHERE cpgs_w_padding IS NOT NULL ")
nb_regions_w_data_and_asm_flagged=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${ML_DATASET}.features_wo_hmm WHERE cpgs_w_padding IS NOT NULL AND asm IS NOT NULL")
echo "Features were prepared for ${nb_regions} regions. Among these, we could extract features for ${nb_regions_w_data} regions. Among these, we have ${nb_regions_w_data_and_asm_flagged} regions with features and ASM."

echo "Exporting the dataset with features (excluding HMM) to the bucket"
bq extract --destination_format=NEWLINE_DELIMITED_JSON "${PROJECT}:${ML_DATASET}.features_wo_hmm" gs://"${BUCKET}"/"${DATA_PATH}"/features_wo_hmm/features_wo_hmm_*.json



#---------------------------------------------
echo "Fitting an HMM model on the training set and infering the states-based features for all datasets"











#--------------------------------------------------------------------------
# Export ENCODE training set to Cloud Storage
#--------------------------------------------------------------------------

# Delete previous JSON files
gsutil -m rm gs://${OUTPUT_B}/${GENOMIC_INTERVAL}bp/encode_training_data_with_genomic_picture/encode_training-*.json

# Export the table into several JSON files (~15 for the 12 ENCODE samples)
bq extract \
    --location=US \
    --destination_format=NEWLINE_DELIMITED_JSON \
    ${DATASET_PRED}.all_samples_all_info_${GENOMIC_INTERVAL}bp_for_notebook \
    gs://${OUTPUT_B}/${GENOMIC_INTERVAL}bp/encode_training_data_with_genomic_picture/encode_training-*.json




#--------------------------------------------------------------------------
# Key numbers about the number of regions evaluated
#--------------------------------------------------------------------------

# Number of distinct regions in the ref genome (with 3 CpGs)
# 3,790,920 (250 bp)
# 2,587,039 (1000 bp)

bq query --use_legacy_sql=false \
    "
    WITH DISTINCT_REGIONS AS (
        SELECT DISTINCT chr, region_inf, region_sup
        FROM ${DATASET_EPI}.hg19_cpg_regions_${GENOMIC_INTERVAL}bp_clean_annotated
    )
    SELECT COUNT(*) FROM DISTINCT_REGIONS
    "

# Number of distinct regions evaluated by CLOUDASM across all ENCODE samples
# 1,419,549 (37% of all regions, 250 bp)
# 1,192,312 (1000 bp)

bq query --use_legacy_sql=false \
    "
    WITH DISTINCT_REGIONS AS (
        SELECT DISTINCT chr, region_inf, region_sup
        FROM ${DATASET_PRED}.all_samples_${GENOMIC_INTERVAL}bp
        WHERE asm_snp IS NOT NULL
    )
    SELECT COUNT(*) FROM DISTINCT_REGIONS
    "

# Number of distinct regions evaluated by CLOUDASM and DEEPASM by sample
# CloudASM evaluated about 10% of all regions with potential ASM, 10-20% when using 1000 bp
bq query --use_legacy_sql=false \
    "
    WITH DISTINCT_REGIONS AS (
        SELECT sample, chr, region_inf, region_sup, asm_snp
        FROM ${DATASET_PRED}.all_samples_${GENOMIC_INTERVAL}bp
    ),
    NO_SNP_DS AS (
        SELECT sample, asm_snp, COUNT(*) AS no_snp 
        FROM DISTINCT_REGIONS
        WHERE asm_snp IS NULL
        GROUP BY sample, asm_snp
    ),
    ASM_DS AS (
        SELECT sample, asm_snp, COUNT(*) AS asm
        FROM DISTINCT_REGIONS
        WHERE asm_snp = true
        GROUP BY sample, asm_snp
    ),
    NO_ASM_DS AS (
        SELECT sample, asm_snp, COUNT(*) AS no_asm
        FROM DISTINCT_REGIONS
        WHERE asm_snp = false
        GROUP BY sample, asm_snp
    ),
    ALL_JOIN AS (
        SELECT
            t1.sample,
            t1.no_snp,
            t2.asm,
            t3.no_asm
        FROM NO_SNP_DS t1
        INNER JOIN ASM_DS t2
        ON t1.sample = t2.sample
        INNER JOIN NO_ASM_DS t3
        ON t1.sample = t3.sample
    )
    SELECT 
        *, 
        ROUND(100*asm/(asm+no_asm),3) AS asm_perc, 
        ROUND(100*(asm+no_asm)/(asm+no_asm+no_snp)) AS cloudasm_cov 
    FROM ALL_JOIN
    "




#--------------------------------------------------------------------------
# Export ASM predictions to BigQuery
#--------------------------------------------------------------------------

# ENCODE table.
bq --location=US load \
               --replace=true \
               --autodetect \
               --source_format=CSV \
               --skip_leading_rows 1 \
               deepasm_june2020.all_encode_with_pred \
               gs://deepasm/encode_asm_prob.csv 




while read SAMPLE ; do 
    echo "Processing sample" ${SAMPLE}
    bq --location=US load \
                --replace=false \
                --autodetect \
                --source_format=CSV \
                --skip_leading_rows 1 \
                tcells_2020.all_tcells_with_pred_encode_model \
                gs://deepasm/${SAMPLE}_asm_prob_encode_model.csv 
done < sample_id.txt

bq --location=US load \
    --replace=false \
    --autodetect \
    --source_format=CSV \
    --skip_leading_rows 1 \
    tcells_2020.tb10206R_chr15_cloudasm_regions_encode_model \
    gs://deepasm/tb10206R_chr15_cloudasm_regions_encode_model.csv 

bq --location=US load \
    --replace=false \
    --autodetect \
    --source_format=CSV \
    --skip_leading_rows 1 \
    tcells_2020.tb10206R_allchr_cloudasm_regions_encode_model \
    gs://deepasm/tb10206R_allchr_cloudasm_regions_encode_model.csv 

bq --location=US load \
    --replace=false \
    --autodetect \
    --source_format=CSV \
    --skip_leading_rows 1 \
    tcells_2020.tb10206R_allchr_cloudasm_regions_encode_model_2 \
    gs://deepasm/tb10206R_allchr_cloudasm_regions_encode_model_2.csv 


8/27

bq --location=US load \
    --replace=false \
    --autodetect \
    --source_format=CSV \
    --skip_leading_rows 1 \
    tcells_2020.tcells_2020_08_26_encode_v3 \
    gs://deepasm/tcells_2020_08_26_encode_v3.csv 


bq --location=US load \
    --replace=false \
    --autodetect \
    --source_format=CSV \
    --skip_leading_rows 1 \
    deepasm_june2020.encode_2020_08_26_encode_v3 \
    gs://deepasm/encode_2020_08_26_encode_v3.csv 

bq --location=US load \
    --replace=true \
    --autodetect \
    --source_format=CSV \
    --skip_leading_rows 1 \
    deepasm_june2020.encode_tcell_2020_08_26_encode_v3 \
    gs://deepasm/encode_tcell_2020_08_26_encode_v3.csv 


bq --location=US load \
    --replace=true \
    --autodetect \
    --source_format=CSV \
    --skip_leading_rows 1 \
    deepasm_june2020.encode_tcell_all_regions_2020_08_26_encode_v3 \
    gs://deepasm/encode_tcell_all_regions_2020_08_26_encode_v3.csv 



#---------- DONE

DATASET_ID="tcells_2020" # "deepasm_june2020"

bq rm -f -t ${DATASET_ID}.all_regions_2020_08_26_encode_v3

bq --location=US load \
    --replace=false \
    --autodetect \
    --source_format=CSV \
    --skip_leading_rows 1 \
    ${DATASET_ID}.all_regions_2020_08_26_encode_v3 \
    gs://deepasm/tb10206R_all_regions_2020_08_26_encode_v3.csv 

bq --location=US load \
    --replace=false \
    --autodetect \
    --source_format=CSV \
    --skip_leading_rows 1 \
    ${DATASET_ID}.all_regions_2020_08_26_encode_v3 \
    gs://deepasm/tb6878R_all_regions_2020_08_26_encode_v3.csv 

bq --location=US load \
    --replace=false \
    --autodetect \
    --source_format=CSV \
    --skip_leading_rows 1 \
    ${DATASET_ID}.all_regions_2020_08_26_encode_v3 \
    gs://deepasm/tb6883R_all_regions_2020_08_26_encode_v3.csv 


#--------------------------------------------------------------------------
# Model evaluation
#--------------------------------------------------------------------------

DATASET_ID="tcells_2020" # "deepasm_june2020"
TABLE="all_regions_2020_08_26_encode_v3" # "deepasm_june2020.all_encode_with_pred"

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${TABLE}_model_summary \
    --replace=true \
    "
    WITH 
        FALSE_POSITIVE AS (
            SELECT 
                sample, 
                COUNT(*) AS fp
            FROM ${DATASET_ID}.${TABLE}
            WHERE asm_probability > 0.5 AND asm_snp = 0
            GROUP BY sample
        ),
        FALSE_NEGATIVE AS (
            SELECT 
                sample, 
                COUNT(*) AS fn
            FROM ${DATASET_ID}.${TABLE}
            WHERE asm_probability < 0.5 AND asm_snp = 1
            GROUP BY sample
        ),
        TRUE_NEGATIVE AS (
            SELECT 
                sample, 
                COUNT(*) AS tn
            FROM ${DATASET_ID}.${TABLE}
            WHERE asm_probability < 0.5 AND asm_snp = 0
            GROUP BY sample
        ),
        TRUE_POSITIVE AS (
            SELECT 
                sample, 
                COUNT(*) AS tp
            FROM ${DATASET_ID}.${TABLE}
            WHERE asm_probability > 0.5 AND asm_snp = 1
            GROUP BY sample
        ),
        ASM_PER AS (
            SELECT 
                sample, 
                COUNT(*) AS total_asm
            FROM ${DATASET_ID}.${TABLE}
            WHERE asm_snp = 1
            GROUP BY sample
        ),
        TOTAL_REGIONS AS (
            SELECT
                sample,
                global_cpg_fm,
                sample_category,
                COUNT(*) AS total_regions
            FROM ${DATASET_ID}.${TABLE}
            WHERE asm_snp > -1
            GROUP BY sample, global_cpg_fm, sample_category
        )
        SELECT 
            t1.sample, 
            t1.sample_category,
            t1.global_cpg_fm,
            ROUND(100*t6.total_asm/t1.total_regions,3) AS asm_perc,
            t2.fp,
            t3.fn,
            t4.tn,
            t5.tp,
            t6.total_asm,
            t1.total_regions,
            ROUND((t5.tp/(t5.tp + t3.fn)),3) AS sensitivity,
            ROUND((t4.tn/(t2.fp + t4.tn)),3) AS specificity,
            ROUND(((t4.tn+t5.tp)/(t2.fp + t3.fn + t4.tn + t5.tp)),3) AS accuracy
        FROM TOTAL_REGIONS t1 
        INNER JOIN FALSE_POSITIVE t2 ON t1.sample = t2.sample
        INNER JOIN FALSE_NEGATIVE t3 ON t1.sample = t3.sample
        INNER JOIN TRUE_NEGATIVE t4 ON t1.sample = t4.sample
        INNER JOIN TRUE_POSITIVE t5 ON t1.sample = t5.sample
        INNER JOIN ASM_PER t6 ON t1.sample = t6.sample
    "



#--------------------------------------------------------------------------
# Comparison with ASM found by mQTL (only for T-cells)
#--------------------------------------------------------------------------

# Table of ASM found by mQTL in T-cells:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4863666/bin/mmc2.xlsx

# After uploading the file to a bucket:
bq --location=US load \
    --replace=true \
    --autodetect \
    --source_format=CSV \
    --skip_leading_rows 1 \
    tcells_2020.tcells_mQTL \
    gs://deepasm/mQTL_tcells.csv


DATASET_ID="tcells_2020" # "deepasm_june2020"
TABLE="all_regions_2020_08_26_encode_v3" # "deepasm_june2020.all_encode_with_pred"

# Table with all the mQTL regions and their intersect with DeepASM
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.tcells_mQTL_DeepASM \
    --replace=true \
    "
    WITH LONG_TABLE AS(
        SELECT
            t2.rank,
            t2.r_square,
            t2.probe_coor,
            t2.strongest_snp,
            t2.dist_cpg_snp,
            t2.chr,
            t1.asm_probability,
            t1.asm_snp,
            t1.sample,
            t1.region_inf,
            t1.region_sup        
        FROM ${DATASET_ID}.${TABLE} t1
        RIGHT JOIN ${DATASET_ID}.tcells_mQTL t2
        ON 
            t1.region_inf <= t2.probe_coor AND 
            t1.region_sup >= t2. probe_coor AND
            t1.chr = t2.chr
    )
    SELECT
        rank,
        r_square,
        probe_coor,
        strongest_snp,
        dist_cpg_snp,
        ANY_VALUE(region_inf) AS region_inf,
        ANY_VALUE(region_sup) AS region_sup,
        MAX(asm_probability) AS max_asm_probability,
        MAX(asm_snp) AS cloudasm_asm,
    FROM LONG_TABLE
    GROUP BY rank, probe_coor, strongest_snp, dist_cpg_snp, r_square
    "


# Metrics to compare mQTL, DeepASM, and CloudASM


# Total number of CpGs evaluated by mQTL and found to have ASM at the single CpG level
# 1,440
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM ${DATASET_ID}.tcells_mQTL_DeepASM
    "

#-------- CLOUDASM
# Number of regions evaluated by CloudASM 
# 526 (36% of all possible regions with ASM)
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM ${DATASET_PRED}.tcells_mQTL_DeepASM
    WHERE cloudasm_asm > -1 
    "

# Number of ASM found by CloudASM 
# 169 (32% of all regions evaluated by CloudASM)
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM ${DATASET_PRED}.tcells_mQTL_DeepASM
    WHERE cloudasm_asm = 1 
    "

#-------- DEEPASM
# Number of regions evaluated by DeepASM 
# 1,268 (88%)
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM ${DATASET_PRED}.tcells_mQTL_DeepASM
    WHERE max_asm_probability IS NOT NULL
    "

# Number of regions evaluated by DeepASM and NOT CloudASM
# 742 (51%)
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM ${DATASET_PRED}.tcells_mQTL_DeepASM
    WHERE max_asm_probability IS NOT NULL AND cloudasm_asm = -1
    "

# Within these 742, regions with ASM found by DeepASM
# 301 (40%)
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM ${DATASET_PRED}.tcells_mQTL_DeepASM
    WHERE cloudasm_asm = -1 AND max_asm_probability > 0.5
    "







