#!/bin/bash

# Here, the instructions to annotate the windows with ASM information
# for the ENCODE samples

# Requires the following files from CloudASM:
# - ${DATASET_CONTEXT}.${SAMPLE}_cpg_asm
# - ${DATASET_CONTEXT}.${SAMPLE}_cpg_read_genotype


#--------------------------------------------------------------------------
# Variables
#--------------------------------------------------------------------------

# Where scripts are located
SCRIPTS="/Users/em/code/DeepASM/asm_annotation"

# BQ dataset where the epigenetic windows are defined
DATASET_EPI="hg19"

# Size of genomic regions:
GENOMIC_INTERVAL="250" # must be the same that in hg19_preparation.sh

# BQ dataset where the output of CloudASM is located
DATASET_PRED="deepasm_feb2022" # T-cells: "tcells_2020" ENCODE: "deepasm_june2020"
 
# BQ dataset where the sample's context files are located (naming defined by CloudASM)
DATASET_CONTEXT="cloudasm_encode_2019" # T-cells: "tcells_2020" ENCODE: "cloudasm_encode_2019" 

# Bucket where to put the txt files for Python analysis
OUTPUT_B="deepasm"

# Cloud Storage location of the logs
LOG="gs://cloudasm-encode/logging/deepasm"

# Docker file required to run the scripts
DOCKER_GCP="google/cloud-sdk:255.0.0"

# Light-weight python Docker image with statistical packages.
DOCKER_PYTHON="gcr.io/hackensack-tyco/python"

# GCP global variables
PROJECT_ID="hackensack-tyco"
REGION_ID="us-central1"
ZONE_ID="us-central1-b"


#--------------------------------------------------------------------------
# ASM Variables
#--------------------------------------------------------------------------

# Minimum number of CpGs we require near a SNP for it to be considered for an ASM region
CPG_PER_ASM_REGION="3"

# Max CpG coverage (to get rid off abherent regions with dozens of thousands of reads overlap.)
MAX_CPG_COV="200"

# p-value cut-off used in all tests for significance
P_VALUE="0.05"

# Benjamin-Hochberg threshold
BH_THRESHOLD="0.05"

# Effect size required at the ASM region level.
ASM_REGION_EFFECT="0.2"

# In an ASM region, minimum bumber of CpGs with significant ASM in the same direction
CPG_SAME_DIRECTION_ASM="3"

# Number of consecutive CpGs with significant ASM in the same direction (among all well-covered CpGs)
CONSECUTIVE_CPG="2" 


#--------------------------------------------------------------------------
# Evaluate ASM in every region overlapping a SNP
#--------------------------------------------------------------------------

# Create a table of regions with CpGs where ASM was calculated (Fisher's test)
# This creates a row for each (snp id, region) combination with at least 3 CpGs 
# All CpGs have 10x coverage and were evaluated for ASM

echo -e "--env SAMPLE\t--env CHR" > all_chr.tsv

# Create a file of job parameters for finding SNPs and their reads.
while read SAMPLE ; do
    for CHR in `seq 1 22` X Y ; do
        echo -e "${SAMPLE}\t${CHR}" >> all_chr.tsv
    done
done < sample_id.txt

# Launch the 288 jobs in parallele (100 at most at the same time)
dsub \
    --project $PROJECT_ID \
    --zones $ZONE_ID \
    --image ${DOCKER_GCP} \
    --logging $LOG \
    --env DATASET_PRED="${DATASET_PRED}" \
    --env DATASET_EPI="${DATASET_EPI}" \
    --env DATASET_CONTEXT="${DATASET_CONTEXT}" \
    --env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
    --env CPG_PER_ASM_REGION="${CPG_PER_ASM_REGION}" \
    --env P_VALUE="${P_VALUE}" \
    --env MAX_CPG_COV="${MAX_CPG_COV}" \
    --script ${SCRIPTS}/cpg_asm.sh \
    --tasks all_chr.tsv \
    --wait

# Delete previous tables
while read SAMPLE ; do
    echo "Deleting the table for sample " ${SAMPLE}
    bq rm -f -t ${DATASET_PRED}.${SAMPLE}_cpg_asm_${GENOMIC_INTERVAL}bp
done < sample_id.txt

# Append to a new table for each sample
{ read
while read SAMPLE CHR ; do 
    echo "Sample is:" ${SAMPLE} ", Chromosome is " ${CHR}
    bq cp --append_table \
        ${DATASET_PRED}.${SAMPLE}_cpg_asm_${CHR} \
        ${DATASET_PRED}.${SAMPLE}_cpg_asm_${GENOMIC_INTERVAL}bp
done 
} < all_chr.tsv

# Check tables
bq query \
    --use_legacy_sql=false \
    "
    SELECT 
        TABLE_NAME, 
        ROUND(TOTAL_ROWS/1000,0) AS thousand_rows 
    FROM ${DATASET_PRED}.INFORMATION_SCHEMA.PARTITIONS
    WHERE TABLE_NAME LIKE '%cpg_asm%_${GENOMIC_INTERVAL}bp'
    "

# Erase intermediary files.
{ read
while read SAMPLE CHR ; do 
    echo "Sample is:" ${SAMPLE} ", Chromosome is " ${CHR}
    bq rm -f -t ${DATASET_PRED}.${SAMPLE}_cpg_asm_${CHR}
done 
} < all_chr.tsv


#### Create a dataset of (region, snp_ip) and arrays of fractional methylation of reads

dsub \
    --project $PROJECT_ID \
    --zones $ZONE_ID \
    --image ${DOCKER_GCP} \
    --logging $LOG \
    --env DATASET_PRED="${DATASET_PRED}" \
    --env DATASET_CONTEXT="${DATASET_CONTEXT}" \
    --env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
    --script ${SCRIPTS}/read_asm.sh \
    --tasks all_samples.tsv \
    --wait

##### Merge the 2 datasets of CpG array and read array

dsub \
    --project $PROJECT_ID \
    --zones $ZONE_ID \
    --image ${DOCKER_GCP} \
    --logging $LOG \
    --env DATASET_PRED="${DATASET_PRED}" \
    --env OUTPUT_B="${OUTPUT_B}" \
    --env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
    --script ${SCRIPTS}/cpg_read_asm.sh \
    --tasks all_samples.tsv \
    --wait


###### Calculate Wilcoxon p-value for regions

# Prepare TSV file
echo -e "--input ASM_REGION\t--output ASM_REGION_PVALUE" > asm_regions.tsv

while read SAMPLE ; do
    echo -e "gs://$OUTPUT_B/${GENOMIC_INTERVAL}bp/$SAMPLE/asm/${SAMPLE}_snp_for_asm_region.json\tgs://$OUTPUT_B/${GENOMIC_INTERVAL}bp/$SAMPLE/asm/${SAMPLE}_asm_region_pvalue.json" >> asm_regions.tsv
done < sample_id.txt


dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --disk-size 100 \
  --machine-type n1-highmem-4 \
  --image ${DOCKER_PYTHON} \
  --logging $LOG \
  --env P_VALUE="${P_VALUE}" \
  --env BH_THRESHOLD="${BH_THRESHOLD}" \
  --script ${SCRIPTS}/asm_region.py \
  --tasks asm_regions.tsv \
  --wait


###### Identify ASM

dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_ID="${DATASET_PRED}" \
  --env OUTPUT_B="${OUTPUT_B}" \
  --env ASM_REGION_EFFECT="${ASM_REGION_EFFECT}" \
  --env CPG_SAME_DIRECTION_ASM="${CPG_SAME_DIRECTION_ASM}" \
  --env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
  --env P_VALUE="${P_VALUE}" \
  --env CONSECUTIVE_CPG="${CONSECUTIVE_CPG}" \
  --script ${SCRIPTS}/summary_asm.sh \
  --tasks all_samples.tsv \
  --wait


#--------------------------------------------------------------------------
# Prepare files for DeepASM
#--------------------------------------------------------------------------


######## Aggregate ASM information with cpg and read FM, and annotation.

######## Add information about the sample (whether it was transformed or not)

# Prepare TSV file for ENCODE samples
echo -e "--env SAMPLE\t--env SAMPLE_STATUS" > sample_status.tsv
echo -e "A549\tmodified" >> sample_status.tsv
echo -e "CD14\tnot_modified" >> sample_status.tsv
echo -e "CD34\tnot_modified" >> sample_status.tsv
echo -e "HeLa_S3\tmodified" >> sample_status.tsv
echo -e "HepG2\tmodified" >> sample_status.tsv
echo -e "fibroblast\tmodified" >> sample_status.tsv
echo -e "mammary_epithelial\tnot_modified" >> sample_status.tsv
echo -e "right_lobe_liver\tnot_modified" >> sample_status.tsv
echo -e "sk_n_sh\tmodified" >> sample_status.tsv
echo -e "spleen_female_adult\tnot_modified" >> sample_status.tsv
echo -e "t_cell_male_adult\tnot_modified" >> sample_status.tsv
echo -e "gm12878\tmodified" >> sample_status.tsv

# Prepare TSV file for t-cells samples
echo -e "--env SAMPLE\t--env SAMPLE_STATUS" > sample_status.tsv
echo -e "tb10206R\tnot_modified" >> sample_status.tsv
echo -e "tb6878R\tnot_modified" >> sample_status.tsv
echo -e "tb6883R\tnot_modified" >> sample_status.tsv

dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_PRED="${DATASET_PRED}" \
  --env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
  --script ${SCRIPTS}/asm_notebook.sh \
  --tasks sample_status.tsv \
  --wait


############## Add genomic window informations

dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_PRED="${DATASET_PRED}" \
  --env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
  --script ${SCRIPTS}/add_genomic_window.sh \
  --tasks all_samples.tsv \
  --wait



##################
######## Concatenate all files per sample

bq rm -f -t ${DATASET_PRED}.all_samples_all_info_${GENOMIC_INTERVAL}bp

while read SAMPLE ; do 
    echo "Sample:" ${SAMPLE}
    bq cp --append_table \
        ${DATASET_PRED}.${SAMPLE}_all_info_${GENOMIC_INTERVAL}bp \
        ${DATASET_PRED}.all_samples_all_info_${GENOMIC_INTERVAL}bp
done < sample_id.txt




#--------------------------------------------------------------------------
# Create table with genomic regions where we know there is ASM
#--------------------------------------------------------------------------

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.all_samples_all_info_${GENOMIC_INTERVAL}bp_for_notebook \
    --replace=true \
    "
    WITH RENAME AS (
        SELECT 
            asm_snp AS asm_snp_tmp, 
            sample_category AS sample_c, 
            * EXCEPT(asm_snp, 
                     sample_category, 
                     wilcoxon_corr_pvalue, 
                     asm_region_effect, 
                     snp_id, 
                     snp_pos)
        FROM ${DATASET_PRED}.all_samples_all_info_${GENOMIC_INTERVAL}bp 
    )
    SELECT 
        IF(asm_snp_tmp = True, 1, IF(asm_snp_tmp = False, 0, -1)) AS asm_snp,
        IF(sample_c = 'not_modified', 0, 1) AS sample_category,
        * EXCEPT(asm_snp_tmp, read, cpg, sample_c, window_details),
        (SELECT ARRAY 
            (SELECT ROUND(fm, 3) FROM UNNEST(read) 
            )
        ) AS read_fm,
        (SELECT ARRAY -- we order by position the FM
            (SELECT ROUND(fm, 3) FROM UNNEST(cpg) ORDER BY pos
            )
        ) AS cpg_fm,
        (SELECT ARRAY -- we order by position. We use Float because python thinks it'a string otherwise
            (SELECT CAST(pos AS FLOAT64) FROM UNNEST(cpg) ORDER BY pos
            ) 
        ) AS cpg_pos,
        window_details AS genomic_picture
    FROM RENAME
    WHERE (asm_snp_tmp = True OR asm_snp_tmp = False) 
    "

#--------------------------------------------------------------------------
# Export ENCODE training set to Cloud Storage
#--------------------------------------------------------------------------

# Delete previous JSON files
gsutil rm gs://${OUTPUT_B}/${GENOMIC_INTERVAL}bp/encode_training_data_with_genomic_picture/encode_training-*.json

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







