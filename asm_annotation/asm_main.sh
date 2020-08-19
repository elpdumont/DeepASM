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
SCRIPTS="/Users/emmanuel/GITHUB_REPOS/DeepASM/asm_annotation"

# Where CloudASM scripts 
CLOUDASM_SCRIPTS="/Users/emmanuel/GITHUB_REPOS/CloudASM-encode-for-deepasm"

# BQ dataset where the epigenetic windows are defined
DATASET_EPI="hg19"

# Size of genomic regions:
GENOMIC_INTERVAL="250" # must be the same that in hg19_preparation.sh

# BQ dataset where the output of CloudASM is located
DATASET_PRED="deepasm_june2020" # tcells_2020"

# BQ dataset where the sample's context files are located (naming defined by CloudASM)
DATASET_CONTEXT="cloudasm_encode_2019" # "tcells_2020" # 

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

# Launch the jobs in parallele (100 at most at the same time)
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
    --script ${SCRIPTS}/cpg_asm.sh \
    --tasks all_chr.tsv \
    --wait

# 1-99, 100-199, 200 - 289

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
    echo -e "gs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_snp_for_asm_region.json\tgs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_asm_region_pvalue.json" >> asm_regions.tsv
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
  --script ${CLOUDASM_SCRIPTS}/asm_region.py \
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

# Prepare TSV file
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


######## Concatenate all files per sample

bq rm -f -t ${DATASET_PRED}.all_samples_${GENOMIC_INTERVAL}bp

while read SAMPLE ; do 
    echo "Sample:" ${SAMPLE}
    bq cp --append_table \
        ${DATASET_PRED}.${SAMPLE}_cpg_read_asm_${GENOMIC_INTERVAL}bp \
        ${DATASET_PRED}.all_samples_${GENOMIC_INTERVAL}bp
done < sample_id.txt


######## Format for using in the DeepASM notebook

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.all_samples_${GENOMIC_INTERVAL}bp_for_notebook \
    --replace=true \
    "
    WITH RENAME AS (
        SELECT asm_snp AS asm_snp_tmp, * EXCEPT(asm_snp)
        FROM ${DATASET_PRED}.all_samples_${GENOMIC_INTERVAL}bp 
    )
    SELECT 
        IF(asm_snp_tmp = True, 1, IF(asm_snp_tmp = False, 0, -1)) AS asm_snp,
        * EXCEPT(asm_snp_tmp, read, cpg),
        (SELECT ARRAY 
            (SELECT fm FROM UNNEST(read) 
            )
        ) AS read_fm,
        (SELECT ARRAY 
            (SELECT fm FROM UNNEST(cpg) 
            )
        ) AS cpg_fm,
        (SELECT ARRAY 
            (SELECT cov FROM UNNEST(cpg) 
            )
        ) AS cpg_cov,
        (SELECT ARRAY 
            (SELECT pos FROM UNNEST(cpg) ORDER BY pos
            ) 
        ) AS cpg_pos 
    FROM RENAME
    "


#--------------------------------------------------------------------------
# Key numbers about the number of regions evaluated
#--------------------------------------------------------------------------

# Number of distinct regions in the ref genome (with 3 CpGs)
# 3,790,920

bq query --use_legacy_sql=false \
    "
    WITH DISTINCT_REGIONS AS (
        SELECT DISTINCT chr, region_inf, region_sup
        FROM ${DATASET_EPI}.hg19_cpg_regions_${GENOMIC_INTERVAL}bp_annotated
    )
    SELECT COUNT(*) FROM DISTINCT_REGIONS
    "

# Number of distinct regions evaluated by CLOUDASM across all ENCODE samples
# 1,419,549 (37% of all regions)

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
# CloudASM evaluated about 10% of all regions with potential ASM
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
               gs://deepasm/csv_encode_with_asm_prob.csv 


bq rm -f -t tcells_2020.all_tcells_with_pred

while read SAMPLE ; do 
    echo "Processing sample" ${SAMPLE}
    bq --location=US load \
                --replace=false \
                --source_format=CSV \
                --skip_leading_rows 1 \
                tcells_2020.all_tcells_with_pred \
                gs://deepasm/csv_${SAMPLE}_with_asm_prob.csv \
                asm_probability:FLOAT64,asm_snp:INTEGER,sample:STRING,chr:STRING,region_inf:INT64,region_sup:INT64,snp_id:STRING,snp_pos:FLOAT64,region_nb_cpg:INT64,nb_cpg_found:INT64,dnase:INT64,encode_ChiP_V2:INT64,tf_motifs:INT64,read_fm:STRING,cpg_fm:STRING,cpg_cov:STRING,cpg_pos:STRING


done < sample_id.txt




#--------------------------------------------------------------------------
# Comparison with ASM found by DeepASM and CloudASM in the ENCODE t-cell
#--------------------------------------------------------------------------

# Universe: 365,034
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM deepasm_june2020.all_encode_with_pred
    WHERE sample = 't_cell_male_adult'
    "

# False negative: 48,728
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM deepasm_june2020.all_encode_with_pred
    WHERE asm_probability > 0.5 AND asm_snp = 0 AND sample = 't_cell_male_adult'
    "

# True positive: 146
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM deepasm_june2020.all_encode_with_pred
    WHERE asm_probability > 0.5 AND asm_snp = 1 AND sample = 't_cell_male_adult'
    "

# True negative: 315,129
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM deepasm_june2020.all_encode_with_pred
    WHERE asm_probability < 0.5 AND asm_snp = 0 AND sample = 't_cell_male_adult'
    "

# False negative: 1031
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM deepasm_june2020.all_encode_with_pred
    WHERE asm_probability < 0.5 AND asm_snp = 1 AND sample = 't_cell_male_adult'
    "

# Accuracy: (80+202,316)/212,186 = 95%
# Precision: 80/(80+7407) = 1%
# Sensibility: 
# VPN = TN/(T)
# Recall:
# AUC: 

#--------------------------------------------------------------------------
# Comparison with ASM found by DeepASM and CloudASM in the 3 T-cells
#--------------------------------------------------------------------------

# Universe: 212,186
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM ${DATASET_PRED}.all_tcells_with_pred
    WHERE asm_snp > -1
    "

# False negative: 7,407
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM ${DATASET_PRED}.all_tcells_with_pred
    WHERE asm_probability > 0.5 AND asm_snp = 0
    "

# True positive: 80
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM ${DATASET_PRED}.all_tcells_with_pred
    WHERE asm_probability > 0.5 AND asm_snp = 1
    "

# True negative: 202,316
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM ${DATASET_PRED}.all_tcells_with_pred
    WHERE asm_probability < 0.5 AND asm_snp = 0
    "

# False negative: 2,383
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM ${DATASET_PRED}.all_tcells_with_pred
    WHERE asm_probability < 0.5 AND asm_snp = 1
    "

# Accuracy: (80+202,316)/212,186 = 95%
# Precision: 80/(80+7407) = 1%
# Sensibility: 
# VPN = TN/(T)
# Recall:
# AUC: 


#--------------------------------------------------------------------------
# Comparison with ASM found by mQTL (only for T-cells)
#--------------------------------------------------------------------------

# Table of ASM found by mQTL in T-cells:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4863666/bin/mmc2.xlsx


bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.tcells_mQTL_DeepASM \
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
        FROM ${DATASET_PRED}.all_tcells_with_pred t1
        RIGHT JOIN ${DATASET_PRED}.tcells_mQTL t2
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

# Total number of CpGs evaluated by mQTL and found to have ASM at the single CpG level
# 1,440
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM ${DATASET_PRED}.tcells_mQTL_DeepASM
    "

# Number of regions evaluated by DeepASM 
# 1,269
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM ${DATASET_PRED}.tcells_mQTL_DeepASM
    WHERE max_asm_probability IS NOT NULL
    "

# Number of regions evaluated by CloudASM within DeepASM regions
# 526 
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM ${DATASET_PRED}.tcells_mQTL_DeepASM
    WHERE cloudasm_asm > -1 AND max_asm_probability IS NOT NULL
    "

# Number of regions evaluated by DeepASM only
# 743
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM ${DATASET_PRED}.tcells_mQTL_DeepASM
    WHERE max_asm_probability IS NOT NULL AND cloudasm_asm = -1
    "


# Within these 743, regions with ASM found by DeepASM
# 64 (8.6%)
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM ${DATASET_PRED}.tcells_mQTL_DeepASM
    WHERE cloudasm_asm = -1 AND max_asm_probability > 0.5
    "

# Number of regions with ASM found by DeepASM within regions evaluated by CloudASM
# 64/526 (12%)
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM ${DATASET_PRED}.tcells_mQTL_DeepASM
    WHERE cloudasm_asm > -1 AND max_asm_probability > 0.5
    "

# Number of ASM = 1 found by CloudASM 
# 169 (32% of )
bq query --use_legacy_sql=false \
    "
    SELECT COUNT(*)
    FROM ${DATASET_PRED}.tcells_mQTL_DeepASM
    WHERE cloudasm_asm = 1 
    "






