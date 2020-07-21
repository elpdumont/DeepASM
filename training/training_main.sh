#!/bin/bash

# Here, the instructions to annotate the windows with ASM information
# for the ENCODE samples


#--------------------------------------------------------------------------
# Variables
#--------------------------------------------------------------------------

# Where scripts are located
SCRIPTS="/Users/emmanuel/GITHUB_REPOS/DeepASM/training"

# Where CloudASM scripts 
CLOUDASM_SCRIPTS="/Users/emmanuel/GITHUB_REPOS/CloudASM-encode-for-deepasm"

# BQ dataset where the epigenetic windows are defined
DATASET_EPI="hg19"

# Size of genomic regions:
GENOMIC_INTERVAL="250" # must be the same that in hg19_preparation.sh

# BQ dataset where the output of CloudASM is located
DATASET_PRED="deepasm_june2020"

# BQ dataset where the sample's context files are located (naming defined by CloudASM)
DATASET_CONTEXT="cloudasm_encode_2019"

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
    bq rm -f -t ${DATASET_PRED}.${SAMPLE}_cpg_asm
done < sample_id.txt

# Append to a new table for each sample
{ read
while read SAMPLE CHR ; do 
    echo "Sample is:" ${SAMPLE} ", Chromosome is " ${CHR}
    bq cp --append_table \
        ${DATASET_PRED}.${SAMPLE}_cpg_asm_${CHR} \
        ${DATASET_PRED}.${SAMPLE}_cpg_asm
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

dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_PRED="${DATASET_PRED}" \
  --env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
  --script ${SCRIPTS}/asm_annotation.sh \
  --tasks all_samples.tsv \
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
