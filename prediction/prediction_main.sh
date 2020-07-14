
#--------------------------------------------------------------------------
# Variables
#--------------------------------------------------------------------------

# Where scripts are located
SCRIPTS="/Users/emmanuel/GITHUB_REPOS/DeepASM/prediction"

# BQ dataset where the epigenetic windows are defined
DATASET_EPI="hg19"

# Size of genomic regions:
GENOMIC_INTERVAL="250" # must be the same that in hg19_preparation.sh

# BQ dataset where the output of CloudASM is located
DATASET_PRED="deepasm_june2020"

# BQ dataset where the sample's context files are located (naming defined by CloudASM)
DATASET_CONTEXT="cloudasm_encode_2019"

# Cloud Storage location of the logs
LOG="gs://cloudasm-encode/logging/deepasm"

# Docker file required to run the scripts
DOCKER_GCP="google/cloud-sdk:255.0.0"

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



#--------------------------------------------------------------------------
# Samples to be prepared for prediction
#--------------------------------------------------------------------------

# Prepare TSV file with just the samples (used for most jobs)
echo -e "--env SAMPLE" > all_samples.tsv

while read SAMPLE ; do
    echo -e "${SAMPLE}" >> all_samples.tsv
done < sample_id.txt


#--------------------------------------------------------------------------
# Create CpG regions to be evaluated by DeepASM
#--------------------------------------------------------------------------

# Create genomic regions used to split jobs per chromosome 
INTERVAL="82000000"

# Prepare TSV file per chromosome (used for many jobs)
echo -e "SAMPLE\tCHR\tLOWER_B\tUPPER_B" > chr_split.tsv

# Create the windows with $INTERVALS bp in it
while read SAMPLE ; do 
    for CHR in `seq 1 22` X Y ; do
        echo "Processing chromosome" ${CHR} " for sample" ${SAMPLE}
        NUCLEOTIDES_IN_CHR=$(awk -v CHR="${CHR}" -F"\t" '{ if ($1 == CHR) print $2}' chr_length.txt)
        INF="1"
        SUP=$(( $NUCLEOTIDES_IN_CHR<$INTERVAL ? $NUCLEOTIDES_IN_CHR: $INTERVAL ))
        echo -e "${SAMPLE}\t${CHR}\t$INF\t$SUP" >> chr_split.tsv # for jobs
        while [ $NUCLEOTIDES_IN_CHR -gt $SUP ] ; do
            INCREMENT=$(( $NUCLEOTIDES_IN_CHR-$SUP<$INTERVAL ? $NUCLEOTIDES_IN_CHR-$SUP: $INTERVAL ))
            INF=$(( ${SUP} + 1 ))
            SUP=$(( ${SUP} + $INCREMENT ))
            echo -e "${SAMPLE}\t${CHR}\t$INF\t$SUP" >> chr_split.tsv

        done
    done
done < sample_id.txt

# Looping over the samples. This step requires manual intervention.
dsub \
--project $PROJECT_ID \
--zones $ZONE_ID \
--image ${DOCKER_GCP} \
--logging $LOG \
--env DATASET_PRED="${DATASET_PRED}" \
--env DATASET_EPI="${DATASET_EPI}" \
--env DATASET_CONTEXT="${DATASET_CONTEXT}" \
--env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
--script ${SCRIPTS}/cpg_regions.sh \
--tasks chr_split.tsv \
--wait

# 1-99, 100-198, 199-297, 298-396, 397-496, 497-577

# Delete previous tables
while read SAMPLE ; do
    echo "Deleting the table for sample " ${SAMPLE}
    bq rm -f -t ${DATASET_PRED}.${SAMPLE}_cpg_regions_${GENOMIC_INTERVAL}bp
done < sample_id.txt

# Append to a new table for each sample
{ read
while IFS=$'\t' read SAMPLE CHR LOWER_B UPPER_B ; do 
    echo "Sample is:" ${SAMPLE} ", Chromosome is " ${CHR} ", lower:" ${LOWER_B} ", and upper:" ${UPPER_B}
    bq cp --append_table \
        ${DATASET_PRED}.${SAMPLE}_cpg_regions_${CHR}_${LOWER_B}_${UPPER_B} \
        ${DATASET_PRED}.${SAMPLE}_cpg_regions_${GENOMIC_INTERVAL}bp
done 
} < chr_split.tsv



# Erase intermediary files.
{ read
while read SAMPLE CHR LOWER_B UPPER_B ; do 
    echo "Sample is:" ${SAMPLE} ", Chromosome is " ${CHR} ", lower:" ${LOWER_B} ", and upper:" ${UPPER_B}
    bq rm -f -t ${DATASET_PRED}.${SAMPLE}_cpg_regions_${CHR}_${LOWER_B}_${UPPER_B}
done 
} < chr_split.tsv



#--------------------------------------------------------------------------
# Calculate fractional methylation of each CpG in each region.
#--------------------------------------------------------------------------

# We request that each CpG is covered at least 10x
# We request that there are at least 3 CpGs that are covered at least 10x in each region

dsub \
    --project $PROJECT_ID \
    --zones $ZONE_ID \
    --image ${DOCKER_GCP} \
    --logging $LOG \
    --env DATASET_PRED="${DATASET_PRED}" \
    --env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
    --script ${SCRIPTS}/cpg_fm.sh \
    --tasks all_samples.tsv \
    --wait

#--------------------------------------------------------------------------
# Calculate fractional methylation of each read in each region.
#--------------------------------------------------------------------------

dsub \
    --project $PROJECT_ID \
    --zones $ZONE_ID \
    --image ${DOCKER_GCP} \
    --logging $LOG \
    --env DATASET_PRED="${DATASET_PRED}" \
    --env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
    --script ${SCRIPTS}/read_fm.sh \
    --tasks all_samples.tsv \
    --wait

#--------------------------------------------------------------------------
# Create a table of regions with fractional methylation of CpGs and reads
#--------------------------------------------------------------------------

dsub \
    --project $PROJECT_ID \
    --zones $ZONE_ID \
    --image ${DOCKER_GCP} \
    --logging $LOG \
    --env DATASET_PRED="${DATASET_PRED}" \
    --env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
    --script ${SCRIPTS}/combine_read_cpg_fm.sh \
    --tasks all_samples.tsv \
    --wait


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
    --script ${SCRIPTS}/cpg_regions_asm.sh \
    --tasks all_chr.tsv 10 \
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



for CHR in `seq 1 22` X Y ; do
    bq rm -f -t ${DATASET_PRED}.${SAMPLE}_cpg_asm_${CHR}
done
