
#--------------------------------------------------------------------------
# Variables
#--------------------------------------------------------------------------

# Where scripts are located
SCRIPTS="/Users/emmanuel/GITHUB_REPOS/DeepASM/prediction"

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
# Format the table for prediction by DeepASM
#--------------------------------------------------------------------------

# If running the ENCODE samples for training and validating the model,
# move on the to the training/training_man.sh script and do not 
# execute the commands below.


######## Concatenate all files per sample

bq rm -f -t ${DATASET_PRED}.tcells_${GENOMIC_INTERVAL}bp

while read SAMPLE ; do 
    echo "Sample:" ${SAMPLE}
    bq cp --append_table \
        ${DATASET_PRED}.${SAMPLE}_cpg_read_${GENOMIC_INTERVAL}bp \
        ${DATASET_PRED}.tcells_${GENOMIC_INTERVAL}bp
done < sample_id.txt

######## Extract the different arrays from structures.

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.tcells_${GENOMIC_INTERVAL}bp_for_notebook \
    --replace=true \
    "
    SELECT 
        '${SAMPLE}' AS sample,
        * EXCEPT(read, cpg),
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
    FROM ${DATASET_PRED}.all_samples_${GENOMIC_INTERVAL}bp 
    "


