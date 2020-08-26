
#--------------------------------------------------------------------------
# Variables
#--------------------------------------------------------------------------

# Where scripts are located
SCRIPTS="/Users/emmanuel/GITHUB_REPOS/DeepASM/sample_preparation"

# BQ dataset where the epigenetic windows are defined
DATASET_EPI="hg19"

# Size of genomic regions:
GENOMIC_INTERVAL="250" # must be the same that in hg19_preparation.sh

# Min CpG coverage
MIN_CPG_COV="10"

# Max CpG coverage (to get rid off abherent regions with dozens of thousands of reads overlap.)
MAX_CPG_COV="200"

# Min number of CpGs mapped per genomic region
MIN_NB_CPG="3"

# BQ dataset where the output of CloudASM is located
DATASET_PRED="deepasm_june2020" # For t-cells: "tcells_2020" # For ENCODE: "deepasm_june2020"

# BQ dataset where the sample's context files are located (naming defined by CloudASM)
DATASET_CONTEXT="cloudasm_encode_2019" # For T-cells: "tcells_2020" # For ENCODE: "cloudasm_encode_2019"

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

# Number of nucleotides in each chromosome
echo -e "chr\tlength" > chr_length.txt
echo -e "1\t249250621" > chr_length.txt && echo -e "2\t243199373" >> chr_length.txt && echo -e "3\t198022430" >> chr_length.txt \
&& echo -e "4\t191154276" >> chr_length.txt && echo -e "5\t180915260" >> chr_length.txt && echo -e "6\t171115067" >> chr_length.txt \
&& echo -e "7\t159138663" >> chr_length.txt && echo -e "8\t146364022" >> chr_length.txt && echo -e "9\t141213431" >> chr_length.txt \
&& echo -e "10\t135534747" >> chr_length.txt && echo -e "11\t135006516" >> chr_length.txt && echo -e "12\t133851895" >> chr_length.txt \
&& echo -e "13\t115169878" >> chr_length.txt && echo -e "14\t107349540" >> chr_length.txt && echo -e "15\t102531392" >> chr_length.txt \
&& echo -e "16\t90354753" >> chr_length.txt && echo -e "17\t81195210" >> chr_length.txt && echo -e "18\t78077248" >> chr_length.txt \
&& echo -e "19\t59128983" >> chr_length.txt && echo -e "20\t63025520" >> chr_length.txt && echo -e "21\t48129895" >> chr_length.txt \
&& echo -e "22\t51304566" >> chr_length.txt && echo -e "X\t155270560 " >> chr_length.txt && echo -e "Y\t59373566" >> chr_length.txt

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
--tasks chr_split.tsv 496-577 \
--wait

# For ENCODE
# 1-99.   dstat --provider google-v2 --project hackensack-tyco --jobs 'cpg-region--emmanuel--200825-121649-59' --users 'emmanuel' --status '*'
# 100-198   dstat --provider google-v2 --project hackensack-tyco --jobs 'cpg-region--emmanuel--200825-144156-62' --users 'emmanuel' --status '*'
# 199-297   dstat --provider google-v2 --project hackensack-tyco --jobs 'cpg-region--emmanuel--200825-181153-59' --users 'emmanuel' --status '*'
# 298-396   dstat --provider google-v2 --project hackensack-tyco --jobs 'cpg-region--emmanuel--200825-212440-04' --users 'emmanuel' --status '*'
# 397-495   dstat --provider google-v2 --project hackensack-tyco --jobs 'cpg-region--emmanuel--200825-222608-91' --users 'emmanuel' --status '*'
# 496-577   dstat --provider google-v2 --project hackensack-tyco --jobs 'cpg-region--emmanuel--200826-080156-50' --users 'emmanuel' --status '*'


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
# Keep the "good" CpGs
#--------------------------------------------------------------------------

# We create a table where CpGs satisfying both MAX_CPG_COV and MIN_CPG_COV

dsub \
    --project $PROJECT_ID \
    --zones $ZONE_ID \
    --image ${DOCKER_GCP} \
    --logging $LOG \
    --env DATASET_PRED="${DATASET_PRED}" \
    --env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
    --env MIN_CPG_COV="${MIN_CPG_COV}" \
    --env MAX_CPG_COV="${MAX_CPG_COV}" \
    --script ${SCRIPTS}/clean_cpg.sh \
    --tasks all_samples.tsv \
    --wait


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

# ENCODE: 

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

# Request that the region has at least 3 CpGs.

dsub \
    --project $PROJECT_ID \
    --zones $ZONE_ID \
    --image ${DOCKER_GCP} \
    --logging $LOG \
    --env DATASET_PRED="${DATASET_PRED}" \
    --env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
    --env MIN_NB_CPG="${MIN_NB_CPG}" \
    --script ${SCRIPTS}/combine_read_cpg_fm.sh \
    --tasks all_samples.tsv \
    --wait


#--------------------------------------------------------------------------
# Format the table for prediction by DeepASM (DO NOT RUN IF YOU DO THE ASM ANNOTATION)
#--------------------------------------------------------------------------

# If running the ENCODE samples for training and validating the model,
# move on the to the asm_annotation/asm_main.sh script and do not 
# execute the commands below.


######## Add sample name as a column

dsub \
    --project $PROJECT_ID \
    --zones $ZONE_ID \
    --image ${DOCKER_GCP} \
    --logging $LOG \
    --env DATASET_PRED="${DATASET_PRED}" \
    --env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
    --script ${SCRIPTS}/summary.sh \
    --tasks all_samples.tsv \
    --wait


######## Concatenate all files per sample

bq rm -f -t ${DATASET_PRED}.${DATASET_PRED}_${GENOMIC_INTERVAL}bp

while read SAMPLE ; do 
    echo "Sample:" ${SAMPLE}
    bq cp --append_table \
        ${DATASET_PRED}.${SAMPLE}_cpg_read_${GENOMIC_INTERVAL}bp_name \
        ${DATASET_PRED}.${DATASET_PRED}_${GENOMIC_INTERVAL}bp
done < sample_id.txt

######## Extract the different arrays from structures.

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${DATASET_PRED}_${GENOMIC_INTERVAL}bp_for_notebook \
    --replace=true \
    "
    SELECT 
        -100 AS asm_snp, 
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
    FROM ${DATASET_PRED}.${DATASET_PRED}_${GENOMIC_INTERVAL}bp
    "


