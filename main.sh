#-----------------------------------------------
# Variables

# Where scripts are located
SCRIPTS="/Users/emmanuel/GITHUB_REPOS/DeepASM"

# BQ dataset where the output of CloudASM is located
DATASET_IN="cloudasm_encode_2019"

# BQ dataset where the data will be generated
DATASET_OUT="deepasm_encode"

# Cloud Storage location of the logs
LOG="gs://cloudasm-encode/logging/deepasm"

# Docker file required to run the scripts
DOCKER_GCP="google/cloud-sdk:255.0.0"

# GCP global variables
PROJECT_ID="hackensack-tyco"
REGION_ID="us-central1"
ZONE_ID="us-central1-b"


#-----------------------------------------------

# Prepare TSV file with just the samples (used for most jobs)
echo -e "--env SAMPLE" > all_samples.tsv

while read SAMPLE ; do
    echo -e "${SAMPLE}" >> all_samples.tsv
done < sample_id.txt


#-----------------------------------------------

# Make a file of CpG coverage and 
# fractional methylation for all samples

# Delete the existing file in the dataset
bq rm -f -t ${PROJECT_ID}:${DATASET_OUT}.cpgs

dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_IN="${DATASET_IN}" \
  --env DATASET_OUT="${DATASET_OUT}" \
  --script ${SCRIPTS}/cpg.sh \
  --tasks all_samples.tsv \
  --wait

#-----------------------------------------------

# Make a file of regions evaluated for ASM for all samples

# Delete existing file
bq rm -f -t ${PROJECT_ID}:${DATASET_OUT}.asm

# Append all samples through a while loop so that we do not exceed the rate
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_IN="${DATASET_IN}" \
  --env DATASET_OUT="${DATASET_OUT}" \
  --script ${SCRIPTS}/asm_region.sh \
  --tasks all_samples.tsv 12 \
  --wait


#-----------------------------------------------

# Combine with asm_snp table.

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_OUT}.asm_cpg_array \
    --replace=true \
    "
    WITH ASM_SNP AS (
        SELECT IF(asm_snp=True, 1, 0) as asm_snp_value, *
        FROM ${DATASET_OUT}.asm
        ),
    CPG_DETAILS AS (
        SELECT snp_id AS snp_id_tmp, cpg
        FROM ${DATASET_OUT}.cpgs
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