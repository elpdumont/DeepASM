#--------------------------------------------------------------------------
# Variables
#--------------------------------------------------------------------------

# Where scripts are located
SCRIPTS="/Users/emmanuel/GITHUB_REPOS/DeepASM/training"

# Where the scripts for enrichment are located
ENRICH_SCRIPTS="/Users/emmanuel/GITHUB_REPOS/DeepASM/enrichment"

# BQ dataset where the output of CloudASM is located
DATASET_IN="cloudasm_encode_2019"

# BQ dataset where the data will be generated
DATASET_OUT="deepasm_encode"

# BQ dataset for epigenetic motifs
DATASET_EPI="hg19"

# Region within which we look for epigenetic signals
EPI_REGION="250"

# Cloud Storage location of the logs
LOG="gs://cloudasm-encode/logging/deepasm"

# Docker file required to run the scripts
DOCKER_GCP="google/cloud-sdk:255.0.0"

# Bucket where CloudASM files are stored
CLOUDASM_BUCKET="cloudasm-encode"

# GCP global variables
PROJECT_ID="hackensack-tyco"
REGION_ID="us-central1"
ZONE_ID="us-central1-b"

#--------------------------------------------------------------------------
# Samples evaluated by CloudASM
#--------------------------------------------------------------------------


# Prepare TSV file with just the samples (used for most jobs)
echo -e "--env SAMPLE" > all_samples.tsv

while read SAMPLE ; do
    echo -e "${SAMPLE}" >> all_samples.tsv
done < sample_id.txt


# Prepare TSV file per chromosome (used for many jobs)
echo -e "--env CHR" > all_chr.tsv

# Create a file of job parameters for finding SNPs and their reads.
for CHR in `seq 1 22` X Y ; do
      echo -e "${CHR}" >> all_chr.tsv
done

#--------------------------------------------------------------------------
# CpG fractional methylation and coverage for all snp_id and all samples
#--------------------------------------------------------------------------

# Delete the existing file in the dataset
bq rm -f -t ${PROJECT_ID}:${DATASET_OUT}.cpgs

# Append all CpG informations to a single table ("cpgs")
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


#--------------------------------------------------------------------------
# Make a table of regions evaluated for ASM for all samples
#--------------------------------------------------------------------------

# Delete existing file
bq rm -f -t ${PROJECT_ID}:${DATASET_OUT}.asm

# Append all samples 
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_IN="${DATASET_IN}" \
  --env DATASET_OUT="${DATASET_OUT}" \
  --script ${SCRIPTS}/asm_region.sh \
  --tasks all_samples.tsv \
  --wait


#--------------------------------------------------------------------------
# Combine CpG arrays and ASM regions together
#--------------------------------------------------------------------------


bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_OUT}.asm_cpg_array \
    --replace=true \
    "
    WITH ASM_SNP AS (
        SELECT 
            IF(asm_snp=True, 1, 0) as asm_snp_value, 
            *
        FROM ${DATASET_OUT}.asm
        ),
    CPG_DETAILS AS (
        SELECT 
            snp_id AS snp_id_tmp,
            sample AS sample_tmp,
            cpg
        FROM ${DATASET_OUT}.cpgs
    ),
    TOGETHER AS (
        SELECT * EXCEPT(asm_snp)
        FROM ASM_SNP 
        INNER JOIN CPG_DETAILS 
        ON snp_id = snp_id_tmp AND sample = sample_tmp
    )
    SELECT 
        asm_snp_value AS asm_snp, 
        * EXCEPT (snp_id_tmp, sample_tmp, asm_snp_value) 
    FROM TOGETHER 
    "


#--------------------------------------------------------------------------
# Create a table of SNPs and their arrays of fractional methylation per genotype
#--------------------------------------------------------------------------

# Import required files created by CloudASM back into the CloudASM dataset.

dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_IN="${DATASET_IN}" \
  --env CLOUDASM_BUCKET="${CLOUDASM_BUCKET}" \
  --command '
        bq --location=US load \
            --autodetect \
            --replace=true \
            --source_format=NEWLINE_DELIMITED_JSON \
            ${DATASET_IN}.${SAMPLE}_asm_region_pvalue \
            gs://${CLOUDASM_BUCKET}/$SAMPLE/asm/${SAMPLE}_asm_region_pvalue.json
        ' \
  --tasks all_samples.tsv \
  --wait


# Concatenate all samples into a single table ("reads")

# Delete existing file
bq rm -f -t ${PROJECT_ID}:${DATASET_OUT}.reads

# Append all samples 
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_IN="${DATASET_IN}" \
  --env DATASET_OUT="${DATASET_OUT}" \
  --script ${SCRIPTS}/reads.sh \
  --tasks all_samples.tsv \
  --wait


#--------------------------------------------------------------------------
# Table with fractional methylation of reads (NOT by genotype)
#--------------------------------------------------------------------------

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_OUT}.asm_read_array \
    --replace=true \
    "
    WITH 
    SNP_READS AS (
        SELECT * 
        FROM ${DATASET_OUT}.reads
    ),
    ASM_RESULTS AS (
        SELECT 
            IF(asm_snp=True, 1, 0) as asm_snp_value, 
            snp_id AS snp_id_results,
            sample AS sample_results
        FROM ${DATASET_OUT}.asm
    ),
    JOINED_ARRAY AS (
        SELECT * 
        FROM SNP_READS INNER JOIN ASM_RESULTS 
        ON snp_id = snp_id_results AND sample = sample_results
    )
    SELECT 
        asm_snp_value AS asm_snp, 
        sample,
        snp_id, 
        snp_pos, 
        chr, 
        ref_reads, 
        alt_reads, 
        ref, 
        alt, 
        (SELECT ARRAY 
            (SELECT methyl_perc FROM UNNEST(REF) 
                UNION ALL SELECT methyl_perc FROM UNNEST(ALT)
            )
        ) AS read_fm
    FROM JOINED_ARRAY
    "


#--------------------------------------------------------------------------
# Table with ASM info, CpG array fractional methyl, and read fractional methyl
#--------------------------------------------------------------------------

# We defined the boundaries of the region as the min and max
# of all CpG positions in the region. 
# These boundaries are used for enrichment with DNASE

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_OUT}.asm_read_cpg_arrays \
    --replace=true \
    "
    WITH
    ASM_READS AS (
        SELECT *
        FROM ${DATASET_OUT}.asm_read_array
    ),
    ASM_CPG AS (
        SELECT 
            snp_id AS snp_id_cpg, 
            sample AS sample_cpg, 
            nb_ref_reads + nb_alt_reads AS nb_reads, 
            nb_cpg, 
            (
                (SELECT max(pos) FROM UNNEST(cpg)) - (SELECT min(pos) FROM UNNEST(cpg))
             ) AS region_length,
             (SELECT min(pos) FROM UNNEST(cpg)) AS region_inf,
             (SELECT max(pos) FROM UNNEST(cpg)) AS region_sup,
            cpg
        FROM ${DATASET_OUT}.asm_cpg_array
    ),
    ASM_READS_CPG_RAW AS (
        SELECT * FROM ASM_READS 
        INNER JOIN ASM_CPG 
        ON snp_id = snp_id_cpg AND sample = sample_cpg
    )
    SELECT 
        asm_snp, 
        sample, 
        snp_id, 
        chr, 
        nb_reads, 
        nb_cpg, 
        region_inf,
        region_sup,
        region_length, 
        read_fm, 
        (SELECT ARRAY 
            (SELECT frac_methyl FROM UNNEST(cpg))
            ) AS cpg_fm,
        (SELECT ARRAY 
            (SELECT cov FROM UNNEST(cpg))
            ) AS cpg_cov,
        (SELECT ARRAY 
            (SELECT pos FROM UNNEST(cpg))
            ) AS cpg_pos
    FROM ASM_READS_CPG_RAW
    "


#--------------------------------------------------------------------------
# DNASE track
#--------------------------------------------------------------------------

# We look for DNAse motifs within 250 bp of the region that 
# was evaluated by CloudASM (boundaries are CpG)

# Combined DNAse data with ASM
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_OUT="${DATASET_OUT}" \
  --env DATASET_EPI="${DATASET_EPI}" \
  --env EPI_REGION="${EPI_REGION}" \
  --script ${SCRIPTS}/dnase.sh \
  --tasks all_chr.tsv \
  --wait


# Concatenate the files
bq rm -f -t ${DATASET_OUT}.asm_read_cpg_dnase

for CHR in `seq 1 22` X Y ; do
    echo "Chromosome is:" ${CHR}
    bq cp --append_table \
        ${DATASET_OUT}.asm_read_cpg_dnase_${CHR} \
        ${DATASET_OUT}.asm_read_cpg_dnase
done

for CHR in `seq 1 22` X Y ; do
    echo "Chromosome is:" ${CHR}
    bq rm -f -t ${DATASET_OUT}.asm_read_cpg_dnase_${CHR}
done


# Gather all DNASE scores under a structure for a given (sample, snp_id) combination
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_OUT}.asm_read_cpg_dnase_struct \
    --replace=true \
    "
    WITH 
        DNASE_AGG AS ( -- we group the DNASe scores together
            SELECT 
                sample AS sample_dnase, 
                snp_id AS snp_id_dnase,
                chr AS chr_dnase,
                region_inf AS region_inf_dnase,
                region_sup AS region_sup_dnase,
                ARRAY_AGG(STRUCT(score_dnase)) AS dnase
            FROM ${DATASET_OUT}.asm_read_cpg_dnase
            GROUP BY 
                sample, 
                snp_id, 
                chr, 
                region_inf, 
                region_sup
        ),
        OTHER_INFO AS (
            SELECT * 
            FROM ${DATASET_OUT}.asm_read_cpg_arrays
        ),
        COMBINED AS (
            SELECT * FROM OTHER_INFO LEFT JOIN DNASE_AGG
            ON 
                sample_dnase = sample AND 
                snp_id_dnase = snp_id AND 
                chr_dnase = chr AND 
                region_inf = region_inf_dnase AND 
                region_sup = region_sup_dnase
        )
        SELECT 
            asm_snp, 
            sample, 
            snp_id, 
            chr, 
            nb_reads, 
            nb_cpg, 
            region_inf, 
            region_sup, 
            region_length, 
            read_fm, 
            cpg_fm, 
            cpg_cov, 
            cpg_pos, 
            -- the command below takes care of the case if there is no dnase score in the array
            (SELECT ARRAY (SELECT score_dnase FROM UNNEST(dnase) WHERE score_dnase IS NOT NULL)) AS dnase_scores
        FROM COMBINED
    "

#--------------------------------------------------------------------------
# TF BINDING FROM CHIP-SEQ DATA
#--------------------------------------------------------------------------

# We look for TF that have measured binding within 250 bp of the CpG window

dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_OUT="${DATASET_OUT}" \
  --env DATASET_EPI="${DATASET_EPI}" \
  --env EPI_REGION="${EPI_REGION}" \
  --script ${SCRIPTS}/tf.sh \
  --tasks all_chr.tsv \
  --wait
        

# Concatenate the files
bq rm -f -t ${DATASET_OUT}.asm_read_cpg_tf 

for CHR in `seq 1 22` X Y ; do
    echo "Chromosome is:" ${CHR}
    bq cp --append_table \
        ${DATASET_OUT}.asm_read_cpg_tf_${CHR} \
        ${DATASET_OUT}.asm_read_cpg_tf
done

for CHR in `seq 1 22` X Y ; do
    echo "Chromosome is:" ${CHR}
    bq rm -f -t ${DATASET_OUT}.asm_read_cpg_tf_${CHR}
done

# Gather all Chip-Seq under a structure for a given (sample, snp_id) combination
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_OUT}.asm_read_cpg_tf_struct \
    --replace=true \
    "
    WITH 
        TF_AGG AS (
            SELECT 
                sample AS sample_tf, 
                snp_id AS snp_id_tf,
                chr AS chr_tf,
                region_inf AS region_inf_tf,
                region_sup AS region_sup_tf,
                ARRAY_AGG(STRUCT(tf_name)) AS tf
            FROM ${DATASET_OUT}.asm_read_cpg_tf
            GROUP BY
                sample, 
                snp_id, 
                chr, 
                region_inf, 
                region_sup
        ),
        OTHER_INFO AS (
            SELECT * 
            FROM ${DATASET_OUT}.asm_read_cpg_arrays
        ),
        COMBINED AS (
            SELECT * FROM OTHER_INFO LEFT JOIN TF_AGG
            ON sample_tf = sample AND 
               snp_id_tf = snp_id AND
               chr_tf = chr AND
               region_inf_tf = region_inf AND
               region_sup_tf = region_sup
        )
        SELECT 
            asm_snp, 
            sample, 
            snp_id, 
            chr, 
            nb_reads, 
            nb_cpg, 
            region_inf, 
            region_sup, 
            region_length, 
            read_fm, 
            cpg_fm, 
            cpg_cov, 
            cpg_pos, 
            (SELECT ARRAY (SELECT tf_name FROM UNNEST(tf) WHERE tf_name IS NOT NULL)) AS tf
        FROM COMBINED
    "

#--------------------------------------------------------------------------
# TF MOTIFS
#--------------------------------------------------------------------------


# Combined ASM motifs with ASM hits
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_OUT="${DATASET_OUT}" \
  --env DATASET_EPI="${DATASET_EPI}" \
  --env EPI_REGION="${EPI_REGION}" \
  --script ${SCRIPTS}/motifs.sh \
  --tasks all_chr.tsv \
  --wait


# Concatenate the files (one per chromosome)
bq rm -f -t ${DATASET_OUT}.asm_read_cpg_motifs

for CHR in `seq 1 22` X Y ; do
    echo "Chromosome is:" ${CHR}
    bq cp --append_table ${DATASET_OUT}.asm_read_cpg_motifs_${CHR} ${DATASET_OUT}.asm_read_cpg_motifs
done

for CHR in `seq 1 22` X Y ; do
    echo "Chromosome is:" ${CHR}
    bq rm -f -t ${DATASET_OUT}.asm_read_cpg_motifs_${CHR}
done


# Gather all motifs under a structure for a given (sample, snp_id) combination
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_OUT}.asm_read_cpg_motifs_struct \
    --replace=true \
    "
    WITH 
        MOTIFS_AGG AS (
            SELECT sample AS sample_motif, snp_id AS snp_id_motif, ARRAY_AGG(STRUCT(motif)) AS tf
            FROM ${DATASET_OUT}.asm_read_cpg_motifs
            GROUP BY sample, snp_id
        ),
        OTHER_INFO AS (
            SELECT * 
            FROM ${DATASET_OUT}.asm_read_cpg_arrays
        ),
        COMBINED AS (
            SELECT * FROM OTHER_INFO LEFT JOIN MOTIFS_AGG
            ON sample_motif = sample AND snp_id_motif = snp_id
        )
        SELECT 
            asm_snp, 
            sample, 
            snp_id, 
            chr, 
            nb_reads, 
            nb_cpg, 
            region_inf, 
            region_sup, 
            region_length, 
            read_fm, 
            cpg_fm, 
            cpg_cov, 
            cpg_pos, 
            (SELECT ARRAY (SELECT motif FROM UNNEST(tf) WHERE motif IS NOT NULL)) AS motifs
        FROM COMBINED
    "



###################################################
# Create a file to export to Notebook
##################################################


# Create a column with a variable to indicate if a motif was found
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_OUT}.asm_for_bq \
    --replace=true \
    "
    WITH 
    MOTIFS AS (
        SELECT * EXCEPT(motifs), 
        (IF(ARRAY_LENGTH(motifs) = 0, 0, 1)) AS motifs_bool 
    FROM ${DATASET_OUT}.asm_read_cpg_motifs_struct
    ),
    CHIPSEQ AS (
        SELECT 
            sample AS sample_chip, 
            snp_id AS snp_id_chip,
            ARRAY_LENGTH(tf) AS nb_tf
        FROM ${DATASET_OUT}.asm_read_cpg_tf_struct
    ),
    DNASE AS (
        SELECT
            sample AS sample_dnase, 
            snp_id AS snp_id_dnase,
            (IF(ARRAY_LENGTH(dnase_scores) = 0, 0, 1)) AS dnase_bool 
        FROM ${DATASET_OUT}.asm_read_cpg_dnase_struct
    ),
    MOTIFS_CHIPSEQ AS (
    SELECT * EXCEPT(sample_chip, snp_id_chip) 
    FROM MOTIFS 
    INNER JOIN CHIPSEQ
    ON sample = sample_chip AND snp_id = snp_id_chip
    )
    SELECT * EXCEPT(sample_dnase, snp_id_dnase) 
    FROM MOTIFS_CHIPSEQ
    INNER JOIN DNASE
    ON sample = sample_dnase AND snp_id = snp_id_dnase
    "

