#--------------------------------------------------------------------------
# Variables
#--------------------------------------------------------------------------

# Where scripts are located
SCRIPTS="/Users/emmanuel/GITHUB_REPOS/DeepASM/training"

# Where the scripts for annotation are located
ANNOTATE_SCRIPTS="/Users/emmanuel/GITHUB_REPOS/DeepASM/annotation"

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
# These boundaries are used for annotation with DNASE

# Created an index for 

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
            (SELECT min(pos) FROM UNNEST(cpg)) AS region_inf,
            -- We add 1 because it's a CpG...
            (1+(SELECT max(pos)FROM UNNEST(cpg))) AS region_sup,
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
        CAST(FLOOR((region_inf + region_sup)/2) AS INT64) AS annotate_ref,
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
# Annotate the CpG windows with additional epigenetic signals.
#--------------------------------------------------------------------------


# Prepare TSV file per chromosome (used for many jobs)
echo -e "--env TABLE\t--env EPI_SIGNAL\t--env CHR" > all_chr.tsv
for SIGNAL in "dnase" "encode_ChiP_V2" "tf_motifs" ; do
    for CHR in `seq 1 22` X Y ; do
    echo -e "asm_read_cpg_arrays\t${SIGNAL}\t${CHR}" >> all_chr.tsv
    done
done


# Combined DNAse data with ASM
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_EPI="${DATASET_EPI}" \
  --env EPI_REGION="${EPI_REGION}" \
  --env DATASET="${DATASET_OUT}" \
  --script ${ANNOTATE_SCRIPTS}/annotation.sh \
  --tasks all_chr.tsv \
  --wait


# Concatenate the files
for EPI_SIGNAL in "dnase" "encode_ChiP_V2" "tf_motifs" ; do
        echo "Processing the signal " ${EPI_SIGNAL}
        bq rm -f -t ${DATASET_OUT}.asm_read_cpg_arrays_${EPI_SIGNAL}_all

        for CHR in `seq 1 22` X Y ; do
            echo "Chromosome is:" ${CHR}
            bq cp --append_table \
                ${DATASET_OUT}.asm_read_cpg_arrays_${EPI_SIGNAL}_${CHR} \
                ${DATASET_OUT}.asm_read_cpg_arrays_${EPI_SIGNAL}_all
        done
    done

# Delete the intermediate files if concatenation was successful
for EPI_SIGNAL in "dnase" "encode_ChiP_V2" "tf_motifs" ; do
    echo "Processing the signal " ${EPI_SIGNAL} "for sample " ${SAMPLE}
    for CHR in `seq 1 22` X Y ; do
        echo "Chromosome is:" ${CHR}
        bq rm -f -t ${DATASET_OUT}.asm_read_cpg_arrays_regions_${EPI_SIGNAL}_${CHR}
    done
done


# Gather all scores in a single array per region. This creates as many tables as there
# are epigenetic signals for annotation.
for EPI_SIGNAL in "dnase" "encode_ChiP_V2" "tf_motifs" ; do
    echo "Processing the signal " ${EPI_SIGNAL} 
    bq query \
        --use_legacy_sql=false \
        --destination_table ${DATASET_OUT}.asm_read_cpg_arrays_${EPI_SIGNAL} \
        --replace=true \
        "
        WITH 
            EPI_AGG AS ( -- we group the DNASe scores together
                SELECT 
                    * EXCEPT(score, read_fm, cpg_fm, cpg_cov, cpg_pos),
                    ARRAY_AGG(STRUCT(score)) AS epi
                FROM ${DATASET_OUT}.asm_read_cpg_arrays_${EPI_SIGNAL}_all
                GROUP BY 
                    chr,
                    sample,
                    snp_id,
                    asm_snp,
                    nb_reads,
                    region_inf,
                    region_sup,
                    annotate_ref,
                    nb_cpg
            )
            SELECT 
                * EXCEPT(epi),
                -- the command below takes care of the case if there is no  score in the array
                IF(
                    ARRAY_LENGTH(
                        (SELECT ARRAY (
                            SELECT score 
                            FROM UNNEST(epi) 
                            WHERE score IS NOT NULL
                            )
                        )) > 0, 1, 0
                ) AS ${EPI_SIGNAL}
            FROM EPI_AGG
        "
done

# Delete the previous tables

for EPI_SIGNAL in "dnase" "encode_ChiP_V2" "tf_motifs" ; do
    echo "Processing the signal " ${EPI_SIGNAL} "for sample " ${SAMPLE}
    bq rm -f -t ${DATASET_OUT}.asm_read_cpg_arrays_${EPI_SIGNAL}_all
done


# Merge the original table (with the CpG and read arrays) with the epigenetic signals.
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_OUT}.asm_annotated \
    --replace=true \
    "
    SELECT
        t1.sample AS sample,
        t1.chr AS chr,
        t1.snp_id AS snp_id,
        t1.asm_snp AS asm_snp,
        t1.nb_reads AS nb_reads,
        t1.region_inf AS region_inf,
        t1.region_sup AS region_sup,
        t1.nb_cpg AS nb_cpg,
        t1.dnase AS dnase,
        t2.encode_ChiP_V2 AS encode_ChiP_V2,
        t3.tf_motifs AS tf_motifs,
        t4.cpg_fm AS cpg_fm,
        t4.read_fm AS read_fm,
        t4.cpg_cov AS cpg_cov,
        t4.cpg_pos AS cpg_pos
    FROM ${DATASET_OUT}.asm_read_cpg_arrays_dnase t1
    JOIN ${DATASET_OUT}.asm_read_cpg_arrays_encode_ChiP_V2 t2 
    ON 
        t1.chr = t2.chr AND 
        t1.sample = t2.sample AND
        t1.snp_id = t2.snp_id AND
        t1.asm_snp = t2.asm_snp AND
        t1.nb_reads = t2.nb_reads AND
        t1.nb_cpg = t2.nb_cpg AND 
        t1.region_inf = t2.region_inf AND 
        t1.region_sup = t2.region_sup
    JOIN ${DATASET_OUT}.asm_read_cpg_arrays_tf_motifs t3 
    ON 
        t1.chr = t3.chr AND 
        t1.sample = t3.sample AND
        t1.snp_id = t3.snp_id AND
        t1.asm_snp = t3.asm_snp AND
        t1.nb_reads = t3.nb_reads AND
        t1.nb_cpg = t3.nb_cpg AND 
        t1.region_inf = t3.region_inf AND 
        t1.region_sup = t3.region_sup
    JOIN ${DATASET_OUT}.asm_read_cpg_arrays t4 
    ON 
        t1.chr = t4.chr AND 
        t1.sample = t4.sample AND
        t1.snp_id = t4.snp_id AND
        t1.asm_snp = t4.asm_snp AND
        t1.nb_reads = t4.nb_reads AND
        t1.nb_cpg = t4.nb_cpg AND 
        t1.region_inf = t4.region_inf AND 
        t1.region_sup = t4.region_sup
    "