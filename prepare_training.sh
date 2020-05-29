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

# Bucket where CloudASM files are stored
CLOUDASM_BUCKET="cloudasm-encode"

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



# Prepare TSV file per chromosome (used for many jobs)
echo -e "--env CHR" > all_chr.tsv

# Create a file of job parameters for finding SNPs and their reads.
for CHR in `seq 1 22` X Y ; do
      echo -e "${CHR}" >> all_chr.tsv
done




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


#-----------------------------------------------
# Create a table with ASMs and arrays of CpGs.

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


#-----------------------------------------------
# Create a table of SNPs and their arrays of fractional methylation per read

# Import required files created by CloudASM

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


# The table asm_region_pvalue has the fractional methylations.

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

# Combine arrays of fractional methylation arrays to regions evaluated for ASM
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

# Built a table with reads array and cpg arrays for all regions evaluated for ASM
# We also evaluate the width of the region (based on first and last CpG)
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



###################################################
# DNASE track
##################################################


# URL to download from
# https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=830774571_pra4VNR81N6YjQ3NUyzCQSqI7hiT&clade=mammal&org=Human&db=hg19&hgta_group=regulation&hgta_track=wgEncodeRegDnaseClustered&hgta_table=0&hgta_regionType=genome&position=chr21%3A23%2C031%2C598-43%2C031%2C597&hgta_outputType=wigData&hgta_outFileName=dnase.txt

# Do a bash command to remove "chr":
sed -i 's|chr||g' dnase.txt

# Upload to bucket
gsutil cp dnase.txt gs://${CLOUDASM_BUCKET}/dnase.txt

# Push DNASe track to BigQuery
bq --location=US load \
    --replace=true \
    --source_format=CSV \
    --field_delimiter "\t" \
    --skip_leading_rows 1 \
    encode_tracks.dnase \
    gs://${CLOUDASM_BUCKET}/dnase.txt \
    bin:INT64,chr:STRING,chr_start:INT64,chr_end:INT64,name:INT64,score:INT64,source_count:FLOAT,source_id:STRING,source_score:STRING


# Combined DNAse data with ASM
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_OUT="${DATASET_OUT}" \
  --script ${SCRIPTS}/dnase.sh \
  --tasks all_chr.tsv \
  --wait
        

# Concatenate the files
bq rm -f -t ${DATASET_OUT}.asm_read_cpg_dnase

for CHR in `seq 1 22` X Y ; do
    bq cp --append_table \
        ${DATASET_OUT}.asm_read_cpg_dnase_${CHR} \
        ${DATASET_OUT}.asm_read_cpg_dnase
done

for CHR in `seq 1 22` X Y ; do
    bq rm -f -t ${DATASET_OUT}.asm_read_cpg_dnase_${CHR}
done



# Gather all DNASE scores under a structure for a given (sample, snp_id) combination
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_OUT}.asm_read_cpg_dnase_struct \
    --replace=true \
    "
    WITH 
        DNASE_AGG AS (
            SELECT 
                sample AS sample_dnase, 
                snp_id AS snp_id_dnase, 
                ARRAY_AGG(STRUCT(score_dnase)) AS dnase
            FROM ${DATASET_OUT}.asm_read_cpg_dnase
            GROUP BY sample, snp_id
        ),
        OTHER_INFO AS (
            SELECT * 
            FROM ${DATASET_OUT}.asm_read_cpg_arrays
        ),
        COMBINED AS (
            SELECT * FROM OTHER_INFO LEFT JOIN DNASE_AGG
            ON sample_dnase = sample AND snp_id_dnase = snp_id
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
            (SELECT ARRAY (SELECT score_dnase FROM UNNEST(dnase) WHERE score_dnase IS NOT NULL)) AS dnase_scores
        FROM COMBINED
    "


######################################
# ENCODE 3 TFBS

# https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=830774571_pra4VNR81N6YjQ3NUyzCQSqI7hiT&clade=mammal&org=Human&db=hg19&hgta_group=regulation&hgta_track=encTfChipPk&hgta_table=0&hgta_regionType=genome&position=chr21%3A23%2C031%2C598-43%2C031%2C597&hgta_outputType=primaryTable&hgta_outFileName=encode_3_tfbs.txt


# Do a bash command to remove "chr":
sed -i 's|chr||g' encode_3_tfbs.txt

# Upload to bucket
gsutil cp encode_3_tfbs.txt gs://${CLOUDASM_BUCKET}/encode_3_tfbs.txt

bq --location=US load \
    --replace=true \
    --source_format=CSV \
    --field_delimiter "\t" \
    --skip_leading_rows 1 \
    ${DATASET_OUT}.encode_3_tfbs \
    gs://${CLOUDASM_BUCKET}/encode_3_tfbs.txt \
    bin:INT64,chr:STRING,chr_start:INT64,chr_end:INT64,name:STRING,score:INT64,strand:STRING,value:FLOAT,pvalue:FLOAT,qvalue:FLOAT,peak:INT64


# Not enough data to be of interest

###################################################
# TF BINDING FROM CHIP-SEQ DATA
##################################################

# Link of the public dataset
# https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=830774571_pra4VNR81N6YjQ3NUyzCQSqI7hiT&clade=mammal&org=Human&db=hg19&hgta_group=regulation&hgta_track=wgEncodeRegTfbsClusteredV2&hgta_table=0&hgta_regionType=genome&position=chr21%3A23%2C031%2C598-43%2C031%2C597&hgta_outputType=primaryTable&hgta_outFileName=encode_ChiP_V2.txt


# Do a bash command to remove "chr":
sed -i 's|chr||g' encode_ChiP_V2.txt

# Upload to bucket
gsutil cp encode_ChiP_V2.txt gs://${CLOUDASM_BUCKET}/encode_ChiP_V2.txt

# Transfer to BigQuery
bq --location=US load \
    --replace=true \
    --source_format=CSV \
    --field_delimiter "\t" \
    --skip_leading_rows 1 \
    ${DATASET_OUT}.encode_ChiP_V2 \
    gs://${CLOUDASM_BUCKET}/encode_ChiP_V2.txt \
    bin:INT64,chr:STRING,chr_start:INT64,chr_end:INT64,name:STRING,score:INT64,strand:STRING,thick_start:INT64,thick_end:INT64,reserved:INT64,block_count:INT64,block_size:INT64,chrom_start:INT64,exp_count:INT64,exp_id:STRING,exp_score:STRING

# Combined previous ASM data with ChiP-seq data (within +/- 250 bp of first and last CpG)
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_OUT="${DATASET_OUT}" \
  --script ${SCRIPTS}/tf.sh \
  --tasks all_chr.tsv \
  --wait
        

# Concatenate the files
bq rm -f -t ${DATASET_OUT}.asm_read_cpg_tf 

for CHR in `seq 1 22` X Y ; do
    bq cp --append_table \
        ${DATASET_OUT}.asm_read_cpg_tf_${CHR} \
        ${DATASET_OUT}.asm_read_cpg_tf
done

for CHR in `seq 1 22` X Y ; do
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
                ARRAY_AGG(STRUCT(tf_name)) AS tf
            FROM ${DATASET_OUT}.asm_read_cpg_tf
            GROUP BY sample, snp_id
        ),
        OTHER_INFO AS (
            SELECT * 
            FROM ${DATASET_OUT}.asm_read_cpg_arrays
        ),
        COMBINED AS (
            SELECT * FROM OTHER_INFO LEFT JOIN TF_AGG
            ON sample_tf = sample AND snp_id_tf = snp_id
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



###################################################
# TF BINDING MOTIFS
##################################################

# Motifs provided by Catherine Do.

# Clean the database of motifs
mv kherad_tf_sorted.bed kherad_tf_sorted.txt
sed -i 's|chr||g' kherad_tf_sorted.txt

# Upload database to bucket
gsutil cp kherad_tf_sorted.txt gs://${CLOUDASM_BUCKET}/kherad_tf_sorted.txt

# Transfer bucket -> BigQuery
bq --location=US load \
    --replace=true \
    --source_format=CSV \
    --field_delimiter "\t" \
    --skip_leading_rows 0 \
    ${DATASET_OUT}.kherad_tf_sorted \
    gs://${CLOUDASM_BUCKET}/kherad_tf_sorted.txt \
    chr:STRING,motif_start:INT64,motif_end:INT64,motif:STRING

# Motifs known to correlate with ASM (from bioRiv publication)
gsutil cp asm_motifs.txt gs://${CLOUDASM_BUCKET}/asm_motifs.txt

# Upload known ASM motifs to BigQuery
bq --location=US load \
    --replace=true \
    --source_format=CSV \
    --field_delimiter "\t" \
    --skip_leading_rows 0 \
    ${DATASET_OUT}.asm_motifs \
    gs://${CLOUDASM_BUCKET}/asm_motifs.txt \
    asm_motif:STRING

# Keep the motifs known to correlate with ASM
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_OUT}.kherad_tf_sorted_asm_motifs \
    --replace=true \
    "
    WITH 
        ASM_MOTIFS AS (
            SELECT * 
            FROM ${DATASET_OUT}.asm_motifs
        ),
        KHERAD AS (
            SELECT * FROM ${DATASET_OUT}.kherad_tf_sorted
        )
        SELECT chr, motif, motif_start, motif_end FROM KHERAD
        INNER JOIN ASM_MOTIFS
        ON asm_motif = motif
    "


# Combined ASM motifs with ASM hits
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_OUT="${DATASET_OUT}" \
  --script ${SCRIPTS}/motifs.sh \
  --tasks all_chr.tsv \
  --wait
        

# Concatenate the files (one per chromosome)
bq rm -f -t ${DATASET_OUT}.asm_read_cpg_motifs

for CHR in `seq 1 22` X Y ; do
    bq cp --append_table ${DATASET_OUT}.asm_read_cpg_motifs_${CHR} ${DATASET_OUT}.asm_read_cpg_motifs
done

for CHR in `seq 1 22` X Y ; do
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

