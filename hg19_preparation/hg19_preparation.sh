#!/bin/bash

# Here we prepare the ref genome with where at least 3 CpGs are present.
# Then we annotate the windows with at least 3 CpGs present.

# Requirements
# CpG bed file in bucket/ref_genomes/reference_genome


export LC_NUMERIC="en_US.UTF-8"
function format_number_with_comma() {
    local number_string=$1
    # Remove existing non-numeric characters (e.g., commas) to ensure the input is a valid integer
    local clean_number=$(echo "$number_string" | sed 's/[^0-9]//g')
    printf "%'d\n" "${clean_number}"
}

source src/env_variables.sh
export FOLDER="${REFERENCE_GENOME}_${GENOMIC_LENGTH}bp_refgenome"
echo "Folder name: ${FOLDER}"

echo "Creating the dataset with an expiration"
if bq ls --project_id="${PROJECT_ID}" | grep -w "${FOLDER}"; then
	echo "Dataset ${FOLDER} already exists in project ${PROJECT_ID}."
else
	# Create the dataset since it does not exist
	bq mk --project_id="${PROJECT_ID}" --dataset --default_table_expiration="${REF_GENOME_DATASET_EXPIRATION_SEC}" "${PROJECT_ID}:${FOLDER}"
	echo "Dataset ${FOLDER} created in project ${PROJECT_ID}."
fi


chmod +x src/make_chr_info_for_ref_genome.sh
src/make_chr_info_for_ref_genome.sh

gsutil cp chr_length.txt gs://"${BUCKET_NAME}"/"${FOLDER}"/chr_length.txt

bq --location=US load \
               --replace=true \
               --source_format=CSV \
               --field_delimiter "\t" \
               --skip_leading_rows 1 \
                "${FOLDER}".chr_length \
               gs://"${BUCKET_NAME}"/"${FOLDER}"/chr_length.txt \
               chr:INT64,chr_length:INT64,absolute_nucleotide_pos:INT64

#--------------------------------------------------------------------------
# Upload a file of the CpG positions in the reference genome
#--------------------------------------------------------------------------

# Transfer the CpG positions to the dataset
bq --location=US load \
               --replace=true \
               --source_format=CSV \
               --field_delimiter "\t" \
               --skip_leading_rows 1 \
                "${FOLDER}".CpG_pos_imported \
               gs://"${BUCKET_NAME}"/ref_genomes/"${REFERENCE_GENOME}"/"${REFERENCE_GENOME}"_CpG_pos.bed \
               chr:STRING,region_inf:INT64,region_sup:INT64

# Converting chr to INT for easy partitioning
bq query \
    --use_legacy_sql=false \
    "
    CREATE OR REPLACE TABLE ${FOLDER}.CpG_pos
    CLUSTER BY clustering_index AS
    SELECT
        CAST(p.chr AS INT64) AS chr,
        p.region_inf,
        p.region_sup,
        CAST(FLOOR((c.absolute_nucleotide_pos + p.region_inf) / ${NB_NUCLEOTIDES_PER_CLUSTER}) AS INT64) AS clustering_index
    FROM ${FOLDER}.CpG_pos_imported p
    JOIN ${FOLDER}.chr_length c ON SAFE_CAST(p.chr AS INT64) = c.chr
    WHERE REGEXP_CONTAINS(p.chr, r'^\\d+$')
    "

# Construct the query to count rows in the table
QUERY="SELECT COUNT(*) as row_count FROM \`${PROJECT_ID}.${FOLDER}.CpG_pos\`"

# Run the query and store the output
OUTPUT=$(bq query --use_legacy_sql=false --format=json "${QUERY}")

# Parse the JSON output to extract the number of rows
ROW_COUNT=$(echo ${OUTPUT} | jq '.[0].row_count')

echo "Number of CpGs in reference genome ${TABLE}: ${ROW_COUNT}"

#--------------------------------------------------------------------------
# Prepare some files related to chromosomes
#--------------------------------------------------------------------------

chmod +x src/make_geomic_regions_for_ref_genome.sh
src/make_geomic_regions_for_ref_genome.sh

echo "Uploading chromosome regions to the bucket"
gzip -f chr_regions.txt
gsutil cp chr_regions.txt.gz gs://"${BUCKET_NAME}"/"${FOLDER}"/chr_regions.txt.gz

# Import the file in BigQuery
bq --location=US load \
               --replace=true \
               --source_format=CSV \
               --field_delimiter "\t" \
               --skip_leading_rows 1 \
                "${FOLDER}".genomic_regions_imported \
               gs://"${BUCKET_NAME}"/"${FOLDER}"/chr_regions.txt.gz \
               chr:INT64,region_inf:INT64,region_sup:INT64,region_center:INT64,clustering_index:INT64


# Cluster the table
bq query --use_legacy_sql=false \
    --location=US \
    "CREATE OR REPLACE TABLE ${FOLDER}.genomic_regions \
    CLUSTER BY clustering_index AS \
    SELECT * FROM ${FOLDER}.genomic_regions_imported"


#--------------------------------------------------------------------------
# Keep regions in the ref genome that have CpGs
#--------------------------------------------------------------------------

echo "Creating the genomic regions with a min number of CpGs for the ref genome"
chmod +x src/find_regions_w_cpg.sh
src/find_regions_in_ref_genome_w_cpg.sh


# Construct the query to count rows in the table
QUERY="SELECT SUM(region_nb_cpg) FROM \`${PROJECT_ID}.${FOLDER}.regions_w_cpg\`"

# Run the query and store the output
OUTPUT=$(bq query --use_legacy_sql=false --format=json "${QUERY}")

# Parse the JSON output to extract the number of rows
NB_CPG=$(echo "${OUTPUT}" | jq '.[0].f0_')

F_NB_CPG=$(format_number_with_comma "${NB_CPG}")
echo "Number of CpGs in reference genome when each region must have ${MIN_NB_CPG_PER_REGION}: ${F_NB_CPG}"

# 20M CpG found.


#--------------------------------------------------------------------------
# Number of CpGs in the 250bp windows.
#--------------------------------------------------------------------------

# 3.7M regions (nb of CpGs >=3), 8.9M (nb of CpGs > 0) vs 12.4M (CpG or not)
# 28.2M CpGs in hg19, 20.8M CpGs (74%) in 250bp windows with at least 3 CpGs.  


#--------------------------------------------------------------------------
# Remove ENCODE's blacklisted regions 
#--------------------------------------------------------------------------

# Link to the bedgraph of the regions:

# https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz

# Unzip the file
gunzip ENCFF001TDO.bed.gz

# Do a bash command to remove "chr":
sed -i 's|chr||g' ENCFF001TDO.bed

# Upload to bucket
gsutil cp ENCFF001TDO.bed gs://${BUCKET}/encode_blacklist.txt

# Transfer to BigQuery
bq --location=US load \
    --replace=true \
    --source_format=CSV \
    --field_delimiter "\t" \
    --skip_leading_rows 1 \
    ${DATASET_EPI}.encode_blacklist \
    gs://${BUCKET}/encode_blacklist.txt \
    chr:STRING,chr_start:INT64,chr_end:INT64,reason:STRING,name1:STRING,name2:STRING

# Remove the genomic regions overlapping with ENCODE's blacklist
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_EPI}.hg19_cpg_regions_${GENOMIC_INTERVAL}bp_black_regions \
    --replace=true \
    "
    SELECT DISTINCT
        t1.chr AS chr,
        t2.chr AS chr_black,
        t1.region_inf,
        t1.region_sup,
        t1.annotate_ref,
        t1.region_nb_cpg
    FROM ${DATASET_EPI}.hg19_cpg_regions_${GENOMIC_INTERVAL}bp t1
    INNER JOIN ${DATASET_EPI}.encode_blacklist t2
    ON t1.chr = t2.chr AND 
       ((chr_start >= t1.region_inf AND chr_start <= t1.region_sup) OR
        (chr_end >= t1.region_inf AND chr_end <= t1.region_sup))
    "

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_EPI}.hg19_cpg_regions_${GENOMIC_INTERVAL}bp_clean \
    --replace=true \
    "
    WITH ALL_DATA AS (
        SELECT 
            t1.chr AS chr,
            t2.chr AS chr_black,
            t1.region_inf,
            t1.region_sup,
            t1.annotate_ref,
            t1.region_nb_cpg
        FROM ${DATASET_EPI}.hg19_cpg_regions_${GENOMIC_INTERVAL}bp t1
        LEFT JOIN ${DATASET_EPI}.hg19_cpg_regions_${GENOMIC_INTERVAL}bp_black_regions t2
        ON t1.region_inf = t2.region_inf AND t1.region_sup = t2.region_sup AND t1.chr = t2.chr
        )
    SELECT chr, region_inf, region_sup, annotate_ref, region_nb_cpg FROM ALL_DATA WHERE chr_black IS NULL
    "

#--------------------------------------------------------------------------
# DNASE track
#--------------------------------------------------------------------------


# URL to download from
# https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=830774571_pra4VNR81N6YjQ3NUyzCQSqI7hiT&clade=mammal&org=Human&db=hg19&hgta_group=regulation&hgta_track=wgEncodeRegDnaseClustered&hgta_table=0&hgta_regionType=genome&position=chr21%3A23%2C031%2C598-43%2C031%2C597&hgta_outputType=wigData&hgta_outFileName=dnase.txt

# Do a bash command to remove "chr":
sed -i 's|chr||g' dnase.txt

# Upload to bucket
gsutil cp dnase.txt gs://${BUCKET}/dnase.txt

# Push DNASe track to BigQuery
bq --location=US load \
    --replace=true \
    --source_format=CSV \
    --field_delimiter "\t" \
    --skip_leading_rows 1 \
    ${DATASET_EPI}.dnase_raw \
    gs://${BUCKET}/dnase.txt \
    bin:INT64,chr:STRING,chr_start:INT64,chr_end:INT64,name:INT64,score:INT64,source_count:FLOAT,source_id:STRING,source_score:STRING

# Clean the database
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_EPI}.dnase \
    --replace=true \
    "
    SELECT 
        chr AS signal_chr, 
        chr_start AS signal_start, 
        chr_end AS signal_end, 
        score
    FROM ${DATASET_EPI}.dnase_raw
    "


#--------------------------------------------------------------------------
# TF BINDING FROM CHIP-SEQ DATA
#--------------------------------------------------------------------------

# Link of the public dataset
# https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=830774571_pra4VNR81N6YjQ3NUyzCQSqI7hiT&clade=mammal&org=Human&db=hg19&hgta_group=regulation&hgta_track=wgEncodeRegTfbsClusteredV2&hgta_table=0&hgta_regionType=genome&position=chr21%3A23%2C031%2C598-43%2C031%2C597&hgta_outputType=primaryTable&hgta_outFileName=encode_ChiP_V2.txt


# Do a bash command to remove "chr":
sed -i 's|chr||g' encode_ChiP_V2.txt

# Upload to bucket
gsutil cp encode_ChiP_V2.txt gs://${BUCKET}/encode_ChiP_V2.txt

# Transfer to BigQuery
bq --location=US load \
    --replace=true \
    --source_format=CSV \
    --field_delimiter "\t" \
    --skip_leading_rows 1 \
    ${DATASET_EPI}.encode_ChiP_V2_raw \
    gs://${BUCKET}/encode_ChiP_V2.txt \
    bin:INT64,chr:STRING,chr_start:INT64,chr_end:INT64,name:STRING,score:INT64,strand:STRING,thick_start:INT64,thick_end:INT64,reserved:INT64,block_count:INT64,block_size:INT64,chrom_start:INT64,exp_count:INT64,exp_id:STRING,exp_score:STRING


bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_EPI}.encode_ChiP_V2 \
    --replace=true \
    "
    SELECT 
        chr AS signal_chr, 
        chr_start AS signal_start, 
        chr_end AS signal_end, 
        score
    FROM ${DATASET_EPI}.encode_ChiP_V2_raw
    "


#--------------------------------------------------------------------------
# TF BINDING MOTIFS known to correlate with ASM
#--------------------------------------------------------------------------

# Provided by table S7 in BiorXiv publication
# Publication link: https://www.biorxiv.org/content/10.1101/815605v3
# Table link: https://www.biorxiv.org/content/biorxiv/early/2020/04/07/815605/DC8/embed/media-8.xlsx?download=true

# Save the XLS file into a txt file.

# Motifs known to correlate with ASM (from bioRiv publication)
gsutil cp asm_motifs.txt gs://${BUCKET}/asm_motifs.txt


# Upload known ASM motifs to BigQuery
bq --location=US load \
    --replace=true \
    --source_format=CSV \
    --field_delimiter "\t" \
    --skip_leading_rows 1 \
    ${DATASET_EPI}.asm_motifs_raw \
    gs://${BUCKET}/asm_motifs.txt \
    motif:STRING,n_tot:INT64,n_asm:INT64,n_no_asm:INT64,odds_ratio:FLOAT,p_value:FLOAT,fdr:FLOAT

# We keep the motifs where the OR > 1
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_EPI}.asm_motifs \
    --replace=true \
    "
    SELECT motif AS asm_motif 
    FROM ${DATASET_EPI}.asm_motifs_raw
    WHERE odds_ratio > 1
    "

#--------------------------------------------------------------------------
# TF BINDING MOTIFS
#--------------------------------------------------------------------------

# Provided by Catherine.
# Originally obtained at http://compbio.mit.edu/encode-motifs/

# Clean the database of motifs
mv kherad_tf_sorted.bed kherad_tf_sorted.txt
sed -i 's|chr||g' kherad_tf_sorted.txt

# Upload database to bucket
gsutil cp kherad_tf_sorted.txt gs://${BUCKET}/kherad_tf_sorted.txt

# Transfer bucket -> BigQuery
bq --location=US load \
    --replace=true \
    --source_format=CSV \
    --field_delimiter "\t" \
    --skip_leading_rows 0 \
    ${DATASET_EPI}.kherad_tf_sorted \
    gs://${BUCKET}/kherad_tf_sorted.txt \
    chr:STRING,chr_start:INT64,chr_end:INT64,motif:STRING

# Keep the motifs known to correlate with ASM
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_EPI}.tf_motifs \
    --replace=true \
    "
    WITH 
        ASM_MOTIFS AS (
            SELECT * 
            FROM ${DATASET_EPI}.asm_motifs
        ),
        KHERAD AS (
            SELECT * FROM ${DATASET_EPI}.kherad_tf_sorted
        )
        SELECT 
            chr AS signal_chr, 
            chr_start AS signal_start, 
            chr_end AS signal_end, 
            motif AS score 
        FROM KHERAD
        INNER JOIN ASM_MOTIFS
        ON asm_motif = motif
    "


#--------------------------------------------------------------------------
# Annotate the reference genome for the various epigenetic tracks
#--------------------------------------------------------------------------

# To use this script, you need to use a table where the chr is 'chr' and 
# the middle of the region as 'annotate_ref'

SAMPLE="hg19_cpg_regions_"${GENOMIC_INTERVAL}"bp_clean"

# Create genomic regions used to split jobs per chromosome 
INTERVAL="60000000"

# Prepare TSV file per chromosome (used for many jobs)
echo -e "--env TABLE\t--env EPI_SIGNAL\tCHR\tLOWER_B\tUPPER_B" > chr_split_epi.tsv

for SIGNAL in "dnase" "encode_ChiP_V2" "tf_motifs" ; do 
    for CHR in `seq 1 22` X Y ; do
            echo "Processing chromosome" ${CHR}
            NUCLEOTIDES_IN_CHR=$(awk -v CHR="${CHR}" -F"\t" '{ if ($1 == CHR) print $2}' chr_length.txt)
            INF="1"
            SUP=$(( $NUCLEOTIDES_IN_CHR<$INTERVAL ? $NUCLEOTIDES_IN_CHR: $INTERVAL ))
            echo -e "${SAMPLE}\t${SIGNAL}\t${CHR}\t$INF\t$SUP" >> chr_split_epi.tsv # for jobs
            while [ $NUCLEOTIDES_IN_CHR -gt $SUP ] ; do
                INCREMENT=$(( $NUCLEOTIDES_IN_CHR-$SUP<$INTERVAL ? $NUCLEOTIDES_IN_CHR-$SUP: $INTERVAL ))
                INF=$(( ${SUP} + 1 ))
                SUP=$(( ${SUP} + $INCREMENT ))
                echo -e "${SAMPLE}\t${SIGNAL}\t${CHR}\t$INF\t$SUP" >> chr_split_epi.tsv

            done
    done
done

# Run the 192 jobs. Only 100 at a time (BQ limit)
dsub \
    --project $PROJECT_ID \
    --zones $ZONE_ID \
    --image ${DOCKER_GCP} \
    --logging $LOG \
    --env DATASET_EPI="${DATASET_EPI}" \
    --env EPI_REGION="${EPI_REGION}" \
    --env DATASET="${DATASET_EPI}" \
    --script ${SCRIPTS}/annotation.sh \
    --tasks chr_split_epi.tsv 1-100 \
    --wait

dsub \
    --project $PROJECT_ID \
    --zones $ZONE_ID \
    --image ${DOCKER_GCP} \
    --logging $LOG \
    --env DATASET_EPI="${DATASET_EPI}" \
    --env EPI_REGION="${EPI_REGION}" \
    --env DATASET="${DATASET_EPI}" \
    --script ${SCRIPTS}/annotation.sh \
    --tasks chr_split_epi.tsv 101-192 \
    --wait


# Delete previous files in case
for EPI_SIGNAL in "dnase" "encode_ChiP_V2" "tf_motifs" ; do
    echo "epi signal: " ${EPI_SIGNAL}
    bq rm -f -t ${DATASET_EPI}.${SAMPLE}_${EPI_SIGNAL}_all
done

# Concatenate files
{ read
while read SAMPLE SIGNAL CHR LOWER_B UPPER_B ; do 
    echo "Sample: " ${SAMPLE} ", signal: " ${SIGNAL} ", chr:" ${CHR} ", lower:" ${LOWER_B} "and upper:" ${UPPER_B}
    bq cp --append_table \
        ${DATASET_EPI}.${SAMPLE}_${SIGNAL}_${CHR}_${LOWER_B}_${UPPER_B} \
        ${DATASET_EPI}.${SAMPLE}_${SIGNAL}_all
done 
} < chr_split_epi.tsv


# Delete intermediary files
{ read
while read SAMPLE SIGNAL CHR LOWER_B UPPER_B ; do 
    echo "Sample: " ${SAMPLE} ", signal: " ${SIGNAL} ", chr:" ${CHR} ", lower:" ${LOWER_B} "and upper:" ${UPPER_B}
    bq rm -f -t ${DATASET_EPI}.${SAMPLE}_${SIGNAL}_${CHR}_${LOWER_B}_${UPPER_B} 
done 
} < chr_split_epi.tsv

#---------------------------------

SAMPLE="hg19_cpg_regions_"${GENOMIC_INTERVAL}"bp_clean"


#--------------------
# The following code estimates if there is any of the epigenetic marks at all
# Will label with 0 and 1 each genomic region for each epigenetic mark

# Gather all scores in a single array per region. This creates as many tables as there
# are epigenetic signals for annotation.
# for EPI_SIGNAL in "dnase" "encode_ChiP_V2" "tf_motifs" ; do
#     echo "Processing the signal " ${EPI_SIGNAL}
#     bq query \
#         --use_legacy_sql=false \
#         --destination_table ${DATASET_EPI}.${SAMPLE}_${EPI_SIGNAL}_binary \
#         --replace=true \
#         "
#         WITH 
#             EPI_AGG AS ( -- we group the DNASe scores together
#                 SELECT 
#                     * EXCEPT(score),
#                     ARRAY_AGG(STRUCT(score)) AS epi
#                 FROM ${DATASET_EPI}.${SAMPLE}_${EPI_SIGNAL}_all
#                 GROUP BY 
#                     chr,
#                     region_inf,
#                     region_sup,
#                     annotate_ref,
#                     region_nb_cpg
#             )
#             SELECT 
#                 * EXCEPT(epi),
#                 -- the command below takes care of the case if there is no  score in the array
#                 IF(
#                     ARRAY_LENGTH(
#                         (SELECT ARRAY (
#                             SELECT score 
#                             FROM UNNEST(epi) 
#                             WHERE score IS NOT NULL
#                             )
#                         )) > 0, 1, 0
#                 ) AS ${EPI_SIGNAL}
#             FROM EPI_AGG
#         "
# done


#--------------------
# The following code estimates the number of epigenetic marks per genomic region per epigenetic mark
# Will label with 0 and 1 each genomic region for each epigenetic mark
# The table of binary (yes/no) can be created from that table.


for EPI_SIGNAL in "dnase" "encode_ChiP_V2" "tf_motifs" ; do
    echo "Processing the signal " ${EPI_SIGNAL}
    bq query \
        --use_legacy_sql=false \
        --destination_table ${DATASET_EPI}.${SAMPLE}_${EPI_SIGNAL} \
        --replace=true \
        "
        SELECT 
            * EXCEPT(score),
            COUNT(score) AS ${EPI_SIGNAL}
        FROM ${DATASET_EPI}.${SAMPLE}_${EPI_SIGNAL}_all
        GROUP BY 
            chr,
            region_inf,
            region_sup,
            annotate_ref,
            region_nb_cpg
        "
done

## Look at the data
bq query \
    --use_legacy_sql=false \
    "
    SELECT DISTINCT(dnase) 
    FROM ${DATASET_EPI}.${SAMPLE}_dnase
    LIMIT 10
    "


# Merge the tables of epigenetic signals.

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_EPI}.${SAMPLE}_annotated \
    --replace=true \
    "
    SELECT
        t1.chr AS chr,
        t1.region_inf AS region_inf,
        t1.region_sup AS region_sup,
        t1.region_nb_cpg AS region_nb_cpg,
        t1.dnase AS dnase,
        t2.encode_ChiP_V2 AS encode_ChiP_V2,
        t3.tf_motifs AS tf_motifs
    FROM ${DATASET_EPI}.${SAMPLE}_dnase t1
    JOIN ${DATASET_EPI}.${SAMPLE}_encode_ChiP_V2 t2 
    ON t1.chr = t2.chr AND 
        t1.region_inf = t2.region_inf AND 
        t1.region_sup = t2.region_sup
    JOIN ${DATASET_EPI}.${SAMPLE}_tf_motifs t3 
    ON t1.chr = t3.chr AND 
    t1.region_inf = t3.region_inf AND 
    t1.region_sup = t3.region_sup
    "

## Look at the data
bq query \
    --use_legacy_sql=false \
    "
    SELECT *
    FROM ${DATASET_EPI}.${SAMPLE}_annotated
    LIMIT 10
    "