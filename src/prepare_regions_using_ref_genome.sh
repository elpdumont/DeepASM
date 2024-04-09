#!/bin/bash

# Requirements
# CpG bed file in bucket/ref_genomes/reference_genome

echo "Creating the dataset with an expiration"
if bq ls --project_id="${PROJECT_ID}" | grep -w "${REFG_FOLDER}"; then
	echo "Dataset ${REFG_FOLDER} already exists in project ${PROJECT_ID}."
else
	# Create the dataset since it does not exist
	bq mk --project_id="${PROJECT_ID}" --dataset --default_table_expiration="${REF_GENOME_DATASET_EXPIRATION_SEC}" "${PROJECT_ID}:${REFG_FOLDER}"
	echo "Dataset ${REFG_FOLDER} created in project ${PROJECT_ID}."
fi

chmod +x src/make_chr_info_for_ref_genome.sh
src/make_chr_info_for_ref_genome.sh

gsutil cp chr_length.txt gs://"${BUCKET_NAME}"/"${REFG_FOLDER}"/chr_length.txt

bq --location=US load \
               --replace=true \
               --source_format=CSV \
               --field_delimiter "\t" \
               --skip_leading_rows 1 \
                "${REFG_FOLDER}".chr_length \
               gs://"${BUCKET_NAME}"/"${REFG_FOLDER}"/chr_length.txt \
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
                "${REFG_FOLDER}".CpG_pos_imported \
               gs://"${BUCKET_NAME}"/ref_genomes/"${REFERENCE_GENOME}"/"${REFERENCE_GENOME}"_CpG_pos.bed \
               chr:STRING,region_inf:INT64,region_sup:INT64

# Converting chr to INT for easy partitioning
bq query \
    --use_legacy_sql=false \
    "
    CREATE OR REPLACE TABLE ${REFG_FOLDER}.CpG_pos
    CLUSTER BY clustering_index AS
    SELECT
        CAST(p.chr AS INT64) AS chr,
        p.region_inf,
        p.region_sup,
        CAST(FLOOR((c.absolute_nucleotide_pos + p.region_inf) / ${NB_NUCLEOTIDES_PER_CLUSTER}) AS INT64) AS clustering_index
    FROM ${REFG_FOLDER}.CpG_pos_imported p
    JOIN ${REFG_FOLDER}.chr_length c ON SAFE_CAST(p.chr AS INT64) = c.chr
    WHERE REGEXP_CONTAINS(p.chr, r'^\\d+$')
    "

# Construct the query to count rows in the table
QUERY="SELECT COUNT(*) as row_count FROM \`${PROJECT_ID}.${REFG_FOLDER}.CpG_pos\`"

# Run the query and store the output
OUTPUT=$(bq query --use_legacy_sql=false --format=json "${QUERY}")

# Parse the JSON output to extract the number of rows
ROW_COUNT=$(echo ${OUTPUT} | jq '.[0].row_count')
F_ROW_COUNT=$(format_number_with_comma "${ROW_COUNT}")

echo "Number of CpGs in reference genome ${TABLE}: ${F_ROW_COUNT}"

#--------------------------------------------------------------------------
# Prepare some files related to chromosomes
#--------------------------------------------------------------------------

chmod +x src/make_genomic_regions_for_ref_genome.sh
src/make_genomic_regions_for_ref_genome.sh

echo "Uploading chromosome regions to the bucket"
gzip -f chr_regions.txt
gsutil cp chr_regions.txt.gz gs://"${BUCKET_NAME}"/"${REFG_FOLDER}"/chr_regions.txt.gz

# Import the file in BigQuery
bq --location=US load \
               --replace=true \
               --source_format=CSV \
               --field_delimiter "\t" \
               --skip_leading_rows 1 \
                "${REFG_FOLDER}".genomic_regions_imported \
               gs://"${BUCKET_NAME}"/"${REFG_FOLDER}"/chr_regions.txt.gz \
               chr:INT64,region_inf:INT64,region_sup:INT64,region_center:INT64,clustering_index:INT64


# Cluster the table
bq query --use_legacy_sql=false \
    --location=US \
    "CREATE OR REPLACE TABLE ${REFG_FOLDER}.genomic_regions \
    CLUSTER BY clustering_index AS \
    SELECT * FROM ${REFG_FOLDER}.genomic_regions_imported"

NB_REGIONS=$(execute_query "SELECT COUNT(*) FROM ${PROJECT_ID}.${REFG_FOLDER}.genomic_regions_imported")

echo "BEFORE FILTERING, found ${NB_REGIONS} regions with a min of ${MIN_NB_CPG_PER_REGION_IN_REF_GENOME} CpGse"


#--------------------------------------------------------------------------
# Keep regions in the ref genome that have CpGs
#--------------------------------------------------------------------------

echo "Creating the genomic regions with a min number of CpGs for the ref genome"
chmod +x src/find_regions_w_cpg.sh
src/find_regions_in_ref_genome_w_cpg.sh

# Construct the query to count rows in the table
NB_CPG=$(execute_query "SELECT SUM(region_nb_cpg) FROM ${PROJECT_ID}.${REFG_FOLDER}.regions_w_cpg")
NB_REGIONS=$(execute_query "SELECT COUNT(*) FROM ${PROJECT_ID}.${REFG_FOLDER}.regions_w_cpg")

echo "Found:${NB_REGIONS} regions with a min of ${MIN_NB_CPG_PER_REGION_IN_REF_GENOME} CpGs, totalling ${NB_CPG} across the whole genome"


echo "Number of regions in reference genome when each region must have ${MIN_NB_CPG_PER_REGION_IN_REF_GENOME}: ${F_NB_CPG}"

#--------------------------------------------------------------------------
# Remove ENCODE's blacklisted regions 
#--------------------------------------------------------------------------

echo "Downloading the ENCODE blacklisted regions"
wget -O encode_blacklist_regions.bed.gz "${ENCODE_BLACKLIST_REGIONS_URL}"
gunzip encode_blacklist_regions.bed.gz

# Do a bash command to remove "chr":
sed -i '' 's|chr||g' encode_blacklist_regions.bed

echo "Uploading blacklisted regions to Cloud Storage"
gsutil cp encode_blacklist_regions.bed gs://${BUCKET_NAME}/${REFG_FOLDER}/

echo "Transfering these regions to BigQuery"
bq --location=US load \
    --replace=true \
    --source_format=CSV \
    --field_delimiter "\t" \
    --skip_leading_rows 1 \
    ${REFG_FOLDER}.encode_blacklist_regions_imported \
    gs://${BUCKET_NAME}/${REFG_FOLDER}/encode_blacklist_regions.bed \
    chr:STRING,chr_start:INT64,chr_end:INT64,reason:STRING,name1:STRING,name2:STRING

# Enforce chr as INT65
bq query \
    --use_legacy_sql=false \
    --destination_table "${REFG_FOLDER}".encode_blacklist_regions \
    --replace=true \
    "
    SELECT
        CAST(chr AS INT64) AS chr,
        chr_start,
        chr_end,
        reason
    FROM ${REFG_FOLDER}.encode_blacklist_regions_imported
    WHERE SAFE_CAST(chr AS INT64) IS NOT NULL
    "

echo "Identify regions of the ref genome that overlap with ENCODE blacklisted regions"
bq query \
    --use_legacy_sql=false \
    --destination_table "${REFG_FOLDER}".regions_w_cpg_in_blacklisted_regions \
    --replace=true \
    "
    SELECT DISTINCT
        t1.chr AS chr,
        t2.chr AS chr_black,
        t1.region_inf,
        t1.region_sup,
        t1.region_center,
        t1.region_nb_cpg,
        t1.clustering_index
    FROM ${REFG_FOLDER}.regions_w_cpg t1
    INNER JOIN ${REFG_FOLDER}.encode_blacklist_regions t2
    ON t1.chr = t2.chr AND 
       ((t2.chr_start >= t1.region_inf AND t2.chr_start <= t1.region_sup) OR
        (t2.chr_end >= t1.region_inf AND t2.chr_end <= t1.region_sup))
    "

NB_REGIONS=$(execute_query "SELECT COUNT(*) FROM ${PROJECT_ID}.${REFG_FOLDER}.regions_w_cpg_in_blacklisted_regions")
echo "There are ${NB_REGIONS} genomic regions in blacklisted regions, which we will remove"


echo "Create a table of genomic regions that do not overlap with ENCODE blacklisted regions"
bq query \
    --use_legacy_sql=false \
    --destination_table "${REFG_FOLDER}".regions_w_cpg_wo_blacklisted_regions \
    --replace=true \
    "
    WITH ALL_DATA AS (
        SELECT 
            t1.chr AS chr,
            t2.chr AS chr_black,
            t1.region_inf,
            t1.region_sup,
            t1.region_center,
            t1.region_nb_cpg,
            t1.clustering_index
        FROM ${REFG_FOLDER}.regions_w_cpg t1
        LEFT JOIN ${REFG_FOLDER}.regions_w_cpg_in_blacklisted_regions t2
        ON t1.region_inf = t2.region_inf AND t1.region_sup = t2.region_sup AND t1.chr = t2.chr
        )
    SELECT chr, region_inf, region_sup, region_center, region_nb_cpg, clustering_index
    FROM ALL_DATA
    WHERE chr_black IS NULL
    "

# Construct the query to count rows in the table
NB_CPG=$(execute_query "SELECT SUM(region_nb_cpg) FROM ${PROJECT_ID}.${REFG_FOLDER}.regions_w_cpg_wo_blacklisted_regions")
NB_REGIONS=$(execute_query "SELECT COUNT(*) FROM ${PROJECT_ID}.${REFG_FOLDER}.regions_w_cpg_wo_blacklisted_regions")

echo "Found:${NB_REGIONS} regions with a min of ${MIN_NB_CPG_PER_REGION_IN_REF_GENOME} CpGs, totalling ${NB_CPG} across the whole genome"