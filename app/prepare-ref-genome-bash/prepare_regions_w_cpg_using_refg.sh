#!/bin/bash

# Requirements
# CpG bed file in bucket/ref_genomes/reference_genome

script_folder="app/prepare-ref-genome-bash"


echo "Create a file with all the lengths of all chromosomes"
"${script_folder}"/make_chr_info_for_ref_genome.sh

echo "Copying this file to Cloud Storage"
gsutil cp chr_length.txt gs://"${BUCKET_NAME}"/"${REFG_DATASET}"/chr_length.txt

echo "Importing this file to BigQuery"
bq --location="${REGION}" load \
               --replace=true \
               --source_format=CSV \
               --field_delimiter "\t" \
               --skip_leading_rows 1 \
                "${PROJECT}":"${REFG_DATASET}".chr_length \
               gs://"${BUCKET_NAME}"/"${REFG_DATASET}"/chr_length.txt \
               chr:INT64,chr_length:INT64,absolute_nucleotide_pos:INT64

#--------------------------------------------------------------------------
# Upload a file of the CpG positions in the reference genome
#--------------------------------------------------------------------------

echo "Transfer the CpG positions from ref genome to the dataset"
bq --location="${REGION}" load \
               --replace=true \
               --source_format=CSV \
               --field_delimiter "\t" \
               --skip_leading_rows 1 \
                "${PROJECT}":"${REFG_DATASET}".CpG_pos_imported \
               gs://"${BUCKET_NAME}"/ref_genomes/"${REFERENCE_GENOME}"/"${REFERENCE_GENOME}"_CpG_pos.bed \
               chr:STRING,region_inf:INT64,region_sup:INT64

echo "Converting chr to INT for easy partitioning"
bq query \
    --destination_table="${PROJECT}":"${REFG_DATASET}".CpG_pos \
    --replace=true \
    --use_legacy_sql=false \
    --clustering_fields=chr \
    --range_partitioning=clustering_index,0,4000,1 \
    "
    SELECT
        CAST(p.chr AS INT64) AS chr,
        p.region_inf,
        p.region_sup,
        CAST(FLOOR((c.absolute_nucleotide_pos + p.region_inf) / ${NB_NUCLEOTIDES_PER_CLUSTER}) AS INT64) AS clustering_index
    FROM ${REFG_DATASET}.CpG_pos_imported p
    JOIN ${REFG_DATASET}.chr_length c ON SAFE_CAST(p.chr AS INT64) = c.chr
    WHERE REGEXP_CONTAINS(p.chr, r'^\\d+$')
    "

ROW_COUNT=$(execute_query "SELECT COUNT(*) as row_count FROM \`${PROJECT}.${REFG_DATASET}.CpG_pos\`")

echo "Number of CpGs in reference genome ${TABLE}: ${ROW_COUNT}"

#--------------------------------------------------------------------------
# Prepare some files related to chromosomes
#--------------------------------------------------------------------------

"${script_folder}"/make_regions_for_ref_genome.sh

echo "Uploading chromosome regions to the bucket"
gzip -f chr_regions.txt
gsutil cp chr_regions.txt.gz gs://"${BUCKET_NAME}"/"${REFG_DATASET}"/chr_regions.txt.gz

echo "Import the file of chr regions in BigQuery"

bq --location="${REGION}" load \
               --replace=true \
               --source_format=CSV \
               --field_delimiter "\t" \
               --skip_leading_rows 1 \
                "${PROJECT}":"${REFG_DATASET}".regions_imported \
               gs://"${BUCKET_NAME}"/"${REFG_DATASET}"/chr_regions.txt.gz \
               chr:INT64,region_inf:INT64,region_sup:INT64,region_center:INT64,clustering_index:INT64


# Cluster the table
bq query \
    --destination_table="${REFG_DATASET}".regions \
    --replace=true \
    --use_legacy_sql=false \
    --clustering_fields=chr \
    --range_partitioning=clustering_index,0,4000,1 \
    "
    SELECT * FROM ${REFG_DATASET}.regions_imported
    "


NB_REGIONS=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${REFG_DATASET}.regions_imported")

echo "BEFORE FILTERING, found ${NB_REGIONS} regions with a min of ${MIN_NB_CPG_PER_REGION_IN_REF_GENOME} CpGs"


#--------------------------------------------------------------------------
# Keep regions in the ref genome that have CpGs
#--------------------------------------------------------------------------

echo "Creating the genomic regions with a min number of CpGs for the ref genome"
"${script_folder}"/find_regions_in_ref_genome_w_cpg.sh

# Construct the query to count rows in the table
NB_CPG=$(execute_query "SELECT SUM(region_nb_cpg) FROM ${PROJECT}.${REFG_DATASET}.regions_w_cpg")
NB_REGIONS=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${REFG_DATASET}.regions_w_cpg")

echo "Found:${NB_REGIONS} regions with a min of ${MIN_NB_CPG_PER_REGION_IN_REF_GENOME} CpGs, totalling ${NB_CPG} across the whole genome"


#--------------------------------------------------------------------------
# Remove ENCODE's blacklisted regions 
#--------------------------------------------------------------------------

echo "Downloading the ENCODE blacklisted regions"
wget -O encode_blacklist_regions.bed.gz "${ENCODE_BLACKLIST_REGIONS_URL}"
gunzip encode_blacklist_regions.bed.gz

# Do a bash command to remove "chr":
sed -i '' 's|chr||g' encode_blacklist_regions.bed

echo "Uploading blacklisted regions to Cloud Storage"
gsutil cp encode_blacklist_regions.bed gs://"${BUCKET_NAME}"/"${REFG_DATASET}"/

echo "Transfering these regions to BigQuery"
bq --location="${REGION}" load \
    --replace=true \
    --source_format=CSV \
    --field_delimiter "\t" \
    --skip_leading_rows 1 \
    "${PROJECT}":"${REFG_DATASET}".encode_blacklist_regions_imported \
    gs://"${BUCKET_NAME}"/"${REFG_DATASET}"/encode_blacklist_regions.bed \
    chr:STRING,chr_start:INT64,chr_end:INT64,reason:STRING,name1:STRING,name2:STRING

# Enforce chr as INT65
bq query \
    --use_legacy_sql=false \
    --destination_table "${PROJECT}":"${REFG_DATASET}".encode_blacklist_regions \
    --replace=true \
    "
    SELECT
        CAST(chr AS INT64) AS chr,
        chr_start,
        chr_end,
        reason
    FROM ${REFG_DATASET}.encode_blacklist_regions_imported
    WHERE SAFE_CAST(chr AS INT64) IS NOT NULL
    "

echo "Identify regions of the ref genome that overlap with ENCODE blacklisted regions"
bq query \
    --use_legacy_sql=false \
    --destination_table "${PROJECT}":"${REFG_DATASET}".regions_w_cpg_no_blacklist \
    --replace=true \
    --clustering_fields=chr \
    --range_partitioning=clustering_index,0,4000,1 \
    "
    WITH REGIONS_BLACK AS (
        SELECT DISTINCT
            t1.chr AS chr,
            t2.chr AS chr_black,
            t1.region_inf,
            t1.region_sup,
            t1.region_center,
            t1.region_nb_cpg,
            t1.clustering_index
        FROM ${REFG_DATASET}.regions_w_cpg t1
        INNER JOIN ${REFG_DATASET}.encode_blacklist_regions t2
        ON t1.chr = t2.chr AND
        ((t2.chr_start >= t1.region_inf AND t2.chr_start <= t1.region_sup) OR
            (t2.chr_end >= t1.region_inf AND t2.chr_end <= t1.region_sup))
        ),
    ALL_DATA AS (
        SELECT 
            t1.chr AS chr,
            t2.chr AS chr_black,
            t1.region_inf,
            t1.region_sup,
            t1.region_center,
            t1.region_nb_cpg,
            t1.clustering_index
        FROM ${REFG_DATASET}.regions_w_cpg t1
        LEFT JOIN REGIONS_BLACK t2
        ON t1.region_inf = t2.region_inf AND t1.region_sup = t2.region_sup AND t1.chr = t2.chr
        )
    SELECT chr, region_inf, region_sup, region_center, region_nb_cpg, clustering_index
    FROM ALL_DATA
    WHERE chr_black IS NULL
    "


# Construct the query to count rows in the table
NB_CPG=$(execute_query "SELECT SUM(region_nb_cpg) FROM ${PROJECT}.${REFG_DATASET}.regions_w_cpg_no_blacklist")
NB_REGIONS=$(execute_query "SELECT COUNT(*) FROM ${PROJECT}.${REFG_DATASET}.regions_w_cpg_no_blacklist")

echo "Found:${NB_REGIONS} regions with a min of ${MIN_NB_CPG_PER_REGION_IN_REF_GENOME} CpGs, totalling ${NB_CPG} across the whole genome"