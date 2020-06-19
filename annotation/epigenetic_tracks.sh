#--------------------------------------------------------------------------
# DNASE track
#--------------------------------------------------------------------------


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
    ${DATASET_EPI}.dnase_raw \
    gs://${CLOUDASM_BUCKET}/dnase.txt \
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
gsutil cp encode_ChiP_V2.txt gs://${CLOUDASM_BUCKET}/encode_ChiP_V2.txt

# Transfer to BigQuery
bq --location=US load \
    --replace=true \
    --source_format=CSV \
    --field_delimiter "\t" \
    --skip_leading_rows 1 \
    ${DATASET_EPI}.encode_ChiP_V2_raw \
    gs://${CLOUDASM_BUCKET}/encode_ChiP_V2.txt \
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
gsutil cp asm_motifs.txt gs://${CLOUDASM_BUCKET}/asm_motifs.txt


# Upload known ASM motifs to BigQuery
bq --location=US load \
    --replace=true \
    --source_format=CSV \
    --field_delimiter "\t" \
    --skip_leading_rows 1 \
    ${DATASET_EPI}.asm_motifs_raw \
    gs://${CLOUDASM_BUCKET}/asm_motifs.txt \
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
gsutil cp kherad_tf_sorted.txt gs://${CLOUDASM_BUCKET}/kherad_tf_sorted.txt

# Transfer bucket -> BigQuery
bq --location=US load \
    --replace=true \
    --source_format=CSV \
    --field_delimiter "\t" \
    --skip_leading_rows 0 \
    ${DATASET_EPI}.kherad_tf_sorted \
    gs://${CLOUDASM_BUCKET}/kherad_tf_sorted.txt \
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

