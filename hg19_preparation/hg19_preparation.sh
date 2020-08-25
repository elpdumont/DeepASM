# Here we prepare the ref genome with where at least 3 CpGs are present.
# Then we annotate the windows with at least 3 CpGs present.

#--------------------------------------------------------------------------
# Variables
#--------------------------------------------------------------------------

SCRIPTS="/Users/emmanuel/GITHUB_REPOS/DeepASM/hg19_preparation"

# Dataset where the ref genome and the epigenetic signals will be located
DATASET_EPI="hg19"

# Size of genomic regions:
GENOMIC_INTERVAL="250"

# Desired window for annotation analysis (needs to be half of INTERVAL)
EPI_REGION=$(( ${GENOMIC_INTERVAL} / 2 ))

# Cloud Storage location of the logs
LOG="gs://cloudasm-encode/logging/deepasm"

# Storage bucket for hg19-related files.
BUCKET="cloudasm-encode"

# Docker file required to run the scripts
DOCKER_GCP="google/cloud-sdk:255.0.0"

# GCP global variables
PROJECT_ID="hackensack-tyco"
REGION_ID="us-central1"
ZONE_ID="us-central1-b"


#--------------------------------------------------------------------------
# Upload a file of the CpG positions in the reference genome
#--------------------------------------------------------------------------

# Provided by Tycko lab @ gs://ref_genomes/grc37/hg19_CpG_pos.bed

# Transfer the CpG positions to the dataset
bq --location=US load \
               --replace=true \
               --source_format=CSV \
               --field_delimiter "\t" \
               --skip_leading_rows 1 \
                ${DATASET_EPI}.hg19_CpG_pos \
               gs://ref_genomes/grc37/hg19_CpG_pos.bed \
               chr_region:STRING,region_inf:INT64,region_sup:INT64


#--------------------------------------------------------------------------
# Split the human genome into windows with the same number of basepairs.
#--------------------------------------------------------------------------

# We enter all the chromosome lengths in base pairs for humans

echo -e "chr\tlength" > chr_length.txt
echo -e "1\t249250621" > chr_length.txt && echo -e "2\t243199373" >> chr_length.txt && echo -e "3\t198022430" >> chr_length.txt \
&& echo -e "4\t191154276" >> chr_length.txt && echo -e "5\t180915260" >> chr_length.txt && echo -e "6\t171115067" >> chr_length.txt \
&& echo -e "7\t159138663" >> chr_length.txt && echo -e "8\t146364022" >> chr_length.txt && echo -e "9\t141213431" >> chr_length.txt \
&& echo -e "10\t135534747" >> chr_length.txt && echo -e "11\t135006516" >> chr_length.txt && echo -e "12\t133851895" >> chr_length.txt \
&& echo -e "13\t115169878" >> chr_length.txt && echo -e "14\t107349540" >> chr_length.txt && echo -e "15\t102531392" >> chr_length.txt \
&& echo -e "16\t90354753" >> chr_length.txt && echo -e "17\t81195210" >> chr_length.txt && echo -e "18\t78077248" >> chr_length.txt \
&& echo -e "19\t59128983" >> chr_length.txt && echo -e "20\t63025520" >> chr_length.txt && echo -e "21\t48129895" >> chr_length.txt \
&& echo -e "22\t51304566" >> chr_length.txt && echo -e "X\t155270560 " >> chr_length.txt && echo -e "Y\t59373566" >> chr_length.txt


# Prepare TSV file per chromosome (used for many jobs)
echo -e "chr\tregion_inf\tregion_sup\tannotate_ref" > chr_regions.txt

# Create the windows with $GENOMIC_INTERVALS base pairs in it
for CHR in `seq 1 22` X Y ; do
    echo "Processing chromosome" ${CHR}
    NUCLEOTIDES_IN_CHR=$(awk -v CHR="${CHR}" -F"\t" '{ if ($1 == CHR) print $2}' chr_length.txt)
    INF="1"
    SUP=$(( $NUCLEOTIDES_IN_CHR<$GENOMIC_INTERVAL ? $NUCLEOTIDES_IN_CHR: $GENOMIC_INTERVAL ))
    ANNOTATE_REF=$(( ($INF + $SUP)/2 ))
    echo -e "${CHR}\t$INF\t$SUP\t$ANNOTATE_REF" >> chr_regions.txt # for jobs
    while [ $NUCLEOTIDES_IN_CHR -gt $SUP ] ; do
        INCREMENT=$(( $NUCLEOTIDES_IN_CHR-$SUP<$GENOMIC_INTERVAL ? $NUCLEOTIDES_IN_CHR-$SUP: $GENOMIC_INTERVAL ))
        INF=$(( ${SUP} + 1 ))
        SUP=$(( ${SUP} + $INCREMENT ))
        ANNOTATE_REF=$(( ($INF + $SUP)/2 ))
        echo -e "${CHR}\t$INF\t$SUP\t$ANNOTATE_REF" >> chr_regions.txt

    done
done

# Upload the file to a bucket
gsutil cp chr_regions.txt gs://${BUCKET}/chr_regions_${GENOMIC_INTERVAL}bp.txt

# Import the file in BigQuery
bq --location=US load \
               --replace=true \
               --source_format=CSV \
               --field_delimiter "\t" \
               --skip_leading_rows 1 \
                ${DATASET_EPI}.hg19_regions_${GENOMIC_INTERVAL}bp \
               gs://${BUCKET}/chr_regions_${GENOMIC_INTERVAL}bp.txt \
               chr_region:STRING,region_inf:INT64,region_sup:INT64,annotate_ref:INT64


#--------------------------------------------------------------------------
# Keep regions that have CpGs (we use hg19 reference genome)
#--------------------------------------------------------------------------

# Create genomic regions used to split jobs per chromosome 
INTERVAL="100000000"

# Prepare TSV file per chromosome (used for many jobs)
echo -e "CHR\tLOWER_B\tUPPER_B" > chr_split_hg19.tsv

for CHR in `seq 1 22` X Y ; do
        echo "Processing chromosome" ${CHR}
        NUCLEOTIDES_IN_CHR=$(awk -v CHR="${CHR}" -F"\t" '{ if ($1 == CHR) print $2}' chr_length.txt)
        INF="1"
        SUP=$(( $NUCLEOTIDES_IN_CHR<$INTERVAL ? $NUCLEOTIDES_IN_CHR: $INTERVAL ))
        echo -e "${CHR}\t$INF\t$SUP" >> chr_split_hg19.tsv # for jobs
        while [ $NUCLEOTIDES_IN_CHR -gt $SUP ] ; do
            INCREMENT=$(( $NUCLEOTIDES_IN_CHR-$SUP<$INTERVAL ? $NUCLEOTIDES_IN_CHR-$SUP: $INTERVAL ))
            INF=$(( ${SUP} + 1 ))
            SUP=$(( ${SUP} + $INCREMENT ))
            echo -e "${CHR}\t$INF\t$SUP" >> chr_split_hg19.tsv

        done
    done

# Launch parallel jobs to establish hg19 regions with CpGs
dsub \
--project $PROJECT_ID \
--zones $ZONE_ID \
--image ${DOCKER_GCP} \
--logging $LOG \
--env DATASET_EPI="${DATASET_EPI}" \
--env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
--script ${SCRIPTS}/hg19_cpg_regions.sh \
--tasks chr_split_hg19.tsv \
--wait

# Concatenate in a single table all hg19 regions with CpGs
bq rm -f -t ${DATASET_EPI}.hg19_cpg_regions_${GENOMIC_INTERVAL}bp

{ read
while read CHR LOWER_B UPPER_B ; do 
    echo "Chromosome is " ${CHR} ", lower:" ${LOWER_B} ", and upper:" ${UPPER_B}
    bq cp --append_table \
        ${DATASET_EPI}.hg19_cpg_regions_${CHR}_${LOWER_B}_${UPPER_B} \
        ${DATASET_EPI}.hg19_cpg_regions_${GENOMIC_INTERVAL}bp
done 
} < chr_split_hg19.tsv

# Delete intermediary tables
{ read
while read CHR LOWER_B UPPER_B ; do 
    echo "Chromosome is " ${CHR} ", lower:" ${LOWER_B} ", and upper:" ${UPPER_B}
    bq rm -f -t ${DATASET_EPI}.hg19_cpg_regions_${CHR}_${LOWER_B}_${UPPER_B} 
done 
} < chr_split_hg19.tsv


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

SAMPLE="hg19_cpg_regions_250bp_clean"

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

# Run the jobs. Only 100 at a time (BQ limit)
dsub \
    --project $PROJECT_ID \
    --zones $ZONE_ID \
    --image ${DOCKER_GCP} \
    --logging $LOG \
    --env DATASET_EPI="${DATASET_EPI}" \
    --env EPI_REGION="${EPI_REGION}" \
    --env DATASET="${DATASET_EPI}" \
    --script ${SCRIPTS}/annotation.sh \
    --tasks chr_split_epi.tsv \
    --wait

# Delete previous files in case
for EPI_SIGNAL in "dnase" "encode_ChiP_V2" "tf_motifs" ; do
    echo "epi signal: " ${EPI_SIGNAL}
    bq rm -f -t ${DATASET_EPI}.hg19_cpg_regions_250bp_${EPI_SIGNAL}_all
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



# Gather all scores in a single array per region. This creates as many tables as there
# are epigenetic signals for annotation.
SAMPLE="hg19_cpg_regions_250bp"
for EPI_SIGNAL in "dnase" "encode_ChiP_V2" "tf_motifs" ; do
    echo "Processing the signal " ${EPI_SIGNAL}
    bq query \
        --use_legacy_sql=false \
        --destination_table ${DATASET_EPI}.${SAMPLE}_${EPI_SIGNAL} \
        --replace=true \
        "
        WITH 
            EPI_AGG AS ( -- we group the DNASe scores together
                SELECT 
                    * EXCEPT(score),
                    ARRAY_AGG(STRUCT(score)) AS epi
                FROM ${DATASET_EPI}.${SAMPLE}_${EPI_SIGNAL}_all
                GROUP BY 
                    chr,
                    region_inf,
                    region_sup,
                    annotate_ref,
                    region_nb_cpg
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