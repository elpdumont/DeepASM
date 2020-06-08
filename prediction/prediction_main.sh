
#--------------------------------------------------------------------------
# Variables
#--------------------------------------------------------------------------

# Where scripts are located
SCRIPTS="/Users/emmanuel/GITHUB_REPOS/DeepASM/prediction"

# Where the scripts for enrichment are located
ENRICH_SCRIPTS="/Users/emmanuel/GITHUB_REPOS/DeepASM/enrichment"

# Desired window for enrichment analysis
EPI_REGION="250"

# BQ dataset where the epigenetic windows are defined
DATASET_EPI="hg19"

# BQ dataset where the output of CloudASM is located
DATASET_PRED="deepasm_prediction"

# BQ dataset where the sample's context files are located (naming defined by CloudASM)
DATASET_CONTEXT="cloudasm_encode_2019"

# Cloud Storage location of the logs
LOG="gs://cloudasm-encode/logging/deepasm"

# Docker file required to run the scripts
DOCKER_GCP="google/cloud-sdk:255.0.0"

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
# Split the human genome into 500bp windows.
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

# Create genomic regions with the same interval:
INTERVAL="500"

# Prepare TSV file per chromosome (used for many jobs)
echo -e "chr\tregion_inf\tregion_sup" > chr_regions.txt

# Create the windows with $INTERVALS base pairs in it
for CHR in `seq 1 22` X Y ; do
    echo "Processing chromosome" ${CHR}
    NUCLEOTIDES_IN_CHR=$(awk -v CHR="${CHR}" -F"\t" '{ if ($1 == CHR) print $2}' chr_length.txt)
    INF="1"
    SUP=$(( $NUCLEOTIDES_IN_CHR<$INTERVAL ? $NUCLEOTIDES_IN_CHR: $INTERVAL ))
    echo -e "${CHR}\t$INF\t$SUP" >> chr_regions.txt # for jobs
    while [ $NUCLEOTIDES_IN_CHR -gt $SUP ] ; do
        INCREMENT=$(( $NUCLEOTIDES_IN_CHR-$SUP<$INTERVAL ? $NUCLEOTIDES_IN_CHR-$SUP: $INTERVAL ))
        INF=$(( ${SUP} + 1 ))
        SUP=$(( ${SUP} + $INCREMENT ))
        echo -e "${CHR}\t$INF\t$SUP" >> chr_regions.txt

    done
done

# Upload the file to a bucket
gsutil cp chr_regions.txt gs://cloudasm-encode/chr_regions.txt

# Import the file in BigQuery
bq --location=US load \
               --replace=true \
               --source_format=CSV \
               --field_delimiter "\t" \
               --skip_leading_rows 1 \
                ${DATASET_PRED}.regions \
               gs://cloudasm-encode/chr_regions.txt \
               chr_region:STRING,region_inf:INT64,region_sup:INT64


#--------------------------------------------------------------------------
# Create CpG regions to be evaluated by DeepASM
#--------------------------------------------------------------------------

# Create genomic regions used to split jobs per chromosome 
# and per interval. We picked the interval to have less than
# 100 queries running at the same time (BQ limit)
INTERVAL="36000000"

# Prepare TSV file per chromosome (used for many jobs)
echo -e "CHR\tLOWER_B\tUPPER_B" > chr_split.tsv

# Create the windows with $INTERVALS bp in it
for CHR in `seq 1 22` X Y ; do
    echo "Processing chromosome" ${CHR}
    NUCLEOTIDES_IN_CHR=$(awk -v CHR="${CHR}" -F"\t" '{ if ($1 == CHR) print $2}' chr_length.txt)
    INF="1"
    SUP=$(( $NUCLEOTIDES_IN_CHR<$INTERVAL ? $NUCLEOTIDES_IN_CHR: $INTERVAL ))
    echo -e "${CHR}\t$INF\t$SUP" >> chr_split.tsv # for jobs
    while [ $NUCLEOTIDES_IN_CHR -gt $SUP ] ; do
        INCREMENT=$(( $NUCLEOTIDES_IN_CHR-$SUP<$INTERVAL ? $NUCLEOTIDES_IN_CHR-$SUP: $INTERVAL ))
        INF=$(( ${SUP} + 1 ))
        SUP=$(( ${SUP} + $INCREMENT ))
        echo -e "${CHR}\t$INF\t$SUP" >> chr_split.tsv

    done
done

# Looping over the samples. This step requires manual intervention.
while read SAMPLE ; do
    dsub \
    --project $PROJECT_ID \
    --zones $ZONE_ID \
    --image ${DOCKER_GCP} \
    --logging $LOG \
    --env SAMPLE="${SAMPLE}" \
    --env DATASET_PRED="${DATASET_PRED}" \
    --env DATASET_CONTEXT="${DATASET_CONTEXT}" \
    --script ${SCRIPTS}/cpg_regions.sh \
    --tasks chr_split.tsv \
    --wait
done < sample_id.txt


# Append all files into a single file.
while read SAMPLE ; do
    echo "Processing sample " ${SAMPLE}
    bq rm -f -t ${DATASET_PRED}.${SAMPLE}_cpg_regions
    { read
    while read CHR LOWER_B UPPER_B ; do 
        echo "Chromosome is " ${CHR} "--" ${LOWER_B} "---" ${UPPER_B}
        bq cp --append_table \
            ${DATASET_PRED}.${SAMPLE}_cpg_regions_${CHR}_${LOWER_B}_${UPPER_B} \
            ${DATASET_PRED}.${SAMPLE}_cpg_regions
    done 
    } < chr_split.tsv
done < sample_id.txt

# Erase intermediary files.
while read SAMPLE ; do
    echo "Processing sample " ${SAMPLE}
    { read
    while read CHR LOWER_B UPPER_B ; do 
        echo "Chromosome is " ${CHR} "--" ${LOWER_B} "---" ${UPPER_B}
        bq rm -f -t ${DATASET_PRED}.${SAMPLE}_cpg_regions_${CHR}_${LOWER_B}_${UPPER_B}
    done 
    } < chr_split.tsv
done < sample_id.txt



#--------------------------------------------------------------------------
# Calculate fractional methylation of each CpG in each region.
#--------------------------------------------------------------------------

# We request that each CpG is covered at least 10x
# We request that there are at least 3 CpGs in each region

dsub \
    --project $PROJECT_ID \
    --zones $ZONE_ID \
    --image ${DOCKER_GCP} \
    --logging $LOG \
    --env DATASET_PRED="${DATASET_PRED}" \
    --script ${SCRIPTS}/cpg_fm.sh \
    --tasks all_samples.tsv \
    --wait


#--------------------------------------------------------------------------
# Calculate fractional methylation of each read in each region.
#--------------------------------------------------------------------------

dsub \
    --project $PROJECT_ID \
    --zones $ZONE_ID \
    --image ${DOCKER_GCP} \
    --logging $LOG \
    --env DATASET_PRED="${DATASET_PRED}" \
    --script ${SCRIPTS}/read_fm.sh \
    --tasks all_samples.tsv \
    --wait

#--------------------------------------------------------------------------
# Create a table of regions with fractional methylation of CpGs and reads
#--------------------------------------------------------------------------


dsub \
    --project $PROJECT_ID \
    --zones $ZONE_ID \
    --image ${DOCKER_GCP} \
    --logging $LOG \
    --env DATASET_PRED="${DATASET_PRED}" \
    --script ${SCRIPTS}/combine_read_cpg_fm.sh \
    --tasks all_samples.tsv \
    --wait



#--------------------------------------------------------------------------
# Compare the number of CpGs and regions evaluated by DeepASM and CloudASM
#--------------------------------------------------------------------------


while read SAMPLE ; do
    echo " "
    echo "********************"
    echo "********************"
    echo "Sample is " ${SAMPLE}
    echo " "
    echo "********************"
    echo "Number of CpGs evaluated by CloudASM"
    bq query \
        --use_legacy_sql=false \
        "
        WITH 
            UNIQUE_REGIONS AS ( 
                SELECT nb_cpg 
                FROM ${DATASET_CONTEXT}.${SAMPLE}_asm_snp
                GROUP BY chr, nb_ref_reads, nb_alt_reads, 
                    asm_region_effect, wilcoxon_corr_pvalue, 
                    nb_cpg, nb_sig_cpg, nb_pos_sig_cpg, 
                    nb_neg_sig_cpg, nb_consec_pos_sig_asm, 
                    nb_consec_neg_sig_asm
            ) 
        SELECT SUM(nb_cpg) FROM UNIQUE_REGIONS
        "

    echo " "
    echo "********************"
    echo "Number of CpGs evaluated by DeepASM"

    bq query \
        --use_legacy_sql=false \
        "
        SELECT SUM(nb_cpg_found) 
        FROM ${DATASET_PRED}.${SAMPLE}_regions
        "

done < sample_id.txt


#--------------------------------------------------------------------------
# Overlap the results with the CloudASM results
#--------------------------------------------------------------------------

while read SAMPLE ; do
    echo "********************"
    echo "Sample is " ${SAMPLE}
    echo "********************"
    bq query \
        --use_legacy_sql=false \
        --destination_table ${DATASET_PRED}.${SAMPLE}_regions_cloudasm \
        --replace=true \
        "
        WITH 
            CPG_REGIONS AS (
                SELECT *
                FROM ${DATASET_PRED}.${SAMPLE}_regions
            ),
            ASM_REGIONS AS (
                SELECT * EXCEPT(region_inf, region_sup, chr, region_length),
                    chr AS chr_asm, 
                    region_length AS region_length_asm,
                    region_inf AS region_inf_asm, 
                    region_sup AS region_sup_asm
                FROM deepasm_encode.asm_for_bq
                WHERE sample = '${SAMPLE}'
            ),
            CPG_ASM_REGIONS AS (
                SELECT * FROM CPG_REGIONS
                INNER JOIN ASM_REGIONS
                ON (region_sup_asm > region_inf AND region_inf_asm < region_inf)
                OR (region_inf_asm < region_sup AND region_sup_asm > region_sup) 
                OR (region_inf_asm > region_inf AND region_sup_asm < region_sup)
            )
            SELECT 
                chr,
                region_inf, 
                region_sup,
                ANY_VALUE(nb_cpg_found) AS nb_cpg_found,
                ANY_VALUE(cpg) AS cpg,
                ANY_VALUE(read) AS read,
                ARRAY_AGG (
                    STRUCT (
                        snp_id,
                        asm_snp,
                        nb_cpg,
                        region_length_asm,
                        region_inf_asm,
                        region_sup_asm
                        )
                    ) AS v
            FROM CPG_ASM_REGIONS
            GROUP BY region_inf, region_sup, chr
        "
done < sample_id.txt


#--------------------------------------------------------------------------
# Enrich the CpG windows with additional epigenetic signals.
#--------------------------------------------------------------------------

# Prepare TSV file per chromosome (used for many jobs)
echo -e "--env TABLE\t--env EPI_SIGNAL\t--env CHR" > all_chr.tsv
while read SAMPLE ; do
    for SIGNAL in "dnase" "encode_ChiP_V2" "tf_motifs" ; do
        for CHR in `seq 1 22` X Y ; do
        echo -e "${SAMPLE}_regions\t${SIGNAL}\t${CHR}" >> all_chr.tsv
        done
    done
done < sample_id.txt

# Run the jobs. Only 100 at a time (BQ limit)
dsub \
--project $PROJECT_ID \
--zones $ZONE_ID \
--image ${DOCKER_GCP} \
--logging $LOG \
--env DATASET_EPI="${DATASET_EPI}" \
--env EPI_REGION="${EPI_REGION}" \
--env DATASET="${DATASET_PRED}" \
--script ${ENRICH_SCRIPTS}/enrichment.sh \
--tasks all_chr.tsv 1-99 \
--wait



# Concatenate the files
while read SAMPLE ; do
    for EPI_SIGNAL in "dnase" "encode_ChiP_V2" "tf_motifs" ; do
        echo "Processing the signal " ${EPI_SIGNAL} "for sample " ${SAMPLE}
        bq rm -f -t ${DATASET_PRED}.${SAMPLE}_regions_${EPI_SIGNAL}_all

        for CHR in `seq 1 22` X Y ; do
            echo "Chromosome is:" ${CHR}
            bq cp --append_table \
                ${DATASET_PRED}.${SAMPLE}_regions_${EPI_SIGNAL}_${CHR} \
                ${DATASET_PRED}.${SAMPLE}_regions_${EPI_SIGNAL}_all
        done
    done
done < sample_id.txt

# Delete the intermediate files if concatenation was successful
while read SAMPLE ; do
    for EPI_SIGNAL in "dnase" "encode_ChiP_V2" "tf_motifs" ; do
        echo "Processing the signal " ${EPI_SIGNAL} "for sample " ${SAMPLE}
        for CHR in `seq 1 22` X Y ; do
            echo "Chromosome is:" ${CHR}
            bq rm -f -t ${DATASET_PRED}.${SAMPLE}_regions_${EPI_SIGNAL}_${CHR}
        done
    done
done < sample_id.txt



# Gather all scores in a single array per region. This creates as many tables as there
# are epigenetic signals for enrichment.
while read SAMPLE ; do
    for EPI_SIGNAL in "dnase" "encode_ChiP_V2" "tf_motifs" ; do
        echo "Processing the signal " ${EPI_SIGNAL} "for sample " ${SAMPLE}
        bq query \
            --use_legacy_sql=false \
            --destination_table ${DATASET_PRED}.${SAMPLE}_regions_${EPI_SIGNAL} \
            --replace=true \
            "
            WITH 
                EPI_AGG AS ( -- we group the DNASe scores together
                    SELECT 
                        * EXCEPT(score, cpg, read),
                        ARRAY_AGG(STRUCT(score)) AS epi
                    FROM ${DATASET_PRED}.${SAMPLE}_regions_${EPI_SIGNAL}_all
                    GROUP BY 
                        chr,
                        region_inf,
                        region_sup,
                        enrich_ref,
                        nb_cpg_found
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
done < sample_id.txt

# Delete the previous tables.
while read SAMPLE ; do
    for EPI_SIGNAL in "dnase" "encode_ChiP_V2" "tf_motifs" ; do
        echo "Processing the signal " ${EPI_SIGNAL} "for sample " ${SAMPLE}
        bq rm -f -t ${DATASET_PRED}.${SAMPLE}_regions_${EPI_SIGNAL}_all
    done
done < sample_id.txt

# Merge the original table (with the CpG and read arrays) with the epigenetic signals.
while read SAMPLE ; do
    echo "Processing the sample " ${SAMPLE}
    bq query \
        --use_legacy_sql=false \
        --destination_table ${DATASET_PRED}.${SAMPLE}_regions_enriched \
        --replace=true \
        "
        SELECT
            ${SAMPLE} AS sample,
            t1.chr AS chr,
            t1.region_inf AS region_inf,
            t1.region_sup AS region_sup,
            t1.nb_cpg_found AS nb_cpg_found,
            t1.dnase AS dnase,
            t2.encode_ChiP_V2 AS encode_ChiP_V2,
            t3.tf_motifs AS tf_motifs,
            t4.cpg AS cpg,
            t4.read AS read
        FROM ${DATASET_PRED}.${SAMPLE}_regions_dnase t1
        JOIN ${DATASET_PRED}.${SAMPLE}_regions_encode_ChiP_V2 t2 
        ON t1.chr = t2.chr AND 
            t1.region_inf = t2.region_inf AND 
            t1.region_sup = t2.region_sup
        JOIN ${DATASET_PRED}.${SAMPLE}_regions_tf_motifs t3 
        ON t1.chr = t3.chr AND 
        t1.region_inf = t3.region_inf AND 
        t1.region_sup = t3.region_sup
        JOIN ${DATASET_PRED}.${SAMPLE}_regions t4 
        ON t1.chr = t4.chr AND 
        t1.region_inf = t4.region_inf AND 
        t1.region_sup = t4.region_sup
        "
done < sample_id.txt

# Delete the individual enrichment tables
while read SAMPLE ; do
    for EPI_SIGNAL in "dnase" "encode_ChiP_V2" "tf_motifs" ; do
        bq rm -f -t ${DATASET_PRED}.${SAMPLE}_regions_${EPI_SIGNAL}
    done
done < sample_id.txt




#--------------------------------------------------------------------------
# Combine all samples in a single table to be given to the model
#--------------------------------------------------------------------------

bq rm -f -t ${DATASET_PRED}.data_for_model

while read SAMPLE ; do
    echo "Processing sample " ${SAMPLE}
    bq rm -f -t ${DATASET_PRED}.${SAMPLE}_cpg_regions
    { read
    while read CHR LOWER_B UPPER_B ; do 
        echo "Chromosome is " ${CHR} "--" ${LOWER_B} "---" ${UPPER_B}
        bq cp --append_table \
            ${DATASET_PRED}.${SAMPLE}_cpg_regions_${CHR}_${LOWER_B}_${UPPER_B} \
            ${DATASET_PRED}.${SAMPLE}_cpg_regions
    done 
    } < chr_split.tsv
done < sample_id.txt