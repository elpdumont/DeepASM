
#--------------------------------------------------------------------------
# Variables
#--------------------------------------------------------------------------

# Where scripts are located
SCRIPTS="/Users/emmanuel/GITHUB_REPOS/DeepASM/prediction"

# BQ dataset where the epigenetic windows are defined
DATASET_EPI="hg19"

# Size of genomic regions:
GENOMIC_INTERVAL="250" # must be the same that in hg19_preparation.sh

# BQ dataset where the output of CloudASM is located
DATASET_PRED="deepasm_june2020"

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
# Create CpG regions to be evaluated by DeepASM
#--------------------------------------------------------------------------

# Create genomic regions used to split jobs per chromosome 
INTERVAL="82000000"

# Prepare TSV file per chromosome (used for many jobs)
echo -e "SAMPLE\tCHR\tLOWER_B\tUPPER_B" > chr_split.tsv

# Create the windows with $INTERVALS bp in it
while read SAMPLE ; do 
    for CHR in `seq 1 22` X Y ; do
        echo "Processing chromosome" ${CHR} " for sample" ${SAMPLE}
        NUCLEOTIDES_IN_CHR=$(awk -v CHR="${CHR}" -F"\t" '{ if ($1 == CHR) print $2}' chr_length.txt)
        INF="1"
        SUP=$(( $NUCLEOTIDES_IN_CHR<$INTERVAL ? $NUCLEOTIDES_IN_CHR: $INTERVAL ))
        echo -e "${SAMPLE}\t${CHR}\t$INF\t$SUP" >> chr_split.tsv # for jobs
        while [ $NUCLEOTIDES_IN_CHR -gt $SUP ] ; do
            INCREMENT=$(( $NUCLEOTIDES_IN_CHR-$SUP<$INTERVAL ? $NUCLEOTIDES_IN_CHR-$SUP: $INTERVAL ))
            INF=$(( ${SUP} + 1 ))
            SUP=$(( ${SUP} + $INCREMENT ))
            echo -e "${SAMPLE}\t${CHR}\t$INF\t$SUP" >> chr_split.tsv

        done
    done
done < sample_id.txt

# Looping over the samples. This step requires manual intervention.
dsub \
--project $PROJECT_ID \
--zones $ZONE_ID \
--image ${DOCKER_GCP} \
--logging $LOG \
--env DATASET_PRED="${DATASET_PRED}" \
--env DATASET_EPI="${DATASET_EPI}" \
--env DATASET_CONTEXT="${DATASET_CONTEXT}" \
--env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
--script ${SCRIPTS}/cpg_regions.sh \
--tasks chr_split.tsv \
--wait

# 1-99, 100-198, 199-297, 298-396, 397-496, 497-577

# Delete previous tables
while read SAMPLE ; do
    echo "Deleting the table for sample " ${SAMPLE}
    bq rm -f -t ${DATASET_PRED}.${SAMPLE}_cpg_regions_${GENOMIC_INTERVAL}bp
done < sample_id.txt

# Append to a new table for each sample
{ read
while IFS=$'\t' read SAMPLE CHR LOWER_B UPPER_B ; do 
    echo "Sample is:" ${SAMPLE} ", Chromosome is " ${CHR} ", lower:" ${LOWER_B} ", and upper:" ${UPPER_B}
    bq cp --append_table \
        ${DATASET_PRED}.${SAMPLE}_cpg_regions_${CHR}_${LOWER_B}_${UPPER_B} \
        ${DATASET_PRED}.${SAMPLE}_cpg_regions_${GENOMIC_INTERVAL}bp
done 
} < chr_split.tsv



# Erase intermediary files.
{ read
while read SAMPLE CHR LOWER_B UPPER_B ; do 
    echo "Sample is:" ${SAMPLE} ", Chromosome is " ${CHR} ", lower:" ${LOWER_B} ", and upper:" ${UPPER_B}
    bq rm -f -t ${DATASET_PRED}.${SAMPLE}_cpg_regions_${CHR}_${LOWER_B}_${UPPER_B}
done 
} < chr_split.tsv



#--------------------------------------------------------------------------
# Calculate fractional methylation of each CpG in each region.
#--------------------------------------------------------------------------

# We request that each CpG is covered at least 10x
# We request that there are at least 3 CpGs that are covered at least 10x in each region

dsub \
    --project $PROJECT_ID \
    --zones $ZONE_ID \
    --image ${DOCKER_GCP} \
    --logging $LOG \
    --env DATASET_PRED="${DATASET_PRED}" \
    --env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
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
    --env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
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
    --env GENOMIC_INTERVAL="${GENOMIC_INTERVAL}" \
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
        FROM ${DATASET_PRED}.${SAMPLE}_regions_${GENOMIC_INTERVAL}bp
        "

done < sample_id.txt


#--------------------------------------------------------------------------
# Combine all samples in a single table to be given to the model
#--------------------------------------------------------------------------

bq rm -f -t ${DATASET_PRED}.data_for_model

while read SAMPLE ; do
    echo "Processing sample " ${SAMPLE}
    bq cp --append_table \
            ${DATASET_PRED}.${SAMPLE}_regions_annotated \
            ${DATASET_PRED}.data_for_model
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

