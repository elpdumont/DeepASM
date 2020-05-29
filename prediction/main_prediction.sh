
#--------------------------------------------------------------------------
# Variables
#--------------------------------------------------------------------------

# Where scripts are located
SCRIPTS="/Users/emmanuel/GITHUB_REPOS/DeepASM/prediction"

# BQ dataset where the output of CloudASM is located
DATASET_PRED="deepasm_prediction_gm12878"

# Name of the sample for which ASM prediction is sought
SAMPLE="gm12878"

# BQ dataset where the sample's context files are located (naming defined by CloudASM)
DATASET_CONTEXT="cloudasm_gm12878"

# Cloud Storage location of the logs
LOG="gs://cloudasm-encode/logging/deepasm"

# Docker file required to run the scripts
DOCKER_GCP="google/cloud-sdk:255.0.0"

# GCP global variables
PROJECT_ID="hackensack-tyco"
REGION_ID="us-central1"
ZONE_ID="us-central1-b"


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

# Create genomic regions used to split jobs per chromosome and per interval
INTERVAL="10000000"

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

# Be careful to split the jobs in max 100 concomittant jobs (BQ limit) or update the limit in GCP
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env SAMPLE="${SAMPLE}" \
  --env DATASET_PRED="${DATASET_PRED}" \
  --env DATASET_CONTEXT="${DATASET_CONTEXT}" \
  --script ${SCRIPTS}/cpg_regions.sh \
  --tasks chr_split.tsv 1-99 \
  --wait

# Erase previous file if there is one.
bq rm -f -t ${DATASET_PRED}.${SAMPLE}_cpg_regions

# Append all files into a single file.
{ read
while read CHR LOWER_B UPPER_B ; do 
    echo "Chromosome is " ${CHR} ". Boundaries are " ${LOWER_B} "and " ${UPPER_B}
    # bq cp --append_table \
    #     ${DATASET_PRED}.${SAMPLE}_cpg_regions_${CHR}_${LOWER_B}_${UPPER_B} \
    #     ${DATASET_PRED}.${SAMPLE}_cpg_regions
    bq rm -f -t ${DATASET_PRED}.${SAMPLE}_cpg_regions_${CHR}_${LOWER_B}_${UPPER_B}
done 
} < chr_split.tsv



#--------------------------------------------------------------------------
# Calculate fractional methylation of each CpG in each region.
#--------------------------------------------------------------------------

# We request that each CpG is covered at least 10x
# We request that there are at least 3 CpGs in each region

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_cpg_fm \
    --replace=true \
    "
    WITH 
        DATASETS_JOINED AS (
            SELECT * 
            FROM ${DATASET_PRED}.${SAMPLE}_cpg_regions
        ),
        CPG_FRAC_METHYL AS (
        -- Creates a list of all CpG sites with their respective fractional
        -- methylation and their CpG region
            SELECT 
                chr, 
                pos, 
                SUM(cov) AS cov,
                ROUND(SUM(meth)/SUM(cov),3) AS fm,
                region_inf,
                region_sup
            FROM DATASETS_JOINED
            GROUP BY chr, pos, region_inf, region_sup
        ),
        GROUP_CPG_INFO_BY_REGION AS (
        SELECT
            region_inf,
            region_sup,
            chr,
            COUNT(*) AS nb_cpg_found,
            ARRAY_AGG(
                STRUCT(fm, cov, pos)
                ) AS cpg
        FROM CPG_FRAC_METHYL
        WHERE cov >= 10
        GROUP BY region_inf, region_sup, chr
        )
        SELECT * FROM GROUP_CPG_INFO_BY_REGION
        WHERE nb_cpg_found >= 3
        "

#--------------------------------------------------------------------------
# Calculate fractional methylation of each read in each region.
#--------------------------------------------------------------------------

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_read_fm \
    --replace=true \
    "
    WITH 
        DATASETS_JOINED AS (
            SELECT * 
            FROM ${DATASET_PRED}.${SAMPLE}_cpg_regions
        ),
        READ_FRAC_METHYL AS (
            SELECT 
                ROUND(SUM(meth)/SUM(cov),3) AS fm,
                chr,
                region_inf,
                region_sup
            FROM DATASETS_JOINED
            GROUP BY read_id, chr, region_inf, region_sup
        )
        SELECT
            chr,
            region_inf,
            region_sup,
            ARRAY_AGG (STRUCT (fm)) AS read
        FROM READ_FRAC_METHYL
        GROUP BY chr, region_inf, region_sup
    "

#--------------------------------------------------------------------------
# Create a table of regions with fractional methylation of CpGs and reads
#--------------------------------------------------------------------------

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_regions \
    --replace=true \
    "
    WITH 
        CPG_INFO AS (
            SELECT * 
            FROM ${DATASET_PRED}.${SAMPLE}_cpg_fm
        ),
        READ_INFO AS (
            SELECT 
                chr AS chr_read,
                region_inf AS region_inf_read,
                region_sup AS region_sup_read,
                read
            FROM ${DATASET_PRED}.${SAMPLE}_read_fm
        )
        SELECT * EXCEPT(
                    chr_read,
                    region_inf_read, 
                    region_sup_read
                    ) 
        FROM CPG_INFO
        INNER JOIN READ_INFO
        ON 
            region_inf =  region_inf_read 
            AND region_sup = region_sup_read
            AND chr_read = chr
    "

#--------------------------------------------------------------------------
# Compare the number of CpGs and regions evaluated by DeepASM and CloudASM
#--------------------------------------------------------------------------

# Number of CpGs evaluated by CloudASM. Should return 2,159,485 CpGs
bq query \
    --use_legacy_sql=false \
    "
    WITH 
        UNIQUE_REGIONS AS ( 
            SELECT nb_cpg 
            FROM cloudasm_gm12878.gm12878_asm_snp
            GROUP BY chr, nb_ref_reads, nb_alt_reads, 
                asm_region_effect, wilcoxon_corr_pvalue, 
                nb_cpg, nb_sig_cpg, nb_pos_sig_cpg, 
                nb_neg_sig_cpg, nb_consec_pos_sig_asm, 
                nb_consec_neg_sig_asm
        ) 
    SELECT SUM(nb_cpg) FROM UNIQUE_REGIONS
    "

-- 39,189,346 Too many!!

# The following command should return 22,370,184 CpG (10x coverage, 3 CpG per region)
bq query \
    --use_legacy_sql=false \
    "
    SELECT SUM(nb_cpg_found) 
    FROM ${DATASET_PRED}.${SAMPLE}_regions
    "


# There are 461,346 regions that were evalulated by CloudASM
# DeepASm will evaluate 4,644,819 regions
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_cpg_regions_cloudasm \
    --replace=true \
    "
    WITH 
        CPG_REGIONS AS (
            SELECT *
            FROM ${DATASET_PRED}.${SAMPLE}_cpg_regions
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
            ANY_VALUE(nb_cpg_found),
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


###################################################
# Enrich the CpG windows with additional epigenetic signals.
##################################################

