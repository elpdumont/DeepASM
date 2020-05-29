
#--------------------------------------------------------------------------
# Variables
#--------------------------------------------------------------------------

# Where scripts are located
SCRIPTS="/Users/emmanuel/GITHUB_REPOS/DeepASM"

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
# Datasets for jobs
#--------------------------------------------------------------------------



#--------------------------------------------------------------------------
# Split the reference genome into 500bp windows.
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
                hg19.hg19_windows \
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
  --script ${SCRIPTS}/reads_cpg_regions.sh \
  --tasks chr_split.tsv 1-100 \
  --wait



# Erase previous file if there is one.

bq rm -f -t ${DATASET_PRED}.${SAMPLE}_reads_CpG_regions

{ read
while read CHR LOWER_B UPPER_B ; do 
    echo "Chromosome is " ${CHR}
    bq cp --append_table \
        ${DATASET_PRED}.${SAMPLE}_reads_CpG_regions_${CHR}_${LOWER_B}_${UPPER_B} \
        ${DATASET_PRED}.${SAMPLE}reads_CpG_regions
    bq rm -f -t ${DATASET_PRED}.${SAMPLE}_reads_CpG_regions_${CHR}_${LOWER_B}_${UPPER_B}
done 
} < chr_split.tsv

bq rm -f -t ${DATASET_PRED}.${SAMPLE}_reads_CpG_regions_*


# Need to calculate fractional methylation of reads, 
#fractional methylation of CpGs, coverage of each CpG, all CpG positions. 
# Then, need to aggregate with other signals.

# We first create a table by intersecting context files 
# and the reads overlapping CpG regions

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_cpg_reads \
    --replace=true \
    "
    WITH 
        CONTEXT AS (
            SELECT * 
            FROM ${DATASET_CONTEXT}.${SAMPLE}_context_filtered
        ),
        READS_CPG_REGIONS AS (
            SELECT * EXCEPT(read_id, chr, read_start, read_end), read_id AS read_id_identified 
            FROM ${DATASET_PRED}.${SAMPLE}_reads_CpG_regions
        ),
        COMBINED AS (
        SELECT
            chr, 
            pos,
            meth,
            cov,
            read_id,
            region_inf,
            region_sup
        FROM CONTEXT
        INNER JOIN READS_CPG_REGIONS 
        ON read_id = read_id_identified
        )
        -- We keep the rows where the CpG is within the region.
        SELECT * FROM COMBINED WHERE pos >= region_inf AND pos <= region_sup
    "


# We create a table of CpG information per CpG region
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
            FROM ${DATASET_PRED}.${SAMPLE}_cpg_reads
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

# We create a table of read information per CpG region
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_read_fm \
    --replace=true \
    "
    WITH 
        DATASETS_JOINED AS (
            SELECT * 
            FROM ${DATASET_PRED}.${SAMPLE}_cpg_reads
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

# We now join the 2 informations (CpG and read)

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_PRED}.${SAMPLE}_cpg_regions \
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

# We check which regions were evaluated by CloudASM
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

##### Number of CpGs evaluated by both methods 
# Number of CpGs evaluated by DeepAS.
# 22,424,388 --> impossible there are only 28M in the genome (80%)
# 2,159,485 CpGs were evaluated by CloudASM (7.5% of all CpGs)


###################################################
# Enrich the CpG windows with additional epigenetic signals.
##################################################

