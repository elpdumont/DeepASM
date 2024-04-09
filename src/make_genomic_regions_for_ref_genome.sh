#!/bin/bash

# Prepare file
echo -e "chr\tregion_inf\tregion_sup\tannotate_ref\tclustering_index" > chr_regions.txt
CLUSTERING_INDEX=0

# Create the windows with $GENOMIC_LENGTH base pairs in it
for CHR in $(seq 1 22) ; do
    echo "Processing chromosome ${CHR}"
    NUCLEOTIDES_IN_CHR=$(awk -v CHR="${CHR}" -F"\t" '{ if ($1 == CHR) print $2}' chr_length.txt)
    CUMULATIVE_NB_NUCLEOTIDES=$(awk -v CHR="${CHR}" -F"\t" '{ if ($1 == CHR) print $3}' chr_length.txt)
    echo "Cumulative number of nucleotides: ${CUMULATIVE_NB_NUCLEOTIDES}"

    INF=1
    SUP=${GENOMIC_LENGTH}
    ANNOTATE_REF=$(( (INF + SUP) / 2 ))

    # Append to chr_regions.txt
    echo -e "${CHR}\t${INF}\t${SUP}\t${ANNOTATE_REF}\t${CLUSTERING_INDEX}" >> chr_regions.txt

    # Continue creating windows as long as NUCLEOTIDES_IN_CHR is greater than SUP
    while [[ ${NUCLEOTIDES_IN_CHR} -gt ${SUP} ]] ; do
        # Determine the increment
        if [[ $((NUCLEOTIDES_IN_CHR - SUP)) -lt ${GENOMIC_LENGTH} ]]; then
            INCREMENT=$((NUCLEOTIDES_IN_CHR - SUP))
        else
            INCREMENT=${GENOMIC_LENGTH}
        fi

        # Update INF and SUP for the next interval
        INF=$((SUP + 1))
        SUP=$((SUP + INCREMENT))
        ANNOTATE_REF=$(( (INF + SUP) / 2 ))
        CLUSTERING_INDEX=$(((CUMULATIVE_NB_NUCLEOTIDES + ANNOTATE_REF) / NB_NUCLEOTIDES_PER_CLUSTER))

        # Append to chr_regions.txt
        echo -e "${CHR}\t${INF}\t${SUP}\t${ANNOTATE_REF}\t${CLUSTERING_INDEX}" >> chr_regions.txt
    done
done



# X\t155270560\n\
# Y\t59373566"

