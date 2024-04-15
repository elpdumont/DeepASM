#!/bin/bash

# Environmental variables: SAMPLE_LIST, DATA_LIST
# Given by the job: BATCH_TASK_INDEX, BATCH_TASK_COUNT

# Initialize a counter for the mapping index
index=0

# Convert string to array
read -a sample_list <<< "${SAMPLE_LIST}"
read -a table_list <<< "${CLOUDASM_TABLES}"


# Loop through each sample
for sample in "${sample_list[@]}"; do
    # Loop through each data item
    for table in "${table_list[@]}"; do
        # Increment the index for each combination
        ((index++))
e
        # Check if the current index matches the BATCH_TASK_INDEX
        if [[ ${index} -eq ${BATCH_TASK_INDEX} ]]; then
            echo "Sample and Table for BATCH_TASK_INDEX ${BATCH_TASK_INDEX}:"
            echo "Sample: ${sample}"
            echo "Table: ${table}"
            break 2  # Exit both loops
        fi
    done
done

# Check if BATCH_TASK_INDEX was out of range
if [[ ${index} -lt ${BATCH_TASK_INDEX} ]]; then
    echo "BATCH_TASK_INDEX ${BATCH_TASK_INDEX} is out of range. Maximum index is ${index}."
fi

TEMP_TABLE="${SAMPLES_DATASET}.${sample}_${table}_temp"
TEMP_TABLE_2="${SAMPLES_DATASET}.${sample}_${table}_temp2"

# Create and populate the temporary table, adding the 'sample' column
bq query \
    --use_legacy_sql=false \
    --destination_table="${TEMP_TABLE}" \
    --replace=true \
    --clustering_fields=sample,chr \
    --range_partitioning=clustering_index,0,4000,1 \
    "
    SELECT
        '${sample}' AS sample,
        CAST(FLOOR((c.absolute_nucleotide_pos + p.pos) / ${NB_NUCLEOTIDES_PER_CLUSTER}) AS INT64) AS clustering_index,
        CAST(p.chr AS INT64) AS chr,
        p.* EXCEPT(chr)
        FROM ${PROJECT}.${CLOUDASM_DATASET}.${sample}_${dataset} p
        JOIN ${REFG_FOLDER}.chr_length c
        ON SAFE_CAST(p.chr AS INT64) = c.chr
        WHERE REGEXP_CONTAINS(p.chr, r'^\\d+$')
    "

# Inner join the table with the regions of the reference genome
bq query \
    --use_legacy_sql=false \
    --destination_table="${TEMP_TABLE_2}" \
    --replace=true \
    --clustering_fields=sample,chr \
    --range_partitioning=clustering_index,0,4000,1 \
    "
    SELECT c.region_inf, c.region_sup, c.region_nb_cpg, p.*
    FROM ${TEMP_TABLE} p
    INNER JOIN ${REFG_FOLDER}.regions_w_cpg_no_blacklist c
    ON p.chr = c.chr AND p.clustering_index = c.clustering_index AND pos >= region_inf AND pos < region_sup
    "


# Maximum number of retries
max_retries=5
# Initial delay in seconds (e.g., 1 second)
delay=1

for (( i=0; i<max_retries; i++ )); do
    # Attempt to append the table
    bq cp --append_table \
        "${TEMP_TABLE_2}" \
        "${PROJECT}:${SAMPLES_DATASET}.${table}" && break

    # If the append fails, print an error message and wait
    echo "Attempt $((i+1)) failed! Retrying in ${delay} seconds..."
    sleep "${delay}"

    # Exponentially increase the delay
    delay=$((delay * 2))
done

# Check if all attempts failed
if [[ ${i} -eq max_retries ]]; then
    echo "All attempts to append the table have failed!"
fi


# Append table
# bq cp --append_table \
#     "${TEMP_TABLE_2}" \
#     "${PROJECT}":"${SAMPLES_DATASET}"."${table}"

# Delete temp table
bq rm -f -t "${TEMP_TABLE}"
bq rm -f -t "${TEMP_TABLE_2}"