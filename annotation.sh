#!/bin/bash

# Requires a reference "annotate_ref" in the center of the window.

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET}.${TABLE}_${EPI_SIGNAL}_${CHR}_${LOWER_B}_${UPPER_B} \
    --replace=true \
    "
    WITH 
    ASM AS (
        SELECT *
        FROM ${DATASET}.${TABLE}
        WHERE 
            chr = '${CHR}' AND 
            region_sup <= ${UPPER_B} AND
            region_inf >= ${LOWER_B}
    ),
    EPI_DATA AS (
        SELECT *
        FROM ${DATASET_EPI}.${EPI_SIGNAL}
        WHERE 
            signal_chr = '${CHR}' AND
            signal_end <= ${UPPER_B} AND
            signal_start >= ${LOWER_B}     
    ),
    COMBINED AS (
        SELECT * FROM ASM 
        LEFT JOIN EPI_DATA
        ON
            (signal_end <= annotate_ref + ${EPI_REGION}) AND
            (signal_start >= annotate_ref - ${EPI_REGION}) AND
            chr = signal_chr
        )
    SELECT * EXCEPT(signal_start, signal_end, signal_chr) FROM COMBINED
    "


