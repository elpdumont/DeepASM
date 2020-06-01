#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET}.${TABLE}_${EPI_SIGNAL}_${CHR} \
    --replace=true \
    "
    WITH 
    ASM AS (
        SELECT *
        FROM ${DATASET}.${TABLE}
        WHERE chr = '${CHR}'
    ),
    EPI_DATA AS (
        SELECT *
        FROM ${DATASET_EPI}.${EPI_SIGNAL}
        WHERE signal_chr = '${CHR}'
    ),
    COMBINED AS (
        SELECT * FROM ASM 
        LEFT JOIN EPI_DATA
        ON
            (signal_end <= enrich_ref + ${EPI_REGION}) AND
            (signal_start >= enrich_ref - ${EPI_REGION}) AND
            chr = signal_chr
        )
    SELECT * EXCEPT(signal_start, signal_end, signal_chr) FROM COMBINED
    "
