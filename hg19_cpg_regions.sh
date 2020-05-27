#!/bin/bash

bq query \
        --use_legacy_sql=false \
        --destination_table hg19.hg19_CpG_regions_${CHR} \
        --replace=true \
        "
        WITH 
        hg19_CpG AS (
            SELECT chr, inf 
            FROM hg19.hg19_CpG_pos
            WHERE chr = '${CHR}'
        ),
        CpG_windows AS (
            SELECT chr, inf AS CpG_pos, inf + ${CPG_DISTANCE} AS CpG_next
            FROM hg19_CpG
            ORDER by chr, CpG_pos
        ),
        CpG_regions AS (
        SELECT chr, CpG_pos, CpG_next,
            COUNTIF(newRange) OVER(ORDER BY chr, CpG_pos) AS region_id
        FROM (
            SELECT *, 
            CpG_pos >= MAX(CpG_next) OVER(ORDER BY CpG_pos ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) AS newRange
            FROM CpG_windows
            )
        ),
        combined AS (
            SELECT
                ANY_VALUE(chr) AS chr_region,
                min(CpG_pos)-1 AS region_inf, 
                max(CpG_pos)+1 AS region_sup,
                max(CpG_pos) - min(CpG_pos) + 2 AS region_length,
                COUNT(*) AS nb_CpGs
            FROM CpG_regions 
            GROUP BY region_id
        )
        SELECT * FROM combined WHERE nb_CpGs >= ${MIN_CPG}
        "