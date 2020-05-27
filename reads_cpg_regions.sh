#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table deepasm_prediction_gm12878.gm12878_reads_CpG_regions_${CHR}_${LOWER_B}_${UPPER_B} \
    --replace=true \
    "
    WITH 
        cpg_regions AS (
            SELECT * 
            FROM hg19.hg19_CpG_regions
            WHERE chr_region = '${CHR}'
                AND region_sup <= ${UPPER_B}
                AND region_inf >= ${LOWER_B}
        ),
        sam_sample AS (
            SELECT read_id, chr, read_start, read_end
            FROM cloudasm_gm12878.gm12878_recal_sam
            WHERE chr = '${CHR}'
                AND read_end <= ${UPPER_B}
                AND read_start >= ${LOWER_B}
        )
        SELECT * EXCEPT(chr_region) FROM cpg_regions
        INNER JOIN sam_sample
        ON (read_end > region_inf AND read_start < region_inf)
        OR (read_start < region_sup AND read_end > region_sup) 
        OR (read_start >= region_inf AND read_end <= region_sup)
    "