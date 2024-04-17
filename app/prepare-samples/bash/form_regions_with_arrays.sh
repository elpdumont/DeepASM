#!/bin/bash

# Here we ensure that all CpGs coverages are within the range defined in the config file

# We also ensure to keep all regions where at least a min number of CpGs were found (defined in the config file)

# Note, it is possible to have two duplicate regions when more than one SNP was available.

bq query \
    --use_legacy_sql=false \
    --destination_table "${SAMPLES_DATASET}".regions_w_arrays \
    --replace=true \
    --clustering_fields=sample,chr \
    --range_partitioning=clustering_index,0,4000,1 \
    "
    WITH ALL_READ_ARRAY_TMP AS (
        SELECT
            chr,
            region_inf,
            region_sup,
            sample,
            region_nb_cpg,
            clustering_index,
            read_id AS id,
            ROUND(SUM(meth)/SUM(cov),3) AS fm,
            COUNT(*) AS nb_cpg,
            ARRAY_AGG (pos-region_inf + 1) AS cpg_pos,
            ARRAY_AGG (meth) AS cpg_meth
        FROM ${SAMPLES_DATASET}.all_cpgs_flagged_w_regions
        GROUP BY sample, clustering_index, chr, region_inf, region_sup, region_nb_cpg, read_id
        ),
    ALL_READ_ARRAY AS (
        SELECT
            chr,
            region_inf,
            region_sup,
            region_nb_cpg,
            sample,
            clustering_index,
            ARRAY_AGG(
                STRUCT(id, fm, nb_cpg, cpg_pos, cpg_meth)
            ) AS all_reads
        FROM ALL_READ_ARRAY_TMP
        GROUP BY sample, clustering_index, chr, region_inf, region_sup, region_nb_cpg
        ),
    CPG_ARRAY_TMP AS (
        SELECT
            chr,
            pos,
            region_inf,
            region_sup,
            sample,
            clustering_index,
            SUM(cov) AS cov,
            ROUND(SUM(meth)/SUM(cov),3) AS fm
        FROM ${SAMPLES_DATASET}.all_cpgs_flagged_w_regions
        GROUP BY sample, clustering_index, chr, region_inf, region_sup, pos
        HAVING cov >= ${MIN_CPG_COV} AND cov <= ${MAX_CPG_COV}
    ),
    CPG_ARRAY AS (
        SELECT
            chr,
            region_inf,
            region_sup,
            sample,
            clustering_index,
            COUNT(*) AS nb_cpg_found,
            ARRAY_AGG(
                STRUCT(fm, cov, pos) ORDER by pos
                ) AS cpg
        FROM CPG_ARRAY_TMP        
        GROUP BY sample, clustering_index, chr, region_inf, region_sup
        HAVING nb_cpg_found >= ${MIN_NB_CPG_FOR_ASM}
    ),
    READ_ARRAY_TMP AS (
        SELECT
            chr,
            region_inf,
            region_sup,
            sample,
            clustering_index,
            ROUND(SUM(meth)/SUM(cov),3) AS fm
        FROM ${SAMPLES_DATASET}.all_cpgs_flagged_w_regions
        GROUP BY sample, clustering_index, chr, region_inf, region_sup, read_id
    ),
    READ_ARRAY AS (
        SELECT
            chr,
            region_inf,
            region_sup,
            sample,
            clustering_index,
            ARRAY_AGG(
                STRUCT(fm) ORDER BY fm
                ) AS reads
        FROM READ_ARRAY_TMP
        GROUP BY sample, clustering_index, chr, region_inf, region_sup
    )
    SELECT
        p.chr,
        p.region_inf,
        p.region_sup,
        p.region_nb_cpg,
        q.nb_cpg_found,
        p.sample,
        p.clustering_index,
        q.cpg,
        c.reads,
        p.all_reads
    FROM ALL_READ_ARRAY p
    INNER JOIN CPG_ARRAY q
    ON p.chr = q.chr AND p.region_inf = q.region_inf AND p.region_sup = q.region_sup AND p.clustering_index = q.clustering_index AND p.sample = q.sample
    INNER JOIN READ_ARRAY c
    ON p.chr = c.chr AND p.region_inf = c.region_inf AND p.region_sup = c.region_sup AND p.clustering_index = c.clustering_index AND p.sample = c.sample
    "
