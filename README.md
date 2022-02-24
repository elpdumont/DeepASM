# DeepASM

Last update: XXX XX, 2022.

## PRe-requesites

you need to have Python, [dsub](https://github.com/DataBiosphere/dsub), and Docker installed on your machine. Dsub is used to execute jobs in parallele using Docker images. You also need an account on [Google Cloud Platform](https://cloud.google.com/) (GCP).

## Overview

DeepASM evaluates the imbalance in methylation between alleles in genomic sequences, without separating the reads per allele using a single nucleotide polymorphism (SNP) in the region.


## Preparation of the reference genome with annotations

The first step is to divide the reference genome in genomic sequences of equal length (e.g. 250 bp, 500 bp, 1000 bp) and annotate the regions for epigenetic marks (DNase, ChiP-Seq, TF binding sites). To do that, the scripts are in `hg19_preparation` where the code in `hg19_preparation.sh` need to be excuted sequentially.

The final output is a BigQuery table called `hg19_cpg_regions_XXXbp_clean_annotated` where XXX is the length of the genomic region selected for analysis.

| Row |	chr	| region_inf | region_sup | region_nb_cpg | dnase | encode_ChiP_V2 | tf_motifs |
|-----|-----|------------|------------|---------------|-------|----------------|-----------|
| 1|1|25856001|25857000|31|1|0|1|
|2|1|156194001|156195000|27|0|0|1|
|3|1|1470001|1471000|64|1|0|1|


## Calculate allele-specific methylation in pre-defined genomic regions for several ENCODE samples


