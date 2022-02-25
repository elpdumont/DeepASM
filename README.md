# DeepASM

Last update: XXX XX, 2022.

## PRe-requesites

you need to have Python, [dsub](https://github.com/DataBiosphere/dsub), and Docker installed on your machine. Dsub is used to execute jobs in parallele using Docker images. You also need an account on [Google Cloud Platform](https://cloud.google.com/) (GCP).

## Overview

DeepASM evaluates the imbalance in methylation between alleles in genomic sequences, without separating the reads per allele using a single nucleotide polymorphism (SNP) in the region.


## Preparation of the reference genome with annotations (hg19_preparation folder)

Execute the bash scripts in `hg19_preparation/hg19_preparation.sh`.

The first step is to divide the reference genome in genomic sequences of equal length (e.g. 250 bp, 500 bp, 1000 bp) and annotate the regions for epigenetic marks (DNase, ChiP-Seq, TF binding sites). 

The final output is a BigQuery table called `hg19_cpg_regions_XXXbp_clean_annotated` where XXX is the length of the genomic region selected for analysis.

| Row |	chr	| region_inf | region_sup | region_nb_cpg | dnase | encode_ChiP_V2 | tf_motifs |
|-----|-----|------------|------------|---------------|-------|----------------|-----------|
| 1|1|25856001|25857000|31|1|0|1|
|2|1|156194001|156195000|27|0|0|1|
|3|1|1470001|1471000|64|1|0|1|


## Prepare samples by integrating read fm, cpg fm, and cpg cov into the reference genome's genomic regions

Execute the bash scripts in `sample_preparation/sample_prep_main.sh`.

The first step labels the data in the context files with their associated genomic regions defined above.


### Annotate genomic regions with ASM information

Execute the bash scripts in `asm_annotation/asm_main.sh`. This pipeline was adapted from our published pipeline [CloudASM](https://academic.oup.com/bioinformatics/article/36/11/3558/5771329?login=false).







