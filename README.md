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

|chr  | region_inf | region_sup | region_nb_cpg | dnase | encode_ChiP_V2 | tf_motifs |
|-----|------------|------------|---------------|-------|----------------|-----------|
| 1   |   87804001 |   87805000 |             6 |     3 |              1 |        39 |
| 1   |  208718001 |  208719000 |             8 |     1 |              1 |        33 |
| 1   |   38460001 |   38461000 |            15 |     2 |              1 |        25 |
| 1   |  152634001 |  152635000 |             7 |     0 |              0 |        25 |


## Prepare samples by integrating read fm, cpg fm, and cpg cov into the reference genome's genomic regions

Execute the bash scripts in `sample_preparation/sample_prep_main.sh`.

The first step labels the data in the context files with their associated genomic regions defined above.


### Annotate genomic regions with ASM information

Execute the bash scripts in `asm_annotation/asm_main.sh`. This pipeline was adapted from our published pipeline [CloudASM](https://academic.oup.com/bioinformatics/article/36/11/3558/5771329?login=false).







