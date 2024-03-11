# DeepASM

Last update: XXX XX, 2022.

## PRe-requesites

You need an account on [Google Cloud Platform](https://cloud.google.com/) (GCP).

## Overview

DeepASM evaluates the imbalance in methylation between alleles in genomic sequences, without separating the reads per allele using a single nucleotide polymorphism (SNP) in the region.

## Prepare GCP containers

```
gcloud config set project hmh-em-deepasm
gcloud config set run/region us-east1
```

gcloud run jobs deploy process-json \
 --image us-east1-docker.pkg.dev/hmh-em-deepasm/docker-repo/python:latest \
 --args="python,/app/process_json.py" \
 --tasks 50 \
 --set-env-vars BUCKET_NAME="hmh_deepasm" \
 --set-env-vars INPUT_DATA_FOLDER_PATH="bq_tables/250bp_asm_labelled/" \
 --set-env-vars OUTPUT_DATA_FOLDER_PATH="ml_datasets/" \
 --set-env-vars GENOMIC_LENGTH="250" \
 --set-env-vars MIN_CPG_COV="20" \
 --set-env-vars KERNEL_FM_NB_VALUES="10" \
 --set-env-vars KERNEL_FM_BANDWIDTH="0.1" \
 --set-env-vars KERNEL_COV_NB_MAX="200" \
 --set-env-vars KERNEL_COV_NB_STEP="10" \
 --set-env-vars KERNEL_COV_BANDWIDTH="5" \
 --max-retries 0 \
 --cpu 4 \
 --memory 16Gi \
 --region us-east1 \
 --project=hmh-em-deepasm

gcloud run jobs execute process-json

## Preparation of the reference genome with annotations (hg19_preparation folder)

Execute the bash scripts in `hg19_preparation/hg19_preparation.sh`.

The first step is to divide the reference genome in genomic sequences of equal length (e.g. 250 bp, 500 bp, 1000 bp) and annotate the regions for epigenetic marks (DNase, ChiP-Seq, TF binding sites).

The final output is a BigQuery table called `hg19_cpg_regions_XXXbp_clean_annotated` where XXX is the length of the genomic region selected for analysis.

| chr | region_inf | region_sup | region_nb_cpg | dnase | encode_ChiP_V2 | tf_motifs |
| --- | ---------- | ---------- | ------------- | ----- | -------------- | --------- |
| 1   | 87804001   | 87805000   | 6             | 3     | 1              | 39        |
| 1   | 208718001  | 208719000  | 8             | 1     | 1              | 33        |
| 1   | 38460001   | 38461000   | 15            | 2     | 1              | 25        |
| 1   | 152634001  | 152635000  | 7             | 0     | 0              | 25        |

## Prepare samples by integrating read fm, cpg fm, and cpg cov into the reference genome's genomic regions

Execute the bash scripts in `sample_preparation/sample_prep_main.sh`.

The first step labels the data in the context files with their associated genomic regions defined above.

### Annotate genomic regions with ASM information

Execute the bash scripts in `asm_annotation/asm_main.sh`. This pipeline was adapted from our published pipeline [CloudASM](https://academic.oup.com/bioinformatics/article/36/11/3558/5771329?login=false).
