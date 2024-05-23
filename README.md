# DeepASM

Last update: May 23, 2024.

Pre-print: 

## PRe-requesites

1. An account on [Google Cloud Platform](https://cloud.google.com/) (GCP).

2. For each whole genome, having processed the CpG context files, flag CpGs with the neareast SNP when possible, and flag each sequencing read with REF or ALT (parental chromosome) when possible, using the pipeline [CloudASM](https://github.com/TyckoLab/CloudASM) and store them on BigQuery in a dataset. 

3. Configure the `GCP` variables in `config.yaml` accordingly.


## Preparation of the reference genome with annotations (hg19_preparation folder)

Execute `main.sh`. 