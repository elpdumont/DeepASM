# DeepASM

Last update: May 23, 2024.

Pre-print: [link on biorxiv](https://www.biorxiv.org/content/10.1101/2024.05.23.595536v1.full.pdf+html).

## Pre-requesites

1. An account on [Google Cloud Platform](https://cloud.google.com/) (GCP).

2. For each whole genome, you need to output the CpG context files, flag CpGs with the neareast SNP when possible, and flag each sequencing read with REF or ALT (parental chromosome) when possible, using the pipeline [CloudASM](https://github.com/TyckoLab/CloudASM) and store them on BigQuery in a dataset. 

3. Configure the `GCP` variables in `config.yaml` accordingly. 


## Script execution

All the steps are in `main.sh`, which can be executed. The main steps are the following:
1. Create a python-based and bash-based image on GCP's artifact repository using the docker folder in the repository.
2. define environmental variables in bash using `config.yaml`
3. Prepare the reference genome in genomic regions of 250 nucleotides.
4. Evaluate ASM using CloudASM in the samples (here, ENCODE).
5. Machine Learning. See our pre-print for a detailed description of the steps.