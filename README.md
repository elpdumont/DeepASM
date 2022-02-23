# DeepASM

Last update: XXX XX, 2022.


## Overview

DeepASM evaluates the imbalance in methylation between alleles in genomic sequences, without separating the reads per allele using a single nucleotide polymorphism (SNP) in the region.

It is designed to run on Google [Cloud Platform](https://cloud.google.com/) (GCP).

The first step is to divide the reference genome in genomic sequences of equal length (e.g. 250 bp, 500 bp, 1000 bp).

