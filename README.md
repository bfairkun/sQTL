# Snakemake workflow: Chimp_sQTL_pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.3.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/sQTL_pipeline.svg?branch=master)](https://travis-ci.org/snakemake-workflows/sQTL_pipeline)


## Authors

* Benjamin Fair (@bfairkun)

## Overview
This pipeline includes read mapping (STAR), preparation of a phenotype table of splicing traits (leafcutter), and sQTL calling (MatrixEQTL calculate nominal associations, and run permutations, saving the best P-value for each intron for each permutation). You will have to make your own script to actually calculate the intron-level, cluster-level, or gene-level (one of those might make more sense than the others depending on your downstream question) P-value from the permutation results. I suppose you are interested in gene-level P-values, could always calculate intron level P-values for the permutation test and then simply ask if the minimum, but then I think it makes sense to do multiple test correction on all introns which might be overly conservative.

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, clone the [latest release](https://github.com/bfairkun/sQTL_pipeline).
If you intend to modify and further develop this workflow, fork this repository. Please consider providing any generally applicable modifications via a pull request.

#### Install dependencies with conda:

`conda env create --file environment.yaml`

Other dependencies that I could not include on conda include the scripts for leafcutter . I have my own fork with small modifications that are required for this pipeline to work:

[leafcutter](https://github.com/bfairkun/leafcutter): modified script to allow nonconventional chromosome names (eg: 2A)

Clone my forks linked above, and add the necessary scripts to $PATH by appending the following to .bashrc:
```
export PATH=$PATH:PathToLeacutterClonedRepo/scripts
export PATH=$PATH:PathToLeacutterClonedRepo/clustering
```

re-source the .bashrc:
`source ~/.bashrc`

Make sure tidyverse and MatrixEQTL are installed for R... I have been using RCC's R/3.4.3 (`module load R/3.4.3`), and installed these with `install.packages()` once in R.

activate the conda environment:
`conda activate my_Chimp_EQTL_env` 

and create rule-specic environments:
`snakemake --use-conda --create-envs-only`

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`. Configure cluster settings in `cluster-config.json`

### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --cluster --cluster-config cluster-config.json --cluster "sbatch --partition={cluster.partition} --job-name={cluster.name} --output=/dev/null --job-name={cluster.name} --nodes={cluster.n} --mem={cluster.mem}"

or by executing the included sbatch script to execute the snakemake process from a cluster

    sbatch snakemake.sbatch

See the [Snakemake documentation](https://snakemake.readthedocs.io) for further details.
