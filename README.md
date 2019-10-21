# Snakemake workflow: Chimp_sQTL_pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.3.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/sQTL_pipeline.svg?branch=master)](https://travis-ci.org/snakemake-workflows/sQTL_pipeline)


## Authors

* Benjamin Fair (@bfairkun)

## Overview
This pipeline includes read mapping (STAR), preparation of a phenotype table of splicing traits (leafcutter), and sQTL calling (MatrixEQTL calculate nominal associations, eigenMT to obtain intron level P-values).

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, clone the [latest release](https://github.com/bfairkun/sQTL_pipeline).
If you intend to modify and further develop this workflow, fork this repository. Please consider providing any generally applicable modifications via a pull request.

Install dependencies with conda:
`conda env create --file environment.yaml`

Other dependencies that I could not include on conda include the scripts for leafcutter . I have my own fork with small modifications that are required for this pipeline to work:

[leafcutter](https://github.com/bfairkun/leafcutter): modified script to allow nonconventional chromosome names (eg: 2A)

Add the necessary scripts to $PATH by appending the following to .bashrc:
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
