# VCF Merge Sequential Pipeline

## Overview

This repository contains a pipeline designed to process VCF files by splitting them into manageable subsets, merging them, and then creating a single merged VCF file using `bcftools`.

## Contents

- `merge_sequential.smk`: The main Snakemake script for the pipeline.
- `config.yaml`: Configuration file that stores parameters and settings for the pipeline.
- `run_merge_sequential.sh`: A scheduler shell script to run the Snakemake pipeline on a Slurm scheduler.

## Setup

### Configuration

Before running the pipeline, you need to set up the `config.yaml` file:

```yaml
vcf_folder: "/path/to/vcf/files"
output_folder: "/desired/output/folder"
vcfs_per_batch: 500
conda_env: "bcftools"
```

- `vcf_folder`: The path to the directory containing your `.vcf.gz` files.
- `output_folder`: The desired directory where the results will be saved.
- `vcfs_per_batch`: The number of VCFs processed in each batch.
- `conda_env`: The name of the Conda environment to be used, which must have `bcftools` installed.

### Conda Environment

Ensure that the specified Conda environment in the `config.yaml` file contains `bcftools`. 

## Running the Pipeline

### Directly with Snakemake

If you have Snakemake installed, you can run the pipeline directly:

```bash
snakemake -s merge_sequential.smk --use-conda
```

### With Slurm Scheduler

To run the pipeline on a system with the Slurm scheduler, use the provided shell script:

```bash
sbatch run_merge_sequential.sh
```

The script sets up necessary temporary directories, logs, and other parameters before invoking Snakemake. The logs for the Snakemake jobs are saved in the `slurm_logs` directory.

## Output

The pipeline will generate merged VCF files, with intermediate lists and logs saved in subdirectories defined in the `config.yaml`.