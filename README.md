# VCF Merge Sequential Pipeline

## Overview

This repository contains a pipeline designed to process VCF files by splitting them into manageable subsets, merging them, and then creating a single merged VCF file using `bcftools`.

## Contents

- `merge_sequential.smk`: The main Snakemake script for the pipeline.
- `config_dummy.yaml`: Configuration file dummy that stores parameters and settings for the pipeline.
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

## Logging

The script logs the start and end times of each rule to facilitate performance profiling and troubleshooting. Log files are stored in the directory specified under `log_subfolder` in the configuration file.

## Contribution

Feel free to fork the repository and submit pull requests for any enhancements or bug fixes. Contributions to improve the script or documentation are welcome.

## TODO

- Implement md5sum calculation for all files to verify data integrity.
- Remove intermediate files to save storage space.
- Ensure proper sequence of index and VCF file creation to maintain organization.
- Add error handling to manage potential issues gracefully, such as missing input files or unsuccessful command executions.
- Consider making the file extensions configurable to allow for more flexible input.
- Explore options for dynamic memory allocation based on the number of threads, possibly through a configuration setting or automatic calculation.
