# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This repository contains two independent Snakemake pipelines for multi-sample VCF merging. Both are designed to run on HPC clusters (Slurm) with Conda environments.

## Running the Pipelines

### bcftools pipeline (direct VCF merge)

```bash
cd bcftools
snakemake -s merge_sequential.smk --use-conda
```

On Slurm:
```bash
cd bcftools
sbatch run_merge_sequential.sh
```

### GATK pipeline (GVCF joint genotyping)

```bash
cd gatk
snakemake -s reblock_genomicsdb_gatk_pipeline.smk --use-conda
```

Both pipelines require a `config.yaml` in their respective directories. Use the provided `config_dummy.yaml` (bcftools) as a template.

## Architecture

### bcftools pipeline (`bcftools/`)

Processes large sets of per-sample VCF files via a 4-stage sequential approach:

1. **normalize_vcf** — normalizes each input VCF with `bcftools norm` (left-align, split multiallelic)
2. **split_vcfs** — splits normalized VCFs into batches of `vcfs_per_batch` (workaround for bcftools' 1021-file limit)
3. **merge_vcfs** — merges each batch with `bcftools merge -m none` (avoids re-creating multiallelic sites)
4. **final_merge** — merges all batch VCFs into a single output; applies `bcftools +fill-tags`

MD5 checksums are computed at normalization, batch merge, and final merge stages.

Key config parameters (`bcftools/config.yaml`):
- `vcf_list_file`: text file with one `.vcf.gz` path per line
- `vcfs_per_batch`: batch size (default 100; max ~1021 due to bcftools limit)
- `final_filter_logic`: `"x"` = PASS if any input passes; `"+"` = PASS only if all inputs pass
- `info_rules`: comma-separated `FIELD:operation` pairs for bcftools merge INFO aggregation

### GATK pipeline (`gatk/`)

Processes per-sample GVCF files via joint genotyping:

1. **reblock_gvcfs** — runs `ReblockGVCF` on each sample to reduce file size
2. **scatter_intervals_by_ns** / **split_intervals** — splits the genome into `scatter_count` intervals for parallelism
3. **genomicsdb_import** — imports all samples into a GenomicsDB workspace per interval
4. **genotype_gvcfs** — runs `GenotypeGVCFs` per interval
5. **merge_vcfs** — merges interval VCFs with Picard `MergeVcfs`

Supports `design: "genome"` (uses ScatterIntervalsByNs) or `design: "exome"` (uses a BED file). Memory allocation doubles with each Snakemake retry attempt (up to 160 GB).

## Key Conventions

- Both pipelines use `TMPDIR` from the environment (set by Slurm) for scratch space, falling back to `/tmp`.
- Sample names are derived by stripping `.vcf.gz` / `.g.vcf.gz` suffixes from input file basenames.
- Output directories are created at pipeline startup (not by Snakemake rules).
- All rules log start/end timestamps and append stderr to log files.
- Conda environments (`bcftools`, `gatk`) must be pre-created with the required tools installed.
