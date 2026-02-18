# merge-multisample-vcf

Snakemake 8 pipeline for normalising and merging large multi-sample VCF cohorts
using **bcftools**.  Handles cohorts too large for a single `bcftools merge` call
(> ~1 000 files) by batching in two stages.

## Pipeline overview

```
input VCFs (one per sample)
       │
       ▼
 normalize_vcf          — bcftools norm: left-align, split multiallelic, index
       │
       ▼  (batches of vcfs_per_batch)
  merge_vcfs             — bcftools merge | +fill-tags | view (per batch)
       │
       ▼
  final_merge            — bcftools merge | +fill-tags | view (all batches)
       │
       ▼
 {output_folder}/final/{final_output_name}  +  .tbi  +  .md5
```

Intermediate VCFs are automatically removed (`temp()`) once consumed.
Run with `--notemp` to keep them for debugging.

## Quick start

**1. Copy and edit the config template:**

```bash
cp config/config_dummy.yaml config/config.yaml
# Edit config/config.yaml — at minimum set vcf_list_file and reference_fasta
```

**2. Create the VCF list file** (one absolute path per line, blank lines and `#` comments ignored):

```
/data/cohort/sample1.vcf.gz
/data/cohort/sample2.vcf.gz
...
```

**3. Run locally (dry-run first):**

```bash
snakemake -s workflow/Snakefile \
          --configfile config/config.yaml \
          --workflow-profile profiles/default \
          --profile profiles/local \
          --dry-run

# Remove --dry-run to execute
```

**4. Run on an HPC cluster (BIH/Charité via Slurm):**

```bash
sbatch scripts/run_snakemake.sh                         # uses config/config.yaml
sbatch scripts/run_snakemake.sh config/my_config.yaml   # custom config
```

The script auto-detects the cluster (BIH or Charité) and loads the matching profile.

## Configuration

Copy `config/config_dummy.yaml` to `config/config.yaml` (gitignored) and set:

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `vcf_list_file` | ✓ | — | Path to file listing input VCFs (one per line) |
| `reference_fasta` | ✓ | — | Reference FASTA for `bcftools norm` |
| `output_folder` | ✓ | — | Root output directory |
| `vcfs_per_batch` | ✓ | — | Samples per intermediate batch (keep ≤ 1 000) |
| `vcf_suffix` | | `.vcf.gz` | Extension stripped to derive sample names |
| `final_output_name` | | `all_merged.vcf.gz` | Final output filename |
| `final_filter_logic` | | `x` | `x` = PASS if any sample passes; `+` = all must pass |
| `info_rules` | | *(GATK defaults)* | `FIELD:OPERATION` pairs for `bcftools merge -i` |

See `config/config_dummy.yaml` for full documentation and the default `info_rules` string.

## Output structure

```
{output_folder}/
├── normalized_vcfs/        per-sample normalised VCFs (temp — auto-deleted)
├── merged_vcfs/            per-batch merged VCFs      (temp — auto-deleted)
├── final/
│   ├── {final_output_name}       final cohort VCF
│   ├── {final_output_name}.tbi   tabix index
│   └── {final_output_name}.md5   MD5 checksum
└── logs/
    ├── normalize_vcf.{sample}.log
    ├── merge_vcfs.{idx}.log
    ├── final_merge.log
    ├── {sample}.normalized.vcf.gz.md5   per-sample checksums
    ├── merge.{idx}.vcf.gz.md5           per-batch checksums
    └── benchmarks/                       Snakemake timing TSVs
```

## Keeping intermediate files

By default all intermediate VCFs are deleted after use.
To keep them (e.g. to re-run only the final merge stage):

```bash
snakemake ... --notemp
```

## Profiles

| Profile | Purpose |
|---------|---------|
| `profiles/default/` | Resource defaults — always loaded via `--workflow-profile` |
| `profiles/bih/` | BIH HPC SLURM settings (auto-selected by launcher) |
| `profiles/charite/` | Charité HPC SLURM settings (auto-selected by launcher) |
| `profiles/local/` | Local execution (low parallelism, short latency-wait) |

Cluster profiles require `snakemake-executor-plugin-slurm`:
```bash
pip install snakemake-executor-plugin-slurm
```

## Conda environment

The pipeline uses `workflow/envs/bcftools.yaml` (bcftools + htslib ≥ 1.18).
Snakemake creates and caches the environment automatically when
`software-deployment-method: conda` is set in the profile (default).

## Archived pipelines

The GATK joint-genotyping pipeline has been archived to `deprecated/gatk/`
due to unresolved indexing bugs; see `deprecated/gatk/DEPRECATED.md`.
