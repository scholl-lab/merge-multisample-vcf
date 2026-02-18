# merge-multisample-vcf

Snakemake 8 pipeline for normalising and merging large multi-sample VCF cohorts
with **bcftools**. Two-stage batching handles cohorts > ~1 000 samples.

```
input VCFs  →  normalize_vcf  →  merge_vcfs (per batch)  →  final_merge
                bcftools norm     bcftools merge|+fill-tags    bcftools merge|+fill-tags
```

Intermediate files are auto-deleted (`temp()`). Use `--notemp` to keep them.

## Quick start

**Generate config interactively:**
```bash
python scripts/generate_config.py
```
Or with flags (see `--help` for all options):
```bash
python scripts/generate_config.py --vcf-folder /data/cohort/vcfs --ref /ref/GRCh38.fa
```

**Dry-run:**
```bash
snakemake -s workflow/Snakefile --configfile config/config.yaml \
          --workflow-profile profiles/default --profile profiles/local --dry-run
```

**HPC (BIH / Charité — auto-detected):**
```bash
sbatch scripts/run_snakemake.sh                        # default config
sbatch scripts/run_snakemake.sh config/my_config.yaml  # custom config
```

## Configuration

| Parameter | Required | Default | Description |
|-----------|:--------:|---------|-------------|
| `vcf_list_file` | ✓ | — | File listing input VCFs, one path per line |
| `reference_fasta` | ✓ | — | Reference FASTA for `bcftools norm` |
| `output_folder` | ✓ | — | Root output directory |
| `vcfs_per_batch` | ✓ | — | VCFs per intermediate batch (keep ≤ 1 000) |
| `vcf_suffix` | | `.vcf.gz` | Suffix stripped to derive sample names |
| `final_output_name` | | `all_merged.vcf.gz` | Final output filename |
| `final_filter_logic` | | `x` | `x` = PASS if any sample passes; `+` = all must pass |
| `info_rules` | | *(GATK defaults)* | `FIELD:OP` pairs for `bcftools merge -i` |

See `config/config_dummy.yaml` for the full `info_rules` default and inline docs.

## Output

```
{output_folder}/
├── normalized_vcfs/   per-sample normalised VCFs  (temp)
├── merged_vcfs/       per-batch merged VCFs        (temp)
├── final/
│   ├── {final_output_name}
│   ├── {final_output_name}.tbi
│   └── {final_output_name}.md5
└── logs/              per-rule logs, md5 files, benchmarks/
```

## Dev tooling

```bash
pip install ruff mypy snakefmt shellcheck-py pytest
make lint     # ruff + snakefmt + shellcheck + mypy
make format   # auto-format all files
make test     # pytest
```

## Requirements

- Snakemake ≥ 8, conda
- `snakemake-executor-plugin-slurm` for HPC cluster execution
- bcftools ≥ 1.18 (managed automatically via `workflow/envs/bcftools.yaml`)
