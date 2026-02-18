# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the pipeline

**Dry-run (always do this first):**
```bash
snakemake -s workflow/Snakefile \
          --configfile config/config.yaml \
          --workflow-profile profiles/default \
          --profile profiles/local \
          --dry-run
```

**Local execution:**
```bash
snakemake -s workflow/Snakefile \
          --configfile config/config.yaml \
          --workflow-profile profiles/default \
          --profile profiles/local
```

**HPC (Slurm):**
```bash
sbatch scripts/run_snakemake.sh                        # uses config/config.yaml
sbatch scripts/run_snakemake.sh config/custom.yaml     # custom config
```

**Keep intermediate files (disable temp() cleanup):**
```bash
snakemake ... --notemp
```

## Repository structure

```
workflow/
├── Snakefile                  Entry point (Snakemake 8, min_version guard)
├── rules/
│   ├── common.smk             Config extraction, dir creation, input helpers
│   ├── merge.smk              All three pipeline rules
│   └── helpers.py             Pure Python helpers (no Snakemake imports)
└── envs/
    └── bcftools.yaml          Pinned conda env (bcftools+htslib >=1.18)
config/
├── config_dummy.yaml          Committed template — copy to config.yaml
└── config.yaml                Live config (gitignored)
profiles/
├── default/config.yaml        Workflow profile: resources, retries, conda
├── bih/config.yaml            BIH HPC SLURM
├── charite/config.yaml        Charité HPC SLURM
└── local/config.yaml          Local execution
scripts/
└── run_snakemake.sh           Cluster-auto-detecting Slurm launcher
bcftools/
└── merge_sequential.smk       Legacy single-file pipeline (Snakemake <8)
deprecated/gatk/               Archived GATK pipeline (non-functional)
.planning/PLAN.md              Implementation plan
```

## Architecture

The pipeline has three rules, each fixing a specific biological concern:

1. **`normalize_vcf`** — `bcftools norm` per sample: left-align indels, split multiallelic sites, produce bgzipped+tabix-indexed output. `temp()` output deleted after its batch merge completes.

2. **`merge_vcfs`** — `bcftools merge | +fill-tags | view` for each batch of `vcfs_per_batch` samples. Batch membership computed by `get_batch_vcfs(idx)` / `get_batch_tbis(idx)` lambdas, so Snakemake has full DAG visibility and `temp()` cleanup is safe. `-m none` preserves the normalized representation.

3. **`final_merge`** — same `bcftools merge | +fill-tags | view` pipeline over all batch VCFs. `expand()` over `BATCH_INDICES` ensures explicit dependency tracking.

**Key design decisions:**
- No list-file helper rules (`split_vcfs`, `list_merged_vcfs` from legacy pipeline). Explicit Snakemake inputs replace them, enabling safe `temp()`.
- All VCF outputs carry a `-W=tbi` tabix index, declared as a named output alongside the VCF.
- MD5 checksums are written inline at the end of each shell block (not in separate rules).
- Memory doubles on each Snakemake retry via `mem_mb(base)` from `helpers.py`.
- `set -euo pipefail` in every shell block.

## Key files to understand first

When making changes, read in this order:
1. `workflow/rules/helpers.py` — pure Python logic, independently testable
2. `workflow/rules/common.smk` — all config variables and helper wrappers
3. `workflow/rules/merge.smk` — the three rules
4. `profiles/default/config.yaml` — resource defaults and thread counts

## Config keys

| Key | Type | Description |
|-----|------|-------------|
| `vcf_list_file` | str | Path to newline-separated list of input VCFs |
| `vcfs_per_batch` | int | Max samples per batch (≤ 1 000) |
| `output_folder` | str | Root output directory |
| `reference_fasta` | str | Reference for bcftools norm |
| `vcf_suffix` | str | Suffix stripped to derive sample name (default `.vcf.gz`) |
| `final_output_name` | str | Filename for final VCF in `final/` |
| `final_filter_logic` | str | `x` or `+` (bcftools merge -F) |
| `info_rules` | str | Comma-separated `FIELD:OP` pairs (bcftools merge -i) |

## Deprecated / legacy

- `bcftools/merge_sequential.smk` — original monolithic Snakefile, kept for backward compatibility. Do not add features here.
- `deprecated/gatk/` — GATK joint-genotyping pipeline, archived 2026-02-18 due to fatal `.tbi` indexing bugs. See `deprecated/gatk/DEPRECATED.md`.
