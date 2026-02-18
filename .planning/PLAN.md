# Modernisation Plan — merge-multisample-vcf (bcftools pipeline)

**Created:** 2026-02-18
**Scope:** Fix all bugs, resolve all README TODOs, modernise the bcftools pipeline
**Out of scope:** GATK pipeline (archived in `deprecated/gatk/`)

---

## Background

### Bugs identified in review (2026-02-18)

| # | Severity | Description | Location |
|---|----------|-------------|----------|
| B1 | FATAL | `normalize_vcf` never indexes its output — bcftools merge requires `.tbi` for every input in a list | `merge_sequential.smk:102–127` |
| B2 | moderate | `merge_vcfs` and `final_merge` create `.tbi` via `-W=tbi` but don't declare it as output — Snakemake can't track or clean it | `smk:172, 242` |
| B3 | moderate | `md5_normalized_vcf` and `md5_merged_vcf` are orphan rules — never wired into `rule all`, so they never execute | `smk:130, 202` |

### README TODOs (from `README.md`)

| # | TODO |
|---|------|
| T1 | Implement MD5 checksum calculation for **all** files (currently only final output is verified) |
| T2 | Remove intermediate files to save storage space |
| T3 | Ensure proper sequence of index and VCF file creation |
| T4 | Add error handling for missing input files or failed commands |
| T5 | Make file extensions configurable (currently hard-coded `.vcf.gz`) |
| T6 | Explore dynamic memory allocation based on thread count |

---

## Phases

### Phase 1 — Critical bug fixes (pipeline correctness)

**Goal:** Make the pipeline runnable end-to-end without crashing.

#### 1.1 — Fix B1: index normalized VCFs

In `normalize_vcf`:
- The `bcftools norm` command already uses `-Oz`; bcftools ≥ 1.16 supports `-W=tbi` on norm as well. Add `-W=tbi` to the norm call **or** add an explicit `tabix {output.normalized}` line.
- Declare the index as a named secondary output: `index=os.path.join(NORMALIZED_DIR, "{sample}.normalized.vcf.gz.tbi")`

Alternatively use `multiext()` to declare both in one line:
```python
output:
    multiext(os.path.join(NORMALIZED_DIR, "{sample}.normalized"), ".vcf.gz", ".vcf.gz.tbi")
```

#### 1.2 — Fix B2: declare `.tbi` outputs in `merge_vcfs` and `final_merge`

Both rules already produce `.tbi` via `-W=tbi`. Add the index to `output:` so Snakemake tracks it:
```python
output:
    vcf=os.path.join(MERGE_DIR, "merge.{idx}.vcf.gz"),
    tbi=os.path.join(MERGE_DIR, "merge.{idx}.vcf.gz.tbi")
```

#### 1.3 — Fix B3: wire intermediate MD5 rules into DAG

Option A (simpler): add the MD5 outputs to `rule all` so they are always computed.
Option B (cleaner): integrate `md5sum` calls into the end of each rule's shell block and remove the standalone MD5 rules.

**Decision:** Use Option B — fold MD5 into the rule that produces the file. This avoids unnecessary rule proliferation and ensures checksums are always current.

---

### Phase 2 — Address README TODOs

#### 2.1 — T2: `temp()` wrappers for intermediate files

Wrap intermediate outputs that are no longer needed after the next stage:
- `normalized_vcfs/{sample}.normalized.vcf.gz` + `.tbi` → `temp()` after `split_vcfs` produces the list
- `merged_vcfs/merge.{idx}.vcf.gz` + `.tbi` → `temp()` after `final_merge` completes

Benefit: large cohorts (thousands of samples) produce hundreds of GB of intermediate VCFs.

**Caution:** `temp()` deletes files even on partial runs. Evaluate whether users re-run from checkpoints. Consider making this opt-in via a config flag `keep_intermediates: true`.

#### 2.2 — T5: configurable input file extensions

Currently the pipeline hard-codes `.vcf.gz` in several places:
- `samples = [os.path.basename(f).replace('.vcf.gz', '') for f in vcf_list]`
- Normalized output path: `{sample}.normalized.vcf.gz`

Add a `vcf_suffix` config parameter (default `".vcf.gz"`) used everywhere file extension stripping occurs.

#### 2.3 — T6: dynamic memory allocation

`get_mem_from_threads` currently returns `threads * 1000`. This is too low for large merges.

Replace with a proper retry-escalating function matching the GATK pattern:
```python
def get_mem_mb(base_mb):
    def mem(wildcards, attempt):
        return min(base_mb * (2 ** (attempt - 1)), 128000)
    return mem
```

Each resource-heavy rule (`normalize_vcf`, `merge_vcfs`, `final_merge`) gets a sensible base:
- `normalize_vcf`: 4 000 MB base
- `merge_vcfs`: 8 000 MB base
- `final_merge`: 16 000 MB base

Also expose `--retries` in the Slurm script (currently absent for bcftools).

#### 2.4 — T4: error handling

Add `set -euo pipefail` to all shell blocks (currently only GATK rules use `set -e`).
The `normalize_vcf` rule already has a manual check for missing input — consolidate this pattern across all rules.

---

### Phase 3 — Modernisation

#### 3.1 — Conda environment YAML files

Replace the string `conda_env: "bcftools"` config approach with a pinned YAML file at `envs/bcftools.yaml`:
```yaml
channels:
  - bioconda
  - conda-forge
dependencies:
  - bcftools >=1.18
```

This makes the environment reproducible and removes the requirement for users to pre-create named environments.

#### 3.2 — Snakemake profile support

Add a `profiles/` directory with a `config.yaml` for local + cluster execution, so users don't need to pass `--profile=cubi-v1` manually.

#### 3.3 — README overhaul

Rewrite `README.md` to:
- Document both execution modes (local dry-run + Slurm)
- Document all config parameters (currently only 4 of 8 are documented)
- Remove the TODO section (todos resolved or tracked in `.planning/`)
- Add a "Quick start" section with the minimal config needed

#### 3.4 — Update CLAUDE.md

Reflect all changes: new `envs/` dir, updated config keys, resolved TODOs, deprecated GATK pipeline.

---

## Execution order

```
Phase 1 (bugs)     →  Phase 2 (TODOs)  →  Phase 3 (modernisation)
  1.1 norm index         2.1 temp()            3.1 conda YAML
  1.2 tbi outputs        2.2 extensions        3.2 profiles
  1.3 MD5 inline         2.3 memory            3.3 README
                         2.4 error handling    3.4 CLAUDE.md
```

Phases 2 and 3 tasks are mostly independent and can be done in any order once Phase 1 is complete.

---

## Files to be modified

| File | Changes |
|------|---------|
| `bcftools/merge_sequential.smk` | All Phase 1 + 2 fixes; main target |
| `bcftools/config_dummy.yaml` | Add `vcf_suffix`, document all params |
| `bcftools/run_merge_sequential.sh` | Add `--retries` flag |
| `envs/bcftools.yaml` | New file — pinned conda environment |
| `README.md` | Full rewrite (Phase 3.3) |
| `CLAUDE.md` | Update to reflect new structure |

## Files NOT to be modified

| File | Reason |
|------|--------|
| `deprecated/gatk/*` | Archived; not touched |
| `bcftools/config.yaml` | Live config; not committed |
