# Modernisation Plan — merge-multisample-vcf (bcftools pipeline)

**Created:** 2026-02-18
**Updated:** 2026-02-18 (sibling repo analysis: sm-alignment, sm-vcf-annotation)
**Scope:** Fix all bugs, resolve all README TODOs, modernise to match scholl-lab Snakemake 8 conventions
**Out of scope:** GATK pipeline (archived in `deprecated/gatk/`)

---

## Reference repositories analysed

| Repo | Path | Key takeaways |
|------|------|---------------|
| sm-alignment | `C:\development\scholl-lab\sm-alignment` | Modern Snakemake 8 layout, profiles, conda YAMLs, temp() everywhere, cluster-auto-detecting launcher |
| sm-vcf-annotation | `C:\development\scholl-lab\sm-vcf-annotation` | benchmark directive, bcftools index -t pattern, scatter/gather, helpers.py, wildcard_constraints |

---

## Bugs to fix (identified in review 2026-02-18)

| # | Severity | Description | Location |
|---|----------|-------------|----------|
| B1 | **FATAL** | `normalize_vcf` never indexes output — `bcftools merge -l` requires `.tbi` for every listed file | `smk:102–127` |
| B2 | moderate | `merge_vcfs` and `final_merge` create `.tbi` via `-W=tbi` but don't declare it as Snakemake output | `smk:172, 242` |
| B3 | moderate | `md5_normalized_vcf` and `md5_merged_vcf` are orphan rules — not wired into `rule all`, never execute | `smk:130, 202` |

## README TODOs to resolve

| # | TODO |
|---|------|
| T1 | MD5 checksums for **all** files (currently only final output) |
| T2 | Remove intermediate files to save storage |
| T3 | Proper sequence of index and VCF file creation |
| T4 | Error handling for missing inputs / failed commands |
| T5 | Configurable file extensions (hard-coded `.vcf.gz` throughout) |
| T6 | Dynamic memory allocation by thread count |

---

## Target structure (matching sibling repos)

```
merge-multisample-vcf/
├── workflow/
│   ├── Snakefile                   # Entry point (Snakemake 8+)
│   ├── rules/
│   │   ├── common.smk              # Config extraction, helper imports, dir creation
│   │   └── merge.smk               # All pipeline rules
│   ├── envs/
│   │   └── bcftools.yaml           # Pinned conda env (replaces conda_env config key)
│   └── rules/
│       └── helpers.py              # Pure Python helpers (sample extraction, etc.)
├── config/
│   ├── config.yaml                 # Active config (gitignored)
│   └── config_dummy.yaml           # Template config (committed)
├── profiles/
│   ├── default/
│   │   └── config.yaml             # Resource defaults, retries, latency-wait
│   ├── bih/
│   │   └── config.yaml             # BIH/CUBI SLURM executor settings
│   └── charite/
│       └── config.yaml             # Charité SLURM settings
├── scripts/
│   └── run_snakemake.sh            # Cluster-auto-detecting launcher (replaces bcftools/run_merge_sequential.sh)
├── deprecated/
│   └── gatk/                       # Archived GATK pipeline
├── .planning/
│   └── PLAN.md
├── CLAUDE.md
└── README.md
```

---

## Phase 1 — Critical bug fixes

**Goal:** Pipeline runs end-to-end without crashing.

### 1.1  Fix B1 — index normalized VCFs

In `normalize_vcf`, `bcftools norm` outputs a compressed VCF but never indexes it.
Fix: add `-W=tbi` to the `bcftools norm` command (supported since bcftools 1.16) **and** declare the index as a named output.

Use `multiext()` to declare both files (pattern from sm-vcf-annotation):
```python
output:
    multiext(
        os.path.join(NORMALIZED_DIR, "{sample}.normalized"),
        ".vcf.gz",
        ".vcf.gz.tbi",
    )
```

Fallback if the bcftools version doesn't support `-W` on norm: append `tabix {output[0]}` after the norm call.

### 1.2  Fix B2 — declare `.tbi` outputs in merge_vcfs and final_merge

Both rules already write the index via `bcftools view -W=tbi`. Add it to `output:`:
```python
output:
    vcf=os.path.join(MERGE_DIR, "merge.{idx}.vcf.gz"),
    tbi=os.path.join(MERGE_DIR, "merge.{idx}.vcf.gz.tbi"),
```

### 1.3  Fix B3 — wire MD5 rules or fold them inline

**Decision:** Fold MD5 into each producing rule (matches sm-alignment `utilities.smk` pattern of keeping checksums close to their file). Remove standalone `md5_normalized_vcf` and `md5_merged_vcf` rules.

- In `normalize_vcf` shell block: append `md5sum {output[0]} > {output[0]}.md5`
- In `merge_vcfs` shell block: append `md5sum {output.vcf} > {output.vcf}.md5`
- In `final_merge` shell block: append `md5sum {output.vcf} > {output.vcf}.md5`
- Declare `.md5` files as named outputs so Snakemake tracks them
- Remove `rule md5_normalized_vcf`, `rule md5_merged_vcf`, `rule md5_final_merge`
- Update `rule all` to only request the final VCF and its `.md5`

---

## Phase 2 — README TODOs

### 2.1  T2 — `temp()` wrappers for intermediate files

Wrap files that are no longer needed once the next stage starts:

```python
# normalized VCFs consumed after split_vcfs writes the list
output:
    multiext(temp(os.path.join(NORMALIZED_DIR, "{sample}.normalized")),
             ".vcf.gz", ".vcf.gz.tbi")

# batch merge VCFs consumed after final_merge
output:
    vcf=temp(os.path.join(MERGE_DIR, "merge.{idx}.vcf.gz")),
    tbi=temp(os.path.join(MERGE_DIR, "merge.{idx}.vcf.gz.tbi")),
```

Make `temp()` opt-out via config:
```yaml
keep_intermediates: false   # set true to retain normalized_vcfs/ and merged_vcfs/
```

### 2.2  T5 — configurable input file extension

Add `vcf_suffix` config key (default `".vcf.gz"`). Use it everywhere a suffix is stripped or appended:
```python
VCF_SUFFIX = config.get("vcf_suffix", ".vcf.gz")
samples = [os.path.basename(f).replace(VCF_SUFFIX, '') for f in vcf_list]
```

### 2.3  T6 — retry-escalating memory (from sm-alignment pattern)

Replace `get_mem_from_threads` with a proper retry function:
```python
def get_mem_mb(base_mb):
    def mem(wildcards, attempt):
        return min(base_mb * (2 ** (attempt - 1)), 128_000)
    return mem
```

Per-rule bases:
- `normalize_vcf`: 4 000 MB
- `merge_vcfs`: 8 000 MB
- `final_merge`: 16 000 MB

Add `retries: 2` in the default profile (not in the Snakefile).

### 2.4  T4 — error handling

- Add `set -euo pipefail` to all shell blocks (currently only the launcher uses this)
- The existing manual missing-file check in `normalize_vcf` is good; keep it

---

## Phase 3 — Modernisation (Snakemake 8 + sibling conventions)

### 3.1  Snakemake 8 entry point

Create `workflow/Snakefile`:
```python
from snakemake.utils import min_version, validate
min_version("8.0")

configfile: "config/config.yaml"

wildcard_constraints:
    sample="[A-Za-z0-9_.+-]+",
    idx="[0-9]+",

include: "rules/common.smk"
include: "rules/merge.smk"
```

Keep `bcftools/merge_sequential.smk` as the single-file legacy option during transition; deprecate once the new layout is working.

### 3.2  Split rules into logical files

`workflow/rules/common.smk` — config extraction, directory creation, sample loading, helper imports
`workflow/rules/merge.smk` — all pipeline rules (normalize → split → batch merge → final merge)

`workflow/rules/helpers.py` — pure Python functions:
- `load_vcf_list(path, suffix)` → `(samples, sample_to_vcf)`
- `get_mem_mb(base_mb)` → retry memory callable
- `n_batches(total, batch_size)` → ceiling division

### 3.3  Conda environment YAML (replaces string env name)

`workflow/envs/bcftools.yaml`:
```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - bcftools >=1.18
  - htslib >=1.18
```

Change all `conda: CONDA_ENV` to `conda: "../envs/bcftools.yaml"` (relative to rules file).
Remove `conda_env` from config (no longer needed).

### 3.4  Profiles

`profiles/default/config.yaml` (workflow profile — ships with repo):
```yaml
latency-wait: 60
rerun-incomplete: true
retries: 2
printshellcmds: true
software-deployment-method:
  - conda

default-resources:
  mem_mb: 4000
  runtime: 120

set-threads:
  normalize_vcf: 4
  merge_vcfs: 4
  final_merge: 4

set-resources:
  normalize_vcf:
    mem_mb: 4000
    runtime: 120
  merge_vcfs:
    mem_mb: 8000
    runtime: 1440
  final_merge:
    mem_mb: 16000
    runtime: 1440
```

`profiles/bih/config.yaml` (BIH/CUBI SLURM):
```yaml
executor: slurm
default-resources:
  slurm_partition: medium
  slurm_account: <account>
jobs: 200
```

`profiles/charite/config.yaml` (Charité SLURM):
```yaml
executor: slurm
default-resources:
  slurm_partition: compute
  slurm_account: sc-users
jobs: 200
```

### 3.5  Cluster-auto-detecting launcher

`scripts/run_snakemake.sh` — replaces `bcftools/run_merge_sequential.sh`.
Copies the cluster detection pattern verbatim from sm-alignment:
```bash
detect_cluster() { ... }  # BIH vs Charite vs local
```
Uses `--workflow-profile profiles/default` + `--profile profiles/${CLUSTER}`.
Add `--retries 2` passthrough (already handled in default profile).

### 3.6  benchmark directive

Add to all shell rules:
```python
benchmark:
    os.path.join(LOG_DIR, "benchmarks/{rule_name}.{wildcard}.tsv")
```
(Pattern from sm-vcf-annotation — all major rules have benchmarks.)

### 3.7  localrule on rule all

```python
rule all:
    localrule: True
    input:
        os.path.join(FINAL_DIR, FINAL_OUTPUT_NAME),
        os.path.join(FINAL_DIR, f"{FINAL_OUTPUT_NAME}.md5"),
```

### 3.8  README overhaul

- Quick start section with minimal config
- Document all 8 config parameters (currently only 4 documented)
- Separate sections for local run and HPC (Slurm)
- Remove TODO section (resolved)
- Add section on output directory structure

### 3.9  Update CLAUDE.md

Reflect new `workflow/` layout, `profiles/`, `scripts/`, resolved TODOs, deprecated GATK.

---

## Execution order

```
Phase 1 (must do first — pipeline broken without these)
  1.1 norm index fix
  1.2 tbi output declaration
  1.3 MD5 inline + orphan rule removal

Phase 2 (independent of each other, do after Phase 1)
  2.1 temp() wrappers + keep_intermediates config
  2.2 vcf_suffix config
  2.3 retry memory
  2.4 set -euo pipefail

Phase 3 (modernisation — can overlap with Phase 2)
  3.1 Snakemake 8 entry point + workflow/ layout
  3.2 Split rules files
  3.3 Conda YAML
  3.4 Profiles
  3.5 New launcher script
  3.6 benchmark directive
  3.7 localrule all
  3.8 README
  3.9 CLAUDE.md
```

Phases 1 and 2 can be done on the existing `bcftools/merge_sequential.smk` directly.
Phase 3 creates the new `workflow/` layout; the old `bcftools/merge_sequential.smk` remains as a legacy fallback until the new layout is validated.

---

## Files to be created / modified

| File | Action | Phase |
|------|--------|-------|
| `bcftools/merge_sequential.smk` | Fix B1–B3, add temp(), set -euo pipefail, retry memory | 1, 2 |
| `bcftools/config_dummy.yaml` | Add vcf_suffix, keep_intermediates; document all params | 2 |
| `workflow/Snakefile` | New — Snakemake 8 entry point | 3.1 |
| `workflow/rules/common.smk` | New — config + dir setup | 3.2 |
| `workflow/rules/merge.smk` | New — all rules (ported from merge_sequential.smk) | 3.2 |
| `workflow/rules/helpers.py` | New — pure Python helpers | 3.2 |
| `workflow/envs/bcftools.yaml` | New — pinned conda env | 3.3 |
| `profiles/default/config.yaml` | New — resource defaults | 3.4 |
| `profiles/bih/config.yaml` | New — BIH SLURM | 3.4 |
| `profiles/charite/config.yaml` | New — Charité SLURM | 3.4 |
| `scripts/run_snakemake.sh` | New — replaces bcftools/run_merge_sequential.sh | 3.5 |
| `README.md` | Rewrite | 3.8 |
| `CLAUDE.md` | Update | 3.9 |

## Files NOT to be modified

| File | Reason |
|------|--------|
| `deprecated/gatk/*` | Archived |
| `bcftools/config.yaml` | Live config; gitignored |
| `bcftools/run_merge_sequential.sh` | Kept as legacy; superseded by scripts/run_snakemake.sh |
