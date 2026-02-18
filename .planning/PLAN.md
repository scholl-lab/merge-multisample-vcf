# Modernisation Plan — merge-multisample-vcf (bcftools pipeline)

**Created:** 2026-02-18
**Implemented:** 2026-02-18
**Status:** ✅ Complete — merged to `feat/modernise-bcftools-pipeline`, PR #5
**Baseline tag:** `v0.1.0` (pre-modernisation snapshot)

---

## Reference repositories analysed

| Repo | Path | Key takeaways applied |
|------|------|-----------------------|
| sm-alignment | `C:\development\scholl-lab\sm-alignment` | Snakemake 8 layout, profiles, conda YAMLs, `temp()`, cluster-auto-detecting launcher |
| sm-vcf-annotation | `C:\development\scholl-lab\sm-vcf-annotation` | `benchmark:` directive, `bcftools index -t`, `helpers.py`, `wildcard_constraints`, `localrule: True` |

---

## Bugs fixed

| # | Severity | Description | Resolution | Commit |
|---|----------|-------------|------------|--------|
| B1 | **FATAL** | `normalize_vcf` never indexed output — `bcftools merge` requires `.tbi` | Added `-W=tbi` to `bcftools norm`; declared `tbi` as named `temp()` output | `a3fb0bf` |
| B2 | moderate | `merge_vcfs` / `final_merge` created `.tbi` via `-W=tbi` but never declared it | Added named `tbi=` output in both rules | `a3fb0bf` |
| B3 | moderate | `md5_normalized_vcf` and `md5_merged_vcf` orphan rules — never executed | Removed all 3 standalone `md5_*` rules; folded `md5sum` inline at end of each producing rule | `a3fb0bf` |

## README TODOs resolved

| # | TODO | Resolution | Commit |
|---|------|------------|--------|
| T1 | MD5 checksums for all files | `md5sum` inline in all 3 rules; `.md5` declared as named output | `a3fb0bf` |
| T2 | Remove intermediate files | `temp()` on all normalized + batch VCFs; `--notemp` to override | `a3fb0bf` |
| T3 | Proper index/VCF creation sequence | `-W=tbi` on every bcftools output; `.tbi` always declared before any downstream rule | `a3fb0bf` |
| T4 | Error handling | `set -euo pipefail` in every shell block | `a3fb0bf` |
| T5 | Configurable file extensions | `vcf_suffix` config key (default `.vcf.gz`) used in `load_vcf_list()` | `a3fb0bf` / `daf77a9` |
| T6 | Dynamic memory allocation | `mem_mb(base)` retry-doubling callable in `helpers.py`; `retries: 2` in default profile | `4b72237` / `a3fb0bf` |

---

## Implemented structure

```
merge-multisample-vcf/
├── workflow/
│   ├── Snakefile                   ✅ Snakemake 8, min_version("8.0"), wildcard_constraints,
│   │                                  localrule all, get_final_outputs()
│   ├── rules/
│   │   ├── common.smk              ✅ Config extraction, dir creation, sample loading,
│   │   │                              input helper lambdas, parse-time batch-count guard
│   │   ├── merge.smk               ✅ 3 rules: normalize_vcf, merge_vcfs, final_merge
│   │   └── helpers.py              ✅ Pure Python: load_vcf_list, n_batches, batch_indices,
│   │                                  get_batch_samples, mem_mb — no Snakemake imports
│   └── envs/
│       └── bcftools.yaml           ✅ Pinned: bcftools+htslib >=1.18
├── config/
│   └── config_dummy.yaml           ✅ All 8 params documented; no conda_env key
├── profiles/
│   ├── default/config.yaml         ✅ latency-wait, retries=2, rerun-incomplete, conda,
│   │                                  shared-fs-usage, per-rule threads+resources+cpus_per_task
│   ├── bih/config.yaml             ✅ SLURM executor, medium partition, scholl-lab account
│   ├── charite/config.yaml         ✅ SLURM executor, compute partition, sc-users account
│   └── local/config.yaml           ✅ jobs=4, latency-wait=5
├── scripts/
│   └── run_snakemake.sh            ✅ detect_cluster() BIH/Charité/local, TMPDIR+EXIT trap,
│                                      --workflow-profile + --profile, "$@" passthrough
├── deprecated/gatk/                ✅ Archived (moved from gatk/), DEPRECATED.md added
├── bcftools/merge_sequential.smk   ✅ Kept unchanged as legacy fallback
├── .gitignore                      ✅ Fixed: path-specific ignores + !profiles/**/config.yaml
├── CLAUDE.md                       ✅ Rewritten for new layout
└── README.md                       ✅ Rewritten: overview, quick-start, config table,
                                       output layout, --notemp docs, profiles table
```

---

## Deviations from original plan — and why

### Deviation 1: No `multiext()` for normalize_vcf output

**Planned:** use `multiext(os.path.join(NORMALIZED_DIR, "{sample}.normalized"), ".vcf.gz", ".vcf.gz.tbi")`

**Implemented:** two named outputs (`vcf=`, `tbi=`)

**Reason:** `multiext()` makes downstream lambda references more verbose (indexing by position rather than by name). Named outputs (`output.vcf`, `output.tbi`) are clearer at the call site in `merge.smk` lambdas and match the sm-alignment style used in all sibling rules with multiple outputs.

---

### Deviation 2: `split_vcfs` and `list_merged_vcfs` eliminated entirely

**Planned:** fix bugs on existing rules, then potentially replace list-file pattern.

**Implemented:** both helper rules removed. Batch membership is computed by `get_batch_vcfs(idx)` / `get_batch_tbis(idx)` lambda functions in `common.smk`, which resolve to explicit file lists at DAG construction time.

**Reason:** the list-file pattern made `temp()` fundamentally unsafe — Snakemake cannot track dependencies that exist only inside text files. Removing these two helper rules (20 lines of pipeline code) unlocks `temp()` on all intermediate VCFs, which is the most important TODO. It also reduces rule count from 7 to 3, satisfying KISS.

The tradeoff: batch VCFs are now passed as positional arguments to `bcftools merge` instead of via `-l`. This limits each invocation to ~1021 files. A parse-time guard in `common.smk` exits with a clear error if `n_batches > 1000`, so the constraint is surfaced immediately rather than silently at runtime.

---

### Deviation 3: `keep_intermediates` config toggle not added

**Planned:** `keep_intermediates: false` config key; conditionally apply `temp()`.

**Implemented:** `temp()` always applied; `--notemp` Snakemake flag documented in README and config_dummy.yaml.

**Reason:** conditional `temp()` in Python requires passing a boolean through to every output declaration, adding noise to every rule. `--notemp` is the standard Snakemake mechanism for exactly this use case and is already well-known to practitioners. KISS.

---

### Deviation 4: Manual input-file existence check removed from `normalize_vcf`

**Planned:** keep the existing `if [ ! -f {input.vcf} ]; then ... exit 1; fi` guard.

**Implemented:** removed. `set -euo pipefail` covers this — bcftools will exit non-zero if the file is missing, and Snakemake will report the failure with the log path. The manual check was redundant and slightly misleading (it ran before Snakemake had verified inputs).

---

### Deviation 5: `validate()` / JSON schema not added

**Planned** (implied by sibling repo pattern): add `snakemake.utils.validate()` with a config schema.

**Not implemented** this cycle — the config is small enough that a comment-documented `config_dummy.yaml` covers the same ground without the schema maintenance burden. Can be added in a future cycle if the config grows.

---

### Deviation 6: `profiles/local/config.yaml` added (not in original plan)

**Not in plan.** Added because local dry-run is the recommended first step in the README quick-start, and having a named local profile makes the invocation consistent with cluster usage.

---

## Commits on `feat/modernise-bcftools-pipeline`

| Commit | Subject | Files |
|--------|---------|-------|
| `4b72237` | feat: add workflow/ skeleton — Snakefile, helpers, conda env | `workflow/Snakefile`, `workflow/rules/helpers.py`, `workflow/envs/bcftools.yaml` |
| `a3fb0bf` | fix: add common.smk + merge.smk — resolve all bugs and TODOs | `workflow/rules/common.smk`, `workflow/rules/merge.smk` |
| `daf77a9` | feat: add config, profiles, and cluster-aware launcher script | `config/config_dummy.yaml`, `scripts/run_snakemake.sh` |
| `ae54c5d` | fix: update .gitignore to track profiles/*/config.yaml | `.gitignore`, `profiles/*/config.yaml` (×4) |
| `3abd6da` | docs: rewrite README and CLAUDE.md for modernised pipeline | `README.md`, `CLAUDE.md` |

Earlier preparation commits (on `main` before branch):

| Commit | Subject |
|--------|---------|
| `aa36e85` | chore: deprecate GATK pipeline, add CLAUDE.md and planning docs |
| `0335865` | docs: update plan with sibling repo patterns |

---

## Files not modified (as planned)

| File | Status |
|------|--------|
| `bcftools/merge_sequential.smk` | Untouched — legacy fallback |
| `bcftools/run_merge_sequential.sh` | Untouched — legacy fallback |
| `bcftools/config_dummy.yaml` | Untouched — superseded by `config/config_dummy.yaml` |
| `deprecated/gatk/*` | Untouched — archived |

---

## Known follow-up items (future cycles)

| Item | Notes |
|------|-------|
| Config JSON schema | Add `workflow/schemas/config.schema.yaml` + `validate()` in Snakefile |
| End-to-end test | Small synthetic VCF fixture + CI action to run `--dry-run` on every PR |
| `bcftools/merge_sequential.smk` cleanup | Apply the B1–B3 fixes to the legacy file too, or formally deprecate it |
| `bcftools norm -W=tbi` version gate | bcftools < 1.16 does not support `-W` on norm; add a version check or fallback `tabix` call if older versions must be supported |
