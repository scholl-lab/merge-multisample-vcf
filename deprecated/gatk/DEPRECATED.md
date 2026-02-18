# GATK Joint-Genotyping Pipeline — DEPRECATED

**Status:** Archived / Not functional
**Archived:** 2026-02-18
**Reason:** Multiple blocking bugs; superseded by the bcftools pipeline

---

## Why this pipeline was deprecated

A review on 2026-02-18 identified two independently fatal bugs that prevent the pipeline from completing any run:

### Fatal Bug 1 — `reblock_gvcfs` never indexes its output (smk line 176)

`ReblockGVCF` produces a `.g.vcf.gz` file but the rule never calls `tabix` or `gatk IndexFeatureFile` on it. The sample map written to `lists/{sample}.sample_map` points `GenomicsDBImport` to these un-indexed files. GATK refuses to import un-indexed GVCFs — **the pipeline always crashes at `genomicsdb_import`**.

### Fatal Bug 2 — `genotype_gvcfs` never indexes its output (smk line 255)

Each `GenotypeGVCFs` rule produces `genotyped_variants.interval_{id}.vcf.gz` with no `.tbi`. Picard `MergeVcfs` requires a `.tbi` alongside every input — **the final merge always crashes**.

### Additional issues found

| Severity | Issue |
|----------|-------|
| Moderate | `merge_vcfs` (Picard) uses hard-coded `20000 MB` with no retry-escalation, unlike the other GATK rules |
| Moderate | Sample name is extracted from column 10 of the GVCF header — fragile for non-standard GVCFs |
| Design | All 400 per-interval intermediate VCFs are placed in `FINAL_DIR` alongside the final output |
| Design | No `temp()` wrappers; intermediate GenomicsDB workspaces and per-interval VCFs consume large disk space |

---

## If you need to revive this pipeline

The intended workflow is:

1. `reblock_gvcfs` → add `tabix {output.reblocked_gvcf}` + declare `.tbi` as output
2. `genomicsdb_import` → unchanged
3. `genotype_gvcfs` → add `gatk IndexFeatureFile -I {output}` + declare `.tbi` as output
4. `merge_vcfs` (Picard) → add retry-escalating memory; optionally replace with `bcftools concat` for speed

Files preserved as-is for reference:
- `reblock_genomicsdb_gatk_pipeline.smk`
- `config.yaml`
- `run_reblock_genomicsdb_gatk_pipeline.sh`
