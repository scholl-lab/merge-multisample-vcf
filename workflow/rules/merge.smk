"""Pipeline rules: normalize → batch-merge → final merge.

Design notes
------------
* All intermediate VCFs are marked temp() so Snakemake removes them once the
  last downstream rule consuming them has finished.  Run with --notemp to keep
  all files (useful for debugging or re-running individual stages).
* MD5 checksums are written inline at the end of each rule that produces a
  VCF, so they are always in sync with the file.
* Tabix (.tbi) indexes are declared as Snakemake outputs so the DAG can track
  them; bcftools creates them atomically via the -W=tbi flag.
* Batch VCFs are passed as positional args (not via -l list files), which
  makes Snakemake aware of every dependency and enables safe temp() cleanup.
  The bcftools ~1021-file limit is validated at parse time in common.smk.
"""


# ---------------------------------------------------------------------------
# Rule: normalize_vcf
# Fixes: B1 — index declared + produced; B3 — MD5 inline; T4 — pipefail
# ---------------------------------------------------------------------------

rule normalize_vcf:
    """Left-align, split multiallelic sites, and index each input VCF."""
    input:
        vcf=lambda wildcards: SAMPLE_TO_VCF[wildcards.sample],
        fasta=REFERENCE_FASTA,
    output:
        vcf=temp(os.path.join(NORMALIZED_DIR, "{sample}.normalized.vcf.gz")),
        tbi=temp(os.path.join(NORMALIZED_DIR, "{sample}.normalized.vcf.gz.tbi")),
        md5=os.path.join(LOG_DIR, "{sample}.normalized.vcf.gz.md5"),
    log:
        os.path.join(LOG_DIR, "normalize_vcf.{sample}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "normalize_vcf.{sample}.tsv"),
    threads: 4
    resources:
        mem_mb=mem_mb(4_000),
        tmpdir=SCRATCH_DIR,
    conda:
        "../envs/bcftools.yaml"
    shell:
        r"""
        set -euo pipefail
        echo "Starting normalize_vcf for {wildcards.sample} at: $(date)" > {log}

        bcftools norm \
            --threads {threads} \
            -m-any --force -a --atom-overlaps . \
            -f {input.fasta} \
            -Oz -o {output.vcf} \
            -W=tbi \
            {input.vcf} >> {log} 2>&1

        md5sum {output.vcf} > {output.md5} 2>> {log}
        echo "Finished normalize_vcf for {wildcards.sample} at: $(date)" >> {log}
        """


# ---------------------------------------------------------------------------
# Rule: merge_vcfs
# Fixes: B2 — tbi declared; B3 — MD5 inline; T2 — temp(); T4 — pipefail
# Removes: split_vcfs helper rule (batch membership computed via lambda)
# ---------------------------------------------------------------------------

rule merge_vcfs:
    """Merge one batch of normalised VCFs into a single multi-sample VCF.

    -m none  prevents bcftools from re-creating multiallelic sites so that
    the output stays normalised and does not require a second normalisation pass.
    """
    input:
        vcfs=lambda wc: get_batch_vcfs(int(wc.idx)),
        tbis=lambda wc: get_batch_tbis(int(wc.idx)),
    output:
        vcf=temp(os.path.join(MERGE_DIR, "merge.{idx}.vcf.gz")),
        tbi=temp(os.path.join(MERGE_DIR, "merge.{idx}.vcf.gz.tbi")),
        md5=os.path.join(LOG_DIR, "merge.{idx}.vcf.gz.md5"),
    log:
        os.path.join(LOG_DIR, "merge_vcfs.{idx}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "merge_vcfs.{idx}.tsv"),
    threads: 4
    resources:
        mem_mb=mem_mb(8_000),
        tmpdir=SCRATCH_DIR,
    params:
        info_rules=INFO_RULES,
        filter_logic=FINAL_FILTER_LOGIC,
    conda:
        "../envs/bcftools.yaml"
    shell:
        r"""
        set -euo pipefail
        echo "Starting merge_vcfs batch {wildcards.idx} at: $(date)" > {log}

        (
          bcftools merge \
            --threads {threads} \
            -F {params.filter_logic} \
            -0 -m none \
            -i '{params.info_rules}' \
            {input.vcfs} \
          | bcftools +fill-tags \
          | bcftools view --threads {threads} -Oz -o {output.vcf} -W=tbi
        ) >> {log} 2>&1

        md5sum {output.vcf} > {output.md5} 2>> {log}
        echo "Finished merge_vcfs batch {wildcards.idx} at: $(date)" >> {log}
        """


# ---------------------------------------------------------------------------
# Rule: final_merge
# Fixes: B2 — tbi declared; B3 — MD5 inline; T4 — pipefail
# Removes: list_merged_vcfs helper rule (explicit expand keeps DAG correct)
# ---------------------------------------------------------------------------

rule final_merge:
    """Merge all batch VCFs into the final cohort-wide VCF."""
    input:
        vcfs=expand(
            os.path.join(MERGE_DIR, "merge.{idx}.vcf.gz"),
            idx=BATCH_INDICES,
        ),
        tbis=expand(
            os.path.join(MERGE_DIR, "merge.{idx}.vcf.gz.tbi"),
            idx=BATCH_INDICES,
        ),
    output:
        vcf=os.path.join(FINAL_DIR, FINAL_OUTPUT_NAME),
        tbi=os.path.join(FINAL_DIR, f"{FINAL_OUTPUT_NAME}.tbi"),
        md5=os.path.join(FINAL_DIR, f"{FINAL_OUTPUT_NAME}.md5"),
    log:
        os.path.join(LOG_DIR, "final_merge.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "final_merge.tsv"),
    threads: 4
    resources:
        mem_mb=mem_mb(16_000),
        tmpdir=SCRATCH_DIR,
    params:
        info_rules=INFO_RULES,
        filter_logic=FINAL_FILTER_LOGIC,
    conda:
        "../envs/bcftools.yaml"
    shell:
        r"""
        set -euo pipefail
        echo "Starting final_merge at: $(date)" > {log}

        (
          bcftools merge \
            --threads {threads} \
            -F {params.filter_logic} \
            -0 -m none \
            -i '{params.info_rules}' \
            {input.vcfs} \
          | bcftools +fill-tags \
          | bcftools view --threads {threads} -Oz -o {output.vcf} -W=tbi
        ) >> {log} 2>&1

        md5sum {output.vcf} > {output.md5} 2>> {log}
        echo "Finished final_merge at: $(date)" >> {log}
        """
