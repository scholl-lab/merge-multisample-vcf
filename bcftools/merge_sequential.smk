# merge_sequential.smk

import os
import functools
from pathlib import Path

# ----------------------------------------------------------------------------------- #
# SCRIPT DESCRIPTION:
# This Snakemake script orchestrates a pipeline to normalize, split, and merge large sets
# of VCF files. It normalizes each VCF, splits them into smaller subsets, merges them in batches, 
# and then merges all batches into a single VCF file. `bcftools` is utilized for 
# VCF processing tasks, which has a limitation of handling up to 1021 files at a time.
# Additionally, MD5 checksums are computed for key output files to verify data integrity.
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# ENVIRONMENT SETUP:
SCRATCH_DIR = os.environ.get('TMPDIR', '/tmp')  # Default to /tmp if TMPDIR not set
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# CONFIGURATION:
configfile: "config.yaml"
VCF_LIST_FILE = config["vcf_list_file"]          # Path to the file listing input VCFs
VCFS_PER_BATCH = config["vcfs_per_batch"]        # Number of VCFs per merge batch
OUTPUT_FOLDER = config["output_folder"]          # Output directory
CONDA_ENV = config["conda_env"]                  # Conda environment
REFERENCE_FASTA = config["reference_fasta"]      # Reference FASTA for normalization

# New Parameter: Final Output Name
# Set a default value if 'final_output_name' is not provided
FINAL_OUTPUT_NAME = config.get("final_output_name", "all_merged.vcf.gz")
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# OUTPUT ORGANIZATION:
prefix_results = functools.partial(os.path.join, OUTPUT_FOLDER)

# Define subdirectories for different stages of processing.
NORMALIZED_DIR = prefix_results('normalized_vcfs')  # Directory for normalized VCFs
LISTS_DIR = prefix_results('lists')                # Directory for list files
MERGE_DIR = prefix_results('merged_vcfs')          # Directory for merged VCFs
FINAL_DIR = prefix_results('final')                # Directory for final merged VCF
LOG_DIR = prefix_results('logs')                    # Directory for logs

# Ensure all output directories exist.
for output_dir in [NORMALIZED_DIR, LISTS_DIR, MERGE_DIR, FINAL_DIR, LOG_DIR]:
    os.makedirs(output_dir, exist_ok=True)
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# SAMPLE PROCESSING:
# Read the list of VCF files and extract sample names.
# ----------------------------------------------------------------------------------- #

# Read the list of VCF files from the input file.
with open(VCF_LIST_FILE, 'r') as f:
    vcf_list = [line.strip() for line in f if line.strip()]

# Extract sample names from the basenames of the input files.
samples = [os.path.basename(f).replace('.vcf.gz', '') for f in vcf_list]
samples = [s[:-4] if s.endswith('.vcf') else s for s in samples]  # Remove '.vcf' if present

# Create a mapping from sample names to VCF file paths.
sample_to_vcf = dict(zip(samples, vcf_list))

# Calculate the number of subsets needed for merging.
total_vcfs = len(vcf_list)
needed_number_of_subsets = -(-total_vcfs // VCFS_PER_BATCH)  # Ceiling division
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# HELPER FUNCTIONS:
def get_mem_from_threads(wildcards, threads):
    return threads * 1000  # Simple memory allocation per thread (in MB)
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# PIPELINE RULES:

rule all:
    input:
        # Final merged VCF with desired name
        os.path.join(FINAL_DIR, FINAL_OUTPUT_NAME),
        # MD5 checksum for the final merged VCF
        os.path.join(FINAL_DIR, f"{FINAL_OUTPUT_NAME}.md5")
        # You can add more checksum files here if needed

# Rule to normalize each VCF file
rule normalize_vcf:
    input:
        vcf=lambda wildcards: sample_to_vcf.get(wildcards.sample, None),
        fasta=REFERENCE_FASTA
    output:
        normalized=os.path.join(NORMALIZED_DIR, "{sample}.normalized.vcf.gz")
    log:
        os.path.join(LOG_DIR, "normalize.{sample}.log")
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time = '02:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        CONDA_ENV
    shell:
        """
        if [ ! -f {input.vcf} ]; then
            echo "VCF file {input.vcf} not found!" >&2
            exit 1
        fi

        echo "Starting normalization for {wildcards.sample} at $(date)" > {log}
        bcftools norm -m-any --force -a --atom-overlaps . -W --threads {threads} -f {input.fasta} -Oz -o {output.normalized} {input.vcf} &>> {log}
        echo "Finished normalization for {wildcards.sample} at $(date)" >> {log}
        """

# Rule to compute MD5 checksum for normalized VCFs
rule md5_normalized_vcf:
    input:
        normalized=os.path.join(NORMALIZED_DIR, "{sample}.normalized.vcf.gz")
    output:
        md5=os.path.join(LOG_DIR, "{sample}.normalized.vcf.gz.md5")
    log:
        log_file=os.path.join(LOG_DIR, "md5_normalized.{sample}.log")
    # Removed 'conda' directive as 'md5sum' is a standard utility
    shell:
        """
        echo "Computing MD5 for {input.normalized} at $(date)" > {log_file}
        md5sum {input.normalized} > {output.md5} 2>> {log_file}
        echo "Finished MD5 for {input.normalized} at $(date)" >> {log_file}
        """

# Rule to split the normalized VCF files into subsets
rule split_vcfs:
    input:
        vcf_files=expand(os.path.join(NORMALIZED_DIR, "{sample}.normalized.vcf.gz"), sample=samples)
    output:
        subset_files=expand(os.path.join(LISTS_DIR, "subset_vcfs.{idx}"), idx=range(needed_number_of_subsets))
    log:
        subset_logs=expand(os.path.join(LOG_DIR, "split_vcfs.{idx}.log"), idx=range(needed_number_of_subsets))
    run:
        from datetime import datetime

        for idx, subset_file, subset_log in zip(range(needed_number_of_subsets), output.subset_files, log.subset_logs):
            with open(subset_log, 'a') as log_file:
                log_file.write(f"Starting split_vcfs subset {idx} at: {datetime.now()}\n")

            subset = input.vcf_files[idx * VCFS_PER_BATCH : (idx + 1) * VCFS_PER_BATCH]
            with open(subset_file, 'w') as f:
                f.write("\n".join(subset) + "\n")

            with open(subset_log, 'a') as log_file:
                log_file.write(f"Finished split_vcfs subset {idx} at: {datetime.now()}\n")

# Rule to merge subsets of normalized VCFs
rule merge_vcfs:
    input:
        os.path.join(LISTS_DIR, "subset_vcfs.{idx}")
    output:
        os.path.join(MERGE_DIR, "merge.{idx}.vcf.gz")
    log:
        os.path.join(LOG_DIR, "merge.{idx}.log")
    threads: 8
    resources:
        mem_mb = get_mem_from_threads,
        time = '24:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        CONDA_ENV
    shell:
        """
        echo "Starting merge_vcfs {wildcards.idx} at: $(date)" >> {log}
        bcftools merge --threads {threads} -0 -l -m none -i BaseQRankSum:avg,ExcessHet:avg,FS:avg,MQ:avg,MQRankSum:avg,QD:avg,ReadPosRankSum:avg,SOR:avg,DP:avg,AF:sum,AS_BaseQRankSum:avg,AS_FS:avg,AS_MQ:avg,AS_MQRankSum:avg,AS_QD:avg,AS_ReadPosRankSum:avg,AS_SOR:avg,AS_UNIQ_ALT_READ_COUNT:avg,MLEAC:avg,MLEAF:avg,AN:sum,AC:sum {input} -Oz -o {output} &>> {log}
        bcftools index --threads {threads} {output} &>> {log}
        echo "Finished merge_vcfs {wildcards.idx} at: $(date)" >> {log}
        """

# Rule to compute MD5 checksum for merged subset VCFs
rule md5_merged_vcf:
    input:
        merged=os.path.join(MERGE_DIR, "merge.{idx}.vcf.gz")
    output:
        md5=os.path.join(LOG_DIR, "merge.{idx}.vcf.gz.md5")
    log:
        log_file=os.path.join(LOG_DIR, "md5_merge.{idx}.log")
    # Removed 'conda' directive as 'md5sum' is a standard utility
    shell:
        """
        echo "Computing MD5 for {input.merged} at $(date)" > {log_file}
        md5sum {input.merged} > {output.md5} 2>> {log_file}
        echo "Finished MD5 for {input.merged} at $(date)" >> {log_file}
        """

# Rule to create a list of all the merged VCF files
rule list_merged_vcfs:
    input:
        merged_files=expand(os.path.join(MERGE_DIR, "merge.{idx}.vcf.gz"), idx=range(needed_number_of_subsets))
    output:
        list_file=os.path.join(LISTS_DIR, "merge.txt")
    log:
        log_file=os.path.join(LOG_DIR, "list_merged_vcfs.log")
    run:
        from datetime import datetime

        with open(log.log_file, 'a') as log_file:
            log_file.write(f"Starting list_merged_vcfs at: {datetime.now()}\n")

        with open(output.list_file, 'w') as f:
            for file in input.merged_files:
                f.write(file + "\n")

        with open(log.log_file, 'a') as log_file:
            log_file.write(f"Finished list_merged_vcfs at: {datetime.now()}\n")

# Rule to perform the final merge of all batches into a single VCF file
rule final_merge:
    input:
        os.path.join(LISTS_DIR, "merge.txt")
    output:
        os.path.join(FINAL_DIR, FINAL_OUTPUT_NAME)
    log:
        os.path.join(LOG_DIR, "final_merge.log")
    threads: 8
    resources:
        mem_mb = get_mem_from_threads,
        time = '24:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        CONDA_ENV
    shell:
        """
        echo "Starting final_merge at: $(date)" >> {log}
        bcftools merge --threads {threads} -l -m none -i BaseQRankSum:avg,ExcessHet:avg,FS:avg,MQ:avg,MQRankSum:avg,QD:avg,ReadPosRankSum:avg,SOR:avg,DP:avg,AF:sum,AS_BaseQRankSum:avg,AS_FS:avg,AS_MQ:avg,AS_MQRankSum:avg,AS_QD:avg,AS_ReadPosRankSum:avg,AS_SOR:avg,AS_UNIQ_ALT_READ_COUNT:avg,MLEAC:avg,MLEAF:avg,AN:sum,AC:sum {input} -0 -Oz -o {output} &>> {log}
        bcftools index -t -f --threads {threads} {output} &>> {log}
        echo "Finished final_merge at: $(date)" >> {log}
        """

# Rule to compute MD5 checksum for the final merged VCF
rule md5_final_merge:
    input:
        merged=os.path.join(FINAL_DIR, FINAL_OUTPUT_NAME)
    output:
        md5=os.path.join(FINAL_DIR, f"{FINAL_OUTPUT_NAME}.md5")
    log:
        os.path.join(LOG_DIR, "md5_final_merge.log")
    # Removed 'conda' directive as 'md5sum' is a standard utility
    shell:
        """
        echo "Computing MD5 for {input.merged} at $(date)" > {log}
        md5sum {input.merged} > {output.md5} 2>> {log}
        echo "Finished MD5 for {input.merged} at $(date)" >> {log}
        """
# ----------------------------------------------------------------------------------- #
