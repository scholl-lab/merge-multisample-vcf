import os
import glob
import functools

# ----------------------------------------------------------------------------------- #
# SCRIPT DESCRIPTION:
# This Snakemake script orchestrates a pipeline to manage large sets of VCF files.
# It splits a large number of VCF files into smaller subsets, merges them in batches, 
# and then merges all batches into a single VCF file. `bcftools` is utilized for 
# VCF processing tasks, which has a limitation of handling up to 1021 files at a time.
# Reference: https://shicheng-guo.github.io/bioinformatics/1923/02/28/bcftools-merge
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# ENVIRONMENT SETUP:
# Define the temporary directory using an environment variable typically set by
# the cluster scheduler.
SCRATCH_DIR = os.environ.get('TMPDIR')
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# CONFIGURATION:
# Extract the necessary parameters from the configuration file, including input 
# directory, batch size for VCF processing, and the Conda environment to be used.
configfile: "config.yaml"
INPUT_DIR = config["vcf_folder"]
VCFS_PER_BATCH = config["vcfs_per_batch"]
CONDA_ENV = config["conda_env"]
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# OUTPUT ORGANIZATION:
# Define result directories based on the provided output folder in the config.
# Using functools.partial to streamline the path joining process.
prefix_results = functools.partial(os.path.join, config['output_folder'])

# Define subdirectories for different stages of processing.
LISTS_DIR = prefix_results('lists')
MERGE_DIR = prefix_results('merged_vcfs')
FINAL_DIR = prefix_results('final')
LOG_DIR = prefix_results('logs')

# Ensure all output directories, including the log directory, exist.
for output_dir in [LISTS_DIR, MERGE_DIR, FINAL_DIR, LOG_DIR]:
    os.makedirs(output_dir, exist_ok=True)
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# HELPER FUNCTIONS:
# This function calculates the required memory based on the number of threads.
def get_mem_from_threads(wildcards, threads):
    return threads * 1000
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# SUBSET CALCULATION:
# Determine the number of subsets required based on the total number of VCFs
# and the desired batch size.
total_vcfs = len(glob.glob("{}/{}.vcf.gz".format(INPUT_DIR, '*')))
needed_number_of_subsets = -(-total_vcfs // VCFS_PER_BATCH)  # Ceiling division.
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# PIPELINE RULES:
# Define the rules for the Snakemake pipeline, guiding the workflow through various 
# stages including splitting, merging, and creating a list of merged VCF files.
# Each rule has clearly defined inputs, outputs, and log files to maintain a structured
# workflow.
# TODO:
# - Implement md5sum calculation for all files to verify data integrity.
# - Remove intermediate files to save storage space.
# - Ensure proper sequence of index and VCF file creation to maintain organization.

# Rule to initiate the pipeline with the final merged VCF file as the target output.
rule all:
    input:
        os.path.join(FINAL_DIR, "all_merged.vcf.gz")

# Rule to split the VCF files into smaller subsets for batch processing. It uses 
# Python's datetime module to log the start and end times of the process.
rule split_vcfs:
    input:
        vcf_files=glob.glob("{}/{}.vcf.gz".format(INPUT_DIR, '*'))
    output:
        expand(os.path.join(LISTS_DIR, "subset_vcfs.{idx}"), idx=range(needed_number_of_subsets))
    log:
        expand(os.path.join(LOG_DIR, "split_vcfs.{idx}.log"), idx=range(needed_number_of_subsets))
    run:
        from datetime import datetime
        
        # Write start timestamp to the log
        with open(log[0], 'a') as log_file:
            log_file.write(f"Starting split_vcfs at: {datetime.now()}\n")

        subsets = [input.vcf_files[i:i+VCFS_PER_BATCH] for i in range(0, len(input.vcf_files), VCFS_PER_BATCH)]
        for idx, subset in enumerate(subsets):
            with open(os.path.join(LISTS_DIR, "subset_vcfs." + str(idx)), 'w') as f:
                f.write("\n".join(subset))
                
        # Write end timestamp to the log
        with open(log[0], 'a') as log_file:
            log_file.write(f"Finished split_vcfs at: {datetime.now()}\n")

# Rule to merge the VCF files based on the subsets created in the previous rule, using
# bcftools for the merging process.
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
        echo "Starting merge_vcfs at: $(date)" >> {log}
        bcftools merge --threads {threads} -0 -l {input} -Oz -o {output} &>> {log}
        bcftools index --threads {threads} -t {output} &>> {log}
        echo "Finished merge_vcfs at: $(date)" >> {log}
        """

# Rule to create a list of all the merged VCF files, aiding in the final merging process.
rule list_merged_vcfs:
    input:
        merged_files=expand(os.path.join(MERGE_DIR, "merge.{idx}.vcf.gz"), idx=range(needed_number_of_subsets))
    output:
        os.path.join(LISTS_DIR, "merge.txt")
    log:
        os.path.join(LOG_DIR, "list_merged_vcfs.log")
    run:
        from datetime import datetime
        
        # Write start timestamp to the log
        with open(log[0], 'a') as log_file:
            log_file.write(f"Starting list_merged_vcfs at: {datetime.now()}\n")

        with open(output[0], 'w') as f:
            for file in input.merged_files:
                f.write(file + "\n")
                
        # Write end timestamp to the log
        with open(log[0], 'a') as log_file:
            log_file.write(f"Finished list_merged_vcfs at: {datetime.now()}\n")


# Rule to merge all the intermediate merged VCF files into a single final VCF file using 
# bcftools, marking the end of the pipeline.
rule final_merge:
    input:
        os.path.join(LISTS_DIR, "merge.txt")
    output:
        os.path.join(FINAL_DIR, "all_merged.vcf.gz")
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
        bcftools merge --threads {threads} -l {input} -0 -Oz -o {output} &>> {log}
        bcftools index --threads {threads} -t {output} &>> {log}
        echo "Finished final_merge at: $(date)" >> {log}
        """
# ----------------------------------------------------------------------------------- #
