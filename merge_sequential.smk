import os
import glob
import functools

# ----------------------------------------------------------------------------------- #
# SCRIPT DESCRIPTION:
# This Snakemake script processes VCF files by first splitting them into manageable
# subsets, then merging them, and finally creating a single merged VCF file.
# The script uses `bcftools` for VCF processing.
# based on the idea here: https://shicheng-guo.github.io/bioinformatics/1923/02/28/bcftools-merge
# bcftools seems to only handle up to 1021 files at a time, so we need to split the files
# ----------------------------------------------------------------------------------- #

# Read the configuration file
configfile: "config.yaml"

# ----------------------------------------------------------------------------------- #
# ENVIRONMENT SETUP:
# Define temporary directory using an environment variable 
# (usually set by the cluster scheduler)
SCRATCH_DIR = os.environ.get('TMPDIR')
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# CONFIGURATION:
# Extract user-defined input and output directories from the configuration.
# VCFS_PER_BATCH determines how many VCFs are processed in each batch.
# CONDA_ENV specifies the Conda environment to be used.
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
# TODO: fix logging (currently it overwrites amd the run rules dont work)
# TODO: add md5sum calculation for all input, intermediate and output files
# TODO: remove intermediate files
# TODO: make sure that the index is created after the vcf file

# Rule to find all .vcf.gz files and initiate the pipeline.
rule all:
    input:
        os.path.join(FINAL_DIR, "all_merged.vcf.gz")

# Rule to split the vcf files into subsets.
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

# Rule to merge the vcf files based on subsets.
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
        bcftools merge --threads {threads} -0 -l {input} -Oz -o {output} &> {log}
        bcftools index --threads {threads} -t {output} &> {log}
        echo "Finished merge_vcfs at: $(date)" >> {log}
        """

# Rule to create a list of all merged vcf files.
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


# Rule to merge all merged vcf files into one final file.
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
        bcftools merge --threads {threads} -l {input} -0 -Oz -o {output} &> {log}
        bcftools index --threads {threads} -t {output} &> {log}
        echo "Finished final_merge at: $(date)" >> {log}
        """
# ----------------------------------------------------------------------------------- #
