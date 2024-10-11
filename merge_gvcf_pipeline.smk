import os
import functools

# ----------------------------------------------------------------------------------- #
# SCRIPT DESCRIPTION:
# This Snakemake script orchestrates a pipeline to merge large sets of GVCF files.
# It splits a large number of GVCF files into smaller subsets, merges them in batches, 
# and then merges all batches into a single GVCF file. GATK's CombineGVCFs is utilized 
# for GVCF processing tasks, which has memory limitations for large numbers of samples.
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
# file, batch size for GVCF processing, and the Conda environment to be used.
configfile: "config.yaml"
INPUT_FILE = config["gvcf_list_file"]
GVCFS_PER_BATCH = config["gvcfs_per_batch"]
CONDA_ENV = config["conda_env"]
REFERENCE_GENOME = config["reference_genome"]
JAVA_MEM = config.get("java_memory", "20g")
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# OUTPUT ORGANIZATION:
# Define result directories based on the provided output folder in the config.
prefix_results = functools.partial(os.path.join, config['output_folder'])

# Define subdirectories for different stages of processing.
LISTS_DIR = prefix_results('lists')
MERGE_DIR = prefix_results('merged_gvcfs')
FINAL_DIR = prefix_results('final')
LOG_DIR = prefix_results('logs')

# Ensure all output directories, including the log directory, exist.
for output_dir in [LISTS_DIR, MERGE_DIR, FINAL_DIR, LOG_DIR]:
    os.makedirs(output_dir, exist_ok=True)
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# SUBSET CALCULATION:
# Determine the number of subsets required based on the total number of GVCFs
# and the desired batch size.
with open(INPUT_FILE, 'r') as f:
    gvcf_files = [line.strip() for line in f if line.strip()]

total_gvcfs = len(gvcf_files)
needed_number_of_subsets = -(-total_gvcfs // GVCFS_PER_BATCH)  # Ceiling division.
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# PIPELINE RULES:
# Define the rules for the Snakemake pipeline, guiding the workflow through various 
# stages including splitting, merging, and creating a list of merged GVCF files.
# ----------------------------------------------------------------------------------- #

# Rule to initiate the pipeline with the final merged GVCF file as the target output.
rule all:
    input:
        os.path.join(FINAL_DIR, "all_merged.g.vcf.gz")

# Rule to split the GVCF files into smaller subsets for batch processing.
rule split_input_list:
    input:
        INPUT_FILE
    output: expand(os.path.join(LISTS_DIR, "subset_gvcfs.{idx}.list"), idx=range(needed_number_of_subsets))
    log:
        expand(os.path.join(LOG_DIR, "split_gvcfs.{idx}.log"), idx=range(needed_number_of_subsets))
    run:
        subsets = [gvcf_files[i:i+GVCFS_PER_BATCH] for i in range(0, len(gvcf_files), GVCFS_PER_BATCH)]
        for idx, subset in enumerate(subsets):
            with open(os.path.join(LISTS_DIR, "subset_gvcfs." + str(idx) + ".list"), 'w') as f:
                f.write("\n".join(subset))

# Rule to merge the GVCF files based on the subsets created in the previous rule, using GATK's CombineGVCFs.
rule merge_gvcfs:
    input: os.path.join(LISTS_DIR, "subset_gvcfs.{idx}.list")
    output:
        os.path.join(MERGE_DIR, "merge.{idx}.g.vcf.gz")
    log:
        os.path.join(LOG_DIR, "merge_gvcfs.{idx}.log")
    threads: 8
    resources:
        mem_mb = lambda wildcards, threads: threads * 2500,
        time = '24:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        CONDA_ENV
    shell:
        """
        echo "Starting merge_gvcfs at: $(date)" >> {log}
        gatk --java-options "-Xmx{resources.mem_mb}m" CombineGVCFs \
          -R {REFERENCE_GENOME} \
          -O {output} \
          -G StandardAnnotation -G AS_StandardAnnotation \
          -V {input} &>> {log}
        echo "Finished merge_gvcfs at: $(date)" >> {log}
        """

# Rule to create a list of all the merged GVCF files, aiding in the final merging process.
rule list_merged_gvcfs:
    input:
        merged_files=expand(os.path.join(MERGE_DIR, "merge.{idx}.g.vcf.gz"), idx=range(needed_number_of_subsets))
    output: os.path.join(LISTS_DIR, "merge.list")
    log:
        os.path.join(LOG_DIR, "list_merged_gvcfs.log")
    run:
        with open(output[0], 'w') as f:
            for file in input.merged_files:
                f.write(file + "\n")

# Rule to merge all the intermediate merged GVCF files into a single final GVCF file using GATK's CombineGVCFs.
rule final_merge:
    input: os.path.join(LISTS_DIR, "merge.list")
    output:
        os.path.join(FINAL_DIR, "all_merged.g.vcf.gz")
    log:
        os.path.join(LOG_DIR, "final_merge.log")
    threads: 8
    resources:
        mem_mb = lambda wildcards, threads: threads * 2500,
        time = '24:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        CONDA_ENV
    shell:
        """
        echo "Starting final_merge at: $(date)" >> {log}
        gatk --java-options "-Xmx{resources.mem_mb}m" CombineGVCFs \
          -R {REFERENCE_GENOME} \
          -O {output} \
          -G StandardAnnotation -G AS_StandardAnnotation \
          -V {input} &>> {log}
        echo "Finished final_merge at: $(date)" >> {log}
        """
# ----------------------------------------------------------------------------------- #