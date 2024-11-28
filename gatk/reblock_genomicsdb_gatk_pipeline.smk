import os
import functools

# ----------------------------------------------------------------------------------- #
# SCRIPT DESCRIPTION:
# This Snakemake script orchestrates a pipeline to reblock input GVCF files and then
# uses GenomicsDBImport and GenotypeGVCFs to create multi-sample VCFs over equal intervals,
# followed by merging these VCFs into a single consolidated VCF using Picard's MergeVcfs.
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# ENVIRONMENT SETUP:
# Define the temporary directory using an environment variable typically set by
# the cluster scheduler.
SCRATCH_DIR = os.environ.get('TMPDIR', '/tmp')
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# CONFIGURATION:
# Extract the necessary parameters from the configuration file, including input 
# file, the Conda environments to be used, and reference genome.
configfile: "config.yaml"
INPUT_FILE = config["gvcf_list_file"]
GATK_ENV = config["gatk_env"]
REFERENCE_GENOME = config["reference_genome"]
SCATTER_COUNT = int(config.get("scatter_count", 400))  # Number of intervals to split into
FINAL_VCF_NAME = config.get("final_vcf_name", "merged_variants.vcf.gz")  # Final VCF filename
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# OUTPUT ORGANIZATION:
# Define result directories based on the provided output folder in the config.
prefix_results = functools.partial(os.path.join, config['output_folder'])

# Define subdirectories for different stages of processing.
LISTS_DIR = prefix_results('lists')
REBLOCK_DIR = prefix_results('reblocked_gvcfs')
GENOMICS_DB_DIR = prefix_results('genomicsdb')
FINAL_DIR = prefix_results('final')
LOG_DIR = prefix_results('logs')
INTERVALS_DIR = prefix_results('intervals')

# Ensure all output directories, including the log directory, exist.
for output_dir in [LISTS_DIR, REBLOCK_DIR, GENOMICS_DB_DIR, FINAL_DIR, LOG_DIR, INTERVALS_DIR]:
    os.makedirs(output_dir, exist_ok=True)
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# SAMPLE PROCESSING:
# Read the list of GVCF files and extract sample names.
# ----------------------------------------------------------------------------------- #

# Read the list of GVCF files from the input file.
with open(INPUT_FILE, 'r') as f:
    gvcf_list = [line.strip() for line in f if line.strip()]

# Extract sample names from the basenames of the input files.
samples = [os.path.basename(f).replace('.g.vcf.gz', '') for f in gvcf_list]

# Create a mapping from sample names to GVCF file paths.
sample_to_gvcf = dict(zip(samples, gvcf_list))
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# INTERVAL IDENTIFICATION:
# Define interval IDs based on SCATTER_COUNT.
# ----------------------------------------------------------------------------------- #

# Generate interval IDs from 1 to SCATTER_COUNT.
interval_ids = [str(i) for i in range(1, SCATTER_COUNT + 1)]
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# PIPELINE RULES:
# Define the rules for the Snakemake pipeline.
# ----------------------------------------------------------------------------------- #

# Rule to initiate the pipeline with the final merged VCF as the target output.
rule all:
    input:
        os.path.join(FINAL_DIR, FINAL_VCF_NAME)

# Rule to generate intervals using ScatterIntervalsByNs.
rule scatter_intervals_by_ns:
    input:
        reference=REFERENCE_GENOME
    output:
        intervals_list=os.path.join(INTERVALS_DIR, "intervals.interval_list")
    log:
        os.path.join(LOG_DIR, "scatter_intervals_by_ns.log")
    conda:
        GATK_ENV
    shell:
        """
        set -e
        echo "Starting ScatterIntervalsByNs at: $(date)" > {log}
        gatk ScatterIntervalsByNs \
            -R {input.reference} \
            -O {output.intervals_list} &>> {log}
        echo "Finished ScatterIntervalsByNs at: $(date)" >> {log}
        """

# Rule to split intervals into equal parts using SplitIntervals.
rule split_intervals:
    input:
        intervals_list=os.path.join(INTERVALS_DIR, "intervals.interval_list"),
        reference=REFERENCE_GENOME
    output:
        interval_files=expand(os.path.join(INTERVALS_DIR, "scattered.interval_{interval_id}_of_{scatter_count}.interval_list"),
                              interval_id=interval_ids,
                              scatter_count=str(SCATTER_COUNT))
    params:
        scatter_count=SCATTER_COUNT
    log:
        os.path.join(LOG_DIR, "split_intervals.log")
    conda:
        GATK_ENV
    shell:
        """
        set -e
        echo "Starting SplitIntervals at: $(date)" > {log}
        gatk SplitIntervals \
            -R {input.reference} \
            -L {input.intervals_list} \
            --scatter-count {params.scatter_count} \
            -O {INTERVALS_DIR} &>> {log}
        echo "Finished SplitIntervals at: $(date)" >> {log}
        """

# Rule to reblock the GVCF files and create sample map entries.
rule reblock_gvcfs:
    input:
        gvcf=lambda wildcards: sample_to_gvcf[wildcards.sample]
    output:
        reblocked_gvcf=os.path.join(REBLOCK_DIR, "{sample}.rb.g.vcf.gz"),
        sample_map_entry=os.path.join(LISTS_DIR, "{sample}.sample_map")
    log:
        os.path.join(LOG_DIR, "reblock_gvcf.{sample}.log")
    threads: 2
    resources:
        mem_mb=8000,
        time='24:00:00',
        tmpdir=SCRATCH_DIR
    conda:
        GATK_ENV
    shell:
        """
        set -e
        echo "Starting ReblockGVCF for sample {wildcards.sample} at: $(date)" > {log}
        gatk --java-options "-Xmx{resources.mem_mb}m" ReblockGVCF \
          -R {REFERENCE_GENOME} \
          -V {input.gvcf} \
          -O {output.reblocked_gvcf} &>> {log}
        
        # Extract sample name from the reblocked GVCF
        sample_name=$(zgrep "^#CHROM" {output.reblocked_gvcf} | head -n1 | cut -f10)
        echo -e "$sample_name\t{output.reblocked_gvcf}" > {output.sample_map_entry}
        
        echo "Finished ReblockGVCF for sample {wildcards.sample} at: $(date)" >> {log}
        """

# Rule to create a unified sample map for all reblocked GVCF files.
rule create_sample_map:
    input:
        sample_map_entries=expand(os.path.join(LISTS_DIR, "{sample}.sample_map"), sample=samples)
    output:
        os.path.join(LISTS_DIR, "all_samples.sample_map")
    log:
        os.path.join(LOG_DIR, "create_sample_map.log")
    run:
        with open(output[0], 'w') as outfile:
            for sample_map_file in input.sample_map_entries:
                with open(sample_map_file, 'r') as infile:
                    outfile.write(infile.read())

# Rule to import GVCFs into GenomicsDB workspaces per interval.
rule genomicsdb_import:
    input:
        sample_map=os.path.join(LISTS_DIR, "all_samples.sample_map"),
        interval_list=lambda wildcards: os.path.join(INTERVALS_DIR, f'scattered.interval_{wildcards.interval_id}_of_{SCATTER_COUNT}.interval_list')
    output:
        db=directory(os.path.join(GENOMICS_DB_DIR, "cohort_db_{interval_id}"))
    log:
        os.path.join(LOG_DIR, "genomicsdb_import.{interval_id}.log")
    threads: 2
    resources:
        mem_mb=20000,
        time='240:00:00',
        tmpdir=SCRATCH_DIR
    conda:
        GATK_ENV
    shell:
        """
        set -e
        echo "Starting GenomicsDBImport for interval {wildcards.interval_id} at: $(date)" > {log}
        gatk --java-options "-Xmx{resources.mem_mb}m -XX:ParallelGCThreads=2" GenomicsDBImport \
          --sample-name-map {input.sample_map} \
          --genomicsdb-workspace-path {output.db} \
          --tmp-dir {resources.tmpdir} \
          --reader-threads {threads} \
          --batch-size 50 \
          --overwrite-existing-genomicsdb-workspace \
          -L {input.interval_list} &>> {log}
        echo "Finished GenomicsDBImport for interval {wildcards.interval_id} at: $(date)" >> {log}
        """

# Rule to perform joint genotyping using GenotypeGVCFs per interval.
rule genotype_gvcfs:
    input:
        db=os.path.join(GENOMICS_DB_DIR, "cohort_db_{interval_id}"),
        interval_list=lambda wildcards: os.path.join(INTERVALS_DIR, f'scattered.interval_{wildcards.interval_id}_of_{SCATTER_COUNT}.interval_list')
    output:
        os.path.join(FINAL_DIR, "genotyped_variants.interval_{interval_id}.vcf.gz")
    log:
        os.path.join(LOG_DIR, "genotype_gvcfs.{interval_id}.log")
    threads: 2
    resources:
        mem_mb=40000,
        time='72:00:00',
        tmpdir=SCRATCH_DIR
    conda:
        GATK_ENV
    shell:
        """
        set -e
        echo "Starting GenotypeGVCFs for interval {wildcards.interval_id} at: $(date)" > {log}
        gatk --java-options "-Xmx{resources.mem_mb}m -XX:ParallelGCThreads=2" GenotypeGVCFs \
          -R {REFERENCE_GENOME} \
          --only-output-calls-starting-in-intervals \
          --use-new-qual-calculator \
          -V gendb://{input.db} \
          -O {output} \
          -L {input.interval_list} &>> {log}
        echo "Finished GenotypeGVCFs for interval {wildcards.interval_id} at: $(date)" >> {log}
        """

# Rule to create a list of all VCFs for merging.
rule create_vcf_list:
    input:
        genotype_vcfs=expand(os.path.join(FINAL_DIR, "genotyped_variants.interval_{interval_id}.vcf.gz"),
                             interval_id=interval_ids)
    output:
        list_file=os.path.join(FINAL_DIR, "merged_vcfs.list")
    log:
        os.path.join(LOG_DIR, "create_vcf_list.log")
    run:
        with open(output.list_file, 'w') as outfile:
            for vcf in input.genotype_vcfs:
                outfile.write(f"I={vcf}\n")

# Rule to merge all VCFs into a single VCF using Picard's MergeVcfs.
rule merge_vcfs:
    input:
        list_file=os.path.join(FINAL_DIR, "merged_vcfs.list")
    output:
        merged_vcf=os.path.join(FINAL_DIR, FINAL_VCF_NAME)
    log:
        os.path.join(LOG_DIR, "merge_vcfs.log")
    threads: 2
    resources:
        mem_mb=20000,
        time='72:00:00',
        tmpdir=SCRATCH_DIR
    conda:
        GATK_ENV
    shell:
        """
        set -e
        echo "Starting MergeVcfs at: $(date)" > {log}
        picard MergeVcfs \
            $(cat {input.list_file}) \
            O={output.merged_vcf} &>> {log}
        echo "Finished MergeVcfs at: $(date)" >> {log}
        """
# ----------------------------------------------------------------------------------- #
