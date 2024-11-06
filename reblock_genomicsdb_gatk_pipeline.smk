import os
import functools

# ----------------------------------------------------------------------------------- #
# SCRIPT DESCRIPTION:
# This Snakemake script orchestrates a pipeline to reblock input GVCF files and then
# use GenomicsDBImport and SelectVariants to create per-contig multi-sample VCFs.
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

# Ensure all output directories, including the log directory, exist.
for output_dir in [LISTS_DIR, REBLOCK_DIR, GENOMICS_DB_DIR, FINAL_DIR, LOG_DIR]:
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
# CONTIG IDENTIFICATION:
# Extract the list of contigs from the reference genome's .fai file.
# ----------------------------------------------------------------------------------- #

# Ensure that the reference genome index (.fai) exists.
fai_file = REFERENCE_GENOME + ".fai"
if not os.path.exists(fai_file):
    raise FileNotFoundError(f"Reference genome index file not found: {fai_file}")

# Extract contig names from the .fai file.
contigs = []
with open(fai_file, 'r') as fai:
    for line in fai:
        contig = line.split()[0]
        contigs.append(contig)
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# PIPELINE RULES:
# Define the rules for the Snakemake pipeline.
# ----------------------------------------------------------------------------------- #

# Rule to initiate the pipeline with all per-contig selected VCFs as the target outputs.
rule all:
    input:
        expand(os.path.join(FINAL_DIR, "selected_variants.{contig}.m.vcf.gz"), contig=contigs)
    # Note: No need to include 'all_samples.sample_map' as it's an intermediate step.

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
        mem_mb = 8000,
        time = '24:00:00',
        tmpdir = SCRATCH_DIR
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

# Rule to import GVCFs into separate GenomicsDB workspaces per contig.
rule genomicsdb_import:
    input:
        sample_map=os.path.join(LISTS_DIR, "all_samples.sample_map")
    output:
        db=os.path.join(GENOMICS_DB_DIR, "cohort_db_{contig}")
    log:
        os.path.join(LOG_DIR, "genomicsdb_import.{contig}.log")
    threads: 8
    resources:
        mem_mb = 40000,
        time = '240:00:00',
        tmpdir = SCRATCH_DIR
    params:
        contig=lambda wildcards: wildcards.contig
    conda:
        GATK_ENV
    shell:
        """
        set -e
        echo "Starting GenomicsDBImport for contig {wildcards.contig} at: $(date)" > {log}
        gatk --java-options "-Xmx{resources.mem_mb}m" GenomicsDBImport \
          --sample-name-map {input.sample_map} \
          --genomicsdb-workspace-path {output.db} \
          --tmp-dir {resources.tmpdir} \
          --reader-threads {threads} \
          --batch-size 50 \
          --overwrite-existing-genomicsdb-workspace \
          -L {wildcards.contig} &>> {log}
        echo "Finished GenomicsDBImport for contig {wildcards.contig} at: $(date)" >> {log}
        """

# Rule to select variants from each GenomicsDB workspace per contig.
rule select_variants:
    input:
        db=os.path.join(GENOMICS_DB_DIR, "cohort_db_{contig}")
    output:
        os.path.join(FINAL_DIR, "selected_variants.{contig}.m.vcf.gz")
    log:
        os.path.join(LOG_DIR, "select_variants.{contig}.log")
    threads: 8
    resources:
        mem_mb = 40000,
        time = '72:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        GATK_ENV
    shell:
        """
        set -e
        echo "Starting SelectVariants for contig {wildcards.contig} at: $(date)" > {log}
        gatk --java-options "-Xmx{resources.mem_mb}m" SelectVariants \
          -R {REFERENCE_GENOME} \
          -V gendb://{input.db} \
          -O {output} &>> {log}
        echo "Finished SelectVariants for contig {wildcards.contig} at: $(date)" >> {log}
        """
# ----------------------------------------------------------------------------------- #
