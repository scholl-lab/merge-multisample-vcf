import os
import functools

# ----------------------------------------------------------------------------------- #
# SCRIPT DESCRIPTION:
# This Snakemake script orchestrates a pipeline to reblock input GVCF files and then
# use GenomicsDBImport and SelectVariants to create a final merged GVCF.
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
# file, the Conda environment to be used.
configfile: "config.yaml"
INPUT_FILE = config["gvcf_list_file"]
GATK_ENV = config["gatk_env"]
BCFTOOLS_ENV = config["bcftools_env"]
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

# Read the list of GVCF files from the input file.
with open(INPUT_FILE, 'r') as f:
    gvcf_list = [line.strip() for line in f if line.strip()]

# Create a mapping from indices to GVCF file paths.
gvcf_dict = { idx: filepath for idx, filepath in enumerate(gvcf_list) }

# Total number of GVCF files.
NUM_GVCFS = len(gvcf_list)
# ----------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------- #
# PIPELINE RULES:
# Define the rules for the Snakemake pipeline, guiding the workflow through various 
# stages including reblocking, GenomicsDBImport, and creating a final merged GVCF.
# ----------------------------------------------------------------------------------- #

# Rule to initiate the pipeline with the final merged GVCF file as the target output.
rule all:
    input:
        os.path.join(FINAL_DIR, "selected_variants.g.vcf.gz")

# Rule to reblock the GVCF files.
rule reblock_gvcfs:
    input:
        gvcf=lambda wildcards: gvcf_dict[int(wildcards.idx)]
    output:
        reblocked_gvcf=os.path.join(REBLOCK_DIR, "{idx}.rb.g.vcf.gz")
    log:
        os.path.join(LOG_DIR, "reblock_gvcf.{idx}.log")
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
        echo "Starting ReblockGVCF at: $(date)" > {log}
        gatk --java-options "-Xmx{resources.mem_mb}m" ReblockGVCF \
          -R {REFERENCE_GENOME} \
          -V {input.gvcf} \
          -O {output.reblocked_gvcf} &>> {log}
        echo "Finished ReblockGVCF at: $(date)" >> {log}
        """

# Rule to get the sample name from the reblocked GVCF and write a sample map file.
rule get_sample_name:
    input:
        reblocked_gvcf=os.path.join(REBLOCK_DIR, "{idx}.rb.g.vcf.gz")
    output:
        sample_map=os.path.join(LISTS_DIR, "{idx}.sample_map")
    log:
        os.path.join(LOG_DIR, "get_sample_name.{idx}.log")
    conda:
        BCFTOOLS_ENV
    shell:
        """
        set -e
        echo "Starting get_sample_name at: $(date)" > {log}
        sample_name=$(bcftools query -l {input.reblocked_gvcf})
        echo -e "$sample_name\t{input.reblocked_gvcf}" > {output.sample_map}
        echo "Finished get_sample_name at: $(date)" >> {log}
        """

# Rule to create a sample map for all reblocked GVCF files.
rule create_sample_map:
    input:
        sample_maps=expand(os.path.join(LISTS_DIR, "{idx}.sample_map"), idx=range(NUM_GVCFS))
    output:
        os.path.join(LISTS_DIR, "all_samples.sample_map")
    log:
        os.path.join(LOG_DIR, "create_sample_map.log")
    run:
        with open(output[0], 'w') as outfile:
            for sample_map_file in input.sample_maps:
                with open(sample_map_file, 'r') as infile:
                    outfile.write(infile.read())

# Rule to import GVCFs into GenomicsDB.
rule genomicsdb_import:
    input:
        sample_map=os.path.join(LISTS_DIR, "all_samples.sample_map")
    output:
        db=os.path.join(GENOMICS_DB_DIR, "cohort_db")
    log:
        os.path.join(LOG_DIR, "genomicsdb_import.log")
    threads: 8
    resources:
        mem_mb = 20000,
        time = '48:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        GATK_ENV
    shell:
        """
        set -e
        echo "Starting GenomicsDBImport at: $(date)" > {log}
        gatk --java-options "-Xmx{resources.mem_mb}m" GenomicsDBImport \
          --sample-name-map {input.sample_map} \
          --genomicsdb-workspace-path {output.db} \
          --tmp-dir {resources.tmpdir} \
          --reader-threads {threads} \
          --batch-size 50 \
          --overwrite-existing-genomicsdb-workspace \
          $(awk '{{print "-L " $1}}' {REFERENCE_GENOME}.fai) &>> {log}
        echo "Finished GenomicsDBImport at: $(date)" >> {log}
        """

# Rule to select variants from the imported GenomicsDB.
rule select_variants:
    input:
        db=os.path.join(GENOMICS_DB_DIR, "cohort_db")
    output:
        os.path.join(FINAL_DIR, "selected_variants.g.vcf.gz")
    log:
        os.path.join(LOG_DIR, "select_variants.log")
    threads: 8
    resources:
        mem_mb = 20000,
        time = '24:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        GATK_ENV
    shell:
        """
        set -e
        echo "Starting SelectVariants at: $(date)" > {log}
        gatk --java-options "-Xmx{resources.mem_mb}m" SelectVariants \
          -R {REFERENCE_GENOME} \
          -V gendb://{input.db} \
          -O {output} &>> {log}
        echo "Finished SelectVariants at: $(date)" >> {log}
        """

# ----------------------------------------------------------------------------------- #
