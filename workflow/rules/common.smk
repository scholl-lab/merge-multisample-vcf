"""Shared config extraction, directory setup, and input-helper functions."""

import os
import sys

from rules.helpers import (
    load_vcf_list as _load_vcf_list,
    batch_indices as _batch_indices,
    get_batch_samples as _get_batch_samples,
    mem_mb,
    n_batches,
)


# ---------------------------------------------------------------------------
# Config extraction
# ---------------------------------------------------------------------------

VCF_LIST_FILE = config["vcf_list_file"]
VCFS_PER_BATCH = int(config["vcfs_per_batch"])
OUTPUT_FOLDER = config["output_folder"]
REFERENCE_FASTA = config["reference_fasta"]

VCF_SUFFIX = config.get("vcf_suffix", ".vcf.gz")
FINAL_OUTPUT_NAME = config.get("final_output_name", "all_merged.vcf.gz")
FINAL_FILTER_LOGIC = config.get("final_filter_logic", "x")

# Validate config values at parse time for clear, early error messages.
if not os.path.isfile(VCF_LIST_FILE):
    sys.exit(
        f"ERROR: vcf_list_file={VCF_LIST_FILE!r} does not exist. " "Update config/config.yaml."
    )
if not os.path.isfile(REFERENCE_FASTA):
    sys.exit(
        f"ERROR: reference_fasta={REFERENCE_FASTA!r} does not exist. " "Update config/config.yaml."
    )
if VCFS_PER_BATCH <= 0:
    sys.exit(
        f"ERROR: vcfs_per_batch={VCFS_PER_BATCH!r} must be a positive integer. "
        "Update config/config.yaml."
    )
_VALID_FILTER_LOGIC = {"x", "+"}
if FINAL_FILTER_LOGIC not in _VALID_FILTER_LOGIC:
    sys.exit(
        f"ERROR: final_filter_logic={FINAL_FILTER_LOGIC!r} is not valid; "
        f"allowed values: {sorted(_VALID_FILTER_LOGIC)}. "
        "Update config/config.yaml."
    )

INFO_RULES = config.get(
    "info_rules",
    (
        "BaseQRankSum:avg,ExcessHet:avg,FS:avg,MQ:avg,MQRankSum:avg,QD:avg,"
        "ReadPosRankSum:avg,SOR:avg,DP:avg,AF:sum,AS_BaseQRankSum:avg,"
        "AS_FS:avg,AS_MQ:avg,AS_MQRankSum:avg,AS_QD:avg,AS_ReadPosRankSum:avg,"
        "AS_SOR:avg,AS_UNIQ_ALT_READ_COUNT:avg,MLEAC:avg,MLEAF:avg,AN:sum,AC:sum"
    ),
)

SCRATCH_DIR = os.environ.get("TMPDIR", "/tmp")


# ---------------------------------------------------------------------------
# Output directory layout
# ---------------------------------------------------------------------------

_p = lambda *parts: os.path.join(OUTPUT_FOLDER, *parts)  # noqa: E731

NORMALIZED_DIR = _p("normalized_vcfs")
MERGE_DIR = _p("merged_vcfs")
FINAL_DIR = _p("final")
LOG_DIR = _p("logs")

for _d in [NORMALIZED_DIR, MERGE_DIR, FINAL_DIR, LOG_DIR, os.path.join(LOG_DIR, "benchmarks")]:
    os.makedirs(_d, exist_ok=True)


# ---------------------------------------------------------------------------
# Sample loading & batch partitioning
# ---------------------------------------------------------------------------

SAMPLES, SAMPLE_TO_VCF = _load_vcf_list(VCF_LIST_FILE, VCF_SUFFIX)
TOTAL_VCFS = len(SAMPLES)
BATCH_INDICES = _batch_indices(TOTAL_VCFS, VCFS_PER_BATCH)

# Guard: bcftools merge accepts at most ~1021 positional file arguments.
# The final merge receives one file per batch, so validate here at parse time.
_n = n_batches(TOTAL_VCFS, VCFS_PER_BATCH)
if _n > 1000:
    sys.exit(
        f"ERROR: {_n} batches would exceed bcftools' ~1021-file limit at the "
        "final merge step. Increase vcfs_per_batch in config/config.yaml "
        f"(currently {VCFS_PER_BATCH}) so that ceil({TOTAL_VCFS} / vcfs_per_batch) <= 1000."
    )


# ---------------------------------------------------------------------------
# Input helper functions (used as lambdas in rule inputs)
# ---------------------------------------------------------------------------


def get_batch_vcfs(idx: int) -> list:
    """Normalized VCF paths for batch *idx*."""
    return [
        os.path.join(NORMALIZED_DIR, f"{s}.normalized.vcf.gz")
        for s in _get_batch_samples(idx, SAMPLES, VCFS_PER_BATCH)
    ]


def get_batch_tbis(idx: int) -> list:
    """Tabix index paths for batch *idx* (ensures Snakemake tracks them)."""
    return [
        os.path.join(NORMALIZED_DIR, f"{s}.normalized.vcf.gz.tbi")
        for s in _get_batch_samples(idx, SAMPLES, VCFS_PER_BATCH)
    ]


def get_final_outputs() -> list:
    """Targets for rule all: final VCF, its tabix index, and MD5."""
    return [
        os.path.join(FINAL_DIR, FINAL_OUTPUT_NAME),
        os.path.join(FINAL_DIR, f"{FINAL_OUTPUT_NAME}.tbi"),
        os.path.join(FINAL_DIR, f"{FINAL_OUTPUT_NAME}.md5"),
    ]
