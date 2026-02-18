"""Unit tests for workflow/rules/helpers.py."""

import sys
from pathlib import Path

import pytest

# Make helpers importable without a Snakemake installation.
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "rules"))

from helpers import (
    _strip_vcf_suffix,
    batch_indices,
    get_batch_samples,
    load_vcf_list,
    mem_mb,
    n_batches,
)

# ---------------------------------------------------------------------------
# load_vcf_list
# ---------------------------------------------------------------------------


def test_load_vcf_list_basic(tmp_path):
    lst = tmp_path / "samples.list"
    lst.write_text("/data/sampleA.vcf.gz\n/data/sampleB.vcf.gz\n")
    samples, s2p = load_vcf_list(str(lst))
    assert samples == ["sampleA", "sampleB"]
    assert s2p["sampleA"] == "/data/sampleA.vcf.gz"
    assert s2p["sampleB"] == "/data/sampleB.vcf.gz"


def test_load_vcf_list_blank_and_comments_skipped(tmp_path):
    lst = tmp_path / "samples.list"
    lst.write_text("# comment\n\n  # indented comment\n/data/sampleA.vcf.gz\n")
    samples, _ = load_vcf_list(str(lst))
    assert samples == ["sampleA"]


def test_load_vcf_list_empty_raises(tmp_path):
    lst = tmp_path / "empty.list"
    lst.write_text("# only a comment\n\n")
    with pytest.raises(ValueError, match="no valid entries"):
        load_vcf_list(str(lst))


def test_load_vcf_list_duplicate_raises(tmp_path):
    lst = tmp_path / "dup.list"
    lst.write_text("/a/sampleA.vcf.gz\n/b/sampleA.vcf.gz\n")
    with pytest.raises(ValueError, match="Duplicate sample name"):
        load_vcf_list(str(lst))


def test_load_vcf_list_custom_suffix(tmp_path):
    lst = tmp_path / "gvcf.list"
    lst.write_text("/data/sampleA.g.vcf.gz\n/data/sampleB.g.vcf.gz\n")
    samples, _ = load_vcf_list(str(lst), suffix=".g.vcf.gz")
    assert samples == ["sampleA", "sampleB"]


# ---------------------------------------------------------------------------
# _strip_vcf_suffix
# ---------------------------------------------------------------------------


def test_strip_vcf_suffix_gz():
    assert _strip_vcf_suffix("sampleA.vcf.gz", ".vcf.gz") == "sampleA"


def test_strip_vcf_suffix_plain_vcf():
    assert _strip_vcf_suffix("sampleA.vcf", ".vcf.gz") == "sampleA"


def test_strip_vcf_suffix_gvcf():
    assert _strip_vcf_suffix("sampleA.g.vcf.gz", ".g.vcf.gz") == "sampleA"


def test_strip_vcf_suffix_no_match():
    # When the suffix doesn't match, name is returned as-is (minus .vcf if present).
    result = _strip_vcf_suffix("sampleA.bam", ".vcf.gz")
    assert result == "sampleA.bam"


# ---------------------------------------------------------------------------
# n_batches / batch_indices
# ---------------------------------------------------------------------------


def test_n_batches_exact():
    assert n_batches(100, 50) == 2


def test_n_batches_ceiling():
    assert n_batches(101, 50) == 3


def test_n_batches_single():
    assert n_batches(1, 100) == 1


def test_n_batches_zero_total():
    assert n_batches(0, 50) == 0


def test_n_batches_invalid_batch_size():
    with pytest.raises(ValueError, match="positive"):
        n_batches(10, 0)


def test_batch_indices():
    assert batch_indices(10, 3) == [0, 1, 2, 3]
    assert batch_indices(9, 3) == [0, 1, 2]


# ---------------------------------------------------------------------------
# get_batch_samples
# ---------------------------------------------------------------------------


def test_get_batch_samples_first():
    samples = list("ABCDE")
    assert get_batch_samples(0, samples, 2) == ["A", "B"]


def test_get_batch_samples_last():
    samples = list("ABCDE")
    assert get_batch_samples(2, samples, 2) == ["E"]


def test_get_batch_samples_middle():
    samples = list("ABCDE")
    assert get_batch_samples(1, samples, 2) == ["C", "D"]


# ---------------------------------------------------------------------------
# mem_mb
# ---------------------------------------------------------------------------


def test_mem_mb_attempt_1():
    fn = mem_mb(4_000)
    assert fn(None, 1) == 4_000


def test_mem_mb_doubles_on_retry():
    fn = mem_mb(4_000)
    assert fn(None, 2) == 8_000
    assert fn(None, 3) == 16_000


def test_mem_mb_capped_at_max():
    fn = mem_mb(64_000, max_mb=100_000)
    assert fn(None, 2) == 100_000  # would be 128 000, capped at 100 000


def test_mem_mb_returns_int():
    fn = mem_mb(4_000)
    assert isinstance(fn(None, 1), int)
