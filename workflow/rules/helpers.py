"""Pure Python helpers for the merge-multisample-vcf pipeline.

No Snakemake imports — functions are independently unit-testable.
"""

from __future__ import annotations

import os
from typing import Any

# ---------------------------------------------------------------------------
# VCF list loading
# ---------------------------------------------------------------------------


def load_vcf_list(
    path: str,
    suffix: str = ".vcf.gz",
) -> tuple[list[str], dict[str, str]]:
    """Read a VCF list file and return ``(samples, sample_to_path)``.

    Sample names are derived by stripping *suffix* from each file basename.
    Blank lines and lines starting with ``#`` are ignored.

    Raises ``ValueError`` if any two entries resolve to the same sample name.
    """
    paths: list[str] = []
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            paths.append(stripped)

    samples: list[str] = []
    sample_to_path: dict[str, str] = {}

    for p in paths:
        name = _strip_vcf_suffix(os.path.basename(p), suffix)
        if name in sample_to_path:
            raise ValueError(
                f"Duplicate sample name '{name}' from files:\n  {sample_to_path[name]}\n  {p}"
            )
        samples.append(name)
        sample_to_path[name] = p

    return samples, sample_to_path


def _strip_vcf_suffix(basename: str, suffix: str) -> str:
    """Remove *suffix* (and a trailing '.vcf' if present) from *basename*."""
    name = basename.removesuffix(suffix)
    # Handle plain .vcf files that were listed without compression
    if name.endswith(".vcf"):
        name = name[:-4]
    return name


# ---------------------------------------------------------------------------
# Batch arithmetic
# ---------------------------------------------------------------------------


def n_batches(total: int, batch_size: int) -> int:
    """Return the number of batches needed (ceiling division)."""
    if batch_size <= 0:
        raise ValueError(f"batch_size must be positive, got {batch_size}")
    return -(-total // batch_size)


def batch_indices(total: int, batch_size: int) -> list[int]:
    """Return zero-based batch index list ``[0, 1, ..., n_batches-1]``."""
    return list(range(n_batches(total, batch_size)))


def get_batch_samples(
    idx: int,
    samples: list[str],
    batch_size: int,
) -> list[str]:
    """Return the slice of sample names belonging to batch *idx*."""
    start = idx * batch_size
    return samples[start : start + batch_size]


# ---------------------------------------------------------------------------
# Resource helpers
# ---------------------------------------------------------------------------


def mem_mb(base_mb: int, max_mb: int = 128_000):
    """Return a Snakemake-compatible callable that doubles memory on each retry.

    Usage::

        resources:
            mem_mb=mem_mb(4_000),
    """

    def _mem(wildcards: Any, attempt: int) -> int:
        return int(min(base_mb * (2 ** (attempt - 1)), max_mb))

    return _mem
