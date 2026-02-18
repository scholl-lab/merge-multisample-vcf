#!/usr/bin/env python3
"""Generate config/config.yaml (and optionally a VCF list file) for merge-multisample-vcf.

Scans a VCF directory or accepts a pre-existing VCF list file, auto-detects a
reference genome FASTA from common project and HPC shared locations, and writes
a ready-to-edit pipeline config.

Two modes:
  Interactive:  python scripts/generate_config.py              (guided wizard)
  Flags:        python scripts/generate_config.py --vcf-folder /path/to/vcfs

Usage:
    python scripts/generate_config.py
    python scripts/generate_config.py --vcf-folder /data/cohort/vcfs
    python scripts/generate_config.py --vcf-list /data/cohort/samples.list
    python scripts/generate_config.py --vcf-folder /data/cohort/vcfs --ref /ref/GRCh38.fa
    python scripts/generate_config.py --vcf-folder /data/cohort/vcfs --dry-run
    python scripts/generate_config.py --vcf-folder /data/cohort/vcfs \\
        --output-folder results/cohort_A --vcfs-per-batch 200
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

GENOME_EXTENSIONS = ("*.fna", "*.fa", "*.fasta")

# Search directories for reference FASTA (relative to project root)
REF_SEARCH_DIRS = [
    "resources/ref/GRCh38",
    "resources/ref/GRCh37",
    "resources/ref",
    "analysis/ref/GRCh38",
    "analysis/ref",
    "../resources/ref/GRCh38",
    "../resources/ref",
]

# Well-known shared locations on BIH and Charité HPC
SHARED_REF_DIRS = [
    "/data/cephfs-1/work/groups/scholl/shared/ref/GRCh38",
    "/data/cephfs-1/work/groups/scholl/shared/ref/GRCh37",
    "/data/cephfs-1/work/groups/scholl/shared/ref",
]

# Default INFO rules (GATK joint-calling conventions)
DEFAULT_INFO_RULES = (
    "BaseQRankSum:avg,ExcessHet:avg,FS:avg,MQ:avg,MQRankSum:avg,QD:avg,"
    "ReadPosRankSum:avg,SOR:avg,DP:avg,AF:sum,AS_BaseQRankSum:avg,"
    "AS_FS:avg,AS_MQ:avg,AS_MQRankSum:avg,AS_QD:avg,AS_ReadPosRankSum:avg,"
    "AS_SOR:avg,AS_UNIQ_ALT_READ_COUNT:avg,MLEAC:avg,MLEAF:avg,AN:sum,AC:sum"
)


# ---------------------------------------------------------------------------
# Section A: VCF discovery and list handling
# ---------------------------------------------------------------------------


def discover_vcf_files(vcf_folder: str | Path) -> list[dict[str, Any]]:
    """Scan a directory for *.vcf.gz files.

    Returns a list of dicts with keys: path, basename, size_bytes.
    """
    folder = Path(vcf_folder)
    if not folder.is_dir():
        sys.exit(f"Error: VCF directory does not exist: {folder}")

    entries = []
    for vcf in sorted(folder.glob("*.vcf.gz")):
        if vcf.is_file():
            entries.append(
                {
                    "path": str(vcf.resolve()),
                    "basename": vcf.name.removesuffix(".vcf.gz"),
                    "size_bytes": vcf.stat().st_size,
                }
            )
    return entries


def load_vcf_list(list_path: str | Path) -> list[str]:
    """Read a VCF list file (one path per line, blank lines and # ignored).

    Returns a list of path strings. Warns about missing files.
    """
    path = Path(list_path)
    if not path.is_file():
        sys.exit(f"Error: VCF list file does not exist: {path}")

    paths: list[str] = []
    missing = 0
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if not Path(stripped).is_file():
                print(
                    f"  Warning: line {lineno}: file not found: {stripped}",
                    file=sys.stderr,
                )
                missing += 1
            paths.append(stripped)

    if missing:
        print(
            f"Warning: {missing} of {len(paths)} VCF path(s) in the list do not exist.",
            file=sys.stderr,
        )
    return paths


def write_vcf_list(
    vcf_paths: list[str],
    output_path: Path,
    *,
    dry_run: bool = False,
    overwrite: bool = False,
) -> None:
    """Write VCF paths to a list file, one per line.

    Does nothing (prints a message) when dry_run is True.
    """
    if dry_run:
        print(f"[dry-run] Would write {len(vcf_paths)} path(s) to: {output_path}")
        return

    if output_path.is_file() and not overwrite:
        response = input(f"File already exists: {output_path}. Overwrite? [y/N] ")
        if response.strip().lower() not in ("y", "yes"):
            print(f"Skipped writing {output_path}")
            return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        for p in vcf_paths:
            fh.write(p + "\n")
    print(f"Written: {output_path}  ({len(vcf_paths)} VCF paths)")


# ---------------------------------------------------------------------------
# Section B: Reference FASTA discovery
# ---------------------------------------------------------------------------


def _find_genome_fastas(search_dir: Path) -> list[dict[str, Any]]:
    """Find FASTA files in a directory.

    Returns list of dicts with keys: path, has_fai, has_dict, name.
    """
    results: list[dict[str, Any]] = []
    if not search_dir.is_dir():
        return results

    for ext in GENOME_EXTENSIONS:
        for fasta in sorted(search_dir.glob(ext)):
            if not fasta.is_file() or fasta.name.startswith("."):
                continue
            has_fai = Path(str(fasta) + ".fai").is_file()
            has_dict = fasta.with_suffix(".dict").is_file() or Path(
                str(fasta) + ".dict"
            ).is_file()
            results.append(
                {
                    "path": str(fasta),
                    "has_fai": has_fai,
                    "has_dict": has_dict,
                    "name": fasta.name,
                }
            )
    return results


def discover_reference_fasta(
    ref_dir: Path | None = None,
    project_root: Path | None = None,
) -> dict[str, Any]:
    """Scan for a reference genome FASTA.

    Search order: explicit ref_dir -> relative project paths -> HPC shared dirs.

    Returns dict with keys: path (str or None), search_log (list[str]).
    """
    search_log: list[str] = []
    if project_root is None:
        project_root = Path.cwd()

    search_dirs: list[Path] = []
    if ref_dir is not None:
        search_dirs.append(ref_dir)
    for rel in REF_SEARCH_DIRS:
        search_dirs.append(project_root / rel)
    for shared in SHARED_REF_DIRS:
        search_dirs.append(Path(shared))

    for sd in search_dirs:
        fastas = _find_genome_fastas(sd)
        if not fastas:
            if sd.is_dir():
                search_log.append(f"No FASTA in: {sd}")
            else:
                search_log.append(f"Dir not found: {sd}")
            continue

        # Prefer indexed (fai) over unindexed
        indexed = [f for f in fastas if f["has_fai"]]
        best = indexed[0] if indexed else fastas[0]
        status = []
        if best["has_fai"]:
            status.append(".fai")
        if best["has_dict"]:
            status.append(".dict")
        label = f"  [{', '.join(status)}]" if status else "  [no index]"
        search_log.append(f"Found: {best['path']}{label}")
        return {"path": str(Path(best["path"]).resolve()), "search_log": search_log}

    search_log.append("No reference genome found in any search directory.")
    return {"path": None, "search_log": search_log}


# ---------------------------------------------------------------------------
# Section C: Config template generation
# ---------------------------------------------------------------------------


def _format_info_rules_yaml(rules: str) -> str:
    """Format the info_rules string as a YAML block scalar."""
    # Wrap the comma-separated string at ~80 chars for readability
    chunks: list[str] = []
    current = ""
    for part in rules.split(","):
        part = part.strip()
        candidate = f"{current},{part}" if current else part
        if len(candidate) > 78 and current:
            chunks.append(current + ",")
            current = part
        else:
            current = candidate
    if current:
        chunks.append(current)

    lines = "\n  ".join(chunks)
    return f">-\n  {lines}"


def generate_config_yaml(
    vcf_list_file: str,
    reference_fasta: str,
    output_folder: str = "results",
    vcfs_per_batch: int = 100,
    vcf_suffix: str = ".vcf.gz",
    final_output_name: str = "all_merged.vcf.gz",
    final_filter_logic: str = "x",
    info_rules: str = DEFAULT_INFO_RULES,
) -> str:
    """Build the config.yaml content string from pipeline parameters."""
    info_rules_yaml = _format_info_rules_yaml(info_rules)

    return f"""\
# merge-multisample-vcf — pipeline configuration
# Generated by: python scripts/generate_config.py
# Review all paths before running the pipeline.

# ---------------------------------------------------------------------------
# Required
# ---------------------------------------------------------------------------

# Path to a plain-text file listing input VCFs, one absolute path per line.
# Blank lines and lines starting with '#' are ignored.
vcf_list_file: "{vcf_list_file}"

# Number of VCFs merged in each intermediate batch.
# Keep below ~1000 (bcftools merge positional-argument limit).
vcfs_per_batch: {vcfs_per_batch}

# Root directory for all pipeline outputs.
output_folder: "{output_folder}"

# Reference FASTA for bcftools norm (bgzipped+indexed or plain FASTA).
reference_fasta: "{reference_fasta}"

# ---------------------------------------------------------------------------
# Optional — defaults shown
# ---------------------------------------------------------------------------

# File extension stripped when deriving sample names from VCF basenames.
vcf_suffix: "{vcf_suffix}"

# Filename for the final merged VCF placed in {{output_folder}}/final/.
final_output_name: "{final_output_name}"

# FILTER field merging logic for bcftools merge -F:
#   "x"  — PASS if *any* input sample passes  (default, more permissive)
#   "+"  — PASS only if *all* input samples pass (more strict)
final_filter_logic: "{final_filter_logic}"

# Comma-separated INFO field aggregation rules for bcftools merge -i.
info_rules: {info_rules_yaml}
"""


def write_config(
    content: str,
    output_path: Path,
    *,
    dry_run: bool = False,
    overwrite: bool = False,
) -> None:
    """Write the config YAML to disk.

    Prints the content and returns when dry_run is True.
    """
    if dry_run:
        print(f"\n--- config.yaml ({output_path}) ---")
        print(content)
        return

    if output_path.is_file() and not overwrite:
        response = input(f"File already exists: {output_path}. Overwrite? [y/N] ")
        if response.strip().lower() not in ("y", "yes"):
            print(f"Skipped writing {output_path}")
            return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        fh.write(content)
    print(f"Written: {output_path}")


# ---------------------------------------------------------------------------
# Section D: Prompt helpers
# ---------------------------------------------------------------------------


def _prompt(prompt: str, default: str = "") -> str:
    """Prompt for text with an optional default."""
    if default:
        result = input(f"{prompt} [{default}]: ").strip()
        return result if result else default
    return input(f"{prompt}: ").strip()


def _prompt_yn(prompt: str, default: bool = True) -> bool:
    """Prompt for yes/no with a default."""
    suffix = "[Y/n]" if default else "[y/N]"
    result = input(f"{prompt} {suffix}: ").strip().lower()
    if not result:
        return default
    return result in ("y", "yes")


def _prompt_path(prompt: str, default: str = "", *, must_exist: bool = True) -> str:
    """Prompt for a filesystem path, re-prompting if validation fails."""
    while True:
        raw = _prompt(prompt, default)
        if not raw:
            if not must_exist:
                return ""
            print("  Path cannot be empty. Please try again.")
            continue
        p = Path(raw).expanduser()
        if must_exist and not p.exists():
            print(f"  Path does not exist: {p}")
            if not _prompt_yn("  Try again?", default=True):
                return str(p)
            continue
        return str(p)


def _prompt_int(prompt: str, default: int, min_val: int = 1, max_val: int = 1000) -> int:
    """Prompt for an integer within a range."""
    while True:
        raw = _prompt(prompt, str(default))
        try:
            value = int(raw)
        except ValueError:
            print(f"  Please enter a number between {min_val} and {max_val}.")
            continue
        if not min_val <= value <= max_val:
            print(f"  Value must be between {min_val} and {max_val}.")
            continue
        return value


def _prompt_choice(prompt: str, choices: list[str], default: str) -> str:
    """Prompt with numbered choices."""
    for i, choice in enumerate(choices, 1):
        marker = " <--" if choice == default else ""
        print(f"    [{i}] {choice}{marker}")
    while True:
        raw = _prompt(f"  {prompt}", default=default)
        if raw.isdigit():
            idx = int(raw) - 1
            if 0 <= idx < len(choices):
                return choices[idx]
        if raw in choices:
            return raw
        print(f"  Invalid — enter 1–{len(choices)} or a value from the list.")


def format_size(size_bytes: int) -> str:
    """Convert bytes to a human-readable string."""
    value = float(size_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024:
            return f"{value:.1f} {unit}"
        value /= 1024
    return f"{value:.1f} PB"


# ---------------------------------------------------------------------------
# Section E: Interactive wizard
# ---------------------------------------------------------------------------


def _interactive_wizard() -> None:
    """Guided wizard for generating pipeline config files."""
    print()
    print("=" * 60)
    print("  merge-multisample-vcf — Config Generator")
    print("=" * 60)
    print()

    project_root = Path.cwd()

    # ---- Step 1: VCF input ----
    print("Step 1: VCF input")
    print("  [1] Scan a directory for *.vcf.gz files")
    print("  [2] Use an existing VCF list file")
    choice = _prompt("  Choice", default="1")

    vcf_paths: list[str] = []
    vcf_list_file: str = ""
    write_list = False

    if choice == "2":
        # Existing list file
        list_path_str = _prompt_path("  Path to VCF list file")
        vcf_paths = load_vcf_list(list_path_str)
        vcf_list_file = str(Path(list_path_str).resolve())
        print(f"\n  Loaded {len(vcf_paths)} path(s) from {vcf_list_file}")
    else:
        # Scan directory
        vcf_folder_str = _prompt_path("  VCF directory to scan")
        entries = discover_vcf_files(vcf_folder_str)
        if not entries:
            print(f"\n  No *.vcf.gz files found in: {vcf_folder_str}")
            if not _prompt_yn("  Continue anyway?", default=False):
                print("Aborted.")
                return
        else:
            print(f"\n  Found {len(entries)} VCF file(s):")
            for e in entries[:10]:
                print(f"    {e['basename']}.vcf.gz  ({format_size(e['size_bytes'])})")
            if len(entries) > 10:
                print(f"    ... and {len(entries) - 10} more")
        print()

        vcf_paths = [e["path"] for e in entries]
        default_list = "input/input.list"
        list_out = _prompt("  Save VCF paths to list file", default=default_list)
        vcf_list_file = str(Path(list_out).resolve())
        write_list = True

    print()

    # ---- Step 2: Output folder ----
    print("Step 2: Pipeline output folder")
    output_folder = _prompt("  Output folder", default="results")
    print()

    # ---- Step 3: Batch size ----
    print("Step 3: Batch size")
    print("  How many VCFs to merge in each intermediate batch.")
    print("  Keep ≤ 1000 (bcftools hard limit). Larger = more RAM per batch.")
    vcfs_per_batch = _prompt_int("  vcfs_per_batch", default=100, min_val=1, max_val=1000)
    print()

    # ---- Step 4: Reference FASTA ----
    print("Step 4: Reference FASTA")
    print("  Scanning for reference genome...")
    ref_result = discover_reference_fasta(project_root=project_root)
    for line in ref_result["search_log"]:
        print(f"  {line}")

    reference_fasta: str
    if ref_result["path"]:
        print(f"\n  Auto-detected: {ref_result['path']}")
        if _prompt_yn("  Use this reference?", default=True):
            reference_fasta = ref_result["path"]
        else:
            reference_fasta = _prompt_path("  Path to reference FASTA", must_exist=False)
    else:
        print("  No reference genome found automatically.")
        reference_fasta = _prompt_path("  Path to reference FASTA (or press Enter to skip)", must_exist=False)
        if not reference_fasta:
            reference_fasta = "EDIT_ME: /path/to/reference.fa"
    print()

    # ---- Step 5: Optional settings ----
    print("Step 5: Optional settings (press Enter to accept defaults)")
    vcf_suffix = _prompt("  vcf_suffix", default=".vcf.gz")
    final_output_name = _prompt("  final_output_name", default="all_merged.vcf.gz")
    print("  final_filter_logic:")
    final_filter_logic = _prompt_choice(
        "filter logic", choices=["x", "+"], default="x"
    )
    print()

    # ---- Step 6: Summary + write ----
    print("=" * 60)
    print("  Summary")
    print("=" * 60)
    print(f"  VCF count:          {len(vcf_paths)}")
    if write_list:
        print(f"  VCF list file:      {vcf_list_file}")
    print(f"  Output folder:      {output_folder}")
    print(f"  Batch size:         {vcfs_per_batch}")
    print(f"  Reference FASTA:    {reference_fasta}")
    print(f"  VCF suffix:         {vcf_suffix}")
    print(f"  Final output name:  {final_output_name}")
    print(f"  Filter logic:       {final_filter_logic}")
    print()

    dry_run = not _prompt_yn("Write files now?", default=True)
    overwrite = False
    if not dry_run:
        overwrite = _prompt_yn("Overwrite existing files without asking?", default=False)
    print()

    # Write VCF list
    if write_list:
        write_vcf_list(
            vcf_paths,
            Path(vcf_list_file),
            dry_run=dry_run,
            overwrite=overwrite,
        )

    # Generate and write config
    content = generate_config_yaml(
        vcf_list_file=vcf_list_file,
        reference_fasta=reference_fasta,
        output_folder=output_folder,
        vcfs_per_batch=vcfs_per_batch,
        vcf_suffix=vcf_suffix,
        final_output_name=final_output_name,
        final_filter_logic=final_filter_logic,
    )
    write_config(content, Path("config/config.yaml"), dry_run=dry_run, overwrite=overwrite)

    if "EDIT_ME" in content:
        print("\nDone! Search for EDIT_ME in config/config.yaml to fill in missing paths.")
    else:
        print("\nDone! Review config/config.yaml before running the pipeline.")


# ---------------------------------------------------------------------------
# Section F: CLI entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """Main entry point: parse arguments, discover files, generate output."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate config/config.yaml for merge-multisample-vcf.\n"
            "Run without arguments for an interactive guided wizard."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  %(prog)s                                          # interactive wizard
  %(prog)s --vcf-folder /data/cohort/vcfs           # scan directory
  %(prog)s --vcf-list /data/cohort/vcfs.list        # existing list
  %(prog)s --vcf-folder /data/cohort/vcfs --dry-run
  %(prog)s --vcf-folder /data/cohort/vcfs --ref /ref/GRCh38.fa --output-folder results/A5297
""",
    )

    # VCF input (mutually exclusive)
    vcf_group = parser.add_mutually_exclusive_group()
    vcf_group.add_argument(
        "--vcf-folder",
        metavar="DIR",
        help="Scan this directory for *.vcf.gz files and generate a VCF list file.",
    )
    vcf_group.add_argument(
        "--vcf-list",
        metavar="FILE",
        help="Use a pre-existing VCF list file (one path per line).",
    )

    parser.add_argument(
        "--list-output",
        default="input/input.list",
        metavar="FILE",
        help="Path for the generated VCF list file (default: input/input.list). "
        "Only used with --vcf-folder.",
    )
    parser.add_argument(
        "--ref",
        metavar="FASTA",
        help="Reference FASTA for bcftools norm. Auto-detected if not provided.",
    )
    parser.add_argument(
        "--ref-dir",
        metavar="DIR",
        help="Directory to search for reference FASTA (supplemental to auto-detection).",
    )
    parser.add_argument(
        "--output-folder",
        default="results",
        metavar="DIR",
        help="Pipeline output folder (default: results).",
    )
    parser.add_argument(
        "--vcfs-per-batch",
        type=int,
        default=100,
        metavar="N",
        help="VCFs per intermediate batch, keep ≤ 1000 (default: 100).",
    )
    parser.add_argument(
        "--vcf-suffix",
        default=".vcf.gz",
        metavar="EXT",
        help="Suffix stripped to derive sample names (default: .vcf.gz).",
    )
    parser.add_argument(
        "--final-output-name",
        default="all_merged.vcf.gz",
        metavar="NAME",
        help="Final merged VCF filename (default: all_merged.vcf.gz).",
    )
    parser.add_argument(
        "--filter-logic",
        default="x",
        choices=["x", "+"],
        dest="filter_logic",
        help="bcftools merge -F: 'x' = PASS if any sample passes (default), "
        "'+' = all samples must pass.",
    )
    parser.add_argument(
        "--config-output",
        default="config/config.yaml",
        metavar="FILE",
        help="Output config.yaml path (default: config/config.yaml).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be generated without writing any files.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files without asking.",
    )

    args = parser.parse_args()

    # No VCF input specified → interactive wizard
    if args.vcf_folder is None and args.vcf_list is None:
        _interactive_wizard()
        return

    # ---- CLI mode ----
    project_root = Path.cwd()

    # Resolve VCF paths
    vcf_paths: list[str]
    vcf_list_file: str

    if args.vcf_folder:
        entries = discover_vcf_files(args.vcf_folder)
        if not entries:
            sys.exit(f"Error: No *.vcf.gz files found in: {args.vcf_folder}")
        vcf_paths = [e["path"] for e in entries]
        vcf_list_file = str(Path(args.list_output).resolve())

        print(f"Found {len(entries)} VCF file(s) in {args.vcf_folder}")
        for e in entries[:5]:
            print(f"  {e['basename']}.vcf.gz  ({format_size(e['size_bytes'])})")
        if len(entries) > 5:
            print(f"  ... and {len(entries) - 5} more")

        write_vcf_list(
            vcf_paths,
            Path(args.list_output),
            dry_run=args.dry_run,
            overwrite=args.overwrite,
        )
    else:
        vcf_paths = load_vcf_list(args.vcf_list)
        vcf_list_file = str(Path(args.vcf_list).resolve())
        print(f"Loaded {len(vcf_paths)} path(s) from {vcf_list_file}")

    # Resolve reference FASTA
    reference_fasta: str
    if args.ref:
        reference_fasta = str(Path(args.ref).resolve())
        print(f"Reference: {reference_fasta}")
    else:
        ref_dir = Path(args.ref_dir).resolve() if args.ref_dir else None
        ref_result = discover_reference_fasta(ref_dir=ref_dir, project_root=project_root)
        for line in ref_result["search_log"]:
            print(f"  {line}")
        if ref_result["path"]:
            reference_fasta = ref_result["path"]
            print(f"Reference: {reference_fasta}")
        else:
            reference_fasta = "EDIT_ME: /path/to/reference.fa"
            print("Warning: No reference FASTA found — using EDIT_ME placeholder.")

    # Validate batch size
    if not 1 <= args.vcfs_per_batch <= 1000:
        sys.exit(f"Error: --vcfs-per-batch must be between 1 and 1000, got {args.vcfs_per_batch}")

    # Generate config
    content = generate_config_yaml(
        vcf_list_file=vcf_list_file,
        reference_fasta=reference_fasta,
        output_folder=args.output_folder,
        vcfs_per_batch=args.vcfs_per_batch,
        vcf_suffix=args.vcf_suffix,
        final_output_name=args.final_output_name,
        final_filter_logic=args.filter_logic,
    )

    write_config(content, Path(args.config_output), dry_run=args.dry_run, overwrite=args.overwrite)

    if "EDIT_ME" in content:
        print("\nDone! Search for EDIT_ME in config/config.yaml to fill in missing paths.")
    else:
        print("\nDone! Review config/config.yaml before running the pipeline.")


if __name__ == "__main__":
    main()
