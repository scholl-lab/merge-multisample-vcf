"""Run both real merge shell blocks against synthetic VCFs (requires bcftools)."""

import hashlib
import re
import shlex
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

pytestmark = pytest.mark.skipif(
    not shutil.which("bcftools") or not shutil.which("bash"),
    reason="requires bcftools >=1.18 with +fill-tags, and bash",
)


@pytest.mark.parametrize("rule", ["merge_vcfs", "final_merge"])
@pytest.mark.parametrize(
    "tag_number", ["1", ".", None], ids=["scalar-tag", "vector-tag", "no-tag"]
)
def test_merge_recomputes_missingness_from_genotypes(tmp_path, rule, tag_number):
    """Inherited F_MISSING must not shadow the genotype-based expression."""
    header = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1,length=100>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    )
    if tag_number:
        header += (
            f'##INFO=<ID=F_MISSING,Number={tag_number},Type=Float,Description="Old missingness">\n'
        )
    info = "F_MISSING=0.9" if tag_number else "."
    inputs = []
    for sample, genotype in [("sample_a", "./."), ("sample_b", "0/1")]:
        vcf = tmp_path / f"{sample}.vcf"
        vcf.write_text(
            header
            + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
            + f"1\t10\t.\tA\tC\t.\tPASS\t{info}\tGT\t{genotype}\n"
            + f"1\t20\t.\tG\tT\t.\tPASS\t{info}\tGT\t0/1\n",
            encoding="utf-8",
        )
        compressed = vcf.with_suffix(".vcf.gz")
        subprocess.run(
            ["bcftools", "view", "-Oz", "-o", str(compressed), "-W=tbi", str(vcf)],
            check=True,
        )
        inputs.append(shlex.quote(str(compressed)))

    # Execute the repository's shell block so a regression in either rule is caught.
    source = (Path(__file__).parents[1] / "workflow/rules/merge.smk").read_text(encoding="utf-8")
    block = re.search(rf"rule {rule}:.*?shell:\s*r\"\"\"(.*?)\"\"\"", source, re.S)
    assert block is not None
    output = tmp_path / "merged.vcf.gz"
    checksum = tmp_path / "merged.vcf.gz.md5"
    log = tmp_path / "merge.log"
    command = (
        block[1]
        .replace("{params.info_rules:q}", "'-'")
        .format(
            threads=1,
            params=SimpleNamespace(filter_logic="x"),
            input=SimpleNamespace(vcfs=" ".join(inputs)),
            output=SimpleNamespace(vcf=shlex.quote(str(output)), md5=shlex.quote(str(checksum))),
            wildcards=SimpleNamespace(idx=0),
            log=shlex.quote(str(log)),
        )
    )
    result = subprocess.run(["bash", "-c", command], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr + log.read_text(encoding="utf-8")
    rows = subprocess.check_output(
        ["bcftools", "query", "-f", "%POS\t%INFO/F_MISSING[\t%GT]\n", str(output)],
        text=True,
    ).splitlines()
    assert rows == ["10\t0.5\t./.\t0/1", "20\t0\t0/1\t0/1"]
    assert output.with_suffix(".gz.tbi").is_file()
    assert checksum.read_text().split()[0] == hashlib.md5(output.read_bytes()).hexdigest()
