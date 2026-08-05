#!/usr/bin/env python3
"""Chromosome-name helpers behind the resolution in parse_workflow.

Runpeng Luo (2026-08-05)

parse_workflow inlines the resolution itself from get_chr_sizes and
strip_chr_prefix; its behaviour is covered end to end by tests/test_chrstyle.py
and the absent-chromosome case in tests/test_dryrun.py.

Dependencies:
  pytest; workflow/scripts/script_utils of this repo, which imports pandas.

Usage:
  pytest tests/test_io_utils.py

Notes/References:
  Naming rules per file: .claude/chr-QA.md
"""

import os
import sys

import pytest

_REPO = os.path.join(os.path.dirname(__file__), "..")
for _sub in ("config", "workflow/scripts/script_utils"):
    sys.path.insert(0, os.path.join(_REPO, _sub))
io_utils = pytest.importorskip("io_utils")
utils = pytest.importorskip("utils")

HG38 = os.path.join(_REPO, "resources", "data", "hg38.chrom.sizes")
MM10 = os.path.join(_REPO, "resources", "data", "mm10.chrom.sizes")


@pytest.mark.parametrize(
    "name,core",
    [("chr1", "1"), ("1", "1"), ("chrX", "X"), ("X", "X"), ("CHR22", "22"), (22, "22")],
)
def test_strip_chr_prefix(name, core):
    """Either spelling, any case, str or int, reduces to the same core name."""
    assert utils.strip_chr_prefix(name) == core


def test_get_chr_sizes_keeps_file_names_and_order(tmp_path):
    """Names come back as the file spells them, in file order."""
    path = tmp_path / "genome_size.txt"
    path.write_text("chr2\t20\nchr1\t10\nchrX\t5\n")
    sizes = io_utils.get_chr_sizes(str(path))
    assert list(sizes) == ["chr2", "chr1", "chrX"]
    assert sizes["chr1"] == 10


def test_bare_contig_size_file(tmp_path):
    """An Ensembl-style file has no chr prefix; the cores still match."""
    path = tmp_path / "genome_size.txt"
    path.write_text("1\t10\n2\t20\nX\t5\n")
    sizes = io_utils.get_chr_sizes(str(path))
    assert list(sizes) == ["1", "2", "X"]
    assert [utils.strip_chr_prefix(c) for c in sizes] == ["1", "2", "X"]


def test_bundled_size_files_cover_their_chromosomes():
    """The shipped size files hold the chromosome set each genome has."""
    hg38 = io_utils.get_chr_sizes(HG38)
    mm10 = io_utils.get_chr_sizes(MM10)
    assert {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]} <= set(hg38)
    assert {f"chr{c}" for c in list(range(1, 20)) + ["X", "Y"]} <= set(mm10)
    assert "chr20" not in mm10


@pytest.mark.parametrize(
    "names,want",
    [
        (["22", "X"], ["chr22", "chrX"]),
        (["chr22", "chrX"], ["chr22", "chrX"]),
        ([], []),
    ],
)
def test_add_chr_prefix(names, want):
    """Readers that take a file in the genome's naming normalize on ingest."""
    pd = pytest.importorskip("pandas")
    assert utils.add_chr_prefix(pd.Series(names, dtype=str)).tolist() == want


@pytest.mark.parametrize(
    "first_line,input_nochr,ok",
    [
        ("chr22\t0\t100\t0.9\n", False, True),
        ("22\t0\t100\t0.9\n", True, True),
        ("chr22\t0\t100\t0.9\n", True, False),
        ("22\t0\t100\t0.9\n", False, False),
    ],
)
def test_require_matching_chr_style(tmp_path, first_line, input_nochr, ok):
    """A BED bedtools resolves against genome_size must use its naming."""
    bed = tmp_path / "x.bed"
    bed.write_text("track name=x\n" + first_line)
    if ok:
        utils.require_matching_chr_style(str(bed), input_nochr, "mappability_bed")
    else:
        with pytest.raises(ValueError, match="chr prefix"):
            utils.require_matching_chr_style(str(bed), input_nochr, "mappability_bed")


def test_require_matching_chr_style_reads_gzip(tmp_path):
    """Mappability tracks ship gzipped."""
    import gzip

    bed = tmp_path / "x.bed.gz"
    with gzip.open(bed, "wt") as fh:
        fh.write("22\t0\t100\t0.9\n")
    utils.require_matching_chr_style(str(bed), True, "mappability_bed")
    with pytest.raises(ValueError):
        utils.require_matching_chr_style(str(bed), False, "mappability_bed")
