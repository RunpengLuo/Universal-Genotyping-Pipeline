#!/usr/bin/env python3
"""Chromosome-name helpers behind the resolution in parse_workflow.

Runpeng Luo
Last update: 2026-08-06

parse_workflow inlines the resolution itself from read_chrom_sizes and
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
    sizes = io_utils.read_chrom_sizes(str(path))
    assert list(sizes) == ["chr2", "chr1", "chrX"]
    assert sizes["chr1"] == 10


def test_bare_contig_size_file(tmp_path):
    """An Ensembl-style file has no chr prefix; the cores still match."""
    path = tmp_path / "genome_size.txt"
    path.write_text("1\t10\n2\t20\nX\t5\n")
    sizes = io_utils.read_chrom_sizes(str(path))
    assert list(sizes) == ["1", "2", "X"]
    assert [utils.strip_chr_prefix(c) for c in sizes] == ["1", "2", "X"]


def test_bundled_size_files_cover_their_chromosomes():
    """The shipped size files hold the chromosome set each genome has."""
    hg38 = io_utils.read_chrom_sizes(HG38)
    mm10 = io_utils.read_chrom_sizes(MM10)
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
    "name,input_nochr,want",
    [
        ("22", False, "chr22"),
        ("chr22", False, "chr22"),
        ("22", True, "22"),
        ("chr22", True, "22"),
        ("chrX", True, "X"),
    ],
)
def test_match_chr_style(name, input_nochr, want):
    """Intervals handed to bedtools take the naming genome_size declares."""
    assert utils.match_chr_style(name, input_nochr) == want


def test_add_chr_prefix_is_per_row():
    """A GTF can start on a scaffold, so the decision cannot be made per file."""
    pd = pytest.importorskip("pandas")
    got = utils.add_chr_prefix(pd.Series(["GL000009.2", "chr1", "2"]))
    assert got.tolist() == ["chrGL000009.2", "chr1", "chr2"]
