#!/usr/bin/env python3
"""io_utils helpers that need no scientific stack.

Runpeng Luo (2026-08-05)

Dependencies:
  pytest; workflow/scripts/script_utils/io_utils.py of this repo. Its module imports
  pandas/numpy, so the tests skip where those are absent.

Usage:
  pytest tests/test_io_utils.py

Notes/References:
  Window build: workflow/scripts/build_window_bed.py
"""

import os
import sys

import pytest

_REPO = os.path.join(os.path.dirname(__file__), "..")
for _sub in ("config", "workflow/scripts/script_utils"):
    sys.path.insert(0, os.path.join(_REPO, _sub))
io_utils = pytest.importorskip("io_utils")

HG38 = os.path.join(_REPO, "resources", "data", "hg38.chrom.sizes")
MM10 = os.path.join(_REPO, "resources", "data", "mm10.chrom.sizes")


def _sizes(tmp_path, names):
    path = tmp_path / "genome_size.txt"
    path.write_text("".join(f"{n}\t1000\n" for n in names))
    return str(path)


def test_contigs_follow_requested_order(tmp_path):
    """Output order is the requested order, not the size file's."""
    sizes = _sizes(tmp_path, ["chr1", "chr2", "chrX"])
    found, absent = io_utils.select_contigs(sizes, ["X", 2, 1])
    assert found == ["chrX", "chr2", "chr1"]
    assert absent == []


def test_bare_contig_names_resolve(tmp_path):
    """An Ensembl-style size file has no chr prefix; requests still match."""
    sizes = _sizes(tmp_path, ["1", "2", "X"])
    found, absent = io_utils.select_contigs(sizes, [1, 2, "X"])
    assert found == ["1", "2", "X"]
    assert absent == []


def test_absent_chromosome_is_reported(tmp_path):
    """A requested chromosome with no contig is returned, not silently dropped."""
    sizes = _sizes(tmp_path, ["chr1", "chr2"])
    found, absent = io_utils.select_contigs(sizes, [1, 2, "X"])
    assert found == ["chr1", "chr2"]
    assert absent == ["X"]


def test_more_than_22_autosomes(tmp_path):
    """No human autosome cap: a 29-autosome genome keeps all of them."""
    sizes = _sizes(tmp_path, [f"chr{c}" for c in list(range(1, 30)) + ["X"]])
    found, absent = io_utils.select_contigs(sizes, list(range(1, 30)) + ["X"])
    assert len(found) == 30
    assert found[-3:] == ["chr28", "chr29", "chrX"]
    assert absent == []


def test_bundled_human_and_mouse_size_files():
    """The shipped size files resolve their own chromosome sets."""
    found, absent = io_utils.select_contigs(HG38, list(range(1, 23)) + ["X"])
    assert found[0] == "chr1" and found[-1] == "chrX" and absent == []
    found, absent = io_utils.select_contigs(MM10, list(range(1, 20)) + ["X"])
    assert len(found) == 20 and absent == []
    # mouse has 19 autosomes; asking for 20-22 reports them rather than dropping
    _, absent = io_utils.select_contigs(MM10, list(range(1, 23)) + ["X"])
    assert absent == ["20", "21", "22"]
