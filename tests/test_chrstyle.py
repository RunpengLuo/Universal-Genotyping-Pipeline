#!/usr/bin/env python3
"""Chromosome naming style: region strings follow the genome-size file.

Runpeng Luo
Last update: 2026-08-06

The genome-size file declares how the raw inputs name their contigs. Region strings
handed to bcftools/mosdepth must use that spelling, or the tool silently matches
nothing. These build the DAG only.

Dependencies:
  pytest; snakemake on PATH.

Usage:
  pytest tests/test_chrstyle.py

Notes/References:
  Naming style: docs/reference.md, `genome_size`
"""

import os
import sys

import pytest

from conftest import dryrun

sys.path.insert(
    0,
    os.path.join(
        os.path.dirname(__file__), "..", "workflow", "scripts", "script_utils"
    ),
)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "config"))
io_utils = pytest.importorskip("io_utils")


def _bare_sizes(workspace):
    """A genome-size file with Ensembl-style bare contig names."""
    path = os.path.join(workspace["root"], "bare.chrom.sizes")
    with open(path, "w") as fh:
        fh.write("".join(f"{c}\t50818468\n" for c in list(range(1, 23)) + ["X"]))
    return path


def test_chr_prefixed_genome_keeps_chr_regions(workspace):
    """The bundled style: contigs are chr-prefixed, so regions are too."""
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=["remote_mode=stream"],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "input has chr-prefix=True" in proc.stdout
    assert "-r chr22" in proc.stdout
    assert "-c chr22" in proc.stdout
    assert "sed 's/^chr//'" not in proc.stdout


def test_bare_contig_genome_drops_the_prefix(workspace):
    """A bare-contig genome yields bare region strings, not chr-prefixed ones."""
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=[f"genome_size={_bare_sizes(workspace)}", "remote_mode=stream"],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "input has chr-prefix=False" in proc.stdout
    assert "-r 22" in proc.stdout and "-r chr22" not in proc.stdout
    assert "-c 22" in proc.stdout and "-c chr22" not in proc.stdout
    assert "sed 's/^chr//'" in proc.stdout, "windows.3col must drop the prefix"


@pytest.mark.parametrize("style", ["chr22", "22"])
def test_segment_bed_is_read_in_one_style(tmp_path, style):
    """Tiling and region_id assignment must see the same contig names.

    build_window_bed read the segment BED twice, raw for tiling and chr-normalized
    for region_id; on a bare-contig genome the two never matched and every window
    was dropped as off-segment.
    """
    seg = tmp_path / "segment.bed"
    seg.write_text(f"{style}\t0\t1000\t22p\t22p#0\n")
    regions = io_utils.read_BED(str(seg))
    contigs = [style]
    chroms = [c if c.lower().startswith("chr") else f"chr{c}" for c in contigs]
    tiled = regions[regions["#CHR"].isin(chroms)]
    assert len(tiled) == 1, f"{style}: segment BED row was filtered out"
    assert set(regions["#CHR"]) == {"chr22"}, "internal frames are chr-prefixed"
