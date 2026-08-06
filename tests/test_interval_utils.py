#!/usr/bin/env python3
"""Unit tests for the shared interval-assignment primitives.

Runpeng Luo (2026-08-06)

``interval_utils`` replaced five near-copies of the same "which interval contains
this?" loop, so these tests pin the behaviours the callers rely on: the vectorized
non-overlapping path, the scan fallback for overlapping references, the many-hit
join, and the largest-overlap rule.

Dependencies:
  pytest, numpy, pandas.

Usage:
  pytest tests/test_interval_utils.py

Notes/References:
  Callers: aggregation_utils (annotate_feature_type, adaptive_segmentation),
  rna_utils (feature_to_blocks), combine_counts_utils (_windows_to_bins,
  assign_snp_bounderies).
"""

import os
import sys

import pytest

_REPO = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(_REPO, "workflow", "scripts", "script_utils"))
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
iv = pytest.importorskip("interval_utils")


def _ref(rows, id_col="rid"):
    """Reference intervals from ``(chrom, start, end, id)`` tuples."""
    return pd.DataFrame(rows, columns=["#CHR", "START", "END", id_col])


DISJOINT = _ref([("chr1", 0, 10, "a"), ("chr1", 10, 20, "b"), ("chr2", 0, 10, "c")])


def test_pos_half_open_boundaries():
    """START is inside the interval, END is not."""
    qry = pd.DataFrame({"#CHR": ["chr1"] * 4, "POS0": [0, 9, 10, 20]})
    got = iv.assign_pos_to_range(qry, DISJOINT, ref_id="rid")
    assert got["rid"].tolist()[:3] == ["a", "a", "b"]
    assert pd.isna(got["rid"].iloc[3])


def test_pos_unmatched_chromosome_is_na():
    """A query on a contig the reference lacks stays unassigned."""
    qry = pd.DataFrame({"#CHR": ["chr9"], "POS0": [5]})
    got = iv.assign_pos_to_range(qry, DISJOINT, ref_id="rid")
    assert got["rid"].isna().all()


def test_pos_empty_reference_chromosome():
    """An empty reference assigns nothing and raises nothing."""
    got = iv.assign_pos_to_range(
        pd.DataFrame({"#CHR": ["chr1"], "POS0": [5]}),
        _ref([]),
        ref_id="rid",
    )
    assert got["rid"].isna().all()


def test_pos_overlapping_reference_takes_first_by_start():
    """With overlapping intervals the scan fallback keeps the earliest-starting hit."""
    ref = _ref([("chr1", 0, 100, "outer"), ("chr1", 40, 60, "inner")])
    qry = pd.DataFrame({"#CHR": ["chr1", "chr1"], "POS0": [50, 10]})
    got = iv.assign_pos_to_range(qry, ref, ref_id="rid")
    assert got["rid"].tolist() == ["outer", "outer"]


def test_pos_unsorted_reference():
    """Reference order does not matter; intervals are start-sorted internally."""
    ref = _ref([("chr1", 10, 20, "b"), ("chr1", 0, 10, "a")])
    qry = pd.DataFrame({"#CHR": ["chr1", "chr1"], "POS0": [5, 15]})
    got = iv.assign_pos_to_range(qry, ref, ref_id="rid")
    assert got["rid"].tolist() == ["a", "b"]


def test_interval_takes_largest_overlap():
    """A query interval straddling two references picks the larger overlap."""
    ref = _ref([("chr1", 0, 10, "a"), ("chr1", 10, 30, "b")])
    qry = pd.DataFrame({"#CHR": ["chr1"], "START": [8], "END": [20]})
    got = iv.assign_interval_to_range(qry, ref, "rid")
    assert got["rid"].tolist() == ["b"]


def test_interval_touching_edge_does_not_overlap():
    """Half-open intervals that only touch share no bases."""
    ref = _ref([("chr1", 10, 20, "a")])
    qry = pd.DataFrame({"#CHR": ["chr1", "chr1"], "START": [0, 20], "END": [10, 30]})
    got = iv.assign_interval_to_range(qry, ref, "rid")
    assert got["rid"].isna().all()


def test_interval_does_not_mutate_input():
    """assign_interval_to_range copies; assign_pos_to_range does not."""
    qry = pd.DataFrame({"#CHR": ["chr1"], "START": [0], "END": [5]})
    iv.assign_interval_to_range(qry, _ref([("chr1", 0, 10, "a")]), "rid")
    assert "rid" not in qry.columns


def test_all_features_joins_every_overlap():
    """Nested genes are all reported, joined in tier order; misses get the default."""
    ref = _ref([("chr1", 0, 100, "G1"), ("chr1", 40, 60, "G2")], id_col="gene_id")
    snps = pd.DataFrame({"#CHR": ["chr1", "chr1", "chr1"], "POS0": [50, 10, 500]})
    got = iv.assign_all_features(snps, ref, id_col="gene_id")
    assert set(got.iloc[0].split(";")) == {"G1", "G2"}
    assert got.iloc[1] == "G1"
    assert got.iloc[2] == "intergenic"


def test_all_features_empty_reference():
    """No annotation at all means every SNP is intergenic."""
    snps = pd.DataFrame({"#CHR": ["chr1"], "POS0": [5]})
    got = iv.assign_all_features(snps, _ref([], id_col="gene_id"), id_col="gene_id")
    assert got.tolist() == ["intergenic"]


def test_overlaps_any_range_matches_assignment():
    """The boolean variant agrees with assign_pos_to_range on nested intervals."""
    ref = _ref([("chr1", 0, 100, "a"), ("chr1", 40, 60, "b"), ("chr2", 0, 10, "c")])
    qry = pd.DataFrame(
        {"#CHR": ["chr1"] * 3 + ["chr2", "chr9"], "POS0": [50, 150, 0, 5, 5]}
    )
    got = iv.overlaps_any_range(qry, ref)
    assert got.tolist() == [True, False, True, True, False]
    assigned = iv.assign_pos_to_range(qry.copy(), ref, ref_id="rid")["rid"].notna()
    assert got.tolist() == assigned.tolist()


def test_overlaps_any_range_empty_reference():
    """No reference intervals means nothing overlaps."""
    qry = pd.DataFrame({"#CHR": ["chr1"], "POS0": [5]})
    assert not iv.overlaps_any_range(qry, _ref([])).any()


@pytest.mark.parametrize("n_ref", [1, 2, 5])
def test_interval_tiers_are_non_overlapping(n_ref):
    """Every tier is internally disjoint, which is what the fast path assumes."""
    starts = np.arange(n_ref) * 5
    ends = starts + 12
    for tier in iv._interval_tiers(starts, ends):
        s, e = starts[tier], ends[tier]
        order = np.argsort(s)
        s, e = s[order], e[order]
        assert np.all(s[1:] >= e[:-1])
