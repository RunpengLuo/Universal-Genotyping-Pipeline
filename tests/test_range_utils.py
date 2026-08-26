#!/usr/bin/env python3
"""Unit tests for the shared range-assignment primitives.

Last update: 2026-08-06

Covers:
- boundaries: the half-open rule, including the last base of a range
- contract: na_idx, dropna, fillna, out_col, and no mutation of qry
- rules: max_overlap, midpoint and contained, including split segments
- invariants: a 1-based pos_col and degenerate ranges are rejected
- trim_range_by_range: splits, full cover, column carry, and row order
"""

import os
import sys

import pytest

_REPO = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(_REPO, "workflow", "scripts", "script_utils"))
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
iv = pytest.importorskip("range_utils")


def _ref(rows, id_col="rid"):
    """Reference ranges from ``(chrom, start, end, id)`` tuples."""
    return pd.DataFrame(rows, columns=["#CHR", "START", "END", id_col])


DISJOINT = _ref([("chr1", 0, 10, "a"), ("chr1", 10, 20, "b"), ("chr2", 0, 10, "c")])


def test_pos_half_open_boundaries():
    """START is inside the range, END is not."""
    qry = pd.DataFrame({"#CHR": ["chr1"] * 4, "POS0": [0, 9, 10, 20]})
    got, _ = iv.assign_pos_to_range(qry, DISJOINT, ref_id="rid")
    assert got["rid"].tolist()[:3] == ["a", "a", "b"]
    assert pd.isna(got["rid"].iloc[3])


def test_pos_unmatched_chromosome_is_na():
    """A query on a contig the reference lacks stays unassigned."""
    qry = pd.DataFrame({"#CHR": ["chr9"], "POS0": [5]})
    got, _ = iv.assign_pos_to_range(qry, DISJOINT, ref_id="rid")
    assert got["rid"].isna().all()


def test_pos_empty_reference_chromosome():
    """An empty reference assigns nothing and raises nothing."""
    got, _ = iv.assign_pos_to_range(
        pd.DataFrame({"#CHR": ["chr1"], "POS0": [5]}),
        _ref([]),
        ref_id="rid",
    )
    assert got["rid"].isna().all()


def test_pos_overlapping_reference_takes_first_by_start():
    """With overlapping ranges the scan fallback keeps the earliest-starting hit."""
    ref = _ref([("chr1", 0, 100, "outer"), ("chr1", 40, 60, "inner")])
    qry = pd.DataFrame({"#CHR": ["chr1", "chr1"], "POS0": [50, 10]})
    got, _ = iv.assign_pos_to_range(qry, ref, ref_id="rid")
    assert got["rid"].tolist() == ["outer", "outer"]


def test_pos_unsorted_reference():
    """Reference order does not matter; ranges are start-sorted internally."""
    ref = _ref([("chr1", 10, 20, "b"), ("chr1", 0, 10, "a")])
    qry = pd.DataFrame({"#CHR": ["chr1", "chr1"], "POS0": [5, 15]})
    got, _ = iv.assign_pos_to_range(qry, ref, ref_id="rid")
    assert got["rid"].tolist() == ["a", "b"]


def test_range_takes_largest_overlap():
    """A query range straddling two references picks the larger overlap."""
    ref = _ref([("chr1", 0, 10, "a"), ("chr1", 10, 30, "b")])
    qry = pd.DataFrame({"#CHR": ["chr1"], "START": [8], "END": [20]})
    got, _ = iv.assign_range_to_range(qry, ref, "rid")
    assert got["rid"].tolist() == ["b"]


def test_range_touching_edge_does_not_overlap():
    """Half-open ranges that only touch share no bases."""
    ref = _ref([("chr1", 10, 20, "a")])
    qry = pd.DataFrame({"#CHR": ["chr1", "chr1"], "START": [0, 20], "END": [10, 30]})
    got, _ = iv.assign_range_to_range(qry, ref, "rid")
    assert got["rid"].isna().all()


def test_assign_never_mutates_the_input():
    """Both primitives copy: the caller must bind the returned frame."""
    ref = _ref([("chr1", 0, 10, "a")])
    rng = pd.DataFrame({"#CHR": ["chr1"], "START": [0], "END": [5]})
    iv.assign_range_to_range(rng, ref, "rid")
    assert "rid" not in rng.columns
    pos = pd.DataFrame({"#CHR": ["chr1"], "POS0": [5]})
    iv.assign_pos_to_range(pos, ref, ref_id="rid")
    assert "rid" not in pos.columns


def test_all_features_joins_every_overlap():
    """Nested genes are all reported, joined in cluster order; misses get the default."""
    ref = _ref([("chr1", 0, 100, "G1"), ("chr1", 40, 60, "G2")], id_col="gene_id")
    snps = pd.DataFrame({"#CHR": ["chr1", "chr1", "chr1"], "POS0": [50, 10, 500]})
    out, na_idx = iv.assign_pos_to_range_ovlp(
        snps, ref, ref_id="gene_id", fillna="intergenic"
    )
    got = out["gene_id"]
    assert set(got.iloc[0].split(";")) == {"G1", "G2"}
    assert got.iloc[1] == "G1"
    assert got.iloc[2] == "intergenic"
    assert na_idx.tolist() == [2]


def test_all_features_empty_reference():
    """No annotation at all means every SNP is intergenic."""
    snps = pd.DataFrame({"#CHR": ["chr1"], "POS0": [5]})
    out, _ = iv.assign_pos_to_range_ovlp(
        snps, _ref([], id_col="gene_id"), ref_id="gene_id", fillna="intergenic"
    )
    assert out["gene_id"].tolist() == ["intergenic"]


def test_overlaps_any_range_matches_assignment():
    """The boolean variant agrees with assign_pos_to_range on nested ranges."""
    ref = _ref([("chr1", 0, 100, "a"), ("chr1", 40, 60, "b"), ("chr2", 0, 10, "c")])
    qry = pd.DataFrame(
        {"#CHR": ["chr1"] * 3 + ["chr2", "chr9"], "POS0": [50, 150, 0, 5, 5]}
    )
    got = iv.overlaps_any_range(qry, ref)
    assert got.tolist() == [True, False, True, True, False]
    assigned = iv.assign_pos_to_range(qry, ref, ref_id="rid")[0]["rid"].notna()
    assert got.tolist() == assigned.tolist()


def test_overlaps_any_range_empty_reference():
    """No reference ranges means nothing overlaps."""
    qry = pd.DataFrame({"#CHR": ["chr1"], "POS0": [5]})
    assert not iv.overlaps_any_range(qry, _ref([])).any()


@pytest.mark.parametrize("n_ref", [1, 2, 5])
def test_range_clusters_are_non_overlapping(n_ref):
    """Every cluster is internally disjoint, which is what the fast path assumes."""
    starts = np.arange(n_ref) * 5
    ends = starts + 12
    for cluster in iv._range_clusters(starts, ends):
        s, e = starts[cluster], ends[cluster]
        order = np.argsort(s)
        s, e = s[order], e[order]
        assert np.all(s[1:] >= e[:-1])


def test_pos_out_col_names_the_added_column():
    """out_col decouples the added column from the reference id column."""
    qry = pd.DataFrame({"#CHR": ["chr1", "chr1"], "POS0": [5, 15]})
    got, _ = iv.assign_pos_to_range(
        qry, DISJOINT, ref_id="rid", pos_col="POS0", out_col="at_start"
    )
    assert got["at_start"].tolist() == ["a", "b"]
    assert "rid" not in got.columns


@pytest.mark.parametrize(
    "start,end,expect",
    [
        (0, 10, "a"),  # contained, exactly one range
        (2, 8, "a"),  # contained, interior
        (10, 20, "b"),  # contained, starts on a shared edge
        (5, 15, None),  # straddles two adjacent ranges
        (18, 25, None),  # crosses out of the last range
        (30, 40, None),  # in no range at all
    ],
)
def test_contained_needs_both_ends_in_one_id(start, end, expect):
    """rule="contained" yields an id only when both ends land in the same one.

    build_segment_bed relies on this to assert each segment sits inside one arm.
    """
    qry = pd.DataFrame({"#CHR": ["chr1"], "START": [start], "END": [end]})
    got, na_idx = iv.assign_range_to_range(qry, DISJOINT, "rid", rule="contained")
    val = None if pd.isna(got["rid"].iloc[0]) else got["rid"].iloc[0]
    assert val == expect
    assert na_idx.tolist() == ([] if expect else [0])


@pytest.mark.parametrize(
    "n,ranges,expect_cuts",
    [
        (5, [], [1, 2, 3, 4]),  # no range: every item is its own cluster
        (5, [(1, 4)], [1, 4]),  # items 1,2,3 glued; cuts only at its edges
        (5, [(0, 5)], []),  # one range over everything: no cut anywhere
        (5, [(0, 2), (1, 4)], [4]),  # overlapping ranges merge transitively
        (5, [(0, 2), (3, 5)], [2, 3]),  # disjoint ranges stay separate clusters
        (5, [(2, 3)], [1, 2, 3, 4]),  # single-item range blocks nothing
        (5, [(2, 2)], [1, 2, 3, 4]),  # empty range blocks nothing
    ],
)
def test_index_clusters_never_cut_inside_a_range(n, ranges, expect_cuts):
    """merge_ranges_to_clusters: a cut lands only where no range spans the boundary."""
    labels = iv.merge_ranges_to_clusters(n, ranges)
    cuts = [i for i in range(1, n) if labels[i] != labels[i - 1]]
    assert cuts == expect_cuts
    assert labels[0] == 0
    assert (np.diff(labels) >= 0).all(), "cluster ids are not run-length contiguous"


def test_index_clusters_on_an_empty_ordering():
    """Zero items yields an empty int64 array, not an error."""
    labels = iv.merge_ranges_to_clusters(0, [(0, 3)])
    assert labels.shape == (0,)
    assert labels.dtype == np.int64


# --- the uniform (qry, na_idx) contract ------------------------------------------


def test_na_idx_is_positional_into_the_input():
    """na_idx indexes the INPUT frame, so a parallel matrix subsets the same way."""
    qry = pd.DataFrame({"#CHR": ["chr1"] * 4, "POS0": [5, 500, 15, 900]})
    got, na_idx = iv.assign_pos_to_range(qry, DISJOINT, ref_id="rid")
    assert na_idx.tolist() == [1, 3]
    assert got["rid"].isna().to_numpy().nonzero()[0].tolist() == na_idx.tolist()


def test_dropna_filters_reindexes_and_casts():
    """dropna returns the survivors reindexed from 0, cast to the reference dtype."""
    ref = _ref([("chr1", 0, 10, 7), ("chr1", 10, 20, 8)])
    qry = pd.DataFrame({"#CHR": ["chr1"] * 3, "POS0": [5, 900, 15]}, index=[9, 8, 7])
    got, na_idx = iv.assign_pos_to_range(qry, ref, ref_id="rid", dropna=True)
    assert len(got) + len(na_idx) == 3
    assert got.index.tolist() == [0, 1]
    assert got["rid"].tolist() == [7, 8]
    assert got["rid"].dtype == ref["rid"].dtype


def test_fillna_fills_and_does_not_cast():
    """fillna replaces the misses but leaves the column object-typed."""
    ref = _ref([("chr1", 0, 10, 7)])
    qry = pd.DataFrame({"#CHR": ["chr1"] * 2, "POS0": [5, 900]})
    got, na_idx = iv.assign_pos_to_range(qry, ref, ref_id="rid", fillna=-1)
    assert got["rid"].tolist() == [7, -1]
    assert got["rid"].dtype == object
    assert na_idx.tolist() == [1]


def test_fillna_and_dropna_are_exclusive():
    """Filling then dropping is a contradiction, not a silent no-op."""
    qry = pd.DataFrame({"#CHR": ["chr1"], "POS0": [5]})
    with pytest.raises(AssertionError):
        iv.assign_pos_to_range(qry, DISJOINT, ref_id="rid", fillna="x", dropna=True)


def test_duplicate_index_qry_is_annotated_positionally():
    """A duplicate index must not scatter the assignment across matching labels."""
    qry = pd.DataFrame({"#CHR": ["chr1"] * 3, "POS0": [5, 15, 900]}, index=[0, 0, 0])
    got, _ = iv.assign_pos_to_range(qry, DISJOINT, ref_id="rid")
    assert got["rid"].tolist()[:2] == ["a", "b"]
    assert pd.isna(got["rid"].iloc[2])


# --- the 0-based invariant --------------------------------------------------------


def test_pos_col_named_POS_is_rejected():
    """POS is the 1-based VCF column; only POS0 may be assigned."""
    qry = pd.DataFrame({"#CHR": ["chr1"], "POS": [5]})
    with pytest.raises(AssertionError):
        iv.assign_pos_to_range(qry, DISJOINT, ref_id="rid", pos_col="POS")
    with pytest.raises(AssertionError):
        iv.overlaps_any_range(qry, DISJOINT, pos_col="POS")


@pytest.mark.parametrize("start,end", [(10, 10), (10, 5)])
def test_degenerate_reference_range_is_rejected(start, end):
    """A closed-interval or inverted BED must fail loudly, not assign silently."""
    ref = _ref([("chr1", start, end, "a")])
    qry = pd.DataFrame({"#CHR": ["chr1"], "POS0": [10]})
    with pytest.raises(AssertionError):
        iv.assign_pos_to_range(qry, ref, ref_id="rid")


def test_degenerate_query_range_is_rejected():
    """rule="contained" probes END-1, so a zero-length query is meaningless."""
    qry = pd.DataFrame({"#CHR": ["chr1"], "START": [10], "END": [10]})
    with pytest.raises(AssertionError):
        iv.assign_range_to_range(qry, DISJOINT, "rid", rule="contained")


def test_last_base_of_a_range_is_inside_it():
    """END-1 is in the range and END is not - the half-open boundary itself."""
    qry = pd.DataFrame({"#CHR": ["chr1"] * 2, "POS0": [9, 10]})
    got, _ = iv.assign_pos_to_range(qry, _ref([("chr1", 0, 10, "a")]), ref_id="rid")
    assert got["rid"].iloc[0] == "a"
    assert pd.isna(got["rid"].iloc[1])


# --- assign_range_to_range rules --------------------------------------------------


def test_contained_spans_two_reference_rows_sharing_one_id():
    """A segment split by a blacklist hole is one id, so a window over it is contained."""
    ref = _ref([("chr1", 0, 100, "seg1"), ("chr1", 120, 300, "seg1")])
    qry = pd.DataFrame({"#CHR": ["chr1"], "START": [50], "END": [200]})
    got, na_idx = iv.assign_range_to_range(qry, ref, "rid", rule="contained")
    assert got["rid"].iloc[0] == "seg1"
    assert na_idx.tolist() == []


def test_midpoint_matches_a_hand_built_midpoint_assignment():
    """rule="midpoint" is exactly assign_pos_to_range on (START+END)//2."""
    ref = _ref([("chr1", 0, 10, "a"), ("chr1", 10, 30, "b")])
    qry = pd.DataFrame({"#CHR": ["chr1"] * 3, "START": [0, 8, 40], "END": [4, 20, 50]})
    got, na_idx = iv.assign_range_to_range(qry, ref, "rid", rule="midpoint")
    mids = pd.DataFrame({"#CHR": qry["#CHR"], "POS0": (qry["START"] + qry["END"]) // 2})
    want, want_na = iv.assign_pos_to_range(mids, ref, ref_id="rid")
    assert got["rid"].tolist() == want["rid"].tolist()
    assert na_idx.tolist() == want_na.tolist()
    assert got["rid"].tolist()[:2] == ["a", "b"]


def test_unknown_rule_is_rejected():
    """A typo in rule must not silently fall through to max_overlap."""
    qry = pd.DataFrame({"#CHR": ["chr1"], "START": [0], "END": [5]})
    with pytest.raises(AssertionError):
        iv.assign_range_to_range(qry, DISJOINT, "rid", rule="nearest")


# --- trim_range_by_range: interval difference ------------------------------------


def _q(rows, extra=None):
    """Query ranges from (chrom, start, end) tuples, with one carried id column."""
    df = pd.DataFrame(rows, columns=["#CHR", "START", "END"])
    df["qid"] = extra if extra is not None else [f"q{i}" for i in range(len(df))]
    return df


@pytest.mark.parametrize(
    "qry,ref,want",
    [
        # no overlap: passes through
        ([("chr1", 0, 100)], [("chr1", 200, 300)], [(0, 100)]),
        # a bite out of the middle splits it in two
        ([("chr1", 0, 100)], [("chr1", 40, 60)], [(0, 40), (60, 100)]),
        # overlapping the left edge trims the start
        ([("chr1", 0, 100)], [("chr1", 0, 40)], [(40, 100)]),
        # overlapping the right edge trims the end
        ([("chr1", 0, 100)], [("chr1", 60, 100)], [(0, 60)]),
        # fully covered: the range disappears
        ([("chr1", 0, 100)], [("chr1", 0, 100)], []),
        # ref wider than qry on both sides: also disappears
        ([("chr1", 10, 90)], [("chr1", 0, 100)], []),
        # abutting, not overlapping: half-open means no cut
        ([("chr1", 0, 100)], [("chr1", 100, 200)], [(0, 100)]),
        # two bites -> three pieces
        (
            [("chr1", 0, 100)],
            [("chr1", 20, 30), ("chr1", 60, 70)],
            [(0, 20), (30, 60), (70, 100)],
        ),
        # OVERLAPPING refs are merged first, so one gap not two
        (
            [("chr1", 0, 100)],
            [("chr1", 20, 50), ("chr1", 40, 70)],
            [(0, 20), (70, 100)],
        ),
        # NESTED refs likewise
        (
            [("chr1", 0, 100)],
            [("chr1", 20, 70), ("chr1", 30, 40)],
            [(0, 20), (70, 100)],
        ),
    ],
)
def test_trim_range_cases(qry, ref, want):
    """Interval difference, including the overlapping/nested reference cases."""
    out = iv.trim_range_by_range(_q(qry), _q(ref)[["#CHR", "START", "END"]])
    assert list(zip(out["START"], out["END"])) == want


def test_trim_range_carries_every_column_and_order():
    """Any non-coordinate column rides along onto each surviving piece."""
    q = pd.DataFrame(
        {
            "#CHR": ["chr1"],
            "START": [0],
            "END": [100],
            "gene": ["A"],
            "score": [1.5],
        }
    )
    ref = pd.DataFrame({"#CHR": ["chr1"], "START": [40], "END": [60]})
    out = iv.trim_range_by_range(q, ref)
    assert list(out.columns) == list(q.columns)
    assert out["gene"].tolist() == ["A", "A"]
    assert out["score"].tolist() == [1.5, 1.5]


def test_trim_range_untouched_chromosome_passes_through():
    """A query on a contig the reference never mentions is returned as-is."""
    q = _q([("chr1", 0, 100), ("chr9", 0, 100)])
    ref = pd.DataFrame({"#CHR": ["chr1"], "START": [0], "END": [50]})
    out = iv.trim_range_by_range(q, ref)
    assert list(zip(out["#CHR"], out["START"], out["END"])) == [
        ("chr1", 50, 100),
        ("chr9", 0, 100),
    ]


def test_trim_range_empty_reference_is_identity():
    """No reference ranges means nothing is cut."""
    q = _q([("chr1", 0, 100), ("chr2", 0, 50)])
    empty = pd.DataFrame({"#CHR": [], "START": [], "END": []})
    out = iv.trim_range_by_range(q, empty)
    assert len(out) == 2 and out["END"].tolist() == [100, 50]


def test_trim_range_rejects_degenerate_ranges():
    """The module's 0-based half-open invariant applies here too."""
    bad = _q([("chr1", 50, 50)])
    ref = pd.DataFrame({"#CHR": ["chr1"], "START": [0], "END": [10]})
    with pytest.raises(AssertionError):
        iv.trim_range_by_range(bad, ref)


# --- split_range_at_pos: cut at SV extremities -----------------------------------


def _p(rows):
    """Cut positions from (chrom, pos0) tuples."""
    return pd.DataFrame(rows, columns=["#CHR", "POS0"])


@pytest.mark.parametrize(
    "qry,pos,want",
    [
        # one interior cut -> two pieces, no base lost
        ([("chr1", 0, 100)], [("chr1", 40)], [(0, 40), (40, 100)]),
        # two cuts -> three pieces
        (
            [("chr1", 0, 100)],
            [("chr1", 20), ("chr1", 70)],
            [(0, 20), (20, 70), (70, 100)],
        ),
        # unsorted input positions cut the same way
        (
            [("chr1", 0, 100)],
            [("chr1", 70), ("chr1", 20)],
            [(0, 20), (20, 70), (70, 100)],
        ),
        # a repeated position cuts once
        ([("chr1", 0, 100)], [("chr1", 40), ("chr1", 40)], [(0, 40), (40, 100)]),
        # on START: no-op, the position is already a bound
        ([("chr1", 0, 100)], [("chr1", 0)], [(0, 100)]),
        # on END: no-op, END is outside the half-open range
        ([("chr1", 0, 100)], [("chr1", 100)], [(0, 100)]),
        # the last base is inside, so it cuts
        ([("chr1", 0, 100)], [("chr1", 99)], [(0, 99), (99, 100)]),
        # outside every range: no-op
        ([("chr1", 0, 100)], [("chr1", 150)], [(0, 100)]),
        # a contig the ranges never mention: no-op
        ([("chr1", 0, 100)], [("chr9", 40)], [(0, 100)]),
        # each range is cut independently
        (
            [("chr1", 0, 100), ("chr2", 0, 100)],
            [("chr1", 40)],
            [(0, 40), (40, 100), (0, 100)],
        ),
    ],
)
def test_split_range_cases(qry, pos, want):
    """A cut lands on a bound, so pieces tile the input range exactly."""
    out = iv.split_range_at_pos(_q(qry), _p(pos))
    assert list(zip(out["START"], out["END"])) == want


def test_split_range_carries_every_column_and_order():
    """Any non-coordinate column rides along onto each piece, input order kept."""
    q = pd.DataFrame(
        {
            "#CHR": ["chr1", "chr1"],
            "START": [0, 100],
            "END": [100, 200],
            "region_id": ["1p", "1q"],
        }
    )
    out = iv.split_range_at_pos(q, _p([("chr1", 40)]))
    assert list(out.columns) == list(q.columns)
    assert out["region_id"].tolist() == ["1p", "1p", "1q"]
    assert list(zip(out["START"], out["END"])) == [(0, 40), (40, 100), (100, 200)]


def test_split_range_empty_positions_is_identity():
    """No cut positions means the ranges pass through."""
    q = _q([("chr1", 0, 100), ("chr2", 0, 50)])
    empty = pd.DataFrame({"#CHR": [], "POS0": []})
    out = iv.split_range_at_pos(q, empty)
    assert len(out) == 2 and out["END"].tolist() == [100, 50]


def test_split_range_rejects_one_based_pos_col():
    """`POS` is 1-based; the module takes `POS0` only."""
    q = _q([("chr1", 0, 100)])
    pos = pd.DataFrame({"#CHR": ["chr1"], "POS": [40]})
    with pytest.raises(AssertionError):
        iv.split_range_at_pos(q, pos, pos_col="POS")


def test_split_range_rejects_degenerate_ranges():
    """The module's 0-based half-open invariant applies here too."""
    bad = _q([("chr1", 50, 50)])
    with pytest.raises(AssertionError):
        iv.split_range_at_pos(bad, _p([("chr1", 40)]))


def test_overlap_bp_partial_and_full():
    """Overlapping bp per query range, counting the half-open intersection."""
    qry = _q([("chr1", 0, 100), ("chr1", 100, 200), ("chr1", 500, 600)])
    ref = _ref([("chr1", 50, 150, "t1")])
    assert iv.range_overlap_bp(qry, ref).tolist() == [50, 50, 0]


def test_overlap_bp_unions_overlapping_references():
    """A base covered by two probes counts once, not twice."""
    qry = _q([("chr1", 0, 100)])
    ref = _ref([("chr1", 10, 60, "a"), ("chr1", 40, 80, "b")])
    assert iv.range_overlap_bp(qry, ref).tolist() == [70]


def test_overlap_bp_reference_inside_query():
    """Several disjoint references inside one query sum."""
    qry = _q([("chr1", 0, 1000)])
    ref = _ref(
        [("chr1", 10, 20, "a"), ("chr1", 100, 130, "b"), ("chr1", 900, 905, "c")]
    )
    assert iv.range_overlap_bp(qry, ref).tolist() == [45]


def test_overlap_bp_other_chromosome_and_empty():
    """A contig the reference never mentions, and an empty reference, are both zero."""
    qry = _q([("chr1", 0, 100), ("chr9", 0, 100)])
    ref = _ref([("chr1", 0, 100, "a")])
    assert iv.range_overlap_bp(qry, ref).tolist() == [100, 0]
    assert iv.range_overlap_bp(qry, _ref([])).tolist() == [0, 0]


def test_overlap_bp_preserves_query_row_order():
    """Output is positional in the query, whatever order the reference arrives in."""
    qry = _q([("chr2", 0, 100), ("chr1", 0, 100)])
    ref = _ref([("chr1", 0, 30, "a"), ("chr2", 0, 70, "b")])
    assert iv.range_overlap_bp(qry, ref).tolist() == [70, 30]
