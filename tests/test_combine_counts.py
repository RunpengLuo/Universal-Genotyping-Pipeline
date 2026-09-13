"""Unit tests for the binning core and the bb-level allele aggregation.

Last update: 2026-08-11

Covers:
- depth: fixed-bin depth aggregated per bb, NaN masked per dataset column
- summation: bb BAF is correct only when SNPs are haplotype-consistent
- phase flips: the split preserves an LOH that naive summation cancels
- binning: cluster keys, SNP-free bins and segments, join-able bbs
- fragments: ATAC counts route through the windows, holes dropped
- tiling: build_fixedwidth_bins matches the loop it replaced
"""

import os
import sys

import numpy as np
import pandas as pd
import pytest

_REPO = os.path.join(os.path.dirname(__file__), "..")
_SCRIPTS = os.path.join(_REPO, "workflow", "scripts")
sys.path.insert(0, os.path.join(_REPO, "config"))
for _sub in ("", "script_utils", "plot"):
    sys.path.insert(0, os.path.join(_SCRIPTS, _sub))
ph = pytest.importorskip("phasing_utils")
ru = pytest.importorskip("range_utils")


@pytest.fixture
def seg():
    """segmentation_utils: the summation kernel and the binning kernel (numba)."""
    pytest.importorskip("numba")
    return pytest.importorskip("segmentation_utils")


@pytest.fixture
def feat():
    """feature_utils; pandas/numpy/scipy only."""
    return pytest.importorskip("feature_utils")


@pytest.fixture
def ccu():
    """combine_counts_utils: the bulk depth and RDR summaries."""
    return pytest.importorskip("combine_counts_utils")


def _bin_snps(snps, bins):
    """Assign SNPs to fixed bins, dropping the ones that land in none."""
    binned, _ = ru.assign_pos_to_range(snps, bins, ref_id="bin_id", dropna=True)
    return binned


def _snps(n, region="r", ps=1):
    """SNP frame with the columns detect_phase_flips expects (incl. #CHR/POS0), index 0..n-1."""
    return pd.DataFrame(
        {
            "#CHR": ["chr1"] * n,
            "POS0": np.arange(n) * 1000,
            "region_id": [region] * n,
            "PS": [ps] * n,
        },
        index=range(n),
    )


def _dense(x):
    return np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x)


def test_sum_features_to_bbs_sums_per_bb(seg):
    x = np.array([[3.0, 0.0], [1.0, 0.0], [2.0, 5.0]])
    out = _dense(seg.sum_features_to_bbs(x, np.array([0, 0, 1]), 2))
    assert out.tolist() == [[4.0, 0.0], [2.0, 5.0]]


def test_same_orientation_bb_baf_no_cancellation(seg):
    # 3 LOH SNPs, all B-skewed (BAF 0.8) -> bin BAF stays 0.8
    b = np.array([[8.0], [7.0], [9.0]])
    tot = b + np.array([[2.0], [3.0], [1.0]])
    bins = np.array([0, 0, 0])
    baf = (
        _dense(seg.sum_features_to_bbs(b, bins, 1))[0, 0]
        / _dense(seg.sum_features_to_bbs(tot, bins, 1))[0, 0]
    )
    assert abs(baf - 0.8) < 1e-9


def test_detect_flip_at_orientation_switch():
    # high depth -> tight CIs; SNPs 0-2 BAF~0.8, SNPs 3-5 BAF~0.2 (orientation flip)
    b = np.array([[80], [80], [80], [20], [20], [20]], float)
    a = np.array([[20], [20], [20], [80], [80], [80]], float)
    pg = ph.detect_phase_flips(
        _snps(6), a, b, cluster_cols=["region_id", "PS"]
    ).to_numpy()
    assert int((np.diff(pg) != 0).sum()) == 1  # exactly one boundary
    assert pg[2] != pg[3]  # ...between SNP 2 and 3
    assert pg[0] == pg[1] == pg[2] and pg[3] == pg[4] == pg[5]


def test_no_flip_when_balanced():
    # BAF ~0.5 everywhere -> no orientation to switch -> one phase cluster
    b = np.full((6, 1), 50.0)
    a = np.full((6, 1), 50.0)
    pg = ph.detect_phase_flips(
        _snps(6), a, b, cluster_cols=["region_id", "PS"]
    ).to_numpy()
    assert len(np.unique(pg)) == 1


def test_loh_cancelled_when_merged_but_preserved_when_split(seg):
    # LOH with a mid-segment orientation flip: SNPs 0-2 B=0.8, SNPs 3-5 B=0.2
    b = np.array([[80], [80], [80], [20], [20], [20]], float)
    a = np.array([[20], [20], [20], [80], [80], [80]], float)
    tot = a + b

    # naive: one bin ignoring phase -> B and A cancel -> BAF 0.5 (LOH lost)
    merged = (
        _dense(seg.sum_features_to_bbs(b, np.zeros(6, int), 1))[0, 0]
        / _dense(seg.sum_features_to_bbs(tot, np.zeros(6, int), 1))[0, 0]
    )
    assert abs(merged - 0.5) < 1e-9

    # phase-aware: split at the detected flip -> each bin keeps the LOH (folded 0.2)
    pg = ph.detect_phase_flips(
        _snps(6), a, b, cluster_cols=["region_id", "PS"]
    ).to_numpy()
    bins, _ = pd.factorize(pg)
    k = int(bins.max()) + 1
    bsum = _dense(seg.sum_features_to_bbs(b, bins, k))[:, 0]
    tsum = _dense(seg.sum_features_to_bbs(tot, bins, k))[:, 0]
    folded = np.minimum(bsum / tsum, 1 - bsum / tsum)
    assert np.allclose(folded, 0.2, atol=1e-9)


def test_bbs_frame_joins_on_bb_id(seg):
    """`bbs` must be join-able on bb_id: the index name must not shadow the column.

    build_adaptive_bins groups the fixed bins by bb_id and then adds a bb_id column.
    If the grouped index keeps that name, pandas rejects `join(on="bb_id")` with
    "both an index level and a column label", which is how interp_cM_between_bbs
    consumes the frame.
    """
    n_bins = 40
    bins = pd.DataFrame(
        {
            "#CHR": ["chr1"] * n_bins,
            "START": np.arange(n_bins) * 1000,
            "END": (np.arange(n_bins) + 1) * 1000,
            "region_id": ["r"] * n_bins,
        }
    )
    bins["bin_id"] = np.arange(n_bins)
    snps = pd.DataFrame(
        {
            "#CHR": ["chr1"] * 80,
            "POS0": np.arange(80) * 500,
            "POS": np.arange(80) * 500 + 1,
            "region_id": ["r"] * 80,
        }
    )
    tot = np.full((80, 2), 10.0)
    bbs, snps_bb = seg.build_adaptive_bins(
        bins,
        _bin_snps(snps, bins),
        tot,
        50,
        2,
        cluster_cols=["region_id"],
    )
    assert bbs.index.name is None, "bb_id index name shadows the bb_id column"
    assert "bb_id" in bbs.columns and "bb_id" in snps_bb.columns
    # the exact operation interp_cM_between_bbs performs
    span = snps_bb.groupby("bb_id", sort=False)["POS"].agg(lo="min", hi="max")
    joined = bbs.join(span, on="bb_id")
    assert len(joined) == len(bbs) and joined["lo"].notna().all()


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"PASS {name}")
    print("all combine_counts unit tests passed")


# ---------------------------------------------------------------------------
# Single-cell binning on the window grid: the fixed bins are windows, so bins
# exist where no SNP does and none lies in a blacklist hole.
# ---------------------------------------------------------------------------


def _windows(n, chrom="chr1", size=1000, start=0):
    """n contiguous fixed-size windows carrying one region_id/seg_id."""
    starts = start + np.arange(n) * size
    win = pd.DataFrame(
        {
            "#CHR": [chrom] * n,
            "START": starts,
            "END": starts + size,
            "region_id": ["arm"] * n,
            "seg_id": ["seg"] * n,
        }
    )
    win["bin_id"] = np.arange(n)
    return win


def _snps_in(windows, rows):
    """One SNP at the midpoint of each named window row."""
    mids = ((windows["START"] + windows["END"]) // 2).to_numpy()[rows]
    return pd.DataFrame(
        {
            "#CHR": ["chr1"] * len(rows),
            "POS0": mids,
            "POS": mids + 1,
            "region_id": ["arm"] * len(rows),
            "seg_id": ["seg"] * len(rows),
        }
    )


def test_leading_snp_free_bin_keeps_a_cluster_key(seg):
    """A bin with no SNP must still carry every cluster key.

    build_adaptive_bins groups the bins by cluster_cols and pandas drops null keys, so a
    null-keyed bin would never enter the merge and would keep the initialized bb_id of 0.
    With a 1 kb grid the first bin almost never holds a SNP, so this is the common case.
    """
    win = _windows(10)
    snps = _snps_in(win, [5, 6, 7])
    snps["PS"] = 1
    tot = np.full((len(snps), 2), 100.0)
    binned = _bin_snps(snps, win)
    modal = binned.groupby("bin_id")["PS"].agg(lambda x: x.mode().iloc[0])
    win["PS"] = win["bin_id"].map(modal).ffill().bfill().fillna(1)

    assert win["PS"].notna().all(), "leading SNP-free bins lost their cluster key"
    bbs, _ = seg.build_adaptive_bins(
        win,
        binned,
        tot,
        50,
        1,
        cluster_cols=["region_id", "seg_id", "PS"],
        gene_aware=False,
    )
    assert win["bb_id"].between(0, len(bbs) - 1).all()


def test_snp_free_segment_gets_its_own_bb(seg):
    """Windows exist without SNPs, so a SNP-free segment still yields a bb."""
    win = pd.concat(
        [_windows(5), _windows(5, start=100_000).assign(seg_id="seg2")],
        ignore_index=True,
    )
    win["bin_id"] = np.arange(len(win))
    snps = _snps_in(win, [1, 2, 3])  # every SNP lands in seg, none in seg2
    snps["PS"] = 1
    tot = np.full((len(snps), 2), 100.0)
    binned = _bin_snps(snps, win)
    modal = binned.groupby("bin_id")["PS"].agg(lambda x: x.mode().iloc[0])
    win["PS"] = win["bin_id"].map(modal).ffill().bfill().fillna(1)
    bbs, _ = seg.build_adaptive_bins(
        win,
        binned,
        tot,
        50,
        1,
        cluster_cols=["region_id", "seg_id", "PS"],
        gene_aware=False,
    )
    assert (bbs["#SNPS"] == 0).any(), "the SNP-free segment produced no bb"
    assert sorted(bbs["bb_id"]) == list(range(len(bbs))), "bb_ids are not contiguous"


def test_switchprobs_finite_across_a_snp_free_bb(seg):
    """A bb holding no SNP must not interpolate to NaN, nor poison the next bb."""
    bbs = pd.DataFrame(
        {
            "#CHR": ["chr1"] * 3,
            "START": [0, 10_000, 20_000],
            "END": [10_000, 20_000, 30_000],
            "bb_id": [0, 1, 2],
        }
    )
    # bb 1 holds no SNP
    snp_info = pd.DataFrame(
        {"#CHR": ["chr1"] * 2, "POS": [5_000, 25_000], "bb_id": [0, 2]}
    )
    gmap = pd.DataFrame(
        {
            "#CHR": ["chr1"] * 4,
            "POS": [0, 10_000, 20_000, 30_000],
            "cM": [0.0, 1.0, 2.0, 3.0],
        }
    )
    dist = ph.interp_cM_between_bbs(bbs, snp_info, gmap, bb_id_col="bb_id")
    assert np.isfinite(dist).all(), "a SNP-free bb produced NaN cM distances"
    probs = ph.estimate_switchprobs_cM(dist)
    assert np.isfinite(probs).all()


def test_gene_cluster_over_bins_never_splits_a_gene(seg, feat):
    """A gene's whole span of bins is one cluster, including its SNP-free interior."""
    win = _windows(10)
    snps = _snps_in(win, [3, 7])  # one gene, SNPs only at its two ends
    snps["feature_id"] = ["G", "G"]
    binned = _bin_snps(snps, win)
    spans = (
        feat.explode_feature_ids(binned, cols=["bin_id"])
        .groupby("feature_id")["bin_id"]
        .agg(["min", "max"])
    )
    win["gene_cluster"] = ru.merge_ranges_to_clusters(
        len(win), zip(spans["min"].to_numpy(), spans["max"].to_numpy() + 1)
    )
    inside = win.loc[3:7, "gene_cluster"].to_numpy()
    assert len(set(inside)) == 1, "bins 3..7 of one gene fell into several clusters"
    assert win.loc[2, "gene_cluster"] != inside[0]


def test_window_frame_routes_fragments_to_the_owning_bb(feat, tmp_path):
    """ATAC counts through the windows: holes are dropped, sibling windows accumulate."""
    frag = tmp_path / "atac_fragments.tsv.gz"
    import gzip

    rows = [
        ("chr1", 100, 200, "AAA"),  # window 0 -> bb 0
        ("chr1", 1100, 1200, "AAA"),  # window 1 -> bb 0
        ("chr1", 2100, 2200, "AAA"),  # the HOLE -> dropped
        ("chr1", 3100, 3200, "BBB"),  # window 2 -> bb 1
    ]
    with gzip.open(frag, "wt") as fh:
        for c, s, e, bc in rows:
            fh.write(f"{c}\t{s}\t{e}\t{bc}\t1\n")

    # windows 0,1 own bb 0; window 2 owns bb 1; 2000-3000 is a blacklist hole
    win = pd.DataFrame(
        {
            "#CHR": ["chr1"] * 3,
            "START": [0, 1000, 3000],
            "END": [1000, 2000, 4000],
            "bb_id": [0, 0, 1],
        }
    )
    barcodes = pd.DataFrame(
        {
            "raw": ["AAA", "BBB"],
            "dataset_id": ["R1", "R1"],
            "assay_type": ["scATAC", "scATAC"],
            "BARCODE": ["AAA_R1_scATAC", "BBB_R1_scATAC"],
        }
    )
    mtx = feat.sum_atac_fragments_to_bins([str(frag)], ["R1"], barcodes, win, 2)
    dense = mtx.toarray()
    assert dense.shape == (2, 2)
    assert dense[0, 0] == 2, "two windows of one bb did not accumulate"
    assert dense[1, 1] == 1
    assert dense.sum() == 3, "the fragment in the blacklist hole was counted"


# ---------------------------------------------------------------------------
# build_fixedwidth_bins: the vectorized segment tiler
# ---------------------------------------------------------------------------


def _tile_loop(start, end, size):
    """The pre-vectorization reference implementation."""
    rows, pos = [], start
    while pos < end:
        e = min(pos + size, end)
        rows.append([pos, e])
        pos = e
    if len(rows) > 1 and (rows[-1][1] - rows[-1][0]) < size // 2:
        rows[-2][1] = rows[-1][1]
        rows.pop()
    return rows


@pytest.mark.parametrize(
    "start,end",
    [
        (0, 3000),  # exact multiple
        (0, 2600),  # remainder >= size//2 -> its own bin
        (0, 2400),  # remainder <  size//2 -> absorbed by its predecessor
        (0, 2500),  # remainder == size//2, the boundary
        (0, 400),  # shorter than one bin
        (0, 1000),  # exactly one bin
        (0, 1001),  # one bin + a 1 bp tail
        (0, 1),  # single base
        (7_000_000, 7_002_450),  # nonzero offset
    ],
)
def test_fixedwidth_bins_match_the_tiling_loop(seg, start, end):
    """Vectorized tiling reproduces the loop it replaced, remainder rule included."""
    segments = pd.DataFrame(
        {"#CHR": ["chr1"], "START": [start], "END": [end], "seg_id": ["s"]}
    )
    got = seg.build_fixedwidth_bins(segments, 1000)[["START", "END"]].values.tolist()
    assert got == _tile_loop(start, end, 1000)


def test_fixedwidth_bins_carry_segment_ids(seg):
    """Each bin inherits its source segment's ids, so no assignment pass is needed."""
    segments = pd.DataFrame(
        {
            "#CHR": ["chr1", "chr1", "chr2"],
            "START": [0, 5000, 0],
            "END": [2400, 7600, 1500],
            "region_id": ["1p", "1q", "2p"],
            "seg_id": ["a", "b", "c"],
        }
    )
    out = seg.build_fixedwidth_bins(segments, 1000, chroms=["chr1"])
    assert out["#CHR"].unique().tolist() == ["chr1"]
    assert out["seg_id"].tolist() == ["a", "a", "b", "b", "b"]
    assert out["region_id"].tolist() == ["1p", "1p", "1q", "1q", "1q"]
    # every bin lies inside the segment it came from
    span = segments.set_index("seg_id")
    for _, r in out.iterrows():
        assert span.loc[r["seg_id"], "START"] <= r["START"]
        assert r["END"] <= span.loc[r["seg_id"], "END"]


def test_fixedwidth_bins_on_empty_and_zero_length(seg):
    """A zero-length segment yields no bin; an empty frame yields an empty frame."""
    zero = pd.DataFrame({"#CHR": ["chr1"], "START": [5], "END": [5], "seg_id": ["s"]})
    assert len(seg.build_fixedwidth_bins(zero, 1000)) == 0
    empty = pd.DataFrame({"#CHR": [], "START": [], "END": [], "seg_id": []})
    assert len(seg.build_fixedwidth_bins(empty, 1000)) == 0


def test_summarize_read_depth_bb_masks_nan_per_column(ccu):
    """Length-weighted mean per bb; a NaN window is dropped for that column only.

    The whole point of one joint depth matrix: rd_correct keeps NaN in place instead of
    deleting the row, so a window the correction lost for one dataset still counts for
    every other dataset of the run.
    """
    bins = pd.DataFrame(
        {
            "#CHR": ["chr1"] * 4,
            "START": [0, 1000, 2000, 3000],
            "END": [1000, 2000, 3000, 5000],
            "bb_id": [0, 0, 1, 1],
        }
    )
    dp = np.array(
        [
            [10.0, 10.0],
            [20.0, np.nan],
            [30.0, 30.0],
            [40.0, 40.0],
        ],
        dtype=np.float32,
    )
    bb_dp, bb_bases = ccu.summarize_read_depth_bb(bins, dp, 2)

    assert bb_dp[0, 0] == pytest.approx(15.0)
    assert bb_dp[0, 1] == pytest.approx(10.0)
    assert bb_dp[1, 0] == pytest.approx((30 * 1000 + 40 * 2000) / 3000)
    assert bb_dp[1, 1] == pytest.approx(bb_dp[1, 0])
    assert bb_bases[0, 0] == pytest.approx(30000.0)
    assert bb_bases[0, 1] == pytest.approx(10000.0)


def test_summarize_read_depth_bb_all_nan_bb_is_nan(ccu):
    """A bb whose every window is NaN for a dataset yields NaN, not a divide-by-zero."""
    bins = pd.DataFrame(
        {"#CHR": ["chr1"] * 2, "START": [0, 1000], "END": [1000, 2000], "bb_id": [0, 0]}
    )
    dp = np.array([[np.nan, 5.0], [np.nan, 5.0]], dtype=np.float32)
    bb_dp, bb_bases = ccu.summarize_read_depth_bb(bins, dp, 1)
    assert np.isnan(bb_dp[0, 0])
    assert bb_dp[0, 1] == pytest.approx(5.0)
    assert bb_bases[0, 0] == 0.0


##################################################
# the close rule: min_snp_reads AND min_total_reads, no span cap


def _bin_all(
    seg, win, snps, tot, msr, min_snp_per_bin=1, rdcount=None, mtr=0, is_loh_col=None
):
    """build_adaptive_bins over one window frame, SNPs already binned."""
    return seg.build_adaptive_bins(
        win,
        _bin_snps(snps, win),
        tot,
        msr,
        min_snp_per_bin,
        cluster_cols=["region_id", "seg_id"]
        + ([is_loh_col.replace("is_loh", "loh_id")] if is_loh_col else []),
        gene_aware=False,
        bin_total_reads=rdcount,
        min_total_reads=mtr,
        is_loh_col=is_loh_col,
    )


def test_msr_alone_closes_when_mtr_is_off(seg):
    """With min_total_reads=0 the rule is the old one: SNP reads decide."""
    win = _windows(10)
    snps = _snps_in(win, list(range(10)))
    tot = np.full((len(snps), 1), 50.0)  # 50 reads per window, MSR 100 -> 2 windows/bb
    bbs, _ = _bin_all(seg, win, snps, tot, 100)
    assert len(bbs) == 5, bbs
    assert (bbs["BLOCKSIZE"] == 2000).all()


def test_mtr_forces_a_larger_bb_than_msr_alone(seg):
    """A read-start floor the SNP reads do not reach keeps merging past the MSR cut."""
    win = _windows(12)
    snps = _snps_in(win, list(range(12)))
    tot = np.full((len(snps), 1), 50.0)
    rd = np.full((12, 1), 100.0)  # MTR 400 needs 4 windows, MSR 100 only 2
    bbs, _ = _bin_all(seg, win, snps, tot, 100, rdcount=rd, mtr=400)
    assert len(bbs) == 3, bbs
    assert (bbs["BLOCKSIZE"] == 4000).all()


def test_mtr_covers_every_column_not_just_tumors(seg):
    """min_total_reads is per observation; the thinnest column sets the bb."""
    win = _windows(12)
    snps = _snps_in(win, list(range(12)))
    tot = np.full((len(snps), 1), 500.0)
    rd = np.column_stack([np.full(12, 400.0), np.full(12, 100.0)])
    bbs, _ = _bin_all(seg, win, snps, tot, 100, rdcount=rd, mtr=400)
    assert (bbs["BLOCKSIZE"] == 4000).all(), "the 100-read column must drive the merge"


def test_loh_cluster_ignores_the_snp_criterion(seg):
    """No germline het survives clonal LOH, so read starts alone close the bb."""
    win = _windows(8)
    win["is_loh"] = True
    win["loh_id"] = 1
    snps = _snps_in(win, [0])  # a single SNP, far below any MSR
    tot = np.full((len(snps), 1), 1.0)
    rd = np.full((8, 1), 100.0)
    bbs, _ = _bin_all(
        seg, win, snps, tot, 10_000, rdcount=rd, mtr=200, is_loh_col="is_loh"
    )
    assert len(bbs) == 4, bbs
    assert (bbs["BLOCKSIZE"] == 2000).all()


def test_loh_boundary_never_falls_inside_a_bb(seg):
    """loh_id in cluster_cols is what stops a bb straddling the boundary."""
    win = _windows(12)
    win["is_loh"] = np.r_[np.zeros(6, bool), np.ones(6, bool)]
    win["loh_id"] = np.r_[np.ones(6, int), np.full(6, 2)]
    snps = _snps_in(win, list(range(12)))
    tot = np.full((len(snps), 1), 50.0)
    rd = np.full((12, 1), 100.0)
    bbs, _ = _bin_all(
        seg, win, snps, tot, 100, rdcount=rd, mtr=200, is_loh_col="is_loh"
    )
    loh_of_bb = win.groupby("bb_id")["is_loh"].nunique()
    assert (loh_of_bb == 1).all(), "a bb spans the LOH boundary"
    assert bbs["END"].isin([6000]).any(), "no bb ends at the LOH boundary"


def test_cluster_end_merges_the_short_trailing_run(seg):
    """The remainder joins the last bb rather than forming an under-threshold one."""
    win = _windows(5)
    snps = _snps_in(win, list(range(5)))
    tot = np.full((len(snps), 1), 50.0)  # MSR 100 -> 2 windows per bb, 1 left over
    bbs, _ = _bin_all(seg, win, snps, tot, 100)
    assert len(bbs) == 2
    assert bbs["BLOCKSIZE"].tolist() == [2000, 3000], "the tail did not merge back"


def test_a_cluster_that_meets_nothing_stays_one_bb(seg):
    """HATCHet2's one-bin-per-arm case: no span cap can cut it, so it is one bb."""
    win = _windows(6)
    snps = _snps_in(win, [0])
    tot = np.full((len(snps), 1), 1.0)
    rd = np.full((6, 1), 1.0)
    bbs, _ = _bin_all(seg, win, snps, tot, 10_000, rdcount=rd, mtr=10_000)
    assert len(bbs) == 1
    assert bbs["BLOCKSIZE"].tolist() == [6000]


##################################################
# clonal-LOH detection: density per tile, not allele fraction per site


@pytest.fixture
def loh():
    """loh_utils: the tile builder, the NB fit and the two-state decode."""
    return pytest.importorskip("loh_utils")


def _tiles(n, het_per_tile, tile_kb=100, chrom="chr1"):
    """n contiguous tiles of full exposure, each holding the given het count."""
    starts = np.arange(n) * tile_kb * 1000
    return pd.DataFrame(
        {
            "#CHR": [chrom] * n,
            "START": starts,
            "END": starts + tile_kb * 1000,
            "region_id": ["arm"] * n,
            "n_het": np.asarray(het_per_tile),
            "exposure_bp": tile_kb * 1000,
            "mid": starts + tile_kb * 500,
        }
    )


def _decode(loh, tiles, ratio=0.02, t=1e-6):
    e = tiles["exposure_bp"].to_numpy() / 1e6
    rate, size = loh.fit_snp_rate(tiles["n_het"].to_numpy(), e)
    return loh.viterbi_loh(
        tiles["n_het"].to_numpy(),
        e,
        tiles["mid"].to_numpy(),
        pd.factorize(tiles["region_id"])[0],
        rate,
        ratio * rate,
        size,
        t,
    )


def test_fit_snp_rate_recovers_a_planted_rate(loh):
    """The rate is per Mb, so a tile's count scales with its assayable span."""
    rng = np.random.default_rng(0)
    e = rng.uniform(0.2, 1.0, 400)
    x = rng.poisson(600 * e)
    rate, size = loh.fit_snp_rate(x, e)
    assert 540 < rate < 660, rate
    assert size > 0


def test_fit_snp_rate_ignores_tiles_without_exposure(loh):
    """A window with no coverage contributes no exposure and must not skew the rate."""
    e = np.r_[np.ones(50), np.zeros(50)]
    x = np.r_[np.full(50, 600), np.zeros(50, dtype=int)]
    rate, _ = loh.fit_snp_rate(x, e)
    assert 590 < rate < 610, rate


def test_decode_finds_a_planted_depleted_block(loh):
    """The signal this exists for: het density collapsing over a run of tiles."""
    rng = np.random.default_rng(1)
    het = rng.poisson(600, 120)
    het[40:80] = rng.poisson(3, 40)
    is_loh = _decode(loh, _tiles(120, het))
    assert is_loh[40:80].mean() > 0.9
    assert is_loh[:40].sum() == 0 and is_loh[80:].sum() == 0


def test_hom_runs_inside_a_het_rich_tile_are_not_loh(loh):
    """The per-site chain's failure, as a unit test.

    Hom-looking SNPs cluster hard - on this data 30% of them sit in runs of 30+ - and a
    per-site chain flips on ~48 of them in a row. A tile still holding most of its hets
    is not depleted, whatever the local runs look like.
    """
    rng = np.random.default_rng(2)
    het = rng.poisson(600, 120)
    het[40:80] = rng.poisson(400, 40)  # a third of the hets gone, not all of them
    assert _decode(loh, _tiles(120, het)).sum() == 0


def test_decode_restarts_at_every_arm(loh):
    """A depleted arm must not drag the state across the centromere."""
    tiles = _tiles(80, np.r_[np.full(40, 3), np.full(40, 600)])
    tiles["region_id"] = ["chr1p"] * 40 + ["chr1q"] * 40
    is_loh = _decode(loh, tiles)
    assert is_loh[:40].all() and not is_loh[40:].any()


def test_build_loh_tiles_never_spans_a_segment(loh):
    """seg_id bounds cut a tile, so a tile never straddles a blacklist hole."""
    win = pd.concat(
        [_windows(120), _windows(120, start=120_000).assign(seg_id="seg2")],
        ignore_index=True,
    )
    win["bin_id"] = np.arange(len(win))
    tiles = loh.build_loh_tiles(
        win, np.array([], dtype=int), np.ones(len(win), bool), 100_000
    )
    assert len(tiles) == 4, tiles  # 2 segments x 2 tile boundaries
    assert (tiles["END"] - tiles["START"] <= 100_000).all()


def test_build_loh_tiles_drops_exposure_where_unassayable(loh):
    """An uncovered window keeps its hets but lends no exposure, so a hole is not LOH."""
    win = _windows(200)
    assayable = np.ones(200, bool)
    assayable[:100] = False
    tiles = loh.build_loh_tiles(win, np.array([0, 1, 150]), assayable, 100_000)
    assert tiles["exposure_bp"].tolist() == [0, 100 * 1000]
    assert tiles["n_het"].tolist() == [2, 1], "hets are kept where exposure is not"


def test_loh_intervals_report_a_single_tile_run(loh):
    """No length filter: the tile size is the floor, so one tile is one interval."""
    tiles = _tiles(20, np.full(20, 600))
    is_loh = np.zeros(20, bool)
    is_loh[5] = True
    out = loh.loh_intervals(tiles, is_loh)
    assert len(out) == 1
    assert (out["END"] - out["START"]).tolist() == [100_000]


def test_loh_intervals_do_not_bridge_a_neutral_gap(loh):
    """A gap is evidence against; merging over it is what chained spurious runs."""
    tiles = _tiles(20, np.full(20, 600))
    is_loh = np.zeros(20, bool)
    is_loh[2:6] = True
    is_loh[8:12] = True
    out = loh.loh_intervals(tiles, is_loh)
    assert len(out) == 2, out
    assert (out["END"] - out["START"]).tolist() == [400_000, 400_000]


def test_loh_intervals_empty_without_loh(loh):
    tiles = _tiles(10, np.full(10, 600))
    out = loh.loh_intervals(tiles, np.zeros(10, bool))
    assert len(out) == 0
    assert list(out.columns) == ["#CHR", "START", "END"]
