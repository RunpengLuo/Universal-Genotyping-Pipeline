"""Unit tests for the bulk combine_counts allele/BAF aggregation core.

Targets the two functions that decide bb-level BAF correctness:

- `matrix_utils.sum_features_to_bbs`: sum SNP-by-observation counts into
  bb-by-observation counts.
- `phasing_utils.detect_phase_flips`: split SNPs into phase clusters at
  haplotype-orientation switches so a bb never sums across a flip.

bb BAF is `sum(B) / sum(total)` over a bb's SNPs, which is only correct when the
SNPs are haplotype-consistent. These tests show the phase-flip split preserves an LOH
signal that naive (phase-blind) summation cancels to 0.5.

Both modules need only numpy/pandas/scipy, so every test here runs wherever that stack
exists; each skips on its own if it does not.

Run in the base.yaml env for full coverage:
  python -m pytest tests/test_combine_counts.py -v
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


@pytest.fixture
def mat():
    """matrix_utils, the module holding the summation kernel."""
    return pytest.importorskip("matrix_utils")


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


def test_sum_features_to_bbs_sums_per_bb(mat):
    x = np.array([[3.0, 0.0], [1.0, 0.0], [2.0, 5.0]])
    out = _dense(mat.sum_features_to_bbs(x, np.array([0, 0, 1]), 2))
    assert out.tolist() == [[4.0, 0.0], [2.0, 5.0]]


def test_same_orientation_bb_baf_no_cancellation(mat):
    # 3 LOH SNPs, all B-skewed (BAF 0.8) -> bin BAF stays 0.8
    b = np.array([[8.0], [7.0], [9.0]])
    tot = b + np.array([[2.0], [3.0], [1.0]])
    bins = np.array([0, 0, 0])
    baf = (
        _dense(mat.sum_features_to_bbs(b, bins, 1))[0, 0]
        / _dense(mat.sum_features_to_bbs(tot, bins, 1))[0, 0]
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


def test_loh_cancelled_when_merged_but_preserved_when_split(mat):
    # LOH with a mid-segment orientation flip: SNPs 0-2 B=0.8, SNPs 3-5 B=0.2
    b = np.array([[80], [80], [80], [20], [20], [20]], float)
    a = np.array([[20], [20], [20], [80], [80], [80]], float)
    tot = a + b

    # naive: one bin ignoring phase -> B and A cancel -> BAF 0.5 (LOH lost)
    merged = (
        _dense(mat.sum_features_to_bbs(b, np.zeros(6, int), 1))[0, 0]
        / _dense(mat.sum_features_to_bbs(tot, np.zeros(6, int), 1))[0, 0]
    )
    assert abs(merged - 0.5) < 1e-9

    # phase-aware: split at the detected flip -> each bin keeps the LOH (folded 0.2)
    pg = ph.detect_phase_flips(
        _snps(6), a, b, cluster_cols=["region_id", "PS"]
    ).to_numpy()
    bins, _ = pd.factorize(pg)
    k = int(bins.max()) + 1
    bsum = _dense(mat.sum_features_to_bbs(b, bins, k))[:, 0]
    tsum = _dense(mat.sum_features_to_bbs(tot, bins, k))[:, 0]
    folded = np.minimum(bsum / tsum, 1 - bsum / tsum)
    assert np.allclose(folded, 0.2, atol=1e-9)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"PASS {name}")
    print("all combine_counts unit tests passed")
