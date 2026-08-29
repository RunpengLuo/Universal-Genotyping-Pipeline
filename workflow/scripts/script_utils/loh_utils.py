"""Detect clonal-LOH regions from het-SNP density along the genome.

Last update: 2026-08-28

Under clonal LOH one haplotype is gone, so every germline het there is called hom and the
observed het-SNP density collapses. Density is the statistic, not the per-site allele
fraction: a run of hom-looking sites happens for panel and mappability reasons inside
perfectly heterozygous sequence, and only an aggregate over a whole tile can tell the two
apart.

Functions:
- build_loh_tiles: windows -> fixed-width tiles carrying het counts and assayable span
- fit_snp_rate: negative-binomial MLE of the neutral rate and its overdispersion
- viterbi_loh: two-state decode, neutral against a depleted rate
- loh_intervals: runs of the LOH state as a BED3 frame

References:
- Numbat ``detect_clonal_loh`` / ``fit_snp_rate`` / ``viterbi_loh``
  (kharchenkolab/numbat, R/utils.R and R/hmm.R), doi:10.1038/s41587-022-01468-y.
  Two differences, both deliberate: the transition is per bp rather than per unit, as
  everywhere else in this pipeline, and Numbat's second emission (expression at
  ``phi = 0.5``, making it a clonal-DELETION detector) is dropped - copy-neutral LOH and
  LOH on a gained arm are common in bulk DNA and would be discarded by it.
"""

import logging

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import nbinom


def _nb_logpmf(x, mu, size):
    """``log NB(x; mu, size)``, guarding the zero-mean tile."""
    mu = np.maximum(mu, 1e-9)
    return nbinom.logpmf(x, size, size / (size + mu))


def build_loh_tiles(bins, snp_bin_ids, assayable, tile_size, cluster_cols=("seg_id",)):
    """Aggregate windows into fixed-width tiles carrying het counts and assayable span.

    Tiles are cut on *cluster_cols* as well as on *tile_size*, so none spans a segment
    boundary; the windows themselves never do, since they are tiled from ``segment.bed``.

    Args:
        bins: Window frame with ``#CHR``, ``START``, ``END``, ``region_id`` and
            *cluster_cols*, in genomic order.
        snp_bin_ids: Window index of every het SNP; one entry per SNP.
        assayable: Boolean per window. A window that is not assayable contributes its
            het count but no exposure, so a coverage hole cannot read as depletion.
        tile_size: Tile width in bp.
        cluster_cols: Columns a tile may not span.

    Returns:
        Tile frame with ``#CHR``, ``START``, ``END``, ``region_id``, ``n_het``,
        ``exposure_bp`` and ``mid``, sorted by position.
    """
    n_het = np.bincount(np.asarray(snp_bin_ids), minlength=len(bins))
    span = (bins["END"] - bins["START"]).to_numpy()
    frame = pd.DataFrame(
        {
            "#CHR": bins["#CHR"].to_numpy(),
            "START": bins["START"].to_numpy(),
            "END": bins["END"].to_numpy(),
            "region_id": bins["region_id"].to_numpy(),
            "n_het": n_het,
            "exposure_bp": np.where(np.asarray(assayable), span, 0),
            "tile": bins["START"].to_numpy() // int(tile_size),
        }
    )
    keys = ["#CHR", *cluster_cols, "tile"]
    for col in cluster_cols:
        frame[col] = bins[col].to_numpy()
    tiles = (
        frame.groupby(keys, sort=False)
        .agg(
            START=("START", "min"),
            END=("END", "max"),
            region_id=("region_id", "first"),
            n_het=("n_het", "sum"),
            exposure_bp=("exposure_bp", "sum"),
        )
        .reset_index()
        .drop(columns=["tile", *cluster_cols])
        .sort_values(["#CHR", "START"], kind="mergesort")
        .reset_index(drop=True)
    )
    tiles["mid"] = (tiles["START"] + tiles["END"]) // 2
    logging.info(
        f"LOH tiles: {len(tiles)} of {tile_size / 1e3:.0f} kb over {len(bins)} windows, "
        f"{int((tiles['exposure_bp'] == 0).sum())} with no assayable span"
    )
    return tiles


def fit_snp_rate(n_het, exposure_mb, init=(500.0, 2.0)):
    """MLE of the het rate per Mb and its negative-binomial overdispersion.

    Fitted over EVERY tile, LOH included, as Numbat does. Refitting on the tiles a first
    decode called neutral is self-reinforcing and must not be added: truncating the
    neutral set raises ``rate`` and tightens ``size``, which calls more LOH, which
    truncates further. The LOH tiles' near-zero counts are what give ``size`` the heavy
    tail that keeps a genuinely quiet neutral tile from flipping.

    Args:
        n_het: Het count per tile.
        exposure_mb: Assayable span per tile, in Mb; zero-exposure tiles are dropped.
        init: Starting ``(rate, size)``.

    Returns:
        ``(rate, size)``; ``size`` is the NB shape, small values meaning clumpier counts
        than Poisson.
    """
    keep = np.asarray(exposure_mb) > 0
    x, e = np.asarray(n_het)[keep], np.asarray(exposure_mb)[keep]
    assert x.size, "fit_snp_rate: no tile has assayable span"

    def nll(par):
        rate, size = par
        if rate <= 0 or size <= 0:
            return np.inf
        return -_nb_logpmf(x, rate * e, size).sum()

    rate, size = minimize(nll, init, method="Nelder-Mead").x
    logging.info(
        f"fit_snp_rate: {rate:.1f} het/Mb, NB size {size:.2f}, over {x.size} tiles "
        f"(observed median {np.median(x / e):.1f} het/Mb)"
    )
    return float(rate), float(size)


def viterbi_loh(n_het, exposure_mb, mid, group_id, rate_ref, rate_loh, size, t):
    """Two-state Viterbi over the tiles: neutral rate against a depleted one.

    The chain restarts at every group, so it never runs across a centromere. Transitions
    are Poisson per bp on the gap between tile midpoints, matching the rest of this
    pipeline; a tile with no assayable span emits nothing and only propagates.

    Args:
        n_het: Het count per tile.
        exposure_mb: Assayable span per tile, in Mb.
        mid: Tile midpoint, for the transition distance.
        group_id: Chain group per tile (the chromosome arm).
        rate_ref, rate_loh: Het rate per Mb under each state.
        size: NB overdispersion, shared by both states.
        t: Breakpoint rate per bp.

    Returns:
        Boolean array, True where the MAP state is LOH.
    """
    x = np.asarray(n_het)
    e = np.asarray(exposure_mb, dtype=float)
    mid = np.asarray(mid, dtype=np.int64)
    group_id = np.asarray(group_id)
    n = x.size

    log_e = np.zeros((n, 2))
    obs = e > 0
    log_e[obs, 0] = _nb_logpmf(x[obs], rate_ref * e[obs], size)
    log_e[obs, 1] = _nb_logpmf(x[obs], rate_loh * e[obs], size)

    nu = np.full((n, 2), -np.inf)
    ptr = np.zeros((n, 2), dtype=np.int8)
    nu[0] = np.log(0.5) + log_e[0]
    for i in range(1, n):
        if group_id[i] != group_id[i - 1]:
            nu[i] = np.log(0.5) + log_e[i]
            continue
        d = max(int(mid[i] - mid[i - 1]), 1)
        log_stay = -t * d
        log_switch = np.log1p(-np.exp(log_stay)) if log_stay > -30 else np.log(t * d)
        for j in range(2):
            cand = nu[i - 1] + np.array(
                [log_stay if k == j else log_switch for k in range(2)]
            )
            ptr[i, j] = int(np.argmax(cand))
            nu[i, j] = cand[ptr[i, j]] + log_e[i, j]

    state = np.zeros(n, dtype=np.int8)
    state[-1] = int(np.argmax(nu[-1]))
    for i in range(n - 2, -1, -1):
        state[i] = (
            int(np.argmax(nu[i]))
            if group_id[i] != group_id[i + 1]
            else ptr[i + 1, state[i + 1]]
        )
    return state.astype(bool)


def loh_intervals(tiles, is_loh):
    """Runs of the LOH state as a BED3 frame.

    Consecutive LOH tiles within one contig merge. Nothing is merged across a non-LOH
    gap - a gap is evidence against, and bridging it is what let short spurious runs
    chain together in the earlier design.

    Every run is reported. A run is at minimum one tile, so ``tile_size`` already floors
    the interval length, and the transition prior in :func:`viterbi_loh` is where short
    runs are resisted. A length filter on top of both cuts real and spurious regions
    alike, so none is applied here, as in Numbat.

    Args:
        tiles: Tile frame from :func:`build_loh_tiles`, position-sorted.
        is_loh: Boolean per tile.

    Returns:
        DataFrame with ``#CHR``, ``START``, ``END``, 0-based half-open.
    """
    is_loh = np.asarray(is_loh, dtype=bool)
    if not is_loh.any():
        logging.info("clonal-LOH intervals: none")
        return pd.DataFrame(columns=["#CHR", "START", "END"])

    chrom = tiles["#CHR"].to_numpy()
    new_run = np.r_[True, (chrom[1:] != chrom[:-1]) | (is_loh[1:] != is_loh[:-1])]
    run_id = np.cumsum(new_run)
    df = (
        pd.DataFrame(
            {
                "#CHR": chrom,
                "START": tiles["START"].to_numpy(),
                "END": tiles["END"].to_numpy(),
                "run": run_id,
                "is_loh": is_loh,
            }
        )
        .groupby("run", sort=True)
        .agg(
            **{
                "#CHR": ("#CHR", "first"),
                "START": ("START", "min"),
                "END": ("END", "max"),
                "is_loh": ("is_loh", "first"),
            }
        )
    )
    df = df[df["is_loh"]].drop(columns="is_loh").reset_index(drop=True)
    spans = (df["END"] - df["START"]).to_numpy()
    logging.info(
        f"clonal-LOH intervals: {len(df)}, {int(spans.sum()) / 1e6:.1f} Mbp, "
        f"shortest {int(spans.min()) / 1e3:.0f} kb"
    )
    return df
