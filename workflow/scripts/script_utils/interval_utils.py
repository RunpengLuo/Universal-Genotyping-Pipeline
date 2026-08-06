"""Assign query positions or intervals to a reference interval set.

Two public entry points, both keyed on ``#CHR`` and both returning *qry* with a
``ref_id`` column added (``pd.NA`` where nothing overlaps):

- ``assign_pos_to_range`` - a query POSITION lands in the interval containing it.
- ``assign_interval_to_range`` - a query INTERVAL takes the reference interval it
  overlaps most.

``assign_all_features`` is the many-hit variant of the first: every overlapping id,
joined. All three share ``_searchsorted_assign``; the scan fallback runs only for a
chromosome whose reference intervals overlap each other.

Reference intervals are 0-based half-open ``[START, END)``.
"""

import heapq

import numpy as np
import pandas as pd


def _searchsorted_assign(starts, ends, positions):
    """Index of the non-overlapping, start-sorted interval containing each position.

    ``starts``/``ends`` are 0-based half-open (``START <= pos < END``), sorted by
    ``starts`` with no overlaps. Returns ``(idx, valid)``: ``idx[k]`` is the interval
    index for ``positions[k]`` where ``valid[k]``, undefined otherwise.
    """
    if len(starts) == 0:
        n = len(positions)
        return np.zeros(n, dtype=np.int64), np.zeros(n, dtype=bool)
    idx = np.searchsorted(starts, positions, side="right") - 1
    safe_idx = idx.clip(min=0)
    valid = (idx >= 0) & (positions < ends[safe_idx])
    return idx, valid


def _sorted_ref(ref_chrom, ref_id):
    """One chromosome's reference intervals as start-sorted ``(starts, ends, ids)``."""
    starts = ref_chrom["START"].to_numpy()
    ends = ref_chrom["END"].to_numpy()
    ids = ref_chrom[ref_id].to_numpy()
    order = np.argsort(starts, kind="stable")
    return starts[order], ends[order], ids[order]


def _interval_tiers(starts, ends):
    """Partition intervals into the minimum number of non-overlapping tiers.

    Greedy earliest-finishing assignment (sort by start, reuse the tier whose last
    END fits before the next START, else open a new tier). Tier count equals the max
    overlap depth. Members within a tier come out start-sorted. Returns a list of
    index arrays into the input.
    """
    order = np.argsort(starts, kind="stable")
    heap = []  # (last_end, tier_id)
    tier_members = []
    for i in order:
        s, e = int(starts[i]), int(ends[i])
        if heap and heap[0][0] <= s:
            _, t = heapq.heappop(heap)
            tier_members[t].append(i)
            heapq.heappush(heap, (e, t))
        else:
            t = len(tier_members)
            tier_members.append([i])
            heapq.heappush(heap, (e, t))
    return [np.array(m, dtype=np.int64) for m in tier_members]


def assign_pos_to_range(qry, ref, ref_id="region_id", pos_col="POS0"):
    """Assign each query position to the reference interval containing it.

    Vectorized per chromosome via ``np.searchsorted``; a chromosome whose reference
    intervals overlap each other falls back to a scan that keeps the first hit.

    Args:
        qry: Query frame with ``#CHR`` and *pos_col*. Modified in place.
        ref: Reference intervals with ``#CHR``, ``START``, ``END`` and *ref_id*.
        ref_id: Reference identifier column, added to *qry*.
        pos_col: 0-based position column of *qry*.

    Returns:
        *qry* with *ref_id* added; ``pd.NA`` where the position is in no interval.
    """

    def _first_overlap_scan(starts, ends, ids, positions):
        """First overlapping interval id per position, for a chromosome with overlaps.

        ``starts``/``ends``/``ids`` must be start-sorted. Returns an object array holding
        ``pd.NA`` where a position overlaps nothing.
        """
        right_bounds = np.searchsorted(starts, positions, side="right")
        out = np.empty(len(positions), dtype=object)
        out[:] = pd.NA
        for i, pos in enumerate(positions):
            cands = slice(0, right_bounds[i])
            mask = ends[cands] > pos
            if mask.any():
                out[i] = ids[cands][mask][0]
        return out

    qry[ref_id] = pd.NA
    for chrom, ref_c in ref.groupby("#CHR", sort=False):
        qmask = (qry["#CHR"] == chrom).to_numpy()
        if not qmask.any():
            continue
        positions = qry.loc[qmask, pos_col].to_numpy()
        qidx = qry.index[qmask]
        starts, ends, ids = _sorted_ref(ref_c, ref_id)

        if len(starts) > 1 and np.any(starts[1:] < ends[:-1]):
            qry.loc[qidx, ref_id] = _first_overlap_scan(starts, ends, ids, positions)
        else:
            idx, valid = _searchsorted_assign(starts, ends, positions)
            qry.loc[qidx[valid], ref_id] = ids[idx[valid]]
    return qry


def assign_interval_to_range(qry, ref, ref_id):
    """Assign each query interval to the reference interval it overlaps most.

    Args:
        qry: Query intervals with ``#CHR``, ``START``, ``END``. Copied, not modified.
        ref: Reference intervals with ``#CHR``, ``START``, ``END`` and *ref_id*.
        ref_id: Reference identifier column, added to the returned frame.

    Returns:
        A copy of *qry* with *ref_id* added; ``pd.NA`` where nothing overlaps.
    """
    qry = qry.copy()
    qry[ref_id] = pd.NA
    for chrom, ref_c in ref.groupby("#CHR", sort=False):
        qmask = (qry["#CHR"] == chrom).to_numpy()
        if not qmask.any():
            continue
        q_starts = qry.loc[qmask, "START"].to_numpy()
        q_ends = qry.loc[qmask, "END"].to_numpy()
        starts, ends, ids = _sorted_ref(ref_c, ref_id)

        right_bounds = np.searchsorted(starts, q_ends, side="left")
        best = np.empty(len(q_starts), dtype=object)
        best[:] = pd.NA
        for i in range(len(q_starts)):
            qs, qe = q_starts[i], q_ends[i]
            cands = slice(0, right_bounds[i])
            mask = ends[cands] > qs
            if not mask.any():
                continue
            c_starts, c_ends = starts[cands][mask], ends[cands][mask]
            overlap = np.minimum(qe, c_ends) - np.maximum(qs, c_starts)
            best[i] = ids[cands][mask][np.argmax(overlap)]
        qry.loc[qmask, ref_id] = best
    return qry


def overlaps_any_range(qry, ref, pos_col="POS0"):
    """Boolean mask: does each query position fall in any reference interval?

    The membership-only variant of ``assign_pos_to_range`` - no id is carried, so
    overlapping references need no tie-break and every chromosome takes the
    vectorized path (one ``_searchsorted_assign`` per overlap tier).

    Args:
        qry: Query frame with ``#CHR`` and *pos_col*. Not modified.
        ref: Reference intervals with ``#CHR``, ``START``, ``END``.
        pos_col: 0-based position column of *qry*.

    Returns:
        A bool array aligned to *qry* rows.
    """
    hit = np.zeros(len(qry), dtype=bool)
    if len(ref) == 0:
        return hit
    for chrom, ref_c in ref.groupby("#CHR", sort=False):
        qmask = (qry["#CHR"] == chrom).to_numpy()
        if not qmask.any():
            continue
        positions = qry.loc[qmask, pos_col].to_numpy()
        starts = ref_c["START"].to_numpy()
        ends = ref_c["END"].to_numpy()
        found = np.zeros(len(positions), dtype=bool)
        for tier in _interval_tiers(starts, ends):
            _, valid = _searchsorted_assign(starts[tier], ends[tier], positions)
            found |= valid
        hit[qmask] = found
    return hit


def assign_all_features(
    snps, ref, id_col, pos_col="POS0", sep=";", default="intergenic"
):
    """All overlapping ``ref`` ids per SNP, ``sep``-joined (``default`` when none).

    Vectorized: ``ref`` intervals are split into non-overlapping tiers so each tier
    uses the ``_searchsorted_assign`` fast path. The first (densest) tier is assigned
    fully vectorized; only the rare SNPs that also hit a higher tier are joined.
    Returns a Series aligned to ``snps.index``.
    """
    result = pd.Series(default, index=snps.index, dtype=object)
    if len(ref) == 0:
        return result
    for chrom, ref_c in ref.groupby("#CHR", sort=False):
        qmask = (snps["#CHR"] == chrom).to_numpy()
        if not qmask.any():
            continue
        positions = snps.loc[qmask, pos_col].to_numpy()
        qidx = snps.index[qmask]
        starts = ref_c["START"].to_numpy()
        ends = ref_c["END"].to_numpy()
        ids = ref_c[id_col].to_numpy().astype(str)
        tiers = _interval_tiers(starts, ends)

        t0 = tiers[0]
        idx0, valid0 = _searchsorted_assign(starts[t0], ends[t0], positions)
        joined = np.where(valid0, ids[t0][idx0.clip(min=0)], "").astype(object)
        for tier in tiers[1:]:
            idx, valid = _searchsorted_assign(starts[tier], ends[tier], positions)
            for k in np.nonzero(valid)[0]:
                gid = ids[tier][idx[k]]
                joined[k] = gid if joined[k] == "" else f"{joined[k]}{sep}{gid}"
        joined = np.where(joined == "", default, joined)
        result.loc[qidx] = joined
    return result
