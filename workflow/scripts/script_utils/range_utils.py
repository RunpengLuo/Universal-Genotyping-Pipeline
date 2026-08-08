"""Assign query positions or ranges to a reference range set.

Coordinates are 0-based throughout: a pos is a 0-based offset (``POS0``, never the
1-based VCF ``POS``) and a range is a 0-based half-open ``[START, END)`` pair, on either
side. Both are asserted on every call. The 1-based boundaries are named in
``docs/DEVELOPER.md``.

Three frame-level entry points, all keyed on ``#CHR``, all returning
``(annotated_qry, na_idx)``:

- ``assign_pos_to_range`` - a query POSITION lands in the range containing it.
- ``assign_pos_to_range_ovlp`` - the many-hit variant: every overlapping id, joined.
- ``assign_range_to_range`` - a query RANGE, by ``rule``: largest overlap, midpoint, or
  contained (both ends in one reference id).

They never mutate *qry*; the out-column is built as one positional array, so a duplicate
index is harmless. ``na_idx`` holds the positional indices of the unassigned rows, which
is what makes ``dropna=True`` usable when a parallel matrix has to be subset the same way.

``overlaps_any_range`` is the boolean variant and ``merge_ranges_to_clusters`` the
array-level one: its ranges index into an ordered sequence rather than a coordinate.
"""

import heapq
import logging

import numpy as np
import pandas as pd

RULES = ("max_overlap", "midpoint", "contained")


def _searchsorted_assign(starts, ends, positions):
    """Index of the non-overlapping, start-sorted range containing each position.

    ``starts``/``ends`` are 0-based half-open (``START <= pos < END``), sorted by
    ``starts`` with no overlaps. Returns ``(idx, valid)``: ``idx[k]`` is the range
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
    """One chromosome's reference ranges as start-sorted ``(starts, ends, ids)``."""
    starts = ref_chrom["START"].to_numpy()
    ends = ref_chrom["END"].to_numpy()
    ids = ref_chrom[ref_id].to_numpy()
    order = np.argsort(starts, kind="stable")
    return starts[order], ends[order], ids[order]


def _range_clusters(starts, ends):
    """Partition ranges into the minimum number of non-overlapping clusters.

    Greedy earliest-finishing assignment (sort by start, reuse the cluster whose last
    END fits before the next START, else open a new cluster). Cluster count equals the
    max overlap depth. Members within a cluster come out start-sorted. Returns a list
    of index arrays into the input.
    """
    order = np.argsort(starts, kind="stable")
    heap = []  # (last_end, cluster_id)
    cluster_members = []
    for i in order:
        s, e = int(starts[i]), int(ends[i])
        if heap and heap[0][0] <= s:
            _, c = heapq.heappop(heap)
            cluster_members[c].append(i)
            heapq.heappush(heap, (e, c))
        else:
            c = len(cluster_members)
            cluster_members.append([i])
            heapq.heappush(heap, (e, c))
    return [np.array(m, dtype=np.int64) for m in cluster_members]


def _first_overlap_scan(starts, ends, ids, positions):
    """First overlapping range id per position, for a chromosome with overlaps.

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


def _check_ranges(df, what):
    """Assert *df* holds 0-based half-open ranges."""
    assert (df["END"].to_numpy() > df["START"].to_numpy()).all(), (
        f"{what}, range is not 0-based [START, END)"
    )


def _check_pos_col(pos_col):
    """Assert *pos_col* is not the 1-based VCF column."""
    assert pos_col != "POS", "pos_col, `POS` is 1-based; pass `POS0`"


def _na_array(n):
    """Object array of *n* ``pd.NA``."""
    out = np.empty(n, dtype=object)
    out[:] = pd.NA
    return out


def _pos_values(chroms, positions, ref, ref_id):
    """Reference id per position, as a positional object array.

    Vectorized per chromosome via ``np.searchsorted``; a chromosome whose reference
    ranges overlap each other falls back to a scan that keeps the first hit.
    """
    out = _na_array(len(chroms))
    for chrom, ref_c in ref.groupby("#CHR", sort=False):
        qmask = chroms == chrom
        if not qmask.any():
            continue
        pos_c = positions[qmask]
        starts, ends, ids = _sorted_ref(ref_c, ref_id)
        # start-sorted, so an adjacent-pair test detects any overlap in the set
        if len(starts) > 1 and np.any(starts[1:] < ends[:-1]):
            out[qmask] = _first_overlap_scan(starts, ends, ids, pos_c)
        else:
            idx, valid = _searchsorted_assign(starts, ends, pos_c)
            sub = _na_array(len(pos_c))
            sub[valid] = ids[idx[valid]]
            out[qmask] = sub
    return out


def _max_overlap_values(qry, ref, ref_id):
    """Reference id of the largest overlap per query range, as a positional array."""
    chroms = qry["#CHR"].to_numpy()
    qs_all = qry["START"].to_numpy()
    qe_all = qry["END"].to_numpy()
    out = _na_array(len(qry))
    for chrom, ref_c in ref.groupby("#CHR", sort=False):
        qmask = chroms == chrom
        if not qmask.any():
            continue
        q_starts, q_ends = qs_all[qmask], qe_all[qmask]
        starts, ends, ids = _sorted_ref(ref_c, ref_id)
        right_bounds = np.searchsorted(starts, q_ends, side="left")
        best = _na_array(len(q_starts))
        for i in range(len(q_starts)):
            qs, qe = q_starts[i], q_ends[i]
            cands = slice(0, right_bounds[i])
            mask = ends[cands] > qs
            if not mask.any():
                continue
            c_starts, c_ends = starts[cands][mask], ends[cands][mask]
            overlap = np.minimum(qe, c_ends) - np.maximum(qs, c_starts)
            best[i] = ids[cands][mask][np.argmax(overlap)]
        out[qmask] = best
    return out


def _finalize(qry, values, out_col, ref, ref_id, fillna, dropna, cast=True):
    """Attach *values* to a copy of *qry*, then apply *fillna* or *dropna*.

    Args:
        qry: Query frame. Copied, never modified.
        values: Positional object array of ids, ``pd.NA`` where unassigned.
        out_col: Column to write *values* into.
        ref: Reference frame, for the *dropna* dtype cast.
        ref_id: Identifier column of *ref*, for the same.
        fillna: Value for unassigned rows; no dtype cast is applied.
        dropna: Drop unassigned rows, reindex from 0, and cast *out_col*.
        cast: Whether *dropna* casts to ``ref[ref_id].dtype``. False for a
            joined-string out-column, whose dtype is not the reference's.

    Returns:
        ``(qry, na_idx)`` - the annotated frame and the positional indices of the
        unassigned rows of the INPUT frame.
    """
    assert not (fillna is not None and dropna), (
        "assign, fillna and dropna are exclusive"
    )
    na = pd.isna(values)
    na_idx = np.flatnonzero(na)
    if fillna is not None and len(na_idx):
        values = values.copy()
        values[na] = fillna

    qry = qry.copy()
    qry[out_col] = values
    if dropna:
        if len(na_idx):
            logging.info(f"{out_col}: dropped {len(na_idx)}/{len(qry)} unassigned rows")
            qry = qry.loc[~na]
        qry = qry.reset_index(drop=True)
        if cast:
            qry[out_col] = qry[out_col].astype(ref[ref_id].dtype)
    return qry, na_idx


def assign_pos_to_range(
    qry,
    ref,
    ref_id="region_id",
    pos_col="POS0",
    out_col=None,
    fillna=None,
    dropna=False,
):
    """Assign each query position to the reference range containing it.

    Args:
        qry: Query frame with ``#CHR`` and *pos_col*. Copied, not modified.
        ref: Reference ranges with ``#CHR``, ``START``, ``END`` and *ref_id*.
        ref_id: Reference identifier column of *ref*.
        pos_col: 0-based position column of *qry*.
        out_col: Column added to *qry*; defaults to *ref_id*.
        fillna: Value for positions in no range, instead of ``pd.NA``.
        dropna: Drop those positions and cast *out_col* to the reference dtype.

    Returns:
        ``(qry, na_idx)``; see :func:`_finalize`.
    """
    _check_pos_col(pos_col)
    _check_ranges(ref, "ref")
    out_col = out_col or ref_id
    values = _pos_values(qry["#CHR"].to_numpy(), qry[pos_col].to_numpy(), ref, ref_id)
    return _finalize(qry, values, out_col, ref, ref_id, fillna, dropna)


def assign_range_to_range(
    qry, ref, ref_id, rule="max_overlap", out_col=None, fillna=None, dropna=False
):
    """Assign each query range to a reference range, by *rule*.

    ``contained`` compares reference id VALUES, not reference rows, so a query spanning
    the gap between two rows that share an id (a segment split by a blacklist hole) is
    contained, not straddling. ``midpoint`` and ``contained`` run the vectorized
    position kernel; only ``max_overlap`` walks the queries one at a time.

    Args:
        qry: Query ranges with ``#CHR``, ``START``, ``END``. Copied, not modified.
        ref: Reference ranges with ``#CHR``, ``START``, ``END`` and *ref_id*.
        ref_id: Reference identifier column of *ref*.
        rule: ``max_overlap`` (largest shared span), ``midpoint`` (the range holding
            ``(START + END) // 2``), or ``contained`` (an id only when ``START`` and
            ``END - 1`` land in the same one).
        out_col: Column added to *qry*; defaults to *ref_id*.
        fillna: Value for ranges with no assignment, instead of ``pd.NA``.
        dropna: Drop those ranges and cast *out_col* to the reference dtype.

    Returns:
        ``(qry, na_idx)``; see :func:`_finalize`.
    """
    assert rule in RULES, f"rule, unknown: {rule}"
    _check_ranges(ref, "ref")
    _check_ranges(qry, "qry")
    out_col = out_col or ref_id

    chroms = qry["#CHR"].to_numpy()
    starts, ends = qry["START"].to_numpy(), qry["END"].to_numpy()
    if rule == "midpoint":
        values = _pos_values(chroms, (starts + ends) // 2, ref, ref_id)
    elif rule == "contained":
        at_start = _pos_values(chroms, starts, ref, ref_id)
        at_end = _pos_values(chroms, ends - 1, ref, ref_id)
        both = ~pd.isna(at_start) & ~pd.isna(at_end)
        same = np.zeros(len(qry), dtype=bool)
        same[both] = at_start[both] == at_end[both]
        values = _na_array(len(qry))
        values[same] = at_start[same]
    else:
        values = _max_overlap_values(qry, ref, ref_id)
    return _finalize(qry, values, out_col, ref, ref_id, fillna, dropna)


def assign_pos_to_range_ovlp(
    qry,
    ref,
    ref_id,
    pos_col="POS0",
    out_col=None,
    sep=";",
    fillna=None,
    dropna=False,
):
    """Assign each query position EVERY overlapping reference id, *sep*-joined.

    The many-hit variant of :func:`assign_pos_to_range`, for a reference whose ranges
    nest (GTF genes). Vectorized: *ref* is split into non-overlapping clusters so each
    takes the ``_searchsorted_assign`` fast path; the first (densest) cluster is assigned
    wholesale and only positions that also hit a later cluster are joined.

    *out_col* holds a joined string, not a reference id, so *dropna* does not cast it.

    Args:
        qry: Query frame with ``#CHR`` and *pos_col*. Copied, not modified.
        ref: Reference ranges with ``#CHR``, ``START``, ``END`` and *ref_id*.
        ref_id: Reference identifier column of *ref*; values are stringified.
        pos_col: 0-based position column of *qry*.
        out_col: Column added to *qry*; defaults to *ref_id*.
        sep: Separator between the ids of one position.
        fillna: Value for positions in no range (e.g. ``"intergenic"``).
        dropna: Drop those positions.

    Returns:
        ``(qry, na_idx)``; see :func:`_finalize`.
    """
    _check_pos_col(pos_col)
    _check_ranges(ref, "ref")
    out_col = out_col or ref_id

    values = _na_array(len(qry))
    chroms = qry["#CHR"].to_numpy()
    pos_all = qry[pos_col].to_numpy()
    for chrom, ref_c in ref.groupby("#CHR", sort=False):
        qmask = chroms == chrom
        if not qmask.any():
            continue
        pos_c = pos_all[qmask]
        starts = ref_c["START"].to_numpy()
        ends = ref_c["END"].to_numpy()
        ids = ref_c[ref_id].to_numpy().astype(str)
        clusters = _range_clusters(starts, ends)

        c0 = clusters[0]
        idx0, valid0 = _searchsorted_assign(starts[c0], ends[c0], pos_c)
        joined = np.where(valid0, ids[c0][idx0.clip(min=0)], "").astype(object)
        for cluster in clusters[1:]:
            idx, valid = _searchsorted_assign(starts[cluster], ends[cluster], pos_c)
            for k in np.nonzero(valid)[0]:
                rid = ids[cluster][idx[k]]
                joined[k] = rid if joined[k] == "" else f"{joined[k]}{sep}{rid}"
        sub = _na_array(len(pos_c))
        hit = joined != ""
        sub[hit] = joined[hit]
        values[qmask] = sub
    return _finalize(qry, values, out_col, ref, ref_id, fillna, dropna, cast=False)


def overlaps_any_range(qry, ref, pos_col="POS0"):
    """Boolean mask: does each query position fall in any reference range?

    The membership-only variant of :func:`assign_pos_to_range` - no id is carried, so
    overlapping references need no tie-break and every chromosome takes the
    vectorized path (one ``_searchsorted_assign`` per overlap cluster).

    Args:
        qry: Query frame with ``#CHR`` and *pos_col*. Not modified.
        ref: Reference ranges with ``#CHR``, ``START``, ``END``.
        pos_col: 0-based position column of *qry*.

    Returns:
        A bool array aligned to the *qry* features.
    """
    _check_pos_col(pos_col)
    hit = np.zeros(len(qry), dtype=bool)
    if len(ref) == 0:
        return hit
    _check_ranges(ref, "ref")
    chroms = qry["#CHR"].to_numpy()
    pos_all = qry[pos_col].to_numpy()
    for chrom, ref_c in ref.groupby("#CHR", sort=False):
        qmask = chroms == chrom
        if not qmask.any():
            continue
        pos_c = pos_all[qmask]
        starts = ref_c["START"].to_numpy()
        ends = ref_c["END"].to_numpy()
        found = np.zeros(len(pos_c), dtype=bool)
        for cluster in _range_clusters(starts, ends):
            _, valid = _searchsorted_assign(starts[cluster], ends[cluster], pos_c)
            found |= valid
        hit[qmask] = found
    return hit


def merge_ranges_to_clusters(n_items, ranges):
    """Cluster id per ordered item, merging every item a range spans into one cluster.

    Here a range is a ``[lo, hi)`` pair of INDICES into the item ordering, not a
    coordinate. Every boundary interior to a range is made non-cuttable, so ranges
    sharing an item merge transitively and an item in no range is its own cluster.
    A caller uses this to keep an indivisible span whole - a gene over the fixed bins
    holding its SNPs, so no bb splits it.

    Args:
        n_items: Number of ordered items.
        ranges: Iterable of ``(lo, hi)`` half-open index pairs.

    Returns:
        int64 array of length *n_items*, run-length contiguous; the boundary between
        items ``i - 1`` and ``i`` is a cut point iff ``labels[i] != labels[i - 1]``.
    """
    if n_items == 0:
        return np.zeros(0, dtype=np.int64)
    blocked = np.zeros(n_items - 1, dtype=bool)  # boundary between item i and i+1
    for lo, hi in ranges:
        if hi - 1 > lo:
            blocked[lo : hi - 1] = True
    labels = np.empty(n_items, dtype=np.int64)
    labels[0] = 0
    if n_items > 1:
        labels[1:] = np.cumsum(~blocked)
    return labels
