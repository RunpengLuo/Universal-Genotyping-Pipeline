"""scATAC helpers: per-bb fragment counts from 10x fragment files.

Fragment files are read in chunks because a single sample runs to hundreds of
millions of records; each fragment is counted once, by its midpoint.
"""

import logging

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from range_utils import assign_pos_to_range
from utils import add_chr_prefix


def count_atac_fragments_to_bbs(
    frag_files, reps, barcodes_full, bb_df, num_bbs, chunksize=5_000_000
):
    """Count deduped ATAC fragments per bb per cell from 10x fragment files.

    Each record of a 10x ``atac_fragments.tsv.gz`` is one deduplicated fragment
    (``chrom, start, end, barcode, readSupport``); the readSupport field is IGNORED.
    Every fragment is counted once, assigned to the bb containing its midpoint, so
    each observation sums to that cell's in-bb fragment count.

    Parameters
    ----------
    frag_files, reps : parallel lists
        ``frag_files[i]`` is the fragment file for replicate ``reps[i]``.
    barcodes_full : pd.DataFrame
        Columns ``REP_ID``, ``BARCODE`` (``BARCODE`` = ``"{raw}_{rep}"``) giving the observation
        order (identical to that assay's ``bb.*allele.npz`` observations).
    bb_df : pd.DataFrame
        bbs with ``#CHR``, ``START``, ``END`` (0-based half-open) and ``bb_id``.
    num_bbs : int
        Number of bbs (output features).

    Returns
    -------
    scipy.sparse.csr_matrix, shape ``(num_bbs, n_cells)``, dtype int32.
    """
    n_cells = len(barcodes_full)
    bc_rep = barcodes_full["REP_ID"].to_numpy().astype(str)
    bc_full = barcodes_full["BARCODE"].to_numpy().astype(str)
    # global observation index keyed by (rep, raw_barcode); strip the "_{rep}" suffix
    obs_of = {}
    for i in range(n_cells):
        rep, raw = bc_rep[i], bc_full[i]
        sfx = "_" + rep
        if raw.endswith(sfx):
            raw = raw[: -len(sfx)]
        obs_of[(rep, raw)] = i

    bb_all, obs_all = [], []
    for frag_file, rep in zip(frag_files, reps):
        rep_map = {raw: c for (r, raw), c in obs_of.items() if r == rep}
        if not rep_map or frag_file is None:
            continue
        n_frag = 0
        for chunk in pd.read_csv(
            frag_file,
            sep="\t",
            comment="#",
            header=None,
            usecols=[0, 1, 2, 3],
            names=["#CHR", "start", "end", "BC"],
            dtype={0: str, 1: np.int64, 2: np.int64, 3: str},
            chunksize=chunksize,
        ):
            obs_vals = chunk["BC"].map(rep_map).to_numpy()
            m = ~pd.isna(obs_vals)
            if not m.any():
                continue
            sub = chunk.loc[m]
            sub = sub.assign(**{"#CHR": add_chr_prefix(sub["#CHR"])})
            mid = (sub["start"].to_numpy() + sub["end"].to_numpy()) // 2
            frag = pd.DataFrame({"#CHR": sub["#CHR"].to_numpy(), "POS0": mid})
            frag = assign_pos_to_range(frag, bb_df, ref_id="bb_id", pos_col="POS0")
            keep = frag["bb_id"].notna().to_numpy()
            if not keep.any():
                continue
            bb_all.append(frag.loc[keep, "bb_id"].to_numpy().astype(np.int64))
            obs_all.append(obs_vals[m][keep].astype(np.int64))
            n_frag += int(keep.sum())
        logging.info(f"  ATAC {rep}: {n_frag} in-bb fragments counted")

    if bb_all:
        bbs = np.concatenate(bb_all)
        obs = np.concatenate(obs_all)
    else:
        bbs = np.zeros(0, dtype=np.int64)
        obs = np.zeros(0, dtype=np.int64)
    data = np.ones(len(bbs), dtype=np.int32)
    return csr_matrix((data, (bbs, obs)), shape=(num_bbs, n_cells), dtype=np.int32)
