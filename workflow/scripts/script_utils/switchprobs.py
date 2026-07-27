import numpy as np
import pandas as pd


def interp_cM_blocks(
    blocks: pd.DataFrame,
    snp_info: pd.DataFrame,
    genetic_map: pd.DataFrame,
    block_id_col: str = "bbc_id",
):
    """Interpolate centimorgan distances between consecutive blocks using a genetic map.

    Parameters
    ----------
    blocks : pd.DataFrame
        Block-level DataFrame with a block ID column and ``#CHR`` column.
    snp_info : pd.DataFrame
        SNP DataFrame with a matching block ID column and ``POS`` column.
    genetic_map : pd.DataFrame
        Genetic map with ``#CHR``, ``POS``, and ``cM`` columns.
    block_id_col : str
        Name of the block ID column in both *blocks* and *snp_info*.

    Returns
    -------
    np.ndarray
        Inter-block cM distances (first block per chromosome gets 0).
    """
    blocks = blocks.copy(deep=True)
    blocks["dist_cM"] = 0.0

    hb_pos = snp_info.groupby(block_id_col, sort=False)["POS"].agg(
        snp_start="min", snp_end="max"
    )
    blocks = blocks.join(hb_pos, on=block_id_col)

    genetic_map_chrs = genetic_map.groupby(by="#CHR", sort=False, observed=True)
    for ch, ch_blocks in blocks.groupby(by="#CHR", sort=False, observed=True):
        ch_map = genetic_map_chrs.get_group(ch)
        start_cMs = np.interp(
            ch_blocks["snp_start"].to_numpy(),
            ch_map["POS"].to_numpy(),
            ch_map["cM"].to_numpy(),
        )
        end_cMs = np.interp(
            ch_blocks["snp_end"].to_numpy(),
            ch_map["POS"].to_numpy(),
            ch_map["cM"].to_numpy(),
        )

        dist_cM = np.zeros(len(ch_blocks), dtype=np.float32)
        dist_cM[1:] = start_cMs[1:] - end_cMs[:-1]
        blocks.loc[ch_blocks.index, "dist_cM"] = np.maximum(dist_cM, 0.0)
    return blocks["dist_cM"].to_numpy()


def estimate_switchprobs_cM(dist_cms: np.ndarray, nu=1, min_switchprob=1e-6):
    """Convert cM distances to haplotype switch probabilities using the Haldane mapping function.

    Computes ``(1 - exp(-2 * nu * d)) / 2``, clipped at *min_switchprob*.

    Parameters
    ----------
    dist_cms : np.ndarray
        Inter-SNP or inter-block centimorgan distances.
    nu : float
        Scaling factor for the Haldane function.
    min_switchprob : float
        Minimum switch probability (floor).

    Returns
    -------
    np.ndarray
        Switch probabilities.
    """
    switchprobs = (1 - np.exp(-2 * nu * dist_cms)) / 2.0
    return np.clip(switchprobs, a_min=min_switchprob, a_max=None)


def estimate_switchprobs_PS(blocks: pd.DataFrame, switchprob_ps=0.05):
    """Assign switch probabilities based on PS phaseset membership.

    Within the same phaseset, the probability is *switchprob_ps*; across
    different phasesets it is approximately 0.5.

    Parameters
    ----------
    blocks : pd.DataFrame
        DataFrame with a ``PS`` column indicating phaseset IDs.
    switchprob_ps : float
        Switch probability within the same phaseset.

    Returns
    -------
    np.ndarray
        Switch probabilities per block.
    """
    switch_bias = 1e-4
    same_block = blocks["PS"] == blocks["PS"].shift(1).fillna(False)
    switchprobs = np.where(
        same_block,
        switchprob_ps,  # within same PS phase block
        0.5 - switch_bias,  # across different PS phase block
    )
    return switchprobs
