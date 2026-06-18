"""SNP-informed adaptive binning across all bulk assays + depth aggregation + RDR.

All bulk assays present in the run (e.g. bulkWGS + bulkWGS-lr) are segmented on ONE
shared bin grid: ``adaptive_segmentation`` requires ``min_snp_reads`` in every tumor
sample, so stacking all assays' tumor columns yields bins that jointly satisfy every
bulk sample. Allele counts are aggregated per bin across all samples; read depth and
RDR are computed per assay, normalizing each assay's tumors by that assay's own normal.

Inputs are lists (one entry per bulk assay, index-aligned to ``params.bulk_assays``).
With a single bulk assay the lists have length 1. Outputs are written flat under
``bb_dir/`` (e.g. ``bb_dir/bb.tsv.gz``); matrix columns are the bulk samples.
"""

import os
import logging

t = int(getattr(snakemake, "threads", 1))
os.environ["OMP_NUM_THREADS"] = str(t)
os.environ["OPENBLAS_NUM_THREADS"] = str(t)
os.environ["MKL_NUM_THREADS"] = str(t)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(t)
os.environ["NUMEXPR_NUM_THREADS"] = str(t)

import numpy as np
import pandas as pd

from utils import setup_logging, maybe_path, qc_path, sort_df_chr
from combine_counts_utils import scatter_counts_to_shared_snps
from aggregation_utils import (
    adaptive_segmentation,
    assign_pos_to_range,
    count_split_genes,
    detect_phase_flips,
    gene_block_labels,
    matrix_segmentation,
)
from matplotlib.backends.backend_pdf import PdfPages
from plot_utils import plot_allele_freqs, plot_rdr_baf, plot_segmentation_qc
from switchprobs import (
    interp_cM_blocks,
    estimate_switchprobs_cM,
    estimate_switchprobs_PS,
)

SNP_KEY = ["#CHR", "POS0"]
WIN_KEY = ["#CHR", "START", "END"]

setup_logging(snakemake.log[0])

snp_info_files = list(snakemake.input["snp_info"])
tot_files = list(snakemake.input["tot_mtx_snp"])
a_files = list(snakemake.input["a_mtx_snp"])
b_files = list(snakemake.input["b_mtx_snp"])
dp_files = list(snakemake.input["dp_corrected"])
win_files = list(snakemake.input["window_df"])
sample_files = list(snakemake.input["sample_file"])

gmap_file = maybe_path(snakemake.input["gmap_file"])
region_bed = snakemake.input["region_bed"]
blacklist_bed = maybe_path(snakemake.input.get("blacklist_bed", None))
genome_size = snakemake.input["genome_size"]
gtf_file = maybe_path(snakemake.input["gtf_file"])

qc_dir = snakemake.params["qc_dir"]
qc_prefix = "combine_counts"
os.makedirs(qc_dir, exist_ok=True)
run_id = snakemake.params["run_id"]
# QC files are <qc_prefix>.<name>.bulk.<run_id>.<ext>: fold "bulk"+run_id into the stamp
qc_stamp = ".".join(p for p in ("bulk", run_id) if p)
bulk_assays = list(snakemake.params["bulk_assays"])
median_normalization = bool(snakemake.params["median_normalization"])
phase_flip_test = bool(snakemake.params["phase_flip_test"])
n_assays = len(bulk_assays)

sids_list = [pd.read_table(f) for f in sample_files]
snps_list = [pd.read_table(f, sep="\t") for f in snp_info_files]
tot_list = [np.load(f)["mat"].astype(np.int32) for f in tot_files]
a_list = [np.load(f)["mat"].astype(np.int32) for f in a_files]
b_list = [np.load(f)["mat"].astype(np.int32) for f in b_files]
dp_list = [np.load(f)["mat"] for f in dp_files]
win_list = [pd.read_table(f, sep="\t") for f in win_files]

sample_name = sids_list[0]["SAMPLE"].iloc[0]
logging.info(f"combine_counts: sample={sample_name}, bulk_assays={bulk_assays}")

has_ps = all("PS" in s.columns for s in snps_list)
has_feature = all("feature_id" in s.columns for s in snps_list)
annot_cols = (
    ["#CHR", "POS", "POS0", "region_id"]
    + (["PS"] if has_ps else [])
    + (["feature_id"] if has_feature else [])
)
snps = (
    pd.concat([s[annot_cols] for s in snps_list], ignore_index=True)
    .drop_duplicates(SNP_KEY)
)
snps = sort_df_chr(snps, ch="#CHR", pos="POS0").reset_index(drop=True)
snps["snp_row"] = np.arange(len(snps))
n_snps = len(snps)
logging.info(f"shared SNP set (union): {n_snps} SNPs across {n_assays} assays")

total_samples = sum(len(s) for s in sids_list)
tot_mtx = np.zeros((n_snps, total_samples), dtype=np.int32)
a_mtx = np.zeros((n_snps, total_samples), dtype=np.int32)
b_mtx = np.zeros((n_snps, total_samples), dtype=np.int32)
col_assay, col_repid = [], []
assay_blocks = []
offset = 0
for k in range(n_assays):
    shared_row = (
        snps_list[k][SNP_KEY]
        .merge(snps[SNP_KEY + ["snp_row"]], on=SNP_KEY, how="left")["snp_row"]
        .to_numpy()
    )
    scatter_counts_to_shared_snps(tot_mtx, tot_list[k], shared_row, offset)
    scatter_counts_to_shared_snps(a_mtx, a_list[k], shared_row, offset)
    scatter_counts_to_shared_snps(b_mtx, b_list[k], shared_row, offset)

    stypes = sids_list[k]["sample_type"].tolist()
    n_k = len(stypes)
    has_normal_k = "normal" in stypes
    tumor_cols = [offset + i for i, st in enumerate(stypes) if st == "tumor"]
    assay_blocks.append(
        {
            "assay": bulk_assays[k],
            "offset": offset,
            "n": n_k,
            "has_normal": has_normal_k,
            "normal_col": offset if has_normal_k else None,
            "tumor_cols": tumor_cols,
        }
    )
    col_assay += [bulk_assays[k]] * n_k
    col_repid += sids_list[k]["REP_ID"].tolist()
    offset += n_k

tumor_cols_all = [c for blk in assay_blocks for c in blk["tumor_cols"]]
logging.info(f"{total_samples} bulk samples, {len(tumor_cols_all)} tumor columns")

assert "region_id" in snps.columns, "invalid SNP file"
grp_cols = ["region_id"]
if not has_ps:
    logging.info("PS not in SNP columns, setting PS=1 for all SNPs")
    snps["PS"] = 1
grp_cols.append("PS")
logging.info(f"#phaseset={snps['PS'].nunique()}")

if phase_flip_test:
    snps["phase_group"] = detect_phase_flips(
        snps,
        a_mtx[:, tumor_cols_all],
        b_mtx[:, tumor_cols_all],
        grp_cols=grp_cols,
        tumor_sidx=0,
        epsilon=float(snakemake.params["phase_flip_epsilon"]),
        alpha=float(snakemake.params["phase_flip_alpha"]),
    )
    grp_cols.append("phase_group")

window_df = (
    pd.concat([w[["#CHR", "START", "END", "region_id"]] for w in win_list], ignore_index=True)
    .drop_duplicates(WIN_KEY)
)
window_df = sort_df_chr(window_df, ch="#CHR", pos="START").reset_index(drop=True)

gene_aware_binning = bool(snakemake.params["gene_aware_binning"]) and has_feature
window_df["win_idx"] = np.arange(len(window_df))
snp_window_cols = ["#CHR", "POS0", "PS"]
if phase_flip_test:
    snp_window_cols.append("phase_group")
if gene_aware_binning:
    snp_window_cols.append("feature_id")
_snps_tmp = snps[snp_window_cols].copy()
_snps_tmp = assign_pos_to_range(_snps_tmp, window_df, ref_id="win_idx", pos_col="POS0")
_snps_tmp = _snps_tmp.dropna(subset=["win_idx"])
_snps_tmp["win_idx"] = _snps_tmp["win_idx"].astype(np.int64)
snps_per_win = _snps_tmp.groupby("win_idx").size()
logging.info(
    f"SNPs per window: {len(snps_per_win)}/{len(window_df)} windows have SNPs, "
    f"mean={snps_per_win.mean():.1f}, median={snps_per_win.median():.1f}"
)
win_ps = _snps_tmp.groupby("win_idx")["PS"].agg(lambda x: x.mode().iloc[0])
window_df["PS"] = window_df["win_idx"].map(win_ps)
if window_df["PS"].isna().any():
    window_df["PS"] = window_df["PS"].ffill()

if phase_flip_test:
    win_pg = _snps_tmp.groupby("win_idx")["phase_group"].agg(lambda x: x.mode().iloc[0])
    window_df["phase_group"] = window_df["win_idx"].map(win_pg)
    if window_df["phase_group"].isna().any():
        window_df["phase_group"] = window_df["phase_group"].ffill()

if gene_aware_binning:
    # Gene blocks: every window holding a gene's SNPs -- and any window between the
    # first and last such window -- is glued into one indivisible block, so a bin can
    # span several whole genes but never a partial gene. A window shared by two genes
    # merges their blocks (inseparable at window resolution). Windows holding no gene
    # SNP are singleton blocks, keeping native window granularity.
    genic = _snps_tmp[
        _snps_tmp["feature_id"].notna() & (_snps_tmp["feature_id"] != "intergenic")
    ]
    rng = genic.groupby("feature_id")["win_idx"].agg(["min", "max"])
    window_df["gene_block"] = gene_block_labels(
        len(window_df), zip(rng["min"].to_numpy(), rng["max"].to_numpy())
    )
    logging.info(
        f"gene-aware binning: {len(rng)} genes over {len(window_df)} windows -> "
        f"{window_df['gene_block'].nunique()} gene/intergenic blocks (bins never split a gene)"
    )

max_blocksize = int(snakemake.params["max_blocksize"])
tot_tumor = np.ascontiguousarray(tot_mtx[:, tumor_cols_all], dtype=np.float64)
bbs, snps = adaptive_segmentation(
    window_df,
    snps,
    tot_tumor,
    int(snakemake.params["min_snp_reads"]),
    int(snakemake.params["min_snp_per_block"]),
    grp_cols=grp_cols,
    tumor_sidx=0,
    max_blocksize=max_blocksize,
    gene_aware=gene_aware_binning,
)
num_bbs = len(bbs)

_split = count_split_genes(snps, grp_cols)
if _split is not None:
    logging.info(
        f"gene-split sanity: {_split[0]}/{_split[1]} genes have SNPs crossing a bin "
        f"boundary (gene_aware_binning={gene_aware_binning})"
    )

bb_ids = snps["bb_id"].to_numpy()

snp_orig_idx = snps["_orig_idx"].to_numpy()
tot_mtx = tot_mtx[snp_orig_idx]
a_mtx = a_mtx[snp_orig_idx]
b_mtx = b_mtx[snp_orig_idx]

a_mtx_bb = matrix_segmentation(a_mtx, bb_ids, num_bbs)
b_mtx_bb = matrix_segmentation(b_mtx, bb_ids, num_bbs)
tot_mtx_bb = matrix_segmentation(tot_mtx, bb_ids, num_bbs)

baf_mtx_bb = np.divide(
    b_mtx_bb,
    tot_mtx_bb,
    where=tot_mtx_bb > 0,
    out=np.full_like(b_mtx_bb, np.nan, dtype=np.float32),
)

logging.info("aggregating corrected window depth into adaptive bins (per assay)")
bb_dp = np.full((num_bbs, total_samples), np.nan, dtype=np.float32)
# total aligned bases per bin = sum_w(depth_w * win_len_w); a per-segment total
# (vs bb_dp which is the length-weighted mean depth), used for the segmentation QC.
bb_bases = np.zeros((num_bbs, total_samples), dtype=np.float64)
for blk, win_a, dp_a in zip(assay_blocks, win_list, dp_list):
    w = win_a.merge(window_df[WIN_KEY + ["bin_id"]], on=WIN_KEY, how="left")
    assert w["bin_id"].notna().all(), f"{blk['assay']} windows missing from scaffold"
    bin_ids = w["bin_id"].to_numpy().astype(np.int64)
    win_lengths = (w["END"] - w["START"]).to_numpy(dtype=np.float64)
    total_len_per_bin = np.bincount(bin_ids, weights=win_lengths, minlength=num_bbs)
    for s in range(blk["n"]):
        weighted_sums = np.bincount(
            bin_ids, weights=dp_a[:, s] * win_lengths, minlength=num_bbs
        )
        bb_bases[:, blk["offset"] + s] = weighted_sums
        with np.errstate(invalid="ignore"):
            bb_dp[:, blk["offset"] + s] = weighted_sums / total_len_per_bin

logging.info(f"compute bb RDR, median_normalization={median_normalization}")
bb_rdr = np.full((num_bbs, len(tumor_cols_all)), np.nan, dtype=np.float32)
rdr_pos = {c: i for i, c in enumerate(tumor_cols_all)}
for blk, win_a, dp_a in zip(assay_blocks, win_list, dp_list):
    if blk["has_normal"] and not median_normalization:
        win_sizes = (win_a["END"] - win_a["START"]).to_numpy(dtype=np.float64)
        total_bases = np.nansum(dp_a * win_sizes[:, None], axis=0)
        library_correction = total_bases[0] / total_bases
        logging.info(f"  {blk['assay']} library factor: {library_correction}")
        normal_bb_dp = bb_dp[:, blk["normal_col"]]
        for c in blk["tumor_cols"]:
            with np.errstate(invalid="ignore", divide="ignore"):
                bb_rdr[:, rdr_pos[c]] = (
                    bb_dp[:, c] / normal_bb_dp * library_correction[c - blk["offset"]]
                )
    else:
        for c in blk["tumor_cols"]:
            col = bb_dp[:, c]
            valid_i = np.isfinite(col) & (col > 0)
            if valid_i.any():
                med = np.median(col[valid_i])
                logging.info(f"  bb median-centering {col_repid[c]}: median={med:.4f}")
                with np.errstate(invalid="ignore", divide="ignore"):
                    bb_rdr[valid_i, rdr_pos[c]] = col[valid_i] / med

rdr_outlier_quantile = float(snakemake.params["rdr_outlier_quantile"])
if rdr_outlier_quantile > 0:
    rdr_upper = np.nanquantile(bb_rdr, 1 - rdr_outlier_quantile)
    n_outlier = int(np.nansum(bb_rdr > rdr_upper))
    logging.info(
        f"RDR outlier filter: quantile={rdr_outlier_quantile}, "
        f"threshold={rdr_upper:.4f}, {n_outlier} entries set to NaN"
    )
    bb_rdr[bb_rdr > rdr_upper] = np.nan

sample_labels = [f"{col_assay[i]}:{col_repid[i]}" for i in range(total_samples)]
tumor_labels = [sample_labels[c] for c in tumor_cols_all]
rdr_ylim = (np.round(np.nanquantile(bb_rdr, 0.99)).astype(int) + 1) * 1.1

joint_sids = pd.concat(
    [sids_list[k].assign(assay_type=bulk_assays[k]) for k in range(n_assays)],
    ignore_index=True,
)
if gene_aware_binning and "feature_id" in snps.columns:
    _genic = snps[snps["feature_id"].notna() & (snps["feature_id"] != "intergenic")]
    bb_gene_count = (
        _genic.groupby("bb_id")["feature_id"].nunique()
        .reindex(range(num_bbs)).fillna(0).to_numpy()
    )
else:
    bb_gene_count = None

pdf_path = qc_path(qc_dir, qc_prefix, "combine_counts.pdf", qc_stamp)
with PdfPages(pdf_path) as pdf:
    plot_segmentation_qc(
        bbs,
        joint_sids,
        bb_bases,
        b_mtx_bb,
        tot_mtx_bb,
        pdf=pdf,
        gene_count=bb_gene_count,
    )
    plot_allele_freqs(
        bbs,
        sample_labels,
        tot_mtx_bb,
        b_mtx_bb,
        genome_size,
        qc_dir,
        apply_pseudobulk=False,
        allele="B",
        unit="bb",
        region_bed=region_bed,
        blacklist_bed=blacklist_bed,
        run_id=run_id,
        pdf=pdf,
    )
    plot_rdr_baf(
        bbs,
        bb_rdr,
        baf_mtx_bb[:, tumor_cols_all],
        tumor_labels,
        genome_size,
        pdf_path,
        unit="bb",
        rdr_ylim=rdr_ylim,
        region_bed=region_bed,
        blacklist_bed=blacklist_bed,
        pdf=pdf,
    )
logging.info(f"saved QC PDF to {pdf_path}")

nan_mask = (
    np.isnan(baf_mtx_bb).any(axis=1)
    | np.isnan(bb_dp).any(axis=1)
    | np.isnan(bb_rdr).any(axis=1)
)
n_nan_rows = int(nan_mask.sum())
n_valid = num_bbs - n_nan_rows
logging.info(
    f"NaN row filter: {n_nan_rows}/{num_bbs} bins have NaN, "
    f"keeping {n_valid} ({n_valid / max(num_bbs, 1) * 100:.1f}%)"
)

if n_nan_rows > 0:
    valid = ~nan_mask
    bbs = bbs.loc[valid].reset_index(drop=True)
    tot_mtx_bb = tot_mtx_bb[valid]
    a_mtx_bb = a_mtx_bb[valid]
    b_mtx_bb = b_mtx_bb[valid]
    baf_mtx_bb = baf_mtx_bb[valid]
    bb_dp = bb_dp[valid]
    bb_rdr = bb_rdr[valid]

kept_bb_ids = np.where(~nan_mask)[0] if n_nan_rows > 0 else np.arange(num_bbs)
old_to_new = {old: new for new, old in enumerate(kept_bb_ids)}
snps_valid = snps[snps["bb_id"].isin(old_to_new)].copy()
snps_valid["bb_id"] = snps_valid["bb_id"].map(old_to_new)
bbs["bb_id"] = np.arange(len(bbs))

logging.info("estimate bin-level switchprobs")
if gmap_file is not None:
    genetic_map = pd.read_table(gmap_file, sep="\t")
    dist_cms = interp_cM_blocks(bbs, snps_valid, genetic_map, block_id_col="bb_id")
    bbs["switchprobs"] = estimate_switchprobs_cM(
        dist_cms,
        nu=float(snakemake.params["nu"]),
        min_switchprob=float(snakemake.params["min_switchprob"]),
    )
else:
    switchprob_ps = float(snakemake.params["switchprob_ps"])
    bbs["switchprobs"] = estimate_switchprobs_PS(bbs, switchprob_ps)

bbs[["#CHR", "START", "END", "#SNPS", "region_id", "switchprobs"]].to_csv(
    snakemake.output["bb_file"], sep="\t", header=True, index=False
)
np.savez_compressed(snakemake.output["tot_mtx_bb"], mat=tot_mtx_bb)
np.savez_compressed(snakemake.output["a_mtx_bb"], mat=a_mtx_bb)
np.savez_compressed(snakemake.output["b_mtx_bb"], mat=b_mtx_bb)
np.savez_compressed(snakemake.output["baf_mtx_bb"], mat=baf_mtx_bb)
np.savez_compressed(snakemake.output["dp_mtx_bb"], mat=bb_dp)
np.savez_compressed(snakemake.output["rdr_mtx_bb"], mat=bb_rdr)

joint_sids.to_csv(snakemake.output["sample_file"], sep="\t", index=False)
logging.info("finished combine_counts.")
