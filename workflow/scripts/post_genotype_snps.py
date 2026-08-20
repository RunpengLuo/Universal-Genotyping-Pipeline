"""Post-process called SNPs into the run's per-chromosome het/hom-alt VCFs.

One branch per workflow mode, selected by ``params.mode``:
- ``bulk``: under ``apply_clonal_loh_hmm`` the calls came off a high-purity tumor, where
  clonal LOH collapses a gHET onto one allele, so every site is re-genotyped by a latent
  LOH-state chain over the arms. Off it, the calls are germline genotypes already and the
  outputs symlink to them.
- ``nonbulk``: sum the per-modality cellsnp-lite DP/AD/OTH per site and keep the sites
  whose allele balance reads as heterozygous (or hom-alt). There is no matched normal to
  call against, so genotype comes from the counts rather than a caller.

Last update: 2026-08-16

Inputs:
- [bulk] snp_dir/raw/chr{chrname}.vcf.gz: the bcftools calls, every chromosome
- [bulk] snp_panel: the population SNP VCF the calls were made over
- [bulk] region_bed: chromosome arms; region_id is the HMM's chain group
- [nonbulk] snp_dir/pseudobulk_{modality}/cellSNP.base.vcf.gz: cellsnp-lite base calls
- [nonbulk] genome_size: chrom sizes TSV, for the ##contig header
Outputs:
- snp_dir/chr{chrname}.vcf.gz: bi-allelic het or hom-alt SNPs
- snp_dir/chr{chrname}.vcf.gz.tbi: tabix index of the above
- qc_dir/post_genotype_snps.{bulk,nonbulk}.pdf: SNP allele-freq, one page per grouping:
  genotype, then LOH state (clonal-LOH HMM runs only)
"""

import logging

snakemake_handle = snakemake

from utils import (
    add_chr_prefix,
    maybe_list,
    set_omp_threads,
    setup_logging,
    strip_chr_prefix,
)

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from genotype_loh_hmm import GT_DTYPE, clonal_loh_hmm
from io_utils import (
    read_BED,
    read_chrom_sizes,
    read_VCF,
    symlink_files,
    write_VCF,
)
from plot_alleles import plot_allele_freqs
from range_utils import assign_pos_to_range

REFINE_PLOT_MAX_SNPS = (
    300_000  # QC scatter resolves no more than this; sampled above it
)
# blue vs red is het against hom, the call this step exists to make; grey is no-call
GT_PLOT_COLORS = {"0/0": "red", "0/1": "blue", "1/1": "red", "./.": "gray"}
# the same points regrouped onto a second QC page; grey stays the class with no call
LOH_PLOT_COLORS = {"LOH": "red", "non-LOH": "blue", "unassigned": "gray"}
LOH_COL = "LOH_STATE"
# VCF fixed fields this step does not compute, and what to write when the input lacks them
VCF_DEFAULT_FALLBACKS = {"ID": ".", "QUAL": ".", "FILTER": "PASS"}


def write_per_chrom(
    snps, chrom_col, ad, dp, oth, gt, chroms, snp_vcfs, input_nochr, genome_size
):
    """Assemble the output VCF records and write one bgzipped, indexed file per contig.

    ``ID``/``QUAL``/``FILTER`` come from *snps* wherever it carries them: this step
    re-genotypes rather than re-calls, so the caller's own fields survive it.

    Args:
        snps: Rows to write, carrying ``POS``, ``REF``, ``ALT`` and *chrom_col*.
        chrom_col: Column of *snps* holding the chr-prefixed contig.
        ad, dp, oth: ALT, total and other-allele counts, aligned to *snps*.
        gt: Genotype string per row, aligned to *snps*.
        chroms: Contigs, chr-prefixed, aligned to *snp_vcfs*.
        snp_vcfs: Output paths.
        input_nochr: Spell output contigs without the ``chr`` prefix.
        genome_size: Chrom sizes TSV, for the ``##contig`` header.
    """
    out = snps[["POS", "REF", "ALT"]].copy()
    out["#CHROM"] = snps[chrom_col].to_numpy()
    for col, placeholder in VCF_DEFAULT_FALLBACKS.items():
        out[col] = snps[col].to_numpy() if col in snps.columns else placeholder
    out["INFO"] = (
        "AD="
        + pd.Series(np.asarray(ad), index=out.index).astype(str)
        + ";DP="
        + pd.Series(np.asarray(dp), index=out.index).astype(str)
        + ";OTH="
        + pd.Series(np.asarray(oth), index=out.index).astype(str)
    )
    out["FORMAT"] = "GT"
    out["SAMPLE"] = gt

    chrom_sizes = {
        f"chr{strip_chr_prefix(k)}": v for k, v in read_chrom_sizes(genome_size).items()
    }
    by_chrom = out.groupby("#CHROM", sort=False)
    for chrom, snp_vcf in zip(chroms, snp_vcfs):
        rows = by_chrom.get_group(chrom) if chrom in by_chrom.groups else out.iloc[:0]
        write_VCF(
            rows,
            snp_vcf,
            strip_chr_prefix(chrom) if input_nochr else chrom,
            chrom_length=chrom_sizes[chrom],
        )


params = snakemake_handle.params
mode = params["mode"]
apply_clonal_loh_hmm = bool(params["apply_clonal_loh_hmm"])
chroms = list(params["chroms"])
input_nochr = params["input_nochr"]
min_dp = int(params["min_dp"])

raw_snp_vcfs = maybe_list(snakemake_handle.input["raw_snp_vcfs"])
snp_vcfs = maybe_list(snakemake_handle.output["snp_vcfs"])
genome_size = snakemake_handle.input["genome_size"]

logging.info(
    f"start post_genotype_snps, mode={mode}, "
    f"apply_clonal_loh_hmm={apply_clonal_loh_hmm}"
)

if mode == "bulk":
    # =========================================================================
    # bulk: bcftools already called these off one alignment set. Under
    # apply_clonal_loh_hmm they came off a high-purity tumor and are re-genotyped
    # by the chain; otherwise they are germline genotypes and pass through.
    # =========================================================================
    assert len(raw_snp_vcfs) == len(snp_vcfs), (
        f"bulk: {len(raw_snp_vcfs)} input VCF(s) for {len(snp_vcfs)} output(s)"
    )
    hmm_segments = snakemake_handle.output["hmm_segments"]
    hmm_params = snakemake_handle.output["hmm_params"]

    frames = []
    for raw_snp_vcf in raw_snp_vcfs:
        snps = read_VCF(
            raw_snp_vcf,
            addchr=False,
            snps_presorted=True,
            read_AD=True,
            required_cols=["GT", "AD", "DP"],
        )
        if snps is None:
            logging.warning(f"{raw_snp_vcf}: no records")
            continue
        frames.append(snps)
    assert frames, "bulk: every input VCF is empty"
    snps = pd.concat(frames, ignore_index=True)
    del frames
    snps["#CHR"] = add_chr_prefix(snps["#CHROM"].astype(str))
    if apply_clonal_loh_hmm:
        snp_panel = snakemake_handle.input["snp_panel"]
        region_bed = snakemake_handle.input["region_bed"]
        tau = params["tau"]
        logging.info(
            f"clonal-LOH HMM over {len(raw_snp_vcfs)} VCF(s), panel={snp_panel}"
        )

        rc = snps["REF_COUNT"].to_numpy()
        ac = snps["ALT_COUNT"].to_numpy()
        depth = pd.to_numeric(snps["DP"], errors="coerce").fillna(0).to_numpy()
        other = np.maximum(depth.astype(np.int64) - rc - ac, 0)

        # the chain must break wherever adjacency stops meaning anything, above all
        # across a centromere, so it is grouped by the run's own chromosome arms
        arms = read_BED(region_bed, col_id="region_id")
        qry = pd.DataFrame(
            {"#CHR": snps["#CHR"].to_numpy(), "POS0": snps["POS"].to_numpy() - 1}
        )
        qry, na_idx = assign_pos_to_range(qry, arms, ref_id="region_id", pos_col="POS0")
        group_id = pd.factorize(qry["region_id"])[0].astype(np.int64)
        logging.info(
            f"arms: {len(arms)} in {region_bed}, "
            f"{int((group_id >= 0).sum())} SNPs assigned, "
            f"{len(na_idx)} outside every arm"
        )
        fitted = group_id >= 0
        fit = clonal_loh_hmm(
            rc[fitted],
            ac[fitted],
            snps["POS"].to_numpy()[fitted],
            group_id[fitted],
            hom_laf=float(params["hom_laf"]),
            loh_laf=float(params["loh_laf"]),
            pi=float(params["pi_het"]),
            t=float(params["breakpoint_rate"]),
            # NB: YAML has no infinity, so null is the binomial
            tau=np.inf if tau is None else float(tau),
            n_retained=int(params["n_retained"]),
            n_iter=int(params["n_iter"]),
            em_tol=float(params["em_tol"]),
            margin=float(params["margin"]),
            learn_pi=bool(params["learn_pi"]),
            loh_min=float(params["loh_min"]),
            p_het_min=float(params["p_het_min"]),
            min_dp=min_dp,
        )
        gt = np.full(len(snps), "./.", dtype=GT_DTYPE)
        gt[fitted] = fit["gt"]

        # NB: state 0 is the clonal-LOH state, so this matches hmm_segments row for row
        loh_state = np.full(len(snps), "unassigned", dtype=object)
        loh_state[fitted] = np.where(fit["state"] == 0, "LOH", "non-LOH")
        snps[LOH_COL] = loh_state

        # what the caller said, so the rescue can be counted in both directions
        called = snps["GT"].astype(str).str.replace("|", "/", regex=False).to_numpy()
        logging.info(
            f"re-called het from 0/0={int(((called == '0/0') & (gt == '0/1')).sum())}, "
            f"from 1/1={int(((called == '1/1') & (gt == '0/1')).sum())}"
        )

        # the Viterbi path as intervals: the model's own allelic-imbalance segmentation.
        # A new segment starts wherever the state changes or the arm does.
        arm = group_id[fitted]
        state = fit["state"]
        cut = np.r_[True, (state[1:] != state[:-1]) | (arm[1:] != arm[:-1])]
        seg = pd.DataFrame(
            {
                "#CHR": snps.loc[fitted, "#CHR"].to_numpy(),
                "POS": snps.loc[fitted, "POS"].to_numpy(),
                "state": state,
                "maf": fit["maf"],
                "is_het": fit["gt"] == "0/1",
                "seg": np.cumsum(cut),
            }
        )
        segments = seg.groupby("seg", sort=True).agg(
            START=("POS", "min"),
            END=("POS", "max"),
            state=("state", "first"),
            n_snps=("POS", "size"),
            n_het=("is_het", "sum"),
            mean_maf=("maf", "mean"),
        )
        segments.insert(0, "#CHR", seg.groupby("seg", sort=True)["#CHR"].first())
        segments["START"] -= 1
        segments.insert(4, "b_s", fit["means"][segments["state"].to_numpy()])
        segments.to_csv(hmm_segments, sep="\t", index=False, float_format="%.4f")
        pd.DataFrame(
            [
                {
                    "means": ",".join(f"{b:.4f}" for b in fit["means"]),
                    "pi": float(np.mean(fit["pi"])),
                    "loglik": fit["loglik"],
                    "n_iter": fit["n_iter"],
                    "converged": fit["converged"],
                    "n_sites": int(fitted.sum()),
                    "n_arms": len(arms),
                    "n_segments": len(segments),
                }
            ]
        ).to_csv(hmm_params, sep="\t", index=False, float_format="%.6g")
        logging.info(f"{len(segments)} HMM segments -> {hmm_segments}")

        # from here GT is this step's call, not the caller's
        snps["GT"] = gt
        keep = (gt == "0/1") | (gt == "1/1")
        logging.info(f"#kept SNPs={int(keep.sum())}/{len(snps)}")
        write_per_chrom(
            snps.loc[keep],
            "#CHR",
            ac[keep],
            (rc + ac)[keep],
            other[keep],
            snps.loc[keep, "GT"].to_numpy(),
            chroms,
            snp_vcfs,
            input_nochr,
            genome_size,
        )
        genotyped = snps
    else:
        # no fit ran, so both aux tables are headers only
        pd.DataFrame(
            columns=[
                "#CHR",
                "START",
                "END",
                "state",
                "b_s",
                "n_snps",
                "n_het",
                "mean_maf",
            ]
        ).to_csv(hmm_segments, sep="\t", index=False)
        pd.DataFrame(
            columns=[
                "means",
                "pi",
                "loglik",
                "n_iter",
                "converged",
                "n_sites",
                "n_arms",
                "n_segments",
            ]
        ).to_csv(hmm_params, sep="\t", index=False)
        logging.info(f"linking {len(raw_snp_vcfs)} VCF(s) and indexes")
        symlink_files(raw_snp_vcfs, snp_vcfs)
        symlink_files(
            maybe_list(snakemake_handle.input["raw_snp_vcfs_tbi"]),
            maybe_list(snakemake_handle.output["snp_vcfs_tbi"]),
        )
        # the calls are kept as they are, so their own GT is what the plot shows
        genotyped = snps
else:
    # =========================================================================
    # single-cell: no matched normal to call against, so genotype comes from the
    # per-modality cellsnp-lite counts summed per site.
    # =========================================================================
    filter_nz_OTH = params["filter_nz_OTH"]
    filter_hom_ALT = params["filter_hom_ALT"]
    min_het_reads = int(params["min_het_reads"])
    min_vaf_thres = float(params["min_vaf_thres"])
    modalities = list(params["modalities"])
    logging.info(
        f"filter_nz_OTH={filter_nz_OTH}, filter_hom_ALT={filter_hom_ALT}, "
        f"min_het_reads={min_het_reads}, min_dp={min_dp}, "
        f"min_vaf_thres={min_vaf_thres}"
    )

    KEY = ["#CHROM", "POS", "REF", "ALT"]
    CNT = ["DP", "AD", "OTH"]
    CARRY = list(VCF_DEFAULT_FALLBACKS)
    raw_snps_list = []
    for idx, modality in enumerate(modalities):
        raw_snps = read_VCF(raw_snp_vcfs[idx], addkey=True, required_cols=KEY + CNT)
        for cnt in CNT:
            raw_snps[f"{cnt}{idx}"] = raw_snps[cnt].astype(np.int64)

        raw_snps = raw_snps[raw_snps["#CHROM"].astype(str).isin(chroms)]

        n_dup_rows = raw_snps.duplicated(subset="KEY", keep=False).sum()
        if n_dup_rows > 0:
            n_dup_keys = raw_snps.loc[
                raw_snps.duplicated(subset="KEY", keep=False), "KEY"
            ].nunique()
            logging.warning(
                f"{modality} have {n_dup_rows} rows with duplicated positions={n_dup_keys}"
            )
            logging.warning("drop duplicated rows.")
            raw_snps = raw_snps.drop_duplicates(subset="KEY", keep="first").reset_index(
                drop=True
            )
        raw_snps_list.append(raw_snps)

    base_snps = pd.concat(
        [
            raw_snps[
                ["KEY", "#CHROM", "POS", "REF", "ALT"]
                + [c for c in CARRY if c in raw_snps.columns]
            ]
            for raw_snps in raw_snps_list
        ],
        ignore_index=True,
    ).drop_duplicates(subset=["#CHROM", "POS", "REF", "ALT"], keep="first")

    dup_mask = base_snps.duplicated(subset="KEY", keep=False)
    n_dup_rows = int(dup_mask.sum())
    if n_dup_rows > 0:
        n_dup_keys = int(base_snps.loc[dup_mask, "KEY"].drop_duplicates().shape[0])
        logging.warning(
            f"merged df have {n_dup_rows} rows with duplicated positions={n_dup_keys}"
        )
        logging.warning("drop duplicated rows.")
        base_snps = base_snps.loc[~dup_mask, :]
    base_snps = base_snps.reset_index(drop=True)
    nsnps_before_genotyping = len(base_snps)

    base_snps = base_snps.sort_values(["#CHROM", "POS"], kind="mergesort")
    for cnt in CNT:
        base_snps[cnt] = 0
    for idx, raw_snps in enumerate(raw_snps_list):
        right_cols = ["KEY"] + [f"{cnt}{idx}" for cnt in CNT]
        base_snps = pd.merge(
            left=base_snps,
            right=raw_snps[right_cols],
            on="KEY",
            how="left",
            validate="one_to_one",
            sort=False,
        )
        for cnt in CNT:
            base_snps[cnt] += base_snps[f"{cnt}{idx}"].fillna(0).astype(np.int64)
        base_snps = base_snps.drop(columns=[f"{cnt}{idx}" for cnt in CNT])

    ref_count = base_snps["DP"] - base_snps["AD"]
    allele_ratio = base_snps["AD"] / base_snps["DP"]
    called = base_snps["DP"] >= min_dp
    is_het = (
        called
        & (np.minimum(base_snps["AD"], ref_count) >= min_het_reads)
        & allele_ratio.between(min_vaf_thres, 1 - min_vaf_thres)
    )
    # NB: hom-ref never survives the keep below, so it is not genotyped
    is_hom_alt = called & (base_snps["AD"] == base_snps["DP"])

    # the same four classes the bulk branch reports, so one QC plot serves both
    gt = np.full(len(base_snps), "./.", dtype=GT_DTYPE)
    gt[called & (base_snps["AD"] == 0)] = "0/0"
    gt[is_hom_alt] = "1/1"
    gt[is_het] = "0/1"

    keep = is_het if filter_hom_ALT else (is_het | is_hom_alt)
    if filter_nz_OTH:
        logging.info("SNPs with nonzero OTHs are filtered.")
        keep &= base_snps["OTH"] == 0
    logging.info(
        f"#nz-OTH SNPs={int((base_snps['OTH'] > 0).sum())}/{nsnps_before_genotyping}"
    )

    final_snps = base_snps[keep]
    logging.info(f"#kept SNPs={len(final_snps)}/{nsnps_before_genotyping}")
    write_per_chrom(
        final_snps,
        "#CHROM",
        final_snps["AD"].to_numpy(),
        final_snps["DP"].to_numpy(),
        final_snps["OTH"].to_numpy(),
        gt[keep],
        chroms,
        snp_vcfs,
        input_nochr,
        genome_size,
    )

    genotyped = base_snps.assign(
        REF_COUNT=base_snps["DP"] - base_snps["AD"],
        ALT_COUNT=base_snps["AD"],
        GT=gt,
        **{"#CHR": base_snps["#CHROM"].astype(str)},
    )

if len(genotyped) > REFINE_PLOT_MAX_SNPS:
    sampled = np.random.default_rng(0).choice(
        len(genotyped), REFINE_PLOT_MAX_SNPS, replace=False
    )
    logging.info(
        f"QC plot: sampled {REFINE_PLOT_MAX_SNPS} of {len(genotyped)} sites, "
        f"the scatter would not resolve more"
    )
    genotyped = genotyped.iloc[np.sort(sampled)].reset_index(drop=True)

ref_counts = genotyped["REF_COUNT"].to_numpy()
depths = ref_counts + genotyped["ALT_COUNT"].to_numpy()

# NB: the pass-through branch carries the caller's GT, which may be phased
genotyped["GT_PLOT"] = genotyped["GT"].astype(str).str.replace("|", "/", regex=False)
# one page per grouping of the same points; the LOH page needs the chain to have run
qc_pages = [("GT_PLOT", GT_PLOT_COLORS, "")]
if LOH_COL in genotyped.columns:
    qc_pages.append((LOH_COL, LOH_PLOT_COLORS, " - LOH state"))
logging.info(f"QC pages: {', '.join(col for col, _, _ in qc_pages)}")

with PdfPages(snakemake_handle.output["qc_pdf"]) as qc_pdf:
    for group_col, group_colors, title_suffix in qc_pages:
        plot_allele_freqs(
            genotyped,
            [],
            [],
            [],
            depths.reshape(-1, 1),
            ref_counts.reshape(-1, 1),
            genome_size,
            params["qc_dir"],
            apply_pseudobulk=True,
            cell_dataset_ids=None,
            allele="ref",
            feature_label="SNP",
            snp_groups=genotyped[group_col].to_numpy(),
            group_colors=group_colors,
            run_id=params["run_id"],
            sample_id=params["sample_id"],
            name_prefix="post_genotype_snps",
            pdf=qc_pdf,
            title_suffix=title_suffix,
        )
logging.info("finished post_genotype_snps")
