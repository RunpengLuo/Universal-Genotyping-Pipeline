##################################################
# SNP-informed adaptive binning + window→bin depth aggregation
# Bulk: combine_counts (with corrected window depth)
# Non-bulk: combine_counts_nonbulk + cnv_segmentation
##################################################


def _bulk_inputs(subdir, fname):
    return [f"{subdir}/{at}/{fname}" for at in BULK_PRESENT]


def _nonbulk_inputs(subdir, fname):
    return [f"{subdir}/{at}/{fname}" for at in NONBULK_PRESENT]


def _nonbulk_outputs(prefix, suffix):
    return [config["bb_dir"] + f"/{prefix}.{at}.{suffix}" for at in NONBULK_PRESENT]


rule combine_counts:
    input:
        dp_corrected=_bulk_inputs(config["pileup_dir"], "window.dp.npz"),
        window_df=_bulk_inputs(config["pileup_dir"], "window.tsv.gz"),
        snp_info=_bulk_inputs(config["allele_dir"], "snps.tsv.gz"),
        tot_mtx_snp=_bulk_inputs(config["allele_dir"], "snp.Tallele.npz"),
        a_mtx_snp=_bulk_inputs(config["allele_dir"], "snp.Aallele.npz"),
        b_mtx_snp=_bulk_inputs(config["allele_dir"], "snp.Ballele.npz"),
        sample_file=_bulk_inputs(config["allele_dir"], "sample_ids.tsv"),
        gmap_file=lambda wc: (
            config["phase_dir"] + "/genetic_map.tsv.gz" if require_genetic_map else []
        ),
        region_bed=config["region_bed"],
        blacklist_bed=config["blacklist_bed"] or [],
        genome_size=config["genome_size"],
        gtf_file=config["gtf_file"],
    output:
        bb_file=config["bb_dir"] + "/bb.tsv.gz",
        tot_mtx_bb=config["bb_dir"] + "/bb.Tallele.npz",
        a_mtx_bb=config["bb_dir"] + "/bb.Aallele.npz",
        b_mtx_bb=config["bb_dir"] + "/bb.Ballele.npz",
        baf_mtx_bb=config["bb_dir"] + "/bb.baf.npz",
        dp_mtx_bb=config["bb_dir"] + "/bb.depth.npz",
        rdr_mtx_bb=config["bb_dir"] + "/bb.rdr.npz",
        sample_file=config["bb_dir"] + "/sample_ids.tsv",
    params:
        qc_dir=config["qc_dir"],
        bulk_assays=BULK_PRESENT,
        nu=config["params_combine_counts"]["nu"],
        min_switchprob=config["params_combine_counts"]["min_switchprob"],
        switchprob_ps=config["params_combine_counts"]["switchprob_ps"],
        min_snp_reads=config["params_combine_counts"]["min_snp_reads"],
        min_snp_per_block=config["params_combine_counts"]["min_snp_per_block"],
        gene_aware_binning=config["params_combine_counts"]["gene_aware_binning"],
        median_normalization=config["params_combine_counts"]["median_normalization"],
        rdr_outlier_quantile=config["params_combine_counts"]["rdr_outlier_quantile"],
        max_blocksize=config["params_combine_counts"]["max_blocksize"],
        phase_flip_test=config["params_combine_counts"]["phase_flip_test"],
        phase_flip_epsilon=config["params_combine_counts"]["phase_flip_epsilon"],
        phase_flip_alpha=config["params_combine_counts"]["phase_flip_alpha"],
        run_id=_run_id,
    threads: 1
    log:
        config["log_dir"] + f"/combine_counts.{_run_id}.log",
    conda:
        "../envs/base.yaml"
    script:
        """../scripts/combine_counts.py"""


rule combine_counts_nonbulk:
    input:
        snp_info=_nonbulk_inputs(config["allele_dir"], "snps.tsv.gz"),
        tot_mtx_snp=_nonbulk_inputs(config["allele_dir"], "snp.Tallele.npz"),
        a_mtx_snp=_nonbulk_inputs(config["allele_dir"], "snp.Aallele.npz"),
        b_mtx_snp=_nonbulk_inputs(config["allele_dir"], "snp.Ballele.npz"),
        sample_file=_nonbulk_inputs(config["allele_dir"], "sample_ids.tsv"),
        all_barcodes=_nonbulk_inputs(config["allele_dir"], "barcodes.tsv.gz"),
        barcodes_full=_nonbulk_inputs(config["allele_dir"], "barcodes.full.tsv.gz"),
        gmap_file=(
            config["phase_dir"] + "/genetic_map.tsv.gz" if require_genetic_map else []
        ),
        region_bed=config["region_bed"],
        genome_size=config["genome_size"],
        gtf_file=config["gtf_file"],
    output:
        bb_file=config["bb_dir"] + "/bb.tsv.gz",
        sample_file=config["bb_dir"] + "/sample_ids.tsv",
        tot_mtx_bb=_nonbulk_outputs("bb", "Tallele.npz"),
        a_mtx_bb=_nonbulk_outputs("bb", "Aallele.npz"),
        b_mtx_bb=_nonbulk_outputs("bb", "Ballele.npz"),
        baf_mtx_bb=_nonbulk_outputs("bb", "baf.npz"),
        multi_snp_file=_nonbulk_outputs("multi_snp", "tsv.gz"),
        tot_mtx_multi=_nonbulk_outputs("multi_snp", "Tallele.npz"),
        a_mtx_multi=_nonbulk_outputs("multi_snp", "Aallele.npz"),
        b_mtx_multi=_nonbulk_outputs("multi_snp", "Ballele.npz"),
        all_barcodes=_nonbulk_outputs("barcodes", "tsv.gz"),
        barcodes_full=_nonbulk_outputs("barcodes", "full.tsv.gz"),
    params:
        qc_dir=config["qc_dir"],
        nonbulk_assays=NONBULK_PRESENT,
        nu=config["params_combine_counts"]["nu"],
        min_switchprob=config["params_combine_counts"]["min_switchprob"],
        switchprob_ps=config["params_combine_counts"]["switchprob_ps"],
        nsnp_multi=config["params_combine_counts"]["nsnp_multi"],
        min_snp_reads=config["params_combine_counts"]["min_snp_reads"],
        min_snp_per_block=config["params_combine_counts"]["min_snp_per_block"],
        gene_aware_binning=config["params_combine_counts"]["gene_aware_binning"],
        run_id=_run_id,
    threads: 1
    log:
        config["log_dir"] + f"/combine_counts_nonbulk.{_run_id}.log",
    conda:
        "../envs/base.yaml"
    script:
        """../scripts/combine_counts_nonbulk.py"""


rule cnv_segmentation:
    input:
        snp_info=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snps.tsv.gz",
        tot_mtx_snp=lambda wc: config["allele_dir"]
        + f"/{wc.assay_type}/snp.Tallele.npz",
        a_mtx_snp=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snp.Aallele.npz",
        b_mtx_snp=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snp.Ballele.npz",
        sample_file=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/sample_ids.tsv",
        all_barcodes=config["allele_dir"] + "/{assay_type}/barcodes.tsv.gz",
        barcodes_full=config["allele_dir"] + "/{assay_type}/barcodes.full.tsv.gz",
        h5ad_file=config["bb_dir"] + "/{assay_type}/{assay_type}.h5ad",
        region_bed=lambda wc: config["region_bed"],
        genome_size=lambda wc: config["genome_size"],
        gtf_file=lambda wc: config["gtf_file"],
        bb_file=lambda wc: config["bb_file"] or [],
    output:
        cnv_segments=config["bb_dir"] + "/{assay_type}/cnv_segments.tsv",
        x_count=config["bb_dir"] + "/{assay_type}/bb.Xcount.npz",
        tot_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.Tallele.npz",
        a_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.Aallele.npz",
        b_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.Ballele.npz",
        barcodes_out=config["bb_dir"] + "/{assay_type}/barcodes.tsv.gz",
        barcodes_full_out=config["bb_dir"] + "/{assay_type}/barcodes.full.tsv.gz",
        sample_file=config["bb_dir"] + "/{assay_type}/sample_ids.tsv",
    wildcard_constraints:
        assay_type="(scRNA|scATAC|VISIUM|VISIUM3prime)",
    params:
        qc_dir=config["qc_dir"],
        sample_name=SAMPLE_ID,
        assay_type=lambda wc: wc.assay_type,
        feature_type=lambda wc: ASSAY_TYPE2FEATURE_TYPE[wc.assay_type],
        run_id=_run_id,
    threads: 1
    log:
        config["log_dir"] + f"/cnv_segmentation.{{assay_type}}.{_run_id}.log",
    conda:
        "../envs/base.yaml"
    script:
        """../scripts/cnv_segmentation.py"""
