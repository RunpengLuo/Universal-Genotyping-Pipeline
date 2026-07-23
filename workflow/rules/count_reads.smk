##################################################
# Read depth computation and bias correction (bulk only)
#
# bulkWGS/bulkWES: mosdepth fixed-window depth + HMMcopy-style LOWESS correction
#
# Produces window.dp.npz + window.tsv.gz consumed by combine_counts.
##################################################

_rdr_cfg = config["params_count_reads"]


rule run_mosdepth:
    input:
        alignment=lambda wc: alignment_input(get_data[(wc.assay_type, wc.dataset_id)]),
        alignment_index=lambda wc: alignment_index_input(
            get_data[(wc.assay_type, wc.dataset_id)]
        ),
        windows_bed=config["pileup_dir"] + "/{assay_type}/windows.bed.gz",
    output:
        mosdepth_file=config["pileup_dir"]
        + "/{assay_type}/out_mosdepth/{dataset_id}.regions.bed.gz",
    log:
        config["log_dir"]
        + f"/run_mosdepth/run_mosdepth.{{assay_type}}_{{dataset_id}}.{_run_id}.log",
    wildcard_constraints:
        assay_type="(bulkWGS|bulkWGS-lr|bulkWES)",
    conda:
        "../envs/tools.yaml"
    threads: config["threads"]["mosdepth"]
    resources:
        downloads=lambda wc: download_slots(get_data[(wc.assay_type, wc.dataset_id)]),
    params:
        out_prefix=config["pileup_dir"] + "/{assay_type}/out_mosdepth/{dataset_id}",
        read_quality=config["params_mosdepth"]["read_quality"],
        extra_params=config["params_mosdepth"]["extra_params"],
    shell:
        r"""
        mosdepth \
            -t {threads} \
            -Q {params.read_quality} \
            --by {input.windows_bed} \
            {params.extra_params} \
            {params.out_prefix} {input.alignment} > {log} 2>&1
        """


rule rd_correct:
    """Per-window LOWESS bias correction (GC/mappability/replication timing)."""
    input:
        mosdepth_files=lambda wc: [
            config["pileup_dir"]
            + f"/{wc.assay_type}/out_mosdepth/{dataset_id}.regions.bed.gz"
            for dataset_id in assay2dataset_ids[wc.assay_type]
        ],
        window_bed=lambda wc: get_assay_window_bed(wc.assay_type),
        genome_size=config["genome_size"],
        region_bed=segment_bed,
        blacklist_bed=config["blacklist_bed"] or [],
    output:
        dp_corrected=config["pileup_dir"] + "/{assay_type}/window.dp.npz",
        window_df=config["pileup_dir"] + "/{assay_type}/window.tsv.gz",
        depth_stats=report(
            config["pileup_dir"] + "/{assay_type}/depth_statistics.tsv",
            category="QC stats",
            subcategory="depth",
            labels={"table": "depth statistics", "assay": "{assay_type}"},
        ),
        qc_pdf=report(
            config["qc_dir"] + "/rd_correction.{assay_type}.pdf",
            category="QC plots",
            subcategory="read-depth correction",
            labels={"assay": "{assay_type}"},
        ),
    log:
        config["log_dir"] + f"/rd_correct.{{assay_type}}.{_run_id}.log",
    wildcard_constraints:
        assay_type="(bulkWGS|bulkWGS-lr|bulkWES)",
    conda:
        "../envs/base.yaml"
    params:
        qc_dir=config["qc_dir"],
        sample_name=sample_id,
        dataset_ids=lambda wc: assay2dataset_ids[wc.assay_type],
        sample_ids=lambda wc: [
            f"{sample_id}_{dataset_id}"
            for dataset_id in assay2dataset_ids[wc.assay_type]
        ],
        mosdepth_dir=lambda wc: config["pileup_dir"] + f"/{wc.assay_type}/out_mosdepth",
        chromosomes=config["chromosomes"],
        samplesize=_rdr_cfg["samplesize"],
        routlier=_rdr_cfg["routlier"],
        doutlier=_rdr_cfg["doutlier"],
        min_mappability=_rdr_cfg["min_mappability"],
        gc_correct=_rdr_cfg["gc_correct"],
        gc_correct_method=_rdr_cfg["gc_correct_method"],
        rt_correct=_rdr_cfg["rt_correct"],
        assay_type=lambda wc: wc.assay_type,
        run_id=_run_id,
    script:
        """../scripts/rd_correct.py"""
