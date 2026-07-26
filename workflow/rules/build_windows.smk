##################################################
# Bulk window construction (gated on workflow_mode == "bulk_genotyping").
# Rule flow (region.bed stays the raw arm-level file; downloads/tool calls are
# Snakemake rules, only binning is Python):
#
#   . build_segment_bed                [always for bulk]
#       in : region_bed (raw arm BED) + blacklist_bed + bedpe
#       out: segment_bed  (#CHR START END region_id[arm] seg_id[chunk]); also read by
#            phase_and_concat (SNP region/seg assignment), rd_correct + combine_counts
#            (QC region overlay) -- so it runs even when windows are pre-built.
#
# The window BED is EITHER consumed pre-built OR built (see the if/else below):
#
#   use_prebuilt_windows [window_bed set, no BEDPE]:
#       config["window_bed"] is read directly by the consumers (window_bed_path);
#       nothing is built and the Repli-seq fetch is skipped. rd_correct filters it to
#       config["chromosomes"], so extra contigs are harmless.
#
#   else -- build the window BED (all inside the `if not use_prebuilt_windows` block):
#     . repliseq_bigwig_to_bedgraph    [do_repliseq only, per bigWig {name}]
#         in : bigWig URL (storage)            out: {name}.hg19.bedGraph
#     . repliseq_liftover              [do_repliseq + hg38 target only, per {name}]
#         in : {name}.hg19.bedGraph + chain URL (storage)
#         out: {name}.hg38.bedGraph (cached under aux); hg19 runs skip this
#     . build_window_bed              [one grid for every bulk assay: WGS/WGS-lr/WES]
#         in : segment_bed + reference + genome_size, and the optional
#              mappability_bed, {reference_version} bedgraphs
#         out: windows.bed.gz  (#CHR START END region_id seg_id GC [MAP] [REPLI])
#              + qc_pdf (segment- and window-length histograms)
#
#   . window_bed_to_3bed               [per assay_type]
#       in : windows.bed.gz           out: pileup/{assay}/windows.bed.gz (mosdepth --by)
#
# Globals from parse_workflow: segment_bed, bedpe_files, has_breakpoints,
# use_prebuilt_windows, do_repliseq, window_size.
##################################################

# One window BED for every bulk assay (WGS/WGS-lr/WES share the segment.bed grid).
window_bed_path = (
    config["window_bed"]
    if use_prebuilt_windows
    else config["aux_dir"] + "/windows.bed.gz"
)


if workflow_mode == "bulk_genotyping":

    rule build_segment_bed:
        """Segment BED: region_id (arm) + seg_id (breakpoint chunk).

        Subtracts the blacklist from the raw region.bed and splits each arm at the
        union of all datasets' BEDPE breakpoints. Always runs for bulk; with no BEDPE
        every arm is one segment (seg_id one-per-arm).
        """
        input:
            region_bed=config["region_bed"],
            blacklist_bed=config["blacklist_bed"] or [],
            bedpe=[file_input(f) for f in bedpe_files],
        output:
            segment_bed=segment_bed,
        log:
            config["log_dir"] + f"/build_segment_bed/build_segment_bed.{_run_id}.log",
        benchmark:
            config["bench_dir"] + f"/build_segment_bed/build_segment_bed.{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        script:
            "../scripts/build_segment_bed.py"

    # ------------------------------------------------------------------------
    # Window BED: EITHER consumed pre-built, OR built from scratch (window_bed_path).
    #   use_prebuilt_windows (window_bed set, no BEDPE) -> config["window_bed"] is read
    #     directly; nothing is built here, and the Repli-seq fetch/convert is skipped
    #     (only build_window_bed consumes it).
    #   else -> build_window_bed builds the one shared window BED (all bulk assays),
    #     fed by the Repli-seq bedGraphs staged just below (when do_repliseq).
    # ------------------------------------------------------------------------
    if not use_prebuilt_windows:
        if do_repliseq:
            _repli_cache = config["aux_dir"] + "/repliseq"
            _repli_names = [f[: -len(".bigWig")] for f in REPLISEQ_BIGWIG_FILES]
            # Repli-seq bigWigs are hg19; lift to hg38 only when the run is hg38.
            _repli_target = config["reference_version"]
            _repli_lift = _repli_target != "hg19"

            rule repliseq_bigwig_to_bedgraph:
                """Fetch an ENCODE Repli-seq bigWig (hg19) and convert to bedGraph."""
                input:
                    bigwig=lambda wc: file_input(
                        f"{UCSC_REPLISEQ_BASE}/{wc.name}.bigWig"
                    ),
                output:
                    (
                        temp(_repli_cache + "/{name}.hg19.bedGraph")
                        if _repli_lift
                        else _repli_cache + "/{name}.hg19.bedGraph"
                    ),
                log:
                    config["log_dir"]
                    + f"/repliseq_bigwig_to_bedgraph/repliseq_bigwig_to_bedgraph.{{name}}.{_run_id}.log",
                benchmark:
                    config["bench_dir"]
                    + f"/repliseq_bigwig_to_bedgraph/repliseq_bigwig_to_bedgraph.{{name}}.{_run_id}.tsv"
                wildcard_constraints:
                    name="[A-Za-z0-9]+",
                conda:
                    "../envs/ucsc.yaml"
                resources:
                    downloads=1,
                shell:
                    "bigWigToBedGraph {input.bigwig} {output} 2> {log}"

            if _repli_lift:

                rule repliseq_liftover:
                    """liftOver an hg19 Repli-seq bedGraph to hg38 (cached under aux)."""
                    input:
                        bedgraph=_repli_cache + "/{name}.hg19.bedGraph",
                        chain=file_input(LIFTOVER_CHAIN_URL),
                    output:
                        bedgraph=_repli_cache + "/{name}.hg38.bedGraph",
                        unmapped=temp(_repli_cache + "/{name}.unmapped"),
                    log:
                        config["log_dir"]
                        + f"/repliseq_liftover/repliseq_liftover.{{name}}.{_run_id}.log",
                    benchmark:
                        config["bench_dir"]
                        + f"/repliseq_liftover/repliseq_liftover.{{name}}.{_run_id}.tsv"
                    wildcard_constraints:
                        name="[A-Za-z0-9]+",
                    conda:
                        "../envs/ucsc.yaml"
                    resources:
                        downloads=1,
                    shell:
                        "liftOver {input.bedgraph} {input.chain} {output.bedgraph} "
                        "{output.unmapped} 2> {log}"

        rule build_window_bed:
            """Build the window BED in one pass -> config["aux_dir"]/windows.bed.gz.

            Tiles segment_bed, assigns region_id + seg_id, then annotates GC (always),
            MAP (when mappability_bed is set), and REPLI (when Repli-seq bedGraphs are
            available). Optional inputs are empty ([]) when absent, and the script skips
            the covariate whose input is empty. One grid for every bulk assay.
            """
            input:
                region_bed=segment_bed,
                reference=config["reference"],
                genome_size=config["genome_size"],
                mappability_bed=config.get("mappability_bed") or [],
                bedgraphs=(
                    [
                        _repli_cache + f"/{n}.{_repli_target}.bedGraph"
                        for n in _repli_names
                    ]
                    if do_repliseq
                    else []
                ),
            output:
                window_bed=config["aux_dir"] + "/windows.bed.gz",
                qc_pdf=report(
                    config["qc_dir"] + "/build_window_bed.pdf",
                    category="QC plots",
                    subcategory="window build",
                ),
            log:
                config["log_dir"] + f"/build_window_bed/build_window_bed.{_run_id}.log",
            benchmark:
                config["bench_dir"]
                + f"/build_window_bed/build_window_bed.{_run_id}.tsv"
            conda:
                "../envs/base.yaml"
            params:
                reference_version=config["reference_version"],
                chromosomes=config["chromosomes"],
                window_size=window_size,
            script:
                "../scripts/build_window_bed.py"


rule window_bed_to_3bed:
    """Per-assay headerless 3-column BED (#CHR/START/END) for mosdepth --by."""
    input:
        window_bed=window_bed_path,
    output:
        mosdepth_bed=temp(config["pileup_dir"] + "/{assay_type}/windows.bed.gz"),
    log:
        config["log_dir"]
        + f"/window_bed_to_3bed/window_bed_to_3bed.{{assay_type}}.{_run_id}.log",
    benchmark:
        config["bench_dir"]
        + f"/window_bed_to_3bed/window_bed_to_3bed.{{assay_type}}.{_run_id}.tsv"
    wildcard_constraints:
        assay_type="(bulkWGS|bulkWGS-lr|bulkWES)",
    shell:
        "gzip -dc {input.window_bed} | tail -n +2 | cut -f1-3 | gzip -c "
        "> {output.mosdepth_bed} 2> {log}"
