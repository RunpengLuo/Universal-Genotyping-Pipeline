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
#       config["window_bed"] is read directly by the consumers for every stream
#       (get_assay_window_bed); nothing is built and the Repli-seq fetch is skipped.
#       rd_correct filters it to config["chromosomes"], so extra contigs are harmless.
#
#   else -- build the window BED (all inside the `if not use_prebuilt_windows` block):
#     . repliseq_bigwig_to_bedgraph    [do_repliseq only, per bigWig {name}]
#         in : bigWig URL (storage)            out: {name}.hg19.bedGraph
#     . repliseq_liftover              [do_repliseq + hg38 target only, per {name}]
#         in : {name}.hg19.bedGraph + chain URL (storage)
#         out: {name}.hg38.bedGraph (cached under aux); hg19 runs skip this
#     . build_window_bed              [per stream = wgs / wes]
#         in : segment_bed + reference + genome_size, and the optional
#              wes_targets (wes), mappability_bed, {reference_version} bedgraphs
#         out: {stream}_windows.bed.gz  (#CHR START END region_id seg_id GC [MAP] [REPLI])
#              + qc_pdf (segment- and window-length histograms)
#
#   . window_bed_to_3bed               [per assay_type]
#       in : {stream}_windows.bed.gz           out: pileup/{assay}/windows.bed.gz (mosdepth --by)
#
# Globals from parse_workflow: segment_bed, bedpe_files, wes_targets_files,
# has_breakpoints, use_prebuilt_windows, window_streams, do_repliseq,
# window_size_wgs, window_size_wes.
##################################################


def get_assay_window_bed(assay_type):
    """Per-assay window BED path (bulkWES -> wes stream, WGS family -> wgs stream).

    With use_prebuilt_windows, the pre-built window_bed is consumed directly for every
    stream (assumed to fit both wgs and wes): rd_correct filters it to
    config["chromosomes"] and window_bed_to_3bed only cuts columns, so extra contigs
    are harmless and nothing is built.
    """
    if use_prebuilt_windows:
        return config["window_bed"]
    stream = "wes" if assay_type == "bulkWES" else "wgs"
    return config["aux_dir"] + f"/{stream}_windows.bed.gz"


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
            config["log_dir"] + f"/build_segment_bed/{_run_id}.log",
        benchmark:
            config["bench_dir"] + f"/build_segment_bed/{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        script:
            "../scripts/build_segment_bed.py"

    # ------------------------------------------------------------------------
    # Window BED: EITHER consumed pre-built, OR built from scratch.
    #   use_prebuilt_windows (window_bed set, no BEDPE) -> config["window_bed"] is read
    #     directly for every stream (get_assay_window_bed); nothing is built here, and
    #     the Repli-seq fetch/convert is skipped (only build_window_bed consumes it).
    #   else -> build_window_bed builds the window BED for every stream (wgs and wes),
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
                    + f"/repliseq_bigwig_to_bedgraph/{{name}}.{_run_id}.log",
                benchmark:
                    config["bench_dir"]
                    + f"/repliseq_bigwig_to_bedgraph/{{name}}.{_run_id}.tsv"
                wildcard_constraints:
                    name="[A-Za-z0-9]+",
                conda:
                    "../envs/tools.yaml"
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
                        config["log_dir"] + f"/repliseq_liftover/{{name}}.{_run_id}.log",
                    benchmark:
                        config["bench_dir"]
                        + f"/repliseq_liftover/{{name}}.{_run_id}.tsv"
                    wildcard_constraints:
                        name="[A-Za-z0-9]+",
                    conda:
                        "../envs/tools.yaml"
                    resources:
                        downloads=1,
                    shell:
                        "liftOver {input.bedgraph} {input.chain} {output.bedgraph} "
                        "{output.unmapped} 2> {log}"

        rule build_window_bed:
            """Build a stream's window BED in one pass -> config["aux_dir"]/{stream}_windows.bed.gz.

            Tiles segment_bed (wgs) or wes_targets (wes), assigns region_id + seg_id, then
            annotates GC (always), MAP (when mappability_bed is set), and REPLI (when
            Repli-seq bedGraphs are available). Optional inputs are empty ([]) when absent,
            and the script skips the covariate whose input is empty.
            """
            input:
                region_bed=segment_bed,
                reference=config["reference"],
                genome_size=config["genome_size"],
                wes_targets_bed=lambda wc: (
                    [file_input(f) for f in wes_targets_files]
                    if wc.stream == "wes"
                    else []
                ),
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
                window_bed=config["aux_dir"] + "/{stream}_windows.bed.gz",
                qc_pdf=report(
                    config["qc_dir"] + "/build_window_bed.{stream}.pdf",
                    category="QC plots",
                    subcategory="window build",
                    labels={"stream": "{stream}"},
                ),
            log:
                config["log_dir"] + f"/build_window_bed/{{stream}}.{_run_id}.log",
            benchmark:
                config["bench_dir"] + f"/build_window_bed/{{stream}}.{_run_id}.tsv"
            wildcard_constraints:
                stream="(wgs|wes)",
            conda:
                "../envs/base.yaml"
            params:
                mode=lambda wc: wc.stream,
                reference_version=config["reference_version"],
                chromosomes=config["chromosomes"],
                window_size=lambda wc: (
                    window_size_wes if wc.stream == "wes" else window_size_wgs
                ),
            script:
                "../scripts/build_window_bed.py"


rule window_bed_to_3bed:
    """Per-assay headerless 3-column BED (#CHR/START/END) for mosdepth --by."""
    input:
        window_bed=lambda wc: get_assay_window_bed(wc.assay_type),
    output:
        mosdepth_bed=temp(config["pileup_dir"] + "/{assay_type}/windows.bed.gz"),
    log:
        config["log_dir"] + f"/window_bed_to_3bed/{{assay_type}}.{_run_id}.log",
    benchmark:
        config["bench_dir"] + f"/window_bed_to_3bed/{{assay_type}}.{_run_id}.tsv"
    wildcard_constraints:
        assay_type="(bulkWGS|bulkWGS-lr|bulkWES)",
    shell:
        "gzip -dc {input.window_bed} | tail -n +2 | cut -f1-3 | gzip -c "
        "> {output.mosdepth_bed} 2> {log}"
