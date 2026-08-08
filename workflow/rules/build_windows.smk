##################################################
# Segment BED + window BED, both built in EVERY mode.
# Bulk counts the windows (mosdepth -> rd_correct); single-cell uses them as the
# fixed-bin skeleton for binning. The build is identical in both: GC and MAP are
# annotated whenever their config input is set. Only the Repli-seq FETCH is bulk-only,
# because do_repliseq gates it.
# Rule flow (downloads/tool calls are Snakemake rules, only binning is Python):
#
#   . build_segment_bed                [always, every mode]
#       in : segment_bed (configured segmentation) + region_bed (arms) + blacklist_bed
#       out: segment_bed  (#CHR START END region_id[arm] seg_id[segment]); also read by
#            phase_and_concat (SNP region/seg assignment), rd_correct + combine_counts
#            (QC region overlay) -- so it runs even when windows are pre-built.
#
# The window BED is EITHER built from the segments OR consumed pre-built as-is:
#
#   build_windows [no window_bed]:
#     . repliseq_bigwig_to_bedgraph    [do_repliseq (bulk only), per bigWig {name}]
#         in : bigWig URL (storage)            out: {name}.hg19.bedGraph
#     . repliseq_liftover              [do_repliseq + hg38 target only, per {name}]
#         in : {name}.hg19.bedGraph + chain URL (storage)
#         out: {name}.hg38.bedGraph (cached under aux); hg19 runs skip this
#     . build_window_bed              [one bin set for every assay of the run]
#         in : segment_bed + genome_size + reference + mappability_bed
#              + {reference_version} bedgraphs
#         out: windows.bed.gz  (#CHR START END region_id seg_id [GC] [MAP] [REPLI])
#              segment- and window-length histograms are logged, not plotted
#         Tiling is per segment row, so no window spans two segments.
#
#   else [window_bed set]: no rule at all. The configured file is read directly by the
#   consumers; its region_id/seg_id are assumed to be the segment_bed ids, and
#   rd_correct filters it to config["chromosomes"], so extra contigs are harmless.
#
#   . window_bed_to_3bed               [one shared file for every bulk assay]
#       in : windows.bed.gz           out: aux/windows.3col.bed.gz (mosdepth --by)
#
# Globals from parse_workflow: segment_bed, window_bed (built or configured),
# build_windows, do_repliseq, window_size.
##################################################

rule build_segment_bed:
    """Segment BED: region_id (arm) + seg_id (the configured segment).

    Stamps each configured segment with the arm it sits in and subtracts the
    blacklist. With segment_bed == region_bed there is one segment per arm.
    """
    input:
        segments=config["segment_bed"],
        region_bed=config["region_bed"],
        blacklist_bed=config["blacklist_bed"] or [],
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
# Window BED: EITHER built from the segments, OR supplied and used as-is.
# The Repli-seq fetch/convert only feeds build_window_bed.
# ------------------------------------------------------------------------
if build_windows:
    if do_repliseq:
        _repli_cache = config["aux_dir"] + "/repliseq"
        _repli_names = [f[: -len(".bigWig")] for f in REPLISEQ_BIGWIG_FILES]
        # Repli-seq bigWigs are hg19; lift to hg38 only when the run is hg38.
        _repli_target = reference_version
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
        """Build the window BED in one pass -> aux_dir/windows.bed.gz.

        Tiles segment_bed and assigns region_id + seg_id, then annotates the
        bias-correction covariates GC / MAP / REPLI. Each is gated on its own config
        input being non-empty, the same way in every mode. One bin set for every assay
        of the run; rd_correct (bulk) is the only consumer of the covariates.
        """
        input:
            segment_bed=segment_bed,
            reference=config["reference"],
            genome_size=config["genome_size"],
            mappability_bed=config["mappability_bed"] or [],
            bedgraphs=(
                [
                    _repli_cache + f"/{n}.{_repli_target}.bedGraph"
                    for n in _repli_names
                ]
                if do_repliseq
                else []
            ),
        output:
            window_bed=window_bed,
        log:
            config["log_dir"] + f"/build_window_bed/build_window_bed.{_run_id}.log",
        benchmark:
            config["bench_dir"]
            + f"/build_window_bed/build_window_bed.{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        params:
            chroms=chroms,
            input_nochr=input_nochr,
            window_size=window_size,
        script:
            "../scripts/build_window_bed.py"



rule window_bed_to_3bed:
    """Headerless 3-column BED (#CHR/START/END) for mosdepth --by; one bin set, all bulk assays."""
    input:
        window_bed=window_bed,
    output:
        mosdepth_bed=temp(config["aux_dir"] + "/windows.3col.bed.gz"),
    log:
        config["log_dir"] + f"/window_bed_to_3bed/window_bed_to_3bed.{_run_id}.log",
    benchmark:
        config["bench_dir"] + f"/window_bed_to_3bed/window_bed_to_3bed.{_run_id}.tsv"
    params:
        strip_chr_prefix="sed 's/^chr//' | " if input_nochr else "",
    shell:
        "gzip -dc {input.window_bed} | tail -n +2 | cut -f1-3 | {params.strip_chr_prefix}"
        "gzip -c > {output.mosdepth_bed} 2> {log}"
