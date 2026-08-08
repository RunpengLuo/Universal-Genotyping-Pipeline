##################################################
# Segment BED + window BED, both built in EVERY mode.
# Bulk counts the windows (mosdepth -> rd_correct); single-cell uses them as the
# fixed-bin skeleton for binning. Only bulk needs the GC/MAP/REPLI covariates.
# Rule flow (downloads/tool calls are Snakemake rules, only binning is Python):
#
#   . build_segment_bed                [always, every mode]
#       in : segment_bed (configured segmentation) + region_bed (arms) + blacklist_bed
#       out: segment_bed  (#CHR START END region_id[arm] seg_id[segment]); also read by
#            phase_and_concat (SNP region/seg assignment), rd_correct + combine_counts
#            (QC region overlay) -- so it runs even when windows are pre-built.
#
# The window BED is EITHER built from the segments OR consumed pre-built and verified:
#
#   build_windows [no window_bed]:
#     . repliseq_bigwig_to_bedgraph    [do_repliseq (bulk only), per bigWig {name}]
#         in : bigWig URL (storage)            out: {name}.hg19.bedGraph
#     . repliseq_liftover              [do_repliseq + hg38 target only, per {name}]
#         in : {name}.hg19.bedGraph + chain URL (storage)
#         out: {name}.hg38.bedGraph (cached under aux); hg19 runs skip this
#     . build_window_bed              [one bin set for every assay of the run]
#         in : segment_bed + genome_size, and (bulk only) reference, mappability_bed,
#              {reference_version} bedgraphs
#         out: windows.bed.gz  (#CHR START END region_id seg_id [GC] [MAP] [REPLI])
#              + qc_pdf (segment- and window-length histograms)
#         Tiling is per segment row, so no window spans two segments.
#
#   else [window_bed set]:
#     . verify_window_bed             [gates the RD path]
#         in : config["window_bed"] + segment_bed
#         out: aux/window_bed.checked -- fails when a window crosses a segment bound.
#         The file itself is read directly by the consumers (window_bed_path); its
#         region_id/seg_id are assumed to be the segment_bed ids. rd_correct filters
#         it to config["chromosomes"], so extra contigs are harmless.
#
#   . window_bed_to_3bed               [one shared file for every bulk assay]
#       in : windows.bed.gz           out: aux/windows.3col.bed.gz (mosdepth --by)
#
# Globals from parse_workflow: segment_bed, build_windows, do_repliseq, window_size.
##################################################

# GC/MAP/REPLI are inputs to rd_correct, which only bulk runs; a non-bulk window BED
# carries just the intervals + region_id/seg_id.
_rd_covariates = workflow_mode == "bulk_genotyping"

# One window BED for every assay of the run (they share the segment.bed tiling).
window_bed_path = (
    config["aux_dir"] + "/windows.bed.gz" if build_windows else config["window_bed"]
)
# gate on verify_window_bed when the window BED is supplied rather than built
window_bed_checked = (
    [] if build_windows else config["aux_dir"] + "/window_bed.checked"
)


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
# Window BED: EITHER built from the segments, OR supplied and verified against
# them (window_bed_path). The Repli-seq fetch/convert only feeds build_window_bed.
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
        """Build the window BED in one pass -> config["aux_dir"]/windows.bed.gz.

        Tiles segment_bed and assigns region_id + seg_id, then annotates the
        bias-correction covariates GC / MAP / REPLI. Every covariate is gated on its
        input being non-empty, and they are all fed only in bulk: rd_correct is their
        one consumer, so a non-bulk run gets the plain intervals. One bin set for every
        assay of the run.
        """
        input:
            region_bed=segment_bed,
            reference=config["reference"] if _rd_covariates else [],
            genome_size=config["genome_size"],
            mappability_bed=(
                (config.get("mappability_bed") or []) if _rd_covariates else []
            ),
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
            chroms=chroms,
            input_nochr=input_nochr,
            window_size=window_size,
        script:
            "../scripts/build_window_bed.py"

else:

    rule verify_window_bed:
        """Fail unless every supplied window sits inside one segment_bed segment."""
        input:
            window_bed=config["window_bed"],
            segment_bed=segment_bed,
        output:
            checked=config["aux_dir"] + "/window_bed.checked",
        log:
            config["log_dir"] + f"/verify_window_bed/verify_window_bed.{_run_id}.log",
        benchmark:
            config["bench_dir"]
            + f"/verify_window_bed/verify_window_bed.{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        params:
            chroms=chroms,
        script:
            "../scripts/verify_window_bed.py"


rule window_bed_to_3bed:
    """Headerless 3-column BED (#CHR/START/END) for mosdepth --by; one bin set, all bulk assays."""
    input:
        window_bed=window_bed_path,
        checked=window_bed_checked,
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
