"""Segment BED and window BED, built in every mode.

Last update: 2026-08-11

Rules:
- build_segment_bed: arms cut at the SV extremities, blacklist-subtracted
- [optional] repliseq_bigwig_to_bedgraph: fetch one ENCODE Repli-seq bigWig
- [optional] repliseq_liftover: lift an hg19 Repli-seq bedGraph to the run build
- [optional] build_window_bed: tile the segments, annotate GC, MAP, REPLI
- [optional] annotate_window_targets: per-window capture-target fraction
- window_bed_to_3bed: headerless 3-column BED for mosdepth --by
Globals:
- extremity_tsv: the configured SV breakpoints, empty when unset
- segment_bed, window_bed: paths, built or configured
- target_bed, window_target: the capture kit's targets and the per-window fraction,
  both empty when target_bed is unset
- build_windows, do_repliseq, window_size: build switches and the tile size
"""


rule build_segment_bed:
    """Cut the arms at the SV extremities, subtract the blacklist."""
    input:
        region_bed=region_bed,
        extremity_tsv=extremity_tsv,
        blacklist_bed=blacklist_bed,
    output:
        segment_bed=segment_bed,
    log:
        log_dir + f"/build_segment_bed.{_run_id}.log",
    benchmark:
        bench_dir + f"/build_segment_bed.{_run_id}.tsv"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/build_segment_bed.py"


if build_windows:
    if do_repliseq:

        rule repliseq_bigwig_to_bedgraph:
            """Fetch an ENCODE Repli-seq bigWig (hg19) and convert to bedGraph."""
            input:
                bigwig=lambda wc: file_input(f"{UCSC_REPLISEQ_BASE}/{wc.name}.bigWig"),
            output:
                (
                    temp(aux_dir + "/repliseq/{name}.hg19.bedGraph")
                    if reference_version in REPLI_LIFTOVER
                    else aux_dir + "/repliseq/{name}.hg19.bedGraph"
                ),
            log:
                log_dir
                + f"/repliseq_bigwig_to_bedgraph/repliseq_bigwig_to_bedgraph.{{name}}.{_run_id}.log",
            benchmark:
                bench_dir
                + f"/repliseq_bigwig_to_bedgraph/repliseq_bigwig_to_bedgraph.{{name}}.{_run_id}.tsv"
            wildcard_constraints:
                name="[A-Za-z0-9]+",
            conda:
                "../envs/ucsc.yaml"
            resources:
                downloads=1,
            shell:
                "bigWigToBedGraph {input.bigwig} {output} 2> {log}"

        if reference_version in REPLI_LIFTOVER:

            rule repliseq_liftover:
                """liftOver an hg19 Repli-seq bedGraph to the run's build."""
                input:
                    bedgraph=aux_dir + "/repliseq/{name}.hg19.bedGraph",
                    chain=file_input(LIFTOVER_CHAIN_URLS[reference_version]),
                output:
                    bedgraph=aux_dir
                    + f"/repliseq/{{name}}.{reference_version}.bedGraph",
                    unmapped=temp(
                        aux_dir + f"/repliseq/{{name}}.{reference_version}.unmapped"
                    ),
                log:
                    log_dir
                    + f"/repliseq_liftover/repliseq_liftover.{{name}}.{_run_id}.log",
                benchmark:
                    bench_dir
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
        """Tile segment_bed, then annotate GC, MAP and REPLI where configured."""
        input:
            segment_bed=segment_bed,
            reference=reference,
            genome_size=genome_size,
            mappability_bed=mappability_bed,
            bedgraphs=(
                [
                    aux_dir
                    + f"/repliseq/{f.removesuffix('.bigWig')}.{reference_version}.bedGraph"
                    for f in REPLISEQ_BIGWIG_FILES
                ]
                if do_repliseq
                else []
            ),
        output:
            window_bed=window_bed,
        log:
            log_dir + f"/build_window_bed.{_run_id}.log",
        benchmark:
            bench_dir + f"/build_window_bed.{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        params:
            chroms=chr_chromosomes,
            input_nochr=input_nochr,
            window_size=window_size,
        script:
            "../scripts/build_window_bed.py"


if target_bed:

    rule annotate_window_targets:
        """Per-window capture-target bp fraction, the on/off-target split."""
        input:
            window_bed=window_bed,
            target_bed=target_bed,
        output:
            window_target=window_target,
        log:
            log_dir + f"/annotate_window_targets.{_run_id}.log",
        benchmark:
            bench_dir + f"/annotate_window_targets.{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        params:
            chroms=chr_chromosomes,
        script:
            "../scripts/annotate_window_targets.py"


rule window_bed_to_3bed:
    """Headerless 3-column BED (#CHR/START/END) for mosdepth --by; one bin set, all bulk assays."""
    input:
        window_bed=window_bed,
    output:
        mosdepth_bed=temp(aux_dir + "/windows.3col.bed.gz"),
    log:
        log_dir + f"/window_bed_to_3bed.{_run_id}.log",
    benchmark:
        bench_dir + f"/window_bed_to_3bed.{_run_id}.tsv"
    params:
        strip_chr_prefix="sed 's/^chr//' | " if input_nochr else "",
    shell:
        "gzip -dc {input.window_bed} | tail -n +2 | cut -f1-3 | {params.strip_chr_prefix}"
        "gzip -c > {output.mosdepth_bed} 2> {log}"
