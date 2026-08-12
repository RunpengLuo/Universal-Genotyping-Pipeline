"""Segment BED and window BED, built in every mode.

Last update: 2026-08-11

Rules:
- build_segment_bed: arm-stamped, blacklist-subtracted segments
- [optional] repliseq_bigwig_to_bedgraph: fetch one ENCODE Repli-seq bigWig
- [optional] repliseq_liftover: lift an hg19 Repli-seq bedGraph to hg38
- [optional] build_window_bed: tile the segments, annotate GC, MAP, REPLI
- window_bed_to_3bed: headerless 3-column BED for mosdepth --by
Globals:
- input_segment_bed: the configured segmentation, region_bed when unset
- segment_bed, window_bed: paths, built or configured
- build_windows, do_repliseq, window_size: build switches and the tile size
"""


rule build_segment_bed:
    """Stamp each segment with its arm, subtract the blacklist."""
    input:
        segments=input_segment_bed,
        region_bed=region_bed,
        blacklist_bed=blacklist_bed,
    output:
        segment_bed=segment_bed,
    log:
        log_dir + f"/build_segment_bed/build_segment_bed.{_run_id}.log",
    benchmark:
        bench_dir + f"/build_segment_bed/build_segment_bed.{_run_id}.tsv"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/build_segment_bed.py"


if build_windows:
    if do_repliseq:
        _repli_cache = aux_dir + "/repliseq"
        _repli_names = [f[: -len(".bigWig")] for f in REPLISEQ_BIGWIG_FILES]
        _repli_target = reference_version
        _repli_lift = _repli_target != "hg19"

        rule repliseq_bigwig_to_bedgraph:
            """Fetch an ENCODE Repli-seq bigWig (hg19) and convert to bedGraph."""
            input:
                bigwig=lambda wc: file_input(f"{UCSC_REPLISEQ_BASE}/{wc.name}.bigWig"),
            output:
                (
                    temp(_repli_cache + "/{name}.hg19.bedGraph")
                    if _repli_lift
                    else _repli_cache + "/{name}.hg19.bedGraph"
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
                [_repli_cache + f"/{n}.{_repli_target}.bedGraph" for n in _repli_names]
                if do_repliseq
                else []
            ),
        output:
            window_bed=window_bed,
        log:
            log_dir + f"/build_window_bed/build_window_bed.{_run_id}.log",
        benchmark:
            bench_dir + f"/build_window_bed/build_window_bed.{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        params:
            chroms=chr_chromosomes,
            input_nochr=input_nochr,
            window_size=window_size,
        script:
            "../scripts/build_window_bed.py"


rule window_bed_to_3bed:
    """Headerless 3-column BED (#CHR/START/END) for mosdepth --by; one bin set, all bulk assays."""
    input:
        window_bed=window_bed,
    output:
        mosdepth_bed=temp(aux_dir + "/windows.3col.bed.gz"),
    log:
        log_dir + f"/window_bed_to_3bed/window_bed_to_3bed.{_run_id}.log",
    benchmark:
        bench_dir + f"/window_bed_to_3bed/window_bed_to_3bed.{_run_id}.tsv"
    params:
        strip_chr_prefix="sed 's/^chr//' | " if input_nochr else "",
    shell:
        "gzip -dc {input.window_bed} | tail -n +2 | cut -f1-3 | {params.strip_chr_prefix}"
        "gzip -c > {output.mosdepth_bed} 2> {log}"
