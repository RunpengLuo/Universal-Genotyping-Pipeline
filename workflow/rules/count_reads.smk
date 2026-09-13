"""Fixed-bin read depth and its bias correction, bulk only.

Last update: 2026-08-27

Rules:
- [optional] run_mosdepth: whole-genome per-bin depth for one replicate
- [optional] run_mosdepth_chrom, merge_mosdepth: the remote_mode stream variants
- count_read_starts_chrom, merge_read_starts: per-window read counts by alignment start
- rd_correct: GC, mappability and replication-timing correction, once for all bulk
Outputs:
- pileup_dir/bulk/window.raw.dp.npz: raw depth, windows x all bulk datasets
- pileup_dir/bulk/window.dp.npz: corrected depth, windows x all bulk datasets
- pileup_dir/bulk/depth_statistics.tsv: per-dataset depth summary
- pileup_dir/{assay_type}/{dataset_id}.rdcount.bed.gz: read-start counts, mosdepth's
  4-column shape. No rule reads it yet; combine_counts only holds it in the DAG.
"""

if remote_mode != "stream":

    rule run_mosdepth:
        input:
            alignment=lambda wc: alignment_input(
                get_data[(wc.assay_type, wc.dataset_id)]
            ),
            alignment_index=lambda wc: alignment_index_input(
                get_data[(wc.assay_type, wc.dataset_id)]
            ),
            windows_bed=aux_dir + "/windows.3col.bed.gz",
        output:
            mosdepth_file=pileup_dir
            + "/{assay_type}/out_mosdepth/{dataset_id}.regions.bed.gz",
        log:
            log_dir
            + f"/run_mosdepth/run_mosdepth.{{assay_type}}_{{dataset_id}}.{_run_id}.log",
        benchmark:
            bench_dir
            + f"/run_mosdepth/run_mosdepth.{{assay_type}}_{{dataset_id}}.{_run_id}.tsv"
        conda:
            "../envs/mosdepth.yaml"
        threads: config["threads"]["mosdepth"]
        resources:
            downloads=lambda wc: download_slots(
                get_data[(wc.assay_type, wc.dataset_id)]
            ),
        params:
            out_prefix=pileup_dir + "/{assay_type}/out_mosdepth/{dataset_id}",
            read_quality=config["params_count_reads"]["read_quality"],
            exclude_flags=config["params_count_reads"]["exclude_flags"],
            extra_params=config["params_count_reads"]["mosdepth_extra_params"],
        shell:
            r"""
            mosdepth \
                -t {threads} \
                -Q {params.read_quality} \
                -F {params.exclude_flags} \
                --by {input.windows_bed} \
                {params.extra_params} \
                {params.out_prefix} {input.alignment} > {log} 2>&1
            """

else:

    rule run_mosdepth_chrom:
        """Stream mode: per-chromosome mosdepth (index jump), merged by merge_mosdepth."""
        input:
            alignment=lambda wc: bam_stream_input(
                get_data[(wc.assay_type, wc.dataset_id)]
            ),
            alignment_index=lambda wc: bam_stream_index_input(
                get_data[(wc.assay_type, wc.dataset_id)]
            ),
            windows_bed=aux_dir + "/windows.3col.bed.gz",
            reference=reference,
        output:
            mosdepth_file=temp(
                pileup_dir
                + "/{assay_type}/out_mosdepth/{dataset_id}.chr{chrname}.regions.bed.gz"
            ),
        log:
            log_dir
            + f"/run_mosdepth/run_mosdepth.{{assay_type}}_{{dataset_id}}.chr{{chrname}}.{_run_id}.log",
        benchmark:
            bench_dir
            + f"/run_mosdepth/run_mosdepth.{{assay_type}}_{{dataset_id}}.chr{{chrname}}.{_run_id}.tsv"
        conda:
            "../envs/mosdepth.yaml"
        threads: config["threads"]["mosdepth"]
        resources:
            downloads=lambda wc: download_slots(
                get_data[(wc.assay_type, wc.dataset_id)]
            ),
        params:
            out_prefix=pileup_dir
            + "/{assay_type}/out_mosdepth/{dataset_id}.chr{chrname}",
            chrom=lambda wc: input_chrom(wc.chrname),
            read_quality=config["params_count_reads"]["read_quality"],
            exclude_flags=config["params_count_reads"]["exclude_flags"],
            extra_params=config["params_count_reads"]["mosdepth_extra_params"],
            bam_arg=lambda wc: bam_stream_arg(get_data[(wc.assay_type, wc.dataset_id)]),
        shell:
            r"""
            set -euo pipefail
            ALN="{input.alignment}"; [ -z "$ALN" ] && ALN="{params.bam_arg}"
            mosdepth \
                -t {threads} \
                -Q {params.read_quality} \
                -F {params.exclude_flags} \
                -c {params.chrom} \
                -f {input.reference} \
                --by {input.windows_bed} \
                {params.extra_params} \
                {params.out_prefix} "$ALN" > {log} 2>&1
            """

    rule merge_mosdepth:
        """Concat per-chrom regions (config-chrom order) -> the file rd_correct reads."""
        input:
            per_chrom=lambda wc: [
                pileup_dir
                + f"/{wc.assay_type}/out_mosdepth/{wc.dataset_id}.chr{c}.regions.bed.gz"
                for c in nochr_chromosomes
            ],
        output:
            mosdepth_file=pileup_dir
            + "/{assay_type}/out_mosdepth/{dataset_id}.regions.bed.gz",
        log:
            log_dir
            + f"/run_mosdepth/merge_mosdepth.{{assay_type}}_{{dataset_id}}.{_run_id}.log",
        shell:
            r"""
            cat {input.per_chrom} > {output.mosdepth_file} 2> {log}
            """


rule count_read_starts_chrom:
    """Reads whose leftmost mapped base falls in each window, one chromosome (index jump).

    A read is charged to the single window holding its SAM POS, so it is counted once
    however many windows it spans - unlike mosdepth, which spreads its bases over all of
    them. Per chromosome in every remote_mode: `bedtools -sorted` sweeps, so A and B must
    agree on chromosome order, and the window BED's genomic sort need not match the
    alignment's @SQ order. One chromosome each side makes that constraint vacuous.
    """
    input:
        alignment=lambda wc: bam_stream_input(get_data[(wc.assay_type, wc.dataset_id)]),
        alignment_index=lambda wc: bam_stream_index_input(
            get_data[(wc.assay_type, wc.dataset_id)]
        ),
        chrom_bed=aux_dir + "/windows.3col.chr{chrname}.bed.gz",
        reference=reference,
    output:
        rdcount=temp(
            pileup_dir + "/{assay_type}/out_rdcount/{dataset_id}.chr{chrname}.bed.gz"
        ),
    log:
        log_dir
        + f"/count_read_starts/count_read_starts.{{assay_type}}_{{dataset_id}}.chr{{chrname}}.{_run_id}.log",
    benchmark:
        bench_dir
        + f"/count_read_starts/count_read_starts.{{assay_type}}_{{dataset_id}}.chr{{chrname}}.{_run_id}.tsv"
    conda:
        "../envs/bedtools.yaml"
    threads: config["threads"]["mosdepth"]
    resources:
        downloads=lambda wc: download_slots(get_data[(wc.assay_type, wc.dataset_id)]),
    params:
        chrom=lambda wc: input_chrom(wc.chrname),
        out_chrom=lambda wc: f"chr{wc.chrname}",
        read_quality=config["params_count_reads"]["read_quality"],
        exclude_flags=config["params_count_reads"]["exclude_flags"],
        bam_arg=lambda wc: bam_stream_arg(get_data[(wc.assay_type, wc.dataset_id)]),
    shell:
        r"""
        set -euo pipefail
        ALN="{input.alignment}"; [ -z "$ALN" ] && ALN="{params.bam_arg}"
        (
          samtools view -@ {threads} \
              -T "{input.reference}" \
              -q {params.read_quality} \
              -F {params.exclude_flags} \
              "$ALN" "{params.chrom}" \
          | cut -f4 \
          | awk -v OFS='\t' -v C="{params.out_chrom}" '{{print C, $1 - 1, $1}}' \
          | bedtools intersect -a "{input.chrom_bed}" -b stdin -c -sorted \
          | gzip -c > "{output.rdcount}"
        ) 2> {log}
        """


rule merge_read_starts:
    """Concat per-chrom read-start counts (config-chrom order), one file per dataset."""
    input:
        per_chrom=lambda wc: [
            pileup_dir + f"/{wc.assay_type}/out_rdcount/{wc.dataset_id}.chr{c}.bed.gz"
            for c in nochr_chromosomes
        ],
    output:
        rdcount=pileup_dir + "/{assay_type}/{dataset_id}.rdcount.bed.gz",
    log:
        log_dir
        + f"/count_read_starts/merge_read_starts.{{assay_type}}_{{dataset_id}}.{_run_id}.log",
    threads: 1
    shell:
        r"""
        cat {input.per_chrom} > {output.rdcount} 2> {log}
        """


rule rd_correct:
    """Per-window sequencing bias correction (GC/mappability/replication timing)."""
    input:
        mosdepth_files=[
            pileup_dir + f"/{at}/out_mosdepth/{rid}.regions.bed.gz"
            for at in assay_types
            for rid in assay2dataset_ids[at]
        ],
        window_bed=window_bed,
        window_target=window_target,
        genome_size=genome_size,
        region_bed=segment_bed,
        blacklist_bed=blacklist_bed,
    output:
        dp_raw=pileup_dir + "/bulk/window.raw.dp.npz",
        dp_corrected=pileup_dir + "/bulk/window.dp.npz",
        depth_stats=report(
            pileup_dir + "/bulk/depth_statistics.tsv",
            category="QC stats",
            subcategory="depth",
            labels={"table": "depth statistics"},
        ),
        qc_pdf=report(
            qc_dir + "/rd_correction.bulk.pdf",
            category="QC plots",
            subcategory="read-depth correction",
        ),
    log:
        log_dir + f"/rd_correct.bulk.{_run_id}.log",
    benchmark:
        bench_dir + f"/rd_correct.bulk.{_run_id}.tsv"
    conda:
        "../envs/base.yaml"
    params:
        qc_dir=qc_dir,
        sample_id=sample_id,
        dataset_ids=[rid for at in assay_types for rid in assay2dataset_ids[at]],
        dataset_assays=[at for at in assay_types for rid in assay2dataset_ids[at]],
        sample_types=[st for at in assay_types for st in assay2sample_types[at]],
        chroms=chr_chromosomes,
        samplesize=config["params_count_reads"]["samplesize"],
        routlier=config["params_count_reads"]["routlier"],
        doutlier=config["params_count_reads"]["doutlier"],
        min_mappability=config["params_count_reads"]["min_mappability"],
        gc_correct=config["params_count_reads"]["gc_correct"],
        rd_correct_method=config["params_count_reads"]["rd_correct_method"],
        rt_correct=config["params_count_reads"]["rt_correct"],
        run_id=_run_id,
    script:
        """../scripts/rd_correct.py"""
