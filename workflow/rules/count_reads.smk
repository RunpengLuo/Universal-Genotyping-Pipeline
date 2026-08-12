"""Fixed-bin read depth and its bias correction, bulk only.

Last update: 2026-08-12

Rules:
- [optional] run_mosdepth: whole-genome per-bin depth for one replicate
- [optional] run_mosdepth_chrom, merge_mosdepth: the remote_mode stream variants
- rd_correct: GC, mappability and replication-timing correction, once for all bulk
Outputs:
- pileup_dir/bulk/window.dp.npz: corrected depth, windows x all bulk datasets
- pileup_dir/bulk/depth_statistics.tsv: per-dataset depth summary
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
            read_quality=config["params_mosdepth"]["read_quality"],
            extra_params=config["params_mosdepth"]["extra_params"],
            bam_arg=lambda wc: bam_stream_arg(get_data[(wc.assay_type, wc.dataset_id)]),
        shell:
            r"""
            set -euo pipefail
            ALN="{input.alignment}"; [ -z "$ALN" ] && ALN="{params.bam_arg}"
            mosdepth \
                -t {threads} \
                -Q {params.read_quality} \
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


rule rd_correct:
    """Per-window sequencing bias correction (GC/mappability/replication timing)."""
    input:
        mosdepth_files=[
            pileup_dir + f"/{at}/out_mosdepth/{rid}.regions.bed.gz"
            for at in assay_types
            for rid in assay2dataset_ids[at]
        ],
        window_bed=window_bed,
        genome_size=genome_size,
        region_bed=segment_bed,
        blacklist_bed=blacklist_bed,
    output:
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
        log_dir + f"/rd_correct/rd_correct.bulk.{_run_id}.log",
    benchmark:
        bench_dir + f"/rd_correct/rd_correct.bulk.{_run_id}.tsv"
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
        gc_correct_method=config["params_count_reads"]["gc_correct_method"],
        rt_correct=config["params_count_reads"]["rt_correct"],
        run_id=_run_id,
    script:
        """../scripts/rd_correct.py"""
