"""Count reads at the phased het SNPs, once per replicate.

Last update: 2026-08-11

Rules:
- pileup_snps_bulk_bcftools_chrom: bcftools REF and ALT depths, one chromosome of one
  bulk replicate, so a replicate fans out instead of walking the whole alignment once
- merge_pileup_counts: concat the per-chrom counts in config-chrom order. mpileup emits
  regions in that order, so the merge is byte-identical to a single whole-genome job
- pileup_snps_nonbulk_mode1a: cellsnp-lite per-cell counts per replicate
Outputs:
- pileup_dir/{assay}_{dataset_id}/bcftools.counts.tsv.gz: bulk allele depths
- pileup_dir/{assay}_{dataset_id}/cellSNP.*: single-cell counts and barcodes
"""


rule pileup_snps_bulk_bcftools_chrom:
    """One chromosome of a bulk replicate's het-SNP REF/ALT depths (index jump)."""
    input:
        alignment=lambda wc: bam_stream_input(get_data[(wc.assay_type, wc.dataset_id)]),
        alignment_index=lambda wc: bam_stream_index_input(
            get_data[(wc.assay_type, wc.dataset_id)]
        ),
        snp_vcf=phased_snp_vcf,
        reference=reference,
    output:
        counts=temp(
            pileup_dir
            + "/{assay_type}_{dataset_id}/bcftools.counts.chr{chrname}.tsv.gz"
        ),
    log:
        log_dir
        + f"/pileup_snps/bcftools.{{assay_type}}_{{dataset_id}}.chr{{chrname}}.{_run_id}.log",
    benchmark:
        bench_dir
        + f"/pileup_snps/bcftools.{{assay_type}}_{{dataset_id}}.chr{{chrname}}.{_run_id}.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: config["threads"]["pileup"]
    resources:
        downloads=lambda wc: download_slots(get_data[(wc.assay_type, wc.dataset_id)]),
    params:
        min_mapq=config["params_bcftools"]["min_mapq"],
        min_baseq=config["params_bcftools"]["min_baseq"],
        max_depth=config["params_bcftools"]["max_depth"],
        extra_params=config["params_bcftools"]["extra_params"],
        bam_arg=lambda wc: bam_stream_arg(get_data[(wc.assay_type, wc.dataset_id)]),
        chrom=lambda wc: input_chrom(wc.chrname),
    shell:
        r"""
        set -euo pipefail
        ALN="{input.alignment}"; [ -z "$ALN" ] && ALN="{params.bam_arg}"
        (
          bcftools mpileup "$ALN" \
              --fasta-ref "{input.reference}" \
              --output-type u \
              --annotate FORMAT/AD \
              --skip-indels \
              --min-MQ {params.min_mapq} \
              --min-BQ {params.min_baseq} \
              --max-depth {params.max_depth} \
              {params.extra_params} \
              --regions {params.chrom} \
              --targets-file "{input.snp_vcf}" \
          | bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' \
          | bgzip -@ {threads} -c > {output.counts}
        ) 2> {log}
        """


rule merge_pileup_counts:
    """Concat per-chrom counts (config-chrom order) -> the file phase_and_concat reads."""
    input:
        per_chrom=lambda wc: [
            pileup_dir
            + f"/{wc.assay_type}_{wc.dataset_id}/bcftools.counts.chr{c}.tsv.gz"
            for c in nochr_chromosomes
        ],
    output:
        counts=pileup_dir + "/{assay_type}_{dataset_id}/bcftools.counts.tsv.gz",
    log:
        log_dir
        + f"/pileup_snps/merge_pileup_counts.{{assay_type}}_{{dataset_id}}.{_run_id}.log",
    threads: 1
    shell:
        r"""
        cat {input.per_chrom} > {output.counts} 2> {log}
        """


rule pileup_snps_nonbulk_mode1a:
    input:
        barcode=lambda wc: file_input(
            get_data[(wc.assay_type, wc.dataset_id)]["barcodes"]
        ),
        alignment=lambda wc: alignment_input(get_data[(wc.assay_type, wc.dataset_id)]),
        alignment_index=lambda wc: alignment_index_input(
            get_data[(wc.assay_type, wc.dataset_id)]
        ),
        snp_vcf=phased_snp_vcf,
    output:
        out_dir=directory(pileup_dir + "/{assay_type}_{dataset_id}/"),
        out_vcf=pileup_dir + "/{assay_type}_{dataset_id}/cellSNP.base.vcf.gz",
        out_tsv=pileup_dir + "/{assay_type}_{dataset_id}/cellSNP.samples.tsv",
        out_dp=pileup_dir + "/{assay_type}_{dataset_id}/cellSNP.tag.DP.mtx",
        out_ad=pileup_dir + "/{assay_type}_{dataset_id}/cellSNP.tag.AD.mtx",
    log:
        log_dir
        + f"/pileup_snps/cellsnp_lite_mode1a.{{assay_type}}_{{dataset_id}}.{_run_id}.log",
    benchmark:
        bench_dir
        + f"/pileup_snps/cellsnp_lite_mode1a.{{assay_type}}_{{dataset_id}}.{_run_id}.tsv"
    wildcard_constraints:
        assay_type="(scRNA|scATAC|VISIUM|VISIUM3prime)",
    conda:
        "../envs/cellsnp.yaml"
    threads: config["threads"]["pileup"]
    resources:
        downloads=lambda wc: download_slots(get_data[(wc.assay_type, wc.dataset_id)]),
    params:
        UMItag=lambda wc: branch(
            wc.assay_type == "scATAC",
            then="None",
            otherwise=config["params_cellsnp_lite"]["UMItag"],
        ),
        cellTAG=config["params_cellsnp_lite"]["cellTAG"],
        minMAF=config["params_cellsnp_lite"]["minMAF_pileup"],
        minCOUNT=config["params_cellsnp_lite"]["minCOUNT_pileup"],
    shell:
        r"""
        cellsnp-lite \
            -b "{input.barcode}" \
            -s "{input.alignment}" \
            -R "{input.snp_vcf}" \
            -O "{output.out_dir}" \
            -p {threads} \
            --minMAF {params.minMAF} \
            --minCOUNT {params.minCOUNT} \
            --UMItag {params.UMItag} \
            --cellTAG {params.cellTAG} \
            --gzip > {log} 2>&1
        """
