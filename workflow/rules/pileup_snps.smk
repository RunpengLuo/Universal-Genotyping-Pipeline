rule pileup_snps_bulk_bcftools:
    """Bulk het-SNP read counting with bcftools (REF/ALT allele depths at the phased loci)."""
    input:
        alignment=lambda wc: bam_stream_input(get_data[(wc.assay_type, wc.dataset_id)]),
        alignment_index=lambda wc: bam_stream_index_input(
            get_data[(wc.assay_type, wc.dataset_id)]
        ),
        snp_vcf=phased_snp_vcf,
        reference=config["reference"],
    output:
        counts=config["pileup_dir"]
        + "/{assay_type}_{dataset_id}/bcftools.counts.tsv.gz",
    log:
        config["log_dir"]
        + f"/pileup_snps_bulk_bcftools/pileup_snps_bulk_bcftools.{{assay_type}}_{{dataset_id}}.{_run_id}.log",
    benchmark:
        config["bench_dir"]
        + f"/pileup_snps_bulk_bcftools/pileup_snps_bulk_bcftools.{{assay_type}}_{{dataset_id}}.{_run_id}.tsv"
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
        region_arg=(
            "-r " + ",".join(raw_chrom(c) for c in config["chromosomes"])
            if remote_stream
            else ""
        ),
    shell:
        r"""
        set -euo pipefail
        ALN="{input.alignment}"; [ -z "$ALN" ] && ALN="{params.bam_arg}"
        (
          bcftools mpileup "$ALN" \
              -f "{input.reference}" \
              -Ou \
              --threads {threads} \
              -a FORMAT/AD \
              --skip-indels \
              -q {params.min_mapq} \
              -Q {params.min_baseq} \
              -d {params.max_depth} \
              {params.extra_params} \
              {params.region_arg} \
              -T "{input.snp_vcf}" \
          | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' \
          | bgzip -c > {output.counts}
        ) 2> {log}
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
        out_dir=directory(config["pileup_dir"] + "/{assay_type}_{dataset_id}/"),
        out_vcf=config["pileup_dir"] + "/{assay_type}_{dataset_id}/cellSNP.base.vcf.gz",
        out_tsv=config["pileup_dir"] + "/{assay_type}_{dataset_id}/cellSNP.samples.tsv",
        out_dp=config["pileup_dir"] + "/{assay_type}_{dataset_id}/cellSNP.tag.DP.mtx",
        out_ad=config["pileup_dir"] + "/{assay_type}_{dataset_id}/cellSNP.tag.AD.mtx",
    log:
        config["log_dir"]
        + f"/pileup_snps_nonbulk_mode1a/pileup_snps_nonbulk_mode1a.{{assay_type}}_{{dataset_id}}.{_run_id}.log",
    benchmark:
        config["bench_dir"]
        + f"/pileup_snps_nonbulk_mode1a/pileup_snps_nonbulk_mode1a.{{assay_type}}_{{dataset_id}}.{_run_id}.tsv"
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
