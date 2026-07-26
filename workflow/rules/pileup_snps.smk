rule pileup_snps_bulk_mode1b:
    input:
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
        + f"/pileup_snps_bulk_mode1b/pileup_snps_bulk_mode1b.{{assay_type}}_{{dataset_id}}.{_run_id}.log",
    benchmark:
        config["bench_dir"]
        + f"/pileup_snps_bulk_mode1b/pileup_snps_bulk_mode1b.{{assay_type}}_{{dataset_id}}.{_run_id}.tsv"
    wildcard_constraints:
        assay_type="(bulkWGS|bulkWGS-lr|bulkWES)",
    conda:
        "../envs/cellsnp.yaml"
    threads: config["threads"]["pileup"]
    resources:
        downloads=lambda wc: download_slots(get_data[(wc.assay_type, wc.dataset_id)]),
    params:
        minMAF=config["params_cellsnp_lite"]["minMAF_pileup"],
        minCOUNT=config["params_cellsnp_lite"]["minCOUNT_pileup"],
    shell:
        r"""
        cellsnp-lite \
            -s "{input.alignment}" \
            -R "{input.snp_vcf}" \
            -O "{output.out_dir}" \
            -p {threads} \
            --minMAF {params.minMAF} \
            --minCOUNT {params.minCOUNT} \
            --UMItag None \
            --cellTAG None \
            --gzip > {log} 2>&1
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
