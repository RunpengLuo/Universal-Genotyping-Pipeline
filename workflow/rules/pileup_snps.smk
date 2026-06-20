rule pileup_snps_bulk_mode1b:
    input:
        bam=lambda wc: get_data[(wc.assay_type, wc.rep_id)][1],
        snp_vcf=lambda wc: branch(
            run_genotyping,
            then=config["phase_dir"] + "/phased_het_snps.vcf.gz",
            otherwise=config["het_snp_vcf"],
        ),
    output:
        out_dir=directory(config["pileup_dir"] + "/{assay_type}_{rep_id}/"),
        out_vcf=config["pileup_dir"] + "/{assay_type}_{rep_id}/cellSNP.base.vcf.gz",
        out_tsv=config["pileup_dir"] + "/{assay_type}_{rep_id}/cellSNP.samples.tsv",
        out_dp=config["pileup_dir"] + "/{assay_type}_{rep_id}/cellSNP.tag.DP.mtx",
        out_ad=config["pileup_dir"] + "/{assay_type}_{rep_id}/cellSNP.tag.AD.mtx",
    wildcard_constraints:
        assay_type="(bulkWGS|bulkWGS-lr|bulkWES)",
    threads: config["threads"]["pileup"]
    params:
        minMAF=config["params_cellsnp_lite"]["minMAF_pileup"],
        minCOUNT=config["params_cellsnp_lite"]["minCOUNT_pileup"],
    log:
        config["log_dir"]
        + f"/pileup_snps_bulk_mode1b/pileup_snps.{{assay_type}}_{{rep_id}}.{_run_id}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        r"""
        cellsnp-lite \
            -s "{input.bam}" \
            -R "{input.snp_vcf}" \
            -O "{output.out_dir}" \
            -p {threads} \
            --minMAF {params.minMAF} \
            --minCOUNT {params.minCOUNT} \
            --UMItag None \
            --cellTAG None \
            --gzip > {log} 2>&1
        """


rule pileup_snps_single_cell_mode1a:
    input:
        barcode=lambda wc: get_data[(wc.assay_type, wc.rep_id)][0],
        bam=lambda wc: get_data[(wc.assay_type, wc.rep_id)][1],
        snp_vcf=lambda wc: branch(
            run_genotyping,
            then=config["phase_dir"] + "/phased_het_snps.vcf.gz",
            otherwise=config["het_snp_vcf"],
        ),
    output:
        out_dir=directory(config["pileup_dir"] + "/{assay_type}_{rep_id}/"),
        out_vcf=config["pileup_dir"] + "/{assay_type}_{rep_id}/cellSNP.base.vcf.gz",
        out_tsv=config["pileup_dir"] + "/{assay_type}_{rep_id}/cellSNP.samples.tsv",
        out_dp=config["pileup_dir"] + "/{assay_type}_{rep_id}/cellSNP.tag.DP.mtx",
        out_ad=config["pileup_dir"] + "/{assay_type}_{rep_id}/cellSNP.tag.AD.mtx",
    wildcard_constraints:
        assay_type="(scRNA|scATAC|VISIUM|VISIUM3prime)",
    threads: config["threads"]["pileup"]
    params:
        UMItag=lambda wc: branch(
            wc.assay_type == "scATAC",
            then="None",
            otherwise=config["params_cellsnp_lite"]["UMItag"],
        ),
        cellTAG=config["params_cellsnp_lite"]["cellTAG"],
        minMAF=config["params_cellsnp_lite"]["minMAF_pileup"],
        minCOUNT=config["params_cellsnp_lite"]["minCOUNT_pileup"],
    log:
        config["log_dir"]
        + f"/pileup_snps_single_cell_mode1a/pileup_snps.{{assay_type}}_{{rep_id}}.{_run_id}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        r"""
        cellsnp-lite \
            -b "{input.barcode}" \
            -s "{input.bam}" \
            -R "{input.snp_vcf}" \
            -O "{output.out_dir}" \
            -p {threads} \
            --minMAF {params.minMAF} \
            --minCOUNT {params.minCOUNT} \
            --UMItag {params.UMItag} \
            --cellTAG {params.cellTAG} \
            --gzip > {log} 2>&1
        """
