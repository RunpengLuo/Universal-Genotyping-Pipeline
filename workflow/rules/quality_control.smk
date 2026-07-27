##################################################
# Quality-control rules
##################################################


if workflow_mode == "bulk_genotyping" and run_genotyping:

    rule plot_genotype_qc:
        """Het vs hom-alt allele-fraction QC from the genotyped VCFs (genotyping sample)."""
        input:
            vcfs=expand(
                config["snp_dir"] + "/chr{chrname}.vcf.gz",
                chrname=config["chromosomes"],
            ),
            genome_size=config["genome_size"],
        output:
            qc_pdf=report(
                config["qc_dir"] + "/genotype_snp_qc.pdf",
                category="QC plots",
                subcategory="genotyping",
                labels={"plot": "genotype het/hom-alt AF"},
            ),
        log:
            config["log_dir"] + f"/plot_genotype_qc/plot_genotype_qc.{_run_id}.log",
        benchmark:
            config["bench_dir"] + f"/plot_genotype_qc/plot_genotype_qc.{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        threads: 1
        script:
            "../scripts/plot/plot_genotype_qc.py"
