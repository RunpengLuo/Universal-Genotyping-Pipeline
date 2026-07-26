if run_phasing and config["phaser"] == "shapeit":

    rule phase_snps_shapeit:
        input:
            snp_vcf=lambda wc: config["snp_dir"] + f"/chr{wc.chrname}.vcf.gz",
            phasing_panel_file=lambda wc: get_phasing_panel(wc.chrname),
            gmap_file=lambda wc: get_genetic_map(wc.chrname),
        output:
            phased_file=config["phase_dir"] + "/chr{chrname}.vcf.gz",
            bcf_file=temp(config["phase_dir"] + "/chr{chrname}.bcf"),
            bcf_file_csi=temp(config["phase_dir"] + "/chr{chrname}.bcf.csi"),
        log:
            config["log_dir"]
            + f"/phase_snps_shapeit/phase_snps.chr{{chrname}}.{_run_id}.log",
        benchmark:
            config["bench_dir"]
            + f"/phase_snps_shapeit/phase_snps.chr{{chrname}}.{_run_id}.tsv"
        conda:
            "../envs/shapeit.yaml"
        threads: config["threads"]["phase"]
        params:
            chrom="chr{chrname}",
        shell:
            r"""
            phase_common \
                --input "{input.snp_vcf}" \
                --map "{input.gmap_file}" \
                --reference "{input.phasing_panel_file}" \
                --region "{params.chrom}" \
                --thread "{threads}" \
                --output "{output.bcf_file}" > {log} 2>&1

            bcftools view -Ov "{output.bcf_file}" | bgzip > "{output.phased_file}"
            tabix -f -p vcf "{output.phased_file}"
            """


if run_phasing and config["phaser"] == "eagle":

    rule phase_snps_eagle:
        input:
            snp_vcf=lambda wc: config["snp_dir"] + f"/chr{wc.chrname}.vcf.gz",
            phasing_panel_file=lambda wc: get_phasing_panel(wc.chrname),
            gmap_file=lambda wc: get_genetic_map(wc.chrname),
        output:
            phased_file=config["phase_dir"] + "/chr{chrname}.vcf.gz",
        log:
            config["log_dir"]
            + f"/phase_snps_eagle/phase_snps.chr{{chrname}}.{_run_id}.log",
        benchmark:
            config["bench_dir"]
            + f"/phase_snps_eagle/phase_snps.chr{{chrname}}.{_run_id}.tsv"
        conda:
            "../envs/eagle.yaml"
        threads: config["threads"]["phase"]
        params:
            chrom="chr{chrname}",
            out_prefix=config["phase_dir"] + "/chr{chrname}",
        shell:
            r"""
            eagle \
                --vcfTarget "{input.snp_vcf}" \
                --geneticMapFile "{input.gmap_file}" \
                --vcfRef "{input.phasing_panel_file}" \
                --vcfOutFormat z \
                --numThreads "{threads}" \
                --outPrefix {params.out_prefix} > {log} 2>&1
            tabix -f -p vcf "{output.phased_file}"
            """


if run_phasing and config["phaser"] == "longphase":

    rule phase_snps_longphase:
        input:
            snp_vcf=lambda wc: config["snp_dir"] + f"/chr{wc.chrname}.vcf.gz",
            alignment=lambda wc: alignment_input(phase_files),
            alignment_index=lambda wc: alignment_index_input(phase_files),
            reference=lambda wc: config["reference"],
        output:
            phased_file=config["phase_dir"] + "/chr{chrname}.vcf.gz",
        log:
            config["log_dir"]
            + f"/phase_snps_longphase/phase_snps.chr{{chrname}}.{_run_id}.log",
        benchmark:
            config["bench_dir"]
            + f"/phase_snps_longphase/phase_snps.chr{{chrname}}.{_run_id}.tsv"
        conda:
            "../envs/longphase.yaml"
        threads: config["threads"]["phase"]
        resources:
            downloads=download_slots(phase_files),
        params:
            chrom="chr{chrname}",
            min_mapq=config["params_longphase"]["min_mapq"],
            extra_params=config["params_longphase"]["extra_params"],
            out_prefix=config["phase_dir"] + "/chr{chrname}",
        shell:
            r"""
            longphase phase \
                --bam-file={input.alignment} \
                --reference={input.reference} \
                --snp-file={input.snp_vcf} \
                --mappingQuality={params.min_mapq} \
                --out-prefix={params.out_prefix} \
                --threads={threads} \
                {params.extra_params} > {log} 2>&1
            bgzip -f "{params.out_prefix}.vcf"
            tabix -f -p vcf "{output.phased_file}"
            """


if run_phasing:

    rule concat_and_extract_phased_het_snps:
        input:
            vcf_files=expand(
                config["phase_dir"] + "/chr{chrname}.vcf.gz",
                chrname=config["chromosomes"],
            ),
        output:
            phased_vcf=config["phase_dir"] + "/phased_het_snps.vcf.gz",
            phased_vcf_tbi=config["phase_dir"] + "/phased_het_snps.vcf.gz.tbi",
            snp_stats=report(
                config["phase_dir"] + "/germline_snp_statistics.tsv",
                category="QC stats",
                subcategory="phasing",
                labels={"table": "germline SNP statistics"},
            ),
            lst_file=temp(config["phase_dir"] + "/phased_snps.lst"),
        log:
            config["log_dir"] + f"/concat_and_extract_phased_het_snps/{_run_id}.log",
        benchmark:
            config["bench_dir"] + f"/concat_and_extract_phased_het_snps/{_run_id}.tsv"
        conda:
            "../envs/bcftools.yaml"
        threads: 1
        shell:
            r"""
            printf "#CHR\ttotal\thet_phased\thet_unphased\thom_alt\thom_ref\n" > "{output.snp_stats}"
            for vcf in {input.vcf_files}; do
                chr=$(basename "$vcf" .vcf.gz)
                n_het_phased=$(bcftools view -H -i 'GT="0|1" || GT="1|0"' "$vcf" 2>/dev/null | wc -l)
                n_het_unphased=$(bcftools view -H -i 'GT="0/1"' "$vcf" 2>/dev/null | wc -l)
                n_hom_alt=$(bcftools view -H -i 'GT="1|1" || GT="1/1"' "$vcf" 2>/dev/null | wc -l)
                n_hom_ref=$(bcftools view -H -i 'GT="0|0" || GT="0/0"' "$vcf" 2>/dev/null | wc -l)
                n_total=$(bcftools view -H "$vcf" 2>/dev/null | wc -l)
                printf "%s\t%d\t%d\t%d\t%d\t%d\n" "$chr" "$n_total" "$n_het_phased" "$n_het_unphased" "$n_hom_alt" "$n_hom_ref" >> "{output.snp_stats}"
            done
            printf "%s\n" {input.vcf_files} > "{output.lst_file}"
            bcftools concat -f "{output.lst_file}" -Ou \
            | bcftools view -Oz -m2 -M2 -i 'GT="0|1" || GT="1|0"' \
                -o "{output.phased_vcf}" 2> "{log}"
            tabix -f -p vcf "{output.phased_vcf}" 2>> "{log}"
            """

    rule parse_genetic_map:
        input:
            gmap_files=lambda wc: [get_genetic_map(c) for c in config["chromosomes"]],
        output:
            gmap_tsv=config["phase_dir"] + "/genetic_map.tsv.gz",
        log:
            config["log_dir"] + f"/parse_genetic_map/{_run_id}.log",
        benchmark:
            config["bench_dir"] + f"/parse_genetic_map/{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        threads: 1
        params:
            chrnames=config["chromosomes"],
            phaser=config["phaser"],
            reference_version=config["reference_version"],
        script:
            "../scripts/parse_genetic_map.py"
