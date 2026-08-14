"""Phase the called SNPs, then extract one het-SNP VCF.

Last update: 2026-08-12

Rules:
- [optional] phase_snps_shapeit, phase_snps_eagle: panel phasing, one chromosome each
- [optional] phase_snps_longphase: read-based phasing for long-read bulk
- concat_and_extract_phased_het_snps: concatenate, keep the phased hets
- [optional] parse_genetic_map: the phaser maps into one chr-prefixed table
Outputs:
- phase_dir/phased_het_snps.vcf.gz: the parent SNP set for every pileup
"""

if phaser == "shapeit":

    rule phase_snps_shapeit:
        input:
            snp_vcf=lambda wc: snp_dir + f"/chr{wc.chrname}.vcf.gz",
            phasing_panel_file=lambda wc: get_phasing_panel(wc.chrname),
            gmap_file=lambda wc: get_genetic_map(wc.chrname),
        output:
            phased_file=phase_dir + "/chr{chrname}.vcf.gz",
            bcf_file=temp(phase_dir + "/chr{chrname}.bcf"),
            bcf_file_csi=temp(phase_dir + "/chr{chrname}.bcf.csi"),
        log:
            log_dir
            + f"/phase_snps_shapeit/phase_snps_shapeit.chr{{chrname}}.{_run_id}.log",
        benchmark:
            bench_dir
            + f"/phase_snps_shapeit/phase_snps_shapeit.chr{{chrname}}.{_run_id}.tsv"
        conda:
            "../envs/shapeit.yaml"
        threads: config["threads"]["phase"]
        params:
            chrom=lambda wc: input_chrom(wc.chrname),
        shell:
            r"""
            SHAPEIT5_phase_common \
                --input "{input.snp_vcf}" \
                --map "{input.gmap_file}" \
                --reference "{input.phasing_panel_file}" \
                --region "{params.chrom}" \
                --thread "{threads}" \
                --output "{output.bcf_file}" > {log} 2>&1

            bcftools view --output-type v "{output.bcf_file}" | bgzip > "{output.phased_file}"
            tabix -f -p vcf "{output.phased_file}"
            """


if phaser == "eagle":

    rule phase_snps_eagle:
        input:
            snp_vcf=lambda wc: snp_dir + f"/chr{wc.chrname}.vcf.gz",
            phasing_panel_file=lambda wc: get_phasing_panel(wc.chrname),
            gmap_file=lambda wc: get_genetic_map(wc.chrname),
        output:
            phased_file=phase_dir + "/chr{chrname}.vcf.gz",
        log:
            log_dir + f"/phase_snps_eagle/phase_snps_eagle.chr{{chrname}}.{_run_id}.log",
        benchmark:
            bench_dir
            + f"/phase_snps_eagle/phase_snps_eagle.chr{{chrname}}.{_run_id}.tsv"
        conda:
            "../envs/eagle.yaml"
        threads: config["threads"]["phase"]
        params:
            out_prefix=phase_dir + "/chr{chrname}",
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


if phaser == "longphase":

    rule phase_snps_longphase:
        input:
            snp_vcf=lambda wc: snp_dir + f"/chr{wc.chrname}.vcf.gz",
            alignment=lambda wc: bam_stream_input(phase_files),
            alignment_index=lambda wc: bam_stream_index_input(phase_files),
            reference=reference,
        output:
            phased_file=phase_dir + "/chr{chrname}.vcf.gz",
        log:
            log_dir
            + f"/phase_snps_longphase/phase_snps_longphase.chr{{chrname}}.{_run_id}.log",
        benchmark:
            bench_dir
            + f"/phase_snps_longphase/phase_snps_longphase.chr{{chrname}}.{_run_id}.tsv"
        conda:
            "../envs/longphase.yaml"
        threads: config["threads"]["phase"]
        resources:
            downloads=download_slots(phase_files),
        params:
            min_mapq=config["params_longphase"]["min_mapq"],
            extra_params=config["params_longphase"]["extra_params"],
            out_prefix=phase_dir + "/chr{chrname}",
            bam_args=lambda wc, input: " ".join(
                [f"--bam-file={p}" for p in input.alignment]
                + [f"--bam-file={t}" for t in bam_stream_arg(phase_files).split()]
            ),
        shell:
            r"""
            longphase phase \
                {params.bam_args} \
                --reference={input.reference} \
                --snp-file={input.snp_vcf} \
                --mappingQuality={params.min_mapq} \
                --out-prefix={params.out_prefix} \
                --threads={threads} \
                {params.extra_params} > {log} 2>&1
            bgzip -f "{params.out_prefix}.vcf"
            tabix -f -p vcf "{output.phased_file}"
            """


rule concat_and_extract_phased_het_snps:
    input:
        vcf_files=expand(
            phase_dir + "/chr{chrname}.vcf.gz",
            chrname=nochr_chromosomes,
        ),
    output:
        phased_vcf=phase_dir + "/phased_het_snps.vcf.gz",
        phased_vcf_tbi=phase_dir + "/phased_het_snps.vcf.gz.tbi",
        snp_stats=report(
            phase_dir + "/germline_snp_statistics.tsv",
            category="QC stats",
            subcategory="phasing",
            labels={"table": "germline SNP statistics"},
        ),
        lst_file=temp(phase_dir + "/phased_snps.lst"),
    log:
        log_dir
        + f"/concat_and_extract_phased_het_snps/concat_and_extract_phased_het_snps.{_run_id}.log",
    benchmark:
        bench_dir
        + f"/concat_and_extract_phased_het_snps/concat_and_extract_phased_het_snps.{_run_id}.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    shell:
        r"""
        printf "#CHR\ttotal\thet_phased\thet_unphased\thom_alt\thom_ref\n" > "{output.snp_stats}"
        for vcf in {input.vcf_files}; do
            chr=$(basename "$vcf" .vcf.gz)
            n_het_phased=$(bcftools view --no-header --include 'GT="0|1" || GT="1|0"' "$vcf" 2>/dev/null | wc -l)
            n_het_unphased=$(bcftools view --no-header --include 'GT="0/1"' "$vcf" 2>/dev/null | wc -l)
            n_hom_alt=$(bcftools view --no-header --include 'GT="1|1" || GT="1/1"' "$vcf" 2>/dev/null | wc -l)
            n_hom_ref=$(bcftools view --no-header --include 'GT="0|0" || GT="0/0"' "$vcf" 2>/dev/null | wc -l)
            n_total=$(bcftools view --no-header "$vcf" 2>/dev/null | wc -l)
            printf "%s\t%d\t%d\t%d\t%d\t%d\n" "$chr" "$n_total" "$n_het_phased" "$n_het_unphased" "$n_hom_alt" "$n_hom_ref" >> "{output.snp_stats}"
        done
        printf "%s\n" {input.vcf_files} > "{output.lst_file}"
        bcftools concat --file-list "{output.lst_file}" --output-type u \
        | bcftools view --output-type z --min-alleles 2 --max-alleles 2 \
            --include 'GT="0|1" || GT="1|0"' \
            --output "{output.phased_vcf}" 2> "{log}"
        tabix -f -p vcf "{output.phased_vcf}" 2>> "{log}"
        """


if gmap_file:

    rule parse_genetic_map:
        input:
            gmap_files=[get_genetic_map(c) for c in nochr_chromosomes],
        output:
            gmap_file=gmap_file,
        log:
            log_dir + f"/parse_genetic_map/parse_genetic_map.{_run_id}.log",
        benchmark:
            bench_dir + f"/parse_genetic_map/parse_genetic_map.{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        threads: 1
        params:
            chroms=nochr_chromosomes,
            phaser=phaser,
            species=species,
        script:
            "../scripts/parse_genetic_map.py"
