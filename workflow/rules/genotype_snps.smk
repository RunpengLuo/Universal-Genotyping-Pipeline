"""Call het and hom-alt SNPs, or split a given VCF per chromosome.

Last update: 2026-08-11

Rules:
- [bulk] genotype_snps_bulk: bcftools calls one chromosome from the alignments, over
  the snp_panel positions. By default `-T` reads CHROM/POS only, so REF comes from the
  reference and ALT from the reads; `fix_panel_allele` instead constrains the
  call to the panel's REF/ALT, discarding reads that carry any other allele
- [single-cell] genotype_snps_pseudobulk_mode1b: cellsnp-lite calls one modality
- [single-cell] genotype_snps_no_normal: genotype het/hom-alt from the pseudobulk counts,
  there being no matched normal to call against
- [optional] split_het_snp_vcf: split a given het_snp_vcf per chromosome
Outputs:
- snp_dir/chr{chrname}.vcf.gz: bi-allelic het and hom-alt SNPs
"""

if workflow_mode == "bulk_genotyping" and run_genotyping:

    rule genotype_snps_bulk:
        input:
            alignment=bam_stream_input(genotype_files),
            alignment_index=bam_stream_index_input(genotype_files),
            snp_panel=snp_panel,
            reference=reference,
        output:
            snp_vcf=snp_dir + "/chr{chrname}.vcf.gz",
            unfiltered_vcf=temp(snp_dir + "/chr{chrname}.unfiltered.vcf.gz"),
        log:
            log_dir
            + f"/genotype_snps_bulk/genotype_snps_bulk.chr{{chrname}}.{_run_id}.log",
        benchmark:
            bench_dir
            + f"/genotype_snps_bulk/genotype_snps_bulk.chr{{chrname}}.{_run_id}.tsv"
        conda:
            "../envs/bcftools.yaml"
        threads: config["threads"]["genotype"]
        resources:
            downloads=download_slots(genotype_files),
        params:
            min_mapq=config["params_bcftools"]["min_mapq"],
            min_baseq=config["params_bcftools"]["min_baseq"],
            min_dp=config["params_bcftools"]["min_dp"],
            max_depth=config["params_bcftools"]["max_depth"],
            min_qual=config["params_bcftools"]["min_qual"],
            extra_params=config["params_bcftools"]["extra_params"],
            bam_arg=bam_stream_arg(genotype_files),
            chrom=lambda wc: input_chrom(wc.chrname),
            alleles_arg=lambda wc, input: (
                "--constrain alleles --targets-file "
                f"<(bcftools query --regions {input_chrom(wc.chrname)} "
                f"--format '%CHROM\\t%POS\\t%REF,%ALT\\n' {input.snp_panel})"
                if fix_panel_allele
                else ""
            ),
        shell:
            r"""
            ALN="{input.alignment}"; [ -z "$ALN" ] && ALN="{params.bam_arg}"
            bcftools mpileup $ALN \
                --fasta-ref "{input.reference}" \
                --output-type u \
                --annotate INFO/AD,AD,DP \
                --skip-indels \
                --min-MQ {params.min_mapq} \
                --min-BQ {params.min_baseq} \
                --max-depth {params.max_depth} \
                {params.extra_params} \
                --regions {params.chrom} \
                --targets-file "{input.snp_panel}" \
            | bcftools call --multiallelic-caller --variants-only \
                {params.alleles_arg} \
                --threads {threads} \
                --output-type z --output {output.unfiltered_vcf} 2> {log}

            NSAMPLE=$(bcftools query --list-samples {output.unfiltered_vcf} | wc -l | tr -d ' ')
            if [ "$NSAMPLE" -ne 1 ]; then
                echo "ERROR: genotyping produced $NSAMPLE samples; expected 1. Pooled alignments must share one @RG SM tag (config genotype_dataset_ids)." >> {log}
                exit 1
            fi

            TOTAL=$(bcftools query --format '\n' {output.unfiltered_vcf} | wc -l | tr -d ' ')

            bcftools view {output.unfiltered_vcf} \
                --types snps --min-alleles 2 --max-alleles 2 \
                --include 'QUAL>={params.min_qual} && GT="alt" && FMT/DP>={params.min_dp}' \
                --threads {threads} \
                --output-type z --output {output.snp_vcf} 2>> {log}

            PASS=$(bcftools query --format '\n' {output.snp_vcf} | wc -l | tr -d ' ')
            echo "Variant sites called: $TOTAL, Passed filters: $PASS, Filtered: $((TOTAL - PASS))" >> {log}

            tabix -p vcf {output.snp_vcf}
            """


if workflow_mode == "single_cell_genotyping" and run_genotyping:

    rule genotype_snps_pseudobulk_mode1b:
        input:
            alignments=lambda wc: alignment_input(modality2files[wc.modality]),
            alignment_indexes=lambda wc: alignment_index_input(
                modality2files[wc.modality]
            ),
            snp_panel=snp_panel,
        output:
            out_dir=directory(snp_dir + "/pseudobulk_{modality}"),
            out_vcf=snp_dir + "/pseudobulk_{modality}/cellSNP.base.vcf.gz",
            out_tsv=snp_dir + "/pseudobulk_{modality}/cellSNP.samples.tsv",
            out_dp=snp_dir + "/pseudobulk_{modality}/cellSNP.tag.DP.mtx",
            out_ad=snp_dir + "/pseudobulk_{modality}/cellSNP.tag.AD.mtx",
            bam_lst=temp("tmp/bams.{modality}.lst"),
        log:
            log_dir
            + f"/genotype_snps_pseudobulk/genotype_snps_pseudobulk.{{modality}}.{_run_id}.log",
        benchmark:
            bench_dir
            + f"/genotype_snps_pseudobulk/genotype_snps_pseudobulk.{{modality}}.{_run_id}.tsv"
        conda:
            "../envs/cellsnp.yaml"
        threads: config["threads"]["genotype"]
        resources:
            downloads=lambda wc: download_slots(modality2files[wc.modality]),
        params:
            UMItag=lambda wc: branch(
                wc.modality == "RNA",
                then=config["params_cellsnp_lite"]["UMItag"],
                otherwise="None",
            ),
            minMAF=config["params_cellsnp_lite"]["minMAF_genotype"],
            minCOUNT=config["params_cellsnp_lite"]["minCOUNT_genotype"],
        shell:
            r"""
            printf "%s\n" {input.alignments} > "{output.bam_lst}"
            cellsnp-lite \
                -S "{output.bam_lst}" \
                -R "{input.snp_panel}" \
                -O "{output.out_dir}" \
                -p {threads} \
                --minMAF {params.minMAF} \
                --minCOUNT {params.minCOUNT} \
                --UMItag {params.UMItag} \
                --cellTAG None \
                --gzip > {log} 2>&1
            """

    rule genotype_snps_no_normal:
        input:
            raw_snp_vcfs=[
                snp_dir + f"/pseudobulk_{modality}/cellSNP.base.vcf.gz"
                for modality in modalities
            ],
            genome_size=genome_size,
        output:
            snp_vcfs=expand(
                snp_dir + "/chr{chrname}.vcf.gz",
                chrname=nochr_chromosomes,
            ),
            snp_vcfs_tbi=expand(
                snp_dir + "/chr{chrname}.vcf.gz.tbi",
                chrname=nochr_chromosomes,
            ),
        log:
            log_dir + f"/genotype_snps_no_normal.{_run_id}.log",
        benchmark:
            bench_dir + f"/genotype_snps_no_normal.{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        threads: 1
        params:
            chroms=chr_chromosomes,
            input_nochr=input_nochr,
            modalities=modalities,
            min_het_reads=config["params_annotate_snps"]["min_het_reads"],
            min_hom_dp=config["params_annotate_snps"]["min_hom_dp"],
            min_vaf_thres=config["params_annotate_snps"]["min_vaf_thres"],
            filter_nz_OTH=config["params_annotate_snps"]["filter_nz_OTH"],
            filter_hom_ALT=config["params_annotate_snps"]["filter_hom_ALT"],
        script:
            "../scripts/genotype_snps_no_normal.py"


if not run_genotyping and run_phasing:

    rule split_het_snp_vcf:
        """Per-chromosome SNPs from a supplied unphased het_snp_vcf, for the phaser."""
        input:
            het_snp_vcf=het_snp_vcf,
        output:
            snp_vcf=snp_dir + "/chr{chrname}.vcf.gz",
            snp_vcf_tbi=snp_dir + "/chr{chrname}.vcf.gz.tbi",
        log:
            log_dir
            + f"/split_het_snp_vcf/split_het_snp_vcf.chr{{chrname}}.{_run_id}.log",
        benchmark:
            bench_dir
            + f"/split_het_snp_vcf/split_het_snp_vcf.chr{{chrname}}.{_run_id}.tsv"
        conda:
            "../envs/bcftools.yaml"
        threads: 1
        params:
            chrom=lambda wc: input_chrom(wc.chrname),
        shell:
            r"""
            if [ ! -f "{input.het_snp_vcf}.tbi" ] && [ ! -f "{input.het_snp_vcf}.csi" ]; then
                tabix -f -p vcf "{input.het_snp_vcf}" 2> {log}
            fi
            bcftools view "{input.het_snp_vcf}" --regions "{params.chrom}" \
                --output-type z --output "{output.snp_vcf}" 2>> {log}
            tabix -f -p vcf "{output.snp_vcf}" 2>> {log}
            """
