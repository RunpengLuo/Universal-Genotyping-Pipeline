"""
Inputs
1. BAM files
2. SNP panels
3. reference genome

Outputs: bi-allelic hom-alt and het ref/alt SNPs, per chromosome.
snps/<chrom>.vcf.gz
"""

if workflow_mode == "bulk_genotyping" and run_genotyping:

    rule genotype_snps_bulk:
        input:
            alignment=alignment_input(genotype_files),
            alignment_index=alignment_index_input(genotype_files),
            target_pos=lambda wc: config["snp_targets"] + "/target.chr{chrname}.pos.gz",
            reference=config["reference"],
        output:
            snp_vcf=config["snp_dir"] + "/chr{chrname}.vcf.gz",
            unfiltered_vcf=temp(config["snp_dir"] + "/chr{chrname}.unfiltered.vcf.gz"),
        log:
            config["log_dir"] + f"/genotype_snps_bulk/chr{{chrname}}.{_run_id}.log",
        benchmark:
            config["bench_dir"] + f"/genotype_snps_bulk/chr{{chrname}}.{_run_id}.tsv"
        conda:
            "../envs/bcftools.yaml"
        threads: config["threads"]["genotype"]
        resources:
            downloads=download_slots(genotype_files),
        params:
            chrom="chr{chrname}",
            min_mapq=config["params_bcftools"]["min_mapq"],
            min_baseq=config["params_bcftools"]["min_baseq"],
            min_dp=config["params_bcftools"]["min_dp"],
            max_depth=config["params_bcftools"]["max_depth"],
            min_qual=config["params_bcftools"]["min_qual"],
        shell:
            r"""
            bcftools mpileup {input.alignment} \
                -f "{input.reference}" \
                -Ou \
                --threads {threads} \
                -a INFO/AD,AD,DP \
                --skip-indels \
                -q {params.min_mapq} \
                -Q {params.min_baseq} \
                -d {params.max_depth} \
                -T {input.target_pos} \
            | bcftools call -m \
                -Oz -o {output.unfiltered_vcf} 2> {log}

            NSAMPLE=$(bcftools query -l {output.unfiltered_vcf} | wc -l | tr -d ' ')
            if [ "$NSAMPLE" -ne 1 ]; then
                echo "ERROR: genotyping produced $NSAMPLE samples; expected 1. Pooled alignments must share one @RG SM tag (config genotype_dataset_ids)." >> {log}
                exit 1
            fi

            TOTAL=$(bcftools view -H {output.unfiltered_vcf} | wc -l | tr -d ' ')

            bcftools view {output.unfiltered_vcf} -v snps -m2 -M2 \
                -i 'QUAL>={params.min_qual} && GT="alt" && FMT/DP>={params.min_dp}' \
                -Oz -o {output.snp_vcf} 2>> {log}

            PASS=$(bcftools view -H {output.snp_vcf} | wc -l | tr -d ' ')
            echo "Total called: $TOTAL, Passed filters: $PASS, Filtered: $((TOTAL - PASS))" >> {log}

            tabix -p vcf {output.snp_vcf}
            """


if workflow_mode == "single_cell_genotyping" and run_genotyping:

    rule genotype_snps_pseudobulk_mode1b:
        input:
            alignments=lambda wc: alignment_input(modality2files[wc.modality]),
            alignment_indexes=lambda wc: alignment_index_input(
                modality2files[wc.modality]
            ),
            snp_panel=config["snp_panel"],
        output:
            out_dir=directory(config["snp_dir"] + "/pseudobulk_{modality}"),
            out_vcf=config["snp_dir"] + "/pseudobulk_{modality}/cellSNP.base.vcf.gz",
            out_tsv=config["snp_dir"] + "/pseudobulk_{modality}/cellSNP.samples.tsv",
            out_dp=config["snp_dir"] + "/pseudobulk_{modality}/cellSNP.tag.DP.mtx",
            out_ad=config["snp_dir"] + "/pseudobulk_{modality}/cellSNP.tag.AD.mtx",
            bam_lst=temp("tmp/bams.{modality}.lst"),
        log:
            config["log_dir"] + f"/genotype_snps_pseudobulk/{{modality}}.{_run_id}.log",
        benchmark:
            config["bench_dir"]
            + f"/genotype_snps_pseudobulk/{{modality}}.{_run_id}.tsv"
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

    rule annotate_snps_pseudobulk:
        input:
            raw_snp_vcfs=[
                config["snp_dir"] + f"/pseudobulk_{modality}/cellSNP.base.vcf.gz"
                for modality in modalities
            ],
            genome_size=config["genome_size"],
        output:
            snp_vcfs=expand(
                config["snp_dir"] + "/chr{chrname}.vcf.gz",
                chrname=config["chromosomes"],
            ),
            snp_vcfs_tbi=expand(
                config["snp_dir"] + "/chr{chrname}.vcf.gz.tbi",
                chrname=config["chromosomes"],
            ),
            snp_stats=report(
                config["snp_dir"] + "/pseudobulk_snp_statistics.tsv",
                category="QC stats",
                subcategory="genotyping",
                labels={"table": "pseudobulk SNP statistics"},
            ),
        log:
            config["log_dir"] + f"/annotate_snps_pseudobulk/{_run_id}.log",
        benchmark:
            config["bench_dir"] + f"/annotate_snps_pseudobulk/{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        threads: 1
        params:
            modalities=modalities,
            min_het_reads=config["params_annotate_snps"]["min_het_reads"],
            min_hom_dp=config["params_annotate_snps"]["min_hom_dp"],
            min_vaf_thres=config["params_annotate_snps"]["min_vaf_thres"],
            filter_nz_OTH=config["params_annotate_snps"]["filter_nz_OTH"],
            filter_hom_ALT=config["params_annotate_snps"]["filter_hom_ALT"],
        script:
            "../scripts/annotate_snps_pseudobulk.py"


if not run_genotyping and run_phasing:

    rule split_het_snp_vcf:
        """Per-chromosome SNPs from a supplied unphased het_snp_vcf, for the phaser."""
        input:
            het_snp_vcf=config["het_snp_vcf"],
        output:
            snp_vcf=config["snp_dir"] + "/chr{chrname}.vcf.gz",
            snp_vcf_tbi=config["snp_dir"] + "/chr{chrname}.vcf.gz.tbi",
        log:
            config["log_dir"] + f"/split_het_snp_vcf/chr{{chrname}}.{_run_id}.log",
        benchmark:
            config["bench_dir"] + f"/split_het_snp_vcf/chr{{chrname}}.{_run_id}.tsv"
        conda:
            "../envs/bcftools.yaml"
        threads: 1
        params:
            chrom="chr{chrname}",
        shell:
            r"""
            if [ ! -f "{input.het_snp_vcf}.tbi" ] && [ ! -f "{input.het_snp_vcf}.csi" ]; then
                tabix -f -p vcf "{input.het_snp_vcf}" 2> {log}
            fi
            bcftools view "{input.het_snp_vcf}" -r "{params.chrom}" \
                -Oz -o "{output.snp_vcf}" 2>> {log}
            tabix -f -p vcf "{output.snp_vcf}" 2>> {log}
            """
