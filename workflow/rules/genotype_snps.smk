"""
Inputs
1. BAM files
2. SNP panels
3. reference genome

Outputs: bi-allelic hom-alt and het ref/alt SNPs, per chromosome.
snps/<chrom>.vcf.gz
"""

if workflow_mode == "bulk_genotyping":

    rule genotype_snps_bulk:
        input:
            bams=genotype_bam,
            target_pos=lambda wc: config["snp_targets"] + "/target.chr{chrname}.pos.gz",
            reference=config["reference"],
        output:
            snp_vcf=config["snp_dir"] + "/chr{chrname}.vcf.gz",
            unfiltered_vcf=temp(config["snp_dir"] + "/chr{chrname}.unfiltered.vcf.gz"),
        log:
            config["log_dir"] + f"/genotype_snps_bulk/chr{{chrname}}.{_run_id}.log",
        threads: config["threads"]["genotype"]
        params:
            chrom="chr{chrname}",
            min_mapq=config["params_bcftools"]["min_mapq"],
            min_baseq=config["params_bcftools"]["min_baseq"],
            min_dp=config["params_bcftools"]["min_dp"],
            max_depth=config["params_bcftools"]["max_depth"],
            min_qual=config["params_bcftools"]["min_qual"],
        conda:
            "../envs/tools.yaml"
        shell:
            r"""
            bcftools mpileup {input.bams} \
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
                echo "ERROR: genotyping BAM produced $NSAMPLE samples; expected 1 (BAM has multiple @RG SM tags?)" >> {log}
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


if workflow_mode == "single_cell_genotyping":

    rule genotype_snps_pseudobulk_mode1b:
        input:
            bams=lambda wc: modality2bams[wc.modality],
            snp_panel=config["snp_panel"],
        output:
            out_dir=directory(config["snp_dir"] + "/pseudobulk_{modality}"),
            out_vcf=config["snp_dir"] + "/pseudobulk_{modality}/cellSNP.base.vcf.gz",
            out_tsv=config["snp_dir"] + "/pseudobulk_{modality}/cellSNP.samples.tsv",
            out_dp=config["snp_dir"] + "/pseudobulk_{modality}/cellSNP.tag.DP.mtx",
            out_ad=config["snp_dir"] + "/pseudobulk_{modality}/cellSNP.tag.AD.mtx",
            bam_lst=temp("tmp/bams.{modality}.lst"),
        log:
            config["log_dir"] + f"/genotype_snps_pseudobulk.{{modality}}.{_run_id}.log",
        threads: config["threads"]["genotype"]
        params:
            UMItag=lambda wc: branch(
                wc.modality == "RNA",
                then=config["params_cellsnp_lite"]["UMItag"],
                otherwise="None",
            ),
            minMAF=config["params_cellsnp_lite"]["minMAF_genotype"],
            minCOUNT=config["params_cellsnp_lite"]["minCOUNT_genotype"],
        conda:
            "../envs/tools.yaml"
        shell:
            r"""
            printf "%s\n" {input.bams} > "{output.bam_lst}"
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
            snp_stats=config["snp_dir"] + "/pseudobulk_snp_statistics.tsv",
        params:
            modalities=modalities,
            min_het_reads=config["params_annotate_snps"]["min_het_reads"],
            min_hom_dp=config["params_annotate_snps"]["min_hom_dp"],
            min_vaf_thres=config["params_annotate_snps"]["min_vaf_thres"],
            filter_nz_OTH=config["params_annotate_snps"]["filter_nz_OTH"],
            filter_hom_ALT=config["params_annotate_snps"]["filter_hom_ALT"],
        threads: 1
        log:
            config["log_dir"] + f"/annotate_snps_pseudobulk.{_run_id}.log",
        conda:
            "../envs/base.yaml"
        script:
            "../scripts/annotate_snps_pseudobulk.py"
