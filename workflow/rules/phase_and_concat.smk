##################################################


# all bulk assays jointly processed in one job -> one joint matrix
rule phase_and_concat_bulk:
    input:
        counts=[
            config["pileup_dir"] + f"/{at}_{rid}/bcftools.counts.tsv.gz"
            for at in assay_types
            for rid in assay2dataset_ids[at]
        ],
        snp_vcf=phased_snp_vcf,
        region_bed=segment_bed,
        genome_size=config["genome_size"],
        gtf_file=config["gtf_file"],
        blacklist_bed=config["blacklist_bed"] or [],
    output:
        snp_info=config["allele_dir"] + "/bulk/snps.tsv.gz",
        tot_mtx_snp=config["allele_dir"] + "/bulk/snp.Tallele.npz",
        a_mtx_snp=config["allele_dir"] + "/bulk/snp.Aallele.npz",
        b_mtx_snp=config["allele_dir"] + "/bulk/snp.Ballele.npz",
        sample_file=config["allele_dir"] + "/bulk/sample_ids.tsv",
        qc_pdf=report(
            config["qc_dir"] + "/phase_and_concat.bulk.pdf",
            category="QC plots",
            subcategory="phasing / allele freq (bulk)",
        ),
    log:
        config["log_dir"] + f"/phase_and_concat/phase_and_concat.bulk.{_run_id}.log",
    benchmark:
        config["bench_dir"] + f"/phase_and_concat/phase_and_concat.bulk.{_run_id}.tsv"
    conda:
        "../envs/base.yaml"
    params:
        qc_dir=config["qc_dir"],
        sample_id=sample_id,
        obs_assays=[at for at in assay_types for rid in assay2dataset_ids[at]],
        obs_reps=[rid for at in assay_types for rid in assay2dataset_ids[at]],
        obs_sample_types=[st for at in assay_types for st in assay2sample_types[at]],
        obs_base_reps=[br for at in assay_types for br in assay2base_reps[at]],
        min_depth=config["params_phase_and_concat"]["min_depth"],
        gamma=config["params_phase_and_concat"]["gamma"],
        exon_only=config["params_phase_and_concat"]["exon_only"],
        run_id=_run_id,
    script:
        """../scripts/phase_and_concat_bulk.py"""


##################################################


rule phase_and_concat_nonbulk:
    input:
        vcfs=lambda wc: [
            config["pileup_dir"] + f"/{wc.assay_type}_{dataset_id}/cellSNP.base.vcf.gz"
            for dataset_id in assay2dataset_ids[wc.assay_type]
        ],
        sample_tsvs=lambda wc: [
            config["pileup_dir"] + f"/{wc.assay_type}_{dataset_id}/cellSNP.samples.tsv"
            for dataset_id in assay2dataset_ids[wc.assay_type]
        ],
        tot_mtxs=lambda wc: [
            config["pileup_dir"] + f"/{wc.assay_type}_{dataset_id}/cellSNP.tag.DP.mtx"
            for dataset_id in assay2dataset_ids[wc.assay_type]
        ],
        ad_mtxs=lambda wc: [
            config["pileup_dir"] + f"/{wc.assay_type}_{dataset_id}/cellSNP.tag.AD.mtx"
            for dataset_id in assay2dataset_ids[wc.assay_type]
        ],
        snp_vcf=phased_snp_vcf,
        h5ad_file=lambda wc: (
            config["bb_dir"] + f"/{wc.assay_type}.h5ad"
            if ASSAY_TYPE2MODALITY[wc.assay_type] == "RNA"
            else []
        ),
        region_bed=lambda wc: config["region_bed"],
        genome_size=lambda wc: config["genome_size"],
        gtf_file=lambda wc: config["gtf_file"],
        blacklist_bed=lambda wc: branch(
            config["blacklist_bed"] is None,
            then=[],
            otherwise=config["blacklist_bed"],
        ),
    output:
        all_barcodes=config["allele_dir"] + "/{assay_type}/barcodes.tsv.gz",
        barcodes_full=config["allele_dir"] + "/{assay_type}/barcodes.full.tsv.gz",
        snp_info=config["allele_dir"] + "/{assay_type}/snps.tsv.gz",
        tot_mtx_snp=config["allele_dir"] + "/{assay_type}/snp.Tallele.npz",
        a_mtx_snp=config["allele_dir"] + "/{assay_type}/snp.Aallele.npz",
        b_mtx_snp=config["allele_dir"] + "/{assay_type}/snp.Ballele.npz",
        sample_file=config["allele_dir"] + "/{assay_type}/sample_ids.tsv",
        unique_snp_ids=config["allele_dir"] + "/{assay_type}/unique_snp_ids.npy",
        qc_pdf=report(
            config["qc_dir"] + "/phase_and_concat.{assay_type}.pdf",
            category="QC plots",
            subcategory="phasing / allele freq",
            labels={"assay": "{assay_type}"},
        ),
    log:
        config["log_dir"]
        + f"/phase_and_concat/phase_and_concat.{{assay_type}}.{_run_id}.log",
    benchmark:
        config["bench_dir"]
        + f"/phase_and_concat/phase_and_concat.{{assay_type}}.{_run_id}.tsv"
    wildcard_constraints:
        assay_type="(scRNA|scATAC|VISIUM|VISIUM3prime)",
    conda:
        "../envs/base.yaml"
    params:
        qc_dir=config["qc_dir"],
        sample_id=sample_id,
        assay_type=lambda wc: wc.assay_type,
        dataset_ids=lambda wc: assay2dataset_ids[wc.assay_type],
        sample_types=lambda wc: assay2sample_types[wc.assay_type],
        exon_only=config["params_phase_and_concat"]["exon_only"],
        run_id=_run_id,
    script:
        """../scripts/phase_and_concat_nonbulk.py"""
