##################################################


# all bulk assays jointly processed in one job -> one joint matrix
rule phase_and_concat_bulk:
    input:
        vcfs=[
            config["pileup_dir"] + f"/{at}_{rid}/cellSNP.base.vcf.gz"
            for at in assay_types
            for rid in assay2rep_ids[at]
        ],
        sample_tsvs=[
            config["pileup_dir"] + f"/{at}_{rid}/cellSNP.samples.tsv"
            for at in assay_types
            for rid in assay2rep_ids[at]
        ],
        tot_mtxs=[
            config["pileup_dir"] + f"/{at}_{rid}/cellSNP.tag.DP.mtx"
            for at in assay_types
            for rid in assay2rep_ids[at]
        ],
        ad_mtxs=[
            config["pileup_dir"] + f"/{at}_{rid}/cellSNP.tag.AD.mtx"
            for at in assay_types
            for rid in assay2rep_ids[at]
        ],
        snp_vcf=(
            config["phase_dir"] + "/phased_het_snps.vcf.gz"
            if run_genotyping
            else config["het_snp_vcf"]
        ),
        region_bed=config["region_bed"],
        genome_size=config["genome_size"],
        gtf_file=config["gtf_file"],
        blacklist_bed=config["blacklist_bed"] or [],
    output:
        snp_info=config["allele_dir"] + f"/{bulk_stream}/snps.tsv.gz",
        tot_mtx_snp=config["allele_dir"] + f"/{bulk_stream}/snp.Tallele.npz",
        a_mtx_snp=config["allele_dir"] + f"/{bulk_stream}/snp.Aallele.npz",
        b_mtx_snp=config["allele_dir"] + f"/{bulk_stream}/snp.Ballele.npz",
        sample_file=config["allele_dir"] + f"/{bulk_stream}/sample_ids.tsv",
    params:
        qc_dir=config["qc_dir"],
        sample_name=SAMPLE_ID,
        col_assays=[at for at in assay_types for rid in assay2rep_ids[at]],
        col_reps=[rid for at in assay_types for rid in assay2rep_ids[at]],
        col_sample_types=[st for at in assay_types for st in assay2sample_types[at]],
        min_depth=config["params_phase_and_concat"]["min_depth"],
        gamma=config["params_phase_and_concat"]["gamma"],
        exon_only=config["params_phase_and_concat"]["exon_only"],
        run_id=_run_id,
    log:
        config["log_dir"] + f"/phase_and_concat.bulk.{_run_id}.log",
    conda:
        "../envs/base.yaml"
    script:
        """../scripts/phase_and_concat_bulk.py"""


##################################################


rule phase_and_concat_single_cell:
    input:
        vcfs=lambda wc: [
            config["pileup_dir"] + f"/{wc.assay_type}_{rep_id}/cellSNP.base.vcf.gz"
            for rep_id in assay2rep_ids[wc.assay_type]
        ],
        sample_tsvs=lambda wc: [
            config["pileup_dir"] + f"/{wc.assay_type}_{rep_id}/cellSNP.samples.tsv"
            for rep_id in assay2rep_ids[wc.assay_type]
        ],
        tot_mtxs=lambda wc: [
            config["pileup_dir"] + f"/{wc.assay_type}_{rep_id}/cellSNP.tag.DP.mtx"
            for rep_id in assay2rep_ids[wc.assay_type]
        ],
        ad_mtxs=lambda wc: [
            config["pileup_dir"] + f"/{wc.assay_type}_{rep_id}/cellSNP.tag.AD.mtx"
            for rep_id in assay2rep_ids[wc.assay_type]
        ],
        snp_vcf=lambda wc: branch(
            run_genotyping,
            then=config["phase_dir"] + "/phased_het_snps.vcf.gz",
            otherwise=config["het_snp_vcf"],
        ),
        h5ad_file=lambda wc: (
            config["bb_dir"] + f"/{wc.assay_type}/{wc.assay_type}.h5ad"
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
        cell_snp_Aallele=config["allele_dir"] + "/{assay_type}/cell_snp_Aallele.npz",
        cell_snp_Ballele=config["allele_dir"] + "/{assay_type}/cell_snp_Ballele.npz",
    wildcard_constraints:
        assay_type="(scRNA|scATAC|VISIUM|VISIUM3prime)",
    params:
        qc_dir=config["qc_dir"],
        sample_name=SAMPLE_ID,
        assay_type=lambda wc: wc.assay_type,
        rep_ids=lambda wc: assay2rep_ids[wc.assay_type],
        sample_types=lambda wc: assay2sample_types[wc.assay_type],
        exon_only=config["params_phase_and_concat"]["exon_only"],
        run_id=_run_id,
    log:
        config["log_dir"] + f"/phase_and_concat.{{assay_type}}.{_run_id}.log",
    conda:
        "../envs/base.yaml"
    script:
        """../scripts/phase_and_concat_single_cell.py"""
