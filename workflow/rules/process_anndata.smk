##################################################
# Build AnnData objects from single-cell data (scRNA / scATAC / VISIUM)
##################################################


rule process_rna_anndata:
    input:
        barcodes=lambda wc: [
            get_data[(wc.assay_type, rid)][0] for rid in assay2rep_ids[wc.assay_type]
        ],
        ranger_dirs=lambda wc: [
            get_data[(wc.assay_type, rid)][2] for rid in assay2rep_ids[wc.assay_type]
        ],
        region_bed=lambda wc: config["region_bed"],
        genome_size=lambda wc: config["genome_size"],
        gene_blacklist_file=lambda wc: branch(
            config["gene_blacklist_file"] is None,
            then=[],
            otherwise=config["gene_blacklist_file"],
        ),
        gtf_file=lambda wc: config["gtf_file"],
    output:
        h5ad_file=config["bb_dir"] + "/{assay_type}/{assay_type}.h5ad",
    params:
        assay_type=lambda wc: wc.assay_type,
        rep_ids=lambda wc: assay2rep_ids[wc.assay_type],
        sample_types=lambda wc: assay2sample_types[wc.assay_type],
        min_frac_barcodes=lambda wc: config["params_process_anndata"][
            "min_frac_barcodes"
        ],
        gene_id_colname=lambda wc: config["params_process_anndata"]["gene_id_colname"],
    wildcard_constraints:
        assay_type="(scRNA|VISIUM|VISIUM3prime)",
    log:
        config["log_dir"] + f"/process_rna_anndata.{{assay_type}}.{_run_id}.log",
    conda:
        "../envs/base.yaml"
    script:
        """../scripts/process_rna_anndata.py"""
