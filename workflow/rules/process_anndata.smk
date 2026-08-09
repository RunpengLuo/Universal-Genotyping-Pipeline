##################################################
# Build AnnData objects from single-cell data (scRNA / scATAC / VISIUM)
##################################################


rule process_rna_anndata:
    input:
        barcodes=lambda wc: file_input(
            [
                get_data[(wc.assay_type, rid)]["barcodes"]
                for rid in assay2dataset_ids[wc.assay_type]
            ]
        ),
        matrix_h5=lambda wc: file_input(
            [
                get_data[(wc.assay_type, rid)]["matrix_h5"]
                for rid in assay2dataset_ids[wc.assay_type]
            ]
        ),
        # spatial only; staged into a Space Ranger layout for squidpy
        spatial_files=lambda wc: (
            file_input(
                spatial_layout(wc.assay_type, assay2dataset_ids[wc.assay_type])[1]
            )
            if wc.assay_type in SPATIAL_ASSAYS
            else []
        ),
        region_bed=lambda wc: config["region_bed"],
        gene_blacklist_file=lambda wc: branch(
            config["gene_blacklist_file"] is None,
            then=[],
            otherwise=config["gene_blacklist_file"],
        ),
        gtf_file=lambda wc: config["gtf_file"],
    output:
        h5ad_file=config["bb_dir"] + "/{assay_type}.h5ad",
    log:
        config["log_dir"]
        + f"/process_rna_anndata/process_rna_anndata.{{assay_type}}.{_run_id}.log",
    benchmark:
        config["bench_dir"]
        + f"/process_rna_anndata/process_rna_anndata.{{assay_type}}.{_run_id}.tsv"
    wildcard_constraints:
        assay_type="(scRNA|VISIUM|VISIUM3prime)",
    conda:
        "../envs/base.yaml"
    params:
        assay_type=lambda wc: wc.assay_type,
        dataset_ids=lambda wc: assay2dataset_ids[wc.assay_type],
        # per-dataset spatial/ filenames, aligned with input.spatial_files
        spatial_names=lambda wc: (
            spatial_layout(wc.assay_type, assay2dataset_ids[wc.assay_type])[0]
            if wc.assay_type in SPATIAL_ASSAYS
            else []
        ),
        min_frac_barcodes=lambda wc: config["params_process_anndata"][
            "min_frac_barcodes"
        ],
        gene_id_colname=lambda wc: config["params_process_anndata"]["gene_id_colname"],
    script:
        """../scripts/process_rna_anndata.py"""
