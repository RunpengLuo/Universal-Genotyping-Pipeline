##################################################
# SNP-informed adaptive binning + window→bin depth aggregation
# Bulk: combine_counts (with corrected window depth)
# Non-bulk: combine_counts_nonbulk + cnv_segmentation
##################################################

# binning-parameter tag for QC plot filenames (distinguishes MSR/MSPB sweeps that
# share one qc_dir)
_cc_tag = (
    f"MSR{config['params_combine_counts']['min_snp_reads']}"
    f"_MSPB{config['params_combine_counts']['min_snp_per_block']}"
)


if workflow_mode == "bulk_genotyping":

    rule combine_counts:
        input:
            # depth/window stay per-assay; allele matrices are one joint set
            dp_corrected=[
                config["pileup_dir"] + f"/{at}/window.dp.npz" for at in assay_types
            ],
            window_df=[
                config["pileup_dir"] + f"/{at}/window.tsv.gz" for at in assay_types
            ],
            snp_info=config["allele_dir"] + f"/{bulk_stream}/snps.tsv.gz",
            tot_mtx_snp=config["allele_dir"] + f"/{bulk_stream}/snp.Tallele.npz",
            a_mtx_snp=config["allele_dir"] + f"/{bulk_stream}/snp.Aallele.npz",
            b_mtx_snp=config["allele_dir"] + f"/{bulk_stream}/snp.Ballele.npz",
            sample_file=config["allele_dir"] + f"/{bulk_stream}/sample_ids.tsv",
            gmap_file=lambda wc: (
                config["phase_dir"] + "/genetic_map.tsv.gz"
                if require_genetic_map
                else []
            ),
            region_bed=config["region_bed"],
            blacklist_bed=config["blacklist_bed"] or [],
            genome_size=config["genome_size"],
        output:
            bb_file=config["bb_dir"] + f"/{bulk_stream}/bb.tsv.gz",
            tot_mtx_bb=config["bb_dir"] + f"/{bulk_stream}/bb.Tallele.npz",
            a_mtx_bb=config["bb_dir"] + f"/{bulk_stream}/bb.Aallele.npz",
            b_mtx_bb=config["bb_dir"] + f"/{bulk_stream}/bb.Ballele.npz",
            dp_mtx_bb=config["bb_dir"] + f"/{bulk_stream}/bb.depth.npz",
            rdr_mtx_bb=config["bb_dir"] + f"/{bulk_stream}/bb.rdr.npz",
            sample_file=config["bb_dir"] + f"/{bulk_stream}/sample_ids.tsv",
            qc_pdf=report(
                config["qc_dir"] + f"/combine_counts.{bulk_stream}.{_cc_tag}.pdf",
                category="QC plots",
                subcategory="bulk binning",
                labels={"stream": bulk_stream, "binning": _cc_tag},
            ),
        params:
            qc_dir=config["qc_dir"],
            bulk_assays=assay_types,
            nu=config["params_combine_counts"]["nu"],
            min_switchprob=config["params_combine_counts"]["min_switchprob"],
            switchprob_ps=config["params_combine_counts"]["switchprob_ps"],
            min_snp_reads=config["params_combine_counts"]["min_snp_reads"],
            min_snp_per_block=config["params_combine_counts"]["min_snp_per_block"],
            gene_aware_binning=config["params_combine_counts"]["gene_aware_binning"],
            median_normalization=config["params_combine_counts"]["median_normalization"],
            rdr_outlier_quantile=config["params_combine_counts"]["rdr_outlier_quantile"],
            max_blocksize=config["params_combine_counts"]["max_blocksize"],
            phase_flip_test=config["params_combine_counts"]["phase_flip_test"],
            phase_flip_epsilon=config["params_combine_counts"]["phase_flip_epsilon"],
            phase_flip_alpha=config["params_combine_counts"]["phase_flip_alpha"],
            run_id=_run_id,
        threads: 1
        log:
            config["log_dir"] + f"/combine_counts.{_run_id}.log",
        conda:
            "../envs/base.yaml"
        script:
            """../scripts/combine_counts.py"""

elif workflow_mode == "single_cell_genotyping":

    rule combine_counts_nonbulk:
        input:
            snp_info=[config["allele_dir"] + f"/{at}/snps.tsv.gz" for at in assay_types],
            tot_mtx_snp=[
                config["allele_dir"] + f"/{at}/snp.Tallele.npz" for at in assay_types
            ],
            a_mtx_snp=[
                config["allele_dir"] + f"/{at}/snp.Aallele.npz" for at in assay_types
            ],
            b_mtx_snp=[
                config["allele_dir"] + f"/{at}/snp.Ballele.npz" for at in assay_types
            ],
            sample_file=[
                config["allele_dir"] + f"/{at}/sample_ids.tsv" for at in assay_types
            ],
            all_barcodes=[
                config["allele_dir"] + f"/{at}/barcodes.tsv.gz" for at in assay_types
            ],
            barcodes_full=[
                config["allele_dir"] + f"/{at}/barcodes.full.tsv.gz"
                for at in assay_types
            ],
            ranger_dirs=[
                get_data[(at, rid)][2]
                for at in assay_types
                for rid in assay2rep_ids[at]
            ],
            h5ad_files=[
                config["bb_dir"] + f"/{at}/{at}.h5ad"
                for at in assay_types
                if ASSAY_TYPE2MODALITY[at] == "RNA"
            ],
            gmap_file=(
                config["phase_dir"] + "/genetic_map.tsv.gz"
                if require_genetic_map
                else []
            ),
            region_bed=config["region_bed"],
            genome_size=config["genome_size"],
        output:
            bb_file=[config["bb_dir"] + f"/{at}/bb.tsv.gz" for at in assay_types],
            sample_file=[
                config["bb_dir"] + f"/{at}/sample_ids.tsv" for at in assay_types
            ],
            tot_mtx_bb=[
                config["bb_dir"] + f"/{at}/bb.Tallele.npz" for at in assay_types
            ],
            a_mtx_bb=[config["bb_dir"] + f"/{at}/bb.Aallele.npz" for at in assay_types],
            b_mtx_bb=[config["bb_dir"] + f"/{at}/bb.Ballele.npz" for at in assay_types],
            multi_snp_file=[
                config["bb_dir"] + f"/{at}/multi_snp.tsv.gz" for at in assay_types
            ],
            tot_mtx_multi=[
                config["bb_dir"] + f"/{at}/multi_snp.Tallele.npz" for at in assay_types
            ],
            a_mtx_multi=[
                config["bb_dir"] + f"/{at}/multi_snp.Aallele.npz" for at in assay_types
            ],
            b_mtx_multi=[
                config["bb_dir"] + f"/{at}/multi_snp.Ballele.npz" for at in assay_types
            ],
            all_barcodes=[
                config["bb_dir"] + f"/{at}/barcodes.tsv.gz" for at in assay_types
            ],
            barcodes_full=[
                config["bb_dir"] + f"/{at}/barcodes.full.tsv.gz" for at in assay_types
            ],
            x_count=[config["bb_dir"] + f"/{at}/bb.Xcount.npz" for at in assay_types],
            qc_pdf=report(
                [
                    config["qc_dir"] + f"/combine_counts.{at}.{_cc_tag}.pdf"
                    for at in assay_types
                ],
                category="QC plots",
                subcategory="single-cell binning",
            ),
        params:
            qc_dir=config["qc_dir"],
            nonbulk_assays=assay_types,
            ranger_assays=[at for at in assay_types for rid in assay2rep_ids[at]],
            ranger_reps=[rid for at in assay_types for rid in assay2rep_ids[at]],
            h5ad_assays=[at for at in assay_types if ASSAY_TYPE2MODALITY[at] == "RNA"],
            nu=config["params_combine_counts"]["nu"],
            min_switchprob=config["params_combine_counts"]["min_switchprob"],
            switchprob_ps=config["params_combine_counts"]["switchprob_ps"],
            nsnp_multi=config["params_combine_counts"]["nsnp_multi"],
            min_snp_reads=config["params_combine_counts"]["min_snp_reads"],
            min_snp_per_block=config["params_combine_counts"]["min_snp_per_block"],
            gene_aware_binning=config["params_combine_counts"]["gene_aware_binning"],
            run_id=_run_id,
        threads: 1
        log:
            config["log_dir"] + f"/combine_counts_nonbulk.{_run_id}.log",
        conda:
            "../envs/base.yaml"
        script:
            """../scripts/combine_counts_nonbulk.py"""

elif workflow_mode == "copytyping_preprocess":

    rule cnv_segmentation:
        input:
            snp_info=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snps.tsv.gz",
            tot_mtx_snp=lambda wc: config["allele_dir"]
            + f"/{wc.assay_type}/snp.Tallele.npz",
            a_mtx_snp=lambda wc: config["allele_dir"]
            + f"/{wc.assay_type}/snp.Aallele.npz",
            b_mtx_snp=lambda wc: config["allele_dir"]
            + f"/{wc.assay_type}/snp.Ballele.npz",
            sample_file=lambda wc: config["allele_dir"]
            + f"/{wc.assay_type}/sample_ids.tsv",
            all_barcodes=config["allele_dir"] + "/{assay_type}/barcodes.tsv.gz",
            barcodes_full=config["allele_dir"] + "/{assay_type}/barcodes.full.tsv.gz",
            h5ad_file=lambda wc: (
                config["bb_dir"] + f"/{wc.assay_type}/{wc.assay_type}.h5ad"
                if ASSAY_TYPE2MODALITY[wc.assay_type] == "RNA"
                else []
            ),
            ranger_dirs=lambda wc: [
                get_data[(wc.assay_type, rid)][2]
                for rid in assay2rep_ids[wc.assay_type]
            ],
            region_bed=lambda wc: config["region_bed"],
            genome_size=lambda wc: config["genome_size"],
            bb_file=lambda wc: config["bb_file"] or [],
        output:
            cnv_segments=config["bb_dir"] + "/{assay_type}/cnv_segments.tsv",
            x_count=config["bb_dir"] + "/{assay_type}/bb.Xcount.npz",
            tot_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.Tallele.npz",
            a_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.Aallele.npz",
            b_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.Ballele.npz",
            barcodes_out=config["bb_dir"] + "/{assay_type}/barcodes.tsv.gz",
            barcodes_full_out=config["bb_dir"] + "/{assay_type}/barcodes.full.tsv.gz",
            sample_file=config["bb_dir"] + "/{assay_type}/sample_ids.tsv",
            qc_pdf=report(
                config["qc_dir"] + "/cnv_segmentation.{assay_type}.pdf",
                category="QC plots",
                subcategory="CNV segmentation",
                labels={"assay": "{assay_type}"},
            ),
        wildcard_constraints:
            assay_type="(scRNA|scATAC|VISIUM|VISIUM3prime)",
        params:
            qc_dir=config["qc_dir"],
            sample_name=SAMPLE_ID,
            assay_type=lambda wc: wc.assay_type,
            feature_type=lambda wc: ASSAY_TYPE2FEATURE_TYPE[wc.assay_type],
            run_id=_run_id,
        threads: 1
        log:
            config["log_dir"] + f"/cnv_segmentation.{{assay_type}}.{_run_id}.log",
        conda:
            "../envs/base.yaml"
        script:
            """../scripts/cnv_segmentation.py"""
