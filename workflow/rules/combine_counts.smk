##################################################
# SNP-informed adaptive binning + window→bin depth aggregation
# Bulk: combine_counts (with corrected window depth)
# Non-bulk: combine_counts_nonbulk (derived bins) + combine_counts_fixed_bins (given bins)
##################################################

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
            snp_info=config["allele_dir"] + "/bulk/snps.tsv.gz",
            tot_mtx_snp=config["allele_dir"] + "/bulk/snp.Tallele.npz",
            a_mtx_snp=config["allele_dir"] + "/bulk/snp.Aallele.npz",
            b_mtx_snp=config["allele_dir"] + "/bulk/snp.Ballele.npz",
            sample_file=config["allele_dir"] + "/bulk/sample_ids.tsv",
            gmap_file=lambda wc: (
                config["phase_dir"] + "/genetic_map.tsv.gz"
                if require_genetic_map
                else []
            ),
            region_bed=segment_bed,
            blacklist_bed=config["blacklist_bed"] or [],
            genome_size=config["genome_size"],
        output:
            bb_file=expand(
                config["bb_dir"] + f"/MSR{{msr}}/bulk/bb.tsv.gz", msr=msr_list
            ),
            tot_mtx_bb=expand(
                config["bb_dir"] + f"/MSR{{msr}}/bulk/bb.Tallele.npz",
                msr=msr_list,
            ),
            a_mtx_bb=expand(
                config["bb_dir"] + f"/MSR{{msr}}/bulk/bb.Aallele.npz",
                msr=msr_list,
            ),
            b_mtx_bb=expand(
                config["bb_dir"] + f"/MSR{{msr}}/bulk/bb.Ballele.npz",
                msr=msr_list,
            ),
            dp_mtx_bb=expand(
                config["bb_dir"] + f"/MSR{{msr}}/bulk/bb.depth.npz",
                msr=msr_list,
            ),
            rdr_mtx_bb=expand(
                config["bb_dir"] + f"/MSR{{msr}}/bulk/bb.rdr.npz",
                msr=msr_list,
            ),
            sample_file=expand(
                config["bb_dir"] + f"/MSR{{msr}}/bulk/sample_ids.tsv",
                msr=msr_list,
            ),
            qc_pdf=report(
                expand(
                    config["qc_dir"] + f"/combine_counts.bulk.MSR{{msr}}.pdf",
                    msr=msr_list,
                ),
                category="QC plots",
                subcategory="bulk binning",
            ),
        log:
            config["log_dir"] + f"/combine_counts/combine_counts.bulk.{_run_id}.log",
        benchmark:
            config["bench_dir"] + f"/combine_counts/combine_counts.bulk.{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        threads: 1
        params:
            qc_dir=config["qc_dir"],
            bulk_assays=assay_types,
            nu=config["params_combine_counts"]["nu"],
            min_switchprob=config["params_combine_counts"]["min_switchprob"],
            switchprob_ps=config["params_combine_counts"]["switchprob_ps"],
            min_snp_reads=msr_list,
            min_snp_per_bin=config["params_combine_counts"]["min_snp_per_bin"],
            gene_aware_binning=config["params_combine_counts"]["gene_aware_binning"],
            rdr_outlier_quantile=config["params_combine_counts"]["rdr_outlier_quantile"],
            max_blocksize=config["params_combine_counts"]["max_blocksize"],
            phase_flip_test=config["params_combine_counts"]["phase_flip_test"],
            phase_flip_epsilon=config["params_combine_counts"]["phase_flip_epsilon"],
            phase_flip_alpha=config["params_combine_counts"]["phase_flip_alpha"],
            run_id=_run_id,
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
            frag_files=file_input(
                [
                    get_data[("scATAC", rid)]["fragments"]
                    for rid in assay2dataset_ids["scATAC"]
                ]
            ),
            h5ad_files=[
                config["bb_dir"] + f"/{at}.h5ad"
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
            bb_file=[
                config["bb_dir"] + f"/MSR{msr}/{at}/bb.tsv.gz"
                for at in assay_types
                for msr in msr_list
            ],
            sample_file=[
                config["bb_dir"] + f"/MSR{msr}/{at}/sample_ids.tsv"
                for at in assay_types
                for msr in msr_list
            ],
            tot_mtx_bb=[
                config["bb_dir"] + f"/MSR{msr}/{at}/bb.Tallele.npz"
                for at in assay_types
                for msr in msr_list
            ],
            a_mtx_bb=[
                config["bb_dir"] + f"/MSR{msr}/{at}/bb.Aallele.npz"
                for at in assay_types
                for msr in msr_list
            ],
            b_mtx_bb=[
                config["bb_dir"] + f"/MSR{msr}/{at}/bb.Ballele.npz"
                for at in assay_types
                for msr in msr_list
            ],
            multi_snp_file=[
                config["bb_dir"] + f"/MSR{msr}/{at}/multi_snp.tsv.gz"
                for at in assay_types
                for msr in msr_list
            ],
            tot_mtx_multi=[
                config["bb_dir"] + f"/MSR{msr}/{at}/multi_snp.Tallele.npz"
                for at in assay_types
                for msr in msr_list
            ],
            a_mtx_multi=[
                config["bb_dir"] + f"/MSR{msr}/{at}/multi_snp.Aallele.npz"
                for at in assay_types
                for msr in msr_list
            ],
            b_mtx_multi=[
                config["bb_dir"] + f"/MSR{msr}/{at}/multi_snp.Ballele.npz"
                for at in assay_types
                for msr in msr_list
            ],
            all_barcodes=[
                config["bb_dir"] + f"/MSR{msr}/{at}/barcodes.tsv.gz"
                for at in assay_types
                for msr in msr_list
            ],
            barcodes_full=[
                config["bb_dir"] + f"/MSR{msr}/{at}/barcodes.full.tsv.gz"
                for at in assay_types
                for msr in msr_list
            ],
            x_count=[
                config["bb_dir"] + f"/MSR{msr}/{at}/bb.Xcount.npz"
                for at in assay_types
                for msr in msr_list
            ],
            qc_pdf=report(
                [
                    config["qc_dir"] + f"/combine_counts.{at}.MSR{msr}.pdf"
                    for at in assay_types
                    for msr in msr_list
                ],
                category="QC plots",
                subcategory="single-cell binning",
            ),
        log:
            config["log_dir"]
            + f"/combine_counts_nonbulk/combine_counts_nonbulk.{_run_id}.log",
        benchmark:
            config["bench_dir"]
            + f"/combine_counts_nonbulk/combine_counts_nonbulk.{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        threads: 1
        params:
            qc_dir=config["qc_dir"],
            nonbulk_assays=assay_types,
            frag_reps=assay2dataset_ids["scATAC"],
            h5ad_assays=[at for at in assay_types if ASSAY_TYPE2MODALITY[at] == "RNA"],
            nu=config["params_combine_counts"]["nu"],
            min_switchprob=config["params_combine_counts"]["min_switchprob"],
            switchprob_ps=config["params_combine_counts"]["switchprob_ps"],
            nsnp_multi=config["params_combine_counts"]["nsnp_multi"],
            min_snp_reads=msr_list,
            min_snp_per_bin=config["params_combine_counts"]["min_snp_per_bin"],
            gene_aware_binning=config["params_combine_counts"]["gene_aware_binning"],
            run_id=_run_id,
        script:
            """../scripts/combine_counts_nonbulk.py"""

elif workflow_mode == "copytyping_preprocess":

    rule combine_counts_fixed_bins:
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
                config["bb_dir"] + f"/{wc.assay_type}.h5ad"
                if ASSAY_TYPE2MODALITY[wc.assay_type] == "RNA"
                else []
            ),
            frag_files=lambda wc: file_input(
                [
                    get_data[(wc.assay_type, rid)]["fragments"]
                    for rid in assay2dataset_ids[wc.assay_type]
                ]
                if wc.assay_type == "scATAC"
                else []
            ),
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
                config["qc_dir"] + "/combine_counts_fixed_bins.{assay_type}.pdf",
                category="QC plots",
                subcategory="fixed-bin aggregation",
                labels={"assay": "{assay_type}"},
            ),
        log:
            config["log_dir"]
            + f"/combine_counts_fixed_bins/combine_counts_fixed_bins.{{assay_type}}.{_run_id}.log",
        benchmark:
            config["bench_dir"]
            + f"/combine_counts_fixed_bins/combine_counts_fixed_bins.{{assay_type}}.{_run_id}.tsv"
        wildcard_constraints:
            assay_type="(scRNA|scATAC|VISIUM|VISIUM3prime)",
        conda:
            "../envs/base.yaml"
        threads: 1
        params:
            qc_dir=config["qc_dir"],
            sample_id=sample_id,
            assay_type=lambda wc: wc.assay_type,
            run_id=_run_id,
        script:
            """../scripts/combine_counts_fixed_bins.py"""
