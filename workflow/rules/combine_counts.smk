"""Bin the SNPs and depth into bbs, the pipeline's output unit.

Last update: 2026-08-11

Rules:
- [bulk] combine_counts: adaptive binning, depth aggregation and RDR
- [single-cell] combine_counts_nonbulk: adaptive binning over every assay
- [copytyping] combine_counts_fixed_bins: counts onto pre-computed bbs
Outputs:
- bb_dir/MSR{msr}/bulk/: the bulk bbs, one subdir per min_snp_reads
- bb_dir/MSR{msr}/{assay}/: the single-cell bbs, sliced per assay
- bb_dir/{assay}/: the copytyping bbs, no binning so no MSR level
- bb_dir/unit/{bulk,assay}/: the un-binned SNP, window and gene levels binning consumes
"""

if workflow_mode == "bulk_genotyping":

    rule combine_counts:
        input:
            dp_corrected=pileup_dir + "/bulk/window.dp.npz",
            window_bed=window_bed,
            snp_info=allele_dir + "/snps.tsv.gz",
            tot_mtx_snp=allele_dir + "/snp.Tallele.npz",
            a_mtx_snp=allele_dir + "/snp.Aallele.npz",
            b_mtx_snp=allele_dir + "/snp.Ballele.npz",
            sample_file=allele_dir + "/sample_ids.tsv",
            gmap_file=gmap_file,
            region_bed=segment_bed,
            blacklist_bed=blacklist_bed,
            genome_size=genome_size,
        output:
            bb_file=expand(bb_dir + f"/MSR{{msr}}/bulk/bb.tsv.gz", msr=msr_list),
            tot_mtx_bb=expand(
                bb_dir + f"/MSR{{msr}}/bulk/bb.Tallele.npz",
                msr=msr_list,
            ),
            a_mtx_bb=expand(
                bb_dir + f"/MSR{{msr}}/bulk/bb.Aallele.npz",
                msr=msr_list,
            ),
            b_mtx_bb=expand(
                bb_dir + f"/MSR{{msr}}/bulk/bb.Ballele.npz",
                msr=msr_list,
            ),
            dp_mtx_bb=expand(
                bb_dir + f"/MSR{{msr}}/bulk/bb.depth.npz",
                msr=msr_list,
            ),
            rdr_mtx_bb=expand(
                bb_dir + f"/MSR{{msr}}/bulk/bb.rdr.npz",
                msr=msr_list,
            ),
            sample_file=expand(
                bb_dir + f"/MSR{{msr}}/bulk/sample_ids.tsv",
                msr=msr_list,
            ),
            unit_snp_file=bb_dir + "/unit/bulk/snp.tsv.gz",
            unit_tot_mtx=bb_dir + "/unit/bulk/snp.Tallele.npz",
            unit_a_mtx=bb_dir + "/unit/bulk/snp.Aallele.npz",
            unit_b_mtx=bb_dir + "/unit/bulk/snp.Ballele.npz",
            unit_window_file=bb_dir + "/unit/bulk/window.tsv.gz",
            unit_dp_mtx=bb_dir + "/unit/bulk/window.depth.npz",
            unit_sample_file=bb_dir + "/unit/bulk/sample_ids.tsv",
            multi_bb_file=bb_dir + "/multi_snp/bulk/bb.tsv.gz",
            multi_tot_mtx=bb_dir + "/multi_snp/bulk/bb.Tallele.npz",
            multi_a_mtx=bb_dir + "/multi_snp/bulk/bb.Aallele.npz",
            multi_b_mtx=bb_dir + "/multi_snp/bulk/bb.Ballele.npz",
            multi_dp_mtx=bb_dir + "/multi_snp/bulk/bb.depth.npz",
            multi_rdr_mtx=bb_dir + "/multi_snp/bulk/bb.rdr.npz",
            multi_sample_file=bb_dir + "/multi_snp/bulk/sample_ids.tsv",
            qc_pdf=report(
                expand(
                    qc_dir + f"/combine_counts.bulk.MSR{{msr}}.pdf",
                    msr=msr_list,
                ),
                category="QC plots",
                subcategory="bulk binning",
            ),
        log:
            log_dir + f"/combine_counts.bulk.{_run_id}.log",
        benchmark:
            bench_dir + f"/combine_counts.bulk.{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        threads: 1
        params:
            qc_dir=qc_dir,
            sample_id=sample_id,
            assay_types=assay_types,
            dataset_ids=[rid for at in assay_types for rid in assay2dataset_ids[at]],
            dataset_assays=[at for at in assay_types for rid in assay2dataset_ids[at]],
            chroms=chr_chromosomes,
            nu=config["params_combine_counts"]["nu"],
            min_switchprob=config["params_combine_counts"]["min_switchprob"],
            switchprob_ps=config["params_combine_counts"]["switchprob_ps"],
            min_snp_reads=msr_list,
            min_snp_per_bin=config["params_combine_counts"]["min_snp_per_bin"],
            nsnp_multi=config["params_combine_counts"]["nsnp_multi"],
            gene_aware_binning=config["params_combine_counts"]["gene_aware_binning"],
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
            snp_info=allele_dir + "/snps.tsv.gz",
            tot_mtx_snp=allele_dir + "/snp.Tallele.npz",
            a_mtx_snp=allele_dir + "/snp.Aallele.npz",
            b_mtx_snp=allele_dir + "/snp.Ballele.npz",
            sample_file=allele_dir + "/sample_ids.tsv",
            all_barcodes=allele_dir + "/barcodes.tsv.gz",
            frag_files=file_input(
                [
                    get_data[("scATAC", rid)]["fragments"]
                    for rid in assay2dataset_ids["scATAC"]
                ]
            ),
            h5ad_files=[
                bb_dir + f"/{at}.h5ad"
                for at in assay_types
                if ASSAY_TYPE2MODALITY[at] == "RNA"
            ],
            gmap_file=gmap_file,
            window_bed=window_bed,
            genome_size=genome_size,
        output:
            bb_file=[
                bb_dir + f"/MSR{msr}/{at}/bb.tsv.gz"
                for at in assay_types
                for msr in msr_list
            ],
            sample_file=[
                bb_dir + f"/MSR{msr}/{at}/sample_ids.tsv"
                for at in assay_types
                for msr in msr_list
            ],
            tot_mtx_bb=[
                bb_dir + f"/MSR{msr}/{at}/bb.Tallele.npz"
                for at in assay_types
                for msr in msr_list
            ],
            a_mtx_bb=[
                bb_dir + f"/MSR{msr}/{at}/bb.Aallele.npz"
                for at in assay_types
                for msr in msr_list
            ],
            b_mtx_bb=[
                bb_dir + f"/MSR{msr}/{at}/bb.Ballele.npz"
                for at in assay_types
                for msr in msr_list
            ],
            unit_snp_file=[bb_dir + f"/unit/{at}/snp.tsv.gz" for at in assay_types],
            unit_tot_mtx=[bb_dir + f"/unit/{at}/snp.Tallele.npz" for at in assay_types],
            unit_a_mtx=[bb_dir + f"/unit/{at}/snp.Aallele.npz" for at in assay_types],
            unit_b_mtx=[bb_dir + f"/unit/{at}/snp.Ballele.npz" for at in assay_types],
            unit_barcodes=[bb_dir + f"/unit/{at}/barcodes.tsv.gz" for at in assay_types],
            unit_sample_file=[bb_dir + f"/unit/{at}/sample_ids.tsv" for at in assay_types],
            unit_window_file=[
                bb_dir + f"/unit/{at}/window.tsv.gz"
                for at in assay_types
                if at == "scATAC"
            ],
            unit_window_x=[
                bb_dir + f"/unit/{at}/window.Xcount.npz"
                for at in assay_types
                if at == "scATAC"
            ],
            unit_gene_file=[
                bb_dir + f"/unit/{at}/gene.tsv.gz"
                for at in assay_types
                if ASSAY_TYPE2MODALITY[at] == "RNA"
            ],
            unit_gene_x=[
                bb_dir + f"/unit/{at}/gene.Xcount.npz"
                for at in assay_types
                if ASSAY_TYPE2MODALITY[at] == "RNA"
            ],
            multi_snp_file=[bb_dir + f"/multi_snp/{at}/bb.tsv.gz" for at in assay_types],
            tot_mtx_multi=[
                bb_dir + f"/multi_snp/{at}/bb.Tallele.npz" for at in assay_types
            ],
            a_mtx_multi=[
                bb_dir + f"/multi_snp/{at}/bb.Aallele.npz" for at in assay_types
            ],
            b_mtx_multi=[
                bb_dir + f"/multi_snp/{at}/bb.Ballele.npz" for at in assay_types
            ],
            all_barcodes=[
                bb_dir + f"/MSR{msr}/{at}/barcodes.tsv.gz"
                for at in assay_types
                for msr in msr_list
            ],
            x_count=[
                bb_dir + f"/MSR{msr}/{at}/bb.Xcount.npz"
                for at in assay_types
                for msr in msr_list
            ],
            qc_pdf=report(
                [
                    qc_dir + f"/combine_counts.{at}.MSR{msr}.pdf"
                    for at in assay_types
                    for msr in msr_list
                ],
                category="QC plots",
                subcategory="single-cell binning",
            ),
        log:
            log_dir + f"/combine_counts_nonbulk.{_run_id}.log",
        benchmark:
            bench_dir + f"/combine_counts_nonbulk.{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        threads: 1
        params:
            qc_dir=qc_dir,
            sample_id=sample_id,
            assay_types=assay_types,
            chroms=chr_chromosomes,
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
            snp_info=allele_dir + "/snps.tsv.gz",
            tot_mtx_snp=allele_dir + "/snp.Tallele.npz",
            a_mtx_snp=allele_dir + "/snp.Aallele.npz",
            b_mtx_snp=allele_dir + "/snp.Ballele.npz",
            sample_file=allele_dir + "/sample_ids.tsv",
            all_barcodes=allele_dir + "/barcodes.tsv.gz",
            h5ad_file=lambda wc: (
                bb_dir + f"/{wc.assay_type}.h5ad"
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
            genome_size=genome_size,
            bb_file=bb_file,
        output:
            bb_file=bb_dir + "/{assay_type}/bb.tsv.gz",
            x_count=bb_dir + "/{assay_type}/bb.Xcount.npz",
            tot_mtx_bb=bb_dir + "/{assay_type}/bb.Tallele.npz",
            a_mtx_bb=bb_dir + "/{assay_type}/bb.Aallele.npz",
            b_mtx_bb=bb_dir + "/{assay_type}/bb.Ballele.npz",
            barcodes_out=bb_dir + "/{assay_type}/barcodes.tsv.gz",
            sample_file=bb_dir + "/{assay_type}/sample_ids.tsv",
            qc_pdf=report(
                qc_dir + "/combine_counts_fixed_bins.{assay_type}.pdf",
                category="QC plots",
                subcategory="fixed-bin aggregation",
                labels={"assay": "{assay_type}"},
            ),
        log:
            log_dir
            + f"/combine_counts_fixed_bins/combine_counts_fixed_bins.{{assay_type}}.{_run_id}.log",
        benchmark:
            bench_dir
            + f"/combine_counts_fixed_bins/combine_counts_fixed_bins.{{assay_type}}.{_run_id}.tsv"
        wildcard_constraints:
            assay_type="(scRNA|scATAC|VISIUM|VISIUM3prime)",
        conda:
            "../envs/base.yaml"
        threads: 1
        params:
            qc_dir=qc_dir,
            sample_id=sample_id,
            assay_type=lambda wc: wc.assay_type,
            run_id=_run_id,
        script:
            """../scripts/combine_counts_fixed_bins.py"""
