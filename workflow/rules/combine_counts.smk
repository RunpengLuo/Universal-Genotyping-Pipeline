"""Bin the SNPs and depth into bbs, the pipeline's output unit.

Last update: 2026-08-28

Rules:
- [bulk] combine_counts: clonal-LOH detection
  (params_combine_counts.detect_loh_tumor_cell_line), then adaptive binning, depth, RDR
  and read-start counts
- [single-cell] combine_counts_nonbulk: adaptive binning over every assay
- [copytyping] combine_counts_fixed_bins{,_rna}: counts onto pre-computed bbs
Outputs:
- bb_dir/MSR{msr}/bulk/: the bulk bbs, one subdir per min_snp_reads
- bb_dir/MSR{msr}/{assay}/: the single-cell bbs, sliced per assay
- bb_dir/{assay}/: the copytyping bbs, no binning so no MSR level
- bb_dir/unit/{bulk,assay}/: the un-binned SNP, window and gene levels binning consumes
- qc_dir/detect_loh.pdf: the het-density decode; only a detect_loh_tumor_cell_line
  run declares it. The regions are the is_loh column of the window and bb tables
"""

if workflow_mode == "bulk_genotyping":

    _detect_loh_cl = config["params_combine_counts"]["detect_loh_tumor_cell_line"]
    # no other run has clonal LOH to find, so no other run declares this
    _loh_out = (
        {
            "loh_pdf": report(
                qc_dir + "/detect_loh.pdf",
                category="QC plots",
                subcategory="bulk binning",
                labels={"plot": "clonal-LOH density"},
            ),
        }
        if _detect_loh_cl
        else {}
    )

    rule combine_counts:
        input:
            dp_corrected=pileup_dir + "/bulk/window.dp.npz",
            rdcount_files=[
                pileup_dir + f"/{at}/{rid}.rdcount.bed.gz"
                for at in assay_types
                for rid in assay2dataset_ids[at]
            ],
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
            **_loh_out,
            bb_file=expand(bb_dir + f"/MSR{{msr}}/bulk/bb.tsv.gz", msr=msr_list),
            flip_tsv_bb=temp(
                expand(
                    bb_dir + f"/MSR{{msr}}/bulk/snp.flip.tsv",
                    msr=msr_list,
                )
            ),
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
            rdcount_mtx_bb=expand(
                bb_dir + f"/MSR{{msr}}/bulk/bb.rdcount.npz",
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
            unit_rdcount_mtx=bb_dir + "/unit/bulk/window.rdcount.npz",
            unit_sample_file=bb_dir + "/unit/bulk/sample_ids.tsv",
            multi_bb_file=bb_dir + "/multi_snp/bulk/bb.tsv.gz",
            multi_tot_mtx=bb_dir + "/multi_snp/bulk/bb.Tallele.npz",
            multi_a_mtx=bb_dir + "/multi_snp/bulk/bb.Aallele.npz",
            multi_b_mtx=bb_dir + "/multi_snp/bulk/bb.Ballele.npz",
            multi_dp_mtx=bb_dir + "/multi_snp/bulk/bb.depth.npz",
            multi_rdr_mtx=bb_dir + "/multi_snp/bulk/bb.rdr.npz",
            multi_rdcount_mtx=bb_dir + "/multi_snp/bulk/bb.rdcount.npz",
            multi_sample_file=bb_dir + "/multi_snp/bulk/sample_ids.tsv",
            qc_stats_pdf=report(
                expand(
                    qc_dir + f"/combine_counts.stats.bulk.MSR{{msr}}.pdf",
                    msr=msr_list,
                ),
                category="QC plots",
                subcategory="bulk binning",
                labels={"plot": "bb statistics"},
            ),
            qc_1d2d_pdf=report(
                expand(
                    qc_dir + f"/combine_counts.1d2d.bulk.MSR{{msr}}.pdf",
                    msr=msr_list,
                ),
                category="QC plots",
                subcategory="bulk binning",
                labels={"plot": "genome-wide RDR/BAF"},
            ),
        log:
            log_dir + f"/combine_counts.bulk.{_run_id}.log",
        benchmark:
            bench_dir + f"/combine_counts.bulk.{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        threads: config["threads"]["combine_counts"]
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
            phaseset_aware_binning=config["params_combine_counts"][
                "phaseset_aware_binning"
            ],
            min_total_reads=config["params_combine_counts"]["min_total_reads"],
            detect_loh_tumor_cell_line=_detect_loh_cl,
            loh_tile_size=config["params_combine_counts"]["loh_tile_size"],
            loh_rate_ratio=config["params_combine_counts"]["loh_rate_ratio"],
            loh_breakpoint_rate=config["params_combine_counts"]["loh_breakpoint_rate"],
            phase_em=config["params_combine_counts"]["phase_em"],
            phase_em_tau=config["params_combine_counts"]["phase_em_tau"],
            phase_em_min_llr=config["params_combine_counts"]["phase_em_min_llr"],
            phase_em_n_grid=config["params_combine_counts"]["phase_em_n_grid"],
            phase_em_grid_eps=config["params_combine_counts"]["phase_em_grid_eps"],
            phase_em_n_restarts=config["params_combine_counts"]["phase_em_n_restarts"],
            phase_em_max_iter=config["params_combine_counts"]["phase_em_max_iter"],
            phase_em_tol=config["params_combine_counts"]["phase_em_tol"],
            phase_em_min_snps=config["params_combine_counts"]["phase_em_min_snps"],
            run_id=_run_id,
        script:
            """../scripts/combine_counts.py"""

    rule apply_snp_em_phase:
        """Re-phase the het SNP VCF with one bb level's within-bb EM orientation.

        `phase_em` decides each SNP's orientation per bb, so the corrected phasing is a
        property of the (SNP grid x bb grid) pair and differs between `min_snp_reads`
        levels. Writing one VCF per level lets a downstream mode that consumes a phased
        VCF - `copytyping_preprocess`, via `het_snp_vcf` - read the same orientation the
        bbs were summed in, by being handed the VCF and the `bb.tsv.gz` of one level.

        The GT swap is a text rewrite, so it stays in the shell rather than reaching for
        `bcftools annotate`, which cannot populate FORMAT/GT from anything but a VCF.
        That is also why `combine_counts` hands the sites over as a file: this rule cannot
        read the flip vector out of the EM. The file is `temp()`, so it is removed once
        this rule has consumed it and does not outlive the run.
        """
        input:
            vcf=phased_snp_vcf,
            flip_sites=bb_dir + "/MSR{msr}/bulk/snp.flip.tsv",
        output:
            vcf=bb_dir + "/MSR{msr}/bulk/phased_het_snps.phase_em.vcf.gz",
            vcf_tbi=bb_dir + "/MSR{msr}/bulk/phased_het_snps.phase_em.vcf.gz.tbi",
        log:
            log_dir + f"/apply_snp_em_phase.MSR{{msr}}.{_run_id}.log",
        benchmark:
            bench_dir + f"/apply_snp_em_phase.MSR{{msr}}.{_run_id}.tsv"
        conda:
            "../envs/bcftools.yaml"
        threads: 1
        shell:
            r"""
            exec > >(tee -a "{log}") 2>&1
            n_flip=$(wc -l < "{input.flip_sites}")
            echo "re-phasing $(basename {input.vcf}) at $n_flip sites -> {output.vcf}"
            bcftools view "{input.vcf}" \
              | awk -F'\t' -v OFS='\t' '
                  NR==FNR {{ f[$1 FS $2]=1; next }}
                  /^#/    {{ print; next }}
                  {{ if (($1 FS $2) in f) {{
                         n=split($10,a,":")
                         if (a[1]=="1|0") a[1]="0|1"; else if (a[1]=="0|1") a[1]="1|0"
                         s=a[1]; for(i=2;i<=n;i++) s=s":"a[i]; $10=s; c++
                     }}
                     print }}
                  END {{ printf("swapped %d genotypes\n", c) > "/dev/stderr" }}
                ' "{input.flip_sites}" - \
              | bgzip -c > "{output.vcf}"
            tabix -f -p vcf "{output.vcf}"
            """

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
            unit_sample_file=[
                bb_dir + f"/unit/{at}/sample_ids.tsv" for at in assay_types
            ],
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
            unit_stats_pdf=report(
                [qc_dir + f"/combine_counts.stats.{at}.pdf" for at in assay_types],
                category="QC plots",
                subcategory="single-cell binning",
                labels={"plot": "unit-level counts"},
            ),
            unit_stats_tsv=report(
                [qc_dir + f"/combine_counts.stats.{at}.tsv" for at in assay_types],
                category="QC stats",
                subcategory="single-cell binning",
                labels={"table": "unit-level counts"},
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
            phaseset_aware_binning=config["params_combine_counts"][
                "phaseset_aware_binning"
            ],
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
            window_bed=window_bed,
            bb_file=bb_file,
        output:
            bb_file=bb_dir + "/{assay_type}/bb.tsv.gz",
            x_count=bb_dir + "/{assay_type}/bb.Xcount.npz",
            tot_mtx_bb=bb_dir + "/{assay_type}/bb.Tallele.npz",
            a_mtx_bb=bb_dir + "/{assay_type}/bb.Aallele.npz",
            b_mtx_bb=bb_dir + "/{assay_type}/bb.Ballele.npz",
            barcodes_out=bb_dir + "/{assay_type}/barcodes.tsv.gz",
            sample_file=bb_dir + "/{assay_type}/sample_ids.tsv",
            unit_snp_file=bb_dir + "/unit/{assay_type}/snp.tsv.gz",
            unit_tot_mtx=bb_dir + "/unit/{assay_type}/snp.Tallele.npz",
            unit_a_mtx=bb_dir + "/unit/{assay_type}/snp.Aallele.npz",
            unit_b_mtx=bb_dir + "/unit/{assay_type}/snp.Ballele.npz",
            unit_barcodes=bb_dir + "/unit/{assay_type}/barcodes.tsv.gz",
            unit_sample_file=bb_dir + "/unit/{assay_type}/sample_ids.tsv",
            unit_window_file=bb_dir + "/unit/{assay_type}/window.tsv.gz",
            unit_window_x=bb_dir + "/unit/{assay_type}/window.Xcount.npz",
            qc_pdf=report(
                qc_dir + "/combine_counts.{assay_type}.pdf",
                category="QC plots",
                subcategory="fixed-bin aggregation",
                labels={"assay": "{assay_type}"},
            ),
            unit_stats_pdf=report(
                qc_dir + "/combine_counts.stats.{assay_type}.pdf",
                category="QC plots",
                subcategory="fixed-bin aggregation",
                labels={"assay": "{assay_type}", "plot": "unit-level counts"},
            ),
            unit_stats_tsv=report(
                qc_dir + "/combine_counts.stats.{assay_type}.tsv",
                category="QC stats",
                subcategory="fixed-bin aggregation",
                labels={"assay": "{assay_type}", "table": "unit-level counts"},
            ),
        log:
            log_dir
            + f"/combine_counts_fixed_bins/combine_counts_fixed_bins.{{assay_type}}.{_run_id}.log",
        benchmark:
            bench_dir
            + f"/combine_counts_fixed_bins/combine_counts_fixed_bins.{{assay_type}}.{_run_id}.tsv"
        wildcard_constraints:
            assay_type="scATAC",
        conda:
            "../envs/base.yaml"
        threads: 1
        params:
            qc_dir=qc_dir,
            sample_id=sample_id,
            assay_type=lambda wc: wc.assay_type,
            chroms=chr_chromosomes,
            run_id=_run_id,
        script:
            """../scripts/combine_counts_fixed_bins.py"""

    # NB: same script; a gene is indivisible, so RNA's unit is the gene, not the window
    use rule combine_counts_fixed_bins as combine_counts_fixed_bins_rna with:
        output:
            bb_file=bb_dir + "/{assay_type}/bb.tsv.gz",
            x_count=bb_dir + "/{assay_type}/bb.Xcount.npz",
            tot_mtx_bb=bb_dir + "/{assay_type}/bb.Tallele.npz",
            a_mtx_bb=bb_dir + "/{assay_type}/bb.Aallele.npz",
            b_mtx_bb=bb_dir + "/{assay_type}/bb.Ballele.npz",
            barcodes_out=bb_dir + "/{assay_type}/barcodes.tsv.gz",
            sample_file=bb_dir + "/{assay_type}/sample_ids.tsv",
            unit_snp_file=bb_dir + "/unit/{assay_type}/snp.tsv.gz",
            unit_tot_mtx=bb_dir + "/unit/{assay_type}/snp.Tallele.npz",
            unit_a_mtx=bb_dir + "/unit/{assay_type}/snp.Aallele.npz",
            unit_b_mtx=bb_dir + "/unit/{assay_type}/snp.Ballele.npz",
            unit_barcodes=bb_dir + "/unit/{assay_type}/barcodes.tsv.gz",
            unit_sample_file=bb_dir + "/unit/{assay_type}/sample_ids.tsv",
            unit_gene_file=bb_dir + "/unit/{assay_type}/gene.tsv.gz",
            unit_gene_x=bb_dir + "/unit/{assay_type}/gene.Xcount.npz",
            qc_pdf=report(
                qc_dir + "/combine_counts.{assay_type}.pdf",
                category="QC plots",
                subcategory="fixed-bin aggregation",
                labels={"assay": "{assay_type}"},
            ),
            unit_stats_pdf=report(
                qc_dir + "/combine_counts.stats.{assay_type}.pdf",
                category="QC plots",
                subcategory="fixed-bin aggregation",
                labels={"assay": "{assay_type}", "plot": "unit-level counts"},
            ),
            unit_stats_tsv=report(
                qc_dir + "/combine_counts.stats.{assay_type}.tsv",
                category="QC stats",
                subcategory="fixed-bin aggregation",
                labels={"assay": "{assay_type}", "table": "unit-level counts"},
            ),
        wildcard_constraints:
            assay_type="(scRNA|VISIUM|VISIUM3prime)",
