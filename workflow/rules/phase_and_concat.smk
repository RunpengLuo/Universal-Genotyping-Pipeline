"""Merge the per-replicate pileups onto one shared SNP set, and phase them.

Last update: 2026-08-11

Rules:
- [bulk] phase_and_concat_bulk: one joint matrix over every bulk replicate
- [non-bulk] phase_and_concat_nonbulk: one joint matrix over every cell
Outputs:
- allele_dir/snps.tsv.gz: the run's SNP grid, the matrix rows
- allele_dir/snp.{T,A,B}allele.npz: one phased allele matrix
- allele_dir/sample_ids.tsv: the observation roster
- allele_dir/barcodes.tsv.gz: single-cell only, the matrix column axis
"""

if workflow_mode == "bulk_genotyping":

    rule phase_and_concat_bulk:
        input:
            counts=[
                pileup_dir + f"/{at}_{rid}/bcftools.counts.tsv.gz"
                for at in assay_types
                for rid in assay2dataset_ids[at]
            ],
            snp_vcf=phased_snp_vcf,
            region_bed=segment_bed,
            genome_size=genome_size,
            gtf_file=gtf_file,
            blacklist_bed=blacklist_bed,
        output:
            snp_info=allele_dir + "/snps.tsv.gz",
            tot_mtx_snp=allele_dir + "/snp.Tallele.npz",
            a_mtx_snp=allele_dir + "/snp.Aallele.npz",
            b_mtx_snp=allele_dir + "/snp.Ballele.npz",
            sample_file=allele_dir + "/sample_ids.tsv",
            qc_pdf=report(
                qc_dir + "/phase_and_concat.bulk.pdf",
                category="QC plots",
                subcategory="phasing / allele freq (bulk)",
            ),
        log:
            log_dir + f"/phase_and_concat.bulk.{_run_id}.log",
        benchmark:
            bench_dir + f"/phase_and_concat.bulk.{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        params:
            qc_dir=qc_dir,
            sample_id=sample_id,
            dataset_assays=[at for at in assay_types for rid in assay2dataset_ids[at]],
            dataset_ids=[rid for at in assay_types for rid in assay2dataset_ids[at]],
            sample_types=[st for at in assay_types for st in assay2sample_types[at]],
            base_dataset_ids=[
                br for at in assay_types for br in assay2base_dataset_ids[at]
            ],
            min_depth=config["params_phase_and_concat"]["min_depth"],
            gamma=config["params_phase_and_concat"]["gamma"],
            exon_only=config["params_phase_and_concat"]["exon_only"],
            run_id=_run_id,
        script:
            """../scripts/phase_and_concat_bulk.py"""

else:

    rule phase_and_concat_nonbulk:
        input:
            vcfs=[
                pileup_dir + f"/{at}_{rid}/cellSNP.base.vcf.gz"
                for at in assay_types
                for rid in assay2dataset_ids[at]
            ],
            sample_tsvs=[
                pileup_dir + f"/{at}_{rid}/cellSNP.samples.tsv"
                for at in assay_types
                for rid in assay2dataset_ids[at]
            ],
            tot_mtxs=[
                pileup_dir + f"/{at}_{rid}/cellSNP.tag.DP.mtx"
                for at in assay_types
                for rid in assay2dataset_ids[at]
            ],
            ad_mtxs=[
                pileup_dir + f"/{at}_{rid}/cellSNP.tag.AD.mtx"
                for at in assay_types
                for rid in assay2dataset_ids[at]
            ],
            snp_vcf=phased_snp_vcf,
            h5ad_files=[
                bb_dir + f"/{at}.h5ad"
                for at in assay_types
                if ASSAY_TYPE2MODALITY[at] == "RNA"
            ],
            region_bed=segment_bed,
            genome_size=genome_size,
            gtf_file=gtf_file,
            blacklist_bed=blacklist_bed,
        output:
            snp_info=allele_dir + "/snps.tsv.gz",
            all_barcodes=allele_dir + "/barcodes.tsv.gz",
            tot_mtx_snp=allele_dir + "/snp.Tallele.npz",
            a_mtx_snp=allele_dir + "/snp.Aallele.npz",
            b_mtx_snp=allele_dir + "/snp.Ballele.npz",
            sample_file=allele_dir + "/sample_ids.tsv",
            qc_pdf=report(
                [qc_dir + f"/phase_and_concat.{at}.pdf" for at in assay_types],
                category="QC plots",
                subcategory="phasing / allele freq",
            ),
        log:
            log_dir + f"/phase_and_concat.nonbulk.{_run_id}.log",
        benchmark:
            bench_dir + f"/phase_and_concat.nonbulk.{_run_id}.tsv"
        conda:
            "../envs/base.yaml"
        params:
            qc_dir=qc_dir,
            sample_id=sample_id,
            assay_types=assay_types,
            dataset_assays=[at for at in assay_types for rid in assay2dataset_ids[at]],
            dataset_ids=[rid for at in assay_types for rid in assay2dataset_ids[at]],
            sample_types=[st for at in assay_types for st in assay2sample_types[at]],
            exon_only=config["params_phase_and_concat"]["exon_only"],
            run_id=_run_id,
        script:
            """../scripts/phase_and_concat_nonbulk.py"""
