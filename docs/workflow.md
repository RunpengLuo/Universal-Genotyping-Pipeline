# Workflow

Three modes, set by `workflow_mode`. Each runs a subset of the stages below and writes its final bins under `bb_dir/`.

Config keys and every output file: [reference.md](reference.md). Sample-file schema: [sample_sheet.md](sample_sheet.md).

`{chrname}`, `{assay_type}`, `{dataset_id}` and `{modality}` are wildcards — one job each. `{msr}` is not: one binning job writes every `MSR{msr}/` subdirectory of the sweep. Bulk allele/bb outputs are one joint set under `bulk/` (WGS and WES mixed).

In **any** mode, setting `het_snp_vcf` skips genotyping, and phasing too unless `het_snp_vcf_phased: false`.

---

## `bulk_genotyping`

Assays: `bulkWGS`, `bulkWGS-lr`, `bulkWES` (WGS and WES may be mixed). Needs bulk records (normal + tumor); per-stream window BEDs are built automatically off `build_segment_bed`.

| Step | Rule | Output |
|------|------|--------|
| Genotype SNPs | `genotype_snps_bulk` | `snp_dir/chr{chrname}.vcf.gz` |
| Phase SNPs | `phase_snps_{eagle,shapeit,longphase}` | `phase_dir/chr{chrname}.vcf.gz` |
| Concat phased VCFs | `concat_and_extract_phased_het_snps` | `phase_dir/phased_het_snps.vcf.gz` |
| Parse genetic map | `parse_genetic_map` | `phase_dir/genetic_map.tsv.gz` |
| Pileup at het SNPs | `pileup_snps_bulk_mode1b` | `pileup_dir/{assay_type}_{dataset_id}/` |
| Phase and concat (joint, all bulk assays) | `phase_and_concat_bulk` | `allele_dir/bulk/` |
| Compute read depth | `run_mosdepth` | `pileup_dir/{assay_type}/out_mosdepth/{dataset_id}.regions.bed.gz` |
| Bias correction | `rd_correct` | `pileup_dir/{assay_type}/window.{dp.npz,tsv.gz}` |
| Adaptive binning + RDR | `combine_counts` | `bb_dir/MSR{msr}/bulk/` |

Every bulk replicate is piled up against one shared phased het-SNP VCF, so `phase_and_concat_bulk` runs **once** and builds a single joint allele matrix — one pseudobulk column per replicate, all bulk assays segmented together. Depth and bias correction stay per-assay. RDR divides each tumor by `params_combine_counts.rdr_normalization` (`auto` \| `median` \| `normal`).

---

## `single_cell_genotyping`

Assays: `scRNA`, `scATAC`, `VISIUM`, `VISIUM3prime`. A multiome pair is one `dataset_id` shared by an `scRNA` and an `scATAC` record.

| Step | Rule | Output |
|------|------|--------|
| Pseudobulk genotyping | `genotype_snps_pseudobulk_mode1b` | `snp_dir/pseudobulk_{modality}/` |
| Annotate SNPs | `annotate_snps_pseudobulk` | `snp_dir/chr{chrname}.vcf.gz` |
| Phase SNPs | `phase_snps_{eagle,shapeit,longphase}` | `phase_dir/chr{chrname}.vcf.gz` |
| Concat phased VCFs | `concat_and_extract_phased_het_snps` | `phase_dir/phased_het_snps.vcf.gz` |
| Parse genetic map | `parse_genetic_map` | `phase_dir/genetic_map.tsv.gz` |
| Single-cell pileup | `pileup_snps_nonbulk_mode1a` | `pileup_dir/{assay_type}_{dataset_id}/` |
| Build RNA AnnData | `process_rna_anndata` | `bb_dir/{assay_type}.h5ad` |
| Phase and concat | `phase_and_concat_nonbulk` | `allele_dir/{assay_type}/` |
| Adaptive binning + fragment/UMI counting | `combine_counts_nonbulk` | `bb_dir/MSR{msr}/{assay_type}/` |

All of the sample's non-bulk assays are segmented on **one shared bin grid** (one pseudobulk column per replicate x assay), duplicated into each per-assay subdirectory. `bb.Xcount.npz` counts native signal per bin: scATAC from the raw fragments, RNA from the h5ad.

To reuse a bulk run's SNPs, point `het_snp_vcf` at its `phase/phased_het_snps.vcf.gz`.

---

## `copytyping_preprocess`

Assays: `scRNA`, `scATAC`, `VISIUM`, `VISIUM3prime`. Never genotypes or phases: `het_snp_vcf` (already phased) and `bb_file` are **required**.

| Step | Rule | Output |
|------|------|--------|
| Single-cell pileup | `pileup_snps_nonbulk_mode1a` | `pileup_dir/{assay_type}_{dataset_id}/` |
| Build AnnData (RNA only) | `process_rna_anndata` | `bb_dir/{assay_type}.h5ad` |
| Phase and concat | `phase_and_concat_nonbulk` | `allele_dir/{assay_type}/` |
| Aggregate onto fixed bins | `combine_counts_fixed_bins` | `bb_dir/{assay_type}/` |

No binning rule runs here: `combine_counts_fixed_bins` aggregates `allele_dir` onto the given `bb_file` blocks, so the output is flat — there is no `MSR{msr}/` layer.
