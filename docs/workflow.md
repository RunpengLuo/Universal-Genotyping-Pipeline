# Workflow

Three workflow modes, set via `config["workflow_mode"]`. Each mode runs a subset of the pipeline stages.

---

## `bulk_genotyping`

Assays: `bulkWGS`, `bulkWGS-lr`, `bulkWES`

| Step | Rule | Script / Tool |
|------|------|---------------|
| 1. Genotype SNPs | `genotype_snps_bulk` | bcftools mpileup + call |
| 2. Phase SNPs (per chr) | `phase_snps_{eagle,shapeit,longphase}` | eagle2 / shapeit5 / longphase |
| 3. Concat phased VCFs | `concat_and_extract_phased_het_snps` | bcftools concat + view |
| 4. Parse genetic map | `parse_genetic_map` | `scripts/parse_genetic_map.py` |
| 5. Pileup at het SNPs | `pileup_snps_bulk_mode1b` | cellsnp-lite |
| 6. Phase and concat (joint, all bulk assays) | `phase_and_concat_bulk` | `scripts/phase_and_concat_bulk.py` |
| 7. Compute read depth | `run_mosdepth` | mosdepth |
| 8. Bias correction | `rd_correct` | `scripts/rd_correct.py` |
| 9. Adaptive binning | `combine_counts` | `scripts/combine_counts.py` |

All bulk replicates are piled up against one shared phased het-SNP VCF, so step 6 runs **once** and builds a single joint allele matrix (under `allele_dir/{stream}/` where stream is `bulkWGS` for the WGS family or `bulkWES`, one pseudobulk column per replicate). Step 9 reads it directly — no per-assay union. Depth/RDR (steps 7–8) stay per-assay.

**Outputs** (`bb_dir/{stream}/MSR{msr}/`, one subdir per `min_snp_reads` value — the sweep runs in a single job with shared preprocessing): `bb.tsv.gz`, `bb.{Tallele,Aallele,Ballele,rdr,depth}.npz`, `sample_ids.tsv`

---

## `single_cell_genotyping`

Assays: `scRNA`, `scATAC`, `VISIUM`, `VISIUM3prime`

| Step | Rule | Script / Tool |
|------|------|---------------|
| 1. Pseudobulk genotyping | `genotype_snps_pseudobulk_mode1b` | cellsnp-lite |
| 2. Annotate SNPs | `annotate_snps_pseudobulk` | `scripts/annotate_snps_pseudobulk.py` |
| 3. Phase SNPs (per chr) | `phase_snps_{eagle,shapeit,longphase}` | eagle2 / shapeit5 / longphase |
| 4. Concat phased VCFs | `concat_and_extract_phased_het_snps` | bcftools concat + view |
| 5. Parse genetic map | `parse_genetic_map` | `scripts/parse_genetic_map.py` |
| 6. Single-cell pileup | `pileup_snps_nonbulk_mode1a` | cellsnp-lite |
| 7. Build RNA AnnData (RNA-family only; feeds the RNA `bb.Xcount.npz`) | `process_rna_anndata` | `scripts/process_rna_anndata.py` |
| 8. Phase and concat | `phase_and_concat_nonbulk` | `scripts/phase_and_concat_nonbulk.py` |
| 9. Adaptive binning (+ per-assay `bb.Xcount.npz`: scATAC from raw fragments, RNA from the h5ad) | `combine_counts_nonbulk` | `scripts/combine_counts_nonbulk.py` |

**Outputs** (all per-assay under `bb_dir/{assay}/MSR{msr}/`, one subdir per `min_snp_reads` value — the sweep runs in a single job with shared preprocessing): `bb.tsv.gz` (the one shared grid, joint across the sample's assays, duplicated into each sub-dir), `sample_ids.tsv` (likewise duplicated), `bb.{Tallele,Aallele,Ballele}.npz` (bins × cells), `bb.Xcount.npz` (per-cell native counts per bb bin — scATAC from raw fragments, scRNA/VISIUM UMIs from the h5ad), `multi_snp.*` (MSR-independent), `barcodes{,.full}.tsv.gz`

---

## `copytyping_preprocess`

Assays: `scRNA`, `scATAC`, `VISIUM`, `VISIUM3prime` (requires pre-computed het SNP VCF)

| Step | Rule | Script / Tool |
|------|------|---------------|
| 1. Single-cell pileup | `pileup_snps_nonbulk_mode1a` | cellsnp-lite |
| 2. Build AnnData (RNA only) | `process_rna_anndata` | `scripts/process_rna_anndata.py` |
| 3. Phase and concat | `phase_and_concat_nonbulk` | `scripts/phase_and_concat_nonbulk.py` |
| 4. Adaptive binning | `combine_counts_nonbulk` | `scripts/combine_counts_nonbulk.py` |
| 5. CNV segmentation | `cnv_segmentation` | `scripts/cnv_segmentation.py` |

**Outputs** (`bb_dir/{assay_type}/`): `cnv_segments.tsv`, `bb.{Tallele,Aallele,Ballele,Xcount}.npz` (scATAC `Xcount` from raw fragments, RNA from the h5ad), `barcodes{,.full}.tsv.gz`, `sample_ids.tsv`
