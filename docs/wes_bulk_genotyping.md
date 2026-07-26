# WES-only bulk genotyping, end to end

How a `bulk_genotyping` run with `assay_types: ["bulkWES"]` (genotyping on, i.e.
`het_snp_vcf` unset) flows from BAMs to `bb.tsv.gz`.

## Stage overview

```
genotype -> phase -> pileup -> build windows -> count reads -> combine counts
 (SNPs)     (0|1)    (allele    (267bp WES       (mosdepth +    (bin on WES
                      counts)    target grid)     LOWESS)        windows -> bb)
```

## Parse (`workflow/scripts/parse_workflow_args.py`)

- `window_streams = ["wes"]` (bulkWES maps to the `wes` stream).
- `wes_targets_bed` on a bulkWES record is **required** (asserted at parse) unless a
  pre-built `window_bed` is supplied.
- `bulk_genotyping_strategy: auto` -> `bcftools` when a `normal` WES dataset exists,
  else `tumor_hmm` (beta-binomial het caller).
- `segment_bed = aux/segment.bed`: arm `region_id` + `seg_id` chunk, from
  `region_bed` - blacklist - `breakpoint_bedpe`.

## 1. Genotype -- `workflow/rules/genotype_snps.smk` (`genotype_snps_bulk`, `bcftools` env)

`bcftools mpileup` on the selected dataset's WES BAM, restricted to the population panel
positions `snp_targets/target.chr{c}.pos.gz` (`--skip-indels`, `-a INFO/AD,AD,DP`) ->
`bcftools call -m` -> keep biallelic SNPs with `QUAL >= min_qual`, `GT="alt"`,
`FMT/DP >= min_dp`.

Off-target panel positions have ~0 WES coverage and are dropped by the `DP` filter, so
het SNPs effectively land in captured regions. Output: `snps/chr{c}.vcf.gz`.

> [!NOTE]
> With no matched normal, the `tumor_hmm` strategy pileups the tumor(s) over `snp_panel`
> and calls het SNPs with a multivariate beta-binomial HMM instead of `bcftools call`.

## 2. Phase -- `workflow/rules/phase_snps.smk` (`eagle` / `shapeit` / `longphase` env)

Assay-agnostic: phase `snps/chr{c}.vcf.gz` against the reference panel, then
`concat_and_extract_phased_het_snps` -> `phase/phased_het_snps.vcf.gz` (het `0|1`/`1|0`).

## 3. Build windows -- `workflow/rules/build_windows.smk` (`base` / `ucsc` envs)

- `build_segment_bed` -> `aux/segment.bed`.
- `build_window_bed` (wes stream): `generate_wes_windows` merges the capture targets and
  **tiles them at `window_size_wes` (default 267 bp)**, assigns `region_id`/`seg_id` by
  window midpoint (windows whose midpoint falls off-segment/blacklist are dropped), and
  annotates GC (always) + MAP (if `mappability_bed`) + REPLI (if Repli-seq) ->
  `aux/wes_windows.bed.gz`.
- `get_assay_window_bed("bulkWES")` returns `aux/wes_windows.bed.gz`.

## 4. Pileup -- `workflow/rules/pileup_snps.smk` (`pileup_snps_bulk_mode1b`, `cellsnp` env)

`cellsnp-lite` over each WES dataset's BAM at the phased het SNPs -> per-dataset AD/DP,
consolidated into `allele_dir/bulk/snp.{Tallele,Aallele,Ballele}.npz` + `snps.tsv.gz` +
`sample_ids.tsv` (one joint SNP grid over all bulk datasets).

## 5. Count reads -- `workflow/rules/count_reads.smk` (`mosdepth` / `base` envs)

- `run_mosdepth --by wes_windows` per dataset -> per-window depth.
- `rd_correct`: LOWESS GC / MAP / REPLI correction **on the WES window grid**, filtered
  to `chromosomes` -> `pileup/bulkWES/window.dp.npz` + `window.tsv.gz` + depth stats + QC.

## 6. Combine counts -- `workflow/scripts/combine_counts.py` (`base` env)

The WES-only branch (`combine_counts.py`, ~line 138):

```python
grid_idxs = [i for i, at in enumerate(bulk_assays) if at != "bulkWES"]
if not grid_idxs:                      # WES-only: no WGS grid exists
    grid_idxs = list(range(len(window_df_list)))   # bin ON the WES windows
```

So the **WES 267 bp windows are the binning grid** -- there is no midpoint projection
(projection only applies to WES *columns* in a mixed WGS+WES run). Then:

- `adaptive_segmentation` groups windows by `seg_id` and grows each bin until every tumor
  column meets its read target -- WES columns use `min_snp_reads_wes` (ANDed with the
  swept `min_snp_reads`) -- capped by `max_blocksize`; gene-aware if enabled.
- Allele counts summed per bin -> BAF; corrected WES-window depth aggregated per bin
  (identity mapping here); RDR = tumor bin depth / base (`rdr_base_dataset_id` or median).

Output: `bb_dir/MSR{msr}/bulk/bb.tsv.gz` + `bb.{Tallele,Aallele,Ballele,depth,rdr}.npz` +
`sample_ids.tsv` + `combine_counts.bulk.MSR{msr}.pdf`.

## Things worth knowing

> [!IMPORTANT]
> WES bins on capture-target-tiled 267 bp windows, not the 1 kbp WGS grid. `window_size_wes`
> and `min_snp_reads_wes` are the knobs that tune WES bin sizes.

> [!WARNING]
> If you set `window_bed` (a pre-built WGS 1 kbp grid) for a WES-only run,
> `use_prebuilt_windows` makes `get_assay_window_bed` return that WGS grid for the wes
> stream too -- so it bins on 1 kbp WGS windows and loses the capture-target resolution.
> For a proper WES run, leave `window_bed` unset so the wes windows build from
> `wes_targets_bed`.
