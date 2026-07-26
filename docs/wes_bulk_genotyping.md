# WES bulk genotyping, end to end

WES is processed **identically to WGS**: the same `segment.bed`/`window_size` grid, no
capture-target file, and the same `min_snp_reads`. WGS/WGS-lr/WES of one patient bin
together on one shared grid. This page notes only what a WES reader should keep in mind;
the full stage flow is in [`workflow.md`](workflow.md).

## Stage overview

```
genotype -> phase -> pileup -> build windows -> count reads -> combine counts
 (SNPs)     (0|1)    (allele    (segment.bed     (mosdepth +    (bin on the
                     counts)    tiled grid)      LOWESS)        shared grid -> bb)
```

## What is the same as WGS

- **Windows:** `build_window_bed` tiles `aux/segment.bed` at
  `params_build_windows.window_size` (default 1 kbp) into one `aux/windows.bed.gz`, shared
  by every bulk assay. No exon tiling, no `wes_targets_bed`.
- **Binning:** `adaptive_segmentation` grows each bin until every tumor column (WGS/WGS-lr/WES
  alike) meets `min_snp_reads`, capped by `max_blocksize`; gene-aware if enabled.
- **Genotype/phase/pileup:** assay-agnostic. With a matched normal WES dataset,
  `bulk_genotyping_strategy: auto` -> `bcftools`; without one -> `tumor_hmm`.

## What to keep in mind

> [!NOTE]
> WES reads concentrate on capture targets, but ~10-50% are off-target, so genome-wide
> 1 kbp bins are not empty and adaptive binning still reaches the read target per bin. WES
> contributes BAF from het SNPs (mostly exonic) and depth/RDR from all its reads.

> [!WARNING]
> RD-correction caveat: WES depth carries capture-enrichment structure, so its
> post-correction RDR is noisier than WGS. The per-sample LOWESS GC/MAP/REPLI fit self-adapts
> and `params_count_reads.routlier`/`doutlier` trim spikes, but inspect
> `qc/rd_correction.bulkWES.pdf` before trusting WES RDR. This is documented, not
> special-cased -- WES uses the same `rd_correct` code as WGS.

> [!TIP]
> A pre-built `window_bed` (e.g. `resources/data/windows.1kbp.hg38.bed.gz`) is honored for
> WES exactly as for WGS: when set (and no `breakpoint_bedpe`), `build_window_bed` is skipped
> and the file is read directly.
