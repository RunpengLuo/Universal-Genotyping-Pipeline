# TODO

- Distinguish germline Het from Hom-alt SNPs from high purity tumor sample without matched-normal sample. Adapt https://github.com/raphael-group/hetdetect.
    - retrieve population ALT frequency as prior genotype info. high ALT freq indicates likely hom-alt 
- Streaming remote data rather than downloading them. Currently `cellsnp-lite` don't allow URL inputs.
- Test suite. end-to-end workflow unit testing.

## Unify WGS/WES windows and allow mixing in bulk genotyping (IMPLEMENTED)

Done: `build_segment_bed` (region_id arm + seg_id chunk) + per-stream window BEDs
(`{wgs,wes}_windows.bed.gz`); `combine_counts` bins on the WGS grid grouped by `seg_id`,
projects WES depth by midpoint, and uses a per-column read threshold (`min_snp_reads` WGS,
`min_snp_reads_wes` WES) with `max_blocksize` gated behind it; the `validate_mode` guard is
lifted; one `bb_dir/MSR{msr}/bulk/`; per-stream `rdr_base_dataset_id` validation. Verified
by DAG tests only (incl. a mixed WGS+WES case). Follow-ups:

- Real-data validation: no execution-level test exists yet; run a real mixed WGS+WES bulk
  sample end-to-end and inspect the projected WES RDR/BAF per bin.
- Residual NaN edge: a `max_blocksize` cut in a target-free desert (no WES window in the
  span) can still yield a WES-empty bin; the existing NaN-row drop removes it. Revisit only
  if it discards useful WGS-covered bins.
