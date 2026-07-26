# TODO

## Test pipeline

Goal: after a user clones + installs the pipeline, one command runs a real end-to-end
test that exercises a full mode, not a per-rule dry run.

Two layers:
- **Dry-run (DAG) tests** (`tests/test_dryrun.py`, `tests/test_remote.py`): built, all
  modes x JSON/TSV; execute no rule, need no data. Keep as the fast default.
- **End-to-end tests** (`tests/data/<case>/`): a real `snakemake` run to final
  `bb_dir/` outputs. One directory per case: `samples.json` + `config.yaml` (+ reference
  staging notes). Later cases: single_cell_genotyping, copytyping_preprocess, spatial.

Cases:
- [x] `hatchet2_chr22_simulation`: HATCHet demo-complete data, hg19/chr22, 1 matched
      normal + 3 tumors (alignments streamed from Zenodo 4046906), phasing panel from
      Zenodo 6709541, pre-built hg19 window BED (GC + MAP + REPLI) + blacklist + gtf.
      Source dir holds only `README.md` + wrapper scripts + `config.yaml` +
      `samples.json`; `prepare_refs.sh` stages references and `run_pipeline.sh` /
      `check_outputs.sh` run + verify. Runtime (staged refs + outputs) under gitignored
      `.test-run/`. CI: `.github/workflows/it-hatchet2_chr22_simulation.yml` (run locally
      via `act`).
- [ ] `single_cell_genotyping`, `copytyping_preprocess`, spatial cases.

Harness: GitHub Actions workflow (`.github/workflows/`) that installs the env
(conda/mamba), stages references, runs the case, and asserts outputs; gate the heavy
end-to-end job behind a manual/scheduled trigger, keep dry-run tests on every push.

## Others
- Distinguish germline Het from Hom-alt SNPs from high purity tumor sample without matched-normal sample. Adapt https://github.com/raphael-group/hetdetect.
    - retrieve population ALT frequency as prior genotype info. high ALT freq indicates likely hom-alt 
- Streaming remote data rather than downloading them. Currently `cellsnp-lite` don't allow URL inputs.

## RD bias correction (potential over-correction)

Observed on `hatchet2_chr22_simulation` (chr22-only, matched normal + 3 tumors): GC/RT
correction lowers the GC/RT correlation but raises RD dispersion (T2 GC MAD 10.80 -> 11.58),
i.e. removes weak bias while injecting the covariate shape (the RT replication wave).
Contributing factors: weak biases (GC r ~ 0.12, RT r ~ 0.1-0.24), per-sample quadratic fit
over CNA-contaminated bins (`rd_correct_utils.py`; highest-CNA sample hurt most), the
`RT + RT**2` term imprinting the replication-timing wave, and a single-chromosome fit.

- [ ] `gc_correct`, `rt_correct`: options `[true, false, auto]`, resolved per dataset, not
      one global switch.
- [ ] `auto`: choose each covariate from the data - explained-variance / correlation cutoff
      per covariate, matched-normal presence (tumor/normal RDR already cancels shared GC/RT
      bias), min bin count / genome coverage - and drop a covariate below threshold.
- [ ] Fit the correction on CN-neutral bins, or derive it from the normal and apply to the
      tumors, instead of a per-tumor fit over CNA-altered bins.
- [ ] QC gate: flag or revert a correction when it increases a sample's RD dispersion
      (MAD after > MAD before).
- [ ] Fit genome-wide (full GC/RT range, enough bins); warn when a single chromosome or low
      bin count makes the quadratic fit unreliable.
- [ ] Lower-order or regularized model when a covariate's bias is weak.

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

## streaming panel
`https://ftp.ncbi.nih.gov/snp/organisms/`, `tabix-streams just that region from the remote file — never downloading the ~16 GB whole thing (genotype_snps.py:232): bcftools query -f '%CHROM\t%POS\n' -r chr22 <URL> -> target_chr22.pos.gz.`
what is usual time cost for tabix snp positions from local/online snp panel -> genotyping?