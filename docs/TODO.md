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
- Streaming remote data rather than downloading them. Bulk pileup now uses `bcftools`
  (`pileup_snps_bulk_bcftools`), which is htslib/URL-capable; single-cell still uses
  `cellsnp-lite`, which rejects URL inputs (its `access(F_OK)` guard). Whole-file `storage()`
  download is still the default (streaming = whole-file transfer + fragility tradeoffs).

## Within-bin BAF phasing (`phase_hmm.py`)

Bulk `combine_counts` collapses CNA/LOH BAF toward 0.5 (seen on `hatchet2_chr22_simulation`
dbSNP151; panel/phaser ruled out). Two causes: (1) `detect_phase_flips` fragments SNPs into
thousands of phase-groups, and the per-group bin-count floor in `_bin_windows_numba` makes
~4-SNP bins so `min_snp_reads` never binds; (2) `apply_phase_to_mat` orients A/B by one
per-SNP `PHASE` bit with no within-bin re-orientation, so unfolded bin BAF is exactly 0.500.
HATCHet2 gets ~220 SNPs/bin and BAF ~0.20 in LOH via a per-bin EM (phase latent shared across
samples).

Fix: a new helper `workflow/scripts/phase_hmm.py`, a within-bin multi-sample beta-binomial
phase HMM (K=1 per bin), gated by `params_combine_counts.phase_correction`
(`none|flip_split|bin_hmm`). Under `bin_hmm`, drop the phase-group split so bins reach the
read target, then a per-SNP phase HMM (latent shared across samples, LD-derived switch/stay
transitions, tau calibrated on the matched normal) re-orients SNPs inside each bin and stores
unfolded phased-frame counts. Optional layer 2 (`cross_bin_phasing: dp`) is a 2-state Viterbi
over bins per region for cross-bin orientation. Emission/recursion port HATCHet3
(`cluster_bins/hmm`); per-bin-EM structure follows HATCHet2. Full design, model equations,
combine_counts flow, and validation: `~/.claude/plans/phase-hmm-within-bin.md`.

- [ ] Phase A: `phase_hmm.py` + `snp_switchprobs`; wire `phase_correction` into bulk
      `combine_counts`; config/const/parser; unit test; A/B comparison vs HATCHet2 `bb`; docs.
- [ ] Phase B: single-cell `combine_counts_nonbulk` reuse; make `bin_hmm` default after A/B.

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

## Unified bulk grid: WES treated exactly like WGS (IMPLEMENTED)

Done: `build_segment_bed` (region_id arm + seg_id chunk) + one shared window BED
(`aux/windows.bed.gz`, tiled from `segment.bed` at `window_size`). Every bulk assay
(WGS/WGS-lr/WES) bins on that one grid grouped by `seg_id`, with a single read target
`min_snp_reads` and `max_blocksize` gated behind it; one `bb_dir/MSR{msr}/bulk/`. The
earlier WES stream (per-record `wes_targets_bed`, 267 bp exon tiling, WES-onto-WGS depth
projection, `min_snp_reads_wes`) was removed. Verified by DAG tests only. Follow-ups:

- Real-data validation: no execution-level test exists yet; run a real WGS+WES bulk sample
  end-to-end and inspect per-bin WES RDR/BAF.
- WES RD-correction: WES depth carries capture-enrichment structure; confirm the per-sample
  LOWESS fit + `routlier`/`doutlier` handle it, or flag over-correction (see the RD bias
  correction section above).

## streaming panel
`https://ftp.ncbi.nih.gov/snp/organisms/`, `tabix-streams just that region from the remote file — never downloading the ~16 GB whole thing (genotype_snps.py:232): bcftools query -f '%CHROM\t%POS\n' -r chr22 <URL> -> target_chr22.pos.gz.`
what is usual time cost for tabix snp positions from local/online snp panel -> genotyping?