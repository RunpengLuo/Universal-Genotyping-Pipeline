# TODO

## MSR vector to support joint whole-genome / targeted segmentation

Replace the scalar `min_snp_reads` with a per-dataset vector derived from a target BAF standard
error, so assays of unequal read supply (WGS + WES, WGS + CRISPR-targeted) can share one bb grid.
Design, model, and implementation steps: `.claude/msr_vector.md`.

## slurm support

## Test pipeline

Goal: after a user clones + installs the pipeline, one command runs a real end-to-end
test that exercises a full mode, not a per-rule dry run.

Two layers:
- **Dry-run (DAG) tests** (`tests/test_dryrun.py`, `tests/test_remote.py`): built, all
  modes x JSON/TSV; execute no rule, need no data. This is what CI runs.
- **End-to-end tests**: a real `snakemake` run to final `bb_dir/` outputs. Not in the
  repo. Sample sheets and configs for real cases live outside it, since the inputs are
  whole-genome remote alignments and no hosted runner can execute them.

- [ ] Decide where end-to-end case definitions live (separate repo, or a gitignored
      working dir) and how references are staged.
- [ ] Cases to cover: `single_cell_genotyping`, `copytyping_preprocess`, spatial.
- [ ] One bulk WGS+WES case run end to end, inspecting per-bin WES RDR/BAF: the shared
      bulk grid and the on/off-target depth fit are both verified by DAG tests only.

Harness: if end-to-end returns, gate it behind a manual/scheduled trigger and keep the
dry-run tests on every push.

## Tumor-only genotyping (no matched normal)

Today every mode needs a normal to call germline SNPs: `genotype_dataset_ids` defaults to the
first normal, and a tumor-only run genotypes off tumor reads with a WARN.

- [ ] Distinguish germline het from hom-alt SNPs in a high-purity tumor without a matched
      normal. Adapt https://github.com/raphael-group/hetdetect.
    - retrieve population ALT frequency as prior genotype info. high ALT freq indicates likely hom-alt

## Segmentation parameter selection

`min_snp_reads` (and `min_total_reads`) are chosen by hand: a list sweeps the grid and the
user picks a point off the QC PDFs. The `select_segmentation` script that scored a
`(min_snp_reads, max_blocksize)` grid is not in the tree, and `max_blocksize` itself is gone.

- [ ] Automated model selection for the segmentation parameters against sequencing coverage
      and segmentation variance, i.e. recommend one grid point instead of a sweep.
- [ ] Concretely: pick a default MSR at the elbow of lag-1 RDR/BAF dispersion vs bin count,
      and record the pick alongside the outputs (was a TODO comment in `combine_counts.py`,
      pointing at a `docs/combine_counts_pseudocode.md` that is not in the repo).

## Within-bb EM phasing

Implemented as `params_combine_counts.phase_em`; model and parameters in
`workflow/scripts/script_utils/phase_em.py` and
[reference](reference.md#params_combine_counts).

> [!IMPORTANT]
> Two designs are superseded and should not be revived without new evidence: the HMM
> (`.claude/phase-hmm-within-bin.md`) and the later per-SNP transition prior from the
> genetic map, with forward-backward and Viterbi. Both were built and measured against the
> independent per-SNP latent and lost. A run prior makes the flips coherent, and a
> coherent run of flips manufactures more false imbalance in a balanced bb than the same
> number of scattered flips.

- [ ] `phaser: "none"` - genotype, skip panel and long-read phasing, carry the unphased
      het SNPs into pileup, for runs where no phasing prior exists. Design in
      `.claude/skip_phasing_and_bb_em.md` (its Part 2/3, the `phase_correction` enum, is
      superseded by `phase_em`). Includes the `estimate_switchprobs_PS` `KeyError: 'PS'`
      fix, which is still live: `build_adaptive_bins` only creates `bbs["PS"]` when the
      SNP frame has one.
- [ ] A/B comparison of the EM against HATCHet2 `bb` on a shared input.
- [ ] Re-derive `phase_em_tau` and `phase_em_min_llr` on a second dataset; both are
      currently set from one.

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

## Remote input

`remote_mode: stream` is implemented for bulk: bcftools/mosdepth/longphase read remote
BAM/CRAM directly and fetch only the config `chromosomes` via index jumps; the default
stays `storage` (whole-file download). Single-cell and copytyping cannot stream, since
`cellsnp-lite` rejects URLs in its `access(F_OK)` guard.

- [ ] Validate `##idx##` remote-index support and numeric parity on a real URL BAM.
- [ ] Stream the SNP panel the same way, rather than fetching the whole file:
      `bcftools query -f '%CHROM\t%POS\n' -r chr22 <URL>` against
      `https://ftp.ncbi.nih.gov/snp/organisms/`. Measure what tabix over a remote panel
      costs against a local one before making it a default.
