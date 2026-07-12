- Streaming remote alignments is not possible today: `cellsnp-lite` rejects a URL before htslib
  opens it (`[E::check_args] '<url>' does not exist`), and it runs the pileup stage in every mode, so
  every alignment must be downloaded anyway. `samtools`, `mosdepth`, and `bcftools mpileup` all read
  remote BAMs fine (verified against a 94 GB GIAB BAM: a 10 kb region in 1.4 s via range requests).
  Revisit if cellsnp-lite gains URL support.
- Test suite: only dry-run (DAG) tests exist today, under `tests/`. They cover sample-file
  parsing/validation, rule wiring, storage() handling of remote inputs, and each mode's targets,
  from both the JSON and the legacy TSV format. No rule is ever executed. Still to build, in order
  of cost:
  - Parser unit tests: `workflow/scripts/parse_workflow_args.py` is stdlib-only and importable; test the
    JSON/TSV record equivalence, every `validate_records` error, `expand_ranger_dir`, and
    `get_visium_layout` directly (milliseconds, no snakemake).
  - Real execution on tiny fixtures: commit a small synthetic reference plus normal/tumor
    alignments (~100 KB; slice a region with `samtools view -b ref chr22:1-100000`, or see
    https://github.com/nf-core/test-datasets). Bulk can run `--until run_mosdepth` to avoid needing
    a phasing panel; `copytyping_preprocess` runs end to end because it takes a pre-computed
    `het_snp_vcf` + `bb_file` and skips genotyping and phasing. Single-cell needs a 10x-shaped
    `filtered_feature_bc_matrix.h5` and `atac_fragments.tsv.gz`.
  - Per-rule tests: once one real run succeeds, `snakemake --generate-unit-tests` writes a pytest
    case per rule into `.tests/unit/`. It requires a prior successful run and small data:
    https://snakemake.readthedocs.io/en/stable/snakefiles/testing.html
