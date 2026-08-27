"""tumor_genotyping_mode picks how a tumor is genotyped with no matched normal.

Two strategies, both tumor-only: `vaf_cutoff` (the default) thresholds depth, minor-allele
reads and VAF; `clonal_loh_hmm` re-genotypes each site from its neighbourhood, since under
clonal LOH a gHET collapses to one allele and reads like a gHOM.

Neither applies when no genotyped bulk dataset is a tumor: the caller's germline GT passes
through. That rule is bulk-only, cellsnp-lite emitting no GT of its own.
"""

from conftest import dryrun

# the key nests under params_genotype_snps, so --config nests it too
HMM = "params_genotype_snps={'tumor_genotyping_mode': 'clonal_loh_hmm'}"
VAF = "params_genotype_snps={'tumor_genotyping_mode': 'vaf_cutoff'}"
LEGACY = "params_genotype_snps={'apply_clonal_loh_hmm': True}"
BOGUS = "params_genotype_snps={'tumor_genotyping_mode': 'nonsense'}"


def bulk(workspace, extra=()):
    return dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=extra,
    )


def test_no_tumor_genotyped_passes_the_calls_through(workspace):
    """The default auto-picks the normal, so neither mode applies and GT survives."""
    proc = bulk(workspace)
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "tumor_genotyping_mode=none" in proc.stdout
    assert "--variants-only" in proc.stdout
    assert 'GT="alt"' in proc.stdout
    assert "--keep-alts" not in proc.stdout


def test_vaf_cutoff_applies_when_a_tumor_is_genotyped(workspace):
    """Naming the tumor engages vaf_cutoff; the calling flags do not change."""
    proc = bulk(workspace, ['genotype_dataset_ids=["D1"]', VAF])
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "tumor_genotyping_mode=vaf_cutoff" in proc.stdout
    assert "--variants-only" in proc.stdout
    assert 'GT="alt"' in proc.stdout
    assert "--keep-alts" not in proc.stdout


def test_hmm_keeps_every_callable_site(workspace):
    """clonal_loh_hmm swaps --variants-only for --keep-alts and drops the GT filter."""
    proc = bulk(workspace, [HMM])
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "tumor_genotyping_mode=clonal_loh_hmm" in proc.stdout
    assert "--keep-alts" in proc.stdout
    assert "--variants-only" not in proc.stdout
    assert 'GT="alt"' not in proc.stdout


def test_only_the_hmm_declares_its_aux_tables(workspace):
    """No rule reads them, so a run that fits no chain must not promise them either."""
    off = bulk(workspace, ['genotype_dataset_ids=["D1"]', VAF])
    on = bulk(workspace, [HMM])
    assert off.returncode == 0 and on.returncode == 0
    assert "clonal_loh_hmm.segments.tsv" not in off.stdout
    assert "clonal_loh_hmm.segments.tsv" in on.stdout


def test_panel_alleles_are_always_constrained(workspace):
    """The panel's REF/ALT are fixed in every mode; there is no opting out."""
    for extra in ((), [HMM], ['genotype_dataset_ids=["D1"]', VAF]):
        proc = bulk(workspace, extra)
        assert proc.returncode == 0, proc.stderr[-1500:]
        assert "--constrain alleles" in proc.stdout
        assert "%CHROM\\t%POS\\t%REF,%ALT" in proc.stdout


def test_hmm_requires_tumor_datasets(workspace):
    """The rescue is tumor-only, so pointing it at a normal is a parse error."""
    proc = bulk(workspace, ['genotype_dataset_ids=["N1"]', HMM])
    assert proc.returncode != 0
    combined = proc.stdout + proc.stderr
    assert "clonal_loh_hmm" in combined
    assert "N1" in combined


def test_hmm_is_rejected_off_bulk(workspace):
    """cellsnp-lite pseudobulk has no per-arm depth to fit, so the key is refused."""
    proc = dryrun(
        workspace,
        workspace["sc_json"],
        "S1",
        "single_cell_genotyping",
        ["scRNA"],
        extra=[HMM],
    )
    assert proc.returncode != 0
    combined = proc.stdout + proc.stderr
    assert "clonal_loh_hmm" in combined
    assert "single_cell_genotyping" in combined


def test_single_cell_always_thresholds_counts(workspace):
    """No caller GT exists off bulk, so vaf_cutoff runs whatever the sample_type."""
    proc = dryrun(
        workspace,
        workspace["sc_json"],
        "S1",
        "single_cell_genotyping",
        ["scRNA"],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "tumor_genotyping_mode=vaf_cutoff" in proc.stdout


def test_legacy_key_is_a_parse_error(workspace):
    """apply_clonal_loh_hmm is gone; carrying it must fail loudly, not be ignored."""
    proc = bulk(workspace, [LEGACY])
    assert proc.returncode != 0
    combined = proc.stdout + proc.stderr
    assert "apply_clonal_loh_hmm" in combined
    assert "tumor_genotyping_mode" in combined


def test_unknown_mode_is_a_parse_error(workspace):
    """The value is closed; a typo names the allowed set rather than falling through."""
    proc = bulk(workspace, [BOGUS])
    assert proc.returncode != 0
    combined = proc.stdout + proc.stderr
    assert "tumor_genotyping_mode" in combined
    assert "nonsense" in combined
