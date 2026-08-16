"""apply_clonal_loh_hmm decides how bulk calls are made and post-processed.

false (default) calls variant-only with the `GT="alt"` filter and post_genotype_snps
symlinks the calls through. true keeps every callable panel site with `--keep-alts` and
re-genotypes it with the clonal-LOH HMM, since under clonal LOH a gHET collapses to one
allele and reads like a gHOM. It is a tumor-only rescue, so it requires tumor datasets;
the converse is not an error, a tumor may be genotyped variant-only.
"""

from conftest import dryrun

# apply_clonal_loh_hmm lives inside params_genotype_snps, so --config nests it
HMM_ON = "params_genotype_snps={'apply_clonal_loh_hmm': True}"


def bulk(workspace, extra=()):
    return dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=extra,
    )


def test_variant_only_calling_is_the_default(workspace):
    """Without the HMM the calls are variant-only, GT-filtered, and passed through."""
    proc = bulk(workspace)
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "apply_clonal_loh_hmm=False" in proc.stdout
    assert "--variants-only" in proc.stdout
    assert 'GT="alt"' in proc.stdout
    assert "--keep-alts" not in proc.stdout


def test_hmm_keeps_every_callable_site(workspace):
    """apply_clonal_loh_hmm swaps --variants-only for --keep-alts, drops the GT filter."""
    proc = bulk(workspace, [HMM_ON, 'genotype_dataset_ids=["D1"]'])
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "apply_clonal_loh_hmm=True" in proc.stdout
    assert "--keep-alts" in proc.stdout
    assert "--variants-only" not in proc.stdout
    assert 'GT="alt"' not in proc.stdout


def test_panel_alleles_are_always_constrained(workspace):
    """The panel's REF/ALT are fixed in both modes; there is no opting out."""
    for extra in ((), [HMM_ON, 'genotype_dataset_ids=["D1"]']):
        proc = bulk(workspace, extra)
        assert proc.returncode == 0, proc.stderr[-1500:]
        assert "--constrain alleles" in proc.stdout
        assert "%CHROM\\t%POS\\t%REF,%ALT" in proc.stdout


def test_hmm_requires_tumor_datasets(workspace):
    """The rescue is tumor-only, so declaring it over a normal is a parse error."""
    proc = bulk(workspace, ['genotype_dataset_ids=["N1"]', HMM_ON])
    assert proc.returncode != 0
    combined = proc.stdout + proc.stderr
    assert "apply_clonal_loh_hmm=True" in combined
    assert "N1" in combined


def test_tumor_may_be_genotyped_without_the_hmm(workspace):
    """The check is one-directional: variant-only calling off a tumor is legitimate."""
    proc = bulk(workspace, ['genotype_dataset_ids=["D1"]'])
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "apply_clonal_loh_hmm=False" in proc.stdout
    assert "--variants-only" in proc.stdout


def test_hmm_flag_is_ignored_off_bulk(workspace):
    """cellsnp-lite genotypes from counts either way, so the key never reaches it."""
    proc = dryrun(
        workspace,
        workspace["sc_json"],
        "S1",
        "single_cell_genotyping",
        ["scRNA"],
        extra=[HMM_ON],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "--keep-alts" not in proc.stdout
