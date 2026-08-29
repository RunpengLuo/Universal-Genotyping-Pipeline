"""Genotyping a tumor, and the clonal-LOH detection that a tumor-only run turns on.

The auto-pick prefers a normal, whose germline GT passes through; it falls through to a
tumor when the sample carries none, and the VAF cutoff then genotypes it.
`params_combine_counts.detect_loh_tumor_cell_line` is independent of all that: it is a binning decision,
made in combine_counts from het-SNP density, since one site's counts cannot tell a gHET
inside LOH from a gHOM.
"""

from conftest import dryrun

DETECT_LOH = "params_combine_counts={detect_loh_tumor_cell_line: True}"


def bulk(workspace, extra=()):
    return dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=extra,
    )


def test_default_genotypes_the_normal_and_passes_the_calls_through(workspace):
    """A sample with a normal auto-picks it, and GT survives untouched."""
    proc = bulk(workspace)
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "tumor_genotyping_mode=none" in proc.stdout


def test_vaf_cutoff_applies_whenever_a_tumor_is_genotyped(workspace):
    """Naming the tumor engages the cutoff; nothing else has to be set."""
    proc = bulk(workspace, ['genotype_dataset_ids=["D1"]'])
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "tumor_genotyping_mode=vaf_cutoff" in proc.stdout


def test_detect_loh_is_independent_of_what_is_genotyped(workspace):
    """It is a binning switch, so it declares the QC page over the normal's calls too."""
    proc = bulk(workspace, [DETECT_LOH])
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "tumor_genotyping_mode=none" in proc.stdout
    assert "detect_loh.bulk.pdf" in proc.stdout


def test_only_a_detect_loh_run_declares_the_qc_page(workspace):
    """A run with no clonal LOH to find must not promise the QC page."""
    off = bulk(workspace)
    on = bulk(workspace, [DETECT_LOH])
    assert off.returncode == 0 and on.returncode == 0
    assert "detect_loh.bulk.pdf" not in off.stdout
    assert "detect_loh.bulk.pdf" in on.stdout


def test_panel_alleles_are_always_constrained(workspace):
    """The panel's REF/ALT are fixed in every mode; there is no opting out."""
    for extra in ((), [DETECT_LOH], ['genotype_dataset_ids=["D1"]']):
        proc = bulk(workspace, extra)
        assert proc.returncode == 0, proc.stderr[-1500:]
        assert "--constrain alleles" in proc.stdout
        assert "%CHROM\\t%POS\\t%REF,%ALT" in proc.stdout


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
