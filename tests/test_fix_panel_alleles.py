"""fix_panel_allele constrains bulk genotyping to the panel's REF/ALT.

On by default: `bcftools call` gets `--constrain alleles` with the panel's alleles piped
in, so reads carrying any other allele are discarded and a site can only be genotyped as
the panel spells it. Off, `-T` gives positions only, REF comes from the reference and ALT
from the reads.
"""

from conftest import dryrun


def test_panel_alleles_constrain_the_call_by_default(workspace):
    """Default run renders --constrain alleles, fed by a stream off the panel."""
    proc = dryrun(
        workspace, workspace["bulk_json"], "T1", "bulk_genotyping", ["bulkWGS"]
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "--constrain alleles" in proc.stdout
    # the alleles are streamed from the panel, keyed to this job's chromosome
    assert "%CHROM\\t%POS\\t%REF,%ALT" in proc.stdout
    assert "--regions chr22" in proc.stdout


def test_panel_alleles_can_be_disabled(workspace):
    """fix_panel_allele=False falls back to calling ALT from the reads."""
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=["fix_panel_allele=False"],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "--constrain alleles" not in proc.stdout


def test_no_matched_normal_forces_panel_alleles(workspace):
    """Genotyping a tumor with no normal must not let a somatic allele become ALT."""
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=['genotype_dataset_ids=["D1"]', "fix_panel_allele=False"],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "--constrain alleles" in proc.stdout, "the explicit False must be overridden"
    assert "no matched normal to genotype" in proc.stdout + proc.stderr


def test_panel_alleles_ignored_off_bulk(workspace):
    """cellsnp-lite already fixes the panel alleles, so the flag never reaches it."""
    proc = dryrun(
        workspace,
        workspace["sc_json"],
        "S1",
        "single_cell_genotyping",
        ["scRNA"],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "--constrain alleles" not in proc.stdout
