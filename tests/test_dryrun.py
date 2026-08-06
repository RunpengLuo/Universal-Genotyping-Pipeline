#!/usr/bin/env python3
"""Dry-run the DAG for every workflow mode, from both sample-file formats.

Runpeng Luo (2026-07-12)

These tests build the DAG only (`snakemake -n`); no rule is executed and no real
data is needed. They cover sample-file parsing and validation, rule wiring, the
storage() wrapping of remote inputs, and the final targets of each mode. A JSON
sample file and the equivalent TSV must yield the same DAG.

Executing the rules on real data is not covered; see docs/TODO.md.

Dependencies:
  pytest; snakemake on PATH.

Usage:
  pytest tests/                 # dry-run tests (run from the repo root)
  pytest tests/ -m network      # also hit the GIAB URLs

Notes/References:
  Snakemake's own per-rule test generator (`--generate-unit-tests`) needs a prior
  successful run, so it is only usable once real fixtures exist: docs/TODO.md.
"""

import json
import os

import pytest

from conftest import dryrun, job_counts

# (id, sample-file key, sample_id, workflow_mode, assay_types)
BULK_CASES = [
    ("bulk", "bulk", "T1", "bulk_genotyping", ["bulkWGS"]),
]
SC_CASES = [
    ("multiome", "sc", "S1", "single_cell_genotyping", ["scRNA", "scATAC"]),
    ("visium", "sc", "V1", "single_cell_genotyping", ["VISIUM"]),
    ("copytyping", "sc", "S1", "copytyping_preprocess", ["scATAC"]),
]
ALL_CASES = BULK_CASES + SC_CASES


@pytest.mark.parametrize("fmt", ["json", "tsv"])
@pytest.mark.parametrize(
    "case_id,sheet,sample_id,mode,assays",
    ALL_CASES,
    ids=[c[0] for c in ALL_CASES],
)
def test_dag_builds(workspace, fmt, case_id, sheet, sample_id, mode, assays):
    """Every mode builds a non-empty DAG from either sample-file format."""
    proc = dryrun(workspace, workspace[f"{sheet}_{fmt}"], sample_id, mode, assays)
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    assert counts, f"no jobs planned:\n{proc.stdout[-2000:]}"
    assert "all" in counts


@pytest.mark.parametrize(
    "case_id,sheet,sample_id,mode,assays",
    ALL_CASES,
    ids=[c[0] for c in ALL_CASES],
)
def test_json_and_tsv_agree(workspace, case_id, sheet, sample_id, mode, assays):
    """A TSV sheet plans exactly the same jobs as the equivalent JSON."""
    a = dryrun(workspace, workspace[f"{sheet}_json"], sample_id, mode, assays)
    b = dryrun(workspace, workspace[f"{sheet}_tsv"], sample_id, mode, assays)
    assert a.returncode == 0 and b.returncode == 0
    assert job_counts(a.stdout) == job_counts(b.stdout)


def test_bulk_rules(workspace):
    """Bulk genotyping plans the genotype/phase/pileup/depth chain."""
    proc = dryrun(
        workspace, workspace["bulk_json"], "T1", "bulk_genotyping", ["bulkWGS"]
    )
    counts = job_counts(proc.stdout)
    assert {"genotype_snps_bulk", "phase_snps_eagle", "combine_counts"} <= set(counts)
    # one bcftools pileup and one mosdepth per dataset (normal + tumor)
    assert counts["pileup_snps_bulk_bcftools"] == 2
    assert counts["run_mosdepth"] == 2


def test_bulk_stream_mode(workspace):
    """remote_mode=stream: bcftools rules carry -r and mosdepth is per-chrom + merged."""
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=("remote_mode=stream",),
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    # depth becomes per-chrom mosdepth + a merge; the whole-file rule is gone
    assert "run_mosdepth_chrom" in counts and "merge_mosdepth" in counts
    assert "run_mosdepth" not in counts
    # genotype/pileup restrict to the config chroms via index jumps
    assert "-r chr22" in proc.stdout


def test_stream_rejected_for_single_cell(workspace):
    """remote_mode=stream errors for non-bulk modes (cellsnp-lite cannot read URLs)."""
    proc = dryrun(
        workspace,
        workspace["sc_json"],
        "S1",
        "single_cell_genotyping",
        ["scRNA", "scATAC"],
        extra=("remote_mode=stream",),
    )
    assert proc.returncode != 0
    assert "only supported for bulk_genotyping" in (proc.stdout + proc.stderr)


def test_breakpoint_presegmentation(workspace):
    """build_segment_bed + per-stream window build always run for bulk; bedpe feeds the segment BED."""
    base = dryrun(
        workspace, workspace["bulk_json"], "T1", "bulk_genotyping", ["bulkWGS"]
    )
    bedpe = dryrun(
        workspace, workspace["bulk_bedpe_json"], "B1", "bulk_genotyping", ["bulkWGS"]
    )
    assert base.returncode == 0 and bedpe.returncode == 0, bedpe.stderr[-1500:]
    base_counts, bedpe_counts = job_counts(base.stdout), job_counts(bedpe.stdout)
    # the segment BED + per-stream window build are always-on for bulk (both runs)
    for rule in (
        "build_segment_bed",
        "build_window_bed",
    ):
        assert rule in base_counts, f"{rule} missing (base):\n{base.stdout[-2000:]}"
        assert rule in bedpe_counts, f"{rule} missing (bedpe):\n{bedpe.stdout[-2000:]}"
    # a breakpoint_bedpe only feeds build_segment_bed when present
    assert "sv.bedpe" in bedpe.stdout
    assert "sv.bedpe" not in base.stdout


def test_mixed_wgs_wes(workspace):
    """bulkWGS + bulkWES share one window grid and one joint bulk bb dir (WES == WGS)."""
    proc = dryrun(
        workspace,
        workspace["bulk_mixed_json"],
        "MX",
        "bulk_genotyping",
        ["bulkWGS", "bulkWES"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    # segment BED + exactly one shared window BED (no per-stream / wes_targets)
    assert "build_segment_bed" in counts
    assert counts.get("build_window_bed", 0) == 1, proc.stdout[-2000:]
    assert "/windows.bed.gz" in proc.stdout
    assert "wes_targets" not in proc.stdout
    assert "wgs_windows.bed.gz" not in proc.stdout
    assert "wes_windows.bed.gz" not in proc.stdout
    # one joint binning into a single bb/bulk dir (no per-stream subdir)
    assert "combine_counts" in counts
    assert "/bulk/bb.tsv.gz" in proc.stdout
    assert "bulkWES/bb.tsv.gz" not in proc.stdout


def test_prebuilt_windows_skip_build(workspace):
    """A pre-built window_bed (WGS-only, no BEDPE) is consumed directly, nothing built."""
    ref = workspace["ref"]
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=[f"window_bed={ref}/window.bed"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    # build_window_bed is not planned; the prebuilt window_bed is read directly
    assert "build_window_bed" not in counts
    assert "build_segment_bed" in counts
    assert f"{ref}/window.bed" in proc.stdout


def test_prebuilt_windows_ignored_with_bedpe(workspace):
    """A BEDPE re-tiles the arms, so a pre-built window_bed is ignored and rebuilt."""
    ref = workspace["ref"]
    proc = dryrun(
        workspace,
        workspace["bulk_bedpe_json"],
        "B1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=[f"window_bed={ref}/window.bed"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    assert "build_window_bed" in counts
    assert "subset_prebuilt_window_bed" not in counts


def test_prebuilt_windows_cover_all_assays(workspace):
    """A pre-built window_bed serves WGS + WES alike, so nothing is built."""
    ref = workspace["ref"]
    proc = dryrun(
        workspace,
        workspace["bulk_mixed_json"],
        "MX",
        "bulk_genotyping",
        ["bulkWGS", "bulkWES"],
        extra=[f"window_bed={ref}/window.bed"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    assert "build_window_bed" not in counts
    assert f"{ref}/window.bed" in proc.stdout


def test_prebuilt_windows_skip_repliseq(workspace):
    """A pre-built window_bed skips the Repli-seq fetch (do_repliseq active on hg38).

    Network-free: the Repli-seq rules are the only URL-storage inputs here, and they
    are gated inside `if not use_prebuilt_windows`, so nothing queries a remote host.
    """
    ref = workspace["ref"]
    sheet = _sheet_with_refvers(workspace, "repliseq.json", ["hg38", "hg38"])
    proc = dryrun(
        workspace,
        sheet,
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=["reference_version=hg38", f"window_bed={ref}/window.bed"],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    counts = job_counts(proc.stdout)
    assert "repliseq_bigwig_to_bedgraph" not in counts
    assert "build_window_bed" not in counts


def test_scatac_fragments_are_tracked(workspace):
    """The ATAC fragments file is a tracked input, not resolved inside a script."""
    proc = dryrun(
        workspace,
        workspace["sc_json"],
        "S1",
        "single_cell_genotyping",
        ["scRNA", "scATAC"],
    )
    assert proc.returncode == 0
    assert "atac_fragments.tsv.gz" in proc.stdout


def test_visium_spatial_files_are_tracked(workspace):
    """Each Space Ranger spatial/ file is a tracked input of process_rna_anndata."""
    proc = dryrun(
        workspace, workspace["sc_json"], "V1", "single_cell_genotyping", ["VISIUM"]
    )
    assert proc.returncode == 0
    for name in (
        "tissue_positions.csv",
        "scalefactors_json.json",
        "tissue_hires_image.png",
        "tissue_lowres_image.png",
    ):
        assert name in proc.stdout, f"{name} is not a tracked input"


def test_null_files_value_is_dropped(workspace):
    """A null optional input is absent, not the literal path "None"."""
    sheet = os.path.join(workspace["root"], "null_file.json")
    doc = json.loads(open(workspace["bulk_json"]).read())
    doc["samples"][1]["files"]["breakpoint_bedpe"] = None
    with open(sheet, "w") as fh:
        json.dump(doc, fh)
    proc = dryrun(workspace, sheet, "T1", "bulk_genotyping", ["bulkWGS"])
    base = dryrun(
        workspace, workspace["bulk_json"], "T1", "bulk_genotyping", ["bulkWGS"]
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert job_counts(proc.stdout) == job_counts(base.stdout)


def test_tsv_missing_required_column_fails(workspace):
    """TSV columns are the record keys; a missing one is named, not silently empty."""
    sheet = os.path.join(workspace["root"], "no_col.tsv")
    lines = open(workspace["bulk_tsv"]).read().splitlines()
    header = lines[0].split("\t")
    drop = header.index("reference_version")
    rows = [
        "\t".join(c for i, c in enumerate(ln.split("\t")) if i != drop) for ln in lines
    ]
    with open(sheet, "w") as fh:
        fh.write("\n".join(rows) + "\n")
    proc = dryrun(workspace, sheet, "T1", "bulk_genotyping", ["bulkWGS"])
    assert proc.returncode != 0
    out = proc.stdout + proc.stderr
    assert "missing required column(s)" in out and "reference_version" in out


def test_record_without_reference_version_fails(workspace):
    """reference_version is required on every record, not optional provenance."""
    sheet = os.path.join(workspace["root"], "no_refver.json")
    doc = json.loads(open(workspace["bulk_json"]).read())
    del doc["samples"][1]["reference_version"]
    with open(sheet, "w") as fh:
        json.dump(doc, fh)
    proc = dryrun(workspace, sheet, "T1", "bulk_genotyping", ["bulkWGS"])
    assert proc.returncode != 0
    out = proc.stdout + proc.stderr
    assert "missing required key(s)" in out and "reference_version" in out


def _sheet_with_refvers(workspace, name, refvers):
    """Copy bulk_json, setting each record's reference_version from *refvers*."""
    sheet = os.path.join(workspace["root"], name)
    doc = json.loads(open(workspace["bulk_json"]).read())
    for rec, refver in zip(doc["samples"], refvers):
        rec["reference_version"] = refver
    with open(sheet, "w") as fh:
        json.dump(doc, fh)
    return sheet


def test_chromosome_absent_from_genome_size_fails(workspace):
    """Config chromosomes are checked against the genome once, at DAG build."""
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=["chromosomes=[22,99]"],
    )
    assert proc.returncode != 0
    out = proc.stdout + proc.stderr
    assert "have no contig in" in out and "'99'" in out


def test_config_species_required(workspace):
    """An unset config species is an error, not a silent fallback to human."""
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=["species="],
    )
    assert proc.returncode != 0
    assert "species is required" in proc.stdout + proc.stderr


def test_species_reaches_parse_genetic_map(workspace):
    """The shipped config says human, and that is what parse_genetic_map is given."""
    proc = dryrun(
        workspace, workspace["bulk_json"], "T1", "bulk_genotyping", ["bulkWGS"]
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "parse_genetic_map" in job_counts(proc.stdout)
    mouse = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=["species=mouse"],
    )
    assert mouse.returncode == 0, mouse.stderr[-1500:]


def test_unsupported_species_warns(workspace):
    """An unknown species warns but runs; it only matters if a gmap needs relabeling."""
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=["species=zebrafish"],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "species='zebrafish' is not natively supported" in proc.stdout


def test_config_reference_version_required(workspace):
    """An unset config reference_version is an error, not a silent no-filter."""
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=["reference_version="],
    )
    assert proc.returncode != 0
    assert "reference_version is required" in proc.stdout + proc.stderr


def test_record_reference_version_alias_matches(workspace):
    """Config and records may spell the build differently; both are canonicalized."""
    sheet = _sheet_with_refvers(workspace, "alias.json", ["T2T-CHM13v2.0", "chm13v2.0"])
    proc = dryrun(
        workspace,
        sheet,
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=["reference_version=CHM13"],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]


def test_records_of_another_build_are_dropped(workspace):
    """A record on a different build is not selected, and the error names the builds."""
    sheet = _sheet_with_refvers(workspace, "otherbuild.json", ["hg19", "hg19"])
    proc = dryrun(workspace, sheet, "T1", "bulk_genotyping", ["bulkWGS"])
    assert proc.returncode != 0
    out = proc.stdout + proc.stderr
    assert "reference_version='chm13v2'" in out
    assert "hg19 (2)" in out


def test_unrecognized_reference_version_still_selects(workspace):
    """An unsupported build warns but still runs, so long as records match it."""
    sheet = _sheet_with_refvers(
        workspace, "giabv3.json", ["GRCh38-GIABv3", "GRCh38-GIABv3"]
    )
    proc = dryrun(
        workspace,
        sheet,
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=["reference_version=GRCh38-GIABv3"],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "is not natively supported" in proc.stdout
    assert "grch38-giabv3" in proc.stdout


def test_unknown_sample_id_fails(workspace):
    """An unmatched sample_id is an error, not an empty DAG."""
    proc = dryrun(
        workspace, workspace["bulk_json"], "nope", "bulk_genotyping", ["bulkWGS"]
    )
    assert proc.returncode != 0
    assert "no records for sample_id" in proc.stdout + proc.stderr


@pytest.mark.parametrize("mode", ["bulk_genotyping", "single_cell_genotyping"])
def test_het_snp_vcf_skips_genotyping(workspace, mode):
    """A supplied het_snp_vcf replaces genotyping in every mode."""
    sheet, sample_id, assays = (
        ("bulk_json", "T1", ["bulkWGS"])
        if mode == "bulk_genotyping"
        else ("sc_json", "S1", ["scRNA", "scATAC"])
    )
    ref = workspace["ref"]
    base = dryrun(workspace, workspace[sheet], sample_id, mode, assays)
    with_vcf = dryrun(
        workspace,
        workspace[sheet],
        sample_id,
        mode,
        assays,
        extra=[f"het_snp_vcf={ref}/het_snps.vcf.gz", "het_snp_vcf_phased=False"],
    )
    assert base.returncode == 0 and with_vcf.returncode == 0, with_vcf.stderr[-1500:]
    genotype_rules = {
        r for r in job_counts(base.stdout) if r.startswith("genotype_snps")
    }
    assert genotype_rules, "baseline should genotype"
    assert not genotype_rules & set(job_counts(with_vcf.stdout))
    # het_snp_vcf_phased=false: split per chromosome, then phased
    counts = job_counts(with_vcf.stdout)
    assert "split_het_snp_vcf" in counts
    assert any(r.startswith("phase_snps") for r in counts)


def test_phased_het_snp_vcf_skips_phasing(workspace):
    """het_snp_vcf_phased=true replaces genotyping and phasing."""
    ref = workspace["ref"]
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=[f"het_snp_vcf={ref}/het_snps.vcf.gz", "het_snp_vcf_phased=True"],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    counts = job_counts(proc.stdout)
    for rule in counts:
        assert not rule.startswith(("genotype_snps", "phase_snps", "split_het_snp_vcf"))
    assert "concat_and_extract_phased_het_snps" not in counts
    assert "pileup_snps_bulk_bcftools" in counts


def test_copytyping_requires_phased_vcf(workspace):
    """copytyping_preprocess never phases, so its het_snp_vcf must be declared phased."""
    ref = workspace["ref"]
    proc = dryrun(
        workspace,
        workspace["sc_json"],
        "S1",
        "copytyping_preprocess",
        ["scATAC"],
        extra=[
            f"het_snp_vcf={ref}/het_snps.vcf.gz",
            "het_snp_vcf_phased=False",
            f"bb_file={ref}/bb.tsv.gz",
        ],
    )
    assert proc.returncode != 0
    assert "het_snp_vcf must be phased" in proc.stdout + proc.stderr
