#!/usr/bin/env python3
"""Dry-run the DAG for every workflow mode, from both sample-file formats.

Last update: 2026-08-06

Covers:
- modes: each of the three builds a non-empty DAG
- formats: a TSV sheet plans exactly the same jobs as JSON
- validation: bad sample sheets and configs fail at DAG build
- wiring: window build, repliseq, streaming, and the short-circuit VCFs
Notes:
- cost: no rule runs and no real data is needed
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


@pytest.mark.parametrize(
    "sheet,sample_id,mode,assays",
    [
        ("bulk_json", "T1", "bulk_genotyping", ["bulkWGS"]),
        ("sc_json", "S1", "single_cell_genotyping", ["scRNA", "scATAC"]),
    ],
    ids=["bulk", "single_cell"],
)
def test_windows_are_built_from_the_segments(workspace, sheet, sample_id, mode, assays):
    """With no window_bed, every mode builds the segment BED and tiles it."""
    ref = workspace["ref"]
    proc = dryrun(workspace, workspace[sheet], sample_id, mode, assays)
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    assert "build_segment_bed" in counts
    assert "build_window_bed" in counts
    # the configured segmentation feeds build_segment_bed, which feeds the tiling
    assert f"{ref}/segment.bed" in proc.stdout
    assert "/windows.bed.gz" in proc.stdout


@pytest.mark.parametrize(
    "mode,assays",
    [
        ("bulk_genotyping", ["bulkWGS"]),
        ("single_cell_genotyping", ["scRNA", "scATAC"]),
        ("copytyping_preprocess", ["scATAC"]),
    ],
)
def test_segment_bed_defaults_to_region_bed(workspace, mode, assays):
    """An unset segment_bed falls back to region_bed: one segment per arm."""
    sheet = (
        workspace["bulk_json"] if mode == "bulk_genotyping" else workspace["sc_json"]
    )
    sample_id = "T1" if mode == "bulk_genotyping" else "S1"
    proc = dryrun(workspace, sheet, sample_id, mode, assays, extra=["segment_bed="])
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "segment_bed unset, one segment per region_bed arm" in proc.stdout
    assert job_counts(proc.stdout)["build_segment_bed"] == 1
    assert f"{workspace['ref']}/region.bed" in proc.stdout


@pytest.mark.parametrize(
    "mode,assays",
    [
        ("bulk_genotyping", ["bulkWGS"]),
        ("single_cell_genotyping", ["scRNA", "scATAC"]),
        ("copytyping_preprocess", ["scATAC"]),
    ],
)
def test_segment_bed_is_built_in_every_mode(workspace, mode, assays):
    """Every mode reads aux/segment.bed, so build_segment_bed is always planned."""
    sheet = (
        workspace["bulk_json"] if mode == "bulk_genotyping" else workspace["sc_json"]
    )
    sample_id = "T1" if mode == "bulk_genotyping" else "S1"
    proc = dryrun(workspace, sheet, sample_id, mode, assays)
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "build_segment_bed" in job_counts(proc.stdout)


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


@pytest.mark.parametrize(
    "sheet,sample_id,mode,assays",
    [
        ("bulk_json", "T1", "bulk_genotyping", ["bulkWGS"]),
        ("sc_json", "S1", "single_cell_genotyping", ["scRNA", "scATAC"]),
    ],
    ids=["bulk", "single_cell"],
)
def test_prebuilt_windows_are_used_not_built(workspace, sheet, sample_id, mode, assays):
    """A pre-built window_bed is consumed as-is; nothing is tiled."""
    ref = workspace["ref"]
    proc = dryrun(
        workspace,
        workspace[sheet],
        sample_id,
        mode,
        assays,
        extra=[f"window_bed={ref}/window.bed"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    # nothing is tiled; the supplied grid is read straight from its configured path
    assert "build_window_bed" not in counts
    # build_segment_bed still runs: the SNP side and the RD QC overlay both read it
    assert "build_segment_bed" in counts
    assert f"{ref}/window.bed" in proc.stdout


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
    are gated inside `if build_windows`, so nothing queries a remote host.
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


def test_unread_files_key_is_dropped(workspace):
    """A files key the assay never reads is dropped, null or not, and plans the same DAG."""
    sheet = os.path.join(workspace["root"], "null_file.json")
    doc = json.loads(open(workspace["bulk_json"]).read())
    doc["samples"][1]["files"]["fragments"] = (
        None  # scATAC-only key on a bulkWGS record
    )
    doc["samples"][0]["files"]["barcodes"] = "/path/to/nowhere.tsv.gz"
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
    assert "missing required key(s)" in out and "reference_version" in out


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


def _sheet_with_refvers(workspace, name, refvers, source="bulk_json"):
    """Copy a sample sheet, setting each record's reference_version from *refvers*."""
    sheet = os.path.join(workspace["root"], name)
    doc = json.loads(open(workspace[source]).read())
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
    assert "are not found in" in out and "'99'" in out


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
    """A record on a different build is not selected, leaving nothing to run."""
    sheet = _sheet_with_refvers(workspace, "otherbuild.json", ["hg19", "hg19"])
    proc = dryrun(workspace, sheet, "T1", "bulk_genotyping", ["bulkWGS"])
    assert proc.returncode != 0
    out = proc.stdout + proc.stderr
    assert "no datasets exist after selection" in out


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


def test_unknown_sample_id_fails(workspace):
    """An unmatched sample_id is an error, not an empty DAG."""
    proc = dryrun(
        workspace, workspace["bulk_json"], "nope", "bulk_genotyping", ["bulkWGS"]
    )
    assert proc.returncode != 0
    assert "no datasets exist after selection" in proc.stdout + proc.stderr


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


def _longread_sheet(workspace, name):
    """Bulk sheet with a short-read normal, a long-read normal and a long-read tumor."""
    ref = workspace["ref"]

    def files(stem):
        return {
            "alignment": f"{ref}/{stem}.bam",
            "alignment_index": f"{ref}/{stem}.bam.bai",
        }

    doc = {
        "version": 1,
        "samples": [
            {
                "sample_id": "T1",
                "dataset_id": "N1",
                "assay_type": "bulkWGS",
                "sample_type": "normal",
                "reference_version": "chm13v2",
                "files": files("normal"),
            },
            {
                "sample_id": "T1",
                "dataset_id": "L1",
                "assay_type": "bulkWGS-lr",
                "sample_type": "normal",
                "reference_version": "chm13v2",
                "files": files("normal"),
            },
            {
                "sample_id": "T1",
                "dataset_id": "D1",
                "assay_type": "bulkWGS-lr",
                "sample_type": "tumor",
                "reference_version": "chm13v2",
                "files": files("tumor"),
            },
        ],
    }
    sheet = os.path.join(workspace["root"], name)
    with open(sheet, "w") as fh:
        json.dump(doc, fh)
    return sheet


def test_longphase_auto_picks_the_long_read_normal(workspace):
    """phaser=longphase with no phase_dataset_ids co-phases the long-read normal."""
    sheet = _longread_sheet(workspace, "longread.json")
    proc = dryrun(
        workspace,
        sheet,
        "T1",
        "bulk_genotyping",
        ["bulkWGS", "bulkWGS-lr"],
        extra=["phaser=longphase"],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "phase_dataset_ids: ['L1']" in proc.stdout
    counts = job_counts(proc.stdout)
    assert "phase_snps_longphase" in counts
    assert "phase_snps_eagle" not in counts and "parse_genetic_map" not in counts


def test_phase_dataset_ids_must_be_long_read(workspace):
    """A short-read dataset named in phase_dataset_ids is rejected, not silently kept."""
    sheet = _longread_sheet(workspace, "longread_shortread.json")
    proc = dryrun(
        workspace,
        sheet,
        "T1",
        "bulk_genotyping",
        ["bulkWGS", "bulkWGS-lr"],
        extra=["phaser=longphase", 'phase_dataset_ids=["N1"]'],
    )
    assert proc.returncode != 0
    assert "needs long reads" in proc.stdout + proc.stderr


def test_phase_dataset_ids_outside_selection_fails(workspace):
    """An id in the sheet but dropped by assay_types is an error, not an auto-pick."""
    sheet = _longread_sheet(workspace, "longread_selection.json")
    proc = dryrun(
        workspace,
        sheet,
        "T1",
        "bulk_genotyping",
        ["bulkWGS-lr"],
        extra=["phaser=longphase", 'phase_dataset_ids=["N1"]'],
    )
    assert proc.returncode != 0
    assert "phase_dataset_ids not found in the records" in proc.stdout + proc.stderr


def test_genotype_dataset_ids_outside_selection_fails(workspace):
    """bulkWES E1 is in the sheet but dropped by assay_types, so naming it fails."""
    proc = dryrun(
        workspace,
        workspace["bulk_mixed_json"],
        "MX",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=['genotype_dataset_ids=["E1"]'],
    )
    assert proc.returncode != 0
    assert "genotype_dataset_ids not found in the records" in proc.stdout + proc.stderr


def test_genotype_dataset_ids_in_selection_is_used(workspace):
    """The same id check accepts a dataset the selection kept."""
    proc = dryrun(
        workspace,
        workspace["bulk_mixed_json"],
        "MX",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=['genotype_dataset_ids=["N1"]'],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "genotype_dataset_ids: ['N1']" in proc.stdout


def test_single_cell_skips_repliseq(workspace):
    """`do_repliseq` is the one thing build_window_bed still gates on the mode.

    GC and MAP are annotated identically in every mode, but the ENCODE bigWig fetch +
    liftOver must stay out of a single-cell DAG: they are 15 UCSC downloads feeding a
    covariate only rd_correct (bulk) reads.
    """
    sheet = _sheet_with_refvers(
        workspace, "sc_hg38.json", ["hg38"] * 4, source="sc_json"
    )
    proc = dryrun(
        workspace,
        sheet,
        "S1",
        "single_cell_genotyping",
        ["scRNA", "scATAC"],
        extra=["reference_version=hg38"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    assert "build_window_bed" in counts
    assert "repliseq_bigwig_to_bedgraph" not in counts
    assert "repliseq_liftover" not in counts


def test_copytyping_does_not_build_windows(workspace):
    """copytyping_preprocess bins onto its own bb_file, so no window grid is needed."""
    proc = dryrun(
        workspace, workspace["sc_json"], "S1", "copytyping_preprocess", ["scATAC"]
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    assert "build_window_bed" not in counts


def test_single_cell_binning_reads_the_window_grid(workspace):
    """combine_counts_nonbulk takes the window BED as its fixed bins."""
    proc = dryrun(
        workspace,
        workspace["sc_json"],
        "S1",
        "single_cell_genotyping",
        ["scRNA", "scATAC"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "combine_counts_nonbulk" in job_counts(proc.stdout)
    assert "/windows.bed.gz" in proc.stdout
