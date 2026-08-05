"""Shared fixtures: a stub workspace the DAG can be built against.

Runpeng Luo (2026-07-12)

Dry-run tests only build the DAG, so every reference asset and input file can be
an empty stub; nothing is read. Sample files are written in both formats from one
description, so the JSON and TSV paths are compared on identical data.

Dependencies:
  pytest; snakemake on PATH.

Inputs
  none: every fixture is generated under pytest's tmp_path_factory
Outputs:
  workspace: paths of the stub reference assets and sample files
Notes/References:
  Sample-file format: docs/sample_sheet.md
"""

import gzip
import json
import os
import subprocess

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAKEFILE = os.path.join(REPO, "workflow", "Snakefile")
CONFIGFILE = os.path.join(REPO, "config", "config.yaml")

# empty stubs; the DAG only needs these paths to exist
REF_FILES = (
    "genome.fa",
    "genome_size.txt",
    "region.bed",
    "window.bed",
    "genes.gtf",
    "gmap.txt.gz",
    "snp_panel.vcf.gz",
    "bb.tsv.gz",
)

VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr22,length=100000>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
)
RANGER_FILES = (
    "gex_possorted_bam.bam",
    "gex_possorted_bam.bam.bai",
    "atac_possorted_bam.bam",
    "atac_possorted_bam.bam.bai",
    "possorted_genome_bam.bam",
    "possorted_genome_bam.bam.bai",
    "atac_fragments.tsv.gz",
    "filtered_feature_bc_matrix.h5",
    "filtered_feature_bc_matrix/barcodes.tsv.gz",
    "spatial/tissue_positions.csv",
    "spatial/scalefactors_json.json",
    "spatial/tissue_hires_image.png",
    "spatial/tissue_lowres_image.png",
)


def _touch(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "a").close()
    return path


@pytest.fixture(scope="session")
def workspace(tmp_path_factory):
    """Stub reference assets plus one bulk and one single-cell sample file per format."""
    root = tmp_path_factory.mktemp("ws")
    ref = root / "ref"
    for name in REF_FILES:
        _touch(str(ref / name))

    # parse_workflow resolves config chromosomes against this at DAG build, so it
    # needs real contigs (the other reference assets stay empty stubs)
    (ref / "genome_size.txt").write_text(
        "".join(f"chr{c}\t50818468\n" for c in list(range(1, 23)) + ["X", "Y"])
    )

    rows = "".join(
        f"chr22\t{pos}\t.\tA\tG\t60\tPASS\t.\tGT\t0|1\n" for pos in (1000, 2000, 3000)
    )
    with gzip.open(str(ref / "het_snps.vcf.gz"), "wt") as fh:
        fh.write(VCF_HEADER + rows)
    for chrom in ("22",):
        _touch(str(ref / "targets" / f"target.chr{chrom}.pos.gz"))
        _touch(str(ref / "panel" / f"chr{chrom}.genotypes.bcf"))

    outs = root / "outs"
    for name in RANGER_FILES:
        _touch(str(outs / name))
    for name in ("normal.bam", "tumor.bam"):
        _touch(str(ref / name))
        _touch(str(ref / f"{name}.bai"))

    # SV breakpoints BEDPE (0-based, like BED): two junctions on chr22
    (ref / "sv.bedpe").write_text(
        "chr22\t16000000\t16000001\tchr22\t16500000\t16500001\tsv1\t60\t+\t-\n"
        "chr22\t20000000\t20000001\tchr22\t30000000\t30000001\tsv2\t42\t-\t+\n"
    )
    barcodes = str(outs / "filtered_feature_bc_matrix" / "barcodes.tsv.gz")
    bulk_json = {
        "version": 1,
        "samples": [
            {
                "sample_id": "T1",
                "dataset_id": "N1",
                "assay_type": "bulkWGS",
                "sample_type": "normal",
                "reference_version": "chm13v2",
                "files": {
                    "alignment": str(ref / "normal.bam"),
                    "alignment_index": str(ref / "normal.bam.bai"),
                },
            },
            {
                "sample_id": "T1",
                "dataset_id": "D1",
                "rdr_base_dataset_id": "N1",
                "assay_type": "bulkWGS",
                "sample_type": "tumor",
                "reference_version": "chm13v2",
                "files": {
                    "alignment": str(ref / "tumor.bam"),
                    "alignment_index": str(ref / "tumor.bam.bai"),
                },
            },
        ],
    }
    # bulk with an SV BEDPE on the tumor -> breakpoint-aware pre-segmentation
    bulk_bedpe_json = {
        "version": 1,
        "samples": [
            {
                "sample_id": "B1",
                "dataset_id": "N1",
                "assay_type": "bulkWGS",
                "sample_type": "normal",
                "reference_version": "chm13v2",
                "files": {
                    "alignment": str(ref / "normal.bam"),
                    "alignment_index": str(ref / "normal.bam.bai"),
                },
            },
            {
                "sample_id": "B1",
                "dataset_id": "D1",
                "assay_type": "bulkWGS",
                "sample_type": "tumor",
                "reference_version": "chm13v2",
                "files": {
                    "alignment": str(ref / "tumor.bam"),
                    "alignment_index": str(ref / "tumor.bam.bai"),
                    "breakpoint_bedpe": str(ref / "sv.bedpe"),
                },
            },
        ],
    }
    # WGS + WES on one individual: joint SNPs/phasing, one shared window grid
    bulk_mixed_json = {
        "version": 1,
        "samples": [
            {
                "sample_id": "MX",
                "dataset_id": "N1",
                "assay_type": "bulkWGS",
                "sample_type": "normal",
                "reference_version": "chm13v2",
                "files": {
                    "alignment": str(ref / "normal.bam"),
                    "alignment_index": str(ref / "normal.bam.bai"),
                },
            },
            {
                "sample_id": "MX",
                "dataset_id": "D1",
                "rdr_base_dataset_id": "N1",
                "assay_type": "bulkWGS",
                "sample_type": "tumor",
                "reference_version": "chm13v2",
                "files": {
                    "alignment": str(ref / "tumor.bam"),
                    "alignment_index": str(ref / "tumor.bam.bai"),
                },
            },
            {
                "sample_id": "MX",
                "dataset_id": "E1",
                "assay_type": "bulkWES",
                "sample_type": "tumor",
                "reference_version": "chm13v2",
                "files": {
                    "alignment": str(ref / "tumor.bam"),
                    "alignment_index": str(ref / "tumor.bam.bai"),
                },
            },
        ],
    }
    sc_json = {
        "version": 1,
        "samples": [
            {
                "sample_id": "S1",
                "dataset_id": "U1",
                "assay_type": "scRNA",
                "sample_type": "tumor",
                "reference_version": "chm13v2",
                "files": {
                    "alignment": str(outs / "gex_possorted_bam.bam"),
                    "alignment_index": str(outs / "gex_possorted_bam.bam.bai"),
                    "barcodes": barcodes,
                    "matrix_h5": str(outs / "filtered_feature_bc_matrix.h5"),
                },
            },
            {
                "sample_id": "S1",
                "dataset_id": "U1",
                "assay_type": "scATAC",
                "sample_type": "tumor",
                "reference_version": "chm13v2",
                "files": {
                    "alignment": str(outs / "atac_possorted_bam.bam"),
                    "alignment_index": str(outs / "atac_possorted_bam.bam.bai"),
                    "barcodes": barcodes,
                    "fragments": str(outs / "atac_fragments.tsv.gz"),
                },
            },
            {
                "sample_id": "V1",
                "dataset_id": "W1",
                "assay_type": "VISIUM",
                "sample_type": "tumor",
                "reference_version": "chm13v2",
                "files": {
                    "alignment": str(outs / "possorted_genome_bam.bam"),
                    "alignment_index": str(outs / "possorted_genome_bam.bam.bai"),
                    "barcodes": barcodes,
                    "matrix_h5": str(outs / "filtered_feature_bc_matrix.h5"),
                    "tissue_positions": str(outs / "spatial" / "tissue_positions.csv"),
                    "scalefactors": str(outs / "spatial" / "scalefactors_json.json"),
                    "image_hires": str(outs / "spatial" / "tissue_hires_image.png"),
                    "image_lowres": str(outs / "spatial" / "tissue_lowres_image.png"),
                },
            },
        ],
    }

    # the TSV sheets encode the same records as the JSON above, one files.<key> column
    # per input; the DAGs must come out identical (test_json_and_tsv_agree)
    bulk_tsv = [
        "sample_id\tdataset_id\trdr_base_dataset_id\tassay_type\tsample_type"
        "\treference_version\tfiles.alignment\tfiles.alignment_index",
        f"T1\tN1\t\tbulkWGS\tnormal\tchm13v2\t{ref / 'normal.bam'}\t{ref / 'normal.bam.bai'}",
        f"T1\tD1\tN1\tbulkWGS\ttumor\tchm13v2\t{ref / 'tumor.bam'}\t{ref / 'tumor.bam.bai'}",
    ]
    sc_tsv = [
        "sample_id\tdataset_id\tassay_type\tsample_type\treference_version"
        "\tfiles.alignment\tfiles.alignment_index\tfiles.barcodes\tfiles.fragments"
        "\tfiles.matrix_h5\tfiles.tissue_positions\tfiles.scalefactors"
        "\tfiles.image_hires\tfiles.image_lowres",
        f"S1\tU1\tscRNA\ttumor\tchm13v2\t{outs / 'gex_possorted_bam.bam'}"
        f"\t{outs / 'gex_possorted_bam.bam.bai'}\t{barcodes}\t"
        f"\t{outs / 'filtered_feature_bc_matrix.h5'}\t\t\t\t",
        f"S1\tU1\tscATAC\ttumor\tchm13v2\t{outs / 'atac_possorted_bam.bam'}"
        f"\t{outs / 'atac_possorted_bam.bam.bai'}\t{barcodes}"
        f"\t{outs / 'atac_fragments.tsv.gz'}\t\t\t\t\t",
        f"V1\tW1\tVISIUM\ttumor\tchm13v2\t{outs / 'possorted_genome_bam.bam'}"
        f"\t{outs / 'possorted_genome_bam.bam.bai'}\t{barcodes}\t"
        f"\t{outs / 'filtered_feature_bc_matrix.h5'}"
        f"\t{outs / 'spatial' / 'tissue_positions.csv'}"
        f"\t{outs / 'spatial' / 'scalefactors_json.json'}"
        f"\t{outs / 'spatial' / 'tissue_hires_image.png'}"
        f"\t{outs / 'spatial' / 'tissue_lowres_image.png'}",
    ]

    paths = {}
    for name, doc in (
        ("bulk", bulk_json),
        ("bulk_bedpe", bulk_bedpe_json),
        ("bulk_mixed", bulk_mixed_json),
        ("sc", sc_json),
    ):
        p = root / f"{name}.json"
        p.write_text(json.dumps(doc, indent=1))
        paths[f"{name}_json"] = str(p)
    for name, lines in (("bulk", bulk_tsv), ("sc", sc_tsv)):
        p = root / f"{name}.tsv"
        p.write_text("\n".join(lines) + "\n")
        paths[f"{name}_tsv"] = str(p)

    paths["root"] = str(root)
    paths["ref"] = str(ref)
    paths["outs"] = str(outs)
    return paths


def dryrun(workspace, sample_file, sample_id, workflow_mode, assay_types, extra=()):
    """Run `snakemake -n` against the stub workspace and return its CompletedProcess."""
    ref = workspace["ref"]
    out = os.path.join(workspace["root"], f"out_{workflow_mode}_{sample_id}")
    os.makedirs(out, exist_ok=True)
    cmd = [
        "snakemake",
        "-n",
        "-p",
        "--cores",
        "1",
        "-s",
        SNAKEFILE,
        "--configfile",
        CONFIGFILE,
        "--directory",
        out,
        "--config",
        f"sample_file={sample_file}",
        f"sample_id={sample_id}",
        f"workflow_mode={workflow_mode}",
        f"assay_types={json.dumps(assay_types)}",
        "chromosomes=[22]",
        # chm13v2 is outside REPLISEQ_REFVERS, so do_repliseq stays off and no rule
        # takes a URL input: the DAG builds offline. Must match the fixture records.
        "reference_version=chm13v2",
        f"reference={ref}/genome.fa",
        f"genome_size={ref}/genome_size.txt",
        f"region_bed={ref}/region.bed",
        f"gtf_file={ref}/genes.gtf",
        f"gmap_path={ref}/gmap.txt.gz",
        f"snp_panel={ref}/snp_panel.vcf.gz",
        f"snp_targets={ref}/targets",
        f"phasing_panel={ref}/panel",
        *extra,
    ]
    # het_snp_vcf now skips genotyping in every mode, so only pass it where required
    if workflow_mode == "copytyping_preprocess" and not any(
        e.startswith("het_snp_vcf=") for e in extra
    ):
        cmd += [
            f"het_snp_vcf={ref}/het_snps.vcf.gz",
            "het_snp_vcf_phased=True",
            f"bb_file={ref}/bb.tsv.gz",
        ]
    return subprocess.run(cmd, capture_output=True, text=True)


def job_counts(stdout):
    """Parse the `Job stats` table of a dry run into {rule: count}."""
    counts = {}
    in_table = False
    for line in stdout.splitlines():
        if line.startswith("Job stats:"):
            in_table = True
            continue
        if in_table:
            parts = line.split()
            if len(parts) == 2 and parts[1].isdigit():
                counts[parts[0]] = int(parts[1])
            elif line.startswith("total"):
                break
    counts.pop("job", None)
    return counts
