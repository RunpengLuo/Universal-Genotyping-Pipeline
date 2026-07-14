"""Assay-type, data-type, and reference constants shared by Snakefile and scripts."""

import os

NONBULK_ASSAYS = {"scATAC", "scRNA", "VISIUM", "VISIUM3prime"}
BULK_ASSAYS = {"bulkWGS", "bulkWGS-lr", "bulkWES"}
# Long-read bulk assays (identified by the "-lr" suffix convention).
LONGREAD_ASSAYS = {"bulkWGS-lr"}
ALLOWED_ASSAY_TYPES = list(BULK_ASSAYS) + list(NONBULK_ASSAYS)
SPATIAL_ASSAYS = {"VISIUM", "VISIUM3prime"}

# Phasers grouped by the evidence they consume. Panel-based phasers need a
# population reference panel + genetic map; the long-read phaser instead reads
# haplotype evidence directly from a long-read BAM.
LONGREAD_PHASER = {"longphase"}
PANEL_PHASER = {"eagle", "shapeit"}

# Reference genome versions that natively supported
REFVERS = ["hg19", "hg38", "chm13v2", "mm10"]

# Sex chromosomes -> their integer label in an Eagle map (X/Y follow the
# autosomes: human X=23,Y=24; mm10 X=20,Y=21). See parse_genetic_map.py.
REFVER2SEXCHROM = {
    "hg19": {"X": 23, "Y": 24},
    "hg38": {"X": 23, "Y": 24},
    "chm13v2": {"X": 23, "Y": 24},
    "mm10": {"X": 20, "Y": 21},
}

ASSAY_TYPE2MODALITY = {
    "bulkWGS": "DNA",
    "bulkWGS-lr": "DNA",
    "bulkWES": "DNA",
    "scATAC": "DNA",
    "scRNA": "RNA",
    "VISIUM": "RNA",
    "VISIUM3prime": "RNA",
}

ASSAY_TYPE2FEATURE_TYPE = {
    "bulkWGS": "dna",
    "bulkWGS-lr": "dna",
    "bulkWES": "exon",
    "scATAC": "tile",
    "scRNA": "gene",
    "VISIUM": "gene",
    "VISIUM3prime": "gene",
}

# 10x Cell/Space Ranger outs/ filenames; tuples list alternate spellings, newest first.
RANGER_MATRIX_H5 = ("filtered_feature_bc_matrix.h5",)
RANGER_BARCODES = "filtered_feature_bc_matrix/barcodes.tsv.gz"
RANGER_ATAC_FRAGMENTS = ("atac_fragments.tsv.gz", "fragments.tsv.gz")
RANGER_BAMS = {
    "scRNA": "gex_possorted_bam.bam",
    "scATAC": "atac_possorted_bam.bam",
    "VISIUM": "possorted_genome_bam.bam",
    "VISIUM3prime": "possorted_genome_bam.bam",
}
RANGER_SPATIAL_DIR = "spatial"
RANGER_TISSUE_POSITIONS = ("tissue_positions.csv", "tissue_positions_list.csv")
RANGER_SCALEFACTORS = ("scalefactors_json.json",)
RANGER_IMAGE_HIRES = ("tissue_hires_image.png",)
RANGER_IMAGE_LOWRES = ("tissue_lowres_image.png",)

# files key -> (Ranger name(s), under spatial/)
RANGER_LAYOUT = {
    "matrix_h5": (RANGER_MATRIX_H5, False),
    "fragments": (RANGER_ATAC_FRAGMENTS, False),
    "tissue_positions": (RANGER_TISSUE_POSITIONS, True),
    "scalefactors": (RANGER_SCALEFACTORS, True),
    "image_hires": (RANGER_IMAGE_HIRES, True),
    "image_lowres": (RANGER_IMAGE_LOWRES, True),
}

GTF_COLUMNS = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attributes",
]


WORKFLOW_MODES = ("bulk_genotyping", "single_cell_genotyping", "copytyping_preprocess")

URL_SCHEMES = ("http://", "https://", "ftp://", "s3://")

# sample-file schema; see docs/sample_sheet.md. Records may carry any other key;
# it is kept as-is and never read.
REQUIRED_RECORD_KEYS = ("sample_id", "dataset_id", "assay_type", "sample_type", "files")
# optional, read by a rule
OPTIONAL_RECORD_KEYS = ("rdr_base_dataset_id",)
# provenance-only; coerced to str, never read by a rule
PROVENANCE_KEYS = ("passage_id", "platform", "reference_version", "cancer_type")

ALIGNMENT_FILES = {"alignment", "alignment_index"}

# files read per assay type; the union is the whole `files` vocabulary, and any other
# key in a record's `files` is ignored. Non-alignment entries apply to single-cell modes only
REQUIRED_FILES = {
    "bulkWGS": ALIGNMENT_FILES,
    "bulkWGS-lr": ALIGNMENT_FILES,
    "bulkWES": ALIGNMENT_FILES,
    "scRNA": ALIGNMENT_FILES | {"barcodes", "matrix_h5"},
    "scATAC": ALIGNMENT_FILES | {"barcodes", "fragments"},
    "VISIUM": ALIGNMENT_FILES
    | {
        "barcodes",
        "matrix_h5",
        "tissue_positions",
        "scalefactors",
        "image_hires",
        "image_lowres",
    },
    # squidpy cannot load tissue images for 3' data
    "VISIUM3prime": ALIGNMENT_FILES
    | {"barcodes", "matrix_h5", "tissue_positions", "scalefactors"},
}

# LEGACY TSV sheet; JSON is canonical (docs/sample_sheet.md)
TSV_REQUIRED_COLUMNS = ("SAMPLE", "REP_ID", "assay_type", "sample_type", "PATH_to_bam")

# final outputs per mode, under bb_dir
BB_ALLELES = ("bb.Tallele.npz", "bb.Aallele.npz", "bb.Ballele.npz")
BB_GRID = ("bb.tsv.gz", "sample_ids.tsv")
BULK_TARGETS = BB_GRID + BB_ALLELES + ("bb.depth.npz", "bb.rdr.npz")
SINGLE_CELL_TARGETS = BB_GRID + BB_ALLELES + ("bb.Xcount.npz", "barcodes.tsv.gz")
COPYTYPING_TARGETS = BB_ALLELES + ("cnv_segments.tsv", "bb.Xcount.npz")


def is_url(path):
    """True if a sample-file path is a remote URL rather than a local path."""
    return str(path).startswith(URL_SCHEMES)


def get_alignment_index_path(alignment):
    """Return the co-located index path of an alignment: .crai for CRAM, else .bai."""
    return alignment + (".crai" if alignment.endswith(".cram") else ".bai")


def get_phasing_panel_path(phasing_panel):
    """Return a per-chromosome phasing panel path function."""
    return lambda chrname: os.path.join(phasing_panel, f"chr{chrname}.genotypes.bcf")


def get_genetic_map_path(gmap_path):
    """Return a per-chromosome genetic map path function.

    ``gmap_path`` is a full path with optional ``{chrname}`` placeholder:
      SHAPEIT5: ``/path/to/maps/chr{chrname}.mm10.gmap.gz``
      Eagle2:   ``/path/to/tables/genetic_map_mm10_withX.txt.gz`` (no placeholder)
    """
    return lambda chrname: gmap_path.format(chrname=chrname)
