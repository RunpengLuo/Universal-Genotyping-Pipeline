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
