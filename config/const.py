"""Assay-type, data-type, and reference constants shared by Snakefile and scripts."""

import os

NONBULK_ASSAYS = {"scATAC", "scRNA", "VISIUM", "VISIUM3prime"}
BULK_ASSAYS = {"bulkWGS", "bulkWGS-lr", "bulkWES"}
ALLOWED_ASSAY_TYPES = list(BULK_ASSAYS) + list(NONBULK_ASSAYS)
SPATIAL_ASSAYS = {"VISIUM", "VISIUM3prime"}
CHROM_ORDER = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]

# Reference genome versions that natively supported
REFVERS = ["hg19", "hg38", "chm13v2", "mm10"]

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
