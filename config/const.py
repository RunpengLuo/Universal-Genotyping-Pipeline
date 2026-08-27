"""Assay, reference and sample-file constants shared by the Snakefile and scripts.

Last update: 2026-08-12

Constants:
- WORKFLOW_MODES, REMOTE_MODES: run modes and remote input handling
- ALLOWED_ASSAY_TYPES, MULTIOME_ASSAYS, PSEUDOBULK_ASSAYS, ASSAY_TYPE2MODALITY:
  assay vocabulary
- SAMPLE_FILE_EXTS, SAMPLE_TYPES, RECORD_ID_KEYS, RECORD_ID_PATTERN,
  REQUIRED_RECORD_KEYS, REQUIRED_FILES: sample-file schema
- TUMOR_GENOTYPING_MODES: tumor-only genotyping strategies
- SNP_PANEL_EXTS: snp_panel extensions bcftools --targets-file can parse. It picks the
  parser from the filename suffix, so .bcf and .vcf.bgz reach the tab parser and match
  nothing without erroring (samtools/bcftools#690, #1368)
- BULK_TARGETS, SINGLE_CELL_TARGETS, COPYTYPING_TARGETS: the per-mode final outputs
- REFVERS, REFVERS_ALIAS, SPECIES2SEXCHROM: reference builds and their spellings
- RANGER_*, RANGER_LAYOUT: 10x Cell and Space Ranger filenames
- REPLISEQ_*, REPLI_LIFTOVER, LIFTOVER_CHAIN_URLS: Repli-seq sources and chains
- URL_SCHEMES: remote input schemes a sample-file value may use
- VCF_COLUMNS, VCF_SAMPLE_COLUMNS, GTF_COLUMNS: the parsed file schemas
Functions:
- canonical_refver, is_known_refver: fold a reference_version spelling
- get_phasing_panel_path, get_genetic_map_path: per-chromosome path builders
"""

import os

##################################################
# Pipeline modes
WORKFLOW_MODES = ("bulk_genotyping", "single_cell_genotyping", "copytyping_preprocess")

# Remote input fetch mode
REMOTE_MODES = ("storage", "stream")

# Remote input schemes
URL_SCHEMES = ("http://", "https://", "ftp://", "s3://")

# params_combine_counts.rdr_normalization mode
RDR_NORMALIZATIONS = ("auto", "median", "normal")

# params_genotype_snps.tumor_genotyping_mode
TUMOR_GENOTYPING_MODES = ("vaf_cutoff", "clonal_loh_hmm")

# snp_panel file extensions
SNP_PANEL_EXTS = (".vcf.gz",)

##################################################
# sample-file schema
SAMPLE_FILE_EXTS = (".json", ".tsv", ".txt")
SAMPLE_TYPES = ("normal", "tumor")
RECORD_ID_KEYS = ("sample_id", "dataset_id")
RECORD_ID_PATTERN = r"[A-Za-z0-9_-]+"
REQUIRED_RECORD_KEYS = (
    "sample_id",
    "dataset_id",
    "assay_type",
    "sample_type",
    "reference_version",
    "files",
)
OPTIONAL_RECORD_KEYS = ("rdr_base_dataset_id",)
PROVENANCE_KEYS = ("passage_id", "platform", "cancer_type")
SCALAR_RECORD_KEYS = (
    tuple(k for k in REQUIRED_RECORD_KEYS if k != "files")
    + OPTIONAL_RECORD_KEYS
    + PROVENANCE_KEYS
)
FILES_COLUMN_PREFIX = "files."

##################################################
# Input/output filenames
ALIGNMENT_FILES = {"alignment", "alignment_index"}
REQUIRED_FILES = {
    "bulkWGS": ALIGNMENT_FILES,
    "bulkWGS-lr": ALIGNMENT_FILES,
    "bulkWES": ALIGNMENT_FILES,
    "scDNA": ALIGNMENT_FILES,
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
# Outputs under <bb_dir>
BB_ALLELES = ("bb.Tallele.npz", "bb.Aallele.npz", "bb.Ballele.npz")
BB_GRID = ("bb.tsv.gz", "sample_ids.tsv")
BULK_TARGETS = BB_GRID + BB_ALLELES + ("bb.depth.npz", "bb.rdr.npz")
SINGLE_CELL_TARGETS = BB_GRID + BB_ALLELES + ("bb.Xcount.npz", "barcodes.tsv.gz")
COPYTYPING_TARGETS = BB_ALLELES + ("bb.tsv.gz", "bb.Xcount.npz")

##################################################
# Supported phasing softwares
LONGREAD_PHASER = {"longphase"}
PANEL_PHASER = {"eagle", "shapeit"}


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


##################################################
# Supported sequencing assays
# single-cell libraries the pipeline processes pooled, through the bulk path only
PSEUDOBULK_ASSAYS = {"scDNA"}
BULK_ASSAYS = {"bulkWGS", "bulkWGS-lr", "bulkWES"} | PSEUDOBULK_ASSAYS
BULK_LR_ASSAYS = {"bulkWGS-lr"}
NONBULK_ASSAYS = {"scATAC", "scRNA", "VISIUM", "VISIUM3prime"}
SPATIAL_ASSAYS = {"VISIUM", "VISIUM3prime"}
ALLOWED_ASSAY_TYPES = list(BULK_ASSAYS) + list(NONBULK_ASSAYS)

MULTIOME_ASSAYS = {"scRNA", "scATAC"}

# bulk-genotyping assay preference
GT_ASSAY_ORD = {"bulkWGS": 0, "bulkWGS-lr": 1, "bulkWES": 2, "scDNA": 3}

ASSAY_TYPE2MODALITY = {
    "bulkWGS": "DNA",
    "bulkWGS-lr": "DNA",
    "bulkWES": "DNA",
    "scDNA": "DNA",
    "scATAC": "DNA",
    "scRNA": "RNA",
    "VISIUM": "RNA",
    "VISIUM3prime": "RNA",
}

##################################################
# Native supported reference genome versions
REFVERS = ["hg19", "hg38", "chm13v2", "mm10"]

# Species
SPECIES = ("human", "mouse")
SPECIES2SEXCHROM = {
    "human": {"X": 23, "Y": 24},
    "mouse": {"X": 20, "Y": 21},
}

REFVERS_ALIAS = {
    "hg19": ["GRCh37", "b37", "hs37", "hs37d5"],
    "hg38": ["GRCh38", "hs38", "hs38DH", "GRCh38.p13", "GRCh38_no_alt"],
    "chm13v2": ["CHM13v2.0", "CHM13", "T2T", "T2T-CHM13v2", "T2T-CHM13v2.0"],
    "mm10": ["GRCm38", "MGSCv38"],
}


def canonical_refver(value):
    """Convert reference version to canonical form if supported."""
    _ALIAS2REFVER = {
        alias.lower(): canon
        for canon, aliases in REFVERS_ALIAS.items()
        for alias in (canon, *aliases)
    }
    key = str(value).strip().lower()
    return _ALIAS2REFVER.get(key, key)


def is_known_refver(value):
    """True when ``value`` names a reference version the pipeline natively supports."""
    return canonical_refver(value) in REFVERS


##################################################
# Repli-seq for replication timing RD correction.
# Liftover hg19 to other reference versions.
REPLISEQ_REFVERS = ("hg19", "hg38", "chm13v2")
REPLI_LIFTOVER = ("hg38", "chm13v2")
LIFTOVER_CHAIN_URLS = {
    "hg38": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
    "chm13v2": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHs1.over.chain.gz",
}

# hg19 Repli-seq data
UCSC_REPLISEQ_BASE = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq"
)
REPLISEQ_BIGWIG_FILES = (
    "wgEncodeUwRepliSeqBg02esWaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqBjWaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqBjWaveSignalRep2.bigWig",
    "wgEncodeUwRepliSeqGm06990WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqGm12801WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqGm12812WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqGm12813WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqHelas3WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqImr90WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqSknshWaveSignalRep1.bigWig",
)

##################################################
# 10x {Cell,Space}-Ranger output filenames
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

##################################################
# VCF fixed columns, then the pair a single-sample genotyped VCF adds
VCF_COLUMNS = [
    "#CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
]
VCF_SAMPLE_COLUMNS = ["FORMAT", "SAMPLE"]

##################################################
# Gene annotation GTF file required columns
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
