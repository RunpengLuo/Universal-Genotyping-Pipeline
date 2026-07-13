#!/usr/bin/env python3
"""Parse and validate everything workflow/Snakefile needs before the DAG is built.

Runpeng Luo (2026-07-11)

``parse_workflow(config)`` is the only entry point the Snakefile calls; it returns
every name the rules read. ``parse_sample_file`` dispatches on the file extension
to the JSON or the legacy TSV parser, both emitting the same record schema.

Dependencies:
  const.py (config/, placed on sys.path by the Snakefile); stdlib only. The sample-file
  schema and the 10x layout live in const.py; this module holds only the logic.

Inputs
  config: the Snakemake config dict (keys: docs/reference.md)
  sample file: JSON, or a legacy TSV sheet (docs/sample_sheet.md)
Outputs:
  parse_workflow -> dict of workflow globals (keys listed in its docstring)
Notes/References:
  Sample-file schema, remote inputs, validation rules: docs/sample_sheet.md
  Config keys and output layout: docs/reference.md
"""

import csv
import json
import os
import sys

from const import (
    ALIGNMENT_FILES,
    ALLOWED_ASSAY_TYPES,
    ASSAY_TYPE2MODALITY,
    BULK_ASSAYS,
    BULK_TARGETS,
    COPYTYPING_TARGETS,
    FILE_KEYS,
    LONGREAD_ASSAYS,
    LONGREAD_PHASER,
    NONBULK_ASSAYS,
    PANEL_PHASER,
    PROVENANCE_KEYS,
    RANGER_LAYOUT,
    RANGER_SPATIAL_DIR,
    RECORD_KEYS,
    REFVERS,
    REQUIRED_FILES,
    REQUIRED_RECORD_KEYS,
    SINGLE_CELL_TARGETS,
    SPATIAL_ASSAYS,
    TSV_REQUIRED_COLUMNS,
    WORKFLOW_MODES,
    get_alignment_index_path,
    get_genetic_map_path,
    get_phasing_panel_path,
    is_url,
)


def _anchor(path, idx, rec):
    """Error prefix naming the offending record."""
    return (
        f"{path}: record {idx} (sample_id={rec.get('sample_id')!r}, "
        f"dataset_id={rec.get('dataset_id')!r}, assay_type={rec.get('assay_type')!r})"
    )


def parse_sample_file_json(path):
    """Read the JSON sample file and return its records.

    Args:
        path: Path to the JSON sample file.

    Returns:
        List of record dicts, with scalar fields coerced to str.

    Raises:
        ValueError: The file is not a records object/list, or a record is not an
            object.
    """
    with open(path) as fh:
        doc = json.load(fh)

    if isinstance(doc, dict):
        if "samples" not in doc:
            raise ValueError(f"{path}: object must hold a 'samples' list")
        records = doc["samples"]
    else:
        records = doc
    if not isinstance(records, list):
        raise ValueError(f"{path}: 'samples' must be a list of records")

    out = []
    for idx, rec in enumerate(records):
        if not isinstance(rec, dict):
            raise ValueError(f"{path}: record {idx} is not an object")
        norm = dict(rec)
        for key in ("sample_id", "dataset_id", "rdr_base_dataset_id", *PROVENANCE_KEYS):
            if key in norm and norm[key] is not None:
                norm[key] = str(norm[key])
        files = norm.get("files")
        if isinstance(files, dict):
            norm["files"] = {k: str(v) for k, v in files.items()}
        out.append(norm)
    return out


def expand_ranger_dir(assay_type, ranger_dir, files):
    """LEGACY (TSV only): fill a record's missing files from a local 10x outs/ dir.

    Used only by the legacy TSV path, which cannot name the files individually.
    Keys already in `files` win; alternate Ranger spellings (const.py RANGER_*)
    are probed on disk. See docs/sample_sheet.md, "Legacy TSV".

    Args:
        assay_type: Assay of the record.
        ranger_dir: Local Cell/Space Ranger outs/ directory.
        files: Files map parsed so far (bam, barcodes).

    Returns:
        The files map with the derivable keys filled in.

    Raises:
        ValueError: ranger_dir is a URL (a directory cannot be fetched).
    """
    if is_url(ranger_dir):
        raise ValueError(
            f"PATH_to_10x_ranger must be a local directory, got a URL: {ranger_dir}. "
            "Name each file explicitly in a JSON sample file to use remote inputs."
        )

    def _probe(names, in_spatial):
        prefix = RANGER_SPATIAL_DIR if in_spatial else ""
        paths = [os.path.join(ranger_dir, prefix, name) for name in names]
        return next((p for p in paths if os.path.exists(p)), paths[0])

    # keys the assay consumes, minus those the TSV names or derives itself
    derived = {
        key: _probe(*RANGER_LAYOUT[key])
        for key in REQUIRED_FILES[assay_type] - ALIGNMENT_FILES - {"barcodes"}
    }
    return {**derived, **files}


def parse_sample_file_tsv(path):
    """LEGACY: read a TSV sample sheet and return records in the JSON schema.

    Columns and their limits: docs/sample_sheet.md, "Legacy TSV".

    Args:
        path: Path to the TSV sample sheet.

    Returns:
        List of record dicts in the same schema as parse_sample_file_json.

    Raises:
        ValueError: A required column is missing.
    """
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"{path}: no rows")
    missing = [c for c in TSV_REQUIRED_COLUMNS if c not in rows[0]]
    if missing:
        raise ValueError(f"{path}: missing required column(s) {missing}")

    def _get(row, col):
        val = (row.get(col) or "").strip()
        return val or None

    records = []
    for row in rows:
        assay_type = row["assay_type"]
        alignment = row["PATH_to_bam"]
        files = {
            "alignment": alignment,
            "alignment_index": get_alignment_index_path(alignment),
        }
        barcodes = _get(row, "PATH_to_barcodes")
        if barcodes:
            files["barcodes"] = barcodes
        ranger_dir = _get(row, "PATH_to_10x_ranger")
        if ranger_dir and assay_type in NONBULK_ASSAYS:
            files = expand_ranger_dir(assay_type, ranger_dir, files)

        rec = {
            "sample_id": row["SAMPLE"],
            "dataset_id": row["REP_ID"],
            "assay_type": assay_type,
            "sample_type": row["sample_type"],
            "files": files,
        }
        base = _get(row, "RDR_BASE_REP_ID")
        if base:
            rec["rdr_base_dataset_id"] = base
        passage = _get(row, "passage")
        if passage:
            rec["passage_id"] = passage
        records.append(rec)
    return records


def parse_sample_file(path):
    """Parse a sample file, dispatching on extension (.json or .tsv).

    Args:
        path: Path to a JSON sample file or a legacy TSV sample sheet.

    Returns:
        List of record dicts. See docs/sample_sheet.md.

    Raises:
        ValueError: The extension is neither .json nor .tsv/.txt.
    """
    ext = os.path.splitext(path)[1].lower()
    if ext == ".json":
        return parse_sample_file_json(path)
    if ext in (".tsv", ".txt"):
        print(
            f"NOTE: {path} is a LEGACY TSV sheet; single-cell files are derived from "
            "PATH_to_10x_ranger and must be local. See docs/sample_sheet.md.",
            file=sys.stderr,
        )
        return parse_sample_file_tsv(path)
    raise ValueError(f"{path}: sample file must be .json or .tsv, got {ext!r}")


def validate_records(records, path, workflow_mode, sample_id, configured_assay_types):
    """Validate records against the spec, then the subset selected for this run.

    Args:
        records: Records from load_sample_records.
        path: Sample file path, for error messages.
        workflow_mode: bulk_genotyping | single_cell_genotyping | copytyping_preprocess.
        sample_id: The sample_id being processed.
        configured_assay_types: Assay types enabled for this run.

    Raises:
        ValueError: Any record violates the spec, or the selected records
            violate a mode/replicate rule.
    """
    single_cell = workflow_mode in ("single_cell_genotyping", "copytyping_preprocess")

    for idx, rec in enumerate(records):
        at = _anchor(path, idx, rec)

        unknown = set(rec) - RECORD_KEYS
        if unknown:
            raise ValueError(
                f"{at}: unknown key(s) {sorted(unknown)}; allowed: {sorted(RECORD_KEYS)}"
            )
        missing = [k for k in REQUIRED_RECORD_KEYS if rec.get(k) in (None, "")]
        if missing:
            raise ValueError(f"{at}: missing required key(s) {missing}")

        assay_type = rec["assay_type"]
        if assay_type not in ALLOWED_ASSAY_TYPES:
            raise ValueError(
                f"{at}: assay_type must be one of {sorted(ALLOWED_ASSAY_TYPES)}"
            )
        if rec["sample_type"] not in ("normal", "tumor"):
            raise ValueError(f"{at}: sample_type must be 'normal' or 'tumor'")
        if not isinstance(rec["files"], dict):
            raise ValueError(f"{at}: files must be an object")
        if "meta" in rec and not isinstance(rec["meta"], dict):
            raise ValueError(f"{at}: meta must be an object")

        files = rec["files"]
        unknown = set(files) - set(FILE_KEYS)
        if unknown:
            raise ValueError(
                f"{at}: unknown files key(s) {sorted(unknown)}; "
                f"allowed: {sorted(FILE_KEYS)}"
            )
        unused = set(files) - REQUIRED_FILES[assay_type]
        if unused:
            raise ValueError(
                f"{at}: files key(s) {sorted(unused)} are not used by {assay_type}"
            )
        required = REQUIRED_FILES[assay_type] if single_cell else ALIGNMENT_FILES
        for key in sorted(required):
            if not files.get(key):
                raise ValueError(
                    f"{at}: files.{key} is required for {assay_type}; "
                    f"got {sorted(files)}"
                )

    selected = [
        r
        for r in records
        if r["sample_id"] == sample_id and r["assay_type"] in configured_assay_types
    ]
    if not selected:
        raise ValueError(
            f"{path}: no records for sample_id={sample_id!r} with assay_type in "
            f"{sorted(configured_assay_types)}"
        )

    seen = {}
    for idx, rec in enumerate(selected):
        key = (rec["dataset_id"], rec["assay_type"])
        if key in seen:
            raise ValueError(
                f"{_anchor(path, idx, rec)}: duplicate (dataset_id, assay_type) {key}"
            )
        seen[key] = idx

    dataset2assays = {}
    for rec in selected:
        dataset2assays.setdefault(rec["dataset_id"], []).append(rec["assay_type"])
    for dataset_id, assays in dataset2assays.items():
        if len(assays) == 1:
            continue
        if any(a in BULK_ASSAYS for a in assays):
            raise ValueError(
                f"{path}: dataset_id={dataset_id!r} is reused across bulk assays {assays}; "
                "bulk dataset_ids must be unique"
            )
        if len(assays) > 2 or set(assays) != {"scRNA", "scATAC"}:
            raise ValueError(
                f"{path}: dataset_id={dataset_id!r} has assays {assays}; a shared dataset_id is "
                "only allowed for an scRNA + scATAC multiome pair"
            )

    dataset_ids = {r["dataset_id"] for r in selected}
    for idx, rec in enumerate(selected):
        base = rec.get("rdr_base_dataset_id")
        if not base:
            continue
        at = _anchor(path, idx, rec)
        if rec["sample_type"] != "tumor":
            raise ValueError(f"{at}: rdr_base_dataset_id is set on a non-tumor record")
        if base == rec["dataset_id"]:
            raise ValueError(f"{at}: rdr_base_dataset_id={base!r} is the record itself")
        if base not in dataset_ids:
            raise ValueError(
                f"{at}: rdr_base_dataset_id={base!r} is not a dataset_id in sample_id={sample_id!r}"
            )


def get_visium_layout(files, assay_type):
    """Map a spatial record's files onto their Space Ranger `spatial/` names.

    squidpy.read.visium() takes a directory, so process_rna_anndata stages these
    under <tmp>/spatial/ before reading. Names come from RANGER_LAYOUT (const.py);
    the canonical (first) spelling is used, since the staged dir is ours to name.

    Args:
        files: The record's files map.
        assay_type: VISIUM or VISIUM3prime.

    Returns:
        Dict of Space Ranger `spatial/` filename -> the record's path.
    """
    assert assay_type in SPATIAL_ASSAYS, f"not a spatial assay: {assay_type}"
    layout = {}
    for key, (names, in_spatial) in RANGER_LAYOUT.items():
        if in_spatial and files.get(key):
            layout[names[0]] = files[key]
    return layout


def select_records(path, workflow_mode, sample_id, configured_assay_types):
    """Parse + validate a sample file, then return the records this run processes.

    Args:
        path: JSON sample file (or legacy TSV sheet).
        workflow_mode: bulk_genotyping | single_cell_genotyping | copytyping_preprocess.
        sample_id: The sample_id to process.
        configured_assay_types: Assay types enabled for this run.

    Returns:
        List of the selected records, each with a `modality` key added.
    """
    records = parse_sample_file(path)
    validate_records(records, path, workflow_mode, sample_id, configured_assay_types)

    selected = []
    for rec in records:
        if (
            rec["sample_id"] != sample_id
            or rec["assay_type"] not in configured_assay_types
        ):
            continue
        rec = {**rec, "modality": ASSAY_TYPE2MODALITY[rec["assay_type"]]}
        if rec["sample_type"] == "tumor" and not rec.get("rdr_base_dataset_id"):
            print(
                f"NOTE: tumor dataset_id={rec['dataset_id']!r} has no "
                "rdr_base_dataset_id; RDR uses median normalization"
            )
        selected.append(rec)
    return selected


def get_sample_context(records):
    """Build the per-assay lookups the rules consume.

    Args:
        records: Output of select_records.

    Returns:
        Dict with:
          get_data: {(assay_type, dataset_id): files map}
          modality2files: {modality: [files map, ...]}
          assay2dataset_ids / assay2sample_types / assay2base_reps: per assay
            lists, aligned with each other; bulk assays are ordered normal-first
            so column 0 is the matched normal. Every allowed assay is present,
            defaulting to [].
    """
    modality2files = {}
    for rec in records:
        modality2files.setdefault(rec["modality"], []).append(rec["files"])

    by_assay = {}
    for rec in records:
        by_assay.setdefault(rec["assay_type"], []).append(rec)

    assay2dataset_ids, assay2sample_types, assay2base_reps = {}, {}, {}
    for assay_type in ALLOWED_ASSAY_TYPES:
        rows = by_assay.get(assay_type, [])
        if assay_type in BULK_ASSAYS:
            rows = sorted(rows, key=lambda r: r["sample_type"] != "normal")
        assay2dataset_ids[assay_type] = [r["dataset_id"] for r in rows]
        assay2sample_types[assay_type] = [r["sample_type"] for r in rows]
        assay2base_reps[assay_type] = [r.get("rdr_base_dataset_id", "") for r in rows]

    return {
        "get_data": {(r["assay_type"], r["dataset_id"]): r["files"] for r in records},
        "modality2files": modality2files,
        "assay2dataset_ids": assay2dataset_ids,
        "assay2sample_types": assay2sample_types,
        "assay2base_reps": assay2base_reps,
    }


def select_files(records, sample_type=None, longread_only=False):
    """Return the files maps of the records matching a sample_type / long-read filter."""
    return [
        r["files"]
        for r in records
        if (sample_type is None or r["sample_type"] == sample_type)
        and (not longread_only or r["assay_type"] in LONGREAD_ASSAYS)
    ]


def parse_assay_types(config, workflow_mode):
    """Validate config assay_types and keep only those the mode can run.

    Raises:
        ValueError: An assay type is not one of ALLOWED_ASSAY_TYPES.
    """
    assay_types = config["assay_types"]
    invalid = [a for a in assay_types if a not in ALLOWED_ASSAY_TYPES]
    if invalid:
        raise ValueError(
            f"invalid assay_types={invalid}; allowed: {sorted(ALLOWED_ASSAY_TYPES)}"
        )
    allowed = BULK_ASSAYS if workflow_mode == "bulk_genotyping" else NONBULK_ASSAYS
    return [a for a in assay_types if a in allowed]


def validate_mode(config, workflow_mode, assay_types):
    """Check the mode-specific config requirements.

    Raises:
        ValueError: copytyping_preprocess is missing a required input or its
            het_snp_vcf is not declared phased, or bulkWES is mixed with another
            bulk assay.
    """
    if workflow_mode == "copytyping_preprocess":
        for key in ("het_snp_vcf", "bb_file"):
            if config.get(key) is None:
                raise ValueError(f"{key} is required for copytyping_preprocess")
        if not config.get("het_snp_vcf_phased", True):
            raise ValueError(
                "copytyping_preprocess requires a phased het_snp_vcf: set "
                "het_snp_vcf_phased=true (it never genotypes or phases)"
            )
    elif workflow_mode == "bulk_genotyping":
        if "bulkWES" in assay_types and len(assay_types) > 1:
            raise ValueError(
                f"mixing bulkWES with other bulk assays is unsupported: {assay_types}"
            )


def select_datasets(records, dataset_ids, config_key, sample_id):
    """Return the records named by a config list of dataset_ids, in that order.

    Args:
        records: Output of select_records.
        dataset_ids: dataset_ids from config; each must be in this run.
        config_key: Config key being resolved, for error messages.
        sample_id: The sample_id being processed, for error messages.

    Returns:
        List of records, one per requested dataset_id.

    Raises:
        ValueError: A dataset_id is duplicated or absent from the run.
    """
    if len(set(dataset_ids)) != len(dataset_ids):
        raise ValueError(f"{config_key} has duplicate dataset_ids: {dataset_ids}")
    by_id = {r["dataset_id"]: r for r in records}
    missing = [d for d in dataset_ids if d not in by_id]
    if missing:
        raise ValueError(
            f"{config_key}={missing} not a dataset_id of sample_id={sample_id!r}; "
            f"available: {sorted(by_id)}"
        )
    return [by_id[d] for d in dataset_ids]


def get_genotype_records(config, records, sample_id):
    """Pick the records whose alignments are piled up to call germline SNPs.

    Uses config genotype_dataset_ids when set. Otherwise prefers a normal over a
    tumor, and a short-read assay over a long-read one: bcftools mpileup defaults
    are tuned for short reads, so a long-read alignment is the last resort even
    when one is present.

    Returns:
        List of records (one unless the config names several).

    Raises:
        ValueError: A named dataset_id is unknown, or the run has no records.
    """
    named = config.get("genotype_dataset_ids") or []
    if named:
        chosen = select_datasets(records, named, "genotype_dataset_ids", sample_id)
        non_normal = [r["dataset_id"] for r in chosen if r["sample_type"] != "normal"]
        if non_normal:
            print(
                f"WARN: genotype_dataset_ids includes non-normal dataset(s) "
                f"{non_normal}; germline SNPs may carry somatic signal"
            )
        return chosen

    def rank(rec):
        return (
            rec["sample_type"] != "normal",
            rec["assay_type"] in LONGREAD_ASSAYS,
        )

    chosen = sorted(records, key=rank)[:1]
    if not chosen:
        raise ValueError(f"no records to genotype for sample_id={sample_id!r}")
    rec = chosen[0]
    print(
        f"NOTE: genotype_dataset_ids unset; genotyping {rec['dataset_id']!r} "
        f"({rec['sample_type']}, {rec['assay_type']})"
    )
    return chosen


def get_phase_records(config, records, sample_id):
    """Pick the long-read records longphase reads.

    Uses config phase_dataset_ids when set, else a long-read normal, else a
    long-read tumor.

    Raises:
        ValueError: A named dataset is not long-read, or the run has none.
    """
    named = config.get("phase_dataset_ids") or []
    if named:
        if len(named) > 1:
            raise ValueError(
                f"phase_dataset_ids takes one dataset_id; longphase reads a single "
                f"alignment, got {named}"
            )
        chosen = select_datasets(records, named, "phase_dataset_ids", sample_id)
        short = [
            r["dataset_id"] for r in chosen if r["assay_type"] not in LONGREAD_ASSAYS
        ]
        if short:
            raise ValueError(
                f"phase_dataset_ids={short} are not long-read assays "
                f"({sorted(LONGREAD_ASSAYS)}); longphase needs long reads"
            )
        return chosen

    lr = [r for r in records if r["assay_type"] in LONGREAD_ASSAYS]
    lr = sorted(lr, key=lambda r: r["sample_type"] != "normal")
    if not lr:
        raise ValueError(
            f"phaser=longphase requires at least one long-read bulk assay "
            f"({sorted(LONGREAD_ASSAYS)}) in the sample file for "
            f"sample_id={sample_id!r}"
        )
    print(f"NOTE: phase_dataset_ids unset; longphase reads {lr[0]['dataset_id']!r}")
    return lr[:1]


def parse_phasing(config, records, sample_id, phaser, run_genotyping, run_phasing):
    """Resolve the genotyping alignment(s) and the phaser's reference inputs.

    Panel phasers (eagle/shapeit) need a genetic map and a reference panel, both
    checked here so a missing file fails at parse time rather than mid-DAG. The
    long-read phaser instead reads haplotypes from a long-read alignment.

    Returns:
        Dict with genotype_files, phase_files, get_genetic_map, get_phasing_panel.

    Raises:
        AssertionError: A gmap or panel file is missing.
        ValueError: The phaser is unknown, or a dataset selection is invalid.
    """
    out = {
        "genotype_files": None,
        "phase_files": None,
        "get_genetic_map": None,
        "get_phasing_panel": None,
    }
    if run_genotyping:
        out["genotype_files"] = [
            r["files"] for r in get_genotype_records(config, records, sample_id)
        ]
    if not run_phasing:
        return out

    if phaser in PANEL_PHASER:
        gmap_path = config.get("gmap_path")
        assert gmap_path, f"gmap_path required for {phaser}"
        get_genetic_map = get_genetic_map_path(gmap_path)
        missing_gmaps = [
            get_genetic_map(c)
            for c in config["chromosomes"]
            if not os.path.exists(get_genetic_map(c))
        ]
        assert not missing_gmaps, f"failed to locate gmap files: {missing_gmaps[:3]}"

        phasing_panel = config["phasing_panel"]
        assert os.path.isdir(phasing_panel), (
            f"failed to locate phasing panel: {phasing_panel}"
        )
        get_phasing_panel = get_phasing_panel_path(phasing_panel)
        missing_panels = [
            get_phasing_panel(c)
            for c in config["chromosomes"]
            if not os.path.exists(get_phasing_panel(c))
        ]
        assert not missing_panels, f"failed to locate panel files: {missing_panels[:3]}"

        out["get_genetic_map"] = get_genetic_map
        out["get_phasing_panel"] = get_phasing_panel
    elif phaser in LONGREAD_PHASER:
        out["phase_files"] = [
            r["files"] for r in get_phase_records(config, records, sample_id)
        ]
    else:
        raise ValueError(f"unknown phaser: {phaser}")
    return out


def build_final_targets(config, workflow_mode, assay_types, bulk_stream, msr_list):
    """Return the files `rule all` requests for this mode."""
    bb_dir = config["bb_dir"]
    targets = []
    if workflow_mode == "bulk_genotyping":
        targets += [
            f"{bb_dir}/MSR{msr}/{bulk_stream}/{f}"
            for msr in msr_list
            for f in BULK_TARGETS
        ]
    elif workflow_mode == "single_cell_genotyping":
        targets += [
            f"{bb_dir}/MSR{msr}/{assay_type}/{f}"
            for assay_type in assay_types
            for msr in msr_list
            for f in SINGLE_CELL_TARGETS
        ]
    elif workflow_mode == "copytyping_preprocess":
        targets += [
            f"{bb_dir}/{assay_type}/{f}"
            for assay_type in assay_types
            for f in COPYTYPING_TARGETS
        ]
    return targets


def parse_workflow(config):
    """Parse and validate the config + sample file into the workflow's globals.

    Args:
        config: The Snakemake config dict.

    Returns:
        Dict of the names workflow/Snakefile unpacks and the rules then read:
          workflow_mode, sample_id, assay_types, modalities, bulk_stream,
          msr_list, get_data, modality2files,
          assay2dataset_ids, assay2sample_types, assay2base_reps, phaser,
          run_genotyping, run_phasing, het_snp_vcf, phased_snp_vcf,
          require_genetic_map, genotype_files, phase_files,
          get_genetic_map, get_phasing_panel, final_targets.

    Raises:
        ValueError: The mode, assay types, sample file, or phaser is invalid.
    """
    workflow_mode = config["workflow_mode"]
    if workflow_mode not in WORKFLOW_MODES:
        raise ValueError(f"workflow_mode must be one of {list(WORKFLOW_MODES)}")

    sample_id = config["sample_id"]
    assay_types = parse_assay_types(config, workflow_mode)
    records = select_records(
        config["sample_file"], workflow_mode, sample_id, assay_types
    )
    # narrow to the assays present in the selected records
    assay_types = list(dict.fromkeys(r["assay_type"] for r in records))
    validate_mode(config, workflow_mode, assay_types)

    refvers = config.get("reference_version")
    if refvers not in REFVERS:
        print(
            f"WARNING: reference_version={refvers!r} is not natively supported: "
            f"{REFVERS}."
        )

    msr = config["params_combine_counts"]["min_snp_reads"]
    msr_list = [int(m) for m in (msr if isinstance(msr, list) else [msr])]
    # bulk allele/combine outputs live under allele_dir/{stream}/ (WGS family or WES)
    bulk_stream = "bulkWES" if "bulkWES" in assay_types else "bulkWGS"

    # A supplied het_snp_vcf replaces genotyping in every mode; het_snp_vcf_phased
    # declares whether it also replaces phasing. It is read only when a VCF is given.
    het_snp_vcf = config.get("het_snp_vcf")
    run_genotyping = het_snp_vcf is None
    run_phasing = True
    if het_snp_vcf is not None:
        if not os.path.exists(het_snp_vcf):
            raise ValueError(f"het_snp_vcf does not exist: {het_snp_vcf}")
        run_phasing = not bool(config.get("het_snp_vcf_phased", True))
        print(
            "NOTE: het_snp_vcf given -> skipping genotyping; "
            + (
                "het_snp_vcf_phased=false -> the VCF will be phased"
                if run_phasing
                else "het_snp_vcf_phased=true -> skipping phasing too"
            )
        )
    phased_snp_vcf = (
        config["phase_dir"] + "/phased_het_snps.vcf.gz" if run_phasing else het_snp_vcf
    )

    phaser = config.get("phaser", "undefined")
    phasing = {
        "genotype_files": None,
        "phase_files": None,
        "get_genetic_map": None,
        "get_phasing_panel": None,
    }
    if run_genotyping or run_phasing:
        phasing = parse_phasing(
            config, records, sample_id, phaser, run_genotyping, run_phasing
        )

    return {
        "workflow_mode": workflow_mode,
        "sample_id": sample_id,
        "assay_types": assay_types,
        "modalities": list(dict.fromkeys(r["modality"] for r in records)),
        "bulk_stream": bulk_stream,
        "msr_list": msr_list,
        "phaser": phaser,
        "run_genotyping": run_genotyping,
        "run_phasing": run_phasing,
        "het_snp_vcf": het_snp_vcf,
        "phased_snp_vcf": phased_snp_vcf,
        "require_genetic_map": run_phasing and phaser in PANEL_PHASER,
        "final_targets": build_final_targets(
            config, workflow_mode, assay_types, bulk_stream, msr_list
        ),
        **get_sample_context(records),
        **phasing,
    }
