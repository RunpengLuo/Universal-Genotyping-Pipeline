#!/usr/bin/env python3
"""Parse and validate everything workflow/Snakefile needs before the DAG is built.

Runpeng Luo (2026-07-11)

``parse_workflow(config)`` is the only entry point the Snakefile calls; it returns
every name the rules read, and inlines the whole flow as ``# ===`` blocks. Three
module-level helpers are shared: ``parse_sample_file_json`` / ``parse_sample_file_tsv``
(dispatched on the sample-file extension, both emitting the same record schema) and
``validate_records`` (reused by resources/scripts/validate_sample_file.py).

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
    LONGREAD_ASSAYS,
    LONGREAD_PHASER,
    NONBULK_ASSAYS,
    OPTIONAL_FILES,
    OPTIONAL_RECORD_KEYS,
    PANEL_PHASER,
    REPLISEQ_REFVERS,
    PROVENANCE_KEYS,
    RANGER_LAYOUT,
    RANGER_SPATIAL_DIR,
    RDR_NORMALIZATIONS,
    REFVERS,
    REQUIRED_FILES,
    REQUIRED_RECORD_KEYS,
    SINGLE_CELL_TARGETS,
    TSV_REQUIRED_COLUMNS,
    WORKFLOW_MODES,
    get_genetic_map_path,
    get_phasing_panel_path,
    is_url,
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
        for key in ("sample_id", "dataset_id", *OPTIONAL_RECORD_KEYS, *PROVENANCE_KEYS):
            if key in norm and norm[key] is not None:
                norm[key] = str(norm[key])
        files = norm.get("files")
        if isinstance(files, dict):
            norm["files"] = {k: str(v) for k, v in files.items()}
        out.append(norm)
    return out


def parse_sample_file_tsv(path):
    """LEGACY: read a TSV sample sheet and return records in the JSON schema.

    Columns and their limits: docs/sample_sheet.md, "Legacy TSV". Single-cell files
    are derived from a local PATH_to_10x_ranger directory via the nested
    ``expand_ranger_dir`` (alternate Ranger spellings, const.py RANGER_*, probed on
    disk; keys already named win).

    Args:
        path: Path to the TSV sample sheet.

    Returns:
        List of record dicts in the same schema as parse_sample_file_json.

    Raises:
        ValueError: A required column is missing, or a PATH_to_10x_ranger is a URL.
    """

    def expand_ranger_dir(assay_type, ranger_dir, files):
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

    def _get(row, col):
        val = (row.get(col) or "").strip()
        return val or None

    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"{path}: no rows")
    missing = [c for c in TSV_REQUIRED_COLUMNS if c not in rows[0]]
    if missing:
        raise ValueError(f"{path}: missing required column(s) {missing}")

    records = []
    for row in rows:
        assay_type = row["assay_type"]
        alignment = row["PATH_to_bam"]
        files = {
            "alignment": alignment,
            "alignment_index": alignment
            + (".crai" if alignment.endswith(".cram") else ".bai"),
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


def validate_records(records, path, workflow_mode, sample_id, configured_assay_types):
    """Validate records against the spec, then the subset selected for this run.

    Mutates each record's ``files`` map in place to keep only the keys the assay
    reads. Shared by ``parse_workflow`` (DAG build) and the standalone
    ``resources/scripts/validate_sample_file.py``.

    Args:
        records: Records from parse_sample_file_{json,tsv}.
        path: Sample file path, for error messages.
        workflow_mode: bulk_genotyping | single_cell_genotyping | copytyping_preprocess.
        sample_id: The sample_id being processed.
        configured_assay_types: Assay types enabled for this run.

    Raises:
        ValueError: Any record violates the spec, or the selected records violate a
            mode/replicate rule.
    """

    def _anchor(idx, rec):
        """Error prefix naming the offending record."""
        return (
            f"{path}: record {idx} (sample_id={rec.get('sample_id')!r}, "
            f"dataset_id={rec.get('dataset_id')!r}, assay_type={rec.get('assay_type')!r})"
        )

    single_cell = workflow_mode in ("single_cell_genotyping", "copytyping_preprocess")

    for idx, rec in enumerate(records):
        at = _anchor(idx, rec)

        req_missing = [k for k in REQUIRED_RECORD_KEYS if rec.get(k) in (None, "")]
        if req_missing:
            raise ValueError(f"{at}: missing required key(s) {req_missing}")

        assay_type = rec["assay_type"]
        if assay_type not in ALLOWED_ASSAY_TYPES:
            raise ValueError(
                f"{at}: assay_type must be one of {sorted(ALLOWED_ASSAY_TYPES)}"
            )
        if rec["sample_type"] not in ("normal", "tumor"):
            raise ValueError(f"{at}: sample_type must be 'normal' or 'tumor'")
        if not isinstance(rec["files"], dict):
            raise ValueError(f"{at}: files must be an object")

        readable = REQUIRED_FILES[assay_type] | OPTIONAL_FILES.get(assay_type, set())
        files = {k: v for k, v in rec["files"].items() if k in readable}
        ignored = set(rec["files"]) - set(files)
        if ignored:
            print(
                f"NOTE: {at}: ignoring files key(s) {sorted(ignored)}; "
                f"{assay_type} reads {sorted(readable)}",
                file=sys.stderr,
            )
        rec["files"] = files

        required_files = REQUIRED_FILES[assay_type] if single_cell else ALIGNMENT_FILES
        for key in sorted(required_files):
            if not files.get(key):
                raise ValueError(
                    f"{at}: files.{key} is required for {assay_type}; got {sorted(files)}"
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
                f"{_anchor(idx, rec)}: duplicate (dataset_id, assay_type) {key}"
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
        at = _anchor(idx, rec)
        if rec["sample_type"] != "tumor":
            raise ValueError(f"{at}: rdr_base_dataset_id is set on a non-tumor record")
        if base == rec["dataset_id"]:
            raise ValueError(f"{at}: rdr_base_dataset_id={base!r} is the record itself")
        if base not in dataset_ids:
            raise ValueError(
                f"{at}: rdr_base_dataset_id={base!r} is not a dataset_id in sample_id={sample_id!r}"
            )


def parse_workflow(config):
    """Parse and validate the config + sample file into the workflow's globals.

    Args:
        config: The Snakemake config dict.

    Returns:
        Dict of the names workflow/Snakefile unpacks and the rules then read:
          workflow_mode, sample_id, assay_types, modalities, msr_list, phaser,
          run_genotyping, run_phasing, het_snp_vcf, phased_snp_vcf,
          require_genetic_map, final_targets, get_data, modality2files,
          assay2dataset_ids, assay2sample_types, assay2base_reps, genotype_files,
          phase_files, get_genetic_map, get_phasing_panel, segment_bed, bedpe_files,
          wes_targets_files, has_breakpoints, window_streams, do_repliseq,
          window_size_wgs, window_size_wes.

    Raises:
        ValueError: The mode, assay types, sample file, or phaser is invalid.
    """
    path = config["sample_file"]

    def select_datasets(records, dataset_ids, config_key):
        """Records named by a config dataset_id list, in that order."""
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

    # === workflow mode + sample_id ===
    workflow_mode = config["workflow_mode"]
    if workflow_mode not in WORKFLOW_MODES:
        raise ValueError(f"workflow_mode must be one of {list(WORKFLOW_MODES)}")
    sample_id = config["sample_id"]

    # === assay_types: validate against the schema, keep those this mode runs ===
    configured = config["assay_types"]
    invalid = [a for a in configured if a not in ALLOWED_ASSAY_TYPES]
    if invalid:
        raise ValueError(
            f"invalid assay_types={invalid}; allowed: {sorted(ALLOWED_ASSAY_TYPES)}"
        )
    allowed = BULK_ASSAYS if workflow_mode == "bulk_genotyping" else NONBULK_ASSAYS
    configured = [a for a in configured if a in allowed]

    # === load the sample file (json or legacy tsv) ===
    ext = os.path.splitext(path)[1].lower()
    if ext == ".json":
        records = parse_sample_file_json(path)
    elif ext in (".tsv", ".txt"):
        print(
            f"NOTE: {path} is a LEGACY TSV sheet; single-cell files are derived from "
            "PATH_to_10x_ranger and must be local. See docs/sample_sheet.md.",
            file=sys.stderr,
        )
        records = parse_sample_file_tsv(path)
    else:
        raise ValueError(f"{path}: sample file must be .json or .tsv, got {ext!r}")

    # === validate against the spec + selection rules (mutates files in place) ===
    validate_records(records, path, workflow_mode, sample_id, configured)

    # === select this run's records + add modality ===
    records = [
        {**rec, "modality": ASSAY_TYPE2MODALITY[rec["assay_type"]]}
        for rec in records
        if rec["sample_id"] == sample_id and rec["assay_type"] in configured
    ]

    # === RDR normalization policy: drop/keep each bulk tumor's rdr_base ===
    rdr_normalization = config["params_combine_counts"].get("rdr_normalization", "auto")
    if rdr_normalization not in RDR_NORMALIZATIONS:
        raise ValueError(
            f"params_combine_counts.rdr_normalization must be one of "
            f"{list(RDR_NORMALIZATIONS)}, got {rdr_normalization!r}"
        )
    unbased, ignored = [], []
    for i, rec in enumerate(records):
        is_bulk_tumor = (
            rec["assay_type"] in BULK_ASSAYS and rec["sample_type"] == "tumor"
        )
        base = rec.get("rdr_base_dataset_id")
        if is_bulk_tumor:
            if rdr_normalization == "median" and base:
                ignored.append(rec["dataset_id"])
                records[i] = {
                    k: v for k, v in rec.items() if k != "rdr_base_dataset_id"
                }
            elif rdr_normalization != "median" and not base:
                unbased.append(rec["dataset_id"])
    if rdr_normalization == "normal" and unbased:
        raise ValueError(
            f"rdr_normalization='normal' requires rdr_base_dataset_id on every bulk "
            f"tumor, missing on: {sorted(unbased)}. Set it in the sample file, or use "
            "rdr_normalization='auto' to median-normalize these."
        )
    if ignored:
        print(
            f"NOTE: rdr_normalization='median' -> ignoring rdr_base_dataset_id on "
            f"{len(ignored)} tumor(s): {sorted(ignored)}"
        )
    if unbased:
        print(
            f"NOTE: {len(unbased)} tumor(s) have no rdr_base_dataset_id; RDR uses "
            f"median normalization: {sorted(unbased)}"
        )

    # assay types actually present in the selected records
    assay_types = list(dict.fromkeys(r["assay_type"] for r in records))

    # === copytyping_preprocess requirements ===
    if workflow_mode == "copytyping_preprocess":
        assert config.get("het_snp_vcf") is not None, (
            "het_snp_vcf is required for copytyping_preprocess"
        )
        assert config.get("het_snp_vcf_phased", True), (
            "het_snp_vcf must be phased for copytyping_preprocess"
        )
        assert config.get("bb_file") is not None, (
            "bb_file is required for copytyping_preprocess"
        )

    # === reference_version sanity ===
    refvers = config.get("reference_version")
    if refvers not in REFVERS:
        print(
            f"WARNING: reference_version={refvers!r} is not natively supported: "
            f"{REFVERS}."
        )

    # === min_snp_reads sweep (one MSR{msr}/ subdir per value) ===
    msr = config["params_combine_counts"]["min_snp_reads"]
    msr_list = [int(m) for m in (msr if isinstance(msr, list) else [msr])]

    # === segment BED + per-stream window build ===
    # region.bed (arm) stays config["region_bed"]; bulk builds aux/segment.bed
    # (region_id + seg_id) and the per-stream window BEDs, read by the bulk rules.
    bedpe_files = list(
        dict.fromkeys(
            r["files"]["breakpoint_bedpe"]
            for r in records
            if "breakpoint_bedpe" in r["files"]
        )
    )
    wes_targets_files = list(
        dict.fromkeys(
            r["files"]["wes_targets_bed"]
            for r in records
            if "wes_targets_bed" in r["files"]
        )
    )
    has_breakpoints = len(bedpe_files) > 0
    is_bulk = workflow_mode == "bulk_genotyping"
    window_streams = sorted({"wes" if at == "bulkWES" else "wgs" for at in assay_types})
    do_repliseq = config.get("reference_version") in REPLISEQ_REFVERS
    pp = config.get("params_build_windows") or {}
    window_size_wgs = int(pp.get("window_size_wgs") or 1000)
    window_size_wes = int(pp.get("window_size_wes") or 267)
    if is_bulk and "wes" in window_streams and not wes_targets_files:
        raise ValueError(
            "building the WES window BED needs files.wes_targets_bed on a bulkWES record"
        )
    segment_bed = (
        config["aux_dir"] + "/segment.bed" if is_bulk else config["region_bed"]
    )
    if is_bulk:
        print(
            f"NOTE: bulk -> segment BED ({segment_bed}); "
            f"{len(bedpe_files)} breakpoint_bedpe file(s), region_id=arm, seg_id=chunk; "
            f"window streams: {window_streams}"
        )

    # === genotyping / phasing switches (a het_snp_vcf short-circuits the front) ===
    het_snp_vcf = config.get("het_snp_vcf")
    run_genotyping = het_snp_vcf is None
    run_phasing = True
    phased_snp_vcf = config["phase_dir"] + "/phased_het_snps.vcf.gz"
    if het_snp_vcf is not None:
        assert os.path.exists(het_snp_vcf), f"het_snp_vcf does not exist: {het_snp_vcf}"
        # default het_snp_vcf_phased=true -> the VCF is taken as phased, phasing skipped
        run_phasing = not bool(config.get("het_snp_vcf_phased", True))
        if not run_phasing:
            phased_snp_vcf = het_snp_vcf

    # === phaser reference inputs (genotype + phase record selection) ===
    phaser = config.get("phaser", "undefined")
    genotype_files = None
    phase_files = None
    get_genetic_map = None
    get_phasing_panel = None
    if run_genotyping:
        named = config.get("genotype_dataset_ids") or []
        if named:
            chosen = select_datasets(records, named, "genotype_dataset_ids")
            non_normal = [
                r["dataset_id"] for r in chosen if r["sample_type"] != "normal"
            ]
            if non_normal:
                print(
                    f"WARN: genotype_dataset_ids includes non-normal dataset(s) "
                    f"{non_normal}; germline SNPs may carry somatic signal"
                )
        else:
            chosen = sorted(
                records,
                key=lambda r: (
                    r["sample_type"] != "normal",
                    r["assay_type"] in LONGREAD_ASSAYS,
                ),
            )[:1]
            if not chosen:
                raise ValueError(f"no records to genotype for sample_id={sample_id!r}")
            r0 = chosen[0]
            print(
                f"NOTE: genotype_dataset_ids unset; genotyping {r0['dataset_id']!r} "
                f"({r0['sample_type']}, {r0['assay_type']})"
            )
        genotype_files = [r["files"] for r in chosen]
    if run_phasing:
        if phaser in PANEL_PHASER:
            gmap_path = config.get("gmap_path")
            assert gmap_path, f"gmap_path required for {phaser}"
            get_genetic_map = get_genetic_map_path(gmap_path)
            missing_gmaps = [
                get_genetic_map(c)
                for c in config["chromosomes"]
                if not os.path.exists(get_genetic_map(c))
            ]
            assert not missing_gmaps, (
                f"failed to locate gmap files: {missing_gmaps[:3]}"
            )

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
            assert not missing_panels, (
                f"failed to locate panel files: {missing_panels[:3]}"
            )
        elif phaser in LONGREAD_PHASER:
            named = config.get("phase_dataset_ids") or []
            if named:
                if len(named) > 1:
                    raise ValueError(
                        f"phase_dataset_ids takes one dataset_id; longphase reads a single "
                        f"alignment, got {named}"
                    )
                chosen = select_datasets(records, named, "phase_dataset_ids")
                short = [
                    r["dataset_id"]
                    for r in chosen
                    if r["assay_type"] not in LONGREAD_ASSAYS
                ]
                if short:
                    raise ValueError(
                        f"phase_dataset_ids={short} are not long-read assays "
                        f"({sorted(LONGREAD_ASSAYS)}); longphase needs long reads"
                    )
            else:
                lr = [r for r in records if r["assay_type"] in LONGREAD_ASSAYS]
                lr = sorted(lr, key=lambda r: r["sample_type"] != "normal")
                if not lr:
                    raise ValueError(
                        f"phaser=longphase requires at least one long-read bulk assay "
                        f"({sorted(LONGREAD_ASSAYS)}) in the sample file for "
                        f"sample_id={sample_id!r}"
                    )
                print(
                    f"NOTE: phase_dataset_ids unset; longphase reads {lr[0]['dataset_id']!r}"
                )
                chosen = lr[:1]
            phase_files = [r["files"] for r in chosen]
        else:
            raise ValueError(f"unknown phaser: {phaser}")

    # === per-assay lookups the rules consume (bulk ordered normal-first) ===
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
    get_data = {(r["assay_type"], r["dataset_id"]): r["files"] for r in records}

    # === final targets `rule all` requests for this mode ===
    bb_dir = config["bb_dir"]
    if workflow_mode == "bulk_genotyping":
        final_targets = [
            f"{bb_dir}/MSR{m}/bulk/{f}" for m in msr_list for f in BULK_TARGETS
        ]
    elif workflow_mode == "single_cell_genotyping":
        final_targets = [
            f"{bb_dir}/MSR{m}/{at}/{f}"
            for at in assay_types
            for m in msr_list
            for f in SINGLE_CELL_TARGETS
        ]
    else:  # copytyping_preprocess
        final_targets = [
            f"{bb_dir}/{at}/{f}" for at in assay_types for f in COPYTYPING_TARGETS
        ]

    # === assemble the globals the Snakefile unpacks ===
    return {
        "workflow_mode": workflow_mode,
        "sample_id": sample_id,
        "assay_types": assay_types,
        "modalities": list(dict.fromkeys(r["modality"] for r in records)),
        "msr_list": msr_list,
        "phaser": phaser,
        "run_genotyping": run_genotyping,
        "run_phasing": run_phasing,
        "het_snp_vcf": het_snp_vcf,
        "phased_snp_vcf": phased_snp_vcf,
        "require_genetic_map": run_phasing and phaser in PANEL_PHASER,
        "final_targets": final_targets,
        "get_data": get_data,
        "modality2files": modality2files,
        "assay2dataset_ids": assay2dataset_ids,
        "assay2sample_types": assay2sample_types,
        "assay2base_reps": assay2base_reps,
        "genotype_files": genotype_files,
        "phase_files": phase_files,
        "get_genetic_map": get_genetic_map,
        "get_phasing_panel": get_phasing_panel,
        "segment_bed": segment_bed,
        "bedpe_files": bedpe_files,
        "wes_targets_files": wes_targets_files,
        "has_breakpoints": has_breakpoints,
        "window_streams": window_streams,
        "do_repliseq": do_repliseq,
        "window_size_wgs": window_size_wgs,
        "window_size_wes": window_size_wes,
    }
