#!/usr/bin/env python3
"""Parse and validate everything workflow/Snakefile needs before the DAG is built.

Runpeng Luo
Last update: 2026-08-06

Dependencies:
  const.py (sample-file schema, 10x layout) and script_utils (get_chr_sizes,
  strip_chr_prefix), both placed on sys.path by the Snakefile.

Inputs
  config: the Snakemake config dict (keys: docs/reference.md)
  sample file: JSON or TSV encoding (docs/sample_sheet.md)
Outputs:
  parse_workflow -> dict of workflow globals (keys listed in its docstring)
Notes/References:
  Sample-file schema, remote inputs, validation rules: docs/sample_sheet.md
  Config keys and output layout: docs/reference.md
"""

import csv
import json
import os

from const import (
    ALIGNMENT_FILES,
    ALLOWED_ASSAY_TYPES,
    ASSAY_TYPE2MODALITY,
    BULK_ASSAYS,
    BULK_TARGETS,
    COPYTYPING_TARGETS,
    FILES_COLUMN_PREFIX,
    LONGREAD_ASSAYS,
    LONGREAD_PHASER,
    NONBULK_ASSAYS,
    OPTIONAL_FILES,
    PANEL_PHASER,
    REPLISEQ_REFVERS,
    RDR_NORMALIZATIONS,
    REFVERS,
    REQUIRED_FILES,
    REQUIRED_RECORD_KEYS,
    SCALAR_RECORD_KEYS,
    SINGLE_CELL_TARGETS,
    SPECIES,
    WORKFLOW_MODES,
    canonical_refver,
    get_genetic_map_path,
    is_known_refver,
    get_phasing_panel_path,
    is_url,
)
from io_utils import get_chr_sizes
from utils import logging_snakemake, strip_chr_prefix


def parse_sample_file_json(path):
    """Read the JSON sample file and return its records.

    Args:
        path: Path to the JSON sample file.

    Returns:
        List of record dicts, with scalar fields coerced to str.

    Raises:
        AssertionError: The file is not a records object/list, or a record is not
            an object.
    """
    with open(path) as fh:
        doc = json.load(fh)

    if isinstance(doc, dict):
        assert "samples" in doc, f"{path}: object must hold a 'samples' list"
        records = doc["samples"]
    else:
        records = doc
    assert isinstance(records, list), f"{path}: 'samples' must be a list"

    out = []
    for idx, rec in enumerate(records):
        assert isinstance(rec, dict), f"{path}: record {idx} is not an object"
        norm = dict(rec)
        for key in SCALAR_RECORD_KEYS:
            if key in norm and norm[key] is not None:
                norm[key] = str(norm[key])
        files = norm.get("files")
        if isinstance(files, dict):
            norm["files"] = {k: str(v) for k, v in files.items() if v is not None}
        out.append(norm)
    return out


def parse_sample_file_tsv(path):
    """Read a TSV sample sheet and return records in the JSON schema.

    The TSV is a flat encoding of the same schema, not a reduced one: columns are
    the record keys, and each input is its own ``files.<key>`` column. An empty cell
    omits the key. Spec: docs/sample_sheet.md, "TSV".

    Args:
        path: Path to the TSV sample sheet.

    Returns:
        List of record dicts in the same schema as parse_sample_file_json.

    Raises:
        AssertionError: The file has no rows, or a required column is missing.
    """
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows, f"{path}: no rows"

    required_columns = [k for k in REQUIRED_RECORD_KEYS if k != "files"]
    missing = [c for c in required_columns if c not in rows[0]]
    assert not missing, f"{path}: missing required column(s) {missing}"
    assert any(c.startswith(FILES_COLUMN_PREFIX) for c in rows[0]), (
        f"{path}: no {FILES_COLUMN_PREFIX}* column"
    )

    records = []
    for row in rows:
        rec, files = {}, {}
        for col, val in row.items():
            val = (val or "").strip()
            if not val or col is None:
                continue
            if col.startswith(FILES_COLUMN_PREFIX):
                files[col[len(FILES_COLUMN_PREFIX) :]] = val
            else:
                rec[col] = val
        rec["files"] = files
        records.append(rec)
    return records


def select_records(records, sample_id, configured_assay_types, reference_version):
    """This run's records: one sample_id, the configured assays, one genome build.

    The single definition of "selected", used by both ``validate_records`` (whose
    replicate rules must see exactly the run's records) and ``parse_workflow``.

    Args:
        records: Records from parse_sample_file_{json,tsv}.
        sample_id: The sample_id being processed.
        configured_assay_types: Assay types enabled for this run.
        reference_version: Canonical build to keep; records of any other are dropped.

    Returns:
        The matching records, in file order.
    """
    return [
        r
        for r in records
        if r["sample_id"] == sample_id
        and r["assay_type"] in configured_assay_types
        and canonical_refver(r["reference_version"]) == reference_version
    ]


def _anchor(path, idx, rec):
    """Error prefix naming the offending record."""
    return (
        f"{path}: record {idx} (sample_id={rec.get('sample_id')!r}, "
        f"dataset_id={rec.get('dataset_id')!r}, assay_type={rec.get('assay_type')!r})"
    )


def _refver_counts(records, sample_id):
    """`refver (n)` summary of one sample_id's records, for error messages."""
    present = {}
    for rec in records:
        if rec["sample_id"] == sample_id:
            refver = canonical_refver(rec["reference_version"])
            present[refver] = present.get(refver, 0) + 1
    return ", ".join(f"{rv} ({n})" for rv, n in sorted(present.items())) or "no records"


def require_record_keys(records, path):
    """Every record carries every REQUIRED_RECORD_KEYS entry.

    Runs before anything reads a record by key, so a malformed sample file fails
    naming the record rather than with a bare KeyError from a later subset.

    Args:
        records: Records from parse_sample_file_{json,tsv}.
        path: Sample file path, for error messages.

    Raises:
        AssertionError: A record is missing a required key or leaves it empty.
    """
    for idx, rec in enumerate(records):
        missing = [k for k in REQUIRED_RECORD_KEYS if rec.get(k) in (None, "")]
        assert not missing, (
            f"{_anchor(path, idx, rec)}: missing required key(s) {missing}"
        )


def validate_records(
    records, path, workflow_mode, sample_id, configured_assay_types, reference_version
):
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
        reference_version: Canonical build to keep; records of any other are dropped.

    Raises:
        AssertionError: Any record violates the spec, or the selected records violate
            a mode/replicate rule.
    """
    require_record_keys(records, path)
    single_cell = workflow_mode in ("single_cell_genotyping", "copytyping_preprocess")

    for idx, rec in enumerate(records):
        at = _anchor(path, idx, rec)
        assay_type = rec["assay_type"]
        assert assay_type in ALLOWED_ASSAY_TYPES, (
            f"{at}: assay_type must be one of {sorted(ALLOWED_ASSAY_TYPES)}"
        )
        assert rec["sample_type"] in ("normal", "tumor"), (
            f"{at}: sample_type must be 'normal' or 'tumor'"
        )
        did = rec["dataset_id"]
        assert did and all(c.isalnum() or c in "_-" for c in did), (
            f"{at}: dataset_id must match [A-Za-z0-9_-], got {did!r}"
        )
        assert isinstance(rec["files"], dict), f"{at}: files must be an object"

        readable = REQUIRED_FILES[assay_type] | OPTIONAL_FILES.get(assay_type, set())
        files = {k: v for k, v in rec["files"].items() if k in readable}
        ignored = set(rec["files"]) - set(files)
        if ignored:
            logging_snakemake(
                f"NOTE: {at}: ignoring files key(s) {sorted(ignored)}; "
                f"{assay_type} reads {sorted(readable)}",
            )
        rec["files"] = files

        required_files = REQUIRED_FILES[assay_type] if single_cell else ALIGNMENT_FILES
        for key in sorted(required_files):
            assert files.get(key), f"{at}: files.{key} is required for {assay_type}"

    selected = select_records(
        records, sample_id, configured_assay_types, reference_version
    )
    assert selected, (
        f"{path}: no records for sample_id={sample_id!r} with assay_type in "
        f"{sorted(configured_assay_types)} and reference_version={reference_version!r}; "
        f"sample_id has: {_refver_counts(records, sample_id)}"
    )

    seen = set()
    for idx, rec in enumerate(selected):
        key = (rec["dataset_id"], rec["assay_type"])
        assert key not in seen, (
            f"{_anchor(path, idx, rec)}: duplicate (dataset_id, assay_type) {key}"
        )
        seen.add(key)

    dataset2assays = {}
    for rec in selected:
        dataset2assays.setdefault(rec["dataset_id"], []).append(rec["assay_type"])
    for dataset_id, assays in dataset2assays.items():
        if len(assays) == 1:
            continue
        assert not any(a in BULK_ASSAYS for a in assays), (
            f"{path}: dataset_id={dataset_id!r} is reused across bulk assays {assays}"
        )
        assert len(assays) == 2 and set(assays) == {"scRNA", "scATAC"}, (
            f"{path}: dataset_id={dataset_id!r} has assays {assays}; only an "
            "scRNA + scATAC pair may share one"
        )

    dataset_ids = {r["dataset_id"] for r in selected}
    for idx, rec in enumerate(selected):
        base = rec.get("rdr_base_dataset_id")
        if not base:
            continue
        at = _anchor(path, idx, rec)
        assert rec["sample_type"] == "tumor", (
            f"{at}: rdr_base_dataset_id is set on a non-tumor record"
        )
        assert base != rec["dataset_id"], (
            f"{at}: rdr_base_dataset_id={base!r} is the record itself"
        )
        assert base in dataset_ids, (
            f"{at}: rdr_base_dataset_id={base!r} is not a dataset_id of {sample_id!r}"
        )


def parse_workflow(config):
    """Parse and validate the config + sample file into the workflow's globals.

    Args:
        config: The Snakemake config dict.

    Returns:
        Dict of the names workflow/Snakefile unpacks and the rules then read:
          workflow_mode, sample_id, remote_mode, reference_version, species,
          assay_types, modalities, msr_list, phaser,
          run_genotyping, run_phasing, het_snp_vcf, phased_snp_vcf,
          require_genetic_map, final_targets, get_data, modality2files,
          assay2dataset_ids, assay2sample_types, assay2base_reps, genotype_files,
          phase_files, get_genetic_map, get_phasing_panel, segment_bed, bedpe_files,
          has_breakpoints, use_prebuilt_windows, do_repliseq, window_size.

    Raises:
        AssertionError: The mode, assay types, sample file, or phaser is invalid.
    """

    def select_datasets(records, dataset_ids, config_key):
        """Records named by a config dataset_id list, in that order."""
        assert len(set(dataset_ids)) == len(dataset_ids), (
            f"{config_key} has duplicate dataset_ids: {dataset_ids}"
        )
        by_id = {r["dataset_id"]: r for r in records}
        missing = [d for d in dataset_ids if d not in by_id]
        assert not missing, (
            f"{config_key}={missing} not a dataset_id of sample_id={sample_id!r}; "
            f"available: {sorted(by_id)}"
        )
        return [by_id[d] for d in dataset_ids]

    # === workflow mode + sample_id ===
    sample_file = config["sample_file"]
    workflow_mode = config["workflow_mode"]
    assert workflow_mode in WORKFLOW_MODES, (
        f"workflow_mode must be one of {list(WORKFLOW_MODES)}"
    )
    sample_id = config["sample_id"]

    # === remote input mode: whole-file storage() download vs direct URL streaming ===
    remote_mode = config["remote_mode"]
    assert remote_mode in ("storage", "stream"), (
        f"remote_mode must be 'storage' or 'stream', got {remote_mode!r}"
    )
    assert remote_mode == "storage" or workflow_mode == "bulk_genotyping", (
        "remote_mode='stream' is only supported for bulk_genotyping"
    )

    # === assay_types requested: validate against the schema, keep this mode's ===
    config_assay_types = config["assay_types"]
    invalid = [a for a in config_assay_types if a not in ALLOWED_ASSAY_TYPES]
    assert not invalid, (
        f"invalid assay_types={invalid}; allowed: {sorted(ALLOWED_ASSAY_TYPES)}"
    )
    allowed = BULK_ASSAYS if workflow_mode == "bulk_genotyping" else NONBULK_ASSAYS
    config_assay_types = [a for a in config_assay_types if a in allowed]

    # === load the sample file (json or legacy tsv) ===
    ext = os.path.splitext(sample_file)[1].lower()
    assert ext in (".json", ".tsv", ".txt"), (
        f"{sample_file}: sample file must be .json or .tsv, got {ext!r}"
    )
    records = (
        parse_sample_file_json(sample_file)
        if ext == ".json"
        else parse_sample_file_tsv(sample_file)
    )

    require_record_keys(records, sample_file)

    # === chromosomes: must exist in genome_size ===
    genome_size = config["genome_size"]
    assert genome_size, "genome_size is required (two-column chrom<TAB>size file)"
    by_core = {}
    for name in get_chr_sizes(genome_size):
        by_core.setdefault(strip_chr_prefix(name), name)
    wanted = [strip_chr_prefix(c) for c in config["chromosomes"]]
    assert wanted, "chromosomes is empty"
    absent_chroms = [c for c in wanted if c not in by_core]
    assert not absent_chroms, (
        f"chromosomes {absent_chroms} have no contig in {genome_size}"
    )
    chroms = [f"chr{c}" for c in wanted]
    input_nochr = not by_core[wanted[0]].lower().startswith("chr")
    logging_snakemake(f"chromosomes: {chroms[:3]}... input_nochr={input_nochr}")

    # === species ===
    species = config["species"]
    assert species, f"species is required in the config; one of {list(SPECIES)}"
    if species not in SPECIES:
        logging_snakemake(
            f"WARNING: species={species!r} is not natively supported ({list(SPECIES)})."
        )

    # === reference version: canonicalize, then filter ===
    raw_refver = config["reference_version"]
    assert raw_refver, (
        f"reference_version is required in the config; one of {REFVERS} (or an alias)"
    )
    reference_version = canonical_refver(raw_refver)
    if not is_known_refver(raw_refver):
        logging_snakemake(
            f"WARNING: reference_version={raw_refver!r} is not natively supported "
            f"({REFVERS})."
        )
    matched = {}
    for rec in records:
        if canonical_refver(rec["reference_version"]) != reference_version:
            continue
        spelling = rec["reference_version"]
        n, ids = matched.get(spelling, (0, set()))
        matched[spelling] = (n + 1, ids | {rec["sample_id"]})
    logging_snakemake(
        f"reference_version: config={raw_refver!r} -> {reference_version}"
    )
    for spelling, (n_rec, ids) in sorted(matched.items()):
        logging_snakemake(
            f"  {spelling:<20} {n_rec:5d} record(s) {len(ids):4d} sample_id(s)"
        )

    # === assay_types supplied: this sample_id's records on this build ===
    sample_sheet_assay_types = list(
        dict.fromkeys(
            r["assay_type"]
            for r in records
            if r["sample_id"] == sample_id
            and canonical_refver(r["reference_version"]) == reference_version
        )
    )

    # === validate against the spec + selection rules (mutates files in place) ===
    validate_records(
        records,
        sample_file,
        workflow_mode,
        sample_id,
        config_assay_types,
        reference_version,
    )

    # === select this run's records + add modality ===
    records = [
        {**rec, "modality": ASSAY_TYPE2MODALITY[rec["assay_type"]]}
        for rec in select_records(
            records, sample_id, config_assay_types, reference_version
        )
    ]

    # === assay_types this run processes: requested and supplied ===
    assay_types = list(dict.fromkeys(r["assay_type"] for r in records))
    logging_snakemake(
        f"assay_types: config={config_assay_types} "
        f"sample_sheet={sample_sheet_assay_types} -> run={assay_types}"
    )

    # === RDR normalization policy: drop/keep each bulk tumor's rdr_base ===
    rdr_normalization = config["params_combine_counts"]["rdr_normalization"]
    assert rdr_normalization in RDR_NORMALIZATIONS, (
        f"rdr_normalization must be one of {list(RDR_NORMALIZATIONS)}, "
        f"got {rdr_normalization!r}"
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
    assert not (rdr_normalization == "normal" and unbased), (
        f"rdr_normalization='normal' requires rdr_base_dataset_id on every bulk "
        f"tumor, missing on: {sorted(unbased)}"
    )
    if ignored:
        logging_snakemake(
            f"NOTE: rdr_normalization='median' -> ignoring rdr_base_dataset_id on "
            f"{len(ignored)} tumor(s): {sorted(ignored)}"
        )
    if unbased:
        logging_snakemake(
            f"NOTE: {len(unbased)} tumor(s) have no rdr_base_dataset_id; RDR uses "
            f"median normalization: {sorted(unbased)}"
        )

    # === copytyping_preprocess requirements ===
    if workflow_mode == "copytyping_preprocess":
        assert config["het_snp_vcf"] is not None, (
            "het_snp_vcf is required for copytyping_preprocess"
        )
        assert config["het_snp_vcf_phased"], (
            "het_snp_vcf must be phased for copytyping_preprocess"
        )
        assert config["bb_file"] is not None, (
            "bb_file is required for copytyping_preprocess"
        )

    # === gtf_file is a required reference input (gene/exon annotation) ===
    assert config["gtf_file"], "gtf_file is required (gene/exon annotation GTF)"

    # === min_snp_reads sweep (one MSR{msr}/ subdir per value) ===
    msr = config["params_combine_counts"]["min_snp_reads"]
    msr_list = [int(m) for m in (msr if isinstance(msr, list) else [msr])]

    # === segment BED + bin BED build (one bin set for every bulk assay: WGS/WGS-lr/WES) ===
    bedpe_files = list(
        dict.fromkeys(
            r["files"]["breakpoint_bedpe"]
            for r in records
            if "breakpoint_bedpe" in r["files"]
        )
    )
    has_breakpoints = len(bedpe_files) > 0
    is_bulk = workflow_mode == "bulk_genotyping"
    do_repliseq = reference_version in REPLISEQ_REFVERS
    window_size = int(config["params_build_windows"]["window_size"])

    # skip window build if pre-built window bed is provided & no breakpoints
    window_bed = config["window_bed"]
    if window_bed is not None and not is_url(window_bed):
        assert os.path.exists(window_bed), f"window_bed does not exist: {window_bed}"
    use_prebuilt_windows = is_bulk and window_bed is not None and not has_breakpoints
    if is_bulk and window_bed is not None and has_breakpoints:
        logging_snakemake(
            f"NOTE: window_bed ignored ({window_bed}); {len(bedpe_files)} "
            "breakpoint_bedpe file(s) re-tile the arms, so windows are built"
        )
    segment_bed = (
        config["aux_dir"] + "/segment.bed" if is_bulk else config["region_bed"]
    )
    if is_bulk:
        logging_snakemake(
            f"NOTE: bulk -> segment BED ({segment_bed}); "
            f"{len(bedpe_files)} breakpoint_bedpe file(s), region_id=arm, seg_id=chunk; "
            f"windows: {'pre-built ' + window_bed if use_prebuilt_windows else 'built'}"
        )

    # === genotyping / phasing switches (a het_snp_vcf short-circuits the front) ===
    het_snp_vcf = config["het_snp_vcf"]
    run_genotyping = het_snp_vcf is None
    run_phasing = True
    phased_snp_vcf = config["phase_dir"] + "/phased_het_snps.vcf.gz"
    if het_snp_vcf is not None:
        assert os.path.exists(het_snp_vcf), f"het_snp_vcf does not exist: {het_snp_vcf}"
        # het_snp_vcf_phased=true -> the VCF is taken as phased, phasing skipped
        run_phasing = not bool(config["het_snp_vcf_phased"])
        if not run_phasing:
            phased_snp_vcf = het_snp_vcf

    # === phaser reference inputs (genotype + phase record selection) ===
    phaser = config["phaser"]
    genotype_files = None
    phase_files = None
    get_genetic_map = None
    get_phasing_panel = None
    if run_genotyping:
        named = config["genotype_dataset_ids"]
        if named:
            chosen = select_datasets(records, named, "genotype_dataset_ids")
            non_normal = [
                r["dataset_id"] for r in chosen if r["sample_type"] != "normal"
            ]
            if non_normal:
                logging_snakemake(
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
            assert chosen, f"no records to genotype for sample_id={sample_id!r}"
            r0 = chosen[0]
            logging_snakemake(
                f"NOTE: genotype_dataset_ids unset; genotyping {r0['dataset_id']!r} "
                f"({r0['sample_type']}, {r0['assay_type']})"
            )
        genotype_files = [r["files"] for r in chosen]
    if run_phasing:
        assert phaser in PANEL_PHASER or phaser in LONGREAD_PHASER, (
            f"unknown phaser: {phaser}"
        )
        if phaser in PANEL_PHASER:
            gmap_path = config["gmap_path"]
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
        else:
            named = config["phase_dataset_ids"]
            if named:
                chosen = select_datasets(records, named, "phase_dataset_ids")
                short = [
                    r["dataset_id"]
                    for r in chosen
                    if r["assay_type"] not in LONGREAD_ASSAYS
                ]
                assert not short, (
                    f"phase_dataset_ids={short} are not long-read assays "
                    f"({sorted(LONGREAD_ASSAYS)})"
                )
            else:
                lr = [r for r in records if r["assay_type"] in LONGREAD_ASSAYS]
                assert lr, (
                    f"{phaser} requires a long-read assay "
                    f"({sorted(LONGREAD_ASSAYS)}) for sample_id={sample_id!r}"
                )
                normals = [r for r in lr if r["sample_type"] == "normal"]
                chosen = normals or lr
                kind = "normal" if normals else "tumor (no normal)"
                logging_snakemake(
                    f"NOTE: phase_dataset_ids unset; longphase co-phases "
                    f"{[r['dataset_id'] for r in chosen]} ({kind})"
                )
            phase_files = [r["files"] for r in chosen]

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
        qc_genotype = config["qc_genotype_snps"]
        if isinstance(qc_genotype, str):  # --config passes bools as strings
            qc_genotype = qc_genotype.strip().lower() in ("true", "1", "yes")
        if run_genotyping and qc_genotype:
            final_targets.append(config["qc_dir"] + "/genotype_snp_qc.pdf")
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
        "remote_mode": remote_mode,
        "reference_version": reference_version,
        "species": species,
        "chroms": chroms,
        "input_nochr": input_nochr,
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
        "has_breakpoints": has_breakpoints,
        "use_prebuilt_windows": use_prebuilt_windows,
        "do_repliseq": do_repliseq,
        "window_size": window_size,
    }
