def file_input(paths):
    """Rule input for sample-file path(s); accepts one path or a list.

    Local paths pass through unchanged. URLs are wrapped in storage(), which
    downloads the file once into .snakemake/storage/ and deletes the local copy
    once no remaining job needs it (default --keep-storage-local-copies=False).
    """
    if isinstance(paths, (list, tuple)):
        return [file_input(p) for p in paths]
    return storage(str(paths)) if is_url(paths) else str(paths)


# Remote bulk alignments are region-subset to `chromosomes` once (subset_remote_alignment),
# so every consuming rule reads a small local BAM instead of storage()-ing the whole file.
# Maps the remote alignment URL -> its local subset BAM path; local or single-cell inputs
# are absent and keep the storage()/local path.
remote_bulk_subset = {
    files["alignment"]: config["aux_dir"]
    + f"/remote_subset/{assay_type}_{dataset_id}.bam"
    for (assay_type, dataset_id), files in get_data.items()
    if assay_type in BULK_ASSAYS and is_url(files["alignment"])
}


def alignment_input(files):
    """Rule input for the .bam/.cram of one files map, or a list of files maps.

    A remote bulk alignment is redirected to its local region-subset BAM
    (subset_remote_alignment); any other URL falls back to whole-file storage().
    """
    if isinstance(files, (list, tuple)):
        return [alignment_input(f) for f in files]
    subset = remote_bulk_subset.get(files["alignment"])
    return subset if subset else file_input(files["alignment"])


def alignment_index_input(files):
    """Rule input for the alignment index of one files map, or a list of files maps.

    For a remote bulk alignment the index is the subset BAM's local ``.bai``; any
    other URL falls back to whole-file storage() of the provided index.
    """
    if isinstance(files, (list, tuple)):
        return [alignment_index_input(f) for f in files]
    subset = remote_bulk_subset.get(files["alignment"])
    return subset + ".bai" if subset else file_input(files["alignment_index"])


def download_slots(files):
    """Storage-retrieval cost of a job: 1 when any of its files is a URL, else 0.

    Gated by the `downloads` resource (profile/config.yaml), so at most N jobs that
    must fetch remote inputs run at once. That caps both concurrent HTTP requests
    (NCBI returns 503 under heavy parallelism) and peak local disk, which holds
    roughly `downloads` x the largest alignment.
    """
    if isinstance(files, (list, tuple)):
        return int(any(download_slots(f) for f in files))
    return int(any(is_url(v) for v in files.values()))


def spatial_layout(assay_type, dataset_ids):
    """Space Ranger `spatial/` members of each rep, for staging before squidpy.

    squidpy.read.visium() takes a directory, so process_rna_anndata stages each
    rep's files under <tmp>/spatial/ using their canonical Space Ranger names
    (RANGER_LAYOUT, const.py). Returns (names, paths): names[i] lists the spatial/
    filenames of rep i, and paths is those files flattened in the same order, so a
    rule can pass paths as an input list and names as a param and still pair them up.
    """
    names, paths = [], []
    for dataset_id in dataset_ids:
        files = get_data[(assay_type, dataset_id)]
        layout = {
            ranger_names[0]: files[key]
            for key, (ranger_names, in_spatial) in RANGER_LAYOUT.items()
            if in_spatial and files.get(key)
        }
        names.append(list(layout.keys()))
        paths.extend(layout.values())
    return names, paths


def cli_flag(params_dict, key, flag_name, is_bool=False):
    """Return CLI flag list for a nullable config param, or [] if null.

    For non-bool params: returns [flag_name, str(value)] if value is not None.
    For bool params: returns [flag_name] if value is truthy, [] otherwise.
    """
    val = params_dict.get(key)
    if val is None:
        return []
    if is_bool:
        return [flag_name] if val else []
    return [flag_name, str(val)]


def cli_flags_str(params_dict, *specs):
    """Build a shell-ready string of optional CLI flags.

    Each spec is (key, flag_name) or (key, flag_name, True) for bool flags.
    """
    parts = []
    for spec in specs:
        is_bool = len(spec) > 2 and spec[2]
        parts.extend(cli_flag(params_dict, spec[0], spec[1], is_bool=is_bool))
    return " ".join(parts)


# When a bulk alignment is a remote URL, fetch only `chromosomes` once into a local BAM
# via the remote index (htslib HTTP range reads), so genotype/phase/pileup/mosdepth all
# read the subset rather than downloading the whole file. BAM or CRAM in (the reference
# decodes CRAM, ignored for BAM); output is always BAM. Single-cell not covered.
if remote_bulk_subset:

    rule subset_remote_alignment:
        input:
            reference=config["reference"],
        output:
            bam=temp(config["aux_dir"] + "/remote_subset/{assay_type}_{dataset_id}.bam"),
            bai=temp(
                config["aux_dir"] + "/remote_subset/{assay_type}_{dataset_id}.bam.bai"
            ),
        log:
            config["log_dir"]
            + f"/subset_remote_alignment/subset_remote_alignment.{{assay_type}}_{{dataset_id}}.{_run_id}.log",
        benchmark:
            config["bench_dir"]
            + f"/subset_remote_alignment/subset_remote_alignment.{{assay_type}}_{{dataset_id}}.{_run_id}.tsv"
        wildcard_constraints:
            assay_type="(bulkWGS|bulkWGS-lr|bulkWES)",
        conda:
            "../envs/samtools.yaml"
        threads: config["threads"]["pileup"]
        resources:
            downloads=1,
        params:
            aln_url=lambda wc: get_data[(wc.assay_type, wc.dataset_id)]["alignment"],
            idx_url=lambda wc: get_data[(wc.assay_type, wc.dataset_id)][
                "alignment_index"
            ],
            regions=" ".join(f"chr{c}" for c in config["chromosomes"]),
        shell:
            r"""
            samtools view -b -@ {threads} -T "{input.reference}" -o "{output.bam}" \
                -X "{params.aln_url}" "{params.idx_url}" {params.regions} 2> {log}
            samtools index -@ {threads} "{output.bam}" 2>> {log}
            """
