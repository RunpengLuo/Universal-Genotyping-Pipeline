"""Input helpers shared by every rule file.

Last update: 2026-08-11

Functions:
- file_input: wrap a sample-file path in storage() when remote
- alignment_input, alignment_index_input: the .bam/.cram pair of one record
- bam_stream_input, bam_stream_index_input: the remote_mode stream variants
- bam_stream_arg: the htslib url##idx##idxurl argument for longphase
- download_slots: the downloads resource that throttles storage() retrieval
"""


def file_input(paths):
    """Rule input for sample-file path(s); accepts one path or a list.

    Local paths pass through unchanged. URLs are wrapped in storage(), which
    downloads the file once into .snakemake/storage/ and deletes the local copy
    once no remaining job needs it (default --keep-storage-local-copies=False).
    """
    if isinstance(paths, (list, tuple)):
        return [file_input(p) for p in paths]
    return storage(str(paths)) if is_url(paths) else str(paths)


def alignment_input(files):
    """Rule input for the .bam/.cram of one files map, or a list of files maps."""
    if isinstance(files, (list, tuple)):
        return [alignment_input(f) for f in files]
    return file_input(files["alignment"])


def alignment_index_input(files):
    """Rule input for the .bai/.crai of one files map, or a list of files maps."""
    if isinstance(files, (list, tuple)):
        return [alignment_index_input(f) for f in files]
    return file_input(files["alignment_index"])


def bam_stream_input(files):
    """Tracked `input:` for a bulk alignment, honoring `remote_mode`.

    Local paths pass through. A URL is wrapped in storage() (whole-file download) under
    `remote_mode: storage`, or dropped ([]) under `remote_mode: stream` so the URL is
    read directly by htslib via bam_stream_arg() instead of being staged.
    """
    if isinstance(files, (list, tuple)):
        out = [bam_stream_input(f) for f in files]
        return [x for x in out if x != []]
    aln = files["alignment"]
    if is_url(aln):
        return [] if remote_mode == "stream" else storage(aln)
    return str(aln)


def bam_stream_index_input(files):
    """Tracked `input:` for a bulk alignment index; [] for a streamed URL (see above)."""
    if isinstance(files, (list, tuple)):
        out = [bam_stream_index_input(f) for f in files]
        return [x for x in out if x != []]
    idx = files["alignment_index"]
    if is_url(idx):
        return [] if remote_mode == "stream" else storage(idx)
    return str(idx)


def bam_stream_arg(files):
    """Shell token for a streamed remote alignment (`url##idx##idxurl`), else ''.

    Non-empty only under `remote_mode: stream` for a URL alignment; rules use it as the
    fallback when the staged `{input.alignment}` is empty. For a list, joins the tokens.
    """
    if isinstance(files, (list, tuple)):
        return " ".join(filter(None, (bam_stream_arg(f) for f in files)))
    aln = files["alignment"]
    if remote_mode == "stream" and is_url(aln):
        return f"{aln}##idx##{files['alignment_index']}"
    return ""


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
    """Space Ranger `spatial/` members of each dataset_id, for staging before squidpy.

    squidpy.read.visium() takes a directory, so process_rna_anndata stages each
    dataset_id's files under <tmp>/spatial/ using their canonical Space Ranger names
    (RANGER_LAYOUT, const.py). Returns (names, paths): names[i] lists the spatial/
    filenames of dataset_id i, and paths is those files flattened in the same order, so a
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
