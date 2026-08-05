#!/usr/bin/env python3
"""Remote sample-file inputs: URLs are wrapped in storage() and fetched with their index.

Runpeng Luo (2026-07-12)

The local test serves a stub alignment over http://127.0.0.1 and dry-runs against
it, so no network is needed. The GIAB test resolves real URLs and is marked
`network`; it is deselected unless `-m network` is given.

Dependencies:
  pytest; snakemake with snakemake-storage-plugin-http.

Usage:
  pytest tests/test_remote.py                # local http server only
  pytest tests/test_remote.py -m network     # also hit the GIAB FTP host

Notes/References:
  Remote inputs: docs/sample_sheet.md
"""

import functools
import http.server
import json
import os
import threading

import pytest

from conftest import dryrun

GIAB_BAM = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/"
    "Liss_lab/Element_AVITI_20240118/HG008-N-D_Element-StdInsert_61x_GRCh38-GIABv3.bam"
)


@pytest.fixture(scope="module")
def http_server(tmp_path_factory):
    """Serve a directory holding a stub alignment + index; yields its base URL."""
    served = tmp_path_factory.mktemp("served")
    for name in ("remote.bam", "remote.bam.bai"):
        (served / name).write_bytes(b"")

    handler = functools.partial(
        http.server.SimpleHTTPRequestHandler, directory=str(served)
    )
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_port}"
    server.shutdown()


def _remote_sheet(path, alignment, index):
    doc = {
        "version": 1,
        "samples": [
            {
                "sample_id": "T1",
                "dataset_id": "N1",
                "assay_type": "bulkWGS",
                "sample_type": "normal",
                "reference_version": "chm13v2",
                "files": {"alignment": alignment, "alignment_index": index},
            }
        ],
    }
    with open(path, "w") as fh:
        json.dump(doc, fh)
    return path


def test_local_url_is_retrieved_from_storage(workspace, http_server):
    """An http(s) alignment and its index are both planned as storage retrievals."""
    sheet = _remote_sheet(
        os.path.join(workspace["root"], "remote.json"),
        f"{http_server}/remote.bam",
        f"{http_server}/remote.bam.bai",
    )
    proc = dryrun(workspace, sheet, "T1", "bulk_genotyping", ["bulkWGS"])
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert proc.stdout.count("retrieve from storage") >= 2
    assert "remote.bam.bai" in proc.stdout


def test_local_and_remote_mix(workspace, http_server):
    """One sheet may mix a remote alignment with local ones."""
    ref = workspace["ref"]
    sheet = os.path.join(workspace["root"], "mixed.json")
    doc = {
        "version": 1,
        "samples": [
            {
                "sample_id": "T1",
                "dataset_id": "N1",
                "assay_type": "bulkWGS",
                "sample_type": "normal",
                "reference_version": "chm13v2",
                "files": {
                    "alignment": f"{http_server}/remote.bam",
                    "alignment_index": f"{http_server}/remote.bam.bai",
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
                    "alignment": f"{ref}/tumor.bam",
                    "alignment_index": f"{ref}/tumor.bam.bai",
                },
            },
        ],
    }
    with open(sheet, "w") as fh:
        json.dump(doc, fh)

    proc = dryrun(workspace, sheet, "T1", "bulk_genotyping", ["bulkWGS"])
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "retrieve from storage" in proc.stdout
    assert f"{ref}/tumor.bam" in proc.stdout


def test_tsv_sheet_accepts_urls(workspace, http_server):
    """The TSV encoding names every input, so it can carry remote URLs like the JSON."""
    sheet = os.path.join(workspace["root"], "remote.tsv")
    with open(sheet, "w") as fh:
        fh.write(
            "sample_id\tdataset_id\tassay_type\tsample_type\treference_version"
            "\tfiles.alignment\tfiles.alignment_index\n"
            f"T1\tN1\tbulkWGS\tnormal\tchm13v2"
            f"\t{http_server}/remote.bam\t{http_server}/remote.bam.bai\n"
        )
    proc = dryrun(workspace, sheet, "T1", "bulk_genotyping", ["bulkWGS"])
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert proc.stdout.count("retrieve from storage") >= 2


@pytest.mark.network
def test_giab_url_resolves(workspace):
    """A real GIAB alignment URL plans a storage retrieval (no download in -n)."""
    sheet = _remote_sheet(
        os.path.join(workspace["root"], "giab.json"), GIAB_BAM, GIAB_BAM + ".bai"
    )
    proc = dryrun(workspace, sheet, "T1", "bulk_genotyping", ["bulkWGS"])
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "retrieve from storage" in proc.stdout
