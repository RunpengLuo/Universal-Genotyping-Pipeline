#!/usr/bin/env python3
"""Unit tests for the sample-sheet loader's record checks.

Last update: 2026-08-13

Covers:
- ids: sample_id and dataset_id must match RECORD_ID_PATTERN, both encodings
- rejection: whitespace, path and shell separators, non-ASCII letters
- timing: the charset is checked at load, before any selection
"""

import json
import os
import sys

import pytest

_REPO = os.path.join(os.path.dirname(__file__), "..")
for _sub in ("config", "workflow/scripts", "workflow/scripts/script_utils"):
    sys.path.insert(0, os.path.join(_REPO, _sub))
const = pytest.importorskip("const")
parse_workflow_args = pytest.importorskip("parse_workflow_args")

RECORD = {
    "sample_id": "HT001",
    "dataset_id": "N1",
    "assay_type": "bulkWGS",
    "sample_type": "normal",
    "reference_version": "hg38",
    "files": {"alignment": "/d/n.bam", "alignment_index": "/d/n.bam.bai"},
}

BAD_IDS = ["N 1", "N\t1", "N.1", "a/b", "HT001;rm", "café", "N+1", ""]


def write_sheet(tmp_path, record, ext):
    """Write one record as a JSON or TSV sample file, return its path."""
    path = tmp_path / f"samples{ext}"
    if ext == ".json":
        path.write_text(json.dumps({"samples": [record]}))
        return str(path)
    cols = [k for k in record if k != "files"] + [
        f"{const.FILES_COLUMN_PREFIX}{k}" for k in record["files"]
    ]
    vals = [str(record[k]) for k in record if k != "files"] + list(
        record["files"].values()
    )
    path.write_text("\t".join(cols) + "\n" + "\t".join(vals) + "\n")
    return str(path)


@pytest.mark.parametrize("ext", [".json", ".tsv"])
@pytest.mark.parametrize("dataset_id", ["N1", "N_1", "N-1", "N_1-a", "0"])
def test_id_charset_accepts_alnum_underscore_dash(tmp_path, ext, dataset_id):
    """Letters, digits, underscore and dash pass in either encoding."""
    record = {**RECORD, "dataset_id": dataset_id, "sample_id": "HT_001-p2"}
    records = parse_workflow_args.read_sample_sheet(write_sheet(tmp_path, record, ext))
    assert records[0]["dataset_id"] == dataset_id


@pytest.mark.parametrize("key", const.RECORD_ID_KEYS)
@pytest.mark.parametrize("bad", BAD_IDS)
def test_id_charset_rejects_separators_and_non_ascii(tmp_path, key, bad):
    """A separator, whitespace or non-ASCII letter fails the load, naming the key."""
    record = {**RECORD, key: bad}
    with pytest.raises(AssertionError, match=key):
        parse_workflow_args.read_sample_sheet(write_sheet(tmp_path, record, ".json"))


def test_id_charset_is_checked_before_selection(tmp_path):
    """A record no run would select still fails: the check is at load, not selection."""
    record = {**RECORD, "sample_id": "OTHER PATIENT"}
    with pytest.raises(AssertionError, match="sample_id"):
        parse_workflow_args.read_sample_sheet(write_sheet(tmp_path, record, ".json"))
