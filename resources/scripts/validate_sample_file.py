#!/usr/bin/env python3
"""Validate a JSON (or legacy TSV) sample file without running the workflow.

Runpeng Luo (2026-07-12)

Runs the same parser and the same checks the Snakefile runs at DAG build, so a
malformed sample file is caught in a second rather than mid-run. Every sample_id
is validated in each workflow mode its assay types belong to, unless --sample-id
or --workflow-mode narrows it. Without --check-files, only the schema is checked.

--check-files resolves every input a rule would read: a local path is stat'ed, a
remote URL gets a 2-byte ranged GET (the same request Snakemake's HTTP storage
plugin issues), so a 404, a typo, or a host that forbids listing is caught before
the run. Nothing is downloaded. --skip-remote checks local paths only.

Only keys the record's assay type consumes are checked: validation drops the rest
(e.g. record-only fastq_r1), so they are never stat'ed or fetched.

Dependencies:
  Python 3 standard library; config/const.py and workflow/scripts/parse_workflow_args.py
  of this repo (located relative to this script).

Usage:
  python resources/scripts/validate_sample_file.py samples.json
  python resources/scripts/validate_sample_file.py samples.tsv --check-files
  python resources/scripts/validate_sample_file.py samples.json --sample-id HG008
      sample_file       # .json sample file, or a legacy .tsv sheet
      --sample-id       # only this sample_id (default: every one in the file)
      --workflow-mode   # only this mode (default: every mode the assays allow)
      --check-files     # resolve every path: stat local files, range-GET remote URLs
      --skip-remote     # with --check-files, do not touch the network
      --jobs            # concurrent URL checks (default 8)
      --timeout         # seconds per URL check (default 30)

Inputs
  sample file: docs/sample_sheet.md
Outputs:
  a per-(sample_id, mode) report on stdout; exit status 1 if anything failed
Notes/References:
  Schema and validation rules: docs/sample_sheet.md
"""

import argparse
import collections
import concurrent.futures
import os
import sys
import urllib.error
import urllib.request

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "config"))
sys.path.insert(0, os.path.join(REPO, "workflow", "scripts"))

from const import BULK_ASSAYS, NONBULK_ASSAYS, is_url  # noqa: E402
from parse_workflow_args import parse_sample_file, validate_records  # noqa: E402

# a mode can only run the assay types it supports
MODE_ASSAYS = {
    "bulk_genotyping": BULK_ASSAYS,
    "single_cell_genotyping": NONBULK_ASSAYS,
    "copytyping_preprocess": NONBULK_ASSAYS,
}


def check_url(url, timeout):
    """Resolve one URL with a 2-byte ranged GET; return None if OK, else the reason.

    A ranged GET is what Snakemake's HTTP storage plugin issues, and unlike HEAD it
    is answered by every host here. 200 means the server ignored the Range header and
    would have sent the whole file, which still proves the URL resolves.
    """
    req = urllib.request.Request(url, headers={"Range": "bytes=0-1"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            if resp.status in (200, 206):
                return None
            return f"HTTP {resp.status}"
    except urllib.error.HTTPError as exc:
        return f"HTTP {exc.code}"
    except Exception as exc:  # noqa: BLE001 - any failure to resolve is a failure
        return str(exc)


def check_files(records, skip_remote=False, jobs=8, timeout=30):
    """Resolve every input of every record: stat local paths, range-GET remote URLs.

    Args:
        records: Validated records; their files maps hold only the keys the assay reads.
        skip_remote: Do not touch the network; URLs are counted, not checked.
        jobs: Concurrent URL checks.
        timeout: Seconds per URL check.

    Returns:
        (bad, n_local, n_remote, n_skipped): the inputs that did not resolve, and how
        many local paths, URLs, and skipped URLs were seen.
    """
    bad, n_local, n_skipped = [], 0, 0
    urls = []
    for rec in records:
        for key, path in rec["files"].items():
            where = f"{rec['sample_id']}/{rec['dataset_id']}/{rec['assay_type']} files.{key}"
            if is_url(path):
                if skip_remote:
                    n_skipped += 1
                    continue
                urls.append((where, path))
            else:
                n_local += 1
                if not os.path.exists(path):
                    bad.append(f"{where}: {path} (no such file)")

    if urls:
        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
            reasons = pool.map(lambda u: check_url(u[1], timeout), urls)
            for (where, url), reason in zip(urls, reasons):
                if reason:
                    bad.append(f"{where}: {url} ({reason})")
    return bad, n_local, len(urls), n_skipped


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("sample_file")
    ap.add_argument("--sample-id", default=None)
    ap.add_argument("--workflow-mode", default=None, choices=sorted(MODE_ASSAYS))
    ap.add_argument("--check-files", action="store_true")
    ap.add_argument("--skip-remote", action="store_true")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=30)
    args = ap.parse_args()

    try:
        records = parse_sample_file(args.sample_file)
    except (ValueError, OSError) as e:
        print(f"FAIL parse: {e}")
        return 1

    sample_ids = sorted({r["sample_id"] for r in records})
    assays = collections.Counter(r["assay_type"] for r in records)
    print(f"{args.sample_file}")
    print(f"  {len(records)} records, {len(sample_ids)} sample_id(s)")
    print(f"  assays: {dict(assays)}")

    if args.sample_id:
        if args.sample_id not in sample_ids:
            print(f"FAIL: sample_id={args.sample_id!r} is not in the file")
            return 1
        sample_ids = [args.sample_id]
    modes = [args.workflow_mode] if args.workflow_mode else sorted(MODE_ASSAYS)

    ok = failed = 0
    for sample_id in sample_ids:
        for mode in modes:
            present = sorted(
                {
                    r["assay_type"]
                    for r in records
                    if r["sample_id"] == sample_id
                    and r["assay_type"] in MODE_ASSAYS[mode]
                }
            )
            if not present:
                continue
            try:
                validate_records(records, args.sample_file, mode, sample_id, present)
                ok += 1
            except ValueError as e:
                failed += 1
                print(f"  FAIL {sample_id} [{mode}]: {e}")

    print(f"\nvalidated: {ok} (sample_id, mode) combo(s) OK, {failed} failed")

    bad = []
    if args.check_files:
        scoped = [r for r in records if r["sample_id"] in set(sample_ids)]
        bad, n_local, n_remote, n_skipped = check_files(
            scoped,
            skip_remote=args.skip_remote,
            jobs=args.jobs,
            timeout=args.timeout,
        )
        print(
            f"\nfiles: {n_local} local stat'ed, {n_remote} remote resolved"
            + (f", {n_skipped} remote skipped" if n_skipped else "")
        )
        if bad:
            print(f"{len(bad)} unresolved:")
            for b in bad[:20]:
                print(f"  {b}")
            if len(bad) > 20:
                print(f"  ... and {len(bad) - 20} more")
        elif n_local or n_remote:
            print("  all inputs resolve")

    return 1 if (failed or bad) else 0


if __name__ == "__main__":
    sys.exit(main())
