#!/usr/bin/env python3
"""Compare a rerun of the manuscript's runs with the runs it repeats, file by file and record by record.

A run directory holds the scored files (the concordance matrix, the firing rates, the rule distribution, the
discrepancies, the per-chromosome and per-consequence tables, the summary) and, for a run whose annotation the harness
keeps, the annotated VCF. Two runs agree when every scored file has the same content and, where both keep it, the VCF
has the same records. Not compared: the provenance file and the log, which differ by design (the time, the path and the
hash of the binary), and the header of a VCF, which names the version of the binary and the date. A scored file is
compared with its line ends normalised, because the scored files of the first runs of the manuscript end their lines
with CR LF and the scoring script now writes LF.

A per-variant table is compared by column name: the tables of the first runs were made by earlier versions of
`12_per_variant_table.py` that wrote fewer columns, so a table agrees when every column the earlier table has is the
same in every row, and the columns it lacks are reported.

    16_compare_runs.py PREV_ROOT NEW_ROOT [--names a,b] [--tables] [--pair PREV_DIR:NEW_DIR:LABEL ...] [--csv OUT]

PREV_ROOT and NEW_ROOT hold one directory per run, with the same names. `--pair` adds a run that is kept under
different names (`output_v37_prev` and `output_v37`), `--tables` adds the per-variant tables `<run>_per_variant.tsv.gz`
that both roots hold, and `--csv` writes one row per comparison for `manuscript/acmg/data/rerun_comparison.csv`, which
`verify_numbers.py` reads. The exit status is 1 when anything differs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import sys
from pathlib import Path

SKIP = {"run_provenance.txt", "run.log"}


def digest_text(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_vcf(path: Path) -> tuple[str, int]:
    h = hashlib.sha256()
    n = 0
    with gzip.open(path, "rb") as f:
        for line in f:
            if not line.startswith(b"#"):
                h.update(line)
                n += 1
    return h.hexdigest(), n


def compare_run(prev: Path, new: Path) -> tuple[list[str], int]:
    """The files that differ between two run directories, and how many VCFs were compared by record."""
    diffs, vcfs = [], 0
    for name in sorted({p.name for p in prev.iterdir()} | {p.name for p in new.iterdir()}):
        if name in SKIP:
            continue
        a, b = prev / name, new / name
        if name.endswith(".vcf.gz"):
            # A VCF kept by one run and not by the other is not a difference between the runs.
            if a.exists() and b.exists():
                vcfs += 1
                (da, na), (db, nb) = digest_vcf(a), digest_vcf(b)
                if da != db:
                    diffs.append(f"{name}: {na} and {nb} records, not the same")
        elif not a.exists():
            diffs.append(f"{name} only in the rerun")
        elif not b.exists():
            diffs.append(f"{name} only in the first run")
        elif digest_text(a) != digest_text(b):
            diffs.append(name)
    return diffs, vcfs


def compare_table(prev: Path, new: Path) -> tuple[list[str], list[str]]:
    """What differs between two per-variant tables by column name, and the columns the earlier table lacks."""
    problems: list[str] = []
    with gzip.open(prev, "rt") as fa, gzip.open(new, "rt") as fb:
        ha, hb = fa.readline().rstrip("\n").split("\t"), fb.readline().rstrip("\n").split("\t")
        missing = [c for c in ha if c not in hb]
        added = [c for c in hb if c not in ha]
        if missing:
            problems.append(f"the rerun lacks {missing}")
        shared = [c for c in ha if c in hb]
        ia, ib = [ha.index(c) for c in shared], [hb.index(c) for c in shared]
        rows = differing = 0
        for la, lb in zip(fa, fb):
            rows += 1
            a, b = la.rstrip("\n").split("\t"), lb.rstrip("\n").split("\t")
            if [a[i] for i in ia] != [b[i] for i in ib]:
                differing += 1
        if fa.readline() or fb.readline():
            problems.append("different numbers of rows")
        if differing:
            problems.append(f"{differing} of {rows} rows differ")
    return problems, added


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prev")
    ap.add_argument("new")
    ap.add_argument("--names", default="", help="comma-separated run names; the default is every directory of NEW_ROOT")
    ap.add_argument("--tables", action="store_true", help="also compare the per-variant tables both roots hold")
    ap.add_argument("--pair", action="append", default=[], help="PREV_DIR:NEW_DIR:LABEL for a run kept under two names")
    ap.add_argument("--csv", default="", help="write one row per comparison here")
    args = ap.parse_args()
    prev_root, new_root = Path(args.prev), Path(args.new)
    wanted = set(filter(None, args.names.split(",")))

    jobs = [(p.name, prev_root / p.name, p) for p in sorted(new_root.iterdir()) if p.is_dir() and (not wanted or p.name in wanted)]
    for spec in args.pair:
        a, b, label = spec.rsplit(":", 2)
        jobs.append((label, Path(a), Path(b)))

    rows: list[dict] = []
    for label, a, b in jobs:
        if not a.is_dir():
            rows.append({"kind": "run", "name": label, "identical": "no", "vcfs_compared": 0, "detail": "no first run to compare with"})
            continue
        diffs, vcfs = compare_run(a, b)
        rows.append({"kind": "run", "name": label, "identical": "no" if diffs else "yes", "vcfs_compared": vcfs, "detail": "; ".join(diffs)})
    if args.tables:
        for p in sorted(new_root.glob("*_per_variant.tsv.gz")):
            q = prev_root / p.name
            if wanted and p.name.removesuffix("_per_variant.tsv.gz") not in wanted:
                continue
            if not q.exists():
                continue
            problems, added = compare_table(q, p)
            rows.append({"kind": "table", "name": p.name.removesuffix("_per_variant.tsv.gz"), "identical": "no" if problems else "yes",
                         "vcfs_compared": 0, "detail": "; ".join(problems + ([f"columns the first table lacks: {added}"] if added else []))})

    runs = [r for r in rows if r["kind"] == "run"]
    tables = [r for r in rows if r["kind"] == "table"]
    print(f"{len(runs)} runs compared: {sum(r['identical'] == 'yes' for r in runs)} identical in every scored file "
          f"({sum(int(r['vcfs_compared']) for r in runs)} annotated VCFs compared by record)")
    if args.tables:
        print(f"{len(tables)} per-variant tables compared: {sum(r['identical'] == 'yes' for r in tables)} identical in every column the first table has")
    for r in rows:
        if r["identical"] != "yes":
            print(f"DIFFERENT {r['kind']} {r['name']}: {r['detail']}")
    if args.csv:
        with open(args.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["kind", "name", "identical", "vcfs_compared", "detail"], lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    return 1 if any(r["identical"] != "yes" for r in rows) else 0


if __name__ == "__main__":
    sys.exit(main())
