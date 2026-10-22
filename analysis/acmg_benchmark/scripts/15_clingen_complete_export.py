#!/usr/bin/env python3
"""Put back the ClinGen curations that the 30 April 2026 export did not list.

The benchmark's gene-disease source is ClinGen's gene-disease validity export of
30 April 2026. Checking the PVS1 gap of the ACMG manuscript found that CFTR, which the
current export lists as Definitive (autosomal recessive, classified 1 June 2022), is not
in it, and neither are 103 other curations that ClinGen dates before that day. The export
is a download, and we did not establish why those rows were missing.

This builds the file the benchmark would have read had the export been complete: the
April rows, plus every row of a later export that is dated on or before the day of the
April export and whose (gene, disease) pair the April export did not list. Rows dated
after 30 April are left out, so nothing ClinGen learned later reaches the benchmark, and
a pair that both exports list keeps its April classification.

The rows it adds are not all curations the April export lacked, and the report says which
are not. A row is a re-label when its assertion identifier (the UUID in the report URL) is
one that the April export lists for the same gene under another disease identifier, which
is how ClinGen's later change of a MONDO term shows up: it adds nothing, and 7 of the 112
rows are of this kind. The one row dated 30 April itself (GATM, 04:00 UTC) is later than
the time at which the April file was downloaded (02:19 UTC, from the file's modification
time), and the manuscript checks that it changes no call. What remains, 104 rows in 100
genes, are curations dated before the day of the export that the export lacked.

Usage

    # 1. the current export (a public file; the date of the download is in its header)
    curl -L -o data/benchmark/sa_sources/clingen_gene_validity_20261011.csv \\
        https://search.clinicalgenome.org/kb/gene-validity/download

    # 2. the merged export, in the same layout as the original
    python3 15_clingen_complete_export.py \\
        data/benchmark/sa_sources/clingen_gene_validity.csv \\
        data/benchmark/sa_sources/clingen_gene_validity_20261011.csv \\
        data/benchmark/sa_sources/clingen_gene_validity_complete.csv

    # 3. the usual conversion, then the database 11_manuscript_runs.py reads
    python3 sa_sources/clingen_gdv_to_oga.py \\
        data/benchmark/sa_sources/clingen_gene_validity_complete.csv \\
        data/benchmark/sa_sources/clingen_gdv_complete.tsv
    mkdir -p data/benchmark/sa_replace
    target/release/fastvep sa-build --source omim \\
        -i data/benchmark/sa_sources/clingen_gdv_complete.tsv \\
        -o data/benchmark/sa_replace/omim_complete

It prints what it added and writes `manuscript/acmg/data/clingen_missing_curations.csv`, with the status of each added
row (lacked, relabel or export day).
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
EXPORT_DAY = "2026-04-30"
KEY = ("GENE SYMBOL", "DISEASE ID (MONDO)")
ESTABLISHED = {"Definitive", "Strong", "Moderate"}
# The report URL ends in the assertion identifier and, for the current form, the classification date;
# older curations carry a short number instead of a UUID (CGGCIEX:assertion_6238).
ASSERTION = re.compile(r"assertion_(.+?)(?:-\d{4}-\d{2}-\d{2}T\d+\.\d+Z)?$")


def read(path: Path) -> tuple[list[list[str]], list[str], list[list[str]]]:
    """The preamble lines, the column header, and the data rows of an export."""
    rows = list(csv.reader(path.open(encoding="utf-8")))
    h = next(i for i, r in enumerate(rows) if r and r[0] == "GENE SYMBOL")
    # The line after the header is a row of plus signs, which is part of the preamble.
    body = [r for r in rows[h + 1 :] if r and not r[0].startswith("+++")]
    return rows[: h + 1] + rows[h + 1 : h + 2], rows[h], body


def main() -> None:
    if len(sys.argv) != 4:
        sys.exit("usage: 15_clingen_complete_export.py <april export> <current export> <out csv>")
    april_path, current_path, out_path = map(Path, sys.argv[1:])
    pre, head, april = read(april_path)
    _, head2, current = read(current_path)
    if head != head2:
        sys.exit("the two exports do not have the same columns")
    ix = {name: head.index(name) for name in head}
    key = lambda r: tuple(r[ix[k]] for k in KEY)  # noqa: E731
    have = {key(r) for r in april}
    date = lambda r: r[ix["CLASSIFICATION DATE"]][:10]  # noqa: E731
    gene = lambda r: r[ix["GENE SYMBOL"]]  # noqa: E731

    def assertion(r: list[str]) -> str:
        m = ASSERTION.search(r[ix["ONLINE REPORT"]])
        if m is None:
            sys.exit(f"no assertion identifier in the report URL of {gene(r)}: {r[ix['ONLINE REPORT']]}")
        return m.group(1)

    april_assertions = {(gene(r), assertion(r)) for r in april}
    added = [r for r in current if key(r) not in have and date(r) <= EXPORT_DAY]

    def status(r: list[str]) -> str:
        if (gene(r), assertion(r)) in april_assertions:
            return "relabel"  # the April export lists this assertion under another disease identifier
        if date(r) == EXPORT_DAY:
            return "export day"  # classified on the day of the export, after the file was made
        return "lacked"

    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, quoting=csv.QUOTE_ALL, lineterminator="\n")
        w.writerows(pre)
        w.writerows(april)
        w.writerows(added)
    cols = ["GENE SYMBOL", "DISEASE LABEL", "MOI", "CLASSIFICATION", "CLASSIFICATION DATE"]
    report = ROOT / "manuscript/acmg/data/clingen_missing_curations.csv"
    with report.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["gene", "disease", "moi", "classification", "classified", "panel", "established", "status"])
        for r in sorted(added, key=lambda r: (gene(r), r[ix["DISEASE LABEL"]])):
            w.writerow([r[ix[c]] for c in cols] + [r[ix["GCEP"]], r[ix["CLASSIFICATION"]] in ESTABLISHED, status(r)])
    kinds = {k: [r for r in added if status(r) == k] for k in ("lacked", "relabel", "export day")}
    print(f"April export: {len(april)} rows; current export: {len(current)} rows")
    print(f"added {len(added)} rows dated on or before {EXPORT_DAY}: " + ", ".join(f"{len(v)} {k}" for k, v in kinds.items()))
    print(f"the {len(kinds['lacked'])} curations that the export lacked are in {len({gene(r) for r in kinds['lacked']})} genes")
    print(f"wrote {out_path} and {report}")


if __name__ == "__main__":
    main()
