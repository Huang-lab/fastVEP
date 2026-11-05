#!/usr/bin/env python3
"""Cut a fixture of real ClinGen Gene-Disease Validity curations for the classifier's parsing test.

    make_clingen_parse_fixture.py <clingen_gene_validity.csv> <out_dir>

Writes into <out_dir> (`crates/fastvep-cli/fixtures/clingen/`):

  genemap.tsv   what `clingen_gdv_to_oga.py` makes of the chosen genes' curations, produced by that
                script and not re-implemented here: the file `fastvep sa-build --source omim` reads
  expected.tsv  what the classifier must read from each gene, worked out from the CSV rows alone,
                by the rules the criteria follow, and not from `genemap.tsv`

The test (`crates/fastvep-cli/tests/clingen_real_rows.rs`) builds an `.oga` from `genemap.tsv` with the
CLI's own `sa-build`, reads each gene back, and compares what `OmimData` makes of it with
`expected.tsv`. So a change to the converter's output format, to the `.oga` builder (which splits a
gene's phenotypes on `;`) or to the typed reading of the tag shows up as a gene that no longer reads as
the curation says it should.

The genes are chosen to cover what makes the reading hard, all of it taken from the real file: every
mode of inheritance the file uses (AD, AR, XL, SD, MT, UD), a gene with an autosomal dominant and an
autosomal recessive disease, a gene whose curations are all weak, the three labels that contain a `;`
(which the `.oga` builder splits into pieces, one without the tag), the three X-linked disorders whose
label says "dominant", and labels with parentheses, commas and slashes. A fixed number of genes of
each plain kind fill it out, taken in file order so that the selection is reproducible.

The rules, which are the criteria's and are repeated here on purpose so that the expectation is
independent of the Rust:
  established   the curation's classification is Definitive, Strong or Moderate
  dominant      an established curation tagged AD or SD
  recessive     an established curation tagged AR or SD
  no valid      the gene has curations and none is established
  contexts      the established curations grouped by tag, in the order AD AR XL SD MT UD
"""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from collections import OrderedDict, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ESTABLISHED = {"Definitive", "Strong", "Moderate"}
ORDER = ["AD", "AR", "XL", "SD", "MT", "UD"]
PLAIN_PER_KIND = 8


def read(path: Path) -> list[dict]:
    rows = []
    with path.open(newline="") as f:
        for r in csv.reader(f):
            if len(r) < 9 or r[0] == "GENE SYMBOL" or r[0].startswith(("+", "CLINGEN", "FILE", "WEBPAGE")):
                continue
            rows.append({"gene": r[0], "label": r[2], "mondo": r[3], "moi": r[4], "cls": r[6], "raw": r})
    return rows


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, out = Path(sys.argv[1]), Path(sys.argv[2])
    rows = read(src)
    by_gene: "OrderedDict[str, list[dict]]" = OrderedDict()
    for r in rows:
        by_gene.setdefault(r["gene"], []).append(r)

    chosen: "OrderedDict[str, None]" = OrderedDict()

    def take(genes):
        for g in genes:
            chosen.setdefault(g, None)

    est = lambda rs: [r for r in rs if r["cls"] in ESTABLISHED]
    moi_set = lambda rs: {r["moi"] for r in est(rs)}

    # the hard cases, all of them by rule from the file
    take(g for g, rs in by_gene.items() if any(";" in r["label"] for r in rs))
    take(g for g, rs in by_gene.items() if any(r["moi"] == "UD" for r in rs))
    take(g for g, rs in by_gene.items() if any(r["moi"] == "XL" and "dominant" in r["label"].lower() for r in rs))
    take(g for g, rs in by_gene.items() if any(r["moi"] == "SD" for r in rs))
    take(g for g, rs in by_gene.items() if any(r["moi"] == "MT" for r in rs))
    punct = [g for g, rs in by_gene.items() if any(("(" in r["label"] or "/" in r["label"]) for r in rs)]
    take(punct[:PLAIN_PER_KIND])
    dual = [g for g, rs in by_gene.items() if {"AD", "AR"} <= moi_set(rs)]
    take(dual[:PLAIN_PER_KIND])
    weak_only = [g for g, rs in by_gene.items() if rs and not est(rs)]
    take(weak_only[:PLAIN_PER_KIND])
    take([g for g, rs in by_gene.items() if moi_set(rs) == {"AD"}][:PLAIN_PER_KIND])
    take([g for g, rs in by_gene.items() if moi_set(rs) == {"AR"}][:PLAIN_PER_KIND])
    take([g for g, rs in by_gene.items() if moi_set(rs) == {"XL"}][:PLAIN_PER_KIND])
    # the genes the manuscript and the review name
    take(g for g in ("PAH", "KMT2D", "BRCA1", "ABCC8", "TTR", "NEK1", "DDX41", "SRCAP", "PROKR2", "VWF") if g in by_gene)

    # the converter, run on exactly these curations
    header = [
        ["CLINGEN GENE DISEASE VALIDITY CURATIONS"] + [""] * 9,
        ["GENE SYMBOL", "GENE ID (HGNC)", "DISEASE LABEL", "DISEASE ID (MONDO)", "MOI", "SOP", "CLASSIFICATION", "ONLINE REPORT", "CLASSIFICATION DATE", "GCEP"],
    ]
    with tempfile.TemporaryDirectory() as tmp:
        subset = Path(tmp) / "subset.csv"
        with subset.open("w", newline="") as f:
            w = csv.writer(f, quoting=csv.QUOTE_ALL)
            w.writerows(header)
            for g in chosen:
                for r in by_gene[g]:
                    w.writerow(r["raw"])
        out.mkdir(parents=True, exist_ok=True)
        subprocess.run([sys.executable, str(HERE / "clingen_gdv_to_oga.py"), str(subset), str(out / "genemap.tsv")], check=True)

    with (out / "expected.tsv").open("w") as f:
        f.write("gene\testablished\tno_valid\tdominant\trecessive\tcontexts\n")
        for g in sorted(chosen):
            rs = by_gene[g]
            e = est(rs)
            contexts: dict[str, int] = defaultdict(int)
            for r in e:
                contexts[r["moi"]] += 1
            dominant = any(r["moi"] in ("AD", "SD") for r in e)
            recessive = any(r["moi"] in ("AR", "SD") for r in e)
            ctx = "|".join(f"{m}:{contexts[m]}" for m in ORDER if m in contexts)
            f.write("\t".join([g, str(bool(e)).lower(), str(bool(rs) and not e).lower(), str(dominant).lower(), str(recessive).lower(), ctx]) + "\n")
    print(f"{len(chosen)} genes, {sum(len(by_gene[g]) for g in chosen)} curations -> {out}")


if __name__ == "__main__":
    main()
