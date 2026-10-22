#!/usr/bin/env python3
"""One row per ClinVar 2-star+ variant: ClinVar's classification, fastVEP's call and
the evidence behind it, read from an annotated run.

`03_evaluate_concordance.py` reduces a run to counts. Anything that needs to ask a
question of individual variants - which criteria sit behind a discordant call, how
the evidence points are distributed by ClinVar class, what a REVEL band holds -
needs the variants themselves. This writes them once, so every later analysis reads
one table instead of re-parsing a 250 MB bgzip.

The transcript each variant is scored on is chosen by `03`'s own `pick_csq`, imported
and not re-implemented, and the script stops if its calls and criteria disagree with
`03`'s on any variant, or if the resulting matrix is not the run's own
`concordance_matrix.csv`. A table that cannot reproduce the published matrix is not
a table of the published run.

ClinVar's 2-star+ classifications are a reference and not ground truth.

Usage
  12_per_variant_table.py --run data/benchmark/output_v37 --out data/benchmark/manuscript_runs/v37_per_variant.tsv.gz
"""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]

_spec = importlib.util.spec_from_file_location("evaluate_concordance", HERE / "03_evaluate_concordance.py")
ev = importlib.util.module_from_spec(_spec)
sys.modules["evaluate_concordance"] = ev
_spec.loader.exec_module(ev)

COLUMNS = [
    "chrom", "pos", "ref", "alt", "gene", "stars", "clinvar", "clnsig", "clnvc",
    "call", "criteria", "consequence", "impact", "biotype", "transcript", "mane_select",
    "hgvsc", "hgvsp", "revel", "spliceai_max", "phylop",
    "gnomad_af", "gnomad_ac", "gnomad_an", "gnomad_hc", "gnomad_popmax_af",
    "pli", "loeuf", "mis_z", "clingen_tags", "symbol",
]

POP_AF_FIELDS = range(5, 15)  # AFR..SAS in FV_GNOMAD, after ALLELE|ALL_AF|ALL_AC|ALL_AN|ALL_HC


def info_dict(info: str) -> dict[str, str]:
    out = {}
    for piece in info.split(";"):
        k, sep, v = piece.partition("=")
        if sep:
            out[k] = v
    return out


def fv_fields(info: dict[str, str], key: str, alleles: list[str]) -> list[str] | None:
    """The `|`-split fields of an FV_ annotation for one allele, or None.

    ClinVar's VCF is biallelic, so there is one entry; where there are several
    (comma-separated) the one whose first field is one of the allele's spellings is
    taken. The annotation names an allele the way VEP's CSQ does (`T` for an inserted
    base, `-` for a deletion), so that spelling is tried before the VCF's own.
    """
    raw = info.get(key)
    if not raw:
        return None
    entries = [e.split("|") for e in raw.split(",")]
    for allele in alleles:
        for parts in entries:
            if parts and parts[0] == allele:
                return parts
    return None


def num(x: str | None) -> str:
    if x is None or x == "" or x == ".":
        return ""
    try:
        float(x)
    except ValueError:
        return ""
    return x


def clingen_tags(omim: str) -> str:
    """The inheritance tags of a gene's established ClinGen curations, joined by `+`.

    Reads the `(ClinGen Definitive/AD, MONDO:...)` the annotation carries, counting
    only Definitive, Strong and Moderate curations, and returns `none` for a gene with
    no such curation. This describes the data for the manuscript's tables; the
    classifier does not read the tag (limitation 20 in docs/ACMG.md).
    """
    tags: set[str] = set()
    for item in omim.replace("%20", " ").replace("%2C", ",").split("&"):
        if "ClinGen" not in item:
            continue
        cls = item.split("ClinGen", 1)[1].strip().split("/", 1)
        if len(cls) < 2 or cls[0].strip() not in ("Definitive", "Strong", "Moderate"):
            continue
        tags.add(cls[1].split(",", 1)[0].split(")", 1)[0].strip())
    return "+".join(sorted(tags)) if tags else "none"


def pick_full(entries, csq_idx):
    """`03`'s pick_csq, returning the whole chosen CSQ row as well."""
    if not entries:
        return None, None
    acmg_i = csq_idx.get("ACMG", -1)
    can_i = csq_idx.get("CANONICAL", -1)

    def acmg_of(parts):
        return parts[acmg_i] if 0 <= acmg_i < len(parts) else ""

    populated = [p for p in entries if acmg_of(p)]
    pool = populated or entries
    canon = [p for p in pool if can_i >= 0 and len(p) > can_i and p[can_i] == "YES"]
    return (canon[0] if canon else pool[0]), bool(populated)


def records(vcf: str):
    """Yield (chrom, pos, ref, alt, csq_idx, entries_for_alt, info) per ALT allele."""
    csq_idx = None
    with ev.open_text(vcf) as f:
        for line in f:
            if line.startswith("##"):
                if "ID=CSQ" in line:
                    csq_idx = {n: i for i, n in enumerate(ev.parse_csq_format(line))}
                continue
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8 or csq_idx is None:
                continue
            chrom = cols[0].removeprefix("chr")
            info = info_dict(cols[7])
            csq = info.get("CSQ")
            if csq is None:
                continue
            ref = cols[3]
            by_allele = defaultdict(list)
            acmg_i = csq_idx.get("ACMG", -1)
            for ent in csq.split(","):
                parts = ent.split("|")
                if len(parts) <= acmg_i:
                    continue
                by_allele[parts[0]].append(parts)
            for alt in cols[4].split(","):
                yield chrom, cols[1], ref, alt, csq_idx, by_allele.get(ev.vep_allele(ref, alt), []), info


def build(run: Path, labels_path: Path, out: Path) -> None:
    labels = ev.load_labels(labels_path)
    vcf = str(run / "clinvar_2star.fastvep.vcf.gz")
    n = 0
    matrix = {t: {p: 0 for p in ev.CLASSES + ["NoCall"]} for t in ev.CLASSES}
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        for chrom, pos, ref, alt, csq_idx, entries, info in records(vcf):
            key = (chrom, pos, ref, alt)
            lab = labels.get(key)
            if lab is None:
                continue
            acmg, criteria, top, _ = ev.pick_csq(entries, csq_idx)
            row, populated = pick_full(entries, csq_idx)
            if row is not None:
                assert (row[csq_idx["ACMG"]] if populated else "") == (acmg or ""), key
            call = ev.class_label(acmg) if acmg else None
            matrix[lab["normalized_class"]][call or "NoCall"] += 1

            def f(name):
                i = csq_idx.get(name, -1)
                return row[i] if row is not None and 0 <= i < len(row) else ""

            alleles = [ev.vep_allele(ref, alt), alt]
            gn = fv_fields(info, "FV_GNOMAD", alleles)
            pop = max(
                (float(gn[i]) for i in POP_AF_FIELDS if gn and i < len(gn) and num(gn[i])),
                default=None,
            )
            sp = fv_fields(info, "SpliceAI", alleles)
            sp_max = ""
            if sp and len(sp) >= 6:
                vals = [float(v) for v in sp[2:6] if num(v)]
                sp_max = f"{max(vals):.2f}" if vals else ""
            rv = fv_fields(info, "FV_REVEL", alleles)
            ph = fv_fields(info, "FV_PHYLOP", alleles)
            gg = (info.get("FV_GNOMAD_GENE") or "").split("|")
            w.writerow(
                [
                    chrom, pos, ref, alt, lab["gene"], lab["review_stars"], lab["normalized_class"],
                    lab["clnsig"], info.get("CLNVC", ""),
                    call or "NoCall", criteria or "", top, f("IMPACT"), f("BIOTYPE"), f("Feature"),
                    f("MANE_SELECT"), f("HGVSc"), f("HGVSp"),
                    num(rv[1]) if rv and len(rv) > 1 else "",
                    sp_max,
                    num(ph[1]) if ph and len(ph) > 1 else "",
                    num(gn[1]) if gn and len(gn) > 1 else "",
                    num(gn[2]) if gn and len(gn) > 2 else "",
                    num(gn[3]) if gn and len(gn) > 3 else "",
                    num(gn[4]) if gn and len(gn) > 4 else "",
                    "" if pop is None else repr(pop),
                    num(gg[1]) if len(gg) > 1 else "",
                    num(gg[2]) if len(gg) > 2 else "",
                    num(gg[3]) if len(gg) > 3 else "",
                    clingen_tags(info.get("FV_OMIM", "")),
                    # The symbol of the picked transcript, which is not always ClinVar's first GENEINFO symbol
                    # (`gene`, the column the gene-level counts of the manuscript use).
                    f("SYMBOL"),
                ]
            )
            n += 1
            if n % 100000 == 0:
                print(f"  {n:,} variants", file=sys.stderr, flush=True)

    # The table must reproduce the run's own matrix, or it is not a table of that run.
    want = {}
    with (run / "concordance_matrix.csv").open() as f:
        rows = list(csv.reader(f))
    for r in rows[1:]:
        want[r[0]] = {c: int(v) for c, v in zip(rows[0][1:], r[1:])}
    if matrix != want:
        diff = [(t, p, matrix[t][p], want[t][p]) for t in matrix for p in matrix[t] if matrix[t][p] != want[t][p]]
        sys.exit(f"the table does not reproduce {run}/concordance_matrix.csv: {diff[:10]}")
    print(f"wrote {n:,} variants to {out}; matrix identical to {run.name}/concordance_matrix.csv")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, type=Path)
    ap.add_argument("--labels", type=Path, default=ROOT / "data/benchmark/clinvar_2star_labels.tsv")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    build(args.run, args.labels, args.out)


if __name__ == "__main__":
    main()
