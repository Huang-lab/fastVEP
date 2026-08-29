#!/usr/bin/env python3
"""Diff two benchmark runs' ACMG calls, variant by variant.

The v23-to-v26 table was built by hand, which is why there is no script for it.
This is that table's generator, so the next run's delta is reproducible.

The picked transcript comes from `03_evaluate_concordance.py`'s own `pick_csq`,
imported rather than reimplemented: a call diff that picked a different
transcript than the concordance report would not describe the same run.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
EVAL = HERE / "03_evaluate_concordance.py"


def load_eval_module(path: Path):
    spec = importlib.util.spec_from_file_location("acmg_eval", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


FIELDS = [
    ("Consequence", "consequence"),
    ("IMPACT", "impact"),
    ("HGVSc", "hgvsc"),
    ("HGVSp", "hgvsp"),
    ("Feature", "transcript"),
    ("SYMBOL", "gene"),
]


def calls(ev, vcf_path):
    """key -> dict of the picked entry's fields, for every variant in the VCF."""
    out = {}
    for chrom, pos, ref, alt, idx, entries in ev.variant_records(vcf_path):
        acmg, crit, _top, _canon = ev.pick_csq(entries, idx)
        row = {"acmg": acmg or "", "criteria": crit or ""}
        chosen = None
        if entries:
            acmg_i = idx.get("ACMG", -1)
            populated = [p for p in entries if 0 <= acmg_i < len(p) and p[acmg_i]]
            pool = populated or entries
            can_i = idx.get("CANONICAL", -1)
            canon = [p for p in pool if can_i >= 0 and len(p) > can_i and p[can_i] == "YES"]
            chosen = canon[0] if canon else pool[0]
        for name, key in FIELDS:
            i = idx.get(name, -1)
            row[key] = chosen[i] if chosen and 0 <= i < len(chosen) else ""
        out[(chrom, pos, ref, alt)] = row
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--old-name", default="old")
    ap.add_argument("--new-name", default="new")
    ap.add_argument("--truth")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ev = load_eval_module(EVAL)

    sys.stderr.write(f"reading {args.old}\n")
    old = calls(ev, args.old)
    sys.stderr.write(f"  {len(old):,} variants\n")
    sys.stderr.write(f"reading {args.new}\n")
    new = calls(ev, args.new)
    sys.stderr.write(f"  {len(new):,} variants\n")

    truth = {}
    if args.truth:
        with open(args.truth) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                truth[(row["chrom"], row["pos"], row["ref"], row["alt"])] = row

    only_old = set(old) - set(new)
    only_new = set(new) - set(old)

    o, n = args.old_name, args.new_name
    changed_call = []
    changed_other = Counter()
    for key in sorted(set(old) & set(new)):
        a, b = old[key], new[key]
        if a["acmg"] != b["acmg"]:
            changed_call.append(key)
        else:
            for _name, k in FIELDS + [("", "criteria")]:
                if a[k] != b[k]:
                    changed_other[k] += 1

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["gene", "chrom", "pos", "ref", "alt", "clinvar_clnsig", "clinvar_stars",
            f"acmg_{o}", f"acmg_{n}", f"criteria_{o}", f"criteria_{n}",
            f"consequence_{o}", f"consequence_{n}", f"impact_{o}", f"impact_{n}",
            f"hgvsc_{o}", f"hgvsc_{n}", f"hgvsp_{o}", f"hgvsp_{n}",
            f"transcript_{o}", f"transcript_{n}"]
    with out.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for key in changed_call:
            a, b = old[key], new[key]
            t = truth.get(key, {})
            w.writerow([a["gene"] or b["gene"], key[0], key[1], key[2], key[3],
                        t.get("clnsig", ""), t.get("review_stars", ""),
                        a["acmg"], b["acmg"], a["criteria"], b["criteria"],
                        a["consequence"], b["consequence"], a["impact"], b["impact"],
                        a["hgvsc"], b["hgvsc"], a["hgvsp"], b["hgvsp"],
                        a["transcript"], b["transcript"]])

    print(f"{o} variants:            {len(old):,}")
    print(f"{n} variants:            {len(new):,}")
    print(f"only in {o}:             {len(only_old):,}")
    print(f"only in {n}:             {len(only_new):,}")
    print(f"ACMG call changed:       {len(changed_call):,}   -> {out}")
    print("same call, field changed:")
    for k, v in changed_other.most_common():
        print(f"  {k:<12} {v:,}")


if __name__ == "__main__":
    main()
