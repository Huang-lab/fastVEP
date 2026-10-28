#!/usr/bin/env python3
"""Does an allele get the same cells in a record with several ALTs as alone in a record of its own?

usage: compare_alone.py MULTI.vcf MULTI.tab ALONE.tab [--show N]

MULTI.tab is a tool's output for MULTI.vcf, and ALONE.tab its output for the same file passed
through `split_alts.py`.
Each row of a multi-allelic record is paired with the row of the same allele and the same
transcript in the other file, and every column but the ones that name the allele or its span
(Allele, Location, REF_ALLELE, and the key) is compared, VARIANT_CLASS too being left out:
both tools write the class of the whole record there, so an allele's differs when its record's does.
The allele is found by ALLELE_NUM when the file has one (Ensembl VEP with --allele_number),
and by the Allele column otherwise, which for a record with several ALTs is the ALT with the
base they all share taken off, as both tools write it.

A difference is not a defect on either side.
It says the tool describes an allele differently according to what else the record holds, and
how much: that is what `docs/VEP_DIVERGENCE.md` means by Ensembl 116 no longer splitting them.
"""
import argparse
import collections


def read(path):
    header, rows = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line[1:].split("\t")
                continue
            rows.append(dict(zip(header, line.split("\t"))))
    return header, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("multi_vcf")
    ap.add_argument("multi_tab")
    ap.add_argument("alone_tab")
    ap.add_argument("--show", type=int, default=3)
    a = ap.parse_args()

    records = {}
    with open(a.multi_vcf) as f:
        for n, line in enumerate(f):
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            alts = c[4].split(",")
            if len(alts) < 2:
                continue
            name = c[2].split(";")[0]
            shown = name if name != "." else f"{c[0]}_{c[1]}_{c[3]}/{'/'.join(alts)}"
            base = name if name != "." else f"{c[0]}_{c[1]}_{c[3]}"
            records[shown] = (n, base, c[3], alts)

    mh, multi = read(a.multi_tab)
    ah, alone = read(a.alone_tab)
    skip = {"Uploaded_variation", "Allele", "Location", "REF_ALLELE", "ALLELE_NUM", "SOURCE", "VARIANT_CLASS"}
    columns = [c for c in mh if c in ah and c not in skip and not c.endswith(".gz")]
    by_alone = collections.defaultdict(list)
    for r in alone:
        by_alone[(r["Uploaded_variation"], r["Feature"])].append(r)

    paired = 0
    cells = collections.Counter()
    shapes = collections.Counter()
    shown = {}
    for r in multi:
        rec = records.get(r["Uploaded_variation"])
        if rec is None:
            continue
        n, base, ref, alts = rec
        if "ALLELE_NUM" in r:
            i = int(r["ALLELE_NUM"])
        else:
            chop = len({x[0] for x in [ref] + alts}) == 1
            i = next(
                (j for j, x in enumerate(alts, 1) if ((x[1:] if chop else x) or "-") == r["Allele"]),
                None,
            )
            if i is None:
                continue
        other = by_alone.get((f"{base}_n{n}_alt{i}", r["Feature"]), [])
        if len(other) != 1:
            continue
        paired += 1
        d = tuple(c for c in columns if r[c] != other[0][c])
        for c in d:
            cells[c] += 1
        if d:
            shapes[d] += 1
            shown.setdefault(d, (r["Uploaded_variation"], i, r["Feature"], {c: (r[c], other[0][c]) for c in d}))
    print(f"{paired} rows of records with several ALTs, paired with the same allele alone;"
          f" {sum(shapes.values())} differ in at least one column")
    print("  cells that differ, by column:", dict(cells.most_common()))
    for d, count in shapes.most_common(a.show):
        name, i, feature, diff = shown[d]
        print(f"  {count:5d}  {', '.join(d)}   e.g. {name} ALT {i} {feature}")
        for c, (m, o) in diff.items():
            print(f"            {c}: with its record {m!r} | alone {o!r}")


main()
