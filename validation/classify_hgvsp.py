#!/usr/bin/env python3
"""Group the cells of an HGVSp diff by the shape of what each tool wrote.

usage: classify_hgvsp.py diff_HGVSp.tsv [N]

`diff_HGVSp.tsv` is what `compare_tab.py --out DIR` writes for the column. A cell is reduced to its
shape, with every number `N` and every run of three-letter residues a single `R` (two or more
`R2`, `R3`...), so `p.Ter583delextTer23` is `RNdelextRN` and `p.Leu89TyrfsTer42` is `RNRfsRN`, and
the cells are counted by the pair of shapes with the first one seen as the example. A few
hundred thousand differing cells are a few dozen classes, and a class is a cause or a few.
"""
import collections
import re
import sys

RESIDUE = r"(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Ter|Xaa|Sec|Pyl)"


def shape(cell):
    if cell == "-":
        return "-"
    cell = cell.split(":p.", 1)[-1].replace("%3D", "=")
    cell = re.sub(r"\d+", "N", cell)
    cell = re.sub(RESIDUE + "+", lambda m: "R" if len(m.group(0)) == 3 else f"R{len(m.group(0)) // 3}", cell)
    return cell


def main():
    path = sys.argv[1]
    top = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    groups = collections.defaultdict(list)
    with open(path) as f:
        next(f)
        for line in f:
            variant, allele, feature, vep, ours = line.rstrip("\n").split("\t")
            groups[(shape(vep), shape(ours))].append((variant, feature, vep, ours))
    print(sum(len(v) for v in groups.values()), "cells in", len(groups), "classes")
    for (sv, so), items in sorted(groups.items(), key=lambda kv: -len(kv[1]))[:top]:
        variant, feature, vep, ours = items[0]
        print(f"{len(items):7d}  VEP {sv:24s} fastVEP {so:24s} e.g. {variant} {feature}")
        print(f"         VEP    : {vep}")
        print(f"         fastVEP: {ours}")


main()
