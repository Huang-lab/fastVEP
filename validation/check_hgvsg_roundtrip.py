#!/usr/bin/env python3
"""Does each tool's HGVSg, applied to the genome, give the allele the VCF gave?

usage: check_hgvsg_roundtrip.py WORLD_DIR DIFF_HGVSg.tsv

WORLD_DIR holds genome.fa and in.vcf, DIFF_HGVSg.tsv is the file `compare_tab.py --out` writes for the
HGVSg column (Uploaded_variation, Allele, Feature, VEP, fastVEP).
For every (record, allele) in it the allele is applied to the genome as the VCF says, and each tool's
description as HGVS says.
Two spellings of one change give the same edited sequence however they are written (shifted or not,
dup or ins), so equal sequences mean the descriptions differ only in how they spell the change, and
unequal ones mean a description is not of this variant.
"""
import collections
import re
import sys


def main():
    world, diff = sys.argv[1:3]
    genome, name = {}, None
    for line in open(f"{world}/genome.fa"):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            genome[name] = []
        else:
            genome[name].append(line.upper())
    genome = {k: "".join(v) for k, v in genome.items()}

    records = {}
    for line in open(f"{world}/in.vcf"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        alts = f[4].split(",")
        shown = f[2].split(";")[0] if f[2] != "." else f"{f[0]}_{f[1]}_{f[3]}/{'/'.join(alts)}"
        records[shown] = (f[0], int(f[1]), f[3], alts)

    def apply_vcf(chrom, pos, ref, alt):
        g = genome[chrom]
        return g[: pos - 1] + alt + g[pos - 1 + len(ref):]

    def apply_description(desc):
        m = re.match(r"^(?P<c>[^:]+):[gm]\.(?P<body>.+)$", desc)
        c, body = m["c"], m["body"]
        g = genome[c]
        sub = re.match(r"^(\d+)([ACGT])>([ACGT])$", body)
        if sub:
            p = int(sub[1])
            assert g[p - 1] == sub[2], f"{desc}: reference base is {g[p - 1]}"
            return g[: p - 1] + sub[3] + g[p:]
        r = re.match(r"^(\d+)(?:_(\d+))?(del|dup|inv|ins([ACGT]+)|delins([ACGT]+))$", body)
        if not r:
            raise ValueError(f"cannot read {desc}")
        s = int(r[1])
        e = int(r[2]) if r[2] else s
        kind = r[3]
        if kind == "del":
            return g[: s - 1] + g[e:]
        if kind == "dup":
            return g[:e] + g[s - 1:e] + g[e:]
        if kind == "inv":
            comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
            return g[: s - 1] + "".join(comp[b] for b in reversed(g[s - 1:e])) + g[e:]
        if kind.startswith("delins"):
            return g[: s - 1] + r[5] + g[e:]
        assert e == s + 1, desc
        return g[:s] + r[4] + g[s:]

    def pick_alt(ref, alts, allele):
        if len(alts) == 1:
            return alts[0]
        chop = len({a[0] for a in [ref] + alts}) == 1
        return next((a for a in alts if ((a[1:] if chop else a) or "-") == allele), None)

    tally, shown, seen = collections.Counter(), collections.defaultdict(list), set()
    for line in list(open(diff))[1:]:
        up, allele, _feature, vep, ours = line.rstrip("\n").split("\t")
        if (up, allele) in seen or up not in records:
            continue
        seen.add((up, allele))
        chrom, pos, ref, alts = records[up]
        alt = pick_alt(ref, alts, allele)
        if alt is None:
            tally["allele not found"] += 1
            continue
        want = apply_vcf(chrom, pos, ref, alt)
        verdict = []
        for who, d in (("VEP", vep), ("fastVEP", ours)):
            if d == "-":
                verdict.append(f"{who} writes none")
            else:
                verdict.append(f"{who} {'reproduces the allele' if apply_description(d) == want else 'does NOT reproduce the allele'}")
        key = ", ".join(verdict)
        tally[key] += 1
        if len(shown[key]) < 3:
            shown[key].append((up, allele, vep, ours))
    print(f"{sum(tally.values())} (record, allele) pairs whose HGVSg differs:")
    for key, n in tally.most_common():
        print(f"  {n:5d}  {key}")
        for up, allele, vep, ours in shown[key]:
            print(f"           {up} allele {allele}: VEP {vep} | fastVEP {ours}")


main()
