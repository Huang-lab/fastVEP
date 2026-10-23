#!/usr/bin/env python3
"""Write a VCF in which every ALT of a multi-allelic record is a record of its own.

usage: split_alts.py IN.vcf OUT.vcf

A record with one ALT is copied as it is.
A record with several becomes one record per ALT, named `<name>_n<line>_alt<i>`: <name> is its
first ID, or CHROM_POS_REF when it has none, <line> is its line in the input (so that two
records at one site stay apart) and <i> counts the ALTs from 1, which is the ALLELE_NUM Ensembl
VEP gives the allele in the original record.
`compare_alone.py` reads the names back.
"""
import sys


def main():
    src, dst = sys.argv[1:3]
    names = set()
    with open(src) as fin, open(dst, "w") as out:
        for n, line in enumerate(fin):
            if line.startswith("#"):
                out.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            alts = f[4].split(",")
            if len(alts) == 1:
                out.write(line)
                continue
            name = f[2].split(";")[0]
            if name == ".":
                name = f"{f[0]}_{f[1]}_{f[3]}"
            for i, alt in enumerate(alts, 1):
                g = list(f)
                g[2] = f"{name}_n{n}_alt{i}"
                assert g[2] not in names, g[2]
                names.add(g[2])
                g[4] = alt
                out.write("\t".join(g) + "\n")


main()
