#!/usr/bin/env python3
"""Cut a reproducible sample out of a VCF, for the VEP divergence comparison.

Two shapes, because the two ends of the range in
[`docs/VEP_DIVERGENCE.md`](../docs/VEP_DIVERGENCE.md) need different samples:

    --stratified   the hard shapes. Quotas per variant class - 3,000 SNV /
                   1,800 deletion / 1,200 insertion / 400 MNV / 200 delins -
                   so 54.5 % of the 6,600 are not SNVs, against 7.2 % of the
                   ClinVar 2-star+ set they are drawn from.

    --every N      an ordinary callset, unweighted: every Nth record, which is
                   what the genome-wide HG002 row is measured on.

Both walk the file in order and take every k-th record of a class, so the
sample is a function of the input alone. Nothing here is random, and nothing
here is a seed that has to be carried around with the numbers: re-running this
on the same VCF reproduces the sample that produced them.
"""

import argparse
import sys

QUOTAS = {"snv": 3000, "del": 1800, "ins": 1200, "mnv": 400, "delins": 200}


def variant_class(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(ref) == 1:
        return "ins"
    if len(alt) == 1:
        return "del"
    if len(ref) == len(alt):
        return "mnv"
    return "delins"


def contig_order(chrom: str) -> int:
    named = {"X": 23, "Y": 24, "MT": 25, "M": 25}
    bare = chrom[3:] if chrom.startswith("chr") else chrom
    if bare.isdigit():
        return int(bare)
    return named.get(bare, 99)


def sort_key(line: str):
    f = line.split("\t", 5)
    return (contig_order(f[0]), int(f[1]), f[3], f[4])


def stratified(src, quotas):
    header, pools = [], {k: [] for k in quotas}
    for line in src:
        if line.startswith("#"):
            header.append(line)
            continue
        f = line.split("\t", 5)
        # A multi-allelic record is one variant to the quota and several to the
        # comparison, which would weight it twice. The genome-wide sample is
        # where those are measured.
        if "," in f[4]:
            continue
        pools[variant_class(f[3], f[4])].append(line)

    picked = []
    for name, want in quotas.items():
        pool = pools[name]
        step = max(1, len(pool) // want)
        taken = pool[::step][:want]
        print(
            f"{name:8s} pool={len(pool):8d} step={step:5d} taken={len(taken)}",
            file=sys.stderr,
        )
        picked.extend(taken)
    picked.sort(key=sort_key)
    return header, picked


def systematic(src, n):
    header, picked, seen = [], [], 0
    for line in src:
        if line.startswith("#"):
            header.append(line)
            continue
        seen += 1
        if seen % n == 0:
            picked.append(line)
    print(f"1-in-{n} of {seen} records: {len(picked)}", file=sys.stderr)
    return header, picked


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input", help="VCF to sample")
    ap.add_argument("output", help="where to write the sample")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--stratified", action="store_true", help="quota sample of the hard shapes"
    )
    mode.add_argument("--every", type=int, metavar="N", help="every Nth record")
    args = ap.parse_args()

    with open(args.input) as src:
        if args.stratified:
            header, picked = stratified(src, QUOTAS)
        else:
            header, picked = systematic(src, args.every)

    with open(args.output, "w") as out:
        out.writelines(header)
        out.writelines(picked)
    print(f"wrote {len(picked)} records to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
