#!/usr/bin/env python3
"""Allele strings for the differential test of VARIANT_CLASS.

One Ensembl allele string per line ("REF/ALT1/ALT2", "-" for an empty allele).
What each one is classed as is decided by VEP's own `SO_variation_class` (see
regen.sh), not here. The strings are built to reach every branch of that
function: single bases, longer alleles, deletions and insertions, mixed
multi-allelic sites, the odd tokens of the older dbSNP imports, and strings that
should fall through to `sequence_alteration`.

Deterministic for a given seed.
"""
import argparse
import random

BASES = "ACGT"


def seq(rng, lo=1, hi=4):
    return "".join(rng.choice(BASES) for _ in range(rng.randint(lo, hi)))


ODD = [
    "(LARGEDELETION)", "(LARGEINSERTION)", "(CAG)8", "(CAG)9", "(2345 BP DELETION)",
    "(308 BP DELETION)", "(5 BP INSERTION)", "DEL", "INS", "N", "R", "Y", "a", "ACGTN", "",
]


def handwritten():
    return [
        "A/G", "A/G/T", "C/T", "AA/TT", "AC/GT", "AA/TTTT", "A/AT", "AT/A", "A/-", "ACGT/-",
        "-/A", "-/ACGT", "-/A/T", "A/-/T", "-/A/-", "A/G/ACG", "AA/A/AAA", "AA/-/TT", "N/A", "A/N",
        "R/Y", "(CAG)8/(CAG)9", "(CAG)8/(CAG)7", "-/(LARGEDELETION)", "-/(LARGEINSERTION)",
        "(LARGEDELETION)/-", "(2345 BP DELETION)/-", "(LARGEINSERTION)", "(308 BP DELETION)",
        "A", "A/", "/A", "A//G", "a/g", "A/g", "-", "-/-", "A/A", "A/G/G", "ACGT/ACGT",
        "A/(LARGEDELETION)", "A/(5 BP INSERTION)", "-/A/(LARGEINSERTION)", "AC/-/A",
    ]


def random_string(rng):
    kind = rng.random()
    if kind < 0.25:  # SNV-like, possibly multi-allelic
        n = rng.randint(2, 4)
        return "/".join(rng.choice("ACGTNRY") for _ in range(n))
    if kind < 0.45:  # deletion shapes
        ref = seq(rng, 1, 5)
        alts = ["-"] * rng.randint(1, 2) if rng.random() < 0.7 else ["-", seq(rng, 1, 3)]
        return "/".join([ref] + alts)
    if kind < 0.65:  # insertion shapes
        alts = [seq(rng, 1, 4) for _ in range(rng.randint(1, 3))]
        if rng.random() < 0.2:
            alts.append("-")
        return "/".join(["-"] + alts)
    if kind < 0.85:  # substitutions and indels
        ref = seq(rng, 1, 4)
        alts = [seq(rng, 1, 5) for _ in range(rng.randint(1, 3))]
        return "/".join([ref] + alts)
    parts = [rng.choice(ODD + [seq(rng), "-"]) for _ in range(rng.randint(1, 4))]
    return "/".join(parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=139)
    ap.add_argument("--count", type=int, default=60000)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    seen = set()
    for s in handwritten() + [random_string(rng) for _ in range(args.count)]:
        if s not in seen:
            seen.add(s)
            print(s)


if __name__ == "__main__":
    main()
