#!/usr/bin/env python3
"""Candidate (input, known variant) pairs for the differential test of the
known-variant allele matcher.

Each output line is one pair, tab separated, in the shape VEP's
`get_matched_variant_alleles` takes it:

    a_allele_string  a_pos  a_strand  b_allele_string  b_pos  b_strand

`a` is the input site, `b` the known variant from the cache. Positions are the
Ensembl start (an insertion sits at the base it is inserted before).

The pairs are built to hit what the matcher has to get right, not to be
realistic: the same event written minimally, anchored on the left, anchored on
the right and anchored on both sides; events that differ by a shift inside a
repeat; multi-allelic sites on either side; the other strand; and plain random
alleles. Whether any pair matches is decided by VEP's own code (see regen.sh),
not here.

Deterministic for a given seed: nothing here reads the clock or the machine.
"""
import argparse
import itertools
import random

# A reference with a period-6 repeat, homopolymers and a dinucleotide repeat, so
# that equivalent events written at different positions exist.
REF = "GATTACAACCCTAACCCTAACCCTAAGGGTTTTCCCAAAGTCAGCACACACATGCAAAAAGGCCGT"
ORIGIN = 1000

INSERTS = ["A", "C", "T", "CC", "TA", "CCT", "ACA"]


def event_forms(p, d, ins):
    """All spellings of one event: delete REF[p:p+d], insert `ins` there.

    Returns {name: (ref, alt, pos)}; an empty allele is written "-"."""
    deleted = REF[p : p + d]
    forms = {"minimal": (deleted or "-", ins or "-", ORIGIN + p)}
    if p >= 1:
        a = REF[p - 1]
        forms["left"] = (a + deleted, a + ins, ORIGIN + p - 1)
    if p + d < len(REF):
        b = REF[p + d]
        forms["right"] = (deleted + b, ins + b, ORIGIN + p)
    if p >= 1 and p + d < len(REF):
        a, b = REF[p - 1], REF[p + d]
        forms["both"] = (a + deleted + b, a + ins + b, ORIGIN + p - 1)
    return forms


def events(rng):
    """(p, d, ins) events: SNVs, deletions, insertions and delins."""
    out = []
    for p in range(1, len(REF) - 6):
        base = REF[p]
        for alt in "ACGT":
            if alt != base:
                out.append((p, 1, alt))
        for d in (1, 2, 3, 4):
            out.append((p, d, ""))
        for ins in INSERTS:
            out.append((p, 0, ins))
        out.append((p, 2, rng.choice(INSERTS)))  # a delins
    return out


def allele_string(ref, alts):
    return "/".join([ref] + alts)


def pair_line(a, b, b_strand=1):
    (a_ref, a_alts, a_pos), (b_ref, b_alts, b_pos) = a, b
    return "\t".join(
        [allele_string(a_ref, a_alts), str(a_pos), "1", allele_string(b_ref, b_alts), str(b_pos), str(b_strand)]
    )


def single(form):
    ref, alt, pos = form
    return (ref, [alt], pos)


def structured(rng, n):
    """Pairs of events near each other, in every combination of spellings."""
    evs = events(rng)
    by_p = {}
    for e in evs:
        by_p.setdefault(e[0], []).append(e)
    names = ["minimal", "left", "right", "both"]
    lines = []
    while len(lines) < n:
        x = rng.choice(evs)
        near = [e for q in range(x[0] - 3, x[0] + 4) for e in by_p.get(q, [])]
        y = x if rng.random() < 0.3 else rng.choice(near)
        fx, fy = event_forms(*x), event_forms(*y)
        nx, ny = rng.choice([k for k in names if k in fx]), rng.choice([k for k in names if k in fy])
        lines.append(pair_line(single(fx[nx]), single(fy[ny]), rng.choice([1, 1, 1, 1, -1])))
    return lines


def multiallelic_group(rng):
    """Several events with one (p, d) share a reference allele: a multi-allelic site."""
    p = rng.randrange(1, len(REF) - 6)
    d = rng.choice([0, 1, 1, 2])
    pool = [i for i in INSERTS if i != REF[p : p + d]] + ["G", "T", "C", "A"]
    inserts = rng.sample(sorted(set(pool)), rng.randint(2, 3))
    inserts = [i for i in inserts if i != REF[p : p + d]] or ["A"]
    return p, d, inserts


def site(rng, p, d, inserts, spelling):
    """The multi-allelic site (p, d, inserts) in one spelling, which is
    (ref, alts, pos). Every allele gets the same flanking bases, so they share
    one reference allele."""
    forms = [event_forms(p, d, ins) for ins in inserts]
    if spelling not in forms[0]:
        spelling = "minimal"
    ref, _, pos = forms[0][spelling]
    return (ref, [f[spelling][1] for f in forms], pos)


def multiallelic(rng, n):
    lines = []
    while len(lines) < n:
        p, d, inserts = multiallelic_group(rng)
        names = ["minimal", "left", "right", "both"]
        a = site(rng, p, d, inserts, rng.choice(names))
        # b: the same group, a subset of it, or a neighbouring group
        kind = rng.choice(["same", "subset", "single", "near"])
        if kind == "near":
            q, e, more = multiallelic_group(rng)
            q = max(1, min(len(REF) - 7, p + rng.randint(-2, 2)))
            b = site(rng, q, e, more, rng.choice(names))
        elif kind == "single":
            forms = event_forms(p, d, rng.choice(inserts))
            b = single(forms[rng.choice([k for k in names if k in forms])])
        elif kind == "subset":
            b = site(rng, p, d, inserts[: max(1, len(inserts) - 1)], rng.choice(names))
        else:
            b = site(rng, p, d, inserts, rng.choice(names))
        lines.append(pair_line(a, b, rng.choice([1, 1, 1, -1])))
    return lines


def random_allele(rng, lo=0, hi=4):
    return "".join(rng.choice("ACGT") for _ in range(rng.randint(lo, hi)))


def random_pairs(rng, n):
    lines = []
    while len(lines) < n:
        def site_(pos):
            ref = random_allele(rng) or "-"
            alts = []
            for _ in range(rng.choice([1, 1, 1, 2, 3])):
                alt = random_allele(rng) or "-"
                if alt != ref:
                    alts.append(alt)
            return (ref, alts, pos) if alts else None

        a_pos = rng.randint(100, 110)
        a = site_(a_pos)
        b = site_(a_pos + rng.randint(-3, 3))
        if a and b:
            lines.append(pair_line(a, b, rng.choice([1, 1, 1, 1, -1])))
    return lines


def handwritten():
    """The shapes issue #139 reports, and the corners of the Perl."""
    H = [
        # the reporter's anchored deletion and insertion, and their minimal forms
        "ACCCT/-\t10401\t1\tACCCTA/A\t10401\t1",
        "ACCCT/-\t10401\t1\tAACCC/-\t10400\t1",
        "AACCC/-\t10400\t1\tAACCC/-\t10400\t1",
        "-/T\t10409\t1\t-/T\t10409\t1",
        "-/A\t10421\t1\tG/GA\t10420\t1",
        "-/A\t10421\t1\tGA/G\t10420\t1",
        # SNVs, multi-allelic on either side
        "C/A\t10414\t1\tC/A\t10414\t1",
        "C/A\t10414\t1\tC/G\t10414\t1",
        "C/A\t10414\t1\tC/A/G\t10414\t1",
        "C/A/G\t10414\t1\tC/G\t10414\t1",
        "C/A/G\t10414\t1\tC/A/G\t10414\t1",
        "C/G/A\t10414\t1\tC/A/G\t10414\t1",
        # a mixed multi-allelic site the VCF parser could not strip an anchor from
        "A/G/ACG\t100\t1\t-/CG\t101\t1",
        "A/G/ACG\t100\t1\tA/ACG\t100\t1",
        "A/G/ACG\t100\t1\tA/G\t100\t1",
        # an allele string with no alleles in it
        "C/A\t10414\t1\tCOSMIC_MUTATION\t10414\t1",
        # other strand
        "C/T\t100\t1\tG/A\t100\t-1",
        "C/T\t100\t1\tG/A\t100\t1",
        "ACG/-\t100\t1\tCGT/-\t100\t-1",
        # same sequence, different position
        "C/A\t100\t1\tC/A\t101\t1",
        "-/T\t100\t1\t-/T\t101\t1",
        # an MNV that trims to an SNV only from the left
        "AAC/AAT\t100\t1\tC/T\t102\t1",
        "C/T\t102\t1\tAAC/AAT\t100\t1",
        "AAC/AAT\t100\t1\tCAC/CAT\t100\t1",
        # right trim moves no position; left trim does
        "TCG/TG\t100\t1\tCG/G\t101\t1",
        "TCG/TG\t100\t1\tTC/T\t100\t1",
        "CCC/CC\t100\t1\tCC/C\t100\t1",
        "CCC/CC\t100\t1\tCC/C\t101\t1",
        "CCC/CC\t100\t1\tC/-\t100\t1",
        "CCC/CC\t100\t1\tC/-\t102\t1",
        # delins on both sides
        "ACGT/TT\t100\t1\tACGT/TT\t100\t1",
        "ACGT/TT\t100\t1\tCG/T\t101\t1",
        "ACGT/TT\t100\t1\tACG/T\t100\t1",
        # single base alleles, where only one trim order is tried
        "A/-\t100\t1\tA/-\t100\t1",
        "-/A\t100\t1\t-/A\t100\t1",
        "A/-\t100\t1\tGA/G\t99\t1",
        "A/-\t100\t1\tAG/G\t100\t1",
    ]
    return H


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=139)
    ap.add_argument("--structured", type=int, default=40000)
    ap.add_argument("--multiallelic", type=int, default=15000)
    ap.add_argument("--random", type=int, default=15000)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    generated = structured(rng, args.structured)
    generated += multiallelic(rng, args.multiallelic)
    generated += random_pairs(rng, args.random)
    # Shuffled, so that taking the first N matches (regen.sh does) samples every
    # family instead of only the first one generated.
    rng.shuffle(generated)
    seen = set()
    for line in handwritten() + generated:
        if line not in seen:
            seen.add(line)
            print(line)


if __name__ == "__main__":
    main()
