#!/usr/bin/env python3
"""Every shape of change at every position around every boundary of the world's transcripts.

    gen_edges.py OUT.vcf [SAMPLE] [RADIUS]

The gene models and the genome are `../vep_tab_world`'s (`gen_world.py`, which leaves out
what happens at an exon, a coding or a transcript edge because its tab output is where the
two tools were known to differ; this is that part). The boundaries are the edges of every exon
and of every CDS segment, so they include the transcript ends and the start and stop codons,
on both strands and in coding and non-coding genes.

At each position within RADIUS (default 10) of a boundary one record is written per shape,
every shape its own single-ALT record so that no multi-allelic rule stands between the tools:

  - an SNV, and MNVs of 2 and 3 bases;
  - deletions of 1 to 6, 9, 12, 15, 20 and 30 bases;
  - insertions of 1 to 6 and 9 bases, random, a copy of the bases before and a copy of the
    bases after (the last two are duplications);
  - the delins that change length and share no base at either end (2 -> 3, 3 -> 2, 3 -> 4,
    4 -> 2);
  - whole-exon deletions, from a few bases before an exon to a few bases after it, and the
    same span replaced by four bases;
  - deletions that run from inside one exon, across the whole intron, into the next, which
    is a change over two exons (`EXON` `1-2/3`), and the same span replaced by four bases.

Every record is minimal and anchored by one base, so what differs is what the transcript
does with the change and not how the change was written.

That is about 17,700 records. SAMPLE (default 8) keeps every SAMPLE-th of them, in position order,
and all the whole-exon and over-an-intron records, which is the fixture: the full set is what the comparison in
docs/VEP_DIVERGENCE.md ("Part 3") was run on (`SAMPLE` 1).
"""
import os
import random
import sys

args = [a for a in sys.argv[1:] if not a.startswith("--")]
OUT = args[0]
SAMPLE = int(args[1]) if len(args) > 1 else 8
RADIUS = int(args[2]) if len(args) > 2 else 10
WORLD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "vep_tab_world")
rng = random.Random(7)

genome = "".join(l.strip() for l in open(f"{WORLD}/genome.fa") if not l.startswith(">"))
LENGTH = len(genome)

exons, cds = [], []
by_transcript = {}
for line in open(f"{WORLD}/genes.gff3"):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9:
        continue
    if f[2] == "exon":
        exons.append((int(f[3]), int(f[4])))
        parent = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)["Parent"]
        by_transcript.setdefault(parent, []).append((int(f[3]), int(f[4])))
    elif f[2] == "CDS":
        cds.append((int(f[3]), int(f[4])))

edges = set()
for s, e in exons + cds:
    edges.update([s, e + 1])
edges = sorted(edges)


def base(p):
    return genome[p - 1]


def other(b, k=0):
    return [x for x in "ACGT" if x != b][k % 3]


def rand_seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def shapes(p):
    """(ref, alt) pairs, in VCF form, for a change anchored at p."""
    out = []
    r = base(p)
    out.append((r, other(r, p)))
    for n in (2, 3):
        ref = genome[p - 1:p - 1 + n]
        out.append((ref, "".join(other(c, p + i) for i, c in enumerate(ref))))
    for n in (1, 2, 3, 4, 5, 6, 9, 12, 15, 20, 30):
        ref = genome[p - 1:p + n]
        out.append((ref, ref[0]))
    for n in (1, 2, 3, 4, 5, 6, 9):
        out.append((r, r + rand_seq(n)))
        if p - 1 >= n and p + n <= LENGTH:
            out.append((r, r + genome[p:p + n]))  # a copy of what follows
            out.append((r, r + genome[p - n:p]))  # a copy of what precedes
    for ra, la in ((2, 3), (3, 2), (3, 4), (4, 2)):
        ref = genome[p - 1:p - 1 + ra]
        while True:
            alt = rand_seq(la)
            if alt[0] != ref[0] and alt[-1] != ref[-1]:
                break
        out.append((ref, alt))
    return out


records = set()
for e in edges:
    for p in range(e - RADIUS, e + RADIUS + 1):
        if not (5 <= p <= LENGTH - 40):
            continue
        for ref, alt in shapes(p):
            if p - 1 + len(ref) <= LENGTH and ref != alt and set(ref + alt) <= set("ACGT"):
                records.add((p, ref, alt))

whole_exon = set()
for s, e in sorted(set(exons)):
    for a in (3, 8):
        for b in (3, 8):
            p = s - a - 1
            n = (e + b) - p
            if p >= 5 and p + n <= LENGTH - 2:
                ref = genome[p - 1:p + n]
                whole_exon.add((p, ref, ref[0]))
                alt = rand_seq(4)
                if alt[0] != ref[0] and alt[-1] != ref[-1]:
                    whole_exon.add((p, ref, alt))

# from inside one exon, over the whole intron, into the next
over_intron = set()
for exs in by_transcript.values():
    exs = sorted(exs)
    for (s1, e1), (s2, e2) in zip(exs, exs[1:]):
        for a, b in ((2, 2), (8, 8)):
            p = e1 - a
            n = (s2 + b) - p
            if p >= 5 and p + n <= LENGTH - 2:
                ref = genome[p - 1:p + n]
                over_intron.add((p, ref, ref[0]))
                alt = rand_seq(4)
                if alt[0] != ref[0] and alt[-1] != ref[-1]:
                    over_intron.add((p, ref, alt))

ordered = sorted(records)
kept = set(ordered[::SAMPLE]) | whole_exon | over_intron
with open(OUT, "w") as v:
    v.write(f"##fileformat=VCFv4.2\n##contig=<ID=chr10,length={LENGTH}>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for p, ref, alt in sorted(kept):
        v.write(f"chr10\t{p}\t.\t{ref}\t{alt}\t.\t.\t.\n")
print(f"{len(kept)} of {len(records) + len(whole_exon) + len(over_intron)} records at {len(edges)} edges")
