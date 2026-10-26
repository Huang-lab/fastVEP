#!/usr/bin/env python3
"""Lone indels written with bases the REF and the ALT share beyond the one-base anchor.

    gen_lone.py OUT.vcf [N [SEED]]

The gene models and the genome are `../vep_tab_world`'s. Every record has one ALT that is an
indel, and the pair is written as a caller that does not trim might write it: a prefix both
alleles share (one to three bases; the first is the anchor VCF asks for, the others are what
a caller that does not left-align or trim leaves in), the bases that change, and a suffix they
share (none to three bases):

    ref = prefix + removed + suffix        alt = prefix + added + suffix

with the removed and added bases differing at both of their ends, so that what is left once
the shared bases are taken off is exactly `removed` and `added`. The three shapes are a
deletion (nothing added), an insertion (nothing removed) and a replacement of one length by
another.

Ensembl VEP minimises such a record from both ends before it annotates it
(`Parser::post_process_vfs`, `minimise_alleles`, `trim_sequences`), which changes the
`Allele`, the `Location`, the `REF_ALLELE` and the `VARIANT_CLASS` of every row, and what the
consequence step sees. A record with several ALTs is not minimised in 116 and is not here.

The positions are spread over the genes, their edges, and the intergenic sequence between.
"""
import os
import random
import sys

OUT = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 1500
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 139
WORLD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "vep_tab_world")
rng = random.Random(SEED)

genome = "".join(l.strip() for l in open(f"{WORLD}/genome.fa") if not l.startswith(">"))
LENGTH = len(genome)

exons, cds = [], []
for line in open(f"{WORLD}/genes.gff3"):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9:
        continue
    if f[2] == "exon":
        exons.append((int(f[3]), int(f[4])))
    elif f[2] == "CDS":
        cds.append((int(f[3]), int(f[4])))
edges = sorted({e for s, t in exons + cds for e in (s, t + 1)})


def rand_seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def position():
    while True:
        r = rng.random()
        if r < 0.45:
            p = rng.choice(edges) + rng.randint(-14, 14)
        elif r < 0.85:
            s, e = rng.choice(exons)
            p = rng.randint(s, e)
        else:
            p = rng.randint(5, LENGTH - 60)
        if 5 <= p <= LENGTH - 60:
            return p


records = set()
guard = 0
while len(records) < N and guard < N * 50:
    guard += 1
    p = position()
    pre = rng.choice([1, 1, 2, 2, 3])          # the anchor, and what a caller leaves beyond it
    suf = rng.choice([0, 1, 1, 2, 3])
    shape = rng.choice(["del", "ins", "indel"])
    removed_n = {"del": rng.randint(1, 8), "ins": 0, "indel": rng.randint(1, 6)}[shape]
    added_n = {"del": 0, "ins": rng.randint(1, 8), "indel": rng.randint(1, 6)}[shape]
    if shape == "indel" and removed_n == added_n:
        added_n += 1
    prefix = genome[p - 1:p - 1 + pre]
    removed = genome[p - 1 + pre:p - 1 + pre + removed_n]
    after = p - 1 + pre + removed_n
    suffix = genome[after:after + suf]
    added = rand_seq(added_n)
    if len(prefix) != pre or len(removed) != removed_n or len(suffix) != suf:
        continue
    # nothing shared at the ends of what changes, or the pair would trim further than meant
    if removed and added and (removed[0] == added[0] or removed[-1] == added[-1]):
        continue
    ref = prefix + removed + suffix
    alt = prefix + added + suffix
    if ref == alt or "N" in ref + alt:
        continue
    records.add((p, ref, alt))

with open(OUT, "w") as v:
    v.write(f"##fileformat=VCFv4.2\n##contig=<ID=chr10,length={LENGTH}>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for p, ref, alt in sorted(records):
        v.write(f"chr10\t{p}\t.\t{ref}\t{alt}\t.\t.\t.\n")
print(f"{len(records)} records")
