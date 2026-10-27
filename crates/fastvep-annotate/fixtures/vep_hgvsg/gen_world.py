#!/usr/bin/env python3
"""A genome of repeats and edges, and variants placed to make HGVS g. hard, for holding
`hgvsg_vep` to what Ensembl's VEP writes.

    gen_world.py OUTDIR [SEED]

`chr10` has a homopolymer, a dinucleotide, a trinucleotide and a hexamer repeat, and
homopolymers at its very start and its very end, so that a shift runs into the edge of the
contig. `MT` is the mitochondrion, whose numbering letter is `m`. Each contig has a small
gene, because VEP's --gff mode skips a contig no feature sits on.

Variants: every kind of change at every position of every repeat (a unit inserted or
deleted, a base inserted, an SNV, a delins), at the first, the middle and the last base of
the run and one either side, near both ends of the contig, and records of several alts.
A record with an allele that trims to a substitution is left out, as in the tab world
(see gen_world.py there): VEP writes a delins of the untrimmed alt over the trimmed span.

Writes OUTDIR/{genome.fa, genes.gff3, in.vcf}. Every record has an ID, `v<N>`, and the ones at
the very edges of a contig come last and are `e<N>`.
"""
import random
import sys

OUT = sys.argv[1]
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 139
rng = random.Random(SEED)

LEN = {"chr10": 1500, "MT": 1000}
seq = {c: [rng.choice("ACGT") for _ in range(n)] for c, n in LEN.items()}


def put(contig, pos, text):
    for i, ch in enumerate(text):
        seq[contig][pos - 1 + i] = ch


repeats = {  # contig: [(start, unit, copies)]
    "chr10": [(200, "A", 11), (300, "CA", 10), (400, "GAT", 6), (500, "CTGAAC", 4), (600, "T", 7),
              (1, "G", 12), (LEN["chr10"] - 11, "T", 12)],
    "MT": [(300, "C", 9), (500, "AT", 8), (700, "GGT", 5)],
}
for contig, runs in repeats.items():
    for start, unit, copies in runs:
        put(contig, start, unit * copies)
genome = {c: "".join(b) for c, b in seq.items()}

with open(f"{OUT}/genome.fa", "w") as f:
    for c, g in genome.items():
        f.write(f">{c}\n")
        for i in range(0, len(g), 60):
            f.write(g[i:i + 60] + "\n")
with open(f"{OUT}/genes.gff3", "w") as f:
    f.write("##gff-version 3\n")
    for c, (s, e) in (("chr10", (1000, 1100)), ("MT", (950, 990))):
        f.write(f"{c}\ttest\tgene\t{s}\t{e}\t.\t+\t.\tID=gene:G{c};Name=G{c};biotype=lncRNA\n")
        f.write(f"{c}\ttest\tlnc_RNA\t{s}\t{e}\t.\t+\t.\tID=transcript:T{c};Parent=gene:G{c};biotype=lncRNA\n")
        f.write(f"{c}\ttest\texon\t{s}\t{e}\t.\t+\t.\tID=exon:E{c};Parent=transcript:T{c};rank=1\n")

comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
other = lambda b: rng.choice([x for x in "ACGT" if x != b])
rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))


def trims_to_a_substitution(ref, alt):
    r, a = ref, alt
    while r and a and r[0] == a[0]:
        r, a = r[1:], a[1:]
    while r and a and r[-1] == a[-1]:
        r, a = r[:-1], a[:-1]
    return (r, a) != (ref, alt) and bool(r) and bool(a)


def events(contig, pos):
    """Alts at pos, as (ref, alt) in VCF form."""
    g = genome[contig]
    b = g[pos - 1]
    out = [(b, other(b))]  # SNV
    for n in (1, 2, 3, 6):
        if pos + n <= len(g):
            out.append((g[pos - 1:pos + n], b))  # a deletion of n bases after the anchor
        out.append((b, b + g[pos:pos + n]))  # an insertion that copies what follows
        if pos - n >= 1:
            out.append((b, b + g[pos - n:pos]))  # an insertion that copies what precedes
    out.append((b, b + rand(rng.randint(1, 3))))
    if pos + 2 <= len(g):
        ref = g[pos - 1:pos + 2]
        out.append((ref, other(ref[0]) + rand(rng.randint(1, 4))))  # a delins
        out.append((ref, "".join(other(c) for c in ref)))  # an MNV
    return out


records = []
for contig, runs in repeats.items():
    spots = set()
    for start, unit, copies in runs:
        last = start + len(unit) * copies - 1
        spots.update({start - 1, start, start + len(unit), (start + last) // 2, last - 1, last, last + 1})
    spots.update({1, 2, 3, 5, LEN[contig] - 4, LEN[contig] - 2, LEN[contig] - 1} if contig == "chr10" else {100, 150, 850, 900})
    for pos in sorted(p for p in spots if 1 <= p <= LEN[contig] - 1):
        evs = events(contig, pos)
        for ref, alt in evs:
            if pos + len(ref) - 1 > LEN[contig]:
                continue
            records.append((contig, pos, ref, [alt]))
        # records of several alts: two events at one position that share one reference
        for _ in range(3):
            (r1, a1), (r2, a2) = rng.sample(evs, 2)
            if len(r1) > len(r2):
                (r1, a1), (r2, a2) = (r2, a2), (r1, a1)
            if r2[:len(r1)] != r1:
                continue
            ref = r2
            alts = [a1 + r2[len(r1):], a2]
            if len(set(alts)) < 2 or ref in alts:
                continue
            if len({len(a) != len(ref) for a in alts}) == 2 and sum(len(a) != len(ref) for a in alts) > 1:
                continue
            if any(trims_to_a_substitution(ref, a) for a in alts):
                continue
            if pos + len(ref) - 1 > LEN[contig]:
                continue
            records.append((contig, pos, ref, alts))

records = [r for r in records if r[1] >= 1]
seen, uniq = set(), []
for r in records:
    key = (r[0], r[1], r[2], tuple(r[3]))
    if key not in seen:
        seen.add(key)
        uniq.append(r)
order = {"chr10": 0, "MT": 1}
uniq.sort(key=lambda r: (order[r[0]], r[1], r[2], r[3]))
# The very edges of each contig, which the events above leave out because they stay inside it: a
# change *before* the first base, one *after* the last, and a reference allele that runs off the
# end. Ensembl VEP writes no HGVSg for them (`hgvs_genomic` returns nothing for a variation
# feature that is not inside the slice), except for a record of several alts, which it describes
# with a position the contig does not have. They are `e<N>` and come last, with no random choice
# made for them, so that the ids and the alleles of the records above do not move.
nxt = lambda b, k=0: [x for x in "ACGT" if x != b][k]
edges = []
for contig, n in LEN.items():
    g = genome[contig]
    first, last, before_last = g[0], g[-1], g[-2]
    edges += [
        (contig, 1, first, [nxt(first) + first]),  # a base inserted before the first
        (contig, 1, first, [nxt(first, 1) + nxt(first) + first]),  # two bases
        (contig, 1, first, ["ACGTACGTA" + first]),  # nine bases
        (contig, 1, first, [first + first]),  # the first base duplicated: inside
        (contig, 1, first, [nxt(first) + first, nxt(first, 2)]),  # several alts: before, and an SNV
        (contig, 1, first, [nxt(first) + first, first + nxt(first)]),  # several alts: before, and after
        (contig, n, last, [last + nxt(last)]),  # a base inserted after the last
        (contig, n, last, [last + nxt(last) + nxt(last, 1)]),  # two bases
        (contig, n, last, [last + "ACGTACGTA"]),  # nine bases
        (contig, n, last, [last + last]),  # the last base duplicated: g.<n>dup, which VEP does not write
        (contig, n, last, [last + nxt(last), nxt(last, 2)]),  # several alts: after, and an SNV
        (contig, n, last, [nxt(last)]),  # an SNV at the last base
        (contig, n - 1, before_last + last, [before_last]),  # the last base deleted: inside
        (contig, n - 1, before_last + last + "A", [before_last]),  # a reference that runs one base off the end
        (contig, n - 1, before_last + last + "ACG", [before_last + last]),  # three bases off the end
        (contig, n - 1, before_last, [before_last + nxt(before_last)]),  # a base inserted before the last: inside
    ]

with open(f"{OUT}/in.vcf", "w") as v:
    v.write("##fileformat=VCFv4.2\n")
    for c, n in LEN.items():
        v.write(f"##contig=<ID={c},length={n}>\n")
    v.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for i, (contig, pos, ref, alts) in enumerate(uniq, 1):
        v.write(f"{contig}\t{pos}\tv{i}\t{ref}\t{','.join(alts)}\t.\t.\t.\n")
    for i, (contig, pos, ref, alts) in enumerate(edges, 1):
        v.write(f"{contig}\t{pos}\te{i}\t{ref}\t{','.join(alts)}\t.\t.\t.\n")
print(f"{len(uniq)} records and {len(edges)} at the edges of the contigs")
