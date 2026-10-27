#!/usr/bin/env python3
"""A world of small transcripts and the deletions that cover them, for comparing fastVEP's tab
output with Ensembl VEP's (issue #139).

    gen_world.py OUTDIR [SEED]

Ensembl calls a deletion (an alternate shorter than the reference) whose span covers a whole
transcript `transcript_ablation`, HIGH, and gives it no other term. A small-variant predictor that
asks the regions the span reaches one by one says `stop_retained_variant,5_prime_UTR_variant,
3_prime_UTR_variant` for it, which is LOW.

The transcripts are 21 to 130 bases, on both strands, coding and not, in one exon and in three. For
each, a deletion at every combination of start from 12 bases before it to 3 bases in and end from 3
bases before its last base to 12 after, so that spans that cover it exactly, by a base, by many, and
those that stop a base short are all there; and, over each transcript, a replacement by something
longer, a replacement by something as long (a substitution of every base), an insertion and a
shorter replacement. Each is a single-ALT record.

Writes OUTDIR/{genome.fa, genes.gff3, in.vcf}.
"""
import random
import sys

OUT = sys.argv[1]
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 139
rng = random.Random(SEED)
CHROM = "chr1"
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
STOPS = ["TAA", "TAG", "TGA"]
SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
         if a + b + c not in STOPS + ["ATG"]]


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


# id, strand, coding, exon lengths in transcript order, intron length, 5' UTR length
MODELS = [
    ("A1", "+", True, [60], 0, 9),
    ("A2", "-", True, [60], 0, 9),
    ("A3", "+", False, [40], 0, 0),
    ("A4", "-", False, [40], 0, 0),
    ("A5", "+", True, [40, 30, 30], 20, 9),
    ("A6", "-", True, [40, 30, 30], 20, 9),
    ("A7", "+", False, [21], 0, 0),
    ("A8", "-", True, [21], 0, 3),
]
SPACING = 400
LENGTH = SPACING * (len(MODELS) + 1)
genome = [rng.choice("ACGT") for _ in range(LENGTH)]
gff, spans, exon_n = [], [], 0
for k, (tid, strand, coding, lens, intron, utr5) in enumerate(MODELS):
    start = SPACING * (k + 1)
    exons, pos = [], start
    for n in lens:
        exons.append((pos, pos + n - 1))
        pos += n + intron
    total = sum(lens)
    if coding:
        n_cds = (total - utr5 - 3) // 3 * 3
        cds = "ATG" + "".join(rng.choice(SENSE) for _ in range(n_cds // 3 - 2)) + rng.choice(STOPS)
        mrna = "".join(rng.choice("CT") for _ in range(utr5)) + cds
        mrna += "".join(rng.choice("GC") for _ in range(total - len(mrna)))
        cds_lo, cds_hi = utr5, utr5 + len(cds)
    else:
        mrna = "".join(rng.choice("ACGT") for _ in range(total))
        cds_lo = cds_hi = 0
    plus = mrna if strand == "+" else revcomp(mrna)
    i = 0
    for s, e in exons:
        for p in range(s, e + 1):
            genome[p - 1] = plus[i]
            i += 1
    lo, hi = exons[0][0], exons[-1][1]
    spans.append((tid, lo, hi))
    gid, ens_t, ens_p = (f"ENSG00000000{k + 301}", f"ENST00000000{k + 301}", f"ENSP00000000{k + 301}")
    biotype = "protein_coding" if coding else "lncRNA"
    gff.append(f"{CHROM}\ttest\tgene\t{lo}\t{hi}\t.\t{strand}\t.\t"
               f"ID=gene:{gid};Name={tid};biotype={biotype};gene_id={gid};version=1")
    gff.append(f"{CHROM}\ttest\t{'mRNA' if coding else 'lnc_RNA'}\t{lo}\t{hi}\t.\t{strand}\t.\t"
               f"ID=transcript:{ens_t};Parent=gene:{gid};Name={tid}-001;biotype={biotype};"
               f"tag=Ensembl_canonical;transcript_id={ens_t};version=1")
    order = exons if strand == "+" else list(reversed(exons))
    for rank, (s, e) in enumerate(order, 1):
        exon_n += 1
        eid = f"ENSE{exon_n:011d}"
        gff.append(f"{CHROM}\ttest\texon\t{s}\t{e}\t.\t{strand}\t.\t"
                   f"Parent=transcript:{ens_t};Name={eid};exon_id={eid};rank={rank};version=1")
    if coding:
        done, offset = 0, 0
        for s, e in order:
            n = e - s + 1
            a, b = max(cds_lo, offset), min(cds_hi, offset + n)
            if a < b:
                g_lo = s + (a - offset) if strand == "+" else e - (b - offset) + 1
                g_hi = s + (b - offset) - 1 if strand == "+" else e - (a - offset)
                gff.append(f"{CHROM}\ttest\tCDS\t{g_lo}\t{g_hi}\t.\t{strand}\t{(3 - done % 3) % 3}\t"
                           f"ID=CDS:{ens_p};Parent=transcript:{ens_t};protein_id={ens_p};version=1")
                done += b - a
            offset += n

G = "".join(genome)
rows = set()
for tid, lo, hi in spans:
    for ds in (-12, -3, -1, 0, 1, 3):
        for de in (-3, -1, 0, 1, 3, 12):
            a, b = lo + ds, hi + de
            if a < 2 or b <= a:
                continue
            # a deletion of [a, b], written with the base in front of it
            rows.add((a - 1, G[a - 2:b], G[a - 2]))
    # over the whole transcript: longer, as long, shorter, and an insertion at its start
    a, b = lo - 3, hi + 3
    ref = G[a - 2:b]
    rows.add((a - 1, ref, G[a - 2] + "".join(rng.choice("ACGT") for _ in range(len(ref) + 5))))
    rows.add((a - 1, ref, G[a - 2] + "".join(rng.choice("ACGT") for _ in range(len(ref) - 1))))
    rows.add((a - 1, ref, G[a - 2] + "".join(rng.choice("ACGT") for _ in range(len(ref) // 2))))
    rows.add((lo, G[lo - 1], G[lo - 1] + "ACGTA"))
    for p in range(lo, hi + 1, 7):
        rows.add((p, G[p - 1], {"A": "C", "C": "G", "G": "T", "T": "A"}[G[p - 1]]))
kept = sorted(rows)

with open(f"{OUT}/genome.fa", "w") as f:
    f.write(f">{CHROM}\n")
    for i in range(0, LENGTH, 60):
        f.write(G[i:i + 60] + "\n")
with open(f"{OUT}/genes.gff3", "w") as g:
    g.write("##gff-version 3\n")
    g.write(f"##sequence-region   {CHROM} 1 {LENGTH}\n")
    g.write("\n".join(gff) + "\n")
with open(f"{OUT}/in.vcf", "w") as o:
    o.write("##fileformat=VCFv4.2\n")
    o.write(f"##contig=<ID={CHROM},length={LENGTH}>\n")
    o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for p, r, a in kept:
        o.write(f"{CHROM}\t{p}\t.\t{r}\t{a}\t.\t.\t.\n")
print(f"{len(kept)} records")
