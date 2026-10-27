#!/usr/bin/env python3
"""A world of transcripts with deliberately built start codons, and every kind of change around
them, for comparing fastVEP's tab output with Ensembl VEP's (issue #139).

    gen_world.py OUTDIR [SEED]

What decides `start_lost` is the sequence the edit leaves at the start of the coding sequence, and
Ensembl VEP decides it from three tests that each fail in a way a real transcript shows:

  - the edit-based test (`_ins_del_start_altered`) asks "is the ATG still at the coding start"
    only for a transcript that has a 5' UTR, so on one that begins at its initiator (T03, T04,
    T05) an insertion that keeps the ATG is a start loss;
  - the peptide test compares the window of the change with its own alternate by prefix and
    suffix, which a deletion cannot satisfy, so the loss of codon 2 is a start loss beside a
    `start_retained_variant` (T01, T02, T05, where the second codon repeats its neighbours so that
    the deletion has a leftmost form that touches the start codon);
  - and what a start loss is *written* as, `p.Met1?`, becomes a range or two residue names when
    the window is longer than a residue.

The transcripts (cDNA = 5' UTR + CDS + 3' UTR, all in one exon unless it says otherwise):

  T01 +  UTR 10  ATG GCG GCT GCT ...        codon 2 repeats the last base of ATG
  T02 -  UTR 10  the same on the reverse strand
  T03 +  UTR 0   ATG GCT GCT ...            a transcript that begins at its initiator
  T04 -  UTR 0   the same on the reverse strand
  T05 +  UTR 0   ATG GCG GCT ...            both
  T06 +  UTR 10  ATG ATG GCT ...            the initiator repeated
  T07 -  UTR 10  ATG ATG GCT ...
  T08 +  UTR 2   ATG GCT ...                a short UTR
  T09 +  two exons, the first ending on the A of ATG (the start codon is split by an intron)
  T10 -  the same on the reverse strand

At every position from 12 bases before the coding start to 14 bases after it: three
substitutions, deletions of 1 to 6, 9 and 12 bases, insertions of 1 to 6 bases (a duplication of
what precedes, and bases chosen at random) and replacements of two and three bases, each a
single-ALT record.

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


def codons(first, second, n=10):
    """The CDS: `first` and `second` codons, then random sense codons, then a stop."""
    body = [first, second] + [rng.choice(SENSE) for _ in range(n)]
    return "".join(body) + rng.choice(STOPS)


def utr(n):
    return "".join(rng.choice("CT") for _ in range(n))  # no ATG can be made of it


# id, strand, 5' UTR length, CDS, exon lengths in transcript order (None: one exon)
MODELS = [
    ("T01", "+", 10, codons("ATG", "GCG"), None),
    ("T02", "-", 10, codons("ATG", "GCG"), None),
    ("T03", "+", 0, codons("ATG", "GCT"), None),
    ("T04", "-", 0, codons("ATG", "GCT"), None),
    ("T05", "+", 0, codons("ATG", "GCG"), None),
    ("T06", "+", 10, codons("ATG", "ATG"), None),
    ("T07", "-", 10, codons("ATG", "ATG"), None),
    ("T08", "+", 2, codons("ATG", "GCT"), None),
    ("T09", "+", 9, codons("ATG", "GCT"), [10, None]),
    ("T10", "-", 9, codons("ATG", "GCT"), [10, None]),
]
SPACING = 1000
LENGTH = SPACING * (len(MODELS) + 1)
genome = [rng.choice("ACGT") for _ in range(LENGTH)]


def place(intervals, strand, tseq):
    plus = tseq if strand == "+" else revcomp(tseq)
    i = 0
    for s, e in intervals:
        for p in range(s, e + 1):
            genome[p - 1] = plus[i]
            i += 1
    assert i == len(plus)


records = []  # (transcript, position, ref, alt)
gff = []
exon_n = 0
for k, (tid, strand, utr5, cds, exon_lens) in enumerate(MODELS):
    start = SPACING * (k + 1)
    mrna = utr(utr5) + cds + "".join(rng.choice("GC") for _ in range(40))
    if exon_lens is None:
        exons = [(start, start + len(mrna) - 1)]
    else:
        first = exon_lens[0]
        intron = 50
        exons = [(start, start + first - 1),
                 (start + first + intron, start + first + intron + len(mrna) - first - 1)]
    place(exons, strand, mrna)
    # the CDS in genomic coordinates
    tx_exons = exons if strand == "+" else list(reversed(exons))
    cds_lo, cds_hi = utr5, utr5 + len(cds)
    cds_iv, offset = [], 0
    for s, e in tx_exons:
        n = e - s + 1
        lo, hi = max(cds_lo, offset), min(cds_hi, offset + n)
        if lo < hi:
            if strand == "+":
                cds_iv.append((s + (lo - offset), s + (hi - offset) - 1))
            else:
                cds_iv.append((e - (hi - offset) + 1, e - (lo - offset)))
        offset += n
    cds_iv.sort()
    gid = f"ENSG00000000{k + 201}"
    ens_t = f"ENST00000000{k + 201}"
    ens_p = f"ENSP00000000{k + 201}"
    lo, hi = exons[0][0], exons[-1][1]
    gff.append(f"{CHROM}\ttest\tgene\t{lo}\t{hi}\t.\t{strand}\t.\t"
               f"ID=gene:{gid};Name={tid};biotype=protein_coding;gene_id={gid};version=1")
    gff.append(f"{CHROM}\ttest\tmRNA\t{lo}\t{hi}\t.\t{strand}\t.\t"
               f"ID=transcript:{ens_t};Parent=gene:{gid};Name={tid}-001;biotype=protein_coding;"
               f"tag=Ensembl_canonical;transcript_id={ens_t};version=1")
    order = exons if strand == "+" else list(reversed(exons))
    for rank, (s, e) in enumerate(order, 1):
        exon_n += 1
        eid = f"ENSE{exon_n:011d}"
        gff.append(f"{CHROM}\ttest\texon\t{s}\t{e}\t.\t{strand}\t.\t"
                   f"Parent=transcript:{ens_t};Name={eid};exon_id={eid};rank={rank};version=1")
    done = 0
    for s, e in (cds_iv if strand == "+" else list(reversed(cds_iv))):
        gff.append(f"{CHROM}\ttest\tCDS\t{s}\t{e}\t.\t{strand}\t{(3 - done % 3) % 3}\t"
                   f"ID=CDS:{ens_p};Parent=transcript:{ens_t};protein_id={ens_p};version=1")
        done += e - s + 1
    # the genomic position of the first base of the CDS, and the window of changes around it
    first_cds = cds_iv[0][0] if strand == "+" else cds_iv[-1][1]
    step = 1 if strand == "+" else -1
    window = sorted({first_cds + step * d for d in range(-12, 15)})
    records.append((tid, window, exons))

other = {"A": "C", "C": "G", "G": "T", "T": "A"}
rows = set()
G = "".join(genome)
for tid, window, exons in records:
    inside = lambda p, n=1: any(s <= p and p + n - 1 <= e for s, e in exons)
    for p in window:
        ref = G[p - 1]
        follow = G[p:p + 15]
        for a in "ACGT":
            if a != ref:
                rows.add((p, ref, a))
        for n in (1, 2, 3, 4, 5, 6, 9, 12):
            if len(follow) >= n:
                rows.add((p, ref + follow[:n], ref))
        for n in (1, 2, 3, 4, 5, 6):
            rows.add((p, ref, ref + "".join(rng.choice("ACGT") for _ in range(n))))
            if p - n >= 1:
                rows.add((p, ref, ref + G[p - n:p]))      # a duplication of what precedes
        for n in (2, 3):
            if len(follow) >= n - 1:
                new = "".join(rng.choice("ACGT") for _ in range(n))
                old = ref + follow[:n - 1]
                if new != old:
                    rows.add((p, old, new))
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
