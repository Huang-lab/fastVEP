#!/usr/bin/env python3
"""A world of transcripts with very short introns, and changes around them, for comparing
fastVEP's tab output with Ensembl VEP's (issue #139).

    gen_world.py OUTDIR [SAMPLE] [SEED]

Ensembl VEP treats a transcript with an intron of 13 bases or fewer (`end - start <= 12`) as
having a "frameshift intron", and that changes what it says about every change near the
transcript's exons, not only about the ones in the intron:

  - a change inside the intron earns no splice term and no `intron_variant`; if the intron is
    within the coding region it is a `coding_sequence_variant`;
  - every exon is 12 bases wider on each side for the "does it touch an exon" test, which takes
    the `splice_polypyrimidine_tract_variant` away from the 3 to 12 bases before every exon
    of that transcript, and makes a change in the intron `coding` for HGVSp.

Intron lengths are chosen around that threshold, on both strands, in coding and non-coding
transcripts, with two controls that have introns of 14 bases and longer only. Every transcript
has a UTR at both ends, so a change that is not a frameshift never reaches the end of it.

The changes are, at every position within RADIUS of an exon edge and in every intron up to 40
bases, a substitution, deletions of 1, 3, 5 and 15 bases and insertions of 1 and 3, each a
single-ALT record. Every substitution is kept, and every SAMPLE-th of the others.

Writes OUTDIR/{genome.fa, genes.gff3, in.vcf}.
"""
import random
import sys

OUT = sys.argv[1]
SAMPLE = int(sys.argv[2]) if len(sys.argv) > 2 else 8
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 139
RADIUS = 18
rng = random.Random(SEED)
CHROM = "chr1"

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
STOPS = {"TAA", "TAG", "TGA"}
SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS | {"ATG"}]


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


LENGTH = 24000
seq = [rng.choice("ACGT") for _ in range(LENGTH)]

# (gene id, symbol, strand, biotype, start, exon lengths in transcript order,
#  intron lengths in transcript order): the 5' UTR is 40 bases and the 3' UTR 30 or more
GENES = [
    ("ENSG00000000101", "FSA", "+", "protein_coding", 1000, [100, 80, 80, 80, 80, 130], [1, 2, 13, 14, 25]),
    ("ENSG00000000102", "FSB", "-", "protein_coding", 4000, [100, 80, 80, 80, 130], [3, 1, 13, 40]),
    ("ENSG00000000103", "FSC", "+", "protein_coding", 7000, [100, 90, 130], [13, 30]),
    ("ENSG00000000104", "FSN", "+", "lncRNA", 10000, [90, 80, 80, 90], [2, 30, 13]),
    ("ENSG00000000105", "FSM", "-", "lncRNA", 13000, [90, 80, 90], [1, 14]),
    ("ENSG00000000106", "CTL", "+", "protein_coding", 16000, [100, 90, 90, 130], [14, 30, 60]),
    ("ENSG00000000107", "CTR", "-", "protein_coding", 19000, [100, 90, 130], [14, 20]),
]
UTR5 = 40


def layout(start, strand, exon_lens, intron_lens, coding):
    """Exons and CDS segments, ascending, as (start, end) pairs, from lengths in transcript order."""
    pos = start
    tx_exons = []
    for i, n in enumerate(exon_lens):
        tx_exons.append((pos, pos + n - 1))
        pos += n + (intron_lens[i] if i < len(intron_lens) else 0)
    total = sum(exon_lens)
    utr3 = 30 + ((total - UTR5 - 30) % 3)
    # the CDS in transcript coordinates is [UTR5, total - utr3), cut out of the exons
    cds_lo, cds_hi = UTR5, total - utr3
    exons_tx_order = tx_exons if strand == "+" else list(reversed(tx_exons))
    cds, offset = [], 0
    for s, e in exons_tx_order:
        n = e - s + 1
        lo, hi = max(cds_lo, offset), min(cds_hi, offset + n)
        if lo < hi:
            if strand == "+":
                cds.append((s + (lo - offset), s + (hi - offset) - 1))
            else:
                cds.append((e - (hi - offset) + 1, e - (lo - offset)))
        offset += n
    return sorted(tx_exons), (sorted(cds) if coding else [])


def place(intervals, strand, tseq):
    plus = tseq if strand == "+" else revcomp(tseq)
    i = 0
    for s, e in intervals:
        for p in range(s, e + 1):
            seq[p - 1] = plus[i]
            i += 1
    assert i == len(plus)


models = []
for gid, sym, strand, biotype, start, exon_lens, intron_lens in GENES:
    coding = biotype == "protein_coding"
    exons, cds = layout(start, strand, exon_lens, intron_lens, coding)
    if coding:
        n = sum(e - s + 1 for s, e in cds)
        assert n % 3 == 0, (sym, n)
        codons = ["ATG"] + [rng.choice(SENSE) for _ in range(n // 3 - 2)] + [rng.choice(sorted(STOPS))]
        place(cds, strand, "".join(codons))
    models.append((gid, sym, strand, biotype, exons, cds))

genome = "".join(seq)
with open(f"{OUT}/genome.fa", "w") as f:
    f.write(f">{CHROM}\n")
    for i in range(0, LENGTH, 60):
        f.write(genome[i:i + 60] + "\n")

with open(f"{OUT}/genes.gff3", "w") as g:
    g.write("##gff-version 3\n")
    g.write(f"##sequence-region   {CHROM} 1 {LENGTH}\n")
    exon_n = 0
    for gid, sym, strand, biotype, exons, cds in models:
        tid = "ENST" + gid[4:]
        lo, hi = exons[0][0], exons[-1][1]
        g.write(f"{CHROM}\ttest\tgene\t{lo}\t{hi}\t.\t{strand}\t.\t"
                f"ID=gene:{gid};Name={sym};biotype={biotype};gene_id={gid};version=1\n")
        ttype = "mRNA" if cds else "lnc_RNA"
        g.write(f"{CHROM}\ttest\t{ttype}\t{lo}\t{hi}\t.\t{strand}\t.\t"
                f"ID=transcript:{tid};Parent=gene:{gid};Name={sym}-001;biotype={biotype};"
                f"tag=Ensembl_canonical;transcript_id={tid};version=1\n")
        order = exons if strand == "+" else list(reversed(exons))
        for rank, (s, e) in enumerate(order, 1):
            exon_n += 1
            eid = f"ENSE{exon_n:011d}"
            g.write(f"{CHROM}\ttest\texon\t{s}\t{e}\t.\t{strand}\t.\t"
                    f"Parent=transcript:{tid};Name={eid};exon_id={eid};rank={rank};version=1\n")
        if cds:
            pid = "ENSP" + gid[4:]
            done = 0
            for s, e in (cds if strand == "+" else list(reversed(cds))):
                g.write(f"{CHROM}\ttest\tCDS\t{s}\t{e}\t.\t{strand}\t{(3 - done % 3) % 3}\t"
                        f"ID=CDS:{pid};Parent=transcript:{tid};protein_id={pid};version=1\n")
                done += e - s + 1

# the changes: around every exon edge, and whole in every short intron
positions = set()
for gid, sym, strand, biotype, exons, cds in models:
    for s, e in exons:
        for edge in (s, e):
            positions.update(range(edge - RADIUS, edge + RADIUS + 1))
    for (s1, e1), (s2, e2) in zip(exons, exons[1:]):
        if s2 - e1 - 1 <= 40:
            positions.update(range(e1 - 2, s2 + 3))

other = {"A": "C", "C": "G", "G": "T", "T": "A"}
rows = []
for p in sorted(positions):
    ref = genome[p - 1]
    follow = genome[p:p + 15]
    rows.append((p, ref, other[ref]))
    for n in (1, 3, 5, 15):
        rows.append((p, ref + follow[:n], ref))
    rows.append((p, ref, ref + "G"))
    rows.append((p, ref, ref + "ACG"))
# every substitution, which is what the polypyrimidine rule is stated over and what the sampling
# would otherwise thin to a handful, and every SAMPLE-th of the rest
kept = [r for i, r in enumerate(rows) if i % SAMPLE == 0 or (len(r[1]) == 1 and len(r[2]) == 1)]

with open(f"{OUT}/in.vcf", "w") as o:
    o.write("##fileformat=VCFv4.2\n")
    o.write(f"##contig=<ID={CHROM},length={LENGTH}>\n")
    o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for p, r, a in kept:
        o.write(f"{CHROM}\t{p}\t.\t{r}\t{a}\t.\t.\t.\n")
print(f"{len(kept)} of {len(rows)} records at {len(positions)} positions")
