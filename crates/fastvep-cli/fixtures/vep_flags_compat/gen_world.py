#!/usr/bin/env python3
"""A RefSeq-style world for the VEP flags that change what text is written (`--no_escape`,
`--total_length`), and changes in and around each of its genes.

    gen_world.py OUTDIR

`chrR` is a random contig of 40 kb. Its genes, in NCBI's GFF3 dress (what Ensembl VEP's `--gff`
reader reads as a RefSeq GFF3):

  NB1      + NM_090001.3: coding, two exons
  NB2      - NM_090002.2: coding, two exons, on the reverse strand
  NB3      + NR_090003.1: non-coding, two exons
  GENE A   + NM_090004.1: coding, one exon, and a symbol with a space in it (`Name=GENE A`)
  (nothing)  from 32,000 to 40,000

The changes are, at every position of a set around and in each transcript (its exon edges, its
exons' interiors, its introns, its flanks up to 3 kb away), a substitution, a deletion of one base
and an insertion of one, each a single-ALT record with an ID of its own, and a synonymous change in
each of the first codons of every coding transcript: the `=` of `p.Val3=` that `--no_escape` is about.

Writes OUTDIR/{genome.fa, genes.gff3, in.vcf}.
"""
import random
import sys

OUT = sys.argv[1]
rng = random.Random(139)
LENGTH = 40000

BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AMINO[16 * i + 4 * j + k]
         for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
SENSE = [c for c, aa in CODON.items() if aa != "*"]
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


contig = [rng.choice("ACGT") for _ in range(LENGTH)]

# (symbol, GeneID, strand, gene_biotype, accession, source, transcript type, exons, coding span, protein)
GENES = [
    ("NB1", 90001, "+", "protein_coding", "NM_090001.3", "BestRefSeq", "mRNA",
     [(1000, 1299), (1600, 1899)], (1100, 1849), "NP_090001.2"),
    ("NB2", 90002, "-", "protein_coding", "NM_090002.2", "BestRefSeq", "mRNA",
     [(6000, 6299), (6600, 6899)], (6100, 6849), "NP_090002.2"),
    ("NB3", 90003, "+", "lncRNA", "NR_090003.1", "BestRefSeq", "lnc_RNA",
     [(20000, 20199), (20400, 20599)], None, None),
    ("GENE A", 90004, "+", "protein_coding", "NM_090004.1", "BestRefSeq", "mRNA",
     [(24000, 24299)], (24100, 24249), "NP_090004.1"),
]


def coding_pieces(exons, cds):
    """The coding span's pieces, ascending: each exon cut to it."""
    lo, hi = cds
    return [(max(s, lo), min(e, hi)) for s, e in exons if max(s, lo) <= min(e, hi)]


def place_cds(pieces, strand):
    """Write a coding sequence into the contig: ATG, sense codons, a stop codon."""
    n = sum(b - a + 1 for a, b in pieces)
    assert n % 3 == 0, (pieces, n)
    codons = ["ATG"] + [rng.choice(SENSE) for _ in range(n // 3 - 2)] + [rng.choice(["TAA", "TAG", "TGA"])]
    plus = "".join(codons) if strand == "+" else revcomp("".join(codons))
    i = 0
    for a, b in pieces:
        for p in range(a, b + 1):
            contig[p - 1] = plus[i]
            i += 1


for symbol, geneid, strand, gbio, acc, source, ttype, exons, cds, prot in GENES:
    if cds:
        place_cds(coding_pieces(exons, cds), strand)

gff = ["##gff-version 3"]
for symbol, geneid, strand, gbio, acc, source, ttype, exons, cds, prot in GENES:
    lo, hi = exons[0][0], exons[-1][1]
    gid = f"gene-{symbol.replace(' ', '_')}"
    gff.append(f"chrR\t{source}\tgene\t{lo}\t{hi}\t.\t{strand}\t.\tID={gid};Dbxref=GeneID:{geneid},HGNC:HGNC:{geneid};"
               f"Name={symbol};description=test gene {geneid};gbkey=Gene;gene={symbol};gene_biotype={gbio}")
    rid = f"rna-{acc}"
    gbkey = "mRNA" if ttype == "mRNA" else "ncRNA"
    gff.append(f"chrR\t{source}\t{ttype}\t{lo}\t{hi}\t.\t{strand}\t.\tID={rid};Parent={gid};"
               f"Dbxref=GeneID:{geneid},GenBank:{acc};Name={acc};gbkey={gbkey};gene={symbol};"
               f"product=test transcript {acc};transcript_id={acc}")
    order = exons if strand == "+" else list(reversed(exons))
    for n, (s, e) in enumerate(order, 1):
        gff.append(f"chrR\t{source}\texon\t{s}\t{e}\t.\t{strand}\t.\tID=exon-{acc}-{n};Parent={rid};"
                   f"Dbxref=GeneID:{geneid},GenBank:{acc};gbkey={gbkey};gene={symbol};transcript_id={acc}")
    if cds:
        pieces = coding_pieces(exons, cds)
        done = 0
        for a, b in (pieces if strand == "+" else list(reversed(pieces))):
            gff.append(f"chrR\t{source}\tCDS\t{a}\t{b}\t.\t{strand}\t{(3 - done % 3) % 3}\tID=cds-{prot};Parent={rid};"
                       f"Dbxref=GeneID:{geneid},GenBank:{prot};Name={prot};gbkey=CDS;gene={symbol};"
                       f"product=test protein;protein_id={prot}")
            done += b - a + 1

with open(f"{OUT}/genes.gff3", "w") as f:
    f.write("\n".join(gff) + "\n")
with open(f"{OUT}/genome.fa", "w") as f:
    f.write(">chrR\n")
    s = "".join(contig)
    for i in range(0, len(s), 60):
        f.write(s[i:i + 60] + "\n")
genome = "".join(contig)

positions = set()
for symbol, geneid, strand, gbio, acc, source, ttype, exons, cds, prot in GENES:
    for s, e in exons:
        for edge in (s, e):
            positions.update(range(edge - 4, edge + 5))
        positions.update(range(s + 10, e - 9, 37))
    for (s1, e1), (s2, e2) in zip(exons, exons[1:]):
        positions.update(range(e1 + 40, s2 - 40, 53))
    positions.update(range(exons[0][0] - 3000, exons[0][0], 400))
    positions.update(range(exons[-1][1] + 100, exons[-1][1] + 3000, 400))
positions.update(range(33000, 39000, 1000))

# `--gff` mode builds the 3' UTR of a reverse-strand transcript from the wrong place (divergence 5
# of docs/VEP_DIVERGENCE.md), so a frameshift near the end of NB2's coding sequence reads on into
# a UTR that is not there. No length-changing change goes in the last 170 bases of it.
NO_INDELS = range(6090, 6260)

other = {"A": "C", "C": "G", "G": "T", "T": "A"}
rows = []
for p in sorted(positions):
    if p < 2 or p > LENGTH - 3:
        continue
    ref = genome[p - 1]
    rows.append((p, ref, other[ref]))
    if p not in NO_INDELS:
        rows.append((p, genome[p - 1:p + 1], ref))
        rows.append((p, ref, ref + "G"))

# a synonymous change in the third base of codons 3, 4 and 20 of every coding transcript: `p.Xxx3=`
for symbol, geneid, strand, gbio, acc, source, ttype, exons, cds, prot in GENES:
    if not cds:
        continue
    pieces = coding_pieces(exons, cds)
    if strand == "-":
        pieces = list(reversed(pieces))
    coding = []
    for a, b in pieces:
        coding.extend(range(a, b + 1) if strand == "+" else range(b, a - 1, -1))
    for codon_index in (3, 4, 20):
        third = coding[codon_index * 3 + 2]
        codon = "".join((genome[q - 1] if strand == "+" else COMP[genome[q - 1]])
                        for q in coding[codon_index * 3:codon_index * 3 + 3])
        for b in "ACGT":
            if b != codon[2] and CODON[codon[:2] + b] == CODON[codon]:
                ref = genome[third - 1]
                alt = b if strand == "+" else COMP[b]
                rows.append((third, ref, alt))
                break

rows = sorted(set(rows))
with open(f"{OUT}/in.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write(f"##contig=<ID=chrR,length={LENGTH}>\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for n, (pos, ref, alt) in enumerate(rows, 1):
        f.write(f"chrR\t{pos}\tc{n}\t{ref}\t{alt}\t.\t.\t.\n")
print(len(rows), "records")
