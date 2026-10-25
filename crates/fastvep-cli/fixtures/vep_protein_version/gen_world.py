#!/usr/bin/env python3
"""Five coding genes whose proteins carry different versions, and missense variants in each.

    gen_world.py OUTDIR

Ensembl's GFF3 writes a protein's version on its CDS records, `protein_id=ENSP00000363557;version=3`,
and one protein in six of Ensembl 115's has one that is not 1. Ensembl VEP 116 reads it into the
accession of every HGVSp (`ENSP00000363557.3:p.Thr266Met`), and fastVEP wrote `.1` whatever it said.
NCBI's GFF3 has no `version` attribute and carries the version in the id (`protein_id=NP_000001.2`),
which VEP writes as it is, and fastVEP wrote `NP_000001.2.1`.

`chrP` is a random contig. Each gene is one exon of 300 bases, `ATG`, 98 codons that are not stops and a
stop, on the forward strand. The proteins:

    ENSP00000000011  version=3
    ENSP00000000012  version=1     written out, which reads the same as no version
    ENSP00000000013  (none)        no attribute: 1
    ENSP00000000014  version=12    two digits
    NP_000000015.2   NCBI-style records (`rna-`, `gene-`), the version in the id

Variants: five missense SNVs in each, at codons spread through the protein, and a synonymous one.

Writes OUTDIR/{genome.fa, genes.gff3, in.vcf}. Every record has an ID, `v<gene>_<codon>`.
"""
import random
import sys

OUT = sys.argv[1]
rng = random.Random(139)

BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AMINO[16 * i + 4 * j + k]
         for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
SENSE = [c for c, aa in CODON.items() if aa != "*"]

contig = [rng.choice("ACGT") for _ in range(6000)]
genes = [
    ("ENSG00000000011", "ENST00000000011", "ENSP00000000011", 1000, 3),
    ("ENSG00000000012", "ENST00000000012", "ENSP00000000012", 2000, 1),
    ("ENSG00000000013", "ENST00000000013", "ENSP00000000013", 3000, None),
    ("ENSG00000000014", "ENST00000000014", "ENSP00000000014", 4000, 12),
    ("gene-NB15", "NM_000015.3", "NP_000015.2", 5000, "ncbi"),
]
gff = ["##gff-version 3", "##sequence-region   chrP 1 6000"]
cds_of = {}
for gid, tid, pid, start, version in genes:
    codons = ["ATG"] + [rng.choice(SENSE) for _ in range(98)] + ["TAA"]
    cds = "".join(codons)
    end = start + len(cds) - 1
    contig[start - 1:end] = cds
    cds_of[gid] = (start, codons)
    if version == "ncbi":
        # as RefSeq's GFF3 writes it: prefixed ids, the transcript and protein versions in the ids
        gff += [
            f"chrP\tBestRefSeq\tgene\t{start}\t{end}\t.\t+\t.\tID={gid};Name=NB15;gbkey=Gene;gene=NB15;gene_biotype=protein_coding",
            f"chrP\tBestRefSeq\tmRNA\t{start}\t{end}\t.\t+\t.\tID=rna-{tid};Parent={gid};gbkey=mRNA;gene=NB15;transcript_id={tid}",
            f"chrP\tBestRefSeq\texon\t{start}\t{end}\t.\t+\t.\tID=exon-{tid}-1;Parent=rna-{tid};gbkey=mRNA;gene=NB15;transcript_id={tid}",
            f"chrP\tBestRefSeq\tCDS\t{start}\t{end}\t.\t+\t0\tID=cds-{pid};Parent=rna-{tid};gbkey=CDS;gene=NB15;protein_id={pid}",
        ]
        continue
    v = f";version={version}" if version is not None else ""
    gff += [
        f"chrP\ttest\tgene\t{start}\t{end}\t.\t+\t.\tID=gene:{gid};Name=GENE{gid[-2:]};biotype=protein_coding;gene_id={gid};version=1",
        f"chrP\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\tID=transcript:{tid};Parent=gene:{gid};Name=GENE{gid[-2:]}-201;biotype=protein_coding;tag=Ensembl_canonical;transcript_id={tid};version=2",
        f"chrP\ttest\texon\t{start}\t{end}\t.\t+\t.\tID=exon:{tid}.1;Parent=transcript:{tid};Name={tid}.1;exon_id={tid}.1;rank=1;version=1",
        f"chrP\ttest\tCDS\t{start}\t{end}\t.\t+\t0\tID=CDS:{pid};Parent=transcript:{tid};protein_id={pid}{v}",
    ]

with open(f"{OUT}/genes.gff3", "w") as f:
    f.write("\n".join(gff) + "\n")
with open(f"{OUT}/genome.fa", "w") as f:
    f.write(">chrP\n")
    s = "".join(contig)
    for i in range(0, len(s), 60):
        f.write(s[i:i + 60] + "\n")

records = []
for gid, tid, pid, start, version in genes:
    _, codons = cds_of[gid]
    for k in (5, 25, 50, 75, 98):
        codon = codons[k]
        alt = None
        for pos in rng.sample(range(3), 3):
            for base in "ACGT":
                if base == codon[pos]:
                    continue
                new = codon[:pos] + base + codon[pos + 1:]
                if CODON[new] not in (CODON[codon], "*"):
                    alt = (pos, base)
                    break
            if alt:
                break
        pos, base = alt
        genomic = start + 3 * k + pos
        records.append((genomic, f"v{gid[-2:]}_{k}", codon[pos], base))
    # a synonymous one: third base of codon 40
    codon = codons[40]
    for base in "ACGT":
        if base != codon[2] and CODON[codon[:2] + base] == CODON[codon]:
            records.append((start + 3 * 40 + 2, f"v{gid[-2:]}_syn", codon[2], base))
            break

records.sort()
with open(f"{OUT}/in.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=chrP,length=6000>\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for pos, rid, ref, alt in records:
        f.write(f"chrP\t{pos}\t{rid}\t{ref}\t{alt}\t.\t.\t.\n")
print(len(records), "records")
