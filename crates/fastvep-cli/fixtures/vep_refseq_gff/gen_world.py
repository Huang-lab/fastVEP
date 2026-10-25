#!/usr/bin/env python3
"""NCBI RefSeq-style records of every kind a human RefSeq GFF3 has, and SNVs in each.

    gen_world.py OUTDIR

What Ensembl VEP 116 does with a RefSeq GFF3 (`--gff`) is not what it does with Ensembl's, and the
difference is in what it calls things (BaseGXF.pm, GFF.pm, release/116):

  - a transcript is named by its `transcript_id` attribute, and by its `ID` (`rna-TRG-GCC1-5`, prefix
    and all) where it has none, as a tRNA record has not;
  - a gene is named by the `GeneID` of its transcript's Dbxref (`100189284`), not by its own `ID`
    (`gene-TRG-GCC1-5`);
  - a transcript's biotype is its `biotype` / `transcript_type` / `transcript_biotype` attribute,
    and failing that: `protein_coding` for an mRNA, the `ncrna_class` of an ncRNA, `miRNA` where the
    gene is described as a microRNA, `IG_<X>_gene` for a gene segment, and the `gbkey` attribute
    for the rest (`ncRNA` for an lnc_RNA, `misc_RNA` for a `transcript` record, `rRNA`, `tRNA`).

`chrR` is a random contig; each transcript is one exon of 240 bases, the mRNAs open with `ATG`, carry
79 sense codons and close with a stop. Writes OUTDIR/{genome.fa, genes.gff3, in.vcf}. Every record has an
ID, `r<N>`.
"""
import random
import sys

OUT = sys.argv[1]
rng = random.Random(116)

BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AMINO[16 * i + 4 * j + k]
         for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
SENSE = [c for c, aa in CODON.items() if aa != "*"]

contig = [rng.choice("ACGT") for _ in range(9000)]


def esc(text):
    return text.replace(",", "%2C").replace(";", "%3B")


# (symbol, GeneID, gene_biotype, gene description, source, transcript type, transcript accession,
#  transcript attributes beyond the common ones, protein accession or None)
GENES = [
    ("NB1", 90001, "protein_coding", "NCBI test gene 1", "BestRefSeq", "mRNA", "NM_090001.3",
     "gbkey=mRNA;product=test protein 1", "NP_090001.2"),
    ("LOC90002", 90002, "protein_coding", "uncharacterized LOC90002", "Gnomon", "mRNA", "XM_090002.1",
     "gbkey=mRNA;product=uncharacterized protein LOC90002", "XP_090002.1"),
    ("NB3", 90003, "lncRNA", "NCBI test long non-coding RNA 3", "BestRefSeq", "lnc_RNA", "NR_090003.1",
     "gbkey=ncRNA;product=NCBI test long non-coding RNA 3", None),
    ("LOC90004", 90004, "lncRNA", "uncharacterized LOC90004", "Gnomon", "transcript", "XR_090004.1",
     "gbkey=misc_RNA;product=uncharacterized LOC90004", None),
    ("MIR90005", 90005, "miRNA", "microRNA 90005", "BestRefSeq", "primary_transcript", "NR_090005.1",
     "gbkey=precursor_RNA;product=microRNA 90005", None),
    ("NB6", 90006, "rRNA", "NCBI test rRNA 6", "BestRefSeq", "rRNA", "NR_090006.1",
     "gbkey=rRNA;product=RNA, test pre-ribosomal 6", None),
    ("LOC90007", 90007, "snRNA", "U2 spliceosomal RNA", "cmsearch", "snRNA", "XR_090007.1",
     "gbkey=ncRNA;product=U2 spliceosomal RNA", None),
    ("NB8", 90008, "other", "NCBI test small Cajal body RNA 8", "BestRefSeq", "ncRNA", "NR_090008.1",
     "gbkey=ncRNA;ncrna_class=scaRNA;product=small Cajal body-specific RNA 8", None),
    # a record with no transcript_id of its own, as RefSeq's tRNA records have none
    ("TRN9", 90009, "tRNA", "tRNA-Gly (anticodon GCC) 9", "tRNAscan-SE", "tRNA", None,
     "gbkey=tRNA;product=tRNA-Gly", None),
]

gff = ["##gff-version 3"]
records = []
for n, (symbol, geneid, gbio, desc, source, ttype, acc, extra, prot) in enumerate(GENES):
    start = 500 + 900 * n
    end = start + 239
    if prot:
        codons = ["ATG"] + [rng.choice(SENSE) for _ in range(78)] + [rng.choice(["TAA", "TAG", "TGA"])]
        seq = "".join(codons)
    else:
        codons = None
        seq = "".join(rng.choice("ACGT") for _ in range(240))
    contig[start - 1:end] = seq
    gid = f"gene-{symbol}"
    rid = f"rna-{acc}" if acc else f"rna-{symbol}"
    tx = f"transcript_id={acc}" if acc else ""
    attrs = lambda *parts: ";".join(p for p in parts if p)
    gff.append(f"chrR\t{source}\tgene\t{start}\t{end}\t.\t+\t.\t" + attrs(
        f"ID={gid}", f"Dbxref=GeneID:{geneid},HGNC:HGNC:{geneid}", f"Name={symbol}",
        f"description={desc}", "gbkey=Gene", f"gene={symbol}", f"gene_biotype={gbio}"))
    gff.append(f"chrR\t{source}\t{ttype}\t{start}\t{end}\t.\t+\t.\t" + attrs(
        f"ID={rid}", f"Parent={gid}",
        f"Dbxref=GeneID:{geneid},GenBank:{acc}" if acc else f"Dbxref=GeneID:{geneid},HGNC:HGNC:{geneid}",
        f"Name={acc}" if acc else "", extra.split(";")[0], f"gene={symbol}", ";".join(extra.split(";")[1:]), tx))
    gff.append(f"chrR\t{source}\texon\t{start}\t{end}\t.\t+\t.\t" + attrs(
        f"ID=exon-{acc or symbol}-1", f"Parent={rid}", f"Dbxref=GeneID:{geneid}", extra.split(";")[0],
        f"gene={symbol}", tx))
    if prot:
        gff.append(f"chrR\t{source}\tCDS\t{start}\t{end}\t.\t+\t0\t" + attrs(
            f"ID=cds-{prot}", f"Parent={rid}", f"Dbxref=GeneID:{geneid},GenBank:{prot}", f"Name={prot}",
            "gbkey=CDS", f"gene={symbol}", f"product=test protein {n}", f"protein_id={prot}"))
    # two SNVs inside the exon (in different codons for the coding ones)
    for k, offset in enumerate((31, 152)):
        pos = start + offset
        ref = contig[pos - 1]
        alt = rng.choice([b for b in "ACGT" if b != ref])
        if codons:
            # a missense change: try the three bases of the codon until the residue changes and is not a stop
            codon_start = start + (offset // 3) * 3
            codon = "".join(contig[codon_start - 1:codon_start + 2])
            for p in range(3):
                for b in "ACGT":
                    if b != codon[p] and CODON[codon[:p] + b + codon[p + 1:]] not in (CODON[codon], "*"):
                        pos, ref, alt = codon_start + p, codon[p], b
                        break
                else:
                    continue
                break
        records.append((pos, f"r{n + 1}_{k + 1}", ref, alt))

# a gene's transcript is at 500, 1400, ...; keep the contig and files
with open(f"{OUT}/genes.gff3", "w") as f:
    f.write("\n".join(gff) + "\n")
with open(f"{OUT}/genome.fa", "w") as f:
    f.write(">chrR\n")
    s = "".join(contig)
    for i in range(0, len(s), 60):
        f.write(s[i:i + 60] + "\n")
records.sort()
with open(f"{OUT}/in.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=chrR,length=9000>\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for pos, rid, ref, alt in records:
        f.write(f"chrR\t{pos}\t{rid}\t{ref}\t{alt}\t.\t.\t.\n")
print(len(records), "records")
