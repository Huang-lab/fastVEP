#!/usr/bin/env python3
"""A VEP-format variation cache and a VCF for an end-to-end comparison of
known-variant matching with Ensembl's own VEP, run as a program.

    gen_world.py OUTDIR N [SEED]

Each (input site, known variant) pair of ../vep_match/gen_cases.py gets a window of
its own on chr10, 200 bp from the next, so that an input meets exactly one cache
neighbourhood. The input becomes a VCF record, anchored on the left when an allele is
empty, the way a caller writes it. The known variant becomes a record of all_vars.gz
with VEP's coordinates (an insertion has end = start - 1). Some cache records are
somatic or failed, some carry no alleles at all (COSMIC_MUTATION), and some windows
get a second, differently named record, so that the ordering of Existing_variation is
exercised too. Every VCF record has an ID of its own (`in<N>`), which is how the
answers are matched back to the input.

Writes OUTDIR/{in.vcf, genes.gff3, genome.fa, genome.fa.fai, cache_10.txt} and the
cache's info.txt and chr_synonyms.txt under OUTDIR/cache/homo_sapiens/116_GRCh38/.
cache_10.txt still has to be bgzipped and indexed: see regen.sh.

Deterministic for a given SEED; nothing here reads the clock or the machine.
"""
import importlib.util
import os
import random
import sys

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("gen_cases", os.path.join(here, "..", "vep_match", "gen_cases.py"))
gc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gc)

OUT = sys.argv[1]
N = int(sys.argv[2])
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 139
rng = random.Random(SEED)

pairs = list(gc.handwritten())
n_struct, n_multi, n_rand = int(N * 0.6), int(N * 0.25), int(N * 0.15)
generated = gc.structured(rng, n_struct) + gc.multiallelic(rng, n_multi) + gc.random_pairs(rng, n_rand)
rng.shuffle(generated)
seen = set(pairs)
for line in generated:
    if line not in seen:
        seen.add(line)
        pairs.append(line)
pairs = pairs[:N]

SPACING = 200
FIRST = 1000
LENGTH = FIRST + SPACING * (len(pairs) + 2)
genome = [rng.choice("ACGT") for _ in range(LENGTH)]


def span_end(start, ref):
    return start - 1 if ref == "-" else start + len(ref) - 1


vcf = []
cache = []  # (start, order, line)
for i, line in enumerate(pairs):
    a_alleles, a_pos, _a_strand, b_alleles, b_pos, b_strand = line.split("\t")
    a_pos, b_pos, b_strand = int(a_pos), int(b_pos), int(b_strand)
    shift = FIRST + SPACING * i + 50 - a_pos
    a_pos, b_pos = a_pos + shift, b_pos + shift

    # the input, as the VCF record a caller would write
    ref, *alts = a_alleles.split("/")
    if ref == "-" or "-" in alts:
        anchor = genome[a_pos - 2]
        vcf_pos = a_pos - 1
        ref = anchor + ("" if ref == "-" else ref)
        alts = [anchor + ("" if x == "-" else x) for x in alts]
    else:
        vcf_pos = a_pos
    if any(a == ref for a in alts) or len(set(alts)) != len(alts):
        continue
    vcf.append((vcf_pos, ref, alts, f"in{i}"))

    # the known variant
    b_ref = b_alleles.split("/")[0]
    start, end = b_pos, span_end(b_pos, b_ref)
    somatic = 1 if rng.random() < 0.10 else 0
    failed = 1 if rng.random() < 0.03 else 0
    allele_string = b_alleles
    if rng.random() < 0.03:
        allele_string, somatic = "COSMIC_MUTATION", 1  # no alleles: matched on start and end
    rest = "\t.\t.\t.\t0\t.\t.\t.\t.\t.\t.\t.\t."
    cache.append((start, len(cache), f"10\trsP{i}\t{failed}\t{somatic}\t{start}\t{end}\t{allele_string}\t{b_strand}{rest}"))
    if rng.random() < 0.12:
        prefix = rng.choice(["rs", "CM", "CI", "CD", "COSV", "COSM", "XYZ", "esv"])
        second_somatic = 1 if prefix.startswith("COS") or rng.random() < 0.05 else 0
        cache.append((start, len(cache), f"10\t{prefix}Q{i}\t0\t{second_somatic}\t{start}\t{end}\t{allele_string}\t{b_strand}{rest}"))

cache_dir = os.path.join(OUT, "cache", "homo_sapiens", "116_GRCh38")
os.makedirs(os.path.join(cache_dir, "10"), exist_ok=True)
with open(os.path.join(OUT, "cache_10.txt"), "w") as f:
    for _start, _order, line in sorted(cache):
        f.write(line + "\n")
with open(os.path.join(cache_dir, "info.txt"), "w") as f:
    f.write("# CACHE UPDATED 2026-04-08 00:00:00\nspecies\thomo_sapiens\nassembly\tGRCh38\ncache_region_size\t1000000\n"
            "variation_cols\tchr,variation_name,failed,somatic,start,end,allele_string,strand,minor_allele,"
            "minor_allele_freq,clin_sig,phenotype_or_disease,clin_sig_allele,pubmed,var_synonyms,AFR,AMR,EAS,EUR,SAS\n"
            "var_type\ttabix\n")
with open(os.path.join(cache_dir, "chr_synonyms.txt"), "w") as f:
    f.write("10\tchr10\tNC_000010.11\n")

vcf.sort(key=lambda r: (r[0], r[1], r[2], r[3]))
with open(os.path.join(OUT, "in.vcf"), "w") as f:
    f.write(f"##fileformat=VCFv4.2\n##contig=<ID=chr10,length={LENGTH}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for pos, ref, alts, vid in vcf:
        f.write(f"chr10\t{pos}\t{vid}\t{ref}\t{','.join(alts)}\t.\t.\t.\n")

with open(os.path.join(OUT, "genome.fa"), "w") as f:
    f.write(">chr10\n")
    seq = "".join(genome)
    for k in range(0, len(seq), 60):
        f.write(seq[k:k + 60] + "\n")
# the .fai VEP wants next to the FASTA: name, length, offset, bases per line, bytes per line
with open(os.path.join(OUT, "genome.fa.fai"), "w") as f:
    f.write(f"chr10\t{LENGTH}\t7\t60\t61\n")

# a gene near the start, so that the contig is a known one to VEP's GFF mode
with open(os.path.join(OUT, "genes.gff3"), "w") as f:
    f.write("##gff-version 3\n"
            "chr10\ttest\tgene\t100\t150\t.\t+\t.\tID=gene:ENSG_X;Name=GENEX;biotype=lncRNA\n"
            "chr10\ttest\tlnc_RNA\t100\t150\t.\t+\t.\tID=transcript:ENST_X;Parent=gene:ENSG_X;biotype=lncRNA\n"
            "chr10\ttest\texon\t100\t150\t.\t+\t.\tID=exon:E1;Parent=transcript:ENST_X;rank=1\n")
print(f"{len(vcf)} VCF records, {len(cache)} cache records, genome {LENGTH} bp")
