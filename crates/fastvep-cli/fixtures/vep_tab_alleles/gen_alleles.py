#!/usr/bin/env python3
"""Alleles that are not plain bases, at places that exercise every kind of row.

    gen_alleles.py OUT.vcf

The gene models and the genome are `../vep_tab_world`'s. At each locus (a coding base on the
forward strand and on the reverse one, a UTR base, an intronic base, a base of a non-coding
transcript, and an intergenic one) the file has records whose ALT, or REF, is not a run of
A, C, G and T:

  - `*`, the spanning deletion: alone, beside a base, beside an indel;
  - `.`, which says no ALT was called;
  - an ambiguity code in the ALT, in an SNV, an MNV, an insertion that is a multiple of three
    bases and one that is not, and a replacement;
  - an ambiguity code in the REF;
  - a record with a real ALT and an ambiguous one.

Left out because VEP does not annotate it the way anything sensible would, and
docs/VEP_DIVERGENCE.md says what fastVEP does instead: `<NON_REF>` (VEP skips the whole
record, even beside a real ALT) and `<*>` (VEP writes a `coding_sequence_variant` row for
it, which says nothing). They are in `non_ref.vcf`, which the test reads on its own terms.
"""
import os
import sys

OUT = sys.argv[1]
WORLD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "vep_tab_world")
genome = "".join(l.strip() for l in open(f"{WORLD}/genome.fa") if not l.startswith(">"))


def base(p):
    return genome[p - 1]


def other(b, k=0):
    return [x for x in "ACGT" if x != b][k % 3]


# ENST00000000001 is forward (CDS 1100..3251), ENST00000000003 reverse (CDS 10099..11300),
# ENST00000000004 a non-coding forward transcript (20000..21500)
loci = [
    ("forward coding", 1150),
    ("forward coding, a codon later", 1171),
    ("reverse coding", 10200),
    ("reverse coding, a codon later", 10231),
    ("5' UTR", 900),
    ("3' UTR", 3400),
    ("intron", 1500),
    ("non-coding exon", 20100),
    ("intergenic", 15000),
]

records = []


def add(pos, ref, alt):
    records.append((pos, ref, alt))


for _, p in loci:
    r = base(p)
    o = other(r)
    # the spanning deletion
    add(p, r, "*")
    add(p + 5, base(p + 5), other(base(p + 5)) + ",*")
    add(p + 7, base(p + 7), "*," + other(base(p + 7), 1))
    add(p + 10, genome[p + 9:p + 11], base(p + 10) + ",*")  # a deletion, and the star
    # no ALT
    add(p + 13, base(p + 13), ".")
    # an ambiguity code in the ALT
    for k, code in enumerate("NRY"):
        q = p + 16 + k * 3
        add(q, base(q), code)
    q = p + 26
    add(q, genome[q - 1:q + 1], "NN")                    # MNV, both unknown
    add(q + 3, genome[q + 2:q + 4], base(q + 3) + "N")   # MNV, one unknown
    add(q + 6, base(q + 6), base(q + 6) + "N")           # a 1-base insertion
    add(q + 9, base(q + 9), base(q + 9) + "NN")          # a 2-base insertion
    add(q + 12, base(q + 12), base(q + 12) + "NNN")      # a 3-base insertion
    add(q + 16, genome[q + 15:q + 18], "NNN")            # a replacement of one codon's length
    add(q + 21, genome[q + 20:q + 23], "NN")             # a replacement that is a frameshift
    # an ambiguity code in the REF
    add(q + 26, "N", other(base(q + 26)))
    add(q + 29, "N" + base(q + 30), "N")                  # a deletion anchored on an N
    add(q + 32, "N", "N" + other(base(q + 32)))           # an insertion anchored on an N
    # a real ALT and an unknown one
    add(q + 36, base(q + 36), other(base(q + 36)) + ",N")

records = sorted(set(records))
with open(OUT, "w") as v:
    v.write("##fileformat=VCFv4.2\n##contig=<ID=chr10,length=30000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for p, ref, alt in records:
        v.write(f"chr10\t{p}\t.\t{ref}\t{alt}\t.\t.\t.\n")
print(f"{len(records)} records")
