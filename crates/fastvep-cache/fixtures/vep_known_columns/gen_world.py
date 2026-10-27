#!/usr/bin/env python3
"""A VEP-format variation cache with the real human 116 columns, a VCF and a
GFF3, to compare the known-variant columns of fastVEP (`AF`, `AFR_AF`, ...,
`MAX_AF`, `CLIN_SIG`, `SOMATIC`, `PHENO`, `PUBMED`) with Ensembl's own VEP run
as a program.

    gen_world.py OUTDIR [N [SEED]]

The cache has the 42 `variation_cols` of the real release-116 cache of Homo
sapiens, in their order, so every population column exists and `clin_sig_allele`,
`clin_sig_ref_allele`, `clinical_impact` and `var_synonyms` sit where VEP finds
them. Two kinds of cases share it:

* `t01` .. `t38` (and their lettered variants), written by hand, one window each between positions
  1600 and about 3500 on contig 10, inside two overlapping transcripts so that every
  allele has a row for each. Each is a corner of VEP's Perl that a plausible
  reading gets wrong (see the comment above it): the `%.4f` of AF, interpolation,
  `MAX_AF` ties, `CLIN_SIG` per allele and its strand flip, the order of the
  known variants. What VEP answers decides what is right; the comments say what
  was being probed, not what is expected.

* N generated pairs (default 400, from 10,000 on, 200 bp apart so that an input
  meets exactly one neighbourhood, outside any transcript so that each allele has
  one intergenic row): the (input site, known variant) pairs of
  ../vep_match/gen_cases.py, which are the spellings known-variant matching
  finds hard, each with a random payload: frequencies for a random subset of the
  42 columns in random formats, clinical significance per variant and per
  allele, PubMed ids, somatic and phenotype flags, failed records, a second
  record in the window.

Deterministic for a given SEED; nothing here reads the clock or the machine.

Writes OUTDIR/{in.vcf, genes.gff3, genome.fa, genome.fa.fai, cache_10.txt} and
OUTDIR/cache/homo_sapiens/116_GRCh38/info.txt. cache_10.txt still has to be
bgzipped and indexed: see regen.sh.
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
N = int(sys.argv[2]) if len(sys.argv) > 2 else 400
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 139
rng = random.Random(SEED)

COLS = ("chr,variation_name,failed,somatic,start,end,allele_string,strand,clin_sig,phenotype_or_disease,"
        "clin_sig_allele,clin_sig_ref_allele,clinical_impact,pubmed,var_synonyms,AF,AFR,AMR,EAS,EUR,SAS,"
        "gnomADe,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_NFE,gnomADe_REMAINING,"
        "gnomADe_SAS,gnomADe_MID,gnomADg,gnomADg_AFR,gnomADg_AMI,gnomADg_AMR,gnomADg_ASJ,gnomADg_EAS,"
        "gnomADg_FIN,gnomADg_MID,gnomADg_NFE,gnomADg_REMAINING,gnomADg_SAS").split(",")
assert len(COLS) == 42
FREQ_COLS = COLS[COLS.index("AF"):]

HAND_FIRST = 1600
HAND_STEP = 40
GEN_FIRST = 10000
GEN_SPACING = 200

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

LENGTH = GEN_FIRST + GEN_SPACING * (len(pairs) + 2)
genome = [rng.choice("ACGT") for _ in range(LENGTH + 1)]  # 1-based: genome[p]
genome[0] = "N"

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def rc(s):
    return "".join(COMP[c] for c in reversed(s))


def base(p):
    return genome[p]


def other(b, k=0):
    return [x for x in "ACGT" if x != b][k]


cache = []  # (start, order, line)
vcf = []  # (pos, ref, alts, id)


def rec(name, start, allele_string, end=None, failed=0, somatic=0, strand=1, **kw):
    f = {c: "." for c in COLS}
    f.update(chr="10", variation_name=name, start=str(start), allele_string=allele_string)
    f["failed"] = str(failed) if failed else "."
    f["somatic"] = str(somatic) if somatic else "."
    f["end"] = str(end) if end is not None else "."
    f["strand"] = str(strand) if strand != 1 else "."
    for k, v in kw.items():
        assert k in COLS, k
        f[k] = str(v)
    cache.append((start, len(cache), "\t".join(f[c] for c in COLS)))


def vrec(pos, ref, alts, vid):
    vcf.append((pos, ref, alts, vid))


# ---------------------------------------------------------------------------
# The hand-written cases
# ---------------------------------------------------------------------------
P = HAND_FIRST


def nxt():
    global P
    p = P
    P += HAND_STEP
    return p


# t01: formatting. AF is %.4f formatted, every other key is printed as stored.
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t01")
rec("rsT01", p, f"{r}/{a}", AF=f"{a}:0.123456", AFR=f"{a}:0.5", AMR=f"{a}:0", EAS=f"{a}:0.00005",
    EUR=f"{a}:1e-05", SAS=f"{a}:0.9999",
    gnomADe=f"{a}:8.291e-06", gnomADe_AFR=f"{a}:0.0001", gnomADe_NFE=f"{a}:0.00123456789",
    gnomADg=f"{a}:0.12", gnomADg_AFR=f"{a}:0.2")

# t02: multi-allelic input and record, per-allele frequencies
p = nxt(); r = base(p); a = other(r, 0); b = other(r, 1)
vrec(p, r, [a, b], "t02")
rec("rsT02", p, f"{r}/{a}/{b}", AF=f"{a}:0.1,{b}:0.2", AFR=f"{a}:0.01,{b}:0.02",
    gnomADg=f"{a}:0.3,{b}:0.4")

# t03: AF interpolation: only the REF allele is stored, AF is stored for 2 alleles
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t03")
rec("rsT03", p, f"{r}/{a}", AF=f"{r}:0.8765", AFR=f"{r}:0.7", gnomADg=f"{r}:0.6")

# t04: 3-allele record, only one alt has a stored AF: no interpolation for the other
p = nxt(); r = base(p); a = other(r, 0); b = other(r, 1)
vrec(p, r, [a, b], "t04")
rec("rsT04", p, f"{r}/{a}/{b}", AF=f"{a}:0.2")

# t05: record on the minus strand
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t05")
rec("rsT05", p, f"{rc(r)}/{rc(a)}", strand=-1, AF=f"{rc(a)}:0.25", AFR=f"{rc(a)}:0.125",
    clin_sig="pathogenic", clin_sig_allele=f"{rc(a)}:pathogenic", clin_sig_ref_allele=rc(r))

# t06: alt not in the record: record is not co-located at all
p = nxt(); r = base(p); a = other(r, 0); b = other(r, 1)
vrec(p, r, [a], "t06")
rec("rsT06", p, f"{r}/{b}", AF=f"{b}:0.4", somatic=0, phenotype_or_disease=1, clin_sig="benign")

# t07: several co-located: file order, ordering, SOMATIC and PHENO lists
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t07")
rec("COSV100", p, f"{r}/{a}", somatic=1, phenotype_or_disease=1)
rec("CM000100", p, "HGMD_MUTATION", phenotype_or_disease=0)
rec("rsT07b", p, f"{r}/{a}", AF=f"{a}:0.05", phenotype_or_disease=1, clin_sig="pathogenic")
rec("rsT07a", p, f"{r}/{a}", AF=f"{a}:0.07", AFR=f"{a}:0.3")
rec("esvT07", p, f"{r}/{a}", AFR=f"{a}:0.4")

# t08: clin_sig without clin_sig_allele
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t08")
rec("rsT08", p, f"{r}/{a}", clin_sig="likely_benign,benign", phenotype_or_disease=1)

# t09: clin_sig_allele per alt, multi-allelic input
p = nxt(); r = base(p); a = other(r, 0); b = other(r, 1)
vrec(p, r, [a, b], "t09")
rec("rsT09", p, f"{r}/{a}/{b}", clin_sig="pathogenic,likely_pathogenic,benign", phenotype_or_disease=1,
    clin_sig_allele=f"{a}:pathogenic;{a}:likely_pathogenic;{b}:benign", clin_sig_ref_allele=r)

# t10: the real chrY rs778145751 shape: a deletion whose clin_sig_allele key is the deleted bases
p = nxt()
for i, c in enumerate("TAAGTAAGT"):
    genome[p + i] = c
vrec(p, "TAAGTAAGT", ["TAAGT"], "t10")
rec("rsT10", p, "TAAGTAAGT/TAAGT", end=p + 8, clin_sig="pathogenic", phenotype_or_disease=1,
    clin_sig_allele="AAGT:pathogenic", clin_sig_ref_allele="AAGT", pubmed="10581029",
    gnomADe="TAAGT:2.632e-05", gnomADg="TAAGT:2.992e-05", gnomADg_NFE="TAAGT:7.423e-05")

# t10b: the same deletion with the key VEP could match ('-')
p = nxt()
for i, c in enumerate("TAAGTAAGT"):
    genome[p + i] = c
vrec(p, "TAAGTAAGT", ["TAAGT"], "t10b")
rec("rsT10b", p, "TAAGTAAGT/TAAGT", end=p + 8, clin_sig="pathogenic", phenotype_or_disease=1,
    clin_sig_allele="-:pathogenic", clin_sig_ref_allele="AAGT", gnomADg="TAAGT:2.992e-05")

# t11: failed record is not used
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t11")
rec("rsT11", p, f"{r}/{a}", failed=1, AF=f"{a}:0.5")

# t12a: MAX_AF with ties; AF, gnomADe, gnomADg do not take part
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t12a")
rec("rsT12a", p, f"{r}/{a}", AF=f"{a}:0.5", AFR=f"{a}:0.01", AMR=f"{a}:0.02", EAS=f"{a}:0.02",
    EUR=f"{a}:0", SAS=f"{a}:0.005", gnomADe=f"{a}:0.9", gnomADe_NFE=f"{a}:0.015", gnomADg=f"{a}:0.8",
    gnomADg_AMR=f"{a}:0.02")

# t12b: every frequency 0
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t12b")
rec("rsT12b", p, f"{r}/{a}", AF=f"{a}:0", AFR=f"{a}:0", gnomADg_NFE=f"{a}:0")

# t12c: known but no frequency at all
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t12c")
rec("rsT12c", p, f"{r}/{a}", phenotype_or_disease=1)

# t12d: two co-located with the same frequency: duplicates in MAX_AF_POPS?
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t12d")
rec("rsT12d1", p, f"{r}/{a}", AFR=f"{a}:0.3", AF=f"{a}:0.11114")
rec("rsT12d2", p, f"{r}/{a}", AFR=f"{a}:0.3", AF=f"{a}:0.22226")

# t12e: two co-located with different AF: order of the AF column
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t12e")
rec("rsT12e2", p, f"{r}/{a}", AF=f"{a}:0.2", EUR=f"{a}:0.2")
rec("COSV12e", p, f"{r}/{a}", somatic=1, AF=f"{a}:0.9", EUR=f"{a}:0.9")
rec("rsT12e1", p, f"{r}/{a}", AF=f"{a}:0.1", EUR=f"{a}:0.1")

# t13: pubmed on two records
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t13")
rec("rsT13a", p, f"{r}/{a}", pubmed="111,222")
rec("rsT13b", p, f"{r}/{a}", pubmed="222,333")

# t14: record without alleles (COSMIC) on a multi-allelic input: reported for every alt
p = nxt(); r = base(p); a = other(r, 0); b = other(r, 1)
vrec(p, r, [a, b], "t14")
rec("COSV14", p, "COSMIC_MUTATION", somatic=1, phenotype_or_disease=1, var_synonyms="COSMIC::COSM5433537")

# t15: gnomAD sub-population holes ('.') and a 2-alt gnomAD column on a 2-allele input
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t15")
rec("rsT15", p, f"{r}/{a}", gnomADe=f"{a}:0.001", gnomADe_AFR=".", gnomADe_NFE=f"{a}:0",
    gnomADg=f"{a}:0.002", gnomADg_AMI=f"{a}:0.5", gnomADg_MID=".")

# t16: an insertion (end = start - 1) and AF with '-'
p = nxt(); r = base(p); x = other(r)
vrec(p, r, [r + x], "t16")  # VCF anchor form; VEP minimises to -/x at p+1
rec("rsT16", p + 1, f"-/{x}", end=p, AF=f"{x}:0.25", gnomADg=f"{x}:0.03")

# t17: record on the neighbour position that the +-1 window still finds, but whose allele differs
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t17")
rec("rsT17", p + 1, f"{base(p + 1)}/{other(base(p + 1))}", AF=f"{other(base(p + 1))}:0.4")

# t18: clin_sig_allele names an allele the record does not have (real: rs1286380745 G/C, csa T:...)
p = nxt(); r = base(p); a = other(r, 0); c = other(r, 1)
vrec(p, r, [a], "t18")
rec("rsT18", p, f"{r}/{a}", clin_sig="uncertain_significance", clin_sig_allele=f"{c}:uncertain_significance",
    clin_sig_ref_allele=r, phenotype_or_disease=1)

# t19: a properly trimmed deletion whose clin_sig_allele key is the deleted bases
p = nxt()
for i, c in enumerate("TAAGT"):
    genome[p + i] = c
vrec(p, "TAAGT", ["T"], "t19")
rec("rsT19", p + 1, "AAGT/-", end=p + 4, clin_sig="pathogenic", phenotype_or_disease=1,
    clin_sig_allele="AAGT:pathogenic", clin_sig_ref_allele="AAGT", gnomADg="-:0.01")

# t19b: the same with the key '-'
p = nxt()
for i, c in enumerate("TAAGT"):
    genome[p + i] = c
vrec(p, "TAAGT", ["T"], "t19b")
rec("rsT19b", p + 1, "AAGT/-", end=p + 4, clin_sig="pathogenic", phenotype_or_disease=1,
    clin_sig_allele="-:pathogenic", clin_sig_ref_allele="AAGT", gnomADg="-:0.01")

# t10m: the real rs778145751 record, un-trimmed as in the real cache, found through a multi-allelic input
# (VEP does not minimise those), so the record's start lies inside the window
p = nxt()
for i, c in enumerate("TAAGTAAGT"):
    genome[p + i] = c
vrec(p, "TAAGTAAGT", ["TAAGT", "T"], "t10m")
rec("rsT10m", p, "TAAGTAAGT/TAAGT", end=p + 8, clin_sig="pathogenic", phenotype_or_disease=1,
    clin_sig_allele="AAGT:pathogenic", clin_sig_ref_allele="AAGT", pubmed="10581029",
    gnomADe="TAAGT:2.632e-05", gnomADg="TAAGT:2.992e-05", gnomADg_NFE="TAAGT:7.423e-05")

# t20: record whose allele_string is lower case on input (VEP upper-cases the input)
p = nxt(); r = base(p); a = other(r)
vrec(p, r.lower(), [a.lower()], "t20")
rec("rsT20", p, f"{r}/{a}", AF=f"{a}:0.33")

# t21: the AF column has both alleles for a biallelic record (as old caches did)
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t21")
rec("rsT21", p, f"{r}/{a}", AF=f"{r}:0.6,{a}:0.4", gnomADg=f"{r}:0.99,{a}:0.01")

# t22*: several co-located records with clin_sig_allele and/or clin_sig
# t22a: two records, both with clin_sig_allele for the alt: distinct strings joined with ';' in hash order
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t22a")
rec("rsU1", p, f"{r}/{a}", clin_sig="pathogenic", clin_sig_allele=f"{a}:pathogenic", clin_sig_ref_allele=r)
rec("rsU2", p, f"{r}/{a}", clin_sig="benign", clin_sig_allele=f"{a}:benign", clin_sig_ref_allele=r)
# t22b: record 1 has clin_sig only, record 2 has clin_sig_allele: the allele-specific value replaces the list
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t22b")
rec("rsV1", p, f"{r}/{a}", clin_sig="likely_benign")
rec("rsV2", p, f"{r}/{a}", clin_sig="pathogenic", clin_sig_allele=f"{a}:pathogenic", clin_sig_ref_allele=r)
# t22c: the other way round in the file
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t22c")
rec("rsW1", p, f"{r}/{a}", clin_sig="pathogenic", clin_sig_allele=f"{a}:pathogenic", clin_sig_ref_allele=r)
rec("rsW2", p, f"{r}/{a}", clin_sig="likely_benign")
# t22d: clin_sig only first, then a record whose clin_sig_allele names another allele
p = nxt(); r = base(p); a = other(r, 0); c = other(r, 1)
vrec(p, r, [a], "t22d")
rec("rsX1", p, f"{r}/{a}", clin_sig="benign")
rec("rsX2", p, f"{r}/{a}", clin_sig="pathogenic", clin_sig_allele=f"{c}:pathogenic", clin_sig_ref_allele=r)
# t22e: the mismatching one first
p = nxt(); r = base(p); a = other(r, 0); c = other(r, 1)
vrec(p, r, [a], "t22e")
rec("rsY1", p, f"{r}/{a}", clin_sig="pathogenic", clin_sig_allele=f"{c}:pathogenic", clin_sig_ref_allele=r)
rec("rsY2", p, f"{r}/{a}", clin_sig="benign")

# The cases below were added with the port, to reach what a reading of the Perl says
# is reachable and the cases above are not.

# t23: a record whose reference is another allele of the input. A reading of the Perl says the
# record is read for every allele of a site it was attached to, so that the interpolated AF
# (1 - 0.3) of its reference could land on an input allele. VEP does not attach it: matching is by
# (reference, alternate) pair, and the record's pair is not one of the input's, so no cell.
p = nxt(); r = base(p); a = other(r, 0); c = other(r, 1)
vrec(p, r, [a, c], "t23")
rec("rsT23", p, f"{c}/{a}", AF=f"{a}:0.3", AFR=f"{a}:0.2")

# t24: the same, with the stored frequency on the record's reference (no cell either)
p = nxt(); r = base(p); a = other(r, 0); c = other(r, 1)
vrec(p, r, [a, c], "t24")
rec("rsT24", p, f"{c}/{a}", AF=f"{c}:0.9", AFR=f"{c}:0.8")

# t25: a multi-allelic input where a record covers two alleles and a second record only one
p = nxt(); r = base(p); a = other(r, 0); b = other(r, 1); c = other(r, 2)
vrec(p, r, [a, b, c], "t25")
rec("rsT25a", p, f"{r}/{a}/{b}", AF=f"{a}:0.11,{b}:0.22", gnomADg=f"{a}:0.01,{b}:0.02",
    clin_sig="benign", clin_sig_allele=f"{a}:benign;{b}:pathogenic", clin_sig_ref_allele=r)
rec("rsT25b", p, f"{r}/{c}", AF=f"{c}:0.33", AFR=f"{c}:0.44", phenotype_or_disease=1)

# t26: minus-strand record and a multi-allelic input, frequencies keyed by the record's own alleles
p = nxt(); r = base(p); a = other(r, 0); b = other(r, 1)
vrec(p, r, [a, b], "t26")
rec("rsT26", p, f"{rc(r)}/{rc(a)}/{rc(b)}", strand=-1, AF=f"{rc(a)}:0.15,{rc(b)}:0.25",
    AMR=f"{rc(a)}:0.35,{rc(b)}:0.45", clin_sig_allele=f"{rc(a)}:benign;{rc(b)}:pathogenic",
    clin_sig_ref_allele=rc(r), clin_sig="benign,pathogenic")

# t27: the record's clin_sig_ref_allele is not the input's reference and the record is on the plus strand:
# VEP reverse-complements the alleles of clin_sig_allele all the same
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t27")
rec("rsT27", p, f"{r}/{a}", clin_sig="pathogenic", clin_sig_allele=f"{rc(a)}:pathogenic",
    clin_sig_ref_allele=rc(r) if rc(r) != r else other(r, 1))

# t28: pubmed that is 0, and one that is empty
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t28")
rec("rsT28", p, f"{r}/{a}", pubmed="0", phenotype_or_disease=0, clin_sig="0")

# t29: a record with no alleles first, then a normal one that has clin_sig_allele
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t29")
rec("COSV29", p, "COSMIC_MUTATION", somatic=1, phenotype_or_disease=1, clin_sig="pathogenic")
rec("rsT29", p, f"{r}/{a}", clin_sig="benign", clin_sig_allele=f"{a}:likely_benign", clin_sig_ref_allele=r)

# t30: a gnomAD sub-population with both alleles of a multi-allelic record, one 0
p = nxt(); r = base(p); a = other(r, 0); b = other(r, 1)
vrec(p, r, [a, b], "t30")
rec("rsT30", p, f"{r}/{a}/{b}", gnomADe=f"{a}:0.5,{b}:0", gnomADe_SAS=f"{a}:0,{b}:0.5",
    gnomADg_FIN=f"{a}:2e-05,{b}:0.0002", AFR=f"{b}:0.7")

# t31: the same frequency on two populations of two different known variants, to see which
# `MAX_AF_POPS` entries VEP repeats
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t31")
rec("rsT31a", p, f"{r}/{a}", EUR=f"{a}:0.25", gnomADg_NFE=f"{a}:0.25")
rec("rsT31b", p, f"{r}/{a}", EUR=f"{a}:0.25", SAS=f"{a}:0.25")

# t32: an AF of an old-style cache, where the record lists the alleles in the opposite order
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t32")
rec("rsT32", p, f"{r}/{a}", AF=f"{a}:0.4,{r}:0.6", AFR=f"{a}:0.4,{r}:0.6")

# t33: AF 0 (an allele nobody has) against a record that also lists the other allele
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t33")
rec("rsT33", p, f"{r}/{a}", AF=f"{a}:0", AFR=f"{a}:0,{r}:1")

# t34: the interpolated AF is 0 (the stored one is 1)
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t34")
rec("rsT34", p, f"{r}/{a}", AF=f"{r}:1")

# t35: a long run of decimals that rounds up at the fourth place, and an exact tie in binary
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t35")
rec("rsT35", p, f"{r}/{a}", AF=f"{a}:0.03125", AFR=f"{a}:0.09375")
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t35b")
rec("rsT35b", p, f"{r}/{a}", AF=f"{a}:0.99995")

# t36: the one way a record the input did not match for an allele can still speak for it. The record is on
# the other strand and matches the first alt through the reverse complement; its own reference allele,
# spelled as the record spells it, is the input's second alt, and AF of the record interpolates a
# frequency for exactly that allele (1 - 0.3). `--af` reads the record for the second alt, `--af_1kg`
# alone does not, so the sets disagree about this row.
p = nxt(); r = base(p); a = [x for x in "ACGT" if x not in (r, rc(r))][0]
vrec(p, r, [a, rc(r)], "t36")
rec("rsT36", p, f"{rc(r)}/{rc(a)}", strand=-1, AF=f"{rc(a)}:0.3", AFR=f"{rc(a)}:0.2")

# t37: the same with the stored frequency on the record's reference
p = nxt(); r = base(p); a = [x for x in "ACGT" if x not in (r, rc(r))][0]
vrec(p, r, [a, rc(r)], "t37")
rec("rsT37", p, f"{rc(r)}/{rc(a)}", strand=-1, AF=f"{rc(r)}:0.6")

# t38: trailing commas, which two records of the real cache have in `pubmed` (the first real-data run
# found them: Perl's `split` drops the empty field at the end, a plain split keeps it) and which the same
# `split` meets in `clin_sig` and in the pairs of a frequency column
p = nxt(); r = base(p); a = other(r)
vrec(p, r, [a], "t38")
rec("rsT38", p, f"{r}/{a}", pubmed="111,222,", clin_sig="benign,likely_benign,", phenotype_or_disease=1,
    AFR=f"{a}:0.1,", AMR=f"{a}:0.2,,", gnomADg_NFE=f"{a}:0.3,")

# ---------------------------------------------------------------------------
# The generated pairs, with random payloads
# ---------------------------------------------------------------------------
CLIN_SIGS = ["pathogenic", "likely_pathogenic", "benign", "likely_benign", "uncertain_significance",
             "risk_factor", "drug_response", "not_provided", "association"]
VALUES = ["0", "0.5", "0.123456", "1e-05", "8.291e-06", "0.0001", "0.9999", "1", "0.03125", "0.25",
          "0.333333", "2.5e-07", "0.00123456789", "0.75", "0.1", "0.9", "0.05"]


def random_value():
    if rng.random() < 0.8:
        return rng.choice(VALUES)
    return repr(round(rng.random() ** 3, rng.randint(1, 6)))


def payload(allele_string, strand):
    """Random columns for a known variant, keyed by the record's own alleles."""
    kw = {}
    alleles = allele_string.split("/") if "/" in allele_string else []
    ref, alts = (alleles[0], alleles[1:]) if alleles else (None, [])
    if rng.random() < 0.30:
        kw["phenotype_or_disease"] = 1
    if rng.random() < 0.35:
        kw["clin_sig"] = ",".join(rng.sample(CLIN_SIGS, rng.randint(1, 3)))
        if alleles and rng.random() < 0.5:
            parts = []
            for _ in range(rng.randint(1, 3)):
                key = rng.choice(alts + [ref] + (["-"] if rng.random() < 0.1 else []))
                parts.append(f"{key}:{rng.choice(CLIN_SIGS)}")
            kw["clin_sig_allele"] = ";".join(parts)
            how = rng.random()
            if how < 0.6:
                kw["clin_sig_ref_allele"] = ref
            elif how < 0.8 and all(c in COMP for c in ref):
                kw["clin_sig_ref_allele"] = rc(ref)
            elif how < 0.9:
                kw["clin_sig_ref_allele"] = rng.choice(["A", "C", "G", "T"])
    if rng.random() < 0.20:
        kw["pubmed"] = ",".join(str(rng.randint(100000, 999999)) for _ in range(rng.randint(1, 3)))
    if alleles:
        for col in FREQ_COLS:
            p_col = 0.5 if col == "AF" else 0.22
            if rng.random() >= p_col:
                continue
            listed = rng.sample(alleles, rng.randint(1, len(alleles))) if rng.random() < 0.5 \
                else [rng.choice(alts)]
            if col == "AF" and rng.random() < 0.5:
                listed = listed[:1]
            kw[col] = ",".join(f"{a}:{random_value()}" for a in listed)
    return kw


def span_end(start, ref):
    return start - 1 if ref == "-" else start + len(ref) - 1


for i, line in enumerate(pairs):
    a_alleles, a_pos, _a_strand, b_alleles, b_pos, b_strand = line.split("\t")
    a_pos, b_pos, b_strand = int(a_pos), int(b_pos), int(b_strand)
    shift = GEN_FIRST + GEN_SPACING * i + 50 - a_pos
    a_pos, b_pos = a_pos + shift, b_pos + shift

    # the input, as the VCF record a caller would write
    ref, *alts = a_alleles.split("/")
    if ref == "-" or "-" in alts:
        anchor = genome[a_pos - 1]
        vcf_pos = a_pos - 1
        ref = anchor + ("" if ref == "-" else ref)
        alts = [anchor + ("" if x == "-" else x) for x in alts]
    else:
        vcf_pos = a_pos
    if any(a == ref for a in alts) or len(set(alts)) != len(alts):
        continue
    vrec(vcf_pos, ref, alts, f"g{i}")

    # the known variant
    b_ref = b_alleles.split("/")[0]
    start, end = b_pos, span_end(b_pos, b_ref)
    somatic = 1 if rng.random() < 0.10 else 0
    failed = 1 if rng.random() < 0.03 else 0
    allele_string = b_alleles
    if rng.random() < 0.03:
        allele_string, somatic = "COSMIC_MUTATION", 1  # no alleles: matched on start and end
    rec(f"rsP{i}", start, allele_string, end=end, failed=failed, somatic=somatic,
        strand=b_strand, **payload(allele_string, b_strand))
    if rng.random() < 0.30:
        prefix = rng.choice(["rs", "CM", "CI", "CD", "COSV", "COSM", "XYZ", "esv"])
        second_somatic = 1 if prefix.startswith("COS") or rng.random() < 0.05 else 0
        rec(f"{prefix}Q{i}", start, allele_string, end=end, somatic=second_somatic,
            strand=b_strand, **payload(allele_string, b_strand))

# ---------------------------------------------------------------------------
cache_dir = os.path.join(OUT, "cache", "homo_sapiens", "116_GRCh38")
os.makedirs(os.path.join(cache_dir, "10"), exist_ok=True)
with open(os.path.join(OUT, "cache_10.txt"), "w") as f:
    for _start, _order, line in sorted(cache):
        f.write(line + "\n")
# info.txt: the real human 116 one, minus the regulatory/cell_types lines
with open(os.path.join(cache_dir, "info.txt"), "w") as f:
    f.write("species\thomo_sapiens\nassembly\tGRCh38\nsift\tb\npolyphen\tb\n"
            "source_sift\t6.2.1\nsource_genebuild\tGENCODE50\nsource_gencode\tGENCODE 50\n"
            "source_assembly\tGRCh38.p14\nsource_polyphen\t2.2.3\n"
            "variation_cols\t" + ",".join(COLS) + "\n"
            "source_dbSNP\t156\nsource_1000genomes\tphase3\nsource_gnomADe\tv4.1\nsource_gnomADg\tv4.1\n"
            "var_type\ttabix\n")

vcf.sort(key=lambda r: (r[0], r[1], r[2], r[3]))
with open(os.path.join(OUT, "in.vcf"), "w") as f:
    f.write(f"##fileformat=VCFv4.2\n##contig=<ID=10,length={LENGTH}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for pos, ref, alts, vid in vcf:
        f.write(f"10\t{pos}\t{vid}\t{ref}\t{','.join(alts)}\t.\t.\t.\n")

with open(os.path.join(OUT, "genome.fa"), "w") as f:
    f.write(">10\n")
    seq = "".join(genome[1:])
    for k in range(0, len(seq), 60):
        f.write(seq[k:k + 60] + "\n")
# the .fai VEP wants next to the FASTA: name, length, offset, bases per line, bytes per line
with open(os.path.join(OUT, "genome.fa.fai"), "w") as f:
    f.write(f"10\t{LENGTH}\t3\t60\t61\n")

# two overlapping transcripts of one gene around the hand-written cases
with open(os.path.join(OUT, "genes.gff3"), "w") as f:
    f.write("##gff-version 3\n"
            "10\ttest\tgene\t1000\t5000\t.\t+\t.\tID=gene:ENSG_X;Name=GENEX;biotype=lncRNA\n"
            "10\ttest\tlnc_RNA\t1000\t3000\t.\t+\t.\tID=transcript:ENST_1;Parent=gene:ENSG_X;biotype=lncRNA\n"
            "10\ttest\texon\t1000\t1500\t.\t+\t.\tID=exon:E1;Parent=transcript:ENST_1;rank=1\n"
            "10\ttest\texon\t2000\t3000\t.\t+\t.\tID=exon:E2;Parent=transcript:ENST_1;rank=2\n"
            "10\ttest\tlnc_RNA\t1500\t5000\t.\t+\t.\tID=transcript:ENST_2;Parent=gene:ENSG_X;biotype=lncRNA\n"
            "10\ttest\texon\t1500\t5000\t.\t+\t.\tID=exon:E3;Parent=transcript:ENST_2;rank=1\n")
print(f"{len(vcf)} VCF records, {len(cache)} cache records, genome {LENGTH} bp, last hand case at {P}")
