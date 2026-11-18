#!/usr/bin/env python3
"""A world of small genes and the structural variants around them, for comparing fastVEP's tab
output with Ensembl VEP's, cell by cell.

    gen_world.py OUTDIR [SEED]

chr1 (14 kb) carries ten transcripts that between them have every shape a structural variant
rule asks about: coding on both strands with a 5' and a 3' UTR, a coding sequence that fills its
exons, a coding sequence that ends in part of a codon, a transcript with a frameshift intron (one of
ten bases), a transcript that skips an exon of its gene, a non-coding transcript of one exon and of two,
a nonsense-mediated decay transcript, a miRNA, and an antisense transcript across the exon of
another. A tandem repeat sits in an intron.

The records are the changes whose answer depends on where their two ends are. For every
transcript, the places that matter (each end of every exon and of the coding sequence and the bases
beside them, the first bases of an intron and the last, a point in the middle of it, the two ends of the
transcript and points a few hundred and a thousand bases beyond them) are taken as the ends of a
variant, in every class VEP has a rule for: `<DEL>`, `<DUP>`, `<INV>`, `<CNV>`, the tandem and mobile
element classes, the copy number alleles, and an `<INS>` at each place and of three lengths from each
of them. Every pair of places that are near each other is made, and a sample of the others.

Then the records that are about how a record is read and not about where it lies: an END that
SVLEN overrules, a type taken from SVTYPE, several ALTs, a breakend with a mate and one without, one
that Manta writes with CHR2, a tandem repeat with the INFO that VEP turns into bases, and the records
VEP skips.

Left out, because VEP's answer for them is not one that can be repeated: a breakend whose mate is on
another contig (which of its rows are written depends on the order of a Perl hash, and differs from one
run to the next), and `<*>` (VEP annotates it as an allele called `<*>`, which fastVEP does not).
`docs/STRUCTURAL_VARIANTS.md` says why.

Writes OUTDIR/{genome.fa, genes.gff3, in.vcf}.
"""
import random
import sys

OUT = sys.argv[1]
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 185
rng = random.Random(SEED)
CHROM = "chr1"
LENGTH = 14000
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
STOPS = {"TAA", "TAG", "TGA"}
SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS | {"ATG"}]


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


seq = [rng.choice("ACGT") for _ in range(LENGTH)]

# (gene id, symbol, strand, biotype, transcripts); a transcript is
# (id, canonical, exons [(start, end)] ascending, cds [(start, end)] ascending or [], biotype)
GENES = [
    ("ENSG00000000401", "SVA", "+", "protein_coding", [
        # exons 3 and 4 are a base and ten bases apart: the intron between them is a frameshift intron
        ("ENST00000000401", True, [(1000, 1129), (1300, 1420), (1600, 1700), (1711, 1900)],
         [(1050, 1129), (1300, 1420), (1600, 1700), (1711, 1849)], "protein_coding"),
        # skips exon 2
        ("ENST00000000402", False, [(1000, 1129), (1600, 1700), (1711, 1900)],
         [(1050, 1129), (1600, 1700), (1711, 1796)], "protein_coding"),
    ]),
    ("ENSG00000000402", "SVB", "-", "protein_coding", [
        ("ENST00000000403", True, [(3000, 3300), (3500, 3560), (3700, 3950)],
         [(3020, 3300), (3500, 3560), (3700, 3921)], "protein_coding"),
    ]),
    ("ENSG00000000403", "SVC", "+", "lncRNA", [
        ("ENST00000000404", True, [(5000, 5200), (5400, 5500)], [], "lncRNA"),
    ]),
    ("ENSG00000000404", "SVD", "-", "lncRNA", [
        ("ENST00000000405", True, [(6000, 6300)], [], "lncRNA"),
    ]),
    ("ENSG00000000405", "SVE", "+", "protein_coding", [
        ("ENST00000000406", True, [(7000, 7200), (7300, 7500)],
         [(7050, 7200), (7300, 7430)], "nonsense_mediated_decay"),
    ]),
    ("ENSG00000000406", "SVF", "-", "protein_coding", [
        # the coding sequence is the exons: no UTR at all
        ("ENST00000000407", True, [(8000, 8029), (8100, 8129), (8200, 8259)],
         [(8000, 8029), (8100, 8129), (8200, 8259)], "protein_coding"),
    ]),
    ("ENSG00000000407", "SVG", "+", "protein_coding", [
        # 182 coding bases, two more than a whole number of codons, and no stop codon
        ("ENST00000000408", True, [(9000, 9100), (9200, 9300)],
         [(9010, 9100), (9200, 9290)], "protein_coding"),
    ]),
    ("ENSG00000000408", "SVH", "+", "miRNA", [
        ("ENST00000000409", True, [(10000, 10090)], [], "miRNA"),
    ]),
    ("ENSG00000000409", "SVI", "-", "lncRNA", [
        # antisense, across the third exon of SVA and the frameshift intron
        ("ENST00000000410", True, [(1450, 1550), (1650, 1850)], [], "lncRNA"),
    ]),
]

# a tandem repeat in the second intron of SVB
REPEAT_UNIT, REPEAT_AT, REPEAT_COPIES = "CA", 3400, 20
for i in range(REPEAT_COPIES):
    seq[REPEAT_AT - 1 + 2 * i] = "C"
    seq[REPEAT_AT + 2 * i] = "A"


def place(intervals, strand, tseq):
    plus = tseq if strand == "+" else revcomp(tseq)
    i = 0
    for s, e in intervals:
        for pos in range(s, e + 1):
            seq[pos - 1] = plus[i]
            i += 1
    assert i == len(plus)


for gid, sym, strand, biotype, trs in GENES:
    tid, canon, exons, cds, tbiotype = trs[0]
    if not cds:
        continue
    n = sum(e - s + 1 for s, e in cds)
    if n % 3 != 0:
        # a coding sequence that ends in part of a codon has no stop codon to place
        codons = ["ATG"] + [rng.choice(SENSE) for _ in range(n // 3 - 1)]
        tseq = "".join(codons) + "".join(rng.choice("ACGT") for _ in range(n % 3))
    else:
        codons = ["ATG"] + [rng.choice(SENSE) for _ in range(n // 3 - 2)] + [rng.choice(sorted(STOPS))]
        tseq = "".join(codons)
    place(cds, strand, tseq)

# The transcript that skips an exon takes the coding sequence of the one it shares its exons with, and
# so its last codon is whatever the bases there are. VEP's `--gff` reader writes no `cds_end_NF` flag
# for a coding sequence that ends in a sense codon, and fastVEP does (it reads the end from the
# sequence), so the transcript ends in a stop of its own: in the other transcript's frame it is the
# `xTA` and `A..` of two sense codons, and no stop.
seq[1793:1796] = list("TAA")

genome = "".join(seq)
with open(f"{OUT}/genome.fa", "w") as f:
    f.write(f">{CHROM}\n")
    for i in range(0, LENGTH, 60):
        f.write(genome[i:i + 60] + "\n")

with open(f"{OUT}/genes.gff3", "w") as g:
    g.write("##gff-version 3\n")
    g.write(f"##sequence-region   {CHROM} 1 {LENGTH}\n")
    exon_n = 0
    for gid, sym, strand, biotype, trs in GENES:
        lo = min(e[0] for t in trs for e in t[2])
        hi = max(e[1] for t in trs for e in t[2])
        g.write(f"{CHROM}\ttest\tgene\t{lo}\t{hi}\t.\t{strand}\t.\t"
                f"ID=gene:{gid};Name={sym};biotype={biotype};gene_id={gid};version=1\n")
        for tid, canon, exons, cds, tbiotype in trs:
            ttype = "mRNA" if cds else "lnc_RNA"
            tag = ";tag=Ensembl_canonical" if canon else ""
            g.write(f"{CHROM}\ttest\t{ttype}\t{exons[0][0]}\t{exons[-1][1]}\t.\t{strand}\t.\t"
                    f"ID=transcript:{tid};Parent=gene:{gid};Name={sym}-{tid[-3:]};biotype={tbiotype}"
                    f"{tag};transcript_id={tid};version=1\n")
            order = exons if strand == "+" else list(reversed(exons))
            for rank, (s, e) in enumerate(order, 1):
                exon_n += 1
                eid = f"ENSE{exon_n:011d}"
                g.write(f"{CHROM}\ttest\texon\t{s}\t{e}\t.\t{strand}\t.\t"
                        f"Parent=transcript:{tid};Name={eid};exon_id={eid};rank={rank};version=1\n")
            if cds:
                pid = "ENSP" + tid[4:]
                done = 0
                for s, e in (cds if strand == "+" else list(reversed(cds))):
                    g.write(f"{CHROM}\ttest\tCDS\t{s}\t{e}\t.\t{strand}\t{(3 - done % 3) % 3}\t"
                            f"ID=CDS:{pid};Parent=transcript:{tid};protein_id={pid};version=1\n")
                    done += e - s + 1


def base(pos):
    return genome[pos - 1]


def landmarks(exons, cds, strand):
    """The places on a transcript where one end of a variant can be on another kind of sequence
    than the other: both sides of every boundary, and a few places that are far from any."""
    marks = set()
    lo, hi = exons[0][0], exons[-1][1]
    for s, e in exons:
        marks.update([s - 1, s, s + 1, e - 1, e, e + 1])
    for (a_s, a_e), (b_s, b_e) in zip(exons, exons[1:]):
        first, last = a_e + 1, b_s - 1
        marks.update([first + 1, first + 2, first + 3, (first + last) // 2, last - 3, last - 2, last - 1])
    if cds:
        cs, ce = cds[0][0], cds[-1][1]
        marks.update([cs - 1, cs, cs + 1, cs + 2, cs + 3, ce - 3, ce - 2, ce - 1, ce, ce + 1])
    marks.update([lo - 1100, lo - 1001, lo - 1000, lo - 500, lo - 1, hi + 1, hi + 500, hi + 1000, hi + 1001, hi + 1100])
    return sorted(m for m in marks if 3 <= m <= LENGTH - 2)


records = []  # (pos, id, ref, alt, info)
seen = set()
names = set()


def add(pos, name, alt, info):
    key = (pos, alt, info)
    if key in seen or pos < 1 or pos > LENGTH:
        return
    seen.add(key)
    # a name is one record's: the same interval written another way is another record
    unique, n = name, 1
    while unique in names:
        n += 1
        unique = f"{name}_{n}"
    names.add(unique)
    records.append((pos, unique, base(pos), alt, info))


# The classes of an interval variant, with how each is written. Every variant has the first form;
# every other one is written the second way, so that both the type from the ALT and the type from
# SVTYPE are in the world.
INTERVAL = [
    ("del", "<DEL>", "DEL"),
    ("dup", "<DUP>", "DUP"),
    ("inv", "<INV>", "INV"),
    ("cnv", "<CNV>", "CNV"),
]
RARE = [
    ("tdup", "<DUP:TANDEM>", "DUP"),
    ("alu_ins", "<INS:ME:ALU>", "INS"),
    ("alu_del", "<DEL:ME:ALU>", "DEL"),
    ("me_ins", "<INS:ME>", "INS"),
    ("cn0", "<CN0>", "CNV"),
    ("cn2", "<CN2>", "CNV"),
    ("cn3", "<CN3>", "CNV"),
]

for gid, sym, strand, biotype, trs in GENES:
    for tid, canon, exons, cds, tbiotype in trs:
        marks = landmarks(exons, cds, strand)
        for i, a in enumerate(marks):
            for j in range(i, len(marks)):
                b = marks[j]
                near = j - i <= 2
                sampled = rng.random() < 0.012
                if not (near or sampled):
                    continue
                # the variant covers `a`..`b` (a >= 3: POS is the base before)
                for k, (kind, alt, svtype) in enumerate(INTERVAL):
                    if not near and rng.random() > 0.5:
                        continue
                    form = (i + j + k) % 3
                    info = f"END={b}" if form == 0 else f"END={b};SVTYPE={svtype}" if form == 1 else f"SVTYPE={svtype};SVLEN={'-' if kind == 'del' else ''}{b - a + 1}"
                    add(a - 1, f"sv_{kind}_{a}_{b}", alt, info)
                if sampled:
                    kind, alt, svtype = RARE[(i * 7 + j) % len(RARE)]
                    if kind in ("alu_ins", "me_ins"):
                        add(a - 1, f"sv_{kind}_{a}_{b}", alt, f"SVLEN={b - a + 1}")
                    else:
                        add(a - 1, f"sv_{kind}_{a}_{b}", alt, f"END={b}")
        # the intervals whose answer turns on more than where their ends are: a whole exon with the
        # intron on both sides of it, an exon and the first bases of the next, a window across each end
        # of the coding sequence (a few bases of the UTR, a few of the CDS), and the whole transcript
        pairs = set()
        for s, e in exons:
            pairs.update([(s - 4, e + 4), (s - 1, e + 1), (s, e), (s - 3, e), (s, e + 3)])
        for (s1, e1), (s2, e2) in zip(exons, exons[1:]):
            pairs.update([(e1 - 5, s2 + 1), (s1 - 4, s2 + 1), (s1, e2)])
        if cds:
            cs, ce = cds[0][0], cds[-1][1]
            pairs.update([(cs - 4, cs + 5), (ce - 5, ce + 4), (cs - 2, cs + 2), (ce - 2, ce + 2)])
        lo, hi = exons[0][0], exons[-1][1]
        whole = [(lo - 1, hi + 1), (lo, hi), (lo - 5, hi + 5), (lo - 1, hi), (lo, hi + 1), (lo + 5, hi - 5)]
        pairs.update(whole)
        for n, (a, b) in enumerate(sorted(pairs)):
            if a < 3 or b > LENGTH - 2 or a > b:
                continue
            for k, (kind, alt, svtype) in enumerate(INTERVAL):
                form = (n + k) % 3
                info = f"END={b}" if form == 0 else f"END={b};SVTYPE={svtype}" if form == 1 else f"SVTYPE={svtype};SVLEN={'-' if kind == 'del' else ''}{b - a + 1}"
                add(a - 1, f"sv_{kind}_{a}_{b}", alt, info)
            if (a, b) in whole:
                for kind, alt, svtype in RARE:
                    if kind in ("alu_ins", "me_ins"):
                        add(a - 1, f"sv_{kind}_{a}_{b}", alt, f"SVLEN={b - a + 1}")
                    else:
                        add(a - 1, f"sv_{kind}_{a}_{b}", alt, f"END={b}")
        # an insertion at each place, of three lengths, and one with END = POS (between two bases)
        for a in marks:
            add(a, f"sv_ins_between_{a}", "<INS>", f"END={a}")
            for n, label in ((1, "one"), (40, "forty"), (400, "four_hundred")):
                add(a - 1, f"sv_ins_{label}_{a}", "<INS>", f"SVLEN={n}")

# how a record is read
first_exon = GENES[0][4][0][2][0]
ex1, ex2 = GENES[0][4][0][2][0], GENES[0][4][0][2][1]
p = ex1[0] + 60
add(p, "read_svlen_over_end", "<DEL>", f"END={p + 900};SVLEN=-25")
add(p, "read_svlen_only", "<DEL>", "SVLEN=-25")
add(p, "read_end_only", "<DEL>", f"END={p + 25}")
add(p, "read_no_extent_del", "<DEL>", ".")
add(p, "read_no_extent_dup", "<DUP>", ".")
add(p, "read_no_extent_ins", "<INS>", ".")
add(p, "read_no_extent_inv", "<INV>", ".")
add(p, "read_end_before_pos", "<INV>", f"END={p - 30}")
add(p, "read_svtype_only", "<X>", f"SVTYPE=DEL;END={p + 25}")
add(p, "read_alt_over_svtype", "<DEL>", f"SVTYPE=DUP;END={p + 25}")
add(p, "read_unknown_cpx", "<CPX>", f"SVTYPE=CPX;END={p + 25}")
add(p, "read_unknown_ctx", "<CTX>", f"SVTYPE=CTX;END={p + 25}")
add(p, "read_unknown_str", "<STR>", f"END={p + 25}")
add(p, "read_nonref", "<NON_REF>", f"END={p + 25}")
add(p, "read_cn1", "<CN1>", f"END={p + 25}")
add(p, "read_cn9", "<CN9>", f"END={p + 25}")
add(p, "read_dup_tandem_noinfo", "<DUP:TANDEM>", f"END={p + 25}")
add(p, "read_del_me", "<DEL:ME>", f"END={p + 25}")
add(p, "read_ins_herv", "<INS:ME:HERV>", "SVLEN=300")
add(p, "read_ins_line1", "<INS:ME:L1>", "SVLEN=300")
add(p, "read_ins_sva", "<INS:ME:SVA>", "SVLEN=300")
add(p, "read_del_line1", "<DEL:ME:LINE1>", f"END={p + 25}")
add(p, "read_del_sva", "<DEL:ME:SVA>", f"END={p + 25}")
add(p, "read_del_herv", "<DEL:ME:HERV>", f"END={p + 25}")
add(p, "read_ins_unknown_me", "<INS:ME:FOO>", "SVLEN=300")
# an SVLEN that is not a length, one for each of several ALTs, and one on a breakend
add(p, "read_svlen_dot_del", "<DEL>", f"END={p + 25};SVLEN=.")
add(p, "read_svlen_dot_inv", "<INV>", f"END={p + 25};SVLEN=.")
add(p, "read_svlen_zero_dup", "<DUP>", f"END={p + 25};SVLEN=0")
add(p, "read_svlen_two_values", "<DEL>,<DEL>", f"SVLEN=-30,-60")
add(p, "read_svlen_on_bnd", f"{base(p)}[{CHROM}:{p + 600}[", "SVTYPE=BND;SVLEN=40")
add(p, "read_svlen_on_bnd_huge", f"{base(p)}[{CHROM}:{p + 600}[", "SVTYPE=BND;SVLEN=20000000")
add(p, "read_svlen_huge_del", "<DEL>", "SVLEN=-20000000")
add(p, "read_end_huge_dup", "<DUP>", "END=20000000")
add(p, "read_svlen_negative_ins", "<INS>", "SVLEN=-40")
add(p, "read_single_bnd_no_svtype_after", f"{base(p)}.", ".")
add(p, "read_single_bnd_no_svtype_before", f".{base(p)}", ".")
add(p, "read_bnd_no_svtype", f"{base(p)}]{CHROM}:{p + 600}]", ".")
# several ALTs
add(p, "multi_del_dup", "<DEL>,<DUP>", f"END={p + 25}")
add(p, "multi_del_del", "<DEL>,<DEL>", f"END={p + 25}")
add(p, "multi_snv_del", "G,<DEL>", f"END={p + 25};SVTYPE=DEL")
add(p, "multi_dup_ins", "<DUP>,<INS>", f"END={p + 25}")
# breakends on this contig: one with a mate and the mate's own record, one without, and Manta's CHR2
for n, (a, b) in enumerate([
    (1050, 1320), (1400, 1650), (1650, 1700), (1160, 1250), (2100, 2500),
    (3100, 3550), (3320, 3320), (3900, 5100), (5100, 5300), (6100, 7100),
    (7250, 10050), (8010, 8210), (9050, 9250), (1500, 1750),
]):
    ra, rb = base(a), base(b)
    add(a, f"bnd_{n}_a", f"{ra}[{CHROM}:{b}[", f"SVTYPE=BND;MATEID=bnd_{n}_b")
    add(b, f"bnd_{n}_b", f"]{CHROM}:{a}]{rb}", f"SVTYPE=BND;MATEID=bnd_{n}_a")
    add(a, f"bnd_{n}_c", f"{ra}]{CHROM}:{b}]", "SVTYPE=BND")
    add(a, f"bnd_{n}_d", f"[{CHROM}:{b}[{ra}", "SVTYPE=BND")
    add(a, f"bnd_{n}_single_after", f"{ra}.", "SVTYPE=BND")
    add(b, f"bnd_{n}_single_before", f".{rb}", "SVTYPE=BND")
    add(a, f"bnd_{n}_symbolic", "<BND>", f"SVTYPE=BND;CHR2={CHROM};END2={b}")
    add(a, f"bnd_{n}_manta", "<BND>", f"SVTYPE=BND;CHR2={CHROM};END={b}")
add(1050, "bnd_two_mates", f"{base(1050)}[{CHROM}:1320[,{base(1050)}]{CHROM}:3100]", "SVTYPE=BND")
add(1050, "bnd_symbolic_alone", "<BND>", "SVTYPE=BND")
add(1050, "bnd_two_mates_one_far", "%s[%s:1320[,%s]%s:12500]" % (base(1050), CHROM, base(1050), CHROM), "SVTYPE=BND")
# a tandem repeat that VEP reads as bases: CA x 20 at 3400..3439, in an intron of SVB
unit, at, copies = REPEAT_UNIT, REPEAT_AT, REPEAT_COPIES
span = 2 * copies
for n, (more, label) in enumerate([(5, "grows"), (-4, "shrinks"), (0, "same"), (1, "one_more")]):
    count = copies + more
    add(at - 1, f"tr_{label}", "<CNV:TR>",
        f"SVTYPE=CNV;END={at + span - 1};SVLEN={span};RN=1;RUS={unit};RUL=2;RUC={count};RB={2 * count}")
add(at - 1, "tr_rb_only", "<CNV:TR>", f"SVTYPE=CNV;END={at + span - 1};SVLEN={span};RN=1;RUS={unit};RUL=2;RB={2 * (copies + 3)}")
add(at - 1, "tr_two_alleles", "<CNV:TR>,<CNV:TR>",
    f"SVTYPE=CNV;END={at + span - 1};SVLEN={span},{span};RN=1,1;RUS={unit},{unit};RUC={copies + 2},{copies - 2}")
add(at - 1, "tr_no_info", "<CNV:TR>", f"SVTYPE=CNV;END={at + span - 1};SVLEN={span}")
add(at - 1, "tr_tandem_dup", "<DUP:TANDEM>", f"END={at + span - 1};SVLEN={span};RUS={unit};RUC={copies + 1}")
a1 = at - 1
tr_base = f"SVTYPE=CNV;END={at + span - 1};SVLEN={span}"
add(a1, "tr_zero", "<CNV:TR>", f"{tr_base};RN=1;RUS={unit};RUC=0")
add(a1, "tr_missing_unit", "<CNV:TR>", f"{tr_base};RN=1;RUS=.;RUC=7")
add(a1, "tr_rb_not_whole", "<CNV:TR>", f"{tr_base};RN=1;RUS={unit};RB=45")
add(a1, "tr_two_units", "<CNV:TR>", f"{tr_base};RN=2;RUS=CA,GT;RUC=3,2")
add(a1, "tr_two_alleles_two_units", "<CNV:TR>,<CNV:TR>", f"{tr_base},{span};RN=2,1;RUS=CA,GT,AC;RUC=3,2,4")
add(a1, "tr_two_alleles_no_rn", "<CNV:TR>,<CNV:TR>", f"{tr_base},{span};RUS=CA,GT;RUC=9,11")
add(a1, "tr_rn_dot", "<CNV:TR>,<CNV:TR>", f"{tr_base},{span};RN=.,1;RUS=CA;RUC=9")
add(a1, "tr_svlen_not_end", "<CNV:TR>", f"SVTYPE=CNV;END={at + span - 1};SVLEN=30;RN=1;RUS={unit};RUC=21")
add(a1, "tr_no_extent", "<CNV:TR>", f"SVTYPE=CNV;RN=1;RUS={unit};RUC=21")
add(a1, "tr_oversized", "<CNV:TR>", f"{tr_base};RN=1;RUS={unit};RUC=6000000")
add(a1, "tr_few_counts", "<CNV:TR>", f"{tr_base};RN=2;RUS=CA,GT;RUC=3")
# the change `tr_one_more` makes, written as the plain record it becomes
anchor, repeat = base(at - 1), genome[at - 1:at - 1 + span]
records.append((at - 1, "plain_one_more", anchor + repeat, anchor + repeat + unit, "."))
# a long deletion written as bases, a tandem duplication of the exon of SVD, and records at both ends of the contig
long_ref = genome[5599:5640]
records.append((5600, "bases_long_deletion", long_ref, long_ref[0], "."))

records.sort(key=lambda r: (r[0], r[1]))
with open(f"{OUT}/in.vcf", "w") as v:
    v.write("##fileformat=VCFv4.2\n")
    v.write(f"##contig=<ID={CHROM},length={LENGTH}>\n")
    v.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n")
    v.write("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Length\">\n")
    v.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n")
    v.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for pos, name, ref, alt, info in records:
        v.write(f"{CHROM}\t{pos}\t{name}\t{ref}\t{alt}\t.\t.\t{info}\n")
print(f"{len(records)} records")
