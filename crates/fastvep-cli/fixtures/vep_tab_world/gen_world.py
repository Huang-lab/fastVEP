#!/usr/bin/env python3
"""A small world for comparing fastVEP's tab output with Ensembl VEP's, cell by cell.

    gen_world.py OUTDIR N [SEED] [--all]

chr10 (30 kb) carries four genes on both strands, a second transcript for one of them, a
non-coding gene and an antisense one that overlaps another. Coding sequence is built from
non-stop codons with a start and a stop codon, and every transcript has a UTR long enough
for a frameshift to end inside it, so that a difference between the two tools is a
difference of rule and not of a broken protein.

The variants are the ones fastVEP is expected to describe exactly as VEP does: SNVs, MNVs
and indels written minimally (one anchor base, nothing shared beyond it), alone or several
to a record, with and without IDs. Left out, because they are known to differ and are
listed in docs/VEP_DIVERGENCE.md ("Cells of the tab output"):

  - a change that crosses a boundary of an exon or of the coding sequence. VEP writes
    `501-?` where only one end is inside, and reports both an EXON and an INTRON number for
    a change that is in both. fastVEP does the same now, and `../vep_tab_edges` is where it
    is checked, change by change; they are left out here so that this world stays the
    baseline of what the rest of the tab output does.
  - an indel in the first or last four codons of a coding sequence, where VEP calls a
    frameshift what fastVEP calls inframe_insertion plus stop_retained, and an indel
    within ten bases of a transcript's end, where VEP writes no HGVSc.
  - an indel whose alleles repeat bases at either end beyond the anchor, which VEP
    minimises before it annotates it and fastVEP does not (a delins such as CA>GGA).
  - a record with several alts in which two are indels, whose second and later alleles VEP
    places differently, and a record in which an allele trims to a substitution, for which
    VEP writes a delins of the untrimmed alt over the trimmed reference span (`g.5559delins
    GGCCA` for `TGCCA>GGCCA`, which is `g.5559T>G`).
  - an MNV that is an inversion only once its shared ends are clipped, which VEP names a
    delins and fastVEP names an inversion, as HGVS does.

With --all nothing is left out, which is how the size of each class above is measured
(validation/run_tab_world.sh).

Writes OUTDIR/{genome.fa, genes.gff3, in.vcf}.
"""
import random
import sys

ALL = "--all" in sys.argv
args = [a for a in sys.argv[1:] if a != "--all"]
OUT = args[0]
N = int(args[1])
SEED = int(args[2]) if len(args) > 2 else 139
rng = random.Random(SEED)

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
STOPS = {"TAA", "TAG", "TGA"}
SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS | {"ATG"}]


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


LENGTH = 30000
seq = [rng.choice("ACGT") for _ in range(LENGTH)]

# (gene id, symbol, strand, biotype, transcripts)
# a transcript: (id, canonical, exons [(s, e)] ascending, cds [(s, e)] ascending or [])
GENES = [
    ("ENSG00000000001", "GENEA", "+", "protein_coding", [
        ("ENST00000000001", True, [(700, 1300), (2000, 2200), (3000, 3600)],
         [(1100, 1300), (2000, 2200), (3000, 3251)]),
        ("ENST00000000002", False, [(700, 1300), (3000, 3600)],
         [(1100, 1300), (3000, 3251)]),
    ]),
    ("ENSG00000000002", "GENEB", "-", "protein_coding", [
        ("ENST00000000003", True, [(9700, 10400), (11000, 11700)],
         [(10099, 10400), (11000, 11300)]),
    ]),
    ("ENSG00000000003", "GENEC", "+", "lncRNA", [
        ("ENST00000000004", True, [(20000, 20500), (21000, 21500)], []),
    ]),
    ("ENSG00000000004", "GENEE", "-", "lncRNA", [
        ("ENST00000000005", True, [(2500, 2600), (3700, 3800)], []),
    ]),
]


def place(intervals, strand, tseq):
    plus = tseq if strand == "+" else revcomp(tseq)
    i = 0
    for s, e in intervals:
        for pos in range(s, e + 1):
            seq[pos - 1] = plus[i]
            i += 1
    assert i == len(plus)


for gid, sym, strand, biotype, trs in GENES:
    cds = trs[0][3]
    if not cds:
        continue
    n = sum(e - s + 1 for s, e in cds)
    assert n % 3 == 0, (sym, n)
    codons = ["ATG"] + [rng.choice(SENSE) for _ in range(n // 3 - 2)] + [rng.choice(sorted(STOPS))]
    place(cds, strand, "".join(codons))

genome = "".join(seq)
with open(f"{OUT}/genome.fa", "w") as f:
    f.write(">chr10\n")
    for i in range(0, LENGTH, 60):
        f.write(genome[i:i + 60] + "\n")

with open(f"{OUT}/genes.gff3", "w") as g:
    g.write("##gff-version 3\n")
    g.write(f"##sequence-region   chr10 1 {LENGTH}\n")
    exon_n = 0
    for gid, sym, strand, biotype, trs in GENES:
        lo = min(e[0] for t in trs for e in t[2])
        hi = max(e[1] for t in trs for e in t[2])
        g.write(f"chr10\ttest\tgene\t{lo}\t{hi}\t.\t{strand}\t.\t"
                f"ID=gene:{gid};Name={sym};biotype={biotype};gene_id={gid};version=1\n")
        for tid, canon, exons, cds in trs:
            ttype = "mRNA" if cds else "lnc_RNA"
            tag = ";tag=Ensembl_canonical" if canon else ""
            g.write(f"chr10\ttest\t{ttype}\t{exons[0][0]}\t{exons[-1][1]}\t.\t{strand}\t.\t"
                    f"ID=transcript:{tid};Parent=gene:{gid};Name={sym}-{tid[-3:]};biotype={biotype}"
                    f"{tag};transcript_id={tid};version=1\n")
            order = exons if strand == "+" else list(reversed(exons))
            for rank, (s, e) in enumerate(order, 1):
                exon_n += 1
                eid = f"ENSE{exon_n:011d}"
                g.write(f"chr10\ttest\texon\t{s}\t{e}\t.\t{strand}\t.\t"
                        f"Parent=transcript:{tid};Name={eid};exon_id={eid};rank={rank};version=1\n")
            if cds:
                pid = "ENSP" + tid[4:]
                done = 0
                for s, e in (cds if strand == "+" else list(reversed(cds))):
                    g.write(f"chr10\ttest\tCDS\t{s}\t{e}\t.\t{strand}\t{(3 - done % 3) % 3}\t"
                            f"ID=CDS:{pid};Parent=transcript:{tid};protein_id={pid};version=1\n")
                    done += e - s + 1

# the places where one end of a change can be on a different kind of sequence than the other
edges = set()
for gid, sym, strand, biotype, trs in GENES:
    for tid, canon, exons, cds in trs:
        for s, e in exons:
            edges.update([s, e + 1])
        for s, e in cds:
            edges.update([s, e + 1])
edges = sorted(edges)

# the first and last codons of the coding span: a change there has rules of its own (start_lost,
# stop_lost, stop_retained, a frameshift read against the stop), and an indel in them is one of the
# cases the two tools do not yet agree on
codon_zones = []
for gid, sym, strand, biotype, trs in GENES:
    for tid, canon, exons, cds in trs:
        if cds:
            # the first and last four codons of the coding span, whichever way the gene runs
            codon_zones.append((cds[0][0], cds[0][0] + 11))
            codon_zones.append((cds[-1][1] - 11, cds[-1][1]))


transcript_ends = []
for gid, sym, strand, biotype, trs in GENES:
    for tid, canon, exons, cds in trs:
        transcript_ends.append((exons[0][0], exons[0][0] + 9))
        transcript_ends.append((exons[-1][1] - 9, exons[-1][1]))


def near_transcript_end(lo, hi):
    return any(lo <= b and a <= hi for a, b in transcript_ends)


def in_codon_zone(lo, hi):
    return any(lo <= b and a <= hi for a, b in codon_zones)


def crosses(lo, hi):
    """Does the span lo..hi (inclusive, for an insertion the two bases it sits between) have
    an edge strictly inside it?  An edge e separates e-1 from e."""
    return any(lo < e <= hi for e in edges)


def base(pos):
    return genome[pos - 1]


def other(b):
    return rng.choice([x for x in "ACGT" if x != b])


def rand_seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def make_alt(pos):
    """One alt at pos as (ref, alt) in VCF form, or None when it is not one to keep."""
    kind = rng.random()
    if kind < 0.50:  # SNV
        r = base(pos)
        return r, other(r)
    if kind < 0.65:  # deletion, anchored on the base before
        n = rng.randint(1, 8)
        ref = genome[pos - 1:pos + n]
        return ref, ref[0]
    if kind < 0.83:  # insertion after pos
        n = rng.randint(1, 8)
        r = base(pos)
        ins = rand_seq(n)
        if rng.random() < 0.35:  # copy what is next to it, so that some are duplications
            ins = genome[pos:pos + n] if rng.random() < 0.5 else genome[pos - n:pos]
        return r, r + ins
    # MNV: the same length, differing somewhere
    n = rng.randint(2, 5)
    ref = genome[pos - 1:pos - 1 + n]
    alt = "".join(c if rng.random() < 0.4 else other(c) for c in ref)
    if alt == ref:
        alt = other(ref[0]) + ref[1:]
    return ref, alt


def span(pos, ref, alts):
    """The bases the record touches, with the neighbours an insertion sits between."""
    lo = pos
    hi = pos + len(ref) - 1
    if any(len(a) > len(ref) for a in alts) or any(len(a) < len(ref) for a in alts):
        lo, hi = pos - 1, hi + 1
    return lo, hi


def positions():
    while True:
        r = rng.random()
        if r < 0.40:
            pos = rng.choice(edges) + rng.randint(-14, 14)
        elif r < 0.65:
            gid, sym, strand, biotype, trs = rng.choice(GENES)
            e = [x for t in trs for x in t[2]]
            pos = rng.randint(min(a for a, b in e), max(b for a, b in e))
        elif r < 0.85:
            gid, sym, strand, biotype, trs = rng.choice(GENES)
            e = [x for t in trs for x in t[2]]
            pos = rng.choice([min(a for a, b in e), max(b for a, b in e)]) + rng.randint(-4800, 4800)
        else:
            pos = rng.randint(5, LENGTH - 40)
        if 5 <= pos <= LENGTH - 40:
            return pos


records = []
seen = set()
guard = 0
while len(records) < N and guard < N * 200:
    guard += 1
    pos = positions()
    multi = rng.random() < 0.12
    first = make_alt(pos)
    ref, alts = first[0], [first[1]]
    if multi:
        for _ in range(rng.randint(1, 2)):
            r2, a2 = make_alt(pos)
            if len(r2) > len(ref):
                ext = r2[len(ref):]
                alts = [a + ext for a in alts]
                ref = r2
            elif len(r2) < len(ref):
                a2 += ref[len(r2):]
            if a2 != ref and a2 not in alts:
                alts.append(a2)
    if pos + len(ref) > LENGTH - 2:
        continue
    lo, hi = span(pos, ref, alts)
    if crosses(lo, hi) and not ALL:  # an MNV can reach across a boundary as well as an indel
        continue
    if any(len(a) != len(ref) for a in alts):
        if in_codon_zone(lo, hi) and not ALL:
            continue
        # every indel stays minimal: the alleles share nothing past the first base
        if not ALL and any(len(a) != len(ref) and (a[-1] == ref[-1] and len(a) > 1 and len(ref) > 1) for a in alts):
            continue
    if any(len(a) != len(ref) for a in alts) and near_transcript_end(lo, hi) and not ALL:
        continue
    if len(alts) > 1 and not ALL:
        # several alleles are one variant in VEP 116, annotated whole. Two things it does
        # with them are not yet matched here, and are listed in docs/VEP_DIVERGENCE.md: it
        # places the second and later insertions of one record differently, and writes
        # a delins of the untrimmed alt over the trimmed reference span for an allele
        # that trims to a substitution.
        if sum(1 for a in alts if len(a) != len(ref)) > 1:
            continue
        def trims_to_a_substitution(r, a):
            t_r, t_a = r, a
            while t_r and t_a and t_r[0] == t_a[0]:
                t_r, t_a = t_r[1:], t_a[1:]
            while t_r and t_a and t_r[-1] == t_a[-1]:
                t_r, t_a = t_r[:-1], t_a[:-1]
            return (t_r, t_a) != (r, a) and t_r and t_a
        if any(trims_to_a_substitution(ref, a) for a in alts):
            continue
    # an MNV that is an inversion only once its shared ends are clipped: VEP says delins
    # for the HGVSc and fastVEP says inv, which is the HGVS name
    if not ALL and len(alts) == 1 and len(alts[0]) == len(ref) and len(ref) > 2:
        t_r, t_a = ref, alts[0]
        while t_r and t_a and t_r[0] == t_a[0]:
            t_r, t_a = t_r[1:], t_a[1:]
        while t_r and t_a and t_r[-1] == t_a[-1]:
            t_r, t_a = t_r[:-1], t_a[:-1]
        if (t_r, t_a) != (ref, alts[0]) and len(t_r) >= 2 and t_a == revcomp(t_r):
            continue
    if not ALL and len(alts) > 1 and len({len(a) for a in alts} | {len(ref)}) > 1:
        # mixed lengths in one record: only when all alleles keep the same first base
        if len({x[0] for x in [ref] + alts}) != 1:
            continue
    key = (pos, ref, tuple(alts))
    if key in seen:
        continue
    seen.add(key)
    r = rng.random()
    vid = "." if r < 0.7 else (f"rs{rng.randint(1, 999999)}" if r < 0.9 else
                               f"rs{rng.randint(1, 999999)};rs{rng.randint(1, 999999)}")
    records.append((pos, vid, ref, alts))

records.sort(key=lambda r: (r[0], r[2], r[3]))
with open(f"{OUT}/in.vcf", "w") as v:
    v.write(f"##fileformat=VCFv4.2\n##contig=<ID=chr10,length={LENGTH}>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for pos, vid, ref, alts in records:
        v.write(f"chr10\t{pos}\t{vid}\t{ref}\t{','.join(alts)}\t.\t.\t.\n")
print(f"{len(records)} records")
