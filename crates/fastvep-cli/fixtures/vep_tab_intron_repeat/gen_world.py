#!/usr/bin/env python3
"""A world of transcripts whose coding repeats continue across an intron, and a twin of it with the
intron cut out, for comparing fastVEP's tab output with Ensembl VEP's (issue #139).

    gen_world.py OUTDIR [SAMPLE] [SEED]

HGVS describes a change in the sequence of the transcript, and says to put it at its 3'-most
position there. Ensembl moves it along the *genome* (`_genomic_shift`) and describes it where that
ends. The two are the same until a repeat goes on across an intron: then the genome carries the
change over bases the transcript does not have.

Every transcript here has two coding exons and an intron, and the genome repeats a short unit
(`A`, `AG`, `GCT`, ...) through the end of the first exon, all of the intron and the start of the
second, so a change inside the repeat that is moved 3' reaches the intron and goes over it. There is one case for each
unit, each intron length (1 to 14 bases, around Ensembl's frameshift-intron limit of 13) and
each strand, and the coding frame is shifted from case to case so that the intron falls in every
codon phase.

Each is built twice. Contig `ra` has the transcript as above; contig `rb` has its twin, the same
spliced transcript in one exon with the intron cut out of the genome. A change inside an exon has
the same place in the transcript in both and so, in an annotator that works in the transcript,
the same description: `rb` is the oracle for `ra`, made by the tool being compared.

The changes are, at every position of the repeat that is in an exon (so that both bases an
insertion sits between, and every base a deletion removes, are in the same exon), insertions of the
unit as it reads there, of one base and of a few random ones, and deletions of one base and of a
unit. A change on `ra` and the same change on `rb` have the IDs `i_N` and `t_N`.

Writes OUTDIR/{genome.fa, genes.gff3, in.vcf}.
"""
import random
import sys

OUT = sys.argv[1]
SAMPLE = int(sys.argv[2]) if len(sys.argv) > 2 else 3
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 139
rng = random.Random(SEED)

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
STOPS = {"TAA", "TAG", "TGA"}
SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS | {"ATG"}]


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


def has_stop(unit):
    """A stop codon in any frame of the unit repeated, which would end the protein inside the repeat."""
    s = unit * 12
    return any(s[i:i + 3] in STOPS for f in range(3) for i in range(f, len(s) - 2, 3))


CANDIDATES = ["A", "C", "AG", "CT", "GCT", "AAG", "GCAT", "AAGC", "GGAAC", "CAGCAT"]
UNITS = [u for u in CANDIDATES if not has_stop(u)]
assert UNITS == CANDIDATES, "a unit reads as a stop codon in some frame"
# 12 and 14 straddle Ensembl's frameshift-intron limit (13 bases or fewer)
INTRON_LENGTHS = [1, 2, 3, 4, 6, 12, 14]
PAD = 60
UTR = 20

contigs = {"ra": [], "rb": []}  # pieces of each contig, in order
size = {"ra": 0, "rb": 0}
gff = []  # (contig, line)
records = []  # (n, pos_a, pos_b, ref, alt_a)


def add(contig, seq):
    """Append bases to a contig and give back the 1-based coordinate of the first."""
    start = size[contig] + 1
    contigs[contig].append(seq)
    size[contig] += len(seq)
    return start


def gff_gene(contig, gid, name, strand, exons, cds):
    """exons and cds are ascending (start, end) lists."""
    tid, pid = "ENST" + gid[4:], "ENSP" + gid[4:]
    lo, hi = exons[0][0], exons[-1][1]
    gff.append((contig, f"{contig}\ttest\tgene\t{lo}\t{hi}\t.\t{strand}\t.\t"
                        f"ID=gene:{gid};Name={name};biotype=protein_coding;gene_id={gid};version=1"))
    gff.append((contig, f"{contig}\ttest\tmRNA\t{lo}\t{hi}\t.\t{strand}\t.\t"
                        f"ID=transcript:{tid};Parent=gene:{gid};Name={name}-001;biotype=protein_coding;"
                        f"tag=Ensembl_canonical;transcript_id={tid};version=1"))
    order = exons if strand == "+" else list(reversed(exons))
    for rank, (s, e) in enumerate(order, 1):
        eid = f"ENSE{gid[4:]}{rank}"
        gff.append((contig, f"{contig}\ttest\texon\t{s}\t{e}\t.\t{strand}\t.\t"
                            f"Parent=transcript:{tid};Name={eid};exon_id={eid};rank={rank};version=1"))
    done = 0
    for s, e in (cds if strand == "+" else list(reversed(cds))):
        gff.append((contig, f"{contig}\ttest\tCDS\t{s}\t{e}\t.\t{strand}\t{(3 - done % 3) % 3}\t"
                            f"ID=CDS:{pid};Parent=transcript:{tid};protein_id={pid};version=1"))
        done += e - s + 1


def span(start, length, strand, lo, hi):
    """Offsets [lo, hi) of a transcript sequence of `length` bases, placed at `start`, as ascending
    genome coordinates (first, last)."""
    if strand == "+":
        return start + lo, start + hi - 1
    return start + length - hi, start + length - 1 - lo


def random_bases(k):
    return "".join(rng.choice("ACGT") for _ in range(k))


def build_transcript(unit, n, m, intron_len):
    """The two exons and the intron of one case, in transcript order, with a protein that has no
    stop but its last codon."""
    for _ in range(1000):
        run = (unit * ((n + intron_len + m) // len(unit) + 2))[:n + intron_len + m]
        pre = "ATG" + "".join(rng.choice(SENSE) for _ in range(5))
        # the spliced repeat is n + m bases, and the protein needs whole codons
        fill = random_bases((-(n + m)) % 3)
        post = "".join(rng.choice(SENSE) for _ in range(5)) + rng.choice(sorted(STOPS))
        exon1 = random_bases(UTR) + pre + run[:n]
        exon2 = run[n + intron_len:] + fill + post + random_bases(UTR)
        cds = (exon1 + exon2)[UTR:-UTR]
        if len(cds) % 3 == 0 and not any(cds[i:i + 3] in STOPS for i in range(0, len(cds) - 3, 3)):
            return exon1, run[n:n + intron_len], exon2
    raise SystemExit(f"no protein without a stop for {unit}, {intron_len}")


# Every transcript is 1,500 bases from the ends of its contig. In an earlier layout of this world
# the last transcript was 60 bases from the end, and three of its changes had an HGVSc from Ensembl
# that is not the change (`c.20dup` for an insertion of a base that the base before it is not),
# which the same world with room around every transcript did not show. Ensembl shifts a change
# along 1,000 bases of the chromosome either side of it; the cause was not looked into, because it
# is not what this world is for.
for contig in ("ra", "rb"):
    add(contig, random_bases(1500))

case = 0
for si, strand in enumerate("+-"):
    for ui, unit in enumerate(UNITS):
        for li, intron_len in enumerate(INTRON_LENGTHS):
            case += 1
            u = len(unit)
            # bases of the repeat before the intron, in the first exon, and after it, in the second;
            # the first varies from case to case, which moves the intron through the codon phases
            n, m = 2 * u + 3 + (ui + li + si) % 3, 2 * u + 3
            exon1, intron, exon2 = build_transcript(unit, n, m, intron_len)

            placed = {}
            for contig, pieces in (("ra", (exon1, intron, exon2)), ("rb", (exon1 + exon2,))):
                add(contig, random_bases(PAD))
                whole = "".join(pieces)
                start = add(contig, whole if strand == "+" else revcomp(whole))
                add(contig, random_bases(PAD))
                if contig == "ra":
                    spans = [(0, len(exon1)), (len(exon1) + len(intron), len(whole))]
                    cds_tx = [(UTR, len(exon1)), (len(exon1) + len(intron), len(whole) - UTR)]
                else:
                    spans = [(0, len(whole))]
                    cds_tx = [(UTR, len(whole) - UTR)]
                exons = sorted(span(start, len(whole), strand, a, b) for a, b in spans)
                cds = sorted(span(start, len(whole), strand, a, b) for a, b in cds_tx)
                gid = f"ENSG{'1' if contig == 'ra' else '2'}{case:07d}"
                gff_gene(contig, gid, f"IR{case}{'a' if contig == 'ra' else 'b'}", strand, exons, cds)
                placed[contig] = (start, whole, exons)

            a_start, a_whole, a_exons = placed["ra"]
            b_start = placed["rb"][0]
            genome_a = a_whole if strand == "+" else revcomp(a_whole)

            def base(p):
                return genome_a[p - a_start]

            def exon_of(p):
                return next((k for k, (s, e) in enumerate(a_exons) if s <= p <= e), None)

            def to_twin(p):
                """Where a base of an exon of ra is on rb: the exon that comes second in the genome
                (the second in transcript order on the plus strand, the first on the minus) is
                that much further along on ra, which has the intron in between."""
                return b_start + (p - a_start) - (intron_len if exon_of(p) == 1 else 0)

            # the exonic bases of the repeat: the end of the first exon and the start of the second
            ends = [span(a_start, len(a_whole), strand, len(exon1) - n, len(exon1)),
                    span(a_start, len(a_whole), strand, len(exon1) + len(intron),
                         len(exon1) + len(intron) + m)]
            positions = sorted({p for lo, hi in ends for p in range(lo, hi + 1)})
            for p in positions:
                k = exon_of(p)
                assert k is not None
                e = a_exons[k][1]
                ref = base(p)
                changes = []
                # insertions between p and p + 1: the unit as it reads there, a base and a few
                if p + u <= e:
                    changes.append((p, ref, ref + "".join(base(q) for q in range(p + 1, p + u + 1))))
                if p + 1 <= e:
                    changes.append((p, ref, ref + rng.choice("ACGT")))
                    changes.append((p, ref, ref + random_bases(rng.randint(2, 5))))
                # deletions of the base after p and of a unit
                for width in sorted({1, u}):
                    if p + width <= e:
                        changes.append((p, ref + "".join(base(q) for q in range(p + 1, p + width + 1)), ref))
                for pos, r, a in changes:
                    records.append((len(records) + 1, pos, to_twin(pos), r, a))

for contig in ("ra", "rb"):
    add(contig, random_bases(1500))

# a record is kept when it is one of every SAMPLE-th
kept = [r for r in records if r[0] % SAMPLE == 0]

ra_seq = "".join(contigs["ra"])
rb_seq = "".join(contigs["rb"])
with open(f"{OUT}/genome.fa", "w") as f:
    for name, seq in (("ra", ra_seq), ("rb", rb_seq)):
        f.write(f">{name}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + "\n")

with open(f"{OUT}/genes.gff3", "w") as g:
    g.write("##gff-version 3\n")
    for name, seq in (("ra", ra_seq), ("rb", rb_seq)):
        g.write(f"##sequence-region   {name} 1 {len(seq)}\n")
    for _, line in gff:
        g.write(line + "\n")

# VEP and fastVEP both ask for the REF to be the genome's, so the twin's bases are checked here
with open(f"{OUT}/in.vcf", "w") as o:
    o.write("##fileformat=VCFv4.2\n")
    o.write(f"##contig=<ID=ra,length={len(ra_seq)}>\n##contig=<ID=rb,length={len(rb_seq)}>\n")
    o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for n, pa, pb, ref, alt in kept:
        assert ra_seq[pa - 1:pa - 1 + len(ref)] == ref, (n, pa, ref)
        assert rb_seq[pb - 1:pb - 1 + len(ref)] == ref, (n, pb, ref)
        o.write(f"ra\t{pa}\ti_{n}\t{ref}\t{alt}\t.\t.\t.\n")
    for n, pa, pb, ref, alt in kept:
        o.write(f"rb\t{pb}\tt_{n}\t{ref}\t{alt}\t.\t.\t.\n")
print(f"{case} transcripts, {len(kept)} of {len(records)} changes, each on both contigs")
