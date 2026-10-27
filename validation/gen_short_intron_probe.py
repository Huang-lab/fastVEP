#!/usr/bin/env python3
"""Variants around every exon edge of the transcripts that have a frameshift intron.

    gen_short_intron_probe.py GFF3 FASTA OUT.vcf RADIUS TRANSCRIPT_ID...

Ensembl VEP treats a transcript with an intron of 13 bases or fewer (`end - start <= 12`) as
having a frameshift intron, which changes what it says near every exon of that transcript
(docs/VEP_DIVERGENCE.md, "Frameshift introns"). `TRANSCRIPT_ID` is the value of the exons'
Parent attribute without a `transcript:` prefix (`ENST00000673971`, `rna-NM_001372044.2`).
FASTA needs a `.fai` beside it. The transcripts of a GFF3 that have such an intron are found with:

    python3 - GFF3 <<'PY'
    import collections, sys
    ex = collections.defaultdict(list)
    for l in open(sys.argv[1]):
        f = l.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "exon":
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            for p in a.get("Parent", "").split(","):
                ex[p].append((int(f[3]), int(f[4])))
    for t, e in ex.items():
        e.sort()
        if any(s <= n and n - s <= 12 for (_, a), (b, _) in zip(e, e[1:]) for s, n in [(a + 1, b - 1)]):
            print(t)
    PY

At every position within RADIUS of an exon edge, and in every intron of 40 bases or fewer, one
record per shape, each a single ALT: a substitution, a deletion of 1 and of 3 bases and an
insertion of 1 and of 3.
"""
import sys

gff, fa, out, radius = sys.argv[1:5]
radius = int(radius)
want = set(sys.argv[5:])
fai = {}
for line in open(fa + ".fai"):
    name, length, offset, line_bases, line_width = line.split("\t")[:5]
    fai[name] = (int(length), int(offset), int(line_bases), int(line_width))
fh = open(fa, "rb")


def sequence(chrom, start, end):
    """Bases `start` to `end`, 1-based inclusive."""
    _, offset, line_bases, line_width = fai[chrom]
    first = offset + ((start - 1) // line_bases) * line_width + (start - 1) % line_bases
    last = offset + ((end - 1) // line_bases) * line_width + (end - 1) % line_bases + 1
    fh.seek(first)
    return fh.read(last - first).decode().replace("\n", "").replace("\r", "").upper()


exons, chroms = {}, {}
for line in open(gff):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "exon":
        continue
    attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    for parent in attrs.get("Parent", "").split(","):
        key = parent.split(":", 1)[1] if parent.startswith("transcript:") else parent
        if key in want:
            exons.setdefault(key, []).append((int(f[3]), int(f[4])))
            chroms[key] = f[0]
if set(exons) != want:
    sys.exit(f"no exons for {sorted(want - set(exons))}")

positions = set()
for transcript, exon_list in exons.items():
    exon_list.sort()
    chrom = chroms[transcript]
    for start, end in exon_list:
        for edge in (start, end):
            positions.update((chrom, p) for p in range(edge - radius, edge + radius + 1))
    for (_, end_a), (start_b, _) in zip(exon_list, exon_list[1:]):
        if start_b - end_a - 1 <= 40:
            positions.update((chrom, p) for p in range(end_a - 2, start_b + 3))

other = {"A": "C", "C": "G", "G": "T", "T": "A"}
rows = []
for chrom, pos in sorted(positions):
    ref = sequence(chrom, pos, pos)
    after = sequence(chrom, pos + 1, pos + 3)
    rows.append((chrom, pos, ref, other[ref]))
    rows.append((chrom, pos, ref + after[:1], ref))
    rows.append((chrom, pos, ref + after[:3], ref))
    rows.append((chrom, pos, ref, ref + "G"))
    rows.append((chrom, pos, ref, ref + "ACG"))

with open(out, "w") as o:
    o.write("##fileformat=VCFv4.2\n")
    for chrom in sorted({r[0] for r in rows}):
        o.write(f"##contig=<ID={chrom}>\n")
    o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for chrom, pos, ref, alt in rows:
        o.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")
print(f"{len(rows)} records at {len(positions)} positions around {len(exons)} transcripts")
