#!/usr/bin/env python3
"""Which transcript Ensembl VEP and fastVEP flag with `--flag_pick*`, and why they differ where they do.

    compare_pick.py --by variant|allele|allele-gene --gff3 GENES.gff3 VEP.tab FASTVEP.tab

Both tables are tab outputs of the same VCF over the same gene models with one flagging form
(`--flag_pick`, `--flag_pick_allele` or `--flag_pick_allele_gene`, which is what `--by` names: the group
that gets one flag is the variant, the variant and allele, or the variant, allele and gene). For each group
that has a transcript row in both, the transcripts whose `PICK` is 1 are compared; a group that
differs is put under the first of these that holds, with the GFF3 read the way Ensembl VEP's own reader
reads it (`AnnotationSource/File/BaseGXF.pm`: `transcript_support_level`, the `MANE_Select`,
`MANE_Plus_Clinical` and `Ensembl_canonical` tags, the biotype; no CCDS, no APPRIS):

  not a row of VEP's    a transcript fastVEP flags that VEP wrote no row for (its type is named: VEP's
                        reader skips `unconfirmed_transcript` and others)
  a criterion before    MANE, canonical, TSL or biotype differs between the two flagged transcripts,
  the rank              which should have decided the pick in both tools
  a better term         the consequences of the two differ: VEP's rank criterion decided it, and
                        the transcript fastVEP flagged had a CCDS or an APPRIS that VEP's reader never
                        sets (so the criterion that fastVEP applied before the rank was not VEP's)
  the longer one        the consequences are the same and VEP's transcript is the longer, which is
                        the ninth criterion of `--pick-order` (`length`): the length of the translateable
                        sequence, with the bases of the start phase, for a transcript that has a
                        translation, and the sum of the lengths of its exons for one that has not
  other                 anything else

It also counts the variants that no transcript is beside (every row of them has no feature): how many each tool has, and
how many of their rows each flags, which is what `docs/VEP_DIVERGENCE.md` divergence 10 quotes.

The counts are what `docs/VEP_DIVERGENCE.md`, Part 2, quotes. Standard library only.
"""

from __future__ import annotations

import argparse
import collections
import sys


def read_table(path: str) -> list[dict[str, str]]:
    header, rows = None, []
    with open(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            cells = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                header = [c.lstrip("#") for c in cells]
            else:
                rows.append(dict(zip(header, cells)))
    if header is None or "PICK" not in header:
        sys.exit(f"{path}: no PICK column (is it the output of a flagging form?)")
    return rows


class Gff3:
    """What `compare_pick.py` needs of a GFF3: for each transcript its type, tags, biotype, TSL, CCDS, and its lengths."""

    def __init__(self, path: str, contigs: set[str]):
        self.kind, self.mane, self.canonical, self.tsl, self.biotype, self.ccds = {}, {}, {}, {}, {}, {}
        exons: dict[str, int] = collections.defaultdict(int)
        cds: dict[str, list[tuple[int, int, str, str]]] = collections.defaultdict(list)
        with open(path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9 or cols[0] not in contigs:
                    continue
                attrs = dict(kv.split("=", 1) for kv in cols[8].split(";") if "=" in kv)
                parent = attrs.get("Parent", "").removeprefix("transcript:")
                if cols[2] == "exon":
                    exons[parent] += int(cols[4]) - int(cols[3]) + 1
                elif cols[2] == "CDS":
                    cds[parent].append((int(cols[3]), int(cols[4]), cols[6], cols[7]))
                elif attrs.get("ID", "").startswith("transcript:"):
                    t = attrs["ID"].removeprefix("transcript:")
                    tags = attrs.get("tag", "").split(",")
                    self.kind[t] = cols[2]
                    self.mane[t] = "MANE_Select" in tags or "MANE_Plus_Clinical" in tags
                    self.canonical[t] = "Ensembl_canonical" in tags
                    tsl = attrs.get("transcript_support_level", "").split(" ")[0]
                    self.tsl[t] = int(tsl) if tsl.isdigit() else 100
                    self.biotype[t] = attrs.get("biotype", "")
                    self.ccds[t] = "ccdsid" in attrs
        self.exons, self.cds = exons, cds

    def length(self, t: str) -> int:
        """VEP's `length` criterion for a transcript: how long it is, as `pick_worst_VariationFeatureOverlapAllele` counts it."""
        if t in self.cds:
            segments = sorted(self.cds[t])
            first = segments[0] if segments[0][2] == "+" else segments[-1]
            phase = int(first[3]) if first[3].isdigit() else 0
            return sum(e - s + 1 for s, e, _, _ in segments) + max(phase, 0)
        return self.exons.get(t, 0)

    def before_the_rank(self, t: str) -> tuple:
        return (not self.mane[t], not self.canonical[t], self.tsl[t], self.biotype[t] != "protein_coding")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--by", choices=["variant", "allele", "allele-gene"], required=True)
    ap.add_argument("--gff3", required=True)
    ap.add_argument("vep")
    ap.add_argument("fastvep")
    args = ap.parse_args()

    group = {
        "variant": lambda r: (r["Uploaded_variation"],),
        "allele": lambda r: (r["Uploaded_variation"], r["Allele"]),
        "allele-gene": lambda r: (r["Uploaded_variation"], r["Allele"], r["Gene"]),
    }[args.by]
    vep, fv = read_table(args.vep), read_table(args.fastvep)
    contigs = {r["Location"].split(":")[0] for r in vep + fv}
    gff = Gff3(args.gff3, contigs)

    def groups(rows):
        g = collections.defaultdict(list)
        for r in rows:
            if r["Feature"] != "-":  # a row with no transcript is not a pick between transcripts
                g[group(r)].append(r)
        return g

    vg, fg = groups(vep), groups(fv)
    total, differing = 0, collections.Counter()
    for key, vrows in vg.items():
        frows = fg.get(key)
        if not frows:
            continue
        total += 1
        vpicked = {r["Feature"] for r in vrows if r["PICK"] == "1"}
        fpicked = {r["Feature"] for r in frows if r["PICK"] == "1"}
        if vpicked == fpicked:
            continue
        if len(vpicked) != 1 or len(fpicked) != 1:
            differing["other: more than one transcript flagged, or none"] += 1
            continue
        (a,), (b,) = vpicked, fpicked
        by_transcript = {r["Feature"]: r for r in vrows}
        if b not in by_transcript:
            differing[f"not a row of VEP's: its type is {gff.kind.get(b, '?')}"] += 1
        elif a not in by_transcript:
            differing["other: VEP's transcript is not in its own group"] += 1
        elif gff.before_the_rank(a) != gff.before_the_rank(b):
            differing["a criterion before the rank differs (MANE, canonical, TSL, biotype)"] += 1
        elif by_transcript[a]["Consequence"] != by_transcript[b]["Consequence"]:
            what = "has a CCDS" if gff.ccds.get(b) and not gff.ccds.get(a) else "has no CCDS that VEP's reader could miss"
            differing[f"a better term: fastVEP's transcript {what}"] += 1
        elif gff.length(a) > gff.length(b):
            differing["the longer one: VEP's transcript is longer (" + ("coding sequence" if a in gff.cds else "exons") + ")"] += 1
        else:
            differing[f"other: the same consequence and VEP's is not longer ({gff.length(a)} against {gff.length(b)})"] += 1

    print(f"--by {args.by}: {total} groups with a transcript row in both tables, {sum(differing.values())} with another transcript flagged")
    for what, n in sorted(differing.items(), key=lambda kv: (-kv[1], kv[0])):
        print(f"  {n:6d}  {what}")

    def featureless(rows):
        """The variants all of whose rows have no feature, their rows, and how many of those rows are flagged."""
        by_variant = collections.defaultdict(list)
        for r in rows:
            by_variant[r["Uploaded_variation"]].append(r)
        free = {v: rs for v, rs in by_variant.items() if all(r["Feature"] == "-" for r in rs)}
        return len(free), sum(len(rs) for rs in free.values()), sum(r["PICK"] == "1" for rs in free.values() for r in rs)

    for name, rows in (("Ensembl VEP", vep), ("fastVEP", fv)):
        variants, n_rows, flagged = featureless(rows)
        print(f"  {name}: {variants} variants with no transcript, {n_rows} rows, {flagged} of them flagged")


if __name__ == "__main__":
    main()
