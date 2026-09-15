#!/usr/bin/env python3
"""Compare fastVEP and Ensembl VEP row by row, and dump every disagreement.

A row is one (variant, allele, transcript). The counts this prints are the ones
in [`docs/VEP_DIVERGENCE.md`](../docs/VEP_DIVERGENCE.md); the per-field TSVs it
writes beside them are what the categories in that document are counted from.

Two things about the matching are worth knowing, because both cost real rows:

* **Rows are paired by ALT ordinal, not by the allele string.** VEP spells a
  deletion `-` and so does fastVEP, which makes every allele of a multi-allelic
  deletion site look alike; keying on the string drops them from the comparison
  instead of disagreeing. Run VEP with `--allele_number` and this reads
  `ALLELE_NUM` to put each entry back on its own ALT - VEP does not emit them in
  ALT order. Without it, file order is used, which is ALT order for both tools
  on a single-ALT record and unreliable on a multi-allelic one.

* **A transcript only one tool annotates is counted, not compared.** VEP's
  `--gff` mode ignores `ncRNA_gene` records, so it has no rows at all for those
  transcripts, and a variant whose only neighbours are ncRNA genes comes out
  `intergenic_variant` there and `upstream_gene_variant` here. Those are
  reported as the two "rows only in" lines rather than mixed into the rates.
"""

import argparse
import re
from collections import defaultdict, Counter

CODING_TERMS = {
    "stop_gained",
    "stop_lost",
    "start_lost",
    "frameshift_variant",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "synonymous_variant",
    "protein_altering_variant",
    "stop_retained_variant",
    "start_retained_variant",
    "incomplete_terminal_codon_variant",
    "coding_sequence_variant",
}

FIELDS = ("cons", "impact", "aa", "codons", "hgvsc", "hgvsp", "splice")


def parse(path):
    """{(chrom, pos, ref, alt): {transcript: [csq dict, ...]}}"""
    columns, records = [], {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                columns = re.search(r"Format: ([^\"]+)", line).group(1).split("|")
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            csq = next(
                (kv[4:] for kv in fields[7].split(";") if kv.startswith("CSQ=")), ""
            )
            by_transcript = defaultdict(list)
            for entry in filter(None, csq.split(",")):
                values = entry.split("|")
                row = {
                    columns[i]: (values[i] if i < len(values) else "")
                    for i in range(len(columns))
                }
                by_transcript[row.get("Feature", "")].append(row)
            records[(fields[0], fields[1], fields[3], fields[4])] = by_transcript
    return records


def pair_up(ours, theirs):
    """One (allele, transcript) row of ours against the same row of theirs."""
    if theirs and theirs[0].get("ALLELE_NUM", ""):
        theirs = sorted(theirs, key=lambda row: int(row["ALLELE_NUM"]))
    if len(ours) == len(theirs):
        return list(zip(ours, theirs)), 0, 0
    # Unequal counts mean one tool split or dropped an allele; fall back to the
    # allele string and count what it cannot place.
    mine = {row.get("Allele", ""): row for row in ours}
    yours = {row.get("Allele", ""): row for row in theirs}
    shared = mine.keys() & yours.keys()
    return (
        [(mine[k], yours[k]) for k in shared],
        len(mine) - len(shared),
        len(yours) - len(shared),
    )


def compare(fastvep, vep, out_prefix):
    matched = coding = only_ours = only_theirs = 0
    counts = Counter()
    dumps = {name: [] for name in FIELDS}

    for key in fastvep.keys() & vep.keys():
        ours, theirs = fastvep[key], vep[key]
        for transcript in ours.keys() | theirs.keys():
            mine, yours = ours.get(transcript, []), theirs.get(transcript, [])
            if not mine:
                only_theirs += len(yours)
                continue
            if not yours:
                only_ours += len(mine)
                continue
            pairs, unpaired_ours, unpaired_theirs = pair_up(mine, yours)
            only_ours += unpaired_ours
            only_theirs += unpaired_theirs

            for ours_row, vep_row in pairs:
                matched += 1
                our_terms = set(filter(None, ours_row["Consequence"].split("&")))
                vep_terms = set(filter(None, vep_row["Consequence"].split("&")))
                is_coding = bool((our_terms | vep_terms) & CODING_TERMS)
                coding += is_coding
                where = [
                    f"{key[0]}:{key[1]} {key[2]}>{key[3]}",
                    transcript,
                    ours_row.get("SYMBOL", ""),
                    ours_row.get("Allele", ""),
                    vep_row.get("Allele", ""),
                ]
                terms = ["&".join(sorted(our_terms)), "&".join(sorted(vep_terms))]

                def note(field, ours_value, vep_value):
                    counts[field] += 1
                    dumps[field].append(where + [ours_value, vep_value] + terms)

                if our_terms != vep_terms:
                    counts["cons_coding"] += is_coding
                    note("cons", terms[0], terms[1])
                our_splice = {t for t in our_terms if "splice" in t}
                vep_splice = {t for t in vep_terms if "splice" in t}
                if our_splice != vep_splice:
                    note("splice", "&".join(sorted(our_splice)), "&".join(sorted(vep_splice)))
                if ours_row["IMPACT"] != vep_row["IMPACT"]:
                    note("impact", ours_row["IMPACT"], vep_row["IMPACT"])
                if is_coding and ours_row["Amino_acids"] != vep_row["Amino_acids"]:
                    note("aa", ours_row["Amino_acids"], vep_row["Amino_acids"])
                if is_coding and ours_row["Codons"] != vep_row["Codons"]:
                    note("codons", ours_row["Codons"], vep_row["Codons"])
                if ours_row["HGVSc"] != vep_row["HGVSc"]:
                    note("hgvsc", ours_row["HGVSc"], vep_row["HGVSc"])
                if ours_row["HGVSp"] != vep_row["HGVSp"]:
                    note("hgvsp", ours_row["HGVSp"], vep_row["HGVSp"])

    print(f"matched rows          {matched}")
    print(f"coding rows           {coding}  ({100.0 * coding / max(matched, 1):.1f} %)")
    print(f"rows only in fastVEP  {only_ours}")
    print(f"rows only in VEP      {only_theirs}")
    print()
    print(f"{'field':24s} {'scope':8s} {'disagreeing':>11s} {'agreement':>11s}")
    rows = [
        ("Consequence terms", "coding", counts["cons_coding"], coding),
        ("Amino_acids", "coding", counts["aa"], coding),
        ("Codons", "coding", counts["codons"], coding),
        ("Splice terms", "all", counts["splice"], matched),
        ("Whole consequence set", "all", counts["cons"], matched),
        ("IMPACT", "all", counts["impact"], matched),
        ("HGVSc", "all", counts["hgvsc"], matched),
        ("HGVSp", "all", counts["hgvsp"], matched),
    ]
    for label, scope, bad, total in rows:
        rate = 100.0 * (total - bad) / total if total else 0.0
        print(f"{label:24s} {scope:8s} {bad:11d} {rate:10.3f} %")

    if out_prefix:
        for name, rows in dumps.items():
            with open(f"{out_prefix}.{name}.tsv", "w") as out:
                out.write(
                    "location\ttranscript\tsymbol\tfastvep_allele\tvep_allele"
                    "\tfastvep\tvep\tfastvep_consequence\tvep_consequence\n"
                )
                for row in rows:
                    out.write("\t".join(row) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("fastvep_vcf")
    ap.add_argument("vep_vcf")
    ap.add_argument(
        "--dump",
        metavar="PREFIX",
        help="write PREFIX.<field>.tsv with every disagreeing row",
    )
    args = ap.parse_args()
    compare(parse(args.fastvep_vcf), parse(args.vep_vcf), args.dump)


if __name__ == "__main__":
    main()
