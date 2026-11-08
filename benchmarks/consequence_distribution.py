#!/usr/bin/env python3
"""Count the consequence calls of an annotated VCF: what Table 8 and Figure 4 of the manuscript are drawn from.

A call is one entry of the CSQ field, one (variant, allele, transcript). Two counts are written, because the table
needs one and the other is what a reader asks next:

    most severe   each call counted once, under the first term of its `Consequence`. fastVEP writes the terms of a call
                  in the order of severity (Ensembl VEP's ranking, `Consequence::rank`), so the first is the most severe.
                  The counts add up to the number of calls. This is Table 8.
    every term    each call counted under every term it carries (`intron_variant&splice_region_variant` is counted
                  under both), so the counts add up to more than the number of calls.

    fastvep annotate -i HG002.vcf --gff3 ... --fasta ... --hgvs -o - | benchmarks/consequence_distribution.py - --out-prefix OUT
    benchmarks/consequence_distribution.py annotated.vcf --out-prefix OUT

writes OUT_most_severe.csv and OUT_every_term.csv (`consequence,count`, most frequent first) and OUT_run.txt, the totals and
the md5 of the VCF it read, so the counts can be tied to the bytes of the run that made them. Reading the VCF from a
pipe means the 10 GB output of a whole genome never has to be kept.

Standard library only.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from collections import Counter


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vcf", help="annotated VCF, or - for standard input")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--csq-tag", default="CSQ", help="the INFO key of the annotation (default CSQ)")
    args = ap.parse_args()

    src = sys.stdin.buffer if args.vcf == "-" else open(args.vcf, "rb")
    digest = hashlib.md5()
    column = None
    first: Counter[str] = Counter()
    every: Counter[str] = Counter()
    variants = calls = without_csq = 0
    key = (args.csq_tag + "=").encode()
    for line in src:
        digest.update(line)
        if line.startswith(b"#"):
            if line.startswith(b"##INFO=<ID=" + args.csq_tag.encode() + b","):
                text = line.decode()
                fields = text.split("Format: ", 1)[1].split('"')[0].split("|")
                column = fields.index("Consequence")
            continue
        if column is None:
            sys.exit("no ##INFO header for the CSQ field before the first record")
        variants += 1
        info = line.split(b"\t", 8)[7]
        start = info.find(key)
        if start < 0 or (start > 0 and info[start - 1:start] != b";"):
            without_csq += 1
            continue
        end = info.find(b";", start)
        value = info[start + len(key):end if end >= 0 else len(info)]
        for entry in value.split(b","):
            terms = entry.split(b"|", column + 1)[column].decode()
            calls += 1
            names = terms.split("&")
            first[names[0]] += 1
            for name in names:
                every[name] += 1
    if args.vcf != "-":
        src.close()

    for suffix, counter in (("most_severe", first), ("every_term", every)):
        with open(f"{args.out_prefix}_{suffix}.csv", "w", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(["consequence", "count"])
            for term, k in sorted(counter.items(), key=lambda kv: (-kv[1], kv[0])):
                w.writerow([term, k])
    summary = {"variants": variants, "calls": calls, "records without a CSQ field": without_csq,
               "distinct terms as the most severe": len(first), "distinct terms in all": len(every), "md5 of the VCF read": digest.hexdigest()}
    for key, value in summary.items():
        print(f"{key}: {value}")
    with open(f"{args.out_prefix}_run.txt", "w") as f:
        for key, value in summary.items():
            f.write(f"{key}: {value}\n")


if __name__ == "__main__":
    main()
