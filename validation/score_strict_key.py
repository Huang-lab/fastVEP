#!/usr/bin/env python3
"""The agreement of two tab outputs by the key of the vep-rs paper (bioRxiv 10.64898/2026.09.22.753614): the F1 of the sets of
distinct tuples (Location, Allele, Feature, Feature_type, the sorted set of Consequence terms), one score for each pair of files.

    validation/score_strict_key.py REFERENCE_DIR OTHER_DIR

REFERENCE_DIR holds the tab files that Ensembl VEP wrote and OTHER_DIR those that fastVEP wrote for the same inputs, under the same
names (`*.tab` or `*.txt`). It prints, for each name, the rows of each, the F1 of the key above and, beside it, the F1 of a key
that leaves out the allele and the position (the same (variant, feature, terms) of each), which tells a different spelling of the
allele or of the span from a different set of terms.

An unmatched row counts against the score in both directions (precision and recall), and a pair of empty files scores 1.0.
`Feature_type` of a row without a feature is taken as `-`, as it is in the paper's script.

Standard library only.
"""

import pathlib
import sys


def read(path):
    """The rows of a tab file as dicts, from its `#Uploaded_variation` header; the `##` lines are not rows.

    `None` for a file that has no such header: Ensembl VEP writes `<output>_warnings.txt` beside its output.
    """
    rows, header = [], None
    for line in open(path):
        if line.startswith("##"):
            continue
        if line.startswith("#Uploaded_variation"):
            header = line.lstrip("#").rstrip("\n").split("\t")
        elif not line.startswith("#") and header:
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows if header else None


def terms(row):
    return ",".join(sorted(set(row["Consequence"].split(","))))


def strict(row):
    feature_type = "-" if row["Feature"] == "-" else row["Feature_type"]
    return (row["Location"], row["Allele"], row["Feature"], feature_type, terms(row))


def without_allele_and_location(row):
    feature_type = "-" if row["Feature"] == "-" else row["Feature_type"]
    return (row["Uploaded_variation"], row["Feature"], feature_type, terms(row))


def f1(reference, other):
    if not reference and not other:
        return 1.0
    shared = len(reference & other)
    precision = shared / len(other) if other else 0.0
    recall = shared / len(reference) if reference else 0.0
    return 2 * precision * recall / (precision + recall) if precision + recall else 0.0


def main():
    reference_dir, other_dir = map(pathlib.Path, sys.argv[1:3])
    names = sorted(p.name for p in reference_dir.iterdir() if p.suffix in (".tab", ".txt"))
    for name in names:
        reference = read(reference_dir / name)
        if reference is None:
            continue  # not a tab output
        other = read(other_dir / name)
        if other is None:
            sys.exit(f"{other_dir / name} is not a tab output: it has no #Uploaded_variation header")
        scores = [
            f1({key(r) for r in reference}, {key(r) for r in other})
            for key in (strict, without_allele_and_location)
        ]
        print(f"{pathlib.Path(name).stem:50s} reference={len(reference):7d} other={len(other):7d} "
              f"strict={scores[0]:.4f} without_allele_and_location={scores[1]:.4f}")


if __name__ == "__main__":
    main()
