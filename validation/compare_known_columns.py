#!/usr/bin/env python3
"""Compare the known-variant columns of a tab file written by Ensembl VEP with one written by fastVEP.

usage: compare_known_columns.py VEP.tab FASTVEP.tab [--show N] [--out DIR] [--unordered COL,COL]

The columns are `Existing_variation`, `AF`, `AFR_AF` ... `SAS_AF`, `gnomADe_*`, `gnomADg_*`,
`MAX_AF`, `MAX_AF_POPS`, `CLIN_SIG`, `SOMATIC`, `PHENO` and `PUBMED`: the ones both files have.

What these columns say is a function of the input variant and its allele, not of the transcript, so
the rows of each file are collapsed to one entry per (Uploaded_variation, Allele): the set of
distinct tuples of these cells over the rows of that pair. The two tools are run with different
transcript sources and write different rows for the same variant, which is the business of
`compare_tab.py`; here a pair is paired whenever both files have a row for it. A file that gives
two different tuples for the rows of one pair is reported as inconsistent: VEP's output hash is
built per row, so that is a finding about the program, not about the comparison.

`--unordered` names columns whose cells are `;`-separated lists that VEP writes in the order of a
Perl hash (`CLIN_SIG` with two distinct allele-specific values), which differs from run to run; their
parts are compared as a set. `CLIN_SIG` is always treated this way.

Prints the pairs compared, those in one file only, and for each column how many pairs differ, with the
commonest differing values. With --out, one TSV per differing column lists every pair that differs.
"""
import argparse
import collections
import os
import sys

KNOWN = ("Existing_variation AF AFR_AF AMR_AF EAS_AF EUR_AF SAS_AF "
         "gnomADe_AF gnomADe_AFR_AF gnomADe_AMR_AF gnomADe_ASJ_AF gnomADe_EAS_AF gnomADe_FIN_AF "
         "gnomADe_MID_AF gnomADe_NFE_AF gnomADe_REMAINING_AF gnomADe_SAS_AF "
         "gnomADg_AF gnomADg_AFR_AF gnomADg_AMI_AF gnomADg_AMR_AF gnomADg_ASJ_AF gnomADg_EAS_AF "
         "gnomADg_FIN_AF gnomADg_MID_AF gnomADg_NFE_AF gnomADg_REMAINING_AF gnomADg_SAS_AF "
         "MAX_AF MAX_AF_POPS CLIN_SIG SOMATIC PHENO PUBMED").split()


def collapse(path, columns, unordered):
    """{(variation, allele): set of tuples of cells}, and the header's columns."""
    header = None
    pairs = collections.defaultdict(set)
    rows = 0
    with open(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            line = line.rstrip("\n")
            if line.startswith("#"):
                header = line[1:].split("\t")
                at = {c: header.index(c) for c in ("Uploaded_variation", "Allele")}
                known = [(c, header.index(c)) for c in columns if c in header]
                continue
            cells = line.split("\t")
            rows += 1
            values = []
            for name, i in known:
                v = cells[i]
                if name in unordered and ";" in v:
                    v = ";".join(sorted(v.split(";")))
                values.append(v)
            # VEP upper-cases the alleles of its input; fastVEP writes them as given
            pairs[(cells[at["Uploaded_variation"]], cells[at["Allele"]].upper())].add(tuple(values))
    return header, [c for c, _ in known], pairs, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vep")
    ap.add_argument("fastvep")
    ap.add_argument("--show", type=int, default=5)
    ap.add_argument("--out", default="")
    ap.add_argument("--unordered", default="")
    a = ap.parse_args()
    unordered = {c for c in a.unordered.split(",") if c} | {"CLIN_SIG"}

    vh, vcols, vpairs, vrows = collapse(a.vep, KNOWN, unordered)
    fh, fcols, fpairs, frows = collapse(a.fastvep, KNOWN, unordered)
    print(f"VEP     : {vrows:8d} rows, {len(vpairs):7d} (variant, allele) pairs, known columns {len(vcols)}")
    print(f"fastVEP : {frows:8d} rows, {len(fpairs):7d} (variant, allele) pairs, known columns {len(fcols)}")
    print("\n== columns ==")
    print("VEP only    :", [c for c in vcols if c not in fcols])
    print("fastVEP only:", [c for c in fcols if c not in vcols])
    shared = [c for c in vcols if c in fcols]
    print("shared columns in the same order:", [c for c in vcols if c in fcols] == [c for c in fcols if c in vcols])
    vi = {c: i for i, c in enumerate(vcols)}
    fi = {c: i for i, c in enumerate(fcols)}

    print("\n== pairs ==")
    both = sorted(set(vpairs) & set(fpairs))
    only_vep = sorted(set(vpairs) - set(fpairs))
    only_ours = sorted(set(fpairs) - set(vpairs))
    print(f"paired: {len(both)}, VEP only: {len(only_vep)}, fastVEP only: {len(only_ours)}")
    for label, keys in (("VEP only", only_vep), ("fastVEP only", only_ours)):
        for k in keys[: a.show]:
            print(f"  {label}: {k}")
    for label, pairs in (("VEP", vpairs), ("fastVEP", fpairs)):
        inconsistent = [k for k, v in pairs.items() if len(v) > 1]
        print(f"{label}: {len(inconsistent)} pairs have rows that disagree among themselves")
        for k in inconsistent[: a.show]:
            print(f"    {k}: {sorted(pairs[k])[:3]}")

    diffs = collections.defaultdict(list)
    cells = 0
    nonempty = 0
    for k in both:
        v, f = sorted(vpairs[k])[0], sorted(fpairs[k])[0]
        for c in shared:
            cells += 1
            nonempty += v[vi[c]] != "-"
            if v[vi[c]] != f[fi[c]]:
                diffs[c].append((k, v[vi[c]], f[fi[c]]))
    total = sum(len(x) for x in diffs.values())
    pairs_differing = len({k for lst in diffs.values() for k, _, _ in lst})
    print(f"\n== cells: {len(both)} paired (variant, allele), {cells} cells compared "
          f"({nonempty} of them not `-` in VEP), {total} differ in {pairs_differing} pairs ==")
    for c in shared:
        if c in diffs:
            print(f"\n[{c}] {len(diffs[c])} differ")
            tally = collections.Counter((v, f) for _, v, f in diffs[c])
            for (v, f), n in tally.most_common(a.show):
                print(f"   VEP {v!r:42.42} fastVEP {f!r:42.42} x{n}")
            for k, v, f in diffs[c][: a.show]:
                print(f"     e.g. {k}")
    if a.out:
        os.makedirs(a.out, exist_ok=True)
        for c, lst in diffs.items():
            with open(os.path.join(a.out, f"diff_{c}.tsv"), "w") as o:
                o.write("Uploaded_variation\tAllele\tVEP\tfastVEP\n")
                for k, v, f in lst:
                    o.write("\t".join([*k, v, f]) + "\n")
    sys.exit(1 if total or only_vep or only_ours else 0)


main()
