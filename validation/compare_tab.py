#!/usr/bin/env python3
"""Compare a tab file written by Ensembl VEP with one written by fastVEP, cell by cell.

usage: compare_tab.py VEP.tab FASTVEP.tab [--ignore COL,COL] [--key COL,COL] [--show N] [--out DIR]

Rows are paired on (Uploaded_variation, Allele, Feature) unless --key says otherwise, so the
input needs a distinct name per record (an ID of its own, or a distinct VCF line, which VEP 116
names it by).
Where every record has one ALT, --key Uploaded_variation,Feature pairs a row whose Allele the
two tools write differently, and the Allele column then reports the difference.
The header is compared first: the columns each has, and whether the shared ones are in
the same order. Then every shared column of every paired row is compared, and a count of
the cells that differ is printed for each column with the commonest differing values.
With --out, one TSV per differing column lists every row that differs.

`validation/run_tab_world.sh` runs both tools over a generated world and calls this.
"""
import argparse
import collections
import os
import sys


def read(path):
    header, rows = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line[1:].split("\t")
                continue
            cells = line.split("\t")
            if header is None:
                sys.exit(f"{path}: data before a header")
            if len(cells) != len(header):
                sys.exit(f"{path}: a row of {len(cells)} cells under {len(header)} columns")
            rows.append(dict(zip(header, cells)))
    return header, rows


def key_of(columns):
    return lambda row: tuple(row[c] for c in columns)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vep")
    ap.add_argument("fastvep")
    ap.add_argument("--ignore", default="")
    ap.add_argument("--key", default="Uploaded_variation,Allele,Feature")
    ap.add_argument("--show", type=int, default=4)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    ignore = {c for c in a.ignore.split(",") if c}
    key = key_of(a.key.split(","))

    vep_header, vep_rows = read(a.vep)
    our_header, our_rows = read(a.fastvep)
    print(f"VEP     : {len(vep_rows):6d} rows, {len(vep_header)} columns")
    print(f"fastVEP : {len(our_rows):6d} rows, {len(our_header)} columns")

    vep_cmp = [c for c in vep_header if c not in ignore]
    print("\n== header ==")
    if vep_cmp == our_header:
        print("identical" + (f" (ignoring {sorted(ignore)})" if ignore else ""))
    else:
        print("VEP only    :", [c for c in vep_cmp if c not in our_header])
        print("fastVEP only:", [c for c in our_header if c not in vep_cmp])
        shared_v = [c for c in vep_cmp if c in our_header]
        shared_f = [c for c in our_header if c in vep_cmp]
        print("shared columns in the same order:", shared_v == shared_f)
    shared = [c for c in our_header if c in vep_header and c not in ignore]

    by_key = {"vep": collections.defaultdict(list), "fastvep": collections.defaultdict(list)}
    for r in vep_rows:
        by_key["vep"][key(r)].append(r)
    for r in our_rows:
        by_key["fastvep"][key(r)].append(r)
    only_vep = sorted(set(by_key["vep"]) - set(by_key["fastvep"]))
    only_ours = sorted(set(by_key["fastvep"]) - set(by_key["vep"]))
    print("\n== rows ==")
    print(f"paired keys: {len(set(by_key['vep']) & set(by_key['fastvep']))}, "
          f"VEP only: {len(only_vep)}, fastVEP only: {len(only_ours)}")
    for label, keys in (("VEP only", only_vep), ("fastVEP only", only_ours)):
        for k in keys[: a.show]:
            print(f"  {label}: {k}")

    diffs = collections.defaultdict(list)
    cells = pairs = 0
    for k in sorted(set(by_key["vep"]) & set(by_key["fastvep"])):
        vs = sorted(by_key["vep"][k], key=lambda r: [r[c] for c in shared])
        fs = sorted(by_key["fastvep"][k], key=lambda r: [r[c] for c in shared])
        for v, f in zip(vs, fs):
            pairs += 1
            for c in shared:
                cells += 1
                if v[c] != f[c]:
                    diffs[c].append((k, v[c], f[c]))
    total = sum(len(x) for x in diffs.values())
    print(f"\n== cells: {pairs} paired rows, {cells} cells compared, {total} differ ==")
    for c in shared:
        if c in diffs:
            print(f"\n[{c}] {len(diffs[c])} differ")
            tally = collections.Counter((v, f) for _, v, f in diffs[c])
            for (v, f), n in tally.most_common(a.show):
                print(f"   VEP {v!r:42.42} fastVEP {f!r:42.42} x{n}")
    if a.out:
        os.makedirs(a.out, exist_ok=True)
        for c, lst in diffs.items():
            with open(os.path.join(a.out, f"diff_{c}.tsv"), "w") as o:
                o.write("Uploaded_variation\tAllele\tFeature\tVEP\tfastVEP\n")
                for k, v, f in lst:
                    o.write("\t".join([*k, v, f]) + "\n")


main()
