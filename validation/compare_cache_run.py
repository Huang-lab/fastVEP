#!/usr/bin/env python3
"""Compare a tab file written by Ensembl VEP from its own cache with one written by fastVEP over the matching GFF3.

usage: compare_cache_run.py VEP.tab FASTVEP.tab OUTDIR
                            [--vep-where COL=VAL ...] [--vep-where-not COL=VAL ...] [--drop-tx-where COL=VAL ...]

Three things are counted, and each is a different question:

  1. rows: how many each file has, how many pair on (Uploaded_variation, Allele, Feature), and what the
     rows that do not pair are (fastVEP's alone are listed by Feature_type and BIOTYPE);
  2. cells: for every column both files have, the paired rows whose cell differs, the commonest
     VEP -> fastVEP values, and one TSV per column in OUTDIR;
  3. transcripts: for the transcript-level columns (Gene, SYMBOL, SYMBOL_SOURCE, HGNC_ID, BIOTYPE, CANONICAL,
     FLAGS, ...), the number of transcripts whose set of values over all of their rows differs, which is not
     dominated by the transcripts that have many rows.

--drop-tx-where drops every row of a transcript that has one row with COL=VAL in VEP's file, in both files:
`--drop-tx-where BAM_EDIT=OK` leaves out the RefSeq transcripts VEP's cache edits with an alignment (`_rna_edit`),
which a GFF3 cannot reproduce. --vep-where and --vep-where-not keep or drop single rows of VEP's file.

`validation/run_cache_oracle.sh` makes both files and calls this.
"""
import collections
import os
import sys

TX_COLS = ["Gene", "Feature_type", "STRAND", "FLAGS", "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "BIOTYPE", "CANONICAL"]


def parse_args(argv):
    pos, where, where_not, drop_tx = [], [], [], []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--vep-where":
            where.append(argv[i + 1]); i += 2
        elif a == "--vep-where-not":
            where_not.append(argv[i + 1]); i += 2
        elif a == "--drop-tx-where":
            drop_tx.append(argv[i + 1]); i += 2
        else:
            pos.append(a); i += 1
    return pos, where, where_not, drop_tx


def rows(path):
    hdr = None
    with open(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            c = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = [h.lstrip("#") for h in c]
                yield hdr
                continue
            yield c


def main():
    pos, where, where_not, drop_tx = parse_args(sys.argv[1:])
    vep_path, fv_path, out = pos[:3]
    os.makedirs(out, exist_ok=True)
    it = rows(vep_path)
    vh = next(it)
    keep = [(vh.index(k), v) for k, _, v in (w.partition("=") for w in where)]
    drop = [(vh.index(k), v) for k, _, v in (w.partition("=") for w in where_not)]
    ki = [vh.index(c) for c in ("Uploaded_variation", "Allele", "Feature")]
    dropped_tx = set()
    if drop_tx:
        specs = [(vh.index(k), v) for k, _, v in (w.partition("=") for w in drop_tx)]
        for c in rows(vep_path):
            if c[0].startswith("Uploaded") or c[0] == "Uploaded_variation":
                continue
            if any(c[i] == v for i, v in specs):
                dropped_tx.add(c[ki[2]])
        print(f"transcripts dropped ({'; '.join(drop_tx)}): {len(dropped_tx)}")
    intern = sys.intern
    vep = {}
    vtx = collections.defaultdict(lambda: collections.defaultdict(set))
    n_vep_all = 0
    for c in it:
        n_vep_all += 1
        if c[ki[2]] in dropped_tx:
            continue
        if any(c[i] != v for i, v in keep) or any(c[i] == v for i, v in drop):
            continue
        key = (c[ki[0]], c[ki[1]], c[ki[2]])
        vep[key] = tuple(intern(x) for x in c)
        for col in TX_COLS:
            if col in vh:
                vtx[c[ki[2]]][col].add(c[vh.index(col)])
    it = rows(fv_path)
    fh = next(it)
    shared = [c for c in vh if c in fh and c not in ("Uploaded_variation", "Allele", "Feature")]
    vi = {c: vh.index(c) for c in shared}
    fi = {c: fh.index(c) for c in shared}
    fki = [fh.index(c) for c in ("Uploaded_variation", "Allele", "Feature")]
    diffs = collections.defaultdict(list)
    ftx = collections.defaultdict(lambda: collections.defaultdict(set))
    only_fv = collections.Counter()
    only_fv_examples = {}
    paired = 0
    n_fv = 0
    seen = set()
    for c in it:
        key = (c[fki[0]], c[fki[1]], c[fki[2]])
        if c[fki[2]] in dropped_tx:
            continue
        n_fv += 1
        v = vep.get(key)
        for col in TX_COLS:
            if col in fh:
                ftx[c[fki[2]]][col].add(c[fh.index(col)])
        if v is None:
            only_fv[(c[fh.index("Feature_type")] if "Feature_type" in fh else "?", c[fh.index("BIOTYPE")] if "BIOTYPE" in fh else "?")] += 1
            only_fv_examples.setdefault(c[fki[2]], key)
            continue
        seen.add(key)
        paired += 1
        for col in shared:
            if v[vi[col]] != c[fi[col]]:
                diffs[col].append((key, v[vi[col]], c[fi[col]]))
    only_vep = len(vep) - len(seen)
    print(f"VEP rows {n_vep_all} ({len(vep)} kept), fastVEP rows {n_fv}, paired {paired}, VEP-only {only_vep}, fastVEP-only {sum(only_fv.values())}")
    if only_fv:
        print("  fastVEP-only rows by (Feature_type, BIOTYPE):", only_fv.most_common(6))
        print("  fastVEP-only transcripts:", len(only_fv_examples), list(only_fv_examples)[:6])
    vep_only_tx = collections.Counter()
    for k in vep:
        if k not in seen:
            vep_only_tx[k[2]] += 1
    if vep_only_tx:
        print("  VEP-only transcripts:", len(vep_only_tx), vep_only_tx.most_common(6))
    print(f"cells compared {paired * len(shared)}; columns differing:")
    for col in shared:
        d = diffs.get(col, [])
        if not d:
            continue
        top = collections.Counter((a, b) for _, a, b in d).most_common(3)
        print(f"  {col:20} {len(d):8}   e.g. {top}")
        with open(os.path.join(out, f"diff_{col}.tsv"), "w") as f:
            f.write("Uploaded_variation\tAllele\tFeature\tVEP\tfastVEP\n")
            for k, a, b in d:
                f.write("\t".join(k) + f"\t{a}\t{b}\n")
    print("transcripts (over the transcripts both files have): how many have a different set of values per column")
    common = [t for t in vtx if t in ftx]
    print(f"  transcripts in both: {len(common)} (VEP {len(vtx)}, fastVEP {len(ftx)})")
    for col in TX_COLS:
        if col not in vh or col not in fh:
            continue
        bad = [t for t in common if vtx[t][col] != ftx[t][col]]
        if not bad:
            continue
        top = collections.Counter((",".join(sorted(vtx[t][col])), ",".join(sorted(ftx[t][col]))) for t in bad).most_common(4)
        print(f"  {col:14} {len(bad):6} transcripts   e.g. {top}")
        with open(os.path.join(out, f"tx_{col}.tsv"), "w") as f:
            f.write("Feature\tVEP\tfastVEP\n")
            for t in bad:
                f.write(f"{t}\t{','.join(sorted(vtx[t][col]))}\t{','.join(sorted(ftx[t][col]))}\n")


main()
