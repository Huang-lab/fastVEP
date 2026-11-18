#!/usr/bin/env python3
"""Ensembl VEP 116.3 against fastVEP on a freshly generated world of structural variants, cell by cell.

    run_sv_world.py --fastvep BIN --workdir DIR [--seeds 11 12 13 ...] [--jobs N]

The golden of `crates/fastvep-cli/fixtures/vep_tab_sv` is one world (the seed 185 of its
`gen_world.py`), and a test that reads one golden cannot tell a rule that is right from a rule that
was fitted to it. This writes a world for each seed given, runs the VEP image that `regen.sh` there
pins over it (in `--gff` mode, with its statistics, as the golden is made), runs `BIN annotate` over
the same files with the same columns, and compares what `vep_tab_sv_golden.rs` compares:

  * the rows: the same (record, allele, feature) keys, and nothing else;
  * their order: the records as the file has them, and for each its rows in VEP's order;
  * the 23 columns of `COMPARED`, cell by cell, on every row.

Every seed gets a directory of its own in DIR (the world, both tables, `differences.tsv`), and DIR is
never emptied or removed here: a run refuses a DIR that already holds a seed's directory.

The cells that are known to differ are those of `KNOWN_IN_THE_GOLDEN`: a tandem repeat with several
alleles is a record with several ALTs, and VEP 116 describes the alleles of such a record as the record
(`docs/VEP_DIVERGENCE.md`, Part 1, section 11). They are listed as known and counted apart; any other
difference is printed and makes the exit status 1.

Needs docker (the image `regen.sh` names, already pulled), python3, bgzip and tabix.
"""
import argparse
import concurrent.futures
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIXTURE = HERE.parent / "crates" / "fastvep-cli" / "fixtures" / "vep_tab_sv"

# The columns the golden test compares. `Existing_variation`, `FLAGS`, `SYMBOL_SOURCE` and `HGNC_ID`
# are not: a `--gff` run of VEP has no variation source, no `cds_start_NF`, and no symbol source.
COMPARED = [
    "Location", "Consequence", "IMPACT", "cDNA_position", "CDS_position", "Protein_position",
    "Amino_acids", "Codons", "DISTANCE", "STRAND", "Gene", "Feature_type", "EXON", "INTRON",
    "VARIANT_CLASS", "REF_ALLELE", "HGVSc", "HGVSp", "HGVSg", "HGVS_OFFSET", "SYMBOL", "BIOTYPE",
    "CANONICAL",
]

# What `regen.sh` gives VEP, and the same columns for fastVEP.
VEP_FLAGS = ("--symbol --biotype --canonical --numbers --variant_class --show_ref_allele "
             "--hgvs --hgvsg --total_length").split()
FASTVEP_FLAGS = ("--symbol --biotype --canonical --numbers --variant-class --show-ref-allele "
                 "--hgvs --hgvsg --total-length").split()
DISTANCE = "1000"  # not VEP's 5000: the world is 14 kb with ten transcripts (see regen.sh)

# A tandem repeat with several alleles is read as a record with several ALTs (divergence 11): the
# columns that describe the change, on the rows of the records that the generator names
# `tr_two_alleles*`.
KNOWN_IN_THE_GOLDEN = {"HGVSc", "HGVS_OFFSET", "HGVSg"}
KNOWN_RECORDS = "tr_two_alleles"


def pinned_image():
    text = (FIXTURE / "regen.sh").read_text()
    found = re.search(r'^IMAGE="([^"]+)"', text, re.M)
    if not found:
        sys.exit("the image is not named in regen.sh")
    return found.group(1)


def read_table(path):
    """(header, rows) of a VEP-style tab file: a row is a dict by column name."""
    header, rows = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith("##"):
            continue
        cells = line.split("\t")
        if line.startswith("#"):
            header = [c.lstrip("#") for c in cells]
            continue
        if len(cells) != len(header):
            sys.exit(f"{path}: a ragged row {line[:120]!r}")
        rows.append(dict(zip(header, cells)))
    if header is None:
        sys.exit(f"{path}: no header")
    return header, rows


def key(row):
    return (row["Uploaded_variation"], row["Allele"], row["Feature"])


def run(cmd, **kw):
    kw.setdefault("text", True)
    return subprocess.run(cmd, check=True, capture_output=True, **kw)


def one_seed(seed, fastvep, workdir, image):
    here = Path(workdir) / f"seed_{seed}"
    here.mkdir(parents=True)  # refuses a seed that was already run: nothing is overwritten
    run([sys.executable, "-I", str(FIXTURE / "gen_world.py"), str(here), str(seed)])
    # VEP's reader of a GFF3 wants it sorted, bgzipped and indexed.
    gff = here / "genes.gff3"
    lines = gff.read_text().splitlines()
    body = sorted((l for l in lines if not l.startswith("#")),
                  key=lambda l: (l.split("\t")[0], int(l.split("\t")[3]), int(l.split("\t")[4])))
    sorted_text = "\n".join([l for l in lines if l.startswith("#")] + body) + "\n"
    with open(here / "genes.sorted.gff3.gz", "wb") as out:
        out.write(run(["bgzip", "-c"], input=sorted_text.encode(), text=False).stdout)
    run(["tabix", "-p", "gff", str(here / "genes.sorted.gff3.gz")])

    # VEP is run as `regen.sh` runs it, with its statistics (the HTML lands beside the table).
    subprocess.run(
        ["docker", "run", "--rm", "-e", "VEP_DATABASE=0", "-v", f"{here}:/data", image, "vep",
         "--input_file", "/data/in.vcf", "--output_file", "/data/vep.tab",
         "--gff", "/data/genes.sorted.gff3.gz", "--fasta", "/data/genome.fa",
         "--distance", DISTANCE, "--force_overwrite", "--tab", *VEP_FLAGS],
        check=True, capture_output=True)
    subprocess.run(
        [fastvep, "annotate", "-i", str(here / "in.vcf"), "--gff3", str(here / "genes.gff3"),
         "--fasta", str(here / "genome.fa"), "--output-format", "tab", "-o", str(here / "fastvep.tab"),
         "--distance", DISTANCE, "--no-progress", *FASTVEP_FLAGS],
        check=True, capture_output=True)

    _, vep = read_table(here / "vep.tab")
    _, fv = read_table(here / "fastvep.tab")
    vep_keys, fv_keys = [key(r) for r in vep], [key(r) for r in fv]
    only_vep = sorted(set(vep_keys) - set(fv_keys))
    only_fv = sorted(set(fv_keys) - set(vep_keys))
    first_out_of_order = next(
        (i for i, (a, b) in enumerate(zip(vep_keys, fv_keys)) if a != b), None
    ) if len(vep_keys) == len(fv_keys) else "the row counts differ"

    by_key = {}
    for r in fv:
        by_key.setdefault(key(r), []).append(r)
    differences, known, compared = [], [], 0
    seen = {}
    for v in vep:
        k = key(v)
        n = seen.get(k, 0)
        seen[k] = n + 1
        f = by_key.get(k, [None] * (n + 1))
        if n >= len(f) or f[n] is None:
            continue
        f = f[n]
        for col in COMPARED:
            compared += 1
            if v[col] != f[col]:
                record = (k[0], k[1], k[2], col, v[col], f[col])
                (known if k[0].startswith(KNOWN_RECORDS) and col in KNOWN_IN_THE_GOLDEN else differences).append(record)
    with open(here / "differences.tsv", "w") as out:
        out.write("kind\trecord\tallele\tfeature\tcolumn\tvep\tfastvep\n")
        for kind, items in (("known", known), ("unexpected", differences)):
            for r in items:
                out.write("\t".join((kind, *r)) + "\n")
    return {
        "seed": seed, "records_in": sum(1 for l in (here / "in.vcf").read_text().splitlines() if not l.startswith("#")),
        "vep_rows": len(vep), "fastvep_rows": len(fv), "only_vep": only_vep, "only_fastvep": only_fv,
        "first_out_of_order": first_out_of_order, "cells": compared,
        "known": len(known), "unexpected": differences,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fastvep", required=True, help="the binary under test")
    ap.add_argument("--workdir", required=True, help="where the seeds' directories go (created if absent)")
    ap.add_argument("--seeds", type=int, nargs="+", default=[11, 23, 37, 41, 59, 73])
    ap.add_argument("--jobs", type=int, default=2, help="seeds run side by side")
    args = ap.parse_args()
    fastvep = str(Path(args.fastvep).resolve())
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    image = pinned_image()

    bad = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(one_seed, s, fastvep, workdir, image): s for s in args.seeds}
        results = []
        for fut in concurrent.futures.as_completed(futures):
            try:
                results.append(fut.result())
            except subprocess.CalledProcessError as e:
                sys.exit(f"seed {futures[fut]}: {' '.join(map(str, e.cmd))[:200]} failed: {(e.stderr or b'')[:400]!r}")
    for r in sorted(results, key=lambda r: r["seed"]):
        wrong = bool(r["only_vep"] or r["only_fastvep"] or r["first_out_of_order"] is not None
                     or r["unexpected"] or r["vep_rows"] != r["fastvep_rows"])
        bad += wrong
        print(f"seed {r['seed']}: {r['records_in']} records, VEP {r['vep_rows']} rows, fastVEP {r['fastvep_rows']} rows, "
              f"{r['cells']} cells compared, {r['known']} known (divergence 11), {len(r['unexpected'])} unexpected, "
              f"{len(r['only_vep'])} rows only VEP writes, {len(r['only_fastvep'])} only fastVEP writes, "
              f"first row out of order: {r['first_out_of_order']}  {'DIFFERS' if wrong else 'same'}")
        for k in r["only_vep"][:5]:
            print("   only VEP:", k)
        for k in r["only_fastvep"][:5]:
            print("   only fastVEP:", k)
        for d in r["unexpected"][:12]:
            print("   cell:", d)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
