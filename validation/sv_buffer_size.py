#!/usr/bin/env python3
"""Does what Ensembl VEP writes for a structural variant depend on the records it reads with it?

    sv_buffer_size.py --fastvep BIN --gff3 GENES.gff3 --vep-gff REF_DIR/GENES.gff3.gz --fasta REF_DIR/GENOME.fa \
                      --image IMAGE --workdir DIR VCF

Ensembl VEP reads its input 5,000 records at a time (`--buffer_size`), and in `--gff` mode it reads the transcripts of a stretch of
records for the span of all of them (`filter_features_by_min_max` of `AnnotationSource.pm`), and pairs a structural variant with a transcript
by the numbers of its interval whatever contig the transcript is on (`_chr_matches_for_vf_overlap` of `InputBuffer.pm`). So what a record gets can
depend on the records buffered with it. This runs fastVEP once, and VEP twice over the same VCF, with its own buffer and with a buffer of one
record, in the way that `regen.sh` and the head-to-head do (statistics on, tab output), and for each of the two VEP runs prints the F1 of the
strict key of the vep-rs paper (`score_strict_key.py`), the tuples that only one tool writes, and those split by the kind of record they
belong to:

  single or symbolic       a breakend with no mate position (`A.`, `.A`), and `<BND>`
  a mate on its own contig `N[chr1:2500[` for a record on chr1
  a mate on another contig `N[chr5:2500[` for a record on chr1

REF_DIR holds the bgzipped, indexed GFF3 and the FASTA (with its .fai) that VEP is given (mounted read-only as /ref), and GENES.gff3 is the same
gene models as fastVEP reads. Needs docker and the image already pulled. WORKDIR must not exist or must be empty; nothing in it is removed.
"""
import argparse
import collections
import re
import shutil
import subprocess
import sys
from pathlib import Path


def read(path):
    header, rows = None, []
    for line in open(path):
        if line.startswith("##"):
            continue
        cells = line.rstrip("\n").split("\t")
        if line.startswith("#"):
            header = [c.lstrip("#") for c in cells]
        elif header:
            rows.append(dict(zip(header, cells)))
    return rows


def key(r):
    """The strict key of the vep-rs paper: location, allele, feature, its type, and the sorted set of terms."""
    feature_type = "-" if r["Feature"] == "-" else r["Feature_type"]
    return (r["Location"], r["Allele"], r["Feature"], feature_type, ",".join(sorted(set(r["Consequence"].split(",")))))


def kinds(vcf):
    """record id -> the kind of its ALT: single or symbolic, a mate on its own contig, a mate on another."""
    out = {}
    for line in open(vcf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        mate = re.search(r"[\[\]]([^:\[\]]+):", f[4])
        if "<BND>" in f[4] or not mate:
            out[f[2]] = "single or symbolic"
        else:
            out[f[2]] = "a mate on its own contig" if mate.group(1) == f[0] else "a mate on another contig"
    return out


def by_record(rows):
    g = collections.defaultdict(set)
    for r in rows:
        g[r["Uploaded_variation"]].add(key(r))
    return g


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fastvep", required=True)
    ap.add_argument("--gff3", required=True, help="the gene models as fastVEP reads them")
    ap.add_argument("--vep-gff", required=True, help="the bgzipped, tabix-indexed GFF3 that VEP reads, in the directory that holds the FASTA")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--image", required=True, help="the Ensembl VEP 116.3 image, by digest")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--runs", type=int, default=1, help="how many times VEP is run with each buffer (its rows for a breakend whose mate is on another contig vary from run to run)")
    ap.add_argument("vcf")
    args = ap.parse_args()

    work = Path(args.workdir).resolve()
    work.mkdir(parents=True, exist_ok=True)
    if any(work.iterdir()):
        sys.exit(f"{work} is not empty")
    vcf = Path(args.vcf).resolve()
    ref_dir = Path(args.vep_gff).resolve().parent
    if Path(args.fasta).resolve().parent != ref_dir:
        sys.exit("the FASTA must be in the directory of the GFF3 that VEP reads")
    (work / "in").mkdir()
    shutil.copyfile(vcf, work / "in" / vcf.name)

    subprocess.run([args.fastvep, "annotate", "-i", str(vcf), "--gff3", args.gff3, "--fasta", args.fasta, "--output-format", "tab",
                    "-o", str(work / "fastvep.tab"), "--no-progress"], check=True, capture_output=True)
    runs = {}
    for label, size in (("default buffer", None), ("--buffer_size 1", 1)):
        for n in range(1, args.runs + 1):
            out = work / f"vep_{'default' if size is None else size}_{n}.tab"
            cmd = ["docker", "run", "--rm", "-e", "VEP_DATABASE=0", "-v", f"{work / 'in'}:/in:ro", "-v", f"{ref_dir}:/ref:ro", "-v", f"{work}:/out",
                   args.image, "vep", "--input_file", f"/in/{vcf.name}", "--gff", f"/ref/{Path(args.vep_gff).name}",
                   "--fasta", f"/ref/{Path(args.fasta).name}", "--output_file", f"/out/{out.name}", "--tab", "--force_overwrite"]
            if size is not None:
                cmd += ["--buffer_size", str(size)]
            subprocess.run(cmd, check=True, capture_output=True)
            runs[f"{label}, run {n}" if args.runs > 1 else label] = out

    fv_rows = read(work / "fastvep.tab")
    fv_by = by_record(fv_rows)
    fv_keys = {k for ks in fv_by.values() for k in ks}
    kind = kinds(vcf)
    print(f"{vcf.name}: {len(kind)} records; fastVEP {len(fv_rows)} rows, {len(fv_keys)} distinct tuples")
    for label, path in runs.items():
        rows = read(path)
        vp_by = by_record(rows)
        vp_keys = {k for ks in vp_by.values() for k in ks}
        common = len(vp_keys & fv_keys)
        f1 = 2 * common / (len(vp_keys) + len(fv_keys)) if (vp_keys or fv_keys) else 1.0
        print(f"  VEP, {label}: {len(rows)} rows, {len(vp_keys)} tuples, {common} in common, strict F1 {f1:.4f}; "
              f"only VEP {len(vp_keys - fv_keys)}, only fastVEP {len(fv_keys - vp_keys)}")
        only_v, only_f = collections.Counter(), collections.Counter()
        for rec in set(vp_by) | set(fv_by):
            k = kind.get(rec, "?")
            only_v[k] += len(vp_by.get(rec, set()) - fv_by.get(rec, set()))
            only_f[k] += len(fv_by.get(rec, set()) - vp_by.get(rec, set()))
        for k in sorted(set(only_v) | set(only_f)):
            print(f"      {k:26s} only VEP {only_v[k]:6d}   only fastVEP {only_f[k]:6d}")


if __name__ == "__main__":
    main()
