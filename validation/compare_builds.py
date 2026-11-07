#!/usr/bin/env python3
"""Does one build of fastVEP write the same bytes as an earlier build, on the inputs of the comparison with Ensembl VEP?

    validation/compare_builds.py --base BIN --base-label 5048fc5 --new BIN --new-label v0.5.0-rc.1 --out DIR

Ensembl VEP's output does not depend on the build of fastVEP that it is compared with, so a comparison measured with one build (the numbers of
Part 3 of docs/VEP_DIVERGENCE.md were measured on 27 October 2026, at 5048fc5) is the comparison of another build whose output on the same inputs
is the same bytes, and only of such a build. The numbers against release 115.1 in the same document were measured at ac2e2b6 and describe that
build, which writes other bytes than the ones after it.

This runs both builds over the inputs of `run_divergence.sh` and `run_validation.sh`, over the same GRCh38 release 115 GFF3 and FASTA, in two modes:
`vcf`, with the flags of those scripts (`--hgvs --symbol --canonical`), and `tab`, with the flags of the comparison with Ensembl VEP 116.3
(`--output-format tab --symbol --biotype --canonical --numbers --variant-class --show-ref-allele --hgvs --hgvsg`), which reach columns and a
writer that the VCF run does not. The two outputs are compared byte for byte. Each build reads a transcript cache of its own, written by itself on
its first run, so the check does not take on trust that a cache written by one build is read the same by the other.

The inputs: the 173 variants of the VEP example (136 on chromosome 22 and 37 on chromosome 21), the stratified ClinVar sample of 6,600 variants,
the 1-in-5 sample of the ClinVar 2-star+ set (134,732), the 1-in-200 sample of GIAB HG002 (20,241), the 400 ClinVar in-frame deletions, and, as a
VCF only, the whole of HG002 (4,048,342 variants). The samples are cut by `sample_variants.py`, which is deterministic.

Writes identity.tsv and environment.tsv into --out. The outputs are hashed and deleted unless --keep is given (the whole genome's is 11 GB).
Exit status 1 if any run differs. Standard library only.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
HUMAN = Path(os.environ.get("HUMAN_DIR", ROOT / "test_data" / "organisms" / "human"))
CLINVAR = Path(os.environ.get("CLINVAR_2STAR", ROOT / "data" / "benchmark" / "clinvar_2star.vcf"))
GFF3 = HUMAN / "Homo_sapiens.GRCh38.115.gff3"
FASTA = HUMAN / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
HG002 = HUMAN / "human_giab_hg002_full.vcf"
# vcf: the flags of run_divergence.sh and run_validation.sh. tab: the flags of the comparison with Ensembl VEP 116.3 (Part 3 of docs/VEP_DIVERGENCE.md),
# which reach columns and a writer that the VCF run does not.
MODES = {
    "vcf": ["--hgvs", "--symbol", "--canonical"],
    "tab": ["--output-format", "tab", "--symbol", "--biotype", "--canonical", "--numbers", "--variant-class", "--show-ref-allele", "--hgvs", "--hgvsg"],
}


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(16 << 20), b""):
            h.update(block)
    return h.hexdigest()


def records(path: Path) -> int:
    with path.open("rb") as f:
        return sum(1 for line in f if not line.startswith(b"#"))


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(16 << 20), b""):
            h.update(block)
    return h.hexdigest()


def gene_model_dir(base: Path, label: str) -> tuple[Path, Path]:
    """Links to the GFF3 and the FASTA in a directory of this build's own, so that the sidecar cache it writes is its own."""
    d = base / f"gene_model_{label}"
    d.mkdir(parents=True, exist_ok=True)
    out = []
    for src in (GFF3, FASTA):
        dst = d / src.name
        if not dst.exists():
            dst.symlink_to(src)
        out.append(dst)
    fai = FASTA.with_name(FASTA.name + ".fai")
    if fai.exists() and not (d / fai.name).exists():
        (d / fai.name).symlink_to(fai)
    return out[0], out[1]


def annotate(binary: Path, gff3: Path, fasta: Path, vcf: Path, out: Path, flags: list[str]) -> float:
    t0 = time.perf_counter()
    r = subprocess.run([str(binary), "annotate", "-i", str(vcf), "--gff3", str(gff3), "--fasta", str(fasta), *flags, "-o", str(out)],
                       capture_output=True, text=True)
    wall = time.perf_counter() - t0
    if r.returncode != 0:
        sys.exit(f"FAILED (exit {r.returncode}): {binary} on {vcf}\n{r.stderr[-1000:]}")
    return wall


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True, type=Path)
    ap.add_argument("--base-label", required=True)
    ap.add_argument("--new", required=True, type=Path)
    ap.add_argument("--new-label", required=True)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--work", type=Path, help="where the samples, the outputs and the caches go (default: OUT/work)")
    ap.add_argument("--keep", action="store_true", help="keep the outputs")
    ap.add_argument("--skip-genome", action="store_true", help="leave out the whole HG002 (11 GB of output from each build)")
    ap.add_argument("--modes", default="vcf,tab", help="comma-separated: vcf, tab (default both; the whole HG002 is run as VCF only)")
    args = ap.parse_args()
    work = args.work or args.out / "work"
    args.out.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    for need in (GFF3, FASTA, HG002, CLINVAR, HERE / "human" / "vep_example_GRCh38.vcf", HERE / "human" / "clinvar_inframe_deletions.vcf"):
        if not need.exists():
            sys.exit(f"missing {need}")

    inputs = {
        "vep_example": HERE / "human" / "vep_example_GRCh38.vcf",
        "inframe_deletions": HERE / "human" / "clinvar_inframe_deletions.vcf",
    }
    clinvar_sample, clinvar_fifth, hg002_sample = work / "clinvar_sample.vcf", work / "clinvar_1in5.vcf", work / "hg002_sample.vcf"
    subprocess.run([sys.executable, str(HERE / "sample_variants.py"), str(CLINVAR), str(clinvar_sample), "--stratified"], check=True, capture_output=True)
    subprocess.run([sys.executable, str(HERE / "sample_variants.py"), str(CLINVAR), str(clinvar_fifth), "--every", "5"], check=True, capture_output=True)
    subprocess.run([sys.executable, str(HERE / "sample_variants.py"), str(HG002), str(hg002_sample), "--every", "200"], check=True, capture_output=True)
    inputs["clinvar_sample"], inputs["clinvar_1in5"], inputs["hg002_sample"] = clinvar_sample, clinvar_fifth, hg002_sample
    if not args.skip_genome:
        inputs["hg002_full"] = HG002

    builds = {args.base_label: args.base, args.new_label: args.new}
    models = {label: gene_model_dir(work, label) for label in builds}
    modes = [m for m in args.modes.split(",") if m]
    rows, differ = [], 0
    for mode in modes:
        for name, vcf in inputs.items():
            if name == "hg002_full" and mode != "vcf":
                continue  # the whole genome as a table is 10 GB from each build, and the VCF run reaches the same predictions
            n = records(vcf)
            digest, size = {}, {}
            for label, binary in builds.items():
                gff3, fasta = models[label]
                out = work / f"{name}.{mode}.{label}.out"
                wall = annotate(binary, gff3, fasta, vcf, out, MODES[mode])
                digest[label], size[label] = md5_of(out), out.stat().st_size
                print(f"{mode} {name:18s} {label:16s} {n:>9,} records  {size[label]:>14,} bytes  {digest[label]}  ({wall:.1f} s)", flush=True)
                if not args.keep:
                    out.unlink()
            same = len(set(digest.values())) == 1
            differ += not same
            rows.append((mode, name, n, digest[args.base_label], digest[args.new_label], size[args.base_label], size[args.new_label], "yes" if same else "NO"))

    with (args.out / "identity.tsv").open("w") as f:
        f.write(f"mode\tinput\trecords\tmd5_{args.base_label}\tmd5_{args.new_label}\tbytes_{args.base_label}\tbytes_{args.new_label}\tidentical\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    with (args.out / "environment.tsv").open("w") as f:
        f.write(f"started_utc\t{started}\n")
        for mode in modes:
            f.write(f"flags_{mode}\t{' '.join(MODES[mode])}\n")
        for label, binary in builds.items():
            version = subprocess.run([str(binary), "--version"], capture_output=True, text=True).stdout.strip()
            f.write(f"binary_{label}\t{binary.name}\t{version}\tsha256 {sha256_of(binary)}\n")
        f.write(f"gff3\t{GFF3.name}\nfasta\t{FASTA.name}\nclinvar_2star\t{CLINVAR.name}\n")
    print(f"{len(rows) - differ} of {len(rows)} runs are byte-identical between {args.base_label} and {args.new_label}")
    sys.exit(1 if differ else 0)


if __name__ == "__main__":
    main()
