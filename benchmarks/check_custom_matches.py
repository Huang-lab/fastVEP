#!/usr/bin/env python3
"""Do fastVEP and Ensembl VEP find the same records in the same ClinVar and gnomAD files?

The annotation conditions of the head-to-head (`run_head_to_head.py`, `run_vep_side.py`) give both tools the same source data: fastVEP as
fastSA databases built from the files, Ensembl VEP as `--custom` sources. A condition is only a comparison of speed if the two do the same work,
and a `--custom` source that VEP cannot match (a chromosome name, an allele it will not equate) is quicker than one that works. So this runs both
tools over one chromosome of HG002 with the joint condition and compares, record by record, whether each found a ClinVar record and a gnomAD record.

    benchmarks/check_custom_matches.py --bin target/release/fastvep --chrom 22 --out benchmarks/results/v0.5.0/with_vep

A record "has" a source for VEP when the column that VEP names after the source (`ClinVar`, `gnomAD`) is filled in any of its CSQ entries, and for
fastVEP when a transcript consequence carries the `clinvar` or `gnomad` object. (Not a field of the source: a ClinVar record can have no CLNSIG, and a
gnomAD record that failed its filters has no AF_joint, and both are matches.) Records are compared by their place in the input, which both tools keep.

Writes custom_matches_chr<N>.tsv (what each tool found, and the records on which they differ) into --out and exits 1 if there is a difference.
Needs Docker, the files that `run_vep_side.py` needs, and the database directories that `run_head_to_head.py` builds.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_vep_side as vs  # noqa: E402  (the VEP side of the harness: its container call and its --custom arguments)


def vep_matches(chrom: str, work: Path) -> tuple[list[bool], list[bool]]:
    custom = vs.prepare_custom(vs.SOURCES / "clinvar.vcf.gz", vs.SOURCES / "gnomad_joint_extracts", [chrom], work / "custom")
    gff = vs.sorted_gff3(vs.HUMAN / "Homo_sapiens.GRCh38.115.gff3")
    fasta = vs.HUMAN / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    vcf = vs.INPUTS / "hg002_by_chr" / f"chr{chrom}.vcf"
    _, _, _, out = vs.docker_vep(vcf, gff, fasta, f"check_chr{chrom}.vcf", work, custom=custom, extra=vs.custom_args("joint", chrom))
    clinvar: list[bool] = []
    gnomad: list[bool] = []
    idx = None
    with out.open() as f:
        for line in f:
            if line.startswith("#"):
                m = re.search(r'ID=CSQ,.*Format: ([^"]+)"', line)
                if m:
                    names = m.group(1).split("|")
                    idx = (names.index("ClinVar"), names.index("gnomAD"))
                continue
            info = line.split("\t", 8)[7]
            entries = [e.split("|") for e in info.split("CSQ=", 1)[1].split(";", 1)[0].split(",")]
            clinvar.append(any(len(e) > idx[0] and e[idx[0]] for e in entries))
            gnomad.append(any(len(e) > idx[1] and e[idx[1]] for e in entries))
    out.unlink()
    return clinvar, gnomad


def fastvep_matches(binary: Path, chrom: str, sa_dir: Path, work: Path) -> tuple[list[bool], list[bool]]:
    vcf = vs.INPUTS / "hg002_by_chr" / f"chr{chrom}.vcf"
    out = work / f"check_chr{chrom}.json"
    subprocess.run([str(binary), "annotate", "-i", str(vcf), "-o", str(out), "--output-format", "json", "--gff3", str(vs.HUMAN / "Homo_sapiens.GRCh38.115.gff3"),
                    "--fasta", str(vs.HUMAN / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"), "--hgvs", "--symbol", "--canonical", "--no-progress",
                    "--sa-dir", str(sa_dir)], check=True, capture_output=True)
    records = json.loads(out.read_text())
    out.unlink()
    has = lambda r, key: any(key in t for t in r.get("transcript_consequences", []))  # noqa: E731
    return [has(r, "clinvar") for r in records], [has(r, "gnomad") for r in records]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", required=True, type=Path)
    ap.add_argument("--chrom", default="22")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--sa-dir", type=Path, default=vs.INPUTS / "sa_joint", help="the directory run_head_to_head.py makes for the joint condition")
    args = ap.parse_args()
    for need in (args.sa_dir, vs.INPUTS / "hg002_by_chr" / f"chr{args.chrom}.vcf"):
        if not need.exists():
            sys.exit(f"missing {need}: run benchmarks/run_head_to_head.py (--only wgs_joint --list) and benchmarks/run_vep_side.py first")
    work = vs.ROOT / "benchmarks" / "output" / "vep_side"
    work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    v_cv, v_gn = vep_matches(args.chrom, work)
    f_cv, f_gn = fastvep_matches(args.bin.resolve(), args.chrom, args.sa_dir, work)
    if not len(v_cv) == len(f_cv):
        sys.exit(f"VEP wrote {len(v_cv):,} records and fastVEP {len(f_cv):,}")
    rows = [["source", "records", "vep_found", "fastvep_found", "both", "only_vep", "only_fastvep"]]
    bad = 0
    for name, v, f in (("ClinVar", v_cv, f_cv), ("gnomAD", v_gn, f_gn)):
        both = sum(1 for a, b in zip(v, f) if a and b)
        only_v = [i for i, (a, b) in enumerate(zip(v, f)) if a and not b]
        only_f = [i for i, (a, b) in enumerate(zip(v, f)) if b and not a]
        rows.append([name, str(len(v)), str(sum(v)), str(sum(f)), str(both), str(len(only_v)), str(len(only_f))])
        bad += len(only_v) + len(only_f)
        print(f"{name:8s} {len(v):>8,} records: VEP found {sum(v):,}, fastVEP found {sum(f):,}, both {both:,}, only VEP {len(only_v)}, only fastVEP {len(only_f)}")
    (args.out / f"custom_matches_chr{args.chrom}.tsv").write_text("".join("\t".join(r) + "\n" for r in rows))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
