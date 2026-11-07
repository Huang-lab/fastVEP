#!/usr/bin/env python3
"""fastSA v1 against v2 on real databases: size, build time and annotation output (Table 6 of the manuscript, the table of the README).

    benchmarks/run_sa_formats.py --bin target/release/fastvep --out benchmarks/results/v0.5.0-rc.1

For each source it builds the v1 (`.osa` and its `.osa.idx`) and the v2 (`.osa2`) database from the same file with `fastvep sa-build --format`,
and records the size on disk and the wall time of the build. ClinVar is the whole database; REVEL is chromosomes 1 and 22, cut out of the
REVEL zip in one pass (`sa-build` reads a plain or gzipped CSV, not a zip). It then annotates the 673,660 variants of the ClinVar 2-star+ set
against a directory that holds the v1 databases and against one that holds the v2, first with ClinVar alone and then with ClinVar and the two
REVEL chromosomes, and compares the two outputs of each pair byte for byte: a format that changed an answer would be a storage optimization
that is not one. The two formats alternate within a round of runs, and the first round, which is not timed, counts the records that carry a
REVEL score, so that an output that is identical because a database was never read cannot pass. Two bugs of the v2 writer were found after the
measurement of July 2026 (it dropped the ClinVar records with N or IUPAC bases, and resolved duplicate keys arbitrarily, #80), so the
byte-for-byte claim is re-checked here and not carried over.

Needs data/benchmark/sa_sources/clinvar.vcf.gz and revel-v1.3_all_chromosomes.zip, data/benchmark/clinvar_2star.vcf, and the human GFF3 and FASTA
of the speed benchmark. Writes sa_formats.csv, sa_formats_annotation.csv and sa_formats_environment.tsv into --out. Standard library only.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SOURCES = ROOT / "data" / "benchmark" / "sa_sources"
CLINVAR_2STAR = ROOT / "data" / "benchmark" / "clinvar_2star.vcf"
HUMAN = ROOT / "test_data" / "organisms" / "human"
WORK = ROOT / "benchmarks" / "output" / "sa_formats"


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def run(cmd: list[str]) -> tuple[float, str]:
    t0 = time.perf_counter()
    r = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    if r.returncode != 0:
        sys.exit(f"FAILED (exit {r.returncode}): {' '.join(cmd)}\n{r.stderr[-800:]}")
    return wall, r.stdout + r.stderr


def records_of(text: str) -> int | None:
    found = [int(m.replace(",", "")) for m in re.findall(r"([\d,]{4,}) (?:records|variants|entries|rows)", text)]
    return max(found) if found else None


def build(binary: Path, source: str, src: Path, base: Path, fmt: str) -> dict:
    for old in base.parent.glob(base.name + ".*"):
        old.unlink()
    wall, out = run([str(binary), "sa-build", "--source", source, "-i", str(src), "-o", str(base), "--assembly", "GRCh38", "--format", fmt, "--no-progress"])
    files = sorted(base.parent.glob(base.name + ".*"))
    size = sum(f.stat().st_size for f in files)
    return {"wall": wall, "bytes": size, "files": [f.name for f in files], "records": records_of(out)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--runs", type=int, default=3, help="annotation runs per format, the median is reported (default 3)")
    args = ap.parse_args()
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    for need in (SOURCES / "clinvar.vcf.gz", SOURCES / "revel-v1.3_all_chromosomes.zip", CLINVAR_2STAR, HUMAN / "Homo_sapiens.GRCh38.115.gff3"):
        if not need.exists():
            sys.exit(f"missing {need}")
    args.out.mkdir(parents=True, exist_ok=True)
    (WORK / "v1").mkdir(parents=True, exist_ok=True)
    (WORK / "v2").mkdir(parents=True, exist_ok=True)

    results = []
    print("ClinVar ...", flush=True)
    cv = {fmt: build(args.bin, "clinvar", SOURCES / "clinvar.vcf.gz", WORK / ("v1" if fmt == "osa" else "v2") / "clinvar", fmt) for fmt in ("osa", "osa2")}
    results.append(("ClinVar (whole database)", cv))

    print("REVEL chromosomes 1 and 22 (one pass over the zip) ...", flush=True)
    awk = ("NR==1 {h=$0; print h > c1; print h > c22} $1==1 {print > c1} $1==22 {print > c22}")
    p1, p22 = WORK / "revel_chr1.csv", WORK / "revel_chr22.csv"
    subprocess.run(f"unzip -p '{SOURCES / 'revel-v1.3_all_chromosomes.zip'}' | awk -F, -v c1='{p1}' -v c22='{p22}' '{awk}'", shell=True, check=True)
    for chrom, csv_path in (("chr1", p1), ("chr22", p22)):
        print(f"REVEL {chrom} ...", flush=True)
        r = {fmt: build(args.bin, "revel", csv_path, WORK / ("v1" if fmt == "osa" else "v2") / f"revel_{chrom}", fmt) for fmt in ("osa", "osa2")}
        results.append((f"REVEL ({chrom})", r))
        csv_path.unlink()

    with (args.out / "sa_formats.csv").open("w") as f:
        f.write("source,records,v1_bytes,v2_bytes,ratio,build_v1_s,build_v2_s\n")
        for name, r in results:
            rec = r["osa"]["records"] or r["osa2"]["records"] or ""
            f.write(f"{name},{rec},{r['osa']['bytes']},{r['osa2']['bytes']},{r['osa2']['bytes'] / r['osa']['bytes']:.2f},{r['osa']['wall']:.1f},{r['osa2']['wall']:.1f}\n")

    print("annotating the ClinVar 2-star+ set against the v1 and the v2 databases ...", flush=True)
    configs = [("clinvar", ["clinvar.*"]), ("clinvar+revel", ["clinvar.*", "revel_chr1.*", "revel_chr22.*"])]
    gff3, fasta = HUMAN / "Homo_sapiens.GRCh38.115.gff3", HUMAN / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    nvar = sum(1 for line in CLINVAR_2STAR.open() if not line.startswith("#"))
    rows, same = [], True
    for name, globs in configs:
        dirs = {}
        for fmt, d in (("osa", "v1"), ("osa2", "v2")):
            sa_dir = WORK / f"sa_{name.replace('+', '_')}_{d}"
            sa_dir.mkdir(exist_ok=True)
            for old in sa_dir.iterdir():
                old.unlink()
            for g in globs:
                for f in (WORK / d).glob(g):
                    (sa_dir / f.name).symlink_to(f)
            dirs[fmt] = sa_dir
        times = {"osa": [], "osa2": []}
        digests = {"osa": set(), "osa2": set()}
        with_revel = {}
        # the two formats alternate within a round, so that what the machine does between two rounds falls on both
        for r in range(args.runs + 1):
            for fmt, d in (("osa", "v1"), ("osa2", "v2")):
                out = WORK / f"annotated_{d}.vcf"
                wall, _ = run([str(args.bin), "annotate", "-i", str(CLINVAR_2STAR), "--gff3", str(gff3), "--fasta", str(fasta),
                               "--sa-dir", str(dirs[fmt]), "--no-progress", "-o", str(out)])
                if r > 0:  # the first run reads every input once
                    times[fmt].append(wall)
                else:
                    # the databases must have been read: count the records that carry a REVEL score
                    with out.open("rb") as fh:
                        with_revel[fmt] = sum(1 for line in fh if b"REVEL" in line)
                digests[fmt].add(md5(out))
                out.unlink()
        for fmt in ("osa", "osa2"):
            one = len(digests[fmt]) == 1
            rows.append((name, "+".join(g.split(".")[0] for g in globs), fmt, nvar, statistics.median(times[fmt]),
                         ";".join(f"{t:.2f}" for t in times[fmt]), next(iter(digests[fmt])) if one else "VARIES", with_revel[fmt]))
        same = same and len(digests["osa"] | digests["osa2"]) == 1
    with (args.out / "sa_formats_annotation.csv").open("w") as f:
        f.write("config,databases,format,variants,median_s,runs_s,md5,records_with_revel\n")
        for name, dbs, fmt, n, med, runs, digest, revel in rows:
            f.write(f"{name},{dbs},{fmt},{n},{med:.2f},{runs},{digest},{revel}\n")
    with (args.out / "sa_formats_environment.tsv").open("w") as f:
        f.write(f"binary\t{args.bin}\nbinary_version\t{subprocess.run([str(args.bin), '--version'], capture_output=True, text=True).stdout.strip()}\n")
        f.write(f"annotation_runs\t{args.runs}\nidentical_output\t{'yes' if same else 'NO'}\nstarted_utc\t{started}\n")
    print(f"identical output from the v1 and the v2 databases: {'yes' if same else 'NO'}")
    if not same:
        sys.exit("the v1 and the v2 databases gave different annotations")

if __name__ == "__main__":
    main()
