#!/usr/bin/env python3
"""Does cutting the benchmark's inputs to five columns change what fastVEP's time says?

    benchmarks/check_lean_inputs.py --bin target/release/fastvep --out benchmarks/results/v0.5.0

`run_head_to_head.py` cuts the chromosome 22 and chromosome 1 inputs of Tables 2, 4 and 5 from the GIAB HG002 VCF to CHROM POS ID REF ALT,
with `.` for the rest, because neither tool reads the INFO of that file to predict a consequence. Ensembl VEP reads the same cut text
(`run_vep_side.py`), so the tables compare like with like; this measures what the cut is worth to fastVEP: all
50,284 chromosome 22 records, cut and as they are, on one thread, with the flags of the head-to-head, the two alternating within a round, one
round that is not reported and then `--runs` that are. It writes lean_inputs.tsv (every run) into --out and prints the medians.
Standard library only.
"""

from __future__ import annotations

import argparse
import statistics
import subprocess
import sys
import time
import os
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_head_to_head as h  # noqa: E402


def as_is(src: Path, dst: Path, chrom: str) -> int:
    n = 0
    with src.open() as fin, dst.open("w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            if line.split("\t", 1)[0] == chrom:
                fout.write(line)
                n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--runs", type=int, default=5)
    args = ap.parse_args()
    d = h.WORK / "inputs"
    d.mkdir(parents=True, exist_ok=True)
    lean, asis = d / "chr22_all.vcf", d / "chr22_all_as_is.vcf"
    if not lean.exists():
        h.lean_vcf(h.HG002, lean, "22", None)
    n = as_is(h.HG002, asis, "22")
    cells = {"lean": lean, "as_is": asis}
    env = {**os.environ, "RAYON_NUM_THREADS": "1"}
    walls = {k: [] for k in cells}
    sizes = {}
    out = h.WORK / "lean_inputs_out.vcf"
    for r in range(args.runs + 1):
        for k, vcf in cells.items():
            t0 = time.perf_counter()
            p = subprocess.run([str(args.bin), "annotate", "-i", str(vcf), "--gff3", str(h.GFF3), "--fasta", str(h.FASTA),
                                "--hgvs", "--symbol", "--canonical", "--no-progress", "-o", str(out)], env=env, capture_output=True, text=True)
            wall = time.perf_counter() - t0
            if p.returncode != 0:
                sys.exit(f"FAILED: {p.stderr[-800:]}")
            sizes[k] = out.stat().st_size
            out.unlink()
            if r > 0:  # the first round reads every input once
                walls[k].append(wall)
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "lean_inputs.tsv").open("w") as f:
        f.write("cell\trecords\tinput_bytes\tout_bytes\tmedian_s\truns_s\n")
        for k, vcf in cells.items():
            f.write(f"{k}\t{n}\t{vcf.stat().st_size}\t{sizes[k]}\t{statistics.median(walls[k]):.2f}\t{';'.join(f'{w:.2f}' for w in walls[k])}\n")
    m = {k: statistics.median(v) for k, v in walls.items()}
    print(f"chromosome 22, {n:,} records, one thread: cut to five columns {m['lean']:.2f} s, as they are {m['as_is']:.2f} s ({(m['as_is'] / m['lean'] - 1) * 100:+.1f}%)")


if __name__ == "__main__":
    main()
