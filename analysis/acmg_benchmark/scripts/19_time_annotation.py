#!/usr/bin/env python3
"""How long the annotation step of the benchmark takes: `annotate --acmg --pick` over the 673,660 variants of the ClinVar 2-star+ set, and the scoring of it.

    19_time_annotation.py --bin $FASTVEP_V38_BIN --out analysis/acmg_benchmark/results/v38_v0.5.0

`run_benchmark.sh` annotates the set with `--acmg --pick --output-format vcf`, writes the VCF to stdout and compresses it with `bgzip`; this
runs that command and times it from outside, on every core (the default) in rounds and on one thread once, and records the peak resident
set size of the `fastvep` process (`wait4`, the child that was waited for). The scoring (`03_evaluate_concordance.py`, which `run_benchmark.sh`
runs next) is timed once, on the output of the last run on every core. The documents of this directory used to say "about 2 minutes" and
"4,103 s (68 min)", the second from the first version of the classifier; this is the number for the release. Nothing else should run
meanwhile: the run takes the whole machine.

Writes annotation_time.tsv (every run, and the scoring) into --out. Standard library only; Linux and macOS (`os.wait4`).
"""

from __future__ import annotations

import argparse
import os
import tempfile
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BENCH = ROOT / "data" / "benchmark"
HUMAN = ROOT / "test_data" / "organisms" / "human"


def one_run(binary: Path, threads: str, dst: Path) -> tuple[float, float, int]:
    env = dict(os.environ)
    if threads != "all":
        env["RAYON_NUM_THREADS"] = threads
    cmd = [str(binary), "annotate", "-i", str(BENCH / "clinvar_2star.vcf"), "-o", "-", "--gff3", str(HUMAN / "Homo_sapiens.GRCh38.115.gff3"),
           "--fasta", str(HUMAN / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"), "--sa-dir", str(BENCH / "sa_db"), "--acmg", "--pick",
           "--output-format", "vcf"]
    t0 = time.perf_counter()
    with dst.open("wb") as out:
        fv = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, env=env)
        bg = subprocess.Popen(["bgzip", "-c"], stdin=fv.stdout, stdout=out)
        fv.stdout.close()
        _, status, usage = os.wait4(fv.pid, 0)
        fv.returncode = os.waitstatus_to_exitcode(status)
        bg.wait()
    wall = time.perf_counter() - t0
    if fv.returncode != 0 or bg.returncode != 0:
        sys.exit(f"FAILED: fastvep exit {fv.returncode}, bgzip exit {bg.returncode}")
    rss = usage.ru_maxrss / (1 << 20) if platform.system() == "Darwin" else usage.ru_maxrss / 1024  # bytes on macOS, KiB on Linux
    return wall, usage.ru_utime + usage.ru_stime, int(rss)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--runs", type=int, default=3, help="timed runs on every core (default 3, after one that is not timed)")
    ap.add_argument("--skip-one-thread", action="store_true")
    ap.add_argument("--work", type=Path, default=Path("/tmp"))
    args = ap.parse_args()
    for need in (args.bin, BENCH / "clinvar_2star.vcf", BENCH / "sa_db", HUMAN / "Homo_sapiens.GRCh38.115.gff3"):
        if not need.exists():
            sys.exit(f"missing {need}")
    args.out.mkdir(parents=True, exist_ok=True)
    dst = args.work / "annotation_time.vcf.gz"
    rows = []
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    one_run(args.bin, "all", dst)  # reads every input once
    for r in range(1, args.runs + 1):
        wall, cpu, rss = one_run(args.bin, "all", dst)
        rows.append(("all", r, wall, cpu, rss))
        print(f"every core, run {r}: {wall:.1f} s wall, {cpu:.0f} s CPU, {rss} MiB", flush=True)
    if not args.skip_one_thread:
        wall, cpu, rss = one_run(args.bin, "1", dst)
        rows.append(("1", 1, wall, cpu, rss))
        print(f"one thread: {wall:.1f} s wall, {cpu:.0f} s CPU, {rss} MiB", flush=True)
    # the scoring of the output of the last run on every core, into a directory that is removed afterwards
    last_all = dst  # the last run was on one thread unless skipped; its output is the same bytes either way
    t0 = time.perf_counter()
    with tempfile.TemporaryDirectory(dir=args.work) as scored:
        subprocess.run([sys.executable, str(Path(__file__).resolve().parent / "03_evaluate_concordance.py"), "--labels", str(BENCH / "clinvar_2star_labels.tsv"),
                        "--predictions", str(last_all), "--out", scored], check=True, capture_output=True)
    score_wall = time.perf_counter() - t0
    print(f"scoring: {score_wall:.1f} s", flush=True)
    rows.append(("score", 1, score_wall, 0.0, 0))
    dst.unlink()
    with (args.out / "annotation_time.tsv").open("w") as f:
        f.write(f"# {args.bin} ({subprocess.run([str(args.bin), '--version'], capture_output=True, text=True).stdout.strip()}), started {started}, {platform.platform()}\n")
        f.write("threads\trun\twall_s\tcpu_s\tpeak_rss_mib\n")
        for t, r, w, c, m in rows:
            f.write(f"{t}\t{r}\t{w:.2f}\t{c:.1f}\t{m}\n")
    walls = [w for t, _, w, _, _ in rows if t == "all"]
    print(f"every core: median {statistics.median(walls):.1f} s over {len(walls)} runs")


if __name__ == "__main__":
    main()
