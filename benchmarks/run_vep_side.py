#!/usr/bin/env python3
"""The Ensembl VEP side of the small head-to-head cells: the chromosome 22 series, start-up and peak memory.

`run_head_to_head.py` measures fastVEP; this measures Ensembl VEP on the same inputs, so that a speed-up in a table is a ratio of two
times taken on one day on one machine. It does not repeat the whole-genome runs: VEP runs out of the Docker default memory above about
2.7 million variants, so its whole-genome time is the sum of 22 per-chromosome runs, each the minimum of two, about 2.7 hours for one
annotation condition (`manuscript/data/vep_baseline_2026-07-29.csv` has the sums of 29 July 2026).

    benchmarks/run_vep_side.py --out benchmarks/results/v0.5.0-rc.1
    benchmarks/run_vep_side.py --out ... --only scaling --runs 3

What it measures, with Ensembl VEP 115.1 in Docker, `--gff` mode against a bgzip and tabix copy of the same GRCh38 Ensembl 115 GFF3 (VEP
rejects an unindexed one), `--vcf --hgvs --symbol --canonical --no_stats`, one thread:

    scaling   the first 1,000, 5,000 and 10,000 chr22 variants of HG002, and all 50,284 (the inputs `run_head_to_head.py` cuts)
    startup   one variant, against the human and the yeast gene model
    memory    peak resident memory (VmHWM of the VEP process, read from /proc inside the container) for 100,000 chr1 variants

The time is the wall-clock time of `docker run`, which includes starting the container: that is what a user of the Docker image waits for,
and the only clean host-side metric (`/usr/bin/time` on `docker run` times the client). The container is started with the options that
`validation/run_validation.sh` uses, `VEP_DATABASE=0` and no `--offline` (see the comment there).

Writes vep_raw.tsv, vep_summary.tsv and vep_resources.csv (the three cells that Table 4 of the manuscript quotes) into --out, and a
`vep_environment.tsv` with the image, Docker's CPUs and memory, and the load. Standard library only; needs Docker, bgzip and tabix.
"""

from __future__ import annotations

import argparse
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HUMAN = ROOT / "test_data" / "organisms" / "human"
YEAST = ROOT / "test_data" / "organisms" / "yeast"
INPUTS = ROOT / "benchmarks" / "output" / "head_to_head" / "inputs"
IMAGE = "ensemblorg/ensembl-vep:release_115.1"
FLAGS = ["--vcf", "--force_overwrite", "--no_stats", "--hgvs", "--symbol", "--canonical"]
PARTS = ("scaling", "startup", "memory")


def sh(*cmd: str, check: bool = True) -> str:
    return subprocess.run(cmd, capture_output=True, text=True, check=check).stdout.strip()


def sorted_gff3(gff3: Path) -> Path:
    """A bgzip + tabix copy, made once next to the GFF3 (VEP reads it through tabix)."""
    out = gff3.with_name(gff3.stem + ".sorted.gff3.gz")
    if not (out.exists() and out.with_name(out.name + ".tbi").exists()):
        print(f"preparing {out.name} (once) ...", flush=True)
        body = subprocess.run(f"grep -v '^#' '{gff3}' | sort -k1,1 -k4,4n -t$'\\t' | bgzip -@ 4 -c > '{out}'", shell=True, executable="/bin/bash")
        if body.returncode != 0:
            sys.exit(f"could not make {out}")
        sh("tabix", "-p", "gff", str(out))
    return out


def docker_vep(vcf: Path, gff: Path, fasta: Path, out_name: str, work: Path, peak: bool = False) -> tuple[float, int | None]:
    """One run of VEP in a fresh container: (wall seconds, peak resident kB or None)."""
    # test_data may be a link to another checkout (a git worktree): resolve every path VEP is given, so that the mount has it.
    vep = ["vep", "--input_file", f"/in/{vcf.name}", "--gff", f"/ref/{gff.name}", "--fasta", f"/ref/{fasta.name}", "--output_file", f"/out/{out_name}", *FLAGS]
    mounts = ["-v", f"{vcf.resolve().parent}:/in:ro", "-v", f"{gff.resolve().parent}:/ref:ro", "-v", f"{work}:/out"]
    if fasta.resolve().parent != gff.resolve().parent:
        mounts += ["-v", f"{fasta.resolve().parent}:/ref2:ro"]
        vep[vep.index("--fasta") + 1] = f"/ref2/{fasta.name}"
    cmd = ["docker", "run", "--rm", "-e", "VEP_DATABASE=0", *mounts, IMAGE]
    if peak:
        # VmHWM is the kernel's high-water mark of the process's resident set; the last value read before it exits is its peak.
        script = ('"$@" >/dev/null 2>&1 & p=$!; m=0; while kill -0 $p 2>/dev/null; do '
                  'v=$(awk \'/VmHWM/ {print $2}\' /proc/$p/status 2>/dev/null); [ -n "$v" ] && m=$v; sleep 0.1; done; wait $p; echo "VmHWM_kB=$m"')
        cmd += ["sh", "-c", script, "sh", *vep]
    else:
        cmd += vep
    t0 = time.perf_counter()
    r = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    if r.returncode != 0:
        sys.exit(f"FAILED (exit {r.returncode}): {' '.join(cmd)}\n{r.stderr[-800:]}")
    out_file = work / out_name
    if not out_file.exists() or out_file.stat().st_size == 0:
        sys.exit(f"VEP wrote nothing to {out_file}: {' '.join(cmd)}\n{r.stdout[-400:]}\n{r.stderr[-400:]}")
    kb = None
    if peak:
        for line in r.stdout.splitlines():
            if line.startswith("VmHWM_kB="):
                kb = int(line.split("=")[1])
    out_file.unlink()
    return wall, kb


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--runs", type=int, default=3, help="timed runs per cell; the median is reported (default 3, as for the VEP times of the manuscript)")
    ap.add_argument("--only", default="", help=f"comma-separated parts: {', '.join(PARTS)}")
    args = ap.parse_args()
    parts = [p for p in (args.only.split(",") if args.only else PARTS) if p]
    for p in parts:
        if p not in PARTS:
            sys.exit(f"unknown part {p!r}")
    for need in (INPUTS / "chr22_all.vcf", INPUTS / "one.vcf", INPUTS / "chr1_100k.vcf", INPUTS / "yeast_one.vcf"):
        if not need.exists():
            sys.exit(f"missing {need}: run benchmarks/run_head_to_head.py first, which cuts the inputs")
    info = sh("docker", "info", "--format", "{{.ServerVersion}} cpus={{.NCPU}} mem={{.MemTotal}}", check=False)
    if not info or "cpus=0" in info:
        sys.exit("Docker's engine is not answering")
    args.out.mkdir(parents=True, exist_ok=True)
    work = ROOT / "benchmarks" / "output" / "vep_side"
    work.mkdir(parents=True, exist_ok=True)
    human_gff, yeast_gff = sorted_gff3(HUMAN / "Homo_sapiens.GRCh38.115.gff3"), sorted_gff3(YEAST / "yeast.gff3")
    human_fa, yeast_fa = HUMAN / "Homo_sapiens.GRCh38.dna.primary_assembly.fa", YEAST / "yeast.fa"
    image_id = sh("docker", "image", "inspect", IMAGE, "--format", "{{.Id}}", check=False)

    cells = []  # (part, name, vcf, gff, fasta, peak)
    if "scaling" in parts:
        cells += [("scaling", f"chr22_{n}", INPUTS / f"chr22_{n}.vcf", human_gff, human_fa, False) for n in ("1000", "5000", "10000", "all")]
    if "startup" in parts:
        cells += [("startup", "human_1_variant", INPUTS / "one.vcf", human_gff, human_fa, False),
                  ("startup", "yeast_1_variant", INPUTS / "yeast_one.vcf", yeast_gff, yeast_fa, False)]
    if "memory" in parts:
        cells += [("memory", "chr1_100k", INPUTS / "chr1_100k.vcf", human_gff, human_fa, True)]

    raw: list[list[str]] = [["part", "cell", "round", "wall_s", "peak_rss_mib", "load1"]]
    print(f"Docker {info}; image {image_id[:19]}; {args.runs} timed runs per cell after one that is not timed", flush=True)
    for part, name, vcf, gff, fasta, peak in cells:
        for r in range(0, args.runs + 1):
            if r == 0 and peak:
                continue  # the 100,000-variant run takes about two minutes; its warm-up is the run itself
            wall, kb = docker_vep(vcf, gff, fasta, f"vep_{name}.vcf", work, peak)
            raw.append([part, name, str(r), f"{wall:.2f}", f"{kb / 1024:.0f}" if kb else "", f"{os.getloadavg()[0]:.2f}"])
            print(f"  {part:8s} {name:18s} round {r}: {wall:8.2f} s" + (f"  {kb / 1024:.0f} MiB" if kb else ""), flush=True)

    with (args.out / "vep_raw.tsv").open("w") as f:
        f.writelines("\t".join(row) + "\n" for row in raw)
    summary = [["part", "cell", "n", "median_wall_s", "min_wall_s", "max_wall_s", "median_peak_rss_mib"]]
    names = list(dict.fromkeys((r[0], r[1]) for r in raw[1:]))
    med: dict[str, float] = {}
    for part, name in names:
        rows = [r for r in raw[1:] if r[0] == part and r[1] == name and (int(r[2]) >= 1 or part == "memory")]
        w = [float(r[3]) for r in rows]
        peaks = [float(r[4]) for r in rows if r[4]]
        summary.append([part, name, str(len(w)), f"{statistics.median(w):.2f}", f"{min(w):.2f}", f"{max(w):.2f}", f"{statistics.median(peaks):.0f}" if peaks else ""])
        med[f"{part}:{name}"] = statistics.median(w)
        if peaks:
            med[f"peak:{name}"] = statistics.median(peaks)
    with (args.out / "vep_summary.tsv").open("w") as f:
        f.writelines("\t".join(row) + "\n" for row in summary)
    with (args.out / "vep_environment.tsv").open("w") as f:
        for k, v in (("image", IMAGE), ("image_id", image_id), ("docker", info), ("flags", " ".join(FLAGS)), ("started_utc", time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())),
                     ("runs", str(args.runs)), ("load_end", " ".join(f"{x:.2f}" for x in os.getloadavg()))):
            f.write(f"{k}\t{v}\n")
    resources = [["metric", "value"]]
    if "startup:yeast_1_variant" in med:
        resources.append(["startup_yeast", f"{med['startup:yeast_1_variant']:.1f} sec"])
        resources.append(["startup_human", f"{med['startup:human_1_variant']:.1f} sec"])
    if "peak:chr1_100k" in med:
        resources.append(["peak_memory", f"{med['peak:chr1_100k'] / 1024:.1f} GiB"])
    if len(resources) > 1:
        with (args.out / "vep_resources.csv").open("w") as f:
            f.writelines(",".join(row) + "\n" for row in resources)
    print(f"\nwrote vep_raw.tsv, vep_summary.tsv, vep_environment.tsv" + (", vep_resources.csv" if len(resources) > 1 else "") + f" in {args.out}")


if __name__ == "__main__":
    main()
