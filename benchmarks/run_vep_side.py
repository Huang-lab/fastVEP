#!/usr/bin/env python3
"""The Ensembl VEP side of the head-to-head: the chromosome 22 series, start-up, peak memory and the whole genome.

`run_head_to_head.py` measures fastVEP; this measures Ensembl VEP on the same inputs on the same machine, so that a speed-up in a table is
a ratio of two times taken on one day.

    benchmarks/run_vep_side.py --out benchmarks/results/v0.5.0/with_vep
    benchmarks/run_vep_side.py --out ... --only scaling --runs 3
    benchmarks/run_vep_side.py --out ... --only genome              # resumes: a (condition, chromosome) already in vep_genome.tsv is skipped

What it measures, with Ensembl VEP 116.3 in Docker (the image is pinned by digest), `--gff` mode against a bgzip and tabix copy of the same
GRCh38 Ensembl 115 GFF3 that fastVEP reads (VEP rejects an unindexed one), `--vcf --hgvs --symbol --canonical --no_stats`, one thread:

    scaling   the first 1,000, 5,000 and 10,000 chr22 variants of HG002, and all 50,284 (the inputs `run_head_to_head.py` cuts)
    startup   one variant, against the human and the yeast gene model
    memory    peak resident memory (VmHWM of the VEP process, read from /proc inside the container) for 100,000 chr1 variants
    genome    the 4,048,342 variants of HG002, one run per chromosome (1 to 22, the chromosomes the file has) and per condition:
                plain    consequences only
                clinvar  and `--custom` ClinVar
                joint    and `--custom` ClinVar and gnomAD v4.1 joint, one file per chromosome
              The whole-genome time is the sum of the 22 runs. A single run over all the variants is not attempted: VEP 115.1 ran out of the
              7.65 GiB that Docker gives a container by default above about 2.7 million variants (29 July 2026), and one chromosome at a time
              is the way VEP completes, and bounds its memory, so the sum is its best case.

The `--custom` files are the very ones the fastSA databases of the same condition were built from (`data/benchmark/sa_sources/clinvar.vcf.gz`
for `clinvar.osa`, and the `gnomad_joint_extracts/` for `sa_db_joint/`), and each is asked for the INFO fields that the database keeps, so that
the two tools read the same source data and write the same number of annotation values. VEP matches the `chr`-prefixed names of the gnomAD
extracts to the unprefixed `1` of the VCF by itself (probed 2026-11-07 on 5,000 chr22 variants, 245 matches with and without the prefix).

The per-chromosome inputs are the sites-only text that the chr22 cells read (CHROM POS ID REF ALT, then three dots), which is VEP's best case:
fastVEP's whole-genome cell reads the 2.6 GB file with its genotype columns.

The time is the wall-clock time of `docker run`, which includes starting the container: that is what a user of the Docker image waits for,
and the only clean host-side metric (`/usr/bin/time` on `docker run` times the client). The CPU time of VEP itself is read inside the container
with the shell's `times`, so that a run that waited (another process on the machine) shows as wall far above CPU. The container is started with
the options that `validation/run_validation.sh` uses, `VEP_DATABASE=0` and no `--offline` (see the comment there).

Writes vep_raw.tsv, vep_summary.tsv and vep_resources.csv (the cells that Table 4 of the manuscript quotes), vep_genome.tsv, vep_genome_sums.csv
and vep_baseline.csv (the sums, in the layout of `manuscript/data/vep_baseline_*.csv`) into --out, and a `vep_environment.tsv` with the image, its
size (pulled and unpacked), Docker's CPUs and memory, the `--custom` arguments and the load. Standard library only; needs Docker, bgzip and tabix.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HUMAN = ROOT / "test_data" / "organisms" / "human"
YEAST = ROOT / "test_data" / "organisms" / "yeast"
HG002 = HUMAN / "human_giab_hg002_full.vcf"
INPUTS = ROOT / "benchmarks" / "output" / "head_to_head" / "inputs"
SOURCES = ROOT / "data" / "benchmark" / "sa_sources"
IMAGE = "ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3"  # 116.3
FLAGS = ["--vcf", "--force_overwrite", "--no_stats", "--hgvs", "--symbol", "--canonical"]
PARTS = ("scaling", "startup", "memory", "genome")
CONDITIONS = ("plain", "clinvar", "joint")
CHROMS = [str(c) for c in range(1, 23)]

# The INFO fields that fastVEP's databases keep (crates/fastvep-sa/src/sources/clinvar.rs and gnomad.rs) and the extract has.
CLINVAR_FIELDS = ["CLNSIG", "CLNREVSTAT", "CLNDN", "CLNVC", "CLNVCSO", "AF_EXAC", "AF_TGP", "AF_ESP"]
POPULATIONS = ["afr", "amr", "asj", "eas", "fin", "mid", "nfe", "remaining", "sas"]
GNOMAD_FIELDS = (["AC_joint", "AN_joint", "AF_joint", "nhomalt_joint", "AC_joint_XY", "AN_joint_XY"]
                 + [f"AF_joint_{p}" for p in POPULATIONS] + ["faf95_joint", "fafmax_faf95_max_joint"])
# The column of the CSQ that is filled when a record found a match: the one that VEP names after the source and fills with its ID or its
# coordinates. Not a field of the source: a ClinVar record can have no CLNSIG, and a gnomAD one that failed its filters has no AF_joint.
MARKERS = {"clinvar": ["ClinVar"], "joint": ["ClinVar", "gnomAD"], "plain": []}
SCALING_VARIANTS = {"chr22_1000": 1000, "chr22_5000": 5000, "chr22_10000": 10000, "chr22_all": 50284}


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


def indexed_copy(src: Path, dst: Path) -> None:
    """`src` linked (or copied, across file systems) to `dst`, with a tabix index beside it: the data directory is not written to."""
    if not dst.exists():
        try:
            os.link(src, dst)
        except OSError:
            shutil.copyfile(src, dst)
    if not dst.with_name(dst.name + ".tbi").exists():
        sh("tabix", "-p", "vcf", str(dst))


def prepare_custom(clinvar: Path, gnomad_dir: Path, chroms: list[str], out: Path) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    indexed_copy(clinvar, out / "clinvar.vcf.gz")
    for c in chroms:
        indexed_copy(gnomad_dir / f"gnomad_chr{c}.vcf.gz", out / f"gnomad_chr{c}.vcf.gz")
    return out


def custom_args(cond: str, chrom: str) -> list[str]:
    """The `--custom` arguments of a condition (paths as the container sees them)."""
    args = []
    if cond in ("clinvar", "joint"):
        args += ["--custom", f"file=/custom/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,fields={'%'.join(CLINVAR_FIELDS)}"]
    if cond == "joint":
        args += ["--custom", f"file=/custom/gnomad_chr{chrom}.vcf.gz,short_name=gnomAD,format=vcf,type=exact,fields={'%'.join(GNOMAD_FIELDS)}"]
    return args


def cut_chromosomes(dst: Path) -> dict[str, int]:
    """One sites-only VCF per chromosome of HG002, made once in one pass over the 2.6 GB file (the same text `lean_vcf` of the harness writes)."""
    counts_file = dst / "counts.tsv"
    if counts_file.exists():
        counts = {k: int(v) for k, v in (line.split("\t") for line in counts_file.read_text().splitlines())}
        if set(counts) == set(CHROMS) and all((dst / f"chr{c}.vcf").stat().st_size > 0 for c in CHROMS):
            return counts
    dst.mkdir(parents=True, exist_ok=True)
    print(f"cutting the chromosomes of {HG002.name} (once) ...", flush=True)
    counts = {c: 0 for c in CHROMS}
    outs = {c: (dst / f"chr{c}.vcf").open("w") for c in CHROMS}
    try:
        for f in outs.values():
            f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        with HG002.open() as fin:
            for line in fin:
                if line.startswith("#"):
                    continue
                f = line.split("\t", 5)
                if f[0] not in outs:
                    sys.exit(f"{HG002.name} has a record on {f[0]!r}, which is not one of chromosomes 1 to 22")
                outs[f[0]].write("\t".join(f[:5]) + "\t.\t.\t.\n")
                counts[f[0]] += 1
    finally:
        for f in outs.values():
            f.close()
    counts_file.write_text("".join(f"{c}\t{n}\n" for c, n in counts.items()))
    return counts


def children_cpu(stdout: str) -> float | None:
    """User + system seconds of the children, from the last line of the shell's `times` ("0m1.2s 0m0.3s" twice, the shell's own, then its children's)."""
    lines = [ln for ln in stdout.splitlines() if re.fullmatch(r"(\d+m[\d.]+s\s*){2}", ln.strip())]
    if len(lines) < 2:
        return None
    return sum(int(m) * 60 + float(s) for m, s in re.findall(r"(\d+)m([\d.]+)s", lines[-1]))


def docker_vep(vcf: Path, gff: Path, fasta: Path, out_name: str, work: Path, *, custom: Path | None = None, extra: list[str] | None = None,
               peak: bool = False) -> tuple[float, float | None, int | None, Path]:
    """One run of VEP in a fresh container: (wall seconds, CPU seconds of VEP, peak resident kB or None, the output file, which the caller removes)."""
    # test_data may be a link to another checkout (a git worktree): resolve every path VEP is given, so that the mount has it.
    vep = ["vep", "--input_file", f"/in/{vcf.name}", "--gff", f"/ref/{gff.name}", "--fasta", f"/ref/{fasta.name}", "--output_file", f"/out/{out_name}", *FLAGS, *(extra or [])]
    mounts = ["-v", f"{vcf.resolve().parent}:/in:ro", "-v", f"{gff.resolve().parent}:/ref:ro", "-v", f"{work}:/out"]
    if fasta.resolve().parent != gff.resolve().parent:
        mounts += ["-v", f"{fasta.resolve().parent}:/ref2:ro"]
        vep[vep.index("--fasta") + 1] = f"/ref2/{fasta.name}"
    if custom is not None:
        mounts += ["-v", f"{custom.resolve()}:/custom:ro"]
    cmd = ["docker", "run", "--rm", "-e", "VEP_DATABASE=0", *mounts, IMAGE]
    if peak:
        # VmHWM is the kernel's high-water mark of the process's resident set; the last value read before it exits is its peak.
        script = ('"$@" >/dev/null 2>&1 & p=$!; m=0; while kill -0 $p 2>/dev/null; do '
                  'v=$(awk \'/VmHWM/ {print $2}\' /proc/$p/status 2>/dev/null); [ -n "$v" ] && m=$v; sleep 0.1; done; wait $p; rc=$?; echo "VmHWM_kB=$m"; times; exit $rc')
    else:
        script = '"$@"; rc=$?; times; exit $rc'
    cmd += ["sh", "-c", script, "sh", *vep]
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
    return wall, children_cpu(r.stdout), kb, out_file


def count_output(out_file: Path, cond: str) -> tuple[int, list[int]]:
    """(records, for each `--custom` source of the condition the records whose first CSQ entry has a value in the column of that source).

    A source that VEP could not read (a wrong path, a name it does not match) leaves its columns out of the CSQ, or empty on every record,
    and is quicker than one that works: the run is refused in the first case and the count shows the second."""
    records = 0
    idx: list[int] = []
    matched = [0] * len(MARKERS[cond])
    with out_file.open() as f:
        for line in f:
            if line.startswith("#"):
                m = re.search(r'ID=CSQ,.*Format: ([^"]+)"', line)
                if m:
                    names = m.group(1).split("|")
                    missing = [c for c in MARKERS[cond] if c not in names]
                    if missing:
                        sys.exit(f"{out_file.name}: no {', '.join(missing)} column in the CSQ of a {cond} run, so the --custom source was not read")
                    idx = [names.index(c) for c in MARKERS[cond]]
                continue
            records += 1
            if idx:
                info = line.split("\t", 8)[7]
                csq = info.split("CSQ=", 1)[1].split(";", 1)[0].split(",", 1)[0].split("|") if "CSQ=" in info else []
                for k, i in enumerate(idx):
                    if len(csq) > i and csq[i]:
                        matched[k] += 1
    return records, matched


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def read_rows(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text().splitlines()] if path.exists() else []


def write_rows(path: Path, rows: list[list[str]]) -> None:
    path.write_text("".join("\t".join(r) + "\n" for r in rows))


def main() -> None:
    global IMAGE
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--runs", type=int, default=3, help="timed runs per cell of scaling, startup and memory; the median is reported (default 3, as for the VEP times of the manuscript)")
    ap.add_argument("--only", default="", help=f"comma-separated parts: {', '.join(PARTS)}")
    ap.add_argument("--conditions", default=",".join(CONDITIONS), help="comma-separated conditions of the genome part (default all three)")
    ap.add_argument("--chroms", default=",".join(CHROMS), help="comma-separated chromosomes of the genome part (default 1 to 22)")
    ap.add_argument("--image", default=IMAGE, help="the VEP image (default: release 116.3 by digest); another release can be timed on the small parts, into a directory of its own")
    ap.add_argument("--clinvar", type=Path, default=SOURCES / "clinvar.vcf.gz", help="the ClinVar VCF that clinvar.osa was built from")
    ap.add_argument("--gnomad-joint", type=Path, default=SOURCES / "gnomad_joint_extracts", help="directory with gnomad_chr*.vcf.gz, the extracts that sa_db_joint was built from")
    args = ap.parse_args()
    IMAGE = args.image
    parts = [p for p in (args.only.split(",") if args.only else PARTS) if p]
    for p in parts:
        if p not in PARTS:
            sys.exit(f"unknown part {p!r}")
    conditions, chroms = args.conditions.split(","), args.chroms.split(",")
    if any(c not in CONDITIONS for c in conditions) or any(c not in CHROMS for c in chroms):
        sys.exit("unknown condition or chromosome")
    needed = [INPUTS / "chr22_all.vcf", INPUTS / "one.vcf", INPUTS / "chr1_100k.vcf", INPUTS / "yeast_one.vcf"]
    if "genome" in parts:
        needed += [HG002, args.clinvar] + [args.gnomad_joint / f"gnomad_chr{c}.vcf.gz" for c in chroms]
    for need in needed:
        if not need.exists():
            sys.exit(f"missing {need}" + (": run benchmarks/run_head_to_head.py first, which cuts the inputs" if need.parent == INPUTS else ""))
    info = sh("docker", "info", "--format", "{{.ServerVersion}} cpus={{.NCPU}} mem={{.MemTotal}}", check=False)
    if not info or "cpus=0" in info:
        sys.exit("Docker's engine is not answering")
    args.out.mkdir(parents=True, exist_ok=True)
    work = ROOT / "benchmarks" / "output" / "vep_side"
    work.mkdir(parents=True, exist_ok=True)
    human_gff, yeast_gff = sorted_gff3(HUMAN / "Homo_sapiens.GRCh38.115.gff3"), sorted_gff3(YEAST / "yeast.gff3")
    human_fa, yeast_fa = HUMAN / "Homo_sapiens.GRCh38.dna.primary_assembly.fa", YEAST / "yeast.fa"
    image_id = sh("docker", "image", "inspect", IMAGE, "--format", "{{.Id}}", check=False)
    # What the image weighs: `.Size` is what is pulled (compressed, with the containerd image store of Docker Desktop), and `docker images` adds that to
    # the unpacked copy; the footprint of an installation is the unpacked root file system, which the container can measure itself.
    image_bytes = sh("docker", "image", "inspect", IMAGE, "--format", "{{.Size}}", check=False)
    image_unpacked_kib = sh("docker", "run", "--rm", "--entrypoint", "sh", IMAGE, "-c", "du -skx / 2>/dev/null | cut -f1", check=False)

    cells = []  # (part, name, vcf, gff, fasta, peak)
    if "scaling" in parts:
        cells += [("scaling", f"chr22_{n}", INPUTS / f"chr22_{n}.vcf", human_gff, human_fa, False) for n in ("1000", "5000", "10000", "all")]
    if "startup" in parts:
        cells += [("startup", "human_1_variant", INPUTS / "one.vcf", human_gff, human_fa, False),
                  ("startup", "yeast_1_variant", INPUTS / "yeast_one.vcf", yeast_gff, yeast_fa, False)]
    if "memory" in parts:
        cells += [("memory", "chr1_100k", INPUTS / "chr1_100k.vcf", human_gff, human_fa, True)]

    raw_path = args.out / "vep_raw.tsv"
    header = ["part", "cell", "round", "wall_s", "cpu_s", "peak_rss_mib", "load1"]
    # Rows of a part that is not run now stay; the rows of a part that is run are replaced.
    raw = [r for r in read_rows(raw_path)[1:] if r[0] not in {c[0] for c in cells}]
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    print(f"Docker {info}; image {image_id[:19]}; {args.runs} timed runs per cell after one that is not timed", flush=True)
    for part, name, vcf, gff, fasta, peak in cells:
        for r in range(0, args.runs + 1):
            if r == 0 and peak:
                continue  # the 100,000-variant run takes about two minutes; its warm-up is the run itself
            # a name of its own for every run: a file that the host has just removed from the shared directory has once been refused to the container
            wall, cpu, kb, out_file = docker_vep(vcf, gff, fasta, f"vep_{name}_r{r}.vcf", work, peak=peak)
            out_file.unlink()
            raw.append([part, name, str(r), f"{wall:.2f}", f"{cpu:.2f}" if cpu is not None else "", f"{kb / 1024:.0f}" if kb else "", f"{os.getloadavg()[0]:.2f}"])
            print(f"  {part:8s} {name:18s} round {r}: {wall:8.2f} s (VEP cpu {cpu if cpu is not None else float('nan'):.2f} s)" + (f"  {kb / 1024:.0f} MiB" if kb else ""), flush=True)
            write_rows(raw_path, [header] + raw)

    genome_path = args.out / "vep_genome.tsv"
    gheader = ["condition", "chrom", "variants", "wall_s", "cpu_s", "out_records", "matched_clinvar", "matched_gnomad", "load1", "finished_utc"]
    custom_dir = None
    if "genome" in parts:
        counts = cut_chromosomes(INPUTS / "hg002_by_chr")
        custom_dir = prepare_custom(args.clinvar, args.gnomad_joint, chroms, work / "custom") if any(c != "plain" for c in conditions) else None
        done = {(r[0], r[1]) for r in read_rows(genome_path)[1:]}
        todo = [(c, cond) for c in chroms for cond in conditions if (cond, c) not in done]
        print(f"genome: {len(todo)} of {len(chroms) * len(conditions)} runs to do ({len(done)} already in {genome_path.name})", flush=True)
        for chrom, cond in todo:
            wall, cpu, _, out_file = docker_vep(INPUTS / "hg002_by_chr" / f"chr{chrom}.vcf", human_gff, human_fa, f"vep_{cond}_chr{chrom}.vcf", work,
                                                custom=custom_dir if cond != "plain" else None, extra=custom_args(cond, chrom))
            records, matched = count_output(out_file, cond)
            out_file.unlink()
            if records != counts[chrom]:
                sys.exit(f"{cond} chr{chrom}: VEP wrote {records:,} records for {counts[chrom]:,} variants")
            rows = read_rows(genome_path) or [gheader]
            by_source = dict(zip(MARKERS[cond], matched))
            rows.append([cond, chrom, str(counts[chrom]), f"{wall:.2f}", f"{cpu:.2f}" if cpu is not None else "", str(records),
                         str(by_source.get("ClinVar", "")), str(by_source.get("gnomAD", "")), f"{os.getloadavg()[0]:.2f}", time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())])
            write_rows(genome_path, rows)
            print(f"  {cond:8s} chr{chrom:<3s} {counts[chrom]:>8,} variants {wall:9.2f} s (VEP cpu {cpu if cpu is not None else float('nan'):.2f} s)"
                  + "".join(f", {n:,} matched {name}" for name, n in by_source.items()), flush=True)

    summary = [["part", "cell", "n", "median_wall_s", "min_wall_s", "max_wall_s", "median_cpu_s", "median_peak_rss_mib"]]
    names = list(dict.fromkeys((r[0], r[1]) for r in raw))
    med: dict[str, float] = {}
    for part, name in names:
        rows = [r for r in raw if r[0] == part and r[1] == name and (int(r[2]) >= 1 or part == "memory")]
        w = [float(r[3]) for r in rows]
        cpus = [float(r[4]) for r in rows if r[4]]
        peaks = [float(r[5]) for r in rows if r[5]]
        summary.append([part, name, str(len(w)), f"{statistics.median(w):.2f}", f"{min(w):.2f}", f"{max(w):.2f}", f"{statistics.median(cpus):.2f}" if cpus else "",
                        f"{statistics.median(peaks):.0f}" if peaks else ""])
        med[f"{part}:{name}"] = statistics.median(w)
        if peaks:
            med[f"peak:{name}"] = statistics.median(peaks)
    if raw:
        write_rows(args.out / "vep_summary.tsv", summary)

    genome = read_rows(genome_path)[1:]
    if genome:
        sums = [["condition", "chromosomes", "variants", "wall_s", "cpu_s", "matched_clinvar", "matched_gnomad"]]
        baseline = [["measurement", "variants", "vep_sec"]]
        for cond in CONDITIONS:
            rows = [r for r in genome if r[0] == cond]
            if rows:
                wall, cpu = sum(float(r[3]) for r in rows), sum(float(r[4]) for r in rows if r[4])
                sums.append([cond, str(len(rows)), str(sum(int(r[2]) for r in rows)), f"{wall:.1f}", f"{cpu:.1f}",
                             str(sum(int(r[6]) for r in rows if r[6])) if cond != "plain" else "", str(sum(int(r[7]) for r in rows if r[7])) if cond == "joint" else ""])
        (args.out / "vep_genome_sums.csv").write_text("".join(",".join(r) + "\n" for r in sums))
        # The layout of manuscript/data/vep_baseline_*.csv; a sum is a whole-genome time only when all 22 chromosomes are in it.
        names = {"plain": "wgs_conseq", "clinvar": "wgs_clinvar", "joint": "wgs_gnomad_clinvar"}
        for row in summary[1:]:
            if row[0] == "scaling":
                baseline.append(["chr22_scaling", str(SCALING_VARIANTS[row[1]]), f"{float(row[3]):.2f}"])
        for r in sums[1:]:
            if int(r[1]) == len(CHROMS):
                baseline.append([names[r[0]], r[2], f"{float(r[3]):.0f}"])
        (args.out / "vep_baseline.csv").write_text("".join(",".join(r) + "\n" for r in baseline))

    # The parts can be run by separate commands: what an earlier one recorded stays unless this one has something to say, and the start is the first.
    env_path = args.out / "vep_environment.tsv"
    env = dict(line.rstrip("\n").split("\t", 1) for line in env_path.open()) if env_path.exists() else {}
    this = {"image": IMAGE, "image_id": image_id, "image_bytes": image_bytes, "image_unpacked_kib": image_unpacked_kib, "docker": info, "flags": " ".join(FLAGS),
            "custom_clinvar": " ".join(custom_args("clinvar", "22")), "custom_gnomad_joint": " ".join(custom_args("joint", "22")[2:]),
            "clinvar_vcf": f"{args.clinvar.name} md5 {md5(args.clinvar)}" if "genome" in parts else "", "runs": str(args.runs),
            "load_end": " ".join(f"{x:.2f}" for x in os.getloadavg())}
    env.update({k: v for k, v in this.items() if v})
    env["started_utc"] = min(env.get("started_utc", started), started)
    env_path.write_text("".join(f"{k}\t{v}\n" for k, v in env.items()))
    resources = [["metric", "value"]]
    if "startup:yeast_1_variant" in med and "startup:human_1_variant" in med:
        resources.append(["startup_yeast", f"{med['startup:yeast_1_variant']:.1f} sec"])
        resources.append(["startup_human", f"{med['startup:human_1_variant']:.1f} sec"])
    if "peak:chr1_100k" in med:
        resources.append(["peak_memory", f"{med['peak:chr1_100k'] / 1024:.1f} GiB"])
    if image_unpacked_kib.isdigit():
        resources.append(["image_size", f"{int(image_unpacked_kib) * 1024 / 1e9:.2f} GB unpacked"])
    if image_bytes.isdigit():
        resources.append(["image_pull_size", f"{int(image_bytes) / 1e9:.2f} GB"])
    if len(resources) > 1:
        (args.out / "vep_resources.csv").write_text("".join(",".join(row) + "\n" for row in resources))
    print(f"\nwrote the files of {args.out} (vep_raw.tsv, vep_summary.tsv, vep_environment.tsv" + (", vep_genome.tsv, vep_genome_sums.csv, vep_baseline.csv" if genome else "")
          + (", vep_resources.csv" if len(resources) > 1 else "") + ")")


if __name__ == "__main__":
    main()
