#!/usr/bin/env python3
"""The fastVEP side of the head-to-head and resource tables, measured in one place.

Tables 2 to 5 of the manuscript and the performance section of the README quote times, peak memory and sizes of one
`fastvep` binary on the human gene model of Ensembl 115. They used to come out of a shell session, so a number in a
document could not be traced to a command, a binary or a day. This script is that command: it writes every run it
made (`raw.tsv`), their medians (`summary.tsv`) and what it ran on (`environment.tsv`) into one directory, and the
documents quote that directory.

    benchmarks/run_head_to_head.py --bin target/release/fastvep --built-from "v0.5.0-rc.1 (c12d9f9)" \
        --out benchmarks/results/v0.5.0-rc.1
    benchmarks/run_head_to_head.py --bin ... --only scaling,memory --runs 3

What is measured, always with the FASTA and `--hgvs --symbol --canonical` (the flags of the head-to-head with Ensembl
VEP) and with the transcript cache already built:

    first-run  parse the GFF3 and write the sidecar cache, once: wall time, peak memory, size of the cache
    startup    one variant, which is what reading the cache costs: one thread, eight and every core; yeast as well
    scaling    Table 2: the first 1,000, 5,000 and 10,000 chr22 variants of GIAB HG002 and all 50,284, one thread
    memory     Table 4: peak memory of 100,000 chr1 variants on 1, 8 and every thread; and of a run without a FASTA
    formats    Table 5: 10,000 chr22 variants written as VCF, tab and JSON, one thread
    wgs        Table 3: the whole of HG002 (4,048,342 variants) on one thread and on every core: consequences only,
               with ClinVar, and with gnomAD and ClinVar
    organisms  one variant against each of the five gene models: how many transcripts each holds (`gene_models.tsv`)

The size of the binary, and the SHA-256 of every input, are in `environment.tsv`.

`benchmarks/run_benchmark.sh` measures the five-organism throughput (Table 7) and is run beside this one.

Why it is built the way it is:

* Rounds, not blocks. Each cell of a part is run once per round and the rounds are repeated, so a slow minute on the
  machine (a backup, a browser tab) falls on every cell and not on one of them. The median of the rounds is reported.
  Round 0 is a warm-up that is not reported (it reads every input once, and is where the outputs are hashed).
* A run that fails is not a run. The script stops and prints the command; it never times how long a failure took.
* One process per run, timed from outside. The peak memory comes from `wait4`, which reports the child that was
  waited for; `getrusage(RUSAGE_CHILDREN)` reports the largest child so far, which would give every later, smaller run
  the peak of an earlier one.
* The output of the runs that must agree is compared, byte for byte, in round 0 (md5): the whole genome on one
  thread against the same on every core, and the same for the 100,000-variant runs. A speed-up that changes an answer
  is not one, and this is where a threading bug would show.
* The script only reads the inputs, the databases and the cache that the repository keeps, and records the SHA-256 of
  the inputs and the md5 of the cache before and after. A run that needs sequences its cache lacks rewrites the
  cache, which would make every later run measure something else. The only thing it deletes is the sidecar cache of
  the gene model, which the first-run part rebuilds at once with the binary being measured: a cache that an older
  build wrote is read by this one but was not cut into frames, so a run on it decodes on one thread.
* The inputs are cut from the GIAB VCF by this script, in one rule, so the cuts are reproducible: the first N
  records of chromosome 22 (or 1), kept to their first five columns with `.` for the rest (the INFO of the benchmark
  VCF is read by neither tool for the consequences), in a directory the repository ignores.

Linux and macOS only (`os.wait4`). Standard library only.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import platform
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ORGANISMS = ROOT / "test_data" / "organisms"
HUMAN = ORGANISMS / "human"
GFF3 = HUMAN / "Homo_sapiens.GRCh38.115.gff3"
FASTA = HUMAN / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
HG002 = HUMAN / "human_giab_hg002_full.vcf"
YEAST = ORGANISMS / "yeast"
DEFAULT_SA_DB = ROOT / "data" / "benchmark" / "sa_db"
WORK = ROOT / "benchmarks" / "output" / "head_to_head"

# --hgvs needs the FASTA; --symbol and --canonical are the two flags the Ensembl VEP side of the comparison is given.
FLAGS = ["--hgvs", "--symbol", "--canonical", "--no-progress"]
ALL_CORES = os.cpu_count() or 1
PARTS = ("first-run", "startup", "scaling", "memory", "formats", "wgs", "organisms")
RSS_UNIT = 1024 * 1024 if sys.platform == "darwin" else 1024  # ru_maxrss is bytes on macOS and KiB on Linux


def sidecar_of(gff3: Path) -> Path:
    # The path as given plus a suffix, never the resolved path: `default_cache_path` in fastvep-cache does the same, which is
    # what lets a symlink to the GFF3 in another directory keep a sidecar of its own.
    return gff3.with_name(gff3.name + ".fastvep.cache")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def sh(*cmd: str) -> str:
    try:
        return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


@dataclass
class Cell:
    part: str
    name: str
    argv: list[str]  # after the binary
    threads: int | None  # None: the default, every core
    out: Path | None = None  # where the run writes, removed after it; its size is recorded
    group: str = ""  # cells of a group must write the same bytes
    runs: list[dict] = field(default_factory=list)


class Recorder:
    def __init__(self, binary: Path, outdir: Path):
        self.bin = binary
        self.outdir = outdir
        self.raw: list[dict] = []
        self.md5s: dict[str, dict[str, str]] = {}  # group -> {cell: md5}
        self.failed = False

    def run(self, cell: Cell, round_no: int) -> None:
        env = {k: v for k, v in os.environ.items() if k != "RAYON_NUM_THREADS"}
        if cell.threads is not None:
            env["RAYON_NUM_THREADS"] = str(cell.threads)
        cmd = [str(self.bin), *cell.argv]
        if cell.out is not None:
            cell.out.parent.mkdir(parents=True, exist_ok=True)
            cell.out.unlink(missing_ok=True)
        with tempfile.TemporaryFile() as err:
            t0 = time.perf_counter()
            proc = subprocess.Popen(cmd, env=env, stdout=subprocess.DEVNULL, stderr=err)
            _, status, ru = os.wait4(proc.pid, 0)
            wall = time.perf_counter() - t0
            proc.returncode = os.waitstatus_to_exitcode(status)
            err.seek(0)
            log = err.read().decode("utf-8", "replace")
            tail = log[-600:]
        if proc.returncode != 0:
            sys.exit(f"FAILED (exit {proc.returncode}): RAYON_NUM_THREADS={cell.threads} {' '.join(cmd)}\n{tail}")
        row = {
            "part": cell.part,
            "cell": cell.name,
            "threads": "all" if cell.threads is None else cell.threads,
            "round": round_no,
            "wall_s": f"{wall:.3f}",
            "user_s": f"{ru.ru_utime:.3f}",
            "sys_s": f"{ru.ru_stime:.3f}",
            "peak_rss_mib": f"{ru.ru_maxrss / RSS_UNIT:.0f}",
            "out_bytes": "",
            "md5": "",
            "load1": f"{os.getloadavg()[0]:.2f}",
            "transcripts": (re.search(r"Loaded (\d+) transcripts", log) or [None, ""])[1],
        }
        if cell.out is not None and cell.out.exists():
            row["out_bytes"] = str(cell.out.stat().st_size)
            if round_no == 0 and (cell.group or cell.part in ("scaling", "formats")):
                row["md5"] = md5(cell.out)
                if cell.group:
                    self.md5s.setdefault(cell.group, {})[cell.name] = row["md5"]
            cell.out.unlink()
        self.raw.append(row)
        print(f"  {cell.part:9s} {cell.name:34s} round {round_no}: {wall:8.2f} s  {row['peak_rss_mib']:>6s} MiB", flush=True)

    def rounds(self, cells: list[Cell], runs: int) -> None:
        # Round 0 is not timed. It warms the file cache, so that the first timed run of a cell does not read the VCF from disk
        # while the others read it from memory, and it is where the outputs are hashed: hashing 10 GB of output between two timed
        # runs would leave the next one with a file cache full of the last one's output.
        for r in range(0, runs + 1):
            for cell in cells:
                self.run(cell, r)

    def identity(self) -> list[tuple[str, bool, dict[str, str]]]:
        out = []
        for group, cells in sorted(self.md5s.items()):
            out.append((group, len(set(cells.values())) == 1, cells))
        return out

    def write(self) -> None:
        cols = ["part", "cell", "threads", "round", "wall_s", "user_s", "sys_s", "peak_rss_mib", "out_bytes", "md5", "load1", "transcripts"]
        with (self.outdir / "raw.tsv").open("w") as f:
            f.write("\t".join(cols) + "\n")
            for r in self.raw:
                f.write("\t".join(str(r[c]) for c in cols) + "\n")
        cells: dict[tuple[str, str, str], list[dict]] = {}
        for r in self.raw:
            if int(r["round"]) >= 1:  # round 0 is the warm-up
                cells.setdefault((r["part"], r["cell"], str(r["threads"])), []).append(r)
        scols = ["part", "cell", "threads", "n", "median_wall_s", "min_wall_s", "max_wall_s", "median_peak_rss_mib", "busy_cores", "out_bytes"]
        with (self.outdir / "summary.tsv").open("w") as f:
            f.write("\t".join(scols) + "\n")
            for (part, name, threads), rows in cells.items():
                walls = [float(r["wall_s"]) for r in rows]
                busy = [(float(r["user_s"]) + float(r["sys_s"])) / float(r["wall_s"]) for r in rows]
                size = next((r["out_bytes"] for r in rows if r["out_bytes"]), "")
                f.write("\t".join([
                    part, name, threads, str(len(rows)), f"{statistics.median(walls):.2f}", f"{min(walls):.2f}", f"{max(walls):.2f}",
                    f"{statistics.median(float(r['peak_rss_mib']) for r in rows):.0f}", f"{statistics.median(busy):.2f}", size,
                ]) + "\n")


def lean_vcf(src: Path, dst: Path, chrom: str, limit: int | None) -> int:
    """The first `limit` records of `chrom`, kept to CHROM POS ID REF ALT and `.` for QUAL FILTER INFO."""
    n = 0
    with src.open() as fin, dst.open("w") as fout:
        fout.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for line in fin:
            if line.startswith("#"):
                continue
            f = line.split("\t", 5)
            if f[0] != chrom:
                continue
            fout.write("\t".join(f[:5]) + "\t.\t.\t.\n")
            n += 1
            if limit is not None and n >= limit:
                break
    return n


def prepare_inputs(sa_db: Path) -> dict[str, Path]:
    """The cuts of the HG002 VCF and the database directories the cells read, made once and reused."""
    d = WORK / "inputs"
    d.mkdir(parents=True, exist_ok=True)
    inputs: dict[str, Path] = {}
    # One pass over the 2.4 GB file for chr22, a second for chr1; each stops when it has enough, which for chr1 is early.
    chr22_all = d / "chr22_all.vcf"
    if not chr22_all.exists():
        print("cutting chr22 ...", flush=True)
        lean_vcf(HG002, chr22_all, "22", None)
    inputs["chr22_all"] = chr22_all
    for n in (1000, 5000, 10000):
        p = d / f"chr22_{n}.vcf"
        if not p.exists():
            with chr22_all.open() as src, p.open("w") as dst:
                kept = 0
                for line in src:
                    if line.startswith("#") or kept < n:
                        dst.write(line)
                        kept += 0 if line.startswith("#") else 1
        inputs[f"chr22_{n}"] = p
    chr1 = d / "chr1_100k.vcf"
    if not chr1.exists():
        print("cutting chr1 ...", flush=True)
        lean_vcf(HG002, chr1, "1", 100_000)
    inputs["chr1_100k"] = chr1
    one = d / "one.vcf"
    if not one.exists():
        lean_vcf(HG002, one, "22", 1)
    inputs["one"] = one
    yeast_one = d / "yeast_one.vcf"
    if not yeast_one.exists():
        with (YEAST / "yeast_ensembl_full.vcf").open() as src, yeast_one.open("w") as dst:
            for line in src:
                dst.write(line)
                if not line.startswith("#"):
                    break
    inputs["yeast_one"] = yeast_one
    # A GFF3 path of its own for the run without a FASTA, so that its cache (built without sequences) is not the main one.
    nofasta = d / "nofasta" / "grch38.gff3"
    nofasta.parent.mkdir(exist_ok=True)
    if not nofasta.is_symlink():
        nofasta.symlink_to(GFF3)
    inputs["nofasta_gff3"] = nofasta
    # The two annotation conditions of Table 3. Symbolic links: nothing here is written to.
    if (sa_db / "clinvar.osa").exists():
        for cond, files in (("sa_clinvar", ["clinvar.osa", "clinvar.osa.idx"]),
                            ("sa_both", ["clinvar.osa", "clinvar.osa.idx", *sorted(p.name for p in sa_db.glob("gnomad_chr*.osa2"))])):
            cd = d / cond
            if cd.exists():
                shutil.rmtree(cd)
            cd.mkdir()
            for name in files:
                (cd / name).symlink_to(sa_db / name)
            inputs[cond] = cd
    return inputs


def annotate(vcf: Path, out: Path, *, fmt: str = "vcf", gff3: Path = GFF3, fasta: Path | None = FASTA, sa: Path | None = None,
             flags: bool = True) -> list[str]:
    """The arguments of one run (after the binary). `flags=False` with no FASTA is the consequence-only run."""
    argv = ["annotate", "-i", str(vcf), "-o", str(out), "--output-format", fmt, "--gff3", str(gff3)]
    if fasta is not None:
        argv += ["--fasta", str(fasta)]
    argv += FLAGS if flags else ["--no-progress"]
    if sa is not None:
        argv += ["--sa-dir", str(sa)]
    return argv


def environment(rec: Recorder, args: argparse.Namespace, started: str, inputs: dict[str, Path]) -> None:
    binary = rec.bin
    rows = [
        ("started_utc", started),
        ("system", f"{platform.uname().system} {platform.uname().release} {platform.uname().machine}"),
        ("macos_product_version", sh("sw_vers", "-productVersion")),
        ("cpu", sh("sysctl", "-n", "machdep.cpu.brand_string") or platform.processor()),
        ("logical_cores", str(ALL_CORES)),
        ("performance_cores", sh("sysctl", "-n", "hw.perflevel0.logicalcpu")),  # Apple silicon only; empty elsewhere
        ("efficiency_cores", sh("sysctl", "-n", "hw.perflevel1.logicalcpu")),
        ("hardware_model", sh("sysctl", "-n", "hw.model")),
        ("swap_in_use", sh("sysctl", "-n", "vm.swapusage")),
        ("memory_gib", f"{int(sh('sysctl', '-n', 'hw.memsize') or 0) / 2**30:.0f}" if sys.platform == "darwin" else ""),
        ("rustc", sh("rustc", "--version")),
        ("python", platform.python_version()),
        ("binary", str(binary)),
        ("binary_version", sh(str(binary), "--version")),
        ("binary_bytes", str(binary.stat().st_size)),
        ("binary_sha256", sha256(binary)),
        ("built_from", args.built_from),
        ("repo_head", sh("git", "-C", str(ROOT), "rev-parse", "HEAD")),
        ("runs_per_cell", str(args.runs)),
        ("parts", ",".join(args.parts)),
        ("loadavg_start", " ".join(f"{x:.2f}" for x in os.getloadavg())),
        ("flags", " ".join(FLAGS)),
        ("gff3_sha256", sha256(GFF3)),
        ("fasta_sha256", sha256(FASTA)),
        ("hg002_vcf_sha256", sha256(HG002)),
        ("hg002_records", str(sum(1 for line in HG002.open() if not line.startswith("#")))),
        ("chr22_records", str(sum(1 for line in inputs["chr22_all"].open() if not line.startswith("#")))),
    ]
    sa = args.sa_db
    if (sa / "clinvar.osa").exists():
        rows.append(("clinvar_osa_sha256", sha256(sa / "clinvar.osa")))
        rows.append(("gnomad_osa2_files", str(len(list(sa.glob("gnomad_chr*.osa2"))))))
    with (rec.outdir / "environment.tsv").open("w") as f:
        for k, v in rows:
            f.write(f"{k}\t{v}\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", required=True, type=Path, help="the fastvep binary to measure")
    ap.add_argument("--out", required=True, type=Path, help="directory for raw.tsv, summary.tsv and environment.tsv")
    ap.add_argument("--runs", type=int, default=5, help="rounds per cell; the median is reported (default 5)")
    ap.add_argument("--only", default="", help=f"comma-separated parts to run: {', '.join(PARTS)} (default all)")
    ap.add_argument("--built-from", default="not stated", help="what the binary was built from, recorded in environment.tsv")
    ap.add_argument("--sa-db", type=Path, default=DEFAULT_SA_DB, help="directory with clinvar.osa and gnomad_chr*.osa2 (default data/benchmark/sa_db)")
    ap.add_argument("--list", action="store_true", help="print the cells and stop")
    args = ap.parse_args()
    args.parts = [p for p in (args.only.split(",") if args.only else PARTS) if p]
    for p in args.parts:
        if p not in PARTS:
            sys.exit(f"unknown part {p!r}; the parts are {', '.join(PARTS)}")
    for f in (args.bin, GFF3, FASTA, HG002):
        if not f.exists():
            sys.exit(f"missing: {f}")
    args.out.mkdir(parents=True, exist_ok=True)
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    inputs = prepare_inputs(args.sa_db)
    tmp = WORK / "tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    rec = Recorder(args.bin.resolve(), args.out)
    sidecar = sidecar_of(GFF3)
    cache_before = md5(sidecar) if sidecar.exists() else "absent"

    one, yeast_one = inputs["one"], inputs["yeast_one"]
    cells: dict[str, list[Cell]] = {p: [] for p in PARTS}
    cells["startup"] = [
        Cell("startup", "human_1_variant_t1", annotate(one, tmp / "o.vcf"), 1, tmp / "o.vcf"),
        Cell("startup", "human_1_variant_t8", annotate(one, tmp / "o.vcf"), 8, tmp / "o.vcf"),
        Cell("startup", "human_1_variant_all", annotate(one, tmp / "o.vcf"), None, tmp / "o.vcf"),
        Cell("startup", "yeast_1_variant_t1", annotate(yeast_one, tmp / "o.vcf", gff3=YEAST / "yeast.gff3", fasta=YEAST / "yeast.fa"), 1, tmp / "o.vcf"),
    ]
    cells["scaling"] = [Cell("scaling", f"chr22_{n}_t1", annotate(inputs[f"chr22_{n}"], tmp / "o.vcf"), 1, tmp / "o.vcf") for n in ("1000", "5000", "10000")]
    cells["scaling"].append(Cell("scaling", "chr22_all_t1", annotate(inputs["chr22_all"], tmp / "o.vcf"), 1, tmp / "o.vcf"))
    mem_vcf = inputs["chr1_100k"]
    cells["memory"] = [Cell("memory", f"chr1_100k_t{t}", annotate(mem_vcf, tmp / "o.vcf"), None if t == "all" else int(t), tmp / "o.vcf", group="memory")
                       for t in ("1", "8", "all")]
    cells["memory"] += [Cell("memory", f"chr1_100k_nofasta_conseq_t{t}", annotate(mem_vcf, tmp / "o.vcf", gff3=inputs["nofasta_gff3"], fasta=None, flags=False),
                             int(t), tmp / "o.vcf", group="memory_nofasta") for t in ("1", "8")]
    cells["formats"] = [Cell("formats", f"chr22_10000_{fmt}_t1", annotate(inputs["chr22_10000"], tmp / f"o.{fmt}", fmt=fmt), 1, tmp / f"o.{fmt}")
                        for fmt in ("vcf", "tab", "json")]
    for cond, sa in (("plain", None), ("clinvar", inputs.get("sa_clinvar")), ("both", inputs.get("sa_both"))):
        if cond != "plain" and sa is None:
            print(f"skipping the {cond} condition: no databases in {args.sa_db}")
            continue
        for t in (1, None):
            label = "t1" if t == 1 else "all"
            cells["wgs"].append(Cell("wgs", f"hg002_{cond}_{label}", annotate(HG002, tmp / "wgs.vcf", sa=sa), t, tmp / "wgs.vcf", group=f"wgs_{cond}"))

    organisms = []
    for name, d, gff3, fasta, vcf in (("yeast", "yeast", "yeast.gff3", "yeast.fa", "yeast_ensembl_full.vcf"),
                                      ("drosophila", "drosophila", "drosophila.gff3", "drosophila.fa", "drosophila_dgrp2_full.vcf"),
                                      ("arabidopsis", "arabidopsis", "arabidopsis.gff3", "arabidopsis.fa", "arabidopsis_1001g_full.vcf"),
                                      ("mouse", "mouse", "mouse.gff3", "mouse.fa", "mouse_cast_eij_full.vcf" if (ORGANISMS / "mouse" / "mouse_cast_eij_full.vcf").exists() else "mouse_ensembl_1m.vcf"),
                                      ("human", "human", GFF3.name, FASTA.name, HG002.name)):
        paths = [ORGANISMS / d / f for f in (gff3, fasta, vcf)]
        if all(p.exists() for p in paths):
            organisms.append((name, *paths))
    for name, gff3, fasta, vcf in organisms:
        first = WORK / "inputs" / f"{name}_first_record.vcf"
        if not first.exists():
            with vcf.open() as src, first.open("w") as dst:
                for line in src:
                    dst.write(line)
                    if not line.startswith("#"):
                        break
        cells["organisms"].append(Cell("organisms", f"{name}_1_variant", annotate(first, tmp / "o.vcf", gff3=gff3, fasta=fasta, flags=False) + ["--hgvs"], None, tmp / "o.vcf"))

    if args.list:
        for p in args.parts:
            for c in cells[p]:
                print(p, c.name, f"threads={c.threads}", " ".join(c.argv))
        return

    environment(rec, args, started, inputs)
    print(f"binary: {rec.bin}  ({sh(str(rec.bin), '--version')})  runs per cell: {args.runs}  load: {os.getloadavg()[0]:.2f}", flush=True)

    if "first-run" in args.parts:
        print("== first run: parse the GFF3 and write the cache", flush=True)
        todo = [(GFF3, "human_with_sequences", FASTA, True), (inputs["nofasta_gff3"], "human_without_sequences", None, False)]
        if (YEAST / "yeast.gff3").exists():
            todo.append((YEAST / "yeast.gff3", "yeast_with_sequences", YEAST / "yeast.fa", True))
        for gff3, name, fasta, flags in todo:
            sidecar_of(gff3).unlink(missing_ok=True)
            vcf = yeast_one if name.startswith("yeast") else one
            cell = Cell("first-run", name, annotate(vcf, tmp / "o.vcf", gff3=gff3, fasta=fasta, flags=flags), None, tmp / "o.vcf")
            rec.run(cell, 1)
            rec.raw[-1]["out_bytes"] = str(sidecar_of(gff3).stat().st_size)  # what the first run leaves behind is the cache
    if "first-run" not in args.parts and not sidecar.exists():
        sys.exit(f"no sidecar cache at {sidecar}: run the first-run part first")
    for part in ("startup", "scaling", "memory", "formats", "wgs"):
        if part in args.parts:
            print(f"== {part}", flush=True)
            rec.rounds(cells[part], args.runs)
    if "organisms" in args.parts:
        print("== organisms: the gene models", flush=True)
        rec.rounds(cells["organisms"], 1)
        with (args.out / "gene_models.tsv").open("w") as f:
            f.write("organism\ttranscripts\tgff3_bytes\tgff3_sha256\tfasta_bytes\tfasta_sha256\tvcf\tvcf_bytes\n")
            by_cell = {r["cell"]: r for r in rec.raw if r["part"] == "organisms"}
            for name, gff3, fasta, vcf in organisms:
                f.write("\t".join([name, by_cell[f"{name}_1_variant"]["transcripts"], str(gff3.stat().st_size), sha256(gff3), str(fasta.stat().st_size),
                                    sha256(fasta), vcf.name, str(vcf.stat().st_size)]) + "\n")
    rec.write()
    cache_after = md5(sidecar) if sidecar.exists() else "absent"
    bad = [g for g, same, _ in rec.identity() if not same]
    with (args.out / "identity.tsv").open("w") as f:
        f.write("group\tidentical\tcell\tmd5\n")
        for group, same, cs in rec.identity():
            for name, digest in sorted(cs.items()):
                f.write(f"{group}\t{'yes' if same else 'NO'}\t{name}\t{digest}\n")
    with (args.out / "environment.tsv").open("a") as f:
        f.write(f"loadavg_end\t{' '.join(f'{x:.2f}' for x in os.getloadavg())}\n")
        f.write(f"sidecar_md5_before\t{cache_before}\n")
        f.write(f"sidecar_md5_after\t{cache_after}\n")
        f.write(f"finished_utc\t{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\n")
    print(f"\nwrote {args.out}/raw.tsv, summary.tsv, identity.tsv, environment.tsv")
    for group, same, _ in rec.identity():
        print(f"  identical output across threads, {group}: {'yes' if same else 'NO'}")
    if "first-run" not in args.parts and cache_before != cache_after:
        print("WARNING: the sidecar cache changed during the runs; a run rewrote it, so later runs did not read the cache the first ones read")
        bad.append("sidecar")
    if bad:
        sys.exit(f"FAILED checks: {', '.join(bad)}")


if __name__ == "__main__":
    main()
