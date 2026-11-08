#!/usr/bin/env python3
"""Write the CSV files of `manuscript/data/` from the runs that measured them, and check that they still agree.

Until this script, the files of `manuscript/data/` had no generator: their numbers were typed from a shell session
(July 2026), and when the code changed nothing said which of them had gone stale. Now each file is a function of
files in a results directory that a command made, and `--check` fails when a file differs from what its inputs say.

    manuscript/build_data.py --results benchmarks/results/v0.5.0            # write the files
    manuscript/build_data.py --results benchmarks/results/v0.5.0 --check    # fail if any file is out of step

The results directory holds what these commands wrote:

    benchmarks/run_head_to_head.py    raw.tsv, gene_models.tsv, environment.tsv (and summary.tsv, identity.tsv); the parts that are compared with
                                      Ensembl VEP (startup, scaling, memory, formats, wgs and wgs_joint) were run again in the session of the VEP runs,
                                      into the subdirectory `with_vep/`, and a cell that is in both is read from there
    benchmarks/run_vep_side.py        with_vep/vep_baseline.csv (Ensembl VEP's times), vep_summary.tsv, vep_environment.tsv (and vep_raw.tsv, vep_genome.tsv)
    validation/measure_vep116.py      vep116_agreement.csv and vep116_hg002_columns.csv (the counts of the comparison with Ensembl VEP 116.3)
    benchmarks/run_benchmark.sh       multi_organism.csv (benchmark_results.csv, copied) and multi_organism_environment.tsv
    benchmarks/consequence_distribution.py   consequence_most_severe.csv and consequence_every_term.csv
    manuscript/code_metrics.py        code_metrics.csv, so_terms.csv, crate_dependencies.csv
    benchmarks/run_sa_formats.py      sa_formats.csv, sa_formats_annotation.csv (Table 6)
    cargo test --workspace            tests_executed.csv (what ran and passed)

The times of Ensembl VEP are read from `with_vep/vep_baseline.csv` of the results directory: VEP 116.3 ran in Docker on the same Apple M4, in the
session of the fastVEP runs of `with_vep/`, and its whole-genome time is the sum over chromosomes 1 to 22 of one run of each (`vep_genome.tsv` has them one
by one, with the CPU time of VEP beside the wall time). A speed-up in the manuscript is that time over the time of fastVEP of the same session: a laptop
runs the same binary at speeds that differ by several per cent from one session to another, which is more than the spread of the rounds of one.
`vep116_provenance.md` accounts for the counts of the comparison with Ensembl VEP 116.3, which `verify_numbers.py` checks against
`docs/VEP_DIVERGENCE.md`.

A file whose input is not in the results directory is left as it is, and `--check` says so.

Standard library only.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

ORGANISMS = [  # (key in the results, label of the CSVs, assembly, variant source)
    ("yeast", "yeast", "R64-1-1", "Ensembl/SGD"),
    ("drosophila", "drosophila", "BDGP6.54", "DGRP2"),
    ("arabidopsis", "arabidopsis", "TAIR10", "1001_Genomes"),
    ("mouse", "mouse", "GRCm39", "MGP_CAST/EiJ"),
    ("human", "human", "GRCh38", "GIAB_HG002"),
]
SCALING_SIZES = ("1000", "5000", "10000", "all")
# (key of the VEP baseline, label, the harness's condition, its part): the third condition reads the joint release of gnomAD, the one
# whose source files Ensembl VEP is given as well (`run_vep_side.py`); the condition of the v4.1 exomes is in raw.tsv and is not tabulated.
CONDITIONS = [("conseq", "Consequence", "plain", "wgs"), ("clinvar", "+ ClinVar", "clinvar", "wgs"), ("gnomad_clinvar", "+ gnomAD + ClinVar", "joint", "wgs_joint")]
# Files that a run of the release made in the results directory and that the manuscript's data directory holds as they are.
COPIED = {"vep_baseline.csv": "with_vep", "vep116_agreement.csv": "", "vep116_hg002_columns.csv": "", "code_metrics.csv": "", "so_terms.csv": "",
          "crate_dependencies.csv": "", "sa_formats.csv": "", "sa_formats_annotation.csv": "", "tests_executed.csv": ""}  # name -> subdirectory of the results


def read_rows(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


class Results:
    def __init__(self, directory: Path):
        self.dir = directory
        raw = directory / "raw.tsv"
        # round 0 is the harness's untimed warm-up and identity round, and is in no median
        self.raw = [r for r in read_rows(raw) if int(r["round"]) >= 1] if raw.exists() else []
        # the later session, the one of the Ensembl VEP runs: a cell that it measured is read from it and not from the first session
        later = directory / "with_vep" / "raw.tsv"
        self.later = [r for r in read_rows(later) if int(r["round"]) >= 1] if later.exists() else []
        self.vep = directory / "with_vep"
        lenv = directory / "with_vep" / "environment.tsv"
        self.env_later = dict(line.rstrip("\n").split("\t", 1) for line in lenv.open()) if lenv.exists() else {}
        gm = directory / "gene_models.tsv"
        self.gene_models = {r["organism"]: r for r in read_rows(gm)} if gm.exists() else {}
        mo = directory / "multi_organism.csv"
        self.multi = {r["dataset"]: r for r in read_rows(mo, ",")} if mo.exists() else {}
        env = directory / "environment.tsv"
        self.env = dict(line.rstrip("\n").split("\t", 1) for line in env.open()) if env.exists() else {}

    def rows(self, part: str, cell: str) -> list[dict[str, str]]:
        mine = [r for r in self.later if r["part"] == part and r["cell"] == cell]
        return mine or [r for r in self.raw if r["part"] == part and r["cell"] == cell]

    def wall(self, part: str, cell: str) -> float | None:
        v = [float(r["wall_s"]) for r in self.rows(part, cell)]
        return statistics.median(v) if v else None

    def rss(self, part: str, cell: str) -> float | None:
        v = [float(r["peak_rss_mib"]) for r in self.rows(part, cell)]
        return statistics.median(v) if v else None

    def out_bytes(self, part: str, cell: str) -> int | None:
        v = [int(r["out_bytes"]) for r in self.rows(part, cell) if r["out_bytes"]]
        return v[0] if v else None


def baseline(res: Results) -> dict[tuple[str, str], float]:
    path = res.vep / "vep_baseline.csv"
    return {(r["measurement"], r["variants"]): float(r["vep_sec"]) for r in read_rows(path, ",")} if path.exists() else {}


def fmt(x: float, digits: int = 2) -> str:
    return f"{x:.{digits}f}"


def speedup(vep: float, fast: float) -> str:
    s = vep / fast
    return f"{s:.2f}" if s < 1 else f"{s:.1f}"


def build(res: Results) -> dict[str, list[list[str]]]:
    """name of the CSV -> its rows, header first. A file with no input here is not in the dict."""
    out: dict[str, list[list[str]]] = {}
    base = baseline(res)

    # Table 2: chr22 scaling, one thread, against VEP's times of the same day.
    sizes = {"1000": 1000, "5000": 5000, "10000": 10000, "all": 50284}
    if all(res.wall("scaling", f"chr22_{k}_t1") for k in SCALING_SIZES) and all(("chr22_scaling", str(n)) in base for n in sizes.values()):
        rows = [["variants", "fastvep_sec", "fastvep_vps", "vep_sec", "vep_vps", "speedup"]]
        for k, n in sizes.items():
            f = res.wall("scaling", f"chr22_{k}_t1")
            v = base[("chr22_scaling", str(n))]
            rows.append([str(n), fmt(f), str(round(n / f)), fmt(v), str(round(n / v)), speedup(v, f)])
        out["vep_comparison.csv"] = rows

    # Table 3: the whole genome, one thread and every core, three annotation conditions.
    n_wgs = 4_048_342
    if all(res.wall(part, f"hg002_{c}_{t}") for _, _, c, part in CONDITIONS for t in ("t1", "all")) and all((f"wgs_{k}", str(n_wgs)) in base for k, _, _, _ in CONDITIONS):
        rows = [["condition", "label", "fastvep_1t_sec", "fastvep_mt_sec", "vep_sec", "fastvep_1t_vps", "fastvep_mt_vps", "vep_vps", "speedup"]]
        for key, label, cond, part in CONDITIONS:
            one, many = res.wall(part, f"hg002_{cond}_t1"), res.wall(part, f"hg002_{cond}_all")
            v = base[(f"wgs_{key}", str(n_wgs))]
            rows.append([key, label, fmt(one), fmt(many), f"{v:.0f}", str(round(n_wgs / one)), str(round(n_wgs / many)), str(round(n_wgs / v)), speedup(v, one)])
        out["vep_genomewide.csv"] = rows

    # Table 5: three output formats, 10,000 chr22 variants, one thread.
    if all(res.wall("formats", f"chr22_10000_{f}_t1") for f in ("vcf", "tab", "json")) and "human" in res.gene_models:
        t = res.gene_models["human"]["transcripts"]
        rows = [["format", "variants", "transcripts", "time_sec", "output_size_bytes"]]
        for f in ("vcf", "tab", "json"):
            rows.append([f, "10000", t, fmt(res.wall("formats", f"chr22_10000_{f}_t1")), str(res.out_bytes("formats", f"chr22_10000_{f}_t1"))])
        out["format_comparison.csv"] = rows

    # Table 7: five organisms, every core, --hgvs, from run_benchmark.sh.
    if res.multi and res.gene_models:
        names = {"yeast": "yeast_full", "drosophila": "drosophila_full", "arabidopsis": "arabidopsis_full", "mouse": "mouse_full", "human": "human_full"}
        comp = [["organism", "assembly", "source", "transcripts", "variants", "time_sec", "variants_per_sec"]]
        scal = [["organism", "variants", "fastvep_time_sec", "fastvep_vps"]]
        for key, label, assembly, source in ORGANISMS:
            r = res.multi.get(names[key])
            if r is None or key not in res.gene_models:
                continue
            n, t = int(r["variants"]), float(r["time_seconds"])
            comp.append([label, assembly, source, res.gene_models[key]["transcripts"], str(n), fmt(t), str(round(n / t))])
            scal.append([label, str(n), fmt(t), str(round(n / t))])
        if len(comp) == 6:
            out["organism_comparison.csv"], out["scaling.csv"] = comp, scal
            h = res.multi["human_full"]
            out["giab_benchmark.csv"] = [["dataset", "build", "source", "variants", "transcripts", "fastvep_time_sec"],
                                         ["HG002_full_wgs", "GRCh38", "GIAB_HG002", h["variants"], res.gene_models["human"]["transcripts"], fmt(float(h["time_seconds"]))]]

    # Table 4: what the binary costs. Start-up is a run of one variant, which is the time to read the transcript cache.
    needed = [("startup", "yeast_1_variant_t1"), ("startup", "human_1_variant_t1"), ("startup", "human_1_variant_t8"), ("memory", "chr1_100k_t8")]
    if all(res.wall(p, c) for p, c in needed) and res.env.get("binary_bytes"):
        rows = [["metric", "value", "unit"],
                ["binary_size", f"{int(res.env['binary_bytes']) / 1e6:.1f}", "MB"],
                ["cache_load_yeast_7k_1_thread", fmt(res.wall("startup", "yeast_1_variant_t1")), "sec"],
                ["cache_load_human_509k_1_thread", fmt(res.wall("startup", "human_1_variant_t1")), "sec"],
                ["cache_load_human_509k_8_threads", fmt(res.wall("startup", "human_1_variant_t8")), "sec"],
                ["peak_memory_100k_variants", f"{res.rss('memory', 'chr1_100k_t8') / 1024:.1f}", "GiB"],
                ["peak_memory_100k_variants_no_fasta_consequence_only", f"{res.rss('memory', 'chr1_100k_nofasta_conseq_t8') / 1024:.1f}", "GiB"],
                ["binary_bytes", res.env["binary_bytes"], "bytes"],
                ["external_dependencies", "0", "count"]]
        # Ensembl VEP's side of the same table, from run_vep_side.py: start-up is a run of one variant, which includes starting the container
        vs = res.vep / "vep_summary.tsv"
        venv = res.vep / "vep_environment.tsv"
        if vs.exists() and venv.exists():
            cells = {(r["part"], r["cell"]): r for r in read_rows(vs)}
            vmeta = dict(line.rstrip("\n").split("\t", 1) for line in venv.open())
            rows += [["vep_image_size_unpacked", f"{int(vmeta['image_unpacked_kib']) * 1024 / 1e9:.2f}", "GB"],
                     ["vep_image_unpacked_bytes", str(int(vmeta['image_unpacked_kib']) * 1024), "bytes"],
                     ["vep_image_size_pulled", f"{int(vmeta['image_bytes']) / 1e9:.2f}", "GB"],
                     ["vep_startup_yeast_1_variant", fmt(float(cells[("startup", "yeast_1_variant")]["median_wall_s"])), "sec"],
                     ["vep_startup_human_1_variant", fmt(float(cells[("startup", "human_1_variant")]["median_wall_s"])), "sec"],
                     ["vep_peak_memory_100k_variants", f"{float(cells[('memory', 'chr1_100k')]['median_peak_rss_mib']) / 1024:.2f}", "GiB"]]
        out["resource_usage.csv"] = rows

    # Table 8: the consequence calls of the whole HG002 genome, each under its most severe term.
    ms = res.dir / "consequence_most_severe.csv"
    if ms.exists():
        out["consequence_distribution.csv"] = [row for row in csv.reader(ms.open())]

    for name, sub in COPIED.items():
        if (res.dir / sub / name).exists():
            out[name] = [row for row in csv.reader((res.dir / sub / name).open())]
    return out


def provenance(res: Results) -> list[list[str]]:
    keys = ["binary_version", "built_from", "binary_sha256", "system", "macos_product_version", "cpu", "logical_cores", "rustc", "runs_per_cell", "started_utc", "finished_utc"]
    rows = [["key", "value"]] + [[k, res.env[k]] for k in keys if k in res.env]
    venv = res.vep / "vep_environment.tsv"
    if venv.exists():
        vmeta = dict(line.rstrip("\n").split("\t", 1) for line in venv.open())
        rows += [["vep_baseline_file", "vep_baseline.csv"],
                 ["vep_baseline_what", "Ensembl VEP 116.3, Docker " + vmeta["image"] + ", --gff against a tabix-indexed GFF3, one thread, --vcf --hgvs --symbol --canonical; --custom ClinVar and gnomAD v4.1 joint for the annotation conditions"],
                 ["vep_baseline_measured", vmeta["started_utc"][:10]],
                 ["vep_baseline_machine", res.env.get("system", "") + ", " + res.env.get("cpu", "")],
                 ["vep_baseline_wgs", "sum over chromosomes 1 to 22 of one run of each"]]
    if res.env_later:
        rows += [["with_vep_started_utc", res.env_later.get("started_utc", "")], ["with_vep_finished_utc", res.env_later.get("finished_utc", "")],
                 ["with_vep_parts", res.env_later.get("parts", "")]]
    rows += [["results_directory", str(res.dir.relative_to(HERE.parent)) if res.dir.is_relative_to(HERE.parent) else str(res.dir)]]
    return rows


def write_csv(path: Path, rows: list[list[str]]) -> None:
    with path.open("w", newline="") as f:
        csv.writer(f, lineterminator="\n").writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True, type=Path)
    ap.add_argument("--check", action="store_true", help="write nothing; exit 1 if a file differs from what the results say")
    args = ap.parse_args()
    res = Results(args.results.resolve())
    files = build(res)
    files["provenance.csv"] = provenance(res)
    bad = []
    for name, rows in sorted(files.items()):
        path = DATA / name
        new = path.with_suffix(".tmp")
        write_csv(new, rows)
        same = path.exists() and new.read_bytes() == path.read_bytes()
        if args.check or same:
            new.unlink()
        else:
            new.replace(path)
        print(f"{'ok     ' if same else ('DIFFERS' if args.check else 'written')} {name}")
        if args.check and not same:
            bad.append(name)
    for name in ("vep_comparison.csv", "vep_genomewide.csv", "format_comparison.csv", "organism_comparison.csv", "scaling.csv", "giab_benchmark.csv",
                 "resource_usage.csv", "consequence_distribution.csv", *COPIED):
        if name not in files:
            print(f"left as it is (no input in the results directory) {name}")
    if bad:
        sys.exit(f"out of step with {args.results}: {', '.join(bad)}")


if __name__ == "__main__":
    main()
