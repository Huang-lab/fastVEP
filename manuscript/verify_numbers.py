#!/usr/bin/env python3
"""Check that the numbers in `manuscript/manuscript.md` (and in the performance section of `README.md`) are the numbers in `manuscript/data/`.

The text is written by hand and the CSV files are written by `build_data.py` from the runs that measured them, so this is the bridge: it
recomputes each quoted number from the CSV files, finds the sentence or table row that quotes it, and fails when that line does not
contain the recomputed value. A line is found by a phrase that occurs in no other line, so an edit that moves or rewords a claim shows up
as an anchor that no longer matches and not as a silent pass.

    manuscript/verify_numbers.py        run every check, print the failures, exit 1 if there are any
    manuscript/verify_numbers.py -v     also print every check that passed

The counts of the comparison with Ensembl VEP 116.3 are written by `validation/measure_vep116.py` (`manuscript/data/vep116_provenance.md` says how), and
the text of `docs/VEP_DIVERGENCE.md` quotes them as well, so each is looked for in its row of that document too.

Not checked: a number that is not a function of the CSV files (a date, a count read from the repository's own documents that is not one
of those above), and the wording of any sentence that passes.

Standard library only.
"""

from __future__ import annotations

import argparse
import collections
import csv
import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
RESULTS = Path(os.environ.get("FASTVEP_RESULTS", HERE.parent / "benchmarks" / "results" / "v0.5.0"))  # the results directory that `build_data.py` was given


def logical_lines(text: str) -> list[str]:
    """A Markdown file as one line to a paragraph, a table row, a heading or a list item: a line that goes on with a paragraph is joined to the one before it.

    The manuscript has a sentence to a line, and a number is looked for in the paragraph that carries it and not in the sentence.
    """
    out: list[str] = []
    joinable = False  # the last line is prose or a list item, which the next line may continue
    for line in text.splitlines():
        stripped = line.strip()
        starts_block = (not stripped or line.startswith(("#", "|", "```", "---")) or re.match(r"(- |\* |\d+\. )", line) is not None)
        if starts_block:
            out.append(line)
            joinable = bool(stripped) and not line.startswith(("#", "|", "```", "---"))
        elif joinable and out:
            out[-1] += " " + stripped
        else:
            out.append(line)
            joinable = True
    return out


PAPER = logical_lines((HERE / "manuscript.md").read_text())
README = logical_lines((HERE.parent / "README.md").read_text())

failures: list[str] = []
passed: list[str] = []


def table(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="") as f:
        return list(csv.DictReader(f))


def table_tsv(relative: str) -> list[dict[str, str]]:
    """A tab-separated file of benchmarks/results/<build>/ (the comparison of builds is not data of the manuscript, but the text quotes it)."""
    with (RESULTS / relative).open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def kv(name: str) -> dict[str, str]:
    return {r[next(iter(r))]: r["value"] for r in table(name)}


def fact(lines: list[str], where: str, anchor: str, *values: str | tuple[str, ...]) -> None:
    """The one line that contains `anchor` must contain each value (a tuple lists acceptable spellings of one)."""
    hits = [i for i, line in enumerate(lines) if anchor in line]
    if len(hits) != 1:
        failures.append(f"{where}: {'no line' if not hits else f'{len(hits)} lines'} contain the anchor {anchor!r}")
        return
    line = lines[hits[0]]
    for v in values:
        alts = v if isinstance(v, tuple) else (v,)
        if any(a in line for a in alts):
            passed.append(f"{where}: {anchor[:48]!r}: {alts[0]}")
        else:
            failures.append(f"{where}, line {hits[0] + 1}: {anchor[:48]!r} does not contain {alts[0]!r}")


def block(lines: list[str], caption: str) -> list[list[str]]:
    """The rows (cells, stripped) of the Markdown table that follows the line starting with `caption`."""
    starts = [i for i, line in enumerate(lines) if line.startswith(caption)]
    if len(starts) != 1:
        failures.append(f"{len(starts)} lines start with {caption!r}")
        return []
    rows, i = [], starts[0] + 1
    while i < len(lines) and not lines[i].startswith("|"):
        i += 1
    while i < len(lines) and lines[i].startswith("|"):
        rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
        i += 1
    return rows[2:]  # the header and the rule


def same(where: str, got: list[list[str]], want: list[list[str]]) -> None:
    if len(got) != len(want):
        failures.append(f"{where}: {len(got)} rows in the text, {len(want)} in the data")
        return
    for g, w in zip(got, want):
        if g != w:
            failures.append(f"{where}: row {g[0]!r} is {g} and the data say {w}")
        else:
            passed.append(f"{where}: row {g[0]!r}")


def sp(x: float) -> str:
    return f"{x:.2f}x" if x < 1 else (f"{x:.1f}x" if x < 10 else f"{x:.0f}x")


# ---- Table 2 and 3: the head-to-head ----------------------------------------------------------------------------------
h2h, gw = table("vep_comparison.csv"), table("vep_genomewide.csv")
same("Table 2", block(PAPER, "**Table 2."), [[f"{int(r['variants']):,}", f"{float(r['fastvep_sec']):.2f}s", f"{float(r['vep_sec']):.2f}s",
                                                  sp(float(r['vep_sec']) / float(r['fastvep_sec'])) + (" (VEP faster)" if float(r['vep_sec']) < float(r['fastvep_sec']) else "")] for r in h2h])
t3 = []
for r in gw:
    vep = float(r["vep_sec"])
    t3.append([r["label"].replace("Consequence", "Consequence only") if r["condition"] == "conseq" else r["label"], f"{float(r['fastvep_1t_sec']):.1f}s", f"{float(r['fastvep_mt_sec']):.1f}s",
               f"{vep:,.0f}s" + (f" (~{vep / 60:.0f} min)" if r["condition"] == "conseq" else ""), f"{vep / float(r['fastvep_1t_sec']):.0f}x"])
same("Table 3", block(PAPER, "**Table 3."), t3)
lo, hi = min(float(r["speedup"]) for r in gw), max(float(r["speedup"]) for r in gw)
fact(PAPER, "abstract", "fastVEP annotates the complete GIAB HG002", f"{lo:.0f}-{hi:.0f}x")
fact(PAPER, "3.2", "On the complete GIAB HG002 WGS, fastVEP is", f"{lo:.0f}-{hi:.0f}x")

# the numbers that the text derives from Table 3 and from the runs behind it
N_WGS = 4_048_342
plain = gw[0]
one, many, vep_plain = float(plain["fastvep_1t_sec"]), float(plain["fastvep_mt_sec"]), float(plain["vep_sec"])
raw_later = [r for r in table_tsv("with_vep/raw.tsv") if int(r["round"]) >= 1]
rounds = lambda cell, part="wgs": [float(r["wall_s"]) for r in raw_later if r["part"] == part and r["cell"] == cell]
with (RESULTS / "with_vep" / "vep_genome_sums.csv").open(newline="") as f:
    vsum = {r["condition"]: r for r in csv.DictReader(f)}
cpu_share = 100 * sum(float(r["cpu_s"]) for r in vsum.values()) / sum(float(r["wall_s"]) for r in vsum.values())
with (RESULTS / "with_vep" / "custom_matches_chr22.tsv").open(newline="") as f:
    cm = {r["source"]: r for r in csv.DictReader(f, delimiter="\t")}
assert all(r["only_vep"] == "0" == r["only_fastvep"] and r["vep_found"] == r["fastvep_found"] for r in cm.values()), "the two tools found different records in the same ClinVar and gnomAD files"
fact(PAPER, "abstract", "fastVEP annotates the complete GIAB HG002", f"in {many:.0f} seconds on ten cores")
fact(PAPER, "abstract", "In single-threaded head-to-head benchmarks", f"{lo:.0f}-{hi:.0f}x", f"approximately {vep_plain / 60:.0f} minutes")
fact(PAPER, "intro", "one VEP process annotates hundreds of variants", f"{N_WGS / vep_plain:,.0f} per second", f"about {vep_plain / 3600:.1f} hours")
fact(PAPER, "intro", "annotates a complete human whole genome", f"{lo:.0f}-{hi:.0f}x")
fact(PAPER, "3.2", "The 5 consequence-only runs behind the medians", f"from {min(rounds('hg002_plain_t1')):.0f} to {max(rounds('hg002_plain_t1')):.0f} seconds on one thread",
     f"from {min(rounds('hg002_plain_all')):.1f} to {max(rounds('hg002_plain_all')):.1f} seconds on ten threads")
fact(PAPER, "3.2", "With its default multi-threaded pipeline", f"{many:.0f} seconds", f"{N_WGS / many:,.0f} variants per second", f"{one / many:.1f}x its single-threaded speed")
fact(PAPER, "3.2", "These whole-genome VEP times are", f"{float(vsum['plain']['cpu_s']):,.0f} s of the {vep_plain:,.0f} s of wall time")
fact(PAPER, "3.2", "In the present comparison both tools consume", f"{lo:.0f}-{hi:.0f}x")
fact(PAPER, "3.2", "Supplementary annotations were supplied to both tools", f"{int(cm['ClinVar']['both']):,} ClinVar records", f"{int(cm['gnomAD']['both']):,} gnomAD records")
fact(PAPER, "4.1", "Benchmarks on the GIAB HG002 gold-standard dataset", f"{lo:.0f}-{hi:.0f}x", f"about {many:.0f} seconds", f"~{N_WGS / many:,.0f} variants per second")
fact(PAPER, "4.1", "The fixed binary-cache start-up", f"{N_WGS / one:,.0f} variants per second")
fact(PAPER, "4.1", "The practical implications are substantial", f"roughly {500_000 / (N_WGS / one):.0f} seconds", f"about {500_000 / (N_WGS / many):.1f} seconds multi-threaded",
     f"approximately {500_000 / (N_WGS / vep_plain) / 60:.0f} minutes", f"{N_WGS / vep_plain:,.0f} variants per second")
fact(PAPER, "Figure 2", "(A) Wall-clock annotation time (log-log)", f"{vep_plain:,.0f} s", f"fastVEP single-thread, {one:.0f} s", f"fastVEP multi-thread, {many:.0f} s")
fact(PAPER, "Figure 2", "(B) Whole-genome single-threaded head-to-head", f"{lo:.0f}-{hi:.0f}x", f"({min(float(r['fastvep_mt_sec']) for r in gw):.0f}-{max(float(r['fastvep_mt_sec']) for r in gw):.0f} s)")
fact(PAPER, "5.2", "The time of a VEP run is the wall-clock time of", f"{cpu_share:.1f}% of the wall time")
fact(PAPER, "5.2", "Supplementary-annotation conditions gave both tools", f"{int(cm['ClinVar']['both']):,} ClinVar records", f"{int(cm['gnomAD']['both']):,} gnomAD records")
first = [r for r in table_tsv("summary.tsv") if r["cell"] == "hg002_plain_t1"][0]
later = [r for r in table_tsv("with_vep/summary.tsv") if r["cell"] == "hg002_plain_t1"][0]
fact(PAPER, "5.2", "The parts of fastVEP's harness that are compared with VEP", f"{float(first['median_wall_s']):.1f} s in an earlier session and {float(later['median_wall_s']):.1f} s in this one")
spot = RESULTS / "with_vep" / "vep_115_1" / "vep_summary.tsv"
if spot.exists():
    s115 = {r["cell"]: float(r["median_wall_s"]) for r in table_tsv("with_vep/vep_115_1/vep_summary.tsv") if r["part"] == "scaling"}
    cells = [("1000", "chr22_1000"), ("5000", "chr22_5000"), ("10000", "chr22_10000"), ("50284", "chr22_all")]
    v116 = {r["variants"]: float(r["vep_sec"]) for r in h2h}
    ratios = [v116[n] / s115[c] for n, c in cells]
    fact(PAPER, "3.2", "Release 116.3 is slower than release 115.1", ", ".join(f"{s115[c]:.2f}" for _, c in cells), ", ".join(f"{v116[n]:.2f}" for n, _ in cells),
         f"{min(ratios):.1f} to {max(ratios):.1f} times as long")
    fact(README, "README", "Release 116.3 takes", f"{min(ratios):.1f} to {max(ratios):.1f} times as long", ", ".join(f"{s115[c]:.2f}" for _, c in cells), ", ".join(f"{v116[n]:.2f}" for n, _ in cells))

# ---- Table 5: formats ---------------------------------------------------------------------------------------------------
fm = table("format_comparison.csv")
want5 = [[{"vcf": "VCF (CSQ, 49 fields)", "tab": "Tab-delimited", "json": "JSON"}[r["format"]], f"{float(r['time_sec']):.2f}", f"{int(r['output_size_bytes']) / 1024:,.0f} KB"] for r in fm]
same("Table 5", block(PAPER, "**Table 5."), want5)
start1 = float(next(r for r in table_tsv("with_vep/summary.tsv") if r["cell"] == "human_1_variant_t1")["median_wall_s"])
adds = {r["format"]: float(r["time_sec"]) - start1 for r in fm}
fact(PAPER, "3.5", "All three output formats perform similarly", f"the {start1:.2f}-second binary-cache start-up", f"add {adds['vcf']:.2f} seconds as VCF", f"{adds['tab']:.2f} as tab-delimited", f"{adds['json']:.2f} as JSON")

# ---- Table 7 and Figure 5: five organisms --------------------------------------------------------------------------------
org = table("organism_comparison.csv")
rows7 = block(PAPER, "**Table 7.")
same("Table 7 (transcripts, variants, time, throughput)", [[r[2], r[3], r[5], r[6]] for r in rows7],
     [[f"{int(o['transcripts']):,}", f"{int(o['variants']):,}", f"{float(o['time_sec']):.2f}", f"{int(o['variants_per_sec']):,}"] for o in org])
vps = [int(o["variants_per_sec"]) for o in org]
fact(PAPER, "3.7", "fastVEP annotates complete variant datasets across all five organisms", f"{min(vps):,}", f"{max(vps):,}")
tx = [int(o["transcripts"]) for o in org]
fact(PAPER, "Figure 5", "**Figure 5.", f"{min(vps) / 1000:.0f}K", f"{max(vps) / 1000:.0f}K", f"{max(tx) / min(tx):.0f}x")

# ---- Table 8: consequences ------------------------------------------------------------------------------------------------
cons = table("consequence_distribution.csv")
calls = sum(int(r["count"]) for r in cons)
fr = lambda c: ("<0.001%" if c / calls * 100 < 0.001 else (f"{c / calls * 100:.3f}%" if c / calls * 100 < 0.01 else (f"{c / calls * 100:.2f}%" if c / calls * 100 < 0.1 else f"{c / calls * 100:.1f}%")))
same("Table 8", block(PAPER, "**Table 8."), [[r["consequence"], f"{int(r['count']):,}", fr(int(r["count"]))] for r in cons])
fact(PAPER, "3.8", "The annotation produced", f"{calls:,}", f"{len(cons)} distinct SO terms")
fact(PAPER, "Figure 4", "**Figure 4.", f"{calls / 1e6:.1f} million")

# ---- Table 1, Table 9 and the code ------------------------------------------------------------------------------------------
m = kv("code_metrics.csv")
crates = ["core", "genome", "cache", "consequence", "hgvs", "io", "filter", "sa", "annotate", "classification", "cli", "web"]
got1 = block(PAPER, "**Table 1.")
same("Table 1 (tests)", [[r[0], r[1]] for r in got1[:-1]], [[f"`fastvep-{c}`", m[f"tests_fastvep-{c}"]] for c in crates])
if got1:
    same("Table 1 (total)", [[got1[-1][1]]], [[f"**{int(m['tests_total']):,}**"]])
fact(PAPER, "2.1", "The twelve crates hold", f"{int(m['total_source_lines']):,}")
fact(PAPER, "2.4", "implemented in `fastvep-consequence`", f"{int(m['loc_fastvep-consequence']):,} lines")
fact(PAPER, "2.2", "The predictors produce", m["so_consequence_types_emitted"], m["so_consequence_types"])
fact(PAPER, "2.6", "a `CSQ` INFO field with", f"{m['csq_output_fields']} pipe-delimited")
fact(PAPER, "2.8", "a single-page application", f"{int(m['web_gui_lines']):,} lines")
fact(PAPER, "3.9", "fastVEP's correctness was validated through", f"{int(m['tests_total']):,}")
so = table("so_terms.csv")
rows_s1 = block(PAPER, "### Supplementary Table S1.")
same("Table S1", [r[:4] for r in rows_s1], [[t["rank"], t["so_term"], t["impact"], t["emitted_by_a_predictor"]] for t in so])

# ---- resource usage -------------------------------------------------------------------------------------------------------
res = kv("resource_usage.csv")
# peak memory by thread count, as the operating system reports it (MiB), in the session of the VEP runs
mem = {r["cell"]: float(r["median_peak_rss_mib"]) for r in table_tsv("with_vep/summary.tsv") if r["part"] == "memory"}
fact(PAPER, "Table 4", "| Binary/install size |", f"{float(res['binary_size']):.1f} MB")
fact(PAPER, "Table 4", "| Startup (7K transcripts, yeast) |", f"{float(res['cache_load_yeast_7k_1_thread']):.2f} sec")
fact(PAPER, "Table 4", "| Startup (510K transcripts, human) |", f"{float(res['cache_load_human_509k_1_thread']):.2f} sec", f"{float(res['cache_load_human_509k_8_threads']):.2f} sec")
fact(PAPER, "Table 4", "| Peak memory (100,000 variants, human) |", f"{float(res['peak_memory_100k_variants']):.1f} GiB")
# Ensembl VEP's side of Table 4 (run_vep_side.py): the image, a run of one variant, and the peak resident set of 100,000 variants
vs = {r["cell"]: r for r in table_tsv("with_vep/vep_summary.tsv")}
vep_yeast, vep_human = float(res["vep_startup_yeast_1_variant"]), float(res["vep_startup_human_1_variant"])
vep_mem = float(res["vep_peak_memory_100k_variants"])
fast_t1, fast_t8 = float(res["cache_load_human_509k_1_thread"]), float(res["cache_load_human_509k_8_threads"])
fast_yeast = float(res["cache_load_yeast_7k_1_thread"])
fact(PAPER, "Table 4", "| Binary/install size |", f"{float(res['vep_image_size_unpacked']):.2f} GB (Docker image, unpacked", f"{float(res['vep_image_size_pulled']):.2f} GB to pull",
     f"{float(res['binary_size']):.1f} MB", f"{int(res['vep_image_unpacked_bytes']) / int(res['binary_bytes']):.0f}x smaller")
fact(PAPER, "Table 4", "| Startup (7K transcripts, yeast) |", f"{vep_yeast:.2f} sec", f"{fast_yeast:.2f} sec", f"{vep_yeast / fast_yeast:.0f}x")
fact(PAPER, "Table 4", "| Startup (510K transcripts, human) |", f"{vep_human:.2f} sec", f"{fast_t1:.2f} sec (1 thread)", f"{fast_t8:.2f} sec (8 threads)", f"VEP {fast_t1 / vep_human:.1f}x faster on one thread")
fact(PAPER, "Table 4", "| Peak memory (100,000 variants, human) |", f"{vep_mem:.2f} GiB", f"VEP {mem['chr1_100k_t8'] / 1024 / vep_mem:.1f}x less")
vep_cpu = {r["cell"]: float(r["median_cpu_s"]) for r in vs.values() if r["part"] == "startup"}
fact(PAPER, "3.4", "Ensembl VEP's run of one variant takes", f"{vep_human:.2f} seconds for either gene model",
     f"{vep_cpu['yeast_1_variant']:.2f} to {vep_cpu['human_1_variant']:.2f} seconds is the CPU time")
assert abs(vep_human - vep_yeast) < 0.005, "the text says that a run of one variant takes the same time against either gene model"
fact(PAPER, "3.4", "Ensembl VEP's process peaks at", f"{vep_mem:.2f} GiB", f"it needs {mem['chr1_100k_t8'] / 1024 / vep_mem:.1f} times less memory")
fact(PAPER, "Table 10", "| Binary size |", f"{float(res['binary_size']):.1f} MB", f"{float(res['vep_image_size_unpacked']):.2f} GB (Docker image, unpacked)")
fact(PAPER, "3.4", "The cache holds the whole gene model with its sequences in memory", f"{mem['chr1_100k_t1'] / 1024:.2f} GiB on one", f"{mem['chr1_100k_t8'] / 1024:.2f} GiB on eight",
     f"{mem['chr1_100k_tall'] / 1024:.2f} GiB on ten", f"{mem['chr1_100k_nofasta_conseq_t8'] / 1024:.1f} GiB")
fact(README, "README", "The peak memory of a run is the gene model", f"{mem['chr1_100k_t1']:,.0f} MiB on one thread", f"{mem['chr1_100k_t8']:,.0f} MiB on eight",
     f"{mem['chr1_100k_nofasta_conseq_t1']:,.0f} MiB and {mem['chr1_100k_nofasta_conseq_t8']:,.0f} MiB")

# ---- Table 6: fastSA v1 against v2 -------------------------------------------------------------------------------------------
sa, saa = table("sa_formats.csv"), table("sa_formats_annotation.csv")
want6 = []
for r in sa:
    label = r["source"].replace(" (whole database)", "")
    want6.append([label, f"{int(r['records']) / 1e6:.2f}M", f"{int(r['v1_bytes']) / 1e6:.1f} MB", f"{int(r['v2_bytes']) / 1e6:.1f} MB",
                  f"{int(r['v2_bytes']) / int(r['v1_bytes']):.2f}x", f"{float(r['build_v1_s']):.1f}s", f"{float(r['build_v2_s']):.1f}s"])
same("Table 6", block(PAPER, "**Table 6."), want6)
A = {(r["config"], r["format"]): r for r in saa}
assert all(A[(c, "osa")]["md5"] == A[(c, "osa2")]["md5"] != "VARIES" for c in ("clinvar", "clinvar+revel")), "the v1 and v2 outputs differ: the text says they are identical"
fact(PAPER, "3.6", "The annotation output is byte-identical between the two formats", f"{int(A[('clinvar', 'osa')]['variants']):,}",
     f"{float(A[('clinvar', 'osa')]['median_s']):.1f} s and {float(A[('clinvar', 'osa2')]['median_s']):.1f} s",
     f"{float(A[('clinvar+revel', 'osa')]['median_s']):.1f} s and {float(A[('clinvar+revel', 'osa2')]['median_s']):.1f} s",
     f"{int(A[('clinvar+revel', 'osa')]['records_with_revel']):,}")

# ---- tests that ran ------------------------------------------------------------------------------------------------------------
ex = kv("tests_executed.csv")
# 1,733 `#[test]` attributes: 1,731 ran and passed on macOS, one is ignored, and one is compiled on Linux only
assert int(ex["passed"]) + int(ex["ignored"]) + int(ex["linux_only"]) == int(m["tests_total"]), "the tests that ran do not add up to the tests that are counted"
fact(PAPER, "3.9", "fastVEP's correctness was validated through", f"{int(ex['passed']):,} of them run and pass", "one is ignored", "Linux only")
fact(PAPER, "5.4", "The release has", f"{int(m['tests_total']):,}", f"{int(ex['passed']):,} ran and passed", "one was ignored", "Linux only")

# ---- accuracy against Ensembl VEP 116.3 -------------------------------------------------------------------------------------------
# The counts are written by validation/measure_vep116.py, and docs/VEP_DIVERGENCE.md, Part 3, quotes them: each is looked for in its own row of that document.
DIV = logical_lines((HERE.parent / "docs" / "VEP_DIVERGENCE.md").read_text())
agree = {r["input"]: r for r in table("vep116_agreement.csv")}
for r in agree.values():
    fact(DIV, "VEP_DIVERGENCE.md", f"| {int(r['records']):,} | {int(r['paired_rows']):,} | {int(r['cells_compared']):,} |", f"{int(r['cells_differ']):,}")
cols = {r["column"]: r for r in table("vep116_hg002_columns.csv")}
fact(DIV, "VEP_DIVERGENCE.md", "On HG002 the", f"{sum(int(c['cells_differ']) for c in cols.values()):,} cells are in four columns", f"`HGVS_OFFSET` {int(cols['HGVS_OFFSET']['cells_differ']):,}",
     f"`HGVSc` {int(cols['HGVSc']['cells_differ']):,}", f"`HGVSp` {cols['HGVSp']['cells_differ']}", f"`HGVSg` {cols['HGVSg']['cells_differ']}")
fact(DIV, "VEP_DIVERGENCE.md", "| `HGVSc` and `HGVS_OFFSET` of an allele in a record with several ALTs |", f"{int(cols['HGVSc']['divergence_11_several_alts']):,} and {int(cols['HGVS_OFFSET']['cells_differ']):,}")
fact(DIV, "VEP_DIVERGENCE.md", "| `HGVSg` of an allele in a record with several ALTs |", cols["HGVSg"]["cells_differ"])
fact(DIV, "VEP_DIVERGENCE.md", "| `HGVSc` of an intronic base after the terminator |", cols["HGVSc"]["divergence_4_terminator_intron"])
fact(DIV, "VEP_DIVERGENCE.md", "| `HGVSp` | 17 (divergence 17) |", cols["HGVSp"]["divergence_17_three_prime_rule"])
hg, c5, c5h, world, edges, lone, ex_set, inframe = (agree[k] for k in ("hg002_1in200", "clinvar_1in5", "clinvar_1in5_hgvs", "world", "edges", "lone_indels", "example", "inframe_deletions"))
assert all(int(c["unclassified"]) == 0 for c in cols.values()), "a cell of the HG002 sample differs for a cause that the document does not record"
share = lambda r: 100.0 * (int(r["cells_compared"]) - int(r["cells_differ"])) / int(r["cells_compared"])
fact(PAPER, "3.3", "On a systematic 1-in-200 sample of the GIAB HG002 call set", f"{int(hg['records']):,}", f"{int(hg['paired_rows']):,}", f"{int(hg['cells_compared']):,}", f"{int(hg['cells_differ']):,}",
     f"{100 - share(hg):.3f}%", f"`HGVS_OFFSET` {int(cols['HGVS_OFFSET']['cells_differ']):,}", f"`HGVSc` {int(cols['HGVSc']['cells_differ']):,}", f"`HGVSg` {cols['HGVSg']['cells_differ']}", f"`HGVSp` {cols['HGVSp']['cells_differ']}")
fact(PAPER, "3.3", "On a systematic 1-in-5 sample of the ClinVar 2-star+ call set", f"{int(c5['records']):,}", f"{int(c5['paired_rows']):,}", f"{int(c5['cells_differ'])} of {int(c5['cells_compared']):,}",
     f"{100 - share(c5):.5f}%", f"{int(c5h['cells_differ']):,} of {int(c5h['cells_compared']):,}", f"{100 - share(c5h):.4f}%")
# the 53 cells of the in-frame deletions, by variant: the text says which are the NT5C2 deletion, which lose the initiator and which are on incomplete transcripts
by_variant: collections.Counter = collections.Counter()
features = collections.defaultdict(set)
for path in sorted((RESULTS / "vep116_diff" / "inframe").glob("diff_*.tsv")):
    with path.open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            by_variant[row["Uploaded_variation"]] += 1
            features[row["Uploaded_variation"]].add(row["Feature"])
assert sum(by_variant.values()) == int(inframe["cells_differ"]), "the differing rows of the in-frame set are not the cells that differ"
nt5c2, initiator = by_variant["10_103089678_TTCCTCC/T"], by_variant["1_110604009_AATGATGTTC/A"]
others = [v for v in by_variant if v not in ("10_103089678_TTCCTCC/T", "1_110604009_AATGATGTTC/A")]
fact(PAPER, "3.3", "On the 173 variants of VEP's own example set", f"{int(ex_set['cells_compared']):,}", f"{int(inframe['cells_differ'])} of {int(inframe['cells_compared']):,}",
     f"{100 * int(inframe['cells_differ']) / int(inframe['cells_compared']):.3f}%", f"{nt5c2} are the NT5C2 deletion", f"{initiator} are a change that loses the initiator",
     f"{sum(by_variant[v] for v in others)} are on {({1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five'})[len(set().union(*(features[v] for v in others)))]} transcripts")
assert int(ex_set["cells_differ"]) == 0, "the text says that no cell of the VEP example set differs"
fact(PAPER, "3.3", "A real call set holds few of the shapes", f"{int(world['records']):,}", f"{int(world['paired_rows']):,}", f"{int(world['cells_differ'])} of {int(world['cells_compared']):,}", f"{100 - share(world):.2f}%",
     f"{int(edges['records']):,}", f"{int(edges['paired_rows']):,}", f"{int(edges['cells_differ'])} of {int(edges['cells_compared']):,}", f"{100 - share(edges):.3f}%",
     f"{int(lone['records']):,}", f"{int(lone['cells_differ'])} of {int(lone['cells_compared']):,}", f"{100 - share(lone):.3f}%")
fact(PAPER, "abstract", "whole-genome call set (", f"{int(hg['paired_rows']):,}", f"{share(hg):.2f}%", f"{int(hg['cells_compared']) / 1e6:.1f} million", f"{int(hg['cells_differ']):,}",
     f"{share(c5h):.4f}%", f"{int(c5h['cells_compared']) / 1e6:.1f} million")
fact(PAPER, "intro", "agreeing with Ensembl VEP 116.3", f"{share(hg):.2f}%", f"{share(c5h):.4f}%")
fact(PAPER, "4.2", "call set agrees on every paired row in the columns that carry", f"{share(hg):.2f}%", f"{int(hg['cells_differ']):,}", f"{share(c5h):.4f}%",
     f"{share(world):.2f}%", f"{share(edges):.2f}%", f"{share(lone):.2f}%")
fact(PAPER, "Table 9", "| **Agreement with Ensembl VEP 116.3** |", f"{share(hg):.2f}%", f"{share(c5h):.4f}%")
fact(PAPER, "Table 10", "| Agreement with Ensembl VEP 116.3 |", f"{share(hg):.2f}%", f"{share(c5h):.4f}%")
ident = [r for r in table_tsv("identity_vs_5048fc5/identity.tsv")]
assert ident and all(r["identical"] == "yes" for r in ident), "the output of 0.5.0 differs from the build of 27 October: the text says it is identical"
fact(PAPER, "3.3", "Ensembl VEP's output does not depend on the build of fastVEP", f"{len(ident)} runs in all")
big = next(r for r in ident if r["input"] == "hg002_full")
fact(PAPER, "3.3", "Ensembl VEP's output does not depend on the build of fastVEP", f"{int(big['records']):,} records", f"{int(big['bytes_v0.5.0-rc.1']) / 1e9:.1f} GB")
fact(PAPER, "Table 9", "| **Golden files of VEP 116.3** |", "26", "20,363")
fact(PAPER, "Figure 3", "**Figure 3.", f"{int(hg['paired_rows']):,}", f"{int(hg['cells_differ']):,}")

fact(PAPER, "4.4", "GFF3 as primary annotation source with binary caching.", f"{float(res['cache_load_yeast_7k_1_thread']):.2f} seconds", f"{float(res['cache_load_human_509k_1_thread']):.2f} seconds", f"{float(res['cache_load_human_509k_8_threads']):.2f} on eight")

# ---- the README quotes the same runs ---------------------------------------------------------------------------------------
plain = gw[0]
fact(README, "README", "| Consequence only |", f"{float(plain['fastvep_1t_sec']):.1f}s", f"{float(plain['vep_sec']):,.0f}s")
for r in h2h:
    fact(README, "README", f"| {int(r['variants']):,} | ", f"{float(r['fastvep_sec']):.2f}s")
for o in org:
    fact(README, "README", f"| {int(o['variants']):,} | {{}}".format("") if False else f"| {int(o['transcripts']):,} | {int(o['variants']):,} |", f"{int(o['variants_per_sec']):,}")

fact(README, "README", "A laptop is a noisy instrument", f"{min(rounds('hg002_plain_t1')):.0f} to {max(rounds('hg002_plain_t1')):.0f} seconds",
     f"{min(rounds('hg002_plain_all')):.1f} to {max(rounds('hg002_plain_all')):.1f}")
fact(README, "README", "With fastVEP's default multi-threading", f"**{many:.1f}s**", f"**{float(gw[-1]['fastvep_mt_sec']):.1f}s**")
fact(README, "README", "Both tools read the same ClinVar and gnomAD files", f"{int(cm['ClinVar']['both']):,} ClinVar records", f"{int(cm['gnomAD']['both']):,} gnomAD records")
for r in gw[1:]:
    fact(README, "README", f"| {r['label']} |", f"{float(r['fastvep_1t_sec']):.1f}s", f"{float(r['vep_sec']):,.0f}s", f"**{float(r['vep_sec']) / float(r['fastvep_1t_sec']):.0f}×**")

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("-v", action="store_true", help="also print every check that passed")
args = ap.parse_args()
if args.v:
    print("\n".join(passed))
print(f"{len(passed)} values found in their lines, {len(failures)} failures")
for f in failures:
    print("FAIL", f)
sys.exit(1 if failures else 0)
