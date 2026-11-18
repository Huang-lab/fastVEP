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
DOC_SV = logical_lines((HERE.parent / "docs" / "STRUCTURAL_VARIANTS.md").read_text())

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
fact(PAPER, "Table 9", "| **Golden files of VEP 116.3** |", "31", "27,055")
fact(PAPER, "Figure 3", "**Figure 3.", f"{int(hg['paired_rows']):,}", f"{int(hg['cells_differ']):,}")

fact(PAPER, "4.4", "GFF3 as primary annotation source with binary caching.", f"{float(res['cache_load_yeast_7k_1_thread']):.2f} seconds", f"{float(res['cache_load_human_509k_1_thread']):.2f} seconds", f"{float(res['cache_load_human_509k_8_threads']):.2f} on eight")

# ---- structural variants (Section 3.11) ----------------------------------------------------------------------------------------
# The three tables are written by hand from sv_agreement.csv (scored by validation/score_strict_key.py), sv_timing.csv
# (benchmarks/run_sv_head_to_head.sh) and small_variant_regression.csv (benchmarks/run_small_variant_regression.sh).
sva = {r["input"]: r for r in table("sv_agreement.csv")}
svt = {(r["input"], r["tool"], r["threads"]): float(r["median_s"]) for r in table("sv_timing.csv")}
svr = table("small_variant_regression.csv")
R1, R2 = "R1_gnomad_sv_v41_chr21_q", "R2_clinvar_chr21_sv"
n_ = lambda x: f"{int(x):,}"
SV_NINE = ["01_snv", "02_mnp_complex", "03_small_indels", "05_symbolic_del_ins", "06_symbolic_dup_inv", "07_cnv_repeat", "08_mobile_elements", "11_multi_allelic", "12_complex_imprecise"]
assert all(sva[k]["strict_f1"] == "1.0000" for k in SV_NINE), "the text says that nine files of the vep-rs set score 1.0000"
assert all(sva[k]["fastvep_rows"] == sva[k]["vep_rows"] for k in SV_NINE), "the text says that the nine files have the rows of Ensembl VEP, row for row"
row_ = lambda label, k: [label, n_(sva[k]["records"]), n_(sva[k]["vep_rows"]), n_(sva[k]["fastvep_rows"]), sva[k]["strict_f1"]]
want11 = [
    row_("Generated world", "world"),
    row_("Generated world, another seed", "world_seed_7"),
    ["vep-rs files 01-03, 05-08, 11, 12 (9 files)", n_(sum(int(sva[k]["records"]) for k in SV_NINE)), n_(sum(int(sva[k]["vep_rows"]) for k in SV_NINE)),
     n_(sum(int(sva[k]["fastvep_rows"]) for k in SV_NINE)), "1.0000"],
    row_("vep-rs file 04 (explicit alleles)", "04_large_explicit"),
    row_("vep-rs file 09 (breakends)", "09_breakends"),
    ["vep-rs files 10 and 13 (`<*>` and VCF 4.5)", n_(int(sva["10_special_alleles"]["records"]) + int(sva["13_vcf45_features"]["records"])),
     n_(int(sva["10_special_alleles"]["vep_rows"]) + int(sva["13_vcf45_features"]["vep_rows"])), n_(int(sva["10_special_alleles"]["fastvep_rows"]) + int(sva["13_vcf45_features"]["fastvep_rows"])),
     f"{sva['10_special_alleles']['strict_f1']} and {sva['13_vcf45_features']['strict_f1']}"],
    row_("ClinVar structural variants, chromosome 21", R2),
    row_("gnomAD-SV v4.1, chromosome 21, 1-in-4 sample", R1),
    row_("the same, the contig spelled `chr21` throughout", R1 + "_one_spelling"),
]
same("Table 11", [r[:5] for r in block(PAPER, "**Table 11.")], want11)
sv_speed = lambda inp, th: svt[(inp, "vep", "1")] / svt[(inp, "fastvep", th)]
want12 = [[name, n_(sva[k]["records"]), f"{svt[(f'{k}.vcf', 'vep', '1')]:.1f} s", f"{svt[(f'{k}.vcf', 'fastvep', '1')]:.2f} s", f"{sv_speed(f'{k}.vcf', '1'):.1f}x",
           f"{svt[(f'{k}.vcf', 'fastvep', 'all')]:.2f} s", f"{sv_speed(f'{k}.vcf', 'all'):.1f}x"]
          for name, k in (("gnomAD-SV v4.1, chromosome 21, 1-in-4 sample", R1), ("ClinVar structural variants, chromosome 21", R2))]
same("Table 12", block(PAPER, "**Table 12."), want12)
rg_ = {(r["workload"], r["threads"], r["build"]): r for r in svr}
# (label, variants, workload, threads) in the order of Table 13; a workload that was not run is not in the table
SV_WORKLOADS = [
    ("HG002, 1-in-20 sample, tab, the flags of Table 3", "202,417", "sample_1_in_20", "all"),
    ("HG002, chromosome 22, `--hgvs`", "50,284", "chr22", "1"),
    ("HG002, 1-in-5 sample, consequences only", "809,668", "genome_1_in_5", "1"),
    ("HG002, 1-in-5 sample, consequences only, tab", "809,668", "genome_1_in_5_tab", "1"),
    ("HG002, 1-in-5 sample, `--hgvs`, JSON", "809,668", "genome_1_in_5_json", "1"),
    ("HG002, 1-in-5 sample, `--hgvs`, with ClinVar", "809,668", "genome_1_in_5_clinvar", "1"),
    ("HG002, 1-in-5 sample, `--hgvs`, with ClinVar and the joint gnomAD", "809,668", "genome_1_in_5_joint", "1"),
    ("HG002, whole genome, consequences only", "4,048,342", "genome", "1"),
    ("HG002, whole genome, consequences only", "4,048,342", "genome", "all"),
    ("HG002, whole genome, `--hgvs` (the workload of Table 3)", "4,048,342", "genome_hgvs", "1"),
    ("HG002, whole genome, `--hgvs` (the workload of Table 3)", "4,048,342", "genome_hgvs", "all"),
    ("ClinVar 2-star+, 1-in-10 sample, `--acmg --pick`", "67,366", "clinvar_1_in_10_acmg", "all"),
    ("Yeast, all variants, `--hgvs`", "260,526", "yeast_full", "all"),
    ("Drosophila, all variants, `--hgvs`", "4,438,427", "drosophila_full", "all"),
    ("Arabidopsis, all variants, `--hgvs`", "12,883,854", "arabidopsis_full", "all"),
    ("Mouse, all variants, `--hgvs`", "26,062,054", "mouse_full", "all"),
]
want13 = []
for label, variants, workload, threads in SV_WORKLOADS:
    if (workload, threads, "new") not in rg_:
        continue
    base, new = float(rg_[(workload, threads, "base")]["median_s"]), float(rg_[(workload, threads, "new")]["median_s"])
    r = rg_[(workload, threads, "new")]
    # the one cell that is not the same bytes is the JSON one, by the 8 characters of `splice_donor_5th_base_variant` over `splice_region_variant`
    # in the `most_severe_consequence` of one record (benchmarks/results/identity_sweep/json_difference.txt says which)
    if r["output_identical"] != "yes":
        assert workload == "genome_1_in_5_json" and int(r["output_bytes"]) - int(r["output_bytes_base"]) == 8, f"the two builds write other bytes for {workload} than the text says"
    fastest = 100 * (float(r["min_s"]) / float(rg_[(workload, threads, "base")]["min_s"]) - 1)
    want13.append([label, variants, "10" if threads == "all" else "1", f"{base:.2f} s", f"{new:.2f} s", f"{100 * (new / base - 1):+.1f}%", f"{fastest:+.1f}%",
                   f"{float(r['instructions_paired_median_pct']):+.2f}%" if threads == "1" else "(threads wait)", f"{r['runs']}",
                   "yes" if r["output_identical"] == "yes" else f"no, {int(r['output_bytes']) - int(r['output_bytes_base'])} bytes"])
same("Table 13", block(PAPER, "**Table 13."), want13)
# the README quotes the same runs
for label, k in (("| gnomAD-SV v4.1, chromosome 21, a 1-in-4 sample |", R1), ("| ClinVar structural variants, chromosome 21 |", R2)):
    fact(README, "README", label, n_(sva[k]["records"]), f"{svt[(f'{k}.vcf', 'vep', '1')]:.1f}s", f"{svt[(f'{k}.vcf', 'fastvep', '1')]:.2f}s", f"**{sv_speed(f'{k}.vcf', '1'):.1f}×**",
         f"{svt[(f'{k}.vcf', 'fastvep', 'all')]:.2f}s", f"**{sv_speed(f'{k}.vcf', 'all'):.1f}×**")
for label, k in (("| A generated world (`crates/fastvep-cli/fixtures/vep_tab_sv`)", "world"), ("| Another world, from another seed |", "world_seed_7"), ("| The structural variants of ClinVar on chromosome 21 |", R2),
                 ("| gnomAD-SV v4.1 on chromosome 21, a 1-in-4 sample, as published |", R1), ("| The same, with the contig spelled `chr21` throughout", R1 + "_one_spelling")):
    fact(README, "README", label, n_(sva[k]["records"]), n_(sva[k]["vep_rows"]), n_(sva[k]["fastvep_rows"]), sva[k]["strict_f1"])
present_ = [x for x in SV_WORKLOADS if (x[2], x[3], "new") in rg_]
med_pct = [100 * (float(rg_[(w, th, "new")]["median_s"]) / float(rg_[(w, th, "base")]["median_s"]) - 1) for _, _, w, th in present_]
fast_pct = [100 * (float(rg_[(w, th, "new")]["min_s"]) / float(rg_[(w, th, "base")]["min_s"]) - 1) for _, _, w, th in present_]
stall_pct = max(max(abs(float(rg_[(w, th, "new")]["paired_min_pct"])), abs(float(rg_[(w, th, "new")]["paired_max_pct"]))) for _, _, w, th in present_)
ins_one = {w: float(rg_[(w, "1", "new")]["instructions_paired_median_pct"]) for _, _, w, th in present_ if th == "1"}
ins_rest = [v for w, v in ins_one.items() if w not in ("genome_1_in_5_tab", "genome_1_in_5_json")]
fastest_range = f"{min(fast_pct):+.1f}% to {max(fast_pct):+.1f}%"
medians_range = f"{min(med_pct):+.1f}% to {max(med_pct):+.1f}%"
work = (f"on one thread the two builds differ by {min(ins_rest):+.2f}% to {max(ins_rest):+.2f}% in {len(ins_rest)} of the {len(ins_one)} workloads, by {ins_one['genome_1_in_5_json']:+.2f}% for JSON "
        f"and by {ins_one['genome_1_in_5_tab']:+.2f}% for tab output, which also takes "
        f"{100 * (float(rg_[('genome_1_in_5_tab', '1', 'new')]['min_s']) / float(rg_[('genome_1_in_5_tab', '1', 'base')]['min_s']) - 1):.1f}% longer in its fastest run")
fact(README, "README", "The two builds were run in turn, the order alternating from round to round, on", f"{len(med_pct)} workloads", f"the fastest run of each build differs by {fastest_range}")
fact(README, "README", "The work that a run does is told by the instructions that it retires", work)
fact(PAPER, "3.11", "Run in turn by the build with the rules and by the build of master (Table 13)", f"{len(med_pct)} workloads")
fact(PAPER, "3.11", "A laptop varies from one round to the next by about 3%, more while someone uses it", f"by up to {stall_pct:.0f}% within a pair", f"differs by {fastest_range}", f"the medians by {medians_range}")
fact(PAPER, "3.11", "The work that a run does is told by the instructions that it retires", work)
fact(PAPER, "4.1", "and the rules for them change the time of the fastest run over small variants by", f"{fastest_range} (Table 13)")
fact(DOC_SV, "docs/STRUCTURAL_VARIANTS.md", "The fastest run of each build, over the rounds of", f"{len(med_pct)} workloads", f"differs by {fastest_range}", f"the medians by {medians_range}")
fact(DOC_SV, "docs/STRUCTURAL_VARIANTS.md", "The work that a run does is told by the instructions that it retires", work)

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

# ---- what Ensembl VEP 116.3 wrote for the runs of validation/, kept in benchmarks/results/vep_oracle/ ---------------------------------
ORACLE = HERE.parent / "benchmarks" / "results" / "vep_oracle"
oracle = lambda name: (ORACLE / name).read_text().splitlines()

# validation/run_sv_world.py: the generated worlds of other seeds, row for row and cell by cell
w_n = w_rows = w_cells = 0
for line in oracle("sv_worlds.txt"):
    m = re.match(r"seed (\d+): (\d+) records, VEP (\d+) rows, fastVEP (\d+) rows, (\d+) cells compared, (\d+) known \(divergence 11\), (\d+) unexpected, "
                 r"(\d+) rows only VEP writes, (\d+) only fastVEP writes, first row out of order: (\w+)", line)
    if not m:
        continue
    seed, _, vep_rows, fv_rows, cells, known, unexpected, only_vep, only_fv, order = m.groups()
    assert vep_rows == fv_rows and unexpected == only_vep == only_fv == "0" and order == "None" and 5 <= int(known) <= 7, f"seed {seed} of sv_worlds.txt is not VEP's"
    if seed != "185":  # the seed of the fixture is the golden of the test suite
        w_n, w_rows, w_cells = w_n + 1, w_rows + int(vep_rows), w_cells + int(cells)
assert w_n == 8, "the text says eight seeds besides the golden's own"
fact(DOC_SV, "docs/STRUCTURAL_VARIANTS.md", "Eight seeds besides the golden's own", f"{w_rows:,} rows and {w_cells:,} cells")
fact(DIV, "VEP_DIVERGENCE.md", "`validation/run_sv_world.py` makes the same comparison on worlds of other seeds", f"{w_rows:,} rows of VEP 116.3")
fact(README, "README", "`validation/run_sv_world.py` makes the first comparison on worlds of other seeds", f"{w_rows:,} rows of Ensembl VEP 116.3")

# validation/compare_pick.py: the transcript that each flagging form flags, over chromosome 22 of HG002
fp: dict[str, dict] = {}
form = None
for line in oracle("flag_pick_chr22.txt"):
    m = re.match(r"--by (\S+): (\d+) groups with a transcript row in both tables, (\d+) with another transcript flagged", line)
    if m:
        form = m.group(1)
        fp[form] = {"groups": int(m.group(2)), "differ": int(m.group(3)), "longer": 0, "ccds": 0, "skipped": 0}
        continue
    m = re.match(r"\s+(\d+)\s+(.*)", line)
    if m and form:
        n, why = int(m.group(1)), m.group(2)
        key_ = "longer" if why.startswith("the longer one") else "ccds" if why.startswith("a better term") else "skipped" if why.startswith("not a row of VEP's") else "other"
        fp[form][key_] = fp[form].get(key_, 0) + n
        continue
    m = re.match(r"\s+(Ensembl VEP|fastVEP): (\d+) variants with no transcript, (\d+) rows, (\d+) of them flagged", line)
    if m and form:
        fp[form][m.group(1)] = tuple(int(x) for x in m.groups()[1:])
assert all(f.get("other", 0) == 0 and f["longer"] + f["ccds"] + f["skipped"] == f["differ"] for f in fp.values()), "every difference of the pick is one of the three causes the text names"
for form, label in (("variant", "| `--flag_pick`, a group is a variant |"), ("allele", "| `--flag_pick_allele`, a variant and an allele |"),
                    ("allele-gene", "| `--flag_pick_allele_gene`, a variant, an allele and a gene |")):
    f = fp[form]
    fact(DIV, "VEP_DIVERGENCE.md", label, f"{f['groups']:,}", f"{f['differ']:,}", f"{f['longer']:,}", f" {f['ccds']} ", f" {f['skipped']} ")
share = [100 * fp[k]["differ"] / fp[k]["groups"] for k in fp]
fact(DIV, "VEP_DIVERGENCE.md", "They flag another transcript for", f"{min(share):.1f} to {max(share):.1f} % of the groups")
fv_variants, fv_rows_, fv_flagged = fp["variant"]["fastVEP"]
vp_variants, vp_rows_, vp_flagged = fp["variant"]["Ensembl VEP"]
assert fv_flagged == fv_rows_ == fp["allele"]["fastVEP"][2], "fastVEP flags every row of a variant that no transcript is beside"
assert fp["allele"]["Ensembl VEP"][2] == fp["allele-gene"]["Ensembl VEP"][2] == vp_rows_, "VEP flags every allele's row for the two allele forms"
fact(DIV, "VEP_DIVERGENCE.md", "A variant that no transcript is beside has only its `intergenic_variant` rows", f"{fv_variants:,} of the 50,284 variants", f"in {fv_rows_:,} rows",
     f"the 61 it leaves", f"{fv_variants:,} of the {fv_rows_:,} rows")
assert fv_rows_ - fv_variants == 61
fact(README, "README", "`--flag-pick*` flags the `intergenic_variant` rows that it left without a flag", f"{fv_variants:,} of the 50,284 variants of chromosome 22 of HG002, in {fv_rows_:,} rows")
fact(DIV, "VEP_DIVERGENCE.md", "| `--flag-pick*` left a variant that no transcript is beside without a flagged row |", f"{fv_rows_:,} of the 718,065 rows", f"the rows of {fv_variants:,} of its 50,284 variants")

# the three JSON records of HG002 whose most_severe_consequence changed: VEP writes the term of 0.6.0, and master wrote the other
js = [l.split("\t") for l in oracle("json_most_severe_consequence.txt") if l and not l.startswith("#")]
assert len(js) == 3 and all(r[3] == r[4] != r[5] for r in js), "VEP and 0.6.0 write the same term, and master another"
fact(DIV, "VEP_DIVERGENCE.md", "| The ranks of `splice_donor_5th_base_variant` and `splice_region_variant`", f"{len(js)} of the 4,048,342 JSON records", "VEP 116.3 writes `splice_donor_5th_base_variant` as the `most_severe_consequence` of all three records")

# validation/sv_buffer_size.py: Ensembl VEP with its own buffer and with a buffer of one record, over file 09 of the vep-rs set (three runs of each) and the two real sets
sections: dict[str, list[str]] = {}
cur_ = None
for line in oracle("sv_buffer_size.txt"):
    m = re.match(r"(\S+\.vcf): (\d+) records;", line)
    if m:
        cur_ = m.group(1)
        sections[cur_] = []
    elif cur_:
        sections[cur_].append(line)


def vep_runs(name: str, buffer_label: str) -> list[dict]:
    out, item = [], None
    for line in sections[name]:
        m = re.match(r"\s+VEP, (default buffer|--buffer_size 1)(?:, run \d+)?: (\d+) rows, (\d+) tuples, (\d+) in common, strict F1 ([\d.]+); only VEP (\d+), only fastVEP (\d+)", line)
        if m:
            item = {"buffer": m.group(1), "tuples": int(m.group(3)), "f1": m.group(5), "only_vep": int(m.group(6)), "only_fv": int(m.group(7)), "kinds": {}}
            out.append(item)
            continue
        m = re.match(r"\s+(single or symbolic|a mate on its own contig|a mate on another contig)\s+only VEP\s+(\d+)\s+only fastVEP\s+(\d+)", line)
        if m and item:
            item["kinds"][m.group(1)] = (int(m.group(2)), int(m.group(3)))
    return [r for r in out if r["buffer"] == buffer_label]


b09 = {k: vep_runs("09_breakends.vcf", k) for k in ("default buffer", "--buffer_size 1")}
assert len(b09["default buffer"]) == len(b09["--buffer_size 1"]) == 3
assert all(r["kinds"]["single or symbolic"] == (40, 863) and r["kinds"]["a mate on its own contig"] == (0, 0) for r in b09["default buffer"]), "the default buffer: 40 tuples only VEP's and 863 only fastVEP's on single breakends, in each run"
assert all(r["kinds"]["single or symbolic"] == (0, 0) and r["kinds"]["a mate on its own contig"] == (0, 0) for r in b09["--buffer_size 1"]), "a buffer of one: no difference on single breakends or mates on their own contig"
other_ = [r["kinds"]["a mate on another contig"] for r in b09["--buffer_size 1"]]
fact(DOC_SV, "docs/STRUCTURAL_VARIANTS.md", "Run one record at a time (`--buffer_size 1`), Ensembl VEP and fastVEP write the same rows", "863 tuples of those are only fastVEP's and 40 only VEP's",
     f"{other_[0][1]:,}, {other_[1][1]:,} and {other_[2][1]:,} tuples that only fastVEP writes", f"{other_[0][0]}, {other_[1][0]} and {other_[2][0]} that only VEP does")
fact(PAPER, "3.11", "of its 2,575 differing tuples", f"({other_[0][1]:,}, {other_[1][1]:,} and {other_[2][1]:,} tuples in three runs)")
for name in ("R1_gnomad_sv_v41_chr21_q.vcf", "R2_clinvar_chr21_sv.vcf"):
    assert len(vep_runs(name, "default buffer")) == len(vep_runs(name, "--buffer_size 1")) == 1
r1d, r1b = vep_runs("R1_gnomad_sv_v41_chr21_q.vcf", "default buffer")[0], vep_runs("R1_gnomad_sv_v41_chr21_q.vcf", "--buffer_size 1")[0]
r2d, r2b = vep_runs("R2_clinvar_chr21_sv.vcf", "default buffer")[0], vep_runs("R2_clinvar_chr21_sv.vcf", "--buffer_size 1")[0]
assert r2d["f1"] == r2b["f1"] == "0.9999", "the ClinVar set scores the same with either buffer"
fact(DOC_SV, "docs/STRUCTURAL_VARIANTS.md", "It is the buffer again, since the other contig's transcripts", f"{r1b['tuples']:,} tuples of VEP", f"{r1b['only_vep']} are not fastVEP's", f"{r1d['tuples']:,} and {r1d['only_vep']:,}",
     f"{r1b['f1']} and not {r1d['f1']}")
fact(README, "README", "its `--gff` mode pairs 905 records with transcripts of other chromosomes", f"the strict F1 is {r1b['f1']})")
fact(PAPER, "3.11", "pairs 905 of the records with transcripts of other chromosomes", f"the strict F1 is {r1b['f1']}.")

# the release binary against the build that was timed: the instructions of one thread, kept in small_variant_regression/
rel_pct, fmt_ = [], None
for line in (HERE.parent / "benchmarks" / "results" / "small_variant_regression" / "instructions_release_vs_timed.txt").read_text().splitlines():
    if line.startswith("== "):
        fmt_ = {}
    m = re.match(r"\s+(final|release)\s+n=\d+ instructions median\s+([\d.]+) G", line)
    if m and fmt_ is not None:
        fmt_[m.group(1)] = float(m.group(2))
        if len(fmt_) == 2:
            rel_pct.append(100 * (fmt_["release"] / fmt_["final"] - 1))
            fmt_ = None
assert len(rel_pct) == 3, "the release binary was counted over VCF, tab and JSON"
rel_text = f"it retires {min(rel_pct):+.1f}% to {max(rel_pct):+.1f}% of the instructions of that build"
fact(README, "README", "The binary of the release does the same work as the build that was timed", rel_text)
fact(DOC_SV, "docs/STRUCTURAL_VARIANTS.md", "The binary of the release does the same work as the build that was timed", rel_text)
fact(PAPER, "3.11", "The binary of the release does the same work as the build that was timed", rel_text)

# the release binary against the build of the structural variant table: instructions on one thread, and the timings of a session in which the machine was in use
import statistics
SVDIR = HERE.parent / "benchmarks" / "results" / "sv_head_to_head"
sv_ins: dict[tuple[str, str], list[float]] = {}
with (SVDIR / "instructions_release_vs_measured.tsv").open(newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        sv_ins.setdefault((r["input"], r["build"]), []).append(float(r["instructions"]))
sv_pct = {k: 100 * (statistics.median(sv_ins[(k, "release")]) / statistics.median(sv_ins[(k, "measured")]) - 1) for k in (R1, R2)}
with (SVDIR / "release_binary_machine_in_use" / "summary.tsv").open(newline="") as f:
    in_use = {(r["input"], r["tool"], r["threads"]): r for r in csv.DictReader(f, delimiter="\t")}
with (SVDIR / "release_binary_machine_in_use" / "environment.tsv").open() as f:
    in_use_load = next(line.split("\t")[1].split()[0] for line in f if line.startswith("load_average_at_start"))
fast_ = [float(in_use[(k + ".vcf", "fastvep", "1")]["min_s"]) for k in (R1, R2)]
med_ = [float(in_use[(k + ".vcf", "fastvep", "1")]["median_s"]) for k in (R1, R2)]
table_med = [svt[(k + ".vcf", "fastvep", "1")] for k in (R1, R2)]
assert all(abs(f / m - 1) < 0.02 for f, m in zip(fast_, table_med)), "the fastest runs of the release binary are within 2% of the medians of the table"
for lines, where, sep in ((README, "README", ""), (DOC_SV, "docs/STRUCTURAL_VARIANTS.md", ""), (PAPER, "3.11", " ")):
    fact(lines, where, "The binary that ships retires", f"{sv_pct[R1]:+.2f}% and {sv_pct[R2]:+.2f}%", f"the load average at the start was {float(in_use_load):.1f}",
         f"{fast_[0]:.2f}{sep}s and {fast_[1]:.2f}{sep}s", f"{med_[0]:.2f}{sep}s and {med_[1]:.2f}{sep}s")

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("-v", action="store_true", help="also print every check that passed")
args = ap.parse_args()
if args.v:
    print("\n".join(passed))
print(f"{len(passed)} values found in their lines, {len(failures)} failures")
for f in failures:
    print("FAIL", f)
sys.exit(1 if failures else 0)
