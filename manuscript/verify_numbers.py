#!/usr/bin/env python3
"""Check that the numbers in `manuscript/manuscript.md` (and in the performance section of `README.md`) are the numbers in `manuscript/data/`.

The text is written by hand and the CSV files are written by `build_data.py` from the runs that measured them, so this is the bridge: it
recomputes each quoted number from the CSV files, finds the sentence or table row that quotes it, and fails when that line does not
contain the recomputed value. A line is found by a phrase that occurs in no other line, so an edit that moves or rewords a claim shows up
as an anchor that no longer matches and not as a silent pass.

    manuscript/verify_numbers.py        run every check, print the failures, exit 1 if there are any
    manuscript/verify_numbers.py -v     also print every check that passed

The counts of the comparison with Ensembl VEP 116.3 are typed (`manuscript/data/vep116_provenance.md` says where they come from), so they
are looked for in the rows of `docs/VEP_DIVERGENCE.md` that they were copied from, and the text is checked against them.

Not checked: a number that is not a function of the CSV files (a date, a count read from the repository's own documents that is not one
of those above), and the wording of any sentence that passes.

Standard library only.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"


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
    results = HERE.parent / "benchmarks" / "results" / "v0.5.0-rc.1"
    with (results / relative).open(newline="") as f:
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

# ---- Table 5: formats ---------------------------------------------------------------------------------------------------
fm = table("format_comparison.csv")
want5 = [[{"vcf": "VCF (CSQ, 49 fields)", "tab": "Tab-delimited", "json": "JSON"}[r["format"]], f"{float(r['time_sec']):.2f}", f"{int(r['output_size_bytes']) / 1024:,.0f} KB"] for r in fm]
same("Table 5", block(PAPER, "**Table 5."), want5)

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
fact(PAPER, "Table 4", "| Binary/install size |", f"{float(res['binary_size']):.1f} MB")
fact(PAPER, "Table 4", "| Startup (7K transcripts, yeast) |", f"{float(res['cache_load_yeast_7k_1_thread']):.2f} sec")
fact(PAPER, "Table 4", "| Startup (510K transcripts, human) |", f"{float(res['cache_load_human_509k_1_thread']):.2f} sec", f"{float(res['cache_load_human_509k_8_threads']):.2f} sec")
fact(PAPER, "Table 4", "| Peak memory (100,000 variants, human) |", f"{float(res['peak_memory_100k_variants']):.1f} GiB")
fact(PAPER, "Table 10", "| Binary size |", f"{float(res['binary_size']):.1f} MB")

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
fact(PAPER, "5.4", "The release candidate has", f"{int(m['tests_total']):,}", f"{int(ex['passed']):,} ran and passed", "one was ignored", "Linux only")

# ---- accuracy against Ensembl VEP 116.3 -------------------------------------------------------------------------------------------
# The counts are typed (docs/VEP_DIVERGENCE.md, Part 3, is where they were measured), so each is looked for in its own row of that document.
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
hg, c5, c5h, world, edges, lone = (agree[k] for k in ("hg002_1in200", "clinvar_1in5", "clinvar_1in5_hgvs", "world", "edges", "lone_indels"))
share = lambda r: 100.0 * (int(r["cells_compared"]) - int(r["cells_differ"])) / int(r["cells_compared"])
fact(PAPER, "3.3", "On a systematic 1-in-200 sample of the GIAB HG002 call set", f"{int(hg['records']):,}", f"{int(hg['paired_rows']):,}", f"{int(hg['cells_compared']):,}", f"{int(hg['cells_differ']):,}",
     f"{100 - share(hg):.3f}%", f"`HGVS_OFFSET` {int(cols['HGVS_OFFSET']['cells_differ']):,}", f"`HGVSc` {int(cols['HGVSc']['cells_differ']):,}", f"`HGVSg` {cols['HGVSg']['cells_differ']}", f"`HGVSp` {cols['HGVSp']['cells_differ']}")
fact(PAPER, "3.3", "On a systematic 1-in-5 sample of the ClinVar 2-star+ call set", f"{int(c5['records']):,}", f"{int(c5['paired_rows']):,}", f"{int(c5['cells_differ'])} of {int(c5['cells_compared']):,}",
     f"{100 - share(c5):.5f}%", f"{int(c5h['cells_differ']):,} of {int(c5h['cells_compared']):,}", f"{100 - share(c5h):.4f}%")
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
assert ident and all(r["identical"] == "yes" for r in ident), "the release candidate's output differs from the build of 27 October: the text says it is identical"
fact(PAPER, "3.3", "The comparisons below were made on 27 October 2026", f"{len(ident)} runs in all")
big = next(r for r in ident if r["input"] == "hg002_full")
fact(PAPER, "3.3", "The comparisons below were made on 27 October 2026", f"{int(big['records']):,} records", f"{int(big['bytes_v0.5.0-rc.1']) / 1e9:.1f} GB")
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

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("-v", action="store_true", help="also print every check that passed")
args = ap.parse_args()
if args.v:
    print("\n".join(passed))
print(f"{len(passed)} values found in their lines, {len(failures)} failures")
for f in failures:
    print("FAIL", f)
sys.exit(1 if failures else 0)
