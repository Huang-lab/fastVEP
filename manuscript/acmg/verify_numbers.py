#!/usr/bin/env python3
"""Check that the numbers in `manuscript.md` are the numbers in `data/`.

The manuscript is written by hand, and every figure in it comes from a table that
`analysis/acmg_benchmark/scripts/13_manuscript_stats.py` wrote. This script is the bridge: it
recomputes each quoted number from those tables, finds the sentence that quotes it (every
sentence is on its own line), and fails if the sentence does not contain the recomputed value.
A sentence is found by a phrase that occurs in no other line, so an edit that moves or rewords
a claim shows up as an anchor that no longer matches and not as a silent pass.

Tables built by `build_tables.py` are checked there (`build_tables.py --check`), so they are
not repeated here. A number that is not a function of the data (a date, a count read from a
document, a reviewer's remark) is not checked here, and each is named in the README.

ClinVar's 2-star+ classifications are a reference and not ground truth.

Two bases. The checks read the files of a base, the classifier whose final run the tables describe: `--base v38` (the default, the
latest classifier, which is what the manuscript says) reads `headline_v38.csv` and the other `_v38` files
(`13_manuscript_stats.py`), and `--base v37` reads their twins for v37, the classifier at the close of the review rounds, with the
runs and the labels that name v38 named v37. The code is spelled for v37 and every name goes through `fname`, `rn` or `lab`, which
turn it into the base's. The text of `manuscript.md` says v38, so under `--base v37` the checks that quote a number that differs
report it: the list is what v38 changed from v37, sentence by sentence (the line, the anchor and the value v37 gives), and it is
not a failure of the data. The sentences that set v38 against v37 read v37's files under both bases. An anchor is a phrase of the
text as it stands, so it is the same under both bases. The files of the history of the versions (the trajectory, the ledger and the
reviewer's rounds) have no twin and are read as they are.

Usage
  verify_numbers.py            run every check, print the failures, exit 1 if there are any
  verify_numbers.py -v         also print every check that passed
  verify_numbers.py --base v37 check the sentences against the files of v37 (the failures are what v38 changed)
  verify_numbers.py --dump FILE write every check (line, anchor, the values the data give) to a JSON file; the dumps of the two
                                bases, set side by side, are the numbers of the text that differ between v37 and v38
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
LINES = (HERE / "manuscript.md").read_text().splitlines()
N_ALL = 673_660


def _base() -> str:
    """The base whose files are read, which the module-level loads below need before `main` parses its arguments."""
    for i, a in enumerate(sys.argv):
        if a == "--base" and i + 1 < len(sys.argv):
            return sys.argv[i + 1]
        if a.startswith("--base="):
            return a.split("=", 1)[1]
    return os.environ.get("ACMG_BASE", "v38")


BASE = _base()
if BASE not in ("v37", "v38"):
    sys.exit(f"unknown base {BASE}; choose v37 or v38")

failures: list[str] = []
passed: list[str] = []
checks_made: list[dict] = []  # every `fact`: the line, the anchor and the values the data give, for `--dump`


def fname(name: str) -> str:
    """The file of a table. Names are spelled for v37, and `_v37` is the base's suffix for another base."""
    return name.replace("_v37", f"_{BASE}")


def rn(name: str) -> str:
    """A run's name in the family of the base: `v37_cut_lp4` is `v38_cut_lp4` under v38."""
    return BASE + name[3:] if name.startswith("v37") else name


def lab(text: str) -> str:
    """A label or a column name that names the classifier, which the code spells v37."""
    return text.replace("v37", BASE)


def csv(name: str) -> pd.DataFrame:
    return pd.read_csv(DATA / fname(name))


def n(x) -> str:
    """A count with thousands separators, rounded half up."""
    return f"{int(float(x) + 0.5):,}"


def p(x, d: int = 1) -> str:
    """A proportion as a percentage with d decimals, rounded half up."""
    scale = 10 ** d
    return f"{int(float(x) * 100 * scale + 0.5) / scale:.{d}f}%"


def pp(x, d: int = 1) -> str:
    """A percentage-point difference (proportion units) as a number with d decimals."""
    return f"{abs(float(x)) * 100:.{d}f}"


WORDS = ["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen",
         "sixteen", "seventeen", "eighteen", "nineteen", "twenty"]


def w(k) -> tuple:
    """A small count as the digits or the word, capitalised or not (the text spells out some counts)."""
    k = int(k)
    return (str(k), WORDS[k], WORDS[k].capitalize()) if 0 <= k < len(WORDS) else (str(k),)


def has(value: str, line: str) -> bool:
    """A value is in the sentence as a whole number or as a phrase, and a value that is not a number (a word, a variant) is found as written.
    A number that only occurs inside something else is a check that passes for no reason, and each of these did: the 1 of `1,795` for a count
    that had become 1, the 2 of `PM2`, the 38 of `v38`, the 5 of `4e-5`, the 7 of the hash `70b7c684` and the 8 of `8%`. So a number is found
    only where it is a number of its own: `24` is not found in `124`, `2.4`, `24,5`, `PM24` or `4e-24`, and `8` is not found in `8%` (a value
    that is a percentage is spelled with its %). A value that is a number with words after it (`6.1 points`, `8 at full strength`) has to
    begin where a number begins in the sentence."""
    if re.fullmatch(r"-?[\d,.]*\d%?", value):
        tail = r"(?![\w%]|,\d|\.\d)" if re.fullmatch(r"-?[\d,]+", value) else r"(?!\w|,\d|\.\d)"
        return re.search(r"(?<![\w.,])(?<!e-)" + re.escape(value) + tail, line) is not None
    if re.match(r"-?\d", value):
        return any(m.start() == 0 or not re.match(r"[\w.,]", line[m.start() - 1]) for m in re.finditer(re.escape(value), line))
    return value in line


def fact(anchor: str, *values) -> None:
    """The one line containing `anchor` must contain each value (a tuple is a set of acceptable spellings)."""
    hits = [i for i, line in enumerate(LINES) if anchor in line]
    if len(hits) != 1:
        failures.append(f"{'no line' if not hits else f'{len(hits)} lines'} contain the anchor {anchor!r}")
        checks_made.append({"anchor": anchor, "line": None, "values": [(v if isinstance(v, tuple) else (str(v),))[0] for v in values]})
        return
    line = LINES[hits[0]]
    checks_made.append({"anchor": anchor, "line": hits[0] + 1, "values": [(v if isinstance(v, tuple) else (str(v),))[0] for v in values]})
    for v in values:
        alts = v if isinstance(v, tuple) else (str(v),)
        if any(has(a, line) for a in alts):
            passed.append(f"{anchor[:50]!r}: {alts[0]}")
        else:
            failures.append(f"line {hits[0] + 1}: expected {alts[0]!r} in the sentence {anchor!r}")


# ---------------------------------------------------------------------------------- the data

runs = csv("runs_v37_variants.csv").set_index("name")
base = runs.loc[BASE]
N_OPP = int(base.opposite)  # the opposite-direction calls of the base: 51 for v38, 52 for v37
N_REV = int(csv("round6_rows_v37.csv").still_opposite.astype(str).eq("True").sum())  # the ones of them that the reviewer saw in the last round: 48
_opp_rows = csv("opposite_direction_v37.csv")
_is_fb = _opp_rows.clinvar.isin(("Pathogenic", "Likely_pathogenic"))
REV_FB, REV_FP = int((_is_fb & (_opp_rows.ruling != "not reviewed")).sum()), int((~_is_fb & (_opp_rows.ruling != "not reviewed")).sum())  # the ones that the reviewer saw, by direction: 24 and 24
# v37's own rows and steps, for the checks of Section 3.8 that set v38 against them (the same objects under v37)
runs37 = pd.read_csv(DATA / "runs_v37_variants.csv").set_index("name")
base37 = runs37.loc["v37"]
rv38 = pd.read_csv(DATA / "runs_v38_variants.csv").set_index("name")  # v38's family, whatever the base: the sentences that set v38 against v37 read it, so they are the same under both bases
b38 = rv38.loc["v38"]
facts38 = pd.read_csv(DATA / "facts_v38.csv").set_index("name")["value"]
steps37 = pd.read_csv(DATA / "threshold_steps_v37.csv")
traj = csv("trajectory.csv").set_index("order")
steps = csv("threshold_steps_v37.csv")
head = csv("headline_v37.csv").set_index("metric", drop=False)
c5 = csv("confusion_5x5_v37.csv").pivot(index="clinvar", columns="fastvep", values="n")
c3 = csv("confusion_3x3_v37.csv")
facts = csv("facts_v37.csv").set_index("name")["value"]
fv = lambda k: float(facts[k]) if k not in ("pvs1_benign_common_genes", "clingen_weak_only_genes_with_semicolon_label") else facts[k]
# How the classifier read inheritance up to v37 (the words of the disease names): the sentences that say what v37 read are about v37
# whichever base is checked, so they read v37's facts; under v37 these are the same numbers as `fv` gives.
facts37 = pd.read_csv(DATA / "facts_v37.csv").set_index("name")["value"]
fv37 = lambda k: float(facts37[k])
stars = csv("by_review_status_v37.csv")
cons = csv("by_consequence_v37.csv")
pts = csv("points_by_clinvar_v37.csv")
rep = csv("combining_rules_replay_v37.csv").set_index("scheme")
val = csv("replay_validation_v37.csv").set_index("run")
pm2dep = csv("pm2_dependence_v37.csv").set_index("rule")
vcep = csv("vcep_bars_v37.csv").set_index("bar")
ledger = csv("change_ledger.csv")
rounds = csv("review_rounds.csv").set_index("round")
plp_pts = csv("plp_to_vus_points_v37.csv").set_index("pathogenic_points")
plp_cons = csv("plp_to_vus_by_consequence_v37.csv")
conc = csv("plp_to_vus_gene_concentration_v37.csv").set_index("top_genes")
pvs1 = csv("pvs1_gap_v37.csv").set_index("measure")["n"]
sigs = csv("signatures_offdiagonal_v37.csv")
gene = csv("gene_level_v37.csv")
res_by_group = csv("opposite_direction_by_group_v37.csv").set_index("group")
r6 = csv("round6_rows_v37.csv")
temporal = csv("temporal_v37.csv")
tclass = csv("temporal_by_class_v37.csv")
trate = csv("temporal_opposite_rates_v37.csv").set_index("set")
tdir = csv("temporal_direction_changes_v37.csv")
aprl = csv("april_opposite_in_october_v37.csv")
rules64 = csv("round6_rules_v37.csv")
baseline = csv("baselines_v37.csv").set_index("rule")
curve = csv("baselines_curve_v37.csv")
mix = csv("temporal_agreement_mix_v37.csv").set_index("set")
hi_af = csv("temporal_high_af_benign_v37.csv").set_index("set")
prop = csv("temporal_proposals_v37.csv").set_index("change")
vaud = csv("vcep_audit_v37.csv").set_index("status")["entries"]
vaud_pm2 = csv("vcep_audit_pm2_v37.csv").set_index("status")["entries"]
new_idx = csv("temporal_new_in_index_v37.csv").set_index("measure")["n"]
blockers = csv("temporal_high_af_blockers_v37.csv")
r2v = csv("baselines_r2_vs_v37.csv")["n"].tolist()
bp1_idx = csv("bp1_index_counts_v37.csv").set_index("gene")
hold = pd.read_csv(HERE.parents[1] / "analysis/acmg_benchmark/results/v37_manuscript/holdout_runs.tsv", sep="\t")
hold37 = hold
if BASE != "v37":  # the seven runs of the base's family on the October release, which share a file with Section 3.8's four
    _h = pd.read_csv(HERE.parents[1] / f"analysis/acmg_benchmark/results/v37_manuscript/holdout_runs_{BASE}.tsv", sep="\t")
    hold = _h[_h.name.isin({rn(x) for x in hold37.name})]
tiers = csv("precision_by_tier_v37.csv").set_index("call")
ratio = csv("temporal_ratio_v37.csv")
stdr = csv("temporal_standardised_recall_v37.csv").set_index("clinvar")
bands = csv("temporal_recall_by_af_band_v37.csv")
gap_reasons = csv("pvs1_gap_reasons_v37.csv")
gap_const = csv("pvs1_gap_constrained_v37.csv").set_index("gene")
closure = csv("pvs1_gap_closure_v37.csv").set_index("run")
closure_fp = csv("pvs1_gap_closure_new_false_pathogenic_v37.csv")
closure_genes = csv("pvs1_gap_closure_genes_v37.csv")
missing = csv("clingen_missing_summary_v37.csv").set_index("measure")["n"]
# the summary keeps three kinds of row apart: curations the April export lacked, rows that only re-label one it lists, and the row dated on the day of the export
LACKED = "curations dated before 30 April that the export lacked"
ROWS = "rows of the 11 October export dated on or before 30 April whose gene and disease pair the April export lacks"
RELABEL = "  of them re-labelled (the April export lists the same assertion under another disease identifier)"
EXPORT_DAY_ROW = "  of them dated 30 April itself, after the April file was downloaded"
DAY_GENE_N = next(k for k in missing.index if k.startswith("variants of "))
DAY_GENE_CALLS = "  whose call differs in the run that reads the completed file"
DAY_GENE_CHANGED = "  whose call or criteria differ in the run that reads the completed file"
DAY_GENE_TAGS = ("  with the established ClinGen tag AR alone in the main run", "  with the established ClinGen tags AD and AR in the run that reads the completed file")
DAY_GENE_LOST = ("  that lose a PM2 criterion in the run that reads the completed file", "  lowest allele frequency in gnomAD among those (0 when there are none)",
                 "  highest allele frequency in gnomAD among those (0 when there are none)")
prov = csv("provenance_v37.csv")
controls = csv("controls_v37.csv").set_index("run")
veto = csv("replay_veto_v37.csv").iloc[0]
hooknew = csv("hook_new_opposite_v37.csv")
crit = csv("criteria_firing_v37.csv")
status = csv("criteria_status.csv").set_index("criterion")  # the 28 criteria; made by 18_criteria_status.py from a sample of the benchmark annotated to JSON
cftr_null = csv("cftr_null_variants_v37.csv")
EARLIER = {rn("v37_bs2off"), rn("v37_bs2prev_3e-5"), rn("v37_bs2prev_1e-4"), rn("v37_bs2prev_1e-3"), rn("v37_bs2prev_3e-3"), rn("v37_bs2ac_20"), rn("v37_bs2ac_100"), rn("v37_bs2ac_1000")}
harness = runs[~runs.index.isin(EARLIER | {BASE})]  # the 89 runs of 11_manuscript_runs.py on the April set
mainline = traj[traj.mainline.astype(str) == "True"]
committed_base = int(sum(c5.loc[t, c] for t in c5.index for c in ("Pathogenic", "Likely_pathogenic", "Likely_benign", "Benign")))  # the calls of one of those four classes
assert committed_base == int(traj.loc[int(BASE[1:])].committed), f"the committed calls of {BASE} are the trajectory's"


def new_york_date(recorded_at: str) -> str:
    """The date of a time that a provenance file records in UTC, in New York, where the dates of the paper are."""
    import datetime as dt
    from zoneinfo import ZoneInfo

    if not recorded_at:
        return ""
    t = dt.datetime.fromisoformat(recorded_at.replace("Z", "+00:00"))
    return t.astimezone(ZoneInfo("America/New_York")).date().isoformat()


def rerun_recorded_at(name: str) -> str:
    """When the repetition of a run by release 0.5.0 was recorded (UTC)."""
    f = HERE.parents[1] / "analysis/acmg_benchmark/results/v38_v0.5.0/provenance" / f"{name}.txt"
    return dict(line.split(": ", 1) for line in f.read_text().splitlines() if ": " in line).get("recorded_at", "") if f.exists() else ""


def same_96() -> bool:
    """The family of the base and the family of v37 are the same 96 runs: the names pair up but for the prefix, 89 on the April set and 7 on the October release."""
    prov37_ = pd.read_csv(DATA / "provenance_v37.csv")
    return ({"v37" + x[3:] for x in prov.run} == set(prov37_.run) and {"v37" + x[3:] for x in hold.name} == set(hold37.name)
            and len(prov37_) == 89 and len(hold37) == 7)


def tier_total() -> int:
    return int(c5.loc["Pathogenic", "Likely_pathogenic"] + c5.loc["Likely_pathogenic", "Pathogenic"] + c5.loc["Likely_benign", "Benign"] + c5.loc["Benign", "Likely_benign"])


def cell(a: str, b: str) -> int:
    x = c3[(c3.clinvar == a) & (c3.fastvep == b)].iloc[0]
    return int(x.n)


def run(name: str):
    return runs.loc[rn(name)]


def star_cell(star: int, a: str, b: str):
    """One cell of the concordance by review status."""
    return stars[(stars.stars == star) & (stars.clinvar == a) & (stars.fastvep == b)].iloc[0]


def sci1(x: float) -> str:
    """A small probability to one digit, as the text writes it: 3.4e-5 is 3e-5."""
    return re.sub(r"e-0*", "e-", f"{float(x):.0e}")


def hl(metric: str):
    """The first row of a metric (for a rate, the row that has an interval; a count that shares the name comes later)."""
    x = head.loc[metric]
    return x.iloc[0] if isinstance(x, pd.DataFrame) else x


def hlc(metric: str):
    """The last row of a metric: the count that shares its name with a rate."""
    x = head.loc[metric]
    return x.iloc[-1] if isinstance(x, pd.DataFrame) else x


# --------------------------------------------------------------------------------- Section 2


def databases() -> None:
    """The first databases and the complete ones (Sections 2.1 and 3.9): what the first lacked, where, and what the runs say under each.

    `21_compare_databases.py` writes the five `databases_*.csv` files from the per-variant tables of the runs on both stacks, and
    they hold both versions whatever the base is.
    """
    pres = pd.read_csv(DATA / "databases_presence.csv")
    cause = pd.read_csv(DATA / "databases_causes.csv").set_index("source")
    blocks = pd.read_csv(DATA / "databases_lost_blocks.csv")
    manifest = pd.read_csv(HERE.parents[1] / "analysis/acmg_benchmark/results/complete_databases/manifest.tsv", sep="\t")
    ap = pres[pres.set == "April set"].set_index("source")
    sa, pp, gn, rv = ap.loc["SpliceAI"], ap.loc["PhyloP"], ap.loc["gnomAD"], ap.loc["REVEL"]
    fact("are read from databases built from the whole of that file, one for each source and chromosome", str(len(manifest)))
    fact("Of the 673,278 variants of the benchmark that have an alternative allele", n(sa["variants"]), n(sa["first"]), n(pp["first"]), n(sa["complete"]), n(pp["complete"]),
         p(sa.only_complete / sa.complete), p(pp.only_complete / pp.complete),
         "no variant had a value only in the first databases" if sa.only_first == 0 and pp.only_first == 0 else "?",
         "where both had one the values were the same" if sa.differ_where_both == 0 and pp.differ_where_both == 0 else "?")
    fact("The first databases lack their records in 14 stretches of SpliceAI and 14 of PhyloP", n(cause.loc["SpliceAI", "batches_with_a_tail"]), n(cause.loc["PhyloP", "batches_with_a_tail"]), "200",
         n(cause.loc["PhyloP", "indels_with_the_second_base_outside_the_regions"]),
         "reaches the end of the batch" if cause.loc["SpliceAI", "unexplained"] == 0 and cause.loc["PhyloP", "unexplained"] == 0 else "?",
         "more indels" if cause.loc["PhyloP", "indels_with_the_second_base_outside_the_regions"] == cause.loc["PhyloP", "variants_lost"] - cause.loc["PhyloP", "in_the_tail_of_a_batch"] else "?")
    fact("On the April set the regional gnomAD databases of the benchmark and the complete ones gave the same call", n(N_ALL),
         "so the counts that the frequency criteria read were not affected" if gn.only_complete == 0 and gn.only_first == 0 and gn.differ_where_both == 0 else "?",
         "REVEL" if rv.only_complete == 0 and rv.only_first == 0 and rv.differ_where_both == 0 else "?")
    # the batches that the section on the databases names are the rows of the file: the same batches for both sources, with the same first lost position
    both = blocks[blocks.spliceai_lost.notna() & blocks.phylop_lost.notna()]
    if len(both) != len(blocks) or len(blocks) != int(cause.loc["SpliceAI", "batches_with_a_tail"]):
        failures.append("databases_lost_blocks.csv: the batches that lose SpliceAI records are not the batches that lose PhyloP records, or not the number that databases_causes.csv counts")

    # Section 3.9: the runs on the first databases and on the complete ones
    ovn = pd.read_csv(DATA / "databases_old_vs_new.csv")
    chg = pd.read_csv(DATA / "databases_opposite_changes.csv", dtype=str, keep_default_na=False)
    a_set, c_set, nw_set, w_set = "April set", "October release, the variants of the first design", "October release, new variants of the first design", "October release, all variants"

    def g(version: str, scope: str, measure: str, db: str) -> float:
        return float(ovn[(ovn.version == version) & (ovn.set == scope) & (ovn.measure == measure)].iloc[0][db])

    def rate(version: str, scope: str, measure: str, db: str) -> str:
        return p(g(version, scope, measure, db) / g(version, scope, "variants", db), 2)

    def both_dbs(version: str, scope: str, measure: str, fmt=n) -> list[str]:
        return [fmt(g(version, scope, measure, "first")), fmt(g(version, scope, measure, "complete"))]

    d_dir = lambda v: g(v, a_set, "same_direction", "complete") - g(v, a_set, "same_direction", "first")  # noqa: E731
    fact("Under v38 agreement in direction rises from 524,380", *both_dbs("v38", a_set, "same_direction"), rate("v38", a_set, "same_direction", "first"), rate("v38", a_set, "same_direction", "complete"),
         *both_dbs("v38", a_set, "exact"), rate("v38", a_set, "exact", "first"), rate("v38", a_set, "exact", "complete"), n(g("v38", a_set, "variants", "first")),
         *both_dbs("v38", a_set, "opposite"), *both_dbs("v38", a_set, "false_benign"),
         "stay at " + n(g("v38", a_set, "false_pathogenic", "complete")) if g("v38", a_set, "false_pathogenic", "first") == g("v38", a_set, "false_pathogenic", "complete") else "?")
    fact("Under v37 agreement in direction rises from 522,886", *both_dbs("v37", a_set, "same_direction"), rate("v37", a_set, "same_direction", "first"), rate("v37", a_set, "same_direction", "complete"),
         *both_dbs("v37", a_set, "opposite"),
         n(g("v38", a_set, "same_direction", "first") - g("v37", a_set, "same_direction", "first")), n(g("v38", a_set, "same_direction", "complete") - g("v37", a_set, "same_direction", "complete")))
    fact("The complete databases add 1.5 points of agreement in direction and three contradictions", f"{d_dir('v38') / g('v38', a_set, 'variants', 'first') * 100:.1f}",
         w(int(g("v38", a_set, "opposite", "complete") - g("v38", a_set, "opposite", "first")))[1],
         "as much to v37" if f"{d_dir('v37') / g('v37', a_set, 'variants', 'first') * 100:.1f}" == f"{d_dir('v38') / g('v38', a_set, 'variants', 'first') * 100:.1f}" or abs(d_dir("v37") - d_dir("v38")) < 300 else "?")
    fact("Recall of ClinVar's benign-side variants rises from 62.9% to 66.4%", p(g("v38", a_set, "blb_recall", "first")), p(g("v38", a_set, "blb_recall", "complete")),
         n(g("v38", a_set, "same_benign", "complete") - g("v38", a_set, "same_benign", "first")), *both_dbs("v38", a_set, "same_benign"),
         p(g("v38", a_set, "plp_recall", "complete")) if g("v38", a_set, "plp_recall", "first") == g("v38", a_set, "plp_recall", "complete") else "?",
         n(g("v38", a_set, "same_pathogenic", "complete")) if g("v38", a_set, "same_pathogenic", "first") == g("v38", a_set, "same_pathogenic", "complete") else "?")
    fact("The records that the first databases lacked are SpliceAI and PhyloP values that BP4, BP7 and PP3 read", *[n(g("v38", a_set, m, "complete") - g("v38", a_set, m, "first")) for m in ("BP4", "BP7", "PP3")],
         *both_dbs("v38", a_set, "BP4"), *both_dbs("v38", a_set, "BP7"), *both_dbs("v38", a_set, "PP3"))
    fact("PVS1 on canonical splice variants, which has to pass a check on SpliceAI", w(int(g("v38", a_set, "PVS1_on_canonical_splice", "complete") - g("v38", a_set, "PVS1_on_canonical_splice", "first")))[1],
         n(g("v38", a_set, "canonical_splice_variants", "complete")), *both_dbs("v38", a_set, "PVS1_on_canonical_splice"))
    # the three variants that enter the opposite-direction calls
    entered = chg[(chg.version == "v38") & (chg.change == "enters the opposite-direction calls")]
    entered37 = chg[(chg.version == "v37") & (chg.change == "enters the opposite-direction calls")]
    if len(chg[(chg.version == "v38") & (chg.change != "enters the opposite-direction calls")]):
        failures.append("databases_opposite_changes.csv: a v38 variant leaves the opposite-direction calls or changes its call inside them, which Section 3.9 does not say")
    pos = lambda gene: entered[entered.gene == gene].iloc[0]  # noqa: E731
    tm, vw, gt = pos("TMEM216"), pos("VWF"), pos("GTPBP3")
    r6 = csv("round6_rows_v37.csv")
    seen = set(zip(r6.chrom.astype(str), r6.pos.astype(str), r6.ref, r6.alt))
    fact("Each is called Likely benign on a SpliceAI value of 0.07 or below with BP4", *[f"{x.gene} {x.chrom}:{x.pos} {x.ref}>{x.alt}" for x in (tm, vw, gt)],
         max(float(x.complete_spliceai_max) for x in (tm, vw, gt)) <= 0.07 and f"{float(tm.complete_spliceai_max):.2f}" or "?",
         f"{float(tm.complete_spliceai_max):.2f}", f"{float(vw.complete_spliceai_max):.2f}", f"{float(gt.complete_spliceai_max):.2f}",
         f"{float(tm.complete_phylop):g}", f"{float(vw.complete_phylop):g}", f"{float(gt.complete_phylop):.3f}",
         "BS1 and BP4" if tm.complete_criteria == "BS1&BP4" else "?", "BP4 and BP7" if vw.complete_criteria == "BP4&BP7" else "?", "PM2, BP4 and BP7" if gt.complete_criteria == "PM2_Supporting&BP4&BP7" else "?",
         "and two with BP7 as well" if sum("BP7" in x.complete_criteria for x in (tm, vw, gt)) == 2 and all("BP4" in x.complete_criteria for x in (tm, vw, gt)) else "?")
    fact("The three contradictions that the complete databases add are ClinVar-Pathogenic variants that the first databases left as VUS", w(len(entered))[1],
         "ClinVar-Pathogenic variants that the first databases left as VUS" if all(x.clinvar == "Pathogenic" for x in (tm, vw, gt)) and all(x.first_call == "VUS" and x.complete_call == "Likely_benign" for x in (tm, vw, gt)) else "?")
    fact("None is a Benign call, none is among the 64 calls that the reviewer saw in the last round", "None is a Benign call" if all(x.complete_call != "Benign" for x in (tm, vw, gt)) else "?",
         "none is among the 64 calls" if not any((x.chrom, x.pos, x.ref, x.alt) in seen for x in (tm, vw, gt)) else "?",
         "the same three enter under v37" if sorted(zip(entered37.gene, entered37.pos)) == sorted(zip(entered.gene, entered.pos)) else "?")
    fact("On the October release the complete databases add 1.7 points of agreement in direction", f"{(g('v38', c_set, 'same_direction', 'complete') - g('v38', c_set, 'same_direction', 'first')) / g('v38', c_set, 'variants', 'first') * 100:.1f}")
    fact("The first drafts scored the 772,585 variants of ClinVar's release of 4 October 2026", n(g("v38", c_set, "variants", "first")), n(g("v38", w_set, "variants", "complete")))
    fact("On the 772,585, agreement in direction rises from 593,343", *both_dbs("v38", c_set, "same_direction"), rate("v38", c_set, "same_direction", "first"), rate("v38", c_set, "same_direction", "complete"),
         *both_dbs("v38", c_set, "opposite"), *both_dbs("v38", c_set, "false_benign"),
         "stay at " + n(g("v38", c_set, "false_pathogenic", "complete")) if g("v38", c_set, "false_pathogenic", "first") == g("v38", c_set, "false_pathogenic", "complete") else "?",
         n(g("v38", nw_set, "variants", "first")), *both_dbs("v38", nw_set, "same_direction"), rate("v38", nw_set, "same_direction", "first"), rate("v38", nw_set, "same_direction", "complete"),
         *both_dbs("v38", nw_set, "opposite"))
    out_n = g("v38", w_set, "variants", "complete") - g("v38", c_set, "variants", "complete")
    out_dir = g("v38", w_set, "same_direction", "complete") - g("v38", c_set, "same_direction", "complete")
    fact("The 8,550 variants that the first drafts left out add one opposite-direction call", n(out_n), w(int(g("v38", w_set, "opposite", "complete") - g("v38", c_set, "opposite", "complete")))[1],
         p(out_dir / out_n), n(out_dir), n(g("v38", w_set, "same_direction", "complete")), rate("v38", w_set, "same_direction", "complete"), n(g("v38", w_set, "opposite", "complete")),
         n(g("v38", w_set, "false_benign", "complete")), n(g("v38", w_set, "false_pathogenic", "complete")))
    sa_share, pp_share = sa["only_complete"] / sa["complete"], pp["only_complete"] / pp["complete"]
    gain = d_dir("v38") / g("v38", a_set, "variants", "first") * 100
    fact("The first SpliceAI and PhyloP databases of the benchmark were a second instance", p(sa_share), p(pp_share))
    fact("Check each input against an independent listing, because a missing record degrades quietly", p(sa_share, 0), p(pp_share, 0) if p(sa_share, 0) == p(pp_share, 0) else "?")
    fact("The first drafts, including the version that the co-authors read, rested on SpliceAI and PhyloP databases that lacked", p(sa_share), p(pp_share))
    traj_ = pd.read_csv(DATA / "trajectory.csv").set_index("order")
    first_rate = {k: float(traj_.loc[k].same_direction_rate) for k in range(13, 37) if k in traj_.index and traj_.loc[k].mainline}
    first_rate[37] = g("v37", a_set, "same_direction", "first") / g("v37", a_set, "variants", "first")
    first_rate[38] = g("v38", a_set, "same_direction", "first") / g("v38", a_set, "variants", "first")
    ks = sorted(first_rate)
    steps_first = {k: abs(first_rate[k] - first_rate[j]) for j, k in zip(ks, ks[1:])}
    biggest = max(steps_first, key=steps_first.get)
    fact("The completion of the databases moved agreement in direction by 1.5 points", f"{d_dir('v38') / g('v38', a_set, 'variants', 'first') * 100:.1f}",
         f"{steps_first[biggest] * 100:.2f} points" if biggest == 38 and d_dir("v38") / g("v38", a_set, "variants", "first") > steps_first[biggest] else "?", "v38's")

    # the sweeps and the effects of the family on both stacks (`21_compare_databases.py` makes them with `20_compare_families.py`'s definitions of a shape and of a change)
    shapes = pd.read_csv(DATA / "databases_sweep_shapes.csv", dtype={"settings_that_change_order": str}).set_index("sweep")
    swp = pd.read_csv(DATA / "databases_sweeps.csv", dtype={"setting": str})
    fam = pd.read_csv(DATA / "databases_family_summary.csv")
    cnt = {r["name"]: r["value"] for _, r in fam[fam.text.isna()].iterrows()}
    kept = [k for k in shapes.index if not shapes.loc[k, "differs"]]
    sd = lambda label, value, db: float(swp[(swp.sweep == label) & (swp.setting == value)].iloc[0][f"same_direction_{db}"])  # noqa: E731
    gaps = []
    for label in shapes.index:
        for pair in str(shapes.loc[label, "settings_that_change_order"]).split("; "):
            if " and " in pair:
                a_, b_ = pair.split(" and ")
                gaps += [abs(sd(label, a_, db) - sd(label, b_, db)) for db in ("first", "complete")]
    fb = lambda value, db: float(swp[(swp.sweep == "BP7 PhyloP bar") & (swp.setting == value)].iloc[0][f"false_benign_{db}"])  # noqa: E731
    floor = shapes.loc["BS2 carrier floor, bar off"]
    nsp = fam[fam["name"] == "effect_that_changed"]
    big_effects = nsp[nsp.value == 1]
    spl = re.findall(r"-?\d+", big_effects.iloc[0].text.split(":", 1)[1]) if len(big_effects) == 1 and "nospliceai opposite" in big_effects.iloc[0].text else []
    fact("Of the 12 sweeps, 4 differ in a detail between the two sets of databases and 8 do not", str(int(cnt["sweeps"])), str(len(shapes[shapes.differs])), str(len(kept)),
         "1e-3" if shapes.loc["BS1 bar", "knee_first"] == shapes.loc["BS1 bar", "knee_complete"] == 0.001 else "?",
         "300 nucleotides" if shapes.loc["BP7 far boundary", "knee_first"] == shapes.loc["BP7 far boundary", "knee_complete"] == 300 else "?",
         "under both" if int(cnt["sweeps_that_differ_in_the_knee"]) == 1 and shapes.loc["BS2 carrier floor, bar off", "differs"] else "?")
    fact("The details are the order of two settings of the BS1 bar", n(max(gaps)),
         "one reference-discordant benign call on the complete databases and none on the first" if fb("3", "complete") - fb("2.0", "complete") == 1 and fb("3", "first") - fb("2.0", "first") == 0 else "?",
         f"{floor.knee_factor_first:.1f}", f"{floor.knee_factor_complete:.1f}", f"{int(floor.knee_complete)} alleles" if pd.isna(floor.knee_first) else "?",
         "BS2 carrier bar" if shapes.loc["BS2 carrier bar, dominant-only genes", "pairs_that_change_order"] > 0 and shapes.loc["BS1 bar", "pairs_that_change_order"] > 0 else "?")
    fact("7 changed sign, appeared or vanished, or changed by more than 20%, and 6 of the 7", str(int(cnt["effects"])), w(int(cnt["effect_measures"]))[1], str(int(cnt["effects_that_changed"])),
         str(int(cnt["effects_that_changed_and_are_small"])), "20%", "6 of the 7" if f"{int(cnt['effects_that_changed_and_are_small'])} of the {int(cnt['effects_that_changed'])}" == "6 of the 7" else "?")
    fact("The one that is larger is the removal of SpliceAI", str(abs(int(spl[1]))) if len(spl) == 2 else "?", str(abs(int(spl[0]))) if len(spl) == 2 else "?")


def methods() -> None:
    # Class sizes and stars of the reference
    cls = {t: int(c5.loc[t].sum()) for t in ("Pathogenic", "Likely_pathogenic", "VUS", "Likely_benign", "Benign")}
    fact("That gives 673,660 records:", n(cls["Pathogenic"]), n(cls["Likely_pathogenic"]), n(cls["VUS"]), n(cls["Likely_benign"]), n(cls["Benign"]))
    s = stars.groupby("stars").row_total.sum()
    tot_by_star = {k: int(stars[stars.stars == k].drop_duplicates(["clinvar"]).row_total.sum()) for k in (2, 3, 4)}
    fact("That gives 673,660 records:", n(tot_by_star[2] + 0), n(tot_by_star[3]), n(tot_by_star[4]))
    fact("The VUS class is the", n(cls["VUS"]), "two RYR1 records" if fv("labels_conflicting_mapped_to_vus") == 2 else "?",
         "fastVEP calls both Likely pathogenic" if fv("labels_conflicting_called_pathogenic_side") == 2 else "?")
    fact("Four records of CFH and CFI", "Four" if fv("labels_lowpen_pair_mapped_to_pathogenic") == 4 else "?")
    fact("We kept both labels so that the 38 runs stay comparable", str(int(traj.index.max())), "by one variant" if fv("labels_lowpen_pair_exact_gain_if_relabelled") == 1 else "?",
         "by two" if fv("labels_conflicting_called_pathogenic_side") == 2 else "?")
    fact("Among the 673,660 records, 382 have `.` as their alternative allele", n(N_ALL), n(hl("variants_no_call").value), n(fv("alt_dot_records")))
    fact("Each of them names the reference sequence itself in its HGVS expression",
         "every one of them as benign-side" if fv("alt_dot_hgvs_is_reference") == fv("alt_dot_clinvar_benign_side") == fv("alt_dot_records") else "?",
         n(c5.loc["Benign", "NoCall"]), n(c5.loc["Likely_benign", "NoCall"]))
    fact("There is no alternative allele to annotate, so they receive no call", n(hl("variants_no_call").value), n(c5.loc[["Likely_benign", "Benign"]].sum().sum()),
         "no call" if fv("alt_dot_records_called") == 0 else "?")
    pr_bl = cls["Pathogenic"] + cls["Likely_pathogenic"] + cls["Likely_benign"] + cls["Benign"]
    fact("Rates use all 673,660 records as the denominator", f"{hl('opposite_rate').value * 100:.4f}%", f"{hl('opposite').value / pr_bl * 100:.3f}%", n(pr_bl),
         f"{hl('opposite').value / committed_base * 100:.3f}%", n(committed_base))
    fact("the reviewer supported ClinVar's class for 13 of the 24 reference-discordant benign calls", str(int(fv("reviewer_supported_clinvar_false_benign"))), str(REV_FB),
         str(int(fv("reviewer_supported_clinvar_false_pathogenic"))), str(REV_FP))
    fact("scored matrices are stored for 32 of the 38 runs", str(int((traj.source == "scored matrix on disk").sum())), str(int(traj.index.max())), "v7 to v38" if traj[traj.source == "scored matrix on disk"].index.min() == 7 else "?")
    # 2.1: the criteria
    # 2.1: the criteria. The status table is the engine's own report on every 40th record of the benchmark, and the script that makes it
    # stops if the engine disagrees with the groups below, so here the sentences are held to the table.
    group = status.group
    codes_of = lambda g: list(status.index[group == g])
    from_variant, genotypes, functional, case, off = (codes_of(g) for g in ("from the variant", "needs genotypes of a proband", "needs functional assays", "needs case-level data", "off by default"))
    n_not_applied = len(status) - len(from_variant)
    fact("The engine implements all 28 ACMG/AMP criteria", str(len(status)), str(len(status) - len(off)))
    fact("From it 15 of the 28 criteria can be applied", str(len(from_variant)), str(len(status)), *from_variant)
    assert sorted(from_variant) == sorted(crit.criterion), "the 15 criteria of the benchmark are the 15 of criteria_firing_v37.csv"
    fact("The other 13 are not applied in this paper, and they fall into four groups", str(n_not_applied), w(len({group[c] for c in status.index if group[c] != "from the variant"}))[1])
    fact("Four need the genotypes of a proband and of its parents", w(len(genotypes)), *genotypes)
    fact("Two need functional assays, PS3 and BS3", w(len(functional)), *functional)
    fact("Five need case-level data that the engine does not read", w(len(case)), *case)
    fact("PP5 and BP6 are off by default.", *off)
    sample = int(status.records.max())
    # the step of the sample is not in `data/`, but the sample is, and 673,660 over it is 40 to the nearest whole number
    fact("we annotated every 40th record of the benchmark", n(sample), f"{round(N_ALL / sample)}th")
    never = [c for c in status.index if status.records[c] and not status.evaluated[c]]
    fact("were not evaluated for any of them", *never)
    assert sorted(never) == sorted(genotypes + functional + case) or sorted(never) == sorted([c for c in genotypes if c != "PM3"] + functional + case), "the criteria never evaluated"
    fact("PM3 was evaluated for 9,918", n(status.evaluated["PM3"]), "met for none" if status.met["PM3"] == 0 else "?")
    ok = all(status.evaluated[c] > 0 and status.met[c] > 0 for c in from_variant)
    fact("Each of the other 15 was evaluated for some records and met for some", "15" if ok and len(from_variant) == 15 else "?")
    prov37 = pd.read_csv(DATA / "provenance_v37.csv")
    # The harness runs were made at once, several at a time on one machine, and are not timings; the timing is release 0.5.0 timed alone on 7 November 2026
    # (`results/v38_v0.5.0/annotation_time.tsv`), with the first databases (the file is the one that the first drafts quoted and was not measured again)
    # The same 96 runs: the runs of the two families have the same names but for the prefix, and v37's family has as many (the base's family is v37's own under v37)
    same_96_runs = same_96()
    # The repetition of the runs by release 0.5.0 (`results/v38_v0.5.0`): which runs, and whether the table that says they agree says so for each
    rerun_dir = HERE.parents[1] / "analysis/acmg_benchmark/results/v38_v0.5.0"
    rerun = pd.read_csv(rerun_dir / "rerun_comparison.csv") if (rerun_dir / "rerun_comparison.csv").exists() else pd.DataFrame(columns=["kind", "name", "identical"])
    rerun_runs = rerun[rerun.kind == "run"]
    # Release 0.5.0, timed alone on 7 November 2026 (`results/v38_v0.5.0/annotation_time.tsv`): the median of three runs on every core, one run on one thread, and the scoring script
    tt = pd.read_csv(HERE.parents[1] / "analysis/acmg_benchmark/results/v38_v0.5.0/annotation_time.tsv", sep="\t", comment="#")
    fact("takes 25.3 seconds to annotate on every core", f"{tt[tt.threads == 'all'].wall_s.median():.1f} seconds", f"{tt[tt.threads == '1'].wall_s.iloc[0]:.0f} seconds",
         f"{tt[tt.threads == 'score'].wall_s.iloc[0]:.0f} seconds")
    fact("Those times are of release 0.5.0", "7 November 2026", "the median of three runs" if len(tt[tt.threads == "all"]) == 3 else "?", "one run" if len(tt[tt.threads == "1"]) == 1 else "?")
    prov37 = pd.read_csv(DATA / "provenance_v37.csv")
    _dir37 = HERE.parents[1] / "analysis/acmg_benchmark/results/v37_manuscript/provenance"
    made_on_all = {new_york_date(dict(ln.split(": ", 1) for ln in (_dir37 / f"{x}.txt").read_text().splitlines() if ": " in ln).get("recorded_at", "")) for x in prov37.run}
    fact("the annotation of one of the 89 runs of v37 took a median of", str(int(prov37.annotation_seconds.median())), str(len(prov37)), "8 November 2026" if made_on_all == {"2026-11-08"} else "?")
    # 2.3
    v1 = traj.loc[1]
    fact("Every rate is over the 673,660 records of the benchmark", n(v1.n), p(v1.same_direction_rate_of_annotated), p(v1.same_direction_rate))
    fact("We also kept a ledger of the changes", str(len(ledger)))
    # 2.5
    fact("so the 51 residual calls of Section 3.6 are in-sample", str(N_OPP), str(int(run("v37_sw_no_exceptions").opposite)), str(int(run("v37_hk_strict").opposite)))
    # 2.6: the runs
    g = harness.groupby("group").size()
    k = harness.groupby("knob").size()
    fact("For this paper we re-ran every sweep of a threshold that no guideline states", str(int(k["pm2_ad_af_threshold"] - 0)).replace("8", "8"),
         str(int(k["pm2_ar_af_threshold"])), str(int(k["bs1_af_threshold"])), str(int(k["bs2_hom_prevalence_threshold"])), str(int(k["bp7_max_intron_offset"])),
         str(int(k["phylop_conserved"])), str(int(k["bp1_max_pathogenic_missense"])), str(int(g["sweep"])))
    fact("Four more values of the PM2 threshold were run under the 2015 combining table", w(int(g["interaction"])),
         str(sum(1 for x in EARLIER if "prev" in x or x.endswith("off"))), str(sum(1 for x in EARLIER if "ac_" in x)), str(len(EARLIER)))
    fact("To measure what each design decision buys", str(int(g["switch"])), str(int(g["source"])))
    fact("so for 9 further runs we built the binary with a patch", str(int(g["cutoff"])))
    fact("Two departures from the sources that we found while reading them again", str(int(g["hook"]) - 3))
    fact("a fifth run with the same build removes together the three adjustments", "fifth" if int(g["hook"]) == 7 else "?")
    fact("two more runs with that build give BP1's veto a count of 1 and a count of 2", "two more runs" if int(g["hook"]) == 7 else "?")
    fact("Three runs change an input", str(int(g["gap"])))
    fact("With nothing set, each measurement build is the default classifier",
         *(["every cell"] if all(controls[lab("cells_equal_to_v37")] == controls.cells) and all(controls[lab("criteria_firing_equal_to_v37")] == controls.criteria) else ["?"]))
    # One binary for each kind of run and a different binary for each kind: the two measurement builds must not be
    # the plain build under another name, which is what a build made into the target directory of the plain one is.
    kinds = {"plain": prov[~prov.run.str.startswith((rn("v37_cut_"), rn("v37_hk_")))], "cutoff": prov[prov.run.str.startswith(rn("v37_cut_"))],
             "hooks": prov[prov.run.str.startswith(rn("v37_hk_"))]}

    def kind_binary(df) -> str:
        """The binary that made most of the runs of a kind (ties broken by the hash, so that the dump does not depend on the order of a set)."""
        return sorted(df.binary_sha256_prefix.value_counts().items(), key=lambda kv: (-kv[1], kv[0]))[0][0]

    kind_bin = {k: kind_binary(v) for k, v in kinds.items()}
    # The family v37 was made with binaries built from commit 5d8004a, the plain one with the hash that the text gives
    k37_plain = prov37[~prov37.run.str.startswith(("v37_cut_", "v37_hk_"))]
    binaries37 = pd.read_csv(HERE.parents[1] / "analysis/acmg_benchmark/results/v37_manuscript/binaries.tsv", sep="\t")
    built37 = binaries37[binaries37.sha256.str.startswith(kind_binary(k37_plain))]
    fact("The same 96 runs were made for v37, with binaries built from commit 5d8004a", kind_binary(k37_plain),
         "5d8004a" if len(built37) == 1 and "commit 5d8004a" in built37.iloc[0].built_from else "?", n(len(prov37) + len(hold37)) if same_96_runs else "?")
    # Every run of a kind is made by the binary of its kind: a run that another binary made is a failure. (The first drafts had 18 sweeps that the
    # binary of the working tree that became v38 made first and the tag's made again; the runs of this paper were all made by the tag's.)
    for kind, df in kinds.items():
        for name, prefix in zip(df.run, df.binary_sha256_prefix):
            if prefix != kind_bin[kind]:
                failures.append(f"{name} was made by {prefix} and not by the {kind} binary {kind_bin[kind]}")
    if len(set(kind_bin.values())) != 3:
        failures.append(f"the three kinds of run are not made by three different binaries: {kind_bin}")
    # The 42 runs of Section 3.8 were made a second time, by the same binary, into a directory of their own (`results/v38_v0.5.0`): the repetition
    # must name the binary of the stored run, `rerun_comparison.csv` must say that every scored file is the same, and the guard does not rely on that
    # table alone: the concordance matrix and the firing counts of the repetition, which are in the repository (`runs_v38.tsv` of the two directories),
    # must be the stored run's, cell for cell
    stored_rows = pd.read_csv(HERE.parents[1] / "analysis/acmg_benchmark/results/v37_manuscript/runs_v38.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("name")
    again_rows = pd.read_csv(rerun_dir / "runs_v38.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("name")
    cells = [c for c in again_rows.columns if c in stored_rows.columns]
    same_cells = lambda name: name in stored_rows.index and name in again_rows.index and len(cells) > 50 and bool((stored_rows.loc[name, cells] == again_rows.loc[name, cells]).all())  # noqa: E731
    stored_dir = HERE.parents[1] / "analysis/acmg_benchmark/results/v37_manuscript/provenance"
    read_prov = lambda f: dict(line.split(": ", 1) for line in f.read_text().splitlines() if ": " in line) if f.exists() else {}  # noqa: E731
    repeated = []
    for name in rerun_runs.name:
        again, stored = read_prov(rerun_dir / "provenance" / f"{name}.txt"), read_prov(stored_dir / f"{name}.txt")
        rows = rerun[(rerun.kind == "run") & (rerun.name == name)]
        if again.get("fastvep_sha256") and again.get("fastvep_sha256") == stored.get("fastvep_sha256") and len(rows) == 1 and rows.iloc[0].identical == "yes" and same_cells(name):
            repeated.append(name)
        else:
            failures.append(f"{name}: its repetition (results/v38_v0.5.0) was not made by the binary of the stored run, or differs from it in a scored file or in a cell of its matrix and firing counts")
    made_on = {new_york_date(read_prov(stored_dir / f"{x}.txt").get("recorded_at", "")) for x in harness.index}
    again_on = {new_york_date(rerun_recorded_at(x)) for x in repeated}
    fact("The 38 runs of Section 3.8 on the April set and the 4 on the October release were then made a second time", str(len(repeated) - sum(1 for x in repeated if x.startswith("v38_oct_"))),
         str(sum(1 for x in repeated if x.startswith("v38_oct_"))), "made a second time" if len(repeated) == 42 else "?",
         "8 November 2026" if made_on == {"2026-11-08"} and again_on == {"2026-11-08"} else "?")
    fact("All 89 runs (Supplementary Table S2) were made by the same harness", str(len(harness)), str(len(kinds["plain"])),
         str(len(kinds["cutoff"])), str(len(kinds["hooks"])), kind_bin["plain"], "8 November 2026" if made_on == {"2026-11-08"} else "?")
    # 2.8
    nw_n = new_idx["new variants of the release"]
    fact("A new variant is not one that ClinVar did not hold in April", p(new_idx["  in the VCF of 4 April 2026, the file of the April set"] / nw_n),
         "at a lower review status" if new_idx["    at a lower review status than two stars"] == new_idx["  in the VCF of 4 April 2026, the file of the April set"] else "?",
         p(new_idx["  with a record in the April variant_summary"] / nw_n))
    fact("Seven runs of the changes that Section 4.5 proposes were made on the same variants", "Seven" if len(hold) == 7 else "?")
    # 2.7
    fact("The replay of the default scheme reproduces v38 for 673,641", n(veto.replay_matches_the_run), n(veto.variants))
    fact("The remaining 19 are variants that the replay calls Likely pathogenic", str(int(veto.held_at_vus_by_the_veto)))
    fact("It reproduces the table run, the control", "3 and 1")


# --------------------------------------------------------------------------------- Section 3.1


def trajectory() -> None:
    t = traj
    m = mainline
    # v37 and v38 were run again with the complete databases, and the rows of the trajectory are those runs; the series that the text calls flat is the
    # one on the first databases, in which v14 to v36 are the rows and v37 and v38 are read from the comparison of the two stacks
    ovn = pd.read_csv(DATA / "databases_old_vs_new.csv")
    one = lambda v, k, db="first": float(ovn[(ovn.version == v) & (ovn.set == "April set") & (ovn.measure == k)].iloc[0][db])  # noqa: E731
    first_rate = {k: float(t.loc[k].same_direction_rate) for k in m.loc[14:36].index}
    first_rate[37] = one("v37", "same_direction") / N_ALL
    first_rate[13] = float(t.loc[13].same_direction_rate)
    first_rate[38] = one("v38", "same_direction") / N_ALL
    ks = sorted(first_rate)
    steps_first = {k: abs(first_rate[k] - first_rate[j]) for j, k in zip(ks, ks[1:])}
    fact("For 37 of the 38 versions (v3 was a partial run", str(len(t)), str(int(t.index.max())), "v3" if 3 not in t.index else "?")
    fact("Same-direction agreement with ClinVar rose from 51.0% in v1", p(t.loc[1].same_direction_rate), p(first_rate[38]), p(t.loc[38].same_direction_rate),
         "on the first databases and to" if (t.databases[[1, 7, 36]] != "complete").all() and t.loc[38].databases == "complete" else "?")
    fact("Agreement rose 8.8 points between v2 and v4", f"{(t.loc[4].same_direction_rate - t.loc[2].same_direction_rate) * 100:.1f}",
         f"{(t.loc[6].same_direction_rate - t.loc[5].same_direction_rate) * 100:.1f}")
    fact("The first review (v7 to v8) added 2.8 points", f"{(t.loc[8].same_direction_rate - t.loc[7].same_direction_rate) * 100:.1f}")
    fact("The next run, v9, lost 2.7 points", f"{(t.loc[8].same_direction_rate - t.loc[9].same_direction_rate) * 100:.1f}")
    fact("Recall of ClinVar's pathogenic-side variants fell from 63.2% at v8", p(t.loc[8].plp_recall), p(t.loc[9].plp_recall), p(t.loc[38].plp_recall),
         f"{(t.loc[8].plp_recall - t.loc[38].plp_recall) * 100:.1f}")
    fact("Agreement then rose by 3.2 points at v10", f"{(t.loc[10].same_direction_rate - t.loc[9].same_direction_rate) * 100:.1f}",
         f"{(t.loc[11].same_direction_rate - t.loc[10].same_direction_rate) * 100:.1f}",
         f"{(t.loc[13].same_direction_rate - t.loc[12].same_direction_rate) * 100:.1f}")
    after = [first_rate[k] for k in ks if 14 <= k <= 37]
    step = max(steps_first[k] for k in ks if 14 <= k <= 37)
    fact("Across the next 23 runs (v14 to v37, without v17)", str(len(after)), p(min(after), 2), p(max(after), 2), f"{step * 100:.2f}", "on the first databases")
    # v38 is the one step after v13 that is larger than the largest of the 23: 0.22 points, to 77.84%, and the complete databases then give v37 and v38 their own rates
    fact("The 38th run, v38, which reads the mode of inheritance that ClinGen curates", p(first_rate[38], 2), f"{(first_rate[38] - first_rate[37]) * 100:.2f} points",
         "the one step after v13 of more than" if sum(1 for k in steps_first if steps_first[k] > step) == 1 and max(steps_first, key=steps_first.get) == 38 else "?", f"{step * 100:.2f}")
    fact("Repeated with the complete databases, v37 and v38 reach 79.07% and 79.32%", p(t.loc[37].same_direction_rate, 2), p(t.loc[38].same_direction_rate, 2),
         f"{(t.loc[38].same_direction_rate - first_rate[38]) * 100:.1f} points more" if f"{(t.loc[38].same_direction_rate - first_rate[38]) * 100:.1f}" == f"{(t.loc[37].same_direction_rate - first_rate[37]) * 100:.1f}" else "?",
         "more than every step since v13" if (t.loc[37].same_direction_rate - first_rate[37]) > max(steps_first.values()) else "?")
    fact("Opposite-direction calls fell from 94 at v13 to 49 at v37 and 48 at v38", str(int(t.loc[13].opposite)), n(one("v37", "opposite")), n(one("v38", "opposite")),
         str(int(t.loc[13].false_benign)), n(one("v37", "false_benign")), n(one("v38", "false_benign")), str(int(t.loc[13].false_pathogenic)), n(one("v37", "false_pathogenic")),
         n(one("v38", "false_pathogenic")), "at both" if one("v37", "false_benign") == one("v38", "false_benign") else "?",
         "on the first databases", str(int(t.loc[37].opposite)), str(int(t.loc[38].opposite)), str(int(t.loc[38].false_benign)),
         "in both" if t.loc[37].false_benign == t.loc[38].false_benign else "?")
    fact("Opposite-direction calls had reached a minimum of 25 at v10", str(int(t.loc[10].opposite)), str(int(t.loc[11].opposite)))
    rate = lambda r: r.opposite / r.committed * 100  # noqa: E731
    fact("The matrix stored for v1 holds 31 opposite-direction calls", str(int(t.loc[1].opposite)), n(t.loc[1].committed), n(t.loc[38].committed),
         f"{rate(t.loc[1]):.3f}%", f"{rate(t.loc[37]):.3f}%", f"{rate(t.loc[38]):.3f}%", f"{min(rate(t.loc[i]) for i in range(7, 39)):.3f}%",
         "in v37 and in v38" if f"{rate(t.loc[37]):.3f}" == f"{rate(t.loc[38]):.3f}" else "?")
    blb_first = [float(t.loc[k].blb_recall) for k in m.loc[13:36].index] + [one("v37", "blb_recall")]
    fact("Recall of ClinVar's pathogenic-side variants in direction has stayed between",
         f"{m.loc[11:38].plp_recall.min() * 100:.1f}", f"{m.loc[11:38].plp_recall.max() * 100:.1f}",
         f"{min(blb_first) * 100:.1f}", f"{max(blb_first) * 100:.1f}", p(one("v38", "blb_recall"), 1), p(t.loc[37].blb_recall, 1), p(t.loc[38].blb_recall, 1))


# --------------------------------------------------------------------------------- Section 3.2


def review() -> None:
    f = ledger.found_by.value_counts()
    fact("The ledger of Section 2.3 holds", str(len(ledger)))
    fact("Expert review first found 30 of the 64", str(f["expert_review"]), str(f["code_audit"]), str(f["discordance_mining"]), str(f["first_run"]), str(f["benchmark_aggregate"]),
         str(f["sweep"]), str(f["vep_comparison"]), str(f["regression"]), f"{f['expert_review'] / len(ledger) * 100:.0f}%")
    c = ledger["class"].value_counts()
    fact("By class, 21 items were a condition of a guideline", str(c["G"]))
    fact("Seventeen were wiring defects", "Seventeen" if c["W"] == 17 else str(c["W"]))
    fact("Eight concerned how far a population frequency", "Eight" if c["F"] == 8 else str(c["F"]),
         "five were annotation" if c["A"] == 5 else str(c["A"]), "five were thresholds" if c["T"] == 5 else str(c["T"]),
         "four concerned scoring" if c["S"] == 4 else str(c["S"]), "two were circularity" if c["C"] == 2 else str(c["C"]),
         "two were combining" if c["R"] == 2 else str(c["R"]))
    r1, r2, r6 = rounds.loc[1], rounds.loc[2], rounds.loc[6]
    fact("In round 1 the workbook held 101", str(int(r1.rows)), str(int(r1.clinvar_supported)), str(int(r1.fastvep_supported)), str(int(r1.neither_supported)), str(int(r1.no_verdict)))
    fact("Round 2 held 122 calls from v8", str(int(r2.rows)), str(int(r2.clinvar_supported)), str(int(r2.fastvep_supported)), str(int(r2.neither_supported)),
         str(int(r2.needs_evidence)), str(int(r2.no_verdict)))
    fact("Seven of the 65", w(int(r2.homology_marked_rows)), w(int(r2.homology_marked_on_fastvep_call)), str(int(r2.clinvar_supported - r2.homology_marked_on_fastvep_call)),
         str(int(r2.no_verdict + r2.homology_marked_on_fastvep_call)), str(int(r2.clinvar_supported)), str(int(r2.no_verdict)))
    fact("Round 6 held 64 calls from v26", str(int(r6.rows)), str(int(r6.clinvar_supported)), str(int(r6.fastvep_supported)), str(int(r6.neither_supported)), str(int(r6.needs_evidence)))
    fact("So when the software and ClinVar disagreed in direction", p(r1.clinvar_share_of_verdicts, 0), p(r2.clinvar_share_of_verdicts, 0), p(r6.clinvar_share_of_verdicts, 0),
         str(int(r1.rows_with_verdict)), str(int(r2.rows_with_verdict)), str(int(r6.rows_with_verdict)))


# --------------------------------------------------------------------------------- Section 3.3


def thresholds() -> None:
    pm2 = runs[(runs.knob == "pm2_ad_af_threshold") & (runs.group == "sweep")]
    lo = min(pm2.plp_recall.min(), base.plp_recall)
    hi = max(pm2.plp_recall.max(), base.plp_recall)
    fact("Across strict absence and every threshold from 1e-6 to 5e-4", p(lo, 2), p(hi, 2), n((hi - lo) * 93_812), "93,812")
    pm2_37 = runs37[(runs37.knob == "pm2_ad_af_threshold") & (runs37.group == "sweep")]
    span37 = max(pm2_37.plp_recall.max(), base37.plp_recall) - min(pm2_37.plp_recall.min(), base37.plp_recall)
    fact("On v38 the same sweep is almost flat, and on v37 it was almost flat as well", "almost flat" if (hi - lo) < 0.01 and span37 < 0.01 else "?",
         "on v37 it was almost flat as well" if span37 < 0.01 else "?")
    fact("It reaches 62.47% at the default 4e-5", p(base.plp_recall, 2),
         str(int(run("v37_sw_pm2ad_0.0001").same_pathogenic - base.same_pathogenic)), str(int(run("v37_sw_pm2ad_0.0002").same_pathogenic - run("v37_sw_pm2ad_0.0001").same_pathogenic)))
    fact("Reference-discordant pathogenic calls rise from 21 at strict absence", str(int(run("v37_sw_pm2ad_0").false_pathogenic)), str(int(base.false_pathogenic)),
         str(int(run("v37_sw_pm2ad_0.0002").false_pathogenic)), n(run("v37_sw_pm2ad_0").same_direction), n(base.same_direction))
    t = {v: run(f"v37_ix_table_pm2ad_{v}") for v in ("0", "1e-06", "1e-05", "0.0001")}
    t4 = run("v37_sw_table_both_sides")
    fact("Recall of pathogenic-side variants is 51.0% at strict absence", p(t["0"].plp_recall), p(t["1e-06"].plp_recall), p(t["1e-05"].plp_recall), p(t4.plp_recall), p(t["0.0001"].plp_recall))
    fact("Raising the threshold from strict absence to 4e-5 gains 5,699", n(t4.same_pathogenic - t["0"].same_pathogenic), n(t4.false_pathogenic - t["0"].false_pathogenic),
         f"{(t4.same_pathogenic - t['0'].same_pathogenic) / 93_812 * 100:.1f} points")
    # v10's sweep under this table was 37.8% at strict absence and 56.8% at 4e-5 (RUN_VERSIONS.md); v37 reproduced it, and v38 has its steps
    # and a third of its rise (the sweep of Section 3.8 finds the knee at 1e-5 under both)
    def rises_every_step(rs, version: str) -> bool:
        """The five settings of the sweep under the 2015 table, in the order of the threshold: each step gains correct pathogenic-side calls."""
        ladder = [rs.loc[f"{version}_ix_table_pm2ad_0"], rs.loc[f"{version}_ix_table_pm2ad_1e-06"], rs.loc[f"{version}_ix_table_pm2ad_1e-05"],
                  rs.loc[f"{version}_sw_table_both_sides"], rs.loc[f"{version}_ix_table_pm2ad_0.0001"]]
        return all(b_.same_pathogenic > a_.same_pathogenic for a_, b_ in zip(ladder, ladder[1:]))

    t37 = {v: runs37.loc[f"v37_ix_table_pm2ad_{v}"] for v in ("0",)}
    t4_37 = runs37.loc["v37_sw_table_both_sides"]
    fact("The sweep that v10 recorded under this table", "37.8%", "56.8%", p(t37["0"].plp_recall), p(t4_37.plp_recall), n(t4_37.same_pathogenic - t37["0"].same_pathogenic),
         f"{(t4_37.same_pathogenic - t37['0'].same_pathogenic) / 93_812 * 100:.1f} points", p(t["0"].plp_recall), p(t4.plp_recall),
         "reproduced by v37" if abs(t37["0"].plp_recall - 0.378) < 0.011 and abs(t4_37.plp_recall - 0.568) < 0.002 else "?",
         n(t4.same_pathogenic - t["0"].same_pathogenic), f"{(t4.same_pathogenic - t['0'].same_pathogenic) / 93_812 * 100:.1f} points",
         "rises at every step as v37 does" if rises_every_step(runs37, "v37") and rises_every_step(rv38, "v38") else "?",
         "by about a third as much" if 0.30 <= (t4.same_pathogenic - t["0"].same_pathogenic) / (t4_37.same_pathogenic - t37["0"].same_pathogenic) <= 0.37 else "?")
    gain_ix, gain_pm = t4.same_pathogenic - t["0"].same_pathogenic, hi - lo
    fact("It is still the steeper sweep: it moves recall 12 times as far", f"{gain_ix / (gain_pm * 93_812):.0f} times as far", n(gain_ix), n(gain_pm * 93_812))
    ix_runs = runs[runs.group == "interaction"]

    def no_ix_run_sets_the_recessive_threshold(ixs) -> bool:
        """The configurations of the interaction runs (the 2015 table with one PM2 threshold) set the threshold of the other genes and not the recessive one."""
        cfgs = [HERE.parents[1] / "analysis/acmg_benchmark/results/v37_manuscript/configs" / f"{name}.toml" for name in ixs.index]
        keys = [set(re.findall(r"^([a-z0-9_]+)\s*=", c.read_text(), re.M)) for c in cfgs if c.exists()]
        return len(ixs) > 0 and len(keys) == len(ixs) and all("pm2_ad_af_threshold" in k and "pm2_ar_af_threshold" not in k for k in keys)
    fact("Why v38 flattens it is a reading that the firing counts support and that no run has tested", n(t["0"].fire_PM2_Supporting), n(runs37.loc["v37_ix_table_pm2ad_0"].fire_PM2_Supporting),
         n(fv("inheritance_read_recessive_only")), "7e-5",
         "no run sweeps both PM2 thresholds under the 2015 table" if no_ix_run_sets_the_recessive_threshold(ix_runs) else "?")
    fact("Of the 51,703 pathogenic-side calls", n(pm2dep.loc["Richards 2015 table", "calls_depending_on_pm2"]), n(pm2dep.loc["Richards 2015 table", "of_which_pvs1_and_pm2_only"]),
         p(pm2dep.loc["Richards 2015 table", "of_which_pvs1_and_pm2_only"] / pm2dep.loc["Richards 2015 table", "calls_depending_on_pm2"], 0))
    fact("The point scale of the default scheme counts a lone PVS1 as 8 points", n(pm2dep.loc["shipped: SVI points on the pathogenic side", "calls_depending_on_pm2"]))

    def step(knob, i):
        sub = steps[(steps.knob == knob) & steps.shipped_between].reset_index(drop=True)
        return sub.iloc[i]

    def rate(row):
        return row.gain_calls / row.added_errors if row.added_errors else float("inf")

    bp7 = steps[steps.knob == "BP7 far boundary"].reset_index(drop=True)
    r = {(a, b): x for a, b, x in zip(bp7.from_value, bp7.to_value, bp7.itertuples())}
    fact("Moving it from 200 to 300 nucleotides buys 163", n(r[(200, 300)].gain_calls), n(r[(300, 500)].gain_calls), n(r[(500, 1000)].gain_calls), n(r[(1000, 1000000)].gain_calls),
         w(r[(1000, 1000000)].added_errors))
    fact("The default 300 is at the corner", n(rate(r[(0, 50)])), n(rate(r[(200, 300)])), n(rate(r[(300, 500)])), n(rate(r[(500, 1000)])))
    fact("Turning the extension off is not free either", n(r[(0, 50)].gain_calls), str(r[(0, 50)].added_errors))
    bs2 = steps[steps.knob == "BS2 homozygote bar"].reset_index(drop=True)
    fact("The steps from 1e-3 to 1e-6 are worth", *[n(rate(x)) for x in bs2.iloc[1:].itertuples()])
    bs1 = steps[steps.knob == "BS1 bar"].reset_index(drop=True)
    fact("Between 2e-2 and 5e-3 (default), lowering the threshold buys", n(bs1.iloc[0].gain_calls + bs1.iloc[1].gain_calls), str(int(bs1.iloc[0].added_errors + bs1.iloc[1].added_errors)))
    fact("Below the default value the exchange rate is 858", n(rate(bs1.iloc[2])), n(rate(bs1.iloc[3])), n(rate(bs1.iloc[4])), n(rate(bs1.iloc[5])))
    b3 = run("v37_sw_bs1_0.001")
    fact("By the knee rule that fixed BP7's bound", p((b3.blb_recall - base.blb_recall), 1).rstrip("%"), n(b3.same_benign - base.same_benign), str(int(b3.false_benign - base.false_benign)),
         str(int(base.false_benign)), str(int(b3.false_benign)))
    fact("In the table of ClinGen VCEP specifications we extracted", f"{vcep.loc['BS1', 'median'] * 1e4:g}e-4".replace("2.5e-4", "2.5e-4"), n(vcep.loc["BS1", "genes_with_a_bar"]),
         n(vcep.loc["BS1", "genes_tighter_than_default"]), n(vcep.loc["BA1", "genes_with_a_bar"]), n(vcep.loc["BA1", "genes_tighter_than_default"]))
    fact("The median PM2 threshold is 4e-5", n(vcep.loc["PM2", "genes_with_a_bar"]), "1.3e-6", "5e-4")
    # The BS2 carrier threshold: a trade on the genes that v38 reads as dominant-only, and flat on the 50 that v37 read (the runs v37_bs2* and their
    # twins are the earlier runs of the carrier settings, paired with the sweeps of Section 3.8)
    cs = {k: run(f"v37_bs2{k}") for k in ("off", "prev_3e-5", "prev_1e-4", "prev_1e-3", "prev_3e-3")}
    cf = {k: run(f"v37_bs2ac_{k}") for k in (20, 100, 1000)}
    c37 = {k: runs37.loc[f"v37_bs2{k}"] for k in ("off", "prev_3e-5", "prev_1e-4", "prev_1e-3", "prev_3e-3")}
    sweep_flat = lambda rs, b: len({(r.false_benign, r.false_pathogenic) for r in (*rs.values(), b)}) == 1  # noqa: E731
    fact("On the 746 genes that v38 reads as dominant-only the threshold is not flat", n(fv("inheritance_read_dominant_only")), "not flat" if not sweep_flat(cs, base) else "?")
    fact("Switching it off, so that an allele count of at least five alone suffices", str(int(cs["off"].false_benign)), str(int(cs["off"].false_pathogenic)),
         str(int(cs["prev_3e-5"].false_benign)), str(int(cs["prev_1e-4"].false_benign)), str(int(base.false_benign)),
         "at 1e-3 and 3e-3" if cs["prev_1e-3"].false_benign == cs["prev_3e-3"].false_benign == base.false_benign else "?")
    above = [cs["prev_1e-3"].same_benign - cs["prev_3e-3"].same_benign, base.same_benign - cs["prev_1e-3"].same_benign]
    below = [(a, b) for a, b in ((base, cs["prev_1e-4"]), (cs["prev_1e-4"], cs["prev_3e-5"]), (cs["prev_3e-5"], cs["off"]))]
    fact("From 3e-3 to the default each step adds benign-side calls and no error", n(above[0]), n(above[1]),
         "no error" if cs["prev_1e-3"].false_benign == cs["prev_3e-3"].false_benign == base.false_benign and cs["prev_1e-3"].false_pathogenic == base.false_pathogenic else "?",
         *[n((b.same_benign - a.same_benign) / (b.false_benign - a.false_benign)) if b.false_benign != a.false_benign else "no added error" for a, b in below])
    fact("The default 3e-4 recovers 439 pathogenic-side calls", n(cs["off"].same_pathogenic), n(base.same_pathogenic), n(base.same_pathogenic - cs["off"].same_pathogenic),
         n(cs["off"].same_benign), n(base.same_benign), n(cs["off"].same_benign - base.same_benign), n(cs["prev_1e-3"].same_pathogenic - base.same_pathogenic),
         n(cs["prev_1e-3"].same_pathogenic))
    fact("With the threshold off, a floor of 20 alleles gives 39", *[str(int(cf[k].false_benign)) for k in (20, 100, 1000)])
    fact("On v37 the same sweep was flat, because the threshold reached the 50 genes", n(fv37("inheritance_read_dominant_only")), str(int(base37.false_benign)), str(int(base37.false_pathogenic)),
         "was flat" if sweep_flat(c37, base37) else "?")
    # the PM2 threshold of recessive-only genes
    pr = runs[(runs.knob == "pm2_ar_af_threshold") & (runs.group == "sweep")]
    rlo, rhi = min(pr.plp_recall.min(), base.plp_recall), max(pr.plp_recall.max(), base.plp_recall)
    fact("The threshold for recessive-only genes, which the code attributed", p(rlo, 2).rstrip("%"), p(rhi, 2).rstrip("%"),
         str(int(max(pr.same_direction.max(), base.same_direction) - min(pr.same_direction.min(), base.same_direction))),
         "do not change" if pr.false_benign.nunique() == 1 and pr.false_pathogenic.nunique() == 1 and pr.false_benign.iloc[0] == base.false_benign else "?")
    # The cell of Table 3 for the same sweep, which `build_tables.py` computes on its own: it left the default out of the span and read 26 where the text said 27
    fact("PM2 threshold, genes read as recessive-only: 7e-5", f"recall {p(rlo, 2)} to {p(rhi, 2)}",
         f"within {int(max(pr.same_direction.max(), base.same_direction) - min(pr.same_direction.min(), base.same_direction))} calls")
    # BP7's conservation threshold
    ph = steps[steps.knob == "BP7 PhyloP bar"].reset_index(drop=True)
    pv = {(a, b): x for a, b, x in zip(ph.from_value, ph.to_value, ph.itertuples())}
    off, wk, sp = run("v37_sw_phylop_1000"), run("v37_hk_bp7_walker"), run("v37_hk_bp7_spliceai01")
    fact("Raising the threshold from 1 to 1.5 gains 5,814", n(pv[(1.0, 1.5)].gain_calls), n(pv[(1.5, 2.0)].gain_calls), str(int(pv[(1.5, 2.0)].added_errors)),
         n(rate(pv[(1.5, 2.0)])))
    fact("Beyond the default value, raising it to 3 gains 7,018", n(pv[(2.0, 3.0)].gain_calls), n(pv[(3.0, 4.0)].gain_calls), n(pv[(4.0, 1000.0)].gain_calls))
    fact("Every step buys at least 1,600 calls per added reference-discordant benign call", n(min(rate(x) for x in pv.values() if x.added_errors > 0) // 100 * 100))
    fact("Switching the test off, which follows the caution of Walker et al.", n(off.same_benign - base.same_benign), n(off.vus_to_blb - base.vus_to_blb),
         n(off.same_direction - base.same_direction), pp((off.same_direction - base.same_direction) / N_ALL, 1), str(int(base.false_benign)), str(int(off.false_benign)))
    # 313 is the SpliceAI threshold alone against v37 and 402 is the same threshold on top of the switched-off conservation test; an earlier
    # draft called 313 an "addition" to the switch, which the numbers pass and the meaning does not, so both are checked here
    no_error_change = (sp.false_benign == base.false_benign and sp.false_pathogenic == base.false_pathogenic
                       and wk.false_benign == off.false_benign and wk.false_pathogenic == off.false_pathogenic)
    fact("Asking that SpliceAI be at most 0.1, as BP4 does", n(base.same_benign - sp.same_benign), n(off.same_benign - wk.same_benign),
         "neither changing an error" if no_error_change else "?", n(wk.same_direction), p(wk.same_direction / N_ALL, 1),
         n(wk.same_direction - base.same_direction), f"{w(int(wk.false_benign - base.false_benign))[1]} additional")
    _walker = hooknew[hooknew.run == rn("v37_hk_bp7_walker")]
    _walker_genes = list(_walker.gene)
    _wk = run("v37_hk_bp7_walker")
    _pm2 = [g_ for g_, c_ in zip(_walker.gene, _walker.criteria_in_run) if "PM2" in c_]
    _bp4 = [g_ for g_, c_ in zip(_walker.gene, _walker.criteria_in_run) if "PM2" not in c_ and "BP4" in c_]
    fact(f"We have not changed the classifier: the {w(len(_walker_genes))[1]} variants", *_walker_genes, n(round(_wk.same_benign - base.same_benign, -3)),
         f"{' and '.join(_pm2)} have PM2 and BP4 beside BP7 and {' and '.join(_bp4)} has BP4" if len(_pm2) + len(_bp4) == len(_walker_genes) and _bp4 else "?",
         "ClinVar-pathogenic" if (_walker.clinvar == "P/LP").all() else "?", f"{w(int(_wk.false_benign - base.false_benign))[1]} reference-discordant benign calls")
    fact("of the 2,924 genes of the ClinGen file it reads 1,310 as recessive-only", n(fv("clingen_genes_in_file")), n(fv("inheritance_read_recessive_only")), n(fv("inheritance_read_dominant_only")),
         n(fv("inheritance_read_both")), n(fv("inheritance_read_neither")), n(fv37("inheritance_read_recessive_only")), n(fv37("inheritance_read_dominant_only")),
         n(fv37("inheritance_read_both")), n(fv37("inheritance_read_neither")))
    fact("The two thresholds therefore act on nearly every gene and on the 1,310 genes read as recessive-only", n(fv("inheritance_read_recessive_only")))
    fact("On v37 the second threshold reached 67 genes, and a flat sweep there said", n(fv37("inheritance_read_recessive_only")))
    fact("That median leaves out the panels that give no single number", n(vaud_pm2.sum()), n(vaud_pm2["ok"]), n(vaud_pm2["absence-required"]), n(vaud_pm2["no-threshold"]),
         n(vaud_pm2["inheritance-conditioned"]), n(vaud_pm2["superscript-notation"]), n(vaud_pm2["conflicting-panels"]),
         n(vaud_pm2["ambiguous"] + vaud_pm2["not-specified"] + vaud_pm2.get("implausible", 0)), "7e-5")
    fact("Across BA1, BS1 and PM2 together the extraction read 411", n(vaud.sum()), n(vaud["ok"]))
    fact(f"rests on {w(int(off.false_benign - base.false_benign))[1]} events", n(round((off.same_benign - base.same_benign) / (off.false_benign - base.false_benign), -2)))
    pw = prop.loc["BP7 as Walker 2023 describes it"]
    fact("On the 116,078 variants that were not in the April set, BP7 as the recommendation describes it adds", n(pw.new_delta_same_direction),
         "one reference-discordant benign call" if pw.new_delta_false_benign == 1 else "?")
    b1c = {k: run(f"v37_sw_bp1veto_{k}") for k in (1, 2, 5, 10, 25)}
    bp1_steps = steps[steps.knob == "BP1 veto count"]
    fact("the veto itself avoids 28 reference-discordant benign calls relative to none", str(int(run("v37_sw_no_bp1_veto").false_benign - base.false_benign)))
    fact("At a count of 1, BP1 is withheld in every gene", n(b1c[1].same_direction - base.same_direction), str(int(base.false_benign)), str(int(b1c[1].false_benign)),
         str(int(b1c[1].false_pathogenic - base.false_pathogenic) or "do not change") if b1c[1].false_pathogenic == base.false_pathogenic else "?",
         n(base.same_benign - b1c[1].same_benign), n(base.vus_to_blb - b1c[1].vus_to_blb))
    rs = csv("bp1_veto_count_resolved_v37.csv")
    both = sorted(rs[rs.run == rn("v37_sw_bp1veto_1")].gene)
    fact("The two reference-discordant benign calls that a count of 1 removes are the BP1 pair of Section 3.6", *(both if both == ["IQSEC1", "PRR12"] else ["?"]),
         "one pathogenic missense record for PRR12" if int(bp1_idx.loc["PRR12", "pathogenic_missense_records"]) == 1 and bp1_idx.loc["PRR12", "own_record_among_them"] else "?",
         "two for IQSEC1" if int(bp1_idx.loc["IQSEC1", "pathogenic_missense_records"]) == 2 and bp1_idx.loc["IQSEC1", "own_record_among_them"] else "?",
         bp1_idx.loc["IQSEC1", "the_others"])
    fact("The rise is those 1,864 variants staying VUS", p(c5.loc["VUS"].sum() / N_ALL, 0), n(base.same_benign - b1c[1].same_benign), n(base.vus_to_blb - b1c[1].vus_to_blb))
    lv1, lv2 = run("v37_hk_bp1_loo_veto1"), run("v37_hk_bp1_loo_veto2")
    rl = csv("bp1_veto_count_resolved_v37.csv")
    fact("With the variant's own record left out of the count as well (Section 3.4)", str(int(lv1.false_benign)), n(lv1.same_direction - base.same_direction), n(base.same_benign - lv1.same_benign),
         str(int(lv2.false_benign)),
         *(["IQSEC1"] if list(rl[rl.run == rn("v37_hk_bp1_loo_veto1")].gene) == ["IQSEC1"] and len(rl[rl.run == rn("v37_hk_bp1_loo_veto2")]) == 0 else ["?"]))
    fact("At a count of 2 without that change the figures are", n(b1c[2].same_direction - base.same_direction), str(int(b1c[2].false_benign)), n(base.same_benign - b1c[2].same_benign),
         *(sorted(rs[rs.run == rn("v37_sw_bp1veto_2")].gene) if list(rs[rs.run == rn("v37_sw_bp1veto_2")].gene) == ["IQSEC1"] else ["?"]))
    fact("a count of 1 would withhold BP1 for both", *both, "two and one" if (int(bp1_idx.loc["IQSEC1", "pathogenic_missense_records"]), int(bp1_idx.loc["PRR12", "pathogenic_missense_records"])) == (2, 1) else "?")
    fact("Above the default value each step buys correct benign-side calls", str(int(b1c[5].false_benign - base.false_benign)), str(int(b1c[10].false_benign - base.false_benign)),
         str(int(b1c[25].false_benign - base.false_benign)))
    fact("No step gives more than 111 calls per added reference-discordant benign call", n(bp1_steps[bp1_steps.added_errors > 0].calls_per_error.max()), "the default count of 3")


def points() -> None:
    ntb = pts[pts.scope == "no benign evidence"].set_index("pathogenic_points")
    fact("Variants with 0 or 1 point", n(ntb.loc[0, "n"]), n(ntb.loc[1, "n"]), p(ntb.loc[0, "plp_fraction_of_all"]), p(ntb.loc[1, "plp_fraction_of_all"]))
    fact("At 2 points the share is 10.7%", p(ntb.loc[2, "plp_fraction_of_all"]), p(ntb.loc[3, "plp_fraction_of_all"]), p(ntb.loc[4, "plp_fraction_of_all"]))
    fact("At 5 points it jumps to 64.1%", p(ntb.loc[5, "plp_fraction_of_all"]), n(ntb.loc[5, "clinvar_plp"]), n(ntb.loc[5, "n"]), p(ntb.loc[5, "clinvar_vus"] / ntb.loc[5, "n"]),
         p(ntb.loc[6, "plp_fraction_of_all"]), p(min(ntb.loc[7:10, "plp_fraction_of_all"])), p(ntb.loc[10, "plp_fraction_of_all"]))
    c5r, c4r, c7r, c8r, c9r = (run(f"v37_cut_lp{i}") for i in (5, 4, 7, 8, 9))
    fact("Moving the Likely pathogenic cutoff from 6 to 5 points", p(base.plp_recall), p(c5r.plp_recall), n(c5r.same_pathogenic - base.same_pathogenic))
    fact("It does so at the price of 6,003 more ClinVar-VUS variants", n(c5r.vus_to_plp - base.vus_to_plp), str(int(c5r.false_pathogenic - base.false_pathogenic)),
         f"{(c5r.same_pathogenic - base.same_pathogenic) / ((c5r.vus_to_plp - base.vus_to_plp) + (c5r.false_pathogenic - base.false_pathogenic)):.1f}")
    fact("Agreement in direction rises by 4,673 calls", n(c5r.same_direction - base.same_direction), n(base.exact - c5r.exact))
    fact("A cutoff of 4 reaches 76.7% recall", p(c4r.plp_recall), n(c4r.vus_to_plp), n(c4r.false_pathogenic))
    fact("Raising the cutoff to 7 or 8 gives up", n(base.same_pathogenic - c7r.same_pathogenic), n(base.same_pathogenic - c8r.same_pathogenic), str(int(base.false_pathogenic - c7r.false_pathogenic)),
         p(c9r.plp_recall), str(int(c9r.false_pathogenic)))
    fact("No cutoff changes the number of reference-discordant benign calls", str(int(base.false_benign)),
         "every run" if all(run(f"v37_cut_{k}").false_benign == base.false_benign for k in ("lp4", "lp5", "lp7", "lp8", "lp9", "p9", "p11", "p12")) else "?")
    fact("Tavtigian's framework places a posterior of 0.90 at 6 points", "0.90", "0.81")
    p9, p11, p12 = run("v37_cut_p9"), run("v37_cut_p11"), run("v37_cut_p12")
    fact("The Pathogenic cutoff of 10 points does not move agreement in direction", str(int(max(abs(x.same_direction - base.same_direction) for x in (p9, p11, p12)))))
    fact("It does move exact agreement", n(p9.exact), n(base.exact))
    fact("Of the 58,201 Likely pathogenic calls", n(fv("lp_calls")), n(fv("lp_calls_at_9_points")), p(fv("lp_calls_at_9_points") / fv("lp_calls"), 0),
         n(fv("lp_at_9_points_pvs1_pm2_only")), n(fv("lp_at_9_points_clinvar_pathogenic")))
    fact("Scoring both sides on the point scale with Likely benign from -1", p(base.blb_recall), p(rep.loc["SVI points, both sides, Likely Benign from -1 (Tavtigian)", "blb_recall"]),
         str(int(base.false_benign)), n(rep.loc["SVI points, both sides, Likely Benign from -1 (Tavtigian)", "false_benign"]),
         n(csv("full_points_false_benign_signatures_v37.csv").iloc[0].n))
    fact("Tightening the floor to -2 gives 50.6% recall", p(rep.loc["SVI points, both sides, Likely Benign from -2", "blb_recall"]), n(rep.loc["SVI points, both sides, Likely Benign from -2", "false_benign"]),
         p(rep.loc["SVI points, both sides, Likely Benign from -3", "blb_recall"]), n(rep.loc["SVI points, both sides, Likely Benign from -3", "false_benign"]))
    m1 = rep.loc["SVI points, both sides, Likely Benign from -1 (Tavtigian)"]
    fact("The table therefore does better than the floors of -2 and -3", str(int(m1.false_benign - base.false_benign)), f"{(m1.blb_recall - base.blb_recall) * 100:.1f}")
    fact("Scoring both sides on the point scale with Likely benign from -1", p(base.vus_recall), p(m1.vus_recall))
    # The three schemes side by side: the default, the 2015 table on both sides (a real run, which its replay reproduces) and the
    # point scale on both sides (a replay only), each from the replay of Appendix Table A2.
    dflt, tbl = rep.loc[lab("shipped (v37)")], rep.loc["Richards 2015 table, both sides"]
    t4s = run("v37_sw_table_both_sides")
    pts_rows = [rep.loc[s] for s in ("SVI points, both sides, Likely Benign from -1 (Tavtigian)", "SVI points, both sides, Likely Benign from -2",
                                     "SVI points, both sides, Likely Benign from -3")]
    fact("The first two are real runs, the second through the setting",
         "the second through the setting" if (dflt.exact == base.exact and dflt.same_direction == base.same_direction and t4s.exact == tbl.exact
                                              and t4s.same_direction == tbl.same_direction and t4s.plp_recall == tbl.plp_recall
                                              and val.loc[rn("v37_sw_table_both_sides"), "cells_differ"] == 0) else "?")
    fact("The default scheme has 431,727 exact calls", n(dflt.exact), n(dflt.same_direction), p(dflt.plp_recall), p(dflt.blb_recall),
         str(int(dflt.false_benign)), str(int(dflt.false_pathogenic)))
    fact("The 2015 table on both sides has 430,652", n(tbl.exact), n(tbl.same_direction), p(tbl.plp_recall), p(tbl.blb_recall), str(int(tbl.false_benign)),
         str(int(tbl.false_pathogenic)), n(dflt.same_direction - tbl.same_direction), pp(dflt.plp_recall - tbl.plp_recall),
         str(int(dflt.false_pathogenic - tbl.false_pathogenic)),
         "on the pathogenic side only" if (dflt.same_benign == tbl.same_benign and dflt.false_benign == tbl.false_benign) else "?")
    fact("These three point scales are replays of the recorded criteria", *[n(r.exact) for r in pts_rows], *[n(r.same_direction) for r in pts_rows],
         "fewer exact calls" if all(r.exact < dflt.exact and r.same_direction < dflt.same_direction for r in pts_rows) else "?")
    fact("Removing it raises same-direction agreement by 684", n(run("v37_sw_no_an_floor").same_direction - base.same_direction), "no opposite-direction call" if run("v37_sw_no_an_floor").opposite == base.opposite else "?")


# --------------------------------------------------------------------------------- Section 3.4


def ablations() -> None:
    d = lambda name, col="same_direction": run(name)[col] - base[col]
    ci, cb = run("v37_sw_clinvar_informed"), run("v37_src_noclinvar_protein")
    fact("Letting it count, as the runs before v9 did", n(d("v37_sw_clinvar_informed")), pp(d("v37_sw_clinvar_informed") / N_ALL, 2), p(base.plp_recall), p(ci.plp_recall),
         str(int(base.opposite - ci.opposite)))
    bl = run("v37_hk_bp1_loo")
    fact("The veto itself counts the variant's own record", str(int(base.false_benign)), str(int(bl.false_benign)), str(int(bl.false_benign - base.false_benign)),
         str(int(run("v37_sw_no_bp1_veto").false_benign - base.false_benign)), *list(hooknew[hooknew.run == rn("v37_hk_bp1_loo")].gene))
    fact("Removing the ClinVar protein and splice index", n(-d("v37_src_noclinvar_protein")), pp(d("v37_src_noclinvar_protein") / N_ALL, 2), f"{(base.plp_recall - cb.plp_recall) * 100:.1f}")
    fact("Removing the allele-level ClinVar record as well", "one further call" if abs(d("v37_src_noclinvar_all") - d("v37_src_noclinvar_protein")) == 1 else "?")
    fact("The dependence on ClinVar-derived inputs is therefore small", str(int(base.false_benign)), str(int(cb.false_benign)), str(int(cb.false_benign - base.false_benign)))
    fact("The dependence on ClinVar-derived inputs is therefore small in agreement and larger in errors", str(int(run("v37_src_noclinvar_all").false_benign)),
         str(int(run("v37_src_noclinvar_all").opposite)), str(int(base.opposite)))
    stc = run("v37_hk_strict")
    alone = [run(x).false_benign - base.false_benign for x in (rn("v37_hk_bp1_loo"), rn("v37_sw_no_lowpen_gate"), rn("v37_sw_no_exceptions"))]
    fact("Two further adjustments also read the variant's own ClinVar record", str(int(alone[1])), str(int(alone[2])))
    fact("Removed together with the veto's own-record count", str(int(base.false_benign)), str(int(stc.false_benign)), str(int(base.opposite)), str(int(stc.opposite)),
         *[str(int(a)) for a in alone], "exactly the sum" if stc.false_benign - base.false_benign == sum(alone) else "?")
    fact("Three of the ten added calls are Benign calls", w(fv("strict_opposite_called_benign")), w(int(stc.opposite - base.opposite)),
         *str(facts["strict_opposite_called_benign_genes"]).replace("G6PD,G6PD", "G6PD").split(","), "two of the 15 curated alleles (both G6PD)"
         if str(facts["strict_opposite_called_benign_genes"]).split(",").count("G6PD") == 2 else "?")
    if fv("strict_opposite_called_pathogenic") != 0:
        failures.append("the strict run has a Pathogenic call among its contradictions, which Section 3.4 and the abstract say it has not")
    sa, ph, gn = run("v37_src_nospliceai"), run("v37_src_nophylop"), run("v37_src_nognomad")
    fact("Removing SpliceAI lowers agreement in direction by 118,039", n(-d("v37_src_nospliceai")), pp(d("v37_src_nospliceai") / N_ALL, 1), p(base.blb_recall), p(sa.blb_recall))
    fact("Removing PhyloP lowers it by 84,759", n(-d("v37_src_nophylop")), pp(d("v37_src_nophylop") / N_ALL, 1), p(ph.blb_recall), n(-d("v37_src_nognomad")), pp(d("v37_src_nognomad") / N_ALL, 1), p(gn.blb_recall))
    fact("Without SpliceAI, 219,372 of the 268,243 BP4 firings are lost", n(base.fire_BP4 - sa.fire_BP4), n(base.fire_BP4), n(base.fire_BP7))
    cg = run("v37_src_noclingen")
    pvs1_any = lambda r_: sum(r_[c] for c in r_.index if c.startswith("fire_PVS1"))
    fact("On the pathogenic side, removing the ClinGen gene-disease file", n(-d("v37_src_noclingen")), pp(d("v37_src_noclingen") / N_ALL, 1), p(base.plp_recall), p(cg.plp_recall),
         n(pvs1_any(base) - pvs1_any(cg)), n(pvs1_any(base)), p((pvs1_any(base) - pvs1_any(cg)) / pvs1_any(base), 0))
    nc = run("v37_src_noconstraint")
    fact("Removing gnomAD gene constraint takes out PP2 and BP1 entirely", n(pvs1_any(base) - pvs1_any(nc)), f"{(base.plp_recall - nc.plp_recall) * 100:.1f}", n(d("v37_src_noconstraint")),
         p(c5.loc["VUS"].sum() / N_ALL, 0))
    bs1s = steps[steps.knob == "BS1 bar"].reset_index(drop=True)
    fact("Their prices are all below the 858 calls", n(bs1s.iloc[2].gain_calls / bs1s.iloc[2].added_errors))
    fact("Removing REVEL, the source specific to missense variants", pp(d("v37_src_norevel") / N_ALL, 2), n(-d("v37_src_norepeat")))
    fact("Removing SpliceAI or PhyloP lowers reference-discordant benign calls from 24 to 4", str(int(base.false_benign)), str(int(sa.false_benign)), str(int(ph.false_benign)))
    fact("BP7 is a criterion of high yield and high cost", n(fv("bp7_firings_clinvar_blb")), n(fv("bp7_firings")), p(fv("bp7_firings_clinvar_blb") / fv("bp7_firings"), 1),
         n(-d("v37_src_nophylop")), str(int(base.false_benign - ph.false_benign)), n(round(-d("v37_src_nophylop") / (base.false_benign - ph.false_benign), -3)))
    nb, ng, nm, ne, nl = (run(k) for k in (rn("v37_sw_no_bp1_veto"), rn("v37_sw_no_gene_validity"), rn("v37_sw_no_mechanism_gate"), rn("v37_sw_no_exceptions"), rn("v37_sw_no_lowpen_gate")))
    fact("Four safeguards remove errors and cost no calls", f"or {int(base.false_benign - run('v37_hk_bp1_loo_veto1').false_benign)} with the own record left out", n(-d("v37_sw_no_bp1_veto")), str(int(nb.false_benign - base.false_benign)), n(-d("v37_sw_no_mechanism_gate")),
         str(int(nm.false_pathogenic - base.false_pathogenic)), str(int(ne.false_benign - base.false_benign)), str(int(nl.false_benign - base.false_benign)))
    fact("Reading absence from gnomAD as rarity for PM2 gains 1,499", n(-d("v37_sw_pm2_needs_record")))
    nh, nf = run("v37_sw_no_homology_list"), run("v37_sw_no_faf")
    fact("Three more remove errors at a price in calls", n(d("v37_sw_no_gene_validity")), str(int(ng.false_pathogenic - base.false_pathogenic)),
         n(d("v37_sw_no_gene_validity") / (ng.false_pathogenic - base.false_pathogenic)), n(d("v37_sw_no_homology_list")), str(int(nh.false_benign - base.false_benign)),
         n(d("v37_sw_no_homology_list") / (nh.false_benign - base.false_benign)), n(d("v37_sw_no_faf")), str(int(nf.false_benign - base.false_benign)),
         str(int(base.false_pathogenic - nf.false_pathogenic)), n(d("v37_sw_no_faf") / (nf.false_benign - base.false_benign)))
    fact("The gnomAD region flags", n(d("v37_sw_no_region_flags")), w(int(run("v37_sw_no_region_flags").opposite - base.opposite)),
         n(d("v37_sw_no_an_floor")), n(d("v37_sw_no_pm1_benign_half")))
    fact("Two opt-in settings change little", n(d("v37_sw_pp3_cap_moderate")), f"{(base.plp_recall - run('v37_sw_pp3_cap_moderate').plp_recall) * 100:.2f}", n(-d("v37_sw_pvs1_nmd50")))
    pm = run("v37_sw_pm2_moderate")
    fact("PM2 at Moderate, as in the 2015 standards", f"{(pm.plp_recall - base.plp_recall) * 100:.1f}", n(pm.same_pathogenic - base.same_pathogenic), n(d("v37_sw_pm2_moderate")),
         str(int(pm.false_pathogenic - base.false_pathogenic)), n(d("v37_sw_pm2_moderate") / (pm.false_pathogenic - base.false_pathogenic)))
    tb = run("v37_sw_table_both_sides")
    fact("The 2015 table on both sides loses 4,900 calls", n(-d("v37_sw_table_both_sides")), f"{(base.plp_recall - tb.plp_recall) * 100:.1f}", str(int(base.false_pathogenic - tb.false_pathogenic)))


# --------------------------------------------------------------------------------- Section 3.5


def discordance() -> None:
    pr_n = int(c5.loc["Pathogenic"].sum() + c5.loc["Likely_pathogenic"].sum())
    bl_n = int(c5.loc["Likely_benign"].sum() + c5.loc["Benign"].sum())
    same = cell("P/LP", "P/LP") + cell("VUS", "VUS") + cell("B/LB", "B/LB")
    ex, sd = hl("exact_match"), head[head.metric == "same_direction"].iloc[0]
    fact("On the benchmark snapshot, v38 reproduces ClinVar's class for", n(hlc("exact").value), p(ex.value), f"{ex.ci_lo * 100:.1f}", f"{ex.ci_hi * 100:.1f}",
         n(hlc("same_direction").value), p(sd.value), f"{sd.ci_lo * 100:.1f}", f"{sd.ci_hi * 100:.1f}", n(hlc("genes").value), n(hlc("genes").value - 1),
         str(int(fv("records_without_gene_symbol"))))
    not_same = N_ALL - same
    fact("Of the 139,325 records that are not in the same direction", n(not_same), n(cell("P/LP", "VUS") + cell("B/LB", "VUS")),
         p((cell("P/LP", "VUS") + cell("B/LB", "VUS")) / not_same), n(cell("VUS", "P/LP") + cell("VUS", "B/LB")), p((cell("VUS", "P/LP") + cell("VUS", "B/LB")) / not_same),
         n(hl("variants_no_call").value), n(hl("opposite").value), p(hl("opposite").value / not_same, 2))
    o = hl("opposite_rate")
    fact("The opposite-direction rate is 0.0076% of records", f"{o.value * 100:.4f}%", f"{o.ci_lo * 100:.4f}", f"{o.ci_hi * 100:.4f}", n(hlc("false_benign").value), n(hlc("false_pathogenic").value),
         f"{hlc('false_benign').value / pr_n:.3%}", f"{hlc('false_pathogenic').value / bl_n:.3%}".replace("0.009%", "0.009%"), n(pr_n), n(bl_n))
    tier = int(c5.loc["Pathogenic", "Likely_pathogenic"] + c5.loc["Likely_pathogenic", "Pathogenic"] + c5.loc["Likely_benign", "Benign"] + c5.loc["Benign", "Likely_benign"])
    fact("A further 102,608 records (15.2%)", n(tier), p(tier / N_ALL), n(c5.loc["Pathogenic", "Likely_pathogenic"] + c5.loc["Benign", "Likely_benign"]), n(c5.loc["Pathogenic", "Likely_pathogenic"]), n(c5.loc["Benign", "Likely_benign"]))
    pp_ = hl("precision_plp_calls")
    pb_ = hl("precision_blb_calls")
    fact("Of the 62,240 variants called pathogenic-side", n(pp_.denominator), p(pp_.value), n(pp_.numerator), f"{pp_.ci_lo * 100:.1f}", f"{pp_.ci_hi * 100:.1f}",
         p(cell("VUS", "P/LP") / pp_.denominator), n(pb_.denominator), p(pb_.value), n(pb_.numerator), f"{pb_.ci_lo * 100:.1f}", f"{pb_.ci_hi * 100:.1f}")
    pr, br, vr = hl("plp_recall"), hl("blb_recall"), hl("vus_recall")
    fact("Recall is much lower than precision on both sides", p(pr.value), f"{pr.ci_lo * 100:.1f}", f"{pr.ci_hi * 100:.1f}", p(br.value), f"{br.ci_lo * 100:.1f}", f"{br.ci_hi * 100:.1f}",
         p(vr.value), f"{vr.ci_lo * 100:.1f}", f"{vr.ci_hi * 100:.1f}")
    fact("Among the 1,273 genes with at least ten", n(fv("genes_with_10_plp")), n(fv("genes_with_10_plp_none_called")), p(fv("genes_with_10_plp_none_called") / fv("genes_with_10_plp"), 0),
         n(fv("genes_with_10_plp_90pct_called")), p(fv("genes_with_10_plp_90pct_called") / fv("genes_with_10_plp"), 0))
    # the baselines and the tiers of the call
    b_all, b1, b2, b3 = (baseline.loc[k] for k in ("all VUS", "R1", "R1 + R2", "R1 + R2 + R3"))
    v37r = baseline.loc[lab("v37 (the classifier)")]
    fact("Calling every variant VUS agrees with ClinVar in direction for 43.8%", p(b_all.same_direction_share))
    r5 = baseline.loc["R1 at 0.005"]
    cv_ = lambda rule, bar: curve[(curve.rule == rule) & (curve.benign_frequency_bar == bar)].iloc[0]
    r1_1e2, r1_2e2 = cv_("R1", 0.01), cv_("R1", 0.02)
    fact("Calling a variant benign-side when its population-maximum allele frequency is 5e-3 or more", p(r5.same_direction_share), n(r5.opposite),
         p(b1.same_direction_share), n(b1.opposite), n(r1_1e2.opposite), p(r1_1e2.same_direction_share), n(r1_2e2.opposite), p(r1_2e2.same_direction_share),
         f"about {(v37r.same_direction_share - r1_1e2.same_direction_share) * 100:.0f} points less")
    fact("Adding a pathogenic-side call for a stop-gained", p(b2.same_direction_share), n(b2.opposite), p(b3.same_direction_share),
         f"{(b3.same_direction_share - v37r.same_direction_share) * 100:.1f}", n(b3.opposite), n(b3.false_benign), n(b3.false_pathogenic), f"{b3.opposite / v37r.opposite:.0f}")
    r3c = curve[curve.rule == "R1 + R2 + R3"]
    r23 = curve[curve.rule == "R2 + R3"].iloc[0]
    fact("Without the frequency rule, the last two rules give", p(r23.same_direction_share), f"{(r23.same_direction_share - v37r.same_direction_share) * 100:.1f} points above {BASE}", n(r23.opposite),
         n(r23.false_benign), n(r23.false_pathogenic), f"{r23.opposite / v37r.opposite:.0f}",
         f"no threshold of the first rule brings the three rules below {int(r23.opposite)} contradictions" if r3c.opposite.min() >= r23.opposite else "?")
    fact("The last rule is not a bargain by the standard of Section 2.5 either", n(b3.blb_calls_that_clinvar_calls_blb - b2.blb_calls_that_clinvar_calls_blb),
         n(b3.false_benign - b2.false_benign), f"{(b3.blb_calls_that_clinvar_calls_blb - b2.blb_calls_that_clinvar_calls_blb) / (b3.false_benign - b2.false_benign):.0f}")
    fact("The population maximum used by these rules is a point estimate over ten groups", str(int(run("v37_sw_no_faf").false_benign - base.false_benign)),
         f"{n(r1_2e2.opposite)}, {n(r1_1e2.opposite)}, {n(r5.opposite)} and {n(b1.opposite)} at 2e-2, 1e-2, 5e-3 and 1e-3")
    fact("The rules are not worse than v38 at everything", p(b2.plp_recall), p(v37r.plp_recall), p(b2.plp_concordant_share), p(v37r.plp_concordant_share),
         n(b2.false_pathogenic), n(v37r.false_pathogenic), n(r2v[1]), n(r2v[3]), n(r2v[2]), n(r2v[5]), n(r2v[4]))
    fact("Its 4,039 Pathogenic calls are ClinVar pathogenic-side", n(tiers.loc["Pathogenic", "n"]), n(tiers.loc["Pathogenic", "clinvar_same_direction"]),
         p(tiers.loc["Pathogenic", "clinvar_same_direction"] / tiers.loc["Pathogenic", "n"]), n(tiers.loc["Pathogenic", "clinvar_vus"]),
         n(tiers.loc["Benign", "n"]), n(tiers.loc["Benign", "clinvar_same_direction"]), p(tiers.loc["Benign", "clinvar_same_direction"] / tiers.loc["Benign", "n"]),
         n(tiers.loc["Benign", "clinvar_vus"]), n(tiers.loc["Likely_pathogenic", "clinvar_opposite"]), n(tiers.loc["Likely_benign", "clinvar_opposite"]))
    st = lambda s, a, b: stars[(stars.stars == s) & (stars.clinvar == a) & (stars.fastvep == b)].iloc[0]
    fact("The 3-star records, which are expert-panel classifications", p(st(3, "P/LP", "P/LP").row_fraction), p(st(2, "P/LP", "P/LP").row_fraction),
         p(st(3, "B/LB", "B/LB").row_fraction), p(st(2, "B/LB", "B/LB").row_fraction), p(st(3, "B/LB", "VUS").row_fraction))
    fact("All 23 four-star records", str(int(st(4, "P/LP", "VUS").row_total)), str(int(st(4, "P/LP", "VUS").n)))
    cg = lambda g, a, b: cons[(cons.consequence_group == g) & (cons.clinvar == a) & (cons.fastvep == b)].iloc[0]
    null_ok = cg("Frameshift or nonsense", "P/LP", "P/LP").n + cg("Canonical splice site", "P/LP", "P/LP").n
    fact("Variants called pathogenic-side are almost all null variants", n(null_ok), p(null_ok / base.same_pathogenic, 1), n(cg("Missense", "P/LP", "P/LP").n))
    fact("Among ClinVar's pathogenic-side variants, 88.4%", p(cg("Canonical splice site", "P/LP", "P/LP").row_fraction), p(cg("Frameshift or nonsense", "P/LP", "P/LP").row_fraction),
         p(cg("Missense", "P/LP", "P/LP").row_fraction), p(cg("In-frame or protein-altering", "P/LP", "P/LP").row_fraction))
    fact("Among ClinVar's benign-side variants the picture is the reverse", p(cg("Intronic", "B/LB", "B/LB").row_fraction), p(cg("Synonymous", "B/LB", "B/LB").row_fraction),
         p(cg("UTR", "B/LB", "B/LB").row_fraction), p(cg("Missense", "B/LB", "B/LB").row_fraction))
    pv = plp_pts["n"]
    fact("By the evidence that was met, 10,694", n(pv[5]), p(pv[5] / pv.sum()), n(pv[0] + pv[1]), p((pv[0] + pv[1]) / pv.sum()), n(pv[2] + pv[3] + pv[4]), str(int(pv[pv.index >= 6].sum())))
    pc = plp_cons.groupby("consequence_group").n.sum()
    tot = int(pv.sum())
    fact("By consequence, 18,992", n(pc["Missense"]), p(pc["Missense"] / tot), n(pc["Frameshift or nonsense"]), p(pc["Frameshift or nonsense"] / tot),
         n(pc["Other (flank, non-coding, unscored)"]), p(pc["Other (flank, non-coding, unscored)"] / tot))
    fact("The concentration in genes is high", p(conc.loc[10, "share_of_plp_to_vus"]), p(conc.loc[100, "share_of_plp_to_vus"]), p(conc.loc[500, "share_of_plp_to_vus"]), n(conc.loc[10, "genes_with_any"]))
    top = gene.sort_values("plp_called_vus", ascending=False).head(10)
    fact("The concentration in genes is high", *[f"{g} {int(v)}" for g, v in zip(top.gene, top.plp_called_vus)])
    fact("Of the 52,584 ClinVar pathogenic-side stop-gained and frameshift variants", n(pvs1["ClinVar P/LP stop-gained and frameshift variants"]), n(pvs1["with a PVS1 criterion (any strength)"]),
         n(pvs1["without PVS1"]), n(pvs1["without PVS1, genes"]))
    fact("Those 5,184 are 57% of the 9,136", p(pvs1["without PVS1"] / fv("trunc_vus"), 0), n(fv("trunc_vus")))
    gr = lambda key: gap_reasons[gap_reasons.reason.str.contains(key, regex=False)].iloc[0]
    fact("For 4,696 of the 5,184 variants", n(gr("no route").variants), p(gr("no route").variants / pvs1["without PVS1"], 0), n(gr("no route").genes))
    fact("For 463 (in 38 genes) the validity gate withholds PVS1", n(gr("validity gate").variants), n(gr("validity gate").genes), n(gr("mechanism gate").variants),
         n(gr("mechanism gate").genes), n(gr("(other)").variants))
    fact("All 76 variants in genes that gnomAD constrains", n(pvs1["without PVS1, variants in genes constrained by gnomAD"]),
         *[f"{g} {int(v)}" for g, v in gap_const.variants.items()], n(gr("validity gate").variants_in_constrained_genes), n(gr("mechanism gate").variants_in_constrained_genes))
    fact("The only ClinGen curation of RASA1 in the file", "Noonan syndrome", "Disputed", f"{int(gap_const.loc['RASA1', 'variants_naming_it'])} of its {int(gap_const.loc['RASA1', 'variants'])}",
         f"{int(gap_const.loc['LRP2', 'variants_naming_it'])} of its {int(gap_const.loc['LRP2', 'variants'])}")
    fact("and 5,184, in 625 genes", n(pvs1["without PVS1, genes"]), n(pvs1["without PVS1, genes by the symbol of the annotation"]),
         n(pvs1["without PVS1, variants whose ClinVar gene differs from the symbol of the annotation"]))
    fact("The second route is nearly empty", n(fv("mechanism_table_entries")), n(fv("mechanism_table_gof")),
         "one is dominant negative" if fv("mechanism_table_dominant_negative") == 1 else "?",
         "five name both" if fv("mechanism_table_lof_and_gof") == 5 else "?",
         "those five are each licensed" if fv("mechanism_lof_genes_licensed_anyway") == fv("mechanism_table_lof_and_gof") else "?")
    fact("The validity gate had one leak that changed no call in v37", n(fv("clingen_weak_only_genes")), facts["clingen_weak_only_genes_with_semicolon_label"],
         "no variant in the benchmark" if fv("clingen_weak_only_semicolon_genes_in_benchmark") == 0 else "?")
    fact("CFTR holds 439 of the 5,184 variants", n(pvs1["without PVS1, in CFTR (stop-gained and frameshift)"]), n(pvs1["CFTR null variants in the benchmark (stop, frameshift, canonical splice)"]),
         str(int(st(4, "P/LP", "VUS").n)), str(int(st(4, "P/LP", "VUS").row_total)))
    c1, c2, c3 = closure.loc[rn("v37_gap_clingen_complete")], closure.loc[rn("v37_gap_orphanet_lof")], closure.loc[rn("v37_gap_both")]
    versions = {h.split(":")[0] for h in cftr_null.hgvsc}
    mane = n(pvs1["  of them null on the MANE Select transcript, the one the pick chose"])
    fact("All 558 are null on CFTR's MANE Select transcript", mane, *versions, "the transcript that the pick chose" if len(versions) == 1 and (cftr_null.on_mane_select == "yes").all() else "?")
    cv = cftr_null.clinvar.value_counts()
    fact("ClinVar classifies 556 of the 558 as pathogenic-side", n(pvs1["  of them ClinVar pathogenic-side"]), n(cftr_null.shape[0]), n(cv["Pathogenic"]), n(cv["Likely_pathogenic"]), n(cv["VUS"]),
         "every one of them" if cftr_null.shape[0] == pvs1["CFTR null variants in the benchmark (stop, frameshift, canonical splice)"] else "?")
    assert (cftr_null[lab("v37_call")] == "VUS").all() and int(cftr_null.complete_clingen_call.isin(["Pathogenic", "Likely_pathogenic"]).sum()) == int(c1.cftr_null_called_plp)
    fact("ClinGen's export of 11 October 2026 lists CFTR as Definitive for cystic fibrosis", n(missing[LACKED]), n(missing["  in genes"]),
         n(missing["  of Definitive, Strong or Moderate strength"]), n(missing["  genes with such a curation"]), n(missing[LACKED] - 1))
    fact("The completed file adds the 112 rows of the October export", n(missing[ROWS]), n(missing[LACKED]), "seven" if missing[RELABEL] == 7 else "?",
         "one row" if missing[EXPORT_DAY_ROW] == 1 else "?", n(missing[DAY_GENE_N]), "changes none of the" if missing[DAY_GENE_CALLS] == 0 else "?")
    fact("It does change the criteria of 6 of them under v38", n(missing[DAY_GENE_CHANGED]), n(missing[DAY_GENE_LOST[0]]),
         "AR alone" if missing[DAY_GENE_TAGS[0]] == missing[DAY_GENE_N] else "?", "AD and AR" if missing[DAY_GENE_TAGS[1]] == missing[DAY_GENE_N] else "?",
         f"{missing[DAY_GENE_LOST[1]] * 1e5:.1f}e-5", f"{missing[DAY_GENE_LOST[2]] * 1e5:.1f}e-5",
         "from 7e-5 to 4e-5", "lose PM2 without changing a call" if missing[DAY_GENE_CALLS] == 0 and missing[DAY_GENE_LOST[0]] == missing[DAY_GENE_CHANGED] else "?")
    gg = closure_genes[closure_genes.run == rn("v37_gap_clingen_complete")].set_index("gene").gap_variants_gaining_pvs1
    fact("Putting all 112 back, in a run that differs from v38 in that file alone", n(missing[ROWS]), n(c1.gap_variants_gaining_pvs1), n(c1.genes_gaining_pvs1), n(c1.gap_variants),
         p(c1.gap_variants_gaining_pvs1 / c1.gap_variants, 0), n(c1.gap_variants_called_plp), *[f"{g} {int(gg[g])}" for g in ("CFTR", "ELP1", "EVC", "SPG7")],
         f"which are {n(sum(gg[g] for g in ('CFTR', 'ELP1', 'EVC', 'SPG7')))} of the {n(c1.gap_variants_gaining_pvs1)}")
    fact("Of CFTR's 558 null variants, 556 then carry PVS1", n(c1.cftr_null_with_pvs1), n(c1.cftr_null_called_plp), n(c1.four_star_called_plp), n(c1.four_star),
         str(int(star_cell(4, "P/LP", "VUS").row_total - star_cell(4, "P/LP", "VUS").n)))
    fp = closure_fp[closure_fp.run == rn("v37_gap_clingen_complete")]
    # the sentence names the variant: one under v38 (ALPK1, ClinVar Benign), where v37 had two (ALPK1 again, ClinVar Benign and Likely benign)
    only = fp.iloc[0] if len(fp) == 1 and len(fp) == int(c1.blb_to_plp) else None
    fact("The run adds 852 correct pathogenic-side calls", n(c1.plp_vus_to_plp), p(base.plp_recall), p(run("v37_gap_clingen_complete").plp_recall), n(c1.vus_to_plp),
         f"{n(c1.blb_to_plp)} reference-discordant pathogenic call, a {only.consequence.replace('_', ' ')} of {only.gene} that ClinVar calls {only.clinvar}, "
         f"at an allele frequency of {only.popmax_af * 1e4:.1f}e-4" if only is not None else "?")
    fp2 = closure_fp[closure_fp.run == rn("v37_gap_orphanet_lof")]
    fact("Licensing PVS1 for the 1,042 genes that Orphanet labels", str(run("v37_gap_orphanet_lof").value).replace(" genes", "").replace("1042", "1,042"), n(c2.gap_variants_gaining_pvs1),
         n(c2.gap_variants_called_plp), n(c2.plp_vus_to_plp), n(c2.vus_to_plp), n(c2.blb_to_plp), *sorted(set(fp2.gene)),
         f"{fp2[fp2.consequence == 'splice_acceptor_variant'].popmax_af.min() * 100:.1f}%")
    bo = run("v37_gap_both")
    fact("With the complete ClinGen file as well", n(c3.gap_variants_gaining_pvs1), n(c3.gap_variants_called_plp), n(bo.same_direction - base.same_direction),
         pp((bo.same_direction - base.same_direction) / N_ALL, 1), p(bo.plp_recall), str(int(bo.false_pathogenic - base.false_pathogenic)),
         n((bo.same_direction - base.same_direction) / (bo.false_pathogenic - base.false_pathogenic)))
    fact("Even so 3,470 of the 5,184 variants", n(c3.gap_variants - c3.gap_variants_gaining_pvs1), p((c3.gap_variants - c3.gap_variants_gaining_pvs1) / c3.gap_variants, 0))
    fact("The other 3,952 truncating variants", n(fv("trunc_vus_with_pvs1_strong") + fv("trunc_vus_with_pvs1_moderate") + fv("trunc_vus_with_pvs1_supporting") + fv("trunc_vus_with_pvs1")),
         n(fv("trunc_vus_with_pvs1_strong")), n(fv("trunc_vus_with_pvs1_moderate")), n(fv("trunc_vus_with_pvs1_supporting")), f"{n(fv('trunc_vus_with_pvs1'))} at full strength",
         n(fv("canonical_splice_vus") - fv("canonical_splice_vus_with_pvs1")), n(fv("canonical_splice_vus")))
    sg = lambda a, b, k: sigs[(sigs.clinvar == a) & (sigs.fastvep == b)].set_index("signature")
    bv = sg("B/LB", "VUS", 0)
    fact("At least 48,291 of these variants (51%) carry PM2", n(bv.loc["PM2_Supporting", "n"] + bv.loc["BP4+PM2_Supporting", "n"]),
         p((bv.loc["PM2_Supporting", "n"] + bv.loc["BP4+PM2_Supporting", "n"]) / cell("B/LB", "VUS"), 0), n(bv.loc["BP4", "n"]), n(bv.loc["BS1", "n"]), n(bv.loc["(none)", "n"]))
    fact("Of the ClinVar benign-side variants that are called VUS, 24,771", n(bv.loc["PM2_Supporting", "n"]), p(bv.loc["PM2_Supporting", "fraction_of_cell"]), n(bv.loc["BP4+PM2_Supporting", "n"]),
         p(bv.loc["BP4+PM2_Supporting", "fraction_of_cell"]), n(bv.loc["BP4", "n"]), p(bv.loc["BP4", "fraction_of_cell"]), n(bv.loc["(none)", "n"]), p(bv.loc["(none)", "fraction_of_cell"]),
         n(bv.loc["BS1", "n"]), p(bv.loc["BS1", "fraction_of_cell"]))
    vp = sg("VUS", "P/LP", 0)
    fact("The pathogenic-side calls are PP2 with PM2", n(vp.loc["PM2_Supporting+PP2+PP3_Strong", "n"]), n(vp.loc["PM2_Supporting+PVS1", "n"]), n(vp.loc["PM1+PM2_Supporting+PP2+PP3_Moderate", "n"]),
         n(vp.loc["PM1+PM2_Supporting+PM5+PP2", "n"]), n(vp.loc["PVS1", "n"]))
    vb = sg("VUS", "B/LB", 0)
    bp1 = sum(sg("VUS", "B/LB", 0).loc[k, "n"] for k in ("BP1+BP4_Moderate+PM2_Supporting", "BP1+BP4+PM2_Supporting", "BP1+BP4_Moderate", "BP1+BP4"))
    vus_blb = int(fv("vus_called_blb"))
    fact("The benign-side calls are chiefly BP1 with BP4", n(bp1), n(vus_blb), n(cell("VUS", "B/LB")))
    # a benign-side call needs two or more benign criteria or one at Very Strong or Stand-alone strength (the 2015 table), and a pathogenic-side one 6 points
    one_vs = int(fv("vus_called_blb_one_benign_criterion_very_strong_or_standalone"))
    fact("Each pathogenic-side call stands on 6 or more points", n(one_vs), n(vus_blb),
         "two or more benign criteria or on one at Very Strong or Stand-alone strength"
         if fv("vus_called_blb_no_benign_criterion") == 0 and fv("vus_called_blb_one_benign_criterion") == one_vs else "?",
         "6 or more points" if vp.pathogenic_points.min() >= 6 else "?")


# --------------------------------------------------------------------------------- Section 3.6


def residual() -> None:
    fact("The last review round (v26) held 64 opposite-direction calls", str(len(r6)), str(int(r6.still_opposite.astype(str).eq("True").sum())))
    fixed = r6[r6.still_opposite.astype(str) == "False"]
    fact("The 16 that were repaired all became VUS", str(len(fixed)), str(int((fixed.ruling == "clinvar").sum())))
    kept = r6[r6.still_opposite.astype(str) == "True"]
    sr = kept.seen_in_round.value_counts()
    fact("Each of the 48 had already been in the table of an earlier round", str(sr["round 3"]), str(sr["round 4"]), str(sr["round 5"]))
    tot = res_by_group[["clinvar", "fastvep", "neither", "needs_evidence", "disorder_specific"]].sum()
    fact("For 29 of the 48 (60%)", str(int(tot.clinvar)), p(tot.clinvar / N_REV, 0), str(int(tot.fastvep)), p(tot.fastvep / N_REV, 0), str(int(tot.neither)), str(int(tot.needs_evidence)), str(int(tot.disorder_specific)))
    fact("So for 29 of the 48 the reviewer judged ClinVar's call the better", str(int(tot.clinvar)), str(N_REV - int(tot.clinvar)))
    r49 = rules64[(rules64.set == f"the {N_REV} still opposite in {BASE}") & (rules64.ruling == "clinvar")].set_index("rule").rows
    fact("Of the 29 rulings for ClinVar, 10 follow from a class", str(int(r49.sum())), str(int(r49["class"])), str(int(r49["stated"])), str(int(r49["evidence"])))
    left = int((aprl.october_status == "left the 2-star+ set").sum())
    fact("Five of the 51 (10%) are no longer in the 2-star+ set", w(left), p(left / N_OPP, 0), p((N_ALL - sum(temporal[temporal.stratum == x].iloc[0].stratum_n for x in ("unchanged", "tier_change", "direction_change"))) / N_ALL, 1))
    fact("Ten of these are a lone PVS1 and 11 are PVS1 with PM2", "Ten" if fv("pvs1_benign_lone_pvs1") == 10 else "?", "11" if fv("pvs1_benign_pvs1_pm2") == 11 else "?")
    fact("**PVS1 on a variant ClinVar calls benign or likely benign (", f"({int(res_by_group.loc['PVS1 on a ClinVar benign-side variant', 'total'])})")
    fact("Five are seen in 611 to 4,892 gnomAD alleles", "Five" if fv("pvs1_benign_common_n") == 5 else "?", n(fv("pvs1_benign_common_ac_min")), n(fv("pvs1_benign_common_ac_max")),
         f"{fv('pvs1_benign_common_af_min') * 1e4:.1f}e-4", f"{fv('pvs1_benign_common_af_max') * 1e3:.1f}e-3", *str(fv("pvs1_benign_common_genes")).split(","),
         p(fv("pvs1_benign_common_popmax_min"), 2), p(fv("pvs1_benign_common_popmax_max"), 1))
    common_genes = str(fv("pvs1_benign_common_genes")).split(",")
    asked = r6[r6.gene.isin(common_genes) & r6.still_opposite.astype(str).eq("True")]
    fact("BS1, which tests the filtering allele frequency against 5e-3", "did not fire for any of the five" if fv("pvs1_benign_common_bs1_fired") == 0 else "?",
         n(fv("pvs1_benign_common_popmax_min") / 2.5e-4), n(fv("pvs1_benign_common_popmax_max") / 2.5e-4),
         "the reviewer supported ClinVar on all five" if len(asked) == len(common_genes) == 5 and (asked.ruling == "clinvar").all() else "?")
    fact("A disease-specific frequency threshold could act on about seven of them", w(int(fv("pvs1_benign_common_n")) + 2)[1] if int(fv("pvs1_benign_common_n")) + 2 == 7 else "?",
         w(int(fv("pvs1_benign_common_n")))[1], "BS1" if int(res_by_group.loc["BS1 on a ClinVar pathogenic-side variant", "total"] - res_by_group.loc["BS1 on a ClinVar pathogenic-side variant", "not reviewed"]) == 2 else "?")
    fact(f"{w(int(fv('bp4bp7_intron_variant')))[2]} are intronic, {w(int(fv('bp4bp7_synonymous_variant')))[1]} synonymous and {w(int(fv('bp4bp7_splice_region_variant')))[1]} in the splice region")


# --------------------------------------------------------------------------------- Section 3.7


def later_release() -> None:
    sub = lambda name: temporal[temporal.stratum == name].iloc[0]
    kept = int(sub("unchanged").stratum_n + sub("tier_change").stratum_n + sub("direction_change").stratum_n)  # the April variants that are still at two stars or better
    fact("ClinVar's release of 4 October 2026 holds 781,135", n(sub("all").stratum_n), n(kept), n(sub("new").stratum_n))
    fact("A further 8,603 of the April variants have left the set", n(N_ALL - kept))
    ni = new_idx
    nn_ = ni["new variants of the release"]
    vcf_in, vcf_single = ni["  in the VCF of 4 April 2026, the file of the April set"], ni["    a single-submitter record in that file"]
    two_star = ni["  with a record of two stars or better in the April variant_summary"]
    fact("Most of the new variants were already in ClinVar in April, at a lower review status", n(vcf_in), n(nn_), p(vcf_in / nn_), n(vcf_single), p(vcf_single / nn_, 0),
         "none of them at two stars or better" if ni["    at a lower review status than two stars"] == vcf_in else "?")
    fact("In the `variant_summary` of 29 April, from which the protein index was built", n(ni["  with a record in the April variant_summary"]), p(ni["  with a record in the April variant_summary"] / nn_),
         n(ni["    a single-submitter record"]), p(ni["    a single-submitter record"] / nn_, 0), n(two_star), p(two_star / nn_),
         n(ni["    of those, in the VCF of 4 April at a lower review status"]), n(ni["    of those, not in the VCF of 4 April"]),
         f"{(__import__('datetime').date(2026, 4, 29) - __import__('datetime').date(2026, 4, 4)).days} days")
    fact("They are variants that reached two stars since", n(ni["  with a pathogenic or likely pathogenic record in the April variant_summary"]))
    bp1_can = "new ClinVar pathogenic or likely pathogenic missense variants in genes where BP1 can fire (pLI at least 0.9, missense z below 2)"
    own_in = "  with their own residue change among the gene's pathogenic missense records in the April protein index"
    apr_can = "April ClinVar pathogenic or likely pathogenic missense variants in genes where BP1 can fire"
    apr_own = "  April set: with their own residue change among the gene's pathogenic missense records in the April protein index"
    fact("so the protein index lacks every record added since April", n(ni[own_in]), n(ni[bp1_can]), n(ni[bp1_can] - ni[own_in]), n(ni[apr_own]), n(ni[apr_can]))
    fa = sub("all")
    fact("Across all 781,135 records of the release", n(fa.stratum_n), p(fa.same_direction / fa.stratum_n), str(int(fa.opposite)), str(int(fa.false_benign)), str(int(fa.false_pathogenic)))
    fu = sub("unchanged")
    fact("On the 662,668 variants whose class did not change", n(fu.stratum_n), p(fu.same_direction / fu.stratum_n), str(int(fu.opposite)))
    t = lambda g, c: tclass[(tclass.group == g) & (tclass.clinvar == c)].iloc[0]
    new = temporal[temporal.stratum == "new"]
    rec = lambda a, b: new[(new.clinvar == a) & (new.fastvep == b)].iloc[0].row_fraction
    fact("On the 116,078 new variants agreement in direction is 72.6%", n(new.iloc[0].stratum_n), p(new.iloc[0].same_direction / new.iloc[0].stratum_n), p(rec("P/LP", "P/LP")),
         p(rec("VUS", "VUS")), p(rec("B/LB", "B/LB")), p(hl("plp_recall").value), p(hl("vus_recall").value), p(hl("blb_recall").value))
    mn, ma = mix.loc["new variants"], mix.loc["April set"]
    fact("The new variants are a different mix of ClinVar classes", p(mn.share_p_lp, 0), p(mn.share_vus, 0), p(mn.share_b_lb, 0), p(ma.share_p_lp, 0), p(ma.share_vus, 0), p(ma.share_b_lb, 0),
         p(mn.agreement_at_april_recall_of_each_direction), f"{(ma.agreement_in_direction - mn.agreement_at_april_recall_of_each_direction) * 100:.1f}",
         f"{(mn.agreement_at_april_recall_of_each_direction - mn.agreement_in_direction) * 100:.1f}", f"{(ma.agreement_in_direction - mn.agreement_in_direction) * 100:.1f}",
         f"{(mn.agreement_at_own_recall_and_april_mix - mn.agreement_in_direction) * 100:.1f}", f"{(ma.agreement_in_direction - mn.agreement_at_own_recall_and_april_mix) * 100:.1f}",
         "taking the mix first")
    fact("New Benign variants are much rarer than the old ones", f"{t('October, new', 'Benign').median_popmax_af * 100:.2f}%", f"{t(lab('April set (v37)'), 'Benign').median_popmax_af * 100:.1f}%",
         p(t("October, new", "Benign")["share_popmax_above_0.5pct"], 0), p(t(lab("April set (v37)"), "Benign")["share_popmax_above_0.5pct"], 0),
         p(t(lab("April set (v37)"), "Benign").recall), p(t("October, new", "Benign").recall), p(t(lab("April set (v37)"), "Likely_benign").recall), p(t("October, new", "Likely_benign").recall))
    fact("Recall of ClinVar Pathogenic variants is 64.2% against 65.0%", p(t("October, new", "Pathogenic").recall), p(t(lab("April set (v37)"), "Pathogenic").recall),
         p(t("October, new", "Likely_pathogenic").recall), p(t(lab("April set (v37)"), "Likely_pathogenic").recall))
    nr, ar = trate.loc["October, new"], trate.loc["April set"]
    r_apr, r_new, r_ratio, r_dir, r_cons, r_cls, r_clscons = (ratio.iloc[i] for i in range(7))
    fact("Among the new variants it is 0.017%", f"{nr.rate * 100:.3f}%", str(int(nr.opposite)), n(nr.records), f"{r_new.ci_lo * 100:.3f}", f"{r_new.ci_hi * 100:.3f}", f"{nr.rate / ar.rate:.1f}",
         f"{ar.rate * 100:.4f}%", f"{r_apr.ci_lo * 100:.3f}", f"{r_apr.ci_hi * 100:.3f}", f"{r_ratio.ci_lo:.1f}", f"{r_ratio.ci_hi:.1f}")
    fact("The mix of ClinVar class and consequence does not explain the excess", f"{r_cons.value:.1f}", f"{r_clscons.value:.1f}", f"{r_cls.value:.1f}",
         f"{r_dir.value:.1f}", str(int(r_cons.k)))
    bn = stdr.loc["Benign"]
    fact("At the April recall of each of nine bands of population-maximum allele frequency", p(bn.new_recall_at_april_recall_per_band), p(bn.share_of_the_difference_explained, 0),
         f"{(bn.april_recall - bn.new_recall_at_april_recall_per_band) * 100:.1f}", f"{(bn.april_recall - bn.new_recall) * 100:.1f}",
         p(stdr.loc["Likely_benign", "new_recall_at_april_recall_per_band"]))
    fact("We made the split into nine bands after a split into six", p(bn.six_band_share_of_the_difference_explained, 0), "nine bands", "six")
    ha, hn = hi_af.loc["April set"], hi_af.loc["new variants"]
    bb = bands[bands.clinvar == "Benign"].copy()
    bb["w"] = bb.n_new / bb.n_new.sum()
    bb["c"] = bb.w * (bb.recall_new - bb.recall_april) * 100
    top2 = -bb[bb.popmax_af_band.isin(["2e-2 to 5e-2", ">=5e-2"])].c.sum()
    # What the mix of bands leaves unexplained is the sum over bands of (April recall less new recall) times the band's share of the new variants.
    tot_c = -bb.c.sum()
    mid_c = -bb[bb.popmax_af_band.isin(["5e-3 to 1e-2", "1e-2 to 2e-2"])].c.sum()
    low_c = -bb[bb.popmax_af_band == "1e-5 to 1e-4"].c.sum()
    gap = ((bb.recall_april - bb.recall_new) * 100).set_axis(bb.popmax_af_band)
    fact("What remains is spread over the bands", f"{top2:.1f}", f"{tot_c:.1f}", p(top2 / tot_c, 0) if top2 / tot_c < 0.5 else "?", f"{mid_c:.1f}", f"{low_c:.1f}")
    fact("The two at 2% or more have the widest gaps between the April recall and the new recall", f"{gap['2e-2 to 5e-2']:.0f}", f"{gap['>=5e-2']:.0f}",
         "have the widest gaps" if set(gap.nlargest(2).index) == {"2e-2 to 5e-2", ">=5e-2"} else "?",
         n(hn["benign_variants_at_2e-2_or_more"]), n(hn.called_benign_side), p(hn.share_called_benign_side, 0),
         n(ha["benign_variants_at_2e-2_or_more"]), p(ha.share_called_benign_side, 0), n(hn.not_called_benign_side), n(hn.of_them_with_bp4_alone))
    bk = blockers[blockers.set == "new variants"].set_index("reason_bs1_did_not_fire").variants
    bka = blockers[blockers.set == "April set"].set_index("reason_bs1_did_not_fire").variants
    withheld = ["gnomAD low-complexity flag", "gnomAD FILTER", "gnomAD segmental-duplication flag", "gene on the homology list"]
    fact("Annotating those 274 again with the details of each criterion shows why", n(hn.not_called_benign_side))
    fact("For 226 (82%) BA1 and BS1 were withheld by design", n(bk[withheld].sum()), p(bk[withheld].sum() / hn.not_called_benign_side, 0), n(bk["gnomAD low-complexity flag"]),
         n(bk["gnomAD FILTER"]), n(bk["gnomAD segmental-duplication flag"]), n(bk["gene on the homology list"]))
    fact("For 20 (7%) the criteria were evaluated", n(bk["evaluated and not met (filtering allele frequency below the bar)"]),
         p(bk["evaluated and not met (filtering allele frequency below the bar)"] / hn.not_called_benign_side, 0))
    fact("Another 26 met BS1 and stay at VUS", n(bk["BS1 met, alone"] + bk["BS1 met, with other criteria"]), n(bk["BS1 met, alone"]), n(bk["BS1 met, with other criteria"]))
    fact("have an allele number below the floor", w(int(bk["allele number below the floor"]))[2] if int(bk["allele number below the floor"]) == 2 else "?")
    fact("The 226 are 30% of the 753 new Benign variants", n(bk[withheld].sum()), p(bk[withheld].sum() / hn["benign_variants_at_2e-2_or_more"], 0), n(hn["benign_variants_at_2e-2_or_more"]),
         p(bka[withheld].sum() / ha["benign_variants_at_2e-2_or_more"], 1), n(bka[withheld].sum()), n(ha["benign_variants_at_2e-2_or_more"]))
    pw, pc = prop.loc["BP7 as Walker 2023 describes it"], prop.loc["the completed ClinGen file"]
    po, pb = prop.loc["PVS1 licensed by Orphanet's loss-of-function label"], prop.loc["both inputs of the PVS1 gap"]
    pl, pa, p0 = prop.loc["BP1's veto without the variant's own record"], prop.loc["no allele-number floor"], prop.loc[lab("none (v37 on the release)")]
    pv1 = prop.loc["BP1's veto at a count of 1"]
    wkr = run("v37_hk_bp7_walker")
    fact("BP7 as the splicing recommendation describes it adds 5,744 calls in direction", n(pw.new_delta_same_direction), p(pw.new_delta_same_direction / pw.new_variants),
         p((wkr.same_direction - base.same_direction) / N_ALL), n(wkr.same_direction - base.same_direction), str(int(wkr.false_benign - base.false_benign)),
         "one reference-discordant benign call" if pw.new_delta_false_benign == 1 else "?")
    fact("The completed ClinGen file adds 52 calls in direction", n(pc.new_delta_same_direction), "no opposite-direction call" if pc.new_delta_false_benign + pc.new_delta_false_pathogenic == 0 else "?",
         n(po.new_delta_same_direction), str(int(po.new_delta_false_pathogenic)), n(pc.april_delta_same_direction), str(int(pc.april_delta_false_pathogenic)),
         n(po.april_delta_same_direction), str(int(po.april_delta_false_pathogenic)), n(pb.new_delta_same_direction), str(int(pb.new_delta_false_pathogenic)),
         n(pb.april_delta_same_direction), str(int(pb.april_delta_false_pathogenic)),
         n(po.new_delta_same_direction / po.new_delta_false_pathogenic), n(pb.new_delta_same_direction / pb.new_delta_false_pathogenic),
         n(po.april_delta_same_direction / po.april_delta_false_pathogenic), n(pb.april_delta_same_direction / pb.april_delta_false_pathogenic), "on 2 events")
    fact("BP1's veto without the variant's own record changes no call", "no call" if pl.new_calls_changed == 0 else "?", "exactly 3", n(new_idx["new ClinVar pathogenic or likely pathogenic missense variants in genes where BP1 can fire (pLI at least 0.9, missense z below 2)"]),
         n(new_idx["  of them in a gene with exactly 3 such records, so that leaving the own record out lifts the veto"]), "NUP85 and RPS29",
         "neither has BP4" if new_idx["    of those, with BP4 (BP1 with BP4 gives Likely benign, BP1 alone gives VUS)"] == 0 else "?")
    fact("The measurement build discounts one record whenever ClinVar's allele-level record", f"up to {n(new_idx[bp1_can] - new_idx[own_in])} of the {n(new_idx[bp1_can])}")
    fact("BP1's veto at a count of 1 gives 137 more calls in direction", n(pv1.new_delta_same_direction), n(-pv1.new_delta_vus_called_blb), n(-pv1.new_delta_blb_correct),
         "no opposite-direction call" if pv1.new_delta_false_benign + pv1.new_delta_false_pathogenic == 0 else "?", str(int(p0.new_false_benign)),
         n(run("v37_sw_bp1veto_1").same_direction - base.same_direction), str(int(base.false_benign - run("v37_sw_bp1veto_1").false_benign)))
    fact("Removing the allele-number floor changes no call", "no call" if pa.new_calls_changed == 0 else "?", n(p0.new_with_allele_number_below_2000), n(p0.new_below_2000_at_or_above_the_bs1_bar), f"at most {int(p0.new_below_2000_at_or_above_the_bs1_bar)} variants")
    fact("These are seven runs on one set of variants", "seven" if len(hold) == 7 else "?")
    dc = tdir
    fact("ClinVar changed the direction of 157 of the 665,057", str(int(dc.n.sum())), n(kept))
    vus_out = dc[(dc.april3 == "VUS")].n.sum()
    fact("The classifier called VUS for 119 of the 134", str(int(dc[(dc.april3 == "VUS") & (dc.fv3 == "VUS")].n.sum())), str(int(vus_out)),
         str(int(dc[dc.cv3 == dc.fv3].n.sum())), str(int(dc[dc.april3.map({"P/LP": "P/LP", "VUS": "VUS", "B/LB": "B/LB"}) == dc.fv3].n.sum())))
    into_vus = dc[dc.cv3 == "VUS"]
    fact("For the 22 moved into VUS it agreed with the new class for 17", str(int(into_vus.n.sum())), str(int(into_vus[into_vus.fv3 == "VUS"].n.sum())))
    fact("Of the 8,483 April variants where the classifier committed against a ClinVar VUS", n(cell("VUS", "P/LP") + cell("VUS", "B/LB")),
         str(int(dc[(dc.april3 == "VUS") & (dc.cv3 == dc.fv3) & (dc.fv3 != "VUS")].n.sum())),
         str(int(dc[(dc.april3 == "VUS") & (dc.cv3 == "B/LB") & (dc.fv3 == "P/LP")].n.sum())))


# --------------------------------------------------------------------------------- Section 4 and others


def discussion() -> None:
    t = traj
    m = mainline
    fact("Agreement in direction with ClinVar reached 77.6% at v13", p(t.loc[13].same_direction_rate), f"{m.loc[14:37].same_direction_rate.min() * 100:.2f}%", f"{m.loc[14:37].same_direction_rate.max() * 100:.2f}%",
         str(int(t.loc[13].opposite)), str(int(t.loc[37].opposite)), str(int((ledger.run.str.extract(r"v(\d+)")[0].astype(int) >= 13).sum())), str(len(m.loc[14:37])))
    fact("A project that watched only agreement would have concluded at v13",
         f"{baseline.loc['R1 + R2 + R3', 'opposite'] / baseline.loc[lab('v37 (the classifier)'), 'opposite']:.0f}")
    fact("The PM2 threshold moved pathogenic-side recall by 19 points", "19", f"{(run('v37_sw_pm2ad_0.0005').plp_recall - run('v37_sw_pm2ad_0').plp_recall) * 100:.1f}",
         f"{(run('v37_sw_table_both_sides').plp_recall - run('v37_ix_table_pm2ad_0').plp_recall) * 100:.1f}",
         f"{(runs37.loc['v37_sw_pm2ad_0.0005'].plp_recall - runs37.loc['v37_sw_pm2ad_0'].plp_recall) * 100:.1f}",
         f"{(runs37.loc['v37_sw_table_both_sides'].plp_recall - runs37.loc['v37_ix_table_pm2ad_0'].plp_recall) * 100:.1f} in v37")
    na = run("v37_src_noclinvar_all")
    fact("Removing every ClinVar-derived input costs 0.78 points", pp((na.same_direction - base.same_direction) / N_ALL, 2), str(int(na.false_benign - base.false_benign)),
         str(int(base.false_benign)), str(int(na.false_benign)))
    fact("Almost none of the discordance is contradiction", n(hl("opposite").value), n(N_ALL), p(hl("opposite").value / N_ALL, 4), n(hl("false_benign").value))
    fact("Almost all the rest is a defined ClinVar class that the classifier leaves as VUS", n(cell("P/LP", "VUS") + cell("B/LB", "VUS")))
    bvs = sigs[(sigs.clinvar == "B/LB") & (sigs.fastvep == "VUS")].set_index("signature")
    fact("Part is design: 10,694 of the 35,179", n(plp_pts["n"][5]), n(cell("P/LP", "VUS")), n(pvs1["without PVS1"]), n(cell("B/LB", "VUS")), n(bvs.loc["BP4", "n"] + bvs.loc["BS1", "n"]),
         n(bvs.loc["PM2_Supporting", "n"]))
    fact("at the wrong tier, and for the pathogenic-side half", n(tier_total()), n(fv("lp_at_9_points_clinvar_pathogenic")), n(c5.loc["Pathogenic", "Likely_pathogenic"]))
    fact("The classifier's calls are therefore best read with their error profile", n(hl("precision_plp_calls").denominator), p(hl("precision_plp_calls").value),
         p(cell("VUS", "P/LP") / hl("precision_plp_calls").denominator), f"{cell('B/LB', 'P/LP') / hl('precision_plp_calls').denominator * 100:.2f}%")
    bl = run("v37_hk_bp1_loo")
    fact("It also counts the variant's own record, and 3 of the 28 rest on it", str(int(bl.false_benign - base.false_benign)), str(int(run("v37_src_noclinvar_protein").false_benign - base.false_benign)),
         str(int(run("v37_hk_strict").false_benign)), str(int(run("v37_hk_strict").opposite)))
    fact("takes the opposite-direction calls from 51 to 79", str(int(base.opposite)), str(int(na.opposite)))
    fact("The ClinGen file that the benchmark read lacked 104 curations, among them CFTR, and the classifier degraded", n(missing[LACKED]))
    fact("That check flags 625 genes", n(pvs1["without PVS1, genes"]), n(closure.loc[rn("v37_gap_clingen_complete"), "genes_gaining_pvs1"]), n(closure.loc[rn("v37_gap_both"), "genes_gaining_pvs1"]),
         n(fv("pvs1_gap_genes_not_in_clingen_file")))
    fact("The design is the defect worth naming", n(fv("benchmark_genes_not_in_clingen_file")), n(fv("benchmark_genes")))
    fact("Nothing in an aggregate showed it", str(int(star_cell(4, "P/LP", "VUS").n)), str(int(star_cell(4, "P/LP", "VUS").row_total)))
    fact("Five of the 51 residual variants (10%) had left the 2-star+ set", p((aprl.october_status == "left the 2-star+ set").sum() / N_OPP, 0))
    pr_n_, bl_n_ = int(c5.loc["Pathogenic"].sum() + c5.loc["Likely_pathogenic"].sum()), int(c5.loc["Likely_benign"].sum() + c5.loc["Benign"].sum())
    fact("It also has no negative controls", p(pr_n_ / (pr_n_ + bl_n_), 0))
    from scipy import stats as _st
    lo2, hi2 = _st.chi2.ppf(0.025, 4) / 2, _st.chi2.ppf(0.975, 6) / 2
    fact("The exchange rates of Section 3.3 are counts of differences", f"{lo2:.1f}", f"{hi2:.1f}")
    cl = closure.loc[rn("v37_gap_clingen_complete")]
    wk = run("v37_hk_bp7_walker")
    pcc = prop.loc["the completed ClinGen file"]
    fact("First, complete the ClinGen input", n(cl.plp_vus_to_plp), n(pcc.april_delta_same_direction), str(int(cl.blb_to_plp)), n(pcc.new_delta_plp_correct), n(pcc.new_delta_same_direction),
         n(fv("benchmark_genes_not_in_clingen_file")))
    fact("Second, a source of established loss-of-function mechanism", n(prop.loc["PVS1 licensed by Orphanet's loss-of-function label"].new_delta_same_direction),
         str(int(prop.loc["PVS1 licensed by Orphanet's loss-of-function label"].new_delta_false_pathogenic)))
    fact("Third, BP7 as the splicing recommendation describes it", n(wk.same_direction - base.same_direction), str(int(wk.false_benign - base.false_benign)),
         n(prop.loc["BP7 as Walker 2023 describes it"].new_delta_same_direction))
    b1_1 = run("v37_sw_bp1veto_1")
    pv1_ = prop.loc["BP1's veto at a count of 1"]
    pv1_ = prop.loc["BP1's veto at a count of 1"]
    lv1_ = run("v37_hk_bp1_loo_veto1")
    fact("Fourth, BP1's veto leaving out the variant's own record", str(int(bl.false_benign - base.false_benign)), n(b1_1.same_direction - base.same_direction),
         n(base.same_benign - b1_1.same_benign), str(int(lv1_.false_benign)), str(int(base.false_benign)), n(pv1_.new_delta_same_direction),
         "avoids no reference-discordant benign call" if pv1_.new_delta_false_benign == 0 else "?",
         "none on the new variants" if prop.loc["BP1's veto without the variant's own record"].new_delta_false_benign == 0 else "?")
    fact("Sixth, review of the allele-number floor", n(run("v37_sw_no_an_floor").same_direction - base.same_direction))
    fact("Second, a source of established loss-of-function mechanism", n(closure.loc[rn("v37_gap_both"), "gap_variants"] - closure.loc[rn("v37_gap_both"), "gap_variants_gaining_pvs1"]))
    fact("Building a software classifier for the ACMG/AMP standards was largely", str(int(ledger["class"].value_counts()["G"])), str(int(ledger.found_by.value_counts()["expert_review"])))
    c1_, c3_ = closure.loc[rn("v37_gap_clingen_complete")], closure.loc[rn("v37_gap_both")]
    pr_bl2 = int(c5.loc["Pathogenic"].sum() + c5.loc["Likely_pathogenic"].sum() + c5.loc["Likely_benign"].sum() + c5.loc["Benign"].sum())
    fact("On the 673,660 ClinVar 2-star+ variants on which every threshold was tuned", p(hl("same_direction").value), n(hlc("opposite").value), f"{hl('opposite_rate').value * 100:.4f}%",
         f"{hl('opposite').value / pr_bl2 * 100:.3f}%", str(int(run("v37_hk_strict").opposite)), n(tiers.loc["Pathogenic", "n"]), n(tiers.loc["Benign", "n"]),
         f"{baseline.loc['R1 + R2 + R3', 'opposite'] / baseline.loc[lab('v37 (the classifier)'), 'opposite']:.0f}", f"{w(fv('strict_opposite_called_benign'))[1]} of the added calls Benign")
    nw = temporal[temporal.stratum == "new"].iloc[0]
    fact("On the 116,078 variants that were not in the April set, agreement in direction is 72.6%", n(nw.stratum_n), p(nw.same_direction / nw.stratum_n),
         f"{trate.loc['October, new'].rate / trate.loc['April set'].rate:.1f}", f"{ratio.iloc[2].ci_lo:.1f}", f"{ratio.iloc[2].ci_hi:.1f}")
    fact("The largest gap, in PVS1, is in part an input of unknown cause", str(int(missing[LACKED] - 1)), p(c1_.gap_variants_gaining_pvs1 / c1_.gap_variants, 0),
         p((c3_.gap_variants_gaining_pvs1 - c1_.gap_variants_gaining_pvs1) / c3_.gap_variants, 0))
    fact("ClinVar moved 157 of 665,057 variants across a direction", str(int(tdir.n.sum())))
    fact("Several lessons are of practical use", "")
    fact("The mode of inheritance that BS2 and PM2 read was the weakest input of v37", n(fv37("clingen_genes_in_file")), n(fv37("inheritance_read_recessive_only")), n(fv37("inheritance_read_dominant_only")),
         n(fv37("inheritance_read_both")), n(fv37("benchmark_genes_established_tags_ad_only")), n(fv37("benchmark_genes_established_tags_ad_only_read_dominant")),
         n(fv37("benchmark_genes_established_tags_ar_only")), n(fv37("benchmark_genes_established_tags_ar_only_read_recessive")), "established ClinGen tags")
    fact("The benchmark's classifier reads the gnomAD exomes only", f"{(base.same_direction - run('v37_sw_pm2_needs_record').same_direction) / N_ALL * 100:.2f}",
         f"{(base.plp_recall - run('v37_sw_pm2_needs_record').plp_recall) * 100:.1f}", p(base.fire_PM2_Supporting / N_ALL, 0))
    fact("The price of an error was not applied uniformly", n((t.loc[11].plp_recall - t.loc[10].plp_recall) * 93_812 / (t.loc[11].false_pathogenic - t.loc[10].false_pathogenic)),
         n((t.loc[13].blb_recall - t.loc[12].blb_recall) * 284_550 / (t.loc[13].false_benign - t.loc[12].false_benign)), n(steps[steps.knob == "BS2 homozygote bar"].reset_index(drop=True).iloc[1].calls_per_error),
         n(steps[steps.knob == "BS1 bar"].reset_index(drop=True).iloc[2].calls_per_error),
         f"on {int(t.loc[11].false_pathogenic - t.loc[10].false_pathogenic)} and {int(t.loc[13].false_benign - t.loc[12].false_benign)} errors",
         f"on {int(steps[steps.knob == 'BS2 homozygote bar'].reset_index(drop=True).iloc[1].added_errors)} and {int(steps[steps.knob == 'BS1 bar'].reset_index(drop=True).iloc[2].added_errors)} errors")


# ------------------------------------------------------------------- abstract, introduction, repository


def abstract_and_intro() -> None:
    import datetime

    t = traj
    nw = temporal[temporal.stratum == "new"].iloc[0]
    ovn = pd.read_csv(DATA / "databases_old_vs_new.csv")
    cs = pd.read_csv(DATA / "databases_causes.csv").set_index("source")
    a_set = "April set"
    one = lambda v, k, db: float(ovn[(ovn.version == v) & (ovn.set == a_set) & (ovn.measure == k)].iloc[0][db])  # noqa: E731
    share = lambda src: cs.loc[src, "variants_lost"] / cs.loc[src, "variants_with_a_value_in_the_complete"]  # noqa: E731
    # the abstract (every number of it is a function of the data)
    fact("The ACMG/AMP criteria are underspecified and fragmented.")
    fact("Thirty-eight versions were scored on 673,660 ClinVar variants", n(N_ALL), "Thirty-eight" if int(t.index.max()) == 38 else "?", "six review rounds" if len(rounds) == 3 else "?")
    fact("whose 116,078 new variants are the only out-of-sample ones", n(nw.stratum_n))
    fact("The input is a sites-only VCF, which allows 15 of 28 criteria", str(len(crit)), str(len(status)))
    fact("The first drafts used SpliceAI and PhyloP databases lacking 8% of their values", p(share("SpliceAI"), 0), p(share("PhyloP"), 0) if p(share("SpliceAI"), 0) == p(share("PhyloP"), 0) else "?")
    fact("On the April set v38 agrees with ClinVar's class for 64.1%", p(hl("exact_match").value), p(hl("same_direction").value), str(N_OPP), f"{hl('opposite_rate').value * 100:.4f}%")
    fact("On the new variants agreement in direction is 72.6%", p(nw.same_direction / nw.stratum_n), f"{trate.loc['October, new'].rate / trate.loc['April set'].rate:.1f}", f"{ratio.iloc[2].ci_lo:.1f}", f"{ratio.iloc[2].ci_hi:.1f}")
    fact("The complete databases added 1.5 points of agreement and three contradictions", f"{(one('v38', 'same_direction', 'complete') - one('v38', 'same_direction', 'first')) / N_ALL * 100:.1f}",
         w(int(one("v38", "opposite", "complete") - one("v38", "opposite", "first")))[1])
    tot_rev = res_by_group[["clinvar", "fastvep", "neither", "needs_evidence", "disorder_specific"]].sum()
    fact("The reviewer supported ClinVar on 29 of 48 contradictions seen", str(int(tot_rev.clinvar)), str(N_REV))
    fact("The agreement is bounded: ClinVar can be wrong, there was one reviewer, and a sites-only input allows 15 of 28 criteria", str(len(crit)), str(len(status)))
    # the introduction
    weeks = (datetime.date(2026, 10, 6) - datetime.date(2026, 4, 21)).days / 7
    fact("That speed allowed 37 numbered benchmark runs in 24 weeks", f"{weeks:.0f} weeks".replace(" weeks", ""), "24 weeks")


def contribution() -> None:
    """The sentences that state the contribution: the problem and the bounds of Section 1, the opening of Section 4 and the lines before the conclusions.

    They quote what other sentences of the paper quote, so each number is tied to the data as the sentence that it summarises is, and each
    claim that the data decide (which setting agrees more, how many calls a threshold moves) is a condition and not a typed word.
    """
    fact("Fourth, an input of variants alone, with no genotype, pedigree or phenotype of any individual, allows", str(len(crit)), str(len(status)))
    # the opening of the Discussion
    m1 = rep.loc["SVI points, both sides, Likely Benign from -1 (Tavtigian)"]
    tbl = run("v37_sw_table_both_sides")
    fact("Scoring both sides by points raises benign-side recall and adds hundreds of reference-discordant benign calls",
         "raises benign-side recall" if m1.blb_recall > base.blb_recall else "?",
         "hundreds of reference-discordant benign calls" if m1.false_benign - base.false_benign >= 100 else "?",
         "gives up pathogenic-side recall" if tbl.plp_recall < base.plp_recall else "?", w(int(base.false_pathogenic - tbl.false_pathogenic)))
    span = lambda rs: max(rs.same_direction.max(), base.same_direction) - min(rs.same_direction.min(), base.same_direction)  # noqa: E731
    cutoffs = [rep.loc[f"shipped, Pathogenic from {k} points"] for k in (9, 10, 11, 12)]
    spread = lambda col: max(getattr(c, col) for c in cutoffs) - min(getattr(c, col) for c in cutoffs)  # noqa: E731
    moves = [span(runs[(runs.knob == "pm2_ad_af_threshold") & (runs.group == "sweep")]), span(runs[(runs.knob == "pm2_ar_af_threshold") & (runs.group == "sweep")]),
             abs(run("v37_sw_no_an_floor").same_direction - base.same_direction)]
    fact("The cutoff for Likely pathogenic, BP7's conservation test and the BS1 and BS2 thresholds matter", n(N_ALL),
         "fewer than a thousand calls" if max(moves) < 1000 else "?",
         "the cutoff for Pathogenic moves exact agreement and not agreement in direction" if spread("same_direction") < 1000 <= spread("exact") else "?")
    fact("Fourth, check each input against an independent listing", f"{(run('v37_sw_clinvar_informed').plp_recall - base.plp_recall) * 100:.1f}",
         "CFTR" if pvs1["CFTR null variants in the benchmark (stop, frameshift, canonical splice)"] > 0 else "?")
    r38 = csv("v38_runs.csv").set_index("run")
    fact("Counting one gnomAD observation as two benign criteria adds", n(r38.loc["v38_alt_unguarded"].d_same_direction - r38.loc["v38"].d_same_direction),
         "two ClinVar-Pathogenic variants Benign" if r38.loc["v38_alt_unguarded"].false_benign - r38.loc["v38"].false_benign == 2 else "?",
         p(baseline.loc["R1 + R2 + R3", "same_direction_share"]), n(baseline.loc["R1 + R2 + R3", "opposite"]), str(N_OPP))
    nw = temporal[temporal.stratum == "new"].iloc[0]
    fact("With these choices v38 agrees with ClinVar's class for", p(hl("exact_match").value, 0), p(hl("same_direction").value, 0), n(N_ALL), str(N_OPP),
         p(hl("opposite").value / N_ALL, 4))
    fact("On the 116,078 variants that were not in that set, agreement in direction is 72.6%", n(nw.stratum_n), p(nw.same_direction / nw.stratum_n),
         f"{trate.loc['October, new'].rate / trate.loc['April set'].rate:.1f}", f"{ratio.iloc[2].ci_lo:.1f}", f"{ratio.iloc[2].ci_hi:.1f}")
    fact("agreement in direction is higher with the PM2 threshold at strict absence than at the default",
         "agreement in direction is higher with the PM2 threshold at strict absence than at the default" if run("v37_sw_pm2ad_0").same_direction > base.same_direction else "?",
         f"avoid {w(int(run('v37_sw_phylop_1000').false_benign - base.false_benign))[1]} contradictions")
    fact("A sites-only input evaluates 15 of the 28 criteria, which caps what any annotator", str(len(crit)), str(len(status)))


def repository() -> None:
    import importlib
    import subprocess

    # The crate as it was at the commit that made the runs (the tag of release 0.5.0 for v38, 5d8004a for v37), not as it is now: comments added
    # since, when a citation was corrected, would otherwise change the line count.
    root = HERE.parents[1]

    def git(*args: str) -> str:
        return subprocess.run(["git", *args], cwd=root, capture_output=True, text=True, check=True).stdout

    lacking = [ref for ref in ("v0.5.0-rc.1", "5d8004a", "5293cc0") if subprocess.run(["git", "cat-file", "-e", f"{ref}^{{commit}}"], cwd=root, capture_output=True).returncode != 0]
    if lacking:
        failures.append(f"the sentences about the crate, the tag and the commits need {', '.join(lacking)}, which this clone lacks (git fetch --tags --unshallow)")
        return

    def crate(ref: str) -> tuple[int, int, str]:
        bodies = [git("show", f"{ref}:{f}") for f in git("ls-tree", "-r", "--name-only", ref, "crates/fastvep-classification/src").split() if f.endswith(".rs")]
        return sum(b.count("\n") for b in bodies), sum(b.count("#[test]") for b in bodies), re.search(r'^version = "([^"]+)"', git("show", f"{ref}:Cargo.toml"), re.M).group(1)

    (l38, t38, ver38), (l37, t37, ver37) = crate("v0.5.0-rc.1"), crate("5d8004a")
    is_ancestor = lambda a, b: subprocess.run(["git", "merge-base", "--is-ancestor", a, b], cwd=root).returncode == 0  # noqa: E731
    fact("The classification engine is the `fastvep-classification` crate", n(l38), n(t38), ver38)
    # the paper calls the candidate "release 0.5.0" for as long as no other tag of 0.5.0 exists, so the text is stale the day that v0.5.0 is tagged
    only_rc = git("tag", "--list", "v0.5*").split() == ["v0.5.0-rc.1"]
    fact("The tag is `v0.5.0-rc.1` (c12d9f9), the first release candidate of fastVEP release 0.5.0", git("rev-parse", "--short=7", "v0.5.0-rc.1^{commit}").strip(),
         "the only tag of it so far" if only_rc else "?",
         "5293cc0" if is_ancestor("5293cc0", "v0.5.0-rc.1") and git("log", "-1", "--format=%s", "5293cc0").startswith("fix(acmg): read the mode of inheritance") else "?")
    fact("v37 (5d8004a) is the version in which the review rounds ended", git("rev-parse", "--short=7", "5d8004a").strip(), n(l37), n(t37), ver37)
    fact("v38 is the classifier of fastVEP release 0.5.0, which for now is the release candidate", "`v0.5.0-rc.1`" if only_rc else "?")
    fact("Every number is for v38, the latest of these versions of the classifier", "v0.5.0-rc.1" if only_rc else "?", "GRCh38")
    v = {m: importlib.import_module(m).__version__ for m in ("pandas", "numpy", "scipy", "matplotlib")}
    fact("Analysis used Python", f"pandas {v['pandas']}", f"NumPy {v['numpy']}", f"SciPy {v['scipy']}", f"matplotlib {v['matplotlib']}", f"Python {sys.version.split()[0]}")
    try:
        rustc = subprocess.run(["rustc", "--version"], capture_output=True, text=True, check=True).stdout.split()[1]
        fact("Analysis used Python", f"Rust {rustc}")
    except (OSError, subprocess.CalledProcessError):
        pass


def availability() -> None:
    """Every file and directory that "Data and code availability" names is in the repository.

    The section tells a reader where each script, patch, table and record of the paper is. A name that the repository does not hold,
    because a script was renamed or a file was never added, is the kind of statement that no number check catches. A bare name is
    looked for by its file name anywhere in the repository, a path by its place, and a directory by a file under it.
    """
    import subprocess

    root = HERE.parents[1]
    # Tracked files and the files that are not yet but are not ignored either, so that a file added with its text passes before it is committed
    # and a file under an ignored path (such as the bulk outputs under data/) does not.
    tracked = subprocess.run(["git", "ls-files", "--cached", "--others", "--exclude-standard"], cwd=root, capture_output=True, text=True, check=True).stdout.split("\n")
    names = {Path(f).name for f in tracked if f}
    section = (HERE / "manuscript.md").read_text().split("\n## Data and code availability\n", 1)[1].split("\n## ", 1)[0]
    found = 0
    for token in dict.fromkeys(re.findall(r"`([^`\s]+)`", section)):
        if not re.search(r"\.(py|sh|patch|csv|tsv|xlsx|md|toml)$|/$", token):
            continue
        if token.endswith("/"):
            ok = any(f.startswith(token) or f.startswith("manuscript/acmg/" + token) for f in tracked) or (root / token).is_dir()
        elif "/" in token:
            ok = token in tracked or ("manuscript/acmg/" + token) in tracked
        else:
            ok = token in names
        found += ok
        if not ok:
            failures.append(f"Data and code availability names {token!r}, which the repository does not hold")
    passed.append(f"Data and code availability: {found} names found in the repository")


def audit() -> None:
    fact("Regrading both by one step would change 2 of 673,660 calls", n(fv("pvs1_regrade_calls_changed")), n(fv("pvs1_regrade_vus_to_lp")), n(fv("pvs1_regrade_b_to_vus")), n(N_ALL))
    fact("Third, BP7 keeps the 2015 test of poor conservation", "0.1", "2.0")
    nm = run("v37_sw_pvs1_nmd50")
    fact("Second, PVS1 predicts NMD by whether the variant lies in the last exon", n(base.same_direction - nm.same_direction), f"{(base.plp_recall - nm.plp_recall) * 100:.2f}",
         w(int(base.false_pathogenic - nm.false_pathogenic)))


def early_trades() -> None:
    t = traj
    n_plp, n_blb = 93_812, 284_550
    fact("Opposite-direction calls had reached a minimum of 25 at v10", n((t.loc[11].plp_recall - t.loc[10].plp_recall) * n_plp),
         str(int(t.loc[11].false_pathogenic - t.loc[10].false_pathogenic)))
    fact("They rose again at v13, when BP7's intronic extension", n((t.loc[13].blb_recall - t.loc[12].blb_recall) * n_blb), str(int(t.loc[13].false_benign - t.loc[12].false_benign)))
    fact("Building BP3's database moved same-direction agreement by 24 calls", n((t.loc[15].same_direction_rate - t.loc[14].same_direction_rate) * N_ALL),
         f"{(t.loc[15].same_direction_rate - t.loc[14].same_direction_rate) * 100:.3f}", "no opposite-direction call" if t.loc[15].opposite == t.loc[14].opposite else "?")
    fact("Switching on PS1's splice path left same-direction agreement unchanged", n((t.loc[16].exact_rate - t.loc[15].exact_rate) * N_ALL),
         f"{(t.loc[16].exact_rate - t.loc[15].exact_rate) * 100:.2f}")
    fact("The workbook of round 1 listed 101 opposite-direction calls, but v7 had 313", str(int(t.loc[7].opposite)))
    fact("so 212 of v7's 313 had never been listed", str(int(t.loc[7].opposite - 101)))


def criteria_consistency() -> None:
    """The per-criterion table, the run table and the facts table count the same firings."""
    cr = csv("criteria_firing_v37.csv").set_index("criterion")
    checks = [
        ("BP7 firings", int(cr.loc["BP7", "variants"]), int(fv("bp7_firings"))),
        ("BP7 firings (run table)", int(cr.loc["BP7", "variants"]), int(base.fire_BP7)),
        ("PVS1 at full strength", int(cr.loc["PVS1", "strengths"].split("default ")[1].split(";")[0].replace(",", "")), int(base.fire_PVS1)),
        ("BA1 firings", int(cr.loc["BA1", "variants"]), int(base.fire_BA1)),
        ("BS1 firings", int(cr.loc["BS1", "variants"]), int(base.fire_BS1)),
        ("PM2 firings", int(cr.loc["PM2", "variants"]), int(base.fire_PM2_Supporting)),
    ]
    for name, a, b in checks:
        if a == b:
            passed.append(f"consistency: {name} = {a}")
        else:
            failures.append(f"{name}: the criteria table says {a} and the other table {b}")


def coverage() -> None:
    """Sentences whose numbers are functions of the data and sit in a line that no other check anchors."""
    t = traj
    f = ledger.found_by.value_counts()
    new = temporal[temporal.stratum == "new"].iloc[0]
    fact("we applied it unchanged to ClinVar's release of 4 October 2026, extracted the same way", n(temporal[temporal.stratum == "all"].iloc[0].stratum_n))
    fact("The output does not record the detail the veto needs", str(int(veto.held_at_vus_by_the_veto)))
    r6row = rounds.loc[6]
    fact("In round 6 the workbook asked for the reviewer's class", str(int(r6row.rows)), str(int(r6row.rows_with_written_reason)))
    w_late = int(((ledger["class"] == "W") & (ledger.run.str.extract(r"v(\d+)")[0].astype(int) > 9)).sum())
    fact("Eight of the 17 wiring defects were found after v9", w(w_late), str(int((ledger["class"] == "W").sum())))
    fact("**Review of the discordant calls first identified 30 of the 64 changes", str(int(f["expert_review"])), str(len(ledger)))
    fact("Seventh, the ClinGen file that the benchmark read lacked", n(missing[LACKED]))
    fact("**Counting a variant's own ClinVar record would add 14.3 points", f"{(run('v37_sw_clinvar_informed').plp_recall - base.plp_recall) * 100:.1f}",
         pp((run("v37_src_noclinvar_all").same_direction - base.same_direction) / N_ALL, 2))
    fact("The tiers of v38 are informative in the same way, with their denominators", n(tiers.loc["Pathogenic", "n"]), f"{tiers.loc['Pathogenic', 'opposite_ci_hi'] * 100:.3f}%",
         n(tiers.loc["Benign", "n"]), f"{tiers.loc['Benign', 'opposite_ci_hi'] * 100:.4f}%", str(int(tiers.loc["Likely_pathogenic", "clinvar_opposite"])),
         n(tiers.loc["Likely_pathogenic", "n"]), f"{tiers.loc['Likely_pathogenic', 'opposite_rate'] * 100:.3f}%", str(int(tiers.loc["Likely_benign", "clinvar_opposite"])),
         n(tiers.loc["Likely_benign", "n"]), f"{tiers.loc['Likely_benign', 'opposite_rate'] * 100:.3f}%",
         "none" if tiers.loc["Pathogenic", "clinvar_opposite"] == 0 and tiers.loc["Benign", "clinvar_opposite"] == 0 else "?")
    fact("The zero for Benign is informative", f"{tiers.loc['Benign', 'n'] * tiers.loc['Likely_benign', 'opposite_rate']:.0f}",
         f"{tiers.loc['Pathogenic', 'n'] * tiers.loc['Likely_pathogenic', 'opposite_rate']:.1f}", "two of the 15 curated alleles (both G6PD)", w(fv("strict_opposite_called_benign"))[1])
    dates = lambda k: missing[k]
    fact("Their dates are recent: 83 of the 104", n(dates("  classified in 2025 or 2026")), n(dates("  classified from January to April 2026")),
         p(dates("  classified in 2025 or 2026") / dates(LACKED), 0), p(dates("  classified from January to April 2026") / dates(LACKED), 0),
         n(dates("curations in the April export classified in 2025 or 2026")), n(dates("curations in the 30 April export")),
         p(dates("curations in the April export classified in 2025 or 2026") / dates("curations in the 30 April export"), 0), n(dates("curations in the April export classified from January to April 2026")),
         p(dates("curations in the April export classified from January to April 2026") / dates("curations in the 30 April export"), 1),
         f"about {dates('  expected in 2025 or 2026 if the omission were independent of the date'):.0f} and {dates('  expected from January to April 2026 if the omission were independent of the date'):.0f}",
         w(dates("  classified before 2024"))[1])
    fact("**Where the 35,179 pathogenic-side variants called VUS come from", n(cell("P/LP", "VUS")))
    fact("**Where the 95,230 benign-side variants called VUS come from", n(cell("B/LB", "VUS")))
    fact("**The 3,610 and 4,873 variants where v38 commits", n(cell("VUS", "P/LP")), n(cell("VUS", "B/LB")))
    pvs1_group = res_by_group.loc["PVS1 on a ClinVar benign-side variant"]
    bp47_group = res_by_group.loc["BP4+BP7 on a ClinVar pathogenic-side variant"]
    fact("**PVS1 on a variant ClinVar calls benign or likely benign (21)", str(int(pvs1_group.total)))
    fact("**BP4 and BP7 on a variant ClinVar calls pathogenic or likely pathogenic (22)", str(int(bp47_group.total)))
    fact("These are the 22 reference-discordant benign calls of Section 3.4 that rest on BP7", str(int(bp47_group.total)), str(int(base.false_benign - run("v37_src_nophylop").false_benign)), "VWF and GTPBP3")
    fact("The reviewer saw 20 of the 22", str(int(bp47_group.total - bp47_group["not reviewed"])), str(int(bp47_group.total)))
    fact("The reviewer supported ClinVar for 10", str(int(bp47_group.clinvar)))
    opp = csv("opposite_direction_v37.csv")
    pol = opp[opp.gene == "POLR3A"].iloc[0]
    fact("POLR3A c.1909+22G>A, at an allele frequency of 2.5e-3", f"{pol.gnomad_af * 1e3:.1f}e-3")
    pad, slc = opp[opp.gene == "PADI3"].iloc[0], opp[opp.gene == "SLC12A3"].iloc[0]
    fact("The BS1 pair (PADI3 and SLC12A3)", n(pad.gnomad_ac), n(slc.gnomad_ac), p(pad.gnomad_popmax_af), p(slc.gnomad_popmax_af))
    tm216 = opp[opp.gene == "TMEM216"].iloc[0]
    fact("A third BS1 call, TMEM216 11:61392563 G>T (upstream of the gene)", f"{tm216.chrom}:{tm216.pos} {tm216.ref}>{tm216.alt}", "BP4" if "BP4" in tm216.criteria and "BS1" in tm216.criteria else "?",
         "0.07" if "0.07" in str(csv("databases_opposite_changes.csv").query("gene == 'TMEM216'").complete_spliceai_max.iloc[0]) else "?", "the reviewer has not seen it" if tm216.ruling == "not reviewed" else "?")
    nopp = csv("temporal_new_opposite_v37.csv")
    fact("the opposite-direction rates and their ratio with its interval, the 20 opposite-direction calls on the new variants", n(len(nopp)))
    kinds = nopp.kind.value_counts()
    fact("The 18 are of three kinds that we had named on the April residual", str(int(kinds["PVS1 on a ClinVar benign-side variant"])),
         str(int(kinds["BP4+BP7 on a ClinVar pathogenic-side variant"])),
         str(int(kinds["missense pathogenic-side call on PP2 and PP3 for a ClinVar benign-side variant"])), *sorted(set(nopp.gene)))
    still = aprl[aprl.october_status == "still 2-star+"]
    fact("The 46 residual calls that are still in the set are unchanged in class", str(len(still)), "unchanged" if (still.october_stratum == "unchanged").all() else "?")
    ng = csv("temporal_new_genes_v37.csv").iloc[0]
    fact("The new variants lie in 4,862 genes", n(ng.genes), n(ng.genes_not_in_the_april_set), p(ng.genes_not_in_the_april_set / ng.genes))
    fact("None of the 20 opposite-direction calls on the new variants rests on BP1, PS1 or PM5", n(ng.opposite_calls),
         "None" if ng.opposite_calls_with_bp1 + ng.opposite_calls_with_ps1 + ng.opposite_calls_with_pm5 == 0 else "?",
         "one rests on PM1" if ng.opposite_calls_with_pm1 == 1 else "?", "at most one" if ng.opposite_calls_with_pm1 == 1 else "?")
    fact("The two largest steps before v13 (8.8 and 5.2 points", f"{(t.loc[4].same_direction_rate - t.loc[2].same_direction_rate) * 100:.1f}",
         f"{(t.loc[6].same_direction_rate - t.loc[5].same_direction_rate) * 100:.1f}")
    fact("Review first identified 30 of the 64 ledger items (47%, by our coding)", str(int(f["expert_review"])), str(len(ledger)), f"{f['expert_review'] / len(ledger) * 100:.0f}%",
         str(int(ledger["class"].value_counts()["G"])))
    fact("The reviewer supported ClinVar's call on 89%, 86% and 67%", p(rounds.loc[1].clinvar_share_of_verdicts, 0), p(rounds.loc[2].clinvar_share_of_verdicts, 0),
         p(rounds.loc[6].clinvar_share_of_verdicts, 0))
    fact("A benchmark cannot find a missing condition when the omission does not change an aggregate", n((t.loc[15].same_direction_rate - t.loc[14].same_direction_rate) * N_ALL))
    c5r = run("v37_cut_lp5")
    pm2 = runs[(runs.knob == "pm2_ad_af_threshold") & (runs.group == "sweep")]
    fact("The cutoffs that define the classes move the pathogenic side more than any numeric threshold", f"{(c5r.plp_recall - base.plp_recall) * 100:.0f}",
         f"{(max(pm2.plp_recall.max(), base.plp_recall) - min(pm2.plp_recall.min(), base.plp_recall)) * 100:.1f}")
    fact("Counting a variant's own ClinVar record would add 14.3 points of pathogenic-side recall, and the classifier excludes it", f"{(run('v37_sw_clinvar_informed').plp_recall - base.plp_recall) * 100:.1f}")
    fact("A ClinVar agreement of 64% exact and 79% in direction", p(hl("exact_match").value, 0), p(hl("same_direction").value, 0))
    tot = res_by_group[["clinvar", "fastvep"]].sum()
    fact("In the last round the reviewer supported ClinVar on 60% of the 48", p(tot.clinvar / N_REV, 0), p(tot.fastvep / N_REV, 0))
    fact("saw only the opposite-direction calls and none of the 35,179", n(cell("P/LP", "VUS")))
    fact("is read from the mode that ClinGen curates for each of the", n(fv("clingen_genes_in_file")))
    fact("The classifier evaluates 15 of 28 criteria", str(len(crit)))
    fact("The ledger of 64 items and its classes are our own coding", str(len(ledger)))
    fact("Each of these can be measured by the harness of this paper before it is adopted", str(int(hl("opposite").value)))
    # Sentences that quote the count of residual opposite-direction calls, a reviewer's ruling or a figure of the main run in a line that no check above anchors
    fact("The reviewer's remarks on the residual 48 opposite-direction calls", str(N_REV))
    fact("Both are among the 51 and were found in the set on which the count was swept", str(N_OPP),
         f"Both are among the {N_OPP}" if {"IQSEC1", "PRR12"} <= set(r6[r6.still_opposite.astype(str).eq("True")].gene) else "?")
    fact("**Counting a variant's own ClinVar record would add 14.3 points of pathogenic-side recall, and removing every ClinVar-derived input costs 0.78 points of agreement",
         str(N_OPP), str(int(run("v37_src_noclinvar_all").opposite)))
    fact("Letting it count, as the runs before v9 did", str(int(base.opposite - run("v37_sw_clinvar_informed").opposite)), "of the " + str(N_OPP))
    fact("Variants called pathogenic-side are almost all null variants", n(base.same_pathogenic))
    fact("Those 5,184 are 57% of the 9,136", n(pvs1["without PVS1"]))
    fact("For 4,696 of the 5,184 variants", n(pvs1["without PVS1"]))
    fact("CFTR holds 439 of the 5,184 variants", n(pvs1["without PVS1"]))
    fact("Even so 3,470 of the 5,184 variants", n(closure.loc[rn("v37_gap_both")].gap_variants))
    fact("### 3.6 The 51 opposite-direction calls that remain", str(N_OPP))
    fact("The last review round (v26) held 64 opposite-direction calls, and 48 of the 51 in v38 are among them", str(N_REV), str(N_OPP), w(N_OPP - N_REV)[1] if N_OPP - N_REV == 3 else "?", "briefly had " + str(int(traj.loc[31].opposite)))
    fact("Each of the 48 had already been in the table of an earlier round", str(N_REV))
    fact("**Table 6. The 51 opposite-direction calls of v38", str(N_OPP))
    fact("Supplementary Table S8 lists the 51 calls with their criteria", str(N_OPP))
    fact("For 29 of the 48 (60%) the reviewer supported ClinVar's call", str(N_REV), str(int(fv("reviewer_supported_clinvar_false_benign"))), str(REV_FB),
         str(int(fv("reviewer_supported_clinvar_false_pathogenic"))), str(REV_FP))
    fact("for the other 19 the reviewer did not side with ClinVar", str(N_REV - int(tot.clinvar)))
    fact("Five of the 51 (10%) are no longer in the 2-star+ set", str(N_OPP))
    fact("The 51 fall into five groups", str(N_OPP))
    fact("A further 8,603 of the April variants have left the set", str(N_OPP), w(int((aprl.october_status == "left the 2-star+ set").sum()))[1])
    newr = temporal[temporal.stratum == "new"].iloc[0]
    fact("**On variants that were never in the benchmark, agreement in direction is", p(newr.same_direction / newr.stratum_n), p(hl("same_direction").value), f"{trate.loc['October, new'].rate / trate.loc['April set'].rate:.1f}")
    fact("Across the 772,585 covered records, agreement in direction is", p(hl("same_direction").value))
    n_blb = int(c5.loc[["Likely_benign", "Benign"]].sum().sum())
    n_blb_called = n_blb - int(c5.loc[["Likely_benign", "Benign"], "NoCall"].sum())
    fact("The recall of ClinVar's benign-side variants is over the 284,168 that have a call", n(n_blb_called), n(n_blb), n(c5.loc[["Likely_benign", "Benign"], "NoCall"].sum()),
         f"{(b38.same_benign / n_blb_called - b38.same_benign / n_blb) * 100:.1f} points higher")
    fact("The reading is as good as the field it reads", n(missing[LACKED]))
    f38v = pd.read_csv(DATA / "v38_facts.csv").set_index("name")["value"]
    fact("(it fired on 7,234 of the 7,328 variants that meet BS1", n(f38v["ad_only_variants_with_bs1_and_bs2"]), n(f38v["ad_only_variants_with_bs1"]))
    fact("The same 96 runs were made for v37", n(len(harness) + len(hold)), "The same 96 runs" if same_96() else "?")
    fact("Agreement in direction with ClinVar reached 77.6% at v13 and stayed between", p(t.loc[38].same_direction_rate, 2), f"{(t.loc[38].same_direction_rate - t.loc[37].same_direction_rate) * 100:.2f} points",
         str(int(t.loc[38].opposite)))
    gix37 = runs37.loc["v37_sw_table_both_sides"].same_pathogenic - runs37.loc["v37_ix_table_pm2ad_0"].same_pathogenic
    fact("The PM2 threshold moved pathogenic-side recall by 19 points when it was set in v10", f"{gix37 / 93_812 * 100:.1f}",
         f"{(run('v37_sw_table_both_sides').same_pathogenic - run('v37_ix_table_pm2ad_0').same_pathogenic) / 93_812 * 100:.1f}")
    fact("Removing every ClinVar-derived input costs 0.78 points of agreement", str(int(run("v37_sw_no_bp1_veto").false_benign - base.false_benign)))
    fact("It also counts the variant's own record, and 3 of the 28 rest on it", str(N_OPP), str(int(run("v37_hk_bp1_loo").false_benign - base.false_benign)), str(int(run("v37_sw_no_bp1_veto").false_benign - base.false_benign)),
         str(int(run("v37_hk_strict").false_benign)), str(int(run("v37_hk_strict").opposite)))
    fact("The classifier's calls are therefore best read with their error profile beside them", n(tiers.loc["Pathogenic", "n"]), n(tiers.loc["Benign", "n"]))
    fact("In the last round the reviewer supported ClinVar on 60% of the 48", str(N_REV))
    fact("Five of the 51 residual variants (10%) had left the 2-star+ set by October", str(N_OPP))
    fact("A discordance with ClinVar is a prompt to look", str(N_REV))
    fact("Count severe errors as a primary measure beside agreement", str(int(t.loc[10].opposite)), str(N_OPP))
    fact("All thresholds were tuned on the April benchmark and its discordances", str(N_OPP), str(int(run("v37_hk_strict").opposite)))


def v38() -> None:
    """Section 2.10 and Section 3.8 (the reading of ClinGen's curated inheritance), and the sentences elsewhere that quote them."""
    res = HERE.parents[1] / "analysis/acmg_benchmark/results/v37_manuscript"
    runs38 = pd.read_csv(res / "runs_v38.tsv", sep="\t")
    hold38 = pd.read_csv(res / "holdout_runs_v38.tsv", sep="\t")
    # The files hold the family v38's runs beside Section 3.8's: its April runs are those of its own groups, and its October
    # runs are the four that the family v37 has no twin of.
    runs38 = runs38[runs38.group.isin(["inherit", "inherit_base", "inherit_alt", "joint"])]
    hold38 = hold38[~hold38.name.map(lambda x: "v37" + x[3:]).isin(set(hold37.name))]
    r = csv("v38_runs.csv").set_index("run")
    ctl = csv("v38_controls.csv").iloc[0]
    fire = csv("v38_criteria_firing.csv").set_index("criterion")
    tr = csv("v38_transitions.csv")
    mv = csv("v38_clinvar_plp_moves.csv")
    f38 = csv("v38_facts.csv").set_index("name")["value"]
    rd = csv("v38_inheritance_reading.csv")
    pdm = csv("v38_per_disorder.csv")
    pdg = csv("v38_per_disorder_genes.csv").set_index("symbol")
    sw = csv("v38_sweeps.csv").set_index("name")
    st38 = csv("v38_threshold_steps.csv")
    tp = csv("v38_triple.csv").set_index("set")
    ho = csv("v38_holdout.csv")
    hmv = csv("v38_holdout_plp_moves.csv")
    jn = csv("v38_joint.csv")
    jc = csv("v38_joint_changes.csv")
    shipped, before, unguarded = r.loc["v38"], r.loc["v38_base"], r.loc["v38_alt_unguarded"]
    label, as_v37, triple = r.loc["v38_alt_label_guard"], r.loc["v38_alt_bs2_as_v37"], r.loc["v38_alt_triple_cap"]
    i = lambda x: int(x)

    def reading(scope: str, kind: str, **match):
        x = rd[(rd.scope == scope) & (rd.kind == kind)]
        for k, v in match.items():
            x = x[x[k] == v]
        assert len(x) == 1, (scope, kind, match, len(x))
        return x.iloc[0]

    def pd_row(measure: str, nth: int = 0):
        return int(pdm[pdm.measure == measure].iloc[nth].n)

    def the_ho(set_: str, name: str):
        return ho[(ho.set == set_) & (ho.run == name)].iloc[0]

    allf, ben = "every gene of the file", "genes with a benchmark variant"
    names_all = reading(allf, "reading", reading="disease names (to v37)")
    field_all = reading(allf, "reading", reading="curated field (v38)")
    names_b = reading(ben, "reading", reading="disease names (to v37)")
    field_b = reading(ben, "reading", reading="curated field (v38)")
    ad_b = reading(ben, "group", group="established tags AD only")
    ar_b = reading(ben, "group", group="established tags AR only")
    both_b = reading(ben, "group", group="established tags AD and AR")
    other = lambda name: i(reading(ben, "other", group=name).genes)

    # ----------------------------------------------------------------------------------- Section 2.10
    fixture = (HERE.parents[1] / "crates/fastvep-cli/fixtures/clingen/genemap.tsv").read_text().splitlines()
    rows = [line for line in fixture if line and not line.startswith("#")]
    fact("A test builds the gene database with the command-line builder from 263 real curations in 185 genes", n(sum(line.count("(ClinGen ") for line in rows)), n(len(rows)),
         "2026" if "11 October 2026" in LINES[[k for k, line in enumerate(LINES) if line.startswith("A test builds the gene database")][0]] else "?")
    fact("The control of that build, with nothing set, reproduces v38 in all 30 cells", n(ctl.cells_equal_to_v38), n(ctl.cells), n(ctl.criteria_firing_equal_to_v38), n(ctl.criteria),
         "all" if ctl.cells_equal_to_v38 == ctl.cells and ctl.criteria_firing_equal_to_v38 == ctl.criteria else "?")
    fact("were swept again on v38 in 26 runs", n(len(sw[sw.group == "sweep"])))
    fact("were applied unchanged to the release of 4 October 2026 (Section 2.8) in 4 runs", n(len(hold38)))
    every_chromosome = {"v38", "v38_noflags", "v38_joint"}
    fact("Three runs on every chromosome compare", w(len(every_chromosome & set(runs38.name))),
         w(sum(runs38.name.str.startswith("v38_chr2122"))), "chromosomes 21 and 22")
    fact("In all there are 38 runs on the April set and 4 on the October release", n(len(runs38)), n(len(hold38)))
    again42 = pd.read_csv(HERE.parents[1] / "analysis/acmg_benchmark/results/v38_v0.5.0/rerun_comparison.csv")
    again42 = again42[again42.kind == "run"]
    fact("All 42 were made on 8 November 2026 with the complete databases", n(len(again42)),
         "8 November 2026" if {new_york_date(rerun_recorded_at(x)) for x in again42.name} == {"2026-11-08"} else "?",
         "agree in every scored file" if set(again42.name) == set(runs38.name) | set(hold38.name) and (again42.identical == "yes").all() else "?")
    shared = [x for x in runs38.name if x.startswith(("v38_sw_pm2ad_", "v38_sw_pm2ar_", "v38_sw_bs2hom_"))]
    fact("The runs that Section 3.8 made on its own, 38 on the April set and 4 on the October release", n(len(runs38)), n(len(hold38)), n(len(shared)),
         n(len(harness) - len(shared)), n(len(hold)))

    # ----------------------------------------------------------------------------------- Section 3.8, v38 against v37 (the same under both bases)
    fam = pd.read_csv(HERE.parents[1] / "analysis/acmg_benchmark/results/v38_family_vs_v37_summary.csv").set_index("name")["value"]
    f38v = pd.read_csv(DATA / "v38_facts.csv").set_index("name")["value"]
    fact("Every run of Sections 2.6 to 2.8 was made with both versions", n(len(harness) + len(hold)))
    d_dir = b38.same_direction - base37.same_direction
    fact("On the benchmark v38 gives 11,449 more calls in direction than v37", n(d_dir), "one opposite-direction call fewer" if base37.opposite - b38.opposite == 1 else "?")
    fact("v38 agrees with ClinVar's class for 431,727 of the 673,660 records", n(b38.exact), n(N_ALL), n(b38.same_direction), p(b38.same_direction / N_ALL, 2), n(base37.exact),
         n(base37.same_direction), p(base37.same_direction / N_ALL, 2), n(b38.exact - base37.exact), n(d_dir), pp(d_dir / N_ALL, 2))
    fact("Of the 1,494, 4 are the commit before the change and 1,490 are the change itself", n(d_dir), n(before.same_direction - base37.same_direction),
         n(shipped.same_direction - before.same_direction), n(shipped.d_exact), n(shipped.d_same_direction))
    fact("The gain is on the benign side: 11,708 more ClinVar benign-side variants", n(b38.same_benign - base37.same_benign), p(base37.blb_recall, 2), p(b38.blb_recall, 2),
         n(base37.same_pathogenic - b38.same_pathogenic), p(base37.plp_recall, 2), p(b38.plp_recall, 2))
    fact("Opposite-direction calls go from 49 to 51", str(int(base37.opposite)), str(int(b38.opposite)), str(int(base37.false_benign)), str(int(base37.false_pathogenic)),
         str(int(b38.false_benign)), str(int(b38.false_pathogenic)), "one call leaves the set" if int(f38v["opposite_left_in_v38"]) == 1 else "?",
         "none enters it" if int(f38v["opposite_entered_in_v38"]) == 0 else "?")
    fact("The call that leaves is RAD21 8:116866762 C>T", f38v["left_opposite_variant"],
         "a splice-acceptor variant" if f38v["left_opposite_consequence"] == "splice_acceptor_variant" else "?",
         "Likely benign at two stars" if f38v["left_opposite_clinvar"] == "Likely_benign" and f38v["left_opposite_stars"] == "2" else "?",
         "Likely pathogenic on PVS1 alone" if f38v["left_opposite_call_v37"] == "Likely_pathogenic" and f38v["left_opposite_criteria_v37"] == "PVS1" else "?",
         "VUS on PVS1 and BS2" if f38v["left_opposite_call_v38"] == "VUS" and f38v["left_opposite_criteria_v38"] == "PVS1&BS2" else "?",
         f"tag of the gene is {f38v['left_opposite_tags_v38']}", n(f38v["left_opposite_gnomad_ac"]), n(f38v["left_opposite_gnomad_an"]), f"{n(f38v['left_opposite_gnomad_hc'])} homozygotes")
    fact("In round 6 the reviewer supported ClinVar's call on this variant", "supported ClinVar" if len(r6[r6.gene == f38v["left_opposite_gene"]]) == 1 and (r6[r6.gene == f38v["left_opposite_gene"]].ruling == "clinvar").all() and int(f38v["left_opposite_in_round6"]) == 1 else "?")
    p37_ar = runs37[(runs37.knob == "pm2_ar_af_threshold") & (runs37.group == "sweep")]
    p38_ar = rv38[(rv38.knob == "pm2_ar_af_threshold") & (rv38.group == "sweep")]
    span_dir = lambda rs, b_: int(max(rs.same_direction.max(), b_.same_direction) - min(rs.same_direction.min(), b_.same_direction))  # noqa: E731
    span_rec = lambda rs, b_: (min(rs.plp_recall.min(), b_.plp_recall), max(rs.plp_recall.max(), b_.plp_recall))  # noqa: E731
    p38_flat = p38_ar.false_benign.nunique() == 1 and p38_ar.false_pathogenic.nunique() == 1 and p38_ar.false_benign.iloc[0] == b38.false_benign and p38_ar.false_pathogenic.iloc[0] == b38.false_pathogenic
    fact("Of the 12 sweeps, 4 differ from v37", n(fam["sweeps"]), n(fam["sweeps_that_differ_in_shape"]), "none differs in the knee" if fam["sweeps_that_differ_in_the_knee"] == 0 else "?",
         sci1(fam["knee_bs1_v38"]), n(fam["knee_bp7_far_boundary_v38"]), "under both versions" if (fam["knee_bs1_v37"], fam["knee_bp7_far_boundary_v37"]) == (fam["knee_bs1_v38"], fam["knee_bp7_far_boundary_v38"]) else "?")
    fact("The BS2 carrier threshold and its floor stop being flat", str(int(rv38.loc["v38_bs2off"].false_benign)), str(int(runs37.loc["v37_bs2off"].false_benign)),
         str(int(rv38.loc["v38_bs2ac_20"].false_benign)), str(int(runs37.loc["v37_bs2ac_20"].false_benign)))
    ix38, ix37 = rv38.loc["v38_ix_table_pm2ad_0"], runs37.loc["v37_ix_table_pm2ad_0"]
    gain_ix38, gain_ix37 = rv38.loc["v38_sw_table_both_sides"].same_pathogenic - ix38.same_pathogenic, runs37.loc["v37_sw_table_both_sides"].same_pathogenic - ix37.same_pathogenic
    fact("The PM2 sweep under the 2015 table loses most of its range", n(gain_ix37), n(gain_ix38), "loses most of its range" if gain_ix38 < gain_ix37 / 2 else "?")
    lo38, hi38 = span_rec(p38_ar, b38)
    fact("The PM2 threshold for genes read as recessive-only reaches 1,310 genes where it reached 67", n(facts38["inheritance_read_recessive_only"]), n(fv37("inheritance_read_recessive_only")),
         p(lo38, 2), p(hi38, 2), str(span_dir(p38_ar, b38)), str(span_dir(p37_ar, base37)), "neither kind of reference-discordant call changes" if p38_flat else "?")
    steps38 = pd.read_csv(DATA / "threshold_steps_v38.csv")
    bs2_38, bs2_37 = (x[x.knob == "BS2 homozygote bar"].reset_index(drop=True) for x in (steps38, steps37))
    fact("The BS2 threshold for individuals with no functional copy keeps the shape it had", *[n(x.gain_calls / x.added_errors) for x in bs2_38.iloc[1:].itertuples()],
         *[n(x.gain_calls / x.added_errors) for x in bs2_37.iloc[1:].itertuples()])
    fact("Of the 245 effects of the single switches", n(fam["effects"]), n(fam["effects_that_changed"]), n(fam["effects_that_changed_and_are_small"]),
         f"on {w(int(fam['effect_measures']))[1]} measures")
    ng38, ng37 = rv38.loc["v38_src_noclingen"], runs37.loc["v37_src_noclingen"]
    fact("The one that is large is the removal of the ClinGen file", n(b38.same_benign - ng38.same_benign), n(base37.same_benign - ng37.same_benign))
    fact("The ClinVar benign-side variants that the removal turns into VUS are 2,007", n(ng38.blb_to_vus - b38.blb_to_vus), n(ng37.blb_to_vus - base37.blb_to_vus),
         n(b38.fire_PM2_Supporting - ng38.fire_PM2_Supporting), n(base37.fire_PM2_Supporting - ng37.fire_PM2_Supporting),
         n(ng38.fire_BS2 - b38.fire_BS2), n(base37.fire_BS2 - ng37.fire_BS2))
    gc38, gc37 = rv38.loc["v38_gap_clingen_complete"], runs37.loc["v37_gap_clingen_complete"]
    fact("The completed ClinGen file changes the correct benign-side calls by 21 under v38", n(gc38.same_benign - b38.same_benign), n(gc37.same_benign - base37.same_benign),
         n(gc38.fire_PM2_Supporting - b38.fire_PM2_Supporting), n(b38.fire_BS2 - gc38.fire_BS2), n(gc37.fire_PM2_Supporting - base37.fire_PM2_Supporting), n(gc37.fire_BS2 - base37.fire_BS2))
    t37, t38 = pd.read_csv(DATA / "temporal_v37.csv"), pd.read_csv(DATA / "temporal_v38.csv")
    ra37, ra38 = pd.read_csv(DATA / "temporal_ratio_v37.csv"), pd.read_csv(DATA / "temporal_ratio_v38.csv")
    a37, a38 = t37[t37.stratum == "all"].iloc[0], t38[t38.stratum == "all"].iloc[0]
    n37, n38 = t37[t37.stratum == "new"].iloc[0], t38[t38.stratum == "new"].iloc[0]
    ao37, ao38 = (pd.read_csv(DATA / f"april_opposite_in_october_{v_}.csv") for v_ in ("v37", "v38"))
    in38 = set(zip(ao38.chrom.astype(str), ao38.pos, ao38.ref, ao38.alt))
    oct_left = ao37[[k_ not in in38 for k_ in zip(ao37.chrom.astype(str), ao37.pos, ao37.ref, ao37.alt)]]
    fact("Against v37, v38 agrees in direction for 2,022 more of the 781,135 covered records", n(a38.same_direction - a37.same_direction), n(a38.stratum_n), p(a37.same_direction / a37.stratum_n, 1),
         p(a38.same_direction / a38.stratum_n, 1), str(int(a38.opposite)), str(int(a37.opposite)), n(n38.stratum_n), p(n37.same_direction / n37.stratum_n, 1), p(n38.same_direction / n38.stratum_n, 1),
         str(int(n38.opposite)), "the same 18" if n37.opposite == n38.opposite else "?",
         f"the {f38v['left_opposite_gene']} variant again" if len(oct_left) == 1 and oct_left.iloc[0].gene == f38v["left_opposite_gene"] and oct_left.iloc[0].october_status == "still 2-star+" else "?",
         f"{ra37[ra37.measure.str.startswith('ratio')].iloc[0].value:.1f} times the April rate under both versions" if f"{ra37[ra37.measure.str.startswith('ratio')].iloc[0].value:.1f}" == f"{ra38[ra38.measure.str.startswith('ratio')].iloc[0].value:.1f}" else "?")

    # ----------------------------------------------------------------------------------- Section 3.8, the baseline and the reading
    fact("they reached 67 and 50 genes", n(names_all.recessive_only), n(names_all.dominant_only))
    fact("Its run differs from v37's in the criteria of 25 indels", n(f38["base_vs_v37_variants_with_other_criteria_or_call"]), n(f38["base_vs_v37_variants_with_another_call"]),
         n(before.exact), n(before.same_direction), n(base37.exact), n(base37.same_direction),
         "indels" if f38["base_vs_v37_differing_variants_that_are_indels"] == f38["base_vs_v37_variants_with_other_criteria_or_call"] else "?",
         "ClinVar VUS variants" if f38["base_vs_v37_differing_variants_called_lp_by_v37_and_vus_in_clinvar"] == f38["base_vs_v37_variants_with_another_call"] else "?")
    fact("Rates are over the 673,660 records, as in the rest of the paper", n(N_ALL), n(pd_row("variants scored")))
    fact("**The field reads 2,056 genes as dominant-only or recessive-only where the names read 117.**",
         n(field_all.dominant_only + field_all.recessive_only), n(names_all.dominant_only + names_all.recessive_only))
    fact("Of the 2,924 genes of the ClinGen file, the names gave 50 dominant-only", n(names_all.genes), n(names_all.dominant_only), n(names_all.recessive_only), n(names_all.both), n(names_all.neither),
         n(field_all.dominant_only), n(field_all.recessive_only), n(field_all.both), n(field_all.neither))
    fact("Among the 2,840 genes with a benchmark variant, the names read as dominant", n(names_b.genes), n(ad_b.genes), n(ad_b.read_dominant_names), n(ar_b.genes), n(ar_b.read_recessive_names),
         n(both_b.genes), n(both_b.read_dominant_names), n(both_b.read_recessive_names))
    fact("Of the other 701 genes, 454 have no established curation", n(reading(ben, "group", group="other or none").genes), n(other("no established curation")), n(other("X-linked only")),
         n(other("mitochondrial only")), n(other("with a semidominant curation")),
         n(sum(other(g) for g in ("no established curation", "X-linked only", "mitochondrial only", "with a semidominant curation", "another combination"))))
    fact("In the file 129 genes, 127 of them with a benchmark variant", n(reading(allf, "contexts").genes), n(reading(ben, "contexts").genes))
    t_ = tr[(tr.run == "v38") & (tr.against == "v38_base") & (tr.clinvar == "B/LB")]
    gained = int(t_[(t_.call_before == "VUS") & t_.call_after.isin(["Benign", "Likely_benign"])].n.sum())
    lost = int(t_[(t_.call_before.isin(["Benign", "Likely_benign"])) & (t_.call_after == "VUS")].n.sum())
    n_blb_all = int(c5.loc[["Likely_benign", "Benign"]].sum().sum())
    fact("v38 gives a benign-side call to 1,924 ClinVar benign-side variants that were VUS", n(gained), n(lost), p(before.blb_recall), p(shipped.blb_recall),
         n(n_blb_all - int(c5.loc[["Likely_benign", "Benign"], "NoCall"].sum())), p(before.same_benign / n_blb_all), p(shipped.same_benign / n_blb_all), n(n_blb_all),
         n(shipped.vus_to_blb - before.vus_to_blb), n(shipped.vus_to_blb), n(before.vus_to_blb), n(shipped.vus_to_plp), "by none" if shipped.vus_to_plp == before.vus_to_plp else "?")
    fact("Reference-discordant benign calls stay at 24, and reference-discordant pathogenic calls fall from 25 to 24", n(before.false_benign), n(before.false_pathogenic), n(shipped.false_pathogenic),
         "stay" if shipped.false_benign == before.false_benign else "?")
    fact("PM2 is met by 9,335 more variants", n(fire.loc["PM2_Supporting", "v38"] - fire.loc["PM2_Supporting", "v38_base"]), n(fire.loc["PM2_Supporting", "v38"]),
         n(fire.loc["PM2_Supporting", "v38_base"]), n(field_all.recessive_only), sci1(float(sw[sw.knob == "pm2_ar_af_threshold"].shipped.iloc[0])),
         sci1(float(sw[sw.knob == "pm2_ad_af_threshold"].shipped.iloc[0])), n(as_v37.d_exact), n(abs(as_v37.d_same_direction)))
    fact("Its carrier arm now reaches the 746 genes read as dominant-only and not 50", n(field_all.dominant_only), n(names_all.dominant_only), n(fire.loc["BS2", "v38"]), n(fire.loc["BS2", "v38_base"]),
         n(fire.loc["BS2", "v38_alt_unguarded"]))

    # ----------------------------------------------------------------------------------- the guard
    u = tr[(tr.run == "v38_alt_unguarded") & (tr.against == "v38")]
    from_lb = int(u[u.call_before == "Likely_benign"].n.sum())
    from_vus = int(u[u.call_before == "VUS"].n.sum())
    by_class = u.groupby("clinvar").n.sum()
    fact("The reading as it first ran, with BS2's carriers counted beside BS1 and BA1, gives 3,432", n(unguarded.exact - shipped.exact), n(unguarded.same_direction - shipped.same_direction))
    fact("It also calls two ClinVar-Pathogenic variants Benign that are VUS in v38", f38["ttr_call_v38_alt_unguarded"], f38["ttr_call_v38"], n(f38["ttr_alleles"]), n(f38["ttr_homozygous"]),
         n(f38["vwf_alleles"]), n(f38["vwf_homozygous"]), "two" if unguarded.false_benign - shipped.false_benign == 2 else "?")
    fact("In genes whose established curations are all autosomal dominant the carrier arm is met on 7,234", n(f38["ad_only_variants_with_bs1_and_bs2"]), n(f38["ad_only_variants_with_bs1"]),
         n(f38["ad_only_variants_with_ba1_and_bs2"]), n(f38["ad_only_variants_with_ba1"]),
         n(int(f38["ad_only_variants_with_bs1_and_bs2"]) + int(f38["ad_only_variants_with_ba1_and_bs2"])),
         "all" if f38["ad_only_variants_with_ba1_and_bs2"] == f38["ad_only_variants_with_ba1"] else "?",
         n(fire.loc["BS2", "v38_alt_unguarded"] - fire.loc["BS2", "v38"]))
    fact("A stand-alone BA1 is Benign already", n(u.n.sum()), n(from_lb), n(from_vus), n(by_class["B/LB"]), n(by_class["VUS"]), n(by_class["P/LP"]))
    fact("It is met by fewer variants than before because the rule that sets the carrier count aside withholds", n(fire.loc["BS2", "v38_alt_unguarded"] - fire.loc["BS2", "v38"]),
         n(fire.loc["BS2", "v38_alt_unguarded"]), n(fire.loc["BS2", "v38"]))
    fact("What the measurement adds is that two of the six have a cause that is not the two variants", n(u.n.sum()))
    fact("Against counting both it gives up 3,432 exact calls", n(unguarded.exact - shipped.exact), n(unguarded.same_direction - shipped.same_direction))
    ns, nu = the_ho("new", "v38_oct_shipped"), the_ho("new", "v38_oct_alt_unguarded")
    fact("On the 116,074 new variants of the October release it gives up 166 exact calls", n(ns.variants), n(nu.exact - ns.exact), n(nu.same_direction - ns.same_direction),
         "avoids no contradiction" if (nu.false_benign, nu.false_pathogenic) == (ns.false_benign, ns.false_pathogenic) else "?")

    # ----------------------------------------------------------------------------------- the four variants
    four = mv[(mv.run == "v38") & (mv.call_after == "VUS")]
    fact("In four loss-of-function variants, NEK1 p.Ser1036Ter", w(len(four))[1], *[f"{g} {h.split(':')[1]}" for g, h in zip(four.gene, four.hgvsp)],
         "PVS1" if set(four.criteria_before) == {"PVS1"} and set(four.criteria_after) == {"PVS1&BS2"} else "?")
    fact("Each is seen in 295 to 324 gnomAD alleles", n(f38["lp_to_vus_alleles_min"]), n(f38["lp_to_vus_alleles_max"]), "none as a homozygote" if int(f38["lp_to_vus_homozygous_max"]) == 0 else "?")
    fact("Withholding the arm for a disorder whose label says susceptibility or predisposition restores two", w(label.same_pathogenic - shipped.same_pathogenic), n(shipped.same_benign - label.same_benign),
         n(shipped.same_direction - label.same_direction))
    oct4 = hmv[(hmv.run == "v38_oct_shipped") & (hmv.call_after == "VUS")]
    fact("The same four variants lose the call on the October release", "same" if sorted(oct4.gene) == sorted(four.gene) else "?")

    # ----------------------------------------------------------------------------------- the cap
    a = tp.loc["April set"]
    o = tp.loc["October, all"]
    nn = tp.loc["October, new variants"]
    fact("**A combined cap on PP2, PM1 and PP3 trades 320 correct pathogenic-side calls", n(a.clinvar_plp_losing_the_call), n(a.clinvar_vus_kept_at_vus), n(a.clinvar_blb_no_longer_called_plp))
    fact("All three are met for 1,592 variants", n(a.variants_with_pp2_pm1_and_pp3_met), n(a.of_them_called_likely_pathogenic))
    fact("Capping the three together at 4 points changes 910 calls", n(a.calls_changed), n(a.clinvar_plp_losing_the_call), n(a.clinvar_vus_kept_at_vus), n(a.clinvar_blb_no_longer_called_plp),
         "all from Likely pathogenic to VUS" if a.calls_changed == a.all_likely_pathogenic_to_vus else "?")
    fact("It raises agreement in direction by 267 calls", n(triple.same_direction - shipped.same_direction), n(triple.exact - shipped.exact), p(shipped.plp_recall), p(triple.plp_recall),
         n(shipped.false_pathogenic - triple.false_pathogenic), n(shipped.false_pathogenic), n(a.plp_calls_lost_per_false_pathogenic_avoided))
    fact("On the October release the same cap changes 1,019 calls", n(o.calls_changed), n(o.clinvar_plp_losing_the_call), n(o.clinvar_blb_no_longer_called_plp), n(nn.calls_changed),
         n(nn.clinvar_plp_losing_the_call), "the one" if nn.clinvar_blb_no_longer_called_plp == 1 else "?")
    fact("We did not adopt it, because it changes the combining rule", n(a.clinvar_vus_kept_at_vus), n(a.plp_calls_lost_per_false_pathogenic_avoided))

    # ----------------------------------------------------------------------------------- the calls for each mode of inheritance
    fact("**A gene with more than one mode of inheritance is classified once for each", f"{pd_row('variants whose modes of inheritance give different classes') / pd_row('variants whose gene has more than one mode of inheritance') * 100:.2f}%")
    fact("Of the 673,278 variants with a call, 87,599 (13.0%) are in 128 genes", n(pd_row("variants scored")), n(pd_row("variants whose gene has more than one mode of inheritance")),
         p(pd_row("variants whose gene has more than one mode of inheritance") / pd_row("variants scored")), n(pd_row("genes among them")))
    different = pd_row("variants whose modes of inheritance meet different criteria")
    by_bs2, by_pm2 = pd_row("  of them, BS2 is met under one mode and not another"), pd_row("  of them, PM2_Supporting is met under one mode and not another")
    fact("The modes meet different criteria for 4,437 of the variants", n(different), n(by_bs2), n(by_pm2), "always" if by_bs2 + by_pm2 == different else "?")
    fact("They give different classes for 492 variants (0.56%)", n(pd_row("variants whose modes of inheritance give different classes")),
         f"{pd_row('variants whose modes of inheritance give different classes') / pd_row('variants whose gene has more than one mode of inheritance') * 100:.2f}%",
         n(pd_row("variants whose modes of inheritance give different directions")),
         f"{pd_row('variants whose modes of inheritance give different directions') / pd_row('variants whose gene has more than one mode of inheritance') * 100:.2f}%")
    c_cls, c_dir = "  of them, modes of inheritance give different classes", "  of them, modes of inheritance give different directions"
    fact("Among the ClinVar classes of such genes the modes give different classes for 15", n(pd_row(c_cls, 0)), n(pd_row("ClinVar P/LP variants in such genes")), n(pd_row(c_cls, 1)),
         n(pd_row("ClinVar VUS variants in such genes")), n(pd_row(c_cls, 2)), n(pd_row("ClinVar B/LB variants in such genes")))
    fact("In no gene do they differ for more than 30 variants", n(pdg.differ_in_class.max()), *[f"{g} {int(pdg.loc[g, 'differ_in_class'])}" for g in ("PKD1", "BRCA2", "RYR1")])

    # ----------------------------------------------------------------------------------- the sweeps
    # The sweeps of the four thresholds that read inheritance are checked in Section 3.3 (`thresholds`) and, against v37, above; what is left here
    # is that the runs of Section 3.8 that made them are the family's own.
    pm2_ad = sw[sw.knob == "pm2_ad_af_threshold"]
    pm2_ar = sw[sw.knob == "pm2_ar_af_threshold"]
    for knob, x in (("pm2_ad_af_threshold", pm2_ad), ("pm2_ar_af_threshold", pm2_ar)):
        if not all(x.index.isin(rv38.index)) or not all((rv38.loc[x.index, "same_direction"] == x.same_direction).values):
            failures.append(f"the sweep of {knob} in v38_sweeps.csv is not the one of the family v38 in runs_v38_variants.csv")

    # ----------------------------------------------------------------------------------- the October release
    ao, bo = the_ho("all", "v38_oct_shipped"), the_ho("all", "v38_oct_base")
    uo = the_ho("all", "v38_oct_alt_unguarded")
    nb = the_ho("new", "v38_oct_base")
    fact("On the 780,749 covered records v38 gives 336 more exact calls", n(ao.variants), n(ao.exact - bo.exact), n(ao.same_direction - bo.same_direction), n(bo.false_benign), n(ao.false_benign),
         n(bo.false_pathogenic), n(ao.false_pathogenic))
    fact("**On the new variants of the October release the effect is", n(ns.same_direction - nb.same_direction))
    fact("with 30 reference-discordant benign calls before and after and 39 reference-discordant pathogenic calls", n((ao.same_direction - bo.same_direction) - (ns.same_direction - nb.same_direction)))
    fact("On the 116,074 new variants with a call it gives 154 more exact calls", n(ns.variants), n(ns.exact - nb.exact), n(ns.same_direction - nb.same_direction), n(nb.false_benign), n(nb.false_pathogenic),
         "before and after" if (ns.false_benign, ns.false_pathogenic) == (nb.false_benign, nb.false_pathogenic) else "?")
    fact("Counting BS2 beside BS1 and BA1 adds 3,943 exact calls", n(uo.exact - bo.exact), n(uo.same_direction - bo.same_direction), "two" if uo.false_benign - bo.false_benign == 2 else "?",
         "both of them old variants" if set(hmv[(hmv.run == "v38_oct_alt_unguarded") & (hmv.call_after == "Benign")].stratum) == {"unchanged"} else "?")
    fact("v38 was developed on the April set, and only the 116,074 new variants", n(ns.variants), n(ns.same_direction - nb.same_direction), n(ao.same_direction - bo.same_direction))

    # ----------------------------------------------------------------------------------- gnomAD's joint release
    from decimal import ROUND_HALF_UP, Decimal

    def pct3(x) -> str:
        """A frequency as a percentage with three decimals, rounded half up from the digits that the file holds."""
        return f"{(Decimal(str(x)) * 100).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP)}%"

    j = jn.set_index("run")
    ex, nf, jt = j.loc["v38"], j.loc["v38_noflags"], j.loc["v38_joint"]
    every = jc[(jc.set == "every chromosome") & (jc.run == "v38_joint")]
    moved = every[every.consequence_group != "(criterion)"]
    by_class = moved.groupby("clinvar").n.sum()
    rec = csv("v38_joint_records.csv").set_index("set").loc["every chromosome"]
    cov = csv("v38_joint_coverage.csv").set_index("group")
    added_all = (jt.false_benign - ex.false_benign) + (jt.false_pathogenic - ex.false_pathogenic)
    added_genomes = (jt.false_benign - nf.false_benign) + (jt.false_pathogenic - nf.false_pathogenic)
    price = lambda gain, cost: n(round(gain / cost, -2))

    fact("On every chromosome (673,278 variants) a database built from it holds a record for every variant that the exomes hold",
         n(rec.variants), n(rec.in_both), n(rec.joint_only), p(rec.joint_only / rec.variants),
         "never lower" if (rec.exomes_only, rec.ac_lower, rec.an_lower, rec.hc_lower) == (0, 0, 0, 0) else "?")
    fact("It changes 13,412 calls (2.0%)", n(moved.n.sum()), p(moved.n.sum() / jt.variants), n(by_class["B/LB"]), n(by_class["VUS"]), n(by_class["P/LP"]),
         n(jt.d_exact), pp(jt.d_exact / jt.variants, 2), n(jt.d_same_direction), pp(jt.d_same_direction / jt.variants, 2),
         n(ex.false_benign), n(jt.false_benign), n(ex.false_pathogenic), n(jt.false_pathogenic),
         "one more reference-discordant benign call" if jt.false_benign - ex.false_benign == 1 else "?",
         "two more reference-discordant pathogenic calls" if jt.false_pathogenic - ex.false_pathogenic == 2 else "?")
    fact("ClinVar benign-side variants called benign-side rise from 188,914 to 198,209", n(ex.same_benign), n(jt.same_benign), n(ex.same_pathogenic), n(jt.same_pathogenic))
    intronic, utr = cov.loc["intronic variants whose call moves"], cov.loc["UTR variants whose call moves"]
    fact("Of the 13,412 calls that move, 4,485 are on variants that the exomes lack", n(cov.loc["calls that move", "variants"]), n(cov.loc["calls that move", "no_exomes_record"]),
         n(intronic.no_exomes_record), n(intronic.variants), n(utr.variants), n(utr.median_an_exomes), n(utr.median_an_joint),
         "of the 6,476" if int(moved[moved.consequence_group == "Intronic"].n.sum()) == intronic.variants else "?",
         "UTR variants" if int(moved[moved.consequence_group == "UTR"].n.sum()) == utr.variants else "?")
    ba1 = every[every.criterion == "BA1"].iloc[0]
    fact("BA1 is gained on 10,115 of the moved calls and lost on 339", n(ba1.criterion_gained), n(ba1.criterion_lost))
    bc = csv("v38_joint_by_chrom.csv").set_index("chrom")
    autosomes_and_x = bc.loc[[str(c) for c in range(1, 23)] + ["X"]]
    share = autosomes_and_x.calls_moved / autosomes_and_x.variants
    fact("Each of the 22 autosomes and the X gains calls", p(share.min()), p(share.max()), n(bc.loc["Y", "variants"]), n(bc.loc["MT", "variants"]),
         "gains calls" if bool(((autosomes_and_x.d_exact > 0) & (autosomes_and_x.d_same_direction > 0)).all()) else "?",
         "nothing moves" if int(bc.loc[["Y", "MT"], "calls_moved"].sum()) == 0 else "?")
    fact("Part of the gain is the flags that the joint release lacks and not the genomes", n(nf.d_exact), n(nf.d_same_direction), n(jt.d_exact - nf.d_exact),
         n(jt.d_same_direction - nf.d_same_direction),
         "one more reference-discordant pathogenic call and no more reference-discordant benign call" if (nf.false_pathogenic - ex.false_pathogenic, nf.false_benign - ex.false_benign) == (1, 0) else "?",
         "one more of each" if (jt.false_benign - nf.false_benign, jt.false_pathogenic - nf.false_pathogenic) == (1, 1) else "?")
    pc = csv("v38_joint_pilot_check.csv")
    fact("A pilot of the same three runs on chromosomes 21 and 22 alone", n(pc.variants.iloc[0]),
         "gives the same call and the same criteria"
         if len(pc) == 3 and (pc.variants == pc.found_in_full).all() and not pc.calls_differing.any() and not pc.criteria_differing.any() else "?")

    je = csv("v38_joint_errors.csv")
    for flag in ("filter_ac0", "filter_vqsr", "lcr", "segdup", "bs2_met", "pm2_met", "bs1_met", "capped_by_veto"):
        je[flag] = je[flag].eq(True)  # a blank (no gnomAD record) is False
    errors = je[~je.kind.str.startswith("pathogenic-side")]
    one = lambda symbol, database: errors[(errors.symbol == symbol) & (errors.database == database)].iloc[0]
    faf = lambda symbol, database: float(one(symbol, database).faf95_max)
    kinds = errors.drop_duplicates(["chrom", "pos"]).kind
    fact("Five variants change their standing as an opposite-direction call, four added and one removed", w(len(kinds)),
         w(int(kinds.str.startswith("no error").sum())), w(int(kinds.str.endswith("no error").sum())))
    bar = float(one("KMT2C", "joint").bs1_bar)
    straddle = {s for s in ("KMT2C", "RAI1", "PIEZO1", "TXNL4A", "NF1") if (faf(s, "exomes") - bar) * (faf(s, "joint") - bar) < 0}
    fact("Four of the five are at the threshold of BS1, a filtering allele frequency of 0.5%", w(len(straddle)), p(bar),
         "the fifth belongs to the missing flags" if straddle == {"KMT2C", "RAI1", "PIEZO1", "TXNL4A"} else "?")
    km, ra = one("KMT2C", "exomes"), one("RAI1", "exomes")
    fact("KMT2C p.Cys391Ter (ClinVar Likely benign) and RAI1 p.Gln280HisfsTer83 (Benign)", km.hgvsp, ra.hgvsp,
         "variant-quality filter" if all(one(s, d).frequency_withheld_by == "FILTER=AS_VQSR" for s in ("KMT2C", "RAI1") for d in ("exomes", "joint")) else "?",
         "at VUS" if (km.call, ra.call) == ("VUS", "VUS") and km.capped_by_veto and ra.capped_by_veto and km.criteria == ra.criteria == "PVS1" else "?",
         "Likely benign" if km.clinvar == "Likely_benign" else "?", "Benign" if ra.clinvar == "Benign" else "?")
    fact("Their filtering allele frequencies are 0.522% and 0.537% in the exomes", pct3(faf("KMT2C", "exomes")), pct3(faf("RAI1", "exomes")), pct3(faf("KMT2C", "joint")),
         pct3(faf("RAI1", "joint")),
         "both are called Likely pathogenic" if all(one(s, "joint").call == "Likely_pathogenic" and not one(s, "joint").capped_by_veto for s in ("KMT2C", "RAI1")) else "?")
    fact("PIEZO1's splice-acceptor deletion (ClinVar Benign) is the same case in the other direction", pct3(faf("PIEZO1", "exomes")), pct3(faf("PIEZO1", "joint")),
         "Likely pathogenic" if one("PIEZO1", "exomes").call == "Likely_pathogenic" else "?", "at VUS" if one("PIEZO1", "joint").call == "VUS" and one("PIEZO1", "joint").capped_by_veto else "?",
         "splice-acceptor" if one("PIEZO1", "exomes").consequence == "splice_acceptor_variant" else "?", "Benign" if one("PIEZO1", "exomes").clinvar == "Benign" else "?")
    tx_e, tx_j = one("TXNL4A", "exomes"), one("TXNL4A", "joint")
    fact("TXNL4A's deletion upstream of the gene (ClinVar Pathogenic)", n(tx_e.gnomad_an), pct3(tx_e.faf95_max), n(tx_j.gnomad_an), pct3(tx_j.faf95_max),
         "BS1 fires" if tx_j.bs1_met and not tx_e.bs1_met else "?", "Likely benign" if tx_j.call == "Likely_benign" else "?",
         "upstream" if tx_e.consequence == "upstream_gene_variant" else "?", "Pathogenic" if tx_e.clinvar == "Pathogenic" else "?")
    nf_e, nf_j = one("NF1", "exomes"), one("NF1", "joint")
    fact("NF1 p.Met895Ile (ClinVar Likely benign) is called Likely pathogenic", nf_e.hgvsp,
         "segmental-duplication flag" if nf_e.frequency_withheld_by == "segmental duplication" and nf_e.segdup else "?",
         "fires" if nf_j.pm2_met and not nf_e.pm2_met and nf_j.call == "Likely_pathogenic" else "?", "Likely benign" if nf_e.clinvar == "Likely_benign" else "?")

    lost_rows = je[je.kind.str.startswith("pathogenic-side")]
    lost_j, lost_e = lost_rows[lost_rows.database == "joint"], lost_rows[lost_rows.database == "exomes"]
    plp = every[every.clinvar == "P/LP"]
    is_p = lambda s: s.isin(["Pathogenic", "Likely_pathogenic"])
    lost_n = int(plp[is_p(plp.call_before) & ~is_p(plp.call_after)].n.sum())
    gained_n = int(plp[~is_p(plp.call_before) & is_p(plp.call_after)].n.sum())
    fact("Thirteen ClinVar pathogenic-side variants lose a pathogenic-side call and 39 gain one", w(lost_n), n(gained_n),
         "lose" if lost_n == len(lost_j) else "?")
    filtered = lost_j[lost_j.frequency_withheld_by.fillna("").str.startswith("FILTER=")]
    in_exomes = lost_e.set_index(["chrom", "pos"]).loc[list(zip(filtered.chrom, filtered.pos))]
    has_record = in_exomes.gnomad_ac.notna()
    fact("Ten of the 13 are missense variants whose record in the joint release carries a filter of gnomAD's quality control", w(len(filtered)), n(len(lost_j)),
         w(int((filtered.filter_ac0 & ~filtered.filter_vqsr).sum())), w(int((~filtered.filter_ac0 & filtered.filter_vqsr).sum())),
         w(int((filtered.filter_ac0 & filtered.filter_vqsr).sum())),
         w(int(has_record.sum())), w(int((~has_record).sum())),
         "one to four alleles" if int(in_exomes[has_record].gnomad_ac.min()) == 1 and int(in_exomes[has_record].gnomad_ac.max()) == 4 else "?",
         "clean record" if not in_exomes[has_record].frequency_withheld_by.fillna("").str.len().any() else "?",
         "missense" if (filtered.consequence == "missense_variant").all() else "?", "withholds PM2" if not filtered.pm2_met.any() else "?")
    low = lost_e[lost_e.frequency_withheld_by.fillna("") == "low-complexity region"]
    fact("Three are frameshifts, MYO6 p.Gln918ThrfsTer24 and PRRT2 p.Arg217ProfsTer8 and p.Arg217GlufsTer12", w(len(low)), *sorted(low.hgvsp),
         "frameshifts" if (low.consequence == "frameshift_variant").all() else "?",
         "BS2" if lost_j[lost_j.hgvsp.isin(low.hgvsp)].bs2_met.all() and not lost_j[lost_j.hgvsp.isin(low.hgvsp)].lcr.any() else "?")

    pm2 = csv("v38_joint_pm2_lost.csv")
    ac0 = pm2[(pm2.reason == "filter AC0") & pm2.joint_count_positive].iloc[0]
    fact("Against the exomes without the flags, 741 variants with a positive joint count lose PM2 to an AC0 in one callset", n(ac0.variants), n(ac0.calls_moved),
         n(ac0.pathogenic_side_calls_moved))
    fact("The change costs about 4,600 calls in direction for each added opposite-direction call", price(jt.d_same_direction, added_all), n(jt.d_same_direction), w(added_all),
         price(jt.d_same_direction - nf.d_same_direction, added_genomes), n(jt.d_same_direction - nf.d_same_direction), w(added_genomes))
    fact("What the joint release gives up is the flags", n(nf.d_same_direction), "one reference-discordant pathogenic call" if nf.false_pathogenic - ex.false_pathogenic == 1 else "?",
         "three pathogenic-side frameshifts" if len(low) == 3 else "?")
    fact("gnomAD's genomes are not read by the benchmark; read through the joint release on every chromosome they give 9,187 more calls in direction",
         n(jt.d_same_direction), w(added_all))
    fact("Seventh, gnomAD's genomes, read through the joint release", n(jt.d_same_direction), w(added_all), price(jt.d_same_direction, added_all))
    fact("Reading gnomAD's genomes through its joint release adds 9,187 calls in direction on every chromosome", n(jt.d_same_direction), w(added_all))
    fact("Reading gnomAD's genomes as well as its exomes, through its joint release, would add 9,321 calls in direction", n(jt.d_same_direction), pp(jt.d_same_direction / jt.variants, 2),
         "one more reference-discordant benign call" if jt.false_benign - ex.false_benign == 1 else "?", "two more reference-discordant pathogenic calls" if jt.false_pathogenic - ex.false_pathogenic == 2 else "?")

    # ----------------------------------------------------------------------------------- the sentences elsewhere that quote it
    fact("Against v37, v38 reads the mode of inheritance that ClinGen curates, which makes 2,056 genes dominant-only or recessive-only", n(names_all.dominant_only + names_all.recessive_only),
         n(field_all.dominant_only + field_all.recessive_only))
    fact("It adds 1,494 calls in direction (0.22 points) with no added reference-discordant benign call", n(b38.same_direction - base37.same_direction),
         pp((b38.same_direction - base37.same_direction) / N_ALL, 2),
         "no added reference-discordant benign call" if b38.false_benign == base37.false_benign else "?",
         "one fewer reference-discordant pathogenic call" if base37.false_pathogenic - b38.false_pathogenic == 1 else "?",
         f"{w(len(four))[1]} ClinVar-Pathogenic loss-of-function variants move from Likely pathogenic to VUS")
    fact("Counting one gnomAD observation as two benign criteria would have added", n(unguarded.d_same_direction - shipped.d_same_direction),
         "two ClinVar-Pathogenic variants Benign" if unguarded.false_benign - shipped.false_benign == 2 else "?")
    fact("**A flat sweep can mean that the setting reaches almost nothing, or that the rule beside it makes the setting irrelevant, and an input read from prose can hide the first.**")
    fact("The sweep of the BS2 carrier threshold was flat in v37 because", n(names_all.dominant_only), n(field_all.dominant_only))
    fact("The sweeps of the PM2 thresholds were flat for the other reason", n(names_all.recessive_only), n(field_all.recessive_only))
    fact("Across the two thresholds the words made", n(names_all.dominant_only + names_all.recessive_only), n(names_all.genes), n(field_all.dominant_only + field_all.recessive_only))
    fact("Read the field that a source curates and not the prose beside it", n(names_all.dominant_only + names_all.recessive_only), n(field_all.dominant_only + field_all.recessive_only), n(len(rows)))
    t37o, t38o = pd.read_csv(DATA / "temporal_v37.csv"), pd.read_csv(DATA / "temporal_v38.csv")
    a37o, a38o = t37o[t37o.stratum == "all"].iloc[0], t38o[t38o.stratum == "all"].iloc[0]
    n37o, n38o = t37o[t37o.stratum == "new"].iloc[0], t38o[t38o.stratum == "new"].iloc[0]
    fact("Against v37, v38 reads the mode of inheritance that ClinGen curates, which v37 took", n(b38.same_direction - base37.same_direction), n(a38o.same_direction - a37o.same_direction),
         n(n38o.same_direction - n37o.same_direction), w(len(four))[1], n(unguarded.d_same_direction - shipped.d_same_direction),
         "no added reference-discordant benign call" if b38.false_benign == base37.false_benign and a38o.false_benign == a37o.false_benign else "?",
         "one fewer reference-discordant pathogenic call" if base37.false_pathogenic - b38.false_pathogenic == 1 else "?")
    fact("The mode of inheritance that BS2 and PM2 read was the weakest input of v37", "")
    fact("The field gives a mode and not a penetrance", w(len(four)))
    # what BS2's arm for homozygotes and hemizygotes still counts beside BS1 and BA1, and the price of the label guard
    fact("**The arm for homozygotes and hemizygotes is not set aside, and the same doubling remains in it for", n(f38["bs2_beside_bs1_without_ba1_variants"]))
    fact("BS2 is met beside BS1 or BA1 through that arm on", n(f38["bs2_beside_frequency_variants"]), n(f38["bs2_beside_frequency_variants_on_x"]))
    fact("BA1 alone is Benign, so the doubling changes a call only beside BS1 without BA1", n(f38["bs2_beside_bs1_without_ba1_variants"]),
         n(f38["bs2_beside_bs1_without_ba1_on_x"]), n(f38["bs2_beside_bs1_without_ba1_clinvar_benign_side"]), n(f38["bs2_beside_bs1_without_ba1_clinvar_vus"]))
    fact("No ClinVar pathogenic-side variant is among them, so the benchmark shows no harm from it",
         "No ClinVar pathogenic-side variant" if int(f38["bs2_beside_bs1_without_ba1_clinvar_pathogenic_side"]) == 0 else "?")
    fact("BS2's arm for homozygotes and hemizygotes still counts", n(f38["bs2_beside_bs1_without_ba1_variants"]))
    per_restored = n(round((shipped.same_benign - label.same_benign) / (label.same_pathogenic - shipped.same_pathogenic)))
    fact("The label guard that restores two of the four gives up about", w(label.same_pathogenic - shipped.same_pathogenic), per_restored)
    fact("and in v38 we declined the label guard that restores two ClinVar-Pathogenic variants at about", per_restored)


def structure() -> None:
    """The figures and tables of the paper: what it cites exists, comes in order, and every data file is in exactly one table.

    The main figures state conclusions and are cited in order (1 to 5); the supplementary figures (S1 to S4) show how the
    classifier was developed and reviewed. There are fewer than ten supplementary tables, and every other table is in the
    appendix, so each file of `data/` must be named by exactly one supplementary or appendix table, and the workbook that
    gathers the tables must hold what those files say.
    """
    text = (HERE / "manuscript.md").read_text()
    body, back = text.split("\n## Figure legends\n", 1)
    legends = re.findall(r"^\*\*Figure (S?\d)\. ", back, re.M)
    sup = re.findall(r"^\*\*Supplementary Table (S\d)\.", back, re.M)
    app = re.findall(r"^\*\*Appendix Table (A\d+)\.", back, re.M)

    def check(ok: bool, what: str) -> None:
        (passed if ok else failures).append(("structure: " if ok else "structure: ") + what)

    check(legends == ["1", "2", "3", "4", "5", "S1", "S2", "S3", "S4"], f"the legends are the main figures 1 to 5 and the supplementary figures S1 to S4 (found {legends})")
    check(len(sup) < 10 and sup == [f"S{i}" for i in range(1, len(sup) + 1)], f"fewer than ten supplementary tables, numbered S1 to S{len(sup)}")
    check(app == [f"A{i}" for i in range(1, len(app) + 1)] and len(app) > 0, f"the appendix tables are numbered A1 to A{len(app)}")
    main_tables = re.findall(r"^\*\*Table (\d+)\.", body, re.M)
    check(main_tables == [str(i) for i in range(1, len(main_tables) + 1)] and len(main_tables) > 0, f"the tables of the text are numbered 1 to {len(main_tables)} in the order they stand")
    prose = "\n".join(line for line in body.splitlines() if not line.startswith("**Table "))
    first_cited = list(dict.fromkeys(re.findall(r"\bTable (\d+)\b", prose)))
    check(first_cited == main_tables, f"the tables of the text are first cited in the order they are numbered (cited {first_cited}, numbered {main_tables})")
    for kind, pattern, defined in (("figure", r"Figure (S?\d)", legends), ("supplementary table", r"Supplementary Table (S\d)", sup), ("appendix table", r"Appendix Table (A\d+)", app)):
        cited = re.findall(pattern, text)
        missing = sorted(set(cited) - set(defined))
        check(not missing, f"every {kind} that the text cites is defined (undefined: {missing})")
        in_body = list(dict.fromkeys(re.findall(pattern, body)))
        check(set(in_body) == set(defined), f"every {kind} is cited in the text before the legends (not cited: {sorted(set(defined) - set(in_body))})")
        groups = [[x for x in in_body if x.startswith("S")], [x for x in in_body if not x.startswith("S")]] if kind == "figure" else [in_body]
        for g in groups:
            check(g == sorted(g, key=lambda x: int(x.lstrip("SA"))), f"{kind}s are first cited in order ({g})")
    entries = [para for para in back.split("\n\n") if para.startswith(("**Supplementary Table", "**Appendix Table"))]
    listed = [f for para in entries for f in sorted(set(re.findall(r"`([\w.\-]+\.(?:csv|tsv))`", para)))]
    files = sorted(f.name for f in DATA.iterdir() if f.suffix in (".csv", ".tsv"))
    if BASE == "v38":
        files = [f for f in files if not ("_v37" in f and (DATA / f.replace("_v37", "_v38")).exists())]  # the files of v37 are in no table of a manuscript that says v38
    else:
        listed = [f.replace("_v38", f"_{BASE}") if (DATA / f.replace("_v38", f"_{BASE}")).exists() else f for f in listed]  # the manuscript names v38's files: their twins stand in
        files = [f for f in files if "_v38" not in f]
    dup = sorted({f for f in listed if listed.count(f) > 1})
    check(not dup, f"no data file is in two tables (twice: {dup})")
    check(sorted(set(listed)) == files, f"every data file is in a table, and every file named exists (not listed: {sorted(set(files) - set(listed))}; not found: {sorted(set(listed) - set(files))})")
    # The workbook gathers those files by table, so a change to a file or to a table entry that it has not followed is a stale table.
    if BASE == "v38":
        import build_supplement

        stale = build_supplement.differences()
        check(not stale, f"supplementary_tables.xlsx holds what data/ and the table entries say (rerun build_supplement.py --write; {stale[:1]})")


def binaries() -> None:
    """`binaries.tsv` and `binaries_v38.tsv` name every binary that made a run of the repository, with the runs it made.

    The provenance file of a run records the SHA-256 of its binary and the commit of the checkout that ran the harness, which is
    not the commit of the binary; the table is where a hash is tied to what it was built from. A hash that is in a provenance
    file and not in the table, a hash in the table that made no run, and a count that is not the number of provenance files
    with the hash are each a failure.
    """
    res = HERE.parents[1] / "analysis/acmg_benchmark/results/v37_manuscript"
    table = pd.concat([pd.read_csv(res / name, sep="\t", dtype=str) for name in ("binaries.tsv", "binaries_v38.tsv") if (res / name).exists()], ignore_index=True)
    seen: dict[str, int] = {}
    for f in sorted((res / "provenance").glob("*.txt")):
        for line in f.read_text().splitlines():
            if line.startswith("fastvep_sha256: "):
                sha = line.split(": ", 1)[1].strip()
                seen[sha] = seen.get(sha, 0) + 1
    listed = dict(zip(table.sha256, table.runs.astype(int)))
    for sha, k in sorted(seen.items()):
        if sha not in listed:
            failures.append(f"binaries.tsv: {sha[:12]} made {k} runs and is not in the table")
        elif listed[sha] != k:
            failures.append(f"binaries.tsv: {sha[:12]} made {k} runs and the table says {listed[sha]}")
    for sha in listed:
        if sha not in seen:
            failures.append(f"binaries.tsv: {sha[:12]} is in the table and made no run")
    passed.append(f"binaries.tsv: {len(seen)} binaries made {sum(seen.values())} runs")


CHECKS = [databases, methods, trajectory, review, thresholds, points, ablations, discordance, residual, later_release, discussion, abstract_and_intro, contribution, repository, availability, audit, early_trades, criteria_consistency, coverage, v38, binaries, structure]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-v", action="store_true", help="print every passed check")
    ap.add_argument("--base", choices=["v37", "v38"], default="v38", help="the classifier whose files are checked (default v38; see the docstring; read before the arguments)")
    ap.add_argument("--dump", default="", help="write every check (the line, the anchor and the values the data give) to this JSON file")
    args = ap.parse_args()
    for check in CHECKS:
        try:
            check()
        except Exception as e:  # a check that cannot run is a failure, not a skip
            failures.append(f"{check.__name__}: {type(e).__name__}: {e}")
    if args.dump:
        import json

        Path(args.dump).write_text(json.dumps(checks_made, indent=1))
    if args.v:
        print("\n".join(passed))
    print(f"{len(passed)} values found in their sentences, {len(failures)} failures")
    for f in failures:
        print("FAIL", f)
    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
