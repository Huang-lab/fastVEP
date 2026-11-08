#!/usr/bin/env python3
"""Set the runs of the family v38 against the runs of the family v37, pair by pair.

The two families are the same list of runs (`11_manuscript_runs.py --family v37` and `--family v38`: the sweeps, switches,
leave-one-source-out ablations, PM2 sweeps under the 2015 table, point cutoffs, hooks and gap runs on the April set, and the
seven changes on the October release), made by the classifier of the review rounds (v37) and by the classifier of release 0.5.0 (v38, the
latest, which is the one that the manuscript describes), so a
pair of runs with the same name apart from the prefix differs in the classifier and in nothing else. This script reads the
scored files of both and writes what moved:

  results/v38_family_vs_v37.csv            one row for every pair, run, number and cell: the 30 cells of the concordance matrix, the
                                           headline metrics and the firing count of every criterion, under v37 and under v38, and
                                           the difference (tidy: kind, item, v37, v38, delta, relative_delta)
  results/v38_family_vs_v37_sweeps.csv     every sweep as a curve, one row for every setting: the metrics under both and the
                                           differences (see `sweep_shape` for the shape of a sweep: its ranges, the sign of each
                                           step, its knee and the order of its settings)
  results/v38_family_vs_v37_effects.csv    every ablation, switch, hook and gap run as an effect (the run less the default under the
                                           same classifier) under both, and whether the effect changed sign or by more than 20%
  results/v38_family_vs_v37_opposite.csv   the opposite-direction calls that are in one run of a pair and not in the other, with
                                           the criteria and the call each has
  results/v38_family_vs_v37_summary.csv    the counts that the manuscript quotes from this comparison and that no table of
                                           `manuscript/acmg/data/` holds (`verify_numbers.py` reads them): the sweeps and how many
                                           differ in shape or knee, the effects and how many changed
  results/v38_family_vs_v37.md             the same in words, one sentence on a line: a summary, then the main pair, the
                                           opposite-direction calls, the sweeps, the ablations, the controls, the criteria that fire
                                           more or less and the October release, and last the statements of the manuscript that read
                                           differently under v37 (Section 9: the claims whose truth or wording differs, with their
                                           line numbers and the numbers that decide them; the statements about how the runs were
                                           made; the tables, workbook and figures that differ; and every number of the text that
                                           `verify_numbers.py` recomputes differently)

Everything is over all 673,660 records of the April set and all 781,135 of the October release, as the paper's tables are: the
reference-agreement records (ALT is a dot; 382 and 385) are no-call under both classifiers, v37's binaries write a record for
them and v38's do not, so the matrices of v38 are filled with them as no-calls (`13_manuscript_stats.py`, `fill_nocall`).

ClinVar's 2-star+ classifications are a reference and not ground truth. Every agreement here is with them and nothing stronger.

Section 9 of the document reads `manuscript/acmg/manuscript.md`, which says v38, and does not change it. Its line numbers are those of
the file as it stands, and it runs `verify_numbers.py --base v37` and `--base v38` with `--dump` and sets the two dumps side by side (the
numbers of the text that differ under v37), builds the workbook and the figures of v37 in a temporary directory, and compares the tables
that `build_tables.py --base v37` makes with the text. The statements that it quotes are those of the manuscript as it stands, and the
wording that they had when the manuscript was written for v37 is in the script. The few
statements of the text above that are typed and not computed (which sweeps have a knee, which four variants lose a pathogenic-side
call, which sweeps are flat, which is the large effect) are asserted by `check_typed_claims`, so the script stops and does not write a
document that the data contradict.

Usage
  20_compare_families.py          write the files (needs the runs of both families, their per-variant tables, and
                                  `manuscript/acmg/data/*_v37.csv` and `*_v38.csv`; takes about a minute)
"""

from __future__ import annotations

import csv
import datetime as dt
import importlib.util
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from zoneinfo import ZoneInfo

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BENCH = ROOT / "data/benchmark"
RES = ROOT / "analysis/acmg_benchmark/results"
RUNS_DIR = RES / "v37_manuscript"
DATA = ROOT / "manuscript/acmg/data"
OUT_CSV = RES / "v38_family_vs_v37.csv"
OUT_SWEEPS = RES / "v38_family_vs_v37_sweeps.csv"
OUT_EFFECTS = RES / "v38_family_vs_v37_effects.csv"
OUT_OPPOSITE = RES / "v38_family_vs_v37_opposite.csv"
OUT_SUMMARY = RES / "v38_family_vs_v37_summary.csv"
OUT_MD = RES / "v38_family_vs_v37.md"

_spec = importlib.util.spec_from_file_location("stats13", HERE / "13_manuscript_stats.py")
st = importlib.util.module_from_spec(_spec)
sys.modules["stats13"] = st
_spec.loader.exec_module(st)
st.set_base("v38")  # `fill_nocall` puts the reference-agreement records back into v38's matrices; v37's already hold them, which it checks

CLASSES, CALLS = st.CLASSES, st.CALLS
METRICS = ["exact", "same_direction", "same_pathogenic", "same_benign", "vus_exact", "false_benign", "false_pathogenic", "opposite", "plp_recall",
           "blb_recall", "vus_recall", "plp_to_vus", "blb_to_vus", "vus_to_plp", "vus_to_blb", "no_call"]
KEY = ["chrom", "pos", "ref", "alt"]
APRIL_N, OCTOBER_N = 673_660, 781_135
TRI = st.TRI


# ------------------------------------------------------------------------------------------------------------- the runs


def matrix_of(row, kind: str) -> dict[str, dict[str, int]]:
    """The 5x6 matrix of a row of `runs.tsv` or `runs_v38.tsv`, with the no-calls of the reference-agreement records filled in."""
    m = {t: {p: int(row[f"cm_{t}__{p}"]) for p in CALLS} for t in CLASSES}
    return st.fill_nocall(m, kind)


def fires_of(row) -> dict[str, int]:
    return {c.removeprefix("fire_"): int(row[c]) for c in row.index if c.startswith("fire_")}


def read_dir(path: Path, kind: str) -> tuple[dict, dict]:
    """The matrix and the firing counts of a run directory (the main runs, the October base runs and the BS2 sweeps of v37)."""
    m = st.read_matrix(path / "concordance_matrix.csv", kind)
    fire = {}
    with (path / "criterion_firing_rates.csv").open() as f:
        for row in csv.DictReader(f):
            fire[row["criterion"]] = sum(int(row[f"{c}_fired"]) for c in CLASSES)
    return m, fire


def opposite_of(path: Path) -> pd.DataFrame:
    """The opposite-direction calls of a run, from its discrepancies file, which holds every one of them."""
    d = pd.read_csv(path / "discrepancies.tsv", sep="\t", dtype=str, keep_default_na=False)
    path_side, ben_side = ("Pathogenic", "Likely_pathogenic"), ("Likely_benign", "Benign")
    keep = (d.clinvar.isin(path_side) & d.predicted.isin(ben_side)) | (d.clinvar.isin(ben_side) & d.predicted.isin(path_side))
    return d[keep].reset_index(drop=True)


class Pair:
    """A run of the family v37 and its twin of the family v38."""

    def __init__(self, v37: str, v38: str, group: str, knob: str, value: str, kind: str, dir37: Path, dir38: Path, m37, m38, f37, f38):
        self.v37, self.v38, self.group, self.knob, self.value, self.kind = v37, v38, group, knob, value, kind
        self.dir37, self.dir38 = dir37, dir38
        self.m37, self.m38, self.f37, self.f38 = m37, m38, f37, f38
        self.x37, self.x38 = st.matrix_metrics(m37), st.matrix_metrics(m38)


def load_pairs() -> tuple[list[Pair], Pair, Pair]:
    """Every pair of the two families, and the two baselines (the main run and the run on the October release)."""
    apr37 = pd.read_csv(RUNS_DIR / "runs.tsv", sep="\t").set_index("name")
    apr38 = pd.read_csv(RUNS_DIR / "runs_v38.tsv", sep="\t").set_index("name")
    oct37 = pd.read_csv(RUNS_DIR / "holdout_runs.tsv", sep="\t").set_index("name")
    oct38 = pd.read_csv(RUNS_DIR / "holdout_runs_v38.tsv", sep="\t").set_index("name")
    runs = BENCH / "manuscript_runs"
    pairs: list[Pair] = []
    missing = []
    for name37, row in apr37.iterrows():
        name38 = "v38" + name37[3:]
        if name38 not in apr38.index:
            missing.append(name38)
            continue
        r38 = apr38.loc[name38]
        pairs.append(Pair(name37, name38, row["group"], row["knob"], str(row["value"]), "april", runs / name37, runs / name38,
                          matrix_of(row, "april"), matrix_of(r38, "april"), fires_of(row), fires_of(r38)))
    for name37, row in oct37.iterrows():
        name38 = "v38" + name37[3:]
        if name38 not in oct38.index:
            missing.append(name38)
            continue
        r38 = oct38.loc[name38]
        pairs.append(Pair(name37, name38, "holdout", row["knob"], str(row["value"]), "october", runs / name37, runs / name38,
                          matrix_of(row, "october"), matrix_of(r38, "october"), fires_of(row), fires_of(r38)))
    if missing:
        sys.exit(f"{len(missing)} runs of the family v38 are not in the files of the runs: {missing[:6]}")
    # The BS2 sweeps of the carrier bar and floor of v37 are runs outside the family (`output_v37_bs2*`), and on v38 they are the runs of Section 3.8 with the same grids.
    earlier = [("v37_bs2off", "v38_sw_bs2ad_0", "bs2_ad_prevalence_threshold", "0"), ("v37_bs2prev_3e-5", "v38_sw_bs2ad_3e-05", "bs2_ad_prevalence_threshold", "3e-05"),
               ("v37_bs2prev_1e-4", "v38_sw_bs2ad_0.0001", "bs2_ad_prevalence_threshold", "0.0001"),
               ("v37_bs2prev_1e-3", "v38_sw_bs2ad_0.001", "bs2_ad_prevalence_threshold", "0.001"),
               ("v37_bs2prev_3e-3", "v38_sw_bs2ad_0.003", "bs2_ad_prevalence_threshold", "0.003"),
               ("v37_bs2ac_20", "v38_sw_bs2ac_20", "bs2_ad_min_ac (bar off)", "20"), ("v37_bs2ac_100", "v38_sw_bs2ac_100", "bs2_ad_min_ac (bar off)", "100"),
               ("v37_bs2ac_1000", "v38_sw_bs2ac_1000", "bs2_ad_min_ac (bar off)", "1000")]
    for name37, name38, knob, value in earlier:
        d37 = BENCH / f"output_{name37}"
        if not (d37 / "concordance_matrix.csv").exists():
            continue
        m37, f37 = read_dir(d37, "april")
        r38 = apr38.loc[name38]
        pairs.append(Pair(name37, name38, "sweep", knob, value, "april", d37, runs / name38, m37, matrix_of(r38, "april"), f37, fires_of(r38)))
    m37, f37 = read_dir(BENCH / "output_v37", "april")
    m38, f38 = read_dir(runs / "v38", "april")
    main = Pair("v37", "v38", "baseline", "", "shipped", "april", BENCH / "output_v37", runs / "v38", m37, m38, f37, f38)
    m37, f37 = read_dir(BENCH / "output_v37_oct2026", "october")
    m38, f38 = read_dir(runs / "v38_oct_shipped", "october")
    october = Pair("v37_oct2026", "v38_oct_shipped", "baseline", "", "shipped", "october", BENCH / "output_v37_oct2026", runs / "v38_oct_shipped", m37, m38, f37, f38)
    return pairs, main, october


# ------------------------------------------------------------------------------------------------------ the tidy table


def rel(a: float, b: float) -> float | str:
    return (b - a) / abs(a) if a else ""


def tidy(pairs: list[Pair]) -> pd.DataFrame:
    rows = []
    for p in pairs:
        base = dict(group=p.group, run_v37=p.v37, run_v38=p.v38, knob=p.knob, value=p.value, set=p.kind)
        for k in METRICS:
            a, b = p.x37[k], p.x38[k]
            rows.append({**base, "kind": "metric", "item": k, "v37": a, "v38": b, "delta": b - a, "relative_delta": rel(a, b)})
        for t in CLASSES:
            for c in CALLS:
                a, b = p.m37[t][c], p.m38[t][c]
                rows.append({**base, "kind": "cell", "item": f"{t}__{c}", "v37": a, "v38": b, "delta": b - a, "relative_delta": rel(a, b)})
        for c in sorted(set(p.f37) | set(p.f38)):
            a, b = p.f37.get(c, 0), p.f38.get(c, 0)
            rows.append({**base, "kind": "criterion", "item": c, "v37": a, "v38": b, "delta": b - a, "relative_delta": rel(a, b)})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------------ sweeps


# (label, knob, the metric the sweep is meant to gain, the metric it costs, the shipped value, the group that holds the runs)
SWEEPS = [
    ("PM2 bar, not recessive-only", "pm2_ad_af_threshold", "same_pathogenic", "false_pathogenic", "4e-05", "sweep"),
    ("PM2 bar, recessive-only", "pm2_ar_af_threshold", "same_pathogenic", "false_pathogenic", "7e-05", "sweep"),
    ("BS1 bar", "bs1_af_threshold", "same_benign", "false_benign", "0.005", "sweep"),
    ("BS2 homozygote bar", "bs2_hom_prevalence_threshold", "same_benign", "false_benign", "0.001", "sweep"),
    ("BS2 carrier bar, dominant-only genes", "bs2_ad_prevalence_threshold", "same_benign", "false_benign", "0.0003", "sweep"),
    ("BS2 carrier floor, bar off", "bs2_ad_min_ac (bar off)", "same_benign", "false_benign", "5", "sweep"),
    ("BP7 far boundary", "bp7_max_intron_offset", "same_benign", "false_benign", "300", "sweep"),
    ("BP7 PhyloP bar", "phylop_conserved", "same_benign", "false_benign", "2.0", "sweep"),
    ("BP1 veto count", "bp1_max_pathogenic_missense", "same_benign", "false_benign", "3", "sweep"),
    ("PM2 bar, not recessive-only, 2015 table both sides", "pm2_ad_af_threshold (Richards 2015 table both sides)", "same_pathogenic", "false_pathogenic", "4e-05", "interaction"),
    ("Likely pathogenic cutoff (points)", "likely_pathogenic_points", "same_pathogenic", "false_pathogenic", "6", "cutoff"),
    ("Pathogenic cutoff (points)", "pathogenic_points", "same_pathogenic", "false_pathogenic", "10", "cutoff"),
]


def sweep_points(pairs: list[Pair], main: Pair, knob: str, shipped: str, group: str) -> list[dict]:
    """The settings of a sweep, the default among them, each with its metrics under both families, ordered by value.

    The default setting is the main run, except in the two sweeps whose curve is not made under the default combination of settings:
    the PM2 sweep under the 2015 table, whose point at 4e-5 is the run with the table and nothing else changed
    (`<family>_sw_table_both_sides`), and the floor of BS2's carrier arm with the bar off, whose point at the default floor of 5 is
    the run with the bar off (`v37_bs2off`, `v38_sw_bs2ad_0`).
    """
    by_name = {p.v37: p for p in pairs}
    out = []
    for p in pairs:
        if p.knob == knob and p.group in (group, "sweep"):
            # the run at 6 points is the control of the cutoff build, which must reproduce the default: the default row stands for it
            if knob == "likely_pathogenic_points" and p.value == "6":
                continue
            out.append({"setting": float(p.value), "value": p.value, "run_v37": p.v37, "run_v38": p.v38, "x37": p.x37, "x38": p.x38, "shipped": False})
    if knob.startswith("pm2_ad_af_threshold (Richards"):
        t = by_name["v37_sw_table_both_sides"]
        out.append({"setting": float(shipped), "value": shipped, "run_v37": t.v37, "run_v38": t.v38, "x37": t.x37, "x38": t.x38, "shipped": True})
    elif knob.startswith("bs2_ad_min_ac"):
        t = by_name["v37_bs2off"]
        out.append({"setting": float(shipped), "value": shipped, "run_v37": t.v37, "run_v38": t.v38, "x37": t.x37, "x38": t.x38, "shipped": True})
    else:
        out.append({"setting": float(shipped), "value": shipped, "run_v37": "v37", "run_v38": "v38", "x37": main.x37, "x38": main.x38, "shipped": True})
    return sorted(out, key=lambda o: o["setting"])


def steps_of(points: list[dict], gain: str, cost: str, which: str) -> list[dict]:
    """The steps of a sweep ordered by what it gains, as `13_manuscript_stats.sweep_steps` orders them: from the setting with the least
    to the one with the most, each pair of neighbours a step with the calls it gains, the errors it adds and the calls per error."""
    order = sorted(points, key=lambda o: (o[which][gain], o["setting"]))
    steps = []
    for a, b in zip(order[:-1], order[1:]):
        dg, dc = b[which][gain] - a[which][gain], b[which][cost] - a[which][cost]
        steps.append({"from": a["setting"], "to": b["setting"], "gain": dg, "cost": dc, "rate": (dg / dc) if dc > 0 else math.inf})
    return steps


KNEE_FACTOR = 3.0  # the fall that makes a knee: a step worth a third or less of the step before it
KNEE_FROM = 150.0  # the step before it is worth hundreds of calls per added error (Section 2.5 of the manuscript)
KNEE_TO = 150.0  # and the step after it about a hundred or fewer (the same Section 2.5); a fall by a factor of 3 to 480 calls per error is not a knee


def largest_fall(steps: list[dict]) -> tuple[float, float, float, float | None]:
    """The two neighbouring steps, in the order of what they gain, between which the exchange rate falls the most, as (the factor, the rate
    of the first, the rate of the second, the setting that they share); the rates of a step that adds no error are infinite and are not compared."""
    best: tuple[float, float, float, float | None] = (1.0, 0.0, 0.0, None)
    for a, b in zip(steps[:-1], steps[1:]):
        if math.isinf(a["rate"]) or math.isinf(b["rate"]) or b["rate"] <= 0:
            continue
        f = a["rate"] / b["rate"]
        if f > best[0]:
            best = (f, a["rate"], b["rate"], a["to"] if a["to"] == b["from"] else b["from"])
    return best


def knee_of(steps: list[dict]) -> tuple[float | None, float]:
    """Where the exchange rate falls the most from one step to the next, along the steps in the order of what they gain, if that is a knee.

    Section 2.5 of the manuscript takes the knee of a curve to be the point past which the exchange rate falls from hundreds of calls per
    added error to about a hundred or fewer. Here it is the setting that the two neighbouring steps of the largest fall share, when the fall
    is by a factor of 3 or more, the first step is worth 150 calls per error or more and the second 150 or fewer; otherwise there is no
    knee. The rule gives the knees that the manuscript names, BS1 at 1e-3 and the far boundary of BP7 at 300, under both versions and on both
    sets of databases, and none for the other sweeps that the manuscript says have none. (The first form of the rule asked only for the
    factor of 3 and the first step; on the complete databases the carrier floor of BS2 with the bar off then fell by a factor of 3.3, from 1,578
    calls per error to 480, which is not the fall to about a hundred that Section 2.5 describes, and the rule now asks for it.)
    Returns the setting and the factor of the largest fall.
    """
    factor, before, after, where = largest_fall(steps)
    if factor < KNEE_FACTOR or before < KNEE_FROM or after > KNEE_TO:
        return None, factor
    return where, factor


def sweep_shape(points: list[dict], gain: str, cost: str) -> dict:
    """What a sweep looks like under each family: its ranges, the sign of each step along the settings and its knee."""
    res = {}
    for tag, which in (("v37", "x37"), ("v38", "x38")):
        ordered = sorted(points, key=lambda o: o["setting"])
        d_gain = [b[which][gain] - a[which][gain] for a, b in zip(ordered[:-1], ordered[1:])]
        d_cost = [b[which][cost] - a[which][cost] for a, b in zip(ordered[:-1], ordered[1:])]
        sign = lambda x: (x > 0) - (x < 0)  # noqa: E731
        steps = steps_of(points, gain, cost, which)
        knee, factor = knee_of(steps)
        _, fall_from, fall_to, _ = largest_fall(steps)
        res[tag] = {
            "range_same_direction": max(o[which]["same_direction"] for o in points) - min(o[which]["same_direction"] for o in points),
            "range_gain": max(o[which][gain] for o in points) - min(o[which][gain] for o in points),
            "range_false_benign": max(o[which]["false_benign"] for o in points) - min(o[which]["false_benign"] for o in points),
            "range_false_pathogenic": max(o[which]["false_pathogenic"] for o in points) - min(o[which]["false_pathogenic"] for o in points),
            "gain_signs": "".join("+" if s > 0 else "-" if s < 0 else "0" for s in map(sign, d_gain)),
            "cost_signs": "".join("+" if s > 0 else "-" if s < 0 else "0" for s in map(sign, d_cost)),
            "knee": knee, "knee_factor": factor, "fall_from_rate": fall_from, "fall_to_rate": fall_to,
            "order_by_same_direction": [o["setting"] for o in sorted(points, key=lambda o: (o[which]["same_direction"], o["setting"]))],
            "max_rate": max([s["rate"] for s in steps if not math.isinf(s["rate"])] or [0]),
            "min_rate": min([s["rate"] for s in steps if not math.isinf(s["rate"])] or [0]),
        }
    # settings whose order by agreement flipped between the families
    flips = []
    for i, a in enumerate(points):
        for b in points[i + 1:]:
            d37 = (a["x37"]["same_direction"] > b["x37"]["same_direction"]) - (a["x37"]["same_direction"] < b["x37"]["same_direction"])
            d38 = (a["x38"]["same_direction"] > b["x38"]["same_direction"]) - (a["x38"]["same_direction"] < b["x38"]["same_direction"])
            if d37 * d38 < 0:
                flips.append((a["value"], b["value"]))
    res["flips"] = flips
    return res


def sweep_table(pairs: list[Pair], main: Pair) -> tuple[pd.DataFrame, list[dict]]:
    rows, summary = [], []
    for label, knob, gain, cost, shipped, group in SWEEPS:
        pts = sweep_points(pairs, main, knob, shipped, group)
        if len(pts) < 3:
            continue
        shape = sweep_shape(pts, gain, cost)
        for o in pts:
            r = {"sweep": label, "knob": knob, "setting": o["value"], "default": o["shipped"], "run_v37": o["run_v37"], "run_v38": o["run_v38"]}
            for k in ("exact", "same_direction", "same_pathogenic", "same_benign", "false_benign", "false_pathogenic", "opposite", "plp_recall", "blb_recall"):
                r[f"{k}_v37"], r[f"{k}_v38"] = o["x37"][k], o["x38"][k]
                r[f"{k}_delta"] = o["x38"][k] - o["x37"][k]
            rows.append(r)
        summary.append({"label": label, "knob": knob, "gain": gain, "cost": cost, "shape": shape, "n": len(pts), "points": pts})
    return pd.DataFrame(rows), summary


# ------------------------------------------------------------------------------------------------------ effects


EFFECT_METRICS = ["same_direction", "same_pathogenic", "same_benign", "false_benign", "false_pathogenic", "opposite", "exact"]


def effects(pairs: list[Pair], main: Pair) -> pd.DataFrame:
    """The effect of every switch, source, hook and gap run, and of every other run on the April set, under both families."""
    rows = []
    for p in pairs:
        if p.kind != "april" or p.group in ("baseline",):
            continue
        for k in EFFECT_METRICS:
            e37, e38 = p.x37[k] - main.x37[k], p.x38[k] - main.x38[k]
            changed_sign = e37 * e38 < 0
            vanished_or_appeared = (e37 == 0) != (e38 == 0)
            rel_ = abs(e38 - e37) / abs(e37) if e37 else math.inf if e38 else 0.0
            rows.append({"group": p.group, "run_v37": p.v37, "run_v38": p.v38, "knob": p.knob, "value": p.value, "metric": k,
                         "effect_v37": e37, "effect_v38": e38, "effect_delta": e38 - e37,
                         "changed_sign": bool(changed_sign), "appeared_or_vanished": bool(vanished_or_appeared),
                         "relative_change": "" if math.isinf(rel_) else rel_, "changed_by_more_than_20_percent": bool(rel_ > 0.2),
                         "small_number": bool(max(abs(e37), abs(e38)) < 10)})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------------ opposite calls


def describe(d: pd.Series) -> str:
    return f"{d.gene} {d.chrom}:{d.pos} {d.ref}>{d.alt}"


def opposite_changes(pairs: list[Pair], main: Pair, october: Pair) -> tuple[pd.DataFrame, dict]:
    """For every pair, the opposite-direction calls that one run has and the other lacks."""
    rows, counts = [], {}
    for p in [main, october] + pairs:
        try:
            a, b = opposite_of(p.dir37), opposite_of(p.dir38)
        except FileNotFoundError:
            continue
        ka = {tuple(r) for r in a[KEY].itertuples(index=False)}
        kb = {tuple(r) for r in b[KEY].itertuples(index=False)}
        counts[p.v37] = (len(a), len(b), len(ka - kb), len(kb - ka))
        for which, other, keys, src, oth in (("left_in_v38", "v38", ka - kb, a, b), ("entered_in_v38", "v37", kb - ka, b, a)):
            for key in sorted(keys):
                m = src[(src.chrom == key[0]) & (src.pos == key[1]) & (src.ref == key[2]) & (src.alt == key[3])].iloc[0]
                rows.append({"pair_v37": p.v37, "pair_v38": p.v38, "change": which, "gene": m.gene, "chrom": m.chrom, "pos": m.pos, "ref": m.ref, "alt": m.alt,
                             "stars": m.stars, "clinvar": m.clinvar,
                             "call_in_the_run_that_has_it": m.predicted, "criteria_in_the_run_that_has_it": m.met_criteria,
                             "consequence": m.consequence})
    return pd.DataFrame(rows), counts


# ------------------------------------------------------------------------------------------------------ the text


class Doc:
    """Markdown with one sentence on a line."""

    def __init__(self) -> None:
        self.out: list[str] = []

    def h(self, level: int, text: str) -> None:
        self.blank()
        self.out += ["#" * level + " " + text, ""]

    def s(self, text: str) -> None:
        self.out.append(text)

    def blank(self) -> None:
        if self.out and self.out[-1] != "":
            self.out.append("")

    def bullets(self, items: list[str]) -> None:
        self.blank()
        self.out += [f"- {t}" for t in items]
        self.blank()

    def table(self, header: list[str], rows: list[list[str]], align: str | None = None) -> None:
        """A table; `align` is one letter for each column, `l` or `r` (the first column on the left and the others on the right by default)."""
        align = align or "l" + "r" * (len(header) - 1)
        self.blank()
        self.out.append("| " + " | ".join(header) + " |")
        self.out.append("|" + "|".join(":---" if a == "l" else "---:" for a in align) + "|")
        self.out += ["| " + " | ".join(r) + " |" for r in rows]
        self.blank()

    def text(self) -> str:
        return "\n".join(self.out).rstrip("\n") + "\n"


def n(x) -> str:
    """A count with thousands separators, rounded half up."""
    return f"{int(float(x) + 0.5):,}" if float(x) >= 0 else f"-{int(abs(float(x)) + 0.5):,}"


def sg(x) -> str:
    """A signed count."""
    v = int(round(float(x)))
    return f"{v:+,}" if v else "0"


def pc(x, d: int = 2) -> str:
    return f"{float(x) * 100:.{d}f}%"


def pa(x, d: int = 2) -> str:
    """A difference of two proportions in percentage points, without its sign (for a sentence that says above or below)."""
    return f"{abs(float(x)) * 100:.{d}f}"


WORDS = ["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve"]


def word(k) -> str:
    """A small count as a word."""
    k = int(k)
    return WORDS[k] if 0 <= k < len(WORDS) else str(k)


def plural(k: int, one: str, many: str) -> str:
    return f"{k} {one if k == 1 else many}"


def setting(v: str, knob: str = "") -> str:
    """A setting as the manuscript writes it: 1e-6 and not 1e-06 or 0.000001, and what a zero means for the knob."""
    x = float(v)
    if 0 < x < 0.1:
        text = f"{x:.0e}".replace("e-0", "e-")
    else:
        text = v
    if x == 0:
        if knob.startswith("pm2_"):
            return "0 (strict absence)"
        if knob.startswith("bs2_ad_prevalence") or knob.startswith("bp7_max_intron"):
            return "0 (off)"
    return text


def short(text: str, width: int = 70) -> str:
    """An anchor cut at a word."""
    if len(text) <= width:
        return text
    return text[:width].rsplit(" ", 1)[0] + " ..."


def pp(x, d: int = 2) -> str:
    """A difference of two proportions in percentage points, signed."""
    v = float(x) * 100
    t = f"{abs(v):.{d}f}"
    return t if float(t) == 0 else ("+" if v > 0 else "-") + t


def prov(name: str) -> dict:
    f = RUNS_DIR / "provenance" / f"{name}.txt"
    return dict(line.split(": ", 1) for line in f.read_text().splitlines() if ": " in line)


def manuscript_lines() -> list[str]:
    return (ROOT / "manuscript/acmg/manuscript.md").read_text().splitlines()


def find_line(lines: list[str], anchor: str) -> int:
    """The number of the one line of the manuscript that holds `anchor`."""
    hits = [i + 1 for i, line in enumerate(lines) if anchor in line]
    if len(hits) != 1:
        raise SystemExit(f"the anchor {anchor!r} is on {len(hits)} lines of the manuscript, not one")
    return hits[0]


# ---------------------------------------------------------------------------------------------- 1. what is compared


def new_york_days(times: list[str]) -> str:
    """The New York date of times that the provenance files record in UTC, which is how the paper dates a run: `8 November 2026`, or `8 to 9 November 2026`."""
    days = sorted({dt.datetime.fromisoformat(t.replace("Z", "+00:00")).astimezone(ZoneInfo("America/New_York")).date() for t in times})
    if len(days) == 1:
        return f"{days[0].day} {days[0].strftime('%B %Y')}"
    return f"{days[0].day} to {days[-1].day} {days[-1].strftime('%B %Y')}"


def provenance_counts(names: list[str]) -> dict[str, int]:
    out: dict[str, int] = {}
    for nm in names:
        k = prov(nm)["fastvep_sha256"][:8]
        out[k] = out.get(k, 0) + 1
    return out


def section_compared(D: Doc, pairs: list[Pair], main: Pair, october: Pair) -> None:
    fam37 = [p.v37 for p in pairs if p.dir37.parent == BENCH / "manuscript_runs"]
    fam38 = [p.v38 for p in pairs if p.dir37.parent == BENCH / "manuscript_runs"]
    shared = [x for x in fam38 if is_shared(x)]
    new38 = fam38  # every run of the family was made once, with the complete databases
    c38 = provenance_counts(new38)
    rec = sorted(prov(x)["recorded_at"] for x in new38)
    secs = sorted(int(prov(x)["annotation_seconds"]) for x in new38)
    groups: dict[str, int] = {}
    for p in pairs:
        g = p.group if p.dir37.parent == BENCH / "manuscript_runs" else "earlier BS2 runs"
        groups[g] = groups.get(g, 0) + 1
    same_def = definitions_agree(pairs)
    D.h(2, "1. What is compared")
    D.s(f"The family v37 is the {sum(1 for p in pairs if p.kind == 'april' and p.dir37.parent == BENCH / 'manuscript_runs')} runs on the April set and the {sum(1 for p in pairs if p.kind == 'october')} runs on the October release that `11_manuscript_runs.py` defines, made on {new_york_days([prov(x)['recorded_at'] for x in fam37])} by binaries built from commit 5d8004a (`binaries.tsv`).")
    D.s(f"The family v38 is the same {len(fam37)} runs made by the classifier of release 0.5.0, the tag `v0.5.0-rc.1`, with the same databases.")
    D.s(f"{len(shared)} of them are runs that Section 3.8 also defines under the same names (the sweeps `v38_sw_pm2ad_*`, `v38_sw_pm2ar_*` and `v38_sw_bs2hom_*`, whose grids are the family's), and the harness makes each of them once.")
    D.s(f"The {len(new38)} were made on {new_york_days(rec)} (New York; the provenance records UTC, {rec[0][:10]}) by binaries built from the tag, {len([x for x in new38 if x.startswith('v38_oct_')])} of them on the October release, several at a time on the ten-core machine (the annotation of one run took {secs[0]} to {secs[-1]} seconds, median {secs[len(secs) // 2]}).")
    D.s(f"The binaries were `fastvep_v38` (SHA-256 `70b7c684...`, {c38.get('70b7c684', 0)} runs), `fastvep_cut38` (`de44db94...`, {c38.get('de44db94', 0)} runs, the point cutoffs) and `fastvep_hook38` (`9886d1b5...`, {c38.get('9886d1b5', 0)} runs, BP1's own-record count and BP7's SpliceAI threshold), each built from a clean export of the tag in a target directory of its own by `build_manuscript_binaries.sh`, and `binaries_v38.tsv` records them.")
    D.s(f"The two families are generated by the same functions of `11_manuscript_runs.py` with the name prefix and the binaries as the only parameters, and {same_def[0]} of the {same_def[1]} pairs have the same knob, value, default, removed and replaced annotation files and environment in `run_matrix.tsv` and `run_matrix_v38.tsv`, so a pair differs in the classifier and in nothing else.")
    earlier = [p for p in pairs if p.dir37.parent == BENCH]
    D.s(f"{len(earlier)} earlier runs of the BS2 carrier settings (`output_v37_bs2off`, `output_v37_bs2prev_*` and `output_v37_bs2ac_*`, made on 8 November 2026 by the binary of v37, like the paper's runs) are paired with the runs of Section 3.8 that have the same grids (`v38_sw_bs2ad_*` and `v38_sw_bs2ac_*`), which makes {len(pairs)} pairs besides the two baselines in the CSV files; the table below says how many pairs each group holds.")
    D.table(["Group", "Pairs"], [[g, str(k)] for g, k in sorted(groups.items())])
    D.s("The main run `v37` is paired with `v38`, and the run of `14_temporal_holdout.py` on the October release (`v37_oct2026`) with `v38_oct_shipped`.")
    D.s("Everything is over all 673,660 records of the April set and all 781,135 of the October release, as the tables of the paper are.")
    D.s("The 382 ClinVar reference-agreement records of the April set (385 in the October subset) have `.` as their alternative allele and get no call under either classifier, but the binaries of v37 write a record for them and the binaries of v38 do not, so the matrices of v38 hold 673,278 records and the scoring script counts the other 382 apart.")
    D.s("The matrices of v38 here are filled with those records as no-calls, as `13_manuscript_stats.py` does for its twin tables (`fill_nocall`), so that both families have the same denominators; the benign-side recall of v38 is therefore over 284,550 variants, as v37's is, and not over 284,168, as in the tables of Section 3.8.")
    D.s("Three kinds of statement are used below.")
    D.bullets([
        "A change is a difference of two counts of the same pair, v38 less v37, and a relative change is that difference over the count under v37.",
        f"The effect of an ablation is the run less the default under the same classifier; it changed sign if the two effects have opposite signs, and by more than 20% if the effect under v38 differs from the one under v37 by more than a fifth of the latter; an effect of fewer than 10 calls under both is called small, because 20% of it is less than two calls.",
        f"The knee of a sweep is where the exchange rate falls the most from one step to the next (Section 2.5 of the manuscript: from hundreds of calls per added error to about a hundred or fewer), taken here as a fall by a factor of {KNEE_FACTOR:g} or more from a step worth {KNEE_FROM:g} calls per error or more; under both classifiers the rule finds the two knees that the manuscript names (BS1 at 1e-3 and BP7's far boundary at 300, under both versions), none in the sweeps that the manuscript says have none (BS2's homozygote bar, BP7's PhyloP bar and the veto count of BP1), and two that it does not state in these terms (the Likely pathogenic cutoff at 5 points, and the PM2 threshold under the 2015 table at 1e-5, a point inside the span from strict absence to 4e-5 whose rise the manuscript says that v10 took for a knee).",
    ])
    D.s("Every number in this file is read from the scored files of the runs and from the twin tables of `manuscript/acmg/data/`, and `20_compare_families.py` writes this file and the five CSV files beside it.")
    D.bullets([
        "`v38_family_vs_v37.csv`: one row for every pair, number and cell (the 30 cells of the concordance matrix, 16 headline metrics and the firing count of every criterion), under both families, with the difference and the relative difference.",
        "`v38_family_vs_v37_sweeps.csv`: every sweep as a curve, one row for every setting.",
        "`v38_family_vs_v37_effects.csv`: every run of the April set as an effect against the default, under both families, with the flags of the second bullet above.",
        "`v38_family_vs_v37_opposite.csv`: every opposite-direction call that one run of a pair has and the other lacks.",
    ])

SHARED_PREFIXES = ("v38_sw_pm2ad_", "v38_sw_pm2ar_", "v38_sw_bs2hom_")


def is_shared(name: str) -> bool:
    """A run of the family v38 that Section 3.8 made under the same name."""
    return name.startswith(SHARED_PREFIXES)


def definitions_agree(pairs: list[Pair]) -> tuple[int, int]:
    """How many pairs have the same setting in the two run matrices (knob, value, default, annotation files and environment)."""
    m37 = pd.read_csv(RUNS_DIR / "run_matrix.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("name")
    m38 = pd.read_csv(RUNS_DIR / "run_matrix_v38.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("name")
    cols = ["knob", "value", "shipped", "removed_from_sa_db", "replaced_in_sa_db", "environment"]
    same = total = 0
    for p in pairs:
        if p.v37 not in m37.index:
            continue  # the earlier runs of the BS2 settings are not in the run matrix
        total += 1
        same += int(all(m37.loc[p.v37, c] == m38.loc[p.v38, c] for c in cols))
    return same, total


# ---------------------------------------------------------------------------------------------- 2. the main run


def three_way_matrix() -> dict:
    """v38_base, the commit before the change of Section 3.8, as a matrix of the same denominators."""
    r = pd.read_csv(RUNS_DIR / "runs_v38.tsv", sep="\t").set_index("name")
    return matrix_of(r.loc["v38_base"], "april")


def per_variant_comparison() -> dict:
    """The per-variant tables of v37, of the commit before the change and of v38 side by side."""
    key = KEY
    v37 = st.read_pv(BENCH / "manuscript_runs/v37_per_variant.tsv.gz")
    v38 = st.read_pv(BENCH / "manuscript_runs/v38_per_variant.tsv.gz")
    base = st.read_pv(BENCH / "manuscript_runs/v38_base_per_variant.tsv.gz")
    cols38 = key + ["call", "criteria", "consequence", "hgvsc", "hgvsp", "clingen_tags", "gnomad_ac", "gnomad_an", "gnomad_hc", "gnomad_popmax_af", "symbol"]
    x = v37.merge(v38[cols38], on=key, suffixes=("", "_38"), validate="one_to_one")
    x = x.merge(base[key + ["call", "criteria"]].rename(columns={"call": "call_base", "criteria": "criteria_base"}), on=key, validate="one_to_one")
    x["cv3"] = x["clinvar"].map(TRI)
    return {"x": x, "n": len(x)}


def code_of(c: str) -> str:
    return c.split("_")[0]


def section_main(D: Doc, main: Pair, cmp_: dict) -> None:
    a, b = main.x37, main.x38
    base_m = st.matrix_metrics(three_way_matrix())
    x = cmp_["x"]
    D.h(2, "2. The main run")
    D.h(3, "2.1 The headline numbers under both classifiers")
    N = APRIL_N
    rows = [
        ["Exact agreement", n(a["exact"]), n(b["exact"]), sg(b["exact"] - a["exact"]), f"{pc(a['exact'] / N)} to {pc(b['exact'] / N)}"],
        ["Same-direction agreement", n(a["same_direction"]), n(b["same_direction"]), sg(b["same_direction"] - a["same_direction"]), f"{pc(a['same_direction'] / N)} to {pc(b['same_direction'] / N)}"],
        ["ClinVar P/LP called P/LP", n(a["same_pathogenic"]), n(b["same_pathogenic"]), sg(b["same_pathogenic"] - a["same_pathogenic"]), f"recall {pc(a['plp_recall'])} to {pc(b['plp_recall'])}"],
        ["ClinVar B/LB called LB/B", n(a["same_benign"]), n(b["same_benign"]), sg(b["same_benign"] - a["same_benign"]), f"recall {pc(a['blb_recall'])} to {pc(b['blb_recall'])}"],
        ["ClinVar VUS called VUS", n(a["vus_exact"]), n(b["vus_exact"]), sg(b["vus_exact"] - a["vus_exact"]), f"recall {pc(a['vus_recall'])} to {pc(b['vus_recall'])}"],
        ["ClinVar P/LP called VUS", n(a["plp_to_vus"]), n(b["plp_to_vus"]), sg(b["plp_to_vus"] - a["plp_to_vus"]), ""],
        ["ClinVar B/LB called VUS", n(a["blb_to_vus"]), n(b["blb_to_vus"]), sg(b["blb_to_vus"] - a["blb_to_vus"]), ""],
        ["ClinVar VUS called P/LP", n(a["vus_to_plp"]), n(b["vus_to_plp"]), sg(b["vus_to_plp"] - a["vus_to_plp"]), ""],
        ["ClinVar VUS called LB/B", n(a["vus_to_blb"]), n(b["vus_to_blb"]), sg(b["vus_to_blb"] - a["vus_to_blb"]), ""],
        ["Reference-discordant benign calls", n(a["false_benign"]), n(b["false_benign"]), sg(b["false_benign"] - a["false_benign"]), ""],
        ["Reference-discordant pathogenic calls", n(a["false_pathogenic"]), n(b["false_pathogenic"]), sg(b["false_pathogenic"] - a["false_pathogenic"]), ""],
        ["Opposite-direction calls", n(a["opposite"]), n(b["opposite"]), sg(b["opposite"] - a["opposite"]), f"{pc(a['opposite'] / N, 4)} to {pc(b['opposite'] / N, 4)} of the records"],
        ["No call", n(a["no_call"]), n(b["no_call"]), sg(b["no_call"] - a["no_call"]), ""],
    ]
    D.table(["Measure", "v37", "v38", "Change", "As a rate"], rows)
    D.s(f"v38 agrees with ClinVar's class for {n(b['exact'] - a['exact'])} more of the {n(N)} records (exact agreement {pc(a['exact'] / N)} to {pc(b['exact'] / N)}) and with its direction for {n(b['same_direction'] - a['same_direction'])} more ({pc(a['same_direction'] / N)} to {pc(b['same_direction'] / N)}).")
    D.s(f"The gain is on the benign side: {n(b['same_benign'] - a['same_benign'])} more ClinVar benign-side variants are called benign-side (recall {pc(a['blb_recall'])} to {pc(b['blb_recall'])}), while {n(abs(b['same_pathogenic'] - a['same_pathogenic']))} fewer ClinVar pathogenic-side variants are called pathogenic-side (recall {pc(a['plp_recall'])} to {pc(b['plp_recall'])}).")
    D.s(f"Opposite-direction calls go from {a['opposite']} to {b['opposite']} ({a['false_benign']} reference-discordant benign calls and {a['false_pathogenic']} reference-discordant pathogenic calls under v37, {b['false_benign']} and {b['false_pathogenic']} under v38).")
    benign_rows = [t for t in CLASSES if t in ("Likely_benign", "Benign")]
    benign_all = sum(sum(main.m38[t].values()) for t in benign_rows)
    benign_called = benign_all - sum(main.m38[t]["NoCall"] for t in benign_rows)
    D.s(f"Recall of ClinVar's benign-side variants is {pc(b['blb_recall'])} under v38 over the same {n(benign_all)} variants as under v37, and {pc(b['same_benign'] / benign_called)} over the {n(benign_called)} that have a call (the denominator of Table 8 of the manuscript).")

    D.h(3, "2.2 Where the difference comes from: v37, the commit before the change, and v38")
    D.s("The commit before the change of Section 3.8 (`v38_base`) has the classification code of v37 and a newer annotation layer, and v38 is the change of Section 3.8 on top of it, so the difference between v37 and v38 is the sum of two steps.")
    rows = []
    for label, m in (("v37", a), ("v38_base, the commit before the change", base_m), ("v38", b)):
        rows.append([label, n(m["exact"]), n(m["same_direction"]), n(m["same_pathogenic"]), n(m["same_benign"]), str(m["false_benign"]), str(m["false_pathogenic"])])
    D.table(["Run", "Exact", "Same direction", "P/LP called P/LP", "B/LB called LB/B", "Ref.-discordant benign", "Ref.-discordant pathogenic"], rows)
    D.s(f"Of the {n(b['same_direction'] - a['same_direction'])} more calls in direction of v38, {n(base_m['same_direction'] - a['same_direction'])} are the commit before the change and {n(b['same_direction'] - base_m['same_direction'])} are the change itself.")
    D.s(f"Of the {n(b['exact'] - a['exact'])} more exact calls, {n(base_m['exact'] - a['exact'])} are the commit before the change and {n(b['exact'] - base_m['exact'])} the change.")
    ka = {tuple(r) for r in opposite_of(main.dir37)[KEY].itertuples(index=False)}
    kb = {tuple(r) for r in opposite_of(BENCH / "manuscript_runs/v38_base")[KEY].itertuples(index=False)}
    kc = {tuple(r) for r in opposite_of(main.dir38)[KEY].itertuples(index=False)}
    D.s(f"The commit before the change has {len(kb)} opposite-direction calls ({len(ka - kb)} of v37's are missing from them and {len(kb - ka)} are new), and v38 has {len(kc)} ({len(kb - kc)} of the commit's are missing from them and {len(kc - kb)} are new; Section 3 below).")

    D.h(3, "2.3 What changed call by call")
    ch = x[x.call != x.call_38]
    chb = x[x.call != x.call_base]
    D.s(f"Of the {n(cmp_['n'])} variants that both runs annotate, v38 gives another call to {n(len(ch))} ({n(len(chb))} of them already in the commit before the change) and other criteria to {n((x.criteria != x.criteria_38).sum())}.")
    tr = ch.groupby(["cv3", "call", "call_38"]).size().sort_values(ascending=False)
    D.s(f"The calls that move, by ClinVar direction and by the call before and after (the largest 14 of the {len(tr)} transitions):")
    D.table(["ClinVar", "v37 call", "v38 call", "Variants"], [[c, a_, b_, n(k)] for (c, a_, b_), k in tr.head(14).items()])
    total_tr = int(tr.sum())
    D.s(f"The 14 largest transitions hold {n(tr.head(14).sum())} of the {n(total_tr)} calls that move ({pc(tr.head(14).sum() / total_tr, 1)}).")
    plp = ch[ch.cv3 == "P/LP"]
    D.s(f"{n(len(plp))} of the calls that move are ClinVar pathogenic-side variants: {n(((plp.call == 'Likely_pathogenic') & (plp.call_38 == 'Pathogenic')).sum())} go from Likely pathogenic to Pathogenic and {n(((plp.call_38 == 'VUS') & plp.call.isin(['Likely_pathogenic', 'Pathogenic'])).sum())} lose their pathogenic-side call to VUS ({', '.join(sorted(set(plp[(plp.call_38 == 'VUS') & plp.call.isin(['Likely_pathogenic', 'Pathogenic'])].symbol)))}, which are the four loss-of-function variants that Section 3.8 names).")
    gl = []
    for c in sorted({code_of(t) for s_ in x.criteria for t in s_.split("&") if t}):
        has = lambda s_, c=c: any(code_of(t) == c for t in s_.split("&") if t)  # noqa: E731
        was, now = x.criteria.map(has), x.criteria_38.map(has)
        gained, lost = int((~was & now).sum()), int((was & ~now).sum())
        if gained or lost:
            gl.append((c, gained, lost))
    gl.sort(key=lambda t: -(t[1] + t[2]))
    D.s("The criteria that a variant gains or loses (a variant counts once for a criterion, whatever its strength):")
    D.table(["Criterion", "Variants gaining it", "Variants losing it"], [[c, n(g), n(l_)] for c, g, l_ in gl])
    big = gl[0]
    D.s(f"{big[0]} is the criterion that changes for the most variants: {n(big[1])} gain it and {n(big[2])} lose it, and the next, {gl[1][0]}, changes for {n(gl[1][1] + gl[1][2])}.")

# ---------------------------------------------------------------------------------------------- 3. opposite-direction calls


def section_opposite(D: Doc, main: Pair, october: Pair, pairs: list[Pair], opp: pd.DataFrame, counts: dict, cmp_: dict) -> None:
    a, b = main.x37, main.x38
    x = cmp_["x"]
    D.h(2, "3. The opposite-direction calls")
    left = opp[(opp.pair_v37 == "v37") & (opp.change == "left_in_v38")]
    entered = opp[(opp.pair_v37 == "v37") & (opp.change == "entered_in_v38")]
    D.s(f"The main run has {a['opposite']} opposite-direction calls under v37 ({a['false_benign']} reference-discordant benign calls and {a['false_pathogenic']} reference-discordant pathogenic calls) and {b['opposite']} under v38 ({b['false_benign']} and {b['false_pathogenic']}); {len(left)} call left the set and {len(entered)} entered it.")
    rulings = pd.read_csv(DATA / "round6_rulings.tsv", sep="\t", dtype=str)
    for r in left.itertuples():
        m = x[(x.chrom == r.chrom) & (x.pos == r.pos) & (x.ref == r.ref) & (x.alt == r.alt)].iloc[0]
        rl = rulings[(rulings.chrom == r.chrom) & (rulings.pos == r.pos) & (rulings.ref == r.ref) & (rulings.alt == r.alt)]
        ruling = f"in round 6 {ruling_text(rl.iloc[0].ruling, rl.iloc[0].rule)}" if len(rl) else "the reviewer never saw it"
        D.s(f"The call that left is {m.gene} {m.chrom}:{m.pos} {m.ref}>{m.alt} ({m.consequence}, {m.hgvsc}), which ClinVar classifies as {m.clinvar.replace('_', ' ')} at {m.stars} stars.")
        same_base = (m.call, m.criteria) == (m.call_base, m.criteria_base)
        D.s(f"v37 called it {m.call.replace('_', ' ')} on {m.criteria.replace('&', ', ')}, the commit before the change {'gave the same call on the same criteria' if same_base else 'called it ' + m.call_base.replace('_', ' ') + ' on ' + m.criteria_base.replace('&', ', ')}, and v38 calls it {m.call_38.replace('_', ' ')} on {m.criteria_38.replace('&', ', ')}.")
        gained = [t for t in m.criteria_38.split("&") if t not in m.criteria.split("&")]
        D.s(f"The established ClinGen curations of the gene carry the tag {m.clingen_tags}, and gnomAD holds the variant in {n(m.gnomad_ac_38)} of {n(m.gnomad_an_38)} alleles with {m.gnomad_hc_38} homozygotes (population-maximum frequency {float(m.gnomad_popmax_af_38) * 100:.2f}%).")
        D.s(f"v38 adds {' and '.join(gained)} beside PVS1, and the combining rule of Section 2.1 of the manuscript (a call on one side stands unless the other side holds evidence of Strong strength or above without reaching a call, and the variant is then a VUS) makes the variant a VUS; {ruling}.")
    keys38 = opposite_of(main.dir38)
    seen = sum(1 for r in keys38.itertuples() if ((rulings.chrom == r.chrom) & (rulings.pos == r.pos) & (rulings.ref == r.ref) & (rulings.alt == r.alt)).any())
    if not len(entered):
        D.s(f"No call entered the set: the opposite-direction calls of v38 are a subset of those of v37, and {seen} of the {len(keys38)} are among the {len(rulings)} rows of `round6_rulings.tsv`, the calls that the reviewer ruled on in round 6.")
    else:
        D.s(f"{len(entered)} call entered the set, and {seen} of the {len(keys38)} opposite-direction calls of v38 are among the {len(rulings)} rows of `round6_rulings.tsv`.")
    oc = opp[opp.pair_v37 == "v37_oct2026"]
    same_oct = len(left) and len(oc) and set(map(tuple, left[KEY].itertuples(index=False))) == set(map(tuple, oc[oc.change == "left_in_v38"][KEY].itertuples(index=False))) and not (oc.change == "entered_in_v38").any()
    D.s(f"On the October release the 2-star+ variants have {october.x37['opposite']} opposite-direction calls under v37 and {october.x38['opposite']} under v38 ({october.x37['false_benign']} and {october.x37['false_pathogenic']}, {october.x38['false_benign']} and {october.x38['false_pathogenic']}); " + ("the same variant leaves it and none enters." if same_oct else f"{len(oc[oc.change == 'left_in_v38'])} leave it and {len(oc[oc.change == 'entered_in_v38'])} enter it (Section 8)."))
    first = f"{left.iloc[0].gene} {left.iloc[0].chrom}:{left.iloc[0].pos}" if len(left) else ""
    gone = opp[(opp.change == "left_in_v38") & (opp.gene == left.iloc[0].gene) & (opp.pos == left.iloc[0].pos)].pair_v37.nunique() if len(left) else 0
    D.s(f"{first} is an opposite-direction call under v37 and not under v38 in {gone} of the {len(counts)} pairs that have a discrepancy file on both sides, so the table below lists, for each other pair, only what differs beyond that one variant.")
    first_keys = {tuple(r) for r in left[KEY].itertuples(index=False)}
    rows, n_same, n_only = [], 0, 0
    considered = [p for p in pairs if p.v37 in counts]
    for p in considered:
        c37, c38, nl, ne = counts[p.v37]
        o = opp[opp.pair_v37 == p.v37]
        lf, en = o[o.change == "left_in_v38"], o[o.change == "entered_in_v38"]
        beyond = lf[[tuple(r) not in first_keys for r in lf[KEY].itertuples(index=False)]]
        if not len(lf) and not len(en):
            n_same += 1
        elif not len(beyond) and not len(en):
            n_only += 1
        else:
            rows.append([p.v37, p.v38, str(c37), str(c38), "; ".join(f"{r.gene} {r.chrom}:{r.pos}" for r in beyond.itertuples()) or "-", "; ".join(f"{r.gene} {r.chrom}:{r.pos}" for r in en.itertuples()) or "-"])
    D.s(f"Of the {len(considered)} other pairs, {n_same} have the same opposite-direction calls under both classifiers, {n_only} differ only in {first}, and {len(rows)} differ in others:")
    D.table(["Run under v37", "Run under v38", "Calls v37", "Calls v38", "Left beyond the one above", "Entered"], rows, align="llrrll")
    n_left = opp[opp.change == "left_in_v38"].drop_duplicates(["gene", "chrom", "pos", "ref", "alt"]).shape[0]
    n_entered = opp[opp.change == "entered_in_v38"].drop_duplicates(["gene", "chrom", "pos", "ref", "alt"]).shape[0]
    D.s(f"Over all the pairs {n_left} distinct variants are opposite-direction calls under v37 and not under v38 in at least one pair, and {n_entered} the other way (`v38_family_vs_v37_opposite.csv` lists each with its criteria).")
    for change, title in (("left_in_v38", "The opposite-direction calls that v38 does not make where v37 makes them (in at least one pair)"),
                          ("entered_in_v38", "The opposite-direction calls that v38 makes where v37 does not (in at least one pair)")):
        d = opp[opp.change == change]
        rows = []
        for key, g in d.groupby(["gene", "chrom", "pos", "ref", "alt"], sort=False):
            first = g.iloc[0]
            who = sorted(set(g.pair_v37), key=lambda t: (len(t), t))
            rows.append([f"{key[0]} {key[1]}:{key[2]} {allele(key[3])}>{allele(key[4])}", f"{first.clinvar.replace('_', ' ')} ({first.stars} stars)",
                         f"{first.call_in_the_run_that_has_it.replace('_', ' ')} on {first.criteria_in_the_run_that_has_it.replace(';', ', ')}", str(len(who)), ", ".join(f"`{w_}`" for w_ in who[:3]) + (f" and {len(who) - 3} more" if len(who) > 3 else "")])
        rows.sort(key=lambda r: (-int(r[3]), r[0]))
        D.s(f"{title}, with the call and the criteria of the run that has it:")
        D.table(["Variant", "ClinVar", "Call", "Pairs", "Pairs in which it differs (run under v37)"], rows, align="lllrl")
    bs2 = opp[opp.pair_v37.isin([p.v37 for p in pairs if p.knob in ("bs2_ad_prevalence_threshold", "bs2_ad_min_ac (bar off)")])]
    en_bs2 = bs2[bs2.change == "entered_in_v38"]
    D.s(f"{len(en_bs2.drop_duplicates(['gene', 'chrom', 'pos', 'ref', 'alt']))} of the {n_entered} that enter are in the pairs of the BS2 carrier threshold and floor, which is where Section 3.8 finds the new reach of the carrier arm: with the threshold off, v38 has {pairs_by_name(pairs, 'v37_bs2off').x38['opposite']} opposite-direction calls where v37 has {pairs_by_name(pairs, 'v37_bs2off').x37['opposite']}.")


RULING_TEXT = {"clinvar": "the reviewer supported ClinVar's call", "fastvep": "the reviewer supported fastVEP's call", "neither": "the reviewer's own call was VUS",
               "needs_evidence": "the reviewer said that more evidence was needed", "disorder_specific": "the reviewer said that the answer depends on the disorder"}
RULE_TEXT = {"class": "by a class that the reviewer stated", "stated": "by the reviewer's statement that ClinVar is correct", "evidence": "by evidence for ClinVar's class without naming one"}


def ruling_text(ruling: str, rule: str) -> str:
    """The reviewer's ruling in the words of Section 3.6 of the manuscript, and the rule by which we coded it."""
    how = f" ({RULE_TEXT[rule]})" if ruling == "clinvar" and rule in RULE_TEXT else f" (coded by us as {rule})"
    return RULING_TEXT.get(ruling, ruling) + how


def allele(a: str) -> str:
    """An allele as it is, or its ends and its length when it is long."""
    return a if len(a) <= 12 else f"{a[:4]}...{a[-3:]} ({len(a)} bp)"


def pairs_by_name(pairs: list[Pair], name37: str) -> Pair:
    return next(p for p in pairs if p.v37 == name37)

# ---------------------------------------------------------------------------------------------- 4. the sweeps


def fmt_setting(v: str) -> str:
    return v


def shape_changes(summary: list[dict]) -> list[tuple[str, str]]:
    """The sweeps whose shape differs between the families: the knee moved, the signs of the steps differ, or settings change order."""
    out = []
    for sw in summary:
        sh = sw["shape"]
        a, b = sh["v37"], sh["v38"]
        reasons = []
        if a["knee"] != b["knee"]:
            reasons.append("the knee moved")
        if a["gain_signs"] != b["gain_signs"] or a["cost_signs"] != b["cost_signs"]:
            reasons.append("the signs of the steps differ")
        if sh["flips"]:
            reasons.append(f"{plural(len(sh['flips']), 'pair of settings changes', 'pairs of settings change')} order, within a range of {n(a['range_same_direction'])} calls under v37 and {n(b['range_same_direction'])} under v38")
        if reasons:
            out.append((sw["label"], "; ".join(reasons)))
    return out


def summary_rows(summary: list[dict], eff: pd.DataFrame) -> list[dict]:
    """The counts that the manuscript quotes from this comparison and that no table of `manuscript/acmg/data/` holds (`verify_numbers.py`
    reads them from here): how many sweeps there are and how many differ, and how many effects there are and how many changed."""
    group_eff = eff[eff.group.isin(["switch", "source", "hook", "gap"])]
    flagged = group_eff[(group_eff.changed_sign == True) | (group_eff.changed_by_more_than_20_percent == True) | (group_eff.appeared_or_vanished == True)]  # noqa: E712
    knee = {s_["label"]: (s_["shape"]["v37"]["knee"], s_["shape"]["v38"]["knee"]) for s_ in summary}
    PM2_2015 = "PM2 bar, not recessive-only, 2015 table both sides"
    fall = {s_["label"]: ((0, s_["shape"]["v37"]["fall_from_rate"], s_["shape"]["v37"]["fall_to_rate"]), (0, s_["shape"]["v38"]["fall_from_rate"], s_["shape"]["v38"]["fall_to_rate"])) for s_ in summary}
    return [
        {"name": "sweeps", "value": len(summary)},
        {"name": "sweeps_that_differ_in_shape", "value": len(shape_changes(summary))},
        {"name": "sweeps_that_differ_in_the_knee", "value": sum(1 for a_, b_ in knee.values() if a_ != b_)},
        {"name": "knee_bs1_v37", "value": knee["BS1 bar"][0]}, {"name": "knee_bs1_v38", "value": knee["BS1 bar"][1]},
        {"name": "knee_bp7_far_boundary_v37", "value": knee["BP7 far boundary"][0]}, {"name": "knee_bp7_far_boundary_v38", "value": knee["BP7 far boundary"][1]},
        {"name": "knee_pm2_2015_v37", "value": knee[PM2_2015][0]}, {"name": "knee_pm2_2015_v38", "value": knee[PM2_2015][1]},
        {"name": "fall_pm2_2015_v37_from", "value": round(fall[PM2_2015][0][1])}, {"name": "fall_pm2_2015_v37_to", "value": round(fall[PM2_2015][0][2])},
        {"name": "fall_pm2_2015_v38_from", "value": round(fall[PM2_2015][1][1])}, {"name": "fall_pm2_2015_v38_to", "value": round(fall[PM2_2015][1][2])},
        {"name": "effects", "value": len(group_eff)},
        {"name": "effect_measures", "value": int(group_eff.metric.nunique())},
        {"name": "effects_that_changed", "value": len(flagged)},
        {"name": "effects_that_changed_and_are_small", "value": int((flagged.small_number == True).sum())},  # noqa: E712
    ]


def knee_sentence(a: dict, b: dict) -> str:
    fa, fb = f"{a['knee_factor']:.1f}", f"{b['knee_factor']:.1f}"
    ka = None if a["knee"] is None else setting(f"{a['knee']:g}")
    kb = None if b["knee"] is None else setting(f"{b['knee']:g}")
    if ka is None and kb is None:
        return f"There is no knee under either classifier (the largest fall in the exchange rate is by a factor of {fa} under v37 and {fb} under v38)."
    if ka == kb:
        return f"The knee is at {ka} under both classifiers (the largest fall in the exchange rate is by a factor of {fa} under v37 and {fb} under v38)."
    return f"The knee is {'absent' if ka is None else 'at ' + ka} under v37 (the largest fall in the exchange rate is by a factor of {fa}) and {'absent' if kb is None else 'at ' + kb} under v38 ({fb}), so it has moved."


def rate_sentence(a: dict, b: dict) -> str:
    def one(x: dict) -> str:
        return "no step adds an error" if x["max_rate"] == 0 and x["min_rate"] == 0 else f"the exchange rate of the steps lies between {n(x['min_rate'])} and {n(x['max_rate'])} calls per added error"
    ra, rb = one(a), one(b)
    if ra == rb:
        return f"Under both classifiers {ra} (steps that add no error are left out of the rates)."
    return f"Under v37 {ra}, and under v38 {rb} (steps that add no error are left out of the rates)."


def section_sweeps(D: Doc, summary: list[dict]) -> None:
    D.h(2, "4. The sweeps")
    D.s("Every sweep is a curve of one threshold, the default among its settings, and each is given under both classifiers.")
    D.s("For each sweep the tables give agreement in direction, the calls that the sweep is meant to gain (correct pathogenic-side calls for the PM2 sweeps and the point cutoffs, correct benign-side calls for the rest) and the two kinds of reference-discordant call.")
    D.s("The settings are in the order of their values, and the row marked default is the point at which the curve passes through the default (see the note under each of the two sweeps that are not made at the default combination of the other settings).")
    for sw in summary:
        sh, pts = sw["shape"], sw["points"]
        D.h(3, f"4.{summary.index(sw) + 1} {sw['label']}")
        gain = sw["gain"]
        gl = "ClinVar P/LP called P/LP" if gain == "same_pathogenic" else "ClinVar B/LB called LB/B"
        rows = []
        for o in pts:
            rows.append([setting(o["value"], sw["knob"]) + (" (default)" if o["shipped"] else ""), n(o["x37"]["same_direction"]), n(o["x38"]["same_direction"]), n(o["x37"][gain]), n(o["x38"][gain]),
                         str(o["x37"]["false_benign"]), str(o["x38"]["false_benign"]), str(o["x37"]["false_pathogenic"]), str(o["x38"]["false_pathogenic"])])
        D.table(["Setting", "Agreement v37", "Agreement v38", f"{gl} v37", f"{gl} v38", "Ref.-discordant benign v37", "v38", "Ref.-discordant pathogenic v37", "v38"], rows)
        if sw["knob"].startswith("pm2_ad_af_threshold (Richards"):
            D.s("The curve is made with the 2015 combining table on both sides, so its point at the default threshold of 4e-5 is the run `sw_table_both_sides` (the table and nothing else changed) and not the main run.")
        if sw["knob"].startswith("bs2_ad_min_ac"):
            D.s("The curve is made with the carrier threshold off, so its point at the default floor of 5 alleles is the run with the threshold off (`v37_bs2off`, `v38_sw_bs2ad_0`) and not the main run.")
        a, b = sh["v37"], sh["v38"]
        D.s(f"Across its settings the sweep moves agreement in direction by {n(a['range_same_direction'])} calls under v37 and {n(b['range_same_direction'])} under v38, the calls it gains by {n(a['range_gain'])} and {n(b['range_gain'])}, reference-discordant benign calls by {a['range_false_benign']} and {b['range_false_benign']}, and reference-discordant pathogenic calls by {a['range_false_pathogenic']} and {b['range_false_pathogenic']}.")
        same_sign = a["gain_signs"] == b["gain_signs"] and a["cost_signs"] == b["cost_signs"]
        D.s(f"The signs of the steps along the settings (the calls gained, then the errors added; + more, - fewer, 0 the same) are {a['gain_signs']} and {a['cost_signs']} under v37 and {b['gain_signs']} and {b['cost_signs']} under v38, so the shape is {'the same' if same_sign else 'not the same'} step for step.")
        D.s(knee_sentence(a, b))
        D.s(rate_sentence(a, b))
        if sh["flips"]:
            fl = ", ".join(f"{setting(u, sw['knob'])} and {setting(v, sw['knob'])}" for u, v in sh["flips"][:6]) + (f" and {len(sh['flips']) - 6} more" if len(sh["flips"]) > 6 else "")
            D.s(f"The order of the settings by agreement in direction differs between the two: {plural(len(sh['flips']), 'pair of settings changes', 'pairs of settings change')} order ({fl}), within a range of {n(a['range_same_direction'])} calls under v37 and {n(b['range_same_direction'])} under v38.")
        else:
            D.s("The order of the settings by agreement in direction is the same under both.")
    D.h(3, f"4.{len(summary) + 1} What changed in shape")
    changed_shape = shape_changes(summary)
    D.s("Three things were compared for every sweep: the sign of each step, the position of the knee, and the order of the settings by agreement in direction.")
    D.s(f"{len(changed_shape)} of the {len(summary)} sweeps differ in at least one of them:")
    D.bullets([f"{label}: {why}." for label, why in changed_shape])
    D.s("A sweep whose whole range is a few dozen calls reorders its settings by a few calls, and that is not a change of shape in the sense that matters; the ranges are given so that it can be told.")
    kn37 = [sw["label"] for sw in summary if sw["shape"]["v37"]["knee"] is not None]
    kn38 = [sw["label"] for sw in summary if sw["shape"]["v38"]["knee"] is not None]
    D.s(f"The sweeps with a knee are the same under both classifiers: {'; '.join(kn38)}." if kn37 == kn38 else f"The sweeps with a knee are {'; '.join(kn37)} under v37 and {'; '.join(kn38)} under v38.")


# ---------------------------------------------------------------------------------------------- 5. the ablations


def short_name(name37: str) -> str:
    return name37.removeprefix("v37_")


def section_effects(D: Doc, pairs: list[Pair], main: Pair, eff: pd.DataFrame) -> None:
    D.h(2, "5. The ablations, the departures from the sources and the inputs that change")
    D.s("The effect of a run is the run less the default under the same classifier, for agreement in direction, for the correct pathogenic-side and benign-side calls and for the two kinds of reference-discordant call.")
    D.s("The tables give the effect under v37 and under v38 for each run of the groups switch, source, hook and gap, and the flags of Section 1 follow them.")
    metrics = [("same_direction", "agreement"), ("same_pathogenic", "P/LP called P/LP"), ("same_benign", "B/LB called LB/B"), ("false_benign", "ref.-discordant benign"), ("false_pathogenic", "ref.-discordant pathogenic")]
    titles = {"switch": "5.1 Single settings switched off (on for the two opt-in settings)", "source": "5.2 One evidence source removed", "hook": "5.3 The departures that the second measurement build exposes",
              "gap": "5.4 The inputs of the PVS1 gap"}
    e = eff[eff.metric.isin([m for m, _ in metrics])]
    for group, title in titles.items():
        D.h(3, title)
        sub = e[e.group == group]
        rows = []
        for run37 in dict.fromkeys(sub.run_v37):
            r = sub[sub.run_v37 == run37]
            row = [short_name(run37)]
            for m, _ in metrics:
                x = r[r.metric == m].iloc[0]
                row += [sg(x.effect_v37), sg(x.effect_v38)]
            rows.append(row)
        D.table(["Run"] + [f"{lbl} {v}" for _, lbl in metrics for v in ("v37", "v38")], rows)
    D.h(3, "5.5 The effects that changed")
    flag = eff[(eff.group.isin(["switch", "source", "hook", "gap"])) & ((eff.changed_sign == True) | (eff.changed_by_more_than_20_percent == True) | (eff.appeared_or_vanished == True))]  # noqa: E712
    allg = eff[eff.group.isin(["switch", "source", "hook", "gap"])]
    D.s(f"Of the {len(allg)} effects of these runs ({allg.run_v37.nunique()} runs and {allg.metric.nunique()} measures), {len(flag)} changed sign, appeared or vanished, or changed by more than 20%, and {int((flag.small_number == True).sum())} of those are small (fewer than 10 calls under both classifiers).")  # noqa: E712
    big = flag[flag.small_number == False]  # noqa: E712
    lines = []
    for r in big.itertuples():
        rel = "" if r.relative_change == "" else f", {float(r.relative_change) * 100:.0f}% of the effect under v37"
        extra = " and changed sign" if r.changed_sign else ""
        lines.append(f"{short_name(r.run_v37)}, {r.metric.replace('_', ' ')}: {sg(r.effect_v37)} under v37 and {sg(r.effect_v38)} under v38 (a difference of {sg(r.effect_delta)}{rel}){extra}.")
    D.bullets(lines)
    sm = flag[flag.small_number == True]  # noqa: E712
    D.s("The small ones, which a single variant makes:")
    D.bullets([f"{short_name(r.run_v37)}, {r.metric.replace('_', ' ')}: {sg(r.effect_v37)} and {sg(r.effect_v38)}." for r in sm.itertuples()])
    P = {p_.v37: p_ for p_ in pairs}

    def effect(name37: str, metric: str) -> tuple[int, int]:
        r_ = eff[(eff.run_v37 == name37) & (eff.metric == metric)].iloc[0]
        return int(r_.effect_v37), int(r_.effect_v38)

    def firing(name37: str, code: str) -> tuple[int, int]:
        q = P[name37]
        return q.f37.get(code, 0) - main.f37.get(code, 0), q.f38.get(code, 0) - main.f38.get(code, 0)

    def benign_to_vus(name37: str) -> tuple[int, int]:
        q = P[name37]
        return (sum(q.m37[t]["VUS"] - main.m37[t]["VUS"] for t in ("Likely_benign", "Benign")), sum(q.m38[t]["VUS"] - main.m38[t]["VUS"] for t in ("Likely_benign", "Benign")))

    nc37, nc38 = effect("v37_src_noclingen", "same_benign")
    D.s(f"The one effect of size on the benign side is the removal of the ClinGen file (`src_noclingen`): it costs {n(abs(nc38))} correct benign-side calls under v38 where it costs {n(abs(nc37))} under v37.")
    pm37, pm38 = firing("v37_src_noclingen", "PM2_Supporting")
    bs37, bs38 = firing("v37_src_noclingen", "BS2")
    others = sorted((abs((q38 - q37)) for q37, q38 in (firing("v37_src_noclingen", c_) for c_ in set(main.f37) | set(main.f38) if c_ not in ("PM2_Supporting", "BS2"))), reverse=True)
    v37_up, v38_up = benign_to_vus("v37_src_noclingen")
    D.s(f"The ClinVar benign-side variants that the removal turns into VUS are {n(v38_up)} under v38 against {n(v37_up)} under v37, PM2 at Supporting strength fires on {n(abs(pm38))} {'fewer' if pm38 < 0 else 'more'} variants under v38 and on {n(abs(pm37))} {'fewer' if pm37 < 0 else 'more'} under v37, and BS2 fires on {n(abs(bs38))} {'more' if bs38 > 0 else 'fewer'} variants under v38 and on {n(abs(bs37))} {'more' if bs37 > 0 else 'fewer'} under v37.")
    D.s(f"No other criterion differs between the two families by more than {others[0]} variants in this pair (`criterion_firing_rates.csv`).")
    D.s("PM2 and BS2 are the two criteria that branch on the mode of inheritance and can be met in the benchmark (Section 3.8 of the manuscript), and the removal of the file takes away the field that v38 reads for them and the names that v37 read; no per-variant table was made for this run, so which of the two accounts for the loss of benign-side calls is not traced here.")
    cc37, cc38 = effect("v37_gap_clingen_complete", "same_benign")
    gb37, gb38 = effect("v37_gap_both", "same_benign")
    cpm37, cpm38 = firing("v37_gap_clingen_complete", "PM2_Supporting")
    cbs37, cbs38 = firing("v37_gap_clingen_complete", "BS2")
    orph_same = all(effect("v37_gap_orphanet_lof", m_)[0] == effect("v37_gap_orphanet_lof", m_)[1] for m_ in ("same_direction", "same_pathogenic", "same_benign", "false_benign", "false_pathogenic"))
    D.s(f"The completed ClinGen file (`gap_clingen_complete`) and the file with Orphanet's label (`gap_both`) change the correct benign-side calls by {sg(cc38)} and {sg(gb38)} under v38 where they changed them by {sg(cc37)} and {sg(gb37)} under v37.")
    D.s(f"Under v38 the completed file changes the firing of PM2 at Supporting strength by {sg(cpm38)} and of BS2 by {sg(cbs38)}, where under v37 it changed them by {sg(cpm37)} and {sg(cbs37)} (PM2 and BS2 are the two criteria that branch on the mode of inheritance).")
    D.s(f"Orphanet's label alone (`gap_orphanet_lof`) has {'the same effect under both classifiers on agreement in direction, on the two recalls and on both kinds of reference-discordant call' if orph_same else 'a different effect under the two classifiers'}.")
    D.s("Every other effect of these runs has the same sign, is within 20% of its value under v37 and neither appeared nor vanished.")
    for group, label in (("cutoff", "the runs of the point cutoffs"), ("interaction", "the runs of the PM2 threshold under the 2015 table"), ("sweep", "the runs of the sweeps")):
        g = eff[(eff.group == group) & eff.metric.isin(["same_direction", "false_benign", "false_pathogenic"])]
        fl = g[(g.changed_sign == True) | (g.changed_by_more_than_20_percent == True)]  # noqa: E712
        D.s(f"For {label} ({g.run_v37.nunique()} runs), {len(fl)} of the {len(g)} effects on agreement and on the two kinds of reference-discordant call changed sign or by more than 20%; the sweeps are examined as curves in Section 4.")

# ---------------------------------------------------------------------------------------------- 6 and 7. controls and criteria


def section_controls(D: Doc, pairs: list[Pair], main: Pair) -> None:
    D.h(2, "6. The controls, and what the second family reproduces")
    c = pd.read_csv(DATA / "controls_v38.csv")
    what = {"v38_cut_lp6": "the control of the cutoff build, with the Likely pathogenic cutoff at its default of 6 points", "v38_hk_control": "the control of the hook build, with nothing set"}
    for r in c.itertuples():
        D.s(f"`{r.run}`, {what[r.run]}, equals the main run `v38` in {r.cells_equal_to_v38} of {r.cells} cells of the concordance matrix and in the firing count of {r.criteria_firing_equal_to_v38} of {r.criteria} criterion codes (`controls_v38.csv`).")
    c37 = pd.read_csv(DATA / "controls_v37.csv")
    D.s(f"Both builds therefore behave as the default classifier apart from the settings they expose, as the controls of the family v37 do (`controls_v37.csv`: {'; '.join(f'`{r.run}` {r.cells_equal_to_v37} of {r.cells} cells and {r.criteria_firing_equal_to_v37} of {r.criteria} codes' for r in c37.itertuples())}).")
    same_files = []
    for ctl in ("v38_cut_lp6", "v38_hk_control"):
        ok = all((BENCH / "manuscript_runs" / ctl / f).read_bytes().replace(b"\r\n", b"\n") == (BENCH / "manuscript_runs/v38" / f).read_bytes().replace(b"\r\n", b"\n")
                 for f in ("concordance_matrix.csv", "criterion_firing_rates.csv", "discrepancies.tsv", "rule_distribution.csv", "concordance_by_chrom.csv", "concordance_by_consequence.csv"))
        same_files.append((ctl, ok))
    D.s("Their scored files are the main run's byte for byte: " + "; ".join(f"`{ctl}` {'identical' if ok else 'DIFFERENT'} in the matrix, the firing rates, the discrepancies, the rule distribution and the two per-group tables" for ctl, ok in same_files) + ".")
    family = [p_ for p_ in pairs if p_.dir37.parent == BENCH / "manuscript_runs"]
    new_runs = [p_.v38 for p_ in family if not is_shared(p_.v38)]
    keys = {k for p_ in family for k in prov(p_.v37)}
    complete = [x for x in new_runs if set(prov(x)) == keys and all(v.strip() for v in prov(x).values())]
    D.s(f"The provenance file of {len(complete)} of the {len(new_runs)} new runs of the family has the {len(keys)} fields that every file of the family v37 has (flags, SHA-256 of the binary and of the configuration, the input, GFF3 and FASTA by name and size, the signature of the database directory, the version of the binary, the commit and the number of uncommitted files of the checkout that ran the harness, the annotation time and the time of the record), and none of them is empty; the other files of the directory are the {len(prov_files()) - len(new_runs)} that existed before.")
    pv = {x: prov(x) for x in new_runs if not x.startswith(("v38_gap_", "v38_src_", "v38_oct_"))}
    same_inputs = {(v["sa_dir_signature"], v["input_bytes"], v["gff3_bytes"], v["fasta_bytes"]) for v in pv.values()}
    D.s(f"The {len(pv)} new April runs that read the unchanged databases agree in the signature of the database directory, the sizes of the input, the GFF3 and the FASTA ({len(same_inputs)} distinct combination{'s' if len(same_inputs) != 1 else ''}), and {'the signature is the one that the main run of Section 3.8 records' if same_inputs and {v[0] for v in same_inputs} == {prov('v38')['sa_dir_signature']} else 'the signature is not the one that the main run of Section 3.8 records'} (`{prov('v38')['sa_dir_signature'][:8]}...`).")
    kept37 = [p_ for p_ in family if list(p_.dir37.glob("*.vcf.gz"))]
    kept38 = [p_ for p_ in family if list(p_.dir38.glob("*.vcf.gz"))]
    same_kept = {p_.v37.removeprefix("v37_") for p_ in kept37} == {p_.v38.removeprefix("v38_") for p_ in kept38}
    D.s(f"The annotated VCF of a run is deleted after scoring unless the run is one whose per-variant table the statistics read: the family v38 keeps it for {sum(1 for p_ in kept38 if p_.kind == 'april')} April runs and {sum(1 for p_ in kept38 if p_.kind == 'october')} October runs, {'the same runs as' if same_kept else 'not the same runs as'} the family v37 ({sum(1 for p_ in kept37 if p_.kind == 'april')} and {sum(1 for p_ in kept37 if p_.kind == 'october')}; `KEEP_VCF` of `11_manuscript_runs.py`).")


def prov_files() -> list[Path]:
    return sorted((RUNS_DIR / "provenance").glob("*.txt"))


def section_criteria(D: Doc, pairs: list[Pair], main: Pair, october: Pair) -> None:
    D.h(2, "7. Criteria whose firing counts moved by more than 5%")
    D.s("A criterion's firing count is the number of variants that carry it at any strength, counted for each code (a strength is a code of its own), as `criterion_firing_rates.csv` of each run holds it.")
    rows = []
    for c in sorted(set(main.f37) | set(main.f38)):
        a, b = main.f37.get(c, 0), main.f38.get(c, 0)
        rel_ = (b - a) / a if a else float("nan")
        rows.append([c, n(a), n(b), sg(b - a), "-" if a == 0 else f"{rel_ * 100:+.2f}%"])
    D.s("The main run, every code:")
    D.table(["Code", "v37", "v38", "Change", "Relative"], rows)
    moved = [c for c in set(main.f37) | set(main.f38) if main.f37.get(c, 0) and abs(main.f38.get(c, 0) - main.f37[c]) / main.f37[c] > 0.05]
    D.s(f"In the main run {len(moved)} code moved by more than 5%: {', '.join(sorted(moved))}" + (f" ({main.f37[moved[0]]:,} to {main.f38[moved[0]]:,}, {(main.f38[moved[0]] - main.f37[moved[0]]) / main.f37[moved[0]] * 100:+.1f}%)" if len(moved) == 1 else "") + ".")
    D.s(f"PM2 at Supporting strength moves by {(main.f38['PM2_Supporting'] - main.f37['PM2_Supporting']) / main.f37['PM2_Supporting'] * 100:+.1f}% ({n(main.f37['PM2_Supporting'])} to {n(main.f38['PM2_Supporting'])}), below the line; Section 3.8 attributes it to the {n(st_recessive())} genes that v38 reads as recessive-only, which use PM2's higher threshold.")
    flagged = []
    for p in [october] + pairs:
        for c in sorted(set(p.f37) | set(p.f38)):
            a, b = p.f37.get(c, 0), p.f38.get(c, 0)
            if a and abs(b - a) / a > 0.05 and c not in moved:
                flagged.append((p.v37, p.v38, c, a, b))
    D.s(f"Across the other {len(pairs) + 1} pairs, {len(flagged)} firing counts of a code that does not move in the main run moved by more than 5%:")
    D.bullets([f"{x37} against {x38}, {c}: {n(a)} to {n(b)} ({(b - a) / a * 100:+.1f}%)." for x37, x38, c, a, b in flagged])
    if flagged:
        pm_runs = all(c_ == "PM2_Supporting" and x37.startswith(("v37_sw_pm2ad_", "v37_sw_pm2ar_", "v37_ix_table_pm2ad_")) for x37, x38, c_, a, b in flagged)
        D.s(f"{'All of them are' if pm_runs else 'They include'} the firing of PM2 at Supporting strength in runs that change PM2's threshold.")
        if pm_runs:
            by_knob = {}
            for p_ in pairs:
                if p_.group == "sweep" and p_.knob in ("pm2_ad_af_threshold", "pm2_ar_af_threshold"):
                    by_knob.setdefault(p_.knob, []).append(p_)
            for knob, label in (("pm2_ad_af_threshold", "genes not read as recessive-only"), ("pm2_ar_af_threshold", "recessive-only genes")):
                run_list = sorted(by_knob[knob], key=lambda q: float(q.value))
                lo, hi = run_list[0], run_list[-1]
                D.s(f"Across the sweep of PM2's threshold for {label} ({lo.value} to {hi.value}) its firing goes from {n(lo.f37['PM2_Supporting'])} to {n(hi.f37['PM2_Supporting'])} under v37 and from {n(lo.f38['PM2_Supporting'])} to {n(hi.f38['PM2_Supporting'])} under v38.")
            D.s(f"Each threshold therefore reaches a different share of the variants under the two classifiers, which is the direction that the move of the genes read as recessive-only from {n(st_recessive_v37())} to {n(st_recessive())} implies, because the classifier uses one PM2 threshold for them and another for the rest (Section 3.3 of the manuscript); the sweeps as curves are in Section 4.")


def st_recessive() -> int:
    f = pd.read_csv(DATA / "facts_v38.csv").set_index("name")["value"]
    return int(float(f["inheritance_read_recessive_only"]))


def st_recessive_v37() -> int:
    f = pd.read_csv(DATA / "facts_v37.csv").set_index("name")["value"]
    return int(float(f["inheritance_read_recessive_only"]))

# ---------------------------------------------------------------------------------------------- 8. the October release


def section_october(D: Doc, pairs: list[Pair], main: Pair, october: Pair, opp: pd.DataFrame) -> None:
    D.h(2, "8. The October release")
    D.s("The classifier of each family was applied unchanged to the 781,135 variants of ClinVar's release of 4 October 2026 at two stars or better (`temporal_v37.csv` and `temporal_v38.csv`).")
    ta, tb = pd.read_csv(DATA / "temporal_v37.csv"), pd.read_csv(DATA / "temporal_v38.csv")
    rows = []
    for stratum in ("all", "unchanged", "tier_change", "direction_change", "new", "new, 2-star", "new, 3-star"):
        x, y = ta[ta.stratum == stratum].iloc[0], tb[tb.stratum == stratum].iloc[0]
        rows.append([stratum, n(x.stratum_n), n(x.exact), n(y.exact), n(x.same_direction), n(y.same_direction), pc(x.same_direction / x.stratum_n, 1), pc(y.same_direction / y.stratum_n, 1),
                     f"{x.opposite} ({x.false_benign} and {x.false_pathogenic})", f"{y.opposite} ({y.false_benign} and {y.false_pathogenic})"])
    D.table(["Stratum", "Variants", "Exact v37", "Exact v38", "In direction v37", "In direction v38", "Rate v37", "Rate v38", "Opposite v37", "Opposite v38"], rows)
    al, nw = ta[ta.stratum == "all"].iloc[0], tb[tb.stratum == "all"].iloc[0]
    D.s(f"On the October release v38 agrees in direction for {n(nw.same_direction)} variants where v37 agrees for {n(al.same_direction)} ({sg(nw.same_direction - al.same_direction)}), with {nw.opposite} opposite-direction calls where v37 has {al.opposite} ({nw.false_benign} reference-discordant benign calls and {nw.false_pathogenic} reference-discordant pathogenic calls, against {al.false_benign} and {al.false_pathogenic}).")
    a2, b2 = ta[ta.stratum == "new"].iloc[0], tb[tb.stratum == "new"].iloc[0]
    same_new_opp = int(a2.opposite) == int(b2.opposite)
    D.s(f"On the {n(a2.stratum_n)} new variants, the ones that are not in the April set, agreement in direction is {pc(a2.same_direction / a2.stratum_n, 1)} under v37 and {pc(b2.same_direction / b2.stratum_n, 1)} under v38, and the opposite-direction calls are {a2.opposite} under v37 ({a2.false_benign} reference-discordant benign calls and {a2.false_pathogenic} reference-discordant pathogenic calls) and {b2.opposite} under v38 ({b2.false_benign} and {b2.false_pathogenic}){', the same number' if same_new_opp else ''}.")
    ra, rb = pd.read_csv(DATA / "temporal_ratio_v37.csv"), pd.read_csv(DATA / "temporal_ratio_v38.csv")
    xa, xb = ra[ra.measure.str.startswith("ratio")].iloc[0], rb[rb.measure.str.startswith("ratio")].iloc[0]
    D.s(f"The rate of opposite-direction calls on the new variants is {xa.value:.2f} times the April rate under v37 (95% interval {xa.ci_lo:.2f} to {xa.ci_hi:.2f}) and {xb.value:.2f} times under v38 ({xb.ci_lo:.2f} to {xb.ci_hi:.2f}), where the April rate is {int(ra.iloc[0].k)} in {n(ra.iloc[0].n)} and {int(rb.iloc[0].k)} in {n(rb.iloc[0].n)}.")
    strata = pd.read_csv(BENCH / "clinvar_20261004/strata.tsv", sep="\t", dtype=str)
    o = opp[opp.pair_v37 == "v37_oct2026"]
    where = []
    for r in o.itertuples():
        st_ = strata[(strata.chrom == r.chrom) & (strata.pos == r.pos) & (strata.ref == r.ref) & (strata.alt == r.alt)].stratum.iloc[0]
        where.append(f"{r.gene} {r.chrom}:{r.pos} ({r.change.replace('_', ' ')}, stratum {st_})")
    D.s(f"The pair `v37_oct2026` and `v38_oct_shipped` differs in {len(o)} opposite-direction call: {'; '.join(where)}, so the {b2.opposite} on the new variants are the same variants under both." if (len(o) == 1 and 'new' not in where[0]) else f"The pair `v37_oct2026` and `v38_oct_shipped` differs in these opposite-direction calls: {'; '.join(where)}.")
    D.h(3, "8.1 The seven changes of Section 4.5 on the new variants")
    pa, pb = pd.read_csv(DATA / "temporal_proposals_v37.csv").set_index("change"), pd.read_csv(DATA / "temporal_proposals_v38.csv").set_index("change")
    rows = []
    for ch in pa.index[1:]:
        x, y = pa.loc[ch], pb.loc[ch]
        rows.append([ch, sg(x.new_delta_same_direction), sg(y.new_delta_same_direction), f"{sg(x.new_delta_false_benign)} / {sg(x.new_delta_false_pathogenic)}", f"{sg(y.new_delta_false_benign)} / {sg(y.new_delta_false_pathogenic)}",
                     sg(x.april_delta_same_direction), sg(y.april_delta_same_direction), f"{sg(x.april_delta_false_benign)} / {sg(x.april_delta_false_pathogenic)}", f"{sg(y.april_delta_false_benign)} / {sg(y.april_delta_false_pathogenic)}"])
    D.table(["Change", "New: agreement v37", "v38", "New: ref.-disc. benign / pathogenic v37", "v38", "April: agreement v37", "v38", "April: benign / pathogenic v37", "v38"], rows)
    same_new = sum(1 for ch in pa.index[1:] if pa.loc[ch].new_delta_false_benign == pb.loc[ch].new_delta_false_benign and pa.loc[ch].new_delta_false_pathogenic == pb.loc[ch].new_delta_false_pathogenic)
    D.s(f"On the new variants the change in the number of reference-discordant calls is the same under both classifiers for {same_new} of the 7 changes, and the change in calls in direction differs by at most {max(abs(pb.loc[ch].new_delta_same_direction - pa.loc[ch].new_delta_same_direction) for ch in pa.index[1:]):.0f} calls.")
    cc = "the completed ClinGen file"
    D.s(f"On the April set the completed ClinGen file adds {n(pa.loc[cc].april_delta_same_direction)} calls in direction and {n(pa.loc[cc].april_delta_false_pathogenic)} reference-discordant pathogenic calls under v37 and {n(pb.loc[cc].april_delta_same_direction)} and {n(pb.loc[cc].april_delta_false_pathogenic)} under v38.")
    ap = opp[opp.pair_v37 == "v37_gap_clingen_complete"]
    D.s("The opposite-direction calls that differ in the pair `v37_gap_clingen_complete` and `v38_gap_clingen_complete` are " + "; ".join(f"{r.gene} {r.chrom}:{r.pos} ({r.change.replace('_', ' ')})" for r in ap.itertuples()) + ".")
    for nm in ("PVS1 licensed by Orphanet's loss-of-function label", "both inputs of the PVS1 gap"):
        x, y = pa.loc[nm], pb.loc[nm]
        D.s(f"{nm[0].upper() + nm[1:]}: {n(x.april_delta_same_direction)} calls in direction for {int(x.april_delta_false_pathogenic)} reference-discordant pathogenic calls under v37 ({x.april_delta_same_direction / x.april_delta_false_pathogenic:.0f} per added error) and {n(y.april_delta_same_direction)} for {int(y.april_delta_false_pathogenic)} under v38 ({y.april_delta_same_direction / y.april_delta_false_pathogenic:.0f} per added error).")

# ---------------------------------------------------------------------------------------------- 9. statements that need rewording


def headings(lines: list[str]) -> list[tuple[int, str]]:
    """The headings of the manuscript with their line numbers, as labels: `Section 3.8`, `Abstract`, `Data and code availability`."""
    out = []
    for i, line in enumerate(lines, 1):
        m = re.match(r"(#{1,3}) (.+)", line)
        if not m:
            continue
        text = m.group(2)
        num = re.match(r"(\d+(?:\.\d+)?)\.? ", text)
        out.append((i, "the title" if len(m.group(1)) == 1 else f"Section {num.group(1)}" if num else text))
    return out


def section_of(line: int, heads: list[tuple[int, str]]) -> str:
    label = heads[0][1]
    for start, name in heads:
        if start <= line:
            label = name
    return label


def run_verify(base: str) -> list[dict]:
    """Every check of `verify_numbers.py` under a base: the line, the anchor and the values the data give."""
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "checks.json"
        r = subprocess.run([sys.executable, str(ROOT / "manuscript/acmg/verify_numbers.py"), "--base", base, "--dump", str(out)], capture_output=True, text=True)
        if not out.exists():
            sys.exit(f"verify_numbers.py --base {base} wrote no dump:\n{r.stdout[-500:]}{r.stderr[-500:]}")
        return json.loads(out.read_text())


def numbers_that_change() -> tuple[list[tuple[int, str, list[tuple[str, str]]]], int, int]:
    """The checks of `verify_numbers.py` whose value differs between the bases: (line, anchor, [(under v37, under v38)])."""
    a, b = run_verify("v37"), run_verify("v38")
    if len(a) != len(b) or any(x["anchor"] != y["anchor"] for x, y in zip(a, b)):
        sys.exit("the checks of the two bases are not the same list")
    out = []
    for x, y in zip(a, b):
        if len(x["values"]) != len(y["values"]):
            # a check that lists a different number of values under the two bases (the genes that a sentence names: six under v37 and five under v38)
            # has values that do not pair by position, so each list is given whole
            d = [(", ".join(x["values"]), ", ".join(y["values"]))]
        else:
            d = [(u, v) for u, v in zip(x["values"], y["values"]) if u != v]
        if d:
            out.append((x["line"], x["anchor"], d))
    return out, len(a), sum(len(x["values"]) for x in a)


def derived_outputs_that_change() -> dict:
    """What the scripts that build the manuscript's tables, workbook and figures make of the files of v37 against what the repository holds (the
    text says v38): the marked tables that differ, the sheets of the workbook and the figures. The workbook and the figures are built into a
    temporary directory, since the ones in the repository are the manuscript's."""
    out: dict = {}
    r = subprocess.run([sys.executable, str(ROOT / "manuscript/acmg/build_tables.py"), "--base", "v37", "--check"], capture_output=True, text=True)
    out["tables"] = [line.split(": ", 1)[1] for line in r.stdout.splitlines() if line.startswith("stale: ")]
    out["n_tables"] = int(re.search(r"(\d+) marked tables", r.stdout).group(1))
    # the cells of Table 3 whose words differ and not only their numbers
    t38 = subprocess.run([sys.executable, str(ROOT / "manuscript/acmg/build_tables.py"), "--print", "thresholds"], capture_output=True, text=True, check=True).stdout.splitlines()
    t37 = subprocess.run([sys.executable, str(ROOT / "manuscript/acmg/build_tables.py"), "--base", "v37", "--print", "thresholds"], capture_output=True, text=True, check=True).stdout.splitlines()
    bare = lambda text: re.sub(r"[0-9][0-9,.]*(e-?\d+)?%?", "#", text)  # noqa: E731
    out["words"] = []
    for x, y in zip(t37, t38):
        cx, cy = [c.strip() for c in x.strip().strip("|").split(" | ")], [c.strip() for c in y.strip().strip("|").split(" | ")]
        for u, v in zip(cx, cy):
            if bare(u) != bare(v):
                out["words"].append((cx[0], u, v))
    with tempfile.TemporaryDirectory() as tmp:
        wb, figs = Path(tmp) / "supplementary_tables_v37.xlsx", Path(tmp) / "figures"
        subprocess.run([sys.executable, str(ROOT / "manuscript/acmg/build_supplement.py"), "--base", "v37", "--write", "--out", str(wb)], check=True, capture_output=True)
        a = pd.read_excel(ROOT / "manuscript/acmg/supplementary_tables.xlsx", sheet_name=None, header=None)
        b = pd.read_excel(wb, sheet_name=None, header=None)
        out["sheets"] = (len(a), [k for k in a if not a[k].equals(b[k])], [k for k in a if a[k].equals(b[k])])
        subprocess.run([sys.executable, str(ROOT / "manuscript/acmg/generate_figures.py"), "--base", "v37", "--out", str(figs)], check=True, capture_output=True)
        pngs = sorted(f.name for f in figs.glob("*.png"))
        same = [f for f in pngs if (ROOT / "manuscript/acmg/figures" / f).read_bytes() == (figs / f).read_bytes()]
        out["figures"] = (len(pngs), [f for f in pngs if f not in same], same)
    return out


def git_out(*args: str) -> str:
    return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True, check=True).stdout


def crate_size(ref: str) -> tuple[int, int, str]:
    """Lines of Rust and unit tests (`#[test]`) of the `src` files of the classification crate at a commit, counted as the manuscript counts them, and its version."""
    lines = tests = 0
    for f in git_out("ls-tree", "-r", "--name-only", ref, "--", "crates/fastvep-classification/src").split():
        if f.endswith(".rs"):
            text = git_out("show", f"{ref}:{f}")
            lines += text.count("\n")
            tests += text.count("#[test]")
    version = re.search(r'^version = "([^"]+)"', git_out("show", f"{ref}:Cargo.toml"), re.M).group(1)
    return lines, tests, version


def regrade_changes(df: pd.DataFrame) -> dict:
    """The variants whose call changes in the replay of Section 3.2 that grades PVS1 one step up (`13_manuscript_stats.py`, `facts`), as {key: (call, call after)}."""
    out = {}
    for key, crit, cons in zip(zip(*[df[k] for k in KEY]), df["criteria"], df["consequence"]):
        cs = [x for x in crit.split("&") if x]
        if not any(x in ("PVS1_Moderate", "PVS1_Supporting") for x in cs):
            continue
        new = []
        for x in cs:
            if x == "PVS1_Moderate" and cons in ("stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant"):
                new.append("PVS1_Strong")
            elif x == "PVS1_Supporting" and cons in ("stop_gained", "frameshift_variant"):
                new.append("PVS1_Moderate")
            else:
                new.append(x)
        a, b = st.replay(cs), st.replay(new)
        if a != b:
            out[key] = (a, b)
    return out


def span(sw: dict, which: str, metric: str, total: int) -> str:
    v = [o[which][metric] for o in sw["points"]]
    return f"{min(v) / total * 100:.2f}% to {max(v) / total * 100:.2f}%"


class Claim:
    """A statement of the manuscript whose truth or wording changes under v38: where it is, what it says and what the data say."""

    def __init__(self, anchors: list[str], says: str, evidence: list[str]) -> None:
        self.anchors, self.says, self.evidence = anchors, says, evidence


def claims_of(pairs: list[Pair], main: Pair, october: Pair, summary: list[dict], cmp_: dict) -> list[Claim]:
    N = APRIL_N
    h37, h38 = main.x37, main.x38
    P = {p.v37: p for p in pairs}
    f37 = pd.read_csv(DATA / "facts_v37.csv").set_index("name")["value"]
    f38 = pd.read_csv(DATA / "facts_v38.csv").set_index("name")["value"]
    sweep = {s["label"]: s for s in summary}
    claims: list[Claim] = []

    # ---- the trajectory
    # v14 to v37 are runs on the first databases (the rows of trajectory.csv), v37 and v38 were run again on the complete ones and their rows are those
    # runs, so the first-database values of the two come from the comparison of the stacks (`21_compare_databases.py`); a band that took v37's row of
    # trajectory.csv would end at 79.07% and not at 77.62%
    traj = pd.read_csv(DATA / "trajectory.csv")
    ovn = pd.read_csv(DATA / "databases_old_vs_new.csv")

    def first_of(version: str, measure: str) -> float:
        return float(ovn[(ovn.version == version) & (ovn.set == "April set") & (ovn.measure == measure)].iloc[0]["first"])

    mainline = traj[traj.mainline.astype(str) == "True"].sort_values("order").set_index("order", drop=False)
    first_rate = {k: float(mainline.loc[k].same_direction_rate) for k in mainline.index if 13 <= k <= 36}
    first_rate[37], first_rate[38] = first_of("v37", "same_direction") / N, first_of("v38", "same_direction") / N
    band_rates = pd.Series([first_rate[k] for k in sorted(first_rate) if 14 <= k <= 37])
    steps = pd.Series([abs(first_rate[k] - first_rate[j]) for j, k in zip(sorted(first_rate), sorted(first_rate)[1:]) if k <= 37])  # from v13's level into the band, as the text counts
    early = (mainline[mainline.order <= 13].same_direction_rate.diff().abs() * 100).dropna().sort_values(ascending=False)
    plp_band = pd.Series([float(mainline.loc[k].plp_recall) for k in mainline.index if 11 <= k <= 36] + [first_of("v37", "plp_recall")])
    blb_band = pd.Series([float(mainline.loc[k].blb_recall) for k in mainline.index if 13 <= k <= 36] + [first_of("v37", "blb_recall")])
    r37, r38 = h37["same_direction"] / N, h38["same_direction"] / N  # the complete databases
    r37f, r38f = first_rate[37], first_rate[38]
    claims.append(Claim(
        ["same-direction agreement on the first databases stayed between", "After v13 the headline stopped moving.", "Agreement in direction with ClinVar reached 77.6% at v13 and stayed between"],
        f"Agreement in direction was flat after the thirteenth version: it stayed between {pc(band_rates.min())} and {pc(band_rates.max())} for the next {len(band_rates)} runs (v14 to v37, without v17), no run moving it by more than {steps.max() * 100:.2f} points; recall of ClinVar's pathogenic-side variants has stayed between {pc(plp_band.min(), 1)} and {pc(plp_band.max(), 1)} since v11, and recall of its benign-side variants between {pc(blb_band.min(), 1)} and {pc(blb_band.max(), 1)} since v13.",
        [f"These are statements about the runs v14 to v37 on the first databases ({len(band_rates)} runs, and the step from v13 into the band: the rows of `trajectory.csv`, and for v37 the first-database run of `databases_old_vs_new.csv`), and they stay true as history.",
         f"v38 on the first databases has {pc(r38f)} agreement in direction, {pa(r38f - r37f)} points above v37 and {pa(r38f - band_rates.max())} points above the top of that band, a step that is {(r38f - r37f) / steps.max():.1f} times the largest step between neighbouring runs of v13 to v37.",
         f"On the complete databases v37 has {pc(r37)} and v38 {pc(r38)} ({n(h38['same_direction'])} calls), {pa(r37 - r37f)} and {pa(r38 - r38f)} points above their runs on the first databases, {(r37 - r37f) / steps.max():.1f} and {(r38 - r38f) / steps.max():.1f} times the largest step of the band.",
         f"Its recall of pathogenic-side variants is {pc(h38['plp_recall'], 1)}, {'inside' if plp_band.min() <= h38['plp_recall'] <= plp_band.max() else 'outside'} the band, and its recall of benign-side variants is {pc(h38['blb_recall'], 1)}, {pa(h38['blb_recall'] - blb_band.max(), 1)} points {'above' if h38['blb_recall'] > blb_band.max() else 'below'} the top of the band, which is of the first databases (the benign-side recall is over the same {n(sum(sum(main.m38[t].values()) for t in ('Likely_benign', 'Benign')))} variants as the trajectory's).",
         f"The manuscript therefore says that the agreement moved once more at v38, by the change of Section 3.8 ({pa(r38f - r37f)} points on the first databases, against {early.iloc[0]:.1f} and {early.iloc[1]:.1f} points for the two largest steps before v13, which line {find_line(manuscript_lines(), 'The two largest steps before v13')} quotes), and that the complete databases moved it by more (Section 3.9).",
         ]))

    # ---- the PM2 sweeps
    pm_ad, pm_ar = sweep["PM2 bar, not recessive-only"], sweep["PM2 bar, recessive-only"]
    rec37, rec38 = int(float(f37["inheritance_read_recessive_only"])), int(float(f38["inheritance_read_recessive_only"]))
    ar_zero = all(pm_ar["shape"][t][k] == 0 for t in ("v37", "v38") for k in ("range_false_benign", "range_false_pathogenic"))
    claims.append(Claim(
        ["On v38 the same sweep is almost flat, and on v37 it was almost flat as well.", "is flat as well: from strict absence to 5e-4", "On v37 the second threshold reached 67 genes, and a flat sweep there said"],
        f"On v37 the sweep of the PM2 threshold for genes not read as recessive-only is almost flat (recall {span(pm_ad, 'x37', 'same_pathogenic', 93_812)}, {n(pm_ad['shape']['v37']['range_gain'])} variants of 93,812), the sweep for recessive-only genes is flat as well (recall {span(pm_ar, 'x37', 'same_pathogenic', 93_812)}, agreement varying by {pm_ar['shape']['v37']['range_same_direction']} calls, no change in either kind of reference-discordant call), and 'the two thresholds therefore acted on nearly every gene and on {rec37} genes, and a flat sweep says what a threshold does to the calls it reaches, not what a threshold would do if inheritance were read from a correct field'.",
        [f"Under v38 the first sweep moves recall across {span(pm_ad, 'x38', 'same_pathogenic', 93_812)} ({n(pm_ad['shape']['v38']['range_gain'])} variants) and reference-discordant pathogenic calls from {min(o['x38']['false_pathogenic'] for o in pm_ad['points'])} to {max(o['x38']['false_pathogenic'] for o in pm_ad['points'])}: it is still almost flat.",
         f"The second moves recall across {span(pm_ar, 'x38', 'same_pathogenic', 93_812)} and agreement by {pm_ar['shape']['v38']['range_same_direction']} calls, and {'changes neither kind of reference-discordant call under either classifier' if ar_zero else 'changes a reference-discordant call under one of the classifiers'}: it is still flat.",
         f"v38 reads {n(rec38)} genes as recessive-only where v37 read {n(rec37)}, so the manuscript gives the count of genes of both versions, and says that the sweep of the recessive-only threshold, which reached {n(rec37)} genes and now reaches {n(rec38)}, stays flat (line {find_line(manuscript_lines(), 'On v37 the second threshold reached 67 genes, and a flat sweep there said')})."]))

    # ---- the PM2 sweep under the 2015 table
    ix = sweep["PM2 bar, not recessive-only, 2015 table both sides"]
    ixp = {o["value"]: o for o in ix["points"]}
    p0, pdef = ixp["0"], next(o for o in ix["points"] if o["shipped"])
    gain = {t: pdef[t]["same_pathogenic"] - p0[t]["same_pathogenic"] for t in ("x37", "x38")}
    cost = {t: pdef[t]["false_pathogenic"] - p0[t]["false_pathogenic"] for t in ("x37", "x38")}
    curve = lambda t: ", ".join(f"{pc(o[t]['plp_recall'], 1)} at {'strict absence' if o['value'] == '0' else setting(o['value'])}" for o in ix["points"])  # noqa: E731
    dep37 = pd.read_csv(DATA / "pm2_dependence_v37.csv").set_index("rule")
    dep38 = pd.read_csv(DATA / "pm2_dependence_v38.csv").set_index("rule")
    tbl, pts_ = "Richards 2015 table", "shipped: SVI points on the pathogenic side"
    ixrun = P["v37_ix_table_pm2ad_0"]
    claims.append(Claim(
        ["Under the 2015 combining table, which is one switch away in the same classifier, the same sweep is steep", "Recall of pathogenic-side variants is 51.0% at strict absence", "The sweep that v10 recorded under this table", "The PM2 threshold moved pathogenic-side recall by 19 points when it was set in v10"],
        f"Under the 2015 combining table the same sweep is steep: recall of pathogenic-side variants is {curve('x37')}; raising the threshold from strict absence to 4e-5 gains {n(gain['x37'])} of 93,812 variants for {cost['x37']} more reference-discordant pathogenic calls, 'which is the knee that v10 recorded: the sweep then gave 37.8% at strict absence and 56.8% at 4e-5'; and the PM2 threshold moved pathogenic-side recall by 19 points when it was set in v10, by {pm_ad['shape']['v37']['range_gain'] / 93_812 * 100:.1f} in v37 and by {gain['x37'] / 93_812 * 100:.1f} again when the 2015 table is restored.",
        [f"Under v38 the recall is {curve('x38')}.",
         f"Raising the threshold from strict absence to 4e-5 gains {n(gain['x38'])} variants ({gain['x38'] / 93_812 * 100:.1f} points) for {cost['x38']} more reference-discordant pathogenic calls under v38, against {n(gain['x37'])} ({gain['x37'] / 93_812 * 100:.1f} points) for {cost['x37']} under v37.",
         f"v10's sweep (37.8% at strict absence and 56.8% at 4e-5) is reproduced by v37 ({pc(p0['x37']['plp_recall'], 1)} and {pc(pdef['x37']['plp_recall'], 1)}) and not by v38 ({pc(p0['x38']['plp_recall'], 1)} and {pc(pdef['x38']['plp_recall'], 1)}), so the manuscript says that v37 reproduced the sweep that v10 recorded and that v38 rises at every step as v37 does and by about a third as much (Section 4.10 below finds the shape the same step for step and the knee at 1e-5 under both), and Section 4.1 quotes {gain['x38'] / 93_812 * 100:.1f} points for v38 and {gain['x37'] / 93_812 * 100:.1f} for v37.",
         f"The sweep is still the steeper one: it moves recall {gain['x38'] / pm_ad['shape']['v38']['range_gain']:.0f} times as far as the same sweep under the default scheme ({n(gain['x38'])} against {n(pm_ad['shape']['v38']['range_gain'])} variants), where under v37 it moves it {gain['x37'] / pm_ad['shape']['v37']['range_gain']:.0f} times as far, and the statement that the combining rule decides how many calls depend on PM2 holds under v38: {n(dep38.loc[tbl, 'calls_depending_on_pm2'])} calls depend on PM2 under the 2015 table ({n(dep38.loc[tbl, 'of_which_pvs1_and_pm2_only'])} on PVS1 and PM2 alone) and {n(dep38.loc[pts_, 'calls_depending_on_pm2'])} under the default scheme (v37: {n(dep37.loc[tbl, 'calls_depending_on_pm2'])}, {n(dep37.loc[tbl, 'of_which_pvs1_and_pm2_only'])} and {n(dep37.loc[pts_, 'calls_depending_on_pm2'])}).",
         f"A reading that the firing counts support, and that no run tests (no run sweeps both PM2 thresholds under the 2015 table): at strict absence PM2 still fires on {n(ixrun.f38['PM2_Supporting'])} variants under v38 against {n(ixrun.f37['PM2_Supporting'])} under v37, because the {n(rec38)} genes that v38 reads as recessive-only keep their own threshold, which this sweep does not change."]))

    # ---- the BS2 carrier threshold
    bc = sweep["BS2 carrier bar, dominant-only genes"]
    bp = {o["value"]: o for o in bc["points"]}
    off, dflt, hi = bp["0"], next(o for o in bc["points"] if o["shipped"]), bp["0.001"]
    fb = lambda t: "/".join(str(o[t]["false_benign"]) for o in bc["points"])  # noqa: E731
    fp = lambda t: "/".join(str(o[t]["false_pathogenic"]) for o in bc["points"])  # noqa: E731
    claims.append(Claim(
        ["**The BS2 carrier threshold for genes read as dominant-only is a trade", "Switching it off, so that an allele count of at least five alone suffices", "On v37 the same sweep was flat, because the threshold reached the 50 genes"],
        f"In v37 the BS2 carrier threshold for genes read as dominant-only is flat, because it reached {int(float(f37['inheritance_read_dominant_only']))} genes: switching it off and setting it at 3e-3 change the numbers of reference-discordant benign and reference-discordant pathogenic calls by nothing ({dflt['x37']['false_benign']} and {dflt['x37']['false_pathogenic']} throughout), and the default 3e-4 recovers {n(dflt['x37']['same_pathogenic'] - off['x37']['same_pathogenic'])} pathogenic-side calls relative to no threshold ({n(off['x37']['same_pathogenic'])} to {n(dflt['x37']['same_pathogenic'])}) and gives up {n(off['x37']['same_benign'] - dflt['x37']['same_benign'])} benign-side calls ({n(off['x37']['same_benign'])} to {n(dflt['x37']['same_benign'])}), 'and nothing above it recovers another'.",
        [f"The manuscript says that the threshold is a trade on the {n(f38['inheritance_read_dominant_only'])} genes that v38 reads as dominant-only, where it was flat on the {int(float(f37['inheritance_read_dominant_only']))} that v37 read (line {find_line(manuscript_lines(), 'On v37 the same sweep was flat, because the threshold reached the 50 genes')}).",
         f"Under v38 the reference-discordant benign calls across the settings (off, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3) are {fb('x38')} and the reference-discordant pathogenic calls {fp('x38')} (under v37 {fb('x37')} and {fp('x37')}).",
         f"The default recovers {n(dflt['x38']['same_pathogenic'] - off['x38']['same_pathogenic'])} pathogenic-side calls relative to no threshold ({n(off['x38']['same_pathogenic'])} to {n(dflt['x38']['same_pathogenic'])}) and gives up {n(off['x38']['same_benign'] - dflt['x38']['same_benign'])} benign-side calls ({n(off['x38']['same_benign'])} to {n(dflt['x38']['same_benign'])}), against {n(dflt['x37']['same_pathogenic'] - off['x37']['same_pathogenic'])} and {n(off['x37']['same_benign'] - dflt['x37']['same_benign'])} under v37.",
         f"Above the default, 1e-3 recovers {n(hi['x38']['same_pathogenic'] - dflt['x38']['same_pathogenic'])} more pathogenic-side calls ({n(dflt['x38']['same_pathogenic'])} to {n(hi['x38']['same_pathogenic'])}), so 'nothing above it recovers another' does not hold under v38 (it holds under v37, where the count stays at {n(dflt['x37']['same_pathogenic'])})."]))

    # ---- the PVS1 regrade
    v37_pv = st.read_pv(BENCH / "manuscript_runs/v37_per_variant.tsv.gz").set_index(KEY)
    v38_pv = st.read_pv(BENCH / "manuscript_runs/v38_per_variant.tsv.gz").set_index(KEY)
    base_pv = st.read_pv(BENCH / "manuscript_runs/v38_base_per_variant.tsv.gz").set_index(KEY)
    g37, g38 = regrade_changes(v37_pv.reset_index()), regrade_changes(v38_pv.reset_index())
    assert (len(g37), len(g38)) == (int(f37["pvs1_regrade_calls_changed"]), int(f38["pvs1_regrade_calls_changed"])), "the replay here is not the one of 13_manuscript_stats.py"
    dropped = [k for k in g37 if k not in g38]
    has_pvs1 = lambda df, k: "PVS1" in df.loc[k].criteria  # noqa: E731
    lost_before = [k for k in dropped if not has_pvs1(base_pv, k)]
    lost_in_change = [k for k in dropped if has_pvs1(base_pv, k) and not has_pvs1(v38_pv, k)]
    other = [k for k in dropped if k not in lost_before and k not in lost_in_change]
    claims.append(Claim(
        ["The same replay on the recorded criteria of v37 changed 7"],
        f"Regrading the two PVS1 grades by one step would change {int(f37['pvs1_regrade_calls_changed'])} of 673,660 calls ({int(f37['pvs1_regrade_vus_to_lp'])} from VUS to Likely pathogenic and {int(f37['pvs1_regrade_b_to_vus'])} from Benign to VUS), by a replay that is an upper bound, because the recorded code does not say which branch gave a PVS1_Moderate; so few change because a grade one step lower plus PM2 seldom reaches 6 points either way.",
        [f"The same replay on the recorded criteria of v38 changes {int(f38['pvs1_regrade_calls_changed'])} calls ({int(f38['pvs1_regrade_vus_to_lp'])} from VUS to Likely pathogenic and {int(f38['pvs1_regrade_b_to_vus'])} from Benign to VUS), still an upper bound.",
         f"{len(dropped)} of the {len(g37)} variants that change under v37 do not under v38: {len(lost_before)} carry no PVS1 at all in the commit before the change, which differs from v37 in the criteria of 25 indels (Section 3.8, line {find_line(manuscript_lines(), 'Its run differs from v37')}), {len(lost_in_change)} lose it in the change itself" + (f", and {len(other)} keep it and no longer change call for another reason." if other else ", and none keeps PVS1 and stops changing."),
         "The reason that the manuscript gives for the small number holds a fortiori."]))

    # ---- the row dated on the day of the export
    gm37 = pd.read_csv(DATA / "clingen_missing_summary_v37.csv").set_index("measure")["n"]
    gm38 = pd.read_csv(DATA / "clingen_missing_summary_v38.csv").set_index("measure")["n"]
    key_day = [k for k in gm37.index if k.startswith("variants of GATM")][0]
    key_diff = [k for k in gm37.index if k.strip().startswith("whose call or criteria differ")][0]
    gdiff, gtags = [], {}
    for fam, pv_main in (("v37", v37_pv), ("v38", v38_pv)):
        gap = st.read_pv(BENCH / f"manuscript_runs/{fam}_gap_clingen_complete_per_variant.tsv.gz").set_index(KEY)
        both = pv_main[pv_main.symbol == "GATM"].join(gap[["call", "criteria", "clingen_tags"]], rsuffix="_gap")
        gdiff.append(both[(both.call != both.call_gap) | (both.criteria != both.criteria_gap)])
        gtags[fam] = (sorted(set(both.clingen_tags)), sorted(set(both.clingen_tags_gap)))
    d37, d38 = gdiff
    lost_codes = sorted({c_ for _, r_ in d38.iterrows() for c_ in set(r_.criteria.split("&")) - set(r_.criteria_gap.split("&"))})
    afs = (d38.gnomad_ac.astype(float) / d38.gnomad_an.astype(float))
    claims.append(Claim(
        ["and changes none of the 163 calls on GATM's variants", "It does change the criteria of 6 of them under v38"],
        f"The completed ClinGen file adds one row (GATM, Moderate) dated 04:00 UTC on 30 April, after the April file was saved, and it 'changes none of the {int(gm37[key_day])} calls on GATM's variants'.",
        [f"Under v38 the call of {int((d38.call != d38.call_gap).sum())} of the {int(gm38[key_day])} variants changes ({int((d37.call != d37.call_gap).sum())} under v37) and the call or criteria of {int(gm38[key_diff])} do ({int(gm37[key_diff])} under v37; `clingen_missing_summary_v38.csv`), so the manuscript says that the file changes none of the calls and then that it changes the criteria of {int(gm38[key_diff])}{'' if int((d38.call != d38.call_gap).sum()) == 0 else ' (the calls of some do change, and that sentence does not hold)'}.",
         f"The {len(d38)} variants lose {' and '.join(lost_codes)}; they are seen in {int(d38.gnomad_ac.astype(float).min())} to {int(d38.gnomad_ac.astype(float).max())} gnomAD alleles, at allele frequencies of {afs.min():.1e} to {afs.max():.1e}, between PM2's two thresholds (4e-5 for genes not read as recessive-only and 7e-5 for recessive-only genes, Section 3.3).",
         f"The established ClinGen tags of GATM's variants are {' and '.join(gtags['v38'][0])} in the main run and {' and '.join(gtags['v38'][1])} in the run that reads the completed file (the same under v37), so under v38, which reads the tag, the gene stops being read as recessive-only and PM2 takes the lower threshold for it."]))

    # ---- PVS1 on a ClinVar benign-side variant
    obg37 = pd.read_csv(DATA / "opposite_direction_by_group_v37.csv").set_index("group")
    obg38 = pd.read_csv(DATA / "opposite_direction_by_group_v38.csv").set_index("group")
    grp = "PVS1 on a ClinVar benign-side variant"
    rulings = pd.read_csv(DATA / "round6_rulings.tsv", sep="\t", dtype=str)
    gone = [r_ for r_ in rulings.itertuples() if r_.gene in set(str(f37["pvs1_benign_common_genes"]).split(",")) - set(str(f38["pvs1_benign_common_genes"]).split(","))]
    claims.append(Claim(
        ["Ten of these are a lone PVS1 and 11 are PVS1 with PM2", "Five are seen in 611 to 4,892 gnomAD alleles", "BS1, which tests the filtering allele frequency against 5e-3", "**PVS1 on a variant ClinVar calls benign or likely benign (21).**"],
        f"Section 3.6 groups {int(obg37.loc[grp, 'total'])} calls of PVS1 on a variant that ClinVar calls benign or likely benign: {int(f37['pvs1_benign_lone_pvs1'])} a lone PVS1 and {int(f37['pvs1_benign_pvs1_pm2'])} PVS1 with PM2; {int(f37['pvs1_benign_common_n'])} are seen in {n(f37['pvs1_benign_common_ac_min'])} to {n(f37['pvs1_benign_common_ac_max'])} gnomAD alleles ({str(f37['pvs1_benign_common_genes']).replace(',', ', ')}) with population-maximum frequencies of {float(f37['pvs1_benign_common_popmax_min']) * 100:.2f}% to {float(f37['pvs1_benign_common_popmax_max']) * 100:.1f}%; BS1 fired for none of them, their population-maximum frequencies being {float(f37['pvs1_benign_common_popmax_min']) / 2.5e-4:.0f} to {float(f37['pvs1_benign_common_popmax_max']) / 2.5e-4:.0f} times the median published threshold of 2.5e-4; the reviewer supported ClinVar on all six, 'remarking for five of them that the variant is too common in gnomAD'.",
        [f"Under v38 the group has {int(obg38.loc[grp, 'total'])} calls ({int(f38['pvs1_benign_lone_pvs1'])} a lone PVS1 and {int(f38['pvs1_benign_pvs1_pm2'])} PVS1 with PM2), because the RAD21 variant of Section 3 above left it: {int(f38['pvs1_benign_common_n'])} are seen in {n(f38['pvs1_benign_common_ac_min'])} to {n(f38['pvs1_benign_common_ac_max'])} gnomAD alleles ({str(f38['pvs1_benign_common_genes']).replace(',', ', ')}), at allele frequencies of {float(f38['pvs1_benign_common_af_min']) * 1e4:.1f}e-4 to {float(f38['pvs1_benign_common_af_max']) * 1e3:.1f}e-3 and population-maximum frequencies of {float(f38['pvs1_benign_common_popmax_min']) * 100:.2f}% to {float(f38['pvs1_benign_common_popmax_max']) * 100:.1f}%, {float(f38['pvs1_benign_common_popmax_min']) / 2.5e-4:.0f} to {float(f38['pvs1_benign_common_popmax_max']) / 2.5e-4:.0f} times the median published threshold, and BS1 fired for {'none' if int(f38['pvs1_benign_common_bs1_fired']) == 0 else str(int(f38['pvs1_benign_common_bs1_fired']))} of them.",
         f"In round 6 {ruling_text(gone[0].ruling, gone[0].rule)} on the RAD21 variant, so the manuscript says that the reviewer supported ClinVar on all {word(f38['pvs1_benign_common_n'])} (the RAD21 variant is {'the only one of the six that is no longer in the group' if len(gone) == 1 else 'one of those that are no longer in the group'}).",
         "The reviewer's remarks on how common each variant is in gnomAD are in the reviewer's own workbook and not in `round6_rulings.tsv`, so neither this script nor `verify_numbers.py` can check them.",
         "When the manuscript was written for v37 it said that the reviewer remarked for five of the six that the variant is too common in gnomAD.",
         "The manuscript now says four of the five (KDM2B and PHKB in those words, EDN3 and PIEZO1 as more than 1,000 carriers) and, for FMR1, its 210 hemizygotes and one homozygote and a question why the population frequencies were not applied.",
         "That is a reading of the returned workbook on the machine that made the paper, made when the manuscript was re-based, and the workbook is not in the repository."]))

    # ---- the reviewer's rulings
    tot37 = obg37[["clinvar", "fastvep", "neither", "needs_evidence", "disorder_specific"]].sum()
    tot38 = obg38[["clinvar", "fastvep", "neither", "needs_evidence", "disorder_specific"]].sum()
    r6a = pd.read_csv(DATA / "round6_rows_v37.csv")
    r6b = pd.read_csv(DATA / "round6_rows_v38.csv")
    rl37 = pd.read_csv(DATA / "round6_rules_v37.csv")
    rl38 = pd.read_csv(DATA / "round6_rules_v38.csv")
    kind = lambda rl, label, rule: int(rl[(rl.set == label) & (rl.rule == rule) & (rl.ruling == "clinvar")].rows.sum())  # noqa: E731
    lab37, lab38 = f"the {int(r6a.still_opposite.sum())} still opposite in v37", f"the {int(r6b.still_opposite.sum())} still opposite in v38"
    sup = lambda rows, side: int(rows[rows.clinvar_class.isin(side) & rows.still_opposite & (rows.ruling == "clinvar")].shape[0])  # noqa: E731
    claims.append(Claim(
        ["the reviewer supported ClinVar's class for 13 of the 24 reference-discordant benign calls", "For 29 of the 48 (60%) the reviewer supported ClinVar's call", "So for 29 of the 48 the reviewer judged ClinVar's call the better", "Of the 29 rulings for ClinVar, 10 follow from a class", "In the last round the reviewer supported ClinVar on 60% of the 48"],
        f"For {int(tot37.clinvar)} of the {int(tot37.sum())} ({pc(tot37.clinvar / tot37.sum(), 0)}) opposite-direction calls the reviewer supported ClinVar's call ({sup(r6a, ('Pathogenic', 'Likely_pathogenic'))} of the {int(h37['false_benign'])} reference-discordant benign calls and {sup(r6a, ('Benign', 'Likely_benign'))} of the {int(h37['false_pathogenic'])} reference-discordant pathogenic calls), for {int(tot37.fastvep)} fastVEP's, for {int(tot37.neither)} the reviewer's own call was VUS, for {int(tot37.needs_evidence)} more evidence was needed and for {int(tot37.disorder_specific)} the answer depends on the disorder; of the {int(tot37.clinvar)} rulings for ClinVar, {kind(rl37, lab37, 'class')} follow from a class that the reviewer stated, {kind(rl37, lab37, 'stated')} from the reviewer's statement that ClinVar is correct and {kind(rl37, lab37, 'evidence')} from evidence for ClinVar's class without naming one.",
        [f"Under v38 the reviewer's rulings on the {int(tot38.sum())} calls give ClinVar {int(tot38.clinvar)} ({pc(tot38.clinvar / tot38.sum(), 0)}; {sup(r6b, ('Pathogenic', 'Likely_pathogenic'))} of the {int(h38['false_benign'])} reference-discordant benign calls and {sup(r6b, ('Benign', 'Likely_benign'))} of the {int(h38['false_pathogenic'])} reference-discordant pathogenic calls), fastVEP {int(tot38.fastvep)}, neither {int(tot38.neither)}, more evidence {int(tot38.needs_evidence)} and the disorder {int(tot38.disorder_specific)}; the rulings for ClinVar are {kind(rl38, lab38, 'class')} from a class, {kind(rl38, lab38, 'stated')} from a statement and {kind(rl38, lab38, 'evidence')} from evidence (`round6_rules_v38.csv`).",
         f"{int((~r6b.still_opposite).sum())} of the 64 rows of round 6 are no longer opposite under v38 ({int((~r6a.still_opposite).sum())} under v37), and the reviewer supported ClinVar's call on {int(((~r6b.still_opposite) & (r6b.ruling == 'clinvar')).sum())} of them ({int(((~r6a.still_opposite) & (r6a.ruling == 'clinvar')).sum())} under v37); the one that was added is the RAD21 variant, which the reviewer first saw in {rulings[(rulings.gene == 'RAD21') & (rulings.pos == '116866762')].iloc[0].seen_in_round}.",
         f"The calls of the 48 were first seen by the reviewer in round 3 for {int((r6b[r6b.still_opposite].seen_in_round == 'round 3').sum())}, round 4 for {int((r6b[r6b.still_opposite].seen_in_round == 'round 4').sum())} and round 5 for {int((r6b[r6b.still_opposite].seen_in_round == 'round 5').sum())} (v37: {int((r6a[r6a.still_opposite].seen_in_round == 'round 3').sum())}, {int((r6a[r6a.still_opposite].seen_in_round == 'round 4').sum())} and {int((r6a[r6a.still_opposite].seen_in_round == 'round 5').sum())})."]))

    # ---- the new variants, which bands carry what the mix does not explain
    b37, b38 = (pd.read_csv(DATA / f"temporal_recall_by_af_band_{b}.csv") for b in ("v37", "v38"))

    def contributions(bands: pd.DataFrame) -> pd.DataFrame:
        bb = bands[bands.clinvar == "Benign"].copy()
        bb["w"] = bb.n_new / bb.n_new.sum()
        bb["c"] = bb.w * (bb.recall_new - bb.recall_april) * 100
        return bb.set_index("popmax_af_band")

    c37_, c38_ = contributions(b37), contributions(b38)
    top = ["2e-2 to 5e-2", ">=5e-2"]
    t37, t38 = -c37_.loc[top].c.sum(), -c38_.loc[top].c.sum()
    s37, s38 = -c37_.c.sum(), -c38_.c.sum()
    rows = [[band, n(c37_.loc[band].n_new), pc(c37_.loc[band].recall_april, 1), pc(c37_.loc[band].recall_new, 1), f"{-c37_.loc[band].c:+.2f}", pc(c38_.loc[band].recall_april, 1), pc(c38_.loc[band].recall_new, 1), f"{-c38_.loc[band].c:+.2f}"] for band in c37_.index]
    table = ["| Band of population-maximum frequency | New Benign variants | April recall v37 | New recall v37 | Points v37 | April recall v38 | New recall v38 | Points v38 |", "|:---|---:|---:|---:|---:|---:|---:|---:|"] + ["| " + " | ".join(r_) + " |" for r_ in rows]
    cl = Claim(
        ["What remains is spread over the bands"],
        f"What remains of the fall in recall of the new Benign variants, once the mix of nine bands of population-maximum frequency is taken as the April recall in each band, 'is mostly in the highest bands (the two at 2% or more carry {t37:.1f} of the {s37:.1f} points that the mix does not explain)'.",
        [f"Under v38 the two bands at 2% or more carry {t38:.1f} of the {s38:.1f} points ({t38 / s38 * 100:.0f}%) where under v37 they carry {t37:.1f} of {s37:.1f} ({t37 / s37 * 100:.0f}%), so the manuscript no longer says 'mostly in the highest bands' and says that what remains is spread over the bands; the unexplained part is larger ({s38:.1f} points against {s37:.1f}).",
         "The positive and negative contributions of each band to the points that the mix does not explain (the points are the band's share of the new Benign variants times the April recall less the new recall, so a negative figure means the new variants are called benign-side more often than the April ones of that band):"])
    cl.table = table
    claims.append(cl)

    # ---- the annotation time and the crate
    prov89 = pd.read_csv(DATA / "provenance_v38.csv")
    prov37 = pd.read_csv(DATA / "provenance_v37.csv")
    claims.append(Claim(
        ["takes 25.3 seconds to annotate on every core", "the annotation of one of the 89 runs of v37 took a median of"],
        "A benchmark of 673,660 ClinVar variants of review status two stars or better 'takes about two minutes to annotate in the runs of v37 (a median of 127 seconds over its 89 runs on the April set, made on 4 November 2026 up to three at a time on one Apple M4 machine with 10 cores and 24 GB that another job was also using) and about ten seconds to score', and release 0.5.0 'was timed alone on the same machine on 7 November 2026: it took 25.3 seconds on every core (the median of three runs), 74 seconds on one thread (one run) and 11 seconds to score'.",
        [f"The {len(prov37)} April runs of v37 made for this paper took a median of {int(prov37.annotation_seconds.median())} seconds ({int(prov37.annotation_seconds.min())} to {int(prov37.annotation_seconds.max())}) and the {len(prov89)} of v38 a median of {int(prov89.annotation_seconds.median())} ({int(prov89.annotation_seconds.min())} to {int(prov89.annotation_seconds.max())}), made on 8 November 2026 several at a time on the same machine with other work beside them, with the complete databases (which hold 8% more records) and by binaries that contain the two speed-ups of 5 November 2026 for v38 only.",
         "The manuscript therefore says that the times of these runs are not timings, and keeps the time of release 0.5.0 measured alone on 7 November 2026 on the first databases (`results/v38_v0.5.0/annotation_time.tsv`), which was not repeated with the complete databases (Section 4.4)."]))
    l37, t37_, v37_ver = crate_size("5d8004a")
    l38, t38_, v38_ver = crate_size("v0.5.0-rc.1")
    claims.append(Claim(
        ["The classification engine is the `fastvep-classification` crate", "v37 (5d8004a) is the version in which the review rounds ended"],
        f"The classification engine is the `fastvep-classification` crate of fastVEP (version 0.4.0): 14,857 lines of Rust with 340 unit tests at the commit analyzed.",
        [f"The method that gives these figures for 5d8004a (`src` files of the crate, `#[test]` attributes) gives {n(l37)} lines and {t37_} tests there, version {v37_ver}, and {n(l38)} lines and {t38_} tests at the tag `v0.5.0-rc.1` (c12d9f9), version {v38_ver}, and the manuscript gives both."]))

    # ---- the veto and the ClinVar-derived inputs
    def false_benign_set(path: Path) -> set:
        d = opposite_of(path)
        d = d[d.clinvar.isin(("Pathogenic", "Likely_pathogenic"))]
        return {tuple(r_) for r_ in d[KEY].itertuples(index=False)}

    sets = {}
    for fam, base_dir in (("v37", BENCH / "output_v37"), ("v38", BENCH / "manuscript_runs/v38")):
        b0 = false_benign_set(base_dir)
        sets[fam] = (len(b0), len(false_benign_set(BENCH / f"manuscript_runs/{fam}_sw_no_bp1_veto") - b0), len(false_benign_set(BENCH / f"manuscript_runs/{fam}_src_noclinvar_protein") - b0),
                     len(false_benign_set(BENCH / f"manuscript_runs/{fam}_src_noclinvar_all") - b0),
                     (false_benign_set(BENCH / f"manuscript_runs/{fam}_sw_no_bp1_veto") - b0) == (false_benign_set(BENCH / f"manuscript_runs/{fam}_src_noclinvar_protein") - b0))
    claims.append(Claim(
        ["The dependence on ClinVar-derived inputs is therefore small in agreement and larger in errors", "Removing every ClinVar-derived input costs 0.78 points"],
        f"Without the ClinVar protein and splice index, reference-discordant benign calls rise from {sets['v37'][0]} to {sets['v37'][0] + sets['v37'][2]}, and from {sets['v37'][0]} to {sets['v37'][0] + sets['v37'][3]} with the allele-level record as well, and 'the {sets['v37'][2]} are the same variants as appear when BP1's veto is switched off'; removing every ClinVar-derived input 'takes the opposite-direction calls from 49 to 74, {sets['v37'][2]} of the reference-discordant benign calls because it removes the BP1 veto'.",
        [f"Under v38 the index removal adds {sets['v38'][2]} reference-discordant benign calls (from {sets['v38'][0]} to {sets['v38'][0] + sets['v38'][2]}) and removing the allele-level record as well adds {sets['v38'][3]} (to {sets['v38'][0] + sets['v38'][3]}), and the {sets['v38'][2]} {'are' if sets['v38'][4] else 'are not'} the same variants as appear when BP1's veto is switched off ({sets['v38'][1]} added, {'the same set' if sets['v38'][4] else 'a different set'}), so the manuscript says {sets['v38'][2]} where the sentence for v37 had {sets['v37'][2]}; the allele-level record accounts for {sets['v38'][3] - sets['v38'][2]} of the additions under v38 where it accounts for {sets['v37'][3] - sets['v37'][2]} under v37."]))
    ms = manuscript_lines()
    claims.sort(key=lambda c_: min(find_line(ms, a_) for a_ in c_.anchors))
    return claims


def section_rewording(D: Doc, pairs: list[Pair], main: Pair, october: Pair, summary: list[dict], cmp_: dict, changes: tuple, claims: list, outs: dict) -> None:
    ms = manuscript_lines()
    heads = headings(ms)
    numbers, n_checks, n_values = changes
    D.h(2, "9. Where the manuscript's statements differ between v38 and v37")
    D.s("The manuscript says v38 and quotes v37 only as the history of the numbered versions and as the before of a difference between the two.")
    D.s("Its line numbers here are those of `manuscript/acmg/manuscript.md` as it stands, and this section lists what would read differently if it were written for v37.")
    D.s("Three kinds of statement differ: a number (the same sentence, another value), a claim whose truth or wording differs, and a statement about how the runs were made.")
    D.s("The first kind is listed in full from `verify_numbers.py`, which recomputes every number of the text from the tables: with `--base v38` (the default) it recomputes them from the twin tables of v38, and with `--base v37` from those of v37, and `--dump` writes the line, the anchor and the value of each check, so the dumps of the two bases, set side by side, say which sentence differs and in what.")
    D.s("The second kind is listed by hand from a reading of the sentences around those numbers, each with the numbers that decide it, computed here from the scored files.")
    D.s("The third kind is listed from the provenance of the runs.")
    D.s("The sentences that say what v37 read (the inheritance read from the words of the disease names: 67 genes recessive-only, 50 dominant-only and 12 both, and the genes the carrier threshold reached) and the trajectory of the numbered versions (`trajectory.csv`, `change_ledger.csv`, `review_rounds.csv`) are the same for both and are not listed as differing.")

    D.h(3, "9.1 Claims whose truth or wording differs")
    for cl in claims:
        lines_ = sorted({find_line(ms, a) for a in cl.anchors})
        where = ", ".join(dict.fromkeys(section_of(x, heads) for x in lines_))
        D.s(f"**{'Line' if len(lines_) == 1 else 'Lines'} {', '.join(str(x) for x in lines_)} ({where}).**")
        D.s(f"Written for v37, the statement read: {cl.says}")
        for e in cl.evidence[:-1] if hasattr(cl, "table") else cl.evidence:
            D.s(e)
        if hasattr(cl, "table"):
            D.s(cl.evidence[-1])
            D.blank()
            D.out += cl.table
        D.blank()

    D.h(3, "9.2 Statements about how the runs were made and what was measured on v38")
    family = [p_ for p_ in pairs if p_.dir37.parent == BENCH / "manuscript_runs"]
    kinds = {"plain": [], "cutoff": [], "hook": []}
    for p_ in family:
        if p_.kind == "april":
            kinds["cutoff" if p_.group == "cutoff" else "hook" if p_.group == "hook" else "plain"].append(prov(p_.v38)["fastvep_sha256"][:8])
    ca = {k: dict(Counter(v)) for k, v in kinds.items()}
    k37 = {k: dict(Counter(prov(p_.v37)["fastvep_sha256"][:8] for p_ in family if p_.kind == "april" and (p_.group == "cutoff") == (k == "cutoff") and (p_.group == "hook") == (k == "hook"))) for k in ("plain", "cutoff", "hook")}
    shared = [p_ for p_ in family if is_shared(p_.v38)]
    new_apr = [p_ for p_ in family if p_.kind == "april" and not is_shared(p_.v38)]
    new_oct = [p_ for p_ in family if p_.kind == "october"]
    all_files = sorted(f.name for f in DATA.iterdir() if f.is_file())
    section38 = [f for f in all_files if f.startswith("v38_")]
    twins = [f for f in all_files if "_v38" in f and not f.startswith("v38_")]
    v37_files = [f for f in all_files if not f.startswith("v38_") and "_v38" not in f]
    no_twin = [f for f in v37_files if "_v37" not in f or f.replace("_v37", "_v38") not in all_files]
    history = [f for f in no_twin if f in ("trajectory.csv", "change_ledger.csv", "review_rounds.csv")]
    independent = {"clingen_missing_curations.csv": "a list of the rows of the ClinGen export", "criteria_status.csv": "read from the engine of v38 already", "tavtigian_posterior.csv": "arithmetic of the framework, which `13_manuscript_stats.py` checks to be the file that is there",
                   **{f: "a comparison of the first and the complete databases (`21_compare_databases.py`), which holds both versions" for f in no_twin if f.startswith("databases_")}}
    assert set(no_twin) == set(history) | {"round6_rulings.tsv"} | set(independent), f"the files without a twin are not the ones that Section 9.2 describes: {no_twin}"

    def where(anchor: str, inline: bool = False) -> str:
        line = find_line(ms, anchor)
        return f"{'line' if inline else 'Line'} {line} ({section_of(line, heads)})"

    D.blank()
    D.s(f"**{where('The analyses of Sections 3.1 to 3.7 were made on v38 and on v37')} and {where('Every number of this paper is v38', inline=True)}.**")
    D.s("The manuscript says that every number is v38's unless it names v37, and that the analyses of Sections 3.1 to 3.7 were made on both versions, except for what needs a reviewer or a run that the repository does not have.")
    D.s(f"Written for v37 it said that the analyses of Sections 3.1 to 3.7 had not been repeated on v38, except for the four thresholds that read inheritance and the runs on gnomAD's joint release (Section 3.8).")
    D.s(f"The family v38 repeats the {len(family)} runs of Section 2.6 ({len([p_ for p_ in family if p_.kind == 'april'])} on the April set and {len(new_oct)} on the October release) with the classifier of release 0.5.0, and `13_manuscript_stats.py --base v38` and `14_temporal_holdout.py --base v38` write {len(twins)} twin tables of `manuscript/acmg/data/` from them (`_v38` where the file of v37 has `_v37` in its name), which is what lets the manuscript say that Sections 3.1 to 3.7 were repeated on v38.")
    why_no = [f"{', '.join(f'`{f}`' for f in history)} are the history of the numbered versions", "`round6_rulings.tsv` is the reviewer's own ruling"] + [f"`{f}` is {why}" for f, why in independent.items()]
    D.s(f"{len(no_twin)} of the {len(v37_files)} files of `manuscript/acmg/data/` that are not v38's own have no twin: {'; '.join(why_no)}.")
    D.s("What was not repeated on v38 is what needs a reviewer or a run that the repository does not have: the review rounds, the ledger, and the runs of round 6, which need builds that are not in the repository.")

    D.blank()
    D.s(f"**{where('All 89 runs (Supplementary Table S2) were made by the same harness')}.**")
    D.s(f"Written for v37, the statement was that all 89 runs were made by the same harness, {sum(k37['plain'].values())} with the binary of v37 (SHA-256 `{next(iter(k37['plain']))}...`), {sum(k37['cutoff'].values())} with the build of the first patch and {sum(k37['hook'].values())} with the build of the second.")
    D.s(f"For the family v38 the split of the 89 April runs is the same, {sum(ca['plain'].values())} plain, {sum(ca['cutoff'].values())} with the cutoff build and {sum(ca['hook'].values())} with the hook build, and each kind has one binary: {', '.join(f'`{k}...` made {v}' for k, v in sorted(ca['plain'].items(), key=lambda kv: -kv[1]))} for the plain runs, `{', '.join(ca['cutoff'])}...` for the cutoff runs and `{', '.join(ca['hook'])}...` for the hook runs (`binaries_v38.tsv`).")
    D.s("`verify_numbers.py` counts a run of another binary than its kind's as a failure, and the 42 runs of Section 3.8 must have been made a second time by the binary of the stored run, with every scored file identical (`results/v38_v0.5.0`).")

    D.blank()
    D.s(f"**{where('The three binaries are builds of the tag')}.**")
    D.s("Written for v37, the statement was that the three binaries are builds of commit 5d8004a, each from a clean export in a target directory of its own (`build_manuscript_binaries.sh`).")
    D.s("The binaries of the family v38 are builds of the tag `v0.5.0-rc.1` (commit c12d9f9), made the same way (`build_manuscript_binaries.sh`, names `v38`, `cut38` and `hook38`); the two patches apply to the tag with every context line matching, the cutoff patch at the line it was written for and the hooks patch with offsets of 2 and 20 lines, the lines that the tree of the tag has gained above its two hunks.")

    D.blank()
    D.s(f"**{where('The 38 runs of Section 3.8 on the April set and the 4 on the October release were then made a second time')}, {where('The same 96 runs were made for v37', inline=True)}, {where('All 42 were made on 8 November 2026', inline=True)}.**")
    D.s("Written for v37, the statement was that every run was made twice, in October 2026 and again on 4 November 2026, and that the two sets agree in every scored file.")
    D.s(f"In the manuscript now every run of the family was made once on {new_york_days([prov(p_.v38)['recorded_at'] for p_ in family])} with the complete databases, and the 42 runs of Section 3.8 were made a second time by the same binaries into a directory of their own, identical in every scored file (`results/v38_v0.5.0/rerun_comparison.csv`), so the check that the others are right is the controls (Section 6 above), which reproduce the main run in every cell and in the firing count of every criterion.")
    commits = Counter(prov(p_.v38)["repo_commit"][:7] for p_ in family)
    D.s(f"The provenance file of each run records the commit of the checkout that ran the harness ({', '.join(f'`{c}` for {k}' for c, k in commits.most_common())}), which is not the commit of the binary, as the manuscript and `results/v37_manuscript/README.md` say for both families.")

    D.blank()
    D.s(f"**{where('The prices of all seven steps are measured on v38')}.**")
    D.s("Written for v37, the statement was that the prices of the first six steps of Section 4.5 were measured on v37 and not repeated on v38, and that the seventh was measured on v38.")
    D.s(f"All seven are measured on v38 in the family (`temporal_proposals_v38.csv`, the April and the October columns), and the numbers of the steps are in 9.4 below (lines {', '.join(str(find_line(ms, a_)) for a_ in ('First, complete the ClinGen input', 'Third, BP7 as the splicing recommendation describes it', 'Fourth, BP1', 'Sixth, review of the allele-number floor'))}).")

    D.blank()
    D.s(f"**{where('The scored concordance matrix, the criterion firing counts and the provenance file of each run made for this paper')}, {where('The results of v38, which is the version of every number', inline=True)}.**")
    D.s("Written for v37, the statement listed the results of v37 as the 89 runs on the April set and the 7 on the October release in `analysis/acmg_benchmark/results/v37_manuscript/` with `binaries.tsv`, whose 74 files in `manuscript/acmg/data/` were the ones whose names do not begin with `v38_`, and the results of v38 as the 42 runs of Section 3.8, their repetition by release 0.5.0 and the 23 files of `manuscript/acmg/data/` whose names begin with `v38_`.")
    D.s(f"The manuscript now lists the results of v38 first and those of v37 as the before: the family v38 adds {len(new_apr)} runs on the April set and {len(new_oct)} on the October release to `runs_v38.tsv`, `holdout_runs_v38.tsv` and `run_matrix_v38.tsv` of the same directory (the {len(shared)} shared runs were already there), with their provenance files and configurations, `binaries_v38.tsv`, the {len(twins)} twin tables whose names have `_v38` in them, each with its twin for v37 ({'the 67 pairs that the manuscript counts' if len(twins) == 67 else 'not the number of pairs that the manuscript counts'}), the {len(section38)} files that begin with `v38_` ({'the 23 that the manuscript counts' if len(section38) == 23 else 'not the number that the manuscript counts'}) and the {len(no_twin)} that belong to neither version ({'the seven that the manuscript counts' if len(no_twin) == 7 else 'not the number that the manuscript counts'}), this comparison with its five CSV files, and the per-variant tables of the new runs in `data/benchmark/manuscript_runs/` (not in the repository).")

    D.h(3, "9.3 The tables, the workbook and the figures")
    D.s(f"`build_tables.py --base v37 --check` finds {len(outs['tables'])} of the {outs['n_tables']} marked tables of the manuscript different from the text, which says v38: {', '.join(outs['tables'])}.")
    caption, in_section = {}, {}
    for i_, line_ in enumerate(ms):
        m_ = re.match(r"<!-- table:(\w+):begin -->", line_)
        if m_:
            caption[m_.group(1)] = next((re.match(r"\*\*(Table \d+)\.", ms[k]).group(1) for k in range(i_, -1, -1) if re.match(r"\*\*Table \d+\.", ms[k])), "")
            in_section[m_.group(1)] = section_of(i_ + 1, heads)
    unchanged = [t for t in caption if t not in outs["tables"]]
    D.s(f"They are {', '.join(f'{caption[t]} (`{t}`)' for t in outs['tables'])}; the ones that do not differ are {', '.join(f'{caption[t]} (`{t}`, in {in_section[t]})' for t in unchanged)}.")
    D.s("`build_tables.py --base v37 --write-to PATH` writes the manuscript with the tables of v37 to a file and leaves `manuscript.md` as it is.")
    for row, u, v in outs["words"]:
        D.s(f"The words of a cell of Table 3 differ and not only its numbers, in the row '{row}': under v37 it reads '{u}' and under v38, as in the manuscript, '{v}'.")
    n_sheets, differ, alike = outs["sheets"]
    D.s(f"The workbook of v37 (`build_supplement.py --base v37 --write --out PATH`) has the same {n_sheets} sheets, and {len(differ)} of them differ from `supplementary_tables.xlsx`: {', '.join(differ)}; the same are {', '.join(alike)}.")
    n_png, differ_png, same_png = outs["figures"]
    D.s(f"`generate_figures.py --base v37 --out DIR` draws the same {n_png} figures as PNG files, and {len(differ_png)} of them differ from the ones in `manuscript/acmg/figures/`: {', '.join(differ_png)}; {', '.join(same_png)} are the same, being the history of the numbered versions, which ends with v38 under both.")
    D.s("The PDF files of the figures are not tracked and carry the time they were made, so they are not compared.")

    D.h(3, "9.4 The numbers of the text that differ")
    changed_values = sum(len(c[2]) for c in numbers)
    D.s(f"`verify_numbers.py` makes {n_checks} checks of {n_values} values; under v37 {len(numbers)} of the checks give another value than under v38, {changed_values} values in all (the sentences that set v38 against v37 contain both, so not every one of these is absent from its sentence).")
    D.s("They are listed below by section of the manuscript, each as the line, the phrase that identifies the sentence, and the values that differ, the value under v37 first and the value under v38, which the manuscript has, second.")
    D.s("A value that is a statement about v37 as history (what v37 read, the trajectory of the numbered versions) does not appear, and a quoted count of `48` is the count of v38 (v37: 49), which the guard reads from the main run.")
    D.s("A value written `?` is a condition that the guard types (`Ten`, `Five`, `none of the`) and that does not hold under v37; the claims above explain each.")
    D.s("The numbers inside the marked tables are not checks of `verify_numbers.py` and are in 9.3.")
    by_section: dict[str, list] = {}
    for line_, anchor, d in numbers:
        by_section.setdefault(section_of(line_, heads), []).append((line_, anchor, d))
    order = list(dict.fromkeys(lbl for _, lbl in heads))
    for sec in order:
        if sec not in by_section:
            continue
        D.blank()
        D.s(f"**{sec[0].upper() + sec[1:]}.**")
        D.bullets([f"Line {line_}, \"{short(anchor)}\": " + "; ".join(f"{('a condition that does not hold' if u == '?' else u)} under v37 becomes {v}" for u, v in d) + "." for line_, anchor, d in sorted(by_section[sec])])


# ------------------------------------------------------------------------------------------------------ main


def check_typed_claims(summary: list[dict], cmp_: dict, eff: pd.DataFrame) -> None:
    """The few statements of the text above that are typed and not computed. They hold for these data, and if the data change the script stops
    here and does not write a document that says something the data no longer say."""
    knees = {s_["label"]: (s_["shape"]["v37"]["knee"], s_["shape"]["v38"]["knee"]) for s_ in summary if s_["shape"]["v37"]["knee"] is not None or s_["shape"]["v38"]["knee"] is not None}
    # the PM2 sweep under the 2015 table has the knee at 1e-5 under v38 only: its exchange rate falls from 880 calls per error to 84 under v38, and from 3,410
    # to 484 under v37, which is a fall by a factor of 7 and not to about a hundred (the first form of the rule, which asked only for the factor, gave it to both)
    expected = {"BS1 bar": (0.001, 0.001), "BP7 far boundary": (300.0, 300.0), "PM2 bar, not recessive-only, 2015 table both sides": (None, 1e-05), "Likely pathogenic cutoff (points)": (5.0, 5.0)}
    assert knees == expected, f"the knees are not the ones that Section 1 of the document names: {knees}"
    sw = {s_["label"]: s_["shape"] for s_ in summary}
    for label in ("BS2 carrier bar, dominant-only genes", "BS2 carrier floor, bar off"):
        assert sw[label]["v37"]["range_false_benign"] == 0 and sw[label]["v38"]["range_false_benign"] > 0, f"{label} is not flat under v37 and not under v38"
    ix = sw["PM2 bar, not recessive-only, 2015 table both sides"]
    assert ix["v38"]["range_gain"] < ix["v37"]["range_gain"] / 2, "the PM2 sweep under the 2015 table keeps most of its range under v38"
    flagged = eff[eff.group.isin(["switch", "source", "hook", "gap"]) & ((eff.changed_sign == True) | (eff.changed_by_more_than_20_percent == True) | (eff.appeared_or_vanished == True)) & (eff.small_number == False)]  # noqa: E712
    big = flagged[flagged.effect_delta.abs() >= 100]
    assert list(big.run_v37) == ["v37_src_noclingen"], f"the one large change of an effect is not the removal of the ClinGen file: {list(big.run_v37)}"
    x = cmp_["x"]
    ch = x[(x.call != x.call_38) & (x.cv3 == "P/LP") & (x.call_38 == "VUS") & x.call.isin(["Likely_pathogenic", "Pathogenic"])]
    assert set(ch.symbol) == {"NEK1", "DDX41", "SRCAP", "PROKR2"}, f"the pathogenic-side variants that lose their call are not the four of Section 3.8: {sorted(ch.symbol)}"


def section_summary(D: Doc, pairs: list[Pair], main_pair: Pair, october: Pair, summary: list[dict], eff: pd.DataFrame, opp: pd.DataFrame, counts: dict, changes: tuple, n_claims: int, outs: dict) -> None:
    a, b = main_pair.x37, main_pair.x38
    N = APRIL_N
    oa, ob = october.x37, october.x38
    numbers, n_checks, n_values = changes
    group_eff = eff[eff.group.isin(["switch", "source", "hook", "gap"])]
    flagged = group_eff[(group_eff.changed_sign == True) | (group_eff.changed_by_more_than_20_percent == True) | (group_eff.appeared_or_vanished == True)]  # noqa: E712
    others = {k: v for k, v in counts.items() if k not in ("v37", "v37_oct2026")}
    more = [(k, c37, c38) for k, (c37, c38, _, _) in others.items() if c38 > c37]
    fewer = [(k, c37, c38) for k, (c37, c38, _, _) in others.items() if c38 < c37]
    left = opp[opp.change == "left_in_v38"].drop_duplicates(["gene", "chrom", "pos", "ref", "alt"])
    entered = opp[opp.change == "entered_in_v38"].drop_duplicates(["gene", "chrom", "pos", "ref", "alt"])
    changed_shape = shape_changes(summary)
    moved = [c for c in set(main_pair.f37) | set(main_pair.f38) if main_pair.f37.get(c, 0) and abs(main_pair.f38.get(c, 0) - main_pair.f37[c]) / main_pair.f37[c] > 0.05]
    ta, tb = pd.read_csv(DATA / "temporal_v37.csv"), pd.read_csv(DATA / "temporal_v38.csv")
    al, nw = ta[ta.stratum == "all"].iloc[0], tb[tb.stratum == "all"].iloc[0]
    na, nb = ta[ta.stratum == "new"].iloc[0], tb[tb.stratum == "new"].iloc[0]
    ra, rb = pd.read_csv(DATA / "temporal_ratio_v37.csv"), pd.read_csv(DATA / "temporal_ratio_v38.csv")
    xa, xb = ra[ra.measure.str.startswith("ratio")].iloc[0], rb[rb.measure.str.startswith("ratio")].iloc[0]
    D.out += ["# The runs of the family v38 against the runs of the family v37", ""]
    D.s("This file sets the classifier of release 0.5.0 (v38), the latest, against v37, the classifier at the close of the review rounds, which the ACMG manuscript described until it was re-based on v38, run for run, on the same inputs and by the same harness.")
    D.s("It is written by `analysis/acmg_benchmark/scripts/20_compare_families.py` from the scored files of the runs and the tables of `manuscript/acmg/data/`, and every number in it is read from them.")
    D.s("ClinVar's 2-star+ classifications are a reference and not ground truth, and every agreement here is with them and nothing stronger.")
    D.h(2, "Summary")
    D.bullets([
        f"On the main run, v38 agrees with ClinVar's class for {n(b['exact'] - a['exact'])} more of the {n(N)} records and with its direction for {n(b['same_direction'] - a['same_direction'])} more ({pc(a['same_direction'] / N)} to {pc(b['same_direction'] / N)}); it calls {n(b['same_benign'] - a['same_benign'])} more ClinVar benign-side variants benign-side (recall {pc(a['blb_recall'])} to {pc(b['blb_recall'])}) and {n(abs(b['same_pathogenic'] - a['same_pathogenic']))} fewer pathogenic-side variants pathogenic-side (recall {pc(a['plp_recall'])} to {pc(b['plp_recall'])}) (Section 2).",
        f"Opposite-direction calls go from {a['opposite']} to {b['opposite']} ({a['false_benign']} reference-discordant benign calls and {a['false_pathogenic']} reference-discordant pathogenic calls, then {b['false_benign']} and {b['false_pathogenic']}): {len(opp[(opp.pair_v37 == 'v37') & (opp.change == 'left_in_v38')])} call leaves the set and {len(opp[(opp.pair_v37 == 'v37') & (opp.change == 'entered_in_v38')])} enter it (Section 3).",
        f"Across the {len(others)} other pairs with a discrepancy file on both sides, {len(fewer)} have fewer opposite-direction calls under v38, {len(others) - len(fewer) - len(more)} the same number and {len(more)} more, and over all the pairs {len(left)} variants leave and {len(entered)} enter in at least one; the pairs with more are {', '.join(f'`{k}` ({c37} to {c38})' for k, c37, c38 in more)} (Section 3).",
        f"Of the {len(summary)} sweeps, {len(changed_shape)} differ in the sign of a step, the position of the knee or the order of their settings, and {plural(sum(1 for s_ in summary if s_['shape']['v37']['knee'] != s_['shape']['v38']['knee']), 'differs', 'differ')} in the knee; the sweeps of the BS2 carrier bar and floor stop being flat and the PM2 sweep under the 2015 table loses most of its range (Section 4).",
        f"Of the {len(group_eff)} effects of the single switches, source removals, hook runs and gap runs on {word(int(group_eff.metric.nunique()))} measures, {len(flagged)} changed sign, appeared or vanished, or changed by more than 20%, and {int((flagged.small_number == True).sum())} of them are of fewer than 10 calls; the one that is large is the removal of the ClinGen file (Section 5).",  # noqa: E712
        f"In the main run {len(moved)} criterion moves its firing count by more than 5% ({', '.join(sorted(moved))}); PM2 at Supporting strength moves by {(main_pair.f38['PM2_Supporting'] - main_pair.f37['PM2_Supporting']) / main_pair.f37['PM2_Supporting'] * 100:+.1f}% (Section 7).",
        f"On the October release agreement in direction goes from {pc(al.same_direction / al.stratum_n, 1)} to {pc(nw.same_direction / nw.stratum_n, 1)} and opposite-direction calls from {al.opposite} to {nw.opposite}; on the {n(na.stratum_n)} new variants agreement goes from {pc(na.same_direction / na.stratum_n, 1)} to {pc(nb.same_direction / nb.stratum_n, 1)} with {na.opposite} and {nb.opposite} opposite-direction calls, and the rate of opposite-direction calls is {xa.value:.1f} and {xb.value:.1f} times the April rate (Section 8).",
        f"The controls reproduce the main run under both classifiers (Section 6).",
        f"{len(numbers)} of the {n_checks} checks of `verify_numbers.py` give another value under v37 than the manuscript's ({sum(len(c[2]) for c in numbers)} of {n_values} values), {n_claims} claims differ in truth or wording, {len(outs['tables'])} of {outs['n_tables']} marked tables, {len(outs['sheets'][1])} of {outs['sheets'][0]} sheets of the workbook and {len(outs['figures'][1])} of {outs['figures'][0]} figures differ (Section 9).",
    ])
    D.s("The sections:")
    D.bullets(["1. what is compared and how", "2. the main run", "3. the opposite-direction calls", "4. the sweeps", "5. the ablations and the inputs", "6. the controls", "7. the criteria that fire more or less", "8. the October release", "9. where the manuscript's statements differ between v38 and v37"])


def write_md(pairs, main_pair, october, summary, eff, opp, counts) -> None:
    cmp_ = per_variant_comparison()
    check_typed_claims(summary, cmp_, eff)
    changes = numbers_that_change()
    claims = claims_of(pairs, main_pair, october, summary, cmp_)
    outs = derived_outputs_that_change()
    top, D = Doc(), Doc()
    section_summary(top, pairs, main_pair, october, summary, eff, opp, counts, changes, len(claims), outs)
    section_compared(D, pairs, main_pair, october)
    section_main(D, main_pair, cmp_)
    section_opposite(D, main_pair, october, pairs, opp, counts, cmp_)
    section_sweeps(D, summary)
    section_effects(D, pairs, main_pair, eff)
    section_controls(D, pairs, main_pair)
    section_criteria(D, pairs, main_pair, october)
    section_october(D, pairs, main_pair, october, opp)
    section_rewording(D, pairs, main_pair, october, summary, cmp_, changes, claims, outs)
    OUT_MD.write_text(top.text() + "\n" + D.text())
    print(f"wrote {OUT_MD.name} ({len(top.out) + len(D.out)} lines)")


def main() -> None:
    pairs, main_pair, october = load_pairs()
    print(f"{len(pairs)} pairs of runs")
    tidy_df = tidy([main_pair, october] + pairs)
    tidy_df.to_csv(OUT_CSV, index=False, lineterminator="\n")
    print(f"wrote {OUT_CSV.name} ({len(tidy_df)} rows)")
    sw, summary = sweep_table(pairs, main_pair)
    sw.to_csv(OUT_SWEEPS, index=False, lineterminator="\n")
    print(f"wrote {OUT_SWEEPS.name} ({len(sw)} rows)")
    eff = effects(pairs, main_pair)
    eff.to_csv(OUT_EFFECTS, index=False, lineterminator="\n")
    print(f"wrote {OUT_EFFECTS.name} ({len(eff)} rows)")
    pd.DataFrame(summary_rows(summary, eff)).to_csv(OUT_SUMMARY, index=False, lineterminator="\n")
    print(f"wrote {OUT_SUMMARY.name}")
    opp, counts = opposite_changes(pairs, main_pair, october)
    opp.to_csv(OUT_OPPOSITE, index=False, lineterminator="\n")
    print(f"wrote {OUT_OPPOSITE.name} ({len(opp)} rows)")
    write_md(pairs, main_pair, october, summary, eff, opp, counts)


if __name__ == "__main__":
    main()
