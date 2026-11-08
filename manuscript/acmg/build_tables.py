#!/usr/bin/env python3
"""Build the manuscript's data tables from `data/` and keep `manuscript.md` in step with them.

A table that carries numbers is written once, here, from the CSV files that
`analysis/acmg_benchmark/scripts/13_manuscript_stats.py` produces, and is placed in
`manuscript.md` between two marker comments:

    <!-- table:cutoffs:begin -->
    ...generated Markdown...
    <!-- table:cutoffs:end -->

`--write` replaces every marked block with the generated text; `--check` exits non-zero if
any block differs from what the data now says, which is how a stale number is caught. Tables
that are text (the guidance ledger, the evidence sources) are not generated and are not marked.

ClinVar's 2-star+ classifications are a reference and not ground truth.

Two bases. The tables are built from the files of a base, the classifier whose final run they describe: `--base v38` (the
default, the latest classifier, which is what the manuscript says) reads `runs_v38_variants.csv` and the other `_v38` files from
`analysis/acmg_benchmark/scripts/13_manuscript_stats.py`, and `--base v37` reads their twins for v37 (`runs_v37_variants.csv`, ...),
the classifier at the close of the review rounds, which is the before of what v38 changes. The code is spelled for v37 and every
name goes through `fname`, `rn` or `lab`, which turn it into the base's. The text of `manuscript.md` says v38, so under `--base v37`
`--check` says what would differ and `--write` is refused: `--write-to PATH` writes the manuscript with the tables of the base to a
file of your own, and the manuscript is not touched. The typed claims of Table 3 (a knee, a range of exchange rates, a flat sweep)
are computed from the steps of the base by the rules in `table_thresholds`, so that they say what the base's data say.

Usage
  build_tables.py --check     compare every marked block with its generated text
  build_tables.py --write     rewrite the marked blocks in place
  build_tables.py --print cutoffs    print one table
  build_tables.py --base v37 --print thresholds     the table for v37
  build_tables.py --base v37 --write-to /tmp/manuscript_v37_tables.md    the manuscript with v37's tables, written elsewhere
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
MANUSCRIPT = HERE / "manuscript.md"

N_ALL = 673_660  # every ClinVar 2-star+ record of the April benchmark


BASE = "v38"  # the classifier whose files the tables are built from (see the docstring)


def fname(name: str) -> str:
    """The file of a table. Names are spelled for v37, and `_v37` is the base's suffix for another base."""
    return name.replace("_v37", f"_{BASE}")


def rn(name: str) -> str:
    """A run's name in the family of the base: `v37_cut_lp4` is `v38_cut_lp4` under v38."""
    return BASE + name[3:] if name.startswith("v37") else name


def lab(text: str) -> str:
    """A label that names the classifier, which the code spells v37."""
    return text.replace("v37", BASE)


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(DATA / fname(name))


def n(x) -> str:
    """A count with thousands separators, rounded half up (Python's round() goes to even)."""
    return f"{int(float(x) + 0.5):,}"


def sci(v: float) -> str:
    """A threshold the way the text writes it: 4e-5, 2.5e-4, 1e-2; a value of 1 or more as is."""
    v = float(v)
    if v == 0:
        return "0"
    if v >= 1:
        return f"{v:g}"
    mantissa, exp = f"{v:e}".split("e")
    return f"{mantissa.rstrip('0').rstrip('.')}e{int(exp)}"


def pct(x, d: int = 1) -> str:
    return f"{float(x) * 100:.{d}f}%"


def md(header: list[str], rows: list[list[str]], align: list[str] | None = None) -> str:
    align = align or ["l"] + ["r"] * (len(header) - 1)
    rule = ["---:" if a == "r" else ":---" for a in align]
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join(rule) + "|"]
    lines += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(lines)


# ----------------------------------------------------------------------------- the tables


def table_cutoffs() -> str:
    """Table 4: the two point cutoffs of the pathogenic side, by real runs."""
    r = read("runs_v37_variants.csv")
    base = r[r.name == BASE].iloc[0]
    order = [
        ("Likely pathogenic from 4 points", rn("v37_cut_lp4")),
        ("Likely pathogenic from 5 points", rn("v37_cut_lp5")),
        ("Likely pathogenic from 6 points (default)", BASE),
        ("Likely pathogenic from 7 points", rn("v37_cut_lp7")),
        ("Likely pathogenic from 8 points", rn("v37_cut_lp8")),
        ("Likely pathogenic from 9 points", rn("v37_cut_lp9")),
        ("Pathogenic from 9 points", rn("v37_cut_p9")),
        ("Pathogenic from 10 points (default)", BASE),
        ("Pathogenic from 11 points", rn("v37_cut_p11")),
        ("Pathogenic from 12 points", rn("v37_cut_p12")),
    ]
    rows = []
    for label, name in order:
        x = r[r.name == name].iloc[0]
        rows.append([
            label,
            f"{n(x.same_pathogenic)} ({pct(x.plp_recall)})",
            n(x.vus_to_plp),
            n(x.false_pathogenic),
            n(x.false_benign),
            f"{n(x.same_direction)} ({pct(float(x.same_direction) / N_ALL)})",
            f"{n(x.exact)} ({pct(float(x.exact) / N_ALL)})",
        ])
    return md(
        ["Setting", "ClinVar P/LP called P/LP (recall)", "ClinVar VUS called P/LP", "ClinVar B/LB called P/LP",
         "ClinVar P/LP called LB/B", "Same-direction agreement", "Exact agreement"],
        rows,
    )



def table_inheritance() -> str:
    """Table 8: the classifier before and after reading ClinGen's curated inheritance, and the alternatives that were measured."""
    r = read("v38_runs.csv").set_index("run")
    order = [
        ("Before: inheritance read from the disease names", "v38_base"),
        ("**After (default)**: curated inheritance; BS2's carriers not counted beside BS1 or BA1", "v38"),
        ("Curated inheritance; BS2's carriers counted beside BS1 and BA1", "v38_alt_unguarded"),
        ("Default, and no BS2 carrier arm for a susceptibility or predisposition disorder", "v38_alt_label_guard"),
        ("Curated inheritance for PM2, PM3 and BP2; BS2 as before", "v38_alt_bs2_as_v37"),
        ("Default, with PP2, PM1 and PP3 capped together at 4 points", "v38_alt_triple_cap"),
    ]
    rows = []
    for label, name in order:
        x = r.loc[name]
        rows.append([
            label,
            f"{n(x.same_pathogenic)} ({pct(x.plp_recall)})",
            f"{n(x.same_benign)} ({pct(x.blb_recall)})",
            n(x.false_benign),
            n(x.false_pathogenic),
            f"{n(x.same_direction)} ({pct(x.same_direction / N_ALL)})",
            f"{n(x.exact)} ({pct(x.exact / N_ALL)})",
        ])
    return md(
        ["Setting", "ClinVar P/LP called P/LP (recall)", "ClinVar B/LB called LB/B (recall)", "Reference-discordant benign calls",
         "Reference-discordant pathogenic calls", "Same-direction agreement", "Exact agreement"],
        rows,
    )


def table_joint() -> str:
    """Table 9: gnomAD's joint exomes and genomes against the exomes, with the region flags separated out."""
    # The pilot on chromosomes 21 and 22 is in the file and in Appendix Table A16; the table is the runs of every chromosome.
    r = read("v38_joint.csv")
    r = r[r.set == "every chromosome"]
    rows = []
    for x in r.itertuples():
        reference = x.database == "gnomAD exomes"

        def change(total, delta):
            return f"{n(total)}" if reference else f"{n(total)} ({int(delta):+,})"

        rows.append([
            x.database,
            n(x.variants),
            n(x.same_pathogenic),
            n(x.same_benign),
            n(x.false_benign),
            n(x.false_pathogenic),
            change(x.same_direction, x.d_same_direction),
            change(x.exact, x.d_exact),
        ])
    return md(
        ["Source of the counts", "Variants", "ClinVar P/LP called P/LP", "ClinVar B/LB called LB/B", "Reference-discordant benign calls",
         "Reference-discordant pathogenic calls", "Same-direction agreement (change)", "Exact agreement (change)"],
        rows,
        align=["l"] + ["r"] * 7,
    )


def _step(steps: pd.DataFrame, knob: str, which: int) -> str:
    """One step of a sweep next to the default value: calls gained, and the error it cost."""
    sub = steps[(steps.knob == knob) & (steps.shipped_between)].reset_index(drop=True)
    if which >= len(sub):
        return "-"  # the sweep is ordered by what it gains, and the default is at the end of that order, with one step beside it
    x = sub.iloc[which]
    cost = int(x.added_errors)
    rate = "no added error" if cost == 0 else f"{n(x.gain_calls / cost)} per added error"
    return f"{sci(x.from_value)} to {sci(x.to_value)}: +{n(x.gain_calls)} calls, +{cost} ({rate})"


def table_thresholds() -> str:
    """Table 3: the thresholds that no guideline states, and what the base (v38 by default) shows."""
    r = read("runs_v37_variants.csv").set_index("name")
    steps = read("threshold_steps_v37.csv")
    vc = read("vcep_bars_v37.csv").set_index("bar")
    base = r.loc[BASE]
    pm2 = r[(r.knob == "pm2_ad_af_threshold") & (r.group == "sweep")]
    lo, hi = min(pm2.plp_recall.min(), base.plp_recall), max(pm2.plp_recall.max(), base.plp_recall)
    pm2r = r[(r.knob == "pm2_ar_af_threshold") & (r.group == "sweep")]
    rlo, rhi = min(pm2r.plp_recall.min(), base.plp_recall), max(pm2r.plp_recall.max(), base.plp_recall)
    # the default is one of the settings of the sweep, as it is for the recall beside it: leaving it out gave 26 calls under v38 where the text said 27
    rspan = int(max(pm2r.same_direction.max(), base.same_direction) - min(pm2r.same_direction.min(), base.same_direction))
    tb0, tb1 = r.loc[rn("v37_ix_table_pm2ad_0")], r.loc[rn("v37_ix_table_pm2ad_0.0001")]
    noan = r.loc[rn("v37_sw_no_an_floor")]
    off = r.loc[rn("v37_sw_phylop_1000")]
    bp1_1 = r.loc[rn("v37_sw_bp1veto_1")]
    bp1_loo1 = r.loc[rn("v37_hk_bp1_loo_veto1")]
    bp1_steps = steps[steps.knob == "BP1 veto count"]

    # The claims that the table made in typed words, computed from the steps of the base. The exchange rate of a step is the calls
    # it gains per error it adds, and a step that adds none has an infinite rate and is left out of every range ("no added error").
    def cpe(knob: str) -> pd.DataFrame:
        x = steps[steps.knob == knob]
        return x[x.calls_per_error != math.inf]

    def rate_from(knob: str, from_value: float) -> str:
        x = steps[(steps.knob == knob) & (abs(steps.from_value - from_value) < 1e-12)]
        if len(x) != 1:
            return "n/a (no step of the sweep starts there)"
        return "no added error" if x.iloc[0].calls_per_error == math.inf else n(x.iloc[0].calls_per_error)

    def joined(values: list[str]) -> str:
        return values[0] if len(values) == 1 else ", ".join(values[:-1]) + " and " + values[-1]

    bs1_below = rate_from("BS1 bar", 0.001)  # the step from the 1e-3 run to the one below it: where the sweep's knee is said to be
    hom_below = steps[(steps.knob == "BS2 homozygote bar") & (steps.from_value <= 0.001)].sort_values("from_value", ascending=False)
    carrier = r[(r.knob == "bs2_ad_prevalence_threshold") & (r.group == "sweep")]
    carrier_fb, carrier_fp = sorted(set(carrier.false_benign) | {base.false_benign}), sorted(set(carrier.false_pathogenic) | {base.false_pathogenic})
    far = cpe("BP7 far boundary")
    far_below, far_above = far[far.to_value <= 300].calls_per_error, far[far.from_value >= 300].calls_per_error
    knee = len(far_below) > 0 and len(far_above) > 0 and far_below.min() > far_above.max()
    phylop_min = cpe("BP7 PhyloP bar").calls_per_error.min()
    # The step v11 to v12 moved the BS1 threshold from 1e-2 to 5e-3; the run record says "3.1 pp of benign recall", which is
    # the Benign class alone, and the figure the text uses everywhere else is the benign side as a whole.
    tr = read("trajectory.csv").set_index("order")
    v12_missed = int(tr.loc[12].false_benign - tr.loc[11].false_benign)
    v12_pts = f"{(tr.loc[12].blb_recall - tr.loc[11].blb_recall) * 100:.1f}"
    v12_calls = (tr.loc[12].blb_recall - tr.loc[11].blb_recall) * 284_550
    rows = [
        ["PM2 threshold, genes not read as recessive-only: 4e-5",
         "Swept in v10 under the 2015 table; the median of the published VCEP thresholds is "
         + f"{sci(vc.loc['PM2', 'median'])} ({int(vc.loc['PM2', 'genes_with_a_bar'])} genes)",
         _step(steps, "PM2 bar, not recessive-only (point system)", 0), _step(steps, "PM2 bar, not recessive-only (point system)", 1),
         f"Recall {pct(lo)} to {pct(hi)} across strict absence to 5e-4 under the default scheme; "
         f"{pct(tb0.plp_recall)} to {pct(tb1.plp_recall)} under the 2015 table on both sides"],
        ["PM2 threshold, genes read as recessive-only: 7e-5",
         "Not from the SVI PM2 recommendation, which the code cited and which gives no frequency; equal to the Hearing Loss expert panel's recessive threshold. First swept here",
         _step(steps, "PM2 bar, recessive-only", 0), _step(steps, "PM2 bar, recessive-only", 1),
         f"Flat: recall {pct(rlo, 2)} to {pct(rhi, 2)} and agreement in direction within {rspan} calls from strict absence to 5e-4"],
        ["BS1 threshold: 5e-3",
         f"Swept on a 1-in-10 sample, where the move from 1e-2 looked free; on the full set it cost {v12_missed} more reference-discordant benign calls for {v12_pts} points of benign-side recall ({n(v12_calls)} calls), and was kept (v12). The median of the published VCEP thresholds is {sci(vc.loc['BS1', 'median'])}",
         _step(steps, "BS1 bar", 0), _step(steps, "BS1 bar", 1),
         f"Knee at 1e-3 ({bs1_below} per added error below it); 5e-3 is on the cautious side of the knee"],
        ["BS2 homozygote threshold, every gene not read as dominant-only: 1e-3",
         "Swept on the full benchmark (12 Aug): no knee, so set on the prevalence of the commonest recessive conditions",
         _step(steps, "BS2 homozygote bar", 0), _step(steps, "BS2 homozygote bar", 1),
         f"No knee: {joined([n(x) if x != math.inf else 'no added error' for x in hom_below.calls_per_error])} calls per added error below the default value"],
        ["BS2 carrier threshold, genes read as dominant-only: 3e-4, with an allele count of at least 5",
         lab("Added after review of one gene (SNRNP200) and swept on the full benchmark (v37)"), "-", "-",
         (f"Flat: {carrier_fb[0]} reference-discordant benign calls and {carrier_fp[0]} reference-discordant pathogenic calls at every setting from off to 3e-3"
          if len(carrier_fb) == 1 and len(carrier_fp) == 1 else
          f"Not flat: {carrier_fb[0]} to {carrier_fb[-1]} reference-discordant benign calls and {carrier_fp[0]} to {carrier_fp[-1]} reference-discordant pathogenic calls from off to 3e-3")],
        ["BP1 veto: 3 pathogenic missense variants in the gene",
         "Added in v9 with the veto, after the round-2 review of BP1 firing in genes with a pathogenic missense spectrum; the count was never swept before this paper",
         _step(steps, "BP1 veto count", 0), _step(steps, "BP1 veto count", 1),
         f"At most {n(bp1_steps.calls_per_error.replace(math.inf, 0).max())} calls per added error at every step, so there is no knee and the default count is on the aggressive side of the stricter ones; "
         f"a count of 1 gives +{n(bp1_1.same_direction - base.same_direction)} calls in direction and {n(base.same_benign - bp1_1.same_benign)} fewer correct benign-side calls, "
         f"and {int(base.false_benign - bp1_1.false_benign)} fewer reference-discordant benign calls ({int(base.false_benign - bp1_loo1.false_benign)} fewer when the veto leaves the variant's own record out)"],
        ["BP7 intronic distance: 300 nt",
         "Set at the knee of a sweep (v14)",
         _step(steps, "BP7 far boundary", 0), _step(steps, "BP7 far boundary", 1),
         (f"Knee confirmed: {n(far_below.min())} to {n(far_below.max())} calls per added error below 300, about {n(far_above.min())} to {n(far_above.max())} above"
          if knee else
          f"No clear knee: {n(far_below.min())} to {n(far_below.max())} calls per added error below 300 and {n(far_above.min())} to {n(far_above.max())} above")],
        ["BP7 PhyloP threshold: 2.0",
         "In the first commit, with no source recorded (the 2015 standards say \"not highly conserved\" and give no number). First swept here",
         _step(steps, "BP7 PhyloP bar", 0), _step(steps, "BP7 PhyloP bar", 1),
         f"No knee: at least {int(phylop_min // 100 * 100):,} calls per added error at every step; switching the test off adds {n(off.same_benign - base.same_benign)} "
         f"benign-side calls and {int(off.false_benign - base.false_benign)} reference-discordant benign calls (Section 3.3)"],
        ["Allele-number floor for BA1 and BS1: 2,000",
         "From the dataset-size condition of the BA1 recommendation, applied here to the site's allele number and to BS1", "-", "-",
         f"Removing it adds {n(noan.same_direction - base.same_direction)} calls in direction and no opposite-direction call"],
        ["Likely pathogenic from 6 points",
         "Fixed by the point system (Tavtigian et al., 2020); varied here as a sensitivity analysis", "-", "-",
         f"5 points: +{n(r.loc[rn('v37_cut_lp5')].same_pathogenic - base.same_pathogenic)} pathogenic-side calls, "
         f"+{n(r.loc[rn('v37_cut_lp5')].vus_to_plp - base.vus_to_plp)} ClinVar-VUS and "
         f"+{n(r.loc[rn('v37_cut_lp5')].false_pathogenic - base.false_pathogenic)} benign-side called pathogenic-side (Table 4)"],
        ["Pathogenic from 10 points",
         "Fixed by the point system", "-", "-",
         f"9 to 12 points move agreement in direction by at most {n(max(abs(r.loc[rn(f'v37_cut_p{c}')].same_direction - base.same_direction) for c in (9, 11, 12)))} calls"],
        ["BA1 at 5%",
         f"Fixed by the guideline; the median of the published VCEP thresholds is {sci(vc.loc['BA1', 'median'])}", "-", "-",
         "Not varied"],
    ]
    return md(
        ["Threshold: default value", "How it was set", "Last step before the default value", "First step beyond it", lab("On v37")],
        rows, align=["l", "l", "l", "l", "l"],
    )


def table_decomposition() -> str:
    """Table 5: every ClinVar 2-star+ record, by how the base's call relates to ClinVar's class."""
    c = read("confusion_5x5_v37.csv").pivot(index="clinvar", columns="fastvep", values="n")
    tot = int(c.to_numpy().sum())
    if tot != N_ALL:
        raise SystemExit(f"the 5x5 matrix holds {tot} records, not {N_ALL}")
    P, LP, V, LB, B = "Pathogenic", "Likely_pathogenic", "VUS", "Likely_benign", "Benign"
    exact = int(sum(c.loc[t, t] for t in (P, LP, V, LB, B)))
    tier = int(c.loc[P, LP] + c.loc[LP, P] + c.loc[LB, B] + c.loc[B, LB])
    p_vus = int(c.loc[P, V] + c.loc[LP, V])
    b_vus = int(c.loc[LB, V] + c.loc[B, V])
    v_p = int(c.loc[V, P] + c.loc[V, LP])
    v_b = int(c.loc[V, LB] + c.loc[V, B])
    fb = int(sum(c.loc[t, x] for t in (P, LP) for x in (LB, B)))
    fp = int(sum(c.loc[t, x] for t in (LB, B) for x in (P, LP)))
    nocall = int(c["NoCall"].sum())
    parts = [exact, tier, p_vus, b_vus, v_p, v_b, fb + fp, nocall]
    if sum(parts) != N_ALL:
        raise SystemExit(f"the kinds sum to {sum(parts)}, not {N_ALL}")
    rows = [
        ["Same class as ClinVar", "the diagonal", n(exact), pct(exact / N_ALL)],
        ["Same direction, different tier", "P and LP, LB and B", n(tier), pct(tier / N_ALL)],
        ["&nbsp;&nbsp;ClinVar Pathogenic called Likely pathogenic", "", n(c.loc[P, LP]), pct(c.loc[P, LP] / N_ALL)],
        ["&nbsp;&nbsp;ClinVar Benign called Likely benign", "", n(c.loc[B, LB]), pct(c.loc[B, LB] / N_ALL)],
        ["ClinVar P/LP called VUS", "a pathogenic-side class, no call made", n(p_vus), pct(p_vus / N_ALL)],
        ["ClinVar LB/B called VUS", "a benign-side class, no call made", n(b_vus), pct(b_vus / N_ALL)],
        ["ClinVar VUS called P/LP", "", n(v_p), pct(v_p / N_ALL)],
        ["ClinVar VUS called LB/B", "", n(v_b), pct(v_b / N_ALL)],
        ["Opposite direction", "", n(fb + fp), pct((fb + fp) / N_ALL, 3)],
        ["&nbsp;&nbsp;ClinVar P/LP called LB/B (reference-discordant benign calls)", "", n(fb), pct(fb / N_ALL, 3)],
        ["&nbsp;&nbsp;ClinVar LB/B called P/LP", "", n(fp), pct(fp / N_ALL, 3)],
        ["No call", "records with no alternative allele (ALT is a dot)", n(nocall), pct(nocall / N_ALL, 2)],
        ["All records", "", n(N_ALL), "100%"],
    ]
    return md(["Kind", "Cells", "Variants", "Share of all records"], rows, align=["l", "l", "r", "r"])


def table_residual() -> str:
    """Table 6: the opposite-direction calls of the base, by the kind of evidence behind the call and the review ruling."""
    g = read("opposite_direction_by_group_v37.csv").set_index("group")
    order = [
        "PVS1 on a ClinVar benign-side variant",
        "BP4+BP7 on a ClinVar pathogenic-side variant",
        "PM1 with PP2 and PP3 on a ClinVar benign-side variant",
        "BS1 on a ClinVar pathogenic-side variant",
        "BP1 on a ClinVar pathogenic-side variant",
    ]
    label = {
        "PVS1 on a ClinVar benign-side variant": "PVS1 on a ClinVar benign or likely benign variant",
        "BP4+BP7 on a ClinVar pathogenic-side variant": "BP4 and BP7 on a ClinVar pathogenic or likely pathogenic variant",
        "PM1 with PP2 and PP3 on a ClinVar benign-side variant": "PM1, PP2 and PP3 on a ClinVar benign or likely benign variant",
        "BS1 on a ClinVar pathogenic-side variant": "BS1 on a ClinVar pathogenic or likely pathogenic variant",
        "BP1 on a ClinVar pathogenic-side variant": "BP1 on a ClinVar pathogenic or likely pathogenic variant",
    }
    cols = ["clinvar", "fastvep", "neither", "needs_evidence", "disorder_specific"]
    heads = ["ClinVar supported", "fastVEP supported", "Neither", "Needs more evidence", "Depends on the disorder"]
    if "not reviewed" in g.columns:  # a later classifier can make a call opposite that the reviewer never saw
        cols, heads = cols + ["not reviewed"], heads + ["Not reviewed"]
    rows = []
    for k in order:
        x = g.loc[k]
        rows.append([label[k], n(x["total"])] + [n(x[c]) for c in cols])
    tot = g.loc[order]
    rows.append(["All", n(tot["total"].sum())] + [n(tot[c].sum()) for c in cols])
    return md(["Evidence behind the call", "Calls"] + heads, rows)


def table_temporal() -> str:
    """Table 7: the base on the April set and on each stratum of the later release."""
    t = read("temporal_v37.csv")
    bc = read("temporal_by_class_v37.csv")
    c3 = read("confusion_3x3_v37.csv")

    def block(label, stratum):
        sub = t[t.stratum == stratum]
        first = sub.iloc[0]
        def rec(a, b):
            x = sub[(sub.clinvar == a) & (sub.fastvep == b)].iloc[0]
            return pct(x.row_fraction) if x.row_total else "-"
        return [label, n(first.stratum_n), pct(first.same_direction / first.stratum_n), rec("P/LP", "P/LP"), rec("B/LB", "B/LB"), rec("VUS", "VUS"),
                f"{int(first.opposite)} ({int(first.false_benign)} and {int(first.false_pathogenic)})"]

    def april_row():
        def cell(a, b):
            return int(c3[(c3.clinvar == a) & (c3.fastvep == b)].iloc[0].n)

        def rec(a, b):
            return pct(c3[(c3.clinvar == a) & (c3.fastvep == b)].iloc[0].row_fraction)

        same = cell("P/LP", "P/LP") + cell("VUS", "VUS") + cell("B/LB", "B/LB")
        fb, fp = cell("P/LP", "B/LB"), cell("B/LB", "P/LP")
        return ["April release, the benchmark set", n(N_ALL), pct(same / N_ALL), rec("P/LP", "P/LP"), rec("B/LB", "B/LB"), rec("VUS", "VUS"),
                f"{fb + fp} ({fb} and {fp})"]

    rows = [
        april_row(),
        block("October release, all 2-star+ variants", "all"),
        block("&nbsp;&nbsp;unchanged since April", "unchanged"),
        block("&nbsp;&nbsp;class moved within its direction", "tier_change"),
        block("&nbsp;&nbsp;class moved across a direction", "direction_change"),
        block("&nbsp;&nbsp;new: not in the April set", "new"),
    ]
    return md(["Set", "Variants", "Same-direction agreement", "ClinVar P/LP called P/LP", "ClinVar LB/B called LB/B", "ClinVar VUS called VUS",
               "Opposite-direction calls (reference-discordant benign and reference-discordant pathogenic)"], rows)

def table_databases() -> str:
    """Table 10: v37 and v38 with the first databases and with the complete ones, from `databases_old_vs_new.csv`.

    The file is made by `21_compare_databases.py` from the per-variant tables of the runs on both stacks and holds both versions
    whatever the base is, so the table is the same under `--base v37`. The whole October release was scored only with the complete
    databases (the first drafts scored the 772,585 variants that a gnomAD extract covered), so its first-database cells are empty.
    """
    d = pd.read_csv(DATA / "databases_old_vs_new.csv")

    def get(version: str, scope: str, measure: str, db: str) -> float:
        x = d[(d.version == version) & (d.set == scope) & (d.measure == measure)]
        return float(x.iloc[0][db])

    def row(label: str, scope: str, measure: str, fmt) -> list[str]:
        cells = []
        for version in ("v37", "v38"):
            for db in ("first", "complete"):
                cells.append("-" if math.isnan(get(version, scope, measure, db)) else fmt(version, scope, measure, db))
        return [label] + cells

    def count_rate(version: str, scope: str, measure: str, db: str) -> str:
        return f"{n(get(version, scope, measure, db))} ({pct(get(version, scope, measure, db) / get(version, scope, 'variants', db), 1)})"

    def opposite(version: str, scope: str, measure: str, db: str) -> str:
        return f"{n(get(version, scope, 'opposite', db))} ({n(get(version, scope, 'false_benign', db))} and {n(get(version, scope, 'false_pathogenic', db))})"

    def plain(version: str, scope: str, measure: str, db: str) -> str:
        return n(get(version, scope, measure, db))

    def share(version: str, scope: str, measure: str, db: str) -> str:
        return pct(get(version, scope, measure, db), 1)

    april, covered, new, whole = "April set", "October release, the variants of the first design", "October release, new variants of the first design", "October release, all variants"
    k = lambda scope: n(get("v38", scope, "variants", "complete"))  # noqa: E731
    rows = [
        row(f"**April set** ({k(april)} variants): exact agreement", april, "exact", count_rate),
        row("Agreement in direction", april, "same_direction", count_rate),
        row("Opposite-direction calls (reference-discordant benign and reference-discordant pathogenic)", april, "opposite", opposite),
        row("Recall of ClinVar pathogenic-side variants", april, "plp_recall", share),
        row("Recall of ClinVar benign-side variants", april, "blb_recall", share),
        row("Variants with BP4", april, "BP4", plain),
        row("Variants with BP7", april, "BP7", plain),
        row("Variants with PP3", april, "PP3", plain),
        row(f"Canonical splice variants with PVS1 (of {n(get('v38', april, 'canonical_splice_variants', 'complete'))})", april, "PVS1_on_canonical_splice", plain),
        row(f"**October release**, the {k(covered)} variants that the first drafts scored: agreement in direction", covered, "same_direction", count_rate),
        row(f"Opposite-direction calls, of the {k(covered)}", covered, "opposite", opposite),
        row(f"The {k(new)} new variants among them: agreement in direction", new, "same_direction", count_rate),
        row(f"Opposite-direction calls, of the {k(new)}", new, "opposite", opposite),
        row(f"The whole release ({k(whole)} variants): agreement in direction", whole, "same_direction", count_rate),
        row(f"Opposite-direction calls, of the {k(whole)}", whole, "opposite", opposite),
    ]
    return md(["Measure", "v37, first databases", "v37, complete databases", "v38, first databases", "v38, complete databases"], rows)


TABLES = {
    "databases": table_databases,
    "inheritance": table_inheritance,
    "joint": table_joint,
    "thresholds": table_thresholds,
    "cutoffs": table_cutoffs,
    "decomposition": table_decomposition,
    "residual": table_residual,
    "temporal": table_temporal,
}


# ---------------------------------------------------------------------------- the driver

BLOCK = re.compile(r"(<!-- table:(?P<name>[a-z0-9_]+):begin -->\n)(?P<body>.*?)(\n<!-- table:(?P=name):end -->)", re.S)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--print", dest="show", default="")
    ap.add_argument("--base", choices=["v37", "v38"], default=os.environ.get("ACMG_BASE", "v38"),
                    help="the classifier whose files the tables are built from (default v38, or ACMG_BASE like the other scripts; see the docstring)")
    ap.add_argument("--write-to", dest="write_to", default="", help="write the manuscript with the tables of the base to this file and leave manuscript.md as it is")
    args = ap.parse_args()
    global BASE
    BASE = args.base
    if BASE not in ("v37", "v38"):
        sys.exit(f"unknown base {BASE}; choose v37 or v38")
    if args.write and BASE != "v38":
        sys.exit(f"manuscript.md says v38, and --write would put the tables of {BASE} in its text: use --write-to PATH")
    if args.show:
        print(TABLES[args.show]())
        return
    text = MANUSCRIPT.read_text()
    stale, seen = [], set()

    def sub(m):
        name = m.group("name")
        if name not in TABLES:
            sys.exit(f"manuscript.md marks a table called {name!r} that this script does not build")
        seen.add(name)
        new = TABLES[name]()
        if m.group("body") != new:
            stale.append(name)
        return m.group(1) + new + m.group(4)

    out = BLOCK.sub(sub, text)
    missing = sorted(set(TABLES) - seen)
    if args.write_to:
        Path(args.write_to).write_text(out)
        print(f"wrote {args.write_to}: {len(seen)} tables ({len(stale)} differ from the text of manuscript.md)")
        return
    if args.write:
        MANUSCRIPT.write_text(out)
        print(f"wrote {len(seen)} tables ({len(stale)} changed)")
    else:
        for name in stale:
            print(f"stale: {name}")
        for name in missing:
            print(f"not in the manuscript: {name}")
        print(f"{len(seen)} marked tables, {len(stale)} stale")
        if stale:
            sys.exit(1)


if __name__ == "__main__":
    main()
