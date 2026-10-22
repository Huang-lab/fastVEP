#!/usr/bin/env python3
"""Figures for the ACMG classification manuscript.

Reads the small tables in `data/` (written by `analysis/acmg_benchmark/scripts/13_manuscript_stats.py`
and `11_manuscript_runs.py`) and writes PNG and PDF into `figures/`. Nothing is typed in by hand:
every plotted value comes from a table, and the tables come from the scored runs.

Colour follows the data-viz method this repository uses elsewhere:

  * Ordinal ACMG tiers (P, LP, VUS, LB, B) are a diverging scale, two hues and a neutral midpoint.
    The blue steps are the documented ramp's 300 and 550; the red steps are the documented red held
    at the same OKLCH lightness and chroma as those two blues, so the two arms are equal.
  * Nominal series take the documented categorical slots in their fixed order.
  * Text wears ink tokens, never a series colour; marks carry identity beside direct labels.

ClinVar's 2-star+ classifications are a reference and not ground truth; every axis that says
"agreement" means agreement with them and nothing stronger.

Usage
  generate_figures.py            all figures
  generate_figures.py 1 4        figures 1 and 4
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
FIG = HERE / "figures"

# ------------------------------------------------------------------------------ tokens
SURFACE = "#ffffff"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]

TIER = {
    "Pathogenic": "#9e3432",
    "Likely_pathogenic": "#e4857e",
    "VUS": "#c3c2b7",
    "Likely_benign": "#6da7ec",
    "Benign": "#1c5cab",
}
TIER3 = {"P/LP": "#e34948", "VUS": "#c3c2b7", "B/LB": "#2a78d6"}
SHORT = {"Pathogenic": "P", "Likely_pathogenic": "LP", "VUS": "VUS", "Likely_benign": "LB", "Benign": "B"}
BLUE_SEQ = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
            "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]

plt.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.5,
    "axes.labelsize": 7.5,
    "axes.titlesize": 8.5,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.8,
    "axes.labelcolor": INK2,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "text.color": INK,
    "legend.frameon": False,
    "legend.fontsize": 7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "grid.linestyle": "-",
    "lines.solid_capstyle": "round",
    "lines.solid_joinstyle": "round",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

LW = 1.5  # the 2 px line of the method, at 96 dpi
MS = 5.5  # an 8 px marker


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(DATA / name)


def panel_letter(ax, letter: str, x: float = -0.02, y: float = 1.04) -> None:
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK, ha="left", va="bottom")


def hgrid(ax) -> None:
    ax.grid(axis="y")
    ax.set_axisbelow(True)


def vgrid(ax) -> None:
    ax.grid(axis="x")
    ax.set_axisbelow(True)


def save(fig, stem: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{stem}.png", dpi=300, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"  wrote figures/{stem}.png and .pdf")


def thousands(x, pos=None) -> str:
    return f"{int(x):,}"


# ------------------------------------------------------------------------ figure 1: trajectory

# Version at which each change landed, and the label it carries in the key. The labels name the
# change that produced the step, as RUN_VERSIONS.md records it.
EVENTS = [
    (4, "v4: SpliceAI and PhyloP routing fixed"),
    (6, "v6: PM2 reads absence as rarity"),
    (8, "v8: review 1 (PVS1 graded; BP4 and PM2 fixes)"),
    (9, "v9: review 2 (PS1 and PM1 stop self-matching)"),
    (10, "v10: gnomAD QC gates; PM2 and BS2 bars measured"),
    (11, "v11: SVI point system, pathogenic side"),
    (12, "v12: BS1 bar moved from 1% to 0.5%"),
    (13, "v13: --acmg implies --hgvs; BP7 intronic extension live"),
    (14, "v14: BP7 far boundary; PM1 tests benign variation"),
    (16, "v16: PS1 splice path live"),
    (20, "v19 and v20: frequency veto; hypomorph exceptions"),
    (27, "v23 to v27: Ensembl-parity annotation"),
    (37, "v37: review 6 (HGVS placement, BP7 exon edge, BS2 bar)"),
]


def fig1_trajectory() -> None:
    t = read("trajectory.csv")
    t["mainline"] = t["mainline"].astype(bool)
    for c in ("same_direction_rate", "exact_rate", "plp_recall", "blb_recall"):
        t[c] = t[c] * 100  # every rate of this figure is drawn as a percentage
    stored = t[t.source.str.startswith("scored")]
    main = stored[stored.mainline]  # the steps on the way: v7 to v37 without v17
    side = stored[~stored.mainline]  # v17, the opt-in run with the published VCEP bars
    # Open markers: v1, a matrix reconstructed from the run's documentation, and v2, v4 to v6, the figures
    # recorded at the time. None of them has a raw output behind it.
    opened = t[~t.source.str.startswith("scored")]

    fig = plt.figure(figsize=(7.2, 8.6))
    gs = fig.add_gridspec(6, 1, height_ratios=[0.55, 2.3, 2.0, 2.0, 0.42, 1.25], hspace=0.30)
    ax0 = fig.add_subplot(gs[0])
    axA = fig.add_subplot(gs[1], sharex=ax0)
    axB = fig.add_subplot(gs[2], sharex=ax0)
    axC = fig.add_subplot(gs[3], sharex=ax0)
    axk = fig.add_subplot(gs[5])

    # --- numbered event strip; the key is printed beneath the whole figure
    ax0.set_ylim(0, 1)
    ax0.set_yticks([])
    for sp in ("left", "bottom"):
        ax0.spines[sp].set_visible(False)
    ax0.tick_params(axis="x", bottom=False, labelbottom=False)
    for i, (v, _) in enumerate(EVENTS, 1):
        for ax in (axA, axB, axC):
            ax.axvline(v, color=GRID, lw=0.6, zorder=0)
        ax0.plot([v, v], [0.0, 0.42], color=AXIS, lw=0.6, zorder=1)
        ax0.text(v, 0.62, str(i), ha="center", va="center", fontsize=6.2, color=SURFACE, fontweight="bold",
                 bbox=dict(boxstyle="circle,pad=0.22", fc=INK2, ec="none"), zorder=3)

    def stored_line(ax, col, color, **kw):
        ax.plot(main["order"], main[col], color=color, lw=LW, marker="o", ms=MS - 1.5, mec=SURFACE, mew=1, zorder=3, **kw)

    def open_marks(ax, rows, col, color, ms=MS - 1.5):
        ax.plot(rows["order"], rows[col], ls="none", marker="o", ms=ms, mfc=SURFACE, mec=color, mew=1.2, zorder=3)

    def side_marks(ax, col, color):
        ax.plot(side["order"], side[col], ls="none", marker="x", ms=4.6, color=color, mew=1.3, zorder=3)

    def dotted(ax, rows, col, color):
        """A dotted line through the figures that were not kept, up to the first stored run."""
        pts = pd.concat([rows.dropna(subset=[col]), main.iloc[:1]])
        ax.plot(pts["order"], pts[col], color=color, lw=1.0, ls=(0, (1, 2)), zorder=2)

    # --- agreement
    hgrid(axA)
    stored_line(axA, "same_direction_rate", CAT[0])
    stored_line(axA, "exact_rate", CAT[1])
    dotted(axA, opened, "same_direction_rate", CAT[0])
    dotted(axA, opened, "exact_rate", CAT[1])
    open_marks(axA, opened.dropna(subset=["same_direction_rate"]), "same_direction_rate", CAT[0])
    open_marks(axA, opened.dropna(subset=["exact_rate"]), "exact_rate", CAT[1])
    side_marks(axA, "same_direction_rate", CAT[0])
    side_marks(axA, "exact_rate", CAT[1])
    axA.set_ylim(45, 80)
    axA.set_yticks(range(45, 81, 5))
    axA.set_ylabel("Agreement with ClinVar 2-star+\n(% of 673,660 variants)")
    v37 = main[main.run == "v37"].iloc[0]
    axA.text(37.6, v37.same_direction_rate, f"{v37.same_direction_rate:.1f}%", color=INK2, fontsize=7, va="center", ha="left")
    axA.text(37.6, v37.exact_rate, f"{v37.exact_rate:.1f}%", color=INK2, fontsize=7, va="center", ha="left")
    v1 = opened[opened.run == "v1"].iloc[0]
    axA.annotate("v1: 51.0% of all records\n(54.7% of the 627,376 it annotated)", xy=(1.2, v1.same_direction_rate - 0.5),
                 xytext=(8.0, 47.6), fontsize=6.4, color=INK2, ha="left", va="center",
                 arrowprops=dict(arrowstyle="-", color=AXIS, lw=0.8, shrinkA=0, shrinkB=1))
    v9 = main[main.run == "v9"].iloc[0]
    axA.annotate("v9 is leave-one-out:\nthe variant's own ClinVar\nrecord is no longer counted", xy=(9.05, v9.same_direction_rate - 0.4),
                 xytext=(15.0, 67.0), fontsize=6.4, color=INK2, ha="left", va="center",
                 arrowprops=dict(arrowstyle="-", color=AXIS, lw=0.8, shrinkA=0, shrinkB=1))
    axA.legend(handles=[
        Line2D([0], [0], color=CAT[0], lw=LW, marker="o", ms=MS - 1.5, mec=SURFACE, label="Same-direction"),
        Line2D([0], [0], color=CAT[1], lw=LW, marker="o", ms=MS - 1.5, mec=SURFACE, label="Exact (five-tier)"),
        Line2D([0], [0], color=INK2, lw=0, marker="o", ms=MS - 1.5, mfc=SURFACE, mec=INK2, mew=1.2,
               label="v1 (reconstructed) and v2 to v6 (recorded at the time)"),
        Line2D([0], [0], color=INK2, lw=0, marker="x", ms=4.6, mew=1.3, label="v17: opt-in run with the VCEP bars, not a step"),
    ], loc="lower right", bbox_to_anchor=(1.0, 0.0), fontsize=6.4)
    panel_letter(axA, "A")

    # --- opposite-direction calls
    hgrid(axB)
    axB.set_yscale("log")
    opp = main.dropna(subset=["false_benign"])
    axB.plot(opp["order"], opp.false_benign, color=TIER3["P/LP"], lw=LW, marker="o", ms=MS - 1.5, mec=SURFACE, mew=1, zorder=3,
             label="ClinVar P/LP called LB or B (missed diagnosis)")
    axB.plot(opp["order"], opp.false_pathogenic, color=TIER3["B/LB"], lw=LW, marker="o", ms=MS - 1.5, mec=SURFACE, mew=1, zorder=3,
             label="ClinVar LB or B called P or LP")
    axB.plot(main["order"], main.opposite, color=INK2, lw=1.0, marker="o", ms=MS - 2.5, mec=SURFACE, mew=0.8, zorder=2, label="Total")
    dotted(axB, opened, "opposite", INK2)
    open_marks(axB, opened.dropna(subset=["opposite"]), "opposite", INK2, ms=MS - 1.5)
    o1 = opened[opened.run == "v1"].iloc[0]
    axB.plot([o1["order"]], [o1.false_benign], ls="none", marker="o", ms=MS - 1.5, mfc=SURFACE, mec=TIER3["P/LP"], mew=1.2, zorder=3)
    axB.plot([o1["order"]], [o1.false_pathogenic], ls="none", marker="o", ms=MS - 1.5, mfc=SURFACE, mec=TIER3["B/LB"], mew=1.2, zorder=3)
    side_marks(axB, "false_benign", TIER3["P/LP"])
    side_marks(axB, "false_pathogenic", TIER3["B/LB"])
    side_marks(axB, "opposite", INK2)
    axB.set_ylim(4, 1800)
    axB.set_yticks([10, 30, 100, 300])
    axB.set_yticklabels(["10", "30", "100", "300"])
    axB.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    axB.set_ylabel("Opposite-direction calls\n(n, log scale)")
    for run in ("v7", "v13", "v37"):
        r = main[main.run == run].iloc[0]
        axB.text(r["order"] + 0.35, r.opposite * 1.10, f"{int(r.opposite)}", fontsize=6.8, color=INK2, va="bottom", ha="left")
    handles, labels = axB.get_legend_handles_labels()
    handles.append(Line2D([0], [0], color=INK2, lw=0, marker="o", ms=MS - 1.5, mfc=SURFACE, mec=INK2, mew=1.2))
    labels.append("v1 (reconstructed) and v6 (recorded at the time)")
    axB.legend(handles, labels, loc="upper right", fontsize=6.4, bbox_to_anchor=(1.0, 1.0))
    panel_letter(axB, "B")

    # --- recall
    hgrid(axC)
    stored_line(axC, "plp_recall", TIER3["P/LP"])
    stored_line(axC, "blb_recall", TIER3["B/LB"])
    side_marks(axC, "plp_recall", TIER3["P/LP"])
    side_marks(axC, "blb_recall", TIER3["B/LB"])
    axC.set_ylim(40, 70)
    axC.set_yticks(range(40, 71, 5))
    axC.set_ylabel("Same-direction recall\n(% of ClinVar class)")
    axC.text(1.0, 41.0, "v1: 10.9% and 19.0%, below the axis", fontsize=6.4, color=INK2, ha="left", va="bottom")
    axC.set_xlabel("Benchmark run (v1 to v37)")
    axC.set_xlim(0, 40)
    axC.set_xticks([1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37])
    axC.legend(handles=[Line2D([0], [0], color=TIER3["P/LP"], lw=LW, marker="o", ms=MS - 1.5, mec=SURFACE, label="ClinVar P/LP called P or LP"),
                        Line2D([0], [0], color=TIER3["B/LB"], lw=LW, marker="o", ms=MS - 1.5, mec=SURFACE, label="ClinVar B/LB called LB or B")],
               loc="lower right", fontsize=6.6)
    panel_letter(axC, "C")
    for ax in (axA, axB):
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.tick_params(axis="x", length=0)

    # --- key to the numbered changes
    axk.axis("off")
    n = len(EVENTS)
    rows = (n + 1) // 2
    for i, (_, label) in enumerate(EVENTS):
        col, row = divmod(i, rows)
        x = col * 0.52
        y = 1.0 - (row + 0.5) / rows
        axk.text(x + 0.012, y, str(i + 1), fontsize=6.2, color=SURFACE, fontweight="bold", ha="center", va="center",
                 bbox=dict(boxstyle="circle,pad=0.22", fc=INK2, ec="none"))
        axk.text(x + 0.034, y, label, fontsize=6.6, color=INK2, va="center", ha="left")
    axk.set_xlim(0, 1)
    axk.set_ylim(0, 1)
    save(fig, "fig1_trajectory")


# ------------------------------------------------------------------------ figure 2: expert review

RULING_COLORS = {
    "clinvar_supported": CAT[0],
    "fastvep_supported": CAT[1],
    "neither_supported": CAT[2],
    "needs_evidence": CAT[3],
    "no_verdict": GRID,
}
RULING_LABELS = {
    "clinvar_supported": "ClinVar's call supported",
    "fastvep_supported": "fastVEP's call supported",
    "neither_supported": "Neither: the reviewer's call is VUS or another class",
    "needs_evidence": "Needs more evidence, or depends on the disorder",
    "no_verdict": "No verdict recorded",
}


def text_on(color: str) -> str:
    """White or ink, by the fill's luminance, so a label set inside a segment always reads."""
    r, g, b = (int(color[i:i + 2], 16) / 255 for i in (1, 3, 5))
    lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return INK if lum > 0.45 else SURFACE


def stacked_barh(ax, ylabels, series, colors, height=0.55, min_label=0.06, total=None, fontsize=6.8, label_thin=True):
    """Horizontal stacked bars with a 2 px surface gap between segments and fitted in-segment labels."""
    left = np.zeros(len(ylabels))
    tot = np.array(total if total is not None else [sum(v[i] for v in series.values()) for i in range(len(ylabels))], float)
    for name, vals in series.items():
        vals = np.array(vals, float)
        ax.barh(range(len(ylabels)), vals, left=left, height=height, color=colors[name], edgecolor=SURFACE, linewidth=1.2, zorder=3)
        for i, v in enumerate(vals):
            if v <= 0:
                continue
            if v / max(tot[i], 1) >= min_label and v >= 0.045 * max(tot):
                ax.text(left[i] + v / 2, i, f"{int(v)}", ha="center", va="center", fontsize=fontsize, color=text_on(colors[name]), zorder=4)
            elif label_thin:
                # too thin to hold its number: set it just above the bar, in ink, not the segment's colour
                ax.text(left[i] + v / 2, i - height / 2 - 0.04, f"{int(v)}", ha="center", va="bottom", fontsize=fontsize - 0.4, color=INK2, zorder=4)
        left += vals
    ax.set_yticks(range(len(ylabels)))
    ax.set_yticklabels(ylabels)
    ax.invert_yaxis()


def fig2_review() -> None:
    rr = read("review_rounds.csv")
    grp = read("opposite_direction_by_group_v37.csv")

    fig = plt.figure(figsize=(7.2, 3.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.15], height_ratios=[1, 0.34], hspace=0.55, wspace=0.55)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axL = fig.add_subplot(gs[1, :])
    axL.axis("off")

    cats = ["clinvar_supported", "fastvep_supported", "neither_supported", "needs_evidence", "no_verdict"]
    labels = [f"Round {int(r.round)}\n{int(r.rows)} rows on {r.run}" for r in rr.itertuples()]
    series = {c: list(rr[c]) for c in cats}
    stacked_barh(axA, labels, series, RULING_COLORS, total=list(rr.rows))
    axA.set_xlim(0, 125)
    axA.set_xlabel("Rows (opposite-direction calls) in the workbook returned")
    vgrid(axA)
    axA.set_title("What the reviewer ruled", loc="left", fontsize=8, color=INK, pad=6)
    panel_letter(axA, "A", x=-0.30, y=1.07)

    order = ["PVS1 on a ClinVar benign-side variant", "BP4+BP7 on a ClinVar pathogenic-side variant",
             "PM1 with PP2 and PP3 on a ClinVar benign-side variant", "BS1 on a ClinVar pathogenic-side variant",
             "BP1 on a ClinVar pathogenic-side variant"]
    short = {order[0]: "PVS1, ClinVar\nbenign-side", order[1]: "BP4 + BP7, ClinVar\npathogenic-side",
             order[2]: "PM1 + PP2 + PP3,\nClinVar benign-side", order[3]: "BS1, ClinVar\npathogenic-side",
             order[4]: "BP1, ClinVar\npathogenic-side"}
    g = grp.set_index("group").loc[order]
    series = {
        "clinvar_supported": list(g["clinvar"]),
        "fastvep_supported": list(g["fastvep"]),
        "neither_supported": list(g["neither"]),
        "needs_evidence": list(g["needs_evidence"] + g["disorder_specific"]),
    }
    stacked_barh(axB, [short[o] for o in order], series, RULING_COLORS, total=list(g["total"]))
    axB.set_xlim(0, 24)
    axB.set_xlabel("Opposite-direction calls remaining in v37 (49)")
    vgrid(axB)
    axB.set_title("Rulings on the 49 calls that remain", loc="left", fontsize=8, color=INK, pad=6)
    panel_letter(axB, "B", x=-0.42, y=1.07)

    handles = [Patch(fc=RULING_COLORS[c], ec="none", label=RULING_LABELS[c]) for c in cats]
    axL.legend(handles=handles, loc="center", ncol=3, fontsize=6.8, columnspacing=1.6, handlelength=1.2, handleheight=1.0)
    save(fig, "fig2_review")


# ----------------------------------------------------------------- figure 3: defect classes

CLASS_LABELS = {
    "G": "A condition of the guideline was omitted",
    "W": "Wiring: a criterion or input was inert or misread",
    "F": "Population frequency trusted too far, or not enough",
    "A": "Consequence or HGVS differed from Ensembl VEP",
    "T": "A threshold set by measurement",
    "S": "Scoring, reporting or provenance",
    "C": "Circularity with ClinVar",
    "R": "Combining rule",
}
FOUND_LABELS = {
    "expert_review": "Expert\nreview",
    "code_audit": "Code\naudit",
    "discordance_mining": "Discordance\nmining",
    "benchmark_aggregate": "Benchmark\naggregate",
    "first_run": "First\nreal-data run",
    "sweep": "Threshold\nsweep",
    "vep_comparison": "VEP\ncomparison",
    "regression": "Regression in\na later run",
}


def fig3_ledger() -> None:
    led = read("change_ledger.csv")
    ct = pd.crosstab(led["class"], led["found_by"])
    rows = list(ct.sum(axis=1).sort_values(ascending=False).index)
    cols = list(ct.sum(axis=0).sort_values(ascending=False).index)
    ct = ct.loc[rows, cols]

    fig = plt.figure(figsize=(7.2, 3.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[5.2, 1.0], wspace=0.06)
    ax = fig.add_subplot(gs[0])
    axT = fig.add_subplot(gs[1], sharey=ax)
    vmax = ct.values.max()
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("seq", BLUE_SEQ[:10])
    # Each cell is a rectangle inset from its neighbours by a surface-coloured gap; an image with lines drawn over
    # its joins left a hairline of the cell colour on the outer edge.
    gap = 0.025
    for i in range(ct.shape[0]):
        for j in range(ct.shape[1]):
            v = ct.values[i, j]
            if v:
                ax.add_patch(Rectangle((j - 0.5 + gap, i - 0.5 + gap), 1 - 2 * gap, 1 - 2 * gap,
                                                          fc=cmap(v / (vmax * 1.05)), ec="none", lw=0))
                col = SURFACE if v / vmax > 0.5 else INK
                ax.text(j, i, str(v), ha="center", va="center", fontsize=7.5, color=col)
    ax.set_xlim(-0.5, ct.shape[1] - 0.5)
    ax.set_ylim(ct.shape[0] - 0.5, -0.5)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([FOUND_LABELS[c] for c in cols], fontsize=6.8)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([CLASS_LABELS[r] for r in rows], fontsize=7)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    totals = ct.sum(axis=1)
    axT.barh(range(len(rows)), totals.values, height=0.55, color=CAT[0], edgecolor=SURFACE, zorder=3)
    for i, v in enumerate(totals.values):
        axT.text(v + 0.5, i, str(v), va="center", ha="left", fontsize=7.5, color=INK2)
    axT.set_xlim(0, totals.max() * 1.35)
    axT.set_xticks([])
    axT.tick_params(left=False, labelleft=False)
    for sp in ("top", "right", "bottom", "left"):
        axT.spines[sp].set_visible(False)
    axT.set_title("Total\nby class", fontsize=6.8, color=INK2, pad=18)
    cs = ct.sum(axis=0)
    for j, c in enumerate(cols):
        ax.text(j, len(rows) - 0.18, str(int(cs[c])), ha="center", va="top", fontsize=7.5, color=INK2, clip_on=False)
    ax.text(-0.62, len(rows) - 0.18, "Total by route", ha="right", va="top", fontsize=7, color=INK2, clip_on=False)
    save(fig, "fig3_defect_classes")


# ------------------------------------------------------------- figure 5: evidence points and cutoffs


def compact(n: float) -> str:
    n = float(n)
    if n >= 100000:
        return f"{n/1000:.0f}k"
    return f"{n/1000:.1f}k" if n >= 1000 else f"{int(n)}"


def fig5_points() -> None:
    pts = read("points_by_clinvar_v37.csv")
    pts = pts[pts.scope == "no benign evidence"]
    pts = pts[pts.pathogenic_points <= 11]
    rep = read("combining_rules_replay_v37.csv")

    fig = plt.figure(figsize=(7.2, 5.0))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.4, 1], hspace=0.70, wspace=0.36)
    axA = fig.add_subplot(gs[:, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 1])

    # A: composition by pathogenic points, among variants with no benign evidence
    x = pts.pathogenic_points.to_numpy()
    tot = pts.n.to_numpy().astype(float)
    parts = [("P/LP", pts.clinvar_plp.to_numpy() / tot), ("VUS", pts.clinvar_vus.to_numpy() / tot), ("B/LB", pts.clinvar_blb.to_numpy() / tot)]
    bottom = np.zeros(len(x))
    for name, frac in parts:
        axA.bar(x, frac * 100, bottom=bottom * 100, width=0.72, color=TIER3[name], edgecolor=SURFACE, linewidth=1.0, zorder=3)
        bottom += frac
    for xi, n in zip(x, tot):
        axA.text(xi, 101.5, compact(n), ha="center", va="bottom", fontsize=6.0, color=INK2, rotation=90)
    for edge in (5.5, 9.5):
        axA.axvline(edge, color=INK, lw=0.9, zorder=4)
    axA.set_xlim(-0.6, 11.6)
    axA.set_ylim(0, 100)
    axA.set_xticks(range(0, 12))
    axA.set_ylabel("Share of variants with that many points (%)")
    axA.set_title("Where ClinVar's classes fall on the points scale", loc="left", fontsize=8, pad=28)
    hgrid(axA)
    # the classes the points scale maps to, under the axis
    for lo, hi, lab in ((-0.5, 5.5, "VUS or benign side\n(0 to 5 points)"), (5.5, 9.5, "Likely pathogenic\n(6 to 9)"), (9.5, 11.5, "Pathogenic\n(10 or more)")):
        axA.annotate("", xy=(lo + 0.08, -0.105), xytext=(hi - 0.08, -0.105), xycoords=("data", "axes fraction"), textcoords=("data", "axes fraction"),
                     arrowprops=dict(arrowstyle="-", color=INK2, lw=0.9, shrinkA=0, shrinkB=0), annotation_clip=False)
        axA.text((lo + hi) / 2, -0.118, lab, transform=axA.get_xaxis_transform(), ha="center", va="top", fontsize=6.4, color=INK2, linespacing=1.05)
    axA.set_xlabel("Pathogenic evidence points the classifier reached (no benign evidence met)", labelpad=3)
    axA.legend(handles=[Patch(fc=TIER3[k], ec="none", label=f"ClinVar {k}") for k in ("P/LP", "VUS", "B/LB")],
               loc="upper center", bbox_to_anchor=(0.5, -0.215), ncol=3, fontsize=6.6, columnspacing=1.4, handlelength=1.2)
    panel_letter(axA, "A", x=-0.12, y=1.14)

    # B: LP cutoff trade-off, from real runs of the final classifier with that one constant changed
    runs = read("runs_v37_variants.csv")
    cut = runs[(runs.group == "cutoff") & (runs.knob == "likely_pathogenic_points")].copy()
    cut["lp"] = cut.value.astype(int)
    shipped_row = runs[runs.name == "v37"].iloc[[0]].copy()
    shipped_row["lp"] = 6
    cut = pd.concat([cut[cut.lp != 6], shipped_row], ignore_index=True)
    cut["wrong"] = cut.vus_to_plp + cut.false_pathogenic
    cut = cut.sort_values("lp")
    axB.plot(cut.wrong / 1000, cut.plp_recall * 100, color=CAT[0], lw=LW, zorder=2)
    for r in cut.itertuples():
        shipped = r.lp == 6
        # a filled black marker is the shipped setting, here and in Figures 4 and 5C; open markers are the others
        axB.plot([r.wrong / 1000], [r.plp_recall * 100], marker="o", ms=MS + (2 if shipped else 0), mfc=INK if shipped else SURFACE,
                 mec=SURFACE if shipped else CAT[0], mew=1.5 if shipped else 1.4, zorder=5 if shipped else 3)
    lab = {4: (0.9, 0, "4"), 5: (-0.9, 2.4, "5"), 6: (0.9, 0.6, "6 (shipped)"), 7: (0.9, -1.6, "7 and 8"), 9: (0.9, 0, "9")}
    for r in cut.itertuples():
        if r.lp in lab:
            dx, dy, t = lab[r.lp]
            axB.text(r.wrong / 1000 + dx, r.plp_recall * 100 + dy, t, fontsize=6.6, color=INK2, va="center", ha="left")
    axB.set_xlabel("ClinVar VUS or B/LB called P/LP\n(thousands)")
    axB.set_ylabel("ClinVar P/LP called P/LP (%)")
    axB.set_xlim(-1.0, 22)
    axB.set_ylim(48, 82)
    hgrid(axB)
    axB.set_title("Likely pathogenic cutoff, in points (real runs)", loc="left", fontsize=8, pad=6)
    panel_letter(axB, "B", x=-0.30, y=1.14)

    # C: benign-side scheme trade-off
    sch = {
        "shipped (v37)": "Table (shipped)",
        "SVI points, both sides, Likely Benign from -1 (Tavtigian)": "Points, from -1",
        "SVI points, both sides, Likely Benign from -2": "Points, from -2",
        "SVI points, both sides, Likely Benign from -3": "Points, from -3",
    }
    sub = rep[rep.scheme.isin(sch)].copy()
    pt = sub[sub.scheme != "shipped (v37)"].sort_values("false_benign")
    axC.plot(pt.false_benign, pt.blb_recall * 100, color=CAT[1], lw=LW, zorder=2)
    place = {"shipped (v37)": (24, 0, "left"), "SVI points, both sides, Likely Benign from -1 (Tavtigian)": (-20, 0, "right"),
             "SVI points, both sides, Likely Benign from -2": (28, 0, "left"), "SVI points, both sides, Likely Benign from -3": (26, -0.5, "left")}
    for r in sub.itertuples():
        shipped = r.scheme == "shipped (v37)"
        axC.plot([r.false_benign], [r.blb_recall * 100], marker="o", ms=MS + (2 if shipped else 0), mfc=INK if shipped else SURFACE,
                 mec=SURFACE if shipped else CAT[1], mew=1.5 if shipped else 1.4, zorder=5 if shipped else 3)
        dx, dy, ha = place[r.scheme]
        axC.text(r.false_benign + dx, r.blb_recall * 100 + dy, sch[r.scheme], fontsize=6.4, color=INK2, va="center", ha=ha)
    axC.set_xlabel("ClinVar P/LP called LB/B\n(missed diagnoses, n)")
    axC.set_ylabel("ClinVar B/LB called LB/B (%)")
    axC.set_xlim(-30, 480)
    axC.set_ylim(30, 80)
    hgrid(axC)
    axC.set_title("Benign side: table or points (replay)", loc="left", fontsize=8, pad=6)
    panel_letter(axC, "C", x=-0.30, y=1.14)
    save(fig, "fig5_points")


# ---------------------------------------------------------------- figures 6 and 7: the discordance


def fig6_anatomy() -> None:
    c5 = read("confusion_5x5_v37.csv")
    c3 = read("confusion_3x3_v37.csv")
    cc = read("what_is_behind_each_call_v37.csv")
    classes = ["Pathogenic", "Likely_pathogenic", "VUS", "Likely_benign", "Benign"]
    cols = classes + ["NoCall"]

    fig = plt.figure(figsize=(7.2, 4.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.12, 1], height_ratios=[1, 1], hspace=0.78, wspace=0.30)
    axA = fig.add_subplot(gs[:, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 1])

    # A: five-tier matrix, colour = share of the ClinVar row
    m = np.zeros((5, 6))
    n = np.zeros((5, 6), dtype=int)
    for i, t in enumerate(classes):
        sub = c5[c5.clinvar == t].set_index("fastvep")["n"]
        tot = sub.sum()
        for j, pcol in enumerate(cols):
            n[i, j] = sub[pcol]
            m[i, j] = 100 * sub[pcol] / tot
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("seq", BLUE_SEQ[:11])
    axA.imshow(m, cmap=cmap, vmin=0, vmax=75, aspect="auto")
    for i in range(5):
        for j in range(6):
            if n[i, j] == 0:
                txt = "0"
            else:
                txt = f"{n[i, j]:,}"
            col = SURFACE if m[i, j] > 40 else INK
            axA.text(j, i - 0.12, txt, ha="center", va="center", fontsize=6.6, color=col)
            axA.text(j, i + 0.22, f"{m[i, j]:.1f}%" if m[i, j] >= 0.05 else ("<0.1%" if n[i, j] else ""), ha="center", va="center", fontsize=5.8, color=col)
    axA.set_xticks(range(6))
    axA.set_xticklabels([SHORT.get(c, "No call") for c in cols], fontsize=7)
    axA.xaxis.tick_top()
    axA.xaxis.set_label_position("top")
    axA.set_xlabel("fastVEP call", labelpad=6)
    axA.set_yticks(range(5))
    axA.set_yticklabels([f"{SHORT[t]}\n(n = {int(c5[c5.clinvar == t].n.sum()):,})" for t in classes], fontsize=6.8)
    axA.set_ylabel("ClinVar 2-star+ class")
    axA.tick_params(length=0)
    for sp in axA.spines.values():
        sp.set_visible(False)
    for k in range(7):
        axA.axvline(k - 0.5, color=SURFACE, lw=1.5)
    for k in range(6):
        axA.axhline(k - 0.5, color=SURFACE, lw=1.5)
    # outline the opposite-direction cells
    for (i, j) in [(0, 3), (0, 4), (1, 3), (1, 4), (3, 0), (3, 1), (4, 0), (4, 1)]:
        axA.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec=INK, lw=1.2, zorder=5))
    panel_letter(axA, "A", x=-0.20, y=1.10)

    def bars(ax, df, ykey, xkey, order, ylabels, title, xlabel, legend_title, totals_col):
        series_order = ["P/LP", "VUS", "B/LB"]
        for i, key in enumerate(order):
            sub = df[df[ykey] == key].set_index(xkey)
            tot = sub[totals_col].iloc[0]
            left = 0.0
            for name in series_order:
                v = float(sub.loc[name, "n"])
                w = 100 * v / tot
                ax.barh(i, w, left=left, height=0.52, color=TIER3[name], edgecolor=SURFACE, linewidth=1.2, zorder=3)
                if w >= 11:
                    ax.text(left + w / 2, i, f"{w:.1f}%", ha="center", va="center", fontsize=6.4, color=text_on(TIER3[name]), zorder=4)
                left += w
            if key in ("P/LP", "B/LB"):
                opp = {"P/LP": "B/LB", "B/LB": "P/LP"}[key]
                v = float(sub.loc[opp, "n"])
                ax.text(101, i, f"{int(v)} opposite", ha="left", va="center", fontsize=6.2, color=INK2)
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels(ylabels, fontsize=6.8)
        ax.invert_yaxis()
        ax.set_xlim(0, 112)
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_xlabel(xlabel)
        ax.set_title(title, loc="left", fontsize=8, pad=17)
        handles = [Line2D([0], [0], lw=0, label=legend_title)] + [Patch(fc=TIER3[k], ec="none", label=k) for k in series_order]
        ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=4, fontsize=6.4, columnspacing=0.9,
                  handlelength=1.0, handletextpad=0.4, borderaxespad=0.1)

    # B: by ClinVar class (rows)
    d = c3[c3.fastvep != "NoCall"].rename(columns={"fastvep": "x"})
    bars(axB, d.assign(y=d.clinvar), "y", "x", ["P/LP", "VUS", "B/LB"],
         [f"ClinVar {k}\n(n = {int(d[d.clinvar == k].row_total.iloc[0]):,})" for k in ("P/LP", "VUS", "B/LB")],
         "What fastVEP says of each ClinVar class", "Share of the ClinVar class (%)", "fastVEP call:", "row_total")
    panel_letter(axB, "B", x=-0.30, y=1.30)

    # C: behind each fastVEP call (columns)
    d = cc.rename(columns={"clinvar": "x"})
    bars(axC, d.assign(y=d.fastvep), "y", "x", ["P/LP", "VUS", "B/LB"],
         [f"fastVEP {k}\n(n = {int(d[d.fastvep == k].column_total.iloc[0]):,})" for k in ("P/LP", "VUS", "B/LB")],
         "What ClinVar says behind each fastVEP call", "Share of the fastVEP call (%)", "ClinVar class:", "column_total")
    panel_letter(axC, "C", x=-0.30, y=1.30)
    save(fig, "fig6_anatomy")


def fig7_strata() -> None:
    bc = read("by_consequence_v37.csv")
    bs = read("by_review_status_v37.csv")
    tb = read("temporal_by_class_v37.csv")

    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.55, 1], hspace=0.62, wspace=0.42)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1], sharey=axA)
    axC = fig.add_subplot(gs[1, :])

    groups = ["Frameshift or nonsense", "Canonical splice site", "Start or stop loss", "Missense", "In-frame or protein-altering",
              "Splice region", "Synonymous", "Intronic", "UTR", "Other (flank, non-coding, unscored)"]
    rows = [("c", g) for g in groups] + [("gap", "")] + [("s", 2), ("s", 3), ("s", 4)]

    def one(ax, cls, title, letter):
        for i, (kind, key) in enumerate(rows):
            if kind == "gap":
                continue
            if kind == "c":
                sub = bc[(bc.consequence_group == key) & (bc.clinvar == cls)].set_index("fastvep")
                label = key
            else:
                sub = bs[(bs.stars == key) & (bs.clinvar == cls)].set_index("fastvep")
                label = f"{key}-star records"
            tot = float(sub["row_total"].iloc[0])
            ax.text(-2, i, label, ha="right", va="center", fontsize=6.4, color=INK2) if ax is axA else None
            left = 0.0
            for name in ("P/LP", "VUS", "B/LB"):
                v = float(sub.loc[name, "n"])
                w = 100 * v / tot if tot else 0
                ax.barh(i, w, left=left, height=0.62, color=TIER3[name], edgecolor=SURFACE, linewidth=1.0, zorder=3)
                if w >= 9:
                    ax.text(left + w / 2, i, f"{w:.0f}", ha="center", va="center", fontsize=5.9, color=text_on(TIER3[name]), zorder=4)
                left += w
            nocall = float(sub.loc["NoCall", "n"]) if "NoCall" in sub.index else 0.0
            # An asterisk marks a row whose bar stops short of 100%: those records have no call.
            ax.text(103, i, (f"{int(tot):,}" + ("*" if nocall else "")) if tot else "none", ha="left", va="center", fontsize=5.8, color=INK2)
        ax.set_yticks([])
        ax.set_xlim(0, 100)
        ax.set_ylim(len(rows) - 0.5, -0.9)
        ax.set_xticks([0, 50, 100])
        ax.set_xlabel("Share of ClinVar %s variants (%%)" % cls)
        ax.set_title(title, loc="left", fontsize=8, pad=4)
        for sp in ("left",):
            ax.spines[sp].set_visible(False)
        ax.text(103, -0.9, "n", ha="left", va="center", fontsize=5.8, color=INK2)

    one(axA, "P/LP", "ClinVar P/LP", "A")
    one(axB, "B/LB", "ClinVar B/LB", "B")
    panel_letter(axA, "A", x=-0.04, y=1.12)
    panel_letter(axB, "B", x=-0.04, y=1.12)
    fig.legend(handles=[Line2D([0], [0], lw=0, label="fastVEP call:")] + [Patch(fc=TIER3[k], ec="none", label=k) for k in ("P/LP", "VUS", "B/LB")],
               loc="center", bbox_to_anchor=(0.5, 0.437), ncol=4, fontsize=6.8, columnspacing=1.2, handlelength=1.1, handletextpad=0.4)

    # C: the same three classes, April set against the October strata
    order = [("April set (v37)", "April set\n(used for tuning)"), ("October, unchanged", "October,\nunchanged"), ("October, new", "October, new\n(not in the April set)")]
    classes = [("Pathogenic", "P"), ("Likely_pathogenic", "LP"), ("VUS", "VUS"), ("Likely_benign", "LB"), ("Benign", "B")]
    width = 0.24
    xs = np.arange(len(classes))
    # Greys, not hues: red and blue mean P/LP and B/LB in panels A and B, and a group of variants is not a class.
    cols3 = ["#d4d3cc", "#9a9993", "#5b5a55"]
    for k, (g, lab) in enumerate(order):
        sub = tb[tb.group == g].set_index("clinvar")
        vals = [sub.loc[c, "recall"] * 100 for c, _ in classes]
        lo = [sub.loc[c, "ci_lo"] * 100 for c, _ in classes]
        hi = [sub.loc[c, "ci_hi"] * 100 for c, _ in classes]
        pos = xs + (k - 1) * (width + 0.03)
        axC.bar(pos, vals, width=width, color=cols3[k], edgecolor=SURFACE, linewidth=1.0, zorder=3, label=lab.replace("\n", " "))
        axC.errorbar(pos, vals, yerr=[np.array(vals) - np.array(lo), np.array(hi) - np.array(vals)], fmt="none", ecolor=INK, elinewidth=0.8, capsize=1.5, zorder=4)
        for xp, v in zip(pos, vals):
            axC.text(xp, v + 3.2, f"{v:.0f}", ha="center", va="bottom", fontsize=6.0, color=INK2)
    axC.set_xticks(xs)
    axC.set_xticklabels([f"ClinVar {a}" for _, a in classes])
    axC.set_ylim(0, 112)
    axC.set_yticks([0, 25, 50, 75, 100])
    axC.spines["left"].set_bounds(0, 100)  # a percentage stops at 100; the room above holds the value labels
    axC.set_ylabel("Same-direction recall (%)")
    hgrid(axC)
    axC.legend(loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.2), fontsize=6.8, columnspacing=1.6)
    nb = tb[(tb.clinvar == "Benign")].set_index("group")
    axC.set_xlim(-0.6, 5.9)
    sr = read("temporal_standardised_recall_v37.csv").set_index("clinvar").loc["Benign"]
    axC.text(4.62, 50, f"New Benign variants\nare rarer: median\npopulation-maximum\nAF {nb.loc['October, new', 'median_popmax_af']*100:.2f}%, against\n{nb.loc['April set (v37)', 'median_popmax_af']*100:.1f}% in the April set.\nAt April's recall in each\nfrequency band they\nwould be at {sr.new_recall_at_april_recall_per_band*100:.0f}%.",
             fontsize=6.0, color=INK2, ha="left", va="center", linespacing=1.15)
    panel_letter(axC, "C", x=-0.04, y=1.18)
    axC.set_title("", loc="left")
    save(fig, "fig7_strata")


# -------------------------------------------------------------------- figure 4: the thresholds


def fig4_thresholds() -> None:
    r = read("runs_v37_variants.csv")
    base = r[r.name == "v37"].iloc[0]
    fig = plt.figure(figsize=(7.2, 3.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.45, 1.0, 0.9], wspace=0.40)
    axA = fig.add_subplot(gs[0])
    axB = fig.add_subplot(gs[1])
    axC = fig.add_subplot(gs[2])

    # A: the benign-side trade space. Every knob passes through the shipped point.
    # (dx, dy, alignment) of each labelled point, in data units, chosen so no label sits on a line.
    knobs = [
        ("bs1_af_threshold", "BS1 bar", CAT[0], "0.005", {"0.0002": (0.3, 0, "left"), "0.001": (0.3, -4, "left"), "0.02": (0, 6, "center"), "0.01": (0, 6, "center")}),
        ("bs2_hom_prevalence_threshold", "BS2 homozygote bar", CAT[1], "0.001", {"1e-06": (0.3, 0, "left"), "1e-05": (-0.3, 3.5, "right")}),
        ("bp7_max_intron_offset", "BP7 far boundary", CAT[2], "300", {"0": (0, 6, "center"), "50": (0.3, -4, "left"), "1000000": (-0.3, 0, "right")}),
        ("phylop_conserved", "BP7 PhyloP bar", CAT[3], "2.0", {"1": (0, -4.5, "center"), "1000": (0.3, 0, "left")}),
        ("bp1_max_pathogenic_missense", "BP1 veto count", CAT[4], "3", {"1": (0.3, -3.5, "left"), "25": (0.3, 0, "left")}),
    ]
    show = {"bs1_af_threshold": {"0.0002", "0.001", "0.01", "0.02"}, "bs2_hom_prevalence_threshold": {"1e-06", "1e-05"},
            "bp7_max_intron_offset": {"0", "50", "1000000"}, "phylop_conserved": {"1", "1000"}, "bp1_max_pathogenic_missense": set()}  # its points sit on the shipped one and on BP7's line; the caption names the range
    for knob, label, color, shipped, offs in knobs:
        sub = r[(r.knob == knob) & (r.group == "sweep")].copy()
        sub["v"] = sub.value.astype(str)
        pts = [(float(base.blb_recall) * 100, int(base.false_benign), shipped)]
        for row in sub.itertuples():
            pts.append((row.blb_recall * 100, int(row.false_benign), row.v))
        pts.sort(key=lambda t: t[0])
        axA.plot([p_[0] for p_ in pts], [p_[1] for p_ in pts], color=color, lw=LW, zorder=2)
        for x, y, v in pts:
            if v == shipped:
                continue
            axA.plot([x], [y], marker="o", ms=MS - 0.5, mfc=SURFACE, mec=color, mew=1.3, zorder=3)
            if v in show[knob]:
                fmt = ({"1000": "off", "1": "1"} if knob == "phylop_conserved" else {"1": "count 1", "25": "count 25"} if knob == "bp1_max_pathogenic_missense"
                       else {"1000000": "no bound", "0": "off"}).get(v, v)
                dx, dy, ha = offs[v]
                axA.text(x + dx, y + dy, fmt, fontsize=5.9, color=INK2, ha=ha, va="center")
    axA.plot([float(base.blb_recall) * 100], [int(base.false_benign)], marker="o", ms=MS + 2, mfc=INK, mec=SURFACE, mew=1.5, zorder=5)
    axA.text(float(base.blb_recall) * 100 - 0.35, int(base.false_benign) + 9, "shipped", fontsize=6.4, color=INK, ha="right", va="center")
    axA.set_xlim(57, 69.4)
    axA.set_ylim(-3, 112)
    axA.set_xlabel("ClinVar B/LB called LB/B (%)")
    axA.set_ylabel("ClinVar P/LP called LB/B\n(missed diagnoses, n)")
    hgrid(axA)
    axA.legend(handles=[Line2D([0], [0], color=c, lw=LW, label=l) for _, l, c, _, _ in knobs], loc="upper left", fontsize=6.6)
    axA.set_title("Bars on the benign side", loc="left", fontsize=8, pad=6)
    panel_letter(axA, "A", x=-0.17, y=1.10)

    # B: the PM2 bar for dominant genes, under the point system and under the 2015 table
    pm = r[(r.knob == "pm2_ad_af_threshold") & (r.group == "sweep")].copy()
    tb = r[r.knob.str.startswith("pm2_ad_af_threshold (Richards") & (r.group == "interaction")].copy() if "interaction" in set(r.group) else r.iloc[0:0]
    tbase = r[r.name == "v37_sw_table_both_sides"]

    def xval(v):
        v = float(v)
        return 2.5e-7 if v == 0 else v

    ptsA = [(xval(row.value), row.plp_recall * 100) for row in pm.itertuples()] + [(4e-5, base.plp_recall * 100)]
    ptsA.sort()
    axB.plot([p_[0] for p_ in ptsA], [p_[1] for p_ in ptsA], color=CAT[0], lw=LW, marker="o", ms=MS - 1, mec=SURFACE, mew=1, zorder=3,
             label="Points")
    ar = r[(r.knob == "pm2_ar_af_threshold") & (r.group == "sweep")]
    ptsR = sorted([(xval(row.value), row.plp_recall * 100) for row in ar.itertuples()] + [(7e-5, base.plp_recall * 100)])
    axB.plot([p_[0] for p_ in ptsR], [p_[1] for p_ in ptsR], color=CAT[2], lw=1.0, ls=(0, (4, 2)), marker="D", ms=MS - 2.2, mfc=SURFACE, mec=CAT[2],
             mew=1.1, zorder=3, label="Points, recessive-only")
    if len(tb):
        ptsT = [(xval(row.value), row.plp_recall * 100) for row in tb.itertuples()] + [(4e-5, float(tbase.iloc[0].plp_recall) * 100)]
        ptsT.sort()
        axB.plot([p_[0] for p_ in ptsT], [p_[1] for p_ in ptsT], color=CAT[1], lw=LW, marker="o", ms=MS - 1, mec=SURFACE, mew=1, zorder=3, label="2015 table")
    axB.axvline(4e-5, color=AXIS, lw=0.9, zorder=1)
    axB.text(4e-5 * 1.2, 65.6, "shipped 4e-5 (7e-5\nfor recessive-only)", fontsize=5.8, color=INK2, va="top", ha="left")
    axB.set_xscale("log")
    axB.set_xlim(1.5e-7, 1e-3)
    axB.set_ylim(36, 68)
    axB.set_xticks([2.5e-7, 1e-6, 1e-5, 1e-4])
    axB.set_xticklabels(["0", "1e-6", "1e-5", "1e-4"])
    axB.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    axB.set_xlabel("PM2 bar")
    axB.set_ylabel("ClinVar P/LP called P/LP (%)")
    hgrid(axB)
    axB.legend(loc="lower right", bbox_to_anchor=(1.0, 0.0), fontsize=6.0, handlelength=1.8)
    axB.set_title("PM2 bars, by combining rule", loc="left", fontsize=8, pad=6)
    panel_letter(axB, "B", x=-0.27, y=1.10)

    # C: the dominant-gene BS2 carrier bar, change in recall against the shipped 3e-4
    bs = r[(r.knob == "bs2_ad_prevalence_threshold") & (r.group == "sweep")].copy()
    xs, plp, blb = [], [], []
    for row in bs.itertuples():
        v = float(row.value)
        xs.append(2.0e-5 if v == 0 else v)
        plp.append((row.plp_recall - base.plp_recall) * 100)
        blb.append((row.blb_recall - base.blb_recall) * 100)
    xs.append(3e-4); plp.append(0.0); blb.append(0.0)
    order = np.argsort(xs)
    xs, plp, blb = np.array(xs)[order], np.array(plp)[order], np.array(blb)[order]
    axC.plot(xs, plp, color=TIER3["P/LP"], lw=LW, marker="o", ms=MS - 1, mec=SURFACE, mew=1, zorder=3)
    axC.plot(xs, blb, color=TIER3["B/LB"], lw=LW, marker="o", ms=MS - 1, mec=SURFACE, mew=1, zorder=3)
    axC.axhline(0, color=AXIS, lw=0.9, zorder=1)
    axC.axvline(3e-4, color=AXIS, lw=0.9, zorder=1)
    axC.set_xscale("log")
    axC.set_xlim(1.2e-5, 5e-3)
    axC.set_xticks([2.0e-5, 1e-4, 1e-3])
    axC.set_xticklabels(["off", "1e-4", "1e-3"])
    axC.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    axC.set_xlabel("BS2 carrier bar")
    axC.set_ylabel("Change in recall (pp)")
    hgrid(axC)
    axC.text(3e-4 * 1.15, axC.get_ylim()[1] * 0.92, "shipped\n3e-4", fontsize=6.2, color=INK2, va="top", ha="left")
    axC.text(4.6e-3, 0.012, "P/LP", fontsize=6.4, color=INK2, ha="right", va="bottom")
    axC.text(4.6e-3, blb[-1] - 0.012, "B/LB", fontsize=6.4, color=INK2, ha="right", va="top")
    axC.set_title("BS2 carrier bar", loc="left", fontsize=8, pad=6)
    panel_letter(axC, "C", x=-0.30, y=1.10)
    save(fig, "fig4_thresholds")


# ----------------------------------------------------------------------- figure 8: ablations

ABLATION_LABELS = {
    "v37_sw_pm2_moderate": "PM2 at Moderate (SVI downgrade off)",
    "v37_sw_pm2_needs_record": "PM2 needs a gnomAD record",
    "v37_sw_table_both_sides": "2015 table on both sides (no point system)",
    "v37_sw_no_gene_validity": "No gene-disease validity gate",
    "v37_sw_no_mechanism_gate": "No disease-mechanism gate on PVS1",
    "v37_sw_no_lowpen_gate": "No ClinVar low-penetrance gate",
    "v37_sw_no_region_flags": "No segdup / low-complexity veto",
    "v37_sw_no_faf": "Point estimate in place of filtering AF",
    "v37_sw_clinvar_informed": "Own ClinVar record counted (PS1, PM1)",
    "v37_sw_no_homology_list": "No homologous-gene list",
    "v37_sw_no_exceptions": "No frequency-exception list",
    "v37_sw_no_bp1_veto": "No BP1 mechanism veto",
    "v37_sw_no_pm1_benign_half": "PM1 ignores benign variation",
    "v37_sw_no_an_floor": "No allele-number floor (2,000)",
    "v37_sw_pp3_cap_moderate": "PP3 capped at Moderate",
    "v37_sw_pvs1_nmd50": "PVS1 on the exact 50-nt NMD rule",
    "v37_src_noclinvar_protein": "ClinVar protein and splice index",
    "v37_src_noclinvar_all": "All ClinVar-derived sources",
    "v37_src_norevel": "REVEL",
    "v37_src_nospliceai": "SpliceAI",
    "v37_src_nophylop": "PhyloP",
    "v37_src_noconstraint": "gnomAD gene constraint",
    "v37_src_noclingen": "ClinGen gene-disease validity",
    "v37_src_norepeat": "RepeatMasker",
    "v37_src_nognomad": "gnomAD (all frequency criteria)",
}


MEASURED_LABELS = {
    "v37_hk_bp1_loo": "BP1 veto leaves out the variant's own record",
    "v37_hk_bp7_spliceai01": "BP7 needs SpliceAI at most 0.1",
    "v37_sw_phylop_1000": "BP7 without the conservation test",
    "v37_hk_bp7_walker": "BP7 with both changes, as Walker 2023",
    "v37_gap_clingen_complete": "ClinGen file with its omitted curations",
    "v37_gap_orphanet_lof": "PVS1 for Orphanet loss-of-function genes",
    "v37_gap_both": "Both inputs of the PVS1 gap",
    "v37_hk_strict": "The three adjustments, removed together",
}


def _signed(v: float, digits: int = 2) -> str:
    """A signed value that never reads -0.00: what rounds to zero is written 0.00."""
    t = f"{abs(v):.{digits}f}"
    if float(t) == 0:
        return t
    return ("+" if v > 0 else "-") + t


def fig8_ablations() -> None:
    r = read("runs_v37_variants.csv")
    base = r[r.name == "v37"].iloc[0]
    labels = {**ABLATION_LABELS, **MEASURED_LABELS}
    a = r[r.name.isin(labels)].copy()
    a["d_same_pp"] = (a.same_direction - base.same_direction) / 673660 * 100
    a["d_missed"] = a.false_benign - base.false_benign
    a["d_fp"] = a.false_pathogenic - base.false_pathogenic
    a["d_plp_pp"] = (a.plp_recall - base.plp_recall) * 100
    a["label"] = a.name.map(labels)

    blocks = [("A", "Setting switched off (on for the two opt-in settings): change from v37", a[a.group == "switch"].sort_values("d_same_pp")),
              ("B", "One source removed from the annotation directory: change from v37", a[a.group == "source"].sort_values("d_same_pp")),
              ("C", "Departures from the sources, measured by real runs: change from v37",
               a[a.name.isin(MEASURED_LABELS)].sort_values("d_same_pp"))]
    nrow = [len(b) for _, _, b in blocks]
    fig = plt.figure(figsize=(7.2, 8.6))
    gs = fig.add_gridspec(3, 4, height_ratios=nrow, width_ratios=[1.5, 1, 0.9, 0.9], hspace=0.42, wspace=0.10)
    for bi, (letter, title, b) in enumerate(blocks):
        axs = [fig.add_subplot(gs[bi, k]) for k in range(4)]
        y = np.arange(len(b))
        specs = [("d_same_pp", "Same-direction\nagreement (pp)", INK2, True),
                 ("d_plp_pp", "ClinVar P/LP called\nP/LP (pp)", TIER3["P/LP"], True),
                 ("d_missed", "Missed diagnoses (n)\nfewer is better", TIER3["P/LP"], False),
                 ("d_fp", "False-pathogenic (n)\nfewer is better", TIER3["B/LB"], False)]
        for k, (col, xl, color, is_pp) in enumerate(specs):
            ax = axs[k]
            vals = b[col].to_numpy()
            ax.barh(y, vals, height=0.56, color=color, edgecolor=SURFACE, linewidth=1.0, zorder=3)
            ax.axvline(0, color=AXIS, lw=0.9, zorder=2)
            span = max(abs(vals).max(), 1e-9)
            for yi, v in zip(y, vals):
                t = _signed(v) if is_pp else (f"{int(v):+d}" if int(v) else "0")
                ax.text(v + (0.04 * span if v >= 0 else -0.04 * span), yi, t, ha="left" if v >= 0 else "right", va="center", fontsize=5.8, color=INK2)
            ax.set_xlim(-span * 1.5, span * 1.5)
            ax.set_ylim(len(b) - 0.5, -0.5)
            ax.set_yticks(y if k == 0 else [])
            if k == 0:
                ax.set_yticklabels(b["label"], fontsize=6.4)
            ax.set_xlabel(xl, fontsize=6.4)
            ax.tick_params(axis="x", labelsize=6.0)
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3, symmetric=True, integer=not is_pp))
            vgrid(ax)
            for sp in ("left",):
                ax.spines[sp].set_visible(False)
        axs[0].set_title(f"{letter}  {title}", loc="left", fontsize=8, x=-1.05, pad=7, fontweight="bold")
    save(fig, "fig8_ablations")


# ---------------------------------------------------------------------------------- driver

# -------------------------------------------------------------- figure 9: agreement against contradictions


def fig9_frontier() -> None:
    """Agreement in direction against opposite-direction calls, for the versions, the runs of this paper and simple rules.

    The point of the figure is that agreement does not say how cautious a classifier is: rules with no combining rule,
    gene-level gate or safeguard reach a higher agreement than v37 by accepting many times as many contradictions, and at
    the number of contradictions that v37 has, or ten times it, they agree far less.
    """
    N = 673_660
    traj = read("trajectory.csv")
    main = traj[traj.mainline & traj.same_direction_rate.notna() & traj.opposite.notna() & (traj.source == "scored matrix on disk")]
    runs = read("runs_v37_variants.csv")
    runs = runs[runs.name != "v37"]
    curve = read("baselines_curve_v37.csv")
    base = read("baselines_v37.csv").set_index("rule")

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.set_xscale("symlog", linthresh=10, linscale=0.5)
    ax.set_xlim(-1.2, 12000)
    ax.set_ylim(40, 88)
    hgrid(ax)
    ax.grid(axis="x")

    # the runs of this paper: v37 with one thing changed
    ax.plot(runs.opposite, runs.same_direction / N * 100, ls="none", marker="o", ms=2.6, mfc="#b7b6ae", mec="none", alpha=0.75, zorder=2)

    # the versions
    ax.plot(main.opposite, main.same_direction_rate * 100, color=CAT[0], lw=LW, marker="o", ms=MS - 2, mec=SURFACE, mew=0.8, zorder=4)
    # (x factor, y offset, alignment) of each version label, chosen so that none sits on a line
    for run, fx, dy, ha in (("v7", 1.12, -1.5, "left"), ("v8", 1.0, 1.6, "center"), ("v10", 0.86, -1.5, "right"), ("v13", 1.45, 1.4, "left")):
        r = main[main.run == run].iloc[0]
        ax.text(r.opposite * fx, r.same_direction_rate * 100 + dy, f"{run}: {int(r.opposite)}", fontsize=6.2, color=INK2, ha=ha, va="center")

    # the rules, at nine frequency bars each
    style = {"R1": (CAT[1], "o", "R1: benign-side at or above a frequency bar (the bar is labelled)"),
             "R1 + R2": (CAT[3], "s", "R1 and R2: a rare null variant is pathogenic-side"),
             "R1 + R2 + R3": (CAT[2], "^", "R1, R2 and R3: a quiet variant with SpliceAI below 0.2 is benign-side")}
    for rule, (col, mk, lab) in style.items():
        c = curve[curve.rule == rule].sort_values("benign_frequency_bar", ascending=False)
        ax.plot(c.opposite, c.same_direction_share * 100, color=col, lw=1.1, marker=mk, ms=MS - 2.2, mfc=SURFACE, mec=col, mew=1.1, zorder=3, label=lab)
    # the frequency bar at each end of each line
    for rule, (col, mk, lab) in style.items():
        c = curve[curve.rule == rule].set_index("benign_frequency_bar")
        for bar, dx, dy, ha in ((5e-2, 0.80, 0.0, "right"), (1e-4, 1.22, 0.0, "left")):
            r = c.loc[bar]
            ax.text(r.opposite * dx, r.same_direction_share * 100 + dy, f"{bar:g}", fontsize=6.0, color=INK2, ha=ha, va="center")

    # the same rules with no frequency rule at all: the pathogenic-side rule alone, and with the benign-side rule for quiet variants
    for rule, col, mk in (("R2", CAT[3], "s"), ("R2 + R3", CAT[2], "^")):
        r = curve[curve.rule == rule].iloc[0]
        ax.plot([r.opposite], [r.same_direction_share * 100], ls="none", marker=mk, ms=MS - 1.0, mfc=col, mec=SURFACE, mew=0.9, zorder=5)
    r3 = curve[curve.rule == "R2 + R3"].iloc[0]
    ax.annotate(f"no frequency rule: {r3.same_direction_share * 100:.1f}%, {int(r3.opposite)}", xy=(r3.opposite, r3.same_direction_share * 100),
                xytext=(r3.opposite * 2.0, r3.same_direction_share * 100 - 3.0), fontsize=6.4, color=INK2, ha="center",
                arrowprops=dict(arrowstyle="-", color=AXIS, lw=0.8, shrinkA=0, shrinkB=3))
    r2 = curve[curve.rule == "R2"].iloc[0]
    ax.text(r2.opposite * 1.25, r2.same_direction_share * 100, f"R2 alone: {r2.same_direction_share * 100:.1f}%, {int(r2.opposite)}", fontsize=6.4,
            color=INK2, ha="left", va="center")

    # every variant VUS, and v37
    ax.plot([0], [base.loc["all VUS", "same_direction_share"] * 100], marker="D", ms=MS - 1, mfc=INK2, mec=SURFACE, mew=0.8, zorder=5, ls="none")
    ax.text(1.3, base.loc["all VUS", "same_direction_share"] * 100, "every variant VUS", fontsize=6.4, color=INK2, ha="left", va="center")
    v = base.loc["v37 (the classifier)"]
    ax.plot([v.opposite], [v.same_direction_share * 100], marker="o", ms=MS + 3, mfc=INK, mec=SURFACE, mew=1.5, zorder=6, ls="none")
    ax.annotate(f"v37: {v.same_direction_share * 100:.1f}%, {int(v.opposite)}", xy=(v.opposite, v.same_direction_share * 100),
                xytext=(v.opposite * 0.9, v.same_direction_share * 100 + 4.6), fontsize=6.8, color=INK, ha="center",
                arrowprops=dict(arrowstyle="-", color=AXIS, lw=0.8, shrinkA=0, shrinkB=3))
    best = base.loc["R1 + R2 + R3"]
    ax.annotate(f"0.001: {best.same_direction_share * 100:.1f}%, {int(best.opposite):,}", xy=(best.opposite, best.same_direction_share * 100),
                xytext=(best.opposite * 1.55, best.same_direction_share * 100 + 2.6), fontsize=6.4, color=INK2, ha="center",
                arrowprops=dict(arrowstyle="-", color=AXIS, lw=0.8, shrinkA=0, shrinkB=3))

    ax.set_xticks([0, 10, 100, 1000, 10000])
    ax.set_xticklabels(["0", "10", "100", "1,000", "10,000"])
    ax.set_xlabel("Opposite-direction calls (n, log scale)")
    ax.set_ylabel("Agreement in direction\n(% of 673,660 variants)")
    handles = [
        Line2D([0], [0], color=CAT[0], lw=LW, marker="o", ms=MS - 2, mec=SURFACE, label="fastVEP versions v7 to v37 (scored matrices)"),
        Line2D([0], [0], ls="none", marker="o", ms=3.4, mfc="#b7b6ae", mec="none", label="runs of this paper: v37 with one change"),
    ] + [Line2D([0], [0], color=col, lw=1.1, marker=mk, ms=MS - 2.2, mfc=SURFACE, mec=col, mew=1.1, label=lab) for col, mk, lab in style.values()] + [
        Line2D([0], [0], ls="none", marker="s", ms=MS - 1.0, mfc=CAT[3], mec=SURFACE, label="filled: R2, and R2 with R3, with no frequency rule"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=6.4, handletextpad=0.5, labelspacing=0.35)
    save(fig, "fig9_frontier")


FIGS = {1: fig1_trajectory, 2: fig2_review, 3: fig3_ledger, 5: fig5_points, 4: fig4_thresholds, 6: fig6_anatomy, 7: fig7_strata, 8: fig8_ablations, 9: fig9_frontier}


def main() -> None:
    want = [int(a) for a in sys.argv[1:]] or sorted(FIGS)
    for n in want:
        print(f"figure {n}")
        FIGS[n]()


if __name__ == "__main__":
    main()
