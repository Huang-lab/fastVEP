#!/usr/bin/env python3
"""Diagnostic figures for the ClinVar 2-star+ benchmark.

Focuses on **which ACMG criteria fire for which classification outcome**,
especially for discordant cases — rather than v1↔v7 baseline comparisons
which are documented in METHODS.md / RUN_VERSIONS.md.

Reads outputs from `03_evaluate_concordance.py` under
`data/benchmark/output_v7/` (override via `argv[1]`) and emits the
following figures under `<output_dir>/figures/`:

  fig_concordance_matrix          row-normalised truth × predicted heatmap
  fig_outcome_breakdown           per-truth-class outcome stack
                                  (same / off / VUS / opposite)
  fig_criterion_fire_heatmap      criterion × truth-class fire rate;
                                  the core "which criteria are diagnostic
                                  for which class" diagnostic
  fig_criterion_signatures_by_class
                                  top criterion combinations per truth
                                  class, colored by predicted outcome
  fig_lost_to_vus_signatures      criterion patterns of P/LP→VUS
                                  variants — the dominant failure mode
                                  driven by missing manual-curation
                                  evidence (PS3 / PP1 / PP4)
  fig_opposite_direction_signatures
                                  for the ~100 opposite-direction cases,
                                  the criterion patterns of
                                  Pathogenic-truth→Benign-pred and the
                                  reverse direction

Usage:
  generate_figures.py                     # uses ../../data/benchmark/output_v7
  generate_figures.py <out_dir>
"""

from __future__ import annotations

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

plt.rcParams.update(
    {
        "font.size": 12,
        "axes.titlesize": 16,
        "axes.labelsize": 14,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 11,
        "legend.title_fontsize": 12,
    }
)

C = {
    "P": "#dc2626",
    "LP": "#f97316",
    "VUS": "#6b7280",
    "LB": "#3b82f6",
    "B": "#10b981",
    "delta_up": "#10b981",
    "delta_down": "#ef4444",
    "neutral": "#94a3b8",
}
CLASSES = ["Pathogenic", "Likely_pathogenic", "VUS", "Likely_benign", "Benign"]
CLASS_SHORT = {
    "Pathogenic": "P",
    "Likely_pathogenic": "LP",
    "VUS": "VUS",
    "Likely_benign": "LB",
    "Benign": "B",
}
CLASS_COLOR = {c: C[CLASS_SHORT[c]] for c in CLASSES}


# ──────────────────────────────────────────────────────────────────────
# data loaders
# ──────────────────────────────────────────────────────────────────────

def read_matrix(out_dir: Path):
    rows = list(csv.reader((out_dir / "concordance_matrix.csv").open()))
    header = rows[0]
    matrix = {}
    for row in rows[1:]:
        truth = row[0]
        for col, val in zip(header[1:], row[1:]):
            matrix[(truth, col)] = int(val)
    return matrix


def read_criterion_fires(out_dir: Path) -> dict[str, dict[str, int]]:
    fires = {}
    with (out_dir / "criterion_firing_rates.csv").open() as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            code = row["criterion"]
            fires[code] = {tcl: int(row.get(f"{tcl}_fired", 0)) for tcl in CLASSES}
    return fires


def read_discrepancies(out_dir: Path) -> list[dict[str, str]]:
    out = []
    with (out_dir / "discrepancies.tsv").open() as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for row in rdr:
            out.append(row)
    return out


def class_totals(matrix) -> dict[str, int]:
    cols = CLASSES + ["NoCall"]
    return {t: sum(matrix.get((t, c), 0) for c in cols) for t in CLASSES}


# ──────────────────────────────────────────────────────────────────────
# figures
# ──────────────────────────────────────────────────────────────────────

def fig_concordance_matrix(matrix, fig_dir: Path):
    fig, ax = plt.subplots(figsize=(11, 7.5))
    cols = CLASSES + ["NoCall"]
    mat = np.zeros((len(CLASSES), len(cols)))
    for i, t in enumerate(CLASSES):
        row_total = sum(matrix.get((t, c), 0) for c in cols)
        for j, c in enumerate(cols):
            mat[i, j] = 100 * matrix.get((t, c), 0) / row_total if row_total else 0

    im = ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(
        [CLASS_SHORT.get(c, c) for c in cols], fontsize=13, fontweight="bold"
    )
    ax.set_yticks(range(len(CLASSES)))
    ax.set_yticklabels([CLASS_SHORT[c] for c in CLASSES], fontsize=13, fontweight="bold")
    ax.set_xlabel("fastVEP predicted", fontweight="bold")
    ax.set_ylabel("ClinVar 2-star+ truth", fontweight="bold")
    ax.set_title("ACMG concordance matrix (row-normalised %)", fontweight="bold")

    for i, t in enumerate(CLASSES):
        row_total = sum(matrix.get((t, c), 0) for c in cols)
        for j, c in enumerate(cols):
            cnt = matrix.get((t, c), 0)
            pct = mat[i, j]
            color = "white" if pct > 50 else "black"
            label = f"{cnt:,}\n({pct:.0f}%)"
            ax.text(
                j, i, label, ha="center", va="center", fontsize=10,
                color=color, fontweight="bold" if cols[j] == t else "normal",
            )
        if i < len(cols):
            ax.add_patch(
                plt.Rectangle(
                    (i - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=C["delta_up"], lw=2.5
                )
            )

    cbar = plt.colorbar(im, ax=ax, shrink=0.85)
    cbar.set_label("Row %")
    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(fig_dir / f"fig_concordance_matrix.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  fig_concordance_matrix")


def fig_outcome_breakdown(matrix, fig_dir: Path):
    fig, ax = plt.subplots(figsize=(11, 5.5))
    n_per = class_totals(matrix)

    def share(truth, mask):
        n = n_per[truth] or 1
        return 100 * sum(matrix.get((truth, c), 0) for c in mask) / n

    same_dir_mask = {
        "Pathogenic": ["Pathogenic", "Likely_pathogenic"],
        "Likely_pathogenic": ["Pathogenic", "Likely_pathogenic"],
        "VUS": ["VUS"],
        "Likely_benign": ["Likely_benign", "Benign"],
        "Benign": ["Likely_benign", "Benign"],
    }
    opp_mask = {
        "Pathogenic": ["Benign", "Likely_benign"],
        "Likely_pathogenic": ["Benign", "Likely_benign"],
        "VUS": [],
        "Likely_benign": ["Pathogenic", "Likely_pathogenic"],
        "Benign": ["Pathogenic", "Likely_pathogenic"],
    }
    same = [share(t, same_dir_mask[t]) for t in CLASSES]
    nocall = [share(t, ["NoCall"]) for t in CLASSES]
    opp = [share(t, opp_mask[t]) for t in CLASSES]
    other = [100 - s - n - o for s, n, o in zip(same, nocall, opp)]

    x = np.arange(len(CLASSES))
    ax.bar(x, same, color=C["delta_up"], label="Same direction", alpha=0.9)
    ax.bar(x, other, bottom=same, color=C["VUS"], alpha=0.7,
           label="VUS / off-direction non-opposite")
    ax.bar(x, nocall, bottom=[s + o for s, o in zip(same, other)],
           color="#fbbf24", alpha=0.85, label="NoCall")
    ax.bar(x, opp,
           bottom=[s + o + n for s, o, n in zip(same, other, nocall)],
           color=C["delta_down"], alpha=0.9, label="Opposite direction")

    ax.set_xticks(x)
    ax.set_xticklabels([CLASS_SHORT[c] for c in CLASSES], fontweight="bold", fontsize=13)
    ax.set_xlabel("ClinVar 2-star+ truth class", fontweight="bold")
    ax.set_ylabel("% of class")
    ax.set_ylim(0, 110)
    ax.set_title("Per-class outcome breakdown", fontweight="bold")
    ax.legend(loc="upper right")
    ax.grid(axis="y", alpha=0.15)

    for i, (s, n) in enumerate(zip(same, [n_per[t] for t in CLASSES])):
        ax.text(i, s + 1.5, f"{s:.0f}%\nn={n:,}", ha="center", fontsize=10,
                fontweight="bold")

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(fig_dir / f"fig_outcome_breakdown.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  fig_outcome_breakdown")


def fig_criterion_fire_heatmap(fires, totals, fig_dir: Path):
    """Criterion × truth-class fire rate heatmap. Rows ordered by total
    fires; cells show % of variants in the truth class where the
    criterion fired. The dominant diagnostic of "which criteria are
    informative for which class".
    """
    # Order by total fire count, drop criteria that never fire
    ordered = sorted(
        [(c, fc, sum(fc.values())) for c, fc in fires.items() if sum(fc.values()) > 0],
        key=lambda kv: -kv[2],
    )
    codes = [c for c, _, _ in ordered]
    mat = np.zeros((len(codes), len(CLASSES)))
    counts = np.zeros((len(codes), len(CLASSES)), dtype=int)
    for i, (_, fc, _) in enumerate(ordered):
        for j, tcl in enumerate(CLASSES):
            n = totals[tcl]
            counts[i, j] = fc.get(tcl, 0)
            mat[i, j] = 100 * counts[i, j] / n if n else 0

    fig, ax = plt.subplots(figsize=(8.5, max(6, 0.32 * len(codes))))
    im = ax.imshow(mat, cmap="YlOrRd", aspect="auto", vmin=0, vmax=mat.max() or 1)
    ax.set_xticks(range(len(CLASSES)))
    ax.set_xticklabels([CLASS_SHORT[c] for c in CLASSES], fontweight="bold", fontsize=12)
    ax.set_yticks(range(len(codes)))
    ax.set_yticklabels(codes, fontfamily="monospace", fontsize=10)
    ax.set_xlabel("ClinVar 2-star+ truth class", fontweight="bold")
    ax.set_title("Criterion fire rate by truth class\n(% of class where criterion was met)",
                 fontweight="bold")

    for i in range(len(codes)):
        for j in range(len(CLASSES)):
            v = mat[i, j]
            if v > 0:
                color = "white" if v > 35 else "black"
                if v >= 1:
                    label = f"{v:.0f}%"
                elif v >= 0.1:
                    label = f"{v:.1f}%"
                else:
                    label = f"{counts[i, j]}"
                ax.text(j, i, label, ha="center", va="center", fontsize=8.5, color=color)

    cbar = plt.colorbar(im, ax=ax, shrink=0.7)
    cbar.set_label("% of class")
    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(fig_dir / f"fig_criterion_fire_heatmap.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  fig_criterion_fire_heatmap")


def fig_criterion_signatures_by_class(disc, fires, totals, fig_dir: Path):
    """For each non-VUS truth class, the top criterion signatures of
    *failed* calls (where fastVEP did not match truth direction). Bars
    are split by the failure mode: → VUS (insufficient evidence) vs
    opposite-direction (real disagreement).

    discrepancies.tsv only logs P/LP→VUS, P/LP→LB/B (opposite), and
    LB/B→P/LP (opposite). Same-direction calls don't appear here.
    """
    by_class: dict[str, Counter] = {tcl: Counter() for tcl in CLASSES}
    by_class_outcome: dict[str, dict[str, int]] = {
        tcl: defaultdict(lambda: defaultdict(int)) for tcl in CLASSES
    }
    for r in disc:
        sig = "+".join(sorted(set(c for c in r["met_criteria"].split(";") if c))) or "(none)"
        by_class[r["truth"]][sig] += 1
        by_class_outcome[r["truth"]][sig][r["predicted"]] += 1

    focus = ["Pathogenic", "Likely_pathogenic", "Likely_benign", "Benign"]
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    axes_flat = axes.flatten()

    for ax_i, tcl in enumerate(focus):
        ax = axes_flat[ax_i]
        sigs = by_class[tcl].most_common(8)
        if not sigs:
            ax.set_visible(False)
            continue

        opp_set = {
            "Pathogenic": {"Benign", "Likely_benign"},
            "Likely_pathogenic": {"Benign", "Likely_benign"},
            "Likely_benign": {"Pathogenic", "Likely_pathogenic"},
            "Benign": {"Pathogenic", "Likely_pathogenic"},
        }[tcl]

        labels: list[str] = []
        vus_counts: list[int] = []
        opp_counts: list[int] = []
        nocall_counts: list[int] = []
        for sig, _total in sigs:
            outcomes = by_class_outcome[tcl][sig]
            labels.append(sig if len(sig) <= 36 else sig[:34] + "…")
            vus_counts.append(outcomes.get("VUS", 0))
            opp_counts.append(sum(v for k, v in outcomes.items() if k in opp_set))
            nocall_counts.append(outcomes.get("NoCall", 0))

        y = np.arange(len(labels))
        ax.barh(y, vus_counts, color=C["VUS"], alpha=0.8, label="→ VUS (insufficient evidence)")
        left = np.array(vus_counts)
        ax.barh(y, opp_counts, left=left, color=C["delta_down"], alpha=0.9,
                label="Opposite direction")
        left = left + np.array(opp_counts)
        ax.barh(y, nocall_counts, left=left, color="#fbbf24", alpha=0.85,
                label="NoCall")

        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontfamily="monospace", fontsize=9.5)
        ax.invert_yaxis()
        n_truth = totals[tcl]
        n_failed = sum(by_class[tcl].values())
        ax.set_title(
            f"{CLASS_SHORT[tcl]} truth (n={n_truth:,};  failed={n_failed:,})\n"
            f"top 8 criterion signatures of failed calls",
            fontweight="bold", fontsize=13,
        )
        ax.set_xlabel("# variants in failure bucket", fontsize=11)
        ax.grid(axis="x", alpha=0.15)
        if ax_i == 0:
            ax.legend(loc="lower right", fontsize=10)

        for yi in range(len(labels)):
            tot = vus_counts[yi] + opp_counts[yi] + nocall_counts[yi]
            ax.text(tot + tot * 0.015, yi, f"{tot:,}", va="center",
                    fontsize=9, color="#374151")

    fig.suptitle(
        "Why fastVEP fails per truth class: criterion signatures of non-matching calls\n"
        "(Same-direction calls excluded — they aren't in discrepancies.tsv)",
        fontweight="bold", fontsize=14, y=1.01,
    )
    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(fig_dir / f"fig_criterion_signatures_by_class.{ext}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  fig_criterion_signatures_by_class")


def fig_lost_to_vus_signatures(disc, fig_dir: Path):
    """Pathogenic-truth → VUS and Likely_pathogenic-truth → VUS — the
    dominant failure mode. Show the top criterion signatures of these
    "lost" cases. Driven by missing manual-curation evidence (PS3
    functional, PP1 segregation, PP4 phenotype-specific) rather than
    classifier disagreement.
    """
    p_to_vus = Counter()
    lp_to_vus = Counter()
    for r in disc:
        if r["predicted"] != "VUS":
            continue
        sig = "+".join(sorted(set(c for c in r["met_criteria"].split(";") if c))) or "(none)"
        if r["truth"] == "Pathogenic":
            p_to_vus[sig] += 1
        elif r["truth"] == "Likely_pathogenic":
            lp_to_vus[sig] += 1

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7.5))
    for ax, ctr, label, title in (
        (ax1, p_to_vus, "P→VUS", f"Pathogenic-truth → VUS  (n={sum(p_to_vus.values()):,})"),
        (ax2, lp_to_vus, "LP→VUS", f"Likely-Pathogenic-truth → VUS  (n={sum(lp_to_vus.values()):,})"),
    ):
        sigs = ctr.most_common(15)
        labels = [s if len(s) <= 35 else s[:32] + "…" for s, _ in sigs]
        counts = [n for _, n in sigs]
        y = np.arange(len(labels))
        ax.barh(y, counts, color=C["P"] if "P→" == label[:2] else C["LP"], alpha=0.85)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontfamily="monospace", fontsize=10)
        ax.invert_yaxis()
        ax.set_xlabel("# variants", fontsize=11)
        ax.set_title(title, fontweight="bold", fontsize=13)
        ax.grid(axis="x", alpha=0.15)
        for yi, n in enumerate(counts):
            ax.text(n + max(counts) * 0.01, yi, f"{n:,}", va="center",
                    fontsize=10, color="#374151")

    fig.suptitle(
        "What criteria fired for the variants that fastVEP could only call VUS?\n"
        "(The dominant failure mode — typically driven by missing PS3 / PP1 / PP4 evidence)",
        fontweight="bold", fontsize=14, y=1.0,
    )
    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(fig_dir / f"fig_lost_to_vus_signatures.{ext}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  fig_lost_to_vus_signatures")


def fig_opposite_direction_signatures(disc, fig_dir: Path):
    """For the ~100 truly opposite-direction cases, show the criterion
    patterns by reversal direction.
    """
    PtoB = Counter()  # Pathogenic-truth → Likely_benign or Benign predicted
    BtoP = Counter()  # Benign-truth → Likely_pathogenic or Pathogenic predicted
    for r in disc:
        truth, pred = r["truth"], r["predicted"]
        sig = "+".join(sorted(set(c for c in r["met_criteria"].split(";") if c))) or "(none)"
        if truth in ("Pathogenic", "Likely_pathogenic") and pred in ("Likely_benign", "Benign"):
            PtoB[sig] += 1
        elif truth in ("Likely_benign", "Benign") and pred in ("Likely_pathogenic", "Pathogenic"):
            BtoP[sig] += 1

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    for ax, ctr, label, color, title in (
        (ax1, PtoB, "P/LP→B/LB", C["delta_down"],
         f"Pathogenic-tier truth → Benign-tier predicted  (n={sum(PtoB.values())})"),
        (ax2, BtoP, "B/LB→P/LP", "#a855f7",
         f"Benign-tier truth → Pathogenic-tier predicted  (n={sum(BtoP.values())})"),
    ):
        sigs = ctr.most_common(15)
        if not sigs:
            ax.text(0.5, 0.5, "(no cases)", ha="center", va="center",
                    transform=ax.transAxes, fontsize=14, color="#666")
            ax.set_title(title, fontweight="bold", fontsize=13)
            continue
        labels = [s if len(s) <= 35 else s[:32] + "…" for s, _ in sigs]
        counts = [n for _, n in sigs]
        y = np.arange(len(labels))
        ax.barh(y, counts, color=color, alpha=0.85)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontfamily="monospace", fontsize=10)
        ax.invert_yaxis()
        ax.set_xlabel("# variants", fontsize=11)
        ax.set_title(title, fontweight="bold", fontsize=13)
        ax.grid(axis="x", alpha=0.15)
        for yi, n in enumerate(counts):
            ax.text(n + max(counts) * 0.02, yi, f"{n}", va="center",
                    fontsize=10, color="#374151")

    fig.suptitle(
        "Opposite-direction discrepancies: criterion patterns by direction\n"
        "(These are the candidate set for medical-geneticist review)",
        fontweight="bold", fontsize=14, y=1.02,
    )
    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(fig_dir / f"fig_opposite_direction_signatures.{ext}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  fig_opposite_direction_signatures")


# ──────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) > 1:
        out_dir = Path(sys.argv[1])
    else:
        out_dir = Path(__file__).resolve().parents[3] / "data/benchmark/output_v7"
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading {out_dir}")
    matrix = read_matrix(out_dir)
    fires = read_criterion_fires(out_dir)
    disc = read_discrepancies(out_dir)
    totals = class_totals(matrix)

    # Wipe stale comparison-style figures from prior versions to keep
    # the figures/ directory clean.
    for stale in (
        "fig_v1_vs_v6_recall.png", "fig_v1_vs_v6_recall.pdf",
        "fig_v1_vs_v7_recall.png", "fig_v1_vs_v7_recall.pdf",
        "fig_headline_v1_vs_v6.png", "fig_headline_v1_vs_v6.pdf",
        "fig_headline_v1_vs_v7.png", "fig_headline_v1_vs_v7.pdf",
        "fig_bp7_pvs1_delta.png", "fig_bp7_pvs1_delta.pdf",
        "fig_recall_by_class.png", "fig_recall_by_class.pdf",
        "fig_criterion_fires.png", "fig_criterion_fires.pdf",
    ):
        try:
            (fig_dir / stale).unlink()
        except FileNotFoundError:
            pass

    print("Generating figures...")
    fig_concordance_matrix(matrix, fig_dir)
    fig_outcome_breakdown(matrix, fig_dir)
    fig_criterion_fire_heatmap(fires, totals, fig_dir)
    fig_criterion_signatures_by_class(disc, fires, totals, fig_dir)
    fig_lost_to_vus_signatures(disc, fig_dir)
    fig_opposite_direction_signatures(disc, fig_dir)

    print(f"\nDone. {len(list(fig_dir.glob('*.png')))} PNG / {len(list(fig_dir.glob('*.pdf')))} PDF in {fig_dir}")


if __name__ == "__main__":
    main()
