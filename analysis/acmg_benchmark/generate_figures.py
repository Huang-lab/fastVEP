#!/usr/bin/env python3
"""
Generate ACMG benchmark figures from concordance analysis output.

Reads CSV files from output/ and generates publication-quality figures.

Usage:
  python generate_figures.py          # uses output/ in same directory
  python generate_figures.py <dir>    # uses specified output directory
"""

import csv
import os
import sys

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import matplotlib.patches as mpatches
    import numpy as np
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("matplotlib not installed. Install with: pip install matplotlib")
    sys.exit(1)

# Global font size bump
plt.rcParams.update({
    'font.size': 14,
    'axes.titlesize': 18,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 13,
    'legend.title_fontsize': 14,
})

# Colors
C = {
    'P': '#ef4444', 'LP': '#f97316', 'VUS': '#6b7280',
    'LB': '#3b82f6', 'B': '#10b981',
    'concordant': '#10b981', 'discordant': '#ef4444', 'downgraded': '#f59e0b',
    'bar': '#6c7aee', 'highlight': '#f5426c',
}
CLASS_COLORS = {
    'Pathogenic': C['P'], 'Likely_pathogenic': C['LP'],
    'VUS': C['VUS'], 'Likely_benign': C['LB'], 'Benign': C['B'],
}
CLASS_SHORT = {
    'Pathogenic': 'P', 'Likely_pathogenic': 'LP',
    'VUS': 'VUS', 'Likely_benign': 'LB', 'Benign': 'B',
}
CLASSES = ['Pathogenic', 'Likely_pathogenic', 'VUS', 'Likely_benign', 'Benign']


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def fig_concordance_matrix(data_dir, fig_dir):
    """Heatmap of concordance matrix."""
    rows = read_csv(os.path.join(data_dir, "concordance_matrix.csv"))

    # Build matrix
    matrix = {}
    for r in rows:
        matrix[(r['ClinVar_class'], r['Predicted_class'])] = int(r['count'])

    fig, ax = plt.subplots(figsize=(10, 8))

    # Create heatmap data
    mat = np.zeros((5, 5))
    for i, cv in enumerate(CLASSES):
        row_total = sum(matrix.get((cv, p), 0) for p in CLASSES)
        for j, pred in enumerate(CLASSES):
            mat[i, j] = 100 * matrix.get((cv, pred), 0) / row_total if row_total > 0 else 0

    im = ax.imshow(mat, cmap='YlOrRd', aspect='auto', vmin=0, vmax=100)

    ax.set_xticks(range(5))
    ax.set_xticklabels([CLASS_SHORT[c] for c in CLASSES], fontsize=16, fontweight='bold')
    ax.set_yticks(range(5))
    ax.set_yticklabels([CLASS_SHORT[c] for c in CLASSES], fontsize=16, fontweight='bold')
    ax.set_xlabel('fastVEP Predicted', fontsize=18, fontweight='bold')
    ax.set_ylabel('ClinVar Gold Standard', fontsize=18, fontweight='bold')
    ax.set_title('ACMG Classification Concordance\n(ClinVar 2-star+ GRCh38 variants)',
                 fontsize=20, fontweight='bold')

    # Annotate cells
    for i in range(5):
        row_total = sum(matrix.get((CLASSES[i], p), 0) for p in CLASSES)
        for j in range(5):
            count = matrix.get((CLASSES[i], CLASSES[j]), 0)
            pct = mat[i, j]
            color = 'white' if pct > 50 else 'black'
            ax.text(j, i, f'{count:,}\n({pct:.0f}%)', ha='center', va='center',
                    fontsize=13, color=color, fontweight='bold' if i == j else 'normal')

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Row percentage (%)', fontsize=15)
    cbar.ax.tick_params(labelsize=13)

    # Highlight diagonal
    for i in range(5):
        ax.add_patch(plt.Rectangle((i - 0.5, i - 0.5), 1, 1,
                                   fill=False, edgecolor='#10b981', lw=3))

    plt.tight_layout()
    for ext in ['png', 'pdf']:
        plt.savefig(os.path.join(fig_dir, f'fig_concordance_matrix.{ext}'),
                    dpi=300, bbox_inches='tight')
    print("  Saved fig_concordance_matrix")
    plt.close()


def fig_summary_metrics(data_dir, fig_dir):
    """Bar chart of key concordance metrics."""
    rows = read_csv(os.path.join(data_dir, "summary_metrics.csv"))
    metrics = {r['metric']: r for r in rows}

    fig, ax = plt.subplots(figsize=(11, 6))

    labels = ['Exact\nmatch', 'Same\ndirection', 'Opposite\ndirection',
              'Pathogenic\nsensitivity', 'Benign\nsensitivity']
    keys = ['exact_match', 'same_direction', 'opposite_direction',
            'pathogenic_sensitivity', 'benign_sensitivity']
    colors = [C['concordant'], C['concordant'], C['discordant'],
              C['P'], C['B']]

    values = []
    for k in keys:
        if k in metrics and metrics[k]['percentage'] != 'N/A':
            values.append(float(metrics[k]['percentage']))
        else:
            values.append(0)

    bars = ax.bar(range(len(labels)), values, color=colors, alpha=0.85, width=0.6)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=15)
    ax.set_ylabel('Percentage (%)', fontsize=17)
    ax.set_title('ACMG Classification Concordance Metrics\n(fastVEP vs ClinVar 2-star+ gold standard)',
                 fontsize=19, fontweight='bold')
    ax.set_ylim(0, 105)
    ax.grid(True, alpha=0.15, axis='y')
    ax.tick_params(axis='y', labelsize=14)

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.5,
                f'{val:.1f}%', ha='center', fontsize=16, fontweight='bold')

    total = int(metrics.get('total_2star_plus', {}).get('value', 0))
    if total:
        ax.text(0.98, 0.02, f'n = {total:,} variants',
                transform=ax.transAxes, ha='right', va='bottom',
                fontsize=14, color='#666', style='italic')

    plt.tight_layout()
    for ext in ['png', 'pdf']:
        plt.savefig(os.path.join(fig_dir, f'fig_summary_metrics.{ext}'),
                    dpi=300, bbox_inches='tight')
    print("  Saved fig_summary_metrics")
    plt.close()


def fig_consequence_concordance(data_dir, fig_dir):
    """Horizontal bar chart: concordance by consequence type."""
    rows = read_csv(os.path.join(data_dir, "concordance_by_consequence.csv"))

    # Sort by N descending
    rows = sorted(rows, key=lambda r: -int(r['n']))[:12]  # top 12

    fig, ax = plt.subplots(figsize=(14, 8))

    consequences = [r['consequence'] for r in rows]
    same_dir = [float(r['same_dir_pct']) for r in rows]
    n_vals = [int(r['n']) for r in rows]
    opposite = [float(r['opposite_pct']) for r in rows]

    y = range(len(consequences))
    colors = [C['concordant'] if sd > 60 else C['downgraded'] if sd > 40 else C['discordant']
              for sd in same_dir]

    bars = ax.barh(y, same_dir, color=colors, alpha=0.85, height=0.6)
    ax.barh(y, [-o for o in opposite], color=C['discordant'], alpha=0.3, height=0.6)

    ax.set_yticks(y)
    ax.set_yticklabels(consequences, fontsize=14, fontfamily='monospace')
    ax.invert_yaxis()
    ax.set_xlabel('Same-direction concordance (%)', fontsize=17)
    ax.set_title('Concordance by Consequence Type\n(fastVEP ACMG vs ClinVar 2-star+)',
                 fontsize=19, fontweight='bold')
    ax.axvline(x=0, color='#333', lw=0.5)
    ax.set_xlim(-10, 110)
    ax.grid(True, alpha=0.15, axis='x')
    ax.tick_params(axis='x', labelsize=14)

    for bar, sd, n in zip(bars, same_dir, n_vals):
        ax.text(bar.get_width() + 1.5, bar.get_y() + bar.get_height() / 2,
                f'{sd:.0f}%  (n={n:,})', va='center', fontsize=13, color='#333')

    legend_el = [
        mpatches.Patch(facecolor=C['concordant'], alpha=0.85, label='High concordance (>60%)'),
        mpatches.Patch(facecolor=C['downgraded'], alpha=0.85, label='Moderate (40-60%)'),
        mpatches.Patch(facecolor=C['discordant'], alpha=0.85, label='Low (<40%)'),
    ]
    ax.legend(handles=legend_el, loc='lower right', fontsize=13)

    plt.tight_layout()
    for ext in ['png', 'pdf']:
        plt.savefig(os.path.join(fig_dir, f'fig_consequence_concordance.{ext}'),
                    dpi=300, bbox_inches='tight')
    print("  Saved fig_consequence_concordance")
    plt.close()


def fig_variant_distribution(data_dir, fig_dir):
    """Stacked bar chart: variant type distribution by ClinVar class."""
    rows = read_csv(os.path.join(data_dir, "variant_counts.csv"))

    # Get top consequence types
    cons_totals = {}
    for r in rows:
        cons = r['consequence']
        cons_totals[cons] = cons_totals.get(cons, 0) + int(r['count'])
    top_cons = sorted(cons_totals.keys(), key=lambda c: -cons_totals[c])[:10]

    # Build data
    data = {c: {cons: 0 for cons in top_cons} for c in CLASSES}
    for r in rows:
        if r['consequence'] in top_cons and r['clinvar_class'] in CLASSES:
            data[r['clinvar_class']][r['consequence']] = int(r['count'])

    fig, ax = plt.subplots(figsize=(13, 7))

    x = np.arange(len(CLASSES))
    width = 0.65
    bottom = np.zeros(len(CLASSES))

    cons_colors = plt.cm.Set3(np.linspace(0, 1, len(top_cons)))

    for i, cons in enumerate(top_cons):
        values = [data[c][cons] for c in CLASSES]
        ax.bar(x, values, width, bottom=bottom, label=cons,
               color=cons_colors[i], edgecolor='white', lw=0.5)
        bottom += values

    ax.set_xticks(x)
    ax.set_xticklabels([CLASS_SHORT[c] for c in CLASSES], fontsize=16, fontweight='bold')
    ax.set_xlabel('ClinVar Classification', fontsize=17)
    ax.set_ylabel('Variant Count', fontsize=17)
    ax.set_title('Variant Consequence Distribution by ClinVar Class\n(2-star+ GRCh38)',
                 fontsize=19, fontweight='bold')
    ax.legend(title='Consequence', bbox_to_anchor=(1.02, 1), loc='upper left',
              fontsize=12, title_fontsize=13)
    ax.grid(True, alpha=0.15, axis='y')
    ax.tick_params(axis='y', labelsize=14)

    plt.tight_layout()
    for ext in ['png', 'pdf']:
        plt.savefig(os.path.join(fig_dir, f'fig_variant_distribution.{ext}'),
                    dpi=300, bbox_inches='tight')
    print("  Saved fig_variant_distribution")
    plt.close()


def fig_criteria_coverage(fig_dir):
    """Horizontal bar chart showing which criteria are automatable."""
    criteria = [
        ("PVS1", "Null variant in LOF gene", True, "Very Strong"),
        ("PS1", "Same AA change pathogenic", True, "Strong"),
        ("PS2", "De novo (confirmed)", False, "Strong"),
        ("PS3", "Functional studies", False, "Strong"),
        ("PS4", "Prevalence in affected", True, "Strong"),
        ("PM1", "Mutational hotspot", True, "Moderate"),
        ("PM2", "Absent from controls", True, "Moderate*"),
        ("PM3", "Compound het (recessive)", True, "Moderate"),
        ("PM4", "Protein length change", True, "Moderate"),
        ("PM5", "Novel missense at known pos", True, "Moderate"),
        ("PM6", "Assumed de novo", True, "Moderate"),
        ("PP1", "Co-segregation", False, "Supporting"),
        ("PP2", "Missense constrained gene", True, "Supporting"),
        ("PP3", "Computational evidence", True, "Supporting+"),
        ("PP4", "Phenotype specific", False, "Supporting"),
        ("PP5", "Reputable source", False, "Supporting"),
        ("BA1", "AF > 5%", True, "Standalone"),
        ("BS1", "AF > expected", True, "Strong"),
        ("BS2", "Healthy adult homozygous", True, "Strong"),
        ("BS3", "Functional (no effect)", False, "Strong"),
        ("BS4", "Lack of segregation", False, "Strong"),
        ("BP1", "Missense in LOF gene", True, "Supporting"),
        ("BP2", "In cis/trans analysis", True, "Supporting"),
        ("BP3", "In-frame in repeat", True, "Supporting"),
        ("BP4", "Computational benign", True, "Supporting+"),
        ("BP5", "Alt molecular basis", False, "Supporting"),
        ("BP6", "Reputable benign", False, "Supporting"),
        ("BP7", "Synonymous no splice", True, "Supporting"),
    ]

    fig, ax = plt.subplots(figsize=(14, 12))

    y_positions = range(len(criteria))
    colors = [C['concordant'] if auto else '#d1d5db' for _, _, auto, _ in criteria]

    bars = ax.barh(y_positions, [1 if auto else 0.5 for _, _, auto, _ in criteria],
                   color=colors, alpha=0.85, height=0.7)

    labels = [f"{code}  {desc}" for code, desc, _, _ in criteria]

    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=13, fontfamily='monospace')
    ax.invert_yaxis()
    ax.set_xlim(0, 1.45)
    ax.set_xticks([])
    ax.set_title('ACMG-AMP Criteria Coverage in fastVEP\n(28 criteria: 18 automatable, 10 require manual curation)',
                 fontsize=19, fontweight='bold')

    for bar, (code, desc, auto, strength) in zip(bars, criteria):
        label = "Automated" if auto else "Manual"
        color = '#065f46' if auto else '#666'
        ax.text(bar.get_width() + 0.02, bar.get_y() + bar.get_height() / 2,
                f'{label}  [{strength}]', va='center', fontsize=12, color=color)

    legend_el = [
        mpatches.Patch(facecolor=C['concordant'], alpha=0.85, label='Automated (18/28)'),
        mpatches.Patch(facecolor='#d1d5db', alpha=0.85, label='Requires manual curation (10/28)'),
    ]
    ax.legend(handles=legend_el, loc='lower right', fontsize=14)

    # PM2 footnote
    ax.text(0.02, 0.02, '* PM2 downgraded to Supporting per ClinGen SVI recommendation\n'
                         '+ PP3/BP4 use ClinGen SVI calibrated REVEL thresholds with strength escalation',
            transform=ax.transAxes, fontsize=11, color='#666', style='italic', va='bottom')

    plt.tight_layout()
    for ext in ['png', 'pdf']:
        plt.savefig(os.path.join(fig_dir, f'fig_criteria_coverage.{ext}'),
                    dpi=300, bbox_inches='tight')
    print("  Saved fig_criteria_coverage")
    plt.close()


def main():
    if len(sys.argv) > 1:
        data_dir = sys.argv[1]
    else:
        data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

    fig_dir = os.path.join(data_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    print("Generating ACMG benchmark figures...\n")

    # Criteria coverage figure (no data dependency)
    fig_criteria_coverage(fig_dir)

    # Data-dependent figures
    matrix_file = os.path.join(data_dir, "concordance_matrix.csv")
    if os.path.exists(matrix_file):
        fig_concordance_matrix(data_dir, fig_dir)
        fig_summary_metrics(data_dir, fig_dir)
        fig_consequence_concordance(data_dir, fig_dir)
        fig_variant_distribution(data_dir, fig_dir)
    else:
        print(f"\n  No concordance data in {data_dir}/")
        print("  Run clinvar_concordance.py first to generate data.")
        print("  Only criteria coverage figure was generated.\n")

    print("\nDone.")


if __name__ == '__main__':
    main()
