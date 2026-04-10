#!/usr/bin/env python3
"""
Generate manuscript figures for fastVEP.
Requires: matplotlib (pip install matplotlib)
"""

import csv
import os

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import matplotlib.patches as mpatches
    from matplotlib.patches import FancyBboxPatch
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("matplotlib not installed. Install with: pip install matplotlib")
    print("Generating text-based summaries instead.\n")

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')
FIG_DIR = os.path.dirname(__file__)

COLORS = {
    'primary': '#6c7aee',
    'vep': '#f5426c',
    'high': '#f5426c',
    'moderate': '#f59e0b',
    'low': '#3b82f6',
    'modifier': '#6b7280',
    'success': '#10b981',
    'cli': '#6c7aee',
    'mid': '#818cf8',
    'data': '#34d399',
    'core': '#f59e0b',
    'sa': '#f472b6',
}

ORG_COLORS = {
    'yeast': '#10b981',
    'drosophila': '#f59e0b',
    'arabidopsis': '#8b5cf6',
    'mouse': '#ef4444',
    'human': '#3b82f6',
}

ORG_LABELS = {
    'yeast': 'Yeast (R64)',
    'drosophila': 'Drosophila (BDGP6)',
    'arabidopsis': 'Arabidopsis (TAIR10)',
    'mouse': 'Mouse (GRCm39)',
    'human': 'Human (GRCh38)',
}

def read_csv(filename):
    path = os.path.join(DATA_DIR, filename)
    with open(path) as f:
        reader = csv.DictReader(f)
        return list(reader)


# ═══════════════════════════════════════════════════════════════════
# Fig 1: Architecture diagram
# ═══════════════════════════════════════════════════════════════════
def fig1_architecture():
    if not HAS_MPL:
        print("Figure 1: Architecture (text-only)")
        return

    fig, ax = plt.subplots(figsize=(14, 9))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 9)
    ax.axis('off')

    def box(x, y, w, h, label, sublabel, color, fontsize=11):
        rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15",
                              facecolor=color, edgecolor='#374151', linewidth=1.5, alpha=0.9)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2 + 0.15, label, ha='center', va='center',
                fontsize=fontsize, fontweight='bold', color='white')
        if sublabel:
            ax.text(x + w/2, y + h/2 - 0.2, sublabel, ha='center', va='center',
                    fontsize=8, color='white', alpha=0.9, style='italic')

    def arrow(x1, y1, x2, y2):
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', color='#6b7280', lw=1.5))

    ax.text(7, 8.6, 'fastVEP Architecture', ha='center', va='center',
            fontsize=16, fontweight='bold', color='#1f2937')
    ax.text(7, 8.25, '10-crate Cargo workspace', ha='center', va='center',
            fontsize=10, color='#6b7280')

    box(1.5, 7.0, 5, 0.9, 'fastvep-cli', 'Pipeline, cache builder, SA builder', COLORS['cli'], 13)
    box(7.5, 7.0, 5, 0.9, 'fastvep-web', 'Production web server (axum)', COLORS['cli'], 13)

    box(0.3, 5.3, 2.6, 0.9, 'fastvep-io', 'VCF/CSQ/JSON I/O', COLORS['mid'])
    box(3.2, 5.3, 2.6, 0.9, 'fastvep-hgvs', 'HGVSg/c/p', COLORS['mid'])
    box(6.1, 5.3, 2.8, 0.9, 'fastvep-consequence', 'SNV/indel/SV engine', COLORS['mid'])
    box(9.2, 5.3, 2.2, 0.9, 'fastvep-filter', 'Filter engine', COLORS['mid'])
    box(11.7, 5.3, 2.0, 0.9, 'fastvep-sa', 'Annotations', COLORS['sa'])

    box(0.5, 3.5, 3.5, 0.9, 'fastvep-genome', 'Transcript, Exon, Gene, CodonTable', COLORS['data'])
    box(4.5, 3.5, 5.0, 0.9, 'fastvep-cache', 'GFF3, FASTA mmap, tabix, binary cache', COLORS['data'])
    box(10.0, 3.5, 3.5, 0.9, 'fastvep-sa (fastSA)', 'ClinVar, gnomAD, REVEL, ...', COLORS['sa'])

    box(3.0, 1.7, 8.0, 0.9, 'fastvep-core', 'GenomicPosition, Consequence (49 SO terms), Allele, Strand, Impact, VariantType', COLORS['core'], 12)

    for x in [1.6, 4.5, 7.5, 10.3, 12.7]:
        arrow(4, 7.0, x, 6.2)
        arrow(10, 7.0, x, 6.2)
    for x in [2.25, 7.0, 11.75]:
        arrow(x, 5.3, x, 4.4)
    for x in [2.25, 7.0, 11.75]:
        arrow(x, 3.5, 7.0, 2.6)

    stats_text = "17,966 LOC  |  173 tests  |  3.2 MB binary (LTO + strip)  |  fastVEP.org"
    ax.text(7, 0.8, stats_text, ha='center', va='center',
            fontsize=10, color='#6b7280',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#f3f4f6', edgecolor='#d1d5db'))

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig1_architecture.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig1_architecture.pdf'), bbox_inches='tight')
    print("Saved fig1_architecture.png/pdf")
    plt.close()


# ═══════════════════════════════════════════════════════════════════
# Fig 2: Multi-organism benchmark (THE money shot)
# Scatter: variant count vs wall-clock time, bubble size = transcripts
# ═══════════════════════════════════════════════════════════════════
def fig2_multi_organism_benchmark():
    data = read_csv('organism_comparison.csv')

    if not HAS_MPL:
        print("Figure 2: Multi-organism benchmark")
        for d in data:
            print(f"  {d['organism']:15s} {int(d['variants']):>12,} variants  {float(d['time_sec']):>6.1f}s  {int(d['variants_per_sec']):>8,} v/s")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    orgs = [d['organism'] for d in data]
    variants = [int(d['variants']) for d in data]
    times = [float(d['time_sec']) for d in data]
    vps = [int(d['variants_per_sec']) for d in data]
    transcripts = [int(d['transcripts']) for d in data]
    colors = [ORG_COLORS.get(o, '#999') for o in orgs]
    sizes = [max(80, t/1000 * 15) for t in transcripts]  # scale bubble by transcript count

    # Panel A: Variants vs Wall-clock time (log-log scatter)
    for o, v, t, c, s, tr in zip(orgs, variants, times, colors, sizes, transcripts):
        ax1.scatter(v, t, c=c, s=s, alpha=0.85, edgecolors='white', linewidth=1.5, zorder=3)
        label = ORG_LABELS.get(o, o)
        offset_y = 0.15 if o != 'mouse' else -0.2
        ax1.annotate(f'{label}\n{v/1e6:.1f}M in {t:.0f}s' if v > 500000 else f'{label}\n{v/1e3:.0f}K in {t:.1f}s',
                     xy=(v, t), xytext=(12, 8),
                     textcoords='offset points', fontsize=8.5, color='#333',
                     fontweight='bold')

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('Variants annotated (gold-standard datasets)', fontsize=12)
    ax1.set_ylabel('Wall-clock time (seconds)', fontsize=12)
    ax1.set_title('A. Annotation Time Across Organisms', fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.2)
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1e6:.0f}M' if x >= 1e6 else f'{x/1e3:.0f}K'))

    # Panel B: Throughput by organism (horizontal bar)
    sorted_data = sorted(zip(orgs, vps, transcripts, colors), key=lambda x: x[1])
    s_orgs, s_vps, s_tr, s_colors = zip(*sorted_data)
    labels = [f"{ORG_LABELS.get(o,o)}\n({tr:,} transcripts)" for o, tr in zip(s_orgs, s_tr)]

    bars = ax2.barh(range(len(s_orgs)), s_vps, color=s_colors, alpha=0.85, height=0.6)
    ax2.set_yticks(range(len(s_orgs)))
    ax2.set_yticklabels(labels, fontsize=10)
    ax2.set_xlabel('Throughput (variants/sec)', fontsize=12)
    ax2.set_title('B. Peak Throughput by Organism', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.2, axis='x')
    ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1000:.0f}K'))

    for bar, v in zip(bars, s_vps):
        ax2.text(bar.get_width() + 2000, bar.get_y() + bar.get_height()/2,
                 f'{v:,} v/s', va='center', fontsize=10, fontweight='bold', color='#333')

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig2_multi_organism.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig2_multi_organism.pdf'), bbox_inches='tight')
    print("Saved fig2_multi_organism.png/pdf")
    plt.close()


# ═══════════════════════════════════════════════════════════════════
# Fig 3: fastVEP vs Ensembl VEP head-to-head (NEW — the striking comparison)
# ═══════════════════════════════════════════════════════════════════
def fig3_vep_comparison():
    data = read_csv('vep_comparison.csv')
    variants = [int(d['variants']) for d in data]
    f_time = [float(d['fastvep_sec']) for d in data]
    v_time = [float(d['vep_sec']) for d in data]
    f_vps = [int(d['fastvep_vps']) for d in data]
    v_vps = [int(d['vep_vps']) for d in data]
    speedups = [float(d['speedup']) for d in data]

    if not HAS_MPL:
        print("Figure 3: VEP Comparison")
        for v, ft, vt, s in zip(variants, f_time, v_time, speedups):
            print(f"  {v:>6,}: fastVEP {ft:.2f}s  VEP {vt:.1f}s  {s:.0f}x")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: Wall-clock time comparison (log-log)
    ax1.plot(variants, f_time, 'o-', color=COLORS['primary'], linewidth=2.5, markersize=10, label='fastVEP', zorder=3)
    ax1.plot(variants, v_time, 's-', color=COLORS['vep'], linewidth=2.5, markersize=10, label='Ensembl VEP v115', zorder=3)

    # Add full-genome extrapolation for VEP (dashed)
    ax1.plot([50000, 4048342], [96.79, 96.79 * (4048342/50000)], 's--', color=COLORS['vep'],
             linewidth=1.5, markersize=8, alpha=0.4, zorder=2)
    ax1.scatter([4048342], [71], c=COLORS['primary'], s=150, marker='*', zorder=4, edgecolors='white', linewidth=2)
    ax1.annotate('fastVEP full WGS\n4.05M variants in 71s', xy=(4048342, 71), xytext=(-80, 30),
                 textcoords='offset points', fontsize=9, fontweight='bold', color=COLORS['primary'],
                 arrowprops=dict(arrowstyle='->', color=COLORS['primary'], lw=1.5))
    ax1.annotate('VEP: cannot complete\n(est. ~2.2 hours)', xy=(4048342, 96.79*(4048342/50000)), xytext=(-80, -30),
                 textcoords='offset points', fontsize=9, color=COLORS['vep'], style='italic',
                 arrowprops=dict(arrowstyle='->', color=COLORS['vep'], lw=1, ls='--'))

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('Number of variants (GIAB HG002 chr22)', fontsize=12)
    ax1.set_ylabel('Wall-clock time (seconds)', fontsize=12)
    ax1.set_title('A. Annotation Time', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11, loc='upper left')
    ax1.grid(True, alpha=0.2)
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1e6:.1f}M' if x >= 1e6 else f'{x/1e3:.0f}K'))

    # Panel B: Throughput comparison
    ax2.plot(variants, f_vps, 'o-', color=COLORS['primary'], linewidth=2.5, markersize=10, label='fastVEP', zorder=3)
    ax2.plot(variants, v_vps, 's-', color=COLORS['vep'], linewidth=2.5, markersize=10, label='Ensembl VEP v115', zorder=3)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('Number of variants', fontsize=12)
    ax2.set_ylabel('Throughput (variants/sec)', fontsize=12)
    ax2.set_title('B. Annotation Throughput', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.2)
    ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1e6:.1f}M' if x >= 1e6 else f'{x/1e3:.0f}K'))
    ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1000:.0f}K' if x >= 1000 else f'{x:.0f}'))

    for v, fv, vv, s in zip(variants, f_vps, v_vps, speedups):
        if s > 1:
            ax2.annotate(f'{s:.0f}x', xy=(v, fv), xytext=(0, 14),
                         textcoords='offset points', ha='center', fontsize=10,
                         color=COLORS['primary'], fontweight='bold')

    # Highlight the divergence
    ax2.fill_between(variants, v_vps, f_vps, alpha=0.1, color=COLORS['primary'])

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig3_vep_comparison.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig3_vep_comparison.pdf'), bbox_inches='tight')
    print("Saved fig3_vep_comparison.png/pdf")
    plt.close()


# ═══════════════════════════════════════════════════════════════════
# Fig 4: VEP concordance
# ═══════════════════════════════════════════════════════════════════
def fig4_vep_concordance():
    data = read_csv('vep_concordance.csv')
    fields = [d['field'] for d in data]
    accuracy = [float(d['accuracy']) for d in data]

    if not HAS_MPL:
        print("Figure 3: VEP Concordance")
        for f, a in zip(fields, accuracy):
            print(f"  {f:25s} {a:.1f}%")
        return

    fig, ax = plt.subplots(figsize=(10, 7))
    colors = [COLORS['success'] if a == 100.0 else COLORS['vep'] for a in accuracy]
    bars = ax.barh(range(len(fields)), accuracy, color=colors, alpha=0.85)
    ax.set_yticks(range(len(fields)))
    ax.set_yticklabels(fields, fontsize=10, fontfamily='monospace')
    ax.invert_yaxis()
    ax.set_xlabel('Concordance with Ensembl VEP v115.1 (%)', fontsize=12)
    ax.set_title('Field-Level Accuracy: fastVEP vs Ensembl VEP\n(2,340 shared transcript-allele pairs, 173 variants)', fontsize=13, fontweight='bold')
    ax.set_xlim(95, 101)
    ax.grid(True, alpha=0.3, axis='x')
    ax.axvline(x=100, color=COLORS['success'], linestyle='--', alpha=0.5, linewidth=1)

    for bar, acc in zip(bars, accuracy):
        ax.text(bar.get_width() - 0.3, bar.get_y() + bar.get_height()/2,
                f'{acc:.1f}%', va='center', ha='right', fontsize=9, color='white', fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig4_vep_concordance.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig4_vep_concordance.pdf'), bbox_inches='tight')
    print("Saved fig4_vep_concordance.png/pdf")
    plt.close()


# ═══════════════════════════════════════════════════════════════════
# Fig 4: Consequence distribution (keep as-is)
# ═══════════════════════════════════════════════════════════════════
def fig5_consequence_distribution():
    data = read_csv('consequence_distribution.csv')
    consequences = [d['consequence'] for d in data]
    counts = [int(d['count']) for d in data]

    impact_map = {
        'splice_acceptor_variant': 'HIGH', 'splice_donor_variant': 'HIGH',
        'stop_gained': 'HIGH', 'frameshift_variant': 'HIGH',
        'missense_variant': 'MODERATE', 'inframe_insertion': 'MODERATE',
        'inframe_deletion': 'MODERATE',
        'splice_region_variant': 'LOW', 'synonymous_variant': 'LOW',
        'splice_polypyrimidine_tract_variant': 'LOW',
        'splice_donor_5th_base_variant': 'LOW', 'splice_donor_region_variant': 'LOW',
    }
    def get_color(c):
        imp = impact_map.get(c, 'MODIFIER')
        return COLORS.get(imp.lower(), COLORS['modifier'])

    if not HAS_MPL:
        print("Figure 4: Consequence Distribution")
        for c, n in zip(consequences, counts):
            print(f"  {c:45s} {n:>6d}")
        return

    fig, ax = plt.subplots(figsize=(11, 7))
    colors = [get_color(c) for c in consequences]
    bars = ax.barh(range(len(consequences)), counts, color=colors, alpha=0.85)
    ax.set_yticks(range(len(consequences)))
    ax.set_yticklabels(consequences, fontsize=10, fontfamily='monospace')
    ax.invert_yaxis()
    ax.set_xlabel('Annotation count', fontsize=12)
    ax.set_title('Predicted Consequence Distribution\n(1,000 real 1KGP chr22 variants, 28,161 annotations)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    ax.set_xscale('log')

    for bar, count in zip(bars, counts):
        ax.text(bar.get_width() * 1.1, bar.get_y() + bar.get_height()/2,
                f'{count:,}', va='center', fontsize=9, color='#555')

    legend_elements = [
        mpatches.Patch(facecolor=COLORS['high'], alpha=0.85, label='HIGH'),
        mpatches.Patch(facecolor=COLORS['moderate'], alpha=0.85, label='MODERATE'),
        mpatches.Patch(facecolor=COLORS['low'], alpha=0.85, label='LOW'),
        mpatches.Patch(facecolor=COLORS['modifier'], alpha=0.85, label='MODIFIER'),
    ]
    ax.legend(handles=legend_elements, title='Impact', loc='lower right', fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig5_consequence_distribution.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig5_consequence_distribution.pdf'), bbox_inches='tight')
    print("Saved fig5_consequence_distribution.png/pdf")
    plt.close()


# ═══════════════════════════════════════════════════════════════════
# Fig 5: Throughput vs Genome Complexity (NEW — replaces old resource usage)
# Shows throughput stays high even as transcript count increases 72x
# ═══════════════════════════════════════════════════════════════════
def fig6_throughput_vs_complexity():
    data = read_csv('organism_comparison.csv')

    if not HAS_MPL:
        print("Figure 5: Throughput vs Complexity")
        for d in data:
            print(f"  {d['organism']:15s} {int(d['transcripts']):>8,} transcripts  {int(d['variants_per_sec']):>8,} v/s")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    orgs = [d['organism'] for d in data]
    transcripts = [int(d['transcripts']) for d in data]
    vps = [int(d['variants_per_sec']) for d in data]
    times = [float(d['time_sec']) for d in data]
    variants = [int(d['variants']) for d in data]
    colors = [ORG_COLORS.get(o, '#999') for o in orgs]

    # Panel A: Throughput vs Transcript count (shows consistency)
    for o, tr, v, c in zip(orgs, transcripts, vps, colors):
        ax1.scatter(tr, v, c=c, s=200, alpha=0.85, edgecolors='white', linewidth=2, zorder=3)
        ax1.annotate(ORG_LABELS.get(o, o), xy=(tr, v), xytext=(8, 8),
                     textcoords='offset points', fontsize=9, fontweight='bold', color='#333')

    ax1.set_xscale('log')
    ax1.set_xlabel('Gene model complexity (number of transcripts)', fontsize=12)
    ax1.set_ylabel('Throughput (variants/sec)', fontsize=12)
    ax1.set_title('A. Throughput vs Genome Complexity', fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.2)
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1000:.0f}K'))
    ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1000:.0f}K'))

    # Add shaded band showing the 57K-177K throughput range
    ax1.axhspan(50000, 185000, alpha=0.08, color=COLORS['primary'])
    ax1.text(300000, 120000, '57K–177K v/s\nacross 72x complexity range',
             ha='center', fontsize=9, color=COLORS['primary'], style='italic')

    # Panel B: Total variants processed (shows scale)
    sorted_data = sorted(zip(orgs, variants, times, colors), key=lambda x: x[1])
    s_orgs, s_vars, s_times, s_colors = zip(*sorted_data)
    labels = [ORG_LABELS.get(o, o) for o in s_orgs]

    bars = ax2.barh(range(len(s_orgs)), [v/1e6 for v in s_vars], color=s_colors, alpha=0.85, height=0.6)
    ax2.set_yticks(range(len(s_orgs)))
    ax2.set_yticklabels(labels, fontsize=10)
    ax2.set_xlabel('Variants annotated (millions)', fontsize=12)
    ax2.set_title('B. Gold-Standard Dataset Scale', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.2, axis='x')

    for bar, v, t in zip(bars, s_vars, s_times):
        ax2.text(bar.get_width() + 0.15, bar.get_y() + bar.get_height()/2,
                 f'{v/1e6:.1f}M variants in {t:.0f}s',
                 va='center', fontsize=10, fontweight='bold', color='#333')

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig6_throughput_complexity.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig6_throughput_complexity.pdf'), bbox_inches='tight')
    print("Saved fig6_throughput_complexity.png/pdf")
    plt.close()


if __name__ == '__main__':
    print("Generating fastVEP manuscript figures...")
    print(f"Data dir: {os.path.abspath(DATA_DIR)}")
    print(f"Output dir: {os.path.abspath(FIG_DIR)}")
    print()
    fig1_architecture()
    fig2_multi_organism_benchmark()
    fig3_vep_comparison()
    fig4_vep_concordance()
    fig5_consequence_distribution()
    fig6_throughput_vs_complexity()
    print("\nDone.")
