#!/usr/bin/env python3
"""
Generate manuscript figures for OxiVEP.
Requires: matplotlib (pip install matplotlib)
"""

import csv
import os

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
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
}

def read_csv(filename):
    path = os.path.join(DATA_DIR, filename)
    with open(path) as f:
        reader = csv.DictReader(f)
        return list(reader)


def fig2_throughput_comparison():
    """Figure 2: OxiVEP vs Ensembl VEP throughput comparison."""
    data = read_csv('scaling.csv')
    variants = [int(d['variants']) for d in data]
    oxi_vps = [float(d['oxivep_vps']) for d in data]
    vep_vps = [float(d['vep_vps']) for d in data]
    oxi_time = [float(d['oxivep_time_sec']) for d in data]
    vep_time = [float(d['vep_time_sec']) for d in data]

    if not HAS_MPL:
        print("Figure 2: Throughput Comparison")
        for v, ot, ovps, vt, vvps in zip(variants, oxi_time, oxi_vps, vep_time, vep_vps):
            print(f"  {v:>8,}v: OxiVEP {ot:.2f}s ({ovps:,.0f} v/s)  VEP {vt:.1f}s ({vvps:,.0f} v/s)  {ovps/vvps:.0f}x")
        print()
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    # Panel A: Wall-clock time
    ax1.plot(variants, oxi_time, 'o-', color=COLORS['primary'], linewidth=2.5, markersize=9, label='OxiVEP', zorder=3)
    ax1.plot(variants, vep_time, 's-', color=COLORS['vep'], linewidth=2.5, markersize=9, label='Ensembl VEP v115', zorder=3)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('Number of variants', fontsize=12)
    ax1.set_ylabel('Wall-clock time (seconds)', fontsize=12)
    ax1.set_title('A. Annotation Time', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1000:.0f}K' if x < 1e6 else f'{x/1e6:.0f}M'))

    # Panel B: Throughput
    ax2.plot(variants, oxi_vps, 'o-', color=COLORS['primary'], linewidth=2.5, markersize=9, label='OxiVEP', zorder=3)
    ax2.plot(variants, vep_vps, 's-', color=COLORS['vep'], linewidth=2.5, markersize=9, label='Ensembl VEP v115', zorder=3)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('Number of variants', fontsize=12)
    ax2.set_ylabel('Throughput (variants/sec)', fontsize=12)
    ax2.set_title('B. Annotation Throughput', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1000:.0f}K' if x < 1e6 else f'{x/1e6:.0f}M'))
    ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1000:.0f}K' if x < 1e6 else f'{x/1e6:.1f}M'))

    # Add speedup annotations
    for v, ovps, vvps in zip(variants, oxi_vps, vep_vps):
        speedup = ovps / vvps
        ax2.annotate(f'{speedup:.0f}x', xy=(v, ovps), xytext=(0, 12),
                     textcoords='offset points', ha='center', fontsize=9,
                     color=COLORS['primary'], fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig2_throughput_scaling.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig2_throughput_scaling.pdf'), bbox_inches='tight')
    print("Saved fig2_throughput_scaling.png/pdf")
    plt.close()


def fig3_vep_concordance():
    """Figure 3: VEP concordance heatmap/bar chart."""
    data = read_csv('vep_concordance.csv')
    fields = [d['field'] for d in data]
    accuracy = [float(d['accuracy']) for d in data]

    if not HAS_MPL:
        print("Figure 3: VEP Concordance")
        for f, a in zip(fields, accuracy):
            print(f"  {f:25s} {a:.1f}%")
        print()
        return

    fig, ax = plt.subplots(figsize=(10, 7))
    colors = [COLORS['success'] if a == 100.0 else COLORS['vep'] for a in accuracy]
    bars = ax.barh(range(len(fields)), accuracy, color=colors, alpha=0.85)
    ax.set_yticks(range(len(fields)))
    ax.set_yticklabels(fields, fontsize=10, fontfamily='monospace')
    ax.invert_yaxis()
    ax.set_xlabel('Concordance with Ensembl VEP v115.1 (%)', fontsize=12)
    ax.set_title('Field-Level Accuracy: OxiVEP vs Ensembl VEP\n(2,340 shared transcript-allele pairs, 173 variants)', fontsize=13, fontweight='bold')
    ax.set_xlim(95, 101)
    ax.grid(True, alpha=0.3, axis='x')
    ax.axvline(x=100, color=COLORS['success'], linestyle='--', alpha=0.5, linewidth=1)

    for bar, acc in zip(bars, accuracy):
        ax.text(bar.get_width() - 0.3, bar.get_y() + bar.get_height()/2,
                f'{acc:.1f}%', va='center', ha='right', fontsize=9, color='white', fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig3_vep_concordance.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig3_vep_concordance.pdf'), bbox_inches='tight')
    print("Saved fig3_vep_concordance.png/pdf")
    plt.close()


def fig4_consequence_distribution():
    """Figure 4: Consequence distribution from real 1KGP data."""
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
            bar = '#' * (n // 200)
            print(f"  {c:45s} {n:>6d}  {bar}")
        print()
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

    # Legend for impact colors
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=COLORS['high'], alpha=0.85, label='HIGH'),
        Patch(facecolor=COLORS['moderate'], alpha=0.85, label='MODERATE'),
        Patch(facecolor=COLORS['low'], alpha=0.85, label='LOW'),
        Patch(facecolor=COLORS['modifier'], alpha=0.85, label='MODIFIER'),
    ]
    ax.legend(handles=legend_elements, title='Impact', loc='lower right', fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig4_consequence_distribution.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig4_consequence_distribution.pdf'), bbox_inches='tight')
    print("Saved fig4_consequence_distribution.png/pdf")
    plt.close()


def fig5_resource_usage():
    """Figure 5: Resource usage comparison."""
    data = read_csv('resource_usage.csv')
    mem_data = {d['metric']: float(d['value']) for d in data}

    if not HAS_MPL:
        print("Figure 5: Resource Usage")
        for k, v in mem_data.items():
            unit = next((d['unit'] for d in data if d['metric'] == k), '')
            print(f"  {k}: {v} {unit}")
        print()
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: Memory comparison
    variants_labels = ['1K', '10K', '100K']
    oxivep_mem = [mem_data.get('peak_memory_1000v', 0),
                  mem_data.get('peak_memory_10000v', 0),
                  mem_data.get('peak_memory_100000v', 0)]
    vep_mem = [500, 500, 600]

    x = range(len(variants_labels))
    w = 0.35
    ax1.bar([i - w/2 for i in x], vep_mem, w, label='Ensembl VEP', color=COLORS['vep'], alpha=0.8)
    ax1.bar([i + w/2 for i in x], oxivep_mem, w, label='OxiVEP', color=COLORS['primary'], alpha=0.8)
    ax1.set_xlabel('Input size (variants)', fontsize=12)
    ax1.set_ylabel('Peak memory (MB)', fontsize=12)
    ax1.set_title('A. Memory Usage', fontsize=13, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(variants_labels)
    ax1.legend()
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3, axis='y')

    # Panel B: Binary size / startup
    metrics = ['Binary Size\n(MB)', 'Startup + GFF3\n(ms)', 'Startup\n(cached, ms)']
    vep_vals = [200, 10000, 10000]
    oxi_vals = [mem_data.get('binary_size', 0),
                mem_data.get('startup_time_gff3_chr22', 0),
                mem_data.get('startup_time_cached', 0)]

    x = range(len(metrics))
    ax2.bar([i - w/2 for i in x], vep_vals, w, label='Ensembl VEP', color=COLORS['vep'], alpha=0.8)
    ax2.bar([i + w/2 for i in x], oxi_vals, w, label='OxiVEP', color=COLORS['primary'], alpha=0.8)
    ax2.set_ylabel('Value', fontsize=12)
    ax2.set_title('B. Size & Startup', fontsize=13, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(metrics, fontsize=10)
    ax2.legend()
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig5_resource_usage.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig5_resource_usage.pdf'), bbox_inches='tight')
    print("Saved fig5_resource_usage.png/pdf")
    plt.close()


def fig1_architecture():
    """Figure 1: Architecture diagram (text-based)."""
    diagram = """
    ┌─────────────────────────────────────────────────────────────────┐
    │                        oxivep-cli                               │
    │    (CLI binary, pipeline, cache builder, embedded web server)    │
    └──────┬──────────┬──────────┬──────────┬──────────┬─────────────┘
           │          │          │          │          │
    ┌──────▼──┐ ┌─────▼─────┐ ┌─▼────────┐ │   ┌──────▼──────┐
    │oxivep-io│ │oxivep-hgvs│ │oxivep-   │ │   │oxivep-filter│
    │(VCF I/O,│ │(HGVSg,    │ │consequen-│ │   │(filtering)  │
    │ 47-field│ │ HGVSc,    │ │ce engine)│ │   └─────────────┘
    │  CSQ)   │ │ HGVSp)    │ └──┬───────┘ │
    └────┬────┘ └─────┬─────┘    │         │
         │            │          │    ┌────▼─────────┐
    ┌────▼────────────▼──────────▼──┐ │oxivep-cache  │
    │        oxivep-genome          │ │(GFF3, FASTA, │
    │  (Transcript, Exon, Gene,     │ │ tabix, mmap, │
    │   CodonTable, coord mapping)  │ │ binary cache)│
    └────────────┬──────────────────┘ └────┬─────────┘
                 │                         │
            ┌────▼─────────────────────────▼──┐
            │          oxivep-core            │
            │  (GenomicPosition, Consequence, │
            │   Allele, Strand, Impact)       │
            └─────────────────────────────────┘
    """
    with open(os.path.join(FIG_DIR, 'fig1_architecture.txt'), 'w') as f:
        f.write(diagram)
    print("Saved fig1_architecture.txt")


if __name__ == '__main__':
    print("Generating OxiVEP manuscript figures...")
    print(f"Data dir: {os.path.abspath(DATA_DIR)}")
    print(f"Output dir: {os.path.abspath(FIG_DIR)}")
    print()
    fig1_architecture()
    fig2_throughput_comparison()
    fig3_vep_concordance()
    fig4_consequence_distribution()
    fig5_resource_usage()
    print("\nDone.")
