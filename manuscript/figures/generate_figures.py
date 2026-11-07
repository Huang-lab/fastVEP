#!/usr/bin/env python3
"""Generate the manuscript figures for fastVEP from `manuscript/data/` and from nothing else.

Every number a figure prints is read from a CSV there. `manuscript/build_data.py` writes the CSVs from the runs that measured them
(`benchmarks/results/<build>/`), so a figure cannot quote a measurement that no run made: the labels that used to be typed here
("4,621 s", "198 s", "93 s", "22/22 fields", the stats strip of Figure 1) are computed, and so are the limits of the axes that they
would have run off when the numbers changed.
"""

import csv
import math
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
    print("matplotlib not installed.")

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')
FIG_DIR = os.path.dirname(__file__)

C = {
    'fast': '#6c7aee', 'vep': '#f5426c', 'success': '#10b981',
    'high': '#f5426c', 'moderate': '#f59e0b', 'low': '#3b82f6', 'modifier': '#6b7280',
    'cli': '#6c7aee', 'mid': '#818cf8', 'data': '#34d399', 'core': '#f59e0b', 'sa': '#f472b6',
}
OC = {'yeast': '#10b981', 'drosophila': '#f59e0b', 'arabidopsis': '#8b5cf6', 'mouse': '#ef4444', 'human': '#3b82f6'}
OL = {'yeast': 'Yeast', 'drosophila': 'Drosophila', 'arabidopsis': 'Arabidopsis', 'mouse': 'Mouse', 'human': 'Human'}


def read_csv(f):
    with open(os.path.join(DATA_DIR, f)) as fh:
        return list(csv.DictReader(fh))


def metrics():
    return {r['metric']: r['value'] for r in read_csv('code_metrics.csv')}


def resources():
    return {r['metric']: r['value'] for r in read_csv('resource_usage.csv')}


# ---------------------------------------------------------------------------------------------------------------------
# Figure 1: the architecture is drawn from the dependencies cargo reports, not from memory.

CRATE_ROLE = {
    'core': 'Consequence (SO terms) · Allele · Impact · VariantType',
    'genome': 'Transcript / Exon / Gene · codon tables',
    'filter': 'filter_vep expressions',
    'classification': 'ACMG/AMP criteria',
    'cache': 'GFF3 · FASTA mmap · transcript cache',
    'consequence': 'SNV / indel / SV predictors',
    'hgvs': 'HGVS g / c / p',
    'io': 'VCF in · CSQ / tab / JSON out',
    'sa': 'fastSA · ClinVar, gnomAD, scores',
    'annotate': 'shared engine · pick · web servers',
    'cli': 'annotate · sa-build · sa-convert · filter · cache · web',
    'web': 'axum production server',
}
CRATE_COLOR = {'cli': C['cli'], 'web': C['cli'], 'annotate': C['mid'], 'io': C['mid'], 'consequence': C['mid'], 'hgvs': C['mid'],
               'filter': C['mid'], 'classification': C['mid'], 'genome': C['data'], 'cache': C['data'], 'sa': C['sa'], 'core': C['core']}


def crate_layers():
    """{crate: set of direct dependencies} with the edges that another path already implies removed, and the layer of each crate."""
    deps = {}
    for r in read_csv('crate_dependencies.csv'):
        deps.setdefault(r['crate'].replace('fastvep-', ''), set()).add(r['depends_on'].replace('fastvep-', ''))
    crates = set(deps) | {d for ds in deps.values() for d in ds}
    for c in crates:
        deps.setdefault(c, set())

    def reach(c, seen=None):
        out = set()
        for d in deps[c]:
            out |= {d} | reach(d)
        return out
    reduced = {c: {d for d in deps[c] if not any(d in reach(o) for o in deps[c] if o != d)} for c in crates}
    layer = {}

    def depth(c):
        if c not in layer:
            layer[c] = 0 if not deps[c] else 1 + max(depth(d) for d in deps[c])
        return layer[c]
    for c in crates:
        depth(c)
    return reduced, layer


def layer_order(reduced, layer):
    """The order of the crates within each layer, by the barycentre of their neighbours, which keeps the arrows from crossing more than they must."""
    n_layers = max(layer.values()) + 1
    order = {k: sorted((c for c in layer if layer[c] == k), key=lambda c: list(CRATE_ROLE).index(c)) for k in range(n_layers)}
    neighbours = {c: set() for c in layer}
    for c, ds in reduced.items():
        for d in ds:
            neighbours[c].add(d)
            neighbours[d].add(c)

    def x_of(c):
        members = order[layer[c]]
        return (members.index(c) + 0.5) / len(members)

    def crossings():
        edges = [(c, d) for c, ds in reduced.items() for d in ds]
        n = 0
        for i, (c1, d1) in enumerate(edges):
            for c2, d2 in edges[i + 1:]:
                if len({c1, d1, c2, d2}) < 4:
                    continue
                # edges between the same two layers cross when their ends are in opposite orders; an edge that skips a layer is judged on its ends alone
                if (x_of(c1) - x_of(c2)) * (x_of(d1) - x_of(d2)) < 0:
                    n += 1
        return n

    best, best_n = {k: list(v) for k, v in order.items()}, crossings()
    for sweep in range(8):
        for k in (range(n_layers - 1, -1, -1) if sweep % 2 == 0 else range(n_layers)):
            order[k] = sorted(order[k], key=lambda c: (sum(x_of(n) for n in neighbours[c]) / len(neighbours[c])) if neighbours[c] else x_of(c))
        n = crossings()
        if n < best_n:
            best, best_n = {k: list(v) for k, v in order.items()}, n
    return best


def fig1_architecture():
    """Architecture: (A) the crate dependency graph with test counts; (B) annotation data flow."""
    m, rs = metrics(), resources()
    reduced, layer = crate_layers()
    if not HAS_MPL: return
    n_layers = max(layer.values()) + 1
    order = layer_order(reduced, layer)
    fig = plt.figure(figsize=(18, 9.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.75, 1], wspace=0.06)
    axA = fig.add_subplot(gs[0]); axB = fig.add_subplot(gs[1])

    # ================= Panel A: crates by layer =================
    axA.set_xlim(0, 14); axA.set_ylim(0.1, 9.3); axA.axis('off')
    H = 1.0
    ys = {k: 0.45 + (8.0 - 0.45) * k / (n_layers - 1) for k in range(n_layers)}
    axA.text(7, 9.15, f"A. Crate architecture: {len(layer)}-crate Cargo workspace (an arrow is a direct dependency)", ha='center', fontsize=13, fontweight='bold', color='#1f2937')
    pos, width = {}, {}
    for k in range(n_layers):
        members = order[k]
        w = min(3.3, 12.6 / len(members) - 0.25)
        total = len(members) * w + (len(members) - 1) * 0.25
        x = 7 - total / 2
        axA.add_patch(plt.Rectangle((0, ys[k] - 0.14), 14, H + 0.28, facecolor='#f8fafc' if k % 2 else '#eef2ff', edgecolor='none', zorder=0))
        for c in members:
            tests = m.get(f'tests_fastvep-{c}', '?')
            axA.add_patch(FancyBboxPatch((x, ys[k]), w, H, boxstyle="round,pad=0.05", facecolor=CRATE_COLOR[c], edgecolor='#374151', lw=1.2, alpha=0.95, zorder=2))
            axA.text(x + w / 2, ys[k] + H - 0.2, f"fastvep-{c}", ha='center', va='center', fontsize=9.6 if w > 2.4 else 8.8, fontweight='bold', color='white', zorder=3)
            # the role breaks at its own " · " separators and not inside a phrase
            role, line = [], ''
            for part in CRATE_ROLE[c].split(' · '):
                if line and len(line) + 3 + len(part) > int(w * 11.5):
                    role.append(line)
                    line = part
                else:
                    line = f"{line} · {part}" if line else part
            role.append(line)
            assert len(role) <= 2, (c, role)
            for i, line in enumerate(role):
                axA.text(x + w / 2, ys[k] + H - 0.46 - 0.18 * i, line, ha='center', va='center', fontsize=6.9, color='white', alpha=0.97, style='italic', zorder=3)
            axA.text(x + w / 2, ys[k] + 0.14, f"{tests} tests", ha='center', va='center', fontsize=7.4, color='white', fontweight='bold', alpha=0.97, zorder=3)
            pos[c], width[c] = (x + w / 2, ys[k]), w
            x += w + 0.25
    # the ends of the arrows are spread along the edges of a box, in the order of the boxes they join, so that they do not meet in one point
    out_edges = {c: sorted(ds, key=lambda d: pos[d][0]) for c, ds in reduced.items()}
    in_edges = {c: sorted((o for o, ds in reduced.items() if c in ds), key=lambda o: pos[o][0]) for c in layer}
    spread = lambda i, n, w: 0.0 if n == 1 else (-0.36 + 0.72 * i / (n - 1)) * w
    for c, ds in out_edges.items():
        for i, d in enumerate(ds):
            x1 = pos[c][0] + spread(i, len(ds), width[c])
            j = in_edges[d].index(c)
            x2 = pos[d][0] + spread(j, len(in_edges[d]), width[d])
            axA.annotate('', xy=(x2, pos[d][1] + H + 0.06), xytext=(x1, pos[c][1] - 0.06),
                         arrowprops=dict(arrowstyle='->', color='#475569', lw=1.1, alpha=0.8, shrinkA=0, shrinkB=0), zorder=1)

    # ================= Panel B: annotation data flow =================
    axB.set_xlim(0, 10); axB.set_ylim(0.1, 9.3); axB.axis('off')
    axB.text(5, 9.15, 'B. Annotation data flow', ha='center', fontsize=13, fontweight='bold', color='#1f2937')
    stages = [('VCF input', 'multi-sample · SV · gzip', '#64748b'),
              ('Parse & normalize', 'fastvep-io', C['mid']),
              ('Transcript overlap', 'fastvep-cache · binary search', C['data']),
              ('Consequence prediction', f"fastvep-consequence · {m['so_consequence_types_emitted']} SO terms emitted", C['mid']),
              ('HGVS g / c / p', 'fastvep-hgvs', C['mid']),
              ('Supplementary annotation', 'fastvep-sa · .osa / .osa2', C['sa']),
              ('Output', f"CSQ ({m['csq_output_fields']} fields) · TSV · JSON", '#64748b')]
    n = len(stages); top = 8.1; bot = 1.5; step = (top - bot) / (n - 1)
    for i, (t, s, col) in enumerate(stages):
        y = top - i * step
        axB.add_patch(FancyBboxPatch((1.1, y - 0.34), 7.8, 0.68, boxstyle="round,pad=0.1", facecolor=col, edgecolor='#374151', lw=1.2, alpha=0.92))
        axB.text(5, y + 0.1, t, ha='center', va='center', fontsize=10, fontweight='bold', color='white')
        axB.text(5, y - 0.17, s, ha='center', va='center', fontsize=7.3, color='white', alpha=0.92, style='italic')
        if i < n - 1: axB.annotate('', xy=(5, y - step + 0.36), xytext=(5, y - 0.38), arrowprops=dict(arrowstyle='->', color='#6b7280', lw=1.8))

    # stats strip across the bottom
    stats = (f"{m['crates']} crates    ·    {int(m['total_source_lines']):,} lines of Rust    ·    {int(m['tests_total']):,} tests    ·    "
             f"{m['so_consequence_types']} SO terms ({m['so_consequence_types_emitted']} emitted)    ·    {m['csq_output_fields']} CSQ fields    ·    "
             f"{rs['binary_size']} MB binary    ·    no run-time dependency but libc    ·    fastVEP.org")
    fig.text(0.5, 0.02, stats, ha='center', fontsize=10.5, color='#374151', bbox=dict(boxstyle='round,pad=0.45', facecolor='#f3f4f6', edgecolor='#d1d5db'))
    plt.subplots_adjust(bottom=0.09, top=0.97)
    for ext in ['png', 'pdf']: plt.savefig(os.path.join(FIG_DIR, f'fig1_architecture.{ext}'), dpi=300, bbox_inches='tight')
    print("Saved fig1_architecture"); plt.close()


# ---------------------------------------------------------------------------------------------------------------------
def fig2_comparison():
    """THE figure: fastVEP vs Ensembl VEP, single-threaded, measured genome-wide."""
    h2h = read_csv('vep_comparison.csv')       # chr22 scaling (single thread)
    gw = read_csv('vep_genomewide.csv')        # genome-wide 3-condition head-to-head
    orgs = read_csv('organism_comparison.csv')
    giab = read_csv('giab_benchmark.csv')[0]
    if not HAS_MPL: return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.5))

    # === Panel A: Wall-clock time (log-log), single thread ===
    hv = [int(d['variants']) for d in h2h]
    ft = [float(d['fastvep_sec']) for d in h2h]
    vt = [float(d['vep_sec']) for d in h2h]
    ax1.plot(hv, ft, 'o-', color=C['fast'], lw=2.5, ms=9, label='fastVEP (chr22, 1 thread)', zorder=3)
    ax1.plot(hv, vt, 's-', color=C['vep'], lw=2.5, ms=9, label='Ensembl VEP (chr22, 1 thread)', zorder=3)

    # the crossover is where the two curves meet, found by interpolating the log of the ratio between the two points around it
    cross = None
    for i in range(1, len(hv)):
        r0, r1 = math.log(vt[i - 1] / ft[i - 1]), math.log(vt[i] / ft[i])
        if r0 < 0 <= r1:
            f = -r0 / (r1 - r0)
            cross = math.exp(math.log(hv[i - 1]) + f * (math.log(hv[i]) - math.log(hv[i - 1])))
            lo, hi = hv[i - 1], hv[i]
            ax1.axvspan(lo, hi, color='#94a3b8', alpha=0.08, zorder=0)
            ax1.annotate(f'crossover\n(~{cross:,.0f} variants)', xy=(math.sqrt(lo * hi), 0.17), fontsize=7.5, color='#475569', ha='center', va='bottom', style='italic')
            break

    # genome-wide MEASURED points (VEP = per-chromosome sum)
    NGW = int(giab['variants'])
    fv1t = float(gw[0]['fastvep_1t_sec']); fvmt = float(gw[0]['fastvep_mt_sec']); vepgw = float(gw[0]['vep_sec'])
    ax1.plot([hv[-1], NGW], [vt[-1], vepgw], color=C['vep'], lw=1, ls=':', alpha=0.5, zorder=2)
    ax1.plot([hv[-1], NGW], [ft[-1], fv1t], color=C['fast'], lw=1, ls=':', alpha=0.5, zorder=2)
    ax1.scatter([NGW], [vepgw], marker='s', s=150, color=C['vep'], edgecolors='white', lw=1.5, zorder=5)
    ax1.scatter([NGW], [fv1t], marker='o', s=150, color=C['fast'], edgecolors='white', lw=1.5, zorder=5)
    ax1.scatter([NGW], [fvmt], marker='*', s=300, color=C['fast'], edgecolors='white', lw=1.2, zorder=6)
    ax1.annotate(f'VEP whole genome\n{vepgw:,.0f} s (per-chr sum)', xy=(NGW, vepgw), xytext=(-12, -4),
                 textcoords='offset points', fontsize=7.5, color=C['vep'], ha='right', va='top', fontweight='bold')
    ax1.annotate(f'fastVEP 1-thread\n{fv1t:.0f} s', xy=(NGW, fv1t), xytext=(-12, 8),
                 textcoords='offset points', fontsize=7.5, color=C['fast'], ha='right')
    ax1.annotate(f'fastVEP multi-thread {fvmt:.0f} s', xy=(NGW, fvmt), xytext=(-12, -16),
                 textcoords='offset points', fontsize=7.5, color=C['fast'], ha='right', fontweight='bold')

    # fastVEP multi-organism points
    for d in orgs:
        o = d['organism']; v, t = int(d['variants']), float(d['time_sec'])
        if o == 'human':
            continue  # shown as the whole-genome markers above
        ax1.scatter(v, t, c=OC.get(o, '#999'), s=110, marker='D', alpha=0.95, edgecolors='white', lw=1.0, zorder=4)
        label = f"{OL.get(o,o)} {v/1e6:.1f}M" if v > 500000 else f"{OL.get(o,o)} {v/1e3:.0f}K"
        near_star = abs(math.log10(v / NGW)) < 0.1
        ax1.annotate(label, xy=(v, t), xytext=(9, -11) if near_star else (7, 5), textcoords='offset points', fontsize=7, color='#333')

    ax1.set_xscale('log'); ax1.set_yscale('log')
    ax1.set_ylim(0.14, vepgw * 2.2)
    ax1.set_xlabel('Variants', fontsize=12)
    ax1.set_ylabel('Wall-clock time (seconds)', fontsize=12)
    ax1.set_title('A. Annotation time: fastVEP vs Ensembl VEP (single-thread)', fontsize=12.5, fontweight='bold')
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1e6:.0f}M' if x >= 1e6 else f'{x/1e3:.0f}K'))
    ax1.grid(True, alpha=0.15)
    legend_h = [
        plt.Line2D([0],[0],marker='o',color=C['fast'],lw=2,ms=8,label='fastVEP (chr22, 1 thread)'),
        plt.Line2D([0],[0],marker='s',color=C['vep'],lw=2,ms=8,label='Ensembl VEP (1 thread)'),
        plt.Line2D([0],[0],marker='*',color=C['fast'],lw=0,ms=13,label='fastVEP multi-thread'),
        plt.Line2D([0],[0],marker='D',color='#999',lw=0,ms=8,label='fastVEP (full organisms)'),
    ]
    ax1.legend(handles=legend_h, fontsize=8.5, loc='upper left')

    # === Panel B: Whole-genome single-thread head-to-head (the headline) ===
    labels = [d['label'] for d in gw]
    fv = [float(d['fastvep_1t_sec']) for d in gw]
    vp = [float(d['vep_sec']) for d in gw]
    sp = [float(d['speedup']) for d in gw]
    x = list(range(len(labels))); w = 0.38
    ax2.bar([i - w/2 for i in x], fv, w, color=C['fast'], label='fastVEP (1 thread)')
    ax2.bar([i + w/2 for i in x], vp, w, color=C['vep'], label='Ensembl VEP (1 thread)')
    ax2.set_yscale('log')
    ax2.set_xticks(x); ax2.set_xticklabels(labels, fontsize=10)
    ax2.set_xlabel('Annotation workload (complete human WGS)', fontsize=12)
    ax2.set_ylabel('Wall-clock time, full WGS (s, log scale)', fontsize=11)
    ax2.set_title(f'B. Whole genome ({NGW/1e6:.2f}M variants): single-thread head-to-head', fontsize=12.5, fontweight='bold')
    ax2.grid(True, alpha=0.15, axis='y')
    ax2.set_ylim(min(fv) / 2.2, max(vp) * 20)
    for i in x:
        ax2.text(i - w/2, fv[i]*1.05, f'{fv[i]:.0f}s', ha='center', va='bottom', fontsize=8.5, color=C['fast'], fontweight='bold')
        ax2.text(i + w/2, vp[i]*1.05, f'{vp[i]:.0f}s', ha='center', va='bottom', fontsize=8.5, color=C['vep'], fontweight='bold')
        ax2.text(i, max(vp[i], fv[i])*2.1, f'{sp[i]:.0f}x' if sp[i] >= 10 else f'{sp[i]:.1f}x', ha='center', va='bottom', fontsize=13, fontweight='bold', color='#1f2937')
    mt = [float(d['fastvep_mt_sec']) for d in gw]
    handles, labs = ax2.get_legend_handles_labels()
    handles.append(plt.Line2D([0], [0], marker='*', color=C['fast'], lw=0, ms=12))
    labs.append(f'fastVEP, all 10 cores: {min(mt):.0f}-{max(mt):.0f} s (not a bar)')
    ax2.legend(handles, labs, fontsize=9, loc='upper right', framealpha=0.95)

    plt.tight_layout()
    for ext in ['png','pdf']: plt.savefig(os.path.join(FIG_DIR, f'fig2_comparison.{ext}'), dpi=300, bbox_inches='tight')
    print("Saved fig2_comparison"); plt.close()


def fig3_concordance():
    """Agreement with Ensembl VEP 116.3: the cells that differ on each input, and, on the HG002 sample, by column and by cause.

    The counts are `manuscript/data/vep116_agreement.csv` and `vep116_hg002_columns.csv`, which are the counts of Part 3 of
    docs/VEP_DIVERGENCE.md (`vep116_provenance.md` says where they come from and why they are those of the release).
    """
    inputs = read_csv('vep116_agreement.csv')
    cols = read_csv('vep116_hg002_columns.csv')
    if not HAS_MPL: return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.5, 7), gridspec_kw={'width_ratios': [1, 1.05]})

    # === Panel A: share of cells that differ, by input (log scale; a dot for each, because a bar on a log axis has the length of its baseline) ===
    label = {'hg002_1in200': 'HG002 call set,\n1-in-200 sample', 'clinvar_1in5': 'ClinVar 2-star+,\n1-in-5, no HGVS columns',
             'clinvar_1in5_hgvs': 'ClinVar 2-star+,\n1-in-5, with HGVS columns', 'world': 'Generated world\n(both strands, multi-ALT)',
             'edges': 'Its edges\n(every shape at every boundary)', 'lone_indels': 'Its lone indels\n(shared bases)'}
    colour = {'hg002_1in200': C['fast'], 'clinvar_1in5': C['fast'], 'clinvar_1in5_hgvs': C['fast'],
              'world': C['moderate'], 'edges': C['moderate'], 'lone_indels': C['moderate']}
    share = [100.0 * int(r['cells_differ']) / int(r['cells_compared']) for r in inputs]
    y = list(range(len(inputs)))
    ax1.hlines(y, 1e-4, share, color=[colour[r['input']] for r in inputs], lw=2, alpha=0.45, zorder=1)
    ax1.scatter(share, y, s=130, color=[colour[r['input']] for r in inputs], edgecolors='white', lw=1.2, zorder=3)
    for yi, r, v in zip(y, inputs, share):
        ax1.text(v * 1.7, yi, f"{int(r['cells_differ']):,} of {int(r['cells_compared']):,} cells ({v:.2g}%)", va='center', fontsize=9, color='#333')
    ax1.set_yticks(y); ax1.set_yticklabels([label[r['input']] for r in inputs], fontsize=9.5)
    ax1.invert_yaxis()
    ax1.set_xscale('log'); ax1.set_xlim(1e-4, 100)
    ax1.set_xlabel('Cells that differ from Ensembl VEP 116.3 (% of the cells compared, log scale)', fontsize=11)
    ax1.set_title('A. Cells that differ, by input', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.2, axis='x')
    ax1.legend(handles=[plt.Line2D([0], [0], marker='o', color=C['fast'], lw=0, ms=9, label='real call sets'),
                        plt.Line2D([0], [0], marker='o', color=C['moderate'], lw=0, ms=9, label='generated worlds, built to be hard')],
               fontsize=9.5, loc='center right', bbox_to_anchor=(1.0, 0.68))

    # === Panel B: the HG002 sample by column, stacked by cause ===
    names = [c['column'] for c in cols]
    yb = list(range(len(names)))
    parts = [('divergence_11_several_alts', C['moderate'], 'alleles of records with several ALTs (divergence 11)'),
             ('divergence_4_terminator_intron', C['high'], 'an intronic base after the terminator (divergence 4)'),
             ('divergence_17_three_prime_rule', C['low'], "a change the 3'-rule carries over an intron (divergence 17)")]
    left = [0] * len(names)
    for key, colr, lab in parts:
        vals = [int(c[key]) for c in cols]
        ax2.barh(yb, vals, left=left, color=colr, alpha=0.9, label=lab)
        left = [l + v for l, v in zip(left, vals)]
    rows_n = int(cols[0]['paired_rows'])
    for yi, tot in zip(yb, left):
        ax2.text(tot + 25 if tot else 25, yi, f'{tot:,}' if tot else '0', va='center', fontsize=9, color='#333' if tot else C['success'], fontweight='bold' if tot == 0 else 'normal')
    ax2.set_yticks(yb); ax2.set_yticklabels(names, fontsize=9.5, fontfamily='monospace')
    ax2.invert_yaxis()
    ax2.set_xlim(0, max(left) * 1.18)
    ax2.set_xlabel(f'Cells that differ (of {rows_n:,} paired rows per column)', fontsize=11)
    total = sum(left)
    ax2.set_title(f'B. HG002 sample, by column: {total:,} cells differ\n(17 of the 21 columns shown agree on every row)', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.2, axis='x')
    ax2.legend(title='every differing cell has a documented cause:', title_fontsize=9, fontsize=9, loc='upper right')

    plt.tight_layout()
    for ext in ['png','pdf']: plt.savefig(os.path.join(FIG_DIR, f'fig3_concordance.{ext}'), dpi=300, bbox_inches='tight')
    print("Saved fig3_concordance"); plt.close()


def fig4_consequences():
    """Consequence distribution (genome-wide, GIAB HG002 full WGS), each call under its most severe term."""
    data = read_csv('consequence_distribution.csv')
    impact = {r['so_term']: r['impact'] for r in read_csv('so_terms.csv')}
    giab = read_csv('giab_benchmark.csv')[0]
    if not HAS_MPL: return
    consequences = [d['consequence'] for d in data]
    counts = [int(d['count']) for d in data]
    colors = [C[impact[c].lower()] for c in consequences]   # a KeyError means a term that so_terms.csv does not know
    fig, ax = plt.subplots(figsize=(11, 8))
    bars = ax.barh(range(len(consequences)), counts, color=colors, alpha=0.85)
    ax.set_yticks(range(len(consequences))); ax.set_yticklabels(consequences, fontsize=9.5, fontfamily='monospace')
    ax.invert_yaxis(); ax.set_xlabel('Consequence calls (log scale)', fontsize=12)
    ax.set_title(f"Predicted Consequence Distribution\n(GIAB HG002 full WGS: {int(giab['variants'])/1e6:.2f}M variants, {sum(counts)/1e6:.1f}M calls, {len(consequences)} SO terms as the most severe)", fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x'); ax.set_xscale('log')
    for bar, count in zip(bars, counts):
        ax.text(bar.get_width()*1.15, bar.get_y()+bar.get_height()/2, f'{count:,}', va='center', fontsize=8.5, color='#555')
    ax.set_xlim(right=ax.get_xlim()[1]*3)
    legend_el = [mpatches.Patch(facecolor=C[k],alpha=0.85,label=l) for k,l in [('high','HIGH'),('moderate','MODERATE'),('low','LOW'),('modifier','MODIFIER')]]
    ax.legend(handles=legend_el, title='Impact', loc='lower right', fontsize=9)
    plt.tight_layout()
    for ext in ['png','pdf']: plt.savefig(os.path.join(FIG_DIR, f'fig4_consequences.{ext}'), dpi=300, bbox_inches='tight')
    print("Saved fig4_consequences"); plt.close()


def fig5_organisms():
    """Multi-organism: throughput consistency + absolute scale handled."""
    data = read_csv('organism_comparison.csv')
    if not HAS_MPL: return

    orgs = [d['organism'] for d in data]
    vps = [int(d['variants_per_sec']) for d in data]
    variants = [int(d['variants']) for d in data]
    times = [float(d['time_sec']) for d in data]
    transcripts = [int(d['transcripts']) for d in data]

    # order by transcript count (genome complexity) for both panels
    idx = sorted(range(len(orgs)), key=lambda i: transcripts[i])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.4))

    # === Panel A: throughput vs genome complexity (consistency) ===
    tx = [transcripts[i] for i in idx]; vp = [vps[i] for i in idx]
    cols = [OC.get(orgs[i], '#999') for i in idx]
    ax1.axhspan(min(vps), max(vps), color='#94a3b8', alpha=0.10, zorder=0)
    ax1.plot(tx, vp, '-', color='#cbd5e1', lw=1.5, zorder=1)
    ax1.scatter(tx, vp, c=cols, s=260, edgecolors='white', lw=1.8, zorder=3)
    # label offsets are given by the order along the curve, so that adjacent organisms do not collide
    shifts = [(6, 16, 'left'), (-42, 16, 'center'), (44, -30, 'center'), (0, 16, 'center'), (-6, 18, 'right')]
    for rank, i in enumerate(idx):
        dx, dy, ha = shifts[rank % len(shifts)]
        ax1.annotate(f"{OL.get(orgs[i],orgs[i])}\n{vps[i]:,} v/s",
                     xy=(transcripts[i], vps[i]), xytext=(dx, dy), textcoords='offset points',
                     ha=ha, fontsize=8.5, fontweight='bold', color='#333')
    ax1.set_xscale('log')
    ax1.set_xlabel('Transcripts (genome complexity, log scale)', fontsize=12)
    ax1.set_ylabel('Throughput (variants/sec)', fontsize=12)
    ax1.set_title(f'A. Throughput across a {max(transcripts)/min(transcripts):.0f}x range\nof genome complexity ({min(transcripts)/1000:.0f}K-{max(transcripts)/1000:.0f}K transcripts)', fontsize=12, fontweight='bold')
    ax1.set_ylim(0, math.ceil(max(vps) * 1.25 / 100000) * 100000)
    ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1000:.0f}K'))
    ax1.grid(True, alpha=0.15)
    ax1.text(0.5, 0.06, f'{min(vps)/1000:.0f}-{max(vps)/1000:.0f}K variants/sec band', transform=ax1.transAxes, ha='center',
             fontsize=8.5, style='italic', color='#64748b')

    # === Panel B: absolute scale: complete datasets, wall-clock ===
    yb = range(len(idx))
    nvar = [variants[i] for i in idx]
    ax2.barh(list(yb), nvar, color=cols, alpha=0.88, height=0.62)
    ax2.set_yticks(list(yb))
    ax2.set_yticklabels([f"{OL.get(orgs[i],orgs[i])}" for i in idx], fontsize=11)
    ax2.set_xscale('log'); ax2.set_xlim(1e5, 5e8)
    ax2.set_xlabel('Variants annotated (complete dataset, log scale)', fontsize=12)
    ax2.set_title('B. Complete datasets annotated\n(wall-clock time per dataset)', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.15, axis='x')
    for i, ii in enumerate(idx):
        t = times[ii]
        tstr = f"{t:.1f}s" if t < 10 else (f"{t:.0f}s" if t < 100 else f"{t/60:.1f} min")
        ax2.text(variants[ii]*1.25, i, f"{variants[ii]/1e6:.2f}M var  ·  {tstr}",
                 va='center', fontsize=9, fontweight='bold', color='#333')
    fig.text(0.5, 0.005, 'The Drosophila and mouse variant files are in the coordinates of dm3 and GRCm38, against gene models of BDGP6.54 and GRCm39: '
             'those two runs measure throughput, not annotation.', ha='center', fontsize=8, style='italic', color='#64748b')

    plt.tight_layout(rect=(0, 0.03, 1, 1))
    for ext in ['png','pdf']: plt.savefig(os.path.join(FIG_DIR, f'fig5_organisms.{ext}'), dpi=300, bbox_inches='tight')
    print("Saved fig5_organisms"); plt.close()


if __name__ == '__main__':
    print("Generating fastVEP manuscript figures...\n")
    fig1_architecture()
    fig2_comparison()
    fig3_concordance()
    fig4_consequences()
    fig5_organisms()
    print("\nDone.")
