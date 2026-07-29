#!/usr/bin/env python3
"""Generate manuscript figures for fastVEP."""

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


def fig1_architecture():
    """Architecture: (A) tiered crate dependency graph with test counts; (B) annotation data flow."""
    if not HAS_MPL: return
    tests = {'core':20,'genome':18,'cache':54,'consequence':27,'hgvs':25,'io':86,'filter':21,
             'sa':153,'annotate':5,'classification':147,'cli':76,'web':9}
    fig = plt.figure(figsize=(18, 8.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.75, 1], wspace=0.06)
    axA = fig.add_subplot(gs[0]); axB = fig.add_subplot(gs[1])

    # ================= Panel A: crate architecture (centered rows) =================
    axA.set_xlim(0, 14); axA.set_ylim(1.0, 9.0); axA.axis('off')
    CX, H = 7.0, 0.85
    tier_y = {'bin':7.0,'lib':5.3,'prov':3.6,'core':1.9}
    for key,lbl,col in [('bin','binaries','#eef2ff'),('lib','libraries','#eef2ff'),
                        ('prov','data providers','#ecfdf5'),('core','core','#fff7ed')]:
        y=tier_y[key]
        axA.add_patch(plt.Rectangle((0,y-0.14),14,H+0.28,facecolor=col,edgecolor='none',zorder=0))
        axA.text(0.42,y+H/2,lbl,ha='center',va='center',fontsize=8,color='#9aa3b2',style='italic',rotation=90)
    def arr(x1,y1,x2,y2,c='#94a3b8'): axA.annotate('',xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle='->',color=c,lw=1.1,alpha=0.85),zorder=1)
    def row(items, y, gap):
        total=sum(w for _,_,_,w,_ in items)+gap*(len(items)-1)
        x=CX-total/2; ctr={}
        for label,sub,color,w,fs in items:
            axA.add_patch(FancyBboxPatch((x,y),w,H,boxstyle="round,pad=0.1",facecolor=color,edgecolor='#374151',lw=1.3,alpha=0.93,zorder=2))
            axA.text(x+w/2,y+H/2+0.17,label,ha='center',va='center',fontsize=fs,fontweight='bold',color='white',zorder=3)
            axA.text(x+w/2,y+H/2-0.19,sub,ha='center',va='center',fontsize=6.8,color='white',alpha=0.93,style='italic',zorder=3)
            ctr[label.split('-')[-1]]=x+w/2; x+=w+gap
        return ctr
    axA.text(CX,8.7,'A. Crate architecture — 12-crate Cargo workspace',ha='center',fontsize=13,fontweight='bold',color='#1f2937')
    b=row([('fastvep-cli','annotate · sa-build · filter · cache  ·  76 tests',C['cli'],5.2,11),
           ('fastvep-web','axum production server  ·  9 tests',C['cli'],5.2,11)], tier_y['bin'], 0.9)
    l=row([('fastvep-io','file I/O  ·  86 tests',C['mid'],2.4,9.3),
           ('fastvep-annotate','shared engine  ·  5 tests',C['mid'],2.4,9),
           ('fastvep-consequence','SNV / indel / SV  ·  27 tests',C['mid'],2.7,8.6),
           ('fastvep-hgvs','HGVS g/c/p  ·  25 tests',C['mid'],2.1,9),
           ('fastvep-filter','filter_vep  ·  21 tests',C['mid'],2.0,9)], tier_y['lib'], 0.3)
    p=row([('fastvep-genome','Transcript / Exon / Gene  ·  18 tests',C['data'],3.4,9.3),
           ('fastvep-cache','GFF3 · FASTA mmap · cache  ·  54 tests',C['data'],4.0,9),
           ('fastvep-sa','fastSA · ClinVar, gnomAD  ·  153 tests',C['sa'],3.4,9)], tier_y['prov'], 0.5)
    c=row([('fastvep-core','Consequence (49 SO) · Allele · Impact · VariantType  ·  20 tests',C['core'],7.6,11)], tier_y['core'], 0)
    # dependency arrows (downward, tier to tier)
    yb,yl,yp,yc = tier_y['bin'],tier_y['lib'],tier_y['prov'],tier_y['core']
    arr(b['cli'], yb, l['annotate']-0.25, yl+H); arr(b['web'], yb, l['annotate']+0.25, yl+H)   # binaries -> shared engine
    arr(l['io'], yl, p['genome'], yp+H); arr(l['consequence'], yl, p['cache'], yp+H)
    arr(l['consequence']-0.35, yl, p['genome']+0.5, yp+H)
    for k in ['genome','cache','sa']: arr(p[k], yp, c['core'], yc+H)                            # providers -> core

    # ================= Panel B: annotation data flow =================
    axB.set_xlim(0,10); axB.set_ylim(1.0,9.0); axB.axis('off')
    axB.text(5,8.7,'B. Annotation data flow',ha='center',fontsize=13,fontweight='bold',color='#1f2937')
    stages=[('VCF input','multi-sample · SV · gzip','#64748b'),
            ('Parse & normalize','fastvep-io',C['mid']),
            ('Transcript overlap','fastvep-cache · binary search',C['data']),
            ('Consequence prediction','fastvep-consequence · 49 SO terms',C['mid']),
            ('HGVS g / c / p','fastvep-hgvs',C['mid']),
            ('Supplementary annotation','fastvep-sa · .osa / .osa2',C['sa']),
            ('Output','CSQ (47 fields) · TSV · JSON','#64748b')]
    n=len(stages); top=7.95; bot=1.5; step=(top-bot)/(n-1)
    for i,(t,s,col) in enumerate(stages):
        y=top-i*step
        axB.add_patch(FancyBboxPatch((1.1,y-0.34),7.8,0.68,boxstyle="round,pad=0.1",facecolor=col,edgecolor='#374151',lw=1.2,alpha=0.92))
        axB.text(5,y+0.1,t,ha='center',va='center',fontsize=10,fontweight='bold',color='white')
        axB.text(5,y-0.17,s,ha='center',va='center',fontsize=7.3,color='white',alpha=0.92,style='italic')
        if i<n-1: axB.annotate('',xy=(5,y-step+0.36),xytext=(5,y-0.38),arrowprops=dict(arrowstyle='->',color='#6b7280',lw=1.8))

    # stats strip across the bottom
    stats='12 crates    ·    ~37,500 LOC    ·    641 tests    ·    49 SO terms    ·    47 CSQ fields    ·    3.3 MB static binary    ·    0 runtime dependencies    ·    fastVEP.org'
    fig.text(0.5,0.02,stats,ha='center',fontsize=10.5,color='#374151',
             bbox=dict(boxstyle='round,pad=0.45',facecolor='#f3f4f6',edgecolor='#d1d5db'))
    plt.subplots_adjust(bottom=0.09, top=0.97)
    for ext in ['png','pdf']: plt.savefig(os.path.join(FIG_DIR, f'fig1_architecture.{ext}'), dpi=300, bbox_inches='tight')
    print("Saved fig1_architecture"); plt.close()


def fig2_comparison():
    """THE figure: fastVEP vs Ensembl VEP, single-threaded, measured genome-wide."""
    h2h = read_csv('vep_comparison.csv')       # chr22 scaling (single thread)
    gw = read_csv('vep_genomewide.csv')        # genome-wide 3-condition head-to-head
    orgs = read_csv('organism_comparison.csv')
    if not HAS_MPL: return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.5))

    # === Panel A: Wall-clock time (log-log), single thread ===
    hv = [int(d['variants']) for d in h2h]
    ft = [float(d['fastvep_sec']) for d in h2h]
    vt = [float(d['vep_sec']) for d in h2h]
    ax1.plot(hv, ft, 'o-', color=C['fast'], lw=2.5, ms=9, label='fastVEP (chr22, 1 thread)', zorder=3)
    ax1.plot(hv, vt, 's-', color=C['vep'], lw=2.5, ms=9, label='Ensembl VEP (chr22, 1 thread)', zorder=3)

    # crossover marker
    ax1.axvspan(1000, 5000, color='#94a3b8', alpha=0.08, zorder=0)
    ax1.annotate('crossover\n(~2-3K variants)', xy=(2200, 1.4), fontsize=7.5, color='#475569', ha='center', style='italic')

    # genome-wide MEASURED points at 4.05M (VEP = per-chromosome sum)
    NGW = 4048342
    fv1t = float(gw[0]['fastvep_1t_sec']); fvmt = float(gw[0]['fastvep_mt_sec']); vepgw = float(gw[0]['vep_sec'])
    ax1.plot([hv[-1], NGW], [vt[-1], vepgw], color=C['vep'], lw=1, ls=':', alpha=0.5, zorder=2)
    ax1.plot([hv[-1], NGW], [ft[-1], fv1t], color=C['fast'], lw=1, ls=':', alpha=0.5, zorder=2)
    ax1.scatter([NGW], [vepgw], marker='s', s=150, color=C['vep'], edgecolors='white', lw=1.5, zorder=5)
    ax1.scatter([NGW], [fv1t], marker='o', s=150, color=C['fast'], edgecolors='white', lw=1.5, zorder=5)
    ax1.scatter([NGW], [fvmt], marker='*', s=340, color=C['fast'], edgecolors='white', lw=1.2, zorder=6)
    ax1.annotate('VEP whole genome\n4,621 s (per-chr sum)', xy=(NGW, vepgw), xytext=(12, 2),
                 textcoords='offset points', fontsize=7.5, color=C['vep'], ha='left', va='center', fontweight='bold')
    ax1.annotate('fastVEP 1-thread\n198 s', xy=(NGW, fv1t), xytext=(-12, 2),
                 textcoords='offset points', fontsize=7.5, color=C['fast'], ha='right')
    ax1.annotate('fastVEP multi-thread 93 s', xy=(NGW, fvmt), xytext=(-12, -16),
                 textcoords='offset points', fontsize=7.5, color=C['fast'], ha='right', fontweight='bold')

    # fastVEP multi-organism points
    for d in orgs:
        o = d['organism']; v, t = int(d['variants']), float(d['time_sec'])
        if o == 'human':
            continue  # shown as the 4.05M genome-wide markers above
        ax1.scatter(v, t, c=OC.get(o, '#999'), s=150, marker='D', alpha=0.9, edgecolors='white', lw=1.2, zorder=4)
        label = f"{OL.get(o,o)} {v/1e6:.1f}M" if v > 500000 else f"{OL.get(o,o)} {v/1e3:.0f}K"
        ax1.annotate(label, xy=(v, t), xytext=(7, 5), textcoords='offset points', fontsize=7, color='#333')

    ax1.set_xscale('log'); ax1.set_yscale('log')
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
    ax2.set_title('B. Whole genome (4.05M variants): single-thread head-to-head', fontsize=12.5, fontweight='bold')
    ax2.grid(True, alpha=0.15, axis='y')
    ax2.set_ylim(80, 12000)
    for i in x:
        ax2.text(i - w/2, fv[i]*1.05, f'{fv[i]:.0f}s', ha='center', va='bottom', fontsize=8.5, color=C['fast'], fontweight='bold')
        ax2.text(i + w/2, vp[i]*1.05, f'{vp[i]:.0f}s', ha='center', va='bottom', fontsize=8.5, color=C['vep'], fontweight='bold')
        ax2.text(i, max(vp[i], fv[i])*2.1, f'{sp[i]:.1f}x', ha='center', va='bottom', fontsize=13, fontweight='bold', color='#1f2937')
    ax2.legend(fontsize=9, loc='lower right')
    ax2.annotate('fastVEP multi-thread (10 cores): 93-105 s', xy=(0.02, 0.03), xycoords='axes fraction',
                 fontsize=8, style='italic', color='#475569')

    plt.tight_layout()
    for ext in ['png','pdf']: plt.savefig(os.path.join(FIG_DIR, f'fig2_comparison.{ext}'), dpi=300, bbox_inches='tight')
    print("Saved fig2_comparison"); plt.close()


def fig3_concordance():
    """VEP concordance: field-level 100% + identical consequence calls."""
    data = read_csv('vep_concordance.csv')
    types = read_csv('vep_type_concordance.csv')
    if not HAS_MPL: return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7.5), gridspec_kw={'width_ratios': [1, 1.12]})

    # === Panel A: 23 fields, all 100% concordant (lollipop) ===
    fields = [d['field'] for d in data]
    acc = [float(d['accuracy']) for d in data]
    y = range(len(fields))
    ax1.hlines(list(y), 0, acc, color=C['success'], lw=2.5, alpha=0.55, zorder=1)
    ax1.scatter(acc, list(y), color=C['success'], s=90, zorder=3, edgecolors='white', lw=1)
    for yi in y:
        ax1.text(101.5, yi, '✓', va='center', ha='center', fontsize=11, color=C['success'], fontweight='bold')
    ax1.set_yticks(list(y)); ax1.set_yticklabels(fields, fontsize=9.5, fontfamily='monospace')
    ax1.invert_yaxis()
    ax1.set_xlim(0, 104); ax1.set_xticks([0, 25, 50, 75, 100])
    ax1.axvline(100, color=C['success'], ls='--', alpha=0.4, lw=1)
    ax1.set_xlabel('Concordance with Ensembl VEP v115.1 (%)', fontsize=11)
    ax1.set_title('A. Field-level concordance: 23/23 fields = 100%\n(2,340 shared transcript-allele pairs)', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.2, axis='x')

    # === Panel B: consequence-type calls, fastVEP vs VEP (identical) ===
    tnames = [d['consequence'] for d in types]
    fv = [int(d['fastvep']) for d in types]
    vp = [int(d['vep']) for d in types]
    yb = range(len(tnames)); h = 0.38
    ax2.barh([i + h/2 for i in yb], fv, h, color=C['fast'], alpha=0.9, label='fastVEP')
    ax2.barh([i - h/2 for i in yb], vp, h, color=C['vep'], alpha=0.9, label='Ensembl VEP')
    ax2.set_yticks(list(yb)); ax2.set_yticklabels(tnames, fontsize=9.5, fontfamily='monospace')
    ax2.invert_yaxis()
    ax2.set_xscale('log'); ax2.set_xlim(1, 2000)
    ax2.set_xlabel('Consequence calls (log scale)', fontsize=11)
    ax2.set_title('B. Consequence calls per type: fastVEP vs VEP\n(identical across all 12 observed types)', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.2, axis='x')
    for i in yb:
        ax2.text(fv[i]*1.15, i, f'{fv[i]:,}', va='center', fontsize=8, color='#444')
    ax2.legend(fontsize=10, loc='lower right')

    plt.tight_layout()
    for ext in ['png','pdf']: plt.savefig(os.path.join(FIG_DIR, f'fig3_concordance.{ext}'), dpi=300, bbox_inches='tight')
    print("Saved fig3_concordance"); plt.close()


def fig4_consequences():
    """Consequence distribution (genome-wide, GIAB HG002 full WGS)."""
    data = read_csv('consequence_distribution.csv')
    if not HAS_MPL: return
    consequences = [d['consequence'] for d in data]
    counts = [int(d['count']) for d in data]
    impact_map = {'splice_acceptor_variant':'HIGH','splice_donor_variant':'HIGH','stop_gained':'HIGH',
                  'frameshift_variant':'HIGH','start_lost':'HIGH','stop_lost':'HIGH',
                  'missense_variant':'MODERATE','inframe_insertion':'MODERATE','inframe_deletion':'MODERATE',
                  'splice_region_variant':'LOW','synonymous_variant':'LOW','splice_polypyrimidine_tract_variant':'LOW',
                  'splice_donor_5th_base_variant':'LOW','splice_donor_region_variant':'LOW','stop_retained_variant':'LOW'}
    colors = [C.get(impact_map.get(c,'MODIFIER').lower(), C['modifier']) for c in consequences]
    fig, ax = plt.subplots(figsize=(11, 8))
    bars = ax.barh(range(len(consequences)), counts, color=colors, alpha=0.85)
    ax.set_yticks(range(len(consequences))); ax.set_yticklabels(consequences, fontsize=9.5, fontfamily='monospace')
    ax.invert_yaxis(); ax.set_xlabel('Annotation count (log scale)', fontsize=12)
    ax.set_title('Predicted Consequence Distribution\n(GIAB HG002 full WGS: 4.05M variants, 50.1M annotations, 26 SO terms)', fontsize=13, fontweight='bold')
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

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # === Panel A: throughput vs genome complexity (consistency) ===
    tx = [transcripts[i] for i in idx]; vp = [vps[i] for i in idx]
    cols = [OC.get(orgs[i], '#999') for i in idx]
    ax1.axhspan(46000, 86000, color='#94a3b8', alpha=0.10, zorder=0)
    ax1.plot(tx, vp, '-', color='#cbd5e1', lw=1.5, zorder=1)
    ax1.scatter(tx, vp, c=cols, s=260, edgecolors='white', lw=1.8, zorder=3)
    # per-organism label offsets to avoid collisions (Drosophila/Arabidopsis are adjacent)
    off = {'yeast': (0, 16, 'center'), 'drosophila': (-42, 16, 'center'),
           'arabidopsis': (44, -30, 'center'), 'mouse': (0, 16, 'center'), 'human': (-6, 18, 'right')}
    for i in idx:
        dx, dy, ha = off.get(orgs[i], (0, 16, 'center'))
        ax1.annotate(f"{OL.get(orgs[i],orgs[i])}\n{vps[i]:,} v/s",
                     xy=(transcripts[i], vps[i]), xytext=(dx, dy), textcoords='offset points',
                     ha=ha, fontsize=8.5, fontweight='bold', color='#333')
    ax1.set_xscale('log')
    ax1.set_xlabel('Transcripts (genome complexity, log scale)', fontsize=12)
    ax1.set_ylabel('Throughput (variants/sec)', fontsize=12)
    ax1.set_title('A. Consistent throughput across a 72x range\nof genome complexity (7K-508K transcripts)', fontsize=12, fontweight='bold')
    ax1.set_ylim(0, 100000)
    ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1000:.0f}K'))
    ax1.grid(True, alpha=0.15)
    ax1.text(0.5, 0.06, '47-86K variants/sec band', transform=ax1.transAxes, ha='center',
             fontsize=8.5, style='italic', color='#64748b')

    # === Panel B: absolute scale — complete datasets, wall-clock ===
    yb = range(len(idx))
    nvar = [variants[i] for i in idx]
    bars = ax2.barh(list(yb), nvar, color=cols, alpha=0.88, height=0.62)
    ax2.set_yticks(list(yb))
    ax2.set_yticklabels([f"{OL.get(orgs[i],orgs[i])}" for i in idx], fontsize=11)
    ax2.set_xscale('log'); ax2.set_xlim(1e5, 1e8)
    ax2.set_xlabel('Variants annotated (complete dataset, log scale)', fontsize=12)
    ax2.set_title('B. Complete gold-standard datasets annotated\n(wall-clock time per dataset)', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.15, axis='x')
    for i, ii in enumerate(idx):
        t = times[ii]
        tstr = f"{t:.0f}s" if t < 100 else f"{t/60:.1f} min"
        ax2.text(variants[ii]*1.25, i, f"{variants[ii]/1e6:.2f}M var  ·  {tstr}",
                 va='center', fontsize=9, fontweight='bold', color='#333')

    plt.tight_layout()
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
