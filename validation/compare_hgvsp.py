#!/usr/bin/env python3
"""Compare fastVEP and Ensembl VEP protein-level HGVS descriptions.

Usage:
    python compare_hgvsp.py <fastvep.vcf> <vep.vcf> [--verbose]

The general field-level comparison in compare_vep.py runs on the VEP example
dataset, which is entirely SNVs. Protein-level 3'-normalisation only engages on
insertions and deletions, so it needs an indel-bearing input to be exercised at
all; validation/human/clinvar_inframe_deletions.vcf is one.

Reports exact-string agreement of HGVSp over every (variant, transcript) pair
both tools annotate, and classifies each disagreement by the shift distance
between the two descriptions.
"""
import argparse
import re
import sys
from collections import Counter

HGVSP = re.compile(r'^p\.([A-Z][a-z]{2})(\d+)(?:_([A-Z][a-z]{2})(\d+))?(del|dup|delins|ins|=)?')


def load(path):
    """Return {(chrom, pos, ref, alt, feature): csq_dict} for rows carrying HGVSp."""
    header, rows = None, {}
    with open(path) as handle:
        for line in handle:
            if line.startswith('##INFO=<ID=CSQ'):
                header = re.search(r'Format: ([^"]+)', line).group(1).split('|')
                continue
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            found = re.search(r'CSQ=([^;\t]*)', fields[7])
            if not found:
                continue
            for entry in found.group(1).split(','):
                csq = dict(zip(header, entry.split('|')))
                csq['HGVSp'] = csq.get('HGVSp', '').replace('%3D', '=')
                if csq['HGVSp']:
                    key = (fields[0], fields[1], fields[3], fields[4], csq.get('Feature', ''))
                    rows[key] = csq
    return rows


def span(hgvsp):
    """(start, end, kind) of a protein description, or None if it is not a range."""
    m = HGVSP.match(hgvsp.split(':')[-1])
    if not m:
        return None
    start = int(m.group(2))
    return start, int(m.group(4)) if m.group(4) else start, m.group(5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('fastvep')
    ap.add_argument('vep')
    ap.add_argument('--verbose', action='store_true')
    args = ap.parse_args()

    ours, theirs = load(args.fastvep), load(args.vep)
    shared = [k for k in theirs if k in ours]
    if not shared:
        sys.exit('no (variant, transcript) pair carries HGVSp in both inputs')

    differ = [k for k in shared if ours[k]['HGVSp'] != theirs[k]['HGVSp']]
    agree = len(shared) - len(differ)
    print(f'HGVSp pairs compared : {len(shared)}')
    print(f'exact agreement      : {agree} ({100 * agree / len(shared):.2f}%)')
    print(f'disagreements        : {len(differ)}')

    sites = {k[:4] for k in shared}
    bad_sites = {k[:4] for k in differ}
    print(f'variant sites        : {len(sites)}, of which {len(bad_sites)} diverge on any transcript')
    if not differ:
        return

    # Ensembl VEP's _shift_3prime scans length(post_seq) - deleted_length
    # positions while consuming one residue per iteration, so for an n-residue
    # change it stops up to n-1 residues short of the protein terminus. Any
    # disagreement outside that envelope is a different phenomenon.
    classes = Counter()
    for key in differ:
        a, b = span(ours[key]['HGVSp']), span(theirs[key]['HGVSp'])
        if a and b and a[2] == b[2] and (a[1] - a[0]) == (b[1] - b[0]):
            n, shift = a[1] - a[0] + 1, a[0] - b[0]
            known = 1 <= shift <= n - 1
            classes[f'{a[2]}: we shift {shift} further, n={n} '
                    f'({"C-terminal shift bound" if known else "UNEXPLAINED"})'] += 1
        else:
            classes['different description shape (UNEXPLAINED)'] += 1
    for label, count in classes.most_common():
        print(f'  {count:5}  {label}')

    unexplained = sum(c for lbl, c in classes.items() if 'UNEXPLAINED' in lbl)
    if args.verbose:
        for key in differ[:20]:
            print(f'    {ours[key].get("SYMBOL", ""):10} {key[4]:18} '
                  f'ours={ours[key]["HGVSp"].split(":")[-1]:24} '
                  f'vep={theirs[key]["HGVSp"].split(":")[-1]}')
    sys.exit(1 if unexplained else 0)


if __name__ == '__main__':
    main()
