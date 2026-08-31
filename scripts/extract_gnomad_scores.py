#!/usr/bin/env python3
"""Distil SpliceAI and PhyloP databases out of the gnomAD v4.1 sites VCF.

Both scores are already in gnomAD's INFO column, which removes the two access
barriers in the ACMG stack: Illumina's precomputed SpliceAI is behind a
BaseSpace account, and UCSC's phyloP ships as 25 per-chromosome wigFix files
totalling about 5 GB. gnomAD carries `spliceai_ds_max` and `phylop` on records
it already publishes, so a run that fetches gnomAD for allele frequencies can
produce all three databases from one pass.

What is lost, precisely:

  * gnomAD stores only SpliceAI's **maximum** delta score, not the four
    acceptor-gain / acceptor-loss / donor-gain / donor-loss components. This
    writes the max into all four DS fields. Every ACMG consumer of SpliceAI in
    fastVEP reads `max_delta_score()` and nothing else (PVS1's splice branch,
    PP3, BP4, BP7), so the classification is unaffected. The DP position fields
    have no gnomAD equivalent and are written as 0, which means the `SpliceAI`
    INFO field fastVEP emits carries placeholder positions. Do not read those
    positions downstream.

  * gnomAD's `phylop` is the Zoonomia 241-placental-mammal score, not UCSC's
    hg38 100-way vertebrate phyloP. They are different scores on different
    scales. fastVEP's `phylop_conserved` default of 2.0 was measured against
    this one, so this is the score the shipped threshold expects.

  * Coverage is gnomAD's coverage. The exomes release covers capture targets
    and their flanks, so a deep-intronic variant has no record and BP7's
    deep-intronic extension has nothing to read there. Use the genomes release
    for whole-genome coverage, at roughly 2.8x the download.

Usage:
    # One chromosome, restricted to regions you care about
    python3 scripts/extract_gnomad_scores.py --chrom 17 --regions targets.bed \
        --out-dir extracts/

    # Then build, per chromosome
    fastvep sa-build --source spliceai -i extracts/spliceai_chr17.vcf \
        -o sa_databases/spliceai_chr17 --assembly GRCh38
    fastvep sa-build --source phylop -i extracts/phylop_chr17.tsv \
        -o sa_databases/phylop_chr17 --assembly GRCh38

Requires `tabix` on PATH. It reads the public bucket over https, so nothing is
downloaded whole.
"""

import argparse
import os
import re
import subprocess
import sys

URL = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/"
       "vcf/{kind}/gnomad.{kind}.v4.1.sites.chr{chrom}.vcf.bgz")

SPLICE_RE = re.compile(r"(?:^|;)spliceai_ds_max=([^;]+)")
PHYLOP_RE = re.compile(r"(?:^|;)phylop=([^;]+)")

# tabix takes regions as argv, so they are fed in batches to stay under the
# platform's argument-length limit. 200 matches the benchmark's build scripts.
REGION_BATCH = 200


def load_regions(path, chrom):
    """BED rows for this chromosome as 1-based inclusive tabix region strings.

    Accepts both `chr17` and `17` in the BED, because a target file written
    against a UCSC reference and one written against an Ensembl reference are
    otherwise identical and it is not worth making the user convert.
    """
    wanted = {f"chr{chrom}", chrom}
    regions = []
    with open(path) as f:
        for line in f:
            if line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] not in wanted:
                continue
            # gnomAD's VCFs are `chr`-prefixed; normalise regardless of input.
            regions.append(f"chr{chrom}:{int(parts[1]) + 1}-{parts[2]}")
    return regions


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chrom", required=True,
                   help="chromosome without the `chr` prefix, e.g. 17 or X")
    p.add_argument("--regions",
                   help="BED of regions to extract. Omitting it reads the "
                        "whole chromosome, which streams the entire file.")
    p.add_argument("--kind", choices=["exomes", "genomes"], default="exomes",
                   help="gnomAD release to read (default: exomes)")
    p.add_argument("--out-dir", default=".")
    p.add_argument("--url", help="override the source URL, e.g. a local .bgz")
    args = p.parse_args()

    url = args.url or URL.format(kind=args.kind, chrom=args.chrom)
    os.makedirs(args.out_dir, exist_ok=True)
    sa_path = os.path.join(args.out_dir, f"spliceai_chr{args.chrom}.vcf")
    pp_path = os.path.join(args.out_dir, f"phylop_chr{args.chrom}.tsv")

    if args.regions:
        regions = load_regions(args.regions, args.chrom)
        if not regions:
            sys.exit(f"no rows for chromosome {args.chrom} in {args.regions}")
        batches = [regions[i:i + REGION_BATCH]
                   for i in range(0, len(regions), REGION_BATCH)]
    else:
        batches = [[f"chr{args.chrom}"]]

    sa_n = pp_n = 0
    with open(sa_path, "w") as sa_out, open(pp_path, "w") as pp_out:
        sa_out.write("##fileformat=VCFv4.2\n")
        sa_out.write('##INFO=<ID=SpliceAI,Number=.,Type=String,'
                     'Description="SpliceAI max delta score from gnomAD '
                     'spliceai_ds_max; the four DS fields carry the same max '
                     'and the DP fields are placeholders">\n')
        sa_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for batch in batches:
            proc = subprocess.Popen(["tabix", url, *batch],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, text=True)
            for line in proc.stdout:
                if not line or line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 8:
                    continue
                chrom, pos, ref, alts, info = (cols[0], cols[1], cols[3],
                                               cols[4], cols[7])

                m = SPLICE_RE.search(info)
                if m and m.group(1) not in (".", ""):
                    ds = m.group(1)
                    for alt in alts.split(","):
                        if alt in (".", "*"):
                            continue
                        sa_out.write(
                            f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t"
                            f"SpliceAI={alt}|gnomad|{ds}|{ds}|{ds}|{ds}|0|0|0|0\n")
                        sa_n += 1

                # PhyloP is per position, not per allele, so one row per record.
                q = PHYLOP_RE.search(info)
                if q and q.group(1) not in (".", ""):
                    pp_out.write(f"{chrom}\t{pos}\t{q.group(1)}\n")
                    pp_n += 1
            proc.wait()
            if proc.returncode != 0:
                sys.exit("tabix failed: " + (proc.stderr.read() or "").strip())

    print(f"chr{args.chrom}: spliceai={sa_n:,} -> {sa_path}", file=sys.stderr)
    print(f"chr{args.chrom}: phylop={pp_n:,} -> {pp_path}", file=sys.stderr)
    if sa_n == 0 and pp_n == 0:
        sys.exit("extracted nothing; check the region file and the chromosome")


if __name__ == "__main__":
    main()
