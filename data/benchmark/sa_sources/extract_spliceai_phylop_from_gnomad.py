#!/usr/bin/env python3
"""Extract SpliceAI + PhyloP from gnomAD v4.1 exomes VCF.

gnomAD v4.1 sites VCFs already include `spliceai_ds_max` and `phylop`
in the INFO field, so we don't need to download separate SpliceAI
(~32 GB, BaseSpace-gated) or PhyloP (~50 GB UCSC bigwig) databases —
we just re-pull the same regions we used for the gnomAD .osa build
and emit two side-files per chromosome:

  spliceai_chr{N}.vcf      one row per allele, INFO=SpliceAI=A|GENE|<ds>|<ds>|<ds>|<ds>|.|.|.|.
                           where <ds> is the max delta score (gnomAD only stores the max,
                           so AG/AL/DG/DL are filled with the same value — lossy but
                           correct for threshold-based criteria like PP3 max_ds ≥ 0.2)
  phylop_chr{N}.tsv        chrom\tpos\tphylop_value (3-column simple format)

Usage:
    extract_spliceai_phylop_from_gnomad.py <chrom>
"""

import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BED = f"{ROOT}/data/benchmark/sa_sources/clinvar_2star_regions.bed"
OUT = f"{ROOT}/data/benchmark/sa_sources/spliceai_phylop_extracts"
URL = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr{}.vcf.bgz"


SPLICE_RE = re.compile(r"spliceai_ds_max=([^;]+)")
PHYLOP_RE = re.compile(r"phylop=([^;]+)")


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: extract_spliceai_phylop_from_gnomad.py <chrom>")
    chrom = sys.argv[1]
    chr_label = f"chr{chrom}"

    regions = []
    with open(BED) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[0] == chr_label:
                regions.append(f"{parts[0]}:{int(parts[1]) + 1}-{parts[2]}")
    if not regions:
        sys.exit(f"no regions for {chr_label} in {BED}")

    os.makedirs(OUT, exist_ok=True)
    spliceai_path = f"{OUT}/spliceai_chr{chrom}.vcf"
    phylop_path = f"{OUT}/phylop_chr{chrom}.tsv"

    sa_n = 0
    pp_n = 0
    with open(spliceai_path, "w") as sa_out, open(phylop_path, "w") as pp_out:
        sa_out.write("##fileformat=VCFv4.2\n")
        sa_out.write(
            "##INFO=<ID=SpliceAI,Number=.,Type=String,Description=\"SpliceAI score (synthesized from gnomAD spliceai_ds_max)\">\n"
        )
        sa_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        # tabix in batches to bound argv size; -n 200 mirrors the gnomad build script.
        # Run a single tabix with all regions piped via xargs for parallelism.
        url = URL.format(chrom)
        # Use one tabix call; tabix accepts many region args
        for batch_start in range(0, len(regions), 200):
            batch = regions[batch_start : batch_start + 200]
            proc = subprocess.Popen(
                ["tabix", url, *batch],
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                text=True,
            )
            assert proc.stdout is not None
            for line in proc.stdout:
                if not line or line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 8:
                    continue
                vcf_chrom = cols[0]
                pos = cols[1]
                ref = cols[3]
                alts = cols[4]
                info = cols[7]

                m = SPLICE_RE.search(info)
                if m and m.group(1) not in (".", ""):
                    ds = m.group(1)
                    for alt in alts.split(","):
                        if alt in (".", "*"):
                            continue
                        # SpliceAI parser expects A|GENE|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL.
                        # gnomAD only stores the max — duplicate it across all four DS fields.
                        # DP fields are int positions; the parser rejects "." here, so we
                        # write 0 as a placeholder (criteria use DS thresholds only).
                        sa_info = f"SpliceAI={alt}|gnomad|{ds}|{ds}|{ds}|{ds}|0|0|0|0"
                        sa_out.write(
                            f"{vcf_chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{sa_info}\n"
                        )
                        sa_n += 1

                p = PHYLOP_RE.search(info)
                if p and p.group(1) not in (".", ""):
                    pp_out.write(f"{vcf_chrom}\t{pos}\t{p.group(1)}\n")
                    pp_n += 1
            proc.wait()

    print(
        f"[chr{chrom}] spliceai={sa_n:,} phylop={pp_n:,} → {spliceai_path}, {phylop_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
