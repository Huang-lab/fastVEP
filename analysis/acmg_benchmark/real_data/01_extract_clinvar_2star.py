#!/usr/bin/env python3
"""
Extract ClinVar 2-star+ variants from clinvar.vcf.gz into a benchmark VCF
plus a TSV of (chrom, pos, ref, alt, clnsig, gene) for evaluation.

ClinVar 2-star+ means CLNREVSTAT contains:
- criteria_provided,_multiple_submitters,_no_conflicts (2 stars)
- reviewed_by_expert_panel (3 stars)
- practice_guideline (4 stars)
"""

import gzip
import re
import sys
from pathlib import Path

CLNVAR_VCF = Path(sys.argv[1] if len(sys.argv) > 1 else
                  "/Users/kuan-lin.huang/Projects/fastVEP/genomes/human/sa_sources/clinvar.vcf.gz")
OUT_DIR = Path(sys.argv[2] if len(sys.argv) > 2 else
               "/Users/kuan-lin.huang/Projects/fastVEP/analysis/acmg_benchmark/real_data")
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_VCF = OUT_DIR / "clinvar_2star.vcf"
OUT_TSV = OUT_DIR / "clinvar_2star_truth.tsv"

# ClinVar review status → star level
def review_stars(clnrevstat: str) -> int:
    s = clnrevstat.lower()
    if "practice_guideline" in s:
        return 4
    if "reviewed_by_expert_panel" in s:
        return 3
    if "criteria_provided,_multiple_submitters,_no_conflicts" in s:
        return 2
    if "criteria_provided,_conflicting_classifications" in s:
        return 1
    if "criteria_provided,_single_submitter" in s:
        return 1
    if "no_assertion_criteria_provided" in s:
        return 0
    if "no_assertion_provided" in s:
        return 0
    return 0


# Map ClinVar CLNSIG to our 5-tier classification.
def normalize_clnsig(clnsig: str):
    s = clnsig.lower().replace(",_", ",").replace(",", " ")
    if "pathogenic/likely_pathogenic" in clnsig.lower():
        return "Pathogenic"
    if "benign/likely_benign" in clnsig.lower():
        return "Benign"
    if "likely_pathogenic" in clnsig.lower():
        return "Likely_pathogenic"
    if "likely_benign" in clnsig.lower():
        return "Likely_benign"
    if clnsig.lower().startswith("pathogenic"):
        return "Pathogenic"
    if clnsig.lower().startswith("benign"):
        return "Benign"
    if "uncertain_significance" in clnsig.lower():
        return "VUS"
    if "conflicting_classifications" in clnsig.lower():
        return "VUS"
    if "conflicting_interpretations" in clnsig.lower():
        return "VUS"
    return None


def parse_info(info: str) -> dict:
    out = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


count_total = 0
count_kept = 0
counts_by_class: dict = {}
counts_by_stars: dict = {}

with gzip.open(CLNVAR_VCF, "rt") as fh, OUT_VCF.open("w") as vcf_out, OUT_TSV.open("w") as tsv_out:
    tsv_out.write("chrom\tpos\tref\talt\tgene\tclnsig\tnormalized_class\treview_stars\trcv\n")
    for line in fh:
        if line.startswith("#"):
            vcf_out.write(line)
            continue
        count_total += 1
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            continue
        chrom, pos, vid, ref, alt, qual, flt, info = fields[:8]
        info_d = parse_info(info)
        clnrevstat = info_d.get("CLNREVSTAT", "")
        clnsig = info_d.get("CLNSIG", "")
        gene = info_d.get("GENEINFO", "").split(":")[0]
        rcv = info_d.get("CLNVI", "") or vid
        stars = review_stars(clnrevstat)
        if stars < 2:
            continue
        norm = normalize_clnsig(clnsig)
        if norm is None:
            continue
        # Skip MNPs / structural variants, keep SNVs and small indels
        if "*" in alt or len(alt) > 100 or len(ref) > 100:
            continue
        # Add chr prefix to match Ensembl FASTA contigs (1, 2, ..., X, Y)
        # — Ensembl uses unprefixed; ClinVar VCF already uses unprefixed. Leave alone.
        vcf_out.write(line)
        tsv_out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t{clnsig}\t{norm}\t{stars}\t{rcv}\n")
        counts_by_class[norm] = counts_by_class.get(norm, 0) + 1
        counts_by_stars[stars] = counts_by_stars.get(stars, 0) + 1
        count_kept += 1

print(f"Total ClinVar records:    {count_total:>10,}")
print(f"Kept (2-star+, classifiable, SNV/small indel): {count_kept:>10,}")
print(f"By normalized class:")
for k, v in sorted(counts_by_class.items(), key=lambda kv: -kv[1]):
    print(f"  {k:<20}{v:>10,}")
print(f"By star level:")
for k, v in sorted(counts_by_stars.items()):
    print(f"  {k}-star{'':<14}{v:>10,}")
print(f"\nOutputs:\n  {OUT_VCF}\n  {OUT_TSV}")
