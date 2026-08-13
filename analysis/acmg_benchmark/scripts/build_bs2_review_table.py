#!/usr/bin/env python3
"""Build the BS2 gene-disease review table.

BS2 has no curated disease list. It decides per variant from gnomAD counts plus
an inheritance string parsed out of the ClinGen Gene-Disease Validity `.oga`,
and the only hard-coded gene list it consults is the 25-gene homology list that
blocks it entirely. The round-2 medical-genetics review asked for the thing that
does not exist yet: a statement, per gene-disease entity, of whether the
disorder's onset, penetrance and prevalence allow BS1/BS2 to be applied at all.

This script produces the seed for that table. It reads an annotated benchmark
VCF, finds every variant where BS2 fired, aggregates to the gene level, and
ranks genes by how much BS2 actually matters there - misfires on
pathogenic truth first, then raw volume. The output has the evidence columns
filled in and the judgement columns left blank for the reviewer.

Usage:
    build_bs2_review_table.py <annotated.vcf.gz> <out.tsv> [--criterion BS2]
"""

import argparse
import gzip
import sys
from collections import defaultdict

# Columns the reviewer fills in. Order matters: these become the right-hand
# side of the spreadsheet, after the evidence.
REVIEW_COLUMNS = [
    "BS2_applicable__Y_N_or_unsure",
    "BS1_applicable__Y_N_or_unsure",
    "onset__congenital_childhood_adult_late",
    "penetrance__complete_incomplete_unknown",
    "variable_expressivity__Y_N",
    "disease_prevalence__eg_1_in_50000",
    "recommended_BA1_threshold",
    "recommended_BS1_threshold",
    "founder_population__if_any",
    "hypomorph_in_trans_mechanism__Y_N",
    "reviewer_notes",
]

# ClinVar significance -> the three-way truth bucket used for ranking.
PATHOGENIC_TERMS = ("pathogenic",)
BENIGN_TERMS = ("benign",)


def truth_bucket(clnsig):
    """Collapse a CLNSIG string to pathogenic / benign / uncertain.

    Checked in this order because `Likely_pathogenic` contains "pathogenic"
    and `Conflicting...` contains neither.
    """
    s = clnsig.lower()
    if "conflicting" in s:
        return "uncertain"
    has_p = any(t in s for t in PATHOGENIC_TERMS)
    has_b = any(t in s for t in BENIGN_TERMS)
    if has_p and not has_b:
        return "pathogenic"
    if has_b and not has_p:
        return "benign"
    return "uncertain"


def parse_info(info):
    out = {}
    for part in info.split(";"):
        k, _, v = part.partition("=")
        out[k] = v
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf")
    ap.add_argument("out")
    ap.add_argument("--criterion", default="BS2",
                    help="criterion code to aggregate (default BS2)")
    ap.add_argument("--max-examples", type=int, default=3)
    ap.add_argument("--misfire-vcf",
                    help="also write the variants where the criterion fired on "
                         "pathogenic truth, as a VCF. Re-annotate it to recover "
                         "the gnomAD counts behind each call.")
    args = ap.parse_args()
    misfires = []

    csq_fields = None
    genes = defaultdict(lambda: {
        "fired": 0,
        "fired_pathogenic": 0,
        "fired_benign": 0,
        "fired_uncertain": 0,
        "diseases": defaultdict(int),
        "examples": [],
        "cosignatures": defaultdict(int),
    })

    opener = gzip.open if args.vcf.endswith(".gz") else open
    with opener(args.vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                if line.startswith("##INFO=<ID=CSQ"):
                    fmt = line.split("Format: ")[1].rstrip('">\n')
                    csq_fields = fmt.split("|")
                continue
            if line.startswith("#"):
                continue
            if csq_fields is None:
                sys.exit("No CSQ header found in the VCF")

            cols = line.rstrip("\n").split("\t")
            info = parse_info(cols[7])
            csq = info.get("CSQ")
            if not csq:
                continue

            # `--pick` leaves one consequence block per variant, but tolerate
            # multiple by taking the first that carries the criterion.
            for block in csq.split(","):
                vals = block.split("|")
                if len(vals) != len(csq_fields):
                    continue
                rec = dict(zip(csq_fields, vals))
                criteria = rec.get("ACMG_CRITERIA", "")
                codes = [c for c in criteria.split("&") if c]
                if args.criterion not in codes:
                    continue

                symbol = rec.get("SYMBOL") or "(no gene)"
                g = genes[symbol]
                g["fired"] += 1
                bucket = truth_bucket(info.get("CLNSIG", ""))
                g["fired_" + bucket] += 1
                for dn in info.get("CLNDN", "").split("|"):
                    if dn and dn not in ("not_provided", "not_specified"):
                        g["diseases"][dn.replace("_", " ")] += 1
                g["cosignatures"]["+".join(sorted(codes))] += 1

                # Keep the harmful cases as examples: the criterion fired on
                # pathogenic truth.
                if bucket == "pathogenic":
                    if len(g["examples"]) < args.max_examples:
                        g["examples"].append(
                            f"{cols[0]}:{cols[1]}{cols[3]}>{cols[4]} "
                            f"({info.get('CLNSIG', '?')}; called {rec.get('ACMG', '?')})"
                        )
                    misfires.append({
                        "chrom": cols[0], "pos": cols[1],
                        "ref": cols[3], "alt": cols[4],
                        "gene": symbol,
                        "clnsig": info.get("CLNSIG", ""),
                        "clndn": info.get("CLNDN", "").replace("_", " "),
                        "hgvs": info.get("CLNHGVS", ""),
                        "called": rec.get("ACMG", ""),
                        "criteria": criteria.replace("&", "+"),
                        "hgvsc": rec.get("HGVSc", ""),
                        "hgvsp": rec.get("HGVSp", ""),
                        "consequence": rec.get("Consequence", ""),
                    })
                break

    evidence_columns = [
        "gene",
        "n_variants_with_" + args.criterion,
        "n_on_pathogenic_truth",
        "n_on_benign_truth",
        "n_on_uncertain_truth",
        "clinvar_diseases_top3",
        "example_misfires",
        "most_common_criteria_signature",
    ]

    rows = []
    for symbol, g in genes.items():
        top_diseases = sorted(g["diseases"].items(), key=lambda kv: -kv[1])[:3]
        top_sig = max(g["cosignatures"].items(), key=lambda kv: kv[1])[0] if g["cosignatures"] else ""
        rows.append({
            "gene": symbol,
            "n_variants_with_" + args.criterion: g["fired"],
            "n_on_pathogenic_truth": g["fired_pathogenic"],
            "n_on_benign_truth": g["fired_benign"],
            "n_on_uncertain_truth": g["fired_uncertain"],
            "clinvar_diseases_top3": "; ".join(d for d, _ in top_diseases),
            "example_misfires": " | ".join(g["examples"]),
            "most_common_criteria_signature": top_sig,
        })

    # Misfires first (those are where the criterion is actively harmful), then
    # raw volume. A gene where BS2 fires 4,000 times and is never wrong needs
    # far less of the reviewer's attention than one where it is wrong twice.
    rows.sort(key=lambda r: (-r["n_on_pathogenic_truth"],
                             -r["n_variants_with_" + args.criterion]))

    header = evidence_columns + REVIEW_COLUMNS
    with open(args.out, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "")) for c in header) + "\n")

    if args.misfire_vcf:
        with open(args.misfire_vcf, "w") as vf:
            vf.write("##fileformat=VCFv4.2\n")
            vf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for m in misfires:
                vf.write(f"{m['chrom']}\t{m['pos']}\t.\t{m['ref']}\t{m['alt']}"
                         f"\t.\t.\t.\n")
        # The per-variant detail rides alongside as a TSV keyed the same way.
        with open(args.misfire_vcf + ".detail.tsv", "w") as df:
            cols_out = ["chrom", "pos", "ref", "alt", "gene", "clnsig", "clndn",
                        "hgvs", "hgvsc", "hgvsp", "consequence", "called", "criteria"]
            df.write("\t".join(cols_out) + "\n")
            for m in misfires:
                df.write("\t".join(str(m.get(c, "")) for c in cols_out) + "\n")
        print(f"wrote {args.misfire_vcf} ({len(misfires)} variants) "
              f"and {args.misfire_vcf}.detail.tsv")

    n_misfire_genes = sum(1 for r in rows if r["n_on_pathogenic_truth"] > 0)
    print(f"{len(rows)} genes with {args.criterion}; "
          f"{n_misfire_genes} have at least one firing on pathogenic truth")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
