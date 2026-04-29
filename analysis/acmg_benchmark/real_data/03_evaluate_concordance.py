#!/usr/bin/env python3
"""
Evaluate ACMG classifier concordance against ClinVar 2-star+ truth.

Reads:
  truth TSV   — chrom\tpos\tref\talt\tgene\tclnsig\tnormalized_class\treview_stars\trcv
  fastvep output (JSON array, --output-format json)

Joins on (chrom, pos, ref, alt) and emits to ./output/:
  concordance_matrix.csv
  concordance_summary.txt
  discrepancies.tsv
  criterion_firing_rates.csv
"""

import argparse
import json
import csv
import sys
from pathlib import Path
from collections import defaultdict, Counter

CLASSES = ["Pathogenic", "Likely_pathogenic", "VUS", "Likely_benign", "Benign"]


def load_truth(path: Path) -> dict:
    truth = {}
    with path.open() as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for row in rdr:
            key = (row["chrom"], row["pos"], row["ref"], row["alt"])
            truth[key] = row
    return truth


def variant_key(rec: dict):
    chrom = str(rec.get("seq_region_name") or rec.get("chromosome") or rec.get("chrom") or "")
    pos = str(rec.get("start") or rec.get("position") or rec.get("pos") or "")
    allele = rec.get("allele_string", "")
    if "/" in allele:
        ref, alts = allele.split("/", 1)
        alt_list = alts.split(",")
        return [(chrom, pos, ref, alt) for alt in alt_list]
    return [(chrom, pos, "", "")]


def class_label(c):
    if c in ("Pathogenic",):
        return "Pathogenic"
    if c in ("Likely_pathogenic", "LikelyPathogenic"):
        return "Likely_pathogenic"
    if c in ("Uncertain_significance", "UncertainSignificance", "VUS"):
        return "VUS"
    if c in ("Likely_benign", "LikelyBenign"):
        return "Likely_benign"
    if c in ("Benign",):
        return "Benign"
    return None


def get_top_consequence(tc):
    cs = tc.get("consequence_terms", [])
    if cs:
        return cs[0]
    return "unknown"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--truth", required=True)
    p.add_argument("--predictions", required=True, help="fastvep JSON output (array)")
    p.add_argument("--out", default="output")
    args = p.parse_args()

    truth = load_truth(Path(args.truth))
    print(f"Loaded {len(truth):,} truth records")

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    cm = {tc: {pc: 0 for pc in CLASSES + ["NoCall"]} for tc in CLASSES}
    cm_consequence = defaultdict(lambda: {tc: {pc: 0 for pc in CLASSES + ["NoCall"]} for tc in CLASSES})
    criterion_fires = {tc: Counter() for tc in CLASSES}
    n_classified = 0
    rule_dist = Counter()
    discrepancies = []
    matched_keys = set()

    # Stream JSON array
    print(f"Streaming predictions from {args.predictions}...")
    with open(args.predictions) as f:
        # Try JSON array
        first_byte = f.read(1)
        if first_byte == "[":
            f.seek(0)
            recs = json.load(f)
        else:
            # JSONL
            f.seek(0)
            recs = (json.loads(line) for line in f if line.strip())
        for rec in recs:
            for k in variant_key(rec):
                if k not in truth:
                    continue
                t = truth[k]
                tc_truth = t["normalized_class"]
                # Take the first transcript_consequence with an acmg block
                acmg_block = None
                top_csq = "unknown"
                for tc in rec.get("transcript_consequences", []) or []:
                    if "acmg" in tc:
                        acmg_block = tc["acmg"]
                        top_csq = get_top_consequence(tc)
                        break
                if acmg_block is None:
                    cm[tc_truth]["NoCall"] += 1
                    continue
                pc = class_label(acmg_block.get("classification")) or "NoCall"
                cm[tc_truth][pc] += 1
                cm_consequence[top_csq][tc_truth][pc] += 1
                rule = acmg_block.get("triggered_rule") or ""
                if rule:
                    rule_dist[rule] += 1
                for c in acmg_block.get("criteria", []) or []:
                    if c.get("met"):
                        criterion_fires[tc_truth][c["code"]] += 1
                n_classified += 1
                matched_keys.add(k)
                # Track discrepancies (opposite-direction or pathogenic→VUS/NoCall)
                if (tc_truth in ("Pathogenic", "Likely_pathogenic") and pc in ("Benign", "Likely_benign", "VUS", "NoCall")) \
                        or (tc_truth in ("Benign", "Likely_benign") and pc in ("Pathogenic", "Likely_pathogenic")):
                    if len(discrepancies) < 10000:
                        met_codes = ";".join(c["code"] for c in (acmg_block.get("criteria") or []) if c.get("met"))
                        discrepancies.append((k[0], k[1], k[2], k[3], t["gene"], t["review_stars"], tc_truth, pc, top_csq, rule, met_codes))
                break

    # Concordance matrix
    matrix_path = out_dir / "concordance_matrix.csv"
    with matrix_path.open("w") as f:
        w = csv.writer(f)
        w.writerow(["truth"] + CLASSES + ["NoCall"])
        for tc in CLASSES:
            w.writerow([tc] + [cm[tc][pc] for pc in CLASSES + ["NoCall"]])

    # By-consequence (top 12 most frequent consequences)
    consq_counts = sorted(cm_consequence.items(), key=lambda kv: -sum(sum(d.values()) for d in kv[1].values()))[:12]
    by_csq_path = out_dir / "concordance_by_consequence.csv"
    with by_csq_path.open("w") as f:
        w = csv.writer(f)
        w.writerow(["consequence", "truth"] + CLASSES + ["NoCall", "n"])
        for csq, mat in consq_counts:
            for tcl in CLASSES:
                row = mat[tcl]
                n = sum(row.values())
                if n == 0:
                    continue
                w.writerow([csq, tcl] + [row[pc] for pc in CLASSES + ["NoCall"]] + [n])

    # Summary text
    summary_path = out_dir / "concordance_summary.txt"
    totals = {"n": 0, "exact": 0, "same_dir": 0, "opp": 0, "no_call": 0}
    with summary_path.open("w") as f:
        f.write("ClinVar 2-star+ concordance against fastvep ACMG classifier (real data)\n")
        f.write("=" * 75 + "\n\n")
        f.write(f"Truth records:       {len(truth):,}\n")
        f.write(f"Classified:          {n_classified:,}\n")
        f.write(f"Truth not annotated: {len(truth) - len(matched_keys):,}\n\n")
        f.write("Per-class breakdown:\n")
        f.write(f"{'truth':<22} {'n':>8} {'exact':>8} {'same_dir':>10} {'opposite':>10} {'no_call':>8}\n")
        for tcl in CLASSES:
            row = cm[tcl]
            n = sum(row.values())
            exact = row[tcl]
            if tcl in ("Pathogenic", "Likely_pathogenic"):
                same_dir = row["Pathogenic"] + row["Likely_pathogenic"]
                opp = row["Benign"] + row["Likely_benign"]
            elif tcl == "VUS":
                same_dir = row["VUS"]
                opp = 0
            else:
                same_dir = row["Benign"] + row["Likely_benign"]
                opp = row["Pathogenic"] + row["Likely_pathogenic"]
            no_call = row["NoCall"]
            f.write(f"{tcl:<22} {n:>8} {exact:>8} {same_dir:>10} {opp:>10} {no_call:>8}\n")
            for k, v in [("n", n), ("exact", exact), ("same_dir", same_dir), ("opp", opp), ("no_call", no_call)]:
                totals[k] += v
        f.write(f"\n{'TOTAL':<22} {totals['n']:>8} {totals['exact']:>8} {totals['same_dir']:>10} {totals['opp']:>10} {totals['no_call']:>8}\n")
        if totals["n"]:
            f.write(f"\nExact-match rate:        {totals['exact']/totals['n']*100:.1f}%\n")
            f.write(f"Same-direction rate:     {totals['same_dir']/totals['n']*100:.1f}%\n")
            f.write(f"Opposite-direction rate: {totals['opp']/totals['n']*100:.1f}%\n")
            f.write(f"NoCall rate:             {totals['no_call']/totals['n']*100:.1f}%\n")

        f.write("\nTop 20 triggered rules:\n")
        for rule, n in rule_dist.most_common(20):
            f.write(f"  {n:>8} {rule}\n")

        f.write("\nCriterion firing rates (% of variants in each truth class where the criterion fired):\n")
        all_codes = sorted({c for cnt in criterion_fires.values() for c in cnt})
        f.write(f"{'criterion':<22}")
        for tcl in CLASSES:
            f.write(f"{tcl:>20}")
        f.write("\n")
        for c in all_codes:
            f.write(f"{c:<22}")
            for tcl in CLASSES:
                tot = sum(cm[tcl].values())
                fired = criterion_fires[tcl].get(c, 0)
                pct = (fired / tot * 100) if tot else 0
                f.write(f"  {fired:>6} ({pct:>5.1f}%) ")
            f.write("\n")

    # Discrepancies TSV
    disc_path = out_dir / "discrepancies.tsv"
    with disc_path.open("w") as f:
        f.write("chrom\tpos\tref\talt\tgene\tstars\ttruth\tpredicted\tconsequence\trule\tmet_criteria\n")
        for d in discrepancies:
            f.write("\t".join(str(x) for x in d) + "\n")

    print("\nOutputs:")
    print(f"  {matrix_path}")
    print(f"  {by_csq_path}")
    print(f"  {summary_path}")
    print(f"  {disc_path}")
    print()
    print(open(summary_path).read())


if __name__ == "__main__":
    main()
