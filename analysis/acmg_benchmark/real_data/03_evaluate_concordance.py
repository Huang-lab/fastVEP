#!/usr/bin/env python3
"""
Evaluate ACMG classifier concordance against ClinVar 2-star+ truth.

Inputs
  --truth        TSV: chrom\tpos\tref\talt\tgene\tclnsig\tnormalized_class\treview_stars\trcv
  --predictions  fastvep --output-format vcf (.vcf.gz / .vcf / .vcf.bgz)

The VCF stores annotations in INFO/CSQ as a comma-separated list of
pipe-separated transcript entries. The CSQ format header line carries
the field order. We pick the first CSQ entry whose ACMG field is
non-empty (with `--pick` this is usually the canonical transcript).

Outputs (under --out)
  concordance_matrix.csv             5×6 (truth × predicted+NoCall)
  concordance_by_chrom.csv           per-chromosome concordance
  concordance_by_consequence.csv     top consequences × class concordance
  concordance_summary.txt            full text report
  criterion_firing_rates.csv         per-criterion fire rates by truth class
  rule_distribution.csv              triggered-rule frequencies
  discrepancies.tsv                  opposite-direction calls (top 10k)
"""

import argparse
import csv
import gzip
from pathlib import Path
from collections import defaultdict, Counter

CLASSES = ["Pathogenic", "Likely_pathogenic", "VUS", "Likely_benign", "Benign"]


def load_truth(path):
    truth = {}
    with open(path) as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for row in rdr:
            key = (row["chrom"], row["pos"], row["ref"], row["alt"])
            truth[key] = row
    return truth


def open_text(path):
    """Open .vcf, .vcf.gz, or .vcf.bgz transparently for text reading."""
    if path.endswith(".gz") or path.endswith(".bgz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_csq_format(header_line: str) -> list[str]:
    """Extract field names from a CSQ INFO header.

    The header form is:
      ##INFO=<ID=CSQ,Number=.,Type=String,Description="... Format: A|B|C|...">
    """
    if "Format: " not in header_line:
        return []
    fmt = header_line.split("Format: ", 1)[1].rstrip('">\n')
    return fmt.split("|")


def class_label(c):
    """Map fastvep ACMG codes (P/LP/VUS/LB/B and long forms) → fixed labels."""
    if not c:
        return None
    c = c.strip()
    if c in ("P", "Pathogenic"):
        return "Pathogenic"
    if c in ("LP", "Likely_pathogenic", "LikelyPathogenic"):
        return "Likely_pathogenic"
    if c in ("VUS", "Uncertain_significance", "UncertainSignificance"):
        return "VUS"
    if c in ("LB", "Likely_benign", "LikelyBenign"):
        return "Likely_benign"
    if c in ("B", "Benign"):
        return "Benign"
    return None


def parse_info_csq(info: str) -> str | None:
    """Extract the CSQ= field value from an INFO column."""
    for piece in info.split(";"):
        if piece.startswith("CSQ="):
            return piece[4:]
    return None


def vep_allele(ref: str, alt: str) -> str:
    """Convert a VCF (REF, ALT) pair to VEP's CSQ-Allele convention.

    VEP strips the leading common prefix between REF and ALT:
      C  → CT      ⇒ "T"        (insertion of T)
      CT → C       ⇒ "-"        (single-base deletion)
      AGGT → A     ⇒ "-"        (multi-base deletion)
      G  → GGGGCC  ⇒ "GGGCC"    (insertion)
      C  → T       ⇒ "T"        (SNV — no shared prefix to strip)
      AT → GC      ⇒ "GC"       (MNV — no shared prefix)

    A pure deletion (alt collapses to empty after stripping) is
    represented as the literal "-" in CSQ.
    """
    if ref == alt:
        return alt
    i = 0
    while i < len(ref) and i < len(alt) and ref[i] == alt[i]:
        i += 1
    new_alt = alt[i:]
    return new_alt if new_alt else "-"


def variant_records(vcf_path):
    """Yield (chrom, pos, ref, alt, csq_field_idx_map, csq_entries_for_alt)
    for each variant in the VCF. Multi-allelic sites are split per ALT.
    Indels are matched on the VEP allele convention (see `vep_allele`).
    """
    csq_idx = None
    with open_text(vcf_path) as f:
        for line in f:
            if line.startswith("##"):
                if "ID=CSQ" in line:
                    csq_fields = parse_csq_format(line)
                    csq_idx = {name: i for i, name in enumerate(csq_fields)}
                continue
            if line.startswith("#"):
                continue
            if csq_idx is None:
                # No CSQ header — bail; nothing to extract.
                return
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom = cols[0].removeprefix("chr")
            pos = cols[1]
            ref = cols[3]
            alts = cols[4].split(",")
            info = cols[7]
            csq_field = parse_info_csq(info)
            if csq_field is None:
                continue
            entries = csq_field.split(",")
            # Group entries by their CSQ Allele (field 0). For indels the
            # VCF ALT (e.g. "CT") doesn't match the CSQ Allele (e.g. "T"
            # for a one-base insertion), so we normalise via vep_allele.
            by_csq_allele: dict[str, list[list[str]]] = defaultdict(list)
            for ent in entries:
                parts = ent.split("|")
                if len(parts) <= csq_idx.get("ACMG", -1):
                    continue
                by_csq_allele[parts[0]].append(parts)
            for alt in alts:
                csq_alt = vep_allele(ref, alt)
                yield (chrom, pos, ref, alt, csq_idx, by_csq_allele.get(csq_alt, []))


def pick_csq(entries: list[list[str]], csq_idx: dict[str, int]):
    """Return (acmg_label, criteria_str, consequence_top, canonical_flag)
    for the picked transcript. Prefers entries with a non-empty ACMG;
    among those, prefers CANONICAL=YES; falls back to first non-empty."""
    if not entries:
        return None, None, "unknown", False
    acmg_i = csq_idx.get("ACMG", -1)
    crit_i = csq_idx.get("ACMG_CRITERIA", -1)
    csq_i = csq_idx.get("Consequence", -1)
    can_i = csq_idx.get("CANONICAL", -1)

    def acmg_of(parts):
        return parts[acmg_i] if 0 <= acmg_i < len(parts) else ""

    populated = [p for p in entries if acmg_of(p)]
    pool = populated or entries
    canon = [p for p in pool if can_i >= 0 and len(p) > can_i and p[can_i] == "YES"]
    chosen = canon[0] if canon else pool[0]

    acmg = acmg_of(chosen) if populated else ""
    crit = chosen[crit_i] if 0 <= crit_i < len(chosen) else ""
    cs = chosen[csq_i] if 0 <= csq_i < len(chosen) else ""
    canonical = (can_i >= 0 and len(chosen) > can_i and chosen[can_i] == "YES")
    # Top consequence is the first listed (CSQ stores them &-separated).
    top = cs.split("&")[0] if cs else "unknown"
    return acmg, crit, top, canonical


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--truth", required=True)
    p.add_argument("--predictions", required=True)
    p.add_argument("--out", default="output")
    args = p.parse_args()

    truth = load_truth(args.truth)
    print(f"Loaded {len(truth):,} truth records")

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    cm = {tc: {pc: 0 for pc in CLASSES + ["NoCall"]} for tc in CLASSES}
    cm_consequence = defaultdict(lambda: {tc: {pc: 0 for pc in CLASSES + ["NoCall"]} for tc in CLASSES})
    cm_chrom = defaultdict(lambda: {tc: {pc: 0 for pc in CLASSES + ["NoCall"]} for tc in CLASSES})
    criterion_fires = {tc: Counter() for tc in CLASSES}
    rule_dist = Counter()
    discrepancies = []
    n_soft_disc = 0
    matched = set()
    n_classified = 0

    for chrom, pos, ref, alt, csq_idx, entries in variant_records(args.predictions):
        key = (chrom, pos, ref, alt)
        if key not in truth:
            continue
        t = truth[key]
        tc_truth = t["normalized_class"]
        acmg, criteria_str, top_csq, _canonical = pick_csq(entries, csq_idx)
        if not acmg:
            cm[tc_truth]["NoCall"] += 1
            cm_chrom[chrom][tc_truth]["NoCall"] += 1
            continue
        pc = class_label(acmg) or "NoCall"
        cm[tc_truth][pc] += 1
        cm_chrom[chrom][tc_truth][pc] += 1
        cm_consequence[top_csq][tc_truth][pc] += 1
        # ACMG_CRITERIA is "&"-joined inside CSQ (because "," is a record
        # separator in INFO and ";" terminates INFO entries). Each token is
        # a met criterion code (e.g. BP4_Moderate, PVS1, PM2_Supporting).
        criteria = [c for c in criteria_str.split("&") if c] if criteria_str else []
        for c in criteria:
            criterion_fires[tc_truth][c] += 1
        # Triggered-rule distribution: derive from the unique combination.
        # We don't get the named rule from VCF output; sort criteria as a
        # surrogate so identical rule signatures group together.
        rule = "+".join(sorted(set(criteria))) if criteria else "(none)"
        rule_dist[rule] += 1
        n_classified += 1
        matched.add(key)
        # An "opposite-direction" call flips benign-side <-> pathogenic-side:
        # these are the review-critical errors and must ALWAYS be logged in
        # full (they feed 04_build_md_review_table.py). "Soft" discrepancies
        # (a directional truth landing on VUS / NoCall) are far more numerous
        # and only the first 10k are kept — but they must never crowd out an
        # opposite-direction call, which the old single 10k cap allowed once
        # the soft calls filled the buffer before the end of the genome.
        is_opposite = (
            tc_truth in ("Pathogenic", "Likely_pathogenic")
            and pc in ("Benign", "Likely_benign")
        ) or (
            tc_truth in ("Benign", "Likely_benign")
            and pc in ("Pathogenic", "Likely_pathogenic")
        )
        is_soft = (
            tc_truth in ("Pathogenic", "Likely_pathogenic")
            and pc in ("VUS", "NoCall")
        )
        if is_opposite or (is_soft and n_soft_disc < 10000):
            if is_soft:
                n_soft_disc += 1
            discrepancies.append((
                chrom, pos, ref, alt, t["gene"], t["review_stars"],
                tc_truth, pc, top_csq, rule, ";".join(criteria)
            ))

    n_truth = len(truth)
    n_unmatched = n_truth - len(matched)

    # ── concordance_matrix.csv ──
    matrix_path = out_dir / "concordance_matrix.csv"
    with matrix_path.open("w") as f:
        w = csv.writer(f)
        w.writerow(["truth"] + CLASSES + ["NoCall"])
        for tcl in CLASSES:
            w.writerow([tcl] + [cm[tcl][pc] for pc in CLASSES + ["NoCall"]])

    # ── concordance_by_chrom.csv ──
    by_chrom_path = out_dir / "concordance_by_chrom.csv"
    with by_chrom_path.open("w") as f:
        w = csv.writer(f)
        w.writerow(["chrom", "truth", "n", "exact", "same_dir", "opposite", "no_call"])
        for chrom in sorted(cm_chrom.keys(), key=lambda x: (len(x), x)):
            for tcl in CLASSES:
                row = cm_chrom[chrom][tcl]
                n = sum(row.values())
                if n == 0:
                    continue
                exact = row[tcl]
                if tcl in ("Pathogenic", "Likely_pathogenic"):
                    same = row["Pathogenic"] + row["Likely_pathogenic"]
                    opp = row["Benign"] + row["Likely_benign"]
                elif tcl == "VUS":
                    same = row["VUS"]
                    opp = 0
                else:
                    same = row["Benign"] + row["Likely_benign"]
                    opp = row["Pathogenic"] + row["Likely_pathogenic"]
                w.writerow([chrom, tcl, n, exact, same, opp, row["NoCall"]])

    # ── concordance_by_consequence.csv ──
    consq_counts = sorted(
        cm_consequence.items(),
        key=lambda kv: -sum(sum(d.values()) for d in kv[1].values()),
    )[:15]
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

    # ── criterion_firing_rates.csv ──
    fr_path = out_dir / "criterion_firing_rates.csv"
    all_codes = sorted({c for cnt in criterion_fires.values() for c in cnt})
    with fr_path.open("w") as f:
        w = csv.writer(f)
        header = ["criterion"] + [f"{tcl}_fired" for tcl in CLASSES]
        w.writerow(header)
        for code in all_codes:
            row = [code] + [criterion_fires[tcl].get(code, 0) for tcl in CLASSES]
            w.writerow(row)

    # ── rule_distribution.csv ──
    rd_path = out_dir / "rule_distribution.csv"
    with rd_path.open("w") as f:
        w = csv.writer(f)
        w.writerow(["criteria_signature", "n"])
        for rule, n in rule_dist.most_common(200):
            w.writerow([rule, n])

    # ── discrepancies.tsv ──
    disc_path = out_dir / "discrepancies.tsv"
    with disc_path.open("w") as f:
        f.write("chrom\tpos\tref\talt\tgene\tstars\ttruth\tpredicted\tconsequence\trule\tmet_criteria\n")
        for d in discrepancies:
            f.write("\t".join(str(x) for x in d) + "\n")

    # ── concordance_summary.txt ──
    totals = {"n": 0, "exact": 0, "same": 0, "opp": 0, "nc": 0}
    summary_path = out_dir / "concordance_summary.txt"
    with summary_path.open("w") as f:
        f.write("ClinVar 2-star+ concordance against fastvep ACMG classifier (real data)\n")
        f.write("=" * 75 + "\n\n")
        f.write(f"Truth records:       {n_truth:,}\n")
        f.write(f"Classified:          {n_classified:,}\n")
        f.write(f"Truth not annotated: {n_unmatched:,}\n\n")
        f.write("Per-class breakdown (entire dataset):\n")
        f.write(f"{'truth':<22} {'n':>8} {'exact':>8} {'same_dir':>10} {'opposite':>10} {'no_call':>8}\n")
        for tcl in CLASSES:
            row = cm[tcl]
            n = sum(row.values())
            exact = row[tcl]
            if tcl in ("Pathogenic", "Likely_pathogenic"):
                same = row["Pathogenic"] + row["Likely_pathogenic"]
                opp = row["Benign"] + row["Likely_benign"]
            elif tcl == "VUS":
                same = row["VUS"]
                opp = 0
            else:
                same = row["Benign"] + row["Likely_benign"]
                opp = row["Pathogenic"] + row["Likely_pathogenic"]
            no_call = row["NoCall"]
            f.write(f"{tcl:<22} {n:>8} {exact:>8} {same:>10} {opp:>10} {no_call:>8}\n")
            for k, v in [("n", n), ("exact", exact), ("same", same), ("opp", opp), ("nc", no_call)]:
                totals[k] += v
        f.write(
            f"\n{'TOTAL':<22} {totals['n']:>8} {totals['exact']:>8} "
            f"{totals['same']:>10} {totals['opp']:>10} {totals['nc']:>8}\n"
        )
        if totals["n"]:
            f.write(f"\nExact-match rate:        {totals['exact']/totals['n']*100:.1f}%\n")
            f.write(f"Same-direction rate:     {totals['same']/totals['n']*100:.1f}%\n")
            f.write(f"Opposite-direction rate: {totals['opp']/totals['n']*100:.1f}%\n")
            f.write(f"NoCall rate:             {totals['nc']/totals['n']*100:.1f}%\n")

        f.write("\nTop 25 criteria signatures (sorted criterion code combos):\n")
        for rule, n in rule_dist.most_common(25):
            f.write(f"  {n:>8}  {rule}\n")

        f.write("\nCriterion fire counts by truth class:\n")
        f.write(f"{'criterion':<22}")
        for tcl in CLASSES:
            f.write(f"  {tcl:>22}")
        f.write("\n")
        for code in all_codes:
            f.write(f"{code:<22}")
            for tcl in CLASSES:
                f.write(f"  {criterion_fires[tcl].get(code, 0):>22}")
            f.write("\n")

    print("\nOutputs:")
    for p in (matrix_path, by_chrom_path, by_csq_path, fr_path, rd_path, disc_path, summary_path):
        print(f"  {p}")
    print()
    print(open(summary_path).read())


if __name__ == "__main__":
    main()
