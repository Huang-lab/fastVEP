#!/usr/bin/env python3
"""Replay a BS2 widening call by call, and ask what a guard would have kept.

`v37_adtag` was a run in which BS2's dominant branch read ClinGen's `/AD` tag
instead of the disease name (limitation 20 in `docs/ACMG.md`). It came from a
locally patched build that is not in the repository, so this script does not
make it again. It reads the annotated VCFs of two runs that exist, the shipped
one and that one, finds every variant where the second carries BS2 and the first
does not (an "addition"), and works out the call each addition would have had if
a guard had withheld its BS2.

That call is the baseline's, because BS2 is the only criterion the two runs
differ in for those variants. The script counts the additions where another
criterion differs as well, and says so, so that the assumption is checked and
not trusted. The replay covers additions only: a variant that loses BS2 in the
widened run keeps the baseline's BS2 here, which is what "baseline plus guarded
additions" means.

The guards are the dominant carrier bar on the half-decade grid the sweeps in
`docs/ACMG.md` use, and "BS1 and BA1 are not already firing", which asks whether
BS2 is counting a frequency the baseline's own BS1 or BA1 already counted. A
guard chosen after seeing which variants a widening harms is fitted to those
variants, and the table says how many that was.

The lower bound is the classifier's own: `poisson_lower_95` in
`crates/fastvep-classification/src/criteria/benign_strong.rs` (a 95% one-sided
bound, tabulated up to 20 events and Wilson-Hilferty above), divided by the
individuals surveyed. Change one and the other has to change. Hemizygote counts
are not in the VCF, so a non-PAR X variant is divided by AN / 2 here, where the
classifier counts its XY individuals once, which makes the bound a little high
there.

The labels are ClinVar's own 2-star-or-better classifications, a reference and
not the truth, so "toward" and "pushed toward benign" mean closer to or further
from ClinVar's label and nothing stronger.

Usage:
  10_replay_bs2_guards.py --baseline v37 --run v37_adtag \\
      --out analysis/acmg_benchmark/results/v37/

Writes, under `--out`:
  bs2_tag_replay.tsv                  one row per guard
  bs2_tag_additions_pathogenic.tsv    the additions ClinVar calls Pathogenic or
                                      Likely pathogenic, with the counts the bound
                                      comes from
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from itertools import zip_longest
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BENCH = ROOT / "data" / "benchmark"
DEFAULT_LABELS = BENCH / "clinvar_2star_labels.tsv"

RANK = {"P": 2, "LP": 1, "VUS": 0, "LB": -1, "B": -2}
LABEL = {"Pathogenic": "P", "Likely_pathogenic": "LP", "VUS": "VUS", "Likely_benign": "LB", "Benign": "B"}

# bs2_ad_prevalence_threshold's default, and the half-decade grid the sweep used.
SHIPPED_BAR = 3e-4
BARS = (3e-4, 1e-3, 3e-3)

CSQ = re.compile(r"(?:^|;)CSQ=([^;]*)")
GNOMAD = re.compile(r"(?:^|;)FV_GNOMAD=([^;]*)")
# The columns this script reads, found by name in each file's own header and not
# by position, so that a column added to CSQ cannot shift them into a plausible
# wrong answer.
CSQ_COLUMNS = ("SYMBOL", "HGVSc", "HGVSp", "ACMG", "ACMG_CRITERIA")
GNOMAD_COLUMNS = ("ALL_AC", "ALL_AN", "ALL_HC")


def poisson_lower_95(observed: int) -> float:
    """Port of `poisson_lower_95` in benign_strong.rs."""
    table = [
        0.0, 0.0513, 0.3554, 0.8177, 1.3663, 1.9702, 2.6130, 3.2853, 3.9808, 4.6952, 5.4254,
        6.1690, 6.9242, 7.6896, 8.4639, 9.2463, 10.0360, 10.8324, 11.6350, 12.4432, 13.2547,
    ]
    if observed < len(table):
        return table[observed]
    k = float(observed)
    term = 1.0 - 1.0 / (9.0 * k) - 1.645 / (3.0 * k**0.5)
    # The product in the order the Rust writes it, so that the two agree to the last bit.
    return max(k * term * term * term, 0.0)


def carrier_bound(ac: int, an: int, hom: int) -> float:
    """The carrier frequency's lower bound, as `evaluate_bs2` computes it for a dominant gene."""
    individuals = max(an // 2, 1)
    return poisson_lower_95(max(ac - hom, 0)) / individuals


def base_code(code: str) -> str:
    """PM1_Moderate -> PM1: a change of strength is not a change of criterion."""
    return re.sub(r"_(Supporting|Moderate|Strong|Very_Strong)$", "", code)


def header_columns(path: Path, header: list[str], info_id: str, wanted: tuple[str, ...]) -> dict[str, int]:
    """Where each wanted column of an INFO field sits, from the field's `Format: A|B|C` header line."""
    for line in header:
        if line.startswith(f"##INFO=<ID={info_id},"):
            m = re.search(r"Format: ([^\">]+)", line)
            if m:
                names = m.group(1).strip().split("|")
                missing = [w for w in wanted if w not in names]
                if missing:
                    sys.exit(f"{path}: the {info_id} header has no column named {', '.join(missing)}")
                return {w: names.index(w) for w in wanted}
    sys.exit(f"{path}: no {info_id} header line with a Format, so its columns cannot be found")


def stream(path: Path):
    """key -> (gene, hgvsc, hgvsp, call, criteria, ac, an, hom), in file order, for the first CSQ entry."""
    header: list[str] = []
    with gzip.open(path, "rt") as f:
        csq_at = gnomad_at = None
        for line in f:
            if line[0] == "#":
                header.append(line)
                continue
            if csq_at is None:
                csq_at = header_columns(path, header, "CSQ", CSQ_COLUMNS)
                gnomad_at = header_columns(path, header, "FV_GNOMAD", GNOMAD_COLUMNS)
            t = line.rstrip("\n").split("\t", 8)
            m = CSQ.search(t[7])
            if not m:
                continue
            csq = m.group(1).split(",")[0].split("|")
            g = GNOMAD.search(t[7])
            counts = (None, None, None)
            if g:
                p = g.group(1).split("|")
                try:
                    counts = tuple(int(p[gnomad_at[c]]) for c in GNOMAD_COLUMNS)
                except (IndexError, ValueError):
                    pass
            crit = csq[csq_at["ACMG_CRITERIA"]]
            yield (
                (t[0], t[1], t[3], t[4]),
                csq[csq_at["SYMBOL"]], csq[csq_at["HGVSc"]], csq[csq_at["HGVSp"]],
                csq[csq_at["ACMG"]], {base_code(c) for c in crit.split("&") if c},
                *counts,
            )


def read_labels(path: Path) -> dict:
    with path.open() as f:
        return {(r["chrom"], r["pos"], r["ref"], r["alt"]): r["normalized_class"] for r in csv.DictReader(f, delimiter="\t")}


def same_direction(call: str, label: str) -> bool:
    return (RANK[call] > 0 and RANK[label] > 0) or (RANK[call] < 0 and RANK[label] < 0) or (RANK[call] == 0 == RANK[label])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baseline", required=True, help="run name under data/benchmark/output_<name>/, the shipped behaviour")
    ap.add_argument("--run", required=True, help="the run whose BS2 was widened")
    ap.add_argument("--labels", type=Path, default=DEFAULT_LABELS, help="ClinVar 2-star+ labels TSV")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    if not args.labels.exists():
        ap.error(f"labels file not found: {args.labels}")
    paths = [BENCH / f"output_{n}" / "clinvar_2star.fastvep.vcf.gz" for n in (args.baseline, args.run)]
    for p in paths:
        if not p.exists():
            ap.error(f"annotated VCF not found: {p}")
    labels = read_labels(args.labels)

    additions, other_changes, lost = [], 0, 0
    totals = {"baseline": [0, 0], "run": [0, 0]}  # exact and same-direction, over the variants both runs call
    for a, b in zip_longest(stream(paths[0]), stream(paths[1])):
        if a is None or b is None:
            sys.exit(f"the two runs do not list the same number of variants: {args.run if b is None else args.baseline} ends first")
        if a[0] != b[0]:
            sys.exit(f"the two runs do not list the same variants in the same order: {a[0]} against {b[0]}")
        (key, gene, hgvsc, hgvsp, call_a, crit_a, ac, an, hom), (_, _, _, _, call_b, crit_b, *_rest) = a, b
        label = LABEL.get(labels.get(key))
        if call_a not in RANK or call_b not in RANK:
            continue
        if label is not None:
            for side, call in (("baseline", call_a), ("run", call_b)):
                totals[side][0] += call == label
                totals[side][1] += same_direction(call, label)
        if "BS2" in crit_b and "BS2" not in crit_a:
            if ac is None or an is None or hom is None:
                continue
            additions.append({
                "key": key, "gene": gene, "hgvsc": hgvsc, "hgvsp": hgvsp, "label": label,
                "base": call_a, "run": call_b, "crit": crit_a, "ac": ac, "an": an, "hom": hom,
                "bound": carrier_bound(ac, an, hom), "other": (crit_a ^ crit_b) - {"BS2"},
            })
        elif call_a != call_b:
            other_changes += 1
            lost += "BS2" in crit_a and "BS2" not in crit_b
    for side, name in (("baseline", args.baseline), ("run", args.run)):
        print(f"{name}: {totals[side][0]:,} exact and {totals[side][1]:,} same-direction matches with ClinVar, "
              "over the variants both runs call")
    mixed = sum(1 for x in additions if x["other"])
    below_bar = sum(1 for x in additions if x["bound"] <= SHIPPED_BAR)
    changed = sum(1 for x in additions if x["base"] != x["run"])
    print(f"{len(additions):,} variants gain BS2 in {args.run} and {changed:,} of them change call; "
          f"{mixed} differ in another criterion too, which the replay does not separate")
    print(f"{below_bar} additions have a bound at or below the shipped bar {SHIPPED_BAR:g} "
          "(0 is expected, and a non-zero count is a non-PAR X variant or a port that has drifted)")
    print(f"{other_changes:,} further calls differ without a BS2 addition ({lost} of them lose BS2), and are not replayed")

    def effect(x):
        """(exact, same_direction, pushed toward benign, resolved) for one addition that changes the call."""
        label = x["label"]
        if label is None or x["base"] == x["run"]:
            return 0, 0, False, False
        exact = int(x["run"] == label) - int(x["base"] == label)
        same = int(same_direction(x["run"], label)) - int(same_direction(x["base"], label))
        pushed = RANK[label] > 0 and RANK[x["run"]] < RANK[x["base"]]
        resolved = RANK[label] < 0 and x["base"] == "VUS" and RANK[x["run"]] < 0
        return exact, same, pushed, resolved

    guards = [("none: BS2 as the widened run made it", lambda x: True)]
    for bar in BARS:
        guards.append((f"carrier bar above {bar:g}", lambda x, bar=bar: x["bound"] > bar))
    guards.append(("BS1 and BA1 not already firing", lambda x: not ({"BS1", "BA1"} & x["crit"])))
    for bar in BARS[1:]:
        guards.append((f"carrier bar above {bar:g} and BS1 and BA1 not already firing",
                       lambda x, bar=bar: x["bound"] > bar and not ({"BS1", "BA1"} & x["crit"])))

    rows = []
    for name, keep in guards:
        kept = [x for x in additions if keep(x)]
        exact = same = resolved = 0
        pushed = []
        for x in kept:
            e, s, p, r = effect(x)
            exact += e
            same += s
            resolved += r
            if p:
                pushed.append(f"{x['gene']} {x['base']}>{x['run']}")
        rows.append({
            "guard": name, "additions_kept": len(kept), "exact_delta": exact, "same_direction_delta": same,
            "benign_labelled_vus_resolved": resolved, "pathogenic_pushed_toward_benign": len(pushed),
            "which": ",".join(pushed) or "-",
        })

    args.out.mkdir(parents=True, exist_ok=True)
    cols = ["guard", "additions_kept", "exact_delta", "same_direction_delta", "benign_labelled_vus_resolved",
            "pathogenic_pushed_toward_benign", "which"]
    with (args.out / "bs2_tag_replay.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    pcols = ["chrom", "pos", "ref", "alt", "gene", "hgvsc", "hgvsp", "clinvar", "baseline_call", "run_call",
             "ac", "homozygotes", "an", "carriers", "carrier_lower_95", "baseline_criteria"]
    with (args.out / "bs2_tag_additions_pathogenic.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(pcols)
        pathogenic = [x for x in additions if x["label"] in ("P", "LP")]
        for x in sorted(pathogenic, key=lambda x: (x["gene"], x["key"])):
            w.writerow([*x["key"], x["gene"], x["hgvsc"], x["hgvsp"], x["label"], x["base"], x["run"],
                        x["ac"], x["hom"], x["an"], x["ac"] - x["hom"], f"{x['bound']:.3g}",
                        "&".join(sorted(x["crit"])) or "-"])

    width = max(len(r["guard"]) for r in rows)
    for r in rows:
        print(f"{r['guard']:<{width}}  kept {r['additions_kept']:>6,}  exact {r['exact_delta']:>+6,}  "
              f"same-direction {r['same_direction_delta']:>+6,}  VUS resolved {r['benign_labelled_vus_resolved']:>5,}  "
              f"Pathogenic or Likely pathogenic pushed toward benign {r['pathogenic_pushed_toward_benign']} "
              f"{'' if r['which'] == '-' else r['which']}")


if __name__ == "__main__":
    main()
