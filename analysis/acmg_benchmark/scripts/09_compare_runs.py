#!/usr/bin/env python3
"""Compare ACMG benchmark runs with a baseline, call by call.

`RUN_VERSIONS.md` quotes numbers for what a change did, and numbers typed in by
hand from console output drift. This is the one place they come from: give it a
baseline run and the runs to set against it and it writes

  run_summary.tsv                    one row per run: the headline agreement
                                     counts, the per-class counts, the two error
                                     types, and how many of the reviewed rows
                                     are still opposite-direction
  attribution_<base>_to_<run>.tsv    every changed call, grouped by the
                                     criteria it lost and gained, with how many
                                     moved toward and away from ClinVar's label
  summary.md                         the same two things as tables

The labels are ClinVar's own 2-star-or-better classifications. They are a
reference, not the truth: ClinVar assigns some of them partly from the same
gnomAD frequency and ClinVar-neighbour evidence the classifier reads, so
agreement with them is partly circular and a threshold fitted to them is partly
fitted to its own input. "Toward" and "away" below mean "closer to ClinVar's
label" and "further from it", nothing stronger.

Usage:
  09_compare_runs.py --baseline v36 --runs v37 v37_bs2off --out analysis/acmg_benchmark/results/v37/
  09_compare_runs.py --baseline v36 --runs v38 --out /tmp/v38 --labels data/benchmark/clinvar_2star_labels.tsv

Runs are directory names under `data/benchmark/output_<run>/`. `--out` is
required and is taken relative to the working directory.

Every changed call is counted once in `changed` and lands in exactly one of
`toward`, `away`, `level` (the same distance from the label before and after,
as when P becomes B against a VUS label) or `unscored` (a NoCall on either side,
or no ClinVar label for the variant), so those four add up to `changed`.

`reviewed_rows_opposite` is -1 in run_summary.tsv, and n/a in summary.md, for a
run whose `discrepancies.tsv` or reviewed-rows file is missing.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BENCH = ROOT / "data" / "benchmark"
DEFAULT_LABELS = BENCH / "clinvar_2star_labels.tsv"
LEGACY_LABELS = BENCH / "clinvar_2star_truth.tsv"
DEFAULT_REVIEWED = HERE.parent / "data" / "round6_reviewed_rows.tsv"

ATTRIBUTION_ROWS_SHOWN = 25
CLASSES = ["Pathogenic", "Likely_pathogenic", "VUS", "Likely_benign", "Benign"]
SHORT = {"P": "Pathogenic", "LP": "Likely_pathogenic", "VUS": "VUS", "LB": "Likely_benign", "B": "Benign"}
# One axis, pathogenic positive: a call is "closer" to a label when the distance
# between their ranks shrinks.
RANK = {"Pathogenic": 2, "Likely_pathogenic": 1, "VUS": 0, "Likely_benign": -1, "Benign": -2}


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def read_labels(path: Path) -> dict[tuple[str, str, str, str], str]:
    out = {}
    with path.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            out[(row["chrom"], row["pos"], row["ref"], row["alt"])] = row["normalized_class"]
    return out


def read_matrix(run_dir: Path) -> dict[str, dict[str, int]]:
    """concordance_matrix.csv: one row per ClinVar label, one column per call."""
    with (run_dir / "concordance_matrix.csv").open() as f:
        rows = list(csv.reader(f))
    header = rows[0][1:]
    return {r[0]: {c: int(v) for c, v in zip(header, r[1:])} for r in rows[1:]}


def metrics(matrix: dict[str, dict[str, int]]) -> dict[str, int]:
    p, lp = matrix["Pathogenic"], matrix["Likely_pathogenic"]
    lb, b = matrix["Likely_benign"], matrix["Benign"]
    vus = matrix["VUS"]
    exact = sum(matrix[c][c] for c in CLASSES)
    same = (
        p["Pathogenic"] + p["Likely_pathogenic"]
        + lp["Pathogenic"] + lp["Likely_pathogenic"]
        + vus["VUS"]
        + lb["Likely_benign"] + lb["Benign"]
        + b["Likely_benign"] + b["Benign"]
    )
    opposite = (
        p["Likely_benign"] + p["Benign"] + lp["Likely_benign"] + lp["Benign"]
        + lb["Pathogenic"] + lb["Likely_pathogenic"] + b["Pathogenic"] + b["Likely_pathogenic"]
    )
    no_call = sum(matrix[c]["NoCall"] for c in CLASSES)
    return {
        "classified": sum(matrix[c][k] for c in CLASSES for k in CLASSES),
        "exact": exact,
        "same_direction": same,
        "opposite": opposite,
        "no_call": no_call,
        "p_same": p["Pathogenic"] + p["Likely_pathogenic"],
        "lp_same": lp["Pathogenic"] + lp["Likely_pathogenic"],
        "vus_exact": vus["VUS"],
        "lb_same": lb["Likely_benign"] + lb["Benign"],
        "b_same": b["Likely_benign"] + b["Benign"],
        # The two errors are not interchangeable: the first is a missed
        # diagnosis, the second sends a healthy person down a clinical pathway.
        "false_benign": p["Likely_benign"] + p["Benign"] + lp["Likely_benign"] + lp["Benign"],
        "false_pathogenic": lb["Pathogenic"] + lb["Likely_pathogenic"] + b["Pathogenic"] + b["Likely_pathogenic"],
    }


def reviewed_still_opposite(run_dir: Path, reviewed: set[tuple[str, str, str, str]]) -> int:
    """How many of the reviewed rows the run still lists as opposite-direction.

    `discrepancies.tsv` carries every opposite-direction call in full (and only
    the first 10,000 of the softer ones), so membership there is exact.
    """
    path = run_dir / "discrepancies.tsv"
    if not path.exists() or not reviewed:
        return -1
    hits = 0
    with path.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pc = row["predicted"]
            tc = row.get("clinvar") or row.get("truth") or row.get("clinvar_class", "")
            opposite = (tc in ("Pathogenic", "Likely_pathogenic") and pc in ("Benign", "Likely_benign")) or (
                tc in ("Benign", "Likely_benign") and pc in ("Pathogenic", "Likely_pathogenic")
            )
            if opposite and (row["chrom"], row["pos"], row["ref"], row["alt"]) in reviewed:
                hits += 1
    return hits


def read_reviewed(path: Path) -> set[tuple[str, str, str, str]]:
    if not path.exists():
        return set()
    with path.open() as f:
        return {(r["chrom"], r["pos"], r["ref"], r["alt"]) for r in csv.DictReader(f, delimiter="\t")}


def base_code(code: str) -> str:
    """PM1_Moderate -> PM1, BP4_Very_Strong -> BP4: a change of strength is not a change of criterion."""
    return re.sub(r"_(Supporting|Moderate|Strong|Very_Strong)$", "", code)


def picked_calls(ev, diff, vcf_path: Path) -> dict:
    """key -> {"acmg": short call, "criteria": set of base codes}, through the same picker the concordance report uses."""
    out = {}
    for key, row in diff.calls(ev, str(vcf_path)).items():
        out[key] = {"acmg": row["acmg"], "criteria": {base_code(c) for c in row["criteria"].split("&") if c}}
    return out


def label_of(call: str) -> str | None:
    return SHORT.get(call) if call in SHORT else call if call in RANK else None


ATTRIBUTION_FIELDS = ["criteria_lost_to_gained", "changed", "toward", "away", "level", "unscored",
                      "label_pathogenic", "label_vus", "label_benign", "opposite_before", "opposite_after"]


def attribute(base_calls: dict, run_calls: dict, labels: dict) -> list[dict]:
    """Group changed calls by the criteria they lost and gained."""
    groups: dict[str, Counter] = defaultdict(Counter)
    for key, a in base_calls.items():
        b = run_calls.get(key)
        if b is None or a["acmg"] == b["acmg"]:
            continue
        lost = ",".join(sorted(a["criteria"] - b["criteria"])) or "-"
        gained = ",".join(sorted(b["criteria"] - a["criteria"])) or "-"
        g = groups[f"{lost} -> {gained}"]
        g["changed"] += 1
        label = labels.get(key)
        if label is None:
            g["unscored"] += 1
            continue
        pathogenic_label = label in ("Pathogenic", "Likely_pathogenic")
        benign_label = label in ("Benign", "Likely_benign")
        g["label_pathogenic"] += pathogenic_label
        g["label_benign"] += benign_label
        g["label_vus"] += label == "VUS"
        # A NoCall is no direction at all, so it is never opposite to anything,
        # and a change into or out of one has no distance to compare.
        ca, cb = label_of(a["acmg"]), label_of(b["acmg"])
        # Opposite-direction before and after, the review-critical count.
        if ca is not None and ((pathogenic_label and RANK[ca] < 0) or (benign_label and RANK[ca] > 0)):
            g["opposite_before"] += 1
        if cb is not None and ((pathogenic_label and RANK[cb] < 0) or (benign_label and RANK[cb] > 0)):
            g["opposite_after"] += 1
        if ca is None or cb is None:
            g["unscored"] += 1
            continue
        before, after = abs(RANK[ca] - RANK[label]), abs(RANK[cb] - RANK[label])
        if after < before:
            g["toward"] += 1
        elif after > before:
            g["away"] += 1
        else:
            g["level"] += 1
    rows = []
    for key, g in groups.items():
        rows.append({"criteria_lost_to_gained": key, **{k: g.get(k, 0) for k in ATTRIBUTION_FIELDS[1:]}})
    rows.sort(key=lambda r: -r["changed"])
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baseline", required=True, help="run name under data/benchmark/output_<name>/")
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--labels", type=Path, default=None, help="ClinVar 2-star+ labels TSV")
    ap.add_argument("--reviewed", type=Path, default=None,
                    help=f"TSV of reviewed rows (chrom pos ref alt); default {DEFAULT_REVIEWED.relative_to(ROOT)}")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--no-attribution", action="store_true", help="skip the per-variant diff (fast; headline only)")
    args = ap.parse_args()

    labels_path = args.labels or (DEFAULT_LABELS if DEFAULT_LABELS.exists() else LEGACY_LABELS)
    if not labels_path.exists():
        ap.error(f"labels file not found: {labels_path}")
    # A reviewed-rows file that was asked for and is missing is an error; one that
    # was not asked for and is missing leaves the column at -1, "not computed".
    if args.reviewed is not None and not args.reviewed.exists():
        ap.error(f"--reviewed file not found: {args.reviewed}")
    reviewed_path = args.reviewed or DEFAULT_REVIEWED
    if not reviewed_path.exists():
        sys.stderr.write(f"warning: no reviewed-rows file at {reviewed_path}; reviewed_rows_opposite will be n/a\n")
    reviewed = read_reviewed(reviewed_path)
    args.out.mkdir(parents=True, exist_ok=True)

    names = [args.baseline] + [r for r in args.runs if r != args.baseline]
    summary = []
    for name in names:
        run_dir = BENCH / f"output_{name}"
        m = metrics(read_matrix(run_dir))
        m["run"] = name
        m["reviewed_rows_opposite"] = reviewed_still_opposite(run_dir, reviewed)
        summary.append(m)

    cols = ["run", "classified", "exact", "same_direction", "opposite", "no_call", "p_same", "lp_same", "vus_exact",
            "lb_same", "b_same", "false_benign", "false_pathogenic", "reviewed_rows_opposite"]
    with (args.out / "run_summary.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in summary:
            w.writerow({c: row[c] for c in cols})

    # One sentence to a line, so that a diff of two of these files shows which
    # sentence changed.
    md = ["# Run comparison", "",
          f"Baseline `{args.baseline}`.",
          f"Labels: ClinVar 2-star+ classifications (`{labels_path.name}`), a reference and not ground truth.", "",
          "| run | classified | exact | same-direction | opposite | P same | LP same | VUS exact | LB same | B same "
          "| P/LP called B/LB | B/LB called P/LP | reviewed rows opposite |",
          "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in summary:
        md.append("| {run} | {classified:,} | {exact:,} | {same_direction:,} | {opposite} | {p_same:,} | {lp_same:,} "
                  "| {vus_exact:,} | {lb_same:,} | {b_same:,} | {false_benign} | {false_pathogenic} "
                  "| {reviewed} |".format(reviewed="n/a" if r["reviewed_rows_opposite"] < 0 else r["reviewed_rows_opposite"], **r))

    if not args.no_attribution:
        ev = load_module("acmg_eval", HERE / "03_evaluate_concordance.py")
        diff = load_module("acmg_diff", HERE / "08_diff_calls.py")
        labels = read_labels(labels_path)
        sys.stderr.write(f"baseline {args.baseline}\n")
        base_calls = picked_calls(ev, diff, BENCH / f"output_{args.baseline}" / "clinvar_2star.fastvep.vcf.gz")
        for name in names[1:]:
            sys.stderr.write(f"run {name}\n")
            run_calls = picked_calls(ev, diff, BENCH / f"output_{name}" / "clinvar_2star.fastvep.vcf.gz")
            rows = attribute(base_calls, run_calls, labels)
            out = args.out / f"attribution_{args.baseline}_to_{name}.tsv"
            with out.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=ATTRIBUTION_FIELDS, delimiter="\t", lineterminator="\n")
                w.writeheader()
                w.writerows(rows)
            total = Counter()
            for r in rows:
                for k in ATTRIBUTION_FIELDS[1:]:
                    total[k] += r[k]
            md += ["", f"## `{args.baseline}` to `{name}`", "",
                   f"{total['changed']:,} calls changed: {total['toward']:,} toward ClinVar's label, "
                   f"{total['away']:,} away, {total['level']:,} the same distance from it, "
                   f"{total['unscored']:,} not scored (a NoCall on either side, or no ClinVar label).",
                   f"Opposite-direction {total['opposite_before']} before, "
                   f"{total['opposite_after']} after, among the changed calls.", "",
                   "| criteria lost to gained | changed | toward | away | level | not scored | ClinVar P/LP | ClinVar VUS | ClinVar B/LB | opposite before | after |",
                   "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
            for r in rows[:ATTRIBUTION_ROWS_SHOWN]:
                md.append("| {criteria_lost_to_gained} | {changed:,} | {toward:,} | {away:,} | {level:,} | {unscored:,} "
                          "| {label_pathogenic:,} | {label_vus:,} | {label_benign:,} | {opposite_before} | {opposite_after} |".format(**r))
            if len(rows) > ATTRIBUTION_ROWS_SHOWN:
                md += ["", f"The {ATTRIBUTION_ROWS_SHOWN} largest of {len(rows):,} groups; all of them, and the totals above, "
                           f"are in `{out.name}`."]

    (args.out / "summary.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
