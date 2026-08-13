#!/usr/bin/env python3
"""Sweep any numeric ACMG config threshold and measure what it costs and buys.

Several ACMG thresholds were set from convention rather than measurement. This
runs the ClinVar 2-star+ benchmark once per value of a config key, holding
everything else at its default, and reports the metrics that decide whether a
threshold is doing its job:

  * the firing count of the criteria that key controls,
  * the effect on the final call in both directions, and
  * separately, the two error types, because they are not interchangeable. A
    false-benign call on a pathogenic variant is a missed diagnosis; a
    false-pathogenic call on a benign one sends a healthy person down a
    clinical pathway. A threshold that trades one for the other at 1:1 is not
    neutral.

Usage:
    sweep_acmg_thresholds.py --sa-dir DIR --out-dir DIR \\
        --sweep "pm2_ad_af_threshold=0.0,1e-6,1e-5,4e-5" \\
        --sweep "bs2_hom_prevalence_threshold=1e-5,1e-3"

Each --sweep is one axis, run independently with every other key left at its
default, so each curve is readable on its own.

What this can and cannot settle
-------------------------------
It can find the best *global default* and show how sharply the answer depends
on the threshold. A flat curve means the knob is not worth curating.

It cannot produce a per-gene-disease threshold, which is what the framework
(Whiffin 2017) actually asks for, because prevalence, penetrance and allelic
heterogeneity are not in this data. And it is fitted to ClinVar labels that
were themselves assigned partly from gnomAD frequency, so a frequency
threshold tuned this way is partly fitted to its own input. Both limits are
why a result here narrows expert review rather than replacing it.
"""

import argparse
import gzip
import json
import os
import subprocess
import sys
import tempfile
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = os.path.dirname(ROOT)  # repo root

BENIGN_CALLS = {"B", "LB"}
PATHOGENIC_CALLS = {"P", "LP"}

# Criteria whose firing counts are worth carrying in the output, so a sweep of
# one key also shows the knock-on effect on the criteria it interacts with.
TRACKED_CRITERIA = (
    "PVS1", "PS1", "PM1", "PM2", "PM2_Supporting", "PP2", "PP3",
    "BA1", "BS1", "BS2", "BP1", "BP4", "BP7",
)


def truth_bucket(clnsig):
    """Collapse a CLNSIG string to pathogenic / benign / uncertain.

    Checked in this order because `Likely_pathogenic` contains "pathogenic"
    and `Conflicting...` contains neither.
    """
    s = clnsig.lower()
    if "conflicting" in s:
        return "uncertain"
    p, b = "pathogenic" in s, "benign" in s
    if p and not b:
        return "pathogenic"
    if b and not p:
        return "benign"
    return "uncertain"


def call_bucket(shorthand):
    if shorthand in PATHOGENIC_CALLS:
        return "pathogenic"
    if shorthand in BENIGN_CALLS:
        return "benign"
    return "uncertain"


def run_one(fastvep, vcf_in, gff3, fasta, sa_dir, config_path, out_vcf):
    cmd = [fastvep, "annotate", "-i", vcf_in, "-o", out_vcf,
           "--gff3", gff3, "--fasta", fasta, "--sa-dir", sa_dir,
           "--acmg", "--pick", "--output-format", "vcf"]
    if config_path:
        cmd += ["--acmg-config", config_path]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)


def score(vcf_path):
    """Read an annotated VCF and return call-level and criterion-level metrics."""
    csq_fields = None
    m = Counter()
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                if line.startswith("##INFO=<ID=CSQ"):
                    csq_fields = line.split("Format: ")[1].rstrip('">\n').split("|")
                continue
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            info = {}
            for part in cols[7].split(";"):
                k, _, v = part.partition("=")
                info[k] = v
            csq = info.get("CSQ")
            if not csq:
                continue
            vals = csq.split(",")[0].split("|")
            if len(vals) != len(csq_fields):
                continue
            rec = dict(zip(csq_fields, vals))

            truth = truth_bucket(info.get("CLNSIG", ""))
            called = call_bucket(rec.get("ACMG", ""))
            codes = set(c for c in rec.get("ACMG_CRITERIA", "").split("&") if c)

            m["n"] += 1
            for code in codes:
                if code in TRACKED_CRITERIA:
                    m["fired_" + code] += 1
            if truth != "uncertain" and called != "uncertain":
                m["same_direction" if truth == called else "opposite_direction"] += 1
            if truth == "pathogenic":
                m["truth_pathogenic"] += 1
                if called == "pathogenic":
                    m["recall_pathogenic"] += 1
                elif called == "benign":
                    m["false_benign"] += 1
            elif truth == "benign":
                m["truth_benign"] += 1
                if called == "benign":
                    m["recall_benign"] += 1
                elif called == "pathogenic":
                    m["false_pathogenic"] += 1
    return m


def parse_sweep(spec):
    """`key=v1,v2,v3` -> (key, [(literal, float), ...]).

    Each value is carried as the caller wrote it as well as parsed. The literal
    is what goes into the generated TOML, so the caller's own notation decides
    the TOML type: `300` is an integer and deserializes into a `u64` field like
    `bp7_max_intron_offset`, while `0.0` is a float and deserializes into an
    `f64` field like `pm2_ad_af_threshold`. Rendering everything from the
    parsed float instead would emit `300.0` and fail every integer-typed key.
    The parsed value is used only for labels and for the JSON summary.
    """
    key, _, values = spec.partition("=")
    if not key or not values:
        raise argparse.ArgumentTypeError(
            f"--sweep must look like key=v1,v2,v3 (got {spec!r})")
    literals = [v.strip() for v in values.split(",")]
    try:
        return key.strip(), [(lit, float(lit)) for lit in literals]
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"--sweep value is not a number: {exc}") from exc


def pct(num, den):
    return 100.0 * num / den if den else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sa-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--sweep", action="append", required=True, type=parse_sweep,
                    help="key=v1,v2,... — repeatable, one axis per flag")
    ap.add_argument("--vcf", default=os.path.join(ROOT, "data/benchmark/clinvar_2star.vcf"))
    ap.add_argument("--gff3", default=os.path.join(ROOT, "test_data/organisms/human/Homo_sapiens.GRCh38.115.gff3"))
    ap.add_argument("--fasta", default=os.path.join(ROOT, "test_data/organisms/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"))
    ap.add_argument("--fastvep", default=os.path.join(ROOT, "target/release/fastvep"))
    ap.add_argument("--name", default="acmg_threshold_sweep",
                    help="basename for the JSON written into --out-dir")
    ap.add_argument("--keep-vcfs", action="store_true",
                    help="retain each run's annotated VCF instead of deleting it")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    results = []
    for key, values in args.sweep:
        for literal, value in values:
            label = f"{key}={literal}"
            with tempfile.NamedTemporaryFile("w", suffix=".toml", delete=False) as cf:
                cf.write(f"{key} = {literal}\n")
                cfg = cf.name
            out_vcf = os.path.join(
                args.out_dir, "sweep_" + label.replace("=", "_") + ".vcf")
            print(f"  running {label} ...", flush=True)
            try:
                run_one(args.fastvep, args.vcf, args.gff3, args.fasta,
                        args.sa_dir, cfg, out_vcf)
                m = score(out_vcf)
            finally:
                os.unlink(cfg)
                if not args.keep_vcfs and os.path.exists(out_vcf):
                    os.unlink(out_vcf)

            row = {"key": key, "value": value, "setting": label, **dict(m)}
            results.append(row)
            print(
                f"    path_recall {pct(m['recall_pathogenic'], m['truth_pathogenic']):5.1f}%"
                f"  benign_recall {pct(m['recall_benign'], m['truth_benign']):5.1f}%"
                f"  false_path {m['false_pathogenic']:>5}"
                f"  false_benign {m['false_benign']:>5}"
                f"  opposite {m['opposite_direction']:>5}",
                flush=True)

    out_json = os.path.join(args.out_dir, args.name + ".json")
    with open(out_json, "w") as fh:
        json.dump(results, fh, indent=1)
    print(f"\nwrote {out_json}")


if __name__ == "__main__":
    sys.exit(main())
