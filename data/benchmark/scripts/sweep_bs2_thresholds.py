#!/usr/bin/env python3
"""Sweep BS2's thresholds and measure what each choice costs and buys.

BS2's two numeric knobs (`bs2_ad_min_ac` for dominant genes,
`bs2_ar_min_hom` plus `bs2_hom_prevalence_threshold` for recessive and
X-linked ones) were set from convention, not measurement. This runs the full
ClinVar 2-star+ benchmark once per setting and reports, for each:

  * how often BS2 fires at all,
  * how often it fires on a variant ClinVar calls pathogenic (the harm), and
  * what that does to the final call, since BS2 only matters when it flips one.

What this can and cannot settle
-------------------------------
It can find the best *global default*, and show how sharply the answer depends
on the threshold. If the curve is flat, the knob is not worth curating.

It cannot produce a per-gene-disease threshold, which is what the ACMG
framework (Whiffin 2017) actually asks for, because prevalence, penetrance and
allelic heterogeneity are not in this data. And it is fitted to ClinVar labels
that were themselves assigned partly from gnomAD frequency, so a threshold
tuned this way is partly fitted to its own input. Both limits are why the
result narrows expert review rather than replacing it.

Usage:
    sweep_bs2_thresholds.py --sa-dir DIR --out-dir DIR [--jobs N]
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


def truth_bucket(clnsig):
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
    """Read an annotated VCF and return BS2 and whole-call metrics."""
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
            if "BS2" in codes:
                m["bs2_fired"] += 1
                m["bs2_on_" + truth] += 1
            if truth != "uncertain" and called != "uncertain":
                if truth == called:
                    m["same_direction"] += 1
                else:
                    m["opposite_direction"] += 1
            if truth == "pathogenic" and called == "benign":
                m["false_benign"] += 1
            if truth == "benign" and called == "benign":
                m["true_benign"] += 1
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sa-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--vcf", default=os.path.join(ROOT, "data/benchmark/clinvar_2star.vcf"))
    ap.add_argument("--gff3", default=os.path.join(ROOT, "test_data/organisms/human/Homo_sapiens.GRCh38.115.gff3"))
    ap.add_argument("--fasta", default=os.path.join(ROOT, "test_data/organisms/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"))
    ap.add_argument("--fastvep", default=os.path.join(ROOT, "target/release/fastvep"))
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # One axis at a time, holding the other at its current default, so each
    # curve is readable on its own.
    settings = []
    for ac in [1, 2, 3, 5, 8, 12, 20, 40, 100]:
        settings.append((f"ad_min_ac={ac}", {"bs2_ad_min_ac": ac}))
    for hom in [1, 2, 3, 5, 10]:
        settings.append((f"ar_min_hom={hom}", {"bs2_ar_min_hom": hom}))
    for prev in [1e-6, 1e-5, 5e-5, 1e-4, 1e-3]:
        settings.append((f"prevalence={prev:.0e}", {"bs2_hom_prevalence_threshold": prev}))

    results = []
    for label, overrides in settings:
        with tempfile.NamedTemporaryFile("w", suffix=".toml", delete=False) as cf:
            for k, v in overrides.items():
                cf.write(f"{k} = {v}\n")
            cfg = cf.name
        out_vcf = os.path.join(args.out_dir, f"sweep_{label.replace('=', '_')}.vcf")
        print(f"  running {label} ...", flush=True)
        run_one(args.fastvep, args.vcf, args.gff3, args.fasta, args.sa_dir, cfg, out_vcf)
        m = score(out_vcf)
        os.unlink(cfg)
        os.unlink(out_vcf)
        row = {"setting": label, **overrides, **dict(m)}
        results.append(row)
        print(f"    BS2 fired {m['bs2_fired']:>7}  "
              f"on-pathogenic {m['bs2_on_pathogenic']:>4}  "
              f"false-benign {m['false_benign']:>4}  "
              f"true-benign {m['true_benign']:>7}", flush=True)

    out_json = os.path.join(args.out_dir, "bs2_threshold_sweep.json")
    with open(out_json, "w") as fh:
        json.dump(results, fh, indent=1)
    print(f"\nwrote {out_json}")


if __name__ == "__main__":
    sys.exit(main())
