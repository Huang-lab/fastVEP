#!/usr/bin/env python3
"""Score the final classifier against a later ClinVar release, and on the variants that were not in the benchmark.

Every threshold, exception and rule in the classifier was chosen against one frozen ClinVar
release (the file dated 2026-04-04, 673,660 2-star+ variants) and the discordant rows it
produced. That makes the benchmark an in-sample measurement for everything that was tuned.
A later release holds variants that no sweep, review round or run used, which is the nearest
thing to a held-out set available without new data:

  new           not in the benchmark set (which holds only variants of two stars or better).
                Almost all were in ClinVar in April at a lower review status (nine in ten as a
                single-submitter record) and have reached two stars since: no version, review
                round or sweep ever scored them
  unchanged     in both releases with the same ClinVar class
  tier_change   in both, class moved within its direction (LB <-> B, LP <-> P)
  direction_change  in both, class moved across a direction (VUS <-> P/LP <-> B/LB)

One restriction is forced by the data and is not a choice. The gnomAD databases the classifier
reads were extracted for the regions of the April release (`clinvar_2star_regions.bed`, built from
the VCF positions of its variants), so a variant outside those regions has no gnomAD record, and
`pm2_absent_when_no_record` would read that as "absent from gnomAD" and fire PM2 on a gap in our
extract. Only variants whose VCF position lies inside a region are kept, which is the rule that
decided which April variants had gnomAD data, applied to every variant alike. How many variants
this drops is reported, not hidden.

The ClinVar-derived inputs the classifier reads (`clinvar_protein.oga`, built from the late-April
variant_summary, and `clinvar.osa`, from the 2026-06-27 release) are not rebuilt: they are what a
laboratory pipeline with a few-months-old index would hold. The protein index is built from every
classified missense record of the variant_summary whatever its review status, so it holds the April
record of each new variant that was then a classified missense variant (105,920 of the 107,528 new
variants have some record there, by `13_manuscript_stats.py`) and lacks the records added or
reclassified since. For a new variant that lowers PS1 and PM5 a little (other variants at its
residue may have been classified since) and lowers or raises PM1 (fewer pathogenic neighbours lower
it, fewer benign neighbours raise it); its own April record, where that record is pathogenic, is in
BP1's count as it is for the April set. A new variant with a pathogenic record in the allele-level
annotation (27 June) has that record left out of the PS1 and PM1 counts although the April index
does not hold it where it was not pathogenic then, which lowers both further. The manuscript counts
how many of the opposite-direction calls on the new variants rest on those criteria, which bounds
what the stale index can have done.

ClinVar's 2-star+ classifications are a reference and not ground truth.

Usage
  14_temporal_holdout.py prepare   build the covered subset and the strata
  14_temporal_holdout.py run       annotate and score it with run_benchmark.sh
  14_temporal_holdout.py table     write the per-variant table of that run
  14_temporal_holdout.py blockers  why BA1 and BS1 did not fire on the uncalled Benign variants of 2% or more
                                   (needs the per-variant tables; writes temporal_high_af_blockers_v37.csv)
  14_temporal_holdout.py all       prepare, run and table
"""

from __future__ import annotations

import argparse
import bisect
import csv
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BENCH = ROOT / "data/benchmark"
APRIL_LABELS = BENCH / "clinvar_2star_labels.tsv"
CURRENT_RAW = BENCH / "sa_sources/clinvar_20261004.vcf.gz"
CURRENT_DIR = BENCH / "clinvar_20261004"
REGIONS = BENCH / "sa_sources/clinvar_2star_regions.bed"
SUBSET_DIR = BENCH / "clinvar_20261004_covered"
RUN_NAME = "v37_oct2026"
RUN_DIR = BENCH / f"output_{RUN_NAME}"
PER_VARIANT = BENCH / "manuscript_runs/v37_oct2026_per_variant.tsv.gz"

DIRECTION = {"Pathogenic": "P/LP", "Likely_pathogenic": "P/LP", "VUS": "VUS", "Likely_benign": "B/LB", "Benign": "B/LB"}


def load_regions() -> dict[str, tuple[list[int], list[int]]]:
    out: dict[str, tuple[list[int], list[int]]] = {}
    with REGIONS.open() as f:
        for line in f:
            c, s, e = line.split()[:3]
            c = c.removeprefix("chr")
            out.setdefault(c, ([], []))
            out[c][0].append(int(s))
            out[c][1].append(int(e))
    return out


def covered(regions, chrom: str, pos: int) -> bool:
    """The VCF position lies inside a region the gnomAD extract fetched (BED start is 0-based)."""
    if chrom not in regions:
        return False
    starts, ends = regions[chrom]
    i = bisect.bisect_right(starts, pos - 1) - 1
    return i >= 0 and pos <= ends[i]


def stratum(april: dict, key: tuple, cls: str) -> str:
    if key not in april:
        return "new"
    before = april[key]
    if before == cls:
        return "unchanged"
    return "tier_change" if DIRECTION[before] == DIRECTION[cls] else "direction_change"


def prepare() -> None:
    if not (CURRENT_DIR / "clinvar_2star_labels.tsv").exists():
        CURRENT_DIR.mkdir(parents=True, exist_ok=True)
        subprocess.run([sys.executable, str(HERE / "01_extract_clinvar_2star.py"), str(CURRENT_RAW), str(CURRENT_DIR)], check=True)
    regions = load_regions()
    april = {}
    with APRIL_LABELS.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            april[(r["chrom"], r["pos"], r["ref"], r["alt"])] = r["normalized_class"]

    keep: dict[tuple, str] = {}
    labels_rows = []
    counts = Counter()
    with (CURRENT_DIR / "clinvar_2star_labels.tsv").open() as f:
        rd = csv.DictReader(f, delimiter="\t")
        header = rd.fieldnames
        for r in rd:
            key = (r["chrom"], r["pos"], r["ref"], r["alt"])
            st = stratum(april, key, r["normalized_class"])
            counts[("all", st)] += 1
            if not covered(regions, r["chrom"], int(r["pos"])):
                counts[("dropped_outside_extract", st)] += 1
                continue
            keep[key] = st
            labels_rows.append(r)
            counts[("kept", st)] += 1

    SUBSET_DIR.mkdir(parents=True, exist_ok=True)
    with (SUBSET_DIR / "clinvar_2star_labels.tsv").open("w") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(labels_rows)
    n_vcf = 0
    with (CURRENT_DIR / "clinvar_2star.vcf").open() as src, (SUBSET_DIR / "clinvar_2star.vcf").open("w") as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue
            c = line.split("\t", 5)
            key = (c[0], c[1], c[3], c[4])
            if key in keep:
                dst.write(line)
                n_vcf += 1
    with (SUBSET_DIR / "strata.tsv").open("w") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["chrom", "pos", "ref", "alt", "stratum"])
        for (c, p, r, a), st in keep.items():
            w.writerow([c, p, r, a, st])
    assert n_vcf == len(keep), (n_vcf, len(keep))
    print(f"kept {len(keep):,} of {sum(v for (k, _), v in counts.items() if k == 'all'):,} current 2-star+ variants")
    for (k, st), v in sorted(counts.items()):
        print(f"  {k:26} {st:17} {v:>9,}")


def run() -> None:
    env = os.environ.copy()
    env.update(INPUT=str(SUBSET_DIR / "clinvar_2star.vcf"), LABELS=str(SUBSET_DIR / "clinvar_2star_labels.tsv"), OUT_DIR=str(RUN_DIR))
    subprocess.run(["bash", str(HERE / "run_benchmark.sh"), RUN_NAME], env=env, check=True, cwd=ROOT)


def table() -> None:
    subprocess.run([sys.executable, str(HERE / "12_per_variant_table.py"), "--run", str(RUN_DIR),
                    "--labels", str(SUBSET_DIR / "clinvar_2star_labels.tsv"), "--out", str(PER_VARIANT)], check=True)


JSON_CLASS = {"Pathogenic": "Pathogenic", "LikelyPathogenic": "Likely_pathogenic", "UncertainSignificance": "VUS",
              "LikelyBenign": "Likely_benign", "Benign": "Benign"}


def bs1_reason(criterion: dict) -> str:
    """Why BS1 (the lower of the two frequency bars) did or did not fire, from the details of its criterion."""
    d = criterion.get("details", {})
    blocked = d.get("frequency_blocked")
    if blocked:
        for text, label in (("low-complexity region", "gnomAD low-complexity flag"),
                            ("segmental duplication", "gnomAD segmental-duplication flag"),
                            ("gnomAD FILTER=", "gnomAD FILTER"),
                            ("paralogue/pseudogene homology", "gene on the homology list"),
                            ("low-penetrance", "ClinVar low-penetrance term")):
            if text in blocked:
                return label
        return "another frequency blocker"
    if d.get("an_below_minimum") is not None or d.get("an_missing"):
        return "allele number below the floor"
    if d.get("exception_match"):
        return "curated frequency exception"
    if criterion.get("met"):
        return "BS1 met"
    return "evaluated and not met (filtering allele frequency below the bar)"


def bs1_met_label(reason: str, criteria: str) -> str:
    """BS1 met is held at VUS either because the 2015 table does not call it alone or because other evidence is beside it."""
    if reason != "BS1 met":
        return reason
    return "BS1 met, alone" if criteria.split("&") == ["BS1"] else "BS1 met, with other criteria"


def blockers() -> None:
    """Why BA1 and BS1 did not fire for the ClinVar Benign variants of 2% or more that are not called benign-side.

    The per-variant table says which criteria fired and not why a frequency criterion did not, so those variants
    (the 238 new ones, and the 7,899 of the April set for comparison) are annotated again with the JSON output,
    which carries the details of every criterion, by the binary and with the flags of the benchmark. The reason
    BS1 gives is counted: BS1 has the lower bar, and the site-level and gene-level blockers withhold BA1 and BS1
    alike. The calls of this second annotation are checked against the table's, one by one.

    Writes manuscript/acmg/data/temporal_high_af_blockers_v37.csv.
    """
    import gzip
    import json
    import tempfile

    out_csv = ROOT / "manuscript/acmg/data/temporal_high_af_blockers_v37.csv"
    fastvep = ROOT / "target/release/fastvep"
    gff3 = ROOT / "test_data/organisms/human/Homo_sapiens.GRCh38.115.gff3"
    fasta = ROOT / "test_data/organisms/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    strata = {}
    with (SUBSET_DIR / "strata.tsv").open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            strata[(r["chrom"], r["pos"], r["ref"], r["alt"])] = r["stratum"]

    def uncalled(path: Path, only_new: bool) -> list[dict]:
        rows = []
        with gzip.open(path, "rt") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["clinvar"] != "Benign" or r["call"] in ("Likely_benign", "Benign"):
                    continue
                if not r["gnomad_popmax_af"] or float(r["gnomad_popmax_af"]) < 2e-2:
                    continue
                if only_new and strata[(r["chrom"], r["pos"], r["ref"], r["alt"])] != "new":
                    continue
                rows.append(r)
        return rows

    sets = (("April set", uncalled(BENCH / "manuscript_runs/v37_per_variant.tsv.gz", False)),
            ("new variants", uncalled(PER_VARIANT, True)))
    out = []
    with tempfile.TemporaryDirectory(dir=BENCH) as tmp:
        for name, rows in sets:
            vcf, js = Path(tmp) / "in.vcf", Path(tmp) / "out.json"
            with vcf.open("w") as f:
                f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                for r in rows:
                    f.write(f"{r['chrom']}\t{r['pos']}\t.\t{r['ref']}\t{r['alt']}\t.\t.\t.\n")
            subprocess.run([str(fastvep), "annotate", "-i", str(vcf), "-o", str(js), "--gff3", str(gff3), "--fasta", str(fasta),
                            "--sa-dir", str(BENCH / "sa_db"), "--acmg", "--pick", "--output-format", "json"],
                           check=True, stderr=subprocess.DEVNULL)
            records = json.loads(js.read_text())
            if len(records) != len(rows):
                sys.exit(f"{name}: {len(rows)} variants in, {len(records)} records out")
            reasons = Counter()
            for rec, r in zip(records, rows):
                tcs = rec["transcript_consequences"]
                if len(tcs) != 1:
                    sys.exit(f"{name}: --pick returned {len(tcs)} transcripts for {r['chrom']}:{r['pos']}")
                acmg = tcs[0]["acmg"]
                if JSON_CLASS[acmg["classification"]] != r["call"]:
                    sys.exit(f"{name}: the call of {r['chrom']}:{r['pos']} differs from the per-variant table's")
                bs1 = next(c for c in acmg["criteria"] if c["code"] == "BS1")
                reasons[bs1_met_label(bs1_reason(bs1), r["criteria"])] += 1
            for reason, k in reasons.most_common():
                out.append({"set": name, "uncalled_benign_variants_at_2e-2_or_more": len(rows), "reason_bs1_did_not_fire": reason,
                            "variants": k, "share": k / len(rows)})
            js.unlink()
    with out_csv.open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(f"wrote {out_csv.name} ({len(out)} rows)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["prepare", "run", "table", "blockers", "all"])
    args = ap.parse_args()
    if args.step in ("prepare", "all"):
        prepare()
    if args.step in ("run", "all"):
        run()
    if args.step in ("table", "all"):
        table()
    if args.step == "blockers":
        blockers()


if __name__ == "__main__":
    main()
