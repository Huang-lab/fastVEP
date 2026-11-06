#!/usr/bin/env python3
"""The calls that gnomAD's joint release turns into an error, or out of a pathogenic-side call, annotated one at a time through each database.

A run of the benchmark records a variant's call, the criteria that made it and the counts of gnomAD it read
(`FV_GNOMAD`), but not the filtering allele frequency that BS1 and the guard on a site that failed gnomAD's quality
control compare with their bar, and not the flags the record carried. A call that moves because one of those moved
cannot be explained from the run's own output.

This script selects the variants that matter most when a source changes,

  - those whose standing as a missed diagnosis (ClinVar pathogenic-side, called benign-side) or as a false-pathogenic
    call (ClinVar benign-side, called pathogenic-side) differs between two runs, and
  - those that ClinVar calls pathogenic-side and the first run calls Pathogenic or Likely pathogenic and the second
    does not,

annotates each of them through the databases of both runs with `--output-format json`, which carries the filtering
allele frequency, the flags and the reason a frequency criterion was not evaluated, and writes one row per variant and
database to `manuscript/acmg/data/v38_joint_errors.csv`. The call that the json gives is checked against the call in
the run's per-variant table, so a row cannot describe something other than what the run did.

It also takes the variants that lose PM2 between the exomes with the region flags withholding nothing (`v38_noflags`) and
the joint release, annotates them through the joint databases, and counts them by the reason PM2 is not met: a filter of
one callset that withholds it, or a frequency above the bar (`manuscript/acmg/data/v38_joint_pm2_lost.csv`).

    17_joint_boundary_cases.py [--exomes v38] [--joint v38_joint]

Needs `FASTVEP_V38_BIN` (the binary that made both runs), the per-variant tables of the two runs
(`12_per_variant_table.py`) and the joint databases in `data/benchmark/sa_db_joint/`.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
BENCH = ROOT / "data/benchmark"
RUNS = BENCH / "manuscript_runs"
SA_DB = BENCH / "sa_db"
SA_JOINT = BENCH / "sa_db_joint"
OUT = ROOT / "manuscript/acmg/data/v38_joint_errors.csv"
OUT_PM2 = ROOT / "manuscript/acmg/data/v38_joint_pm2_lost.csv"
GFF3 = ROOT / "test_data/organisms/human/Homo_sapiens.GRCh38.115.gff3"
FASTA = ROOT / "test_data/organisms/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
KEY = ["chrom", "pos", "ref", "alt"]

SIDE = {"Pathogenic": "P", "Likely_pathogenic": "P", "VUS": "V", "Likely_benign": "B", "Benign": "B"}
JSON_CALL = {"Pathogenic": "Pathogenic", "LikelyPathogenic": "Likely_pathogenic", "UncertainSignificance": "VUS",
             "LikelyBenign": "Likely_benign", "Benign": "Benign"}
# What the guard on a site that failed gnomAD's quality control writes into the rule that it overrides.
VETO_TEXT = "the population-frequency criteria could not be evaluated at all"


def error_kind(clinvar: str, call: str) -> str:
    if SIDE[clinvar] == "P" and SIDE[call] == "B":
        return "missed diagnosis"
    if SIDE[clinvar] == "B" and SIDE[call] == "P":
        return "false-pathogenic call"
    return ""


def standing(clinvar: str, call: str) -> str:
    """What a call is to the reference: an error, a pathogenic-side call on a pathogenic-side record, or neither."""
    return error_kind(clinvar, call) or ("pathogenic-side call" if SIDE[clinvar] == "P" and SIDE[call] == "P" else "")


def per_variant(run: str) -> pd.DataFrame:
    path = RUNS / f"{run}_per_variant.tsv.gz"
    if not path.exists():
        sys.exit(f"{path} does not exist: run 12_per_variant_table.py on {run}")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def staged_joint_dir(parent: Path) -> Path:
    """The databases of the exomes run with the gnomAD files of the joint release in place of the exomes'."""
    target = parent / "sa_joint"
    target.mkdir()
    for f in sorted(SA_DB.iterdir()):
        replacement = SA_JOINT / f.name
        if f.name.startswith("gnomad_chr") and f.suffix == ".osa2":
            if not replacement.exists():
                sys.exit(f"{replacement} does not exist")
            (target / f.name).symlink_to(replacement)
        else:
            (target / f.name).symlink_to(f)
    return target


def annotate(binary: str, vcf: Path, sa_dir: Path, out: Path) -> list[dict]:
    subprocess.run(
        [binary, "annotate", "-i", str(vcf), "-o", str(out), "--gff3", str(GFF3), "--fasta", str(FASTA),
         "--sa-dir", str(sa_dir), "--acmg", "--pick", "--output-format", "json"],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    return json.loads(out.read_text())


def records_of(variants: pd.DataFrame) -> tuple[list[str], list[str]]:
    """The benchmark's own VCF header and the records of `variants`, in their order, so that the json can be read back by position."""
    wanted = {tuple(r): i for i, r in enumerate(variants[KEY].itertuples(index=False, name=None))}
    header, found = [], {}
    with (BENCH / "clinvar_2star.vcf").open() as f:
        for line in f:
            if line.startswith("#"):
                header.append(line)
                continue
            fields = line.split("\t", 5)
            key = (fields[0], fields[1], fields[3], fields[4])
            if key in wanted:
                found[wanted[key]] = line
    if len(found) != len(wanted):
        sys.exit(f"{len(wanted) - len(found)} of the variants are not in the benchmark's VCF")
    return header, [found[i] for i in range(len(variants))]


def pm2_lost(binary: str, joint_dir: Path, tmp: Path, joint_run: str) -> pd.DataFrame:
    """The variants that lose PM2 between the exomes with the flags withholding nothing and the joint release, by the reason."""
    nf, jt = per_variant("v38_noflags"), per_variant(joint_run)
    m = nf[KEY + ["clinvar", "call", "criteria"]].merge(jt[KEY + ["call", "criteria"]], on=KEY, suffixes=("_nf", "_j"), validate="one_to_one")
    has = lambda s: s.str.contains("PM2_Supporting", regex=False)
    lost = m[has(m.criteria_nf) & ~has(m.criteria_j)].reset_index(drop=True)
    header, lines = records_of(lost)
    vcf = tmp / "pm2_lost.vcf"
    vcf.write_text("".join(header) + "".join(lines))
    rows = []
    for x, rec in zip(lost.itertuples(index=False), annotate(binary, vcf, joint_dir, tmp / "pm2_lost.json")):
        tc = rec["transcript_consequences"][0]
        g = tc.get("gnomad") or {}
        pm2 = next(c for c in tc["acmg"]["criteria"] if c["code"].startswith("PM2"))
        if pm2["evaluated"]:
            reason = "frequency above the bar"
        else:
            flags = [name for name, key in (("AC0", "filterAc0"), ("AS_VQSR", "filterVqsr"), ("InbreedingCoeff", "filterInbreeding")) if g.get(key)]
            reason = "filter " + " and ".join(flags) if flags else "withheld, no filter"
        rows.append({"reason": reason, "joint_count_positive": (g.get("allAc") or 0) > 0, "calls_moved": x.call_nf != x.call_j,
                     "pathogenic_side_calls_moved": x.call_nf != x.call_j and SIDE[x.clinvar] == "P"})
    out = pd.DataFrame(rows).groupby(["reason", "joint_count_positive"], as_index=False).agg(
        variants=("calls_moved", "size"), calls_moved=("calls_moved", "sum"), pathogenic_side_calls_moved=("pathogenic_side_calls_moved", "sum"))
    return out.sort_values("variants", ascending=False)


def blocked_by(criteria: list[dict]) -> str:
    """Why a frequency criterion was not evaluated at this record, in a few words; empty when none was withheld."""
    for c in criteria:
        if c["code"].split("_")[0] in ("BA1", "BS1", "BS2", "PM2") and not c["evaluated"]:
            why = c["details"].get("frequency_blocked") if c.get("details") else None
            if why is None:
                continue
            if why.startswith("gnomAD FILTER="):
                return why.split(";")[0].replace("gnomAD ", "")
            if "low-complexity" in why:
                return "low-complexity region"
            if "segmental duplication" in why:
                return "segmental duplication"
            return why[:60]
    return ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exomes", default="v38", help="the run on the exomes")
    ap.add_argument("--joint", default="v38_joint", help="the run on the joint release")
    args = ap.parse_args()
    binary = os.environ.get("FASTVEP_V38_BIN", "")
    if not binary or not Path(binary).is_file():
        sys.exit("FASTVEP_V38_BIN must name the binary that made the two runs")

    runs = {"exomes": per_variant(args.exomes), "joint": per_variant(args.joint)}
    cols = KEY + ["symbol", "clinvar", "call", "criteria", "consequence", "hgvsp"]
    both = runs["exomes"][cols].merge(runs["joint"][cols], on=KEY, suffixes=("_exomes", "_joint"), validate="one_to_one")
    both["standing_exomes"] = [standing(c, k) for c, k in zip(both.clinvar_exomes, both.call_exomes)]
    both["standing_joint"] = [standing(c, k) for c, k in zip(both.clinvar_joint, both.call_joint)]
    errors = [error_kind(c, k) for c, k in zip(both.clinvar_exomes, both.call_exomes)]
    errors_joint = [error_kind(c, k) for c, k in zip(both.clinvar_joint, both.call_joint)]
    lost = (both.standing_exomes == "pathogenic-side call") & (both.standing_joint != "pathogenic-side call")
    changed = both[(pd.Series(errors, index=both.index) != pd.Series(errors_joint, index=both.index)) | lost].reset_index(drop=True)
    changed["kind"] = [f"{a or 'no error'} -> {b or 'no error'}" if (a in ('missed diagnosis', 'false-pathogenic call') or b in ('missed diagnosis', 'false-pathogenic call'))
                       else "pathogenic-side call -> no pathogenic-side call"
                       for a, b in zip(changed.standing_exomes, changed.standing_joint)]
    print(f"{len(changed)} variants selected: {changed.kind.value_counts().to_dict()}")

    header, lines = records_of(changed)

    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        vcf = tmp / "variants.vcf"
        vcf.write_text("".join(header) + "".join(lines))
        sa = {"exomes": SA_DB, "joint": staged_joint_dir(tmp)}
        for database in ("exomes", "joint"):
            records = annotate(binary, vcf, sa[database], tmp / f"{database}.json")
            if len(records) != len(changed):
                sys.exit(f"{database}: {len(records)} records came back for {len(changed)} variants")
            for x, rec in zip(changed.itertuples(index=False), records):
                tc = rec["transcript_consequences"][0]
                acmg = tc["acmg"]
                g = tc.get("gnomad") or {}
                call = JSON_CALL[acmg["classification"]]
                in_table = getattr(x, f"call_{database}")
                if call != in_table:
                    sys.exit(f"{x.chrom}:{x.pos} {database}: the json says {call} and the run's table says {in_table}")
                by_code = {c["code"].split("_")[0]: c for c in acmg["criteria"]}
                rows.append({
                    "kind": x.kind,
                    "chrom": x.chrom, "pos": x.pos, "ref": x.ref, "alt": x.alt, "symbol": x.symbol_exomes,
                    "clinvar": x.clinvar_exomes, "database": database, "call": call,
                    "criteria": getattr(x, f"criteria_{database}"), "consequence": x.consequence_exomes,
                    "hgvsp": x.hgvsp_exomes.split(":")[-1],
                    "gnomad_af": g.get("allAf"), "gnomad_ac": g.get("allAc"), "gnomad_an": g.get("allAn"),
                    "faf95_max": g.get("faf95Max"), "bs1_bar": by_code["BS1"]["details"].get("af_threshold"),
                    "filter_ac0": g.get("filterAc0"), "filter_vqsr": g.get("filterVqsr"),
                    "filter_inbreeding": g.get("filterInbreeding"), "lcr": g.get("lcr"), "segdup": g.get("segdup"),
                    "bs1_evaluated": by_code["BS1"]["evaluated"], "bs1_met": by_code["BS1"]["met"],
                    "bs2_met": by_code["BS2"]["met"], "pm2_met": by_code["PM2"]["met"],
                    "frequency_withheld_by": blocked_by(acmg["criteria"]),
                    "capped_by_veto": VETO_TEXT in (acmg.get("triggered_rule") or ""),
                })
        lost = pm2_lost(binary, sa["joint"], tmp, args.joint)
    out = pd.DataFrame(rows).sort_values(["kind", "symbol", "pos", "database"], key=lambda s: s.map(str))
    out.to_csv(OUT, index=False, lineterminator="\n")
    print(f"wrote {OUT.relative_to(ROOT)} ({len(out)} rows)")
    lost.to_csv(OUT_PM2, index=False, lineterminator="\n")
    print(f"wrote {OUT_PM2.relative_to(ROOT)} ({len(lost)} rows)")


if __name__ == "__main__":
    main()
