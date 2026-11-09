#!/usr/bin/env python3
"""Which of the 28 ACMG/AMP criteria the benchmark can apply, and what the engine reports for each of them on a sample of it.

The manuscript says that the benchmark, a VCF of variants with no genotype, pedigree or phenotype of anyone, can
apply 15 of the 28 criteria, and names the other 13 with the reason each cannot be applied. That is a claim about
what the engine does with its input, and the annotated VCFs of the runs cannot show it: they list only the criteria
that were met. The JSON output lists every criterion that the engine considered, with whether it could be evaluated
(`evaluated`) and whether it was met (`met`), so this script annotates every 40th record of the benchmark to JSON,
counts the two for each criterion, and checks the counts against the status that the manuscript states.

The status, what each criterion needs and how the engine reads it where it can are written here and not read from
the engine, because they describe the criteria and not the code. What the engine reports is read from its output,
and the script stops, without writing, if the two disagree:

  - the 15 criteria of the benchmark are evaluated for some records and met for some;
  - no other criterion is met, and none but PM3 is evaluated for any record (PM3 is evaluated in genes that are not
    read as recessive, where it is not applicable, and not evaluated in the others, where it needs the genotype
    of a proband);
  - PP5 and BP6 do not appear, because they are off by default.

Writes `manuscript/acmg/data/criteria_status.csv`, one row for each of the 28 criteria.

    18_criteria_status.py [--every 40]

Needs `FASTVEP_V38_BIN` (the binary of the shipped classifier) and the inputs of `run_benchmark.sh` in `data/benchmark/`.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
BENCH = ROOT / "data/benchmark"
OUT = ROOT / "manuscript/acmg/data/criteria_status.csv"
GFF3 = ROOT / "test_data/organisms/human/Homo_sapiens.GRCh38.115.gff3"
FASTA = ROOT / "test_data/organisms/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

VARIANT = "from the variant"
GENOTYPES = "needs genotypes of a proband"
FUNCTIONAL = "needs functional assays"
CASE = "needs case-level data"
OFF = "off by default"

# code: (group, what it needs, how the engine reads it where it can)
CRITERIA: dict[str, tuple[str, str, str]] = {
    "PVS1": (VARIANT, "consequence; gene constraint; ClinGen gene-disease validity and mechanism", ""),
    "PS1": (VARIANT, "ClinVar protein and splice index", ""),
    "PS2": (GENOTYPES, "de novo status with confirmed parentage", "a trio VCF (--proband --mother --father)"),
    "PS3": (FUNCTIONAL, "functional studies showing damage", "a curated file (--functional-evidence)"),
    "PS4": (CASE, "prevalence in affected individuals against controls", "not read; a proxy from ClinVar review stars is invalid and off"),
    "PM1": (VARIANT, "ClinVar protein index; ClinGen gene-disease validity", ""),
    "PM2": (VARIANT, "gnomAD; ClinVar frequencies; inheritance", ""),
    "PM3": (GENOTYPES, "in trans with a pathogenic variant for a recessive disorder", "the proband's genotype and its other variants in the gene"),
    "PM4": (VARIANT, "consequence", ""),
    "PM5": (VARIANT, "ClinVar protein index", ""),
    "PM6": (GENOTYPES, "assumed de novo, parentage unconfirmed", "a trio VCF (--proband --mother --father)"),
    "PP1": (CASE, "co-segregation in a family", "not read"),
    "PP2": (VARIANT, "gene constraint; ClinGen gene-disease validity", ""),
    "PP3": (VARIANT, "REVEL; SpliceAI", ""),
    "PP4": (CASE, "a phenotype specific to the gene's disease", "not read"),
    "PP5": (OFF, "a reputable source reports the variant as pathogenic", "use_pp5_bp6 in the configuration"),
    "BA1": (VARIANT, "gnomAD", ""),
    "BS1": (VARIANT, "gnomAD", ""),
    "BS2": (VARIANT, "gnomAD; inheritance", ""),
    "BS3": (FUNCTIONAL, "functional studies showing no damage", "a curated file (--functional-evidence)"),
    "BS4": (CASE, "lack of segregation in a family", "not read"),
    "BP1": (VARIANT, "gene constraint; ClinVar protein index", ""),
    "BP2": (GENOTYPES, "in trans with a pathogenic variant for a dominant disorder or in cis with one", "the proband's genotype and its other variants in the gene"),
    "BP3": (VARIANT, "consequence; RepeatMasker", ""),
    "BP4": (VARIANT, "REVEL; SpliceAI", ""),
    "BP5": (CASE, "an alternative molecular basis for the disease in the patient", "not read"),
    "BP6": (OFF, "a reputable source reports the variant as benign", "use_pp5_bp6 in the configuration"),
    "BP7": (VARIANT, "consequence; SpliceAI; PhyloP; exon position", ""),
}
ORDER = ["PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5", "PM6", "PP1", "PP2", "PP3", "PP4", "PP5",
         "BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"]
# The two criteria that a default run does not report, so that the engine's output cannot give their strength: Richards et al. (2015).
NOT_REPORTED_STRENGTH = {"PP5": "Supporting", "BP6": "Supporting"}


def sample_vcf(every: int, path: Path) -> int:
    """Every `every`th record of the benchmark VCF, with its header; the number of records written."""
    n = kept = 0
    with (BENCH / "clinvar_2star.vcf").open() as src, path.open("w") as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue
            n += 1
            if n % every == 0:
                dst.write(line)
                kept += 1
    return kept


def top_level_objects(path: Path):
    """The objects of the JSON array that `fastvep annotate --output-format json` writes, one at a time.

    The writer puts the brackets of the array on lines of their own and each variant between a line that is `{` and
    one that is `}` or `},`, so the file is read line by line and not parsed whole: it is some 800 MB for a sample of
    this size. The caller checks the count of objects against the number of records written, which would not hold if
    the layout were ever to change.
    """
    buf: list[str] = []
    with path.open() as f:
        for line in f:
            stripped = line.rstrip("\n")
            if stripped == "{":
                buf = [line]
            elif stripped in ("}", "},"):
                buf.append("}\n")
                yield json.loads("".join(buf))
                buf = []
            elif buf:
                buf.append(line)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--every", type=int, default=40, help="take every Nth record of the benchmark")
    args = ap.parse_args()
    binary = os.environ.get("FASTVEP_V38_BIN")
    if not binary or not Path(binary).exists():
        sys.exit("FASTVEP_V38_BIN must name the binary of the shipped classifier")

    seen, evaluated, met = Counter(), Counter(), Counter()
    default_strength: dict[str, str] = {}
    direction: dict[str, str] = {}
    variants = records = 0
    with tempfile.TemporaryDirectory() as tmp:
        vcf, out = Path(tmp) / "sample.vcf", Path(tmp) / "sample.json"
        written = sample_vcf(args.every, vcf)
        subprocess.run(
            [binary, "annotate", "-i", str(vcf), "-o", str(out), "--gff3", str(GFF3), "--fasta", str(FASTA),
             "--sa-dir", str(BENCH / os.environ.get("BENCH_SA_DB", "sa_db_complete")), "--acmg", "--pick", "--output-format", "json"],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        for variant in top_level_objects(out):
            variants += 1
            for consequence in variant.get("transcript_consequences") or []:
                result = consequence.get("acmg")
                if not result:
                    continue
                records += 1
                for c in result["criteria"]:
                    code = c["code"].split("_")[0]
                    seen[code] += 1
                    evaluated[code] += bool(c["evaluated"])
                    met[code] += bool(c["met"])
                    default_strength.setdefault(code, c["default_strength"])
                    direction.setdefault(code, c["direction"])
    if variants != written:
        sys.exit(f"the json holds {variants} variants and the sample {written}: the layout that top_level_objects reads has changed")

    problems = []
    if set(seen) | {"PP5", "BP6"} != set(ORDER):
        problems.append(f"the engine reports {sorted(seen)}, which is not the 26 criteria of a default run")
    for code in ORDER:
        group = CRITERIA[code][0]
        if group == VARIANT:
            if not evaluated[code] or not met[code]:
                problems.append(f"{code} is declared evaluable from the variant but was evaluated {evaluated[code]} and met {met[code]} times")
        elif group == OFF:
            if seen[code]:
                problems.append(f"{code} is declared off by default but appears {seen[code]} times")
        else:
            if met[code]:
                problems.append(f"{code} is declared not applicable to the benchmark but was met {met[code]} times")
            if code != "PM3" and evaluated[code]:
                problems.append(f"{code} is declared not evaluable in the benchmark but was evaluated {evaluated[code]} times")
    if not 0 < evaluated["PM3"] < records:
        problems.append(f"PM3 should be evaluated for some of the {records} records and not for all, and was for {evaluated['PM3']}")
    if problems:
        sys.exit("the engine's report disagrees with the status stated here, nothing written:\n  " + "\n  ".join(problems))

    rows = [
        {"criterion": code, "direction": direction.get(code, "Pathogenic" if code.startswith("P") else "Benign").lower(),
         "default_strength": default_strength.get(code) or NOT_REPORTED_STRENGTH[code], "group": CRITERIA[code][0], "needs": CRITERIA[code][1],
         "how_the_engine_reads_it": CRITERIA[code][2], "records": seen[code], "evaluated": evaluated[code], "met": met[code]}
        for code in ORDER
    ]
    pd.DataFrame(rows).to_csv(OUT, index=False, lineterminator="\n")
    print(f"wrote {OUT.name} ({len(rows)} rows) from {records} records of {variants} variants (every {args.every}th of the benchmark)")


if __name__ == "__main__":
    main()
