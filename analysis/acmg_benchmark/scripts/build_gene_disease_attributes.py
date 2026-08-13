#!/usr/bin/env python3
"""Assemble the per-gene-disease attribute table the frequency criteria need.

BS2 asks whether a disorder has "full penetrance expected at an early age".
BA1, BS1 and PM2 all ask what frequency is "greater than expected for the
disorder". Both questions are about the *disorder*, and fastVEP answers them
today with one global number each, because it has no table of disorder
attributes to consult.

This builds the fetchable part of that table, and marks the rest as missing
rather than inventing it.

Sources, all CC-BY and all public:

  * Orphanet `en_product6.xml`       gene to disorder, with the association type
  * Orphanet `en_product9_prev.xml`  prevalence, as a class and where stated a value
  * Orphanet `en_product9_ages.xml`  average age of onset, and inheritance
  * `vcep_thresholds_audit.tsv`      published ClinGen VCEP frequency bars, from
                                     build_vcep_thresholds.py

What is *not* sourced, and is emitted blank:

  * **penetrance** - Orphanet has no penetrance field, and no other public
    resource publishes it per gene-disease pair in machine-readable form. This
    is the field BS2 most wants, and it needs a clinician.
  * **allelic contribution** and **genetic heterogeneity** - the other two
    inputs to a Whiffin 2017 maximum-credible-AF calculation. VCEPs compute
    them per gene, but publish only the resulting bar, which is why the VCEP
    columns here carry the bar rather than its ingredients.

Usage:

    python3 build_gene_disease_attributes.py \\
        --orphanet-dir data/benchmark/sa_sources/orphanet \\
        --vcep-audit analysis/acmg_benchmark/data/vcep_thresholds_audit.tsv \\
        --out analysis/acmg_benchmark/data/gene_disease_attributes.tsv

Download the three XML files into `--orphanet-dir` first:

    for f in en_product6 en_product9_prev en_product9_ages; do
      curl -L -o <dir>/$f.xml https://www.orphadata.com/data/xml/$f.xml
    done
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

# Only germline disease-causing associations. Orphanet also records modifiers,
# susceptibility factors and candidate genes, none of which is the "this gene
# causes this disorder" relation the frequency criteria are reasoning about.
CAUSAL_ASSOCIATIONS = {
    "Disease-causing germline mutation(s) in",
    "Disease-causing germline mutation(s) (loss of function) in",
    "Disease-causing germline mutation(s) (gain of function) in",
}

# Orphanet inheritance terms, mapped onto the vocabulary AcmgConfig uses.
INHERITANCE = {
    "Autosomal dominant": "AD",
    "Autosomal recessive": "AR",
    "X-linked dominant": "AD",
    "X-linked recessive": "AR",
    "X-linked": "AR",
    "Mitochondrial inheritance": "MT",
    "Semi-dominant": "AD_AR",
}

# Orphanet's onset classes, ordered youngest first. BS2's precondition is
# "early age", so the earliest onset a disorder is recorded with is the one
# that matters - a disorder that can present neonatally is not late-onset just
# because it can also present in adulthood.
ONSET_ORDER = [
    "Antenatal",
    "Neonatal",
    "Infancy",
    "Childhood",
    "Adolescent",
    "Adult",
    "Elderly",
]
# Neither of these names an age. "All ages" *spans* the scale rather than
# sitting on it, so appending it to the ordering above would make a disorder
# recorded as {Adult, All ages} come out as adult-onset - understating exactly
# the thing BS2 asks about. Both are excluded from the earliest-onset choice
# and reported only when nothing specific was recorded.
ONSET_UNRANKED = ("All ages", "No data available")

# Orphanet prevalence classes, as an upper bound on the fraction of individuals
# affected. "1-9 / 100 000" becomes 9e-5. These feed BS2's prevalence bar,
# which is a *maximum credible* prevalence, so the top of the class is the
# right end to take.
PREVALENCE_CLASS_MAX = {
    ">1 / 1000": 1e-2,
    "6-9 / 10 000": 9e-4,
    "1-5 / 10 000": 5e-4,
    "1-9 / 100 000": 9e-5,
    "1-9 / 1 000 000": 9e-6,
    "<1 / 1 000 000": 1e-6,
    "Unknown": None,
    "Not yet documented": None,
}


def text_of(node, path: str) -> str:
    found = node.find(path)
    return (found.text or "").strip() if found is not None else ""


def parse_disorders(path: str, handler) -> None:
    """Stream `<Disorder>` elements through `handler`, freeing each as we go."""
    for _, element in ET.iterparse(path, events=("end",)):
        if element.tag != "Disorder":
            continue
        handler(element)
        element.clear()


def load_gene_associations(path: str) -> dict[str, list[tuple[str, str, str]]]:
    """OrphaCode -> [(gene symbol, association type, disorder name)]."""
    associations: dict[str, list[tuple[str, str, str]]] = defaultdict(list)

    def handle(disorder) -> None:
        code = text_of(disorder, "OrphaCode")
        name = text_of(disorder, "Name")
        for link in disorder.iter("DisorderGeneAssociation"):
            kind = text_of(link, "DisorderGeneAssociationType/Name")
            if kind not in CAUSAL_ASSOCIATIONS:
                continue
            symbol = text_of(link, "Gene/Symbol")
            if symbol:
                associations[code].append((symbol, kind, name))

    parse_disorders(path, handle)
    return associations


def load_prevalence(path: str) -> dict[str, dict[str, str]]:
    """OrphaCode -> the worldwide point prevalence, preferring validated rows.

    Orphanet records several prevalence estimates per disorder - by geography,
    by type (point, birth, lifetime, cases/families) and by validation status.
    A single row has to be picked, and the one that answers "how common is this
    disorder" is the validated worldwide point prevalence. Where that is
    absent, any validated point prevalence is taken and its geography recorded,
    so a reader can see the estimate is regional.
    """
    result: dict[str, dict[str, str]] = {}

    def handle(disorder) -> None:
        code = text_of(disorder, "OrphaCode")
        best = None
        best_rank = 99
        for row in disorder.iter("Prevalence"):
            kind = text_of(row, "PrevalenceType/Name")
            if kind not in ("Point prevalence", "Prevalence at birth"):
                continue
            geography = text_of(row, "PrevalenceGeographic/Name")
            validated = text_of(row, "PrevalenceValidationStatus/Name") == "Validated"
            rank = (
                0 if (kind == "Point prevalence" and geography == "Worldwide" and validated)
                else 1 if (kind == "Point prevalence" and validated)
                else 2 if validated
                else 3
            )
            if rank < best_rank:
                best_rank, best = rank, (row, kind, geography)
        if best is None:
            return
        row, kind, geography = best
        klass = text_of(row, "PrevalenceClass/Name")
        upper = PREVALENCE_CLASS_MAX.get(klass)
        value = text_of(row, "ValMoy")
        result[code] = {
            "prevalence_class": klass,
            "prevalence_max": f"{upper:g}" if upper is not None else "",
            "prevalence_type": kind,
            "prevalence_geographic": geography,
            "prevalence_value": value if value not in ("", "0.0") else "",
        }

    parse_disorders(path, handle)
    return result


def load_onset_and_inheritance(path: str) -> dict[str, dict[str, str]]:
    """OrphaCode -> earliest recorded onset, all onsets, and inheritance."""
    result: dict[str, dict[str, str]] = {}

    def handle(disorder) -> None:
        code = text_of(disorder, "OrphaCode")
        onsets = [text_of(o, "Name") for o in disorder.iter("AverageAgeOfOnset")]
        onsets = [o for o in onsets if o]
        ranked = [o for o in onsets if o in ONSET_ORDER]
        if ranked:
            earliest = min(ranked, key=ONSET_ORDER.index)
        else:
            spanning = [o for o in onsets if o in ONSET_UNRANKED]
            earliest = spanning[0] if spanning else ""
        modes = [text_of(m, "Name") for m in disorder.iter("TypeOfInheritance")]
        mapped = sorted({INHERITANCE[m] for m in modes if m in INHERITANCE})
        result[code] = {
            "onset_earliest": earliest,
            "onset_all": "|".join(onsets),
            "inheritance": "|".join(mapped),
            "inheritance_raw": "|".join(m for m in modes if m),
        }

    parse_disorders(path, handle)
    return result


def load_vcep(path: str) -> dict[str, dict[str, str]]:
    """Gene -> its usable published VCEP bars, from the audit table."""
    if not os.path.exists(path):
        return {}
    bars: dict[str, dict[str, str]] = defaultdict(dict)
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["status"] != "ok":
                continue
            bars[row["gene"]][f"vcep_{row['criterion'].lower()}"] = row["value"]
            bars[row["gene"]]["vcep_panel"] = row["vcep"]
            bars[row["gene"]]["vcep_spec"] = row["spec"]
    return bars


COLUMNS = [
    "gene", "orpha_code", "disorder", "association", "inheritance",
    "inheritance_raw", "onset_earliest", "onset_all", "prevalence_class",
    "prevalence_max", "prevalence_type", "prevalence_geographic",
    "prevalence_value", "penetrance", "vcep_ba1", "vcep_bs1", "vcep_pm2",
    "vcep_panel", "vcep_spec",
]


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--orphanet-dir", required=True)
    parser.add_argument("--vcep-audit", default="")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    paths = {
        name: os.path.join(args.orphanet_dir, f"{name}.xml")
        for name in ("en_product6", "en_product9_prev", "en_product9_ages")
    }
    missing = [p for p in paths.values() if not os.path.exists(p)]
    if missing:
        print("missing Orphanet XML: " + ", ".join(missing), file=sys.stderr)
        print("see this script's docstring for the download commands", file=sys.stderr)
        return 1

    associations = load_gene_associations(paths["en_product6"])
    prevalence = load_prevalence(paths["en_product9_prev"])
    natural_history = load_onset_and_inheritance(paths["en_product9_ages"])
    vcep = load_vcep(args.vcep_audit)

    rows = []
    for code, members in associations.items():
        for symbol, kind, disorder in members:
            row = {c: "" for c in COLUMNS}
            row.update(gene=symbol, orpha_code=code, disorder=disorder, association=kind)
            row.update(prevalence.get(code, {}))
            row.update(natural_history.get(code, {}))
            row.update(vcep.get(symbol, {}))
            # Not sourced anywhere machine-readable. Left blank on purpose;
            # see the module docstring.
            row["penetrance"] = ""
            rows.append(row)

    rows.sort(key=lambda r: (r["gene"], r["orpha_code"]))
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    genes = {r["gene"] for r in rows}
    filled = {c: sum(1 for r in rows if r[c]) for c in COLUMNS}
    print(f"{len(rows)} gene-disorder pairs across {len(genes)} genes -> {args.out}")
    print("column coverage:")
    for column in COLUMNS:
        share = filled[column] / len(rows) * 100 if rows else 0
        note = "  <- not sourced" if column == "penetrance" else ""
        print(f"  {column:<22} {filled[column]:>6}  {share:5.1f}%{note}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
