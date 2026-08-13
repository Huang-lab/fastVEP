#!/usr/bin/env python3
"""Build a per-gene frequency-threshold table from the ClinGen CSpec Registry.

ClinGen Variant Curation Expert Panels publish gene- and disease-specific
specifications of the ACMG/AMP criteria, and the ones that matter most to
fastVEP are the frequency bars: BA1, BS1 and PM2. A VCEP has done the work
fastVEP's global defaults cannot - it has looked up the disorder's prevalence,
its penetrance and its allelic heterogeneity, and derived a bar for that gene.
Where a VCEP has published one, it outranks anything measured across all genes.

This script fetches those specifications from the CSpec Registry's JSON-LD API
and emits two files:

  * a TOML fragment of `[gene_overrides.<GENE>]` blocks, loadable with
    `--acmg-config`;
  * a TSV audit table carrying, for every gene and criterion, the verbatim
    sentence the number was read out of, so a curator can check the parse
    rather than trust it.

Nothing here guesses. A threshold is emitted only when the source sentence
states it unambiguously, and every rejection is recorded in the audit table
with the reason. See `--help` and the module docstring in
`extract_threshold` for what "unambiguously" means.

Usage:

    python3 build_vcep_thresholds.py \\
        --cache data/benchmark/sa_sources/cspec \\
        --out-toml analysis/acmg_benchmark/data/vcep_thresholds.toml \\
        --out-audit analysis/acmg_benchmark/data/vcep_thresholds_audit.tsv

The cache directory is reused across runs; delete it to re-fetch.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from typing import Iterable

API = "https://cspec.genome.network/cspec/api"

# Criteria whose specifications this script understands. BA1 and BS1 are
# "frequency at or above X is benign evidence"; PM2 is "frequency at or below X
# is pathogenic evidence". The direction matters for which of several numbers in
# a sentence is the threshold.
BENIGN_CRITERIA = ("BA1", "BS1")
PATHOGENIC_CRITERIA = ("PM2",)

# HPO mode-of-inheritance terms, mapped onto the vocabulary AcmgConfig uses.
INHERITANCE = {
    "autosomal dominant inheritance": "AD",
    "autosomal recessive inheritance": "AR",
    "x-linked dominant inheritance": "AD",
    "x-linked recessive inheritance": "AR",
    "x-linked inheritance": "AR",
    "semidominant inheritance": "AD_AR",
    "mitochondrial inheritance": "MT",
}


@dataclass
class Row:
    """One extracted (or rejected) threshold, for the audit table."""

    gene: str
    criterion: str
    strength: str
    value: float | None
    status: str
    note: str
    vcep: str
    spec: str
    spec_status: str
    diseases: str
    inheritance: str
    flag: str
    source_text: str


@dataclass
class GeneEntry:
    vcep: str = ""
    spec: str = ""
    spec_status: str = ""
    diseases: list[str] = field(default_factory=list)
    inheritance: set[str] = field(default_factory=set)
    thresholds: dict[str, float] = field(default_factory=dict)
    flags: dict[str, str] = field(default_factory=dict)


def fetch(url: str, cache_path: str) -> dict:
    """GET `url` as JSON, caching the raw bytes at `cache_path`."""
    if os.path.exists(cache_path) and os.path.getsize(cache_path) > 100:
        with open(cache_path, "rb") as handle:
            return json.loads(handle.read())
    with urllib.request.urlopen(url, timeout=60) as response:
        raw = response.read()
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    with open(cache_path, "wb") as handle:
        handle.write(raw)
    return json.loads(raw)


def strip_markdown(text: str) -> str:
    """Reduce a CSpec description to plain prose.

    The descriptions are Markdown with emphasis, footnote links and escaped
    punctuation, all of which sit between a number and its operator often
    enough to defeat a naive regex.
    """
    text = re.sub(r"\[[^\]]*\]\([^)]*\)", " ", text)  # links, incl. footnotes
    text = re.sub(r"<[^>]+>", " ", text)  # inline HTML
    text = text.replace("\\", "")
    text = re.sub(r"[*_`]+", "", text)
    text = text.replace(" ", " ")
    return re.sub(r"\s+", " ", text).strip()


# A frequency written as a percentage, a decimal, or in scientific notation.
#
# The percentage alternative comes FIRST and that ordering is load-bearing:
# Python's alternation is leftmost-first, not longest-match, so putting the
# decimal pattern first makes "0.1%" parse as 0.1 rather than 0.001. That is a
# thousand-fold error in the benign direction, and it is silent.
NUMBER = r"(\d*\.?\d+(?:[eE][-+]?\d+)?%|\d*\.\d+(?:[eE][-+]?\d+)?|\d+(?:\.\d+)?[eE][-+]?\d+)"
GE = r"(?:≥|>=|&ge;|greater (?:than )?or equal to|at or above|of at least)"
GT = r"(?:>|&gt;|greater than|above|exceeds)"
LE = r"(?:≤|<=|&le;|less than or equal to|at or below|no more than|no greater than)"
LT = r"(?:<|&lt;|less than|below)"
OPERATOR = re.compile(rf"({GE}|{GT}|{LE}|{LT})\s*(?:than\s+)?{NUMBER}", re.IGNORECASE)
# "MAF cutoff of 0.2%" states a bar without an inequality. Narrow on purpose:
# the noun has to be a frequency, so a "prevalence cutoff" cannot match.
BARE_CUTOFF = re.compile(
    rf"(?:MAF|minor allele frequency|allele frequency|AF|FAF)\s*"
    rf"(?:cut[- ]?off|threshold)\s*(?:of|:|is|=)?\s*{NUMBER}",
    re.IGNORECASE,
)
PERCENT_GLOSS = re.compile(rf"{NUMBER}\s*\(\s*(?:FAF|MAF|AF)?\s*[<>≤≥=]*\s*([\d.]+)\s*%\s*\)")

# Words that make a nearby number an allele frequency rather than something
# else the specification happens to quantify - a disease prevalence, an allele
# count, a read depth, a share of the variant spectrum.
FREQUENCY_CONTEXT = re.compile(
    r"(?:MAF|minor allele frequency|allele frequency|frequency|FAF|popmax|grpmax|"
    r"gnomAD|population database|filtering allele)",
    re.IGNORECASE,
)
CONTEXT_WINDOW = 80
# Descriptions shorter than this are the whole of what the VCEP wrote for a
# frequency criterion, so the criterion itself supplies the context. ETHE1's
# entire BA1 description is ">0.001 (>0.1%)".
SELF_EVIDENT_LENGTH = 60

# Scientific notation written with a superscript exponent - "8 x 10<sup>-3</sup>"
# - loses its minus sign and its binding when the markup is stripped, leaving
# "8 x 10 -3" for a regex to misread as the number 3 or 10. Refuse the whole
# description rather than parse around it.
SUPERSCRIPT_NOTATION = re.compile(r"\d\s*[x×]\s*10\s*-?\s*\d", re.IGNORECASE)

# Ranges an allele-frequency bar for each criterion must fall in to be
# believable. These are not the criteria's own semantics; they are a net for
# mis-parses, and every catch so far has been one. BA1 above 50 % is not a
# frequency bar, and a PM2 bar of 2 % came from "the most frequent pathogenic
# variant accounts for no more than 2% of..." - a statement about variant
# spectrum that happens to sit beside the word frequency.
PLAUSIBLE = {"BA1": (1e-5, 0.5), "BS1": (1e-6, 0.1), "PM2": (0.0, 0.01)}

# fastVEP's global defaults, for flagging outliers. A VCEP bar looser than the
# default deserves a curator's eye before it ships, because it makes the gene
# easier to call benign than the tool's baseline - ABCA4's published BA1 of
# 0.163 is a real example, faithfully parsed. Flagging is not rejection: the
# panel's number is the authority, and the flag only says where to look.
DEFAULTS = {"BA1": 0.05, "BS1": 0.005, "PM2": 0.0001}


def outlier_flag(criterion: str, value: float) -> str:
    """A short note when a bar sits far from fastVEP's global default."""
    default = DEFAULTS[criterion]
    if criterion in BENIGN_CRITERIA and value > default:
        return f"looser than the {default:g} default"
    if criterion in PATHOGENIC_CRITERIA and value > default:
        return f"looser than the {default:g} default"
    if value > 0 and default / value >= 100:
        return f"{default / value:.0f}x tighter than the {default:g} default"
    return ""

# Specifications that state one bar for dominant and another for recessive
# disease. Real content, not a parse failure, but `gene_overrides` carries one
# number per gene so it cannot be represented.
INHERITANCE_SPLIT = re.compile(
    r"autosomal (?:dominant|recessive)|X-linked|for (?:AD|AR)\b", re.IGNORECASE
)
# PM2 specifications that ask for absence rather than a bar.
ABSENCE_ONLY = re.compile(r"\babsent\b", re.IGNORECASE)


def as_fraction(token: str) -> float:
    """Parse one frequency token into a fraction of 1."""
    if token.endswith("%"):
        return float(token[:-1]) / 100.0
    return float(token)


def has_frequency_context(plain: str, at: int, end: int) -> bool:
    """Whether the text around offsets `at`..`end` is talking about a frequency.

    Both directions matter: specifications write "Allele frequency ≥0.1%" and
    "≥0.000156 (0.0156%) GroupMax Filtering Allele Frequency" about equally
    often. A description short enough to be nothing but the bar carries its
    context in the criterion label instead.
    """
    if len(plain) <= SELF_EVIDENT_LENGTH:
        return True
    window = plain[max(0, at - CONTEXT_WINDOW) : end + CONTEXT_WINDOW]
    return FREQUENCY_CONTEXT.search(window) is not None


def extract_threshold(criterion: str, text: str) -> tuple[float | None, str, str]:
    """Read the frequency bar out of one criterion description.

    Returns `(value, status, note)`. `status` is `ok` when the text states a
    single unambiguous bar, and otherwise names why nothing was emitted.

    The rules are deliberately strict, because a wrong threshold silently
    changes calls in a clinical direction and is far worse than a missing one:

    * A number counts only when an inequality operator is attached to it, or
      when it is named as a cutoff. Specifications quantify plenty of other
      things in passing.
    * The number must sit in frequency-talking prose. Without that test,
      "the most frequent pathogenic variant accounts for no more than 2% of..."
      reads as a 2 % allele-frequency bar.
    * The operator has to point the right way for the criterion. BA1 and BS1
      fire *above* their bar and PM2 *below* its bar, so an operator pointing
      the other way means the sentence is about something else.
    * Where several qualifying numbers appear, they must agree, and a
      disagreement caused by a dominant/recessive split is reported as such
      rather than as a parse failure.
    * Where the specification glosses its own number as a percentage -
      "0.005 (0.5%)" - the two must agree. This catches transcription errors in
      the source and mis-parses here alike.
    """
    plain = strip_markdown(text)
    if not plain:
        return None, "empty", "no description text"
    if SUPERSCRIPT_NOTATION.search(plain):
        return None, "superscript-notation", "exponent lost when markup was stripped"

    for value_token, percent_token in PERCENT_GLOSS.findall(plain):
        value = as_fraction(value_token)
        gloss = float(percent_token) / 100.0
        if abs(value - gloss) > max(value, gloss) * 1e-6:
            return (
                None,
                "gloss-mismatch",
                f"stated {value_token} but glossed as {percent_token}%",
            )

    wants_upper_bound = criterion in PATHOGENIC_CRITERIA
    candidates: list[float] = []
    wrong_way = 0
    out_of_context = 0
    for match in OPERATOR.finditer(plain):
        operator, token = match.group(1), match.group(2)
        if not has_frequency_context(plain, match.start(), match.end()):
            out_of_context += 1
            continue
        is_upper_bound = operator[0] in "≤<" or operator.lower().startswith(
            ("less", "at or below", "no more", "no greater", "<=", "&le;", "&lt;", "below")
        )
        if is_upper_bound != wants_upper_bound:
            wrong_way += 1
            continue
        candidates.append(as_fraction(token))

    for match in BARE_CUTOFF.finditer(plain):
        candidates.append(as_fraction(match.group(1)))

    if not candidates:
        if criterion == "PM2" and ABSENCE_ONLY.search(plain):
            return None, "absence-required", "specifies absence from controls, not a bar"
        if wrong_way:
            return None, "wrong-direction", f"{wrong_way} bound(s), none the right way"
        if out_of_context:
            return None, "no-threshold", f"{out_of_context} number(s), none an allele frequency"
        return None, "no-threshold", "no inequality attached to a frequency"

    unique = sorted(set(candidates))
    low, high = PLAUSIBLE[criterion]
    implausible = [v for v in unique if not low <= v <= high]
    unique = [v for v in unique if low <= v <= high]
    if not unique:
        formatted = ", ".join(f"{v:g}" for v in implausible)
        return None, "implausible", f"outside the believable range for {criterion}: {formatted}"
    if len(unique) > 1:
        formatted = ", ".join(f"{v:g}" for v in unique)
        if INHERITANCE_SPLIT.search(plain):
            return None, "inheritance-conditioned", f"separate bars by inheritance: {formatted}"
        return None, "ambiguous", f"{len(unique)} disagreeing bars: {formatted}"
    note = ""
    if implausible:
        formatted = ", ".join(f"{v:g}" for v in implausible)
        note = f"ignored as out of range: {formatted}"
    return unique[0], "ok", note


def strongest_applicable(criterion_code: dict, criterion: str) -> tuple[str, str] | None:
    """The evidence strength this VCEP actually applies, and its description.

    A specification lists every strength and marks each Applicable or not. For
    BS1 some VCEPs apply both Strong and Supporting; fastVEP carries one number
    per gene, so the criterion's own default strength wins - `BS1` means the
    Strong row, and the Supporting row is a separate, weaker code fastVEP does
    not model per gene.
    """
    preferred = {"BA1": "Stand Alone", "BS1": "Strong", "PM2": None}[criterion]
    applicable = [
        (es.get("label", ""), es.get("description") or "")
        for es in criterion_code.get("evidenceStrengths", [])
        if (es.get("applicability") or "").lower().startswith("applicable")
    ]
    if not applicable:
        return None
    if preferred:
        for label, description in applicable:
            if label == preferred:
                return label, description
        return None
    return applicable[0]


def collect(cache_dir: str, released_only: bool) -> tuple[dict[str, GeneEntry], list[Row]]:
    registry = fetch(f"{API}/svis", os.path.join(cache_dir, "svis.json"))["data"]
    genes: dict[str, GeneEntry] = {}
    rows: list[Row] = []

    for spec in registry:
        status = spec.get("status", "")
        if released_only and status != "Released":
            continue
        vcep = (spec.get("affiliation") or {}).get("label", "")
        spec_id = spec.get("@id", "").rsplit("/", 1)[-1]

        for rule_set in spec.get("ruleSets", []):
            rule_set_id = rule_set["@id"].rsplit("/", 1)[-1]
            members = []
            for gene in rule_set.get("genes", []):
                diseases = [d.get("label", "") for d in gene.get("diseases", [])]
                inheritance = {
                    INHERITANCE[moi.get("@label", "").lower()]
                    for d in gene.get("diseases", [])
                    for moi in d.get("modeOfInheritance", [])
                    if moi.get("@label", "").lower() in INHERITANCE
                }
                members.append((gene.get("label", ""), diseases, inheritance))
            if not members:
                continue

            document = fetch(
                f"{API}/RuleSet/id/{rule_set_id}",
                os.path.join(cache_dir, f"{rule_set_id}.json"),
            )
            by_label = {c.get("label"): c for c in document.get("criteriaCodes", [])}

            for criterion in BENIGN_CRITERIA + PATHOGENIC_CRITERIA:
                code = by_label.get(criterion)
                chosen = strongest_applicable(code, criterion) if code else None
                if chosen is None:
                    for symbol, diseases, inheritance in members:
                        rows.append(
                            Row(symbol, criterion, "", None, "not-specified",
                                "VCEP marks this criterion Not Applicable, or omits it",
                                vcep, spec_id, status, "|".join(diseases),
                                "|".join(sorted(inheritance)), "", "")
                        )
                    continue
                strength, description = chosen
                value, parse_status, note = extract_threshold(criterion, description)
                excerpt = strip_markdown(description)[:400]
                for symbol, diseases, inheritance in members:
                    rows.append(
                        Row(symbol, criterion, strength, value, parse_status, note,
                            vcep, spec_id, status, "|".join(diseases),
                            "|".join(sorted(inheritance)),
                            outlier_flag(criterion, value) if value is not None else "",
                            excerpt)
                    )
                    if value is None:
                        continue
                    entry = genes.setdefault(symbol, GeneEntry())
                    entry.vcep = entry.vcep or vcep
                    entry.spec = entry.spec or spec_id
                    entry.spec_status = entry.spec_status or status
                    for disease in diseases:
                        if disease not in entry.diseases:
                            entry.diseases.append(disease)
                    entry.inheritance |= inheritance
                    entry.thresholds[criterion] = value
                    flag = outlier_flag(criterion, value)
                    if flag:
                        entry.flags[criterion] = flag

    return genes, rows


def misordered_genes(genes: dict[str, GeneEntry]) -> dict[str, str]:
    """Genes whose extracted bars do not nest the way the criteria require.

    BA1 ≥ BS1 ≥ PM2 holds by the meaning of the criteria: a frequency high
    enough to be standalone-benign is at least as high as one that is strong
    benign evidence, which is at least as high as one rare enough to support
    pathogenicity. A published specification will not violate that, so a
    violation is this script mis-reading one of them.

    This is the check that earns its keep. It is what caught PM2 bars of 2 %
    read out of "the most frequent pathogenic variant accounts for no more than
    2% of..." across eight cardiomyopathy genes, and BA1 bars of 4 and 6 read
    out of superscript exponents. Both are now blocked upstream; the check
    stays because the next mis-parse will be one nobody predicted.
    """
    bad: dict[str, str] = {}
    for symbol, entry in genes.items():
        ba1 = entry.thresholds.get("BA1")
        bs1 = entry.thresholds.get("BS1")
        pm2 = entry.thresholds.get("PM2")
        for higher, lower, label in (
            (ba1, bs1, "BA1 < BS1"),
            (bs1, pm2, "BS1 < PM2"),
            (ba1, pm2, "BA1 < PM2"),
        ):
            if higher is not None and lower is not None and higher < lower:
                bad[symbol] = f"{label} ({higher:g} < {lower:g})"
                break
    return bad


def conflicting_genes(rows: Iterable[Row]) -> set[str]:
    """Genes where two specifications state different bars for one criterion.

    A gene can appear in more than one VCEP's scope - most often because two
    panels curate it for two different diseases - and their bars need not
    agree. `gene_overrides` is keyed by gene alone, so there is no honest way
    to pick one. These are dropped from the TOML and left in the audit table.
    """
    seen: dict[tuple[str, str], set[float]] = {}
    for row in rows:
        if row.value is None:
            continue
        seen.setdefault((row.gene, row.criterion), set()).add(row.value)
    return {gene for (gene, _), values in seen.items() if len(values) > 1}


def write_toml(path: str, genes: dict[str, GeneEntry], dropped: set[str]) -> int:
    lines = [
        "# ClinGen VCEP frequency thresholds, per gene.",
        "#",
        "# GENERATED FILE - do not edit by hand.",
        "# Regenerate with analysis/acmg_benchmark/scripts/build_vcep_thresholds.py",
        "#",
        "# Every value below is a published ClinGen Variant Curation Expert Panel",
        "# specification of BA1, BS1 or PM2 for that gene, read from the CSpec",
        "# Registry. Where a VCEP has set a bar it outranks fastVEP's global",
        "# default, which is measured across all genes and cannot know a specific",
        "# disorder's prevalence or penetrance.",
        "#",
        "# Load with:  fastvep annotate --acmg --acmg-config <this file>",
        "#",
        "# What is NOT here, and why:",
        "#  - genes whose specification states no numeric bar for a criterion, or",
        "#    states it in prose this generator will not guess at;",
        "#  - genes curated by two panels that disagree, since `gene_overrides` is",
        "#    keyed by gene alone and picking one would be arbitrary;",
        "#  - prevalence, penetrance and age of onset, which no VCEP publishes in",
        "#    machine-readable form. Those remain a curation gap.",
        "#",
        "# The `disorders` blocks are informational: they record which MONDO",
        "# disorders the thresholds were specified for, and their inheritance.",
        "# fastVEP has no active-disorder selection yet, so the gene-level values",
        "# are what the classifier reads.",
        "",
    ]
    written = 0
    for symbol in sorted(genes):
        entry = genes[symbol]
        if symbol in dropped or not entry.thresholds:
            continue
        written += 1
        lines.append(f"# {entry.vcep} ({entry.spec}, {entry.spec_status})")
        for criterion, flag in sorted(entry.flags.items()):
            lines.append(f"# REVIEW: {criterion} is {flag}")
        lines.append(f"[gene_overrides.{symbol}]")
        for criterion, key in (
            ("BA1", "ba1_af_threshold"),
            ("BS1", "bs1_af_threshold"),
            ("PM2", "pm2_af_threshold"),
        ):
            if criterion in entry.thresholds:
                lines.append(f"{key} = {entry.thresholds[criterion]:g}")
        for disease in entry.diseases:
            key = disease.replace('"', "")
            lines.append(f'[gene_overrides.{symbol}.disorders."{key}"]')
            if len(entry.inheritance) == 1:
                lines.append(f'inheritance = "{next(iter(entry.inheritance))}"')
            for criterion, field_name in (("BS1", "bs1_af_threshold"), ("PM2", "pm2_af_threshold")):
                if criterion in entry.thresholds:
                    lines.append(f"{field_name} = {entry.thresholds[criterion]:g}")
        lines.append("")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as handle:
        handle.write("\n".join(lines))
    return written


def write_audit(
    path: str, rows: list[Row], dropped: set[str], misordered: dict[str, str]
) -> None:
    header = [
        "gene", "criterion", "strength", "value", "status", "note", "vcep",
        "spec", "spec_status", "diseases", "inheritance", "flag", "source_text",
    ]
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in sorted(rows, key=lambda r: (r.gene, r.criterion, r.spec)):
            status, note = row.status, row.note
            if row.value is not None and row.gene in misordered:
                status, note = "misordered", misordered[row.gene]
            elif row.value is not None and row.gene in dropped:
                status = "conflicting-panels"
            handle.write(
                "\t".join(
                    [
                        row.gene, row.criterion, row.strength,
                        f"{row.value:g}" if row.value is not None else "",
                        status, note, row.vcep, row.spec, row.spec_status,
                        row.diseases, row.inheritance, row.flag,
                        row.source_text.replace("\t", " "),
                    ]
                )
                + "\n"
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cache", required=True, help="directory for fetched JSON")
    parser.add_argument("--out-toml", required=True)
    parser.add_argument("--out-audit", required=True)
    parser.add_argument(
        "--include-unreleased",
        action="store_true",
        help="also read specifications still in preparation (default: Released only)",
    )
    args = parser.parse_args()

    try:
        genes, rows = collect(args.cache, released_only=not args.include_unreleased)
    except urllib.error.URLError as error:
        print(f"CSpec Registry unreachable: {error}", file=sys.stderr)
        return 1

    misordered = misordered_genes(genes)
    dropped = conflicting_genes(rows) | set(misordered)
    written = write_toml(args.out_toml, genes, dropped)
    write_audit(args.out_audit, rows, dropped, misordered)

    counts: dict[str, int] = {}
    for row in rows:
        counts[row.status] = counts.get(row.status, 0) + 1
    print(f"genes with at least one usable threshold: {written}")
    print(f"genes dropped for disagreeing panels:     {len(dropped) - len(misordered)}")
    print(f"genes dropped for out-of-order bars:      {len(misordered)}")
    for symbol, reason in sorted(misordered.items()):
        print(f"    {symbol}: {reason}")
    print("per-criterion parse outcomes:")
    for status in sorted(counts):
        print(f"  {status:<16} {counts[status]}")
    print(f"wrote {args.out_toml}")
    print(f"wrote {args.out_audit}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
