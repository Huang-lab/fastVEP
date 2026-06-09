#!/usr/bin/env python3
"""Build the enriched medical-geneticist review TSV.

Joins the 112 opposite-direction discrepancies (where fastVEP
commits to a directional call that contradicts ClinVar 2-star+)
with full annotation fields from three sources:

- ClinVar VCF INFO: CLNHGVS, CLNDN, CLNDISDB, CLNREVSTAT, CLNSIG,
  CLNSIGSCV, CLNSIGCONF, ALLELEID, CLNVC, MC, GENEINFO, ORIGIN,
  CLNVI (RCV identifiers), AF_EXAC, AF_TGP, AF_ESP
- fastVEP CSQ (picked transcript, from `--output-format vcf`):
  HGVSc, HGVSp, BIOTYPE, EXON, INTRON, MANE_SELECT, CANONICAL,
  ENSP, CCDS, HGNC_ID, SIFT, PolyPhen, ACMG, ACMG_CRITERIA
- fastVEP JSON (picked transcript, from `--output-format json`)
  for the full per-variant score panel that the VCF CSQ doesn't
  carry:
    REVEL score
    SpliceAI dsAg / dsAl / dsDg / dsDl + max
    PhyloP, GERP
    gnomAD allAf, allAc, allAn, allHc, max-pop AF
    Gene-level: gnomAD constraints (pLI, LOEUF, misZ, synZ)
    Gene-level: ClinGen GDV / OMIM phenotypes
    Gene-level: ClinVar-protein pathogenic neighbours (count + list)
    ACMG triggered combination rule

Also computes a priority score for ranking:
  - 3-star > 2-star (panel reliability)
  - Pathogenic↔Benign (extreme reversal) > LP↔LB
  - More criteria fired = higher classifier confidence in disagreement

Outputs `data/benchmark/output_v7/discrepancies_for_md_review.tsv`
sorted highest-priority-first.

Usage:
  python3 04_build_md_review_table.py
"""

from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = ROOT / "data/benchmark/output_v7"
DISC = OUT_DIR / "discrepancies.tsv"
TRUTH = ROOT / "data/benchmark/clinvar_2star_truth.tsv"
CLINVAR_VCF = ROOT / "data/benchmark/clinvar_2star.vcf"
FASTVEP_VCF = OUT_DIR / "clinvar_2star.fastvep.vcf.gz"
# Override: when present, prefer a smaller VCF re-annotated with `--hgvs`
# so the review TSV gets populated HGVSc / HGVSp fields. The full-benchmark
# VCF.gz omits HGVS to save space; we recompute just for the ~100
# discrepancy variants.
FASTVEP_HGVS_VCF = OUT_DIR / "opposite_direction.fastvep.vcf"
# JSON re-annotation of the same 112 variants. Carries the full per-variant
# SA score panel (REVEL, SpliceAI components, PhyloP, gnomAD AFs and
# constraints, OMIM/ClinGen GDV phenotypes, ClinVar-protein pathogenic
# neighbours) that the VCF CSQ format doesn't expose.
FASTVEP_JSON = OUT_DIR / "opposite_direction.fastvep.json.gz"
OUT_TSV = OUT_DIR / "discrepancies_for_md_review.tsv"

# Fields we'll lift verbatim from each source. ClinVar field set covers
# everything an MD curator would consult before re-classifying. The fastVEP
# CSQ subset includes both transcript context (MANE_SELECT, CCDS, ENSP) and
# the prediction-tool fields needed to evaluate the criteria signature.
CLINVAR_FIELDS = [
    "ALLELEID", "CLNHGVS", "CLNDN", "CLNDISDB", "CLNREVSTAT",
    "CLNSIG", "CLNSIGCONF", "CLNSIGSCV", "CLNVC", "CLNVI",
    "GENEINFO", "MC", "ORIGIN", "AF_EXAC", "AF_TGP", "AF_ESP",
]
FASTVEP_FIELDS = [
    "HGVSc", "HGVSp", "BIOTYPE", "EXON", "INTRON", "MANE_SELECT",
    "CANONICAL", "ENSP", "CCDS", "HGNC_ID", "SIFT", "PolyPhen",
    "ACMG", "ACMG_CRITERIA",
]
# Per-variant score / annotation panel from the JSON output. Each
# entry is `(output_column_suffix, JSONPath-style accessor)`. The
# accessor is a list/tuple where strings index dicts and the special
# value (`SPLICEAI_MAX`, `MAX_POP_AF`, `PROTEIN_VARIANT_COUNT`,
# `JOIN`) trigger small computed transforms.
FASTVEP_SCORE_FIELDS = [
    ("revel_score",         ["revel", "score"]),
    ("phylop",              ["phylop"]),
    ("gerp",                ["gerp"]),
    ("spliceai_dsAg",       ["spliceAI", "dsAg"]),
    ("spliceai_dsAl",       ["spliceAI", "dsAl"]),
    ("spliceai_dsDg",       ["spliceAI", "dsDg"]),
    ("spliceai_dsDl",       ["spliceAI", "dsDl"]),
    ("spliceai_max_ds",     ["spliceAI", "SPLICEAI_MAX"]),
    ("spliceai_gene",       ["spliceAI", "gene"]),
    ("gnomad_allAf",        ["gnomad", "allAf"]),
    ("gnomad_allAc",        ["gnomad", "allAc"]),
    ("gnomad_allAn",        ["gnomad", "allAn"]),
    ("gnomad_allHc",        ["gnomad", "allHc"]),
    ("gnomad_afrAf",        ["gnomad", "afrAf"]),
    ("gnomad_amrAf",        ["gnomad", "amrAf"]),
    ("gnomad_asjAf",        ["gnomad", "asjAf"]),
    ("gnomad_easAf",        ["gnomad", "easAf"]),
    ("gnomad_finAf",        ["gnomad", "finAf"]),
    ("gnomad_midAf",        ["gnomad", "midAf"]),
    ("gnomad_nfeAf",        ["gnomad", "nfeAf"]),
    ("gnomad_remainingAf",  ["gnomad", "remainingAf"]),
    ("gnomad_sasAf",        ["gnomad", "sasAf"]),
    ("gnomad_max_pop_af",   ["gnomad", "MAX_POP_AF"]),
    ("acmg_classification", ["acmg", "classification"]),
    ("acmg_triggered_rule", ["acmg", "triggered_rule"]),
]
# Gene-level fields from `record.genes[gene]`. Same accessor convention.
FASTVEP_GENE_FIELDS = [
    ("gene_pLI",                    ["gnomad_genes", "pLI"]),
    ("gene_LOEUF",                  ["gnomad_genes", "loeuf"]),
    ("gene_misZ",                   ["gnomad_genes", "misZ"]),
    ("gene_synZ",                   ["gnomad_genes", "synZ"]),
    ("gene_omim_phenotypes",        ["omim", "phenotypes", "JOIN"]),
    ("gene_clinvar_protein_n",      ["clinvar_protein", "PROTEIN_VARIANT_COUNT"]),
]

OPPOSITE_PAIRS = {
    ("Pathogenic", "Likely_benign"), ("Pathogenic", "Benign"),
    ("Likely_pathogenic", "Likely_benign"), ("Likely_pathogenic", "Benign"),
    ("Likely_benign", "Pathogenic"), ("Likely_benign", "Likely_pathogenic"),
    ("Benign", "Pathogenic"), ("Benign", "Likely_pathogenic"),
}
EXTREME_PAIRS = {
    ("Pathogenic", "Benign"), ("Benign", "Pathogenic"),
    ("Likely_benign", "Pathogenic"), ("Pathogenic", "Likely_benign"),
}


def vep_allele(ref: str, alt: str) -> str:
    """VCF (REF, ALT) → VEP CSQ Allele convention."""
    if ref == alt:
        return alt
    i = 0
    while i < len(ref) and i < len(alt) and ref[i] == alt[i]:
        i += 1
    new_alt = alt[i:]
    return new_alt if new_alt else "-"


def vep_ref_alt(ref: str, alt: str) -> tuple[str, str]:
    """Convert a VCF (REF, ALT) pair to the (ref, alt) form VEP emits in
    the JSON `allele_string` field. Strips the leading common prefix; an
    empty side becomes the literal `-`."""
    if ref == alt:
        return ref, alt
    i = 0
    while i < len(ref) and i < len(alt) and ref[i] == alt[i]:
        i += 1
    nr = ref[i:] or "-"
    na = alt[i:] or "-"
    return nr, na


def vep_pos_ref_alt(pos: str, ref: str, alt: str) -> tuple[str, str, str]:
    """Convert a VCF (pos, REF, ALT) trio to the (pos, ref, alt) form
    VEP emits in the JSON. VEP not only strips the leading common
    prefix but also advances the start position by the strip count —
    for an insertion `chr2:47806652 G>GTAAC` the JSON record's
    `start` is `47806653` and `allele_string` is `-/TAAC`."""
    if ref == alt:
        return pos, ref, alt
    i = 0
    while i < len(ref) and i < len(alt) and ref[i] == alt[i]:
        i += 1
    nr = ref[i:] or "-"
    na = alt[i:] or "-"
    return str(int(pos) + i), nr, na


def parse_info_field(info: str, key: str) -> str:
    """Pull `key=...` value from a VCF INFO column. Empty if absent."""
    prefix = key + "="
    for piece in info.split(";"):
        if piece.startswith(prefix):
            return piece[len(prefix):]
    return ""


def load_clinvar_index() -> dict[tuple[str, str, str, str], dict[str, str]]:
    """Index ClinVar INFO fields by (chrom, pos, ref, alt)."""
    out: dict[tuple[str, str, str, str], dict[str, str]] = {}
    with open(CLINVAR_VCF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom = cols[0].removeprefix("chr")
            pos = cols[1]
            ref = cols[3]
            alts = cols[4].split(",")
            info = cols[7]
            ann = {k: parse_info_field(info, k) for k in CLINVAR_FIELDS}
            for alt in alts:
                out[(chrom, pos, ref, alt)] = ann
    return out


def _open_vcf(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def load_fastvep_csq(path: Path = FASTVEP_VCF) -> tuple[
    dict[tuple[str, str, str, str], dict[str, str]],
    dict[str, int],
]:
    """Stream a fastVEP VCF (gzipped or plain) and pick the canonical /
    first-with-ACMG CSQ entry per (chrom, pos, ref, alt)."""
    csq_idx: dict[str, int] | None = None
    out: dict[tuple[str, str, str, str], dict[str, str]] = {}
    with _open_vcf(path) as f:
        for line in f:
            if line.startswith("##"):
                if "ID=CSQ" in line:
                    fmt = line.split("Format: ", 1)[1].rstrip('">\n')
                    csq_idx = {n: i for i, n in enumerate(fmt.split("|"))}
                continue
            if line.startswith("#"):
                continue
            assert csq_idx is not None
            cols = line.rstrip("\n").split("\t")
            chrom = cols[0].removeprefix("chr")
            pos = cols[1]
            ref = cols[3]
            alts = cols[4].split(",")
            info = cols[7]
            csq_str = parse_info_field(info, "CSQ")
            if not csq_str:
                continue
            entries = [e.split("|") for e in csq_str.split(",")]
            ACMG = csq_idx["ACMG"]
            CAN = csq_idx.get("CANONICAL", -1)
            by_alt: dict[str, list[list[str]]] = defaultdict(list)
            for parts in entries:
                if len(parts) <= ACMG:
                    continue
                by_alt[parts[0]].append(parts)
            for alt in alts:
                csq_alt = vep_allele(ref, alt)
                pool = by_alt.get(csq_alt, [])
                if not pool:
                    continue
                populated = [p for p in pool if p[ACMG]] or pool
                canon = [
                    p for p in populated
                    if CAN >= 0 and len(p) > CAN and p[CAN] == "YES"
                ]
                chosen = canon[0] if canon else populated[0]
                out[(chrom, pos, ref, alt)] = {
                    name: chosen[csq_idx[name]] if csq_idx[name] < len(chosen) else ""
                    for name in FASTVEP_FIELDS
                }
    return out, csq_idx or {}


def _accessor(obj: object, path: list[str]):
    """Navigate `path` through nested dicts; handle computed transforms."""
    cur: object = obj
    for i, key in enumerate(path):
        if cur is None:
            return ""
        if key == "SPLICEAI_MAX":
            # SpliceAI max delta score across the four prediction heads.
            if not isinstance(cur, dict):
                return ""
            ds = [cur.get(k) for k in ("dsAg", "dsAl", "dsDg", "dsDl")]
            ds = [v for v in ds if isinstance(v, (int, float))]
            return max(ds) if ds else ""
        if key == "MAX_POP_AF":
            # Max across all populations the parser knows about. Mirrors
            # `GnomadData::max_pop_af()` in fastvep_classification.
            if not isinstance(cur, dict):
                return ""
            keys = [
                "allAf", "afrAf", "amrAf", "asjAf", "easAf", "finAf",
                "midAf", "nfeAf", "othAf", "remainingAf", "sasAf",
            ]
            af = [cur.get(k) for k in keys]
            af = [v for v in af if isinstance(v, (int, float))]
            return max(af) if af else ""
        if key == "JOIN":
            # Stringify a list of phenotype names with "; ".
            if isinstance(cur, list):
                return "; ".join(str(x) for x in cur)
            return str(cur)
        if key == "PROTEIN_VARIANT_COUNT":
            # ClinVar-protein index: count of pathogenic protein variants
            # registered for this gene (PS1 / PM5 / PM1 driver source).
            if isinstance(cur, dict):
                pv = cur.get("proteinVariants") or []
                return len(pv)
            return 0
        if not isinstance(cur, dict):
            return ""
        cur = cur.get(key)
    return cur if cur is not None else ""


def _format_value(v: object) -> str:
    if v is None or v == "":
        return ""
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def load_fastvep_json(path: Path) -> dict[
    tuple[str, str, str, str], dict[str, str]
]:
    """Pull the picked-tc score panel for each variant out of the JSON.

    JSON keys variants by VEP-form (chrom, pos, vep_ref, vep_alt) so we
    return the index in the same form. Callers convert their VCF-form
    (REF, ALT) via `vep_ref_alt(...)` before lookup.
    """
    import json as _json
    out: dict[tuple[str, str, str, str], dict[str, str]] = {}
    with _open_vcf(path) as f:
        recs = _json.load(f)
    for r in recs:
        chrom = str(r.get("seq_region_name", ""))
        pos = str(r.get("start", ""))
        allele = r.get("allele_string") or ""
        if "/" in allele:
            ref, alts = allele.split("/", 1)
        else:
            # No "/" — happens for "*" / "." alt placeholders. Skip.
            continue
        # Pick the same transcript the VCF picked: prefer canonical+ACMG.
        tcs = r.get("transcript_consequences") or []
        populated = [tc for tc in tcs if "acmg" in tc and tc["acmg"].get("classification")]
        pool = populated or [tc for tc in tcs if "acmg" in tc] or tcs
        canon = [tc for tc in pool if tc.get("canonical")]
        chosen = canon[0] if canon else (pool[0] if pool else {})
        gene_sym = chosen.get("gene_symbol") or ""
        gene_obj = (r.get("genes") or {}).get(gene_sym) or {}

        row: dict[str, str] = {}
        for name, path_ in FASTVEP_SCORE_FIELDS:
            row[name] = _format_value(_accessor(chosen, path_))
        for name, path_ in FASTVEP_GENE_FIELDS:
            row[name] = _format_value(_accessor(gene_obj, path_))

        for alt in alts.split(","):
            out[(chrom, pos, ref, alt)] = row
    return out


def priority_score(stars: str, truth: str, predicted: str, n_crit: int) -> int:
    star_w = 100 if stars == "3" else 10
    extremity = 500 if (truth, predicted) in EXTREME_PAIRS else 100
    return star_w * extremity + n_crit


def main() -> None:
    print("Loading truth + ClinVar INFO + fastVEP CSQ ...")
    with open(TRUTH) as f:
        truth = {
            (r["chrom"], r["pos"], r["ref"], r["alt"]): r
            for r in csv.DictReader(f, delimiter="\t")
        }
    clinvar_idx = load_clinvar_index()
    fastvep_idx, _ = load_fastvep_csq()
    if FASTVEP_HGVS_VCF.exists():
        # Overlay the HGVS-populated re-annotation for the discrepancy
        # subset. Other 673,548 variants keep their full-benchmark CSQ
        # entry (HGVS fields empty there but those rows aren't in the
        # output anyway).
        hgvs_idx, _ = load_fastvep_csq(FASTVEP_HGVS_VCF)
        fastvep_idx.update(hgvs_idx)
        print(f"  overlaid {len(hgvs_idx):,} HGVS-annotated variants")
    json_idx: dict[tuple[str, str, str, str], dict[str, str]] = {}
    if FASTVEP_JSON.exists():
        json_idx = load_fastvep_json(FASTVEP_JSON)
        print(f"  loaded JSON score panel for {len(json_idx):,} variants")
    print(f"  truth: {len(truth):,}  clinvar: {len(clinvar_idx):,}  "
          f"fastvep: {len(fastvep_idx):,}")

    print("Filtering opposite-direction discrepancies...")
    rows: list[dict[str, object]] = []
    with DISC.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if (r["truth"], r["predicted"]) not in OPPOSITE_PAIRS:
                continue
            rows.append(r)
    print(f"  {len(rows)} opposite-direction discrepancies")

    score_field_names = [name for name, _ in FASTVEP_SCORE_FIELDS]
    gene_field_names = [name for name, _ in FASTVEP_GENE_FIELDS]
    header = (
        ["priority_score"]
        + ["chrom", "pos", "ref", "alt", "gene", "stars"]
        + ["truth_class", "fastvep_class", "n_criteria_met", "fastvep_met_criteria"]
        + ["consequence_top"]
        + [f"clinvar_{k}" for k in CLINVAR_FIELDS]
        + [f"fastvep_{k}" for k in FASTVEP_FIELDS]
        + [f"fastvep_{k}" for k in score_field_names]
        + [f"fastvep_{k}" for k in gene_field_names]
        + ["review_question"]
    )
    out_rows = []
    for r in rows:
        key = (r["chrom"], r["pos"], r["ref"], r["alt"])
        n_crit = len([c for c in r["met_criteria"].split(";") if c])
        score = priority_score(r["stars"], r["truth"], r["predicted"], n_crit)
        cv = clinvar_idx.get(key, {})
        fv = fastvep_idx.get(key, {})
        # JSON index is keyed by VEP-form (pos, ref, alt). For indels
        # VEP shifts pos by the prefix-strip count and uses `-` for
        # the collapsed side.
        v_pos, v_ref, v_alt = vep_pos_ref_alt(r["pos"], r["ref"], r["alt"])
        scores = json_idx.get((r["chrom"], v_pos, v_ref, v_alt), {})
        review_q = (
            f"Why does fastVEP call {r['predicted']} when ClinVar "
            f"({r['stars']}-star) says {r['truth']}? "
            f"Inspect: {fv.get('HGVSp') or fv.get('HGVSc') or 'no HGVS'}; "
            f"criteria fired = {r['met_criteria'] or '(none)'}"
        )
        out_rows.append({
            "priority_score": score,
            "chrom": r["chrom"], "pos": r["pos"],
            "ref": r["ref"], "alt": r["alt"],
            "gene": r["gene"], "stars": r["stars"],
            "truth_class": r["truth"],
            "fastvep_class": r["predicted"],
            "n_criteria_met": n_crit,
            "fastvep_met_criteria": r["met_criteria"],
            "consequence_top": r["consequence"],
            **{f"clinvar_{k}": cv.get(k, "") for k in CLINVAR_FIELDS},
            **{f"fastvep_{k}": fv.get(k, "") for k in FASTVEP_FIELDS},
            **{f"fastvep_{name}": scores.get(name, "") for name in score_field_names},
            **{f"fastvep_{name}": scores.get(name, "") for name in gene_field_names},
            "review_question": review_q,
        })
    out_rows.sort(key=lambda d: -d["priority_score"])

    with OUT_TSV.open("w") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        w.writeheader()
        for row in out_rows:
            w.writerow(row)
    print(f"\nWrote {OUT_TSV} ({len(out_rows)} rows, {len(header)} columns)")
    print()
    print("Top 5 review cases:")
    for r in out_rows[:5]:
        hgvs = r["fastvep_HGVSp"] or r["fastvep_HGVSc"]
        cv_dn = (r["clinvar_CLNDN"] or "").split("|")[0].replace("_", " ")[:60]
        print(f"  {r['gene']:<10} {r['stars']}*  "
              f"{r['truth_class']:<18} → {r['fastvep_class']:<18}  "
              f"{hgvs}  ({cv_dn})")


if __name__ == "__main__":
    main()
