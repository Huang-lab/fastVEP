#!/usr/bin/env python3
"""Every number and table the ACMG manuscript reports about the final classifier.

Reads the per-variant table of the final run (`12_per_variant_table.py`), the
concordance matrix of every other run on disk, and the sweeps and ablations of
`11_manuscript_runs.py`, and writes small CSV/TSV files to `manuscript/acmg/data/`.
Nothing in the manuscript's tables or figures is typed in by hand: they are read
from those files, and each file is written by one function below.

ClinVar's 2-star+ classifications are a reference and not ground truth. Every
rate here is agreement with them, and says nothing stronger.

Definitions (the same as `03_evaluate_concordance.py`, which scores every run)
  exact           the call is ClinVar's class, on the five-tier scale
  same-direction  both pathogenic-side (P or LP), both benign-side (LB or B), or both VUS
  opposite        ClinVar P/LP and fastVEP LB/B, or the reverse
  false-benign    ClinVar P/LP called LB/B (a missed diagnosis)
  false-pathogenic ClinVar LB/B called P/LP
ClinVar's `Pathogenic/Likely pathogenic` is scored as `Pathogenic` and its
`Benign/Likely benign` as `Benign`, so `exact` is a strict measure and
`same-direction` the fair one.

Usage
  13_manuscript_stats.py            write every table
  13_manuscript_stats.py --only headline,stars
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BENCH = ROOT / "data/benchmark"
OUT = ROOT / "manuscript/acmg/data"
RES = ROOT / "analysis/acmg_benchmark/results"
PER_VARIANT = BENCH / "manuscript_runs/v37_per_variant.tsv.gz"
RUNS = RES / "v37_manuscript/runs.tsv"

CLASSES = ["Pathogenic", "Likely_pathogenic", "VUS", "Likely_benign", "Benign"]
CALLS = CLASSES + ["NoCall"]
SHORT = {"Pathogenic": "P", "Likely_pathogenic": "LP", "VUS": "VUS", "Likely_benign": "LB", "Benign": "B", "NoCall": "NoCall"}
TRI = {"Pathogenic": "P/LP", "Likely_pathogenic": "P/LP", "VUS": "VUS", "Likely_benign": "B/LB", "Benign": "B/LB", "NoCall": "NoCall"}
PATH = ("Pathogenic", "Likely_pathogenic")
BEN = ("Likely_benign", "Benign")


# --------------------------------------------------------------------------- helpers


def write(name: str, rows: list[dict] | pd.DataFrame) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
    df.to_csv(OUT / name, index=False, lineterminator="\n")
    print(f"  wrote {name} ({len(df)} rows)")


def load() -> pd.DataFrame:
    df = pd.read_csv(PER_VARIANT, sep="\t", dtype=str, keep_default_na=False)
    df["cv3"] = df["clinvar"].map(TRI)
    df["fv3"] = df["call"].map(TRI)
    df["stars"] = df["stars"].astype(int)
    return df


def read_matrix(path: Path) -> dict[str, dict[str, int]]:
    rows = list(csv.reader(path.open()))
    hdr = rows[0][1:]
    m = {r[0]: {c: int(v) for c, v in zip(hdr, r[1:])} for r in rows[1:]}
    for t in CLASSES:
        m[t].setdefault("NoCall", 0)
    return m


def matrix_metrics(m: dict[str, dict[str, int]]) -> dict[str, float]:
    n = sum(sum(m[t].values()) for t in CLASSES)
    nocall = sum(m[t]["NoCall"] for t in CLASSES)
    exact = sum(m[t][t] for t in CLASSES)
    same_p = sum(m[t][p] for t in PATH for p in PATH)
    same_b = sum(m[t][p] for t in BEN for p in BEN)
    same = same_p + same_b + m["VUS"]["VUS"]
    fb = sum(m[t][p] for t in PATH for p in BEN)
    fp = sum(m[t][p] for t in BEN for p in PATH)
    n_p = sum(sum(m[t].values()) for t in PATH)
    n_b = sum(sum(m[t].values()) for t in BEN)
    n_v = sum(m["VUS"].values())
    p_to_vus = sum(m[t]["VUS"] for t in PATH)
    b_to_vus = sum(m[t]["VUS"] for t in BEN)
    v_to_p = sum(m["VUS"][p] for p in PATH)
    v_to_b = sum(m["VUS"][p] for p in BEN)
    return {
        "n": n, "no_call": nocall, "exact": exact, "same_direction": same,
        "same_pathogenic": same_p, "same_benign": same_b, "vus_exact": m["VUS"]["VUS"],
        "false_benign": fb, "false_pathogenic": fp, "opposite": fb + fp,
        "n_clinvar_plp": n_p, "n_clinvar_blb": n_b, "n_clinvar_vus": n_v,
        "plp_to_vus": p_to_vus, "blb_to_vus": b_to_vus, "vus_to_plp": v_to_p, "vus_to_blb": v_to_b,
        "plp_recall": same_p / n_p, "blb_recall": same_b / n_b, "vus_recall": m["VUS"]["VUS"] / n_v,
        "exact_rate": exact / n, "same_direction_rate": same / n,
    }


def clopper_pearson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    if n == 0:
        return (math.nan, math.nan)
    lo = 0.0 if k == 0 else stats.beta.ppf(alpha / 2, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(1 - alpha / 2, k + 1, n - k)
    return float(lo), float(hi)


AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
       "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}


def read_protein_index(path: Path) -> dict[str, list[dict]]:
    """The classified missense variants of each gene in the ClinVar protein index the classifier reads.

    An `.oga` file is the magic `FSTGA_01`, a 2-byte schema, an 8-byte length and a bincode `GeneIndex`: a header of
    five fields (a 2-byte schema and four strings) and a map from gene symbol to JSON strings, with every integer
    little-endian and every string and sequence preceded by its 8-byte length (`crates/fastvep-sa/src/gene.rs`).
    """
    import json
    import struct

    raw = path.read_bytes()
    if raw[:8] != b"FSTGA_01":
        sys.exit(f"{path} is not an .oga file")
    body, o = raw[18:18 + struct.unpack("<Q", raw[10:18])[0]], 0

    def u64() -> int:
        nonlocal o
        o += 8
        return struct.unpack("<Q", body[o - 8:o])[0]

    def text() -> str:
        nonlocal o
        n = u64()
        o += n
        return body[o - n:o].decode()

    o += 2
    for _ in range(4):
        text()
    out = {}
    for _ in range(u64()):
        gene = text()
        records = [text() for _ in range(u64())]
        out[gene] = json.loads(records[0])["proteinVariants"]
    if o != len(body):
        sys.exit(f"{path}: {len(body) - o} bytes left over after the index")
    return out


def _index_row(index: dict[str, list[dict]], gene: str, hgvsp: str, kind: str) -> dict:
    """The pathogenic missense records of `gene` in the protein index, and whether the variant's own residue change is among them."""
    m = re.search(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})", hgvsp)
    own = (int(m.group(2)), AA3[m.group(1)], AA3[m.group(3)])
    records = [v for v in index[gene] if "pathogenic" in v["sig"].lower()]
    has_own = any((v["pos"], v["refAa"], v["altAa"]) == own for v in records)
    aa1 = {v: k for k, v in AA3.items()}
    others = [f"p.{aa1.get(v['refAa'], v['refAa'])}{v['pos']}{aa1.get(v['altAa'], v['altAa'])}" for v in records if (v["pos"], v["refAa"], v["altAa"]) != own]
    return {"gene": gene, "variant": hgvsp.split(":")[-1], "kind": kind, "pathogenic_missense_records": len(records),
            "own_record_among_them": has_own, "records_other_than_own": len(records) - int(has_own), "the_others": ";".join(others)}


# ------------------------------------------------------------- consequence groups


def consequence_group(c: str) -> str:
    if c == "missense_variant":
        return "Missense"
    if c in ("synonymous_variant", "stop_retained_variant"):
        return "Synonymous"
    if c in ("frameshift_variant", "stop_gained"):
        return "Frameshift or nonsense"
    if c in ("splice_donor_variant", "splice_acceptor_variant"):
        return "Canonical splice site"
    if c in ("splice_region_variant", "splice_polypyrimidine_tract_variant",
             "splice_donor_region_variant", "splice_donor_5th_base_variant"):
        return "Splice region"
    if c in ("inframe_deletion", "inframe_insertion", "protein_altering_variant", "coding_sequence_variant"):
        return "In-frame or protein-altering"
    if c in ("start_lost", "stop_lost"):
        return "Start or stop loss"
    if c == "intron_variant":
        return "Intronic"
    if c in ("3_prime_UTR_variant", "5_prime_UTR_variant"):
        return "UTR"
    return "Other (flank, non-coding, unscored)"


# ------------------------------------------------------------ evidence points


def code_points(code: str) -> tuple[str, int]:
    """Direction and ClinGen SVI points of one met criterion, from the code the output carries.

    The output names a strength only where it differs from the criterion's default
    (`PP3_Moderate`, `PM2_Supporting`, `PVS1_Strong`), so a bare code takes the
    default of Richards 2015: PVS1 very strong (8), PS strong (4), PM moderate (2),
    PP supporting (1), BA1 stand-alone (8), BS strong (4), BP supporting (1).
    """
    base = code.split("_")[0]
    if code == "BA1":
        return "B", 8
    suffix = code[len(base) + 1:] if "_" in code else ""
    direction = "P" if base.startswith("P") else "B"
    if suffix == "Supporting":
        pts = 1
    elif suffix == "Moderate":
        pts = 2
    elif suffix == "Strong":
        pts = 4
    elif suffix == "Very_Strong":
        pts = 8
    elif base == "PVS1":
        pts = 8
    elif base.startswith(("PS", "BS")):
        pts = 4
    elif base.startswith("PM"):
        pts = 2
    else:
        pts = 1
    return direction, pts


def strength_counts(codes: list[str]):
    """The strength counts `EvidenceCounts::from_criteria` makes, and the pathogenic points."""
    pvs = ps = pm = pp = ba = bs = bp = 0
    pts = 0
    for c in codes:
        d, p = code_points(c)
        if d == "P":
            pts += p
            pvs += p == 8
            ps += p == 4
            pm += p == 2
            pp += p == 1
        elif c == "BA1":
            ba += 1
        elif p == 4 or p == 2:
            bs += 1  # Benign Moderate is counted as Strong
        elif p == 8:
            bs += 2  # Benign Very Strong is two Strong
        else:
            bp += 1
    return pvs, ps, pm, pp, ba, bs, bp, pts


def table_pathogenic(pvs, ps, pm, pp):
    if pvs >= 1 and ps >= 1:
        return "P"
    if pvs >= 1 and pm >= 2:
        return "P"
    if pvs >= 1 and pm >= 1 and pp >= 1:
        return "P"
    if pvs >= 1 and pp >= 2:
        return "P"
    if ps >= 2:
        return "P"
    if ps >= 1 and pm >= 3:
        return "P"
    if ps >= 1 and pm >= 2 and pp >= 2:
        return "P"
    if ps >= 1 and pm >= 1 and pp >= 4:
        return "P"
    if pvs >= 1 and pm >= 1:
        return "LP"
    if pvs >= 1 and pp >= 1:
        return "LP"
    if ps >= 1 and 1 <= pm <= 2:
        return "LP"
    if ps >= 1 and pp >= 2:
        return "LP"
    if pm >= 3:
        return "LP"
    if pm >= 2 and pp >= 2:
        return "LP"
    if pm >= 1 and pp >= 4:
        return "LP"
    return None


def table_benign(ba, bs, bp):
    if ba >= 1 or bs >= 2:
        return "B"
    if (bs >= 1 and bp >= 1) or bp >= 2:
        return "LB"
    return None


def replay(codes: list[str], *, scheme: str = "shipped", lp_cut: int = 6, p_cut: int = 10,
           lb_floor: int | None = None) -> str:
    """Re-derive a call from the criteria a variant met, under a combining scheme.

    `shipped` is v37's own: SVI points on the pathogenic side, the Richards 2015 table
    on the benign side, and the rule that Strong or stronger evidence on one side
    blocks a definite call on the other. `table` is the 2015 table on both sides.
    `points` scores both sides on the SVI scale; `lb_floor` is the benign-points
    edge of Likely Benign (-1 is Tavtigian's, -2 the stricter one tried in v10).

    This is a replay of `combiner.rs` from the codes in the output. It cannot see the
    one thing the output omits, the site-level frequency veto, which takes a lone
    Likely Pathogenic to VUS where gnomAD rejected the site; the manuscript reports
    how many calls that is.
    """
    pvs, ps, pm, pp, ba, bs, bp, pts = strength_counts(codes)
    benign_pts = 0
    for c in codes:
        d, p = code_points(c)
        if d == "B":
            benign_pts += p
    if scheme == "table":
        pc = table_pathogenic(pvs, ps, pm, pp)
    else:
        pc = "P" if pts >= p_cut else "LP" if pts >= lp_cut else None
    if scheme == "points":
        net = pts - benign_pts
        floor = 1 if lb_floor is None else lb_floor
        if net >= p_cut:
            return "P"
        if net >= lp_cut:
            return "LP"
        if net <= -7 or ba >= 1:
            return "B"
        if net <= -floor:
            return "LB"
        return "VUS"
    bc = table_benign(ba, bs, bp)
    strong_p = pvs >= 1 or ps >= 1
    strong_b = ba >= 1 or bs >= 1
    if bc and strong_p and pc is None:
        return "VUS"
    if pc and strong_b and bc is None:
        return "VUS"
    if pc and bc:
        return "VUS"
    return pc or bc or "VUS"


# ---------------------------------------------------------------- the analyses


def headline(df: pd.DataFrame) -> None:
    m = read_matrix(BENCH / "output_v37/concordance_matrix.csv")
    mm = matrix_metrics(m)
    scored = mm["n"] - mm["no_call"]
    if mm["n"] != len(df):
        sys.exit(f"the per-variant table has {len(df)} rows and the matrix {mm['n']}")

    # Every record is in the denominator, as in the scorer: a reference-agreement record (ALT is
    # '.') produces no call and counts as not concordant. Gene-cluster bootstrap: variants in one
    # gene are not independent, so the interval comes from resampling genes, not variants.
    sub = df.copy()
    sub["exact_f"] = (sub["clinvar"] == sub["call"]).astype(int)
    sub["same_f"] = (sub["cv3"] == sub["fv3"]).astype(int)
    sub["plp"] = (sub["cv3"] == "P/LP").astype(int)
    sub["blb"] = (sub["cv3"] == "B/LB").astype(int)
    sub["vus"] = (sub["cv3"] == "VUS").astype(int)
    sub["plp_same"] = ((sub["cv3"] == "P/LP") & (sub["fv3"] == "P/LP")).astype(int)
    sub["blb_same"] = ((sub["cv3"] == "B/LB") & (sub["fv3"] == "B/LB")).astype(int)
    sub["vus_same"] = ((sub["cv3"] == "VUS") & (sub["fv3"] == "VUS")).astype(int)
    sub["opp"] = (((sub["cv3"] == "P/LP") & (sub["fv3"] == "B/LB")) | ((sub["cv3"] == "B/LB") & (sub["fv3"] == "P/LP"))).astype(int)
    sub["called_plp"] = (sub["fv3"] == "P/LP").astype(int)
    sub["called_blb"] = (sub["fv3"] == "B/LB").astype(int)
    sub["plp_when_called_plp"] = ((sub["fv3"] == "P/LP") & (sub["cv3"] == "P/LP")).astype(int)
    sub["blb_when_called_blb"] = ((sub["fv3"] == "B/LB") & (sub["cv3"] == "B/LB")).astype(int)
    sub["one"] = 1
    g = sub.groupby("gene")[["one", "exact_f", "same_f", "plp", "blb", "vus", "plp_same", "blb_same", "vus_same",
                             "opp", "called_plp", "called_blb", "plp_when_called_plp", "blb_when_called_blb"]].sum()
    arr = g.to_numpy(dtype=float)
    cols = {c: i for i, c in enumerate(g.columns)}
    rng = np.random.default_rng(20261011)
    B = 2000
    idx = rng.integers(0, len(arr), size=(B, len(arr)))
    boot = np.stack([arr[i].sum(axis=0) for i in idx])

    def ratio(num, den):
        return boot[:, cols[num]] / boot[:, cols[den]]

    def ci(v):
        return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))

    tot = arr.sum(axis=0)

    def point(num, den):
        return tot[cols[num]] / tot[cols[den]]

    spec = [
        ("exact_match", "exact_f", "one", "Exact match (five-tier), of all variants"),
        ("same_direction", "same_f", "one", "Same-direction (P/LP, VUS, B/LB), of all variants"),
        ("plp_recall", "plp_same", "plp", "ClinVar P/LP called P or LP"),
        ("blb_recall", "blb_same", "blb", "ClinVar B/LB called LB or B"),
        ("vus_recall", "vus_same", "vus", "ClinVar VUS called VUS"),
        ("opposite_rate", "opp", "one", "Opposite-direction, of all variants"),
        ("precision_plp_calls", "plp_when_called_plp", "called_plp", "fastVEP P/LP calls that ClinVar also calls P/LP"),
        ("precision_blb_calls", "blb_when_called_blb", "called_blb", "fastVEP B/LB calls that ClinVar also calls B/LB"),
    ]
    rows = []
    for key, num, den, label in spec:
        v = point(num, den)
        lo, hi = ci(ratio(num, den))
        rows.append({"metric": key, "label": label, "value": v, "ci_lo": lo, "ci_hi": hi,
                     "numerator": int(tot[cols[num]]), "denominator": int(tot[cols[den]]),
                     "ci": "gene-cluster bootstrap, 2000 resamples, 95% percentile"})
    k, n = mm["opposite"], mm["n"]
    lo, hi = clopper_pearson(k, n)
    rows.append({"metric": "opposite_rate_exact_binomial", "label": "Opposite-direction, of all variants (exact binomial, treating variants as independent)",
                 "value": k / n, "ci_lo": lo, "ci_hi": hi, "numerator": k, "denominator": n,
                 "ci": "Clopper-Pearson 95%"})
    n_genes = len(g)
    if int(tot[cols["exact_f"]]) != mm["exact"] or int(tot[cols["same_f"]]) != mm["same_direction"]:
        sys.exit("the per-variant table and the scorer's matrix disagree on exact or same-direction counts")
    for key, val, label in [
        ("variants_total", mm["n"], "ClinVar 2-star+ variants (SNVs and small indels, GRCh38)"),
        ("variants_scored", scored, "Scored (an ACMG call was made)"),
        ("variants_no_call", mm["no_call"], "No call (ClinVar reference-agreement records, ALT is '.')"),
        ("genes", n_genes, "Distinct genes among the variants"),
        ("exact", mm["exact"], "Exact matches"),
        ("same_direction", mm["same_direction"], "Same-direction calls"),
        ("opposite", mm["opposite"], "Opposite-direction calls"),
        ("false_benign", mm["false_benign"], "ClinVar P/LP called LB/B"),
        ("false_pathogenic", mm["false_pathogenic"], "ClinVar LB/B called P/LP"),
    ]:
        rows.append({"metric": key, "label": label, "value": val, "ci_lo": "", "ci_hi": "",
                     "numerator": "", "denominator": "", "ci": ""})
    write("headline_v37.csv", rows)


def confusion(df: pd.DataFrame) -> None:
    m = read_matrix(BENCH / "output_v37/concordance_matrix.csv")
    rows = []
    for t in CLASSES:
        for p in CALLS:
            rows.append({"clinvar": t, "fastvep": p, "n": m[t][p]})
    write("confusion_5x5_v37.csv", rows)

    ct = pd.crosstab(df["cv3"], df["fv3"])
    rows = []
    for a in ("P/LP", "VUS", "B/LB"):
        tot = int(ct.loc[a].sum())
        for b in ("P/LP", "VUS", "B/LB", "NoCall"):
            k = int(ct.loc[a, b]) if b in ct.columns else 0
            lo, hi = clopper_pearson(k, tot)
            rows.append({"clinvar": a, "fastvep": b, "n": k, "row_total": tot, "row_fraction": k / tot,
                         "ci_lo": lo, "ci_hi": hi})
    write("confusion_3x3_v37.csv", rows)

    # Column view: what is behind each of fastVEP's own calls.
    rows = []
    for b in ("P/LP", "VUS", "B/LB"):
        tot = int(ct[b].sum())
        for a in ("P/LP", "VUS", "B/LB"):
            k = int(ct.loc[a, b])
            lo, hi = clopper_pearson(k, tot)
            rows.append({"fastvep": b, "clinvar": a, "n": k, "column_total": tot, "column_fraction": k / tot,
                         "ci_lo": lo, "ci_hi": hi})
    write("what_is_behind_each_call_v37.csv", rows)


def by_stars(df: pd.DataFrame) -> None:
    rows = []
    for s in (2, 3, 4):
        sub = df[df.stars == s]
        for a in ("P/LP", "VUS", "B/LB"):
            sa = sub[sub.cv3 == a]
            tot = len(sa)
            for b in ("P/LP", "VUS", "B/LB", "NoCall"):
                k = int((sa.fv3 == b).sum())
                lo, hi = clopper_pearson(k, tot) if tot else (math.nan, math.nan)
                rows.append({"stars": s, "clinvar": a, "fastvep": b, "n": k, "row_total": tot,
                             "row_fraction": k / tot if tot else math.nan, "ci_lo": lo, "ci_hi": hi})
    write("by_review_status_v37.csv", rows)

    # The 4-star records, individually: there are 23 of them.
    s4 = df[df.stars == 4][["gene", "chrom", "pos", "ref", "alt", "clinvar", "call", "criteria", "consequence", "hgvsc", "hgvsp"]]
    write("four_star_variants_v37.csv", s4)


def by_consequence(df: pd.DataFrame) -> None:
    df = df.copy()
    df["cgrp"] = df["consequence"].map(consequence_group)
    rows = []
    for grp, sub in df.groupby("cgrp"):
        for a in ("P/LP", "VUS", "B/LB"):
            sa = sub[sub.cv3 == a]
            tot = len(sa)
            for b in ("P/LP", "VUS", "B/LB", "NoCall"):
                k = int((sa.fv3 == b).sum())
                rows.append({"consequence_group": grp, "clinvar": a, "fastvep": b, "n": k, "row_total": tot,
                             "row_fraction": k / tot if tot else math.nan})
    write("by_consequence_v37.csv", rows)


def signature(s: str) -> str:
    return "+".join(sorted(c for c in s.split("&") if c)) or "(none)"


def offdiagonal(df: pd.DataFrame) -> None:
    df = df.copy()
    df["sig"] = df["criteria"].map(signature)
    pb = df["criteria"].map(lambda s: strength_counts([c for c in s.split("&") if c])[7])
    df["path_points"] = pb
    rows = []
    for a, b in [("P/LP", "VUS"), ("B/LB", "VUS"), ("VUS", "P/LP"), ("VUS", "B/LB"), ("P/LP", "B/LB"), ("B/LB", "P/LP")]:
        sub = df[(df.cv3 == a) & (df.fv3 == b)]
        vc = sub["sig"].value_counts()
        for sig, k in vc.head(25).items():
            pts = int(sub[sub.sig == sig]["path_points"].iloc[0])
            rows.append({"clinvar": a, "fastvep": b, "signature": sig, "n": int(k),
                         "fraction_of_cell": k / len(sub), "pathogenic_points": pts, "cell_total": len(sub)})
    write("signatures_offdiagonal_v37.csv", rows)

    # Where in ClinVar P/LP -> VUS the evidence stops: pathogenic points reached.
    sub = df[(df.cv3 == "P/LP") & (df.fv3 == "VUS")]
    c = sub["path_points"].value_counts().sort_index()
    rows = [{"pathogenic_points": int(k), "n": int(v), "fraction": v / len(sub)} for k, v in c.items()]
    write("plp_to_vus_points_v37.csv", rows)

    # Which consequence classes the ClinVar P/LP -> VUS belong to, and whether PVS1 was in play.
    sub = sub.copy()
    sub["cgrp"] = sub["consequence"].map(consequence_group)
    sub["has_pvs1"] = sub["criteria"].str.contains("PVS1")
    rows = [{"consequence_group": g, "has_pvs1_family": bool(h), "n": int(len(x))}
            for (g, h), x in sub.groupby(["cgrp", "has_pvs1"])]
    write("plp_to_vus_by_consequence_v37.csv", rows)


def points_tables(df: pd.DataFrame) -> None:
    df = df.copy()
    codes = df["criteria"].map(lambda s: [c for c in s.split("&") if c])
    sc = codes.map(strength_counts)
    df["pp"] = [x[7] for x in sc]
    df["bp_pts"] = [sum(code_points(c)[1] for c in cs if code_points(c)[0] == "B") for cs in codes]
    sub = df[df.cv3 != "NoCall"]
    rows = []
    for ben_only in (False, True):
        s = sub[sub.bp_pts == 0] if ben_only else sub
        ct = pd.crosstab(s["pp"], s["cv3"])
        for pts in range(0, 13):
            if pts not in ct.index:
                continue
            r = ct.loc[pts]
            tot = int(r.sum())
            defined = int(r.get("P/LP", 0) + r.get("B/LB", 0))
            rows.append({
                "scope": "no benign evidence" if ben_only else "all variants",
                "pathogenic_points": pts, "n": tot,
                "clinvar_plp": int(r.get("P/LP", 0)), "clinvar_vus": int(r.get("VUS", 0)), "clinvar_blb": int(r.get("B/LB", 0)),
                "plp_fraction_of_all": r.get("P/LP", 0) / tot,
                "plp_fraction_of_defined": r.get("P/LP", 0) / defined if defined else math.nan,
            })
    write("points_by_clinvar_v37.csv", rows)

    # Tavtigian 2018 / 2020 expected probability of pathogenicity from points, for the same x axis,
    # at the 0.10 prior the framework is built on and at this benchmark's own prior among defined classes.
    n_plp = int((sub.cv3 == "P/LP").sum())
    n_blb = int((sub.cv3 == "B/LB").sum())
    prior_bench = n_plp / (n_plp + n_blb)
    rows = []
    for pts in range(-8, 14):
        for name, prior in (("framework prior 0.10", 0.10), ("benchmark prior", prior_bench)):
            odds = (350 ** (pts / 8)) * prior / (1 - prior)
            rows.append({"points": pts, "prior": name, "prior_value": prior, "posterior": odds / (1 + odds)})
    write("tavtigian_posterior.csv", rows)


def replay_tables(df: pd.DataFrame) -> None:
    codes = df["criteria"].map(lambda s: [c for c in s.split("&") if c])
    nocall = (df["call"] == "NoCall").to_numpy()
    actual = df["call"].map(SHORT).to_numpy()
    base = np.array([replay(c) if not nc else "NoCall" for c, nc in zip(codes, nocall)])
    differs = base != actual
    mism = int(differs.sum())
    print(f"  replay of the shipped scheme differs from the run on {mism} of {len(df)} variants")
    # The only difference there should be is the site-level frequency veto, which takes a lone
    # Likely Pathogenic to VUS and is not recoverable from the output. Check that it is, and
    # then hold those variants at VUS wherever a scheme would call them Likely Pathogenic.
    d = pd.crosstab(pd.Series(actual[differs], name="run"), pd.Series(base[differs], name="replay"))
    print(d.to_string())
    veto = differs & (actual == "VUS") & (base == "LP")
    if int(veto.sum()) != mism:
        sys.exit("the replay differs from the run in a way the site-level veto does not explain")
    print(f"  all {mism} are the site-level frequency veto")
    write("replay_veto_v37.csv", [{"variants": len(df), "replay_matches_the_run": len(df) - mism, "held_at_vus_by_the_veto": mism}])

    def call(c, nc, i, **kw):
        if nc:
            return "NoCall"
        r = replay(c, **kw)
        return "VUS" if (r == "LP" and veto[i]) else r

    long = {"P": "Pathogenic", "LP": "Likely_pathogenic", "VUS": "VUS", "LB": "Likely_benign", "B": "Benign", "NoCall": "NoCall"}
    truth = df["clinvar"].to_numpy()

    def matrix_of(pred):
        m = {t: {p: 0 for p in CALLS} for t in CLASSES}
        for t, p in zip(truth, pred):
            m[t][long[p]] += 1
        return m

    schemes = [
        ("shipped (v37)", {}),
        ("Richards 2015 table, both sides", {"scheme": "table"}),
        ("SVI points, both sides, Likely Benign from -1 (Tavtigian)", {"scheme": "points", "lb_floor": 1}),
        ("SVI points, both sides, Likely Benign from -2", {"scheme": "points", "lb_floor": 2}),
        ("SVI points, both sides, Likely Benign from -3", {"scheme": "points", "lb_floor": 3}),
    ]
    rows = []
    for lp in (4, 5, 6, 7, 8, 9):
        schemes.append((f"shipped, Likely Pathogenic from {lp} points", {"lp_cut": lp}))
    for pc in (9, 10, 11, 12):
        schemes.append((f"shipped, Pathogenic from {pc} points", {"p_cut": pc}))
    for name, kw in schemes:
        pred = np.array([call(c, nc, i, **kw) for i, (c, nc) in enumerate(zip(codes, nocall))])
        mm = matrix_metrics(matrix_of(pred))
        rows.append({"scheme": name, **{k: mm[k] for k in
                     ("exact", "same_direction", "same_pathogenic", "same_benign", "opposite", "false_benign", "false_pathogenic",
                      "plp_recall", "blb_recall", "vus_recall", "plp_to_vus", "blb_to_vus", "vus_to_plp", "vus_to_blb")}})
    write("combining_rules_replay_v37.csv", rows)
    validate_replay(call, codes, nocall, matrix_of)

    # A lone BP4 under the full point system: the class of call v10 measured at 22 of 36.
    pred = np.array([call(c, nc, i, scheme="points", lb_floor=1) for i, (c, nc) in enumerate(zip(codes, nocall))])
    sig = df["criteria"].map(signature).to_numpy()
    fbmask = np.array([(t in PATH) and (p in ("LB", "B")) for t, p in zip(truth, pred)])
    c = Counter(sig[fbmask])
    write("full_points_false_benign_signatures_v37.csv",
          [{"signature": s, "n": n} for s, n in c.most_common(15)])


def validate_replay(call, codes, nocall, matrix_of) -> None:
    """Check the replay against real runs of the same configuration, cell by cell.

    A replay is a model of `combiner.rs`, and a model that has only been checked against the
    configuration it was fitted on has not been checked. Each real run here differs from v37 in
    the combining step only, so its 5x6 concordance matrix is what the replay must reproduce:
    the Richards 2015 table on both sides (a shipped switch) and the point cutoffs (runs of a
    patched binary, `11_manuscript_runs.py`, group `cutoff`).

    One thing the replay cannot see is the site-level frequency veto, which takes a Likely
    Pathogenic call to VUS where gnomAD rejected the site. It is held at VUS for the 19 variants
    v37 shows it on. A different cutoff makes other variants Likely Pathogenic, and the veto may
    apply to some of them, which shows as a replayed Likely Pathogenic that the run called VUS.
    That is the only difference a faithful replay should show, so the check is that a difference
    is of that shape and not that there is none.
    """
    checks = [("v37", BENCH / "output_v37/concordance_matrix.csv", {}),
              ("v37_sw_table_both_sides", BENCH / "manuscript_runs/v37_sw_table_both_sides/concordance_matrix.csv", {"scheme": "table"})]
    # lp6 is the shipped cutoff run on the patched binary: the control that the patch changes nothing else.
    checks += [(f"v37_cut_lp{n}", BENCH / f"manuscript_runs/v37_cut_lp{n}/concordance_matrix.csv", {"lp_cut": n}) for n in (4, 5, 6, 7, 8, 9)]
    checks += [(f"v37_cut_p{n}", BENCH / f"manuscript_runs/v37_cut_p{n}/concordance_matrix.csv", {"p_cut": n}) for n in (9, 11, 12)]
    rows = []
    for name, path, kw in checks:
        if not path.exists():
            print(f"  (no real run {name}: not validated)")
            continue
        real = read_matrix(path)
        pred = np.array([call(c, nc, i, **kw) for i, (c, nc) in enumerate(zip(codes, nocall))])
        rep_m = matrix_of(pred)
        diff = {t: {p: real[t][p] - rep_m[t][p] for p in CALLS} for t in CLASSES}
        cells = [(t, p) for t in CLASSES for p in CALLS if diff[t][p] != 0]
        moved = sum(max(diff[t][p], 0) for t in CLASSES for p in CALLS)
        # The veto's shape: within each ClinVar class, replayed Likely Pathogenic -> real VUS and nothing else.
        veto_shaped = all(
            diff[t][p] == 0 or p in ("Likely_pathogenic", "VUS") for t in CLASSES for p in CALLS
        ) and all(diff[t]["Likely_pathogenic"] <= 0 and diff[t]["VUS"] == -diff[t]["Likely_pathogenic"] for t in CLASSES)
        mr, mp = matrix_metrics(real), matrix_metrics(rep_m)
        rows.append({
            "run": name, "scheme": ", ".join(f"{k}={v}" for k, v in kw.items()) or "shipped",
            "cells_compared": len(CLASSES) * len(CALLS), "cells_differ": len(cells), "variants_differ": moved,
            "difference_is_veto_shaped": bool(veto_shaped),
            "same_direction_real": mr["same_direction"], "same_direction_replay": mp["same_direction"],
            "opposite_real": mr["opposite"], "opposite_replay": mp["opposite"],
            "plp_recall_real": mr["plp_recall"], "plp_recall_replay": mp["plp_recall"],
        })
        print(f"  replay vs {name}: {len(cells)} of 30 cells differ, {moved} variants, veto-shaped={veto_shaped}")
        if not veto_shaped:
            sys.exit(f"the replay of {name} differs from the real run in a way the site-level veto does not explain: {cells}")
        if name in ("v37", "v37_sw_table_both_sides", "v37_cut_lp6") and cells:
            sys.exit(f"the replay of {name} must reproduce the real run exactly, and does not: {cells}")
    write("replay_validation_v37.csv", rows)


def criteria_table(df: pd.DataFrame) -> None:
    """How often each computable criterion fires in v37, at which strengths, and on which ClinVar classes.

    A criterion's firings are counted per variant (a variant that carries PM2 once counts once). The
    share in its own direction is the fraction of the variants it fired on whose ClinVar class is
    pathogenic-side (for a pathogenic criterion) or benign-side (for a benign one); the share in the
    opposite direction is the fraction on the other side. A criterion firing on a ClinVar VUS is not
    an error, since a VUS is a class that submitters left open.
    """
    pathogenic = {"PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5", "PM6", "PP1", "PP2", "PP3", "PP4", "PP5"}
    rows = []
    split = df["criteria"].map(lambda c: [x for x in c.split("&") if x])
    long = pd.DataFrame({"code": [x for c in split for x in c], "cv3": [cv for c, cv in zip(split, df["cv3"]) for _ in c]})
    long["base"] = long["code"].str.split("_").str[0]
    for base, sub in long.groupby("base"):
        strengths = sub["code"].map(lambda x: x[len(base) + 1:] or "default").value_counts()
        n = len(sub)
        is_p = base in pathogenic
        own, opp = ("P/LP", "B/LB") if is_p else ("B/LB", "P/LP")
        rows.append({
            "criterion": base, "direction": "pathogenic" if is_p else "benign", "variants": n,
            "strengths": "; ".join(f"{k} {int(v):,}" for k, v in strengths.items()),
            "clinvar_plp": int((sub.cv3 == "P/LP").sum()), "clinvar_vus": int((sub.cv3 == "VUS").sum()), "clinvar_blb": int((sub.cv3 == "B/LB").sum()),
            "share_in_own_direction": float((sub.cv3 == own).sum() / n), "share_in_opposite_direction": float((sub.cv3 == opp).sum() / n),
        })
    out = pd.DataFrame(rows)
    order = {c: i for i, c in enumerate(["PVS1", "PS1", "PM1", "PM2", "PM4", "PM5", "PP2", "PP3", "BA1", "BS1", "BS2", "BP1", "BP3", "BP4", "BP7"])}
    out["_o"] = out["criterion"].map(lambda c: order.get(c, 99))
    write("criteria_firing_v37.csv", out.sort_values("_o").drop(columns="_o"))


def pm2_dependence(df: pd.DataFrame) -> None:
    """How many pathogenic-side calls stand on PM2 alone, under each combining rule.

    A call depends on PM2 if it is Pathogenic or Likely Pathogenic with PM2 at Supporting and
    would not be without it. This is what makes the PM2 bar matter: a bar can only move calls
    that need the criterion it gates.
    """
    codes = df["criteria"].map(lambda s: [c for c in s.split("&") if c])
    rows = []
    for name, kw in (("shipped: SVI points on the pathogenic side", {}), ("Richards 2015 table", {"scheme": "table"})):
        veto_free = 0
        dep = Counter()
        pvs_only = 0
        for c, cv in zip(codes, df["clinvar"]):
            if "PM2_Supporting" not in c:
                continue
            a = replay(c, **kw)
            b = replay([x for x in c if x != "PM2_Supporting"], **kw)
            if a in ("P", "LP") and b not in ("P", "LP"):
                dep[TRI[cv]] += 1
                pvs_only += sorted(c) == ["PM2_Supporting", "PVS1"]
        n = sum(dep.values())
        rows.append({"rule": name, "calls_depending_on_pm2": n, "of_which_pvs1_and_pm2_only": pvs_only,
                     "clinvar_plp": dep["P/LP"], "clinvar_vus": dep["VUS"], "clinvar_blb": dep["B/LB"]})
    write("pm2_dependence_v37.csv", rows)


def opposite_cases(df: pd.DataFrame) -> None:
    rev = pd.read_csv(OUT / "round6_rulings.tsv", sep="\t", dtype=str)
    opp = df[((df.cv3 == "P/LP") & (df.fv3 == "B/LB")) | ((df.cv3 == "B/LB") & (df.fv3 == "P/LP"))].copy()
    opp["sig"] = opp["criteria"].map(signature)

    def group(r):
        s = r["sig"]
        if "PVS1" in s:
            return "PVS1 on a ClinVar benign-side variant"
        if "BP7" in s:
            return "BP4+BP7 on a ClinVar pathogenic-side variant"
        if "PM1" in s and r["cv3"] == "B/LB":
            return "PM1 with PP2 and PP3 on a ClinVar benign-side variant"
        if "BS1" in s:
            return "BS1 on a ClinVar pathogenic-side variant"
        if "BP1" in s:
            return "BP1 on a ClinVar pathogenic-side variant"
        return "other"

    opp["group"] = opp.apply(group, axis=1)
    m = opp.merge(rev[["chrom", "pos", "ref", "alt", "ruling", "row"]], on=["chrom", "pos", "ref", "alt"], how="left")
    if m["ruling"].isna().any():
        sys.exit("an opposite-direction call is not one of the reviewer's 64 rows")
    cols = ["group", "gene", "chrom", "pos", "ref", "alt", "stars", "clinvar", "call", "criteria", "consequence", "hgvsc", "hgvsp",
            "gnomad_ac", "gnomad_af", "gnomad_popmax_af", "ruling", "row"]
    write("opposite_direction_v37.csv", m[cols].sort_values(["group", "gene"]))
    t = pd.crosstab(m["group"], m["ruling"])
    t["total"] = t.sum(axis=1)
    t = t.reset_index()
    write("opposite_direction_by_group_v37.csv", t)

    # All 64 rows: what became of them.
    all64 = rev.merge(df[["chrom", "pos", "ref", "alt", "call", "criteria", "cv3", "fv3"]], on=["chrom", "pos", "ref", "alt"])
    all64["still_opposite"] = ((all64.cv3 == "P/LP") & (all64.fv3 == "B/LB")) | ((all64.cv3 == "B/LB") & (all64.fv3 == "P/LP"))
    all64["v37_call"] = all64["call"]
    write("round6_rows_v37.csv", all64[["row", "gene", "chrom", "pos", "ref", "alt", "clinvar_class", "stars",
                                       "fastvep_call_at_round6", "v37_call", "still_opposite", "ruling", "rule", "seen_in_round"]])

    # Which coding rule produced each ruling, for all 64 rows and for the 49 that are still opposite.
    rows = []
    for label, sub in (("all 64 rows of round 6", all64), ("the 49 still opposite in v37", all64[all64.still_opposite])):
        t = pd.crosstab(sub["rule"], sub["ruling"])
        for rule in t.index:
            for ruling in t.columns:
                if t.loc[rule, ruling]:
                    rows.append({"set": label, "rule": rule, "ruling": ruling, "rows": int(t.loc[rule, ruling])})
    write("round6_rules_v37.csv", rows)


def gene_level(df: pd.DataFrame) -> None:
    df = df.copy()
    g = df.groupby("gene")
    rows = []
    for gene, sub in g:
        n_plp = int((sub.cv3 == "P/LP").sum())
        n_blb = int((sub.cv3 == "B/LB").sum())
        rows.append({
            "gene": gene, "n": len(sub), "clinvar_plp": n_plp, "clinvar_blb": n_blb,
            "plp_called_plp": int(((sub.cv3 == "P/LP") & (sub.fv3 == "P/LP")).sum()),
            "plp_called_vus": int(((sub.cv3 == "P/LP") & (sub.fv3 == "VUS")).sum()),
            "blb_called_blb": int(((sub.cv3 == "B/LB") & (sub.fv3 == "B/LB")).sum()),
            "blb_called_vus": int(((sub.cv3 == "B/LB") & (sub.fv3 == "VUS")).sum()),
            "opposite": int((((sub.cv3 == "P/LP") & (sub.fv3 == "B/LB")) | ((sub.cv3 == "B/LB") & (sub.fv3 == "P/LP"))).sum()),
            "clingen_tags": sub["clingen_tags"].iloc[0],
        })
    gdf = pd.DataFrame(rows).sort_values("n", ascending=False)
    write("gene_level_v37.csv", gdf)

    # Concentration: how much of ClinVar P/LP -> VUS sits in how many genes.
    gl = gdf.sort_values("plp_called_vus", ascending=False)
    tot = gl["plp_called_vus"].sum()
    cum = gl["plp_called_vus"].cumsum() / tot
    rows = []
    for k in (10, 25, 50, 100, 250, 500, 1000):
        rows.append({"top_genes": k, "share_of_plp_to_vus": float(cum.iloc[k - 1]),
                     "genes_with_any": int((gl["plp_called_vus"] > 0).sum())})
    write("plp_to_vus_gene_concentration_v37.csv", rows)

    # Per-gene concordance heterogeneity among genes with enough P/LP variants to say anything.
    big = gdf[gdf.clinvar_plp >= 100].copy()
    big["plp_recall"] = big["plp_called_plp"] / big["clinvar_plp"]
    write("plp_recall_by_gene_v37.csv", big.sort_values("plp_recall"))


def read_clingen(path: Path) -> pd.DataFrame:
    """A ClinGen gene-disease validity export, without its preamble and separator lines."""
    import io

    lines = path.read_text(encoding="utf-8").splitlines()
    head = next(i for i, ln in enumerate(lines) if ln.startswith('"GENE SYMBOL"'))
    return pd.read_csv(io.StringIO("\n".join(ln for ln in lines[head:] if not ln.startswith('"+++'))))


def shipped_mechanisms() -> dict[str, str]:
    """The curated disease mechanisms the classifier ships, read from the source of `config.rs`."""
    src = (ROOT / "crates/fastvep-classification/src/config.rs").read_text()
    body = src[src.index("fn default_gene_mechanisms"):]
    body = body[: body.index("collect()")]
    return dict(re.findall(r'\(\s*"([A-Z0-9-]+)"\s*,\s*"([A-Za-z_]+)"\s*\)', body))


def pvs1_gap(df: pd.DataFrame) -> None:
    """ClinVar P/LP nonsense and frameshift variants that never reach PVS1, and the genes they sit in.

    PVS1 needs an established loss-of-function mechanism, which the classifier reads from gnomAD
    constraint (pLI at least 0.9 or LOEUF at most 0.35), a curated mechanism, or a ClinGen
    gene-disease curation. A recessive gene that tolerates heterozygous loss in gnomAD and is
    absent from the ClinGen file meets none of them, so a stop-gained variant in it carries no
    PVS1 at all.
    """
    nl = df[df.consequence.isin(["stop_gained", "frameshift_variant"]) & (df.cv3 == "P/LP")].copy()
    nl["pvs1"] = nl["criteria"].str.contains("PVS1")
    nl["pli_f"] = pd.to_numeric(nl["pli"], errors="coerce")
    nl["loeuf_f"] = pd.to_numeric(nl["loeuf"], errors="coerce")
    nl["constrained"] = (nl.pli_f >= 0.9) | (nl.loeuf_f <= 0.35)
    g = nl.groupby("gene").agg(n=("call", "size"), with_pvs1=("pvs1", "sum"), pli=("pli_f", "first"), loeuf=("loeuf_f", "first"),
                               clingen_tags=("clingen_tags", "first"))
    g["without_pvs1"] = g.n - g.with_pvs1
    write("pvs1_gap_genes_v37.csv", g.reset_index().sort_values("without_pvs1", ascending=False).head(40))
    miss = nl[~nl.pvs1]
    rows = [
        {"measure": "ClinVar P/LP stop-gained and frameshift variants", "n": len(nl)},
        {"measure": "with a PVS1 criterion (any strength)", "n": int(nl.pvs1.sum())},
        {"measure": "without PVS1", "n": len(miss)},
        {"measure": "without PVS1, genes", "n": int(miss.gene.nunique())},
        # `gene` is the first symbol of ClinVar's GENEINFO; the annotation names the gene of the picked transcript
        {"measure": "without PVS1, genes by the symbol of the annotation", "n": int(miss.symbol.nunique())},
        {"measure": "without PVS1, variants whose ClinVar gene differs from the symbol of the annotation", "n": int((miss.gene != miss.symbol).sum())},
        {"measure": "without PVS1, variants in genes constrained by gnomAD", "n": int(miss.constrained.sum())},
        {"measure": "without PVS1, variants in genes with no established ClinGen curation", "n": int((miss.clingen_tags == "none").sum())},
        {"measure": "without PVS1, in the 15 genes with most such variants", "n": int(g.sort_values("without_pvs1", ascending=False).head(15).without_pvs1.sum())},
        {"measure": "without PVS1, in CFTR (stop-gained and frameshift)", "n": int((miss.gene == "CFTR").sum())},
    ]
    nul = df[df.consequence.isin(["stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant"])]
    c = nul[nul.gene == "CFTR"]
    rows += [{"measure": "CFTR null variants in the benchmark (stop, frameshift, canonical splice)", "n": len(c)},
             {"measure": "  of them called VUS", "n": int((c.fv3 == "VUS").sum())},
             {"measure": "  of them carrying any PVS1", "n": int(c.criteria.str.contains("PVS1").sum())}]

    # Why PVS1 is absent, in the order PVS1 asks its questions (`criteria/pvs1.rs`): the validity gate
    # (every ClinGen curation of the gene is Limited, Disputed, Refuted or No Known Disease
    # Relationship), the mechanism gate (the shipped table says the gene acts by gain of function or
    # dominant negatively), then the three routes to "loss of function is the mechanism" (gnomAD
    # constraint, a curated loss-of-function mechanism, a ClinGen curation of Moderate or better).
    cg = BENCH / "sa_sources/clingen_gene_validity.csv"
    weak = {"Limited", "Disputed", "Refuted", "No Known Disease Relationship"}
    classes = read_clingen(cg).groupby("GENE SYMBOL")["CLASSIFICATION"].agg(set)
    weak_only = {gn for gn, k in classes.items() if k <= weak}
    established = {gn for gn, k in classes.items() if k - weak}
    mech = shipped_mechanisms()
    non_lof = {gn for gn, m in mech.items() if "LOF" not in m.upper()}
    lof = {gn for gn, m in mech.items() if "LOF" in m.upper()}

    def reason(r):
        if r.gene in weak_only:
            return "every ClinGen curation of the gene is weak (validity gate)"
        if r.gene in non_lof:
            return "the curated mechanism is gain of function or dominant negative (mechanism gate)"
        if r.constrained or r.gene in lof or r.gene in established:
            return "a route to loss of function is met but PVS1 is not (other)"
        return "no route to loss of function as the mechanism"

    miss = miss.copy()
    miss["reason"] = miss.apply(reason, axis=1)
    out = []
    for why, sub in miss.groupby("reason"):
        out.append({"reason": why, "variants": len(sub), "genes": int(sub.gene.nunique()), "variants_in_constrained_genes": int(sub.constrained.sum())})
    write("pvs1_gap_reasons_v37.csv", sorted(out, key=lambda r: -r["variants"]))
    cons = miss[miss.constrained]
    # The conditions ClinVar's own records name for these variants, from the April VCF (CLNDN).
    wanted = {(r.chrom, r.pos, r.ref, r.alt): r.gene for r in cons.itertuples()}
    names: dict[str, Counter] = {}
    with (BENCH / "clinvar_2star.vcf").open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            gn = wanted.get((c[0].removeprefix("chr"), c[1], c[3], c[4])) or wanted.get((c[0], c[1], c[3], c[4]))
            if gn:
                info = dict(kv.split("=", 1) for kv in c[7].split(";") if "=" in kv)
                for d in set(info.get("CLNDN", "").split("|")):
                    names.setdefault(gn, Counter())[d] += 1

    def top_condition(gn):
        cnt = {d: k for d, k in names.get(gn, Counter()).items() if d != "not_provided"}
        d = max(cnt, key=cnt.get) if cnt else ""
        return d.replace("_", " "), cnt.get(d, 0)

    write("pvs1_gap_constrained_v37.csv", [
        {"gene": gn, "variants": int(k), "reason": cons[cons.gene == gn].reason.iloc[0],
         "clingen_curations": "; ".join(sorted(f"{d} ({c})" for d, c in
                                               read_clingen(cg).query("`GENE SYMBOL` == @gn")[["DISEASE LABEL", "CLASSIFICATION"]].itertuples(index=False))),
         "clinvar_most_named_condition": top_condition(gn)[0], "variants_naming_it": top_condition(gn)[1]}
        for gn, k in cons.gene.value_counts().items()])
    rows += [{"measure": "without PVS1, constrained genes", "n": int(cons.gene.nunique())},
             {"measure": "without PVS1, validity gate", "n": int((miss.reason.str.contains("validity")).sum())},
             {"measure": "without PVS1, mechanism gate", "n": int((miss.reason.str.contains("mechanism gate")).sum())},
             {"measure": "without PVS1, no route to loss of function", "n": int((miss.reason.str.startswith("no route")).sum())},
             {"measure": "without PVS1, no route, genes", "n": int(miss[miss.reason.str.startswith("no route")].gene.nunique())}]
    write("pvs1_gap_v37.csv", rows)


def baselines(df: pd.DataFrame) -> None:
    """What simple rules agree with ClinVar on, to put the agreement figure in context.

    The rules read the inputs of some ACMG criteria (population frequency, consequence class and SpliceAI)
    and apply no combining rule, gene-level gate or safeguard. Four rules, each added to the one before and
    none tuned on this benchmark. The bars are round
    numbers: 1e-3 is the median BA1 bar of the published VCEP specifications (Section 3.3), 1e-4 is
    a rarity bar of the same order as the PM2 bars, and 0.2 is SpliceAI's pathogenic bar. A rule is
    applied only to a variant that no earlier rule has called, and a reference-agreement record
    (ALT is a dot) is never concordant, as in the scorer. The population maximum is the highest allele
    frequency of the ten groups the annotation carries, "other" and "remaining" among them, which is a
    point estimate and not the filtering allele frequency that BA1 and BS1 test.

      all VUS   call every variant VUS
      R1        population-maximum allele frequency of 1e-3 or more: benign-side
      R2        a stop-gained, frameshift, canonical splice or start-lost variant that is absent from
                gnomAD or below 1e-4: pathogenic-side
      R3        a synonymous, intronic, UTR, flanking or non-coding variant with a SpliceAI score
                below 0.2: benign-side. Intronic includes the polypyrimidine tract and flanking means
                upstream or downstream of a gene, because the consequence column holds the most severe
                term and VEP ranks those two above intron_variant and below the coding terms
    """
    af = pd.to_numeric(df["gnomad_popmax_af"].replace("", np.nan), errors="coerce")
    ds = pd.to_numeric(df["spliceai_max"].replace("", np.nan), errors="coerce")
    null = df["consequence"].isin(["stop_gained", "frameshift_variant", "splice_donor_variant",
                                   "splice_acceptor_variant", "start_lost"])
    quiet = df["consequence"].isin(["synonymous_variant", "intron_variant", "3_prime_UTR_variant",
                                    "5_prime_UTR_variant", "upstream_gene_variant", "downstream_gene_variant",
                                    "non_coding_transcript_exon_variant", "splice_polypyrimidine_tract_variant"])
    nocall = df["call"] == "NoCall"
    cv = df["cv3"]

    def score(name: str, rule: str, pred: pd.Series) -> dict:
        pred = pred.where(~nocall, "NoCall")
        same = int((pred == cv).sum())
        missed = int(((cv == "P/LP") & (pred == "B/LB")).sum())
        fp = int(((cv == "B/LB") & (pred == "P/LP")).sum())
        cp, cb = int((pred == "P/LP").sum()), int((pred == "B/LB").sum())
        return {"rule": name, "definition": rule, "same_direction": same, "same_direction_share": same / len(df),
                "opposite": missed + fp, "false_benign": missed, "false_pathogenic": fp,
                "called_plp": cp, "called_blb": cb,
                "plp_calls_that_clinvar_calls_plp": int(((pred == "P/LP") & (cv == "P/LP")).sum()),
                "blb_calls_that_clinvar_calls_blb": int(((pred == "B/LB") & (cv == "B/LB")).sum()),
                "plp_recall": ((pred == "P/LP") & (cv == "P/LP")).sum() / (cv == "P/LP").sum(),
                "plp_concordant_share": ((pred == "P/LP") & (cv == "P/LP")).sum() / cp if cp else "",
                "blb_recall": ((pred == "B/LB") & (cv == "B/LB")).sum() / (cv == "B/LB").sum(),
                "blb_concordant_share": ((pred == "B/LB") & (cv == "B/LB")).sum() / cb if cb else ""}

    v = pd.Series("VUS", index=df.index)
    rows = [score("v37 (the classifier)", "the ACMG/AMP classification of this paper", df["fv3"])]
    rows.append(score("all VUS", "every variant VUS", v))
    r1 = v.copy()
    r1[af >= 1e-3] = "B/LB"
    rows.append(score("R1", "popmax AF at least 1e-3: benign-side", r1))
    r2 = r1.copy()
    r2[null & (af.isna() | (af < 1e-4)) & (r2 == "VUS")] = "P/LP"
    rows.append(score("R1 + R2", "R1; a null variant absent from gnomAD or below 1e-4: pathogenic-side", r2))
    r3 = r2.copy()
    r3[quiet & (ds < 0.2) & (r3 == "VUS")] = "B/LB"
    rows.append(score("R1 + R2 + R3", "R1, R2; a synonymous, intronic (polypyrimidine tract included), untranslated, flanking or non-coding variant with SpliceAI below 0.2: benign-side", r3))
    for bar in (5e-3, 1e-4):
        x = v.copy()
        x[af >= bar] = "B/LB"
        rows.append(score(f"R1 at {bar:g}", f"popmax AF at least {bar:g}: benign-side", x))
    if rows[0]["same_direction"] != int((df["cv3"] == df["fv3"]).sum()):
        sys.exit("the baseline scorer does not reproduce the classifier's agreement")
    write("baselines_v37.csv", rows)

    # The same three rules at other frequency bars, so that the family is a curve and not three points:
    # what agreement the rules reach at each number of contradictions, against v37's. The bars run from 5e-2, where the first
    # rule has fewer contradictions than v37, to 1e-4. The last two rules are also given with no frequency rule at all (the bar
    # is empty), which is the least that the third can have to do with a population frequency.
    curve = []

    def curve_row(label: str, bar, pred: pd.Series) -> dict:
        r = score(label, "", pred)
        return {"rule": label, "benign_frequency_bar": bar, "same_direction_share": r["same_direction_share"],
                "opposite": r["opposite"], "false_benign": r["false_benign"], "false_pathogenic": r["false_pathogenic"],
                "plp_recall": r["plp_recall"], "plp_concordant_share": r["plp_concordant_share"],
                "blb_recall": r["blb_recall"], "blb_concordant_share": r["blb_concordant_share"]}

    def family(bar):
        a = v.copy()
        if bar is not None:
            a[af >= bar] = "B/LB"
        b = a.copy()
        b[null & (af.isna() | (af < 1e-4)) & (b == "VUS")] = "P/LP"
        c = b.copy()
        c[quiet & (ds < 0.2) & (c == "VUS")] = "B/LB"
        return a, b, c

    for bar in (5e-2, 2e-2, 1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4):
        a, b, c = family(bar)
        for label, pred in (("R1", a), ("R1 + R2", b), ("R1 + R2 + R3", c)):
            curve.append(curve_row(label, bar, pred))
    _, b0, c0 = family(None)
    curve.append(curve_row("R2", "", b0))
    curve.append(curve_row("R2 + R3", "", c0))
    write("baselines_curve_v37.csv", curve)

    # The second rule against the PVS1 gap: its pathogenic-side calls that v37 leaves as VUS, and whether v37 gave them PVS1.
    pvs1 = df["criteria"].str.contains("PVS1")
    r2_correct = (b0 == "P/LP") & (cv == "P/LP")
    r2_wrong = (b0 == "P/LP") & (cv == "B/LB")
    left = r2_correct & (df["fv3"] == "VUS")
    wrong_left = r2_wrong & (df["fv3"] == "VUS")
    write("baselines_r2_vs_v37.csv", [
        {"measure": "pathogenic-side calls of the second rule that ClinVar calls pathogenic-side", "n": int(r2_correct.sum())},
        {"measure": "  that v37 leaves as VUS", "n": int(left.sum())},
        {"measure": "    with a PVS1 criterion in v37, at a lower strength or held back", "n": int((left & pvs1).sum())},
        {"measure": "    with no PVS1 criterion in v37", "n": int((left & ~pvs1).sum())},
        {"measure": "pathogenic-side calls of the second rule that ClinVar calls benign-side", "n": int(r2_wrong.sum())},
        {"measure": "  that v37 leaves as VUS", "n": int(wrong_left.sum())},
        {"measure": "    with no PVS1 criterion in v37", "n": int((wrong_left & ~pvs1).sum())}])

    # What each tier of the call is worth: of the variants called at each of the five tiers, how many
    # ClinVar puts in the same direction, in VUS, and in the opposite direction.
    out = []
    for call in CLASSES:
        sub = df[df["call"] == call]
        d = TRI[call]
        opp = "B/LB" if d == "P/LP" else "P/LP" if d == "B/LB" else None
        k_opp = int((sub["cv3"] == opp).sum()) if opp else None
        lo_o, hi_o = clopper_pearson(k_opp, len(sub)) if opp else ("", "")
        out.append({"call": call, "n": len(sub),
                    "clinvar_same_direction": int((sub["cv3"] == d).sum()),
                    "clinvar_vus": int((sub["cv3"] == "VUS").sum()) if d != "VUS" else "",
                    "clinvar_opposite": k_opp if opp else "",
                    "clinvar_exact": int((sub["clinvar"] == call).sum()),
                    "opposite_rate": k_opp / len(sub) if opp else "", "opposite_ci_lo": lo_o, "opposite_ci_hi": hi_o})
    write("precision_by_tier_v37.csv", out)


def gap_runs() -> None:
    """What the gap runs of `11_manuscript_runs.py` do to the variants that carry no PVS1 in v37.

    Three runs, each a full benchmark of the final classifier with one input changed: the ClinGen file
    with the curations the 30 April export lacked put back (`15_clingen_complete_export.py`), PVS1
    licensed by Orphanet's loss-of-function label, and both. For each: how many of the 5,184
    stop-gained and frameshift variants of ClinVar P/LP that lack PVS1 in v37 gain it, how many of
    those are then called pathogenic-side, and every change of call that the run makes, by ClinVar
    class, with the variants that become false-pathogenic one by one.
    """
    base = pd.read_csv(PER_VARIANT, sep="\t", dtype=str, keep_default_na=False)
    key = ["chrom", "pos", "ref", "alt"]
    for d in (base,):
        d["cv3"] = d["clinvar"].map(TRI)
        d["fv3"] = d["call"].map(TRI)
        d["pvs1"] = d["criteria"].str.contains("PVS1")
    gap = base[base.consequence.isin(["stop_gained", "frameshift_variant"]) & (base.cv3 == "P/LP") & ~base.pvs1]
    names = {"v37_gap_clingen_complete": "ClinGen file with the omitted curations put back",
             "v37_gap_orphanet_lof": "PVS1 licensed by Orphanet's loss-of-function label",
             "v37_gap_both": "both"}
    rows, new_fp, per_gene = [], [], []
    for run, label in names.items():
        x = pd.read_csv(BENCH / f"manuscript_runs/{run}_per_variant.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
        x["cv3"] = x["clinvar"].map(TRI)
        x["fv3"] = x["call"].map(TRI)
        x["pvs1"] = x["criteria"].str.contains("PVS1")
        g = gap[key + ["gene"]].merge(x[key + ["pvs1", "fv3"]], on=key)
        b = base[key + ["gene", "cv3", "fv3", "call", "pvs1", "criteria", "clinvar", "consequence", "gnomad_popmax_af"]].merge(
            x[key + ["call", "fv3", "pvs1", "criteria"]], on=key, suffixes=("", "_new"))
        gained = b[~b.pvs1 & b.pvs1_new]
        changed = b[b.call != b.call_new]
        cftr = x[(x.gene == "CFTR") & x.consequence.isin(["stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant"])]
        four = x[x.stars == "4"]
        rows.append({
            "run": run, "label": label,
            "gap_variants": len(g), "gap_variants_gaining_pvs1": int(g.pvs1.sum()),
            "gap_variants_called_plp": int((g.fv3 == "P/LP").sum()),
            "genes_gaining_pvs1": int(g[g.pvs1].gene.nunique()),
            "pvs1_gained_clinvar_plp": int((gained.cv3 == "P/LP").sum()),
            "pvs1_gained_clinvar_vus": int((gained.cv3 == "VUS").sum()),
            "pvs1_gained_clinvar_blb": int((gained.cv3 == "B/LB").sum()),
            "calls_changed": len(changed),
            "plp_vus_to_plp": int(((changed.cv3 == "P/LP") & (changed.fv3 == "VUS") & (changed.fv3_new == "P/LP")).sum()),
            "vus_to_plp": int(((changed.cv3 == "VUS") & (changed.fv3_new == "P/LP")).sum()),
            "blb_to_plp": int(((changed.cv3 == "B/LB") & (changed.fv3_new == "P/LP")).sum()),
            "blb_lost_a_call": int(((changed.cv3 == "B/LB") & (changed.fv3 == "B/LB") & (changed.fv3_new == "VUS")).sum()),
            "cftr_null_variants": len(cftr), "cftr_null_with_pvs1": int(cftr.pvs1.sum()),
            "cftr_null_called_plp": int((cftr.fv3 == "P/LP").sum()),
            "four_star_called_plp": int((four.fv3 == "P/LP").sum()), "four_star": len(four),
        })
        fp = changed[(changed.cv3 == "B/LB") & (changed.fv3_new == "P/LP")]
        for r in fp.itertuples():
            new_fp.append({"run": run, "gene": r.gene, "clinvar": r.clinvar, "call_in_v37": r.call, "call_in_run": r.call_new,
                           "criteria_in_run": r.criteria_new, "consequence": r.consequence, "popmax_af": r.gnomad_popmax_af})
        top = g[g.pvs1].gene.value_counts().head(12)
        for gn, k in top.items():
            per_gene.append({"run": run, "gene": gn, "gap_variants_gaining_pvs1": int(k)})
    write("pvs1_gap_closure_v37.csv", rows)
    write("pvs1_gap_closure_new_false_pathogenic_v37.csv", new_fp)
    write("pvs1_gap_closure_genes_v37.csv", per_gene)
    # How many curations the April export lacked, from the file `15_clingen_complete_export.py` wrote. The rows it
    # added are not all curations the April export lacked: some re-label an assertion the export does list, and one
    # is dated after the export was downloaded, so the summary keeps the three apart.
    cur = pd.read_csv(OUT / "clingen_missing_curations.csv")
    lacked = cur[cur.status == "lacked"]
    april = read_clingen(BENCH / "sa_sources/clingen_gene_validity.csv")
    # The one row dated on the day of the export is in a gene whose calls are checked here to be the same with and
    # without the completed file, because the runs that read the completed file carry it.
    day_gene = cur[cur.status == "export day"].gene.iloc[0]
    v37 = pd.read_csv(PER_VARIANT, sep="\t", dtype=str, keep_default_na=False)
    done = pd.read_csv(BENCH / "manuscript_runs/v37_gap_clingen_complete_per_variant.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    both = v37[v37.gene == day_gene].merge(done[key + ["call", "criteria"]], on=key, suffixes=("", "_run"))
    write("clingen_missing_summary_v37.csv", [
        {"measure": "curations in the 30 April export", "n": len(april)},
        {"measure": "rows of the 11 October export dated on or before 30 April whose gene and disease pair the April export lacks", "n": len(cur)},
        {"measure": "  of them re-labelled (the April export lists the same assertion under another disease identifier)", "n": int((cur.status == "relabel").sum())},
        {"measure": "  of them dated 30 April itself, after the April file was downloaded", "n": int((cur.status == "export day").sum())},
        {"measure": "curations dated before 30 April that the export lacked", "n": len(lacked)},
        {"measure": "  in genes", "n": int(lacked.gene.nunique())},
        {"measure": "  of Definitive, Strong or Moderate strength", "n": int(lacked.established.sum())},
        {"measure": "  genes with such a curation", "n": int(lacked[lacked.established].gene.nunique())},
        {"measure": "  genes with such a curation and no other curation in the April export",
         "n": int(len(set(lacked[lacked.established].gene) - set(april["GENE SYMBOL"])))},
        {"measure": "  classified in 2025 or 2026", "n": int((lacked.classified.str[:4].astype(int) >= 2025).sum())},
        {"measure": "  classified from January to April 2026", "n": int((lacked.classified.str[:7] >= "2026-01").sum())},
        {"measure": "  classified before 2024", "n": int((lacked.classified.str[:4].astype(int) < 2024).sum())},
        {"measure": "curations in the April export classified in 2025 or 2026", "n": int((april["CLASSIFICATION DATE"].str[:4].astype(int) >= 2025).sum())},
        {"measure": "curations in the April export classified from January to April 2026", "n": int((april["CLASSIFICATION DATE"].str[:7] >= "2026-01").sum())},
        {"measure": "curations in the April export classified before 2024", "n": int((april["CLASSIFICATION DATE"].str[:4].astype(int) < 2024).sum())},
        # What a random omission would leave: the lacked curations with the dates of the ones the April export lists.
        {"measure": "  expected in 2025 or 2026 if the omission were independent of the date",
         "n": len(lacked) * int((april["CLASSIFICATION DATE"].str[:4].astype(int) >= 2025).sum()) / len(april)},
        {"measure": "  expected from January to April 2026 if the omission were independent of the date",
         "n": len(lacked) * int((april["CLASSIFICATION DATE"].str[:7] >= "2026-01").sum()) / len(april)},
        {"measure": "  expected before 2024 if the omission were independent of the date",
         "n": len(lacked) * int((april["CLASSIFICATION DATE"].str[:4].astype(int) < 2024).sum()) / len(april)},
        {"measure": f"variants of {day_gene}, the gene of the row dated 30 April itself, in the benchmark", "n": len(both)},
        {"measure": "  whose call or criteria differ in the run that reads the completed file",
         "n": int(((both.call != both.call_run) | (both.criteria != both.criteria_run)).sum())}])


def hook_runs() -> None:
    """The variants whose call becomes opposite-direction in the runs that test BP1's veto and BP7.

    Each is a run of the final classifier that differs from v37 in one setting, with its per-variant
    table made by `12_per_variant_table.py`. A call that was already opposite in v37 is not listed.
    """
    base = pd.read_csv(PER_VARIANT, sep="\t", dtype=str, keep_default_na=False)
    key = ["chrom", "pos", "ref", "alt"]
    base["cv3"] = base["clinvar"].map(TRI)
    rows = []
    for run in ("v37_hk_bp1_loo", "v37_hk_bp1_loo_veto1", "v37_hk_bp1_loo_veto2", "v37_sw_phylop_1000", "v37_hk_bp7_walker"):
        x = pd.read_csv(BENCH / f"manuscript_runs/{run}_per_variant.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
        x["fv3"] = x["call"].map(TRI)
        b = base[key + ["gene", "cv3", "call", "criteria", "consequence", "hgvsc", "gnomad_popmax_af"]].merge(
            x[key + ["call", "fv3", "criteria"]], on=key, suffixes=("", "_run"))
        was_opp = ((b.cv3 == "P/LP") & (b.call.map(TRI) == "B/LB")) | ((b.cv3 == "B/LB") & (b.call.map(TRI) == "P/LP"))
        now_opp = ((b.cv3 == "P/LP") & (b.fv3 == "B/LB")) | ((b.cv3 == "B/LB") & (b.fv3 == "P/LP"))
        for r in b[now_opp & ~was_opp].itertuples():
            rows.append({"run": run, "gene": r.gene, "clinvar": TRI[r.cv3] if r.cv3 in TRI else r.cv3, "call_in_v37": r.call, "call_in_run": r.call_run,
                         "criteria_in_v37": r.criteria, "criteria_in_run": r.criteria_run, "consequence": r.consequence,
                         "popmax_af": r.gnomad_popmax_af})
    write("hook_new_opposite_v37.csv", rows)

    # The opposite-direction calls of v37 that a stricter count for BP1's veto removes: the two BP1 calls of Section 3.6 and no other,
    # with and without the variant's own record in the count.
    resolved = []
    for run in ("v37_sw_bp1veto_1", "v37_sw_bp1veto_2", "v37_hk_bp1_loo_veto1", "v37_hk_bp1_loo_veto2"):
        x = pd.read_csv(BENCH / f"manuscript_runs/{run}_per_variant.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
        b = base[key + ["gene", "cv3", "call", "criteria"]].merge(x[key + ["call", "criteria"]], on=key, suffixes=("", "_run"))
        b["fv3"] = b["call_run"].map(TRI)
        was_opp = ((b.cv3 == "P/LP") & (b.call.map(TRI) == "B/LB")) | ((b.cv3 == "B/LB") & (b.call.map(TRI) == "P/LP"))
        now_opp = ((b.cv3 == "P/LP") & (b.fv3 == "B/LB")) | ((b.cv3 == "B/LB") & (b.fv3 == "P/LP"))
        for r in b[was_opp & ~now_opp].itertuples():
            resolved.append({"run": run, "gene": r.gene, "clinvar": r.cv3, "call_in_v37": r.call, "call_in_run": r.call_run,
                             "criteria_in_v37": r.criteria, "criteria_in_run": r.criteria_run})
        if (now_opp & ~was_opp).any():
            sys.exit(f"{run} makes a call opposite-direction that was not")
    write("bp1_veto_count_resolved_v37.csv", resolved)

    # What the veto counts for the calls that depend on it: the pathogenic missense records that the protein index holds for the
    # gene, and whether the variant's own is among them. The two BP1 calls of v37 are calls that the veto would have prevented at
    # a stricter count, and the three that the measurement build with the own record left out adds are calls that the own
    # record was keeping from BP1.
    index = read_protein_index(BENCH / "sa_db/clinvar_protein.oga")
    own_rows = []
    for r in base[(base.cv3 == "P/LP") & (base.call.map(TRI) == "B/LB") & base.criteria.str.contains("BP1")].itertuples():
        own_rows.append(_index_row(index, r.gene, r.hgvsp, "opposite-direction in v37, with BP1"))
    loo = pd.read_csv(BENCH / "manuscript_runs/v37_hk_bp1_loo_per_variant.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    x = loo.merge(base[key + ["call", "cv3"]], on=key, suffixes=("", "_v37"))
    for r in x[(x.cv3 == "P/LP") & (x.call.map(TRI) == "B/LB") & (x.call_v37.map(TRI) != "B/LB")].itertuples():
        own_rows.append(_index_row(index, r.gene, r.hgvsp, "opposite-direction only when the veto leaves its own record out"))
    write("bp1_index_counts_v37.csv", own_rows)


def provenance_summary() -> None:
    """What produced each run of this paper, read from the provenance files in the repository."""
    rows = []
    # the runs on the April set; the held-out runs on the October release have provenance files too, and are not among them
    april_runs = set(pd.read_csv(RES / "v37_manuscript/runs.tsv", sep="\t", usecols=["name"])["name"])
    for f in sorted((RES / "v37_manuscript/provenance").glob("*.txt")):
        if f.stem not in april_runs:
            continue
        kv = dict(line.split(": ", 1) for line in f.read_text().splitlines() if ": " in line)
        rows.append({"run": f.stem, "binary_sha256_prefix": kv.get("fastvep_sha256", "")[:8], "repo_commit": kv.get("repo_commit", "")[:7],
                     "repo_dirty_files": kv.get("repo_dirty_files", ""), "annotation_seconds": int(kv.get("annotation_seconds", 0))})
    write("provenance_v37.csv", rows)


def trajectory() -> None:
    """One row per run. v7 and later are read from the run's own matrix, v1 from the matrix stored for it,
    and v2 to v6 are the figures METHODS.md recorded at the time, because their raw outputs were overwritten.

    Every rate is over the 673,660 records of the benchmark. v1 annotated only 627,376 of them (its
    annotation phase wrote JSON, and 46,285 records were not annotated), and the 54.7% recorded at the
    time is over those; over all records it is 51.0%. The matrix stored for v1 is itself reconstructed from
    the documentation of the run (`results/v1/README.md`), and holds 31 opposite-direction calls where the
    run record says 36. v17 is not a step on the way: it is the opt-in run with the published VCEP bars,
    which v18 does not carry, so it is marked and left out of every range the text quotes.
    """
    rows = []

    def committed(m):
        """Variants the run called Pathogenic, Likely pathogenic, Likely benign or Benign."""
        return sum(m[t][c] for t in CLASSES for c in ("Pathogenic", "Likely_pathogenic", "Likely_benign", "Benign"))

    v1 = read_matrix(RES / "v1/concordance_matrix.csv")
    m1 = matrix_metrics(v1)
    n_all = 673_660
    plp_total = sum(sum(v1[t].values()) for t in ("Pathogenic", "Likely_pathogenic"))
    early = [
        # name, same-direction % as recorded, exact % as recorded, opposite count as recorded (None = not recorded), source
        ("v2", 54.9, None, None, "METHODS.md: PhyloP, SpliceAI, ClinGen loaded, two wiring bugs latent"),
        ("v4", 63.7, None, None, "METHODS.md: SpliceAI key and PhyloP route fixed, VCF output"),
        ("v5", 65.1, None, None, "METHODS.md: indel allele matching fixed in the scorer"),
        ("v6", 70.3, 56.8, 338, "RUN_VERSIONS.md: PM2 when absent, BP4 splice gate"),
    ]
    rows.append({"run": "v1", "order": 1, "source": "matrix stored for the run, reconstructed from its documentation",
                 "note": "RUN_VERSIONS.md: baseline, no PhyloP, SpliceAI or ClinGen; annotated 627,376 of the 673,660 records, "
                         "54.7% of those in the same direction; the run record says 36 opposite-direction calls, the stored matrix 31",
                 "same_direction_rate": m1["same_direction"] / n_all, "exact_rate": m1["exact"] / n_all,
                 "opposite": m1["opposite"], "false_benign": m1["false_benign"], "false_pathogenic": m1["false_pathogenic"],
                 "plp_recall": m1["same_pathogenic"] / (79_823 + 13_989), "blb_recall": m1["same_benign"] / (128_038 + 156_512),
                 "n": m1["n"], "denominator": n_all, "mainline": True,
                 "same_direction_rate_of_annotated": m1["same_direction"] / m1["n"], "committed": committed(v1)})
    for name, same, exact, opp, src in early:
        rows.append({"run": name, "order": int(name[1:]), "source": "recorded at the time (raw output not preserved)",
                     "note": src, "same_direction_rate": same / 100, "exact_rate": exact / 100 if exact else "",
                     "opposite": opp if opp is not None else "", "false_benign": "", "false_pathogenic": "",
                     "plp_recall": "", "blb_recall": "", "n": "", "denominator": "", "mainline": True,
                     "same_direction_rate_of_annotated": "", "committed": ""})
    paths = {7: RES / "v7/concordance_matrix.csv", 8: RES / "v8/concordance_matrix.csv"}
    for v in range(9, 38):
        p = BENCH / f"output_v{v}/concordance_matrix.csv"
        if p.exists():
            paths[v] = p
    # v9 is reported leave-one-out, which is what the run directory holds.
    for v, p in sorted(paths.items()):
        mat = read_matrix(p)
        mm = matrix_metrics(mat)
        rows.append({
            "run": f"v{v}", "order": v, "source": "scored matrix on disk",
            "note": "opt-in run with the published VCEP bars, not a step on the way (v18 does not carry it)" if v == 17 else "",
            "same_direction_rate": mm["same_direction"] / mm["n"], "exact_rate": mm["exact"] / mm["n"],
            "opposite": mm["opposite"], "false_benign": mm["false_benign"], "false_pathogenic": mm["false_pathogenic"],
            "plp_recall": mm["plp_recall"], "blb_recall": mm["blb_recall"], "n": mm["n"], "denominator": mm["n"],
            "mainline": v != 17, "same_direction_rate_of_annotated": "", "committed": committed(mat),
        })
    write("trajectory.csv", sorted(rows, key=lambda r: r["order"]))


def runs_tables() -> None:
    if not RUNS.exists():
        print("  (no runs.tsv yet: skipping sweeps and ablations)")
        return
    r = pd.read_csv(RUNS, sep="\t")
    rows = []

    def mat_of(row):
        return {t: {p: int(row[f"cm_{t}__{p}"]) for p in CALLS} for t in CLASSES}

    base = read_matrix(BENCH / "output_v37/concordance_matrix.csv")
    firing = {}
    # Criterion firing counts for the baseline, from its own file.
    with (BENCH / "output_v37/criterion_firing_rates.csv").open() as f:
        for row in csv.DictReader(f):
            firing[row["criterion"]] = sum(int(row[f"{c}_fired"]) for c in CLASSES)

    def line(name, group, knob, value, shipped, m, fire=None):
        mm = matrix_metrics(m)
        return {"name": name, "group": group, "knob": knob, "value": value, "shipped": shipped,
                **{k: mm[k] for k in ("exact", "same_direction", "same_pathogenic", "same_benign", "opposite",
                                      "false_benign", "false_pathogenic", "plp_recall", "blb_recall", "vus_recall",
                                      "plp_to_vus", "blb_to_vus", "vus_to_plp", "vus_to_blb")}, **(fire or {})}

    out = [line("v37", "baseline", "", "shipped", "", base, {f"fire_{c}": v for c, v in firing.items()})]
    for _, row in r.iterrows():
        fire = {c: int(row[c]) for c in r.columns if c.startswith("fire_")}
        out.append(line(row["name"], row["group"], row["knob"], str(row["value"]), str(row["shipped"]), mat_of(row), fire))
    # The dominant-gene BS2 bar: the sweep that was already run on this stack (v37 variants).
    bs2 = [("v37_bs2off", "0"), ("v37_bs2prev_3e-5", "3e-05"), ("v37_bs2prev_1e-4", "0.0001"),
           ("v37_bs2prev_1e-3", "0.001"), ("v37_bs2prev_3e-3", "0.003")]
    for name, val in bs2:
        p = BENCH / f"output_{name}/concordance_matrix.csv"
        if p.exists():
            out.append(line(name, "sweep", "bs2_ad_prevalence_threshold", val, "0.0003", read_matrix(p)))
    for name, val in [("v37_bs2ac_20", "20"), ("v37_bs2ac_100", "100"), ("v37_bs2ac_1000", "1000")]:
        p = BENCH / f"output_{name}/concordance_matrix.csv"
        if p.exists():
            out.append(line(name, "sweep", "bs2_ad_min_ac (bar off)", val, "5 (with bar 3e-4)", read_matrix(p)))
    write("runs_v37_variants.csv", out)

    # The controls of the two measurement builds: with nothing set, each must be the shipped classifier.
    rows = []
    for name in ("v37_cut_lp6", "v37_hk_control"):
        row = r[r["name"] == name]
        if row.empty:
            continue
        got = mat_of(row.iloc[0])
        cells = sum(got[t][pp] == base[t][pp] for t in CLASSES for pp in CALLS)
        fires = [c for c in r.columns if c.startswith("fire_")]
        same_fire = sum(int(row.iloc[0][c]) == firing.get(c.removeprefix("fire_"), 0) for c in fires)
        rows.append({"run": name, "cells_equal_to_v37": cells, "cells": len(CLASSES) * len(CALLS),
                     "criteria_firing_equal_to_v37": same_fire, "criteria": len(fires)})
    write("controls_v37.csv", rows)


def threshold_steps() -> None:
    """What each step along a sweep buys and costs, on the final classifier.

    Within one knob the runs are ordered from the setting that gives the fewest benign-side calls to
    the one that gives the most (for the pathogenic-side PM2 bar, from the fewest pathogenic-side
    calls to the most), and each adjacent pair is compared. The exchange rate is the extra
    correct benign-side (or pathogenic-side) calls bought per extra opposite-direction call of the
    kind that setting trades against (missed diagnoses for the benign side, false-pathogenic calls
    for the pathogenic side). A single call is one variant, and counts of 20 to 60 move by several
    from run to run only if the classifier changes, not by chance: every run is deterministic.
    """
    f = DATA_RUNS()
    if f is None:
        return
    base = f[f.name == "v37"].iloc[0]
    knobs = [
        ("bs1_af_threshold", "same_benign", "false_benign", "BS1 bar", "0.005"),
        ("bs2_hom_prevalence_threshold", "same_benign", "false_benign", "BS2 homozygote bar", "0.001"),
        ("bp7_max_intron_offset", "same_benign", "false_benign", "BP7 far boundary", "300"),
        ("phylop_conserved", "same_benign", "false_benign", "BP7 PhyloP bar", "2.0"),
        ("bp1_max_pathogenic_missense", "same_benign", "false_benign", "BP1 veto count", "3"),
        ("pm2_ad_af_threshold", "same_pathogenic", "false_pathogenic", "PM2 bar, not recessive-only (point system)", "4e-05"),
        ("pm2_ar_af_threshold", "same_pathogenic", "false_pathogenic", "PM2 bar, recessive-only", "7e-05"),
    ]
    rows = []
    for knob, gain, cost, label, shipped in knobs:
        sub = f[(f.knob == knob) & (f.group == "sweep")].copy()
        sub["x"] = sub["value"].astype(float)
        sub = pd.concat([sub, pd.DataFrame([{**base.to_dict(), "value": shipped, "x": float(shipped)}])], ignore_index=True)
        sub = sub.sort_values(gain).reset_index(drop=True)
        for a, b in zip(sub.index[:-1], sub.index[1:]):
            ra, rb = sub.loc[a], sub.loc[b]
            d_gain = int(rb[gain] - ra[gain])
            d_cost = int(rb[cost] - ra[cost])
            rows.append({"knob": label, "from_value": ra["value"], "to_value": rb["value"], "gain_calls": d_gain, "added_errors": d_cost,
                         "calls_per_error": (d_gain / d_cost) if d_cost > 0 else math.inf,
                         "shipped_between": shipped in (ra["value"], rb["value"])})
    write("threshold_steps_v37.csv", rows)


def vcep_bars() -> None:
    """The frequency bars ClinGen Variant Curation Expert Panels have published, from the generated table.

    `analysis/acmg_benchmark/data/vcep_thresholds.toml` is written by `build_vcep_thresholds.py` from
    the CSpec Registry and keeps, for every bar, the sentence it was read from in a sibling audit file.
    Only the gene-level blocks are read here (a `[gene_overrides.GENE]` header and the key = value lines
    under it); the `disorders` sub-blocks repeat a subset of them for one disorder.
    """
    path = ROOT / "analysis/acmg_benchmark/data/vcep_thresholds.toml"
    bars: dict[str, dict[str, float]] = {}
    gene = None
    for line in path.read_text().splitlines():
        m = re.match(r"^\[gene_overrides\.([^.\]]+)\]$", line)
        if m:
            gene = m.group(1)
            bars[gene] = {}
            continue
        if line.startswith("["):
            gene = None
            continue
        m = re.match(r"^(ba1_af_threshold|bs1_af_threshold|pm2_af_threshold)\s*=\s*([0-9.eE+-]+)\s*$", line)
        if m and gene:
            bars[gene][m.group(1)] = float(m.group(2))
    defaults = {"ba1_af_threshold": 0.05, "bs1_af_threshold": 0.005, "pm2_af_threshold": 4e-5}
    rows = []
    for key, default in defaults.items():
        v = np.array(sorted(b[key] for b in bars.values() if key in b))
        rows.append({
            "bar": key.split("_")[0].upper(), "genes_with_a_bar": len(v), "median": float(np.median(v)),
            "q1": float(np.percentile(v, 25)), "q3": float(np.percentile(v, 75)),
            "min": float(v.min()), "max": float(v.max()), "fastvep_default": default,
            "genes_tighter_than_default": int((v < default).sum()),
            "genes_looser_than_default": int((v > default).sum()),
        })
    write("vcep_bars_v37.csv", rows)
    n_bars = sum(len(b) for b in bars.values())
    print(f"  {len(bars)} genes with at least one gene-level block, {n_bars} bars")

    # What the extraction did with every criterion entry it read: emitted, or rejected and why.
    audit = pd.read_csv(ROOT / "analysis/acmg_benchmark/data/vcep_thresholds_audit.tsv", sep="\t", dtype=str, keep_default_na=False)
    counts = audit["status"].value_counts()
    if int(counts.get("ok", 0)) != n_bars:
        sys.exit(f"the audit table has {counts.get('ok', 0)} usable bars and the TOML {n_bars}")
    write("vcep_audit_v37.csv", [{"status": k, "entries": int(v)} for k, v in counts.items()])
    # The same for PM2 alone, whose median bar the text quotes: the 411 entries are BA1, BS1 and PM2 together (137 each).
    pm2_counts = audit[audit["criterion"] == "PM2"]["status"].value_counts()
    write("vcep_audit_pm2_v37.csv", [{"status": k, "entries": int(v)} for k, v in pm2_counts.items()])


def facts(df: pd.DataFrame) -> None:
    """Aggregates the text quotes that no other table holds, one name and one value per row.

    `manuscript/acmg/verify_numbers.py` reads this file and the other tables, so that it never has to
    open the per-variant table, which is not in the repository.
    """
    out: list[tuple[str, float | int | str]] = []
    codes = df["criteria"].map(lambda s_: [c for c in s_.split("&") if c])
    pp = codes.map(lambda c: strength_counts(c)[7])
    sig = df["criteria"].map(signature)

    lp = df[df.call == "Likely_pathogenic"]
    lp9 = lp[pp[lp.index] == 9]
    out += [("lp_calls", len(lp)), ("lp_calls_at_9_points", len(lp9)),
            ("lp_at_9_points_pvs1_pm2_only", int((sig[lp9.index] == "PM2_Supporting+PVS1").sum())),
            ("lp_at_9_points_clinvar_pathogenic", int((lp9.clinvar == "Pathogenic").sum()))]

    bp7 = df[df["criteria"].str.contains(r"(?:^|&)BP7(?:&|$)")]
    called = bp7[bp7.fv3 == "B/LB"]
    out += [("bp7_firings", len(bp7)), ("bp7_firings_clinvar_blb", int((bp7.cv3 == "B/LB").sum())),
            ("bp7_called_blb", len(called)), ("bp7_called_blb_clinvar_blb", int((called.cv3 == "B/LB").sum())),
            ("bp7_called_blb_clinvar_plp", int((called.cv3 == "P/LP").sum())), ("bp7_called_blb_clinvar_vus", int((called.cv3 == "VUS").sum()))]

    g = df.groupby("gene")
    gp = pd.DataFrame({"plp": g.apply(lambda x: int((x.cv3 == "P/LP").sum())),
                       "called": g.apply(lambda x: int(((x.cv3 == "P/LP") & (x.fv3 == "P/LP")).sum()))})
    g10 = gp[gp.plp >= 10]
    out += [("genes_with_10_plp", len(g10)), ("genes_with_10_plp_none_called", int((g10.called == 0).sum())),
            ("genes_with_10_plp_90pct_called", int((g10.called / g10.plp >= 0.9).sum()))]

    trunc = df[(df.cv3 == "P/LP") & (df.fv3 == "VUS") & df.consequence.isin(["stop_gained", "frameshift_variant"])]
    has = trunc[trunc.criteria.str.contains("PVS1")]

    def pvs1_strength(c):
        m = re.search(r"(?:^|&)(PVS1(?:_[A-Za-z]+)?)(?:&|$)", c)
        return m.group(1) if m else ""

    for k, v in has.criteria.map(pvs1_strength).value_counts().items():
        out.append((f"trunc_vus_with_{k.lower()}", int(v)))
    out.append(("trunc_vus", len(trunc)))
    splice = df[(df.cv3 == "P/LP") & (df.fv3 == "VUS") & df.consequence.isin(["splice_donor_variant", "splice_acceptor_variant"])]
    out += [("canonical_splice_vus", len(splice)), ("canonical_splice_vus_with_pvs1", int(splice.criteria.str.contains("PVS1").sum()))]

    # The residual calls: gnomAD counts of the PVS1-on-benign-side group, and the consequence classes of the BP4+BP7 group.
    res = pd.read_csv(OUT / "opposite_direction_v37.csv", dtype=str)
    res["ac"] = pd.to_numeric(res["gnomad_ac"], errors="coerce")
    res["af"] = pd.to_numeric(res["gnomad_af"], errors="coerce")
    res["popmax"] = pd.to_numeric(res["gnomad_popmax_af"], errors="coerce")
    pv = res[res.group == "PVS1 on a ClinVar benign-side variant"]
    common = pv[pv.ac >= 600]
    out += [("records_without_gene_symbol", int((df.gene == "").sum()))]
    # The reviewer's rulings on the 49, by kind of contradiction (our coding; Section 3.6): "missed diagnosis" and "false-pathogenic"
    # name the direction of the error against ClinVar and are not a finding that ClinVar is right.
    rul = pd.read_csv(OUT / "opposite_direction_v37.csv", dtype=str)
    is_fb = rul["clinvar"].isin(["Pathogenic", "Likely_pathogenic"])
    out += [("reviewer_supported_clinvar_false_benign", int(((rul.ruling == "clinvar") & is_fb).sum())),
            ("reviewer_supported_clinvar_false_pathogenic", int(((rul.ruling == "clinvar") & ~is_fb).sum()))]
    # Two quirks of the label script (`01_extract_clinvar_2star.py`) that the class counts of Section 2.2 carry. It maps
    # "Conflicting_classifications_of_pathogenicity" to VUS, which only an expert-panel record reaches the 2-star filter
    # with, and its test for "pathogenic/likely_pathogenic" also matches inside
    # "Likely_pathogenic/Likely_pathogenic,_low_penetrance", so those records are labelled Pathogenic.
    pl = {"Pathogenic", "Likely_pathogenic"}
    conf = df[df.clnsig.str.contains("onflict") & (df.clinvar == "VUS")]
    lowpen = df[df.clnsig.str.startswith("Likely_pathogenic/Likely_pathogenic") & (df.clinvar == "Pathogenic")]
    out += [("labels_conflicting_mapped_to_vus", len(conf)),
            ("labels_conflicting_called_pathogenic_side", int(conf.call.isin(pl).sum())),
            ("labels_lowpen_pair_mapped_to_pathogenic", len(lowpen)),
            ("labels_lowpen_pair_exact_gain_if_relabelled", int((lowpen.call == "Likely_pathogenic").sum() - (lowpen.call == "Pathogenic").sum()))]
    # The validity gate withholds PVS1 from a gene whose ClinGen curations are all weak; a label that holds a semicolon
    # defeats it, because the builder splits phenotypes on ';' and the untagged fragment counts as established.
    cg = read_clingen(BENCH / "sa_sources/clingen_gene_validity.csv")
    weak_cls = {"Limited", "Disputed", "Refuted", "No Known Disease Relationship"}
    weak_only = {g for g, k in cg.groupby("GENE SYMBOL")["CLASSIFICATION"].agg(set).items() if k <= weak_cls}
    semicolon = sorted(set(cg[cg["DISEASE LABEL"].str.contains(";") & cg["GENE SYMBOL"].isin(weak_only)]["GENE SYMBOL"]))
    out += [("clingen_weak_only_genes", len(weak_only)), ("clingen_weak_only_genes_with_semicolon_label", ",".join(semicolon)),
            ("clingen_weak_only_semicolon_genes_in_benchmark", int(df.gene.isin(semicolon).sum()))]
    # How the classifier reads inheritance, which PM2's two bars and BS2's two arms branch on: from the disease names in the
    # ClinGen file, by the substring rules of `OmimData::has_dominant_inheritance` and `has_recessive_inheritance`
    # (`sa_extract.rs`), over the phenotypes of the file the benchmark's database was built from (the genemap parser splits
    # them on ';'). The MOI column is not read.
    def weak_tag(p_: str) -> bool:
        parts = p_.split("(ClinGen ")
        return len(parts) > 1 and parts[1].split("/")[0].strip() in weak_cls

    phen: dict[str, list[str]] = {}
    for line in (BENCH / "sa_sources/clingen_gdv.tsv").read_text(encoding="utf-8").splitlines():
        f_ = line.split("\t")
        if line.startswith("#") or len(f_) < 13 or not f_[8].strip() or f_[12] in ("", "."):
            continue
        phen.setdefault(f_[8].split(",")[0].strip(), []).extend(x.strip() for x in f_[12].split(";"))

    def reads(ps: list[str], kind: str) -> bool:
        for p_ in ps:
            if weak_tag(p_):
                continue
            lo = p_.lower()
            if kind == "ad" and ("autosomal dominant" in lo or "{ad}" in lo or ("dominant" in lo and "recessive" not in lo)):
                return True
            if kind == "ar" and ("autosomal recessive" in lo or "{ar}" in lo or ("recessive" in lo and "dominant" not in lo)):
                return True
        return False

    kinds_ = Counter()
    for ps in phen.values():
        d_, r_ = reads(ps, "ad"), reads(ps, "ar")
        kinds_["both" if d_ and r_ else "dominant_only" if d_ else "recessive_only" if r_ else "neither"] += 1
    out += [("clingen_genes_in_file", len(phen))] + [(f"inheritance_read_{k}", int(kinds_[k])) for k in ("recessive_only", "dominant_only", "both", "neither")]
    # How many benchmark genes the file does not hold at all, which is what a flag for "a gene that is not in the file" would list.
    bench_genes = set(df.gene) - {""}
    trunc_plp = df[(df.cv3 == "P/LP") & df.consequence.isin(["stop_gained", "frameshift_variant"])]
    gap_genes = set(trunc_plp[~trunc_plp.criteria.str.contains("PVS1")].gene)
    out += [("benchmark_genes", len(bench_genes)), ("benchmark_genes_not_in_clingen_file", len(bench_genes - set(phen))),
            ("truncating_plp_genes", trunc_plp.gene.nunique()), ("truncating_plp_genes_not_in_clingen_file", len(set(trunc_plp.gene) - set(phen))),
            ("pvs1_gap_genes", len(gap_genes)), ("pvs1_gap_genes_not_in_clingen_file", len(gap_genes - set(phen)))]
    # The inheritance the file states against the one the classifier reads: of the benchmark genes whose established curations (not
    # Limited, Disputed, Refuted or No Known Disease Relationship) all carry /AD, how many the classifier reads as dominant, and the same
    # for /AR (RUN_VERSIONS.md, "Reading ClinGen's inheritance tag").
    established_moi: dict[str, set[str]] = {}
    for g_, ps in phen.items():
        for p_ in ps:
            m_ = re.search(r"\(ClinGen ([^/]+)/([^,)]+)", p_)
            if m_ and m_.group(1).strip() not in weak_cls:
                established_moi.setdefault(g_, set()).add(m_.group(2).strip())
    ad_only = [g_ for g_ in bench_genes if established_moi.get(g_) == {"AD"}]
    ar_only = [g_ for g_ in bench_genes if established_moi.get(g_) == {"AR"}]
    out += [("benchmark_genes_established_tags_ad_only", len(ad_only)), ("benchmark_genes_established_tags_ad_only_read_dominant", sum(reads(phen[g_], "ad") for g_ in ad_only)),
            ("benchmark_genes_established_tags_ar_only", len(ar_only)), ("benchmark_genes_established_tags_ar_only_read_recessive", sum(reads(phen[g_], "ar") for g_ in ar_only))]
    # The strict reading of Section 3.4 (v37_hk_strict) and what its added contradictions are called: the three adjustments it removes
    # are 15 curated alleles, ClinVar's low-penetrance gate and the veto's own-record count, and a Benign call among them contradicts
    # ClinVar's pathogenic-side class at the tier that the abstract says has none.
    strict = pd.read_csv(BENCH / "manuscript_runs/v37_hk_strict/discrepancies.tsv", sep="\t", dtype=str, keep_default_na=False)
    s_opp = strict[(strict.clinvar.isin(PATH) & strict.predicted.isin(BEN)) | (strict.clinvar.isin(BEN) & strict.predicted.isin(PATH))]
    out += [("strict_opposite_calls", len(s_opp)), ("strict_opposite_called_benign", int((s_opp.predicted == "Benign").sum())),
            ("strict_opposite_called_pathogenic", int((s_opp.predicted == "Pathogenic").sum())),
            ("strict_opposite_called_benign_genes", ",".join(sorted(s_opp[s_opp.predicted == "Benign"].gene)))]
    # The curated-mechanism route to loss of function: how many genes the shipped table holds, and whether the ones that
    # name loss of function would be licensed anyway, by gnomAD constraint or by a ClinGen curation that is not weak.
    mech = shipped_mechanisms()
    lof_genes = sorted(g for g, m in mech.items() if "LOF" in m.upper())
    established_cg = {g for g, k in cg.groupby("GENE SYMBOL")["CLASSIFICATION"].agg(set).items() if k - weak_cls}
    pli_ = pd.to_numeric(df.pli, errors="coerce"); loeuf_ = pd.to_numeric(df.loeuf, errors="coerce")
    constrained_genes = set(df[(pli_ >= 0.9) | (loeuf_ <= 0.35)].gene)
    out += [("mechanism_table_entries", len(mech)),
            ("mechanism_table_gof", sum(m.upper() == "GOF" for m in mech.values())),
            ("mechanism_table_dominant_negative", sum(m.upper() == "DOMINANT_NEGATIVE" for m in mech.values())),
            ("mechanism_table_lof_and_gof", len(lof_genes)),
            ("mechanism_lof_genes_licensed_anyway", sum(g in constrained_genes or g in established_cg for g in lof_genes))]
    out += [("pvs1_benign_lone_pvs1", int((pv.criteria == "PVS1").sum())), ("pvs1_benign_pvs1_pm2", int((pv.criteria == "PVS1&PM2_Supporting").sum())),
            ("pvs1_benign_common_n", len(common)), ("pvs1_benign_common_ac_min", int(common.ac.min())), ("pvs1_benign_common_ac_max", int(common.ac.max())),
            ("pvs1_benign_common_af_min", float(common.af.min())), ("pvs1_benign_common_af_max", float(common.af.max())),
            ("pvs1_benign_common_genes", ",".join(sorted(common.gene))),
            ("pvs1_benign_common_popmax_min", float(common.popmax.min())), ("pvs1_benign_common_popmax_max", float(common.popmax.max())),
            ("pvs1_benign_common_bs1_fired", int(common.criteria.str.contains("BS1").sum()))]
    bb = res[res.group == "BP4+BP7 on a ClinVar pathogenic-side variant"]
    for k, v in bb.consequence.value_counts().items():
        out.append((f"bp4bp7_{k}", int(v)))
    # PVS1 graded as the published tree grades it. The implementation grades an NMD-escaping truncation
    # outside a critical region one step below the tree (Moderate and Supporting where the tree gives
    # Strong and Moderate), and a canonical splice variant in the last exon as Moderate whatever the
    # fraction removed. The recorded code does not say which branch gave PVS1_Moderate, so this takes
    # every PVS1_Moderate on a stop, frameshift or canonical splice variant up to Strong and every
    # PVS1_Supporting on a stop or frameshift up to Moderate: an upper bound on what a regrade changes.
    changed = Counter()
    for c, cons in zip(codes, df["consequence"]):
        if not any(x in ("PVS1_Moderate", "PVS1_Supporting") for x in c):
            continue
        new = []
        for x in c:
            if x == "PVS1_Moderate" and cons in ("stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant"):
                new.append("PVS1_Strong")
            elif x == "PVS1_Supporting" and cons in ("stop_gained", "frameshift_variant"):
                new.append("PVS1_Moderate")
            else:
                new.append(x)
        a, b = replay(c), replay(new)
        if a != b:
            changed[(a, b)] += 1
    out.append(("pvs1_regrade_calls_changed", sum(changed.values())))
    out.append(("pvs1_regrade_vus_to_lp", changed[("VUS", "LP")]))
    out.append(("pvs1_regrade_b_to_vus", changed[("B", "VUS")]))
    out.append(("pvs1_regrade_other", sum(changed.values()) - changed[("VUS", "LP")] - changed[("B", "VUS")]))
    write("facts_v37.csv", pd.DataFrame(out, columns=["name", "value"]))


def DATA_RUNS():
    p = OUT / "runs_v37_variants.csv"
    return pd.read_csv(p) if p.exists() else None


def review_rounds() -> None:
    """Rows and verdicts of every workbook the reviewer returned, decoded from the files.

    Rounds 1 and 2 recorded a verdict as a cell fill, decoded here against each workbook's own key.
    Round 6 recorded it as text, coded by us (`round6_rulings.tsv`); the reviewer did not assign
    the categories. All three are mapped onto one set so they can be read together.
    """
    import os

    import openpyxl

    def fill_of(cell):
        if cell.fill and cell.fill.fill_type:
            fg = cell.fill.fgColor
            return fg.rgb if fg.type == "rgb" else ("theme", fg.theme, round(fg.tint or 0, 2))
        return None

    rows = []
    # Round 1: the v7 queue. Key: theme 6 tint 0.8 = correct call, theme 8 tint 0.4 = both discordant.
    p = Path(os.path.expanduser("~/Downloads/FastVEP-ClinVar-discordance-MG.xlsx"))
    if p.exists():
        ws = openpyxl.load_workbook(p, data_only=True)["Sheet1"]
        hdr = [c.value for c in ws[1]]
        ix = {h: i for i, h in enumerate(hdr)}
        tally, n, reason = Counter(), 0, 0
        for r in ws.iter_rows(min_row=2):
            if r[ix["chrom"]].value is None:
                continue
            n += 1
            reason += bool(r[ix["REASON FOR DISCORDANCE"]].value)
            t, f = fill_of(r[ix["truth_class"]]), fill_of(r[ix["fastvep_class"]])
            if t == ("theme", 6, 0.8):
                tally["clinvar_supported"] += 1
            elif t == ("theme", 8, 0.4):
                tally["neither_supported"] += 1
            elif f in ("FF92D050", ("theme", 6, 0.4)):
                tally["fastvep_supported"] += 1
            else:
                tally["no_verdict"] += 1
        rows.append({"round": 1, "run": "v7", "rows": n, "rows_with_written_reason": reason, **tally})
    # Round 2: the v8 queue. Key: green = correct call, theme 8 = both discordant, theme 5 = needs literature.
    p = Path(os.path.expanduser("~/Downloads/discordant 122.xlsx"))
    if p.exists():
        ws = openpyxl.load_workbook(p, data_only=True)["Sheet1"]
        hdr = [c.value for c in ws[1]]
        ix = {h: i for i, h in enumerate(hdr)}
        tally, n = Counter(), 0
        homology, homology_on_call = 0, 0
        for r in ws.iter_rows(min_row=2):
            if r[ix["chrom"]].value is None:
                continue
            n += 1
            t, f = fill_of(r[ix["truth_class"]]), fill_of(r[ix["fastvep_class"]])
            # The orange marker for "a gene with a homology problem" is on the identifying cells of the row.
            if fill_of(r[ix["gene"]]) == "FFFFC000":
                homology += 1
                homology_on_call += f == "FFFFC000"
            if t == "FF92D050":
                tally["clinvar_supported"] += 1
            elif f == "FF92D050":
                tally["fastvep_supported"] += 1
            elif t == ("theme", 8, 0.4):
                tally["neither_supported"] += 1
            elif t == ("theme", 5, 0.4):
                tally["needs_evidence"] += 1
            else:
                tally["no_verdict"] += 1
        rows.append({"round": 2, "run": "v8", "rows": n, **tally, "homology_marked_rows": homology,
                     "homology_marked_on_fastvep_call": homology_on_call})
    # Round 6: the v26 queue, coded by us from the returned text.
    p = OUT / "round6_rulings.tsv"
    if p.exists():
        rv = pd.read_csv(p, sep="\t", dtype=str)
        c = Counter(rv["ruling"])
        notes = ""
        wb6 = Path(os.path.expanduser("~/Downloads/discordant_64_round6-1.xlsx"))
        if wb6.exists():
            ws6 = openpyxl.load_workbook(wb6, data_only=True)["review"]
            h6 = [c_.value for c_ in ws6[1]]
            ix6 = {h: i for i, h in enumerate(h6) if h}
            # The two columns the workbook asked the reviewer to fill came back empty; the verdict column holds the text.
            notes = sum(bool(r[ix6["md_note"]].value) for r in ws6.iter_rows(min_row=2) if r[ix6["gene"]].value is not None)
            empty_cols = sum(not any(r[ix6[k]].value for r in ws6.iter_rows(min_row=2) if r[ix6["gene"]].value is not None) for k in ("MD_CALL", "WHO_IS_RIGHT"))
            if empty_cols != 2:
                sys.exit("the class and who-is-right columns of round 6 are not empty")
        rows.append({"round": 6, "run": "v26", "rows": len(rv), "rows_with_written_reason": notes,
                     "clinvar_supported": c["clinvar"], "fastvep_supported": c["fastvep"],
                     "neither_supported": c["neither"],
                     "needs_evidence": c["needs_evidence"] + c["disorder_specific"]})
    df = pd.DataFrame(rows)
    cats = ["clinvar_supported", "fastvep_supported", "neither_supported", "needs_evidence", "no_verdict"]
    for col in cats:
        if col not in df:
            df[col] = 0
    df[cats] = df[cats].fillna(0).astype(int)
    assert (df[cats].sum(axis=1) == df["rows"]).all(), df
    # Shares are taken over the rows that carry a recorded verdict, the same denominator in every round.
    df["rows_with_verdict"] = df["rows"] - df["no_verdict"]
    df["clinvar_share_of_verdicts"] = df["clinvar_supported"] / df["rows_with_verdict"]
    write("review_rounds.csv", df)


def holdout_rates(apr: pd.DataFrame, octo: pd.DataFrame) -> None:
    """Is the higher rate of opposite-direction calls on new variants a difference in kind or in mix?

    The exact binomial intervals of `temporal_opposite_rates_v37.csv` treat variants as independent, and
    the new variants lie in the genes of the April set. Three analyses that do not:

      ratio      the rate on the new variants over the rate on the April set, with a gene-cluster
                 bootstrap that resamples genes jointly across the two sets (2,000 resamples, 95%
                 percentile, the seed of `headline`)
      expected   the number of opposite-direction calls the new variants would have at the April rate of
                 their own stratum (ClinVar direction by consequence group, ClinVar direction alone, and the same
                 two with ClinVar's five classes in place of its three directions), with the Poisson probability
                 of at least the observed number
      recall     the recall in direction of each ClinVar class on the new variants, against what the April
                 recall of each population-maximum allele-frequency band would give them at their own mix
                 of bands
    """
    new = octo[octo["stratum"] == "new"].copy()
    apr = apr.copy()
    for d in (apr, new):
        d["opp"] = (((d.cv3 == "P/LP") & (d.fv3 == "B/LB")) | ((d.cv3 == "B/LB") & (d.fv3 == "P/LP"))).astype(int)
        d["cgrp"] = d["consequence"].map(consequence_group)

    # ratio, by gene cluster
    ga = apr.groupby("gene").agg(n_apr=("opp", "size"), o_apr=("opp", "sum"))
    gn = new.groupby("gene").agg(n_new=("opp", "size"), o_new=("opp", "sum"))
    g = ga.join(gn, how="outer").fillna(0)
    arr = g[["n_apr", "o_apr", "n_new", "o_new"]].to_numpy(dtype=float)
    rng = np.random.default_rng(20261011)
    idx = rng.integers(0, len(arr), size=(2000, len(arr)))
    boot = np.stack([arr[i].sum(axis=0) for i in idx])
    r_apr, r_new = boot[:, 1] / boot[:, 0], boot[:, 3] / boot[:, 2]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = r_new / r_apr
    ratio = ratio[np.isfinite(ratio)]
    tot = arr.sum(axis=0)
    rows = [{"measure": "opposite-direction rate, April set", "value": tot[1] / tot[0],
             "ci_lo": np.percentile(r_apr, 2.5), "ci_hi": np.percentile(r_apr, 97.5), "n": int(tot[0]), "k": int(tot[1])},
            {"measure": "opposite-direction rate, new variants", "value": tot[3] / tot[2],
             "ci_lo": np.percentile(r_new, 2.5), "ci_hi": np.percentile(r_new, 97.5), "n": int(tot[2]), "k": int(tot[3])},
            {"measure": "ratio of the new rate to the April rate", "value": (tot[3] / tot[2]) / (tot[1] / tot[0]),
             "ci_lo": np.percentile(ratio, 2.5), "ci_hi": np.percentile(ratio, 97.5), "n": len(g), "k": ""}]
    # expected number at the April rate of the same stratum
    # The first two stratifications are the ones the text leads with; the two that follow use ClinVar's five classes, which
    # `expected` above names in its docstring as "ClinVar direction" only, so that the wording "class" is tested as well.
    for label, keys in (("ClinVar direction", ["cv3"]), ("ClinVar direction by consequence group", ["cv3", "cgrp"]),
                        ("ClinVar class", ["clinvar"]), ("ClinVar class by consequence group", ["clinvar", "cgrp"])):
        ra = apr.groupby(keys)["opp"].agg(["sum", "size"])
        rate = (ra["sum"] / ra["size"]).rename("rate")
        nn = new.groupby(keys).size().rename("n_new")
        j = pd.concat([nn, rate], axis=1).dropna()
        unseen = int(nn.sum() - j["n_new"].sum())
        expected = float((j["n_new"] * j["rate"]).sum())
        k = int(new["opp"].sum())
        p = float(stats.poisson.sf(k - 1, expected))
        rows.append({"measure": f"expected opposite-direction calls among the new variants at the April rate by {label}",
                     "value": expected, "ci_lo": "", "ci_hi": "", "n": int(j["n_new"].sum()), "k": k,
                     "poisson_p_at_least_observed": p, "new_variants_in_strata_without_april_rate": unseen})
    write("temporal_ratio_v37.csv", rows)

    # recall in direction, standardised for the allele-frequency mix
    # Nine bands. Six (the last two joined, and the three from 5e-3 up as one) left a fall of 5.8 points unexplained in the
    # commonest band, which a finer split of that band explains: the new Benign variants sit at its low end.
    def band(a):
        if a == "":
            return "absent"
        x = float(a)
        for hi, name in ((1e-5, "<1e-5"), (1e-4, "1e-5 to 1e-4"), (1e-3, "1e-4 to 1e-3"), (5e-3, "1e-3 to 5e-3"), (1e-2, "5e-3 to 1e-2"),
                         (2e-2, "1e-2 to 2e-2"), (5e-2, "2e-2 to 5e-2")):
            if x < hi:
                return name
        return ">=5e-2"

    order = ["absent", "<1e-5", "1e-5 to 1e-4", "1e-4 to 1e-3", "1e-3 to 5e-3", "5e-3 to 1e-2", "1e-2 to 2e-2", "2e-2 to 5e-2", ">=5e-2"]
    for d in (apr, new):
        d["afband"] = d["gnomad_popmax_af"].map(band)
        d["same"] = (d.cv3 == d.fv3).astype(int)
    rows = []
    for cls in CLASSES[:2] + CLASSES[3:]:
        a, n_ = apr[apr.clinvar == cls], new[new.clinvar == cls]
        rec_a = a.groupby("afband")["same"].mean()
        w = n_.groupby("afband").size() / len(n_)
        exp = float((w * rec_a.reindex(w.index)).sum())
        # The six bands tried first, before the finer split above 1e-2 was made: the two lowest above 1e-4 joined, and the three
        # from 1e-2 up as one.
        coarse = {"1e-3 to 5e-3": "1e-3 to 1e-2", "5e-3 to 1e-2": "1e-3 to 1e-2", "1e-2 to 2e-2": ">=1e-2", "2e-2 to 5e-2": ">=1e-2", ">=5e-2": ">=1e-2"}
        rec6 = a.groupby(a["afband"].replace(coarse))["same"].mean()
        w6 = n_.groupby(n_["afband"].replace(coarse)).size() / len(n_)
        exp6 = float((w6 * rec6.reindex(w6.index)).sum())
        rows.append({"clinvar": cls, "april_recall": a["same"].mean(), "new_recall": n_["same"].mean(),
                     "new_recall_at_april_recall_per_band": exp, "n_april": len(a), "n_new": len(n_),
                     "share_of_the_difference_explained": (a["same"].mean() - exp) / (a["same"].mean() - n_["same"].mean())
                     if a["same"].mean() != n_["same"].mean() else math.nan,
                     "six_band_new_recall_at_april_recall_per_band": exp6,
                     "six_band_share_of_the_difference_explained": (a["same"].mean() - exp6) / (a["same"].mean() - n_["same"].mean())
                     if a["same"].mean() != n_["same"].mean() else math.nan})
    write("temporal_standardised_recall_v37.csv", rows)
    rows = []
    for cls in ("Benign", "Likely_benign"):
        for b in order:
            a = apr[(apr.clinvar == cls) & (apr.afband == b)]
            n_ = new[(new.clinvar == cls) & (new.afband == b)]
            rows.append({"clinvar": cls, "popmax_af_band": b, "n_april": len(a), "recall_april": a["same"].mean() if len(a) else math.nan,
                         "n_new": len(n_), "recall_new": n_["same"].mean() if len(n_) else math.nan})
    write("temporal_recall_by_af_band_v37.csv", rows)

    # The highest bands, where BA1 and BS1 should decide: the new Benign variants at 2% or more are called benign-side far less
    # often than the April ones, and the criteria of the ones that are not show what did not fire.
    rows = []
    for name, d in (("April set", apr), ("new variants", new)):
        hi = d[(d.clinvar == "Benign") & (pd.to_numeric(d["gnomad_popmax_af"].replace("", np.nan), errors="coerce") >= 2e-2)]
        nb = hi[~hi.call.isin(["Likely_benign", "Benign"])]
        toks = nb["criteria"].map(lambda c: [t.split("_")[0] for t in c.split("&") if t])
        rows.append({"set": name, "benign_variants_at_2e-2_or_more": len(hi), "called_benign_side": len(hi) - len(nb),
                     "share_called_benign_side": (len(hi) - len(nb)) / len(hi), "not_called_benign_side": len(nb),
                     "of_them_with_bp4_alone": int(toks.map(lambda t: t == ["BP4"]).sum()),
                     "of_them_with_ba1_or_bs1": int(toks.map(lambda t: "BA1" in t or "BS1" in t).sum())})
    write("temporal_high_af_benign_v37.csv", rows)

    # Agreement in direction on the new variants against April's, with the mix of ClinVar directions held at the April mix of
    # each set, so that the part of the fall that is a smaller share of VUS (97% concordant) is separate from the recalls.
    rec_dir = apr.groupby("cv3")["same"].mean()
    rows = []
    for name, d in (("April set", apr), ("new variants", new)):
        mix = d.groupby("cv3").size() / len(d)
        own_rec = d.groupby("cv3")["same"].mean()
        mix_a = apr.groupby("cv3").size() / len(apr)
        rows.append({"set": name, "variants": len(d), "share_p_lp": mix.get("P/LP", 0), "share_vus": mix.get("VUS", 0),
                     "share_b_lb": mix.get("B/LB", 0), "agreement_in_direction": d["same"].mean(),
                     "agreement_at_april_recall_of_each_direction": float((mix * rec_dir.reindex(mix.index)).sum()),
                     # the other order: the set's own recalls at the April mix, which splits the fall into a recall part first
                     "agreement_at_own_recall_and_april_mix": float((mix_a * own_rec.reindex(mix_a.index)).sum())})
    write("temporal_agreement_mix_v37.csv", rows)


def temporal() -> None:
    """The final classifier against ClinVar 2026-10-04: the whole covered set and each stratum of it.

    Strata are defined in `14_temporal_holdout.py`: `new` variants were not in the April benchmark set
    at all, so no iteration, review round or sweep ever saw them.
    """
    pv = BENCH / "manuscript_runs/v37_oct2026_per_variant.tsv.gz"
    st = BENCH / "clinvar_20261004_covered/strata.tsv"
    if not (pv.exists() and st.exists()):
        print("  (no October run yet: skipping)")
        return
    df = pd.read_csv(pv, sep="\t", dtype=str, keep_default_na=False)
    df["cv3"] = df["clinvar"].map(TRI)
    df["fv3"] = df["call"].map(TRI)
    df["stars"] = df["stars"].astype(int)
    strata = pd.read_csv(st, sep="\t", dtype=str)
    df = df.merge(strata, on=["chrom", "pos", "ref", "alt"], how="left")
    if df["stratum"].isna().any():
        sys.exit("a scored October variant has no stratum")
    april = pd.read_csv(BENCH / "clinvar_2star_labels.tsv", sep="\t", dtype=str)
    april = april.rename(columns={"normalized_class": "april_class"})[["chrom", "pos", "ref", "alt", "april_class"]]
    df = df.merge(april, on=["chrom", "pos", "ref", "alt"], how="left")

    def block(sub: pd.DataFrame, label: str) -> list[dict]:
        out = []
        n = len(sub)
        exact = int((sub.clinvar == sub.call).sum())
        same = int((sub.cv3 == sub.fv3).sum())
        fb = int(((sub.cv3 == "P/LP") & (sub.fv3 == "B/LB")).sum())
        fp = int(((sub.cv3 == "B/LB") & (sub.fv3 == "P/LP")).sum())
        for a in ("P/LP", "VUS", "B/LB"):
            sa = sub[sub.cv3 == a]
            for b in ("P/LP", "VUS", "B/LB", "NoCall"):
                k = int((sa.fv3 == b).sum())
                lo, hi = clopper_pearson(k, len(sa)) if len(sa) else (math.nan, math.nan)
                out.append({"stratum": label, "clinvar": a, "fastvep": b, "n": k, "row_total": len(sa),
                            "row_fraction": k / len(sa) if len(sa) else math.nan, "ci_lo": lo, "ci_hi": hi,
                            "stratum_n": n, "exact": exact, "same_direction": same,
                            "opposite": fb + fp, "false_benign": fb, "false_pathogenic": fp})
        return out

    rows = block(df, "all covered")
    for name in ("unchanged", "tier_change", "direction_change", "new"):
        rows += block(df[df.stratum == name], name)
    # By review status within the new variants (the 3-star ones are expert-panel calls).
    new = df[df.stratum == "new"]
    for s_ in (2, 3):
        rows += block(new[new.stars == s_], f"new, {s_}-star")
    write("temporal_v37.csv", rows)

    # The 157 variants whose class moved across a direction: what fastVEP said, against the April and the October class.
    dc = df[df.stratum == "direction_change"].copy()
    dc["april3"] = dc["april_class"].map(TRI)
    t = (dc.groupby(["april3", "cv3", "fv3"]).size().reset_index(name="n"))
    write("temporal_direction_changes_v37.csv", t)
    rows = []
    for _, r in dc.iterrows():
        rows.append({"gene": r["gene"], "chrom": r["chrom"], "pos": r["pos"], "ref": r["ref"], "alt": r["alt"],
                     "april_class": r["april_class"], "october_class": r["clinvar"], "fastvep": r["call"],
                     "criteria": r["criteria"], "consequence": r["consequence"]})
    write("temporal_direction_change_variants_v37.csv", rows)
    # The composition of the new variants: gene concentration and consequence classes.
    new = df[df.stratum == "new"].copy()
    new["cgrp"] = new["consequence"].map(consequence_group)
    rows = [{"consequence_group": g, "clinvar": a, "n": int(len(x))} for (g, a), x in new.groupby(["cgrp", "cv3"])]
    write("temporal_new_composition_v37.csv", rows)
    # Class-wise same-direction recall, April set against the October strata, with the allele-frequency
    # context that explains the one class that moved (novel Benign variants are rarer than the old ones).
    apr = load()  # every record, NoCall included, as in the October strata: a reference-agreement record is not a concordant call
    groups = [("April set (v37)", apr), ("October, unchanged", df[df.stratum == "unchanged"]),
              ("October, new", df[df.stratum == "new"])]
    rows = []
    for name, sub in groups:
        for cls in CLASSES:
            sc = sub[sub.clinvar == cls]
            k = int((sc.cv3 == sc.fv3).sum())
            lo, hi = clopper_pearson(k, len(sc)) if len(sc) else (math.nan, math.nan)
            af = pd.to_numeric(sc["gnomad_popmax_af"], errors="coerce")
            rows.append({"group": name, "clinvar": cls, "n": len(sc), "same_direction": k,
                         "recall": k / len(sc) if len(sc) else math.nan, "ci_lo": lo, "ci_hi": hi,
                         "has_gnomad_popmax": int(af.notna().sum()),
                         "median_popmax_af": float(af.median()) if af.notna().any() else math.nan,
                         "share_popmax_above_0.5pct": float((af > 0.005).sum() / len(sc)) if len(sc) else math.nan})
    write("temporal_by_class_v37.csv", rows)
    top = new["gene"].value_counts().head(15)
    write("temporal_new_top_genes_v37.csv", [{"gene": g, "n": int(v)} for g, v in top.items()])

    # The opposite-direction calls among the new variants, one row each.
    newopp = new[((new.cv3 == "P/LP") & (new.fv3 == "B/LB")) | ((new.cv3 == "B/LB") & (new.fv3 == "P/LP"))].copy()

    def kind(r):
        sig_ = signature(r["criteria"])
        if "PVS1" in sig_ and r["cv3"] == "B/LB":
            return "PVS1 on a ClinVar benign-side variant"
        if "BP7" in sig_ and r["cv3"] == "P/LP":
            return "BP4+BP7 on a ClinVar pathogenic-side variant"
        if r["cv3"] == "B/LB" and r["consequence"] == "missense_variant":
            return "missense pathogenic-side call on PP2 and PP3 for a ClinVar benign-side variant"
        return "other"

    newopp["kind"] = newopp.apply(kind, axis=1)
    write("temporal_new_opposite_v37.csv", newopp[["kind", "gene", "chrom", "pos", "ref", "alt", "stars", "clinvar", "call", "criteria", "consequence"]].sort_values(["kind", "gene"]))

    # How many October records the gnomAD extract covers, by stratum, and how many the rule drops.
    import importlib.util

    spec = importlib.util.spec_from_file_location("temporal_holdout", HERE / "14_temporal_holdout.py")
    th = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(th)
    regions = th.load_regions()
    apr_cls = {(r["chrom"], r["pos"], r["ref"], r["alt"]): r["normalized_class"]
               for r in csv.DictReader((BENCH / "clinvar_2star_labels.tsv").open(), delimiter="\t")}
    cov = Counter()
    with (BENCH / "clinvar_20261004/clinvar_2star_labels.tsv").open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            key = (r["chrom"], r["pos"], r["ref"], r["alt"])
            st_ = th.stratum(apr_cls, key, r["normalized_class"])
            cov[(st_, "covered" if th.covered(regions, r["chrom"], int(r["pos"])) else "outside the extract")] += 1
    rows = []
    for name in ("unchanged", "tier_change", "direction_change", "new"):
        k, d = cov[(name, "covered")], cov[(name, "outside the extract")]
        rows.append({"stratum": name, "records": k + d, "covered": k, "outside_the_extract": d})
    rows.append({"stratum": "all", "records": sum(cov.values()), "covered": sum(v for (a, b), v in cov.items() if b == "covered"),
                 "outside_the_extract": sum(v for (a, b), v in cov.items() if b != "covered")})
    write("temporal_coverage_v37.csv", rows)
    if rows[-1]["covered"] != len(df):
        sys.exit("the covered count does not match the October run")

    # The rate of opposite-direction calls, with an exact interval, for the April set and each October stratum.
    rows = []
    for name, sub in (("April set", load()), ("October, unchanged", df[df.stratum == "unchanged"]), ("October, new", df[df.stratum == "new"]),
                      ("October, all covered", df)):
        k = int((((sub.cv3 == "P/LP") & (sub.fv3 == "B/LB")) | ((sub.cv3 == "B/LB") & (sub.fv3 == "P/LP"))).sum())
        lo, hi = clopper_pearson(k, len(sub))
        rows.append({"set": name, "opposite": k, "records": len(sub), "rate": k / len(sub), "ci_lo": lo, "ci_hi": hi})
    write("temporal_opposite_rates_v37.csv", rows)
    holdout_rates(load(), df)

    # What became of the 49 residual opposite-direction calls of the April set in the later release.
    # A call that is still in the 2-star+ set keeps its class (or moves); one that left it was, in
    # the raw release, either downgraded in review status or reclassified to Conflicting.
    res = pd.read_csv(OUT / "opposite_direction_v37.csv", dtype=str)
    in_oct = df.set_index(["chrom", "pos", "ref", "alt"])
    missing = {}
    rows = []
    for _, r in res.iterrows():
        key = (r["chrom"], r["pos"], r["ref"], r["alt"])
        if key in in_oct.index:
            o = in_oct.loc[key]
            rows.append({"gene": r["gene"], "chrom": r["chrom"], "pos": r["pos"], "ref": r["ref"], "alt": r["alt"],
                         "april_clinvar": r["clinvar"], "fastvep_call": r["call"], "october_status": "still 2-star+",
                         "october_clinvar": o["clinvar"], "october_stratum": o["stratum"], "october_review_status": ""})
        else:
            missing[(r["chrom"].removeprefix("chr"), r["pos"], r["ref"], r["alt"])] = r
    if missing:
        import gzip
        with gzip.open(BENCH / "sa_sources/clinvar_20261004.vcf.gz", "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.split("\t", 8)
                key = (c[0].removeprefix("chr"), c[1], c[3], c[4])
                if key in missing:
                    info = dict(kv.split("=", 1) for kv in c[7].split(";") if "=" in kv)
                    r = missing.pop(key)
                    rows.append({"gene": r["gene"], "chrom": r["chrom"], "pos": r["pos"], "ref": r["ref"], "alt": r["alt"],
                                 "april_clinvar": r["clinvar"], "fastvep_call": r["call"], "october_status": "left the 2-star+ set",
                                 "october_clinvar": info.get("CLNSIG", ""), "october_stratum": "", "october_review_status": info.get("CLNREVSTAT", "")})
                    if not missing:
                        break
    if missing:
        sys.exit(f"{len(missing)} of the April opposite-direction calls are neither in the October covered set nor in the October release")
    write("april_opposite_in_october_v37.csv", rows)
    april_genes = set(pd.read_csv(PER_VARIANT, sep="\t", usecols=["gene"], dtype=str)["gene"])
    new_genes = set(new["gene"])
    print(f"  new variants: {len(new):,} in {len(new_genes):,} genes, {len(new_genes - april_genes):,} of them not in the April set")
    # The protein index that PS1, PM1, PM5 and the BP1 veto read was not rebuilt, so it lacks the records added since April.
    # How many of the opposite-direction calls on the new variants rest on one of those criteria bounds what that can do.
    uses = lambda code: int(newopp["criteria"].map(lambda c: any(t.split("_")[0] == code for t in c.split("&"))).sum())  # noqa: E731
    write("temporal_new_genes_v37.csv", [{"new_variants": len(new), "genes": len(new_genes), "genes_not_in_the_april_set": len(new_genes - april_genes),
                                          "opposite_calls": len(newopp), "opposite_calls_with_bp1": uses("BP1"), "opposite_calls_with_ps1": uses("PS1"),
                                          "opposite_calls_with_pm5": uses("PM5"), "opposite_calls_with_pm1": uses("PM1")}])


def new_in_index() -> None:
    """How the new variants of the October release stand in the ClinVar inputs that the April classifier reads.

    The protein index is built from every classified record of the late-April `variant_summary`, whatever its review status,
    so a variant that reached two stars after the benchmark set was drawn is mostly in it, with its April record. Counted
    here: the April record of each new variant (matched on GRCh38 VCF coordinates), and what that leaves of the test that
    the held-out run of BP1's veto without the variant's own record makes. The veto counts the gene's pathogenic missense
    records in the index, the variant's own among them, so the run can only change a call where the gene's count is
    exactly at the bar and the variant's own record is in it.
    """
    pv = BENCH / "manuscript_runs/v37_oct2026_per_variant.tsv.gz"
    st = BENCH / "clinvar_20261004_covered/strata.tsv"
    if not (pv.exists() and st.exists()):
        print("  (no October run yet: skipping)")
        return
    key = ["chrom", "pos", "ref", "alt"]
    new = pd.read_csv(pv, sep="\t", dtype=str, keep_default_na=False).merge(pd.read_csv(st, sep="\t", dtype=str), on=key, how="left")
    new = new[new.stratum == "new"].copy()
    summary = pd.read_csv(BENCH / "sa_sources/variant_summary.txt.gz", sep="\t", dtype=str, low_memory=False,
                          usecols=["Assembly", "Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF", "ReviewStatus", "ClinicalSignificance"])
    summary = summary[(summary.Assembly == "GRCh38") & (summary.PositionVCF != "-1")]
    summary["k"] = summary.Chromosome + ":" + summary.PositionVCF + ":" + summary.ReferenceAlleleVCF + ":" + summary.AlternateAlleleVCF
    summary = summary.drop_duplicates("k").set_index("k")
    new["k"] = new.chrom + ":" + new.pos + ":" + new.ref + ":" + new.alt
    rec = new.join(summary[["ReviewStatus", "ClinicalSignificance"]], on="k")
    present = rec.ReviewStatus.notna()
    status = rec.ReviewStatus.fillna("")
    sig = rec.ClinicalSignificance.fillna("").str.lower()
    plp = sig.isin(["pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"])
    single = status == "criteria provided, single submitter"
    multiple = status == "criteria provided, multiple submitters, no conflicts"
    expert = status.isin(["reviewed by expert panel", "practice guideline"])
    rows = [("new variants of the covered release", len(new)),
            ("  with a record in the April variant_summary", int(present.sum())),
            ("    a single-submitter record", int(single.sum())),
            ("    a record with multiple submitters and no conflict", int(multiple.sum())),
            ("    a record reviewed by an expert panel", int(expert.sum())),
            ("    a record of another review status", int((present & ~single & ~multiple & ~expert).sum())),
            ("  without a record in the April variant_summary", int((~present).sum())),
            ("  with a pathogenic or likely pathogenic record in the April variant_summary", int(plp.sum()))]
    # The test of the held-out run of BP1's veto without the variant's own record.
    index = read_protein_index(BENCH / "sa_db/clinvar_protein.oga")
    can = new[(new.clinvar.isin(PATH)) & (new.consequence == "missense_variant")
              & (pd.to_numeric(new.pli, errors="coerce") >= 0.9) & (pd.to_numeric(new.mis_z, errors="coerce") < 2)]
    in_index, at_bar, at_bar_bp4 = 0, 0, 0
    for r in can.itertuples():
        row = _index_row(index, r.gene, r.hgvsp, "")
        in_index += row["own_record_among_them"]
        if row["own_record_among_them"] and row["pathogenic_missense_records"] == 3:
            at_bar += 1
            at_bar_bp4 += "BP4" in r.criteria
    # the same count for the April set, where the index was built from the records of the variants themselves
    april = pd.read_csv(PER_VARIANT, sep="\t", dtype=str, keep_default_na=False)
    a_can = april[(april.clinvar.isin(PATH)) & (april.consequence == "missense_variant")
                  & (pd.to_numeric(april.pli, errors="coerce") >= 0.9) & (pd.to_numeric(april.mis_z, errors="coerce") < 2)]
    a_in_index = sum(_index_row(index, r.gene, r.hgvsp, "")["own_record_among_them"] for r in a_can.itertuples())
    rows += [("new ClinVar pathogenic or likely pathogenic missense variants in genes where BP1 can fire (pLI at least 0.9, missense z below 2)", len(can)),
             ("  with their own residue change among the gene's pathogenic missense records in the April protein index", in_index),
             ("  of them in a gene with exactly 3 such records, so that leaving the own record out lifts the veto", at_bar),
             ("    of those, with BP4 (BP1 with BP4 gives Likely benign, BP1 alone gives VUS)", at_bar_bp4),
             ("April ClinVar pathogenic or likely pathogenic missense variants in genes where BP1 can fire", len(a_can)),
             ("  April set: with their own residue change among the gene's pathogenic missense records in the April protein index", a_in_index)]
    write("temporal_new_in_index_v37.csv", pd.DataFrame(rows, columns=["measure", "n"]))


def proposals() -> None:
    """The changes that Section 4.5 proposes, on the variants of the October release that were not in the April set.

    Each run of the `holdout` group of `11_manuscript_runs.py` differs from `v37_oct2026` (the shipped classifier on that release)
    in one change that was measured on the April set, where it was found by reading that set's discordances and the sources. The
    `new` stratum was not scored by any sweep, review round or run, so the change a proposal makes there is the nearest thing to
    a test of it.
    The per-variant table of each run is made by `12_per_variant_table.py` from the run's annotated VCF.
    """
    st = BENCH / "clinvar_20261004_covered/strata.tsv"
    base_path = BENCH / "manuscript_runs/v37_oct2026_per_variant.tsv.gz"
    changes = [
        ("the completed ClinGen file", "v37_oct_clingen_complete", "v37_gap_clingen_complete"),
        ("PVS1 licensed by Orphanet's loss-of-function label", "v37_oct_orphanet_lof", "v37_gap_orphanet_lof"),
        ("both inputs of the PVS1 gap", "v37_oct_both", "v37_gap_both"),
        ("BP7 as Walker 2023 describes it", "v37_oct_bp7_walker", "v37_hk_bp7_walker"),
        ("BP1's veto without the variant's own record", "v37_oct_bp1_loo", "v37_hk_bp1_loo"),
        ("BP1's veto at a count of 1", "v37_oct_bp1_veto1", "v37_sw_bp1veto_1"),
        ("no allele-number floor", "v37_oct_no_an_floor", "v37_sw_no_an_floor"),
    ]
    if not (st.exists() and base_path.exists()) or not all((BENCH / f"manuscript_runs/{h}_per_variant.tsv.gz").exists() for _, h, _ in changes):
        print("  (the held-out runs have no per-variant tables yet: skipping)")
        return
    key = ["chrom", "pos", "ref", "alt"]
    strata = pd.read_csv(st, sep="\t", dtype=str)

    def load(path: Path) -> pd.DataFrame:
        x = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        x = x.merge(strata, on=key, how="left")
        x = x[x.stratum == "new"].copy()
        x["cv3"] = x["clinvar"].map(TRI)
        x["fv3"] = x["call"].map(TRI)
        return x.set_index(key)

    def metrics(x: pd.DataFrame) -> dict:
        return {"same": int((x.cv3 == x.fv3).sum()), "missed": int(((x.cv3 == "P/LP") & (x.fv3 == "B/LB")).sum()),
                "fp": int(((x.cv3 == "B/LB") & (x.fv3 == "P/LP")).sum()),
                "plp_correct": int(((x.cv3 == "P/LP") & (x.fv3 == "P/LP")).sum()), "blb_correct": int(((x.cv3 == "B/LB") & (x.fv3 == "B/LB")).sum()),
                "vus_called_plp": int(((x.cv3 == "VUS") & (x.fv3 == "P/LP")).sum()), "vus_called_blb": int(((x.cv3 == "VUS") & (x.fv3 == "B/LB")).sum())}

    base = load(base_path)
    mb = metrics(base)
    runs = pd.read_csv(OUT / "runs_v37_variants.csv").set_index("name")
    an = pd.to_numeric(base["gnomad_an"].replace("", np.nan), errors="coerce")
    popmax = pd.to_numeric(base["gnomad_popmax_af"].replace("", np.nan), errors="coerce")
    # The floor only matters at a site with a record (an absent record has no frequency to withhold) and a frequency at or above the
    # lowest bar that it gates, BS1's 5e-3, so the run tests it on the variants that are both.
    rows = [{"change": "none (v37 on the release)", "holdout_run": "v37_oct2026", "april_run": "v37", "new_variants": len(base),
             "new_same_direction": mb["same"], "new_opposite": mb["missed"] + mb["fp"], "new_false_benign": mb["missed"], "new_false_pathogenic": mb["fp"],
             "new_with_allele_number_below_2000": int((an < 2000).sum()),
             "new_below_2000_at_or_above_the_bs1_bar": int(((an < 2000) & (popmax >= 5e-3)).sum())}]
    for label, hold, apr in changes:
        x = load(BENCH / f"manuscript_runs/{hold}_per_variant.tsv.gz").reindex(base.index)
        if x["call"].isna().any():
            sys.exit(f"{hold} lacks calls for some new variants")
        m = metrics(x)
        a, a0 = runs.loc[apr], runs.loc["v37"]
        rows.append({"change": label, "holdout_run": hold, "april_run": apr, "new_variants": len(x),
                     "april_delta_same_direction": int(a.same_direction - a0.same_direction),
                     "april_delta_false_benign": int(a.false_benign - a0.false_benign),
                     "april_delta_false_pathogenic": int(a.false_pathogenic - a0.false_pathogenic),
                     "new_same_direction": m["same"], "new_delta_same_direction": m["same"] - mb["same"],
                     "new_opposite": m["missed"] + m["fp"], "new_false_benign": m["missed"], "new_false_pathogenic": m["fp"],
                     "new_delta_false_benign": m["missed"] - mb["missed"], "new_delta_false_pathogenic": m["fp"] - mb["fp"],
                     "new_delta_plp_correct": m["plp_correct"] - mb["plp_correct"], "new_delta_blb_correct": m["blb_correct"] - mb["blb_correct"],
                     "new_delta_vus_called_plp": m["vus_called_plp"] - mb["vus_called_plp"],
                     "new_delta_vus_called_blb": m["vus_called_blb"] - mb["vus_called_blb"],
                     "new_calls_changed": int((x["call"] != base["call"]).sum())})
    write("temporal_proposals_v37.csv", rows)


STEPS = {
    "headline": headline, "confusion": confusion, "stars": by_stars, "consequence": by_consequence,
    "offdiagonal": offdiagonal, "points": points_tables, "replay": replay_tables, "pm2dep": pm2_dependence, "criteria": criteria_table, "facts": facts,
    "opposite": opposite_cases, "genes": gene_level, "pvs1gap": pvs1_gap, "baselines": baselines,
}
NO_DF = {"vcep": vcep_bars, "gaprun": gap_runs, "hooks": hook_runs, "prov": provenance_summary, "trajectory": trajectory, "runs": runs_tables, "steps": threshold_steps, "rounds": review_rounds, "temporal": temporal, "newindex": new_in_index, "proposals": proposals}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", default="", help="comma-separated step names")
    args = ap.parse_args()
    want = [s for s in args.only.split(",") if s] or list(STEPS) + list(NO_DF)
    df = None
    for name in want:
        print(name)
        if name in NO_DF:
            NO_DF[name]()
        elif name in STEPS:
            if df is None:
                df = load()
            STEPS[name](df)
        else:
            sys.exit(f"unknown step {name}; choose from {', '.join(list(STEPS) + list(NO_DF))}")


if __name__ == "__main__":
    main()
