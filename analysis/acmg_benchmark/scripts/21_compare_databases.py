#!/usr/bin/env python3
"""What the complete databases changed: the runs with the first databases set against the same runs with the complete ones.

The benchmark was first built with SpliceAI and PhyloP databases that lack records which the gnomAD file they were distilled
from carries (`results/complete_databases/README.md`). Every run of the paper was made again with databases built from the
whole of that file, by the same classifiers on the same inputs, and the earlier runs were kept under the names
`manuscript_runs_april_databases/` and `output_v37_april_databases/`. This script reads the per-variant tables of the two sets of
runs and writes four things, so that the text states no number of the change that a file does not hold.

  databases_presence.csv           for each source (SpliceAI, PhyloP, gnomAD, REVEL), how many variants have a value under each
                                   stack, how many only under the complete one, and in how many the two values differ
  databases_causes.csv             for SpliceAI and PhyloP, the variants that only the complete databases have, by what explains
                                   them: the tail of a batch of the extraction that ended early, or an indel whose second base
                                   lay just outside the extracted regions (see below)
  databases_lost_blocks.csv        one row for each batch of the extraction in which the first databases stop short: the variants
                                   of each source that have a value before the break, the variants lost after it, and where
  databases_old_vs_new.csv         v37 and v38, on the April set and on the October release, with the first and with the
                                   complete databases: the headline counts and the firing of BP4, BP7, PP3 and of PVS1 on
                                   canonical splice variants (long format: version, set, measure, first, complete)
  databases_opposite_changes.csv   the opposite-direction calls that differ between the two stacks, with the criteria on each side

and the report `results/complete_vs_april_databases.md`, which is made from those files.

Why the first databases lack records. `extract_spliceai_phylop_from_gnomad.py` read gnomAD in batches of 200 regions of
`clinvar_2star_regions.bed` (the regions of each chromosome, in file order), ignored the exit status of `tabix` and sent its messages
to /dev/null. A batch whose stream ended early would therefore have left the databases without every record after the break and
no message. That is what the first databases show: every variant that only the complete databases hold lies in a run of one
batch that reaches the end of the batch, apart from indels whose second base (position + 1) is outside every region. The pattern
is what a stream that ends early gives, and this script does not show that a stream did: the extraction was not repeated.

Inputs (all under `data/benchmark/`): `manuscript_runs/` and `manuscript_runs_april_databases/` (the per-variant tables of `v37`,
`v38`, `v37_oct2026` and `v38_oct_shipped`), `clinvar_2star_labels.tsv`, `clinvar_20261004/clinvar_2star_labels.tsv` and
`clinvar_20261004_covered/clinvar_2star_labels.tsv`, `clinvar_20261004/strata.tsv`, and `sa_sources/clinvar_2star_regions.bed`.

Usage
  21_compare_databases.py            write the CSV files and the report
"""

from __future__ import annotations

import bisect
import importlib.util
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BENCH = ROOT / "data/benchmark"
OUT = ROOT / "manuscript/acmg/data"
REPORT = ROOT / "analysis/acmg_benchmark/results/complete_vs_april_databases.md"
NEW_RUNS = BENCH / "manuscript_runs"
OLD_RUNS = BENCH / "manuscript_runs_april_databases"
BED = BENCH / "sa_sources/clinvar_2star_regions.bed"
BATCH = 200  # extract_spliceai_phylop_from_gnomad.py: `for batch_start in range(0, len(regions), 200)`

# `13_manuscript_stats.py` holds the matrix and the metrics of the paper, and this script uses them as they are.
_spec = importlib.util.spec_from_file_location("manuscript_stats", HERE / "13_manuscript_stats.py")
st = importlib.util.module_from_spec(_spec)
sys.modules["manuscript_stats"] = st
_spec.loader.exec_module(st)

KEY = ["chrom", "pos", "ref", "alt"]
COLS = KEY + ["gene", "clinvar", "call", "criteria", "consequence", "spliceai_max", "phylop", "revel", "gnomad_af", "gnomad_ac", "gnomad_an", "gnomad_hc", "gnomad_popmax_af"]
SOURCES = {  # the name, the column that says whether a variant has a value, and the columns whose values are compared
    "SpliceAI": ("spliceai_max", ["spliceai_max"]),
    "PhyloP": ("phylop", ["phylop"]),
    "gnomAD": ("gnomad_an", ["gnomad_af", "gnomad_ac", "gnomad_an", "gnomad_hc", "gnomad_popmax_af"]),
    "REVEL": ("revel", ["revel"]),
}
SPLICE_CANONICAL = ("splice_donor_variant", "splice_acceptor_variant")
LABELS = {  # the labels of each set of variants: the ClinVar classes that the reference-agreement records (no call) are counted in
    "april": BENCH / "clinvar_2star_labels.tsv",
    "october": BENCH / "clinvar_20261004/clinvar_2star_labels.tsv",
    "october_covered": BENCH / "clinvar_20261004_covered/clinvar_2star_labels.tsv",
}
RUNS = {  # the per-variant table of the main run of each version on each set
    ("v37", "april"): "v37_per_variant.tsv.gz",
    ("v38", "april"): "v38_per_variant.tsv.gz",
    ("v37", "october"): "v37_oct2026_per_variant.tsv.gz",
    ("v38", "october"): "v38_oct_shipped_per_variant.tsv.gz",
}


def read(root: Path, name: str) -> pd.DataFrame:
    x = pd.read_csv(root / name, sep="\t", dtype=str, keep_default_na=False, usecols=COLS)
    if x.duplicated(KEY).any():
        sys.exit(f"{root / name} holds a variant twice")
    return x.set_index(KEY, drop=False)


def has_criterion(criteria: pd.Series, code: str) -> pd.Series:
    """The variants whose criteria include `code` at any strength (`PVS1`, `PVS1_Strong` and `PVS1_Moderate` are PVS1)."""
    return criteria.str.contains(rf"(?:^|&){code}(?:_|&|$)", regex=True)


def matrix(x: pd.DataFrame, labels: Path, stratum: str | None = None) -> dict[str, dict[str, int]]:
    """The five-tier matrix of a per-variant table with the reference-agreement records that a run did not annotate counted as no-calls.

    The records that a run did not annotate are counted from the labels of the set, or of one stratum of the release
    (`clinvar_20261004/strata.tsv`), so that the matrix has the denominator that the tables of the paper have.
    """
    m = {t: {c: 0 for c in st.CALLS} for t in st.CLASSES}
    for (t, c), k in x.groupby(["clinvar", "call"]).size().items():
        m[t][c] += int(k)
    lab = pd.read_csv(labels, sep="\t", dtype=str, usecols=KEY + ["normalized_class"])
    if stratum is not None:
        strata = pd.read_csv(BENCH / "clinvar_20261004/strata.tsv", sep="\t", dtype=str)
        lab = lab.merge(strata, on=KEY, how="left")
        lab = lab[lab["stratum"] == stratum]
    counts = lab["normalized_class"].value_counts()
    for t in st.CLASSES:
        short = int(counts[t]) - sum(m[t].values())
        if short < 0:
            sys.exit(f"the table has more {t} records than the labels {labels}")
        m[t]["NoCall"] += short
    return m


def headline(x: pd.DataFrame, labels: Path, stratum: str | None = None) -> dict[str, float]:
    mm = st.matrix_metrics(matrix(x, labels, stratum))
    return {k: mm[k] for k in ("n",) + tuple(MEASURES[1:])} | {"variants": mm["n"]}


def firing(x: pd.DataFrame) -> dict[str, int]:
    crit = x["criteria"]
    canonical = x["consequence"].isin(SPLICE_CANONICAL)
    pvs1 = has_criterion(crit, "PVS1")
    return {"BP4": int(has_criterion(crit, "BP4").sum()), "BP7": int(has_criterion(crit, "BP7").sum()), "PP3": int(has_criterion(crit, "PP3").sum()),
            "PVS1": int(pvs1.sum()), "canonical_splice_variants": int(canonical.sum()), "PVS1_on_canonical_splice": int((pvs1 & canonical).sum())}


# ------------------------------------------------------------------------------------------------ the sets of variants


def sets_of_runs() -> dict[tuple[str, str], dict]:
    """The tables of both stacks for each version and set, each restricted to the variants of one set of variants.

    The October runs of the first stack hold the 772,585 variants that a gnomAD extract covered (`clinvar_20261004_covered`), and
    the runs of the complete stack hold all 781,135. A comparison of the stacks is over the variants that both hold; the complete
    stack's other 8,550 are counted apart.
    """
    out = {}
    for (version, kind), table in RUNS.items():
        if not (NEW_RUNS / table).exists() or not (OLD_RUNS / table).exists():
            print(f"  ({version}, {kind}): a table is missing, skipped ({table})")
            continue
        out[(version, kind)] = {"first": read(OLD_RUNS, table), "complete": read(NEW_RUNS, table)}
    return out


# ------------------------------------------------------------------------------------------------ presence


def with_alt(x: pd.DataFrame) -> pd.DataFrame:
    """The records that have an alternative allele.

    The 382 reference-agreement records (ALT is a dot) get no call under either stack, and the table of v37 holds a row for each
    where the table of v38 holds none, so a count of the variants that have a value for a source depends on which table is read
    unless they are left out: PhyloP, which is looked up by position, has a value for 312 of them under the first databases and 348
    under the complete ones in the table of v37, and SpliceAI, which is looked up by allele, for none.
    """
    return x[x["alt"] != "."]


def any_pair(runs: dict, kind: str) -> dict | None:
    """The tables of both stacks of a set of variants, from either version (which databases a run reads does not depend on the version)."""
    return runs.get(("v37", kind)) or runs.get(("v38", kind))


def presence(runs: dict) -> list[dict]:
    rows = []
    for kind, scope, key_set in (("april", "April set", "april"), ("october", "October release, the variants of the first design", "covered")):
        pair = any_pair(runs, kind)
        if pair is None:
            continue
        first, complete = with_alt(pair["first"]), with_alt(pair["complete"])
        if key_set == "covered":
            complete = complete.loc[complete.index.isin(first.index)]
        if not first.index.equals(complete.index):
            complete = complete.reindex(first.index)
        for name, (present, compare) in SOURCES.items():
            f, c = first[present] != "", complete[present] != ""
            both = f & c
            differ = pd.Series(False, index=first.index)
            for col in compare:
                differ |= both & (first[col] != complete[col])
            rows.append({"set": scope, "source": name, "variants": len(first), "first": int(f.sum()), "complete": int(c.sum()), "only_complete": int((c & ~f).sum()),
                         "only_first": int((f & ~c).sum()), "both": int(both.sum()), "differ_where_both": int(differ.sum())})
    pair = any_pair(runs, "october")
    if pair is not None:
        outside = with_alt(pair["complete"]).loc[~with_alt(pair["complete"]).index.isin(pair["first"].index)]
        for name, (present, _) in SOURCES.items():
            rows.append({"set": "October release, the variants outside the first design", "source": name, "variants": len(outside), "first": "", "complete": int((outside[present] != "").sum()),
                         "only_complete": "", "only_first": "", "both": "", "differ_where_both": ""})
    return rows


# ------------------------------------------------------------------------------------------------ the causes


def load_regions() -> dict[str, list[tuple[int, int]]]:
    regions: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for line in BED.read_text().splitlines():
        c, s, e = line.split("\t")[:3]
        regions[c].append((int(s) + 1, int(e)))  # the extraction asks tabix for `chrom:start+1-end`
    for c, rs in regions.items():
        if rs != sorted(rs) or any(a[1] >= b[0] for a, b in zip(rs, rs[1:])):
            sys.exit(f"the regions of {c} are not sorted and disjoint, which the batches below assume")
    return regions


def causes(runs: dict) -> tuple[list[dict], list[dict]]:
    """For SpliceAI and PhyloP, why the variants that only the complete databases hold are missing from the first ones.

    A batch of the extraction has a tail when the first database has no value for any variant after some point of the batch while
    the complete database has values for variants after it: the variants of the batch after the last one that the first database
    has a value for, among those that the complete database has a value for, are the tail. The two sources are written from one
    `tabix` stream per batch by one script, so a stream that ended early would leave the same tail in both.
    """
    pair = runs[("v37", "april")]
    first = with_alt(pair["first"])
    complete = with_alt(pair["complete"]).reindex(first.index)
    regions = load_regions()
    starts = {c: [r[0] for r in rs] for c, rs in regions.items()}

    def region_of(chrom: str, pos: int) -> int:
        rs = regions["chr" + chrom]
        i = bisect.bisect_right(starts["chr" + chrom], pos) - 1
        return i if i >= 0 and rs[i][0] <= pos <= rs[i][1] else -1

    base = pd.DataFrame({"chrom": first["chrom"].values, "pos": first["pos"].astype(int).values, "ref": first["ref"].values, "alt": first["alt"].values})
    base["region"] = [region_of(c, p) for c, p in zip(base.chrom, base.pos)]
    base["batch"] = base.region // BATCH
    base["indel"] = base.ref.str.len() != base.alt.str.len()
    base["second_base_outside"] = [region_of(c, p + 1) < 0 for c, p in zip(base.chrom, base.pos)]
    rows = []
    per_batch: dict[tuple[str, int], dict] = {}
    for source, (present, _) in list(SOURCES.items())[:2]:
        f, c = (first[present] != "").values, (complete[present] != "").values
        lost = c & ~f
        tail = pd.Series(False, index=range(len(base)))
        n_tails = 0
        for (chrom, batch), part in base[(base.region >= 0)].groupby(["chrom", "batch"]):
            idx = part.index.values
            if not lost[idx].any():
                continue
            kept = part.pos.values[f[idx]]
            last_kept = int(kept.max()) if len(kept) else 0
            after = part.pos.values > last_kept
            in_tail = lost[idx] & after
            if not in_tail.any():
                continue
            tail.iloc[idx[in_tail]] = True
            n_tails += 1
            nreg = len(regions["chr" + chrom])
            first_region = int(batch) * BATCH
            entry = per_batch.setdefault((chrom, int(batch) + 1), {"chrom": chrom, "batch": int(batch) + 1, "regions_in_batch": min(first_region + BATCH, nreg) - first_region})
            entry[f"{source.lower()}_with_a_value_before_the_break"] = int(f[idx][part.pos.values <= last_kept].sum())
            entry[f"{source.lower()}_lost"] = int(in_tail.sum())
            entry[f"{source.lower()}_first_position_lost"] = int(part.pos.values[in_tail].min())
            entry[f"{source.lower()}_last_position_lost"] = int(part.pos.values[in_tail].max())
        tail_v = tail.values
        edge = lost & ~tail_v & base.indel.values & base.second_base_outside.values
        other = lost & ~tail_v & ~edge
        rows.append({"source": source, "variants_lost": int(lost.sum()), "in_the_tail_of_a_batch": int((lost & tail_v).sum()), "batches_with_a_tail": n_tails,
                     "indels_with_the_second_base_outside_the_regions": int(edge.sum()), "unexplained": int(other.sum()),
                     "indels_among_the_unexplained": int((other & base.indel.values).sum()),
                     "variants_with_a_value_in_the_complete": int(c.sum()), "share_lost_of_the_complete": float(lost.sum() / c.sum())})
    blocks = [per_batch[k] for k in sorted(per_batch, key=lambda k: (k[0].rjust(2, "0"), k[1]))]
    cols = ["chrom", "batch", "regions_in_batch"] + [f"{src}_{what}" for src in ("spliceai", "phylop") for what in ("with_a_value_before_the_break", "lost", "first_position_lost", "last_position_lost")]
    blocks = [{c: b.get(c, "") for c in cols} for b in blocks]
    return rows, blocks


# ------------------------------------------------------------------------------------------------ old against new


MEASURES = ["variants", "exact", "same_direction", "same_pathogenic", "same_benign", "vus_exact", "opposite", "false_benign", "false_pathogenic", "plp_recall", "blb_recall", "vus_recall", "plp_to_vus", "blb_to_vus"]
FIRING = ["BP4", "BP7", "PP3", "PVS1", "canonical_splice_variants", "PVS1_on_canonical_splice"]


def old_vs_new(runs: dict) -> list[dict]:
    rows = []
    strata = pd.read_csv(BENCH / "clinvar_20261004/strata.tsv", sep="\t", dtype=str).set_index(KEY)["stratum"] if (BENCH / "clinvar_20261004/strata.tsv").exists() else None
    for (version, kind), pair in runs.items():
        first, complete = pair["first"], pair["complete"]
        if kind == "april":
            groups = [("April set", first, complete, LABELS["april"], LABELS["april"], None)]
        else:
            if strata is None:
                continue
            inside = complete.loc[complete.index.isin(first.index)]
            new_first = first.loc[strata.reindex(first.index).values == "new"]
            new_inside = inside.loc[strata.reindex(inside.index).values == "new"]
            groups = [("October release, the variants of the first design", first, inside, LABELS["october_covered"], LABELS["october_covered"], None),
                      ("October release, new variants of the first design", new_first, new_inside, LABELS["october_covered"], LABELS["october_covered"], "new"),
                      ("October release, all variants", None, complete, None, LABELS["october"], None)]
        for scope, f, c, lab_f, lab_c, stratum in groups:
            vals = {}
            for db, x, lab in (("first", f, lab_f), ("complete", c, lab_c)):
                if x is None:
                    vals[db] = None
                    continue
                h = headline(x, lab, stratum)
                h.update(firing(x))
                vals[db] = h
            for measure in MEASURES + FIRING:
                a, b = (vals[db].get(measure) if vals[db] is not None else None for db in ("first", "complete"))
                if a is None and b is None:
                    continue
                rows.append({"version": version, "set": scope, "measure": measure, "first": "" if a is None else a, "complete": "" if b is None else b})
    return rows


def opposite_changes(runs: dict) -> list[dict]:
    rows = []
    for version in ("v37", "v38"):
        pair = runs.get((version, "april"))
        if pair is None:
            continue
        first, complete = pair["first"], pair["complete"].reindex(pair["first"].index)

        def opposite(x):
            cv, fv = x["clinvar"].map(st.TRI), x["call"].map(st.TRI)
            return ((cv == "P/LP") & (fv == "B/LB")) | ((cv == "B/LB") & (fv == "P/LP"))

        of, oc = opposite(first), opposite(complete)
        for key in first.index[of | oc]:
            if of[key] and oc[key] and first.at[key, "call"] == complete.at[key, "call"]:
                continue
            a, b = first.loc[key], complete.loc[key]
            rows.append({"version": version, "gene": b["gene"], "chrom": key[0], "pos": key[1], "ref": key[2], "alt": key[3], "clinvar": b["clinvar"], "consequence": b["consequence"],
                         "first_call": a["call"], "complete_call": b["call"], "first_criteria": a["criteria"], "complete_criteria": b["criteria"],
                         "change": "enters the opposite-direction calls" if oc[key] and not of[key] else ("leaves them" if of[key] and not oc[key] else "stays, call differs"),
                         "first_spliceai_max": a["spliceai_max"], "complete_spliceai_max": b["spliceai_max"], "first_phylop": a["phylop"], "complete_phylop": b["phylop"]})
    return rows


# ------------------------------------------------------------------------------------------------ the report


def n(x) -> str:
    return f"{int(round(float(x))):,}"


def pct(x: float, d: int = 2) -> str:
    return f"{float(x) * 100:.{d}f}%"


def signed(x: float) -> str:
    return f"{int(x):+,}"


def report(pres: list[dict], cause_rows: list[dict], blocks: list[dict], ovn: list[dict], opp: list[dict]) -> None:
    pres_df, cause_df, ovn_df, blocks_df = pd.DataFrame(pres), pd.DataFrame(cause_rows), pd.DataFrame(ovn), pd.DataFrame(blocks)
    lines: list[str] = []

    def s(text: str = "") -> None:
        lines.append(text)

    def get(version: str, scope: str, measure: str, db: str) -> float:
        v = ovn_df[(ovn_df.version == version) & (ovn_df.set == scope) & (ovn_df.measure == measure)][db]
        return float(v.iloc[0])

    s("# What the complete databases changed")
    s()
    s("The benchmark was first built with SpliceAI and PhyloP databases that lack records which the gnomAD file they were distilled from carries (`results/complete_databases/README.md`).")
    s("Every run of the paper was then made again with databases built from the whole of that file, by the same classifiers on the same inputs, and the runs on the first databases were kept.")
    s("This file sets the two side by side.")
    s("It is written by `scripts/21_compare_databases.py` from the per-variant tables of the runs, and the numbers are in `manuscript/acmg/data/databases_*.csv`.")
    s("REVEL, ClinVar, ClinGen, the ClinVar protein index and RepeatMasker are the same files under both stacks.")
    s()
    s("## What the first databases lacked")
    s()
    april = pres_df[pres_df.set == "April set"]
    s("| Source | Variants with a value under the first databases | Under the complete databases | Only under the complete | Only under the first | Differing where both have one |")
    s("|:---|---:|---:|---:|---:|---:|")
    for r in april.itertuples():
        s(f"| {r.source} | {n(r.first)} | {n(r.complete)} | {n(r.only_complete)} | {n(r.only_first)} | {n(r.differ_where_both)} |")
    s()
    s(f"The table is over the {n(april.iloc[0].variants)} variants of the April set that have an alternative allele (the other 382 are reference-agreement records, which get no call).")
    sa, pp = (cause_df[cause_df.source == x].iloc[0] for x in ("SpliceAI", "PhyloP"))
    s(f"The complete SpliceAI database has a value for {n(sa.variants_with_a_value_in_the_complete)} of them and the first for {n(sa.variants_with_a_value_in_the_complete - sa.variants_lost)}, so the first lacks {pct(sa.share_lost_of_the_complete, 1)} of the records that exist; for PhyloP the shares are {n(pp.variants_with_a_value_in_the_complete)} and {n(pp.variants_with_a_value_in_the_complete - pp.variants_lost)}, and {pct(pp.share_lost_of_the_complete, 1)}.")
    s("No variant has a value only under the first databases, and where both have one the values are the same, so the complete databases add records and change none.")
    s("gnomAD and REVEL have the same variants with a value, and the same values, under both stacks.")
    s()
    s("## Where the missing records are")
    s()
    s(f"The extraction of April read gnomAD in batches of {BATCH} regions of `clinvar_2star_regions.bed`, a chromosome at a time, ignored the exit status of `tabix` and sent its messages to /dev/null.")
    s(f"The {n(sa.variants_lost)} SpliceAI records that only the complete database holds lie in {int(sa.batches_with_a_tail)} stretches, each inside one batch and reaching its end: in each of those batches the first database has values up to a point and none after it, and the complete database has values after it.")
    rest = int(pp.variants_lost - pp.in_the_tail_of_a_batch)
    edge = int(pp.indels_with_the_second_base_outside_the_regions)
    s(f"For PhyloP the same is true of {n(pp.in_the_tail_of_a_batch)} of the {n(pp.variants_lost)} records, in {int(pp.batches_with_a_tail)} batches.")
    if rest == edge and not pp.unexplained:
        s(f"The other {n(rest)} are all indels, and the second base (position + 1) of each lies outside every region of the extraction.")
    else:
        s(f"The other {n(rest)} are indels or not, and for {n(edge)} of them the second base (position + 1) lies outside every region of the extraction, which leaves {n(pp.unexplained)} that neither explanation covers.")
    both = blocks_df[(blocks_df.spliceai_lost != "") & (blocks_df.phylop_lost != "")]
    same_start = int((both.spliceai_first_position_lost == both.phylop_first_position_lost).sum())
    s(f"The batches are the same for the two sources ({len(both)} batches lose records in both of {len(blocks_df)}), and the first lost position is the same in {same_start} of them, as it would be if one stream per batch fed both files, which is how the extraction script works (one `tabix` stream is read once and written to both).")
    s("The pattern is what a stream that ends early gives, and it is not a record of a stream that did: the extraction was not repeated, and the script threw its messages away.")
    s("`databases_lost_blocks.csv` lists each batch, with the variants of each source that have a value before the break and those lost after it.")
    s()
    s("## The April set")
    s()
    scope = "April set"
    s("| Measure | v37, first databases | v37, complete | v38, first databases | v38, complete |")
    s("|:---|---:|---:|---:|---:|")

    def line(label: str, measure: str, fmt) -> None:
        cells = [fmt(get(v, scope, measure, db)) for v in ("v37", "v38") for db in ("first", "complete")]
        s(f"| {label} | " + " | ".join(cells) + " |")

    line("Exact agreement", "exact", n)
    line("Agreement in direction", "same_direction", n)
    line("Opposite-direction calls", "opposite", n)
    line("Reference-discordant benign calls", "false_benign", n)
    line("Reference-discordant pathogenic calls", "false_pathogenic", n)
    line("Recall of ClinVar pathogenic-side variants", "plp_recall", lambda x: pct(x, 1))
    line("Recall of ClinVar benign-side variants", "blb_recall", lambda x: pct(x, 1))
    line("Recall of ClinVar VUS", "vus_recall", lambda x: pct(x, 1))
    line("BP4 fires on", "BP4", n)
    line("BP7 fires on", "BP7", n)
    line("PP3 fires on", "PP3", n)
    line("PVS1 fires on", "PVS1", n)
    line("PVS1 fires on canonical splice variants", "PVS1_on_canonical_splice", n)
    s()
    s(f"There are {n(get('v38', scope, 'canonical_splice_variants', 'complete'))} canonical splice donor and acceptor variants in the set.")
    for v in ("v38", "v37"):
        d_dir = get(v, scope, "same_direction", "complete") - get(v, scope, "same_direction", "first")
        d_exact = get(v, scope, "exact", "complete") - get(v, scope, "exact", "first")
        n_all = get(v, scope, "variants", "first")
        s(f"Under {v} the complete databases give {signed(d_exact)} exact calls ({d_exact / n_all * 100:+.2f} points) and {signed(d_dir)} calls in direction ({d_dir / n_all * 100:+.2f} points), and the opposite-direction calls go from {n(get(v, scope, 'opposite', 'first'))} to {n(get(v, scope, 'opposite', 'complete'))}.")
        s(f"Recall of ClinVar benign-side variants goes from {pct(get(v, scope, 'blb_recall', 'first'), 1)} to {pct(get(v, scope, 'blb_recall', 'complete'), 1)} and recall of pathogenic-side variants from {pct(get(v, scope, 'plp_recall', 'first'), 1)} to {pct(get(v, scope, 'plp_recall', 'complete'), 1)}.")
    s()
    changed = pd.DataFrame(opp)
    if len(changed):
        s("### The opposite-direction calls that differ")
        s()
        for v in ("v38", "v37"):
            part = changed[changed.version == v]
            enter, leave, stay = (part[part.change == k] for k in ("enters the opposite-direction calls", "leaves them", "stays, call differs"))
            s(f"Under {v}, {len(enter)} variants enter the opposite-direction calls, {len(leave)} leave them and {len(stay)} stay in them with another call.")
            s()
            for r in enter.itertuples():
                s(f"- {r.gene} {r.chrom}:{r.pos} {r.ref}>{r.alt} (ClinVar {r.clinvar}, {r.consequence}): {r.first_call} ({r.first_criteria or 'no criterion'}) becomes {r.complete_call} ({r.complete_criteria or 'no criterion'}); SpliceAI {r.first_spliceai_max or 'none'} to {r.complete_spliceai_max or 'none'}, PhyloP {r.first_phylop or 'none'} to {r.complete_phylop or 'none'}.")
            for r in leave.itertuples():
                s(f"- {r.gene} {r.chrom}:{r.pos} {r.ref}>{r.alt} (ClinVar {r.clinvar}, {r.consequence}) leaves: {r.first_call} ({r.first_criteria or 'no criterion'}) becomes {r.complete_call} ({r.complete_criteria or 'no criterion'}).")
            s()
        s()
    october = ovn_df[ovn_df.set.str.startswith("October")]
    if len(october):
        s("## The October release")
        s()
        s("The first design kept the 772,585 variants of the release that a gnomAD extract covered, and the complete databases need no such rule, so the release is scored whole (781,135 variants at two stars or better).")
        s("The rows below set the two stacks side by side over the variants that both hold, and the last rows give the complete stack on the whole release.")
        s()
        s("| Set | Measure | v37, first | v37, complete | v38, first | v38, complete |")
        s("|:---|:---|---:|---:|---:|---:|")
        for sc in october.set.drop_duplicates():
            for measure, label, fmt in (("variants", "Variants", n), ("same_direction", "Agreement in direction", n), ("opposite", "Opposite-direction calls", n)):
                part = october[(october.set == sc) & (october.measure == measure)]
                if part.empty:
                    continue
                cells = []
                for v in ("v37", "v38"):
                    row = part[part.version == v]
                    for db in ("first", "complete"):
                        val = row.iloc[0][db] if len(row) else ""
                        cells.append(fmt(val) if val != "" else "-")
                s(f"| {sc} | {label} | " + " | ".join(cells) + " |")
        s()
    s("## What this file does not say")
    s()
    s("The first databases were not rebuilt, and the run of the extraction that left them short was not repeated, so the cause above is read from the pattern of the records and not observed.")
    s("The earlier numbered versions v1 to v36 were not repeated on the complete databases, and their results are on the first ones.")
    s("A comparison of two stacks says what the missing records changed in these runs, and it does not say which stack is closer to the truth, which only ClinVar's classifications can be compared with.")
    text = []
    for line in lines:  # one blank line between blocks, never two
        if line == "" and text and text[-1] == "":
            continue
        text.append(line)
    REPORT.write_text("\n".join(text) + "\n")
    print(f"  wrote {REPORT.relative_to(ROOT)} ({len(text)} lines)")


# ------------------------------------------------------------------------------------------------ main


def write(name: str, rows: list[dict]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT / name, index=False, lineterminator="\n")
    print(f"  wrote {name} ({len(rows)} rows)")


def main() -> None:
    runs = sets_of_runs()
    if ("v37", "april") not in runs:
        sys.exit("the per-variant tables of v37 on the April set are missing from one of the two roots")
    pres = presence(runs)
    cause_rows, blocks = causes(runs)
    ovn = old_vs_new(runs)
    opp = opposite_changes(runs)
    write("databases_presence.csv", pres)
    write("databases_causes.csv", cause_rows)
    write("databases_lost_blocks.csv", blocks)
    write("databases_old_vs_new.csv", ovn)
    write("databases_opposite_changes.csv", opp)
    report(pres, cause_rows, blocks, ovn, opp)


if __name__ == "__main__":
    main()
