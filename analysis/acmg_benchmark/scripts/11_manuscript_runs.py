#!/usr/bin/env python3
"""Run the threshold sweeps and leave-one-out ablations the ACMG manuscript reports.

Every run here is a variant of `v37`, the final classifier, made by
`run_benchmark.sh` so that it is scored by `03_evaluate_concordance.py` and records
what produced it. A variant of a run is a run of its own, so each one has its own
version name, `v37_<knob>_<value>`.

Three groups:

  sweep   one numeric threshold varied alone, everything else the shipped default.
          The shipped value is `v37` itself and is not re-run.
  switch  one boolean or list setting turned off (or, for the two opt-in settings,
          on), everything else the shipped default.
  source  one supplementary-annotation source removed from `--sa-dir`, by symlinking
          every other file into a directory of its own. The classifier degrades
          quietly when a source is missing (a criterion with no data reports
          `evaluated: false`), so the run shows what the source was buying.
  interaction  a sweep repeated under a different combining rule, to show whether the
          threshold's importance belongs to the threshold or to the rule beside it.
  cutoff  the SVI point cutoffs of the pathogenic side (Likely pathogenic from 6, Pathogenic
          from 10), which are constants in the code and not settings. These runs need a
          binary built from `manuscript/acmg/replay_validation/combiner_cutoff_env.patch`,
          which reads them from the environment; name it in FASTVEP_CUTOFF_BIN. They exist
          to check the offline replay of `13_manuscript_stats.py` against real runs.
  hook    two departures found while reading the sources again, measured by real runs of a
          binary built from `manuscript/acmg/replay_validation/bp_hooks_env.patch` (BP1's veto
          counting a variant's own ClinVar record, and BP7's SpliceAI bar); name it in
          FASTVEP_HOOKS_BIN. With no environment variable set it is the shipped classifier,
          which the control run `v37_hk_control` checks cell for cell.
  gap     what the PVS1 gap of Section 3.5 is made of: the ClinGen file with the curations the
          April export lacked put back, and PVS1 licensed by Orphanet's loss-of-function label
          for genes that no other route reaches. Needs the files named in `gap_inputs()`.
  holdout the changes that Section 4.5 proposes, run on the variants of the ClinVar release of
          4 October 2026 (the covered subset that `14_temporal_holdout.py prepare` writes), so
          that they are tried on variants that none of them was tuned on. They are run with
          `--holdout`, and are not among the runs of `runs.tsv` (which are all on the April set):
          their matrices are in `holdout_runs.tsv`, and `13_manuscript_stats.py` scores them on
          the new variants from the per-variant tables of `12_per_variant_table.py`.

Outputs
  data/benchmark/manuscript_runs/<name>/        bulk output, gitignored. The annotated
                                                VCF is deleted after scoring unless the
                                                run is in KEEP_VCF: the scored files and
                                                the provenance are what the manuscript reads.
  analysis/acmg_benchmark/results/v37_manuscript/
      run_matrix.tsv       what each run is
      configs/<name>.toml  the --acmg-config of each run that needs one
      provenance/<name>.txt
      runs.tsv             one row per run: the 5x6 concordance matrix and the
                           firing count of every criterion, read back from the scored files

Usage
  11_manuscript_runs.py                 run everything that has not been run
  11_manuscript_runs.py --workers 2     two runs at a time
  11_manuscript_runs.py --only pm2ad    only runs whose name contains the text
  11_manuscript_runs.py --list          print the matrix and stop
  11_manuscript_runs.py --collect       rebuild runs.tsv from the scored files and stop
  11_manuscript_runs.py --holdout       the runs of the `holdout` group (the others are not run)

  FASTVEP_CUTOFF_BIN=/path/to/patched/fastvep 11_manuscript_runs.py --only cut

Every sweep and ablation is fitted to, and scored against, ClinVar's 2-star+
classifications, which are a reference and not ground truth.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BENCH = ROOT / "data/benchmark"
OUT_ROOT = BENCH / "manuscript_runs"
SA_DB = BENCH / "sa_db"
SA_ABLATE = BENCH / "sa_db_ablate"
RES = ROOT / "analysis/acmg_benchmark/results/v37_manuscript"
RUNNER = ROOT / "analysis/acmg_benchmark/scripts/run_benchmark.sh"

CLASSES = ["Pathogenic", "Likely_pathogenic", "VUS", "Likely_benign", "Benign"]
COLS = CLASSES + ["NoCall"]

# Runs whose annotated VCF is kept, because a later analysis reads it call by call.
KEEP_VCF = {
    "v37_src_noclinvar_protein",
    "v37_src_noclinvar_all",
    "v37_sw_clinvar_informed",
    "v37_sw_no_lowpen_gate",
    "v37_sw_table_both_sides",
    "v37_sw_phylop_1000",
    "v37_hk_bp1_loo",
    "v37_hk_bp7_walker",
    "v37_sw_bp1veto_1",
    "v37_sw_bp1veto_2",
    "v37_hk_bp1_loo_veto1",
    "v37_hk_bp1_loo_veto2",
    "v37_gap_clingen_complete",
    "v37_gap_orphanet_lof",
    "v37_gap_both",
}


@dataclass
class Run:
    name: str
    group: str  # sweep | switch | source | interaction | cutoff
    knob: str  # the setting varied, or the source removed
    value: str  # the value it was set to, as text
    shipped: str  # the shipped value, as text
    description: str
    toml: list[str] = field(default_factory=list)  # lines of the config, or empty
    sa_remove: list[str] = field(default_factory=list)  # glob patterns removed from sa_db
    sa_replace: dict[str, str] = field(default_factory=dict)  # sa_db file name -> file used instead
    env: dict[str, str] = field(default_factory=dict)  # environment the run needs (the cutoff runs)
    inputs: dict[str, str] = field(default_factory=dict)  # INPUT and LABELS of a run that is not on the April set
    binary: str = ""  # the patched build the run needs ("hook"), where its group does not say so


def sweeps() -> list[Run]:
    grids = [
        (
            "pm2_ad_af_threshold",
            "pm2ad",
            [0.0, 1e-6, 5e-6, 1e-5, 2e-5, 1e-4, 2e-4, 5e-4],
            "4e-05",
            "PM2 allele-frequency bar for dominant or unknown-inheritance genes",
        ),
        (
            "bs1_af_threshold",
            "bs1",
            [2e-4, 5e-4, 1e-3, 2e-3, 1e-2, 2e-2],
            "0.005",
            "BS1 allele-frequency bar",
        ),
        (
            "bs2_hom_prevalence_threshold",
            "bs2hom",
            [1e-6, 1e-5, 1e-4, 1e-2],
            "0.001",
            "BS2 prevalence bar for the individuals with no functional copy (every gene not read as dominant-only)",
        ),
        (
            "pm2_ar_af_threshold",
            "pm2ar",
            [0.0, 1e-6, 1e-5, 2e-5, 2e-4, 5e-4],
            "7e-05",
            "PM2 allele-frequency bar for genes read as recessive-only",
        ),
        (
            "bp1_max_pathogenic_missense",
            "bp1veto",
            [1, 2, 5, 10, 25],
            "3",
            "number of ClinVar pathogenic missense variants in a gene at which BP1's veto applies (no veto at all is v37_sw_no_bp1_veto)",
        ),
        (
            "phylop_conserved",
            "phylop",
            [1.0, 1.5, 3.0, 4.0, 1000.0],
            "2.0",
            "PhyloP score at or above which a nucleotide counts as conserved for BP7 (1000 switches the test off)",
        ),
    ]
    out: list[Run] = []
    for knob, short, values, shipped, what in grids:
        for v in values:
            label = f"{v:g}"
            out.append(
                Run(
                    name=f"v37_sw_{short}_{label}",
                    group="sweep",
                    knob=knob,
                    value=label,
                    shipped=shipped,
                    description=f"{what}, set to {label}",
                    toml=[f"{knob} = {v!r}"],
                )
            )
    for v in [0, 50, 100, 200, 500, 1000, 1_000_000]:
        out.append(
            Run(
                name=f"v37_sw_bp7far_{v}",
                group="sweep",
                knob="bp7_max_intron_offset",
                value=str(v),
                shipped="300",
                description=(
                    "BP7 far intron boundary, set to "
                    + str(v)
                    + (
                        " (intronic extension off)"
                        if v == 0
                        else " (effectively unbounded)"
                        if v == 1_000_000
                        else ""
                    )
                ),
                toml=[f"bp7_max_intron_offset = {v}"],
            )
        )
    return out


def switches() -> list[Run]:
    def sw(short, knob, value, shipped, what, toml=None):
        return Run(
            name=f"v37_sw_{short}",
            group="switch",
            knob=knob,
            value=value,
            shipped=shipped,
            description=what,
            toml=toml if toml is not None else [f"{knob} = {value}"],
        )

    return [
        sw(
            "pm2_moderate",
            "pm2_downgrade_to_supporting",
            "false",
            "true",
            "PM2 kept at Moderate, as in Richards 2015 (the ClinGen SVI September 2020 downgrade to Supporting off)",
        ),
        sw(
            "pm2_needs_record",
            "pm2_absent_when_no_record",
            "false",
            "true",
            "PM2 not evaluated when gnomAD holds no record for the variant (absence is not read as rarity)",
        ),
        sw(
            "table_both_sides",
            "use_point_system",
            "false",
            "true",
            "Richards 2015 combining table on both sides, in place of the SVI point system on the pathogenic side",
        ),
        sw(
            "no_gene_validity",
            "require_gene_disease_validity",
            "false",
            "true",
            "PVS1, PP2 and PM1 no longer withheld where every ClinGen curation of the gene is Limited, Disputed or Refuted",
        ),
        sw(
            "no_mechanism_gate",
            "mechanism_gates_pvs1",
            "false",
            "true",
            "PVS1 no longer withheld for genes curated as gain-of-function or dominant-negative",
        ),
        sw(
            "no_lowpen_gate",
            "clinvar_low_penetrance_blocks_benign_frequency",
            "false",
            "true",
            "ClinVar's low-penetrance and risk-allele terms no longer withhold BS1 and BS2 (a ClinVar-derived gate)",
        ),
        sw(
            "no_region_flags",
            "gnomad_region_flags_block_frequency",
            "false",
            "true",
            "gnomAD's segmental-duplication and low-complexity flags no longer withhold the frequency criteria",
        ),
        sw(
            "no_faf",
            "use_filtering_af",
            "false",
            "true",
            "BA1 and BS1 test the population-maximum point estimate, not the filtering allele frequency",
        ),
        sw(
            "clinvar_informed",
            "exclude_self_from_clinvar_evidence",
            "false",
            "true",
            "PS1 and PM1 count the variant's own ClinVar record (the ClinVar-informed mode)",
        ),
        sw(
            "no_homology_list",
            "homology_unreliable_genes",
            "[]",
            "29 genes",
            "no gene is treated as having unreliable population frequencies (Mandelker 2016 list emptied)",
        ),
        sw(
            "no_exceptions",
            "ba1_exceptions",
            "[]",
            "15 alleles",
            "the curated frequency-exception list emptied (Ghosh 2018's nine, three hypomorphs, three reviewer-named)",
        ),
        sw(
            "no_bp1_veto",
            "bp1_max_pathogenic_missense",
            "1000000",
            "3",
            "BP1 no longer vetoed in genes with an established pathogenic-missense spectrum",
        ),
        sw(
            "no_pm1_benign_half",
            "pm1_max_benign_in_window",
            "1000000",
            "0",
            "PM1 no longer requires the absence of benign variation in the window",
        ),
        sw(
            "no_an_floor",
            "min_an_for_frequency_criteria",
            "0",
            "2000",
            "BA1 and BS1 no longer require an allele number of at least 2,000 at the variant's site",
        ),
        sw(
            "pp3_cap_moderate",
            "pp3_max_strength",
            '"Moderate"',
            "uncapped",
            "PP3 capped at Moderate (the stricter convention the reviewer described)",
        ),
        sw(
            "pvs1_nmd50",
            "pvs1_nmd_50nt_rule",
            "true",
            "false",
            "PVS1 uses the exact 50-nucleotide NMD rule of Abou Tayoun 2018 (opt-in)",
        ),
    ]


def sources() -> list[Run]:
    def src(short, remove, what, toml=None):
        return Run(
            name=f"v37_src_{short}",
            group="source",
            knob=short,
            value="removed",
            shipped="loaded",
            description=what,
            toml=toml or [],
            sa_remove=remove,
        )

    return [
        src(
            "noclinvar_protein",
            ["clinvar_protein.oga"],
            "ClinVar protein and splice-position index removed (PS1, PM1, PM5, BP1 veto, PVS1 critical region)",
        ),
        src(
            "noclinvar_all",
            ["clinvar_protein.oga", "clinvar.osa", "clinvar.osa.idx"],
            "every ClinVar-derived source removed (the protein index and the allele-level record)",
        ),
        src("norevel", ["revel_chr*"], "REVEL removed (PP3 and BP4 for missense)"),
        src("nospliceai", ["spliceai_chr*"], "SpliceAI removed (PP3 and BP4 for splicing, BP7)"),
        src("nophylop", ["phylop_chr*"], "PhyloP removed (BP7's conservation test)"),
        src(
            "noconstraint",
            ["gnomad_genes.oga"],
            "gnomAD gene constraint removed (PVS1's loss-of-function-intolerance gate, PP2, BP1)",
        ),
        src(
            "noclingen",
            ["omim.oga"],
            "ClinGen gene-disease validity and inheritance removed (PVS1 disease gate, inheritance for BS2, PM2 and PM3)",
        ),
        src("norepeat", ["repeatmasker.osi"], "RepeatMasker removed (BP3)"),
        src(
            "nognomad",
            ["gnomad_chr*"],
            "gnomAD removed, with PM2 told not to read absence as rarity, so no frequency criterion can fire",
            toml=["pm2_absent_when_no_record = false"],
        ),
    ]


def interactions() -> list[Run]:
    """The PM2 sweep again, with the Richards 2015 table on both sides.

    PM2 was set at 4e-5 in v10, under the 2015 table, where a lone PVS1 matched no row and
    PM2 was the criterion that took a null variant to Likely pathogenic. The point system of
    v11 changed that. These runs hold the table fixed and move the PM2 bar, which is the
    sweep v10 made, on the final code.
    """
    out = []
    for v in (0.0, 1e-6, 1e-5, 1e-4):
        label = f"{v:g}"
        out.append(
            Run(
                name=f"v37_ix_table_pm2ad_{label}",
                group="interaction",
                knob="pm2_ad_af_threshold (Richards 2015 table both sides)",
                value=label,
                shipped="4e-05",
                description=f"PM2 bar for dominant genes set to {label}, with the Richards 2015 table on both sides",
                toml=["use_point_system = false", f"pm2_ad_af_threshold = {v!r}"],
            )
        )
    return out


def cutoffs() -> list[Run]:
    """The point cutoffs of the pathogenic side, by real runs of a patched binary.

    Likely pathogenic from 6 points and Pathogenic from 10 are constants in `combiner.rs`. The
    patch makes them read `FASTVEP_LP_CUT` and `FASTVEP_P_CUT`, so that what a different cutoff
    does is measured and not only replayed from the criteria a variant met.
    """
    out = []
    # 6 is the shipped cutoff: that run is the control, which must reproduce v37 exactly, because
    # it shows that the patched binary behaves as the shipped one does apart from the two cutoffs.
    for lp in (4, 5, 6, 7, 8, 9):
        out.append(
            Run(
                name=f"v37_cut_lp{lp}",
                group="cutoff",
                knob="likely_pathogenic_points",
                value=str(lp),
                shipped="6",
                description=(
                    f"Likely pathogenic from {lp} pathogenic points"
                    + (" (the shipped value: a control that must reproduce v37)" if lp == 6 else " (6 is shipped)")
                ),
                env={"FASTVEP_LP_CUT": str(lp)},
            )
        )
    for p in (9, 11, 12):
        out.append(
            Run(
                name=f"v37_cut_p{p}",
                group="cutoff",
                knob="pathogenic_points",
                value=str(p),
                shipped="10",
                description=f"Pathogenic from {p} pathogenic points (10 is shipped)",
                env={"FASTVEP_P_CUT": str(p)},
            )
        )
    return out


def hooks() -> list[Run]:
    """BP1 and BP7 measured by real runs of a build that exposes one setting of each.

    `v37_hk_control` has nothing set, and must reproduce v37 in every cell: that is what
    licenses reading the other three as the shipped classifier with one thing changed.
    """
    return [
        Run(
            name="v37_hk_control",
            group="hook",
            knob="measurement build, nothing set",
            value="-",
            shipped="-",
            description="the measurement build with no hook set (a control that must reproduce v37)",
        ),
        Run(
            name="v37_hk_bp1_loo",
            group="hook",
            knob="BP1 veto, own ClinVar record excluded",
            value="on",
            shipped="off",
            description=(
                "BP1's veto, which counts a gene's pathogenic missense records in the ClinVar protein "
                "index, no longer counts the variant's own record when it is itself pathogenic or likely "
                "pathogenic (PS1 and PM1 already exclude it)"
            ),
            env={"FASTVEP_BP1_LOO": "1"},
        ),
        Run(
            name="v37_hk_bp7_spliceai01",
            group="hook",
            knob="BP7 SpliceAI bar",
            value="0.1 or below",
            shipped="below 0.2",
            description="BP7 requires SpliceAI at most 0.1 (BP4's bar, as Walker 2023 requires) in place of below 0.2",
            env={"FASTVEP_BP7_SPLICEAI_BENIGN": "1"},
        ),
        Run(
            name="v37_hk_bp7_walker",
            group="hook",
            knob="BP7 as Walker 2023 describes it",
            value="SpliceAI at most 0.1, no conservation test",
            shipped="SpliceAI below 0.2 and PhyloP below 2.0",
            description="BP7 requires SpliceAI at most 0.1 and no longer tests conservation (the PhyloP bar set to 1000)",
            toml=["phylop_conserved = 1000.0"],
            env={"FASTVEP_BP7_SPLICEAI_BENIGN": "1"},
        ),
        Run(
            name="v37_hk_strict",
            group="hook",
            knob="the adjustments that read the variant's own ClinVar record or followed the reviewer's rows, removed together",
            value="BP1's veto leaves out the own record; no low-penetrance gate; no frequency-exception list",
            shipped="all three in place",
            description=(
                "BP1's veto no longer counts the variant's own ClinVar record, ClinVar's low-penetrance and risk-allele terms "
                "no longer withhold BS1 and BS2, and the curated frequency-exception list is emptied (each of the three was "
                "measured alone; none had been run with the others)"
            ),
            toml=["clinvar_low_penetrance_blocks_benign_frequency = false", "ba1_exceptions = []"],
            env={"FASTVEP_BP1_LOO": "1"},
        ),
        # The stricter BP1 count of Section 3.3 and the own-record count of Section 3.4 are proposed together
        # (Section 4.5), and the two missed diagnoses that a count of 1 avoids owe the veto to their own record,
        # so the pair is measured and not added up from the two runs that change one each.
        *[
            Run(
                name=f"v37_hk_bp1_loo_veto{count}",
                group="hook",
                knob="BP1 veto count, own ClinVar record excluded",
                value=str(count),
                shipped="3, own record counted",
                description=(
                    f"BP1's veto at {count} pathogenic missense variant{'s' if count > 1 else ''} in the gene, and no longer "
                    "counting the variant's own record when it is itself pathogenic or likely pathogenic"
                ),
                toml=[f"bp1_max_pathogenic_missense = {count}"],
                env={"FASTVEP_BP1_LOO": "1"},
            )
            for count in (1, 2)
        ],
    ]


def gap_inputs() -> dict[str, Path]:
    """Files the gap runs read, all under data/benchmark and none of them in the repository."""
    return {
        "complete_oga": BENCH / "sa_replace" / "omim_complete.oga",
        "orphanet_table": ROOT / "analysis/acmg_benchmark/data/gene_disease_attributes.tsv",
        "config_rs": ROOT / "crates/fastvep-classification/src/config.rs",
    }


def orphanet_mechanisms() -> dict[str, str]:
    """Genes whose disease mechanism Orphanet labels, and that the shipped table does not already name.

    A gene with a loss-of-function association and no gain-of-function one is `LOF`; with both it is
    `LOF_and_GOF`, which the classifier already reads as "a null allele is still pathogenic for the
    first disease". Orphanet's label is per gene-disorder association and is independent of ClinVar.
    """
    import re

    attrs = gap_inputs()["orphanet_table"]
    src = gap_inputs()["config_rs"].read_text()
    body = src[src.index("fn default_gene_mechanisms") :]
    body = body[: body.index("collect()")]
    shipped = set(re.findall(r'\(\s*"([A-Z0-9-]+)"\s*,\s*"[A-Za-z_]+"\s*\)', body))
    kinds: dict[str, set[str]] = {}
    with attrs.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            k = kinds.setdefault(row["gene"], set())
            if "(loss of function)" in row["association"]:
                k.add("LOF")
            if "(gain of function)" in row["association"]:
                k.add("GOF")
    out = {}
    for gene, k in sorted(kinds.items()):
        if gene in shipped or "LOF" not in k:
            continue
        out[gene] = "LOF_and_GOF" if "GOF" in k else "LOF"
    return out


def gaps() -> list[Run]:
    mech = orphanet_mechanisms()
    toml: list[str] = []
    for gene, m in mech.items():
        toml += [f"[gene_overrides.{gene}]", f'mechanism = "{m}"']
    complete = {"omim.oga": str(gap_inputs()["complete_oga"])}
    return [
        Run(
            name="v37_gap_clingen_complete",
            group="gap",
            knob="ClinGen gene-disease file",
            value="the 30 April export plus the curations it lacked",
            shipped="the 30 April export",
            description=(
                "the ClinGen file with the 112 rows put back that the 11 October 2026 export dates on or before 30 April "
                "and whose gene and disease pair the export of that day did not list (104 curations, CFTR among them, "
                "7 re-labelled rows and 1 row dated that day)"
            ),
            sa_replace=complete,
        ),
        Run(
            name="v37_gap_orphanet_lof",
            group="gap",
            knob="PVS1 licensed by Orphanet's loss-of-function label",
            value=f"{len(mech)} genes",
            shipped="none",
            description=(
                f"PVS1 licensed for {len(mech)} genes that Orphanet labels as a loss-of-function disease "
                "and that the shipped mechanism table does not name"
            ),
            toml=toml,
        ),
        Run(
            name="v37_gap_both",
            group="gap",
            knob="the complete ClinGen file and Orphanet's label together",
            value=f"{len(mech)} genes",
            shipped="none",
            description="both of the above together",
            toml=toml,
            sa_replace=complete,
        ),
    ]


def holdout_runs() -> list[Run]:
    """The changes that Section 4.5 proposes, on the variants of the 4 October 2026 release.

    Each differs from `v37_oct2026` (the shipped classifier on that release) in one change that
    the manuscript measured on the April set. The release's variants that were not in the April
    set were not seen by any sweep, review round or run, so these are the nearest thing to a test
    of a change that was found by reading the April discordances.
    """
    subset = BENCH / "clinvar_20261004_covered"
    inputs = {"INPUT": str(subset / "clinvar_2star.vcf"), "LABELS": str(subset / "clinvar_2star_labels.tsv")}
    mech = orphanet_mechanisms()
    toml: list[str] = []
    for gene, m in mech.items():
        toml += [f"[gene_overrides.{gene}]", f'mechanism = "{m}"']
    complete = {"omim.oga": str(gap_inputs()["complete_oga"])}

    def hold(short, what, **kw):
        return Run(name=f"v37_oct_{short}", group="holdout", knob=what, value="-", shipped="-",
                   description=f"{what}, on the variants of the 4 October 2026 release", inputs=inputs, **kw)

    return [
        hold("clingen_complete", "the ClinGen file completed with the rows its 30 April export lacked", sa_replace=complete),
        hold("orphanet_lof", f"PVS1 licensed for {len(mech)} genes that Orphanet labels as a loss-of-function disease", toml=toml),
        hold("both", "the completed ClinGen file and Orphanet's label together", toml=toml, sa_replace=complete),
        hold("bp7_walker", "BP7 with SpliceAI at most 0.1 and no conservation test", toml=["phylop_conserved = 1000.0"],
             env={"FASTVEP_BP7_SPLICEAI_BENIGN": "1"}, binary="hook"),
        hold("bp1_loo", "BP1's veto leaving out the variant's own ClinVar record", env={"FASTVEP_BP1_LOO": "1"}, binary="hook"),
        hold("no_an_floor", "BA1 and BS1 without the allele-number floor of 2,000", toml=["min_an_for_frequency_criteria = 0"]),
        # Not run with the own record left out as well: the hook discounts one record whenever ClinVar's allele-level
        # record (27 June) is pathogenic, which for a new variant that the April index does not hold as pathogenic
        # discounts a record that is not in its count.
        hold("bp1_veto1", "BP1's veto at a count of 1 pathogenic missense variant in the gene", toml=["bp1_max_pathogenic_missense = 1"]),
    ]


def all_runs() -> list[Run]:
    return sweeps() + switches() + sources() + interactions() + cutoffs() + hooks() + gaps()


def write_config(run: Run) -> Path | None:
    if not run.toml:
        return None
    path = RES / "configs" / f"{run.name}.toml"
    path.parent.mkdir(parents=True, exist_ok=True)
    body = [f"# {run.name}: {run.description}.", "# Everything else is the shipped default."]
    body += run.toml
    path.write_text("\n".join(body) + "\n")
    return path


def make_sa_dir(run: Run) -> Path:
    """A directory of symlinks to every sa_db file except the ones this run removes."""
    import fnmatch

    target = SA_ABLATE / run.name
    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True)
    removed = 0
    replaced = set()
    for f in sorted(SA_DB.iterdir()):
        if any(fnmatch.fnmatch(f.name, pat) for pat in run.sa_remove):
            removed += 1
            continue
        if f.name in run.sa_replace:
            replacement = Path(run.sa_replace[f.name])
            if not replacement.exists():
                raise SystemExit(f"{run.name}: {replacement} does not exist (see the module docstring)")
            (target / f.name).symlink_to(replacement)
            replaced.add(f.name)
            continue
        (target / f.name).symlink_to(f)
    if run.sa_remove and removed == 0:
        raise SystemExit(f"{run.name}: nothing in {SA_DB} matched {run.sa_remove}")
    if set(run.sa_replace) != replaced:
        raise SystemExit(f"{run.name}: {sorted(set(run.sa_replace) - replaced)} is not in {SA_DB}")
    return target


def read_matrix(path: Path) -> dict[str, dict[str, int]]:
    rows = list(csv.reader(path.open()))
    header = rows[0][1:]
    return {r[0]: {c: int(v) for c, v in zip(header, r[1:])} for r in rows[1:]}


def read_firing(path: Path) -> dict[str, dict[str, int]]:
    out: dict[str, dict[str, int]] = {}
    with path.open() as f:
        for row in csv.DictReader(f):
            out[row["criterion"]] = {c: int(row[f"{c}_fired"]) for c in CLASSES}
    return out


def collect(runs: list[Run]) -> None:
    """Rebuild runs.tsv (the April set) and holdout_runs.tsv (the October release) from the scored files."""
    RES.mkdir(parents=True, exist_ok=True)
    (RES / "provenance").mkdir(exist_ok=True)
    scored = []
    for run in runs:
        out = OUT_ROOT / run.name
        m = out / "concordance_matrix.csv"
        if not m.exists():
            continue
        prov = out / "run_provenance.txt"
        if prov.exists():
            shutil.copyfile(prov, RES / "provenance" / f"{run.name}.txt")
        scored.append((run, read_matrix(m), read_firing(out / "criterion_firing_rates.csv")))

    def write_runs(path: Path, rows) -> None:
        crit_cols = sorted({c for _, _, fire in rows for c in fire})
        with path.open("w") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(
                ["name", "group", "knob", "value", "shipped"]
                + [f"cm_{t}__{p}" for t in CLASSES for p in COLS]
                + [f"fire_{c}" for c in crit_cols]
            )
            for run, mat, fire in rows:
                w.writerow(
                    [run.name, run.group, run.knob, run.value, run.shipped]
                    + [mat[t][p] for t in CLASSES for p in COLS]
                    + [sum(fire.get(c, {}).values()) for c in crit_cols]
                )

    write_runs(RES / "runs.tsv", [r for r in scored if r[0].group != "holdout"])
    write_runs(RES / "holdout_runs.tsv", [r for r in scored if r[0].group == "holdout"])
    with (RES / "run_matrix.tsv").open("w") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["name", "group", "knob", "value", "shipped", "description", "removed_from_sa_db",
                    "replaced_in_sa_db", "environment"])
        for run in runs:
            w.writerow(
                [run.name, run.group, run.knob, run.value, run.shipped, run.description,
                 ",".join(run.sa_remove), ",".join(sorted(run.sa_replace)),
                 " ".join(f"{k}={v}" for k, v in sorted(run.env.items()))]
            )
    print(f"collected {len([r for r in scored if r[0].group != 'holdout'])} of {len([r for r in runs if r.group != 'holdout'])} runs "
          f"into {RES / 'runs.tsv'} and {len([r for r in scored if r[0].group == 'holdout'])} of "
          f"{len([r for r in runs if r.group == 'holdout'])} held-out runs into {RES / 'holdout_runs.tsv'}")


_print_lock = threading.Lock()


def execute(run: Run) -> tuple[str, bool]:
    out = OUT_ROOT / run.name
    if (out / "concordance_matrix.csv").exists():
        return run.name, True
    env = os.environ.copy()
    patched = {"cutoff": "FASTVEP_CUTOFF_BIN", "hook": "FASTVEP_HOOKS_BIN"}
    kind = run.binary or run.group
    if kind in patched:
        binary = env.get(patched[kind], "")
        if not binary:
            with _print_lock:
                print(f"skip {run.name}: {patched[kind]} is not set (see the module docstring)", flush=True)
            return run.name, True
        env["FASTVEP"] = binary
        env.update(run.env)
    env.update(run.inputs)
    out.mkdir(parents=True, exist_ok=True)
    env["OUT_DIR"] = str(out)
    cfg = write_config(run)
    if cfg:
        env["ACMG_CONFIG"] = str(cfg)
    if run.sa_remove or run.sa_replace:
        env["SA_DIR"] = str(make_sa_dir(run))
    log = out / "run.log"
    with log.open("w") as lf:
        rc = subprocess.call(
            ["bash", str(RUNNER), run.name], env=env, stdout=lf, stderr=subprocess.STDOUT, cwd=ROOT
        )
    ok = rc == 0 and (out / "concordance_matrix.csv").exists()
    if ok and run.name not in KEEP_VCF and run.group != "holdout":
        # the VCFs of the held-out runs are read by 12_per_variant_table.py, which is run on them by hand
        vcf = out / "clinvar_2star.fastvep.vcf.gz"
        if vcf.exists():
            vcf.unlink()
    if ok and (run.sa_remove or run.sa_replace):
        shutil.rmtree(SA_ABLATE / run.name, ignore_errors=True)
    with _print_lock:
        print(f"{'ok  ' if ok else 'FAIL'} {run.name} (exit {rc})", flush=True)
    return run.name, ok


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--only", default="", help="run only runs whose name contains this text")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--collect", action="store_true")
    ap.add_argument("--holdout", action="store_true", help="run the held-out group (on the 4 October 2026 release) and nothing else")
    args = ap.parse_args()

    runs = all_runs() + holdout_runs()
    names = [r.name for r in runs]
    if len(set(names)) != len(names):
        sys.exit("duplicate run names in the matrix")
    if args.list:
        for r in runs:
            print(f"{r.group:7} {r.name:34} {r.description}")
        print(f"{len(runs)} runs")
        return
    if args.collect:
        collect(runs)
        return

    todo = [r for r in runs if args.only in r.name and (r.group == "holdout") == args.holdout]
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    failed = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for name, ok in pool.map(execute, todo):
            if not ok:
                failed.append(name)
    collect(runs)
    if failed:
        sys.exit(f"failed: {', '.join(failed)}")


if __name__ == "__main__":
    main()
