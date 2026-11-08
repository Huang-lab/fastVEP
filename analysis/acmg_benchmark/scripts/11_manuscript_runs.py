#!/usr/bin/env python3
"""Run the threshold sweeps and leave-one-out ablations the ACMG manuscript reports.

Every run here is a variant of `v37`, the final classifier, made by
`run_benchmark.sh` so that it is scored by `03_evaluate_concordance.py` and records
what produced it. A variant of a run is a run of its own, so each one has its own
version name, `v37_<knob>_<value>`.

The runs are defined once, for a family: `v37` (the classifier of Sections 2 to 3.7, names `v37_<...>`) and `v38` (the
classifier of release 0.5.0, the tag `v0.5.0-rc.1`, names `v38_<...>`). The two families are made by the same functions
below with the name prefix, the binaries and the names of two controls as their only parameters, so that they can differ in
the classifier and in nothing else (`--family`). The family `v37` is what this script has always made. The family `v38` is
the same list made again by the later classifier, and it takes in the runs of Section 3.8 that are members of it
(`v38_sw_pm2ad_*`, `v38_sw_pm2ar_*` and `v38_sw_bs2hom_*`, which have the grids of the family): they are not made again,
and `main` stops if Section 3.8 defines one of them differently from the family. Section 3.8's own runs (`--v38`) are
unchanged, and so are their files, to which the runs of the family `v38` add rows.

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
  v38     what changed after v37 (Section 3.8): the classifier reading ClinGen's curated inheritance,
          classifying a gene against each of its modes of inheritance, and not counting BS2's carriers
          beside BS1 or BA1. Four kinds of run, kept apart from the rest (their names start `v38`, they
          are written to `runs_v38.tsv`, `holdout_runs_v38.tsv` and `run_matrix_v38.tsv`, and the
          other files are not touched by them), run with `--v38`:
            inherit       the shipped code; the sweeps of the four thresholds that read inheritance,
                          re-run now that the field is read (binary in FASTVEP_V38_BIN)
            inherit_base  the commit before the change, for the comparison (FASTVEP_V38_BASE_BIN)
            inherit_alt   the alternatives that were measured and not shipped, and one that was not
                          tried until now, by a build of `replay_validation/inheritance_alternatives_env.patch`
                          (FASTVEP_V38_ALT_BIN)
            joint         the benchmark's gnomAD exomes replaced by gnomAD's joint exomes and genomes
                          callset on every chromosome, once its 24 databases are built
                          (`build_gnomad_release_per_chrom.sh`), against the same variants with the exomes;
                          and the pilot of the same on chr21 and chr22 alone, made first and kept as a
                          check that the runs of every chromosome give its calls
          and each of them also on the variants of the October 2026 release (`--holdout --v38`).
  holdout the changes that Section 4.5 proposes, run on the variants of the ClinVar release of
          4 October 2026 (the covered subset that `14_temporal_holdout.py prepare` writes), so
          that they are tried on variants that none of them was tuned on. They are run with
          `--holdout`, and are not among the runs of `runs.tsv` (which are all on the April set):
          their matrices are in `holdout_runs.tsv`, and `13_manuscript_stats.py` scores them on
          the new variants from the per-variant tables of `12_per_variant_table.py`.

Outputs (MANUSCRIPT_RUNS_ROOT and MANUSCRIPT_RESULTS_DIR name another root and another results directory, for a rerun)
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
      run_matrix_v38.tsv, runs_v38.tsv, holdout_runs_v38.tsv
                           the same for every run whose name starts `v38`: Section 3.8's, then the family v38's

Usage
  11_manuscript_runs.py                 run everything that has not been run (the family v37)
  11_manuscript_runs.py --workers 2     two runs at a time
  11_manuscript_runs.py --only pm2ad    only runs whose name contains the text
  11_manuscript_runs.py --list          print the matrix and stop
  11_manuscript_runs.py --collect       rebuild runs.tsv from the scored files and stop
  11_manuscript_runs.py --holdout       the runs of the `holdout` group (the others are not run)
  11_manuscript_runs.py --family v38    the family v38 in place of v37 (with --list: print its matrix; with --holdout: its
                                        seven runs on the October release); the runs that Section 3.8 made are found
                                        scored and not made again, so a second invocation makes nothing
  11_manuscript_runs.py --v38           the runs of Section 3.8 (not the family v38)

  FASTVEP_V37_BIN=/path/to/v37/fastvep 11_manuscript_runs.py --workers 3
  FASTVEP_CUTOFF_BIN=/path/to/patched/fastvep 11_manuscript_runs.py --only cut
  export FASTVEP_V38_BIN=/path/to/fastvep_v38 FASTVEP_CUTOFF38_BIN=/path/to/fastvep_cut38 FASTVEP_HOOKS38_BIN=/path/to/fastvep_hook38
  11_manuscript_runs.py --family v38 --workers 3        (and again with --holdout)

The binary of every run is named by an environment variable, and a run whose variable is not set is not run:
FASTVEP_V37_BIN for the runs of the family v37 that are v37 with one setting changed, and FASTVEP_V38_BIN for those of
the family v38 that are v38 with one setting changed (the script stops if it is missing, because the checkout's own
build is the classifier of the checkout and not v37 or v38 once the classifier has moved on), and the variable of its
group for a patched build or a build of Section 3.8 (FASTVEP_CUTOFF_BIN and FASTVEP_HOOKS_BIN, built from v37, and
FASTVEP_CUTOFF38_BIN and FASTVEP_HOOKS38_BIN, built from v38, all four by `build_manuscript_binaries.sh`;
FASTVEP_V38_BIN, FASTVEP_V38_BASE_BIN, FASTVEP_V38_ALT_BIN for Section 3.8), which is skipped where it is not set.

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
# A rerun is made into a root and a results directory of its own (MANUSCRIPT_RUNS_ROOT, MANUSCRIPT_RESULTS_DIR), so that it
# neither finds the first run's directories already scored and skips them nor overwrites the tables the manuscript reads;
# `16_compare_runs.py` then compares the two roots.
OUT_ROOT = Path(os.environ.get("MANUSCRIPT_RUNS_ROOT", BENCH / "manuscript_runs"))
SA_DB = BENCH / "sa_db"
SA_ABLATE = BENCH / "sa_db_ablate"
RES = Path(os.environ.get("MANUSCRIPT_RESULTS_DIR", ROOT / "analysis/acmg_benchmark/results/v37_manuscript"))
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
    # Section 3.8: the comparisons that are read call by call
    "v38_base",
    "v38",
    "v38_alt_unguarded",
    "v38_alt_label_guard",
    "v38_alt_bs2_as_v37",
    "v38_alt_triple_cap",
    "v38_noflags",
    "v38_chr2122_exomes",
    "v38_chr2122_exomes_noflags",
    "v38_chr2122_joint",
    "v38_joint",
}
# The family v38 keeps the annotated VCF of every run whose v37 twin is kept, so that the same call-by-call analyses can be made of it.
KEEP_VCF |= {"v38_" + n.removeprefix("v37_") for n in KEEP_VCF if n.startswith("v37_")}

# The two families of runs (see the module docstring): the classifier of Sections 2 to 3.7 and the one of release 0.5.0.
FAMILIES = ("v37", "v38")


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


def sweeps(fam: str = "v37") -> list[Run]:
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
            f"number of ClinVar pathogenic missense variants in a gene at which BP1's veto applies (no veto at all is {fam}_sw_no_bp1_veto)",
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
                    name=f"{fam}_sw_{short}_{label}",
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
                name=f"{fam}_sw_bp7far_{v}",
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


def switches(fam: str = "v37") -> list[Run]:
    def sw(short, knob, value, shipped, what, toml=None):
        return Run(
            name=f"{fam}_sw_{short}",
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


def sources(fam: str = "v37") -> list[Run]:
    def src(short, remove, what, toml=None):
        return Run(
            name=f"{fam}_src_{short}",
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


def interactions(fam: str = "v37") -> list[Run]:
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
                name=f"{fam}_ix_table_pm2ad_{label}",
                group="interaction",
                knob="pm2_ad_af_threshold (Richards 2015 table both sides)",
                value=label,
                shipped="4e-05",
                description=f"PM2 bar for dominant genes set to {label}, with the Richards 2015 table on both sides",
                toml=["use_point_system = false", f"pm2_ad_af_threshold = {v!r}"],
            )
        )
    return out


def cutoffs(fam: str = "v37") -> list[Run]:
    """The point cutoffs of the pathogenic side, by real runs of a patched binary.

    Likely pathogenic from 6 points and Pathogenic from 10 are constants in `combiner.rs`. The
    patch makes them read `FASTVEP_LP_CUT` and `FASTVEP_P_CUT`, so that what a different cutoff
    does is measured and not only replayed from the criteria a variant met.
    """
    out = []
    # 6 is the shipped cutoff: that run is the control, which must reproduce the family's own run (v37 or v38) exactly,
    # because it shows that the patched binary behaves as the shipped one does apart from the two cutoffs.
    for lp in (4, 5, 6, 7, 8, 9):
        out.append(
            Run(
                name=f"{fam}_cut_lp{lp}",
                group="cutoff",
                knob="likely_pathogenic_points",
                value=str(lp),
                shipped="6",
                description=(
                    f"Likely pathogenic from {lp} pathogenic points"
                    + (f" (the shipped value: a control that must reproduce {fam})" if lp == 6 else " (6 is shipped)")
                ),
                env={"FASTVEP_LP_CUT": str(lp)},
            )
        )
    for p in (9, 11, 12):
        out.append(
            Run(
                name=f"{fam}_cut_p{p}",
                group="cutoff",
                knob="pathogenic_points",
                value=str(p),
                shipped="10",
                description=f"Pathogenic from {p} pathogenic points (10 is shipped)",
                env={"FASTVEP_P_CUT": str(p)},
            )
        )
    return out


def hooks(fam: str = "v37") -> list[Run]:
    """BP1 and BP7 measured by real runs of a build that exposes one setting of each.

    `<family>_hk_control` has nothing set, and must reproduce the family's own run (v37 or v38) in every cell: that is
    what licenses reading the other three as the shipped classifier with one thing changed.
    """
    return [
        Run(
            name=f"{fam}_hk_control",
            group="hook",
            knob="measurement build, nothing set",
            value="-",
            shipped="-",
            description=f"the measurement build with no hook set (a control that must reproduce {fam})",
        ),
        Run(
            name=f"{fam}_hk_bp1_loo",
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
            name=f"{fam}_hk_bp7_spliceai01",
            group="hook",
            knob="BP7 SpliceAI bar",
            value="0.1 or below",
            shipped="below 0.2",
            description="BP7 requires SpliceAI at most 0.1 (BP4's bar, as Walker 2023 requires) in place of below 0.2",
            env={"FASTVEP_BP7_SPLICEAI_BENIGN": "1"},
        ),
        Run(
            name=f"{fam}_hk_bp7_walker",
            group="hook",
            knob="BP7 as Walker 2023 describes it",
            value="SpliceAI at most 0.1, no conservation test",
            shipped="SpliceAI below 0.2 and PhyloP below 2.0",
            description="BP7 requires SpliceAI at most 0.1 and no longer tests conservation (the PhyloP bar set to 1000)",
            toml=["phylop_conserved = 1000.0"],
            env={"FASTVEP_BP7_SPLICEAI_BENIGN": "1"},
        ),
        Run(
            name=f"{fam}_hk_strict",
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
                name=f"{fam}_hk_bp1_loo_veto{count}",
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


def gaps(fam: str = "v37") -> list[Run]:
    mech = orphanet_mechanisms()
    toml: list[str] = []
    for gene, m in mech.items():
        toml += [f"[gene_overrides.{gene}]", f'mechanism = "{m}"']
    complete = {"omim.oga": str(gap_inputs()["complete_oga"])}
    return [
        Run(
            name=f"{fam}_gap_clingen_complete",
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
            name=f"{fam}_gap_orphanet_lof",
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
            name=f"{fam}_gap_both",
            group="gap",
            knob="the complete ClinGen file and Orphanet's label together",
            value=f"{len(mech)} genes",
            shipped="none",
            description="both of the above together",
            toml=toml,
            sa_replace=complete,
        ),
    ]


def holdout_runs(fam: str = "v37") -> list[Run]:
    """The changes that Section 4.5 proposes, on the variants of the 4 October 2026 release.

    Each differs from the family's run of the shipped classifier on that release (`v37_oct2026`, made by
    `14_temporal_holdout.py`, for v37; `v38_oct_shipped` for v38) in one change that
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
        return Run(name=f"{fam}_oct_{short}", group="holdout", knob=what, value="-", shipped="-",
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


V38_BINARIES = {"inherit": "FASTVEP_V38_BIN", "inherit_base": "FASTVEP_V38_BASE_BIN", "inherit_alt": "FASTVEP_V38_ALT_BIN"}
# The binary of every run that is not made by a patched build or a build of Section 3.8: the shipped classifier of v37,
# with one setting changed. It is named, and not taken from the checkout's `target/release`, because the checkout's
# shipped classifier is v37 only until the classifier moves on, and from then on a run called `v37_...` made by the
# checkout's build would be a run of something else under v37's name.
V37_BINARY = "FASTVEP_V37_BIN"
PATCHED_BINARIES = {"cutoff": "FASTVEP_CUTOFF_BIN", "hook": "FASTVEP_HOOKS_BIN", **V38_BINARIES}
# The same three kinds of binary for each family. The family v38's plain binary is the one that Section 3.8's `inherit`
# group is made by (the release's own, the tag v0.5.0-rc.1), and its two patched builds are the patches of the family v37
# applied to that tree (`build_manuscript_binaries.sh`, names `cut38` and `hook38`).
PLAIN_BINARY = {"v37": V37_BINARY, "v38": V38_BINARIES["inherit"]}
CUTOFF_BINARY = {"v37": "FASTVEP_CUTOFF_BIN", "v38": "FASTVEP_CUTOFF38_BIN"}
HOOK_BINARY = {"v37": "FASTVEP_HOOKS_BIN", "v38": "FASTVEP_HOOKS38_BIN"}
PILOT = BENCH / "pilot_chr2122"
SA_JOINT = BENCH / "sa_db_joint"
CHROMS = [str(c) for c in range(1, 23)] + ["X", "Y"]


def is_v38(run: Run) -> bool:
    return run.name.startswith("v38")


def family_of(run: Run) -> str:
    """The family a run belongs to: by its name, as `is_v38` says which of the output files it goes to."""
    return "v38" if is_v38(run) else "v37"


def binary_variable(run: Run) -> str:
    """The environment variable that names the binary this run is made by."""
    kind = run.binary or run.group
    if kind == "joint":
        kind = "inherit"
    fam = family_of(run)
    return {"cutoff": CUTOFF_BINARY[fam], "hook": HOOK_BINARY[fam], **V38_BINARIES}.get(kind, PLAIN_BINARY[fam])


def must_have_plain_binary(run: Run) -> bool:
    """Whether the script stops when the plain binary of the run's family is not named, in place of skipping the run.

    That is every run of a family that is its classifier with one setting changed. A patched build and the runs of
    Section 3.8 are skipped where their variable is not set, as they always have been.
    """
    kind = run.binary or run.group
    return kind not in ("cutoff", "hook", "joint", *V38_BINARIES)


def ensure_pilot_inputs() -> None:
    """The April benchmark's variants and labels on chromosomes 21 and 22, for the pilot of gnomAD's joint release.

    A filter of the benchmark's own two files with their headers, so that nothing about the subset is a choice.
    The pilot was made before the joint release had been extracted for the other chromosomes, where a variant would
    have read as absent from it, and it stays because the runs of every chromosome must give its calls on these
    variants (`13_manuscript_stats.py` checks that). Written whole each time, so that a file left over from a
    different benchmark cannot be the input.
    """
    PILOT.mkdir(parents=True, exist_ok=True)
    for name, header_prefix in (("clinvar_2star.vcf", "#"), ("clinvar_2star_labels.tsv", "chrom\t")):
        src, dst = BENCH / name, PILOT / name
        tmp = dst.with_suffix(dst.suffix + ".partial")
        with src.open() as f, tmp.open("w") as out:
            for line in f:
                if line.startswith(header_prefix) or line.split("\t", 1)[0] in ("21", "22"):
                    out.write(line)
        tmp.replace(dst)


def v38_runs() -> list[Run]:
    """What changed after v37, on the April set. See the module docstring.

    The alternatives are named for what each one does to the shipped code, not for the experiment:
    `unguarded` reads the tag and counts BS2's carriers beside BS1 and BA1 (the one-line change of
    reading the curated field and nothing else), `label_guard` also withholds BS2's carrier arm from a
    disorder that ClinGen's label calls a susceptibility or a predisposition, `bs2_as_v37` reads the
    tag for PM2, PM3 and BP2 and leaves BS2 as it was, and `triple_cap` caps PP2, PM1 and PP3 together
    at 4 points, which no version of the classifier does.
    """
    def run(name, group, what, **kw):
        return Run(name=name, group=group, knob=what, value="-", shipped="-", description=what, **kw)

    out = [
        run("v38_base", "inherit_base", "the commit before the change: inheritance read from the disease names"),
        run("v38", "inherit", "the shipped classifier: inheritance read from ClinGen's curated field, a gene classified against each "
            "of its modes of inheritance, BS2's carriers not counted beside BS1 or BA1"),
        run("v38_alt_control", "inherit_alt", "the alternatives build with nothing set (a control that must reproduce v38 in every cell)"),
        run("v38_alt_unguarded", "inherit_alt", "ClinGen's inheritance read, BS2's carriers counted beside BS1 and BA1",
            env={"FASTVEP_ALT": "g1_off"}),
        run("v38_alt_label_guard", "inherit_alt", "as shipped, and BS2's carrier arm withheld from a disorder that ClinGen's label calls "
            "a susceptibility or a predisposition", env={"FASTVEP_ALT": "g3"}),
        run("v38_alt_bs2_as_v37", "inherit_alt", "ClinGen's inheritance read for PM2, PM3 and BP2, and BS2 as before: dominant by "
            "disease name, carriers counted beside BS1 and BA1", env={"FASTVEP_ALT": "bs2_names,g1_off"}),
        run("v38_alt_triple_cap", "inherit_alt", "as shipped, with PP2, PM1 and PP3 capped together at 4 points",
            env={"FASTVEP_ALT": "triple_cap"}),
        run("v38_noflags", "inherit", "as shipped, with gnomAD's segmental-duplication and low-complexity flags no longer withholding "
            "the frequency criteria", toml=["gnomad_region_flags_block_frequency = false"]),
    ]

    # The four thresholds that read inheritance, re-swept now that the field is read. The grids are the ones
    # swept on v37, and the shipped value is `v38` itself.
    grids = [
        ("pm2_ad_af_threshold", "pm2ad", [0.0, 1e-6, 5e-6, 1e-5, 2e-5, 1e-4, 2e-4, 5e-4], "4e-05",
         "PM2 allele-frequency bar for genes not read as recessive-only"),
        ("pm2_ar_af_threshold", "pm2ar", [0.0, 1e-6, 1e-5, 2e-5, 2e-4, 5e-4], "7e-05",
         "PM2 allele-frequency bar for genes read as recessive-only"),
        ("bs2_hom_prevalence_threshold", "bs2hom", [1e-6, 1e-5, 1e-4, 1e-2], "0.001",
         "BS2 prevalence bar for the individuals with no functional copy"),
        ("bs2_ad_prevalence_threshold", "bs2ad", [0.0, 3e-5, 1e-4, 1e-3, 3e-3], "0.0003",
         "BS2 carrier bar for genes read as dominant-only (0 is no bar, so that the allele count of at least 5 alone decides)"),
    ]
    for knob, short, values, shipped, what in grids:
        for v in values:
            label = f"{v:g}"
            out.append(Run(name=f"v38_sw_{short}_{label}", group="inherit", knob=knob, value=label, shipped=shipped,
                           description=f"{what}, set to {label}", toml=[f"{knob} = {v!r}"]))
    for v in (20, 100, 1000):
        out.append(Run(name=f"v38_sw_bs2ac_{v}", group="inherit", knob="bs2_ad_min_ac (bar off)", value=str(v), shipped="5 (with bar 3e-4)",
                       description=f"BS2 carrier floor for genes read as dominant-only, set to {v} alleles with the carrier bar off",
                       toml=["bs2_ad_prevalence_threshold = 0.0", f"bs2_ad_min_ac = {v}"]))

    # The pilot of gnomAD's joint exomes and genomes against the exomes, on chr21 and chr22, which were the
    # two chromosomes extracted first, with the same variants and the same binary. The runs of every chromosome
    # below must give its calls (`13_manuscript_stats.py` checks that), so it stays as a check.
    inputs = {"INPUT": str(PILOT / "clinvar_2star.vcf"), "LABELS": str(PILOT / "clinvar_2star_labels.tsv")}
    pilot_joint = {f"gnomad_chr{c}.osa2": str(SA_JOINT / f"gnomad_chr{c}.osa2") for c in ("21", "22")}
    if all(Path(v).exists() for v in pilot_joint.values()) and (BENCH / "clinvar_2star.vcf").exists():
        out += [
            run("v38_chr2122_exomes", "joint", "chr21 and chr22 variants, gnomAD v4.1 exomes", inputs=inputs, binary="inherit"),
            run("v38_chr2122_exomes_noflags", "joint", "chr21 and chr22 variants, gnomAD v4.1 exomes, with its region flags no longer "
                "withholding the frequency criteria", inputs=inputs, binary="inherit", toml=["gnomad_region_flags_block_frequency = false"]),
            run("v38_chr2122_joint", "joint", "chr21 and chr22 variants, gnomAD v4.1 joint exomes and genomes", inputs=inputs,
                binary="inherit", sa_replace=pilot_joint),
        ]
    # and on every chromosome, where the joint database has been built for all of them
    full_joint = {f"gnomad_chr{c}.osa2": str(SA_JOINT / f"gnomad_chr{c}.osa2") for c in CHROMS}
    if all(Path(v).exists() for v in full_joint.values()):
        out.append(run("v38_joint", "joint", "gnomAD v4.1 joint exomes and genomes in place of the exomes", binary="inherit",
                       sa_replace=full_joint))
    return out


def v38_holdout_runs() -> list[Run]:
    """The same comparisons on the variants of the 4 October 2026 release (the covered subset)."""
    subset = BENCH / "clinvar_20261004_covered"
    inputs = {"INPUT": str(subset / "clinvar_2star.vcf"), "LABELS": str(subset / "clinvar_2star_labels.tsv")}

    def hold(short, group, what, **kw):
        return Run(name=f"v38_oct_{short}", group="holdout", knob=what, value="-", shipped="-",
                   description=f"{what}, on the variants of the 4 October 2026 release", inputs=inputs, binary=group, **kw)

    return [
        hold("base", "inherit_base", "the commit before the change"),
        hold("shipped", "inherit", "the shipped classifier"),
        hold("alt_unguarded", "inherit_alt", "ClinGen's inheritance read, BS2's carriers counted beside BS1 and BA1",
             env={"FASTVEP_ALT": "g1_off"}),
        hold("alt_triple_cap", "inherit_alt", "as shipped, with PP2, PM1 and PP3 capped together at 4 points",
             env={"FASTVEP_ALT": "triple_cap"}),
    ]


def all_runs(fam: str = "v37") -> list[Run]:
    return sweeps(fam) + switches(fam) + sources(fam) + interactions(fam) + cutoffs(fam) + hooks(fam) + gaps(fam)


def v38_family(section38: list[Run]) -> tuple[list[Run], list[Run]]:
    """The family v38 (what `all_runs("v38") + holdout_runs("v38")` define) and, of it, the runs that Section 3.8 has not made.

    Section 3.8 swept three of the four thresholds that read inheritance on the grids of the family
    (`v38_sw_pm2ad_*`, `v38_sw_pm2ar_*` and `v38_sw_bs2hom_*`), so those are members of the family under the same names and
    are not made again. They keep Section 3.8's definition (its group, `inherit`, and its description). A member is only
    the same run if the setting is the same, so the two definitions are compared and the script stops at a difference:
    a run that was made with another grid would otherwise stand in the family's sweep under its name.

    Returns the whole family, in the family's order and with Section 3.8's definition for the shared runs, and the runs
    of it that the matrix does not hold yet.
    """
    have = {r.name: r for r in section38}
    fam, new = [], []
    for r in all_runs("v38") + holdout_runs("v38"):
        old = have.get(r.name)
        if old is None:
            fam.append(r)
            new.append(r)
            continue
        key = lambda x: (x.knob, x.value, x.shipped, x.toml, x.sa_remove, x.sa_replace, x.env, x.inputs, x.binary)  # noqa: E731
        if key(old) != key(r):
            sys.exit(f"{r.name}: Section 3.8 defines it as {key(old)} and the family as {key(r)}, so it is not a member of the family")
        fam.append(old)
    return fam, new


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

    # The runs of Section 3.8 have files of their own, so that the tables of the other sections, which read
    # `runs.tsv` and `holdout_runs.tsv` whole, do not change when they are made.
    apr = [r for r in scored if r[0].group != "holdout" and not is_v38(r[0])]
    hold = [r for r in scored if r[0].group == "holdout" and not is_v38(r[0])]
    write_runs(RES / "runs.tsv", apr)
    write_runs(RES / "holdout_runs.tsv", hold)
    if any(is_v38(r[0]) for r in scored):
        write_runs(RES / "runs_v38.tsv", [r for r in scored if is_v38(r[0]) and r[0].group != "holdout"])
        write_runs(RES / "holdout_runs_v38.tsv", [r for r in scored if is_v38(r[0]) and r[0].group == "holdout"])
        with (RES / "run_matrix_v38.tsv").open("w") as f38:
            w38 = csv.writer(f38, delimiter="\t", lineterminator="\n")
            w38.writerow(["name", "group", "knob", "value", "shipped", "description", "removed_from_sa_db", "replaced_in_sa_db", "environment"])
            for run in runs:
                if is_v38(run):
                    w38.writerow([run.name, run.group, run.knob, run.value, run.shipped, run.description, ",".join(run.sa_remove),
                                  ",".join(sorted(run.sa_replace)), " ".join(f"{k}={v}" for k, v in sorted(run.env.items()))])
    runs = [r for r in runs if not is_v38(r)]
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
    print(f"collected {len(apr)} of {len([r for r in runs if r.group != 'holdout'])} runs "
          f"into {RES / 'runs.tsv'} and {len(hold)} of "
          f"{len([r for r in runs if r.group == 'holdout'])} held-out runs into {RES / 'holdout_runs.tsv'}"
          + (f"; {len([r for r in scored if is_v38(r[0])])} runs of Section 3.8 into the files with `v38` in their names"
             if any(is_v38(r[0]) for r in scored) else ""))


_print_lock = threading.Lock()


def execute(run: Run) -> tuple[str, bool]:
    out = OUT_ROOT / run.name
    if (out / "concordance_matrix.csv").exists():
        return run.name, True
    env = os.environ.copy()
    variable = binary_variable(run)
    binary = env.get(variable, "")
    if not binary:
        # `main` has already stopped the run of a script that lacks the binary of v37, so this is a patched build.
        with _print_lock:
            print(f"skip {run.name}: {variable} is not set (see the module docstring)", flush=True)
        return run.name, True
    if not Path(binary).is_file():
        with _print_lock:
            print(f"FAIL {run.name}: {variable} names {binary}, which is not a file", flush=True)
        return run.name, False
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
    ap.add_argument("--v38", action="store_true", help="run the runs of Section 3.8 (what changed after v37) in place of the others")
    ap.add_argument("--family", choices=FAMILIES, default="v37",
                    help="the family of runs to list and to run: v37, the runs of Sections 2 to 3.7 (the default, which lists and runs "
                         "what this script has always listed and run), or v38, the same runs made by the classifier of release 0.5.0")
    args = ap.parse_args()

    section38 = v38_runs() + v38_holdout_runs()
    family38, new38 = v38_family(section38)
    # Every run, in the order that decides the order of the rows of the files: the family v37, Section 3.8, then the rest of the family v38.
    runs = all_runs() + holdout_runs() + section38 + new38
    names = [r.name for r in runs]
    if len(set(names)) != len(names):
        sys.exit("duplicate run names in the matrix")
    if args.list:
        shown = runs[: len(runs) - len(new38)] if args.family == "v37" else family38
        for r in shown:
            print(f"{r.group:7} {r.name:34} {r.description}")
        print(f"{len(shown)} runs")
        return
    if args.collect:
        collect(runs)
        return

    if args.v38:
        pool = section38
    elif args.family == "v38":
        pool = family38
    else:
        pool = [r for r in runs if not is_v38(r)]
    todo = [r for r in pool if args.only in r.name and (r.group == "holdout") == args.holdout]
    for fam in FAMILIES:
        variable = PLAIN_BINARY[fam]
        if not os.environ.get(variable) and any(
            family_of(r) == fam and must_have_plain_binary(r) and not (OUT_ROOT / r.name / "concordance_matrix.csv").exists() for r in todo
        ):
            what = ("built at the v37 commit (5d8004a)" if fam == "v37" else
                    "of release 0.5.0 (built from the tag v0.5.0-rc.1)")
            sys.exit(f"{variable} is not set: it names the binary {what}, which makes every run of the family {fam} that is not made "
                     "by a patched build or a build of Section 3.8 (see the module docstring)")
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    if any(str(PILOT) in v for r in todo for v in r.inputs.values()):
        ensure_pilot_inputs()
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
