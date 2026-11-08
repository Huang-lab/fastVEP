# Runs made for the ACMG manuscript

The runs here are of two families, each a list of 96 runs made by `run_benchmark.sh` and scored by `03_evaluate_concordance.py`: the family v38, the latest classifier (release 0.5.0, the tag `v0.5.0-rc.1`), which is the classifier that the manuscript describes, and the family v37, the classifier at the close of the review rounds (commit `5d8004a`), which is the before of what v38 changes (Section 3.8).
Each run differs from its family's main run (`v38` or `v37`) in one setting, one removed evidence source, one point cutoff, one patched setting of a measurement build, or one input, and the two families are the same list, so that a pair `v37_x` and `v38_x` differs in the classifier and in nothing else.
`scripts/11_manuscript_runs.py` defines them (`--family v38`, the default, and `--family v37`) and writes this directory, which keeps the name it had when it held the runs of v37 alone.
The first section below is the family v37, whose files carry no suffix; the family v38 and the runs of Section 3.8 are the `_v38` files, described in the sections after it.

ClinVar's 2-star+ classifications are a reference and not ground truth.

## The family v37: the runs of the classifier at the close of the review rounds

| File | What it holds |
|---|---|
| `run_matrix.tsv` | One row per run: its group (sweep, switch, source, interaction, cutoff, hook, gap or holdout), the setting, the default value (the column `shipped`), a description, the files replaced in the annotation directory and the environment the run needs |
| `runs.tsv` | One row per run on the April set: the 5x6 concordance matrix and the firing count of every criterion, read back from the scored files |
| `holdout_runs.tsv` | The same for the seven runs on the variants of the 4 October 2026 release (the `holdout` group of `run_matrix.tsv`) |
| `configs/<run>.toml` | The `--acmg-config` of every run that needs one |
| `provenance/<run>.txt` | What produced each run's annotation: the flags, the SHA-256 of the binary and of the configuration, the inputs read and the commit |
| `binaries.tsv` | Every binary that made a run: its SHA-256, what it was built from, the patch, the environment variables it reads and the number of runs it made |
| `rerun_comparison.csv` | The comparison of each run and of each per-variant table of the repetition of 4 November 2026 with the first run it repeats: identical or not, the annotated VCFs compared by record, and what differs. The outputs of the first runs were deleted once this comparison showed that the repetition reproduces them, so this file is what remains of them |

89 runs on the April set: 41 sweeps, 16 single-switch ablations, 9 leave-one-source-out ablations, 4 PM2 sweeps under the 2015 table, 9 runs of the point cutoffs, 7 runs of the second measurement build (`hook`) and 3 runs that change an input (`gap`).
Seven more runs (`holdout`) apply one change each to the variants of the October 2026 release, and their matrices cover all 781,135 variants of that release; `13_manuscript_stats.py` scores them on the 116,078 new ones from their per-variant tables.

Things to know when reading the provenance files.

- The runs were made on 8 November 2026, with the complete databases, by three binaries built from clean exports of commit `5d8004a`, the merge that made `v37`, each in a target directory of its own.
  78 runs used the plain build (SHA-256 `1b8e7c46...`).
  The 9 cutoff runs (`v37_cut_*`) used the build with `manuscript/acmg/replay_validation/combiner_cutoff_env.patch` (`d79ad4f5...`), which differs from the code of v37 in eight lines (six added, two removed).
  9 runs used the build with `bp_hooks_env.patch` (`823086d3...`), which reads BP1's own-record count and BP7's SpliceAI threshold from the environment: the 7 `hook` runs and the two held-out runs that change those settings.
  `binaries.tsv` lists every binary that made a run, with what it was built from and how many runs it made.
- Each of these runs was made once on 8 November 2026, with the complete databases, and none is a repetition of the others.
  The same runs had been made twice before with the first databases (on 6, 11 and 22 October 2026 by binaries built at earlier commits of the branch, and again on 4 November 2026 by the binaries above), and those two sets agreed in every scored file of every run and in every record of the annotated VCFs that both kept (`rerun_comparison.csv` is that comparison and was not made again; `RUN_VERSIONS.md`, "The full rerun of 4 November 2026").
  The runs of 8 November are not a third repetition of them, because the databases differ.
  The scored CSV files of the runs of 6 October ended their lines with CR LF, and the scoring script now writes LF.
- The run at the default cutoff (`v37_cut_lp6`) is a control of the cutoff build and reproduces `v37` in every cell.
  `v37_hk_control`, with nothing set, is the control of the hook build and reproduces `v37` in every cell and in the firing count of every criterion.
- The 7 `holdout` runs used the plain build, except `v37_oct_bp7_walker` and `v37_oct_bp1_loo`, which used the hook build; the three that change the ClinGen file or Orphanet's label need the files of the `gap` runs.
- The 3 `gap` runs need files that are not in the repository: the ClinGen export with its omitted curations put back (`scripts/15_clingen_complete_export.py` builds it) and, for two of them, a configuration of 1,042 gene overrides that `11_manuscript_runs.py` writes from `analysis/acmg_benchmark/data/gene_disease_attributes.tsv`.
- `repo_commit` is the commit of the checkout that ran the harness, which moved during 8 November 2026 as the branch received commits (`3d8c331` to `54c9477`), and it is not the commit the binary was built from: `binaries.tsv` says that.
  `repo_dirty_files` counts the files that were modified or untracked when a run was made, which include the results of the runs themselves, so it differs from run to run and says nothing about the classifier of `v37`, which was built from a clean export.
  The `fastvep` line is the path of the binary when the run was made, in a temporary directory.
- The environment variable that sets a run of a patched build (`FASTVEP_LP_CUT`, `FASTVEP_P_CUT`, `FASTVEP_BP1_LOO` or `FASTVEP_BP7_SPLICEAI_BENIGN`) is not part of what the harness compares and is not in the provenance file.
  It is in the `environment` column of `run_matrix.tsv`, which is the only record of it.

## The runs of Section 3.8 (v38)

42 more runs, made after `v37`, measure the change that reads the mode of inheritance ClinGen curates (Section 3.8 of the manuscript).
Their files carry the suffix `_v38`, and the files above are `v37`'s and were not touched.

| File | What it holds |
|---|---|
| `run_matrix_v38.tsv` | One row per run, as `run_matrix.tsv` has: the group, the setting, the default value (the column `shipped`), the description and the environment the run needs |
| `runs_v38.tsv` | The 38 runs on the April set: the 5x6 concordance matrix and the firing count of every criterion |
| `holdout_runs_v38.tsv` | The 4 runs on the variants of the 4 October 2026 release |
| `configs/v38_*.toml` | The `--acmg-config` of each run that needs one |
| `provenance/v38_*.txt` | What produced each run's annotation |

The 38 runs on the April set are the commit before the change (`v38_base`), the adopted change (`v38`), the control of the measurement build with nothing set (`v38_alt_control`) and four alternatives (`v38_alt_unguarded`, `v38_alt_label_guard`, `v38_alt_bs2_as_v37`, `v38_alt_triple_cap`), the adopted change with gnomAD's region flags withholding nothing (`v38_noflags`), the 26 sweeps of the four thresholds that read inheritance (`v38_sw_*`) the three runs of the pilot on chromosomes 21 and 22 (`v38_chr2122_*`) and the joint release on every chromosome (`v38_joint`, made on 5 November 2026, once the databases of all 24 chromosomes existed).
The 4 on the October release are `v38_oct_base`, `v38_oct_shipped`, `v38_oct_alt_unguarded` and `v38_oct_alt_triple_cap`.
`v38_joint` is defined in `11_manuscript_runs.py` and is made only when the database of all 24 chromosomes exists; it was made on 8 November 2026 like the others.

Things to know when reading these provenance files.

- Three binaries made the runs, all on 8 November 2026, with the complete databases.
  `v38_base` and `v38_oct_base` used the binary of the commit before the change, `9d84be0` (SHA-256 `6d30419c...`), whose classification code is `v37`'s (the two differ in comments only) and whose annotation layer is newer: its run differs from `v37`'s in the criteria of 25 indels and in the call of 4 of them.
  7 runs, the alternatives and their control among them, used a build of `manuscript/acmg/replay_validation/inheritance_alternatives_env.patch` on the tag (`87ebbc9a...`).
  The other 33, the adopted change among them, used the plain build of the tag `v0.5.0-rc.1` (`70b7c684...`).
  `binaries_v38.tsv` says what each was built from.
- Each of the 42 was made a second time by the same binary, into `data/benchmark/manuscript_runs_v0.5.0/`, and the two sets agree in every scored file and in the 15 annotated VCFs that both kept (`../v38_v0.5.0/rerun_comparison.csv`).
  Every one of them had been made before with the first databases, by builds of the working tree that became `v38` (`f184277d...`, `a2bb20ee...`, `f6e7f4ca...`) and by builds of the tag, and those sets agreed (`RUN_VERSIONS.md`, "v38" and "Release 0.5.0 reproduces v38"); the runs of this directory are not repetitions of those, because the databases differ.
- `repo_commit` is the commit of the checkout that ran the harness (the commits of the branch of 8 November 2026) and not the one a binary was built from; `repo_dirty_files` counts the modified and untracked files, which include the results of the runs themselves, and so differs from run to run.
- The environment variable that sets a run of the alternatives build (`FASTVEP_ALT`, with the values `g1_off`, `g3`, `bs2_names,g1_off` and `triple_cap`) is not in the provenance file.
  It is in the `environment` column of `run_matrix_v38.tsv`, which is the only record of it.
- The runs of the pilot read the databases of gnomAD's exomes and joint release for chromosomes 21 and 22, and `v38_joint` those of the joint release for all 24, that `build_gnomad_release_per_chrom.sh` builds; they are not in the repository.

## The family v38: the runs of Sections 2 to 3.7 made again by release 0.5.0

The 96 runs of `v37`'s family (the 89 above on the April set and the 7 of the `holdout` group) were made again by the classifier of release 0.5.0 (the tag `v0.5.0-rc.1`, `c12d9f9`), under the same names with `v38` for `v37`, so that the manuscript, which says v38, can set the latest classifier against the version in which the review rounds ended, run for run (`RUN_VERSIONS.md`, "The family v38").
The files are the `_v38` files above, which now hold the runs of Section 3.8 and the family together: `run_matrix_v38.tsv` 120 rows, `runs_v38.tsv` 109 rows (the April set) and `holdout_runs_v38.tsv` 11 rows (the October release), with `provenance/v38_*.txt` and `configs/v38_*.toml` for each, and `binaries_v38.tsv` for the binaries.
`11_manuscript_runs.py --family v38` writes the 96; `--v38` still writes the 42 of Section 3.8, and the 18 sweeps that both definitions name (`v38_sw_pm2ad_*`, `v38_sw_pm2ar_*`, `v38_sw_bs2hom_*`) are the same runs, made once by the first and used by the second.
Each pair of runs, `v37_x` and `v38_x`, has the same knob, value, default, replaced annotation files and environment in `run_matrix.tsv` and `run_matrix_v38.tsv`; `analysis/acmg_benchmark/results/v38_family_vs_v37.md` sets each pair against its twin, and its four CSV files hold the numbers.

| File | What it holds |
|---|---|
| `binaries_v38.tsv` | The five binaries that made the runs of v38: `v38` (SHA-256 `70b7c684...`, 93 runs), `cut38` (`de44db94...`, 9 runs, the point cutoffs), `hook38` (`9886d1b5...`, 9 runs, BP1's own-record count and BP7's SpliceAI threshold), `base` (`6d30419c...`, 2 runs, the commit before the change) and `alt` (`87ebbc9a...`, 7 runs, the alternatives), each built from a clean export by `scripts/build_manuscript_binaries.sh` |

Things to know when reading the provenance files of the family v38.

- All 96 runs were made on 8 November 2026, with the complete databases, by the three binaries of the family (`v38`, `cut38` and `hook38`), built from the tag by `scripts/build_manuscript_binaries.sh` (`V38_REF=v0.5.0-rc.1 NAMES="cut38 hook38"` for the two patched builds); the patches apply to the tag as they are, the cutoff patch with no offset and the hooks patch with offsets of 2 and 20 lines.
  The 18 that Section 3.8 also defines were made once and are used by both.
- The controls (`v38_cut_lp6`, `v38_hk_control`) reproduce `v38` in every cell of the concordance matrix and in the firing count of every criterion, and in every scored file byte for byte.
- `repo_commit` is the commit of the checkout that ran the harness, which moved during the day (`3d8c331` to `54c9477`), and not the one the binary was built from.
- Each of the 96 was made once, and the 42 of Section 3.8 were made a second time (`../v38_v0.5.0/`).
- The annotated VCF is kept for the same 15 April runs and 7 October runs as the family v37 keeps (`KEEP_VCF`), and the per-variant tables are in `data/benchmark/manuscript_runs/`, which is not in the repository.
- The 382 ClinVar reference-agreement records (385 on the October release) have `.` as their alternative allele and no call under either classifier; the binaries of v37 write a record for them and the binaries of v38 do not, so the scored matrices of the family v38 hold 673,278 records (the harness counts the other 382 apart), and `13_manuscript_stats.py --base v38` and `20_compare_families.py` put them back as no-calls.
