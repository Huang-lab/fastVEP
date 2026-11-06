# Runs made for the ACMG manuscript

Every run here is a variant of `v37`, the classifier of Sections 2 to 3.7 of the manuscript, made by `run_benchmark.sh` and scored by `03_evaluate_concordance.py`.
Each differs from `v37` in one setting, one removed evidence source, one point cutoff, one patched setting of a measurement build, or one input.
`scripts/11_manuscript_runs.py` defines them and writes this directory.

ClinVar's 2-star+ classifications are a reference and not ground truth.

| File | What it holds |
|---|---|
| `run_matrix.tsv` | One row per run: its group (sweep, switch, source, interaction, cutoff, hook, gap or holdout), the setting, the shipped value, a description, the files replaced in the annotation directory and the environment the run needs |
| `runs.tsv` | One row per run on the April set: the 5x6 concordance matrix and the firing count of every criterion, read back from the scored files |
| `holdout_runs.tsv` | The same for the seven runs on the variants of the 4 October 2026 release (the `holdout` group of `run_matrix.tsv`) |
| `configs/<run>.toml` | The `--acmg-config` of every run that needs one |
| `provenance/<run>.txt` | What produced each run's annotation: the flags, the SHA-256 of the binary and of the configuration, the inputs read and the commit |
| `binaries.tsv` | Every binary that made a run: its SHA-256, what it was built from, the patch, the environment variables it reads and the number of runs it made |
| `rerun_comparison.csv` | The comparison of each run and of each per-variant table of the repetition of 4 November 2026 with the first run it repeats: identical or not, the annotated VCFs compared by record, and what differs. The outputs of the first runs were deleted once this comparison showed that the repetition reproduces them, so this file is what remains of them |

89 runs on the April set: 41 sweeps, 16 single-switch ablations, 9 leave-one-source-out ablations, 4 PM2 sweeps under the 2015 table, 9 runs of the point cutoffs, 7 runs of the second measurement build (`hook`) and 3 runs that change an input (`gap`).
Seven more runs (`holdout`) apply one change each to the variants of the October 2026 release, and their matrices cover all 772,585 variants of that release's covered subset; `13_manuscript_stats.py` scores them on the 107,528 new ones from their per-variant tables.

Things to know when reading the provenance files.

- The runs were made on 4 November 2026 by three binaries built from clean exports of commit `5d8004a`, the merge that made `v37`, each in a target directory of its own.
  78 runs used the plain build (SHA-256 `1b8e7c46...`).
  The 9 cutoff runs (`v37_cut_*`) used the build with `manuscript/acmg/replay_validation/combiner_cutoff_env.patch` (`d79ad4f5...`), which differs from the shipped code in eight lines (six added, two removed).
  9 runs used the build with `bp_hooks_env.patch` (`823086d3...`), which reads BP1's own-record count and BP7's SpliceAI bar from the environment: the 7 `hook` runs and the two held-out runs that change those settings.
  `binaries.tsv` lists every binary that made a run, with what it was built from and how many runs it made.
- Each of these runs is a repetition.
  The same runs were first made on 6, 11 and 22 October 2026 with binaries built at earlier commits of the same branch (the plain build at `5d2521f`, SHA-256 `d12c0d8e...`, and the two patched builds, `dc275d1c...` and `2bccb937...`), and the repetition agrees with them in every scored file of every run and in every record of the annotated VCFs that both keep (`rerun_comparison.csv`; `RUN_VERSIONS.md`, "The full rerun of 4 November 2026").
  The scored CSV files of the runs of 6 October end their lines with CR LF, and the scoring script now writes LF.
- The run at the shipped cutoff (`v37_cut_lp6`) is a control of the cutoff build and reproduces `v37` in every cell.
  `v37_hk_control`, with nothing set, is the control of the hook build and reproduces `v37` in every cell and in the firing count of every criterion.
- The 7 `holdout` runs used the plain build, except `v37_oct_bp7_walker` and `v37_oct_bp1_loo`, which used the hook build; the three that change the ClinGen file or Orphanet's label need the files of the `gap` runs.
- The 3 `gap` runs need files that are not in the repository: the ClinGen export with its omitted curations put back (`scripts/15_clingen_complete_export.py` builds it) and, for two of them, a configuration of 1,042 gene overrides that `11_manuscript_runs.py` writes from `analysis/acmg_benchmark/data/gene_disease_attributes.tsv`.
- `repo_commit` is the commit of the checkout that ran the harness, which is `9d84be0` for every run because the change of Section 3.8 was not committed when they were made, and it is not the commit the binary was built from: `binaries.tsv` says that.
  `repo_dirty_files` counts the files that were modified or untracked when a run was made, which include that change and the results of the runs themselves, so it differs from run to run and says nothing about the classifier of `v37`, which was built from a clean export.
  The `fastvep` line is the path of the binary when the run was made, in a temporary directory.
- The environment variable that sets a run of a patched build (`FASTVEP_LP_CUT`, `FASTVEP_P_CUT`, `FASTVEP_BP1_LOO` or `FASTVEP_BP7_SPLICEAI_BENIGN`) is not part of what the harness compares and is not in the provenance file.
  It is in the `environment` column of `run_matrix.tsv`, which is the only record of it.

## The runs of Section 3.8 (v38)

41 more runs, made after `v37`, measure the change that reads the mode of inheritance ClinGen curates (Section 3.8 of the manuscript).
Their files carry the suffix `_v38`, and the files above are `v37`'s and were not touched.

| File | What it holds |
|---|---|
| `run_matrix_v38.tsv` | One row per run, as `run_matrix.tsv` has: the group, the setting, the shipped value, the description and the environment the run needs |
| `runs_v38.tsv` | The 38 runs on the April set: the 5x6 concordance matrix and the firing count of every criterion |
| `holdout_runs_v38.tsv` | The 4 runs on the variants of the 4 October 2026 release |
| `configs/v38_*.toml` | The `--acmg-config` of each run that needs one |
| `provenance/v38_*.txt` | What produced each run's annotation |

The 38 runs on the April set are the commit before the change (`v38_base`), the shipped change (`v38`), the control of the measurement build with nothing set (`v38_alt_control`) and four alternatives (`v38_alt_unguarded`, `v38_alt_label_guard`, `v38_alt_bs2_as_v37`, `v38_alt_triple_cap`), the shipped change with gnomAD's region flags withholding nothing (`v38_noflags`), the 26 sweeps of the four thresholds that read inheritance (`v38_sw_*`) the three runs of the pilot on chromosomes 21 and 22 (`v38_chr2122_*`) and the joint release on every chromosome (`v38_joint`, made on 5 November 2026, once the databases of all 24 chromosomes existed).
The 4 on the October release are `v38_oct_base`, `v38_oct_shipped`, `v38_oct_alt_unguarded` and `v38_oct_alt_triple_cap`.
`v38_joint` is defined in `11_manuscript_runs.py` and is made only when the database of all 24 chromosomes exists.

Things to know when reading these provenance files.

- Three binaries made the runs, all on 4 November 2026 but `v38_joint` (5 November).
  `v38_base` and `v38_oct_base` used the binary of the commit before the change, `9d84be0` (SHA-256 `f6e7f4ca...`), whose classification code is `v37`'s (the two differ in comments only) and whose annotation layer is newer: its run differs from `v37`'s in the criteria of 25 indels and in the call of 4 of them.
  33 runs, the shipped change among them, used the binary of the working tree that became `v38` (`f184277d...`); `v38_joint` is the only one of them made on 5 November 2026, and the only one made once.
  7 runs, the alternatives and their control among them, used a build of `manuscript/acmg/replay_validation/inheritance_alternatives_env.patch` on the same tree (`a2bb20ee...`).
  The last two were built by `scripts/build_manuscript_binaries.sh` after the last edit of the tree.
- The same 41 runs (all but `v38_joint`) were first made with binaries built from earlier states of the tree (SHA-256 `336fea23...`, `a83a5143...` and `0982a653...`), before its last edits: two fixes of allocations, comments, tests, and the changes that followed two independent reviews (`RUN_VERSIONS.md`, "v38").
  The repetition agrees with them in every scored file of every run and in every record of the annotated VCFs that both keep (`rerun_comparison.csv`).
- `repo_commit` is `9d84be0`, the commit the change was made on, for every run but `v38_joint`, because the change was not committed when they were made: `repo_dirty_files` counts the modified and untracked files, which include the change and the results of the runs themselves, and so differs from run to run. `binaries.tsv` says what each binary was built from.
- The environment variable that sets a run of the alternatives build (`FASTVEP_ALT`, with the values `g1_off`, `g3`, `bs2_names,g1_off` and `triple_cap`) is not in the provenance file.
  It is in the `environment` column of `run_matrix_v38.tsv`, which is the only record of it.
- The runs of the pilot read the databases of gnomAD's exomes and joint release for chromosomes 21 and 22, and `v38_joint` those of the joint release for all 24, that `build_gnomad_release_per_chrom.sh` builds; they are not in the repository.
- `v38_joint`'s provenance says `repo_commit` 4078116 with 6 uncommitted files, which are the documentation and the scripts of this change; the binary is the one of the other `v38` runs.
