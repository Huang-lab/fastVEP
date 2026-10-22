# Runs made for the ACMG manuscript

Every run here is a variant of `v37`, the final classifier, made by `run_benchmark.sh` and scored by `03_evaluate_concordance.py`.
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

89 runs on the April set: 41 sweeps, 16 single-switch ablations, 9 leave-one-source-out ablations, 4 PM2 sweeps under the 2015 table, 9 runs of the point cutoffs, 7 runs of the second measurement build (`hook`) and 3 runs that change an input (`gap`).
Seven more runs (`holdout`) apply one change each to the variants of the October 2026 release, and their matrices cover all 772,585 variants of that release's covered subset; `13_manuscript_stats.py` scores them on the 107,528 new ones from their per-variant tables.

Things to know when reading the provenance files.

- 73 runs used the binary of `v37` (SHA-256 `d12c0d8e...`), built at commit `5d2521f`; the commits above it up to the head of that branch add tests, comments and two guards that no benchmark input reaches, and a binary built from the head produced an annotation byte-identical to `v37`'s (`RUN_VERSIONS.md`).
- The 9 cutoff runs (`v37_cut_*`) used a binary built from `manuscript/acmg/replay_validation/combiner_cutoff_env.patch`, which differs from the shipped one in eight lines (six added, two removed), so its SHA-256 (`dc275d1c...`) is not that of the other runs.
  The run at the shipped cutoff (`v37_cut_lp6`) is a control and reproduces `v37` in every cell.
- The 7 `hook` runs used a binary built from `manuscript/acmg/replay_validation/bp_hooks_env.patch` (SHA-256 `2bccb937...`), which reads BP1's own-record count and BP7's SpliceAI bar from the environment.
  `v37_hk_control`, with nothing set, is its control and reproduces `v37` in every cell and in the firing count of every criterion.
- The 7 `holdout` runs used the binary of `v37`, except `v37_oct_bp7_walker` and `v37_oct_bp1_loo`, which used the second measurement build; the three that change the ClinGen file or Orphanet's label need the files of the `gap` runs.
- The 3 `gap` runs need files that are not in the repository: the ClinGen export with its omitted curations put back (`scripts/15_clingen_complete_export.py` builds it) and, for two of them, a configuration of 1,042 gene overrides that `11_manuscript_runs.py` writes from `analysis/acmg_benchmark/data/gene_disease_attributes.tsv`.
- `repo_commit` is the commit of the checkout that ran the harness, which is not always the commit the binary was built from; the binary is identified by its hash.
  63 of the 89 files record `5d8004a` and 26 record `1dfa488`, the head of the branch that was squash-merged as `5d8004a`, whose tree is the same.
- The environment variable that sets a run of a patched build (`FASTVEP_LP_CUT`, `FASTVEP_P_CUT`, `FASTVEP_BP1_LOO` or `FASTVEP_BP7_SPLICEAI_BENIGN`) is not part of what the harness compares and is not in the provenance file.
  It is in the `environment` column of `run_matrix.tsv`, which is the only record of it.
- `repo_dirty_files` counts the files that were modified or untracked when a run was made.
  These runs were made while the manuscript, its scripts and these results were untracked, so the count is not zero and says nothing about the classifier, whose crate was at the recorded commit.
