# Results: v36 to v37, the round-6 review

This is the benchmark record for the change that answers `discordant_64_round6.xlsx`.
What each change did, and what it cost, is in [`RUN_VERSIONS.md`](../../RUN_VERSIONS.md), section "v36 to v37".
The note sent back to the medical geneticist is [`MD_REVIEW_ROUND7.md`](MD_REVIEW_ROUND7.md).

All comparisons are with ClinVar's 2-star+ classifications, which are a reference and not ground truth.

| File | Content |
|---|---|
| `concordance_summary.txt`, `concordance_matrix.csv`, `criterion_firing_rates.csv`, `rule_distribution.csv` | The headline agreement and the criterion statistics of the final run, `v37`. |
| `provenance/<run>.txt` | What produced each run's annotation: the flags, the hash of the binary and of every `--acmg-config`, which inputs were read, and the commit. |
| `call_changes_v36_to_v37.tsv` | Every one of the 2,052 calls that changed between master (`v36`) and `v37`, with the criteria before and after. |
| `run_summary.tsv`, `summary.md` | `09_compare_runs.py` over `v36`, `v37` and its variants. The headline counts per run, the two error types, how many of the reviewed rows are still opposite-direction, and every changed call grouped by the criteria it lost and gained. |
| `attribution_v36_to_*.tsv` | The same grouping, one file per run. |
| `configs/` | The `--acmg-config` of each BS2 sweep run, one TOML named for the run. |

## The runs

`v36` is master at the start of the round, re-run.
It predates `run_provenance.txt`, so it has none.
Every run below was made by `run_benchmark.sh` at commit `5d2521f` of this branch, with the binary built from it (hash in each provenance file).
The two commits above it change documentation and results, and the line endings of the tables `03` and `08` write.
The runs that need a configuration or a flag are made with the `ACMG_CONFIG` or `EXTRA_ARGS` in the last column, from the repository root.

| Run | What differs from `v37` | How to make it |
|---|---|---|
| `v37_bs2off` | `bs2_ad_prevalence_threshold = 0`: the old dominant-gene BS2 floor alone. | `ACMG_CONFIG=analysis/acmg_benchmark/results/v37/configs/bs2off.toml` |
| `v37_bs2prev_3e-5`, `_1e-4`, `_1e-3`, `_3e-3` | The dominant carrier bar at that value. The shipped bar, 3e-4, is `v37` itself. | `ACMG_CONFIG=.../configs/bs2prev_<value>.toml` |
| `v37_bs2ac_20`, `_100`, `_1000` | The bar off and the observed-carrier floor raised from 5. | `ACMG_CONFIG=.../configs/bs2ac_<value>.toml` |
| `v37_mpc` | `--pick-order mane_plus_clinical,mane_select,canonical,appris,tsl,biotype,ccds,rank`. | `EXTRA_ARGS="--pick-order mane_plus_clinical,mane_select,canonical,appris,tsl,biotype,ccds,rank"` |
| `v37_rank` | `--pick-order rank,mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds`. | `EXTRA_ARGS="--pick-order rank,mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds"` |
| `v37_adtag` | BS2 reads ClinGen's `/AD` and `/AR` tags. An experiment, not shipped; see limitation 20 in `docs/ACMG.md`. | A locally patched build of `fastvep` that is not in this repository: no configuration or flag does it. |
| `v37_moi_global` | Every criterion reads the tags. An experiment, not shipped. | The same. |

`v37_adtag` and `v37_moi_global` were run on the stack as it stood before PS1's splice path, PM1 and BP3 learned to follow the HGVS placement, so they have no provenance file here.
That change altered the criteria of 43 variants and the call of none, so the two runs compare with `v37` call for call.

To make a run, give it a version name and name its configuration or flags:

```bash
ACMG_CONFIG=analysis/acmg_benchmark/results/v37/configs/bs2prev_1e-4.toml \
  bash analysis/acmg_benchmark/scripts/run_benchmark.sh v37_bs2prev_1e-4
```

`run_benchmark.sh` refuses to reuse an annotation that something other than the current binary, inputs, flags and configuration produced, and it records what it ran in the run's `run_provenance.txt`.
The BS2 bar and floor sweeps are tabulated in [`docs/ACMG.md`](../../../../docs/ACMG.md#the-dominant-bar) and in `RUN_VERSIONS.md`.
