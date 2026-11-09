# ACMG concordance benchmark

How well does `fastvep annotate --acmg` agree with expert-curated ClinVar calls?

This directory is the single home for that benchmark's **code, documentation and
results**.
Everything it reads and writes at scale - input VCFs, the ~30 GB of supplementary
annotation databases, and the full per-run output directories - lives under
`data/benchmark/`, which is entirely gitignored.

That split is the point of the layout, and it is worth stating why, because the
repository spent several months with the opposite arrangement.
Git cannot re-include a file whose parent directory is excluded.
While `.gitignore` excluded `/data/` and tracked files lived beneath it, every
one of them needed `git add -f`, and three were silently missed - a benchmark
runner and two analysis write-ups that simply never entered the repository.
Nothing under `data/` is tracked now, and nothing needs to be.

**If it is code, documentation or a small result, it belongs here. If it is an
input, a database or a bulk output, it belongs in `data/benchmark/`.**

## Layout

| Path | What |
|---|---|
| [`METHODS.md`](METHODS.md) | What the benchmark measures, how the ClinVar 2-star+ set is built, what the numbers mean (ClinVar's classifications are a reference, not ground truth) |
| [`RUN_VERSIONS.md`](RUN_VERSIONS.md) | Per-run history: what changed between v1 and v38 (the latest version, release 0.5.0), with measured effect, and the record of the full rerun of 4 November 2026 and of the family v38 |
| [`STATUS.md`](STATUS.md) | On-disk inventory - what is downloaded and built locally |
| [`scripts/`](scripts/) | Everything executable, plus [its own README](scripts/README.md) |
| [`scripts/sa_sources/`](scripts/sa_sources/) | Building the supplementary annotation databases from public sources |
| [`results/`](results/) | Small tracked artifacts per run: concordance matrices, summaries, figures |
| [`data/`](data/) | Small tracked *derived* data: curated tables generated from public sources, not bulk inputs |
| [`notes/`](notes/) | Analysis write-ups that are not tied to a single run |

For the *speed* and VEP-equivalence benchmarks, see [`benchmarks/`](../../benchmarks/)
instead. This directory is only about clinical concordance.

## Running it

```bash
# One-time: fetch sources (~30 GB, into data/benchmark/sa_sources/) and build databases
bash analysis/acmg_benchmark/scripts/sa_sources/download_sa_sources.sh

# One-time: extract ClinVar's 2-star+ classifications
python3 analysis/acmg_benchmark/scripts/01_extract_clinvar_2star.py

# Each run - the version string names the output directory
bash analysis/acmg_benchmark/scripts/run_benchmark.sh v15

# Or with a config, to measure one change against the same ClinVar 2-star+ set.
# A variant of a run is a run of its own, so it gets its own version name.
ACMG_CONFIG=analysis/acmg_benchmark/data/vcep_thresholds.toml \
  bash analysis/acmg_benchmark/scripts/run_benchmark.sh v15_vcep

# Any other `fastvep annotate` flag, the same way
EXTRA_ARGS="--pick-order mane_plus_clinical,mane_select,canonical,appris,tsl,biotype,ccds,rank" \
  bash analysis/acmg_benchmark/scripts/run_benchmark.sh v15_mpc

# Any number of runs against a baseline, call by call
python3 analysis/acmg_benchmark/scripts/09_compare_runs.py \
  --baseline v15 --runs v15_vcep v15_mpc --out /tmp/v15_comparison
```

The run script reads the databases of `data/benchmark/sa_db_complete/` (`SA_DIR` names another directory): links to every file of `sa_db/` but its SpliceAI, PhyloP and gnomAD files, which are replaced by 72 files built from the whole of gnomAD's v4.1 exomes file (`results/complete_databases/` holds their checksums and says how they were built and checked, and `scripts/sa_sources/make_complete_stack.sh` makes the directory and stops if a file is not the one that the checksums name).
The SpliceAI, PhyloP and gnomAD files of `sa_db/` itself are not read by the benchmark (`results/complete_databases/README.md`); the runs made with them are kept under the names `manuscript_runs_april_databases/` and `output_v37_april_databases/`, and `results/complete_vs_april_databases.md` sets them against the runs of the paper.
The run script takes about 35 seconds on the ten-core machine of the runs (measured alone on 7 November 2026, with a stack of databases that is smaller than the final one): 25 seconds to annotate 673,660 variants (74 seconds on one thread; `scripts/19_time_annotation.py`, 7 November 2026) and 11 to score, and writes to `data/benchmark/output_v15/`.
It records what produced the annotation in `run_provenance.txt` beside it: the
flags, the hash of the binary and of every `--acmg-config` and
`--functional-evidence` file, which inputs were read (by name and size, not
content) and the commit.
It reuses an existing annotated VCF only when the same binary, flags,
configuration and inputs made it, and stops with exit code 3 otherwise, so a
second run under a name already used cannot be scored against the first run's
annotation.
Record what changed, and its measured effect, in [`RUN_VERSIONS.md`](RUN_VERSIONS.md).

To vary one threshold instead of one version, use the sweep harness rather than a
hand-edited config:

```bash
python3 analysis/acmg_benchmark/scripts/sweep_acmg_thresholds.py \
  --sa-dir data/benchmark/sa_db --out-dir data/benchmark/sweeps \
  --sweep "bp7_max_intron_offset=50,100,300,1000"
```

## The ACMG manuscript

[`manuscript/acmg/`](../../manuscript/acmg/) holds a manuscript on how the classifier was built, reviewed and measured, and on its discordance with ClinVar 2-star+ classifications in its latest version, v38 (release 0.5.0), set against v37, the version in which the review rounds ended.
v37 and v38 are versions of the classifier, and every run of either is on GRCh38.
Every number of the manuscript is v38's unless the text names v37, and the numbers, tables, workbook and figures of v37 are kept as the twins of v38's, so that what v38 changed can be regenerated and checked (`--base v37`).
Its tables are produced from the scored runs by the scripts below, in this order, and its text is checked against them.
Every run is made by a binary that an environment variable names, and not by the checkout's own `target/release/fastvep`.
The checkout's classifier is v38, and a run called `v37_...` made by it would be a run of v38 under v37's name, which `11_manuscript_runs.py` and `14_temporal_holdout.py` refuse to do.
The provenance file of a run records the commit and the number of uncommitted files of the checkout that ran the harness, which is not the commit of the binary; `results/v37_manuscript/binaries_v38.tsv` and `binaries.tsv` say what each binary was built from.

```bash
# 0a. The databases of the runs: the 72 complete files of SpliceAI, PhyloP and gnomAD are copied into data/benchmark/sa_sources/complete/
#     (they were built on another machine, from the whole of the gnomAD file: results/complete_databases/README.md), and the stack that
#     the harness reads is made from them and from the files of sa_db/ that did not change
bash analysis/acmg_benchmark/scripts/sa_sources/make_complete_stack.sh

# 0. The eight binaries, each built from a clean export of its commit in a target directory of its own (the script stops if a
#    patched build does not contain its patch): for the family v38, v38 (the tag v0.5.0-rc.1, release 0.5.0) and cut38 and
#    hook38 (the tag with replay_validation/combiner_cutoff_env.patch and bp_hooks_env.patch); for Section 3.8, base (the commit
#    before the change) and alt (the tag with replay_validation/inheritance_alternatives_env.patch); for the family v37, v37
#    (commit 5d8004a) and cut and hook (the same two patches on it). NAMES="v38 cut38 hook38" builds some of them;
#    V38_REF=<commit> builds v38 at another commit, and V38_TREE=. builds the checkout as it is, uncommitted files included.
#    It prints the eight `export` lines that the steps below read their binaries from.
analysis/acmg_benchmark/scripts/build_manuscript_binaries.sh /path/to/binaries
export FASTVEP_V38_BIN=/path/to/binaries/fastvep_v38 FASTVEP_CUTOFF38_BIN=/path/to/binaries/fastvep_cut38 \
       FASTVEP_HOOKS38_BIN=/path/to/binaries/fastvep_hook38 FASTVEP_V38_BASE_BIN=/path/to/binaries/fastvep_base \
       FASTVEP_V38_ALT_BIN=/path/to/binaries/fastvep_alt FASTVEP_V37_BIN=/path/to/binaries/fastvep_v37 \
       FASTVEP_CUTOFF_BIN=/path/to/binaries/fastvep_cut FASTVEP_HOOKS_BIN=/path/to/binaries/fastvep_hook

# 0b. v37 itself, which every run of the family v37 differs from by one setting, and the earlier runs of v37 that Section 3.8 pairs
#     with its own sweeps (the configurations are in results/v37/configs)
FASTVEP="$FASTVEP_V37_BIN" bash analysis/acmg_benchmark/scripts/run_benchmark.sh v37
for n in bs2off bs2prev_3e-5 bs2prev_1e-4 bs2prev_1e-3 bs2prev_3e-3 bs2ac_20 bs2ac_100 bs2ac_1000; do
  ACMG_CONFIG=analysis/acmg_benchmark/results/v37/configs/$n.toml FASTVEP="$FASTVEP_V37_BIN" \
    bash analysis/acmg_benchmark/scripts/run_benchmark.sh v37_$n
done
FASTVEP="$FASTVEP_V37_BIN" EXTRA_ARGS="--pick-order mane_plus_clinical,mane_select,canonical,appris,tsl,biotype,ccds,rank" \
  bash analysis/acmg_benchmark/scripts/run_benchmark.sh v37_mpc
FASTVEP="$FASTVEP_V37_BIN" EXTRA_ARGS="--pick-order rank,mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds" \
  bash analysis/acmg_benchmark/scripts/run_benchmark.sh v37_rank

# 1. The three runs that change an input need the ClinGen file completed first
#    (the commands are in the docstring of 15_clingen_complete_export.py)

# 2. The runs of Section 3.8 (what v38 changed: ClinGen's curated mode of inheritance is read), which need the binaries v38, base
#    and alt. The databases of gnomAD's joint exomes and genomes release come first, for all 24 chromosomes: about 6 hours from
#    gnomAD's public files, limited by the network and not by the script (JOBS=8 runs eight chromosomes at a time, the default is
#    four), and resumable, so that a run that was stopped continues after its last finished batch. 11_manuscript_runs.py makes
#    `v38_joint` once all 24 exist; the pilot on chr21 and chr22 that was made first reads the same two databases, and its input,
#    the benchmark's variants on those two chromosomes, is written by 11_manuscript_runs.py from the benchmark's own files
FASTVEP="$FASTVEP_V38_BIN" JOBS=8 analysis/acmg_benchmark/scripts/sa_sources/build_gnomad_release_per_chrom.sh joint
python3 analysis/acmg_benchmark/scripts/11_manuscript_runs.py --v38 --workers 3
python3 analysis/acmg_benchmark/scripts/11_manuscript_runs.py --v38 --holdout --workers 3

# 3. The family v38, which is the default: the sweeps, ablations, interaction and gap runs and the runs of the two measurement
#    builds (about 3 minutes each; --workers 3 is comfortable), and the seven runs of the changes that Section 4.5 of the
#    manuscript proposes, made on the variants of the later ClinVar release. 18 of the 96 are runs of step 2 that have the same
#    names and are not made again, and the script stops if step 2 defines one of them differently
python3 analysis/acmg_benchmark/scripts/11_manuscript_runs.py --workers 3
python3 analysis/acmg_benchmark/scripts/11_manuscript_runs.py --holdout --workers 3

# 4. The same 96 runs by the classifier of v37, which are what Section 3.8 sets v38 against, and the later ClinVar release by v37
#    (prepare writes the strata of the release, which the two versions share, and the 2-star+ variants of the whole release are scored)
python3 analysis/acmg_benchmark/scripts/11_manuscript_runs.py --family v37 --workers 3
python3 analysis/acmg_benchmark/scripts/11_manuscript_runs.py --family v37 --holdout --workers 3
python3 analysis/acmg_benchmark/scripts/14_temporal_holdout.py --base v37 all

# 5. The per-variant table of v37 and of v38 and of each run whose variants the statistics read one by one
python3 analysis/acmg_benchmark/scripts/12_per_variant_table.py \
  --run data/benchmark/output_v37 --out data/benchmark/manuscript_runs/v37_per_variant.tsv.gz
for fam in v37 v38; do
  for run in ${fam}_gap_clingen_complete ${fam}_gap_orphanet_lof ${fam}_gap_both ${fam}_hk_bp1_loo ${fam}_sw_phylop_1000 \
             ${fam}_hk_bp7_walker ${fam}_sw_bp1veto_1 ${fam}_sw_bp1veto_2 ${fam}_hk_bp1_loo_veto1 ${fam}_hk_bp1_loo_veto2; do
    python3 analysis/acmg_benchmark/scripts/12_per_variant_table.py \
      --run data/benchmark/manuscript_runs/$run --out data/benchmark/manuscript_runs/${run}_per_variant.tsv.gz
  done
  for run in ${fam}_oct_clingen_complete ${fam}_oct_orphanet_lof ${fam}_oct_both ${fam}_oct_bp7_walker ${fam}_oct_bp1_loo \
             ${fam}_oct_bp1_veto1 ${fam}_oct_no_an_floor; do
    python3 analysis/acmg_benchmark/scripts/12_per_variant_table.py \
      --run data/benchmark/manuscript_runs/$run --labels data/benchmark/clinvar_20261004/clinvar_2star_labels.tsv \
      --out data/benchmark/manuscript_runs/${run}_per_variant.tsv.gz
  done
done
# and of the runs of Section 3.8 (v38 itself and its October run are in the loops of the family)
for run in v38 v38_base v38_alt_unguarded v38_alt_label_guard v38_alt_bs2_as_v37 v38_alt_triple_cap \
           v38_noflags v38_joint v38_chr2122_exomes v38_chr2122_exomes_noflags v38_chr2122_joint; do
  python3 analysis/acmg_benchmark/scripts/12_per_variant_table.py \
    --run data/benchmark/manuscript_runs/$run --out data/benchmark/manuscript_runs/${run}_per_variant.tsv.gz
done
for run in v38_oct_base v38_oct_shipped v38_oct_alt_unguarded v38_oct_alt_triple_cap; do
  python3 analysis/acmg_benchmark/scripts/12_per_variant_table.py \
    --run data/benchmark/manuscript_runs/$run --labels data/benchmark/clinvar_20261004/clinvar_2star_labels.tsv \
    --out data/benchmark/manuscript_runs/${run}_per_variant.tsv.gz
done
python3 analysis/acmg_benchmark/scripts/14_temporal_holdout.py table

# 6. Why BA1 and BS1 did not fire on the Benign variants of 2% or more that are not called benign-side (needs the tables above),
#    for each version
FASTVEP_V38_BIN="$FASTVEP_V38_BIN" python3 analysis/acmg_benchmark/scripts/14_temporal_holdout.py blockers
FASTVEP_V37_BIN="$FASTVEP_V37_BIN" python3 analysis/acmg_benchmark/scripts/14_temporal_holdout.py --base v37 blockers

# 6b. The variants that the joint release turns into errors, or out of a pathogenic-side call, annotated one at a time
#     through both databases (the filtering allele frequency and the flags are in the json and not in the VCF)
FASTVEP_V38_BIN="$FASTVEP_V38_BIN" python3 analysis/acmg_benchmark/scripts/17_joint_boundary_cases.py

# 6c. The status of the 28 criteria against the engine's json output, on every 40th record of the benchmark (the VCF lists
#     only the criteria that were met; the json says whether each could be evaluated). It stops, writing nothing, if the engine
#     disagrees with the groups of criteria that the script states
FASTVEP_V38_BIN="$FASTVEP_V38_BIN" python3 analysis/acmg_benchmark/scripts/18_criteria_status.py

# 7. The runs of the paper against the same runs made with the databases of April 2026 (results/complete_vs_april_databases.md and the
#    databases_*.csv files of results/complete_vs_april_databases/: an internal record, which no table of the paper holds; the trajectory
#    reads databases_old_vs_new.csv, so this comes first): it reads the per-variant tables of the runs kept under the names
#    manuscript_runs_april_databases/ and output_v37_april_databases/
python3 analysis/acmg_benchmark/scripts/21_compare_databases.py
#    Then every table (the files of v38, which the manuscript reads, and then their twins for v37), then the figures, then the checks
python3 analysis/acmg_benchmark/scripts/13_manuscript_stats.py
python3 analysis/acmg_benchmark/scripts/13_manuscript_stats.py --base v37
python3 manuscript/acmg/generate_figures.py
python3 manuscript/acmg/build_tables.py --write
python3 manuscript/acmg/build_tables.py --check
python3 manuscript/acmg/build_supplement.py --write
python3 manuscript/acmg/build_supplement.py --check
python3 manuscript/acmg/verify_numbers.py

# 8. The comparison of the two versions, run for run (results/v38_family_vs_v37.md and its five CSV files; its Section 9 says
#    which statements and numbers of the manuscript would read differently under v37)
python3 analysis/acmg_benchmark/scripts/20_compare_families.py

# 9. What v38 changed, from the twins: the sentences whose numbers differ under v37, and the tables, workbook and figures of v37,
#    which are written outside the repository because the ones in it are the manuscript's
python3 manuscript/acmg/verify_numbers.py --base v37 --dump /tmp/checks_v37.json
python3 manuscript/acmg/build_tables.py --base v37 --check
python3 manuscript/acmg/build_supplement.py --base v37 --write --out /tmp/supplementary_tables_v37.xlsx
python3 manuscript/acmg/generate_figures.py --base v37 --out /tmp/figures_v37
```

Run without `--base`, `13_manuscript_stats.py`, `14_temporal_holdout.py`, `build_tables.py`, `build_supplement.py`, `generate_figures.py` and `verify_numbers.py` are about v38, the files of which the manuscript reads, and the three checks pass on the manuscript as it stands.
With `--base v37` they are about v37: `13_manuscript_stats.py` writes the twins of the tables (`_v37` where the file of v38 has `_v38`), and writes them as they were when v37 was the base but for the rows that the text quotes (`facts_v37.csv` and `clingen_missing_summary_v37.csv` have more rows), the other scripts write the tables, workbook and figures of v37 outside the repository, and `verify_numbers.py` fails on the sentences whose numbers v38 changed, which is its list of them and not a fault of the data.
Seven files have no twin: `trajectory.csv`, `change_ledger.csv`, `review_rounds.csv`, `round6_rulings.tsv`, `clingen_missing_curations.csv`, `criteria_status.csv` and `tavtigian_posterior.csv`, which are the record of the numbered versions (the trajectory ends with v38), the reviewer's rulings or independent of the classifier.

`run_benchmark.sh` and the scripts that start it read the environment, so a run that is made by hand needs `FASTVEP` set the same way.
The tables and the checks are deterministic: after the statistics script had run under each of two hash seeds (`PYTHONHASHSEED=0` and `1`) the files of `manuscript/acmg/data/` that it writes were the same, and the figures are byte for byte the same on a rerun.
Some files are not written by the statistics script: `change_ledger.csv` and `round6_rulings.tsv` are typed from the project's record, `clingen_missing_curations.csv` comes from `15_clingen_complete_export.py`, `temporal_high_af_blockers_v38.csv` from `14_temporal_holdout.py blockers` (and `_v37` from `--base v37 blockers`), `v38_joint_errors.csv` and `v38_joint_pm2_lost.csv` from `17_joint_boundary_cases.py`, and `criteria_status.csv` from `18_criteria_status.py`.

The inputs are dated and none of them is in the repository, and the download scripts fetch the current files and not the dated ones: the ClinVar weekly releases of 4 April 2026 (the benchmark and its labels) and of 4 October 2026, the ClinGen Gene-Disease Validity exports of 30 April 2026 (the file the benchmark reads) and of 11 October 2026, and the ClinVar `variant_summary` and gnomAD files of Section 2.2 of the manuscript.
A rerun on other releases of them is a rerun of the method and not of the paper's numbers.
The run records are in [`results/v37_manuscript/`](results/v37_manuscript/) (the directory keeps the name it had when it held the runs of v37 alone; the runs of v38 are in the files whose names end in `_v38`).
`scripts/16_compare_runs.py` compares a rerun with the runs it repeats, file by file and record by record, and writes `results/v37_manuscript/rerun_comparison.csv`; use it when any step above is made again.

## The classifier itself

Thresholds, criteria and the reasoning behind every value set by measurement are
documented in [`docs/ACMG.md`](../../docs/ACMG.md), not here.
This directory measures the classifier; it does not define it.

The medical-genetics review that drove most of the criteria changes is in
[`docs/ACMG_EXPERT_REVIEW_ROUND2.md`](../../docs/ACMG_EXPERT_REVIEW_ROUND2.md).
What each subsequent run did about it is in
[`RUN_VERSIONS.md`](RUN_VERSIONS.md); the runs between reviews are ours, not hers.

## Curated tables in `data/`

Two tables are generated from public sources and tracked, because they are small,
because they are inputs to classification rather than outputs of a run, and
because their provenance needs to be reviewable line by line.

| File | What | Regenerate with |
|---|---|---|
| [`data/vcep_thresholds.toml`](data/vcep_thresholds.toml) | Published ClinGen VCEP BA1/BS1/PM2 thresholds for 117 genes, as `--acmg-config` input | [`scripts/build_vcep_thresholds.py`](scripts/build_vcep_thresholds.py) |
| [`data/vcep_thresholds_audit.tsv`](data/vcep_thresholds_audit.tsv) | Every threshold and every rejection, with the verbatim sentence it came from | same script |
| [`data/gene_disease_attributes.tsv`](data/gene_disease_attributes.tsv) | Prevalence, onset and inheritance per gene-disorder pair, joined to the VCEP thresholds | [`scripts/build_gene_disease_attributes.py`](scripts/build_gene_disease_attributes.py) |

The audit table is the point. A threshold read out of prose by a regex is a
claim about a clinical guideline, and the only way to make that claim checkable
is to keep the sentence next to the number. Every rejection is in there too,
with its reason, so "this gene has no threshold" and "this generator would not guess"
stay distinguishable.

`vcep_thresholds.toml` is **not** on by default. It is a measured improvement
that also introduces four new reference-discordant benign calls, two of which are founder
variants the panels excluded in prose the generator cannot read. See the v17
section of [`RUN_VERSIONS.md`](RUN_VERSIONS.md) for the full trade.
