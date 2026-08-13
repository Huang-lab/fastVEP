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
| [`METHODS.md`](METHODS.md) | What the benchmark measures, how the truth set is built, what the numbers mean |
| [`RUN_VERSIONS.md`](RUN_VERSIONS.md) | Per-run history: what changed between v1 and v14, with measured effect |
| [`STATUS.md`](STATUS.md) | On-disk inventory - what is downloaded and built locally |
| [`scripts/`](scripts/) | Everything executable, plus [its own README](scripts/README.md) |
| [`scripts/sa_sources/`](scripts/sa_sources/) | Building the supplementary annotation databases from public sources |
| [`results/`](results/) | Small tracked artifacts per run: concordance matrices, summaries, figures |
| [`notes/`](notes/) | Analysis write-ups that are not tied to a single run |

For the *speed* and VEP-equivalence benchmarks, see [`benchmarks/`](../../benchmarks/)
instead. This directory is only about clinical concordance.

## Running it

```bash
# One-time: fetch sources (~30 GB, into data/benchmark/sa_sources/) and build databases
bash analysis/acmg_benchmark/scripts/sa_sources/download_sa_sources.sh

# One-time: build the truth set
python3 analysis/acmg_benchmark/scripts/01_extract_clinvar_2star.py

# Each run - the version string names the output directory
bash analysis/acmg_benchmark/scripts/run_benchmark.sh v15
```

The run script takes about 90 seconds to annotate 673,660 variants plus a couple
of minutes to score, reuses an existing annotated VCF if one is present, and
writes to `data/benchmark/output_v15/`.
Record what changed, and its measured effect, in [`RUN_VERSIONS.md`](RUN_VERSIONS.md).

To vary one threshold instead of one version, use the sweep harness rather than a
hand-edited config:

```bash
python3 analysis/acmg_benchmark/scripts/sweep_acmg_thresholds.py \
  --sa-dir data/benchmark/sa_db --out-dir data/benchmark/sweeps \
  --sweep "bp7_max_intron_offset=50,100,300,1000"
```

## The classifier itself

Thresholds, criteria and the reasoning behind every value set by measurement are
documented in [`docs/ACMG.md`](../../docs/ACMG.md), not here.
This directory measures the classifier; it does not define it.

The expert-review rounds that drove most of the criteria changes are in
[`docs/ACMG_EXPERT_REVIEW_ROUND2.md`](../../docs/ACMG_EXPERT_REVIEW_ROUND2.md) and
[`docs/REVIEWER_REPLY_ROUND3.md`](../../docs/REVIEWER_REPLY_ROUND3.md).
