#!/usr/bin/env bash
# Two builds of fastVEP on small variants: do they write the same bytes, and do they take the same time?
#
#   benchmarks/run_small_variant_regression.sh RESULTS_DIR BASE_BINARY NEW_BINARY
#
# Read from the environment: GFF3 and FASTA (as for run_head_to_head.py), HG002_VCF (the whole call set
# of GIAB HG002, with its genotype columns).
#
# Six workloads, each run once by each binary per round, the two in turn, the order alternating from one
# round to the next, so that a slow minute of the machine, and the heat of the run before, fall on both
# (with the same order every round, the second run of a pair of 30-second runs was 1 to 2% slower for it),
# after one round that is not timed and where the outputs are compared byte for byte (`cmp`):
#   sample_1_in_20    every 20th record of the call set (202,417), the flags of the head-to-head with Ensembl VEP
#                     (`--hgvs --hgvsg --symbol --canonical --numbers --variant-class --show-ref-allele --total-length`),
#                     tab output, every core, 12 rounds
#   chr22             chromosome 22 (50,284 variants), `--hgvs --symbol --canonical`, VCF output, one thread, 5 rounds
#   genome            the whole call set (4,048,342), consequences only (no FASTA), VCF output, one thread (5 rounds)
#                     and every core (8 rounds: the whole genome on every core is the cell that varies most from round to round)
#   genome_hgvs       the same with the FASTA and `--hgvs --symbol --canonical`, which is what Table 3 of the manuscript measures
#                     as "consequence only" (it has no database): one thread, 3 rounds, and every core, 5 rounds
#   genome_1_in_5     every 5th record of the call set (809,668), the flags of the whole genome, VCF output, one
#                     thread, 10 rounds: a difference per row shows in the whole-genome run, which is too long
#                     to take enough rounds of
#   genome_1_in_5_tab the same, tab output
#
# More, each run only when the directory or the file it reads is named in the environment, for the conditions of the
# head-to-head and the organisms that the six above do not reach (SA_CLINVAR and SA_JOINT are the annotation directories that
# run_head_to_head.py makes: ClinVar's `clinvar.osa`, with the joint release of gnomAD's `gnomad_chr*.osa2` for the second;
# SA_COMPLETE is the directory of the ACMG benchmark, CLINVAR_2STAR its input):
#   genome_1_in_5_clinvar, genome_1_in_5_joint   the flags of the whole genome and the databases, VCF output, one thread, 10 rounds
#   yeast_full, drosophila_full, arabidopsis_full, mouse_full
#                                                 the whole variant file of each organism of Table 7 with `--hgvs`, VCF output, every
#                                                 core (ORGANISMS_DIR is test_data/organisms: its `<organism>/<organism>.gff3` and `.fa`)
#   genome_1_in_5_json                            the same without databases, JSON output, one thread, 5 rounds
#   clinvar_1_in_10_acmg                          every 10th record of the ClinVar 2-star+ set, `--acmg --pick` over SA_COMPLETE,
#                                                 every core, 5 rounds
#
# BASE_LABEL and NEW_LABEL (default `base` and `new`) say in environment.tsv what the two binaries are built from, for instance
# `master ed17877` and `structural variants, uncommitted tree on ed17877`.
#
# The output of a run is removed before the next one starts. The outputs of the whole-genome and the larger organism runs are 10 to
# 22 GB; written over the same two files round after round they left the disk in a state that moved the time of a run on every core
# by up to a factor of two (9.9 s to 20 s for one command), and in one series put one build 2% to 5% behind the other in eight pairs of
# eight. Four builds in turn, each writing a file that did not exist and removing it, took 9.80, 9.79, 9.88 and 9.81 s (medians of six).
#
# Writes RESULTS_DIR/timings.tsv (every run: its wall time, the instructions it retired and its peak resident memory, which
# `/usr/bin/time -l` reports: the instructions of a run on one thread are the same to 0.1% from one run to the next, which no wall time
# on a laptop is; on every core the threads that wait move them by 0.5 to 0.8%),
# summary.tsv (the medians), identity.tsv (the comparisons) and environment.tsv (the machine, the binaries with their SHA-256 and the
# commit of the checkout this script is in).
set -euo pipefail
export LC_ALL=C

# One `caffeinate -i` around the whole script, which keeps a laptop from idle-sleeping between the runs. It is not put in front of
# each run: as the parent of `/usr/bin/time -l` it makes that report its own counters (20 million instructions for a run of 24 billion).
if [ -z "${SVR_UNDER_CAFFEINATE:-}" ] && command -v caffeinate >/dev/null 2>&1; then
  export SVR_UNDER_CAFFEINATE=1
  exec caffeinate -i "$BASH" "$0" "$@"
fi

results="${1:?usage: run_small_variant_regression.sh RESULTS_DIR BASE_BINARY NEW_BINARY}"
base="${2:?}"
new="${3:?}"
: "${GFF3:?set GFF3}" "${FASTA:?set FASTA}" "${HG002_VCF:?set HG002_VCF}"
mkdir -p "$results/work"
results="$(cd "$results" && pwd)"
work="$results/work"

{
  printf 'started_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'system\t%s\n' "$(uname -srm)"
  printf 'cpu\t%s\n' "$(sysctl -n machdep.cpu.brand_string 2>/dev/null || true)"
  printf 'load_at_start\t%s\n' "$(uptime | sed 's/.*load averages*: //')"
  # `caffeinate -i` is what keeps a laptop from idle-sleeping between the runs; a shell that cannot take the assertion (some
  # sandboxes) says so on stderr and the runs go on without it, so it is recorded
  awake="$(caffeinate -i -t 1 2>&1 || true)"
  case "$awake" in *Failed*) printf 'caffeinate_holds_the_machine_awake\tno\n' ;; *) printf 'caffeinate_holds_the_machine_awake\tyes\n' ;; esac
  for which in base new; do
    bin="$base"; label="${BASE_LABEL:-base}"
    [ "$which" = new ] && { bin="$new"; label="${NEW_LABEL:-new}"; }
    printf 'binary_%s\t%s\t%s\tsha256 %s\n' "$which" "$label" "$("$bin" --version 2>&1 | head -1)" "$(shasum -a 256 "$bin" | cut -d' ' -f1)"
  done
  repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
  printf 'repo_head\t%s\n' "$(git -C "$repo" rev-parse HEAD 2>/dev/null || true)"
  printf 'repo_dirty_files\t%s\n' "$(git -C "$repo" status --porcelain 2>/dev/null | wc -l | tr -d ' ')"
  printf 'gff3\t%s\nfasta\t%s\nhg002_vcf\t%s\n' "$GFF3" "$FASTA" "$HG002_VCF"
  # both builds read the sidecar cache beside the GFF3, and a run that needs sequences it lacks rewrites it: it is hashed before and after
  printf 'sidecar_cache_md5_before\t%s\n' "$(md5 -q "$GFF3.fastvep.cache" 2>/dev/null || echo absent)"
} > "$results/environment.tsv"

seconds() {
  local TIMEFORMAT=%R
  { time "$@" >/dev/null 2>&1; } 2>&1
}
timed=(/usr/bin/time -l -o "$work/last.time")
# what `/usr/bin/time -l` said of the last run: its instructions retired and its peak resident memory
counter() { awk -v what="$1" '$0 ~ what { print $1; exit }' "$work/last.time"; }

[ -s "$work/sample_1_in_20.vcf" ] || awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} (++n % 20 == 0)' "$HG002_VCF" > "$work/sample_1_in_20.vcf"
[ -s "$work/genome_1_in_5.vcf" ] || awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} (++n % 5 == 0)' "$HG002_VCF" > "$work/genome_1_in_5.vcf"
[ -z "${CLINVAR_2STAR:-}" ] || [ -s "$work/clinvar_1_in_10.vcf" ] || awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} (++n % 10 == 0)' "$CLINVAR_2STAR" > "$work/clinvar_1_in_10.vcf"
[ -s "$work/chr22.vcf" ] || awk 'BEGIN{FS=OFS="\t"} /^#/ || $1 == "chr22" || $1 == "22"' "$HG002_VCF" > "$work/chr22.vcf"

# WORKLOAD BINARY THREADS OUT
run() {
  local workload="$1" bin="$2" threads="$3" out="$4"
  if [ "$threads" = 1 ]; then export RAYON_NUM_THREADS=1; else unset RAYON_NUM_THREADS; fi
  case "$workload" in
    sample_1_in_20)
      "${timed[@]}" "$bin" annotate -i "$work/sample_1_in_20.vcf" --gff3 "$GFF3" --fasta "$FASTA" \
        --hgvs --hgvsg --symbol --canonical --numbers --variant-class --show-ref-allele --total-length \
        --output-format tab --no-progress -o "$out" ;;
    chr22)
      "${timed[@]}" "$bin" annotate -i "$work/chr22.vcf" --gff3 "$GFF3" --fasta "$FASTA" \
        --hgvs --symbol --canonical --output-format vcf --no-progress -o "$out" ;;
    genome)
      "${timed[@]}" "$bin" annotate -i "$HG002_VCF" --gff3 "$GFF3" \
        --output-format vcf --no-progress -o "$out" ;;
    genome_hgvs)
      "${timed[@]}" "$bin" annotate -i "$HG002_VCF" --gff3 "$GFF3" --fasta "$FASTA" \
        --hgvs --symbol --canonical --output-format vcf --no-progress -o "$out" ;;
    genome_1_in_5)
      "${timed[@]}" "$bin" annotate -i "$work/genome_1_in_5.vcf" --gff3 "$GFF3" \
        --output-format vcf --no-progress -o "$out" ;;
    genome_1_in_5_tab)
      "${timed[@]}" "$bin" annotate -i "$work/genome_1_in_5.vcf" --gff3 "$GFF3" \
        --output-format tab --no-progress -o "$out" ;;
    genome_1_in_5_clinvar)
      "${timed[@]}" "$bin" annotate -i "$work/genome_1_in_5.vcf" --gff3 "$GFF3" --fasta "$FASTA" \
        --hgvs --symbol --canonical --sa-dir "$SA_CLINVAR" --output-format vcf --no-progress -o "$out" ;;
    genome_1_in_5_joint)
      "${timed[@]}" "$bin" annotate -i "$work/genome_1_in_5.vcf" --gff3 "$GFF3" --fasta "$FASTA" \
        --hgvs --symbol --canonical --sa-dir "$SA_JOINT" --output-format vcf --no-progress -o "$out" ;;
    yeast_full|drosophila_full|arabidopsis_full|mouse_full)
      local organism="${workload%_full}" vcf
      case "$organism" in
        yeast) vcf="$ORGANISMS_DIR/yeast/yeast_ensembl_full.vcf" ;;
        drosophila) vcf="$ORGANISMS_DIR/drosophila/drosophila_dgrp2_full.vcf" ;;
        arabidopsis) vcf="$ORGANISMS_DIR/arabidopsis/arabidopsis_1001g_full.vcf" ;;
        mouse) vcf="$ORGANISMS_DIR/mouse/mouse_cast_eij_full.vcf" ;;
      esac
      "${timed[@]}" "$bin" annotate -i "$vcf" --gff3 "$ORGANISMS_DIR/$organism/$organism.gff3" --fasta "$ORGANISMS_DIR/$organism/$organism.fa" \
        --hgvs --output-format vcf --no-progress -o "$out" ;;
    genome_1_in_5_json)
      "${timed[@]}" "$bin" annotate -i "$work/genome_1_in_5.vcf" --gff3 "$GFF3" --fasta "$FASTA" \
        --hgvs --symbol --canonical --output-format json --no-progress -o "$out" ;;
    clinvar_1_in_10_acmg)
      "${timed[@]}" "$bin" annotate -i "$work/clinvar_1_in_10.vcf" --gff3 "$GFF3" --fasta "$FASTA" \
        --sa-dir "$SA_COMPLETE" --acmg --pick --output-format vcf --no-progress -o "$out" ;;
  esac
}

printf 'workload\tbuild\tthreads\tround\twall_s\tinstructions\tmax_rss_bytes\n' > "$results/timings.tsv"
printf 'workload\tthreads\tidentical\tbytes\tbytes_base\n' > "$results/identity.tsv"

# workload threads rounds
cells=("sample_1_in_20 all 12" "chr22 1 5" "genome_1_in_5 1 10" "genome_1_in_5_tab 1 10" "genome 1 5" "genome all 8" "genome_hgvs 1 3" "genome_hgvs all 5")
[ -z "${SA_CLINVAR:-}" ] || cells+=("genome_1_in_5_clinvar 1 10")
[ -z "${SA_JOINT:-}" ] || cells+=("genome_1_in_5_joint 1 10")
cells+=("genome_1_in_5_json 1 5")
[ -z "${ORGANISMS_DIR:-}" ] || cells+=("yeast_full all 5" "drosophila_full all 5" "arabidopsis_full all 3" "mouse_full all 3")
[ -z "${SA_COMPLETE:-}" ] || [ -z "${CLINVAR_2STAR:-}" ] || cells+=("clinvar_1_in_10_acmg all 5")
for cell in "${cells[@]}"; do
  read -r workload threads rounds <<< "$cell"
  echo "$workload $threads: warm-up and comparison" >&2
  run "$workload" "$base" "$threads" "$work/base.out" >/dev/null 2>&1
  run "$workload" "$new" "$threads" "$work/new.out" >/dev/null 2>&1
  # two empty outputs are equal, and say nothing: both must have something in them
  if [ -s "$work/base.out" ] && [ -s "$work/new.out" ] && cmp -s "$work/base.out" "$work/new.out"; then same=yes; else same=no; fi
  printf '%s\t%s\t%s\t%s\t%s\n' "$workload" "$threads" "$same" "$(stat -f %z "$work/new.out")" "$(stat -f %z "$work/base.out")" >> "$results/identity.tsv"
  rm -f "$work/base.out" "$work/new.out"
  for round in $(seq 1 "$rounds"); do
    echo "$workload $threads: round $round" >&2
    if [ $((round % 2)) = 1 ]; then order="base new"; else order="new base"; fi
    for build in $order; do
      bin="$base"; [ "$build" = new ] && bin="$new"
      t="$(seconds run "$workload" "$bin" "$threads" "$work/$build.out")"
      rm -f "$work/$build.out"
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$workload" "$build" "$threads" "$round" "$t" "$(counter 'instructions retired')" "$(counter 'maximum resident set size')" >> "$results/timings.tsv"
    done
  done
  rm -f "$work/base.out" "$work/new.out"
done

printf 'sidecar_cache_md5_after\t%s\n' "$(md5 -q "$GFF3.fastvep.cache" 2>/dev/null || echo absent)" >> "$results/environment.tsv"

python3 - "$results" <<'PY'
import collections, statistics, sys
results = sys.argv[1]
cells = collections.defaultdict(list)
instructions = collections.defaultdict(list)
for line in list(open(f"{results}/timings.tsv"))[1:]:
    workload, build, threads, _, wall, ins, _rss = line.rstrip("\n").split("\t")
    cells[(workload, threads, build)].append(float(wall))
    instructions[(workload, threads, build)].append(float(ins))
with open(f"{results}/summary.tsv", "w") as out:
    out.write("workload\tthreads\tbuild\truns\tmedian_s\tmin_s\tmax_s\tmedian_instructions_g\n")
    for (workload, threads, build), v in sorted(cells.items()):
        out.write(f"{workload}\t{threads}\t{build}\t{len(v)}\t{statistics.median(v):.2f}\t{min(v):.2f}\t{max(v):.2f}\t{statistics.median(instructions[(workload, threads, build)]) / 1e9:.1f}\n")
print(open(f"{results}/summary.tsv").read())
print(open(f"{results}/identity.tsv").read())
PY
