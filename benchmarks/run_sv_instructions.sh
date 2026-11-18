#!/usr/bin/env bash
# Do two builds of fastVEP do the same work over structural variant VCFs? The instructions that each retires on one thread, in turn,
# the order alternating from one round to the next.
#
#   benchmarks/run_sv_instructions.sh OUT.tsv ROUNDS BINARY_A BINARY_B VCF [VCF...]
#
# Read from the environment: GFF3 and FASTA (as for run_sv_head_to_head.sh), LABEL_A and LABEL_B (default `a` and `b`).
#
# The wall time of a run on a laptop that someone is using moves by 5% or more from one session to the next, and the head-to-head with
# Ensembl VEP (run_sv_head_to_head.sh) takes its times in one session for both tools. A build made after that session answers
# whether its time is the same by the instructions it retires on one thread (`/usr/bin/time -l`), which are the same to 0.1% from
# one run to the next and do not depend on what else the machine is doing. Every run is checked for an output that is not empty.
#
# Writes OUT.tsv: input, build, round, wall seconds, instructions retired, cycles elapsed.
set -euo pipefail
export LC_ALL=C

out="${1:?usage: run_sv_instructions.sh OUT.tsv ROUNDS BINARY_A BINARY_B VCF [VCF...]}"
rounds="${2:?}"
bin_a="${3:?}"
bin_b="${4:?}"
shift 4
[ "$#" -ge 1 ] || { echo "no VCF given" >&2; exit 2; }
: "${GFF3:?set GFF3}" "${FASTA:?set FASTA}"
label_a="${LABEL_A:-a}"
label_b="${LABEL_B:-b}"
scratch="$(mktemp -d)"
trap 'rm -rf "${scratch:?}"' EXIT

printf 'input\tbuild\tround\twall_s\tinstructions\tcycles\n' > "$out"
export RAYON_NUM_THREADS=1
for round in $(seq 1 "$rounds"); do
  for vcf in "$@"; do
    if [ $((round % 2)) = 1 ]; then order=(a b); else order=(b a); fi
    for which in "${order[@]}"; do
      bin="$bin_a"; label="$label_a"
      if [ "$which" = b ]; then bin="$bin_b"; label="$label_b"; fi
      /usr/bin/time -l -o "$scratch/time.txt" "$bin" annotate -i "$vcf" --gff3 "$GFF3" --fasta "$FASTA" \
        --output-format tab --no-progress -o "$scratch/out.tab" > /dev/null 2>&1 || { echo "run failed: $label $vcf" >&2; exit 1; }
      [ -s "$scratch/out.tab" ] || { echo "empty output: $label $vcf" >&2; exit 1; }
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$(basename "$vcf" .vcf)" "$label" "$round" \
        "$(awk '/ real /{print $1}' "$scratch/time.txt")" \
        "$(awk '/instructions retired/{print $1}' "$scratch/time.txt")" \
        "$(awk '/cycles elapsed/{print $1}' "$scratch/time.txt")" >> "$out"
      rm -f "$scratch/out.tab"
    done
  done
done
python3 - "$out" <<'PY'
import collections, csv, statistics, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
by = collections.defaultdict(list)
for r in rows:
    by[(r["input"], r["build"])].append(float(r["instructions"]))
inputs = list(dict.fromkeys(r["input"] for r in rows))
builds = list(dict.fromkeys(r["build"] for r in rows))
for inp in inputs:
    med = {b: statistics.median(by[(inp, b)]) for b in builds}
    print(f"{inp}: " + ", ".join(f"{b} {med[b] / 1e9:.3f} G" for b in builds) + f"; {builds[1]} against {builds[0]}: {100 * (med[builds[1]] / med[builds[0]] - 1):+.3f}%")
PY
