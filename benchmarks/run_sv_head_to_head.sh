#!/usr/bin/env bash
# fastVEP against Ensembl VEP 116.3 on structural variant VCFs: the wall-clock time of annotating
# the whole file, on one thread each, and of fastVEP on every core.
#
#   benchmarks/run_sv_head_to_head.sh RESULTS_DIR VCF [VCF...]
#
# Read from the environment:
#   FASTVEP       the binary under test
#   GFF3          the gene models of the fastVEP side; its binary transcript cache is read when it is
#                 beside the file (the sidecar `<gff3>.fastvep.cache`)
#   FASTA         the reference, with its .fai
#   VEP_GFF3_GZ   the same gene models, bgzipped and tabix-indexed, which Ensembl VEP needs in
#                 `--gff` mode, with FASTA in the same directory (it is mounted as /ref)
#   VEP_IMAGE     the image of Ensembl VEP 116.3, by digest
#
# Each input is run once on each side first and not timed (it reads every file once). Then the
# rounds: 5 of fastVEP (one thread, every core) and 3 of Ensembl VEP, one run at a time, with the
# machine held awake, and nothing else running. The one-thread fastVEP run is `RAYON_NUM_THREADS=1`,
# the same as for the head-to-head of the manuscript (benchmarks/run_head_to_head.py), and VEP is
# one process in Docker, timed from `docker run` to its exit, which includes the start of the
# container, with `--no_stats` as there.
#
# Writes RESULTS_DIR/timings.tsv (every run), environment.tsv (the machine, the versions and the
# SHA-256 of the binary and of every input) and summary.tsv (the medians).
set -euo pipefail
export LC_ALL=C

results="${1:?usage: run_sv_head_to_head.sh RESULTS_DIR VCF [VCF...]}"
shift
[ "$#" -ge 1 ] || { echo "no VCF given" >&2; exit 2; }
: "${FASTVEP:?set FASTVEP to the binary}"
: "${GFF3:?set GFF3}"
: "${FASTA:?set FASTA}"
: "${VEP_GFF3_GZ:?set VEP_GFF3_GZ}"
: "${VEP_IMAGE:?set VEP_IMAGE}"

mkdir -p "$results/out"
results="$(cd "$results" && pwd)"
ref_dir="$(cd "$(dirname "$VEP_GFF3_GZ")" && pwd)"
gff_name="$(basename "$VEP_GFF3_GZ")"
# the reference by its own name, which is the one in the directory that Ensembl VEP is given
fasta_name="$(basename "$(readlink -f "$FASTA")")"
[ -e "$ref_dir/$fasta_name" ] || { echo "$fasta_name is not in $ref_dir, which holds VEP_GFF3_GZ" >&2; exit 2; }
rounds_fastvep=5
rounds_vep=3

# the wall-clock seconds of a command or a function of this script
seconds() {
  local TIMEFORMAT=%R
  { time "$@" >/dev/null 2>&1; } 2>&1
}

run_fastvep() { # THREADS VCF
  local threads="$1" vcf="$2"
  if [ "$threads" = 1 ]; then
    RAYON_NUM_THREADS=1 caffeinate -i "$FASTVEP" annotate -i "$vcf" --gff3 "$GFF3" --fasta "$FASTA" \
      --output-format tab -o "$results/out/fastvep.tab"
  else
    caffeinate -i "$FASTVEP" annotate -i "$vcf" --gff3 "$GFF3" --fasta "$FASTA" \
      --output-format tab -o "$results/out/fastvep.tab"
  fi
}

run_vep() { # VCF
  local vcf="$1"
  caffeinate -i docker run --rm -e VEP_DATABASE=0 \
    -v "$(cd "$(dirname "$vcf")" && pwd):/in:ro" -v "$ref_dir:/ref:ro" -v "$results/out:/out" \
    "$VEP_IMAGE" vep --input_file "/in/$(basename "$vcf")" --gff "/ref/$gff_name" \
    --fasta "/ref/$fasta_name" --output_file /out/vep.tab --tab --no_stats --force_overwrite
}

{
  echo "key	value"
  echo "date	$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "machine	$(sysctl -n machdep.cpu.brand_string), $(sysctl -n hw.ncpu) cores, $(( $(sysctl -n hw.memsize) / 1073741824 )) GiB"
  echo "os	macOS $(sw_vers -productVersion), Darwin $(uname -r)"
  echo "load_average_at_start	$(uptime | sed 's/.*load averages*: //')"
  echo "docker	$(docker info --format '{{.ServerVersion}}, {{.NCPU}} CPUs, {{.MemTotal}} bytes')"
  echo "vep_image	$VEP_IMAGE"
  echo "fastvep_version	$("$FASTVEP" --version)"
  echo "fastvep_sha256	$(shasum -a 256 "$FASTVEP" | cut -d' ' -f1)"
  echo "gff3_sha256	$(shasum -a 256 "$(readlink -f "$GFF3" 2>/dev/null || echo "$GFF3")" | cut -d' ' -f1)"
  for vcf in "$@"; do
    echo "input	$(basename "$vcf")	$(grep -vc '^#' "$vcf") records	sha256 $(shasum -a 256 "$vcf" | cut -d' ' -f1)"
  done
} > "$results/environment.tsv"

printf 'tool\tinput\tthreads\tround\twall_s\n' > "$results/timings.tsv"

for vcf in "$@"; do
  name="$(basename "$vcf")"
  echo "warm-up $name" >&2
  run_fastvep 1 "$vcf" >/dev/null 2>&1
  run_vep "$vcf" >/dev/null 2>&1
done

for round in $(seq 1 "$rounds_fastvep"); do
  for vcf in "$@"; do
    name="$(basename "$vcf")"
    echo "round $round $name" >&2
    t="$(seconds run_fastvep 1 "$vcf")"
    printf 'fastvep\t%s\t1\t%s\t%s\n' "$name" "$round" "$t" >> "$results/timings.tsv"
    t="$(seconds run_fastvep 10 "$vcf")"
    printf 'fastvep\t%s\tall\t%s\t%s\n' "$name" "$round" "$t" >> "$results/timings.tsv"
    if [ "$round" -le "$rounds_vep" ]; then
      t="$(seconds run_vep "$vcf")"
      printf 'vep\t%s\t1\t%s\t%s\n' "$name" "$round" "$t" >> "$results/timings.tsv"
    fi
  done
done

python3 - "$results" <<'PY'
import collections, statistics, sys
results = sys.argv[1]
cells = collections.defaultdict(list)
for line in list(open(f"{results}/timings.tsv"))[1:]:
    tool, name, threads, _, wall = line.rstrip("\n").split("\t")
    cells[(name, tool, threads)].append(float(wall))
with open(f"{results}/summary.tsv", "w") as out:
    out.write("input\ttool\tthreads\truns\tmedian_s\tmin_s\tmax_s\n")
    for (name, tool, threads), v in sorted(cells.items()):
        out.write(f"{name}\t{tool}\t{threads}\t{len(v)}\t{statistics.median(v):.2f}\t{min(v):.2f}\t{max(v):.2f}\n")
print(open(f"{results}/summary.tsv").read())
PY
