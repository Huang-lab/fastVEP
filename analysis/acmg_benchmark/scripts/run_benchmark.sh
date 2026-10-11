#!/usr/bin/env bash
# End-to-end ClinVar 2-star+ ACMG concordance benchmark.
#
# Replaces the run_full_benchmark_v9.sh ... _v14.sh family, which were seven
# byte-identical copies differing only in the version string. Pass the version
# instead:
#
#     bash analysis/acmg_benchmark/scripts/run_benchmark.sh v15
#
# Code and documentation live under analysis/acmg_benchmark/ and are tracked;
# everything this script reads or writes lives under data/benchmark/ and is not.
# That split is the whole reason for the layout - see ../README.md.
#
# Inputs (all under data/benchmark/):
#   clinvar_2star.vcf         ClinVar 2-star+ VCF, ~673k SNV/small-indel
#   clinvar_2star_labels.tsv  ClinVar's 2-star+ classifications, the reference for concordance
#                             (a reference, not ground truth; the file was called
#                             clinvar_2star_truth.tsv, which is still read if the new name is absent)
#   sa_db/                    built .osa/.osa2/.oga databases
#   plus test_data/organisms/human/ for the FASTA + GFF3 + cache
#
# Outputs (data/benchmark/output_<version>/):
#   clinvar_2star.fastvep.vcf.gz   annotation + ACMG, bgzipped
#   concordance_summary.txt        headline metrics
#   concordance_matrix.csv         ClinVar classification x predicted
#   concordance_by_chrom.csv, concordance_by_consequence.csv
#   criterion_firing_rates.csv, rule_distribution.csv
#   discrepancies.tsv              every non-exact call
#   run_provenance.txt             what produced the annotation: the flags, the
#                                  hashes of the binary and of every --acmg-config
#                                  and --functional-evidence, which inputs were
#                                  read (by name and size, not content: the FASTA
#                                  and the SA databases are too large to hash on
#                                  every run), and the commit
#
# Output is VCF-bgzipped (~10x smaller than pretty-printed JSON and ~100x
# smaller than tab format with all fields). VCF is the only format that carries
# ACMG/ACMG_CRITERIA in the per-transcript CSQ entry.
#
# An existing annotated VCF is reused rather than regenerated, so re-running to
# recompute concordance after a scoring change is cheap - but only when
# run_provenance.txt says the same binary, flags, configuration and inputs made
# it. Anything else stops the run (exit 3): delete the directory, or give the
# run its own version name.
#
# Variants of a run, under their own version names:
#   ACMG_CONFIG=path.toml  an extra --acmg-config
#   EXTRA_ARGS="..."       any other `fastvep annotate` flags, split on whitespace
#   LABELS=path.tsv        the ClinVar 2-star+ labels (TRUTH is the old name)

set -euo pipefail

VERSION="${1:-}"
if [ -z "$VERSION" ]; then
  echo "usage: $(basename "$0") <version>    e.g. $(basename "$0") v15" >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

INPUT="${INPUT:-$ROOT/data/benchmark/clinvar_2star.vcf}"
# TRUTH is the old name of this variable and is still honoured. The old file name
# is the fallback only when neither variable is set: a path that was asked for and
# is missing is an error, and silently reading a different file instead would
# score against the wrong labels.
LABELS="${LABELS:-${TRUTH:-}}"
if [ -z "$LABELS" ]; then
  LABELS="$ROOT/data/benchmark/clinvar_2star_labels.tsv"
  if [ ! -e "$LABELS" ] && [ -e "$ROOT/data/benchmark/clinvar_2star_truth.tsv" ]; then
    LABELS="$ROOT/data/benchmark/clinvar_2star_truth.tsv"
  fi
fi
GFF3="${GFF3:-$ROOT/test_data/organisms/human/Homo_sapiens.GRCh38.115.gff3}"
FASTA="${FASTA:-$ROOT/test_data/organisms/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa}"
SA_DIR="${SA_DIR:-$ROOT/data/benchmark/sa_db}"
FASTVEP="${FASTVEP:-$ROOT/target/release/fastvep}"
OUT_DIR="${OUT_DIR:-$ROOT/data/benchmark/output_$VERSION}"

# An extra --acmg-config, to run a variant of the same version.
#
# Expanded below as ${CONFIG_ARGS[@]+...} rather than plain "${CONFIG_ARGS[@]}":
# under `set -u`, bash 3.2 - which is what macOS ships - treats an empty array
# as an unbound variable and aborts. An `if` rather than `[ ... ] && ...` for
# the same class of reason: a failing test as the last command in an && list
# is a `set -e` trap.
ACMG_CONFIG="${ACMG_CONFIG:-}"
CONFIG_ARGS=()
if [ -n "$ACMG_CONFIG" ]; then
  CONFIG_ARGS=(--acmg-config "$ACMG_CONFIG")
fi

# Any other `fastvep annotate` flags, for a variant of the same version - for
# instance `EXTRA_ARGS="--pick-order rank,mane_select"`. Split on whitespace,
# newlines included, and never globbed, so a flag and its value stay two words.
EXTRA_ARGS="${EXTRA_ARGS:-}"
EXTRA_ARR=()
if [ -n "$EXTRA_ARGS" ]; then
  set -f
  # shellcheck disable=SC2206 # splitting is the point
  EXTRA_ARR=($EXTRA_ARGS)
  set +f
fi

for f in "$INPUT" "$LABELS" "$GFF3" "$FASTA" "$FASTVEP"; do
  [ -e "$f" ] || { echo "missing required input: $f" >&2; exit 1; }
done
[ -d "$SA_DIR" ] || { echo "missing SA database directory: $SA_DIR" >&2; exit 1; }
# `sa_db/` and `sa_db` are one directory, and the signature below strips this
# prefix from every path `find` prints.
SA_DIR="${SA_DIR%/}"

mkdir -p "$OUT_DIR"
VCFGZ="$OUT_DIR/clinvar_2star.fastvep.vcf.gz"
PROVENANCE="$OUT_DIR/run_provenance.txt"

# ---- What produced the annotation, recorded beside it --------------------------
#
# The version name is the only part of the output's key, so without this a second
# run under the same name with a different config, flags, binary or inputs would
# find the first run's annotation and score it as if it were its own. What is
# compared is the content that decides the annotation - the flags, the config
# files' hashes, the binary's hash, and which inputs were read - and not how a
# path happened to be spelled.
sha256() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum | cut -d' ' -f1; else shasum -a 256 | cut -d' ' -f1; fi
}
abspath() {
  if [ -d "$1" ]; then (cd "$1" && pwd -P); else printf '%s/%s' "$(cd "$(dirname "$1")" && pwd -P)" "$(basename "$1")"; fi
}
bytes() { wc -c < "$1" | tr -d ' '; }

# Every --acmg-config the run reads, whether it came in through ACMG_CONFIG or
# through EXTRA_ARGS, so editing either one is seen.
CONFIG_FILES=()
if [ -n "$ACMG_CONFIG" ]; then
  CONFIG_FILES+=("$ACMG_CONFIG")
fi
prev=""
for a in ${EXTRA_ARR[@]+"${EXTRA_ARR[@]}"}; do
  if [ "$prev" = "--acmg-config" ]; then
    CONFIG_FILES+=("$a")
  fi
  case "$a" in
    --acmg-config=*) CONFIG_FILES+=("${a#--acmg-config=}") ;;
  esac
  prev="$a"
done
CONFIG_HASHES=""
for c in ${CONFIG_FILES[@]+"${CONFIG_FILES[@]}"}; do
  [ -f "$c" ] || { echo "missing --acmg-config file: $c" >&2; exit 1; }
  CONFIG_HASHES="${CONFIG_HASHES:+$CONFIG_HASHES,}$(sha256 < "$c")"
done
CONFIG_HASHES="${CONFIG_HASHES:-none}"

# The flags, with each config's path replaced by a placeholder: its hash says
# what was in it. A `--functional-evidence` file supplies PS3 and BS3, so it is
# read for its content the same way, and its path is not part of the key.
ARGS_NORM="annotate --acmg --pick"
prev=""
for a in ${CONFIG_ARGS[@]+"${CONFIG_ARGS[@]}"} ${EXTRA_ARR[@]+"${EXTRA_ARR[@]}"}; do
  if [ "$prev" = "--acmg-config" ]; then
    ARGS_NORM="$ARGS_NORM <config>"
  elif [ "$prev" = "--functional-evidence" ] && [ -f "$a" ]; then
    ARGS_NORM="$ARGS_NORM <functional-evidence:$(sha256 < "$a")>"
  else
    case "$a" in
      --acmg-config=*) ARGS_NORM="$ARGS_NORM --acmg-config=<config>" ;;
      --functional-evidence=*)
        fe="${a#--functional-evidence=}"
        if [ -f "$fe" ]; then
          ARGS_NORM="$ARGS_NORM --functional-evidence=<functional-evidence:$(sha256 < "$fe")>"
        else
          ARGS_NORM="$ARGS_NORM $a"
        fi
        ;;
      *) ARGS_NORM="$ARGS_NORM $a" ;;
    esac
  fi
  prev="$a"
done

# Each file the loader opens, by name under the directory and by size, in a fixed
# order. `-L` follows a symlinked directory or file, which `find` would otherwise
# read as no files at all and sign as the hash of nothing; the extensions are the
# ones `--sa-dir` reads (and `.idx`, the index of a `.osa`), so a `.DS_Store` or an
# editor's backup is not an input.
SA_SIGNATURE="$(find -L "$SA_DIR" -type f \( -name '*.osa' -o -name '*.osa2' -o -name '*.osi' -o -name '*.oga' -o -name '*.idx' \) | LC_ALL=C sort | while read -r f; do printf '%s %s\n' "${f#"$SA_DIR"/}" "$(bytes "$f")"; done | sha256)"
# Read now and not after the annotation, so that a commit made while it runs is
# not recorded as the one that made it. `unknown` and not a number when git cannot
# say: a count of zero would read as a clean tree.
if GIT_COMMIT="$(git -C "$ROOT" rev-parse HEAD 2>/dev/null)" \
  && GIT_DIRTY="$(git -C "$ROOT" status --porcelain 2>/dev/null | wc -l | tr -d ' ')"; then
  :
else
  GIT_COMMIT="unknown"
  GIT_DIRTY="unknown"
fi
CURRENT="provenance_version: 2
args: $ARGS_NORM
acmg_config_sha256: $CONFIG_HASHES
fastvep_sha256: $(sha256 < "$FASTVEP")
input: $(abspath "$INPUT")
input_bytes: $(bytes "$INPUT")
gff3: $(abspath "$GFF3")
gff3_bytes: $(bytes "$GFF3")
fasta: $(abspath "$FASTA")
fasta_bytes: $(bytes "$FASTA")
sa_dir: $(abspath "$SA_DIR")
sa_dir_signature: $SA_SIGNATURE"

# A bgzip with a header and no records, or with nothing in it at all, passes `-s`
# and used to be reused and scored as zero classified variants.
have_records() { [ -s "$1" ] && [ -n "$(gzip -dc "$1" 2>/dev/null | grep -v '^#' | head -n 1)" ]; }
HAVE_ANNOTATION=0
if have_records "$VCFGZ"; then HAVE_ANNOTATION=1; fi

# One line per key, whitespace-trimmed, with a recorded config path folded into
# its placeholder the way the current one is.
trim() { sed -E 's/[[:space:]]+$//; s/[[:space:]]+/ /g; s/--acmg-config [^ ]+/--acmg-config <config>/'; }
recorded_keys() { grep -E "^($1):" "$PROVENANCE" | trim || true; }

if [ "$HAVE_ANNOTATION" = 1 ] && [ -f "$PROVENANCE" ]; then
  if grep -q '^provenance_version: 2$' "$PROVENANCE"; then
    KEYS='provenance_version|args|acmg_config_sha256|fastvep_sha256|input|input_bytes|gff3|gff3_bytes|fasta|fasta_bytes|sa_dir|sa_dir_signature'
    WANT="$CURRENT"
  else
    # Written by an earlier version of this script, which recorded the flags and
    # the config's hash and nothing else.
    KEYS='args|acmg_config_sha256'
    WANT="$(printf '%s\n' "$CURRENT" | grep -E '^(args|acmg_config_sha256):')"
  fi
  RECORDED="$(recorded_keys "$KEYS")"
  if [ "$RECORDED" != "$(printf '%s\n' "$WANT" | trim)" ]; then
    echo "error: $VCFGZ was produced by something other than what this run asks for." >&2
    echo "--- recorded ($PROVENANCE)" >&2
    echo "$RECORDED" >&2
    echo "--- requested" >&2
    printf '%s\n' "$WANT" | trim >&2
    echo "Delete $OUT_DIR, or give this run its own version name." >&2
    exit 3
  fi
elif [ "$HAVE_ANNOTATION" = 1 ] && { [ -n "$ACMG_CONFIG" ] || [ -n "$EXTRA_ARGS" ]; }; then
  # An annotation from before provenance was recorded, and a run that asks for
  # something other than the defaults: there is nothing to check it against.
  echo "error: $VCFGZ predates run_provenance.txt, so it cannot be checked against ACMG_CONFIG or EXTRA_ARGS." >&2
  echo "Delete $OUT_DIR, or give this run its own version name." >&2
  exit 3
fi

if [ "$HAVE_ANNOTATION" = 0 ]; then
  echo "==> Annotating $(grep -vc '^#' "$INPUT") variants with --acmg ($VERSION)..."
  T0=$(date +%s)
  # Written under a temporary name and moved into place only once the whole
  # pipeline has succeeded. A failed annotate - a mistyped flag, say - used to
  # leave a valid, empty bgzip under the final name, which the next run found,
  # reused and scored as zero classified variants.
  rm -f "$VCFGZ.partial"
  "$FASTVEP" annotate \
    -i "$INPUT" \
    -o - \
    --gff3 "$GFF3" \
    --fasta "$FASTA" \
    --sa-dir "$SA_DIR" \
    --acmg --pick \
    ${CONFIG_ARGS[@]+"${CONFIG_ARGS[@]}"} \
    ${EXTRA_ARR[@]+"${EXTRA_ARR[@]}"} \
    --output-format vcf \
    | bgzip -c > "$VCFGZ.partial"
  T1=$(date +%s)
  echo "==> Annotation took $((T1-T0)) seconds"
  if ! have_records "$VCFGZ.partial"; then
    echo "error: fastvep wrote no variant records to $VCFGZ.partial, so there is nothing to record or score." >&2
    exit 1
  fi
  {
    echo "$CURRENT"
    echo "fastvep: $FASTVEP"
    echo "fastvep_version: $("$FASTVEP" --version 2>&1 | head -1)"
    echo "repo_commit: $GIT_COMMIT"
    echo "repo_dirty_files: $GIT_DIRTY"
    echo "annotation_seconds: $((T1-T0))"
    echo "recorded_at: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  } > "$PROVENANCE.partial"
  # The record goes first. Interrupted between the two moves, the directory holds
  # a record and no annotation, which the next run simply redoes; the other order
  # would leave an annotation with no record of what made it.
  mv "$PROVENANCE.partial" "$PROVENANCE"
  mv "$VCFGZ.partial" "$VCFGZ"
elif [ ! -f "$PROVENANCE" ]; then
  echo "==> Reusing existing annotation at $VCFGZ, skipping. It predates run_provenance.txt, so what produced it is not recorded."
else
  echo "==> Reusing existing annotation at $VCFGZ, skipping. Delete it to re-annotate."
fi
ls -la "$VCFGZ"

echo "==> Computing concordance..."
python3 "$SCRIPT_DIR/03_evaluate_concordance.py" \
  --labels "$LABELS" \
  --predictions "$VCFGZ" \
  --out "$OUT_DIR"

echo "==> Done. See $OUT_DIR/"
ls -la "$OUT_DIR/"
echo
echo "Record what changed in analysis/acmg_benchmark/RUN_VERSIONS.md."
