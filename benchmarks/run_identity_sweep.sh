#!/usr/bin/env bash
# Do two or more builds of fastVEP write the same bytes, and say the same on stderr, on the inputs and with the flags of
# every benchmark of this repository that is not about structural variants?
#
#   benchmarks/run_identity_sweep.sh RESULTS_DIR label=BINARY [label=BINARY ...]
#
# The first label is the base, and each case is compared with it. CASES='regex' runs only the cases whose name matches.
#
# What it runs (the inputs of benchmarks/run_head_to_head.py, benchmarks/run_benchmark.sh and
# analysis/acmg_benchmark/scripts/run_benchmark.sh, and the flags of each):
#   hg002_*        the whole of GIAB HG002 as VCF with `--hgvs --symbol --canonical` plain, with ClinVar, with ClinVar and the
#                  v4.1 exomes of gnomAD, and with ClinVar and the joint release; as VCF with no FASTA (consequences only); as the
#                  table of the comparison with Ensembl VEP 116.3
#   hg002_json     the whole of HG002 as JSON (the writer of a record is its own, and a structural variant's `Allele` and
#                  `total_length` are written there)
#   chr22_*        chromosome 22 of HG002 as VCF, tab and JSON, and with each of the pick options
#   <organism>_full  the whole variant file of the yeast, drosophila, arabidopsis and mouse benchmarks, with `--hgvs`
#   clinvar_*      the three ClinVar 2-star+ sets (April 2026, October 2026, October 2026 covered), and a 1-in-5 sample with `--flag-pick`
#   acmg_*         the three ClinVar sets with `--acmg --pick` over the databases of the ACMG benchmark
#
# Each build reads a transcript cache of its own (links to the GFF3 and the FASTA in a directory per label, so the sidecar
# cache it writes is its own, which also compares the caches: they are hashed at the end). Output goes through a pipe into a
# hash, so nothing of the 10 GB outputs is kept. A run on every core, at low priority.
#
# Writes RESULTS_DIR/identity.tsv (one row per case and label: exit status, bytes, MD5 of stdout, MD5 of stderr with the numbers
# of times and sizes and the label's directory taken out, seconds), summary.tsv (one row per case: whether every label wrote the
# base's bytes and its warnings) and environment.tsv (the binaries, the inputs, the caches).
#
# Expected: every case identical between two builds that differ in nothing a small variant meets. Two builds differ in
# what they do meet:
#   - one that flags the `intergenic_variant` row of a variant that no transcript is beside (`--flag-pick*`, the two `*flagpick*`
#     cases) differs from one that does not in exactly those rows;
#   - one that ranks `splice_donor_5th_base_variant` before `splice_region_variant`, as release 116.3 does, differs from one that
#     ranks them the other way round in the `most_severe_consequence` of a JSON record (`hg002_json`) whose transcripts carry
#     one term each, and in nothing else: a row's terms are in rank order, and a row with the one has not the other.
# The paper's `fastvep_v38` of the ACMG benchmark differs from a later build in the `acmg_*` cases (the per-variant coverage
# rule of PM2: 550 mitochondrial variants of the 673,660 of April).
set -u
export LC_ALL=C

out="${1:?usage: run_identity_sweep.sh RESULTS_DIR label=BINARY [label=BINARY ...]}"; shift
[ "$#" -ge 2 ] || { echo "name at least two builds, label=BINARY" >&2; exit 2; }
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
human="${HUMAN_DIR:-$root/test_data/organisms/human}"
organisms="${ORGANISMS_DIR:-$root/test_data/organisms}"
bench="${BENCHMARK_DIR:-$root/data/benchmark}"
hg002="${HG002_VCF:-$human/human_giab_hg002_full.vcf}"
sa_db="${SA_DB:-$bench/sa_db}" sa_joint="${SA_DB_JOINT:-$bench/sa_db_joint}" sa_complete="${SA_DB_COMPLETE:-$bench/sa_db_complete}"
cv_april="${CLINVAR_2STAR:-$bench/clinvar_2star.vcf}"
cv_oct="${CLINVAR_OCT:-$bench/clinvar_20261004/clinvar_2star.vcf}"
cv_octcov="${CLINVAR_OCT_COVERED:-$bench/clinvar_20261004_covered/clinvar_2star.vcf}"
mkdir -p "$out/res" "$out/sa" "$out/in"
out="$(cd "$out" && pwd)"
if command -v md5 >/dev/null 2>&1; then md5_of() { md5 -q "$@"; }; else md5_of() { md5sum "$@" | cut -d' ' -f1; }; fi

for need in "$hg002" "$human/Homo_sapiens.GRCh38.115.gff3" "$human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"; do
  [ -e "$need" ] || { echo "missing $need" >&2; exit 2; }
done

# inputs that have to be cut
[ -s "$out/in/chr22.vcf" ] || awk -F'\t' '/^#/ || $1 == "22" || $1 == "chr22"' "$hg002" > "$out/in/chr22.vcf"
[ -s "$out/in/clinvar_1in5.vcf" ] || awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} (++n % 5 == 0)' "$cv_april" > "$out/in/clinvar_1in5.vcf"
# the annotation conditions of the head-to-head, as links (run_head_to_head.py makes the same directories)
for cond in clinvar both joint; do rm -rf "$out/sa/$cond"; mkdir -p "$out/sa/$cond"; done
for f in clinvar.osa clinvar.osa.idx; do for cond in clinvar both joint; do ln -s "$sa_db/$f" "$out/sa/$cond/$f"; done; done
for f in "$sa_db"/gnomad_chr*.osa2; do ln -s "$f" "$out/sa/both/$(basename "$f")"; done
for f in "$sa_joint"/gnomad_chr*.osa2; do ln -s "$f" "$out/sa/joint/$(basename "$f")"; done

labels=(); bins=()
for pair in "$@"; do labels+=("${pair%%=*}"); bins+=("${pair#*=}"); done

gene_model() { # LABEL NAME GFF3 FASTA: a directory of links of this label's own, so that its sidecar cache is its own
  local d="$out/gm/$1/$2"; mkdir -p "$d"
  [ -e "$d/m.gff3" ] || ln -s "$3" "$d/m.gff3"
  [ -e "$d/m.fa" ] || ln -s "$4" "$d/m.fa"
  [ -e "$d/m.fa.fai" ] || { [ -e "$4.fai" ] && ln -s "$4.fai" "$d/m.fa.fai"; }
  echo "$d"
}

hashing='import hashlib, sys
h, n = hashlib.md5(), 0
for block in iter(lambda: sys.stdin.buffer.read(1 << 24), b""):
    h.update(block); n += len(block)
print(n, h.hexdigest())'

# name | input | gene model | arguments
cases=(
"hg002_plain|$hg002|human|--hgvs --symbol --canonical --output-format vcf"
"hg002_clinvar|$hg002|human|--hgvs --symbol --canonical --output-format vcf --sa-dir $out/sa/clinvar"
"hg002_exomes|$hg002|human|--hgvs --symbol --canonical --output-format vcf --sa-dir $out/sa/both"
"hg002_joint|$hg002|human|--hgvs --symbol --canonical --output-format vcf --sa-dir $out/sa/joint"
"hg002_conseq_only|$hg002|human_nofasta|--output-format vcf"
"hg002_tab116|$hg002|human|--output-format tab --symbol --biotype --canonical --numbers --variant-class --show-ref-allele --hgvs --hgvsg"
"hg002_json|$hg002|human|--hgvs --symbol --canonical --output-format json"
"chr22_vcf|$out/in/chr22.vcf|human|--hgvs --symbol --canonical --output-format vcf"
"chr22_tab|$out/in/chr22.vcf|human|--hgvs --symbol --canonical --output-format tab"
"chr22_json|$out/in/chr22.vcf|human|--hgvs --symbol --canonical --output-format json"
"chr22_pick|$out/in/chr22.vcf|human|--hgvs --symbol --canonical --pick --output-format tab"
"chr22_pick_allele|$out/in/chr22.vcf|human|--hgvs --symbol --canonical --pick-allele --output-format tab"
"chr22_flagpick|$out/in/chr22.vcf|human|--hgvs --symbol --canonical --flag-pick --output-format tab"
"chr22_flagpick_allele|$out/in/chr22.vcf|human|--hgvs --symbol --canonical --flag-pick-allele --output-format tab"
"yeast_full|$organisms/yeast/yeast_ensembl_full.vcf|yeast|--hgvs --output-format vcf"
"drosophila_full|$organisms/drosophila/drosophila_dgrp2_full.vcf|drosophila|--hgvs --output-format vcf"
"arabidopsis_full|$organisms/arabidopsis/arabidopsis_1001g_full.vcf|arabidopsis|--hgvs --output-format vcf"
"mouse_full|$organisms/mouse/mouse_cast_eij_full.vcf|mouse|--hgvs --output-format vcf"
"clinvar_april_plain|$cv_april|human|--hgvs --symbol --canonical --output-format vcf"
"clinvar_oct_plain|$cv_oct|human|--hgvs --symbol --canonical --output-format vcf"
"clinvar_octcov_plain|$cv_octcov|human|--hgvs --symbol --canonical --output-format vcf"
"clinvar_1in5_flagpick|$out/in/clinvar_1in5.vcf|human|--hgvs --symbol --canonical --flag-pick --output-format tab"
"acmg_april|$cv_april|human|--sa-dir $sa_complete --acmg --pick --output-format vcf"
"acmg_oct|$cv_oct|human|--sa-dir $sa_complete --acmg --pick --output-format vcf"
"acmg_octcov|$cv_octcov|human|--sa-dir $sa_complete --acmg --pick --output-format vcf"
)

printf 'case\tlabel\texit\tbytes\tmd5_stdout\tmd5_stderr_normalised\tseconds\n' > "$out/identity.tsv"
for entry in "${cases[@]}"; do
  IFS='|' read -r name input model args <<< "$entry"
  if [ -n "${CASES:-}" ] && ! [[ "$name" =~ $CASES ]]; then continue; fi
  if [ ! -e "$input" ]; then echo "$name: no input $input, skipped" >&2; continue; fi
  for i in "${!labels[@]}"; do
    label="${labels[$i]}"; bin="${bins[$i]}"
    case "$model" in
      human)         d="$(gene_model "$label" human "$human/Homo_sapiens.GRCh38.115.gff3" "$human/Homo_sapiens.GRCh38.dna.primary_assembly.fa")"; fasta="--fasta $d/m.fa";;
      human_nofasta) d="$(gene_model "$label" human_nofasta "$human/Homo_sapiens.GRCh38.115.gff3" "$human/Homo_sapiens.GRCh38.dna.primary_assembly.fa")"; fasta="";;
      *)             d="$(gene_model "$label" "$model" "$organisms/$model/$model.gff3" "$organisms/$model/$model.fa")"; fasta="--fasta $d/m.fa";;
    esac
    t0=$(date +%s)
    # shellcheck disable=SC2086
    result="$(nice -n 10 "$bin" annotate -i "$input" --gff3 "$d/m.gff3" $fasta $args --no-progress -o - 2> "$out/res/$name.$label.err" | python3 -c "$hashing"; exit "${PIPESTATUS[0]}")"
    status=$?
    t1=$(date +%s)
    [ -n "$result" ] || result="0 none"
    # warnings with the label's own directory, the times and the sizes taken out
    # (`\b` is not a word boundary in the `sed -E` of macOS, so the end of a number with its unit is "followed by something that is not a letter, or by the end of the line")
    norm="$(sed -E 's#/gm/[^/]+/#/gm/L/#g; s/[0-9]+(\.[0-9]+)? ?(s|ms|seconds|MB|MiB|GB|GiB|KB)([^A-Za-z0-9_]|$)/N\3/g; s/[0-9]+(\.[0-9]+)?x([^A-Za-z0-9_]|$)/Nx\2/g' "$out/res/$name.$label.err" | grep -v '^$' | md5_of)"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$label" "$status" "${result% *}" "${result#* }" "$norm" "$((t1 - t0))" >> "$out/identity.tsv"
    echo "$name $label exit=$status bytes=${result% *} md5=${result#* } ($((t1 - t0)) s)" >&2
  done
done

python3 - "$out" "${labels[0]}" <<'PY'
import collections, csv, sys
out, base = sys.argv[1], sys.argv[2]
cases = collections.OrderedDict()
for r in csv.DictReader(open(f"{out}/identity.tsv"), delimiter="\t"):
    cases.setdefault(r["case"], {})[r["label"]] = r
with open(f"{out}/summary.tsv", "w") as f:
    labels = list(next(iter(cases.values())))
    f.write("case\tbytes_" + base + "\t" + "\t".join(f"same_bytes_as_{base}:{l}" for l in labels if l != base) + "\t" + "\t".join(f"same_stderr_as_{base}:{l}" for l in labels if l != base) + "\n")
    for case, d in cases.items():
        b = d[base]
        same = ["yes" if d[l]["md5_stdout"] == b["md5_stdout"] and d[l]["exit"] == b["exit"] else "NO" for l in labels if l != base]
        err = ["yes" if d[l]["md5_stderr_normalised"] == b["md5_stderr_normalised"] else "NO" for l in labels if l != base]
        f.write(case + "\t" + b["bytes"] + "\t" + "\t".join(same) + "\t" + "\t".join(err) + "\n")
print(open(f"{out}/summary.tsv").read())
PY

{
  printf 'started_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  for i in "${!labels[@]}"; do printf 'binary_%s\t%s\t%s\tsha256 %s\n' "${labels[$i]}" "$(basename "${bins[$i]}")" "$("${bins[$i]}" --version 2>&1 | head -1)" "$(shasum -a 256 "${bins[$i]}" | cut -d' ' -f1)"; done
  for i in "${!labels[@]}"; do
    for cache in "$out"/gm/"${labels[$i]}"/*/m.gff3.fastvep.cache; do
      [ -e "$cache" ] && printf 'cache_%s_%s\t%s bytes\tmd5 %s\n' "${labels[$i]}" "$(basename "$(dirname "$cache")")" "$(stat -f %z "$cache" 2>/dev/null || stat -c %s "$cache")" "$(md5_of "$cache")"
    done
  done
  printf 'hg002\t%s\n' "$(basename "$hg002")"
} > "$out/environment.tsv"
