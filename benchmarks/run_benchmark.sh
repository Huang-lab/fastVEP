#!/usr/bin/env bash
#
# OxiVEP Benchmark Suite
# Measures annotation performance across organisms and output formats.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
OXIVEP="$PROJECT_DIR/target/release/oxivep"
OUTPUT_DIR="$SCRIPT_DIR/output"

# Colors for terminal output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

print_header() {
    echo ""
    echo -e "${BOLD}========================================${NC}"
    echo -e "${BOLD}  OxiVEP Benchmark Suite${NC}"
    echo -e "${BOLD}========================================${NC}"
    echo ""
}

print_section() {
    echo ""
    echo -e "${CYAN}--- $1 ---${NC}"
}

# Build release binary
build_release() {
    print_section "Building OxiVEP (release mode)"
    cd "$PROJECT_DIR"
    cargo build --release 2>&1
    if [ ! -f "$OXIVEP" ]; then
        echo -e "${RED}ERROR: Build failed. $OXIVEP not found.${NC}"
        exit 1
    fi
    echo -e "${GREEN}Build successful.${NC}"
}

# Count variants in a VCF (excluding header lines)
count_variants() {
    grep -c -v '^#' "$1"
}

# Run a single benchmark and return elapsed time in seconds (float)
# Usage: run_single <label> <input_vcf> <gff3> <output_format> <output_file> [fasta]
run_single() {
    local label="$1"
    local input_vcf="$2"
    local gff3="$3"
    local fmt="$4"
    local outfile="$5"
    local fasta="${6:-}"

    local fasta_args=""
    if [[ -n "$fasta" && -f "$fasta" ]]; then
        fasta_args="--fasta $fasta"
    fi

    local start end elapsed
    start=$(date +%s%N 2>/dev/null || python3 -c 'import time; print(int(time.time()*1e9))')
    "$OXIVEP" annotate \
        --input "$input_vcf" \
        --gff3 "$gff3" \
        $fasta_args \
        --output "$outfile" \
        --output-format "$fmt" \
        --symbol --hgvs --canonical \
        2>/dev/null || true
    end=$(date +%s%N 2>/dev/null || python3 -c 'import time; print(int(time.time()*1e9))')

    elapsed=$(python3 -c "print(f'{($end - $start) / 1e9:.4f}')")
    echo "$elapsed"
}

# Main benchmark logic
run_benchmarks() {
    mkdir -p "$OUTPUT_DIR"

    # Datasets: name, vcf, gff3
    # Synthetic datasets (small, bundled GFF3)
    declare -a NAMES=("human_chr21" "mouse_chr19" "zebrafish_chr5")
    declare -a VCFS=(
        "$SCRIPT_DIR/human_chr21.vcf"
        "$SCRIPT_DIR/mouse_chr19.vcf"
        "$SCRIPT_DIR/zebrafish_chr5.vcf"
    )
    declare -a GFF3S=(
        "$SCRIPT_DIR/human_chr21.gff3"
        "$SCRIPT_DIR/mouse_chr19.gff3"
        "$SCRIPT_DIR/zebrafish_chr5.gff3"
    )
    declare -a FASTAS=("" "" "")

    # Real-data datasets (require Ensembl GFF3/FASTA in test_data/)
    local test_data="$PROJECT_DIR/test_data"
    if [[ -f "$test_data/chr22.gff3" && -f "$test_data/chr22.fa" ]]; then
        # VEP's own example: 173 chr22 variants with full Ensembl annotations
        NAMES+=("vep_example_chr22")
        VCFS+=("$SCRIPT_DIR/vep_example_chr22.vcf")
        GFF3S+=("$test_data/chr22.gff3")
        FASTAS+=("$test_data/chr22.fa")

        # 1KGP chr22: 1000 real variants with full Ensembl annotations
        if [[ -f "$SCRIPT_DIR/human_chr22_1kgp.vcf" ]]; then
            NAMES+=("human_chr22_1kgp")
            VCFS+=("$SCRIPT_DIR/human_chr22_1kgp.vcf")
            GFF3S+=("$test_data/chr22.gff3")
            FASTAS+=("$test_data/chr22.fa")
        fi
    fi
    declare -a FORMATS=("vcf" "tab" "json")

    # Collect results for summary table
    # results[name_format] = "elapsed variants vps"
    declare -A RESULTS
    declare -A VARIANT_COUNTS

    for i in "${!NAMES[@]}"; do
        local name="${NAMES[$i]}"
        local vcf="${VCFS[$i]}"
        local gff="${GFF3S[$i]}"
        local fasta="${FASTAS[$i]:-}"
        local nvar
        nvar=$(count_variants "$vcf")
        VARIANT_COUNTS["$name"]="$nvar"

        print_section "Benchmarking: $name ($nvar variants)"

        for fmt in "${FORMATS[@]}"; do
            local ext="$fmt"
            [ "$ext" = "tab" ] && ext="tsv"
            local outfile="$OUTPUT_DIR/${name}.annotated.${ext}"
            local label="${name} / ${fmt}"

            printf "  %-35s ... " "$label"
            local elapsed
            elapsed=$(run_single "$label" "$vcf" "$gff" "$fmt" "$outfile" "$fasta")

            local vps
            vps=$(python3 -c "
e = float('$elapsed')
n = int('$nvar')
print(f'{n/e:.1f}' if e > 0 else 'inf')
")
            printf "${GREEN}%8s sec${NC}  (%s variants/sec)\n" "$elapsed" "$vps"
            RESULTS["${name}_${fmt}"]="${elapsed} ${nvar} ${vps}"
        done
    done

    # Print summary table
    echo ""
    echo -e "${BOLD}========================================${NC}"
    echo -e "${BOLD}  Benchmark Summary${NC}"
    echo -e "${BOLD}========================================${NC}"
    echo ""
    printf "${BOLD}%-20s %-10s %10s %10s %15s${NC}\n" "Dataset" "Format" "Variants" "Time (s)" "Variants/sec"
    printf "%-20s %-10s %10s %10s %15s\n"   "--------------------" "----------" "----------" "----------" "---------------"

    for name in "${NAMES[@]}"; do
        for fmt in "${FORMATS[@]}"; do
            local key="${name}_${fmt}"
            if [ -n "${RESULTS[$key]+x}" ]; then
                local parts=(${RESULTS[$key]})
                printf "%-20s %-10s %10s %10s %15s\n" "$name" "$fmt" "${parts[1]}" "${parts[0]}" "${parts[2]}"
            fi
        done
    done

    echo ""
    echo -e "${BOLD}Output files written to:${NC} $OUTPUT_DIR/"
    echo ""

    # Write machine-readable results
    local csv_file="$OUTPUT_DIR/benchmark_results.csv"
    echo "dataset,format,variants,time_seconds,variants_per_second" > "$csv_file"
    for name in "${NAMES[@]}"; do
        for fmt in "${FORMATS[@]}"; do
            local key="${name}_${fmt}"
            if [ -n "${RESULTS[$key]+x}" ]; then
                local parts=(${RESULTS[$key]})
                echo "${name},${fmt},${parts[1]},${parts[0]},${parts[2]}" >> "$csv_file"
            fi
        done
    done
    echo -e "${GREEN}CSV results written to:${NC} $csv_file"
}

# Entry point
print_header

echo "Date:     $(date '+%Y-%m-%d %H:%M:%S')"
echo "System:   $(uname -s) $(uname -m)"
echo "Rust:     $(rustc --version 2>/dev/null || echo 'not found')"
echo ""

build_release
run_benchmarks

echo ""
echo -e "${GREEN}Benchmark complete.${NC}"
