#!/usr/bin/env bash
#
# Download gold-standard benchmark data for OxiVEP
#
# Downloads reference genomes (GFF3 + FASTA) and gold-standard variant call sets
# from authoritative sources for each model organism.
#
# Usage: ./download_data.sh [--all | --human | --mouse | --yeast | ...]
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
TEST_DATA="$PROJECT_DIR/test_data"
ORG_DATA="$TEST_DATA/organisms"
VCF_DATA="$TEST_DATA/benchmark_vcfs"

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

mkdir -p "$TEST_DATA" "$ORG_DATA" "$VCF_DATA"

download() {
    local url="$1" dest="$2"
    if [[ -f "$dest" ]]; then
        echo -e "  ${GREEN}Already exists:${NC} $(basename "$dest")"
        return 0
    fi
    echo -e "  ${YELLOW}Downloading:${NC} $(basename "$dest")"
    curl -L --progress-bar -o "$dest" "$url"
}

decompress() {
    local gz="$1"
    local out="${gz%.gz}"
    if [[ -f "$out" ]]; then
        echo -e "  ${GREEN}Already decompressed:${NC} $(basename "$out")"
        return 0
    fi
    echo -e "  Decompressing $(basename "$gz")..."
    gunzip -k "$gz"
}

index_fasta() {
    local fa="$1"
    if [[ -f "${fa}.fai" ]]; then return 0; fi
    if command -v samtools &>/dev/null; then
        echo -e "  Indexing $(basename "$fa")..."
        samtools faidx "$fa"
    else
        echo -e "  ${YELLOW}samtools not found — skipping FASTA index (install for mmap support)${NC}"
    fi
}

# ═══════════════════════════════════════════════════════════════
# Reference data: GFF3 + FASTA per organism
# ═══════════════════════════════════════════════════════════════
ENSEMBL_RELEASE=115
ENSEMBL_FTP="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}"

download_human_ref() {
    echo -e "\n${BOLD}Human (GRCh38, Ensembl ${ENSEMBL_RELEASE})${NC}"
    download "${ENSEMBL_FTP}/gff3/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gff3.gz" \
        "$TEST_DATA/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gff3.gz"
    decompress "$TEST_DATA/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gff3.gz"

    download "${ENSEMBL_FTP}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
        "$TEST_DATA/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    decompress "$TEST_DATA/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    index_fasta "$TEST_DATA/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
}

download_mouse_ref() {
    echo -e "\n${BOLD}Mouse (GRCm39, Ensembl ${ENSEMBL_RELEASE})${NC}"
    download "${ENSEMBL_FTP}/gff3/mus_musculus/Mus_musculus.GRCm39.${ENSEMBL_RELEASE}.gff3.gz" \
        "$ORG_DATA/mouse.gff3.gz"
    decompress "$ORG_DATA/mouse.gff3.gz"
    # Rename for consistency
    [[ -f "$ORG_DATA/mouse.gff3" ]] || mv "$ORG_DATA/Mus_musculus.GRCm39.${ENSEMBL_RELEASE}.gff3" "$ORG_DATA/mouse.gff3" 2>/dev/null || true

    download "${ENSEMBL_FTP}/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz" \
        "$ORG_DATA/mouse.fa.gz"
    decompress "$ORG_DATA/mouse.fa.gz"
    index_fasta "$ORG_DATA/mouse.fa"
}

download_yeast_ref() {
    echo -e "\n${BOLD}Yeast (R64, Ensembl ${ENSEMBL_RELEASE})${NC}"
    download "${ENSEMBL_FTP}/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.${ENSEMBL_RELEASE}.gff3.gz" \
        "$ORG_DATA/yeast.gff3.gz"
    decompress "$ORG_DATA/yeast.gff3.gz"
    [[ -f "$ORG_DATA/yeast.gff3" ]] || mv "$ORG_DATA/Saccharomyces_cerevisiae.R64-1-1.${ENSEMBL_RELEASE}.gff3" "$ORG_DATA/yeast.gff3" 2>/dev/null || true

    download "${ENSEMBL_FTP}/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz" \
        "$ORG_DATA/yeast.fa.gz"
    decompress "$ORG_DATA/yeast.fa.gz"
    index_fasta "$ORG_DATA/yeast.fa"
}

download_drosophila_ref() {
    echo -e "\n${BOLD}Drosophila (BDGP6, Ensembl ${ENSEMBL_RELEASE})${NC}"
    download "${ENSEMBL_FTP}/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.${ENSEMBL_RELEASE}.gff3.gz" \
        "$ORG_DATA/drosophila.gff3.gz"
    decompress "$ORG_DATA/drosophila.gff3.gz"
    [[ -f "$ORG_DATA/drosophila.gff3" ]] || mv "$ORG_DATA/Drosophila_melanogaster.BDGP6.46.${ENSEMBL_RELEASE}.gff3" "$ORG_DATA/drosophila.gff3" 2>/dev/null || true

    download "${ENSEMBL_FTP}/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz" \
        "$ORG_DATA/drosophila.fa.gz"
    decompress "$ORG_DATA/drosophila.fa.gz"
    index_fasta "$ORG_DATA/drosophila.fa"
}

download_elegans_ref() {
    echo -e "\n${BOLD}C. elegans (WBcel235, Ensembl ${ENSEMBL_RELEASE})${NC}"
    download "${ENSEMBL_FTP}/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.${ENSEMBL_RELEASE}.gff3.gz" \
        "$ORG_DATA/elegans.gff3.gz"
    decompress "$ORG_DATA/elegans.gff3.gz"
    [[ -f "$ORG_DATA/elegans.gff3" ]] || mv "$ORG_DATA/Caenorhabditis_elegans.WBcel235.${ENSEMBL_RELEASE}.gff3" "$ORG_DATA/elegans.gff3" 2>/dev/null || true

    download "${ENSEMBL_FTP}/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz" \
        "$ORG_DATA/elegans.fa.gz"
    decompress "$ORG_DATA/elegans.fa.gz"
    index_fasta "$ORG_DATA/elegans.fa"
}

download_arabidopsis_ref() {
    echo -e "\n${BOLD}Arabidopsis (TAIR10, Ensembl Plants ${ENSEMBL_RELEASE})${NC}"
    local PLANTS_FTP="https://ftp.ensemblgenomes.org/pub/plants/release-${ENSEMBL_RELEASE}"
    download "${PLANTS_FTP}/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.${ENSEMBL_RELEASE}.gff3.gz" \
        "$ORG_DATA/arabidopsis.gff3.gz"
    decompress "$ORG_DATA/arabidopsis.gff3.gz"
    [[ -f "$ORG_DATA/arabidopsis.gff3" ]] || mv "$ORG_DATA/Arabidopsis_thaliana.TAIR10.${ENSEMBL_RELEASE}.gff3" "$ORG_DATA/arabidopsis.gff3" 2>/dev/null || true

    download "${PLANTS_FTP}/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz" \
        "$ORG_DATA/arabidopsis.fa.gz"
    decompress "$ORG_DATA/arabidopsis.fa.gz"
    index_fasta "$ORG_DATA/arabidopsis.fa"
}

# ═══════════════════════════════════════════════════════════════
# Gold-standard benchmark VCFs
# ═══════════════════════════════════════════════════════════════

download_human_vcf() {
    echo -e "\n${BOLD}Human VCF: GIAB HG002 (GRCh38)${NC}"
    download "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
        "$VCF_DATA/human_giab_HG002.vcf.gz"
    download "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" \
        "$VCF_DATA/human_giab_HG002.vcf.gz.tbi"
}

download_mouse_vcf() {
    echo -e "\n${BOLD}Mouse VCF: Mouse Genomes Project (GRCm39)${NC}"
    download "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz" \
        "$VCF_DATA/mouse_mgp.vcf.gz"
}

download_yeast_vcf() {
    echo -e "\n${BOLD}Yeast VCF: 1002 Yeast Genomes Project (R64)${NC}"
    download "http://1002genomes.u-strasbg.fr/files/1011Matrix.gvcf.gz" \
        "$VCF_DATA/yeast_1002g.gvcf.gz"
}

download_drosophila_vcf() {
    echo -e "\n${BOLD}Drosophila VCF: DGRP2 (dm6)${NC}"
    echo -e "  ${YELLOW}Note: DGRP2 uses dm6 coordinates. For Ensembl BDGP6 benchmarking,${NC}"
    echo -e "  ${YELLOW}use Ensembl variation VCF instead:${NC}"
    download "${ENSEMBL_FTP}/variation/vcf/drosophila_melanogaster/drosophila_melanogaster.vcf.gz" \
        "$VCF_DATA/drosophila_ensembl.vcf.gz"
}

download_elegans_vcf() {
    echo -e "\n${BOLD}C. elegans VCF: CaeNDR (WBcel235)${NC}"
    echo -e "  ${YELLOW}Note: CaeNDR uses WormBase coordinates. Using Ensembl variation instead:${NC}"
    download "${ENSEMBL_FTP}/variation/vcf/caenorhabditis_elegans/caenorhabditis_elegans.vcf.gz" \
        "$VCF_DATA/elegans_ensembl.vcf.gz"
}

download_arabidopsis_vcf() {
    echo -e "\n${BOLD}Arabidopsis VCF: 1001 Genomes Project (TAIR10)${NC}"
    local PLANTS_FTP="https://ftp.ensemblgenomes.org/pub/plants/release-${ENSEMBL_RELEASE}"
    download "${PLANTS_FTP}/variation/vcf/arabidopsis_thaliana/arabidopsis_thaliana.vcf.gz" \
        "$VCF_DATA/arabidopsis_ensembl.vcf.gz"
}

# ═══════════════════════════════════════════════════════════════
# VCF subsetting: extract N variants from a VCF for benchmarking
# ═══════════════════════════════════════════════════════════════

subset_vcf() {
    local input="$1" output="$2" n="$3"
    if [[ -f "$output" ]]; then
        echo -e "  ${GREEN}Already exists:${NC} $(basename "$output")"
        return 0
    fi
    echo -e "  Extracting $n variants -> $(basename "$output")"
    if [[ "$input" == *.gz ]]; then
        {
            zcat "$input" | grep '^#'
            zcat "$input" | grep -v '^#' | head -n "$n"
        } > "$output"
    else
        {
            grep '^#' "$input"
            grep -v '^#' "$input" | head -n "$n"
        } > "$output"
    fi
}

prepare_benchmark_vcfs() {
    echo -e "\n${BOLD}Preparing benchmark VCF subsets...${NC}"

    # Human: extract chr22 and full subsets from GIAB
    if [[ -f "$VCF_DATA/human_giab_HG002.vcf.gz" ]]; then
        subset_vcf "$VCF_DATA/human_giab_HG002.vcf.gz" "$VCF_DATA/human_100k.vcf" 100000
        subset_vcf "$VCF_DATA/human_giab_HG002.vcf.gz" "$VCF_DATA/human_500k.vcf" 500000
        # Full VCF (decompress all)
        if [[ ! -f "$VCF_DATA/human_giab_full.vcf" ]]; then
            echo -e "  Decompressing full GIAB VCF..."
            zcat "$VCF_DATA/human_giab_HG002.vcf.gz" > "$VCF_DATA/human_giab_full.vcf"
        fi
    fi

    # Mouse
    if [[ -f "$VCF_DATA/mouse_mgp.vcf.gz" ]]; then
        subset_vcf "$VCF_DATA/mouse_mgp.vcf.gz" "$VCF_DATA/mouse_100k.vcf" 100000
        subset_vcf "$VCF_DATA/mouse_mgp.vcf.gz" "$VCF_DATA/mouse_500k.vcf" 500000
    fi

    # Yeast
    if [[ -f "$VCF_DATA/yeast_1002g.gvcf.gz" ]]; then
        subset_vcf "$VCF_DATA/yeast_1002g.gvcf.gz" "$VCF_DATA/yeast_100k.vcf" 100000
    fi

    # Drosophila
    if [[ -f "$VCF_DATA/drosophila_ensembl.vcf.gz" ]]; then
        subset_vcf "$VCF_DATA/drosophila_ensembl.vcf.gz" "$VCF_DATA/drosophila_100k.vcf" 100000
    fi

    # C. elegans
    if [[ -f "$VCF_DATA/elegans_ensembl.vcf.gz" ]]; then
        subset_vcf "$VCF_DATA/elegans_ensembl.vcf.gz" "$VCF_DATA/elegans_100k.vcf" 100000
    fi

    # Arabidopsis
    if [[ -f "$VCF_DATA/arabidopsis_ensembl.vcf.gz" ]]; then
        subset_vcf "$VCF_DATA/arabidopsis_ensembl.vcf.gz" "$VCF_DATA/arabidopsis_100k.vcf" 100000
        subset_vcf "$VCF_DATA/arabidopsis_ensembl.vcf.gz" "$VCF_DATA/arabidopsis_500k.vcf" 500000
    fi
}

# ═══════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════

usage() {
    echo "Usage: $0 [--all | --refs | --vcfs | --human | --mouse | --yeast | --drosophila | --elegans | --arabidopsis]"
    echo ""
    echo "  --all          Download everything (references + VCFs + subsets)"
    echo "  --refs         Download reference genomes only (GFF3 + FASTA)"
    echo "  --vcfs         Download gold-standard VCFs only"
    echo "  --human        Human only (GIAB HG002)"
    echo "  --mouse        Mouse only (MGP)"
    echo "  --yeast        Yeast only (1002 Genomes)"
    echo "  --drosophila   Drosophila only (Ensembl variation)"
    echo "  --elegans      C. elegans only (Ensembl variation)"
    echo "  --arabidopsis  Arabidopsis only (Ensembl variation)"
    echo ""
    echo "Data is saved to: $TEST_DATA/"
}

if [[ $# -eq 0 ]]; then
    usage
    exit 0
fi

for arg in "$@"; do
    case "$arg" in
        --all)
            download_human_ref; download_mouse_ref; download_yeast_ref
            download_drosophila_ref; download_elegans_ref; download_arabidopsis_ref
            download_human_vcf; download_mouse_vcf; download_yeast_vcf
            download_drosophila_vcf; download_elegans_vcf; download_arabidopsis_vcf
            prepare_benchmark_vcfs
            ;;
        --refs)
            download_human_ref; download_mouse_ref; download_yeast_ref
            download_drosophila_ref; download_elegans_ref; download_arabidopsis_ref
            ;;
        --vcfs)
            download_human_vcf; download_mouse_vcf; download_yeast_vcf
            download_drosophila_vcf; download_elegans_vcf; download_arabidopsis_vcf
            prepare_benchmark_vcfs
            ;;
        --human)      download_human_ref; download_human_vcf ;;
        --mouse)      download_mouse_ref; download_mouse_vcf ;;
        --yeast)      download_yeast_ref; download_yeast_vcf ;;
        --drosophila) download_drosophila_ref; download_drosophila_vcf ;;
        --elegans)    download_elegans_ref; download_elegans_vcf ;;
        --arabidopsis) download_arabidopsis_ref; download_arabidopsis_vcf ;;
        --help|-h)    usage; exit 0 ;;
        *)            echo "Unknown option: $arg"; usage; exit 1 ;;
    esac
done

echo -e "\n${GREEN}Done.${NC} Run 'benchmarks/run_benchmark.sh' to benchmark."
