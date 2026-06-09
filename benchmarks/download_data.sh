#!/usr/bin/env bash
#
# Download reference genomes + benchmark VCFs for the fastVEP speed benchmark
# (multi-organism throughput; see benchmarks/run_benchmark.sh).
#
# Layout produced (matches what run_benchmark.sh expects):
#
#   test_data/organisms/yeast/
#       yeast.gff3            (Ensembl 115, R64-1-1)
#       yeast.fa              (Ensembl 115, R64-1-1 toplevel)
#       yeast.fa.fai
#       yeast_ensembl_full.vcf      (Ensembl 115 variation VCF, full genome)
#   test_data/organisms/drosophila/
#       drosophila.gff3       (Ensembl 115, BDGP6.54)
#       drosophila.fa         (Ensembl 115, BDGP6.54 toplevel)
#       drosophila.fa.fai
#       drosophila_dgrp2_full.vcf   (DGRP2 freeze-2 dm6-lifted, Aerts Lab Zenodo deposit)
#   test_data/organisms/arabidopsis/
#       arabidopsis.gff3      (Ensembl Plants 60, TAIR10)
#       arabidopsis.fa        (Ensembl Plants 60, TAIR10 toplevel)
#       arabidopsis.fa.fai
#       arabidopsis_1001g_full.vcf  (1001 Genomes — see notes below)
#   test_data/organisms/mouse/
#       mouse.gff3            (Ensembl 115, GRCm39)
#       mouse.fa              (Ensembl 115, GRCm39 primary assembly)
#       mouse.fa.fai
#       mouse_ensembl_1m.vcf  (MGP REL-2112 snps, head -n 1,000,000)
#   test_data/organisms/elegans/  (optional; opt in with `elegans`)
#       elegans.gff3          (Ensembl 115, WBcel235)
#       elegans.fa            (Ensembl 115, WBcel235 toplevel)
#       elegans.fa.fai
#       elegans_caendr_full.vcf  (CaeNDR WI.20250625 hard-filter isotype)
#   test_data/organisms/human/
#       Homo_sapiens.GRCh38.115.gff3                    (Ensembl 115)
#       Homo_sapiens.GRCh38.dna.primary_assembly.fa     (Ensembl 115)
#       Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
#       human_giab_hg002_full.vcf                       (GIAB HG002 v4.2.1)
#       HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed   (high-confidence regions)
#
# All URLs are verified working at the time of this script. Run
#   `bash benchmarks/check_urls.sh`
# to re-verify before a fresh download (see benchmarks/check_urls.sh).
#
# Notes on data substitutions:
#
#   * Mouse MGP CAST/EiJ-specific VCF (~26M variants): the per-strain
#     v8 file at ftp.ebi.ac.uk/.../CAST_EiJ.mgp.v8.snps.dbSNP142.vcf.gz
#     is 404 as of 2026, and ftp.mousegenomes.org/pub/release/current_snps/
#     is unreachable. We instead download the aggregated REL-2112
#     mgp_REL2021_snps.vcf.gz on the same EBI mirror (~80M variants;
#     subset to 1M for the benchmark slot).
#
#   * 1001 Genomes Arabidopsis (12.9M variants): the original
#     1001genomes.org SNPs_all_methods/*.vcf.gz path is 404, and the
#     site is now a SPA. The intersection_snp_short_indel_vcf/ dir
#     contains per-strain VCFs only, not a combined VCF. We use the
#     Ensembl Plants 60 variation VCF as a working equivalent on the
#     same TAIR10 coordinates.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
TEST_DATA="$PROJECT_DIR/test_data"
ORG_DATA="$TEST_DATA/organisms"

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BOLD='\033[1m'
NC='\033[0m'

mkdir -p "$ORG_DATA"

# ═══════════════════════════════════════════════════════════════
# Source URL constants (kept here so they're easy to audit / update)
# ═══════════════════════════════════════════════════════════════
ENSEMBL_RELEASE=115
ENSEMBL_FTP="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}"
# Ensembl Plants uses its own release numbering (NOT main Ensembl's).
# Verified working at release-60 (matches Ensembl main 115 era data freeze).
PLANTS_RELEASE=60
PLANTS_FTP="http://ftp.ensemblgenomes.org/pub/plants/release-${PLANTS_RELEASE}"

# ═══════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════
download() {
    local url="$1" dest="$2"
    if [[ -f "$dest" ]]; then
        echo -e "  ${GREEN}Already exists:${NC} $(basename "$dest")"
        return 0
    fi
    mkdir -p "$(dirname "$dest")"
    echo -e "  ${YELLOW}Downloading:${NC} $(basename "$dest")"
    curl -L --progress-bar --fail -o "$dest" "$url"
}

decompress() {
    local gz="$1" out="${1%.gz}"
    if [[ -f "$out" ]]; then return 0; fi
    echo -e "  Decompressing $(basename "$gz")..."
    gunzip -k "$gz"
}

# Stream-decompress a .gz to a different output path. Uses `gunzip -c` because
# macOS `zcat` only handles `.Z` files, not `.gz`.
gunzip_to() {
    local gz="$1" out="$2"
    [[ -f "$out" ]] && return 0
    echo -e "  Decompressing $(basename "$gz") -> $(basename "$out")..."
    gunzip -c "$gz" > "$out"
}

index_fasta() {
    local fa="$1"
    [[ -f "${fa}.fai" ]] && return 0
    if command -v samtools &>/dev/null; then
        echo -e "  Indexing $(basename "$fa")..."
        samtools faidx "$fa"
    else
        echo -e "  ${YELLOW}samtools not found — skipping FASTA index (install for mmap support)${NC}"
    fi
}

# Subset a VCF.gz to the first N variant rows, plus the full header.
subset_vcf() {
    local input="$1" output="$2" n="$3"
    [[ -f "$output" ]] && return 0
    echo -e "  Extracting $n variants -> $(basename "$output")"
    if [[ "$input" == *.gz || "$input" == *.bgz ]]; then
        { gunzip -c "$input" | grep '^#'; gunzip -c "$input" | grep -v '^#' | head -n "$n"; } > "$output"
    else
        { grep '^#' "$input"; grep -v '^#' "$input" | head -n "$n"; } > "$output"
    fi
}

# ═══════════════════════════════════════════════════════════════
# Yeast — Ensembl 115, S. cerevisiae R64-1-1
# ═══════════════════════════════════════════════════════════════
download_yeast() {
    echo -e "\n${BOLD}Yeast (R64-1-1, Ensembl ${ENSEMBL_RELEASE})${NC}"
    local D="$ORG_DATA/yeast"

    download "${ENSEMBL_FTP}/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.${ENSEMBL_RELEASE}.gff3.gz" "$D/yeast.gff3.gz"
    decompress "$D/yeast.gff3.gz"

    download "${ENSEMBL_FTP}/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz" "$D/yeast.fa.gz"
    decompress "$D/yeast.fa.gz"
    index_fasta "$D/yeast.fa"

    # Variation: Ensembl 115 saccharomyces_cerevisiae variation VCF (full genome, ~260K variants)
    download "${ENSEMBL_FTP}/variation/vcf/saccharomyces_cerevisiae/saccharomyces_cerevisiae.vcf.gz" "$D/yeast_ensembl.vcf.gz"
    gunzip_to "$D/yeast_ensembl.vcf.gz" "$D/yeast_ensembl_full.vcf"
}

# ═══════════════════════════════════════════════════════════════
# Drosophila — Ensembl 115, BDGP6.54
# ═══════════════════════════════════════════════════════════════
download_drosophila() {
    echo -e "\n${BOLD}Drosophila (BDGP6.54, Ensembl ${ENSEMBL_RELEASE})${NC}"
    local D="$ORG_DATA/drosophila"

    download "${ENSEMBL_FTP}/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.54.${ENSEMBL_RELEASE}.gff3.gz" "$D/drosophila.gff3.gz"
    decompress "$D/drosophila.gff3.gz"

    download "${ENSEMBL_FTP}/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa.gz" "$D/drosophila.fa.gz"
    decompress "$D/drosophila.fa.gz"
    index_fasta "$D/drosophila.fa"

    # DGRP2 freeze-2 (~4.4M variants), dm6-lifted with dbSNP IDs.
    # NCSU's original http://dgrp2.gnets.ncsu.edu/ is offline; this is the
    # Aerts Lab Zenodo deposit (Zenodo DOI 10.5281/zenodo.155396).
    # dm6 == BDGP6 coordinates, so this matches the Ensembl GFF3 / FASTA above.
    download "https://zenodo.org/record/155396/files/dgrp2_dm6_dbSNP.vcf.gz" "$D/drosophila_dgrp2.vcf.gz"
    gunzip_to "$D/drosophila_dgrp2.vcf.gz" "$D/drosophila_dgrp2_full.vcf"
}

# ═══════════════════════════════════════════════════════════════
# Arabidopsis — Ensembl Plants 60, TAIR10
# ═══════════════════════════════════════════════════════════════
download_arabidopsis() {
    echo -e "\n${BOLD}Arabidopsis (TAIR10, Ensembl Plants ${PLANTS_RELEASE})${NC}"
    local D="$ORG_DATA/arabidopsis"

    download "${PLANTS_FTP}/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.${PLANTS_RELEASE}.gff3.gz" "$D/arabidopsis.gff3.gz"
    decompress "$D/arabidopsis.gff3.gz"

    download "${PLANTS_FTP}/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz" "$D/arabidopsis.fa.gz"
    decompress "$D/arabidopsis.fa.gz"
    index_fasta "$D/arabidopsis.fa"

    # 1001 Genomes direct download no longer works; use Ensembl Plants
    # variation VCF as a working equivalent on the same TAIR10 coords.
    download "${PLANTS_FTP}/variation/vcf/arabidopsis_thaliana/arabidopsis_thaliana.vcf.gz" "$D/arabidopsis_ensembl.vcf.gz"
    gunzip_to "$D/arabidopsis_ensembl.vcf.gz" "$D/arabidopsis_1001g_full.vcf"
}

# ═══════════════════════════════════════════════════════════════
# Mouse — Ensembl 115, GRCm39
# ═══════════════════════════════════════════════════════════════
download_mouse() {
    echo -e "\n${BOLD}Mouse (GRCm39, Ensembl ${ENSEMBL_RELEASE})${NC}"
    local D="$ORG_DATA/mouse"

    download "${ENSEMBL_FTP}/gff3/mus_musculus/Mus_musculus.GRCm39.${ENSEMBL_RELEASE}.gff3.gz" "$D/mouse.gff3.gz"
    decompress "$D/mouse.gff3.gz"

    download "${ENSEMBL_FTP}/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz" "$D/mouse.fa.gz"
    decompress "$D/mouse.fa.gz"
    index_fasta "$D/mouse.fa"

    # Mouse Genomes Project REL-2112 aggregated SNPs across all strains.
    # The strain-specific CAST_EiJ v8 file (used in the manuscript)
    # is no longer hosted; aggregated REL2021 is the closest working
    # equivalent.
    download "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz" "$D/mouse_mgp.vcf.gz"
    # Subset to ~1M variants for the "mouse_1m" benchmark slot.
    subset_vcf "$D/mouse_mgp.vcf.gz" "$D/mouse_ensembl_1m.vcf" 1000000
}

# ═══════════════════════════════════════════════════════════════
# C. elegans — Ensembl 115, WBcel235 + CaeNDR WI.20250625
# (optional; not part of --all, opt in with `elegans`)
# ═══════════════════════════════════════════════════════════════
download_elegans() {
    echo -e "\n${BOLD}C. elegans (WBcel235, Ensembl ${ENSEMBL_RELEASE})${NC}"
    local D="$ORG_DATA/elegans"

    download "${ENSEMBL_FTP}/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.${ENSEMBL_RELEASE}.gff3.gz" "$D/elegans.gff3.gz"
    decompress "$D/elegans.gff3.gz"

    download "${ENSEMBL_FTP}/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz" "$D/elegans.fa.gz"
    decompress "$D/elegans.fa.gz"
    index_fasta "$D/elegans.fa"

    # CaeNDR — hard-filter isotype VCF, latest release as of 2026.
    # ~11 GB compressed; if you only want a smoke test set the env var
    # CAENDR_SUBSET=1 (any non-empty value) to subset to 1M variants instead.
    download "https://caendr-open-access-data-bucket.s3.us-east-2.amazonaws.com/dataset_release/c_elegans/20250625/variation/WI.20250625.hard-filter.isotype.vcf.gz" \
        "$D/elegans_caendr.vcf.gz"
    if [[ -n "${CAENDR_SUBSET:-}" ]]; then
        subset_vcf "$D/elegans_caendr.vcf.gz" "$D/elegans_caendr_full.vcf" 1000000
    else
        gunzip_to "$D/elegans_caendr.vcf.gz" "$D/elegans_caendr_full.vcf"
    fi
}

# ═══════════════════════════════════════════════════════════════
# Human — Ensembl 115, GRCh38 + GIAB HG002 v4.2.1
# ═══════════════════════════════════════════════════════════════
download_human() {
    echo -e "\n${BOLD}Human (GRCh38, Ensembl ${ENSEMBL_RELEASE})${NC}"
    local D="$ORG_DATA/human"

    download "${ENSEMBL_FTP}/gff3/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gff3.gz" "$D/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gff3.gz"
    decompress "$D/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gff3.gz"

    download "${ENSEMBL_FTP}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" "$D/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    decompress "$D/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    index_fasta "$D/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

    # GIAB HG002 v4.2.1 benchmark VCF (~4.05M variants) + high-confidence BED
    local GIAB_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38"
    download "${GIAB_BASE}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"     "$D/human_giab_hg002.vcf.gz"
    download "${GIAB_BASE}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" "$D/human_giab_hg002.vcf.gz.tbi"
    # High-confidence regions BED (use with `bcftools view -R` to restrict
    # benchmark scoring to truth-set callable regions, per GIAB convention).
    download "${GIAB_BASE}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
        "$D/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
    gunzip_to "$D/human_giab_hg002.vcf.gz" "$D/human_giab_hg002_full.vcf"
}

# ═══════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════
usage() {
    cat <<EOF
Usage: $0 [--all | <organism>...]

Organisms: yeast, drosophila, arabidopsis, mouse, human, elegans

  --all      Download the 5 manuscript organisms (~12 GB raw, ~30 GB decompressed):
             yeast, drosophila, arabidopsis, mouse, human.
             Does NOT include elegans by default — the CaeNDR VCF is 11 GB.
  elegans    Add C. elegans (CaeNDR WI.20250625 hard-filter isotype). Optional;
             set CAENDR_SUBSET=1 to subset to 1M variants instead of the full ~11 GB.

Data is placed under: $ORG_DATA/<organism>/
Then run: bash benchmarks/run_benchmark.sh
EOF
}

[[ $# -eq 0 ]] && { usage; exit 0; }

for arg in "$@"; do
    case "$arg" in
        --all)
            download_yeast; download_drosophila; download_arabidopsis
            download_mouse; download_human
            ;;
        yeast)        download_yeast ;;
        drosophila)   download_drosophila ;;
        arabidopsis)  download_arabidopsis ;;
        mouse)        download_mouse ;;
        elegans)      download_elegans ;;
        human)        download_human ;;
        --help|-h)    usage; exit 0 ;;
        *)            echo "Unknown option: $arg"; usage; exit 1 ;;
    esac
done

echo -e "\n${GREEN}Done.${NC} Run 'bash benchmarks/run_benchmark.sh' to benchmark."
