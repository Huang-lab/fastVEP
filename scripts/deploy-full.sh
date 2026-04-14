#!/usr/bin/env bash
# ============================================================================
# fastVEP Deployment: FULL RESEARCH (~55GB data for hg38 only)
# ============================================================================
# Best for: Production clinical genomics, research institutions
# Server:   Hetzner CCX33 (8 vCPU, 32GB RAM, 160GB SSD) — $55/month
#           OR CCX23 + 200GB Hetzner Volume ($30 + $10/month)
# Genomes:  Human GRCh38 + GRCh37 + model organisms
# SA:       ClinVar + dbSNP + gnomAD (downloaded one chr at a time to save disk)
#
# IMPORTANT: gnomAD raw VCFs are 20-40GB *per chromosome*. This script
# downloads one at a time, builds the index, then deletes the VCF.
# You need ~40GB of TEMPORARY free space during the gnomAD build.
# ============================================================================
set -euo pipefail

FASTVEP_BIN="${FASTVEP_BIN:-/root/fastVEP/target/release}"
DATA_DIR="/opt/fastvep/data"
BIN_DIR="/opt/fastvep/bin"

echo "=== fastVEP Full Research Deployment ==="
echo ""

# --- Check disk space ---
AVAIL_GB=$(df --output=avail "$DATA_DIR" 2>/dev/null | tail -1 | awk '{printf "%.0f", $1/1024/1024}')
echo "Available disk space: ${AVAIL_GB}GB"
if [ "$AVAIL_GB" -lt 60 ]; then
    echo "WARNING: You have less than 60GB free. gnomAD alone needs ~40GB temp space."
    echo "Consider upgrading to CCX33 (160GB) or attaching a Hetzner Volume."
    read -p "Continue anyway? [y/N] " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# --- 1. Setup directories and copy binaries ---
echo "[1/8] Setting up directories and binaries..."
mkdir -p "$BIN_DIR" "$DATA_DIR"
cp "$FASTVEP_BIN/fastvep-web" "$BIN_DIR/fastvep-web"
cp "$FASTVEP_BIN/fastvep"     "$BIN_DIR/fastvep"
chmod +x "$BIN_DIR/fastvep-web" "$BIN_DIR/fastvep"

# --- 2. Install samtools if missing ---
if ! command -v samtools &>/dev/null; then
    echo "[1.5/8] Installing samtools..."
    apt-get update -qq && apt-get install -y -qq samtools
fi

# --- 3. Download Human GRCh38 ---
echo "[2/8] Downloading Human GRCh38 genome (~4GB)..."
mkdir -p "$DATA_DIR/hg38_ensembl_115/sa"
cd "$DATA_DIR/hg38_ensembl_115"

if [ ! -f Homo_sapiens.GRCh38.115.gff3 ]; then
    wget -q --show-progress https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz
    gunzip Homo_sapiens.GRCh38.115.gff3.gz
fi

if [ ! -f Homo_sapiens.GRCh38.dna.primary_assembly.fa ]; then
    wget -q --show-progress https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
fi

# --- 4. Build SA databases for GRCh38 ---
cd "$DATA_DIR/hg38_ensembl_115/sa"

echo "[3/8] Building ClinVar (GRCh38)..."
if [ ! -f clinvar.osa ]; then
    wget -q --show-progress https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
    "$BIN_DIR/fastvep" sa-build --source clinvar -i clinvar.vcf.gz -o clinvar --assembly GRCh38
    rm -f clinvar.vcf.gz
fi

echo "[4/8] Building dbSNP (GRCh38, ~20GB download → ~5GB .osa)..."
if [ ! -f dbsnp.osa ]; then
    wget -q --show-progress https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
    "$BIN_DIR/fastvep" sa-build --source dbsnp -i GCF_000001405.40.gz -o dbsnp --assembly GRCh38
    # Clean up the 20GB download immediately
    rm -f GCF_000001405.40.gz
fi

echo "[5/8] Building gnomAD v4 (GRCh38) — downloading one chromosome at a time..."
echo "       This will take several hours. Each chr is 7-25GB download."
if [ ! -f gnomad.osa ]; then
    for chr in {1..22} X; do
        VCF="gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz"
        echo "  → Downloading chr${chr}..."
        wget -q --show-progress \
            "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/${VCF}"
        echo "  → Indexing chr${chr}..."
        "$BIN_DIR/fastvep" sa-build --source gnomad -i "$VCF" -o gnomad --assembly GRCh38
        echo "  → Cleaning up chr${chr} VCF..."
        rm -f "$VCF"
        echo "  → chr${chr} done. Disk free: $(df -h --output=avail "$DATA_DIR" | tail -1 | xargs)"
    done
fi

# --- 5. Download Human GRCh37 (hg19) ---
echo "[6/8] Downloading Human GRCh37 genome (~4GB)..."
mkdir -p "$DATA_DIR/hg19_ensembl_113/sa"
cd "$DATA_DIR/hg19_ensembl_113"

if [ ! -f Homo_sapiens.GRCh37.87.gff3 ]; then
    wget -q --show-progress https://ftp.ensembl.org/pub/grch37/release-113/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz
    gunzip Homo_sapiens.GRCh37.87.gff3.gz
fi

if [ ! -f Homo_sapiens.GRCh37.dna.primary_assembly.fa ]; then
    wget -q --show-progress https://ftp.ensembl.org/pub/grch37/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
    gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
    samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa
fi

# ClinVar for GRCh37
cd "$DATA_DIR/hg19_ensembl_113/sa"
echo "  → Building ClinVar (GRCh37)..."
if [ ! -f clinvar.osa ]; then
    wget -q --show-progress https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
    "$BIN_DIR/fastvep" sa-build --source clinvar -i clinvar.vcf.gz -o clinvar --assembly GRCh37
    rm -f clinvar.vcf.gz
fi

# --- 6. Download model organisms ---
echo "[7/8] Downloading model organisms..."

download_organism() {
    local name="$1" dir="$2" gff_url="$3" fa_url="$4"
    echo "  → ${name}..."
    mkdir -p "$DATA_DIR/$dir"
    cd "$DATA_DIR/$dir"
    local gff_file=$(basename "$gff_url" .gz)
    local fa_file=$(basename "$fa_url" .gz)
    if [ ! -f "$gff_file" ]; then
        wget -q --show-progress "$gff_url"
        gunzip "$(basename "$gff_url")"
    fi
    if [ ! -f "$fa_file" ]; then
        wget -q --show-progress "$fa_url"
        gunzip "$(basename "$fa_url")"
        samtools faidx "$fa_file"
    fi
}

download_organism "Mouse (GRCm39)" "mouse_grcm39_ensembl_115" \
    "https://ftp.ensembl.org/pub/release-115/gff3/mus_musculus/Mus_musculus.GRCm39.115.gff3.gz" \
    "https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"

download_organism "Zebrafish (GRCz11)" "zebrafish_grcz11_ensembl_115" \
    "https://ftp.ensembl.org/pub/release-115/gff3/danio_rerio/Danio_rerio.GRCz11.115.gff3.gz" \
    "https://ftp.ensembl.org/pub/release-115/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz"

download_organism "Drosophila (BDGP6)" "drosophila_bdgp6_ensembl_115" \
    "https://ftp.ensembl.org/pub/release-115/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.115.gff3.gz" \
    "https://ftp.ensembl.org/pub/release-115/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz"

download_organism "C. elegans (WBcel235)" "celegans_wbcel235_ensembl_115" \
    "https://ftp.ensembl.org/pub/release-115/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.115.gff3.gz" \
    "https://ftp.ensembl.org/pub/release-115/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz"

download_organism "Yeast (R64-1-1)" "yeast_r64_ensembl_115" \
    "https://ftp.ensembl.org/pub/release-115/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz" \
    "https://ftp.ensembl.org/pub/release-115/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"

# --- 7. Set ownership ---
echo "[8/8] Setting permissions..."
chown -R fastvep:fastvep /opt/fastvep 2>/dev/null || true

echo ""
echo "=== Full Research deployment complete ==="
echo ""
echo "Disk usage:"
du -sh "$DATA_DIR"/*/
echo ""
echo "Total:"
du -sh "$DATA_DIR"
echo ""
echo "What you get:"
echo "  ✓ Full variant consequence prediction"
echo "  ✓ HGVS notation"
echo "  ✓ ClinVar clinical significance (GRCh38 + GRCh37)"
echo "  ✓ dbSNP rs IDs (GRCh38)"
echo "  ✓ gnomAD v4 population frequencies (GRCh38)"
echo "  ✓ Human GRCh38 + GRCh37"
echo "  ✓ Mouse, Zebrafish, Drosophila, C. elegans, Yeast"
echo ""
echo "Next: run 'systemctl start fastvep' or see DEPLOYMENT.md Section 4"
