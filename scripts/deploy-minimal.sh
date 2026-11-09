#!/usr/bin/env bash
# ============================================================================
# fastVEP Deployment: MINIMAL (~5GB data)
# ============================================================================
# Best for: Demos, testing, budget-conscious deployments
# Server:   Hetzner CCX23 (4 vCPU, 16GB RAM, 160GB SSD) - $30/month
# Genomes:  Human GRCh38 only
# SA:       ClinVar only (free, tiny, most clinically useful)
# ============================================================================
set -euo pipefail

FASTVEP_BIN="${FASTVEP_BIN:-/root/fastVEP/target/release}"
DATA_DIR="/opt/fastvep/data"
BIN_DIR="/opt/fastvep/bin"

echo "=== fastVEP Minimal Deployment ==="
echo "Estimated disk usage: ~5GB"
echo ""

# --- 1. Setup directories and copy binaries ---
echo "[1/4] Setting up directories and binaries..."
mkdir -p "$BIN_DIR" "$DATA_DIR"
cp "$FASTVEP_BIN/fastvep-web" "$BIN_DIR/fastvep-web"
cp "$FASTVEP_BIN/fastvep"     "$BIN_DIR/fastvep"
chmod +x "$BIN_DIR/fastvep-web" "$BIN_DIR/fastvep"

# --- 2. Install samtools if missing ---
if ! command -v samtools &>/dev/null; then
    echo "[1.5/4] Installing samtools..."
    apt-get update -qq && apt-get install -y -qq samtools
fi

# --- 3. Download Human GRCh38 (Ensembl 115) ---
echo "[2/4] Downloading Human GRCh38 genome (~4GB)..."
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

# --- 4. Build ClinVar SA database ---
echo "[3/4] Building ClinVar database (~40MB)..."
cd "$DATA_DIR/hg38_ensembl_115/sa"
# sa-build writes clinvar.osa2 (clinvar.osa in the v1 format); either one means it is built.
if [ ! -s clinvar.osa2 ] && [ ! -s clinvar.osa ]; then
    wget -q --show-progress https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
    "$BIN_DIR/fastvep" sa-build --source clinvar -i clinvar.vcf.gz -o clinvar --assembly GRCh38
    rm -f clinvar.vcf.gz
fi

# --- 5. Set ownership ---
echo "[4/4] Setting permissions..."
chown -R fastvep:fastvep /opt/fastvep 2>/dev/null || true

echo ""
echo "=== Minimal deployment complete ==="
echo "Data:     $DATA_DIR/hg38_ensembl_115/"
echo "Binaries: $BIN_DIR/"
echo ""
echo "Disk usage:"
du -sh "$DATA_DIR/hg38_ensembl_115/"
echo ""
echo "What you get:"
echo "  ✓ Full variant consequence prediction (508k transcripts)"
echo "  ✓ HGVS notation"
echo "  ✓ ClinVar clinical significance"
echo "  ✗ No population frequencies (gnomAD)"
echo "  ✗ No dbSNP rs IDs"
echo ""
echo "Next: run 'systemctl start fastvep' or see DEPLOYMENT.md Section 4"
