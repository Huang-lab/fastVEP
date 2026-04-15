#!/usr/bin/env bash
# ============================================================================
# fastVEP Deployment: CLINICAL (~20GB data)
# ============================================================================
# Best for: Clinical variant interpretation, research labs
# Server:   Hetzner CCX23 (4 vCPU, 16GB RAM, 80GB SSD) — $30/month
# Genomes:  Human GRCh38 + 2-3 model organisms
# SA:       ClinVar + dbSNP (rs IDs for every known variant)
# ============================================================================
set -euo pipefail

FASTVEP_BIN="${FASTVEP_BIN:-/root/fastVEP/target/release}"
DATA_DIR="/opt/fastvep/data"
BIN_DIR="/opt/fastvep/bin"

echo "=== fastVEP Clinical Deployment ==="
echo "Estimated disk usage: ~20GB"
echo ""

# --- 1. Setup directories and copy binaries ---
echo "[1/6] Setting up directories and binaries..."
mkdir -p "$BIN_DIR" "$DATA_DIR"
cp "$FASTVEP_BIN/fastvep-web" "$BIN_DIR/fastvep-web"
cp "$FASTVEP_BIN/fastvep"     "$BIN_DIR/fastvep"
chmod +x "$BIN_DIR/fastvep-web" "$BIN_DIR/fastvep"

# --- 2. Install samtools if missing ---
if ! command -v samtools &>/dev/null; then
    echo "[1.5/6] Installing samtools..."
    apt-get update -qq && apt-get install -y -qq samtools
fi

# --- 3. Download Human GRCh38 (Ensembl 115) ---
echo "[2/6] Downloading Human GRCh38 genome (~4GB)..."
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

# --- 4. Build SA databases (ClinVar + dbSNP) ---
echo "[3/6] Building ClinVar database (~40MB)..."
cd "$DATA_DIR/hg38_ensembl_115/sa"
if [ ! -f clinvar.osa ]; then
    wget -q --show-progress https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
    "$BIN_DIR/fastvep" sa-build --source clinvar -i clinvar.vcf.gz -o clinvar --assembly GRCh38
    rm -f clinvar.vcf.gz
fi

echo "[4/6] Building dbSNP database (~5GB .osa, ~20GB download — this takes a while)..."
if [ ! -f dbsnp.osa ]; then
    wget -q --show-progress https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
    "$BIN_DIR/fastvep" sa-build --source dbsnp -i GCF_000001405.40.gz -o dbsnp --assembly GRCh38
    rm -f GCF_000001405.40.gz
fi

# --- 5. Download model organisms ---
echo "[5/6] Downloading model organisms..."

# Mouse (~3GB)
echo "  → Mouse (GRCm39)..."
mkdir -p "$DATA_DIR/mouse_grcm39_ensembl_115"
cd "$DATA_DIR/mouse_grcm39_ensembl_115"
if [ ! -f Mus_musculus.GRCm39.115.gff3 ]; then
    wget -q --show-progress https://ftp.ensembl.org/pub/release-115/gff3/mus_musculus/Mus_musculus.GRCm39.115.gff3.gz
    gunzip Mus_musculus.GRCm39.115.gff3.gz
fi
if [ ! -f Mus_musculus.GRCm39.dna.primary_assembly.fa ]; then
    wget -q --show-progress https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
    gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
    samtools faidx Mus_musculus.GRCm39.dna.primary_assembly.fa
fi

# Zebrafish (~2GB)
echo "  → Zebrafish (GRCz11)..."
mkdir -p "$DATA_DIR/zebrafish_grcz11_ensembl_115"
cd "$DATA_DIR/zebrafish_grcz11_ensembl_115"
if [ ! -f Danio_rerio.GRCz11.115.gff3 ]; then
    wget -q --show-progress https://ftp.ensembl.org/pub/release-115/gff3/danio_rerio/Danio_rerio.GRCz11.115.gff3.gz
    gunzip Danio_rerio.GRCz11.115.gff3.gz
fi
if [ ! -f Danio_rerio.GRCz11.dna.primary_assembly.fa ]; then
    wget -q --show-progress https://ftp.ensembl.org/pub/release-115/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
    gunzip Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
    samtools faidx Danio_rerio.GRCz11.dna.primary_assembly.fa
fi

# Drosophila (~200MB)
echo "  → Drosophila (BDGP6)..."
mkdir -p "$DATA_DIR/drosophila_bdgp6_ensembl_115"
cd "$DATA_DIR/drosophila_bdgp6_ensembl_115"
if [ ! -f Drosophila_melanogaster.BDGP6.46.115.gff3 ]; then
    wget -q --show-progress https://ftp.ensembl.org/pub/release-115/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.115.gff3.gz
    gunzip Drosophila_melanogaster.BDGP6.46.115.gff3.gz
fi
if [ ! -f Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa ]; then
    wget -q --show-progress https://ftp.ensembl.org/pub/release-115/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz
    gunzip Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz
    samtools faidx Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa
fi

# C. elegans (~200MB)
echo "  → C. elegans (WBcel235)..."
mkdir -p "$DATA_DIR/celegans_wbcel235_ensembl_115"
cd "$DATA_DIR/celegans_wbcel235_ensembl_115"
if [ ! -f Caenorhabditis_elegans.WBcel235.115.gff3 ]; then
    wget -q --show-progress https://ftp.ensembl.org/pub/release-115/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.115.gff3.gz
    gunzip Caenorhabditis_elegans.WBcel235.115.gff3.gz
fi
if [ ! -f Caenorhabditis_elegans.WBcel235.dna.toplevel.fa ]; then
    wget -q --show-progress https://ftp.ensembl.org/pub/release-115/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
    gunzip Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
    samtools faidx Caenorhabditis_elegans.WBcel235.dna.toplevel.fa
fi

# Yeast (~20MB)
echo "  → Yeast (R64)..."
mkdir -p "$DATA_DIR/yeast_r64_ensembl_115"
cd "$DATA_DIR/yeast_r64_ensembl_115"
if [ ! -f Saccharomyces_cerevisiae.R64-1-1.115.gff3 ]; then
    wget -q --show-progress https://ftp.ensembl.org/pub/release-115/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz
    gunzip Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz
fi
if [ ! -f Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa ]; then
    wget -q --show-progress https://ftp.ensembl.org/pub/release-115/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
    gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
    samtools faidx Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
fi

# Rat (~2GB)
echo "  → Rat (mRatBN7.2)..."
mkdir -p "$DATA_DIR/rat_mratbn72_ensembl_115"
cd "$DATA_DIR/rat_mratbn72_ensembl_115"
if [ ! -f Rattus_norvegicus.mRatBN7.2.115.gff3 ]; then
    wget -q --show-progress https://ftp.ensembl.org/pub/release-115/gff3/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.115.gff3.gz
    gunzip Rattus_norvegicus.mRatBN7.2.115.gff3.gz
fi
if [ ! -f Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa ]; then
    wget -q --show-progress https://ftp.ensembl.org/pub/release-115/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz
    gunzip Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz
    samtools faidx Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa
fi

# Arabidopsis (~200MB)
echo "  → Arabidopsis (TAIR10)..."
mkdir -p "$DATA_DIR/arabidopsis_tair10_ensembl_115"
cd "$DATA_DIR/arabidopsis_tair10_ensembl_115"
if [ ! -f Arabidopsis_thaliana.TAIR10.58.gff3 ]; then
    wget -q --show-progress https://ftp.ensemblgenomes.org/pub/plants/release-58/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.58.gff3.gz
    gunzip Arabidopsis_thaliana.TAIR10.58.gff3.gz
fi
if [ ! -f Arabidopsis_thaliana.TAIR10.dna.toplevel.fa ]; then
    wget -q --show-progress https://ftp.ensemblgenomes.org/pub/plants/release-58/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
    gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
    samtools faidx Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
fi

# --- 6. Set ownership ---
echo "[6/6] Setting permissions..."
chown -R fastvep:fastvep /opt/fastvep 2>/dev/null || true

echo ""
echo "=== Clinical deployment complete ==="
echo ""
echo "Disk usage:"
du -sh "$DATA_DIR"/*/
echo ""
echo "What you get:"
echo "  ✓ Full variant consequence prediction (508k transcripts)"
echo "  ✓ HGVS notation"
echo "  ✓ ClinVar clinical significance"
echo "  ✓ dbSNP rs IDs for all known variants"
echo "  ✓ Mouse, Zebrafish, Drosophila genomes"
echo "  ✗ No population frequencies (gnomAD) — add with deploy-full.sh"
echo ""
echo "Next: run 'systemctl start fastvep' or see DEPLOYMENT.md Section 4"
