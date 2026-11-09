#!/usr/bin/env bash
# ============================================================================
# fastVEP Deployment: FULL RESEARCH (~26GB of data, plus gnomAD)
# ============================================================================
# Best for: Production clinical genomics, research institutions
# Server:   Hetzner CCX23 (4 vCPU, 16GB RAM, 160GB SSD) - $30/month
# Genomes:  Human GRCh38 + GRCh37 + seven other organisms
# SA:       ClinVar (both builds) + dbSNP (GRCh38) + gnomAD v4.1 genomes (GRCh38)
#
# IMPORTANT: gnomAD v4.1 genomes VCFs are 0.6 to 47GB *per chromosome* (chr2 is
# the largest; 563GB in all). This script downloads one at a time, builds its
# database, then deletes the VCF. You need ~47GB of TEMPORARY free space during
# the gnomAD build, on top of the ~26GB of data.
#
# This does not build the rest of the ACMG reference stack (REVEL, SpliceAI,
# PhyloP, gnomAD gene constraints, ClinGen, RepeatMasker). Without it the web
# server refuses ACMG requests for the genome; see docs/ACMG_SETUP.md.
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
if [ "$AVAIL_GB" -lt 80 ]; then
    echo "WARNING: You have less than 80GB free. The data take ~26GB and the gnomAD build needs ~47GB of temporary space."
    echo "Consider a bigger disk or a Hetzner Volume."
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
# sa-build writes .osa2 (.osa in the v1 format); either one means the database is built.
if [ ! -s clinvar.osa2 ] && [ ! -s clinvar.osa ]; then
    wget -q --show-progress https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
    "$BIN_DIR/fastvep" sa-build --source clinvar -i clinvar.vcf.gz -o clinvar --assembly GRCh38
    rm -f clinvar.vcf.gz
fi

echo "[4/8] Building dbSNP (GRCh38, ~30GB download → ~9.4GB .osa2)..."
if [ ! -s dbsnp.osa2 ] && [ ! -s dbsnp.osa ]; then
    wget -q --show-progress https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
    # 16-bit chunks: dbSNP is dense, and at the default width a scattered lookup costs 23 ms (docs/SUPPLEMENTARY_ANNOTATIONS.md)
    "$BIN_DIR/fastvep" sa-build --source dbsnp -i GCF_000001405.40.gz -o dbsnp --assembly GRCh38 --chunk-bits 16
    # Clean up the 30GB download immediately
    rm -f GCF_000001405.40.gz
fi

echo "[5/8] Building gnomAD v4.1 genomes (GRCh38), one chromosome at a time..."
echo "       This will take several hours. Each chr is 0.6-47GB to download."
# sa-build writes the database under the name it is given, so each chromosome gets a name of its own.
# One name for all of them leaves the database of the last chromosome alone: a gnomAD that says
# "absent" for every variant of chromosomes 1 to 22.
# A chromosome is built outside sa/ and moved in when it is whole, so that a build killed halfway
# is not mistaken for a finished one (and not loaded by the server) when this script is run again.
SA_DIR="$DATA_DIR/hg38_ensembl_115/sa"
BUILD_DIR="$DATA_DIR/hg38_ensembl_115/build.tmp"
mkdir -p "$BUILD_DIR"
for chr in {1..22} X; do
    if [ -s "$SA_DIR/gnomad_chr${chr}.osa2" ] || [ -s "$SA_DIR/gnomad_chr${chr}.osa" ]; then
        echo "  → chr${chr}: already built"
        continue
    fi
    VCF="gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz"
    echo "  → Downloading chr${chr}..."
    wget -q --show-progress -O "$BUILD_DIR/$VCF" \
        "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/${VCF}"
    echo "  → Indexing chr${chr}..."
    "$BIN_DIR/fastvep" sa-build --source gnomad -i "$BUILD_DIR/$VCF" -o "$BUILD_DIR/gnomad_chr${chr}" --assembly GRCh38 --chunk-bits 16
    mv "$BUILD_DIR"/gnomad_chr${chr}.osa* "$SA_DIR/"
    echo "  → Cleaning up chr${chr} VCF..."
    rm -f "$BUILD_DIR/$VCF"
    echo "  → chr${chr} done. Disk free: $(df -h --output=avail "$DATA_DIR" | tail -1 | xargs)"
done
rmdir "$BUILD_DIR" 2>/dev/null || true
BUILT=$(ls "$SA_DIR"/gnomad_chr*.osa2 "$SA_DIR"/gnomad_chr*.osa 2>/dev/null | wc -l)
if [ "$BUILT" -ne 23 ]; then
    echo "WARNING: $BUILT gnomAD databases in $SA_DIR, expected 23 (chr1-22 and X)."
fi

# --- 5. Download Human GRCh37 ---
# The gene set of GRCh37 is Ensembl 87: Homo_sapiens.GRCh37.87.gff3.gz is the same file (size and date)
# under release-87, release-113 and current, so the directory is named after the gene set.
echo "[6/8] Downloading Human GRCh37 genome (~4GB)..."
mkdir -p "$DATA_DIR/hg37_ensembl_87/sa"
cd "$DATA_DIR/hg37_ensembl_87"

if [ ! -f Homo_sapiens.GRCh37.87.gff3 ]; then
    wget -q --show-progress https://ftp.ensembl.org/pub/grch37/release-87/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz
    gunzip Homo_sapiens.GRCh37.87.gff3.gz
fi

if [ ! -f Homo_sapiens.GRCh37.dna.primary_assembly.fa ]; then
    wget -q --show-progress https://ftp.ensembl.org/pub/grch37/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
    gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
    samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa
fi

# ClinVar for GRCh37
cd "$DATA_DIR/hg37_ensembl_87/sa"
echo "  → Building ClinVar (GRCh37)..."
if [ ! -s clinvar.osa2 ] && [ ! -s clinvar.osa ]; then
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
    "https://ftp.ensembl.org/pub/release-115/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.54.115.gff3.gz" \
    "https://ftp.ensembl.org/pub/release-115/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa.gz"

download_organism "C. elegans (WBcel235)" "celegans_wbcel235_ensembl_115" \
    "https://ftp.ensembl.org/pub/release-115/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.115.gff3.gz" \
    "https://ftp.ensembl.org/pub/release-115/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz"

download_organism "Yeast (R64-1-1)" "yeast_r64_ensembl_115" \
    "https://ftp.ensembl.org/pub/release-115/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz" \
    "https://ftp.ensembl.org/pub/release-115/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"

download_organism "Rat (GRCr8)" "rat_grcr8_ensembl_115" \
    "https://ftp.ensembl.org/pub/release-115/gff3/rattus_norvegicus/Rattus_norvegicus.GRCr8.115.gff3.gz" \
    "https://ftp.ensembl.org/pub/release-115/fasta/rattus_norvegicus/dna/Rattus_norvegicus.GRCr8.dna.toplevel.fa.gz"

# Ensembl Plants has its own release numbers (63 here, against 115 of Ensembl), and the directory says so.
download_organism "Arabidopsis (TAIR10)" "arabidopsis_tair10_ensemblplants_63" \
    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.63.gff3.gz" \
    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"

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
echo "  ✓ gnomAD v4.1 genomes population frequencies (GRCh38)"
echo "  ✓ Human GRCh38 + GRCh37"
echo "  ✓ Mouse, rat, zebrafish, fruit fly, C. elegans, yeast, Arabidopsis"
echo "  ✗ Not the rest of the ACMG reference stack: see docs/ACMG_SETUP.md"
echo ""
echo "Next: run 'systemctl start fastvep' or see DEPLOYMENT.md Section 4"
