#!/usr/bin/env bash
#
# Download the source files needed to build the ACMG concordance benchmark's
# supplementary annotation (SA) databases — ClinVar, gnomAD, REVEL, etc.
#
# After this script finishes, run the per-source build helpers in this
# directory:
#   - fastvep sa-build --source clinvar ...          (clinvar.osa)
#   - fastvep sa-build --source clinvar_protein ...  (clinvar_protein.oga)
#   - fastvep sa-build --source gnomad_gene ...      (gnomad_genes.oga)
#   - fastvep sa-build --source revel ...            (per-chrom revel_chrN.osa)
#   - build_gnomad_per_chrom.sh                      (per-chrom gnomad_chrN.osa)
#   - build_spliceai_phylop.sh                       (per-chrom spliceai/phylop)
#   - clingen_gdv_to_oga.py + fastvep sa-build --source omim ...  (omim.oga)
#
# URL audit: run `bash benchmarks/check_urls.sh acmg` first.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

download() {
    local url="$1" dest="$2"
    if [[ -f "$dest" ]]; then
        echo -e "  ${GREEN}Already exists:${NC} $(basename "$dest")"
        return 0
    fi
    echo -e "  ${YELLOW}Downloading:${NC} $(basename "$dest")"
    curl -L --progress-bar --fail -o "$dest" "$url"
}

cd "$SCRIPT_DIR"

echo "== ClinVar =="
# ClinVar variant VCF on GRCh38 — drives clinvar.osa (per-allele clinical significance)
download "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"          "clinvar.vcf.gz"
# ClinVar variant_summary.txt.gz — drives clinvar_protein.oga (per-protein-position pathogenic catalog for PS1/PM1/PM5)
download "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz" "variant_summary.txt.gz"

echo ""
echo "== gnomAD v4.1 =="
# Gene-level constraint metrics (pLI, LOEUF, mis-Z) — drives gnomad_genes.oga
download "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv" \
    "gnomad.v4.1.constraint_metrics.tsv"
# Allele-level per-chrom exomes VCFs are NOT downloaded in bulk here — the
# benchmark only needs the ClinVar 2-star+ regions. See build_gnomad_per_chrom.sh
# which tabix-extracts only the merged ClinVar 2-star+ regions (clinvar_2star_regions.bed)
# from gs://gcp-public-data--gnomad/release/4.1/vcf/exomes/.

echo ""
echo "== REVEL v1.3 =="
download "https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip" \
    "revel-v1.3_all_chromosomes.zip"

echo ""
echo "== ClinGen Gene-Disease Validity =="
# Used as PVS1 disease-gene fallback (preferred over OMIM per Abou Tayoun 2018).
# Note: returns a redirect to a presigned S3 URL; the CSV is hand-saved as
# clingen_gene_validity.csv and then converted via clingen_gdv_to_oga.py.
download "https://search.clinicalgenome.org/kb/gene-validity/download" \
    "clingen_gene_validity.csv"

echo ""
echo -e "${GREEN}Done.${NC} Next steps:"
echo "  1) Build clinvar/clinvar_protein/revel/gnomad_genes .osa/.oga via:"
echo "       fastvep sa-build --source <name> -i <path> -o <prefix>"
echo "  2) Per-chrom gnomAD exomes:  bash $SCRIPT_DIR/build_gnomad_per_chrom.sh"
echo "  3) SpliceAI + PhyloP (distilled from gnomAD):"
echo "       bash $SCRIPT_DIR/build_spliceai_phylop.sh"
echo "  4) ClinGen GDV -> omim.oga:"
echo "       python3 $SCRIPT_DIR/clingen_gdv_to_oga.py < clingen_gene_validity.csv > clingen_gdv.tsv"
echo "       fastvep sa-build --source omim -i clingen_gdv.tsv -o ../sa_db/omim"
