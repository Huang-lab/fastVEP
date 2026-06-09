#!/usr/bin/env bash
#
# Verifies every external download URL used by the fastVEP benchmark suite.
# Run this before a fresh download so you can fail fast on broken upstream
# mirrors (Ensembl release rotations, host renames, etc.).
#
# Usage:
#   bash benchmarks/check_urls.sh           # check all
#   bash benchmarks/check_urls.sh speed     # speed benchmark only
#   bash benchmarks/check_urls.sh acmg      # ACMG concordance only
#
set -uo pipefail

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

ok=0
fail=0
failed_urls=()

check() {
    local label="$1" url="$2"
    local code
    # -L to follow redirects (Zenodo issues 301; S3 occasionally redirects).
    code=$(curl -s -o /dev/null -w "%{http_code}" --max-time 25 -L -I "$url" 2>/dev/null)
    if [[ "$code" != "200" ]]; then
        # Some FTP-over-HTTP gateways disallow HEAD; retry with a tiny range GET
        code=$(curl -s -o /dev/null -w "%{http_code}" --max-time 25 -L -r 0-511 "$url" 2>/dev/null)
    fi
    if [[ "$code" == "200" || "$code" == "206" ]]; then
        echo -e "  ${GREEN}OK${NC}    [$code] $label"
        ok=$((ok+1))
    else
        echo -e "  ${RED}FAIL${NC}  [$code] $label  ${YELLOW}$url${NC}"
        fail=$((fail+1))
        failed_urls+=("$label  $url")
    fi
}

check_speed() {
    echo -e "\n== Speed benchmark (benchmarks/download_data.sh) =="
    check "Yeast GFF3 (Ensembl 115)"        "https://ftp.ensembl.org/pub/release-115/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz"
    check "Yeast FASTA (Ensembl 115)"       "https://ftp.ensembl.org/pub/release-115/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
    check "Yeast variation VCF"             "https://ftp.ensembl.org/pub/release-115/variation/vcf/saccharomyces_cerevisiae/saccharomyces_cerevisiae.vcf.gz"
    check "Drosophila GFF3 (BDGP6.54)"      "https://ftp.ensembl.org/pub/release-115/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.54.115.gff3.gz"
    check "Drosophila FASTA (BDGP6.54)"     "https://ftp.ensembl.org/pub/release-115/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa.gz"
    check "Drosophila DGRP2 dm6 (Zenodo)"   "https://zenodo.org/record/155396/files/dgrp2_dm6_dbSNP.vcf.gz"
    check "Arabidopsis GFF3 (Plants 60)"    "http://ftp.ensemblgenomes.org/pub/plants/release-60/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.60.gff3.gz"
    check "Arabidopsis FASTA (Plants 60)"   "http://ftp.ensemblgenomes.org/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
    check "Arabidopsis variation VCF"       "http://ftp.ensemblgenomes.org/pub/plants/release-60/variation/vcf/arabidopsis_thaliana/arabidopsis_thaliana.vcf.gz"
    check "Mouse GFF3 (GRCm39)"             "https://ftp.ensembl.org/pub/release-115/gff3/mus_musculus/Mus_musculus.GRCm39.115.gff3.gz"
    check "Mouse FASTA (GRCm39)"            "https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
    check "Mouse MGP REL-2112 snps"         "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz"
    check "C. elegans GFF3 (WBcel235)"      "https://ftp.ensembl.org/pub/release-115/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.115.gff3.gz"
    check "C. elegans FASTA (WBcel235)"     "https://ftp.ensembl.org/pub/release-115/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz"
    check "CaeNDR WI.20250625 hard-filter"  "https://caendr-open-access-data-bucket.s3.us-east-2.amazonaws.com/dataset_release/c_elegans/20250625/variation/WI.20250625.hard-filter.isotype.vcf.gz"
    check "Human GFF3 (Ensembl 115)"        "https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz"
    check "Human FASTA (Ensembl 115)"       "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    check "GIAB HG002 v4.2.1 benchmark VCF" "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    check "GIAB HG002 v4.2.1 high-conf BED" "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
}

check_acmg() {
    echo -e "\n== ACMG concordance benchmark (data/benchmark) =="
    check "ClinVar GRCh38 VCF (NCBI)"        "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
    check "ClinVar variant_summary (NCBI)"   "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
    check "gnomAD v4.1 exomes chr22 (gs)"    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr22.vcf.bgz"
    check "gnomAD v4.1 gene constraint TSV"  "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
    check "REVEL v1.3 all chromosomes"       "https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip"
    check "ClinGen Gene-Disease Validity"    "https://search.clinicalgenome.org/kb/gene-validity/download"
}

case "${1:-all}" in
    speed) check_speed ;;
    acmg)  check_acmg ;;
    all)   check_speed; check_acmg ;;
    *)     echo "Usage: $0 [all|speed|acmg]"; exit 1 ;;
esac

echo ""
if [[ $fail -eq 0 ]]; then
    echo -e "${GREEN}All $ok URLs OK.${NC}"
    exit 0
fi
echo -e "${RED}${fail} of $((ok+fail)) URLs FAILED${NC} — update benchmarks/download_data.sh or the ACMG download script."
for u in "${failed_urls[@]}"; do echo "  - $u"; done
exit 1
