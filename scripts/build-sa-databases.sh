#!/usr/bin/env bash
# Download the publicly available supplementary-annotation sources for one
# assembly and convert each with `fastvep sa-build`.
#
#   scripts/build-sa-databases.sh --assembly GRCh38 --out sa_databases
#   scripts/build-sa-databases.sh --assembly GRCh37 --out sa37 --sources small
#   scripts/build-sa-databases.sh --assembly GRCh38 --sources gnomad,revel --chroms 22
#
# docs/SA_DATABASES.md is the manifest this script implements: it lists every
# source, which of them need an account (and so are not here), and what each
# download costs. Read it before running anything other than `small`.
#
# `sa-build` is a converter. A truncated download builds an empty database
# without an error, so every file is checked (gzip -t, then a record count that
# must be non-zero) before the build is trusted.

set -euo pipefail

ASSEMBLY=""
OUT="sa_databases"
WORK=""
SOURCES="small"
CHROMS=""
REGIONS=""
GNOMAD_RELEASE="exomes"
FASTVEP="${FASTVEP:-fastvep}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

SMALL="clinvar clinvar_protein gnomad_genes clingen repeatmasker"
LARGE="gnomad onekg alphamissense revel phylop dbsnp"

usage() {
    cat <<'EOF'
Usage: build-sa-databases.sh --assembly GRCh37|GRCh38 [options]

  --assembly A        GRCh37 or GRCh38 (required)
  --out DIR           where the finished databases go (default: sa_databases)
  --work DIR          where downloads go (default: DIR/../sa_sources)
  --sources LIST      comma-separated names, or `small` (default), `large`, `all`
                        small: clinvar clinvar_protein gnomad_genes clingen repeatmasker
                        large: gnomad onekg alphamissense revel phylop dbsnp
  --chroms LIST       restrict the per-chromosome sources (gnomad, onekg, revel,
                      phylop) to e.g. 21,22,X. Default: 1-22 and X
  --regions BED       gnomad only: fetch just these regions with tabix instead of
                      downloading whole chromosomes
  --gnomad-release R  exomes (default) or genomes. GRCh38 is v4.1, GRCh37 is v2.1.1
  --fastvep PATH      the fastvep binary (default: fastvep on PATH, or $FASTVEP)
  -h, --help          this text

alphamissense, dbsnp and clingen are whole-file sources: --chroms does not apply.
Sources that need an account are not here; see docs/SA_DATABASES.md.
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --assembly) ASSEMBLY="$2"; shift 2 ;;
        --out) OUT="$2"; shift 2 ;;
        --work) WORK="$2"; shift 2 ;;
        --sources) SOURCES="$2"; shift 2 ;;
        --chroms) CHROMS="$2"; shift 2 ;;
        --regions) REGIONS="$2"; shift 2 ;;
        --gnomad-release) GNOMAD_RELEASE="$2"; shift 2 ;;
        --fastvep) FASTVEP="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

case "$ASSEMBLY" in
    GRCh37|GRCh38) ;;
    *) echo "--assembly GRCh37 or GRCh38 is required" >&2; exit 2 ;;
esac
case "$GNOMAD_RELEASE" in exomes|genomes) ;; *) echo "--gnomad-release: exomes or genomes" >&2; exit 2 ;; esac
case "$SOURCES" in
    small) SOURCES="$SMALL" ;;
    large) SOURCES="$LARGE" ;;
    all) SOURCES="$SMALL $LARGE" ;;
    *) SOURCES="${SOURCES//,/ }" ;;
esac
[ -n "$CHROMS" ] || CHROMS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X"
CHROMS="${CHROMS//,/ }"

command -v "$FASTVEP" >/dev/null || { echo "fastvep not found: $FASTVEP" >&2; exit 1; }
mkdir -p "$OUT"
[ -n "$WORK" ] || WORK="$(cd "$OUT" && cd .. && pwd)/sa_sources"
mkdir -p "$WORK"
OUT="$(cd "$OUT" && pwd)"
cd "$WORK"

BUILT=""
say() { printf '\n== %s\n' "$*"; }

# fetch URL FILE: resumable, skipped when already present, gzip-checked.
fetch() {
    local url="$1" file="$2"
    if [ ! -s "$file" ]; then
        curl -fL -sS --retry 3 -C - -o "$file" "$url"
    fi
    case "$file" in
        *.gz|*.bgz) gzip -t "$file" || { echo "truncated download: $file (delete it and re-run)" >&2; exit 1; } ;;
    esac
}

# build SOURCE INPUT OUTNAME [extra sa-build args]: refuses a database of 0 records.
build() {
    local source="$1" input="$2" name="$3"; shift 3
    local log n
    log="$("$FASTVEP" sa-build --source "$source" -i "$input" -o "$OUT/$name" \
        --assembly "$ASSEMBLY" --no-progress "$@" 2>&1)" || { echo "$log" >&2; exit 1; }
    # Allele sources report "(N records)"; interval sources report "Parsed N".
    n="$(printf '%s\n' "$log" | grep -Eo '\(([0-9,]+) (records|genes)\)|Parsed [0-9,]+' | tail -1 | grep -Eo '[0-9][0-9,]*' | tr -d , || true)"
    if [ -z "$n" ] || [ "$n" -eq 0 ]; then
        echo "$name built 0 records from $input: the input did not match --assembly $ASSEMBLY, or it is empty" >&2
        printf '%s\n' "$log" >&2
        exit 1
    fi
    echo "   $name: $n records"
    BUILT="$BUILT $name"
}

# Per-release naming. gnomAD 2.1.1 (GRCh37) uses bare contig names, 4.1 uses chr*.
if [ "$ASSEMBLY" = GRCh38 ]; then
    GN_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/$GNOMAD_RELEASE/gnomad.$GNOMAD_RELEASE.v4.1.sites"
    GN_CHR_PREFIX="chr"
    GN_CONSTRAINT="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
    CLINVAR_DIR="vcf_GRCh38"
    UCSC="hg38"
    AM="hg38"
    DBSNP="GCF_000001405.40"
    PP_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.100way.phyloP100way"
else
    GN_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/$GNOMAD_RELEASE/gnomad.$GNOMAD_RELEASE.r2.1.1.sites"
    GN_CHR_PREFIX=""
    GN_CONSTRAINT="https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
    CLINVAR_DIR="vcf_GRCh37"
    UCSC="hg19"
    AM="hg19"
    DBSNP="GCF_000001405.25"
    PP_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way"
fi

src_clinvar() {
    say "ClinVar ($ASSEMBLY)"
    local f="clinvar_$ASSEMBLY.vcf.gz"
    fetch "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/$CLINVAR_DIR/clinvar.vcf.gz" "$f"
    # The PM2 frequency backstop reads these. A stripped ClinVar lacks them.
    # `head` closes the pipe early, which makes gzip fail under pipefail, hence `|| true`.
    [ "$(gzip -dc "$f" 2>/dev/null | head -n 500 | grep -c '^##INFO=<ID=AF_EXAC' || true)" -gt 0 ] \
        || { echo "$f has no AF_EXAC header: not the NCBI release" >&2; exit 1; }
    build clinvar "$f" clinvar
}

src_clinvar_protein() {
    say "ClinVar protein and splice index ($ASSEMBLY)"
    fetch "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz" variant_summary.txt.gz
    build clinvar_protein variant_summary.txt.gz clinvar_protein
}

src_gnomad_genes() {
    say "gnomAD gene constraints ($ASSEMBLY)"
    local f="gnomad_constraint_$ASSEMBLY.${GN_CONSTRAINT##*.}"
    fetch "$GN_CONSTRAINT" "$f"
    # The GRCh37 file is BGZF-compressed text; sa-build reads gzip by magic bytes.
    build gnomad_genes "$f" gnomad_genes
}

src_clingen() {
    say "ClinGen Gene-Disease Validity (json_key omim)"
    fetch "https://search.clinicalgenome.org/kb/gene-validity/download" clingen_gene_validity.csv
    python3 "$REPO/analysis/acmg_benchmark/scripts/sa_sources/clingen_gdv_to_oga.py" \
        clingen_gene_validity.csv clingen_gdv.tsv
    build omim clingen_gdv.tsv omim
}

src_repeatmasker() {
    say "RepeatMasker ($UCSC)"
    fetch "https://hgdownload.soe.ucsc.edu/goldenPath/$UCSC/database/rmsk.txt.gz" "rmsk_$UCSC.txt.gz"
    python3 "$REPO/analysis/acmg_benchmark/scripts/sa_sources/repeatmasker_to_bed.py" \
        "rmsk_$UCSC.txt.gz" > "repeatmasker_$UCSC.bed"
    # The name is load-bearing: the classifier finds this track by `repeat`.
    build custom_bed "repeatmasker_$UCSC.bed" repeatmasker --name repeatmasker
}

src_gnomad() {
    say "gnomAD $GNOMAD_RELEASE ($ASSEMBLY)"
    local c url f
    for c in $CHROMS; do
        url="$GN_BASE.chr$c.vcf.bgz"
        [ "$ASSEMBLY" = GRCh38 ] || url="$GN_BASE.$c.vcf.bgz"
        if [ -n "$REGIONS" ]; then
            command -v tabix >/dev/null || { echo "--regions needs tabix" >&2; exit 1; }
            f="gnomad_${GNOMAD_RELEASE}_${ASSEMBLY}_chr$c.region.vcf"
            tabix -h "$url" "${GN_CHR_PREFIX}$c:1-1" > "$f"
            # BED is 0-based; accept either contig spelling.
            awk -v a="chr$c" -v b="$c" -v p="$GN_CHR_PREFIX" '$1==a||$1==b{print p b":"$2+1"-"$3}' "$REGIONS" \
                | xargs -n 200 tabix "$url" >> "$f"
        else
            f="gnomad_${GNOMAD_RELEASE}_${ASSEMBLY}_chr$c.vcf.bgz"
            fetch "$url" "$f"
        fi
        build gnomad "$f" "gnomad_chr$c"
    done
}

src_onekg() {
    say "1000 Genomes ($ASSEMBLY)"
    local c url f
    for c in $CHROMS; do
        [ "$c" = X ] && { echo "   skipping chrX: the 1000G chrX files use a different name; see docs/SA_DATABASES.md"; continue; }
        if [ "$ASSEMBLY" = GRCh38 ]; then
            url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr$c.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        else
            url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$c.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        fi
        f="onekg_${ASSEMBLY}_chr$c.sites.vcf"
        # Only columns 1-8 are used, so the genotype columns are dropped on the
        # way down instead of being stored.
        if [ ! -s "$f" ]; then
            curl -fL -sS --retry 3 "$url" | gzip -dc | cut -f1-8 > "$f.part"
            mv "$f.part" "$f"
        fi
        build onekg "$f" "onekg_chr$c"
    done
}

src_alphamissense() {
    say "AlphaMissense ($AM)"
    fetch "https://zenodo.org/records/8208688/files/AlphaMissense_$AM.tsv.gz" "AlphaMissense_$AM.tsv.gz"
    build alphamissense "AlphaMissense_$AM.tsv.gz" alphamissense
}

src_revel() {
    say "REVEL v1.3 ($ASSEMBLY)"
    fetch "https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip" revel-v1.3_all_chromosomes.zip
    local c
    for c in $CHROMS; do
        # One database per chromosome bounds peak memory. The position column
        # follows --assembly; the file carries hg19 and GRCh38 side by side.
        unzip -p revel-v1.3_all_chromosomes.zip | awk -F, -v c="$c" 'NR==1||$1==c' > "revel_chr$c.csv"
        build revel "revel_chr$c.csv" "revel_chr$c"
        rm -f "revel_chr$c.csv"
    done
}

src_phylop() {
    say "PhyloP 100-way ($UCSC)"
    local c
    for c in $CHROMS; do
        fetch "$PP_URL/chr$c.phyloP100way.wigFix.gz" "phylop_${UCSC}_chr$c.wigFix.gz"
        build phylop "phylop_${UCSC}_chr$c.wigFix.gz" "phylop_chr$c"
    done
}

src_dbsnp() {
    say "dbSNP $DBSNP ($ASSEMBLY): about 28-30 GB"
    fetch "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/$DBSNP.gz" "$DBSNP.gz"
    build dbsnp "$DBSNP.gz" dbsnp
}

for s in $SOURCES; do
    case "$s" in
        clinvar|clinvar_protein|gnomad_genes|clingen|repeatmasker|gnomad|onekg|alphamissense|revel|phylop|dbsnp)
            "src_$s" ;;
        *) echo "unknown source: $s" >&2; exit 2 ;;
    esac
done

printf '\nBuilt into %s:%s\n' "$OUT" "$BUILT"
echo "Next: confirm each source answers a query (docs/SA_DATABASES.md, 'Check it answers')."
