#!/usr/bin/env bash
# Run fastvep with --acmg on the ClinVar 2-star+ VCF.
# Output JSON includes ACMG classification per variant.
set -euo pipefail

ROOT=/Users/kuan-lin.huang/Projects/fastVEP
GENOME=$ROOT/genomes/human
INPUT=$ROOT/analysis/acmg_benchmark/real_data/clinvar_2star.vcf
OUTPUT=$ROOT/analysis/acmg_benchmark/real_data/clinvar_2star.fastvep.jsonl

mkdir -p "$(dirname "$OUTPUT")"

if [ ! -f "$GENOME/Homo_sapiens.GRCh38.115.gff3" ]; then
  echo "Decompressing GFF3..."
  gunzip -k "$GENOME/Homo_sapiens.GRCh38.115.gff3.gz"
fi
if [ ! -f "$GENOME/Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]; then
  echo "Decompressing FASTA..."
  gunzip -k "$GENOME/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fi
if [ ! -f "$GENOME/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai" ]; then
  echo "Indexing FASTA..."
  samtools faidx "$GENOME/Homo_sapiens.GRCh38.dna.primary_assembly.fa" 2>/dev/null \
    || python3 -c "
import sys
fa = '$GENOME/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
fai = fa + '.fai'
with open(fa, 'rb') as f, open(fai, 'w') as o:
    name, ln, off, bp_line, bytes_line = None, 0, 0, 0, 0
    pos = 0
    for line in f:
        if line.startswith(b'>'):
            if name:
                o.write(f'{name}\t{ln}\t{name_off}\t{bp_line}\t{bytes_line}\n')
            name = line[1:].split()[0].decode()
            ln = 0
            name_off = pos + len(line)
            bp_line = 0
            bytes_line = 0
        else:
            stripped = line.rstrip(b'\n').rstrip(b'\r')
            ln += len(stripped)
            if bp_line == 0:
                bp_line = len(stripped)
                bytes_line = len(line)
        pos += len(line)
    if name:
        o.write(f'{name}\t{ln}\t{name_off}\t{bp_line}\t{bytes_line}\n')
"
fi

SA_DIR=$GENOME/sa_db
mkdir -p "$SA_DIR"

echo "==> Running fastvep annotate --acmg"
$ROOT/target/release/fastvep annotate \
  -i "$INPUT" \
  -o "$OUTPUT" \
  --gff3 "$GENOME/Homo_sapiens.GRCh38.115.gff3" \
  --fasta "$GENOME/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
  --sa-dir "$SA_DIR" \
  --acmg \
  --hgvs \
  --symbol \
  --canonical \
  --pick \
  --output-format json 2>&1 | tail -30

echo "Done. Output: $OUTPUT"
ls -la "$OUTPUT"
