# Changelog

All notable changes to fastVEP will be documented in this file. Dates are
ISO 8601. Format loosely follows [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Added

- `fastvep cache --synonyms <chr_synonyms.txt>`: VEP-style chromosome
  synonym support, so a merged Ensembl + RefSeq cache built against a
  single FASTA reconciles accession seqids (`NC_000017.11`) with the
  FASTA's contig names (`17`). Transcripts are canonicalized to the FASTA
  naming at build time, so the merged cache uses one consistent scheme and
  `annotate` matches a VCF regardless of which GFF3 a transcript came from.
  Resolves issue #47.

### Changed

- Cache build and FASTA lookups now resolve `chr` ↔ bare and
  mitochondrial (`M`/`MT`/`chrM`/`chrMT`) contig aliases automatically
  (no synonyms file needed). `IndexedTranscriptProvider` applies the same
  aliasing at query time, so a `chr17` VCF matches a `17` cache.
- The "Chromosome … not found in FASTA index" error now suggests
  `--synonyms` when the missing name looks like a RefSeq accession.
- `chrom_aliases()` moved to `fastvep-core` (re-exported from
  `fastvep-sa::common`) so the cache builder and SA readers share one
  implementation.

## [0.2.0] — 2026-06-10

This release accumulates ~55 commits since v0.1.0, headlined by an
ACMG-AMP classification engine, custom annotation sources, VEP-style
merged-cache support, and a ~900× faster supplementary-annotation path.

### Added

#### ACMG-AMP variant classification (new `fastvep-classification` crate)

- `--acmg` flag on `fastvep annotate` runs full ACMG-AMP classification
  per Richards et al. 2015, with ClinGen-SVI–aligned criteria.
- Pathogenic criteria: PVS1 (Abou Tayoun 2018 decision tree),
  PS1 (incl. splice-RNA path), PS2 (de novo), PS3, PS4, PM1, PM2
  (inheritance-aware, ClinGen SVI v1.0), PM3 (v1.0 points-based scoring),
  PM5, PM6, PP1, PP2, PP3 (ClinGen SVI; Pejaver 2022 + Walker 2023, with
  anti-double-count against PVS1/PS1/PM5/PM1).
- Benign criteria: BA1 (Ghosh 2018 exception list), BS1, BS2, BS3, BS4,
  BP1, BP2, BP3, BP4 (splice gating), BP5, BP6, BP7 (Walker 2023
  exon-edge exclusion + deep-intronic extension).
- Trio / compound-het analysis: `--proband`, `--mother`, `--father`
  flags wire PS2 / PM6 / PM3 / BP2.
- Configurable thresholds via `--acmg-config <toml>`.
- ClinVar 2-star+ benchmark suite (`benchmarks/`); recall against P/LP
  reached 64% in v6 of the iteration and continues to improve.

#### Supplementary annotation (fastSA) sources & format

- Custom user-supplied annotations: `sa-build --source custom_vcf` /
  `--source custom_bed` / `--source custom` (auto-detect from input
  extension), with `--name` controlling the JSON-key / column name and
  `--info-fields` selecting which VCF INFO fields to extract. Custom
  BEDs produce a `.osi` interval-level database that is loaded
  alongside `.osa` / `.osa2` via `--sa-dir`. (#46, closes #43)
- gnomAD v4.1 *joint* VCF support (#41) — both per-chromosome and
  combined `joint` releases supported.
- Multi-allelic INFO splitting per VCF `Number=A` / `Number=R`
  semantics (custom_vcf); bi-allelic categoricals are kept whole.
- Gene-level annotations (`.oga`): wire-up for OMIM (ClinGen GDV),
  gnomAD constraint metrics, and ClinVar protein-position indices
  (#20).
- `--sa-only` mode emits only supplementary annotations, skipping the
  default CSQ pipeline — useful for re-annotation of already-annotated
  VCFs (#34).
- VCF-compatible INFO projections: `FV_CLINVAR`, `FV_GNOMAD`,
  `FV_DBSNP`, `FV_REVEL`, `FV_OMIM`, plus standard `SpliceAI` (#25).
- Supplementary annotation columns flow through tab output (#31).
- Reader hardening: refuse malformed/malicious `.osa.idx` / `.osi` /
  `.oga` payloads with bounded-size limits (#28).
- ~900× faster SA annotation via byte-budgeted LRU block cache plus
  per-variant deduplication (#33). Override budget via
  `FASTVEP_SA_CACHE_BYTES_PER_READER`.

#### Annotation pipeline

- VEP `--merged`-style cache: `--gff3` on `annotate` and `cache` is
  repeatable, supports `LABEL=path` syntax (auto-detects Ensembl /
  RefSeq from filenames; `gencode` → Ensembl; `GCF_*` / `refseq` →
  RefSeq), and emits per-transcript SOURCE labels through the CSQ /
  JSON / tab outputs side-by-side. (#46, closes #44)
- `fastvep cache` accepts multiple `--gff3` to pre-build a merged
  binary cache; `--transcript-cache` round-trips per-source labels.
- Gzipped VCF input for `fastvep annotate` (#21).

#### Output & CLI ergonomics

- Gene panel filter: `--gene-list <file>` keeps only tab rows whose
  transcript's gene_id or gene_symbol is on the list (#42).
- Explicit REF column for tab output via `--explicit-alleles` (#42).
- QC class column for tab output via `--qc-rules <toml>` (#42).
- `--pick` selection fixes: pre-filter to the surviving transcript
  before SA / ACMG passes (#40).

### Changed

- Cache builds are now deterministic — bit-for-bit reproducible
  across runs given the same GFF3 + FASTA inputs (#40).
- Custom-VCF INFO key iteration uses `BTreeMap` for stable JSON
  output across runs (content-hash reproducibility, #46).
- gnomAD annotations no longer drop records on `chr*`-style VCF input
  (#37/#38).
- ACMG criteria spec/impl alignment from per-criterion audit (#22).
- ACMG combiner defers conflict gating until rule resolution (#14).

### Fixed

- Path-traversal vulnerability in `resolve_genome_paths` for the web
  server (#5).
- Custom VCF parser handles CRLF line endings, JSON-special characters
  in INFO values (via `serde_json::Map` end-to-end), and flag-only
  INFO entries (stored as `"true"`) (#46).
- Custom BED parser tolerates CRLF, saturates `start = u32::MAX`
  instead of panicking, and skips malformed `end < start` records
  (#46).
- `OsiReader` resolves chromosome aliases (`chrM` / `M` / `MT` /
  `chrMT`) for both BED build and VCF query side (#46).
- Data corruption surfaced explicitly instead of silently masked as
  false negatives in SA readers (#35).
- Four real-data bugs in the ACMG classifier surfaced by the ClinVar
  2-star+ benchmark (#24).
- SA-build accepts gzipped inputs across all sources (#28).

### Documentation

- `docs/ACMG.md` — full ACMG-AMP methods writeup, ClinGen SVI–aligned.
- `docs/ACMG_SETUP.md` — per-source setup guide (REVEL, SpliceAI,
  PhyloP, dbNSFP, OMIM, ClinVar protein index, gnomAD gene constraint).
- `docs/SUPPLEMENTARY_ANNOTATIONS.md` — per-source FV_* / tab column /
  JSON-key schema.
- README rewrites for ACMG, multi-organism setup, merged cache, and
  custom annotation sources.
- Benchmarks reorganised under `benchmarks/`; URLs checked, scripts
  regrouped (#45).

### Internal

- New `fastvep-classification` crate (ACMG-AMP engine).
- New `fastvep-annotate` crate hosts the shared annotate pipeline used
  by both `fastvep annotate` (batch) and `fastvep-web` (per-request).
- CI workflow added for branch protection.
- 515 workspace tests at release (up from 233 at v0.1.0).

## [0.1.0] — 2026-04-23

Initial release. CLI (`fastvep`) and web server (`fastvep-web`); GFF3
gene-model loading; consequence prediction across 49 SO terms (incl.
SVs); HGVSg / HGVSc / HGVSp; allele-level supplementary annotations
(ClinVar, gnomAD, dbSNP, COSMIC, 1000 Genomes, TOPMed, MitoMap,
PhyloP, GERP, DANN, REVEL, SpliceAI, PrimateAI, dbNSFP); filter
engine; VCF / tab / JSON output.
