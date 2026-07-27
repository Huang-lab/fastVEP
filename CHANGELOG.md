# Changelog

All notable changes to fastVEP will be documented in this file. Dates are
ISO 8601. Format loosely follows [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Added

- **fastvep-sa / CLI**: `sa-build --format osa2` builds the v2 `.osa2` format
  for the numeric-payload sources — `--source gnomad`, `onekg` (`1000g`),
  `topmed`, and `alphamissense`. The v2 format (chunked ZIP, u32 value arrays,
  Var32 binary search) is smaller at genome scale and faster to query on the
  sparse, scattered lookups a real VCF generates — `bench_v1_vs_v2` measures
  v2 query throughput at ~3.8–4.5× v1 at 10M records on a strided workload,
  with ~33% smaller files (v2 build is a one-time ~4× slower). The gnomAD
  encoder reuses the v1 INFO field-name detection (v2/v3/v4/joint releases,
  hyphen-separated locals); 1000G and TOPMed reuse their existing v1 parsers
  through a generic `osa2_record_from_v1` bridge that re-encodes each record's
  flat-object JSON into the parallel u32 layout. All v2 builds stream one chunk
  at a time, so peak memory is bounded by the densest ~1 MB window rather than
  the whole file, and v1/v2 output parity is verified by test.

- **fastvep-sa / CLI**: new `--source alphamissense` builds AlphaMissense
  pathogenicity predictions (Cheng et al. 2023; Zenodo 8208688) from the
  genome-coordinate TSV releases into both `.osa` (v1) and `.osa2` (v2). Each
  allele carries a pathogenicity score plus a three-level class
  (`likely_benign` / `ambiguous` / `likely_pathogenic`) — a numeric-plus-small-
  categorical payload that builds natively into v2 (one u32 score column plus a
  u32 class-index column against a 3-entry string table), the first v2 source
  to use categorical string tables. AlphaMissense is first-class in every
  output: JSON (`alphaMissense`), tab, and VCF (`FV_ALPHAMISSENSE`, format
  `ALLELE|PATHOGENICITY|CLASS`). The v1 JSON is built through the very same
  `Field`/`format_value` code the v2 reader reconstructs with, so the two
  formats emit byte-identical annotations (verified by test).

- **fastvep-sa / CLI**: `--format osa2` now also supports the string/array
  sources that don't fit the numeric u32 layout — `--source dbsnp`,
  `--source cosmic`, and `--source clinvar`. Their whole-record JSON is stored
  as one opaque blob per variant (`raw_json_blob_fields` — a single JsonBlob
  field with an empty alias that the reader emits verbatim), so v2 output is
  byte-identical to v1 while v2's chunk-level zstd of the blob column shrinks
  the database sharply (~0.30× at genome scale per `bench_shapes`). dbSNP
  streams through the existing v1 parser via a `bridge_v1_raw_blobs` adapter;
  COSMIC and ClinVar are buffered and sorted by their v1 parsers, then bridged.
  ClinVar's nested significance/phenotype arrays and its `is_array` metadata
  survive intact. v1/v2 output parity is verified by test for all three.

- **fastvep-sa**: `bench_shapes` example measures v1 `.osa` vs v2 `.osa2`
  on-disk size across the payload *shapes* fastVEP sources carry (numeric,
  score+categorical, opaque id-string, array/blob). It answers "is v2 smaller
  for everything?" empirically: **no, not universally** — the answer is shape-
  and scale-dependent. At 2M records v2 is *larger* for numeric/score payloads
  (1.07–1.12×, fixed per-chunk/ZIP overhead dominates) but already much smaller
  for JSON-blob payloads (0.40×). At 10M records v2 is smaller or comparable
  across *all* shapes: numeric 0.67×, score 0.93×, and the blob shapes
  (dbSNP-/ClinVar-like) 0.30–0.34× — because v2 zstd-compresses a whole chunk's
  JSON blobs together, exploiting cross-record redundancy v1's per-block scheme
  can't. So v2 is a clear size win at genome scale, dramatically so for
  blob-heavy sources; it is not a win for small inputs.

### Fixed

- **fastvep-sa**: v2 (`.osa2`) reader re-opened the ZIP file and re-parsed the
  entire central directory on *every* chunk load, making genome-scale queries
  hundreds of times slower than v1. The archive is now parsed once at
  `open()` and reused, restoring the format's intended performance.
- **fastvep-sa**: v2 (`.osa2`) long variants (indels, ref+alt > 4 bases)
  returned the wrong values or were unreadable. The writer keyed each
  `LongVariant` to its input-order position while value columns held only
  short variants in Var32-sorted order, and `Chunk::is_empty` ignored long
  variants entirely (long-only chunks were skipped). Value columns now use a
  combined short-then-long layout with each long variant's slot recorded, and
  `is_empty` accounts for both. gnomAD is full of indels, so this is required
  for correct v2 annotation.

- **fastvep-web**: stored XSS via gene/transcript metadata (symbol, IDs,
  HGVS strings, supplementary-annotation values) rendered unescaped into
  the results table and ACMG modal; a crafted GFF3 `Name`/`ID` or
  supplementary-annotation string could execute in every viewer's browser.
  All such fields are now HTML-escaped before interpolation.
- **fastvep-web**: `/api/annotate` responses leaked internal error detail
  (file paths, parse internals) to clients; errors are now logged
  server-side only and clients get a generic message.
- **fastvep-io**: VCF lines with an empty REF field caused an integer
  underflow in end-coordinate calculation (panic in debug, silent
  corruption in release). Now rejected with a clear parse error.
- **fastvep-cache**: GFF3 lines with `start == 0` or `start > end` parsed
  as valid u64 coordinates but are invalid 1-based GFF3, risking
  downstream underflow in exon/CDS offset math. Now skipped with a
  warning, matching the existing non-numeric-coordinate guard.

### Changed

- **fastvep-web**: `/api/annotate` no longer holds a write lock on the
  shared `AnnotationContext` for the duration of annotation (previously
  needed only to toggle `acmg_config` per request); annotation now takes
  a read lock and passes the ACMG config as a call argument
  (`annotate_vcf_text_with_acmg`), so concurrent requests — including
  unrelated `/api/status` reads — no longer serialize behind whichever
  annotation is running, and one request's ACMG toggle can no longer
  clobber another's mid-flight.
- **fastvep-web**: stats persistence (`save_stats`) now runs on the
  blocking thread pool instead of the async request path, so it no
  longer stalls a tokio worker thread on disk I/O for every annotate call.
- **fastvep-web**: CORS now scopes `allow_methods`/`allow_headers` to
  what the API actually uses (GET/POST, `Content-Type`) instead of `Any`
  on every axis.

### Added

- Unit tests for `fastvep-annotate` and `fastvep-web` (previously zero in
  both crates despite being the code every request flows through), plus
  regression tests for the fixes above.

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
