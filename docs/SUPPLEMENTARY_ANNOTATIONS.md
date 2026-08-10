# Supplementary annotation (fastSA) output contract

This document is the authoritative schema for every supplementary annotation
source produced by `fastvep sa-build` and emitted by `fastvep annotate`. The
same per-source pipe format is used by the **VCF** `FV_*` INFO fields and by
the **tab** `FV_*` columns; the **JSON** output carries the same data as a
structured object under the source's JSON key.

All identifiers prefixed with `FV_` are owned by fastVEP. When you annotate
an input VCF that already declares one of these IDs, fastVEP strips the
input's `##INFO=<ID=FV_*>` headers and any existing `FV_*` values from each
record's INFO column before writing its own. Non-fastVEP INFO fields and the
input's other headers pass through unchanged.

## Identifiers across output formats

| Source            | `sa-build --source` | JSON key          | VCF INFO ID         | Tab column          | Scope        |
|-------------------|---------------------|-------------------|---------------------|---------------------|--------------|
| ClinVar           | `clinvar`           | `clinvar`         | `FV_CLINVAR`        | `FV_CLINVAR`        | Allele       |
| gnomAD            | `gnomad`            | `gnomad`          | `FV_GNOMAD`         | `FV_GNOMAD`         | Allele       |
| dbSNP             | `dbsnp`             | `dbsnp`           | `FV_DBSNP`          | `FV_DBSNP`          | Allele       |
| COSMIC            | `cosmic`            | `cosmic`          | `FV_COSMIC`         | `FV_COSMIC`         | Allele       |
| 1000 Genomes      | `onekg`             | `oneKg`           | `FV_1KG`            | `FV_1KG`            | Allele       |
| TOPMed            | `topmed`            | `topmed`          | `FV_TOPMED`         | `FV_TOPMED`         | Allele       |
| MitoMap           | `mitomap`           | `mitomap`         | `FV_MITOMAP`        | `FV_MITOMAP`        | Allele       |
| PhyloP            | `phylop`            | `phylop`          | `FV_PHYLOP`         | `FV_PHYLOP`         | Allele       |
| GERP              | `gerp`              | `gerp`            | `FV_GERP`           | `FV_GERP`           | Allele       |
| DANN              | `dann`              | `dann`            | `FV_DANN`           | `FV_DANN`           | Allele       |
| REVEL             | `revel`             | `revel`           | `FV_REVEL`          | `FV_REVEL`          | Allele       |
| PrimateAI         | `primateai`         | `primateAI`       | `FV_PRIMATEAI`      | `FV_PRIMATEAI`      | Allele       |
| dbNSFP            | `dbnsfp`            | `dbnsfp`          | `FV_DBNSFP`         | `FV_DBNSFP`         | Allele       |
| AlphaMissense     | `alphamissense`     | `alphaMissense`   | `FV_ALPHAMISSENSE`  | `FV_ALPHAMISSENSE`  | Allele       |
| SpliceAI          | `spliceai`          | `spliceAI`        | `SpliceAI`          | `SpliceAI`          | Allele       |
| OMIM / ClinGen GDV| `omim`              | `omim`            | `FV_OMIM`           | `FV_OMIM`           | Gene         |
| gnomAD constraint | `gnomad_genes`      | `gnomad_genes`    | `FV_GNOMAD_GENE`    | `FV_GNOMAD_GENE`    | Gene         |
| ClinVar protein   | `clinvar_protein`   | `clinvar_protein` | `FV_CLINVAR_PROTEIN`| `FV_CLINVAR_PROTEIN`| Gene         |
| Custom VCF        | `custom_vcf`        | `<--name>`        | `FV_<--NAME>`*      | `FV_<--NAME>`*      | Allele       |
| Custom BED        | `custom_bed`        | `<--name>`        | `FV_<--NAME>`*      | `FV_<--NAME>`*      | Interval     |
| Custom (auto)     | `custom`            | `<--name>`        | `FV_<--NAME>`*      | `FV_<--NAME>`*      | depends on input |

`SpliceAI` is intentionally **not** namespaced under `FV_*` to remain
compatible with the standard SpliceAI INFO contract that downstream tools
already parse.

\* For `custom_vcf` / `custom_bed` / `custom`, the JSON key and `FV_*`
INFO ID derive from the `--name` flag at build time (or the input
filename if `--name` is omitted). For example, `sa-build --source
custom_vcf --name clinical -i my.vcf -o my` produces a `.osa` whose
JSON key is `clinical` and whose VCF projection is emitted under
`FV_CLINICAL`.

## On-disk file formats

`sa-build` writes one of three binary container types depending on the
source scope:

| Extension | Scope        | Producers                                             | Loaded by `--sa-dir`? |
|-----------|--------------|-------------------------------------------------------|-----------------------|
| `.osa2`   | Allele-level | Every allele-level source (the `--format auto` default) | Yes (`Osa2Reader`)  |
| `.osa`    | Allele-level | Any allele-level source, via explicit `--format osa`  | Yes (`SaReader`)      |
| `.osi`    | Interval     | `custom_bed` (and any future structural-variant DB)   | Yes (`OsiReader`)     |
| `.oga`    | Gene         | `omim`, `gnomad_genes`, `clinvar_protein`             | Yes (`GeneIndex`)     |

All readers refuse malformed/oversized index payloads (the `.osa.idx`,
`.osi`, and `.oga` headers carry a `schema_version` and the payload
size is bounded — see `crates/fastvep-sa/src/common.rs::MAX_INDEX_PAYLOAD`).
A directory passed to `--sa-dir` is scanned for any of the four
extensions; mismatched files (e.g. a `.tsv` left in place) are silently
skipped.

## Choosing v1 vs v2 (`--format`)

Allele-level sources can build to either the v1 `.osa` or v2 `.osa2`
container. **You normally don't have to choose:** `sa-build --format` defaults
to `auto`, and every allele-level source now has a v2 encoder, so `auto` builds
`.osa2` for all of them.
The annotate side loads either transparently from `--sa-dir`, so downstream
usage is identical.

Rule of thumb:

- **`auto` (default)** — recommended. Builds `.osa2` for every allele-level
  source; the gene- and interval-level sources have no v2 form and build
  `.oga`/`.osi` as always.
- **`osa` (force v1)** — the escape hatch. Use for a faster one-time build, or
  when a downstream tool specifically expects the `.osa`/`.osa.idx` file pair.
- **`osa2` (force v2)** — same as `auto` for allele-level sources; errors out
  for the gene- and interval-level sources rather than silently building
  something else.

Why v2 is the higher-quality default where available: at genome scale it is
**smaller on disk** (≈0.67× for numeric payloads, as low as ≈0.30× for the
JSON-blob sources) and **faster to query** on the sparse, scattered lookups a
real VCF produces (≈3.8–4.5× at 10M records). The trade-off is a one-time
build that is ≈4× slower; v2 is also not a size win for very small inputs,
where its fixed per-chunk overhead dominates. Output is byte-identical between
the two formats for the blob-backed sources (dbSNP, COSMIC, ClinVar, REVEL,
PrimateAI, dbNSFP, MitoMap, `custom_vcf`, and the positional scores) and for
AlphaMissense and SpliceAI. For the
frequency sources (gnomAD, 1000G, TOPMed) the integer counts are exact and
frequencies match to a fixed 5e-7 resolution — identical in practice except
for the very rarest gnomAD v4 singletons, whose AF falls below that floor.

AlphaMissense and SpliceAI reach that byte-identity by having their **v1**
builder render each record through the same field-encode/format code the v2
reader reconstructs with.
For SpliceAI that changed the v1 `.osa` payload's number formatting: delta
scores now render in the shared scientific-notation form (`8.500000e-1` rather
than `0.85`).
The JSON keys and their order are unchanged, and every consumer parses these
with a JSON parser, so nothing downstream is affected — but a `.osa` built
before this change will differ textually from one built after.
The delta positions (`dpAg`…`dpDl`) are unaffected and still render as plain
signed integers.

Every allele-level and positional source has a v2 encoder: `gnomad`, `onekg`
(`1000g`), `topmed`, `alphamissense`, `spliceai`, `dbsnp`, `cosmic`, `clinvar`,
`revel`, `primateai`, `dbnsfp`, `mitomap`, `custom_vcf`, and the positional
per-base scores `phylop`, `gerp`, `dann`.
The only sources that ignore `--format` are the ones v2 structurally cannot
represent: gene-keyed sources (`omim`, `gnomad_genes`, `clinvar_protein`) build
`.oga`, and interval-keyed `custom_bed` builds `.osi`.

Two encodings sit behind `.osa2`.
**Column-oriented** sources (`gnomad`, `onekg`/`1000g`, `topmed`,
`alphamissense`, `spliceai`) decompose each record into parallel u32 value
columns, plus a categorical string table where one applies (AlphaMissense's
three-level class, SpliceAI's gene symbol).
The rest store one **whole-record JSON blob** per variant, because their
payloads are high-cardinality ID strings or nested arrays that the numeric
layout cannot represent; they still gain from v2's chunk-level zstd, which
compresses a whole chunk's records together.

Positional scores get an especially large v2 win: their per-base coordinates
delta-encode to almost nothing and the score column compresses well, so a
dense per-base database is roughly **0.23× the size** of the v1 `.osa`
(measured via `bench_shapes` on realistic-entropy synthetic scores).

Measured on real releases, whole-database, with every record queried through both
readers and compared:

| Source | Records | v1 (`.osa` + `.idx`) | v2 (`.osa2`) | Ratio |
|--------|---------|----------------------|--------------|-------|
| ClinVar (full release) | 4,438,232 | 50.5 MB | 40.9 MB | 0.81× |
| SpliceAI (chr22) | 827,827 | 4.8 MB | 2.7 MB | 0.56× |
| REVEL (chr22) | 1,776,286 | 8.7 MB | 3.8 MB | 0.44× |
| gnomAD (chr22) | 852,255 | 15.4 MB | 12.9 MB | 0.84× |
| PhyloP (chr22) | 852,255 | 3.1 MB | 1.4 MB | 0.45× |

Annotating 308,740 real chr22 variants against all five of those databases at
once takes **12.7 s from the `.osa` set and 8.5 s from the `.osa2` set (1.50×)**,
with md5-identical output.

The SpliceAI ratio is data-dependent and can be much better than the 0.56×
above: that shard comes from gnomAD's summarised SpliceAI fields, where the four
delta scores are equal per record and all four delta positions are zero, so v1's
JSON is unusually compressible. On a fixture with independently varying scores,
positions, and 200 gene symbols the ratio is **0.14×**.

### v1/v2 record fidelity

Two ways v2 used to disagree with v1 about *which records exist* were found by
running the full real ClinVar release and a real per-chromosome `--sa-dir`
through both formats and comparing every record. Both are fixed; both are worth
knowing about if you have `.osa2` files built before this change.

**Alleles containing `N` or IUPAC codes were dropped.** Var32 and kmer16 are both
two bits per base, so neither can encode them, and the writer skipped such
records with a `log::warn!` that the CLI installs no logger to display. Real
ClinVar carries 668 of them among 4,438,232 records (0.015%) — `Microsatellite`
and long `Insertion` entries. They are now kept in a third per-chunk bucket that
stores the allele bytes verbatim, costing 0.20% on the ClinVar archive
(82 KB of 41 MB).

**Duplicate keys resolved arbitrarily.** Several sources hold more than one
record for the same `(position, ref, alt)`: REVEL has one per transcript or
protein change (111,270 duplicate keys on chr22 alone, 6.4% of its records),
SpliceAI one per overlapping gene. The v1 reader scans forward from a position's
first entry and so always returns the first such record, but a v2 chunk index is
a sorted key array resolved by `binary_search`, which on a run of equal keys
returns an arbitrary member — so the two formats could return different records
for the same query, which showed up as differing REVEL scores in a real chr22
annotate run. The writer now keeps only the first record of each duplicate run,
which reproduces v1 exactly and shrinks the archive, since the dropped records
were unreachable through either reader.

If you have `.osa2` files built before this change, rebuild or `sa-convert` them.

### Converting an existing v1 database (`sa-convert`)

If you already have `.osa` files from before every source defaulted to v2, you
do not have to re-download and re-parse the upstream release:

```bash
fastvep sa-convert -i sa_databases/dbsnp.osa          # writes dbsnp.osa2
fastvep sa-convert -i old/spliceai.osa -o new/spliceai # explicit output base
```

The conversion preserves the record set and every record's JSON payload byte
for byte, along with the database's identity and matching semantics (`json_key`,
name, version, description, assembly, `match_by_allele`, `is_array`,
`is_positional`), so a converted database answers every query exactly as the
`.osa` did.

One caveat, which the tool prints: a `.osa` retains no field schema, so a
conversion can only produce the JSON-blob encoding. For the blob-backed sources
that is exactly what a native v2 build produces anyway. For the five
column-oriented sources it is not, so if you still have the upstream release
prefer `sa-build --format osa2` to get the column layout. Which of the two ends
up smaller is data-dependent, so compare if size matters.

`sa-convert` refuses `.osi` and `.oga` inputs with an explanation rather than
failing obscurely: v2 is a variant-level container and has no equivalent form
for interval- or gene-keyed data.

## Pipe formats

Each value is a pipe-delimited string. Multiple values for the same record
(for example, multiple alt alleles or multiple gene entries) are separated by
`,` in VCF and by `,` within the same tab cell. Empty fields render as the
empty string between two pipes (`A||C`).

Allele-level sources lead with the **uploaded ALT allele** (preserving the
original REF/ALT of the input VCF, especially for indels); gene-level sources
lead with the **gene symbol**.

### Allele-level

- `FV_CLINVAR`: `ALLELE|SIGNIFICANCE|REVIEW_STATUS|PHENOTYPES|VARIANT_CLASS|SO_ACCESSION`
- `FV_GNOMAD`: `ALLELE|ALL_AF|ALL_AC|ALL_AN|ALL_HC|AFR_AF|AMR_AF|ASJ_AF|EAS_AF|FIN_AF|MID_AF|NFE_AF|OTH_AF|REMAINING_AF|SAS_AF`
- `FV_DBSNP`: `ALLELE|ID|GLOBAL_MAF`
- `FV_COSMIC`: `ALLELE|ID|GENE|COUNT`
- `FV_1KG`: `ALLELE|ALL_AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF`
- `FV_TOPMED`: `ALLELE|ALL_AF|ALL_AC|ALL_AN`
- `FV_MITOMAP`: `ALLELE|DISEASE|STATUS`
- `FV_PHYLOP`: `ALLELE|SCORE`
- `FV_GERP`: `ALLELE|SCORE`
- `FV_DANN`: `ALLELE|SCORE`
- `FV_REVEL`: `ALLELE|SCORE`
- `FV_PRIMATEAI`: `ALLELE|SCORE`
- `FV_DBNSFP`: `ALLELE|SIFT|POLYPHEN`
- `FV_ALPHAMISSENSE`: `ALLELE|PATHOGENICITY|CLASS`
- `SpliceAI`: `ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL`

### Gene-level

- `FV_OMIM`: `SYMBOL|MIM_NUMBER|PHENOTYPES`
- `FV_GNOMAD_GENE`: `SYMBOL|PLI|LOEUF|MIS_Z|SYN_Z`
- `FV_CLINVAR_PROTEIN`: `SYMBOL|PROTEIN_VARIANTS` — the `PROTEIN_VARIANTS`
  segment is itself a `&`-joined list of `pos:ref>alt:significance` records.

## Escaping inside pipe fields

To keep `FV_*` values parseable by `bcftools` and similar tools without
double-decoding, fastVEP percent-encodes the following characters within any
pipe field:

| Character | Replacement |
|-----------|-------------|
| `:`       | `%3A`       |
| `;`       | `%3B`       |
| `=`       | `%3D`       |
| `%`       | `%25`       |
| `,`       | `%2C`       |
| `\r`      | `%0D`       |
| `\n`      | `%0A`       |
| `\t`      | `%09`       |
| ` ` (space)| `%20`       |
| `"`       | `%22`       |
| `\|`      | `%7C`       |
| `&`       | `%26`       |

Lists within a single pipe field (for example, multiple ClinVar
significances) are joined with `&` *after* per-element escaping, so the
delimiter cannot collide with payload content.

`bcftools query -f '%INFO/FV_CLINVAR\n'` returns the raw escaped value; a
single percent-decode pass recovers the original. JSON output is **not**
escaped this way — it carries the original strings unmodified inside a
structured object.

## Output-format behavior

- **VCF**: each loaded source emits one `##INFO=<ID=FV_*,Number=.,Type=String,Description="...">`
  header line and one `FV_*=<pipe value>` entry per record (omitted when the
  variant has no annotation). The header `Description` carries the exact
  pipe format above.
- **Tab**: each loaded source appends one column to the row, after the 17
  built-in columns. The file prologue contains one
  `## COLUMN=<ID=FV_*,Description="...">` line per loaded source documenting
  the pipe format. Empty cells render as `-`.
- **JSON**: each loaded source places its full structured payload under its
  JSON key (`clinvar`, `dbnsfp`, …) on the relevant transcript consequence
  (allele-level) or under the variant's `genes` map (gene-level). JSON
  output is the richest projection; VCF and tab are flattened views.

## Adding a new source

`crates/fastvep-io/src/output.rs` defines the single `VCF_PROJECTION_SPECS`
constant that drives VCF emission, tab columns, header documentation, and
input-VCF conflict stripping. Adding an entry there automatically:

- declares the `##INFO=<ID=FV_NEW,...>` header,
- appends the `FV_NEW` tab column,
- strips any pre-existing `FV_NEW` INFO field from the input VCF,
- and the unit test `tab_supplementary_column_names_match_vcf_header_order`
  asserts the column appears in tab output.

When you add a spec, update this document with the new row in the table
above and the new pipe format in the per-source list. The doc-coverage check
in CI fails if the schema in `output.rs` documents an INFO ID that is missing
from this file.
