# Supplementary annotation (fastSA) output contract

To build the databases themselves (URLs, sizes, GRCh37 and GRCh38), see [SA_DATABASES.md](SA_DATABASES.md).

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

\* For `custom_vcf` / `custom_bed` / `custom`, the JSON key and `FV_*` INFO ID derive from the `--name` flag at build time (or the input filename if `--name` is omitted).
For example, `sa-build --source custom_vcf --name clinical -i my.vcf -o my` produces a `.osa` whose JSON key is `clinical` and whose VCF projection is emitted under `FV_CLINICAL`.

A VCF INFO ID must match `[A-Za-z_][0-9A-Za-z_.]*` and `--name` is free text, so every character outside that set becomes `_` and the name is upper-cased.
Distinct names can therefore sanitise alike (`my-panel` and `my_panel` both give `FV_MY_PANEL`), and a derived ID can land on a built-in one (`--name gnomad_gene` gives `FV_GNOMAD_GENE`, which `gnomad_genes` owns).
Two INFO fields sharing an ID is a malformed VCF, so a collision takes a `_2`, `_3`, … suffix in load order; the built-in keeps its own ID, and each header line's `Description` names the source it came from.
The JSON key is always the unmodified `--name`.

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
| REVEL (chr22) | 1,776,286 | 9.1 MB | 4.0 MB | 0.44× |
| gnomAD (chr22) | 852,255 | 15.4 MB | 12.9 MB | 0.84× |
| PhyloP (chr22) | 852,255 | 3.1 MB | 1.4 MB | 0.45× |

The ClinVar and REVEL rows were measured again on 2026-11-07 with 0.5.0 (`benchmarks/run_sa_formats.py`), in megabytes of 10^6 bytes; the SpliceAI, gnomAD and PhyloP rows come from an earlier measurement, whose sources are not among the files that the repository's scripts fetch, and their unit was not recorded (the REVEL chr22 files of this run measure 8.7 and 3.8 in mebibytes, the numbers that row had, so the older rows are probably mebibytes).

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

### Databases that hold nothing

A database with no record opens without an error and answers every query with "no record", which for a frequency source reads as "absent from the population".
The dbSNP database that was meant to serve GRCh37 on fastvep.org was ten bytes of header: it had been built from a release whose contigs are RefSeq accessions, none matched the assembly, and the builder wrote it and exited with success.
So the built-in sources refuse to write one.
`sa-build` (and `sa-convert`, which would carry an empty `.osa` over) stops with an error that names the usual causes, removes what it had begun to write, and exits non-zero; `--allow-empty` writes the database anyway, for the case where nothing is the right answer, such as the shard of a chromosome that has no record in the regions being built.
`custom_vcf` in its default format (`.osa2`) refuses like the others; the older behaviour of a warning and an empty file remains for `custom_bed` and for `custom_vcf` under `--format osa`, so that a user can iterate on a filter.

A loader refuses to serve one that already exists, whoever built it.
`fastvep annotate` prints `warning: <file> was not loaded: it holds no records ...` on stderr for each database file that is empty or cannot be opened (a truncated copy is the common one), because `tracing` has no subscriber in the command line and the file would otherwise be a source that is silently missing from every row.
`fastvep-web` logs the same and reports it in the ACMG verdict, where the source counts as missing.

#### Chunk width for dense sources (`--chunk-bits`)

A chunk is the unit a query loads, and it covers a fixed genomic width -- 2^20
bases, 1 Mb, by default.
A blob-encoded chunk stores its records' JSON one after another in a single
compressed stream, so reading one record decompresses the whole chunk's worth,
and a column-encoded one inflates and indexes every record of the chunk to
answer for one.
The width is therefore the price of a lookup that misses the cache, and how
many records a megabase holds decides whether that price is small.

It is small for the sparse sources.
ClinVar and COSMIC put a few thousand records in a megabase, and the wide chunk
is what keeps them compact.
It is not small for a source with hundreds of thousands, and a source does not
have to score every base to be one.
dbSNP is dense, 9.4 GB for the 3,100 megabases of GRCh38, which is about 3 MB of
compressed records in each megabase: on the hosted site, 1,000 scattered GRCh38
variants took 24 s with it loaded and 0.5 s without it (23 ms a lookup, and the
same on a second pass, because the shared chunk cache is 24 MiB for each core
and holds a few dozen chunks of that size of the 3,100).
SpliceAI is the extreme: three alternates at essentially every position of a
gene body, which is hundreds of megabytes of JSON in a single 1 Mb chunk, so
that queries slow to a crawl and, past 256 MiB, the chunk exceeds the blob
decompression limit outright and *every* lookup in that megabase fails.

Pass `--chunk-bits 16` (65 kb) for such a source, to `sa-build` or, for a
database that is already built as a v1 `.osa`, to `sa-convert`:

```bash
fastvep sa-build --source dbsnp -i GCF_000001405.40.gz -o dbsnp --assembly GRCh38 --chunk-bits 16
fastvep sa-convert -i old/spliceai.osa -o new/spliceai --chunk-bits 16
```

`sa-build` applies it to every `.osa2` it writes and says so, with a warning,
when the build writes something that has no chunks (a v1 file, an interval or a
gene database); without the option each source keeps its own width.
`scripts/build-sa-databases.sh --chunk-bits 16` passes it to the allele-level
sources.

Measured on gnomAD v4.1 exomes, chromosome 21 (2,188,842 records), 3,041
variants in random order under a cache too small to hold one chunk, so that
every lookup loads its chunk:

| `--chunk-bits` | file size | build | 3,041 lookups | CPU of the lookups |
|---|---:|---:|---:|---:|
| 20 (default) | 20.7 MB | 153 s | 13.2 s | 42 s |
| 16 | 22.7 MB (+9.7%) | 156 s | 4.2 s | 7.3 s |
| no database | | | 1.9 s | |

The two databases return the same annotation for every variant (the sorted
outputs have the same checksum).
With the variants in position order, and so one chunk loaded for many, the two
took the same time as no database at all: the width matters when the variants
are scattered over many chunks, which is what a genome-wide VCF is.

Measured on dbSNP build 157 for GRCh38 (1,391,478,698 records) on the hosted site, 4 vCPUs:
the file at 16 bits is 9.19 GB, and 4,000 scattered variants went through a whole `annotate` run, process start included, in 18 s,
against 106 s through the file of 1 Mb chunks that the builder of July had written from the release that was then the latest (9.40 GB), which is 5.9 times faster.
The two files return the same RS number for 3,883 of the 3,884 variants that have a record.
The other variant has two current RS numbers in dbSNP for one allele; the older file returned one of them arbitrarily, and the writer now keeps the first record of the source for any run of equal keys, at every width.

Measured on a source built to SpliceAI's density, 300 scattered queries:

| `--chunk-bits` | file size | 300 queries | records annotated |
|---|---:|---:|---:|
| 20 (default) | 0.26 MB | 9.7 s | 5 / 300 |
| 18 | 0.27 MB | 0.83 s | 300 / 300 |
| 16 | 0.29 MB | 0.51 s | 300 / 300 |
| 14 | 0.35 MB | 0.47 s | 300 / 300 |

For comparison, the same queries against the native column encoding took 0.41 s,
so a converted database at `--chunk-bits 16` is close to it.

Do not lower it for the sparse sources: narrower chunks mean more of them, and
real ClinVar grows from 40.9 MB to 53.6 MB (and from 7,566 to 54,911 archive
entries) between 20 and 16.
`chunk_bits` is recorded per file, so a database built at any width is read
correctly without telling the reader anything.

A run that hits the limit now says so, once per source, rather than returning
records with the column quietly absent:

```
warning: spliceAI could not be read for 295 variant lookup(s); those variants
carry no spliceAI annotation, which is not the same as having none. First error:
SpliceAI: the JSON blobs for chunk chr1/1 decompress past the 256 MiB limit ...
```

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

### Custom sources

A custom source has no fixed pipe layout, because it has no fixed field set: with `--info-fields` unset, `sa-build` stores every INFO key each record happens to carry, so the fields vary from record to record.
Its projection is therefore **self-describing** rather than positional:

- `FV_<NAME>`: `ALLELE|KEY=VALUE&KEY=VALUE` for an allele-level source (`custom_vcf`, `custom_bed`)
- `FV_<NAME>`: `SYMBOL|KEY=VALUE&KEY=VALUE` for a gene-level source

The pairs are key-sorted, so a source lists its fields in the same order on every record.
Keys and values are escaped by the table below before assembly, so neither can contain `=`, `&`, `|` or `,` and the list splits back apart exactly: split the value on `|` for the subject and the pair list, the pair list on `&`, and each pair on its first `=`.

An array value repeats its key (`AF=0.1&AF=0.2`) rather than joining the elements, so no element is left without one.
Nested objects are dropped - neither custom builder emits one.
A record with no fields still emits the subject and an empty pair list (`G|`), because for a three-column BED "the allele is in a listed region" is the whole annotation and dropping the entry would be indistinguishable from a miss.

A frequency panel is the common case, and needs nothing special: build it with `--source custom_vcf` and the `AC` / `AN` / `AF` INFO fields come through as `FV_<NAME>=G|AC=7&AF=0.005814&AN=1204`, readable with `bcftools query -f '%INFO/FV_<NAME>\n'`.
Pass `--info-fields AC,AN,AF` at build time to keep only those three and get the same three pairs on every record.

## Escaping inside pipe fields

To keep `FV_*` values parseable by `bcftools` and similar tools without
double-decoding, fastVEP percent-encodes the following characters within any
pipe field:

| Character | Replacement | Why                                            |
|-----------|-------------|------------------------------------------------|
| `;`       | `%3B`       | ends an INFO field                             |
| `=`       | `%3D`       | separates an INFO key from its value           |
| `,`       | `%2C`       | separates the values of one INFO key           |
| `%`       | `%25`       | introduces an escape                           |
| `\|`      | `%7C`       | separates the pipe field's own segments        |
| `&`       | `%26`       | separates the elements of a list segment       |
| ` ` (space)| `%20`      | keeps the INFO column free of whitespace       |
| `\"`      | `%22`       | opens a quoted field in `read.delim` / `csv`   |
| `\r`      | `%0D`       | would end the record                           |
| `\n`      | `%0A`       | would end the record                           |
| `\t`      | `%09`       | would end the column                           |

Nothing else is encoded.
VCF 4.3 (data lines, fixed fields) restricts only `;`, `=` and `,` inside an INFO value; `|` and `&` are the delimiters these fields are assembled from; `%` has to be encoded for a single decode pass to be exact; and the last three would end the value, the column or the record.

The one character that used to be encoded here and is not any more is the colon.
It is reserved by nothing, Ensembl VEP writes it through literally, and fastVEP's own `CSQ` column already did too - every `HGVSc` carries one - so encoding it here contradicted the same record's `CSQ` field.
`validation/run_escaping_probe.sh` measures what VEP does with each of these characters against a real VEP 115.1, so the rules above can be rechecked rather than taken on trust.

Space and `"` are encoded although the spec permits both, and for the same reason in each case: there is no `FV_*` field in VEP to be diffed against, so matching VEP buys nothing here, while a reader of the column loses something real.
`awk` and `cut -d' '` split a record on the space; `R`'s `read.delim` and Python's `csv` both treat `"` as opening a quoted field and swallow the rest of the line.
`%20` is also reversible, where VEP's `_` substitution is not: `Breast_cancer` cannot be told back from `Breast cancer`.

Lists within a single pipe field (for example, multiple ClinVar
significances) are joined with `&` *after* per-element escaping, so the
delimiter cannot collide with payload content.
A field that introduces a delimiter of its own encodes that delimiter in its own leaves for the same reason.
`FV_CLINVAR_PROTEIN` is the only one that does: its `:` and `>` are literal between the leaves of a `pos:ref>alt:significance` record, and encoded (`%3A`, `%3E`) inside any leaf that contains one.

`bcftools query -f '%INFO/FV_CLINVAR\n'` returns the raw escaped value; a
single percent-decode pass recovers the original. JSON output is **not**
escaped this way — it carries the original strings unmodified inside a
structured object.

## Output-format behavior

- **VCF**: each loaded source emits one `##INFO=<ID=FV_*,Number=.,Type=String,Description="...">`
  header line and one `FV_*=<pipe value>` entry per record (omitted when the
  variant has no annotation). The header `Description` carries the exact
  pipe format above. The declared headers, the owned-ID list that a
  re-annotation replaces, and the emitted values are all derived from one
  loaded-source lookup, so a run cannot emit an `FV_*` field it did not
  declare.
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
