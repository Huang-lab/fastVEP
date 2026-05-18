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
| SpliceAI          | `spliceai`          | `spliceAI`        | `SpliceAI`          | `SpliceAI`          | Allele       |
| OMIM / ClinGen GDV| `omim`              | `omim`            | `FV_OMIM`           | `FV_OMIM`           | Gene         |
| gnomAD constraint | `gnomad_genes`      | `gnomad_genes`    | `FV_GNOMAD_GENE`    | `FV_GNOMAD_GENE`    | Gene         |
| ClinVar protein   | `clinvar_protein`   | `clinvar_protein` | `FV_CLINVAR_PROTEIN`| `FV_CLINVAR_PROTEIN`| Gene         |

`SpliceAI` is intentionally **not** namespaced under `FV_*` to remain
compatible with the standard SpliceAI INFO contract that downstream tools
already parse.

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
