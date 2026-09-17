# fastVEP REST API

`fastvep-web` exposes a small JSON HTTP API over the same annotation engine the CLI uses.
It is meant to be **self-hosted** on your own workstation, lab server, or local network.

## Please self-host

The hosted instance at [fastVEP.org](https://fastVEP.org) is sized for interactive browser use: someone pastes a handful of variants and reads the table.
It is a single small machine with no rate limiting, so it is not a public API endpoint and there is no API key to request.
Please do not point scripts, pipelines, or services at it.

Self-hosting is the better option for you regardless:

- Your variant data never leaves your network, which for patient data is usually a requirement rather than a preference.
- No shared instance to compete with, and no dependency on our uptime.
- You choose the gene model release and which supplementary databases are loaded, instead of inheriting ours.

## What the API is for

The server loads the gene model, reference FASTA, and supplementary annotation databases once at startup and keeps them resident.
A full human gene model is over 500,000 transcripts; annotating one variant by invoking the CLI pays that load cost on every call, while the server pays it once.

That makes the API the right tool for **small numbers of variants, queried repeatedly**: a curation interface, a LIMS integration, an internal service, a notebook.

It is the wrong tool for bulk.
fastVEP's throughput advantage is in whole-genome files, and you do not POST those over HTTP.
The default request body limit is 10 MiB, and raising it is not the intended path.
For exome- or genome-scale VCFs use `fastvep annotate` directly.

## Setup

### 1. Install the server

```bash
cargo install --path crates/fastvep-web
```

### 2. Get reference data

The API needs the same inputs as the CLI: a GFF3 gene model, and a reference FASTA to resolve coding consequences.
Follow [Step 1 of the Local Setup Guide](../README.md#step-1-download-reference-data) in the README, then come back here.

### 3. Run the server

For a single user on one machine, bind to loopback so nothing else on the network can reach it:

```bash
fastvep-web \
  --bind 127.0.0.1 \
  --port 8080 \
  --gff3 data/Homo_sapiens.GRCh38.115.gff3 \
  --fasta data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --sa-dir data/sa/
```

To serve a lab or local network, bind to all interfaces (this is the default):

```bash
fastvep-web --bind 0.0.0.0 --port 8080 --gff3 ... --fasta ...
```

> **There is no authentication.**
> `--bind` is the only access control the server has, and it defaults to `0.0.0.0`.
> On a trusted local network that is fine, which is the deployment this API is designed for.
> On a machine with a public IP it means the server is open to the internet, so put it behind a reverse proxy with auth, or a firewall rule, before exposing it.
> [DEPLOYMENT.md](../DEPLOYMENT.md) covers the nginx and TLS setup.

A firewall does not cover the other risk, so the server closes that one itself.
`POST /api/upload-gff3` and `POST /api/load-genome` replace the gene model for the whole process, which on a server several colleagues share means one client's experiment changes what every other client is answered against.
Both clients are authorised and neither is doing anything wrong, which is exactly why being inside the firewall does not help.
So both endpoints are **off by default** and return `403` unless the server was started with `--allow-model-replacement`.
Turn it on for a single-user desktop run, where the browser interface's own genome switching needs it.
Leave it off for anything shared, and the gene model becomes a property of how the server was started.

Every `/api/annotate` response also reports the model that produced it, so a client that does not control the server can still tell.
See [`POST /api/annotate`](#post-apiannotate).

### 4. Verify it is up

```bash
curl -s http://localhost:8080/api/status | python3 -m json.tool
```

```json
{
    "allow_model_replacement": false,
    "backend": true,
    "gene_model_generation": 0,
    "gff3_source": "Homo_sapiens.GRCh38.115.gff3",
    "has_fasta": true,
    "sa_sources": ["ClinVar", "COSMIC", "dbSNP", "gnomAD", "1000 Genomes", "SpliceAI"],
    "status": "ok",
    "total_genomes": 0,
    "total_variants": 0,
    "transcripts": 509650,
    "version": "0.3.0"
}
```

Check `transcripts`, `has_fasta`, and `sa_sources` here before you trust any annotation.
Each of the three can be missing without the server saying so: a zero transcript count makes every variant `intergenic_variant`, and the other two are the subject of the next section.

`allow_model_replacement` says whether a request can change the gene model on this server; `gene_model_generation` counts how many times one has.
A UI can read the first to decide whether to offer its model-switching controls at all, rather than finding out with a `403` a user cannot interpret.

## Always pass `--fasta`

`--fasta` is optional to the argument parser and effectively required for correct answers.
Without reference sequence the server cannot read the codon a variant falls in, so a coding change degrades to a generic term instead of failing.
Annotating `17:43124027 G>A` against the same gene model, with and without `--fasta`:

| `--fasta` | `most_severe_consequence` |
| --- | --- |
| absent | `coding_sequence_variant` |
| present | `synonymous_variant` |

Neither response is an error, and both look plausible.
That is the failure mode worth guarding against here, so confirm `"has_fasta": true` in `/api/status` after every deployment and gene-model change.

The same applies to `sa_sources`: an empty list means ClinVar, gnomAD, and the rest were not loaded, and any downstream logic reading those fields silently sees nothing.

One wrinkle when you are testing this: `fastvep-web` writes a sidecar transcript cache next to the GFF3, and a run that had `--fasta` saves the built sequences into it.
A later run of the same gene model *without* `--fasta` then loads those sequences from the cache and returns correct coding terms while still reporting `"has_fasta": false`.
So a server that looks right on a machine you have already used a FASTA on can degrade on a fresh deployment, where no such cache exists.
Judge a deployment by what it returns on the target machine, not by what the same flags produced on yours.

## Endpoints

All request and response bodies are JSON, except `POST /api/upload-gff3`, which takes a raw GFF3 body.

| Method | Path | Purpose |
| --- | --- | --- |
| `GET` | `/api/status` | Server version and what is loaded |
| `POST` | `/api/annotate` | Annotate VCF records |
| `GET` | `/api/genomes` | List genomes available under `--data-dir` |
| `POST` | `/api/load-genome` | Switch the active genome - needs `--allow-model-replacement` |
| `POST` | `/api/upload-gff3` | Replace the active gene model with a posted GFF3 - needs `--allow-model-replacement` |

The two marked endpoints replace the gene model for every client of the server, so they return `403` unless it was started with `--allow-model-replacement`.

`GET /` serves the browser interface.

### `POST /api/annotate`

| Field | Type | Default | Meaning |
| --- | --- | --- | --- |
| `vcf` | string | required | VCF text. Header lines are optional; records are tab-separated. |
| `pick` | bool | `false` | Reduce each variant to one transcript consequence (see below). |
| `pick_allele` | bool | `false` | Reduce to one consequence per (variant, allele). VEP's `--pick_allele`. |
| `pick_allele_gene` | bool | `false` | Reduce to one consequence per (variant, allele, gene). VEP's `--pick_allele_gene`. |
| `flag_pick` | bool | `false` | Return every consequence, with `"pick": 1` on the one `pick` would have kept. VEP's `--flag_pick`. |
| `flag_pick_allele` | bool | `false` | As `flag_pick`, marking the pick per allele. |
| `flag_pick_allele_gene` | bool | `false` | As `flag_pick`, marking the pick per allele and gene. |
| `acmg` | bool | `false` | Attach ACMG-AMP classification to each consequence. |

At most one of the six pick switches may be set; the request is read in the order listed.

```bash
curl -s -X POST http://localhost:8080/api/annotate \
  -H 'Content-Type: application/json' \
  -d '{"vcf": "17\t43124027\t.\tG\tA\t50\tPASS\t.", "acmg": false}'
```

Every response also reports the gene model that produced it:

```json
"gene_model": {"gff3_source": "Homo_sapiens.GRCh38.115.gff3", "transcripts": 509650, "generation": 0}
```

This is the field to assert on if you are running a batch through a server you do not control.
Without it, a client whose gene model was replaced mid-batch got a well-formed `200` computed against the wrong model, with nothing in the response to say so: a BRCA1 coding variant coming back `intergenic_variant`, plausibly.
`/api/status` showed the swap, but only to a client that polled it between every request.

`generation` is the field to compare, not the other two.
It counts replacements, so it changes on every one; `gff3_source` reads `user-upload` for every upload, and a transcript count can coincide across two different models.

Multiple records in one request are annotated in one pass, and an `ID` column is carried through to `id` in the response:

```json
{
  "count": 2,
  "time_ms": 0,
  "results": [
    {
      "seq_region_name": "17",
      "start": 43124027,
      "end": 43124027,
      "allele_string": "G/A",
      "id": null,
      "strand": 1,
      "most_severe_consequence": "synonymous_variant",
      "transcript_consequences": [
        {
          "transcript_id": "ENST00000357654",
          "gene_id": "ENSG00000012048",
          "gene_symbol": "BRCA1",
          "biotype": "protein_coding",
          "consequence_terms": ["synonymous_variant"],
          "impact": "LOW",
          "hgvsc": "ENST00000357654.9:c.70C>T",
          "hgvsp": "ENSP00000350283.1:p.Cys24=",
          "hgvsg": "17:g.43124027G>A",
          "strand": -1,
          "variant_allele": "A",
          "canonical": 1
        }
      ]
    },
    {
      "seq_region_name": "17",
      "start": 43104121,
      "end": 43104121,
      "allele_string": "G/T",
      "id": "rs80357064",
      "strand": 1,
      "most_severe_consequence": "splice_donor_variant",
      "transcript_consequences": [
        {
          "transcript_id": "ENST00000357654",
          "gene_id": "ENSG00000012048",
          "gene_symbol": "BRCA1",
          "biotype": "protein_coding",
          "consequence_terms": ["splice_donor_variant"],
          "impact": "HIGH",
          "hgvsc": "ENST00000357654.9:c.441+1C>A",
          "hgvsp": null,
          "hgvsg": "17:g.43104121G>T",
          "strand": -1,
          "variant_allele": "T",
          "canonical": 1
        }
      ]
    }
  ]
}
```

Both `transcript_consequences` arrays are elided above: against a full human gene model those two variants return 58 and 41 entries respectively, one per overlapping transcript.
Expect responses to be large, and select the transcript you care about rather than assuming the first entry is the interesting one.

The per-variant shape follows Ensembl VEP's JSON output, so clients written against VEP's `transcript_consequences` mostly port over.
`time_ms` is server-side annotation time, excluding transport.

`"pick": true` reduces each variant to one transcript, chosen by the same hierarchy as `fastvep annotate --pick`: MANE Select, then MANE Plus Clinical, canonical, APPRIS, TSL, biotype, CCDS, and only then consequence severity.
That is Ensembl VEP's default `--pick_order`, so a default run of either tool picks the same transcript.
It reduces to one *transcript*, not to one entry: `transcript_consequences` carries one entry per (transcript, allele), so a biallelic site still returns two.
A variant that overlaps no transcript at all is left untouched, because there its entries are one per alt allele rather than competing transcripts.

`"pick_allele": true` is the one that reduces to one *entry* per allele, which is what Ensembl VEP's `--pick` does; `"pick_allele_gene": true` keeps one per allele per gene.
The three `flag_pick*` switches return everything and mark the winner with `"pick": 1` on the chosen `transcript_consequences` entries, so a client can apply its own rule and still see what the hierarchy chose.
The key is absent, not `0`, on the entries that were not picked and on every entry of a request that picked nothing.

Which of these agrees with VEP, measured, is in [docs/VEP_DIVERGENCE.md](VEP_DIVERGENCE.md#10---pick-keeps-every-allele-of-the-transcript-it-picks): the four allele-scoped switches agree exactly, and `pick` / `flag_pick` keep every allele of the winning transcript where VEP would drop all but one.

> **`pick` does not mean "most severe".**
> Severity is the *last* tie-break in that order, so where genes overlap, a neighbouring gene's MANE transcript outranks a non-MANE transcript the variant actually disrupts, and the reported consequence can be `upstream_gene_variant` on the neighbour rather than the damaging term on the real target.
> Stock VEP behaves the same way; it is still the wrong answer for clinical reporting.
> The server does not expose `--pick-order`, so if severity must come first, use the CLI with `--pick-order rank,mane_select,...` -- [docs/ACMG.md](ACMG.md#which-transcript---pick-reports) has the measured effect.
> `most_severe_consequence` is computed from the transcripts actually reported, so with `"pick": true` it describes the picked transcript, not the whole locus.

Setting `"acmg": true` adds an `acmg` object to each transcript consequence with the classification, the per-criterion verdicts, and the evidence counts.
Classification quality depends on which supplementary databases are loaded, so read [docs/ACMG_SETUP.md](ACMG_SETUP.md) before relying on it.

### `GET /api/genomes` and `POST /api/load-genome`

Start the server with `--data-dir` pointing at a directory of genome subdirectories (the layout is described under [Multi-organism setup](../README.md#multi-organism-setup) in the README) and the server can switch between them at runtime.

```bash
curl -s http://localhost:8080/api/genomes
```

```json
{"genomes": [{"name": "human_grch38", "has_fasta": true, "has_sa": true}]}
```

```bash
curl -s -X POST http://localhost:8080/api/load-genome \
  -H 'Content-Type: application/json' \
  -d '{"name": "human_grch38"}'
```

```json
{"name": "human_grch38", "transcripts": 509650, "sa_sources": ["ClinVar"], "generation": 1, "time_ms": 6}
```

`name` must be a plain directory name; anything containing a path separator or `..` is rejected with 400.

Loading a genome swaps the model out from under every client of that server, so it is an admin operation rather than something a request handler calls, and the server now treats it as one.
`POST /api/load-genome` returns `403` unless the server was started with `--allow-model-replacement`; `GET /api/genomes` is a read and always answers.
The `generation` in the response is the value `/api/annotate` will report for answers computed against this model.

## Errors

Failures return the matching HTTP status and a single-field body:

```json
{"error": "No VCF data provided"}
```

`400` carries a caller-facing message, as above.
`403` means the request was understood and refused by how the server was started, so changing the request will not help. Currently only the two gene-model-replacing endpoints return it, and the message names the flag that lifts it.
`500` is always the literal string `Internal server error`; the detail is written to the server log instead, so that paths and parser internals are not returned to clients.
Check the server's stderr when you get one.

## Options

| Flag | Env | Default | Meaning |
| --- | --- | --- | --- |
| `--bind` | `FASTVEP_BIND` | `0.0.0.0` | Interface to listen on |
| `--port` | `FASTVEP_PORT` | `8080` | Port to listen on |
| `--gff3` | `FASTVEP_GFF3` | none | Gene model |
| `--fasta` | `FASTVEP_FASTA` | none | Reference FASTA, needs a `.fai` |
| `--sa-dir` | `FASTVEP_SA_DIR` | none | Directory of `.osa`/`.osa2`/`.osi`/`.oga` files |
| `--data-dir` | `FASTVEP_DATA_DIR` | none | Directory of switchable genomes |
| `--distance` | | `5000` | Upstream/downstream distance in bp |
| `--max-body-size` | | `10485760` | Request body cap in bytes |
| `--max-concurrent` | | `64` | Concurrent annotation requests |
| `--stats-file` | `FASTVEP_STATS_FILE` | `stats.json` | Where the served-variant counters persist |
| `--allow-model-replacement` | `FASTVEP_ALLOW_MODEL_REPLACEMENT` | off | Let `/api/upload-gff3` and `/api/load-genome` replace the gene model |

Three things to know about the defaults.
`--gff3` has none, and a server started without one does not refuse to start: it comes up with `"status": "ok"`, reports `"transcripts": 0`, and answers every request with `intergenic_variant` at HTTP 200.
`--stats-file` is relative, so a bare `fastvep-web` writes `stats.json` into whatever directory you launched it from; pass an absolute path under a service account's own state directory when running it as a daemon.
Requests beyond `--max-concurrent` queue rather than fail, which is why raising it is rarely the fix for a slow server.

## Which binary

Use `fastvep-web`, the axum server documented here.

The older `fastvep web` subcommand still works, but it is single-threaded, and it serves only `/api/status`, `/api/annotate`, and `/api/upload-gff3` (no `/api/genomes` or `/api/load-genome`).
