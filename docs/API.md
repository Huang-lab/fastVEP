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

To serve several genomes from one process, give it a directory with one subdirectory per genome.
A request then names the genome it wants, as described under [Choosing the genome](#choosing-the-genome):

```bash
fastvep-web \
  --data-dir genomes/ \
  --gff3 genomes/hg38_ensembl_115/Homo_sapiens.GRCh38.115.gff3 \
  --fasta genomes/hg38_ensembl_115/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --sa-dir genomes/hg38_ensembl_115/sa \
  --preload-genomes
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
Turn it on for a single-user desktop run that wants to upload its own GFF3 or change the default gene model.
Leave it off for anything shared, and the gene model is a property of how the server was started.

Serving more than one genome does not need it.
With `--data-dir`, a request names the genome it wants and is answered from that genome's own model.
Nothing is replaced, so no client can change what another is answered against.
See [Choosing the genome](#choosing-the-genome).

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
    "default_genome": "hg38_ensembl_115",
    "gene_model_generation": 0,
    "gff3_source": "Homo_sapiens.GRCh38.115.gff3",
    "has_fasta": true,
    "limits": {"max_body_bytes": 10485760, "max_response_bytes": 134217728, "max_request_seconds": 90, "max_parallel": 4, "max_alt_alleles": 256},
    "sa_sources": ["ClinVar", "dbSNP"],
    "sa_info": [
        {"key": "clinvar", "name": "ClinVar", "version": "2026-11-02", "assembly": "GRCh38", "files": 1},
        {"key": "dbsnp", "name": "dbSNP", "version": "latest", "assembly": "GRCh38", "files": 1}
    ],
    "sa_skipped": [],
    "acmg": {
        "ready": false,
        "available": false,
        "loaded": ["clinvar"],
        "missing": [
            {"id": "gnomad", "label": "gnomAD allele frequencies", "scope": "allele", "drives": ["BA1", "BS1", "BS2", "PM2"], "without": "no population frequency: ..."}
        ]
    },
    "status": "ok",
    "total_genomes": 0,
    "total_variants": 0,
    "transcripts": 509650,
    "version": "0.5.0"
}
```

Check `transcripts`, `has_fasta`, and `sa_sources` here before you trust any annotation.
Each of the three can be missing without the server saying so: a zero transcript count makes every variant `intergenic_variant`, and the other two are the subject of the next section.

`sa_info` says what each loaded source is: its `key`, its `name`, its `version`, which is the date of the release when the release dated itself (a VCF with a `##fileDate=` line, as ClinVar's has) and `latest` when it did not, its `assembly`, and how many files hold it (`files`; a source kept in a shard for each chromosome is one source).
`sa_skipped` lists the database files that were found in the directory and not loaded, by name and by class of reason, `empty` (a build that found nothing to write) or `unopenable` (a truncated copy, a file of another format), so that a source that is missing can be told from one that was never provided; the reason in words, which can carry a path, is in the server's log.
Each entry of `GET /api/genomes` carries `sa_skipped` for its own directory too.

`limits` is what one request may cost, so that a client can size its requests instead of finding out with a `413`; see [What one request may cost](#what-one-request-may-cost).

`acmg` is the server's verdict on whether `"acmg": true` will be answered for the default genome, and which of the nine sources the classifier reads are missing if it will not (abridged above: `missing` has one entry for each source that is not loaded); see [ACMG needs the reference stack](#acmg-needs-the-reference-stack).

`allow_model_replacement` says whether a request can change the gene model on this server; `gene_model_generation` counts how many times one has.
A UI can read the first to decide whether to offer its model-replacing controls at all, rather than finding out with a `403` a user cannot interpret.

`default_genome` is the name of the genome that a request naming none is answered from, and the model fields above describe that genome.
A server that sends it can also answer from any other genome of its `--data-dir` when a request names one.
A server that does not send it, such as the older `fastvep web`, answers only from the gene model it started with.

`total_genomes` counts the analysis runs that named a genome, plus the `/api/load-genome` calls of a server that allows them.

## Always pass `--fasta`

`--fasta` is optional to the argument parser and effectively required for correct answers.
Without reference sequence the server cannot read the codon a variant falls in, so a coding change degrades to a generic term instead of failing.
Annotating BRCA1 `17:43124025 A>G` against the same GRCh38 gene model, with and without `--fasta` (measured on 8 November 2026):

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
| `GET` | `/api/genomes` | List the genomes a request can name, and which are loaded |
| `POST` | `/api/load-genome` | Replace the default genome with another one - needs `--allow-model-replacement` |
| `POST` | `/api/upload-gff3` | Replace the active gene model with a posted GFF3 - needs `--allow-model-replacement` |

The two marked endpoints replace the gene model for every client of the server, so they return `403` unless it was started with `--allow-model-replacement`.

`GET /` serves the browser interface.

### `POST /api/annotate`

| Field | Type | Default | Meaning |
| --- | --- | --- | --- |
| `vcf` | string | required | VCF text. Header lines are optional; records are tab-separated. |
| `genome` | string | the default genome | The genome to answer from, by a name that `GET /api/genomes` lists. An unknown name is refused with `400`; it is never answered from the default instead. |
| `pick` | bool | `false` | Reduce each variant to one transcript consequence (see below). |
| `pick_allele` | bool | `false` | Reduce to one consequence per (variant, allele). VEP's `--pick_allele`. |
| `pick_allele_gene` | bool | `false` | Reduce to one consequence per (variant, allele, gene). VEP's `--pick_allele_gene`. |
| `flag_pick` | bool | `false` | Return every consequence, with `"pick": 1` on the one `pick` would have kept. VEP's `--flag_pick`. |
| `flag_pick_allele` | bool | `false` | As `flag_pick`, marking the pick per allele. |
| `flag_pick_allele_gene` | bool | `false` | As `flag_pick`, marking the pick per allele and gene. |
| `acmg` | bool | `false` | Attach ACMG-AMP classification to each consequence. Refused with `422` for a genome whose database directory lacks part of the reference stack; see [ACMG needs the reference stack](#acmg-needs-the-reference-stack). |

At most one of the six pick switches may be set; the request is read in the order listed.

```bash
curl -s -X POST http://localhost:8080/api/annotate \
  -H 'Content-Type: application/json' \
  -d '{"vcf": "17\t43124025\t.\tA\tG\t50\tPASS\t.", "acmg": false}'
```

Every response also reports the gene model that produced it:

```json
"gene_model": {"genome": "hg38_ensembl_115", "gff3_source": "Homo_sapiens.GRCh38.115.gff3", "transcripts": 509650, "generation": 0}
```

This is the field to assert on if you are running a batch through a server you do not control.
Without it, a client whose gene model was replaced mid-batch got a well-formed `200` computed against the wrong model, with nothing in the response to say so: a BRCA1 coding variant coming back `intergenic_variant`, plausibly.
`/api/status` showed the swap, but only to a client that polled it between every request.

`genome` names the genome that answered, in the form a request uses to ask for one.
A client that chose a genome should check that it got the one it asked for.

`generation` is the field to compare for the default genome, not the other two.
It counts replacements, so it changes on every one; `gff3_source` reads `user-upload` for every upload, and a transcript count can coincide across two different models.
A genome that a request chose is never replaced, so its `generation` is always `0`.

Multiple records in one request are annotated in one pass, and an `ID` column is carried through to `id` in the response:

```json
{
  "count": 2,
  "time_ms": 0,
  "results": [
    {
      "seq_region_name": "17",
      "start": 43124025,
      "end": 43124025,
      "allele_string": "A/G",
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
          "hgvsc": "ENST00000357654.9:c.72T>C",
          "hgvsp": "ENSP00000350283.3:p.Cys24=",
          "hgvsg": "17:g.43124025A>G",
          "strand": -1,
          "variant_allele": "G",
          "canonical": 1
        }
      ]
    },
    {
      "seq_region_name": "17",
      "start": 43104121,
      "end": 43104121,
      "allele_string": "C/A",
      "id": "rs397509172",
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
          "hgvsc": "ENST00000357654.9:c.441+1G>T",
          "hgvsp": null,
          "hgvsg": "17:g.43104121C>A",
          "strand": -1,
          "variant_allele": "A",
          "canonical": 1
        }
      ]
    }
  ]
}
```

Both `transcript_consequences` arrays are elided above: against the full GRCh38 gene model those two variants return 58 and 41 entries respectively, one per overlapping transcript, and each record also carries the ClinVar and dbSNP fields of a transcript entry when those databases are loaded and hold the variant.
The REF alleles are the genome's: `A` at 17:43124025 and `C` at 17:43104121.
fastVEP writes HGVS against the reference base, as Ensembl does, so a REF that is not the genome's changes the HGVS you get back; the CLI counts such records and says so at the end of a run.
Expect responses to be large, and select the transcript you care about rather than assuming the first entry is the interesting one.

The per-variant shape follows Ensembl VEP's JSON output, so clients written against VEP's `transcript_consequences` mostly port over.
`time_ms` is server-side annotation time, excluding transport.

`warnings`, present only when there is something to say, lists the databases that could not be read for some of the lookups while the request was answered (`gnomad could not be read for 3 variant lookup(s) ...`).
A database that opens and cannot inflate one chunk answers "no record" at the positions of that chunk, which is the same in the annotation as a variant it does not hold, and for gnomAD that reads as absence from the population; the warning is how a client tells the two apart.
The count is the server's, so a lookup of a concurrent request that failed is among them, and `GET /api/status` carries the total since the server started as `sa_lookup_errors` (`{"gnomad": 3}`, empty when every lookup was answered).

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
What it can say depends on which supplementary databases are loaded, which the next section is about.

### ACMG needs the reference stack

The classifier reads nine supplementary sources: gnomAD, ClinVar, REVEL, SpliceAI, PhyloP, gnomAD gene constraint, ClinGen gene-disease validity, ClinVar's protein and splice index, and a repeat track ([docs/ACMG_SETUP.md](ACMG_SETUP.md)).
A criterion with none of its source reports `evaluated: false` and nothing else.
A server that holds only some of the nine therefore answers `"acmg": true` with a classification that rests on less evidence than it appears to, almost always VUS, and nothing in a VUS says that the evidence was never looked at.
So the server decides before it classifies:

- A genome whose database directory holds all nine is *ready*, and `"acmg": true` is answered, provided its gnomAD databases hold every autosome and X that the gene model has transcripts on.
  A gnomAD that holds only some of them (one shard of a per-chromosome build, a chromosome that failed to build) says "no record" for the rest, which reads as absence from the population, so it counts as missing: its entry in `missing` is labelled `none for chromosomes 1-22` (the ones it lacks) and says that PM2 would fire there.
  The response carries `"acmg_stack": {"complete": true, "missing": []}` beside the results.
- Any other genome is refused with HTTP `422` before anything is annotated, and the body says what is missing:

  ```json
  {
    "error": "ACMG classification is not available for genome 'hg38_ensembl_115': 8 of the 9 sources it reads are not loaded, and a classification made without them would rest on less evidence than it appears to. Missing: gnomAD allele frequencies (BA1, BS1, BS2, PM2); ...",
    "acmg": {"ready": false, "available": false, "loaded": ["clinvar"], "missing": [{"id": "gnomad", "label": "gnomAD allele frequencies", "scope": "allele", "drives": ["BA1", "BS1", "BS2", "PM2"], "without": "..."}]}
  }
  ```

- `--acmg-allow-partial` (`FASTVEP_ACMG_ALLOW_PARTIAL`) answers anyway, for a run where an incomplete stack is deliberate.
  Every response then carries `"acmg_stack": {"complete": false, "missing": ["gnomad", ...]}`, so that a client that stores or shows the classifications keeps what they rest on with them.
  A criterion that reads a missing source is `evaluated: false`; PM2 in particular is not evaluated when no gnomAD database is loaded, where it used to be met for every variant.

`GET /api/status`, for the default genome, and each entry of `GET /api/genomes` carry the same verdict in `acmg`, so that a client can offer classification only where it will be answered and say why where it will not.
`ready` means all nine are loaded; `available` means a request for ACMG will be answered, and differs from `ready` only under `--acmg-allow-partial`.
`loaded` and `missing` use the `json_key` of each source as its `id`, and each missing source says which criteria read it (`drives`) and what is lost without it (`without`).
For a genome that is not loaded yet the verdict comes from what its `sa/` directory holds, found once by opening its files, so that a client does not have to load a genome to ask; restart the server after adding databases to a directory.
The chromosomes of a gene model are known once it is loaded, so the check of what gnomAD holds is added then: the server preloads its genomes at start, and a client that asks for a genome before it has loaded is told what its files hold.

### What one request may cost

A response is built in memory, and on a 16 GB host that is the one thing a single client can use to take the server down.
Measured on 8 November 2026 against GRCh38, with the reference stack of fastvep.org, 1,000 ClinVar variants came to about 43 KB a variant, 18 KB with `pick`, and 27 KB with `pick` and `acmg`, and variants spread over the genome to about 18 KB with `pick`, so that the 128 MiB limit holds about 7,300 of them.
Without `pick`, an earlier measurement of `acmg` was 258 KB a variant, and 10,000 variants with ACMG came to several gigabytes in the process.
The 10 MiB the server accepts as a body holds some 440,000 variants.
So a request is bounded, and told how when it is not met:

- **Malformed input is a `400`.**
  Every record is read before any is annotated, and a line that cannot be read is refused with its number and the parser's reason (`VCF line 2 cannot be read: Invalid POS: foo`).
  It used to be an internal error, met at the line that failed, with nothing said.
  A record with more than 256 alternate alleles (`limits.max_alt_alleles`) is refused before it is parsed, because the cost of annotating one grows with the square of their number (200,000 alleles in a 2.4 MB line ran 108 s), and a position or END beyond 4,294,967,295 is refused because the arithmetic of the annotation overflows there.
- **The answer has a size limit.**
  The annotation is made a batch of variants at a time and written out as it goes.
  A batch is at most 100 variants and, past the first, at most 4,000 rows, a row being a transcript that a variant overlaps (for each of its alternate alleles): a variant makes a row for every transcript, so an SNV is a few rows and a deletion of a chromosome arm thousands, and a hundred of those held 2.4 GB before the first byte went out.
  A variant whose answer alone would pass the limit (it overlaps more than about 190,000 transcripts, or 26,000 with `acmg`, which are 700 and 5,000 bytes a row at least) is refused before it is annotated, with its line and the transcripts it overlaps, and the same request with `pick` is answered.
  A request whose answer passes `--max-response-bytes` (128 MiB by default) is refused with `413` at the batch that passes it, saying how many variants were done and what to change: fewer variants per request, `pick`, or the command-line tool.
  The batches are also what bound the memory, since the answer as values is five to ten times the size of its text.
- **The work has a time limit.**
  A request that runs past `--max-request-seconds` (90 by default) is refused with `413`, with what was done.
  The time counts from the start of the work, once the request has a place and its genome is loaded, and not from its arrival: a request that waited 700 ms for a place with 500 ms allowed was refused before it began.
  Keep it below the read timeout of the proxy in front of the server.
- **Annotations run in at most `--max-parallel` places** (one for each core by default), because one annotation is one thread and more of them than cores only make each slower and the process larger.
  A request that finds every place taken waits `--queue-timeout-seconds` (30) for one, and is then turned away with `503` and `Retry-After: 5`.
  Reads (`/api/status`, `/api/genomes`) do not need a place.
- **A client that goes away stops its annotation** between batches, instead of finishing an answer nobody will read, and its place comes back.

A client with a large file should send it in batches and halve the batch on a `413`: a thousand variants with `pick` is about 14 MB, which is well inside the default.

## Choosing the genome

Start the server with `--data-dir` pointing at a directory of genome subdirectories (the layout is described under [Multi-organism setup](../README.md#multi-organism-setup) in the README), and a request can name the genome it wants.
`GET /api/genomes` lists what it can name:

```bash
curl -s http://localhost:8080/api/genomes
```

```json
{
  "default": "hg38_ensembl_115",
  "genomes": [
    {"name": "hg38_ensembl_115", "gff3": "Homo_sapiens.GRCh38.115.gff3", "has_fasta": true, "has_sa": true, "loaded": true, "default": true, "transcripts": 509650},
    {"name": "mouse_grcm39_ensembl_115", "gff3": "Mus_musculus.GRCm39.115.gff3", "has_fasta": true, "has_sa": false, "loaded": false, "default": false, "transcripts": null}
  ]
}
```

`name` is what a request passes as `genome`.
`gff3` is the file the genome is read from, which Ensembl names `<Species>.<Assembly>.<release>.gff3`, so a client can say what a genome is without the server describing it.
`loaded` says whether a request for it is answered at once or waits for it to load, and `transcripts` is what it holds once it has.
The genome the server started with is listed even when it is not under `--data-dir`, and is always `loaded`.

```bash
curl -s -X POST http://localhost:8080/api/annotate \
  -H 'Content-Type: application/json' \
  -d '{"vcf": "1\t4878171\t.\tC\tT\t.\t.\t.", "genome": "mouse_grcm39_ensembl_115"}'
```

The response says which genome answered in `gene_model.genome`, as described under [`POST /api/annotate`](#post-apiannotate).
Choosing a genome replaces nothing and needs no `--allow-model-replacement`: the answer comes from that genome's own model, and no other client's answers move.

How a chosen genome is served:

- A genome is loaded the first time something asks for it and then kept, so memory grows with the genomes in use, bounded by the number of directories under `--data-dir`.
  Concurrent first requests for one genome share a single load.
  `--preload-genomes` loads all of them in the background at startup instead, so that no request waits.
  On the hosted site's nine genomes that took about 70 seconds the first time, while the transcript caches were being built, and the genome the server started with answered from the first second.
- Each genome is read from its own GFF3, FASTA and `sa/` directory and from nothing else.
  The `--sa-dir` the server was started with belongs to the default genome and is never lent to another: on a server that holds human ClinVar and dbSNP there, a mouse variant that shares a coordinate with a human record would otherwise be annotated with the human record.
- A name that is not a genome is refused with `400`, and is never answered from the default genome instead.
  A misspelling that quietly returned human consequences for a mouse variant would look exactly like a correct answer.
  Names are matched against what is on disk letter for letter, on filesystems that ignore case as well.
- Every answer is the same as the one the server gives when it is started on that genome alone, and the options of `/api/annotate` (`pick`, `acmg` and the rest) apply to a chosen genome as they do to the default.
  ACMG-AMP classification is defined for human variants, and a genome with no supplementary annotations holds none of the sources it reads, so the server refuses it there ([ACMG needs the reference stack](#acmg-needs-the-reference-stack)).
- A genome is loaded by a task of its own, and a request only waits for the outcome, so a client that gives up while a genome loads neither stops the load nor starts a second one.
  A load that fails is not remembered: the next request tries again.

### Replacing the default genome

`POST /api/load-genome` is the other way to change what a request that names no genome is answered from: it replaces the default gene model, for every client, with another genome of `--data-dir`.

```bash
curl -s -X POST http://localhost:8080/api/load-genome \
  -H 'Content-Type: application/json' \
  -d '{"name": "human_grch38"}'
```

```json
{"name": "human_grch38", "transcripts": 509650, "sa_sources": ["ClinVar"], "generation": 1, "time_ms": 6}
```

`name` must be a plain directory name; anything containing a path separator or `..` is rejected with 400.

Loading a genome this way swaps the model out from under every client of that server, so it is an admin operation rather than something a request handler calls, and the server treats it as one.
`POST /api/load-genome` returns `403` unless the server was started with `--allow-model-replacement`; `GET /api/genomes` is a read and always answers.
The `generation` in the response is the value `/api/annotate` will report for answers computed against this model, and the default genome is renamed to the one just loaded, so a later request that names the genome it replaced is answered from that genome's own model and not from the new default.
Unlike a genome chosen per request, this keeps the older behaviour of falling back to `--sa-dir` for a genome that has no `sa/` of its own, which a single-user run with one SA directory and several gene models relies on.

## Errors

Failures return the matching HTTP status and a single-field body:

```json
{"error": "No VCF data provided"}
```

`400` carries a caller-facing message, as above; an unknown or empty `genome` is one, and so is a VCF line that cannot be read.
`403` means the request was understood and refused by how the server was started, so changing the request will not help. Currently only the two gene-model-replacing endpoints return it, and the message names the flag that lifts it.
`413` means the request is larger than this server answers in one go: its answer would pass `--max-response-bytes`, or it ran past `--max-request-seconds`. The message says how many variants were done and what to change.
`503` means every place for an annotation was taken and none came free in time; it carries `Retry-After`, and the request is worth repeating.
`422` means the request is well formed and the server could answer it, but not faithfully: ACMG classification was asked of a genome that lacks sources the classifier reads.
The body carries the verdict under `acmg` beside `error`, and the request is not counted as a run.
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
| `--data-dir` | `FASTVEP_DATA_DIR` | none | Directory of genomes, one subdirectory each, that a request can name |
| `--preload-genomes` | `FASTVEP_PRELOAD_GENOMES` | off | Load every genome of `--data-dir` in the background at startup instead of on its first request |
| `--distance` | | `5000` | Upstream/downstream distance in bp |
| `--max-body-size` | | `10485760` | Request body cap in bytes |
| `--max-concurrent` | | `64` | Concurrent requests to each route (every route has its own limit); beyond it a request waits |
| `--max-parallel` | `FASTVEP_MAX_PARALLEL` | cores | Annotations that run at once; a request finds a place or waits `--queue-timeout-seconds` and is turned away with `503` |
| `--queue-timeout-seconds` | `FASTVEP_QUEUE_TIMEOUT_SECONDS` | `30` | How long a request waits for a place |
| `--max-response-bytes` | `FASTVEP_MAX_RESPONSE_BYTES` | `134217728` | Largest JSON answer one request may produce; beyond it `413` |
| `--max-request-seconds` | `FASTVEP_MAX_REQUEST_SECONDS` | `90` | Longest one annotation may run; beyond it `413` |
| `--stats-file` | `FASTVEP_STATS_FILE` | `stats.json` | Where the served-variant counters persist |
| `--allow-model-replacement` | `FASTVEP_ALLOW_MODEL_REPLACEMENT` | off | Let `/api/upload-gff3` and `/api/load-genome` replace the gene model; choosing a genome per request does not need it |
| `--acmg-allow-partial` | `FASTVEP_ACMG_ALLOW_PARTIAL` | off | Answer `"acmg": true` for a genome whose database directory lacks part of the reference stack, saying so in `acmg_stack`, instead of refusing with `422` |

The three switches take `true`, `1`, `yes` or `on` from their environment variable.

Three things to know about the defaults.
`--gff3` has none, and a server started without one does not refuse to start: it comes up with `"status": "ok"`, reports `"transcripts": 0`, and answers every request with `intergenic_variant` at HTTP 200.
`--stats-file` is relative, so a bare `fastvep-web` writes `stats.json` into whatever directory you launched it from; pass an absolute path under a service account's own state directory when running it as a daemon.
Requests beyond `--max-concurrent` queue rather than fail, which is why raising it is rarely the fix for a slow server; `--max-parallel` is the number of annotations that run at once, and a request that cannot get a place in `--queue-timeout-seconds` is turned away with `503`.

## Which binary

Use `fastvep-web`, the axum server documented here.

The older `fastvep web` subcommand still works, but it is single-threaded, and it serves only `/api/status`, `/api/annotate`, and `/api/upload-gff3` (no `/api/genomes` or `/api/load-genome`).
