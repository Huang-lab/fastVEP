# Working in this repository

fastVEP predicts the functional consequence of genomic variants. The output is
a prediction a clinician or researcher may act on, which sets the bar for
everything below: **the failure mode that matters here is not a crash, it is a
plausible-looking wrong answer.** Most of the fixes in `git log` are exactly
that — a variant annotated against the wrong transcript, a gene that lost its
symbol, a cache trusted after the parser that wrote it was fixed. None of them
would have shown up as an error message.

## Layout

Twelve crates, layered so that nothing below depends on anything above it:

```
fastvep-core            Consequence (49 SO terms), Allele, VariantType, Impact
  fastvep-genome        Transcript, Exon, Gene, CodonTable
    fastvep-cache       GFF3 parser, FASTA reader, transcript cache, providers
    fastvep-consequence consequence prediction (small variants + SVs)
    fastvep-hgvs        HGVS g. / c. / p. nomenclature
    fastvep-io          VCF parsing, CSQ / tab / JSON output
  fastvep-filter        filter_vep-compatible expression engine
  fastvep-classification ACMG-AMP criteria (Richards 2015 + ClinGen SVI)
      fastvep-sa        supplementary annotations (.osa/.osa2/.osi/.oga)
        fastvep-annotate  the shared annotation pipeline
          fastvep-cli     the `fastvep` binary
          fastvep-web     the axum production server
```

`fastvep-annotate` and `fastvep-cli` both drive annotation and their per-variant
loops have drifted apart before. A change to one is usually a change to both —
check the other before you call the work done.

## The loop

```bash
cargo fmt --all                                        # before every commit
cargo clippy --workspace --all-targets -- -D warnings  # what CI gates on
cargo test --workspace                                 # 1,731 tests, all must pass
```

CI runs exactly these three and nothing else, in this order. The lint set lives
in `[workspace.lints]` in the root `Cargo.toml`, so a local `cargo clippy`
already sees what CI sees — the only difference is the `-D`.

Keep that set narrow. Every lint in it is one the workspace currently
satisfies, which is what makes a warning worth reading. Adding a lint that
fires on hundreds of existing lines does not raise the bar; it trains everyone
to scroll past the output.

## Conventions the history already follows

**Commits.** Conventional commits with a scope, and a subject that names the
consequence rather than the edit:

```
fix(gff): read ncRNA_gene records, so non-coding genes keep their symbol (#98)
fix(cache): reject the pre-#90 cache format instead of trusting it (#95)
```

`cliff.toml` drafts release notes from these, but `CHANGELOG.md` is
hand-curated: git-cliff writes the unreleased block, a human writes the prose.
Released sections are never regenerated. Bug-fix commits do not touch
`CHANGELOG.md`; it is updated when a release is cut.

**Silencing a lint.** `#[allow(...)]` carries a one-line reason at the site.
An `#[allow]` with no explanation is a TODO nobody wrote down.

**Comments explain why, not what.** The code says what it does. Comments in
this repo exist to record the reasoning that is not recoverable from reading
it — which provider's file format forced a branch, which issue proved an
assumption wrong, why the obvious simplification is unsafe. When you fix a bug,
leave behind the fact that made it a bug.

**Formatting-only commits** go in alone and get their hash added to
`.git-blame-ignore-revs`. Never mix a reformat with a behaviour change; it
makes the behaviour change invisible in review.

## Two hazards worth knowing before you start

**A parser fix does not reach users on its own.** `fastvep annotate` writes a
sidecar transcript cache next to the GFF3 and reuses it whenever its mtime is
newer than the source. A user who fixes their fastVEP but not their GFF3 keeps
reading the cache the *old* parser wrote, and sees no change at all. So a fix
to `crates/fastvep-cache/src/gff.rs` that changes what gets parsed needs the
cache magic in `transcript_cache.rs` bumped alongside it, with the old magic
moved to the rejected list and a `SUMMARY_PRE_NN` / `DETAIL_PRE_NN` pair saying
what was wrong with what it holds. `#95` and `#98` are the worked examples.

The sidecar is also a sequence of zstd frames of about 8 MB with a table of them at its end (a
skippable frame, in the layout of zstd's seekable format), which a run with more than two threads
reads to decode the frames side by side. The bytes of the stream the frames make are what the one
frame of earlier builds held, and a decoder passes over the table, so a reader that takes the file
as one stream, as every build before #166 does, still reads it. A change to how a cache is cut has
to keep that, or bump the magic. A cache has to be read to its end, where zstd checks the checksum
of each frame: a loader that stopped at the last transcript let a flipped byte, or a second cache
joined on behind the first, load without a word. And nothing is allocated on the word of what is
read before that checksum has been checked: every `bincode` call on a cache has a limit
(`bounded`), because a damaged byte that decodes to a length of a terabyte is otherwise an abort,
with no word on which cache it was and no rebuild from the GFF3; so `save_cache` refuses a
transcript that the limit would refuse at every load. The table is read only by a run that uses
it (more than two threads), so a cache whose table is damaged loads on two threads and is refused
on three: the frames have their own checksums, and the table is only where they are.

The pair is two lengths of the same fact, because the two callers are not in
the same situation. The sidecar path rebuilds from the GFF3 by itself, so it
prints the summary and stops; the explicit `--transcript-cache` path has to
refuse the run, so it owes the user the detail. Recovery advice belongs at the
call site, not in the error: the sidecar has already rebuilt by the time it
reports, and deleting a cache named on the command line produces no sidecar to
replace it.

Rejecting is deliberate, and the reason generalises: a stale cache is
*intact*, so nothing in it distinguishes a gene that lost its symbol from one
that never had a symbol. There is no safe way to read it, only a safe way to
refuse it.

**Region-restricted transcript sets must never be cached.** Reading a GFF3
through the tabix path returns only the features overlapping *this* VCF's
variants. That set is correct for this run and wrong for every other one, so
`run_annotate` tracks `region_restricted` and skips the cache save. If you add
another partial-load path, it has to set that flag.

## Performance

The release profile is already tuned (`lto = "fat"`, `codegen-units = 1`).
What remains is allocation. A callgrind profile of the annotate path over a
gene-dense locus, with the system allocator, put about a quarter of its
instructions in `malloc`/`free`, and output formatting, not consequence
prediction, accounted for most of the rest. When you touch the per-variant
loop, the useful question is not "is this branch fast" but "does this allocate
per variant, per transcript, or per allele": the loop runs once per (variant ×
overlapping transcript × allele), so a `String` built there is built millions of
times.

Three shapes account for most of what has been found there so far. One is work
done before deciding whether it is needed: a value built on entry as a fallback
the common path discards. The other is a per-run cost paid per item: resolving
the same column names, or rescanning the same list, once for every record
rather than once. The third is a choice made after the work it was meant to
save: the pick (`--pick*`) used to run after a row had been built for every
transcript of the window, which at `--distance 500000` is hundreds a variant,
and then deleted all but one per gene. It is made now on what the predictor
said, before any row exists (`crates/fastvep-annotate/src/pick.rs`), and
`crates/fastvep-cli/tests/pick_allocations.rs` holds that the allocations of a
record follow the transcripts that are kept and not the ones that are around.

Threads. Free what a thread allocated on that thread: one thread freeing a whole
batch that the others had built was 11 to 13% of the batch loop at 8 threads, so
a record is annotated, rendered and freed by the thread that takes it
(`render_batch`), and each batch is written on a thread of its own. Peak memory
should not grow with the thread count either: look at 1, 8 and the cores of the
machine (the parallel decode of the transcript cache is capped at 8 decoders
because its peak did).

The `fastvep` binary uses mimalloc (`crates/fastvep-cli/src/main.rs`). The
libraries and their tests, and the `fastvep-web` binary, keep the system
allocator; `fastvep web`, a subcommand of the same binary, does not. It does
less work per allocation, and it avoids something glibc does to a thread that
builds a large structure: glibc grows the arena of each thread a page at a time,
one `mprotect` for every 4 KB, so decoding the transcript cache on 8 threads
made 231,689 of them where mimalloc made 62, and under Rosetta 2 that left the
run no faster than master's one thread. The price is peak memory on Linux, 4 to
9% (100 to 230 MB) over glibc on one to 32 threads and 12% on 64, and about 1.2 GiB more address space under
`ulimit -v`. The counting-allocator tests wrap the system allocator, since the
library crates they test never link mimalloc: they count allocations, and the
count does not depend on which allocator serves them.

Measuring. Interleave the binaries and take medians, and compare their output
byte for byte in the same runs: a speed-up that changes an answer is not one.
`benchmarks/run_head_to_head.py` does that for the fastVEP side of the tables in the README and the manuscript and writes every run into
`benchmarks/results/<build>/`; run nothing else while it does (in the run behind those tables the five whole-genome runs on one thread took 43 to 55
seconds, and the 55 had the CPU time of the others and spent ten seconds waiting), and hold a laptop awake with `caffeinate -i`, because macOS
idle-sleeps between commands and suspends a background job while it does (one run of the five organisms took 18.7 seconds for 12.4 that way). A
run that needs sequences its sidecar cache lacks (`--hgvs` over a cache built
without `--fasta`) rewrites the cache, so checksum it before and after a series.
A comparison with Ensembl VEP that was measured with one build describes another build only if the two write the same bytes on its inputs, which
`validation/compare_builds.py` checks: the numbers of `docs/VEP_DIVERGENCE.md` against release 115.1 had gone stale that way, silently, when the
output changed to follow release 116.
An `amd64` container on Apple silicon runs under Rosetta 2, which tells you that
the code runs and gives ratios at best: it aborts on a thread stack of 1 << 50
bytes, where a kernel answers an error, and it showed the glibc behaviour above
at its worst.

`crates/fastvep-io/src/output.rs` shows the shape to copy: `format_csq_entry_into`
writes fields into a caller-owned `String` instead of returning one.
`crates/fastvep-sa/tests/query_path_allocations.rs` shows how to hold the line —
a counting global allocator asserting a query path stays within a budget. That
is the right way to make a performance property permanent; a benchmark in a
commit message is not.

## Claims

Do not assert a provider's file format, a VEP behaviour, or a biological fact
you have not checked in this session against a fixture in `tests/` or a real
file. `#98` shipped once with a rationale that was confidently wrong about
RefSeq GFF3 (it parents non-transcript records to genes, which the fix had
assumed impossible), and the review caught it, not the tests. If you cannot
verify a claim, say so where you make it rather than writing it as fact.
