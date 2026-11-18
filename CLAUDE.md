# Working in this repository

fastVEP predicts the functional consequence of genomic variants.
A clinician or researcher may act on the output, so **the failure mode that matters is not a crash, it is a plausible-looking wrong answer.**
Most fixes in `git log` are exactly that: a variant annotated against the wrong transcript, a gene that lost its symbol, a cache trusted after the parser that wrote it was fixed.
None of them showed up as an error message.

## Layout

Twelve crates, layered so that nothing below depends on anything above it:

```
fastvep-core            Consequence (41 SO terms), Allele, VariantType, Impact
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

`fastvep-annotate` and `fastvep-cli` both drive annotation, and their per-variant loops have drifted apart before.
A change to one is usually a change to both, so check the other before calling the work done.

## The loop

```bash
cargo fmt --all                                        # before every commit
cargo clippy --workspace --all-targets -- -D warnings  # what CI gates on
cargo test --workspace                                 # all must pass
```

CI runs exactly these three, in this order.
The lint set lives in `[workspace.lints]` in the root `Cargo.toml`, so a local `cargo clippy` sees what CI sees; the only difference is `-D`.

Keep that set narrow.
Every lint in it is one the workspace currently satisfies, which is what makes a warning worth reading.
A lint that fires on hundreds of existing lines trains everyone to scroll past the output.

## Conventions

**Commits.**
Conventional commits with a scope, and a subject that names the consequence rather than the edit:

```
fix(gff): read ncRNA_gene records, so non-coding genes keep their symbol (#98)
fix(cache): reject the pre-#90 cache format instead of trusting it (#95)
```

`cliff.toml` drafts release notes from these.
`CHANGELOG.md` is hand-curated: git-cliff writes the unreleased block, a human writes the prose, and released sections are never regenerated.
Bug-fix commits do not touch it; it is updated when a release is cut.

**Silencing a lint.**
`#[allow(...)]` carries a one-line reason at the site.

**Comments explain why, not what.**
Record the reasoning that is not recoverable from the code: which provider's file format forced a branch, which issue proved an assumption wrong, why the obvious simplification is unsafe.
When you fix a bug, leave behind the fact that made it a bug.

**Formatting-only commits** go in alone and get their hash added to `.git-blame-ignore-revs`.
Never mix a reformat with a behaviour change; it hides the behaviour change in review.

## Two hazards

**A parser fix does not reach users on its own.**
`fastvep annotate` writes a sidecar transcript cache next to the GFF3 and reuses it whenever its mtime is newer than the source.
A user who updates fastVEP but not their GFF3 keeps reading what the old parser wrote.
So a change to `crates/fastvep-cache/src/gff.rs` that changes what gets parsed needs the cache magic in `transcript_cache.rs` bumped.
Move the old magic to the rejected list with a `SUMMARY_PRE_NN` / `DETAIL_PRE_NN` pair saying what was wrong with what it holds (`#95` and `#98` are the worked examples).

The summary and the detail are two lengths of one fact, because the callers differ.
The sidecar path rebuilds from the GFF3 by itself, so it prints the summary and stops.
The explicit `--transcript-cache` path has to refuse the run, so it owes the user the detail.
Recovery advice belongs at the call site, not in the error: the sidecar has already rebuilt by the time it reports, and deleting a cache named on the command line produces no sidecar to replace it.

Rejecting is deliberate.
A stale cache is *intact*, so nothing in it distinguishes a gene that lost its symbol from one that never had one.
There is no safe way to read it, only a safe way to refuse it.

The cache format has its own invariants:

- The sidecar is a sequence of zstd frames of about 8 MB, with a table of them at the end (a skippable frame, in zstd's seekable layout).
  A run with more than two threads reads the table to decode frames side by side.
- The frames concatenate to the byte stream that earlier builds wrote as one frame, and a decoder skips the table, so a reader that takes the file as one stream (every build before #166) still reads it.
  A change to how a cache is cut keeps that property or bumps the magic.
- A cache is read to its end, where zstd checks each frame's checksum.
  A loader that stops at the last transcript lets a flipped byte, or a second cache appended behind the first, load silently.
- Nothing is allocated on the word of bytes read before their checksum is checked.
  Every `bincode` call on a cache has a limit (`bounded`), because a damaged byte that decodes to a length of a terabyte is otherwise an abort, with no word on which cache it was and no rebuild.
  So `save_cache` refuses a transcript that the limit would refuse at every load.
- The table is read only by a run that uses it, so a cache with a damaged table loads on two threads and is refused on three.

**Region-restricted transcript sets must never be cached.**
Reading a GFF3 through the tabix path returns only the features overlapping *this* VCF's variants.
That set is correct for this run and wrong for every other one, so `run_annotate` tracks `region_restricted` and skips the cache save.
Any new partial-load path has to set that flag.

## Performance

The release profile is tuned (`lto = "fat"`, `codegen-units = 1`).
What remains is allocation: a callgrind profile of the annotate path over a gene-dense locus put about a quarter of its instructions in `malloc`/`free`, and output formatting, not consequence prediction, accounted for most of the rest.
When you touch the per-variant loop, ask whether it allocates per variant, per transcript, or per allele.
The loop runs once per (variant x overlapping transcript x allele), so a `String` built there is built millions of times.

Three shapes account for most of what has been found:

- Work done before deciding whether it is needed, such as a fallback value built on entry that the common path discards.
- A per-run cost paid per item, such as resolving the same column names or rescanning the same list for every record.
- A choice made after the work it was meant to save.
  The pick (`--pick*`) used to run after a row was built for every transcript of the window (hundreds per variant at `--distance 500000`) and then deleted all but one per gene.
  It now runs on what the predictor said, before any row exists (`crates/fastvep-annotate/src/pick.rs`), and `crates/fastvep-cli/tests/pick_allocations.rs` holds that allocations follow the transcripts kept, not the ones around.

**Threads.**
Free what a thread allocated on that thread.
One thread freeing a batch that the others built cost 11 to 13% of the batch loop at 8 threads, so a record is annotated, rendered and freed by the thread that takes it (`render_batch`), and each batch is written on a thread of its own.
Peak memory should not grow with the thread count: check 1, 8 and the machine's cores (parallel cache decode is capped at 8 decoders because its peak grew).

**Allocator.**
The `fastvep` binary, including `fastvep web`, uses mimalloc (`crates/fastvep-cli/src/main.rs`).
The libraries, their tests and the `fastvep-web` binary keep the system allocator.
mimalloc does less work per allocation and avoids glibc's per-thread arena growth, one `mprotect` per 4 KB, which dominated decoding the transcript cache on 8 threads.
The price is some peak memory on Linux and about 1.2 GiB more address space under `ulimit -v`.
The counting-allocator tests wrap the system allocator, since the library crates never link mimalloc; allocation counts do not depend on which allocator serves them.

**Measuring.**
- Interleave the binaries, take medians, and compare output byte for byte in the same runs: a speed-up that changes an answer is not one.
- `benchmarks/run_head_to_head.py` does this for the fastVEP side of the README and manuscript tables and writes every run into `benchmarks/results/<build>/`.
- Run nothing else while it does, and hold a laptop awake with `caffeinate -i`: macOS idle-sleeps between commands and suspends a background job.
- A run that needs sequences its sidecar cache lacks (`--hgvs` over a cache built without `--fasta`) rewrites the cache, so checksum it before and after a series.
- A VEP comparison measured with one build describes another build only if the two write the same bytes on its inputs, which `validation/compare_builds.py` checks.
  Numbers in `docs/VEP_DIVERGENCE.md` went stale this way when the output changed to follow release 116.
- An `amd64` container on Apple silicon runs under Rosetta 2.
  It tells you the code runs and gives ratios at best, and it aborts on a thread stack of 1 << 50 bytes where a kernel answers an error.

`crates/fastvep-io/src/output.rs` shows the shape to copy: `format_csq_entry_into` writes fields into a caller-owned `String` instead of returning one.
`crates/fastvep-sa/tests/query_path_allocations.rs` shows how to hold the line: a counting global allocator asserting a query path stays within a budget.
That is how to make a performance property permanent; a benchmark in a commit message is not.

## Claims

Do not assert a provider's file format, a VEP behaviour, or a biological fact that you have not checked in this session against a fixture in `tests/` or a real file.
`#98` shipped with a rationale that was confidently wrong about RefSeq GFF3 (it parents non-transcript records to genes, which the fix had assumed impossible), and review caught it, not the tests.
If you cannot verify a claim, say so where you make it rather than writing it as fact.
