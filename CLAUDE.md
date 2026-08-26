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
cargo test --workspace                                 # ~1000 tests, all must pass
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
moved to the rejected list and a `WHY_PRE_NN` string saying what was wrong with
what it holds. `#95` and `#98` are the worked examples.

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
gene-dense locus puts about a quarter of its instructions in `malloc`/`free`,
and output formatting — not consequence prediction — accounts for most of the
rest. When you touch the per-variant loop, the useful question is not "is this
branch fast" but "does this allocate per variant, per transcript, or per
allele": the loop runs once per (variant × overlapping transcript × allele), so
a `String` built there is built millions of times.

Two shapes account for most of what has been found there so far. One is work
done before deciding whether it is needed — a value built on entry as a fallback
the common path discards. The other is a per-run cost paid per item: resolving
the same column names, or rescanning the same list, once for every record
rather than once.

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
