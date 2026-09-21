//! Differential test: `IndexedTranscriptProvider` must return exactly the set a
//! naive overlap filter returns, for every query, on real gene models.
//!
//! Point the test at one or more GFF3 files with `FASTVEP_TEST_GFF3`
//! (colon-separated). Without it the test is a no-op, so CI stays hermetic.
//!
//! The reference answer comes from a sweep line rather than from re-filtering
//! every transcript per query: the query points are dense (hundreds of
//! thousands per contig) and an O(n) reference per query would not finish. The
//! sweep is itself cross-checked against `MemoryTranscriptProvider`'s plain
//! filter on a subsample, so a bug in the reference cannot hide a bug in the
//! index.

use std::cmp::Reverse;
use std::collections::{BTreeSet, BinaryHeap};
use std::fs::File;
use std::io::BufReader;

use fastvep_cache::gff::parse_gff3;
use fastvep_cache::providers::{
    IndexedTranscriptProvider, MemoryTranscriptProvider, TranscriptProvider,
};
use fastvep_genome::Transcript;

type Key = (String, u64, u64);

fn load(path: &str) -> Vec<Transcript> {
    let file = File::open(path).unwrap_or_else(|e| panic!("open {path}: {e}"));
    if path.ends_with(".gz") {
        parse_gff3(flate2::read::MultiGzDecoder::new(BufReader::new(file)))
            .unwrap_or_else(|e| panic!("parse {path}: {e}"))
    } else {
        parse_gff3(BufReader::new(file)).unwrap_or_else(|e| panic!("parse {path}: {e}"))
    }
}

fn key(t: &Transcript) -> Key {
    (t.stable_id.to_string(), t.start, t.end)
}

/// Every transcript boundary, nudged by +/-2, plus a coarse sweep of the whole
/// contig extended past its last base. Boundaries are where an off-by-one in
/// the scan window shows; the sweep catches the rest, including the region past
/// the last transcript start where issue #126 lives.
fn query_points(trs: &[&Transcript], contig_len: u64) -> Vec<u64> {
    let mut pts: BTreeSet<u64> = BTreeSet::new();
    for t in trs {
        for anchor in [t.start, t.end] {
            for delta in [-2i64, -1, 0, 1, 2] {
                let p = anchor as i64 + delta;
                if p >= 1 {
                    pts.insert(p as u64);
                }
            }
        }
    }
    let step = (contig_len / 20_000).max(1);
    let mut p = 1;
    while p <= contig_len + 10_000 {
        pts.insert(p);
        p += step;
    }
    pts.into_iter().collect()
}

/// Expected overlap sets for `points` (ascending) at a fixed query width, by a
/// sweep line over transcripts sorted by start. Everything in the heap at
/// point `p` satisfies `start <= p + width` (it was pushed) and `end >= p`
/// (it was not popped) — exactly the overlap predicate.
fn expected_by_sweep(sorted: &[&Transcript], points: &[u64], width: u64) -> Vec<BTreeSet<Key>> {
    let mut out = Vec::with_capacity(points.len());
    let mut next = 0usize;
    let mut active: BinaryHeap<Reverse<(u64, usize)>> = BinaryHeap::new();
    for &p in points {
        while next < sorted.len() && sorted[next].start <= p + width {
            active.push(Reverse((sorted[next].end, next)));
            next += 1;
        }
        while let Some(&Reverse((e, _))) = active.peek() {
            if e < p {
                active.pop();
            } else {
                break;
            }
        }
        out.push(
            active
                .iter()
                .map(|Reverse((_, i))| key(sorted[*i]))
                .collect(),
        );
    }
    out
}

#[test]
fn indexed_provider_matches_naive_overlap_on_real_gene_models() {
    let Ok(paths) = std::env::var("FASTVEP_TEST_GFF3") else {
        eprintln!("FASTVEP_TEST_GFF3 unset; skipping");
        return;
    };

    let mut total_queries = 0usize;
    let mut bad_queries = 0usize;
    let mut dropped_rows = 0usize;
    let mut extra_rows = 0usize;
    let mut bad_sites: BTreeSet<(String, u64)> = BTreeSet::new();
    let mut examples: Vec<String> = Vec::new();

    for path in paths.split(':').filter(|p| !p.is_empty()) {
        let transcripts = load(path);
        assert!(!transcripts.is_empty(), "{path} yielded no transcripts");
        let indexed = IndexedTranscriptProvider::new(transcripts.clone());
        let linear = MemoryTranscriptProvider::new(transcripts.clone());

        let mut chroms: BTreeSet<String> = BTreeSet::new();
        for t in &transcripts {
            chroms.insert(t.chromosome.to_string());
        }

        for chrom in &chroms {
            let mut sorted = linear.get_transcripts_by_chrom(chrom).unwrap();
            sorted.sort_by_key(|t| t.start);
            let contig_len = sorted.iter().map(|t| t.end).max().unwrap_or(0);
            let points = query_points(&sorted, contig_len);

            for width in [0u64, 1, 50, 1_000] {
                let expected = expected_by_sweep(&sorted, &points, width);

                // Cross-check the sweep itself against the plain filter, on a
                // subsample (the plain filter is O(n) per query).
                for (i, &p) in points.iter().enumerate().step_by(points.len() / 40 + 1) {
                    let naive: BTreeSet<Key> = linear
                        .get_transcripts(chrom, p, p + width)
                        .unwrap()
                        .into_iter()
                        .map(key)
                        .collect();
                    assert_eq!(
                        naive, expected[i],
                        "sweep reference disagrees with the plain filter at \
                         {chrom}:{p}+{width} — the test itself is wrong"
                    );
                }

                for (i, &p) in points.iter().enumerate() {
                    let actual: BTreeSet<Key> = indexed
                        .get_transcripts(chrom, p, p + width)
                        .unwrap()
                        .into_iter()
                        .map(key)
                        .collect();
                    total_queries += 1;
                    if actual != expected[i] {
                        bad_queries += 1;
                        dropped_rows += expected[i].difference(&actual).count();
                        extra_rows += actual.difference(&expected[i]).count();
                        bad_sites.insert((chrom.clone(), p));
                        if examples.len() < 10 {
                            examples.push(format!(
                                "{chrom}:{p}+{width} expected {} got {} missing {:?}",
                                expected[i].len(),
                                actual.len(),
                                expected[i]
                                    .difference(&actual)
                                    .map(|(id, s, e)| format!("{id}[{s}-{e}]"))
                                    .take(3)
                                    .collect::<Vec<_>>()
                            ));
                        }
                    }
                }
            }
        }
    }

    eprintln!(
        "queries={total_queries} wrong={bad_queries} sites={} dropped_rows={dropped_rows} \
         extra_rows={extra_rows}",
        bad_sites.len()
    );
    for e in &examples {
        eprintln!("  {e}");
    }
    assert_eq!(
        bad_queries,
        0,
        "{bad_queries}/{total_queries} queries disagree with a naive overlap scan at {} sites, \
         dropping {dropped_rows} transcript rows",
        bad_sites.len()
    );
}

/// The shape from issue #126, taken from the real model: query just past the
/// end of the last transcript by start, where a longer earlier transcript is
/// still open.
#[test]
fn nested_long_transcript_survives_backward_scan() {
    let Ok(paths) = std::env::var("FASTVEP_TEST_GFF3") else {
        return;
    };
    for path in paths.split(':').filter(|p| !p.is_empty()) {
        let transcripts = load(path);
        let indexed = IndexedTranscriptProvider::new(transcripts.clone());
        let linear = MemoryTranscriptProvider::new(transcripts.clone());
        let mut chroms: BTreeSet<String> = BTreeSet::new();
        for t in &transcripts {
            chroms.insert(t.chromosome.to_string());
        }
        for chrom in &chroms {
            let mut sorted = linear.get_transcripts_by_chrom(chrom).unwrap();
            sorted.sort_by_key(|t| t.start);
            let Some(last) = sorted.last() else { continue };
            let max_end = sorted.iter().map(|t| t.end).max().unwrap();
            if last.end >= max_end {
                continue; // no transcript reaches past the last one by start
            }
            for pos in [last.end + 1, (last.end + max_end) / 2, max_end] {
                let expected: BTreeSet<Key> = linear
                    .get_transcripts(chrom, pos, pos)
                    .unwrap()
                    .into_iter()
                    .map(key)
                    .collect();
                let actual: BTreeSet<Key> = indexed
                    .get_transcripts(chrom, pos, pos)
                    .unwrap()
                    .into_iter()
                    .map(key)
                    .collect();
                assert_eq!(expected, actual, "{path} {chrom}:{pos}");
            }
        }
    }
}
