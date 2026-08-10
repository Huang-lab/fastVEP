//! `AnnotationContext::sa_providers` must allow genuinely concurrent lookups.
//!
//! The providers used to be stored as `Vec<Mutex<Box<dyn AnnotationProvider>>>`
//! and the library/web annotate path locked each one per allele per variant, so
//! concurrent `/annotate` requests serialized on every source in turn - even
//! though every `AnnotationProvider` method takes `&self` and the trait requires
//! `Send + Sync`. (The CLI pipeline sidestepped it by unwrapping the mutexes
//! immediately after loading, which is the tell that they were never needed.)
//!
//! This is proved without timing: the provider below blocks inside
//! `annotate_position` until every worker has entered it. Under an external
//! per-provider mutex the first worker in would wait for peers that cannot
//! acquire the lock, so the rendezvous could never complete.

use anyhow::Result;
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue, SaMetadata};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

const WORKERS: usize = 8;

/// How long a worker waits for its peers before giving up. Only reached when the
/// providers serialize, in which case the whole test costs `WORKERS` × this and
/// then fails.
const RENDEZVOUS_TIMEOUT: Duration = Duration::from_secs(2);

#[derive(Default)]
struct Counters {
    /// Monotonic: how many callers have ever entered. Drives the rendezvous, so
    /// a worker that has already been released cannot un-release its peers.
    arrived: AtomicUsize,
    /// Up on entry, down on exit: current occupancy.
    inside: AtomicUsize,
    /// High-water mark of `inside` - the actual measure of simultaneity.
    peak_inside: AtomicUsize,
}

/// Counts how many callers are inside `annotate_position` at once and refuses to
/// return until all of them are, or the deadline passes.
struct RendezvousProvider {
    metadata: SaMetadata,
    counters: Arc<Counters>,
}

impl RendezvousProvider {
    fn new(counters: Arc<Counters>) -> Self {
        Self {
            metadata: SaMetadata {
                name: "Rendezvous".into(),
                version: "1".into(),
                description: "test".into(),
                assembly: "GRCh38".into(),
                json_key: "rendezvous".into(),
                match_by_allele: true,
                is_array: false,
                is_positional: false,
            },
            counters,
        }
    }
}

impl AnnotationProvider for RendezvousProvider {
    fn name(&self) -> &str {
        &self.metadata.name
    }
    fn json_key(&self) -> &str {
        &self.metadata.json_key
    }
    fn metadata(&self) -> &SaMetadata {
        &self.metadata
    }

    fn annotate_position(
        &self,
        _chrom: &str,
        _pos: u64,
        _ref_allele: &str,
        _alt_allele: &str,
    ) -> Result<Option<AnnotationValue>> {
        self.counters.arrived.fetch_add(1, Ordering::SeqCst);
        let now = self.counters.inside.fetch_add(1, Ordering::SeqCst) + 1;
        self.counters.peak_inside.fetch_max(now, Ordering::SeqCst);

        // Bounded rendezvous: wait for everyone, but give up after
        // `RENDEZVOUS_TIMEOUT` so a regression fails the assertion below instead
        // of hanging the suite. `inside` is decremented on the way out, so
        // `peak_inside` measures true simultaneity - serialized callers each see
        // a count of 1 no matter how many of them eventually run.
        let deadline = Instant::now() + RENDEZVOUS_TIMEOUT;
        while self.counters.arrived.load(Ordering::SeqCst) < WORKERS && Instant::now() < deadline {
            std::thread::yield_now();
        }
        self.counters.inside.fetch_sub(1, Ordering::SeqCst);
        Ok(None)
    }
}

#[test]
fn all_providers_can_be_queried_at_once() {
    let counters = Arc::new(Counters::default());
    // Exactly the storage shape `AnnotationContext::sa_providers` uses.
    let providers: Vec<Box<dyn AnnotationProvider>> =
        vec![Box::new(RendezvousProvider::new(Arc::clone(&counters)))];

    std::thread::scope(|scope| {
        for _ in 0..WORKERS {
            scope.spawn(|| {
                for sa in &providers {
                    sa.annotate_position("chr1", 100, "A", "G").unwrap();
                }
            });
        }
    });

    let peak = counters.peak_inside.load(Ordering::SeqCst);
    assert_eq!(
        peak, WORKERS,
        "expected all {WORKERS} workers inside annotate_position simultaneously, \
         peaked at {peak}; providers are being serialized"
    );
}
