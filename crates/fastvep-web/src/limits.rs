//! What one request is allowed to cost.
//!
//! A response is built in memory: the annotation of every variant of the request, as JSON, before the
//! first byte goes out. On the 16 GB host of fastvep.org that is the one thing a single client
//! can use to take the server down. Measured on ClinVar variants against GRCh38 (a response is about
//! 38 KB a variant without `pick`, 14 KB with it, and 258 KB a variant with `acmg` and no `pick`), 10,000
//! variants with ACMG came to several gigabytes in the process, and a body of the size the server
//! accepts, 10 MiB, holds some 440,000 variants. The limits below turn that into a refusal that names
//! what to change, and keep a busy server answering instead of swapping.
//!
//! * **Response size.** The annotation is made in batches and serialised as it goes, and a request
//!   whose answer would pass [`Limits::max_response_bytes`] is refused with `413` the moment it does,
//!   with the work of one batch wasted and not the whole request. The batches are what bound the
//!   memory: a request holds the bytes of its answer and the values of one batch, where the whole
//!   answer as values is five to ten times the size of its text.
//! * **Time.** A request that runs past [`Limits::max_request_time`] is refused, also with `413`:
//!   what is wrong with it is its size.
//! * **Parallelism.** At most [`Limits::max_parallel`] annotations run at once, one for each core by
//!   default, because one annotation is one thread and more of them than cores only make each
//!   slower and hold more memory at once. A request that has waited [`Limits::queue_timeout`] for a
//!   place is refused with `503` and a `Retry-After`, so that a burst is told to come back and not
//!   held open until the proxy gives up on it.
//! * **A client that goes away** stops its work between batches, instead of finishing an answer that
//!   nobody will read.

use std::time::Duration;

/// How many variants are annotated, and serialised, at a time.
pub const DEFAULT_BATCH_VARIANTS: usize = 250;

#[derive(Debug, Clone)]
pub struct Limits {
    /// Largest request body, in bytes. Enforced by the router (`--max-body-size`); kept here so
    /// that `GET /api/status` can report all the limits in one place.
    pub max_body_bytes: usize,
    /// Largest JSON answer a request may produce, in bytes.
    pub max_response_bytes: usize,
    /// Longest an annotation may run.
    pub max_request_time: Duration,
    /// Annotations that may run at once.
    pub max_parallel: usize,
    /// How long a request waits for a place before it is turned away.
    pub queue_timeout: Duration,
    /// Variants annotated and serialised at a time.
    pub batch_variants: usize,
}

impl Default for Limits {
    fn default() -> Self {
        Self {
            max_body_bytes: 10 * 1024 * 1024,
            max_response_bytes: 64 * 1024 * 1024,
            max_request_time: Duration::from_secs(90),
            max_parallel: std::thread::available_parallelism().map_or(2, |n| n.get()),
            queue_timeout: Duration::from_secs(30),
            batch_variants: DEFAULT_BATCH_VARIANTS,
        }
    }
}

impl Limits {
    /// What `GET /api/status` tells a client, so that it can size its requests.
    pub fn describe(&self) -> serde_json::Value {
        serde_json::json!({
            "max_body_bytes": self.max_body_bytes,
            "max_response_bytes": self.max_response_bytes,
            "max_request_seconds": self.max_request_time.as_secs(),
            "max_parallel": self.max_parallel,
        })
    }
}
