//! The work behind each `fastvep` subcommand.
//!
//! One module per command, so a change to `sa-build` does not mean scrolling
//! past the annotate loop. Everything here was one 5,154-line file; this module
//! keeps only what more than one command needs.

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, Read};

pub mod annotate;
pub mod cache_build;
pub mod custom;
pub mod filter;
pub mod pick;
pub mod sa_build;

pub use annotate::{run_annotate, AnnotateConfig};
pub use cache_build::{parse_gff3_arg, run_cache_build, Gff3Spec};
pub use custom::run_oga_build;
pub use filter::run_filter;
pub use pick::{parse_pick_order, PickCriterion, DEFAULT_PICK_ORDER};
pub use sa_build::{
    run_sa_build, run_sa_build_format, run_sa_build_v2, run_sa_convert, source_from_json_key,
    source_has_decomposed_osa2, source_supports_osa2, OSA2_DECOMPOSED_SOURCES,
    OSA2_SUPPORTED_SOURCES,
};

pub(crate) fn open_vcf_input_reader(input: &str) -> Result<Box<dyn io::Read>> {
    let reader: Box<dyn io::Read> = if input == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(input).with_context(|| format!("Opening input file: {}", input))?)
    };

    wrap_maybe_gzip_reader(reader, input)
}

pub(crate) fn wrap_maybe_gzip_reader(
    mut reader: Box<dyn io::Read>,
    source: &str,
) -> Result<Box<dyn io::Read>> {
    let mut prefix = [0u8; 2];
    let bytes_read = reader.read(&mut prefix)?;
    let looks_like_gzip = bytes_read == 2 && prefix == [0x1f, 0x8b];

    let replay = io::Cursor::new(prefix[..bytes_read].to_vec()).chain(reader);
    if looks_like_gzip || (source != "-" && source.ends_with(".gz")) {
        Ok(Box::new(MultiGzDecoder::new(replay)))
    } else {
        Ok(Box::new(replay))
    }
}
