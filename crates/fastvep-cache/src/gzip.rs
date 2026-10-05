//! Sequential reads that may be gzip-compressed.
//!
//! Same rule as the VCF opener in `fastvep-cli`: the gzip magic bytes decide,
//! and a `.gz` / `.bgz` suffix is the fallback. The two peeked bytes are put
//! back in front of the stream, so a plain-text file is unchanged.

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, Read};
use std::path::Path;

/// Open `path` for a sequential read, decompressing when the file is gzip.
pub fn open_maybe_gzip(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path).with_context(|| format!("Opening {}", path.display()))?;
    wrap_maybe_gzip(file, &path.to_string_lossy())
}

/// Decompress `reader` when `source` is gzip, otherwise return it unchanged.
///
/// `source` is the path the bytes came from, or `"-"` for stdin. A stdin
/// stream has no filename, so only the magic bytes apply there.
pub fn wrap_maybe_gzip(mut reader: impl Read + 'static, source: &str) -> Result<Box<dyn Read>> {
    let mut prefix = [0u8; 2];
    let bytes_read = reader.read(&mut prefix)?;
    let looks_like_gzip = bytes_read == 2 && prefix == [0x1f, 0x8b];

    let replay = io::Cursor::new(prefix[..bytes_read].to_vec()).chain(reader);
    if looks_like_gzip || source_claims_gzip(source) {
        Ok(Box::new(MultiGzDecoder::new(replay)))
    } else {
        Ok(Box::new(replay))
    }
}

/// Whether `path` is gzip, by magic bytes or by a `.gz` / `.bgz` suffix.
///
/// A `.fai` stores offsets into the uncompressed FASTA. Applying it to a
/// gzip member reads compressed bytes, so callers use this to refuse the
/// memory-mapped path.
pub fn is_gzip_compressed(path: &Path) -> Result<bool> {
    let mut file = File::open(path).with_context(|| format!("Opening {}", path.display()))?;
    let mut prefix = [0u8; 2];
    let n = file.read(&mut prefix)?;
    Ok((n == 2 && prefix == [0x1f, 0x8b]) || source_claims_gzip(&path.to_string_lossy()))
}

fn source_claims_gzip(source: &str) -> bool {
    source != "-" && (source.ends_with(".gz") || source.ends_with(".bgz"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Read;

    fn gzip(data: &str) -> Vec<u8> {
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        std::io::Write::write_all(&mut enc, data.as_bytes()).unwrap();
        enc.finish().unwrap()
    }

    fn read_all(path: &Path) -> String {
        let mut reader = open_maybe_gzip(path).unwrap();
        let mut s = String::new();
        reader.read_to_string(&mut s).unwrap();
        s
    }

    #[test]
    fn plain_text_keeps_its_first_bytes() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("genes.gff3");
        std::fs::write(&path, "##gff-version 3\n").unwrap();
        assert_eq!(read_all(&path), "##gff-version 3\n");
        assert!(!is_gzip_compressed(&path).unwrap());
    }

    #[test]
    fn magic_bytes_decompress_a_file_that_is_not_named_gz() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("genes.gff3");
        std::fs::write(&path, gzip("##gff-version 3\n")).unwrap();
        assert!(is_gzip_compressed(&path).unwrap());
        assert_eq!(read_all(&path), "##gff-version 3\n");
    }

    #[test]
    fn a_bgz_suffix_is_decompressed_too() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("genes.gff3.bgz");
        std::fs::write(&path, gzip(">chr1\nACGT\n")).unwrap();
        assert_eq!(read_all(&path), ">chr1\nACGT\n");
    }
}
