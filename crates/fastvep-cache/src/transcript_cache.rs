//! Binary transcript cache for fast startup.
//!
//! Serializes fully-built `Vec<Transcript>` (including spliced sequences)
//! to a compact binary format using bincode + zstd compression.
//! Subsequent loads skip GFF3 parsing, FASTA loading, and sequence building.

use anyhow::{Context, Result};
use fastvep_genome::Transcript;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;
use std::time::SystemTime;

/// Magic header for zstd-compressed caches (current format).
const CACHE_MAGIC_V2: &[u8; 8] = b"FSTVEP02";
/// Magic header for legacy gzip-compressed caches (read-only support).
const CACHE_MAGIC_V1: &[u8; 8] = b"FSTVEP01";

/// Save transcripts to a binary cache file (bincode + zstd).
pub fn save_cache(transcripts: &[Transcript], path: &Path) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Creating cache file: {}", path.display()))?;
    let writer = BufWriter::new(file);
    // zstd level 1: fast compression, still much better decompression than gzip
    let mut zst = zstd::Encoder::new(writer, 1)?;

    // Write magic header
    use std::io::Write;
    zst.write_all(CACHE_MAGIC_V2)?;

    // Serialize with bincode
    bincode::serialize_into(&mut zst, transcripts)
        .with_context(|| "Serializing transcripts to cache")?;

    zst.finish()?;
    Ok(())
}

/// Load transcripts from a binary cache file.
/// Supports both zstd (v2) and legacy gzip (v1) formats.
pub fn load_cache(path: &Path) -> Result<Vec<Transcript>> {
    let file = File::open(path)
        .with_context(|| format!("Opening cache file: {}", path.display()))?;
    let mut reader = BufReader::new(file);

    // Peek at the first bytes to detect format.
    // zstd frames start with 0x28B52FFD; gzip starts with 0x1F8B.
    use std::io::Read;
    let mut peek = [0u8; 4];
    reader.read_exact(&mut peek)
        .with_context(|| "Reading cache header")?;

    // Rewind so the decompressor sees the full stream
    use std::io::Seek;
    reader.seek(std::io::SeekFrom::Start(0))?;

    if peek[0..2] == [0x1F, 0x8B] {
        // Legacy gzip format (v1)
        load_cache_gzip(reader)
    } else {
        // zstd format (v2, or future)
        load_cache_zstd(reader)
    }
}

fn load_cache_zstd<R: std::io::Read>(reader: R) -> Result<Vec<Transcript>> {
    let mut zst = zstd::Decoder::new(reader)?;

    use std::io::Read;
    let mut magic = [0u8; 8];
    zst.read_exact(&mut magic)
        .with_context(|| "Reading cache header")?;
    if &magic != CACHE_MAGIC_V2 {
        anyhow::bail!("Invalid cache file (wrong magic header, expected FSTVEP02)");
    }

    let transcripts: Vec<Transcript> = bincode::deserialize_from(&mut zst)
        .with_context(|| "Deserializing transcripts from cache")?;
    Ok(transcripts)
}

fn load_cache_gzip<R: std::io::Read>(reader: R) -> Result<Vec<Transcript>> {
    use flate2::read::GzDecoder;

    let mut gz = GzDecoder::new(reader);

    use std::io::Read;
    let mut magic = [0u8; 8];
    gz.read_exact(&mut magic)
        .with_context(|| "Reading cache header")?;
    if &magic != CACHE_MAGIC_V1 {
        anyhow::bail!("Invalid cache file (wrong magic header, expected FSTVEP01)");
    }

    let transcripts: Vec<Transcript> = bincode::deserialize_from(&mut gz)
        .with_context(|| "Deserializing transcripts from cache")?;
    Ok(transcripts)
}

/// Check if cache file is newer than source file.
pub fn cache_is_fresh(cache_path: &Path, source_path: &Path) -> bool {
    let cache_mtime = cache_path
        .metadata()
        .and_then(|m| m.modified())
        .unwrap_or(SystemTime::UNIX_EPOCH);
    let source_mtime = source_path
        .metadata()
        .and_then(|m| m.modified())
        .unwrap_or(SystemTime::now());
    cache_mtime > source_mtime
}

/// Get the default cache path for a given GFF3 path.
pub fn default_cache_path(gff3_path: &Path) -> std::path::PathBuf {
    let mut cache_path = gff3_path.to_path_buf();
    let name = cache_path
        .file_name()
        .map(|n| {
            let s = n.to_string_lossy();
            if s.ends_with(".fastvep.cache") {
                s.to_string()
            } else {
                format!("{}.fastvep.cache", s)
            }
        })
        .unwrap_or_else(|| "transcripts.fastvep.cache".to_string());
    cache_path.set_file_name(name);
    cache_path
}

#[cfg(test)]
mod tests {
    use super::*;
    use fastvep_core::Strand;
    use fastvep_genome::{Exon, Gene, Transcript};
    use std::sync::Arc;
    use tempfile::NamedTempFile;

    fn make_test_transcript() -> Transcript {
        Transcript {
            stable_id: Arc::from("ENST00000001"),
            version: Some(1),
            gene: Gene {
                stable_id: Arc::from("ENSG00000001"),
                symbol: Some(Arc::from("TEST")),
                symbol_source: None,
                hgnc_id: None,
                biotype: Arc::from("protein_coding"),
                chromosome: Arc::from("1"),
                start: 1000,
                end: 5000,
                strand: Strand::Forward,
            },
            biotype: Arc::from("protein_coding"),
            chromosome: Arc::from("1"),
            start: 1000,
            end: 5000,
            strand: Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "ENSE001".into(),
                    start: 1000,
                    end: 1200,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: -1,
                    rank: 1,
                },
            ],
            translation: None,
            cdna_coding_start: Some(1),
            cdna_coding_end: Some(200),
            coding_region_start: Some(1000),
            coding_region_end: Some(1200),
            spliced_seq: Some("ACGTACGT".into()),
            translateable_seq: Some("ACGT".into()),
            peptide: Some("T".into()),
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: Some(1),
            appris: Some("P1".into()),
            ccds: None,
            protein_id: Some("ENSP001".into()),
            protein_version: Some(1),
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: None,
            gencode_primary: false,
            flags: vec![],
            codon_table_start_phase: 0,
        }
    }

    #[test]
    fn test_cache_roundtrip() {
        let transcripts = vec![make_test_transcript()];
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();

        save_cache(&transcripts, path).unwrap();
        let loaded = load_cache(path).unwrap();

        assert_eq!(loaded.len(), 1);
        assert_eq!(&*loaded[0].stable_id, "ENST00000001");
        assert_eq!(&**loaded[0].gene.symbol.as_ref().unwrap(), "TEST");
        assert_eq!(loaded[0].spliced_seq.as_deref(), Some("ACGTACGT"));
        assert_eq!(loaded[0].canonical, true);
        assert_eq!(loaded[0].tsl, Some(1));
    }

    #[test]
    fn test_invalid_magic() {
        let tmp = NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), b"NOTVALID").unwrap();
        assert!(load_cache(tmp.path()).is_err());
    }

    #[test]
    fn test_legacy_gzip_cache_loads() {
        // Create a legacy gzip cache and verify it still loads
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;

        let transcripts = vec![make_test_transcript()];
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();

        let file = File::create(path).unwrap();
        let writer = BufWriter::new(file);
        let mut gz = GzEncoder::new(writer, Compression::fast());
        gz.write_all(CACHE_MAGIC_V1).unwrap();
        bincode::serialize_into(&mut gz, &transcripts).unwrap();
        gz.finish().unwrap();

        let loaded = load_cache(path).unwrap();
        assert_eq!(loaded.len(), 1);
        assert_eq!(&*loaded[0].stable_id, "ENST00000001");
    }
}
