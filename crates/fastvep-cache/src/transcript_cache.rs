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
///
/// The write is atomic: the stream goes to a temporary file in the same
/// directory and is renamed into place once complete. `File::create` on the
/// destination would publish the path the instant writing began, and both the
/// magic header and the zstd frame header are written first, so a reader
/// arriving mid-write got a file that passed every check we make and then
/// failed deep inside deserialization. Two annotations started back to back
/// against a fresh GFF3 is enough to hit it, and the loser degrades to a
/// plausible-looking, wrong annotation rather than an error. `rename` within
/// one directory is atomic, so a reader now sees either the previous cache or
/// the complete new one.
pub fn save_cache(transcripts: &[Transcript], path: &Path) -> Result<()> {
    // The temp file has to live in the destination's own directory: `rename` is
    // only atomic within a filesystem, and /tmp is frequently a different one.
    let dir = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));

    let tmp = tempfile::Builder::new()
        .prefix(".fastvep-cache-")
        .suffix(".partial")
        .tempfile_in(dir)
        .with_context(|| format!("Creating temporary cache file in {}", dir.display()))?;

    let writer = BufWriter::new(tmp.reopen().with_context(|| {
        format!(
            "Reopening temporary cache file {}",
            tmp.path().to_path_buf().display()
        )
    })?);
    // zstd level 1: fast compression, still much better decompression than gzip.
    let mut zst = zstd::Encoder::new(writer, 1)?;
    // Record a frame checksum so a cache corrupted in place - a bad sector, a
    // half-flushed page - fails loudly at load instead of deserializing into
    // something subtly wrong.
    zst.include_checksum(true)?;

    // Write magic header
    use std::io::Write;
    zst.write_all(CACHE_MAGIC_V2)?;

    // Serialize with bincode
    bincode::serialize_into(&mut zst, transcripts)
        .with_context(|| "Serializing transcripts to cache")?;

    let mut writer = zst.finish()?;
    writer.flush().with_context(|| "Flushing cache writer")?;
    // fsync before the rename. Without it a crash between the two can leave the
    // renamed path pointing at unwritten blocks, which is the same silent
    // partial-cache failure by a slower route.
    writer
        .get_ref()
        .sync_all()
        .with_context(|| "Syncing cache file to disk")?;
    drop(writer);

    // tempfile creates at 0600. `File::create` produced 0666 & ~umask, so
    // publishing the temp file as-is would quietly make a cache unreadable to
    // everyone but its author - which matters where one user builds a cache in
    // a shared reference directory for a whole group. Keep the mode the
    // destination already had, and otherwise fall back to the 0644 a standard
    // umask would have given.
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mode = std::fs::metadata(path)
            .map(|m| m.permissions().mode() & 0o7777)
            .unwrap_or(0o644);
        std::fs::set_permissions(tmp.path(), std::fs::Permissions::from_mode(mode)).with_context(
            || {
                format!(
                    "Setting permissions on the new cache for {}",
                    path.display()
                )
            },
        )?;
    }

    tmp.persist(path)
        .map_err(|e| anyhow::anyhow!("Publishing cache file {}: {}", path.display(), e.error))?;
    Ok(())
}

/// Load transcripts from a binary cache file.
/// Supports both zstd (v2) and legacy gzip (v1) formats.
pub fn load_cache(path: &Path) -> Result<Vec<Transcript>> {
    let file =
        File::open(path).with_context(|| format!("Opening cache file: {}", path.display()))?;
    let mut reader = BufReader::new(file);

    // Peek at the first bytes to detect format.
    // zstd frames start with 0x28B52FFD; gzip starts with 0x1F8B.
    use std::io::Read;
    let mut peek = [0u8; 4];
    reader
        .read_exact(&mut peek)
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
            exons: vec![Exon {
                stable_id: "ENSE001".into(),
                start: 1000,
                end: 1200,
                strand: Strand::Forward,
                phase: 0,
                end_phase: -1,
                rank: 1,
            }],
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
        assert!(loaded[0].canonical);
        assert_eq!(loaded[0].tsl, Some(1));
    }

    /// A poorly-compressible transcript set, large enough that the cache write
    /// takes long enough for a concurrent reader to be scheduled during it.
    fn bulky_transcripts(n: usize) -> Vec<Transcript> {
        (0..n)
            .map(|i| {
                let mut tr = make_test_transcript();
                // Vary the sequence per transcript, or zstd folds the whole set
                // into a few kilobytes and there is no write window at all.
                let mut seq = String::with_capacity(4096);
                let mut x = (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
                for _ in 0..4096 {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    seq.push(match x % 4 {
                        0 => 'A',
                        1 => 'C',
                        2 => 'G',
                        _ => 'T',
                    });
                }
                tr.stable_id = Arc::from(format!("ENST{i:08}").as_str());
                tr.spliced_seq = Some(seq);
                tr
            })
            .collect()
    }

    #[test]
    fn a_truncated_cache_is_rejected_rather_than_read_short() {
        let transcripts = bulky_transcripts(64);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        save_cache(&transcripts, &path).unwrap();

        let full = std::fs::read(&path).unwrap();
        for fraction in [1, 2, 10, 50, 90, 99] {
            let cut = full.len() * fraction / 100;
            std::fs::write(&path, &full[..cut]).unwrap();
            let loaded = load_cache(&path);
            assert!(
                loaded.is_err(),
                "a cache truncated to {fraction}% loaded anyway, as {:?} transcripts",
                loaded.map(|t| t.len())
            );
        }

        // And the intact file still loads, so the check is not just refusing
        // everything.
        std::fs::write(&path, &full).unwrap();
        assert_eq!(load_cache(&path).unwrap().len(), transcripts.len());
    }

    #[test]
    fn a_cache_write_is_never_visible_half_finished() {
        // The failure this guards: `File::create` publishes the destination the
        // instant writing starts, and the magic header goes out first, so a
        // second annotation launched while the first was still writing read a
        // prefix that passed the header check and then died inside
        // deserialization - or, with a cache it could not rebuild, degraded to
        // an all-intergenic answer at exit code 0.
        let transcripts = bulky_transcripts(4_000);
        let expected = transcripts.len();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("concurrent.cache");

        let writer_path = path.clone();
        let writer = std::thread::spawn(move || {
            save_cache(&transcripts, &writer_path).unwrap();
        });

        let mut absent = 0usize;
        let mut complete = 0usize;
        while !writer.is_finished() {
            // Decide "published or not" from this one open, not from a second
            // stat: the file can appear between the two, and then a legitimate
            // not-yet-there reads as a corrupt cache.
            match File::open(&path) {
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => absent += 1,
                Err(e) => panic!("unexpected error opening {}: {e}", path.display()),
                Ok(_) => match load_cache(&path) {
                    Ok(trs) => {
                        assert_eq!(
                            trs.len(),
                            expected,
                            "observed a published cache holding {} of {expected} transcripts",
                            trs.len()
                        );
                        complete += 1;
                    }
                    Err(e) => panic!(
                        "a cache file was visible at {} but would not load: {e}",
                        path.display()
                    ),
                },
            }
        }
        writer.join().unwrap();

        assert_eq!(load_cache(&path).unwrap().len(), expected);
        assert!(
            absent + complete > 0,
            "the poll loop never ran, so this proved nothing"
        );
        // No temp file may be left behind on the happy path.
        let leftovers: Vec<_> = std::fs::read_dir(dir.path())
            .unwrap()
            .filter_map(|e| e.ok())
            .map(|e| e.file_name().to_string_lossy().into_owned())
            .filter(|n| n.ends_with(".partial"))
            .collect();
        assert!(
            leftovers.is_empty(),
            "left temp files behind: {leftovers:?}"
        );
    }

    #[test]
    #[cfg(unix)]
    fn a_published_cache_keeps_the_mode_the_old_writer_gave_it() {
        // A cache in a shared reference directory has to stay readable by the
        // group that did not build it. tempfile creates at 0600, so publishing
        // the temp file unmodified would have silently locked everyone else out.
        use std::os::unix::fs::PermissionsExt;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("modes.cache");

        save_cache(&make_transcripts(), &path).unwrap();
        let mode = std::fs::metadata(&path).unwrap().permissions().mode() & 0o777;
        assert_eq!(mode, 0o644, "fresh cache published as {mode:o}, want 644");

        // Replacing an existing cache preserves whatever mode it had.
        std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o640)).unwrap();
        save_cache(&make_transcripts(), &path).unwrap();
        let mode = std::fs::metadata(&path).unwrap().permissions().mode() & 0o777;
        assert_eq!(
            mode, 0o640,
            "replacement cache published as {mode:o}, want 640"
        );
    }

    fn make_transcripts() -> Vec<Transcript> {
        vec![make_test_transcript()]
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
