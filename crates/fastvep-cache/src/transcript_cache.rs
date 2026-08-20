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
///
/// V3 differs from V2 in nothing but the guarantee it carries: every V3 cache
/// was written by a build that cannot persist a region-restricted transcript
/// set. See [`load_cache_zstd`] for why that has to be a format bump rather
/// than a note in the changelog.
const CACHE_MAGIC_V3: &[u8; 8] = b"FSTVEP03";
/// Magic header for zstd caches written before #90 (rejected, see
/// [`load_cache_zstd`]).
const CACHE_MAGIC_V2: &[u8; 8] = b"FSTVEP02";
/// Magic header for legacy gzip-compressed caches (rejected, same reason).
const CACHE_MAGIC_V1: &[u8; 8] = b"FSTVEP01";

/// What a rejected pre-#90 cache says, in both the sidecar and the explicit
/// `--transcript-cache` path. Written once so the two read the same.
fn stale_format_error(found: &str) -> anyhow::Error {
    anyhow::anyhow!(
        "Transcript cache is in the {found} format, which predates the #90 fix and \
         cannot be trusted. Caches written between 2026-06-10 (the tabix read path) \
         and 2026-08-18 (#90) may hold only the transcripts overlapping one input \
         VCF's variants while looking like a whole-file cache to every later run - \
         the #87 failure, which reported 47,196 of 47,196 chr17 variants as \
         intergenic at exit code 0. Nothing in the file distinguishes the two, so \
         it is rejected rather than read. Rebuild with `fastvep cache`, or delete \
         it and let the sidecar rebuild itself."
    )
}

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
    zst.write_all(CACHE_MAGIC_V3)?;

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

/// Read a zstd cache, accepting only the current format.
///
/// A `FSTVEP02` cache is rejected, and this is deliberate even though the bytes
/// after the header are identical. #90 stopped *writing* region-restricted
/// transcript sets to the sidecar path, but every cache already on disk from a
/// build between 2026-06-10 and 2026-08-18 may be one, and there is no way to
/// tell from the file: the transcript count is plausible either way, the path is
/// keyed on the GFF3 alone, and `cache_is_fresh` only compares mtimes, so a
/// poisoned cache written after its GFF3 was downloaded passes every check we
/// make. The failure it produces is the one #87 measured - a complete,
/// well-formed VCF calling every variant intergenic at exit code 0.
///
/// So the fix has to be a format bump. Rejecting costs one rebuild, which the
/// sidecar path does automatically; reading costs a silently wrong annotation
/// that nothing downstream can detect. Same trade as #88's decision to error
/// rather than annotate against an empty transcript set.
fn load_cache_zstd<R: std::io::Read>(reader: R) -> Result<Vec<Transcript>> {
    let mut zst = zstd::Decoder::new(reader)?;

    use std::io::Read;
    let mut magic = [0u8; 8];
    zst.read_exact(&mut magic)
        .with_context(|| "Reading cache header")?;
    if &magic == CACHE_MAGIC_V2 {
        return Err(stale_format_error("FSTVEP02"));
    }
    if &magic != CACHE_MAGIC_V3 {
        anyhow::bail!("Invalid cache file (wrong magic header, expected FSTVEP03)");
    }

    let transcripts: Vec<Transcript> = bincode::deserialize_from(&mut zst)
        .with_context(|| "Deserializing transcripts from cache")?;
    Ok(transcripts)
}

/// Legacy gzip caches are recognised and rejected.
///
/// V1 predates the zstd format and therefore also predates #90, so it carries
/// the same ambiguity as V2 and gets the same answer. Recognising the magic
/// rather than falling through to "invalid cache" is the point: a user with a
/// years-old cache should be told to rebuild, not told their file is corrupt.
fn load_cache_gzip<R: std::io::Read>(reader: R) -> Result<Vec<Transcript>> {
    use flate2::read::GzDecoder;

    let mut gz = GzDecoder::new(reader);

    use std::io::Read;
    let mut magic = [0u8; 8];
    gz.read_exact(&mut magic)
        .with_context(|| "Reading cache header")?;
    if &magic == CACHE_MAGIC_V1 {
        return Err(stale_format_error("FSTVEP01 (gzip)"));
    }
    anyhow::bail!("Invalid cache file (wrong magic header, expected FSTVEP01)");
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

    /// Write a cache in an older format, byte-for-byte as the builds of the day
    /// wrote it. Both are readable: rejecting them is a decision, not a
    /// limitation, which is exactly why it needs a test.
    fn write_v2_cache(path: &Path, transcripts: &[Transcript]) {
        use std::io::Write;
        let file = File::create(path).unwrap();
        let mut zst = zstd::Encoder::new(BufWriter::new(file), 1).unwrap();
        zst.write_all(CACHE_MAGIC_V2).unwrap();
        bincode::serialize_into(&mut zst, transcripts).unwrap();
        zst.finish().unwrap().flush().unwrap();
    }

    fn write_v1_cache(path: &Path, transcripts: &[Transcript]) {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;
        let file = File::create(path).unwrap();
        let mut gz = GzEncoder::new(BufWriter::new(file), Compression::fast());
        gz.write_all(CACHE_MAGIC_V1).unwrap();
        bincode::serialize_into(&mut gz, transcripts).unwrap();
        gz.finish().unwrap();
    }

    #[test]
    fn a_pre_90_cache_is_rejected_rather_than_read() {
        // A well-formed FSTVEP02 cache holding a plausible transcript set. The
        // point is that this is indistinguishable from a region-restricted one
        // written by the same build: same magic, same bincode, a transcript
        // count that looks fine, and a path keyed on the GFF3 alone. #90 stopped
        // writing them; it could do nothing about the ones already on disk.
        let tmp = NamedTempFile::new().unwrap();
        write_v2_cache(tmp.path(), &make_transcripts());

        let err = load_cache(tmp.path()).expect_err("FSTVEP02 must not load");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("FSTVEP02"),
            "the error should name the format found, got: {msg}"
        );
        assert!(
            msg.contains("fastvep cache"),
            "the error should say how to recover, got: {msg}"
        );
        assert!(
            msg.contains("intergenic"),
            "the error should say what goes wrong if it were read, got: {msg}"
        );
    }

    #[test]
    fn a_legacy_gzip_cache_is_rejected_with_the_same_reason() {
        // V1 predates the zstd format and therefore predates the tabix read
        // path too, so it carries the same ambiguity. Recognised specifically,
        // rather than falling through to "invalid cache": a user with an old
        // cache should be told to rebuild, not told their file is corrupt.
        let tmp = NamedTempFile::new().unwrap();
        write_v1_cache(tmp.path(), &make_transcripts());

        let err = load_cache(tmp.path()).expect_err("FSTVEP01 must not load");
        let msg = format!("{err:#}");
        assert!(msg.contains("FSTVEP01"), "got: {msg}");
        assert!(msg.contains("fastvep cache"), "got: {msg}");
    }

    #[test]
    fn a_v3_cache_written_now_round_trips() {
        // The other half of the bump: rejecting V2 is only acceptable if the
        // current writer produces something the current reader accepts.
        let tmp = NamedTempFile::new().unwrap();
        let transcripts = make_transcripts();
        save_cache(&transcripts, tmp.path()).unwrap();

        let loaded = load_cache(tmp.path()).unwrap();
        assert_eq!(loaded.len(), transcripts.len());
        assert_eq!(&*loaded[0].stable_id, "ENST00000001");
    }

    #[test]
    fn the_published_magic_is_v3() {
        // Pin the on-disk byte, so a future change to the writer that forgets
        // the reader shows up here rather than in someone's annotation.
        use std::io::Read;
        let tmp = NamedTempFile::new().unwrap();
        save_cache(&make_transcripts(), tmp.path()).unwrap();

        let mut zst = zstd::Decoder::new(BufReader::new(File::open(tmp.path()).unwrap())).unwrap();
        let mut magic = [0u8; 8];
        zst.read_exact(&mut magic).unwrap();
        assert_eq!(&magic, CACHE_MAGIC_V3, "published magic drifted from V3");
    }
}
