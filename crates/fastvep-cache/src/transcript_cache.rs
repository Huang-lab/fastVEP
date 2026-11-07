//! Binary transcript cache for fast startup.
//!
//! Serializes fully-built `Vec<Transcript>` (including spliced sequences)
//! to a compact binary format using bincode + zstd compression.
//! Subsequent loads skip GFF3 parsing, FASTA loading, and sequence building.
//!
//! The stream (the magic, the length of the list, the transcripts) is written as zstd frames of
//! about 8 MB, each with a checksum, and a table of them in a skippable frame at the end (the
//! layout of zstd's "seekable format"). A decoder that reads one frame after another decompresses
//! the frames to the same bytes as the one frame it used to be, and skips the table, so the file
//! is read as it always was by anything that takes it for a stream. A run with more than two
//! threads takes the frames from the table and decodes them side by side; any other run, and a
//! cache that an earlier build wrote as one frame, which has no table, is read as a stream with
//! the decoder on a thread of its own. Either way every frame is read to its end, where the
//! decoder checks its checksum: a cache that fails one, or whose stream goes on after its list, is
//! an error.

use anyhow::{Context, Result};
use fastvep_genome::Transcript;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;
use std::time::SystemTime;

/// Magic header for zstd-compressed caches (current format).
///
/// Like V3 to V9 before it, V10 differs from its predecessor in nothing but the
/// guarantee it carries: every V10 cache was written by a build that does not read the CDS flags
/// and the initiator `Met` of an Ensembl export off the codons of an organelle genome with the
/// standard genetic code, as V9's did. See [`load_cache_zstd`] for why that has to be a format bump
/// rather than a note in the changelog.
const CACHE_MAGIC_V10: &[u8; 8] = b"FSTVEP10";
/// Magic header for zstd caches written before the organelle rule (rejected, see
/// [`load_cache_zstd`]).
const CACHE_MAGIC_V9: &[u8; 8] = b"FSTVEP09";
/// Magic header for zstd caches written before the GFF3 phase and exon checks (rejected, see
/// [`load_cache_zstd`]).
const CACHE_MAGIC_V8: &[u8; 8] = b"FSTVEP08";
/// Magic header for zstd caches written before the VEP-cache parity fix (rejected,
/// see [`load_cache_zstd`]).
const CACHE_MAGIC_V7: &[u8; 8] = b"FSTVEP07";
/// Magic header for zstd caches written before the RefSeq naming fix (rejected,
/// see [`load_cache_zstd`]).
const CACHE_MAGIC_V6: &[u8; 8] = b"FSTVEP06";
/// Magic header for zstd caches written before the protein-version fix
/// (rejected, see [`load_cache_zstd`]).
const CACHE_MAGIC_V5: &[u8; 8] = b"FSTVEP05";
/// Magic header for zstd caches written before the APPRIS fix (rejected, see
/// [`load_cache_zstd`]).
const CACHE_MAGIC_V4: &[u8; 8] = b"FSTVEP04";
/// Magic header for zstd caches written before #98 (rejected, see
/// [`load_cache_zstd`]).
const CACHE_MAGIC_V3: &[u8; 8] = b"FSTVEP03";
/// Magic header for zstd caches written before #90 (rejected, see
/// [`load_cache_zstd`]).
const CACHE_MAGIC_V2: &[u8; 8] = b"FSTVEP02";
/// Magic header for legacy gzip-compressed caches (rejected, same reason).
const CACHE_MAGIC_V1: &[u8; 8] = b"FSTVEP01";

/// One line naming what is wrong with a pre-#90 cache.
const SUMMARY_PRE_90: &str = "predates the #90 fix and may hold only part of the \
     transcript set";

/// The evidence behind [`SUMMARY_PRE_90`], for the caller that has to stop.
const DETAIL_PRE_90: &str = "Caches written between 2026-06-10 (the tabix read \
     path) and 2026-08-18 (#90) may hold only the transcripts overlapping one \
     input VCF's variants while looking like a whole-file cache to every later \
     run - the #87 failure, which reported 47,196 of 47,196 chr17 variants as \
     intergenic at exit code 0. Nothing in the file distinguishes the two, so it \
     is rejected rather than read.";

/// One line naming what is wrong with a pre-organelle-rule cache.
const SUMMARY_PRE_ORGANELLE: &str = "predates the organelle rule and may hold CDS flags read off \
     the codons of a mitochondrial or plastid contig with the standard genetic code";

/// The evidence behind [`SUMMARY_PRE_ORGANELLE`], for the caller that has to stop.
const DETAIL_PRE_ORGANELLE: &str = "The reader that wrote it set `cds_start_NF`, `cds_end_NF` and \
     the initiator `Met` of an Ensembl GFF3 export from the codons of every contig but the \
     vertebrate mitochondrial ones it knew by name (`MT`, `chrM`), with the standard genetic \
     code. A yeast `Mito` gene that starts with `ATA`, which that genome's code reads as a start \
     and the standard one as Ile, was flagged `cds_start_NF`, an incomplete transcript, and a \
     change in its first codon was a `missense_variant` and not a `start_lost`; the contigs \
     `mitochondrion_genome`, `Mt` and `Pt` were read the same way. A cache built from a GFF3 \
     with no such contig is unchanged by the rebuild, and nothing in the file says which it \
     is, so it is rejected rather than read.";

/// One line naming what is wrong with a pre-GFF3-checks cache.
const SUMMARY_PRE_GFF_CHECKS: &str = "predates the GFF3 phase and exon checks and may hold \
     a transcript read from a CDS with no phase, or from exons that overlap";

/// The evidence behind [`SUMMARY_PRE_GFF_CHECKS`], for the caller that has to stop.
const DETAIL_PRE_GFF_CHECKS: &str = "The GFF3 parser that wrote it read a CDS phase of `.`, \
     or of anything that is not a number, as -1 and then as the number of bases to pad the \
     first codon with, which is `u64::MAX` of them: a run with a FASTA stopped in \
     `capacity overflow`, one without in an overflow, and a build without overflow checks, \
     which is the release profile, wraps instead. It also kept a transcript whose exons \
     overlap, which stopped a run with a FASTA in `byte range starts at 350 but ends at 200` \
     and gave a run without one cDNA positions for a base that two exons share. Ensembl VEP \
     reads `.` as no phase, stops on any other value, and leaves out a transcript whose exons \
     overlap (or, where they carry a `rank`, keeps it and stops under `--hgvs`). A cache built \
     from a GFF3 with neither is unchanged by the rebuild, and nothing in the file says which \
     it is, so it is rejected rather than read.";

/// One line naming what is wrong with a pre-parity cache.
const SUMMARY_PRE_CACHE_PARITY: &str = "predates the VEP-cache parity fix and \
     names no symbol source, canonical or MANE transcript for a RefSeq gene, no \
     CDS flags for an Ensembl one, and keeps the transcripts of mature miRNAs";

/// The evidence behind [`SUMMARY_PRE_CACHE_PARITY`], for the caller that has to
/// stop.
const DETAIL_PRE_CACHE_PARITY: &str = "Read against Ensembl VEP 116.3's own caches \
     for chr21 and chr22, the GFF3 parser that wrote it differed in what it made \
     of a transcript. For a RefSeq GFF3 it gave a gene no symbol source (the cache \
     writes `EntrezGene`), took the `gbkey` of a non-coding transcript for its \
     biotype (`ncRNA` where the cache has `lncRNA`, `misc_RNA` for \
     `transcribed_pseudogene`), marked no transcript canonical, did not read RefSeq's \
     `tag=MANE Select` (it has a space), kept a transcript for each mature miRNA of \
     a precursor, and did not read `transl_except ... aa:Met` (a GTG start is `Val` \
     where the cache has `Met`). For an Ensembl GFF3 it gave a gene neither the \
     symbol source nor the HGNC ID its description names, and wrote no \
     `cds_start_NF` or `cds_end_NF` (a GFF3 export has none, and the cache's are \
     read from the sequence once a run has a FASTA; a cache keeps the sequences it \
     was built with) or `initial_met`. A cache built from a GFF3 whose genes have \
     none of these is unchanged by the rebuild; nothing in the file says which \
     source it came from, so it is rejected rather than read.";

/// One line naming what is wrong with a pre-RefSeq-naming cache.
const SUMMARY_PRE_REFSEQ_NAMES: &str = "predates the RefSeq naming fix and names \
     RefSeq transcripts, genes and biotypes as the GFF3's ids and SO terms";

/// The evidence behind [`SUMMARY_PRE_REFSEQ_NAMES`], for the caller that has to
/// stop.
const DETAIL_PRE_REFSEQ_NAMES: &str = "The GFF3 parser that wrote it named a \
     transcript by its ID and a gene by its own ID, and took a transcript's \
     biotype from its SO term when the record states none. For a RefSeq GFF3 that \
     is `rna-NM_000001.3` where Ensembl VEP writes `NM_000001.3`, `gene-TP53` \
     where it writes the GeneID `7157`, and `mRNA` or `lnc_RNA` where it writes \
     `protein_coding` or `ncRNA`, in the Feature, Gene and BIOTYPE of every row \
     and the transcript of every HGVSc. An Ensembl GFF3 is named the same by \
     either rule, so a cache built from one is unchanged by the rebuild; nothing \
     in the file says which source it came from, so it is rejected rather than \
     read.";

/// One line naming what is wrong with a pre-protein-version cache.
const SUMMARY_PRE_PROTEIN_VERSION: &str = "predates the protein-version fix and \
     names every protein version 1";

/// The evidence behind [`SUMMARY_PRE_PROTEIN_VERSION`], for the caller that has to
/// stop.
const DETAIL_PRE_PROTEIN_VERSION: &str = "The GFF3 parser that wrote it parsed the \
     `version` of each CDS record and then gave every translation version 1, so \
     every HGVSp the cache produces names `ENSP….1` for a protein that is \
     `ENSP….3` or `ENSP….12`, which is one protein in six of Ensembl 115's. \
     Ensembl VEP 116 reads the version (and 115 did not), so the two disagree on \
     the accession of those rows and on nothing else. A cache built from a GFF3 \
     whose proteins are all version 1, or that has no `version` attribute (NCBI's), \
     is unchanged by the rebuild; nothing in the file says which source it came \
     from, so it is rejected rather than read.";

/// One line naming what is wrong with a pre-APPRIS cache.
const SUMMARY_PRE_APPRIS: &str = "predates the APPRIS fix and carries no APPRIS \
     tier for any transcript";

/// The evidence behind [`SUMMARY_PRE_APPRIS`], for the caller that has to stop.
const DETAIL_PRE_APPRIS: &str = "The GFF3 parser that wrote it dropped the \
     `appris_principal_*` / `appris_alternative_*` tags GENCODE carries in the \
     same `tag=` list as MANE_Select, so every transcript in the file has an \
     empty APPRIS column and `--pick` reaches the APPRIS tier to find every \
     candidate equal - at a locus where two genes' canonical transcripts \
     overlap and neither is MANE, the pick then falls through to TSL and \
     finally to alphabetical transcript ID. Ensembl's own GFF3 has no APPRIS \
     tag to lose, so a cache built from one is unchanged by the rebuild; \
     nothing in the file says which source it came from, so it is rejected \
     rather than read.";

/// One line naming what is wrong with a pre-#98 cache.
const SUMMARY_PRE_98: &str = "predates the #98 fix and holds unnamed non-coding \
     genes";

/// The evidence behind [`SUMMARY_PRE_98`], for the caller that has to stop.
const DETAIL_PRE_98: &str = "The GFF3 parser that wrote it took only `gene` and \
     `pseudogene` as gene records, so every locus Ensembl types `ncRNA_gene` - \
     the lncRNA, miRNA and snRNA genes - was cached with no symbol and biotype \
     \"unknown\". SYMBOL comes out empty where VEP prints a name, and every \
     symbol-keyed annotation (OMIM, gnomAD constraint, ClinVar protein) silently \
     misses those genes. The file is intact, and nothing in it separates a gene \
     that lost its name from one that never had one, so it is rejected rather \
     than read.";

/// A cache whose format carries a guarantee this build requires and the file
/// cannot supply: it may be region-restricted (pre-#90) or hold unnamed
/// non-coding genes (pre-#98).
///
/// A distinct type rather than a bare message because the two callers need to
/// say different things about it, at different lengths.
///
/// The sidecar path rebuilds from the GFF3 on its own, so nothing is being
/// asked of the user and one line is the whole story - [`Self::summary`]. The
/// explicit `--transcript-cache` path has to stop, so it owes the user the
/// evidence and a way out; it uses [`Display`], which adds [`Self::detail`].
/// Neither carries recovery advice, because the two paths recover differently
/// and a fixed clause was wrong on both: the sidecar had already rebuilt by the
/// time it printed "let the sidecar rebuild itself", and deleting a cache named
/// by `--transcript-cache` conjures no sidecar to replace it.
///
/// The generic read-failure wording at the explicit call site ("most likely
/// truncated or corrupt") is also the wrong diagnosis here - the file is
/// intact, it is the guarantee that is missing. Callers tell the two apart with
/// `anyhow::Error::downcast_ref`.
///
/// [`Display`]: std::fmt::Display
#[derive(Debug)]
pub struct StaleCacheFormat {
    /// The format found on disk, named as the user would see it in the file.
    pub found: &'static str,
    /// One clause completing "the {found} format, which ...". Per-format
    /// because the formats are stale for different reasons: V1/V2 may be
    /// region-restricted (#90), V3 is whole-file but has unnamed non-coding
    /// genes (#98). Telling a user the wrong one sends them looking for the
    /// wrong problem.
    pub summary: &'static str,
    /// The evidence behind `summary`, for the caller that has to stop. Kept out
    /// of the sidecar warning, which recovers by itself and would otherwise
    /// spend a paragraph on a problem the user does not have to act on.
    pub detail: &'static str,
}

impl StaleCacheFormat {
    /// The one-line form, for a caller that recovers on its own.
    ///
    /// A bare predicate with no subject and no trailing stop, so a call site can
    /// put the file in front of it - "sidecar cache {path} {summary}; rebuilding
    /// from GFF3" - and punctuate to suit.
    pub fn summary(&self) -> String {
        format!("is in the {} format, which {}", self.found, self.summary)
    }
}

impl std::fmt::Display for StaleCacheFormat {
    /// The long form: the summary plus the evidence for it. Supplies the
    /// pronoun subject that [`Self::summary`] leaves out, so it reads after a
    /// caller's "cannot be used:" and can be followed by whatever recovery
    /// applies on that path.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "it {}. {}", self.summary(), self.detail)
    }
}

impl std::error::Error for StaleCacheFormat {}

/// About how many bytes of the decompressed stream go in one frame of a cache.
///
/// Big enough that zstd, which starts each frame with an empty window, loses next to nothing of its
/// ratio, and small enough that a cache of the transcripts of Ensembl's GFF3, about 1.6 GB
/// decompressed, is some two hundred frames for the cores of a machine to share.
const FRAME_BYTES: usize = 8 << 20;

/// The magic of the skippable frame that holds the table of frames, and the last four bytes of the
/// file, which say that a table ends it: both from the seekable format of zstd
/// (`contrib/seekable_format` in the zstd repository), so that anything that reads that finds the
/// frames in it too.
const FRAME_TABLE_MAGIC: u32 = 0x184D_2A5E;
const FRAME_TABLE_FOOTER_MAGIC: u32 = 0x8F92_EAB1;

/// Counts the bytes that go through it, so that the length of a frame is known when it is done.
struct Counted<W> {
    inner: W,
    written: u64,
}

impl<W: std::io::Write> std::io::Write for Counted<W> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let n = self.inner.write(buf)?;
        self.written += n as u64;
        Ok(n)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

/// A writer of zstd frames, each with its own checksum, and then the table of them. The caller says
/// where a frame may end, between two transcripts, so that a frame holds whole ones.
struct Frames<W: std::io::Write> {
    encoder: Option<zstd::Encoder<'static, Counted<W>>>,
    /// The bytes of the frame that is open, and where in the output it began.
    in_frame: usize,
    began_at: u64,
    /// The bytes of the stream so far, which tells how long a transcript was once it is written.
    total: u64,
    /// How many bytes of the stream go in a frame before the next transcript begins another.
    frame_bytes: usize,
    /// The length of each frame that is done, as it is in the file and as it is decompressed.
    done: Vec<(u64, u64)>,
}

impl<W: std::io::Write> Frames<W> {
    fn new(writer: W, frame_bytes: usize) -> std::io::Result<Self> {
        let writer = Counted {
            inner: writer,
            written: 0,
        };
        Ok(Self {
            encoder: Some(Self::encoder(writer)?),
            in_frame: 0,
            began_at: 0,
            total: 0,
            frame_bytes,
            done: Vec::new(),
        })
    }

    fn encoder(writer: Counted<W>) -> std::io::Result<zstd::Encoder<'static, Counted<W>>> {
        // zstd level 1: fast compression, still much better decompression than gzip.
        let mut encoder = zstd::Encoder::new(writer, 1)?;
        // Record a frame checksum so a cache corrupted in place - a bad sector, a
        // half-flushed page - fails loudly at load instead of deserializing into
        // something subtly wrong.
        encoder.include_checksum(true)?;
        Ok(encoder)
    }

    /// Ends the frame that is open and notes how long it is.
    fn end_the_frame(&mut self) -> std::io::Result<Counted<W>> {
        let writer = self.encoder.take().expect("a frame is open").finish()?;
        self.done
            .push((writer.written - self.began_at, self.in_frame as u64));
        Ok(writer)
    }

    /// Ends the frame, and begins another, if this one has `frame_bytes` in it. Called between two
    /// transcripts, so that a frame holds whole ones.
    fn begin_a_frame_if_this_one_is_full(&mut self) -> std::io::Result<()> {
        if self.in_frame >= self.frame_bytes {
            let writer = self.end_the_frame()?;
            self.began_at = writer.written;
            self.encoder = Some(Self::encoder(writer)?);
            self.in_frame = 0;
        }
        Ok(())
    }

    /// Ends the last frame, and writes the table of them: a skippable frame, which a decoder that
    /// reads the file as a stream passes over, with for each frame its length in the file and
    /// decompressed (four bytes each, as the seekable format has them), and then how many frames
    /// there are, a byte of flags (none: no checksum of its own for an entry), and the magic.
    fn finish(mut self) -> std::io::Result<W> {
        let mut writer = self.end_the_frame()?;
        let too_long = |_| std::io::Error::other("a frame of the cache is longer than 4 GiB");
        let too_many = |_| std::io::Error::other("a cache of more than 500 million frames");
        let entries = u32::try_from(self.done.len()).map_err(too_many)?;
        let size = u32::try_from(8 * u64::from(entries) + 9).map_err(too_many)?;
        let mut table = Vec::with_capacity(8 * self.done.len() + 17);
        table.extend_from_slice(&FRAME_TABLE_MAGIC.to_le_bytes());
        table.extend_from_slice(&size.to_le_bytes());
        for &(in_file, decompressed) in &self.done {
            table.extend_from_slice(&u32::try_from(in_file).map_err(too_long)?.to_le_bytes());
            table.extend_from_slice(&u32::try_from(decompressed).map_err(too_long)?.to_le_bytes());
        }
        table.extend_from_slice(&entries.to_le_bytes());
        table.push(0);
        table.extend_from_slice(&FRAME_TABLE_FOOTER_MAGIC.to_le_bytes());
        std::io::Write::write_all(&mut writer, &table)?;
        Ok(writer.inner)
    }
}

impl<W: std::io::Write> std::io::Write for Frames<W> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let n = self.encoder.as_mut().expect("a frame is open").write(buf)?;
        self.in_frame += n;
        self.total += n as u64;
        Ok(n)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.encoder.as_mut().expect("a frame is open").flush()
    }
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
    save_cache_with(transcripts, path, FRAME_BYTES)
}

/// [`save_cache`] with the size of a frame as a parameter, which a test makes small, so that a
/// handful of transcripts are several frames.
fn save_cache_with(transcripts: &[Transcript], path: &Path, frame_bytes: usize) -> Result<()> {
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
    use std::io::Write;
    // The stream is the magic, the length of the list and then the transcripts, one after the
    // other: the bytes `bincode::serialize_into(&mut zst, transcripts)` gave after the magic when
    // the whole of it was one frame. It is cut into frames between two transcripts, each of about
    // `frame_bytes`, so that a reader can decode the frames side by side, and the table of them
    // ends the file (`Frames::finish`). A reader that takes the file as one stream of frames, as
    // every reader before this one did, reads the same bytes and passes over the table.
    let mut frames = Frames::new(writer, frame_bytes)?;
    frames.write_all(CACHE_MAGIC_V10)?;
    bincode::serialize_into(&mut frames, &(transcripts.len() as u64))
        .with_context(|| "Serializing transcripts to cache")?;
    for transcript in transcripts {
        frames.begin_a_frame_if_this_one_is_full()?;
        let before = frames.total;
        bincode::serialize_into(&mut frames, transcript)
            .with_context(|| "Serializing transcripts to cache")?;
        // A cache is read with a limit on the bytes of one transcript (`bounded`), and a cache with
        // one that is over it would be refused at every load, and written again by every run that
        // rebuilds it.
        let size = frames.total - before;
        anyhow::ensure!(
            size <= LARGEST_VALUE,
            "the transcript {} is {size} bytes in a cache, and a cache is read with a limit of \
             {LARGEST_VALUE} bytes for one",
            transcript.stable_id
        );
    }
    let mut writer = frames.finish()?;
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
    use std::io::{Read, Seek};
    let mut file =
        File::open(path).with_context(|| format!("Opening cache file: {}", path.display()))?;

    // Peek at the first bytes to detect format.
    // zstd frames start with 0x28B52FFD; gzip starts with 0x1F8B.
    let mut peek = [0u8; 4];
    file.read_exact(&mut peek)
        .with_context(|| "Reading cache header")?;

    // Rewind so the decompressor sees the full stream
    file.seek(std::io::SeekFrom::Start(0))?;

    if peek[0..2] == [0x1F, 0x8B] {
        // Legacy gzip format (v1)
        load_cache_gzip(BufReader::new(file))
    } else {
        // zstd format (v2, or future)
        load_cache_zstd(file)
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
///
/// `FSTVEP09` is rejected on the same grounds: the reader that wrote it flagged the genes of a
/// yeast `Mito` contig that start with `ATA` `cds_start_NF`, and the sidecar path reuses a cache
/// whenever it is newer than its GFF3, so without a bump an upgrade leaves the flag in every row of
/// them. A V9 cache loads cleanly; it just has the flag, which is indistinguishable from an
/// annotation that really called the transcript incomplete.
///
/// `FSTVEP08` is rejected on the same grounds: the parser that wrote it read a CDS phase of `.`
/// as `u64::MAX` bases of padding and kept a transcript whose exons overlap, and the sidecar
/// path reuses a cache whenever it is newer than its GFF3, so without a bump an upgrade leaves
/// both in the cache that a run then annotates from. A V8 cache loads cleanly; it just holds
/// them, which is indistinguishable from a GFF3 that really has such a CDS.
///
/// `FSTVEP07` is rejected on the same grounds: the parser that wrote it left the symbol
/// source of a RefSeq gene empty, called a long non-coding RNA `ncRNA`, marked no RefSeq
/// transcript canonical and read no `MANE Select`, and (with a FASTA) wrote no
/// `cds_start_NF` or `cds_end_NF` for an Ensembl transcript. The sidecar path reuses a cache
/// whenever it is newer than its GFF3, so without a bump an upgrade leaves every one of them in
/// every row. A V7 cache loads cleanly; it just has none of these, which is indistinguishable from
/// a GFF3 that never had them.
///
/// `FSTVEP06` is rejected on the same grounds: the parser that wrote it named a RefSeq
/// transcript `rna-NM_000001.3`, its gene `gene-TP53` and its biotype `mRNA`, and the
/// sidecar path reuses a cache whenever it is newer than its GFF3, so without a bump
/// an upgrade leaves those names in every row. A V6 cache loads cleanly; it just has
/// the old names, which are indistinguishable from a GFF3 whose ids really are those.
///
/// `FSTVEP05` is rejected on the same grounds: the parser that wrote it gave every
/// protein version 1, and the sidecar path reuses a cache whenever it is newer than
/// its GFF3, so without a bump an upgrade would leave `ENSP….1` in front of every
/// `HGVSp` and the user would re-run the identical command to the identical
/// answer. A V5 cache loads cleanly; it just has no protein version but 1, which
/// is indistinguishable from a GFF3 whose proteins really are version 1.
///
/// `FSTVEP04` is rejected on the same grounds as `FSTVEP03`: the parser that
/// wrote it could not see GENCODE's `appris_*` tags, and the sidecar path
/// reuses a cache whenever it is newer than its GFF3, so without a bump the
/// APPRIS column stays empty and `--pick`'s APPRIS tier stays inert after the
/// upgrade. A V4 cache loads cleanly; it just has no APPRIS call for any
/// transcript, which is indistinguishable from an Ensembl source that never
/// carried one.
///
/// `FSTVEP03` is rejected for the same reason at one remove. #98 fixed the GFF3
/// parser, not the caches it had already written, and the sidecar path reuses a
/// cache whenever it is newer than its GFF3 - so without a bump the user who
/// reported #98 would upgrade, re-run the identical command, and still get an
/// empty SYMBOL, because the parser never runs again. A V3 cache is whole-file
/// and loads cleanly; it just has no name for any `ncRNA_gene` locus, which is
/// indistinguishable from a genuinely unnamed gene.
fn load_cache_zstd(mut file: File) -> Result<Vec<Transcript>> {
    use std::io::Seek;
    // A cache of several frames has them decoded side by side, when there are cores to do it on.
    // The frames are found from the table at the end of the file, without reading any of them. The
    // table is only where they are: the stream has every one, each with its checksum, and a run
    // that reads the file as one does not look at the table. Where it is used, a table that does
    // not describe its file is a cache that is damaged, and is refused as one.
    if let Some(workers) = side_by_side_workers(rayon::current_num_threads()) {
        if let Some(frames) = read_frame_table(&file)? {
            if frames.len() > 1 {
                return decode_frames_side_by_side(&file, &frames, workers)
                    .with_context(|| "Deserializing transcripts from cache");
            }
        }
    }

    // A cache of one frame, which is what every build before the one that cut them wrote, or a run
    // with no cores to spare, is read as a stream: the decompression goes ahead of the
    // deserialization on a thread of its own. (A positioned read leaves the cursor of the file
    // where it was on Unix, and does not on Windows.)
    file.seek(std::io::SeekFrom::Start(0))?;
    let mut zst = zstd::Decoder::new(file)?;
    check_magic(&read_magic(&mut zst)?)?;
    deserialize_decoded(zst).with_context(|| "Deserializing transcripts from cache")
}

/// The transcripts of the frames of a cache, in order, with the frames decoded side by side.
///
/// A frame is appended to the list as soon as it and every frame before it are decoded, and let
/// go: the list is the one copy of the transcripts that is held, and what is read of the file at
/// any moment is the few frames that are being decoded, and not the file.
fn decode_frames_side_by_side(
    file: &File,
    frames: &[FrameEntry],
    workers: usize,
) -> Result<Vec<Transcript>> {
    let mut transcripts: Vec<Transcript> = Vec::new();
    let mut count: Option<u64> = None;
    decode_in_order(
        frames,
        workers,
        |i, frame| {
            read_frame(file, frame)
                .and_then(|bytes| decode_frame(&bytes, i == 0, frame.decompressed))
        },
        |(mut frame, header): (Vec<Transcript>, Option<u64>)| {
            if let Some(n) = header {
                count = Some(n);
                // a hint, and not a promise: the length is the cache's own word, so it may be refused
                let _ = transcripts.try_reserve_exact(usize::try_from(n).unwrap_or(0).min(1 << 24));
            }
            transcripts.append(&mut frame);
            Ok(())
        },
    )?;
    let count = count.expect("the first frame holds the length of the list");
    anyhow::ensure!(
        transcripts.len() as u64 == count,
        "the cache says it holds {count} transcripts and its {} frames hold {}",
        frames.len(),
        transcripts.len()
    );
    Ok(transcripts)
}

/// The most frames decoded at once. A decoder has some 10 MB in flight, the frame and the
/// transcripts built from it until they are appended to the list, so a decoder for every thread of a
/// large machine raised the peak of a run with its cores: 3.5 GB on 64 threads, against 2.7 GB on 8
/// and on 64 with this, in the same time.
const MOST_DECODERS: usize = 8;

/// How many frames to decode at once in a pool of `threads`: one thread of them merges, and the rest
/// decode, up to `MOST_DECODERS` of them.
fn decoders_on(threads: usize) -> usize {
    threads.saturating_sub(1).clamp(1, MOST_DECODERS)
}

/// How many frames a run of `threads` decodes at once, and `None` for a run that reads the cache as
/// a stream: two threads or fewer have no core to spare for it, and a platform that has neither
/// `read_exact_at` nor `seek_read` has no frame to read apart from the others.
fn side_by_side_workers(threads: usize) -> Option<usize> {
    (cfg!(any(unix, windows)) && threads > 2).then(|| decoders_on(threads))
}

/// Runs `decode` on every item, on `workers` threads at most at once, and gives each result to
/// `merge` in the order of the items, as soon as it and every item before it are done. The first
/// error, from a decode or from the merge, is returned, and the items not yet begun are left alone.
///
/// A worker takes the next item that nobody has begun, so that the items are begun in order and a
/// slow thread, an efficiency core of a machine that has both, holds up no more than the item it is
/// on. Items of about the same cost finish in about the order they were begun, and few of them wait
/// for one before them; a slow one lets the others run ahead of it, and what they decoded waits in
/// memory, as pieces of the list, until it is done.
///
/// The workers are threads of their own, and the calling thread merges. They are not tasks of the
/// rayon pool: a call from a task of the pool, which holds a thread of it while it waits for the
/// results, would need the others to be free, and two such calls at once in a pool of two would
/// wait for each other for ever.
fn decode_in_order<I, T>(
    items: &[I],
    workers: usize,
    decode: impl Fn(usize, &I) -> Result<T> + Sync,
    merge: impl FnMut(T) -> Result<()>,
) -> Result<()>
where
    I: Sync,
    T: Send,
{
    decode_in_order_with(items, workers, || Ok(()), decode, merge)
}

/// Sets the flag that stops the items that are not yet begun if the thread that holds it unwinds: a
/// panic in a decode or in the merge would otherwise be reported once every item had been decoded.
struct StopOnPanic<'a>(&'a std::sync::atomic::AtomicBool);

impl Drop for StopOnPanic<'_> {
    fn drop(&mut self) {
        if std::thread::panicking() {
            self.0.store(true, std::sync::atomic::Ordering::Relaxed);
        }
    }
}

/// For a test to tell how many workers a load had, and so that it was side by side at all: the
/// calling thread, which is the one a test knows, and the workers that were given a thread.
#[cfg(test)]
static WORKERS_OF_THE_LOADS: std::sync::Mutex<Vec<(std::thread::ThreadId, usize)>> =
    std::sync::Mutex::new(Vec::new());

/// [`decode_in_order`] with a request for each worker's thread that comes before the thread is
/// spawned, and that a test makes fail: a host that has reached its limit on threads refuses the
/// spawn. The workers that were given a thread do all of it, and with none the items are decoded
/// and merged one after the other, on the calling thread, as a run with no cores to spare would.
fn decode_in_order_with<I, T>(
    items: &[I],
    workers: usize,
    ask_for_a_thread: impl Fn() -> std::io::Result<()>,
    decode: impl Fn(usize, &I) -> Result<T> + Sync,
    mut merge: impl FnMut(T) -> Result<()>,
) -> Result<()>
where
    I: Sync,
    T: Send,
{
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
    let (sender, decoded) = std::sync::mpsc::channel::<(usize, Result<T>)>();
    // Set by the first failure, so that the items not yet begun are left alone.
    let failed = AtomicBool::new(false);
    // The next item that nobody has begun.
    let unbegun = AtomicUsize::new(0);
    // the results of items that were decoded before an item that comes earlier was
    let mut waiting: std::collections::BTreeMap<usize, T> = Default::default();
    let mut next = 0;
    std::thread::scope(|scope| -> Result<()> {
        let (failed, unbegun, decode) = (&failed, &unbegun, &decode);
        let _stop_on_panic = StopOnPanic(failed);
        let mut started = 0;
        for _ in 0..workers.max(1) {
            let sender = sender.clone();
            let spawned = ask_for_a_thread().and_then(|()| {
                std::thread::Builder::new()
                    .name("cache-frames".into())
                    .spawn_scoped(scope, move || {
                        let _stop_on_panic = StopOnPanic(failed);
                        loop {
                            if failed.load(Ordering::Relaxed) {
                                break;
                            }
                            let i = unbegun.fetch_add(1, Ordering::Relaxed);
                            if i >= items.len() {
                                break;
                            }
                            let _ = sender.send((i, decode(i, &items[i])));
                        }
                    })
            });
            if spawned.is_err() {
                break;
            }
            started += 1;
        }
        #[cfg(test)]
        WORKERS_OF_THE_LOADS
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .push((std::thread::current().id(), started));
        // the loop below ends when the last of the clones is gone, and not before
        drop(sender);
        if started == 0 {
            for (i, item) in items.iter().enumerate() {
                merge(decode(i, item)?)?;
                next += 1;
            }
            return Ok(());
        }
        for (i, outcome) in decoded {
            match outcome {
                Ok(item) => {
                    waiting.insert(i, item);
                }
                Err(e) => {
                    failed.store(true, Ordering::Relaxed);
                    return Err(e);
                }
            }
            while let Some(item) = waiting.remove(&next) {
                if let Err(e) = merge(item) {
                    failed.store(true, Ordering::Relaxed);
                    return Err(e);
                }
                next += 1;
            }
        }
        Ok(())
    })?;
    anyhow::ensure!(
        next == items.len(),
        "{next} of {} items were decoded",
        items.len()
    );
    Ok(())
}

/// The first eight bytes of the stream of a cache.
fn read_magic(stream: &mut impl std::io::Read) -> Result<[u8; 8]> {
    let mut magic = [0u8; 8];
    stream
        .read_exact(&mut magic)
        .with_context(|| "Reading cache header")?;
    Ok(magic)
}

/// What a cache has to say of its format before anything is read of the rest of it: the formats that
/// were written by builds that are stale, with why, and nothing else but the current one.
fn check_magic(magic: &[u8; 8]) -> Result<()> {
    if magic == CACHE_MAGIC_V2 {
        return Err(StaleCacheFormat {
            found: "FSTVEP02",
            summary: SUMMARY_PRE_90,
            detail: DETAIL_PRE_90,
        }
        .into());
    }
    if magic == CACHE_MAGIC_V3 {
        return Err(StaleCacheFormat {
            found: "FSTVEP03",
            summary: SUMMARY_PRE_98,
            detail: DETAIL_PRE_98,
        }
        .into());
    }
    if magic == CACHE_MAGIC_V4 {
        return Err(StaleCacheFormat {
            found: "FSTVEP04",
            summary: SUMMARY_PRE_APPRIS,
            detail: DETAIL_PRE_APPRIS,
        }
        .into());
    }
    if magic == CACHE_MAGIC_V5 {
        return Err(StaleCacheFormat {
            found: "FSTVEP05",
            summary: SUMMARY_PRE_PROTEIN_VERSION,
            detail: DETAIL_PRE_PROTEIN_VERSION,
        }
        .into());
    }
    if magic == CACHE_MAGIC_V6 {
        return Err(StaleCacheFormat {
            found: "FSTVEP06",
            summary: SUMMARY_PRE_REFSEQ_NAMES,
            detail: DETAIL_PRE_REFSEQ_NAMES,
        }
        .into());
    }
    if magic == CACHE_MAGIC_V7 {
        return Err(StaleCacheFormat {
            found: "FSTVEP07",
            summary: SUMMARY_PRE_CACHE_PARITY,
            detail: DETAIL_PRE_CACHE_PARITY,
        }
        .into());
    }
    if magic == CACHE_MAGIC_V8 {
        return Err(StaleCacheFormat {
            found: "FSTVEP08",
            summary: SUMMARY_PRE_GFF_CHECKS,
            detail: DETAIL_PRE_GFF_CHECKS,
        }
        .into());
    }
    if magic == CACHE_MAGIC_V9 {
        return Err(StaleCacheFormat {
            found: "FSTVEP09",
            summary: SUMMARY_PRE_ORGANELLE,
            detail: DETAIL_PRE_ORGANELLE,
        }
        .into());
    }
    if magic != CACHE_MAGIC_V10 {
        anyhow::bail!("Invalid cache file (wrong magic header, expected FSTVEP10)");
    }

    Ok(())
}

/// A frame of a cache as its table has it: where it begins in the file, how long it is there, and how
/// long it is decompressed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct FrameEntry {
    at: u64,
    len: u64,
    decompressed: u64,
}

/// Where each frame of a cache begins, and how long it is, from the table that ends the file that
/// `save_cache` wrote: the frames one after the other, and then the table, with nothing left over.
/// `None` for a file with no table, which is every cache that an earlier build wrote, and is read
/// as a stream.
///
/// A file that ends as a table does, and whose table does not describe it, is an error here: it is not
/// what `save_cache` wrote. It is only here that the table is read. A run that reads the file as a
/// stream (one or two threads) has the frames without it, each with its own checksum, and loads what a
/// table that was right would have found, so the same file loads on two threads and is refused on
/// three: the sidecar is then rebuilt, with a table that describes it, and a cache named on the
/// command line stops the run. What the table says is not taken on trust for what is in the frames
/// either: each is decoded to its end, where its own checksum is checked, to the length that its entry
/// says, and the transcripts they hold have to be as many as the cache says. One cache in 2^32 of
/// those that earlier builds wrote, one whose last four bytes (its checksum) are those of the footer,
/// is refused here for want of a table, and rebuilt.
fn read_frame_table(file: &File) -> Result<Option<Vec<FrameEntry>>> {
    let file_len = file.metadata()?.len();
    let mut footer = [0u8; 9];
    if file_len < footer.len() as u64 {
        return Ok(None);
    }
    read_at(file, &mut footer, file_len - 9).with_context(|| "Reading the end of the cache")?;
    if u32::from_le_bytes([footer[5], footer[6], footer[7], footer[8]]) != FRAME_TABLE_FOOTER_MAGIC
    {
        return Ok(None);
    }
    let entries = u64::from(u32::from_le_bytes([
        footer[0], footer[1], footer[2], footer[3],
    ]));
    anyhow::ensure!(
        footer[4] == 0,
        "the table of frames at the end of the cache has flags that this does not read ({:#04x})",
        footer[4]
    );
    // the header of the skippable frame, the entries, and the footer
    let table_len = 8 + 8 * entries + 9;
    anyhow::ensure!(
        entries > 0 && table_len <= file_len,
        "the table of frames at the end of the cache says {entries} frames, in a file of {file_len} bytes"
    );
    let table_at = file_len - table_len;
    let mut table = vec![0u8; usize::try_from(table_len - 9)?];
    read_at(file, &mut table, table_at).with_context(|| "Reading the table of frames")?;
    let word =
        |at: usize| u32::from_le_bytes([table[at], table[at + 1], table[at + 2], table[at + 3]]);
    anyhow::ensure!(
        word(0) == FRAME_TABLE_MAGIC && u64::from(word(4)) == 8 * entries + 9,
        "the table of frames at the end of the cache is not in a skippable frame of its size"
    );
    let mut frames = Vec::with_capacity(usize::try_from(entries)?);
    let mut at = 0u64;
    for entry in 0..entries as usize {
        let len = u64::from(word(8 + 8 * entry));
        anyhow::ensure!(len > 0, "the table of frames has a frame of no bytes");
        frames.push(FrameEntry {
            at,
            len,
            decompressed: u64::from(word(12 + 8 * entry)),
        });
        at += len;
    }
    anyhow::ensure!(
        at == table_at,
        "the frames in the table of the cache are {at} bytes, and the file holds {table_at} before it"
    );
    Ok(Some(frames))
}

/// `buf.len()` bytes of `file` from `offset`, on the one handle, which is what every thread that
/// decodes a frame reads through: a file that is replaced while the cache is read is the file it was
/// when it was opened, to all of them.
#[cfg(unix)]
fn read_at(file: &File, buf: &mut [u8], offset: u64) -> std::io::Result<()> {
    use std::os::unix::fs::FileExt;
    file.read_exact_at(buf, offset)
}

#[cfg(windows)]
fn read_at(file: &File, mut buf: &mut [u8], mut offset: u64) -> std::io::Result<()> {
    use std::os::windows::fs::FileExt;
    while !buf.is_empty() {
        match file.seek_read(buf, offset)? {
            0 => return Err(std::io::ErrorKind::UnexpectedEof.into()),
            n => {
                buf = &mut buf[n..];
                offset += n as u64;
            }
        }
    }
    Ok(())
}

/// A platform with neither has no frame to read apart from the others: the cache is one stream.
#[cfg(not(any(unix, windows)))]
fn read_at(_: &File, _: &mut [u8], _: u64) -> std::io::Result<()> {
    Err(std::io::ErrorKind::Unsupported.into())
}

fn read_frame(file: &File, frame: &FrameEntry) -> Result<Vec<u8>> {
    let mut bytes = vec![0u8; usize::try_from(frame.len)?];
    read_at(file, &mut bytes, frame.at)
        .with_context(|| format!("Reading the frame at byte {}", frame.at))?;
    Ok(bytes)
}

/// Counts the bytes that are read through it.
struct CountedRead<R> {
    inner: R,
    read: u64,
}

impl<R: std::io::Read> std::io::Read for CountedRead<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.read += n as u64;
        Ok(n)
    }
}

/// The transcripts of one frame of a cache of several, and, for the first, the length of the list
/// that the frames hold between them.
///
/// The frame is read as a stream, to its end, which is where the decoder checks its checksum, so a
/// frame that is damaged is an error and none of what was read of it is used. Until then what is
/// read is not known to be what was written, which is why every call of `bincode` has the limit of
/// [`bounded`]. The first frame has the magic, which says whether the cache is one that this reads
/// and why not if it is not, before anything is made of the rest. The frame is the bytes of its entry
/// in the table and no others, and decompresses to as many as the table says: a table that were
/// wrong about a frame would leave some of its bytes unread here, or fail to decode them, or make
/// more or fewer of them.
fn decode_frame(
    frame: &[u8],
    first: bool,
    decompressed: u64,
) -> Result<(Vec<Transcript>, Option<u64>)> {
    use bincode::Options;
    use std::io::BufRead;
    let mut unread = frame;
    let mut stream = BufReader::with_capacity(
        DECODE_CHUNK / 4,
        CountedRead {
            inner: zstd::Decoder::with_buffer(&mut unread)?.single_frame(),
            read: 0,
        },
    );
    let count = if first {
        check_magic(&read_magic(&mut stream)?)?;
        Some(bounded().deserialize_from::<_, u64>(&mut stream)?)
    } else {
        None
    };
    let mut transcripts = Vec::new();
    while !stream.fill_buf()?.is_empty() {
        transcripts.push(bounded().deserialize_from::<_, Transcript>(&mut stream)?);
    }
    let read = stream.get_ref().read;
    drop(stream);
    anyhow::ensure!(
        unread.is_empty(),
        "a frame of the cache is followed, inside its entry in the table, by {} more bytes",
        unread.len()
    );
    anyhow::ensure!(
        read == decompressed,
        "a frame of the cache decompresses to {read} bytes, and the table of frames says {decompressed}"
    );
    Ok((transcripts, count))
}

/// How much of the decompressed stream goes from the thread that decodes it to the thread that
/// deserializes it in one piece, and how many pieces the decoder may be ahead by.
const DECODE_CHUNK: usize = 1 << 20;
const DECODE_AHEAD: usize = 4;

/// Deserialize the rest of a cache whose header has been read, decompressing it on a thread of its own.
///
/// A cache of the 500,000 transcripts of Ensembl's GFF3 is a second or two of one thread, about half of
/// it zstd and half of it `bincode` building the transcripts, and the two were one after the other
/// in the thread that wanted the answer: `bincode::deserialize_from` reads a field at a time
/// and each read went into the decoder. The decoder now runs ahead of it, in pieces of a megabyte, and
/// the wait for the cache is the longer of the two and not the sum. Every cache reads this way
/// when the run has two threads or fewer, and a cache of one frame, which is what every build before
/// the frames wrote, whatever the threads.
///
/// The stream is then read to its end, and has to end where the list does. The end is where the
/// decoder checks the checksum `save_cache` records in the frame "so a cache corrupted in place
/// fails loudly at load instead of deserializing into something subtly wrong", and a decoder that
/// is not asked for another byte does not reliably get there: the deserializer stops at the last
/// byte of the last transcript. With the checksum flipped, or a byte in the middle of the cache, or
/// bytes after the frame, the loader of v0.4.0 read caches of several sizes as the transcripts they
/// had become. What follows the list is not nothing either: two caches joined into one file by `cat`
/// loaded as the first alone, with nothing to say that the transcripts of the second had been
/// left out.
fn deserialize_decoded<R: std::io::Read + Send>(decoder: R) -> Result<Vec<Transcript>> {
    deserialize_decoded_with(decoder, || Ok(()))
}

/// [`deserialize_decoded`] with a request for the decoder's thread that comes before the thread is
/// spawned, and that a test makes fail: a host that has reached its limit on threads, as a
/// container with a low limit on its processes has, refuses the spawn and the load has to be made
/// all the same. A stack that no system can give is not a refusal that every system returns: where
/// a kernel answers an error, Rosetta (x86_64 images on Apple silicon) aborts the process with
/// "could not find free space for allocation" for a stack of 1 << 50 bytes, and the test that asked
/// for one ended the run of every test after it.
fn deserialize_decoded_with<R: std::io::Read + Send>(
    decoder: R,
    ask_for_a_thread: impl FnOnce() -> std::io::Result<()>,
) -> Result<Vec<Transcript>> {
    use std::sync::mpsc;
    let (pieces, decoded) = mpsc::sync_channel(DECODE_AHEAD);
    let (spare_sender, spare) = mpsc::channel();
    // The decoder waits here until a thread has it. A closure that held it would take it where it
    // goes when the thread cannot be spawned, which is nowhere: a refused request, like
    // `spawn_scoped` when it fails, drops the closure it is given.
    let waiting = std::sync::Mutex::new(Some(decoder));
    std::thread::scope(|scope| {
        let thread = std::thread::Builder::new().name("cache-decode".into());
        let waiting = &waiting;
        let spawned = ask_for_a_thread().and_then(|()| {
            thread.spawn_scoped(scope, move || {
                let decoder = waiting
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner)
                    .take()
                    .expect("the decoder is taken once, by the thread that decodes");
                decode_ahead(decoder, pieces, spare)
            })
        });
        match spawned {
            Ok(_) => deserialize_to_the_end(Decoded {
                pieces: decoded,
                spare: spare_sender,
                piece: Piece {
                    bytes: Vec::new(),
                    filled: 0,
                },
                at: 0,
            }),
            // No thread to be had: read it all in this one, as it was done before.
            Err(_) => {
                let decoder = waiting
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner)
                    .take()
                    .expect("a decoder that no thread took");
                deserialize_to_the_end(BufReader::with_capacity(DECODE_CHUNK, decoder))
            }
        }
    })
}

/// The most that `bincode` is allowed to take for one transcript, or for the length of the list: a
/// transcript of a human gene is a few hundred kilobytes with its sequences.
const LARGEST_VALUE: u64 = 64 << 20;

/// `bincode` as `save_cache` writes it, with a limit on the bytes that one call may read.
///
/// What a loader reads is not known to be what was written until the decoder has checked the
/// checksum, at the end of the frame, and `bincode` takes the length of a string or of a list from
/// what it reads and allocates it. A byte that is damaged inside the compressed data can decode to a
/// length of a terabyte, and with no limit that is an allocation no machine can give: the process
/// aborts, with no word on which cache it was and no rebuild from the GFF3. With the limit, which is
/// far above the largest transcript and far below what a corrupt length asks for, the call fails and
/// the cache is refused like any other that is damaged. It is a limit on one call, so it is one
/// transcript at a time and not the list, which is hundreds of megabytes.
fn bounded() -> impl bincode::Options {
    use bincode::Options;
    bincode::DefaultOptions::new()
        .with_fixint_encoding()
        .allow_trailing_bytes()
        .with_limit(LARGEST_VALUE)
}

/// The transcripts of `stream`, and then the rest of it, which is nothing for a cache that
/// `save_cache` wrote. Reading it is where a decoder finds a bad checksum.
///
/// The list is read as `save_cache` wrote it, its length and then the transcripts one after
/// another, each with a limit of its own.
fn deserialize_to_the_end(mut stream: impl std::io::Read) -> Result<Vec<Transcript>> {
    use bincode::Options;
    let count: u64 = bounded().deserialize_from(&mut stream)?;
    // not `with_capacity(count)`: the count is not known to be what was written either
    let mut transcripts = Vec::new();
    for _ in 0..count {
        transcripts.push(bounded().deserialize_from::<_, Transcript>(&mut stream)?);
    }
    let after = std::io::copy(&mut stream, &mut std::io::sink())
        .with_context(|| "Reading the rest of the cache after its transcripts")?;
    anyhow::ensure!(
        after == 0,
        "the cache holds {after} bytes after its {} transcripts, which `save_cache` does not write \
         (two caches joined into one file?)",
        transcripts.len()
    );
    Ok(transcripts)
}

/// A piece of the decompressed stream: `bytes[..filled]` is what the decoder put in it.
struct Piece {
    bytes: Vec<u8>,
    filled: usize,
}

/// Read `decoder` to its end, a piece at a time, into `pieces`; an error it meets is the last piece.
///
/// The pieces go round: what the deserializer is done with comes back through `spare` to be filled
/// again, so that a cache of a gigabyte or two goes through half a dozen buffers and not through a
/// thousand, each a fresh megabyte that the allocator has to find and the kernel to page in.
fn decode_ahead<R: std::io::Read>(
    mut decoder: R,
    pieces: std::sync::mpsc::SyncSender<std::io::Result<Piece>>,
    spare: std::sync::mpsc::Receiver<Vec<u8>>,
) {
    loop {
        let mut bytes = spare.try_recv().unwrap_or_else(|_| vec![0u8; DECODE_CHUNK]);
        let mut filled = 0;
        while filled < bytes.len() {
            match decoder.read(&mut bytes[filled..]) {
                // the end of the stream
                Ok(0) => break,
                Ok(n) => filled += n,
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                Err(e) => {
                    // what was read before the error is delivered before it, as a stream does
                    if filled > 0 && pieces.send(Ok(Piece { bytes, filled })).is_err() {
                        return;
                    }
                    let _ = pieces.send(Err(e));
                    return;
                }
            }
        }
        let last = filled < bytes.len();
        if filled > 0 && pieces.send(Ok(Piece { bytes, filled })).is_err() {
            // the deserializer has what it wanted, or has failed
            return;
        }
        if last {
            return;
        }
    }
}

/// The decompressed stream, as the pieces the decoder's thread sent.
struct Decoded {
    pieces: std::sync::mpsc::Receiver<std::io::Result<Piece>>,
    /// Where a piece that has been read goes, to be filled again.
    spare: std::sync::mpsc::Sender<Vec<u8>>,
    piece: Piece,
    at: usize,
}

impl std::io::Read for Decoded {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        while self.at == self.piece.filled {
            match self.pieces.recv() {
                Ok(Ok(piece)) => {
                    let done = std::mem::replace(&mut self.piece, piece);
                    self.at = 0;
                    // not to the decoder if it has gone
                    if !done.bytes.is_empty() {
                        let _ = self.spare.send(done.bytes);
                    }
                }
                Ok(Err(e)) => return Err(e),
                // the thread has sent the whole stream and gone
                Err(_) => return Ok(0),
            }
        }
        let n = out.len().min(self.piece.filled - self.at);
        out[..n].copy_from_slice(&self.piece.bytes[self.at..self.at + n]);
        self.at += n;
        Ok(n)
    }
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
        return Err(StaleCacheFormat {
            found: "FSTVEP01 (gzip)",
            summary: SUMMARY_PRE_90,
            detail: DETAIL_PRE_90,
        }
        .into());
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
        .unwrap_or_else(|_| SystemTime::now());
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
    use std::fs::File;
    use std::sync::Arc;
    use std::time::Duration;
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
            initiator_met: false,
            cds_flags_from_sequence: false,
            readthrough: false,
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
    fn a_cache_of_several_pieces_loads_whole_and_in_order() {
        // 1,200 transcripts of 4 KB that do not compress: the decoder's thread sends several pieces,
        // and a field of a transcript now and then lies across the end of one
        let transcripts = bulky_transcripts(1_200);
        let decompressed: usize = transcripts
            .iter()
            .map(|t| t.spliced_seq.as_ref().map_or(0, |seq| seq.len()))
            .sum();
        assert!(
            decompressed > 3 * DECODE_CHUNK,
            "the cache is meant to be several pieces long"
        );
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        save_cache(&transcripts, &path).unwrap();

        let loaded = load_cache(&path).unwrap();
        assert_eq!(loaded.len(), transcripts.len());
        for (read, written) in loaded.iter().zip(&transcripts) {
            assert_eq!(read.stable_id, written.stable_id);
            assert_eq!(read.spliced_seq, written.spliced_seq);
        }
    }

    /// `load_cache` in a pool of `threads`. With one or two a cache is read as a stream, with the
    /// decompression ahead of the deserialization; with more, a cache of several frames has them
    /// decoded side by side.
    fn load_on(threads: usize, path: &Path) -> Result<Vec<Transcript>> {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| load_cache(path))
    }

    /// [`load_on`], and how many workers decoded the frames if they were decoded side by side, and
    /// `None` if the cache was read as a stream.
    fn load_observed(threads: usize, path: &Path) -> (Result<Vec<Transcript>>, Option<usize>) {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| {
                // the load is made on this thread, which no other load has
                let me = std::thread::current().id();
                let loaded = load_cache(path);
                let workers = WORKERS_OF_THE_LOADS
                    .lock()
                    .unwrap()
                    .iter()
                    .find(|(id, _)| *id == me)
                    .map(|(_, workers)| *workers);
                (loaded, workers)
            })
    }

    /// Whether `error` is a stream that ended before the list did, which is the same whatever the
    /// text that the standard library gives it.
    fn is_a_stream_that_ended_early(error: &anyhow::Error) -> bool {
        error.chain().any(|cause| {
            matches!(
                cause.downcast_ref::<bincode::Error>().map(|e| &**e),
                Some(bincode::ErrorKind::Io(io)) if io.kind() == std::io::ErrorKind::UnexpectedEof
            )
        })
    }

    /// Where the frames of the file at `path` begin, and how long they are, from its table.
    fn frames_of(path: &Path) -> Vec<(u64, u64)> {
        read_frame_table(&File::open(path).unwrap())
            .unwrap()
            .expect("the cache has a table")
            .into_iter()
            .map(|frame| (frame.at, frame.len))
            .collect()
    }

    /// `save_cache` with frames of 16 KB, so that a few dozen transcripts are several frames.
    fn save_small_frames(transcripts: &[Transcript], path: &Path) {
        save_cache_with(transcripts, path, 16 << 10).unwrap();
    }

    /// What each frame of the cache at `path` holds, decompressed.
    fn parts_of(path: &Path) -> Vec<Vec<u8>> {
        let bytes = std::fs::read(path).unwrap();
        frames_of(path)
            .into_iter()
            .map(|(at, len)| {
                zstd::stream::decode_all(&bytes[at as usize..(at + len) as usize]).unwrap()
            })
            .collect()
    }

    /// A cache of `parts`, one frame each, with the table that `save_cache` would have written for
    /// them: the way to make a cache with the bytes in it that a test chooses.
    fn write_cache_of_parts(path: &Path, parts: &[Vec<u8>]) {
        use std::io::Write;
        let mut frames = Frames::new(BufWriter::new(File::create(path).unwrap()), 1).unwrap();
        for part in parts {
            frames.begin_a_frame_if_this_one_is_full().unwrap();
            frames.write_all(part).unwrap();
        }
        frames.finish().unwrap().flush().unwrap();
    }

    /// What `load_cache` says of a cache that `change` has made of the one `save_cache` wrote, as
    /// the whole chain of its message, on `threads` threads.
    fn what_a_changed_cache_says_on(
        threads: usize,
        transcripts: &[Transcript],
        change: impl FnOnce(&mut Vec<u8>),
    ) -> String {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        save_cache(transcripts, &path).unwrap();
        let mut bytes = std::fs::read(&path).unwrap();
        change(&mut bytes);
        std::fs::write(&path, &bytes).unwrap();
        match load_on(threads, &path) {
            Ok(loaded) => panic!("the changed cache loaded, as {} transcripts", loaded.len()),
            Err(e) => format!("{e:#}"),
        }
    }

    /// [`what_a_changed_cache_says_on`] as a cache is read by a run of one thread, which is as a
    /// stream.
    fn what_a_changed_cache_says(
        transcripts: &[Transcript],
        change: impl FnOnce(&mut Vec<u8>),
    ) -> String {
        what_a_changed_cache_says_on(1, transcripts, change)
    }

    #[test]
    fn a_cache_with_a_flipped_checksum_is_rejected_for_its_checksum() {
        // `save_cache` records a checksum for each frame so that a cache corrupted in place fails
        // loudly. Reading a frame to its end is what reads it, and it is the last four bytes of
        // the frame, and not of the file, which ends with the table.
        for threads in [1, 4] {
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("t.cache");
            save_small_frames(&bulky_transcripts(64), &path);
            let (at, len) = *frames_of(&path).last().unwrap();
            let mut bytes = std::fs::read(&path).unwrap();
            bytes[(at + len) as usize - 1] ^= 0xFF;
            std::fs::write(&path, &bytes).unwrap();
            let said = format!("{:#}", load_on(threads, &path).unwrap_err());
            assert!(said.contains("checksum"), "{threads} threads: {said}");
        }
    }

    #[test]
    fn a_cache_with_bytes_after_it_is_rejected() {
        // not the file `save_cache` wrote, whatever else about it reads. The table ends the file, so
        // a file that goes on after it has none to be found and is read as a stream, in which the
        // bytes are not a frame, however many threads there are.
        let junk = |bytes: &mut Vec<u8>| bytes.extend_from_slice(b"not a zstd frame");
        for threads in [1, 4] {
            let said = what_a_changed_cache_says_on(threads, &bulky_transcripts(8), junk);
            assert!(
                said.contains("Reading the rest of the cache after its transcripts"),
                "{threads} threads: {said}"
            );
        }
    }

    #[test]
    fn two_caches_joined_into_one_file_are_not_read_as_the_first_alone() {
        // What `cat a.cache b.cache` makes. The decoder reads on from the first frame into the
        // second, the deserializer stops at the end of the first list, and v0.4.0 loaded the
        // first and said nothing of the second.
        let transcripts = bulky_transcripts(8);
        let join = |bytes: &mut Vec<u8>| {
            let again = bytes.clone();
            bytes.extend_from_slice(&again);
        };
        // two caches that have a table each: read as a stream the second is what follows the first
        // list, and the table at the end of the file, which the frames path reads, is the
        // second's, and says what is in that one and not what the file holds
        let said = what_a_changed_cache_says_on(1, &transcripts, join);
        assert!(said.contains("bytes after its 8 transcripts"), "{said}");
        let said = what_a_changed_cache_says_on(4, &transcripts, join);
        assert!(
            said.contains("the frames in the table of the cache are"),
            "{said}"
        );
        // two caches as every build before the table wrote them have no table, and the second is
        // what follows the first list in a stream
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        write_a_cache_of_one_frame(&path, &transcripts);
        let mut bytes = std::fs::read(&path).unwrap();
        join(&mut bytes);
        std::fs::write(&path, &bytes).unwrap();
        for threads in [1, 4] {
            let said = format!("{:#}", load_on(threads, &path).unwrap_err());
            assert!(
                said.contains("bytes after its 8 transcripts"),
                "{threads} threads: the message names what is wrong with it: {said}"
            );
        }
    }

    #[test]
    fn a_byte_flipped_in_the_middle_of_a_cache_is_rejected() {
        // v0.4.0 loaded this as the transcripts it had become
        what_a_changed_cache_says(&bulky_transcripts(64), |bytes| {
            let middle = bytes.len() / 2;
            bytes[middle] ^= 0xFF;
        });
    }

    /// Reads from `inner`, and remembers the threads that read it and the buffers it was asked to
    /// fill whole.
    struct Spy {
        inner: std::io::Cursor<Vec<u8>>,
        threads: Arc<std::sync::Mutex<Vec<std::thread::ThreadId>>>,
        buffers: Arc<std::sync::Mutex<Vec<usize>>>,
    }

    impl Spy {
        #[allow(clippy::type_complexity)] // the three things a test reads back, named where they are made
        fn over(
            bytes: Vec<u8>,
        ) -> (
            Self,
            Arc<std::sync::Mutex<Vec<std::thread::ThreadId>>>,
            Arc<std::sync::Mutex<Vec<usize>>>,
        ) {
            let threads = Arc::new(std::sync::Mutex::new(Vec::new()));
            let buffers = Arc::new(std::sync::Mutex::new(Vec::new()));
            let spy = Self {
                inner: std::io::Cursor::new(bytes),
                threads: threads.clone(),
                buffers: buffers.clone(),
            };
            (spy, threads, buffers)
        }
    }

    impl std::io::Read for Spy {
        fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
            self.threads
                .lock()
                .unwrap()
                .push(std::thread::current().id());
            if out.len() == DECODE_CHUNK {
                self.buffers.lock().unwrap().push(out.as_ptr() as usize);
            }
            self.inner.read(out)
        }
    }

    /// Runs `f`, which has to answer in `seconds`: a test of something that would hang if it were
    /// wrong fails, and does not stop the run.
    fn within<T: Send + 'static>(seconds: u64, f: impl FnOnce() -> T + Send + 'static) -> T {
        let (answer, answered) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            let _ = answer.send(f());
        });
        match answered.recv_timeout(std::time::Duration::from_secs(seconds)) {
            Ok(value) => value,
            Err(std::sync::mpsc::RecvTimeoutError::Timeout) => {
                panic!("no answer in {seconds} seconds: it is waiting for something that will not come")
            }
            Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => {
                panic!("the thread that was to answer panicked")
            }
        }
    }

    #[test]
    fn the_decoder_is_read_on_a_thread_of_its_own() {
        // What makes the load faster is that the decompression and the deserialization overlap,
        // and nothing else of a test would notice if they did not: the answers are the same. So
        // this holds the property itself, that the decoder is read by a thread other than the one
        // that asked for the transcripts.
        let transcripts = bulky_transcripts(600);
        let stream = bincode::serialize(&transcripts).unwrap();
        assert!(stream.len() > 2 * DECODE_CHUNK, "several pieces long");
        let (spy, threads, _) = Spy::over(stream);
        let loaded = deserialize_decoded(spy).unwrap();
        assert_eq!(loaded.len(), transcripts.len());
        let threads = threads.lock().unwrap();
        assert!(threads.len() >= 3, "{} reads", threads.len());
        let this = std::thread::current().id();
        assert!(
            threads.iter().all(|t| *t != this),
            "the decoder was read by the thread that deserializes"
        );
    }

    #[test]
    fn the_pieces_go_round_through_a_few_buffers() {
        // 64 pieces and the buffers they went through: the deserializer returns each when it is
        // done with it, and the decoder fills it again
        let transcripts = bulky_transcripts(16_000);
        let stream = bincode::serialize(&transcripts).unwrap();
        let pieces = stream.len().div_ceil(DECODE_CHUNK);
        assert!(pieces >= 50, "{pieces} pieces");
        let (spy, _, buffers) = Spy::over(stream);
        let loaded = deserialize_decoded(spy).unwrap();
        assert_eq!(loaded.len(), transcripts.len());
        let mut distinct = buffers.lock().unwrap().clone();
        distinct.sort_unstable();
        distinct.dedup();
        assert!(
            distinct.len() <= DECODE_AHEAD + 3,
            "{} buffers for {pieces} pieces",
            distinct.len()
        );
    }

    #[test]
    fn a_load_that_cannot_have_a_thread_decodes_on_the_calling_one() {
        // The host refuses the decoder's thread. The load decodes in the thread that asked for it,
        // as it did before there was a thread, and the answer is the same. The decoder is not in
        // what the refused request drops, or there would be nothing to decode.
        let transcripts = bulky_transcripts(40);
        let stream = bincode::serialize(&transcripts).unwrap();
        let (spy, threads, _) = Spy::over(stream);
        let loaded = deserialize_decoded_with(spy, || {
            Err(std::io::Error::from(std::io::ErrorKind::WouldBlock))
        })
        .unwrap();
        assert_eq!(loaded.len(), transcripts.len());
        for (read, written) in loaded.iter().zip(&transcripts) {
            assert_eq!(read.stable_id, written.stable_id);
            assert_eq!(read.spliced_seq, written.spliced_seq);
        }
        let threads = threads.lock().unwrap();
        assert!(!threads.is_empty());
        let this = std::thread::current().id();
        assert!(
            threads.iter().all(|t| *t == this),
            "no thread was to be had, and the decoder was read by one"
        );
    }

    #[test]
    fn a_stream_that_is_a_whole_number_of_pieces_long_loads() {
        // One more piece would be empty: the decoder's thread has to see the end of the stream
        // in the read that returns nothing, and not wait for a piece that is not coming
        let mut transcripts = bulky_transcripts(300);
        let size = bincode::serialized_size(&transcripts).unwrap() as usize;
        let pad = DECODE_CHUNK - size % DECODE_CHUNK;
        transcripts
            .last_mut()
            .unwrap()
            .spliced_seq
            .as_mut()
            .unwrap()
            .push_str(&"A".repeat(pad));
        let stream = bincode::serialize(&transcripts).unwrap();
        assert_eq!(stream.len() % DECODE_CHUNK, 0);
        let loaded = within(60, move || {
            deserialize_decoded(std::io::Cursor::new(stream))
        })
        .unwrap();
        assert_eq!(loaded.len(), 300);
    }

    /// Fails with `Interrupted` before every other read, as a read does when a signal arrives.
    struct Interrupted {
        inner: std::io::Cursor<Vec<u8>>,
        calls: usize,
    }

    impl std::io::Read for Interrupted {
        fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
            self.calls += 1;
            if self.calls % 2 == 1 {
                return Err(std::io::Error::from(std::io::ErrorKind::Interrupted));
            }
            // and short reads, which a decoder gives
            let n = out.len().min(100_000);
            self.inner.read(&mut out[..n])
        }
    }

    #[test]
    fn an_interrupted_read_of_the_decoder_is_read_again() {
        let transcripts = bulky_transcripts(300);
        let stream = bincode::serialize(&transcripts).unwrap();
        let reader = Interrupted {
            inner: std::io::Cursor::new(stream),
            calls: 0,
        };
        let loaded = deserialize_decoded(reader).unwrap();
        assert_eq!(loaded.len(), transcripts.len());
        assert_eq!(loaded[299].spliced_seq, transcripts[299].spliced_seq);
    }

    /// Reads `good` and then fails, as a decoder does on a corrupt stream.
    struct FailsAfter {
        good: usize,
    }

    impl std::io::Read for FailsAfter {
        fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
            if self.good == 0 {
                return Err(std::io::Error::other("corrupt on purpose"));
            }
            let n = out.len().min(self.good);
            out[..n].fill(0);
            self.good -= n;
            Ok(n)
        }
    }

    #[test]
    fn an_error_of_the_decoder_reaches_the_deserializer_after_the_bytes_before_it() {
        use std::io::Read;
        let (pieces, decoded) = std::sync::mpsc::sync_channel(DECODE_AHEAD);
        let (spare_sender, spare) = std::sync::mpsc::channel();
        let thread = std::thread::spawn(move || {
            decode_ahead(
                FailsAfter {
                    good: 3 * DECODE_CHUNK + 5,
                },
                pieces,
                spare,
            )
        });
        let mut reader = Decoded {
            pieces: decoded,
            spare: spare_sender,
            piece: Piece {
                bytes: Vec::new(),
                filled: 0,
            },
            at: 0,
        };
        let mut got = Vec::new();
        let error = reader.read_to_end(&mut got).unwrap_err();
        thread.join().unwrap();
        assert_eq!(
            got.len(),
            3 * DECODE_CHUNK + 5,
            "every byte before the error is delivered"
        );
        assert_eq!(error.to_string(), "corrupt on purpose");
    }

    #[test]
    fn a_stream_longer_than_its_list_is_read_over_and_is_an_error() {
        // An empty list, and then 32 MiB of zeros that decode cleanly: the decoder's thread is
        // blocked on a full channel when the list is done, and is read to the end all the same,
        // for what it finds there
        use std::io::Read;
        let stream = std::io::Cursor::new(0u64.to_le_bytes().to_vec())
            .chain(std::io::repeat(0).take(32 << 20));
        let said = within(60, move || deserialize_decoded(stream)).unwrap_err();
        assert!(
            format!("{said:#}").contains("33554432 bytes after its 0 transcripts"),
            "{said:#}"
        );
    }

    #[test]
    fn a_length_that_is_not_one_is_an_error_and_not_an_allocation() {
        // The decoder checks the checksum of a frame at its end, after the deserializer has taken
        // what the frame says, and a damaged byte inside the compressed data can decode to a length
        // of a gigabyte or an exabyte. The length of the id of the first transcript is made 1 GiB
        // here: with no limit that is allocated (and filled with zeros) before the missing bytes
        // are found, and an exabyte would not be an error at all but an abort.
        let mut stream = bincode::serialize(&bulky_transcripts(4)).unwrap();
        stream[8..16].copy_from_slice(&(1u64 << 30).to_le_bytes());
        let said = within(60, move || {
            deserialize_decoded(std::io::Cursor::new(stream))
        })
        .unwrap_err();
        assert!(format!("{said:#}").contains("size limit"), "{said:#}");
    }

    #[test]
    fn a_deserializer_that_fails_does_not_wait_for_a_decoder_that_has_no_end() {
        use std::io::Read;
        // A list of one transcript whose first field is a one-byte string that is not UTF-8,
        // and after it a stream with no end: the deserializer fails on the first byte of the
        // string, and the thread, blocked on a full channel, has to be let go for this to return
        // at all.
        let mut prefix = Vec::new();
        prefix.extend_from_slice(&1u64.to_le_bytes());
        prefix.extend_from_slice(&1u64.to_le_bytes());
        prefix.push(0xFF);
        let stream = std::io::Cursor::new(prefix).chain(std::io::repeat(0));
        assert!(within(60, move || deserialize_decoded(stream)).is_err());
    }

    /// Panics at its first read.
    struct Panics;

    impl std::io::Read for Panics {
        fn read(&mut self, _: &mut [u8]) -> std::io::Result<usize> {
            panic!("the decoder panicked on purpose")
        }
    }

    #[test]
    fn a_panic_of_the_decoder_is_the_loads_panic_and_not_a_short_cache() {
        // To the deserializer a thread that has gone is the end of the stream, which for a stream
        // that stops short of its list is an error, and for one that stops after it would be a
        // cache with the rest of it missing. The scope the thread runs in re-raises its panic, so
        // that neither is what the caller gets.
        let outcome = within(60, || {
            std::panic::catch_unwind(|| deserialize_decoded(Panics))
        });
        assert!(outcome.is_err(), "the load returned {outcome:?}");
    }

    /// A cache as every build before the one that cut it into frames wrote it: the whole of it one
    /// frame, and no table.
    fn write_a_cache_of_one_frame(path: &Path, transcripts: &[Transcript]) {
        use std::io::Write;
        let mut zst = zstd::Encoder::new(BufWriter::new(File::create(path).unwrap()), 1).unwrap();
        zst.include_checksum(true).unwrap();
        zst.write_all(CACHE_MAGIC_V10).unwrap();
        bincode::serialize_into(&mut zst, transcripts).unwrap();
        zst.finish().unwrap().flush().unwrap();
    }

    fn same_transcripts(read: &[Transcript], written: &[Transcript]) {
        assert_eq!(read.len(), written.len());
        for (read, written) in read.iter().zip(written) {
            assert_eq!(read.stable_id, written.stable_id);
            assert_eq!(read.spliced_seq, written.spliced_seq);
        }
    }

    #[test]
    fn a_big_cache_is_several_frames_and_the_bytes_of_one_stream() {
        // 5,000 transcripts of 4 KB that do not compress: 20 MB, which is more than two frames
        let transcripts = bulky_transcripts(5_000);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        save_cache(&transcripts, &path).unwrap();
        let bytes = std::fs::read(&path).unwrap();

        let frames = frames_of(&path);
        assert!(frames.len() >= 3, "{} frames", frames.len());
        // the frames, and then the table: 8 bytes for the header of its frame, 8 for each frame,
        // and 9 for the footer
        let table = 8 + 8 * frames.len() as u64 + 9;
        assert_eq!(
            frames.iter().map(|f| f.1).sum::<u64>() + table,
            bytes.len() as u64
        );
        // The file, read as one stream, is the magic and what `bincode` writes of the list, which
        // is what the one frame held: a reader of one stream reads the same bytes, and passes
        // over the table.
        let mut stream = CACHE_MAGIC_V10.to_vec();
        stream.extend_from_slice(&bincode::serialize(&transcripts).unwrap());
        assert!(
            zstd::stream::decode_all(&bytes[..]).unwrap() == stream,
            "the frames of the cache, read as one stream, are not the stream of the one frame"
        );
    }

    #[test]
    fn the_table_names_each_frame_of_the_file_and_what_it_decompresses_to() {
        // none, one and many transcripts, in frames of a byte (a frame for each), of a hundred, of a
        // few kilobytes, and one that is never full
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        for (n, frame_bytes) in [
            (0, 1),
            (1, 1),
            (7, 100),
            (7, 1 << 30),
            (60, 16 << 10),
            (60, 1),
        ] {
            let transcripts = bulky_transcripts(n);
            save_cache_with(&transcripts, &path, frame_bytes).unwrap();
            let bytes = std::fs::read(&path).unwrap();
            let frames = frames_of(&path);
            let (last_at, last_len) = *frames.last().unwrap();
            let table = &bytes[(last_at + last_len) as usize..];
            let word = |at: usize| u32::from_le_bytes(table[at..at + 4].try_into().unwrap());
            assert_eq!(word(0), 0x184D_2A5E, "the magic of a skippable frame");
            assert_eq!(word(4) as usize, table.len() - 8, "the size of that frame");
            assert_eq!(
                word(table.len() - 9) as usize,
                frames.len(),
                "{n} {frame_bytes}"
            );
            assert_eq!(table[table.len() - 5], 0, "no flags");
            assert_eq!(
                word(table.len() - 4),
                0x8F92_EAB1,
                "the magic of the footer"
            );

            let mut stream = Vec::new();
            for (k, &(at, len)) in frames.iter().enumerate() {
                let frame = &bytes[at as usize..(at + len) as usize];
                // an entry is one frame, whole: the library, which is given the bytes, says where a
                // frame that begins here ends
                assert_eq!(
                    zstd::zstd_safe::find_frame_compressed_size(frame).unwrap(),
                    frame.len(),
                    "frame {k} of {n} transcripts in frames of {frame_bytes}"
                );
                let decompressed = zstd::stream::decode_all(frame).unwrap();
                assert_eq!(word(8 + 8 * k), len as u32);
                assert_eq!(word(12 + 8 * k) as usize, decompressed.len());
                stream.extend(decompressed);
            }
            let mut expected = CACHE_MAGIC_V10.to_vec();
            expected.extend_from_slice(&bincode::serialize(&transcripts).unwrap());
            assert!(
                stream == expected,
                "{n} transcripts in frames of {frame_bytes}"
            );
            if frame_bytes == 1 {
                assert_eq!(
                    frames.len(),
                    n + 1,
                    "the header is a frame, and then one for each transcript"
                );
            }
            for threads in [1, 4] {
                same_transcripts(&load_on(threads, &path).unwrap(), &transcripts);
            }
        }
    }

    #[test]
    fn a_cache_of_several_frames_loads_the_same_on_any_number_of_threads() {
        let transcripts = bulky_transcripts(200);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        save_small_frames(&transcripts, &path);
        assert!(frames_of(&path).len() >= 20);
        // 32 is more threads than there are decoders
        for threads in [1, 2, 3, 8, 32] {
            same_transcripts(&load_on(threads, &path).unwrap(), &transcripts);
        }
    }

    #[test]
    fn at_most_eight_frames_are_decoded_at_once_and_one_thread_of_the_pool_merges() {
        assert_eq!(decoders_on(1), 1);
        assert_eq!(decoders_on(3), 2);
        assert_eq!(decoders_on(8), 7);
        assert_eq!(decoders_on(9), MOST_DECODERS);
        assert_eq!(decoders_on(64), MOST_DECODERS);
        // and two threads or fewer have no cores to spare for them
        assert_eq!(side_by_side_workers(1), None);
        assert_eq!(side_by_side_workers(2), None);
        assert_eq!(side_by_side_workers(3), Some(2));
        assert_eq!(side_by_side_workers(64), Some(MOST_DECODERS));
    }

    #[test]
    fn a_load_takes_the_frames_from_the_table_on_three_threads_or_more_and_with_as_many_workers_as_that_says(
    ) {
        let transcripts = bulky_transcripts(200);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        save_small_frames(&transcripts, &path);
        // the table is what is read to decode side by side, so with two threads or fewer it is not
        // read, and is not what is decoded with
        for (threads, workers) in [
            (1, None),
            (2, None),
            (3, Some(2)),
            (8, Some(7)),
            (9, Some(8)),
            (32, Some(8)),
        ] {
            let (loaded, observed) = load_observed(threads, &path);
            same_transcripts(&loaded.unwrap(), &transcripts);
            assert_eq!(observed, workers, "{threads} threads");
        }
        // a cache of one frame, and one with no table, are read as a stream at any number of them
        write_a_cache_of_one_frame(&path, &transcripts);
        assert_eq!(load_observed(8, &path).1, None);
        save_small_frames(&transcripts, &path);
        let bytes = std::fs::read(&path).unwrap();
        let (at, len) = *frames_of(&path).last().unwrap();
        std::fs::write(&path, &bytes[..(at + len) as usize]).unwrap();
        assert_eq!(load_observed(8, &path).1, None);
    }

    #[test]
    fn items_are_merged_in_order_whatever_order_they_are_decoded_in() {
        let merged = within(60, || {
            let items: Vec<u64> = (0..40).collect();
            let mut merged = Vec::new();
            decode_in_order(
                &items,
                4,
                |i, item| {
                    // the early items are the slow ones, so that the later ones are done first
                    let slow = if i < 8 { 40 - 4 * i as u64 } else { 0 };
                    std::thread::sleep(Duration::from_millis(slow));
                    Ok(*item * 10)
                },
                |value| {
                    merged.push(value);
                    Ok(())
                },
            )
            .unwrap();
            merged
        });
        assert_eq!(merged, (0..40).map(|i| i * 10).collect::<Vec<u64>>());
    }

    #[test]
    fn the_workers_decode_side_by_side_and_no_more_items_are_decoded_at_once_than_there_are_workers(
    ) {
        use std::sync::atomic::{AtomicUsize, Ordering};
        use std::sync::{Arc, Barrier};
        let most = within(60, || {
            let workers = 3;
            let barrier = Arc::new(Barrier::new(workers));
            let running = Arc::new(AtomicUsize::new(0));
            let most = Arc::new(AtomicUsize::new(0));
            let items: Vec<usize> = (0..60).collect();
            let (b, r, m) = (barrier.clone(), running.clone(), most.clone());
            decode_in_order(
                &items,
                workers,
                move |i, _| {
                    let now = r.fetch_add(1, Ordering::SeqCst) + 1;
                    m.fetch_max(now, Ordering::SeqCst);
                    // the first item of each worker waits for the others, so they have to be
                    // decoding at the same time, or this does not return
                    if i < workers {
                        b.wait();
                    }
                    std::thread::sleep(Duration::from_millis(2));
                    r.fetch_sub(1, Ordering::SeqCst);
                    Ok(())
                },
                |_| Ok(()),
            )
            .unwrap();
            most.load(Ordering::SeqCst)
        });
        assert_eq!(most, 3);
    }

    #[test]
    fn the_first_error_of_a_decode_is_returned_and_the_items_not_yet_begun_are_left_alone() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        use std::sync::Arc;
        let (error, begun) = within(60, || {
            let begun = Arc::new(AtomicUsize::new(0));
            let counted = begun.clone();
            let items: Vec<u64> = (0..400).collect();
            let result = decode_in_order(
                &items,
                4,
                move |i, _| {
                    counted.fetch_add(1, Ordering::SeqCst);
                    if i == 5 {
                        return Err(anyhow::anyhow!("item 5 is bad"));
                    }
                    std::thread::sleep(Duration::from_millis(20));
                    Ok(i)
                },
                |_| Ok(()),
            );
            (
                format!("{:#}", result.unwrap_err()),
                begun.load(Ordering::SeqCst),
            )
        });
        assert!(error.contains("item 5 is bad"), "{error}");
        // four workers begin four items, and four more when those are done, which is when item 5
        // has failed; every worker then looks at the flag before it begins another
        assert!(begun < 200, "{begun} of 400 items were begun");
    }

    #[test]
    fn an_error_of_the_merge_is_returned_and_the_items_not_yet_begun_are_left_alone() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        use std::sync::Arc;
        let (error, begun) = within(60, || {
            let begun = Arc::new(AtomicUsize::new(0));
            let counted = begun.clone();
            let items: Vec<u64> = (0..400).collect();
            let result = decode_in_order(
                &items,
                4,
                move |i, _| {
                    counted.fetch_add(1, Ordering::SeqCst);
                    std::thread::sleep(Duration::from_millis(20));
                    Ok(i)
                },
                |value| {
                    anyhow::ensure!(value != 3, "item 3 cannot be merged");
                    Ok(())
                },
            );
            (
                format!("{:#}", result.unwrap_err()),
                begun.load(Ordering::SeqCst),
            )
        });
        assert!(error.contains("item 3 cannot be merged"), "{error}");
        assert!(begun < 200, "{begun} of 400 items were begun");
    }

    #[test]
    fn there_may_be_more_workers_than_items_and_no_items_at_all() {
        for (n, workers) in [(0, 3), (1, 3), (2, 8), (5, 1), (7, 0)] {
            let items: Vec<usize> = (0..n).collect();
            let mut merged = Vec::new();
            decode_in_order(
                &items,
                workers,
                |_, item| Ok(*item),
                |value| {
                    merged.push(value);
                    Ok(())
                },
            )
            .unwrap();
            assert_eq!(merged, items, "{n} items on {workers} workers");
        }
    }

    #[test]
    fn workers_that_cannot_have_a_thread_leave_the_items_to_those_that_can_and_to_the_caller() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let items: Vec<usize> = (0..30).collect();
        // The host refuses every thread, or every one after the first two: the same items are
        // merged in the same order, on the threads that there are, and on this one with none.
        for allowed in [0, 2] {
            let asked = AtomicUsize::new(0);
            let mut merged = Vec::new();
            decode_in_order_with(
                &items,
                5,
                || {
                    if asked.fetch_add(1, Ordering::SeqCst) < allowed {
                        Ok(())
                    } else {
                        Err(std::io::ErrorKind::WouldBlock.into())
                    }
                },
                |_, item| Ok(*item),
                |value| {
                    merged.push(value);
                    Ok(())
                },
            )
            .unwrap();
            assert_eq!(merged, items, "{allowed} threads were given");
        }
        // and an error is the error of the load there too
        let said = decode_in_order_with(
            &items,
            5,
            || Err(std::io::ErrorKind::WouldBlock.into()),
            |i, _| {
                anyhow::ensure!(i != 4, "item 4 is bad");
                Ok(())
            },
            |_| Ok(()),
        )
        .unwrap_err();
        assert!(format!("{said:#}").contains("item 4 is bad"), "{said:#}");
    }

    #[test]
    fn a_panic_of_a_decode_or_of_the_merge_is_the_panic_of_the_load_and_stops_the_items_not_yet_begun(
    ) {
        use std::sync::atomic::{AtomicUsize, Ordering};
        use std::sync::Arc;
        for panics_in in ["the decode", "the merge"] {
            let begun = Arc::new(AtomicUsize::new(0));
            let counted = begun.clone();
            let outcome = within(60, move || {
                let items: Vec<usize> = (0..400).collect();
                std::panic::catch_unwind(|| {
                    let _ = decode_in_order(
                        &items,
                        4,
                        |i, _| {
                            counted.fetch_add(1, Ordering::SeqCst);
                            assert!(!(panics_in == "the decode" && i == 7), "item 7 panics");
                            std::thread::sleep(Duration::from_millis(20));
                            Ok(i)
                        },
                        |value| {
                            assert!(!(panics_in == "the merge" && value == 3), "item 3 panics");
                            Ok(())
                        },
                    );
                })
            });
            assert!(outcome.is_err(), "{panics_in}: the load returned");
            // a panic that is not told to the workers is found once every item has been decoded
            let begun = begun.load(Ordering::SeqCst);
            assert!(begun < 200, "{panics_in}: {begun} of 400 items were begun");
        }
    }

    #[test]
    fn loads_at_once_in_a_pool_that_has_no_thread_to_spare_do_not_wait_for_one_another() {
        // The merge of a load holds the thread it is called from while it waits for the decoders.
        // As tasks of the pool they would wait for threads that are all held by loads, and this
        // would never return (it did, when they were).
        within(60, || {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(2)
                .build()
                .unwrap();
            let items: Vec<usize> = (0..20).collect();
            pool.scope(|scope| {
                for _ in 0..4 {
                    scope.spawn(|_| {
                        let mut merged = Vec::new();
                        decode_in_order(
                            &items,
                            3,
                            |_, item| {
                                std::thread::sleep(Duration::from_millis(2));
                                Ok(*item)
                            },
                            |value| {
                                merged.push(value);
                                Ok(())
                            },
                        )
                        .unwrap();
                        assert_eq!(merged, items);
                    });
                }
            });
        });
    }

    #[test]
    fn a_loader_that_reads_one_stream_reads_a_cache_of_several_frames() {
        // what every build before the one that cut the frames did with the file
        use std::io::Read;
        let transcripts = bulky_transcripts(200);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        save_small_frames(&transcripts, &path);
        let mut zst = zstd::Decoder::new(File::open(&path).unwrap()).unwrap();
        let mut magic = [0u8; 8];
        zst.read_exact(&mut magic).unwrap();
        assert_eq!(&magic, CACHE_MAGIC_V10);
        let loaded: Vec<Transcript> = bincode::deserialize_from(zst).unwrap();
        same_transcripts(&loaded, &transcripts);
    }

    #[test]
    fn a_cache_of_one_frame_that_an_earlier_build_wrote_still_loads() {
        let transcripts = bulky_transcripts(1_200);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        write_a_cache_of_one_frame(&path, &transcripts);
        assert!(read_frame_table(&File::open(&path).unwrap())
            .unwrap()
            .is_none());
        for threads in [1, 4] {
            same_transcripts(&load_on(threads, &path).unwrap(), &transcripts);
        }
    }

    #[test]
    fn frames_that_end_where_a_block_or_a_buffer_does_load_on_every_path() {
        // zstd blocks are 128 KiB, the buffer a frame is read through is 256 KiB, and a piece of
        // the stream of the other path is a MiB: a frame of that many bytes, or one more or fewer,
        // ends where one of them does. A cache of the transcripts of Ensembl's GFF3 has about one
        // such frame in six hundred.
        let transcript_of = |i: usize, bytes: usize| {
            let mut transcript = make_test_transcript();
            transcript.stable_id = Arc::from(format!("ENST{i:08}").as_str());
            transcript.spliced_seq = Some(String::new());
            let empty = bincode::serialized_size(&transcript).unwrap() as usize;
            // pseudo-random bases, which zstd makes a quarter of and not nothing of
            let mut x = (i as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
            transcript.spliced_seq = Some(
                (0..bytes - empty)
                    .map(|_| {
                        x ^= x << 13;
                        x ^= x >> 7;
                        x ^= x << 17;
                        ['A', 'C', 'G', 'T'][(x % 4) as usize]
                    })
                    .collect(),
            );
            assert_eq!(
                bincode::serialized_size(&transcript).unwrap() as usize,
                bytes
            );
            transcript
        };
        let sizes: Vec<usize> = [128 << 10, 256 << 10, 1 << 20]
            .into_iter()
            .flat_map(|size| [size - 1, size, size + 1])
            .collect();
        let transcripts: Vec<Transcript> = sizes
            .iter()
            .enumerate()
            // the first frame has the magic and the length of the list before its transcript
            .map(|(i, &size)| transcript_of(i, if i == 0 { size - 16 } else { size }))
            .collect();
        let mut parts: Vec<Vec<u8>> = transcripts
            .iter()
            .map(|t| bincode::serialize(t).unwrap())
            .collect();
        let mut first = CACHE_MAGIC_V10.to_vec();
        first.extend_from_slice(&(transcripts.len() as u64).to_le_bytes());
        first.append(&mut parts[0]);
        parts[0] = first;
        assert_eq!(parts.iter().map(Vec::len).collect::<Vec<_>>(), sizes);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        write_cache_of_parts(&path, &parts);
        for threads in [1, 3, 4, 8] {
            same_transcripts(&load_on(threads, &path).unwrap(), &transcripts);
        }
    }

    #[test]
    fn a_transcript_that_a_cache_could_not_be_read_with_is_not_written() {
        // A cache is read with a limit on the bytes of one transcript. One over it would be written,
        // refused at every load, and written again by every run that rebuilds it.
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        let mut big = make_test_transcript();
        big.spliced_seq = Some("A".repeat(LARGEST_VALUE as usize));
        let said = format!("{:#}", save_cache(&[big.clone()], &path).unwrap_err());
        assert!(said.contains("a cache is read with a limit of"), "{said}");
        assert!(
            std::fs::read_dir(dir.path()).unwrap().next().is_none(),
            "nothing is left, neither the cache nor a part of one"
        );
        // and one under it is written and read
        big.spliced_seq = Some("A".repeat(LARGEST_VALUE as usize - 4096));
        save_cache(&[big.clone()], &path).unwrap();
        same_transcripts(&load_cache(&path).unwrap(), &[big]);
    }

    #[test]
    fn a_cache_without_a_table_is_read_as_a_stream_of_the_frames_it_has() {
        // The table is where the frames are, and not what is in them: a file that has lost it (cut
        // off after its last frame, or with a footer that no longer says it is one) is still the
        // frames, each with its checksum, and loads as every cache did before there was a table.
        let transcripts = bulky_transcripts(100);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        save_small_frames(&transcripts, &path);
        let bytes = std::fs::read(&path).unwrap();
        let (at, len) = *frames_of(&path).last().unwrap();
        let without_table = &bytes[..(at + len) as usize];
        let mut footer_flipped = bytes.clone();
        *footer_flipped.last_mut().unwrap() ^= 0xFF;
        for changed in [without_table, &footer_flipped[..]] {
            std::fs::write(&path, changed).unwrap();
            assert!(read_frame_table(&File::open(&path).unwrap())
                .unwrap()
                .is_none());
            for threads in [1, 4] {
                same_transcripts(&load_on(threads, &path).unwrap(), &transcripts);
            }
        }
    }

    #[test]
    fn a_table_that_does_not_describe_the_file_is_refused_where_the_frames_are_taken_from_it() {
        let transcripts = bulky_transcripts(60);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        save_small_frames(&transcripts, &path);
        let bytes = std::fs::read(&path).unwrap();
        let frames = frames_of(&path);
        assert!(frames.len() >= 5, "{} frames", frames.len());
        let (last_at, last_len) = *frames.last().unwrap();
        let table_at = (last_at + last_len) as usize;
        let n = frames.len();
        let put = |bytes: &mut Vec<u8>, at: usize, value: u32| {
            bytes[at..at + 4].copy_from_slice(&value.to_le_bytes())
        };
        let get =
            |bytes: &Vec<u8>, at: usize| u32::from_le_bytes(bytes[at..at + 4].try_into().unwrap());
        // what each change to the table is refused for, on a run that has the cores to use it
        type Change = Box<dyn Fn(&mut Vec<u8>)>;
        let cases: Vec<(&str, Change, &str)> = vec![
            (
                "the first frame one byte longer than it is",
                Box::new(move |b| put(b, table_at + 8, get(b, table_at + 8) + 1)),
                "the frames in the table of the cache are",
            ),
            (
                "the first frame one byte shorter than it is",
                Box::new(move |b| put(b, table_at + 8, get(b, table_at + 8) - 1)),
                "the frames in the table of the cache are",
            ),
            (
                "a byte moved from the second frame to the first, so that they still add up",
                Box::new(move |b| {
                    put(b, table_at + 8, get(b, table_at + 8) + 1);
                    put(b, table_at + 16, get(b, table_at + 16) - 1);
                }),
                "Deserializing transcripts from cache",
            ),
            (
                "two frames under one entry, the frames in it more than the one it is for",
                Box::new(move |b| {
                    let entry = |b: &Vec<u8>, k: usize| {
                        (get(b, table_at + 8 + 8 * k), get(b, table_at + 12 + 8 * k))
                    };
                    let mut table = 0x184D_2A5Eu32.to_le_bytes().to_vec();
                    table.extend_from_slice(&(8 * (n as u32 - 1) + 9).to_le_bytes());
                    let (c0, d0) = entry(b, 0);
                    let (c1, d1) = entry(b, 1);
                    table.extend_from_slice(&(c0 + c1).to_le_bytes());
                    table.extend_from_slice(&(d0 + d1).to_le_bytes());
                    for k in 2..n {
                        let (c, d) = entry(b, k);
                        table.extend_from_slice(&c.to_le_bytes());
                        table.extend_from_slice(&d.to_le_bytes());
                    }
                    table.extend_from_slice(&(n as u32 - 1).to_le_bytes());
                    table.push(0);
                    table.extend_from_slice(&0x8F92_EAB1u32.to_le_bytes());
                    b.truncate(table_at);
                    b.extend_from_slice(&table);
                }),
                "inside its entry in the table",
            ),
            (
                "the first frame said to decompress to one byte more than it does",
                Box::new(move |b| put(b, table_at + 12, get(b, table_at + 12) + 1)),
                "decompresses to",
            ),
            (
                "the last frame said to decompress to one byte fewer than it does",
                Box::new(move |b| {
                    let at = table_at + 12 + 8 * (n - 1);
                    put(b, at, get(b, at) - 1)
                }),
                "decompresses to",
            ),
            (
                "a frame of no bytes",
                Box::new(move |b| {
                    put(
                        b,
                        table_at + 16,
                        get(b, table_at + 16) + get(b, table_at + 8),
                    );
                    put(b, table_at + 8, 0);
                }),
                "a frame of no bytes",
            ),
            (
                "one frame more than the table has entries for",
                Box::new(move |b| {
                    let footer = table_at + 8 + 8 * n;
                    put(b, footer, n as u32 + 1)
                }),
                "not in a skippable frame of its size",
            ),
            (
                "more frames than the file can hold",
                Box::new(move |b| {
                    let footer = table_at + 8 + 8 * n;
                    put(b, footer, u32::MAX)
                }),
                "in a file of",
            ),
            (
                "no frames",
                Box::new(move |b| {
                    let footer = table_at + 8 + 8 * n;
                    put(b, footer, 0)
                }),
                "in a file of",
            ),
            (
                "flags that this does not read",
                Box::new(move |b| {
                    let footer = table_at + 8 + 8 * n;
                    b[footer + 4] = 0x80;
                }),
                "flags that this does not read",
            ),
            (
                "the magic of the skippable frame",
                Box::new(move |b| b[table_at] ^= 0x01),
                "not in a skippable frame of its size",
            ),
            (
                "the size of the skippable frame",
                Box::new(move |b| put(b, table_at + 4, 8 * n as u32 + 10)),
                "not in a skippable frame of its size",
            ),
        ];
        for (what, change, said) in &cases {
            let mut damaged = bytes.clone();
            change(&mut damaged);
            std::fs::write(&path, &damaged).unwrap();
            let got = format!("{:#}", load_on(4, &path).unwrap_err());
            assert!(got.contains(said), "{what}: {got}");
            // A run that does not take the frames from the table does not read it, and the frames
            // are intact, each with its own checksum. The skippable frame is part of the stream,
            // though, and a decoder that is told by its header that it is longer than the file is
            // refuses that.
            let as_a_stream = load_on(1, &path);
            if *what == "the size of the skippable frame" {
                assert!(as_a_stream.is_err());
            } else {
                same_transcripts(&as_a_stream.unwrap(), &transcripts);
            }
        }
    }

    #[test]
    fn a_stale_cache_of_several_frames_is_refused_as_stale_and_not_loaded() {
        // The next bump of the magic depends on this: the magic is in the first frame, and the
        // frames are decoded side by side, where a cache that is only stale is whole in every
        // other way and loads.
        let transcripts = bulky_transcripts(60);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        save_small_frames(&transcripts, &path);
        let mut parts = parts_of(&path);
        assert!(parts.len() >= 5);
        parts[0][..8].copy_from_slice(CACHE_MAGIC_V9);
        write_cache_of_parts(&path, &parts);
        for threads in [1, 4] {
            let error = load_on(threads, &path).unwrap_err();
            let stale = error
                .downcast_ref::<StaleCacheFormat>()
                .unwrap_or_else(|| panic!("{threads} threads: not a stale format: {error:#}"));
            assert_eq!(stale.found, "FSTVEP09");
        }
        // and a magic that is not one of ours at all
        parts[0][..8].copy_from_slice(b"FSTVEP99");
        write_cache_of_parts(&path, &parts);
        for threads in [1, 4] {
            let said = format!("{:#}", load_on(threads, &path).unwrap_err());
            assert!(
                said.contains("wrong magic header"),
                "{threads} threads: {said}"
            );
        }
    }

    #[test]
    fn a_length_that_is_not_one_is_an_error_in_any_frame_and_not_an_allocation() {
        // A frame whose bytes are not what was written, with the checksum of the bytes it holds:
        // what a damaged byte inside the compressed data can decode to, and what a writer with a
        // bug makes. The id of the first transcript of a frame is made 1 GiB long, in the first
        // frame and in a later one: with no limit it is allocated, and with an exabyte the
        // process aborts.
        let transcripts = bulky_transcripts(60);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        save_small_frames(&transcripts, &path);
        let parts = parts_of(&path);
        assert!(parts.len() >= 5);
        // the magic, the length of the list, and then the length of the first id; in a later
        // frame the first transcript begins at its first byte
        for (frame, at) in [(0, 16), (3, 0)] {
            let mut damaged = parts.clone();
            damaged[frame][at..at + 8].copy_from_slice(&(1u64 << 30).to_le_bytes());
            write_cache_of_parts(&path, &damaged);
            for threads in [1, 4] {
                let said = format!("{:#}", load_on(threads, &path).unwrap_err());
                assert!(
                    said.contains("size limit"),
                    "frame {frame}, {threads} threads: {said}"
                );
            }
        }
    }

    #[test]
    fn a_count_that_the_frames_do_not_make_is_refused() {
        let transcripts = bulky_transcripts(60);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        save_small_frames(&transcripts, &path);
        let parts = parts_of(&path);
        // the cache says one more than it has, and one fewer
        for (count, frames_said, stream_said) in [
            (61u64, "the cache says it holds 61 transcripts", None),
            (
                59,
                "the cache says it holds 59 transcripts",
                Some("bytes after its 59 transcripts"),
            ),
        ] {
            let mut damaged = parts.clone();
            damaged[0][8..16].copy_from_slice(&count.to_le_bytes());
            write_cache_of_parts(&path, &damaged);
            let said = format!("{:#}", load_on(4, &path).unwrap_err());
            assert!(said.contains(frames_said), "{said}");
            // read as a stream, a list that says more than there is ends early, and one that says
            // less is followed by what it left out
            let error = load_on(1, &path).unwrap_err();
            match stream_said {
                Some(said) => assert!(format!("{error:#}").contains(said), "{error:#}"),
                None => assert!(is_a_stream_that_ended_early(&error), "{error:#}"),
            }
        }
    }

    #[test]
    fn a_byte_flipped_in_any_frame_of_a_cache_is_rejected_on_every_path() {
        let transcripts = bulky_transcripts(60);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        save_small_frames(&transcripts, &path);
        let bytes = std::fs::read(&path).unwrap();
        let frames: Vec<(usize, usize)> = frames_of(&path)
            .into_iter()
            .map(|(at, len)| (at as usize, len as usize))
            .collect();
        assert!(frames.len() >= 5);
        for (k, (at, len)) in frames.iter().enumerate() {
            for position in [at + len / 3, at + len - 1] {
                let mut damaged = bytes.clone();
                damaged[position] ^= 0xFF;
                std::fs::write(&path, &damaged).unwrap();
                for threads in [1, 4] {
                    assert!(
                        load_on(threads, &path).is_err(),
                        "frame {k}, byte {position}, {threads} threads: loaded"
                    );
                }
            }
        }
    }

    #[test]
    fn a_cache_missing_a_frame_or_cut_short_is_rejected_on_every_path() {
        let transcripts = bulky_transcripts(60);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.cache");
        save_small_frames(&transcripts, &path);
        let bytes = std::fs::read(&path).unwrap();
        let frames = frames_of(&path);
        assert!(frames.len() >= 5);
        let (at, len) = (frames[1].0 as usize, frames[1].1 as usize);

        // the second frame taken out, and the table left as it was: every frame that is left is
        // whole, and the table is of a file that is no longer this one
        let mut without = bytes[..at].to_vec();
        without.extend_from_slice(&bytes[at + len..]);
        std::fs::write(&path, &without).unwrap();
        let said = format!("{:#}", load_on(4, &path).unwrap_err());
        assert!(
            said.contains("the frames in the table of the cache are"),
            "{said}"
        );
        // and read as a stream, which does not look at the table, the list is short
        let error = load_on(1, &path).unwrap_err();
        assert!(is_a_stream_that_ended_early(&error), "{error:#}");

        // the second frame taken out of the cache, and a table made for the frames that are left:
        // nothing is wrong with the file but the list, which is short
        std::fs::write(&path, &bytes).unwrap();
        let mut parts = parts_of(&path);
        parts.remove(1);
        write_cache_of_parts(&path, &parts);
        let said = format!("{:#}", load_on(4, &path).unwrap_err());
        assert!(
            said.contains("the cache says it holds 60 transcripts and its"),
            "{said}"
        );
        let error = load_on(1, &path).unwrap_err();
        assert!(is_a_stream_that_ended_early(&error), "{error:#}");

        // the last frame cut short, which takes the table with it
        std::fs::write(&path, &bytes[..bytes.len() - 100]).unwrap();
        for threads in [1, 4] {
            assert!(load_on(threads, &path).is_err(), "{threads} threads");
        }
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

    fn write_v3_cache(path: &Path, transcripts: &[Transcript]) {
        use std::io::Write;
        let file = File::create(path).unwrap();
        let mut zst = zstd::Encoder::new(BufWriter::new(file), 1).unwrap();
        zst.write_all(CACHE_MAGIC_V3).unwrap();
        bincode::serialize_into(&mut zst, transcripts).unwrap();
        zst.finish().unwrap().flush().unwrap();
    }

    fn write_v9_cache(path: &Path, transcripts: &[Transcript]) {
        use std::io::Write;
        let file = File::create(path).unwrap();
        let mut zst = zstd::Encoder::new(BufWriter::new(file), 1).unwrap();
        zst.write_all(CACHE_MAGIC_V9).unwrap();
        bincode::serialize_into(&mut zst, transcripts).unwrap();
        zst.finish().unwrap().flush().unwrap();
    }

    fn write_v8_cache(path: &Path, transcripts: &[Transcript]) {
        use std::io::Write;
        let file = File::create(path).unwrap();
        let mut zst = zstd::Encoder::new(BufWriter::new(file), 1).unwrap();
        zst.write_all(CACHE_MAGIC_V8).unwrap();
        bincode::serialize_into(&mut zst, transcripts).unwrap();
        zst.finish().unwrap().flush().unwrap();
    }

    fn write_v7_cache(path: &Path, transcripts: &[Transcript]) {
        use std::io::Write;
        let file = File::create(path).unwrap();
        let mut zst = zstd::Encoder::new(BufWriter::new(file), 1).unwrap();
        zst.write_all(CACHE_MAGIC_V7).unwrap();
        bincode::serialize_into(&mut zst, transcripts).unwrap();
        zst.finish().unwrap().flush().unwrap();
    }

    fn write_v6_cache(path: &Path, transcripts: &[Transcript]) {
        use std::io::Write;
        let file = File::create(path).unwrap();
        let mut zst = zstd::Encoder::new(BufWriter::new(file), 1).unwrap();
        zst.write_all(CACHE_MAGIC_V6).unwrap();
        bincode::serialize_into(&mut zst, transcripts).unwrap();
        zst.finish().unwrap().flush().unwrap();
    }

    fn write_v5_cache(path: &Path, transcripts: &[Transcript]) {
        use std::io::Write;
        let file = File::create(path).unwrap();
        let mut zst = zstd::Encoder::new(BufWriter::new(file), 1).unwrap();
        zst.write_all(CACHE_MAGIC_V5).unwrap();
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
            msg.contains("intergenic"),
            "the error should say what goes wrong if it were read, got: {msg}"
        );
        // Recovery advice is the caller's to add: it differs between the
        // sidecar path, which has already rebuilt by the time it reports, and
        // the explicit path, where deleting the named file produces no sidecar.
        assert!(
            !msg.contains("fastvep cache"),
            "the type should not prescribe a recovery its callers disagree \
             about, got: {msg}"
        );
    }

    #[test]
    fn a_pre_98_cache_is_rejected_rather_than_read() {
        // A well-formed FSTVEP03 cache: whole-file, loads cleanly, transcript
        // count fine. What it cannot show is that every `ncRNA_gene` locus in it
        // lost its symbol to the pre-#98 parser. The sidecar path reuses a cache
        // that is newer than its GFF3, so without this rejection upgrading to the
        // fix would change nothing for anyone who had already run the tool.
        let tmp = NamedTempFile::new().unwrap();
        write_v3_cache(tmp.path(), &make_transcripts());

        let err = load_cache(tmp.path()).expect_err("FSTVEP03 must not load");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("FSTVEP03"),
            "the error should name the format found, got: {msg}"
        );
        assert!(
            msg.contains("SYMBOL"),
            "the error should say what goes wrong if it were read, got: {msg}"
        );
        assert!(
            !msg.contains("intergenic"),
            "V3 is whole-file; the #90 reason would send the user after the wrong \
             problem, got: {msg}"
        );
    }

    #[test]
    fn a_pre_protein_version_cache_is_rejected_rather_than_read() {
        // A well-formed FSTVEP05 cache: whole-file, loads cleanly, transcript count fine.
        // What it cannot show is that the parser which wrote it gave every protein
        // version 1, so every HGVSp it produces says `ENSP....1` for a protein that is
        // `.3` or `.12`. The sidecar path reuses a cache that is newer than its GFF3, so
        // without this rejection upgrading to the fix would change nothing for anyone who
        // had already run the tool.
        let tmp = NamedTempFile::new().unwrap();
        write_v5_cache(tmp.path(), &make_transcripts());

        let err = load_cache(tmp.path()).expect_err("FSTVEP05 must not load");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("FSTVEP05"),
            "the error should name the format found, got: {msg}"
        );
        assert!(
            msg.contains("protein") && msg.contains("version"),
            "the error should say what goes wrong if it were read, got: {msg}"
        );
        assert!(
            !msg.contains("intergenic"),
            "V5 is whole-file; the #90 reason would send the user after the wrong \
             problem, got: {msg}"
        );
    }

    #[test]
    fn a_pre_refseq_names_cache_is_rejected_rather_than_read() {
        // A well-formed FSTVEP06 cache: whole-file, loads cleanly, transcript count fine.
        // What it cannot show is that the parser which wrote it named a RefSeq transcript
        // `rna-NM_000001.3`, its gene `gene-TP53` and its biotype `mRNA`. The sidecar path
        // reuses a cache that is newer than its GFF3, so without this rejection upgrading
        // to the fix would change nothing for anyone who had already run the tool.
        let tmp = NamedTempFile::new().unwrap();
        write_v6_cache(tmp.path(), &make_transcripts());

        let err = load_cache(tmp.path()).expect_err("FSTVEP06 must not load");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("FSTVEP06"),
            "the error should name the format found, got: {msg}"
        );
        assert!(
            msg.contains("rna-NM_000001.3") && msg.contains("BIOTYPE"),
            "the error should say what goes wrong if it were read, got: {msg}"
        );
        assert!(
            !msg.contains("intergenic"),
            "V6 is whole-file; the #90 reason would send the user after the wrong \
             problem, got: {msg}"
        );
    }

    #[test]
    fn a_pre_cache_parity_cache_is_rejected_rather_than_read() {
        // A well-formed FSTVEP07 cache: whole-file, transcript count fine. What it cannot show is
        // that the parser which wrote it left a RefSeq gene with no symbol source, no canonical
        // transcript and no MANE Select, and an Ensembl transcript with no CDS flags. The sidecar
        // path reuses a cache that is newer than its GFF3, so without this rejection upgrading
        // would change nothing for anyone who had already run the tool.
        let tmp = NamedTempFile::new().unwrap();
        write_v7_cache(tmp.path(), &make_transcripts());

        let err = load_cache(tmp.path()).expect_err("FSTVEP07 must not load");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("FSTVEP07"),
            "the error should name the format found, got: {msg}"
        );
        assert!(
            msg.contains("EntrezGene") && msg.contains("canonical") && msg.contains("cds_start_NF"),
            "the error should say what goes wrong if it were read, got: {msg}"
        );
        assert!(
            !msg.contains("intergenic"),
            "V7 is whole-file; the #90 reason would send the user after the wrong \
             problem, got: {msg}"
        );
    }

    #[test]
    fn a_pre_organelle_cache_is_rejected_rather_than_read() {
        // A well-formed FSTVEP09 cache. What it cannot show is that the reader which wrote it
        // flagged the genes of an organelle contig incomplete by the standard genetic code. The
        // sidecar path reuses a cache that is newer than its GFF3, so without this rejection
        // upgrading would change nothing for anyone who had already run the tool on such a file.
        let tmp = NamedTempFile::new().unwrap();
        write_v9_cache(tmp.path(), &make_transcripts());

        let err = load_cache(tmp.path()).expect_err("FSTVEP09 must not load");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("FSTVEP09"),
            "the error should name the format found, got: {msg}"
        );
        assert!(
            msg.contains("Mito") && msg.contains("cds_start_NF"),
            "the error should say what goes wrong if it were read, got: {msg}"
        );
        assert!(
            !msg.contains("intergenic"),
            "V9 is whole-file; the #90 reason would send the user after the wrong \
             problem, got: {msg}"
        );
    }

    #[test]
    fn a_pre_gff_checks_cache_is_rejected_rather_than_read() {
        // A well-formed FSTVEP08 cache: whole-file, transcript count fine. What it cannot show is
        // that the parser which wrote it padded a codon with `u64::MAX` bases for a CDS with no
        // phase, or kept a transcript whose exons overlap. The sidecar path reuses a cache that
        // is newer than its GFF3, so without this rejection upgrading would change nothing for
        // anyone who had already run the tool on such a file.
        let tmp = NamedTempFile::new().unwrap();
        write_v8_cache(tmp.path(), &make_transcripts());

        let err = load_cache(tmp.path()).expect_err("FSTVEP08 must not load");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("FSTVEP08"),
            "the error should name the format found, got: {msg}"
        );
        assert!(
            msg.contains("no phase") && msg.contains("exons that overlap"),
            "the error should say what goes wrong if it were read, got: {msg}"
        );
        assert!(
            !msg.contains("intergenic"),
            "V8 is whole-file; the #90 reason would send the user after the wrong \
             problem, got: {msg}"
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
        assert!(
            msg.contains("#90"),
            "it should carry the pre-#90 reason, got: {msg}"
        );
    }

    #[test]
    fn the_summary_is_one_sentence_and_the_detail_is_only_in_the_long_form() {
        // The split exists so the sidecar path, which recovers by itself, can
        // report in one line while the explicit path still owes the user the
        // evidence. If `summary` ever grows into the paragraph, the sidecar
        // warning silently becomes one again.
        let tmp = NamedTempFile::new().unwrap();
        write_v3_cache(tmp.path(), &make_transcripts());
        let err = load_cache(tmp.path()).expect_err("FSTVEP03 must not load");
        let stale = err
            .downcast_ref::<StaleCacheFormat>()
            .expect("a rejected format should downcast to StaleCacheFormat");

        let summary = stale.summary();
        assert!(
            summary.contains("FSTVEP03") && summary.contains("#98"),
            "the one-liner must still say which format and why, got: {summary}"
        );
        assert!(
            !summary.contains('.'),
            "the one-liner is a clause, not a paragraph; the caller punctuates \
             it, got: {summary}"
        );
        assert!(
            summary.starts_with("is in the "),
            "the one-liner must leave its subject to the caller, which puts the \
             file in front of it, got: {summary}"
        );
        assert!(
            !summary.contains("SYMBOL"),
            "the evidence belongs in `detail`, got: {summary}"
        );
        assert!(
            stale.to_string().contains("SYMBOL"),
            "the long form must still carry the evidence"
        );
    }

    #[test]
    fn a_v4_cache_written_now_round_trips() {
        // The other half of the bump: rejecting V2 and V3 is only acceptable if
        // the current writer produces something the current reader accepts.
        let tmp = NamedTempFile::new().unwrap();
        let transcripts = make_transcripts();
        save_cache(&transcripts, tmp.path()).unwrap();

        let loaded = load_cache(tmp.path()).unwrap();
        assert_eq!(loaded.len(), transcripts.len());
        assert_eq!(&*loaded[0].stable_id, "ENST00000001");
    }

    #[test]
    fn the_published_magic_is_v10() {
        // Pin the on-disk byte, so a future change to the writer that forgets
        // the reader shows up here rather than in someone's annotation.
        use std::io::Read;
        let tmp = NamedTempFile::new().unwrap();
        save_cache(&make_transcripts(), tmp.path()).unwrap();

        let mut zst = zstd::Decoder::new(BufReader::new(File::open(tmp.path()).unwrap())).unwrap();
        let mut magic = [0u8; 8];
        zst.read_exact(&mut magic).unwrap();
        assert_eq!(&magic, CACHE_MAGIC_V10, "published magic drifted from V10");
    }
}
