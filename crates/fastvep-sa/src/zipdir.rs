//! Minimal ZIP central-directory reader over a memory-mapped archive.
//!
//! **Why this exists (issue #78).** `.osa2` files are ZIP archives whose entry
//! table the reader needs at `open`. Walking that table with
//! `zip::ZipArchive::by_index` costs one *seek + read of the entry's local file
//! header* per entry, because the crate resolves each entry's data offset
//! eagerly. Local headers sit next to their data, so those reads are scattered
//! across the whole (multi-GB) file: a per-chromosome gnomAD shard has a few
//! thousand entries, and a 24-shard `--sa-dir` therefore paid ~100k random
//! reads before annotation could start - 25+ minutes on rotational/network
//! storage, and invisible because it happens before the progress meter.
//!
//! The central directory is a single contiguous region at the end of the file
//! and holds everything the reader needs except each entry's local-header
//! padding. Parsing it here, straight from the mmap, makes `open` a sequential
//! read of that one region; the local header is then resolved lazily, on the
//! first read of that entry, when its page is about to be touched anyway.
//!
//! Only what `.osa2` needs is implemented: no encryption, no multi-disk
//! archives, no data descriptors. ZIP64 is supported because gnomAD-scale
//! shards cross the 4 GiB offset boundary.

use anyhow::{bail, Context, Result};

const EOCD_SIG: u32 = 0x0605_4b50;
const ZIP64_EOCD_LOCATOR_SIG: u32 = 0x0706_4b50;
const ZIP64_EOCD_SIG: u32 = 0x0606_4b50;
const CENTRAL_HEADER_SIG: u32 = 0x0201_4b50;
pub(crate) const LOCAL_HEADER_SIG: u32 = 0x0403_4b50;

const EOCD_FIXED_LEN: usize = 22;
const ZIP64_EOCD_LOCATOR_LEN: usize = 20;
const CENTRAL_HEADER_FIXED_LEN: usize = 46;
pub(crate) const LOCAL_HEADER_FIXED_LEN: u64 = 30;

/// The ZIP32 "value lives in the ZIP64 extra field" sentinel.
const U32_SENTINEL: u32 = u32::MAX;
const U16_SENTINEL: u16 = u16::MAX;

/// One entry as described by the central directory. `data_start` is
/// deliberately absent: deriving it requires the entry's *local* header, which
/// is the random read this module exists to avoid.
#[derive(Debug)]
pub(crate) struct CentralEntry {
    pub name: String,
    pub header_start: u64,
    pub comp_size: u64,
    /// Raw ZIP compression method id (0 = stored, 8 = deflated).
    pub method: u16,
}

/// Fixed-width little-endian read with an explicit bounds check. Offsets here
/// can come from the file's own (possibly corrupt) fields, so the addition is
/// checked rather than trusted.
fn read_le<const N: usize>(buf: &[u8], at: usize) -> Result<[u8; N]> {
    let end = at
        .checked_add(N)
        .filter(|&end| end <= buf.len())
        .ok_or_else(|| anyhow::anyhow!("truncated ZIP structure at byte {}", at))?;
    let mut b = [0u8; N];
    b.copy_from_slice(&buf[at..end]);
    Ok(b)
}

fn read_u16(buf: &[u8], at: usize) -> Result<u16> {
    read_le::<2>(buf, at).map(u16::from_le_bytes)
}

fn read_u32(buf: &[u8], at: usize) -> Result<u32> {
    read_le::<4>(buf, at).map(u32::from_le_bytes)
}

fn read_u64(buf: &[u8], at: usize) -> Result<u64> {
    read_le::<8>(buf, at).map(u64::from_le_bytes)
}

/// Locate the end-of-central-directory record. The ZIP comment may be up to
/// 64 KiB, so the record can sit that far from the end; scan backwards and
/// accept the first candidate whose comment length lands exactly on EOF, which
/// rules out a stray signature inside compressed data.
fn find_eocd(mmap: &[u8]) -> Result<usize> {
    if mmap.len() < EOCD_FIXED_LEN {
        bail!(
            "file is too small to be a ZIP archive ({} bytes)",
            mmap.len()
        );
    }
    let max_back = (u16::MAX as usize + EOCD_FIXED_LEN).min(mmap.len());
    let search_start = mmap.len() - max_back;
    for pos in (search_start..=mmap.len() - EOCD_FIXED_LEN).rev() {
        if read_u32(mmap, pos)? != EOCD_SIG {
            continue;
        }
        let comment_len = read_u16(mmap, pos + 20)? as usize;
        if pos + EOCD_FIXED_LEN + comment_len == mmap.len() {
            return Ok(pos);
        }
    }
    bail!("not a ZIP archive: no end-of-central-directory record found")
}

/// Where the central directory lives and how many records it holds.
struct CentralDirectoryLocator {
    offset: u64,
    size: u64,
    count: u64,
    /// Whether `count` is the true record count. The ZIP32 footer clamps the
    /// count to 65535, so an archive with exactly that many entries reports a
    /// count that happens to equal the sentinel without being one.
    count_is_exact: bool,
    /// Byte position of whichever record terminates the central directory
    /// (the ZIP64 EOCD if present, else the EOCD). Used to recover the real
    /// directory position when the archive has prepended data.
    directory_end: u64,
}

fn locate_central_directory(mmap: &[u8], eocd: usize) -> Result<CentralDirectoryLocator> {
    let count32 = read_u16(mmap, eocd + 10)?;
    let size32 = read_u32(mmap, eocd + 12)?;
    let offset32 = read_u32(mmap, eocd + 16)?;

    // The ZIP64 footer, when present, is authoritative - probe for it rather
    // than inferring from the ZIP32 fields. Writers clamp all three of those
    // to their maximum whether or not ZIP64 is in play, so a clamped value is
    // not by itself proof of anything: an archive with exactly 65535 entries
    // reports count32 == 0xFFFF and has no ZIP64 footer at all.
    if let Some(z64) = find_zip64_eocd(mmap, eocd)? {
        return Ok(CentralDirectoryLocator {
            offset: read_u64(mmap, z64 + 48)?,
            size: read_u64(mmap, z64 + 40)?,
            count: read_u64(mmap, z64 + 32)?,
            count_is_exact: true,
            directory_end: z64 as u64,
        });
    }

    // No ZIP64 footer, so the ZIP32 fields must be real. A clamped size or
    // offset here means the archive genuinely needed ZIP64 and does not have
    // it; a clamped count is only ambiguous, and `count` is advisory anyway.
    if size32 == U32_SENTINEL || offset32 == U32_SENTINEL {
        bail!("archive needs a ZIP64 end-of-central-directory record but has none");
    }
    Ok(CentralDirectoryLocator {
        offset: offset32 as u64,
        size: size32 as u64,
        count: count32 as u64,
        count_is_exact: count32 != U16_SENTINEL,
        directory_end: eocd as u64,
    })
}

/// Offset of the ZIP64 end-of-central-directory record, if this archive has
/// one. The locator sits immediately before the ZIP32 EOCD and points at it.
fn find_zip64_eocd(mmap: &[u8], eocd: usize) -> Result<Option<usize>> {
    let Some(locator) = eocd.checked_sub(ZIP64_EOCD_LOCATOR_LEN) else {
        return Ok(None);
    };
    if read_u32(mmap, locator)? != ZIP64_EOCD_LOCATOR_SIG {
        return Ok(None);
    }
    let z64 = read_u64(mmap, locator + 8)?;
    let z64: usize = z64
        .try_into()
        .map_err(|_| anyhow::anyhow!("ZIP64 EOCD offset {} exceeds usize", z64))?;
    if z64 >= mmap.len() || read_u32(mmap, z64)? != ZIP64_EOCD_SIG {
        bail!(
            "ZIP64 end-of-central-directory record not found at offset {}",
            z64
        );
    }
    Ok(Some(z64))
}

/// Byte range of the central directory within `mmap`.
///
/// The recorded offset is relative to the start of the ZIP data, which is not
/// the start of the file for an archive with prepended bytes (a
/// self-extracting stub, or anything concatenated in front). When the recorded
/// offset does not land on a directory record, fall back to deriving the
/// position from where the directory ends - the same recovery the `zip` crate
/// performs, kept so such archives keep opening.
fn central_directory_range(mmap: &[u8], loc: &CentralDirectoryLocator) -> Result<(usize, usize)> {
    let size: usize = loc
        .size
        .try_into()
        .map_err(|_| anyhow::anyhow!("central directory size {} exceeds usize", loc.size))?;

    let starts_a_record = |start: u64| -> bool {
        let Ok(start) = usize::try_from(start) else {
            return false;
        };
        match start.checked_add(size) {
            Some(end) if end <= mmap.len() => {
                size == 0
                    || read_u32(mmap, start)
                        .map(|s| s == CENTRAL_HEADER_SIG)
                        .unwrap_or(false)
            }
            _ => false,
        }
    };

    let candidates = [Some(loc.offset), loc.directory_end.checked_sub(loc.size)];
    for start in candidates.into_iter().flatten() {
        if starts_a_record(start) {
            let start = start as usize;
            return Ok((start, start + size));
        }
    }
    bail!(
        "central directory not found at offset {} ({} bytes) in a {}-byte archive",
        loc.offset,
        loc.size,
        mmap.len()
    )
}

/// Parse every central-directory record in `mmap`.
///
/// Sequential over one contiguous region: no random access, no syscalls, and
/// no local file headers touched.
pub(crate) fn parse_central_directory(mmap: &[u8]) -> Result<Vec<CentralEntry>> {
    let eocd = find_eocd(mmap)?;
    let loc = locate_central_directory(mmap, eocd)?;
    let (start, end) = central_directory_range(mmap, &loc)?;
    let cd = &mmap[start..end];

    // Every offset the directory records - its own and each local header's -
    // is relative to the start of the ZIP data. Shift them all by however far
    // that is from the start of the file (zero for anything this crate wrote).
    let archive_offset = (start as u64).saturating_sub(loc.offset);

    // `count` comes from the file, so it only pre-sizes the Vec - the loop
    // below is bounded by the directory's byte length, not by this number.
    let mut entries = Vec::with_capacity((loc.count as usize).min(1 << 16));
    let mut at = 0usize;
    while at + CENTRAL_HEADER_FIXED_LEN <= cd.len() {
        if read_u32(cd, at)? != CENTRAL_HEADER_SIG {
            bail!(
                "bad central-directory record signature at offset {}",
                start + at
            );
        }
        let method = read_u16(cd, at + 10)?;
        let comp_size32 = read_u32(cd, at + 20)?;
        let uncomp_size32 = read_u32(cd, at + 24)?;
        let name_len = read_u16(cd, at + 28)? as usize;
        let extra_len = read_u16(cd, at + 30)? as usize;
        let comment_len = read_u16(cd, at + 32)? as usize;
        let header_start32 = read_u32(cd, at + 42)?;

        let name_at = at + CENTRAL_HEADER_FIXED_LEN;
        let extra_at = name_at + name_len;
        let record_end = extra_at + extra_len + comment_len;
        if record_end > cd.len() {
            bail!(
                "central-directory record at offset {} is truncated",
                start + at
            );
        }

        let name = std::str::from_utf8(&cd[name_at..extra_at])
            .map_err(|_| {
                anyhow::anyhow!(
                    "central-directory entry at offset {} has a non-UTF-8 name",
                    start + at
                )
            })?
            .to_string();

        let (comp_size, header_start) = resolve_zip64(
            &cd[extra_at..extra_at + extra_len],
            uncomp_size32,
            comp_size32,
            header_start32,
        )
        .with_context(|| format!("central-directory entry '{}'", name))?;

        entries.push(CentralEntry {
            name,
            header_start: header_start.saturating_add(archive_offset),
            comp_size,
            method,
        });
        at = record_end;
    }

    // A short or otherwise wrong directory would silently truncate the entry
    // table, and a missing entry reads back as "this chunk has no data for
    // this field" - annotations would quietly disappear instead of the load
    // failing. Cross-check against the count the footer declares.
    if loc.count_is_exact && entries.len() as u64 != loc.count {
        bail!(
            "central directory declares {} entries but {} were readable in its {} bytes",
            loc.count,
            entries.len(),
            end - start
        );
    }
    if at != cd.len() {
        bail!(
            "central directory has {} trailing bytes after its last record",
            cd.len() - at
        );
    }

    Ok(entries)
}

/// Pull the full-width compressed size / local-header offset out of the ZIP64
/// extended-information extra field.
///
/// The field packs values in a fixed order - uncompressed size, compressed
/// size, local header offset, disk number - but which of them are present is
/// not stated anywhere in the field itself. The spec says only the values that
/// overflowed their ZIP32 slot appear, yet `zip`'s own writer emits all three
/// u64s whenever `large_file` is set, sentinel or not. Reading purely by
/// sentinel therefore mis-aligns on those archives and hands back an
/// uncompressed size where a header offset belongs. Mirror `zip`'s reader
/// exactly: a body of 24 bytes or more carries all three regardless.
fn resolve_zip64(
    extra: &[u8],
    uncomp_size32: u32,
    comp_size32: u32,
    header_start32: u32,
) -> Result<(u64, u64)> {
    let mut comp_size = comp_size32 as u64;
    let mut header_start = header_start32 as u64;
    if comp_size32 != U32_SENTINEL && header_start32 != U32_SENTINEL {
        return Ok((comp_size, header_start));
    }

    let mut at = 0usize;
    while at + 4 <= extra.len() {
        let id = read_u16(extra, at)?;
        let len = read_u16(extra, at + 2)? as usize;
        let body_at = at + 4;
        let body_end = body_at
            .checked_add(len)
            .filter(|&e| e <= extra.len())
            .ok_or_else(|| anyhow::anyhow!("ZIP extra field at offset {} is truncated", at))?;
        if id == 0x0001 {
            let body = &extra[body_at..body_end];
            let all_present = len >= 24;
            let mut cursor = 0usize;
            if all_present || uncomp_size32 == U32_SENTINEL {
                read_u64(body, cursor)?;
                cursor += 8;
            }
            if all_present || comp_size32 == U32_SENTINEL {
                comp_size = read_u64(body, cursor)?;
                cursor += 8;
            }
            if all_present || header_start32 == U32_SENTINEL {
                header_start = read_u64(body, cursor)?;
            }
            return Ok((comp_size, header_start));
        }
        at = body_end;
    }

    bail!("entry needs a ZIP64 extended-information field but none is present")
}

/// Offset of an entry's data, derived from its *local* file header.
///
/// Separate from the central directory on purpose: the local header's name and
/// extra-field lengths can differ from the central copy, so this is the one
/// piece that cannot be answered without touching the entry itself. Callers
/// resolve it lazily, on first read of the entry.
pub(crate) fn local_data_start(mmap: &[u8], header_start: u64) -> Result<u64> {
    let h: usize = header_start
        .try_into()
        .map_err(|_| anyhow::anyhow!("ZIP local header offset {} exceeds usize", header_start))?;
    let fixed_end = h
        .checked_add(LOCAL_HEADER_FIXED_LEN as usize)
        .ok_or_else(|| anyhow::anyhow!("ZIP local header offset overflow"))?;
    if fixed_end > mmap.len() {
        bail!("ZIP local header at {} extends beyond file", header_start);
    }
    if read_u32(mmap, h)? != LOCAL_HEADER_SIG {
        bail!("bad ZIP local header signature at {}", header_start);
    }
    let name_len = read_u16(mmap, h + 26)? as u64;
    let extra_len = read_u16(mmap, h + 28)? as u64;
    Ok(header_start + LOCAL_HEADER_FIXED_LEN + name_len + extra_len)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Cursor, Write};

    /// Build a small in-memory ZIP with the same writer `.osa2` uses, then
    /// check the hand-rolled parser agrees with the `zip` crate on every entry.
    fn round_trip(names_and_bodies: &[(&str, Vec<u8>)], large: bool) -> Vec<u8> {
        let mut buf = Vec::new();
        {
            let mut zw = zip::ZipWriter::new(Cursor::new(&mut buf));
            let opts = zip::write::SimpleFileOptions::default()
                .compression_method(zip::CompressionMethod::Deflated)
                .large_file(large);
            for (name, body) in names_and_bodies {
                zw.start_file(*name, opts).unwrap();
                zw.write_all(body).unwrap();
            }
            zw.finish().unwrap();
        }
        buf
    }

    fn assert_parity(buf: &[u8]) {
        let parsed = parse_central_directory(buf).unwrap();
        let mut archive = zip::ZipArchive::new(Cursor::new(buf)).unwrap();
        assert_eq!(parsed.len(), archive.len(), "entry count");
        for (i, entry) in parsed.iter().enumerate() {
            let zf = archive.by_index(i).unwrap();
            assert_eq!(entry.name, zf.name(), "name at {}", i);
            assert_eq!(
                entry.comp_size,
                zf.compressed_size(),
                "comp_size of {}",
                entry.name
            );
            assert_eq!(
                entry.header_start,
                zf.header_start(),
                "header_start of {}",
                entry.name
            );
            let expected_method = match zf.compression() {
                zip::CompressionMethod::Stored => 0u16,
                zip::CompressionMethod::Deflated => 8u16,
                other => panic!("unexpected compression {other:?}"),
            };
            assert_eq!(entry.method, expected_method, "method of {}", entry.name);
            let expected_start = zf.data_start();
            drop(zf);
            assert_eq!(
                local_data_start(buf, entry.header_start).unwrap(),
                expected_start,
                "data_start of {}",
                entry.name
            );
        }
    }

    #[test]
    fn parses_a_plain_archive() {
        let buf = round_trip(
            &[
                ("fastsa/metadata.json", b"{\"a\":1}".to_vec()),
                ("fastsa/chr1/0/var32.bin", vec![7u8; 4096]),
                ("fastsa/chr1/0/allAf.bin", vec![9u8; 1024]),
            ],
            false,
        );
        assert_parity(&buf);
    }

    #[test]
    fn parses_zip64_entries() {
        // `large_file(true)` forces ZIP64 extra fields on every entry, which is
        // the layout a multi-GB gnomAD shard produces.
        let buf = round_trip(
            &[
                ("fastsa/metadata.json", b"{}".to_vec()),
                ("fastsa/chr1/0/var32.bin", vec![3u8; 8192]),
            ],
            true,
        );
        assert_parity(&buf);
    }

    #[test]
    fn parses_archive_with_trailing_comment() {
        let mut buf = Vec::new();
        {
            let mut zw = zip::ZipWriter::new(Cursor::new(&mut buf));
            zw.set_comment("a comment that follows the EOCD record");
            zw.start_file(
                "fastsa/metadata.json",
                zip::write::SimpleFileOptions::default(),
            )
            .unwrap();
            zw.write_all(b"{}").unwrap();
            zw.finish().unwrap();
        }
        assert_parity(&buf);
    }

    /// An archive with bytes glued in front of it records offsets relative to
    /// the ZIP data, not the file. The `zip` crate recovers from this, so the
    /// hand-rolled parser must too or such files would stop opening.
    #[test]
    fn parses_an_archive_with_prepended_data() {
        let inner = round_trip(
            &[
                ("fastsa/metadata.json", b"{}".to_vec()),
                ("fastsa/chr1/0/var32.bin", vec![5u8; 2048]),
            ],
            false,
        );
        let mut buf = vec![0xAAu8; 4096];
        buf.extend_from_slice(&inner);

        let parsed = parse_central_directory(&buf).unwrap();
        assert_eq!(parsed.len(), 2);
        for entry in &parsed {
            // Offsets must be file-relative, so each one lands on a real local
            // header past the prepended block.
            assert!(entry.header_start >= 4096, "{} not shifted", entry.name);
            local_data_start(&buf, entry.header_start).unwrap();
        }
    }

    /// The ZIP32 footer clamps its entry count to 65535, so an archive with
    /// exactly that many entries reports the sentinel while having no ZIP64
    /// footer at all. Reading the clamp as "this is ZIP64" rejected the file -
    /// and `load_sa_providers` downgrades an open failure to a warning, so the
    /// shard would have been dropped silently.
    #[test]
    fn parses_archives_around_the_zip32_entry_clamp() {
        for n in [65534usize, 65535, 65536] {
            let mut buf = Vec::new();
            {
                let mut zw = zip::ZipWriter::new(Cursor::new(&mut buf));
                let opts = zip::write::SimpleFileOptions::default()
                    .compression_method(zip::CompressionMethod::Stored);
                for i in 0..n {
                    zw.start_file(format!("e{i}"), opts).unwrap();
                }
                zw.finish().unwrap();
            }
            let parsed = parse_central_directory(&buf)
                .unwrap_or_else(|e| panic!("{n}-entry archive failed to parse: {e}"));
            assert_eq!(parsed.len(), n, "entry count for a {n}-entry archive");
        }
    }

    /// `zip`'s writer emits all three ZIP64 u64s whenever `large_file` is set,
    /// even for values that fit in 32 bits. Picking fields by sentinel alone
    /// mis-aligns on those and returns the uncompressed size as the header
    /// offset.
    #[test]
    fn zip64_extra_field_of_full_length_holds_all_three_values() {
        let mut extra = Vec::new();
        extra.extend_from_slice(&7u64.to_le_bytes()); // uncompressed size
        extra.extend_from_slice(&11u64.to_le_bytes()); // compressed size
        extra.extend_from_slice(&0x1_0000_0000u64.to_le_bytes()); // header start
        let field = {
            let mut f = vec![0x01, 0x00, 24, 0x00];
            f.extend_from_slice(&extra);
            f
        };
        // Sizes fit in 32 bits; only the header offset is a sentinel.
        let (comp_size, header_start) = resolve_zip64(&field, 7, 11, U32_SENTINEL).unwrap();
        assert_eq!(comp_size, 11);
        assert_eq!(header_start, 0x1_0000_0000);
    }

    /// The spec-minimal form: only the overflowed value is present.
    #[test]
    fn zip64_extra_field_of_minimal_length_holds_only_the_overflowed_value() {
        let mut field = vec![0x01, 0x00, 8, 0x00];
        field.extend_from_slice(&0x2_0000_0000u64.to_le_bytes());
        let (comp_size, header_start) = resolve_zip64(&field, 7, 11, U32_SENTINEL).unwrap();
        assert_eq!(comp_size, 11);
        assert_eq!(header_start, 0x2_0000_0000);
    }

    /// A directory that parses into fewer records than the footer declares
    /// must fail the load, not yield a short entry table whose missing entries
    /// read back as "this field has no data".
    #[test]
    fn rejects_a_directory_shorter_than_its_declared_count() {
        let mut buf = round_trip(
            &[
                ("a.bin", vec![1u8; 64]),
                ("b.bin", vec![2u8; 64]),
                ("c.bin", vec![3u8; 64]),
            ],
            false,
        );
        assert_eq!(parse_central_directory(&buf).unwrap().len(), 3);

        // Shrink the recorded directory size so the last record falls outside.
        let eocd = find_eocd(&buf).unwrap();
        let size = read_u32(&buf, eocd + 12).unwrap();
        buf[eocd + 12..eocd + 16].copy_from_slice(&(size / 2).to_le_bytes());
        let err = parse_central_directory(&buf).unwrap_err().to_string();
        assert!(
            err.contains("declares 3 entries") || err.contains("trailing bytes"),
            "{err}"
        );
    }

    #[test]
    fn rejects_a_non_zip_file() {
        let err = parse_central_directory(&[0u8; 512])
            .unwrap_err()
            .to_string();
        assert!(err.contains("no end-of-central-directory"), "{}", err);
    }

    #[test]
    fn rejects_a_truncated_file() {
        let err = parse_central_directory(b"PK").unwrap_err().to_string();
        assert!(err.contains("too small"), "{}", err);
    }

    #[test]
    fn rejects_a_corrupt_central_directory() {
        let mut buf = round_trip(&[("a.bin", vec![1u8; 64])], false);
        // Corrupt the first central-directory record signature.
        let eocd = find_eocd(&buf).unwrap();
        let cd_off = read_u32(&buf, eocd + 16).unwrap() as usize;
        buf[cd_off] = 0x00;
        assert!(parse_central_directory(&buf).is_err());
    }
}
