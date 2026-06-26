//! Dump records from a .osa file to stdout as tab-separated values.
//! Usage: cargo run --release -p fastvep-sa --example dump_sa -- <file.osa> [max_records]
//!
//! If the .osa.idx index file does not exist (e.g. the build is still in
//! progress), the file is scanned sequentially and chromosome names are
//! reported as "?".

use anyhow::Result;
use fastvep_sa::block::SaBlock;
use fastvep_sa::index::SaIndex;
use memmap2::Mmap;
use std::env;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;

fn main() -> Result<()> {
    let mut args = env::args().skip(1);
    let path: PathBuf = args
        .next()
        .expect("usage: dump_sa <file.osa> [max_records]")
        .into();
    let max_records: Option<u64> = args.next().and_then(|s| s.parse().ok());

    let idx_path = path.with_extension("osa.idx");

    println!("chrom\tpos\tref\talt\tjson");

    if idx_path.exists() {
        dump_with_index(&path, &idx_path, max_records)
    } else {
        eprintln!("note: no index found — scanning sequentially, chromosome names will be '?'");
        dump_sequential(&path, max_records)
    }
}

fn dump_with_index(path: &PathBuf, idx_path: &PathBuf, max_records: Option<u64>) -> Result<()> {
    let mut idx_file = File::open(idx_path)?;
    let idx = SaIndex::read_from(&mut idx_file)?;
    let data_file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&data_file)? };

    let mut chroms: Vec<&String> = idx.chromosomes.keys().collect();
    chroms.sort_by(|a, b| chrom_sort_key(a).cmp(&chrom_sort_key(b)));

    let mut printed = 0u64;
    'outer: for chrom in chroms {
        for br in &idx.chromosomes[chrom] {
            let off = br.file_offset as usize;
            let data = &mmap[off + 4..off + 4 + br.compressed_len as usize];
            let entries = SaBlock::decompress(data)?;
            for entry in entries {
                println!(
                    "{}\t{}\t{}\t{}\t{}",
                    chrom, entry.position, entry.ref_allele, entry.alt_allele, entry.json
                );
                printed += 1;
                if max_records.map_or(false, |m| printed >= m) {
                    break 'outer;
                }
            }
        }
    }
    Ok(())
}

fn dump_sequential(path: &PathBuf, max_records: Option<u64>) -> Result<()> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);

    // Skip 8-byte magic + 2-byte schema version
    let mut header = [0u8; 10];
    reader.read_exact(&mut header)?;

    let mut printed = 0u64;
    let mut len_buf = [0u8; 4];
    loop {
        match reader.read_exact(&mut len_buf) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
            Err(e) => return Err(e.into()),
        }
        let compressed_len = u32::from_le_bytes(len_buf) as usize;
        let mut compressed = vec![0u8; compressed_len];
        reader.read_exact(&mut compressed)?;

        let entries = SaBlock::decompress(&compressed)?;
        for entry in entries {
            println!(
                "?\t{}\t{}\t{}\t{}",
                entry.position, entry.ref_allele, entry.alt_allele, entry.json
            );
            printed += 1;
            if max_records.map_or(false, |m| printed >= m) {
                return Ok(());
            }
        }
    }
    Ok(())
}

fn chrom_sort_key(c: &str) -> (u8, u32, &str) {
    let name = c.strip_prefix("chr").unwrap_or(c);
    if let Ok(n) = name.parse::<u32>() {
        (0, n, "")
    } else {
        let order: u8 = match name {
            "X" => 1,
            "Y" => 2,
            "M" | "MT" => 3,
            _ => 4,
        };
        (order, 0, name)
    }
}
