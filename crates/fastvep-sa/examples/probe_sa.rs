//! Probe a .osa file to print block sizes (compressed + decompressed).
//! Usage: cargo run --release -p fastvep-sa --example probe_sa -- <path/to/file.osa>

use anyhow::Result;
use fastvep_sa::block::SaBlock;
use fastvep_sa::index::SaIndex;
use memmap2::Mmap;
use std::env;
use std::fs::File;
use std::path::PathBuf;

fn main() -> Result<()> {
    let path: PathBuf = env::args().nth(1).expect("usage: probe_sa <file.osa>").into();
    let idx_path = path.with_extension("osa.idx");
    let mut f = File::open(&idx_path)?;
    let idx = SaIndex::read_from(&mut f)?;
    let data_file = File::open(&path)?;
    let mmap = unsafe { Mmap::map(&data_file)? };

    let mut total_comp = 0u64;
    let mut total_decomp = 0u64;
    let mut count = 0u64;
    let mut max_decomp = 0usize;
    let mut decomp_struct_overhead = 0u64;

    for (_chrom, blocks) in &idx.chromosomes {
        for br in blocks {
            count += 1;
            let off = br.file_offset as usize;
            let comp_len = br.compressed_len as usize;
            // skip the 4-byte length prefix
            let data = &mmap[off + 4..off + 4 + comp_len];
            let entries = SaBlock::decompress(data)?;
            total_comp += comp_len as u64;
            // raw bytes: u32 + 3 strings' data
            let raw: usize = entries.iter().map(|e| 4 + e.ref_allele.len() + e.alt_allele.len() + e.json.len()).sum();
            total_decomp += raw as u64;
            // struct overhead: Vec header (24) + per-entry (size_of::<BlockEntry>() = ~76)
            let overhead = 24 + entries.len() * std::mem::size_of::<fastvep_sa::block::BlockEntry>();
            decomp_struct_overhead += overhead as u64;
            let block_inmem = raw + overhead;
            max_decomp = max_decomp.max(block_inmem);
        }
    }
    println!("blocks: {}", count);
    println!("avg compressed: {} KB", total_comp / count / 1024);
    println!("avg raw bytes (strings only): {} KB", total_decomp / count / 1024);
    println!("avg struct overhead: {} KB", decomp_struct_overhead / count / 1024);
    println!("total compressed: {} MB", total_comp / (1024 * 1024));
    println!("total decompressed (strings+overhead): {} MB", (total_decomp + decomp_struct_overhead) / (1024 * 1024));
    println!("max single in-mem block: {} MB", max_decomp / (1024 * 1024));
    println!("with cap=4: max footprint ≈ {} MB", (max_decomp * 4) / (1024 * 1024));
    Ok(())
}
