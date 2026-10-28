//! The same RefSeq GFF3 gives a transcript the same `CANONICAL` flag whether it is read as a plain
//! file or through tabix.
//!
//! RefSeq's GFF3 has no `Ensembl_canonical` tag, so fastVEP works the flag out: the MANE Select
//! transcript of a gene, else the transcript that spans the most of the genome, else the first in
//! the file among equals (`derive_refseq_canonical`). That needs every transcript of the gene. A
//! tabix read has the ones that overlap its regions, and used to judge none of them but a MANE
//! Select one, so a run over the sorted `.gff3.gz` of chromosomes 21 and 22 of NCBI's GRCh38 file
//! gave 4,266 of 212,240 rows no `CANONICAL` that the same run over the same lines as a plain file
//! gave `YES` (38 transcripts, and with them `--pick` and `--canonical` of anything that reads the
//! flag).
//!
//! The genes here are in separate 16 kb windows of the tabix index, so that a region reads one
//! transcript of a gene and not the other, and each test says which it expects to have read: a
//! layout that read everything would pass without trying anything.

use fastvep_cache::gff::{parse_gff3, parse_gff3_indexed};
use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::binning_index::index::header::format::CoordinateSystem;
use noodles_csi::binning_index::index::header::Format;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_csi::binning_index::index::Header;
use noodles_csi::binning_index::Indexer;
use noodles_tabix as tabix;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

/// A RefSeq gene (`Dbxref=GeneID:`), which is what makes the flag be worked out.
fn gene(id: &str, geneid: u32, start: u64, end: u64) -> String {
    format!(
        "1\tBestRefSeq\tgene\t{start}\t{end}\t.\t+\t.\tID=gene-{id};Dbxref=GeneID:{geneid};Name={id};gbkey=Gene;gene={id};gene_biotype=protein_coding"
    )
}

/// A transcript of it and its one exon, which is the transcript.
fn transcript(
    accession: &str,
    gene: &str,
    geneid: u32,
    start: u64,
    end: u64,
    mane: bool,
) -> Vec<String> {
    let tag = if mane { ";tag=MANE Select" } else { "" };
    let common = format!(
        "Dbxref=GeneID:{geneid},GenBank:{accession};gbkey=mRNA;gene={gene};transcript_id={accession}"
    );
    vec![
        format!(
            "1\tBestRefSeq\tmRNA\t{start}\t{end}\t.\t+\t.\tID=rna-{accession};Parent=gene-{gene};{common};Name={accession}{tag}"
        ),
        format!(
            "1\tBestRefSeq\texon\t{start}\t{end}\t.\t+\t.\tID=exon-{accession}-1;Parent=rna-{accession};{common}"
        ),
    ]
}

/// The lines, in the order `sort -s -k1,1 -k4,4n` leaves them: by start, and the ones with the same
/// start in the order they were written. (Not the order of a plain `sort -k1,1 -k4,4n`, which puts
/// equal starts in the order of their text; see `docs/VEP_DIVERGENCE.md`.)
fn world() -> Vec<String> {
    let mut lines = Vec::new();
    // MANE Select settles it, though the other transcript is first in the file
    lines.push(gene("MANEGENE", 101, 1_000, 9_000));
    lines.extend(transcript(
        "NM_000002.1",
        "MANEGENE",
        101,
        1_000,
        9_000,
        false,
    ));
    lines.extend(transcript(
        "NM_000001.1",
        "MANEGENE",
        101,
        1_000,
        9_000,
        true,
    ));
    // no MANE Select: the longest span, NM_000004.1, which is neither first nor last
    lines.push(gene("LONGEST", 102, 20_000, 30_000));
    lines.extend(transcript(
        "NM_000003.1",
        "LONGEST",
        102,
        20_000,
        26_000,
        false,
    ));
    lines.extend(transcript(
        "NM_000004.1",
        "LONGEST",
        102,
        20_000,
        30_000,
        false,
    ));
    lines.extend(transcript(
        "NM_000005.1",
        "LONGEST",
        102,
        22_000,
        30_000,
        false,
    ));
    // two of one span: the first in the file
    lines.push(gene("TIE", 103, 40_000, 45_000));
    lines.extend(transcript("NM_000006.1", "TIE", 103, 40_000, 45_000, false));
    lines.extend(transcript("NM_000007.1", "TIE", 103, 40_000, 45_000, false));
    // a gene of two transcripts in windows of their own: the longer is the first, and a region over
    // one of them has not read the other
    lines.push(gene("SPLIT", 104, 50_000, 90_000));
    lines.extend(transcript(
        "NM_000008.1",
        "SPLIT",
        104,
        50_000,
        60_000,
        false,
    ));
    lines.extend(transcript(
        "NM_000009.1",
        "SPLIT",
        104,
        85_000,
        90_000,
        false,
    ));
    // the MANE Select transcript is the short one, in a window of its own, and a region over the
    // long one has not read it
    lines.push(gene("FARMANE", 105, 100_000, 130_000));
    lines.extend(transcript(
        "NM_000010.1",
        "FARMANE",
        105,
        100_000,
        105_000,
        true,
    ));
    lines.extend(transcript(
        "NM_000011.1",
        "FARMANE",
        105,
        115_000,
        130_000,
        false,
    ));
    sort_by_start(lines)
}

fn sort_by_start(mut lines: Vec<String>) -> Vec<String> {
    // `sort_by_key` is stable
    lines.sort_by_key(|l| l.split('\t').nth(3).unwrap().parse::<u64>().unwrap());
    lines
}

/// Write a bgzipped, tabix-indexed GFF3 of `lines`, one chunk to a line, so that a region reads
/// the lines of the bins it touches and no others.
fn write_tabix_gff3(dir: &Path, lines: &[String]) -> PathBuf {
    let gz_path = dir.join("refseq.gff3.gz");
    let mut writer = bgzf::io::Writer::new(File::create(&gz_path).unwrap());
    let mut indexer = Indexer::default();
    for line in lines {
        let start_vp = writer.virtual_position();
        writer.write_all(line.as_bytes()).unwrap();
        writer.write_all(b"\n").unwrap();
        let end_vp = writer.virtual_position();
        let f: Vec<&str> = line.split('\t').collect();
        let start = Position::try_from(f[3].parse::<usize>().unwrap()).unwrap();
        let end = Position::try_from(f[4].parse::<usize>().unwrap()).unwrap();
        indexer
            .add_record(Some((0, start, end, true)), Chunk::new(start_vp, end_vp))
            .unwrap();
    }
    writer.finish().unwrap();
    let header = Header::builder()
        .set_format(Format::Generic(CoordinateSystem::Gff))
        .set_reference_sequence_name_index(0)
        .set_start_position_index(3)
        .set_end_position_index(Some(4))
        .set_line_comment_prefix(b'#')
        .set_line_skip_count(0)
        .set_reference_sequence_names(std::iter::once("1").map(Into::into).collect())
        .build();
    let index = indexer.set_header(header).build(1);
    tabix::fs::write(dir.join("refseq.gff3.gz.tbi"), &index).unwrap();
    gz_path
}

/// transcript -> canonical, for what a read returned
fn flags(transcripts: &[fastvep_genome::Transcript]) -> BTreeMap<String, bool> {
    transcripts
        .iter()
        .map(|t| (t.stable_id.to_string(), t.canonical))
        .collect()
}

fn read(gz: &Path, start: u64, end: u64) -> BTreeMap<String, bool> {
    flags(&parse_gff3_indexed(gz, &[("1".to_string(), start, end)]).unwrap())
}

/// Every transcript a tabix read returns has the flag a read of the whole file gives it.
fn assert_same_as_whole_file(
    whole: &BTreeMap<String, bool>,
    tabix_read: &BTreeMap<String, bool>,
    what: &str,
) {
    assert!(!tabix_read.is_empty(), "{what}: nothing was read");
    for (id, canonical) in tabix_read {
        assert_eq!(
            Some(canonical),
            whole.get(id),
            "{what}: {id} is canonical={canonical} through tabix and canonical={:?} in a whole file; \
             the tabix read has {tabix_read:?}",
            whole.get(id)
        );
    }
}

#[test]
fn a_tabix_read_gives_each_transcript_the_canonical_flag_a_whole_file_read_gives_it() {
    let lines = world();
    let whole = flags(&parse_gff3(lines.join("\n").as_bytes()).unwrap());
    // what the rule says, so that "the same as the whole file" is not "the same wrong thing"
    let canonical: Vec<&str> = whole
        .iter()
        .filter(|(_, c)| **c)
        .map(|(id, _)| id.as_str())
        .collect();
    assert_eq!(
        canonical,
        [
            "NM_000001.1", // MANE Select
            "NM_000004.1", // the longest span
            "NM_000006.1", // the first of two of one span
            "NM_000008.1", // the longer of two
            "NM_000010.1", // MANE Select, short as it is
        ]
    );

    let dir = TempDir::new().unwrap();
    let gz = write_tabix_gff3(dir.path(), &lines);

    // every gene of the file is in the interval: nothing is partial
    assert_same_as_whole_file(&whole, &read(&gz, 1, 200_000), "the whole contig");

    // one gene, and the interval holds it
    let tie = read(&gz, 39_000, 46_000);
    assert_eq!(
        tie.keys().map(String::as_str).collect::<Vec<_>>(),
        ["NM_000006.1", "NM_000007.1"],
        "a region over the gene TIE reads its two transcripts"
    );
    assert_same_as_whole_file(&whole, &tie, "a region holding TIE whole");

    // a gene across the edge of the interval, every transcript of it overlapping the interval
    let longest = read(&gz, 21_000, 23_000);
    assert_eq!(longest.len(), 3, "{longest:?}");
    assert_same_as_whole_file(&whole, &longest, "a region inside LONGEST");
    assert!(longest["NM_000004.1"]);
}

/// The transcript the rule picks is read and the other one is not: nothing in what was read says
/// it is the longest, which only the rest of the gene does.
#[test]
fn a_gene_read_in_part_is_judged_by_all_of_its_transcripts() {
    let lines = world();
    let whole = flags(&parse_gff3(lines.join("\n").as_bytes()).unwrap());
    let dir = TempDir::new().unwrap();
    let gz = write_tabix_gff3(dir.path(), &lines);

    // SPLIT: NM_000008.1 (50,000-60,000) is read and NM_000009.1 (85,000-90,000) is not
    let split = read(&gz, 52_000, 53_000);
    assert!(split.contains_key("NM_000008.1"), "{split:?}");
    assert!(
        !split.contains_key("NM_000009.1"),
        "the layout no longer keeps the two apart: {split:?}"
    );
    assert!(
        split["NM_000008.1"],
        "it is the longest transcript of the gene, which what was read cannot say"
    );
    assert_same_as_whole_file(
        &whole,
        &split,
        "a region over the longer of two transcripts",
    );

    // ... and over the shorter one: it is not canonical, though it is the longest of what was read
    let short = read(&gz, 86_000, 87_000);
    assert!(
        short.contains_key("NM_000009.1") && !short.contains_key("NM_000008.1"),
        "{short:?}"
    );
    assert!(!short["NM_000009.1"]);
}

/// The MANE Select transcript is not among those that were read, and one that is longer than it
/// is: the longest of what was read is not the gene's canonical transcript.
#[test]
fn a_longer_transcript_is_not_canonical_because_the_mane_select_one_was_not_read() {
    let lines = world();
    let whole = flags(&parse_gff3(lines.join("\n").as_bytes()).unwrap());
    let dir = TempDir::new().unwrap();
    let gz = write_tabix_gff3(dir.path(), &lines);

    let far = read(&gz, 120_000, 125_000);
    assert!(far.contains_key("NM_000011.1"), "{far:?}");
    assert!(
        !far.contains_key("NM_000010.1"),
        "the layout no longer keeps the MANE Select transcript out of the read: {far:?}"
    );
    assert!(!far["NM_000011.1"]);
    assert_same_as_whole_file(
        &whole,
        &far,
        "a region over the longer, non-MANE transcript",
    );

    // and over the MANE Select one: it is
    let mane = read(&gz, 101_000, 102_000);
    assert!(mane["NM_000010.1"]);
    assert!(!mane.contains_key("NM_000011.1"), "{mane:?}");
}

/// Reading the rest of a gene is for the flag and for nothing else: its other transcripts are not
/// returned, and a record in them that a whole-file read stops on is not one the run has any use
/// for (the first read of the file did not meet it).
#[test]
fn the_rest_of_a_gene_is_read_for_its_flag_and_does_not_stop_the_run() {
    let mut lines = world();
    // a gene with a long transcript in one window and a short one in another, whose CDS has a
    // phase that is not one
    lines.push(gene("DAMAGED", 106, 140_000, 180_000));
    lines.extend(transcript(
        "NM_000012.1",
        "DAMAGED",
        106,
        140_000,
        147_000,
        false,
    ));
    lines.extend(transcript(
        "NM_000013.1",
        "DAMAGED",
        106,
        170_000,
        175_000,
        false,
    ));
    lines.push(
        "1\tBestRefSeq\tCDS\t170000\t175000\t.\t+\tx\tID=cds-NP_000013.1;Parent=rna-NM_000013.1;gbkey=CDS;protein_id=NP_000013.1"
            .to_string(),
    );
    let lines = sort_by_start(lines);
    assert!(
        parse_gff3(lines.join("\n").as_bytes()).is_err(),
        "a whole-file read stops on the phase"
    );
    let dir = TempDir::new().unwrap();
    let gz = write_tabix_gff3(dir.path(), &lines);

    let read = read(&gz, 141_000, 142_000);
    assert_eq!(
        read.keys().map(String::as_str).collect::<Vec<_>>(),
        ["NM_000012.1"],
        "only what overlaps the region is returned"
    );
    assert!(read["NM_000012.1"], "the longer of the gene's two");
}
