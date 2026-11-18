//! A tabix (`.gz` + `.tbi`) GFF3 is read region-restricted: only features
//! overlapping the input VCF's variants come back. That set is correct for the
//! VCF that produced it and wrong for every other one, so it must never be
//! persisted as a `<gff3>.fastvep.cache` sidecar - the sidecar is keyed on the
//! GFF3 path alone and looks fresh to every later run.
//!
//! Before the gate in `run_annotate`, annotating one VCF and then a second one
//! against the same `--gff3` made the second run load the first run's
//! neighbourhood and report its variants as intergenic. On
//! validation/human/vep_example_GRCh38.vcf that was 171 of 173 variants wrong,
//! including a missense call reduced to `intergenic_variant`, with exit code 0
//! and no warning.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, PickFlags};
use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::binning_index::index::header::format::CoordinateSystem;
use noodles_csi::binning_index::index::header::Format;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_csi::binning_index::index::Header;
use noodles_csi::binning_index::Indexer;
use noodles_tabix as tabix;
use tempfile::TempDir;

/// Two genes far apart on the same contig, each a single-exon coding
/// transcript. A VCF hitting only one of them must not poison the other.
const GENE_A_CDS: u64 = 1_001; // ATG TGG CGG TAA at 1001..1012
const GENE_B_CDS: u64 = 900_001;

fn gff3_lines() -> Vec<String> {
    let mut v = Vec::new();
    for (tag, cds) in [("A", GENE_A_CDS), ("B", GENE_B_CDS)] {
        let end = cds + 11;
        v.push(format!(
            "1\ttest\tgene\t{cds}\t{end}\t.\t+\t.\tID=gene:ENSG_{tag};biotype=protein_coding"
        ));
        v.push(format!(
            "1\ttest\tmRNA\t{cds}\t{end}\t.\t+\t.\tID=transcript:ENST_{tag};Parent=gene:ENSG_{tag};biotype=protein_coding;tag=Ensembl_canonical"
        ));
        v.push(format!(
            "1\ttest\texon\t{cds}\t{end}\t.\t+\t.\tID=exon:ENSE_{tag};Parent=transcript:ENST_{tag};rank=1"
        ));
        v.push(format!(
            "1\ttest\tCDS\t{cds}\t{end}\t.\t+\t0\tID=CDS:ENSP_{tag};Parent=transcript:ENST_{tag}"
        ));
    }
    v
}

/// Write a bgzipped, tabix-indexed GFF3. Records are emitted in coordinate
/// order, which is what tabix requires and what `sort -k1,1 -k4,4n -k5,5n`
/// produces.
fn write_tabix_gff3(dir: &Path) -> PathBuf {
    let gz_path = dir.join("ann.gff3.gz");
    let mut writer = bgzf::io::Writer::new(File::create(&gz_path).unwrap());
    let mut indexer = Indexer::default();

    let mut lines = gff3_lines();
    lines.sort_by_key(|l| {
        let f: Vec<&str> = l.split('\t').collect();
        (
            f[0].to_string(),
            f[3].parse::<u64>().unwrap(),
            f[4].parse::<u64>().unwrap(),
        )
    });

    for line in &lines {
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
    tabix::fs::write(dir.join("ann.gff3.gz.tbi"), &index).unwrap();

    gz_path
}

fn write_fasta(dir: &Path) -> PathBuf {
    // Two coding stretches at the gene positions, N everywhere else.
    let mut seq = vec![b'N'; 900_020];
    for cds in [GENE_A_CDS, GENE_B_CDS] {
        let off = (cds - 1) as usize;
        seq[off..off + 12].copy_from_slice(b"ATGTGGCGGTAA");
    }
    let path = dir.join("ref.fa");
    let mut f = File::create(&path).unwrap();
    writeln!(f, ">1").unwrap();
    for chunk in seq.chunks(60) {
        f.write_all(chunk).unwrap();
        f.write_all(b"\n").unwrap();
    }
    path
}

fn write_vcf(dir: &Path, name: &str, pos: u64) -> PathBuf {
    let path = dir.join(name);
    let mut f = File::create(&path).unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    // Middle base of the Arg codon (CGG -> CAG): a missense call, distinct from
    // both the synonymous third base and the stop the Trp codon would give.
    writeln!(f, "1\t{pos}\t.\tG\tA\t.\t.\t.").unwrap();
    path
}

fn annotate(input: &Path, gff3: &Path, fasta: &Path, out: &Path) {
    run_annotate(AnnotateConfig {
        input: input.to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        output_format: "vcf".into(),
        pick: PickFlags::default(),
        hgvs: true,
        distance: 5000,
        cache_dir: None,
        transcript_cache: None,
        sa_dir: None,
        sa_only: false,
        acmg: false,
        acmg_config: None,
        pick_order: None,
        functional_evidence: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
        show_progress: false,
        compat: Default::default(),
        columns: Default::default(),
        exclude: Default::default(),
    })
    .unwrap();
}

fn csq_of(out: &Path) -> String {
    std::fs::read_to_string(out)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn a_region_restricted_tabix_load_is_never_persisted_as_a_sidecar_cache() {
    let dir = TempDir::new().unwrap();
    let gff3 = write_tabix_gff3(dir.path());
    let fasta = write_fasta(dir.path());
    let sidecar = dir.path().join("ann.gff3.gz.fastvep.cache");

    // Run 1 hits gene A only. Its region-restricted transcript set knows
    // nothing about gene B, 900 kb away.
    let vcf_a = write_vcf(dir.path(), "a.vcf", GENE_A_CDS + 7);
    let out_a = dir.path().join("a.out.vcf");
    annotate(&vcf_a, &gff3, &fasta, &out_a);
    assert!(
        csq_of(&out_a).contains("missense_variant"),
        "run 1 should annotate its own gene: {}",
        csq_of(&out_a)
    );
    assert!(
        !sidecar.exists(),
        "a region-restricted load must not be written to {}; \
         the sidecar is keyed on the GFF3 path and would look fresh to every later run",
        sidecar.display()
    );

    // Run 2 hits gene B, against the same --gff3. Before the gate it loaded
    // run 1's cache and called this intergenic.
    let vcf_b = write_vcf(dir.path(), "b.vcf", GENE_B_CDS + 7);
    let out_b = dir.path().join("b.out.vcf");
    annotate(&vcf_b, &gff3, &fasta, &out_b);
    let csq_b = csq_of(&out_b);
    assert!(
        csq_b.contains("missense_variant"),
        "run 2 must see its own gene, not run 1's neighbourhood: {csq_b}"
    );
    assert!(
        !csq_b.contains("intergenic_variant"),
        "run 2 fell back to a stale region-restricted cache: {csq_b}"
    );
    assert!(
        csq_b.contains("ENST_B"),
        "run 2 should hit gene B's transcript: {csq_b}"
    );
}

#[test]
fn a_whole_file_gff3_load_still_writes_its_sidecar_cache() {
    // The gate must not cost the optimisation it protects: a plain (non-tabix)
    // GFF3 is parsed in full, so its cache is valid for any input.
    let dir = TempDir::new().unwrap();
    let plain = dir.path().join("ann.gff3");
    std::fs::write(dir.path().join("ann.gff3"), gff3_lines().join("\n") + "\n").unwrap();
    let fasta = write_fasta(dir.path());
    let sidecar = dir.path().join("ann.gff3.fastvep.cache");

    let vcf = write_vcf(dir.path(), "a.vcf", GENE_A_CDS + 7);
    let out = dir.path().join("a.out.vcf");
    annotate(&vcf, &plain, &fasta, &out);

    assert!(
        sidecar.exists(),
        "a whole-file load should still be cached at {}",
        sidecar.display()
    );

    // And the cached set must serve a different region correctly.
    let vcf_b = write_vcf(dir.path(), "b.vcf", GENE_B_CDS + 7);
    let out_b = dir.path().join("b.out.vcf");
    annotate(&vcf_b, &plain, &fasta, &out_b);
    let csq_b = csq_of(&out_b);
    assert!(
        csq_b.contains("missense_variant") && csq_b.contains("ENST_B"),
        "cached whole-file load should annotate gene B: {csq_b}"
    );
}

/// A VCF of one record, written as it is.
fn write_vcf_record(dir: &Path, name: &str, record: &str) -> PathBuf {
    let path = dir.join(name);
    let mut f = File::create(&path).unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    writeln!(f, "{record}").unwrap();
    path
}

/// A structural variant is read from a tabix GFF3 for the whole of what it reaches, not for the
/// `POS +- distance` of a small variant: a deletion of 900 kb covers the transcript 900 kb away
/// from its position, and a breakend's mate is paired with the transcript around the mate. Before the
/// prescan read the records it loaded only the transcripts around `POS`, and the rows of the others
/// were missing from the answer with exit code 0 and nothing said.
#[test]
fn a_structural_variant_loads_the_transcripts_its_span_and_its_mates_reach() {
    let dir = TempDir::new().unwrap();
    let gff3 = write_tabix_gff3(dir.path());
    let fasta = write_fasta(dir.path());

    let deletion = write_vcf_record(
        dir.path(),
        "del.vcf",
        "1\t900\tdel\tN\t<DEL>\t.\t.\tEND=950000",
    );
    let out = dir.path().join("del.out.vcf");
    annotate(&deletion, &gff3, &fasta, &out);
    let csq = csq_of(&out);
    for transcript in ["ENST_A", "ENST_B"] {
        assert!(
            csq.contains(&format!("|{transcript}|")),
            "the deletion covers {transcript}, which is {} bases from its position: {csq}",
            if transcript == "ENST_A" { 100 } else { 900_000 }
        );
    }
    assert!(
        csq.contains("transcript_ablation"),
        "a deletion that covers a transcript is its ablation: {csq}"
    );

    let breakend = write_vcf_record(
        dir.path(),
        "bnd.vcf",
        &format!("1\t2000\tbnd\tN\tN[1:{}[\t.\t.\tSVTYPE=BND", GENE_B_CDS + 5),
    );
    let out = dir.path().join("bnd.out.vcf");
    annotate(&breakend, &gff3, &fasta, &out);
    let csq = csq_of(&out);
    assert!(
        csq.contains("|ENST_B|"),
        "the mate of the breakend is in ENST_B: {csq}"
    );
    assert!(
        csq.contains("feature_truncation"),
        "a breakend inside a transcript is its truncation: {csq}"
    );
}

/// An index addresses positions up to `(1 << (min_shift + 3 * depth)) - 1` (536,870,911 for a
/// `.tbi`) and refuses an interval that reaches past it, as a whole. One record whose reach widened
/// the region past it (a deletion with `END=2000000000`, which the run then skips as larger than
/// `--max_sv_size`) left the contig with no transcript at all, and every variant on it intergenic,
/// with exit code 0 and no warning: the interval is cut to what the index covers.
#[test]
fn a_record_that_reaches_past_what_the_index_addresses_does_not_void_its_contig() {
    let dir = TempDir::new().unwrap();
    let gff3 = write_tabix_gff3(dir.path());
    let fasta = write_fasta(dir.path());
    let vcf = write_vcf_record(
        dir.path(),
        "past.vcf",
        &format!(
            "1\t{}\t.\tG\tA\t.\t.\t.\n1\t1200\t.\tN\t<DEL>\t.\t.\tEND=2000000000",
            GENE_A_CDS + 7
        ),
    );
    let out = dir.path().join("past.out.vcf");
    annotate(&vcf, &gff3, &fasta, &out);
    let csq = csq_of(&out);
    assert!(
        csq.contains("|ENST_A|") && csq.contains("missense_variant"),
        "the SNV in gene A is annotated against it: {csq}"
    );
}
