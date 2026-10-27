//! The protein notation of a change in a repeat that goes on across an intron, against what
//! Ensembl VEP 116.3 writes for the same change (issue #139), and against itself.
//!
//! HGVS describes a change in the sequence of the transcript and puts it at its 3'-most position
//! there. Ensembl moves the change along the *genome* and describes where that ends. The two are
//! the same sequence until a repeat goes on over an intron, and then the genome carries the change
//! over bases the transcript does not have: it ends in the intron, or in the next exon as a
//! different change of the transcript.
//!
//! `fixtures/vep_tab_intron_repeat` has 140 transcripts with two coding exons and an intron of 1
//! to 14 bases, on both strands, in a repeat of 1 to 6 bases that runs through the end of the
//! first exon, the whole of the intron and the start of the second, and the same 140 with the
//! intron cut out of the genome: the transcript is the same sequence and so, to an annotator that
//! works in the transcript, the same change is the same description. 3,001 changes (insertions of
//! the repeat's unit, of a base and of a few; deletions of a base and of a unit) at every exonic
//! position of the repeat are made on both: `i_N` on the contig with the intron, `t_N` on its
//! twin. `regen.sh` there runs `vep` from the 116.3 image in `--gff` mode over them and keeps what
//! it wrote. Nothing in the golden was worked out here.
//!
//! What it shows, and what is asserted:
//!
//! - **fastVEP describes a change the same way whether or not the repeat goes over an intron**:
//!   consequence, residues, codons, positions and `HGVSp` are equal on every pair, 3,001 of them.
//! - **VEP does not**: it differs from its own twin on 565 of the 3,001 pairs, in `HGVSp` only:
//!   the change at another residue (346), `Ter?` where the stop is there to be counted (211), a
//!   different residue after the same position (5), and no notation at all where the twin has one
//!   (3).
//! - **On the contig with the intron fastVEP is VEP's in every cell but `HGVSp`**, and that
//!   differs only where VEP differs from its own twin, or where VEP is wrong on the twin as well.
//! - **On the twin the two agree** but for one thing, which is not about introns:
//!   an insertion on a reverse-strand transcript that the 3'-rule moves further than it is long is
//!   described with the wrong string (`shift_feature_seqs` turns it `seq_length - shift_length`
//!   times, a negative number, and so not at all). `vep_tab_edges_golden.rs` lists the same
//!   defect as `Why::Rotation`.
//!
//! `HGVSc` is not part of any of this: it is VEP's on every row of both contigs. It is a
//! description of the genome with the transcript's numbering, and it is the same genome.

use std::collections::{BTreeMap, BTreeSet};
use std::io::Read;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, PickFlags};
use flate2::read::GzDecoder;
use tempfile::TempDir;

fn here(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_tab_intron_repeat")
        .join(name)
}

type Row = BTreeMap<String, String>;
/// One row for each change: `--distance 0` and one transcript at each place.
type Table = BTreeMap<String, Row>;

fn table(text: &str) -> Table {
    let mut lines = text.lines().filter(|l| !l.starts_with("##"));
    let header: Vec<String> = lines
        .next()
        .expect("a header")
        .trim_start_matches('#')
        .split('\t')
        .map(String::from)
        .collect();
    let mut rows = Table::new();
    for line in lines {
        let cells: Vec<&str> = line.split('\t').collect();
        assert_eq!(cells.len(), header.len(), "ragged row {line:?}");
        let row: Row = header
            .iter()
            .cloned()
            .zip(cells.into_iter().map(String::from))
            .collect();
        let id = row["Uploaded_variation"].clone();
        assert!(rows.insert(id.clone(), row).is_none(), "two rows for {id}");
    }
    rows
}

fn golden() -> Table {
    let mut text = String::new();
    GzDecoder::new(std::fs::File::open(here("golden_flags.tab.gz")).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    table(&text)
}

fn ours() -> Table {
    let dir = TempDir::new().unwrap();
    // Copied, because the transcript cache a run writes goes beside its GFF3 and the index
    // of the FASTA beside it, and neither belongs in the source tree.
    let (gff3, fasta) = (dir.path().join("genes.gff3"), dir.path().join("genome.fa"));
    std::fs::copy(here("genes.gff3"), &gff3).unwrap();
    std::fs::copy(here("genome.fa"), &fasta).unwrap();
    let out = dir.path().join("out.tab");
    run_annotate(AnnotateConfig {
        input: here("in.vcf").to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        output_format: "tab".into(),
        pick: PickFlags::default(),
        hgvs: true,
        columns: ColumnFlags {
            show_ref_allele: true,
            variant_class: true,
            symbol: true,
            biotype: true,
            canonical: true,
            numbers: true,
            hgvsg: true,
            ..Default::default()
        },
        distance: 0,
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
        exclude: Default::default(),
    })
    .expect("annotation should succeed");
    table(&std::fs::read_to_string(&out).unwrap())
}

/// The `i_N` and `t_N` rows of a table, paired by `N`.
fn pairs(t: &Table) -> Vec<(&Row, &Row)> {
    t.iter()
        .filter_map(|(id, row)| {
            let n = id.strip_prefix("i_")?;
            Some((row, t.get(&format!("t_{n}")).expect("a twin")))
        })
        .collect()
}

/// A description without the name of the transcript or protein it is of, which differs between
/// the two contigs.
fn bare(cell: &str) -> &str {
    cell.split_once(':').map_or(cell, |(_, d)| d)
}

/// The terms of a consequence that have nothing to say about the repeat: where the intron is, the
/// contig without one has no splice site to be near.
fn without_splice_terms(consequence: &str) -> Vec<&str> {
    consequence
        .split(',')
        .filter(|t| {
            !matches!(
                *t,
                "splice_region_variant"
                    | "splice_acceptor_variant"
                    | "splice_donor_variant"
                    | "splice_donor_5th_base_variant"
                    | "splice_donor_region_variant"
                    | "splice_polypyrimidine_tract_variant"
                    | "intron_variant"
            )
        })
        .collect()
}

/// Columns that say where the change is in the transcript and what it does to the protein: the
/// same on both contigs for a change of one transcript.
const TRANSCRIPT_COLUMNS: &[&str] = &[
    "Amino_acids",
    "Codons",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
];

#[test]
fn both_tools_write_a_row_for_every_change_on_both_contigs() {
    let (vep, fv) = (golden(), ours());
    assert_eq!(vep.len(), 6002);
    assert!(
        vep.keys().eq(fv.keys()),
        "the changes the tools wrote a row for differ"
    );
    assert_eq!(pairs(&vep).len(), 3001);
}

/// fastVEP describes a change of a transcript the same way whether or not the repeat it is in goes
/// on over an intron.
#[test]
fn fastvep_describes_a_change_the_same_way_with_or_without_an_intron_in_the_repeat() {
    let fv = ours();
    let mut differ = Vec::new();
    for (i, t) in pairs(&fv) {
        let same = without_splice_terms(&i["Consequence"])
            == without_splice_terms(&t["Consequence"])
            && TRANSCRIPT_COLUMNS.iter().all(|c| i[*c] == t[*c])
            && bare(&i["HGVSp"]) == bare(&t["HGVSp"]);
        if !same {
            differ.push(format!(
                "{}: {} {} | {} {}",
                i["Uploaded_variation"], i["HGVSp"], i["Consequence"], t["HGVSp"], t["Consequence"]
            ));
        }
    }
    assert!(
        differ.is_empty(),
        "{} of 3001 pairs differ, e.g.\n{}",
        differ.len(),
        differ[..differ.len().min(10)].join("\n")
    );
}

/// What the golden shows of VEP alone, so that the comparisons below are about VEP and not about
/// the fixture.
#[test]
fn vep_describes_the_same_change_differently_when_the_repeat_goes_over_an_intron() {
    let vep = golden();
    let mut differ = 0usize;
    let (mut unwritten, mut unknown_stop, mut other) = (0usize, 0usize, 0usize);
    for (i, t) in pairs(&vep) {
        let (a, b) = (bare(&i["HGVSp"]), bare(&t["HGVSp"]));
        if a == b {
            continue;
        }
        differ += 1;
        if a == "-" {
            unwritten += 1;
        } else if a.ends_with("Ter?") && b.contains("fsTer") && !b.ends_with("Ter?") {
            unknown_stop += 1;
        } else {
            other += 1;
        }
        // the consequence, the residues and the positions are the same sequence and so the same
        assert_eq!(
            without_splice_terms(&i["Consequence"]),
            without_splice_terms(&t["Consequence"]),
            "{}",
            i["Uploaded_variation"]
        );
        for c in TRANSCRIPT_COLUMNS {
            assert_eq!(i[*c], t[*c], "{} {c}", i["Uploaded_variation"]);
        }
    }
    assert!(differ > 500, "{differ} of 3001 pairs");
    assert!(unwritten > 0, "{unwritten} changes with no HGVSp in VEP");
    assert!(
        unknown_stop > 150,
        "{unknown_stop} with `Ter?` for a stop that is there"
    );
    assert!(other > 250, "{other} with another change in the protein");
}

/// An insertion on a reverse-strand transcript that the 3'-rule moved further than it is long: the
/// string it inserts is not the one VEP describes (`Why::Rotation` of `vep_tab_edges_golden.rs`).
fn is_the_rotation_defect(vep: &Row) -> bool {
    let shifted = vep["HGVS_OFFSET"].trim_start_matches('-').parse::<usize>();
    vep["STRAND"] == "-1"
        && vep["VARIANT_CLASS"] == "insertion"
        && shifted.is_ok_and(|n| n > vep["Allele"].len())
}

/// Without an intron in the repeat the tools agree on every column, and on `HGVSp` but for the one
/// defect of VEP's that is not about introns.
#[test]
fn on_the_contig_without_the_intron_fastvep_is_vep_s_but_for_the_rotation_defect() {
    let (vep, fv) = (golden(), ours());
    let mut rotation = 0usize;
    for (id, v) in vep.iter().filter(|(id, _)| id.starts_with("t_")) {
        let f = &fv[id];
        for c in [
            "Location",
            "Consequence",
            "IMPACT",
            "Amino_acids",
            "Codons",
            "cDNA_position",
            "CDS_position",
            "Protein_position",
            "HGVS_OFFSET",
        ] {
            assert_eq!(v[c], f[c], "{id} {c}");
        }
        assert_eq!(bare(&v["HGVSc"]), bare(&f["HGVSc"]), "{id} HGVSc");
        if bare(&v["HGVSp"]) != bare(&f["HGVSp"]) {
            assert!(
                is_the_rotation_defect(v),
                "{id}: VEP {} fastVEP {} for {} {} (offset {})",
                v["HGVSp"],
                f["HGVSp"],
                v["STRAND"],
                v["VARIANT_CLASS"],
                v["HGVS_OFFSET"]
            );
            rotation += 1;
        }
    }
    assert!(rotation > 50, "{rotation} rows with the rotation defect");
}

/// With the intron in the repeat fastVEP is VEP's in every cell but `HGVSp`, and that differs only
/// where VEP is not the same as it is without the intron (the walk over the genome) or is wrong
/// there too.
#[test]
fn with_the_intron_fastvep_is_vep_s_in_every_cell_but_the_hgvsp_the_genome_walk_changes() {
    let (vep, fv) = (golden(), ours());
    let (mut walk, mut rotation) = (BTreeSet::new(), BTreeSet::new());
    for (id, v) in vep.iter().filter(|(id, _)| id.starts_with("i_")) {
        let f = &fv[id];
        for c in [
            "Location",
            "Consequence",
            "IMPACT",
            "Amino_acids",
            "Codons",
            "cDNA_position",
            "CDS_position",
            "Protein_position",
            "HGVS_OFFSET",
        ] {
            assert_eq!(v[c], f[c], "{id} {c}");
        }
        assert_eq!(bare(&v["HGVSc"]), bare(&f["HGVSc"]), "{id} HGVSc");
        if bare(&v["HGVSp"]) == bare(&f["HGVSp"]) {
            continue;
        }
        let twin = format!("t_{}", &id[2..]);
        if bare(&v["HGVSp"]) != bare(&vep[&twin]["HGVSp"]) {
            // VEP's own description changes with the intron, and fastVEP's does not
            walk.insert(id.clone());
        } else {
            assert!(
                is_the_rotation_defect(&vep[&twin]),
                "{id}: VEP {} fastVEP {}",
                v["HGVSp"],
                f["HGVSp"]
            );
            rotation.insert(id.clone());
        }
    }
    assert!(
        walk.len() > 400,
        "{} cells from the genome walk",
        walk.len()
    );
    assert!(rotation.len() > 20, "{} from the rotation", rotation.len());
}
