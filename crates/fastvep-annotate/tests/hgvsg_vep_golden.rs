//! `hgvsg_vep` against what Ensembl VEP 116.3 writes as `HGVSg`, for every allele of
//! every record of `fixtures/vep_hgvsg` (issue #139).
//!
//! The world is a genome of homopolymer, dinucleotide, trinucleotide and hexamer repeats,
//! with a run at each end of the contig and a mitochondrion, and variants chosen to make
//! the 3'-rule matter: every kind of change at the first, middle and last base of each run
//! and one either side, at both ends of the contig, and records of several alts.
//! `regen.sh` there runs `vep --hgvsg --allele_number` from the pinned 116.3 image over it
//! and keeps what it wrote for each (record, allele): that is `golden.tsv`. Nothing in it
//! was worked out here.
//!
//! What that covers, and a reading of VEP's source would not: where the 3'-rule stops at
//! the edge of the flank, an insertion that is a duplication on either side of itself,
//! a deletion that shifts all the way to the end of the contig, the `m.` of the
//! mitochondrion (VEP numbers a contig named `MT` that way, and `chrM` as `g.`), and a
//! reference allele in the VCF that is not the one in the genome. The records `e<N>` are at the
//! very edges of each contig: a change before the first base, one after the last, and a reference
//! allele that runs off the end, for which VEP writes nothing (its `hgvs_genomic` returns nothing
//! for a variation feature that is not inside the slice).

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::{anyhow, Result};
use fastvep_annotate::hgvsg_per_alt;
use fastvep_cache::providers::SequenceProvider;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_hgvsg")
        .join(name)
}

/// A FASTA read whole. Past the end of a contig it returns what there is, as Ensembl's
/// `subseq` does, rather than failing or wrapping round a circular contig.
struct Genome(HashMap<String, Vec<u8>>);

impl Genome {
    fn read(path: &Path) -> Self {
        let mut contigs: HashMap<String, Vec<u8>> = HashMap::new();
        let mut name = String::new();
        for line in std::fs::read_to_string(path).unwrap().lines() {
            if let Some(n) = line.strip_prefix('>') {
                name = n.to_string();
                contigs.entry(name.clone()).or_default();
            } else {
                contigs
                    .get_mut(&name)
                    .unwrap()
                    .extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
            }
        }
        Genome(contigs)
    }
}

impl SequenceProvider for Genome {
    fn fetch_sequence(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        let contig = self
            .0
            .get(chrom)
            .ok_or_else(|| anyhow!("no contig {chrom}"))?;
        if start == 0 || start as usize > contig.len() {
            return Err(anyhow!("{chrom}:{start} is outside the contig"));
        }
        let last = (end as usize).min(contig.len());
        Ok(contig[start as usize - 1..last].to_vec())
    }

    fn sequence_length(&self, chrom: &str) -> Option<u64> {
        self.0.get(chrom).map(|contig| contig.len() as u64)
    }
}

/// Alleles that Ensembl VEP 116.3 describes with a base the contig does not have, and that this
/// does not: `(record, ALLELE_NUM, what VEP writes)`. Every one is an allele of a record of
/// several alts at the very edge of a contig. VEP's check that a variant is inside the contig is
/// on the record as a whole (`$vf_start < 1 || $vf_end < 1 || $vf_start > $ref_length ||
/// $vf_end > $ref_length`), which such a record passes whatever its alleles are, so an insertion
/// before the first base comes out as `g.0dup` and one after the last as `g.1500_1501insA`. Neither
/// is a description of anything, and the same VEP writes nothing for the allele when it is the
/// only one. This writes nothing for them, as it does for the others.
const NAMES_A_BASE_THE_CONTIG_LACKS: [(&str, usize, &str); 6] = [
    ("e5", 1, "chr10:g.0dup"),
    ("e6", 1, "chr10:g.0dup"),
    ("e11", 1, "chr10:g.1500_1501insA"),
    ("e21", 1, "MT:m.0_1insA"),
    ("e22", 1, "MT:m.0_1insA"),
    ("e27", 1, "MT:m.1000_1001insA"),
];

#[test]
fn every_allele_gets_the_hgvsg_vep_wrote_for_it() {
    let genome = Genome::read(&fixture("genome.fa"));

    // (record, ALLELE_NUM) -> HGVSg
    let golden: HashMap<(String, usize), String> = std::fs::read_to_string(fixture("golden.tsv"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            ((c[0].to_string(), c[1].parse().unwrap()), c[2].to_string())
        })
        .collect();

    let mut compared = 0;
    let mut edges = 0;
    let mut wrong = Vec::new();
    for line in std::fs::read_to_string(fixture("in.vcf")).unwrap().lines() {
        if line.starts_with('#') {
            continue;
        }
        let vf = fastvep_io::vcf::parse_vcf_line(line).unwrap();
        let id = line.split('\t').nth(2).unwrap().to_string();
        let ours = hgvsg_per_alt(
            Some(&genome as &dyn SequenceProvider),
            &vf.position.chromosome,
            vf.position.start,
            vf.position.end,
            &vf.ref_allele,
            &vf.alt_alleles,
        );
        for (i, got) in ours.iter().enumerate() {
            let want = &golden[&(id.clone(), i + 1)];
            compared += 1;
            edges += usize::from(id.starts_with('e'));
            let got = got.as_deref().unwrap_or("-");
            if let Some((_, _, vep)) = NAMES_A_BASE_THE_CONTIG_LACKS
                .iter()
                .find(|(record, allele, _)| *record == id && *allele == i + 1)
            {
                // the list says what VEP wrote, so it cannot go stale unnoticed, and what is here
                if want != vep || got != "-" {
                    wrong.push(format!(
                        "{id} allele {} ({line}): VEP {want} (listed as {vep}), here {got}",
                        i + 1
                    ));
                }
                continue;
            }
            if got != want {
                wrong.push(format!(
                    "{id} allele {} ({line}): VEP {want}, here {got}",
                    i + 1
                ));
            }
        }
    }
    assert!(compared > 1400, "{compared} alleles compared");
    assert!(
        edges >= 38,
        "{edges} alleles at the edges of a contig compared"
    );
    assert!(
        wrong.is_empty(),
        "{} of {compared} alleles differ from VEP's, first ten:\n{}",
        wrong.len(),
        wrong
            .iter()
            .take(10)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}
