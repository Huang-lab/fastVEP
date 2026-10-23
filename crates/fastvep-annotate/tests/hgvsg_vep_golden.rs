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
//! reference allele in the VCF that is not the one in the genome.

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
}

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
            let got = got.as_deref().unwrap_or("-");
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
