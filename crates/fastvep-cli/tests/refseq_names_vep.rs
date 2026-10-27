//! A RefSeq GFF3 is named the way Ensembl VEP 116 names it: its transcripts, its genes and
//! their biotypes (issue #139, whose reporter runs `--refseq`).
//!
//! NCBI's GFF3 gives a transcript the ID `rna-NM_000001.3` and its `transcript_id` is
//! `NM_000001.3`; a gene the ID `gene-TP53` and a transcript the Dbxref `GeneID:7157`; and it
//! states no biotype, so there is nothing to read one from but the SO term and a few other
//! attributes. VEP's GFF reader (`BaseGXF.pm`, `GFF.pm`) names the transcript by its
//! `transcript_id`, the gene by the GeneID, and guesses the biotype: `protein_coding` for an
//! mRNA, the `ncrna_class` of an ncRNA, `miRNA` for a transcript of a gene described as a
//! microRNA, and the `gbkey` for the rest. fastVEP named them by ID, by the gene's own ID and by
//! SO term: on chr21 and chr22 of GRCh38.p14 and 527 variants of GIAB HG002, `rna-NM_...` in 2,550
//! rows where VEP has `NM_...`, `gene-TRG-GCC1-5` where it has `100189284`, and `mRNA` in 2,130 rows
//! where it has `protein_coding`.
//!
//! `fixtures/vep_refseq_gff` is a RefSeq-style GFF3 with one transcript of each kind a human
//! one has (an mRNA, a predicted mRNA, an lnc_RNA, a `transcript` record, a microRNA precursor, an
//! rRNA, an snRNA, an ncRNA with a class and a tRNA with no `transcript_id`), 18 variants, and what
//! Ensembl VEP 116.3 wrote for every CSQ entry (`regen.sh`, from the image pinned by digest). This
//! test runs fastVEP over the same files and holds the names and what depends on them to VEP's.
//!
//! The golden is VEP's `--gff` reader, and the reporter runs its *cache*, which names two kinds of
//! transcript differently: the biotype of an `lnc_RNA` of an lncRNA gene is `lncRNA` there and the
//! `gbkey`, `ncRNA`, to the GFF reader, and an `snRNA` of an snRNA gene is `snRNA` and `ncRNA`
//! (`CACHE_BIOTYPES`; measured on chr21 and chr22 of VEP 116.3's RefSeq cache, where all 503 and all
//! 13 + 12 + 1 of them read so). fastVEP follows the cache, so those two transcripts are listed one
//! by one, each asserted to still differ from the golden.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, PickFlags};
use tempfile::TempDir;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_refseq_gff")
        .join(name)
}

/// fastVEP over the fixture, in a directory of its own so that neither the transcript cache it
/// writes beside the GFF3 nor the FASTA index lands in the fixture.
fn annotate(dir: &Path) -> String {
    let (gff3, fasta, out) = (
        dir.join("genes.gff3"),
        dir.join("genome.fa"),
        dir.join("out.vcf"),
    );
    std::fs::copy(fixture("genes.gff3"), &gff3).unwrap();
    std::fs::copy(fixture("genome.fa"), &fasta).unwrap();
    run_annotate(AnnotateConfig {
        input: fixture("in.vcf").to_string_lossy().into(),
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
    .expect("annotation should succeed");
    std::fs::read_to_string(&out).unwrap()
}

type Key = (String, String, String);
type Cells = (String, String, String, String, String, String);

/// (record, allele, feature) -> (Gene, SYMBOL, BIOTYPE, Consequence, HGVSc, HGVSp)
fn entries(annotated: &str) -> BTreeMap<Key, Cells> {
    let at = |name: &str| {
        fastvep_io::output::DEFAULT_CSQ_FIELDS
            .iter()
            .position(|f| *f == name)
            .unwrap_or_else(|| panic!("{name} is not a CSQ field"))
    };
    let (allele, feature, gene, symbol, biotype, consequence, hgvsc, hgvsp) = (
        at("Allele"),
        at("Feature"),
        at("Gene"),
        at("SYMBOL"),
        at("BIOTYPE"),
        at("Consequence"),
        at("HGVSc"),
        at("HGVSp"),
    );
    let mut rows = BTreeMap::new();
    for line in annotated.lines().filter(|l| !l.starts_with('#')) {
        let cols: Vec<&str> = line.split('\t').collect();
        let csq = cols[7]
            .split(';')
            .find_map(|kv| kv.strip_prefix("CSQ="))
            .unwrap_or_else(|| panic!("no CSQ on:\n{line}"));
        for entry in csq.split(',') {
            let f: Vec<&str> = entry.split('|').collect();
            rows.insert(
                (
                    cols[2].to_string(),
                    f[allele].to_string(),
                    f[feature].to_string(),
                ),
                (
                    f[gene].to_string(),
                    f[symbol].to_string(),
                    f[biotype].to_string(),
                    f[consequence].to_string(),
                    f[hgvsc].to_string(),
                    f[hgvsp].to_string(),
                ),
            );
        }
    }
    rows
}

/// The transcripts whose biotype Ensembl VEP's RefSeq *cache* writes differently from its `--gff`
/// reader, with what the cache writes. Why: [`gene_biotype`](fastvep_cache) of the gene, as the
/// module comment says; `Transcript::biotype` has the rule (`refseq_biotype` in `gff.rs`).
const CACHE_BIOTYPES: [(&str, &str); 2] = [("NR_090003.1", "lncRNA"), ("XR_090007.1", "snRNA")];

fn golden() -> BTreeMap<Key, Cells> {
    std::fs::read_to_string(fixture("golden.tsv"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (
                (c[0].to_string(), c[1].to_string(), c[2].to_string()),
                (
                    c[3].to_string(),
                    c[4].to_string(),
                    c[5].to_string(),
                    c[6].to_string(),
                    c[7].to_string(),
                    c[8].to_string(),
                ),
            )
        })
        .collect()
}

#[test]
fn a_refseq_gff3_is_named_as_vep_names_it() {
    let dir = TempDir::new().unwrap();
    let ours = entries(&annotate(dir.path()));
    let want = golden();

    // the fixture has what it says it has
    let names: std::collections::BTreeSet<&str> = want.keys().map(|k| k.2.as_str()).collect();
    for name in [
        "NM_090001.3",
        "NR_090003.1",
        "XR_090004.1",
        "rna-TRN9",
        "NR_090005.1",
    ] {
        assert!(names.contains(name), "VEP names no {name}: {names:?}");
    }
    let biotypes: std::collections::BTreeSet<&str> = want.values().map(|c| c.2.as_str()).collect();
    for biotype in [
        "protein_coding",
        "ncRNA",
        "misc_RNA",
        "miRNA",
        "rRNA",
        "scaRNA",
        "tRNA",
    ] {
        assert!(
            biotypes.contains(biotype),
            "VEP writes no {biotype}: {biotypes:?}"
        );
    }

    // the cache's names for the two kinds of transcript it names differently
    let mut want = want;
    let on_gff_reader: Vec<Key> = want
        .keys()
        .filter(|k| CACHE_BIOTYPES.iter().any(|(t, _)| *t == k.2))
        .cloned()
        .collect();
    assert!(
        !on_gff_reader.is_empty(),
        "the fixture has the transcripts the cache names differently"
    );
    for key in on_gff_reader {
        let cells = want.get_mut(&key).unwrap();
        let (_, cache) = CACHE_BIOTYPES.iter().find(|(t, _)| *t == key.2).unwrap();
        assert_eq!(cells.2, "ncRNA", "VEP's `--gff` reader guesses the gbkey");
        cells.2 = cache.to_string();
    }

    let wrong: Vec<String> = want
        .iter()
        .filter(|(k, v)| ours.get(*k) != Some(*v))
        .map(|(k, v)| format!("{} {} {}: VEP {v:?}, here {:?}", k.0, k.1, k.2, ours.get(k)))
        .collect();
    assert!(
        wrong.is_empty(),
        "{} of {} entries differ from VEP 116.3's, first ten:\n{}",
        wrong.len(),
        want.len(),
        wrong
            .iter()
            .take(10)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
    let extra: Vec<&Key> = ours.keys().filter(|k| !want.contains_key(*k)).collect();
    assert!(extra.is_empty(), "entries VEP does not write: {extra:?}");
}
