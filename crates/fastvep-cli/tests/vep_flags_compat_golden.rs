//! `--total_length` and `--no_escape` against what Ensembl VEP 116.3 writes for the same input
//! (issue #139).
//!
//! `fixtures/vep_flags_compat` holds what VEP wrote, with a `regen.sh` that runs it in Docker from the
//! image pinned there; nothing in the goldens was worked out here.
//!
//! - `golden_total_length.tsv.gz`: the three position cells of every row of `vep --tab --total_length`
//!   over the 2,271 changes at every exon and CDS edge of `fixtures/vep_tab_world` (4,753 rows). They
//!   come in every form `format_coords` has, with the transcript's length after each, and `-/1403`
//!   where a change touches an exon and has no coordinate in it: a deletion that covers an exon from
//!   intron to intron, an insertion before the initiator.
//! - `golden_no_escape.tsv.gz`: the `SYMBOL` and `HGVSp` of every entry of `vep --vcf` over a
//!   RefSeq-style world, as it is and with `--no_escape`: `p.Val3%3D` and `GENE%20A` are `p.Val3=` and
//!   `GENE_A`.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, CompatFlags, PickFlags};
use flate2::read::GzDecoder;
use tempfile::TempDir;

fn fixtures(dir: &str, name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures")
        .join(dir)
        .join(name)
}

fn golden(name: &str) -> Vec<Vec<String>> {
    let mut text = String::new();
    GzDecoder::new(std::fs::File::open(fixtures("vep_flags_compat", name)).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').map(String::from).collect())
        .collect()
}

/// A run of `fastvep annotate` over `input` against the genes and sequence of a fixture directory.
fn run(world: &str, input: &Path, format: &str, compat: CompatFlags, hgvs: bool) -> String {
    let dir = TempDir::new().unwrap();
    // Copied, because the transcript cache a run writes goes beside its GFF3, and neither it nor
    // the FASTA's index belongs in the source tree.
    let (gff3, fasta) = (dir.path().join("genes.gff3"), dir.path().join("genome.fa"));
    std::fs::copy(fixtures(world, "genes.gff3"), &gff3).unwrap();
    std::fs::copy(fixtures(world, "genome.fa"), &fasta).unwrap();
    let out = dir.path().join(format!("out.{format}"));
    run_annotate(AnnotateConfig {
        input: input.to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        output_format: format.into(),
        pick: PickFlags::default(),
        hgvs,
        columns: ColumnFlags {
            symbol: true,
            ..ColumnFlags::default()
        },
        compat,
        exclude: Default::default(),
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
    })
    .expect("annotation should succeed");
    std::fs::read_to_string(&out).unwrap()
}

fn header_index(header: &[&str], name: &str) -> usize {
    header
        .iter()
        .position(|h| *h == name)
        .unwrap_or_else(|| panic!("no column {name}"))
}

/// The CSQ entries of a VCF as maps by field name, with the name of the record that holds each:
/// the ID if it has one, otherwise the line VEP names it by, `{chr}_{pos}_{ref}/{alt}`.
fn csq_entries(vcf: &str) -> Vec<(String, BTreeMap<String, String>)> {
    let format: Vec<String> = vcf
        .lines()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .and_then(|l| l.split("Format: ").nth(1))
        .map(|f| {
            f.trim_end_matches("\">")
                .split('|')
                .map(String::from)
                .collect()
        })
        .expect("a CSQ header");
    let mut out = Vec::new();
    for line in vcf.lines().filter(|l| !l.starts_with('#')) {
        let c: Vec<&str> = line.split('\t').collect();
        let name = if c[2] == "." {
            format!("{}_{}_{}/{}", c[0], c[1], c[3], c[4].replace(',', "/"))
        } else {
            c[2].to_string()
        };
        let Some(csq) = c[7].split(';').find_map(|kv| kv.strip_prefix("CSQ=")) else {
            continue;
        };
        for entry in csq.split(',') {
            out.push((
                name.clone(),
                format
                    .iter()
                    .cloned()
                    .zip(entry.split('|').map(String::from))
                    .collect(),
            ));
        }
    }
    out
}

// --- --total_length ------------------------------------------------------------------------

/// The golden's rows by `(Uploaded_variation, Allele, Feature)`.
fn totals() -> BTreeMap<(String, String, String), [String; 3]> {
    golden("golden_total_length.tsv.gz")
        .into_iter()
        .map(|r| {
            (
                (r[0].clone(), r[1].clone(), r[2].clone()),
                [r[3].clone(), r[4].clone(), r[5].clone()],
            )
        })
        .collect()
}

fn edge_input() -> PathBuf {
    fixtures("vep_tab_edges", "in.vcf")
}

#[test]
fn the_tab_position_cells_with_total_length_are_vep_s_on_every_row() {
    let compat = CompatFlags {
        total_length: true,
        ..CompatFlags::default()
    };
    let text = run("vep_tab_world", &edge_input(), "tab", compat, false);
    let mut lines = text.lines().filter(|l| !l.starts_with("##"));
    let header: Vec<&str> = lines
        .next()
        .unwrap()
        .trim_start_matches('#')
        .split('\t')
        .collect();
    let at = |n: &str| header_index(&header, n);
    let (uploaded, allele, feature) = (at("Uploaded_variation"), at("Allele"), at("Feature"));
    let positions = [
        at("cDNA_position"),
        at("CDS_position"),
        at("Protein_position"),
    ];

    let vep = totals();
    let mut ours = BTreeMap::new();
    for line in lines {
        let c: Vec<&str> = line.split('\t').collect();
        ours.insert(
            (
                c[uploaded].to_string(),
                c[allele].to_string(),
                c[feature].to_string(),
            ),
            positions.map(|i| c[i].to_string()),
        );
    }
    assert_eq!(
        ours.keys().collect::<Vec<_>>(),
        vep.keys().collect::<Vec<_>>(),
        "the rows"
    );
    let wrong: Vec<String> = vep
        .iter()
        .filter(|(k, v)| &ours[*k] != *v)
        .take(10)
        .map(|(k, v)| {
            format!(
                "{} {}: VEP {v:?}, fastVEP {:?}",
                k.0.chars().take(40).collect::<String>(),
                k.2,
                ours[k]
            )
        })
        .collect();
    assert!(wrong.is_empty(), "{wrong:#?}");
}

#[test]
fn the_golden_has_every_form_of_a_cell_with_a_total() {
    // a cell with every run of digits made one `9`: `451-?/1403` is `9-?/9`
    let shape = |s: &str| -> String {
        let mut out = String::new();
        for c in s.chars() {
            if c.is_ascii_digit() {
                if !out.ends_with('9') {
                    out.push('9');
                }
            } else {
                out.push(c);
            }
        }
        out
    };
    let mut shapes = std::collections::BTreeSet::new();
    for cells in totals().values() {
        for cell in cells {
            shapes.insert(shape(cell));
        }
    }
    for wanted in ["9-9/9", "9-?/9", "?-9/9", "9/9", "-/9", "-"] {
        assert!(shapes.contains(wanted), "no `{wanted}` in {shapes:?}");
    }
}

#[test]
fn the_csq_position_fields_with_total_length_are_vep_s_with_the_dash_left_blank() {
    let compat = CompatFlags {
        total_length: true,
        ..CompatFlags::default()
    };
    let vcf = run("vep_tab_world", &edge_input(), "vcf", compat, false);
    let vep = totals();
    let mut seen = 0;
    for (record, entry) in csq_entries(&vcf) {
        // the CSQ leaves the `Feature` of an intergenic entry empty, where the tab says `-`
        let feature = match entry["Feature"].as_str() {
            "" => "-".to_string(),
            f => f.to_string(),
        };
        let key = (record.clone(), entry["Allele"].clone(), feature);
        let expected = &vep[&key];
        let ours = [
            entry["cDNA_position"].as_str(),
            entry["CDS_position"].as_str(),
            entry["Protein_position"].as_str(),
        ];
        for (theirs, ours) in expected.iter().zip(ours) {
            // VEP's VCF writer blanks a cell that is a bare `-` (`$data = '' if $data eq '-'`),
            // and `-/1403` is not one
            let theirs = if theirs == "-" { "" } else { theirs.as_str() };
            assert_eq!(theirs, ours, "{key:?}");
        }
        seen += 1;
    }
    assert_eq!(seen, vep.len());
}

// --- --no_escape ---------------------------------------------------------------------------

/// Entries where fastVEP and stock VEP 116.3 differ and VEP is the one that is wrong:
/// `(record, allele, transcript)` and the HGVSp fastVEP writes where VEP writes none. An
/// insertion at the end of an exon that the intron repeats is carried into the intron by the
/// 3'-rule's walk over the genome, and VEP then has no protein change for a frameshift; fastVEP
/// walks the sequence of the transcript (`docs/VEP_DIVERGENCE.md`).
const KNOWN: &[((&str, &str, &str), &str)] = &[(
    ("c323", "G", "NM_090002.2"),
    "NP_090002.2:p.Ile84ThrfsTer31",
)];

#[test]
fn the_symbol_and_hgvsp_of_every_csq_entry_are_vep_s_as_they_are_and_with_no_escape() {
    let input = fixtures("vep_flags_compat", "in.vcf");
    let golden = golden("golden_no_escape.tsv.gz");
    let mut checked = 0;
    for (mode, no_escape) in [("escaped", false), ("no_escape", true)] {
        let compat = CompatFlags {
            no_escape,
            ..CompatFlags::default()
        };
        let vcf = run("vep_flags_compat", &input, "vcf", compat, true);
        let ours: BTreeMap<(String, String, String), (String, String)> = csq_entries(&vcf)
            .into_iter()
            .map(|(record, e)| {
                (
                    (record, e["Allele"].clone(), e["Feature"].clone()),
                    (e["SYMBOL"].clone(), e["HGVSp"].clone()),
                )
            })
            .collect();
        let vep: BTreeMap<_, _> = golden
            .iter()
            .filter(|r| r[0] == mode)
            .map(|r| {
                (
                    (r[1].clone(), r[2].clone(), r[3].clone()),
                    (r[4].clone(), r[5].clone()),
                )
            })
            .collect();
        assert_eq!(ours.len(), vep.len(), "{mode}: the entries");
        let known = |k: &(String, String, String), v: &(String, String)| {
            KNOWN.iter().any(|((record, allele, transcript), hgvsp)| {
                (record, allele, transcript) == (&k.0.as_str(), &k.1.as_str(), &k.2.as_str())
                    && v.1.is_empty()
                    && ours.get(k).is_some_and(|o| o.0 == v.0 && o.1 == *hgvsp)
            })
        };
        let differing_known = vep
            .iter()
            .filter(|(k, v)| ours.get(*k) != Some(*v) && known(k, v))
            .count();
        assert_eq!(differing_known, KNOWN.len(), "{mode}: the known entries");
        let wrong: Vec<String> = vep
            .iter()
            .filter(|(k, v)| ours.get(*k) != Some(*v) && !known(k, v))
            .take(10)
            .map(|(k, v)| format!("{k:?}: VEP {v:?}, fastVEP {:?}", ours.get(k)))
            .collect();
        assert!(wrong.is_empty(), "{mode}: {wrong:#?}");
        checked += vep.len();
    }
    assert!(checked > 2000, "{checked} entries");
}

#[test]
fn the_golden_has_what_the_flag_is_about() {
    let golden = golden("golden_no_escape.tsv.gz");
    let cells = |mode: &str, column: usize| -> Vec<&str> {
        golden
            .iter()
            .filter(|r| r[0] == mode)
            .map(|r| r[column].as_str())
            .collect()
    };
    // an `=` escaped by default and not with the flag, and a space as `%20` and as `_`
    assert!(cells("escaped", 5).iter().any(|h| h.ends_with("%3D")));
    assert!(!cells("no_escape", 5).iter().any(|h| h.contains("%3D")));
    assert!(cells("no_escape", 5).iter().any(|h| h.ends_with('=')));
    assert!(cells("escaped", 4).contains(&"GENE%20A"));
    assert!(cells("no_escape", 4).contains(&"GENE_A"));
}
