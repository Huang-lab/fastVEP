//! `--pick`: choosing one consequence per variant, in VEP's order.
//!
//! This lives in `fastvep-annotate` rather than in the CLI because both
//! drivers need it and they had already drifted: the CLI ran this hierarchy
//! while `annotate_vcf_text` (the web server's path) kept "canonical, or the
//! first transcript seen", which returns two rows for a single-gene variant
//! and puts a non-canonical one first. One implementation, one behaviour.

use anyhow::Result;
use fastvep_io::variant::TranscriptVariation;

/// One tier of the `--pick-order` hierarchy, in VEP's vocabulary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PickCriterion {
    ManeSelect,
    ManePlusClinical,
    Canonical,
    Appris,
    Tsl,
    Biotype,
    Ccds,
    Rank,
}

/// VEP's default `--pick_order`, which fastVEP matches exactly so that a
/// default run of each tool picks the same transcript.
///
/// Note where `Rank` sits: **last**. Transcript status outranks consequence
/// severity, so at a locus where a MANE transcript of one gene merely
/// neighbours the variant while a non-MANE transcript of another is disrupted
/// by it, both tools report the neighbour. That is correct VEP behaviour and
/// wrong clinical reporting - see `--pick-order` in docs/ACMG.md.
pub const DEFAULT_PICK_ORDER: &[PickCriterion] = &[
    PickCriterion::ManeSelect,
    PickCriterion::ManePlusClinical,
    PickCriterion::Canonical,
    PickCriterion::Appris,
    PickCriterion::Tsl,
    PickCriterion::Biotype,
    PickCriterion::Ccds,
    PickCriterion::Rank,
];

/// Parse a VEP-style `--pick_order` string, e.g. `rank,mane_select,canonical`.
pub fn parse_pick_order(spec: &str) -> Result<Vec<PickCriterion>> {
    let mut out = Vec::new();
    for raw in spec.split(',') {
        let name = raw.trim().to_ascii_lowercase();
        if name.is_empty() {
            continue;
        }
        let c = match name.as_str() {
            "mane_select" | "mane" => PickCriterion::ManeSelect,
            "mane_plus_clinical" => PickCriterion::ManePlusClinical,
            "canonical" => PickCriterion::Canonical,
            "appris" => PickCriterion::Appris,
            "tsl" => PickCriterion::Tsl,
            "biotype" => PickCriterion::Biotype,
            "ccds" => PickCriterion::Ccds,
            "rank" => PickCriterion::Rank,
            "length" => {
                // VEP's final tie-break. fastVEP's TranscriptVariation does not
                // carry transcript length, so honouring it would mean silently
                // doing nothing - worse than saying so.
                return Err(anyhow::anyhow!(
                    "--pick-order: 'length' is not supported; transcript length is not carried on \
                     the annotation record. Ties beyond the criteria you list are broken by \
                     transcript ID, which is deterministic."
                ));
            }
            other => {
                return Err(anyhow::anyhow!(
                    "--pick-order: unknown criterion {:?}. Valid: mane_select, mane_plus_clinical, \
                     canonical, appris, tsl, biotype, ccds, rank",
                    other
                ))
            }
        };
        if out.contains(&c) {
            return Err(anyhow::anyhow!(
                "--pick-order: {:?} listed more than once",
                name
            ));
        }
        out.push(c);
    }
    if out.is_empty() {
        return Err(anyhow::anyhow!("--pick-order: no criteria given"));
    }
    Ok(out)
}

/// Whether `--pick` has competing transcripts to choose between.
///
/// The placeholder rows `annotate_intergenic` and `annotate_sa_only_scaffold`
/// emit carry `transcript_id` `-`, and they are built one per *alt allele*
/// rather than one per transcript. Running the hierarchy over those keeps a
/// single row and silently drops every other alt from the output - no error,
/// just a missing allele. There is no transcript to pick at such a site, so
/// pick does not apply there.
///
/// A mix cannot currently occur, because those scaffolds run only when nothing
/// overlaps the variant; `any` is the deliberate reading if one ever does, so
/// that a real transcript still gets picked.
pub fn has_transcripts_to_pick(tvs: &[TranscriptVariation]) -> bool {
    tvs.len() > 1 && tvs.iter().any(|tv| tv.transcript_id.as_ref() != "-")
}

/// Index of the best transcript variation under the given `--pick-order`
/// hierarchy, with transcript_id alphabetical order as a final deterministic
/// tie-breaker.
///
/// Criteria omitted from `order` are not consulted at all, which is what lets
/// a caller drop a tier rather than only reorder it.
pub fn pick_best_transcript_idx_with(
    tvs: &[TranscriptVariation],
    order: &[PickCriterion],
) -> Option<usize> {
    (0..tvs.len())
        .min_by(|&a, &b| pick_key_with(&tvs[a], order).cmp(&pick_key_with(&tvs[b], order)))
}

/// Index of the best transcript variation under VEP's default `--pick_order`.
///
/// The production path resolves the order from `--pick-order` and calls
/// [`pick_best_transcript_idx_with`] directly; this is the tests' shorthand for
/// "what a default run does", which is the property most of them assert.
#[cfg(test)]
fn pick_best_transcript_idx(tvs: &[TranscriptVariation]) -> Option<usize> {
    pick_best_transcript_idx_with(tvs, DEFAULT_PICK_ORDER)
}

/// Score one transcript on one criterion. Lower is better, uniformly, so the
/// tiers compose by plain lexicographic comparison however they are ordered.
fn pick_score(tv: &TranscriptVariation, c: PickCriterion) -> u32 {
    match c {
        PickCriterion::ManeSelect => tv.mane_select.is_none() as u32,
        PickCriterion::ManePlusClinical => tv.mane_plus_clinical.is_none() as u32,
        PickCriterion::Canonical => !tv.canonical as u32,
        PickCriterion::Appris => appris_rank(tv.appris.as_deref()) as u32,
        PickCriterion::Tsl => tv.tsl.unwrap_or(u8::MAX) as u32,
        PickCriterion::Biotype => u32::from(tv.biotype.as_ref() != "protein_coding"),
        PickCriterion::Ccds => tv.ccds.is_none() as u32,
        PickCriterion::Rank => tv
            .allele_annotations
            .iter()
            .flat_map(|aa| aa.consequences.iter())
            .map(|c| c.rank())
            .min()
            .unwrap_or(u32::MAX),
    }
}

/// Upper bound on `--pick-order` length: there are eight criteria and
/// `parse_pick_order` rejects repeats, so no order can be longer.
const MAX_PICK_CRITERIA: usize = 8;

/// Score a transcript across the configured order.
///
/// Returns a fixed-size array rather than a `Vec` because this sits inside
/// `min_by`, which evaluates the key twice per comparison: a `Vec` here is two
/// heap allocations for every pair of transcripts considered, on every variant.
/// Unused slots stay zero and compare equal, which is harmless since every key
/// in a given comparison is built from the same `order`.
fn pick_key_with<'a>(
    tv: &'a TranscriptVariation,
    order: &[PickCriterion],
) -> ([u32; MAX_PICK_CRITERIA], &'a str) {
    let mut key = [0u32; MAX_PICK_CRITERIA];
    for (slot, &c) in key.iter_mut().zip(order.iter()) {
        *slot = pick_score(tv, c);
    }
    (key, tv.transcript_id.as_ref())
}

/// Map an APPRIS tag (`P1`/`principal1`, ..., `A1`/`alternative1`, ...) to a
/// rank where lower is better, matching VEP's `--pick_order` APPRIS tier:
/// principal1 < principal2 < ... < alternative1 < alternative2 < absent.
fn appris_rank(appris: Option<&str>) -> u8 {
    let Some(s) = appris else { return u8::MAX };
    let lower = s.to_ascii_lowercase();
    let (is_alt, digits) = if let Some(d) = lower.strip_prefix("principal") {
        (false, d)
    } else if let Some(d) = lower.strip_prefix("alternative") {
        (true, d)
    } else if let Some(d) = lower.strip_prefix('p') {
        (false, d)
    } else if let Some(d) = lower.strip_prefix('a') {
        (true, d)
    } else {
        // Present but unrecognised: still better than absent, worse than any
        // recognised principal/alternative.
        return u8::MAX - 1;
    };
    let n: u8 = digits.parse().unwrap_or(9);
    if is_alt {
        5u8.saturating_add(n)
    } else {
        n
    }
}

#[cfg(test)]
mod pick_tests {
    use super::*;
    use fastvep_core::{Allele, Consequence, Impact, Strand};
    use fastvep_io::variant::AlleleAnnotation;
    use std::sync::Arc;

    // Each argument is an independent coordinate or flag with no natural
    // grouping; bundling them would only move the list to the call site.
    #[allow(clippy::too_many_arguments)]
    fn make_tv(
        transcript_id: &str,
        canonical: bool,
        biotype: &str,
        consequences: Vec<Consequence>,
        mane_select: Option<&str>,
        mane_plus_clinical: Option<&str>,
        appris: Option<&str>,
        tsl: Option<u8>,
        ccds: Option<&str>,
    ) -> TranscriptVariation {
        TranscriptVariation {
            transcript_id: Arc::from(transcript_id),
            gene_id: Arc::from("GENE"),
            gene_symbol: Some(Arc::from("GENE")),
            biotype: Arc::from(biotype),
            allele_annotations: vec![AlleleAnnotation {
                allele: Allele::from_str("A"),
                consequences,
                impact: Impact::Modifier,
                cdna_position: None,
                cds_position: None,
                protein_position: None,
                amino_acids: None,
                codons: None,
                exon: None,
                intron: None,
                distance: None,
                protein_length: None,
                escapes_nmd: None,
                hgvsc: None,
                hgvsp: None,
                hgvsg: None,
                hgvs_offset: None,
                existing_variation: Vec::new(),
                sift: None,
                polyphen: None,
                supplementary: Vec::new(),
                acmg_classification: None,
            }],
            canonical,
            strand: Strand::Forward,
            source: None,
            protein_id: None,
            mane_select: mane_select.map(String::from),
            mane_plus_clinical: mane_plus_clinical.map(String::from),
            tsl,
            appris: appris.map(String::from),
            ccds: ccds.map(String::from),
            gencode_primary: false,
            symbol_source: None,
            hgnc_id: None,
            flags: Vec::new(),
        }
    }

    #[test]
    fn intergenic_placeholders_are_not_pickable() {
        // Regression test: a multi-allelic intergenic site arrives as one
        // placeholder row per alt allele. Picking among them drops alleles.
        let rows: Vec<TranscriptVariation> = ["-", "-"]
            .iter()
            .map(|id| {
                make_tv(
                    id,
                    false,
                    "-",
                    vec![Consequence::IntergenicVariant],
                    None,
                    None,
                    None,
                    None,
                    None,
                )
            })
            .collect();
        assert!(!has_transcripts_to_pick(&rows));
    }

    #[test]
    fn real_transcripts_are_pickable() {
        let rows = vec![
            make_tv(
                "ENST1",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                None,
                None,
                None,
            ),
            make_tv(
                "ENST2",
                true,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                None,
                None,
                None,
            ),
        ];
        assert!(has_transcripts_to_pick(&rows));
        // One transcript is not a choice either.
        assert!(!has_transcripts_to_pick(&rows[..1]));
    }

    #[test]
    fn pick_prefers_mane_select_over_canonical() {
        let tvs = vec![
            make_tv(
                "TX_CANON",
                true,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                None,
                None,
                None,
            ),
            make_tv(
                "TX_MANE",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                Some("TX_MANE.1"),
                None,
                None,
                None,
                None,
            ),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_prefers_mane_plus_clinical_over_canonical() {
        let tvs = vec![
            make_tv(
                "TX_CANON",
                true,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                None,
                None,
                None,
            ),
            make_tv(
                "TX_MANE_PC",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                Some("TX_MANE_PC.1"),
                None,
                None,
                None,
            ),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_prefers_mane_select_over_mane_plus_clinical() {
        let tvs = vec![
            make_tv(
                "TX_MANE_PC",
                true,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                Some("TX_MANE_PC.1"),
                None,
                None,
                None,
            ),
            make_tv(
                "TX_MANE",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                Some("TX_MANE.1"),
                None,
                None,
                None,
                None,
            ),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_falls_back_to_canonical_when_no_mane() {
        let tvs = vec![
            make_tv(
                "TX_NONCAN",
                false,
                "protein_coding",
                vec![Consequence::StopGained],
                None,
                None,
                None,
                None,
                None,
            ),
            make_tv(
                "TX_CANON",
                true,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                None,
                None,
                None,
            ),
        ];
        // Canonical wins even though TX_NONCAN has a more severe consequence.
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_prefers_protein_coding_biotype() {
        let tvs = vec![
            make_tv(
                "TX_NONCODING",
                false,
                "lncRNA",
                vec![Consequence::MissenseVariant],
                None,
                None,
                None,
                None,
                None,
            ),
            make_tv(
                "TX_PC",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                None,
                None,
                None,
            ),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_uses_severity_when_other_fields_equal() {
        let tvs = vec![
            make_tv(
                "TX_A",
                false,
                "protein_coding",
                vec![Consequence::SynonymousVariant],
                None,
                None,
                None,
                None,
                None,
            ),
            make_tv(
                "TX_B",
                false,
                "protein_coding",
                vec![Consequence::StopGained],
                None,
                None,
                None,
                None,
                None,
            ),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_tie_breaks_alphabetically_on_transcript_id() {
        let tvs = vec![
            make_tv(
                "TX_Z",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                None,
                None,
                None,
            ),
            make_tv(
                "TX_A",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                None,
                None,
                None,
            ),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_prefers_lower_tsl() {
        let tvs = vec![
            make_tv(
                "TX_TSL5",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                None,
                Some(5),
                None,
            ),
            make_tv(
                "TX_TSL1",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                None,
                Some(1),
                None,
            ),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_prefers_lower_appris_principal() {
        // P1 should beat P3 even though both are APPRIS-tagged — would fail
        // if APPRIS were compared by presence-only.
        let tvs = vec![
            make_tv(
                "TX_P3",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                Some("P3"),
                None,
                None,
            ),
            make_tv(
                "TX_P1",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                Some("P1"),
                None,
                None,
            ),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_prefers_principal_over_alternative_appris() {
        let tvs = vec![
            make_tv(
                "TX_A1",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                Some("A1"),
                None,
                None,
            ),
            make_tv(
                "TX_P5",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                Some("P5"),
                None,
                None,
            ),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_accepts_long_form_appris_tags() {
        // Ensembl GFF3 sometimes uses "principal1" / "alternative2".
        let tvs = vec![
            make_tv(
                "TX_ALT",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                Some("alternative2"),
                None,
                None,
            ),
            make_tv(
                "TX_PRINC",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                Some("principal1"),
                None,
                None,
            ),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_returns_none_for_empty_input() {
        let tvs: Vec<TranscriptVariation> = vec![];
        assert_eq!(pick_best_transcript_idx(&tvs), None);
    }

    // ── C5: configurable --pick-order ────────────────────────────────────

    #[test]
    fn pick_order_default_matches_vep_and_prefers_status_over_severity() {
        // The behaviour the round-2 review flagged, pinned as a fact rather
        // than left implicit: under VEP's default order a MANE transcript the
        // variant merely neighbours outranks a non-MANE one it disrupts.
        // CYP21A2 variants coming out on C4B `downstream_gene_variant`, and
        // STRC-region ones on TIMM9 `upstream_gene_variant`, are this rule.
        let tvs = vec![
            make_tv(
                "TX_DISRUPTED",
                false,
                "protein_coding",
                vec![Consequence::StartLost],
                None,
                None,
                None,
                None,
                None,
            ),
            make_tv(
                "TX_MANE_NEIGHBOUR",
                false,
                "protein_coding",
                vec![Consequence::UpstreamGeneVariant],
                Some("TX_MANE_NEIGHBOUR.1"),
                None,
                None,
                None,
                None,
            ),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_order_with_rank_first_reports_the_disrupted_transcript() {
        // The clinical order. This is what makes the CYP21A2 and KIAA0586
        // rows report the gene the variant actually hits.
        let tvs = vec![
            make_tv(
                "TX_DISRUPTED",
                false,
                "protein_coding",
                vec![Consequence::StartLost],
                None,
                None,
                None,
                None,
                None,
            ),
            make_tv(
                "TX_MANE_NEIGHBOUR",
                false,
                "protein_coding",
                vec![Consequence::UpstreamGeneVariant],
                Some("TX_MANE_NEIGHBOUR.1"),
                None,
                None,
                None,
                None,
            ),
        ];
        let order = parse_pick_order("rank,mane_select,canonical").unwrap();
        assert_eq!(pick_best_transcript_idx_with(&tvs, &order), Some(0));
    }

    #[test]
    fn pick_order_still_breaks_ties_by_the_later_criteria() {
        // Equal severity must fall through to MANE, or putting rank first
        // would turn every same-consequence choice into a transcript-ID sort.
        let tvs = vec![
            make_tv(
                "TX_PLAIN",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                None,
                None,
                None,
            ),
            make_tv(
                "TX_MANE",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                Some("TX_MANE.1"),
                None,
                None,
                None,
                None,
            ),
        ];
        let order = parse_pick_order("rank,mane_select,canonical").unwrap();
        assert_eq!(pick_best_transcript_idx_with(&tvs, &order), Some(1));
    }

    #[test]
    fn pick_order_parses_the_vep_default_spelling() {
        let order = parse_pick_order(
            "mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank",
        )
        .unwrap();
        assert_eq!(order, DEFAULT_PICK_ORDER.to_vec());
    }

    #[test]
    fn pick_order_omitted_criteria_are_not_consulted() {
        // Listing a subset drops the rest rather than appending them, which is
        // what lets a caller say "severity, then MANE, and nothing else".
        let tvs = vec![
            make_tv(
                "TX_A",
                true,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                None,
                None,
                None,
            ),
            make_tv(
                "TX_B",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                None,
                None,
                None,
            ),
        ];
        // Canonical would pick index 0; with only `rank` the tie falls to the
        // transcript-ID tie-break, which is TX_A anyway - so use ccds to show
        // an omitted criterion really is ignored.
        let order = parse_pick_order("rank").unwrap();
        assert_eq!(pick_best_transcript_idx_with(&tvs, &order), Some(0));
        let order = parse_pick_order("canonical,rank").unwrap();
        assert_eq!(pick_best_transcript_idx_with(&tvs, &order), Some(0));
    }

    #[test]
    fn pick_order_rejects_bad_input_rather_than_guessing() {
        for (spec, expect) in [
            ("rank,notacriterion", "unknown criterion"),
            ("rank,rank", "more than once"),
            ("", "no criteria"),
            ("rank,length", "not supported"),
        ] {
            let err = parse_pick_order(spec).expect_err("must reject").to_string();
            assert!(
                err.contains(expect),
                "{spec:?} gave {err:?}, wanted {expect:?}"
            );
        }
    }
}
