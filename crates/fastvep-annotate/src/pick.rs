//! `--pick`: choosing one consequence per variant, in VEP's order.
//!
//! This lives in `fastvep-annotate` rather than in the CLI because both
//! drivers need it and they had already drifted: the CLI ran this hierarchy
//! while `annotate_vcf_text` (the web server's path) kept "canonical, or the
//! first transcript seen", which returns two rows for a single-gene variant
//! and puts a non-canonical one first. One implementation, one behaviour.
//!
//! The scoring is that one implementation: it is written against [`PickRow`], and everything that
//! picks goes through it. What a pick *keeps* is decided by [`pick_winners`], which the library
//! applies to finished rows ([`apply_pick`]) and the CLI, when it reduces, applies to what the
//! predictor said before any row is built (see [`PredictedRow`]). That is two ways of reaching one
//! answer, so `tests/loop_parity.rs` holds the CLI to the library's and `tests/pick_before_rows.rs`
//! holds a reducing run to the rows its flagging counterpart marks.

use std::collections::HashMap;

use anyhow::Result;
use fastvep_consequence::{AlleleConsequenceResult, Flank, TranscriptConsequence};
use fastvep_core::Allele;
use fastvep_genome::Transcript;
use fastvep_io::variant::TranscriptVariation;
use foldhash::fast::FixedState;
use serde::Deserialize;

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
/// "Matches exactly" is a claim about the tiers *and* about what feeds them.
/// It was false for two years on the APPRIS tier: `Transcript::appris` came
/// from the GFF3 parser, which never read the `appris_*` tags, so the tier
/// scored every candidate the same and the pick fell through to the next one
/// that could tell them apart. The three unit tests below covered the tier and
/// passed throughout, because they build `TranscriptVariation` by hand. A tier
/// added here needs a test that reaches it through the parser.
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

/// What one pick is chosen *per*, matching VEP's option family.
///
/// Measured against Ensembl VEP 115.1 on a two-alt site in TP53, which
/// produces 76 consequence entries unpicked: `--pick` leaves 1,
/// `--pick_allele` leaves 2, `--pick_allele_gene` leaves 2 (one gene at that
/// locus), and `--flag_pick_allele_gene` leaves all 76 with 2 flagged.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PickScope {
    /// One pick for the whole variant.
    ///
    /// **This is where fastVEP diverges from VEP, deliberately.** VEP's
    /// `--pick` reduces that TP53 site to a single entry, which means one of
    /// the two alt alleles is simply absent from the output - no error, no
    /// note. fastVEP keeps the winning transcript's annotation for *every*
    /// allele, so `--pick` there leaves 2 entries rather than 1. A dropped alt
    /// is indistinguishable from an alt that had no consequence, and this is a
    /// file a clinician reads. `--pick-allele` is the option for one entry per
    /// allele, and it agrees with VEP exactly. See docs/VEP_DIVERGENCE.md.
    #[default]
    Variant,
    /// One pick per (variant, allele).
    Allele,
    /// One pick per (variant, allele, gene).
    AlleleGene,
}

/// What happens to the (transcript, allele) pairs a pick did not choose.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PickMode {
    /// Delete them, leaving only the winners. VEP's `--pick*`.
    #[default]
    Reduce,
    /// Keep them, and set `pick` on the winners. VEP's `--flag_pick*`.
    ///
    /// Costs what it sounds like: every pass that runs after the pick - SA
    /// lookup, gene annotation, ACMG classification - runs on every
    /// transcript, because that is what "retains others" means. Reduce mode
    /// exists partly to avoid exactly that.
    Flag,
}

/// A resolved pick request: what to pick per, what to do with the losers, and
/// the hierarchy to decide with.
#[derive(Debug, Clone, Copy)]
pub struct PickPlan<'a> {
    pub scope: PickScope,
    pub mode: PickMode,
    pub order: &'a [PickCriterion],
}

impl<'a> PickPlan<'a> {
    /// A plan equivalent to the `--pick` fastVEP has always had.
    pub fn reduce_per_variant(order: &'a [PickCriterion]) -> Self {
        Self {
            scope: PickScope::Variant,
            mode: PickMode::Reduce,
            order,
        }
    }
}

/// Read the six switches out of a JSON request body.
///
/// Both HTTP entry points go through this rather than naming the keys
/// themselves, so the wire contract has one definition; the `Deserialize` on
/// [`PickFlags`] is what actually spells the names.
///
/// Takes the body by *reference*. `serde_json::from_value` consumes its
/// argument, so reading six bools out of a request meant cloning the whole
/// thing - and the body of the endpoint that carries these is a VCF, which is
/// the one field in it that can be megabytes.
///
/// Returns the error rather than defaulting: all six are read in one pass, so
/// a wrong type on any one of them would otherwise discard a correctly spelled
/// sibling and leave the caller annotating with no pick at all.
pub fn pick_flags_from_json(body: &serde_json::Value) -> Result<PickFlags, serde_json::Error> {
    PickFlags::deserialize(body)
}

/// A pick request, before a `--pick-order` has been resolved to attach to it.
///
/// Separate from [`PickPlan`] only because the plan borrows the order, and a
/// run's configuration is built before the order is parsed - a bad
/// `--pick-order` should fail the run at startup, not per variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PickRequest {
    pub scope: PickScope,
    pub mode: PickMode,
}

impl PickRequest {
    /// Attach a resolved hierarchy, giving something [`apply_pick`] can run.
    pub fn plan<'a>(&self, order: &'a [PickCriterion]) -> PickPlan<'a> {
        PickPlan {
            scope: self.scope,
            mode: self.mode,
            order,
        }
    }

    /// Whether the output needs a `PICK` column. Only flagging does: a
    /// reducing pick answers the question by what it leaves behind.
    pub fn needs_pick_column(&self) -> bool {
        self.mode == PickMode::Flag
    }
}

/// The six mutually exclusive pick switches as the command line spells them.
///
/// Named fields rather than six positional `bool`s, because the call site
/// reads `--flag-pick-allele` and `--flag-pick-allele-gene` next to each
/// other and a transposition between them would be invisible.
///
/// `Deserialize` so that the two HTTP entry points - fastvep-web's typed
/// handler and the CLI's bundled server - name these switches from one
/// definition rather than each spelling six JSON keys. They are two
/// implementations of one wire contract, and the annotation path behind them
/// has already drifted once (see the module header). Missing keys default to
/// `false`, so a body that mentions none of them asks for no pick.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Deserialize)]
#[serde(default)]
pub struct PickFlags {
    pub pick: bool,
    pub pick_allele: bool,
    pub pick_allele_gene: bool,
    pub flag_pick: bool,
    pub flag_pick_allele: bool,
    pub flag_pick_allele_gene: bool,
}

impl PickFlags {
    /// The request these switches describe, or `None` for a run that picks
    /// nothing.
    ///
    /// The CLI declares all six in one clap group, so at most one can be set
    /// from a command line. A caller that sets more than one by hand gets the
    /// first match in this order, which is deterministic rather than correct -
    /// there is no correct answer to two conflicting requests.
    pub fn resolve(self) -> Option<PickRequest> {
        let (scope, mode) = if self.pick {
            (PickScope::Variant, PickMode::Reduce)
        } else if self.pick_allele {
            (PickScope::Allele, PickMode::Reduce)
        } else if self.pick_allele_gene {
            (PickScope::AlleleGene, PickMode::Reduce)
        } else if self.flag_pick {
            (PickScope::Variant, PickMode::Flag)
        } else if self.flag_pick_allele {
            (PickScope::Allele, PickMode::Flag)
        } else if self.flag_pick_allele_gene {
            (PickScope::AlleleGene, PickMode::Flag)
        } else {
            return None;
        };
        Some(PickRequest { scope, mode })
    }

    /// The switch a resolved request came from, for a message that has to name
    /// the option the user actually typed.
    pub fn requested_option(self) -> Option<&'static str> {
        [
            (self.pick, "--pick"),
            (self.pick_allele, "--pick-allele"),
            (self.pick_allele_gene, "--pick-allele-gene"),
            (self.flag_pick, "--flag-pick"),
            (self.flag_pick_allele, "--flag-pick-allele"),
            (self.flag_pick_allele_gene, "--flag-pick-allele-gene"),
        ]
        .into_iter()
        .find_map(|(set, name)| set.then_some(name))
    }
}

/// What a pick reads of one transcript's row.
///
/// Two things are a row. [`TranscriptVariation`] is the finished one, which [`apply_pick`]
/// reduces or flags. [`PredictedRow`] is the predictor's answer for one transcript, borrowed and
/// not yet a row, so that a run can settle the pick *before* it builds any rows and build only
/// the ones it keeps. A variant of a run with `--distance 500000` has hundreds of transcripts
/// around it and the pick keeps a few tens; building the rest to delete them was most of what
/// such a run did. Both go through the one scoring below, so what is picked cannot depend on
/// which of the two it was asked of.
pub trait PickRow {
    fn transcript_id(&self) -> &str;
    fn gene_id(&self) -> &str;
    fn canonical(&self) -> bool;
    fn has_mane_select(&self) -> bool;
    fn has_mane_plus_clinical(&self) -> bool;
    fn appris(&self) -> Option<&str>;
    fn tsl(&self) -> Option<u8>;
    fn biotype(&self) -> &str;
    fn has_ccds(&self) -> bool;
    /// How many alleles the row has. They are numbered `0..allele_count()`, in the order the row
    /// writes them.
    fn allele_count(&self) -> usize;
    /// The allele numbered `index`, which has to be below [`allele_count`](Self::allele_count).
    fn allele(&self, index: usize) -> &Allele;
    /// The lowest consequence rank among the alleles equal to `allele`, or among all of them for
    /// `None`; `u32::MAX` when there is none.
    fn best_rank(&self, allele: Option<&Allele>) -> u32;
    /// A placeholder row, which stands for no transcript (see [`TranscriptVariation::NO_FEATURE`]):
    /// no pick chooses it and none deletes it.
    fn is_featureless(&self) -> bool {
        self.transcript_id() == TranscriptVariation::NO_FEATURE
    }
}

impl PickRow for TranscriptVariation {
    fn transcript_id(&self) -> &str {
        &self.transcript_id
    }
    fn gene_id(&self) -> &str {
        &self.gene_id
    }
    fn canonical(&self) -> bool {
        self.canonical
    }
    fn has_mane_select(&self) -> bool {
        self.mane_select.is_some()
    }
    fn has_mane_plus_clinical(&self) -> bool {
        self.mane_plus_clinical.is_some()
    }
    fn appris(&self) -> Option<&str> {
        self.appris.as_deref()
    }
    fn tsl(&self) -> Option<u8> {
        self.tsl
    }
    fn biotype(&self) -> &str {
        &self.biotype
    }
    fn has_ccds(&self) -> bool {
        self.ccds.is_some()
    }
    fn allele_count(&self) -> usize {
        self.allele_annotations.len()
    }
    fn allele(&self, index: usize) -> &Allele {
        &self.allele_annotations[index].allele
    }
    fn best_rank(&self, allele: Option<&Allele>) -> u32 {
        self.allele_annotations
            .iter()
            .filter(|aa| allele.is_none_or(|a| &aa.allele == a))
            .flat_map(|aa| aa.consequences.iter())
            .map(|c| c.rank())
            .min()
            .unwrap_or(u32::MAX)
    }
}

/// Whether an allele makes a row of a transcript: `*` (a deletion that began upstream) is not a
/// variant of the transcript, and VEP writes no row for it.
///
/// The one definition of it, for every place that has to number the alleles of a row the way the
/// pick does: [`row_alleles`], and the caller that numbers the alleles of a transcript the change
/// is beside before it has a prediction to read them from. If two of them disagreed the pick's
/// allele numbers would point at other alleles than the ones the rows have.
pub fn makes_row(allele: &Allele) -> bool {
    *allele != Allele::Missing
}

/// The alleles of a prediction that make a row, in the order the row writes them.
pub fn row_alleles(
    prediction: &TranscriptConsequence,
) -> impl Iterator<Item = &AlleleConsequenceResult> {
    prediction
        .allele_consequences
        .iter()
        .filter(|ac| makes_row(&ac.allele))
}

/// What the predictor said of one transcript, or, for a transcript the change is beside, what it
/// would say, read as the row a pick can choose from before that row is built.
///
/// It answers every question of [`PickRow`] as the [`TranscriptVariation`] built from the same
/// prediction does: the ids, the canonical flag and the biotype are the prediction's, the MANE,
/// APPRIS, TSL and CCDS attributes are the transcript's, and a prediction whose transcript could
/// not be found has none of them.
pub enum PredictedRow<'a> {
    /// A transcript predicted in full.
    Built {
        prediction: &'a TranscriptConsequence,
        transcript: Option<&'a Transcript>,
    },
    /// A transcript the change is beside: every allele in `alleles` has the one answer `flank`
    /// for it, which is all a pick can ask, so nothing has been built for it. The ids, the
    /// canonical flag and the biotype are read from the transcript, which is where the
    /// predictor takes them from.
    Flank {
        transcript: &'a Transcript,
        flank: Flank,
        alleles: &'a [&'a Allele],
    },
}

impl PredictedRow<'_> {
    fn transcript(&self) -> Option<&Transcript> {
        match self {
            Self::Built { transcript, .. } => *transcript,
            Self::Flank { transcript, .. } => Some(transcript),
        }
    }
}

impl PickRow for PredictedRow<'_> {
    fn transcript_id(&self) -> &str {
        match self {
            Self::Built { prediction, .. } => &prediction.transcript_id,
            Self::Flank { transcript, .. } => &transcript.stable_id,
        }
    }
    fn gene_id(&self) -> &str {
        match self {
            Self::Built { prediction, .. } => &prediction.gene_id,
            Self::Flank { transcript, .. } => &transcript.gene.stable_id,
        }
    }
    fn canonical(&self) -> bool {
        match self {
            Self::Built { prediction, .. } => prediction.canonical,
            Self::Flank { transcript, .. } => transcript.canonical,
        }
    }
    fn has_mane_select(&self) -> bool {
        self.transcript().is_some_and(|t| t.mane_select.is_some())
    }
    fn has_mane_plus_clinical(&self) -> bool {
        self.transcript()
            .is_some_and(|t| t.mane_plus_clinical.is_some())
    }
    fn appris(&self) -> Option<&str> {
        self.transcript().and_then(|t| t.appris.as_deref())
    }
    fn tsl(&self) -> Option<u8> {
        self.transcript().and_then(|t| t.tsl)
    }
    fn biotype(&self) -> &str {
        match self {
            Self::Built { prediction, .. } => &prediction.biotype,
            Self::Flank { transcript, .. } => &transcript.biotype,
        }
    }
    fn has_ccds(&self) -> bool {
        self.transcript().is_some_and(|t| t.ccds.is_some())
    }
    fn allele_count(&self) -> usize {
        match self {
            Self::Built { prediction, .. } => row_alleles(prediction).count(),
            Self::Flank { alleles, .. } => alleles.len(),
        }
    }
    fn allele(&self, index: usize) -> &Allele {
        match self {
            Self::Built { prediction, .. } => {
                &row_alleles(prediction)
                    .nth(index)
                    .expect("an allele number below `allele_count`")
                    .allele
            }
            Self::Flank { alleles, .. } => alleles[index],
        }
    }
    fn best_rank(&self, allele: Option<&Allele>) -> u32 {
        match self {
            Self::Built { prediction, .. } => row_alleles(prediction)
                .filter(|ac| allele.is_none_or(|a| &ac.allele == a))
                .flat_map(|ac| ac.consequences.iter())
                .map(|c| c.rank())
                .min()
                .unwrap_or(u32::MAX),
            // each allele has the one consequence, so the best of them is that one if any allele
            // is asked about
            Self::Flank { flank, alleles, .. } => {
                if alleles
                    .iter()
                    .any(|a| allele.is_none_or(|asked| *a == asked))
                {
                    flank.consequence.rank()
                } else {
                    u32::MAX
                }
            }
        }
    }
}

/// Run a pick over one variant's transcript variations, in place.
///
/// The entry point for finished rows: the library's loop, and the CLI's for a flagging pick, which
/// keeps every row. A reducing pick in the CLI is settled before the rows exist, with
/// [`pick_winners`] on [`PredictedRow`]s, and applied by what it builds; the two are held to one
/// answer by the tests named in the module header.
///
/// There used to be a `has_transcripts_to_pick` gate here, returning early
/// unless the record held more than one row and at least one real transcript.
/// It was written for the reducing modes, where collapsing a set of
/// per-allele scaffold rows loses an alt, and it was wrong for the flagging
/// ones in a way that a `PICK is 1` filter turns into missing variants: a
/// variant overlapping exactly *one* transcript failed `len > 1`, so
/// `--flag-pick*` declared the column and left it empty on the only entry
/// there was. Ensembl VEP 115.1 flags that entry (measured on a
/// single-transcript GFF3: `PICK=1` for all three `--flag_pick*` forms), and
/// so does this now. The scaffold hazard the gate existed for is handled
/// where it belongs, in [`pick_winners`], for all three scopes at once.
pub fn apply_pick(tvs: &mut Vec<TranscriptVariation>, plan: &PickPlan<'_>) {
    let winners = pick_winners(tvs, plan);
    match plan.mode {
        PickMode::Reduce => {
            // Nothing to keep means nothing was a candidate, which is not the
            // same as "every row lost": a site whose only rows are scaffolds
            // must come through untouched.
            if !winners.is_empty() {
                retain_winners(tvs, &winners);
            }
        }
        PickMode::Flag => {
            for &(ti, ai) in &winners {
                tvs[ti].allele_annotations[ai].pick = true;
            }
        }
    }
}

/// The `(row index, allele index)` pairs one pick chooses.
///
/// Candidates are the rows carrying a real transcript, at every scope. A
/// scaffold row (`transcript_id` `-`) is never a candidate, so it is never
/// flagged - it has no transcript, so "the transcript the pick chose" is not a
/// thing it can be - and never deleted, because it carries an alt allele that
/// nothing else in the record reports.
///
/// Public for the caller that settles the pick on [`PredictedRow`]s, before the rows exist, and
/// then builds only the winners.
pub fn pick_winners<R: PickRow>(rows: &[R], plan: &PickPlan<'_>) -> Vec<(usize, usize)> {
    match plan.scope {
        PickScope::Variant => {
            let Some(best) = best_real_transcript_idx(rows, plan.order) else {
                return Vec::new();
            };
            // Every allele of the winning transcript, which is what makes
            // `--pick` keep an alt that VEP's `--pick` drops. See `PickScope`.
            (0..rows[best].allele_count())
                .map(|ai| (best, ai))
                .collect()
        }
        PickScope::Allele | PickScope::AlleleGene => grouped_winners(rows, plan),
    }
}

/// Index of the best row carrying a real transcript.
///
/// `pick_best_transcript_idx_with` scores whatever it is given, and `"-"` sorts
/// before every real transcript ID in the final tie-break - so at a site
/// holding both a scaffold and a transcript, an order whose listed criteria all
/// tie would pick the scaffold, and a reducing pick would then delete the real
/// annotation. `--pick-order mane_select` is enough to produce that tie.
fn best_real_transcript_idx<R: PickRow>(rows: &[R], order: &[PickCriterion]) -> Option<usize> {
    (0..rows.len())
        .filter(|&i| !rows[i].is_featureless())
        .min_by(|&a, &b| {
            pick_key_with(&rows[a], order, None).cmp(&pick_key_with(&rows[b], order, None))
        })
}

/// The best `(row index, allele index)` pair in each group, where a group is
/// one allele or one (allele, gene) pair.
///
/// One pass. A row finds its group through numbers rather than comparisons: its gene is hashed to a
/// dense number once, its allele is looked up among the few distinct alleles of the site, and the
/// pair of them finds the group. Groups used to be found by comparing the row's allele and gene
/// against every group before it, which is short work while a locus has a handful of genes. At
/// `--distance 500000` there are dozens of genes around a variant and each row compared its gene's
/// text with each of their groups in turn: the comparisons and the `memcmp` under them were about a
/// fifth of the CPU time of a run over HG002 chr1 with the flags of #166, in a sampling profile.
/// Nothing here builds a `String`: the keys hold borrows of the rows.
///
/// The tables are sized for the rows up front, so that none of them grows: this runs once per
/// variant, and a table that starts small and doubles is five allocations where one would do.
fn grouped_winners<R: PickRow>(rows: &[R], plan: &PickPlan<'_>) -> Vec<(usize, usize)> {
    struct Group<'a> {
        best: (usize, usize),
        key: PickKey<'a>,
    }

    let by_gene = plan.scope == PickScope::AlleleGene;
    // A gene's number is the order in which the rows first met it. Not used when the pick is not
    // per gene, which puts every row in the one gene `0`.
    let mut gene_numbers: HashMap<&str, usize, FixedState> = HashMap::with_capacity_and_hasher(
        if by_gene { rows.len() } else { 0 },
        FixedState::default(),
    );
    let mut alleles: Vec<&Allele> = Vec::new();
    // The number of the group a (gene number, allele number) pair is in.
    let mut group_of: HashMap<(usize, usize), usize, FixedState> =
        HashMap::with_capacity_and_hasher(rows.len(), FixedState::default());
    let mut groups: Vec<Group<'_>> = Vec::with_capacity(rows.len());

    for (ti, row) in rows.iter().enumerate() {
        // A scaffold row carries no transcript, so it cannot win a pick and
        // must not be able to lose one either - it would take an allele with
        // it under Reduce.
        if row.is_featureless() {
            continue;
        }
        let gene = if by_gene {
            let next = gene_numbers.len();
            *gene_numbers.entry(row.gene_id()).or_insert(next)
        } else {
            0
        };
        for ai in 0..row.allele_count() {
            let allele = row.allele(ai);
            let key = pick_key_with(row, plan.order, Some(allele));
            let allele_number = match alleles.iter().position(|a| *a == allele) {
                Some(number) => number,
                None => {
                    alleles.push(allele);
                    alleles.len() - 1
                }
            };
            match group_of.entry((gene, allele_number)) {
                std::collections::hash_map::Entry::Occupied(entry) => {
                    let group = &mut groups[*entry.get()];
                    if key < group.key {
                        group.best = (ti, ai);
                        group.key = key;
                    }
                }
                std::collections::hash_map::Entry::Vacant(entry) => {
                    entry.insert(groups.len());
                    groups.push(Group {
                        best: (ti, ai),
                        key,
                    });
                }
            }
        }
    }

    groups.into_iter().map(|g| g.best).collect()
}

/// Drop every (transcript, allele) pair that is not a winner, keeping the
/// original transcript and allele order.
///
/// A transcript left with no alleles is dropped, and a scaffold row is kept
/// untouched: it was never in a group, and deleting it would delete an allele
/// nothing else reports.
fn retain_winners(tvs: &mut Vec<TranscriptVariation>, winners: &[(usize, usize)]) {
    // Sorted, so that asking whether a pair won is a search, and not a scan of every winner for
    // every pair of every row.
    let mut winners = winners.to_vec();
    winners.sort_unstable();
    let mut out: Vec<TranscriptVariation> = Vec::with_capacity(winners.len());
    for (ti, mut tv) in std::mem::take(tvs).into_iter().enumerate() {
        if tv.is_featureless() {
            out.push(tv);
            continue;
        }
        let anns = std::mem::take(&mut tv.allele_annotations);
        let kept: Vec<_> = anns
            .into_iter()
            .enumerate()
            .filter(|(ai, _)| winners.binary_search(&(ti, *ai)).is_ok())
            .map(|(_, aa)| aa)
            .collect();
        if kept.is_empty() {
            continue;
        }
        tv.allele_annotations = kept;
        out.push(tv);
    }
    *tvs = out;
}

/// Index of the best transcript variation under the given `--pick-order`
/// hierarchy, with transcript_id alphabetical order as a final deterministic
/// tie-breaker.
///
/// Criteria omitted from `order` are not consulted at all, which is what lets
/// a caller drop a tier rather than only reorder it.
pub fn pick_best_transcript_idx_with<R: PickRow>(
    rows: &[R],
    order: &[PickCriterion],
) -> Option<usize> {
    (0..rows.len()).min_by(|&a, &b| {
        pick_key_with(&rows[a], order, None).cmp(&pick_key_with(&rows[b], order, None))
    })
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
///
/// `allele` scopes the `Rank` tier. Every other criterion is a property of the
/// transcript and ignores it, but severity is a property of the *change*: a
/// transcript that a site's first alt makes nonsense and its second alt makes
/// synonymous has two ranks, and an allele-scoped pick has to compare the one
/// belonging to the allele it is picking for. `None` means "every allele on
/// this transcript", which is what a variant-scoped pick wants.
fn pick_score<R: PickRow>(row: &R, c: PickCriterion, allele: Option<&Allele>) -> u32 {
    match c {
        PickCriterion::ManeSelect => !row.has_mane_select() as u32,
        PickCriterion::ManePlusClinical => !row.has_mane_plus_clinical() as u32,
        PickCriterion::Canonical => !row.canonical() as u32,
        PickCriterion::Appris => appris_rank(row.appris()),
        PickCriterion::Tsl => row.tsl().unwrap_or(u8::MAX) as u32,
        PickCriterion::Biotype => u32::from(row.biotype() != "protein_coding"),
        PickCriterion::Ccds => !row.has_ccds() as u32,
        PickCriterion::Rank => row.best_rank(allele),
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
fn pick_key_with<'a, R: PickRow>(
    row: &'a R,
    order: &[PickCriterion],
    allele: Option<&Allele>,
) -> PickKey<'a> {
    let mut key = [0u32; MAX_PICK_CRITERIA];
    for (slot, &c) in key.iter_mut().zip(order.iter()) {
        *slot = pick_score(row, c, allele);
    }
    (key, row.transcript_id())
}

/// The scores in `--pick-order` order, then the transcript ID as the final
/// deterministic tie-break.
type PickKey<'a> = ([u32; MAX_PICK_CRITERIA], &'a str);

/// Every APPRIS spelling reaches here, so the bands are spaced rather than
/// adjacent: a tier number is added to its band's base, and an un-numbered
/// call sits at the top of its own band.
///
/// APPRIS itself only issues principal 1-5 and alternative 1-2, so the 99
/// slots per band are slack, not a claim about the vocabulary.
const APPRIS_PRINCIPAL_BASE: u32 = 0;
const APPRIS_PRINCIPAL_UNNUMBERED: u32 = 100;
const APPRIS_ALTERNATIVE_BASE: u32 = 200;
const APPRIS_ALTERNATIVE_UNNUMBERED: u32 = 300;
/// Present but in no spelling this understands. Worse than any call it does
/// understand, better than no call at all - the transcript was annotated.
const APPRIS_UNRECOGNISED: u32 = u32::MAX - 1;
/// No APPRIS annotation. Ensembl's own GFF3 carries none for any transcript, so
/// on that source this is every transcript's rank and the tier decides nothing,
/// correctly, because there is nothing to decide with. GENCODE's GFF3 does
/// carry it, and before the parser read those tags this was every transcript's
/// rank on every source.
const APPRIS_ABSENT: u32 = u32::MAX;

/// Map an APPRIS tag to a rank where lower is better, matching VEP's
/// `--pick_order` APPRIS tier:
/// principal1 < ... < principal5 < principal < alternative1 < alternative2 <
/// alternative < unrecognised < absent.
///
/// Four spellings of the same thing reach this function, because three sources
/// write it differently and the separator moved between GENCODE releases:
/// GENCODE's own tag (`appris_principal_1`), the tag with the `appris_` prefix
/// already stripped (`principal_1`, `principal1`), and VEP's short form (`P1`,
/// `ALT1`). Handling only some of them is how `appris_principal_1` used to rank
/// *below* every alternative: it fell through to the bare-`a` arm.
fn appris_rank(appris: Option<&str>) -> u32 {
    let Some(s) = appris else {
        return APPRIS_ABSENT;
    };
    // Compared without regard to case, and without building a lower-cased copy: this is called for
    // every row of every allele of every variant, and a `String` for each of them is an allocation
    // that a run with GENCODE's tags and the default `--pick-order` made millions of times.
    let s = s.trim();
    let body = strip_prefix_ignore_ascii_case(s, "appris_").unwrap_or(s);
    // `alt` before the bare `a`, and `alternative` before `alt`: the shorter
    // prefixes would otherwise swallow the longer spellings' digits.
    let (base, unnumbered, digits) =
        if let Some(d) = strip_prefix_ignore_ascii_case(body, "principal") {
            (APPRIS_PRINCIPAL_BASE, APPRIS_PRINCIPAL_UNNUMBERED, d)
        } else if let Some(d) = strip_prefix_ignore_ascii_case(body, "alternative") {
            (APPRIS_ALTERNATIVE_BASE, APPRIS_ALTERNATIVE_UNNUMBERED, d)
        } else if let Some(d) = strip_prefix_ignore_ascii_case(body, "alt") {
            (APPRIS_ALTERNATIVE_BASE, APPRIS_ALTERNATIVE_UNNUMBERED, d)
        } else if let Some(d) = strip_prefix_ignore_ascii_case(body, "p") {
            (APPRIS_PRINCIPAL_BASE, APPRIS_PRINCIPAL_UNNUMBERED, d)
        } else if let Some(d) = strip_prefix_ignore_ascii_case(body, "a") {
            (APPRIS_ALTERNATIVE_BASE, APPRIS_ALTERNATIVE_UNNUMBERED, d)
        } else {
            return APPRIS_UNRECOGNISED;
        };
    let digits = digits.trim_start_matches(['_', '-']);
    if digits.is_empty() {
        return unnumbered;
    }
    match digits.parse::<u32>() {
        // Clamped so an out-of-vocabulary tier still sorts inside its own
        // band rather than past the next one.
        Ok(n) => base + n.clamp(1, 99),
        // `P1extra`: the class is legible, the tier is not.
        Err(_) => unnumbered,
    }
}

/// `s` without `prefix`, when it starts with it in any case of ASCII letters.
///
/// `get` is `None` for a prefix longer than `s` and for a cut that falls inside a multi-byte
/// character, and neither can be a match, because `prefix` is ASCII.
fn strip_prefix_ignore_ascii_case<'a>(s: &'a str, prefix: &str) -> Option<&'a str> {
    let head = s.get(..prefix.len())?;
    head.eq_ignore_ascii_case(prefix)
        .then(|| &s[prefix.len()..])
}

#[cfg(test)]
mod pick_tests {
    use super::*;
    use fastvep_core::{Allele, Consequence, Impact, Strand, TranscriptLocation};
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
                location: TranscriptLocation::default(),
                protein_position: None,
                amino_acids: None,
                codons: None,
                exon: None,
                splice: Default::default(),
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
                pick: false,
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
    fn a_scaffold_only_site_keeps_every_allele_at_every_scope() {
        // A multi-allelic intergenic site arrives as one placeholder row per
        // alt allele, so a pick that treated them as competing transcripts
        // would keep one and silently lose the other alt - no error, just a
        // missing allele. Nothing here is a candidate, so nothing is picked.
        for scope in [PickScope::Variant, PickScope::Allele, PickScope::AlleleGene] {
            for mode in [PickMode::Reduce, PickMode::Flag] {
                let mut rows: Vec<TranscriptVariation> = ["-", "-"]
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
                apply_pick(
                    &mut rows,
                    &PickPlan {
                        scope,
                        mode,
                        order: DEFAULT_PICK_ORDER,
                    },
                );
                assert_eq!(rows.len(), 2, "{scope:?}/{mode:?} dropped a scaffold row");
                assert!(
                    rows.iter()
                        .all(|tv| tv.allele_annotations.iter().all(|aa| !aa.pick)),
                    "{scope:?}/{mode:?} flagged a row with no transcript"
                );
            }
        }
    }

    #[test]
    fn a_scaffold_cannot_outrank_a_real_transcript() {
        // `"-"` sorts before every real transcript ID, which is the final
        // tie-break, so an order whose listed criteria all tie used to hand
        // the pick to the placeholder - and a reducing pick then deleted the
        // real annotation. `mane_select` alone is such an order: neither row
        // has one.
        let order = parse_pick_order("mane_select").unwrap();
        for mode in [PickMode::Reduce, PickMode::Flag] {
            let mut rows = vec![
                make_tv(
                    "-",
                    false,
                    "-",
                    vec![Consequence::IntergenicVariant],
                    None,
                    None,
                    None,
                    None,
                    None,
                ),
                make_tv(
                    "ENST1",
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
            apply_pick(
                &mut rows,
                &PickPlan {
                    scope: PickScope::Variant,
                    mode,
                    order: &order,
                },
            );
            assert!(
                rows.iter().any(|tv| tv.transcript_id.as_ref() == "ENST1"),
                "{mode:?} deleted the real transcript in favour of a scaffold"
            );
            let flagged: Vec<&str> = rows
                .iter()
                .filter(|tv| tv.allele_annotations.iter().any(|aa| aa.pick))
                .map(|tv| tv.transcript_id.as_ref())
                .collect();
            if mode == PickMode::Flag {
                assert_eq!(flagged, vec!["ENST1"], "the transcript is the pick");
            }
        }
    }

    #[test]
    fn a_lone_transcript_is_still_the_pick() {
        // The bug this replaced a `len > 1` gate for: with one overlapping
        // transcript there is nothing to choose between, but there is still
        // something to *flag*, and a client filtering on `PICK is 1` drops the
        // variant entirely if it is left blank. VEP 115.1 flags it.
        for scope in [PickScope::Variant, PickScope::Allele, PickScope::AlleleGene] {
            let mut rows = vec![tv_in_gene(
                "ENST_ONLY",
                "GENE",
                true,
                &[("A", Consequence::MissenseVariant)],
            )];
            apply_pick(
                &mut rows,
                &PickPlan {
                    scope,
                    mode: PickMode::Flag,
                    order: DEFAULT_PICK_ORDER,
                },
            );
            assert_eq!(
                rendered(&rows),
                vec![("ENST_ONLY".into(), "A".into(), true)],
                "{scope:?} left the only entry unflagged"
            );
        }
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
    fn appris_orders_alternatives_by_their_tier() {
        // ALT1 and ALT2 both reach the ranker as the short form the GFF3
        // parser now writes, and the pre-fix ranker scored them equal: the
        // bare-`a` arm swallowed the `A`, `lt1` and `lt2` both failed to
        // parse, and both landed on the same fallback tier.
        //
        // The IDs are chosen so that a tie is not silently right. `TX_A`
        // carries the *worse* call, so the alphabetical transcript-ID
        // tie-break returns it whenever the APPRIS tier declines to decide -
        // which is what the equal scores used to produce.
        let tvs = vec![
            make_tv(
                "TX_A",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                Some("ALT2"),
                None,
                None,
            ),
            make_tv(
                "TX_Z",
                false,
                "protein_coding",
                vec![Consequence::MissenseVariant],
                None,
                None,
                Some("ALT1"),
                None,
                None,
            ),
        ];
        assert_eq!(
            pick_best_transcript_idx(&tvs),
            Some(1),
            "ALT1 should outrank ALT2; index 0 is what the ID tie-break gives"
        );
    }

    #[test]
    fn appris_reads_every_spelling_of_the_same_call() {
        // Three sources write the same annotation three ways, and the
        // separator moved between GENCODE releases. All four spellings of
        // principal 1 have to rank identically, and all of them ahead of
        // every spelling of alternative 2.
        for principal in [
            "P1",
            "p1",
            "principal1",
            "principal_1",
            "appris_principal_1",
        ] {
            for alternative in [
                "ALT2",
                "A2",
                "alternative2",
                "alternative_2",
                "appris_alternative_2",
            ] {
                assert!(
                    appris_rank(Some(principal)) < appris_rank(Some(alternative)),
                    "{principal} ({}) should outrank {alternative} ({})",
                    appris_rank(Some(principal)),
                    appris_rank(Some(alternative)),
                );
            }
        }
        // And the numbered forms agree with each other exactly.
        let p1 = appris_rank(Some("P1"));
        for spelling in ["p1", "principal1", "principal_1", "appris_principal_1"] {
            assert_eq!(appris_rank(Some(spelling)), p1, "{spelling}");
        }
    }

    /// `appris_rank` as it was written when it lower-cased a copy of its argument. It is the oracle
    /// of the test below, which holds the version that builds no copy to the same answers.
    fn appris_rank_through_a_lower_cased_copy(appris: Option<&str>) -> u32 {
        let Some(s) = appris else {
            return APPRIS_ABSENT;
        };
        let lower = s.trim().to_ascii_lowercase();
        let body = lower.strip_prefix("appris_").unwrap_or(&lower);
        let (base, unnumbered, digits) = if let Some(d) = body.strip_prefix("principal") {
            (APPRIS_PRINCIPAL_BASE, APPRIS_PRINCIPAL_UNNUMBERED, d)
        } else if let Some(d) = body.strip_prefix("alternative") {
            (APPRIS_ALTERNATIVE_BASE, APPRIS_ALTERNATIVE_UNNUMBERED, d)
        } else if let Some(d) = body.strip_prefix("alt") {
            (APPRIS_ALTERNATIVE_BASE, APPRIS_ALTERNATIVE_UNNUMBERED, d)
        } else if let Some(d) = body.strip_prefix('p') {
            (APPRIS_PRINCIPAL_BASE, APPRIS_PRINCIPAL_UNNUMBERED, d)
        } else if let Some(d) = body.strip_prefix('a') {
            (APPRIS_ALTERNATIVE_BASE, APPRIS_ALTERNATIVE_UNNUMBERED, d)
        } else {
            return APPRIS_UNRECOGNISED;
        };
        let digits = digits.trim_start_matches(['_', '-']);
        if digits.is_empty() {
            return unnumbered;
        }
        match digits.parse::<u32>() {
            Ok(n) => base + n.clamp(1, 99),
            Err(_) => unnumbered,
        }
    }

    #[test]
    fn appris_rank_reads_a_tag_in_any_case_as_it_did_through_a_lower_cased_copy() {
        // Every string of up to three of these pieces: the prefixes in every case, the separators,
        // tiers, a sign `parse` accepts, whitespace `trim` removes, and characters that are not
        // ASCII, including ones whose cut would fall inside a character.
        const PIECES: [&str; 26] = [
            "",
            " ",
            "\t",
            "appris_",
            "APPRIS_",
            "Appris_",
            "principal",
            "PRINCIPAL",
            "Principal",
            "alternative",
            "ALTERNATIVE",
            "alt",
            "ALT",
            "a",
            "A",
            "p",
            "P",
            "_",
            "-",
            "1",
            "2",
            "007",
            "+3",
            "\u{e9}",
            "\u{212a}",
            "\u{130}",
        ];
        let mut checked = 0usize;
        let mut each = |tag: &str| {
            assert_eq!(
                appris_rank(Some(tag)),
                appris_rank_through_a_lower_cased_copy(Some(tag)),
                "{tag:?}"
            );
            checked += 1;
        };
        each("");
        for a in PIECES {
            each(a);
            for b in PIECES {
                each(&format!("{a}{b}"));
                for c in PIECES {
                    each(&format!("{a}{b}{c}"));
                }
            }
        }
        assert!(checked > 17_000, "{checked} tags");
        assert_eq!(
            appris_rank(None),
            appris_rank_through_a_lower_cased_copy(None)
        );
    }

    #[test]
    fn appris_ranks_an_unnumbered_call_at_the_end_of_its_own_class() {
        // Older GENCODE releases write `appris_principal` with no tier. It is
        // still a principal call: worse than every numbered principal, better
        // than any alternative.
        let (p1, p5) = (appris_rank(Some("P1")), appris_rank(Some("P5")));
        let (p, alt1) = (appris_rank(Some("P")), appris_rank(Some("ALT1")));
        let (alt2, alt) = (appris_rank(Some("ALT2")), appris_rank(Some("ALT")));
        assert!(p1 < p5, "P1 < P5");
        assert!(p5 < p, "P5 < P");
        assert!(p < alt1, "P < ALT1");
        assert!(alt1 < alt2, "ALT1 < ALT2");
        assert!(alt2 < alt, "ALT2 < ALT");
        assert!(alt < appris_rank(Some("weird")), "ALT < unrecognised");
        assert!(
            appris_rank(Some("weird")) < appris_rank(None),
            "unrecognised < absent"
        );
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

    /// A transcript with an explicit gene and one annotation per allele, for
    /// the scope tests. `make_tv` hardcodes one gene and one allele, which is
    /// all the tier tests need and neither of the things a scope decides on.
    fn tv_in_gene(
        transcript_id: &str,
        gene_id: &str,
        canonical: bool,
        per_allele: &[(&str, Consequence)],
    ) -> TranscriptVariation {
        let mut tv = make_tv(
            transcript_id,
            canonical,
            "protein_coding",
            vec![Consequence::IntronVariant],
            None,
            None,
            None,
            None,
            None,
        );
        tv.gene_id = Arc::from(gene_id);
        tv.gene_symbol = Some(Arc::from(gene_id));
        tv.allele_annotations = per_allele
            .iter()
            .map(|(allele, csq)| {
                let mut aa = tv.allele_annotations[0].clone();
                aa.allele = Allele::from_str(allele);
                aa.consequences = vec![*csq];
                aa
            })
            .collect();
        tv
    }

    /// `(transcript, allele, flagged)` for every annotation left standing, in
    /// output order - which is what every writer iterates.
    fn rendered(tvs: &[TranscriptVariation]) -> Vec<(String, String, bool)> {
        tvs.iter()
            .flat_map(|tv| {
                tv.allele_annotations
                    .iter()
                    .map(move |aa| (tv.transcript_id.to_string(), aa.allele.to_string(), aa.pick))
            })
            .collect()
    }

    /// Two genes at one locus, each with a canonical and a non-canonical
    /// transcript, and two alt alleles whose severity ordering is opposite
    /// between the genes - so a scope that ignores the allele cannot get the
    /// per-allele answers right by accident.
    fn two_gene_locus() -> Vec<TranscriptVariation> {
        vec![
            tv_in_gene(
                "ENST_A_CANON",
                "GENE_A",
                true,
                &[
                    ("A", Consequence::MissenseVariant),
                    ("C", Consequence::SynonymousVariant),
                ],
            ),
            tv_in_gene(
                "ENST_A_OTHER",
                "GENE_A",
                false,
                &[
                    ("A", Consequence::StopGained),
                    ("C", Consequence::StopGained),
                ],
            ),
            tv_in_gene(
                "ENST_B_CANON",
                "GENE_B",
                true,
                &[
                    ("A", Consequence::SynonymousVariant),
                    ("C", Consequence::MissenseVariant),
                ],
            ),
        ]
    }

    #[test]
    fn variant_scope_reduce_keeps_one_transcript_with_all_its_alleles() {
        // fastVEP's deliberate divergence: VEP's `--pick` emits one entry for
        // the record, which silently drops an alt. Measured on VEP 115.1 at a
        // two-alt TP53 site: 1 entry against fastVEP's 2. See `PickScope`.
        let mut tvs = two_gene_locus();
        apply_pick(&mut tvs, &PickPlan::reduce_per_variant(DEFAULT_PICK_ORDER));
        assert_eq!(
            rendered(&tvs),
            vec![
                ("ENST_A_CANON".into(), "A".into(), false),
                ("ENST_A_CANON".into(), "C".into(), false),
            ]
        );
    }

    #[test]
    fn allele_scope_picks_per_allele_across_genes() {
        let mut tvs = two_gene_locus();
        apply_pick(
            &mut tvs,
            &PickPlan {
                scope: PickScope::Allele,
                mode: PickMode::Reduce,
                order: DEFAULT_PICK_ORDER,
            },
        );
        // Canonical eliminates GENE_A's other transcript for both alleles,
        // even though it carries the most severe consequence at the locus.
        // That leaves the two canonical transcripts tied on every tier above
        // `rank`, so `rank` decides - and it decides *per allele*, because
        // severity is a property of the change and not of the transcript.
        // Each allele therefore goes to a different gene, which is the whole
        // difference between this scope and the variant one.
        assert_eq!(
            rendered(&tvs),
            vec![
                ("ENST_A_CANON".into(), "A".into(), false),
                ("ENST_B_CANON".into(), "C".into(), false),
            ]
        );
    }

    #[test]
    fn allele_gene_scope_picks_once_per_allele_and_gene() {
        let mut tvs = two_gene_locus();
        apply_pick(
            &mut tvs,
            &PickPlan {
                scope: PickScope::AlleleGene,
                mode: PickMode::Reduce,
                order: DEFAULT_PICK_ORDER,
            },
        );
        // Four groups: two alleles x two genes. GENE_A's non-canonical
        // transcript loses both of its, even though it carries the most severe
        // consequence at the locus - canonical outranks rank, as in VEP.
        assert_eq!(
            rendered(&tvs),
            vec![
                ("ENST_A_CANON".into(), "A".into(), false),
                ("ENST_A_CANON".into(), "C".into(), false),
                ("ENST_B_CANON".into(), "A".into(), false),
                ("ENST_B_CANON".into(), "C".into(), false),
            ]
        );
    }

    #[test]
    fn flag_mode_marks_the_same_choices_it_would_have_reduced_to() {
        // The property that makes `--flag-pick*` trustworthy: flagging and
        // reducing must not disagree about which entry wins, or the flag means
        // something different from the option it is named after.
        for scope in [PickScope::Variant, PickScope::Allele, PickScope::AlleleGene] {
            let mut reduced = two_gene_locus();
            apply_pick(
                &mut reduced,
                &PickPlan {
                    scope,
                    mode: PickMode::Reduce,
                    order: DEFAULT_PICK_ORDER,
                },
            );
            let mut flagged = two_gene_locus();
            apply_pick(
                &mut flagged,
                &PickPlan {
                    scope,
                    mode: PickMode::Flag,
                    order: DEFAULT_PICK_ORDER,
                },
            );

            let kept: Vec<_> = rendered(&reduced)
                .into_iter()
                .map(|(t, a, _)| (t, a))
                .collect();
            let marked: Vec<_> = rendered(&flagged)
                .into_iter()
                .filter(|(_, _, picked)| *picked)
                .map(|(t, a, _)| (t, a))
                .collect();
            assert_eq!(
                kept, marked,
                "{scope:?} flags a different set than it keeps"
            );

            // And flagging keeps everything, which is the whole point.
            assert_eq!(rendered(&flagged).len(), 6, "{scope:?} dropped an entry");
        }
    }

    #[test]
    fn allele_scope_consults_the_rank_of_the_allele_being_picked_for() {
        // Two transcripts of one gene, neither canonical, so `rank` is the
        // only tier that can separate them - and it points a different way for
        // each allele. Scoring `rank` over every allele at once, as a
        // variant-scoped pick does, would give both alleles the same winner.
        let mut tvs = vec![
            tv_in_gene(
                "ENST_X",
                "GENE",
                false,
                &[
                    ("A", Consequence::StopGained),
                    ("C", Consequence::SynonymousVariant),
                ],
            ),
            tv_in_gene(
                "ENST_Y",
                "GENE",
                false,
                &[
                    ("A", Consequence::SynonymousVariant),
                    ("C", Consequence::StopGained),
                ],
            ),
        ];
        apply_pick(
            &mut tvs,
            &PickPlan {
                scope: PickScope::Allele,
                mode: PickMode::Reduce,
                order: DEFAULT_PICK_ORDER,
            },
        );
        assert_eq!(
            rendered(&tvs),
            vec![
                ("ENST_X".into(), "A".into(), false),
                ("ENST_Y".into(), "C".into(), false),
            ]
        );
    }

    #[test]
    fn a_scaffold_row_keeps_its_allele_through_every_scope() {
        // A placeholder carries an alt that nothing else in the record
        // reports, so a reducing pick at a grouping scope must not be able to
        // delete it while keeping the real transcript beside it.
        for scope in [PickScope::Allele, PickScope::AlleleGene] {
            let mut tvs = vec![
                tv_in_gene(
                    "ENST_REAL",
                    "GENE",
                    true,
                    &[("A", Consequence::MissenseVariant)],
                ),
                tv_in_gene("-", "-", false, &[("C", Consequence::IntergenicVariant)]),
            ];
            apply_pick(
                &mut tvs,
                &PickPlan {
                    scope,
                    mode: PickMode::Reduce,
                    order: DEFAULT_PICK_ORDER,
                },
            );
            let alleles: Vec<String> = rendered(&tvs).into_iter().map(|(_, a, _)| a).collect();
            assert!(
                alleles.contains(&"C".to_string()),
                "{scope:?} dropped the scaffolded allele: {alleles:?}"
            );
        }
    }

    #[test]
    fn pick_flags_resolve_to_the_six_vep_options() {
        assert_eq!(PickFlags::default().resolve(), None);
        for (flags, scope, mode, name) in [
            (
                PickFlags {
                    pick: true,
                    ..Default::default()
                },
                PickScope::Variant,
                PickMode::Reduce,
                "--pick",
            ),
            (
                PickFlags {
                    pick_allele: true,
                    ..Default::default()
                },
                PickScope::Allele,
                PickMode::Reduce,
                "--pick-allele",
            ),
            (
                PickFlags {
                    pick_allele_gene: true,
                    ..Default::default()
                },
                PickScope::AlleleGene,
                PickMode::Reduce,
                "--pick-allele-gene",
            ),
            (
                PickFlags {
                    flag_pick: true,
                    ..Default::default()
                },
                PickScope::Variant,
                PickMode::Flag,
                "--flag-pick",
            ),
            (
                PickFlags {
                    flag_pick_allele: true,
                    ..Default::default()
                },
                PickScope::Allele,
                PickMode::Flag,
                "--flag-pick-allele",
            ),
            (
                PickFlags {
                    flag_pick_allele_gene: true,
                    ..Default::default()
                },
                PickScope::AlleleGene,
                PickMode::Flag,
                "--flag-pick-allele-gene",
            ),
        ] {
            let request = flags.resolve().expect("a switch is set");
            assert_eq!(request.scope, scope, "{name}");
            assert_eq!(request.mode, mode, "{name}");
            assert_eq!(flags.requested_option(), Some(name));
            assert_eq!(
                request.needs_pick_column(),
                mode == PickMode::Flag,
                "{name}"
            );
        }
    }

    #[test]
    fn a_malformed_pick_switch_is_an_error_rather_than_no_pick() {
        // All six are read in one pass, so defaulting on failure would let a
        // wrong type on one key discard a correctly spelled sibling and
        // annotate with no pick at all - a wrong answer that looks right.
        let err = pick_flags_from_json(&serde_json::json!({
            "vcf": "ignored",
            "pick_allele": true,
            "pick": "yes",
        }))
        .expect_err("a non-bool switch must not be silently ignored");
        assert!(
            err.to_string().contains("boolean") || err.to_string().contains("bool"),
            "the error should name the problem: {err}"
        );

        // And a well-formed body still resolves, reading through a borrow.
        let body = serde_json::json!({ "vcf": "ignored", "pick_allele": true });
        assert_eq!(
            pick_flags_from_json(&body).unwrap().requested_option(),
            Some("--pick-allele")
        );
    }

    #[test]
    fn pick_switches_deserialize_under_their_wire_names() {
        // Both HTTP entry points read the request body through this, so the
        // names are part of the API, not an implementation detail.
        let body = serde_json::json!({
            "vcf": "ignored",
            "flag_pick_allele_gene": true,
        });
        let flags: PickFlags = serde_json::from_value(body).expect("unknown keys are ignored");
        assert_eq!(flags.requested_option(), Some("--flag-pick-allele-gene"));

        let empty: PickFlags = serde_json::from_value(serde_json::json!({})).unwrap();
        assert_eq!(empty.resolve(), None);
    }

    /// The grouping as it was written before a row found its group by number: every row compared
    /// with every group so far. Kept as the oracle for the faster one, which has to choose the same
    /// pair in every group, in the same order, on any input, ties included.
    fn reference_grouped_winners(
        tvs: &[TranscriptVariation],
        plan: &PickPlan<'_>,
    ) -> Vec<(usize, usize)> {
        struct Group<'a> {
            allele: &'a Allele,
            gene: Option<&'a str>,
            best: (usize, usize),
            key: PickKey<'a>,
        }
        let by_gene = plan.scope == PickScope::AlleleGene;
        let mut groups: Vec<Group<'_>> = Vec::new();
        for (ti, tv) in tvs.iter().enumerate() {
            if tv.transcript_id.as_ref() == "-" {
                continue;
            }
            let gene = by_gene.then(|| tv.gene_id.as_ref());
            for (ai, aa) in tv.allele_annotations.iter().enumerate() {
                let key = pick_key_with(tv, plan.order, Some(&aa.allele));
                match groups
                    .iter_mut()
                    .find(|g| g.allele == &aa.allele && g.gene == gene)
                {
                    Some(g) => {
                        if key < g.key {
                            g.best = (ti, ai);
                            g.key = key;
                        }
                    }
                    None => groups.push(Group {
                        allele: &aa.allele,
                        gene,
                        best: (ti, ai),
                        key,
                    }),
                }
            }
        }
        groups.into_iter().map(|g| g.best).collect()
    }

    /// A small deterministic generator, so that a failure is the same failure on every machine.
    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }
        fn chance(&mut self, one_in: usize) -> bool {
            self.below(one_in) == 0
        }
    }

    /// A locus of up to 40 transcripts in up to 12 genes, over up to three alt alleles, with every
    /// attribute a pick reads varied, transcript ids drawn from a pool small enough that some
    /// repeat (a tie on every tier), and now and then a scaffold row.
    fn random_locus(rng: &mut Rng) -> Vec<TranscriptVariation> {
        const ALLELES: [&str; 5] = ["A", "C", "G", "-", "<DEL>"];
        const CONSEQUENCES: [Consequence; 6] = [
            Consequence::MissenseVariant,
            Consequence::SynonymousVariant,
            Consequence::StopGained,
            Consequence::IntronVariant,
            Consequence::UpstreamGeneVariant,
            Consequence::DownstreamGeneVariant,
        ];
        const APPRIS: [&str; 6] = [
            "principal1",
            "P3",
            "alternative2",
            "ALT1",
            "appris_principal_1",
            "x",
        ];
        let mut tvs = Vec::new();
        for _ in 0..rng.below(41) {
            let id = if rng.chance(25) {
                "-".to_string()
            } else {
                format!("ENST{:02}", rng.below(30))
            };
            let consequence = CONSEQUENCES[rng.below(CONSEQUENCES.len())];
            let mut tv = make_tv(
                &id,
                rng.chance(2),
                if rng.chance(3) {
                    "lncRNA"
                } else {
                    "protein_coding"
                },
                vec![consequence],
                rng.chance(4).then_some("NM_1"),
                rng.chance(6).then_some("NM_2"),
                rng.chance(2).then(|| APPRIS[rng.below(APPRIS.len())]),
                rng.chance(2).then(|| rng.below(6) as u8),
                rng.chance(3).then_some("CCDS1"),
            );
            tv.gene_id = Arc::from(format!("GENE{:02}", rng.below(12)).as_str());
            let alleles = 1 + rng.below(3);
            let mut annotations = Vec::new();
            for _ in 0..alleles {
                let mut aa = tv.allele_annotations[0].clone();
                aa.allele = Allele::from_str(ALLELES[rng.below(ALLELES.len())]);
                aa.consequences = (0..=rng.below(2))
                    .map(|_| CONSEQUENCES[rng.below(CONSEQUENCES.len())])
                    .collect();
                annotations.push(aa);
            }
            tv.allele_annotations = annotations;
            tvs.push(tv);
        }
        tvs
    }

    fn random_order(rng: &mut Rng) -> Vec<PickCriterion> {
        let mut all = DEFAULT_PICK_ORDER.to_vec();
        // Fisher-Yates, then a prefix: orders of every length and every sequence
        for i in (1..all.len()).rev() {
            all.swap(i, rng.below(i + 1));
        }
        all.truncate(1 + rng.below(all.len()));
        all
    }

    #[test]
    fn the_grouped_pick_chooses_what_the_scan_of_every_group_chose_on_any_locus() {
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        let mut compared = 0usize;
        for case in 0..4_000 {
            let tvs = random_locus(&mut rng);
            let order = random_order(&mut rng);
            for scope in [PickScope::Allele, PickScope::AlleleGene] {
                let plan = PickPlan {
                    scope,
                    mode: PickMode::Reduce,
                    order: &order,
                };
                assert_eq!(
                    grouped_winners(&tvs, &plan),
                    reference_grouped_winners(&tvs, &plan),
                    "case {case}, {scope:?}, order {order:?}"
                );
                compared += 1;
            }
        }
        assert_eq!(compared, 8_000);
    }

    #[test]
    fn a_reducing_pick_keeps_what_the_flagging_pick_flags_on_any_locus() {
        // `retain_winners` finds a pair among the winners by search now, and the flagging mode
        // never searches at all: they have to agree about which pairs won
        let mut rng = Rng(0xD1B5_4A32_D192_ED03);
        for case in 0..2_000 {
            let tvs = random_locus(&mut rng);
            let order = random_order(&mut rng);
            for scope in [PickScope::Variant, PickScope::Allele, PickScope::AlleleGene] {
                let mut reduced = tvs.clone();
                apply_pick(
                    &mut reduced,
                    &PickPlan {
                        scope,
                        mode: PickMode::Reduce,
                        order: &order,
                    },
                );
                let mut flagged = tvs.clone();
                apply_pick(
                    &mut flagged,
                    &PickPlan {
                        scope,
                        mode: PickMode::Flag,
                        order: &order,
                    },
                );
                let kept = |tvs: &[TranscriptVariation], only_flagged: bool| {
                    tvs.iter()
                        .flat_map(|tv| {
                            tv.allele_annotations
                                .iter()
                                .filter(move |aa| !only_flagged || aa.pick)
                                .map(move |aa| {
                                    (tv.transcript_id.to_string(), aa.allele.to_string())
                                })
                        })
                        .collect::<Vec<_>>()
                };
                let winners = pick_winners(
                    &tvs,
                    &PickPlan {
                        scope,
                        mode: PickMode::Reduce,
                        order: &order,
                    },
                );
                if winners.is_empty() {
                    // nothing was a candidate, so a reducing pick leaves the record as it came
                    assert_eq!(
                        kept(&reduced, false),
                        kept(&tvs, false),
                        "case {case}, {scope:?}"
                    );
                } else {
                    // scaffold rows are never candidates, and a reducing pick leaves them where they were
                    let candidates_only = |v: Vec<(String, String)>| {
                        v.into_iter().filter(|(t, _)| t != "-").collect::<Vec<_>>()
                    };
                    assert_eq!(
                        candidates_only(kept(&reduced, false)),
                        candidates_only(kept(&flagged, true)),
                        "case {case}, {scope:?}, order {order:?}"
                    );
                }
            }
        }
    }
}
