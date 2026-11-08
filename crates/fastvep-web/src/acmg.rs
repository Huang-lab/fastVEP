//! Whether a genome can be classified, from the sources it holds.
//!
//! ACMG classification reads nine supplementary sources (`fastvep_classification::REFERENCE_STACK`).
//! A criterion that has no source reports `evaluated: false` and nothing else, so a server whose
//! human database directory holds ClinVar and dbSNP only answers `"acmg": true` with a
//! classification that is almost always VUS, and nothing in it says the classifier was not
//! given what it reads. A clinician cannot tell that from a variant on which the evidence was
//! weighed and found balanced, and a researcher reading the class has no reason to try. So the
//! server decides before it classifies:
//!
//! * a genome whose directory holds all nine sources is **ready**, and `acmg: true` is answered;
//! * any other genome is refused with HTTP 422 and the list of what is missing, unless the
//!   operator started the server with `--acmg-allow-partial`, which answers anyway and says in
//!   every response that the stack was incomplete.
//!
//! `GET /api/status` and `GET /api/genomes` carry the same verdict, so a client can offer
//! ACMG only where it will be answered and say why where it will not.

use fastvep_classification::{LoadedSources, SourceInfo, REFERENCE_STACK};
use serde::Serialize;

/// The verdict for one genome.
#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct AcmgReadiness {
    /// Every source of the reference stack is loaded: what the benchmark of the paper had, and
    /// the only state in which a classification is comparable with it.
    pub ready: bool,
    /// A request with `"acmg": true` for this genome will be answered. Differs from `ready` only
    /// on a server started with `--acmg-allow-partial`.
    pub available: bool,
    /// The sources that are loaded, by the `json_key` of their databases.
    pub loaded: Vec<&'static str>,
    /// The sources that are not, with the criteria that read each and what is lost without it.
    pub missing: Vec<&'static SourceInfo>,
}

impl AcmgReadiness {
    pub fn of(sources: LoadedSources, allow_partial: bool) -> Self {
        let missing: Vec<&'static SourceInfo> = sources.missing().collect();
        let ready = missing.is_empty();
        Self {
            ready,
            available: ready || allow_partial,
            loaded: sources.present().map(|s| s.id).collect(),
            missing,
        }
    }

    /// What a refused request is told, naming the genome and what is missing.
    pub fn refusal(&self, genome: &str) -> String {
        let missing: Vec<String> = self
            .missing
            .iter()
            .map(|s| format!("{} ({})", s.label, s.drives.join(", ")))
            .collect();
        format!(
            "ACMG classification is not available for genome '{genome}': {} of the {} sources it \
             reads are not loaded, and a classification made without them would rest on less \
             evidence than it appears to. Missing: {}.",
            self.missing.len(),
            REFERENCE_STACK.len(),
            missing.join("; ")
        )
    }

    /// The part of the verdict a successful ACMG response carries: whether the stack was
    /// complete, and which sources were missing if it was not.
    pub fn stack(&self) -> AcmgStack {
        AcmgStack {
            complete: self.ready,
            missing: self.missing.iter().map(|s| s.id).collect(),
        }
    }
}

/// What a response that carries classifications says about the sources behind them.
#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct AcmgStack {
    pub complete: bool,
    pub missing: Vec<&'static str>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_genome_with_the_whole_stack_is_ready_and_available() {
        let r = AcmgReadiness::of(LoadedSources::all(), false);
        assert!(r.ready && r.available);
        assert!(r.missing.is_empty());
        assert_eq!(r.loaded.len(), REFERENCE_STACK.len());
        assert_eq!(
            r.stack(),
            AcmgStack {
                complete: true,
                missing: vec![]
            }
        );
    }

    #[test]
    fn an_incomplete_stack_is_refused_unless_the_operator_allowed_it() {
        let loaded = LoadedSources::from_json_keys(["clinvar"]);
        let strict = AcmgReadiness::of(loaded, false);
        assert!(!strict.ready && !strict.available);
        assert_eq!(strict.loaded, ["clinvar"]);
        assert_eq!(strict.missing.len(), 8);

        let lenient = AcmgReadiness::of(loaded, true);
        assert!(!lenient.ready && lenient.available);
        assert!(!lenient.stack().complete);
        assert!(lenient.stack().missing.contains(&"gnomad"));
    }

    #[test]
    fn the_refusal_names_the_genome_the_count_and_each_missing_source() {
        let text = AcmgReadiness::of(LoadedSources::from_json_keys(["clinvar", "gnomad"]), false)
            .refusal("hg38_ensembl_115");
        assert!(text.contains("genome 'hg38_ensembl_115'"), "{text}");
        assert!(text.contains("7 of the 9 sources"), "{text}");
        assert!(text.contains("REVEL missense scores (PP3, BP4)"), "{text}");
        assert!(!text.contains("gnomAD allele frequencies ("), "{text}");
    }

    #[test]
    fn the_verdict_serialises_to_the_shape_the_api_documents() {
        let json = serde_json::to_value(AcmgReadiness::of(
            LoadedSources::from_json_keys(["gnomad"]),
            false,
        ))
        .unwrap();
        assert_eq!(json["ready"], false);
        assert_eq!(json["available"], false);
        assert_eq!(json["loaded"], serde_json::json!(["gnomad"]));
        let first = &json["missing"][0];
        assert_eq!(first["id"], "clinvar");
        assert_eq!(first["scope"], "allele");
        assert!(first["drives"].as_array().unwrap().len() > 1);
        assert!(first.get("source").is_none(), "{first}");
    }
}
