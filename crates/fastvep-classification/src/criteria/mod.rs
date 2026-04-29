pub mod pvs1;
pub mod pathogenic_strong;
pub mod pathogenic_moderate;
pub mod pathogenic_supporting;
pub mod benign_standalone;
pub mod benign_strong;
pub mod benign_supporting;

use crate::config::AcmgConfig;
use crate::sa_extract::ClassificationInput;
use crate::types::EvidenceCriterion;

/// Evaluate all 28 ACMG-AMP criteria and return the full list.
pub fn evaluate_all_criteria(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> Vec<EvidenceCriterion> {
    let gene = input.gene_symbol.as_deref();

    let mut criteria = Vec::with_capacity(28);

    // Pathogenic Very Strong
    let pvs1 = pvs1::evaluate_pvs1(input, config);
    if !is_disabled(gene, &pvs1.code, config) {
        criteria.push(pvs1);
    }

    // Pathogenic Strong
    for c in pathogenic_strong::evaluate_all(input, config) {
        if !is_disabled(gene, &c.code, config) {
            criteria.push(c);
        }
    }

    // Pathogenic Moderate
    for c in pathogenic_moderate::evaluate_all(input, config) {
        if !is_disabled(gene, &c.code, config) {
            criteria.push(c);
        }
    }

    // Pathogenic Supporting
    for c in pathogenic_supporting::evaluate_all(input, config) {
        if !is_disabled(gene, &c.code, config) {
            criteria.push(c);
        }
    }

    // Benign Standalone
    let ba1 = benign_standalone::evaluate_ba1(input, config);
    if !is_disabled(gene, &ba1.code, config) {
        criteria.push(ba1);
    }

    // Benign Strong
    for c in benign_strong::evaluate_all(input, config) {
        if !is_disabled(gene, &c.code, config) {
            criteria.push(c);
        }
    }

    // Benign Supporting
    for c in benign_supporting::evaluate_all(input, config) {
        if !is_disabled(gene, &c.code, config) {
            criteria.push(c);
        }
    }

    criteria
}

fn is_disabled(gene: Option<&str>, code: &str, config: &AcmgConfig) -> bool {
    gene.map_or(false, |g| config.is_criterion_disabled(g, code))
}
