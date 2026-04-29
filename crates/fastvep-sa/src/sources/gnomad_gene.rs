//! gnomAD gene constraint scores parser for building .oga files.
//!
//! Extracts pLI, LOEUF, mis_z, and syn_z per gene.

use crate::common::GeneRecord;
use anyhow::{Context, Result};
use std::io::BufRead;

/// Parse gnomAD constraint metrics TSV into GeneRecords.
///
/// Supports both gnomAD v2.1 column naming (`pLI`, `oe_lof_upper`, `mis_z`,
/// `syn_z`) and the v4.x dotted-namespace naming (`lof.pLI`, `lof.oe_ci.upper`,
/// `mis.z_score`, `syn.z_score`). For v4 we also collapse to a single
/// canonical-transcript row per gene, since v4 emits one row per transcript.
pub fn parse_gnomad_gene_scores<R: BufRead>(reader: R) -> Result<Vec<GeneRecord>> {
    let mut records = Vec::new();
    let mut seen_genes: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut col_indices: Option<GnomadGeneCols> = None;

    for line in reader.lines() {
        let line = line.context("Reading gnomAD gene scores")?;
        if line.starts_with("gene\t") || line.starts_with("#") {
            col_indices = Some(GnomadGeneCols::from_header(&line));
            continue;
        }
        if line.is_empty() {
            continue;
        }

        let cols = match &col_indices {
            Some(c) => c,
            None => continue,
        };

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= cols.max_idx() {
            continue;
        }

        let gene = fields[cols.gene].trim();
        if gene.is_empty() {
            continue;
        }

        // gnomAD v4 emits one row per transcript per gene. Prefer the
        // canonical / MANE_select row when those columns exist; otherwise
        // accept the first row we see for each gene.
        if let Some(idx) = cols.canonical {
            if !is_truthy(fields.get(idx).copied().unwrap_or("")) {
                if let Some(mane_idx) = cols.mane_select {
                    if !is_truthy(fields.get(mane_idx).copied().unwrap_or("")) {
                        continue;
                    }
                } else {
                    continue;
                }
            }
        }
        if !seen_genes.insert(gene.to_string()) {
            continue;
        }

        let mut parts = Vec::new();

        if let Some(idx) = cols.pli {
            if let Ok(v) = fields[idx].parse::<f64>() {
                parts.push(format!("\"pLI\":{:.4}", v));
            }
        }
        if let Some(idx) = cols.loeuf {
            if let Ok(v) = fields[idx].parse::<f64>() {
                parts.push(format!("\"loeuf\":{:.4}", v));
            }
        }
        if let Some(idx) = cols.mis_z {
            if let Ok(v) = fields[idx].parse::<f64>() {
                parts.push(format!("\"misZ\":{:.2}", v));
            }
        }
        if let Some(idx) = cols.syn_z {
            if let Ok(v) = fields[idx].parse::<f64>() {
                parts.push(format!("\"synZ\":{:.2}", v));
            }
        }

        if parts.is_empty() {
            continue;
        }

        records.push(GeneRecord {
            gene_symbol: gene.to_string(),
            json: format!("{{{}}}", parts.join(",")),
        });
    }

    Ok(records)
}

fn is_truthy(s: &str) -> bool {
    matches!(s.trim().to_lowercase().as_str(), "true" | "t" | "1" | "yes")
}

struct GnomadGeneCols {
    gene: usize,
    canonical: Option<usize>,
    mane_select: Option<usize>,
    pli: Option<usize>,
    loeuf: Option<usize>,
    mis_z: Option<usize>,
    syn_z: Option<usize>,
}

impl GnomadGeneCols {
    fn from_header(header: &str) -> Self {
        let fields: Vec<&str> = header.split('\t').collect();
        let find = |needles: &[&str]| {
            for n in needles {
                if let Some(i) = fields
                    .iter()
                    .position(|f| f.eq_ignore_ascii_case(n))
                {
                    return Some(i);
                }
            }
            None
        };

        Self {
            gene: find(&["gene", "gene_symbol"]).unwrap_or(0),
            canonical: find(&["canonical"]),
            mane_select: find(&["mane_select"]),
            // v4: lof.pLI / lof_hc_lc.pLI ; v2.1: pLI
            pli: find(&["lof.pLI", "lof_hc_lc.pLI", "pLI"]),
            // v4: lof.oe_ci.upper ; v2.1: oe_lof_upper / loeuf
            loeuf: find(&["lof.oe_ci.upper", "oe_lof_upper", "loeuf"]),
            // v4: mis.z_score ; v2.1: mis_z
            mis_z: find(&["mis.z_score", "mis_z"]),
            // v4: syn.z_score ; v2.1: syn_z
            syn_z: find(&["syn.z_score", "syn_z"]),
        }
    }

    fn max_idx(&self) -> usize {
        let mut m = self.gene;
        for opt in [
            self.canonical,
            self.mane_select,
            self.pli,
            self.loeuf,
            self.mis_z,
            self.syn_z,
        ] {
            if let Some(i) = opt {
                m = m.max(i);
            }
        }
        m
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_gnomad_gene_scores_v21_format() {
        let data = "\
gene\ttranscript\tobs_lof\texp_lof\toe_lof\toe_lof_upper\tpLI\tmis_z\tsyn_z
BRCA1\tENST00000357654\t0\t50.2\t0.00\t0.03\t1.0000\t3.45\t0.12
TP53\tENST00000269305\t0\t25.1\t0.00\t0.05\t0.9999\t5.67\t-0.34
";
        let records = parse_gnomad_gene_scores(data.as_bytes()).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].gene_symbol, "BRCA1");
        assert!(records[0].json.contains("\"pLI\":1.0000"));
        assert!(records[0].json.contains("\"loeuf\":0.0300"));
        assert!(records[0].json.contains("\"misZ\":3.45"));

        assert_eq!(records[1].gene_symbol, "TP53");
        assert!(records[1].json.contains("\"pLI\":0.9999"));
    }

    #[test]
    fn test_parse_gnomad_gene_scores_v41_format() {
        // gnomAD v4.1 emits one row per transcript with dotted-namespace
        // column names. The parser should resolve `lof.pLI`,
        // `lof.oe_ci.upper`, `mis.z_score`, `syn.z_score` and prefer the
        // canonical row when there's both canonical and non-canonical.
        let data = "\
gene\ttranscript\tcanonical\tmane_select\tlof.pLI\tlof.oe_ci.upper\tmis.z_score\tsyn.z_score
BRCA1\tENST_alt\tfalse\tfalse\t0.5\t0.99\t1.5\t0.1
BRCA1\tENST_canonical\ttrue\ttrue\t1.0\t0.03\t3.45\t0.12
TP53\tENST_canonical\ttrue\ttrue\t0.9999\t0.05\t5.67\t-0.34
";
        let records = parse_gnomad_gene_scores(data.as_bytes()).unwrap();
        assert_eq!(
            records.len(),
            2,
            "expected one row per gene (canonical only); got {:?}",
            records
        );
        assert_eq!(records[0].gene_symbol, "BRCA1");
        assert!(
            records[0].json.contains("\"pLI\":1.0000"),
            "expected canonical pLI=1.0, got {}",
            records[0].json
        );
        assert!(records[0].json.contains("\"loeuf\":0.0300"));
        assert!(records[0].json.contains("\"misZ\":3.45"));
        assert_eq!(records[1].gene_symbol, "TP53");
    }
}
