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
    let mut col_indices: Option<GnomadGeneCols> = None;
    let mut gene_order: Vec<String> = Vec::new();
    let mut gene_best: std::collections::HashMap<String, (u8, String)> =
        std::collections::HashMap::new();

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

        // gnomAD v4 emits one row per transcript per gene. Rank candidate rows
        // so the highest-ranked row wins regardless of input order:
        // canonical=true (2) > mane_select=true (1) > neither (dropped).
        // When neither column exists (v2.1: one row per gene), rank=0 and the
        // first-seen row per gene wins.
        let rank: u8 = if cols.canonical.is_some() || cols.mane_select.is_some() {
            let canonical_truthy = cols
                .canonical
                .map(|i| is_truthy(fields.get(i).copied().unwrap_or("")))
                .unwrap_or(false);
            let mane_truthy = cols
                .mane_select
                .map(|i| is_truthy(fields.get(i).copied().unwrap_or("")))
                .unwrap_or(false);
            if canonical_truthy {
                2
            } else if mane_truthy {
                1
            } else {
                continue;
            }
        } else {
            0
        };

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

        let json = format!("{{{}}}", parts.join(","));
        match gene_best.get(gene) {
            None => {
                gene_order.push(gene.to_string());
                gene_best.insert(gene.to_string(), (rank, json));
            }
            Some((existing_rank, _)) => {
                if rank > *existing_rank {
                    gene_best.insert(gene.to_string(), (rank, json));
                }
            }
        }
    }

    Ok(gene_order
        .into_iter()
        .map(|gene| {
            let (_, json) = gene_best.remove(&gene).expect("gene_best entry missing");
            GeneRecord {
                gene_symbol: gene,
                json,
            }
        })
        .collect())
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
                if let Some(i) = fields.iter().position(|f| f.eq_ignore_ascii_case(n)) {
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
        for i in [
            self.canonical,
            self.mane_select,
            self.pli,
            self.loeuf,
            self.mis_z,
            self.syn_z,
        ]
        .into_iter()
        .flatten()
        {
            m = m.max(i);
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

    #[test]
    fn test_parse_gnomad_gene_scores_v41_canonical_wins_regardless_of_row_order() {
        // Pre-fix the parser accepted the first row whose canonical OR
        // mane_select was truthy and locked the gene out for later rows. So a
        // mane_select-only row appearing before the canonical row would win,
        // and BRCA1 would carry the wrong transcript's metrics.
        let data = "\
gene\ttranscript\tcanonical\tmane_select\tlof.pLI\tlof.oe_ci.upper\tmis.z_score\tsyn.z_score
BRCA1\tENST_mane_only\tfalse\ttrue\t0.5\t0.99\t1.5\t0.1
BRCA1\tENST_canonical\ttrue\ttrue\t1.0\t0.03\t3.45\t0.12
BRCA1\tENST_alt\tfalse\tfalse\t0.2\t0.80\t0.5\t0.0
";
        let records = parse_gnomad_gene_scores(data.as_bytes()).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].gene_symbol, "BRCA1");
        assert!(
            records[0].json.contains("\"pLI\":1.0000"),
            "canonical row should win regardless of order, got {}",
            records[0].json
        );
        assert!(records[0].json.contains("\"loeuf\":0.0300"));
    }
}
