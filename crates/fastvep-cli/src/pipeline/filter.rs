//! `fastvep filter`: filter an annotated VCF by CSQ expression.

use anyhow::{Context, Result};
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write};

/// Filter annotated VCF by CSQ field expressions.
///
/// Reads a VCF with CSQ INFO annotations (produced by `fastvep annotate --output-format vcf`),
/// parses each CSQ entry into fields, evaluates the filter expression, and writes matching lines.
pub fn run_filter(input: &str, output_path: &str, filter_expr: &str) -> Result<()> {
    use fastvep_filter::{Filter, FilterContext};

    let filter = Filter::parse(filter_expr)
        .with_context(|| format!("Parsing filter expression: {}", filter_expr))?;

    // Open input
    let reader: Box<dyn BufRead> = if input == "-" {
        Box::new(io::BufReader::new(io::stdin()))
    } else {
        let f = File::open(input).with_context(|| format!("Opening input: {}", input))?;
        Box::new(io::BufReader::new(f))
    };

    // Open output. 1 MiB buffer same rationale as the main annotation path.
    let mut writer: Box<dyn Write> = if output_path == "-" {
        Box::new(BufWriter::with_capacity(1 << 20, io::stdout()))
    } else {
        let f = File::create(output_path)
            .with_context(|| format!("Creating output: {}", output_path))?;
        Box::new(BufWriter::with_capacity(1 << 20, f))
    };

    // Parse CSQ header to get field names
    let mut csq_fields: Vec<String> = Vec::new();
    let mut kept = 0u64;
    let mut total = 0u64;
    let mut malformed = 0u64;

    for line in reader.lines() {
        let line = line?;

        // Header lines: pass through, extract CSQ format
        if line.starts_with('#') {
            // Parse CSQ format from ##INFO=<ID=CSQ,...,Description="... Format: A|B|C">
            if line.starts_with("##INFO=<ID=CSQ") {
                if let Some(fmt_start) = line.find("Format: ") {
                    let fmt = &line[fmt_start + 8..];
                    let fmt = fmt.trim_end_matches('"').trim_end_matches('>');
                    csq_fields = fmt.split('|').map(|s| s.to_string()).collect();
                    eprintln!(
                        "[filter] CSQ format: {} fields ({})",
                        csq_fields.len(),
                        csq_fields
                            .iter()
                            .take(5)
                            .cloned()
                            .collect::<Vec<_>>()
                            .join(", ")
                    );
                }
            }
            writeln!(writer, "{}", line)?;
            continue;
        }

        total += 1;

        // Parse VCF data line
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            malformed += 1;
            eprintln!(
                "[filter] skipping malformed line {} (fewer than 8 tab-separated fields): {}",
                total, line
            );
            continue;
        }

        let info = fields[7];

        // Extract CSQ entries from INFO field
        let csq_str = info
            .split(';')
            .find(|f| f.starts_with("CSQ="))
            .map(|f| &f[4..]);

        let Some(csq_str) = csq_str else {
            // No CSQ field — skip this line (doesn't match any filter)
            continue;
        };

        // Check if ANY CSQ entry matches the filter
        let mut any_match = false;
        for entry in csq_str.split(',') {
            let values: Vec<&str> = entry.split('|').collect();
            let mut ctx = FilterContext::new();

            for (i, val) in values.iter().enumerate() {
                if let Some(field_name) = csq_fields.get(i) {
                    if !val.is_empty() {
                        ctx.set(field_name, val);
                    }
                }
            }

            // Also set VCF-level fields
            ctx.set("CHROM", fields[0]);
            ctx.set("POS", fields[1]);
            ctx.set("REF", fields[3]);
            ctx.set("ALT", fields[4]);
            ctx.set("FILTER", fields[6]);

            if filter.matches(&ctx) {
                any_match = true;
                break;
            }
        }

        if any_match {
            writeln!(writer, "{}", line)?;
            kept += 1;
        }
    }

    eprintln!(
        "[filter] {} of {} variants passed filter: {}",
        kept, total, filter_expr
    );
    if malformed > 0 {
        eprintln!(
            "[filter] {} of {} lines were skipped as malformed (fewer than 8 tab-separated fields)",
            malformed, total
        );
    }

    Ok(())
}
