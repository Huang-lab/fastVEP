use anyhow::{Context, Result};
use oxivep_cache::fasta::FastaReader;
use oxivep_cache::gff::parse_gff3;
use oxivep_cache::providers::{FastaSequenceProvider, MemoryTranscriptProvider, SequenceProvider, TranscriptProvider};
use oxivep_consequence::ConsequencePredictor;
use oxivep_core::Consequence;
use oxivep_hgvs;
use oxivep_io::output;
use oxivep_io::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature};
use oxivep_io::vcf::VcfParser;
use std::fs::File;
use std::io::{self, BufWriter, Write};

pub struct AnnotateConfig {
    pub input: String,
    pub output: String,
    pub gff3: Option<String>,
    pub fasta: Option<String>,
    pub output_format: String,
    pub pick: bool,
    pub hgvs: bool,
    pub distance: u64,
}

pub fn run_annotate(config: AnnotateConfig) -> Result<()> {
    // Load transcript models from GFF3
    let transcripts = if let Some(ref gff3_path) = config.gff3 {
        let gff_file = File::open(gff3_path)
            .with_context(|| format!("Opening GFF3 file: {}", gff3_path))?;
        let trs = parse_gff3(gff_file)?;
        eprintln!("Loaded {} transcripts from {}", trs.len(), gff3_path);
        trs
    } else {
        eprintln!("Warning: No GFF3 file provided. Only intergenic variants will be annotated.");
        Vec::new()
    };

    let transcript_provider = MemoryTranscriptProvider::new(transcripts);

    // Load FASTA reference
    let seq_provider: Option<FastaSequenceProvider> = if let Some(ref fasta_path) = config.fasta {
        let fasta_file = File::open(fasta_path)
            .with_context(|| format!("Opening FASTA file: {}", fasta_path))?;
        let reader = FastaReader::from_reader(fasta_file)?;
        eprintln!("Loaded reference FASTA from {}", fasta_path);
        Some(FastaSequenceProvider::new(reader))
    } else {
        None
    };

    // Create consequence predictor
    let predictor = ConsequencePredictor::new(config.distance, config.distance);

    // Open input VCF
    let input_reader: Box<dyn io::Read> = if config.input == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(
            File::open(&config.input)
                .with_context(|| format!("Opening input file: {}", config.input))?,
        )
    };
    let mut vcf_parser = VcfParser::new(input_reader)?;

    // Open output
    let output_writer: Box<dyn io::Write> = if config.output == "-" {
        Box::new(io::stdout())
    } else {
        Box::new(
            File::create(&config.output)
                .with_context(|| format!("Creating output file: {}", config.output))?,
        )
    };
    let mut writer = BufWriter::new(output_writer);

    // Write headers based on output format
    match config.output_format.as_str() {
        "vcf" => {
            // Pass through original VCF headers
            for header_line in vcf_parser.header_lines() {
                if header_line.starts_with("#CHROM") {
                    // Insert CSQ header before #CHROM
                    writeln!(writer, "{}", output::csq_header_line(output::DEFAULT_CSQ_FIELDS))?;
                }
                writeln!(writer, "{}", header_line)?;
            }
        }
        "tab" => {
            writeln!(
                writer,
                "## OxiVEP output"
            )?;
            writeln!(
                writer,
                "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation"
            )?;
        }
        "json" => {
            writeln!(writer, "[")?;
        }
        _ => {}
    }

    // Process variants
    let mut count = 0u64;
    let mut first_json = true;

    while let Some(mut vf) = vcf_parser.next_variant()? {
        // Find overlapping transcripts
        let chrom = &vf.position.chromosome;
        let query_start = if vf.position.start > config.distance {
            vf.position.start - config.distance
        } else {
            1
        };
        let query_end = vf.position.end + config.distance;
        let overlapping = transcript_provider.get_transcripts(chrom, query_start, query_end)?;

        if overlapping.is_empty() {
            // Intergenic
            annotate_intergenic(&mut vf);
        } else {
            // Get reference sequence if available
            let ref_seq = seq_provider.as_ref().and_then(|sp| {
                sp.fetch_sequence(chrom, query_start, query_end).ok()
            });

            // Run consequence prediction
            let result = predictor.predict(
                &vf.position,
                &vf.ref_allele,
                &vf.alt_alleles,
                &overlapping,
                ref_seq.as_deref(),
            );

            // Convert prediction results to VariationFeature annotations
            for tc in &result.transcript_consequences {
                let transcript = overlapping.iter().find(|t| t.stable_id == tc.transcript_id);

                let allele_annotations: Vec<AlleleAnnotation> = tc
                    .allele_consequences
                    .iter()
                    .map(|ac| {
                        let mut ann = AlleleAnnotation {
                            allele: ac.allele.clone(),
                            consequences: ac.consequences.clone(),
                            impact: ac.impact,
                            cdna_position: zip_positions(ac.cdna_start, ac.cdna_end),
                            cds_position: zip_positions(ac.cds_start, ac.cds_end),
                            protein_position: zip_positions(ac.protein_start, ac.protein_end),
                            amino_acids: ac.amino_acids.clone(),
                            codons: ac.codons.clone(),
                            exon: ac.exon,
                            intron: ac.intron,
                            distance: ac.distance,
                            hgvsc: None,
                            hgvsp: None,
                            hgvsg: None,
                            existing_variation: vec![],
                            sift: None,
                            polyphen: None,
                        };

                        // Generate HGVS if requested
                        if config.hgvs {
                            ann.hgvsg = Some(oxivep_hgvs::hgvsg(
                                chrom,
                                vf.position.start,
                                vf.position.end,
                                &vf.ref_allele,
                                &ac.allele,
                            ));

                            if let (Some(cs), Some(ce)) = (ac.cdna_start, ac.cdna_end) {
                                if let Some(tr) = transcript {
                                    if let Some(coding_start) = tr.cdna_coding_start {
                                        ann.hgvsc = oxivep_hgvs::hgvsc(
                                            &tc.transcript_id,
                                            cs, ce,
                                            &vf.ref_allele,
                                            &ac.allele,
                                            coding_start,
                                        );
                                    }
                                }
                            }

                            if let (Some(ref aa), Some(ps)) = (&ac.amino_acids, ac.protein_start) {
                                if let Some(tr) = transcript {
                                    if let Some(ref pid) = tr.protein_id {
                                        let is_fs = ac.consequences.contains(&Consequence::FrameshiftVariant);
                                        let ref_aa_byte = aa.0.as_bytes().first().copied().unwrap_or(b'X');
                                        let alt_aa_byte = aa.1.as_bytes().first().copied().unwrap_or(b'X');
                                        ann.hgvsp = oxivep_hgvs::hgvsp(
                                            pid, ps, ref_aa_byte, alt_aa_byte, is_fs,
                                        );
                                    }
                                }
                            }
                        }

                        ann
                    })
                    .collect();

                // Apply pick filter if needed
                let should_include = !config.pick || tc.canonical || vf.transcript_variations.is_empty();

                if should_include {
                    vf.transcript_variations.push(TranscriptVariation {
                        transcript_id: tc.transcript_id.clone(),
                        gene_id: tc.gene_id.clone(),
                        gene_symbol: tc.gene_symbol.clone(),
                        biotype: tc.biotype.clone(),
                        allele_annotations,
                        canonical: tc.canonical,
                        strand: tc.strand,
                        source: None,
                        protein_id: transcript.and_then(|t| t.protein_id.clone()),
                        mane_select: transcript.and_then(|t| t.mane_select.clone()),
                        mane_plus_clinical: transcript.and_then(|t| t.mane_plus_clinical.clone()),
                        tsl: transcript.and_then(|t| t.tsl),
                        appris: transcript.and_then(|t| t.appris.clone()),
                        ccds: transcript.and_then(|t| t.ccds.clone()),
                    });
                }
            }

            vf.compute_most_severe();
        }

        // Write output
        match config.output_format.as_str() {
            "vcf" => write_vcf_line(&mut writer, &vf)?,
            "tab" => {
                for line in output::format_tab_line(&vf) {
                    writeln!(writer, "{}", line)?;
                }
            }
            "json" => {
                if !first_json {
                    writeln!(writer, ",")?;
                }
                first_json = false;
                let json = output::format_json(&vf);
                write!(writer, "{}", serde_json::to_string_pretty(&json)?)?;
            }
            _ => {}
        }

        count += 1;
    }

    // Close JSON array
    if config.output_format == "json" {
        writeln!(writer, "\n]")?;
    }

    writer.flush()?;
    eprintln!("Annotated {} variants", count);

    Ok(())
}

fn annotate_intergenic(vf: &mut VariationFeature) {
    for alt in &vf.alt_alleles {
        vf.transcript_variations.push(TranscriptVariation {
            transcript_id: "-".into(),
            gene_id: "-".into(),
            gene_symbol: None,
            biotype: "-".into(),
            allele_annotations: vec![AlleleAnnotation {
                allele: alt.clone(),
                consequences: vec![Consequence::IntergenicVariant],
                impact: oxivep_core::Impact::Modifier,
                cdna_position: None,
                cds_position: None,
                protein_position: None,
                amino_acids: None,
                codons: None,
                exon: None,
                intron: None,
                distance: None,
                hgvsc: None,
                hgvsp: None,
                hgvsg: None,
                existing_variation: vec![],
                sift: None,
                polyphen: None,
            }],
            canonical: false,
            strand: oxivep_core::Strand::Forward,
            source: None,
            protein_id: None,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
        });
    }
    vf.most_severe_consequence = Some(Consequence::IntergenicVariant);
}

fn write_vcf_line(writer: &mut impl Write, vf: &VariationFeature) -> Result<()> {
    if let Some(ref fields) = vf.vcf_fields {
        let csq = output::format_csq(vf, output::DEFAULT_CSQ_FIELDS);
        let info = if fields.info == "." && !csq.is_empty() {
            format!("CSQ={}", csq)
        } else if !csq.is_empty() {
            format!("{};CSQ={}", fields.info, csq)
        } else {
            fields.info.clone()
        };

        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            fields.chrom, fields.pos, fields.id, fields.ref_allele, fields.alt,
            fields.qual, fields.filter, info
        )?;

        for rest_field in &fields.rest {
            write!(writer, "\t{}", rest_field)?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

fn zip_positions(start: Option<u64>, end: Option<u64>) -> Option<(u64, u64)> {
    match (start, end) {
        (Some(s), Some(e)) => Some((s, e)),
        (Some(s), None) => Some((s, s)),
        (None, Some(e)) => Some((e, e)),
        _ => None,
    }
}

use serde_json;
