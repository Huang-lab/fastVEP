use anyhow::Result;
use oxivep_core::Strand;
use oxivep_genome::{Exon, Gene, Transcript, Translation};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};

/// Parse a GFF3 file into Transcript models.
///
/// Builds gene -> transcript -> exon/CDS hierarchy from GFF3 features.
pub fn parse_gff3<R: Read>(reader: R) -> Result<Vec<Transcript>> {
    let buf = BufReader::new(reader);
    let mut genes: HashMap<String, GffGene> = HashMap::new();
    let mut transcripts: HashMap<String, GffTranscript> = HashMap::new();
    let mut exons: Vec<GffExon> = Vec::new();
    let mut cds_features: Vec<GffCds> = Vec::new();

    for line in buf.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let seqid = fields[0];
        let feature_type = fields[2];
        let start: u64 = fields[3].parse().unwrap_or(0);
        let end: u64 = fields[4].parse().unwrap_or(0);
        let strand = match fields[6] {
            "-" => Strand::Reverse,
            _ => Strand::Forward,
        };
        let phase: i8 = fields[7].parse().unwrap_or(-1);
        let attrs = parse_attributes(fields[8]);

        match feature_type {
            "gene" | "pseudogene" => {
                let gene_id = attrs
                    .get("ID")
                    .or_else(|| attrs.get("gene_id"))
                    .cloned()
                    .unwrap_or_default();
                // Strip "gene:" prefix if present (Ensembl GFF3)
                let gene_id = gene_id.strip_prefix("gene:").unwrap_or(&gene_id).to_string();
                let symbol = attrs.get("Name").cloned();
                let biotype = attrs
                    .get("biotype")
                    .or_else(|| attrs.get("gene_biotype"))
                    .cloned()
                    .unwrap_or_else(|| "unknown".to_string());

                genes.insert(
                    gene_id.clone(),
                    GffGene {
                        id: gene_id,
                        symbol,
                        biotype,
                        chromosome: seqid.to_string(),
                        start,
                        end,
                        strand,
                    },
                );
            }
            "mRNA" | "transcript" | "lnc_RNA" | "miRNA" | "snRNA" | "snoRNA" | "rRNA"
            | "ncRNA" | "tRNA" | "scRNA" | "V_gene_segment" | "D_gene_segment"
            | "J_gene_segment" | "C_gene_segment" | "NMD_transcript_variant"
            | "pseudogenic_transcript" => {
                let transcript_id = attrs
                    .get("ID")
                    .or_else(|| attrs.get("transcript_id"))
                    .cloned()
                    .unwrap_or_default();
                let transcript_id = transcript_id
                    .strip_prefix("transcript:")
                    .unwrap_or(&transcript_id)
                    .to_string();
                let parent = attrs
                    .get("Parent")
                    .cloned()
                    .unwrap_or_default();
                let parent = parent.strip_prefix("gene:").unwrap_or(&parent).to_string();
                let biotype = attrs
                    .get("biotype")
                    .or_else(|| attrs.get("transcript_biotype"))
                    .cloned()
                    .unwrap_or_else(|| feature_type.to_string());
                let canonical = attrs.get("tag").map(|t| t.contains("Ensembl_canonical")).unwrap_or(false);

                transcripts.insert(
                    transcript_id.clone(),
                    GffTranscript {
                        id: transcript_id,
                        parent_gene: parent,
                        biotype,
                        chromosome: seqid.to_string(),
                        start,
                        end,
                        strand,
                        canonical,
                    },
                );
            }
            "exon" => {
                let parent = attrs.get("Parent").cloned().unwrap_or_default();
                let parent = parent
                    .strip_prefix("transcript:")
                    .unwrap_or(&parent)
                    .to_string();
                let exon_id = attrs
                    .get("ID")
                    .or_else(|| attrs.get("exon_id"))
                    .cloned()
                    .unwrap_or_default();
                let exon_id = exon_id.strip_prefix("exon:").unwrap_or(&exon_id).to_string();
                let rank: u32 = attrs
                    .get("rank")
                    .or_else(|| attrs.get("exon_number"))
                    .and_then(|r| r.parse().ok())
                    .unwrap_or(0);

                exons.push(GffExon {
                    id: exon_id,
                    parent_transcript: parent,
                    start,
                    end,
                    strand,
                    phase,
                    rank,
                });
            }
            "CDS" => {
                let parent = attrs.get("Parent").cloned().unwrap_or_default();
                let parent = parent
                    .strip_prefix("transcript:")
                    .unwrap_or(&parent)
                    .to_string();
                let protein_id = attrs
                    .get("ID")
                    .or_else(|| attrs.get("protein_id"))
                    .cloned()
                    .unwrap_or_default();
                let protein_id = protein_id
                    .strip_prefix("CDS:")
                    .unwrap_or(&protein_id)
                    .to_string();

                cds_features.push(GffCds {
                    parent_transcript: parent,
                    protein_id,
                    start,
                    end,
                    strand,
                    phase,
                });
            }
            _ => {}
        }
    }

    // Build transcripts
    let mut result = Vec::new();

    for (tid, gff_tr) in &transcripts {
        let gene = genes.get(&gff_tr.parent_gene);
        let gene_model = gene
            .map(|g| Gene {
                stable_id: g.id.clone(),
                symbol: g.symbol.clone(),
                symbol_source: Some("GFF3".to_string()),
                hgnc_id: None,
                biotype: g.biotype.clone(),
                chromosome: g.chromosome.clone(),
                start: g.start,
                end: g.end,
                strand: g.strand,
            })
            .unwrap_or_else(|| Gene {
                stable_id: gff_tr.parent_gene.clone(),
                symbol: None,
                symbol_source: None,
                hgnc_id: None,
                biotype: "unknown".into(),
                chromosome: gff_tr.chromosome.clone(),
                start: gff_tr.start,
                end: gff_tr.end,
                strand: gff_tr.strand,
            });

        // Collect exons for this transcript
        let mut tr_exons: Vec<Exon> = exons
            .iter()
            .filter(|e| e.parent_transcript == *tid)
            .map(|e| Exon {
                stable_id: e.id.clone(),
                start: e.start,
                end: e.end,
                strand: e.strand,
                phase: e.phase,
                end_phase: -1,
                rank: e.rank,
            })
            .collect();

        // Sort exons by position
        match gff_tr.strand {
            Strand::Forward => tr_exons.sort_by_key(|e| e.start),
            Strand::Reverse => tr_exons.sort_by(|a, b| b.start.cmp(&a.start)),
        }

        // Assign ranks if not set
        for (i, exon) in tr_exons.iter_mut().enumerate() {
            if exon.rank == 0 {
                exon.rank = (i + 1) as u32;
            }
        }

        // Collect CDS features for this transcript
        let tr_cds: Vec<&GffCds> = cds_features
            .iter()
            .filter(|c| c.parent_transcript == *tid)
            .collect();

        let translation = if !tr_cds.is_empty() {
            let protein_id = tr_cds[0].protein_id.clone();
            let cds_start = tr_cds.iter().map(|c| c.start).min().unwrap();
            let cds_end = tr_cds.iter().map(|c| c.end).max().unwrap();

            // genomic_start/end always refer to the min/max genomic coords
            let genomic_start = cds_start;
            let genomic_end = cds_end;

            // For translation: "start" means start of translation in transcript order
            // Forward: translation starts at cds_start, ends at cds_end
            // Reverse: translation starts at cds_end, ends at cds_start
            let (tl_start_pos, tl_end_pos) = match gff_tr.strand {
                Strand::Forward => (cds_start, cds_end),
                Strand::Reverse => (cds_end, cds_start),
            };

            let start_exon_rank = tr_exons
                .iter()
                .find(|e| tl_start_pos >= e.start && tl_start_pos <= e.end)
                .map(|e| e.rank)
                .unwrap_or(1);
            let end_exon_rank = tr_exons
                .iter()
                .find(|e| tl_end_pos >= e.start && tl_end_pos <= e.end)
                .map(|e| e.rank)
                .unwrap_or(1);

            let start_exon = tr_exons.iter().find(|e| e.rank == start_exon_rank);
            let end_exon = tr_exons.iter().find(|e| e.rank == end_exon_rank);

            let start_offset = start_exon
                .map(|e| match gff_tr.strand {
                    Strand::Forward => tl_start_pos.saturating_sub(e.start),
                    Strand::Reverse => e.end.saturating_sub(tl_start_pos),
                })
                .unwrap_or(0);
            let end_offset = end_exon
                .map(|e| match gff_tr.strand {
                    Strand::Forward => tl_end_pos.saturating_sub(e.start),
                    Strand::Reverse => e.end.saturating_sub(tl_end_pos),
                })
                .unwrap_or(0);

            Some(Translation {
                stable_id: protein_id,
                genomic_start,
                genomic_end,
                start_exon_rank,
                start_exon_offset: start_offset,
                end_exon_rank,
                end_exon_offset: end_offset,
            })
        } else {
            None
        };

        // Compute cDNA coding positions
        let (cdna_coding_start, cdna_coding_end) = if let Some(ref tl) = translation {
            let mut cdna_pos = 0u64;
            let mut cs = None;
            let mut ce = None;

            let sorted_exons: Vec<&Exon> = match gff_tr.strand {
                Strand::Forward => {
                    let mut e: Vec<&Exon> = tr_exons.iter().collect();
                    e.sort_by_key(|e| e.start);
                    e
                }
                Strand::Reverse => {
                    let mut e: Vec<&Exon> = tr_exons.iter().collect();
                    e.sort_by(|a, b| b.start.cmp(&a.start));
                    e
                }
            };

            for exon in &sorted_exons {
                let exon_len = exon.end - exon.start + 1;
                match gff_tr.strand {
                    Strand::Forward => {
                        if cs.is_none() && tl.genomic_start >= exon.start && tl.genomic_start <= exon.end {
                            cs = Some(cdna_pos + (tl.genomic_start - exon.start) + 1);
                        }
                        if ce.is_none() && tl.genomic_end >= exon.start && tl.genomic_end <= exon.end {
                            ce = Some(cdna_pos + (tl.genomic_end - exon.start) + 1);
                        }
                    }
                    Strand::Reverse => {
                        if cs.is_none() && tl.genomic_end >= exon.start && tl.genomic_end <= exon.end {
                            cs = Some(cdna_pos + (exon.end - tl.genomic_end) + 1);
                        }
                        if ce.is_none() && tl.genomic_start >= exon.start && tl.genomic_start <= exon.end {
                            ce = Some(cdna_pos + (exon.end - tl.genomic_start) + 1);
                        }
                    }
                }
                cdna_pos += exon_len;
            }
            (cs, ce)
        } else {
            (None, None)
        };

        result.push(Transcript {
            stable_id: tid.clone(),
            gene: gene_model,
            biotype: gff_tr.biotype.clone(),
            chromosome: gff_tr.chromosome.clone(),
            start: gff_tr.start,
            end: gff_tr.end,
            strand: gff_tr.strand,
            exons: tr_exons,
            translation,
            cdna_coding_start,
            cdna_coding_end,
            coding_region_start: cds_features
                .iter()
                .filter(|c| c.parent_transcript == *tid)
                .map(|c| c.start)
                .min(),
            coding_region_end: cds_features
                .iter()
                .filter(|c| c.parent_transcript == *tid)
                .map(|c| c.end)
                .max(),
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: gff_tr.canonical,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: cds_features
                .iter()
                .find(|c| c.parent_transcript == *tid)
                .map(|c| c.protein_id.clone()),
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: Some("GFF3".into()),
            gencode_primary: false,
        });
    }

    Ok(result)
}

fn parse_attributes(attr_str: &str) -> HashMap<String, String> {
    let mut attrs = HashMap::new();
    for part in attr_str.split(';') {
        let part = part.trim();
        if let Some((key, value)) = part.split_once('=') {
            let value = url_decode(value);
            attrs.insert(key.to_string(), value);
        }
    }
    attrs
}

fn url_decode(s: &str) -> String {
    let mut result = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if c == '%' {
            let hex: String = chars.by_ref().take(2).collect();
            if let Ok(byte) = u8::from_str_radix(&hex, 16) {
                result.push(byte as char);
            } else {
                result.push('%');
                result.push_str(&hex);
            }
        } else {
            result.push(c);
        }
    }
    result
}

#[derive(Debug)]
struct GffGene {
    id: String,
    symbol: Option<String>,
    biotype: String,
    chromosome: String,
    start: u64,
    end: u64,
    strand: Strand,
}

#[derive(Debug)]
#[allow(dead_code)]
struct GffTranscript {
    id: String,
    parent_gene: String,
    biotype: String,
    chromosome: String,
    start: u64,
    end: u64,
    strand: Strand,
    canonical: bool,
}

#[derive(Debug)]
struct GffExon {
    id: String,
    parent_transcript: String,
    start: u64,
    end: u64,
    strand: Strand,
    phase: i8,
    rank: u32,
}

#[derive(Debug)]
#[allow(dead_code)]
struct GffCds {
    parent_transcript: String,
    protein_id: String,
    start: u64,
    end: u64,
    strand: Strand,
    phase: i8,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_gff3() -> &'static str {
        "##gff-version 3
chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene:ENSG00000001;Name=TESTGENE;biotype=protein_coding
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST00000001;Parent=gene:ENSG00000001;biotype=protein_coding;tag=Ensembl_canonical
chr1\tensembl\texon\t1000\t1200\t.\t+\t.\tID=exon:ENSE00000001;Parent=transcript:ENST00000001;rank=1
chr1\tensembl\texon\t2000\t2300\t.\t+\t.\tID=exon:ENSE00000002;Parent=transcript:ENST00000001;rank=2
chr1\tensembl\texon\t4000\t5000\t.\t+\t.\tID=exon:ENSE00000003;Parent=transcript:ENST00000001;rank=3
chr1\tensembl\tCDS\t1050\t1200\t.\t+\t0\tID=CDS:ENSP00000001;Parent=transcript:ENST00000001
chr1\tensembl\tCDS\t2000\t2300\t.\t+\t0\tID=CDS:ENSP00000001;Parent=transcript:ENST00000001
chr1\tensembl\tCDS\t4000\t4500\t.\t+\t1\tID=CDS:ENSP00000001;Parent=transcript:ENST00000001"
    }

    #[test]
    fn test_parse_gff3_basic() {
        let transcripts = parse_gff3(sample_gff3().as_bytes()).unwrap();
        assert_eq!(transcripts.len(), 1);

        let tr = &transcripts[0];
        assert_eq!(tr.stable_id, "ENST00000001");
        assert_eq!(tr.gene.stable_id, "ENSG00000001");
        assert_eq!(tr.gene.symbol.as_deref(), Some("TESTGENE"));
        assert_eq!(tr.biotype, "protein_coding");
        assert_eq!(tr.chromosome, "chr1");
        assert_eq!(tr.start, 1000);
        assert_eq!(tr.end, 5000);
        assert_eq!(tr.strand, Strand::Forward);
        assert!(tr.canonical);
    }

    #[test]
    fn test_parse_gff3_exons() {
        let transcripts = parse_gff3(sample_gff3().as_bytes()).unwrap();
        let tr = &transcripts[0];
        assert_eq!(tr.exons.len(), 3);
        assert_eq!(tr.exons[0].start, 1000);
        assert_eq!(tr.exons[0].rank, 1);
        assert_eq!(tr.exons[1].start, 2000);
        assert_eq!(tr.exons[2].start, 4000);
    }

    #[test]
    fn test_parse_gff3_translation() {
        let transcripts = parse_gff3(sample_gff3().as_bytes()).unwrap();
        let tr = &transcripts[0];
        assert!(tr.is_coding());

        let tl = tr.translation.as_ref().unwrap();
        assert_eq!(tl.stable_id, "ENSP00000001");
        assert_eq!(tl.genomic_start, 1050);
        assert_eq!(tl.genomic_end, 4500);

        assert_eq!(tr.coding_region_start, Some(1050));
        assert_eq!(tr.coding_region_end, Some(4500));
        // cDNA coding start: exon1 starts at 1000, CDS starts at 1050 → offset 50 → cDNA pos 51
        assert_eq!(tr.cdna_coding_start, Some(51));
    }

    #[test]
    fn test_parse_gff3_reverse_strand() {
        let gff = "##gff-version 3
chr2\tensembl\tgene\t1000\t3000\t.\t-\t.\tID=gene:ENSG00000002;Name=REVGENE;biotype=protein_coding
chr2\tensembl\tmRNA\t1000\t3000\t.\t-\t.\tID=transcript:ENST00000002;Parent=gene:ENSG00000002;biotype=protein_coding
chr2\tensembl\texon\t2500\t3000\t.\t-\t.\tID=exon:ENSE00000010;Parent=transcript:ENST00000002
chr2\tensembl\texon\t1000\t1500\t.\t-\t.\tID=exon:ENSE00000011;Parent=transcript:ENST00000002
chr2\tensembl\tCDS\t2500\t2900\t.\t-\t0\tID=CDS:ENSP00000002;Parent=transcript:ENST00000002
chr2\tensembl\tCDS\t1100\t1500\t.\t-\t1\tID=CDS:ENSP00000002;Parent=transcript:ENST00000002";

        let transcripts = parse_gff3(gff.as_bytes()).unwrap();
        assert_eq!(transcripts.len(), 1);
        let tr = &transcripts[0];
        assert_eq!(tr.strand, Strand::Reverse);
        assert!(tr.is_coding());
        // For reverse strand: first exon in transcript order is 2500-3000
        assert_eq!(tr.exons[0].start, 2500);
        assert_eq!(tr.exons[1].start, 1000);
    }

    #[test]
    fn test_url_decode() {
        assert_eq!(url_decode("hello%20world"), "hello world");
        assert_eq!(url_decode("100%25"), "100%");
        assert_eq!(url_decode("normal"), "normal");
    }
}
