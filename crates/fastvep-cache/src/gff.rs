use anyhow::Result;
use fastvep_core::Strand;
use fastvep_genome::{Exon, Gene, Transcript, Translation};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::sync::Arc;

/// Parse a GFF3 file into Transcript models.
///
/// Builds gene -> transcript -> exon/CDS hierarchy from GFF3 features.
///
/// Streams lines from `reader` straight into the parser so a multi-GB
/// Ensembl/GENCODE GFF3 doesn't have to be materialised as `Vec<String>`
/// before parsing begins. Per-line IO errors are surfaced via the
/// `Result<String>` iterator with a 1-based line number for diagnostics
/// (earlier code folded IO errors into empty lines via `unwrap_or_default`,
/// which silently produced an empty transcript set on truncated files).
pub fn parse_gff3<R: Read>(reader: R) -> Result<Vec<Transcript>> {
    parse_gff3_with_source(reader, "GFF3")
}

/// Like `parse_gff3`, but stamps every returned transcript with `source`.
///
/// Used by the merged-cache path so transcripts from different GFF3s
/// (e.g. Ensembl + RefSeq) carry their origin through to the SOURCE
/// column of the output. Without per-file tagging the whole call would
/// collapse to a single "GFF3" label and a merged run would be
/// indistinguishable from a single-source one.
pub fn parse_gff3_with_source<R: Read>(reader: R, source: &str) -> Result<Vec<Transcript>> {
    let buf = BufReader::new(reader);
    let lines = buf
        .lines()
        .enumerate()
        .map(|(i, line)| line.map_err(|e| anyhow::anyhow!("Reading GFF3 line {}: {}", i + 1, e)));
    let mut trs = parse_gff3_lines(lines)?;
    for tr in &mut trs {
        tr.source = Some(source.to_string());
    }
    Ok(trs)
}

/// Parse a tabix-indexed GFF3 file, loading only transcripts overlapping given regions.
///
/// Regions are specified as (chrom, start, end) tuples. This loads gene/transcript/exon/CDS
/// features from the indexed file for each region, then assembles transcripts.
pub fn parse_gff3_indexed(
    gff3_gz_path: &Path,
    regions: &[(String, u64, u64)],
) -> Result<Vec<Transcript>> {
    parse_gff3_indexed_with_source(gff3_gz_path, regions, "GFF3")
}

/// Like `parse_gff3_indexed`, but stamps every transcript with `source`.
pub fn parse_gff3_indexed_with_source(
    gff3_gz_path: &Path,
    regions: &[(String, u64, u64)],
    source: &str,
) -> Result<Vec<Transcript>> {
    let mut trs = parse_gff3_indexed_inner(gff3_gz_path, regions)?;
    for tr in &mut trs {
        tr.source = Some(source.to_string());
    }
    Ok(trs)
}

fn parse_gff3_indexed_inner(
    gff3_gz_path: &Path,
    regions: &[(String, u64, u64)],
) -> Result<Vec<Transcript>> {
    use noodles_bgzf as bgzf;
    use noodles_core::region::Interval;
    use noodles_core::Position;
    use noodles_csi::binning_index::BinningIndex;
    use noodles_tabix as tabix;

    let tbi_path = format!("{}.tbi", gff3_gz_path.display());
    let index = tabix::fs::read(&tbi_path)
        .map_err(|e| anyhow::anyhow!("Reading tabix index {}: {}", tbi_path, e))?;

    let header = index
        .header()
        .ok_or_else(|| anyhow::anyhow!("Missing tabix header"))?;
    let ref_names: Vec<String> = header
        .reference_sequence_names()
        .iter()
        .map(|n| n.to_string())
        .collect();

    let mut all_lines: Vec<String> = Vec::new();
    let mut seen_lines: std::collections::HashSet<u64> = std::collections::HashSet::new();

    for (chrom, start, end) in regions {
        // Find reference ID, trying with/without "chr" prefix
        let ref_id = ref_names.iter().position(|n| n == chrom).or_else(|| {
            if chrom.starts_with("chr") {
                ref_names.iter().position(|n| n == &chrom[3..])
            } else {
                let with_chr = format!("chr{}", chrom);
                ref_names.iter().position(|n| *n == with_chr)
            }
        });

        let ref_id = match ref_id {
            Some(id) => id,
            None => continue,
        };

        let pos_start = Position::try_from((*start).max(1) as usize).unwrap_or(Position::MIN);
        let pos_end = Position::try_from(*end as usize).unwrap_or(Position::MIN);
        let query_interval: Interval = (pos_start..=pos_end).into();

        let chunks = match index.query(ref_id, query_interval) {
            Ok(c) => c,
            Err(_) => continue,
        };

        let file = std::fs::File::open(gff3_gz_path)?;
        let mut reader = bgzf::io::Reader::new(file);

        for chunk in &chunks {
            reader.seek(chunk.start())?;
            let mut line = String::new();

            loop {
                line.clear();
                let bytes = reader.read_line(&mut line)?;
                if bytes == 0 {
                    break;
                }

                let trimmed = line.trim();
                if trimmed.is_empty() || trimmed.starts_with('#') {
                    if reader.virtual_position() >= chunk.end() {
                        break;
                    }
                    continue;
                }

                // Deduplicate lines by hashing (regions may overlap in tabix chunks)
                let hash = {
                    use std::hash::{Hash, Hasher};
                    let mut h = std::collections::hash_map::DefaultHasher::new();
                    trimmed.hash(&mut h);
                    h.finish()
                };
                if seen_lines.insert(hash) {
                    all_lines.push(trimmed.to_string());
                }

                if reader.virtual_position() >= chunk.end() {
                    break;
                }
            }
        }
    }

    // Lines were already collected into memory by the tabix-chunk loop above,
    // so wrap them in `Ok` to satisfy the streaming-aware parser signature.
    parse_gff3_lines(all_lines.into_iter().map(Ok))
}

/// Core GFF3 parsing logic: takes an iterator of `Result<String>` lines and
/// builds Transcript models. Each item is either a raw line or a typed IO
/// error from the underlying reader; the latter aborts parsing immediately
/// rather than being silently treated as an empty line.
fn parse_gff3_lines(lines: impl Iterator<Item = Result<String>>) -> Result<Vec<Transcript>> {
    let mut genes: HashMap<String, GffGene> = HashMap::new();
    let mut transcripts: HashMap<String, GffTranscript> = HashMap::new();
    let mut exons: Vec<GffExon> = Vec::new();
    let mut cds_features: Vec<GffCds> = Vec::new();

    for line in lines {
        let line = line?;
        let line = line.trim().to_string();
        let line = line.as_str();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let seqid = fields[0];
        let feature_type = fields[2];
        // Coordinates are required and integral. Earlier these defaulted to 0
        // on parse failure, silently emitting features at position 0 — which
        // then collided with everything during overlap queries. Skip the
        // line with a warning instead so the rest of the file still parses.
        let Ok(start) = fields[3].parse::<u64>() else {
            log::warn!("GFF3: skipping line with non-numeric start: {}", fields[3]);
            continue;
        };
        let Ok(end) = fields[4].parse::<u64>() else {
            log::warn!("GFF3: skipping line with non-numeric end: {}", fields[4]);
            continue;
        };
        // GFF3 coordinates are 1-based and inclusive, so start must be >= 1
        // and <= end. A `start` of 0 (or start > end) parses fine as a u64
        // but is invalid; feeding it through would underflow the `- 1`
        // conversions to 0-based offsets done downstream (e.g. transcript
        // exon/CDS math), the same class of bug the non-numeric check above
        // guards against.
        if start == 0 || start > end {
            log::warn!(
                "GFF3: skipping line with invalid coordinates (start={}, end={})",
                start,
                end
            );
            continue;
        }
        let strand = match fields[6] {
            "-" => Strand::Reverse,
            _ => Strand::Forward,
        };
        // Phase is allowed to be "." in GFF3 (no phase); use -1 sentinel.
        let phase: i8 = fields[7].parse().unwrap_or(-1);
        let attrs = parse_attributes(fields[8]);

        match feature_type {
            t if is_gene_feature(t) => {
                let gene_id = attrs
                    .get("ID")
                    .or_else(|| attrs.get("gene_id"))
                    .cloned()
                    .unwrap_or_default();
                // Strip "gene:" prefix if present (Ensembl GFF3)
                let gene_id = gene_id
                    .strip_prefix("gene:")
                    .unwrap_or(&gene_id)
                    .to_string();
                let symbol = attrs.get("Name").or_else(|| attrs.get("gene")).cloned();
                let biotype = attrs
                    .get("biotype")
                    .or_else(|| attrs.get("gene_biotype"))
                    .cloned()
                    .unwrap_or_else(|| default_gene_biotype(feature_type));

                genes.insert(
                    gene_id.clone(),
                    GffGene {
                        id: gene_id,
                        symbol,
                        biotype,
                        description: attrs.get("description").cloned(),
                        gene_id_attr: attrs.get("gene_id").cloned(),
                        chromosome: seqid.to_string(),
                        start,
                        end,
                        strand,
                    },
                );
            }
            t if is_known_transcript_term(t) || is_provisional_transcript(t, &attrs) => {
                let transcript_id = attrs
                    .get("ID")
                    .or_else(|| attrs.get("transcript_id"))
                    .cloned()
                    .unwrap_or_default();
                let transcript_id = transcript_id
                    .strip_prefix("transcript:")
                    .unwrap_or(&transcript_id)
                    .to_string();
                let parent = attrs.get("Parent").cloned().unwrap_or_default();
                let parent = parent.strip_prefix("gene:").unwrap_or(&parent).to_string();
                let biotype = attrs
                    .get("biotype")
                    .or_else(|| attrs.get("transcript_type"))
                    .or_else(|| attrs.get("transcript_biotype"))
                    .cloned();
                let public_id = public_transcript_id(&attrs, &transcript_id);
                let gene_id_attr = attrs
                    .get("gene_id")
                    .cloned()
                    .or_else(|| gene_id_from_dbxref(&attrs));
                let tag_str = attrs.get("tag").cloned().unwrap_or_default();
                let canonical = tag_str.contains("Ensembl_canonical");
                let mane_select = tag_str.contains("MANE_Select");
                let mane_plus_clinical = tag_str.contains("MANE_Plus_Clinical");
                let gencode_primary = tag_str.contains("gencode_primary");

                // Parse CDS completeness flags from tags
                let mut flags = Vec::new();
                if tag_str.contains("cds_end_NF") {
                    flags.push("cds_end_NF".to_string());
                }
                if tag_str.contains("cds_start_NF") {
                    flags.push("cds_start_NF".to_string());
                }

                let version: Option<u32> = attrs.get("version").and_then(|v| v.parse().ok());
                let ccds = attrs.get("ccdsid").cloned();

                // Parse TSL: "1 (assigned to previous version 2)" → 1
                let tsl: Option<u8> = attrs
                    .get("transcript_support_level")
                    .and_then(|v| v.split_whitespace().next().and_then(|n| n.parse().ok()));

                // APPRIS rides in the same `tag=` list as MANE and canonical,
                // spelled `appris_principal_1` / `appris_alternative_2`. VEP
                // reports the suffix, so that is what is stored.
                let appris = appris_from_tags(&tag_str);

                transcripts.insert(
                    transcript_id.clone(),
                    GffTranscript {
                        id: transcript_id,
                        public_id,
                        parent_gene: parent,
                        gene_id_attr,
                        biotype,
                        feature_type: feature_type.to_string(),
                        gbkey: attrs.get("gbkey").cloned(),
                        ncrna_class: attrs.get("ncrna_class").cloned(),
                        chromosome: seqid.to_string(),
                        start,
                        end,
                        strand,
                        canonical,
                        mane_select,
                        mane_plus_clinical,
                        gencode_primary,
                        ccds,
                        tsl,
                        appris,
                        flags,
                        version,
                        provisional: !is_known_transcript_term(feature_type),
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
                let exon_id = exon_id
                    .strip_prefix("exon:")
                    .unwrap_or(&exon_id)
                    .to_string();
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
                // Protein ID: prefer explicit protein_id attribute (NCBI format),
                // then fall back to ID with prefix stripping (Ensembl format)
                let protein_id = attrs
                    .get("protein_id")
                    .cloned()
                    .or_else(|| {
                        attrs.get("ID").map(|id| {
                            id.strip_prefix("CDS:")
                                .or_else(|| id.strip_prefix("cds-"))
                                .unwrap_or(id)
                                .to_string()
                        })
                    })
                    .unwrap_or_default();

                let protein_version: Option<u32> =
                    attrs.get("version").and_then(|v| v.parse().ok());

                cds_features.push(GffCds {
                    parent_transcript: parent,
                    protein_id,
                    protein_version,
                    start,
                    end,
                    strand,
                    phase,
                });
            }
            _ => {}
        }
    }

    // For NCBI GFF3 format: create implicit transcripts for genes that have
    // CDS children but no transcript/mRNA children (common in bacterial genomes).
    let _genes_with_transcripts: std::collections::HashSet<String> = transcripts
        .values()
        .map(|t| t.parent_gene.clone())
        .collect();

    for cds in &cds_features {
        let parent = &cds.parent_transcript;
        // If the CDS parent is a gene (not a transcript), create an implicit transcript
        if !transcripts.contains_key(parent) && genes.contains_key(parent) {
            // Strip "gene-" prefix if present (NCBI format)
            let gene_id = parent.strip_prefix("gene-").unwrap_or(parent);
            let gene = &genes[parent];
            let implicit_tid = format!("{}_t1", gene_id);
            if !transcripts.contains_key(&implicit_tid) {
                transcripts.insert(
                    implicit_tid.clone(),
                    GffTranscript {
                        public_id: implicit_tid.clone(),
                        id: implicit_tid,
                        parent_gene: parent.clone(),
                        gene_id_attr: None,
                        biotype: Some(gene.biotype.clone()),
                        feature_type: "gene".to_string(),
                        gbkey: None,
                        ncrna_class: None,
                        chromosome: gene.chromosome.clone(),
                        start: gene.start,
                        end: gene.end,
                        strand: gene.strand,
                        canonical: true,
                        mane_select: false,
                        mane_plus_clinical: false,
                        gencode_primary: false,
                        ccds: None,
                        tsl: None,
                        appris: None,
                        flags: vec![],
                        version: None,
                        // Synthesised from a real CDS, so it is confirmed by
                        // construction; the implicit-exon step below gives it
                        // exons too.
                        provisional: false,
                    },
                );
            }
        }
    }

    // Remap CDS features whose parent is a gene to the implicit transcript
    for cds in &mut cds_features {
        let parent = &cds.parent_transcript;
        if !transcripts.contains_key(parent) {
            // Try with gene- prefix stripped
            let gene_id = parent.strip_prefix("gene-").unwrap_or(parent);
            let implicit_tid = format!("{}_t1", gene_id);
            if transcripts.contains_key(&implicit_tid) {
                cds.parent_transcript = implicit_tid;
            }
        }
    }

    // Create implicit exons for CDS features under implicit transcripts
    // that have no corresponding exon (bacterial genomes where CDS = exon)
    {
        let exon_parents: std::collections::HashSet<&str> =
            exons.iter().map(|e| e.parent_transcript.as_str()).collect();
        let implicit_exons: Vec<GffExon> = cds_features
            .iter()
            .filter(|cds| {
                cds.parent_transcript.ends_with("_t1")
                    && !exon_parents.contains(cds.parent_transcript.as_str())
            })
            .map(|cds| GffExon {
                id: format!("exon_{}_{}_{}", cds.parent_transcript, cds.start, cds.end),
                parent_transcript: cds.parent_transcript.clone(),
                start: cds.start,
                end: cds.end,
                strand: cds.strand,
                phase: cds.phase,
                rank: 0,
            })
            .collect();
        exons.extend(implicit_exons);
    }

    // Pre-index exons and CDS by parent transcript for O(1) lookup per transcript.
    // This replaces the O(T*E) nested scan in the assembly loop.
    let mut exons_by_tx: HashMap<String, Vec<GffExon>> = HashMap::new();
    for exon in exons {
        exons_by_tx
            .entry(exon.parent_transcript.clone())
            .or_default()
            .push(exon);
    }
    let mut cds_by_tx: HashMap<String, Vec<GffCds>> = HashMap::new();
    for cds in cds_features {
        cds_by_tx
            .entry(cds.parent_transcript.clone())
            .or_default()
            .push(cds);
    }

    // Confirm or discard the records admitted on structure alone. What makes a
    // record a transcript is that a spliced feature hangs off it, and that is
    // not known until the whole file has been read - hence the two passes.
    //
    // Both exons and CDS count. Requiring exons alone would re-close the
    // allowlist for the case this rule exists to open: a file that names its
    // transcripts with an SO term not in `is_known_transcript_term` *and*
    // describes them with CDS records only. The implicit-exon step above cannot
    // rescue those, because it fires only for the `_t1` transcripts it
    // synthesised itself, so such a transcript would be admitted here and then
    // dropped - taking every consequence call at that locus with it, which is
    // exactly the failure the provisional rule was added to prevent.
    //
    // A record that carries neither is annotation, not a transcript: reading one
    // as a transcript invents a zero-exon phantom that VEP has no counterpart
    // for, which then reaches the output as its own CSQ entry, inflates the
    // "Loaded N transcripts" count, is written to the cache, and competes in
    // `--pick`. Records matched by a known SO term are untouched here, so
    // nothing that parsed before this rule existed can be dropped by it.
    transcripts.retain(|tid, tr| {
        !tr.provisional
            || exons_by_tx.contains_key(tid.as_str())
            || cds_by_tx.contains_key(tid.as_str())
    });

    // Build transcripts
    let mut result = Vec::with_capacity(transcripts.len());
    let empty_exons: Vec<GffExon> = Vec::new();
    let empty_cds: Vec<GffCds> = Vec::new();

    for (tid, gff_tr) in &transcripts {
        let gene = genes.get(&gff_tr.parent_gene);
        let gene_model = gene
            .map(|g| Gene {
                stable_id: Arc::from(gene_public_id(gff_tr, g).as_str()),
                symbol: g.symbol.as_deref().map(Arc::from),
                symbol_source: None,
                hgnc_id: None,
                biotype: Arc::from(g.biotype.as_str()),
                chromosome: Arc::from(g.chromosome.as_str()),
                start: g.start,
                end: g.end,
                strand: g.strand,
            })
            .unwrap_or_else(|| Gene {
                stable_id: Arc::from(
                    gff_tr
                        .gene_id_attr
                        .as_deref()
                        .unwrap_or(gff_tr.parent_gene.as_str()),
                ),
                symbol: None,
                symbol_source: None,
                hgnc_id: None,
                biotype: Arc::from("unknown"),
                chromosome: Arc::from(gff_tr.chromosome.as_str()),
                start: gff_tr.start,
                end: gff_tr.end,
                strand: gff_tr.strand,
            });

        // Collect exons for this transcript (O(1) lookup via pre-indexed map)
        let mut tr_exons: Vec<Exon> = exons_by_tx
            .get(tid.as_str())
            .unwrap_or(&empty_exons)
            .iter()
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
            Strand::Reverse => tr_exons.sort_by_key(|e| std::cmp::Reverse(e.start)),
        }

        // Assign ranks if not set
        for (i, exon) in tr_exons.iter_mut().enumerate() {
            if exon.rank == 0 {
                exon.rank = (i + 1) as u32;
            }
        }

        // Collect CDS features for this transcript (O(1) lookup via pre-indexed map)
        let tr_cds_owned = cds_by_tx.get(tid.as_str()).unwrap_or(&empty_cds);
        let tr_cds: Vec<&GffCds> = tr_cds_owned.iter().collect();

        // Find the phase of the first CDS in transcript order and convert to Ensembl phase.
        // GFF3 phase → Ensembl phase: 0→0, 1→2, 2→1
        // Ensembl phase on the starting exon indicates how many bases from the
        // previous exon are needed to complete the first codon. For incomplete CDS
        // (cds_start_NF), this means the CDS numbering should account for those
        // "missing" bases, effectively shifting cdna_coding_start earlier.
        let first_cds_ensembl_phase: u64 = if !tr_cds.is_empty() {
            let gff_phase = match gff_tr.strand {
                Strand::Forward => tr_cds
                    .iter()
                    .min_by_key(|c| c.start)
                    .map(|c| c.phase)
                    .unwrap_or(0),
                Strand::Reverse => tr_cds
                    .iter()
                    .max_by_key(|c| c.end)
                    .map(|c| c.phase)
                    .unwrap_or(0),
            };
            match gff_phase {
                1 => 2,
                2 => 1,
                other => other as u64,
            }
        } else {
            0
        };

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

            // tr_exons already sorted by position at line 325-328 above
            for exon in &tr_exons {
                let exon_len = exon.end - exon.start + 1;
                match gff_tr.strand {
                    Strand::Forward => {
                        if cs.is_none()
                            && tl.genomic_start >= exon.start
                            && tl.genomic_start <= exon.end
                        {
                            cs = Some(cdna_pos + (tl.genomic_start - exon.start) + 1);
                        }
                        if ce.is_none()
                            && tl.genomic_end >= exon.start
                            && tl.genomic_end <= exon.end
                        {
                            ce = Some(cdna_pos + (tl.genomic_end - exon.start) + 1);
                        }
                    }
                    Strand::Reverse => {
                        if cs.is_none()
                            && tl.genomic_end >= exon.start
                            && tl.genomic_end <= exon.end
                        {
                            cs = Some(cdna_pos + (exon.end - tl.genomic_end) + 1);
                        }
                        if ce.is_none()
                            && tl.genomic_start >= exon.start
                            && tl.genomic_start <= exon.end
                        {
                            ce = Some(cdna_pos + (exon.end - tl.genomic_start) + 1);
                        }
                    }
                }
                cdna_pos += exon_len;
            }
            // Note: first_cds_ensembl_phase is stored on the Transcript and
            // applied in cdna_to_cds() to account for incomplete CDS starts.

            (cs, ce)
        } else {
            (None, None)
        };

        result.push(Transcript {
            stable_id: Arc::from(gff_tr.public_id.as_str()),
            version: gff_tr.version,
            gene: gene_model,
            biotype: Arc::from(transcript_biotype(gff_tr, gene).as_str()),
            chromosome: Arc::from(gff_tr.chromosome.as_str()),
            start: gff_tr.start,
            end: gff_tr.end,
            strand: gff_tr.strand,
            exons: tr_exons,
            translation,
            cdna_coding_start,
            cdna_coding_end,
            coding_region_start: tr_cds.iter().map(|c| c.start).min(),
            coding_region_end: tr_cds.iter().map(|c| c.end).max(),
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: gff_tr.canonical,
            mane_select: if gff_tr.mane_select {
                let vid = match gff_tr.version {
                    Some(v) => format!("{}.{}", gff_tr.public_id, v),
                    None => tid.to_string(),
                };
                Some(vid)
            } else {
                None
            },
            mane_plus_clinical: if gff_tr.mane_plus_clinical {
                let vid = match gff_tr.version {
                    Some(v) => format!("{}.{}", gff_tr.public_id, v),
                    None => tid.to_string(),
                };
                Some(vid)
            } else {
                None
            },
            tsl: gff_tr.tsl,
            appris: gff_tr.appris.clone(),
            ccds: gff_tr.ccds.clone(),
            protein_id: tr_cds.first().map(|c| c.protein_id.clone()),
            // The version of the first CDS record, which is where Ensembl's GFF3 puts a protein's
            // (`protein_id=ENSP00000363557;version=3`), and 1 where it has none or says 0: Ensembl
            // VEP 116 reads it so (`$ordered_cdss->[0]->{attributes}->{version} || 1`). VEP 115
            // gave every translation version 1, and so did this, which put `ENSP….1` in front of
            // the `HGVSp` of one protein in six of Ensembl 115's: 184 of the 844 `HGVSp` rows of
            // GIAB HG002 differed from VEP 116.3 in that suffix and nothing else.
            protein_version: tr_cds
                .first()
                .map(|c| c.protein_version.filter(|&v| v != 0).unwrap_or(1)),
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: Some("GFF3".into()),
            gencode_primary: gff_tr.gencode_primary,
            flags: gff_tr.flags.clone(),
            codon_table_start_phase: first_cds_ensembl_phase,
        });
    }

    // HashMap iteration order is randomized per process, so sort by stable_id
    // to make cache serialization bit-for-bit reproducible across builds.
    result.sort_by(|a, b| {
        a.stable_id
            .cmp(&b.stable_id)
            .then_with(|| a.chromosome.cmp(&b.chromosome))
            .then_with(|| a.start.cmp(&b.start))
            .then_with(|| a.end.cmp(&b.end))
    });

    Ok(result)
}

/// What the output calls a transcript: its `transcript_id` attribute, and its ID where it has none,
/// with a leading `gene:` or `transcript:` taken off. That is how Ensembl VEP names it
/// (`_create_transcript`, `BaseGXF.pm`), and it is what makes a RefSeq record `NM_000001.3` and
/// not `rna-NM_000001.3`. A record with no `transcript_id` keeps its ID as it is: NCBI's tRNA
/// records have none, and VEP calls one `rna-TRG-GCC1-5`.
fn public_transcript_id(attrs: &HashMap<String, String>, id: &str) -> String {
    let named = attrs.get("transcript_id").map_or(id, String::as_str);
    ["gene:", "transcript:"]
        .iter()
        .find_map(|prefix| named.strip_prefix(prefix))
        .unwrap_or(named)
        .to_string()
}

/// The `GeneID` of a record's `Dbxref` (`Dbxref=GeneID:100189284,HGNC:HGNC:34850`), which is how
/// NCBI says what gene a transcript is of, and so what Ensembl VEP calls the gene (`_add_identifiers`).
fn gene_id_from_dbxref(attrs: &HashMap<String, String>) -> Option<String> {
    let dbxref = attrs.get("dbxref").or_else(|| attrs.get("Dbxref"))?;
    dbxref.split(',').find_map(|pair| {
        let (key, value) = pair.split_once(':')?;
        (key == "GeneID").then(|| value.split(':').next().unwrap_or(value).to_string())
    })
}

/// What the output calls a gene, in the order Ensembl VEP looks (`_add_identifiers`): the `gene_id`
/// or the Dbxref `GeneID` of its transcript, then its own record's `gene_id`, then its ID. An NCBI
/// gene is `100189284` and not `gene-TRG-GCC1-5`; an Ensembl one is the same by any of them.
fn gene_public_id(transcript: &GffTranscript, gene: &GffGene) -> String {
    transcript
        .gene_id_attr
        .clone()
        .or_else(|| gene.gene_id_attr.clone())
        .unwrap_or_else(|| gene.id.clone())
}

/// A transcript's biotype as Ensembl VEP gets it (`_record_get_biotype`, `GFF.pm`): the one the
/// record states, and for a record that does not (RefSeq's) a guess from its SO term and what
/// else it says: an mRNA is `protein_coding`, an ncRNA is its `ncrna_class`, a `V_gene_segment`
/// is `IG_V_gene`, a transcript of a gene described as a microRNA is a `miRNA`, and any other
/// is its `gbkey` (`ncRNA` for an lnc_RNA, `misc_RNA`, `rRNA`, `tRNA`).
///
/// Where none of that says anything VEP drops the transcript with a warning. This keeps it, named
/// by its SO term as it always was: a transcript fastVEP annotates and VEP does not is a row
/// more, and a transcript it drops is a consequence lost.
fn transcript_biotype(transcript: &GffTranscript, gene: Option<&GffGene>) -> String {
    if let Some(stated) = &transcript.biotype {
        return stated.clone();
    }
    let kind = transcript.feature_type.to_ascii_lowercase();
    let guessed = if kind == "mrna" {
        Some("protein_coding".to_string())
    } else if kind == "ncrna" {
        transcript.ncrna_class.clone()
    } else if let Some(segment) = kind
        .strip_suffix("_gene_segment")
        .filter(|s| !s.is_empty() && s.bytes().all(|b| b.is_ascii_lowercase()))
    {
        Some(format!("IG_{}_gene", segment.to_ascii_uppercase()))
    } else if gene
        .and_then(|g| g.description.as_deref())
        .is_some_and(|d| d.starts_with("microRNA"))
    {
        Some("miRNA".to_string())
    } else {
        transcript.gbkey.clone()
    };
    guessed.unwrap_or_else(|| transcript.feature_type.clone())
}

/// True when a GFF3 record carries a gene, whichever SO term the file uses.
///
/// Ensembl types non-coding loci `ncRNA_gene` and pseudogene loci
/// `pseudogene`, keeping plain `gene` for protein-coding ones, and adds
/// subtypes as new biotypes appear. Matching only `gene`/`pseudogene` dropped
/// every `ncRNA_gene` record, so lncRNA/miRNA/snRNA transcripts fell back to a
/// gene synthesised from the transcript's `Parent`: right stable ID, no symbol
/// and biotype "unknown". That is why SYMBOL came out empty where VEP printed
/// WASH7P (#98). Match the whole `*_gene` family instead; the V/D/J/C
/// `*_gene_segment` transcript types end in `_segment` and are not caught.
fn is_gene_feature(feature_type: &str) -> bool {
    matches!(feature_type, "gene" | "pseudogene") || feature_type.ends_with("_gene")
}

/// Biotype to assume for a gene record with no `biotype`/`gene_biotype`.
///
/// Only reached when the file omits the attribute; every gene record in
/// Ensembl GRCh38.115 and RefSeq GRCh38.p14 carries one, so this is the path
/// for sparser sources rather than the common case.
///
/// A bare `gene` in a file that omits the attribute (NCBI) is almost always
/// protein-coding. Every other SO term in the `*_gene` family already names its
/// own biotype, so dropping the suffix is both the answer and the reason:
/// `ncRNA_gene` -> `ncRNA`, `miRNA_gene` -> `miRNA`,
/// `transposable_element_gene` -> `transposable_element`. Defaulting those to
/// protein_coding - which is what a fixed fallback did before they were read at
/// all - would mislabel every non-coding locus. `pseudogene` has no `_gene`
/// suffix and stays as it is.
fn default_gene_biotype(feature_type: &str) -> String {
    match feature_type {
        "gene" => "protein_coding".to_string(),
        t => t.strip_suffix("_gene").unwrap_or(t).to_string(),
    }
}

/// SO terms known to name a transcript.
///
/// What Ensembl and RefSeq emit today. A record matching one of these is taken
/// at its word; anything else needs [`is_provisional_transcript`] to earn it.
fn is_known_transcript_term(feature_type: &str) -> bool {
    matches!(
        feature_type,
        "mRNA"
            | "transcript"
            | "lnc_RNA"
            | "miRNA"
            | "snRNA"
            | "snoRNA"
            | "rRNA"
            | "ncRNA"
            | "tRNA"
            | "scRNA"
            | "V_gene_segment"
            | "D_gene_segment"
            | "J_gene_segment"
            | "C_gene_segment"
            | "NMD_transcript_variant"
            | "pseudogenic_transcript"
    )
}

/// True when an unrecognised record is worth *provisionally* reading as a
/// transcript, to be confirmed or discarded once its children are known.
///
/// Both vocabularies grow, and an unlisted term costs more than the
/// `ncRNA_gene` gap in [`is_gene_feature`] did: the transcript is dropped
/// outright, taking every consequence call at that locus with it rather than
/// just its name. Measured against the releases current when this was written,
/// the terms this admits that the explicit list does not carry are:
///
/// - RefSeq GRCh38.p14 (whole genome): `primary_transcript` (2,139),
///   `scaRNA` (54), `antisense_RNA` (42), `vault_RNA` (4), `Y_RNA` (4),
///   `RNase_P_RNA` (2), `telomerase_RNA` (1), `RNase_MRP_RNA` (1)
/// - Ensembl GRCh38.115: `unconfirmed_transcript` (40 on chr1 alone)
/// - FlyBase BDGP6: `transposable_element` (5,896), `pre_miRNA` (262)
///
/// Hanging off a gene is necessary but *not* sufficient, which is why the answer
/// here is only provisional: a GFF3 may hang plain annotation records off a gene
/// as well as transcripts, and reading one as a transcript invents a zero-exon
/// feature that VEP has no counterpart for, which then shows up as an extra CSQ
/// entry and as a candidate for `--pick`. What actually distinguishes a
/// transcript is that a spliced feature hangs off *it*, and that is not known
/// until the file has been read; see the `retain` in [`parse_gff3_lines`].
///
/// That is not hypothetical. RefSeq GRCh38.p14 has six `enhancer` records at
/// the immunoglobulin heavy chain locus - `ID=id-IGH`, `Parent=gene-IGH`,
/// `gbkey=regulatory` - with no exon and no CDS. They are the only six of its
/// 199,249 gene-parented records that carry neither, and without the `retain`
/// each becomes a transcript with biotype `enhancer` that reports
/// `non_coding_transcript_variant` at a locus VEP calls intergenic.
///
/// It is *not* `sequence_feature` that does this, contrary to what an earlier
/// version of this comment said: RefSeq has 1,964 of those and not one carries
/// a `Parent`.
///
/// `exon` and `CDS` are excluded up front because RefSeq attaches those straight
/// to a gene for loci with no RNA record, and the arms below turn them into an
/// implicit transcript instead.
fn is_provisional_transcript(feature_type: &str, attrs: &HashMap<String, String>) -> bool {
    !matches!(feature_type, "exon" | "CDS")
        && parent_is_gene(attrs)
        && (attrs.contains_key("ID") || attrs.contains_key("transcript_id"))
}

/// True when `Parent` names a gene record: `gene:ENSG...` in Ensembl GFF3,
/// `gene-XXXX` in RefSeq. A transcript's own children name it instead
/// (`transcript:`/`rna-`), so this tells the two levels apart without relying on
/// the SO term.
fn parent_is_gene(attrs: &HashMap<String, String>) -> bool {
    attrs
        .get("Parent")
        .is_some_and(|p| p.starts_with("gene:") || p.starts_with("gene-"))
}

/// The APPRIS annotation carried in a GENCODE GFF3 `tag=` list, in the short
/// form VEP reports.
///
/// GENCODE spells it `appris_principal_1` .. `appris_principal_5` and
/// `appris_alternative_1` / `appris_alternative_2`, in the same comma-separated
/// `tag=` attribute that carries `MANE_Select` and `Ensembl_canonical`. VEP
/// reports the short form, `P1`..`P5` and `ALT1`/`ALT2`, and that is what goes
/// in the APPRIS column and into `--pick-order`'s APPRIS tier.
///
/// Older GENCODE releases write the un-numbered `appris_principal` /
/// `appris_alternative`; those become `P` / `ALT`, which `pick`'s APPRIS rank
/// orders after every numbered principal and alternative respectively.
///
/// Ensembl's own GFF3 carries no APPRIS tag at all, so this returns `None`
/// there and the tier stays inert for that source - which is what it was for
/// every source before this function existed.
fn appris_from_tags(tag_str: &str) -> Option<String> {
    for tag in tag_str.split(',') {
        let tag = tag.trim().to_ascii_lowercase();
        let Some(rest) = tag.strip_prefix("appris_") else {
            continue;
        };
        // `principal_1` / `principal1` / `principal`, and the same for
        // `alternative`. The separator has moved between GENCODE releases.
        let (short, digits) = if let Some(d) = rest.strip_prefix("principal") {
            ("P", d)
        } else if let Some(d) = rest.strip_prefix("alternative") {
            ("ALT", d)
        } else {
            continue;
        };
        let digits = digits.trim_start_matches(['_', '-']);
        return Some(if digits.is_empty() {
            short.to_string()
        } else {
            format!("{}{}", short, digits)
        });
    }
    None
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
    /// The record's `description`: Ensembl VEP calls a transcript of a gene described as a
    /// microRNA a `miRNA` when the transcript says nothing of its own ([`transcript_biotype`]).
    description: Option<String>,
    /// The record's `gene_id` attribute, which names the gene in the output when its transcripts
    /// do not ([`gene_public_id`]).
    gene_id_attr: Option<String>,
    chromosome: String,
    start: u64,
    end: u64,
    strand: Strand,
}

#[derive(Debug)]
#[allow(dead_code)]
struct GffTranscript {
    /// What exons and CDS records name as their parent: the record's ID, prefix and all.
    id: String,
    /// What the output calls it ([`public_transcript_id`]).
    public_id: String,
    parent_gene: String,
    /// The `GeneID` of the record's Dbxref, or its `gene_id`: what the output calls its gene.
    gene_id_attr: Option<String>,
    /// The biotype the record states (`biotype`, `transcript_type`, `transcript_biotype`), and
    /// the three things Ensembl VEP guesses one from when it does not ([`transcript_biotype`]).
    biotype: Option<String>,
    feature_type: String,
    gbkey: Option<String>,
    ncrna_class: Option<String>,
    chromosome: String,
    start: u64,
    end: u64,
    strand: Strand,
    canonical: bool,
    mane_select: bool,
    mane_plus_clinical: bool,
    gencode_primary: bool,
    ccds: Option<String>,
    tsl: Option<u8>,
    appris: Option<String>,
    flags: Vec<String>,
    version: Option<u32>,
    /// Admitted only because it hangs off a gene, its SO term being unknown.
    /// Kept only if exons turn out to claim it — see the `retain` in
    /// [`parse_gff3_lines`] and [`is_provisional_transcript`].
    provisional: bool,
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
    protein_version: Option<u32>,
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
        assert_eq!(&*tr.stable_id, "ENST00000001");
        assert_eq!(&*tr.gene.stable_id, "ENSG00000001");
        assert_eq!(tr.gene.symbol.as_deref(), Some("TESTGENE"));
        assert_eq!(&*tr.biotype, "protein_coding");
        assert_eq!(&*tr.chromosome, "chr1");
        assert_eq!(tr.start, 1000);
        assert_eq!(tr.end, 5000);
        assert_eq!(tr.strand, Strand::Forward);
        assert!(tr.canonical);
    }

    /// Ensembl GFF3 types non-coding gene records with SO term `ncRNA_gene`,
    /// not `gene` (issue #98). Dropping those records left every lncRNA /
    /// miRNA / snRNA transcript with a synthesised gene carrying no symbol
    /// and biotype "unknown", so SYMBOL came out empty in the CSQ column
    /// while VEP filled it in.
    /// GENCODE carries APPRIS in the same comma-separated `tag=` list as
    /// `MANE_Select` and `Ensembl_canonical`, and the parser read every other
    /// tag in that list and dropped this one. `Transcript::appris` was
    /// therefore `None` for every transcript from every source, so the APPRIS
    /// column of the CSQ line was always empty and `--pick`'s APPRIS tier
    /// scored every candidate identically.
    #[test]
    fn test_parse_gff3_reads_the_gencode_appris_tag() {
        let gff = "##gff-version 3
1\ttest\tgene\t1000\t2000\t.\t+\t.\tID=gene:ENSG_A;Name=GA;biotype=protein_coding
1\ttest\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:ENST_P;Parent=gene:ENSG_A;biotype=protein_coding;tag=gencode_basic,Ensembl_canonical,appris_principal_1
1\ttest\texon\t1000\t2000\t.\t+\t.\tID=exon:E_P;Parent=transcript:ENST_P;rank=1
1\ttest\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:ENST_A;Parent=gene:ENSG_A;biotype=protein_coding;tag=gencode_basic,appris_alternative_2
1\ttest\texon\t1000\t2000\t.\t+\t.\tID=exon:E_A;Parent=transcript:ENST_A;rank=1
1\ttest\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:ENST_N;Parent=gene:ENSG_A;biotype=protein_coding;tag=gencode_basic
1\ttest\texon\t1000\t2000\t.\t+\t.\tID=exon:E_N;Parent=transcript:ENST_N;rank=1";

        let transcripts = parse_gff3(gff.as_bytes()).unwrap();
        let appris_of = |id: &str| {
            transcripts
                .iter()
                .find(|t| &*t.stable_id == id)
                .unwrap_or_else(|| panic!("{id} missing"))
                .appris
                .clone()
        };
        assert_eq!(appris_of("ENST_P").as_deref(), Some("P1"));
        assert_eq!(appris_of("ENST_A").as_deref(), Some("ALT2"));
        // No APPRIS tag is not "unannotated principal": it stays absent, which
        // is what Ensembl's own GFF3 gives for every transcript.
        assert_eq!(appris_of("ENST_N"), None);
    }

    /// Every spelling the tag has had, and the two ways it can be absent.
    #[test]
    fn test_appris_from_tags_spellings() {
        let cases = [
            ("appris_principal_1", Some("P1")),
            ("appris_principal1", Some("P1")),
            ("appris_principal", Some("P")),
            ("appris_alternative_2", Some("ALT2")),
            ("appris_alternative2", Some("ALT2")),
            ("appris_alternative", Some("ALT")),
            ("gencode_basic,MANE_Select,appris_principal_4", Some("P4")),
            ("APPRIS_PRINCIPAL_3", Some("P3")),
            ("gencode_basic,Ensembl_canonical", None),
            ("", None),
            // Not APPRIS, and close enough to catch a prefix match that is too
            // eager: `appris_` is required, and so is a class after it.
            ("apprisal", None),
            ("appris_something_else", None),
        ];
        for (tags, want) in cases {
            assert_eq!(appris_from_tags(tags).as_deref(), want, "tag list {tags:?}");
        }
    }

    #[test]
    fn test_parse_gff3_ncrna_gene_keeps_symbol() {
        let gff = "##gff-version 3
1\thavana\tncRNA_gene\t14696\t24886\t.\t-\t.\tID=gene:ENSG00000310526;Name=WASH7P;biotype=lncRNA;gene_id=ENSG00000310526
1\thavana\tlnc_RNA\t14696\t24886\t.\t-\t.\tID=transcript:ENST00000488147;Parent=gene:ENSG00000310526;biotype=lncRNA;tag=Ensembl_canonical
1\thavana\texon\t24738\t24886\t.\t-\t.\tID=exon:ENSE00000001;Parent=transcript:ENST00000488147;rank=1
1\thavana\texon\t14696\t14829\t.\t-\t.\tID=exon:ENSE00000002;Parent=transcript:ENST00000488147;rank=2";

        let transcripts = parse_gff3(gff.as_bytes()).unwrap();
        assert_eq!(transcripts.len(), 1);
        let tr = &transcripts[0];
        assert_eq!(&*tr.gene.stable_id, "ENSG00000310526");
        assert_eq!(tr.gene.symbol.as_deref(), Some("WASH7P"));
        assert_eq!(&*tr.gene.biotype, "lncRNA");
        assert_eq!(tr.gene.start, 14696);
        assert_eq!(tr.gene.end, 24886);
        assert_eq!(tr.gene.strand, Strand::Reverse);
    }

    /// `ends_with("_gene")` must not claim the V/D/J/C segment types: those are
    /// transcripts (`*_gene_segment`), and treating one as a gene would drop
    /// every IG/TR transcript.
    #[test]
    fn test_parse_gff3_gene_segment_is_a_transcript() {
        let gff = "##gff-version 3
14\tensembl_havana\tgene\t22438547\t22438554\t.\t+\t.\tID=gene:ENSG00000211821;Name=TRDV2;biotype=IG_V_gene
14\tensembl_havana\tV_gene_segment\t22438547\t22438554\t.\t+\t.\tID=transcript:ENST00000390469;Parent=gene:ENSG00000211821;biotype=IG_V_gene
14\tensembl_havana\texon\t22438547\t22438554\t.\t+\t.\tID=exon:ENSE00001;Parent=transcript:ENST00000390469;rank=1";

        let transcripts = parse_gff3(gff.as_bytes()).unwrap();
        assert_eq!(transcripts.len(), 1);
        let tr = &transcripts[0];
        assert_eq!(&*tr.stable_id, "ENST00000390469");
        assert_eq!(tr.gene.symbol.as_deref(), Some("TRDV2"));
        assert_eq!(tr.exons.len(), 1);
    }

    /// RefSeq names transcripts with SO terms the explicit list does not carry
    /// (`antisense_RNA`, `primary_transcript`, ...). Such a record hangs off a
    /// gene and has its own ID, so keep it rather than dropping the locus.
    #[test]
    fn test_parse_gff3_keeps_unlisted_transcript_term() {
        let gff = "##gff-version 3
NC_000001.11\tBestRefSeq\tgene\t1000\t2000\t.\t+\t.\tID=gene-TESTAS;Name=TESTAS;gene_biotype=lncRNA
NC_000001.11\tBestRefSeq\tantisense_RNA\t1000\t2000\t.\t+\t.\tID=rna-NR_000001.1;Parent=gene-TESTAS
NC_000001.11\tBestRefSeq\texon\t1000\t2000\t.\t+\t.\tID=exon-1;Parent=rna-NR_000001.1";

        let transcripts = parse_gff3(gff.as_bytes()).unwrap();
        assert_eq!(transcripts.len(), 1);
        let tr = &transcripts[0];
        assert_eq!(&*tr.stable_id, "rna-NR_000001.1");
        assert_eq!(tr.gene.symbol.as_deref(), Some("TESTAS"));
        // No biotype/transcript_biotype attribute: fall back to the SO term.
        assert_eq!(&*tr.biotype, "antisense_RNA");
        assert_eq!(tr.exons.len(), 1);
    }

    /// A CDS parented straight to a gene (RefSeq) still becomes an implicit
    /// transcript — the structural transcript test must not claim `exon`/`CDS`
    /// records just because their parent is a gene.
    #[test]
    fn test_parse_gff3_gene_parented_cds_still_implicit() {
        let gff = "##gff-version 3
NC_000001.11\tRefSeq\tgene\t1000\t1300\t.\t+\t.\tID=gene-TESTG;Name=TESTG;gene_biotype=protein_coding
NC_000001.11\tRefSeq\tCDS\t1000\t1300\t.\t+\t0\tID=cds-NP_000001.1;Parent=gene-TESTG;protein_id=NP_000001.1";

        let transcripts = parse_gff3(gff.as_bytes()).unwrap();
        assert_eq!(transcripts.len(), 1);
        let tr = &transcripts[0];
        assert_eq!(&*tr.stable_id, "TESTG_t1");
        assert_eq!(&*tr.gene.stable_id, "gene-TESTG");
        assert_eq!(tr.gene.symbol.as_deref(), Some("TESTG"));
        assert!(tr.is_coding());
        assert_eq!(
            tr.exons.len(),
            1,
            "CDS should have produced an implicit exon"
        );
    }

    /// With no `biotype`/`gene_biotype` attribute the gene biotype comes from
    /// the SO term. `ncRNA_gene` must not inherit the protein_coding default
    /// that bare `gene` gets.
    #[test]
    fn test_gene_biotype_defaults_follow_the_so_term() {
        assert_eq!(default_gene_biotype("gene"), "protein_coding");
        assert_eq!(default_gene_biotype("pseudogene"), "pseudogene");
        assert_eq!(default_gene_biotype("ncRNA_gene"), "ncRNA");
        assert_eq!(default_gene_biotype("miRNA_gene"), "miRNA");
        // FlyBase and SGD, where this family is not hypothetical.
        assert_eq!(
            default_gene_biotype("transposable_element_gene"),
            "transposable_element"
        );

        let gff = "##gff-version 3
1\thavana\tncRNA_gene\t1000\t2000\t.\t+\t.\tID=gene:ENSG00000000001;Name=NCG
1\thavana\tlnc_RNA\t1000\t2000\t.\t+\t.\tID=transcript:ENST00000000001;Parent=gene:ENSG00000000001
1\thavana\texon\t1000\t2000\t.\t+\t.\tID=exon:ENSE00000001;Parent=transcript:ENST00000000001;rank=1";

        let transcripts = parse_gff3(gff.as_bytes()).unwrap();
        assert_eq!(transcripts.len(), 1);
        assert_eq!(&*transcripts[0].gene.biotype, "ncRNA");
    }

    /// A gene-parented record with an ID but no exons and no CDS is annotation,
    /// not a transcript. Admitting one invents a zero-exon feature that VEP has
    /// no counterpart for, and it then shows up as an extra CSQ entry, in the
    /// "Loaded N transcripts" count, in the cache, and among the candidates
    /// `--pick` chooses from.
    ///
    /// The `enhancer` record below is copied from RefSeq GRCh38.p14, where six
    /// of them sit at the immunoglobulin heavy chain locus and are the only six
    /// of its 199,249 gene-parented records with no child of any kind. Without
    /// the confirmation pass this record parses as a transcript of biotype
    /// `enhancer` and reports `non_coding_transcript_variant` where VEP says
    /// intergenic. It is the reason the `provisional` flag exists; see
    /// [`is_provisional_transcript`].
    #[test]
    fn test_parse_gff3_ignores_gene_parented_annotation_records() {
        let gff = "##gff-version 3
NC_000014.9\tBestRefSeq\tgene\t105550000\t105900000\t.\t-\t.\tID=gene-IGH;Name=IGH;gene_biotype=C_region
NC_000014.9\tBestRefSeq\tmRNA\t105860000\t105870000\t.\t-\t.\tID=rna-NM_IGHM;Parent=gene-IGH;transcript_id=NM_IGHM
NC_000014.9\tBestRefSeq\texon\t105860000\t105870000\t.\t-\t.\tID=exon-NM_IGHM-1;Parent=rna-NM_IGHM
NC_000014.9\tCurated Genomic\tenhancer\t105686122\t105686482\t.\t-\t.\tID=id-IGH;Parent=gene-IGH;Note=HS4 enhancer;gbkey=regulatory;regulatory_class=enhancer";

        let transcripts = parse_gff3(gff.as_bytes()).unwrap();
        let ids: Vec<&str> = transcripts.iter().map(|t| &*t.stable_id).collect();
        assert_eq!(
            ids,
            vec!["NM_IGHM"],
            "a gene-parented regulatory record became a phantom transcript"
        );
    }

    /// A transcript may be described by CDS records alone, with no `exon` line.
    /// When its SO term is also one the explicit list does not carry, both new
    /// rules apply at once: it is admitted provisionally, and it must survive
    /// confirmation on its CDS.
    ///
    /// Requiring exons there would have re-closed the allowlist for exactly the
    /// case the provisional rule exists to open - the implicit-exon step cannot
    /// rescue it, because that fires only for the `_t1` transcripts it
    /// synthesised itself - and the locus would lose every consequence call,
    /// not just its name.
    #[test]
    fn test_parse_gff3_keeps_cds_only_transcript_with_unlisted_term() {
        let gff = "##gff-version 3
NC_000001.11\tRefSeq\tgene\t1000\t1300\t.\t+\t.\tID=gene-SRP;Name=SRPGENE;gene_biotype=protein_coding
NC_000001.11\tRefSeq\tSRP_RNA\t1000\t1300\t.\t+\t.\tID=rna-SRP1;Parent=gene-SRP
NC_000001.11\tRefSeq\tCDS\t1000\t1300\t.\t+\t0\tID=cds-SRP1;Parent=rna-SRP1;protein_id=NP_SRP1";

        let transcripts = parse_gff3(gff.as_bytes()).unwrap();
        let ids: Vec<&str> = transcripts.iter().map(|t| &*t.stable_id).collect();
        assert_eq!(
            ids,
            vec!["rna-SRP1"],
            "a CDS-only transcript with an unlisted SO term was dropped"
        );
        assert_eq!(transcripts[0].gene.symbol.as_deref(), Some("SRPGENE"));
        assert!(
            transcripts[0].is_coding(),
            "its CDS should have survived the confirmation pass"
        );
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

    /// A RefSeq GFF3 is named by Ensembl VEP's rules and not by its IDs (`BaseGXF.pm` and `GFF.pm`,
    /// release/116, and a run of VEP 116.3 over `fixtures/vep_refseq_gff`): a transcript is its
    /// `transcript_id` (`NM_090001.3`, where its ID is `rna-NM_090001.3`), or its ID as it is where
    /// it has none (an NCBI tRNA record has none); a gene is the `GeneID` of its transcript's Dbxref
    /// (`90001`, where its ID is `gene-NB1`); and a biotype the record does not state is guessed.
    #[test]
    fn refseq_records_are_named_as_ensembl_vep_names_them() {
        let gff = "##gff-version 3
chr1\tBestRefSeq\tgene\t1000\t1239\t.\t+\t.\tID=gene-NB1;Dbxref=GeneID:90001,HGNC:HGNC:90001;Name=NB1;description=test gene;gbkey=Gene;gene=NB1;gene_biotype=protein_coding
chr1\tBestRefSeq\tmRNA\t1000\t1239\t.\t+\t.\tID=rna-NM_090001.3;Parent=gene-NB1;Dbxref=GeneID:90001,GenBank:NM_090001.3;Name=NM_090001.3;gbkey=mRNA;gene=NB1;transcript_id=NM_090001.3
chr1\tBestRefSeq\texon\t1000\t1239\t.\t+\t.\tID=exon-NM_090001.3-1;Parent=rna-NM_090001.3;gbkey=mRNA;gene=NB1;transcript_id=NM_090001.3
chr1\tBestRefSeq\tCDS\t1000\t1239\t.\t+\t0\tID=cds-NP_090001.2;Parent=rna-NM_090001.3;gbkey=CDS;gene=NB1;protein_id=NP_090001.2
chr1\tBestRefSeq\tgene\t2000\t2239\t.\t+\t.\tID=gene-NB3;Dbxref=GeneID:90003;Name=NB3;gbkey=Gene;gene=NB3;gene_biotype=lncRNA
chr1\tBestRefSeq\tlnc_RNA\t2000\t2239\t.\t+\t.\tID=rna-NR_090003.1;Parent=gene-NB3;Dbxref=GeneID:90003,GenBank:NR_090003.1;gbkey=ncRNA;gene=NB3;transcript_id=NR_090003.1
chr1\tBestRefSeq\texon\t2000\t2239\t.\t+\t.\tID=exon-NR_090003.1-1;Parent=rna-NR_090003.1;transcript_id=NR_090003.1
chr1\tBestRefSeq\tgene\t3000\t3239\t.\t+\t.\tID=gene-MIR5;Dbxref=GeneID:90005;Name=MIR5;description=microRNA 5;gbkey=Gene;gene=MIR5;gene_biotype=miRNA
chr1\tBestRefSeq\tprimary_transcript\t3000\t3239\t.\t+\t.\tID=rna-NR_090005.1;Parent=gene-MIR5;Dbxref=GeneID:90005,GenBank:NR_090005.1;gbkey=precursor_RNA;gene=MIR5;transcript_id=NR_090005.1
chr1\tBestRefSeq\texon\t3000\t3239\t.\t+\t.\tID=exon-NR_090005.1-1;Parent=rna-NR_090005.1;transcript_id=NR_090005.1
chr1\tBestRefSeq\tgene\t4000\t4239\t.\t+\t.\tID=gene-NB8;Dbxref=GeneID:90008;Name=NB8;gbkey=Gene;gene=NB8;gene_biotype=other
chr1\tBestRefSeq\tncRNA\t4000\t4239\t.\t+\t.\tID=rna-NR_090008.1;Parent=gene-NB8;Dbxref=GeneID:90008;gbkey=ncRNA;gene=NB8;ncrna_class=scaRNA;transcript_id=NR_090008.1
chr1\tBestRefSeq\texon\t4000\t4239\t.\t+\t.\tID=exon-NR_090008.1-1;Parent=rna-NR_090008.1;transcript_id=NR_090008.1
chr1\ttRNAscan-SE\tgene\t5000\t5071\t.\t+\t.\tID=gene-TRN9;Dbxref=GeneID:90009;Name=TRN9;gbkey=Gene;gene=TRN9;gene_biotype=tRNA
chr1\ttRNAscan-SE\ttRNA\t5000\t5071\t.\t+\t.\tID=rna-TRN9;Parent=gene-TRN9;Dbxref=GeneID:90009;gbkey=tRNA;gene=TRN9
chr1\ttRNAscan-SE\texon\t5000\t5071\t.\t+\t.\tID=exon-TRN9-1;Parent=rna-TRN9;gbkey=tRNA;gene=TRN9
chr1\tGnomon\tgene\t6000\t6239\t.\t+\t.\tID=gene-LOC90010;Dbxref=GeneID:90010;Name=LOC90010;gbkey=Gene;gene=LOC90010;gene_biotype=C_region
chr1\tGnomon\tV_gene_segment\t6000\t6239\t.\t+\t.\tID=id-LOC90010;Parent=gene-LOC90010;Dbxref=GeneID:90010;gbkey=V_segment;gene=LOC90010
chr1\tGnomon\texon\t6000\t6239\t.\t+\t.\tID=exon-id-LOC90010-1;Parent=id-LOC90010;gbkey=V_segment;gene=LOC90010";

        let transcripts = parse_gff3(gff.as_bytes()).unwrap();
        let named = |id: &str| {
            let t = transcripts
                .iter()
                .find(|t| &*t.stable_id == id)
                .unwrap_or_else(|| {
                    panic!(
                        "no transcript {id}: {:?}",
                        transcripts
                            .iter()
                            .map(|t| &*t.stable_id)
                            .collect::<Vec<_>>()
                    )
                });
            (
                t.gene.stable_id.to_string(),
                t.gene.symbol.as_deref().map(str::to_string),
                t.biotype.to_string(),
            )
        };
        let want = |gene: &str, symbol: &str, biotype: &str| {
            (
                gene.to_string(),
                Some(symbol.to_string()),
                biotype.to_string(),
            )
        };
        assert_eq!(
            named("NM_090001.3"),
            want("90001", "NB1", "protein_coding"),
            "an mRNA"
        );
        assert_eq!(
            named("NR_090003.1"),
            want("90003", "NB3", "ncRNA"),
            "an lnc_RNA is its gbkey"
        );
        assert_eq!(
            named("NR_090005.1"),
            want("90005", "MIR5", "miRNA"),
            "a microRNA's gene says so"
        );
        assert_eq!(
            named("NR_090008.1"),
            want("90008", "NB8", "scaRNA"),
            "an ncRNA is its ncrna_class"
        );
        assert_eq!(
            named("rna-TRN9"),
            want("90009", "TRN9", "tRNA"),
            "no transcript_id: the ID as it is"
        );
        assert_eq!(
            named("id-LOC90010"),
            want("90010", "LOC90010", "IG_V_gene"),
            "a gene segment"
        );
        assert_eq!(transcripts.len(), 6);

        // A biotype the record states wins over every guess, and `transcript_type` (GENCODE's
        // spelling) counts as stating one.
        let stated = "##gff-version 3
chr1\tensembl\tgene\t1000\t1239\t.\t+\t.\tID=gene:ENSG00000001;Name=G;biotype=protein_coding
chr1\tensembl\tmRNA\t1000\t1239\t.\t+\t.\tID=transcript:ENST00000001;Parent=gene:ENSG00000001;transcript_id=ENST00000001;transcript_type=nonsense_mediated_decay
chr1\tensembl\texon\t1000\t1239\t.\t+\t.\tID=exon:E1;Parent=transcript:ENST00000001;rank=1";
        let transcripts = parse_gff3(stated.as_bytes()).unwrap();
        assert_eq!(&*transcripts[0].stable_id, "ENST00000001");
        assert_eq!(&*transcripts[0].gene.stable_id, "ENSG00000001");
        assert_eq!(&*transcripts[0].biotype, "nonsense_mediated_decay");
    }

    /// RefSeq writes a transcript that aligns in two places twice, `rna-NM_x.1` and `rna-NM_x.1-2`,
    /// and both are `transcript_id=NM_x.1`. Ensembl VEP names both so, and the cache has to be
    /// written in one order whichever the hash map yields them in.
    #[test]
    fn two_records_with_one_transcript_id_are_two_transcripts_in_a_fixed_order() {
        let gff = "##gff-version 3
chr1\tBestRefSeq\tgene\t1000\t1239\t.\t+\t.\tID=gene-A;Dbxref=GeneID:1;Name=A;gene_biotype=lncRNA
chr1\tBestRefSeq\tlnc_RNA\t1000\t1239\t.\t+\t.\tID=rna-NR_1.1;Parent=gene-A;Dbxref=GeneID:1;gbkey=ncRNA;transcript_id=NR_1.1
chr1\tBestRefSeq\texon\t1000\t1239\t.\t+\t.\tID=exon-NR_1.1-1;Parent=rna-NR_1.1;transcript_id=NR_1.1
chr1\tBestRefSeq\tgene\t5000\t5239\t.\t+\t.\tID=gene-A-2;Dbxref=GeneID:1;Name=A;gene_biotype=lncRNA
chr1\tBestRefSeq\tlnc_RNA\t5000\t5239\t.\t+\t.\tID=rna-NR_1.1-2;Parent=gene-A-2;Dbxref=GeneID:1;gbkey=ncRNA;transcript_id=NR_1.1
chr1\tBestRefSeq\texon\t5000\t5239\t.\t+\t.\tID=exon-NR_1.1-2-1;Parent=rna-NR_1.1-2;transcript_id=NR_1.1";
        let first: Vec<(String, u64)> = parse_gff3(gff.as_bytes())
            .unwrap()
            .iter()
            .map(|t| (t.stable_id.to_string(), t.start))
            .collect();
        assert_eq!(
            first,
            vec![("NR_1.1".to_string(), 1000), ("NR_1.1".to_string(), 5000)]
        );
        for _ in 0..20 {
            let again: Vec<(String, u64)> = parse_gff3(gff.as_bytes())
                .unwrap()
                .iter()
                .map(|t| (t.stable_id.to_string(), t.start))
                .collect();
            assert_eq!(
                again, first,
                "the order of equal names changed between parses"
            );
        }
    }

    /// The protein version is the CDS records' `version`, as Ensembl VEP 116 reads it, and 1
    /// where they have none or say 0. It was 1 for every protein whatever they said.
    #[test]
    fn the_protein_version_is_the_cds_records_version() {
        let version = |cds_attributes: &str| {
            let gff = format!(
                "##gff-version 3\n\
chr1\tensembl\tgene\t1000\t2000\t.\t+\t.\tID=gene:ENSG00000001;Name=G;biotype=protein_coding\n\
chr1\tensembl\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:ENST00000001;Parent=gene:ENSG00000001;biotype=protein_coding;version=2\n\
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:E1;Parent=transcript:ENST00000001;rank=1\n\
chr1\tensembl\tCDS\t1000\t1998\t.\t+\t0\tID=CDS:ENSP00000001;Parent=transcript:ENST00000001;protein_id=ENSP00000001{cds_attributes}"
            );
            parse_gff3(gff.as_bytes()).unwrap()[0].protein_version
        };
        assert_eq!(version(";version=3"), Some(3));
        assert_eq!(version(";version=12"), Some(12));
        assert_eq!(version(";version=1"), Some(1));
        assert_eq!(version(""), Some(1), "a record with no version is 1");
        assert_eq!(version(";version=0"), Some(1), "VEP's `|| 1` reads 0 as 1");

        // and a transcript with no CDS has no protein, so no version
        let non_coding = "##gff-version 3
chr1\tensembl\tgene\t1000\t2000\t.\t+\t.\tID=gene:ENSG00000002;Name=N;biotype=lncRNA
chr1\tensembl\tlnc_RNA\t1000\t2000\t.\t+\t.\tID=transcript:ENST00000002;Parent=gene:ENSG00000002;biotype=lncRNA
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:E2;Parent=transcript:ENST00000002;rank=1";
        assert_eq!(
            parse_gff3(non_coding.as_bytes()).unwrap()[0].protein_version,
            None
        );
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
    fn test_parse_gff3_mane_and_metadata() {
        let gff = "##gff-version 3
chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene:ENSG00000001;Name=TESTGENE;biotype=protein_coding
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST00000001;Parent=gene:ENSG00000001;biotype=protein_coding;tag=gencode_basic,gencode_primary,Ensembl_canonical,MANE_Select;ccdsid=CCDS30547.2;transcript_support_level=1 (assigned to previous version 2);version=7
chr1\tensembl\texon\t1000\t5000\t.\t+\t.\tID=exon:ENSE00000001;Parent=transcript:ENST00000001;rank=1
chr1\tensembl\tgene\t6000\t9000\t.\t-\t.\tID=gene:ENSG00000002;Name=TESTGENE2;biotype=protein_coding
chr1\tensembl\tmRNA\t6000\t9000\t.\t-\t.\tID=transcript:ENST00000002;Parent=gene:ENSG00000002;biotype=protein_coding;tag=gencode_basic,MANE_Plus_Clinical;transcript_support_level=2;version=3
chr1\tensembl\texon\t6000\t9000\t.\t-\t.\tID=exon:ENSE00000002;Parent=transcript:ENST00000002;rank=1";

        let transcripts = parse_gff3(gff.as_bytes()).unwrap();
        assert_eq!(transcripts.len(), 2);

        // MANE Select transcript
        let ms = transcripts
            .iter()
            .find(|t| &*t.stable_id == "ENST00000001")
            .unwrap();
        assert!(ms.canonical);
        assert_eq!(ms.mane_select.as_deref(), Some("ENST00000001.7"));
        assert!(ms.mane_plus_clinical.is_none());
        assert!(ms.gencode_primary);
        assert_eq!(ms.ccds.as_deref(), Some("CCDS30547.2"));
        assert_eq!(ms.tsl, Some(1));

        // MANE Plus Clinical transcript
        let mpc = transcripts
            .iter()
            .find(|t| &*t.stable_id == "ENST00000002")
            .unwrap();
        assert!(!mpc.canonical);
        assert!(mpc.mane_select.is_none());
        assert_eq!(mpc.mane_plus_clinical.as_deref(), Some("ENST00000002.3"));
        assert!(!mpc.gencode_primary);
        assert_eq!(mpc.tsl, Some(2));
        assert!(mpc.ccds.is_none());
    }

    #[test]
    fn test_url_decode() {
        assert_eq!(url_decode("hello%20world"), "hello world");
        assert_eq!(url_decode("100%25"), "100%");
        assert_eq!(url_decode("normal"), "normal");
    }

    #[test]
    fn feature_order_does_not_change_the_transcript_set() {
        // Issue #87 supposed that `sort -k1,1 -k4,4n -k5,5n` broke the parser by
        // placing children before their parents, dropping ~2% of transcripts.
        // It does not: parsing is two-pass, collecting all features before
        // resolving any parent, so the file may arrive in any order. Pinning
        // that here, because the ordering hypothesis is the natural one to reach
        // for again and the real cause was elsewhere (a region-restricted tabix
        // load being persisted as a whole-file cache).
        let lines: Vec<&str> = sample_gff3().lines().filter(|l| !l.is_empty()).collect();

        let parent_first = lines.join("\n");
        let child_first = {
            let mut v = lines.clone();
            v.reverse();
            v.join("\n")
        };
        // And the coordinate sort a user would actually apply.
        let coordinate_sorted = {
            let mut v: Vec<&str> = lines
                .iter()
                .filter(|l| !l.starts_with('#'))
                .copied()
                .collect();
            v.sort_by_key(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                (
                    f[0].to_string(),
                    f[3].parse::<u64>().unwrap_or(0),
                    f[4].parse::<u64>().unwrap_or(0),
                )
            });
            v.join("\n")
        };

        let fingerprint = |text: &str| {
            let mut trs = parse_gff3(text.as_bytes()).unwrap();
            trs.sort_by(|a, b| a.stable_id.cmp(&b.stable_id));
            trs.iter()
                .map(|t| {
                    let mut exons: Vec<(u64, u64)> =
                        t.exons.iter().map(|e| (e.start, e.end)).collect();
                    exons.sort_unstable();
                    format!(
                        "{}|{}|{}|{:?}|{:?}|{:?}|{exons:?}",
                        t.stable_id,
                        t.start,
                        t.end,
                        t.strand,
                        t.cdna_coding_start,
                        t.cdna_coding_end
                    )
                })
                .collect::<Vec<_>>()
        };

        let expected = fingerprint(&parent_first);
        assert!(!expected.is_empty(), "fixture parsed to nothing");
        assert_eq!(
            fingerprint(&child_first),
            expected,
            "reversing the file changed the transcript set"
        );
        assert_eq!(
            fingerprint(&coordinate_sorted),
            expected,
            "coordinate-sorting the file changed the transcript set"
        );
    }

    #[test]
    fn test_gff3_skips_malformed_coordinate_lines() {
        // Two lines: one valid mRNA, one with a non-numeric start. The malformed
        // line used to silently parse as start=0, end=0 — colliding with every
        // overlap query. It should now be skipped (with a logged warning),
        // leaving the valid line intact.
        let gff = "##gff-version 3
chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene:G1;Name=T;biotype=protein_coding
chr1\tensembl\tmRNA\tNOT_A_NUMBER\t5000\t.\t+\t.\tID=transcript:TX1;Parent=gene:G1;biotype=protein_coding
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:TX2;Parent=gene:G1;biotype=protein_coding;tag=Ensembl_canonical
chr1\tensembl\texon\t1000\t1200\t.\t+\t.\tID=exon:E1;Parent=transcript:TX2;rank=1
chr1\tensembl\tCDS\t1050\t1200\t.\t+\t0\tID=CDS:P1;Parent=transcript:TX2";
        let transcripts = parse_gff3(gff.as_bytes()).unwrap();
        // Only TX2 should appear; TX1 had a bogus start and was skipped.
        assert_eq!(transcripts.len(), 1);
        assert_eq!(&*transcripts[0].stable_id, "TX2");
        assert_eq!(transcripts[0].start, 1000);
    }

    #[test]
    fn test_gff3_skips_zero_and_inverted_coordinate_lines() {
        // Both lines parse fine as u64 (unlike the non-numeric case above)
        // but are invalid GFF3: start=0 is out of the 1-based coordinate
        // system, and start > end is inverted. Left unchecked, either would
        // underflow the `- 1` conversions to 0-based offsets done downstream
        // (e.g. transcript exon/CDS math in fastvep-genome).
        let gff = "##gff-version 3
chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene:G1;Name=T;biotype=protein_coding
chr1\tensembl\tmRNA\t0\t5000\t.\t+\t.\tID=transcript:TX1;Parent=gene:G1;biotype=protein_coding
chr1\tensembl\tmRNA\t5000\t1000\t.\t+\t.\tID=transcript:TX2;Parent=gene:G1;biotype=protein_coding
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:TX3;Parent=gene:G1;biotype=protein_coding;tag=Ensembl_canonical
chr1\tensembl\texon\t1000\t1200\t.\t+\t.\tID=exon:E1;Parent=transcript:TX3;rank=1
chr1\tensembl\tCDS\t1050\t1200\t.\t+\t0\tID=CDS:P1;Parent=transcript:TX3";
        let transcripts = parse_gff3(gff.as_bytes()).unwrap();
        // Only TX3 should appear; TX1 (start=0) and TX2 (start>end) were skipped.
        assert_eq!(transcripts.len(), 1);
        assert_eq!(&*transcripts[0].stable_id, "TX3");
        assert_eq!(transcripts[0].start, 1000);
    }

    #[test]
    fn test_parse_gff3_with_source_tags_every_transcript() {
        // The merged-cache path relies on per-file source tagging — without
        // it, transcripts from Ensembl and RefSeq runs are indistinguishable
        // in the SOURCE column.
        let transcripts = parse_gff3_with_source(sample_gff3().as_bytes(), "RefSeq").unwrap();
        assert!(!transcripts.is_empty());
        for tr in &transcripts {
            assert_eq!(tr.source.as_deref(), Some("RefSeq"));
        }
    }

    #[test]
    fn test_parse_gff3_default_source_is_gff3() {
        // Bare `parse_gff3()` must keep its historical "GFF3" label for
        // backwards compat with any caller (tests, web upload, etc.) that
        // doesn't thread a source through.
        let transcripts = parse_gff3(sample_gff3().as_bytes()).unwrap();
        assert!(!transcripts.is_empty());
        assert_eq!(transcripts[0].source.as_deref(), Some("GFF3"));
    }

    #[test]
    fn test_fasta_rejects_empty_header() {
        // An empty `>` line would silently make every following base part of an
        // unnamed sequence that no downstream query can fetch.
        let bad = ">\nACGT\n";
        let res = crate::fasta::FastaReader::from_reader(bad.as_bytes());
        assert!(res.is_err(), "FASTA with empty header should fail to parse");
    }
}
