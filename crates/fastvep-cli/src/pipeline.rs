use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue};
use fastvep_cache::fasta::FastaReader;
use fastvep_cache::gff::parse_gff3;
use fastvep_cache::info::CacheInfo;
use fastvep_cache::providers::{
    FastaSequenceProvider, IndexedTranscriptProvider, MatchedVariant, SequenceProvider,
    TabixVariationProvider, TranscriptProvider, VariationProvider,
};
use fastvep_consequence::ConsequencePredictor;
use fastvep_core::{Allele, Consequence};
use fastvep_hgvs;
use fastvep_io::output;
use fastvep_io::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature};
use fastvep_io::vcf::VcfParser;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

const BATCH_SIZE: usize = 1024;

fn open_vcf_input_reader(input: &str) -> Result<Box<dyn io::Read>> {
    let reader: Box<dyn io::Read> = if input == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(
            File::open(input)
                .with_context(|| format!("Opening input file: {}", input))?,
        )
    };

    wrap_maybe_gzip_reader(reader, input)
}

fn wrap_maybe_gzip_reader(mut reader: Box<dyn io::Read>, source: &str) -> Result<Box<dyn io::Read>> {
    let mut prefix = [0u8; 2];
    let bytes_read = reader.read(&mut prefix)?;
    let looks_like_gzip = bytes_read == 2 && prefix == [0x1f, 0x8b];

    let replay = io::Cursor::new(prefix[..bytes_read].to_vec()).chain(reader);
    if looks_like_gzip || (source != "-" && source.ends_with(".gz")) {
        Ok(Box::new(MultiGzDecoder::new(replay)))
    } else {
        Ok(Box::new(replay))
    }
}

pub struct AnnotateConfig {
    pub input: String,
    pub output: String,
    pub gff3: Option<String>,
    pub fasta: Option<String>,
    pub output_format: String,
    pub pick: bool,
    pub hgvs: bool,
    pub distance: u64,
    pub cache_dir: Option<String>,
    pub transcript_cache: Option<String>,
    /// Directory containing supplementary annotation files (.osa, .osi, .oga).
    pub sa_dir: Option<String>,
}

pub fn run_annotate(config: AnnotateConfig) -> Result<()> {
    // Extract the GFF3 source name (filename) for the SOURCE field
    let gff3_source: Option<String> = config.gff3.as_ref().map(|p| {
        Path::new(p).file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| p.clone())
    });

    // Load transcript models: try binary cache first, fall back to GFF3
    let cache_path = config.transcript_cache.as_ref().map(|p| Path::new(p).to_path_buf())
        .or_else(|| config.gff3.as_ref().map(|p| fastvep_cache::transcript_cache::default_cache_path(Path::new(p))));

    let mut transcripts = 'load: {
        // Try loading from binary cache
        if let Some(ref cp) = cache_path {
            if cp.exists() {
                let is_fresh = config.gff3.as_ref()
                    .map(|gff| fastvep_cache::transcript_cache::cache_is_fresh(cp, Path::new(gff)))
                    .unwrap_or(true);
                if is_fresh {
                    match fastvep_cache::transcript_cache::load_cache(cp) {
                        Ok(trs) => {
                            eprintln!("Loaded {} transcripts from cache {}", trs.len(), cp.display());
                            break 'load trs;
                        }
                        Err(e) => {
                            eprintln!("Warning: cache load failed ({}), falling back to GFF3", e);
                        }
                    }
                } else {
                    eprintln!("Cache is stale, rebuilding from GFF3");
                }
            }
        }

        // Fall back to GFF3 parsing
        if let Some(ref gff3_path) = config.gff3 {
            let gff_path = Path::new(gff3_path);
            let tbi_path = format!("{}.tbi", gff3_path);

            // Use indexed loading if .gff3.gz + .tbi available
            if gff3_path.ends_with(".gz") && Path::new(&tbi_path).exists() {
                // Pre-scan VCF to collect regions
                let regions = prescan_vcf_regions(&config.input, config.distance)?;
                eprintln!("Pre-scanned {} variant regions from {}", regions.len(), config.input);
                let trs = fastvep_cache::gff::parse_gff3_indexed(gff_path, &regions)?;
                eprintln!("Loaded {} transcripts from indexed {}", trs.len(), gff3_path);
                trs
            } else {
                let gff_file = File::open(gff3_path)
                    .with_context(|| format!("Opening GFF3 file: {}", gff3_path))?;
                let trs = parse_gff3(gff_file)?;
                eprintln!("Loaded {} transcripts from {}", trs.len(), gff3_path);
                trs
            }
        } else {
            eprintln!("Warning: No GFF3 file provided. Only intergenic variants will be annotated.");
            Vec::new()
        }
    };

    // Load FASTA reference (prefer mmap with .fai index, fall back to in-memory)
    let seq_provider: Option<Box<dyn SequenceProvider>> = if let Some(ref fasta_path) = config.fasta {
        let fai_path = format!("{}.fai", fasta_path);
        if Path::new(&fai_path).exists() {
            let reader = fastvep_cache::fasta::MmapFastaReader::open(Path::new(fasta_path))?;
            eprintln!("Memory-mapped reference FASTA from {} (using .fai index)", fasta_path);
            Some(Box::new(fastvep_cache::providers::MmapFastaSequenceProvider::new(reader)))
        } else {
            let fasta_file = File::open(fasta_path)
                .with_context(|| format!("Opening FASTA file: {}", fasta_path))?;
            let reader = FastaReader::from_reader(fasta_file)?;
            eprintln!("Loaded reference FASTA from {}", fasta_path);
            Some(Box::new(FastaSequenceProvider::new(reader)))
        }
    } else {
        None
    };

    // Build sequences for coding transcripts from FASTA (skip if loaded from cache with sequences)
    let needs_seq_build = transcripts.iter().any(|t| t.is_coding() && t.spliced_seq.is_none());
    if needs_seq_build {
        if let Some(ref sp) = seq_provider {
            let built = AtomicUsize::new(0);
            transcripts.par_iter_mut().for_each(|tr| {
                if tr.is_coding() && tr.spliced_seq.is_none() {
                    if let Err(e) = tr.build_sequences(|chrom, start, end| {
                        sp.fetch_sequence(chrom, start, end)
                            .map_err(|e| e.to_string())
                    }) {
                        eprintln!("Warning: could not build sequences for {}: {}", tr.stable_id, e);
                    } else {
                        built.fetch_add(1, Ordering::Relaxed);
                    }
                }
            });
            eprintln!("Built sequences for {} coding transcripts", built.load(Ordering::Relaxed));
        }
    }

    // Save cache after sequence build (only if sequences were built or cache doesn't exist)
    if needs_seq_build {
        if let Some(ref cp) = cache_path {
            if config.gff3.is_some() {
                if let Err(e) = fastvep_cache::transcript_cache::save_cache(&transcripts, cp) {
                    eprintln!("Warning: could not save cache: {}", e);
                } else {
                    eprintln!("Saved transcript cache to {}", cp.display());
                }
            }
        }
    }

    let transcript_provider = IndexedTranscriptProvider::new(transcripts);

    // Initialize variation provider from VEP cache if provided
    let var_provider: Option<TabixVariationProvider> = if let Some(ref dir) = config.cache_dir {
        let info_path = Path::new(dir).join("info.txt");
        let cache_info = CacheInfo::from_file(&info_path)
            .with_context(|| format!("Reading cache info: {}", info_path.display()))?;
        eprintln!(
            "Loaded VEP cache info: species={}, assembly={}, {} variation columns",
            cache_info.species, cache_info.assembly, cache_info.variation_cols.len()
        );
        Some(TabixVariationProvider::new(Path::new(dir), &cache_info)?)
    } else {
        None
    };

    // Load supplementary annotation providers from --sa-dir
    // Shared load_sa_providers returns Mutex-wrapped providers; unwrap for batch pipeline
    // (batch pipeline processes variants sequentially per chunk, no concurrent access).
    let sa_providers: Vec<Box<dyn AnnotationProvider>> = if let Some(ref dir) = config.sa_dir {
        load_sa_providers(Path::new(dir))?
            .into_iter()
            .map(|m| m.into_inner().unwrap())
            .collect()
    } else {
        Vec::new()
    };

    // Create consequence predictor
    let predictor = ConsequencePredictor::new(config.distance, config.distance);

    // Open input VCF (supports plain text or gzipped VCF)
    let input_reader = open_vcf_input_reader(&config.input)?;
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
                "## fastVEP output"
            )?;
            writeln!(
                writer,
                "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS"
            )?;
        }
        "json" => {
            writeln!(writer, "[")?;
        }
        _ => {}
    }

    // Process variants in batches for parallel annotation
    let mut count = 0u64;
    let mut first_json = true;

    loop {
        // Phase 1: Read a batch of variants (sequential - VCF parser is not Sync)
        let mut batch: Vec<(VariationFeature, HashMap<String, Vec<MatchedVariant>>)> = Vec::with_capacity(BATCH_SIZE);
        for _ in 0..BATCH_SIZE {
            match vcf_parser.next_variant()? {
                Some(mut vf) => {
                    // Variation lookup (sequential - TabixVariationProvider is not Sync)
                    let matched_by_allele: HashMap<String, Vec<MatchedVariant>> =
                        if let Some(ref vp) = var_provider {
                            let mut by_allele = HashMap::new();
                            for alt in &vf.alt_alleles {
                                let alt_str = alt.to_string();
                                let ref_str = vf.ref_allele.to_string();
                                let matches = vp.get_matched_variants(
                                    &vf.position.chromosome,
                                    vf.position.start,
                                    vf.position.end,
                                    &ref_str,
                                    &alt_str,
                                ).unwrap_or_default();
                                if !matches.is_empty() {
                                    by_allele.insert(alt_str, matches);
                                }
                            }
                            by_allele
                        } else {
                            HashMap::new()
                        };

                    // Populate existing_variants on the VF for output access
                    for matches in matched_by_allele.values() {
                        for m in matches {
                            if !vf.existing_variants.iter().any(|kv| kv.name == m.name) {
                                vf.existing_variants.push(fastvep_io::variant::KnownVariant {
                                    name: m.name.clone(),
                                    allele_string: None,
                                    minor_allele: m.minor_allele.clone(),
                                    minor_allele_freq: m.minor_allele_freq,
                                    clinical_significance: m.clin_sig.clone(),
                                    somatic: m.somatic,
                                    phenotype_or_disease: m.phenotype_or_disease,
                                    pubmed: m.pubmed.clone(),
                                    frequencies: m.frequencies.clone(),
                                });
                            }
                        }
                    }

                    batch.push((vf, matched_by_allele));
                }
                None => break,
            }
        }

        if batch.is_empty() {
            break;
        }

        // Phase 1.5: Preload SA providers for this batch (sequential)
        if !sa_providers.is_empty() && !batch.is_empty() {
            // Collect all positions in this batch, grouped by chromosome
            let mut chrom_positions: HashMap<&str, Vec<u64>> = HashMap::new();
            for (vf, _) in &batch {
                chrom_positions
                    .entry(&vf.position.chromosome)
                    .or_default()
                    .push(vf.position.start);
            }
            for sa in &sa_providers {
                for (chrom, positions) in &chrom_positions {
                    let _ = sa.preload(chrom, positions);
                }
            }
        }

        // Phase 2: Annotate batch in parallel (transcript lookup + consequence prediction + HGVS)
        batch.par_iter_mut().for_each(|(vf, matched_by_allele)| {
            let chrom = &vf.position.chromosome;
            let query_start = if vf.position.start > config.distance {
                vf.position.start - config.distance
            } else {
                1
            };
            let query_end = vf.position.end + config.distance;
            let overlapping = transcript_provider.get_transcripts(chrom, query_start, query_end)
                .unwrap_or_default();

        if overlapping.is_empty() {
            // Intergenic
            annotate_intergenic(vf);
            // Populate existing_variation on intergenic annotations too
            for tv in &mut vf.transcript_variations {
                for aa in &mut tv.allele_annotations {
                    if let Some(matches) = matched_by_allele.get(&aa.allele.to_string()) {
                        aa.existing_variation = matches.iter().map(|m| m.name.clone()).collect();
                    }
                }
            }
        } else {
            // Get reference sequence if available
            let ref_seq = seq_provider.as_ref().and_then(|sp| {
                sp.fetch_sequence(chrom, query_start, query_end).ok()
            });

            // Run consequence prediction — dispatch SVs to SV predictor
            let transcript_consequences = if vf.variant_type.is_structural() {
                fastvep_consequence::sv_predictor::predict_sv_consequences(
                    chrom,
                    vf.position.start,
                    vf.position.end,
                    vf.variant_type,
                    &vf.alt_alleles,
                    &overlapping,
                    config.distance,
                    config.distance,
                )
            } else {
                let result = predictor.predict(
                    &vf.position,
                    &vf.ref_allele,
                    &vf.alt_alleles,
                    &overlapping,
                    ref_seq.as_deref(),
                );
                result.transcript_consequences
            };

            // Convert prediction results to VariationFeature annotations
            for tc in &transcript_consequences {
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
                            hgvs_offset: None,
                            existing_variation: matched_by_allele
                                .get(&ac.allele.to_string())
                                .map(|matches| matches.iter().map(|m| m.name.clone()).collect())
                                .unwrap_or_default(),
                            sift: None,
                            polyphen: None,
                            supplementary: Vec::new(),
                        };

                        // Generate HGVS if requested
                        if config.hgvs {
                            ann.hgvsg = Some(fastvep_hgvs::hgvsg(
                                chrom,
                                vf.position.start,
                                vf.position.end,
                                &vf.ref_allele,
                                &ac.allele,
                            ));

                            if let Some(tr) = transcript {
                                // Build versioned IDs for HGVS notation
                                let versioned_tid = match tr.version {
                                    Some(v) => format!("{}.{}", tc.transcript_id, v),
                                    None => tc.transcript_id.to_string(),
                                };

                                // Determine alleles for HGVS - complement for minus strand
                                let (hgvs_ref, hgvs_alt) = if tr.strand == fastvep_core::Strand::Reverse {
                                    (complement_allele(&vf.ref_allele), complement_allele(&ac.allele))
                                } else {
                                    (vf.ref_allele.clone(), ac.allele.clone())
                                };

                                if let Some(coding_start) = tr.cdna_coding_start {
                                    if let (Some(cs), Some(ce)) = (ac.cdna_start, ac.cdna_end) {
                                        // Normalize cDNA positions (minus-strand can reverse order)
                                        let (cs, ce) = (cs.min(ce), cs.max(ce));
                                        // Exonic variant: standard HGVSc with 3' shifting
                                        ann.hgvsc = fastvep_hgvs::hgvsc_with_seq(
                                            &versioned_tid,
                                            cs, ce,
                                            &hgvs_ref,
                                            &hgvs_alt,
                                            coding_start,
                                            tr.cdna_coding_end,
                                            tr.spliced_seq.as_deref(),
                                            tr.codon_table_start_phase,
                                        );
                                    } else if ac.intron.is_some() {
                                        // Intronic variant: offset notation
                                        // Note: intronic HGVS uses original coding_start (no phase adjustment)
                                        // Apply HGVS 3' normalization for intronic indels
                                        let (shifted_start, shifted_end) = if let Some(ref sp) = seq_provider {
                                            let is_indel = matches!((&hgvs_ref, &hgvs_alt),
                                                (Allele::Sequence(_), Allele::Deletion) |
                                                (Allele::Deletion, Allele::Sequence(_)));
                                            if is_indel {
                                                if let Some((istart, iend)) = tr.intron_bounds_at(vf.position.start) {
                                                    // Use genomic-strand alleles for ref comparison
                                                    three_prime_shift_intronic(
                                                        &**sp as &dyn SequenceProvider, chrom,
                                                        vf.position.start, vf.position.end,
                                                        &vf.ref_allele, &ac.allele,
                                                        tr.strand, istart, iend,
                                                    )
                                                } else {
                                                    (vf.position.start, vf.position.end)
                                                }
                                            } else {
                                                (vf.position.start, vf.position.end)
                                            }
                                        } else {
                                            (vf.position.start, vf.position.end)
                                        };
                                        // For insertions, build the rotated insertion bases
                                        // after 3' shifting (bases rotate as position shifts)
                                        let shifted_hgvs_alt = if let (Allele::Deletion, Allele::Sequence(ins_bases)) = (&hgvs_ref, &hgvs_alt) {
                                            if shifted_start != vf.position.start && !ins_bases.is_empty() {
                                                // Calculate how many positions we shifted
                                                let shift_amount = if tr.strand == fastvep_core::Strand::Forward {
                                                    (shifted_start as i64 - vf.position.start as i64) as usize
                                                } else {
                                                    (vf.position.start as i64 - shifted_start as i64) as usize
                                                };
                                                // Rotate: for forward strand, each shift moves first base to end
                                                // For reverse strand, each shift moves last base to front
                                                let mut rotated = ins_bases.clone();
                                                let len = rotated.len();
                                                if len > 0 {
                                                    let effective_shift = shift_amount % len;
                                                    match tr.strand {
                                                        fastvep_core::Strand::Forward => {
                                                            rotated.rotate_left(effective_shift);
                                                        }
                                                        fastvep_core::Strand::Reverse => {
                                                            rotated.rotate_right(effective_shift);
                                                        }
                                                    }
                                                }
                                                Allele::Sequence(rotated)
                                            } else {
                                                hgvs_alt.clone()
                                            }
                                        } else {
                                            hgvs_alt.clone()
                                        };

                                        // For insertions, use position before insertion
                                        // for the primary HGVS coordinate (ins is BETWEEN two bases).
                                        // On reverse strand, the insertion is between P and P+1 in
                                        // genomic coords, but P+1 is 5' in transcript order, so we
                                        // use P+1 as the HGVS start coordinate.
                                        let is_insertion = matches!((&hgvs_ref, &shifted_hgvs_alt), (Allele::Deletion, Allele::Sequence(_)));
                                        let hgvs_pos = if is_insertion {
                                            if tr.strand == fastvep_core::Strand::Reverse {
                                                shifted_end + 1
                                            } else {
                                                shifted_end // base before insertion
                                            }
                                        } else {
                                            shifted_start
                                        };
                                        if let Some((cdna_pos, offset)) = tr.genomic_to_intronic_cdna(hgvs_pos) {
                                            // For multi-base variants, compute end position too
                                            let (end_cdna, end_offset) = if shifted_start != shifted_end && hgvs_pos == shifted_start {
                                                tr.genomic_to_intronic_cdna(shifted_end)
                                                    .map(|(c, o)| (Some(c), Some(o)))
                                                    .unwrap_or((None, None))
                                            } else {
                                                (None, None)
                                            };
                                            let mut hgvsc = fastvep_hgvs::hgvsc_intronic_range(
                                                &versioned_tid,
                                                cdna_pos,
                                                offset,
                                                end_cdna,
                                                end_offset,
                                                &hgvs_ref,
                                                &shifted_hgvs_alt,
                                                coding_start,
                                                tr.cdna_coding_end,
                                            );
                                            // For intronic insertions, check if it's a dup.
                                            if let (Some(ref h), Allele::Deletion, Allele::Sequence(_)) =
                                                (&hgvsc, &hgvs_ref, &hgvs_alt)
                                            {
                                                if h.contains("ins") {
                                                    let orig_ins = match &ac.allele {
                                                        Allele::Sequence(b) => b.clone(),
                                                        _ => vec![],
                                                    };
                                                    if !orig_ins.is_empty() {
                                                        if let Some(ref sp) = seq_provider {
                                                            let ins_len = orig_ins.len() as u64;
                                                            // Check dup_before: base(s) before insertion match
                                                            let check_end = vf.position.end;
                                                            let check_start = check_end.saturating_sub(ins_len - 1);
                                                            let dup_before = if let Ok(ref_seq) = sp.fetch_sequence_slice(chrom, check_start, check_end) {
                                                                ref_seq.len() == orig_ins.len()
                                                                    && ref_seq.iter().zip(orig_ins.iter())
                                                                        .all(|(a, b)| a.eq_ignore_ascii_case(b))
                                                            } else { false };
                                                            // Check dup_after: base(s) after insertion match
                                                            let dup_after = if !dup_before {
                                                                let cs = vf.position.start;
                                                                let ce = cs + ins_len - 1;
                                                                if let Ok(ref_seq) = sp.fetch_sequence_slice(chrom, cs, ce) {
                                                                    ref_seq.len() == orig_ins.len()
                                                                        && ref_seq.iter().zip(orig_ins.iter())
                                                                            .all(|(a, b)| a.eq_ignore_ascii_case(b))
                                                                } else { false }
                                                            } else { false };
                                                            if dup_before || dup_after {
                                                                // For dups, determine the dup base position and 3' shift it
                                                                let dup_base_pos = if dup_before {
                                                                    // Dup base is before insertion: position.end
                                                                    vf.position.end
                                                                } else {
                                                                    // Dup base is after insertion: position.start
                                                                    vf.position.start
                                                                };
                                                                // 3' shift the dup position within the intron
                                                                let shifted_dup = if let Some((istart, iend)) = tr.intron_bounds_at(dup_base_pos) {
                                                                    let (sd, _) = three_prime_shift_intronic(
                                                                        &**sp as &dyn SequenceProvider, chrom,
                                                                        dup_base_pos, dup_base_pos,
                                                                        &Allele::Sequence(orig_ins.clone()), &Allele::Deletion,
                                                                        tr.strand, istart, iend,
                                                                    );
                                                                    sd
                                                                } else {
                                                                    dup_base_pos
                                                                };
                                                                // Use shifted_dup (start of dup region) for offset computation
                                                                if let Some((dup_cdna, dup_offset)) = tr.genomic_to_intronic_cdna(shifted_dup) {
                                                                    hgvsc = convert_ins_to_dup(h, dup_offset, ins_len, dup_cdna, coding_start, tr.cdna_coding_end);
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            ann.hgvsc = hgvsc;
                                        }
                                    }
                                } else {
                                    // Non-coding transcript: use n. notation
                                    if let (Some(cs), Some(ce)) = (ac.cdna_start, ac.cdna_end) {
                                        ann.hgvsc = fastvep_hgvs::hgvsc_noncoding(
                                            &versioned_tid,
                                            cs, ce,
                                            &hgvs_ref,
                                            &hgvs_alt,
                                        );
                                    } else if ac.intron.is_some() {
                                        // Apply 3' normalization for non-coding intronic indels
                                        let (nc_shifted_start, nc_shifted_end) = if let Some(ref sp) = seq_provider {
                                            let is_indel = matches!((&hgvs_ref, &hgvs_alt),
                                                (Allele::Sequence(_), Allele::Deletion) |
                                                (Allele::Deletion, Allele::Sequence(_)));
                                            if is_indel {
                                                if let Some((istart, iend)) = tr.intron_bounds_at(vf.position.start) {
                                                    three_prime_shift_intronic(
                                                        &**sp as &dyn SequenceProvider, chrom,
                                                        vf.position.start, vf.position.end,
                                                        &vf.ref_allele, &ac.allele,
                                                        tr.strand, istart, iend,
                                                    )
                                                } else {
                                                    (vf.position.start, vf.position.end)
                                                }
                                            } else {
                                                (vf.position.start, vf.position.end)
                                            }
                                        } else {
                                            (vf.position.start, vf.position.end)
                                        };

                                        // Rotate insertion bases for non-coding
                                        let nc_shifted_hgvs_alt = if let (Allele::Deletion, Allele::Sequence(ins_bases)) = (&hgvs_ref, &hgvs_alt) {
                                            if nc_shifted_start != vf.position.start && !ins_bases.is_empty() {
                                                let shift_amount = if tr.strand == fastvep_core::Strand::Forward {
                                                    (nc_shifted_start as i64 - vf.position.start as i64) as usize
                                                } else {
                                                    (vf.position.start as i64 - nc_shifted_start as i64) as usize
                                                };
                                                let mut rotated = ins_bases.clone();
                                                let len = rotated.len();
                                                if len > 0 {
                                                    let effective_shift = shift_amount % len;
                                                    match tr.strand {
                                                        fastvep_core::Strand::Forward => rotated.rotate_left(effective_shift),
                                                        fastvep_core::Strand::Reverse => rotated.rotate_right(effective_shift),
                                                    }
                                                }
                                                Allele::Sequence(rotated)
                                            } else {
                                                hgvs_alt.clone()
                                            }
                                        } else {
                                            hgvs_alt.clone()
                                        };

                                        if let Some((cdna_pos, offset)) = tr.genomic_to_intronic_cdna(nc_shifted_start) {
                                            let (end_cdna, end_offset) = if nc_shifted_start != nc_shifted_end {
                                                tr.genomic_to_intronic_cdna(nc_shifted_end)
                                                    .map(|(c, o)| (Some(c), Some(o)))
                                                    .unwrap_or((None, None))
                                            } else {
                                                (None, None)
                                            };
                                            let mut hgvsc = fastvep_hgvs::hgvsc_noncoding_intronic_range(
                                                &versioned_tid,
                                                cdna_pos,
                                                offset,
                                                end_cdna,
                                                end_offset,
                                                &hgvs_ref,
                                                &nc_shifted_hgvs_alt,
                                            );
                                            // Dup detection for non-coding intronic insertions
                                            if let (Some(ref h), Allele::Deletion, Allele::Sequence(_)) =
                                                (&hgvsc, &hgvs_ref, &hgvs_alt)
                                            {
                                                if h.contains("ins") {
                                                    let orig_ins = match &ac.allele {
                                                        Allele::Sequence(b) => b.clone(),
                                                        _ => vec![],
                                                    };
                                                    if !orig_ins.is_empty() {
                                                        if let Some(ref sp) = seq_provider {
                                                            let ins_len = orig_ins.len() as u64;
                                                            let check_end = vf.position.end;
                                                            let check_start = check_end.saturating_sub(ins_len - 1);
                                                            let dup_before = if let Ok(ref_seq) = sp.fetch_sequence_slice(chrom, check_start, check_end) {
                                                                ref_seq.len() == orig_ins.len()
                                                                    && ref_seq.iter().zip(orig_ins.iter())
                                                                        .all(|(a, b)| a.eq_ignore_ascii_case(b))
                                                            } else { false };
                                                            let dup_after = if !dup_before {
                                                                let cs = vf.position.start;
                                                                let ce = cs + ins_len - 1;
                                                                if let Ok(ref_seq) = sp.fetch_sequence_slice(chrom, cs, ce) {
                                                                    ref_seq.len() == orig_ins.len()
                                                                        && ref_seq.iter().zip(orig_ins.iter())
                                                                            .all(|(a, b)| a.eq_ignore_ascii_case(b))
                                                                } else { false }
                                                            } else { false };
                                                            if dup_before || dup_after {
                                                                let dup_base_pos = if dup_before { vf.position.end } else { vf.position.start };
                                                                let shifted_dup = if let Some((istart, iend)) = tr.intron_bounds_at(dup_base_pos) {
                                                                    let (sd, _) = three_prime_shift_intronic(
                                                                        &**sp as &dyn SequenceProvider, chrom,
                                                                        dup_base_pos, dup_base_pos,
                                                                        &Allele::Sequence(orig_ins.clone()), &Allele::Deletion,
                                                                        tr.strand, istart, iend,
                                                                    );
                                                                    sd
                                                                } else {
                                                                    dup_base_pos
                                                                };
                                                                if let Some((dup_cdna, dup_offset)) = tr.genomic_to_intronic_cdna(shifted_dup) {
                                                                    hgvsc = convert_ins_to_dup_noncoding(h, dup_offset, ins_len, dup_cdna);
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            ann.hgvsc = hgvsc;
                                        }
                                    }
                                }
                            }

                            if let (Some(ref aa), Some(ps)) = (&ac.amino_acids, ac.protein_start) {
                                if let Some(tr) = transcript {
                                    if let Some(ref pid) = tr.protein_id {
                                        let versioned_pid = match tr.protein_version {
                                            Some(v) => {
                                                let suffix = format!(".{}", v);
                                                if pid.ends_with(&suffix) {
                                                    pid.clone()
                                                } else {
                                                    format!("{}.{}", pid, v)
                                                }
                                            }
                                            None => pid.clone(),
                                        };
                                        let is_fs = ac.consequences.contains(&Consequence::FrameshiftVariant);

                                        if is_fs {
                                            // Frameshift: build alt sequence and scan for first changed AA + new stop
                                            // Use spliced_seq from CDS start onwards (includes 3'UTR for stop codon search)
                                            if let (Some(ref spliced), Some(coding_start), Some(cds_s)) =
                                                (&tr.spliced_seq, tr.cdna_coding_start, ac.cds_start)
                                            {
                                                // Extract from CDS start to end of spliced seq (includes 3'UTR)
                                                let coding_start_idx = (coding_start - 1) as usize;
                                                let ref_from_cds = &spliced.as_bytes()[coding_start_idx..];
                                                let cds_idx = (cds_s - 1) as usize;
                                                let mut alt_from_cds = ref_from_cds.to_vec();

                                                // Apply the indel to build the frameshifted sequence
                                                if ac.allele == Allele::Deletion {
                                                    let del_len = vf.ref_allele.len();
                                                    let end = (cds_idx + del_len).min(alt_from_cds.len());
                                                    alt_from_cds.drain(cds_idx..end);
                                                } else if let Allele::Sequence(ins_bases) = &ac.allele {
                                                    let mut bases = ins_bases.clone();
                                                    if tr.strand == fastvep_core::Strand::Reverse {
                                                        bases = bases.iter().map(|&b| match b {
                                                            b'A' => b'T', b'T' => b'A',
                                                            b'C' => b'G', b'G' => b'C',
                                                            o => o,
                                                        }).collect();
                                                    }
                                                    for (j, &b) in bases.iter().enumerate() {
                                                        if cds_idx + j <= alt_from_cds.len() {
                                                            alt_from_cds.insert(cds_idx + j, b);
                                                        }
                                                    }
                                                }

                                                let codon_start = cds_idx / 3;
                                                ann.hgvsp = fastvep_hgvs::hgvsp_frameshift(
                                                    &versioned_pid,
                                                    ref_from_cds,
                                                    &alt_from_cds,
                                                    codon_start,
                                                );
                                            }
                                        } else {
                                            let ref_aa_byte = aa.0.as_bytes().first().copied().unwrap_or(b'X');
                                            let alt_aa_byte = aa.1.as_bytes().first().copied().unwrap_or(b'X');
                                            ann.hgvsp = fastvep_hgvs::hgvsp(
                                                &versioned_pid, ps, ref_aa_byte, alt_aa_byte, false,
                                            );
                                        }
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
                        source: gff3_source.clone(),
                        protein_id: transcript.and_then(|t| t.protein_id.clone()),
                        mane_select: transcript.and_then(|t| t.mane_select.clone()),
                        mane_plus_clinical: transcript.and_then(|t| t.mane_plus_clinical.clone()),
                        tsl: transcript.and_then(|t| t.tsl),
                        appris: transcript.and_then(|t| t.appris.clone()),
                        ccds: transcript.and_then(|t| t.ccds.clone()),
                        gencode_primary: transcript.map(|t| t.gencode_primary).unwrap_or(false),
                        symbol_source: transcript.and_then(|t| t.gene.symbol_source.clone()),
                        hgnc_id: transcript.and_then(|t| t.gene.hgnc_id.clone()),
                        flags: transcript.map(|t| t.flags.clone()).unwrap_or_default(),
                    });
                }
            }

            // Supplementary annotation: query SA providers for each allele
            if !sa_providers.is_empty() {
                let chrom = &vf.position.chromosome;
                for tv in &mut vf.transcript_variations {
                    for aa in &mut tv.allele_annotations {
                        let alt_str = aa.allele.to_string();
                        let ref_str = vf.ref_allele.to_string();
                        for sa in &sa_providers {
                            if let Ok(Some(ann)) =
                                sa.annotate_position(chrom, vf.position.start, &ref_str, &alt_str)
                            {
                                let json_str = match ann {
                                    AnnotationValue::Json(j) => j,
                                    AnnotationValue::Positional(j) => j,
                                    AnnotationValue::Interval(v) => {
                                        format!("[{}]", v.join(","))
                                    }
                                };
                                aa.supplementary.push((
                                    sa.json_key().to_string(),
                                    json_str,
                                ));
                            }
                        }
                    }
                }
            }

            vf.compute_most_severe();
        }
        }); // end par_iter_mut

        // Phase 3: Write output sequentially (preserves VCF order)
        for (vf, _) in &batch {
            match config.output_format.as_str() {
                "vcf" => write_vcf_line(&mut writer, vf)?,
                "tab" => {
                    for line in output::format_tab_line(vf) {
                        writeln!(writer, "{}", line)?;
                    }
                }
                "json" => {
                    if !first_json {
                        writeln!(writer, ",")?;
                    }
                    first_json = false;
                    let json = output::format_json(vf);
                    write!(writer, "{}", serde_json::to_string_pretty(&json)?)?;
                }
                _ => {}
            }
        }

        count += batch.len() as u64;
    } // end batch loop

    // Close JSON array
    if config.output_format == "json" {
        writeln!(writer, "\n]")?;
    }

    writer.flush()?;
    eprintln!("Annotated {} variants", count);

    Ok(())
}

// Shared annotation utilities from fastvep-annotate (used by batch pipeline).
use fastvep_annotate::{
    annotate_intergenic, complement_allele, convert_ins_to_dup, convert_ins_to_dup_noncoding,
    load_sa_providers, three_prime_shift_intronic, zip_positions,
};



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


use serde_json;

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

    // Open output
    let mut writer: Box<dyn Write> = if output_path == "-" {
        Box::new(BufWriter::new(io::stdout()))
    } else {
        let f = File::create(output_path)
            .with_context(|| format!("Creating output: {}", output_path))?;
        Box::new(BufWriter::new(f))
    };

    // Parse CSQ header to get field names
    let mut csq_fields: Vec<String> = Vec::new();
    let mut kept = 0u64;
    let mut total = 0u64;

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
                        csq_fields.iter().take(5).cloned().collect::<Vec<_>>().join(", ")
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

    Ok(())
}

/// Build a binary transcript cache from GFF3 + optional FASTA.
pub fn run_cache_build(gff3_path: &str, fasta_path: Option<&str>, output_path: &str) -> Result<()> {
    let gff_file = File::open(gff3_path)
        .with_context(|| format!("Opening GFF3 file: {}", gff3_path))?;
    let mut transcripts = parse_gff3(gff_file)?;
    eprintln!("Loaded {} transcripts from {}", transcripts.len(), gff3_path);

    if let Some(fasta) = fasta_path {
        let fasta_file = File::open(fasta)
            .with_context(|| format!("Opening FASTA file: {}", fasta))?;
        let reader = FastaReader::from_reader(fasta_file)?;
        let sp = FastaSequenceProvider::new(reader);
        eprintln!("Loaded reference FASTA from {}", fasta);

        let mut built = 0usize;
        for tr in &mut transcripts {
            if tr.is_coding() {
                if let Err(e) = tr.build_sequences(|chrom, start, end| {
                    sp.fetch_sequence(chrom, start, end)
                        .map_err(|e| e.to_string())
                }) {
                    eprintln!("Warning: could not build sequences for {}: {}", tr.stable_id, e);
                } else {
                    built += 1;
                }
            }
        }
        eprintln!("Built sequences for {} coding transcripts", built);
    }

    fastvep_cache::transcript_cache::save_cache(&transcripts, Path::new(output_path))?;
    eprintln!("Saved transcript cache to {}", output_path);
    Ok(())
}

/// Quick VCF pre-scan to collect variant regions for indexed GFF3 loading.
/// Returns merged (chrom, start, end) regions expanded by the given distance.
fn prescan_vcf_regions(vcf_path: &str, distance: u64) -> Result<Vec<(String, u64, u64)>> {
    let input_reader = open_vcf_input_reader(vcf_path)
        .with_context(|| format!("Pre-scanning VCF: {}", vcf_path))?;
    let reader = io::BufReader::new(input_reader);

    let mut regions: HashMap<String, (u64, u64)> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() { continue; }
        let mut fields = line.split('\t');
        let chrom = match fields.next() { Some(c) => c.to_string(), None => continue };
        let pos: u64 = match fields.next().and_then(|p| p.parse().ok()) { Some(p) => p, None => continue };
        let start = pos.saturating_sub(distance);
        let end = pos + distance;

        let entry = regions.entry(chrom).or_insert((start, end));
        entry.0 = entry.0.min(start);
        entry.1 = entry.1.max(end);
    }

    Ok(regions.into_iter().map(|(chrom, (s, e))| (chrom, s, e)).collect())
}

// =============================================================================
// SA Build: Build supplementary annotation databases from source VCFs
// =============================================================================

/// Standard chromosome ordering for SA builds.
fn standard_chrom_map() -> (Vec<String>, std::collections::HashMap<String, u16>) {
    // Support both "chr1" and "1" naming conventions
    let chroms: Vec<String> = (1..=22)
        .map(|i| i.to_string())
        .chain(["X", "Y", "MT"].iter().map(|s| s.to_string()))
        .collect();
    let mut map: std::collections::HashMap<String, u16> = chroms
        .iter()
        .enumerate()
        .map(|(i, c)| (c.clone(), i as u16))
        .collect();
    // Also map "chr" prefixed names to the same indices
    for (i, c) in chroms.iter().enumerate() {
        map.insert(format!("chr{}", c), i as u16);
    }
    // Common aliases
    map.insert("chrM".to_string(), *map.get("MT").unwrap_or(&24));
    map.insert("M".to_string(), *map.get("MT").unwrap_or(&24));
    (chroms, map)
}

/// Build a supplementary annotation .osa file from a source VCF.
pub fn run_sa_build(source: &str, input: &str, output: &str, assembly: &str) -> Result<()> {
    use fastvep_sa::index::IndexHeader;
    use fastvep_sa::writer::SaWriter;

    let (chrom_list, chrom_map) = standard_chrom_map();

    let header = match source {
        "clinvar" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "clinvar".into(),
            name: "ClinVar".into(),
            version: "latest".into(),
            description: format!("ClinVar annotations for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: true,
            is_positional: false,
        },
        "gnomad" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "gnomad".into(),
            name: "gnomAD".into(),
            version: "latest".into(),
            description: format!("gnomAD population frequencies for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "dbsnp" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "dbsnp".into(),
            name: "dbSNP".into(),
            version: "latest".into(),
            description: format!("dbSNP RS IDs for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "phylop" | "gerp" | "dann" => {
            let json_key = match source {
                "phylop" => "phylopScore",
                "gerp" => "gerpScore",
                "dann" => "dannScore",
                _ => unreachable!(),
            };
            IndexHeader {
                schema_version: fastvep_sa::common::SCHEMA_VERSION,
                json_key: json_key.into(),
                name: source.to_uppercase().into(),
                version: "latest".into(),
                description: format!("{} conservation/prediction scores for {}", source, assembly),
                assembly: assembly.into(),
                match_by_allele: false,
                is_array: false,
                is_positional: true,
            }
        },
        "revel" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "revel".into(),
            name: "REVEL".into(),
            version: "latest".into(),
            description: format!("REVEL missense pathogenicity scores for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "spliceai" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "spliceAI".into(),
            name: "SpliceAI".into(),
            version: "latest".into(),
            description: format!("SpliceAI splice site predictions for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "primateai" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "primateAI".into(),
            name: "PrimateAI".into(),
            version: "latest".into(),
            description: format!("PrimateAI pathogenicity predictions for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "dbnsfp" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "dbnsfp".into(),
            name: "dbNSFP".into(),
            version: "latest".into(),
            description: format!("dbNSFP SIFT/PolyPhen predictions for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "cosmic" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "cosmic".into(),
            name: "COSMIC".into(),
            version: "latest".into(),
            description: format!("COSMIC somatic mutations for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "onekg" | "1000g" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "oneKg".into(),
            name: "1000 Genomes".into(),
            version: "latest".into(),
            description: format!("1000 Genomes population frequencies for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "topmed" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "topmed".into(),
            name: "TOPMed".into(),
            version: "latest".into(),
            description: format!("TOPMed population frequencies for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "mitomap" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "mitomap".into(),
            name: "MitoMap".into(),
            version: "latest".into(),
            description: format!("MitoMap mitochondrial variants for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        _ => anyhow::bail!(
            "Unknown source: {}. Supported: clinvar, gnomad, dbsnp, cosmic, onekg, topmed, mitomap, phylop, gerp, dann, revel, spliceai, primateai, dbnsfp",
            source
        ),
    };

    eprintln!("Building {} .osa from: {}", source, input);

    let file = File::open(input)
        .with_context(|| format!("Opening input file: {}", input))?;
    let reader: Box<dyn io::Read> = if input.ends_with(".gz") || input.ends_with(".bgz") {
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let buf_reader = io::BufReader::new(reader);

    let records = match source {
        "clinvar" => fastvep_sa::sources::clinvar::parse_clinvar_vcf(buf_reader, &chrom_map)?,
        "gnomad" => fastvep_sa::sources::gnomad::parse_gnomad_vcf(buf_reader, &chrom_map)?,
        "dbsnp" => fastvep_sa::sources::dbsnp::parse_dbsnp_vcf(buf_reader, &chrom_map)?,
        "cosmic" => fastvep_sa::sources::cosmic::parse_cosmic_vcf(buf_reader, &chrom_map)?,
        "onekg" | "1000g" => fastvep_sa::sources::onekg::parse_onekg_vcf(buf_reader, &chrom_map)?,
        "topmed" => fastvep_sa::sources::topmed::parse_topmed_vcf(buf_reader, &chrom_map)?,
        "mitomap" => fastvep_sa::sources::mitomap::parse_mitomap(buf_reader, &chrom_map)?,
        "phylop" => fastvep_sa::sources::scores::parse_wigfix(buf_reader, &chrom_map)?,
        "gerp" | "dann" => fastvep_sa::sources::scores::parse_score_tsv(buf_reader, &chrom_map, false)?,
        "revel" => fastvep_sa::sources::revel::parse_revel(buf_reader, &chrom_map, 2)?,
        "spliceai" => fastvep_sa::sources::spliceai::parse_spliceai_vcf(buf_reader, &chrom_map)?,
        "primateai" => fastvep_sa::sources::primateai::parse_primateai(buf_reader, &chrom_map)?,
        "dbnsfp" => fastvep_sa::sources::dbnsfp::parse_dbnsfp(buf_reader, &chrom_map)?,
        _ => unreachable!(),
    };

    eprintln!("Parsed {} records from {}", records.len(), source);

    let output_path = Path::new(output);
    let mut writer = SaWriter::new(header);
    writer.write_to_files(output_path, records.into_iter(), &chrom_list)?;

    eprintln!(
        "Wrote: {} and {}",
        output_path.with_extension("osa").display(),
        output_path.with_extension("osa.idx").display()
    );

    Ok(())
}

// SA provider loading is now in fastvep-annotate::load_sa_providers.
