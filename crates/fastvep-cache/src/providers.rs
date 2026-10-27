use anyhow::{Context, Result};
use fastvep_core::chrom_aliases;
use fastvep_genome::{is_mitochondrial, wrap_position_for, Transcript, MT_LENGTH};
use std::cell::RefCell;
use std::collections::{BTreeMap, HashMap};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

use crate::info::CacheInfo;
use crate::variation::{
    self, InputSite, Matched, VariationCacheRecord, VariationCursor, VariationTabixReader,
};

/// Trait for providing transcript annotations for a genomic region.
pub trait TranscriptProvider: Send + Sync {
    /// Return all transcripts that overlap the given genomic region.
    fn get_transcripts(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<&Transcript>>;

    /// Return all transcripts on a chromosome.
    fn get_transcripts_by_chrom(&self, chrom: &str) -> Result<Vec<&Transcript>>;
}

/// Trait for providing reference sequences.
pub trait SequenceProvider: Send + Sync {
    /// Fetch reference sequence for a region (1-based, inclusive).
    fn fetch_sequence(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>>;

    /// Fetch a reference sequence slice. Default delegates to fetch_sequence.
    fn fetch_sequence_slice(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        self.fetch_sequence(chrom, start, end)
    }

    /// The contigs this provider has sequence for, when it can say.
    fn sequence_names(&self) -> Option<Vec<String>> {
        None
    }
}

/// In-memory transcript provider backed by a Vec<Transcript>.
pub struct MemoryTranscriptProvider {
    transcripts: Vec<Transcript>,
}

impl MemoryTranscriptProvider {
    pub fn new(transcripts: Vec<Transcript>) -> Self {
        Self { transcripts }
    }

    pub fn transcript_count(&self) -> usize {
        self.transcripts.len()
    }
}

impl TranscriptProvider for MemoryTranscriptProvider {
    fn get_transcripts(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<&Transcript>> {
        Ok(self
            .transcripts
            .iter()
            .filter(|t| &*t.chromosome == chrom && t.start <= end && t.end >= start)
            .collect())
    }

    fn get_transcripts_by_chrom(&self, chrom: &str) -> Result<Vec<&Transcript>> {
        Ok(self
            .transcripts
            .iter()
            .filter(|t| &*t.chromosome == chrom)
            .collect())
    }
}

/// High-performance transcript provider using per-chromosome sorted arrays and
/// two binary searches to bound an overlap query.
///
/// Each chromosome's transcripts are sorted by start position. A parallel
/// `prefix_max_end` array stores `max(end)` over `transcripts[..=i]`. Because a
/// running maximum is monotonically non-decreasing, that array can itself be
/// binary-searched for the first transcript able to reach back to a query, so
/// only a bounded window is examined rather than every transcript that begins
/// at or before the query end.
///
/// A *prefix* maximum, not a suffix one. The question the scan has to answer is
/// "can anything at index `i` or earlier still reach this query?", and a
/// transcript nested inside a longer one that began earlier is the normal case
/// in a real gene model. A suffix maximum answers the opposite question, and
/// using it drops every transcript still open past the last transcript to
/// *begin* on the chromosome — 2,186 bp of Ensembl 116 chr22 covered by 29
/// transcripts were reported as intergenic (issue #126).
pub struct IndexedTranscriptProvider {
    /// Transcripts grouped by chromosome, sorted by start position within each group.
    by_chrom: HashMap<Arc<str>, Vec<Transcript>>,
    /// Prefix-max-end arrays: `prefix_max_end[chrom][i]` == max(end) over
    /// `by_chrom[chrom][..=i]`. Built once in [`new`](Self::new) alongside
    /// `by_chrom` and never mutated after, so it is always the same length as
    /// the list it indexes.
    prefix_max_end: HashMap<Arc<str>, Vec<u64>>,
}

/// The order [`IndexedTranscriptProvider`] keeps transcripts in: by contig, then by start, and
/// transcripts that tie as the file had them (the sort is stable, which is what makes the transcript
/// the one that comes first in the file the one a tie-break prefers). A list in this order is one
/// the provider has nothing to sort, and it is the order the sidecar transcript cache is written in.
pub fn sort_for_provider(transcripts: &mut [Transcript]) {
    transcripts.sort_by(|a, b| a.chromosome.cmp(&b.chromosome).then(a.start.cmp(&b.start)));
}

impl IndexedTranscriptProvider {
    pub fn new(mut transcripts: Vec<Transcript>) -> Self {
        let mut by_chrom: HashMap<Arc<str>, Vec<Transcript>> = HashMap::new();
        for tr in transcripts.drain(..) {
            by_chrom
                .entry(Arc::clone(&tr.chromosome))
                .or_default()
                .push(tr);
        }
        // Sort each chromosome's transcripts by start position, unless they are in that order, which
        // a sidecar cache is written in ([`sort_for_provider`]): a pass over them, and not a sort of
        // structs of some 600 bytes each, which was a tenth of the start of a run on the human GFF3
        for trs in by_chrom.values_mut() {
            if !trs.is_sorted_by_key(|t| t.start) {
                trs.sort_by_key(|t| t.start);
            }
        }
        // Build prefix-max-end arrays to bound the overlap window
        let mut prefix_max_end = HashMap::new();
        for (chrom, trs) in &by_chrom {
            let mut running = 0u64;
            let maxima = trs
                .iter()
                .map(|t| {
                    running = running.max(t.end);
                    running
                })
                .collect();
            prefix_max_end.insert(Arc::clone(chrom), maxima);
        }
        Self {
            by_chrom,
            prefix_max_end,
        }
    }

    pub fn transcript_count(&self) -> usize {
        self.by_chrom.values().map(|v| v.len()).sum()
    }

    /// Resolve a query chromosome to the stored key, falling back to chr↔bare
    /// and mitochondrial aliases. Lets a `chr17` VCF match a cache built with
    /// `17` (and vice versa) instead of silently returning no transcripts.
    fn resolve_key(&self, chrom: &str) -> Option<&str> {
        if let Some((k, _)) = self.by_chrom.get_key_value(chrom) {
            return Some(k.as_ref());
        }
        chrom_aliases(chrom).into_iter().find_map(|alias| {
            self.by_chrom
                .get_key_value(alias.as_str())
                .map(|(k, _)| k.as_ref())
        })
    }

    /// The half-open `[lo, hi)` slice of `by_chrom[key]` that can contain a
    /// transcript overlapping `[start, end]`. `key` must already be resolved.
    ///
    /// * `hi` — one past the last transcript that begins at or before `end`,
    ///   by binary search over the start-sorted list. Anything later begins
    ///   after the query and cannot overlap it.
    /// * `lo` — the first transcript that can still reach back to `start`, by
    ///   binary search over the monotone `prefix_max_end`. Every transcript
    ///   before `lo` ends strictly before `start`, so none can overlap.
    ///
    /// Split out from [`get_transcripts`](TranscriptProvider::get_transcripts)
    /// so the window itself is testable and not just the rows it yields: a
    /// pruning bug that only *widens* the window returns correct results while
    /// costing a full scan per query, which is how the suffix-max-end version
    /// of this index came to scan an average of 5,257 entries per chr22 query
    /// where 17 suffice.
    fn candidate_window(&self, key: &str, start: u64, end: u64) -> (usize, usize) {
        let trs = &self.by_chrom[key];
        let maxima = &self.prefix_max_end[key];
        debug_assert_eq!(maxima.len(), trs.len(), "prefix_max_end out of sync");

        let hi = trs.partition_point(|t| t.start <= end);
        let lo = maxima.partition_point(|&max_end| max_end < start);
        (lo.min(hi), hi)
    }
}

impl TranscriptProvider for IndexedTranscriptProvider {
    fn get_transcripts(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<&Transcript>> {
        let Some(key) = self.resolve_key(chrom) else {
            return Ok(Vec::new());
        };
        let (lo, hi) = self.candidate_window(key, start, end);
        Ok(self.by_chrom[key][lo..hi]
            .iter()
            .filter(|t| t.end >= start)
            .collect())
    }

    fn get_transcripts_by_chrom(&self, chrom: &str) -> Result<Vec<&Transcript>> {
        match self.resolve_key(chrom) {
            Some(key) => Ok(self.by_chrom[key].iter().collect()),
            None => Ok(Vec::new()),
        }
    }
}

/// Fetch a genomic region, wrapping around the origin for circular
/// (mitochondrial) contigs when `end` extends past the contig's length.
///
/// A linear FASTA record for `MT`/`chrM` only has bases `1..=length`, so a
/// variant whose reference allele runs past the last base (e.g. a 2bp allele
/// starting at the human rCRS's last position, 16569) logically continues
/// from position 1. Without this, `fetch` on the underlying reader either
/// errors (mmap/.fai readers clamp `end` to `length` and would silently
/// return a truncated, wrong-length slice) or simply omits the wrapped
/// bases — see fastVEP issue #68.
///
/// Only mitochondrial contigs get this treatment; all other chromosomes (and
/// any MT `end` that doesn't actually exceed the contig length) take the
/// untouched fast path straight through to `fetch`.
fn fetch_circular<F>(
    chrom: &str,
    start: u64,
    end: u64,
    contig_len: Option<u64>,
    fetch: F,
) -> Result<Vec<u8>>
where
    F: Fn(&str, u64, u64) -> Result<Vec<u8>>,
{
    if !is_mitochondrial(chrom) {
        return fetch(chrom, start, end);
    }
    let length = contig_len.unwrap_or(MT_LENGTH);
    if length == 0 || end <= length {
        return fetch(chrom, start, end);
    }

    // The region wraps past the origin: take [start, length] then [1, wrapped_end].
    let mut result = fetch(chrom, start, length)?;
    let wrapped_end = wrap_position_for(end, length);
    if wrapped_end > 0 {
        let rest = fetch(chrom, 1, wrapped_end)?;
        result.extend_from_slice(&rest);
    }
    Ok(result)
}

/// Sequence provider backed by a FASTA reader.
pub struct FastaSequenceProvider {
    reader: crate::fasta::FastaReader,
}

impl FastaSequenceProvider {
    pub fn new(reader: crate::fasta::FastaReader) -> Self {
        Self { reader }
    }
}

impl SequenceProvider for FastaSequenceProvider {
    fn fetch_sequence(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        fetch_circular(
            chrom,
            start,
            end,
            self.reader.sequence_length(chrom),
            |c, s, e| self.reader.fetch(c, s, e),
        )
    }

    fn sequence_names(&self) -> Option<Vec<String>> {
        Some(
            self.reader
                .sequence_names()
                .into_iter()
                .map(String::from)
                .collect(),
        )
    }
}

/// Sequence provider backed by a memory-mapped FASTA reader.
/// Uses .fai index for random access without loading the full file into RAM.
pub struct MmapFastaSequenceProvider {
    reader: crate::fasta::MmapFastaReader,
}

impl MmapFastaSequenceProvider {
    pub fn new(reader: crate::fasta::MmapFastaReader) -> Self {
        Self { reader }
    }
}

impl SequenceProvider for MmapFastaSequenceProvider {
    fn fetch_sequence(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        fetch_circular(
            chrom,
            start,
            end,
            self.reader.sequence_length(chrom),
            |c, s, e| self.reader.fetch(c, s, e),
        )
    }

    fn sequence_names(&self) -> Option<Vec<String>> {
        Some(
            self.reader
                .sequence_names()
                .into_iter()
                .map(String::from)
                .collect(),
        )
    }
}

/// A matched known variant with its allele-specific frequency data.
///
/// This is the per-allele view `Existing_variation` and fastVEP's own CSQ
/// summaries read. VEP's frequency columns are not computed from it: they need
/// the whole record and what it matched, which is [`KnownSite`].
#[derive(Debug, Clone)]
pub struct MatchedVariant {
    pub name: String,
    /// The allele as the cache record writes it, the one its per-population
    /// frequencies are keyed by. Empty for a record that carries no alleles.
    pub matched_allele: String,
    pub minor_allele: Option<String>,
    pub minor_allele_freq: Option<f64>,
    pub clin_sig: Option<String>,
    pub somatic: bool,
    pub phenotype_or_disease: bool,
    pub pubmed: Vec<String>,
    /// Population → allele-specific frequency for the matched allele.
    pub frequencies: HashMap<String, f64>,
}

impl MatchedVariant {
    fn from_record(record: &VariationCacheRecord, record_allele: &str) -> Self {
        // Per-allele frequencies for the matched allele
        let mut frequencies = HashMap::new();
        for (pop, freq_str) in &record.frequencies {
            if let Some(f) = variation::get_allele_freq(freq_str, record_allele) {
                frequencies.insert(pop.clone(), f);
            }
        }

        // Also include MAF if minor_allele matches
        if let (Some(ma), Some(maf)) = (&record.minor_allele, record.minor_allele_freq) {
            if ma.eq_ignore_ascii_case(record_allele) {
                frequencies.entry("minor_allele_freq".into()).or_insert(maf);
            }
        }

        MatchedVariant {
            name: record.variation_name.clone(),
            matched_allele: record_allele.to_string(),
            minor_allele: record.minor_allele.clone(),
            minor_allele_freq: record.minor_allele_freq,
            clin_sig: record.clin_sig.clone(),
            somatic: record.somatic,
            phenotype_or_disease: record.phenotype_or_disease,
            pubmed: record.pubmed.clone(),
            frequencies,
        }
    }
}

/// A known variant VEP attached to an input site, and which input alleles it
/// matched.
#[derive(Debug, Clone)]
pub struct KnownRecord {
    /// Shared, so that the alleles of a site and the output rows read one copy.
    pub record: Arc<VariationCacheRecord>,
    pub matched: Matched,
}

/// Every known variant attached to one input site, in the order the cache
/// holds them: `$vf->{existing}` of VEP.
///
/// The per-allele lists of [`MatchedVariant`] are not enough for VEP's frequency
/// columns. `add_colocated_frequency_data` runs over *all* the site's known
/// variants for every allele, including those that matched only another allele
/// of a multi-allelic record, and it reads them in this order, not in the order
/// `Existing_variation` is sorted into.
#[derive(Debug, Clone, Default)]
pub struct KnownSite {
    pub records: Vec<KnownRecord>,
}

impl KnownSite {
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// The known variants of each of the site's `alt_count` alleles, each list in
    /// the order VEP prints `Existing_variation`: germline before somatic, then by
    /// the prefix of the name. The sort is stable, so ties keep cache order.
    pub fn per_alt(&self, alt_count: usize) -> Vec<Vec<MatchedVariant>> {
        let mut per_alt: Vec<Vec<MatchedVariant>> = vec![Vec::new(); alt_count];
        for known in &self.records {
            for (alt_index, list) in per_alt.iter_mut().enumerate() {
                match &known.matched {
                    Matched::AllAlleles => {
                        list.push(MatchedVariant::from_record(&known.record, ""))
                    }
                    Matched::Alleles(matches) => {
                        // The first matching record allele speaks for the input
                        // allele, as the first of `matched_alleles` does.
                        if let Some(m) = matches.iter().find(|m| m.alt_index == alt_index) {
                            list.push(MatchedVariant::from_record(&known.record, &m.record_allele));
                        }
                    }
                }
            }
        }
        for list in &mut per_alt {
            list.sort_by_key(|m| variation::existing_variation_rank(&m.name, m.somatic));
        }
        per_alt
    }
}

/// Trait for providing co-located known variant annotations.
pub trait VariationProvider {
    /// The known variants VEP's `--check_existing` attaches to one input site.
    ///
    /// `start` and `end` are Ensembl coordinates (an insertion has
    /// `end == start - 1`); alleles are written as in an Ensembl allele string,
    /// `-` for an empty one. `AlleleMatch::alt_index` of the result counts into
    /// `alt_alleles`.
    ///
    /// An error means the lookup could not be done, and the caller must not
    /// read it as "no known variants".
    fn known_site(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        ref_allele: &str,
        alt_alleles: &[String],
    ) -> Result<KnownSite>;

    /// [`known_site`](Self::known_site) as one list per alt allele:
    /// `result[i]` belongs to `alt_alleles[i]`, in the order VEP prints
    /// `Existing_variation`.
    fn get_matched_variants(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        ref_allele: &str,
        alt_alleles: &[String],
    ) -> Result<Vec<Vec<MatchedVariant>>> {
        Ok(self
            .known_site(chrom, start, end, ref_allele, alt_alleles)?
            .per_alt(alt_alleles.len()))
    }

    /// What to tell the user once the run is over, when lookups came back
    /// empty for a reason that is theirs to fix.
    fn finish(&self) -> Option<String> {
        None
    }
}

/// Variation provider backed by VEP's tabix-indexed cache files.
pub struct TabixVariationProvider {
    reader: VariationTabixReader,
    cache_dir: std::path::PathBuf,
    /// Sites looked up.
    sites: AtomicU64,
    /// Contigs of the input that the cache does not know, with the number of
    /// sites on each.
    unresolved: Mutex<BTreeMap<String, u64>>,
}

impl TabixVariationProvider {
    /// Create a provider from a VEP cache directory.
    ///
    /// Reads `info.txt` for column definitions and valid chromosomes. Fails
    /// when the directory is not a tabix variation cache, or holds a variation
    /// file that cannot be searched; see [`VariationTabixReader::new`].
    pub fn new(cache_dir: &std::path::Path, cache_info: &CacheInfo) -> Result<Self> {
        let reader = VariationTabixReader::new(cache_dir, &cache_info.variation_cols)
            .with_context(|| {
                format!(
                    "--cache-dir {} (info.txt: var_type {})",
                    cache_dir.display(),
                    cache_info.var_type.as_deref().unwrap_or("absent")
                )
            })?;
        Ok(Self {
            reader,
            cache_dir: cache_dir.to_path_buf(),
            sites: AtomicU64::new(0),
            unresolved: Mutex::new(BTreeMap::new()),
        })
    }

    fn note_unresolved(&self, chrom: &str) {
        // A poisoned lock only loses a count that feeds a warning.
        if let Ok(mut unresolved) = self.unresolved.lock() {
            *unresolved.entry(chrom.to_string()).or_insert(0) += 1;
        }
    }
}

impl TabixVariationProvider {
    /// A place to look the sites of one thread up from, which carries on from its own last site
    /// and does not share a stream of windows with the provider's other callers.
    ///
    /// Run `known_site` on it for sites in input order and a sorted input is read through once,
    /// by each cursor over the stretch it was given.
    pub fn cursor(&self) -> TabixCursor<'_> {
        TabixCursor {
            provider: self,
            cursor: RefCell::new(self.reader.cursor()),
        }
    }

    /// The known variants of one input site, in cache order.
    fn lookup(
        &self,
        query: impl FnOnce(&str, u64, u64) -> Result<Vec<VariationCacheRecord>>,
        chrom: &str,
        site: &InputSite,
    ) -> Result<KnownSite> {
        // VEP asks the index for `start - 1 ..= end + 1` of the input
        // (VariationTabix::_annotate_pm), and the index holds each record at its
        // start.
        let records = query(
            chrom,
            site.start.saturating_sub(1),
            site.end.saturating_add(1),
        )
        .with_context(|| {
            format!(
                "known-variant lookup at {chrom}:{}-{}",
                site.start, site.end
            )
        })?;
        if records.is_empty() {
            return Ok(KnownSite::default());
        }

        let alts: Vec<&str> = site.alts.iter().map(String::as_str).collect();
        let forms = variation::SiteForms::new(&site.reference, &alts, site.start, site.end);
        let mut known: Vec<KnownRecord> = Vec::new();

        for record in records {
            // Failed variants are out unless asked for (`filter_variation`), and
            // VEP skips a record whose name it has already matched for this
            // site, whichever allele that match was for.
            if record.failed
                || known
                    .iter()
                    .any(|k| k.record.variation_name == record.variation_name)
            {
                continue;
            }
            let Some(matched) = forms.compare(&record) else {
                continue;
            };
            known.push(KnownRecord {
                record: Arc::new(record),
                matched,
            });
        }
        Ok(KnownSite { records: known })
    }
}

impl TabixVariationProvider {
    /// `known_site`, reading the windows by `query`.
    fn known_site_by(
        &self,
        query: impl FnOnce(&str, u64, u64) -> Result<Vec<VariationCacheRecord>>,
        chrom: &str,
        start: u64,
        end: u64,
        ref_allele: &str,
        alt_alleles: &[String],
    ) -> Result<KnownSite> {
        self.sites.fetch_add(1, Ordering::Relaxed);

        if self.reader.normalize_chrom(chrom).is_none() {
            self.note_unresolved(chrom);
            return Ok(KnownSite::default());
        }

        // A lone indel alt is looked up minimised, as VEP does it; see
        // `vep_input_site`.
        let alts: Vec<&str> = alt_alleles.iter().map(String::as_str).collect();
        self.lookup(
            query,
            chrom,
            &variation::vep_input_site(ref_allele, &alts, start, end),
        )
    }
}

/// A [`TabixVariationProvider`] read from one thread, with a stream of windows of its own; see
/// [`TabixVariationProvider::cursor`]. It answers as the provider does, and counts the sites it
/// is asked about and the contigs it cannot place in the provider's totals.
pub struct TabixCursor<'a> {
    provider: &'a TabixVariationProvider,
    cursor: RefCell<VariationCursor<'a>>,
}

impl VariationProvider for TabixCursor<'_> {
    fn known_site(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        ref_allele: &str,
        alt_alleles: &[String],
    ) -> Result<KnownSite> {
        self.provider.known_site_by(
            |chrom, lo, hi| self.cursor.borrow_mut().query(chrom, lo, hi),
            chrom,
            start,
            end,
            ref_allele,
            alt_alleles,
        )
    }
}

impl VariationProvider for TabixVariationProvider {
    fn known_site(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        ref_allele: &str,
        alt_alleles: &[String],
    ) -> Result<KnownSite> {
        self.known_site_by(
            |chrom, lo, hi| self.reader.query(chrom, lo, hi),
            chrom,
            start,
            end,
            ref_allele,
            alt_alleles,
        )
    }

    fn finish(&self) -> Option<String> {
        let unresolved = self.unresolved.lock().ok()?;
        if unresolved.is_empty() {
            return None;
        }
        let sites = self.sites.load(Ordering::Relaxed);
        let missed: u64 = unresolved.values().sum();
        let mut examples: Vec<String> = unresolved
            .iter()
            .take(5)
            .map(|(chrom, n)| format!("{chrom} ({n})"))
            .collect();
        if unresolved.len() > 5 {
            examples.push(format!("and {} more", unresolved.len() - 5));
        }
        let which = if missed == sites {
            "none of the input's variants is on a contig".to_string()
        } else {
            format!("{missed} of {sites} variants are on contigs")
        };
        Some(format!(
            "warning: {which} that --cache-dir {} knows, so Existing_variation is empty for \
             {}: {}. The cache has {}. If the names differ, put a chr_synonyms.txt in the cache \
             directory (the file VEP reads, one line of equivalent names per contig).",
            self.cache_dir.display(),
            if missed == sites {
                "all of them"
            } else {
                "those"
            },
            examples.join(", "),
            self.reader.chromosomes().join(", ")
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fastvep_core::Strand;
    use fastvep_genome::{Exon, Gene};

    fn make_transcript(chrom: &str, start: u64, end: u64) -> Transcript {
        Transcript {
            stable_id: Arc::from(format!("ENST_{}", start).as_str()),
            version: None,
            gene: Gene {
                stable_id: "ENSG_1".into(),
                symbol: None,
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: chrom.into(),
                start,
                end,
                strand: Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: chrom.into(),
            start,
            end,
            strand: Strand::Forward,
            exons: vec![Exon {
                stable_id: "ENSE_1".into(),
                start,
                end,
                strand: Strand::Forward,
                phase: 0,
                end_phase: 0,
                rank: 1,
            }],
            translation: None,
            cdna_coding_start: None,
            cdna_coding_end: None,
            coding_region_start: None,
            coding_region_end: None,
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: false,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: None,
            protein_version: None,
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: None,
            gencode_primary: false,
            flags: vec![],
            codon_table_start_phase: 0,
            initiator_met: false,
            cds_flags_from_sequence: false,
            readthrough: false,
        }
    }

    /// A list sorted for the provider answers every query as the list in the order the file had
    /// it does, for transcripts that tie on their start too, on several contigs, and keeps the
    /// order of the ones that tie.
    #[test]
    fn a_list_sorted_for_the_provider_answers_as_the_list_in_file_order_does() {
        let mut seed = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move |bound: u64| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed % bound
        };
        let mut in_file_order: Vec<Transcript> = (0..600)
            .map(|n| {
                let chrom = ["chr1", "chr2", "chrX", "chrM"][next(4) as usize];
                // starts from a small range, so that many tie, and some that reach far
                let start = 1 + next(40) * 100;
                let mut t = make_transcript(chrom, start, start + next(600));
                t.stable_id = Arc::from(format!("T{n:03}").as_str());
                t
            })
            .collect();
        in_file_order.sort_by_key(|t| (t.start % 7, t.stable_id.len())); // not in any useful order
        let mut sorted = in_file_order.clone();
        sort_for_provider(&mut sorted);
        // by contig, then start, and the ties in the order the file had them
        assert!(sorted
            .windows(2)
            .all(|w| (&w[0].chromosome, w[0].start) <= (&w[1].chromosome, w[1].start)));
        for pair in sorted.windows(2) {
            if pair[0].chromosome == pair[1].chromosome && pair[0].start == pair[1].start {
                let at = |t: &Transcript| {
                    in_file_order
                        .iter()
                        .position(|o| o.stable_id == t.stable_id)
                        .unwrap()
                };
                assert!(at(&pair[0]) < at(&pair[1]), "a tie changed places");
            }
        }

        let from_file = IndexedTranscriptProvider::new(in_file_order);
        let from_sorted = IndexedTranscriptProvider::new(sorted);
        assert_eq!(from_file.transcript_count(), from_sorted.transcript_count());
        for chrom in ["chr1", "chr2", "chrX", "chrM", "chr9"] {
            for _ in 0..200 {
                let start = 1 + next(4200);
                let end = start + next(300);
                let ids = |p: &IndexedTranscriptProvider| -> Vec<String> {
                    p.get_transcripts(chrom, start, end)
                        .unwrap()
                        .iter()
                        .map(|t| t.stable_id.to_string())
                        .collect()
                };
                assert_eq!(ids(&from_file), ids(&from_sorted), "{chrom}:{start}-{end}");
            }
            let all = |p: &IndexedTranscriptProvider| -> Vec<String> {
                p.get_transcripts_by_chrom(chrom)
                    .unwrap()
                    .iter()
                    .map(|t| t.stable_id.to_string())
                    .collect()
            };
            assert_eq!(all(&from_file), all(&from_sorted), "{chrom}");
        }
    }

    #[test]
    fn test_memory_transcript_provider() {
        let provider = MemoryTranscriptProvider::new(vec![
            make_transcript("chr1", 1000, 2000),
            make_transcript("chr1", 3000, 4000),
            make_transcript("chr2", 1000, 2000),
        ]);

        // Overlapping query
        let results = provider.get_transcripts("chr1", 1500, 1600).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].start, 1000);

        // Non-overlapping
        let results = provider.get_transcripts("chr1", 2500, 2600).unwrap();
        assert_eq!(results.len(), 0);

        // Different chromosome
        let results = provider.get_transcripts("chr2", 1500, 1600).unwrap();
        assert_eq!(results.len(), 1);

        // By chromosome
        let results = provider.get_transcripts_by_chrom("chr1").unwrap();
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn test_indexed_transcript_provider() {
        let provider = IndexedTranscriptProvider::new(vec![
            make_transcript("chr1", 1000, 2000),
            make_transcript("chr1", 3000, 4000),
            make_transcript("chr2", 1000, 2000),
        ]);

        // Overlapping query
        let results = provider.get_transcripts("chr1", 1500, 1600).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].start, 1000);

        // Non-overlapping
        let results = provider.get_transcripts("chr1", 2500, 2600).unwrap();
        assert_eq!(results.len(), 0);

        // Different chromosome
        let results = provider.get_transcripts("chr2", 1500, 1600).unwrap();
        assert_eq!(results.len(), 1);

        // By chromosome
        let results = provider.get_transcripts_by_chrom("chr1").unwrap();
        assert_eq!(results.len(), 2);

        // Missing chromosome
        let results = provider.get_transcripts("chr99", 1, 100).unwrap();
        assert_eq!(results.len(), 0);

        // Query spanning two transcripts
        let results = provider.get_transcripts("chr1", 1500, 3500).unwrap();
        assert_eq!(results.len(), 2);

        // Query at exact boundaries
        let results = provider.get_transcripts("chr1", 2000, 2000).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].start, 1000);

        // Transcript count
        assert_eq!(provider.transcript_count(), 3);
    }

    /// Intervals as `(start, end)`, in the order the provider reports them.
    fn spans(rows: Vec<&Transcript>) -> Vec<(u64, u64)> {
        rows.into_iter().map(|t| (t.start, t.end)).collect()
    }

    /// A short transcript that both begins and ends after a longer one, so it
    /// sorts last by start while the longer one is still open past its end.
    /// The suffix-max-end version of the scan stopped at the short one and
    /// reported the region past it as having no transcripts at all (#126).
    #[test]
    fn indexed_provider_retains_nested_long_transcript() {
        let provider = IndexedTranscriptProvider::new(vec![
            make_transcript("chr1", 1000, 5000),
            make_transcript("chr1", 3000, 4000),
        ]);
        let results = provider.get_transcripts("chr1", 4500, 4501).unwrap();
        assert_eq!(spans(results), vec![(1000, 5000)]);
    }

    /// Every query against a nested set must agree with a plain filter. The
    /// fixture deliberately mixes containment, shared starts and a transcript
    /// that outlives every later start.
    #[test]
    fn indexed_provider_matches_linear_nested_queries() {
        let transcripts = vec![
            make_transcript("chr1", 3000, 4000),
            make_transcript("chr1", 1000, 6000),
            make_transcript("chr1", 1100, 1600),
            make_transcript("chr1", 1000, 1500),
            make_transcript("chr1", 2000, 5000),
        ];
        let indexed = IndexedTranscriptProvider::new(transcripts.clone());
        let linear = MemoryTranscriptProvider::new(transcripts);

        for start in (900..=6200).step_by(37) {
            for width in [0, 1, 100, 5000] {
                // `MemoryTranscriptProvider` reports in input order and the
                // index in start order, so compare as sorted multisets: the
                // index's own ordering is pinned by the assertion below.
                let mut expected = spans(
                    linear
                        .get_transcripts("chr1", start, start + width)
                        .unwrap(),
                );
                let mut actual = spans(
                    indexed
                        .get_transcripts("chr1", start, start + width)
                        .unwrap(),
                );
                let ordered = actual.clone();
                expected.sort_unstable();
                actual.sort_unstable();
                assert_eq!(actual, expected, "query {start} + {width}");

                let mut by_start = ordered.clone();
                by_start.sort_by_key(|&(s, _)| s);
                assert_eq!(
                    ordered, by_start,
                    "results not start-sorted at {start} + {width}"
                );
            }
        }
    }

    /// Correct results alone cannot tell a working prune from a broken one: a
    /// window that is too wide still yields the right rows. Pin the window.
    #[test]
    fn candidate_window_excludes_transcripts_that_end_before_the_query() {
        let provider = IndexedTranscriptProvider::new(
            (0..100)
                .map(|i| make_transcript("chr1", i * 100 + 1, i * 100 + 50))
                .collect(),
        );

        // Non-overlapping transcripts: a point query needs a window of one.
        let (lo, hi) = provider.candidate_window("chr1", 5001, 5001);
        assert_eq!((lo, hi), (50, 51));

        // A transcript spanning everything forces the window open from 0,
        // because the running maximum never drops below the query.
        let with_spanner = IndexedTranscriptProvider::new(
            std::iter::once(make_transcript("chr1", 1, 10_000))
                .chain((0..100).map(|i| make_transcript("chr1", i * 100 + 1, i * 100 + 50)))
                .collect(),
        );
        let (lo, _) = with_spanner.candidate_window("chr1", 5001, 5001);
        assert_eq!(
            lo, 0,
            "a chromosome-spanning transcript must keep the window open"
        );
    }

    #[test]
    fn test_indexed_provider_overlapping_transcripts() {
        // Test with overlapping transcripts (common in real genomes)
        let provider = IndexedTranscriptProvider::new(vec![
            make_transcript("chr1", 1000, 5000),
            make_transcript("chr1", 2000, 3000),
            make_transcript("chr1", 4000, 6000),
        ]);

        // Query in the overlap region
        let results = provider.get_transcripts("chr1", 2500, 2600).unwrap();
        assert_eq!(results.len(), 2); // Both 1000-5000 and 2000-3000

        // Query spanning all three
        let results = provider.get_transcripts("chr1", 1000, 6000).unwrap();
        assert_eq!(results.len(), 3);
    }

    #[test]
    fn test_fasta_sequence_provider() {
        let fasta = ">chr1\nACGTACGTAAAACCCC\n";
        let reader = crate::fasta::FastaReader::from_reader(fasta.as_bytes()).unwrap();
        let provider = FastaSequenceProvider::new(reader);
        let seq = provider.fetch_sequence("chr1", 1, 4).unwrap();
        assert_eq!(seq, b"ACGT");
    }

    #[test]
    fn test_fasta_fetch_slice_matches_fetch() {
        let fasta = ">chr1\nacgtACGTaaaa\n>chr2\nTTTTgggg\n";
        let reader = crate::fasta::FastaReader::from_reader(fasta.as_bytes()).unwrap();
        let provider = FastaSequenceProvider::new(reader);

        // fetch_slice returns same data as fetch (both uppercase)
        let slice = provider.fetch_sequence_slice("chr1", 1, 4).unwrap();
        let vec = provider.fetch_sequence("chr1", 1, 4).unwrap();
        assert_eq!(slice, vec.as_slice());
        assert_eq!(slice, b"ACGT");

        // Lowercase input is uppercased at load time
        let slice = provider.fetch_sequence_slice("chr2", 5, 8).unwrap();
        assert_eq!(slice, b"GGGG");
    }

    #[test]
    fn test_mt_fetch_wraps_around_origin() {
        // A short synthetic "MT" contig: 10 bases, so a query for [8, 12]
        // should wrap and read positions 8,9,10 then 1,2 -> "HIJAB".
        let fasta = ">MT\nABCDEFGHIJ\n";
        let reader = crate::fasta::FastaReader::from_reader(fasta.as_bytes()).unwrap();
        let provider = FastaSequenceProvider::new(reader);

        let seq = provider.fetch_sequence("MT", 8, 12).unwrap();
        assert_eq!(seq, b"HIJAB");

        // A non-wrapping MT query behaves exactly like a normal fetch.
        let seq = provider.fetch_sequence("MT", 1, 4).unwrap();
        assert_eq!(seq, b"ABCD");
    }

    #[test]
    fn test_mt_fetch_wraps_recognizes_chrm_alias() {
        // Same as above but using the UCSC `chrM` alias for both the FASTA
        // record and the query, to confirm is_mitochondrial + alias
        // resolution compose correctly.
        let fasta = ">chrM\nABCDEFGHIJ\n";
        let reader = crate::fasta::FastaReader::from_reader(fasta.as_bytes()).unwrap();
        let provider = FastaSequenceProvider::new(reader);

        let seq = provider.fetch_sequence("MT", 9, 11).unwrap();
        assert_eq!(seq, b"IJA");
    }

    #[test]
    fn test_non_mt_chrom_is_never_wrapped() {
        // A non-mitochondrial chromosome with `end` past its length must
        // fail like a normal out-of-range fetch, not silently wrap.
        let fasta = ">chr1\nACGT\n";
        let reader = crate::fasta::FastaReader::from_reader(fasta.as_bytes()).unwrap();
        let provider = FastaSequenceProvider::new(reader);
        // fetch_slice-backed fetch clamps end to length rather than erroring;
        // assert it returns the clamped (not wrapped) sequence.
        let seq = provider.fetch_sequence("chr1", 1, 10).unwrap();
        assert_eq!(seq, b"ACGT");
    }

    #[test]
    fn test_mmap_mt_fetch_wraps_around_origin() {
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!(
            "fastvep-cache-test-mt-wrap-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let fasta_path = dir.join("mt.fa");
        std::fs::write(&fasta_path, ">MT\nABCDEFGHIJ\n").unwrap();
        let fai_path = dir.join("mt.fa.fai");
        let mut fai = std::fs::File::create(&fai_path).unwrap();
        // name, length, offset (after ">MT\n" = 4 bytes), line_bases, line_bytes
        writeln!(fai, "MT\t10\t4\t10\t11").unwrap();

        let reader = crate::fasta::MmapFastaReader::open(&fasta_path).unwrap();
        let provider = MmapFastaSequenceProvider::new(reader);

        let seq = provider.fetch_sequence("MT", 8, 12).unwrap();
        assert_eq!(seq, b"HIJAB");

        let _ = std::fs::remove_dir_all(&dir);
    }
}
