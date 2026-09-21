use anyhow::Result;
use fastvep_core::chrom_aliases;
use fastvep_genome::{is_mitochondrial, wrap_position_for, Transcript, MT_LENGTH};
use std::collections::HashMap;
use std::sync::Arc;

use crate::info::CacheInfo;
use crate::variation::{self, VariationTabixReader};

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

impl IndexedTranscriptProvider {
    pub fn new(mut transcripts: Vec<Transcript>) -> Self {
        let mut by_chrom: HashMap<Arc<str>, Vec<Transcript>> = HashMap::new();
        for tr in transcripts.drain(..) {
            by_chrom
                .entry(Arc::clone(&tr.chromosome))
                .or_default()
                .push(tr);
        }
        // Sort each chromosome's transcripts by start position
        for trs in by_chrom.values_mut() {
            trs.sort_by_key(|t| t.start);
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
}

/// A matched known variant with its allele-specific frequency data.
#[derive(Debug, Clone)]
pub struct MatchedVariant {
    pub name: String,
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

/// Trait for providing co-located known variant annotations.
pub trait VariationProvider {
    /// Look up known variants overlapping a position that match the given alleles.
    fn get_matched_variants(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Vec<MatchedVariant>>;
}

/// Variation provider backed by VEP's tabix-indexed cache files.
pub struct TabixVariationProvider {
    reader: VariationTabixReader,
}

impl TabixVariationProvider {
    /// Create a provider from a VEP cache directory.
    ///
    /// Reads `info.txt` for column definitions and valid chromosomes.
    pub fn new(cache_dir: &std::path::Path, cache_info: &CacheInfo) -> Result<Self> {
        let reader = VariationTabixReader::new(
            cache_dir,
            &cache_info.variation_cols,
            &cache_info.valid_chromosomes,
        )?;
        Ok(Self { reader })
    }
}

impl VariationProvider for TabixVariationProvider {
    fn get_matched_variants(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Vec<MatchedVariant>> {
        let records = self.reader.query(chrom, start, end)?;
        let mut matched = Vec::new();

        for record in &records {
            // Skip failed variants
            if record.failed {
                continue;
            }

            // Check allele match
            if let Some(matched_alt) =
                variation::match_alleles(ref_allele, alt_allele, start, end, record)
            {
                // Extract per-allele frequencies for the matched allele
                let mut freqs = HashMap::new();
                for (pop, freq_str) in &record.frequencies {
                    if let Some(f) = variation::get_allele_freq(freq_str, &matched_alt) {
                        freqs.insert(pop.clone(), f);
                    }
                }

                // Also include MAF if minor_allele matches
                if let (Some(ref ma), Some(maf)) = (&record.minor_allele, record.minor_allele_freq)
                {
                    if ma.eq_ignore_ascii_case(&matched_alt) {
                        freqs.entry("minor_allele_freq".into()).or_insert(maf);
                    }
                }

                matched.push(MatchedVariant {
                    name: record.variation_name.clone(),
                    matched_allele: matched_alt,
                    minor_allele: record.minor_allele.clone(),
                    minor_allele_freq: record.minor_allele_freq,
                    clin_sig: record.clin_sig.clone(),
                    somatic: record.somatic,
                    phenotype_or_disease: record.phenotype_or_disease,
                    pubmed: record.pubmed.clone(),
                    frequencies: freqs,
                });
            }
        }

        Ok(matched)
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
