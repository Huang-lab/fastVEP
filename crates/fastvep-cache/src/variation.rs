//! VEP's variation cache, as `convert_cache.pl` leaves it: one bgzipped,
//! tabix-indexed `all_vars.gz` per chromosome directory, holding the known
//! variants (rsIDs and the like) that `--check_existing` reports as
//! `Existing_variation`.
//!
//! This module reads those files and decides which records are the *same
//! variant* as an input site. Both halves reproduce VEP (ensembl-vep
//! release/115) rather than reinterpret it, because a wrong answer here is a
//! plausible-looking one: a known variant reported as novel, or the reverse.
//! The Perl that defines each rule is named next to it.

use anyhow::{anyhow, bail, Context, Result};
use fastvep_core::ChromSynonyms;
use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_core::Position;
use noodles_csi::binning_index::BinningIndex;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

/// The variation file of a chromosome directory.
const DATA_FILE: &str = "all_vars.gz";

/// How many chromosome files stay open at once. VEP keeps five
/// (`VariationTabix::_get_tabix_obj`); a sorted VCF visits one contig at a
/// time, so a handful covers it, and a cache with thousands of scaffolds cannot
/// run the process out of file descriptors.
const MAX_OPEN_FILES: usize = 8;

/// Columns that are not population frequencies. Everything else in
/// `variation_cols` (`AFR`, `gnomADe_NFE`, ...) is a `POP -> "ALLELE:FREQ,..."`
/// column.
const NON_FREQUENCY_COLUMNS: [&str; 16] = [
    "chr",
    "variation_name",
    "failed",
    "somatic",
    "start",
    "end",
    "allele_string",
    "strand",
    "minor_allele",
    "minor_allele_freq",
    "clin_sig",
    "phenotype_or_disease",
    "pubmed",
    "clin_sig_allele",
    "clinical_impact",
    "var_synonyms",
];

/// A parsed record from the VEP variation cache tabix file.
#[derive(Debug, Clone)]
pub struct VariationCacheRecord {
    pub variation_name: String,
    pub failed: bool,
    pub somatic: bool,
    pub start: u64,
    pub end: u64,
    pub allele_string: String,
    pub strand: i8,
    pub minor_allele: Option<String>,
    pub minor_allele_freq: Option<f64>,
    pub clin_sig: Option<String>,
    pub phenotype_or_disease: bool,
    pub pubmed: Vec<String>,
    /// Population name → raw frequency string (e.g., "A:0.001,G:0.999")
    pub frequencies: HashMap<String, String>,
}

/// Where each field sits on a line, resolved once from `variation_cols`.
///
/// Parsing used to build a `HashMap<String, &str>` for every line it looked at,
/// cloning every column name each time. The names are fixed for the whole run.
#[derive(Debug, Clone)]
struct ColumnLayout {
    chr: usize,
    variation_name: Option<usize>,
    failed: Option<usize>,
    somatic: Option<usize>,
    start: usize,
    end: Option<usize>,
    allele_string: usize,
    strand: Option<usize>,
    minor_allele: Option<usize>,
    minor_allele_freq: Option<usize>,
    clin_sig: Option<usize>,
    phenotype_or_disease: Option<usize>,
    pubmed: Option<usize>,
    populations: Vec<(String, usize)>,
}

impl ColumnLayout {
    fn new(columns: &[String]) -> Result<Self> {
        let find = |name: &str| columns.iter().position(|c| c == name);
        let require = |name: &str| {
            find(name).ok_or_else(|| {
                anyhow!(
                    "info.txt `variation_cols` has no `{name}` column. A tabix variation cache \
                     starts with `chr` and carries `start` and `allele_string`; this is not the \
                     column set convert_cache.pl writes"
                )
            })
        };
        let populations = columns
            .iter()
            .enumerate()
            .filter(|(_, name)| !NON_FREQUENCY_COLUMNS.contains(&name.as_str()))
            .map(|(i, name)| (name.clone(), i))
            .collect();
        Ok(Self {
            chr: require("chr")?,
            variation_name: find("variation_name"),
            failed: find("failed"),
            somatic: find("somatic"),
            start: require("start")?,
            end: find("end"),
            allele_string: require("allele_string")?,
            strand: find("strand"),
            minor_allele: find("minor_allele"),
            minor_allele_freq: find("minor_allele_freq"),
            clin_sig: find("clin_sig"),
            phenotype_or_disease: find("phenotype_or_disease"),
            pubmed: find("pubmed"),
            populations,
        })
    }

    /// `start` alone, so a line outside the window is rejected without being
    /// parsed in full.
    fn start_of(&self, line: &str) -> Result<u64> {
        line.split('\t')
            .nth(self.start)
            .and_then(|v| v.parse().ok())
            .ok_or_else(|| anyhow!("variation cache line has no usable start: {:?}", line))
    }

    /// `parse_variation` (`BaseCacheVariation.pm`): `.` is no value, `failed`,
    /// `somatic` and `phenotype_or_disease` default to 0, `end` to `start`, and
    /// `strand` to 1.
    fn parse(&self, line: &str) -> Result<VariationCacheRecord> {
        let fields: Vec<&str> = line.split('\t').collect();
        let get = |col: Option<usize>| {
            col.and_then(|i| fields.get(i))
                .copied()
                .filter(|v| !v.is_empty() && *v != ".")
        };
        let flag = |col: Option<usize>| get(col).is_some_and(|v| v != "0");

        let start = self.start_of(line)?;
        // `||=`: a 0 end is no end either.
        let end = get(self.end)
            .and_then(|v| v.parse().ok())
            .filter(|&e| e != 0)
            .unwrap_or(start);

        let mut frequencies = HashMap::new();
        for (name, i) in &self.populations {
            if let Some(v) = get(Some(*i)) {
                frequencies.insert(name.clone(), v.to_string());
            }
        }

        Ok(VariationCacheRecord {
            variation_name: get(self.variation_name).unwrap_or_default().to_string(),
            failed: flag(self.failed),
            somatic: flag(self.somatic),
            start,
            end,
            allele_string: get(Some(self.allele_string))
                .unwrap_or_default()
                .to_string(),
            // `$v{strand} ||= 1`: a 0 strand is no strand either.
            strand: get(self.strand)
                .and_then(|v| v.parse().ok())
                .filter(|&s: &i8| s != 0)
                .unwrap_or(1),
            minor_allele: get(self.minor_allele).map(str::to_string),
            minor_allele_freq: get(self.minor_allele_freq).and_then(|v| v.parse().ok()),
            clin_sig: get(self.clin_sig).map(str::to_string),
            phenotype_or_disease: flag(self.phenotype_or_disease),
            pubmed: get(self.pubmed)
                .map(|v| v.split(',').map(str::to_string).collect())
                .unwrap_or_default(),
            frequencies,
        })
    }
}

/// One chromosome's `all_vars.gz`, open, with its index in memory.
struct Source {
    index: Box<dyn BinningIndex + Send + Sync>,
    /// Position of this chromosome in the index. `None` when the index lists no
    /// contig at all, which is what an empty file's index looks like.
    ref_id: Option<usize>,
    reader: bgzf::io::Reader<File>,
}

impl Source {
    fn open(dir: &str, data: &Path, index_path: &Path, layout: &ColumnLayout) -> Result<Self> {
        let index = read_index(index_path)?;

        let header = index.header().ok_or_else(|| {
            anyhow!(
                "{} has no tabix header, so it cannot be told which columns it indexes; \
                 rebuild it with `tabix -C -s {} -b {} -e {} {}`",
                index_path.display(),
                layout.chr + 1,
                layout.start + 1,
                layout.start + 1,
                data.display()
            )
        })?;

        // The index only finds a record if it was built on the column the
        // record's position is in. One built on another column answers every
        // query without error and returns the wrong records.
        if header.reference_sequence_name_index() != layout.chr
            || header.start_position_index() != layout.start
        {
            bail!(
                "{} indexes columns (contig {}, start {}) but info.txt puts `chr` at column {} and \
                 `start` at column {} (1-based). Rebuild it with `tabix -C -s {} -b {} -e {} {}`",
                index_path.display(),
                header.reference_sequence_name_index() + 1,
                header.start_position_index() + 1,
                layout.chr + 1,
                layout.start + 1,
                layout.chr + 1,
                layout.start + 1,
                layout.start + 1,
                data.display()
            );
        }

        let names = header.reference_sequence_names();
        let ref_id = if names.is_empty() {
            None
        } else {
            Some(names.iter().position(|name| name == dir).ok_or_else(|| {
                anyhow!(
                    "{} does not list contig `{}` (it lists {} other contig(s)); \
                             the file in directory `{}` holds another chromosome's records",
                    index_path.display(),
                    dir,
                    names.len(),
                    dir
                )
            })?)
        };

        let reader = File::open(data)
            .map(bgzf::io::Reader::new)
            .with_context(|| format!("Opening variation file: {}", data.display()))?;
        Ok(Self {
            index,
            ref_id,
            reader,
        })
    }

    /// The records whose start lies in `lo..=hi`, in file order.
    fn query(
        &mut self,
        lo: u64,
        hi: u64,
        layout: &ColumnLayout,
    ) -> Result<Vec<VariationCacheRecord>> {
        let Some(ref_id) = self.ref_id else {
            return Ok(Vec::new());
        };
        let lo = lo.max(1);
        if hi < lo {
            return Ok(Vec::new());
        }
        let interval: Interval = (position(lo)?..=position(hi)?).into();
        let chunks = self
            .index
            .query(ref_id, interval)
            .context("Tabix query failed")?;

        let mut records = Vec::new();
        let mut line = String::new();

        // The index hands back chunks sorted and merged (`optimize_chunks`), so no
        // record is read twice.
        for chunk in &chunks {
            self.reader.seek(chunk.start())?;
            loop {
                line.clear();
                if self.reader.read_line(&mut line)? == 0 {
                    break;
                }
                let text = line.trim_end_matches(['\n', '\r']);
                if !text.is_empty() {
                    // Sorted by start: the first record past the window ends the chunk.
                    let start = layout.start_of(text)?;
                    if start > hi {
                        break;
                    }
                    if start >= lo {
                        records.push(layout.parse(text)?);
                    }
                }
                if self.reader.virtual_position() >= chunk.end() {
                    break;
                }
            }
        }
        Ok(records)
    }
}

fn position(p: u64) -> Result<Position> {
    usize::try_from(p)
        .ok()
        .and_then(|p| Position::try_from(p).ok())
        .ok_or_else(|| anyhow!("position {p} is outside what an index can address"))
}

/// A `.csi` or a `.tbi`, by extension.
fn read_index(path: &Path) -> Result<Box<dyn BinningIndex + Send + Sync>> {
    let is_csi = path.extension().is_some_and(|e| e == "csi");
    let context = || format!("Reading index: {}", path.display());
    let index: Box<dyn BinningIndex + Send + Sync> = if is_csi {
        Box::new(noodles_csi::fs::read(path).with_context(context)?)
    } else {
        Box::new(noodles_tabix::fs::read(path).with_context(context)?)
    };
    Ok(index)
}

/// The index belonging to `data`: a `.csi` if there is one, else a `.tbi`.
///
/// That is the order the `tabix` command takes (checked with tabix 1.23.1: given
/// a file with both, it answered from the `.csi`), and VEP reads these files
/// through htslib. VEP's `convert_cache.pl` writes a `.csi` (`tabix -C`), so a
/// cache converted by VEP has no `.tbi`; having both is a hand-made directory.
fn find_index(data: &Path, start_column: usize) -> Result<PathBuf> {
    for extension in ["csi", "tbi"] {
        let mut name = data.as_os_str().to_owned();
        name.push(".");
        name.push(extension);
        let candidate = PathBuf::from(name);
        if candidate.is_file() {
            return Ok(candidate);
        }
    }
    bail!(
        "{} has no .csi or .tbi index, so it cannot be searched. Looking up nothing would read \
         exactly like looking up a variant that has no known neighbours. VEP's convert_cache.pl \
         writes this index with `tabix -C -s 1 -b {} -e {} {}`",
        data.display(),
        start_column + 1,
        start_column + 1,
        data.display()
    )
}

/// Open chromosome files, least recently used first.
struct OpenSources {
    cap: usize,
    order: Vec<String>,
    by_dir: HashMap<String, Source>,
}

impl Default for OpenSources {
    fn default() -> Self {
        Self::with_cap(MAX_OPEN_FILES)
    }
}

impl OpenSources {
    /// At least one file stays open, or there would be nothing to evict to.
    fn with_cap(cap: usize) -> Self {
        Self {
            cap: cap.max(1),
            order: Vec::new(),
            by_dir: HashMap::new(),
        }
    }

    fn get_or_open(
        &mut self,
        dir: &str,
        open: impl FnOnce() -> Result<Source>,
    ) -> Result<&mut Source> {
        if self.by_dir.contains_key(dir) {
            if self.order.last().map(String::as_str) != Some(dir) {
                self.order.retain(|d| d != dir);
                self.order.push(dir.to_string());
            }
        } else {
            let source = open()?;
            while self.order.len() >= self.cap {
                let evicted = self.order.remove(0);
                self.by_dir.remove(&evicted);
            }
            self.order.push(dir.to_string());
            self.by_dir.insert(dir.to_string(), source);
        }
        self.by_dir
            .get_mut(dir)
            .ok_or_else(|| anyhow!("chromosome file `{dir}` vanished from the open set"))
    }
}

/// Reader for VEP's tabix-indexed variation cache files.
///
/// Each chromosome has its own `all_vars.gz` and index under `{cache_dir}/{chr}/`.
pub struct VariationTabixReader {
    cache_dir: PathBuf,
    layout: ColumnLayout,
    /// Every accepted spelling of a contig → the name of its directory.
    contigs: HashMap<String, String>,
    /// Directory → index file, for the directories that hold an `all_vars.gz`.
    indexes: HashMap<String, PathBuf>,
    /// The directory names, sorted, for messages.
    chromosomes: Vec<String>,
    open: Mutex<OpenSources>,
}

impl VariationTabixReader {
    /// Open a VEP cache directory.
    ///
    /// `variation_cols` comes from `info.txt` and defines the column order.
    /// `valid_chromosomes` lists the chromosome directories; real caches do not
    /// carry that line (VEP takes the sub-directories, `CacheDir.pm`), so when
    /// it is empty the sub-directories are scanned.
    ///
    /// Fails, here and not once per variant later, when the directory holds no
    /// variation file or a variation file has no index: either would otherwise
    /// answer every lookup with "nothing known".
    ///
    /// A `chr_synonyms.txt` in the directory is read as VEP reads it
    /// (`CacheDir.pm` loads it automatically), so an input that names a contig
    /// differently from the directory (`chrM` and `MT`, or an accession) still
    /// finds it.
    pub fn new(
        cache_dir: &Path,
        variation_cols: &[String],
        valid_chromosomes: &[String],
    ) -> Result<Self> {
        let layout = ColumnLayout::new(variation_cols)?;

        let mut chromosomes: Vec<String> = if valid_chromosomes.is_empty() {
            let mut found = Vec::new();
            let entries = std::fs::read_dir(cache_dir)
                .with_context(|| format!("Reading cache directory: {}", cache_dir.display()))?;
            for entry in entries.flatten() {
                let is_dir = entry.file_type().map(|t| t.is_dir()).unwrap_or(false);
                if is_dir && entry.path().join(DATA_FILE).is_file() {
                    found.push(entry.file_name().to_string_lossy().to_string());
                }
            }
            found
        } else {
            valid_chromosomes.to_vec()
        };
        chromosomes.sort();
        chromosomes.dedup();

        let mut indexes = HashMap::new();
        for dir in &chromosomes {
            let data = cache_dir.join(dir).join(DATA_FILE);
            if data.is_file() {
                let index = find_index(&data, layout.start)?;
                // Open it once now and drop it. An index that cannot be read, or
                // was built on the wrong column, would otherwise surface on the
                // first variant of that chromosome, possibly hours into a run, or
                // not at all if the lookup error were ever swallowed. Indexes are
                // small (a CSI over a whole chromosome is about a megabyte).
                Source::open(dir, &data, &index, &layout)?;
                indexes.insert(dir.clone(), index);
            }
        }
        if indexes.is_empty() {
            bail!(
                "{} holds no `{{chr}}/{DATA_FILE}`. fastVEP reads the tabix form of VEP's variation \
                 cache (`var_type tabix` in info.txt), which is what the indexed cache download \
                 ships and what convert_cache.pl produces; it does not read the unconverted \
                 `{{start}}-{{end}}_var.gz` slices",
                cache_dir.display()
            );
        }

        let synonyms_path = cache_dir.join("chr_synonyms.txt");
        let synonyms = if synonyms_path.is_file() {
            let text = std::fs::read_to_string(&synonyms_path)
                .with_context(|| format!("Reading {}", synonyms_path.display()))?;
            ChromSynonyms::parse(&text)
        } else {
            ChromSynonyms::new()
        };
        let contigs = synonyms.alias_map(&chromosomes);

        Ok(Self {
            cache_dir: cache_dir.to_path_buf(),
            layout,
            contigs,
            indexes,
            chromosomes,
            open: Mutex::new(OpenSources::default()),
        })
    }

    /// The directory a contig name from the input resolves to, if the cache
    /// knows the contig at all.
    pub fn normalize_chrom(&self, chrom: &str) -> Option<&str> {
        self.contigs.get(chrom).map(String::as_str)
    }

    /// The chromosome directories of the cache, sorted.
    pub fn chromosomes(&self) -> &[String] {
        &self.chromosomes
    }

    /// The records whose `start` lies in `lo..=hi`.
    ///
    /// Not the records overlapping the window: `convert_cache.pl` indexes the
    /// start column alone (`-b N -e N`), so a record is found by its start, and
    /// the window VEP asks for is `start - 1 ..= end + 1` of the input
    /// (`VariationTabix::_annotate_pm`). A record's own `end` plays no part,
    /// which is what lets an insertion stored with `end < start` be found.
    ///
    /// A contig the cache does not know yields nothing; ask
    /// [`normalize_chrom`](Self::normalize_chrom) first to tell that from a
    /// window with nothing in it.
    pub fn query(&self, chrom: &str, lo: u64, hi: u64) -> Result<Vec<VariationCacheRecord>> {
        let Some(dir) = self.contigs.get(chrom) else {
            return Ok(Vec::new());
        };
        let Some(index_path) = self.indexes.get(dir) else {
            // The contig is known but the cache holds no variants for it.
            return Ok(Vec::new());
        };
        let mut open = self
            .open
            .lock()
            .map_err(|_| anyhow!("the variation cache reader's lock was poisoned"))?;
        let source = open.get_or_open(dir, || {
            Source::open(
                dir,
                &self.cache_dir.join(dir).join(DATA_FILE),
                index_path,
                &self.layout,
            )
        })?;
        source.query(lo, hi, &self.layout).with_context(|| {
            format!(
                "Reading {}/{DATA_FILE} in {}",
                dir,
                self.cache_dir.display()
            )
        })
    }

    /// Parse a single tab-separated line from the variation cache.
    pub fn parse_line(&self, line: &str) -> Result<VariationCacheRecord> {
        self.layout.parse(line)
    }
}

// -- matching ---------------------------------------------------------------

/// One input allele that a known variant matches.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlleleMatch {
    /// Position of the input allele in the alt list the site was built from.
    pub alt_index: usize,
    /// The allele as the cache record writes it (VEP's `b_allele`). Its
    /// per-population frequencies are keyed by this string.
    pub record_allele: String,
}

/// What `compare_existing` returns for a record that matches.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Matched {
    /// The record carries no alleles (`COSMIC_MUTATION`, `HGMD_MUTATION`): VEP
    /// reports it for every input allele when it has the same start and end,
    /// and keeps no `matched_alleles` for it.
    AllAlleles,
    Alleles(Vec<AlleleMatch>),
}

impl Matched {
    /// Whether the record is reported for the input allele at `alt_index`.
    pub fn covers(&self, alt_index: usize) -> bool {
        match self {
            Matched::AllAlleles => true,
            Matched::Alleles(matches) => matches.iter().any(|m| m.alt_index == alt_index),
        }
    }
}

/// A reference/alternate pair after `trim_sequences`. An allele that trimmed
/// away to nothing is the empty slice, and compares equal to `-`.
#[derive(Debug, Clone, Copy)]
struct Form<'a> {
    reference: &'a [u8],
    alternate: &'a [u8],
    pos: u64,
    alt_index: usize,
}

fn trim_right(reference: &mut &[u8], alternate: &mut &[u8]) {
    while let (Some((r, r_rest)), Some((a, a_rest))) =
        (reference.split_last(), alternate.split_last())
    {
        if !r.eq_ignore_ascii_case(a) {
            break;
        }
        *reference = r_rest;
        *alternate = a_rest;
    }
}

fn trim_left(reference: &mut &[u8], alternate: &mut &[u8], pos: &mut u64) {
    while let (Some((r, r_rest)), Some((a, a_rest))) =
        (reference.split_first(), alternate.split_first())
    {
        if !r.eq_ignore_ascii_case(a) {
            break;
        }
        *reference = r_rest;
        *alternate = a_rest;
        *pos += 1;
    }
}

/// `trim_sequences($ref, $alt, $pos, undef, 1, $right_first)` of
/// `Bio::EnsEMBL::Variation::Utils::Sequence`, as `get_matched_variant_alleles`
/// calls it: no strand, so only a left trim moves the position.
///
/// Case is ignored. VEP compares as written but uppercases the input
/// (`Parser::validate_vf`) and the cache is uppercase.
fn trim<'a>(
    reference: &'a [u8],
    alternate: &'a [u8],
    pos: u64,
    right_first: bool,
) -> (&'a [u8], &'a [u8], u64) {
    let (mut r, mut a, mut pos) = (reference, alternate, pos);
    if right_first {
        trim_right(&mut r, &mut a);
        trim_left(&mut r, &mut a, &mut pos);
    } else {
        trim_left(&mut r, &mut a, &mut pos);
        trim_right(&mut r, &mut a);
    }
    (r, a, pos)
}

/// `_get_trim_directions`: both orders when either allele is longer than one
/// base, since the two can land on different positions; one otherwise.
fn trim_orders(reference: &str, alternate: &str) -> &'static [bool] {
    if reference.len() > 1 || alternate.len() > 1 {
        &[false, true]
    } else {
        &[false]
    }
}

/// An allele that trimmed to nothing is `-`, as VEP's `$empty_to_dash` has it.
fn allele_eq(a: &[u8], b: &[u8]) -> bool {
    let a = if a.is_empty() { b"-" } else { a };
    let b = if b.is_empty() { b"-" } else { b };
    a.eq_ignore_ascii_case(b)
}

/// `reverse_comp` of `Bio::EnsEMBL::Utils::Sequence`.
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'a' => 't',
            'c' => 'g',
            'g' => 'c',
            't' => 'a',
            'r' => 'y',
            'y' => 'r',
            'm' => 'k',
            'k' => 'm',
            'w' => 'w',
            's' => 's',
            'b' => 'v',
            'd' => 'h',
            'h' => 'd',
            'v' => 'b',
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            'R' => 'Y',
            'Y' => 'R',
            'M' => 'K',
            'K' => 'M',
            'B' => 'V',
            'D' => 'H',
            'H' => 'D',
            'V' => 'B',
            other => other,
        })
        .collect()
}

/// Each of the record's alt alleles against the site's reduced forms.
///
/// `get_matched_variant_alleles`, second half: every (ref, alt) pair of the
/// record is trimmed both ways and looked up by (ref, alt, position). The
/// *last* form of the site with that key wins, because VEP stores them in a
/// hash and a later pair overwrites an earlier one.
fn match_record_side(
    forms: &[Form<'_>],
    record_allele_string: &str,
    record_start: u64,
) -> Vec<AlleleMatch> {
    let mut alleles = record_allele_string.split('/');
    let record_ref = alleles.next().unwrap_or("");
    let mut matches = Vec::new();
    for record_alt in alleles {
        for &right_first in trim_orders(record_ref, record_alt) {
            let (r, a, pos) = trim(
                record_ref.as_bytes(),
                record_alt.as_bytes(),
                record_start,
                right_first,
            );
            let hit = forms
                .iter()
                .rev()
                .find(|f| f.pos == pos && allele_eq(f.reference, r) && allele_eq(f.alternate, a));
            if let Some(form) = hit {
                matches.push(AlleleMatch {
                    alt_index: form.alt_index,
                    record_allele: record_alt.to_string(),
                });
                break;
            }
        }
    }
    matches
}

/// An input site, reduced once so that many candidate records can be compared
/// against it without redoing the input's half of the work.
pub struct SiteForms<'a> {
    ref_allele: &'a str,
    alts: Vec<&'a str>,
    start: u64,
    end: u64,
    forms: Vec<Form<'a>>,
}

impl<'a> SiteForms<'a> {
    /// `ref_allele` and `alts` as VEP writes an allele string (`-` for an empty
    /// allele); `start`/`end` are Ensembl coordinates, so an insertion has
    /// `end == start - 1`.
    pub fn new(ref_allele: &'a str, alts: &[&'a str], start: u64, end: u64) -> Self {
        let mut forms = Vec::with_capacity(alts.len() * 2);
        for (alt_index, alt) in alts.iter().enumerate() {
            for &right_first in trim_orders(ref_allele, alt) {
                let (reference, alternate, pos) =
                    trim(ref_allele.as_bytes(), alt.as_bytes(), start, right_first);
                forms.push(Form {
                    reference,
                    alternate,
                    pos,
                    alt_index,
                });
            }
        }
        Self {
            ref_allele,
            alts: alts.to_vec(),
            start,
            end,
            forms,
        }
    }

    /// `get_matched_variant_alleles` for this site against one record's allele
    /// string, start and strand.
    ///
    /// A record on the other strand is compared against the input's reverse
    /// complement with the input's position unchanged, which is what VEP does
    /// (`$a_ref` and `$rev_alt` are reversed, `pos` is not).
    pub fn match_alleles(
        &self,
        record_allele_string: &str,
        record_start: u64,
        record_strand: i8,
    ) -> Vec<AlleleMatch> {
        if record_strand == 1 {
            return match_record_side(&self.forms, record_allele_string, record_start);
        }
        let reference = reverse_complement(self.ref_allele);
        let alts: Vec<String> = self.alts.iter().map(|a| reverse_complement(a)).collect();
        let mut forms = Vec::with_capacity(alts.len() * 2);
        for (alt_index, alt) in alts.iter().enumerate() {
            for &right_first in trim_orders(&reference, alt) {
                let (r, a, pos) = trim(
                    reference.as_bytes(),
                    alt.as_bytes(),
                    self.start,
                    right_first,
                );
                forms.push(Form {
                    reference: r,
                    alternate: a,
                    pos,
                    alt_index,
                });
            }
        }
        match_record_side(&forms, record_allele_string, record_start)
    }

    /// `compare_existing` (`AnnotationType/Variation.pm`): does this record
    /// count as a known variant of the site, and for which alleles?
    ///
    /// A record without a `/` in its allele string carries no alleles to
    /// compare. VEP matches it on identical start and end and nothing else,
    /// and when they differ it does not match: the fall-through to
    /// `get_matched_variant_alleles` has no alt to try.
    pub fn compare(&self, record: &VariationCacheRecord) -> Option<Matched> {
        if !record.allele_string.contains('/') {
            return (record.start == self.start && record.end == self.end)
                .then_some(Matched::AllAlleles);
        }
        let matches = self.match_alleles(&record.allele_string, record.start, record.strand);
        (!matches.is_empty()).then_some(Matched::Alleles(matches))
    }
}

/// VEP's `get_matched_variant_alleles` for an input site against one known
/// variant, given the way its tests and callers give it: reference and
/// alternates of the input, and the record's allele string, start and strand.
pub fn match_alleles(
    ref_allele: &str,
    alts: &[&str],
    start: u64,
    record_allele_string: &str,
    record_start: u64,
    record_strand: i8,
) -> Vec<AlleleMatch> {
    SiteForms::new(ref_allele, alts, start, start).match_alleles(
        record_allele_string,
        record_start,
        record_strand,
    )
}

/// How VEP orders the known variants of one allele in `Existing_variation`
/// (`OutputFactory::add_colocated_variant_info`): germline before somatic, then
/// by the prefix of the name (`rs`, then HGMD's `CM`/`CI`/`CD`, then COSMIC's
/// `CO`, then everything else). The sort is stable, so ties keep file order.
pub fn existing_variation_rank(name: &str, somatic: bool) -> (bool, u8) {
    let prefix = name.as_bytes().get(..2).unwrap_or(&[]);
    let rank = if prefix.eq_ignore_ascii_case(b"rs") {
        1
    } else if ["cm", "ci", "cd"]
        .iter()
        .any(|p| prefix.eq_ignore_ascii_case(p.as_bytes()))
    {
        2
    } else if prefix.eq_ignore_ascii_case(b"co") {
        3
    } else {
        100
    };
    (somatic, rank)
}

/// Parse a VEP frequency string like "A:0.001,G:0.999" into allele-frequency pairs.
pub fn parse_freq_string(freq_str: &str) -> Vec<(String, f64)> {
    let mut result = Vec::new();
    for pair in freq_str.split(',') {
        if let Some((allele, freq_s)) = pair.split_once(':') {
            if let Ok(freq) = freq_s.parse::<f64>() {
                result.push((allele.to_string(), freq));
            }
        }
    }
    result
}

/// Extract the frequency for a specific allele from a population's frequency string.
///
/// Given a freq string like "A:0.001,G:0.999" and allele "A", returns Some(0.001).
pub fn get_allele_freq(freq_str: &str, allele: &str) -> Option<f64> {
    for pair in freq_str.split(',') {
        if let Some((a, f)) = pair.split_once(':') {
            if a.eq_ignore_ascii_case(allele) {
                return f.parse().ok();
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(allele_string: &str, start: u64, end: u64, strand: i8) -> VariationCacheRecord {
        VariationCacheRecord {
            variation_name: "rs1".into(),
            failed: false,
            somatic: false,
            start,
            end,
            allele_string: allele_string.into(),
            strand,
            minor_allele: None,
            minor_allele_freq: None,
            clin_sig: None,
            phenotype_or_disease: false,
            pubmed: vec![],
            frequencies: HashMap::new(),
        }
    }

    /// A reader over no files, for exercising line parsing alone.
    fn parsing_reader(columns: &[String]) -> Result<VariationTabixReader> {
        Ok(VariationTabixReader {
            cache_dir: PathBuf::new(),
            layout: ColumnLayout::new(columns)?,
            contigs: HashMap::new(),
            indexes: HashMap::new(),
            chromosomes: Vec::new(),
            open: Mutex::new(OpenSources::default()),
        })
    }

    /// A reader that may keep one file open still has to answer for every
    /// chromosome: each change of chromosome evicts the open file and reopens
    /// the other. The cap only bites with more than `MAX_OPEN_FILES`
    /// chromosomes, which no fixture here has, so lower it.
    #[test]
    fn eviction_and_reopening_keep_every_answer() {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/vep_cache/csi");
        let info = crate::info::CacheInfo::from_file(&dir.join("info.txt")).unwrap();
        let reader =
            VariationTabixReader::new(&dir, &info.variation_cols, &info.valid_chromosomes).unwrap();
        *reader.open.lock().unwrap() = OpenSources::with_cap(1);

        for _ in 0..3 {
            let on_10 = reader.query("chr10", 10401, 10401).unwrap();
            assert_eq!(on_10.len(), 1);
            assert_eq!(on_10[0].variation_name, "rsFIX_DEL_ANCH");
            let on_mt = reader.query("chrM", 100, 100).unwrap();
            assert_eq!(on_mt.len(), 1);
            assert_eq!(on_mt[0].variation_name, "rsFIX_MT");
            assert_eq!(
                reader.open.lock().unwrap().by_dir.len(),
                1,
                "the cap of one holds"
            );
        }
    }

    fn alt_indexes(matches: &[AlleleMatch]) -> Vec<(usize, &str)> {
        matches
            .iter()
            .map(|m| (m.alt_index, m.record_allele.as_str()))
            .collect()
    }

    #[test]
    fn test_parse_freq_string() {
        let freqs = parse_freq_string("A:0.001,G:0.999");
        assert_eq!(freqs.len(), 2);
        assert_eq!(freqs[0], ("A".to_string(), 0.001));
        assert_eq!(freqs[1], ("G".to_string(), 0.999));

        // Empty/dot
        let freqs = parse_freq_string("");
        assert!(freqs.is_empty());
    }

    #[test]
    fn test_get_allele_freq() {
        assert_eq!(get_allele_freq("A:0.001,G:0.999", "A"), Some(0.001));
        assert_eq!(get_allele_freq("A:0.001,G:0.999", "G"), Some(0.999));
        assert_eq!(get_allele_freq("A:0.001,G:0.999", "T"), None);
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ACGT"), "ACGT");
        assert_eq!(reverse_complement("A"), "T");
        assert_eq!(reverse_complement("ATG"), "CAT");
        // The empty-allele marker passes through, as it does in the Perl.
        assert_eq!(reverse_complement("-"), "-");
    }

    #[test]
    fn test_parse_line() {
        let cols: Vec<String> = "chr,variation_name,failed,somatic,start,end,allele_string,strand,minor_allele,minor_allele_freq,clin_sig,phenotype_or_disease,pubmed,AFR,AMR,EAS,EUR,SAS,AA,EA"
            .split(',')
            .map(String::from)
            .collect();

        let reader = parsing_reader(&cols).expect("a complete column set parses");

        let line = "21\trs559462325\t.\t.\t8522406\t.\tG/A\t.\tA\t0.0002\t.\t.\t.\tA:0\tA:0\tA:0.001\tA:0\tA:0\t.\t.";
        let record = reader.parse_line(line).unwrap();

        assert_eq!(record.variation_name, "rs559462325");
        assert!(!record.failed);
        assert!(!record.somatic);
        assert_eq!(record.start, 8522406);
        assert_eq!(record.end, 8522406); // end defaults to start when "."
        assert_eq!(record.allele_string, "G/A");
        assert_eq!(record.strand, 1); // defaults to 1 when "."
        assert_eq!(record.minor_allele, Some("A".into()));
        assert_eq!(record.minor_allele_freq, Some(0.0002));
        assert!(record.clin_sig.is_none());
        assert!(!record.phenotype_or_disease);

        // Check population frequencies
        assert_eq!(record.frequencies.get("AFR"), Some(&"A:0".to_string()));
        assert_eq!(record.frequencies.get("EAS"), Some(&"A:0.001".to_string()));
    }

    #[test]
    fn a_zero_end_is_no_end() {
        let cols: Vec<String> = "chr,variation_name,start,end,allele_string"
            .split(',')
            .map(String::from)
            .collect();
        let reader = parsing_reader(&cols).unwrap();
        let r = reader.parse_line("1\trs1\t100\t0\tA/G").unwrap();
        assert_eq!(r.end, 100, "`$v{{end}} ||= $v{{start}}` treats 0 as unset");
        // An insertion keeps the end below its start.
        let r = reader.parse_line("1\trs2\t100\t99\t-/T").unwrap();
        assert_eq!((r.start, r.end), (100, 99));
    }

    #[test]
    fn a_line_without_a_start_is_an_error_not_a_record_at_zero() {
        let cols: Vec<String> = "chr,variation_name,start,end,allele_string"
            .split(',')
            .map(String::from)
            .collect();
        let reader = parsing_reader(&cols).unwrap();
        assert!(reader.parse_line("1\trs1\t.\t.\tA/G").is_err());
        assert!(reader.parse_line("1\trs1").is_err());
    }

    #[test]
    fn column_sets_without_the_columns_matching_needs_are_refused() {
        let no_chr: Vec<String> = ["variation_name", "start", "allele_string"]
            .map(String::from)
            .to_vec();
        let err = ColumnLayout::new(&no_chr).unwrap_err().to_string();
        assert!(err.contains("`chr`"), "{err}");
        let no_alleles: Vec<String> = ["chr", "variation_name", "start"]
            .map(String::from)
            .to_vec();
        let err = ColumnLayout::new(&no_alleles).unwrap_err().to_string();
        assert!(err.contains("`allele_string`"), "{err}");
    }

    #[test]
    fn test_match_alleles_snv() {
        let r = record("C/T", 100, 100, 1);
        let site = SiteForms::new("C", &["T"], 100, 100);
        assert_eq!(
            alt_indexes(&site.match_alleles(&r.allele_string, 100, 1)),
            [(0, "T")]
        );
        assert!(SiteForms::new("C", &["G"], 100, 100)
            .match_alleles(&r.allele_string, 100, 1)
            .is_empty());
        // Another position: the key includes it.
        assert!(SiteForms::new("C", &["T"], 200, 200)
            .match_alleles(&r.allele_string, 100, 1)
            .is_empty());
    }

    #[test]
    fn test_match_alleles_reverse_strand() {
        // A record on the minus strand, G/A, is C/T on the plus strand.
        let m = match_alleles("C", &["T"], 100, "G/A", 100, -1);
        assert_eq!(alt_indexes(&m), [(0, "A")]);
    }

    #[test]
    fn test_match_alleles_multiallelic() {
        assert_eq!(
            alt_indexes(&match_alleles("C", &["T"], 100, "C/T/G", 100, 1)),
            [(0, "T")]
        );
        assert_eq!(
            alt_indexes(&match_alleles("C", &["G"], 100, "C/T/G", 100, 1)),
            [(0, "G")]
        );
        assert!(match_alleles("C", &["A"], 100, "C/T/G", 100, 1).is_empty());
        // Both input alts are reported, each with its own record allele.
        let m = match_alleles("C", &["T", "G"], 100, "C/T/G", 100, 1);
        assert_eq!(alt_indexes(&m), [(0, "T"), (1, "G")]);
    }

    /// Issue #139, claim 5. The input `AACCCT>A` is minimised to `ACCCT/-` at
    /// 10401; the cache writes the same deletion anchored, `ACCCTA/A` at 10401.
    /// VEP trims the record from the right and gets `ACCCT/-` at 10401. The
    /// oracle (VEP's own `get_matched_variant_alleles`) says the same.
    #[test]
    fn an_anchored_deletion_record_matches_the_minimised_input() {
        let m = match_alleles("ACCCT", &["-"], 10401, "ACCCTA/A", 10401, 1);
        assert_eq!(alt_indexes(&m), [(0, "A")]);
        // And a minimal record that is a different deletion does not match.
        assert!(match_alleles("ACCCT", &["-"], 10401, "AACCC/-", 10400, 1).is_empty());
    }

    #[test]
    fn an_anchored_insertion_record_matches_the_minimised_input() {
        // Input `G>GA` at 10420 is `-/A` at 10421; the record is `G/GA` at 10420.
        let m = match_alleles("-", &["A"], 10421, "G/GA", 10420, 1);
        assert_eq!(alt_indexes(&m), [(0, "GA")]);
        // The minimal record `-/A` at 10421 matches itself.
        assert_eq!(
            alt_indexes(&match_alleles("-", &["A"], 10421, "-/A", 10421, 1)),
            [(0, "A")]
        );
        // An insertion one base along is another insertion.
        assert!(match_alleles("-", &["A"], 10421, "-/A", 10422, 1).is_empty());
    }

    #[test]
    fn a_record_without_alleles_matches_on_coordinates_alone() {
        let site = SiteForms::new("C", &["A", "G"], 10414, 10414);
        let cosmic = record("COSMIC_MUTATION", 10414, 10414, 1);
        assert_eq!(site.compare(&cosmic), Some(Matched::AllAlleles));
        assert!(Matched::AllAlleles.covers(0) && Matched::AllAlleles.covers(1));
        // A different span does not match, whatever the alleles are.
        assert_eq!(
            site.compare(&record("COSMIC_MUTATION", 10414, 10415, 1)),
            None
        );
        assert_eq!(
            site.compare(&record("COSMIC_MUTATION", 10413, 10413, 1)),
            None
        );
        // An empty allele string is the same thing: VEP's `!~ /\//` is true for it.
        assert_eq!(
            site.compare(&record("", 10414, 10414, 1)),
            Some(Matched::AllAlleles)
        );
    }

    #[test]
    fn compare_reports_which_alleles_matched() {
        let site = SiteForms::new("C", &["A", "G"], 100, 100);
        let matched = site.compare(&record("C/G", 100, 100, 1)).unwrap();
        assert!(!matched.covers(0));
        assert!(matched.covers(1));
        assert_eq!(site.compare(&record("C/T", 100, 100, 1)), None);
    }

    #[test]
    fn existing_variants_sort_the_way_vep_lists_them() {
        let mut names = [
            ("COSV1", true),
            ("CM001", false),
            ("zzz", false),
            ("rs2", false),
            ("co9", false),
            ("rs1", true),
            ("RS3", false),
        ];
        names.sort_by_key(|(name, somatic)| existing_variation_rank(name, *somatic));
        let order: Vec<&str> = names.iter().map(|(n, _)| *n).collect();
        // germline first (rs, then CM, then CO, then the rest), somatic after.
        assert_eq!(order, ["rs2", "RS3", "CM001", "co9", "zzz", "rs1", "COSV1"]);
    }

    #[test]
    fn test_tabix_query_real_data() {
        // Use real VEP cache test data
        let cache_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("test_data/vep_cache/homo_sapiens/84_GRCh38");

        if !cache_dir.join("info.txt").exists() {
            eprintln!("Skipping tabix test: no test data at {:?}", cache_dir);
            return;
        }

        let info = crate::info::CacheInfo::from_file(&cache_dir.join("info.txt")).unwrap();

        let reader =
            VariationTabixReader::new(&cache_dir, &info.variation_cols, &info.valid_chromosomes)
                .unwrap();

        // Query for a known variant: rs879182429 at chr21:8013350
        let records = reader.query("21", 8013350, 8013350).unwrap();
        assert!(
            !records.is_empty(),
            "Expected at least one record at 21:8013350"
        );

        let rs = records.iter().find(|r| r.variation_name == "rs879182429");
        assert!(rs.is_some(), "Expected rs879182429 in results");
        let rs = rs.unwrap();
        assert_eq!(rs.allele_string, "T/G");
        assert_eq!(rs.start, 8013350);

        // Query for a variant with frequency data: rs559462325 at chr21:8522406
        let records = reader.query("21", 8522406, 8522406).unwrap();
        let rs = records.iter().find(|r| r.variation_name == "rs559462325");
        assert!(rs.is_some(), "Expected rs559462325 in results");
        let rs = rs.unwrap();
        assert_eq!(rs.minor_allele, Some("A".to_string()));
        assert_eq!(rs.minor_allele_freq, Some(0.0002));
        assert!(rs.frequencies.contains_key("AFR"));

        // Query with "chr" prefix should also work
        let records = reader.query("chr21", 8013350, 8013350).unwrap();
        assert!(!records.is_empty(), "chr21 prefix should normalize to 21");

        // Query for a non-existent region should return empty
        let records = reader.query("21", 1, 1).unwrap();
        assert!(records.is_empty(), "No variants expected at position 1");
    }
}
