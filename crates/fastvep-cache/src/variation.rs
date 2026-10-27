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
use std::collections::{HashMap, VecDeque};
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
const NON_FREQUENCY_COLUMNS: [&str; 17] = [
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
    "clin_sig_ref_allele",
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
    /// `T:pathogenic;C:benign`: the significance per allele, which VEP prefers
    /// to `clin_sig` where the record has it (`add_colocated_variant_info`).
    pub clin_sig_allele: Option<String>,
    /// The reference allele `clin_sig_allele` is written against. When it is not
    /// the input's reference, VEP reverse-complements the alleles of
    /// `clin_sig_allele` before it looks the input allele up in them.
    pub clin_sig_ref_allele: Option<String>,
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
    clin_sig_allele: Option<usize>,
    clin_sig_ref_allele: Option<usize>,
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
            clin_sig_allele: find("clin_sig_allele"),
            clin_sig_ref_allele: find("clin_sig_ref_allele"),
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
            clin_sig_allele: get(self.clin_sig_allele).map(str::to_string),
            clin_sig_ref_allele: get(self.clin_sig_ref_allele).map(str::to_string),
            phenotype_or_disease: flag(self.phenotype_or_disease),
            // `split(',', $ex->{pubmed})`: Perl drops the empty fields at the end, and the
            // real cache has a trailing comma on a few records (two of the 105,494 of a
            // whole-genome callset on chr21 and chr22 had one in their PUBMED list)
            pubmed: get(self.pubmed)
                .map(|v| {
                    let mut ids: Vec<String> = v.split(',').map(str::to_string).collect();
                    while ids.last().is_some_and(String::is_empty) {
                        ids.pop();
                    }
                    ids
                })
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
    /// Lines read for the last query and not yet past, in file order, each with its
    /// start. A sorted input asks for a window that overlaps the one before it, so
    /// these are served again and the reader goes on from where it stopped,
    /// instead of seeking back to the index's chunk start (which inflates the
    /// block again and splits every line up to `lo` again, for every variant).
    window: VecDeque<(u64, String)>,
    /// Every line of the file that was read and has a start of at least this is in
    /// `window`.
    floor: u64,
    /// `reader` sits right after the last line of `window`.
    positioned: bool,
    /// `reader` has read the last line of the file.
    at_eof: bool,
    /// Lines read from the file, for a test that holds a sorted sweep to one read of each.
    lines_read: u64,
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
            window: VecDeque::new(),
            floor: 0,
            positioned: false,
            at_eof: false,
            lines_read: 0,
        })
    }

    /// The records whose start lies in `lo..=hi`, in file order.
    ///
    /// The file is sorted by start, so what a query needs is the run of lines from
    /// the first with a start of at least `lo` to the first past `hi`. The index
    /// only says where a run that long can begin (`chunks[0].start()`), and a
    /// window that follows the one before it needs no index at all: the lines are
    /// the ones already read, and the reader carries on from the last of them.
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
        // The earliest place a record of the window can be.
        let Some(first) = chunks.first().map(|c| c.start()) else {
            return Ok(Vec::new());
        };

        // Carry on from the last query when everything this one needs that is
        // before the reader's position is still in the window: it starts at or after
        // what the window holds, and the index puts its first line no later than
        // where the reader is. Otherwise (an unsorted input, or a jump forward past
        // what is worth reading through) go back to the index.
        let carry_on =
            self.positioned && lo >= self.floor && first <= self.reader.virtual_position();
        // Not trusted again until this query has read to its end: an error half way
        // leaves the window and the reader out of step.
        self.positioned = false;
        if !carry_on {
            self.window.clear();
            self.reader.seek(first)?;
            self.at_eof = false;
        }
        self.floor = lo;
        while self.window.front().is_some_and(|(start, _)| *start < lo) {
            self.window.pop_front();
        }

        // Read on until a line is past `hi`, keeping the lines that are not before `lo`.
        let mut line = String::new();
        while !self.at_eof && self.window.back().is_none_or(|(start, _)| *start <= hi) {
            line.clear();
            if self.reader.read_line(&mut line)? == 0 {
                self.at_eof = true;
                break;
            }
            self.lines_read += 1;
            let text = line.trim_end_matches(['\n', '\r']);
            if text.is_empty() {
                continue;
            }
            let start = layout.start_of(text)?;
            if start >= lo {
                self.window.push_back((start, text.to_string()));
            }
        }

        self.positioned = true;
        let mut records = Vec::new();
        for (start, text) in &self.window {
            if *start > hi {
                break;
            }
            records.push(layout.parse(text)?);
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
    /// The chromosomes are the sub-directories that hold an `all_vars.gz`, as in
    /// VEP, which takes them from the file system and overwrites whatever
    /// `info.txt` says (`CacheDir::info`); a `valid_chromosomes` line there is not
    /// consulted.
    ///
    /// Fails, here and not once per variant later, when the directory holds no
    /// variation file or a variation file has no index: either would otherwise
    /// answer every lookup with "nothing known".
    ///
    /// A `chr_synonyms.txt` in the directory is read as VEP reads it
    /// (`CacheDir.pm` loads it automatically), so an input that names a contig
    /// differently from the directory (`chrM` and `MT`, or an accession) still
    /// finds it.
    pub fn new(cache_dir: &Path, variation_cols: &[String]) -> Result<Self> {
        let layout = ColumnLayout::new(variation_cols)?;

        let mut chromosomes = Vec::new();
        let entries = std::fs::read_dir(cache_dir)
            .with_context(|| format!("Reading cache directory: {}", cache_dir.display()))?;
        for entry in entries.flatten() {
            let is_dir = entry.file_type().map(|t| t.is_dir()).unwrap_or(false);
            if is_dir && entry.path().join(DATA_FILE).is_file() {
                chromosomes.push(entry.file_name().to_string_lossy().to_string());
            }
        }
        chromosomes.sort();

        let mut indexes = HashMap::new();
        for dir in &chromosomes {
            let data = cache_dir.join(dir).join(DATA_FILE);
            if data.is_file() {
                let index = find_index(&data, layout.start)?;
                // Open it once now and drop it. An index that cannot be read, or
                // was built on the wrong column, would otherwise surface on the
                // first variant of that chromosome, possibly hours into a run, or
                // not at all if the lookup error were ever swallowed. An index is
                // small next to the data file it indexes.
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
pub fn reverse_complement(seq: &str) -> String {
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

/// One input variant as VEP hands it to `compare_existing`: a reference, the
/// alts that share it, and the coordinates they sit at.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InputSite {
    pub reference: String,
    pub alts: Vec<String>,
    pub start: u64,
    pub end: u64,
}

/// `trim_sequences($ref, $alt, $start, $end, 1)`: common bases off the left,
/// moving `start`, then off the right, moving `end`; an allele left empty is
/// `-`.
fn minimise(reference: &str, alternate: &str, start: u64, end: u64) -> (String, String, u64, u64) {
    fn bases(allele: &str) -> &[u8] {
        if allele == "-" {
            &[]
        } else {
            allele.as_bytes()
        }
    }
    let (mut r, mut a) = (bases(reference), bases(alternate));
    let (mut start, mut end) = (start, end);
    while let (Some((x, r_rest)), Some((y, a_rest))) = (r.split_first(), a.split_first()) {
        if !x.eq_ignore_ascii_case(y) {
            break;
        }
        (r, a) = (r_rest, a_rest);
        start += 1;
    }
    while let (Some((x, r_rest)), Some((y, a_rest))) = (r.split_last(), a.split_last()) {
        if !x.eq_ignore_ascii_case(y) {
            break;
        }
        (r, a) = (r_rest, a_rest);
        end = end.saturating_sub(1);
    }
    let text = |s: &[u8]| {
        if s.is_empty() {
            "-".to_string()
        } else {
            String::from_utf8_lossy(s).into_owned()
        }
    };
    (text(r), text(a), start, end)
}

/// The site VEP looks known variants up for, given one VCF record as its
/// parser leaves it (the first base already taken off where `Parser::VCF`
/// does).
///
/// A record with one alt whose length differs from the reference is what VEP
/// calls a non-minimised indel (`Parser::post_process_vfs`), and it replaces
/// the record with the pair trimmed on both sides, left first and then right
/// (`Parser::minimise_alleles`), before anything is looked up. That is not a
/// detail. `get_matched_variant_alleles` tries both trim orders on what it is
/// given, so a pair that is already minimal has only one reading, and an event
/// written in a repeat (`AC>ACCC`, the insertion of `CC` after the `C`) no
/// longer matches the same insertion the cache writes one base to the left.
/// Matching the record as written found a known variant that VEP reports as
/// novel.
///
/// Every other record goes through as it is: SNVs and MNVs, and records with
/// several alts. VEP 115 split a multi-allelic indel into one variant per alt
/// and minimised each (`InputBuffer::split_variants`, switched on for the whole
/// run by the first indel it met), and then lost the known variants of every
/// alt but the first when it rejoined them. VEP 116 no longer switches that on
/// (`InputBuffer::next`, the lines are commented out), so it matches the whole
/// record, and so does this.
pub fn vep_input_site(ref_allele: &str, alts: &[&str], start: u64, end: u64) -> InputSite {
    let length = |s: &str| if s == "-" { 0 } else { s.len() };
    if let [alt] = alts {
        if length(alt) != length(ref_allele) {
            let (reference, alt, start, end) = minimise(ref_allele, alt, start, end);
            return InputSite {
                reference,
                alts: vec![alt],
                start,
                end,
            };
        }
    }
    InputSite {
        reference: ref_allele.to_string(),
        alts: alts.iter().map(|a| a.to_string()).collect(),
        start,
        end,
    }
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
            clin_sig_allele: None,
            clin_sig_ref_allele: None,
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
    /// The cache of `fixtures/vep_cache_window`: 10,040 records on one chromosome in the densities a
    /// real one has (a record every few bases, three hundred at one position, a sparse stretch, a
    /// gap of a megabase) and some that are indels with an end before their start.
    fn window_cache() -> (PathBuf, Vec<String>) {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/vep_cache_window");
        let info = crate::info::CacheInfo::from_file(&dir.join("info.txt")).unwrap();
        (dir, info.variation_cols)
    }

    /// Every record of the file as `(start, name)`, read with no index at all.
    fn window_cache_records() -> Vec<(u64, String)> {
        let (dir, cols) = window_cache();
        let layout = ColumnLayout::new(&cols).unwrap();
        let mut reader = bgzf::io::Reader::new(File::open(dir.join("7/all_vars.gz")).unwrap());
        let mut all = Vec::new();
        let mut line = String::new();
        while reader.read_line(&mut line).unwrap() > 0 {
            let text = line.trim_end();
            let record = layout.parse(text).unwrap();
            all.push((layout.start_of(text).unwrap(), record.variation_name));
            line.clear();
        }
        all
    }

    fn answer(records: &[VariationCacheRecord]) -> Vec<(u64, String)> {
        records
            .iter()
            .map(|r| (r.start, r.variation_name.clone()))
            .collect()
    }

    /// A reader that carries on from the last query answers every window as one that starts from
    /// the index each time, and as a scan of the whole file does, whatever order the windows come
    /// in: a sweep, a sweep with jumps back and forward, the same window again, windows that
    /// start before the first record, in the gap, at the hot site and past the last.
    #[test]
    fn a_reader_that_carries_on_from_the_last_query_answers_as_one_that_does_not() {
        let (dir, cols) = window_cache();
        let all = window_cache_records();
        assert_eq!(all.len(), 10_040);
        let by_scan = |lo: u64, hi: u64| -> Vec<(u64, String)> {
            all.iter()
                .filter(|(s, _)| *s >= lo.max(1) && *s <= hi)
                .cloned()
                .collect()
        };
        let carrying = VariationTabixReader::new(&dir, &cols).unwrap();

        // a deterministic stream of windows
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        let mut next = |n: u64| -> u64 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 33) % n
        };
        let mut windows: Vec<(u64, u64)> = Vec::new();
        // a sweep through the dense region with the windows an input would ask for
        let mut p = 99_990;
        for _ in 0..3_000 {
            p += next(12);
            windows.push((p.saturating_sub(1), p + next(30) + 1));
        }
        // the same windows again, and with jumps among them
        windows.extend_from_slice(&windows.clone()[..200]);
        for _ in 0..1_500 {
            let at = 99_000 + next(1_600_000);
            windows.push((at, at + next(200)));
        }
        // the hot site, the sparse stretch, the gap, before the first record, past the last
        for (lo, hi) in [
            (140_000, 140_000),
            (139_999, 140_001),
            (140_000, 140_000),
            (199_990, 400_000),
            (300_000, 300_100),
            (700_000, 1_400_000),
            (1, 50),
            (0, 10),
            (1_000_000_000, 1_000_000_100),
            (99_000, 99_001),
            (1_500_000, 1_900_000),
        ] {
            windows.push((lo, hi));
        }
        // and a long sorted sweep through the second dense region
        let mut p = 1_499_990;
        for _ in 0..2_500 {
            p += next(3);
            windows.push((p, p + next(8)));
        }

        for (i, (lo, hi)) in windows.iter().enumerate() {
            let carried = answer(&carrying.query("7", *lo, *hi).unwrap());
            let fresh = answer(
                &VariationTabixReader::new(&dir, &cols)
                    .unwrap()
                    .query("7", *lo, *hi)
                    .unwrap(),
            );
            assert_eq!(carried, fresh, "window {i}: {lo}..={hi}");
            assert_eq!(
                carried,
                by_scan(*lo, *hi),
                "window {i}: {lo}..={hi}, against a scan"
            );
        }
    }

    /// What makes the lookup fast is that a sorted input reads each line of the cache about once
    /// and not once for every variant: sweeping the dense region with a window for each of 3,000
    /// variants read the lines from the index's chunk start for every one of them, thousands of
    /// lines a time, before the window carried on.
    #[test]
    fn a_sorted_sweep_reads_each_line_about_once() {
        let (dir, cols) = window_cache();
        let reader = VariationTabixReader::new(&dir, &cols).unwrap();
        let mut p = 100_000u64;
        let mut asked = 0usize;
        while p < 120_000 {
            p += 3;
            reader.query("7", p - 1, p + 1).unwrap();
            asked += 1;
        }
        assert!(asked > 6_000);
        let open = reader.open.lock().unwrap();
        let source = open.by_dir.values().next().expect("an open source");
        // the dense region holds about 5,000 records below 120,000; every one of them, and the
        // line past the last window, and the ones the first seek landed before: a few thousand
        assert!(
            source.lines_read < 12_000,
            "{} lines read for {asked} windows",
            source.lines_read
        );
    }

    #[test]
    fn eviction_and_reopening_keep_every_answer() {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/vep_cache/csi");
        let info = crate::info::CacheInfo::from_file(&dir.join("info.txt")).unwrap();
        let reader = VariationTabixReader::new(&dir, &info.variation_cols).unwrap();
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

    /// `split(',', ...)` of Perl drops the empty fields at the end, so a trailing comma
    /// in the cache's `pubmed` is no PubMed id.
    #[test]
    fn a_trailing_comma_in_pubmed_is_no_id() {
        let cols: Vec<String> = "chr,variation_name,start,allele_string,pubmed"
            .split(',')
            .map(String::from)
            .collect();
        let reader = parsing_reader(&cols).unwrap();
        let r = reader.parse_line("1\trs1\t100\tA/G\t111,222,").unwrap();
        assert_eq!(r.pubmed, ["111", "222"]);
        let r = reader.parse_line("1\trs2\t100\tA/G\t111,,333").unwrap();
        assert_eq!(
            r.pubmed,
            ["111", "", "333"],
            "a field in the middle is kept"
        );
        let r = reader.parse_line("1\trs3\t100\tA/G\t.").unwrap();
        assert!(r.pubmed.is_empty());
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

    /// VEP's own `trim_sequences($ref, $alt, $start, $end, 1)`, run in the VEP
    /// 116.3 image, for pairs of the shape a lone indel alt takes: the left is
    /// trimmed first and moves the start, the right after it and moves the end.
    #[test]
    fn a_lone_indel_alt_is_minimised_as_trim_sequences_does() {
        // ref, alt, start, end => ref, alt, start, end
        let cases = [
            ("C", "CCC", 207651, 207651, "-", "CC", 207652, 207651),
            ("A", "-", 151, 151, "A", "-", 151, 151),
            ("CG", "G", 5851, 5852, "C", "-", 5851, 5851),
            ("-", "GG", 780, 779, "-", "GG", 780, 779),
            ("CA", "G", 100, 101, "CA", "G", 100, 101),
            ("AT", "ATT", 100, 101, "-", "T", 102, 101),
            ("ACCCT", "-", 10401, 10405, "ACCCT", "-", 10401, 10405),
            ("GA", "G", 10899, 10900, "A", "-", 10900, 10900),
            ("TCG", "TG", 100, 102, "C", "-", 101, 101),
            ("CAT", "GT", 202850, 202852, "CA", "G", 202850, 202851),
            ("AAAC", "AC", 50, 53, "AA", "-", 51, 52),
        ];
        for (r, a, s, e, want_r, want_a, want_s, want_e) in cases {
            let site = vep_input_site(r, &[a], s, e);
            assert_eq!(
                (site.reference.as_str(), site.alts, site.start, site.end),
                (want_r, vec![want_a.to_string()], want_s, want_e),
                "{r}>{a} at {s}-{e}"
            );
        }
    }

    /// An event in a repeat is not the same known variant once it is minimised.
    /// `AC>ACCC` reaches the matcher as `C/CCC` at 207651; written as it is, its
    /// right-first reading is the insertion of `CC` at 207651, which is what the
    /// cache has. VEP minimises the lone indel first (`-/CC` at 207652, the
    /// insertion after the second `C`) and calls the cache's record novel.
    #[test]
    fn a_minimised_insertion_in_a_repeat_does_not_match_its_neighbour() {
        let cache_record = "-/CC";
        // matched as written: a match, which VEP does not report
        assert!(!match_alleles("C", &["CCC"], 207651, cache_record, 207651, 1).is_empty());
        // as VEP matches it
        let site = vep_input_site("C", &["CCC"], 207651, 207651);
        let alts: Vec<&str> = site.alts.iter().map(String::as_str).collect();
        let forms = SiteForms::new(&site.reference, &alts, site.start, site.end);
        assert!(forms.match_alleles(cache_record, 207651, 1).is_empty());
        // the same insertion written in the cache where VEP puts it still matches
        assert!(!forms.match_alleles(cache_record, 207652, 1).is_empty());
    }

    /// Only a lone alt is minimised. SNVs and MNVs, and records with several
    /// alts (VEP 116 matches those whole), go through as they are.
    #[test]
    fn everything_but_a_lone_indel_alt_is_looked_up_as_it_is() {
        for (r, alts, s, e) in [
            ("C", vec!["A"], 10414, 10414),
            ("C", vec!["A", "G"], 10414, 10414),
            ("AAC", vec!["AAT"], 100, 102),
            ("CAT", vec!["GT", "TAT"], 202850, 202852),
            ("AT", vec!["G", "GCC"], 100, 101),
            ("-", vec!["A", "CC"], 100, 99),
        ] {
            let site = vep_input_site(r, &alts, s, e);
            assert_eq!(
                (
                    site.reference.as_str(),
                    site.alts.iter().map(String::as_str).collect::<Vec<_>>(),
                    site.start,
                    site.end
                ),
                (r, alts, s, e)
            );
        }
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

        let reader = VariationTabixReader::new(&cache_dir, &info.variation_cols).unwrap();

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
