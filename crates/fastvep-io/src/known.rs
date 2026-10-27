//! The columns VEP fills from the known variants of its variation cache:
//! `AF` and the population frequencies, `MAX_AF`, `CLIN_SIG`, `SOMATIC`,
//! `PHENO` and `PUBMED`.
//!
//! What goes in them is decided per input allele, from the known variants of the
//! whole site (`fastvep_annotate::colocated`), so it is worked out once for each
//! (variant x alt allele) and every transcript row of that allele only borrows
//! it. This module holds what the writers and that port share: the columns and
//! their order, the flags that ask for them, and the cells.
//!
//! Names and order are VEP's (ensembl-vep release/116, `Constants.pm`
//! `@FLAG_FIELDS`), checked against a run of VEP 116.3.

use fastvep_core::Allele;

/// One column VEP fills from the known variants.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KnownColumn {
    // 1000 Genomes
    Af,
    AfrAf,
    AmrAf,
    EasAf,
    EurAf,
    SasAf,
    // gnomAD exomes
    GnomadeAf,
    GnomadeAfrAf,
    GnomadeAmrAf,
    GnomadeAsjAf,
    GnomadeEasAf,
    GnomadeFinAf,
    GnomadeMidAf,
    GnomadeNfeAf,
    GnomadeRemainingAf,
    GnomadeSasAf,
    // gnomAD genomes
    GnomadgAf,
    GnomadgAfrAf,
    GnomadgAmiAf,
    GnomadgAmrAf,
    GnomadgAsjAf,
    GnomadgEasAf,
    GnomadgFinAf,
    GnomadgMidAf,
    GnomadgNfeAf,
    GnomadgRemainingAf,
    GnomadgSasAf,
    MaxAf,
    MaxAfPops,
    ClinSig,
    Somatic,
    Pheno,
    Pubmed,
}

impl KnownColumn {
    /// Every column, in the order VEP writes them when all its flags are given.
    pub const ALL: [KnownColumn; 33] = [
        Self::Af,
        Self::AfrAf,
        Self::AmrAf,
        Self::EasAf,
        Self::EurAf,
        Self::SasAf,
        Self::GnomadeAf,
        Self::GnomadeAfrAf,
        Self::GnomadeAmrAf,
        Self::GnomadeAsjAf,
        Self::GnomadeEasAf,
        Self::GnomadeFinAf,
        Self::GnomadeMidAf,
        Self::GnomadeNfeAf,
        Self::GnomadeRemainingAf,
        Self::GnomadeSasAf,
        Self::GnomadgAf,
        Self::GnomadgAfrAf,
        Self::GnomadgAmiAf,
        Self::GnomadgAmrAf,
        Self::GnomadgAsjAf,
        Self::GnomadgEasAf,
        Self::GnomadgFinAf,
        Self::GnomadgMidAf,
        Self::GnomadgNfeAf,
        Self::GnomadgRemainingAf,
        Self::GnomadgSasAf,
        Self::MaxAf,
        Self::MaxAfPops,
        Self::ClinSig,
        Self::Somatic,
        Self::Pheno,
        Self::Pubmed,
    ];

    /// The header, which is VEP's name for the column.
    pub fn name(self) -> &'static str {
        match self {
            Self::Af => "AF",
            Self::AfrAf => "AFR_AF",
            Self::AmrAf => "AMR_AF",
            Self::EasAf => "EAS_AF",
            Self::EurAf => "EUR_AF",
            Self::SasAf => "SAS_AF",
            Self::GnomadeAf => "gnomADe_AF",
            Self::GnomadeAfrAf => "gnomADe_AFR_AF",
            Self::GnomadeAmrAf => "gnomADe_AMR_AF",
            Self::GnomadeAsjAf => "gnomADe_ASJ_AF",
            Self::GnomadeEasAf => "gnomADe_EAS_AF",
            Self::GnomadeFinAf => "gnomADe_FIN_AF",
            Self::GnomadeMidAf => "gnomADe_MID_AF",
            Self::GnomadeNfeAf => "gnomADe_NFE_AF",
            Self::GnomadeRemainingAf => "gnomADe_REMAINING_AF",
            Self::GnomadeSasAf => "gnomADe_SAS_AF",
            Self::GnomadgAf => "gnomADg_AF",
            Self::GnomadgAfrAf => "gnomADg_AFR_AF",
            Self::GnomadgAmiAf => "gnomADg_AMI_AF",
            Self::GnomadgAmrAf => "gnomADg_AMR_AF",
            Self::GnomadgAsjAf => "gnomADg_ASJ_AF",
            Self::GnomadgEasAf => "gnomADg_EAS_AF",
            Self::GnomadgFinAf => "gnomADg_FIN_AF",
            Self::GnomadgMidAf => "gnomADg_MID_AF",
            Self::GnomadgNfeAf => "gnomADg_NFE_AF",
            Self::GnomadgRemainingAf => "gnomADg_REMAINING_AF",
            Self::GnomadgSasAf => "gnomADg_SAS_AF",
            Self::MaxAf => "MAX_AF",
            Self::MaxAfPops => "MAX_AF_POPS",
            Self::ClinSig => "CLIN_SIG",
            Self::Somatic => "SOMATIC",
            Self::Pheno => "PHENO",
            Self::Pubmed => "PUBMED",
        }
    }

    /// The column a header names, if it is one of these.
    pub fn from_name(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|c| c.name() == name)
    }
}

/// A population column of the cache, and the output column its value goes to.
#[derive(Debug, Clone, Copy)]
pub struct FrequencyKey {
    /// The name of the column in the cache's `variation_cols`.
    pub key: &'static str,
    /// The output column, if VEP names one. `ASN` is a key of `--af_1kg` that has
    /// no column (no 116 cache has the data either), but it still takes part in
    /// `MAX_AF` where a cache does.
    pub column: Option<KnownColumn>,
}

const fn key(key: &'static str, column: KnownColumn) -> FrequencyKey {
    FrequencyKey {
        key,
        column: Some(column),
    }
}

/// A set of populations that one VEP flag asks for.
#[derive(Debug, Clone, Copy)]
pub struct FrequencyGroup {
    /// The flag, in VEP's spelling.
    pub flag: &'static str,
    /// What the populations are, for a message.
    pub description: &'static str,
    /// The cache columns, in the order VEP reads them. That order decides which
    /// populations `MAX_AF_POPS` lists first.
    pub keys: &'static [FrequencyKey],
}

/// VEP's `%FREQUENCY_KEYS` (`OutputFactory.pm`), in the order it walks them:
/// `foreach my $group (sort @keys)`. `--af_gnomad` is the old name of
/// `--af_gnomade` and has the same populations, so it is one group here.
pub const FREQUENCY_GROUPS: [FrequencyGroup; 4] = [
    FrequencyGroup {
        flag: "af",
        description: "the 1000 Genomes `AF`",
        keys: &[key("AF", KnownColumn::Af)],
    },
    FrequencyGroup {
        flag: "af_1kg",
        description: "the 1000 Genomes",
        keys: &[
            key("AF", KnownColumn::Af),
            key("AFR", KnownColumn::AfrAf),
            key("AMR", KnownColumn::AmrAf),
            FrequencyKey {
                key: "ASN",
                column: None,
            },
            key("EAS", KnownColumn::EasAf),
            key("EUR", KnownColumn::EurAf),
            key("SAS", KnownColumn::SasAf),
        ],
    },
    FrequencyGroup {
        flag: "af_gnomade",
        description: "the gnomAD exomes",
        keys: &[
            key("gnomADe", KnownColumn::GnomadeAf),
            key("gnomADe_AFR", KnownColumn::GnomadeAfrAf),
            key("gnomADe_AMR", KnownColumn::GnomadeAmrAf),
            key("gnomADe_ASJ", KnownColumn::GnomadeAsjAf),
            key("gnomADe_EAS", KnownColumn::GnomadeEasAf),
            key("gnomADe_FIN", KnownColumn::GnomadeFinAf),
            key("gnomADe_MID", KnownColumn::GnomadeMidAf),
            key("gnomADe_NFE", KnownColumn::GnomadeNfeAf),
            key("gnomADe_REMAINING", KnownColumn::GnomadeRemainingAf),
            key("gnomADe_SAS", KnownColumn::GnomadeSasAf),
        ],
    },
    FrequencyGroup {
        flag: "af_gnomadg",
        description: "the gnomAD genomes",
        keys: &[
            key("gnomADg", KnownColumn::GnomadgAf),
            key("gnomADg_AFR", KnownColumn::GnomadgAfrAf),
            key("gnomADg_AMI", KnownColumn::GnomadgAmiAf),
            key("gnomADg_AMR", KnownColumn::GnomadgAmrAf),
            key("gnomADg_ASJ", KnownColumn::GnomadgAsjAf),
            key("gnomADg_EAS", KnownColumn::GnomadgEasAf),
            key("gnomADg_FIN", KnownColumn::GnomadgFinAf),
            key("gnomADg_MID", KnownColumn::GnomadgMidAf),
            key("gnomADg_NFE", KnownColumn::GnomadgNfeAf),
            key("gnomADg_REMAINING", KnownColumn::GnomadgRemainingAf),
            key("gnomADg_SAS", KnownColumn::GnomadgSasAf),
        ],
    },
];

/// The known-variant flags of VEP that fastVEP accepts.
///
/// `--af_esp` and `--af_exac` are not here because VEP 116 has no such options
/// (`Config.pm` lists `af`, `af_1kg`, `af_gnomade`, `af_gnomad`, `af_gnomadg`
/// and `max_af`), and a 116 cache has no ESP or ExAC columns.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct KnownFlags {
    /// `--af`: `AF`, the 1000 Genomes combined frequency.
    pub af: bool,
    /// `--af-1kg`: `AF`, `AFR_AF`, `AMR_AF`, `EAS_AF`, `EUR_AF`, `SAS_AF`.
    pub af_1kg: bool,
    /// `--af-gnomade` (VEP's `--af_gnomad` too): the gnomAD exomes columns.
    pub af_gnomade: bool,
    /// `--af-gnomadg`: the gnomAD genomes columns.
    pub af_gnomadg: bool,
    /// `--max-af`: `MAX_AF` and `MAX_AF_POPS`.
    pub max_af: bool,
    /// `--check-existing`: `CLIN_SIG`, `SOMATIC` and `PHENO`.
    pub check_existing: bool,
    /// `--pubmed`: `PUBMED`.
    pub pubmed: bool,
    /// `--clin-sig-allele 0`. VEP's default is 1, which prefers the significance
    /// a record gives for the allele over the one it gives for the variant, so
    /// the flag is the exception and the zero value of this struct is VEP's
    /// default.
    pub no_clin_sig_allele: bool,
}

impl KnownFlags {
    /// No known-variant flag, for a `const`, where `Default` cannot be called.
    pub const NONE: KnownFlags = KnownFlags {
        af: false,
        af_1kg: false,
        af_gnomade: false,
        af_gnomadg: false,
        max_af: false,
        check_existing: false,
        pubmed: false,
        no_clin_sig_allele: false,
    };

    /// What `--everything` turns on among these flags (`OPTION_SETS` in
    /// `Config.pm`): every frequency flag, `--max_af` and `--pubmed`.
    pub fn everything() -> Self {
        Self {
            af: true,
            af_1kg: true,
            af_gnomade: true,
            af_gnomadg: true,
            max_af: true,
            pubmed: true,
            ..Self::default()
        }
    }

    /// Every flag set in either.
    pub fn union(self, other: Self) -> Self {
        Self {
            af: self.af || other.af,
            af_1kg: self.af_1kg || other.af_1kg,
            af_gnomade: self.af_gnomade || other.af_gnomade,
            af_gnomadg: self.af_gnomadg || other.af_gnomadg,
            max_af: self.max_af || other.max_af,
            check_existing: self.check_existing || other.check_existing,
            pubmed: self.pubmed || other.pubmed,
            no_clin_sig_allele: self.no_clin_sig_allele || other.no_clin_sig_allele,
        }
    }

    /// Whether the run looks the known variants up for these columns. Every
    /// frequency flag, `--max_af` and `--pubmed` switch `--check_existing` on in
    /// VEP (`OPTION_SETS`), which is what adds `CLIN_SIG`, `SOMATIC` and `PHENO`.
    pub fn checks_existing(self) -> bool {
        self.check_existing
            || self.af
            || self.af_1kg
            || self.af_gnomade
            || self.af_gnomadg
            || self.max_af
            || self.pubmed
    }

    /// The flags that are set, by the names fastVEP gives them, for messages.
    pub fn requested(self) -> Vec<&'static str> {
        [
            (self.af, "--af"),
            (self.af_1kg, "--af-1kg"),
            (self.af_gnomade, "--af-gnomade"),
            (self.af_gnomadg, "--af-gnomadg"),
            (self.max_af, "--max-af"),
            (self.check_existing, "--check-existing"),
            (self.pubmed, "--pubmed"),
        ]
        .into_iter()
        .filter_map(|(set, name)| set.then_some(name))
        .collect()
    }

    /// Whether the frequency group of this flag is asked for.
    pub fn group(self, flag: &str) -> bool {
        match flag {
            "af" => self.af,
            "af_1kg" => self.af_1kg,
            "af_gnomade" => self.af_gnomade,
            "af_gnomadg" => self.af_gnomadg,
            _ => false,
        }
    }

    /// Whether `AF` is VEP's value for this run (`--af` or `--af-1kg`); without
    /// either, the `AF` of the CSQ is the one fastVEP has always written.
    pub fn writes_af(self) -> bool {
        self.af || self.af_1kg
    }

    /// The columns the flags add, in VEP's order (`@FLAG_FIELDS`: the frequency
    /// flags, `max_af`, `check_existing`, `pubmed`), each once: `AF` is in both
    /// `--af` and `--af_1kg`, and VEP keeps its first place.
    pub fn columns(self) -> Vec<KnownColumn> {
        let mut columns: Vec<KnownColumn> = Vec::new();
        let mut push = |c: KnownColumn| {
            if !columns.contains(&c) {
                columns.push(c);
            }
        };
        if self.af {
            push(KnownColumn::Af);
        }
        for group in &FREQUENCY_GROUPS[1..] {
            if self.group(group.flag) {
                for k in group.keys {
                    if let Some(c) = k.column {
                        push(c);
                    }
                }
            }
        }
        if self.max_af {
            push(KnownColumn::MaxAf);
            push(KnownColumn::MaxAfPops);
        }
        if self.checks_existing() {
            push(KnownColumn::ClinSig);
            push(KnownColumn::Somatic);
            push(KnownColumn::Pheno);
        }
        if self.pubmed {
            push(KnownColumn::Pubmed);
        }
        columns
    }
}

/// The cells of one input allele, as they read in VEP's tab file.
///
/// A list is joined with `,` (the VCF writer turns that into `&`). A column the
/// allele has nothing for has no entry, so a variant with no known variant, which
/// is nearly all of a typical callset, owns no strings at all.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct KnownCells {
    cells: Vec<(KnownColumn, String)>,
}

impl KnownCells {
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            cells: Vec::with_capacity(capacity),
        }
    }

    pub fn get(&self, column: KnownColumn) -> Option<&str> {
        self.cells
            .iter()
            .find(|(c, _)| *c == column)
            .map(|(_, v)| v.as_str())
    }

    /// Set the cell of a column, replacing a value it already has.
    pub fn set(&mut self, column: KnownColumn, value: String) {
        match self.cells.iter_mut().find(|(c, _)| *c == column) {
            Some(slot) => slot.1 = value,
            None => self.cells.push((column, value)),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.cells.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = (KnownColumn, &str)> {
        self.cells.iter().map(|(c, v)| (*c, v.as_str()))
    }
}

/// The cells of one alt allele of a variant.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlleleKnown {
    pub allele: Allele,
    pub cells: KnownCells,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(flags: KnownFlags) -> Vec<&'static str> {
        flags.columns().into_iter().map(KnownColumn::name).collect()
    }

    #[test]
    fn every_column_is_found_again_by_its_name() {
        for c in KnownColumn::ALL {
            assert_eq!(KnownColumn::from_name(c.name()), Some(c));
        }
        assert_eq!(KnownColumn::ALL.len(), 33);
        let mut seen: Vec<&str> = KnownColumn::ALL.iter().map(|c| c.name()).collect();
        seen.sort_unstable();
        seen.dedup();
        assert_eq!(seen.len(), 33, "two columns share a name");
    }

    /// `@FLAG_FIELDS`, `Constants.pm` of VEP 116: `af` has `AF`; `af_1kg` has `AF` first
    /// and then the five populations; `AF` keeps the place `af` gives it.
    #[test]
    fn columns_come_in_the_order_of_vep_flag_fields() {
        assert_eq!(
            names(KnownFlags {
                af_1kg: true,
                ..Default::default()
            }),
            [
                "AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF", "CLIN_SIG", "SOMATIC",
                "PHENO"
            ]
        );
        assert_eq!(
            names(KnownFlags {
                af: true,
                af_1kg: true,
                ..Default::default()
            }),
            names(KnownFlags {
                af_1kg: true,
                ..Default::default()
            }),
            "`--af --af_1kg` is `--af_1kg`: AF is written once"
        );
        let all = names(KnownFlags::everything());
        assert_eq!(all.len(), 33);
        assert_eq!(
            all.iter().position(|n| *n == "MAX_AF"),
            Some(27),
            "the frequency columns, then MAX_AF"
        );
        assert_eq!(
            &all[28..],
            ["MAX_AF_POPS", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED"]
        );
        assert_eq!(
            names(KnownFlags {
                pubmed: true,
                ..Default::default()
            }),
            ["CLIN_SIG", "SOMATIC", "PHENO", "PUBMED"],
            "`--pubmed` switches `--check_existing` on, which adds three columns before it"
        );
        assert!(names(KnownFlags::default()).is_empty());
        assert!(
            names(KnownFlags {
                no_clin_sig_allele: true,
                ..Default::default()
            })
            .is_empty(),
            "`--clin-sig-allele 0` adds no column and looks nothing up"
        );
    }

    #[test]
    fn a_gnomad_flag_adds_ten_columns_and_a_genome_flag_eleven() {
        let exomes = names(KnownFlags {
            af_gnomade: true,
            ..Default::default()
        });
        assert_eq!(exomes[0], "gnomADe_AF");
        assert_eq!(exomes[9], "gnomADe_SAS_AF");
        assert_eq!(exomes.len(), 10 + 3);
        let genomes = names(KnownFlags {
            af_gnomadg: true,
            ..Default::default()
        });
        assert_eq!(genomes[2], "gnomADg_AMI_AF");
        assert_eq!(genomes.len(), 11 + 3);
    }

    #[test]
    fn cells_are_set_once_and_read_back() {
        let mut cells = KnownCells::default();
        assert!(cells.is_empty());
        cells.set(KnownColumn::Af, "0.1000".into());
        cells.set(KnownColumn::Pubmed, "1,2".into());
        cells.set(KnownColumn::Af, "0.2000".into());
        assert_eq!(cells.get(KnownColumn::Af), Some("0.2000"));
        assert_eq!(cells.get(KnownColumn::Pubmed), Some("1,2"));
        assert_eq!(cells.get(KnownColumn::MaxAf), None);
        assert_eq!(cells.iter().count(), 2);
    }
}
