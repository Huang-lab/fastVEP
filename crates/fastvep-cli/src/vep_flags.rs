//! What to tell a user who types a flag of Ensembl VEP that fastVEP does not have.
//!
//! `fastvep annotate` refuses every flag it does not know, with clap's `unexpected argument`, and
//! clap adds a guess at what was meant when one of its own flags is spelled nearly the same. For
//! VEP's flags the guess is wrong more often than not: `--assembly` was "similar" to `--symbol`,
//! `--json` to `--sa-only`, `--bam` to `--acmg`, `--no_stats` to `--no_escape`, `--protein` to
//! `--proband` and `--uniprot` to `--input`, and `--use_transcript_ref` to `--transcript-cache`,
//! which is another thing. A command line copied from a `vep` run is exactly where a user meets them,
//! so each of the flags `docs/VEP_MIGRATION.md` lists as refused says what there is instead, and
//! nothing else is said about a flag this does not know.

/// Flags of VEP and what a user of one wants to be told. `{flag}` is the flag as it was typed.
struct Refused {
    /// The names, as VEP spells them, with no dashes.
    flags: &'static [&'static str],
    tip: &'static str,
}

const SEE: &str = "docs/VEP_MIGRATION.md lists the flags fastVEP has.";

const REFUSED: &[Refused] = &[
    Refused {
        flags: &["cache", "offline", "dir_cache"],
        tip: "fastVEP never opens VEP's transcript cache or a database: the gene models come from a \
              GFF3 (`--gff3`) and the reference from a FASTA (`--fasta`), and the known variants of \
              a VEP variation cache are read from `--cache-dir DIR`.",
    },
    Refused {
        flags: &["refseq", "merged"],
        tip: "the transcripts are the GFF3 you pass: `--gff3 refseq.gff3`, or both of them with \
              `--gff3 Ensembl=a.gff3 --gff3 RefSeq=b.gff3`, which the SOURCE field tells apart.",
    },
    Refused {
        flags: &["tab", "vcf", "json"],
        tip: "the format is `--output-format tab`, `vcf` (the default) or `json`.",
    },
    Refused {
        flags: &["assembly", "species"],
        tip: "there is no such flag: the assembly and the species are whatever the GFF3 and the \
              FASTA are.",
    },
    Refused {
        flags: &["use_transcript_ref", "bam"],
        tip: "fastVEP has no BAM-edited transcripts: it builds each transcript from the GFF3 and the \
              FASTA.",
    },
    Refused {
        flags: &["fields"],
        tip: "the columns of a tab file are the ones the flags ask for, in VEP's order, and there is \
              no way to choose them one by one.",
    },
    Refused {
        flags: &["plugin", "custom"],
        tip: "there is no plugin API and no custom track: `fastvep sa-build` builds a supplementary \
              annotation database from a source such as gnomAD, REVEL, SpliceAI or any VCF or BED, \
              and `--sa-dir` reads it.",
    },
    Refused {
        flags: &["fork"],
        tip: "fastVEP uses all the cores of the machine on its own; set RAYON_NUM_THREADS=N to \
              limit it.",
    },
    Refused {
        flags: &["no_stats"],
        tip: "fastVEP writes no statistics file, so there is nothing to switch off.",
    },
    Refused {
        flags: &["quiet"],
        tip: "there is no {flag}: `--no-progress` turns the progress display off, and the warnings \
              go to stderr.",
    },
    Refused {
        flags: &["mane", "tsl", "appris", "ccds", "protein"],
        tip: "there is no {flag}: the VCF `CSQ` always has the MANE, TSL, APPRIS, CCDS and ENSP \
              fields, filled where the GFF3 says.",
    },
    Refused {
        flags: &[
            "sift",
            "polyphen",
            "domains",
            "uniprot",
            "gene_phenotype",
            "regulatory",
            "nearest",
            "shift_3prime",
            "no_intergenic",
            "coding_only",
            "minimal",
            "allele_number",
            "individual",
            "vcf_info_field",
            "terms",
        ],
        tip: "fastVEP has no counterpart of VEP's {flag}.",
    },
    Refused {
        flags: &["af_esp", "af_exac"],
        tip: "VEP 116 has no {flag} either; its frequency flags are `--af`, `--af_1kg`, \
              `--af_gnomade`, `--af_gnomadg` and `--max_af`, and fastVEP has those.",
    },
];

/// The tip for `flag`, as the user typed it (`--no_stats`, `--no-stats`, `--no_stats=1`), when it is
/// a flag of VEP that fastVEP does not have.
pub fn tip(flag: &str) -> Option<String> {
    let name = flag
        .trim_start_matches('-')
        .split('=')
        .next()
        .unwrap_or_default()
        .replace('-', "_")
        .to_ascii_lowercase();
    let refused = REFUSED.iter().find(|r| r.flags.contains(&name.as_str()))?;
    let flag = flag.split('=').next().unwrap_or(flag);
    Some(format!("{} {SEE}", refused.tip.replace("{flag}", flag)))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The flags `docs/VEP_MIGRATION.md` lists as refused, each of which was run against the binary
    /// when the page was written.
    const LISTED: &[&str] = &[
        "cache",
        "offline",
        "dir_cache",
        "refseq",
        "merged",
        "tab",
        "vcf",
        "json",
        "assembly",
        "species",
        "use_transcript_ref",
        "bam",
        "fields",
        "plugin",
        "custom",
        "fork",
        "no_stats",
        "quiet",
        "sift",
        "polyphen",
        "mane",
        "tsl",
        "appris",
        "ccds",
        "protein",
        "domains",
        "uniprot",
        "gene_phenotype",
        "regulatory",
        "nearest",
        "shift_3prime",
        "no_intergenic",
        "coding_only",
        "minimal",
        "allele_number",
        "individual",
        "vcf_info_field",
        "terms",
        "af_esp",
        "af_exac",
    ];

    #[test]
    fn every_flag_the_migration_page_lists_has_a_tip_in_either_spelling() {
        for name in LISTED {
            for typed in [
                format!("--{name}"),
                format!("--{}", name.replace('_', "-")),
                format!("--{name}=1"),
            ] {
                let tip = tip(&typed).unwrap_or_else(|| panic!("no tip for {typed}"));
                assert!(tip.ends_with(SEE), "{typed}: {tip}");
                assert!(!tip.contains("{flag}"), "{typed}: {tip}");
            }
        }
    }

    #[test]
    fn a_flag_that_is_not_one_of_them_has_none() {
        // a typo of a flag fastVEP has keeps clap's own suggestion, which is right for it
        for typed in [
            "--hgv",
            "--symbols",
            "--nope",
            "--cache-dir",
            "--output-format",
            "-x",
        ] {
            assert_eq!(tip(typed), None, "{typed}");
        }
    }

    #[test]
    fn the_tip_names_what_to_use_and_the_flag_where_it_says_so() {
        assert!(tip("--json").unwrap().contains("--output-format"));
        assert!(tip("--refseq").unwrap().contains("--gff3"));
        assert!(tip("--dir_cache").unwrap().contains("--cache-dir"));
        assert!(tip("--quiet").unwrap().contains("there is no --quiet"));
        assert!(tip("--no-intergenic")
            .unwrap()
            .contains("VEP's --no-intergenic"));
        // the value of `--flag=value` is not part of the flag
        assert!(tip("--sift=b").unwrap().contains("VEP's --sift."));
    }

    #[test]
    fn no_flag_is_in_two_groups() {
        let mut seen = std::collections::BTreeSet::new();
        for group in REFUSED {
            for flag in group.flags {
                assert!(seen.insert(*flag), "{flag} is in two groups");
            }
        }
    }
}
