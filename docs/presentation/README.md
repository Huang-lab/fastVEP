# The fastVEP deck

`fastvep-deck.html` is a presentation of fastVEP and of its ACMG/AMP classification, for people who want to know why to use them and how.
It has 29 slides and 4 slides of appendix.
It is one self-contained HTML file: its styles, scripts and icons are inline, its type is the system's, and it makes no network request of any kind.
It works from a file, offline, on a laptop screen, on a projector and on a phone.
It was tested in Chrome 156 on macOS, and it uses only standard features of current browsers, but it was not tested in Safari or Firefox.

## What it describes, and as of when

- fastVEP is the pre-release v0.5.0-rc.1, and the latest stable release is v0.4.0.
- The speed and agreement-with-VEP numbers were measured on 7 November 2026 on an Apple M4, with a binary built from the tag of that release (v0.5.0-rc.1, c12d9f9; `benchmarks/results/v0.5.0/with_vep/environment.tsv:13-17`), and are the files of `benchmarks/results/v0.5.0/`.
- Every command on a slide was run on 8 November 2026 with the release binary of `fastVEP-v0.5.0-rc.1-aarch64-apple-darwin.tar.gz`, or is shown in the shape that its documentation gives, with its flags checked against the help of that binary (the tables below say which).
- Every ACMG number is v38's, except that the version history on slide 30 plots the values of v1 to v38 as the manuscript records them: each version at its run of the time, and v37 and v38 on the final benchmark databases.
  They are from the manuscript at commit `2a293d0` of the branch `docs/acmg-scheme-and-versions` (`manuscript/acmg/manuscript.md` and the tables of its `data/` directory).
  v38 is the classifier of release 0.5.0, which for now is the release candidate.
- The manuscript is unpublished.
  Its branch `docs/acmg-scheme-and-versions` is not merged to master, and the commit that the footers cite, `2a293d0`, is a commit of that branch, which is opened on that branch once the author has pushed it (`git ls-remote origin` listed the branch at `7814589` on 8 November 2026, and the commits after it, up to `2a293d0`, are pushed by the author with the paper).
  The deck cites no publication for it.
- The deck says what the manuscript says and no more: agreement with ClinVar 2-star+ and with one reviewer is how closely the calls agree with community opinion, and it is not accuracy and not validation.
  The four bounds of slide 28 are the manuscript's, in its words.
- Slides 7, 13, 14, 15, 27 and 29 say something about https://fastvep.org as it was on 8 November 2026.
  The server runs a web build (`0.5.0-rc.1+web.e068377`, read from `/api/status` at 18:58 UTC and again at 19:15 EST that day) that is newer than the tag, so what those slides say about the page, about all nine genomes at once and about ACMG is a statement about the site and not about the tag.
  The `fastvep-web` of the tag serves one genome, has no `acmg` object in `/api/status`, and answers `"acmg": true` with a class made from whatever is loaded.
  The site offers ACMG for GRCh38, and refuses it with HTTP 422 for a genome without the stack (a request for `mouse_grcm39_ensembl_115` at about 19:27 EST); `/api/genomes` marks the other eight genomes as not available.
- The cost of a request on slide 27 is the measured figures of `DEPLOYMENT.md` of the branch `feat/web-genomes` at commit `95e4d2a` (section "What has been measured on the 4 vCPU, 16 GB host", 8 November 2026).
  That branch is on GitHub and is not merged to master.
- Which numbers come from which commit, and what to check again if the manuscript or the site changes, is under "Which numbers come from which commit".

## Present it

Open `fastvep-deck.html` in a browser.

| Do | Keys |
|---|---|
| Next slide | `→`, `Space`, `PageDown` |
| Previous slide | `←`, `Shift+Space`, `PageUp` |
| First and last slide | `Home`, `End` |
| Overview of every slide (arrows move, `Enter` opens) | `O` |
| Speaker notes under the slide | `N` |
| Full screen | `F` |
| Light and dark | `D` |
| List of keys | `?` |

A click on the left or the right third of a slide moves, and so does a swipe.
The address holds the slide number, so `fastvep-deck.html#7` opens slide 7, and `?notes` in the address opens the notes first.
Every slide has speaker notes of a few sentences, with the source of each number.
The page follows the light or dark setting of the system, and `D` overrides it.
It respects the setting that asks for less motion.
On a phone the slide reflows and scrolls inside itself, and the bar at the bottom holds the controls.
Every command block has a Copy button, which uses the clipboard where the page may and a selection where it may not, so it works from a file and in a frame.

## Print to PDF

Use the print dialog of the browser.
The page sets its own page size, 13.333 by 7.5 inches, and no margins, so each slide is one page.
If your browser leaves out backgrounds, turn background graphics on: in Chrome the pages come out the same with it on and off.
Press `N` before printing to put the speaker notes under each slide: the pages are then 13.333 by 10.5 inches.
The printed slides are always on a light background.
The print was checked on the PDF files that Chrome's print engine writes, and not in its print dialog (see What was checked).

## What is on each slide

| # | Title |
|---|---|
| 1 | fastVEP: fast variant annotation and ACMG/AMP classification |
| 2 | Two problems we set out to fix |
| 3 | Ensembl VEP, reimplemented in Rust |
| 4 | 134 to 154 times faster than Ensembl VEP |
| 5 | What it costs: memory and a start-up time |
| 6 | Compared with VEP cell by cell: 0.077% differ |
| 7 | Nine genomes, the databases you use, and a filter |
| 8 | Install: one archive, two binaries |
| 9 | Your first annotation |
| 10 | The flags that matter |
| 11 | Three output formats: VCF, tab and JSON |
| 12 | Build a database once, annotate with it often |
| 13 | No install: fastvep.org in three steps |
| 14 | A REST API for notebooks and internal services |
| 15 | Guards against a plausible wrong answer |
| 16 | The ACMG problem and our contribution |
| 17 | How it was built: implement, review, measure |
| 18 | Nine sources, 15 of 28 criteria |
| 19 | Points on one side, the 2015 table on the other |
| 20 | How closely it agrees with community opinion |
| 21 | What the disagreement is made of |
| 22 | Agreement needs contradictions beside it |
| 23 | What we learned |
| 24 | Build the nine-source stack |
| 25 | Check that every source answers |
| 26 | Run it, then read the evidence |
| 27 | ACMG on the web and the API |
| 28 | What it is not |
| 29 | Where to go next |
| 30 | How agreement moved over 38 versions (appendix) |
| 31 | Coming from VEP: what carries over (appendix) |
| 32 | Tuning ACMG: thresholds, trios, functional evidence (appendix) |
| 33 | Present this deck (appendix) |

## Where every number comes from

Paths are of this repository at the commit that the deck branches from (`9574be1`), unless a row says otherwise.
`summary.tsv` is `benchmarks/results/v0.5.0/with_vep/summary.tsv`, with the line numbers of the file.
`MS` is `manuscript/acmg/manuscript.md` at commit `2a293d0`, with its line numbers, and `MS data` is `manuscript/acmg/data/` at that commit.
`DEPLOYMENT` is `DEPLOYMENT.md` of the branch `feat/web-genomes` at commit `95e4d2a`, with its line numbers.
A date is the day on which the command was run or the request made.

| Slide | What | Source |
|---|---|---|
| 1 | v0.5.0-rc.1 is a pre-release, v0.4.0 is the latest | `gh release view v0.5.0-rc.1` and `gh release list --repo Huang-lab/fastVEP`, 8 Nov 2026 |
| 1 | Apache-2.0 | `README.md:912`, `gh repo view Huang-lab/fastVEP` |
| 1, 4 | 4,048,342 variants of GIAB HG002; 42.2, 44.5 and 54.6 s on one thread; 11.1 and 17.0 s on every core | `summary.tsv:18,20,24` and `:19,25`, `with_vep/environment.tsv:26`, `README.md:797-803` |
| 1, 2, 4 | The call set has chromosomes 1 to 22 only | `benchmarks/run_vep_side.py:17,144`, `with_vep/vep_genome_sums.csv:2-4` |
| 2, 4 | VEP 116.3 took 6,516, 6,662 and 7,334 s (109 minutes for the first) | `with_vep/vep_baseline.csv:6-8`, `README.md:799-801` |
| 2 | The criteria are underspecified and fragmented | MS:20, 68-70 |
| 3 | 49 consequence terms are defined and 38 are emitted | `README.md:13`, `benchmarks/results/v0.5.0/so_terms.csv` |
| 3, 7, 11 | The CSQ has 49 fields (50 with `--acmg`), matching Ensembl VEP's extended format, of which the last two are fastVEP's; formats vcf, tab and json | `README.md:20,655`, `fastvep filter` printed `CSQ format: 49 fields` and `50 fields` on 8 Nov 2026; no file of this deck compares the field order with VEP's own VCF output |
| 3, 10 | The engine is a port of VEP's, it writes VEP's columns, and it accepts the VEP flags it has in VEP's spelling | `docs/VEP_MIGRATION.md:4-6,40`, `fastvep annotate --help` |
| 3, 8 | Two binaries, `fastvep` and `fastvep-web`; the subcommands | the release archive and `fastvep --help`, 8 Nov 2026 |
| 4 | Apple M4, 10 logical cores (4 performance, 6 efficiency), 24 GiB, macOS 26.5.1, medians of 5, `--hgvs --symbol --canonical` | `with_vep/environment.tsv:3-7,10,19,22`, `README.md:767` |
| 4 | The benchmark binary is a build of the tag (v0.5.0-rc.1, c12d9f9), and not the file of the archive, which was not timed | `with_vep/environment.tsv:13-17` |
| 4 | VEP in Docker, `--gff`, per chromosome and summed, the same ClinVar and gnomAD files | `README.md:770,791,792,805` |
| 4 | fastVEP read its binary transcript cache, which was built once in 13 s (12.95 s with the default threads, 1.23 cores busy on average) and is not in the times; 103,525 of the variants have a gnomAD record | `README.md:791`, `benchmarks/results/v0.5.0/summary.tsv:2` (first run, threads all), `with_vep/environment.tsv:33-34` (the sidecar's checksum is the same before and after), `with_vep/vep_genome_sums.csv:4` |
| 4 | The speedups are VEP's time over fastVEP's time: 154, 150 and 134 | 6,516 / 42.19, 6,662 / 44.46, 7,334 / 54.63; `README.md:799-801` |
| 4 | The bars are drawn to one scale, none thinner than 6 px | the CSS of the page (`min-width` of a bar): to scale the fastVEP bars would be 2.2 to 2.9 px wide |
| 5 | Peak memory for 100,000 variants: 959 MiB for VEP and 2,700 MiB for fastVEP (2,702 MiB on eight threads) | `with_vep/vep_summary.tsv:8`, `summary.tsv:10-11`, `README.md:819,826` |
| 5 | Whole-genome peaks of 4,823 MiB and 6,224 MiB on one thread | `summary.tsv:18,24` |
| 5 | Start-up of 0.86 s on one thread and 0.37 s on eight, 509,650 transcripts | `summary.tsv:2-3`, `README.md:807` |
| 5 | The chr22 table; its speedups are VEP's median over fastVEP's unrounded median (0.877, 0.935, 1.048 and 1.544 s), which gives 8.5, 23.0 and 72.9 where a division of the rounded cells of `summary.tsv` and `README.md:813-815` give 8.4, 23 and 73 | `with_vep/raw.tsv` (rounds 1 to 5 of the four `chr22_*_t1` cells), `summary.tsv:6-9`, `with_vep/vep_summary.tsv:2-5` for VEP, `README.md:812-815` |
| 5 | The install: archive of 4,932,188 bytes, `fastvep` of 6,000,304 bytes and `fastvep-web` of 4,394,080 bytes on macOS, 7,472,848 and 5,343,440 bytes on Linux x86_64; VEP's image of 0.97 GB unpacked and 0.31 GB to pull | `ls -l` of the archives and `gh release view`, 8 Nov 2026; `with_vep/vep_resources.csv:5-6` |
| 5 | The binaries need only system libraries: libSystem and libiconv on macOS, libc, libm and libgcc_s and glibc 2.34 or newer on Linux | `otool -L` of the macOS binaries and the symbol versions of the Linux binary, 8 Nov 2026; `README.md:820-821` |
| 6 | 5,023,485, 72,134,785, 84,679,965, 216,390, 80,379 and 73,359 cells; 3,851, 112, 2,999, 954, 0 and 818 differ; VEP 116.3 was run with `--gff` on the same GFF3 | `benchmarks/results/v0.5.0/vep116_agreement.csv` (the RefSeq row is line 10), `README.md:853-859`, `docs/VEP_MIGRATION.md:166` |
| 6 | 27 of the 30 columns are compared: SYMBOL_SOURCE, HGNC_ID and FLAGS cannot be filled from a GFF3 | `validation/measure_vep116.py:24-25`, `docs/VEP_DIVERGENCE.md:690`, `README.md:847` |
| 6 | The 3,851 are in four columns (2,070, 1,729, 35 and 17); the other 23 compared columns agree on every row; 3,822 are multi-ALT | `vep116_hg002_columns.csv`, `README.md:863-865`, `docs/VEP_MIGRATION.md:166` |
| 6 | Against VEP's own cache (`vep --cache`), FLAGS and HGNC_ID, which a GFF3 cannot fill, differ on 35 and 533 of 81,309 HG002 rows of chr21 and chr22 with the Ensembl cache (HGVSc and HGVS_OFFSET differ on 338 and 437) | `docs/VEP_MIGRATION.md:167-168` |
| 6 | 128 rows are in one output only (121 fastVEP's, 7 VEP's), and not all of them were checked | `vep116_agreement.csv:4`, `docs/VEP_DIVERGENCE.md:712-715` |
| 6, 31 | Not a drop-in: no plugins, no `--custom`, no database mode | `docs/VEP_MIGRATION.md:5` |
| 6 | `docs/VEP_DIVERGENCE.md` gives each divergence found and says which tool the project believes is right | `README.md:894-895`, `docs/VEP_MIGRATION.md:169` |
| 7 | The nine genomes | `https://fastvep.org/api/genomes`, 8 Nov 2026 |
| 7 | All nine at once is the web build; the release's `fastvep-web` serves one genome and switches only with `--allow-model-replacement` | the live page and `/api/genomes`, 8 Nov 2026; `README.md:357`, `docs/API.md:264-289`; the release's `fastvep-web --data-dir` run on two genomes (HTTP 200 for the list, HTTP 403 for `POST /api/load-genome`), 8 Nov 2026 |
| 7, 12 | The supplementary sources and their build commands | `README.md:432-458`, `docs/SA_DATABASES.md:45-64` |
| 7 | COSMIC, PrimateAI and TOPMed are marked as account sources | `docs/SA_DATABASES.md:59-61,66` |
| 7 | The filter expression | `README.md:611` (Command Reference); run 8 Nov 2026 |
| 8 | The download recipe, the three platforms, glibc 2.34 | `README.md:62-85` |
| 8 | Six assets of the release | `gh release view v0.5.0-rc.1`, 8 Nov 2026 |
| 9 | The command and its 8 rows, 280 rows without `--pick` (each of the seven BRCA1 variants overlaps 35 to 42 of the 48 transcripts, the TP53 variant one), the warning; one row per allele with `--pick` (all 8 variants have one ALT allele) | `README.md:109-116`, run 8 Nov 2026 on `tests/test.vcf` and `tests/test.gff3`; `fastvep annotate --help` for `--pick`, and a run on a record with two ALT alleles, which gave two rows |
| 9 | `rs_3utr_brca1` lies in the last coding block of the reverse-strand transcript, so `coding_sequence_variant` is right | `tests/test.gff3`: CDS 43045678-43045802 of ENST00000357654 and its 3' UTR 43044295-43045677 |
| 9, 11 | The tab output has 20 columns (4 and 5 are shown) | the header line of the tab output, run 8 Nov 2026 |
| 9 | `tests/` is in the repository and not in the release archive | the file list of the archive (`tar -tzf`), 8 Nov 2026 |
| 10 | The flags and what they do | `fastvep annotate --help` of v0.5.0-rc.1, 8 Nov 2026; `docs/VEP_MIGRATION.md:40,62` |
| 11 | The three outputs of one run on two ClinVar variants; the JSON keys left out are marked with an ellipsis | run 8 Nov 2026; `README.md:20`; `docs/VEP_MIGRATION.md:154-156` |
| 12 | `sa-build` of a custom panel and the two `FV_MYPANEL` values | `README.md:464-512`, run 8 Nov 2026 |
| 12 | No pre-built package, the small tier of about 1 GB, every database of `--sa-dir` is loaded; the script builds the public sources it lists | `docs/SA_DATABASES.md:7,9,18,33,47-57`, `scripts/build-sa-databases.sh:29-30` |
| 12 | ClinVar's 4,614,856 records make a file of 43,242,353 bytes (`.osa2`) | the run of the small tier of 8 Nov 2026 (`clinvar: 4614856 records` and `ls -l` of `clinvar.osa2`) |
| 13 | The three steps, the options (ACMG-AMP classification ticked and enabled for GRCh38), the tabs and downloads, the limits of 10 MiB and 90 seconds, the loaded databases (`sa_sources`: ClinVar, dbSNP, gnomAD, PHYLOP, repeatmasker, REVEL, SpliceAI) | https://fastvep.org (the page loaded in a browser at 19:44 UTC, and its script), `/api/status` at 18:58 UTC, 8 Nov 2026 |
| 13, 14 | The hosted instance is for interactive use | `docs/API.md:6-16`, `README.md:201` |
| 14 | The server command and the request | `docs/API.md:48-54,159-162`; run on a local server, and once on fastvep.org, 8 Nov 2026 |
| 15 | The failure that matters is a plausible-looking wrong answer | `CLAUDE.md:6` |
| 15 | Nine older cache formats are refused: the sidecar is rebuilt and a named cache stops the run | `crates/fastvep-cache/src/transcript_cache.rs:32-58,810-874,1283-1290` (FSTVEP01 to FSTVEP09; FSTVEP10 is current), `crates/fastvep-cli/src/pipeline/annotate.rs:1094-1135` |
| 15 | 11 runs write the same bytes | `README.md:849`, `benchmarks/results/v0.5.0/identity_vs_5048fc5/` |
| 15 | A run without the stack classifies from what is loaded; the stack check exits non-zero; fastvep.org refuses a genome without the stack with HTTP 422 | `docs/ACMG.md:19`, `docs/ACMG_SETUP.md:394`; `https://fastvep.org/api/annotate` answered HTTP 422 to `"acmg": true` for `mouse_grcm39_ensembl_115` at about 19:27 EST on 8 Nov 2026; `crates/fastvep-web/src/acmg.rs` and `errors.rs` of the branch `feat/web-genomes` (not in the tag or in master) |
| 15 | 1,731 tests passed, 0 failed, 1 ignored (macOS) | `cargo test --workspace` on this commit, 8 Nov 2026; `README.md:759` gives 1,732 on Linux |
| 15 | 20,363 records in nine fixtures of VEP 116.3's own output | `README.md:870`, `crates/fastvep-cli/fixtures/vep_tab_*/in.vcf` |
| 16 | The contribution, community opinion, the four measures | MS:20-21, 68-75 |
| 17 | 38 versions, 6 rounds, 673,660 variants | MS:28, 94, 199, 249 |
| 17 | The reviewer read the opposite-direction calls only; most thresholds were swept, three first for the paper, and seven numbers were not varied; the unit of the process is a run, re-run after every change | MS:246 (what the reviewer saw), MS:262 and 270 (the sweeps, the three and the seven numbers), MS:222-223 (the run) |
| 17 | The opposite-direction calls are counted from v6 on | MS data `trajectory.csv` (v1 holds a count reconstructed from its documentation, v2 to v5 none, v6 the first one recorded), MS:235, 393 |
| 18 | The nine sources and what they feed | `scripts/check_acmg_stack.py` output of 8 Nov 2026, MS:173-190 (Table 2), `docs/ACMG_SETUP.md:26-39` |
| 18 | 15 of 28 criteria can be applied to a sites-only input, and why 13 are not | MS:24, 42, 124-137, MS data `criteria_status.csv` |
| 18 | What ClinVar feeds: the PM2 backstop, a gate for BS1 and BS2, the counts of PS1 and PM1 | MS:186 (Table 2); the check lists PM3, BP2, PP5 and BP6 as well (`docs/ACMG_SETUP.md:32`) |
| 19 | The point scale, the 2015 table, and the reason for the split | MS:142-147 |
| 19 | 66.4% to 77.0% and 27 to 399 (a replay, Likely benign from -1 point); 62.5% to 57.1% and 10 fewer (a run) | MS:600-607, MS data `combining_rules_replay_v38.csv` |
| 19 | Recall is the share of a ClinVar group that receives a call in the same direction | MS:208 |
| 20 | 64.1%, 79.3%, 51 and their intervals | MS:668, 693, MS data `headline_v38.csv` |
| 20 | 116,078 new variants, 72.6%, 2.3 times (1.2 to 3.8) | MS:856, 862, 879 (Table 7), MS data `temporal_v38.csv` and `temporal_ratio_v38.csv` |
| 20, 21, 22 | Agreement in direction means the same side (pathogenic, benign or uncertain), so a VUS against a ClinVar VUS counts: 286,815 of the 431,727 exact records; every variant VUS scores 43.8% | MS:213, 702, MS data `confusion_3x3_v38.csv` (VUS to VUS 286,815) |
| 20 | 62 under a strict leave-one-out reading; the new variants lie in 4,862 genes, all but 157 in the April set (a hold-out of variants and almost not of genes); 98.4% had a record in the April variant_summary and 92% of all of them a single-submitter one (114,218 with a record, of which 106,607 single-submitter, 6,505 multiple submitters without conflict, 54 expert panel, 1,052 other) | MS:1198, 1243; MS:839, 903, MS data `temporal_new_genes_v38.csv` and `temporal_new_in_index_v38.csv` |
| 21 | The segments of the bar, 93.6% of 139,325; the share 1.3% of the segment of ClinVar VUS called pathogenic- or benign-side is computed (3,610 and 4,873 of 673,660, which Table 5 gives as 0.5% and 0.7%) | MS:671-692 (Table 5), MS data `confusion_3x3_v38.csv` and `headline_v38.csv` |
| 21 | 10,694, 5,184, 20,072 (17,098 with a lone BP4 and 2,974 with a lone BS1) and the 104 missing curations; the paper does not apportion them | MS:726, 732, 768, 751, 1156-1159 |
| 22 | The four rows and 36 times | MS:702-706, MS data `baselines_v38.csv` |
| 22 | None of the 4,039 Pathogenic or 67,733 Benign calls is among the 51; the zero for Benign is informative (15 expected) and the zero for Pathogenic is not (1.7 expected); 62 and three added Benign calls under the strict reading | MS:710-711, 1243, 633 |
| 23 | Split the combining scheme: points on both sides (a replay, Likely benign from -1 point) gain 10.6 points of benign-side recall for 372 more contradictions; the 2015 table on both sides gives up 5.4 points of pathogenic-side recall to avoid 10 | MS:147, 601-606, 1087-1089, MS data `combining_rules_replay_v38.csv` |
| 23 | 21 of the 64 changes were a condition of a guideline that a criterion did not apply, and review of the opposite-direction calls first identified 30 of the 64 (by our coding) | MS:397-402, 1090, MS data `change_ledger.csv` |
| 23 | The Likely pathogenic cutoff, BP7's conservation test, BS1 and BS2 matter; the PM2 thresholds and the allele-number floor move agreement in direction by under 1,000 calls (534,835 against 534,335 for PM2; 684 for the floor) | MS:1091-1092, 491, 614, and Table 3 (MS:463-483) for the exchange rates |
| 23 | The three rules reach 85.6% with 1,817 contradictions, v38 79.3% with 51 | MS:704, 1095, MS data `baselines_v38.csv` |
| 23 | The ClinGen file lacked 104 curations, CFTR among them; every null variant of CFTR in the benchmark was called VUS; no aggregate showed it | MS:747-751, 1093, 1136-1138, 1187 |
| 23 | Counting BS2 beside BS1 or BA1 adds 638 calls in direction and calls two ClinVar-Pathogenic variants Benign | MS:962-971, 1094-1095 |
| 23 | One tool, one reviewer, one benchmark | MS:1103 |
| 24 | The build commands | `docs/ACMG_SETUP.md:88-312`, `docs/SA_DATABASES.md:20-30`, `scripts/README.md:8-10,38-49`, the usage text of `scripts/build-sa-databases.sh` (`--regions`, `--chroms`, the default chromosomes 1 to 22 and X) |
| 24 | No account, 198.4 GB, 3.4 GB (3,448,659,975 bytes in `data/benchmark/sa_db_complete`, measured on 8 Nov 2026, following its links; the table of `docs/ACMG_SETUP.md` on master lists region-restricted databases and totals 1.3 GB), ClinVar is re-released weekly, pin a dated release | `docs/ACMG_SETUP.md:9,120-121,128,154-159,497` |
| 24 | The small tier took 593 s, 0.9 GB of downloads (911,527,527 bytes), 304 MB of databases (304,092,314 bytes) | `scripts/build-sa-databases.sh` run 8 Nov 2026 (exit 0): the sizes of the files it wrote |
| 24 | The release archive holds `fastvep`, `fastvep-web` and three text files, not the scripts; the extractor needs `tabix` | the file list of the archive, 8 Nov 2026; `scripts/extract_gnomad_scores.py` |
| 25 | The check, its table, and the exit statuses | `docs/ACMG_SETUP.md:381-433`, `scripts/check_acmg_stack.py` run 8 Nov 2026 on the stack of the benchmark (exit 0; gnomAD 7/8, REVEL 3/8, SpliceAI 7/8 and PhyloP 7/8 annotated) and on a ClinVar-only directory (exit 1) |
| 25 | The eight test variants are on chromosomes 1, 2, 7, 10, 14 and 17, and eight of the nine sources pass when they answer at any one of them and RepeatMasker, through BP3, only when BP3 is evaluated on all eight: a stack for chromosome 22 alone fails, one for chromosome 17 alone passes | `scripts/check_acmg_stack.py:46-55,171-193`, and runs of the check on stacks of one chromosome (8 Nov 2026) |
| 26 | The command; the four calls and criteria | `docs/ACMG_SETUP.md:435-446`, run 8 Nov 2026 |
| 26 | A default run reports 26 of the 28 criteria; with `--pick` the call is made on the picked transcript | `docs/ACMG.md:7`, `docs/ACMG.md:776` |
| 26 | ClinVar of 27 June 2026 | MS:186 (Table 2) |
| 27 | The web page and the badges | `docs/ACMG.md:57-63`, `docs/API.md:154,261` |
| 27 | GRCh38 has all nine sources loaded and `/api/status` says available and ready; the other eight genomes do not; `"acmg": true` was answered with HTTP 200 and an `acmg` object on each consequence; the page offers the option ticked and enabled, with the hint that all 9 reference sources are loaded | `https://fastvep.org/api/status`, `/api/genomes` and `/api/annotate` at 18:58 UTC (`/api/status` again at 19:15 EST) and the page in a browser at 19:44 UTC, 8 Nov 2026 |
| 27 | The rule of the gate: all nine sources, with gnomAD for every autosome and X that the gene model has transcripts on, else HTTP 422 and the list of what is missing unless the operator allows partial stacks | `crates/fastvep-web/src/acmg.rs:11-18,64-67`, `crates/fastvep-annotate/src/coverage.rs:34-37` and `errors.rs` of the branch `feat/web-genomes` at e068377, the live build (it is in neither the tag nor master) |
| 27 | A request with `"acmg": true` for a genome without the stack is refused with HTTP 422 and an error that names the missing sources | `https://fastvep.org/api/annotate` for `mouse_grcm39_ensembl_115`, at about 19:27 EST, 8 Nov 2026 |
| 27 | What a request costs: 1,000 ClinVar variants with ACMG and pick took 5.0 s alone and 8.0 s each when four ran at once (500 variants a second in all); six at once were all answered, two after waiting for a permit; 10,000 variants over the genome were refused with HTTP 413 after 7,300, which is what 128 MiB holds (the source does not say whether that request asked for ACMG or pick) | `DEPLOYMENT:932-934` and `DEPLOYMENT:938-939` (a second instance of the same binary, beside the live service, on the host of 4 vCPU and 16 GB), 8 Nov 2026; the limit is `max_response_bytes` 134,217,728 in `/api/status` |
| 27 | The release's `fastvep-web` has no `acmg` object in `/api/status` and answers `"acmg": true` with a classification from whatever is loaded | the release's `fastvep-web` run on a copy of `tests/test.gff3`, 8 Nov 2026 |
| 28 | The research-use sentence | MS:1222 |
| 28 | The four bounds, and the hold-out of variants and almost not of genes | MS:77-81, 903, 1101 |
| 29 | The links, the Slack invitation and the citation | the header and footer of https://fastvep.org, 8 Nov 2026; `README.md:900-909`; `gh repo view` |
| 29 | Binaries for macOS (Apple silicon) and Linux | `README.md:66-70`, the assets of the release |
| 29 | The manuscript is unpublished and its declarations are still to be completed; its branch is not on master, and the commit that the footers cite is `2a293d0` | MS:1353-1369 |
| 30 | The trajectory: each version at its run of the time, and v37 and v38 on the final benchmark databases (column `databases` of `trajectory.csv`); v3 and v17 not plotted; what v1 to v6 are; the line of contradictions starts at v6 and v1's 31 rest on 74,355 variants (v38: 256,054) | MS data `trajectory.csv`, MS:233-237, 370-393 |
| 30 | The sentence about the versions before v37 (v37: 77.62% and 49 contradictions as run at the time, 79.07% and 52 on the final databases) | MS:371, 384, 1430 (the legend of Figure S1), and MS data `trajectory.csv` (the row of v37 is its run on the final databases, and its columns `*_of_the_time` are its run at the time) |
| 30 | 77.6% at v13 and 77.50% to 77.62% for the next 23 runs, 94 to 49, v10's 25 against v38's 51, and 6,632 | MS:381, 384, 386, 1110, 1190 |
| 31 | The flag table | `docs/VEP_MIGRATION.md:25-30,42-62,96-101,201-205`, `fastvep cache --help`, `README.md:825` |
| 32 | The options and the configuration fields; the default of `pm2_ad_af_threshold` is 0.00004 and 0.0 is the strict-absence reading | `docs/ACMG.md:41-54,70-75,83-118`, `docs/ACMG_SETUP.md:428`, `fastvep annotate --help`, MS:463-483 |

## Every command on a slide

All of these were run with the release binary (v0.5.0-rc.1, macOS, Apple silicon) on 8 November 2026 unless the last column says otherwise.
Where a run used a copy of the transcript cache of the Ensembl 115 GFF3 passed with `--transcript-cache`, it was so that the sidecar cache next to the GFF3 stayed as it was, and the result is the same as without the flag.

| Slide | Command | How it was checked |
|---|---|---|
| 7 | `fastvep filter -i annotated.vcf --filter "IMPACT is HIGH and AF < 0.01"` | run on two annotated variants: `[filter] CSQ format: 49 fields` and `0 of 2 variants passed filter` (that output has no AF), exit 0; the expression is `README.md:611` |
| 7 | the release's `fastvep-web --data-dir genomes` with two genome directories, then `GET /api/genomes`, `GET /api/status` and `POST /api/load-genome` | run: both genomes listed (HTTP 200), `allow_model_replacement` false, no `acmg` key, and HTTP 403 for the load (the server was started without `--allow-model-replacement`) |
| 8 | the download recipe with `VERSION=v0.5.0-rc.1` | run: the checksum printed `OK`, `fastvep --version` printed `fastvep 0.5.0-rc.1`, and the archive lists `fastvep`, `fastvep-web`, `LICENSE.md`, `CHANGELOG.md` and `README.md` |
| 9 | `fastvep annotate -i tests/test.vcf --gff3 tests/test.gff3 --hgvs --pick --output-format tab` | run on copies of the two files, identical to the repository's (`cmp`) |
| 11 | `fastvep annotate -i two.vcf --gff3 genes.gff3 --fasta genome.fa --hgvs --pick --output-format tab`, `vcf` and `json` | run on the Ensembl 115 GFF3 and FASTA |
| 12 | `mkdir -p sa`, `fastvep sa-build --source custom_vcf --name mypanel -i panel.vcf -o sa/mypanel` and the `annotate` that follows | run with a `panel.vcf` of two records that carry `AC`, `AN` and `AF`: `Wrote: sa/mypanel.osa2 (2 records)` and two `FV_MYPANEL` values; without `mkdir -p sa` the build stops with `No such file or directory` and exit 1 |
| 14 | `fastvep-web --bind 127.0.0.1 --port 8080 --gff3 genes.gff3 --fasta genome.fa --sa-dir sa/` and the `curl` request | run: a local server returned the same consequence and HGVS fields as fastvep.org for the same request; its supplementary field was the one database of its `--sa-dir` (`mypanel`) |
| 24 | `scripts/build-sa-databases.sh --assembly GRCh38 --out sa_databases` | run: the small tier, exit 0 in 593 s |
| 24 | `scripts/build-sa-databases.sh ... --sources gnomad,revel --chroms 22 --regions targets.bed` | not run, because the gnomAD and REVEL downloads are large; its flags are those of the script's usage text |
| 24 | `python3 scripts/extract_gnomad_scores.py --chrom 22 --regions targets.bed --out-dir extracts/` and the two `fastvep sa-build` commands | run on a 200-base region of chr17 (140 SpliceAI and 140 PhyloP records, two databases written) |
| 25 | `python3 scripts/check_acmg_stack.py --gff3 genes.gff3 --fasta genome.fa --sa-dir sa_databases/` | run on the stack of the benchmark (`data/benchmark/sa_db_complete`): all nine OK, exit 0, with gnomAD 7/8, REVEL 3/8, SpliceAI 7/8 and PhyloP 7/8 annotated, in 2 s; on a ClinVar-only directory: INCOMPLETE, exit 1 (the shell printed the statuses); a stack of chromosome 22 alone ends INCOMPLETE and one of chromosome 17 alone passes |
| 26 | `fastvep annotate -i four.vcf -o four.acmg.vcf --gff3 genes.gff3 --fasta genome.fa --sa-dir sa_databases/ --acmg --pick --output-format vcf` | run on the four variants and the stack of the benchmark (`data/benchmark/sa_db_complete`) |
| 27 | the release's `fastvep-web --gff3 tests/test.gff3`, `GET /api/status` and `POST /api/annotate` with `"acmg": true` | run on a local port: HTTP 200 and no `acmg` key in the status, and HTTP 200 with the class `UncertainSignificance` and 26 criteria for the one variant tried |
| 27 | `GET /api/status`, `GET /api/genomes` and one `POST /api/annotate` with `"acmg": true` to https://fastvep.org | one request each at 18:58 UTC on 8 Nov 2026, and the page loaded once in a browser at 19:44 UTC: `acmg` available and ready for GRCh38, not available for the other genomes, and HTTP 200 with an `acmg` object per consequence; `GET /api/status` again at 19:15 EST, and a `POST /api/annotate` with `"acmg": true` for `mouse_grcm39_ensembl_115` at about 19:27 EST, which answered HTTP 422 and named the missing sources |
| 29 | `fastvep annotate -i your.vcf --gff3 genes.gff3 --fasta genome.fa --hgvs --pick` | the flags were run in slide 11's command and exist in the help |
| 31 | `fastvep cache --gff3 genes.gff3 --fasta genome.fa -o genes.cache` and `--transcript-cache` | run without `--fasta` on the test fixture, then annotated from the cache |
| 32 | the trio and functional-evidence options | not run, because they need a trio VCF and a curated file; the flags exist in the help and the shape is `docs/ACMG.md:44-54` |
| 32 | a configuration file with `bs1_af_threshold`, `pm2_ad_af_threshold` and `bp7_max_intron_offset` | run with the first two on the stack of the benchmark: accepted, and for the four variants of slide 26 it changed no byte of the VCF; the third is in `docs/ACMG.md:118` |

## What was checked

- Layout: every slide at 1920 by 1080, 1280 by 720 and 390 by 844 pixels, in light and dark, was screenshotted and looked at, and a script checked each of them for content outside the slide, content that spills out of its box, a box shorter than its content, a gap between blocks that the style sheet declares and the page lost, overlapping blocks, a heading on more than one line (on the large sizes), a single word alone on the last line of a heading, code in a size under 24 px on the 1080p canvas, and any sideways scroll of the page or the slide.
- Accessibility: axe-core (WCAG 2 A and AA, 2.1 and 2.2 AA, and best practice) found no violation on any slide at the two sizes in both schemes, nor on the overview, the open notes and the list of keys: 144 scans.
- Behaviour: 63 checks of the keys, the address, the clicks on the thirds, the swipe, the overview, the notes, the list of keys, the theme, the copy buttons (the clipboard API and the fallback), the phone layout and the reduced-motion setting.
- Print: Chrome's print engine wrote PDF files of the deck, with and without the notes and from windows of 1280 and 1920 pixels, and every page was checked: 33 pages of 13.333 by 7.5 inches (13.333 by 10.5 with the notes), the text of every page at the scale that the print style sheet means, and the pages rasterized and looked at.
  Width queries of the phone layout are written for the screen only, because Chrome evaluates a width query in print at the paper width of its print settings (816 px for Letter) and an unscoped one turned the printed slides into the phone layout.
  The slides are also contained in print, because Chrome shrank every page of the deck to fit the visually hidden text of slide 18, which runs past its canvas.
- No request: the deck was opened from a clean browser profile with every host unresolvable and a dead proxy, and stepped end to end: two requests, one for its own file and one for its inline logo (a `data:` URI), and none that left the page.
- Facts: independent reviewers audited every claim of every slide and of every note, and this file, against its source, recomputed every number and ran the commands again.
  VEP's 6,662 s of slide 4 is 6,661.46 s before rounding, which the repository's README and CSV round to 6,662, and the deck keeps the repository's figure.
  The ratios of the chr22 table are those of the unrounded medians.
  The ACMG numbers of slides 17 to 23, 28 and 30 were recomputed from the data tables of the manuscript at commit `2a293d0`, the runs of slides 24 to 26 from their outputs, and the cost of a request on slide 27 from `DEPLOYMENT.md` and `/api/status`.
  A further reviewer audited the slides that changed, their notes and this file against their sources and recomputed every number: it found no wrong number, and what it found in wording, scope and traceability was corrected.
  A last, short review read the slides, notes, footers and rows of this file that depend on the manuscript against its final commit, recomputed every plotted point and number of slide 30, and found no wrong number and no wrong slide reference; its eight findings, of wording, citation and one chart marker, were corrected.

## What is not in the deck, and why

- Usage statistics of fastvep.org: the README and the site give them, they change daily, and no claim needs them.
- An author's email address: the deck sends people to GitHub, to Slack and to the lab's page.
- A comparison of fastVEP's ACMG calls with another tool: the manuscript made none.
- The number of lines or of tests of the classifier crate: the manuscript gives them for the tag, and they were not measured for this deck.
- A citation for the ACMG manuscript, which is unpublished.
- How the web host behaves over hours of mixed traffic, or with the nginx rate limit at its limit: `DEPLOYMENT.md` says that these have not been measured.

## Which numbers come from which commit

- The ACMG numbers of slides 17 to 23, 28 and 30 are the manuscript's at commit `2a293d0` of the branch `docs/acmg-scheme-and-versions`: its text, and the `_v38.csv` and other files of its `manuscript/acmg/data/` directory.
- The runs of slides 24 to 26 (the size of the benchmark's stack, the check of the stack and the four calls) were made on 8 November 2026 with the release binary on `data/benchmark/sa_db_complete`, the stack that the benchmark of that commit reads.
- The cost of a request on slide 27 is `DEPLOYMENT.md` of the branch `feat/web-genomes` at commit `95e4d2a`, and the rest of that slide is the live site on 8 November 2026.
- The section, table and line numbers in this file are those of commit `2a293d0`.
  If the manuscript changes, its sections move: read each claim again in the new text and do not follow the number given here.

If the manuscript changes, check these slides against its text and its data tables, and point the footers at the new commit and file names:

| Slide | What to check |
|---|---|
| 16 | The wording of the problem and of the contribution (abstract and Section 1) |
| 17 | The counts of versions, rounds and variants; which thresholds were swept and which numbers were not varied; the description of a run |
| 18 | The nine sources and what they feed (Table 2); the grid of criteria (`criteria_status.csv` is read from the engine; 15 of 28 is by design) |
| 19 | The four pairs of the measured box (recall and contradictions of the two replays and of the run), and the reason for the split |
| 20 | Every number: 64.1%, 79.3%, the intervals, 51, 0.0076%, 62; 116,078, 72.6%, 2.3 times and its interval; 4,862 genes and 157; 98.4% and 92% |
| 21 | Every segment of the bar and its percentage, 93.6% of 139,325, 10,694, 5,184, 20,072, 286,815 of the 431,727; 104 curations is an input |
| 22 | Every row (43.8%, 64.0%, 72.8%, 85.6%, 79.3%; 0, 813, 865, 1,817, 51), 6.3 points, 36 times, 4,039 and 67,733, the expected counts 15 and 1.7, 62 and three |
| 23 | Each of the six lessons and its numbers: 10.6 points and 372; 21 and 30 of 64; the thresholds named and the bound of 1,000 calls; 85.6%, 1,817, 79.3% and 51; 104 curations and CFTR; 638 calls and two variants |
| 24 | "The benchmark's whole stack: 3.4 GB" (measure the directory again) |
| 25 | The table of the check on the stack of the benchmark (7/8, 3/8, 7/8 and 7/8 annotated for gnomAD, REVEL, SpliceAI and PhyloP) and "3.4 GB" in the notes |
| 26 | The four calls, their criteria and points (PP3, BP4 and BP7 read SpliceAI and PhyloP), the ClinVar release date, and "26 criteria" |
| 28 | The four bounds (Introduction) and the hold-out of variants |
| 29 | The commit that the footers cite, and whether the branch has been merged |
| 30 | The points of both charts and the figures quoted beside them (77.6% at v13; 77.50% to 77.62%; 94 to 49; v37 at 77.62% as run at the time and 79.07% on the final databases; 25 against 51; 6,632; 31 and 74,355; 256,054), and which versions are on the final databases |
| 32 | The pointer to Table 3 (how each threshold was set and what the steps on either side buy and cost); the defaults in the file are the classifier's own and do not change |

Slides 1 to 15, 27, 31 and 33 carry no number of the manuscript.
Slide 27 carries the cost of a request from `DEPLOYMENT.md`, which changes if the host or the limits do.

## Slides to look at again if the site changes

- Slides 7, 13, 14, 15 and 27 describe fastvep.org as it was at 18:58 UTC on 8 November 2026 (slide 7: all nine genomes served at once): GRCh38 has all nine sources and offers ACMG-AMP classification, the other eight genomes do not, and the web build is `0.5.0-rc.1+web.e068377`.
  If another genome gets the nine sources, or the build changes, read `/api/status`, `/api/genomes` and the page again, and the right card and the excerpt of slide 27 with them.
- Slide 27 gives the cost of a request from `DEPLOYMENT.md`: if the host, the limits (`/api/status`: 4 parallel requests, 90 s, 128 MiB) or the build change, measure again.
- Slide 29 holds the Slack invitation, which is a link that can expire.

## Keeping it current

The deck is plain HTML: each slide is a `section` with a `div` for its canvas and a `div` that holds its notes, and the styles and the script are in the `style` and `script` elements.
When a release replaces v0.5.0-rc.1, search the file for `rc.1`, and replace the numbers of slides 1 to 6 and 15 from the new `benchmarks/results/` directory.
When the manuscript changes, take the numbers of the slides listed under "Which numbers come from which commit" from its `data/` tables again.
Sizes in the styles are written as `calc(N * var(--u))`, where `--u` is one pixel on the 1920 by 1080 canvas and half a pixel on a phone.
Keep the width queries of the phone layout under `screen` and keep `contain` on the printed slide, or the print breaks (see What was checked).
