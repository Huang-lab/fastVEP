The two CSV files `vep116_agreement.csv` and `vep116_hg002_columns.csv` are not written by `manuscript/build_data.py`: no run of this release made them.
Ensembl VEP 116.3 ran in Docker on 27 October 2026, and the cell counts are those that `docs/VEP_DIVERGENCE.md`, Part 3 ("the tab output, against VEP 116.3"), records for the build at commit 5048fc5.
They are the rows of the first table of that Part and the sentence "On HG002 the 3,851 cells are in four columns" with the table of causes below it.
The document names 21 of the 27 columns compared on HG002; the other six are not named there, and the 3,851 differing cells are all in the four HGVS columns that the document lists, so they agree on every row.

The release candidate v0.5.0-rc.1 writes the same bytes as that build on the 173 variants of the VEP example set, the stratified ClinVar sample of 6,600 variants, the 1-in-5 sample of ClinVar 2-star+, the 1-in-200 sample of HG002, the 400 in-frame deletions and the whole HG002 call set (`validation/compare_builds.py`, `benchmarks/results/v0.5.0-rc.1/identity_vs_5048fc5/identity.tsv`), so these counts are those of its output on those inputs.
The generated worlds were not re-run: their fixtures, with the output that VEP 116.3 wrote for them, are in `crates/fastvep-cli/fixtures/`, and the golden-file tests that compare fastVEP with that output passed in the test run of this release (`benchmarks/results/v0.5.0-rc.1/tests_executed.csv`).
`manuscript/verify_numbers.py` checks that every number of the two files is in `docs/VEP_DIVERGENCE.md`.
To measure them again, run `validation/run_tab_world.sh` and the commands of Part 3 with Docker, and replace the counts.
