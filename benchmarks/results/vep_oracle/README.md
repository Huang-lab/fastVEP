# What Ensembl VEP 116.3 wrote for the checks of this release

These are the outputs of the scripts in `validation/` that run the pinned image of Ensembl VEP 116.3 (`--gff`, the default statistics, tab or JSON output) beside the binary of the release, over the same gene models and reference.
Each file starts with the command that made it.
`manuscript/verify_numbers.py` reads them, and checks that every number that `README.md`, `manuscript/manuscript.md`, `docs/STRUCTURAL_VARIANTS.md` and `docs/VEP_DIVERGENCE.md` quote from them is the one here.

| File | Made by | What it holds |
|------|---------|---------------|
| `sv_worlds.txt` | `validation/run_sv_world.py` | nine generated worlds of structural variants, the seed of the golden fixture and eight others: for each, the rows of VEP and of fastVEP, whether they are the same rows in the same order, and the cells of 23 columns that differ (none but the 5 to 7 of divergence 11 for each) |
| `json_most_severe_consequence.txt` | `fastvep annotate --output-format json` and VEP's `--json` | the three records of the whole-genome HG002 JSON whose `most_severe_consequence` changed with the ranks of release 116.3: VEP writes the donor term for each, as 0.6.0 does and master did not |
| `flag_pick_chr22.txt` | `validation/compare_pick.py` | the three flagging forms of `--pick` over chromosome 22 of HG002: how many picks are another transcript than VEP's, why (the ninth criterion, a CCDS that only fastVEP reads, a transcript type that VEP skips), and the rows of the variants that no transcript is beside |
| `sv_buffer_size.txt` | `validation/sv_buffer_size.py` | VEP with its own buffer of 5,000 records and with `--buffer_size 1` over file 09 of the vep-rs set (three runs each), the structural variants of ClinVar and a 1-in-4 sample of gnomAD-SV: its rows for a record depend on the records read with it |

The VEP side of every one of them is the program run on this machine, in Docker; nothing in these files was worked out by hand.
