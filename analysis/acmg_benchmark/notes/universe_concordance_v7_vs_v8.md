# Universe-wide concordance: before vs after the ACMG fixes

ClinVar 2★+ benchmark, **673,660 variants**, ClinVar = truth.
- **v7** = pre-fix baseline (engine before commit 91e4ceb)
- **v8** = fixed engine (91e4ceb) + rebuilt ClinVar cache with AF backstop

| Metric | v7 (before) | v8 (after) | Change |
|---|---|---|---|
| Exact-match rate | 58.7% | **61.3%** | +2.6 pts (+17,465 variants) |
| Same-direction rate | 70.8% | **73.6%** | +2.8 pts (+18,972 variants) |
| Opposite-direction (worst errors) | 313 (0.046%) | **122 (0.018%)** | −191 (−61%) |
| No-call | 1 | 0 | — |
| Discordant by exact | 277,993 | **260,528** | −17,465 |
| Discordant by same-direction | 196,555 | **177,584** | −18,971 |

**Overall discordance decreased** on every measure. Opposite-direction calls
(benign↔pathogenic, the dangerous ones) fell by 61%.

Per-class exact-match count (v7→v8): Pathogenic +255, Likely_pathogenic +1,394,
VUS +6,467, Likely_benign +7,009, Benign +2,340.

One conservative trade-off: Pathogenic *same-direction* dropped by 582
(~0.7% of pathogenics) — genuine pathogenics that PVS1-downgrade / PM2-removal
moved to VUS. Dwarfed by the gains, and Pathogenic opposite-errors fell 168→67.

---

## Version-matched confirmation (v9)

Re-run with truth **and** annotation both from the ClinVar **2026-06-27**
release (`clinvar_2star.vcf` + truth table regenerated from it), 684,160 variants:

| Metric | v9 (2026-06-27, fixed engine) |
|---|---|
| Exact-match | **61.2%** (418,966) |
| Same-direction | **73.6%** (503,221) |
| Opposite-direction | **130 (0.019%)** |
| No-call | 0 |

Essentially identical to v8 (61.3% / 73.6% / 0.018%), so the fix improvement is
stable across ClinVar versions — not an artifact of the April snapshot.

Caveat: `clinvar.osa`, `clinvar_protein.oga`, and the truth table are fully
refreshed to 2026-06-27, but the per-chrom gnomAD/SpliceAI/PhyloP `.osa` extracts
are still limited to the April `clinvar_2star_regions.bed` (5 kb-padded). The
~10.5k net-new 2★ variants overwhelmingly fall in already-curated loci, so this
did not move the aggregate (opposite stayed ~0.02%, no-call 0); rebuild those
extracts from a fresh region BED for a fully pristine release-matched run.
