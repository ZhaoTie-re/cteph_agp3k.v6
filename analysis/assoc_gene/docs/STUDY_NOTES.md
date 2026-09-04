# Study notes — `assoc_gene`

The measured facts the design rests on. Each was established on the data or in the rvtest source, and
each is enforced somewhere in the pipeline or in `verify.sh`; the pointer is given.

## 1. The phenotype coding defect

`model_inputs/<cohort>/phenotype/pheno.tsv` is 0/1 today; it was 1/2 when an earlier scan ran and was
rewritten in place. rvtest reads a phenotype value of 0 as *missing* (`DataLoader.cpp:1248`): a 0/1
column passes its binary check, maps 0 → −1, drops every negative row, logs `WARN There are no case!`
and **exits 0**. Every downstream p is then a null scan with plausible magnitudes. Fixed in
`prep_pheno_cov.py` (normalise, hand rvtest 1/2) and gated in the RVTEST task, which re-reads rvtest's
"Loaded *n* cases, *m* controls, *k* missing" line against `pheno_coding.json`.

## 2. `--setFile` renames the first column

With `--setFile` rvtest's first output column is `Range` (the set name) and the second `RANGE` (the
consolidated range list, ~30 kB per large gene); with `--geneFile` the first is `Gene`. Read by
position or through a case-folding path, the 30 kB cell lands where the set name belongs and the scan
joins to nothing. `merge_rvtest.py` asserts the header and refuses a `Gene` table — which would also
mean refFlat spans were tested, not the map.

## 3. The span-overlap over-count is exact

A 1-bp tabix range at position q matches a record at p with `len(REF) = L` whenever `p ≤ q ≤ p + L − 1`.
rvtest therefore counts a deletion twice when a later mapped variant sits strictly inside its span
(`p < q ≤ p + L − 1`; two records at the same position do not over-count because rvtest consolidates
identical ranges). On the reference data the predictor matches rvtest's `NumVar` on **81,746 / 81,746**
genes across nine cohort × stratum cells, zero counter-examples. The first formula tried (`p ≤ q`)
matched 3 of 150; strictness is the whole difference. `scan.py`, `verify.sh` §10.

## 4. Map composition and denominators

MinNumVar = 2, per-cohort minAC = 2.

| cohort | N (case / control) | `moderate_high` mapped / dropped | `low_moderate_high` mapped / dropped |
|---|---|---|---|
| narrow_mainland | 2,165 (418 / 1,747) | 10,832 / 3,969 | 14,557 / 3,911 |
| intermediate_mainland | 2,477 (428 / 2,049) | 11,436 / 3,771 | 15,009 / 3,848 |
| full_mainland | 3,068 (438 / 2,630) | 12,378 / 3,437 | 15,691 / 3,707 |

`full_mainland` / `moderate_high`: 75,787 variants, 76,845 (variant, gene) rows, 12,378 sets, Bonferroni
0.05 / 12,378 = 4.04e-06. The `high` stratum (509–740 genes, median 2 variants, no call anywhere) was
dropped from the design.

## 5. The gene families are not nested across cohorts

PPOX, ZBP1 and ZNF30 have ≥ 2 qualifying variants in `narrow_mainland` and fewer in `full_mainland`
(minAC is applied per cohort). A cross-cohort table that read "absent" as "not significant" would be
wrong for them; `evidence.tsv` keeps `not_in_map` as a state, decided from `gene_index_all.tsv`, and
`verify.sh` §17 counts those rows against the index sizes.

## 6. The cohorts are nested; 20 cases and 883 controls separate the extremes

narrow ⊂ intermediate ⊂ full with zero exceptions. `full_mainland` has 20 more cases (438 vs 418) and
883 more controls (2,630 vs 1,747) than `narrow_mainland`, so the three differ almost entirely in
controls. "Called in all three cohorts" therefore means the signal survives adding controls. It is a
sample-selection sensitivity analysis. Every summary caption and figure sidecar states this.

## 7. λ<sub>GC</sub> on a discrete statistic

Over the six cohort × stratum cells, CMC λ is 1.33–1.50 over genes with `n_var = 2` and 1.09–1.14 over
`n_var ≥ 3`; SKAT-O 1.17–1.37 and 1.16–1.23. A null simulation on the real MAC distribution gives
λ ≈ 0.58 for an exact test, and Fisher's exact test on the same counts gives ≈ 0.6. The median-based λ is inflated for two-variant genes by the normal
approximation near p = 0.5, not by association, and the tail is unaffected the same way. λ is reported
overall and split at `n_var = 2 / ≥ 3` in `scan_qc.tsv` and the figure READMEs; the figures draw the
overall value only, and it is never used as a target. `scan.py`, `docs/METHODS.md` §4.4.

## 8. The per-gene counts, verified on STBD1

`full_mainland` / `moderate_high`, 3 variants (chr4:76306990 G>C splice donor, chr4:76309487 C>G and
chr4:76309559 G>T missense): carriers 10 (2.28 %) / 7 (0.27 %), MAC 10 / 7, MAF 0.011 / 0.001, SKAT-O
ρ 0.3, CMC OR 9.60 (3.57–25.85), CMC score p 7.2e-08. The bed decode agrees with `plink2 --freq
counts` on all 151,574 per-variant checks of that cell, the set files are byte-identical to the
previous implementation's, and the MAC columns agree on all 12,378 sets. OR = exp(β) from cmcWald
reproduces the statsmodels refit (9.601) to three figures.

## 9. Multi-site carriers

A sample can carry a minor allele at several sites of one gene. CMC counts it once; the cumulative MAC
counts every site. In `full_mainland` / `moderate_high` this affects 19,300 (gene, sample) pairs in
2,453 genes; GTF3C3 has 10 case carriers holding 21 minor alleles (4 case samples carry 3–4 of its 14
sites). The naive bound `MAC ≤ 2 · carriers` is therefore false; the asserted bound is
`carriers ≤ MAC ≤ 2 · n_var · carriers`. The summary table carries `multi_site_carriers` per block and
`multi_carrier_detail.tsv` lists the samples for the tiered genes (`emit_test_inputs.py`,
`robust_genes.py`, `verify.sh` §13).

## 10. Minor-allele orientation

192 of 75,787 variants in `full_mainland` / `moderate_high` are ALT-major. The flip is decided once on
the pooled counts and applied to both groups, as rvtest's own collapse does; a per-group decision could
make the two groups disagree about which allele is minor. Computed, never assumed: one un-flipped
ALT-major variant would contribute 2N − AC to its gene.

## 11. rvtest facts pinned by the smoke test (`verify.sh` §6)

Output names `<out>.CMC.assoc`, `<out>.CMCWald.assoc`, `<out>.SkatO.assoc`; cmcWald writes `1 + n_covar`
rows per gene with the burden term first and no labels (`Model.h:961`); SkatO's header is `Range RANGE
N_INFORMATIVE NumVar NumPolyVar Q rho Pvalue`; a constant covariate is refused; `--impute` accepts only
`mean`, `hwe`, `drop` and an unrecognised value silently exits 0 (`Main.cpp:1002-1013` has no `else`).

## 12. Tiers on the reference data (BH-called, MinNumVar = 2)

| stratum | Tier 1 | Tier 2 | Tier 3 |
|---|---|---|---|
| moderate_high | OSBP, STBD1 | JAK2, PIK3C2G, RPH3A | NKD2, P2RY2, PES1, SCFD1 |
| low_moderate_high | — | PIK3C2G, RPH3A, STBD1 | CD248, PES1, SCAMP3, SCFD1 |

To be re-read from `results/_comparison/robust_genes/tiers.tsv` after each run.

## 13. Limitations, in order

1. **Platform confound.** All cases DNBSeq / NovaSeq, all controls HiSeq X 15×, zero overlap; unfixable by
   covariate and present in every p-value. Checked elsewhere in the project not to distort allele
   frequencies, but a gene here is a candidate, not a finding.
2. **Not replication.** §6.
3. **Nested strata.** A gene tiered in both strata is one finding tested twice on overlapping sets.
4. **Set-level statements only.** A gene p says nothing about which variant carries it; the per-gene
   figures and `variant_detail.tsv` are descriptive.
5. **λ is a coarse instrument** over ~10<sup>4</sup> non-independent, discrete tests (§7).
6. **cmcWald degeneracy.** A gene with no minor allele in one group has no OR; shown as `—`, noted in
   `tiers.tsv`.
7. **Tier semantics under a `--Cohorts` subset.** Tier 1/2 mean "every listed cohort"; the METHODS
   wording assumes the full three-cohort run.
