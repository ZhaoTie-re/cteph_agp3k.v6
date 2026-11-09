# CTEPH with DVT vs CTEPH without DVT (cases only): association along the HLA genes

**Covers:** `region.allele.png`, `region.residue.png` in `01.region/`

## The question this figure answers

Where along the typed HLA genes do the allele and the amino-acid residue markers differ between the two groups compared (CTEPH with DVT vs CTEPH without DVT (cases only)), and does the picture hold across the three cohorts and both models?

## Panels

**(region.allele.png) rows narrow, intermediate, full mainland; columns fixed model, random model**

Each point is one tested allele, at its place in one marker order shared by all panels: genes in GRCh38 chromosomal order (alternate genes shaded), then residues by IMGT position and residue, alleles by their numeric fields. Height is −log10 P (fixed: plink2 Wald P; random: SAIGE unadjusted score P). Three family-wise thresholds (0.05 over the class's tests in that cohort): solid red = LD-adjusted, 0.05 / Me with Me from GEC v0.2 (primary); dashed grey = exact max-T (validation); dotted = conventional Bonferroni, 0.05 / tests (sensitivity). Red points are significant at the GEC line; the three smallest of them per panel are labelled.

**(region.residue.png) rows narrow, intermediate, full mainland; columns fixed model, random model**

Each point is one tested amino-acid residue, at its place in one marker order shared by all panels: genes in GRCh38 chromosomal order (alternate genes shaded), then residues by IMGT position and residue, alleles by their numeric fields. Height is −log10 P (fixed: plink2 Wald P; random: SAIGE unadjusted score P). Three family-wise thresholds (0.05 over the class's tests in that cohort): solid red = LD-adjusted, 0.05 / Me with Me from GEC v0.2 (primary); dashed grey = exact max-T (validation); dotted = conventional Bonferroni, 0.05 / tests (sensitivity). Red points are significant at the GEC line; the three smallest of them per panel are labelled.

## Interpretation

At the LD-adjusted threshold, at most 0 allele and 0 residue tests are significant in any one scan. A residue position with two residues gives two mirror columns with the same P (one test), and neighbouring positions on one haplotype carry nearly the same signal; the number of points above a line is not the number of independent signals (see _conditional/). Cross-cohort consistency is read in the random model: the tables below count, per threshold, the markers significant in 3/3, 2/3 and 1/3 ("Single cohort") cohorts with the same direction of effect.

## Values in this rendering

| quantity | value |
|---|---|
| region.allele: narrow fixed: markers / tests | 143 / 143 |
| region.allele: narrow fixed: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: narrow fixed: min P | 9.94e-03 (DQB1*05:01) |
| region.allele: narrow random: markers / tests | 143 / 143 |
| region.allele: narrow random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: narrow random: min P | 8.10e-03 (DQB1*05:01) |
| region.allele: intermediate fixed: markers / tests | 144 / 144 |
| region.allele: intermediate fixed: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: intermediate fixed: min P | 6.31e-03 (DQB1*05:01) |
| region.allele: intermediate random: markers / tests | 144 / 144 |
| region.allele: intermediate random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: intermediate random: min P | 4.79e-03 (DQB1*05:01) |
| region.allele: full fixed: markers / tests | 142 / 142 |
| region.allele: full fixed: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: full fixed: min P | 8.45e-03 (DQB1*05:01) |
| region.allele: full random: markers / tests | 142 / 142 |
| region.allele: full random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: full random: min P | 7.01e-03 (DQB1*05:01) |
| region.residue: narrow fixed: markers / tests | 672 / 457 |
| region.residue: narrow fixed: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.residue: narrow fixed: min P | 4.45e-04 (DQB1:-27:A) |
| region.residue: narrow random: markers / tests | 672 / 457 |
| region.residue: narrow random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.residue: narrow random: min P | 3.59e-04 (DQB1:-27:A) |
| region.residue: intermediate fixed: markers / tests | 681 / 463 |
| region.residue: intermediate fixed: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.residue: intermediate fixed: min P | 5.46e-04 (DQB1:-27:A) |
| region.residue: intermediate random: markers / tests | 681 / 463 |
| region.residue: intermediate random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.residue: intermediate random: min P | 4.61e-04 (DQB1:-27:A) |
| region.residue: full fixed: markers / tests | 679 / 462 |
| region.residue: full fixed: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.residue: full fixed: min P | 1.25e-03 (DQB1:-27:A) |
| region.residue: full random: markers / tests | 679 / 462 |
| region.residue: full random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.residue: full random: min P | 1.13e-03 (DQB1:-27:A) |

## Full statistics

**Scan summary (one row per cohort × class × model)**

| class | cohort | model | samples | cases | controls | markers_tested | tests | me_gec | p_gec | p_maxt | p_bonf | tests_gec | tests_maxt | tests_bonf | lambda_gc | min_p | top_marker | top_or |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| allele | narrow | fixed | 290 | 150 | 140 | 143 | 143 | 111.3 | 4.490e-04 | 4.309e-04 | 3.497e-04 | 0 | 0 | 0 | 0.965 | 0.00994 | DQB1*05:01 | 3.566 |
| allele | narrow | random | 290 | 150 | 140 | 143 | 143 | 111.3 | 4.490e-04 | 4.309e-04 | 3.497e-04 | 0 | 0 | 0 | 1.142 | 0.008098 | DQB1*05:01 | 3.159 |
| allele | intermediate | fixed | 296 | 155 | 141 | 144 | 144 | 112.3 | 4.454e-04 | 4.260e-04 | 3.472e-04 | 0 | 0 | 0 | 0.836 | 0.006312 | DQB1*05:01 | 3.794 |
| allele | intermediate | random | 296 | 155 | 141 | 144 | 144 | 112.3 | 4.454e-04 | 4.260e-04 | 3.472e-04 | 0 | 0 | 0 | 1.153 | 0.004794 | DQB1*05:01 | 3.305 |
| allele | full | fixed | 304 | 159 | 145 | 142 | 142 | 110.5 | 4.523e-04 | 4.362e-04 | 3.521e-04 | 0 | 0 | 0 | 1.062 | 0.00845 | DQB1*05:01 | 3.355 |
| allele | full | random | 304 | 159 | 145 | 142 | 142 | 110.5 | 4.523e-04 | 4.362e-04 | 3.521e-04 | 0 | 0 | 0 | 1.12 | 0.007008 | DQB1*05:01 | 3.019 |
| residue | narrow | fixed | 290 | 150 | 140 | 672 | 457 | 185.2 | 2.700e-04 | 2.925e-04 | 1.094e-04 | 0 | 0 | 0 | 1.453 | 4.452e-04 | DQB1:-27:A | 0.315 |
| residue | narrow | random | 290 | 150 | 140 | 672 | 457 | 185.2 | 2.700e-04 | 2.925e-04 | 1.094e-04 | 0 | 0 | 0 | 1.438 | 3.592e-04 | DQB1:-27:A | 0.344 |
| residue | intermediate | fixed | 296 | 155 | 141 | 681 | 463 | 186.6 | 2.680e-04 | 2.943e-04 | 1.080e-04 | 0 | 0 | 0 | 1.599 | 5.460e-04 | DQB1:-27:A | 0.332 |
| residue | intermediate | random | 296 | 155 | 141 | 681 | 463 | 186.6 | 2.680e-04 | 2.943e-04 | 1.080e-04 | 0 | 0 | 0 | 1.618 | 4.607e-04 | DQB1:-27:A | 0.36 |
| residue | full | fixed | 304 | 159 | 145 | 679 | 462 | 186.5 | 2.681e-04 | 2.902e-04 | 1.082e-04 | 0 | 0 | 0 | 1.448 | 0.001246 | DQB1:-27:A | 0.375 |
| residue | full | random | 304 | 159 | 145 | 679 | 462 | 186.5 | 2.681e-04 | 2.902e-04 | 1.082e-04 | 0 | 0 | 0 | 1.529 | 0.001131 | DQB1:-27:A | 0.398 |

**Alleles: markers by cross-cohort consistency (random model)**

| threshold | 3/3 | 2/3 | 1/3 |
|---|---|---|---|
| LD-adjusted (GEC) (primary) | 0 | 0 | 0 |
| Exact (max-T) (validation) | 0 | 0 | 0 |
| Conventional Bonferroni (sensitivity) | 0 | 0 | 0 |

**Amino-acid residues: markers by cross-cohort consistency (random model)**

| threshold | 3/3 | 2/3 | 1/3 |
|---|---|---|---|
| LD-adjusted (GEC) (primary) | 0 | 0 | 0 |
| Exact (max-T) (validation) | 0 | 0 | 0 |
| Conventional Bonferroni (sensitivity) | 0 | 0 | 0 |

## How to read it

1. Find the gene band with the tallest points.
2. Read across a row: does the random model (right) reproduce the fixed model (left)?
3. Read down a column: does the peak hold as controls are added (narrow → full)?
4. Look up the exact P, OR and counts of any marker in _comparison/assoc_all.tsv.

## What this figure does *not* establish

- Which of several correlated markers carries the signal: see the conditional analyses (results/_conditional/).
- The three cohorts share almost all their cases and this set has no controls, so the three cohorts are nearly the same analysis; agreement across them is expected, not a replication.
- Both groups are cases, so the case–control platform split does not enter. But the sequencing platform also differs between cases with and without DVT (docs/STUDY_NOTES.md), so a platform-dependent typing difference would still look like an association; platform is not a covariate here.

---

Methods and rationale: [`METHODS.md`](../../../../docs/METHODS.md)
