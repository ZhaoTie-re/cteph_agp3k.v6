# CTEPH (all cases) vs controls: association along the HLA genes

**Covers:** `region.allele.png`, `region.residue.png` in `01.region/`

## The question this figure answers

Where along the typed HLA genes do the allele and the amino-acid residue markers differ between the two groups compared (CTEPH (all cases) vs controls), and does the picture hold across the three cohorts and both models?

## Panels

**(region.allele.png) rows narrow, intermediate, full mainland; columns fixed model, random model**

Each point is one tested allele, at its place in one marker order shared by all panels: genes in GRCh38 chromosomal order (alternate genes shaded), then residues by IMGT position and residue, alleles by their numeric fields. Height is −log10 P (fixed: plink2 Wald P; random: SAIGE unadjusted score P). Three family-wise thresholds (0.05 over the class's tests in that cohort): solid red = LD-adjusted, 0.05 / Me with Me from GEC v0.2 (primary); dashed grey = exact max-T (validation); dotted = conventional Bonferroni, 0.05 / tests (sensitivity). Red points are significant at the GEC line; the three smallest of them per panel are labelled.

**(region.residue.png) rows narrow, intermediate, full mainland; columns fixed model, random model**

Each point is one tested amino-acid residue, at its place in one marker order shared by all panels: genes in GRCh38 chromosomal order (alternate genes shaded), then residues by IMGT position and residue, alleles by their numeric fields. Height is −log10 P (fixed: plink2 Wald P; random: SAIGE unadjusted score P). Three family-wise thresholds (0.05 over the class's tests in that cohort): solid red = LD-adjusted, 0.05 / Me with Me from GEC v0.2 (primary); dashed grey = exact max-T (validation); dotted = conventional Bonferroni, 0.05 / tests (sensitivity). Red points are significant at the GEC line; the three smallest of them per panel are labelled.

## Interpretation

At the LD-adjusted threshold, at most 1 allele and 3 residue tests are significant in any one scan. A residue position with two residues gives two mirror columns with the same P (one test), and neighbouring positions on one haplotype carry nearly the same signal; the number of points above a line is not the number of independent signals (see _conditional/). Cross-cohort consistency is read in the random model: the tables below count, per threshold, the markers significant in 3/3, 2/3 and 1/3 ("Single cohort") cohorts with the same direction of effect.

## Values in this rendering

| quantity | value |
|---|---|
| region.allele: narrow fixed: markers / tests | 143 / 143 |
| region.allele: narrow fixed: tests significant GEC / max-T / Bonferroni | 1 / 1 / 0 |
| region.allele: narrow fixed: min P | 3.74e-04 (B*40:02) |
| region.allele: narrow random: markers / tests | 143 / 143 |
| region.allele: narrow random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: narrow random: min P | 5.09e-04 (B*40:02) |
| region.allele: intermediate fixed: markers / tests | 144 / 144 |
| region.allele: intermediate fixed: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: intermediate fixed: min P | 1.78e-03 (B*40:02) |
| region.allele: intermediate random: markers / tests | 144 / 144 |
| region.allele: intermediate random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: intermediate random: min P | 2.63e-03 (B*40:02) |
| region.allele: full fixed: markers / tests | 142 / 142 |
| region.allele: full fixed: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: full fixed: min P | 1.10e-03 (A*31:01) |
| region.allele: full random: markers / tests | 142 / 142 |
| region.allele: full random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: full random: min P | 2.31e-03 (A*31:01) |
| region.residue: narrow fixed: markers / tests | 672 / 494 |
| region.residue: narrow fixed: tests significant GEC / max-T / Bonferroni | 1 / 1 / 1 |
| region.residue: narrow fixed: min P | 3.54e-05 (DQB1:57:V) |
| region.residue: narrow random: markers / tests | 672 / 494 |
| region.residue: narrow random: tests significant GEC / max-T / Bonferroni | 1 / 1 / 1 |
| region.residue: narrow random: min P | 4.57e-05 (DQB1:57:V) |
| region.residue: intermediate fixed: markers / tests | 681 / 503 |
| region.residue: intermediate fixed: tests significant GEC / max-T / Bonferroni | 1 / 1 / 0 |
| region.residue: intermediate fixed: min P | 1.48e-04 (DQB1:57:V) |
| region.residue: intermediate random: markers / tests | 681 / 503 |
| region.residue: intermediate random: tests significant GEC / max-T / Bonferroni | 1 / 1 / 0 |
| region.residue: intermediate random: min P | 1.85e-04 (DQB1:57:V) |
| region.residue: full fixed: markers / tests | 679 / 512 |
| region.residue: full fixed: tests significant GEC / max-T / Bonferroni | 3 / 3 / 0 |
| region.residue: full fixed: min P | 1.08e-04 (DQB1:57:V) |
| region.residue: full random: markers / tests | 679 / 512 |
| region.residue: full random: tests significant GEC / max-T / Bonferroni | 1 / 1 / 0 |
| region.residue: full random: min P | 1.74e-04 (DQB1:57:V) |

## Full statistics

**Scan summary (one row per cohort × class × model)**

| class | cohort | model | samples | cases | controls | markers_tested | tests | me_gec | p_gec | p_maxt | p_bonf | tests_gec | tests_maxt | tests_bonf | lambda_gc | min_p | top_marker | top_or |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| allele | narrow | fixed | 2,168 | 419 | 1,749 | 143 | 143 | 114 | 4.388e-04 | 4.144e-04 | 3.497e-04 | 1 | 1 | 0 | 1.759 | 3.743e-04 | B*40:02 | 1.684 |
| allele | narrow | random | 2,195 | 419 | 1,776 | 143 | 143 | 114 | 4.388e-04 | 4.144e-04 | 3.497e-04 | 0 | 0 | 0 | 1.747 | 5.088e-04 | B*40:02 | 1.786 |
| allele | intermediate | fixed | 2,480 | 429 | 2,051 | 144 | 144 | 115 | 4.349e-04 | 4.145e-04 | 3.472e-04 | 0 | 0 | 0 | 1.906 | 0.001778 | B*40:02 | 1.534 |
| allele | intermediate | random | 2,508 | 429 | 2,079 | 144 | 144 | 115 | 4.349e-04 | 4.145e-04 | 3.472e-04 | 0 | 0 | 0 | 1.728 | 0.002633 | B*40:02 | 1.599 |
| allele | full | fixed | 3,071 | 439 | 2,632 | 142 | 142 | 113.6 | 4.402e-04 | 4.134e-04 | 3.521e-04 | 0 | 0 | 0 | 1.538 | 0.001098 | A*31:01 | 1.497 |
| allele | full | random | 3,101 | 439 | 2,662 | 142 | 142 | 113.6 | 4.402e-04 | 4.134e-04 | 3.521e-04 | 0 | 0 | 0 | 1.374 | 0.002306 | A*31:01 | 1.542 |
| residue | narrow | fixed | 2,168 | 419 | 1,749 | 672 | 494 | 201.8 | 2.477e-04 | 2.775e-04 | 1.012e-04 | 1 | 1 | 1 | 2.176 | 3.540e-05 | DQB1:57:V | 0.591 |
| residue | narrow | random | 2,195 | 419 | 1,776 | 672 | 494 | 201.8 | 2.477e-04 | 2.775e-04 | 1.012e-04 | 1 | 1 | 1 | 2.061 | 4.570e-05 | DQB1:57:V | 0.621 |
| residue | intermediate | fixed | 2,480 | 429 | 2,051 | 681 | 503 | 200.8 | 2.491e-04 | 2.730e-04 | 9.940e-05 | 1 | 1 | 0 | 2.067 | 1.482e-04 | DQB1:57:V | 0.626 |
| residue | intermediate | random | 2,508 | 429 | 2,079 | 681 | 503 | 200.8 | 2.491e-04 | 2.730e-04 | 9.940e-05 | 1 | 1 | 0 | 2.003 | 1.847e-04 | DQB1:57:V | 0.652 |
| residue | full | fixed | 3,071 | 439 | 2,632 | 679 | 512 | 200.9 | 2.489e-04 | 2.740e-04 | 9.766e-05 | 3 | 3 | 0 | 2.368 | 1.076e-04 | DQB1:57:V | 0.627 |
| residue | full | random | 3,101 | 439 | 2,662 | 679 | 512 | 200.9 | 2.489e-04 | 2.740e-04 | 9.766e-05 | 1 | 1 | 0 | 2.223 | 1.741e-04 | DQB1:57:V | 0.656 |

**Alleles: markers by cross-cohort consistency (random model)**

| threshold | 3/3 | 2/3 | 1/3 |
|---|---|---|---|
| LD-adjusted (GEC) (primary) | 0 | 0 | 0 |
| Exact (max-T) (validation) | 0 | 0 | 0 |
| Conventional Bonferroni (sensitivity) | 0 | 0 | 0 |

**Amino-acid residues: markers by cross-cohort consistency (random model)**

| threshold | 3/3 | 2/3 | 1/3 |
|---|---|---|---|
| LD-adjusted (GEC) (primary) | 1 | 0 | 0 |
| Exact (max-T) (validation) | 1 | 0 | 0 |
| Conventional Bonferroni (sensitivity) | 0 | 0 | 1 |

## How to read it

1. Find the gene band with the tallest points.
2. Read across a row: does the random model (right) reproduce the fixed model (left)?
3. Read down a column: does the peak hold as controls are added (narrow → full)?
4. Look up the exact P, OR and counts of any marker in _comparison/assoc_all.tsv.

## What this figure does *not* establish

- Which of several correlated markers carries the signal: see the conditional analyses (results/_conditional/).
- The three cohorts are nested and share almost all cases, so agreement across cohorts is a check of the control selection, not a replication.
- All cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X; an association that reflects a typing difference between platforms would look the same as a true one here.

---

Methods and rationale: [`METHODS.md`](../../../../docs/METHODS.md)
