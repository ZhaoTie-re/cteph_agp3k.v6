# CTEPH with DVT vs controls: association along the HLA genes

**Covers:** `region.allele.png`, `region.residue.png` in `01.region/`

## The question this figure answers

Where along the typed HLA genes do the allele and the amino-acid residue markers differ between the two groups compared (CTEPH with DVT vs controls), and does the picture hold across the three cohorts and both models?

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
| region.allele: narrow fixed: min P | 1.57e-03 (C*03:04) |
| region.allele: narrow random: markers / tests | 143 / 143 |
| region.allele: narrow random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: narrow random: min P | 1.57e-03 (C*03:04) |
| region.allele: intermediate fixed: markers / tests | 144 / 144 |
| region.allele: intermediate fixed: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: intermediate fixed: min P | 3.11e-03 (C*03:04) |
| region.allele: intermediate random: markers / tests | 144 / 144 |
| region.allele: intermediate random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: intermediate random: min P | 3.03e-03 (C*03:04) |
| region.allele: full fixed: markers / tests | 142 / 142 |
| region.allele: full fixed: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: full fixed: min P | 4.31e-03 (DQB1*05:02) |
| region.allele: full random: markers / tests | 142 / 142 |
| region.allele: full random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: full random: min P | 4.39e-03 (DQB1*05:02) |
| region.residue: narrow fixed: markers / tests | 672 / 492 |
| region.residue: narrow fixed: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.residue: narrow fixed: min P | 3.21e-03 (B:114:D) |
| region.residue: narrow random: markers / tests | 672 / 492 |
| region.residue: narrow random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.residue: narrow random: min P | 3.27e-03 (C:95:L) |
| region.residue: intermediate fixed: markers / tests | 681 / 501 |
| region.residue: intermediate fixed: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.residue: intermediate fixed: min P | 2.64e-03 (B:114:D) |
| region.residue: intermediate random: markers / tests | 681 / 501 |
| region.residue: intermediate random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.residue: intermediate random: min P | 3.55e-03 (B:114:D) |
| region.residue: full fixed: markers / tests | 679 / 511 |
| region.residue: full fixed: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.residue: full fixed: min P | 1.28e-03 (B:114:D) |
| region.residue: full random: markers / tests | 679 / 511 |
| region.residue: full random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.residue: full random: min P | 1.72e-03 (B:114:D) |

## Full statistics

**Scan summary (one row per cohort × class × model)**

| class | cohort | model | samples | cases | controls | markers_tested | tests | me_gec | p_gec | p_maxt | p_bonf | tests_gec | tests_maxt | tests_bonf | lambda_gc | min_p | top_marker | top_or |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| allele | narrow | fixed | 1,899 | 150 | 1,749 | 143 | 143 | 114.3 | 4.374e-04 | 4.166e-04 | 3.497e-04 | 0 | 0 | 0 | 1.281 | 0.001568 | C*03:04 | 1.738 |
| allele | narrow | random | 1,926 | 150 | 1,776 | 143 | 143 | 114.3 | 4.374e-04 | 4.166e-04 | 3.497e-04 | 0 | 0 | 0 | 1.464 | 0.001569 | C*03:04 | 1.913 |
| allele | intermediate | fixed | 2,206 | 155 | 2,051 | 144 | 144 | 115.3 | 4.337e-04 | 4.063e-04 | 3.472e-04 | 0 | 0 | 0 | 1.455 | 0.003107 | C*03:04 | 1.636 |
| allele | intermediate | random | 2,234 | 155 | 2,079 | 144 | 144 | 115.3 | 4.337e-04 | 4.063e-04 | 3.472e-04 | 0 | 0 | 0 | 1.699 | 0.003032 | C*03:04 | 1.78 |
| allele | full | fixed | 2,791 | 159 | 2,632 | 142 | 142 | 113.8 | 4.394e-04 | 4.146e-04 | 3.521e-04 | 0 | 0 | 0 | 1.34 | 0.004307 | DQB1*05:02 | 2.231 |
| allele | full | random | 2,821 | 159 | 2,662 | 142 | 142 | 113.8 | 4.394e-04 | 4.146e-04 | 3.521e-04 | 0 | 0 | 0 | 1.414 | 0.004392 | DQB1*05:02 | 3.003 |
| residue | narrow | fixed | 1,899 | 150 | 1,749 | 672 | 492 | 200.2 | 2.497e-04 | 2.772e-04 | 1.016e-04 | 0 | 0 | 0 | 2.148 | 0.003206 | B:114:D | 0.673 |
| residue | narrow | random | 1,926 | 150 | 1,776 | 672 | 492 | 200.2 | 2.497e-04 | 2.772e-04 | 1.016e-04 | 0 | 0 | 0 | 2.012 | 0.003266 | C:95:L | 0.667 |
| residue | intermediate | fixed | 2,206 | 155 | 2,051 | 681 | 501 | 200.5 | 2.494e-04 | 2.744e-04 | 9.980e-05 | 0 | 0 | 0 | 1.928 | 0.002638 | B:114:D | 0.675 |
| residue | intermediate | random | 2,234 | 155 | 2,079 | 681 | 501 | 200.5 | 2.494e-04 | 2.744e-04 | 9.980e-05 | 0 | 0 | 0 | 1.891 | 0.003551 | B:114:D | 0.693 |
| residue | full | fixed | 2,791 | 159 | 2,632 | 679 | 511 | 200.6 | 2.493e-04 | 2.758e-04 | 9.785e-05 | 0 | 0 | 0 | 1.727 | 0.001279 | B:114:D | 0.661 |
| residue | full | random | 2,821 | 159 | 2,662 | 679 | 511 | 200.6 | 2.493e-04 | 2.758e-04 | 9.785e-05 | 0 | 0 | 0 | 1.681 | 0.001722 | B:114:D | 0.678 |

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
- The three cohorts are nested and share almost all cases, so agreement across cohorts is a check of the control selection, not a replication.
- All cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X; an association that reflects a typing difference between platforms would look the same as a true one here.

---

Methods and rationale: [`METHODS.md`](../../../../docs/METHODS.md)
