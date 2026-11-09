# CTEPH without DVT vs controls: association along the HLA genes

**Covers:** `region.allele.png`, `region.residue.png` in `01.region/`

## The question this figure answers

Where along the typed HLA genes do the allele and the amino-acid residue markers differ between the two groups compared (CTEPH without DVT vs controls), and does the picture hold across the three cohorts and both models?

## Panels

**(region.allele.png) rows narrow, intermediate, full mainland; columns fixed model, random model**

Each point is one tested allele, at its place in one marker order shared by all panels: genes in GRCh38 chromosomal order (alternate genes shaded), then residues by IMGT position and residue, alleles by their numeric fields. Height is −log10 P (fixed: plink2 Wald P; random: SAIGE unadjusted score P). Three family-wise thresholds (0.05 over the class's tests in that cohort): solid red = LD-adjusted, 0.05 / Me with Me from GEC v0.2 (primary); dashed grey = exact max-T (validation); dotted = conventional Bonferroni, 0.05 / tests (sensitivity). Red points are significant at the GEC line; the three smallest of them per panel are labelled.

**(region.residue.png) rows narrow, intermediate, full mainland; columns fixed model, random model**

Each point is one tested amino-acid residue, at its place in one marker order shared by all panels: genes in GRCh38 chromosomal order (alternate genes shaded), then residues by IMGT position and residue, alleles by their numeric fields. Height is −log10 P (fixed: plink2 Wald P; random: SAIGE unadjusted score P). Three family-wise thresholds (0.05 over the class's tests in that cohort): solid red = LD-adjusted, 0.05 / Me with Me from GEC v0.2 (primary); dashed grey = exact max-T (validation); dotted = conventional Bonferroni, 0.05 / tests (sensitivity). Red points are significant at the GEC line; the three smallest of them per panel are labelled.

## Interpretation

At the LD-adjusted threshold, at most 1 allele and 33 residue tests are significant in any one scan. A residue position with two residues gives two mirror columns with the same P (one test), and neighbouring positions on one haplotype carry nearly the same signal; the number of points above a line is not the number of independent signals (see _conditional/). Cross-cohort consistency is read in the random model: the tables below count, per threshold, the markers significant in 3/3, 2/3 and 1/3 ("Single cohort") cohorts with the same direction of effect.

## Values in this rendering

| quantity | value |
|---|---|
| region.allele: narrow fixed: markers / tests | 143 / 143 |
| region.allele: narrow fixed: tests significant GEC / max-T / Bonferroni | 1 / 1 / 1 |
| region.allele: narrow fixed: min P | 1.91e-04 (DPB1*05:01) |
| region.allele: narrow random: markers / tests | 143 / 143 |
| region.allele: narrow random: tests significant GEC / max-T / Bonferroni | 1 / 1 / 1 |
| region.allele: narrow random: min P | 2.81e-04 (DPB1*05:01) |
| region.allele: intermediate fixed: markers / tests | 144 / 144 |
| region.allele: intermediate fixed: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: intermediate fixed: min P | 6.02e-04 (DPB1*05:01) |
| region.allele: intermediate random: markers / tests | 144 / 144 |
| region.allele: intermediate random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: intermediate random: min P | 8.44e-04 (DPB1*05:01) |
| region.allele: full fixed: markers / tests | 142 / 142 |
| region.allele: full fixed: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: full fixed: min P | 1.03e-03 (DPB1*05:01) |
| region.allele: full random: markers / tests | 142 / 142 |
| region.allele: full random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.allele: full random: min P | 1.41e-03 (DPB1*05:01) |
| region.residue: narrow fixed: markers / tests | 672 / 491 |
| region.residue: narrow fixed: tests significant GEC / max-T / Bonferroni | 30 / 31 / 0 |
| region.residue: narrow fixed: min P | 1.15e-04 (DQA1:11:C) |
| region.residue: narrow random: markers / tests | 672 / 491 |
| region.residue: narrow random: tests significant GEC / max-T / Bonferroni | 33 / 33 / 0 |
| region.residue: narrow random: min P | 1.19e-04 (DQB1:87:Y) |
| region.residue: intermediate fixed: markers / tests | 681 / 500 |
| region.residue: intermediate fixed: tests significant GEC / max-T / Bonferroni | 0 / 1 / 0 |
| region.residue: intermediate fixed: min P | 2.69e-04 (DQA1:218:Q) |
| region.residue: intermediate random: markers / tests | 681 / 500 |
| region.residue: intermediate random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.residue: intermediate random: min P | 3.29e-04 (DQB1:87:Y) |
| region.residue: full fixed: markers / tests | 679 / 510 |
| region.residue: full fixed: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.residue: full fixed: min P | 3.10e-04 (DQB1:87:Y) |
| region.residue: full random: markers / tests | 679 / 510 |
| region.residue: full random: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| region.residue: full random: min P | 3.15e-04 (DQB1:87:Y) |

## Full statistics

**Scan summary (one row per cohort × class × model)**

| class | cohort | model | samples | cases | controls | markers_tested | tests | me_gec | p_gec | p_maxt | p_bonf | tests_gec | tests_maxt | tests_bonf | lambda_gc | min_p | top_marker | top_or |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| allele | narrow | fixed | 1,889 | 140 | 1,749 | 143 | 143 | 114 | 4.385e-04 | 4.179e-04 | 3.497e-04 | 1 | 1 | 1 | 1.725 | 1.908e-04 | DPB1*05:01 | 1.61 |
| allele | narrow | random | 1,916 | 140 | 1,776 | 143 | 143 | 114 | 4.385e-04 | 4.179e-04 | 3.497e-04 | 1 | 1 | 1 | 1.662 | 2.815e-04 | DPB1*05:01 | 1.629 |
| allele | intermediate | fixed | 2,192 | 141 | 2,051 | 144 | 144 | 115 | 4.346e-04 | 4.129e-04 | 3.472e-04 | 0 | 0 | 0 | 1.665 | 6.024e-04 | DPB1*05:01 | 1.538 |
| allele | intermediate | random | 2,220 | 141 | 2,079 | 144 | 144 | 115 | 4.346e-04 | 4.129e-04 | 3.472e-04 | 0 | 0 | 0 | 1.564 | 8.435e-04 | DPB1*05:01 | 1.554 |
| allele | full | fixed | 2,777 | 145 | 2,632 | 142 | 142 | 113.7 | 4.399e-04 | 4.139e-04 | 3.521e-04 | 0 | 0 | 0 | 1.605 | 0.001031 | DPB1*05:01 | 1.496 |
| allele | full | random | 2,807 | 145 | 2,662 | 142 | 142 | 113.7 | 4.399e-04 | 4.139e-04 | 3.521e-04 | 0 | 0 | 0 | 1.549 | 0.001415 | DPB1*05:01 | 1.507 |
| residue | narrow | fixed | 1,889 | 140 | 1,749 | 672 | 491 | 200.1 | 2.499e-04 | 2.770e-04 | 1.018e-04 | 30 | 31 | 0 | 3.279 | 1.152e-04 | DQA1:11:C | 0.596 |
| residue | narrow | random | 1,916 | 140 | 1,776 | 672 | 491 | 200.1 | 2.499e-04 | 2.770e-04 | 1.018e-04 | 33 | 33 | 0 | 3.047 | 1.190e-04 | DQB1:87:Y | 0.543 |
| residue | intermediate | fixed | 2,192 | 141 | 2,051 | 681 | 500 | 200.3 | 2.497e-04 | 2.731e-04 | 1.000e-04 | 0 | 1 | 0 | 2.773 | 2.691e-04 | DQA1:218:Q | 0.617 |
| residue | intermediate | random | 2,220 | 141 | 2,079 | 681 | 500 | 200.3 | 2.497e-04 | 2.731e-04 | 1.000e-04 | 0 | 0 | 0 | 2.584 | 3.290e-04 | DQB1:87:Y | 0.569 |
| residue | full | fixed | 2,777 | 145 | 2,632 | 679 | 510 | 200.4 | 2.495e-04 | 2.767e-04 | 9.804e-05 | 0 | 0 | 0 | 2.497 | 3.104e-04 | DQB1:87:Y | 0.518 |
| residue | full | random | 2,807 | 145 | 2,662 | 679 | 510 | 200.4 | 2.495e-04 | 2.767e-04 | 9.804e-05 | 0 | 0 | 0 | 2.192 | 3.149e-04 | DQB1:87:Y | 0.577 |

**Alleles: markers by cross-cohort consistency (random model)**

| threshold | 3/3 | 2/3 | 1/3 |
|---|---|---|---|
| LD-adjusted (GEC) (primary) | 0 | 0 | 1 |
| Exact (max-T) (validation) | 0 | 0 | 1 |
| Conventional Bonferroni (sensitivity) | 0 | 0 | 1 |

**Amino-acid residues: markers by cross-cohort consistency (random model)**

| threshold | 3/3 | 2/3 | 1/3 |
|---|---|---|---|
| LD-adjusted (GEC) (primary) | 0 | 0 | 51 |
| Exact (max-T) (validation) | 0 | 0 | 51 |
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
