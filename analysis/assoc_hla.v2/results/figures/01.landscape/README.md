# Where in the MHC, and in how many cohorts

**Main figure.** **Covers:** `bars.all.png`, `bars.dvt_case_only.png`, `bars.dvt_no.png`, `bars.dvt_yes.png`, `manhattan.all.png`, `manhattan.dvt_case_only.png`, `manhattan.dvt_no.png`, `manhattan.dvt_yes.png` in `01.landscape/`

## The question this figure answers

Along the typed HLA genes, which allele and residue tests are associated in each case set, and in how many of the three nested cohorts does each hold?

## Panels

**(bars.<case_set>.png a-b) one bar per test, alleles (a) and residues (b)**

For each test (a residue mirror pair is drawn once) the bar runs from −log10 of its largest P to −log10 of its smallest P over the three cohorts, random model: its top is the strongest cohort and its bottom the weakest. Colour and line width give the number of cohorts in which the test is significant at that cohort's LD-adjusted (GEC) threshold with the same direction (3/3 darkest, then 2/3, 1/3 "Single cohort"; grey not significant). Horizontal lines: GEC (solid), exact max-T (dashed) and conventional Bonferroni (dotted), each the median of the three cohorts' thresholds (they differ by < 1 %). The violet triangle marks where the genome-wide MHC lead lies among the genes. Signal leads are labelled with their n/3.

**(manhattan.<case_set>.png a-f) every test in each cohort**

Rows narrow, intermediate, full mainland; columns alleles, residues. One point per test, red if significant at that cohort's own GEC threshold; the three lines are the cohort's own thresholds. The same data as the bars, cohort by cohort.

## Interpretation

Random model. Independent signals (02.signals, 03.conditioning): DQB1 57 V (all cases, residue: GEC 3/3, max-T 3/3, Bonferroni 1/3); DQB1 87 Y (DVT no, residue: GEC 1/3, max-T 1/3, Bonferroni 0/3); DPB1*05:01 (DVT no, allele: GEC 1/3, max-T 1/3, Bonferroni 1/3). In the bar encoding, a test whose whole bar lies above the solid line is significant in all three cohorts; a bar that crosses it with its top only is significant in one. Residues on one haplotype rise together, so a run of coloured bars can be one signal: in DVT no, 51 residue markers are significant in some cohort, and the conditional analysis (03.conditioning) asks how many signals they are.

## Values in this rendering

| quantity | value |
|---|---|
| bars.all: allele: tests drawn | 145 |
| bars.all: allele: tests 3/3 at GEC | 0 |
| bars.all: allele: tests 2/3 at GEC | 0 |
| bars.all: allele: tests 1/3 at GEC | 0 |
| bars.all: allele: thresholds GEC / max-T / Bonferroni (median over cohorts) | 4.39e-04 / 4.14e-04 / 3.50e-04 |
| bars.all: residue: tests drawn | 513 |
| bars.all: residue: tests 3/3 at GEC | 1 |
| bars.all: residue: tests 2/3 at GEC | 0 |
| bars.all: residue: tests 1/3 at GEC | 0 |
| bars.all: residue: thresholds GEC / max-T / Bonferroni (median over cohorts) | 2.49e-04 / 2.74e-04 / 9.94e-05 |
| bars.dvt_case_only: allele: tests drawn | 145 |
| bars.dvt_case_only: allele: tests 3/3 at GEC | 0 |
| bars.dvt_case_only: allele: tests 2/3 at GEC | 0 |
| bars.dvt_case_only: allele: tests 1/3 at GEC | 0 |
| bars.dvt_case_only: allele: thresholds GEC / max-T / Bonferroni (median over cohorts) | 4.49e-04 / 4.31e-04 / 3.50e-04 |
| bars.dvt_case_only: residue: tests drawn | 463 |
| bars.dvt_case_only: residue: tests 3/3 at GEC | 0 |
| bars.dvt_case_only: residue: tests 2/3 at GEC | 0 |
| bars.dvt_case_only: residue: tests 1/3 at GEC | 0 |
| bars.dvt_case_only: residue: thresholds GEC / max-T / Bonferroni (median over cohorts) | 2.68e-04 / 2.93e-04 / 1.08e-04 |
| bars.dvt_no: allele: tests drawn | 145 |
| bars.dvt_no: allele: tests 3/3 at GEC | 0 |
| bars.dvt_no: allele: tests 2/3 at GEC | 0 |
| bars.dvt_no: allele: tests 1/3 at GEC | 1 |
| bars.dvt_no: allele: thresholds GEC / max-T / Bonferroni (median over cohorts) | 4.38e-04 / 4.14e-04 / 3.50e-04 |
| bars.dvt_no: residue: tests drawn | 511 |
| bars.dvt_no: residue: tests 3/3 at GEC | 0 |
| bars.dvt_no: residue: tests 2/3 at GEC | 0 |
| bars.dvt_no: residue: tests 1/3 at GEC | 34 |
| bars.dvt_no: residue: thresholds GEC / max-T / Bonferroni (median over cohorts) | 2.50e-04 / 2.77e-04 / 1.00e-04 |
| bars.dvt_yes: allele: tests drawn | 145 |
| bars.dvt_yes: allele: tests 3/3 at GEC | 0 |
| bars.dvt_yes: allele: tests 2/3 at GEC | 0 |
| bars.dvt_yes: allele: tests 1/3 at GEC | 0 |
| bars.dvt_yes: allele: thresholds GEC / max-T / Bonferroni (median over cohorts) | 4.37e-04 / 4.15e-04 / 3.50e-04 |
| bars.dvt_yes: residue: tests drawn | 512 |
| bars.dvt_yes: residue: tests 3/3 at GEC | 0 |
| bars.dvt_yes: residue: tests 2/3 at GEC | 0 |
| bars.dvt_yes: residue: tests 1/3 at GEC | 0 |
| bars.dvt_yes: residue: thresholds GEC / max-T / Bonferroni (median over cohorts) | 2.49e-04 / 2.76e-04 / 9.98e-05 |
| manhattan.all: allele narrow: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| manhattan.all: allele intermediate: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| manhattan.all: allele full: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| manhattan.all: residue narrow: tests significant GEC / max-T / Bonferroni | 1 / 1 / 1 |
| manhattan.all: residue intermediate: tests significant GEC / max-T / Bonferroni | 1 / 1 / 0 |
| manhattan.all: residue full: tests significant GEC / max-T / Bonferroni | 1 / 1 / 0 |
| manhattan.dvt_case_only: allele narrow: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| manhattan.dvt_case_only: allele intermediate: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| manhattan.dvt_case_only: allele full: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| manhattan.dvt_case_only: residue narrow: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| manhattan.dvt_case_only: residue intermediate: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| manhattan.dvt_case_only: residue full: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| manhattan.dvt_no: allele narrow: tests significant GEC / max-T / Bonferroni | 1 / 1 / 1 |
| manhattan.dvt_no: allele intermediate: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| manhattan.dvt_no: allele full: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| manhattan.dvt_no: residue narrow: tests significant GEC / max-T / Bonferroni | 34 / 34 / 0 |
| manhattan.dvt_no: residue intermediate: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| manhattan.dvt_no: residue full: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| manhattan.dvt_yes: allele narrow: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| manhattan.dvt_yes: allele intermediate: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| manhattan.dvt_yes: allele full: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| manhattan.dvt_yes: residue narrow: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| manhattan.dvt_yes: residue intermediate: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |
| manhattan.dvt_yes: residue full: tests significant GEC / max-T / Bonferroni | 0 / 0 / 0 |

## Full statistics

**Tests by cross-cohort consistency (random model), per case set and class**

| case set | class | tests drawn | gec: 3/3 / 2/3 / 1/3 | maxt: 3/3 / 2/3 / 1/3 | bonf: 3/3 / 2/3 / 1/3 |
|---|---|---|---|---|---|
| all cases | allele | 145 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |
| all cases | residue | 513 | 1 / 0 / 0 | 1 / 0 / 0 | 0 / 0 / 1 |
| DVT yes | allele | 145 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |
| DVT yes | residue | 512 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |
| DVT no | allele | 145 | 0 / 0 / 1 | 0 / 0 / 1 | 0 / 0 / 1 |
| DVT no | residue | 511 | 0 / 0 / 34 | 0 / 0 / 34 | 0 / 0 / 0 |
| DVT yes vs no, cases only | allele | 145 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |
| DVT yes vs no, cases only | residue | 463 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |

## How to read it

1. Start with bars.all.png: find the bars above the solid line and their colour (n/3).
2. Check the dashed and dotted lines: does the call survive the exact and the conventional threshold?
3. Open manhattan.<case_set>.png of the same case set to see which cohort carries a 1/3 call.
4. Compare case sets (bars.dvt_yes / dvt_no / dvt_case_only) with 02.signals for the effect sizes.

## What this figure does *not* establish

- Independence of coloured bars: neighbouring residues and alleles share haplotypes; see 03.conditioning.
- Physical position: the x axis orders genes and markers; it is not a coordinate.
- Replication: the three cohorts are nested and share almost all cases; n/3 is a sample-selection sensitivity analysis, not independent replication.
- Platform: all cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X; a typing difference between platforms would look like an association, and no covariate separates the two. Within cases the platform mix also differs by DVT status.

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
