# CTEPH with DVT vs controls: markers significant in at least one cohort

**Covers:**  in `04.hits/`. Not drawn, because no marker of the class is significant in any cohort: `hits.allele.png`, `hits.residue.png`

## The question this figure answers

For each marker significant at any of the three thresholds in at least one cohort, how large is the effect, and in how many of the three cohorts does it hold?

## Panels

## Interpretation

alleles: 0 markers significant in at least one cohort under some threshold, 0 at the LD-adjusted threshold in 3/3 cohorts; amino-acid residues: 0 markers significant in at least one cohort under some threshold, 0 at the LD-adjusted threshold in 3/3 cohorts. Consistency is read in the random model with the same direction of effect in every significant cohort; the cohorts are nested, so 3/3 is a sample-selection sensitivity analysis, not replication.

## Values in this rendering

| quantity | value |
|---|---|
| hits.allele: markers significant in >= 1 cohort | 0 |
| hits.allele: markers drawn | 0 |
| hits.residue: markers significant in >= 1 cohort | 0 |
| hits.residue: markers drawn | 0 |

## Full statistics

## How to read it

1. Read the consistency columns first: GEC (primary), then max-T and Bonferroni.
2. In the figure, check that the six estimates lie on the same side of 1.
3. In the table, a row with `mirror_of` set and its partner are one test (OR and 1/OR, same P).

## What this figure does *not* establish

- Independence of the markers drawn: many are on one haplotype; see the conditional analyses (results/_conditional/signals_all.tsv).
- Exact effect sizes in the random model: they are SAIGE's one-step score estimates (Tstat / var), an approximation (SAIGE's Firth-corrected estimate was not requested, as in the genome-wide scan); read them for direction and rough size.
- The three cohorts are nested and share almost all cases, so agreement across cohorts is a check of the control selection, not a replication.
- All cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X; an association that reflects a typing difference between platforms would look the same as a true one here.

---

Methods and rationale: [`METHODS.md`](../../../../docs/METHODS.md)
