# CTEPH with DVT vs CTEPH without DVT (cases only): markers significant in at least one cohort

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
- The three cohorts share almost all their cases and this set has no controls, so the three cohorts are nearly the same analysis; agreement across them is expected, not a replication.
- Both groups are cases, so the case–control platform split does not enter. But the sequencing platform also differs between cases with and without DVT (docs/STUDY_NOTES.md), so a platform-dependent typing difference would still look like an association; platform is not a covariate here.

---

Methods and rationale: [`METHODS.md`](../../../../docs/METHODS.md)
