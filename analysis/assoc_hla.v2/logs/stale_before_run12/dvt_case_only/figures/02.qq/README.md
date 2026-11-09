# CTEPH with DVT vs CTEPH without DVT (cases only): distribution of P values against the null

**Covers:** `qq.png` in `02.qq/`

## The question this figure answers

How far does each scan depart from the uniform null, and do the two models depart alike?

## Panels

**(qq.png a–f) rows alleles, residues; columns narrow, intermediate, full mainland**

Observed against expected −log10 P for the fixed (blue) and random (red) model over the markers each returned; expected quantiles are (i − 0.5)/n. The grey band is the pointwise 95 % interval of the uniform null, Beta(i, n − i + 1). The corner text gives each scan's genomic-control factor λGC = median χ²(1) / 0.455.

## Interpretation

λGC spans 0.84–1.62 over the 12 scans. All markers come from one region of chromosome 6 and are strongly correlated (mirror residue pairs, alleles and residues on shared haplotypes), so a modest number of real or artefactual differences between cases and controls lifts the whole curve; λGC here describes the region, not the calibration of the test. It was not applied to any P.

## Values in this rendering

| quantity | value |
|---|---|
| qq: allele narrow fixed: lambda_GC | 0.965 |
| qq: allele narrow random: lambda_GC | 1.142 |
| qq: allele intermediate fixed: lambda_GC | 0.836 |
| qq: allele intermediate random: lambda_GC | 1.153 |
| qq: allele full fixed: lambda_GC | 1.062 |
| qq: allele full random: lambda_GC | 1.12 |
| qq: residue narrow fixed: lambda_GC | 1.453 |
| qq: residue narrow random: lambda_GC | 1.438 |
| qq: residue intermediate fixed: lambda_GC | 1.599 |
| qq: residue intermediate random: lambda_GC | 1.618 |
| qq: residue full fixed: lambda_GC | 1.448 |
| qq: residue full random: lambda_GC | 1.529 |

## Full statistics

**λGC per scan**

| class | cohort | model | markers_tested | lambda_gc |
|---|---|---|---|---|
| allele | narrow | fixed | 143 | 0.965 |
| allele | narrow | random | 143 | 1.142 |
| allele | intermediate | fixed | 144 | 0.836 |
| allele | intermediate | random | 144 | 1.153 |
| allele | full | fixed | 142 | 1.062 |
| allele | full | random | 142 | 1.12 |
| residue | narrow | fixed | 672 | 1.453 |
| residue | narrow | random | 672 | 1.438 |
| residue | intermediate | fixed | 681 | 1.599 |
| residue | intermediate | random | 681 | 1.618 |
| residue | full | fixed | 679 | 1.448 |
| residue | full | random | 679 | 1.529 |

## How to read it

1. Compare the two curves in a panel: the models should lie on top of each other.
2. Compare across a row: the departure as controls are added.
3. Do not read λGC as genome-wide inflation (see Interpretation).

## What this figure does *not* establish

- Whether the departure is biology or artefact: the genome-wide λ of the same models is the calibration reference (analysis/assoc_plink2, assoc_saige).
- The three cohorts share almost all their cases and this set has no controls, so the three cohorts are nearly the same analysis; agreement across them is expected, not a replication.
- Both groups are cases, so the case–control platform split does not enter. But the sequencing platform also differs between cases with and without DVT (docs/STUDY_NOTES.md), so a platform-dependent typing difference would still look like an association; platform is not a covariate here.

---

Methods and rationale: [`METHODS.md`](../../../../docs/METHODS.md)
