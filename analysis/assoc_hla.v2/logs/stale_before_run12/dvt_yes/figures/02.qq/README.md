# CTEPH with DVT vs controls: distribution of P values against the null

**Covers:** `qq.png` in `02.qq/`

## The question this figure answers

How far does each scan depart from the uniform null, and do the two models depart alike?

## Panels

**(qq.png a–f) rows alleles, residues; columns narrow, intermediate, full mainland**

Observed against expected −log10 P for the fixed (blue) and random (red) model over the markers each returned; expected quantiles are (i − 0.5)/n. The grey band is the pointwise 95 % interval of the uniform null, Beta(i, n − i + 1). The corner text gives each scan's genomic-control factor λGC = median χ²(1) / 0.455.

## Interpretation

λGC spans 1.28–2.15 over the 12 scans. All markers come from one region of chromosome 6 and are strongly correlated (mirror residue pairs, alleles and residues on shared haplotypes), so a modest number of real or artefactual differences between cases and controls lifts the whole curve; λGC here describes the region, not the calibration of the test. It was not applied to any P.

## Values in this rendering

| quantity | value |
|---|---|
| qq: allele narrow fixed: lambda_GC | 1.281 |
| qq: allele narrow random: lambda_GC | 1.464 |
| qq: allele intermediate fixed: lambda_GC | 1.455 |
| qq: allele intermediate random: lambda_GC | 1.699 |
| qq: allele full fixed: lambda_GC | 1.34 |
| qq: allele full random: lambda_GC | 1.414 |
| qq: residue narrow fixed: lambda_GC | 2.148 |
| qq: residue narrow random: lambda_GC | 2.012 |
| qq: residue intermediate fixed: lambda_GC | 1.928 |
| qq: residue intermediate random: lambda_GC | 1.891 |
| qq: residue full fixed: lambda_GC | 1.727 |
| qq: residue full random: lambda_GC | 1.681 |

## Full statistics

**λGC per scan**

| class | cohort | model | markers_tested | lambda_gc |
|---|---|---|---|---|
| allele | narrow | fixed | 143 | 1.281 |
| allele | narrow | random | 143 | 1.464 |
| allele | intermediate | fixed | 144 | 1.455 |
| allele | intermediate | random | 144 | 1.699 |
| allele | full | fixed | 142 | 1.34 |
| allele | full | random | 142 | 1.414 |
| residue | narrow | fixed | 672 | 2.148 |
| residue | narrow | random | 672 | 2.012 |
| residue | intermediate | fixed | 681 | 1.928 |
| residue | intermediate | random | 681 | 1.891 |
| residue | full | fixed | 679 | 1.727 |
| residue | full | random | 679 | 1.681 |

## How to read it

1. Compare the two curves in a panel: the models should lie on top of each other.
2. Compare across a row: the departure as controls are added.
3. Do not read λGC as genome-wide inflation (see Interpretation).

## What this figure does *not* establish

- Whether the departure is biology or artefact: the genome-wide λ of the same models is the calibration reference (analysis/assoc_plink2, assoc_saige).
- The three cohorts are nested and share almost all cases, so agreement across cohorts is a check of the control selection, not a replication.
- All cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X; an association that reflects a typing difference between platforms would look the same as a true one here.

---

Methods and rationale: [`METHODS.md`](../../../../docs/METHODS.md)
