# MHC association scan — narrow_mainland, HLA SEX + 10 bbj_mainland PCs, full GRM model

**Figure file:** `scan.random.png`

## The question this figure answers

Where in the classical MHC does this cohort show association under this model, in which GENE does it sit, and is any of it the shape a typing artefact makes rather than a result?

## Panels

**(a) HLA allele markers**

One point per tested HLA allele marker, x its marker position on chr6, y = -log10 P from the association engine. Points are coloured by gene in alternating tones, in chr6 order, so adjacent genes separate without a legend and each point matches its box in (c). The solid red line is this class's own Bonferroni threshold (3e-4 = 0.05 over the allele markers tested), which is the decision rule for this panel. An open red ring marks a significant allele with in_reference == 0 in the marker QC — a two-field call the reference allele panel never carries.

**(b) Amino-acid residue markers**

The same construction over the residue markers, on the SAME y limit as (a) so a taller point is a stronger point in either track, and against its OWN Bonferroni line (7e-5). The two classes are separate analyses with separate multiple-testing burdens — that is why they are separate panels rather than one scatter. Residue markers carry no in_reference flag: the reference panel is an allele panel, so "not in it" is not a statement about a residue, and none is ringed.

**(c) Gene layout**

One box per HLA gene, spanning the chr6 extent of that gene's markers, in the same two-tone colour the points above use. A name sits directly under its own box unless a neighbour is closer than one label width, in which case that run of names — and only that run — is spread at exactly one label width and each displaced name keeps an L-shaped leader back to its box. The boxes are centred on their true extent and drawn at least 12 kb wide: a class I gene is 3-5 kb, under a thousandth of this axis, and would otherwise be invisible.

## Interpretation

1 allele and 1 residue marker(s) clear their own Bonferroni threshold. The two counts are NOT additive evidence: an allele marker and the residue markers it carries are the same chromosomes counted twice, so a peak appearing in both tracks at the same gene is one signal seen two ways, not two. Reading the position off the x-axis alone is what (c) exists to prevent — the marker coordinate is gene start plus index, so it identifies a GENE and no more. 1 significant allele marker(s) are absent from the reference allele panel and are ringed. hla.typing found that the controls carry roughly twice the share of such calls as the cases, differentially by phenotype, so a hit on one of them has the shape a typing artefact makes and cannot be reported as a finding without that check.

## Values in this rendering

| quantity | value |
|---|---|
| cohort | narrow_mainland |
| model | random |
| allele markers drawn | 152 |
| residue markers drawn | 706 |
| allele markers excluded | 0 |
| residue markers excluded | 0 |
| allele markers passing QC | 152 |
| residue markers passing QC | 706 |
| genes drawn | 16 |
| alpha Bonferroni allele | 3.289e-04 |
| alpha Bonferroni residue | 7.082e-05 |
| alpha effective | 4.587e-04 |
| alpha genome-wide | 5.000e-08 |
| allele markers significant | 1 |
| residue markers significant | 1 |
| significant alleles absent from the reference panel | 1 |
| markers named in (a) | 1 |
| markers named in (b) | 1 |
| strongest marker | AA_DQB1_57_V |
| smallest P | 4.866e-05 |
| axis (Mb) | 29.67–33.13 |

## Full statistics

**Strongest allele markers**

| ID | gene | position | A1_FREQ | OBS_CT | OR | L95 | U95 | P | in_reference |
|---|---|---|---|---|---|---|---|---|---|
| HLA_F*01 | F |  | 0.006378 | 2,195 | 2.886 | 1.585 | 5.253 | 2.625e-04 | 0 |
| HLA_B*40:02 | B |  | 0.05809 | 2,195 | 1.721 | 1.255 | 2.359 | 3.705e-04 | 1 |
| HLA_DQA1*01:01 | DQA1 |  | 0.06697 | 2,195 | 0.5561 | 0.3728 | 0.8296 | 0.002018 | 1 |
| HLA_A*31:01 | A |  | 0.07312 | 2,195 | 1.526 | 1.131 | 2.06 | 0.002837 | 1 |
| HLA_DQB1*05:01 | DQB1 |  | 0.06834 | 2,195 | 0.58 | 0.3898 | 0.8629 | 0.003599 | 1 |
| HLA_DPB1*05:01 | DPB1 |  | 0.3478 | 2,195 | 1.274 | 1.067 | 1.522 | 0.00375 | 1 |
| HLA_DRB1*01:01 | DRB1 |  | 0.06697 | 2,194 | 0.5947 | 0.3993 | 0.8858 | 0.005281 | 1 |
| HLA_C*03:04 | C |  | 0.0959 | 2,194 | 1.416 | 1.055 | 1.9 | 0.01018 | 1 |
| HLA_DRB1*13:02 | DRB1 |  | 0.07813 | 2,194 | 0.6829 | 0.4806 | 0.9703 | 0.01667 | 1 |
| HLA_A*33:03 | A |  | 0.08679 | 2,195 | 0.6928 | 0.4916 | 0.9765 | 0.01805 | 1 |
| HLA_A*02:07 | A |  | 0.0344 | 2,195 | 0.594 | 0.3573 | 0.9877 | 0.02234 | 1 |
| HLA_B*51:01 | B |  | 0.07221 | 2,195 | 1.405 | 0.9826 | 2.008 | 0.03117 | 1 |

**Strongest residue markers**

| ID | gene | position | A1_FREQ | OBS_CT | OR | L95 | U95 | P | in_reference |
|---|---|---|---|---|---|---|---|---|---|
| AA_DQB1_57_V | DQB1 | DQB1:57 | 0.1556 | 2,195 | 0.5899 | 0.4524 | 0.7692 | 4.866e-05 | — |
| AA_DQA1_175_E | DQA1 | DQA1:175 | 0.4371 | 2,195 | 1.324 | 1.115 | 1.571 | 6.617e-04 | — |
| AA_B_114_D | B | B:114 | 0.3888 | 2,195 | 0.7487 | 0.6258 | 0.8958 | 7.806e-04 | — |
| AA_B_114_N | B | B:114 | 0.6087 | 2,195 | 1.331 | 1.113 | 1.592 | 8.760e-04 | — |
| AA_DQB1_87_Y | DQB1 | DQB1:87 | 0.2153 | 2,195 | 0.708 | 0.5671 | 0.8839 | 0.001143 | — |
| AA_DRB1_11_L | DRB1 | DRB1:11 | 0.06902 | 2,194 | 0.5723 | 0.3872 | 0.8459 | 0.002561 | — |
| AA_DRB1_30_C | DRB1 | DRB1:30 | 0.06879 | 2,194 | 0.5754 | 0.3891 | 0.8509 | 0.002813 | — |
| AA_DQB1_53_L | DQB1 | DQB1:53 | 0.5089 | 2,195 | 1.284 | 1.073 | 1.536 | 0.003228 | — |
| AA_DQB1_53_Q | DQB1 | DQB1:53 | 0.4911 | 2,195 | 0.7791 | 0.651 | 0.9324 | 0.003228 | — |
| AA_DQB1_85_L | DQB1 | DQB1:85 | 0.5091 | 2,195 | 1.283 | 1.072 | 1.535 | 0.003331 | — |
| AA_DQB1_85_V | DQB1 | DQB1:85 | 0.4909 | 2,195 | 0.7796 | 0.6513 | 0.9332 | 0.003331 | — |
| AA_DQB1_86_E | DQB1 | DQB1:86 | 0.5091 | 2,195 | 1.283 | 1.072 | 1.535 | 0.003331 | — |

## How to read it

1. Read (a) and (b) against their OWN red line. The two classes carry different numbers of tests, so the same P is significant in one panel and not in the other; that is the design, not an inconsistency.
2. Take any peak and read straight down into (c). The gene box under it names the gene the signal is in. Nothing finer than the gene is claimed by this figure.
3. Compare heights ACROSS the two panels freely — they share one y limit and, after the label strips are taken out, one axes height.
4. Check whether a significant allele carries a red ring before quoting it. A ringed marker needs the typing check in hla.typing before it is a result.
5. The dotted genome-wide line is a reference to the rest of the study, not this figure's decision rule. A marker between the Bonferroni line and it has cleared the burden that applies to it.

## What this figure does *not* establish

- It does not fine-map to a base. Marker positions are gene start plus index within the gene, so the axis resolves genes, not variants.
- It cannot separate one signal from another inside the MHC. Long-range LD across the region means a peak at one gene may be tagging a causal allele at a neighbouring one; that is what the conditional rounds and the omnibus figure are for.
- An allele absent from the reference panel is not thereby a wrong call — a small panel cannot sample a rare allele. The ring is a flag for follow-up, never a verdict.
- The two tracks are not independent evidence about each other, so no combined count of "significant markers" is meaningful.

## Symbols

- **genetic model** — the encoding of the genotype in the GLM: ADD counts alt alleles (0/1/2), DOM contrasts carriers against non-carriers, REC contrasts alt-homozygotes against the rest. Three separate genome-wide scans, not one joint test; ADD is the primary.

- **OR** — odds ratio per copy of the A1 (effect) allele from logistic regression, with its 95% CI. Compared on the log scale, so a protective and a risk allele of equal strength are equally far from 1.

- **genome-wide significance** — P<5x10^-8, applied identically to every cohort and model. The three cohorts are nested and the three models correlated, so these are not independent tests and no further multiplicity adjustment is made — stated, not silently assumed.

- **ERRCODE** — plink2's per-variant fit diagnostic. Variants with a non-'.' code (e.g. VIF_INFINITE, SEPARATION) did not fit cleanly and are excluded from lambda_GC and from the hit list rather than silently carried.

---

Methods and rationale: [`METHODS.md`](../docs/METHODS.md)
