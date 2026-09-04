# MHC association scan — intermediate_mainland, HLA SEX + 10 bbj_mainland PCs, relatives removed model

**Figure file:** `scan.fixed.png`

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

1 allele and 0 residue marker(s) clear their own Bonferroni threshold. The two counts are NOT additive evidence: an allele marker and the residue markers it carries are the same chromosomes counted twice, so a peak appearing in both tracks at the same gene is one signal seen two ways, not two. Reading the position off the x-axis alone is what (c) exists to prevent — the marker coordinate is gene start plus index, so it identifies a GENE and no more. 1 significant allele marker(s) are absent from the reference allele panel and are ringed. hla.typing found that the controls carry roughly twice the share of such calls as the cases, differentially by phenotype, so a hit on one of them has the shape a typing artefact makes and cannot be reported as a finding without that check.

## Values in this rendering

| quantity | value |
|---|---|
| cohort | intermediate_mainland |
| model | fixed |
| allele markers drawn | 159 |
| residue markers drawn | 692 |
| allele markers excluded | 0 |
| residue markers excluded | 0 |
| allele markers passing QC | 159 |
| residue markers passing QC | 692 |
| genes drawn | 16 |
| alpha Bonferroni allele | 3.145e-04 |
| alpha Bonferroni residue | 7.225e-05 |
| alpha effective | 4.425e-04 |
| alpha genome-wide | 5.000e-08 |
| allele markers significant | 1 |
| residue markers significant | 0 |
| significant alleles absent from the reference panel | 1 |
| markers named in (a) | 1 |
| markers named in (b) | 3 |
| strongest marker | HLA_F*01 |
| smallest P | 6.066e-05 |
| axis (Mb) | 29.67–33.13 |

## Full statistics

**Strongest allele markers**

| ID | gene | position | A1_FREQ | OBS_CT | OR | L95 | U95 | P | in_reference |
|---|---|---|---|---|---|---|---|---|---|
| HLA_F*01 | F |  | 0.005645 | 2,480 | 3.069 | 1.774 | 5.308 | 6.066e-05 | 0 |
| HLA_B*40:02 | B |  | 0.06149 | 2,480 | 1.55 | 1.187 | 2.024 | 0.0013 | 1 |
| HLA_A*31:01 | A |  | 0.07621 | 2,480 | 1.489 | 1.157 | 1.916 | 0.001999 | 1 |
| HLA_DQA1*01:01 | DQA1 |  | 0.06452 | 2,480 | 0.6016 | 0.4202 | 0.8615 | 0.005542 | 1 |
| HLA_DQB1*05:01 | DQB1 |  | 0.06593 | 2,480 | 0.6214 | 0.4364 | 0.8848 | 0.008313 | 1 |
| HLA_DRB1*01:01 | DRB1 |  | 0.06454 | 2,479 | 0.6394 | 0.45 | 0.9084 | 0.01256 | 1 |
| HLA_A*33:03 | A |  | 0.08347 | 2,480 | 0.6756 | 0.4938 | 0.9243 | 0.01419 | 1 |
| HLA_DPB1*05:01 | DPB1 |  | 0.3516 | 2,480 | 1.208 | 1.038 | 1.405 | 0.01458 | 1 |
| HLA_DRB1*13:02 | DRB1 |  | 0.07443 | 2,479 | 0.6771 | 0.488 | 0.9397 | 0.01969 | 1 |
| HLA_C*03:04 | C |  | 0.1015 | 2,479 | 1.299 | 1.034 | 1.632 | 0.0249 | 1 |
| HLA_F*01:01 | F |  | 0.9643 | 2,480 | 0.6899 | 0.498 | 0.9556 | 0.02556 | 1 |
| HLA_DQB1*04:02 | DQB1 |  | 0.03508 | 2,480 | 1.479 | 1.033 | 2.118 | 0.03262 | 1 |

**Strongest residue markers**

| ID | gene | position | A1_FREQ | OBS_CT | OR | L95 | U95 | P | in_reference |
|---|---|---|---|---|---|---|---|---|---|
| AA_DQB1_57_V | DQB1 | DQB1:57 | 0.1492 | 2,480 | 0.6256 | 0.4912 | 0.7968 | 1.442e-04 | — |
| AA_DQA1_175_E | DQA1 | DQA1:175 | 0.4359 | 2,480 | 1.315 | 1.132 | 1.529 | 3.521e-04 | — |
| AA_B_114_D | B | B:114 | 0.3877 | 2,480 | 0.7537 | 0.6433 | 0.883 | 4.635e-04 | — |
| AA_B_114_N | B | B:114 | 0.6097 | 2,480 | 1.326 | 1.132 | 1.553 | 4.674e-04 | — |
| AA_DQB1_87_Y | DQB1 | DQB1:87 | 0.2115 | 2,480 | 0.7217 | 0.5925 | 0.8789 | 0.001184 | — |
| AA_DQB1_185_I | DQB1 | DQB1:185 | 0.4078 | 2,457 | 1.268 | 1.09 | 1.475 | 0.002138 | — |
| AA_DQB1_185_T | DQB1 | DQB1:185 | 0.5922 | 2,457 | 0.7888 | 0.6779 | 0.9178 | 0.002138 | — |
| AA_DQA1_215_F | DQA1 | DQA1:215 | 0.6065 | 2,479 | 0.7985 | 0.6863 | 0.9291 | 0.003606 | — |
| AA_DQA1_215_L | DQA1 | DQA1:215 | 0.3935 | 2,479 | 1.252 | 1.076 | 1.457 | 0.003606 | — |
| AA_DQA1_50_L | DQA1 | DQA1:50 | 0.3937 | 2,480 | 1.251 | 1.075 | 1.456 | 0.003749 | — |
| AA_DQA1_53_R | DQA1 | DQA1:53 | 0.3937 | 2,480 | 1.251 | 1.075 | 1.456 | 0.003749 | — |
| AA_DQA1_187_A | DQA1 | DQA1:187 | 0.6091 | 2,479 | 0.7987 | 0.686 | 0.9298 | 0.00376 | — |

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
