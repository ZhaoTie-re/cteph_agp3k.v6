# MHC association scan — full_mainland, HLA SEX + 10 bbj_mainland PCs, relatives removed model

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

1 allele and 1 residue marker(s) clear their own Bonferroni threshold. The two counts are NOT additive evidence: an allele marker and the residue markers it carries are the same chromosomes counted twice, so a peak appearing in both tracks at the same gene is one signal seen two ways, not two. Reading the position off the x-axis alone is what (c) exists to prevent — the marker coordinate is gene start plus index, so it identifies a GENE and no more. 1 significant allele marker(s) are absent from the reference allele panel and are ringed. hla.typing found that the controls carry roughly twice the share of such calls as the cases, differentially by phenotype, so a hit on one of them has the shape a typing artefact makes and cannot be reported as a finding without that check.

## Values in this rendering

| quantity | value |
|---|---|
| cohort | full_mainland |
| model | fixed |
| allele markers drawn | 170 |
| residue markers drawn | 705 |
| allele markers excluded | 0 |
| residue markers excluded | 0 |
| allele markers passing QC | 170 |
| residue markers passing QC | 705 |
| genes drawn | 16 |
| alpha Bonferroni allele | 2.941e-04 |
| alpha Bonferroni residue | 7.092e-05 |
| alpha effective | 4.065e-04 |
| alpha genome-wide | 5.000e-08 |
| allele markers significant | 1 |
| residue markers significant | 1 |
| significant alleles absent from the reference panel | 1 |
| markers named in (a) | 1 |
| markers named in (b) | 1 |
| strongest marker | AA_A_178_T |
| smallest P | 7.975e-09 |
| axis (Mb) | 29.67–33.13 |

## Full statistics

**Strongest allele markers**

| ID | gene | position | A1_FREQ | OBS_CT | OR | L95 | U95 | P | in_reference |
|---|---|---|---|---|---|---|---|---|---|
| HLA_F*01 | F |  | 0.004884 | 3,071 | 3.002 | 1.782 | 5.055 | 3.583e-05 | 0 |
| HLA_A*31:01 | A |  | 0.07571 | 3,071 | 1.502 | 1.179 | 1.914 | 9.941e-04 | 1 |
| HLA_DQA1*01:01 | DQA1 |  | 0.06447 | 3,071 | 0.5876 | 0.414 | 0.8339 | 0.002916 | 1 |
| HLA_B*40:02 | B |  | 0.06529 | 3,071 | 1.446 | 1.117 | 1.871 | 0.005069 | 1 |
| HLA_DQB1*05:01 | DQB1 |  | 0.06627 | 3,071 | 0.6238 | 0.4434 | 0.8774 | 0.006709 | 1 |
| HLA_B*51:01 | B |  | 0.07131 | 3,071 | 1.419 | 1.092 | 1.842 | 0.008706 | 1 |
| HLA_A*33:03 | A |  | 0.07962 | 3,071 | 0.6649 | 0.4885 | 0.905 | 0.009469 | 1 |
| HLA_DRB1*01:01 | DRB1 |  | 0.06384 | 3,070 | 0.6357 | 0.4512 | 0.8956 | 0.009585 | 1 |
| HLA_DQB1*04:02 | DQB1 |  | 0.03549 | 3,071 | 1.542 | 1.091 | 2.179 | 0.01409 | 1 |
| HLA_F*01:01 | F |  | 0.9658 | 3,071 | 0.6782 | 0.4942 | 0.9306 | 0.01615 | 1 |
| HLA_C*15:02 | C |  | 0.0241 | 3,070 | 1.666 | 1.095 | 2.536 | 0.01723 | 1 |
| HLA_DRB1*13:02 | DRB1 |  | 0.07068 | 3,070 | 0.6818 | 0.4948 | 0.9394 | 0.01916 | 1 |

**Strongest residue markers**

| ID | gene | position | A1_FREQ | OBS_CT | OR | L95 | U95 | P | in_reference |
|---|---|---|---|---|---|---|---|---|---|
| AA_A_178_T | A | A:178 | 0.9967 | 3,067 | 0.04557 | 0.01596 | 0.1301 | 7.975e-09 | — |
| AA_DQB1_57_V | DQB1 | DQB1:57 | 0.1457 | 3,071 | 0.6267 | 0.495 | 0.7936 | 1.044e-04 | — |
| AA_B_114_N | B | B:114 | 0.6083 | 3,071 | 1.351 | 1.158 | 1.577 | 1.316e-04 | — |
| AA_B_114_D | B | B:114 | 0.3891 | 3,071 | 0.7408 | 0.6347 | 0.8647 | 1.428e-04 | — |
| AA_DQA1_175_E | DQA1 | DQA1:175 | 0.4373 | 3,071 | 1.309 | 1.131 | 1.515 | 3.156e-04 | — |
| AA_DQB1_87_Y | DQB1 | DQB1:87 | 0.2107 | 3,071 | 0.7253 | 0.5985 | 0.8789 | 0.001053 | — |
| AA_DQB1_185_I | DQB1 | DQB1:185 | 0.4086 | 3,037 | 1.26 | 1.088 | 1.46 | 0.00203 | — |
| AA_DQB1_185_T | DQB1 | DQB1:185 | 0.5914 | 3,037 | 0.7934 | 0.6849 | 0.919 | 0.00203 | — |
| AA_A_62_Q | A | A:62 | 0.1897 | 3,070 | 1.306 | 1.099 | 1.553 | 0.0024 | — |
| AA_DQB1_30_H | DQB1 | DQB1:30 | 0.2121 | 3,071 | 0.7562 | 0.6259 | 0.9137 | 0.003793 | — |
| AA_DQA1_215_F | DQA1 | DQA1:215 | 0.6057 | 3,070 | 0.807 | 0.6968 | 0.9347 | 0.004208 | — |
| AA_DQA1_215_L | DQA1 | DQA1:215 | 0.3943 | 3,070 | 1.239 | 1.07 | 1.435 | 0.004208 | — |

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
