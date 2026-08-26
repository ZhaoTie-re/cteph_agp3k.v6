# HLA allele frequencies against the 1000 Genomes JPT panel

**Figure file:** `allele_frequency.png`

## The question this figure answers

Do the alleles this component calls occur at the frequencies published for a Japanese population — or has the typing moved the spectrum?

## Panels

**(a) HLA-A**

Every 2-field allele seen at HLA-A in the controls or in the reference. x is the JPT frequency, y is ours, and the dashed line is equality. Error bars are 95 % Wilson intervals on the observed frequency — a normal-approximation interval goes below zero for the rare alleles that make up most of an HLA table. Named in red: the commonest alleles and the two furthest from the line.

**(b) HLA-B**

Every 2-field allele seen at HLA-B in the controls or in the reference. x is the JPT frequency, y is ours, and the dashed line is equality. Error bars are 95 % Wilson intervals on the observed frequency — a normal-approximation interval goes below zero for the rare alleles that make up most of an HLA table. Named in red: the commonest alleles and the two furthest from the line.

**(c) HLA-C**

Every 2-field allele seen at HLA-C in the controls or in the reference. x is the JPT frequency, y is ours, and the dashed line is equality. Error bars are 95 % Wilson intervals on the observed frequency — a normal-approximation interval goes below zero for the rare alleles that make up most of an HLA table. Named in red: the commonest alleles and the two furthest from the line.

**(d) HLA-DQB1**

Every 2-field allele seen at HLA-DQB1 in the controls or in the reference. x is the JPT frequency, y is ours, and the dashed line is equality. Error bars are 95 % Wilson intervals on the observed frequency — a normal-approximation interval goes below zero for the rare alleles that make up most of an HLA table. Named in red: the commonest alleles and the two furthest from the line.

**(e) HLA-DRB1**

Every 2-field allele seen at HLA-DRB1 in the controls or in the reference. x is the JPT frequency, y is ours, and the dashed line is equality. Error bars are 95 % Wilson intervals on the observed frequency — a normal-approximation interval goes below zero for the rare alleles that make up most of an HLA table. Named in red: the commonest alleles and the two furthest from the line.

**(f) Concordance**

Pearson $r$ and Spearman $\rho$ between the two frequency vectors per locus. $r$ is dominated by the common alleles, $\rho$ weights every allele equally, so the two disagreeing means the disagreement is in the tail.

## Interpretation

Agreement here is evidence that the typing is not systematically mis-assigning alleles — the one thing every other QC output in this component is blind to, because call rate and field depth measure whether a call was produced, not whether it was right. Disagreement at a COMMON allele is the informative failure: it means a frequent haplotype is being read as something else. Disagreement in the tail is expected and mostly reflects the reference's size.

## What this figure does *not* establish

- This is a comparison of POPULATION FREQUENCIES. It cannot say that any given sample was typed correctly, and a set of errors that happens to preserve the frequency spectrum is invisible to it. Only typing samples with known types measures accuracy — see docs/OPEN_QUESTIONS.md §2.
- The reference is 105 samples. The standard error on a frequency near 0.4 is about 0.034, so this bounds gross error and nothing finer.
- Five loci only — A, B, C, DQB1, DRB1. The reference carries no DPB1, no DQA1 and no DRB3/4/5, so the loci with the least reliable typing in this component are exactly the ones it cannot check.
- Cases are excluded from the comparison and shown only in the table. They are a disease series, and the MHC is where a disease series is least expected to match a population reference.

---

Methods and rationale: [`METHODS.md`](../../docs/METHODS.md)
