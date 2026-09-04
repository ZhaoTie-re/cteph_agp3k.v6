# HLA allele frequencies against jMorp 61KJPN-HLA

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

**(d) HLA-DPA1**

Every 2-field allele seen at HLA-DPA1 in the controls or in the reference. x is the JPT frequency, y is ours, and the dashed line is equality. Error bars are 95 % Wilson intervals on the observed frequency — a normal-approximation interval goes below zero for the rare alleles that make up most of an HLA table. Named in red: the commonest alleles and the two furthest from the line.

**(e) HLA-DPB1**

Every 2-field allele seen at HLA-DPB1 in the controls or in the reference. x is the JPT frequency, y is ours, and the dashed line is equality. Error bars are 95 % Wilson intervals on the observed frequency — a normal-approximation interval goes below zero for the rare alleles that make up most of an HLA table. Named in red: the commonest alleles and the two furthest from the line.

**(f) HLA-DQA1**

Every 2-field allele seen at HLA-DQA1 in the controls or in the reference. x is the JPT frequency, y is ours, and the dashed line is equality. Error bars are 95 % Wilson intervals on the observed frequency — a normal-approximation interval goes below zero for the rare alleles that make up most of an HLA table. Named in red: the commonest alleles and the two furthest from the line.

**(g) HLA-DQB1**

Every 2-field allele seen at HLA-DQB1 in the controls or in the reference. x is the JPT frequency, y is ours, and the dashed line is equality. Error bars are 95 % Wilson intervals on the observed frequency — a normal-approximation interval goes below zero for the rare alleles that make up most of an HLA table. Named in red: the commonest alleles and the two furthest from the line.

**(h) HLA-DRB1**

Every 2-field allele seen at HLA-DRB1 in the controls or in the reference. x is the JPT frequency, y is ours, and the dashed line is equality. Error bars are 95 % Wilson intervals on the observed frequency — a normal-approximation interval goes below zero for the rare alleles that make up most of an HLA table. Named in red: the commonest alleles and the two furthest from the line.

**(i) HLA-DRB3**

Every 2-field allele seen at HLA-DRB3 in the controls or in the reference. x is the JPT frequency, y is ours, and the dashed line is equality. Error bars are 95 % Wilson intervals on the observed frequency — a normal-approximation interval goes below zero for the rare alleles that make up most of an HLA table. Named in red: the commonest alleles and the two furthest from the line.

**(j) HLA-DRB4**

Every 2-field allele seen at HLA-DRB4 in the controls or in the reference. x is the JPT frequency, y is ours, and the dashed line is equality. Error bars are 95 % Wilson intervals on the observed frequency — a normal-approximation interval goes below zero for the rare alleles that make up most of an HLA table. Named in red: the commonest alleles and the two furthest from the line.

**(k) HLA-E**

Every 2-field allele seen at HLA-E in the controls or in the reference. x is the JPT frequency, y is ours, and the dashed line is equality. Error bars are 95 % Wilson intervals on the observed frequency — a normal-approximation interval goes below zero for the rare alleles that make up most of an HLA table. Named in red: the commonest alleles and the two furthest from the line.

**(l) HLA-F**

Every 2-field allele seen at HLA-F in the controls or in the reference. x is the JPT frequency, y is ours, and the dashed line is equality. Error bars are 95 % Wilson intervals on the observed frequency — a normal-approximation interval goes below zero for the rare alleles that make up most of an HLA table. Named in red: the commonest alleles and the two furthest from the line.

**(m) HLA-G**

Every 2-field allele seen at HLA-G in the controls or in the reference. x is the JPT frequency, y is ours, and the dashed line is equality. Error bars are 95 % Wilson intervals on the observed frequency — a normal-approximation interval goes below zero for the rare alleles that make up most of an HLA table. Named in red: the commonest alleles and the two furthest from the line.

**(n) Concordance**

Pearson $r$ and Spearman $\rho$ between the two frequency vectors per locus. $r$ is dominated by the common alleles, $\rho$ weights every allele equally, so the two disagreeing means the disagreement is in the tail.

## Interpretation

Agreement here is evidence that the typing is not systematically mis-assigning alleles — the one thing every other QC output in this component is blind to, because call rate and field depth measure whether a call was produced, not whether it was right. Disagreement at a COMMON allele is the informative failure: it means a frequent haplotype is being read as something else. Disagreement in the tail is expected and mostly reflects the reference's size.

## What this figure does *not* establish

- This is a comparison of POPULATION FREQUENCIES. It cannot say that any given sample was typed correctly, and a set of errors that happens to preserve the frequency spectrum is invisible to it. Only typing samples with known types measures accuracy — see docs/OPEN_QUESTIONS.md §2.
- The reference is 61424 samples. The standard error on a frequency near 0.4 is about 0.034, so this bounds gross error and nothing finer.
- Five loci only — A, B, C, DQB1, DRB1. The reference carries no DPB1, no DQA1 and no DRB3/4/5, so the loci with the least reliable typing in this component are exactly the ones it cannot check.
- Cases are excluded from the comparison and shown only in the table. They are a disease series, and the MHC is where a disease series is least expected to match a population reference.

---

Methods and rationale: [`METHODS.md`](../../docs/METHODS.md)
