# HLA allele hits against jMorp 61KJPN-HLA — intermediate_mainland / random

**Figure file:** `freq_qc.random.png`

## The question this figure answers

Are the alleles this scan reports as associated alleles that a large reference population actually carries, at the frequency we call them — or are they typing artefacts and rarity artefacts that an association test cannot distinguish from findings?

## Panels

**(a) Frequency against the reference**

One point per tested allele: the frequency in intermediate_mainland on y against the frequency in jMorp 61KJPN-HLA on x, with the dashed identity line. Alleles clearing the significance threshold are red and named. The narrow left-hand strip holds the alleles the reference has NEVER observed: they have no reference frequency, and drawing them at x = 0 would assert that the reference measured them at zero when it did not measure them at all. Their x positions inside the strip are spread deterministically by rank and carry no meaning.

**(b) Evidence against count**

-log10 P against minor-allele count on a log axis. A hit list that tracks the count floor rather than the data shows up here as evidence concentrating at the left-hand end. Open squares mark alleles absent from the reference, so the two artefact classes can be seen together: a significant allele that is both rare and unknown to the reference is the weakest kind of hit this component can produce.

**(c) Significant and never observed in the reference**

Every allele that clears the threshold and has in_reference = 0, with its frequency in this cohort, its minor-allele count, its call rate, its odds ratio and its P-value. in_reference = 0 means the reference assigned NO chromosome to that allele's P group. When the panel is empty it says so in words: that is a positive result, and a blank panel would read as an analysis that was not run.

## Interpretation

An association test cannot tell a real allele from a mis-called one. If the typing systematically assigns a handful of case chromosomes to an allele they do not carry, the result is a clean odds ratio with a small P-value for an allele that is not there, and no P-value, interval or lambda_GC will reveal it. The only external handle available is population frequency. Panel (a) is that check drawn over every tested allele: agreement with the reference means the typing has not moved the frequency spectrum, and a significant allele sitting on the identity line is one the reference independently agrees exists at that frequency. Panel (b) covers the other artefact class, which is arithmetic rather than biological: at counts in the tens a few chromosomes moving between arms produce a large odds ratio and a small P. Panel (c) is the strongest single caution the component can issue — an allele that 61,424 Japanese individuals have never been observed to carry, called here often enough to be tested and significant, is more likely a typing artefact than a discovery. None of this proves any individual call wrong; it identifies which hits must not be reported without further evidence.

## Values in this rendering

| quantity | value |
|---|---|
| cohort | intermediate_mainland |
| model | random |
| marker class | allele |
| reference panel | jMorp 61KJPN-HLA (61,424 individuals) |
| alleles with a usable fit and a QC row | 159 |
| fits excluded on ERRCODE | 0 |
| alleles in the sumstats file | 159 |
| allele rows in marker_qc.tsv | 2,175 |
| alleles with a reference frequency | 136 |
| alleles the reference has never observed | 23 |
| significance threshold | 4.425e-04 |
| threshold source | given with --alpha |
| significant alleles | 1 |
| significant AND absent from the reference | 1 |

## Full statistics

**Every significant allele, with its reference frequency and count**

| marker | gene | af | freq_reference | in_reference | mac | maf | call_rate | OR | L95 | U95 | P |
|---|---|---|---|---|---|---|---|---|---|---|---|
| HLA_F*01 | F | 0.005582 | — | 0 | 28 | 0.005582 | 1 | 3.043 | 1.678 | 5.52 | 1.247e-04 |

**Significant alleles the reference has never observed**

| marker | gene | af | freq_reference | in_reference | mac | maf | call_rate | OR | L95 | U95 | P |
|---|---|---|---|---|---|---|---|---|---|---|---|
| HLA_F*01 | F | 0.005582 | — | 0 | 28 | 0.005582 | 1 | 3.043 | 1.678 | 5.52 | 1.247e-04 |

## How to read it

1. Start with panel (c). If it names an allele, that allele does not go into a result sentence without independent evidence, whatever panels (a) and (b) show.
2. In (a), read the red points only. A red point on the dashed line is a significant allele whose frequency the reference agrees with; a red point far above it is an allele we call more often than the reference does, which is what a typing bias toward a common allele looks like.
3. In (b), read left to right. Evidence that appears only at the smallest counts is evidence about the count floor.
4. The left-hand strip in (a) is not "frequency zero in the reference". It is "the reference never observed this allele", which is a different statement and the reason the strip is separated from the axis.
5. Nothing here measures per-sample typing accuracy. See the typing component for what the frequency comparison can and cannot establish.

## What this figure does *not* establish

- It compares POPULATION frequencies. A set of typing errors that happens to preserve the frequency spectrum is invisible to it, and it can never say that a particular sample was typed correctly.
- Agreement with the reference is not evidence that an association is real. It only removes one specific way for it to be false.
- The reference is 61,424 individuals from a Japanese panel matched on IPD-IMGT P groups. A P group is coarser than a 2-field allele, so an allele that differs from a reference allele only below P-group resolution is compared at the group level, not at its own.
- in_reference = 0 is a flag, not a verdict. A genuinely rare allele absent from the reference by sampling alone is possible; the flag says the hit needs external evidence, not that it is wrong.
- Residue markers are excluded entirely. The reference is an allele panel, so a blank freq_reference on a residue means "not applicable", not "never observed", and merging the two would manufacture flags.

## Symbols

- **OR** — odds ratio per copy of the A1 (effect) allele from logistic regression, with its 95% CI. Compared on the log scale, so a protective and a risk allele of equal strength are equally far from 1.

- **EAF** — frequency of the EFFECT allele (plink2 A1) among called samples of that group. Case and control EAF are the two numbers the odds ratio is computed from, so plotting them directly shows what drives an effect estimate and whether one group carries the whole difference.

- **ERRCODE** — plink2's per-variant fit diagnostic. Variants with a non-'.' code (e.g. VIF_INFINITE, SEPARATION) did not fit cleanly and are excluded from lambda_GC and from the hit list rather than silently carried.

## Model

```
logit,Pr(case_i) = beta_0 + beta,g_i + gamma_sex,SEX_i + sum_k=1^Kgamma_k,PC_k,i      (g_i = genotype under the stated model)
```

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
