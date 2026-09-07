# HLA allele frequencies against jMorp 61KJPN-HLA

**Figure file:** `allele_frequency.png`

## The question this figure answers

Do the alleles this component calls occur at the frequencies published for a Japanese population — or has the typing moved the spectrum?

## Panels

**(a) All loci pooled**

Every 2-field allele seen in the controls or in the reference, at every locus, on one axis. x is the reference frequency, y is ours, the dashed line is equality. Error bars are 95 % Wilson intervals on the observed frequency. Named in red: the two commonest and the three furthest from the line. The per-locus version is allele_frequency_loci.png.

**(b) Per-locus concordance**

Pearson $r$ between the two frequency vectors, per locus, with a dotted line at 0.98. There is deliberately no rank correlation beside it: most of each vector is an allele only one source carries, so a rank statistic over it measures panel coverage rather than agreement.

**(c) Unconfirmed share by locus**

The fraction of control chromosomes carrying an allele the reference does not list, per locus, on panel (b)'s rows. Pooled over loci it reproduces the cohort figure in typing_confound.png exactly.

## Interpretation

Agreement here is evidence that the typing is not systematically mis-assigning alleles — the one thing every other QC output in this component is blind to, because call rate and field depth measure whether a call was produced, not whether it was right. Disagreement at a COMMON allele is the informative failure: it means a frequent haplotype is being read as something else. Disagreement in the tail is expected and mostly reflects which alleles the reference happens to carry. Panels (b) and (c) are deliberately on the same rows: a locus can have high concordance and still put a large share of chromosomes on alleles the reference never lists, and the two together say which of those is happening.

## Values in this rendering

| quantity | value |
|---|---|
| loci compared | 13 |
| loci with r >= 0.98 | 11 |
| control samples | 2,662 |
| reference individuals | 61,424 |
| lowest r | 0.6352 (DRB3) |
| highest r | 1.0000 |
| worst unconfirmed locus | DRB4 19.2% |
| best unconfirmed locus | DQA1 0.5% |
| alleles compared | 2,452 |
| alleles both sources carry, fewest | 2 |
| alleles both sources carry, most | 114 |

## Full statistics

**Per-locus concordance and unconfirmed share**

| gene | n_alleles | n_chr_ctrl | n_chr_ref | pearson_r | n_shared | max_abs_diff | allele_at_max | n_ref_only | n_obs_only | unconfirmed_share |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 552 | 5,324 | 121,971 | 0.993 | 114 | 0.0299 | A*26:01P | 233 | 205 | 0.06067 |
| B | 546 | 5,324 | 122,813 | 0.9755 | 103 | 0.0313 | B*40:02P | 227 | 216 | 0.0618 |
| C | 442 | 5,322 | 122,800 | 0.9869 | 91 | 0.0347 | C*03:04P | 178 | 173 | 0.04961 |
| DPA1 | 33 | 5,324 | 122,848 | 0.9998 | 6 | 0.0125 | DPA1*01:03P | 5 | 22 | 0.02085 |
| DPB1 | 193 | 5,324 | 122,848 | 0.9966 | 41 | 0.0401 | DPB1*02:01P | 46 | 106 | 0.03343 |
| DQA1 | 27 | 5,324 | 122,847 | 0.9957 | 8 | 0.0254 | DQA1*03:01P | 1 | 18 | 0.004696 |
| DQB1 | 128 | 5,324 | 122,848 | 0.9874 | 34 | 0.0253 | DQB1*06:01P | 52 | 42 | 0.009391 |
| DRB1 | 285 | 5,324 | 122,287 | 0.9828 | 56 | 0.0234 | DRB1*15:13 | 147 | 82 | 0.01972 |
| DRB3 | 101 | 2,540 | 122,848 | 0.6352 | 10 | 0.4852 | DRB3*01:01P | 24 | 67 | 0.1669 |
| DRB4 | 78 | 3,028 | 56,082 | 0.9909 | 3 | 0.2083 | DRB4*01:01P | 2 | 73 | 0.1915 |
| E | 34 | 5,324 | 122,848 | 0.9991 | 2 | 0.0262 | E*01:01P | 0 | 32 | 0.03512 |
| F | 15 | 5,324 | 122,848 | 1 | 2 | 0.015 | F*01:01P | 0 | 13 | 0.01446 |
| G | 18 | 5,324 | 122,564 | 0.9896 | 5 | 0.088 | G*01:01P | 2 | 11 | 0.05992 |

## What this figure does *not* establish

- This is a comparison of POPULATION FREQUENCIES. It cannot say that any given sample was typed correctly, and a set of errors that happens to preserve the frequency spectrum is invisible to it. Only typing samples with known types measures accuracy — see docs/OPEN_QUESTIONS.md §2.
- Only the 13 loci the reference carries can be checked at all. Every other locus this component types has no external reference here, so its calls are unverified rather than verified-and-passing.
- The reference is one Japanese population, not Japanese people in general. A locus where our cohort and the reference panel differ in ancestry will differ in frequency for reasons that have nothing to do with typing.
- Cases are excluded from the comparison and shown only in the table. They are a disease series, and the MHC is where a disease series is least expected to match a population reference.

---

Methods and rationale: [`METHODS.md`](../../docs/METHODS.md)
