# The genome-wide MHC signal and the HLA signal

**Main figure.** **Covers:** `gwas_link.png` in `04.gwas_link/`

## The question this figure answers

Is the genome-wide scan's MHC lead the same signal as an HLA residue or allele?

## Panels

**(gwas_link.png a) the GWAS lead per cohort**

Marginal −log10 P (black), given the first residue signal's lead (red), and given the HLA marker in highest LD with it (open square); dotted: the GWAS suggestive line 1e-5.

**(gwas_link.png b) the HLA signal lead per cohort**

Marginal (filled) and given the GWAS lead (open); tick: the cohort's residue GEC threshold.

**(gwas_link.png c-e) every single-marker conditioning, per cohort**

x: r² of the conditioning marker with the GWAS lead (ALT dosage); y: the lead's −log10 P given that marker. Markers tried: the 10 in highest |r|, the 10 with the smallest marginal P and the signal leads. Dashed: the lead's marginal −log10 P.

## Interpretation

All cases, random model. The GWAS lead chr6:32632861:A:G has marginal P 7.0e-06 / 1.8e-06 / 1.5e-06 here (mean-imputed; the genome-wide scan, with best_guess imputation, reports a slightly smaller P, METHODS §7.1). narrow: GWAS lead P marginal / given DQB1:57:V / given best LD tag (DRB1:10:Q, r2 0.47): 7.0e-06 / 6.2e-03 / 1.5e-06. narrow: DQB1:57:V P marginal / given the GWAS lead: 4.6e-05 / 4.8e-02. intermediate: GWAS lead P marginal / given DQB1:57:V / given best LD tag (DRB1:-25:K, r2 0.49): 1.8e-06 / 9.8e-04 / 7.0e-06. intermediate: DQB1:57:V P marginal / given the GWAS lead: 1.8e-04 / 1.5e-01. full: GWAS lead P marginal / given DQB1:57:V / given best LD tag (DRB1:-16:A, r2 0.49): 1.5e-06 / 8.0e-04 / 3.4e-06. full: DQB1:57:V P marginal / given the GWAS lead: 1.7e-04 / 1.4e-01. A marker that explains the lead pulls its P towards 1; one in LD that does not leaves it near the marginal line.

## Values in this rendering

| quantity | value |
|---|---|
| gwas_link: narrow: GWAS lead P marginal / given DQB1:57:V / given best LD tag (DRB1:10:Q, r2 0.47) | 7.0e-06 / 6.2e-03 / 1.5e-06 |
| gwas_link: narrow: DQB1:57:V P marginal / given the GWAS lead | 4.6e-05 / 4.8e-02 |
| gwas_link: intermediate: GWAS lead P marginal / given DQB1:57:V / given best LD tag (DRB1:-25:K, r2 0.49) | 1.8e-06 / 9.8e-04 / 7.0e-06 |
| gwas_link: intermediate: DQB1:57:V P marginal / given the GWAS lead | 1.8e-04 / 1.5e-01 |
| gwas_link: full: GWAS lead P marginal / given DQB1:57:V / given best LD tag (DRB1:-16:A, r2 0.49) | 1.5e-06 / 8.0e-04 / 3.4e-06 |
| gwas_link: full: DQB1:57:V P marginal / given the GWAS lead | 1.7e-04 / 1.4e-01 |

## Full statistics

**The GWAS lead given each single HLA marker (all cases, random model)**

| conditioning marker | why | r2 with lead (narrow) | lead P given it N / I / F |
|---|---|---|---|
| DRB1:10:Q | ld_tag | 0.469 | 1.5e-06 / 2.2e-06 / 1.9e-06 |
| DQA1:175:E | top_p | 0.467 | - / 8.3e-04 / 5.8e-04 |
| DQA1:175:E | ld_tag,top_p | 0.467 | 2.8e-03 / - / - |
| DRB1:-16:A | ld_tag | 0.467 | 4.8e-06 / 7.0e-06 / 3.4e-06 |
| DRB1:-25:K | ld_tag | 0.467 | 4.8e-06 / 7.0e-06 / 3.4e-06 |
| DRB1:-25:R | ld_tag | 0.467 | 4.8e-06 / 7.0e-06 / 3.4e-06 |
| DRB1:-16:V | ld_tag | 0.462 | 4.2e-06 / 5.5e-06 / 2.6e-06 |
| DRB1:10:Y | ld_tag | 0.457 | 7.4e-07 / 9.6e-07 / 6.9e-07 |
| DRB1:11:S | ld_tag | 0.457 | 7.4e-07 / 9.6e-07 / 6.7e-07 |
| DRB1:12:K | ld_tag | 0.457 | 7.4e-07 / 9.6e-07 / - |
| DRB1:12:T | ld_tag | 0.457 | 7.4e-07 / 9.6e-07 / 6.9e-07 |
| DQB1:87:Y | top_p | 0.412 | 1.6e-03 / 3.1e-04 / 2.3e-04 |
| DQB1:185:T | top_p | 0.399 | 7.0e-04 / 1.9e-04 / 1.1e-04 |
| DQB1:57:V | top_p,signal_lead | 0.314 | 6.2e-03 / 9.8e-04 / 8.0e-04 |
| DQA1*01:01 | top_p | 0.124 | 2.8e-04 / - / 4.3e-05 |
| DRB1:11:L | top_p | 0.124 | 2.2e-04 / - / - |
| B:114:D | top_p | 0.117 | 3.1e-04 / 8.1e-05 / 1.0e-04 |
| B:114:N | top_p | 0.117 | 3.0e-04 / 8.1e-05 / 1.0e-04 |
| B*40:02 | top_p | 0.003 | 1.3e-05 / 4.2e-06 / - |
| A*31:01 | top_p | 0.002 | 1.3e-05 / 3.5e-06 / 2.9e-06 |
| A:62:Q | top_p | — | - / - / 3.6e-06 |
| DQA1:215:F | top_p | — | - / 1.0e-04 / - |
| DQB1:185:I | top_p | — | - / 1.9e-04 / 1.1e-04 |
| DRB1:149:H | ld_tag | — | - / - / 1.4e-06 |
| DRB1:149:Q | ld_tag | — | - / 3.3e-06 / 1.4e-06 |

## How to read it

1. Panel a: does conditioning on the HLA signal pull the lead down in every cohort, and the LD tag not?
2. Panel b: does the HLA signal lose its significance given the lead (the reverse)?
3. Panels c-e: is the marker that explains the lead the one with the highest r²?

## What this figure does *not* establish

- Which variant is causal: two markers that explain each other are one signal statistically.
- The lead is a SNP between HLA-DRB1 and HLA-DQA1; its relation to DQB1 runs through haplotypes.
- Replication: the three cohorts are nested and share almost all cases; n/3 is a sample-selection sensitivity analysis, not independent replication.
- Platform: all cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X; a typing difference between platforms would look like an association, and no covariate separates the two. Within cases the platform mix also differs by DVT status.

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
