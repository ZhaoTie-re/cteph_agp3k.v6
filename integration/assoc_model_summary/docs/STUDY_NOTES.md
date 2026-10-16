# assoc_model_summary — study notes

Numbers of the first build (2026-10-02), from `analysis/assoc_plink2/results` and
`analysis/assoc_saige/results` (SAIGE P = `p.value.NA`). They are reproduced by `reduce_scan.py` and checked
by `verify.sh` §7 against the upstream `scan_qc_all.tsv`.

## Samples (cases / controls)

| cohort | fixed (relatives removed) | random |
|---|---|---|
| narrow mainland | 419 / 1,749 | 419 / 1,776 |
| intermediate mainland | 429 / 2,051 | 429 / 2,079 |
| full mainland | 439 / 2,632 | 439 / 2,662 |

## Calibration

| cohort | λGC fixed | λGC GRM + PCs | λGC GRM only | P < 10⁻⁵ fixed / GRM + PCs / GRM only (≈ 51 expected) |
|---|---|---|---|---|
| narrow | 1.120 | 1.039 | 1.040 | 130 / 87 / 180 |
| intermediate | 1.121 | 1.038 | 1.024 | 216 / 160 / 168 |
| full | 1.132 | 1.021 | 1.017 | 204 / 111 / 121 |

By MAF, from 0.01: fixed 1.03–1.15, GRM + PCs 0.95–1.06.

## Genome-wide loci

- Fixed: 0 / 1 / 2 — *FTO* (intermediate, full) and *MIR3681HG* (full).
- GRM + PCs: 0 / 1 / 0 — *FTO* (intermediate, chr16:53882967:G:A, P 3.8 × 10⁻⁸).
- In full mainland under the random model, the smallest P is 1.8 × 10⁻⁷ at the *FTO* locus and 1.0 × 10⁻⁷ at
  *MIR3681HG*.

## PCs on top of the GRM (question 2; conclusion pending)

- The PCs predict case status: LR χ²(10) = 64.4 / 23.8 / 31.0, P = 5.4 × 10⁻¹⁰ / 8.1 × 10⁻³ / 6.0 × 10⁻⁴.
- Peaks only without PCs: 8 / 2 / 3. They include, in the extended MHC, *DDR1* (narrow, 27 variants),
  *XXbac-BPG308K3.5* (narrow) and *HCG17* (intermediate).
- Peaks only with PCs: 6 / 3 / 6.
- Cross-cohort leads: median |Δ ln OR| 0.014, maximum 0.23.

## Found on the way

- About 0.4–0.9% of variants are kept by the control-MAF ≥ 0.01 filter but have an overall MAF below 0.01:
  22,555–31,347 in the fixed scans, 41,060–47,560 in the random scans.
- These variants are inflated under every model: λGC 2.76–2.84 fixed, 1.46–2.01 random.
- This is consistent with a selection effect of filtering on control frequency alone. It is not plotted and is
  reported for follow-up upstream.

## History

Built 2026-10-02 for the model-choice slides. Question 2's argument is to be agreed with the user.
