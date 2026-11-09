# assoc_hla.v2 — study notes (run of 2026-10-23)

These are measured facts from `results/`. Every number below can be traced to a table named beside it. Run 16: 294 tasks, none failed (run 10 computed the thresholds and conditional analyses; run 12 redrew the figures; run 13 added the benchmark MAXT_BOOT; run 14 the omnibus tests; run 15 removed the gene-level omnibus tests; run 16 added the GWAS lead given whole positions). `verify.sh`: 796 checks passed, 0 failed.

**Reading rules (METHODS §5).**

- **Thresholds.** The primary threshold is LD-adjusted: 0.05 / Me, with Me from GEC v0.2. Exact max-T is the validation threshold and conventional Bonferroni the sensitivity threshold. BH is not used.
- **Model.** Two main analyses: single markers, in the random model (primary for single markers; the fixed model is their comparison); and whole residue positions, in fixed-model omnibus tests at Bonferroni (§10). Genes are not tested as a whole.
- **Consistency.** A marker is reported as significant in n/3 nested cohorts, with the same direction of effect. 1/3 is called "Single cohort".
- **Not replication.** The cohorts share almost all cases, so consistency across them is a sample-selection sensitivity analysis.

## 1. Tests and thresholds

Source: `<case_set>/<cohort>/01.markers/<class>/effective_tests.<class>.tsv`.

Case set `all` (narrow / intermediate / full mainland):

| | alleles | residues |
|---|---|---|
| markers | 143 / 144 / 142 | 672 / 681 / 679 |
| mirror pairs | 0 | 178 / 178 / 167 |
| tests | 143 / 144 / 142 | 494 / 503 / 512 |
| Me (GEC v0.2) | 114.0 / 115.0 / 113.6 | 201.8 / 200.8 / 200.9 |
| **LD-adjusted (GEC)** | **4.39e-4 / 4.35e-4 / 4.40e-4** | **2.48e-4 / 2.49e-4 / 2.49e-4** |
| exact (max-T) | 4.14e-4 / 4.15e-4 / 4.13e-4 | 2.77e-4 / 2.73e-4 / 2.74e-4 |
| conventional Bonferroni | 3.50e-4 / 3.47e-4 / 3.52e-4 | 1.01e-4 / 9.94e-5 / 9.77e-5 |

**Family-wise error rate under max-T.**

- The GEC threshold gives 0.045–0.053 over all 24 (case set, cohort, class) families. That is 0.052–0.053 for alleles and 0.045–0.047 for residues.
- Conventional Bonferroni gives 0.041–0.043 for alleles and 0.019–0.021 for residues, so it is twice as strict as needed for residues.

**Other case sets.**

- The DVT-yes and DVT-no thresholds are within a few percent of these.
- `dvt_case_only` has about 300 cases, so more residue pairs are complements (215–218). It has 457–463 residue tests, Me 185–187 and a GEC threshold of 2.68–2.70e-4.

**Earlier implementation.** GEC as implemented in KGGSEE v1.1 gave smaller Me for case set `all` (alleles about 113, residues about 170). Its residue threshold, about 2.9e-4, lies above the exact max-T threshold. The original v0.2 is used (`tools/gec-0.2/SOURCE.md`).

## 2. Results by case set

Source: `<case_set>/_comparison/hits.tsv`, random model. These are counts of markers by cross-cohort consistency.

| case set | class | LD-adjusted (GEC) | exact (max-T) | conventional Bonferroni |
|---|---|---|---|---|
| all | residue | DQB1 57 V: **3/3** | 3/3 | 1/3 (narrow) |
| DVT yes | both | none | none | none |
| DVT no | residue | 51 markers (33 tests), all 1/3 (narrow): one DQ signal (§4) | the same 51, 1/3 | none |
| DVT yes vs no, cases only | both | none | none | none |

No marker has significant cohorts that disagree in direction.

## 3. All cases

Source: `all/_comparison/assoc_all.tsv`, `marker_summary.tsv`.

**DQB1 position 57 = valine is the only marker significant in the random model, and it is 3/3 at the LD-adjusted threshold.** Its frequency is 10.0–10.3 % in cases against 15.3–16.8 % in controls.

| cohort | random OR (95 % CI) | random P | fixed OR (95 % CI) | fixed P |
|---|---|---|---|---|
| narrow | 0.62 (0.49–0.78) | 4.57e-5 | 0.59 (0.46–0.76) | 3.54e-5 |
| intermediate | 0.65 (0.52–0.82) | 1.85e-4 | 0.63 (0.49–0.80) | 1.48e-4 |
| full | 0.66 (0.53–0.82) | 1.74e-4 | 0.63 (0.50–0.79) | 1.08e-4 |

- **Tested alleles that carry 57 V** are DQB1\*05:01 and DQB1\*06:04 (`_conditional/signals_all.tsv`).
- **No mirror column.** Position 57 has four tested residues (V, A, S, D), so DQB1:57:V is its own test.
- **Significant in the fixed model only, one cohort each.** These do not reach the threshold in the primary model:
  - B\*40:02, narrow: fixed P 3.74e-4, random 5.09e-4, against the allele threshold 4.39e-4;
  - B:114:N and B:114:D, full: fixed P 1.5e-4 / 1.6e-4, random 3.6e-4 / 3.8e-4, against the residue threshold 2.49e-4.

## 4. DVT sets and the case-only test

Source: `_case_sets/dvt_contrast.tsv`, `dvt_no/_comparison/hits.tsv`.

**Cases without DVT carry the DQ signal; cases with DVT do not reach any threshold.**

- **DVT no, residues.**
  - 51 markers (28 DQA1, 23 DQB1; 33 tests) are significant in narrow mainland only, and none in intermediate or full.
  - Conditioning shows they are one signal, led by DQB1 87 Y (§5): random OR 0.54 / 0.57 / 0.58, P 1.19e-4 / 3.29e-4 / 3.15e-4.
  - Intermediate and full miss the GEC threshold (2.50e-4) by a factor of about 1.3.
- **DVT no, alleles.** DPB1\*05:01 is a risk allele at 1/3 under all three thresholds (narrow, P 2.81e-4). It gives OR 1.63 / 1.55 / 1.51 with P 2.81e-4 / 8.44e-4 / 1.41e-3. Its frequency is 45 % in DVT-no cases against 34–35 % in controls.
- **DVT yes.** No marker reaches any threshold. The smallest P is 1.6e-3 for alleles (C\*03:04) and 1.3e-3 for residues (B:114:D, fixed).
- **Case-only (DVT yes vs no).** No marker reaches any threshold. The smallest P is DQB1:-27:A, random P 3.6e-4 / 4.6e-4 / 1.1e-3 against thresholds of about 2.7e-4.

**The leads across case sets** (random model, OR and P narrow / intermediate / full):

| marker | all | DVT yes | DVT no | case-only (DVT yes vs no) |
|---|---|---|---|---|
| DQB1:57:V | 0.62–0.66; 4.6e-5–1.9e-4 | 0.64–0.69; 0.011–0.027 | 0.52–0.55; 3.6e-4–1.0e-3 | 1.43–1.47; 0.22–0.24 |
| DQB1:87:Y | 0.72–0.74; 1.2e-3–1.9e-3 | 0.81–0.84; 0.16–0.22 | 0.54–0.58; 1.2e-4–3.3e-4 | 1.70–1.91; 0.012–0.029 |
| DPB1\*05:01 | 1.17–1.27; 4.2e-3–0.043 | 1.11–1.28; 0.058–0.37 | 1.51–1.63; 2.8e-4–1.4e-3 | 0.75–0.77; 0.088–0.14 |

**What the case-only test shows.**

- DQB1 57 V protects in both DVT sets, so the case-only test gives no difference.
- For DQB1 87 Y and DPB1\*05:01, the case-only OR points the way the DVT-no effect is stronger: Y is rarer and DPB1\*05:01 more common in cases without DVT. Neither reaches a threshold (nominal P 0.012–0.029 and 0.09–0.14).
- The DVT sets are small (140–159 cases). The platform mix differs between cases with and without DVT (METHODS §2.1).

## 5. Independent signals

Source: `_conditional/signals_all.tsv`, `steps_all.tsv`. Signals are selected within a class, jointly over the three cohorts (METHODS §7.4).

| case set | model | class | signal | consistency at entry (GEC / max-T / Bonferroni) | members | tested carrier alleles |
|---|---|---|---|---|---|---|
| all | random | residue | DQB1:57:V | 3/3 / 3/3 / 1/3 | 1 | DQB1\*05:01, DQB1\*06:04 |
| DVT no | random | residue | DQB1:87:Y | 1/3 / 1/3 / 0/3 | 51 | DQB1\*05:01, \*05:02, \*05:03, \*06:04 |
| DVT no | random | allele | DPB1\*05:01 | 1/3 / 1/3 / 1/3 | 1 | — |

**One signal per walk.** After the first signal entered, no marker remained significant in any cohort:

- all cases: the next lead is DPB1:65:I, largest P 0.024;
- DVT-no residues: A:99:Y, 0.0066;
- DVT-no alleles: DQB1\*04:02, 0.0026.

**Fixed-model walks.** The fixed model gives the same residue signals: DQB1:57:V 3/3, and DQB1:87:Y 1/3 with 46 members. It also gives B\*40:02 (1/3) as an allele signal in all cases (§3).

## 6. The GWAS MHC signal is the DQB1 57 V signal

Source: `_conditional/consistency.tsv`, `summary.tsv`. The lead is chr6:32632861:A:G (effect allele G), between HLA-DRB1 and HLA-DQA1. The GWAS P here is mean-imputed, so it is slightly weaker than the genome-wide scan's (METHODS §7.1).

All cases, random model (narrow / intermediate / full):

| quantity | value |
|---|---|
| GWAS lead, marginal P | 7.0e-6 / 1.8e-6 / 1.5e-6 |
| GWAS lead given DQB1:57:V | 6.2e-3 / 9.8e-4 / 8.0e-4 |
| GWAS lead given the whole position 57 (Ala, Ser, Val; Asp left out) | 6.3e-3 / 9.2e-4 / 6.5e-4 |
| DQB1:57:V given the GWAS lead | 0.048 / 0.15 / 0.14 |
| the single HLA marker that explains the lead most | DQB1:57:V in all three cohorts |
| GWAS lead given its best LD tag | 1.5e-6 / 7.0e-6 / 3.4e-6 (DRB1:10:Q / DRB1:-25:K / DRB1:-16:A) |

**Reading.**

- Conditioning on DQB1 57 V lowers the GWAS lead's evidence by about three orders of magnitude, and the reverse removes most of DQB1 57 V's. They are one signal.
- The markers in highest LD with the lead (DRB1) do not explain it. LD alone does not identify the HLA marker behind the GWAS signal.
- **The whole position 57 explains the lead to the same degree as V alone**, so V carries the position's part. Some evidence (P about 1e-3) remains under either; part of the peak may tag other DR–DQ variation that no single HLA marker accounts for.
- **DVT no:** the lead (marginal P 3.9e-4 / 4.5e-4 / 5.1e-4) given DQB1:87:Y is 0.13 / 0.095 / 0.10; given the whole position 87 (Phe, Tyr; Leu left out) 0.29 / 0.20 / 0.21; given the whole position 125, the omnibus lead (Gly, Ser; Ala left out), 0.18 / 0.14 / 0.13.
- **Case-only:** the lead has no association (P 0.55 / 0.74 / 0.86).

## 7. Fixed against random

Source: `<case_set>/_comparison/model_concordance.tsv`.

- **Agreement.** r(−log10 P) is 0.993–0.999 and r(β) is 0.962–0.998 in every case set, cohort and class.
- **Where they differ.** The random model's P is slightly larger. The median −log10 P ratio is 0.91–1.00.
- **Calls at the LD-adjusted threshold.**

  | case set | markers significant, fixed / random / both |
  |---|---|
  | all | 6 / 3 / 3 |
  | DVT no | 47 / 52 / 47 |

  The fixed-only calls are listed in §3.

## 8. λGC is a property of the region, not of the test

λGC (case set `all`) is 1.37–1.91 for alleles and 2.00–2.37 for residues. For reference, the genome-wide fixed scan (`analysis/assoc_plink2`, full mainland, the same covariates, SNPs with MAF ≥ 0.5 %) gives λGC 1.89 in chr6:29.6–33.4 Mb (16,656 SNPs) against 1.13 in the rest of the genome (5,098,423 SNPs). The HLA markers inflate as much as the MHC SNPs of the same scan. λGC is reported and never applied.

## 9. Alleles previously reported in Japanese CTEPH

Source: `all/_comparison/assoc_all.tsv`.

| allele | case / control frequency (full, fixed) | OR range over the 6 cells | P range |
|---|---|---|---|
| B\*52:01 | 12.0 % / 10.8 % | 1.02–1.08 | 0.50–0.87 |
| DPB1\*02:02 | 3.7 % / 4.1 % | 0.86–0.96 | 0.46–0.82 |

Neither is associated in these data.

## 10. Omnibus tests of residue positions (a main analysis; fixed model, Bonferroni)

Source: `_omnibus/omnibus_units.tsv`, `omnibus_signals_all.tsv` (METHODS §8). P values below are narrow / intermediate / full mainland. The Bonferroni thresholds are 1.69–1.72 × 10⁻⁴ (290–296 positions tested per cohort).

| case set | class | entered unit | P | n/3 | members | random-model single-marker lead |
|---|---|---|---|---|---|---|
| all | residue | DQB1 position 57 (A, S, V against D; df 3) | 3.3e-5 / 1.6e-4 / 8.4e-5 | 3/3 | DQB1:57, B:114 | DQB1:57:V, in the unit |
| DVT no | residue | DQB1 position 125 (G, S against A; df 2) | 1.8e-5 / 5.0e-5 / 7.7e-5 | 3/3 | 26 units in DQA1 and DQB1 | DQB1:87:Y, among the members (position 87 is itself 3/3) |

No position is significant after conditioning on the entered one, and DVT yes and the case-only test have no significant position.

**Reading.**

- **Residues agree with the single-marker scans.** The omnibus lead holds the random-model signal lead (DQB1 57 V) or holds it among its members (DQB1 87 Y). Position 57 tested as a whole (P 8.4e-5, full) is about as strong as V alone. Conditioning on the whole position leaves no further signal, so conditioning on V alone did not miss one.
- **Without DVT, the DQ positions are significant in all three cohorts.** Positions 125, 87, −27 and −9 of DQB1 are 3/3, although the single residues reach the GEC threshold in narrow mainland only. The position threshold (1.69–1.72 × 10⁻⁴) is stricter than the single-residue GEC threshold (2.49–2.50 × 10⁻⁴), so the gain comes from the test: at position 87, F and Y both carry a lower risk than L (OR 0.76 and 0.47, full), and the position test uses both. In the same fixed model, Y alone has P 1.3e-4 / 3.2e-4 / 3.1e-4 and position 87 has P 2.2e-5 / 8.3e-5 / 9.4e-5, so the gain is not a change of model.
- **No gene-level omnibus tests (removed 2026-10-25, run 15).** Testing all alleles of a gene jointly is not standard in HLA fine-mapping (Hirata et al. 2019 test positions only). It also depends on the `unlisted` category, the copies of alleles outside the tested matrix, whose share is lower in cases than in controls in every gene, including narrow mainland:

  | gene | unlisted copies, cases / controls (full) | ratio, full | ratio, narrow |
  |---|---|---|---|
  | A | 9.7 % / 14.3 % | 0.68 | 0.69 |
  | B | 13.7 % / 15.7 % | 0.87 | 0.90 |
  | C | 10.1 % / 11.8 % | 0.86 | 0.87 |
  | DRB1 | 2.6 % / 5.1 % | 0.51 | 0.56 |
  | DQA1 | 1.5 % / 2.0 % | 0.76 | 0.63 |
  | DQB1 | 2.5 % / 3.5 % | 0.71 | 0.65 |
  | DPB1 | 3.7 % / 6.2 % | 0.59 | 0.54 |

  Computed on 2026-10-25 from `hla.typing.v2/results/06.assoc_inputs/<cohort>/allele_dosage.tsv` (samples called at every allele of the gene) and the cohort phenotype file. A shift in the same direction across unrelated genes is what a platform difference in typing would produce. The run-14 gene-level results (all cases: HLA-A and DRB1 at 3/3, with no single allele at the GEC threshold; DRB1's `unlisted` OR 0.51) are kept in `logs/stale_before_run15/` and are not reported.
- **Fits.** 31 of 3,517 position × cohort tests separated (17 DVT yes, 14 case-only). Their P is kept but never counts as significant. 35 positions were not tested in a cohort.

## 11. Limits that apply to every result above

- **Platform.** All cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X (15×). A typing difference between platforms would look like an association, and no covariate separates the two. A platform-sensitivity and down-sampling check of these results is planned separately.
- **Platform within cases.** Within cases, the platform mix differs by DVT status (χ² P 9e-7). This affects the DVT sets and the case-only test.
- **Nested cohorts.** "3/3" means a result survives adding controls; it is not replication.
