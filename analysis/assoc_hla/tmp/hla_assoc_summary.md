# HLA association in CTEPH — status report

**Component** `cteph_agp3k.v6/analysis/assoc_hla` · **Date** 2026-08-28 · **Status** analysis complete, not yet written up

---

## 1. Summary

Three findings, in order of confidence.

| # | Finding | Evidence |
|---|---|---|
| **1** | **No HLA allele is associated with CTEPH.** Two previously reported Japanese associations (`B*52:01`, `DPB1*02:02`) do not replicate, at ~100 % power for the reported effect sizes. | Section 3 |
| **2** | **One amino-acid position is associated: DQB1 position 57.** Valine is protective, OR 0.59–0.63, consistent across 3 cohorts × 2 models × 2 test types (6 readings). | Section 4 |
| **3** | **Two typing artefacts were detected and excluded.** Both would have been reported as significant associations by a conventional pipeline; both are sequencing-platform specific and absent from a 61,424-person Japanese reference panel. | Section 6 |

The association is only visible at the residue level — the allele-level scan is null. That asymmetry is the analytical point of the component (Section 5).

---

## 2. Design

HLA dosages from `hla.typing` (HLA-HD, IPD-IMGT/HLA 3.64.0, 2-field, 3,569 samples) encoded as chr6 PLINK markers, so both engines read them unchanged. Round-trip verified exact (10,685,260 genotypes, 0 mismatches).

| Cohort | N | Cases | Controls |
|---|---|---|---|
| `narrow_mainland` | 2,195 | 419 | 1,776 |
| `intermediate_mainland` | 2,508 | 429 | 2,079 |
| `full_mainland` | 3,101 | 439 | 2,662 |

Cohorts are nested and defined by PopGMM without reference to phenotype.

- **Model** logit(p) ~ dosage + SEX + 10 `bbj_mainland` PCs
- **Fixed** plink2 `--glm` (Wald) · **Random** SAIGE full GRM, `--LOCO=TRUE` (chr6 excluded from the null), Firth
- **Layers** allele (203 markers) · residue (705) · amino-acid omnibus per position (300 positions, LRT with *m*−1 df)
- **Threshold** Li & Ji effective-tests, primary; Bonferroni and 5×10⁻⁸ reported alongside
- **Framework** SNP2HLA / Okada 2015 / Hirata 2019, codified in Nat. Protoc. 2023 (doi:10.1038/s41596-023-00853-4)

---

## 3. Allele level — null

**No allele survives conditional analysis in `intermediate_mainland` or `full_mainland`, under either model.**

The single conditional signal anywhere is `HLA_B*40:02` in `narrow_mainland` (OR 1.72, P = 2.0×10⁻⁴). It does not hold in the larger cohorts (`full_mainland`: OR 1.45, P = 5.1×10⁻³). Since narrow and full share 419/439 of the same cases and differ mainly in controls, this pattern points to control heterogeneity rather than a case effect. Not reported as a signal.

### Non-replication of prior Japanese reports

Prior reports were in DVT-negative patients (Kominami 2009 *J Hum Genet* 54:108-114; Tanabe 2005 *Eur Respir J* 25:131-138). DVT status is available for our cases, so the comparison is stratum-matched.

| Marker | Reported OR | Stratum | n cases | Observed OR [95 % CI] | P | Power for reported OR |
|---|---|---|---|---|---|---|
| `B*52:01` | 2.47 | all | 439 | 1.09 [0.86–1.38] | 0.48 | 1.00 |
| `B*52:01` | 2.47 | DVT-negative | 145 | 1.13 [0.77–1.65] | 0.53 | 1.00 |
| `B*52:01` | 2.47 | DVT-positive | 159 | 1.02 [0.70–1.48] | 0.92 | 1.00 |
| `DPB1*02:02` | 5.07 | all | 439 | 1.01 [0.69–1.46] | 0.98 | 1.00 |
| `DPB1*02:02` | 5.07 | DVT-negative | 145 | 0.89 [0.47–1.69] | 0.72 | 1.00 |
| `DPB1*02:02` | 5.07 | DVT-positive | 159 | 0.98 [0.54–1.77] | 0.94 | 1.00 |

The reported ORs fall outside every observed confidence interval. This is a non-replication with power, not an underpowered null.

---

## 4. Residue level — DQB1 position 57

`AA_DQB1_57_V`, six independent readings:

| Cohort | Model | Test | OR | P |
|---|---|---|---|---|
| narrow | fixed | conditional, signal 1 | **0.586** | 2.63×10⁻⁵ |
| narrow | random | conditional, signal 1 | 0.590 | 4.87×10⁻⁵ |
| intermediate | fixed | conditional, signal 1 | 0.621 | 1.14×10⁻⁴ |
| intermediate | random | conditional, signal 1 | 0.625 | 1.90×10⁻⁴ |
| full | fixed | conditional, signal 1 | 0.624 | 8.73×10⁻⁵ |
| full | random | conditional, signal 1 | 0.628 | 1.81×10⁻⁴ |

Single-marker estimates with confidence intervals (fixed model, `--ci 0.95`):

| Cohort | OR [95 % CI] | P |
|---|---|---|
| narrow | 0.591 [0.460–0.758] | 3.45×10⁻⁵ |
| intermediate | 0.626 [0.491–0.797] | 1.44×10⁻⁴ |
| full | 0.627 [0.495–0.794] | 1.04×10⁻⁴ |

Amino-acid omnibus at DQB1:57 (4 residues, 3 df) — tests the position, not one residue:

| Cohort | deviance | P |
|---|---|---|
| narrow (n = 2,195) | 24.37 | **2.09×10⁻⁵** |
| intermediate (n = 2,508) | 21.06 | 1.02×10⁻⁴ |
| full (n = 3,101) | 22.26 | 5.76×10⁻⁵ |

Thresholds (full): residue 4.35×10⁻⁴, omnibus 6.33×10⁻⁴. All six readings clear both. None reaches 5×10⁻⁸.

**Marker QC** — `determined_rate` 1.00, `frac_unconfirmed_allele` 0.014, MAC 907. Clean.

**After conditioning on DQB1:57**, the other MHC candidates fall away: `AA_B_114_N/D` and `AA_DQA1_175_E` are LD-dependent, not independent. DQB1:57 is the single independent MHC signal.

### Why the position matters

DQB1 position 57 sits at the end of pocket 9 of the DQ peptide-binding groove and is the most mechanistically characterised position in the HLA system — Asp57 is the classic type 1 diabetes determinant. The signal therefore lands on a physical property of the molecule, not a haplotype label.

### Independent convergence

The genome-wide SNP scan's strongest MHC variant, chr6:32,632,861 (P = 7.3×10⁻⁷, OR 0.67), lies inside `DQA1`, adjacent to `DQB1`. The two analyses were run independently and converge on the same sub-region.

---

## 5. Why the residue layer was necessary

The allele scan is null; the residue scan is not. This is not a threshold artefact — it is structural.

DQB1 57V is carried by several distinct alleles (`DQB1*06:04`, `DQB1*05:01`, others). Each individually carries too little of the effect to reach significance:

| Allele | 57 | Case freq | Control freq | case/control |
|---|---|---|---|---|
| `DQB1*06:04` | **V** | 0.0547 | 0.0736 | 0.74 |
| `DQB1*05:01` | **V** | 0.0467 | 0.0697 | 0.67 |
| these two combined | **V** | **0.1014** | **0.1433** | **0.71** |

Residue assignment is taken from samples homozygous for the allele, so only these two are confirmed here; the `AA_DQB1_57_V` marker frequency is 0.1462, the small remainder carried by rarer 57V alleles with no homozygote in the cohort.

Neither allele reaches the allele-level threshold (4.07×10⁻⁴) on its own. Pooling by residue recovers the signal (OR 0.63, P = 1.0×10⁻⁴). **An allele-level analysis of this dataset would have concluded "no HLA association in CTEPH" and stopped.**

---

## 6. Typing artefacts detected and excluded

Two markers passed the significance threshold and are **not** results. Both were caught by checks built into the component.

### 6.1 `AA_A_178_T` — P = 8.0×10⁻⁹, the most significant residue in the table

| Property | Value |
|---|---|
| Driver allele | `A*02:783` |
| Occurrences in jMorp 61KJPN (122,848 chromosomes) | **0** |
| P-group assignment in `hla_nom_p.txt` | **singleton** — not grouped with `A*02:01P` over the antigen recognition domain |
| Frequency in our cohort | 20 / 6,194 = 0.32 % |
| 95 % upper bound from 122,848 Japanese chromosomes | 0.0024 % |
| Ratio | **130×** |
| Platform | DNBSeq-T7 specific (15/327 vs 0/125, Fisher P = 0.015) |
| `frac_unconfirmed_allele` | **0.85** |

`A*02:783` differs from `A*02:01` at exactly positions 175 and 178 — the two positions that produced signal. Regional variation within Japan is on the order of 1.5–2×; it cannot produce 130×.

Excluded from conditional analysis by `--max-frac-unconfirmed 0.5`. **Still present in `sumstats.tsv` as the top row — must be flagged explicitly in the manuscript, not omitted.**

### 6.2 `HLA_F*01` — P = 3.6×10⁻⁵, OR 3.00, top allele hit

| Property | Value |
|---|---|
| Resolution | **1 field** — an unresolved call, not an allele |
| `in_reference` | 0 |
| Case chromosomes from DNBSeq-T7 | **14 of 16** |
| Case frequency | 1.82 % (T7 subset 2.19 %, non-T7 cases 0.83 %) |
| Control frequency | 0.26 % |

Global 1-field call rate is depth-dependent (NovaSeq 30× 0.42 %, DNBSeq-T7 30× 0.62 %, HiSeqX 15× 1.03 %, DNBseq-G400RS 15× 2.78 %), confirming that unresolved calls track technical conditions.

**HLA-F typing is unreliable as a whole**, independent of this marker: jMorp reports HLA-F as effectively monomorphic in Japanese (`F*01:01P` 99.9935 %, `F*01:04P` 0.0065 % — two alleles in 122,848 chromosomes), whereas our typing calls 21 distinct 2-field F alleles totalling 3.6 % non-`F*01:01`. HLA-F is a non-classical class Ib gene; HLA-HD performance there is far below its performance at A/B/C/DRB1.

Correctly excluded by the `in_reference` filter — the conditional analysis returned zero allele signals.

---

## 7. Reference-panel validation

Two Japanese panels, with different strengths. Neither alone is sufficient.

| Panel | Genes | Chromosomes | Resolution | SE at f ≈ 0.10 |
|---|---|---|---|---|
| jMorp 61KJPN (Tohoku) | 13 | 122,848 | P group | **±0.0009** |
| 1000 Genomes JPT (Tokyo) | 5 (A, B, C, DQB1, DRB1) | 186 | 2 field | **±0.022** |

**Existence questions** (does this allele exist in Japanese?) — jMorp alone is decisive; JPT at n = 186 has no power. This is what condemns `A*02:783`.

**Frequency questions** (is this frequency plausible?) — both are required, because HLA frequencies vary within Japan and jMorp samples only the northeast.

### Applied to DQB1 57V

```
                    CASE      CTRL      jMorp(Tohoku)   JPT(Tokyo)
DQB1*06:04  (57V)   0.0547    0.0736       0.0549         0.0806
DQB1*05:01  (57V)   0.0467    0.0697       0.0566         0.0645
─────────────────────────────────────────────────────────────────
57V total           0.1014    0.1433       0.1115         0.1451
```

| Comparison | Difference | z |
|---|---|---|
| Controls vs jMorp (Tohoku) | +0.032 | **+6.5** |
| Controls vs JPT (Tokyo) | −0.002 | −0.07 |
| Cases vs jMorp | −0.010 | −1.0 (n.s.) |
| Cases vs JPT | −0.044 | −1.6 (n.s.) |

**A genuine regional gradient exists at this position (Tohoku 0.1115 → Tokyo 0.1451).** Our controls sit at the Tokyo end, consistent with broad mainland recruitment; our cases fall below the entire reference range. Neither case comparison is individually significant — 878 case chromosomes give SE ≈ 0.010 — so only the direct case-control contrast has power.

### Why this is not residual stratification

The signal is **strongest in the most ancestrally homogeneous cohort**, which is the opposite of what stratification predicts:

| Cohort | n | omnibus P |
|---|---|---|
| narrow (most homogeneous) | 2,195 | **2.09×10⁻⁵** |
| intermediate | 2,508 | 1.02×10⁻⁴ |
| full | 3,101 | 5.76×10⁻⁵ |

Tightening the ancestry window strengthens rather than weakens the signal, despite losing 906 samples. This test is only available because PopGMM defines nested cohorts.

---

## 8. Open issues before submission

| # | Issue | Why it matters | Effort |
|---|---|---|---|
| 1 | **λ_GC is 1.47–2.39** (residue layer 2.37) | A reviewer stops here. λ_GC assumes independent, mostly-null tests; Meff is 115 of 705 (≈6× redundancy) and the MHC is one LD block, so λ_GC measures signal density, not confounding. Stating this is not enough — recompute λ excluding markers in LD with DQB1:57, or give a permutation null. | 1 day |
| 2 | **`AA_B_114_N` and `AA_B_114_D` are the same position** (freq 0.609 + 0.389 = 0.998) | The "5 significant residues" in `full_mainland` fixed are really 3 distinct positions plus 1 artefact. Worse, the random-model conditional lists B_114_N as signal 2 and B_114_D as signal 3 — near-tautological. Conditional analysis should condition on the position, not one residue. Same issue at DQB1:185 (I/T, identical P). | 1 day |
| 3 | **Platform confound not yet formally excluded for DQB1:57** | Cases are DNBSeq/NovaSeq 30×, controls HiSeqX 15× — fully confounded with phenotype. The planned leave-one-platform-out must cover this marker. | planned |
| 4 | **5 of 305 omnibus positions failed to fit** (all DRB1: 10, 30, 31, 38, 40) | Likely separation or collinearity. Needs a stated reason in Methods, not a silent NA. | 0.5 day |
| 5 | **No replication cohort** | CTEPH incidence ≈ 1.9 per 100,000 per year; an independent Japanese WGS case series is not obtainable. Must be stated as a limitation, with the three-cohort consistency and cross-method convergence offered in its place. | — |

Items 1, 2 and 4 do not change any conclusion; they change whether a reviewer believes the conclusions.

---

## 9. What to report

> No HLA allele is associated with CTEPH in this cohort, including `B*52:01` and `DPB1*02:02`, previously reported in Japanese patients and here excluded at ~100 % power for the reported effect sizes. Decomposing the signal to amino-acid residues identifies DQB1 position 57 — pocket 9 of the DQ peptide-binding groove — as consistently associated across three nested cohorts and two model classes, with valine protective (OR 0.63, 95 % CI 0.50–0.79). Two markers that passed the significance threshold (P = 8.0×10⁻⁹ and 3.6×10⁻⁵) were shown to be typing artefacts: both were specific to a single sequencing platform and absent from a 61,424-person Japanese reference panel.

The artefact paragraph belongs in Results, not in a supplementary note. For a cohort in which platform is fully confounded with phenotype, demonstrating that the framework detects and removes its own false positives is a load-bearing claim, not a caveat.

---

## 10. Files

```
analysis/assoc_hla/results/
  <cohort>/00.prep/marker_qc.tsv            per-marker QC, in_reference, frac_unconfirmed_allele
  <cohort>/01.assoc/<model>/<class>/sumstats.tsv
  <cohort>/02.omnibus/omnibus.tsv           per-position LRT
  <cohort>/03.signals/significance.*.tsv    thresholds, Meff, lambda_GC
  <cohort>/03.signals/cond.*.signals.tsv    independent signals after conditioning
  <cohort>/04.prior/prior_screen.*.tsv      prior Japanese reports, DVT-stratified
  <cohort>/figures/{01.scan,02.omnibus,03.qc}/
  _comparison/{tables,figures}/             model and cohort comparisons
```

Reference panels: `../../../jMorp_HLA_types/` (61KJPN + `hla_nom_p.txt`, pinned with MD5SUM), `../../../1KG_HLA_types/`.
