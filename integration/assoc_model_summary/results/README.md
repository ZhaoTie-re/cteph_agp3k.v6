# Model choice and main results of the CTEPH GWAS — figures and text for slides

Three questions, one folder each. Every folder has the figures and a `README.md` with a **slide summary** (brief) and a **detailed description**; the slide summaries are repeated here for copying.

| Folder | Question | Figures |
|---|---|---|
| `01.fixed_vs_random/` | Fixed vs random model: is a fixed-effect model enough? | `fixed_vs_random.png`, `calibration.png` |
| `02.pc_effect/` | Random model with and without PCs: what do the PCs add? (conclusion pending) | `grm_pc_vs_grm.png` |
| `03.main_results/` | Main analysis in the three cohorts | `main_manhattan.png`, `main_summary.png` |

---

### 1 · Fixed vs random model

**Take-home: a fixed-effect logistic model remains inflated even after adjusting for 10 PCs and removing close relatives; the logistic mixed model with a full GRM is well calibrated. The main analysis therefore uses the random (mixed) model.**

- Fixed model (PLINK2; SEX + 10 PCs; relatives removed): λGC = 1.120 / 1.121 / 1.132 (narrow / intermediate / full mainland).
- Random model (SAIGE; SEX + 10 PCs + full GRM; relatives retained): λGC = 1.039 / 1.038 / 1.021.
- The excess is genome-wide: in every allele-frequency bin (MAF 0.01–0.50) the fixed model has λGC 1.03–1.15, the random model 0.95–1.06.
- Variants at P < 10⁻⁵ (≈ 51 expected by chance): fixed 130 / 216 / 204; random 87 / 160 / 111.
- Genome-wide loci: fixed 0 / 1 / 2, random 0 / 1 / 0; genome-wide only under the fixed model: *MIR3681HG* in full mainland (P 4.3 × 10⁻⁸ → 1.0 × 10⁻⁷); *FTO* in full mainland (P 2.9 × 10⁻⁸ → 1.8 × 10⁻⁷).

---

### 2 · Random model with and without PCs

**Conclusion: pending — to be agreed. The facts:**

- With the full GRM in the model, adding 10 PCs leaves λGC nearly unchanged: 1.040 / 1.024 / 1.017 without → 1.039 / 1.038 / 1.021 with PCs (narrow / intermediate / full mainland).
- The PCs differ between cases and controls in every cohort: likelihood-ratio test of PC1–10 on case status, P = 5.4 × 10⁻¹⁰ / 0.0081 / 6.0 × 10⁻⁴.
- Variants at P < 10⁻⁵ (≈ 51 expected): 180 / 168 / 121 without PCs, 87 / 160 / 111 with PCs.
- Peaks found only without PCs: 8 / 2 / 3, among them in the extended MHC (chr6 28.4–33.5 Mb): *XXbac-BPG308K3.5* (narrow mainland, 1 variant), *DDR1* (narrow mainland, 27 variants), *HCG17* (intermediate mainland, 1 variant); only with PCs: 6 / 3 / 6.
- Effect sizes of the 48 cross-cohort lead × cohort pairs barely move: median |Δ ln OR| = 0.014 (max 0.23).

---

### 3 · Main analysis in the three cohorts

**Take-home: the main analysis (random model: SAIGE, SEX + 10 PCs + full GRM) is well calibrated in all three nested cohorts (λGC 1.039 / 1.038 / 1.021). *FTO* (rs9934504): OR 1.67 / 1.68 / 1.61, P = 1.4 × 10⁻⁷ / 3.8 × 10⁻⁸ / 1.8 × 10⁻⁷; genome-wide in intermediate mainland. 15 further loci are suggestive in all three cohorts, each with the same direction of effect in every cohort.**

- Cohorts (cases / controls): Narrow mainland 419 / 1,776; Intermediate mainland 429 / 2,079; Full mainland 439 / 2,662.
- Loci, genome-wide + suggestive: narrow 0 + 29; intermediate 1 + 27; full 0 + 32.
- Suggestive in all three cohorts: *ARL13B*, *HS3ST3A1*, *LINC01146*, *GALNTL6*, *ARHGEF26-AS1*, *HLA-DQA1*, *INPP5D*, *MFSD11*, *CDH4*, *SLC44A2*, *PHF20L1*, *ABI3BP*, *ALPI*, *OR13C8*, *ADARB2*.
- The cohorts are nested and share every case: agreement between them shows robustness to the ancestry filter, not replication.

---
