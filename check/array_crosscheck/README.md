# `array_crosscheck` — do the association results depend on the WGS genotypes?

A check, not an analysis of record. Every control is sequenced on HiSeq X and every case on another platform, so a
WGS genotyping difference could look like association, and platform cannot enter the model. SNP arrays give a second
technology for the same people: if a signal were made by WGS calling, the arrays would not reproduce it.

**What it is not.** It is not replication: the people are the same. The arrays are confounded with case status too
(cases on Illumina ASA24v1, controls on JSA24v1), but by a different technology, so an artifact of one would not be
reproduced by the other. Agreement argues against both; it cannot rule out an error shared by both.

## Design

- **People:** unrelated full-mainland genomes (fixed-model set) with array genotypes, 1,683 in all:
  - 299 cases, all DNBSeq-T7 (target 30x) on ASA24v1;
  - 1,384 controls, all HiSeq X (target 15x) on JSA24v1.
  - NovaSeq and DNBSeq-G400RS genomes have no array genotypes, so they are not checked.
- **Array data:** the QC'd set before the ToMMo filter, `cteph_agp3k.v5/array/results/10.variant_qc_filtered`
  (1,750 samples, 409,088 variants: sample call rate, merge, biallelic, variant missingness, autosomes, MAF, HWE).
  The rare drivers come from the per-array sets after sample QC (`02.sample_qc`), because the merge keeps only sites
  on both arrays.
- **WGS data:** `wgs.auto.par/results/12_model_inputs/full_mainland/genotype/fixed_model/…fixed_ready`.
- **Model:** the same on both technologies: logistic regression with sex + 10 BBJ-mainland PCs and Firth fallback
  (plink2 `--glm`). The log check requires "11 covariates loaded".
- **No extra missingness filter** on the shared sites. WGS missingness itself differs by platform (median 3.2% in
  HiSeq X controls vs. 1.3% in DNBSeq-T7 cases; arrays about 0% in both), so filtering on it would select sites by
  platform.

## Steps (SLURM, `--rsc`, queue gr10478b or eb)

| step | script | output |
|---|---|---|
| 1 | `scripts/01_prepare.sh` | `results/01.inputs` (keep lists, shared IDs, both sets), `02.assoc/{wgs,arr}`, `03.concordance` (missingness, HWE, all discordant calls) |
| 2 | `scripts/02_targets.py` | `01.inputs/targets.tsv`, windows (+/-500 kb) |
| 3 | `scripts/03_export.sh` | genotypes near the leads, WGS and array association near the leads, concordance by case status |
| 4 | `scripts/04_compare.py` | `results/04.compare/` (`summary.json`, `signals.tsv`, `signals.direct_on_array.tsv`, `genome_wide.*`, `hla57_tags.tsv`), `results/figures/` |
| 5 | `scripts/05_rare_direct.sh` (+ `.py`) | `results/05.rare_direct/rare_direct.tsv` |
| 6 | `scripts/06_imputed.sh` (+ `.py`) | `results/06.imputed/` (`recovery_by_group.tsv`, `presence_by_carriers.tsv`, `drivers.tsv`, `leads.tsv`, `summary.json`) |

```bash
cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/array_crosscheck
PY=/home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3
R="-t 01:00:00 --rsc p=1:t=8:c=8:m=36568M"
sbatch -p gr10478b $R -o logs/01_prepare.log --wrap "bash scripts/01_prepare.sh"
sbatch -p gr10478b $R -o logs/02_04.log --wrap "$PY -I scripts/02_targets.py && bash scripts/03_export.sh && $PY -I scripts/04_compare.py"
sbatch -p gr10478b $R -o logs/05_rare.log --wrap "bash scripts/05_rare_direct.sh"
sbatch -p gr10478b $R -o logs/06_imputed.log --wrap "bash scripts/06_imputed.sh"
```

Outputs hold counts and estimates only. The plink `.diff` files in `03.concordance` and the `.raw` genotype exports
(`01.inputs`, `05.rare_direct`, `06.imputed`) are keyed by sample ID; they stay internal.

## Results (2026-10-29)

**Genotype concordance** of non-missing calls at the 319,308 shared sites:

| group | concordance |
|---|---|
| cases (DNBSeq-T7 vs. ASA24v1) | 99.9882% |
| controls (HiSeq X vs. JSA24v1) | 99.9837% |

**Genome-wide**, 317,703 shared sites with MAF ≥ 0.01:

| metric | value |
|---|---|
| λ, WGS genotypes | 1.144 |
| λ, array genotypes | 1.143 |
| correlation of z (WGS vs. array) | 0.986 |
| slope of z (array on WGS) | 0.988 |
| sites with WGS P < 10⁻⁴ | 67 |
| … of which, same direction and P < 0.01 on the array | 67 |
| sites with \|Δz\| > 4 | 5 |

The inflation of this fixed-effect model is the same on both technologies, so it is not made by WGS genotyping.
All 5 sites with |Δz| > 4 are array-side calling errors in cases (ASA discordance 3.5–18%), none near a signal.

**Signals** (`signals.tsv`; ORs for the lead's A1, in the 1,683 people):

- Typed directly on the array: the *FTO* lead (WGS OR 1.56, P = 2.7 × 10⁻⁵; array 1.55, P = 3.8 × 10⁻⁵), FGG, F11,
  HLA-DRA, and 11 credible-set members.
- With a proxy of r² ≥ 0.8, where the WGS and array ORs agree to the second decimal: HS3ST3A1, GALNTL6, MHC lead /
  HLA-DQA1, INPP5D, SLC44A2, ABI3BP.
- ABO: the lead rs687289 failed the WGS missingness filter (missing in 5.4% of the target-15x group vs. 1.2% of
  the target-30x group). On the array it gives OR 1.32 for the
  European risk allele A (P = 0.0039 in these 1,683), against 1.35 from WGS in full mainland.
- Low-frequency leads have no array proxy (r² < 0.31) and cannot be checked this way: ARL13B, LINC01146,
  ARHGEF26-AS1, MFSD11, CDH4, PHF20L1, ALPI, ADARB2.

**HLA-DQB1 57 V:**

- Dosage predicted from 1,569 array SNPs (chr6:32.1–33.3 Mb) by ridge regression, trained on the WGS typing of
  controls only. Out-of-fold r² is 0.988 in controls and 0.999 in cases.
- Frequency, WGS typing vs. array prediction: cases 10.4% vs. 10.7%; controls 17.5% vs. 17.8%. The residuals do not
  differ between the groups (P = 0.67).
- Association per copy:
  - WGS typing: OR 0.58 (0.43–0.78), P = 3.3 × 10⁻⁴;
  - array prediction: OR 0.57 (0.42–0.77), P = 2.8 × 10⁻⁴.
- The best single tags show the same pattern against ToMMo, e.g. chr6:32441707 T: cases 16.9%, controls 22.7%,
  ToMMo 18.3% (`hla57_tags.tsv`). The controls' excess over ToMMo is present on the array too, so it is not a WGS
  typing error.

**Rare drivers typed directly** (`05.rare_direct/rare_direct.tsv`; ASA24v1, 299 cases, only 5 controls):

| variant | result |
|---|---|
| *RPH3A* p.Arg544Cys | 10/10 WGS case carriers confirmed; no array-only carrier |
| *OSBP* p.Pro712Ser | 4/4 confirmed; no array-only carrier |
| *JAK2* p.Val617Phe | of 4 WGS case carriers, 2 are called on the array and 2 have no array call (as expected for small clones); 0 carriers in 1,379 controls on both technologies |

The *STBD1* drivers (c.220+1G>C, p.Trp212Cys) and *OSBP* p.Arg462Gln are on neither array.

**Imputed arrays** (step 6; `cteph_agp3k/array_imputed_rev1/results/04.sample_qc`, all variants before the MAF
filter; Minimac4 4.1.4, 2025-10-28). The reference panel is AGP3K whole genomes (the user, 2026-10-30), i.e. the
controls' own sequences. Imputed controls therefore reproduce their WGS calls. The cases' imputed genotypes come from
their array and the panel haplotypes, without their DNBSeq-T7 reads: they are the independent part of this check.

- **The panel holds the controls, not the cases** (rare variants within +/-500 kb of the drivers, `recovery_by_group.tsv`,
  `presence_by_carriers.tsv`):
  - WGS carriers imputed as carriers (dosage >= 0.5): controls 98.8–99.8% (singletons 99.5%); cases 41–83%, rising with
    the number of control carriers;
  - WGS rare variants present in the imputed set: 97% of those carried only by controls, 26% of those carried only by
    cases.
- **Rare drivers** (299 cases, `drivers.tsv`; expected recovery = case recovery of variants with the same number of
  control carriers):

  | variant | WGS case carriers imputed as carriers | expected | cases imputed as carriers without a WGS call |
  |---|---|---|---|
  | *STBD1* c.220+1G>C | 3 of 3 (dosage 1.00 each) | 41% | 0 |
  | *STBD1* p.Trp212Cys | 1 of 2 | 54% | 0 |
  | *RPH3A* p.Arg544Cys | 7 of 10 | 78% | 0 |
  | *OSBP* p.Pro712Ser | 3 of 4 | 73% | 0 |
  | *OSBP* p.Arg462Gln | 0 of 1 | 54% | 0 |
  | *JAK2* p.Val617Phe | not in the imputed set (no carrier in the panel) | — | — |

  The 3 *RPH3A* and 1 *OSBP* carriers that imputation misses are typed directly on the ASA array as carriers (step 5),
  so a miss here is an imputation miss, not evidence against the WGS call. At *STBD1* c.220+1G>C, 67 of the 1,384
  controls have no WGS genotype (2 of 299 cases); one of them is a carrier in the panel. *RPH3A* has one such control
  too. Missing control genotypes can hide carriers (see the missing-genotype share of the gene tests).
- **Common-variant leads** (`leads.tsv`; same model, same 1,683): ORs on imputed dosages agree with WGS, slightly
  smaller, e.g. *ARL13B* 4.42 vs. 4.00, *ARHGEF26-AS1* 3.30 vs. 2.84, *MFSD11* 2.15 vs. 1.82, *ALPI* 1.82 vs. 1.68,
  *LINC01146* 0.60 vs. 0.64 (WGS vs. imputed). Case allele frequencies agree (dosage–genotype r^2 in cases 0.91–1.00).
  Exceptions:
  - *CDH4*: r^2 in cases 0.18 (34 of 299 cases discordant), although the case frequency is similar (5.4% WGS vs.
    4.8% imputed; controls 2.2%) and the OR holds (2.40 vs. 2.16). The lead is not confirmed person by person;
  - *PHF20L1*: r^2 in cases 0.76, OR 1.97 vs. 1.70 (P 0.012 vs. 0.061);
  - *ADARB2* (indel chr10:1747167:CGTG:C): not in the imputed set, cannot be checked.
- In controls r^2 is 1.00 at every lead: that is the panel, not a check.

## Open

- *CDH4* and *PHF20L1* leads: imputation agrees poorly with WGS in cases; *ADARB2* cannot be imputed.
- *STBD1*: Sanger or amplicon sequencing of the carriers, or a read-level review, remains the direct check.
