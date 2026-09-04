# Methods — `assoc_gene`

Gene-based rare-variant association for CTEPH, three nested cohorts, two impact strata, two statistics
(CMC burden and SKAT-O) under one fixed-effect logistic model, with a robustness tier per gene and a
publication summary table for the robust genes.

---

## 1. Samples and genotypes

Genotypes are the `tuning.rv` filtered callset for each cohort
(`02.callset_filter/filtered.{bed,bim,fam}`), i.e. the `12_model_inputs` fixed-model MAF < threshold
set after `tuning.rv` sample and variant QC. Only variants with minor-allele count ≥ the cohort's
`tuning.rv` recommendation (`05.qc_collect/minac_recommendation.tsv`; minAC = 2 in all three cohorts)
enter the map.

| cohort | N | cases | controls |
|---|---|---|---|
| `narrow_mainland` | 2,165 | 418 | 1,747 |
| `intermediate_mainland` | 2,477 | 428 | 2,049 |
| `full_mainland` | 3,068 | 438 | 2,630 |

The cohorts are **strictly nested** (narrow ⊂ intermediate ⊂ full, zero exceptions) and differ almost
entirely in controls: the widest has 883 more controls and only 20 more cases than the narrowest.
Agreement across them therefore shows that a signal **survives adding controls**, and is a
**sample-selection sensitivity analysis, not replication**. Every table and figure sidecar in this
component says so.

The phenotype is normalised in `PREP_PHENO_COV` and handed to rvtest as **1/2** (1 = control,
2 = case). rvtest reads a phenotype value of 0 as *missing* (`DataLoader.cpp:1248`), so a 0/1 column
silently drops every control and runs to completion on cases only; the RVTEST task re-reads rvtest's
own "Loaded *n* cases, *m* controls" line and aborts unless both equal `pheno_coding.json`.

## 2. Gene membership: the snpEff canonical map

Variant → gene assignment is taken from the JHRP snpEff index (`params.SnpEffIndex`), joined on the
exact variant ID `chr<N>:<POS>:<REF>:<ALT>` (the callset `.bim` IDs are strictly of this form; verified
on all 647,142 tested IDs and asserted on every row by `emit_test_inputs.py` and `verify.sh` §20). The
assignment rule is **per-gene most-severe**: one row per (variant, gene) carrying that pair's most
severe consequence, so a variant that is HIGH in one gene and MODERATE in another is tested in both.

The canonical unit is (gene symbol, chromosome). A symbol that occurs on more than one chromosome
(small-RNA families) is named `<symbol>@chr<N>` on every occurrence, because rvtest keys its set table
on the name and would silently merge the two into one cross-chromosome set.

### 2.1 Two nested variant sets

| stratum | snpEff impact classes |
|---|---|
| `moderate_high` | HIGH, MODERATE |
| `low_moderate_high` | HIGH, MODERATE, LOW |

`moderate_high` ⊂ `low_moderate_high` by construction. MODIFIER is never tested: it is what snpEff
gives intergenic and deep-intronic variants, which a gene-based test has no principled way to assign.

### 2.2 `MinNumVar` and the denominator, fixed before any test

A gene enters a stratum's family only if ≥ `MinNumVar` = 2 mapped variants survive in that cohort's
callset. `merge_gene_map.py` applies this filter **first** and writes `n_genes_mapped` and
`bonferroni_threshold = α / n_genes_mapped` to `denominator.<stratum>.tsv` at barrier B1, before a single
association has run. Genes below `MinNumVar` are absent from the map, from the set files and from the
denominator alike. Nothing downstream recomputes the threshold: `scan.py` copies the string verbatim
and `gene_significance.py` derives it again only to assert that both files agree.

| cohort | `moderate_high` genes | threshold | `low_moderate_high` genes | threshold |
|---|---|---|---|---|
| `narrow_mainland` | 10,832 | 4.62e-06 | 14,557 | 3.43e-06 |
| `intermediate_mainland` | 11,436 | 4.37e-06 | 15,009 | 3.33e-06 |
| `full_mainland` | 12,378 | 4.04e-06 | 15,691 | 3.19e-06 |

Because the minAC floor is applied per cohort and the cohorts differ in size, the gene families are
**not** nested across cohorts (a gene can have two qualifying variants in `narrow_mainland` and one in
`full_mainland`, and the wider cohorts map more genes overall). The evidence
table keeps `not_in_map` as a state of its own for that reason.

### 2.3 The map is a file: `--setFile` and the span-overlap over-count

rvtest is driven with `--setFile`, one line per gene, one **1-bp range per variant**
(`<chrom>:<pos>-<pos>`), never `--geneFile`: a refFlat span would pull in every variant between the
first and last mapped position, including the ones the annotation filter excluded. The set file is
written from the same in-memory dictionary as the per-gene counts, over the variant IDs read back from
the `.bim` plink2 actually wrote.

A 1-bp range is still not exact. A tabix interval matches a record over `[POS, POS + len(REF) − 1]`, so
a multi-base-REF record (a deletion) is matched a second time by the range of a later mapped variant
inside its span, and rvtest counts it twice. `scan.py` predicts this over-count **exactly** — the number
of records with `len(REF) > 1` whose span strictly contains another mapped variant's position — and
classifies such genes `explained_by_span_overlap`; the map identity `n_var_map == NumVar` stays an
equality. This is an analytic fact as well as bookkeeping: rvtest genuinely enters the deletion twice.
It is harmless for CMC (a binary collapse) but the SKAT-O kernel double-weights it; the flag is carried
in `map_identity.tsv` so the affected genes can be seen.

## 3. The model

Logistic regression of case status on the gene statistic with covariates **sex + PC1…PC10** (the
`bbj_mainland` ancestry PCs, `pc<k>_avg`), the covariate set every association component of this project uses. `params.NPcs` and
`params.CovarName` change it; the RVTEST task aborts if any named covariate is absent from
`covar_rvt.tsv`, because rvtest itself would exit 0 and fit an uncovaried model.

Missing genotypes are mean-imputed (`--impute mean`). The CMC collapse casts the genotype to `int`
(`Model.cpp:82`), so it is unaffected; the SKAT-O kernel needs a dosage and uses the imputed value.
rvtest accepts only `mean`, `hwe` and `drop`; an unrecognised value leaves the consolidator strategy
uninitialised and still exits 0, which is why the RVTEST task also checks its own log for the strategy
line.

Cases and controls were sequenced on different platforms with zero overlap (all cases DNBSeq or NovaSeq,
all controls HiSeq X 15×). No covariate can remove that confound, and it is present in every p-value
this component produces. **No gene here is established as a CTEPH gene; a robust gene is a candidate
that survived the internal checks below.**

## 4. Statistics

### 4.1 CMC burden — score p and Wald effect size

`--burden cmc,cmcWald`. Per gene, the sample is collapsed to carrier / non-carrier of any minor allele
across the set (`getFlippedToMinorPolymorphicGenotype`). **The published p is the CMC score-test p.**
The odds ratio and 95 % CI come from `cmcWald`, the second model of the same invocation on the same
collapsed genotype and covariates: OR = exp(β), CI = exp(β ± 1.96 SE). The Wald p is not published —
two p-values for one model on one row invites the reader to choose — and `2Φ(−|β/SE|)` will not
reproduce the score p. On STBD1 in `full_mainland` / `moderate_high` this gives OR 9.60 (3.57–25.85)
against p = 7.2e-08, and reproduces the statsmodels refit used elsewhere in the project (9.601) to three
figures.

cmcWald writes `1 + n_covar` rows per gene with no coefficient labels; the burden row is the **first**
of each block, and `merge_rvtest.py` asserts the block height against `--n-covar` before taking it.

A gene with no minor allele in one group gives a Wald fit rvtest resets to NA; the table shows `—`
and `tiers.tsv` carries `cmc_or_na`.

### 4.2 SKAT-O

`--kernel skato`, Beta(1, 25) weights, ρ searched on rvtest's fixed grid 0, 0.1, …, 1. The selected ρ
is kept (`rho` column) and published beside the p: ρ = 0 is pure SKAT (variants may act in opposite
directions), ρ = 1 is a weighted burden.

### 4.3 Significance

Per (cohort, stratum, method), two rules over the same family of `n_genes_mapped` genes:

* **called** — Benjamini–Hochberg q < 0.05, genes rvtest did not return entering at p = 1;
* **Bonferroni** — p < α / `n_genes_mapped`, reported as the stricter tier.

Both use the same α and the same n, so BH's smallest critical value *is* α / n and Bonferroni's
rejection set is always contained in BH's; `scan.py` asserts the containment. No further correction is
applied across the two methods or the three cohorts: they test the same genes on largely the same
samples and are not independent hypotheses. The BH cut is a property of each method's own p-value
distribution and is never compared across methods (`comparable_across_methods = 0` in the significance
tables).

### 4.4 Calibration: λ<sub>GC</sub> on a discrete statistic

λ<sub>GC</sub> is reported per method, overall and split at `n_var = 2` / `≥ 3`, and its value is
**reported, not targeted**. The statistic is a count over a handful of minor alleles; for two-variant
genes the normal approximation is anti-conservative near p = 0.5 (a pile-up just under the median) and
fine in the tail, which inflates a median-based λ for those genes only. On the reference data (six
cohort × stratum cells) CMC gives λ = 1.33–1.50 at `n_var = 2` against 1.09–1.14 at `≥ 3`, SKAT-O
1.17–1.37 against 1.16–1.23; a null simulation on the real MAC distribution gives λ ≈ 0.58 for an
exact test, and Fisher's exact test on the same counts gives ≈ 0.6. The residual 1.1–1.2 at `≥ 3` is
reported as measured; the platform confound (§3) is the obvious candidate and cannot be separated
from it here.
λ is computed through `inv_cdf(p/2)` rather than `inv_cdf(1 − p/2)`, which is exactly 1.0 in binary64
below p ≈ 2.2e-16. The figures draw the overall λ<sub>GC</sub> only; the split by set size is in
`scan_qc.tsv` and in the figure READMEs, where it can be read beside the number it explains.

## 5. Per-gene counts

All counts are computed in `EMIT_TEST_INPUTS` from the same dictionary that wrote the set file, on
exactly the tested samples.

* **Minor allele** — decided **once, on the pooled counts**, and applied to both groups, as rvtest does
  (`CMCWaldTest::fit` collapses the flipped-to-minor genotype). 192 of 75,787 variants in
  `full_mainland` / `moderate_high` are ALT-major and are flipped; `n_flipped` is the audit trail.
* **Cumulative MAC** — sum of minor-allele counts over the set, per group, from `plink2 --freq counts`
  on the case and on the control keep list (`--nonfounders`, so the counted sample set is the tested
  one by construction).
* **Carriers** — samples with ≥ 1 minor allele at any site of the set, per group, from a direct decode of
  `tested.bed`. The decode is **cross-checked per variant against the `.acount`** (ALT count and
  observed-allele count must both agree, or the task aborts); 151,574 checks pass on
  `full_mainland` / `moderate_high`.
* **MAF** — MAC / (2 · N<sub>group</sub>), a fixed denominator, so case and control MAFs are on the same
  footing; using `OBS_CT` instead would differ by the ~2 % genotype missingness.
* **Multi-site carriers** — samples with a minor allele at **≥ 2 distinct sites** of one gene. CMC counts
  such a sample once; the MAC counts every site; the two columns of the summary table can therefore
  diverge, and where they do the `multi_site_carriers` column says by how many samples. The bound the
  pipeline asserts is `carriers ≤ MAC ≤ 2 · n_var · carriers`. In `full_mainland` / `moderate_high`,
  19,300 (gene, sample) pairs in 2,453 genes are multi-site; GTF3C3 has 10 case carriers holding 21
  minor alleles. Every such sample is listed in `multi_carriers.tsv`, and the tiered genes' rows are
  collected in `_comparison/robust_genes/multi_carrier_detail.tsv`.

## 6. Robustness tiers

`robust_genes.tier_of(called, n_cohorts)` is a pure function over the BH calls of one (gene, stratum):

```
any[c]  = called[c].cmc or called[c].skato          (c over the cohorts where the gene is mapped)
Tier 1  = any in EVERY cohort  and  both methods called in ≥ 1 cohort
Tier 2  = any in EVERY cohort
Tier 3  = any in ≥ 2 cohorts
```

Bonferroni never enters the rule; it is reported as `called_bonferroni_anywhere`. A gene not in a
cohort's map cannot count as called there, so Tier 1 and 2 require presence in every map. The summary
table and the per-gene figures carry Tier 1 and 2 genes; Tier 3 appears in `tiers.tsv` and in the
robust-genes figure's evidence matrix.

The rule introduces no new threshold and no new p-value. It asks one further question of a call the scan
already made: does it survive a change of statistic and the addition of controls.

## 7. Scope and limits

* Not replication (nested cohorts, §1). Not free of the platform confound (§3).
* A gene-set p-value says nothing about which variant carries the signal; the per-gene figures and
  `variant_detail.tsv` are the only variant-level statements made, and they are descriptive.
* The two strata are nested, so a gene tiered in both is one finding tested on overlapping sets.
* λ<sub>GC</sub> over ~10<sup>4</sup> non-independent, discrete gene tests is a coarse instrument.

## 8. Software

rvtests 20190205 (`--setFile`, `--burden cmc,cmcWald`, `--kernel skato`, `--impute mean`), plink2
alpha 6, htslib 1.9 tabix, bcftools, snpEff JHRP v6 index, Nextflow DSL2 on SLURM, Python 3 with
numpy / pandas / matplotlib (conda env `cteph_geno_pro`).
