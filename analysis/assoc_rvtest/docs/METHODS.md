# METHODS — assoc_rvtest

Gene-based rare-variant association. This document is about *why* each choice was made; the
output tree and column definitions are in [OUTPUTS.md](OUTPUTS.md), and this study's own
numbers are in [STUDY_NOTES.md](STUDY_NOTES.md).

---

## 1. What a hit from this component means

The design is fully confounded. Cases were sequenced at 30× on DNBSeq/NovaSeq, controls at
15× on HiSeqX, with **zero platform overlap**. A raw rare-variant burden difference between
the groups therefore cannot be separated from a depth and chemistry difference, and no
covariate fixes that, because the confounder and the phenotype are the same variable.

[`tuning.rv`](../../../tuning.rv/README.md) addresses this one layer at a time and returns
the minor-allele-count threshold at which the **adjusted apparent group effect reaches a
statistical null** — the point at which the genome-wide burden artefact is gone. This
component inherits that threshold rather than re-deriving it.

That inheritance is what makes a burden test interpretable here. It does not make any
individual gene real. **A hit is a candidate.** There is no independent cohort for this
phenotype, so nothing in this component is validation.

## 2. Inputs, and why none of the QC is repeated here

Per cohort, from `tuning.rv/results/<cohort>/`:

- `02.callset_filter/filtered.{bed,bim,fam}` — already sample-QC'd (depth-residual outliers
  removed) and variant-QC'd (depth-differential variants above the Kneedle knee excluded).
- `05.qc_collect/minac_recommendation.tsv` — the minAC, and the model it was chosen under.

Re-applying sample or variant QC here would be both wrong and impossible: it is already in
the call set, and the exclusion lists the previous version of this pipeline pointed at no
longer exist. The only filter left is `--mac <minAC>`, applied in `SPLIT_CALLSET`.

`params.MinAC = null` means inherit. The resolved value is read at DAG-build time and written
to `run_manifest.json`, so the run is reconstructable without reading `tuning.rv`'s logs. If
the recommendation file is absent the pipeline **fails at launch naming the path it wanted**;
proceeding with an unfiltered call set would produce a plausible and meaningless result.

### Covariates

`SEX` plus the first 10 `bbj_mainland` PCs, matching the fixed-model components.
`rvtest_prepare_tools.py` lowercases every column when it reformats, so `params.CovarName` is
lowercase (`sex,pc1_avg,…`) while the source header is uppercase (`SEX,PC1_AVG,…`).

**rvtest exits 0 when a `--covar-name` column is absent**, silently fitting an uncovaried
model. That failure mode is invisible in the logs and in the output. `RVTEST` therefore
asserts every named column exists in the prepared covariate file before it runs, and fails
loudly if not. This check is the reason the covariate set can be trusted.

## 3. Impact strata

snpEff `IMPACT` defines three nested strata:

| stratum | includes | genes tested (narrow_mainland) |
| --- | --- | --- |
| `high` | HIGH | 164 |
| `moderate_high` | MODERATE, HIGH | 8,566 |
| `low_moderate_high` | LOW, MODERATE, HIGH | 12,966 |

MODIFIER is excluded from all three — it is the annotation snpEff gives to intergenic and
deep-intronic variants, which a gene-based test has no principled way to assign.

The strata are reported side by side rather than one being chosen. A narrower stratum tests
fewer, more interpretable variants at the cost of power; there is no way to know in advance
which trade is right for a given gene.

## 4. Test methods: SKAT-O and CMC, not Zeggini

One kernel test and one burden test, which is the information the previous three carried.

**Zeggini was dropped on measurement, not on taste.** Against CMC it returned a bit-identical
*P*-value for **63 % / 76 % / 96 %** of genes (`low_moderate_high` / `moderate_high` /
`high`). The two differ only when a sample carries **two or more** rare variants in the same
gene — CMC collapses to a 0/1 indicator, Zeggini sums the count. At this minAC almost no
sample does, so the two statistics coincide. Running both cost a third of the association
time and added no independent evidence.

SKAT-O and CMC are kept because they make genuinely different assumptions about effect
direction: CMC assumes the qualifying variants act the same way, SKAT-O does not. They
disagree in this data (gene-level λ 1.10–1.24 for CMC against 1.18–1.24 for SKAT-O), which
is the point of running both.

## 5. Tiers: Benjamini-Hochberg, per scan

`significant` = FDR < 0.05, `suggestive` = FDR < 0.10, computed over the genes of **one
scan** — one cohort × one stratum × one method.

FDR rather than a fixed *P* because the number of genes tested varies 80-fold between strata.
A fixed cut would impose a different multiple-testing burden on each and make them
incomparable; BH adapts to the number of tests actually performed.

**The price, stated plainly:** the *P* that FDR 0.05 implies is not the same in every
stratum. Measured in `narrow_mainland`:

| stratum | genes | *P* at FDR 0.05 | *P* at FDR 0.10 |
| --- | --- | --- | --- |
| `high` (SKAT-O) | 164 | 2.6 × 10⁻⁴ | 2.6 × 10⁻⁴ |
| `moderate_high` (SKAT-O) | 8,566 | 4.5 × 10⁻⁵ | 1.2 × 10⁻⁴ |
| `low_moderate_high` (SKAT-O) | 12,966 | 1.4 × 10⁻⁵ | 7.5 × 10⁻⁵ |

So a `significant` gene in the HIGH stratum rests on roughly an order of magnitude weaker
evidence than one in LOW+MODERATE+HIGH. Every scan figure draws **its own** BH-implied lines
and prints their values, rather than a shared constant that would imply one rule. **A tier
label is comparable within a stratum, not across strata.**

Genes with fewer than `params.MinNumVar` (3) qualifying variants are removed before BH, and
genes whose test did not converge (a blank or non-numeric *P*) are dropped rather than
counted — a non-result must not dilute every other gene's adjusted *P*.

## 6. Parallelism, and the one barrier it requires

Everything from `SPLIT_CALLSET` to `RVTEST` runs on balanced bins of whole chromosomes.

This is **exact, not approximate**. No gene spans a chromosome boundary, and rvtest tests
each gene independently of every other, so the union of per-bin results is identical to a
whole-genome run. Any grouping of *whole* chromosomes is equally valid; the only question is
how many jobs to pay for.

`chr2` alone carries 9.05 % of variants, so no split can bring the critical path below 9 % of
the serial time. Twelve bins reach exactly that floor while running 216 association jobs
instead of the 396 that one-bin-per-chromosome would need.

**`MERGE_ASSOC` is not optional.** Benjamini-Hochberg is a property of the whole scan, so the
per-bin tables must be concatenated into one before any tier is assigned. It is the single
place where the chromosome split could be silently broken — computing FDR per bin would give
each bin its own, far more lenient, threshold.

Measured effect: the critical path per cohort falls from ~203 min to ~19 min.

## 7. Per-gene follow-up

Every gene reaching a tier gets a deep dive: its variants, their case and control genotype
counts and frequencies, the cumulative allele count the test collapsed, a burden odds ratio
fitted on the same covariates, and each variant's frequency in ToMMo 60KJPN.

This is not optional decoration. A gene-based test reports one *P* for a whole gene, and at
this minAC that *P* can rest on a dozen carriers. *STBD1* — the top gene in every stratum and
method — has a cumulative minor-allele count of 14 (9 case, 5 control) in `narrow_mainland`.
The *P* of 4.1 × 10⁻⁷ is a statement about fourteen alleles, and the follow-up is what makes
that visible.

**The ToMMo comparison is the artefact check that matters under this design.** Cases and
controls were sequenced differently, so the control side is where a technical artefact shows.
A control frequency that departs from the population reference is a genotyping difference,
not a disease association. Agreement with ToMMo is what makes a case excess worth reading.

## 8. Cross-cohort comparison — robustness, not replication

The three cohorts are **nested**: `narrow ⊂ intermediate ⊂ full`. They share their cases
entirely and most of their controls. A gene appearing in all three has therefore **not been
replicated** — a chance hit in the shared core reappears in all three by construction.

What the comparison does measure is robustness to the ancestry filter, because what the wider
cohorts add is precisely the samples the stricter filter rejected. Read instability, not
agreement: a gene called only in `full_mainland` is a candidate artefact of the
ancestry-outlier samples that cohort admits.

`cross_cohort_genes.py` distinguishes `not_a_hit` (tested, reached no tier) from `not_tested`
(never tested in that cohort, too few qualifying variants). Collapsing them would hide
whether a cohort saw anything there at all.

## 9. What would actually validate a hit

Not available in this component, and worth stating so it is not mistaken for something that
is: an independent cohort; or the same gene in an external rare-variant resource for a
related phenotype; or read-level inspection of the carriers' alignments, which is the only
way to rule out that a low-frequency call is a platform artefact rather than a variant.
