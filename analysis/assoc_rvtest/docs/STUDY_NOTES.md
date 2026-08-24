# STUDY_NOTES — assoc_rvtest, CTEPH / AGP3K v6

This study's configuration, its numbers, and the anchors a refactor must not move. The
method itself is in [METHODS.md](METHODS.md).

## Configuration

| | |
| --- | --- |
| `Cohorts` | `narrow_mainland` ⊂ `intermediate_mainland` ⊂ `full_mainland` (nested, narrowest first) |
| call set | `tuning.rv` `02.callset_filter/filtered` — sample- and variant-QC'd |
| `MinAC` | **inherited** from `tuning.rv`: 2 / 2 / 2, `preferred_model = noPC`, `calibrated = true` |
| covariates | `SEX` + `PC1_AVG–PC10_AVG` in the `bbj_mainland` space |
| strata | LOW+MODERATE+HIGH, MODERATE+HIGH, HIGH |
| methods | SKAT-O, CMC (Zeggini dropped — see METHODS §4) |
| tiers | BH FDR < 0.05 significant, < 0.10 suggestive, per scan |
| `MinNumVar` | 3 |
| phenotype | CTEPH case/control |

Cohort sizes after `tuning.rv`'s QC: 2,165 / 2,477 / 3,068 samples and 14.75 M / 16.51 M /
18.95 M variants.

## The confounding this analysis sits on

Cases 30× (DNBSeq/NovaSeq), controls 15× (HiSeqX), **zero platform overlap**. This is not
correctable by covariate — the confounder and the phenotype are the same variable. What
`tuning.rv` establishes is that at minAC 2 the *adjusted* apparent group effect is at a
statistical null (β_group −0.0076 / +0.0023 / −0.0082, *P* 0.32 / 0.75 / 0.19). That removes
the genome-wide burden artefact. It does not make any individual gene real.

**Every result in this component is a candidate.** No independent cohort exists for CTEPH at
this scale, and the three cohorts here are nested, so they cannot replicate each other.

## Sanity anchors — a refactor must not move these

Measured in `narrow_mainland` under the previous (pre-rebuild) configuration, which used the
same call set lineage and covariates:

| gene | stratum | method | *P* | FDR | NumVar |
| --- | --- | --- | --- | --- | --- |
| **STBD1** | MODERATE+HIGH | CMC | 4.10 × 10⁻⁷ | 0.0035 | 3 |
| **STBD1** | MODERATE+HIGH | SKAT-O | 4.74 × 10⁻⁷ | 0.0041 | 3 |
| **STBD1** | LOW+MODERATE+HIGH | SKAT-O | 2.67 × 10⁻⁷ | — | 6 |
| PES1 | MODERATE+HIGH | SKAT-O | 1.13 × 10⁻⁶ | 0.0048 | 4 |

STBD1 is the top gene in **every** stratum × method. Its three MODERATE-impact variants are
exactly `tuning.rv`'s `params.ProtectVariants`, which exist so that variant QC can never
remove them; the pipeline must still find them after the rebuild.

Note STBD1 is absent from the HIGH stratum: its variants are MODERATE, so a HIGH-only scan
has nothing to collapse.

## What STBD1 actually rests on

| | case | control |
| --- | --- | --- |
| cumulative MAC | 9 | 5 |
| cumulative MAF | 0.0110 | 0.0014 |

Burden OR 12.23, 95 % CI [3.80, 39.41]. **Fourteen minor alleles across 2,165 samples.** The
interval spans an order of magnitude, and the point estimate is a description of the carriers
observed rather than a population estimate. This is the reason the per-gene follow-up is part
of the pipeline and not an afterthought.

The three variants' control frequencies agree closely with ToMMo 60KJPN (0.00114 vs 0.00132;
0.00029 vs 0.00041), so the signal is **not** a control-side genotyping artefact — which,
under a design this confounded, is the first thing that had to be ruled out.

## Gene-level calibration

λ_GC over the gene *P*-values, `narrow_mainland`:

| stratum | SKAT-O | CMC |
| --- | --- | --- |
| HIGH | 1.241 | 1.222 |
| MODERATE+HIGH | 1.176 | 1.105 |
| LOW+MODERATE+HIGH | 1.219 | 1.102 |

All above 1, and SKAT-O consistently above CMC. Some of this is real polygenicity, some is
residual confounding that minAC 2 did not remove. It is a reason to read the tiers as
optimistic, and a reason the cross-cohort figure shows λ before it shows any gene.

## Tier thresholds are not comparable across strata

BH adapts to the number of genes tested, so the implied *P* differs (SKAT-O,
`narrow_mainland`):

| stratum | genes tested | *P* at FDR 0.05 |
| --- | --- | --- |
| HIGH | 164 | 2.6 × 10⁻⁴ |
| MODERATE+HIGH | 8,566 | 4.5 × 10⁻⁵ |
| LOW+MODERATE+HIGH | 12,966 | 1.4 × 10⁻⁵ |

A `significant` label in HIGH is roughly an order of magnitude weaker evidence than the same
label in LOW+MODERATE+HIGH. The figures state their own threshold for this reason.

## Zeggini redundancy, measured

Bit-identical *P* to CMC for 96 % / 76 % / 63 % of genes (HIGH / MODERATE+HIGH /
LOW+MODERATE+HIGH). Dropped; see METHODS §4.

## Known limitations, in order of how much they should worry a reader

1. **The design is confounded and no internal analysis fixes it.** minAC 2 removes the
   genome-wide artefact, not the per-gene one.
2. **No replication is possible here.** The cohorts are nested; external validation would
   need a different cohort or an external rare-variant resource.
3. **Effect sizes are uninterpretable at these allele counts.** Report the counts alongside
   any odds ratio, always.
4. **Gene-level λ is 1.10–1.24.** The tiers are optimistic by an amount that has not been
   separated into polygenicity and residual confounding.
5. **Residual platform chemistry is not fixable from within this data set.** Read-level
   inspection of carrier alignments is the outstanding check for any gene taken forward.
