# assoc_hla — method and rationale

## 1. What a result from this component means

An odds ratio here is **per copy of the HLA allele or amino-acid residue**, from a logistic
model with `SEX` and 10 `bbj_mainland` principal components. It is an association in a
Japanese case-control series of chronic thromboembolic pulmonary hypertension, not a causal
claim, and the MHC's linkage disequilibrium means a significant marker is usually one of many
on the same haplotype — which is what §5 and §6 exist to disentangle.

## 2. The framework, and why this one

HLA association has a settled procedure, and this component follows it rather than inventing
one. The reference is *Tutorial: a statistical genetics guide to identifying HLA alleles
driving complex disease*, Nat. Protoc. 2023
([doi:10.1038/s41596-023-00853-4](https://doi.org/10.1038/s41596-023-00853-4)), which
codifies the SNP2HLA / Okada / Raychaudhuri approach and names **PLINK and SAIGE** as the
association engines. The Japanese precedent is Okada et al., Nat. Genet. 2015
([doi:10.1038/ng.3310](https://doi.org/10.1038/ng.3310)) and Hirata et al., Nat. Genet. 2019
([doi:10.1038/s41588-018-0336-0](https://doi.org/10.1038/s41588-018-0336-0)), the latter being
the MHC-wide fine-mapping of the Japanese population this study's cohort belongs to.

**One departure from the protocol, deliberately.** It recommends excluding MAF < 1 % because
imputation is noisy at low frequency. This study does not impute — HLA-HD types the alleles
directly from the reads — so that reason does not apply. See §4 for what replaces it.

## 3. HLA dosages as chr6 markers

`hla.typing` emits two matrices of exact integer counts: allele dosage (how many of a sample's
two chromosomes carry `A*24:02`) and residue dosage (how many carry phenylalanine at HLA-A
position 9). `build_hla_markers.py` writes them as biallelic markers, `REF = A` (absent) and
`ALT = P` (present) — the SNP2HLA convention — on real GRCh38 chr6 coordinates taken from
`EnsDb.Hsapiens.v86` at build time.

Three things follow, and all three are the point:

- **Both engines read them unchanged.** plink2 `--glm` and SAIGE step 2 take PLINK genotypes;
  nothing about them knows these are HLA calls. The comparison between the two models is
  therefore a comparison of *estimators*, with the data held identical.
- **The whole downstream is inherited.** Everything speaks the canonical 16-column plink2
  `--glm` schema the sibling components already use.
- **`--LOCO=TRUE` becomes meaningful.** Because the markers are genuinely on chr6, SAIGE holds
  chr6 out of the null model it tests them against. Without that, the GRM — built from
  genome-wide pruned markers, chr6 among them — would partly contain the signal it exists to
  control for. This is the proximal-contamination problem, and placing the markers on their
  real chromosome solves it for free.

**Positions inside a gene are a plotting convention, not a genomic claim.** Each marker sits at
its gene's start plus its index within the gene. That puts it in the right gene, in a stable
order, on the real MHC axis, which is what the figures need. It does not assert that the allele
is caused by a variant at that base.

**The encoding is verified, not assumed.** `plink2 --export A` reproduces the input dosage
matrix exactly — 10,685,260 genotypes, zero mismatches — and `plink2 --freq` confirms
`ALT_FREQS` equals the frequency of `P` at every marker, so the effect allele is carriage.
`hla_to_sumstats.py` re-asserts `A1 == 'P'` on every run, because an inverted effect allele
would silently invert every odds ratio in the analysis.

## 4. Which markers are tested

### DRB3, DRB4 and DRB5 are held out of the primary analysis

Three independent reasons, none of them a judgement call:

1. **Their dosages are wrong.** `hla.typing` OPEN_QUESTIONS §1: HLA-HD writes `-` for a second
   allele it did not call, and `collect_alleles.py` fills it with the first. At the classical
   loci that is correct — `-` means homozygous. At DRB3/4/5 it means *hemizygous*, because
   these are haplotype-dependent paralogues genuinely absent from many chromosomes. A sample
   with one DRB4 copy gets dosage 2. This affects **305 of 2,544 allele columns and 635 of
   3,328 residue columns**.
2. **DRB3 fails its external check.** OPEN_QUESTIONS §2: against jMorp's 61,424 Japanese, DRB3
   reconciles at *r* ≈ 0.64 while eleven other loci sit at 0.986–1.000, and no mechanism tested
   explains it.
3. **DRB3 and DRB4 have no GRCh38 primary-assembly locus.** They exist only on ALT haplotypes,
   so any coordinate assigned to them would be invented. `hla_gene_coords.R` drops them on its
   own; the exclusion is not something the analyst has to remember.

They run as a separately labelled sensitivity arm (`--DropGenes ''`), never merged with the
primary result.

Genes with dosages but no primary-assembly locus at all — the pseudogenes H, W, Y, T, J, K, L,
V, DRB2 and DRB6–9, 64 allele columns between them — are dropped for reason 3 alone. None
carries residue data.

### The frequency floor is a count

**MAC ≥ 20 is primary; MAF ≥ 1 % is the sensitivity arm.** With 419–439 cases, MAF 1 % is about
44 chromosomes in the narrow cohort and roughly 8 on the case side — too thin for a stable odds
ratio, and the published 1 % floor was chosen against a different failure mode (imputation
noise) that does not exist here. Both arms are run and both are reported.

Call rate ≥ 0.95 as well, which drops nothing in practice: HLA-HD's call rate on this cohort is
0.9995 at worst.

### Residue determination

IMGT's protein alignment writes `*` for an unsequenced residue and `.` for a gap, so at some
positions a chromosome contributes no residue at all and the position's dosages sum to 1 or 0
rather than 2. Over the 1,115 positions at the 16 primary genes the median undetermined rate is
0.03 %, but **157 positions exceed 5 %** and DQA1:56 reaches 24.8 %.

**An undetermined residue is not an absent one, and treating it as one manufactures
associations.** A chromosome with `*` or `.` at a position contributes 0 to every residue
there, which in the dosage matrix is indistinguishable from carrying a different residue. That
is wrong in the worst available direction: typing quality in this cohort is *differential*
between cases and controls (`hla.typing` OPEN_QUESTIONS §3), so any residue whose complement is
partly "unknown" absorbs the typing artefact as an association.

This was not hypothetical. On the first run, before the fix, the strongest residue signal in
the entire study was `AA_A_175_G` at *P* = 2.5 × 10⁻⁸. HLA-A position 175 is 99.67 % glycine;
of the 21 samples not homozygous for it, **five — 24 % — were undetermined rather than carrying
alanine or arginine**. The signal rested on twenty chromosomes, a quarter of which were unknown.

`build_hla_markers.py` therefore sets a sample MISSING at **every** residue of a position whose
dosages do not sum to 2, because at least one of its chromosomes has no determined residue. It
is not that the sample is known to lack those residues; it is that we do not know. On the full
cohort this masks 47,439 (sample, position) pairs, and `AA_A_175_G` falls to MAC 16 and out of
the scan through the ordinary frequency filter. `determined_rate` is recorded **before** the
mask, since afterwards every retained sample sums to 2 by construction and the column would
report 1.0 everywhere.

Single-residue tests are otherwise unaffected — a residue either is or is not on a determined
chromosome. The **omnibus additionally** requires the position to be well determined: its
*m*−1 df assumes the position is observed, so positions below `params.MinDetermined` (0.95) are
excluded from it and flagged while their residues are still tested individually.

### A residue marker inherits the allele artefact, so it inherits the flag too

`in_reference` marks an allele that jMorp's 61,424 Japanese have never observed — the strongest
typing-artefact warning this study has. A residue marker carries no allele name, so at first it
had no such flag. **It needs one, because it is carried BY alleles.**

The demonstration is the run itself. After the undetermined-residue fix, the strongest residue
signal in the study was `AA_A_178_T`, *P* = 8 × 10⁻⁹, OR 0.05. Tracing its twenty minor-state
chromosomes back to the alleles that carry them:

| allele | chromosomes | `in_reference` |
|---|---|---|
| **`A*02:783`** | **15** | **0** |
| `A*02:06` | 5 | 1 |
| others | 1–4 each | mostly 1 |

Fifteen of twenty came from **one allele the reference panel has never seen**, in fifteen case
samples against five control ones. That is not a finding about HLA-A position 178; it is
`A*02:783` — an allele that very likely does not exist — wearing a residue's clothes.

`build_hla_markers.py` therefore computes `frac_unconfirmed_allele` for every residue marker:
the share of its **minor-state** chromosomes whose allele is absent from the reference. The
minor state, because an association at low count is made of it — asking what carries T when
99.6 % of chromosomes do says nothing. On the full cohort the median is **0.02** and 31 of 705
tested residues exceed 0.5; `AA_A_178_T` scores **0.85**.

It is a **warning, not a filter**. A real Japanese allele that even a 61,424-person panel never
sampled lands here too. But a significant residue with a high value should not be reported as a
biological result without saying so.

### The variance ratio is a single one, and that is forced

SAIGE calibrates its score test with a variance ratio estimated from randomly chosen markers in
the **step 1 PLINK file** — here the LD-pruned GRM marker set. A per-MAC-category ratio
(`--isCateVarianceRatio=TRUE`) would suit HLA markers, whose frequencies span a much wider range
than that marker set does, and this component was written to use one.

It cannot. Measured on this cohort, the GRM marker set has **minimum MAC 51** and **zero markers
in SAIGE's default low category of (10, 20.5]** — it is a common-variant, MAF-filtered, pruned
set, which is exactly what a GRM should be built from. There is nothing to estimate a low-MAC
ratio from, so requesting one would either fail or return a single category under another name.
A single ratio, as `assoc_saige` uses, is what these markers can support, and the component says
so rather than carrying an option that reads as more careful than it is.

## 5. The four layers

| layer | fixed model | random model |
|---|---|---|
| single allele | plink2 `--glm`, Wald, `--ci 0.95`, `firth-fallback` | SAIGE step 2, `--is_Firth_beta=TRUE`, SPA |
| single residue | same | same |
| **amino-acid omnibus** | *m*−1 df likelihood-ratio test | *(see below)* |
| conditional | `--condition-list` | SAIGE `--condition` |

`firth-fallback` rather than `no-firth`: HLA alleles are often carried by few cases, and a Wald
estimate at the separation boundary is worse than a Firth one. This differs from
`assoc_plink2`, which scans 5.1 M common variants where separation is rare.

### The omnibus is the layer that makes this fine mapping

At an amino-acid position carrying *m* residues, the question is not whether one residue
associates but whether the position's residue **content** does. `hla_omnibus.py` fits the
covariate-only model and a model carrying *m*−1 residue dosages — one residue is the reference
level, because the dosages at a position sum to its determined chromosome count and the full
set is collinear — and compares them by deviance,

> *D* = 2 (log *L*<sub>full</sub> − log *L*<sub>null</sub>) ~ χ²<sub>*m*−1</sub>

which is exactly Hirata et al.'s test. The commonest residue is the reference, which is
conventional and makes the remaining effects readable.

**Neither engine can do this**, because both test one marker at a time, so the fixed-model
omnibus is fitted directly with `statsmodels`. The random model's counterpart is SAIGE-GENE+'s
group test over the same residue set — that is **SKAT-O, a different test, not an *m*−1 df
likelihood ratio** — and it is reported beside the fixed-model omnibus, never merged with it.

Fits that do not converge, or whose deviance is not finite, are reported with `status != 'ok'`
rather than dropped: a position where the model hits the separation boundary is information
about that position.

## 6. Significance

Three thresholds, all computed and all reported, with the effective-test one primary.

**Why not 5 × 10⁻⁸ alone.** It is calibrated for a genome-wide scan of roughly a million
independent common variants. This scan tests about 150 alleles and 970 residues inside 3.4 Mb
of the most extreme linkage disequilibrium in the genome. Applying it here is not conservative;
it is mis-calibrated, because it answers a question about a different experiment.

**Why not plain Bonferroni alone.** 0.05 / 1,100 assumes 1,100 independent tests. The residues
at one position are near-collinear by construction and whole haplotypes travel together across
genes, so that assumption is badly wrong in the other direction.

**What is used.** Li & Ji (2005, *Heredity* 95:221–227) estimate the effective number of
independent tests from the eigenvalues of the marker correlation matrix,

> *M*<sub>eff</sub> = Σ<sub>*i*</sub> [ 𝟙(λ<sub>*i*</sub> ≥ 1) + (λ<sub>*i*</sub> − ⌊λ<sub>*i*</sub>⌋) ]

and the threshold is 0.05 / *M*<sub>eff</sub>. `effective_tests.py` writes all three thresholds,
the number of markers clearing each, *M*<sub>eff</sub>, and λ_GC, so a reader can apply whichever
convention their journal expects without recomputing anything.

**The omnibus gets its own threshold, computed the same way.** There are about 300
position-level tests, not 700 marker-level ones, and they are correlated differently — residues
within a position are near-collinear by construction, which is precisely what the omnibus
collapses. Judging a position against the single-marker threshold would be judging it against
the wrong experiment. A position has no single genotype vector to correlate, so the effective
count uses **each position's first principal component** over its own residue dosages as a
stand-in; that is an approximation and is labelled as one, with the plain Bonferroni count
printed beside it needing no such assumption. Measured on the full cohort: 300 positions,
*M*<sub>eff</sub> = 79, effective threshold 6.3 × 10⁻⁴ against Bonferroni 1.7 × 10⁻⁴ — and
**both select the same two positions**, so the conclusion does not rest on the choice.

λ_GC is **not** reported for the omnibus. It assumes every test carries 1 df; these carry
*m*−1, and *m* varies by position, so the statistic is undefined there and printing it would
look like a calibration check while being nothing of the kind.

**No correction across cohorts or models.** The three cohorts are nested — narrow ⊂ intermediate
⊂ full — and the two models are fitted on the same data. These are not independent experiments
and Bonferroni across them would be badly conservative. `assoc_plink2/docs/METHODS.md` §11 takes
the same position for the same reason; it is inherited here rather than silently reversed.

## 7. Conditional analysis

A single causal allele drags its whole haplotype over the threshold. The number a reader wants
is how many *separate* things are going on, and the accepted way to get it is forward stepwise
conditioning: add the lead marker's dosage to the covariate set, re-scan, repeat until nothing
clears the threshold.

**Conditioning is global across the MHC, not per gene.** LD here spans genes — DRB1, DQA1 and
DQB1 travel together — so conditioning within a gene would leave the signal standing in its
neighbour and count one thing twice.

**The null model is never refitted.** Conditioning adds a covariate at the association step
only, so every round stays comparable with round 0 and with the other rounds. `assoc_saige`
does the same, for the same reason.

## 8. Screening the published Japanese CTEPH claims

CTEPH already has an HLA literature in this population, and a study of this size should say
whether it holds. Two reports:

- **Kominami et al., *J Hum Genet* 2009;54:108–114** — 160 CTEPH patients (99 without deep vein
  thrombosis, 61 with) against 380 controls, typed by HLA-region microsatellites. Reported
  **HLA-B\*52:01 OR 2.47** and **HLA-DPB1\*02:02 OR 5.07**, and — the load-bearing detail —
  **only in the DVT-negative patients**, explicitly not in the DVT-positive ones.
- **Tanabe et al., *Eur Respir J* 2005;25:131–138** — carrier frequencies 40 % vs 24 %
  (B\*52:01) and 19 % vs 6 % (DPB1\*02:02).

`prior_screen.py` tests exactly those alleles, in exactly that stratum. Three points of method:

1. **The stratum is the hypothesis.** Testing the whole case series against controls tests
   something the paper did not claim and would dilute a real DVT-negative effect. DVT status
   comes from the study's own sample sheet via `prep_subgroups.py`, joined on `WGS_ID` — the
   only id column that matches the genotypes, and it matches all of them. A control has no DVT
   status and is marked `not_applicable` rather than blank, so a stratified comparison cannot
   silently drop the control group.
2. **This is not part of the scan and does not compete with it.** Two pre-specified alleles,
   with a direction and a magnitude stated in advance, tested at α = 0.05. The scan pays for
   asking every marker; this asks two.
3. **Power is on every row, always.** A null is only informative if the effect could have been
   seen. Each row carries the power to detect the *reported* odds ratio, the smallest odds ratio
   the stratum could detect at 80 % power, and whether our confidence interval **excludes** the
   published point estimate — which is a stronger and more useful statement than *P* > 0.05.

## 9. What this component cannot do

- It cannot separate a real HLA effect from the typing artefact of `OPEN_QUESTIONS` §3. Controls
  are typed measurably worse than cases for a technical reason that is perfectly confounded with
  phenotype, and a common allele depleted more in controls reads as enrichment in cases. The
  `in_reference` flag marks the alleles most likely to be affected; it does not remove the
  effect.
- It cannot fine-map below the resolution of a 2-field call. Two alleles identical over the
  antigen recognition domain but different elsewhere are two markers here, but the residue layer
  is what actually separates them.
- The three cohorts are **nested**, so agreement between them is not replication. It is a
  statement about robustness to the ancestry boundary, and the figures say so.
