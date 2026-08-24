# Methods — `assoc_plink2`

What the component does and why, independent of any dataset. Figures report; this document states the
rules behind them. Everything measured on a particular dataset — counts, λ values, named loci,
configured resources — belongs in [`STUDY_NOTES.md`](STUDY_NOTES.md), not here.

Values written as `params.X` are configuration; the pipeline's `SITE CONFIGURATION` block is the only
place they need to be set.

---

## 1. What this component is, and is not

A **fixed-effects scaffold**. It establishes the scale of the signal, produces candidate peaks and
gives the calibration read-out. It is **not** a definitive association analysis — see §7.

## 2. Symbols and notation

Defined once in `_shared/scripts/plot_style.py` (`SYMBOL_DEFS` / `FORMULAS`) and reused verbatim in every
figure caption, every figure document, and this table.

| symbol | meaning |
|---|---|
| −log₁₀*P* | evidence against the null for one variant under one genetic model |
| λ_GC | genomic-control inflation factor, median χ² / 0.4549; **reported, never applied** |
| OR | odds ratio per copy of the A1 (effect) allele, with its 95 % CI |
| EA / OA | effect allele (plink2's `A1`) and the other allele |
| MAF | minor-allele frequency in the analysed samples |
| *r*² | LD with the peak lead, binned on one scale for all three LD sources |
| EAF | frequency of the **effect** allele among called samples of one group |
| N_eff | 4 / (1/N_case + 1/N_ctrl), the balanced-design equivalent sample size |
| PIP | SuSiE posterior inclusion probability |
| ERRCODE | plink2's per-variant fit diagnostic |
| `called_peak` | whether a variant was a peak **in that cohort** — `genome_wide` / `suggestive` / `not_a_peak` |

## 3. Association model

```
logit Pr(case_i) = β₀ + β g_i + γ_sex SEX_i + Σ_{k=1..K} γ_k PC_{k,i}      K = params.NPcs
```

`g_i` is the genotype under the stated model — additive counts alt alleles (0/1/2), dominant contrasts
carriers against non-carriers, recessive contrasts alt-homozygotes against the rest.

```
plink2 --bfile <cohort genotypes>
       --pheno pheno.tsv --pheno-name PHENO1
       --covar <params.CovarFile> --covar-name SEX,PC1_AVG-PC<K>_AVG
       --glm [dominant|recessive] omit-ref no-firth hide-covar --ci 0.95
```

**Covariates are SEX + the first `params.NPcs` principal components**, and the covariate set is fixed
across cohorts and models so the nine scans are comparable. `params.PcLabel` names the space the PCs
were computed in; `params.CovarLabel` is the human-readable form printed on every figure and document,
so a figure can never describe covariates the run did not fit.

**Firth.** `no-firth` is the default. Firth penalisation would stabilise variants near complete
separation, but it changes the estimand — the penalised estimate is shrunk toward zero by a
data-dependent amount, so a Firth OR and an unpenalised OR are not on the same scale and cannot share
a table. Variants that fail to fit are instead flagged, counted and excluded (§5), which is an
explicit, auditable exclusion rather than a silent change of estimator. `params.FirthMode = ''`
restores plink2's Firth fallback.

## 4. One model defines the peaks

`params.PeakModel` (additive by default) defines every peak. The other models are run and are
*reported* — their Manhattan/QQ figure, their row in `scan_qc.tsv`, and their own genome-wide loci with
full annotation (§8b) — but they define no peaks and get no fine-mapping, conditional analysis or
regional plot.

The alternative — union the significant variants across models and report each peak's lead under the
peak model's statistics — produces hit-list rows whose headline *P* is the peak model's value at a peak
that only another coding flagged, and sends fine-mapping summary statistics that do not contain the
evidence which created the peak. Reporting a model's result under another model's statistics is not a
defensible summary.

The non-primary models remain informative about *calibration*: a strongly deflated λ together with a
large `VIF_INFINITE` count is a direct statement about how few alt-homozygotes the sample carries.
That deflation is why their genome-wide loci (§8b) are reported with an explicit caveat rather than
alongside the primary results.

## 5. Which rows count as a result

Two conditions, both necessary (`usable_mask()` in `_shared/scripts/call_peaks.py`, applied identically to the
QC summary, the lead pick and the per-peak statistics):

1. `ERRCODE == '.'` — plink2 itself reports a clean fit.
2. `OR` and `SE` finite and positive, and `0 < P ≤ 1`.

**The second is not redundant.** plink2 leaves `ERRCODE = '.'` on some degenerate fits: at a near-fixed
variant a non-additive coding can give complete separation and return `OR = inf`, an enormous `SE`, and
*P* underflowed to exactly `0`. Because `0 < 5e-8`, such a row satisfies any significance test and
becomes a phantom genome-wide peak. λ_GC is unaffected either way — that computation already requires
*P* > 0. `n_degenerate` is reported per model in `scan_qc.tsv` beside `n_errcode`.

## 6. Peaks — two tiers, one merge

Peaks are called **once**, at the suggestive threshold, then labelled. The thresholds are nested, so
every genome-wide peak is a suggestive peak whose strongest variant also clears `params.PGenomeWide`;
calling twice would count the same physical peak in both tiers.

- **Merging is by distance, not LD clumping.** Variants passing `params.PSuggestive` are joined when
  they lie within 2 × `params.PeakFlank` of one another; the window is the variant span extended by
  `PeakFlank` on each side. A distance rule needs no reference panel, is reproducible from the summary
  statistics alone, and cannot be perturbed by the choice of LD sample — which matters precisely
  because three LD sources are compared downstream and none of them may be allowed to define the peak.
- **Lead** = smallest *P* in the peak, under the peak model.
- **Tier** = `genome_wide` if the lead clears `params.PGenomeWide`, else `suggestive`.

### What each tier is for

| | genome-wide | suggestive |
|---|---|---|
| annotation table (§8) | ✔ | ✔ |
| cross-cohort table (§12) | ✔ | ✔ |
| landscape figure | ✔ | ✔ |
| conditional analysis | ✔ | ✔ — against the **suggestive** threshold |
| SuSiE fine-mapping | ✔ | ✔ |
| three-source regional plot | ✔ | ✔ |

Both tiers now receive the full per-peak follow-up, and their outputs are separated by tier
(`03.peaks/<tier>/`, `figures/0{2,3,4}.*/<tier>/`) so the two are never mixed. **The conditional
analysis is judged against the threshold that defined the peak**: a suggestive peak tested against
5 × 10⁻⁸ would fail at round 0 every time, and the table would report zero signals everywhere —
meaningless rather than empty.

What separates the tiers is now interpretation, not machinery. A suggestive peak is still a
description of the scan's shape rather than a finding, and the follow-up exists so that shape can be
inspected — not so that 128 loci can be reported as results.

**The suggestive tier describes the scan's shape; it is not a list of findings.** Over *M* analysed
variants a threshold α yields ≈ *M*α crossings by chance alone, and an inflated scan yields a multiple
of that. Fine-mapping signals of that strength at a few thousand effective samples returns empty
credible sets, so the tier receives no per-peak follow-up and the figures state the expected count so
a reader cannot mistake it for a hit list.

## 7. Calibration: λ is reported, never applied

λ_GC is a read-out, not a correction. Where relateds are already removed and global PCs already
fitted, residual λ > 1 is fine-scale structure that global PCs do not absorb — and it matters
precisely when case status is correlated with that structure, which happens whenever the sample sets
retain cases and controls at different rates.

Consequences, which hold for any dataset this is run on:

- The peak list is a **candidate** list, not a final answer.
- λ is never applied. Dividing χ² by λ rescales every statistic uniformly, including any true signal,
  and hides the structure rather than fixing it.
- The follow-up for residual structure is a **GRM random-effect model** (SAIGE / REGENIE) on the same
  samples. That is a separate component; this one does not attempt it.

**LDSC is deliberately not part of this component.** It is the standard way to split confounding from
polygenicity and looks like it fits this question, but it needs samples in the tens of thousands; at a
few thousand effective samples the intercept's standard error swamps the estimate. Recorded so it is
not added later by reflex.

## 8. Lead-variant annotation

One row per peak lead, both tiers, in `03.peaks/lead_annotation.tsv`.

- `EA` = plink2's `A1`; `OA` = the other of REF/ALT; `Beta` = log(OR). `OR`, `L95` and `U95` are
  **three numeric columns**, never one packed string — a formatted `1.591 (1.350–1.874)` cannot be
  used without parsing it back apart. Format for display at the point of display.
- Per-group genotype counts, EAF, missing rate and HWE come from **two plink2 calls per cohort** —
  `--keep` the case list / the control list, `--extract` the leads, `--geno-counts --hardy`. Keep lists
  are derived from the `.fam` phenotype column (2 = case, 1 = control). `Case_Genotype_Distribution` is
  `hom-REF/het/hom-ALT` as counts. EAF is computed for the **effect** allele, flipping the ALT-based
  count when plink2 tested REF.
- **rsID** from `params.rsid_vcf`: any tabix-indexed VCF on the same genome build whose ID column
  carries the rsID, looked up by `chr:pos:REF:ALT`. Variants absent from it keep `.`.

  Full-key matching finds the right *record*, which is necessary but not sufficient. A `norm`ed VCF
  splits multi-allelic records by **copying the whole ID field onto every split allele** instead of
  distributing the accessions, so one ID string can end up stamped on several alleles at one
  position. That is only a defect when the string holds **more than one accession**, because RS
  numbers are site-level:

  | at one position | example | reading |
  |---|---|---|
  | one accession, several alleles | `rs790041` on `C>T` and `C>G` | the RS names the site; it holds for whichever ALT we took. **Resolved** |
  | one accession, one indel written from two anchors | `rs75682226` on `GGGCTT>G` and `G>GGGCTT` | same variant, two representations. **Resolved** |
  | several accessions, one allele | a dbSNP merge | both name our allele. **Resolved**, both reported |
  | several accessions, several alleles | `rs1347066655;rs1564211184` on `CGTG>C` and `C>T` | the ID field is the UNION over the pre-split record's alleles. At most one names ours and the file no longer says which. **Unresolved** |

  The key rule therefore extends to the **accession**: several accessions are attributed to a variant
  only when no other allele at that position carries the same ID string. When another does, `rsID`
  stays `.` and the candidates are kept in `rsID_unresolved`, so nothing downstream — a table, a
  Manhattan label, a forest label — can present an accession that may belong to a different allele.
- **Gene** from the Ensembl gene models via `gene_annotate.R`: the overlapping gene if there is one
  (preferring protein-coding, then the longest model), else the nearest, with `Gene_Distance_bp`.
  The same "informative gene" rule as the regional plot (§9) — a variant is never labelled with a
  clone accession.

## 8b. Peaks per model, two tiers

The peak model is the only one that **defines the fan-out**, and nothing downstream reads another
model's peak. But a scan that produces peaks and never names them is not reported, only run — so peaks
are called for every model, in `03.peaks/model_peaks.tsv` (positions) and
`03.peaks/model_peaks_annotation.tsv` (the full field set of §8).

- Same rule as the primary peak caller: `usable_mask` (§5), then *P* < `PSuggestive`, merge within
  `PeakFlank`, then tier on `PGenomeWide`. One row is one locus. For the peak model this reproduces
  `peaks.tsv` exactly, and that duplication is deliberate — `peaks.tsv` must stay byte-identical to
  keep the expensive LD/SuSiE tasks cached.
- **Each scan figure draws its own model's rows.** A figure drawn with another model's peaks says
  nothing about the scan it is drawn from.
- **The two tables are annotated separately.** A variant can be a peak under more than one model, so
  annotating a union keyed by variant ID hands one model's effect to the other. Rows are matched **by
  position, never by variant ID**, and an assertion checks that every row's *P* agrees with the tier it
  claims.

**Read non-primary peaks against their own calibration, not beside the primary results.** A deflated
model is tested against a null that is too narrow, and its peak count is an upper bound: the variants
plink2 could fit are a biased subset of the call set. The figures carry this caveat in their captions;
it is not left to this document alone.

## 9. Gene models

`_shared/scripts/gene_utils.R` is the single definition, shared by the regional-plot gene track and the `Gene`
column, so the two can never disagree.

- **Exon/intron structure**, one row per gene, drawn from a **representative transcript**: the
  transcript with the greatest total exonic length, preferring `protein_coding`. Where the annotation
  release carries no canonical or MANE flag this heuristic stands in for one. Mixing in a newer release
  for the gene track alone would make it inconsistent with the `Gene` column, which is why one release
  serves both.
- **Only genes with an official symbol are drawn.** Ensembl carries many entries named only for the
  sequencing clone they were called from, and in a dense window these can outnumber the named genes.
  They are suppressed by pattern
  (`^(AC|AL|AP|AF|BX|CR|CU|FP|FO|LL|Z)\d{6}\.\d+$`, `^(RP|CTD|CTA|CTB|CTC|LA|XX|KB|DK|WI)\d*-`), along
  with bare `ENSG` names, `TEC` and pseudogene biotypes. HGNC-approved non-coding symbols (`LINC…`,
  `MIR…`, `SNOR…`) are **kept** — they are real names, not clone accessions.

## 10. Per-peak follow-up (every peak, both tiers)

Conditional analysis, SuSiE and the three LD sources run for **every peak in both tiers** — 131 peaks
here (40 / 38 / 53 across the three cohorts), not the 3 genome-wide ones. What `params.MaxSuggestive`
caps is the number of suggestive peaks that get **figures**, not the number that get computed: the
tables are complete, and `lead_variants.tsv` marks the drawn subset with its `figure` column.

The distinction matters because the cross-cohort locus rule reads the conditional analysis of every
peak (§12a) — capping the computation would have made locus identity depend on which peaks happened
to be plotted.

### 10a. Conditional analysis
Stepwise inside the peak window, with the same samples, model and covariates as the genome-wide scan.
Round 0 is the unconditioned fit restricted to the window; each later round adds the previous round's
top variant to `--condition-list` and re-fits. The loop stops when nothing in the window still clears
**the peak's own tier threshold** — `params.PGenomeWide` for a genome-wide peak, `params.PSuggestive`
for a suggestive one (`peakThreshold(peak_id)` in the workflow) — or after `params.MaxCondRounds`.
Judging a suggestive peak against 5 × 10⁻⁸ would have ended every one of them at round 0.

The number of rounds that yielded a signal is the number of independent signals. If plink2 produces no
association file for a round, that is recorded as an outcome — every variant collinear with the
conditioning set — not treated as a crash.

### 10b. SuSiE fine-mapping
`susie_rss(bhat, shat, R, n, L = params.SusieL)`, credible sets at `params.SusieCoverage`,
`min_abs_corr = params.SusieMinAbsCorr`, with an **in-sample LD matrix**
(`plink --r square --keep-allele-order`, signed *r*).

Four requirements, each a real defect in the obvious implementation:

1. **`susie_get_cs(..., Xcorr = R)`, not `X = R`.** In susieR's signature `X` is the genotype matrix
   and `Xcorr` the correlation matrix. Passing an LD matrix as `X` makes susieR read *p* variants as
   *p* "samples", so the reported purity is not credible-set purity and the `min_abs_corr` filter is
   silently inert.
2. **Align `bhat`/`shat` to the LD matrix by variant ID, with an assertion.** Relying on the extract
   list and the `.bim` happening to agree in order is not a guarantee.
3. **`n` = the GWAS analysis N**, not the LD-reference size.
4. **LD from all analysed samples**, not controls only: `R` must describe the sample the summary
   statistics came from.

Signed `--r` for SuSiE, squared `--r2` for plot colouring, kept separate — `r²` discards the sign and
would make every pair of variants look positively correlated.

### 10c. LD from three sources, one scale

Three **roles**, whose configured resources and display names are `params.PopR2Template`,
`params.RefPanelBfile` and `params.LdPanelLabels`:

| role | key | what it is | reporting floor |
|---|---|---|---|
| in-sample | `cohort` | the samples the statistics came from | none |
| population panel | `tommo` | published co-occurrence *r²*, tabix-indexed | typically **r² ≥ 0.2 only** |
| reference panel | `1000g_eas` | an out-of-sample genotype panel | none |

**A reporting floor is a bound, not missing data.** Two states must be told apart and only one is
unknown:

- variant **is** in the panel but returns no pair with the lead → **r² below the floor**, plotted in
  the lowest bin like any other low-LD variant;
- variant is **not in the panel at all** → genuinely unknown, plotted grey.

Panel membership is read off the self-pairs the co-occurrence table already carries, so it costs no
extra query. All three sources share one binned scale (< 0.2 / 0.2–0.4 / 0.4–0.6 / 0.6–0.8 / 0.8–1.0),
which is what makes the truncation harmless: bounded values fall entirely inside the lowest bin.
Treating them as missing instead renders most of the panel grey and makes the source look useless.

**The reference-panel window is cut out with plink2 before plink 1.9 sees it.** Running `--r2` directly
against a panel of tens of millions of variants makes plink 1.9 size its workspace from the whole file
and request tens of GB, which fails.

## 11. Multiplicity

`params.PGenomeWide` is applied identically to every cohort. No further adjustment is made across
cohorts, and this is a decision. Where the sample sets are nested the scans are nowhere near
independent tests; a Bonferroni factor equal to the cohort count would be badly conservative, and a
correct correction needs the effective number of independent tests, which is not identified here.

## 12. Cross-cohort reporting

`_comparison/tables/lead_crosscohort.tsv` reports **every peak lead of every cohort, in every cohort**,
whether or not that cohort called a peak there. One row per (variant × cohort), carrying the same field
set as §8 plus `called_peak`.

The reason is structural. Where the sample sets are nested, a variant that is a peak in one of them
**always has an estimate in the others**; there is no such thing as a missing value, only an unreported
one. A blank cell is indistinguishable from a variant that could not be tested, and dropping such rows
silently removes a cohort from a forest plot. `called_peak` states the distinction explicitly and is
never inferred from a blank. A variant genuinely absent from a cohort's call set is `not_in_call_set`,
which is a different fact again.

**This is description, not replication.** Nested cohorts share most of their samples by construction,
so agreement between them is expected and carries little independent information. What the rows do
support is reading how one estimate behaves as the sample grows: an interval that narrows while the
point estimate holds is a variant gaining sample size, and an estimate that drifts toward 1 as a filter
is relaxed is what a structure-driven signal does. Neither is evidence of replication.

### 12a. Locus identity — what counts as "the same signal"

The union of peak leads is not a list of loci. Three cohorts scan overlapping samples, so the same
signal is often led by a *different* variant in each of them, and a naive union counts one signal two
or three times. `lead_crosscohort.tsv` therefore carries `locus_id` and `locus_rep`, assigned in
`cross_cohort.py` in two steps:

1. **Group by proximity.** Leads on one chromosome within `params.PeakFlank` (250 kb) of the previous
   one become a candidate group.
2. **Split only on evidence.** A candidate group is split back into individual loci **only when some
   cohort's conditional analysis reported two of its members as separate independent signals of one
   peak** (`<peak>.signals.tsv`, §10a). Otherwise the group is one locus.

The representative is the group's smallest-*P* member across cohorts, and it is the variant the
figures draw — `cohort_compare` and `model_compare` plot one row-group per locus, at its
representative, so all three rows compare the *same* variant.

Measured here: **84 lead variants → 75 loci** (9 merged; conditional analysis split 0 candidate
groups) for the fixed-effects component, and **52 → 47** (5 merged) for SAIGE.

Two consequences worth stating:

- **Proximity alone would have been wrong, and so would splitting alone.** *HLA-DQA1* had two leads
  137 bp apart, each the lead in a different cohort; conditioning on either left the window's best
  residual at 5.6 × 10⁻³, so they are one signal and are now drawn once. Conversely a gene *may*
  legitimately appear twice — it just requires conditional evidence, which no group in this data has.
- **`panel` is decided once per locus, on its representative,** then stamped on every member. Without
  that a locus could reach panel (b) through one of its leads and panel (c) through another.

`_comparison/figures/cohort_manhattan.png` reports the same sample sets at scan level: one row per
cohort, Manhattan and QQ on the same 2.5 : 1 grid the scan figure uses, so each Manhattan occupies the
identical 4.08 × 1.77 in box as `scan.<model>.png` panel (a) and a peak has the same shape on both.
The Manhattans sit on **one** cumulative chromosome-offset map, built from the *union* of
per-chromosome extents rather than per cohort, so a genomic column is the same x in every row. They
also share one −log₁₀ *P* limit, one data height and one tick locator, without which a taller column
could be a rescaled axis rather than a smaller *P*. The genome-wide loci are named once, over the union
across cohorts, and a dotted guide drops through every row; a diamond appears only in the cohorts that
called the locus themselves, which is what separates "genome-wide here" from "present here".

## 13. Figures

All figures go through `_shared/scripts/plot_style.py`, a collision check, and a visual pass. Every PNG has a
companion `.md` written by `_shared/scripts/figure_doc.py` giving the full panel-by-panel explanation, the
concrete numbers behind that rendering, how to read it, and what it cannot answer.

- **A per-locus figure names its variant above the crop line.** `regional`, `finemap` and
  `conditional` open with `Gene · CHROM:POS:REF:ALT · rsID` for the lead. The caption names the
  cohort and the locus, but the caption sits *below* `plot_h` — crop the publication figure out and
  it goes too, leaving a locus plot that no longer says which locus it is. The rsID appears only
  when it is resolved to this allele (§8); otherwise the line reads `rsID unresolved` and the
  variant id above it is the identifier that holds. Details in
  [../../_shared/docs/FIGURES.md](../../_shared/docs/FIGURES.md).
- **Nothing about the run is written into the figure code.** The covariate set, the LD panel names and
  the cohort names all arrive as parameters, so a figure describes the run that produced it.
- **No downsampling.** Manhattan and QQ draw every analysed variant; the scatter is rasterised so the
  file stays small while axes and text remain vector. Thinning would distort the QQ's lower arm, which
  is the region λ is read from.
- **λ is shown as a deviation from 1**, on stems anchored at 1.0 — a bar chart with a truncated
  baseline misrepresents the size of the deviation.
- **One figure per scan.** Manhattan, that scan's own peaks on the same genomic axis, QQ and
  effect-against-frequency are one figure, not two. The peak panel is stacked directly under the
  Manhattan and shares its x-axis, so a peak can be traced up into the column it came from; that
  adjacency is the reason to merge them rather than merely co-locate them.
- **Gene labels are capped, and the cap is stated.** Every genome-wide lead is labelled; only the
  `params.MaxSuggestive` smallest-*P* suggestive leads are, because a full suggestive tier needs more
  text than a 6.5 in axis holds. Labels follow the **gwaslab** idiom at `params.AnnoStyle = auto`
  (`expand` and the other gwaslab styles remain selectable): one band, italic text at the least
  intrusive rotation that fits, an L-shaped leader arm whose vertical segment marks the peak, and
  positions repelled symmetrically at `params.RepelForce`.
- **The label band is outside the data area, and there is one band per panel.** The axes box is shrunk
  to make the strip, so the data limits are never inflated for text. Genome-wide and suggestive names
  share one band and one rotation, distinguished by colour and weight.
- **The regional figure states its LD coverage in its document, not in the panels.** The per-source
  counts (`measured` / `below_threshold` / `not_in_panel`) are provenance rather than a finding, and
  in-panel they sat over the data in all three LD panels. They go to the family catalogue as
  columns, one row per locus, which also makes them comparable across loci.
- **The caption carries three things: the finding, one line per panel, one data line.** The
  interpretation, the symbol glossary and the estimator live in the figure's document. These are
  laid out as slides, and a caption that fills most of the canvas makes the plot the smaller half.
- **One regional figure per peak**, three LD panels on a shared x-axis, because LD concordance across
  sources is the point of the figure.
- **Collisions are prevented by measurement, not by tuned offsets.** `plot_style` owns the helpers and
  every figure routes through them: `fit_left_margin` (the left margin is measured from the widest y
  tick label, so a long label cannot be clipped off the canvas), `value_labels` (a point's value label
  grows the axis to make room for itself), `spread_labels` (overlapping text is separated on rendered
  boxes), `equalise_row_heights` / `align_panel_tops` (panels meant to be compared are given the same
  data height) and `place_legend` (above the axes when the panel is short, otherwise in the corner
  holding the fewest data points). A fixed offset works until the data changes shape; a measurement
  does not.
- **Tick labels that would collide are dropped**, measured on the rendered figure
  (`plot_style.thin_tick_labels`), rather than drawn on top of each other.
- **No power or minimum-detectable-effect figure.** Removed by decision. Without replication its only
  load-bearing statement is that a peak sits at the detection boundary, which is arithmetically the
  same statement as "the peak barely cleared the threshold" — already visible from *P*. Per variant,
  `|log OR| / log(MDE)` reduces to `(|z| / 5.451) × (SE_obs / SE_pred)`, so it adds no information
  to *P* that was not already there.

## 14. External resources

The component needs the following, all configured in the pipeline's `SITE CONFIGURATION` block. What
each must *provide* is fixed; which release fills it is not.

| purpose | param | requirement |
|---|---|---|
| PLINK 1.9 | `plink19` | `--r` / `--r2 square` |
| PLINK 2 | `plink2` | `--glm`, `--geno-counts`, `--hardy`, window extraction |
| tabix | `tabix` | indexed access to the resources below |
| rsIDs | `rsid_vcf` | tabix-indexed VCF on the analysis build, rsID in the ID column |
| population LD | `PopR2Template` | per-chromosome co-occurrence *r²* tables, `@@CHROM@@` in the path |
| reference panel | `RefPanelBfile` | plink bfile on the analysis build, out-of-sample |
| recombination | `recomb_bw` | bigWig of recombination rate on the analysis build |
| gene models | conda `conda_env_r` | an Ensembl `EnsDb` package plus `rtracklayer` |
| consequence | `snpeff_index_dir` | pre-computed snpEff index, one tabix TSV per chromosome, keyed `chr:pos:REF:ALT` |

The versions and paths in use for this project are recorded in [`STUDY_NOTES.md`](STUDY_NOTES.md).
