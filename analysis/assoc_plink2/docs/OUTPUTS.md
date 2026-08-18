# Outputs — `assoc_plink2`

Everything lands under `results/`. Per-cohort directories are numbered in the order the pipeline
produces them, and figures are grouped by kind so a directory listing reads as the analysis reads.

```
results/
  <cohort>/
    01.assoc/<model>/        genome-wide statistics (symlink) + plink2 log (copy)
    02.scan/                 scan_qc.tsv — cohort size, calibration, exclusions, all 3 models
    03.peaks/                peaks.tsv · lead_variants.tsv · lead_annotation.tsv
                             model_peaks.tsv · model_peaks_annotation.tsv  (all 3 models)
      <peak_id>/             genome-wide peaks only: LD, SuSiE, conditional
    figures/
      01.scan/               scan.<model>.png + .md                  (3)
      02.regional/           regional.<peak_id>.png + .md            (genome-wide peaks)
      03.finemap/            finemap.<peak_id>.png + .md             (genome-wide peaks)
      04.conditional/        conditional.<peak_id>.png + .md         (genome-wide peaks)
  _comparison/
    tables/                  scan_qc_all.tsv · peaks_all.tsv · lead_annotation_all.tsv
                             model_peaks_all.tsv · cs_variants_all.tsv
                             lead_crosscohort.tsv   <- every lead, in every cohort
    figures/                 cohort_compare.png + .md
                             cohort_manhattan.png + .md
  _run_info/                 trace · report · timeline · dag · run_manifest
```

**Every PNG has a companion `.md` of the same name**, written by `scripts/figure_doc.py`: what the
figure shows, each panel's statistic and how it is computed, axis and colour semantics, the concrete
numbers behind that particular rendering, how to read it in order, and what it cannot establish.

---

## `01.assoc/<model>/`

| file | mode | content |
|---|---|---|
| `<cohort>.<model>.<covar tag>.PHENO1.glm.logistic` | symlink | one row per variant; the covariate tag is `covarTag()` |
| `*.log` | copy | plink2 log; check `covariates loaded` = 11 |

Published as a **symlink** into `work/` — nine copies of 585 MB is 5.3 GB of duplication for files
that are fully reproducible. Small files are copied, so tables and figures survive `work/` being
cleaned.

Columns are plink2's own: `#CHROM POS ID REF ALT PROVISIONAL_REF? A1 OMITTED A1_FREQ TEST OBS_CT OR
LOG(OR)_SE L95 U95 Z_STAT P ERRCODE`. `hide-covar` means exactly one `TEST` value per file.

## `02.scan/scan_qc.tsv`

One row per model. This is the **only** output dominant and recessive contribute to besides their own
Manhattan/QQ figure.

| column | meaning |
|---|---|
| `n_case` / `n_ctrl` | phenotype counts from the `.fam` (2 = case, 1 = control) |
| `n_eff` | 4 / (1/`n_case` + 1/`n_ctrl`) — the balanced-design equivalent size, reported as a cohort descriptor |
| `n_variants` | rows in the GLM file |
| `n_errcode` | rows plink2 flagged with a non-`.` ERRCODE |
| `errcode_breakdown` | `CODE=count;…` |
| `n_degenerate` | rows plink2 left `ERRCODE = '.'` that still carry no usable result — non-finite OR/SE, or P outside (0, 1]. See METHODS §5; these are what produced three phantom loci before the guard existed |
| `n_analysed` | rows surviving both checks |
| `obs_ct_median` | median non-missing genotype count |
| `lambda_gc` | median χ² / 0.4549, computed on analysed rows only |
| `n_genomewide` / `n_suggestive` | analysed variants passing each threshold |
| `min_p` | smallest analysed P |

## `03.peaks/`

| file | content |
|---|---|
| `peaks.tsv` | one row per **additive** peak, both tiers |
| `lead_variants.tsv` | the fan-out key: `cohort, peak_id, tier, lead_id, chrom, pos, start, end` |
| `lead_annotation.tsv` | the table below — **read this one first** |
| `model_peaks.tsv` | **two-tier peaks of all three models**, merged at 250 kb |
| `model_peaks_annotation.tsv` | the same peaks with the full `lead_annotation` field set |
| `<peak_id>.sumstat.tsv` | every analysed variant in the window; genome-wide peaks only |

`peak_id` is `<gw|sg><NNN>_<chrom>_<lead_pos>` — the prefix encodes the tier, so a directory listing
is already sorted by importance. It is the join key for everything downstream.

`peaks.tsv` columns: window (`chrom/start/end`), `n_sig_variants`, `n_genomewide_variants`, and the
lead's `ID/pos/P/OR/SE/L95/U95/A1/REF/ALT/A1_FREQ/OBS_CT`.

### `lead_annotation.tsv`

One row per peak lead, **both tiers**.

| column | how it is derived |
|---|---|
| `rsID` | `params.rsid_vcf`, matched on `chr:pos:REF:ALT`; `.` if absent from that VCF |
| `Gene`, `Gene_Biotype`, `Gene_Distance_bp` | EnsDb v86 — overlapping gene, else nearest. `Gene_Distance_bp` is 0 inside a gene, negative when the variant lies before the gene, positive after |
| `EA` / `OA` | plink2's `A1` / the other of REF, ALT |
| `Beta`, `SE` | `log(OR)` and `LOG(OR)_SE` |
| `OR`, `L95`, `U95` | odds ratio and its 95 % CI, as **three numeric columns**. Never one packed string: a packed `OR (L95–U95)` cannot be used without parsing it back apart, and it makes an en-dash load-bearing. Format for display at the point of display |
| `P` | additive P at the lead |
| `Case_Genotype_Distribution` | `hom-REF/het/hom-ALT` **counts** among cases |
| `Case_EAF` | frequency of the **effect** allele among called cases |
| `Case_Missing_Rate` | missing / (called + missing), cases only |
| `Case_HWE_P` | plink2 `--hardy` on the case subset |
| `Control_*` | the same four, on controls |
| `A1_FREQ`, `MAF`, `OBS_CT`, `N_case`, `N_ctrl` | cohort-level context |

Genotype counts, EAF, missing rate and HWE come from two plink2 calls per cohort (`--keep` the case
list / the control list, `--extract` the leads, `--geno-counts --hardy`).

### `03.peaks/<peak_id>/` — genome-wide peaks only

| file | content |
|---|---|
| `<peak_id>.ld_cohort.tsv` · `.ld_tommo.tsv` · `.ld_1000g_eas.tsv` | `ID, r2, state` per source |
| `<peak_id>.ld_coverage.tsv` | per source: `n_measured`, `n_below_threshold`, `n_not_in_panel`, `pct_informative` |
| `<peak_id>.ld_matrix.tsv` / `.vars` | signed square *r* matrix and its row order — the input SuSiE uses |
| `<peak_id>.exons.tsv` | exon structure of the representative transcript of each named gene |
| `<peak_id>.recomb.tsv` | recombination rate in the window, from `params.recomb_bw` |
| `<peak_id>.pip.tsv` / `.cs.tsv` / `.susie.json` | per-variant PIP, credible sets, run metadata |
| `<peak_id>.cs_variants.tsv` | **one row per credible-set member** — the file downstream annotation reads |
| `<peak_id>.rounds.tsv` / `.signals.tsv` / `.round<N>.tsv` | conditional analysis |

**`state` takes three values and they are not interchangeable:**

| state | meaning |
|---|---|
| `measured` | *r*² with the lead was reported by that source |
| `below_threshold` | the variant **is** in the panel but the source's reporting floor hid the pair — for a co-occurrence panel this is typically *r*² < 0.2, a bound, not missing data |
| `not_in_panel` | absent from that panel entirely — genuinely unknown |

Only `not_in_panel` is grey in the figures. See METHODS §10c.

### `cs_variants.tsv` — the credible-set variant table

`cs.tsv` holds the members only as a comma-joined string and `pip.tsv` only as a `cs` column; neither
is a usable input. This table gives one row per member, with its set context, its posterior weight and
rank, and three annotations.

| column group | columns |
|---|---|
| set context | `cohort, peak_id, cs, cs_size, cs_coverage, cs_purity_min_abs_corr` |
| variant | `variant_id, chrom, pos, EA, OA, A1_FREQ, P, OR, L95, U95, Beta, SE` |
| posterior | `pip, pip_rank, cum_pip, is_lead, is_top_pip` |
| identity | `rsID` (`params.rsid_vcf`), `Gene`, `Gene_Biotype` (Ensembl gene models, nearest/overlapping) |
| consequence | `snpEff_Effect, snpEff_Impact, snpEff_Gene, snpEff_GeneID, snpEff_Feature, snpEff_FeatureID, snpEff_Biotype, snpEff_HGVS_c, snpEff_HGVS_p, snpEff_n_annotations, snpEff_all_effects` |

`pip_rank` and `cum_pip` are **within a set, PIP-descending**, so reading down the rows is reading
down the posterior.

**snpEff comes from a pre-computed index** (`params.snpeff_index_dir`) — one tabix-indexed TSV per
chromosome, keyed `CHROM/POS/REF/ALT`. Nothing here runs snpEff.

A variant carries **one row per transcript** in that index, so the lookup collapses them to the most
severe by `IMPACT` (HIGH > MODERATE > LOW > MODIFIER), preferring a `protein_coding` feature on ties.
The discarded rows are not lost: `snpEff_n_annotations` counts them and `snpEff_all_effects` lists the
distinct effects, so a MODIFIER call with a MODERATE sibling transcript is visible as such.

`snpEff_Gene` and `Gene` answer different questions and are both kept: the first is the transcript the
variant falls in, the second is the nearest or overlapping gene by position.

## `model_peaks_annotation.tsv` — peaks of every model, both tiers

One row per (model × peak), with the **same field set** as `lead_annotation.tsv` so the two stack.
Prefix columns: `cohort, model, peak_id, tier, n_sig_variants, n_genomewide_variants`.

Peaks are formed exactly as the additive ones are (METHODS §8b): usable rows, *P* < 10⁻⁵, merged within
250 kb, then tiered by whether the lead clears 5 × 10⁻⁸. Roughly 29–53 peaks per model per cohort.

`peak_id` is `<add|dom|rec><NNN>_<chrom>_<pos>`. **Nothing fans out from this file** — it exists so each
scan figure can draw its own peaks and so a DOM/REC result is reportable at all. Read the recessive
rows against λ_REC ≈ 0.65 (METHODS §8b).

**The additive rows duplicate `peaks.tsv`, deliberately.** `peaks.tsv` is the fan-out key for the
expensive LD/SuSiE/conditional tasks and must stay byte-identical for those to stay cached, so the
general table is written alongside it rather than replacing it. The two agree row for row; if they ever
disagree, `peaks.tsv` is the one the pipeline acted on.

## `_comparison/`

| file | content |
|---|---|
| `tables/scan_qc_all.tsv` | every cohort × model row, stacked |
| `tables/peaks_all.tsv` | every cohort's peaks, stacked |
| `tables/lead_annotation_all.tsv` | every cohort's annotated leads, stacked (also carries the snpEff columns) |
| `tables/model_peaks_all.tsv` | every cohort's peaks, all three models, both tiers, stacked |
| `tables/cs_variants_all.tsv` | every credible-set member of every genome-wide peak, stacked |
| `tables/lead_crosscohort.tsv` | **every peak lead, in every cohort** — see below |
| `figures/cohort_compare.png` + `.md` | additive calibration + N_eff; case-vs-control EAF at every lead; forest of every genome-wide lead in all three cohorts |
| `figures/cohort_manhattan.png` + `.md` | the three additive scans stacked — Manhattan + QQ per row, on the scan figure's own 4.08 × 1.77 in box; one shared genomic axis, y-limit, data height and QQ limit; header carries cases, controls and N_eff, the QQ carries λ_GC |

### `lead_crosscohort.tsv` — every lead, in every cohort

One row per (**peak lead variant** × **cohort**), for every variant that is a peak in *any* cohort,
both tiers — (number of distinct leads) × (number of cohorts) rows.

| column group | columns |
|---|---|
| variant | `variant_id, chrom, pos, best_tier` |
| cohort context | `cohort, called_peak, peak_id` |
| everything else | the full `lead_annotation.tsv` field set (§8), including the per-group genotype block and the snpEff columns |

`best_tier` is the strongest tier the variant reached in **any** cohort, so sorting on it puts the
genome-wide leads first. `called_peak` is what **this** cohort made of it:

| value | meaning |
|---|---|
| `genome_wide` | this cohort called a genome-wide peak here |
| `suggestive` | this cohort called a peak, but only at the weaker tier |
| `not_a_peak` | this cohort called no peak here — **the estimate is still reported** |
| `not_in_call_set` | the variant is absent from this cohort's call set; the statistics columns are empty |

The cohorts are nested, so a lead always has an estimate in the larger sets. Reporting only the cohort
that called the peak leaves a blank that cannot be told apart from a missing estimate — which is what
the earlier comparison figure did. See METHODS §12 for why this is description and not replication.

## `_run_info/`

| file | content |
|---|---|
| `trace.txt` | one row per task: status, exit, attempt, duration, peak RSS, I/O |
| `report.html` · `timeline.html` · `dag.html` | Nextflow execution report, timeline, DAG |
| `run_manifest.json` | every parameter and resource path the run used |

`trace.txt` is what distinguishes a cached task from a completed one after the fact — the ambiguity
that once made a verification script report 150 cached tasks as failures.
