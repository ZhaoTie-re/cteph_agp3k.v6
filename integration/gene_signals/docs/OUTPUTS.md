# gene_signals — outputs

```
results/
  gene_signals.xlsx                 the workbook (9 sheets)
  workbook.json                     sheet row counts, note widths, Open Targets records listed
  gene_signals.genes.long.tsv       gene x variant set x cohort
  gene_signals.variants.long.tsv    gene x variant x cohort (no sample IDs)
  gene_signals.carriers.tsv         gene x set x cohort x group x subgroup (case platform / control source)
  gene_signals.scans.tsv            cohort x set x statistic
  gene_signals.multisite.tsv        multi-site carriers by sample ID, with sequencing and clinical information
  gene_signals.lovo.tsv             leave-one-variant-out: gene x set x cohort x statistic x left-out variant
  checks.json                       counts and cross-checks (all_passed)
  00.selection/    genes.tsv variants.tsv selection.json
  01.genotypes/    geno.<cohort>.tsv carriers.<cohort>.tsv samples.<cohort>.tsv multisite.<cohort>.tsv
                   (geno / samples / multisite: sample IDs, analyst files)
  02.annotation/   function.tsv
  03.frequencies/  tommo.tsv kg.tsv annovar.tsv
  04.opentargets/  opentargets.tsv opentargets_detail.tsv opentargets_meta.json raw/<variant>.json
  05.rv_qc/        rv_qc.<cohort>.json map_recount.<cohort>.json variant_sets.terms.tsv variant_sets.flow.tsv variant_sets.json
  06.lovo/         lovo.<cohort>.<set>.tsv (the rvtest results of the re-tests) lovo.<cohort>.<set>.<cmc|skato>.log
  _run_info/       run_manifest.json trace.txt report.html timeline.html dag.html
  _superseded_2026-10-09/   outputs of the first version no longer produced (kept, not used)
```

## The workbook

Every sheet except Rare-variant QC has the frame of the common-variant table:

- **rows 1–2:** title and descriptor;
- **rows 3–5:** header (block spanner, sub-spanner, leaf, with footnote letters);
- **then:** the table, a closing rule and the notes in column A, in a fixed order:
  - lettered notes keyed to the header letters, each led by its column name in bold;
  - where they apply, Cohorts (codes, sizes, nesting), Model, Marks and Shading;
  - **Missing values**, one line per placeholder and column (—, –, n/a, empty) with its meaning;
  - Abbreviations.

Cohort codes are N, I, F in the cohort colours (narrow `0571B0`, intermediate `008837`, full `B2182B`).
Variant sets are M/H (moderate/high-impact) and L/M/H (low/moderate/high-impact).

Sheets in the order of the workbook:

| sheet | layout and columns |
|---|---|
| **Rare-variant QC** | The project's rare-variant QC table. Steps Input → Sample QC → Variant QC → Allele count, with Samples (case/ctrl) and Variants per cohort. Then the basis of κ: gap at κ = 1, gap at κ = 2, κ with no significant gap, κ chosen, and the check under the association covariates. Footnotes a–e and references |
| **Variant sets** | One row per step, with its Variants (one column per cohort) and its Genes (one column per cohort) side by side; an empty cell means the step concerns only variants or only genes. Rows: analysed variants and their highest snpEff impact; not mapped to a gene (MODIFIER only; HIGH / MODERATE / LOW annotation without a gene) and mapped to genes (in more than one gene). Per variant set: qualifying, → not tested (gene with 1 qualifying variant), → tested (gene with 2 or more), tested genes with 2 and with ≥ 3 variants, tested variants by impact, variants in more than one tested gene |
| **Consequence terms** | Per impact, a shaded total row ("All HIGH terms" …), then each snpEff term; a closing "All" row. Per cohort: Analysed · Tested genes M/H · Tested genes L/M/H ("–" where the impact is outside the set). Totals asserted equal to sheet Variant sets |
| **Gene-based scans** | One row per cohort × variant set. Two blocks side by side: *Primary analysis: genes with ≥ 2 variants* and *Sensitivity analysis: genes with ≥ 3 variants*. Each block has: Genes tested · Bonferroni threshold · and for SKAT-O and CMC, Smallest P · Significant genes (Bonferroni, BH) · Largest significant P (BH; "–" when none) · λ<sub>GC</sub> of that family. Footnote e gives the λ<sub>GC</sub> ranges of both families |
| **Gene-based signals** | One row per gene × variant set × cohort: per gene (Gene, Chr, Class merged over its six rows), the M/H set and then the L/M/H set (merged), each with its three cohorts. Variant burden: Variants · Carriers n (%) case / control · Cumulative MAC case / control · Cumulative MAF case / control · Multi-site carriers case / control. Gene-based tests (primary, ≥ 2 variants): SKAT-O ρ / P / q · CMC burden OR (95% CI) / P / q · Tier (merged per set). Sensitivity (≥ 3 variants): SKAT-O q · CMC q · Tier. A BH call is shaded and bold; † marks Bonferroni; "—" when the cohort tested the gene on < 3 variants (sensitivity) or the OR is not estimable |
| **Gene-based variants** | One row per gene × variant; genes in the order of Gene-based signals, variants in genomic order, gene label and Class on the gene's first row. Variant: Tested gene · Class · Set (M/H or L/M/H) · Variant (GRCh38 CHROM:POS:REF:ALT) · rsID · Chr · Position · EA/OA (EA = ALT) · Gene (Ensembl 86) · Used in (cohort codes). Functional annotation: Consequence · Impact · Transcript · HGVS · Other consequences. Open Targets: ClinVar · GWAS credible sets · L2G · molQTL credible sets · Enhancer-to-gene. Effect-allele frequency, with data bars: This study CTEPH / AGP3K (`AfCohort`) · ToMMo · 1KG JPT · 1KG EAS · gnomAD EAS · gnomAD NFE. Per cohort: Used in test (`Yes`, `No: MAC k`, `No: no carrier`, `No: common (MAF x)`, `No: depth-gap filter`) · Genotypes CTEPH / AGP3K (OA/OA / OA/EA / EA/EA / missing) · EAF · Missing · HWE P |
| **Leave-one-variant-out** | All 15 genes in both variant sets, side by side. A grey gene row (gene · class · variants each set tested, ‡ where a cohort tested 2) with the published tests; below, one row per variant either set tests in any cohort, indented, in genomic order as in Gene-based variants (asserted). Columns: Gene / HGVS · Variant (GRCh38 CHROM:POS:REF:ALT) · rsID · Consequence · Impact (both as in Gene-based variants, asserted) · Carriers case / control (N, I, F; grey when the cohort's tests did not use the variant, 0 / 0 without a carrier) · per set (M/H, then L/M/H): Tier (grey row: 1–3 green, S3 = sensitivity tier only, – not tiered; a tiered set has its −log10 P in bold) · SKAT-O and CMC −log10 P (N, I, F): grey row with all variants, variant rows without that variant. Shading (key under the table), red for a fall and blue for a rise in −log10 P from the grey row, in steps 1–2, 2–3, ≥ 3 (darker = larger); none below 1. Empty: a low-impact variant in the M/H columns. "—" where the cohort's test did not use the variant |
| **Gene-based Open Targets** | As the common-variant table's Open Targets sheet. A header row per gene (outline level 0), one row per variant (level 1), and its records below (level 2, collapsed). Columns: Gene · Set · Variant · rsID · Source · Type · Trait, condition or target gene · cis / trans · QTL tissue / cell type · Study · Publication · Study N · Cohorts · Population · P · Beta · P / Beta from · Set lead variant · PP · Fine-mapping · L2G gene · L2G score · E2G score · E2G best tissue / cell type · E2G tissues / cell types (n) · E2G interval · Clinical significance · Review status |
| **Multi-site carriers** | A grey header row per gene (samples: cases, controls), then one row per gene × sample, cases first. Columns: Sample ID · Group · Variants carried (HGVS, coloured by impact, het / hom) · Multi-site carrier in M/H and L/M/H (cohort codes) · Sequencing: Project · Platform (verbatim) · Target depth · Observed depth · Clinical: Sex · Age at DNA collection (cases from the project's sample table, controls from the AGP3K age / sex release) · Age at diagnosis · DVT · Max mPAP (mmHg) (cases only; n/a for controls). The only sheet with sample IDs and project names |

## `gene_signals.genes.long.tsv` (gene × set × cohort)

| column | meaning |
|---|---|
| `gene_order`, `set_name`, `gene_symbol`, `chrom`, `class`, `stratum`, `cohort` | identity |
| `tier`, `cohorts_called` | consistency tier of the set (0 = untiered); cohorts with a BH call |
| `n_case`, `n_ctrl`, `bonferroni_threshold`, `n_genes_tested` | cohort sizes; alpha / genes tested |
| `tested`, `not_tested_reason`, `n_var` | whether and on how many variants the cohort tested the gene |
| `carriers_*`, `mac_*`, `multi_*`, `maf_*` (case / ctrl) | the burden the test saw |
| `skato_p`, `skato_q`, `skato_bh`, `skato_bonferroni`, `skato_rho` | SKAT-O |
| `in_sens`, `*_q_sens`, `*_bh_sens`, `*_bonferroni_sens`, `bonferroni_threshold_sens`, `n_genes_tested_sens`, `tier_sens` | the sensitivity family (≥ `SensMinNumVar` variants) |
| `cmc_p`, `cmc_q`, `cmc_bh`, `cmc_bonferroni`, `cmc_or`, `cmc_l95`, `cmc_u95` | CMC: score-test P, Wald OR (95% CI) |

## `gene_signals.variants.long.tsv` (gene × variant × cohort)

| column | meaning |
|---|---|
| `set_name`, `gene_symbol`, `gene_order`, `variant_id`, `chrom`, `pos`, `OA` (REF), `EA` (ALT), `impact`, `effect`, `cohorts_<set>` | identity and gene-map annotation |
| `tx_gene`, `transcript`, `transcript_consequence`, `transcript_biotype`, `hgvs_c`, `hgvs_p`, `other_consequences`, `lof_nmd`, `n_annotations`, `other_genes` | snpEff |
| `Gene_EnsDb`, `Gene_Biotype`, `Gene_Distance_bp` | Ensembl 86 gene |
| `rsID`, `rsID_source` (`ToMMo` / `dbSNP allele` / `dbSNP 151`), `rsID_unresolved` | rsID |
| `tommo_*`, `kg_*`, `gnomad_af`, `gnomad_eas_af`, `gnomad_nfe_af`, `dbsnp_rsid` | reference lookups |
| `ot_*` | Open Targets summary (as the common-variant table) |
| `af_cteph`, `af_agp3k` | ALT frequency in the `AfCohort`'s cases / controls |
| `in_callset`, `callset_reason`, `ctrl_maf_split`, `minor_allele`, `mac` | presence in the cohort call set and why not; the control MAF where common |
| `case_*`, `ctrl_*` | homRef, het, homAlt, missing_n, called, alt_af, missing_rate, hwe_p, mac, carriers |
| `in_set_<set>`, `set_reason_<set>` | used by the cohort's test in that set; reason when not |

## `gene_signals.lovo.tsv` (gene × set × cohort × statistic × left-out variant)

| column | meaning |
|---|---|
| `gene_order`, `set_name`, `gene_symbol`, `stratum`, `cohort`, `method` | identity |
| `primary` | the set is the gene's primary variant set (docs/METHODS.md §5b) |
| `k`, `dropped` | 0 and empty: the test on all variants (baseline); else the k-th variant in genomic order, left out |
| `n_var`, `n_left` | variants the cohort tested the gene on; variants in this test |
| `n_var_engine`, `n_poly_engine` | rvtest's NumVar (asserted = `n_left`) and NumPolyVar |
| `p`, `p_all` | P of this test; P with all variants (= the published P) |
| `delta_log10`, `rank`, `top` | log10(`p` / `p_all`); its rank within gene × set × cohort × statistic; the largest positive |
| `rho`, `beta`, `se` | SKAT-O ρ; cmcWald β and SE (CMC) |

## `05.rv_qc/`

- **`rv_qc.<cohort>.json`.** n_in and n_kept (all / case / control); sample_qc (z threshold, removed);
  variant_qc (knee, n_in, n_removed, n_protected_rescued, n_mono_after_sample_qc, n_mono_removed_by_delta);
  minac (recommendation); sweep (κ 0, 1, chosen); model_check per κ, as `decision` / `association`, each
  with beta, CI, P and % of mean.
- **`map_recount.<cohort>.json`.** The gene map before the gene filter: analysed, mapped, mapped rows,
  multi-gene variants. Per set: qualifying variants (by impact), tested variants (by impact), tested
  multi-gene variants, tested variant–gene pairs, genes with 1 / 2 / ≥ 3 variants, genes tested, and
  variants only in untested genes.
- **`variant_sets.flow.tsv`.** cohort · section · stratum · item · value: the rows of sheet Variant sets.
- **`variant_sets.terms.tsv`.** cohort · impact · effect · analysed · `<set>`: the rows of sheet
  Consequence terms.
