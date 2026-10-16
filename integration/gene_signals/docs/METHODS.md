# gene_signals — methods

Each section is one step: what it computes and what it asserts before writing. A failed assertion stops the
run (`errorStrategy 'finish'`). The workbook follows the common-variant table (`saige_signals`) wherever the
two analyses share a question, and uses the same functions there: the shared `variant_annot.py` and
`gene_annotate.R`, the Open Targets query, and the ToMMo, 1000 Genomes and gnomAD lookups.

## 0. Upstream guard (`gene_signals.nf`, before any submission)

- **Upstream manifest.** `--gene_dir` must hold an `assoc_gene` manifest with these keys: cohorts, minAC,
  strata, min_num_var, covariates, methods, alpha, snpeff_index.
- **Parameters.**
  - `--Cohorts` must follow the upstream order.
  - `--Strata` must be the upstream variant sets.
  - The statistics must be `cmc` and `skato`.
  - minAC, the minimum number of variants and the number of PCs are read from the manifest.
- **Inputs.** Every input must exist before anything is submitted:
  - tuning.rv call sets and QC files;
  - the rare/common split and the controls' frequency table;
  - covariates;
  - per-cohort gene maps.

## 1. Gene selection (`select_genes.py`)

- **Genes.** Every gene with a consistency tier ≤ `MaxTier` (3) in either variant set, plus every gene below
  the Bonferroni threshold in any cohort × set × statistic. Class is the best tier; a gene that is only
  Bonferroni-significant is "Single cohort".
- **Variants.** For each gene, the union over cohorts and sets of the variants each cohort tested it on
  (`variant_stats_all`).
- **Asserts.**
  - The counts equal the gene index.
  - MODERATE / HIGH membership is identical in both sets.
  - No LOW variant is in the moderate/high set.

## 1b. Primary and sensitivity gene families (`gs_common.family_calls`, `tiers_from_calls`)

- **Primary.** The upstream family: the genes with ≥ `min_num_var` (2) qualifying variants. BH q, BH calls,
  Bonferroni calls and tiers are recomputed from `gene_scan_all` and asserted equal to the upstream flags,
  q-values and `tiers.tsv`.
- **Sensitivity** (`SensMinNumVar`, 3). The same tests and P values, with the multiple-testing family
  restricted to the genes with ≥ 3 qualifying variants. A gene's test does not depend on the other genes, so
  only the family changes:
  - BH q over the family;
  - Bonferroni threshold alpha / family size;
  - calls and tiers by the same rules.

  Genes with exactly 2 variants leave the family. Because BH depends on the ranks, a smaller family can add
  calls and can also remove them.
- **Selection.** The union of both families' selections. A gene selected only by the sensitivity family has
  class S ("Sensitivity only").
- **λ<sub>GC</sub> of the sensitivity family.** Equal to the upstream λ over the genes with ≥ 3 variants; this is
  asserted.

## 2. Rare-variant QC table (`rv_qc_stats.py`, per cohort)

The numbers of the project's rare-variant QC summary are recomputed from the tuning.rv outputs.

- **Input.** Samples by case / control come from `sample_metrics`; variants, from the Δ table after sample QC.
- **Sample QC.** The outlier rule is replayed: Theil–Sen line of the burden rate on mean depth, robust Z of
  the residual, and |Z| above the run's threshold (read from the run log). The replay must reproduce the
  removal list. The kept samples must be the call set's `.fam`.
- **Variant QC.**
  - Δ = |VMISS 30× − VMISS 15×|. The count above the knee must equal the removal count plus the protected
    rescues.
  - Monomorphic variants (MinAC 0 after the sample removal) are counted, together with how many of them the
    Δ filter removed anyway.
  - The sweep at κ = 0 must equal the variant-QC output.
  - The sweep at κ = 1 must equal that output minus the monomorphic variants not already removed.
- **Allele count.** The minAC recommendation and its non-significant band are taken from tuning.rv.
- **Burden gap.** The case–control burden-rate gap is computed at κ = 1 and at the chosen κ under two models:
  - the decision model, r ~ case + z(mean depth), whose replay must reproduce the sweep;
  - the association covariates, r ~ case + sex + the PCs.

## 3. Variant sets (`variant_sets.py`)

This step reads each cohort's upstream gene map (`01.map`).

- **Analysed variants.** Each variant is counted once, by its most severe snpEff annotation over all genes
  (`anno_counts`). The sum must equal the analysed variants of step 2 and `map_qc` TOTAL.
- **Assigned to genes** (`map_qc` TOTAL):
  - variants with a high-, moderate- or low-impact consequence in a named gene;
  - variant–gene rows;
  - variants assigned to more than one gene.

  The difference to the variants whose most severe impact is high, moderate or low is the set of such
  consequences that name no gene.
- **Per variant set** (the filtered maps and denominators; the pre-filter counts come from step 3b):
  - genes with a qualifying variant, not tested (fewer than `min_num_var`), and tested;
  - variants of the tested genes. Each is counted once, by its most severe consequence over the tested genes
    it is in, both by impact and by consequence term.
- **Asserts.**
  - Map rows equal the gene index.
  - Distinct variants equal the denominator.
  - Rows per impact equal the gene index's `n_high`, `n_moderate` and `n_low`.
  - Terms sum to the totals.

## 3b. The gene map before the gene filter (`map_recount.py`, per cohort)

- **Why.** The published map holds only the genes with at least `min_num_var` variants. The variants that
  qualify for a set but sit only in untested genes, and the genes with 1, 2 or ≥ 3 qualifying variants, are
  therefore not stored upstream.
- **How.** The map is rebuilt with the upstream script itself (`assoc_gene/scripts/build_gene_map.py`,
  staged):
  - on the cohort's analysed variants (`plink2 --mac <minAC> --write-snplist`, as upstream);
  - per chromosome;
  - with the impacts of the widest set.
- **Asserts.** The rebuild must equal the published map:
  - the summed annotation counts (every impact × term × biotype);
  - map QC TOTAL;
  - per set, the genes tested, the genes dropped, the variants tested, the variant–gene pairs, and the genes
    with 2 and with ≥ 3 variants.
- **Gives, per set:**
  - qualifying variants (a consequence of the set's impacts in at least one gene);
  - those only in untested genes;
  - genes with 1, 2 and ≥ 3 qualifying variants.

## 4. Genotypes, reasons, carriers (`genotype_stats.py`, per cohort)

- **Genotypes.** plink2 `--export A` on the cohort call set the tests used. Per group:
  - REF/REF, REF/ALT, ALT/ALT and no-call counts;
  - ALT frequency (EAF; EA = ALT) and missing share;
  - HWE P, from the shared `group_stats` (plink2 `--hardy`, the function behind the common-variant table).
    Its genotype counts are asserted equal to the export.
- **Minor allele and MAC.** Pooled, decided as the tests do.
- **Why a variant was not used:**

  | code | meaning | established from |
  |---|---|---|
  | `impact` | low-impact variant, outside the moderate/high set | the gene map |
  | `mac` | in the call set, MAC < the cohort's minAC | the genotypes |
  | `gap_filter` | removed by the depth-gap filter | the exclude list |
  | `common` | control MAF ≥ `RareMaf` in this cohort | the split's list; the MAF is read from the split's control frequencies and asserted ≥ the cut |
  | `no_carrier` | absent for none of the above reasons | asserted in step 7: none of its carriers elsewhere is a sample of this cohort |

- **Cross-check.** Every tested variant's counts are asserted equal to `variant_stats_all`.
- **Carriers.** Carriers per gene × set (in total, by case platform, by control source) are asserted equal to
  `set_stats_all`: carriers, MAC, multi-site carriers and number of variants.
- **Multi-site carriers.** Each sample with minor alleles at ≥ 2 of a gene's tested variants is listed with
  the variants and dosages (`multisite.<cohort>.tsv`). One row per counted multi-site carrier is asserted,
  and every row is asserted equal to the upstream `multi_carriers_all.tsv`. ASSEMBLE then merges each
  sample's rows into one per gene (the variants carried, and the sets and cohorts counting it) and adds the
  sample information:
  - platform, project and depths from the sample sheet;
  - sex, and for cases the age at DNA collection, age at diagnosis, DVT and max mPAP,
    from the project's sample table;
  - control ages from the AGP3K age / sex release.

  The sheet's platform and project are asserted equal to those recorded with the genotypes.

## 5. Variant annotation (as the common-variant table)

- **snpEff** (`annotate_function.py`), from the index the gene map was built from:
  - the tested gene's most severe consequence and impact, asserted equal to the map's;
  - its most severe transcript (protein-coding first), giving the transcript and HGVS;
  - other consequences, every other term over all transcripts and features at the variant.
- **Ensembl 86 gene.** The overlapping gene, else the nearest (`gene_annotate.R` via
  `variant_annot.lookup_genes`).
- **rsID.** The ToMMo record's ID. Where it lists several accessions, the one dbSNP 151 (ANNOVAR
  `avsnp151`) assigns to this allele. For variants absent from ToMMo, dbSNP 151's for this allele.
- **Frequencies.**
  - ToMMo 61KJPN: normalised lookup, 60KJPN cross-check.
  - 1000 Genomes 30x: EAS (504) and JPT (104) from genotypes, with the AF_EAS_unrel and LD-panel
    cross-checks.
  - gnomAD v4.1 genomes: AF, EAS and NFE (ANNOVAR). ANNOVAR's coordinates are re-derived and asserted.
- **Open Targets** (`opentargets.py`, unchanged from the common-variant table; a local task). Per variant:
  - ClinVar evidence;
  - GWAS and molQTL credible sets (molQTL classed cis / trans by distance to the target TSS);
  - enhancer-to-gene predictions.

  The summary cell gives the number of records and the top `OtTopN` (3). The second sheet lists every
  record.

## 5b. Leave-one-variant-out (`lovo.py`, per cohort × variant set)

- **Question.** Which of a gene's variants carry its gene-based signal. Each selected gene is tested again in
  each cohort and variant set with each of its tested variants left out in turn, by both statistics.
  Descriptive: no multiple-testing correction, and the P are not compared with the scan thresholds.
- **The upstream model, reproduced.** The inputs and options are those of `assoc_gene` (EMIT_TEST_INPUTS and
  RVTEST):
  - the cohort call set restricted to the set's variant list (`01.map/variants.<set>.txt`) and exported to
    the test VCF with the same plink2 build (`--export vcf id-paste=iid bgz`);
  - the published phenotype and covariate files (`00.prep`; sex and 10 PCs);
  - rvtest with `--setFile`, `--impute mean`, `--burden cmc,cmcWald` and `--kernel skato`;
  - the same gates: every covariate present; the case and control counts rvtest loaded equal the recorded
    ones; mean imputation reported; SKAT-O weights Beta(1, 25).
- **Sets.** Per gene, one set with all its tested variants (the baseline) and one set per variant without it.
  The set file is written as upstream writes it: one 1-bp range per variant, with the VCF's contig string.
- **Asserts.**
  - Every set is returned once, on every sample.
  - rvtest's variant count of each set equals its ranges, so a left-out variant is really out. For the
    selected genes the upstream map identity is `ok` in every cohort and set: no span overlap, no position
    shared with another gene. Two variants of one gene at one position would make removal impossible; this
    aborts.
  - The baseline equals the published result: P of both statistics, SKAT-O ρ, cmcWald β and SE, and the
    variant counts.
  - In ASSEMBLE, the variants left out per gene × set × cohort are exactly the variants the cohort's test
    used (step 7), and the baseline P equals the gene table's.
- **Derived (ASSEMBLE).**
  - Δ = log10(P without the variant) − log10(P with all variants): the orders of magnitude by which P rises
    when the variant is left out.
  - The top contributor per gene × set × cohort × statistic: the largest positive Δ (long table only; the sheet
    shades Δ by magnitude instead, so several contributing variants show alike).
  - Two-variant tests are flagged: each test without one variant is a single-variant test.
- **Primary variant set of a gene.** The set(s) with a consistency tier. For a gene without one, the set of
  its Bonferroni-significant result (class Single cohort) or of its sensitivity tier (class Sensitivity only).
  The sheet shows both sets of every gene side by side and marks the tiered ones (Tier column); `primary`
  stays in the long table.

## 6. Gene × set × cohort (`assemble_table.py`)

- **Burden.** The burden the test saw (carriers, MAC, variants per group) equals the counts of step 4.
- **Calls.**
  - `significant_bh` ⇔ q < alpha.
  - `significant_bonferroni` ⇔ P < alpha / genes tested.
- **Tier.** Recomputed from the calls, and equal to `tiers.tsv`.

## 7. Variant × cohort and the QC links (`assemble_table.py`)

- **Variants used.** Per gene × set × cohort, the variants marked used equal `n_var`, and their MACs sum to
  the gene's cumulative MAC.
- **Absent variants.** A variant absent from a cohort for no listed reason is shown to have no carrier there.
- **Frequencies and HWE P** lie in [0, 1].
- **QC links.**
  - The QC's kept cases and controls equal the upstream cohort samples.
  - Its analysed variants equal the variant-set totals.
  - Its minAC equals the manifest.
  - The variant-set genes and variants equal the denominators.

## 8. The workbook (`write_table.py`)

- **Layout.** XlsxWriter, in the style of the common-variant table. Presentation only: no interpretation is
  written. Every "—" has its reason in the sheet's footnotes.
- **Sheet Rare-variant QC.** Reproduces the project's QC summary table: same layout, rules, colours and
  footnotes; Arial for the table, DengXian for the notes.
- **Re-opened after writing.**
  - Sheet order is checked.
  - Gene order is checked on Gene-based signals and Gene-based variants; on Leave-one-variant-out, every
    gene row and its All-variants row.
  - The Open Targets variant rows must be the variants sheet's rows, in order, with every record at outline
    level 2.
  - Every header letter has its footnote, and every "—" its line under Missing values.
  - There is no `nan` / `None` / `…` or internal name; sample IDs and project names appear only in sheet
    Multi-site carriers (checked against the real sample IDs on every other sheet).
  - There are no inline strings and no unpreserved edge whitespace.

## 9. Deliberately not done

- No test of the scans is re-run, and no P value is corrected (no genomic control). λ<sub>GC</sub> is shown as
  observed. The only tests run here are the leave-one-variant-out re-tests of the selected genes (5b), whose
  baselines equal the published tests.
- No covariate is added. The association model is the upstream one: sex + 10 BBJ-mainland PCs.
- Sample IDs reach the workbook only in sheet Multi-site carriers, at the project's request, with their sample
  information. They are absent from the variant long table. Besides `gene_signals.multisite.tsv`, only `01.genotypes/geno.*.tsv` and
  `samples.*.tsv` hold IDs; they are analyst files.
