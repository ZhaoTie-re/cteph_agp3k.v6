# gene_signals — study notes

## Design decisions

- **Follow the common-variant table.** Variant annotation and the other variant-level information follow
  `saige_signals.xlsx`: its columns, its Open Targets summary and sheet, and its frequency panels. Content
  specific to the rare-variant analysis is added only where the gene-based test needs it:
  - per cohort, whether the test used the variant, and why not;
  - the variant sets by impact and consequence term;
  - the rare-variant QC table.
- **No interpretation in the workbook.** It reports what was observed; reading the signals is left to the
  analysts.
- **Three cohorts side by side.** Every gene is shown in narrow, intermediate and full; none is judged on
  full alone. The nested cohorts share most samples, so agreement across them is not replication.
- **No genomic-control correction.** λ<sub>GC</sub> is reported as observed.
- **The association model is the upstream one** (sex + 10 BBJ-mainland PCs, no depth).
- **A reason for every missing value.** Each "—" is explained in the sheet's footnotes. Each reason a
  cohort's test did not use a variant is established from data.

## What the build of 2026-10-13 showed

Upstream: the assoc_gene run of 2026-10 (gene_scan_all md5 `36f2aa84…`) and the tuning.rv run of 2026-10-08.

### Selection

There are 14 genes and 118 gene × variant pairs (6 HIGH, 60 MODERATE, 52 LOW).

| class | genes |
|---|---|
| Tier 1 | STBD1, OSBP |
| Tier 2 | RPH3A, JAK2 |
| Tier 3 | PIK3C2G, CD248, SCAMP3, NKD2, PES1, SCFD1, P2RY2 |
| Single cohort | FAM217B, TBC1D14, FAS |

### Rare-variant QC

The sheet reproduces the QC summary of 2026-10-08 exactly. For narrow:

| step | samples | variants |
|---|---|---|
| input | 2,168 (419/1,749) | 16,103,343 |
| after sample QC | 2,165 (418/1,747) | 16,058,030 |
| after variant QC | | 15,939,305 |
| analysed (κ = 2) | | 7,431,380 |

### Variant sets

Analysed variants by most severe impact:

| cohort | analysed | HIGH | MODERATE | LOW |
|---|---|---|---|---|
| narrow | 7,431,380 | 3,432 | 61,467 | 83,569 |
| intermediate | 8,050,528 | 3,755 | 67,044 | 90,602 |
| full | 9,161,666 | 4,413 | 77,396 | 103,231 |

- **Assigned to a gene.** 145,225 / 157,886 / 181,090. A further 3,243 / 3,515 / 3,950 high-, moderate- or
  low-impact consequences name no gene, mostly TF binding sites.
- **What each impact contains.**
  - HIGH: `structural_interaction_variant` (1,138 in narrow) is as frequent as `stop_gained` (1,123), and
    `protein_protein_contact` adds 29.
  - MODERATE: `sequence_feature` holds 5,267 variants, beside 54,672 `missense_variant`.
  - LOW: the largest term is `sequence_feature` (36,423), then `synonymous_variant` (33,598).
- **Tested genes, moderate/high-impact set** (N / I / F): 11,192 / 11,720 / 12,554 genes and
  61,153 / 67,257 / 78,558 variants.

### Variants of the table

- **Used in each cohort's test:** 97 (N), 102 (I), 117 (F) of 118.
- **Not used:**
  - 9 cohort × variant cases with no carrier in the cohort;
  - 7 where the variant is common in that cohort (control MAF 0.0101–0.0141). One example is PES1
    p.Thr229Met, which is common in full (0.0141) and used in narrow and intermediate;
  - the rest with MAC below 2.
- **rsID.** 113 from ToMMo, 3 from dbSNP 151 (absent from ToMMo), 2 with none.
- **Open Targets** (release 26.09). 9 of the 118 variants are in the variant index; 526 records (ClinVar
  8 variants, GWAS 2, E2G 3, molQTL 0).

## History

- **2026-10-09** — First version, with a multi-sheet summary, allele balance and carrier composition.
- **2026-10-13** — Rebuilt on request:
  - variant annotation and its Open Targets sheet follow the common-variant table;
  - added the variant-set and consequence-term sheets (impact composition per cohort);
  - the rare-variant QC sheet reproduces the project's QC summary;
  - interpretation removed;
  - the allele-balance step, CADD / REVEL, the ANNOVAR ClinVar columns, and the Summary, Carrier
    composition, All tested genes and Methods sheets dropped.

  Outputs no longer produced were moved to `results/_superseded_2026-10-09/`.
- **2026-10-13, later the same day.**
  - Sheets reordered to follow the analysis: QC → variant sets → consequence terms → scans → signals →
    variants → Open Targets.
  - Variant sets: rewritten around the mapping to genes. It now states:
    - why a variant can carry HIGH / MODERATE / LOW impact yet stay unassigned (TF binding-site motifs
      belong to no gene);
    - why a variant can be tested in several genes (readthrough transcripts such as ARHGAP8 /
      PRR5-ARHGAP8; the UGT1A and protocadherin clusters, up to 22 genes per variant).
  - Consequence terms: a total row per impact and an overall total, asserted equal to Variant sets.
  - Gene-based scans: one row per cohort × set, with the two statistics side by side. The Bonferroni
    threshold (shared) and the BH cut-off (per statistic) can no longer be misaligned; "–" marks a scan
    with no BH-significant gene.

  The other four sheets are unchanged.
- **2026-10-13, evening.**
  - Variant sets split into two tables: A follows the variants, B the genes. Variant and gene counts are no
    longer mixed in one list.
  - The counts the published map does not store are now exact. They come from the map rebuilt before the
    gene filter (`map_recount.py`, running the upstream `build_gene_map.py` and asserted equal to the
    published map):
    - the qualifying variants per set;
    - those only in untested genes;
    - the genes with 1, 2 and ≥ 3 variants.
  - Gene-based scans: λ<sub>GC</sub> is shown over the tested genes with ≥ 2 variants (all tested genes,
    the minimum of the tests) and with ≥ 3 variants, with the number of genes behind each. The upstream λ
    values are recomputed in ASSEMBLE and asserted equal (|Δ| < 5 × 10⁻⁴).
- **2026-10-13, evening, later.**
  - Variant sets reorganised: each step on one row, with its variants and the genes they fall in side by side,
    so the qualifying → not tested / tested split reads the same way in variants and genes.
  - New sheet Multi-site carriers. It lists the 10 people behind the multi-site counts (14 gene × set ×
    carrier rows), with anonymous labels, the variants they carry, and the cohorts counting them. All 30
    cohort records equal the upstream `multi_carriers_all.tsv`.
- **2026-10-13, sensitivity analysis.** The primary analysis stays upstream's gene family (≥ 2 qualifying
  variants). A sensitivity analysis judges the same P values over the genes with ≥ 3 variants
  (`SensMinNumVar 3`). In the code, the primary family's BH q, calls, Bonferroni and tiers are recomputed
  and asserted equal to upstream; the sensitivity family uses the same code.
  - **Family size and threshold.** For example, narrow M/H: 11,192 → 8,372 genes, Bonferroni
    4.47 → 5.97 × 10⁻⁶.
  - **λ<sub>GC</sub>.**
    - CMC: 1.163–1.189 (≥ 2) → 1.096–1.151 (≥ 3).
    - SKAT-O: 1.172–1.243 → 1.169–1.231.
  - **Tiers.**
    - OSBP (M/H, Tier 1) and P2RY2 (M/H, Tier 3) leave the family: they have 2 variants.
    - NKD2: M/H Tier 3 → 2, and L/M/H enters at Tier 3.
    - EFCAB6 (L/M/H) enters at Tier 3. It is added to the workbook as class S, "Sensitivity only".
    - Unchanged: STBD1 (Tier 1), JAK2 and RPH3A (Tier 2), and the other Tier 3 genes.
  - **Bonferroni genes.** OSBP leaves. PIK3C2G (narrow M/H SKAT-O) newly passes.
  - **Ranks.** BH depends on ranks. In intermediate M/H, STBD1's CMC q is 0.042 with OSBP in the family and
    0.064 without it.
  - **Workbook size.** 15 genes, 141 variants; Multi-site carriers now has 13 people (48 records equal to
    upstream).
- **2026-10-13, Multi-site carriers by sample.** At the project's request, the sheet now lists the samples by ID
  instead of anonymous labels, one row per gene × sample, under a header row per gene. It shows:
  - the variants carried, coloured by impact;
  - the sets and cohorts counting each sample;
  - project, platform (verbatim), target and observed depth;
  - sex and age at DNA collection;
  - for cases, age at diagnosis, DVT and max mPAP.

  With EFCAB6 added, there are 13 samples (4 cases, 9 controls). This is the only sheet with sample IDs; the
  checks enforce that.
- **2026-10-13, Gene-based signals layout.**
  - The two variant sets are stacked under each gene (six rows: M/H then L/M/H, each with the three cohorts)
    instead of side by side, which halves the width (42 → 21 columns).
  - Paired numbers share a cell: MAC, MAF and multi-site carriers as case / control.
  - BH-significant P and q are shaded as well as bold, so a call recurring across the three cohorts reads
    as a shaded run.
- **2026-10-13, leave-one-variant-out.** A new step (`LOVO`, 6 tasks, cohort × variant set) and sheet. Every
  selected gene is tested again in each cohort and set with each tested variant left out in turn, by both
  statistics.
  - **Scope.** 3 cohorts × 2 statistics × both variant sets: 1,198 tests. Descriptive only; no BH or
    Bonferroni.
  - **Exactness.** The 180 baselines (gene × set × cohort × statistic) equal the published P, ρ, β and SE. For
    the selected genes the upstream map identity is `ok` in every cell, so leaving a range out removes exactly
    one variant. This is asserted from rvtest's variant count.
  - **Primary variant set.** The set(s) with a tier; for class Single cohort the Bonferroni set; for class
    Sensitivity only the sensitivity-tier set:
    - M/H: STBD1, OSBP, JAK2, NKD2, SCFD1, P2RY2, FAM217B;
    - L/M/H: CD248, SCAMP3, TBC1D14, FAS, EFCAB6;
    - both: RPH3A, PIK3C2G, PES1.

    The other set is shown below as supplementary.
  - **Largest Δ in the primary set** (6 cells = 3 cohorts × 2 statistics):

    | gene (set) | variant with the largest Δ | cells | Δ |
    |---|---|---|---|
    | STBD1 (M/H) | p.Trp212Cys | 6/6 | 2.45–3.90 |
    | OSBP (M/H) | p.Arg462Gln | 6/6 | 2.34–3.20 |
    | RPH3A (M/H; L/M/H) | p.Arg544Cys | 6/6; 6/6 | 2.97–6.08; 2.46–6.01 |
    | JAK2 (M/H) | p.Val617Phe (SKAT-O); p.His608Tyr (CMC, whose P is 0.17–0.34) | 3/6; 3/6 | 4.14–4.97; 0.13–0.17 |
    | PIK3C2G (M/H; L/M/H) | p.Phe186Leu | 6/6; 6/6 | 1.76–4.87; 1.65–4.87 |
    | CD248 (L/M/H) | p.Asp152Asp | 6/6 | 1.39–5.49 |
    | SCAMP3 (L/M/H) | p.Val215Val | 6/6 | 2.50–3.71 |
    | NKD2 (M/H) | p.Gly104Arg | 6/6 | 2.96–4.26 |
    | PES1 (M/H; L/M/H) | p.Thr229Met (N, I; common in F, not tested there) | 4/6; 4/6 | 4.07–4.97; 3.33–4.74 |
    | SCFD1 (M/H) | p.Ile122Val | 6/6 | 0.46–4.09 |
    | P2RY2 (M/H) | p.Gly337Asp; p.Lys373Glu | 4/6; 2/6 | 1.67–2.23; 2.44 |
    | FAM217B (M/H) | p.Gly321Ser | 6/6 | 1.65–4.97 |
    | TBC1D14 (L/M/H) | p.Asn372Asn | 6/6 | 2.24–4.81 |
    | FAS (L/M/H) | p.Gly56Asp | 6/6 | 1.93–2.74 |
    | EFCAB6 (L/M/H) | p.Gln1340fs | 6/6 | 1.46–1.68 |

  - **Two-variant tests.** OSBP, CD248 and FAS (M/H) in all three cohorts, and P2RY2 and TBC1D14 (M/H) in N and
    I. Each test without one variant is a single-variant test; the sheet marks them ‡.
- **2026-10-13, Leave-one-variant-out sheet redesigned** (on request: the first layout was hard to read).
  - Each gene shows its tiered (primary) set(s). The other set sits directly under the gene, in grouped rows
    collapsed under one summary line, instead of a separate section at the end. The separate section had made
    STBD1's five LOW variants look missing.
  - Δ is the visual focus: its own block, shaded by magnitude (≥ 2; 1 to < 2; ≤ −1). This replaces the mark on
    the single largest Δ, which hid a second contributor (STBD1: c.220+1G>C beside p.Trp212Cys).
  - Carriers replace MAC. P after removal sits in a muted block of its own.
  - The cohort-coloured vertical rules and the full-width gene header rows are gone. The gene, its set with its
    tier and the variants tested are in the first column, one line per row.
- **2026-10-13, Leave-one-variant-out sheet, third layout** (the second was still hard to read). One quantity
  per cell: −log10 P, the same unit on the gene row (all variants, as published) and on the variant rows (that
  variant left out), so the drop reads down the column. Shading by the size of the fall (≥ 3 / 2–3 / 1–2; a rise
  ≥ 1 blue) is the only colour in the body. Variants are listed under their gene, sorted by the mean drop, with
  the consequence as text instead of impact colours. Only the tiered set(s) are shown; there are no grouped
  rows. The separate Δ and P blocks are gone: 12 columns instead of 18.
- **2026-10-13, Leave-one-variant-out sheet, both sets** (on request). Every gene is shown in both variant
  sets side by side (3 cohorts × 2 statistics × 2 sets). A Tier column per set marks the tiered ones (1–3; S3 =
  sensitivity tier only; – none), and a tiered set's −log10 P is bold on the gene row. The variant rows are the
  union of both sets; a low-impact variant leaves the M/H columns empty.
- **2026-10-13, sheet order and LOVO colours.**
  - Leave-one-variant-out and Gene-based Open Targets now follow Gene-based signals (… signals → LOVO → Open
    Targets → Multi-site carriers → Gene-based variants).
  - The LOVO shading is a diverging ColorBrewer RdBu scale in steps of 1 (a tenfold change in P). Red: −log10 P
    falls without the variant. Blue: it rises. Darker means 1–2, 2–3, ≥ 3. A colour key sits under the table.
  - The variant rows follow Gene-based variants (genomic order, asserted equal) instead of the mean drop, and
    the colour key moved to the title row, above the frozen panes, so it stays in view.
- **2026-10-13, sheet order.** Each sheet now follows the sheets it refers to: QC → Variant sets →
  Consequence terms → Gene-based scans → Gene-based signals → Gene-based variants → Leave-one-variant-out →
  Gene-based Open Targets → Multi-site carriers. LOVO and Open Targets both build on Gene-based variants: LOVO
  uses its variants and their order, and Open Targets lists the records for its rows. The sample-level sheet
  is last.
  - LOVO's Consequence now carries the full snpEff term as in Gene-based variants, and an Impact column was
    added in the same colours. Both are asserted equal to that sheet, variant by variant.
- **2026-10-13, notes and variant IDs.**
  - Gene-based variants and Leave-one-variant-out now show the Variant (GRCh38 CHROM:POS:REF:ALT) beside the
    rsID; Gene-based Open Targets already did.
  - Every sheet's notes were rewritten to one pattern: bold column names, a shared cohort definition, and a
    Missing values list with one line per placeholder and column. In LOVO, for example, a "—" in −log10 P
    points to the cohort's Carriers cell: grey count, MAC below 2; grey 0 / 0, no carrier; —, not among the
    rare variants. The Rare-variant QC sheet keeps the notes of the project's QC summary, which it reproduces.
