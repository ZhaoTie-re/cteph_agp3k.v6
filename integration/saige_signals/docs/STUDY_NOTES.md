# saige_signals — study notes

Numbers of the current build, from `analysis/assoc_saige/results` (P = `p.value.NA`, peak model `main`,
SuSiE L = 10, coverage 0.95, min |r| 0.5). They are the sanity anchors `verify.sh` §8 checks; a refactor
must not change them. Per-row reference values are frozen in `anchors.no_spa.tsv` (504 rows).

## Selection

- 16 loci, in order: FTO (GW), ARL13B, HS3ST3A1, LINC01146, GALNTL6, ARHGEF26-AS1, HLA-DQA1, INPP5D,
  MFSD11, CDH4, SLC44A2, PHF20L1, ABI3BP, ALPI, OR13C8, ADARB2 (SG).
- 168 variants = 16 leads + 3 alternative leads (ABI3BP `chr3:100830861:G:A`, MFSD11 `chr17:76738241:G:A`,
  SLC44A2 `chr19:10631611:C:T`) + 149 credible-set variants; 504 variant × cohort rows.
- Every locus has a peak in every cohort (48 pairs); every table variant lies in every cohort's window.

## Credible sets

- 48 refits reproduce the stored PIPs (max |ΔPIP| = 5e-16).
- 17 cohort × locus pairs have a set (one each): FTO 5 / 5 / 6, ARL13B 1 / 1 / 1, ARHGEF26-AS1 3 / 3 / 2,
  HLA-DQA1 intermediate 45, INPP5D full 22, HS3ST3A1 intermediate 15 / full 14, MFSD11 narrow 33 /
  intermediate 31, SLC44A2 intermediate 33 / full 32 (narrow / intermediate / full).
- 2 null fits (ADARB2 narrow, intermediate); 29 impure — 20 single-effect, 9 with all ten effects
  diffuse (INPP5D N/I, CDH4 N/I/F, PHF20L1 N/I/F, ALPI N).
- Status of the 504 rows: in CS 252 · below cut 9 · null fit 2 · impure 241. Two below-cut rows are
  near-ties (0.027984 < 0.028025; 0.006061 < 0.006084), hence four significant digits.

## Genotypes and references

- N cases / controls: narrow 419 / 1,776; intermediate 429 / 2,079; full 439 / 2,662.
- snpEff 168 / 168; Ensembl gene 168 / 168; no LoF / NMD prediction; one MODERATE (SLC44A2 set,
  `chr19:10631494:A:G`, p.Gln154Arg).
- ToMMo 167 / 168 (`chr6:32628445:T:A` absent; three matched only after normalisation;
  `chr2:233066247:C:T` FILTER HighlyMultiAllelic); 1KG EAS / JPT 168 / 168; gnomAD 168 / 168.
- r(ToMMo, AGP3K EAF) = 0.999; r(gnomAD EAS, 1KG EAS) = 0.999.
- The project EAS LD panel differs from the NYGC genotypes by one allele at `chr6:32628130:G:A`,
  `chr6:32631268:T:G` and `chr10:1747167:CGTG:C`, and lacks 8 of the 168 variants.

## Open Targets (release 26.09, queried 2026-09-30)

Not anchors: they change with the Open Targets release.

- 97 / 168 variants in the variant index; none of the 45 HLA-DQA1 variants (plausibly because Open
  Targets fine-mapping leaves out the MHC; not checked), nor the ADARB2, ARHGEF26-AS1, ARL13B, CDH4,
  GALNTL6, OR13C8 and PHF20L1 variants.
- GWAS credible sets for 52 variants: FTO 3 / 6 (L2G FTO), HS3ST3A1 14 / 16 (Dupuytren's disease; L2G
  HS3ST3A1), SLC44A2 33 / 33, ALPI and LINC01146 (L2G GALC) 1 / 1 each.
- Thrombosis-related GWAS sets only at SLC44A2: 9 variants, incl. the lead rs4548995 and the alternative
  lead rs1560711 (venous thromboembolism, pulmonary embolism, DVT; ~33–35 sets each; L2G SLC44A2 0.66).
- molQTL 79 variants, E2G 26 (38 variant × gene links), ClinVar 1 (MFSD11 `chr17:76744318:G:A`,
  "evidence only", cancer conditions); 5,847 detail records.
- molQTL cis / trans (1,638 records, 2026-10-05):
  - cis 1,629: eQTL 1,040, sQTL 329, tuQTL 239, sc-eQTL 11, pQTL 10. Median distance to the TSS 16 kb;
    the largest 912 kb (UNC13D eQTL at `chr17:76756425:G:A`).
  - trans 9: one UKB-PPP (EUR) pQTL for PDCD1 (chr2) at each of 9 SLC44A2 variants (chr19). Under their
    variants they rank 12th–24th by PP (0.018–0.035), so none is among the three shown in the main sheet.
  - None at 1–5 Mb. Open Targets' `isTransQtl` agrees with every record.
- Studies behind the 5,803 credible-set records: 554.
  - **GWAS, 310 studies**: 280 list cohorts (the 30 without are GWAS Catalog studies, 379 records); all 310 have a
    population structure. The most common are `NFE 100%` (188), `NFE 96%, EAS 2%, AFR 2%` (20), `FIN 100%` (19,
    FinnGen) and `EAS 100%` (19).
  - **molQTL, 244 studies**: none list cohorts. Only the 2 UKB-PPP studies (19 records) have a population (NFE
    100%); the eQTL Catalogue studies have none.

## History

Built first as a standalone script in `tmp/saige.signals.v6` (2026-09-16; values independently
re-derived from the raw sources by six verifiers), then ported to this component (2026-09-30) with the
same numbers. Same day: the workbook writer moved from openpyxl to XlsxWriter (Excel's repair prompt on
macOS), footnotes became short lines, Open Targets annotation was added, and the Python interpreter was
pinned (`source activate` had left another env's python3 first on PATH; same numbers either way).
2026-10-01: Open Targets summaries no longer truncated ("…"; one record per line in Signals), and the
Open Targets sheet's two biosample columns became "QTL tissue / cell type" and the "E2G" group. Then:
Signals' Open Targets cells back to one unwrapped line with up to 3 records by PP (no thrombosis-first
ranking), the CS footnote spelled out format by format, and the Open Targets sheet laid out as one row per
Signals variant (same order) with its records collapsed below and locus folding; a header row per locus
was then added so that Lead / Alt. lead rows own their records group (a lead row had been the summary row of
both its records and the locus group, so its records could not be opened on their own). Then: cohort labels
"Narrow / Intermediate / Full mainland"; title without "Table."; Gene without distance; gnomAD columns
widened; rsIDs attributed at the two multi-allelic positions: chr19:10670113:A:G → rs2360943 (dbSNP 151 A>G;
rs2067827755 is an insertion) and chr10:1747167:CGTG:C → rs1564211184 (rs1347066655 is the C>T SNV), both
confirmed on Ensembl (GRCh38). Then (2026-10-05): the seven frequency columns share one width;
the sheets are "Common-variant signals" and "Common-variant Open Targets"; and the credible-set P / β are the
variant's own where Open Targets reports them, else the set lead's (grey italic, source and lead in two new
columns). Of 5,803 credible-set records:

| statistic | variant's own | set lead's | none |
|---|---|---|---|
| P | 2,982 | 2,821 | 0 |
| β | 3,770 | 1,776 | 257 |

Before this change, every P and β shown was the lead's.

Later the same day, each molQTL record was classed cis / trans by the variant's distance to the target gene's TSS:
cis ≤ 1 Mb; trans = another chromosome or > 5 Mb; 1–5 Mb neither. Open Targets' set-level `isTransQtl` is
kept as a cross-check. The change adds:
- the `cis-` / `trans-` prefix and "(k trans)" in the main sheet;
- a "cis / trans (to target TSS)" column in sheet `Common-variant Open Targets`.

Every other value is unchanged.

Then sheet `Common-variant Open Targets` changed in two ways; its other cells and the main sheet are unchanged:
- The one-line description on row 3 became lettered notes a–f below the table, in the style of the main sheet,
  with letters on the headers. The header moved up to row 3.
- Two new columns, Cohorts and Population (% of sample).

The CS-column wording was then made plainer; the statuses and every number are unchanged:

| before | after |
|---|---|
| `Out (x < y)` | `Not in CS1 (PIP < y)` |
| `None: null fit` | `No CS: no signal` |
| `None: impure (n = …; min \|r\| …)` | `No CS: low LD (… variants; min \|r\| …)` |
| `None: impure` | `No CS: low LD` |
| `None: impure (10 effects)` | the strongest effect's set, e.g. CDH4 narrow `(592 variants; min \|r\| 0.00)` |

`dropped_set_*` is now filled for those 9 multi-effect fits as well.
