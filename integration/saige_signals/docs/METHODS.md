# saige_signals — methods

## 1. Selection

Loci are the two cross-cohort panels of `assoc_saige` (`lead_crosscohort.tsv`, column `panel`):

- **GW** — lead P < `p_genomewide` in at least one cohort (panel b);
- **SG** — lead P < `p_suggestive` in every cohort and genome-wide in none (panel c).

How the upstream component defines them: per cohort, variants with P < `p_suggestive` within
±`peak_flank_bp` form one peak, whose lead is the smallest P; leads of different cohorts within
`peak_flank_bp` are one locus unless a cohort's stepwise conditional analysis lists them as independent
signals; the representative lead (`locus_rep = 1`) has the smallest P across cohorts. Loci are ordered
as the cross-cohort figure draws them (panel b first, then by that smallest P).

Rows per locus: the representative **Lead**, any **Alt lead** (another lead of the same locus under which
one cohort called its peak), and every **CS variant** — member of a SuSiE credible set of the locus's peak
in at least one cohort (union over cohorts). Each locus × cohort pair is mapped to the peak that cohort
called there, found through every lead of the locus.

## 2. Per-cohort statistics

Inside a cohort's fine-mapping window (the ±`peak_flank_bp` window of its peak) P, OR, CI and PIP come
from that peak's `pip.tsv` / `sumstat.tsv`. A variant outside the window, or a locus without a peak in the
cohort, takes P / OR / CI from the cohort's genome-wide `sumstats.tsv` and has no PIP ("not fine-mapped").
The genome-wide file is read for every variant, and the two sources are asserted equal wherever both
exist. The effect allele is always ALT (`A1 == ALT` asserted).

Genotype counts, effect-allele frequency, missingness and HWE (plink2 `--geno-counts --hardy`, exact
test) are computed per cases and controls on the genotype set the scan ran on (related samples retained),
through `_shared/scripts/variant_annot.group_stats` — the function behind the upstream lead tables — and
asserted equal to those tables for every lead × cohort.

## 3. Credible-set status

`susie_finemap.R` (upstream) calls `susie_rss(bhat = log OR, shat = SE, R = in-sample LD, n = analysed N,
L, estimate_residual_variance = FALSE)` and keeps the 95 % sets of `susie_get_cs(…, min_abs_corr)`.
`susieR` drops every set whose purity (smallest absolute pairwise LD) is below `min_abs_corr`, and
excludes effects whose prior variance V ≤ 1e-9; neither is recorded. `susie_cs_diagnostic.R` therefore
refits every locus × cohort with the stored inputs and the stored arguments (defaults taken from the
peak's `susie.json`), asserts that the stored PIPs (tolerance `SusieTol`), sets, coverage and purity are
reproduced, and lists every effect's 95 % set with exact purity (`susie_get_cs(min_abs_corr = 0,
n_purity = Inf)`). Each variant × cohort then gets one status:

| status | meaning | shown as |
|---|---|---|
| `in_cs` | member of set k | `CSk (n)` |
| `below_cut` | the cohort has a set; the variant's PIP is below y, the set's smallest member PIP | `Not in CSk (PIP < y)` |
| `no_cs_null` | no effect fitted: every V shrank to 0, all PIP = 0 | `No CS: no signal` |
| `no_cs_impure` | effect(s) fitted; every 95 % set failed the purity threshold | `No CS: low LD (n variants; min |r| r)` on lead rows, `No CS: low LD` on CS-variant rows |
| `not_fine_mapped` | no peak in that cohort, or outside its window | `Not fine-mapped` |
| `not_tested` | absent from that cohort's scan | `Not tested` |

Several credible sets per locus × cohort are supported throughout (`CS1`, `CS2`, …). A variant outside all
of them is `Not in CS1, CS2 (PIP < y1, y2)`, with every set's cut. "PIP < y" is asserted on both the PIP and the
α scale.

**The dropped set of a `no_cs_impure` fit** (`dropped_set_size`, `dropped_set_min_r`; the lead-row bracket) is
the 95% set of the fit's strongest effect, the one with the largest prior variance V:
- In a single-effect fit it is the only effect.
- In the 9 fits where all L = 10 effects keep V > 1e-9, it dominates the rest: V is 40–353× the next, and the
  other effects are near-null and spread over the whole window. Showing that one set replaces the former
  "(10 effects)", which read as ten signals.

## 4. Annotation

snpEff: the JHRP v6 pre-computed index; one annotation per transcript, the most severe impact shown
(protein-coding preferred on ties) with the other consequence terms listed; `sequence_feature` (LOW)
ranks above `intron_variant`, so LOW on an intronic variant is not a coding effect. Gene: Ensembl 86
gene overlapping the variant, else nearest within 1 Mb (protein-coding preferred; clone-named models and
pseudogenes excluded); the signed distance stays in the long table (`Gene_Distance_bp`, negative =
upstream of the gene start) but is not shown in the workbook. rsID: the ToMMo 60KJPN ID column matched on
chr:pos:REF:ALT (upstream). At a multi-allelic position that field lists every accession of the position, so
`ASSEMBLE` attributes it with allele-specific dbSNP 151 (ANNOVAR `hg38_avsnp151`, read through its index):
the candidate dbSNP lists for this exact allele, else the only candidate dbSNP does not list for another
allele at the position (`rsID_source`); otherwise it stays unresolved (‡).

## 5. Reference frequencies

All frequencies are of the effect allele. Reference records are split and left-normalised
(`bcftools norm -m -any -f hs38DH.fa`) before matching on chr:pos:REF:ALT, so multi-allelic and
un-trimmed records match.

| source | detail | cross-check |
|---|---|---|
| ToMMo 61KJPN (2025-06-16) | INFO/AF of the genotype-frequency VCFs; FILTER kept | 60KJPN normalised VCF, r > 0.99 |
| 1000 Genomes 30x GRCh38 | AC/AN of the NYGC genotypes subset to the 504 unrelated EAS of the project LD panel and the 104 JPT among them | plink2 `--freq` on the LD panel (within one allele); INFO/AF_EAS_unrel at biallelic sites |
| gnomAD v4.1 genomes | ANNOVAR `hg38_gnomad41_genome`, AF_eas and AF_nfe (4 dp) | ANNOVAR coordinates re-derived from REF/ALT |

The LD panel is not used for the reported values: it keeps only the last ALT at multi-allelic sites and
lacks some variants.

## 6. Open Targets

One GraphQL query per variant to the Open Targets Platform (`chr:pos:REF:ALT` → `chrom_pos_ref_alt`),
all pages of each list (page size 100; larger pages are refused as too expensive), 4 tries with
back-off; the API and data release are recorded once (`meta`). A variant absent from the Open Targets
variant index has no record of any kind there (the index is not exhaustive: common variants such as the
FTO CS variant rs16952623 and all 45 HLA-DQA1 variants are absent) — it is shown as "—", not as zero.

| column | source | summary in the table |
|---|---|---|
| ClinVar | `evidences(datasourceIds: ["eva"])` | number of records; first `OtTopN` clinical significance: condition, alphabetical |
| GWAS credible sets | `credibleSets(studyTypes: [gwas])` with the variant's posterior probability (`locus(variantIds)`) and the set's top L2G gene | number of sets (in brackets, those whose trait matches `OtFocusTraits`: counted, not ranked); first `OtTopN` distinct traits by PP |
| L2G | top `l2GPredictions` of those sets | the most frequent gene (ties: higher median score, then name), median score |
| molQTL credible sets | `credibleSets(studyTypes: [eqtl, pqtl, sqtl, tuqtl, sc*])`, with `isTransQtl` and the target gene's canonical transcript | number of sets (in brackets, those in trans); top `OtTopN` by PP as cis- / trans-type · gene (tissue / cell type, PP) |
| Enhancer-to-gene | `enhancerToGenes` | number of target genes; top `OtTopN` by the best score (best-scoring tissue / cell type, number of tissues / cell types) |

Names are never truncated: the summaries keep full trait, tissue / cell type and condition names (records
separated by `; `, a `;` inside a name becomes `,`). Sheet `Common-variant signals` shows each summary on one
unwrapped line (`OtTopN` = 3 records by PP, E2G by score); text beyond the column is clipped on screen and read
in the formula bar. Sheet `Common-variant Open Targets` ranks the records under each variant with the same keys
and stable sorts on the fetch order, so the records named in the main sheet are the first ones listed there.

**P and β of a credible-set record.** The set-level `pValueMantissa / pValueExponent / beta` of Open Targets
belong to the set's LEAD variant (`credibleSets.variant`), not to the annotated variant. Each statistic is
therefore taken, in this order:
1. the variant's own value in that study (`locus(variantIds)` row: `pValueMantissa / pValueExponent`,
   `beta`), or the set-level value when the variant is itself the lead;
2. otherwise the lead's value, shown in grey italic, with its source in "P / Beta from" and the lead in
   "Set lead variant";
3. otherwise missing.

PICS sets (fine-mapped from summary statistics without per-variant effects) carry only the lead's values.
SuSiE sets usually carry the variant's β, and sometimes its P. Study N, the fine-mapping method and L2G
belong to the study or the set and are labelled so.

**Cohorts and population of a study.** Both come from the study record:
- **Cohorts** = `cohorts`, as listed by Open Targets. Where it lists none, a molQTL record shows the study's
  `projectId` instead (grey italic), e.g. GTEx-v10, BLUEPRINT, UKB_PPP_EUR. A GWAS record shows —.
- **Population** = `ldPopulationStructure`: the share of the study sample assigned to each gnomAD reference
  population, which Open Targets uses for LD. It is written largest first, e.g. `NFE 96%, EAS 2%, AFR 2%`; a share
  below 0.5% is written `<1%`. Open Targets reports it for GWAS studies and UKB-PPP, not for eQTL Catalogue
  studies (—).

**cis / trans of a molQTL record.** The annotated variant is compared with the QTL's target gene. The
gene's transcription start site (TSS) is the start of its Ensembl canonical transcript on the + strand and
its end on the − strand (Open Targets target annotation, `canonicalTranscript`). Without a canonical
transcript the gene's own span is used, and with an unknown strand the nearer end. Each record is classed:
- **cis**: same chromosome, |variant − TSS| ≤ `--cis-window`. The default is 1 Mb, the usual cis window
  (GTEx, eQTL Catalogue).
- **trans**: another chromosome, or more than `--trans-distance` away. The default is 5 Mb, the threshold of
  Open Targets and eQTLGen.
- **neither**: in between; shown with its distance, no call made.
- **unknown**: no gene location.

The class is per variant. Open Targets' `isTransQtl` (gentropy `StudyLocus.flag_trans_qtls`) is per
credible set: the set is trans when every variant in it is on another chromosome or more than 5 Mb from the
TSS (null for GWAS sets). It is kept as a cross-check (`ot_is_trans`). A contradiction (cis with the flag
true, trans with it false) is written after the class in sheet `Common-variant Open Targets` and counted in
`opentargets_meta.json` (`qtl_contradictions`).

The main sheet prefixes each record shown with `cis-` / `trans-` and gives the number of trans sets in
brackets, because the trans records need not be among those shown.

E2G predictions come per gene × biosample (up to ~2,900 for a promoter variant); the detail sheet keeps
one row per target gene with the best score, its best-scoring tissue / cell type and the number of tissues
/ cell types predicting the link (columns prefixed "E2G", apart from the molQTL study's "QTL tissue / cell
type"). Every fetched record count is asserted against the API's `count`.

## 7. Workbook

Written with XlsxWriter (shared strings, rich strings, print titles columns-first as Excel writes them,
column page breaks spanning all rows). An openpyxl-written version triggered Excel's repair prompt on
macOS (print titles rows-first, a row-break extent on column breaks, whitespace-only rich-text runs
without `xml:space="preserve"`); `write_table.py` re-reads the file and asserts none of these is present.
Footnotes are one line per row in column A, wrapped to the width of the frozen columns, so they neither
cross the pane boundary nor depend on row-height estimates.

## 8. Checks

Every step asserts invariants of the data (never fixed counts): each selected variant once per cohort;
one representative lead per locus; at least one peak per locus; one peak maps to one locus; credible-set
rows equal the recorded set sizes; role counts sum; SuSiE fits = pairs with a peak, each reproducing the
stored fit; the status partition covers every row and every non-member row has a reason; recomputed lead
statistics equal `lead_crosscohort.tsv`; credible-set rows equal `cs_variants_all.tsv`;
r(ToMMo, AGP3K) > 0.98 and r(gnomAD EAS, 1KG EAS) > 0.95; every frequency in [0, 1]; one Open Targets
row per variant, and a variant absent from its index has no records. The written workbook is re-opened
and scanned (both sheets) for internal file / script / pipeline names and for the XML Excel repairs. `verify.sh` adds the site rules,
the DAG cardinality and, for a no-SPA run, the frozen anchors of `docs/anchors.no_spa.tsv`.
