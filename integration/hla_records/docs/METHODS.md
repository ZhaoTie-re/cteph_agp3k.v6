# hla_records — methods

The component builds one workbook from finished runs and the data of two decks. It runs no association test;
its only computations are summaries of upstream results (counts, frequencies, λGC, a GEC recomputation) and
the checks. Definitions follow the decks; where a deck defines a quantity, the deck's build script is the
authority and its code was ported, not reinterpreted.

## 0. Upstream guard (`hla_records.nf`, before any submission; `upstream.py`)

Before any job is submitted, the workflow requires: the association run manifest exists and names
`assoc_hla.v2` and the three cohorts; the newest `logs/verify.run*.log` reports `N passed, 0 failed` and is
newer than the manifest; the typing summary (`06.assoc_inputs/summary.json`) and the eight deck files exist.
A stamp of their names and modification times is an input of every task, so a new upstream run or a rebuilt
deck re-runs everything on `-resume`.

`TYPING_VERIFY` runs the typing component's own `verify.sh` (read only: it prints and writes nothing; the
process sets `PYTHONDONTWRITEBYTECODE=1`). `UPSTREAM` then fails unless that log has no `[FAIL]` line, the
md5 of the six association matrices equals the typing summary, and the association checks above hold. It
writes `upstream.json`: the runs and their verification, and every file the table scripts read, with its
sheets, modification time and md5 (sheet *Sources and checks*).

## 1. Typing tables (`typing_tables.py`; sheets T1–T8)

The typing component's own `config.py` and `nomenclature.py` are imported, so thresholds (call rate 95%,
MAF 0.5%, the 15% minor-share and 10-read rules, the outlier rules, tier criteria), gene lists, P-group keys
and release numbers are the typing component's, not copies. Definitions as the typing deck:

- **Calls** are genome × reported gene (3,101 × 19). *Genotype share*: heterozygous + homozygous + present +
  absent over calls. *High confidence*: no weakness in the read evidence.
- **Sample outliers** are listed, never removed: HLA read pairs below the 1st percentile, three or more
  reported genes without a two-field genotype, or a failed gene that every haplotype carries.
- **Platform** is `WGS_Platform` of the sample sheet, never regrouped.
- **DRB3/4/5** copies are compared with those the two DRB1 alleles imply, before and after the typing
  component's reconciliation.
- **ToMMo comparison**: controls only, our two-field alleles mapped to P-group keys as ToMMo names them,
  null alleles dropped on both sides; r over the keys in either table, D = total variation distance, excess
  D beyond the sampling expectation; presence of our keys: in ToMMo, absent, or named after ToMMo's newest
  allele (release 3.21.0). DRB3 and DRB4 are kept but not compared (ToMMo's tables are not comparable).
- **Marker QC**: per cohort, columns that vary; DRB3/4/5 left out; call rate ≥ 95% in cases and in
  controls; MAF ≥ 0.5%. Each column's status is the first rule it fails, in that order.
- **P groups** (T6): tested columns of one P group and their frequencies in the controls (the slide's
  quantity), the other alleles of the group, and ToMMo's entry.

Every count is asserted against the typing component's own record (`qc_record.tsv`, `markers.tsv`,
`summary.json`) before it is written.

## 2. Association tables (`assoc_tables.py`; sheets A1–A12)

Definitions as the association deck (`tmp/ppt/hla_assoc/collect_stats.py`, `make_figures.py`):

- **Tests**: a mirror pair (two residue columns that are complements in every sample) counts once
  (`test_id == marker`). Every count of the workbook is a count of tests.
- **Models**: alleles and residues use the random model (SAIGE, `p = p.value.NA`); the fixed model is a
  sensitivity analysis. Whole positions use the fixed-model omnibus likelihood-ratio test.
- **Thresholds**: single markers GEC = 0.05 / Me (primary) and Bonferroni = 0.05 / tests (sensitivity);
  positions Bonferroni over the positions tested. The exact max-T (parametric bootstrap, 10,000 null
  phenotypes) is the benchmark of the single-marker thresholds; it calls nothing. The pipeline's
  multivariate-normal max-T (`sig_maxt`, `p_maxt`) is not presented.
- **Significant in n/3**: the number of the three nested cohorts in which a test passes; robustness to the
  ancestry definition, not replication.
- **λGC** = median χ²₁.isf(P) / χ²₁.ppf(0.5), each P read as one degree of freedom, over tests with a
  defined P (`errcode == '.'`). Before and after conditioning are computed over the same tests.
- **Stepwise** (A4): the upstream walks. When the omnibus walk enters nothing at a step it records no unit;
  as the residue walks do, the step's strongest position (smallest P in any cohort) is then taken from that
  step's own scan, and asserted not to pass. The allele scan given DQB1 57 V comes from `HLA_GIVEN_SIGNAL`,
  as the deck's summary slide reads it.
- **GWAS lead** (A5): the lead of the genome-wide scan in the MHC kept by the upstream rule (P < 1 × 10⁻⁵ in
  all three cohorts), its P given HLA markers and the reverse, and every scan of the main comparison given
  the lead (as slide 22; the DVT comparisons given the lead are not presented).
- **Genome-wide λGC** (`gwas_lambda.py`, per cohort): the SAIGE scan of the same samples, SNPs with MAF ≥
  0.5%, inside chr6:29.6–33.4 Mb and outside.
- **DQB1 57 V** (A7): both models; copies per person from the residue matrix; each residue at position 57
  with the alleles carrying it (residues through `allele_match_map.tsv` and `residue_reference/DQB1.tsv`);
  valine in ToMMo and 1KG JPT.
- **DVT** (A8): the three DVT comparisons with the same tests and thresholds; the without-DVT HLA-DQ
  signal; the supplementary run given DQB1 57 V (every marginal P asserted equal to the scan's).
- **Reported alleles** (A9): the published values as printed; ours per copy (both models) and as carriers
  on the unrelated samples (Woolf CI, Fisher P); power at the published carrier OR and the OR detectable
  with 80% power (normal approximation on log OR, two-sided α 0.05).
- **GEC recomputation** (A2): Me of the all-cases full-mainland scans recomputed from the genotypes (plink2
  export of the unrelated samples), with GEC's blocks and published polynomial, duplicates |r| ≥ 0.9999 set
  aside.

## 3. The decks (`deck_index.py`, `assets/slide_map.tsv`)

`deck_index.py` reads both `.pptx` (python-pptx): per slide its number, section and title, tied to its key
through the `ORDER` array of the deck's build script. `slide_map.tsv` maps every key to its sheets, or gives
the reason it has none (cover, background, method illustrations, the typing deck's ten-person excerpt). The
writer refuses a slide without an entry, so the Contents sheet follows the decks.

## 4. Deck agreement (`deck_check.py`)

An explicit map, sheet by sheet, of every quantity shared with a slide to the value in the deck's data
(typing `stats.json`; association `stats.json`, `data/results_all.json`, `data/pos57_pgroups.tsv`). Counts
must be equal; floats agree to a relative 1e-6 (the tables hold seven significant digits), shares the typing
deck keeps at four digits to 1e-3, the genome-wide λGC to 1e-5 (its P values are read at another
precision); text must be equal. A mismatch fails the run. It is traced, never copied from the deck.

## 5. Independent recheck (`recheck.py`)

Code of its own (it imports none of the table scripts), other routes:

| sheet | recomputed from |
|---|---|
| T3 | call states and confidence counted from the per-sample genotype and confidence tables |
| T7, T8 | every column in every cohort: status, copies in cases and controls, call rates, MAF, from the dosage matrices of all genomes and each cohort's case status (MAF compared at six decimals, as the typing component prints it) |
| T5 | the stepwise record counted from those statuses |
| A1 | samples of both models counted from the phenotype files of each run |
| A11 | every row against the association tables of the 48 runs (OR, CI, P, SPA P, counts, frequencies) |
| A12 | every row against the omnibus table of each comparison |
| A2 | GEC = 0.05 / Me and Bonferroni = 0.05 / tests from each run's effective-test table |
| A3, A8 | per scan: tests, significant (primary, Bonferroni, SPA), smallest P and its marker, λGC (own χ² code) |
| A2, A10 | mirror pairs: columns of one position that are complements in every genome |
| A7 | copies of valine at DQB1 57 per person from the residue matrices |
| cross-sheet | a number on two sheets is the same: columns kept (T5 = A2), tests (A1 = A2 = A3 = A11/A12), DQB1 57 V (A7 = A8 = A11), position 57 (A7 = A12) |

## 6. The workbook (`write_workbook.py`)

XlsxWriter only. House style of the integration workbooks (`hr_xlsx.py`, copied from gene_signals):
booktabs rules, spanner / sub-spanner / leaf headers, cohort colours `#0571B0` / `#008837` / `#B2182B`, gene
symbols in italics, lettered notes, a *Missing values* block, an abbreviations line. Additions:

- **Numbers as the slides print them.** `deckfmt.py` ports the decks' formatters: JavaScript rounds the exact
  binary value half up (`toFixed`, `toPrecision`), and V8's `10 ** e` is one ulp off for eight exponents
  (−29, −24, −23, −21, −20, −17, −5, −4), which the decks' `pv` divides by; the port reproduces both and
  agrees with node on 20,000 random values. P: two significant digits from 0.001, `m.mm × 10⁻ⁿ` below.
  Detail sheets hold P as numbers for sorting.
- **Missing and not applicable.** `—` is a missing value and `–` not applicable; each sheet that shows one
  states why under *Missing values*.
- **Long text** spans merged columns and wraps, so it does not widen the columns of the other tables.
- A column without a single value stops the build (a wrong key, never data).

The writer refuses to write unless the deck agreement and the recheck passed. It also writes
`workbook.json` (sheets, blocks, rows, columns) and `workbook_cells.json.gz` (what every body cell shows).

## 7. Audit (`audit_workbook.py`)

openpyxl re-reads the written file: every body cell shows what the writer meant, at its row and column; on
the detail sheets every number and P equals its source table row by row; every displayed P has the decks'
format; every header letter has its note and every `—` / `–` its reason; no `nan`, `None`, `True`, `False`
or `…`; no internal names on the content sheets; none of the 3,101 sample identifiers on any sheet; no
inline strings, frozen pane without a split, overlapping merges or untrimmed shared strings. It records the
md5 of the workbook it audited; `verify.sh` §9 requires the published workbook to be that one.

## 8. Deliberately not done

- No test is re-run, and no deck is changed. A deck error is reported, not copied.
- Not presented, because the decks dropped them: the multivariate-normal max-T, Benjamini–Hochberg calls,
  gene-level omnibus tests.
- Nothing per sample (the typing deck's ten-person excerpt is not reproduced); no per-draw bootstrap values;
  no figures.
