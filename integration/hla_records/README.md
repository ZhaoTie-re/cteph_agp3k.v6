# `hla_records` — one workbook recording the HLA typing and the HLA association

An integrative component (`cteph_agp3k.v6/integration/`). It reads four inputs:

- the finished [`hla.typing.v2`](../../hla.typing.v2) run;
- the finished [`analysis/assoc_hla.v2`](../../analysis/assoc_hla.v2) run;
- the supplementary without-DVT run given DQB1 57 V ([`tmp/hla_dvtno_given57v`](../../tmp/hla_dvtno_given57v));
- the two decks and the data their builds read: `hla.typing.v2/deck` and `tmp/ppt/hla_assoc`.

From them it builds one workbook, `results/hla_records.xlsx`. The workbook follows the two decks: the same
order, definitions, terms and numbers. It omits what the decks dropped (the multivariate-normal max-T,
Benjamini–Hochberg, gene-level omnibus tests), anything per sample, and the per-draw bootstrap values. It
re-runs no test.

| sheet | deck · slides | one row per | content |
|---|---|---|---|
| Contents | both | sheet; slide | what each sheet records, and every slide of both decks with its sheet, or why it has none |
| T1 Typing run | typing · 11–12 | setting; group × platform | HLA-HD run; genomes, sex, HLA read pairs, depth, sample outliers by reason, genotype and high-confidence share per platform |
| T2 DRB3, DRB4, DRB5 | typing · 13 | gene | copies against DRB1 before and after reconciliation; two-allele calls |
| T3 Completeness | typing · 14, 16 | gene (19) | calls by result, basis of the allele, read evidence, confidence, locus tier |
| T4 ToMMo comparison | typing · 15, 17–19, 28 | step; gene × panel | harmonization; presence in ToMMo; r, D, excess D against ToMMo 61KJPN and 1KG JPT; DRB1\*15:13; homozygosity |
| T5 Marker QC | typing · 20–27 | step; gene; cohort × class | encoding; call-rate and MAF floors; stepwise record (143/144/142 allele and 672/681/679 residue columns); per gene |
| T6 P groups | typing · 29–30 | cohort; P group | tested columns sharing a P group; DQA1 through the comparison |
| T7 Alleles (detail) | typing · 25 | allele (1,995) | every allele column per cohort: status, copies, call rates, MAF; P group, panels |
| T8 Residues (detail) | typing · 26 | residue column (2,872) | every residue column per cohort |
| A1 Design | association · 4–7 | comparison × cohort; gene | samples of both models, history of DVT, tests per gene, models and settings |
| A2 Thresholds | association · 8, 11–14, 17–19, 39 | scan | Me, GEC and Bonferroni; exact max-T benchmark (score and Wald) with FWER; positions; Me recomputed |
| A3 Q1 Association | association · 20–22 | scan; marker | significant per scan (main comparison); the five strongest alleles, residues and positions |
| A4 Q2 Conditioning | association · 23 | step | stepwise conditioning; what is left given DQB1 57 V; B 114 |
| A5 Q3 GWAS lead | association · 2, 24 | cohort; scan; HLA marker | the lead and DQB1 57 V conditioned on each other; every scan of the main comparison given the lead |
| A6 Inflation | association · 25–26 | scan × conditioning | λGC before and after conditioning; genome-wide SNPs in and outside the MHC |
| A7 Summary | association · 27–29 | cohort; residue | DQB1 57 V in both models; copies per person; residues at 57 and their alleles; panels; Q1–Q3 |
| A8 DVT | association · 30–35 | comparison × cohort; marker | the three DVT comparisons; without DVT, the HLA-DQ signal and the run given DQB1 57 V |
| A9 Reported alleles | association · 36 | study; comparison × cohort | B\*52:01 and DPB1\*02:02: published, ours per copy and as carriers, power |
| A10 Appendix | association · 37–42 | scan; platform; P group | mirror pairs, model agreement, platform by group, position 57 by P group |
| A11 Single markers (detail) | association · 20–21 | test (9,844) | every allele and residue test, random and fixed model side by side |
| A12 Positions (detail) | association · 22 | position (3,552) | every whole-position test |
| Sources and checks | — | run; file; check | the upstream runs and their verification, every file read (md5), every check |

Table captions name their slides through the decks' own slide keys, so they follow the decks when slides move.
Content sheets carry no file, script or internal names and no sample identifiers; the last sheet lists the
sources and checks for the analysts. Every number shared with a slide is printed as the slide prints it.

## Checks of every record

1. **Deck agreement** (`DECK_CHECK`). Every quantity a sheet shares with a slide is compared with the data
   the deck's build read (typing `stats.json`; association `stats.json`, `data/results_all.json`,
   `data/pos57_pgroups.tsv`). Counts must be equal, floats agree to 1e-6. A mismatch fails the run; it is
   never resolved by copying the deck.
2. **Independent recheck** (`RECHECK`). The records are recomputed by another route, with code that imports
   nothing of the table scripts: call states and confidence from the per-sample tables; every column's
   status, copies, call rates and MAF from the dosage matrices; samples from the phenotype files; every
   single-marker and position row against the per-run tables; thresholds as 0.05 / Me and 0.05 / tests;
   significance counts and λGC from the detail rows; mirror pairs and the 57 V genotypes from the matrices;
   and every number that appears on two sheets.
3. **Audit** (`AUDIT`). The written workbook is re-read: every body cell, the P format, the notes behind
   every header letter and every missing value, no artefacts, no internal names, no sample identifiers, and
   the XML Excel is strict about.

The workbook is written only if checks 1 and 2 passed; the run manifest only if check 3 passed.

## Run

```bash
cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/integration/hla_records
./verify.sh                                                       # pre-flight
PATH=/home/b/b37974/anaconda3/envs/dsl2/bin:$PATH nextflow run hla_records.nf
./verify.sh                                                       # post-run: DAG, checks, anchors, workbook, decks
```

- **Jobs.** 11 short SLURM jobs on `gr10478b` and 1 local task (the run manifest). If the jobs sit
  pending with `QOSGrpCpuLimit`, add `--queue eb` (1 h limit; every task fits).
- **Upstream.** The guard refuses to submit unless the association run manifest names `assoc_hla.v2` and
  its newest `verify.run*.log` is clean and newer than the manifest. `TYPING_VERIFY` runs the typing
  component's own `verify.sh`, which writes nothing; every process sets `PYTHONDONTWRITEBYTECODE=1` so that
  importing the typing component's `config.py` leaves no `__pycache__` behind.
- **Resume.** `nextflow run hla_records.nf -resume <session>`, the session id from `nextflow log`. A new
  upstream run or a rebuilt deck changes the stamp every task takes, so all of them re-run.
- **After a deck changes.** Rerun: `verify.sh` §10 fails while the deck data differ from those the run
  checked.
- **Syntax check.** `nextflow lint hla_records.nf nextflow.config`. Do not use `nextflow run … --help`: it
  launches the pipeline.

## Files

```
hla_records.nf         the pipeline (upstream guard, 12 processes)
nextflow.config        SLURM resources (--rsc only), queue, trace / report
verify.sh              pre-flight and post-run checks (§1–§10)
assets/slide_map.tsv   every slide of both decks → its sheets, or why none
scripts/
  typing_tables.py     T1–T8 tables (imports the typing component's config and nomenclature as the authority)
  assoc_tables.py      A1–A12 tables
  gwas_lambda.py       λGC of the genome-wide scan in and outside the MHC, per cohort
  deck_index.py        number, section and title of every slide (python-pptx)
  deck_check.py        agreement with the decks' data
  recheck.py           the independent recomputation
  upstream.py          the upstream guard in detail; every file read with its md5
  write_workbook.py    the workbook (XlsxWriter)
  audit_workbook.py    the audit of the written workbook (openpyxl)
  deckfmt.py           the decks' number formats, ported from their build_slides.js
  hr_common.py         cohorts, comparisons, labels, TSV/JSON I/O
  hr_xlsx.py           house style of the integration workbooks (from gene_signals)
docs/METHODS.md        how every table is built and checked
docs/OUTPUTS.md        the results tree, the workbook, the tables
docs/STUDY_NOTES.md    design decisions, what the build showed, history
```
