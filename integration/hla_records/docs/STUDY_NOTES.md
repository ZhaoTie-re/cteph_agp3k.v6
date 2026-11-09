# hla_records — study notes

## Design decisions

- **One workbook, in the order of the two decks.** Contents, the typing part (T1–T8), the association part
  (A1–A12), Sources and checks. Each sheet names the slides it records, and Contents lists every slide of
  both decks with its sheet or the reason it has none, read from the `.pptx` files at build time, so the map
  follows the decks.
- **The decks are the authority on definitions and numbers.** Counts are counts of tests (a mirror pair once);
  single markers use the random model with GEC (primary) and Bonferroni (sensitivity); positions use the
  fixed-model omnibus test with Bonferroni; the exact max-T is a benchmark only. What the decks dropped is not
  presented (multivariate-normal max-T, Benjamini–Hochberg, gene-level omnibus). P values are printed by a
  port of the decks' own formatter, including V8's inexact powers of ten, so a P reads as on the slide.
- **Two readers.** Content sheets are in the academic house style of the integration workbooks: no file,
  script or internal names and no sample identifiers. The last sheet lists the runs, the files read (with
  md5) and every check, for the analysts.
- **Everything per marker is kept.** Every allele and residue column per cohort (T7, T8) and every test
  of every comparison and cohort, both models side by side (A11, A12), sortable with P as numbers.
- **Every record is checked three times**: against the decks' data, by an independent recomputation from
  the per-sample and per-run files, and by re-reading the written workbook. A failure stops the run; the
  workbook is never written from records that did not pass.
- **Read only upstream.** The typing component's `verify.sh` writes nothing and is run as is; no process
  writes bytecode into the typing component's `scripts/`. The supplementary run given DQB1 57 V was made
  with the association pipeline's scripts unchanged, outside the component.

## What the build of 2026-10-28 showed

- **Upstream.** Association run 17, `verify.run17.log`: 814 checks passed, 0 failed. The typing component's
  `verify.sh`: 76 checks ok, 0 failed. The six association matrices match the typing summary's md5 sums.
- **Deck agreement.** 5,926 values compared over the 16 sheets that share numbers with a slide: 0
  mismatches. Every quantity of the summary slides (DQB1 57 V in both models, position 57, the GWAS lead
  given 57 V, what is left after conditioning, λGC) equals the deck's data.
- **Independent recheck.** 267,108 values recomputed: 0 mismatches. This includes 157,505 values of the
  9,844 single-marker rows against the 48 per-run tables, 21,313 of the 3,552 position rows, and every
  column's status, copies, call rates and MAF in every cohort (87,606 values) from the dosage matrices.
- **Audit.** 507,974 body cells as meant; 1,653 displayed P values in the decks' format; 312,139 numbers of
  the detail sheets equal to their tables; every note and every missing value explained; no internal names,
  no sample identifiers, clean XML.

## Found upstream (reported, not changed here)

1. **The association deck reads the run-16 verify log.** `collect_stats.py` reads `verify.run16.log`, so
   its `stats.json` gives 796 checks; run 17's log gives 814. The workbook reads the newest log.
2. **The association study notes describe run 16** (294 tasks, 796 checks); run 17 is the current run.
3. **Stale fixed-model files.** 24 pairs of `02.assoc/fixed/<class>/fixed.<class>.{json,tsv}` (12 runs × 2
   classes) date from 2026-10-22, before the final runs. No table of this component reads them.
4. **The typing README is stale.** It describes a 24-slide deck (now 30) and the earlier marker floor (163
   alleles and 697 residues under MAC ≥ 20; the floor applied is MAF ≥ 0.5%: 142–144 and 672–681).

## History

- 2026-10-27/28: component created; first build (above).
- 2026-10-28: sheet A5 shows the scans given the GWAS lead for the main comparison only, as slide 22; the
  DVT comparisons given the lead are not presented (the user's decision).
