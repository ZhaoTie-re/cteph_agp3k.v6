# Deck: is any CTEPH association signal a technical artifact?

Summary deck of `check/platform_depth` for the CTEPH research group and supervisors, for a live talk (speaker notes on
every slide) and for reading offline. Narrative, colors and visual grammar follow the HLA association deck
(`tmp/ppt/hla_assoc`): paper order, driven by four questions (Q1 depth, Q2 platform, Q3 technology, Q4 control source),
each answered by one experiment; no roadmap and no part dividers. Every reading is shown in the three nested cohorts
(narrow, intermediate, full mainland). Professional American English; DengXian, conclusion titles, cited rationale. One
message per main slide (about 60 body words, tables of at most six rows); the reasoning that does not fit is in the
speaker notes or the appendix.

**Deliverables:** `platform_depth_summary.pptx`, `platform_depth_summary.pdf` (LibreOffice; DengXian drawn as Liberation Sans).

## Structure (32 slides: 28 main, about 35 min; 4 appendix)

| Section | Slide keys (`ORDER` in `build_slides.js`) |
|---|---|
| Introduction | cover; background (why four checks: the sequencing groups to scale, beside a table of the four features of our data (depth, case platforms, no case on the controls' platform, five control sources), how each could create a false signal, and the question that tests it, each row citing the study it rests on; Decision / Rationale: test each feature, no platform covariate); aims (one aim sentence; a table with a column per question: question, schematic with hypothetical values, approach, and one row per signal type with ● or – and the limit; the order of the talk as stage boxes) |
| Data and design | samples (cohort table, row headers in the cohort colors with *k* components, rationale, cases by platform, controls by source; the signals re-tested (GWAS leads, DQB1:57:V, the Tier 1–2 genes in their M/H sets); decision box); models (each signal type's model, covariates, relatedness, test and threshold; decision box: re-tests reuse them unchanged) |
| Threats | threat (T0–T5: path → prediction → the question that tests it; T2 a calling error of one case platform, tested by the strata; T3 a DNBSeq-T7 calling error against HiSeq X (most of the cases), tested directly by the arrays; genotype filters are not a threat of their own, and an error shared by every case platform is not modeled, as NovaSeq shares the controls' chemistry) |
| Design and fidelity | workflow (the design at a glance and the result of the two fidelity checks: same data and analysis reproduce every result; reads re-processed with none removed leave genotypes and HLA types unchanged; the gates slide merged here in review round 8, its table and the ten pipeline gates G1–G10 in the notes and the appendix) |
| Q1–Q4 | each question's design slide followed at once by its results (user 2026-11-04): dsDesign, dsEmulation, dsSignals (Q1); platDesign (question, three readouts, two layers), looRule (three steps beside the schematic `loo_regions`: removing cases lowers z anyway → compare with z_exp, chance → judge only where a fall can tell: flag / below the bar / not testable), looPlat (by outcome: flagged with z, z_exp and chance; below the bar ordered by chance with a reading; everything else; gene sets read with the CMC burden β), platRes (first why two MAF layers, FTO vs INPP5D, then the five smallest Q P and the threshold basis) (Q2); arrayDesign (question; compact q3_flow; readout table "If DNBSeq-T7 calls made the signal"; control and blind spot), arrayMatch (one marker per signal with Minimac R² and the genotype concordance with WGS in cases and controls; drivers with the WGS case carriers found), hlahdMethod, hibagMethod, hlaCompare (concordance; HLA-HD − HIBAG not differential), arrayRes (GWAS leads and DQB1:57:V: OR per ALT and R²/concordance panels; verdicts N / I / F; CDH4 reproduced at the association level only; genome-wide agreement as the control), arrayGenes (driver case carriers confirmed on the imputed arrays, a description; why describe, not test) (Q3, user 2026-11-05); ctrlDesign (question; two readouts with their flag rule; every case stays; blind spot), looSrc (every signal without each source by signal type: readings / consistent / nt / flag, all flags 0; the RPH3A-without-Nagahama table; gene note), ctrlRes (Q P of all eight common markers, N / I / F, and the threshold basis) (Q4 focused on "no signal changes when the controls are stratified", user 2026-11-04; source-against-rest differences and per-source chances in the notes only); the joint leaving-out slide split by question (user 2026-11-04) |
| Signals | matrix (21 signals × 7 readouts × 3 cohorts: the registry's 7 checks with leaving out a case platform split into the DNBSeq-T7 cases and one small platform, Q and RE in one column; below-the-bar dots; threats per readout; flags with their follow-up fact; signals read per readout; the gene sets' arrays one cell; the rules, their changes after review and the chance calibration in its notes); then the two flags, one slide each (user 2026-11-07: "26～29页重构我们不需要这么复杂，就看看这两个信号，说明其被flag的来源"): adarb2 and linc, each with the same figure (every group in the calls and in the reads; the group the flag comes from shaded) and the same four steps (what is flagged; which cases make it; are their genotypes right; the source of the flag), every number on the slide from full mainland (round 11; the other cohorts in the notes), the verdict in the conclusion; the unflagged signals' key numbers are in the summary notes |
| Summary | summary (round 11: the four questions only, one row each: the question, how many signals it flags or contradicts, large, and a one- or two-sentence answer; the former boxes (flags, not flagged, numbers, caveats and open decisions) are word for word in the notes) |
| Appendix | backup only (three slides): apGates (all ten gates in two tables, grouped as the two fidelity checks: reproduces, rule, result), apDsGenes (gene sets at the controls' depth: P, replicates below Bonferroni, case carriers; reading box), apTbl (the matrix as one table: 21 signals with discovery level and P) |

Grammar (as `tmp/ppt/hla_assoc`). Palette: navy `1F3864` (headers, heads), grey `595959`, ink `212121`, rules `BFBFBF`,
band `DEEAF6` (label column, chips), pale `F4F7FB` (boxes), SIG `FCEFD2` (a flag, a key row); cohorts narrow `0571B0`,
intermediate `008837`, full `B2182B` in every three-cohort table (the N / I / F header cells, or the row header when cohorts
are rows) and figure; cases dark red, controls gray, reference panels black diamonds. No red fill marks a flag: a flag is
a SIG cell in bold (in the evidence matrix a bold navy border and "!", the flagged signal's row on a SIG band), and a flagged dossier's status strip is SIG with navy text (BAND when not flagged). Chips under the
title, right-aligned, name the threats a slide tests ("T1 Depth"); the "if an artifact" prediction is the title of the
results table; figures are lettered a, b, c with short centered titles; text boxes are pale with a navy border and a bold
lead-in ("Decision:", "Rationale:", "What it can detect:", "Verdict:"); every slide ends with an underlined conclusion;
footers give the abbreviations and, on the background slide only, the numbered references ("[1] Author (year) Journal
vol:pages"), cited as "[n]"; every other source (methods, panels, the clonal-hematopoiesis literature) is named in the
speaker notes, so the slides keep their room for the results. There is no reference list slide. The heterogeneity results show the five smallest Q P beside a scan box (counts N / I / F).

Templates. Every design card shows six fields on the slide (question, what changes, who changes N / I / F, decision
rule, experimental controls, blind spot) beside a to-scale schematic, and a box with what the experiment can detect and why it is built
that way; the rest of the card (held constant, comparator, model and statistic, readout, positive and negative controls
in full, processing chain, output tables) is in the notes and `../docs/DESIGN.md`. Every results slide has threat chips,
the figure, the prediction as the table title, a three-cohort table of at most six rows and a verdict; the observations
behind it are in the notes. Every dossier has a status strip, a figure, at most six rows of per-check evidence N / I / F and
a verdict. The follow-up tables (`13.followup`: group frequencies, gene carriers, pre-QC records, the window around a lead,
DQB1:57:V typings) are descriptive counts behind the flags, not a check. Native PowerPoint shapes and tables everywhere
except the diagrams a layout tool draws better (master workflow, numbers flow): `../docs/design/*.png`, written by
`../scripts/design_figures.py` with their sources and placed at their natural size, so their text stays at least 9 pt. The
same design, longer: `../docs/DESIGN.md` (`design_figures.py --check-doc` fails when it is out of date).

## Files

| File | Role |
|---|---|
| `collect_stats.py` | Reads the results (read only), `references.json`, and `numbers()` of `../scripts/design_figures.py` (registry, design numbers, run traces); writes `stats.json`, the only source of numbers for the slides. Checks the final results' consistency (below). |
| `references.json` | The literature, each record checked against Crossref and PubMed (2026-10-31); the background slide's footer cites five, the speaker notes name the rest. |
| `make_figures.py` | Figures redrawn at slide size from the result tables; the reason for each redraw is in its header. |
| `build_slides.js` | Builds `platform_depth_summary.pptx`, `slide_keys.json` and `slide_titles.json` (pptxgenjs 3.12.0, jszip 3.10.1); estimates every table and text height as `check_layout.py` does and stops when an element would run into the conclusion. |
| `deck_check.py` | Every number on a slide, footer or note must be found in `stats.json`; text rules (American spelling, no DiseaseLens, no composite numbers, no internal path on a slide except the verbatim rules slide; SAIGE only there and on the models slide); every title fits one line (width in Liberation Sans Bold 24 pt ≤ 11.9 in, from `slide_titles.json`). |
| `check_layout.py` | Rendered-text collisions (PDF) and shape overlaps (pptx). |
| `render.sh`, `fonts.conf` | PDF via LibreOffice with DengXian mapped to Liberation Sans; optional slide images. |

Figures. The component's figure families (`<results>/figures/0N.*`) are reused as is, except: the evidence matrix, the
platform and control-source forests (full mainland), the ADARB2 and LINC01146 frequency panels (registration-study
split), the gene-carrier panels, the DQB1:57:V typing trace (`v57_trace`, drawn at its slot under the dossier's verdict: the family figure prints its cells at 4.2 pt;
not drawn when no down-sampled case changes its 57 V dosage, and the slide then says so), the leave-one-out scatter
(`loo_scatter`: the family figure's labels are 4.5 pt and overlap), the array evidence per cohort
(`array_panels.full_mainland`: 5–5.5 pt text at slide size), and the three-cohort lead-frequency and gene-carrier panels of
the appendix. Drawn for this deck's grammar: the four question schematics (`fig_q1`–`fig_q4`, hypothetical values, at
the size of their cell), the sequencing groups to scale (`background_depth`), Q1's emulation as before -> after against the controls
(`emulation_dumbbell`: a, each case group's depth, all reads -> the controls' depth, after = fraction kept x depth with
the achieved depth checked by G8; b, the genotypes failing a QC rule), Q1's signals on the |z| scale (`ds_zscale`: one
panel per cohort; each lead's original |z|, the replicates at the controls' depth, and the |z| if depth had made the
deep cases' whole excess, from `ds.gwas.tsv` and `stats.json` `depth_reach`; the drop must equal the matrix reading;
HLA markers and gene sets on the same |z| scale (z from the two-sided P; the test that called each set), with no
threshold of any kind (no significance lines, no flag threshold, no reach) and no range lines: original and the median
of the 3 replicates only, as on slide 10 (medians, no IQR); decluttered: one shared |z| axis title, the gene-set
test named in the caption, larger markers; slides 10–11 decluttered again (user): no value labels on slide 10's figure and no G400RS 15x row (the caption
states it), a three-line checklist on 10 and a three-row table on 11, no verdict box on either (the conclusion line
carries it); user 2026-11-04), DQB1:57:V across the arms with the cohorts in their colors
(`v57_pg`, from `ds.hla_dqb1_57v.tsv`), and the per-genome depth
(`depth_platforms`, from `01.depth`: the family figure colors two platforms in the cohort colors). The typing trace
names a case by its platform (`stats.json` `v57_traced_platform`), never by its sample ID. The per-platform concordance
of Q1 is a table (the family figure's text is about 6 pt).

Numbers computed in `collect_stats.py` from the result tables (not in any table): the depth reach of Q1 (allelic test
with the deep cases set to the controls' frequency), the cohort overlap, exact tests on allele counts (ADARB2), Fisher
and Poisson tests on carrier counts (RPH3A driver, JAK2 by case group, the control deficit against ToMMo), the indel and
SNV λ of the DNBSeq-T7 stratum on the null markers (mean over the MAF bins). Read from outside the results: the number
of mainland ancestry components of each cohort (`PopGMM_output/keep_list_summary.tsv`, `cohort_k`); the thresholds of the
main analyses are parsed from their discovery readings (`thresholds`).

## Results, development tables and strict mode

`collect_stats.py`, `make_figures.py` and `design_figures.py` read the same environment:

| Variable | Meaning |
|---|---|
| `PD_RESULTS` | the results root (default `../results`) |
| `PD_TABLES` | optional: a flat folder of the report stage's tables (`ds.*`, `array.*`, `controls.*`, `followup.*`, `evidence.*`) read instead of their published folders, e.g. `../.dev_report/tabs` |
| `PD_FIGURES` | optional: the figure families (default `<PD_RESULTS>/figures`); `build_slides.js` takes the folder from `stats.json` |

Without `PD_TABLES` the run is strict: every table must hold exactly the arm plan's replicates (r01–r03), the
superseded R = 10 arms must be gone, and the evidence tables may hold only the registry's checks. With `PD_TABLES` the same
rules print warnings, numbers missing from the development tables print as "pending", and the deck says "built from
development tables". Wording that depends on the replicates is computed, not written: the DQB1:57:V title, verdict and
summary line (a GEC flip, single typings or one case's rare allele), the traced cases and what changed in each, the JAK2
depth sentence and the appendix captions.

The imputation tool of the arrays lives in one key, `stats.json` → `imputation` (tool, version, record, panel). It is
parsed from the `imputation` record of `11.array/array.summary.json` in `numbers()` of `../scripts/design_figures.py`
(`N['imputation'] = imputation(asum['imputation'])`, the one line to change to name another tool); every slide and
`docs/DESIGN.md` take it from there.

## Rebuild

```sh
PY=/home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3
NM=<dir with node_modules for pptxgenjs@3.12.0 and jszip@3.10.1>
(cd .. && PYTHONDONTWRITEBYTECODE=1 $PY -B scripts/design_figures.py)            # docs/design figures, docs/DESIGN.md
(cd .. && PYTHONDONTWRITEBYTECODE=1 $PY -B scripts/design_figures.py --check-doc) # expect: no difference
PYTHONDONTWRITEBYTECODE=1 $PY -B -W ignore collect_stats.py                      # strict: no PD_TABLES
PYTHONDONTWRITEBYTECODE=1 $PY -B -W ignore make_figures.py
NODE_PATH=$NM node build_slides.js --apply-theme <pptx skill>/scripts/apply_theme.js
$PY -B deck_check.py                                                             # expect: 0 not found, 0 text-rule issues
bash render.sh [img_dir]                                                         # platform_depth_summary.pdf (+ images)
$PY -B check_layout.py platform_depth_summary.pptx platform_depth_summary.pdf    # expect: 0 and 0
```

Final build 2026-11-04, strict (no `PD_TABLES`), from the results of full.4 and report.8, in the grammar of
`tmp/ppt/hla_assoc` (questions Q1–Q4; no roadmap or dividers): `design_figures.py --check-doc` up to date; deck_check
2,041 numbers on 54 slides, 0 unmatched, 0 text-rule issues; check_layout 0 collisions, 0 shape problems; pptx
`validate.py` passed; `verify.sh` all passed (round 4, text and logic pass: the gates open the results; summary questions
counted on the signals each can rule on; visual pass: Q1's readout drawn as the rule reads it, the P scale in the
appendix; every figure placed at 0.92 or more of its drawn size; then two independent reviewers verified the build
against the result tables and their must-fix items were applied: the 11 below-bar leave-one-out readings named on the
summary and the appendix tables, the 16 array units that could be contradicted spelled out, the summary title stated as
a fact).

Round 5, 2026-11-04 (appendix and citations): the appendix cut from 26 to 13 slides, ordered as the main slides they back.
Removed as duplicates of main slides or notes: the reference list (2 slides), the problem drawn natively (= background),
the numbers flow (= workflow; compute in the notes), the P-scale leads (= the Q1 signals slide and the per-lead table),
the replicate count (an "If asked" note on the Q1 design slide), the narrow and intermediate forests and array panels
(the Q P table and the verdict table give the three cohorts side by side), the leave-one-out table of the two flagged
leads (= the leaving-out slide and the dossiers). Moved into the main text: the coverage strips (after the Q4 design)
and the DQB1:57:V typing trace (into its dossier, under the verdict). The two carrier slides merged into one. Slide
citations kept only for the prior evidence of the background slide (Clayton 2005, Tom 2017, Garner 2011, Thuesen 2022,
Guo 2018); Kim 2021, Lee 2017, Cheng 2014, Kishikawa 2019, Lin & Sullivan 2009, Li 2012, Tadaka 2024, Xie 2014 and
Jaiswal 2014 are named in the notes. Checks: deck_check 1,667 numbers on 42 slides, 0 unmatched, 0 text-rule issues;
check_layout 0 / 0; validate passed; no sample ID in slides or notes; `verify.sh` all passed. The cover follows the group's decks (tmp/ppt/hla_assoc): "Chronic Thromboembolic Pulmonary Hypertension
(CTEPH): Technical Artifact Assessment of Association Signals in JHRP v6", a statement, not a question, with no
subtitle; JHRP v6 is not spelled out there. Threat model revised the same day (user): T3 is a DNBSeq-T7 calling error
against HiSeq X (Q3's question: "Do DNBSeq-T7 calls make the signal?"), the filters threat removed, T5/T6 renumbered T4/T5;
`scripts/design_figures.py` threats() and COVERAGE, DESIGN.md regenerated, stats.json refreshed (design keys only).
Text sweep (user): the overview figure (master_workflow) one idea per line and no unit counts (168 signal variants,
genome-arms, defined/applying readings moved to the notes), drawn to fill the slide; numbers a reader would ask "why"
about removed from the slides (Q3 discordant-site counts, Nagahama counts on the Q4 card, the controls' imputation
panel on the strips, the T7 × 0.965 caption detail, the ADARB2 window counts, "381 readings apply" on the summary).

## Open decisions (for the group) and next steps

The evidence-matrix rules (fixed before the final readings, to be confirmed); ADARB2 and LINC01146 (keep with a flag or
exclude); the DQB1:57:V wording (no flag; its GEC call sits near the threshold); the RPH3A wording (proposed: consistent
in direction; the effect size depends on the control mix); the LINC01146 wording (proposed: flagged when the DNBSeq-T7
cases are left out; no calling artifact; the G400RS 15x cases run the other way). Next steps: a PCR check of the
ADARB2 haplotype in a few DNBSeq-T7 read carriers; a matched-MAF benchmark of case imputation r²; the non-converged
case-only fits.

Q1 rule changed after review (user, 2026-11-04): read down-sampling flags a signal when the replicate-median |z| falls
by more than 1 in any one cohort (was two of three; the other checks keep two of three). Registry
`scripts/checks.py` (`MIN_COHORTS_BY_CHECK`, `min_cohorts()`), `evidence_matrix.py` and `verify_checks.py` read it per
check; report stage re-run as report.9 (8 report processes, 82 cached, guard clean); docs/PIPELINE.md and DESIGN.md
regenerated; flags unchanged (ADARB2, LINC01146: leave one platform out); every depth reading consistent. The slide's
figure (`ds_zscale`) shows no reach and no threshold of any kind; the reach and the DNBSeq-T7 share of each lead's excess are in
the notes only.

All checks flag in any one cohort (user, 2026-11-04; was two of three): `scripts/checks.py` MIN_COHORTS = 1, report
stage re-run as report.10 (8 report processes, 82 cached, guard clean), flags unchanged (ADARB2 3/3, LINC01146 2/3), no
single-cohort flag; PIPELINE.md, METHODS.md, DESIGN.md and catalogue text updated.

Q4 redesigned on the Q2 pattern (user, 2026-11-04: "我们现在的Q2很好，请你据此修改优化Q4"): order ctrlDesign → looSrc →
ctrlRes; the Q4 card fields moved to the notes; new stats `source_strata_example` (FTO 77–502 control ALT copies per
source, OR 1.40–1.72; INPP5D 2–27, OR 1.22–3.79; full mainland, original arm), `source_rare_small` (rare leads keep 1–27
control copies in the three smallest sources) and `gene_loo_src` (16 of 40 gene readings nt, F); looSrc throws if a
source flag or a below-bar reading other than RPH3A without Nagahama appears.
Q4 then given its own logic (user, 2026-11-04: "Q4也有其特殊性，优化，不要一味一样"): every case stays, so leaving out a
source costs little (z_exp / z_pooled 0.89–0.90 without Nagahama, 0.98–1.00 otherwise, 0.55–0.57 without the DNBSeq-T7
cases); almost all the chance of a flag (0.29 of 0.31) sits with Nagahama; a small source is read against the rest and in
Q; direction decides (ABI3BP/ACC, CDH4/HTLV-1, INPP5D/HTLV-1 weaker; ADARB2/ACC stronger, z unchanged without ACC). New
stats `src_loo_power`, `t7_loo_ratio`, `source_pairs`, `gene_loo_src` (weak pooled 15, lost by leaving out 1).
Q4 refocused (user, 2026-11-04: "我们的重点是没有标记因为ctrl分层而改变"): the source-against-rest table and the
direction argument left the slides (notes only); slide 19 counts every reading by signal type (GWAS 240/240 consistent,
DQB1:57:V 15/15, gene sets 75 consistent + 45 nt incl. 3 below the bar; flags 0), new stat `src_loo_kinds`; slide 20
shows Q P of all eight common markers.
Q3 restructured (user, 2026-11-04: "重构优化 Q3"): design slide in the Q2/Q4 form (card fields to the notes); result
slide title "Arrays Contradict No Signal in the Same People" (was "Support 20 of 25, Contradict 0 of 16"); verdicts
grouped by evidence (typed 8 reproduced; imputed 6 confirmed + CDH4, PHF20L1 not confirmed; drivers 6 sets confirmed +
JAK2 partly; ADARB2 not checkable; contradicted 0 of 16); fixed a stale "in at least two cohorts" in the prediction.
Q3 matching slide added (user, 2026-11-04: "每一个marker是如何匹配上的说明白，系统，直接"): arrayMatch between the design
and the results; new stat `array_match` (per signal: strongest match in F, WGS and array marker, proxy r² and case r²
over N–F, other matches, gene drivers typed / imputed); HS3ST3A1's array SNP is a proxy only in F (r² 0.82), a
credible-set member in N and I; MHC added to the glossary. 45 slides (32 main).
Where down-sampling acts (user, 2026-11-05): the Q1 design slide gets a "Where the reads are down-sampled" table
(23 regions, 930 kb around the 168 GWAS signal and 36 gene variants, ± 10 kb; the HLA reads: chr6:29.54–33.42 Mb,
16 ALT and 525 HLA decoy contigs, unmapped; rest of genome unchanged) and "in the reads around the signals and the HLA
reads" in its What-changes row; new appendix slide apRegions lists the 23 regions with their signals and variant
counts. New stats `ds_regions`, `hla_reads`. METHODS §3–§4 state the regions. 46 slides (32 main, 14 appendix).
Variant calls are described on the slides by the 16 GWAS loci and the 4 genes only (user, 2026-11-05): slide 9 row
"the 16 GWAS loci and the 4 genes: each signal or tested variant ± 10 kb"; apRegions lists one row per locus and gene
(region, kb, signal or tested variants) beside the HLA reads; the 23-interval / 930-kb / 2,156-site detail stays in
the notes and METHODS.
DQB1 position 57 left out of the deck (user, 2026-11-05: redundant with DQB1:57:V): collect_stats.py and
make_figures.py drop kind `hla_position` from the evidence tables, so the deck reads 25 signals (525 cells); every
slide text, the models table, the matrix legend ("reported"), apTblOther and the master-workflow box lose it; the
pipeline still tests and reports the position.
Q3's DQB1:57:V check uses HIBAG only (user, 2026-11-05): scripts/array_hibag.py + array_hibag.R in ARRAY_EVIDENCE
(report.12, 8 report processes, 82 cached, guard clean; report.11 failed on an input path, fixed), HIBAG 1.42.0 and its
published Asian model in assets/hibag; evidence label `HIBAG-imputed`; OR typed vs HIBAG N 0.646/0.643, I 0.634/0.632,
F 0.581/0.578, all reproduced. The ridge prediction stays in 04.compare for gate G10 only.
Q3 simplified (user, 2026-11-05: "Q3不要这么复杂"): one marker per signal — the lead typed, else imputed (Minimac4
R² ≥ 0.8, one R² per marker, full mainland), else an LD proxy (r² ≥ 0.8); DQB1:57:V by HIBAG (on the QC'd array SNPs);
gene drivers typed when the whole group is typed, else imputed; cases and controls from the same kind of data; the
association re-run in the same people. array_evidence.py, checks.py (ARRAY_READING), verify_checks.py, figures.py
and catalogue.py changed; report.13 (8 report processes, 82 cached, guard clean). Verdicts: 16 / 16 / 15 of 17
reproduced (PHF20L1 not reproduced in F only, array P 0.06), ADARB2 not checkable, genes RPH3A confirmed, STBD1 /
OSBP / JAK2 partly; none contradicted. Q3 is now four slides: arrayDesign (question + q3_flow schematic: QC'd arrays
1,750 / 409,088 → Minimac4 42,962,013 → one marker per signal → same signal?), hibagMethod (hibag_scheme + HLA-HD
vs HIBAG table; model facts in assets/hibag/DQB1_model.json from scripts/hibag_model_summary.R), arrayMatch (marker
and R² per lead; driver sources), arrayRes (array_panels: OR, imputation R², driver carriers; verdict counts).
47 slides (33 main, 14 appendix).
Q3 restructured again (user, 2026-11-05: "重构16～19，让问题说的更加清晰明白"; HIBAG drawn like the HLA-HD slide; the
HLA-HD slide brought in for the comparison): arrayDesign (question + q3_flow) → arrayMatch (marker per signal) →
hlahdMethod (hla.typing.v2's 'The HLA-HD Typing Algorithm', figure concept_hlahd.png copied unchanged by
make_figures.py, numbers from hla.typing.v2/deck/stats.json) → hibagMethod (hibag_algorithm: worked example, three
steps, illustrative values) → hlaCompare (HLA-HD vs HIBAG table; agreement N / I / F: DQB1 two-field ~95–96%, 57 V
100% cases, 99.6–99.7% controls; OR side by side) → arrayRes. Reference zheng (HIBAG, Pharmacogenomics J 2014) added
to references.json. hibag_scheme retired. 49 slides (35 main, 14 appendix).


Q3 genes: carriers counted, no verdict (user, 2026-11-05: "针对hla和gwas我们查看其信号是否一致。而针对罕见变异的driven
变体，我们看一下我们插补的芯片可以复现多少carrier（只看full就行）"). The arrays cannot reproduce a rare-variant gene
test (same people; the controls are in the imputation panel; imputation misses many rare carriers), so the gene sets get
no array verdict. scripts/array_evidence.py counts, in full mainland only, the WGS case carriers of each driver that
the imputed arrays recover (dosage ≥ 0.5), with the number expected from the genome-wide recovery of case carriers of
rare variants with as many control carriers (06.imputed/recovery_by_group.tsv); typed array calls are no longer used
for the genes; controls listed, not counted. Matrix: gene array cells 'reported' in F, 'not applicable' in N and I
(checks.py ARRAY_GENE_COHORT); 'partly supported' retired. Deck: arrayDesign, arrayMatch (driver table: WGS case
carriers, in the imputed arrays), arrayRes (verdicts for the 17 association signals; gene drivers: 16 of 22 case
carriers recovered, about 14.6 expected; JAK2 V617F not in the imputed arrays), array_panels b = imputed leads only,
c = case carriers with the expected count; matrix legend row for 'reported'; RPH3A, STBD1/OSBP, JAK2 dossiers,
summary, apArrayTab, apTblOther updated.

Q3 marker quality (user, 2026-11-05: "每个变体有一个minimac的R2，同时应该还是有case ctrl的concordance（基因型比对）"):
array_cohorts.py step 06 writes 06.imputed/concordance.tsv (best-guess imputed genotype vs WGS, cases and controls
apart, and among people non-reference on either side; not part of gate G10); array_evidence.py adds conc_cases /
conc_controls to every evidence row (typed: the array call from 04.compare discordance counts; DQB1:57:V: HIBAG's
best-guess 57 V agreement). Slide 17: both tables show R² (Minimac, one per marker) and the concordance in cases and
controls; slide 21 panel b shows R² and both concordances per imputed lead. Controls sit near 100% (panel members);
CDH4 is lowest in cases (about 88%).

Restructure (user, 2026-11-05: "1. 针对Q3重新规划设计，展示的更加直观专业学术 2. 针对25页之后，请你重新规划设计…删除不重要的附件内容
3. 针对个别信号的深入，我们只深入被flag的信号，且case不分新旧研究"):
- Q3 = 7 slides: arrayDesign (question box, compact q3_flow, readout table "If DNBSeq-T7 calls made the signal", control
  and blind-spot box) · arrayMatch · hlahdMethod · hibagMethod · hlaCompare · arrayRes (GWAS leads and DQB1:57:V: OR and
  R²/concordance panels, verdicts N / I / F, genome-wide agreement as the control; appendix apArrayTab folded in) ·
  arrayGenes (new figure array_genes: per-driver case carriers recovered vs expected; why count, not test).
- After Q4: the sample-slices slide removed; blind = "Coverage" (threat × experiment, gaps); ruleDesign = an academic
  table (check, question, statistic, flag rule, signals read, threat) with the changes after review as a footnote (the
  appendix rules slide removed); matrix; only the flagged signals keep a dossier (adarb2, linc); the DQB1:57:V,
  RPH3A, STBD1/OSBP and JAK2 dossiers removed, their key numbers in the summary notes; summary.
- Cases are never split by registration study (freq_ADARB2 / freq_LINC01146 redrawn without the new/old triangles;
  notes cleaned). LINC01146 dossier adds the T7 cases' genotype concordance on the arrays.
- Appendix: apRules, apArrayTab, apLeadFreq1/2, apCarriers removed (figures v57_*, carriers_*, lead_frequencies_cohorts.*,
  gene_carriers_cohorts retired and listed in make_figures.py's stale-file list).
Review round (independent reviewer agent, same day): array ORs per ALT (pipeline, report.16); CDH4 "reproduced at the
association level only" (case non-reference concordance 22%) on arrayRes and summary; slide-17 drivers show the WGS case
carriers found; hlaCompare: concordance wording and the HLA-HD − HIBAG non-differential P; rules retitled "Fixed
Rules", Q1 reach and the source row written out; coverage/rules say "source of bias"; matrix legend "reproduced
(arrays)"; looSrc no longer cites a removed slide; the removed dossiers' numbers in the summary box and notes; apQ removed
(repeats platRes/ctrlRes); master_workflow box without the 26-signal count. Not adopted (user decisions): r² instead of
concordance; typed drivers in the gene readout (notes only, If asked).
Slide 22 (user, 2026-11-05: "22页中我们只要a图就行吧？另外如何计算的expected，你需要说明白（可以用公式）"): array_genes is one
panel (the recovery-by-bin panel b removed); the right column states the formula E = Σ_d n_d × r(c_d) with n_d, c_d and
r(c) defined (5,685 rare untyped variants within 500 kb of the drivers, same bin of control carriers), and a table
computes it per driver (c_d with its bin, r(c_d), n_d, E_d, recovered; total 20 / 13.6 / 14). The bin values stay in
the notes. Earlier the same day the reference was mislabeled "genome-wide"; it is the regional ±500 kb set.
Then (user: "22页说的更加直观，专业，明白一些") the formula box was replaced by two worked steps without subscript notation:
1 · Recovery of comparable variants = a look-up table (control carriers 0 / 1 / 2–3 / 4–10 / ≥ 11; reference variants
585 / 874 / 1,249 / 1,657 / 1,320; case carriers recovered 41 / 54 / 73 / 78 / 83%; stats.json
array_drivers.reference_variants_by_bin added); 2 · Expected for each driver = WGS case carriers × that recovery, a
table "Driver · Control carriers · Carriers × recovery (e.g. 10 × 78%) · Expected · Recovered" with the total. The
controls' 99.6% moved into the "Why count, not test" box; JAK2 V617F as a footnote; conclusion states 14 of 20 against
13.6 expected and no false carrier; figure legend "expected (carriers × recovery)".
Then (user, same night: "好的，给我改成描述。证据矩阵中关于这部分也不用分N I F了，因为反正都是subset，三个格子我人诶可以合并，然后这里都是描述证据"):
the gene drivers are a description, no expected number. Slide 22 = "Imputation Confirms 14 of 20 Gene-Driver
Carriers": figure without the expected tick; confirmed per gene (STBD1 4 of 5, RPH3A 7 of 10, OSBP 3 of 5), none added;
what a confirmed carrier shows (independent of the reads), what a missed one shows (nothing against; context 41–83%
recovery nearby), controls in the panel, JAK2 V617F not counted. Evidence matrix (pipeline report.20 and deck): the gene
sets' array check is ONE cell for N / I / F (evidence.long.tsv: one full-mainland row; verify_checks updated), drawn
with its count (e.g. 4/5); grid cell counts subtract the merged cells. Wording "recovered" → "confirmed" on arrayDesign,
summary, apTblOther, notes; q3_flow box "WGS case carriers the imputed arrays confirm". array.evidence.tsv keeps the
expected count as an unread diagnostic (METHODS §7).
L/M/H removed from the deck (user, 2026-11-05: "我想从ppt中删除所有L/M/H的结果，只要M/H。数据你可以留着，但是判断以及整个deck要删除";
L/M/H is not the main analysis). Pipeline: checks.GENE_SETS_READ = ['moderate_high']; the evidence matrix reads the M/H
sets only (report.21; 21 signals); the L/M/H re-tests stay in the stage tables as data. Deck: no L/M/H reading, count,
threshold, table row or abbreviation; the dagger for tiered sets dropped (every M/H set read is tiered); build_slides.js
throws if a non-M/H gene set reaches it.
Slides 26 (blind, "Coverage") and 27 (ruleDesign, "Fixed Rules") removed (user, same night: "26，27我认为都应该删除"):
their gaps are in the summary's Caveats; the partial coverage and the gaps went to the threat slide's notes (If asked);
the rules, the changes after review, the leave-one-out conditional model with the winner's-curse citations and the
chance calibration went to the matrix slide's notes (If asked). Script before the cut: scratchpad
build_slides.before_cut_26_27.js (session only).
Appendix review (user, same night: "从30也开始（包括30页），仔细看看我们那些需要那些不需要，需要的保留并优化排版和展示方式，不需要的删除"):
8 → 3 slides. Kept and redone: apGates (all ten gates G1–G10, before only five were shown), apDsGenes (wider carriers
column, a reading box: every gene set keeps its Bonferroni call in every replicate — checked in the build — and STBD1
loses at most one case carrier), apTbl (apTblGwas + apTblOther merged, one 21-row table, flagged rows highlighted).
Removed: apDepth (depth and read length already on the background and Q1 slides; per-genome files named in the notes;
figure depth_platforms retired), apRegions (the Q1 design slide states the regions; coordinates in
results/02.targets/hc_intervals.tsv), apDsEmu (the Q1 emulation slide is its positive control; per-platform
concordance, heterozygotes, call rate and missingness were already in its notes, plus the all-cases and G400RS 15x
shares now), apDsLeads (the Q1 signals figure is the same data; replicates with P < 10⁻⁵ added to its notes). Notes
that pointed to removed slides were rewritten. Script before this step: scratchpad build_slides.before_appendix.js.
Q3 rule (user, same night: "reproduced = same direction, array P < 0.05.我们不强制要求，只要flag if opposite in direction or
|zWGS − zarray| > 3, in any cohort就行，另外可视化出来二者z的差别"): array_evidence.py verdicts are consistent /
contradicted / not checkable (array P kept as data; z_wgs and z_array added); the matrix reads 'consistent' (no
'supported' / 'not confirmed'); report.22. arrayRes (slide 21) shows the new figure array_z: a, full mainland, z on WGS
(open) and on the array data (filled), signed to the WGS direction, the opposite direction shaded; b, z_WGS − z_array in
N / I / F with the ±3 flag zone shaded; the table counts consistent / not checkable / contradicted and the largest
|Δz| per cohort. The R²/concordance panel left the slide (it is on the marker slide). "Reproduced" wording replaced by
"consistent" on the Q2 slides, the design card, the LINC01146 dossier, the summary and the appendix table.
Q3 cohorts (user, 2026-11-06: "Q3中的N I F应该和其他问题面对的不一样，注意分辨"): Q3 reads the people with arrays in
each cohort, not the cohorts of Q1, Q2 and Q4. They are written N~arr~ / I~arr~ / F~arr~ (1,398 / 1,460 / 1,683:
DNBSeq-T7 cases 279 / 289 / 299, HiSeq X controls 1,119 / 1,171 / 1,384), defined on the Q3 design slide and in the
footers (abbreviation "N~arr~ / I~arr~ / F~arr~"); cohort tables take `arr: true`; figures array_z and array_genes and
the matrix's SNP-array column label them N_arr, I_arr, F_arr; the summary's Q3 answer, the matrix note and the appendix
table say so too. Notes use N_arr / I_arr / F_arr.
Evidence matrix redesigned (user, 2026-11-06: "优化，专业化我们26页证据矩阵的可视化，更加直观，专业，系统"): one full-width
figure (make_figures.matrix): the seven checks grouped under Q1 Depth / Q2 Platform / Q3 Technology / Q4 Control source
with their cohort letters (the SNP arrays N_arr / I_arr / F_arr); the signals grouped under GWAS leads (16), HLA residue
(1), Gene sets M/H (4); a reading that does not apply is a dash, not an empty box; a right column names each signal's
flags with their cohorts; a bottom row counts the flag readings per check; the legend carries every reading's count.
The slide's legend table and verdict box are gone; one caption line states the cohort cells, that a flag never removes
a signal, and the readings below the testability bar.
JAK2 and the leave-one-out rule (user, 2026-11-06: "我更想给现在的决定一个合理性"; also "skat可以给beta吗" and
"我们没有使用去一的p值来检测这一点，是否合理"): the rule stays (CMC burden β for every set; JAK2 not testable). The
looRule slide carries a box "Gene sets: the CMC burden β" (z_exp needs an effect and its SE; SKAT gives none; JAK2,
SKAT-O with ρ = 0, is nt, not consistent; depth reads each set with its calling test) and two If-asked answers (why not
SKAT-O for JAK2: post hoc and no expected z; why not the P without the group: it is read where |z_exp| ≥ 2.8). looPlat
describes JAK2's SKAT-O P without each platform and its carriers by platform (stats jak2_skato, jak2_platform_carriers);
looSrc's notes describe it without each source; the summary's Q2 row gives the reason; METHODS §10 the full rationale.
hlaCompare (slide 20; user, 2026-11-06: "不用，你就报他们俩的P值就行"): the row "HLA-HD − HIBAG 57 V, cases vs controls: P"
(the typed-minus-imputed residual test) is replaced by "57 V P: HLA-HD (HIBAG)", beside the odds-ratio row; the residual
test stays in the notes.


Review round 7 (2026-11-06; user: two reviewers of numbers and logic, two of presentation and visuals, one integrating
editor, then automatic fixes). Reports, the integrated plan and the edit scripts: session scratchpad `review7/`
(A_numbers, B_logic, C_story_content, D_visual_layout, E_integrated). 131 findings → 63 items applied (11 P0, 16 P1,
36 P2); no rule, check aggregation, evidence table or settled decision changed; no pipeline re-run.
- P0 (facts): summary notes, RPH3A without Nagahama falls below 1.96 in N only (I and F testable and holding); looSrc
  gene box and Say name RPH3A in N beside OSBP; gates notes count the M/H G2 scans only; workflow notes say how the 433
  cells arise (four merged gene-array cells); dsDesign notes give each platform's median fraction; dsSignals notes: no
  57 V dosage changes, one call turns ambiguous in replicate 1 (guarded); dsEmulation Say "in full mainland";
  platDesign/ctrlDesign notes and DESIGN.md no longer say the bar "keeps a small group not testable" (small groups
  always pass it); linc table labels its array row N_arr / I_arr / F_arr; arrayRes source names array_z.png; METHODS
  §10 JAK2 CMC P in N / I / F order.
- P1 (logic and disclosure): no "1.96 plus 0.84"; the bar described as fixed in advance and conservative (new stats
  `loo_below_bar[].sd/.power_cond`; If-asked on looRule, matrix, looSrc notes); a consistent leave-one-out cell defined
  (matrix caption, apTbl legend); the reach of Q1 (6 / 7 / 7 leads, ≤ 23 / 21 / 22 %) and of leaving out DNBSeq-T7
  (testable for 5 / 6 / 5 of 21 signals; new stat `t7_loo_testable`) in the summary answers; random effects described
  as unable to fire for a one-platform effect (table cell, If-asked); platform-against-rest P_diff described for the
  rare leads (new stats `platform_pdiff`, `platform_differs`: INPP5D without G400RS 30x, MFSD11 without G400RS 15x, with
  allele frequencies; LINC01146's P_diff "does not pass a correction"); Q3 re-checks only the cases' calls at imputed
  markers (arrayDesign box, summary caveats) and the array products split by case status (If-asked); T3 tested in part
  by Q2 (threat slide, design_figures.py threats()); the source-check calibration stated as approximate (ctrlDesign,
  ctrlRes, looSrc notes, METHODS §8/§10); JAK2 V617F "not in the QC'd or imputed arrays (a multi-allelic array site,
  removed at array QC)" in sentences (short labels unchanged); arrayGenes: a confirmed carrier is evidence independent
  of the case's own reads, a missed one weak evidence against; summary: JAK2 read by depth only, STBD1/OSBP ToMMo
  shortfall in every cohort (notes give N / I / F and the hidden carriers), all seven below-bar readings named,
  "For the group" in bold, "Next" moved to the notes; slide citations on hlahdMethod/hibagMethod moved to the notes;
  "next slide" cross-references replaced by titles.
- P2 (visuals): figure text ≥ 9 pt at its placed size (background_depth, ds_zscale, q3_flow, hibag_algorithm,
  concept_hlahd slot, evidence_matrix, array_z); loo_regions with a filled observed point and a visible tail;
  loo_scatter: JAK2 label above its box, the below-bar band darker (#DDE1E6), below-bar labels placed by a small
  search for an arrangement without crossings (`_label_below_bar`; DQB1:57:V named by its square marker); forest
  clipped intervals end in arrows, an OR beyond the axis drawn as ">"; array_z one marker meaning per panel (HIBAG
  square, panel b circles); array_genes WGS case carriers in pale case red; matrix cells hold the count only, the
  JAK2 array cell gray (not checkable), the flags column widened; summary answers regular weight with bold leads;
  apTbl highlights every flag cell, "ns" for not significant.
For the user (not applied): U1 keep the 2.8 bar (conservative; a conditional bar would flag HS3ST3A1 I/F and RPH3A N);
U2 mark the hidden below-bar cells; U3 random-effects rule; U4 exact allele test by platform for rare leads; U5
empirical P for the source checks; U6 the short labels "not on the arrays"; U7 reach of the Q3 rule at weak leads;
U8 the Q threshold 0.05 / 17; U9 the cover date (2026-11-04).

Round 7 decisions U1–U9 (2026-11-06; user: "choose the reasonable ones and update"). No rule, reading, flag or
pipeline table changed; `verify.sh` ALL PASSED (frozen pieces unchanged).
- Applied:
  - U2: apTbl marks with "°" each leave-one-out cell that holds a below-bar reading: HS3ST3A1 platforms N/I/F,
    LINC01146 platforms F, DQB1:57:V platforms I/F, RPH3A M/H sources N. The legend defines "°"; a guard ties the marks
    to `loo_below_bar`. The main matrix is unchanged.
  - U3: the platDesign If-asked now says LINC01146 is the only common marker whose random-effects interval includes 1,
    and that a rule "interval includes 1" would flag LINC01146 alone, which leaving out the DNBSeq-T7 cases already
    flags. The pre-specified rule (a reversal) is kept. Guarded.
  - U4: new stats `rare_leads`, `rare_platform_exact` (Fisher's exact test of one case platform's ALT alleles against
    the other cases', every rare lead × platform × cohort, from followup.gwas_af.tsv), `rare_platform_exact_lt05` and
    `rare_platform_exact_pairs`. In the platRes notes:
    - INPP5D G400RS 30x has P 0.020 / 0.018 / 0.019, against P_diff 3 × 10⁻⁴ / 5 × 10⁻⁶ / 2 × 10⁻⁶;
    - MFSD11 G400RS 15x has P 0.076 / 0.076 / 0.073;
    - 4 of the 36 pairs have P < 0.05, all named (ADARB2 DNBSeq-T7 is its flag).
    platDesign card and Say: "an exact test of its allele counts". METHODS §10 updated.
  - U5: new stats `sources_null_q`, `sources_q_emp`, `loo_null_sd` and `loo_null_undefined`. They are computed in
    collect_stats.py from the subsets' SAIGE rows with the pipeline's formulas (platform_stats.rho_shared, maf_bin, GLS Q),
    validated against controls.heterogeneity.tsv (relative difference < 1e-8); no bytecode is written.
    - Source Q at MAF ≥ 5%: analytic P < 0.05 in 8.1 / 9.1 / 10.9% of the null SNPs, and < 0.0029 in
      0.85 / 1.04 / 1.27%.
    - ABI3BP's source Q P, read against its null MAF bin, is 0.042 / 0.047 / 0.107 (analytic 0.022 / 0.020 / 0.049).
    - Leave-one-out model on the null SNPs: SD of (z − z_exp) / SD is 0.93–1.01 for the platforms and 1.10–1.92 for
      the sources. The SE without a source is not above the pooled SE in up to 52% of the null SNPs (Lung cancer, F).
    Where the numbers appear:
    - ctrlRes caption, Say and notes;
    - ctrlDesign negative-control card and Say;
    - looSrc notes;
    - platDesign notes;
    - the matrix answer "Why 2.8?";
    - METHODS §8/§10.
  - U6: the short labels now say "removed at QC": the arrayDesign table ("*JAK2* V617F removed at QC"), the arrayMatch
    driver table, and apTbl ("V617F removed at QC"). Guards confirm V617F is the only uncounted driver.
  - U8: a platRes If-asked on 0.05 / 17 against the 8 markers Q reads. The smallest Q P is 0.013 for the platforms and
    0.020 for the sources, both above 0.05 / 8. Guarded; the threshold is unchanged.
  - U9: cover date 2026-11-06 (`const.deck_date`).
- Kept as they are:
  - U1: the 2.8 bar. It was fixed before the readings; its conservativeness and the readings a conditional bar would
    flag are already disclosed.
  - U7: no new concordance rule. Only CDH4 has a case non-reference concordance below 50% (22%; PHF20L1 is next at
    75%), and arrayRes already names both.

Visual fixes, 2026-11-06. User: "26，元素堆叠，图片左倾" (slide 26: elements stacked, the figure leaning left) and
"24，0那里有悬停的点" (slide 24: an unexplained point at 0).
- evidence_matrix (slide 26)
  - Cause: the legend's one row ran past the table. The figure is saved cropped around its content, so the table sat in
    the left three quarters and its heads collided (Q3 · Technology against Q4 · Control source; N_arr I_arr F_arr
    touching).
  - The axes now fill the figure, and the legend wraps within the table's width (two rows).
  - Spacing: questions stand 2.3 units apart, cells are 0.85 units wide, and the arrays' cells are 0.45 apart with their
    letters at 9.5 pt.
  - The build stops if any text leaves the frame.
- loo_scatter (slide 24)
  - Cause: the points at 0 are OSBP M/H without Nagahama in N and I. With no control carrier left, z_exp is 0.03.
  - Any reading with |z_exp| < `const.void_z` (0.5) is now named under its point. The looSrc caption says why
    (guarded against the data).
  - Labels below the bar now try the later rows before a place beside the point, so RPH3A M/H keeps its leader instead
    of sitting on the 1.96 line.
  - Slide 14's arrangement is unchanged.
- evidence_matrix legend (user 2026-11-06, option chosen from three): "not testable (6) or not checkable (4)" instead of
  "not testable or not checkable (10)". One gray swatch; the 4 are the arrays' cells without a marker (ADARB2 × 3,
  JAK2 V617F × 1). The matrix Say names them, guarded against the readings.

Evidence matrix rebuilt as the evidence chain (user 2026-11-06: "现在的证据矩阵太让人误解了。重新审视一下我们的整个证据链条，
设计一个更加系统，直观，专业的，减少误解"). Rules, readings, flags and the pipeline are unchanged. What misled, and the fix:
- Leaving out a case platform was one cell over four groups: "consistent" when no testable group flags. Without the
  DNBSeq-T7 cases (most of the cases) the rule reads only 5 / 6 / 5 of 21 signals, so most of those blue cells came from
  the three small platforms alone, as if the T7 cases had been tested; 7 below-the-bar readings sat inside blue cells
  (LINC01146 F looked consistent).
  Now two readouts: "Without one small platform" (T2) and "Without DNBSeq-T7 cases" (T3), the latter next to the
  arrays, which re-genotype the same T7 cases. Their union reproduces the registry's lopo reading (collect_stats
  checks every cell).
- Random effects cannot flag a one-platform effect (U3) yet stood as a column of its own. Now "Platform estimates (Q,
  RE)", flagged if either flags; the rules are unchanged.
- Below-the-bar readings are dots in their cells, a white dot on a blue cell (RPH3A M/H without Nagahama, N).
- The cohorts are nested, so three cells are one result seen three times; the caption says so.
- Each readout names the threats it tests (T1–T5) on a "Threats tested" row.
- The flags column carries each flag's readout and cohorts plus its follow-up fact (guarded):
  - ADARB2: "no T7 case carries it; no array marker";
  - LINC01146: "T7 calls agree with the arrays (100%)".
- A bottom row gives the signals each readout could read: 21; 8; 20; 5 / 6 / 5; 16; 8; 20.
- The gene sets' array cells are outlined (a description), not filled like a reading.
- The legend counts every cell: consistent 290; not testable 53 or not checkable 4; flag 5; below-the-bar dots 7;
  does not apply 78; gene carriers 3. Total 433, the registry's grid.
- New stats key `matrix2` (collect_stats.py): columns with question and threats; per signal × readout × cohort the
  reading and the below-the-bar mark; cell counts; signals read; flags per readout. Its guards check:
  - the split and merged readouts reproduce the registry's readings;
  - 433 cells, 5 flags and 7 dots;
  - the T7 reach equals `t7_loo_testable`.
- make_figures.matrix redrawn from `matrix2` (12.7 × 5.4 in; no text outside the frame).
- The matrix slide is retitled "The Evidence Matrix: 21 Signals, Four Questions, Three Cohorts, Two Flags", with a new
  caption (nested cohorts; why two columns), the abbreviations RE and nt, and new notes (why two columns, the cell
  vocabulary, the flags' facts, signals read, the readout ↔ registry mapping).
- apTbl uses the same seven readouts ("Platform Q and RE", "Small platform out", "DNBSeq-T7 out", ..., "Source out"),
  with nt° for a below-the-bar reading that is itself nt.
- The workflow slide's notes say the matrix shows the seven checks as seven readouts.
- Backup of the earlier version in the session scratchpad (u19/backup_matrix_v1).
- Legend completed (user 2026-11-06: the white dot on a blue cell had no legend). The below-the-bar dots now have two
  entries:
  - "below the bar (6): nt, and |z| < 1.96", a dark dot on gray;
  - "consistent, one left-out group below the bar (1)", a white dot on blue: RPH3A M/H without Nagahama, N.
  The matrix Say and Numbers name both, guarded against the readings.

CRAM check of the two flags (user 2026-11-06: "在27和28，我需要你查看cram原文件来判明问题（各自新加页面）"). Two new
slides: adarb2Cram (28) and lincCram (30); the deck now has 34 slides (31 main). Every full-mainland genome was read
(3,101 CRAMs, read only); the scripts, the method and the per-genome table are in cram_check/ (cram_check/README.md).
- ADARB2. The lead lies in a minisatellite of about 390 bp (31-bp units). The carrier haplotype is one unit shorter;
  the caller writes it as a 3-bp and a 28-bp deletion, and the aligner places its reads without a gap.
  - Method: carriers are read by a 16-bp signature, GATCGCATGATGTGTA. It is absent from the reference window, and
    185 of 186 called carriers among the controls (5 of 5 among the other cases) hold it.
  - Result: in the reads, 32 of 319 T7 cases carry the haplotype (10.0%), against 9.4% of controls and 9.3% of all
    cases. None of the 32 is called; the HiSeq X calls also miss 65 control carriers.
  - Mechanism: the T7 carrier reads map well (91% MAPQ >= 20), so the miss lies in calling. T7 is the only 100-bp
    platform; this is likely, not shown.
  - Verdict: a DNBSeq-T7 calling artifact.
- LINC01146.
  - Result: the T7 reads support their calls as well as the controls' reads do. Read genotype = call in 99.4% vs
    99.4%; heterozygote A fraction 0.50 vs 0.50; forward strand 0.51 vs 0.51; A in G/G 0.15% vs 0.10%.
  - The A excess is in the reads (0.70 vs 0.60).
  - Verdict: not a calling artifact; cause open.
- Statements the result settles were updated: "pileups pending" is gone from platDesign (card, If asked), the adarb2
  dossier, the summary and docs/DESIGN.md (scripts/design_figures.py).
  - adarb2 retitled "*ADARB2*: The DNBSeq-T7 Calls Hold No Copy of a Low-Frequency 3-bp Deletion", with a new verdict,
    mechanism and asks.
  - Summary title: "... T7 Calls Missed *ADARB2*'s Carriers", with new flags text, conclusion and next steps.
  - Matrix flags column (ADARB2 "T7 reads carry it; T7 calls miss it (CRAM)"; LINC01146 "... and reads") and the
    matrix Say.
  - Models and threat slide asks.
  - New stats key `cram`. New figures cram_ADARB2 and cram_LINC01146. New abbreviations CRAM and MAPQ.
- LINC01146 re-read (user 2026-11-06: "LINC01146 应该是G400RS 15X的问题，你看看"); lincCram rebuilt.
  - The platform strata (09.platform_assoc only_*, shared-control GLS as scripts/platform_stats.py): ORs T7 1.61,
    NovaSeq 1.32, G400RS 30x 1.12, G400RS 15x 0.73 (F). G400RS 15x is the only stratum departing from all cases
    (−2.5 to −2.7 SE, opposite direction).
  - Without T7, the GLS of the strata gives OR 1.07 (reproduces the pipeline's z 0.46–0.68). Without both T7 and
    G400RS 15x, the other 89 cases give OR 1.23, z 1.30 against z_exp 2.26 (chance 0.07–0.13): they fit the pooled
    effect.
  - Q falls from P 0.010–0.022 to 0.16–0.28 without G400RS 15x. Without T7 it falls likewise, so the tie is broken by
    the evidence: T7 genotypes are confirmed on arrays and in the reads and sit 1.3 SE from all cases.
  - G400RS 15x at the lead: 4 duplicates in 434 reads, homozygotes on 5–20 independent reads, read genotype = call
    100% (cram_check/linc_fragments.tsv). The genotypes are genuine, so this is a small-group deviation, not a calling
    artifact.
  - Genome-wide on 6,000 common null SNVs (cram_check/het_by_sample.tsv, plink2 --het --missing on the QC'd cohort
    genotypes), G400RS 15x genotypes are weaker: F 0.022 vs 0.005 in the controls (P 1e-4), 7.5% missing vs 2.6%, at
    14.2x. Controls at 14–16x show the same direction (F 0.012). T7 and G400RS 30x match the controls.
  - New script cram_check/linc_g400rs15x.py. New stats keys `linc2`, `g15`, `het_genome`. Figure cram_LINC01146 is
    now two panels (strata and combinations; reads per genome).
  - lincCram retitled "*LINC01146*: The Flag Comes From the G400RS 15x Cases, Not From DNBSeq-T7".
  - Updated to match: the linc dossier (verdict, conclusion, Say); the summary (flags text; caveat "weaker 14x
    genotypes"; Say); the matrix flags column ("from G400RS 15x, OR 0.73; T7 calls confirmed") and Say; the models ask.

Review round 8 (2026-11-06; user: two reviewers of numbers and logic, two of presentation, content and visuals with
condensation, one integrating editor, then automatic fixes). Reports, the integrated plan and the edit scripts: session
scratchpad `review8/` (A_numbers, B_logic, C_story_content, D_visual_layout, E_integrated; `apply/`). 74 findings → 37
items applied (4 P0, 11 P1, 22 P2); six questions for the user (U1–U6). No rule, reading, flag or pipeline table
changed; no pipeline re-run.
- P0: adarb2Cram notes, 11 (not 12) controls called non-carriers hold 3 signature reads (keys `cram.adarb2.*.nc_*`);
  adarb2, no read is assigned to the deletion in the T7 records (was "no ALT read", which the reads contradict);
  lincCram, G400RS 15x at 13.8x effective depth (14.2x is the sample sheet's measure, now labeled); lincCram's verdict
  no longer says the other cases fit the pooled effect without G400RS 15x (they exceed it).
- P1: lincCram rebuilt: each end against the other cases (the pipeline's P_diff, signed: G400RS 15x −2.7 to −2.8, T7
  +2.6 to +2.8), not against all cases (73% T7); Q P without each end; the strata stated as symmetric; the attribution
  to G400RS 15x rests on evidence outside them (this replaces the earlier tie-break "T7 … sits 1.3 SE from all cases",
  which was circular); the genotypes the association reads (QC'd; a genotype set missing entered as one A copy) beside
  the reads; new keys `linc2.*.strata.*.z_other/p_other`, `input_af`, `qc_missing`, `af_diff`, `g15.qc`. linc: empirical
  Q P labeled, stale clauses, conclusion de-duplicated. Summary: Q2 and Q3 say what each flag turned out to be;
  LINC01146 proposal and PCR next step; symmetry in the notes. ADARB2: mechanism spoken as likely; the 150-bp platforms
  call 190 of their 259 read carriers, T7 0 of 32 (key `cram.adarb2.*.p150`); the signature spans the junction of the
  two called deletions; nearby SNVs called less often in T7 (keys `adarb2.*.snv_t7_pct`, `snv_ctrl_pct`). platDesign
  If asked: LINC01146's RE interval reflects disagreeing platforms; the positive control's reach. "Signals read" on the
  design slides → "Applies to".
- Condensation: the gates slide merged into the workflow (33 slides, 30 main; its table, gates and questions in the
  workflow notes, every gate on apGates); per-slide timing (`MIN` in `build_slides.js`: 36.75 min, "about 37"; was
  40.5); Say from 5,708 to 4,459 words, detail moved to Numbers and If asked, none dropped; a guard stops the build when
  a main slide's Say exceeds 135 words per minute; section transitions; cross-references in the notes written as
  `'{{title:key}}'` and resolved to the slide's title at build time.
- P2: matrix flag lines spaced (row 1.2), no dash for "no flag", follow-up "(OR < 1)" (the flag fires in N and I),
  legend "nt, |z_exp| ≥ 1.96, |z| < 1.96"; titles of workflow, arrayGenes, matrix and the five signal and summary
  slides shortened, with a width guard in `deck_check.py` (`slide_titles.json`); loo_regions label off the diagonal;
  loo_scatter legend without "flag" when none is drawn; forests print Q P and I² for the common markers only;
  thousands separators; exact P at two digits; signed z in the looPlat and adarb2 tables; dsSignals caption names the
  test per gene set; footers; workflow figure names the matrix readouts (`scripts/design_figures.py`, DESIGN.md
  regenerated); cram_LINC01146 x labels on two lines; cram_check docs (about 390 bp; signature provenance).
- Checks: deck_check 1,793 numbers, 0 not found, 0 text-rule or title-width issues; check_layout 0 collisions, 0 shape
  problems; validate passed; verify.sh ALL PASSED; every Timing line "of about 37 min".
- For the user: U1 the LINC01146 symmetry (the strata alone cannot tell which end deviates); U2 delete hlahdMethod; U3
  merge samples and models; U4 the cohort of looRule's example (narrow; full mainland sits on the bar); U5 LINC01146
  with the QC-missing genotypes left out (one descriptive reading); U6 the source of the "institution" follow-up on the
  linc slide. Decided the same day (below).

Round 8 decisions U1–U6 (2026-11-06; user: "apply U1–U6 by the most reasonable recommendation"). No rule, reading,
flag, slide face or pipeline table changed (all 33 renders pixel-identical to round 8); four note items changed.
- U1, LINC01146 symmetry: kept as rebuilt. The title attributes the flag to the G400RS 15x cases; the notes state that
  the strata alone are symmetric and that the attribution rests on evidence outside them. Every reading agrees there
  is no calling artifact and the flag stays, so nothing decided depends on the tie.
- U2, delete hlahdMethod: kept (0.75 min); the HLA comparison reads better after both algorithms.
- U3, merge samples and models: kept apart; a merged slide would carry two messages.
- U4, looRule's example: narrow mainland kept as an exception to full-mainland examples. In full mainland z_exp is
  2.82, on the 2.8 bar, and the figure could not show the flag region apart (the reason is in the notes).
- U5, LINC01146 without the QC-missing genotypes: one descriptive reading in the lincCram notes (If asked; Numbers).
  Left out, the missing genotypes change the case–control difference in A frequency by less than 0.001 (0.098 / 0.098
  / 0.085 against 0.099 / 0.098 / 0.086 in the input; N / I / F), so their best-guess filling does not move the
  association. Which genotypes QC sets missing matters more: the 8 cases set missing carry 4 A alleles of 16 in their
  pre-QC calls (5 carry none; the cases that pass 0.68), the controls set missing carry A near their group's
  frequency (134 of 220, F); the difference is 0.077 in the pre-QC calls, 0.085 after QC and 0.080 in the reads (F).
  New keys `linc2.*.called_af`, `af_diff.called`, `missing_pre` (collect_stats checks that the counts add up and are
  whole); the build throws if the words stop holding.
- U6, the "institution" answer on linc: sourced, not dropped. The sample sheet (Providing_Institution) gives one
  institution for all but 5 cases in every cohort (434 of 439, F), every DNBSeq case among them; the other 5 are
  NovaSeq cases from 3 institutions, so institution cannot separate the platform groups. Registration study
  (followup.gwas_af.tsv) does not either: within each study the T7 cases carry more A than the G400RS 15x cases (F:
  0.69 against 0.56 in the earlier study, 0.74 against 0.50 in the later). New key `institution` (counts only); the
  build throws if a cohort departs. This reproduces the ad hoc follow-up of 2026-10-31, which had no source in the
  component.
- Checks: deck_check 1,808 numbers, 0 not found, 0 text-rule or title-width issues; check_layout 0 collisions, 0 shape
  problems; validate passed; verify.sh ALL PASSED.

Restructure of the flag slides (2026-11-06; user: "26～29 重构 … 只要能直观，系统地说明我们的问题就行"). Slides 26–29 now
explain the two flags as one system instead of four dossiers; the 29 other slides render pixel-identical.
- flags (26, 1.25 min): both flags through the same four checks (which cases move the signal; do the calls match the
  reads; do the arrays agree; what makes the flag), one table, the reading per flag. Replaces nothing; it is the spine.
- adarb2 (27, 1.5 min): the former adarb2 and adarb2Cram on the question "do the T7 genomes lack the deletion or do the
  T7 calls miss it". Figure cram_ADARB2 redrawn: a, how a carrier is read (the minisatellite to scale, the carrier one
  unit shorter, reads without a gap, the 16-bp junction signature); b, carriers called and in the reads; c, signature
  reads per genome. Table: without T7 z (z_exp), T7 alleles called, carriers called → reads by group, signature in the
  called carriers.
- adarb2Where (28, 1.25 min): new, from the user's question "我们有没有办法确认是calling的问题？". A pipeline diagram
  (native shapes): reads (present) → alignment (one command for every genome; no gap) → variant calling (lost: no T7
  record holds it; Q1's re-call of the same reads at full depth calls none of the 32 again and reproduces every other
  genotype) → QC (not involved) → association (a false protective signal); then why (likely, not shown: the only
  100-bp platform) and beyond (indel-specific λ in the T7 stratum; two experiments that would show the cause, not done).
  New script cram_check/recall_check.py (read only): aligner.tsv (every full-mainland CRAM's @PG, one command: bwa
  0.7.17-r1188 "bwa mem -p hs38DH.fa -", samblaster 0.1.24) and recall_adarb2.tsv (the baseline arm's genotype at the
  lead, 428 genomes: 424 called in both, all identical to the original calls; the 32 T7 read carriers 31 hom-ref and 1
  missing). New stats keys cram.aligner, cram.recall (collect_stats checks one command, identical re-calls, no T7 call).
- linc (29, 1.5 min): the former linc and lincCram. Figure cram_LINC01146 redrawn: a, each group's A frequency in the
  pre-QC calls (95% CI) and in the reads; b, each case platform's odds ratio and the strata combined. Table: without T7
  z (z_exp), reads against calls, arrays, each end against the other cases, without T7 and G400RS 15x. Platform Q,
  random effects and source Q moved from the table to the notes (all consistent readings).
- Keys adarb2Cram and lincCram retired; cross-references in other notes now quote adarb2 and linc. freq_ADARB2 and
  freq_LINC01146 are no longer drawn (stale list). PCR added to the glossary. Every sentence of the former four slides
  was compared with the new text; three details first missed were restored to the notes (the analytic Q P by cohort,
  the pre-QC T7 carriers by cohort, the shared-control covariance of the combined strata). Timing unchanged (5.5 min).
- Checks: deck_check 1,801 numbers, 0 not found, 0 text-rule or title-width issues; check_layout 0 collisions, 0 shape
  problems; validate passed; verify.sh ALL PASSED; no sample ID in slides or notes. Backup of the round-8 state:
  the session scratchpad, review8/backup_restructure_2629.

Flag slides simplified (2026-11-07; user: "26～29页重构我们不需要这么复杂，就看看这两个信号，说明其被flag的来源").
The four flag slides of 2026-11-06 (flags, adarb2, adarb2Where, linc) become two, one per signal, built the same way;
the deck now has 28 main and 3 appendix slides; slides 1–25 and the appendix are unchanged.
- adarb2 (26, 2.0 min): title "The Flag Comes From T7 Calls Missing a Deletion in the Reads". Figure cram_ADARB2
  redrawn: a, carriers per 100 genomes, called and in the reads (since round 9: calls open, reads filled), 95% Wilson CI, for the controls, all
  cases, the DNBSeq-T7 cases (shaded) and the cases without T7 (the reading that raised the flag); b, HaplotypeCaller's
  assembly in the diagnosed genomes, a square per genome (called 150-bp carriers 6 of 6, missed 150-bp carriers 0 of 2,
  T7 carriers 0 of 5). Steps: the flag (z without T7, N / I / F); the T7 alleles called (0 of 590 / 610 / 630); the T7
  reads (32 carriers in every cohort, 10.7% / 10.4% / 10.0%, controls 9.7% / 9.6% / 9.4%); the source (variant calling:
  HaplotypeCaller assembles the deletion from 150-bp reads, not from the T7's 100-bp reads; the re-call of the 319 T7
  genomes calls none of the 32 again).
- linc (27, 1.75 min): title "The Flag Comes From 31 G400RS 15x Cases; Genotypes Hold". Figure cram_LINC01146 redrawn
  in the same design: the A frequency of the controls, all cases, DNBSeq-T7, the cases without T7 (bold) and their three
  platforms (G400RS 15x shaded). Steps: the flag; the G400RS 15x cases (A 0.52 against the controls' 0.60; OR against the
  controls 0.74 / 0.76 / 0.73, corrected in round 9; without T7 the 120 remaining cases come near the controls); the genotypes (reads and
  arrays); the source (a small group, most likely chance; the strata alone symmetric).
- The overview (flags) and the pipeline slide (adarb2Where) are folded into the two slides' notes: the rule and the
  flag cells, "Why these questions?" and "Is a flagged signal removed?", the alignment, re-call, pre-QC, window, MAPQ,
  soft-clip and read-length numbers, the indel-specific λ of the T7 stratum, and the two experiments not done. The former
  figure panels (the minisatellite drawing, the signature reads per genome, the odds-ratio forest) are retired; their
  numbers stay in the notes. Every sentence of the four former slides was compared with the new text.
- New: HaplotypeCaller diagnostics (cram_check/hc_debug; hc_tables.py writes hc_genotypes.tsv, cohort_records.tsv,
  cohort_ga.tsv), stats key cram.hc (collect_stats checks: k = 85 in every assembly over the lead; the deletion
  haplotype in every called 150-bp carrier and in no T7 carrier; a heterozygous G>A in every T7 carrier; the G>A record
  removed by the site QC that the deletion passes; counts equal to cram_scan). The former statement "QC: not involved"
  is made precise in the notes: QC never touched the deletion's genotypes; it removed the G>A the caller wrote instead.
- Checks: deck_check 1,761 numbers, 0 not found, 0 text-rule or title-width issues; check_layout 0 collisions, 0 shape
  problems; validate passed. Backups: the session scratchpad, *.before_r10.*.

Round 9 review (2026-11-07; user: "调用2个agent审核ppt的展示逻辑和数字 / 再调用2个agent优化ppt的文本和可视化方案 / 最后用一个
agent总结，并自动开始优化修改"). Four read-only reviewers (A logic, B numbers, C text, D visual) and one integrator; plan and
log in the session scratchpad, review9/reports/{A_logic,B_numbers,C_text,D_visual,E_integrated,E_applied}.md. No rule,
reading, flag or pipeline table changed; stats.json gained keys only (loo_sd_genes; adarb2.<cohort>.qc_missing;
cram.hc.carriers, cram.hc.cohort.{cohort_genomes, ga_carriers_fm}, cram.hc.groups.*.sig_range).
- Errors fixed on the flag slides: the G400RS 15x odds ratio 0.74 / 0.76 / 0.73 is against the controls, not "against
  the other cases" (against them it is 0.47–0.48); "the other 120 cases" is now "the 120 remaining cases" (the 120
  include the 31); the notes' "the deletion's own genotypes are not touched by QC" is replaced by the counts (the record
  passes the site QC; the genotype filters set 93 control genotypes missing, 4 of them called carriers, and 6 case
  genotypes, no carrier; F); "re-run on 14 genomes" against 13 squares is now "13 carriers (and one T7 non-carrier)".
- Hedge restored on adarb2: step 4 compares like with like (every called carrier reconstructed, 6 of 6; no T7 carrier,
  0 of 5; calls 0 of 32 T7 read carriers against 190 of 259 on the 150-bp platforms) and ends "Likely cause, not
  tested: the T7's 100-bp reads"; the Say says so too. linc step 4 states the symmetry first and the evidence outside
  the strata as reasons, not a test (P_diff 0.006–0.008 does not survive the 64 comparisons, in the Say).
- One system: both slides ask the same four questions (What raised the flag? / Which cases make the flag? / Are their
  genotypes right? / Where does the flag come from?) with a bold answer first; same grid (FLAG_FW 5.9; strip 3.0 / 3.04
  in at 10 pt, the status box over the steps); titles "*ADARB2*: The Flag Comes From T7 Calls; Genotypes Miss 32
  Carriers" / "…; Genotypes Hold"; verdicts end "Kept on record with its source"; conclusions state what the reads
  show. Figures: calls open, reads filled (the deck's open = original, filled = check), key "calls ○ | reads ●", x
  title "Calls and reads, 95% CI", dashed reference in gray 1.0 pt, the SIG band across the whole row (and on panel b's
  T7 row), panel titles "Carriers per 100 genomes" / "Frequency of allele A", panel b "Does HaplotypeCaller reconstruct
  the deletion?". Captions explain the dashed line and the bold row.
- Deck-wide: gene symbols with digits fully italic in the matrix, ds_zscale and array_genes (two-artist labels; mathtext
  kept digits upright); slide 13 title and JAK2 lead-in italic; loo_regions annotation without its gene line; array_z
  panel b title centered (legend inside panel b); fig_q2 tick "NovaSeq"; footers completed per slide (CHB and JPT added
  to the glossary); serial commas; "2-field"; "missing genotypes" on slide 9; "logistic mixed model" in slide 5's Say;
  "Six Possible Sources of Bias, Each With a Prediction and Its Test"; If-asked answers open with the answer; cross-
  references quote titles ("{{title:looPlat}}", "{{title:looSrc}}"); SD ranges of the conditional model now name their
  scope (GWAS leads and DQB1:57:V; gene sets 0.51–0.74 and 0.00–0.40).
- Summary: parallel answers ("Not for any of the 20 signals it can rule on"), ADARB2's verdict once in the table, DQB1:57:V
  "P / GEC threshold", "fewer control carriers than ToMMo", "on record", N / I / F in the footer; next steps add the
  100-bp re-call test; "Done" adds HaplotypeCaller's diagnostics; the T7 indel caveat goes to the summary notes (the
  caveats box has no line left).
- Left for the user (not applied): LINC01146 "no artifact" or "no calling artifact" in the summary (decision 21 wording);
  "ADARB2 can be dropped as an artifact" in the summary's paper proposals; the cover date (2026-11-06); AGP3K's expansion
  for the samples footer.
- Checks: deck_check 1,816 numbers, 0 not found, 0 text-rule or title-width issues; check_layout 0 collisions, 0 shape
  problems; validate passed; verify.sh ALL PASSED; common.py mtime 2026-10-30 17:30:54. Old against new text: every
  number of the old slides 25–28 is still on the slide or in its notes, except the ACC leave-out values now signed to the
  pooled effect (4.35 / 4.26 / 4.31 for −4.35 / −4.26 / −4.31) and the F "weakest source" (HTLV-1), which the ACC
  sentence had mislabeled.

Round 11 (2026-11-07; user: "1. 26，27再直观，简洁一些。只用full mainland来举例。2. 28页直观，简洁一些，就回答4个问题就行").
Text and layout only; no number, rule, reading or figure changed.
- adarb2 and linc: every number on the slide is full mainland (strip "Discovery (F): P …, suggestive"; |z| from all
  cases to without T7, with the expected value; alleles, carriers and concordance in F). Each step shows its question
  on its own line (navy, 11 pt) and the answer below (12.5 pt, bold answer first); the four boxes share one height. The
  verdict boxes are folded into the conclusions ("*ADARB2* is a T7 calling artifact: …"; "*LINC01146* is not a calling
  artifact: …"). linc step 1 says that in F the expected z (2.59) sits just under the testability bar of 2.8 and that
  the flag is in N and I. The three-cohort numbers the steps gave, the arrays' missing marker for the indel and the
  symmetric strata of LINC01146 are in each slide's notes (Numbers, first lines). Guards: the F reading of ADARB2 is a
  flag of a suggestive lead; the F reading of LINC01146 is below the bar, the 89 lean with T7, and the 31's P does not
  survive the 64 comparisons.
- summary: one table, Question | Signals | Answer; the Signals column gives the count large (Q1 0 of 21 flagged, Q2 2 of
  20, Q3 0 of 16 contradicted, Q4 0 of 20; NFS 28 pt, column 2.6 in so check_layout's largest-run estimate fits the
  1.2-in rows). Q3 answers "Only for *ADARB2*" (the arrays contradict none of the 16; the deletion has no array marker,
  and its reads show the T7 calls missed every carrier). The flags / not flagged / numbers / caveats boxes and the
  former answers in full are the first items of the notes' Numbers; the Say follows Q1–Q4, then what holds and the
  evidence levels; the open-decision sentence of the Say is unchanged (user decisions pending).
- Checks: deck_check 1,791 numbers on 31 slides, 0 not found, 0 text-rule issues; check_layout 0 collisions, 0 shape
  problems; validate passed; verify.sh all passed; common.py mtime 2026-10-30 17:30:54.
- Label (2026-11-07, user): adarb2 panel a is titled "Carriers (% of samples)" (was "Carriers per 100 genomes"; the same
  numbers: carriers in a group ÷ its samples × 100); only cram_ADARB2.png redrawn (make_figures.cram_adarb2), and the
  figure line of the adarb2 notes uses the same words. Checks as above (1,791 numbers, 0 not found; layout 0 / 0;
  validate and verify.sh passed).
- Label (2026-11-07, user): adarb2 panel a x title "Calls and reads, 95% Wilson CI" (the interval was always Wilson's,
  make_figures._wilson; linc keeps "95% CI", a Wald interval of the allele frequency). Only cram_ADARB2.png redrawn;
  checks unchanged (1,791 / 0; layout 0 / 0; validate and verify.sh passed).
- Interval (2026-11-07, user): linc now uses Wilson's interval too (make_figures._wilson_af: the allele count p × 2n,
  checked to be whole, then _wilson; the Wald helper is gone); x title "Calls and reads, 95% Wilson CI" on both flag
  figures. The intervals move by at most a few thousandths (G400RS 15x calls 0.392–0.641 → 0.394–0.636); points and
  numbers unchanged. Only cram_LINC01146.png redrawn; checks unchanged (1,791 / 0; layout 0 / 0; validate and
  verify.sh passed).
- Summary Q2 (2026-11-07, user): the lead-in "Not for 18 of 20." removed; the answer opens with "Leaving out the
  DNBSeq-T7 cases flags two: …" (the count column gives 2 of 20). Checks: deck_check 1,790 numbers, 0 not found; layout
  0 / 0; validate and verify.sh passed.

HaplotypeCaller census (2026-11-07; user: "我们有没有什么办法取样取的更加系统有说服力，另外把样本ID也可以表上"; choices: every
carrier plus random non-carriers, no read-length test for now, sample IDs also on the slides).
- cram_check/hc_census: the pilot's HaplotypeCaller command on every full-mainland carrier in the reads or the calls
  (292) and 20 random non-carriers per platform (seed 20261107), 392 genomes. The re-run reproduces every call, and the
  deletion is assembled exactly where it was called (392 of 392): called carriers 190 of 190 (150-bp reads), missed
  150-bp carriers 0 of 69, T7 carriers 0 of 32, non-carriers 0 of 100; k = 85 throughout.
- Slide 26: panel b redrawn from the census (share reconstructed per group, 95% Wilson CI; the T7 row banded); step 4
  "in all 190 called carriers and in none of the 32 T7 carriers"; the caption names the census; the Say, Numbers and the
  If-asked answer on the controls' misses use the census; the pilot's paragraph stays in the notes, marked as a pilot
  chosen by hand. The Say lost one sentence (the 150-bp call rate, still in the notes) to stay within 2 min.
- New appendix slide apHcIds (32 slides: 28 main + 4 appendix): every case carrier by sample ID (32 T7 carriers in two
  tables by signature reads; the 9 case carriers on 150-bp reads), with signature reads, whether the deletion is
  reconstructed, and the G>A written instead (GQ); the control and non-carrier counts below. It reads the IDs from
  cram_check/hc_census/hc_census.tsv (git-ignored), never from stats.json. deck_check.py allows sample IDs on that slide
  only (any ID elsewhere fails) and takes its per-sample numbers from the census table. The deck PDF is now listed as
  never-committed in the repository's .gitignore, since it shows the IDs.
- Checks: deck_check 1,827 numbers on 32 slides, 0 not found, 0 text-rule issues; check_layout 0 / 0; validate passed;
  verify.sh all passed; common.py mtime 2026-10-30 17:30:54; no census ID in any file git would add, except a DQB1:57:V
  trace ID that stats.json already held (v57_traced; reported to the user).

Round 12 review (2026-11-07; four read-only reviewers, A numbers, B logic, C text, D visual, and one integrator; plan
and log in the session scratchpad, review12/reports/{A_numbers,B_logic,C_text,D_visual,E_integrated,E_applied}.md). No
rule, reading, flag, pipeline table or stats.json key changed (collect_stats.py not re-run); every number on the slides
and in the notes is in stats.json (deck_check). Items that need the user are listed in E_applied.md.
- adarb2: the census counts every called carrier (191, 186 of them HiSeq X; one control is called with 3 signature reads,
  below the bar of 4, and is reconstructed too); step 4 compares like with like, the carriers in the reads on both read
  lengths ("190 of 259 carriers on 150-bp reads and in none of the 32 on T7"); panel b's rows are those read carriers
  ("Carriers called / missed, 150-bp reads", T7 carriers, random non-carriers), drawn in ink except the T7 row (case
  red), titled "b Deletion reconstructed (% of samples)" with the x title "HaplotypeCaller re-run, 95% Wilson CI". Title
  "*ADARB2*: The Flag Comes From T7 Calls, Which Miss All 32 T7 Carriers"; step 1 says 2.82 is just above the bar;
  step 2 "the cases without T7". The If-asked on the signature says it is a reference unit with one G>A (a G>A in a
  full-length unit gives the same 16 bases) and gives what points to the deletion. Notes: the census caveats (45 of 100
  non-carriers assembled at the lead; one missed case carrier without an assembly region; k = 85 in every assembly),
  the G>A also in 1 of the 4 missed case carriers, "nearly as much read support" (8 against 9), "described instead".
  Guards: panel b's first two rows equal the 150-bp read carriers; the called carriers equal the census's called plus
  the one below the bar; the full-mainland z_exp sits within 0.1 above the bar.
- linc: step 4 "The 89 … lean toward T7 (OR 1.23); the 31 differ from the rest of the cases at P 0.008, yet Q does not
  flag (P 0.023)": "chance" now rests on the test fixed in advance, and the 64 comparisons are defined in the notes
  (16 leads × 4 platforms); the notes give the 89's 95% CI (0.90–1.67) against T7's 1.61; step 2 "the 120 cases
  without T7"; step 3 "the 299 T7 cases they hold". Guards updated (the 89 below T7 and nearer no effect; Q P above
  0.0029 in every cohort; the GLS P_diff equals the pipeline's).
- Both flag figures: shaded rows in regular weight, the captions key "shaded: the cases the flag comes from; bold: the
  flag's reading"; one SIG band extent per figure (_band_rows); row label "T7 cases"; strip at 10 pt, centered.
- summary: the counts in navy (red is never a flag); Q1 "the largest fall in a signal's replicate-median |z| is 0.22
  … a flag needs a fall of more than 1"; Q3 "the reads show that the T7 calls missed all 32 T7 carriers"; labels
  "Q1 · Depth". The Say follows the slide (no "not for 18 of 20"; the T7 reading's limit, five or six of 21, said
  instead; Q3 and Q4 open with their lead-ins); DQB1:57:V near its GEC threshold in I and F only. Notes: P ÷ GEC; OSBP
  not testable without Nagahama; the third platform caveat removed; "Done: … the HaplotypeCaller census (392 genomes)".
- dsSignals: title "Q1 · Depth: At the Controls' Depth the Largest Fall in Median |z| Is 0.22" (single replicates
  fall further; the rule reads the median), rows "Largest fall in median |z|", conclusion "no median |z| falls by
  anywhere near 1"; ds_zscale panels lettered a–c.
- platDesign: the table has Q4's "Flag if" column (multi-line cells: grid() takes an array of lines); "re-tests the
  remaining cases"; "The pooled set and eight subsets"; Say in the table's order; notes "Hartung–Knapp" and the RE
  95% intervals; gene-set SD "up to 0.40, undefined where …" (slides 11, 12, 25).
- arrayDesign: conclusion "with an array marker"; title "… Are Re-Tested …"; the If-asked on two array products says
  a HiSeq X artifact could recur at the 14 imputed leads (only FTO, typed, and DQB1:57:V, HIBAG, re-check the
  controls); q3_flow without the software name. arrayRes: new If-asked on the array z below the WGS z at all 14
  imputed leads (F_arr); the PHF20L1 sentence; concordance in %; panel b's F_arr filled. arrayGenes: "not the cases'
  reads", "in the reference panel".
- Also: slide 2 "false signal"; slide 4 "Every reading is made in all three cohorts", HLA "mixed model; GEC"; slides 5
  and 29 conclusions say what the gates exclude; ctrlDesign/looSrc "a median 0.89–0.90"; looPlat flag points filled,
  "as in N and I, where it is a flag"; platRes INPP5D MAF and convergence wording and a DNB-chemistry If-asked; forests
  "OR (95% CI)", ADARB2/ACC faded on slide 24; hlaCompare "nearly the same", 57 V rows bold; arrayMatch ADARB2 row SIG,
  CDH4 88.2% bold; apGates "(three cohorts together)", "MD5 checksums", "100.0%"; apTbl flags in slide order; apHcIds
  "186 of 186", non-carriers 20 per platform, the shading and sort keys, conclusion "5 of 9 … and 0 of 32 on T7";
  curly apostrophes on every slide face (typo); "vs.", "re-fit", "Q1 · Depth", "3 of 3 cohorts"; sentence-case axis
  titles; footers completed (T7, nt, z_exp, Q, D, IPD-IMGT/HLA, P group).
- Not applied (reasons in E_applied.md): hand breaks on slide 12 (pptxgenjs cuts "z" from its subscript; the box would
  pass 6.3 in); "on the arrays" and "(all but JAK2)" under the summary counts (the count cell holds 20 characters at
  28 pt for check_layout); the word joiner (missing in the PDF font); the matrix legend (decision 20).
- Checks: deck_check 1,853 numbers on 32 slides, 0 not found, 0 text-rule issues; check_layout 0 / 0; validate passed;
  verify.sh all passed; common.py mtime 2026-10-30 17:30:54; no census ID in any file git would add (351 files), except
  the DQB1:57:V trace ID that stats.json already held (v57_traced; stats.json unchanged this round).

Round 12 user decisions (2026-11-07; user: "1～10 按照你认为最合理的开始修改", the ten open questions of the round 12
review, each decided by the most reasonable option):
1. Summary: no answer opens with a bold lead-in any more ("No." / "Only for *ADARB2*." dropped from Q1, Q3 and Q4, as
   the user had dropped Q2's); the count column gives the yes or no. The Say keeps its spoken answers.
2. Slide 27 If-asked on the recruiting institution: the registration-study claim stays (within each study, in every
   cohort, the T7 cases carry A more often than the G400RS 15x cases; guarded), the per-study frequencies and counts
   are left out (decision 13: cases are not split by study).
3. LINC01146's paper wording stays "no calling artifact" (the Say), as slide 27's conclusion says; no slide says "no
   artifact" any more.
4. The summary Say's proposal for ADARB2: "a DNBSeq-T7 calling artifact; reported, not a candidate locus" (was "can be
   dropped as an artifact"); the opening sentence lost three words to stay within the 3-min Say.
5. Cover date 2026-11-07, the date of this version (`collect_stats.py` CONST deck_date).
6. AGP3K: no file on record expands it (the source deck in JHRPv1 is not readable), so the footer of slide 4 defines it
   by what it is: "the group's whole-genome release of the control projects, all on HiSeq X; also the arrays'
   imputation panel" (glossary entry AGP3K).
7. stats.json holds no sample ID: the DQB1:57:V trace keys (`v57_traced`, `v57_trace`, `v57_traced_platform`,
   `v57_changed`, `v57_ambiguous`) carry the pseudonym "traced case 1"; collect_stats.py now refuses to write
   stats.json if any sample ID of the sample table or cram_scan.tsv occurs in it. A scan of every file git would add
   (352) against all 3,101 IDs finds none.
8. Summary Q2 on the slide: "Leaving out the DNBSeq-T7 cases, testable for five or six of the 21 signals, flags two: …"
   (guarded), so that "2 of 20" is not read as 20 signals re-tested without them.
9. Slides 6–7 terms, in `../scripts/design_figures.py` (not frozen): "Re-test discrepancy", "the re-test", "missing
   genotypes", "people", and the workflow's "the group decides what to do with it"; slide 6 Say "T0 concerns the
   re-test itself"; slide 7 footer defines Tier 1 and Tier 2. `design_figures.py` re-run (only master_workflow.png and
   one DESIGN.md line changed; `--check-doc` up to date).
10. One further analysis run, the cheapest and most direct: `cram_check/tommo_signature.py` (ToMMo 61KJPN, allele counts
    only) → stats key `cram.tommo_sig`. The deletion has AF 4.73% (9.2% carriers expected), the G>A at 1747168 0.036%
    (0.07%): the 9.4%–10.0% of samples with the signature fit the deletion, not the G>A (slide 26 If-asked "Is the
    signature specific?"; guarded). Not run: the read-length test (user: "先不做"), template lengths, HaplotypeCaller on
    the 1–3-read genomes, case-imputation r², person-level intervals, allele balance, a platform-label permutation.
- stats.json changed only in: const.deck_date, cram.tommo_sig (added), design_threats (terms), the v57 pseudonyms.
- Checks: deck_check 1,855 numbers on 32 slides, 0 not found, 0 text-rule issues; check_layout 0 / 0; validate passed;
  verify.sh all passed; common.py mtime 2026-10-30 17:30:54.

Round 13 readability (2026-11-07; user: "调用两个agent，这次两个agent对于展示的排版，图片，文字等进行优化（但是不更改我们的既定主题方式，和其他ppt在整体保持一致），让图片文字中的重要必要内容能让人更加方便看，不会模糊不清，混杂"; two read-only reviewers, F figures and G layout and text, and one integrator; plan and log in the session scratchpad, review13/reports/{F_figures,G_layout,H_integrated,H_applied}.md). Readability only: no claim, number, rule, reading, title, theme element (fonts, palette, title and conclusion style, footer, component look) or user-decided design changed; no slide added, removed or reordered.
- Figures (`make_figures.py`; only these redrawn): forests (14, 24) drawn at their slot's size (placed at 0.985 of their drawn size, was 0.92 / 0.93; Q P, I², legend 9.2 → 9.85 pt), intervals 1.4 pt, every estimate drawn above the pooled diamond, which has a white edge (ms 5.6, its outer size as before); `loo_scatter` (13, 23) gene sets filled (every point is a check, so the only open marks are the readings below the bar), the bar lines and leaders in GREY_ (0.8 / 0.7 pt), and the JAK2 label lifted clear of the 1.96 line; `background_depth` (2) the controls' 17.0x in GREY_ (the mark stays the controls' grey); `loo_regions` (12) and `array_z` (20) upright subscripts, as on the slides, and panel b's ticks −3, −1.5, 0, 1.5, 3; `q3_flow` (15) one type size in the five step-3/4 boxes; `cram_ADARB2` (26) panel b to 110, so the 100% point is clear of "190 of 190".
- Workflow (7): `../scripts/design_figures.py` gives the four "Readouts in the matrix" boxes one body height (text at the top), so their tops and the arrows into the matrix line up; only master_workflow.png and its .dot changed; DESIGN.md unchanged (`--check-doc` up to date). `collect_stats.py` re-run strict: stats.json changed only in `design_figs.master_workflow.height_in` (4.65 → 4.707).
- Slides (`build_slides.js`): subscripts written at PowerPoint's own offset, −25% (pptxgenjs writes −40%, which made every line holding a subscript about 3 pt taller and stepped slide 12's underline down; the group's decks use −25%); slide 2 no-break spaces keep "[1, 2]" and each citation with its last word, and the Q4 cell keeps "(178–1,298 genomes each)" whole; slide 3 even gaps (aim → grid → stages → conclusion 0.10 / 0.13 / 0.14 in); slide 6 threat table at 11 pt (rows 0.72 in); slides 8 and 12 subheads at 11.5 pt; slides 9 and 10 the table 0.1 in under the "If an artifact" line; slide 10 the answer row "0 of 21" bold; slides 11 and 22 the schematic captions as figure captions (10 pt, under the schematic), the shared-group tags 9.5 pt, both boxes 12 pt; slide 12 conclusion "… where the testability bar (2.8) is met …" (same claim, no subscript under the underline); slide 13 "z (z_exp), nt" / "chance …" in two lines in every flagged cell, "below expectation;" / "P_diff 0.11–0.13" broken by hand, Cohorts and Chance centered, one lead-in "JAK2 (burden z 1.2, nt), described:"; slide 16 both tables centered; slide 20 the right column at 10.5 pt; slide 21 "What a confirmed carrier shows" and "Controls" (with its 99.6%) moved to the notes word for word (Numbers, first items), the two other paragraphs at 12 pt, the caption's "(those of N_arr and I_arr are among them)" cut (the box says it); slide 29 both gate tables get row heights from their estimate, so they render where they are placed.
- Not applied (reasons in H_integrated.md): F-09, O-01, O-03, O-05, O-07. For the user: larger sub/superscripts (deck-wide, footers or the N_arr headers; the group's decks keep the default size) and drawing the evidence matrix at h = 5.15 in (text 9.15 → 9.5 pt; decision 20).
- Checks: deck_check 1,855 numbers on 32 slides, 0 not found, 0 text-rule issues; check_layout 0 collisions, 0 shape problems; validate passed; verify.sh all passed; common.py mtime 2026-10-30 17:30:54; no sample ID of the sample table or cram_scan.tsv in any file git would add.
- Lead's check after the round (2026-11-07): slide 21's caption had lost the clause "(those of N_arr and I_arr are
  among them)"; it is restored in the slide's notes ("The F_arr cases include those of N_arr and I_arr (the cohorts are
  nested), …"). Independent rebuild: deck_check 1,855 numbers, 0 not found, 0 text-rule issues; check_layout 0 / 0;
  validate passed; verify.sh all passed; `--check-doc` up to date; stats.json differs from the pre-round state only in
  `design_figs.master_workflow.height_in` (4.65 → 4.707); no sample ID in any of the 352 publishable files.

Slide 14 conclusion (2026-11-07; user: "你只要在下面说明白，虽然我们的规则没有在这里flag这个，但是我们明显发现了G400RS 15x在这个位点上的
不同（之所以没有被flag，是因为这个平台怎么样）"): the bottom line now reads "The rules flag nothing here, but *LINC01146*'s 31
G400RS 15x cases clearly run the other way (OR 0.73; pooled 1.43): so few cases give a wide interval (0.44–1.23), and Q
stays above its threshold (P 0.013–0.023)." (two lines; `conclusion()` takes an optional height). It replaces "Among
common markers no platform stands apart; the rare leads are judged by leaving platforms out.", whose rare-lead clause
the box above already states. Guarded: no platform_Q or platform_RE flag in any cohort, LINC01146's Q P above the
threshold in every cohort, the G400RS 15x OR below 1 with its upper bound below the pooled OR, 31 cases. The Say names
the 31 cases and why Q does not flag them (128 of 135 words). Title unchanged. Checks: deck_check 1,861 numbers, 0 not
found; check_layout 0 / 0; validate passed; verify.sh all passed.
