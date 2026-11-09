# deck — HLA typing in JHRP v6

**One deck: `hla_typing.pptx`, 30 slides.** It follows the order of a paper: the question and aim first, then the
background a reader needs, the method, the quality of the typing, the preparation of the markers, and a summary;
a Notes section after the summary holds the detail a question may need.
The background is written for readers new to HLA, one figure per concept. The quality evidence is kept to two
checks (completeness; concordance of control allele frequencies with ToMMo 61KJPN). The summary closes the deck and
states the methods sentence for a manuscript.

**Outline.** Slide 3 names the five parts of the talk (Background, Method, Quality, Markers, Summary), the aim
of each, and the evidence it presents, without slide numbers; the PowerPoint sections use the same names.

| part | slides |
|---|---|
| Cover | 1 — the CTEPH cover template (disease name, "(CTEPH):", subtitle), author, date |
| Introduction | 2 — HLA associations with disease and the alleles reported in Japanese CTEPH; this study (cases, controls) · 3 — aims and outline: the question, and the aim and evidence of each part |
| Background | 4 — function of the HLA molecules (class I / II, CD8 / CD4 T cells) · 5 — structure of the HLA molecules and the site of allelic variation (crystal structures, panels a–c; variability of the HLA-A groove) · 6 — the HLA region on chromosome 6 · 7 — inheritance of the haplotypes · 8 — allelic diversity (alleles per gene; growth across IPD-IMGT/HLA releases) · 9 — allele nomenclature, 2-field resolution, P groups · 10 — DR haplotypes and *DRB3/4/5* |
| Method | 11 — typing from sequencing reads (a real 48-base stretch of HLA-A exon 2; a look-alike gene) · 12 — the HLA-HD algorithm, and when a 2-field name is written (candidates agree / frequency table / one P group; otherwise none, never an arbitrary pick) · 13 — validation of *DRB3/4/5* copy number against *DRB1*: one real case per outcome, counts, why *DRB1* is decisive · 14 — the two criteria by which typing quality is assessed · 15 — how our controls' calls and the ToMMo table are made comparable, as a table of five steps (comparable population, genes, resolution, frequency definition, dictionary release) with the bias each step prevents (`collect_stats.py` `harmonization`) |
| Quality | 16 — completeness of the typing (results table with denominators) · 17 — presence of the called alleles in ToMMo, in the controls, at P-group resolution: (a) by distinct allele, (b) by copy, side by side with a pooled "All" row (35.9% of distinct alleles, 97.8% of copies) · 18 — concordance of allele frequencies with ToMMo 61KJPN, the one off-diagonal allele explained · 19 — the one dot off the diagonal, DRB1\*15:13: what, why (renamed DRB1\*15:13Q in 3.59.0, no P group, protein = DRB1\*15:01, not in HLA-HD's frequency table), where (the DRB1\*15 group on both sides, `fig_outlier_group`; the DRB1 panel enlarged, `fig_outlier_panel`) |
| Markers | 20 — amino-acid residues, on real HLA-A sequences: (a) why residues are tested as well as alleles (three allele groups share Y at position 9), (b) how a residue column is filled (amino acids read from the IPD-IMGT/HLA protein alignment, not the reads; one column per amino acid seen; 0/1/2 copies), (c) when a cell is blank and why not 0 (a person in the cohort with A\*24:82, unsequenced at 184) · 21 — the 0/1/2 encoding of alleles and residues, a real slice of the full mainland cohort's final matrix · 22 — call rate (the first floor) · 23 — minor allele frequency (the second floor, on what the first leaves) · 24 — stepwise quality-control record, three cohorts (typed → *DRB3–5* excluded → call rate → MAF) · 25 — allele columns retained at each step, by gene, and the share of the typed copies the tested alleles carry (the rest, alleles below MAF 0.5%, enter through the residue columns) · 26 — residue columns, the same (each with its conclusion line) |
| Summary | 27 — four headline results (the only slide with large figures), the two parts side by side with a caveat row (the ToMMo concordance checks P-group resolution only), the methods statement |
| Notes | 28 — how our alleles map onto the ToMMo table at P-group resolution: a worked example (DQA1\*03:01P) and, per gene, ToMMo's P groups and alleles in no P group, our 2-field alleles, how many merge into how many P groups, and how many result (`collect_stats.py` `pgroup_mapping`) · 29 — one gene (DQA1) through the comparison: all 24 alleles of either table at P-group resolution, 8 in both, 15 only in our controls, 1 only in ToMMo, with the 2-field alleles merged into each and a scatter drawn as on the frequency slide (`worked_example`, `fig_worked_example`) · 30 — the tested allele columns that share a P group: ToMMo's entry and frequency beside our alleles in it, with their frequencies in the controls and the sum (`pgroup_splits`, `pgroup_mapping`) |

The further checks (reads behind each call, confidence classes, DRB reconciliation, 1KG JPT, locus tiers) stay in
`../results` and `../results/figures`; they are not in the deck.

## Sources

Background, each statement checked against the full text (PMC) or abstract:
- Medhasi & Chantratita (2022): function and structure of class I / II;
- Trowsdale & Knight (2013): region size, gene content, diversity, disease associations;
- Choo (2007): inheritance, DR / DQ / DP chain pairing, laboratory typing methods, disease examples;
- Chung et al. (2004): carbamazepine and HLA-B\*15:02;
- Barker et al. (2023): the IPD-IMGT/HLA database;
- Kawaguchi et al. (2017): HLA-HD and the HLA-Y pseudogene;
- Marsh et al. (2010): nomenclature;
- Tanabe et al. (2005), Kominami et al. (2009): HLA alleles reported in Japanese CTEPH; numbers checked against the full texts in `../ref` (carrier frequencies and ORs from Table 2 of each, Tanabe's 95% CIs from its Results text, all as reported and none recomputed; Kominami's 160 patients include Tanabe's 80);
- Madden et al. (1993), Stern et al. (1994): the crystal structures on slide 5 (PDB 1HHK, 1DLH);
- Sakaue et al. (2023), full text PMC10786448: allele, residue, omnibus and conditional tests; rare alleles;
- hla.alleles.org: P groups and DRB haplotypes.

Data drawn on the background slides:
- **Structures** (slide 5): `../assets/pdb/1HHK.pdb`, `../assets/pdb/1DLH.pdb` (RCSB, downloaded 2026-10-19), drawn as
  C-alpha traces by `molview.py`. Variability = Shannon entropy of the amino acids at each position across the
  9,175 HLA-A alleles of the IPD-IMGT/HLA 3.64.0 protein alignment (`../results/05.exports/residue_reference/A.tsv`).
- **Matrix slice** (slide 21): ten people drawn at random from the full mainland cohort (fixed seed), the five most frequent
  alleles and two whole residue positions of *A, B, DRB1, DQB1*, among the columns the association tests. After
  the alleles, a grey *other* cell (not a column of the matrix) holds the person's remaining copies at the gene —
  alleles with a column not drawn, or below the MAF floor and without a column — read from the genotype in
  `05.exports`, so the alleles drawn plus *other* add up to 2 for every person and gene. A whole
  position shows every one of its columns, so a person's residue columns there add up to 2 as on slide 20 (b);
  only complete positions are eligible (every amino acid seen there kept a column; at the others a carrier of an
  amino acid below the MAF floor adds up to less than 2), the two most variable that fit in six columns, and
  *HLA-A* always shows position 9 of slide 20. `collect_stats.py` asserts the sums.
  SNPs are deliberately NOT shown beside them: they are different markers of the same region, and placing them
  there would suggest a correspondence with the HLA columns that does not exist.
- **Reads** (slide 11): HLA-A exon 2 of A\*24:02:01:01 and A\*11:01:01:01 from the pinned IPD-IMGT/HLA 3.64.0
  nucleotide alignment (`../assets/imgt/A_nuc.txt`, downloaded 2026-10-20).
- **Gene map** (slide 6): positions from Ensembl (EnsDb v86, GRCh38), the table used by `analysis/assoc_hla`; bands
  from UCSC hg38 (`../assets/ucsc/cytoBand.hg38.chr6.txt`, downloaded 2026-10-19).
- **Residue examples** (slide 20): real rows of the IPD-IMGT/HLA 3.64.0 protein alignment for HLA-A, through the exports' allele-to-alignment match (`../results/05.exports/residue_reference/A.tsv`, `allele_match_map.tsv`); panel c is a real person of the cohort carrying A\*24:02 / A\*24:82, shown without an identifier. Drawn by `make_figures.residue_concept()`.
- **The marked dot** (slide 18): the allele of *DRB1* furthest from the diagonal *as the figure draws it*
  (square-root axes), picked by `freq_outlier()` in `collect_stats.py`; everything said about it is read off the
  dictionaries, not typed in: the name in release 3.64.0 and the name ToMMo's release gives it
  (`../assets/imgt/Allelelist_history.txt`), its P group or the absence of one (`hla_nom_p.txt`, 3.64.0), whether
  HLA-HD's frequency table holds an entry (`<HLA-HD>/freq_data/DRB1_count.txt`), the allele it cannot be told
  apart from in the IPD-IMGT protein alignment, and the 2-field group totals on both sides.
- **Allele counts** (slide 8): the pinned IPD-IMGT/HLA allele history (`../assets/imgt/Allelelist_history.txt`);
  the growth panel starts at release 3.0.0, where HLA-C took its present name (it was "Cw" before).

**The slides and the pipeline apply one rule set in one order** — a column counts while it varies in the
cohort (at least one copy of its less common form); *DRB3/4/5* excluded; call rate ≥ 95 % in cases and in controls (platform and depth differ between them); then MAF ≥ 0.5 % — and
the slides read the pipeline's own record of it, so slides 21–27 carry the numbers of the files the association
reads (`../results/06.assoc_inputs/<cohort>/`). A residue IPD-IMGT/HLA leaves undetermined is blank from the export
on (step 05), so it is not a filter step of its own: it removes no column and lowers its position's call rate. A
loss on slide 22 is counted over the columns that vary in each cohort; slide 23 starts from what the call-rate floor
leaves.

## Conventions

- **Counting:** a marker's frequency is counted in **copies** (each person contributes 0, 1 or 2), never in
  "chromosomes"; the word *chromosome* is kept for the chromosome itself (slides 6, 10).
- **Speaker notes have one structure** (`speak()` in `build_slides.js`): *Say* (the talking points, in order), *Numbers* (figures to have at hand, from `stats.json`), *If asked* (expected questions with answers), *Source* (data files and literature); one bullet per line, section names bold. No decision dates or project history in the notes.
- **Cross-references quote the slide title, never a page number** (pages move as slides are added): the referenced titles live in one table, `T` in `build_slides.js`, used both by `title()` and by `ref()` (slide text) / `refN()` (speaker notes); section names (Background, Method, Notes) may be named as such.
- **Resolution is written *2-field*** (never *two-field*), and "resolved to 2 fields"; *the first two fields* stays where the nomenclature counts fields.
- **Markers are columns.** A count of markers always says *allele columns* (one per 2-field allele) or
  *residue columns* (one per amino acid at a protein position), never bare "alleles" or "residues": "679
  residues" would read as 679 positions, when it is 679 columns at 295 positions.
- **Cohorts are named in full:** *narrow / intermediate / full mainland* (`COH3` in `build_slides.js`, the
  order of every "a / b / c" cell), and *full mainland cohort* wherever a single cohort's number is shown; a
  call-rate header always adds *in cases and in controls* (`CR_HEAD`), unless slide 22's definition is beside it.
- **Text style:** all slide text in DengXian; gene symbols in italics, allele names upright; plain tables with
  a navy header.
- **Type sizes** match `../../tmp/ppt/rare_variant_assoc` (pt): title 24 bold underlined; conclusion line and
  in-slide heads 15; body bullets 11.5; tables 10.5 (the two per-gene tables 9.5, the outline table 11); boxes
  12, dense rationale boxes 11; figure captions 10.5; notes 9.5; footer 7.5; slide number 10; the summary
  callouts 24 with 11 labels; cover 34 (the CTEPH template). Set them in the helpers, not per call.
- **Slide layout:** one underlined conclusion per slide; a figure legend as plain caption text under each data
  figure, the way a paper carries it (no tinted panel, and never a "how to read" label); abbreviations and
  references in the footer.
- **Figures:** one visual system — class I blue, class II green, peptides and what must not be missed red,
  reference data grey; drawn at slide size (12.6 in wide) in Liberation Sans.
- **The two floor slides (22, 23)** are drawn the same way: alleles and residues in separate panels,
  because a floor may act on one kind and not the other, and a broken y-axis (`broken_panel`), because one
  class holds most markers and would otherwise flatten the classes the floor acts in. The break marks are
  placed in figure inches, so both halves get the same slope however the halves are divided.
- **Layout checks.** Before a preview is saved, two checks must be clean: no two text lines of the rendered PDF
  overlap, and no shape (tables at the height their wrapped text gives them) reaches into another or past the
  slide edge. The footer places itself from the bottom up, counting about 205 characters to a line at 7.5 pt;
  `build_slides.js` refuses to build when a footer would not fit. Run `check_layout.py` (below).
- **No headline numbers** except on the summary slide. Elsewhere a result is a row of a table with its denominator, not a large figure in a
  tinted box, and nothing on a slide states a value without what it is a share of.
- **Marker numbers are the pipeline's.** Every count on slides 21–27 is read from `../results/06.assoc_inputs`
  (`markers.tsv`, `qc_record.tsv`, `summary.json`), which step 06 writes together with the matrices the
  association reads; `collect_stats.py` recomputes none of them and stops if the record and the matrices
  disagree. Cohorts and floors come from `../scripts/config.py`.
- **Numbers** come from `stats.json` (`collect_stats.py`, read from `../results`) and `figure_facts.json`
  (`make_figures.py`); none is typed by hand.

## Files

| file | content |
|---|---|
| `hla_typing.pptx`, `hla_typing.preview.pdf` | the deck and its LibreOffice render (DengXian replaced by an Arial-metric font, so text widths are approximate) |
| `collect_stats.py` → `stats.json` | every number on the slides, including the association inputs per cohort |
| `make_figures.py` → `figures/`, `figure_facts.json` | the 18 figures; facts read off the data they draw |
| `molview.py` | pseudo-3D backbone drawings from PDB files (matplotlib) |
| `check_layout.py` | the two layout checks: overlapping text in the rendered PDF; shapes reaching into each other or past the slide edge in the .pptx (exit 1 if either finds anything) |
| `build_slides.js` | the generator (pptxgenjs). Markup: `*italic*`, `^bold^`, `~subscript~`; `\` makes the next character literal |

## Rebuild

```
python3 collect_stats.py
python3 make_figures.py
NODE_PATH=<node_modules with pptxgenjs, jszip> node build_slides.js [--apply-theme <pptx skill>/scripts/apply_theme.js]
soffice --headless --convert-to pdf hla_typing.pptx      # preview; rename to hla_typing.preview.pdf
python3 check_layout.py hla_typing.pptx hla_typing.preview.pdf   # must report 0 and 0
```

Use the `cteph_geno_pro` interpreter for the Python steps. LibreOffice 26.x needs glibc ≥ 2.34; on this EL8 system
use 25.8.
