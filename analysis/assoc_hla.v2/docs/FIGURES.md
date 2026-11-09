# assoc_hla.v2 — figures

The figures follow the questions a reader asks of the analysis, in paper order, not the case sets. All of them are in `results/figures/0N.<family>/`. Each family has one `README.md` with:

- the question;
- what each panel plots;
- an interpretation computed from the tables;
- the values drawn and the full statistics;
- how to read it;
- what it cannot establish.

| # | family | file(s) | question |
|---|---|---|---|
| **main** | | | |
| 01 | `landscape` | `bars.<case_set>.png`, `manhattan.<case_set>.png` (4 case sets each) | Where in the MHC are the associated tests, and in how many cohorts? |
| 02 | `signals` | `signals.png` | How large is each independent signal in every case-set comparison, and at which thresholds is it significant? |
| 03 | `conditioning` | `conditioning.png`, `carriers.png` | Does one marker explain each set of significant markers, and which alleles carry each residue signal? |
| 04 | `gwas_link` | `gwas_link.png` | Is the genome-wide MHC lead the same signal as the HLA signal? |
| 05 | `dvt` | `dvt.png` | Does any test differ between cases with and without DVT? |
| **supplementary** | | | |
| 06 | `thresholds` | `thresholds.png` | How strict is each threshold, and how close is each to a family-wise error of 0.05? |
| 07 | `models` | `models.png` | Do the fixed and random models agree test by test? |
| 08 | `qq` | `qq.png` | How far does each scan depart from the uniform null? |

A family with nothing to draw (for example `signals.png` with no signal) writes no PNG; its README says so.

## Two encodings of the landscape

`01.landscape` draws the same tests in two ways:

- **`bars.<case_set>.png`.** One vertical bar per test. It runs from −log10 of the test's largest P over the three cohorts to −log10 of its smallest, so the top is the strongest cohort and the bottom the weakest. A bar wholly above the solid line is significant in 3/3 cohorts; a bar that crosses it with its top only is significant in 1/3. This mirrors the conditional analysis's rule, where the lead is the marker whose largest P over the cohorts is smallest.
- **`manhattan.<case_set>.png`.** The three cohorts as rows, alleles and residues as columns, each cohort with its own thresholds.

## Encodings shared by every figure

**Model.** The random model (SAIGE) is the primary model and the only one in the main figures. The fixed model appears only in `07.models`.

**Cross-cohort consistency at the LD-adjusted threshold (n/3)** is one hue, light to dark, with line width and marker size growing with n:

| label | colour |
|---|---|
| 3/3 cohorts | `#5E1A18` |
| 2/3 cohorts | `#B53F3B` |
| 1/3 ("Single cohort") | `#E3896B` |
| not significant | `#CDD2D7` (recessive grey) |

The colours were checked with the dataviz palette validator: adjacent steps differ by ΔE ≥ 16 under colour-vision deficiency and ≥ 19 for normal vision, and lightness falls monotonically with n. Signal leads are labelled directly, so identity never rests on colour alone.

**Thresholds** are drawn in the same three styles everywhere:

| threshold | role | line |
|---|---|---|
| LD-adjusted (GEC) | primary | solid red |
| exact (max-T) | validation | dashed grey |
| conventional Bonferroni | sensitivity | dotted dark grey |

**Cohorts** use an ordered ramp from dark (narrow) to light (full mainland), abbreviated N / I / F where space is short.

**Tests, not markers.** A residue mirror pair, two columns that are complements, is one test and is drawn once.

**Marker order.** Genes run in GRCh38 chromosomal order, then residues by IMGT position and letter, and alleles by their numeric fields. Each gene gets a minimum band width so its name stays legible. The x axis is therefore not a physical coordinate. A violet triangle marks where the genome-wide MHC lead falls among the genes.

**Titles are computed.** Every figure title that states a result is built from the tables. A title cannot outlive the data it describes.

**Style.** Style, type scale and caption layout come from `analysis/_shared/scripts/plot_style.py` and `figure_doc.py`. Names, orders, colours and labels come from `scripts/vocab.py`. Each figure carries:

- one bold sentence;
- one line per panel;
- a `Data:` line.

The interpretation lives in the family README, not in the PNG.

**λGC** is drawn only in `08.qq`, and is described as a property of the region, not as a calibration of the test (METHODS §5).
