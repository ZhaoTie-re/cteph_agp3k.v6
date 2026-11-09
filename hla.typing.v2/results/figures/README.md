# hla.typing.v2 — figures

Journal double-column width (7.2 in), 300 dpi. Drawn by `scripts/07_figures.py` from the tables in `results/`.
Colours follow the project's semantic palette:
- blue: the measurement;
- red: what fails a criterion or must not be missed;
- grey: references and thresholds.

| figure | shows | how to read |
|---|---|---|
| `01.typing_overview.png` | Per reported gene, three stacked bars: (a) call state, (b) how the 2-field allele was set, (c) call confidence. | (a) separates real absence (*DRB3/4/5* "absent", grey) from failure (red). In (b), light blue is HLA-HD's frequency table deciding between alleles the reads cannot separate, mostly *DPB1*. In (c), high confidence is ≥ 95 % at 10 of 19 genes, with *A* (91.7 %), *DRB3* (81.9 %) and *DRB4* (79.4 %) lowest. |
| `02.read_evidence.png` | (a) Minor-allele share of unique reads, heterozygotes at obligate genes. (b) The same for HLA-HD's two-allele *DRB3/4/5* calls, split by the copies *DRB1* implies. (c) *DRB1*-implied vs HLA-HD-reported copies per paralogue. | A real heterozygote sits near 0.5. The red tail in (b) is a second paralogue allele where only one copy can exist: 453 of 660 such calls are below 0.15. In (c) the diagonal is agreement (82 %, 82 %, 95 %). |
| `03.panel_frequencies.png` (**main figure**) | Allele frequencies in our controls vs ToMMo 61KJPN, one panel per locus, square-root axes. Red triangle: allele absent from ToMMo although its dictionary had it. Open triangle: allele named after ToMMo's dictionary. | Points on the diagonal agree. "absent" is the share of control chromosomes on red-triangle alleles; it is bold red where it exceeds 2 % (*B*). *DRB3* and *DRB4* frequencies are not compared (reasons in docs/METHODS.md §5). |
| `04.locus_summary.png` | Per locus: (a) absent alleles, all calls (filled) and high-confidence calls (open); (b) heterozygotes reported as homozygotes, with 95 % CI; (c) high-confidence calls; (d) tier. | Red marks a criterion not met; dashed lines are the thresholds (2 %, 3 %, 95 %). Tier = 1 + criteria failed; loci with no panel start at 2. |
| `05.panel_ceiling.png` | Discordance D with 1KG JPT: our controls (dark) and ToMMo (light), with the D expected from sampling alone (black tick). | With 105 people, sampling alone gives D of about 0.1–0.2. Both our controls and ToMMo sit at or below that, so 1KG JPT can rule out gross error only. |
| `06.confidence_vs_panel.png` | Share of control chromosomes on an allele absent from ToMMo, by call confidence, with 95 % CI. | 0.5 % (high), 5.8 % (medium), 16.8 % (low): the read-level grading is borne out by an external panel that played no part in it. |
| `07.sample_qc.png` | (a) Read pairs given to HLA-HD per sample (log axis), with the 1st percentile. (b) Samples by number of reported genes with a low-confidence or no call. | 2,492 of 3,101 samples have none; 13 have three or more. |
| `08.shift_direction.png` | Alleles whose frequency differs from ToMMo by ≥ 1 point at the five 1KG loci: ours − ToMMo (x) against 1KG JPT − ToMMo (y). | Blue quadrants: the independent 1KG JPT sample differs from ToMMo in the same direction as our controls (24 of 33). Shifts of this kind are consistent with differences between cohorts rather than typing error. |
