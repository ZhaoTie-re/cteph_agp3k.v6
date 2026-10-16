# Study notes

## Design decisions (2026-10-16, with the user)

**History.** The first build, `integration/cteph_opentargets`, took everything Open Targets holds on CTEPH: targets, drugs, mouse models, publications. The user redirected it. The graph is to centre on:
- CTEPH and the conditions related to it, every relation with its evidence;
- the genetic studies of both, with their methods;
- other important content.

The first snapshot of that build was deleted on 2026-10-16 at the user's request.

**Relatedness evidence.**
- The sources are Open Targets data plus PubMed-verified literature, kept apart by type: clinical, ontology and colocalisation.
- Open Targets literature "co-occurrence" counts are unusable: `literatureOcurrences(additionalIds)` returns the sum of the two entities' counts, not their overlap.

**Categories.** All four, plus the data-driven shared-signal traits:
- precursor / spectrum;
- risk factors and comorbidities;
- intermediate traits;
- PH group.

**Genetic evidence.**
- GWAS credible sets with L2G;
- curated Mendelian sources;
- the CTEPH genetic literature (PubMed), since Open Targets has a single CTEPH GWAS (Liley 2024).

**Methods.** The user asked that the methods of the genetic studies be summarised too, so they are extracted per study, from the abstract and the open-access Methods sections.

**Other content.**
- East Asian / Japanese evidence is flagged throughout and gathered in one section and one sheet.
- Text-mined CTEPH targets are kept with their source sentences and summarised by Reactome theme.
- Drugs and mouse models are dropped.

**Language model.** Summaries were requested, so the screening, extraction and summaries use an isolated, cached model. Everything else is deterministic.

**External only.** Nothing of this project's own results is read or compared.

## What the seed research showed (PubMed, 2026-10-16)
- **Control choice.** Most classic risk factors were found against non-thromboembolic PH controls (Bonderman 2005, 2009; Jaïs 2005; Lang 2013). Against PE survivors (Dodson 2022), malignancy, thyroid replacement, intravascular devices, IBD and osteomyelitis were not enriched; only splenectomy and non-O blood group trended. The seed keeps these contrary results.
- **IBD and osteomyelitis** rest mainly on one Vienna case-control study (OR 67, 95% CI 7.9–8832).
- **Platelets.** The evidence conflicts: a review lists high platelet counts, while measured cohorts found lower counts with more activation.
- **FXI** has genetic evidence only (the F11 locus).
- **Lp(a)** rests on one small 1998 study.
- **Hereditary spherocytosis** has no supporting abstract and is not in the seed.
- **The 2022 ESC/ERS guidelines** have no PubMed abstract, so a focused review of them (Rosenkranz 2023) is quoted instead.
- **HLA naming.** Older abstracts write HLA-B*5201 / DPB1*0202 (now B*52:01 / DPB1*02:02). PubMed renders some as "B(*)5201", so the grounding check ignores parentheses.

## What the build showed
- **Colocalisation.** It is eCAVIAR only (PICS fine-mapping; no coloc H4). The ABO locus colocalises with many unrelated traits (pancreatic cancer, type 2 diabetes, protein levels), so relations are graded by the number of CTEPH loci shared.
- **PubMed records.** The genetics query returns book records (GeneReviews) as `PubmedBookArticle`; they are parsed too.
- **Extraction smoke test.**
  - 6 records.
  - **Correct:** HLA-B*5201 / DPB1*0202, NFKBIL1, MUC6 and the CHIP genes were extracted, with their sample sizes and effects.
  - **Dropped by grounding:** an HLA-DRB1 the abstract does not name, and a "Mutations of factor II" variant that the abstract does not write.
  - **From the full text:** WES / MuTect2 / Holm / HLA-HD methods were read from the Keio Methods section.
- **Summary grounding.** A gene check against the theme's target list was too strict: sentences name other genes. The check is now that a gene named must occur in the source text.

## History
- 2026-10-16: redesigned from `cteph_opentargets`; seed researched and verified; component built. Then:
  - At the user's choice, a deterministic pre-screen was added: 471 → 98 records, all 16 known genetic studies kept. Open Targets has no literature summaries, and its text mining misses and mis-maps the genetic studies.
  - The full-text fetch failed on PMC1729855, an HTTP 500 from Europe PMC. Such failures now count as "no full text".
  - At the user's request, the methods of studies without a full text now come from their abstract.
  - The user found the 16-sheet workbook "still cluttered, without logic, not intuitive" and chose a question-led workbook whose rows open level by level, with the CTEPH literature grouped by gene / pathway theme. It now has a guide, seven question sheets (answer in brief, then group › record › source as Excel row outlines; quotes in cell notes; links between sheets), the changes, and six flat appendix tables with filters. The knowledge base follows the same order.
  - In the same change, "Mendelian" evidence is shown as **curated gene evidence** (it includes gene burden and UniProt), and the ClinVar review status, which Open Targets returns as the record's confidence, is kept apart from the classifications of the other sources.
  - The first build's snapshot (`cteph_opentargets`) was deleted at the user's request.
