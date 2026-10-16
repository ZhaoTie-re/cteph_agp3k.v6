# Outputs

`results/latest.json` points to the newest complete snapshot (`snapshot`, `dir`, `knowledge`, `workbook`, `release`, all `snapshots`).

`results/snapshots/<FetchDate>/` holds one snapshot:

| Path | Content |
|---|---|
| `00.raw/` | Open Targets responses for CTEPH: `meta.json`, `disease.json`, `associations.json.gz`, `evidence.json.gz` (text-mined and GWAS), `targets.json.gz` (Reactome themes), `gwas.json.gz` (credible sets with colocalisations), `fetch_cteph.json` (counts) |
| `00.raw/diseases/disease.<id>.json.gz` | one related condition: record, GWAS studies and credible sets, curated / Mendelian evidence |
| `00.pubmed/seed_pubmed.json.gz` | the PubMed records the seed cites |
| `00.pubmed/pubmed_search.json.gz` | the records of the CTEPH-genetics query, with the query and count |
| `00.pubmed/fulltext.json.gz` | the Methods text of the open-access full texts |
| `01.conditions/` | `conditions.tsv`, `relations.tsv`, `colocalisation.tsv`, `seed_check.tsv` (every seed quote found in its abstract), `fetch_list.tsv` |
| `02.literature/` | `screening.tsv` (every record: included or not, and why), `studies.tsv` (with `methods_source`: full text or abstract), `findings.tsv`, `methods.tsv` (field, value, quote, source), `extract_log.tsv` (items dropped by grounding), `extract.json` |
| `03.tables/` | the long tables and the graph (below), `tidy.json` (counts) |
| `04.summaries/` | `summaries.tsv` / `.json`: kind, theme, text, grounded, problems, inputs |
| `05.changes/` | `changes.tsv` (category, change, key, item, detail), `changes.json` |
| `06.knowledge/` | `cteph_knowledge_graph.md`, `cteph_knowledge_graph.json`: scope, then one section per question of the workbook (2–8), changes, caveats |
| `cteph_knowledge_graph_<FetchDate>.xlsx`, `workbook.json` | the workbook (below); its layout and rows per sheet |
| `_run_info/` | trace, report, timeline, DAG; `run_manifest.json` (parameters, seed and script md5s, release, counts, changes) |

## 03.tables

| Table | One row per |
|---|---|
| `conditions.tsv` | condition: category, identifiers, evidence sources, genetics mode, clinical items, CTEPH loci shared, GWAS studies (East Asian), credible sets, Mendelian genes |
| `relations.tsv` | piece of relation evidence: type, relation / direction, statement, verbatim quote, design, population, East Asian, PMID, citation, strength (CLPP), CTEPH loci |
| `colocalisation.tsv` | CTEPH credible set × colocalising GWAS credible set: trait, study, ancestry, lead variant and P, CLPP / H4, method |
| `gwas_studies.tsv` | GWAS study of CTEPH or a condition: trait, publication, sample, cases, controls, cohorts, ancestry, EAS flag, summary statistics, credible sets |
| `gwas_loci.tsv` | credible set: lead variant, rsID, P (mantissa, exponent, log10), β, fine-mapping, confidence, top three L2G genes and scores |
| `mendelian.tsv` | condition × gene × source (curated gene evidence): source score, records, classifications (ClinGen / Gene2Phenotype class, PanelApp rating, Orphanet assessment, UniProt confidence), ClinVar review status, clinical significance, specific diseases, PMIDs |
| `mendelian_records.tsv` | curated / Mendelian evidence row |
| `cteph_targets.tsv` | CTEPH target: score, literature score, publications, years, Reactome themes |
| `literature.tsv` | text-mined target × publication: title, matched words, sentences |
| `cteph_studies.tsv`, `cteph_findings.tsv`, `cteph_methods.tsv`, `cteph_screening.tsv` | the literature extraction (as in 02.literature) |
| `gene_condition.tsv` | gene × condition × evidence kind (GWAS, Mendelian, Literature): records, details, East Asian |
| `nodes.tsv`, `edges.tsv`, `graph.json` | the knowledge graph (node: id, label, type, category; edge: source, target, relation, detail, evidence, strength) |

## The workbook

Organised by question. Each question sheet opens with its answer in brief, then a nested table (Excel row outlines): open a row with + at the left, or a whole level with the buttons 1 2 3. Links join the same gene or condition across sheets; cells with a note hold the source passage or the evidence behind a code.

| Sheet | Question | Rows: group › record › source |
|---|---|---|
| Guide | — | the questions with their answers in brief, how to read, the snapshot, caveats |
| 1 Related conditions | Which conditions are related to CTEPH, and on what evidence? | category › condition (evidence counts, CTEPH loci shared, GWAS, curated genes, genes shared with CTEPH) › evidence (quote, ontology relation, shared GWAS signal) |
| 2 Shared genes | Which genes do CTEPH and its related conditions share? | gene × condition map: G GWAS, C curated, L CTEPH literature; genes with CTEPH evidence, in ≥ 3 related conditions, in 2 |
| 3 CTEPH genetics | What is known about the genetics of CTEPH? | the CTEPH GWAS loci; gene theme › gene (studies by result) › finding with its source sentence; generated summary of the findings |
| 4 Study methods | How were the genetic studies of CTEPH done? | design › study, its methods by step (samples, genotyping / sequencing, analysis, replication); tallies; generated summary of the methods |
| 5 Condition genetics | What are the genetics of each related condition? | category › condition › gene (GWAS loci, best P, L2G, curated evidence, in CTEPH) |
| 6 Mechanisms | Which mechanisms does the CTEPH literature describe? | Reactome theme (generated summary) › text-mined target › publication with its sentences |
| 7 East Asian | What evidence comes from East Asian and Japanese populations? | kind › item: clinical evidence, CTEPH genetic studies, GWAS, colocalising signals |
| 8 Changes | What changed since the previous snapshot? | one row per change |
| A1 GWAS studies … A6 Screening | — | flat tables with filters: GWAS studies, GWAS loci, colocalisation, curated genes, CTEPH studies, the screening of every PubMed record |
