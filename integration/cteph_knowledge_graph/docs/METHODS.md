# Methods

## 1. Sources
- **Open Targets Platform.** The GraphQL API, `https://api.platform.opentargets.org/api/v4/graphql`. The release (data and API version) is recorded in every snapshot.
- **PubMed.** NCBI E-utilities (esearch, efetch), without a key; requests are spaced 0.4 s apart.
- **Europe PMC.** REST `fullTextXML`, for the Methods sections of open-access full texts.
- **The curated seed.** `config/related_conditions.yaml` holds the clinical relations, each with a verbatim quote from a PubMed abstract. It was built on 2026-10-16 from 67 PMIDs; all 103 quotes were checked word for word against efetch records.
- **Where it runs.** Every query runs from the launch host, because compute nodes have no internet. Requests are retried with exponential backoff on server errors, HTTP 429 and network failures.

## 2. Which conditions are related, and why (`scripts/resolve.py`)
A condition enters the graph through at least one of three kinds of evidence. Every relation is kept as an evidence row with its source.

1. **Clinical.** The seed lists guidelines, reviews, meta-analyses and cohort or case-control studies whose abstracts state a relation to CTEPH. The categories are:
   - precursor / spectrum;
   - risk factor;
   - comorbidity;
   - intermediate trait;
   - PH group.

   Studies that found **no** association are kept, labelled with their direction. At every build, each quote must occur in the abstract fetched that day (normalised Unicode, dashes and spaces), or the build stops. Evidence from East Asian populations is flagged.
2. **Ontology.** From the Open Targets disease ontology (EFO / MONDO):
   - the parent of CTEPH (pulmonary hypertension);
   - its siblings;
   - the children of pulmonary arterial hypertension.
3. **Shared genetic signal.** A GWAS credible set of the condition colocalises with a CTEPH GWAS credible set: eCAVIAR CLPP ≥ 0.01 or coloc H4 ≥ 0.8. In release 26.09 all colocalisations are eCAVIAR on PICS fine-mapping. A relation is graded by the number of CTEPH loci it shares: **multi-locus** relations are stronger. A trait that shares only the ABO locus mostly reflects ABO pleiotropy.

**Genetics mode.** Each condition gets one:
- `full`: its GWAS and curated gene evidence are fetched;
- `coloc_only`: only the colocalising credible sets are kept (the colocalisation-only traits, and platelet count, whose 125 GWAS are too many);
- `none`: no Open Targets record (procedures such as splenectomy, and CTEPD), or too broad a phenotype (cancer).

## 3. Genetics (`scripts/fetch.py`)
**GWAS.** `studies(diseaseIds)` gives every GWAS study of the condition: its sample size, cases, controls, cohorts, ancestry (`ldPopulationStructure`) and whether summary statistics are available. For each study, every credible set is fetched:
- the lead variant;
- its own P (mantissa and exponent, which avoids float underflow below 10⁻³⁰⁸) and β;
- the fine-mapping method;
- the top three locus-to-gene (L2G) genes.

A study is **East Asian** when its EAS share is ≥ 0.5.

**Curated / Mendelian evidence.** The sources are ClinGen, Gene2Phenotype, Genomics England PanelApp, Orphanet, ClinVar (eva), UniProt and gene burden. For each source:
- the targets it scores ≥ 0.5 for the condition or its ontology descendants are listed (one page ordered by that source's score);
- every evidence row for those targets is then read, with its classification, clinical significance, specific disease and PMIDs. For ClinVar the record's confidence is its review status, so it is kept in its own column (`review_status`) and not mixed with the classifications of the other sources.

The workbook and the knowledge base call this **curated gene evidence**: it is mostly rare-variant, Mendelian evidence, but gene burden and UniProt are not curated gene–disease validity in the ClinGen sense. The tables keep the name `mendelian`.

**CTEPH itself:**
- every associated target;
- its text-mined (Europe PMC) and GWAS evidence;
- the targets' Reactome top-level pathways;
- its GWAS credible sets with all their colocalisations.

**Completeness.** Every list is checked against the count the API reports, and its IDs must be unique. Index/size pages are **not** stable where the sort key ties (2026-10-16: the 615 CTEPH targets paged by 100 differed from one full page in 10 targets). A list longer than its first page is therefore re-read as one page of its full length.

## 4. The genetic studies of CTEPH in the literature (`scripts/pubmed.py`, `extract_genetics.py`)
1. **Search.** PubMed, with the query in `params.PubmedQuery`: CTEPH in title or abstract AND a genetics term. Publication date up to the fetch date.
2. **Pre-screen (deterministic, before any model call).**
   - Records typed Review, Systematic Review, (Practice) Guideline, Editorial, Comment, News or Consensus Development Conference are excluded.
   - So are records whose title and abstract contain no genetic-study term: polymorphism, genotype, allele, mutation, sequencing, exome, genome-wide, GWAS, HLA, SNP, rsID, germline, somatic, clonal haematopoiesis / CHIP, heritability, Mendelian randomisation, copy number, pathogenic / rare variant, missense, haplotype, linkage.
   - On 2026-10-16 this kept 98 of 471 records, including all 16 known CTEPH genetic studies used as a check.
   - Excluded records are listed in `screening.tsv` with the rule that excluded them (`stage = rule`). The user chose this over reading all records (2026-10-16): Open Targets has no literature summaries, and its text mining misses and mis-maps genetic studies.
3. **Screening and extraction (pass A).**
   - An isolated language model reads the remaining abstracts, five per call. It decides whether each record is an original genetic study (or a case report) of CTEPH.
   - For the included ones it extracts:
     - design;
     - population (country, ancestry, Japanese / East Asian);
     - cases and controls, and the control type;
     - genes (as written, plus the HGNC symbol), variants;
     - findings (result, effect, P, and the source sentence);
     - replication;
     - the methods the abstract states.
4. **Methods (pass B).** For every included study, the methods come from the Methods section of its open-access full text (Europe PMC). When no full text can be fetched, they come from its PubMed title and abstract instead (the user, 2026-10-16). Each field records its source. The fields extracted are:
   - cohort and recruitment, diagnosis, controls, DNA source;
   - genotyping platform, sequencing, HLA typing, variant calling and quality control, imputation;
   - statistical test, covariates, multiple testing, replication, power, software.

   Each field comes with its source passage. Europe PMC answers some non-open-access PMCIDs with a persistent HTTP 500; that is treated as "no full text".
5. **Grounding.** Every gene (as written or its symbol), variant, finding quote and sample size must occur in the abstract. Every methods value needs a quote that occurs in the text it came from. "Japanese" or "East Asian" must be stated in the text. Items that fail are dropped and logged (`extract_log.tsv`). Official symbols are kept only when Open Targets knows them by exact name.
6. **Isolation.** Each call is `claude -p --safe-mode --tools "" --strict-mcp-config --no-session-persistence`, with auto-memory disabled, run from `/tmp/cteph_kg_llm`. No project memory, file or instruction reaches the model; `verify.sh` section 5 tests this.
7. **Caching.** Answers are cached by hash of (model, effort, system prompt, prompt) in `cache/llm/`, so later snapshots call the model only for new records.

## 5. Mechanisms and summaries (`scripts/summarise.py`)
- **Mechanism themes.** The text-mined CTEPH targets are grouped by Reactome top-level pathway. For each theme with at least 3 targets, a summary of at most 150 words is written. Its input is up to two sentences from the best-scored publications of each of the theme's 12 best literature targets.
- **Methods and findings.** Summaries of the methods and of the findings of the CTEPH genetic studies (at most 350 words each) are written from the extracted tables and their tallies.
- **Checks.** Every [PMID] a summary cites must be one of its inputs, and every gene it names must occur in its source text. A summary that fails is kept but withheld from the knowledge base and the workbook.

## 6. Tables, graph, changes, knowledge base, workbook
- **`tidy_graph.py`.** Writes the long tables (docs/OUTPUTS.md) and the graph:
  - nodes: conditions, genes, variants, publications;
  - edges: relation evidence; gene–condition genetic evidence (GWAS top L2G ≥ 0.5, Mendelian, literature); colocalisation; genetic study of CTEPH.
- **`diff_snapshots.py`.** Compares stable keys with the newest earlier snapshot.
- **Gene themes.** The findings of the CTEPH genetic studies are grouped by gene theme (`config/gene_themes.yaml`, an editorial grouping reviewed by hand): fibrinogen; ABO / VWF / FVIII / FXI; thrombophilia and haemostasis; HLA and immunity; clonal haematopoiesis and myeloproliferative genes; PAH genes; endothelium, platelets and vascular tone; anticoagulant pharmacogenetics. A gene belongs to the first theme that lists it. A finding without a gene symbol joins the first theme whose genes or words occur in its gene text. The rest are listed under "Other genes" (logged by the build, to be placed in the config) or "No gene named". Per gene, the studies are counted by result: associated (or enriched), not associated, variants identified, other. "Associated in two or more studies" counts studies, not independent samples.
- **Populations.** A literature study's population is its country, else its ancestry, else "not stated". It is Japanese when the extraction says so or it was done in Japan, and East Asian when it is Japanese, flagged East Asian, or done in China, Korea, Taiwan, Hong Kong or Mongolia.
- **`build_knowledge.py`.** Writes the Markdown/JSON knowledge base in the order of the workbook's questions (ten sections), with generated summaries marked.
- **`write_table.py`.** Writes the workbook in the house style, organised by question (docs/OUTPUTS.md):
  - each question sheet opens with its answer in brief (counts computed from the tables, then the generated summary where there is one);
  - then comes a nested table, group › record › source, built with Excel row outlines: groups open, sources closed;
  - source passages sit in cell notes, and links join the same gene or condition across sheets;
  - the appendix holds the flat tables, with filters.

  The workbook is re-opened and checked before publication: sheets, rows, footnote letters, the reason for every "—", outline levels, every link between sheets, and no internal names in cells or notes.
