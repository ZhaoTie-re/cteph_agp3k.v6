# cteph_knowledge_graph — CTEPH, related conditions, and their genetics

A knowledge graph of chronic thromboembolic pulmonary hypertension (CTEPH, `EFO_0022954`). It is rebuilt as dated snapshots that can be repeated, each compared with the one before. External knowledge only: Open Targets, PubMed and Europe PMC; nothing of this project's own results.

**What it holds**
1. **Related conditions, each with its evidence:**
   - clinical relations from guidelines and studies, every quote checked against its PubMed abstract (`config/related_conditions.yaml`);
   - the disease ontology;
   - GWAS signals that colocalise with a CTEPH locus.

   Categories: precursor / spectrum (PE, VTE, DVT, CTEPD), risk factors, comorbidities, intermediate traits (FVIII, VWF, ADAMTS13, fibrinogen, FXI, D-dimer, platelets, Lp(a)), the PH group, and shared-signal traits.
2. **Genetics of CTEPH and of each condition:**
   - GWAS studies and credible sets (lead variant, L2G genes, ancestry);
   - curated gene evidence, mostly Mendelian (ClinGen, PanelApp, Orphanet, ClinVar, UniProt, Gene2Phenotype, gene burden);
   - the genetic studies of CTEPH in PubMed, with design, population, findings and **methods** (abstract and open-access Methods).
3. **Mechanisms.** The text-mined CTEPH targets by Reactome theme, with source sentences and grounded summaries.
4. **East Asian / Japanese evidence**, flagged throughout.

**Outputs**, per snapshot (`results/snapshots/<date>/`):
- knowledge base `06.knowledge/cteph_knowledge_graph.md` / `.json`;
- graph `03.tables/nodes.tsv`, `edges.tsv`, `graph.json`;
- workbook `cteph_knowledge_graph_<date>.xlsx`, organised by question: a guide, seven question sheets whose rows open from group to record to source (CTEPH literature by gene theme), the changes, and six flat appendix tables with filters;
- long tables;
- changes since the previous snapshot.

`results/latest.json` points to the newest snapshot. See docs/OUTPUTS.md.

## First snapshot (2026-10-16; Open Targets 26.09)

**Related conditions.** 95 in all:
- 27 by clinical evidence: 76 PubMed-verified quotes, contrary findings included;
- 11 by the ontology;
- 67 by colocalisation with the three CTEPH loci (ABO, F11, FGA/FGG). VTE, DVT, PE, FXI and VWF share two or more loci.

**Genetics of the related conditions.**
- 514 GWAS studies, 53 of them East Asian, with 8,422 credible sets.
- 600 curated gene links. For example: PAH *BMPR2*, *SOX17*; thrombophilia *PROC*, *PROS1*, *SERPINC1*, *F5*, *F2*; MPN *JAK2*.

**CTEPH literature.**
- 471 PubMed records, pre-screened to 98.
- **51 genetic studies** with 139 checked findings. The main ones:
  - fibrinogen *FGA* Thr312Ala, replicated;
  - Japanese HLA-B\*52:01 / DPB1\*02:02;
  - JAK2 V617F / clonal haematopoiesis;
  - *STAB2* rare variants;
  - Chinese and European GWAS.
- Methods: 266 fields, from the open-access Methods section for 19 studies and from the abstract for 32.

**Mechanisms.** 27 Reactome themes, each with a grounded summary.

**Graph.** 1,664 nodes and 2,658 edges.

## Run (from the login node: compute nodes have no internet)

```bash
cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/integration/cteph_knowledge_graph
PATH=/home/b/b37974/anaconda3/envs/dsl2/bin:$PATH ./verify.sh              # pre-flight (1-5; NO_LLM=1 skips the model probe)
PATH=/home/b/b37974/anaconda3/envs/dsl2/bin:$PATH nextflow run cteph_knowledge_graph.nf --FetchDate $(date +%F) -ansi-log false
PATH=/home/b/b37974/anaconda3/envs/dsl2/bin:$PATH ./verify.sh              # post-run (6-9, newest snapshot)
```

**Running.**
- `--FetchDate` names the snapshot. A new date re-queries Open Targets and PubMed.
- Language-model answers are cached in `cache/llm/`, so a later snapshot only screens new PubMed records.
- For a smoke run, add `--LlmLimit 10 --MaxThemes 2 --results_root <scratch>`.
- Give `-resume` its session ID (`nextflow log`).
- `nextflow lint` is the safe syntax check; `nextflow run --help` would launch the pipeline.

**Updating the evidence.** Edit `config/related_conditions.yaml`. Each new evidence item needs a PMID, a verbatim quote from its title or abstract, and a one-line statement. The build stops if a quote is not in the abstract.

**Gene themes.** `config/gene_themes.yaml` groups the CTEPH literature by gene theme. A new snapshot logs the genes no theme lists ("Other genes"); add them to a theme and rebuild.

## Files
- `cteph_knowledge_graph.nf`, `nextflow.config`, `verify.sh`.
- `config/related_conditions.yaml`: the curated seed.
- `config/gene_themes.yaml`: the gene themes of the CTEPH literature.
- `scripts/`:
  - `fetch.py`: Open Targets;
  - `pubmed.py`: PubMed and Europe PMC;
  - `resolve.py`: conditions and evidence;
  - `extract_genetics.py`: literature screening and extraction;
  - `llm.py`: isolated, cached model calls and grounding;
  - `tidy_graph.py`: tables and graph;
  - `summarise.py`;
  - `diff_snapshots.py`;
  - `build_knowledge.py`;
  - `write_table.py`;
  - `manifest.py`;
  - `ot_common.py`.
- `docs/METHODS.md`, `docs/OUTPUTS.md`, `docs/STUDY_NOTES.md`.
