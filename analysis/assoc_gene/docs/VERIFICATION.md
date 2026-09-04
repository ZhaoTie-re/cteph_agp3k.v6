# Verification — `assoc_gene`

Two layers. **In-task gates** abort the pipeline at the point where a silent failure would otherwise be
born. **`verify.sh`** is read-only and runs before a launch (sections 1–8, 12, 15, 2b) and again after
(9–11, 13–14, 16–21). A section that cannot run yet says `skipped`; nothing is reported as passing that
was not checked.

## In-task gates

| where | what aborts |
|---|---|
| `PREP_PHENO_COV` | a genotyped sample without a phenotype; an empty class after normalisation; no `SEX` column |
| `MERGE_MAP` (B1) | a set on two chromosomes without the `@chr<N>` suffix; a set below `MinNumVar` surviving |
| `EMIT_TEST_INPUTS` | `set(map) ≠ set(bim)`; bim A1 ≠ the ID's ALT; a chromosome the VCF cannot name; a set name rvtest would mis-tokenise; keep lists not partitioning the `.fam` or disagreeing with `pheno_coding.json`; **bed-decoded ALT / observed counts ≠ `.acount` for any variant and group**; `carriers ≤ MAC ≤ 2·n_var·carriers` or `multi ≤ carriers ≤ N` violated |
| `RVTEST` | a named covariate absent from `covar_rvt.tsv`; rvtest's "Loaded *n* cases" ≠ `pheno_coding.json`; any missing phenotype; the `--impute` strategy not confirmed in the log; the declared output file absent |
| `MERGE_RVTEST` (B2) | a `Gene` first column (refFlat, not the map); a method tail that is not the declared method's; a wald part empty for cmc or non-empty for skato; cmcWald block height ≠ `1 + n_covar`; a set returned by the score model and not the Wald model; a returned set with no `set_stats` row |
| `SCAN` (B3) | a merged table off the `gene_scan` schema; a denominator for another cell or disagreeing with the gene index; any set `missing` or `unexplained` in the map identity (unless `--no-strict`); a Bonferroni call that is not a BH call |
| `SIGNIFICANCE` | `n_genes_mapped`, `min_num_var` or `alpha` differing between the denominator and `scan_qc`; a threshold in either file ≠ α / n |
| `ROBUSTNESS` | a set naming different genes in different cohorts; an evidence row count ≠ the design; `alpha` in `denominators.tsv` ≠ `--alpha` |

## `verify.sh`

### Pre-run

| § | check | what it proves |
|---|---|---|
| 1 | `nextflow lint` clean; every script `py_compile`s | syntax |
| 2 | no `cpus`/`memory`, no `task.cpus`, no `--geneFile`; every `m=` a multiple of 4,571; `errorStrategy = 'finish'` | site rules; a figure failure cannot kill 264 rvtest tasks |
| 2b | a `publishDir` that references a process input is a closure | the "No such variable" runtime failure lint cannot see |
| 3 | every flag each process passes is accepted by its script's `--help`; every process running a script is in the map (15 pairs incl. `CATALOGUE` and `SCAN_SCALE`) | `.nf` ↔ script contract |
| 3b | every `stageAs`-expanded input goes to an `nargs='+'` flag | name-right, arity-wrong failures |
| 4 | `../_shared` and the four sibling components unmodified since the session start | isolation |
| 5 | every `plot_*.py` imports `plot_style` / `figure_doc` | the shared style loads |
| 6 | a 20-second real rvtest run on 60 callset variants: output suffixes match `params.RvtestMethods`; cmcWald writes `1 + n_covar` rows per gene; SkatO's header has `rho` at column 7 | the assumptions the shell string makes about rvtest |
| 7 | the six documents exist (incl. `docs/FIGURES.md`) | |
| 8 | `SCAN_COLUMNS` identical in `merge_rvtest.py`, `scan.py`, `robust_genes.py`; `QC_COLUMNS` in `scan.py`, `gene_significance.py`; `DENOM_COLUMNS`; `STATES` in `robust_genes.py` ↔ `vocab.py`; `METHODS` in five scripts incl. `vocab.py`; every reader's `set_stats` / `scan_qc` columns ⊆ the writer's | schema atomicity across scripts |
| 12 | no occurrence of vocabulary from any retired design (the pattern is in `verify.sh` §12) in the `.nf`, config, scripts or documents | the component describes one design |
| 15 | `robust_genes.tier_of` on seven fixture cases | the tier rule |

### Post-run

| § | check | what it proves |
|---|---|---|
| 9 | `trace.txt` has exactly the expected task count per process (397 in total, incl. `CATALOGUE` ×2, `GENE_MODELS` and `SCAN_SCALE`), no unknown process, every task completed | the DAG ran whole; a barrier that emitted nothing is visible |
| 10 | for every `03.assoc/<s>/cmc/merged.tsv`, `n_var_map + span_overlap_excess == n_var_engine` on every gene | the tested sets are the mapped sets, exactly |
| 11 | one Bonferroni value per (cohort, stratum) across the two methods | one denominator judged both |
| 13 | every `set_stats.tsv` row: `mac == case + control`; `carriers ≤ MAC ≤ 2·n_var·carriers`; `multi ≤ carriers ≤ N`; N equals `pheno_coding.json` | the counts are self-consistent and on the right samples |
| 14 | every `gene_scan_all.tsv` row: cmc `or == exp(beta)`, `or_l95 ≤ or ≤ or_u95`, `rho` NA; skato no `or`/`beta`, `rho` present wherever p is | the effect sizes are what METHODS says |
| 16 | every cell of `summary_table.<ReportCohort>.tsv` matches its column's format regex | the table is the manuscript's layout |
| 17 | `evidence.tsv` rows == Σ<sub>s</sub> |genes<sub>s</sub>| × cohorts × 2; `not_in_map` per cell == 2 × (|union| − |index|) | the evidence table is the full design |
| 18 | `map_identity_all.tsv` reasons ⊆ {ok, explained_by_span_overlap, explained_by_collision} | no unexplained set anywhere |
| 19 | every family directory under `results/figures/` has a `README.md` and no other Markdown; a multi-member family's README names every member; `05.genes/index.tsv` rows == PNGs and every named PNG exists; no PNG outside `results/figures/` | one README per figure family, every figure documented once |
| 20 | every `variant_stats.tsv` `variant_id` matches `^chr<N>:<POS>:<REF>:<ALT>$` and agrees with its own `pos`/`ref`/`alt` | the ID invariant the decode and the predictor parse |
| 21 | every PNG under `results/figures/` is 7.20 in wide (4,320 px at 600 dpi) with a plausible height (`robust_genes.png` computed, ≤ 11.5 in + caption) | every figure is authored at the journal width |

### Not covered

Statistical validity of the tests themselves (rvtest's), the platform confound, and the biological
plausibility of any gene. `verify.sh` proves the pipeline computed what METHODS describes, not that
what METHODS describes answers the question.
