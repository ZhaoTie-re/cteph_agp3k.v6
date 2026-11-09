# Refactor proposals — `check/platform_depth` (not applied)

The 2026-10-31 refactor changed only what the cache does not see: the report-stage scripts, the report processes
and block of the workflow, `verify.sh` / `verify_checks.py`, the docs, and new files (`scripts/checks.py`,
`scripts/report_lib.py`, `scripts/frozen.py`, `run.sh`, `docs/PIPELINE.md`, `docs/FROZEN.manifest.tsv`). Every
report output stayed byte-identical (report.4 against report.3). What follows would touch a frozen piece or a
published column; each item names what it would re-run. None is applied.

## What re-runs what

Nextflow's standard cache hashes a task from the session, the process name, its script and every staged file by
path + size + mtime. A new hash re-runs the task; its outputs get a new path and mtime, so every task downstream
re-runs too. Measured cost of the recorded runs (inputs.2, pilot.4–7, full.1, full.3; realtime × allocated cores):
903 core-h, of which HLAHD 639, FIT_NULL 134, JOINT 39, CALL_GVCF 36, DEPTH 28.

| Piece changed | Re-runs | Cost |
|---|---|---|
| `scripts/common.py` (staged by SAMPLE_TABLE, BUILD_TARGETS, DEPTH, … ) | everything | ≈ 900 core-h, 2–3 days of queue |
| `sample_table.py`, `build_targets.py`, MAKE_CHUNKS' inline script | everything (the sample table, targets and chunks feed every task) | ≈ 900 core-h |
| `arm_plan.py`, `--Replicates` (ARM_PLAN interpolates it) | every arm: reads, calls, HLA-HD, re-tests | ≈ 740 core-h |
| `make_subsets.py`, `control_subsets.py` | null fits, step 2, gene tests, omnibus, statistics of the subsets | ≈ 140 core-h |
| `platform_stats.py`; CARRIER_FIRTH's script, or a param it interpolates (`params.Strata`, `params.CasePlatformSlugs`) | PLATFORM_STATS ×3; with CARRIER_FIRTH ×3 | < 1 core-h, but outside the report stage |
| GATE_SUMMARY, GATE_G5 (terminal, local) | themselves | seconds |
| report processes, their scripts, `report_lib.py`, `checks.py`, `catalogue.py` | nothing (`cache false`) | minutes |

`run.sh` refuses to launch when a frozen piece changed (`scripts/frozen.py check`), and its guard stops the driver
when a task that must be cached is submitted.

## A. Code and workflow (would touch frozen pieces)

1. **Split `platform_depth.nf` into modules** (`modules/{inputs,gates,arms,platform,report}.nf`, included by the
   main file). The process names stay (no named sub-workflows: a process called inside one gets `WF:NAME` as its
   name, and the name is hashed), the script texts stay, so the hashes should not change — but it is not verified,
   and a miss on SAMPLE_TABLE re-runs everything. Test protocol before applying: a scratch launch directory with a
   copy of `.nextflow/cache/8eec4a16-…` (the same session UUID, so the same hashes), `-w` pointing at this `work/`,
   `--out_dir <scratch>` (no cached script interpolates `params.out_dir`), `--run inputs` behind the guard: a miss
   then writes only to the copy. Cost if the test passes: none.
2. **One home for the names.** `CONTROL_SOURCES` (in `control_subsets.py`), the variant sets, the platform slugs
   and labels now live in `common.py`, `control_subsets.py`, `report_lib.py` and the workflow's `params.Cohorts`,
   `params.Strata`, `params.CasePlatformSlugs`, `params.ControlSources`. Moving them into `common.py` and deriving
   the params from it re-runs everything (≈ 900 core-h); fold it into the next full re-run if there is one. Until
   then `verify.sh` §17 holds the copies equal.
3. **`stratum` → `variant_set`** for the genes' variant sets (`params.Strata`, `gene_retest.py`'s output column,
   `carrier_firth.py`): the same word means a case-platform stratum elsewhere. Re-runs RETEST_GENE, PLATFORM_GENE,
   CARRIER_FIRTH, PLATFORM_STATS (≈ 1 core-h) but also changes published columns (B.2).
4. **Inline scripts into `scripts/`**: GATE_SUMMARY and GATE_G5 (cheap, terminal); never MAKE_CHUNKS (its chunk
   files feed DEPTH and the fractions, i.e. every arm).
5. **A G7 gate record.** G7 (the swap leaves every other genome byte-identical) is enforced by `bedio.py swap`,
   which aborts otherwise, and is now checked by `verify.sh` from the 66 `swap.*.json`; a GATE_G7 process writing
   `_gates/<label>/G7.json`, and GATE_SUMMARY listing G5–G9 beside G1–G4, would make the gate files complete.
   A new process re-runs nothing, but it edits the workflow outside the report block (re-freeze the manifest with
   `python3 scripts/frozen.py write --force` after a reviewed diff).
6. **`saige_collect.py` column order** (`cohort` in the middle of `platform.strata.tsv`): frozen; re-runs every
   step 2 (≈ 10 core-h) and changes published columns.

## B. Published columns (would change published content; `deck/` and the user read them)

The conventions that already hold everywhere are in [OUTPUTS.md](OUTPUTS.md#table-conventions) and are checked by
`verify.sh` §15. These deviations would need a column change:

1. **z observed:** `z_obs` (`10.platform/<cohort>/platform.lopo.tsv`, `12.controls/controls.loo.tsv`,
   `controls.signals.tsv` `weakest_loo_z_obs`) vs `z_observed` (`14.evidence/evidence.loo.tsv`).
2. **`stratum`:** the genes' variant set in every gene table (`ds.genes`, `platform.genes`, `controls.genes`,
   `followup.gene_*`), a case-platform stratum in `platform.rho.tsv` (`stratum_1`, `stratum_2`).
3. **P column:** `pvalue` in the rvtest-derived gene tables (`ds.genes` `pvalue_<arm>`, `platform.genes`,
   `controls.genes`) vs `p` everywhere else.
4. **Confidence intervals:** `platform.signals.tsv` writes `or_only_<slug>_ci` as text `l-u`; every other table has
   numeric `*_l95` / `*_u95`.
5. **Sample counts:** `cases` / `controls` (`platform.null_fits`), `cases_random` / `controls_random`
   (`controls.counts`), `n_case` / `n_ctrl` (strata), `n_cases` / `n_controls` (`array.evidence`), `n`,
   `samples` (`followup.gene_filters`).
6. **Control frequency:** `af_ctrl` (strata), `control_freq_wgs` (array),
   `freq_control` (re-test tables).
7. **Identifier columns:** a signal's marker is `marker` (`platform.signals`, `controls.signals`, `array.*`,
   `platform.case_only`), `variant_id` (heterogeneity, LOPO, strata), `unit` (`ds.positions`, `followup.*`). One
   pair `signal` (display name) + `marker` (the tested ID) in every per-signal table would join them directly.
8. **`cohort`:** first in `10.platform/_summary` and `12.controls`; after the keys in `ds.*` and `14.evidence`; in
   the middle of `platform.strata.tsv`; absent from the other `10.platform/<cohort>/` tables (the folder names it).
9. **Arm columns:** `ds.hla_n3.tsv` names them bare (`baseline`, `original`, `r01` …) where every other `ds.*` table
   writes `<quantity>_<arm>`; the pivots order them alphabetically (`baseline` before `original`).
10. **`controls.counts.tsv`:** `source` last; the pooled row says `all`.
11. **`evidence.long.tsv` `value`:** its meaning depends on the check (Δ|z|, P, OR, a ratio); a `value_kind`
    column, or the registry's readout per check, would say which.
12. **Figure records:** `03.platform/forest.<cohort>.stats.json` and `05.controls/sources.<cohort>.stats.json` hold
    only the number of markers (identical across cohorts), not the values drawn, which every other `stats.json`
    holds (catalogue.py's promise).

A change of any of these is a report-stage change only (minutes, no cache), but it changes what the deck and the
user read; it needs the user's decision and a synchronized update of `docs/OUTPUTS.md`, `figures.py` and
`evidence_matrix.py` (which read some of these columns).
