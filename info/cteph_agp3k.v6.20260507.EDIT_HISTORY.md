# Edit history — `cteph_agp3k.v6.20260507.xlsx`

Each edit to this sheet is backed up before it is applied. `.bak-<date>` holds the
state **just before** that day's edit.

| file | state |
|---|---|
| `…xlsx.bak-20260716` | original (2026-05-07 data), before any edit |
| `…xlsx.bak-20260729` | after edit 1, before edit 2 |
| `…xlsx.bak-20260824` | after edit 2, before edit 3 |
| `…xlsx` (current) | after edit 3 |

All edits touch data only; no rows added or removed (3,655 rows throughout).

---

## Edit 1 — 2026-07-16 · unify `NovaSeq 6000 30x` → `NovaSeq 30x`

`WGS_Platform` column, 46 rows. Two spellings of one platform were merged; the
label `NovaSeq 6000 30x` no longer exists.

- `NovaSeq 6000 30x` (46) + `NovaSeq 30x` (26) → `NovaSeq 30x` (72)
- Nothing else changed (all 36 other columns identical to `.bak-20260716`).

Why: they are the same platform + target depth, split only by how the instrument
model was written down. Carried through, they would have made two meaningless
platform strata. `WGS_Platform` is the pipeline's authoritative platform label, so
the fix is made once here at the source.

---

## Edit 2 — 2026-07-29 · correct DNBSeq-T7 30x `Observed_Depth` read length

`Observed_Depth` column, 331 rows (every `DNBSeq-T7 30x` sample).

**DNBSeq-T7 30x is a 100 bp library, but the upstream FASTQ depth was computed as
if reads were 150 bp**, overstating its depth by exactly 150/100 = 1.50×. Fix:

```
Observed_Depth(T7) := Observed_Depth(T7) × 100/150
```

- T7 `Observed_Depth` mean 27.71× → **18.47×** (≈19×, its true depth).
- Applied to T7 rows only; no other platform touched (all others are ~150 bp and
  their `Observed_Depth` already agrees with the CRAM to within 2%).
- The clean 100/150 is used, not the empirically measured 1.52×: that extra ~2% is
  the ordinary CRAM-vs-FASTQ offset every platform shows, not the read-length error.

Evidence: `check/down_sampling/observed_depth_audit/` (audit + figure + README).
The audit derives the implied read length from `L_assumed = L_real × Observed/CRAM`
and finds ~150 bp for T7 against its true 100 bp, ~150 bp elsewhere.

---

## Edit 3 — 2026-08-24 · rebuild `Cram_Path` from the current CRAM directories

`Cram_Path` column, 3,620 rows changed. It had stopped describing where the CRAMs are:

- **3,135 rows held `/mnt/ngs8_disc1/pub/ngs/WGS/cram/hg38/`** — a bare directory, not a
  file. That mount is gone, and the value was never a usable path to begin with.
- **457 rows** pointed into `/LARGE1/gr10478/workspace/pipeline/output/2020070…/`, the old
  per-sample pipeline output tree.

Rebuilt in two steps. The column was **cleared outright**, every row; then every
`Flag_JHRPv6 = True` row was filled by matching `ID_JHRPv6` exactly against the CRAM
basename in the two directories the CRAMs now live in:

```
/LARGE1/gr10478/pub/WGS/cram/Pulmonary_Hypertension       712 *.cram
/LARGE1/gr10478/pub/WGS/cram/AGP3K                      3,148 *.cram
```

The written value is `<dir>/<ID_JHRPv6>.cram`.

| | |
|---|---|
| `Flag_JHRPv6 = True` rows | 3,592 (`ID_JHRPv6` missing on 0) |
| resolved to exactly one CRAM | **3,592** |
| basename found in *both* directories (ambiguous) | 0 |
| **not found** | **0** |

Split **3,135** `AGP3K` + **457** `Pulmonary_Hypertension`. Every written path exists on
disk and its basename equals the row's `ID_JHRPv6`. The 268 CRAMs in those directories that
no row uses are samples outside JHRPv6.

The 3,620 changed rows are the 3,592 filled plus **28 `Flag_JHRPv6 = False` rows whose stale
value was cleared and deliberately not refilled** — only `Flag_JHRPv6 = True` carries a path
now, and all 63 `False` rows are blank. The other 35 rows were already blank. Nothing else
changed: all 36 other columns are identical to `.bak-20260824`.

`tmp/cram.v6/cram.v6.summary.csv` answers the same question against an older directory
layout and resolves the symlinks to their targets. **It is stale and was not used as a
source here.**

### Known consequence: the index is not beside the path

All 3,135 `AGP3K` entries are **symlinks** into `Control_Nagahama`, `Control_BBJ`, `HTLV1`,
`Control_ACC`, `Lung_cancer` and `cram4temporary`. The `.crai` sits next to the *real* file,
not next to the symlink, so `<Cram_Path>.crai` exists for **457 of 3,592** rows.

That matters to any consumer that derives the index by appending — as
`check/down_sampling/down_sampling.nf:798` does:

```groovy
file("${row[params.cram_path_col]}.crai")
```

That pipeline is not affected today because it reads its own `cram_info` CSV rather than this
workbook. A future consumer reading `Cram_Path` directly needs `samtools -X`, its own index
map, or `.crai` symlinks placed alongside. Recorded here rather than worked around, because
the path convention is the deliberate choice: the two directories above are the authority.

---

## Reproduce any comparison

```python
import pandas as pd
cur = pd.read_excel("cteph_agp3k.v6.20260507.xlsx")
bak = pd.read_excel("cteph_agp3k.v6.20260507.xlsx.bak-20260729")   # or .bak-20260716
cur = cur.reindex(columns=bak.columns)
for col in bak.columns:
    ne = ~((cur[col] == bak[col]) | (cur[col].isna() & bak[col].isna()))
    if ne.sum():
        print(col, int(ne.sum()))
```

## Note on file size

`.bak-20260716` is 460 KB; `.bak-20260729` and `.bak-20260824` are ~288 KB. That
drop is not a data difference — rewriting the workbook drops redundant style
caches.

The current file is **468 KB**, and that jump *is* data. Before edit 3, 3,135 of
the `Cram_Path` cells held the same short string (`/mnt/…/cram/hg38/`), which the
workbook stores once in its shared-string table. They now hold 3,135 distinct
full paths. Nothing else grew.
