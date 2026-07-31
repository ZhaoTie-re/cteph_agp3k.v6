# Edit history — `cteph_agp3k.v6.20260507.xlsx`

Each edit to this sheet is backed up before it is applied. `.bak-<date>` holds the
state **just before** that day's edit.

| file | state |
|---|---|
| `…xlsx.bak-20260716` | original (2026-05-07 data), before any edit |
| `…xlsx.bak-20260729` | after edit 1, before edit 2 |
| `…xlsx` (current) | after edit 2 |

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

`.bak-20260716` is 460 KB; the later files are ~288 KB. That is not a data
difference — rewriting the workbook drops redundant style caches. The only
data-level differences are the cells listed above.
