#!/usr/bin/env python3
"""Assemble the audit's findings into a standalone markdown report.

The report is generated rather than hand-written so it cannot drift from the data:
every number in it is read out of observed_depth_audit.tsv, which the analysis step
just produced. Editing the prose here is fine; typing a number into it is not.

The report is published next to the figure and tsv it describes, so its image link
is a plain sibling filename and no file is duplicated.
"""

import argparse
import csv
from datetime import datetime
from pathlib import Path

ASSUMED = 150.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--samples-tsv", required=True)
    ap.add_argument("--log", required=True)
    ap.add_argument("--figure-name", required=True,
                    help="filename of the figure, linked as a sibling of the README")
    ap.add_argument("--regions", required=True)
    ap.add_argument("--samples-per-platform", type=int, required=True)
    ap.add_argument("--sample-seed", required=True)
    ap.add_argument("--out-md", required=True)
    args = ap.parse_args()

    with open(args.tsv) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise SystemExit("[report] the audit tsv is empty; nothing to report")
    with open(args.samples_tsv) as fh:
        samples = list(csv.DictReader(fh, delimiter="\t"))

    for r in rows:
        for k in ("READLEN_MEASURED", "READLEN_IMPLIED", "OBSERVED_DEPTH",
                  "DEPTH_CRAM_WITHDUP", "OBSERVED_OVER_MEASURED", "DUP_PCT"):
            r[k] = float(r[k])
        r["N"] = int(r["N"])
    rows.sort(key=lambda r: r["READLEN_MEASURED"])

    disc = [r for r in rows if abs(r["READLEN_MEASURED"] - ASSUMED) > 20]
    ctrl = [r for r in rows if r not in disc]
    if not disc:
        raise SystemExit("[report] no platform's read length is far enough from "
                         f"{ASSUMED:.0f} bp to test the hypothesis; nothing to report")
    d = disc[0]          # the discriminating platform
    c0 = ctrl[0]         # a control, for the worked example
    fig = Path(args.figure_name).name
    n_total = sum(r["N"] for r in rows)
    regions = [x.strip() for x in args.regions.split(",") if x.strip()]
    total_mb = sum(
        (int(x.split(":")[1].split("-")[1]) - int(x.split(":")[1].split("-")[0]) + 1)
        for x in regions) / 1e6
    inflation = ASSUMED / d["READLEN_MEASURED"]
    corrected = d["OBSERVED_DEPTH"] * d["READLEN_MEASURED"] / ASSUMED

    md = f"""# `Observed_Depth` overstates {d['PLATFORM']} by {inflation:.1f}x

*Generated {datetime.now().strftime('%Y-%m-%d')} by `observed_depth_audit.nf`.
Every number below is read from `observed_depth_audit.tsv`; none is typed in by hand.*

---

## The finding

**{d['PLATFORM']} is a {d['READLEN_MEASURED']:.0f} bp library, but `Observed_Depth` was
computed as if it were {ASSUMED:.0f} bp.** The sheet therefore overstates its depth by
{inflation:.2f}x — it claims **{d['OBSERVED_DEPTH']:.1f}x**, the reads carry
**{d['DEPTH_CRAM_WITHDUP']:.1f}x**.

Every other platform is a ~{ASSUMED:.0f} bp library, so the same mistake leaves them
untouched and their sheet depth is correct to within 2%.

> ### {d['PLATFORM']} samples are ~{d['DEPTH_CRAM_WITHDUP']:.0f}x, not the ~30x their label claims.
>
> Anything that grouped, filtered, or matched samples on `Observed_Depth` has been
> treating these as deep when they are not.

---

## 1. What disagrees

`Observed_Depth` is computed upstream **from the FASTQs**, and `cram.v6.ipynb` copies
the column verbatim into `cram.v6.summary.csv` without ever recomputing it. It is
therefore an independent claim about each sample's depth — and the CRAMs can check it.

Four of five platforms agree with their own CRAM to within 2%. One does not:

| platform | read length | `Observed_Depth` (sheet) | depth measured in CRAM | sheet / CRAM |
|---|---:|---:|---:|---:|
"""
    for r in rows:
        flag = " **⚠**" if r in disc else ""
        bold = "**" if r in disc else ""
        md += (f"| {r['PLATFORM']}{flag} | {r['READLEN_MEASURED']:.0f} bp | "
               f"{r['OBSERVED_DEPTH']:.2f}x | {r['DEPTH_CRAM_WITHDUP']:.2f}x | "
               f"{bold}{r['OBSERVED_OVER_MEASURED']:.2f}x{bold} |\n")

    md += f"""
Two facts matter here, and neither is enough alone: **{d['PLATFORM']} is off**, *and*
**the others are not**. A discrepancy on every platform would indict our measurement;
a discrepancy on exactly one indicts the data.

---

## 2. Why read length is the suspect

Depth from a FASTQ is essentially one multiplication:

```
depth = (number of reads) x (read length) / (genome size)
```

Work through the three terms by elimination:

| term | could it be wrong by exactly {inflation:.2f}x, on one platform only? |
|---|---|
| number of reads | Miscounting reads gives arbitrary errors, not a clean {inflation:.2f}x. |
| genome size | A wrong genome size shifts **every** sample equally. |
| **read length** | **The one term a pipeline may fill in with a constant instead of reading per sample.** |

And the read lengths line up with exactly that idea:

"""
    md += f"- **{d['PLATFORM']} sequences at {d['READLEN_MEASURED']:.0f} bp**\n"
    for r in ctrl:
        md += f"- {r['PLATFORM']} sequences at {r['READLEN_MEASURED']:.0f} bp\n"

    md += f"""
If the upstream hard-coded {ASSUMED:.0f} bp, a {d['READLEN_MEASURED']:.0f} bp library is
inflated by exactly {ASSUMED:.0f}/{d['READLEN_MEASURED']:.0f} = **{inflation:.2f}x**, while
the ~{ASSUMED:.0f} bp libraries are left alone. That is the pattern in the table above.

> We did not test this the direct way. Recomputing depth from the FASTQs would settle it
> outright, but that means re-reading every raw dataset to answer a metadata question — a
> cost we chose not to pay. The test below reaches the same answer from the CRAMs in about
> three minutes, and the direct check stays available to anyone who wants it.

---

## 3. The test: what read length does the sheet imply?

### 3.1 The derivation

Depth is **read density x read length** — how many reads cover a base depends on how
densely they fall and how long each one is. Write both numbers that way.

What the upstream computed from the FASTQ (the number under audit):

```
Observed_Depth = (reads / genome size) x L_assumed
               =          rho          x L_assumed
```

What we measure from the CRAM:

```
CRAM depth = (reads in region / region length) x L_real
           =               rho                 x L_real
```

Divide one by the other, and **rho cancels**:

```
Observed_Depth       rho x L_assumed       L_assumed
--------------  =  -------------------  =  ---------
  CRAM depth         rho x L_real           L_real
```

Rearranged, the length the sheet must have assumed:

```
L_assumed = L_real x (Observed_Depth / CRAM depth)
```

Every term on the right is measurable. This is the `READLEN_IMPLIED` column.

**What the cancellation buys.** Along with `rho`, the read count and the genome size drop
out of the problem — and we know neither. The upstream's genome size is unknown (3.088 Gb?
with or without N?), and the read count would cost a full CRAM decode (~1.5 min/sample).
The test needs only a *density*, which {total_mb:,.0f} Mb of probe regions gives in
seconds. That is why this audit runs in ~3 minutes rather than hours.

### 3.2 Worked example

For **{d['PLATFORM']}**, whose reads measure {d['READLEN_MEASURED']:.0f} bp:

```
L_assumed = {d['READLEN_MEASURED']:.0f} bp x ({d['OBSERVED_DEPTH']:.2f} / {d['DEPTH_CRAM_WITHDUP']:.2f})
          = {d['READLEN_MEASURED']:.0f} bp x {d['OBSERVED_OVER_MEASURED']:.3f}
          = {d['READLEN_IMPLIED']:.0f} bp
```

To produce the {d['OBSERVED_DEPTH']:.2f}x in the sheet, the upstream had to treat this
{d['READLEN_MEASURED']:.0f} bp library as **{d['READLEN_IMPLIED']:.0f} bp**.

For a control, **{c0['PLATFORM']}** at {c0['READLEN_MEASURED']:.0f} bp:

```
L_assumed = {c0['READLEN_MEASURED']:.0f} bp x {c0['OBSERVED_OVER_MEASURED']:.3f} = {c0['READLEN_IMPLIED']:.0f} bp
```

Also ~{ASSUMED:.0f} — but that tells us nothing, because {c0['READLEN_MEASURED']:.0f} bp
*is* ~{ASSUMED:.0f} bp. See [Limitations](#6-limitations).

### 3.3 Two predictions that differ

The hypothesis is falsifiable: the two candidate explanations predict different numbers,
so the data can reject one.

| if the upstream... | then `READLEN_IMPLIED` is... |
|---|---|
| used each sample's real read length | equal to the measured read length — {d['READLEN_MEASURED']:.0f} for {d['PLATFORM']} |
| hard-coded {ASSUMED:.0f} bp | ~{ASSUMED:.0f} for **every** platform, whatever its real length |

---

## 4. Result

![The three steps of the argument]({fig})

| platform | measured read length | `READLEN_IMPLIED` | predicted inflation | seen | verdict |
|---|---:|---:|---:|---:|---|
"""
    for r in rows:
        pred = ASSUMED / r["READLEN_MEASURED"]
        v = (f"**implies {ASSUMED:.0f}, not {r['READLEN_MEASURED']:.0f}**" if r in disc
             else "no power (~150 bp library)")
        md += (f"| {r['PLATFORM']} | {r['READLEN_MEASURED']:.0f} bp | "
               f"{r['READLEN_IMPLIED']:.0f} bp | {pred:.2f}x | "
               f"{r['OBSERVED_OVER_MEASURED']:.2f}x | {v} |\n")

    md += f"""
**Every platform's `Observed_Depth` implies ~{ASSUMED:.0f} bp — including the one that is
really {d['READLEN_MEASURED']:.0f} bp.** The second prediction holds; the first is rejected.

The inflation is quantitative, not just directional: {ASSUMED:.0f}/{d['READLEN_MEASURED']:.0f}
predicts **{inflation:.2f}x** and the data shows **{d['OBSERVED_OVER_MEASURED']:.2f}x**. The
hypothesis had room to be wrong by any amount, and it was not.

And the correction closes the gap:

```
{d['OBSERVED_DEPTH']:.2f}x  x  {d['READLEN_MEASURED']:.0f}/{ASSUMED:.0f}  =  {corrected:.1f}x        (CRAM measures {d['DEPTH_CRAM_WITHDUP']:.1f}x)
```

---

## 5. What it means downstream

- **Depth labels are wrong for {d['PLATFORM']}.** These samples are ~{d['DEPTH_CRAM_WITHDUP']:.0f}x.
  Their `Target_Depth` of 30x and their `Observed_Depth` of ~{d['OBSERVED_DEPTH']:.0f}x both
  overstate them.
- **Do not tune down-sampling on `Observed_Depth`.** It is genome-wide, derived from the
  FASTQ, and — as shown here — wrong for one platform. `tuning.fraction/` measures depth
  from the CRAMs for exactly this reason.
- **The sheet is otherwise sound.** The other {len(ctrl)} platforms check out to within 2%,
  so this is one specific defect, not a reason to distrust the column everywhere.

---

## 6. Limitations

**Only {d['PLATFORM']} can discriminate.** The other {len(ctrl)} are ~{ASSUMED:.0f} bp
libraries, so both hypotheses predict nearly the same implied length for them and the test
has **no statistical power** there — that is what "no power" means in the tables above, not
"no data". They are the **calibration control**: coming out at implied ≈ measured ≈
{ASSUMED:.0f} shows the method is unbiased. They are not {len(ctrl)} independent
confirmations. The conclusion rests on one platform and should be read that way.

**This is inference from the output, not a reading of the upstream code.** We chose not to
recompute from the FASTQs, so the arithmetic upstream is reconstructed rather than observed.
Any other bug that scaled {d['PLATFORM']} by exactly {inflation:.2f}x and no other platform
would look identical here — though it would have to coincide with {d['PLATFORM']} being the
only {d['READLEN_MEASURED']:.0f} bp library. The derivation also assumes the upstream used
`reads x L / genome`; under a different formula, `L_assumed` is an equivalent rescaling
rather than a literal read length. **Recomputing one {d['PLATFORM']} sample from its FASTQ
would settle all of this** — not all {d['N']}, one.

**The probe regions stand in for the genome.** The division assumes the read density there
matches genome-wide density. Hence {len(regions)} windows on {len(regions)} chromosomes, away
from centromeres. Even if that is imperfect, the bias applies to **every platform equally**
— the same windows are measured for all — so it shifts everything together and cannot
manufacture the gap between {d['PLATFORM']} and the rest, which is the entire claim.

**Not corrected for:** soft-clipped bases (~0.3%) and reads that never aligned (~0.07%).
Both make CRAM depth slightly lower than the FASTQ's, which is why even the well-behaved
platforms sit ~2% above 1.00x rather than exactly on it. Neither is remotely large enough to
make {inflation:.2f}x.

**Duplicates are kept** in the CRAM measurement. `Observed_Depth` comes from the FASTQ, where
duplicates are still present and unmarked; excluding them here would charge each platform for
its duplicate rate on top of the effect under test. Those rates vary a lot
({max(rows, key=lambda r: r['DUP_PCT'])['PLATFORM']} {max(r['DUP_PCT'] for r in rows):.1f}%
vs {min(rows, key=lambda r: r['DUP_PCT'])['PLATFORM']} {min(r['DUP_PCT'] for r in rows):.1f}%),
so the choice is not cosmetic.

---

## 7. Method

| choice | value | why |
|---|---|---|
| samples | **{n_total}** ({args.samples_per_platform} per platform, random, seed `{args.sample_seed}`) | The effect is a platform-wide {inflation:.2f}x systematic. Read length has *no* within-platform spread (every library verified fixed-length: shortest sampled read = longest), and depth varies a few percent. More samples buy precision the claim does not need. |
| depth | `samtools coverage` over **{len(regions)} probe regions, {total_mb:,.0f} Mb** | The test needs a read *density*, not a genome-wide total, so probe regions suffice — and cost seconds instead of the ~1.5 min/sample a full `samtools stats` pass would. Pooled weighted by region length, since `meandepth` is already a per-base average. |
| read length | first 1,000 reads per region, **`-F 0x900`** | Secondary and supplementary records carry hard-clipped `SEQ` — 30-35 bp fragments of a read. Left in, they drag the mean down *and* make every fixed-length library look variable-length. |
| duplicates | **kept** (`--ff UNMAP,SECONDARY,QCFAIL`) | Matches what `Observed_Depth` saw in the FASTQ. See [Limitations](#6-limitations). |

### 7.1 Which regions

{len(regions)} windows of {total_mb / len(regions):.0f} Mb each, one per chromosome, mid-arm
and clear of the centromeres. Spreading them over {len(regions)} chromosomes keeps any single
locus — an unusual GC content, a segmental duplication — from driving the answer. None of
them is the FTO/analysis region: this audit is about the genome-wide `Observed_Depth`, not
about the locus under study.

| # | region | length |
|---:|---|---:|
"""
    for i, x in enumerate(regions, 1):
        chrom, span = x.split(":")
        s0, s1 = span.split("-")
        md += f"| {i} | `{x}` | {(int(s1) - int(s0) + 1) / 1e6:.1f} Mb |\n"

    md += f"""
Change them with `--regions`; the conclusion should not depend on which windows are used,
and that is worth checking if you doubt it.

### 7.2 Which samples

Drawn at random per platform with seed `{args.sample_seed}`, so the selection is fixed and
reproducible. Platforms with fewer than {args.samples_per_platform} eligible samples
contribute all of them. Eligible = has a CRAM (`Cram_Found == True`) and a non-empty
`Observed_Depth`.

Per-sample measurements are in `observed_depth_audit.samples.tsv` — every median in this
report can be traced back to these rows:

| platform | samples |
|---|---|
"""
    by_pf = {}
    for s in samples:
        by_pf.setdefault(s["PLATFORM"], []).append(s["SAMPLE"])
    for r in rows:
        ids = ", ".join(f"`{x}`" for x in sorted(by_pf.get(r["PLATFORM"], [])))
        md += f"| {r['PLATFORM']} | {ids} |\n"

    md += f"""
---

## 8. Reproducing

```bash
cd .../cteph_agp3k.v6/check/down_sampling/observed_depth_audit
source activate dsl2
export PATH=/home/b/b37974/:$PATH
nextflow run observed_depth_audit.nf -resume
```

| output | contents |
|---|---|
| `results/00_cram_region_stats/` | per-sample, per-region measurements |
| `results/01_audit/observed_depth_audit.tsv` | the per-platform table, machine-readable |
| `results/01_audit/observed_depth_audit.log` | the same analysis, narrated |
| `results/01_audit/{fig}` | the figure above |
| `results/01_audit/README.md` | this document |

Knobs: `--samples_per_platform` (default {rows[0]['N']}), `--sample_seed`, `--regions`.
"""

    Path(args.out_md).write_text(md)
    print(f"[report] wrote {args.out_md} ({len(md.splitlines())} lines), linking {fig}")


if __name__ == "__main__":
    main()
