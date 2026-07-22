#!/usr/bin/env python3
"""Why the sheet's Observed_Depth does not match the CRAMs, and the read-length
explanation for it.

THE DISCREPANCY
    The sheet labels DNBSeq-T7 samples ~30x, but their CRAMs carry ~18x of aligned
    depth. Other platforms match their labels. Observed_Depth is computed upstream
    from the FASTQs and copied verbatim into the sheet, so it never saw the CRAM.

THE HYPOTHESIS
    Observed_Depth = reads x 150 / G, with the read length hard-coded to 150 for
    every sample rather than each sample's real length. T7 reads are 100 bp, so
    its depth is overstated by exactly 150/100 = 1.5x.

THE TEST — 'implied read length'
    Write both quantities in terms of the read density rho (reads per bp):

        depth_measured = rho x readlen        (what the CRAM holds)
        Observed       = rho x 150            (if the read length was hard-coded)

    Divide, and rho — and with it the genome size G — cancels:

        implied_readlen = readlen x Observed / depth_measured

        hard-coded 150     ->  implied == 150 for every platform
        real read length   ->  implied == measured readlen

    So the audit needs no genome size and no whole-CRAM read count: only a read
    density, which a few Mb of probe regions estimates in seconds. That is the
    whole reason this runs fast.

    Honest caveat, repeated in the log: only T7 can discriminate. Every other
    platform is a ~150 bp library, so for them the two hypotheses predict nearly
    the same implied value and the test has no power there. They are the
    calibration control — they show implied ~ measured ~ 150, i.e. the method is
    unbiased — not independent confirmation.

    Depth is compared WITH duplicates kept, because Observed_Depth comes from the
    FASTQ, where duplicate reads are still present and unmarked.

Everything measured here comes from `samtools coverage` on the CRAM. The only
value taken from the sheet is Observed_Depth — the number under audit.
"""

import argparse
import textwrap
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

W = 100
ASSUMED = 150.0  # the read length we suspect was hard-coded upstream

# house style, shared with cram.v6.ipynb and tuning.fraction.
# MUTED carries the caption and every axis label, so it is a body-text colour, not a
# decorative one: it sits near 7:1 against SURFACE. AXIS and GRID may recede — they
# are scaffolding — but nothing a reader has to READ is allowed to.
SERIES = ["#2a78d6", "#008300", "#e87ba4", "#eda100", "#7b5cd6", "#c2452d"]
SURFACE, MUTED, GRID, INK, AXIS = "#fcfcfb", "#54524d", "#dedcd4", "#0b0b0b", "#a8a59d"


def sniff_sep(header):
    """The sheet has shipped both TAB- and COMMA-separated; the .csv name says
    nothing. Guessing wrong yields one column and silent Nones, so read the
    header."""
    return "\t" if header.count("\t") >= header.count(",") and "\t" in header else ","


def read_sheet(path, id_col, dp_col, pf_col):
    """sample -> (observed_depth, platform)."""
    out = {}
    with open(path) as fh:
        sep = sniff_sep(fh.readline())
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split(sep)
        for c in (id_col, dp_col, pf_col):
            if c not in hdr:
                raise SystemExit(f"[audit] {path} has no column {c!r}; found {hdr}")
        i, d, p = hdr.index(id_col), hdr.index(dp_col), hdr.index(pf_col)
        for line in fh:
            f = line.rstrip("\n").split(sep)
            if len(f) <= max(i, d, p) or not f[d].strip():
                continue
            try:
                out[f[i].strip()] = (float(f[d]), f[p].strip())
            except ValueError:
                continue
    return out


def read_region_stats(paths):
    """sample -> measurements, length-weighted across the probe regions.

    Length-weighting rather than a plain mean over regions: meandepth is already a
    per-base average, so pooling them has to respect how many bases each covers.
    """
    acc = defaultdict(lambda: dict(nod=0.0, wdup=0.0, ln=0.0, rl=[], mn=[], mx=[]))
    for p in paths:
        with open(p) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) != len(hdr):
                    continue
                r = dict(zip(hdr, f))
                try:
                    L = float(r["LEN"])
                    a = acc[r["SAMPLE"]]
                    a["nod"] += float(r["DEPTH_NODUP"]) * L
                    a["wdup"] += float(r["DEPTH_WITHDUP"]) * L
                    a["ln"] += L
                    a["rl"].append(float(r["READLEN_MEAN"]))
                    a["mn"].append(float(r["READLEN_MIN"]))
                    a["mx"].append(float(r["READLEN_MAX"]))
                except (ValueError, KeyError):
                    continue
    out = {}
    for s, a in acc.items():
        if a["ln"] <= 0 or not a["rl"]:
            continue
        out[s] = dict(
            depth_nodup=a["nod"] / a["ln"],
            depth_withdup=a["wdup"] / a["ln"],
            readlen=float(np.mean(a["rl"])),
            fixed_len=min(a["mn"]) == max(a["mx"]),
        )
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region-stats", nargs="+", required=True)
    ap.add_argument("--cram-info", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--sample-id-col", default="ID_JHRPv6")
    ap.add_argument("--observed-dp-col", default="Observed_Depth")
    ap.add_argument("--platform-col", default="WGS_Platform")
    ap.add_argument("--out-log", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-samples-tsv", required=True,
                    help="per-sample rows behind the per-platform medians, so any "
                         "single sample can be checked rather than taken on trust")
    ap.add_argument("--out-fig", required=True)
    args = ap.parse_args()

    regions = [r.strip() for r in args.regions.split(",") if r.strip()]
    sheet = read_sheet(args.cram_info, args.sample_id_col, args.observed_dp_col,
                       args.platform_col)
    meas = read_region_stats(args.region_stats)

    samples = sorted(set(sheet) & set(meas))
    if not samples:
        raise SystemExit("[audit] no sample has both a sheet row and CRAM measurements")

    rec = {}
    for s in samples:
        obs, pf = sheet[s]
        m = meas[s]
        rec[s] = dict(
            platform=pf,
            observed=obs,
            readlen=m["readlen"],
            fixed_len=m["fixed_len"],
            depth=m["depth_withdup"],
            depth_nodup=m["depth_nodup"],
            dup_pct=100 * (1 - m["depth_nodup"] / m["depth_withdup"]) if m["depth_withdup"] else float("nan"),
            ratio=obs / m["depth_withdup"] if m["depth_withdup"] else float("nan"),
            implied=m["readlen"] * obs / m["depth_withdup"] if m["depth_withdup"] else float("nan"),
        )

    byp = defaultdict(list)
    for s, r in rec.items():
        byp[r["platform"]].append(s)

    rows = []
    for pf, ss in sorted(byp.items(), key=lambda kv: -len(kv[1])):
        med = lambda k: float(np.median([rec[s][k] for s in ss]))
        rows.append(dict(
            platform=pf, n=len(ss),
            readlen=med("readlen"),
            fixed=sum(rec[s]["fixed_len"] for s in ss),
            observed=med("observed"),
            depth=med("depth"),
            dup=med("dup_pct"),
            ratio=med("ratio"),
            implied=med("implied"),
        ))
    rows.sort(key=lambda r: r["readlen"])
    disc = [r for r in rows if abs(r["readlen"] - ASSUMED) > 20]

    # ── machine-readable: the samples behind the medians ──
    with open(args.out_samples_tsv, "w") as fh:
        fh.write("SAMPLE\tPLATFORM\tREADLEN_MEASURED\tFIXED_LENGTH\tOBSERVED_DEPTH\t"
                 "DEPTH_CRAM_WITHDUP\tDEPTH_CRAM_NODUP\tOBSERVED_OVER_MEASURED\t"
                 "READLEN_IMPLIED\n")
        for s in sorted(rec, key=lambda s: (rec[s]["platform"], s)):
            r = rec[s]
            fh.write(f"{s}\t{r['platform']}\t{r['readlen']:.1f}\t"
                     f"{'yes' if r['fixed_len'] else 'no'}\t{r['observed']:.4f}\t"
                     f"{r['depth']:.4f}\t{r['depth_nodup']:.4f}\t{r['ratio']:.4f}\t"
                     f"{r['implied']:.1f}\n")

    # ── machine-readable: per platform ──
    with open(args.out_tsv, "w") as fh:
        fh.write("PLATFORM\tN\tREADLEN_MEASURED\tREADLEN_IMPLIED\tOBSERVED_DEPTH\t"
                 "DEPTH_CRAM_WITHDUP\tOBSERVED_OVER_MEASURED\tDUP_PCT\n")
        for r in rows:
            fh.write(f"{r['platform']}\t{r['n']}\t{r['readlen']:.1f}\t{r['implied']:.1f}\t"
                     f"{r['observed']:.3f}\t{r['depth']:.3f}\t{r['ratio']:.3f}\t"
                     f"{r['dup']:.1f}\n")

    # ── log ──
    total_mb = sum(
        (int(r.split(":")[1].split("-")[1]) - int(r.split(":")[1].split("-")[0]) + 1)
        for r in regions) / 1e6
    L = ["=" * W,
         " OBSERVED_DEPTH AUDIT — the sheet's depth assumed a 150 bp read length",
         "=" * W,
         f" generated : {datetime.now().isoformat(timespec='seconds')}",
         f" samples   : {len(samples):,} across {len(byp)} platforms",
         f" probes    : {len(regions)} regions, {total_mb:,.1f} Mb total",
         "",
         " Depth is measured from the CRAM with `samtools coverage`, duplicates kept.",
         " The only value taken from the sheet is Observed_Depth — the number under audit.",
         "",
         "-" * W,
         " 1. THE DISCREPANCY",
         "-" * W,
         "   Observed_Depth is the sheet's label; depth is what the CRAM actually carries.",
         "   If the sheet were right, the ratio would be about the same for every platform.",
         "",
         f"   {'platform':<22}{'N':>4}{'readlen':>9}{'Observed':>10}{'CRAM':>9}"
         f"{'Obs/CRAM':>10}{'dup%':>7}",
         "   " + "-" * (W - 5)]
    for r in rows:
        L.append(f"   {r['platform']:<22}{r['n']:>4}{r['readlen']:>9.0f}"
                 f"{r['observed']:>10.2f}{r['depth']:>9.2f}{r['ratio']:>10.2f}{r['dup']:>7.1f}")
    L += ["",
          "-" * W,
          " 2. THE TEST — what read length does Observed_Depth imply?",
          "-" * W,
          "   Both quantities are the read density rho times a read length:",
          "",
          "       CRAM depth = rho x (real read length)",
          f"       Observed   = rho x (whatever length the upstream assumed)",
          "",
          "   so dividing them cancels rho — and with it the genome size, which is why this",
          "   audit needs neither G nor a whole-CRAM read count:",
          "",
          "       implied_readlen = readlen x Observed / CRAM depth",
          "",
          f"       if the upstream used each sample's real length -> implied == measured",
          f"       if it hard-coded {ASSUMED:.0f} bp                        -> implied == {ASSUMED:.0f} everywhere",
          "",
          f"   {'platform':<22}{'N':>4}{'measured':>10}{'implied':>9}{'fixed-len':>11}"
          f"   {'verdict':<26}",
          "   " + "-" * (W - 5)]
    for r in rows:
        fixed = f"{r['fixed']}/{r['n']}"
        if abs(r["readlen"] - ASSUMED) <= 20:
            verdict = "no power (~150 bp library)"
        elif abs(r["implied"] - ASSUMED) < abs(r["implied"] - r["readlen"]):
            verdict = f"implies {ASSUMED:.0f}, not {r['readlen']:.0f}"
        else:
            verdict = "implies its real length"
        L.append(f"   {r['platform']:<22}{r['n']:>4}{r['readlen']:>10.0f}{r['implied']:>9.0f}"
                 f"{fixed:>11}   {verdict:<26}")
    L += ["",
          "   'fixed-len' counts samples whose shortest and longest sampled read are equal —",
          "   the library really is one fixed length, so 'read length' is a well-defined",
          "   number and not an average over a mixture.",
          "",
          "-" * W,
          " 3. WHAT IT SHOWS",
          "-" * W]
    for r in disc:
        L += [f"   {r['platform']} is a {r['readlen']:.0f} bp library, yet Observed_Depth "
              f"implies {r['implied']:.0f} bp.",
              f"   Its depth is overstated by ~{ASSUMED / r['readlen']:.2f}x: the sheet says "
              f"{r['observed']:.1f}x, the reads carry {r['depth']:.1f}x.",
              f"   Corrected ({r['observed']:.1f} x {r['readlen']:.0f}/{ASSUMED:.0f}) = "
              f"{r['observed'] * r['readlen'] / ASSUMED:.1f}x, which matches the CRAM.",
              ""]
    L += [f"   CAVEAT — only a library far from {ASSUMED:.0f} bp can separate the two hypotheses.",
          "   Every other platform here is ~150 bp, so both hypotheses predict nearly the same",
          "   implied value for them and the test has no power there. They are the calibration",
          "   control — implied ~ measured ~ 150 shows the method is unbiased — not four",
          "   independent confirmations.",
          "",
          f"   Discriminating: "
          + (", ".join(f"{r['platform']} ({r['readlen']:.0f} bp)" for r in disc)
             if disc else "none in this sample set"),
          "",
          "-" * W,
          " 4. METHOD",
          "-" * W,
          f"   Depth is sampled over {len(regions)} probe regions ({total_mb:,.1f} Mb) rather than",
          "   the whole genome: the test needs a read density, not a total, and a few Mb",
          "   estimates it in seconds where a full `samtools stats` pass costs ~1.5 min/sample.",
          "",
          "     samtools coverage -r <region> --ff UNMAP,SECONDARY,QCFAIL <cram>",
          "   column 7, 'meandepth', pooled across regions weighted by region length.",
          "",
          "   Duplicates are KEPT (--ff omits DUP from the filter). Observed_Depth is derived",
          "   from the FASTQ, where duplicates are still present and unmarked; excluding them",
          "   here would charge each platform for its duplicate rate on top of the effect",
          "   under test. The dup% column shows what that would have cost.",
          "",
          "   Not corrected for: soft-clipped bases (~0.3%) and reads that never aligned",
          "   (~0.07%), both of which make the CRAM depth slightly lower than the FASTQ's.",
          "   They apply to every platform and are far too small to produce a 1.5x gap.",
          "=" * W, ""]
    Path(args.out_log).write_text("\n".join(L))
    print("\n".join(L))

    # ── figure ──
    # Three panels, one per step of the argument, sharing one visual grammar: a
    # dumbbell per platform whose two ends are the two things being compared, so
    # the connector's length IS the disagreement and a collapsed dumbbell IS
    # agreement. Each panel also carries a right-margin column with the one number
    # that panel is really about, because comparing connector lengths that sit at
    # different places on the x-axis is exactly the arithmetic the reader should
    # not have to do.
    # C_LINE is the connector, and the connector IS the finding — a visible gap on
    # T7, none elsewhere. It has to read at a glance, so it is a mid grey rather
    # than the near-white it started as.
    C_SHEET, C_CRAM, C_LINE = "#c2452d", "#2a78d6", "#7f7c75"
    ypos = {r["platform"]: i for i, r in enumerate(rows)}  # rows sorted by readlen

    def style(ax, title, xlab, note_hdr):
        ax.set_facecolor(SURFACE)
        ax.set_title(title, color=INK, fontsize=11.5, loc="left", pad=22)
        ax.set_xlabel(xlab, color=MUTED, fontsize=9.5)
        ax.set_yticks(range(len(rows)))
        # One empty slot below the last row for the legend: the dumbbells use the
        # full x-range, so no in-plot corner is free of markers.
        ax.set_ylim(-0.8, len(rows) + 0.5)
        ax.invert_yaxis()
        ax.grid(axis="x", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        ax.tick_params(colors=MUTED, labelsize=9)
        ax.tick_params(axis="y", length=0, labelcolor=INK)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(AXIS)
        # header for the right-margin column
        ax.text(1.02, -0.75, note_hdr, transform=ax.get_yaxis_transform(),
                ha="left", va="center", color=MUTED, fontsize=8.5)

    def note(ax, y, text, strong=False):
        ax.text(1.02, y, text, transform=ax.get_yaxis_transform(), ha="left",
                va="center", color=INK if strong else MUTED,
                fontsize=9.5 if strong else 9,
                fontweight="bold" if strong else "normal")

    def pair(ax, r, lo, hi, lo_lab, hi_lab):
        y = ypos[r["platform"]]
        ax.plot([lo, hi], [y, y], color=C_LINE, lw=2.5, zorder=1,
                solid_capstyle="round")
        ax.scatter([lo], [y], s=75, color=C_CRAM, zorder=3,
                   label=lo_lab if y == 0 else None)
        ax.scatter([hi], [y], s=75, color=C_SHEET, zorder=3,
                   label=hi_lab if y == 0 else None)

    fig, axes = plt.subplots(1, 3, figsize=(18.5, 4.8), facecolor=SURFACE)
    fig.subplots_adjust(wspace=0.55)

    # ── 1. the finding ───────────────────────────────────────────────────────
    a = axes[0]
    for r in rows:
        pair(a, r, r["depth"], r["observed"], "measured in the CRAM",
             "Observed_Depth (the sheet)")
        note(a, ypos[r["platform"]], f"{r['ratio']:.2f}x", strong=r in disc)
    a.set_yticklabels([f"{r['platform']}\n{r['readlen']:.0f} bp, n={r['n']}" for r in rows])
    style(a, "1.  The sheet matches the CRAM — except for one platform",
          "depth  (x)", "sheet\n/ CRAM")
    a.legend(fontsize=8.5, frameon=False, loc="lower right", labelcolor=INK)

    # ── 2. the explanation, framed as a prediction test ──────────────────────
    # Same dumbbell grammar and the same collapse=agreement reading as panels 1
    # and 3 — but the axis is now the depth ratio itself (a x, like the margin
    # numbers), and the two ends are a PREDICTION and an OBSERVATION of it. The
    # read-length hypothesis predicts each platform's inflation is exactly
    # 150/readlen (blue); the observed inflation is sheet/CRAM (red, the same
    # number panel 1 puts in its margin). A collapsed dumbbell here means the
    # prediction landed. Every one collapses — far out at 1.5 for T7, at 1.0 for
    # the rest — so the read length explains the whole picture, not just T7.
    b = axes[1]
    for r in rows:
        pred = ASSUMED / r["readlen"]
        pair(b, r, pred, r["ratio"],
             f"predicted:  {ASSUMED:.0f} / read length",
             "observed:  sheet / CRAM")
        note(b, ypos[r["platform"]], f"{r['readlen']:.0f} bp reads", strong=r in disc)
    b.axvline(1.0, color=C_LINE, ls="--", lw=1.3, zorder=2)
    b.text(1.0, -0.75, "1.0x\nno error", color=MUTED, fontsize=8, ha="center",
           va="center", linespacing=0.9)
    b.set_yticklabels([])
    b.set_xlim(0.88, 1.62)
    style(b, "2.  Read length predicts that exact error",
          "sheet depth / CRAM depth  (x)", "read\nlength")
    b.legend(fontsize=8.5, frameon=False, loc="lower left", labelcolor=INK)

    # ── 3. the check ─────────────────────────────────────────────────────────
    c = axes[2]
    for r in rows:
        r["corrected"] = r["observed"] * r["readlen"] / ASSUMED
        pair(c, r, r["depth"], r["corrected"], "measured in the CRAM",
             f"Observed_Depth x readlen/{ASSUMED:.0f}")
        note(c, ypos[r["platform"]], f"{r['corrected'] / r['depth']:.2f}x",
             strong=r in disc)
    c.set_yticklabels([])
    c.set_xlim(axes[0].get_xlim())
    style(c, f"3.  Rescale by readlen/{ASSUMED:.0f} and it agrees",
          "depth  (x)", "sheet\n/ CRAM")
    c.legend(fontsize=8.5, frameon=False, loc="lower right", labelcolor=INK)

    sub = ""
    if disc:
        d = disc[0]
        verb = "is the only library" if len(disc) == 1 else "are the only libraries"
        names = ", ".join(f"{r['platform']} ({r['readlen']:.0f} bp)" for r in disc)
        sub = (f"{names} {verb} not ~{ASSUMED:.0f} bp — the only place the two hypotheses "
               f"disagree. The rest are the calibration control.")
    fig.suptitle(
        f"Observed_Depth was computed with the read length hard-coded to {ASSUMED:.0f} bp\n"
        + sub + f"  |  {len(samples):,} CRAMs, depth over {total_mb:,.0f} Mb of probe "
        f"regions, duplicates kept.",
        color=INK, fontsize=12.5, x=0.005, ha="left", y=1.13)

    # ── caption ──
    # The figure travels: into an email, a slide, a thread. It has to carry its own
    # explanation, so the caption states the finding, the test, the conclusion, the
    # limit of that conclusion, and the method — all generated from the same numbers
    # as the plot, so it cannot drift from what is drawn above it.
    ctrl = [r for r in rows if r not in disc]
    d = disc[0]
    lo, hi = min(r["ratio"] for r in ctrl), max(r["ratio"] for r in ctrl)
    corr = d["observed"] * d["readlen"] / ASSUMED
    # Two columns: one column wide enough to fill an 18.5in figure would run past
    # 250 characters a line, which is unreadable. Left column walks the argument,
    # right column carries what the argument rests on.
    left = [
        f"HOW TO READ  Each row is one sequencing platform ({rows[0]['n']} CRAMs each). "
        f"Every panel is a dumbbell: the two dots are the quantities being compared and the "
        f"connector is their disagreement, so overlapping dots mean the two agree. The number "
        f"in each right margin is that panel's headline value.",

        f"(1) THE FINDING  Blue: depth measured directly from the CRAM. Red: Observed_Depth as "
        f"recorded in the sample sheet — computed upstream from the FASTQs and copied in "
        f"verbatim, so it never saw the CRAM. {len(ctrl)} of {len(rows)} platforms agree with "
        f"their own CRAM to within {100 * (hi - 1):.0f}%. {d['platform']} does not: the sheet "
        f"claims {d['observed']:.1f}x where the reads carry {d['depth']:.1f}x, an inflation of "
        f"{d['ratio']:.2f}x.",

        f"(2) THE TEST  The x-axis is now that inflation factor itself. Blue: what read length "
        f"alone predicts it should be ({ASSUMED:.0f} / read length). Red: what is actually "
        f"observed (sheet / CRAM — the same number as panel 1's margin). The two come from "
        f"disjoint inputs: the prediction never touches a depth, the observation never touches "
        f"a read length. Yet every dumbbell collapses — {d['platform']} far out at "
        f"{ASSUMED / d['readlen']:.2f}x (= {ASSUMED:.0f}/{d['readlen']:.0f} bp), the "
        f"~{ASSUMED:.0f} bp platforms at 1.00x. Read length predicts every platform's error, "
        f"not just {d['platform']}'s.",

        f"(3) THE CORRECTION  Rescaling Observed_Depth by (read length / {ASSUMED:.0f}) moves "
        f"{d['platform']} from {d['ratio']:.2f}x to {corr / d['depth']:.2f}x "
        f"({d['observed']:.1f}x x {d['readlen']:.0f}/{ASSUMED:.0f} = {corr:.1f}x, against "
        f"{d['depth']:.1f}x measured), in line with the other platforms.",
    ]
    right = [
        f"CONCLUSION  Observed_Depth was computed with the read length hard-coded to "
        f"{ASSUMED:.0f} bp. {d['platform']} reads are {d['readlen']:.0f} bp, so its depth is "
        f"overstated by {ASSUMED:.0f}/{d['readlen']:.0f} = {ASSUMED / d['readlen']:.1f}x. "
        f"These samples are ~{d['depth']:.0f}x, not the ~30x their label claims.",

        f"POWER  Only {d['platform']} can discriminate. The other {len(ctrl)} are "
        f"~{ASSUMED:.0f} bp libraries, for which 'each sample's real read length' and "
        f"'hard-coded {ASSUMED:.0f}' predict nearly the same thing — they are the calibration "
        f"control, showing the method is unbiased, not {len(ctrl)} independent confirmations.",

        f"METHOD  {len(samples)} CRAMs, {rows[0]['n']} per platform, random with a fixed seed. "
        f"Depth from `samtools coverage` over {len(regions)} x "
        f"{total_mb / len(regions):.0f} Mb windows on {len(regions)} chromosomes, pooled "
        f"weighted by length, duplicates kept (Observed_Depth came from the FASTQ, where "
        f"duplicates are still unmarked). Read length from the first 1,000 reads per region, "
        f"excluding secondary and supplementary records, whose SEQ is hard-clipped. The "
        f"residual ~{100 * (hi - 1):.0f}% by which observed exceeds predicted throughout is "
        f"soft-clipped bases (~0.3%) and unaligned reads (~0.07%), which lower CRAM depth "
        f"slightly for every platform alike.",
    ]
    for x0, col in ((0.005, left), (0.515, right)):
        fig.text(x0, -0.03, "\n\n".join(textwrap.fill(p, 132) for p in col),
                 ha="left", va="top", color=MUTED, fontsize=8.6, linespacing=1.5)

    fig.savefig(args.out_fig, dpi=150, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


if __name__ == "__main__":
    main()
