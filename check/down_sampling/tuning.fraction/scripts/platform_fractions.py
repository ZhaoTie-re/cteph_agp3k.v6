#!/usr/bin/env python3
"""Per-platform down-sampling fractions, all measured from the CRAMs.

down_sampling.nf keeps a fraction of each sample's reads so that it matches the
depth of the cohort's reference platform. A single fraction for a whole
Target_Depth group is wrong, because a group is not one population: its platforms
sit at very different depths in the regions we analyse. So one fraction is
derived per platform.

    fraction(P) = median CRAM depth of the baseline platform
                  ------------------------------------------
                  median CRAM depth of platform P

Depth is the aligned read depth in the target regions, measured with
`samtools coverage` on the CRAMs — the very read pool `samtools view -s`
subsamples. Nothing here uses the sample sheet's Observed_Depth: that number is
genome-wide, computed upstream from the FASTQs, and is not what the down-sampling
acts on.

fraction > 1 would mean adding reads, which is impossible; such a platform is
already at or below the baseline and is reported with fraction 1.0 (no
down-sampling). Samples without a CRAM cannot be measured and are excluded — the
count is reported so the exclusion is never silent.
"""

import argparse
import glob
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

W = 100


def sniff_sep(header):
    """The sample sheet has shipped both TAB- and COMMA-separated (the .csv name
    says nothing). Splitting on the wrong one silently yields one column and all
    lookups return None, so decide from the header rather than hard-coding."""
    return "\t" if header.count("\t") >= header.count(",") and "\t" in header else ","


def read_meta(csv, id_col, td_col, pf_col):
    """sample -> (target_depth_label, platform). No depth is taken from here."""
    meta = {}
    with open(csv) as fh:
        header = fh.readline().rstrip("\n")
        sep = sniff_sep(header)
        hdr = [c.strip() for c in header.split(sep)]
        missing = [c for c in (id_col, td_col, pf_col) if c not in hdr]
        if missing:
            raise SystemExit(f"[platform_fractions] {csv} (separator {sep!r}) lacks "
                             f"column(s) {missing}. Found: {hdr}")
        i_id, i_td, i_pf = hdr.index(id_col), hdr.index(td_col), hdr.index(pf_col)
        for line in fh:
            f = [x.strip() for x in line.rstrip("\n").split(sep)]
            try:
                meta[f[i_id]] = (f[i_td], f[i_pf])
            except IndexError:
                pass
    if not meta:
        raise SystemExit(f"[platform_fractions] no usable rows parsed from {csv}")
    return meta


def read_cram_dp(paths, wanted):
    """(SAMPLE, REGION, LEN, DEPTH_NODUP, DEPTH_WITHDUP, READLEN) over `wanted`
    regions -> per-sample length-weighted depth / with-dup depth / read length."""
    acc = defaultdict(lambda: [0.0, 0.0, 0, 0, 0.0, 0])
    for p in paths:
        with open(p) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 5 or f[1] not in wanted:
                    continue
                s, ln, nod, wd = f[0], f[2], f[3], f[4]
                rl = f[5] if len(f) > 5 else "NA"
                try:
                    ln = int(ln)
                    acc[s][0] += float(nod) * ln
                    acc[s][2] += ln
                    if wd != "NA":
                        acc[s][1] += float(wd) * ln
                        acc[s][3] += ln
                    if rl != "NA":
                        acc[s][4] += float(rl) * ln
                        acc[s][5] += ln
                except ValueError:
                    pass
    depth = {s: v[0] / v[2] for s, v in acc.items() if v[2]}
    withdup = {s: v[1] / v[3] for s, v in acc.items() if v[3]}
    readlen = {s: v[4] / v[5] for s, v in acc.items() if v[5]}
    return depth, withdup, readlen


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cram-dp", nargs="+", required=True)
    ap.add_argument("--regions", required=True, help="comma-separated target regions")
    ap.add_argument("--cram-info", required=True)
    ap.add_argument("--refined-fam", required=True)
    ap.add_argument("--baseline-platform", required=True)
    ap.add_argument("--sample-id-col", default="ID_JHRPv6")
    ap.add_argument("--target-dp-col", default="Target_Depth")
    ap.add_argument("--platform-col", default="WGS_Platform")
    ap.add_argument("--regions-label", default="target regions")
    ap.add_argument("--out-log", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-fig", required=True)
    args = ap.parse_args()

    regions = {r.strip() for r in args.regions.split(",") if r.strip()}
    keep = {l.split()[1] for l in open(args.refined_fam) if l.strip()}
    meta = read_meta(args.cram_info, args.sample_id_col, args.target_dp_col,
                     args.platform_col)
    depth, withdup, readlen = read_cram_dp(args.cram_dp, regions)

    # group the measured samples by platform, and count what could not be measured
    byp = defaultdict(list)
    for s in keep:
        if s in meta and s in depth:
            byp[meta[s][1]].append(s)
    measured = sum(len(v) for v in byp.values())
    unmeasured = defaultdict(int)
    for s in keep:
        if s in meta and s not in depth:
            unmeasured[meta[s][1]] += 1

    if args.baseline_platform not in byp:
        raise SystemExit(
            f"[platform_fractions] baseline platform {args.baseline_platform!r} has no "
            f"measured sample. Present: {sorted(byp)}")

    base_depth = float(np.median([depth[s] for s in byp[args.baseline_platform]]))
    if base_depth <= 0:
        raise SystemExit("[platform_fractions] baseline depth is zero")

    # one row per platform
    rows = []
    for pf, ss in byp.items():
        d = np.array([depth[s] for s in ss])
        raw = base_depth / float(np.median(d))
        rows.append(dict(
            platform=pf,
            group=meta[ss[0]][0],
            n=len(ss),
            n_missing=unmeasured.get(pf, 0),
            depth_med=float(np.median(d)),
            depth_mean=float(d.mean()),
            readlen=float(np.median([readlen[s] for s in ss if s in readlen])) if any(
                s in readlen for s in ss) else float("nan"),
            dup=float(np.median([100 * (1 - depth[s] / withdup[s]) for s in ss
                                 if s in withdup and withdup[s] > 0])) if any(
                s in withdup for s in ss) else float("nan"),
            raw_fraction=raw,
            fraction=min(raw, 1.0),
            capped=raw > 1.0,
            baseline=pf == args.baseline_platform,
        ))
    rows.sort(key=lambda r: (not r["baseline"], -r["n"]))

    # ── machine-readable table for down_sampling.nf ──
    with open(args.out_tsv, "w") as fh:
        fh.write("PLATFORM\tTARGET_DEPTH\tN\tDEPTH_MEDIAN\tREADLEN\tKEEP_FRACTION\tNOTE\n")
        for r in rows:
            note = ("baseline" if r["baseline"] else
                    "at_or_below_baseline_no_downsampling" if r["capped"] else "")
            fh.write(f"{r['platform']}\t{r['group']}\t{r['n']}\t{r['depth_med']:.4f}\t"
                     f"{r['readlen']:.0f}\t{r['fraction']:.4f}\t{note}\n")

    # ── log ──
    L = ["=" * W,
         " PER-PLATFORM DOWN-SAMPLING FRACTIONS",
         "=" * W,
         f" generated : {datetime.now().isoformat(timespec='seconds')}",
         f" regions   : {', '.join(sorted(regions))}",
         f" baseline  : {args.baseline_platform}   (median CRAM depth {base_depth:.2f})",
         f" samples   : refined_core {len(keep):,};"
         f" measured {measured:,}; without CRAM {sum(unmeasured.values()):,}",
         "",
         " Every depth below is the aligned read depth of the target regions measured on",
         " the CRAM with `samtools coverage` — the read pool `samtools view -s` subsamples.",
         " The sample sheet's Observed_Depth is deliberately not used: it is genome-wide,",
         " derived upstream from the FASTQs, and is not what the down-sampling acts on.",
         "",
         "-" * W,
         " 1. THE FRACTIONS",
         "-" * W,
         "   fraction(P) = median depth(baseline) / median depth(P)",
         "",
         f"   {'platform':<24}{'grp':>5}{'N':>6}{'no CRAM':>9}{'readlen':>9}{'dup%':>7}"
         f"{'depth med':>11}{'depth mean':>12}{'fraction':>10}   note",
         "   " + "-" * (W - 5)]
    for r in rows:
        note = ("<- baseline" if r["baseline"] else
                f"raw {r['raw_fraction']:.3f} -> capped, already <= baseline" if r["capped"]
                else "")
        L.append(f"   {r['platform']:<24}{r['group']:>5}{r['n']:>6}{r['n_missing']:>9}"
                 f"{r['readlen']:>9.0f}{r['dup']:>7.1f}{r['depth_med']:>11.2f}"
                 f"{r['depth_mean']:>12.2f}{r['fraction']:>10.3f}   {note}")

    L += ["",
          "-" * W,
          " 2. HOW TO READ IT",
          "-" * W,
          "   A fraction of 0.50 means: keep half of that platform's reads and it will sit",
          "   at the baseline depth. A fraction of 1.000 marked 'already <= baseline' means",
          "   the platform is at or below the baseline already, so there is nothing to take",
          "   away — reads cannot be added.",
          "",
          "   The fractions differ a lot BETWEEN platforms of the same Target_Depth group.",
          "   That is the point of doing this per platform: a group is not one population,",
          "   and a single group-wide fraction would over-down-sample some platforms while",
          "   under-down-sampling others.",
          "",
          "   'readlen' and 'dup' are reported as context only; they do not enter the",
          "   fraction. Note that `samtools coverage` excludes duplicate-flagged reads, so",
          "   a platform's depth here is already net of its duplicates.",
          "",
          "-" * W,
          " 3. METHOD",
          "-" * W,
          "   Per sample and per region:",
          "       samtools coverage -r <region> --reference <fasta> <cram>",
          "   column 7, 'meandepth' = sum of per-base read depth / region length, i.e. an",
          "   average over every base of the region.",
          "",
          "   At samtools' defaults (verified, not assumed):",
          "     excluded : UNMAP, SECONDARY, QCFAIL, DUP   (the default --ff)",
          "     counted  : SUPPLEMENTARY reads             (~0.1% of depth here)",
          "     no mapping-quality (-q 0) or base-quality (-Q 0) threshold",
          "",
          "   Across regions, per sample, length-weighted:",
          "       depth = sum(meandepth_r * len_r) / sum(len_r)",
          "",
          "   'dup' is a second pass with duplicates kept (--ff UNMAP,SECONDARY,QCFAIL);",
          "   per sample dup = 1 - depth_nodup / depth_withdup, then the platform median.",
          "",
          "   Samples of refined_core without a CRAM cannot be measured and are dropped;",
          "   the 'no CRAM' column shows how many per platform. The medians rest on the",
          "   measured samples only.",
          "=" * W, ""]
    Path(args.out_log).write_text("\n".join(L))
    print("\n".join(L))

    # ── figure ──
    # Two panels on a shared x: the measured depth, and what the fractions do to it.
    # House style follows cram.v6.ipynb — horizontal boxes, platforms direct-labelled
    # on the y-axis, so identity never rests on colour alone.
    SERIES = ["#2a78d6", "#008300", "#e87ba4", "#eda100", "#7b5cd6", "#c2452d"]
    BASE_C = "#0b0b0b"
    # MUTED carries axis labels and the per-row medians, so it is a body-text
    # colour (~7:1 on SURFACE), not a decorative one. GRID and AXIS are scaffolding
    # and may recede; anything a reader has to READ may not.
    SURFACE, MUTED, GRID, INK, AXIS = "#fcfcfb", "#54524d", "#dedcd4", "#0b0b0b", "#a8a59d"

    # Deepest platform on top: the reader's question is "who is furthest above the
    # baseline", and that ordering answers it without a second look.
    prows = sorted(rows, key=lambda r: r["depth_med"], reverse=True)
    pdata = [np.array([depth[s] for s in byp[r["platform"]]]) for r in prows]
    colors = [BASE_C if r["baseline"] else SERIES[j % len(SERIES)]
              for j, r in enumerate(prows)]

    def whisker_span(vals):
        """Where the box and its 1.5xIQR whiskers actually end.

        The x limits are taken from this rather than from the data range: with
        fliers hidden, letting the extremes set the axis would reserve a third of
        the width for marks that are not drawn — in this cohort HiSeqX runs to 45x
        on a handful of samples while its whisker ends near 22x.
        """
        q1, q3 = np.percentile(vals, [25, 75])
        iqr = q3 - q1
        inside = vals[(vals >= q1 - 1.5 * iqr) & (vals <= q3 + 1.5 * iqr)]
        return float(inside.min()), float(inside.max())

    def draw(ax, series, title):
        for i, (r, vals, color) in enumerate(zip(prows, series, colors)):
            # Boxes only, no fliers and no underlying points. The raw points were
            # drawn here once, but n runs from 27 to 1,758 across the rows, so cloud
            # density differed by two orders of magnitude and was read as spread.
            # Quartiles are comparable across rows; ink per sample is not.
            ax.boxplot([vals], positions=[i + 1], vert=False, widths=0.55,
                       showfliers=False, patch_artist=True, zorder=3,
                       boxprops=dict(facecolor=SURFACE, color=color, linewidth=1.8),
                       whiskerprops=dict(color=color, linewidth=1.4),
                       capprops=dict(color=color, linewidth=1.4),
                       medianprops=dict(color=INK, linewidth=2.2))
        ax.axvline(base_depth, color=BASE_C, ls="--", lw=1.4, zorder=1)
        ax.set_yticks(range(1, len(prows) + 1))
        ax.set_ylim(0.45, len(prows) + 0.55)
        ax.invert_yaxis()
        ax.set_xlim(XLO, XHI)
        ax.set_xticks(xticks)
        ax.set_facecolor(SURFACE)
        # No ax.set_title: a title starts at the AXES edge, which here is well right
        # of the panel's y labels, so it floats over the middle of the panel instead
        # of heading it. The titles are drawn as figure text beside the letters.
        ax.set_xlabel(f"CRAM aligned depth in the {args.regions_label}",
                      color=MUTED, fontsize=9.5)
        ax.grid(axis="x", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        ax.tick_params(colors=MUTED, labelsize=9)
        ax.tick_params(axis="y", length=0, labelcolor=INK)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(AXIS)

    def col_header(ax, text):
        """Header for the panel's left-hand column, over the y tick labels."""
        ax.text(-0.02, 1.015, text, transform=ax.transAxes,
                va="bottom", ha="right", color=MUTED, fontsize=8)

    def gutter(ax, series, header):
        """The median as a number, in a column of its own outside the panel.

        A box shows a median's position; the design is stated in depths, so the
        value has to be readable too. Set as a column with a header rather than as
        five floating labels.
        """
        # The header goes above the axes in axes coordinates. In data coordinates it
        # would have to sit beyond the last row, and the y axis is inverted, so
        # "beyond the last row" is the BOTTOM of the panel — under the x label.
        ax.text(1.02, 1.015, header, transform=ax.transAxes,
                va="bottom", ha="left", color=MUTED, fontsize=8)
        for i, vals in enumerate(series):
            ax.text(1.02, i + 1, f"{np.median(vals):.1f}x",
                    transform=ax.get_yaxis_transform(), va="center", ha="left",
                    color=INK, fontsize=8.5)

    # B — the same samples with their platform's fraction applied. Random read
    # subsampling scales depth linearly in expectation, so depth x fraction is the
    # projection; it is a design projection, not a re-measurement.
    bdata = [vals * r["fraction"] for r, vals in zip(prows, pdata)]

    spans = [whisker_span(v) for v in pdata + bdata]
    XLO = min(s[0] for s in spans) - 1.2
    XHI = max(s[1] for s in spans) + 1.2
    xticks = [t for t in range(0, 200, 5) if XLO <= t <= XHI]

    fig, axes = plt.subplots(1, 2, figsize=(13.4, 4.5), facecolor=SURFACE, sharex=True)

    # A — the raw measurement, nothing applied
    draw(axes[0], pdata, "Measured depth, per platform")
    axes[0].set_yticklabels([f"{r['platform']}   (n={r['n']:,})" for r in prows])
    col_header(axes[0], "platform")
    gutter(axes[0], pdata, "median")

    draw(axes[1], bdata, "Projected depth after applying keep_fraction")
    axes[1].set_yticklabels([f"x {r['fraction']:.3f}" + ("  (capped)" if r["capped"] else "")
                             for r in prows])
    col_header(axes[1], "keep_fraction")
    gutter(axes[1], bdata, "median")

    # The headline must not overclaim: a capped platform is NOT brought onto the
    # baseline, it stays where it is. Name the exceptions from the data.
    capped_names = [r["platform"] for r in prows if r["capped"] and not r["baseline"]]
    exception = ""
    if capped_names:
        verb = "sits" if len(capped_names) == 1 else "sit"
        exception = (f"  {', '.join(capped_names)} already {verb} below the baseline and "
                     f"stays there — reads cannot be added.")
    fig.suptitle(
        f"Per-platform down-sampling: every platform above the baseline is brought onto "
        f"{args.baseline_platform} ({base_depth:.2f}x, dashed)\n"
        f"refined_core, {measured:,} samples measured from CRAM"
        + (f"; {sum(unmeasured.values()):,} excluded for having no CRAM."
           if sum(unmeasured.values()) else ".")
        + exception
        + "\nBox, quartiles; whiskers, 1.5 x IQR; individual samples not shown. "
          "A and B share one x axis.",
        color=INK, fontsize=10.5, x=0.006, ha="left", va="top", y=0.985)

    # Placed, not fitted. The median column is fig-level text outside the axes, which
    # tight_layout does not measure and would therefore crop or crowd; the margins
    # here are budgeted for it (0.045 of the width per gutter) and for B's y labels.
    fig.subplots_adjust(left=0.175, right=0.945, top=0.735, bottom=0.155, wspace=0.46)
    # Letter and title on one line at the top-left of the whole panel block — labels,
    # plot and median column — rather than at the axes edge. LGUT is how far each
    # block's y labels reach left of its axes; A's are platform names and B's are
    # fractions, so the two are not the same width and one constant would leave one
    # letter hanging in the gap between the panels.
    LGUT = (0.163, 0.108)
    for ax_, letter, title, lg in zip(
            axes, "AB",
            ("Measured depth, per platform",
             "Projected depth after applying keep_fraction"), LGUT):
        p = ax_.get_position()
        fig.text(p.x0 - lg, p.y1 + 0.075, letter, fontsize=13, fontweight="bold",
                 color=INK, va="bottom", ha="left")
        fig.text(p.x0 - lg + 0.030, p.y1 + 0.078, title, fontsize=10.5,
                 color=INK, va="bottom", ha="left")
    fig.savefig(args.out_fig, dpi=200, facecolor=SURFACE)
    plt.close(fig)


if __name__ == "__main__":
    main()
