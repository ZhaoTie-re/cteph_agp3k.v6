#!/usr/bin/env python3
"""Is the signal an artefact of one sequencing platform?

Down-sampling removes the DEPTH difference between platforms. It does not remove
the platforms: a T7 read at baseline depth is still a 100 bp DNBSEQ read, and a
HiSeqX read a 151 bp Illumina one. So the depth result rules out depth and leaves
read length, chemistry and error profile untouched — and in this cohort those
track the phenotype exactly as depth did, because every control is HiSeqX and no
case is.

No platform holds both cases and controls, so the clean test — hold the platform
fixed, compare cases to controls — cannot be run here at all. That is how the
cohort was collected; no reprocessing repairs it. What is left is weaker:

  only_<P>   the effect estimated from ONE platform's cases (+ all controls).
             These share no case with each other, so they are independent, and
             Cochran's Q / I-squared over them is a real heterogeneity test.

  minus_<P>  the effect with one platform's cases removed.
             These are NESTED — any two share most of their samples — so no
             heterogeneity test may be run over them, and this script does not
             compute one. Worse, three of the four keep ~90% of the cases, which
             forces their estimates onto the full-cohort one whatever the truth
             is. Only the one dropping the dominant platform carries information.
             They are reported because "does it survive losing this platform" is
             a fair question, and labelled because a forest plot of them looks
             like four confirmations and is not.

Every estimate here comes from SAIGE — step 1 refit on the subset, then step 2 —
i.e. the main analysis run on those samples. Nothing here is a crude allele-count
test, so every OR on this page is comparable with the study's own.
"""

import argparse
import math
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2

W = 102
SURFACE, MUTED, GRID, INK, AXIS = "#fcfcfb", "#54524d", "#dedcd4", "#0b0b0b", "#a8a59d"
# Okabe-Ito, as in 03_report.
C_ALL, C_ONLY, C_MINUS = "#0b0b0b", "#D55E00", "#0072B2"
# A minus_* subset keeping more than this share of the cases is forced onto the
# full-cohort estimate by construction. The log's 'informative?' column uses it;
# the figure does not split the rows on it — its title says the rows are nested,
# and the 'cases kept' on every row is the number the reader weighs them by.
FORCED_FRAC = 0.60


def read_saige(path, variant):
    with open(path) as fh:
        hdr = fh.readline().split()
        idx = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            f = line.split()
            if len(f) >= len(hdr) and f[idx["MarkerID"]] == variant:
                return dict(beta=float(f[idx["BETA"]]), se=float(f[idx["SE"]]),
                            p=float(f[idx["p.value.NA"]]),
                            p_spa=float(f[idx["p.value"]]),
                            af=float(f[idx["AF_Allele2"]]),
                            n_case=int(f[idx["N_case"]]), n_ctrl=int(f[idx["N_ctrl"]]),
                            a1=f[idx["Allele1"]], a2=f[idx["Allele2"]])
    return None


def orci(beta, se):
    return math.exp(beta), math.exp(beta - 1.96 * se), math.exp(beta + 1.96 * se)


def heterogeneity(betas, ses):
    """Cochran's Q and I-squared over INDEPENDENT estimates.

    Valid only when the estimates share no samples. Applied to nested subsets it
    would return a near-zero I-squared no matter what, because the estimates are
    forced to agree by construction — a fabricated reassurance. The caller is
    responsible for only passing independent arms; this function cannot check it.
    """
    b = np.asarray(betas, float)
    w = 1.0 / np.asarray(ses, float) ** 2
    if len(b) < 2:
        return None
    bbar = float((w * b).sum() / w.sum())
    q = float((w * (b - bbar) ** 2).sum())
    df = len(b) - 1
    p = float(chi2.sf(q, df))
    i2 = max(0.0, 100.0 * (q - df) / q) if q > 0 else 0.0
    se_bar = math.sqrt(1.0 / w.sum())
    return dict(q=q, df=df, p=p, i2=i2, beta=bbar, se=se_bar)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assoc", action="append", required=True, metavar="TAG:PATH")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--platform-meta", required=True)
    ap.add_argument("--variant", required=True)
    ap.add_argument("--out-log", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-fig", required=True)
    args = ap.parse_args()

    man = {}
    with open(args.manifest) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5:
                man[f[0]] = dict(family=f[1], n_case=int(f[2]), n_ctrl=int(f[3]), note=f[4])

    meta = {}
    with open(args.platform_meta) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            r = dict(zip(hdr, f))
            meta[r["PLATFORM"].replace(" ", "_")] = (float(r["DEPTH_MEDIAN"]), float(r["READLEN"]))

    res = {}
    for spec in args.assoc:
        tag, path = spec.split(":", 1)
        r = read_saige(path, args.variant)
        if r:
            r.update(man.get(tag, {}))
            res[tag] = r

    if "all" not in res:
        raise SystemExit("[platform_report] the 'all' subset is missing; nothing to compare to")
    allr = res["all"]
    only = {t: r for t, r in res.items() if t.startswith("only_")}
    minus = {t: r for t, r in res.items() if t.startswith("minus_")}

    def pf_of(tag):
        return tag.split("_", 1)[1]

    def rl(tag):
        return meta.get(pf_of(tag), (float("nan"), float("nan")))

    het = heterogeneity([r["beta"] for r in only.values()],
                        [r["se"] for r in only.values()]) if len(only) > 1 else None

    # ── tsv ──
    with open(args.out_tsv, "w") as fh:
        fh.write("TAG\tFAMILY\tINDEPENDENT\tN_CASE\tN_CONTROL\tREADLEN\tDEPTH\tAF\t"
                 "BETA\tSE\tOR\tOR_L95\tOR_U95\tP_NONSPA\n")
        for tag in ["all"] + sorted(only) + sorted(minus):
            r = res.get(tag)
            if not r:
                continue
            o, lo, hi = orci(r["beta"], r["se"])
            d, rr = rl(tag) if tag != "all" else (float("nan"), float("nan"))
            ind = {"independent": "yes", "nested": "no"}.get(r.get("family", ""), "NA")
            fh.write(f"{tag}\t{r.get('family','')}\t{ind}\t{r['n_case']}\t{r['n_ctrl']}\t"
                     f"{rr:.0f}\t{d:.2f}\t{r['af']:.6f}\t{r['beta']:.6f}\t{r['se']:.6f}\t"
                     f"{o:.4f}\t{lo:.4f}\t{hi:.4f}\t{r['p']:.6E}\n")

    # ── log ──
    o, lo, hi = orci(allr["beta"], allr["se"])
    L = ["=" * W,
         " Is the signal an artefact of one sequencing platform?",
         "=" * W,
         f" generated : {datetime.now().isoformat(timespec='seconds')}",
         f" variant   : {args.variant}",
         "",
         " Every estimate below is SAIGE — step 1 refit on the subset, then step 2. That is",
         " the main analysis run on those samples, so every OR here is the study's own kind",
         " of OR: exp(BETA), adjusted for SEX + 10 PCs and the GRM. Nothing is a crude",
         " allele-count test.",
         "",
         "-" * W,
         " 1. WHAT THIS COHORT CANNOT ANSWER",
         "-" * W,
         "   No platform holds both cases and controls — every control is HiSeqX, no case is.",
         "   So platform and phenotype are one axis, and the clean test (hold the platform",
         "   fixed, compare cases to controls) cannot be run at all. Down-sampling removed the",
         "   depth difference; it did not remove the 100 bp vs 151 bp reads, or the chemistry.",
         "   Only an external replication, where platform does not track phenotype, closes this.",
         "",
         "-" * W,
         " 2. REFERENCE",
         "-" * W,
         f"   {'subset':<26}{'cases':>7}{'controls':>10}{'OR (95% CI)':>22}{'p (non-SPA)':>14}",
         f"   {'all':<26}{allr['n_case']:>7}{allr['n_ctrl']:>10}"
         f"{f'{o:.2f} ({lo:.2f}-{hi:.2f})':>22}{allr['p']:>14.3E}",
         "",
         "-" * W,
         " 3. EACH PLATFORM ON ITS OWN  —  independent, so heterogeneity is testable",
         "-" * W,
         "   These share no case with one another. An artefact of one chemistry has to live",
         "   in that chemistry; if unrelated chemistries agree, one artefact would have to be",
         "   in all of them at once, at the same size.",
         "",
         f"   {'subset':<26}{'cases':>7}{'readlen':>9}{'depth':>7}{'OR (95% CI)':>22}"
         f"{'p (non-SPA)':>14}"]
    for tag in sorted(only, key=lambda t: -only[t]["n_case"]):
        r = only[tag]
        o_, lo_, hi_ = orci(r["beta"], r["se"])
        d, rr = rl(tag)
        L.append(f"   {tag:<26}{r['n_case']:>7}{rr:>9.0f}{d:>7.1f}"
                 f"{f'{o_:.2f} ({lo_:.2f}-{hi_:.2f})':>22}{r['p']:>14.3E}")
    if het:
        o_, lo_, hi_ = orci(het["beta"], het["se"])
        L += ["",
              f"   Cochran's Q = {het['q']:.2f} on {het['df']} df,  p = {het['p']:.3f}",
              f"   I-squared   = {het['i2']:.1f}%",
              f"   fixed-effect meta of the {len(only)} platforms: "
              f"OR {o_:.2f} ({lo_:.2f}-{hi_:.2f})",
              ""]
        if het["p"] > 0.05:
            L += ["   No detectable heterogeneity: the platforms are consistent with one common",
                  "   effect. That is what an artefact confined to one chemistry would NOT look",
                  "   like. It does not exclude an artefact shared by all of them."]
        else:
            L += ["   HETEROGENEITY DETECTED. The platforms do not agree on one effect, which is",
                  "   what a platform-specific artefact looks like. Read the per-platform rows",
                  "   before reading the pooled estimate."]
        L += ["",
              "   Small platforms carry wide intervals — that is imprecision, not absence. A",
              "   CI spanning 1 with 27 cases says nothing either way; the heterogeneity test",
              "   is what weighs them properly."]

    L += ["",
          "-" * W,
          " 4. DROPPING ONE PLATFORM  —  nested, so NO heterogeneity test",
          "-" * W,
          "   These overlap each other heavily: any two share most of their samples. A",
          "   heterogeneity test over them would return ~0% whatever the truth, because the",
          "   estimates are forced to agree by construction. None is computed.",
          "",
          "   Read the 'kept' column before the OR. A subset that keeps ~90% of the cases",
          "   MUST land on the full-cohort estimate; only a subset that drops a large share",
          "   of them can move, and only that one is evidence.",
          "",
          f"   {'subset':<26}{'cases kept':>11}{'% of all':>9}{'OR (95% CI)':>22}"
          f"{'p (non-SPA)':>14}{'informative?':>14}"]
    for tag in sorted(minus, key=lambda t: minus[t]["n_case"]):
        r = minus[tag]
        o_, lo_, hi_ = orci(r["beta"], r["se"])
        frac = 100.0 * r["n_case"] / allr["n_case"]
        note = "yes" if frac < 100 * FORCED_FRAC else "no — forced"
        L.append(f"   {tag:<26}{r['n_case']:>11}{frac:>8.0f}%"
                 f"{f'{o_:.2f} ({lo_:.2f}-{hi_:.2f})':>22}{r['p']:>14.3E}{note:>14}")

    inf = [t for t in minus if minus[t]["n_case"] / allr["n_case"] < FORCED_FRAC]
    if inf:
        L += [""]
        for t in inf:
            r = minus[t]
            o_, lo_, hi_ = orci(r["beta"], r["se"])
            pf = pf_of(t)
            L += [f"   {t} keeps {r['n_case']} of {allr['n_case']} cases and still gives "
                  f"OR {o_:.2f} ({lo_:.2f}-{hi_:.2f}).",
                  f"   Those cases are on the other platforms, i.e. a different chemistry and a",
                  f"   different read length from {pf}. An artefact of {pf} cannot be in them."]

    L += ["",
          "-" * W,
          " 5. WHAT THIS DOES AND DOES NOT ESTABLISH",
          "-" * W,
          "   RULES OUT   an artefact confined to one platform — including read length, since",
          "               the only 100 bp platform can be dropped and the effect stays.",
          "   DOES NOT    rule out anything shared by every non-HiSeqX platform. That class is",
          "               narrow — it would have to act on unrelated chemistries at once, at",
          "               the same size — but it is not empty, and this cohort cannot test it.",
          "=" * W, ""]
    Path(args.out_log).write_text("\n".join(L))
    print("\n".join(L))

    # ── figure ──
    fig, ax = plt.subplots(1, 2, figsize=(15, 5.4), facecolor=SURFACE,
                           gridspec_kw=dict(width_ratios=[1, 1]))

    def style(a, title, xlab):
        a.set_facecolor(SURFACE)
        a.set_title(title, color=INK, fontsize=11, loc="left", pad=10)
        a.set_xlabel(xlab, color=MUTED, fontsize=9.5)
        a.grid(axis="x", color=GRID, linewidth=0.8)
        a.set_axisbelow(True)
        a.tick_params(colors=MUTED, labelsize=9)
        a.tick_params(axis="y", length=0, labelcolor=INK)
        for s in ("top", "right", "left"):
            a.spines[s].set_visible(False)
        a.spines["bottom"].set_color(AXIS)

    def or_axis(a, lo, hi):
        """Log x, but labelled 1.5 / 2 / 3 rather than 1.5 x 10^0.

        Matplotlib's default log formatter writes every minor tick in scientific
        notation; over the one decade an odds ratio lives in, that is a dozen
        overlapping labels saying nothing. Fixed ticks at the values a reader
        actually looks for, minors off.
        """
        a.set_xscale("log")
        cand = [0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10]
        ticks = [t for t in cand if lo <= t <= hi]
        a.set_xticks(ticks)
        a.set_xticklabels([f"{t:g}" for t in ticks])
        a.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        a.set_xlim(lo, hi)

    def forest(a, tags, src, colour_fn, label_fn):
        """One row per subset. The interval IS the precision.

        An earlier version also scaled the marker area by 1/(CI width), which in
        panel B made the forced rows — the ones that cannot move — the largest dots
        on the page. One unexplained encoding of precision on top of the interval
        that already shows it is noise at best; here it was backwards.
        """
        labs = []
        for i, t in enumerate(tags):
            r = src[t]
            o_, lo_, hi_ = orci(r["beta"], r["se"])
            c = colour_fn(t, r)
            a.plot([lo_, hi_], [i, i], color=c, lw=1.8, alpha=0.9)
            a.scatter([o_], [i], s=95, color=c, zorder=3)
            a.annotate(f"{o_:.2f}", (o_, i), xytext=(0, 10), textcoords="offset points",
                       ha="center", color=INK, fontsize=8.5)
            labs.append(label_fn(t, r))
        return np.arange(len(labs)), labs

    def ref_row(a, yy, o_, lo_, hi_, label):
        a.plot([lo_, hi_], [yy, yy], color=C_ALL, lw=2.2)
        a.scatter([o_], [yy], s=140, color=C_ALL, marker="D", zorder=3)
        a.annotate(f"{o_:.2f}", (o_, yy), xytext=(0, 11), textcoords="offset points",
                   ha="center", color=INK, fontsize=9, fontweight="bold")
        return label

    # One x range for both panels. On separate scales the nested rows of B — narrow
    # because they share their samples, not because they are precise — get stretched
    # across their panel and read as the sharpest evidence on the page.
    ors = [orci(r["beta"], r["se"]) for r in list(only.values()) + list(minus.values())] \
        + [orci(allr["beta"], allr["se"])]
    XLO = min(min(x[1] for x in ors), 1.0) * 0.92
    XHI = max(x[2] for x in ors) * 1.08

    o_all, lo_all, hi_all = orci(allr["beta"], allr["se"])

    def panel(a, tags, src, colour, label_fn, title):
        """Both panels are the same object — subsets of cases against one reference —
        so they are built by one function and can only differ where the data differ."""
        y, labs = forest(a, tags, src, colour, label_fn)
        yy = len(tags)
        labs.append(ref_row(a, yy, o_all, lo_all, hi_all,
                            f"ALL CASES (reference)\n {allr['n_case']} cases"))
        a.axvline(1.0, color=INK, lw=1.2, ls="--")            # no effect
        a.axvline(o_all, color=C_ALL, lw=1, ls=":", alpha=0.5)  # the reference OR
        # The reference is not one of the subsets; a rule says so without a caption.
        a.axhline(yy - 0.5, color=AXIS, lw=0.8)
        a.set_yticks(list(y) + [yy]); a.set_yticklabels(labs, fontsize=8)
        a.set_ylim(-0.7, yy + 0.7); a.invert_yaxis()
        or_axis(a, XLO, XHI)
        style(a, title, "odds ratio (SAIGE, 95% CI, log scale)")

    # The heterogeneity test belongs to A alone: its rows share no case, so Q means
    # something. B's rows are nested and a Q over them would return ~0% by
    # construction, so B carries no such line at all.
    ttl_a = "Each platform on its own"
    if het:
        ttl_a += f"\nCochran's Q p={het['p']:.2f}, I²={het['i2']:.0f}%  →  " \
                 + ("consistent with one effect" if het["p"] > 0.05 else "platforms DISAGREE")
    panel(ax[0], sorted(only, key=lambda t: -only[t]["n_case"]), only,
          lambda t, r: C_ONLY,
          lambda t, r: f"{pf_of(t).replace('_', ' ')}\n"
                       f" {r['n_case']} cases, {rl(t)[1]:.0f} bp, {rl(t)[0]:.1f}x",
          ttl_a)

    panel(ax[1], sorted(minus, key=lambda t: minus[t]["n_case"]), minus,
          lambda t, r: C_MINUS,
          lambda t, r: f"drop {pf_of(t).replace('_', ' ')}\n"
                       f" {r['n_case']} of {allr['n_case']} cases kept",
          "Dropping one platform\n"
          "any two rows share most of their samples; they are not four confirmations")

    fig.suptitle(
        f"Is {args.variant} an artefact of one platform?  Every estimate is SAIGE, "
        "not a crude test\n"
        "No platform carries both cases and controls, so a within-platform comparison is "
        "impossible here — the cohort's limit, not the analysis's.\n"
        "A and B share one x axis: the panels are directly comparable, and B's short "
        "intervals are nestedness, not precision.",
        color=INK, fontsize=11.5, x=0.004, ha="left", va="top", y=0.995)

    # tight_layout reserves the room for the suptitle itself; leave it alone. Two
    # things that look like improvements are not:
    #   * a `rect` — it is applied ON TOP of that reservation, not instead of it, so
    #     it only ever pushes the panels further down.
    #   * measuring the gap and growing the panels into it — a title's height reads
    #     back as less than half its true value (its position is resolved during the
    #     draw, and even a forced draw() reports the same), so the correction
    #     overshoots and lands the titles inside the suptitle.
    # The band of air under a three-line suptitle is the honest cost of three lines.
    fig.tight_layout()
    # The letters come last, from the axes' final positions, which tight_layout set.
    for a_, letter in ((ax[0], "A"), (ax[1], "B")):
        p = a_.get_position()
        fig.text(p.x0 - 0.052, p.y1 + 0.045, letter, fontsize=13, fontweight="bold",
                 color=INK, va="bottom", ha="left")
    fig.savefig(args.out_fig, dpi=200, facecolor=SURFACE)
    plt.close(fig)


if __name__ == "__main__":
    main()
