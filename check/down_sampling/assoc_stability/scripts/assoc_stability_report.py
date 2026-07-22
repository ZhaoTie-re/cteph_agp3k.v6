#!/usr/bin/env python3
"""Does the association survive having the depth difference taken away?

Every control in this cohort was sequenced on HiSeqX and no case was, so platform
and phenotype are one axis: a depth-driven genotyping difference lands exactly
where an association lands, and could manufacture one out of nothing. Each
replicate re-runs the association after every platform was brought down to the
baseline depth and re-called; the baseline arm is the cohort's own genotypes
through this same code.

Three replicates, because the subsampling is random. One draw cannot separate a
robust signal from a lucky one.

EVERYTHING HERE IS POST-GENOTYPE-QC, on both sides. That is the point, not a
detail: genotype QC is the step through which depth reaches an odds ratio — a
shallow genotype is the one QC turns into a no-call — so a stability check run
before QC would be testing the wrong thing.

SNP-based mirrors analysis/assoc_saige and reads p.value.NA, the p WITHOUT the
saddlepoint correction: SPA's own convergence varies run to run, and the genotypes
are what is under test. OR is exp(BETA), the study's own estimand.

Gene-based mirrors analysis/assoc_rvtest on the MODERATE+HIGH impact stratum.

TWO figures, not one. The SNP test rests on 970 variants across 2,193 samples; the
gene test on 3 variants and about a dozen carriers. They are different weights of
evidence, and one frame invites reading the thinner with the confidence of the
thicker.
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

W = 100
BASELINE = "baseline"
SURFACE, MUTED, GRID, INK, AXIS = "#ffffff", "#4d4d4d", "#e6e6e6", "#000000", "#999999"
# Okabe-Ito: distinguishable under every common form of colour blindness, and
# distinguishable in greyscale, which is how a fair number of readers will see it.
C_SIG = "#CC79A7"                 # Okabe-Ito reddish purple: the threshold lines
C_ARM = {"baseline": "#000000",   # the cohort's own calls
         "rep1": "#0072B2",       # blue
         "rep2": "#D55E00",       # vermillion
         "rep3": "#009E73"}       # bluish green
# Must match 07_genotype_concordance and 09_assoc_concordance, which draw the same
# kind of matrix. A reader flipping between them reads colour as meaning; two
# palettes for one kind of plot invents a distinction that is not there. If this
# ever changes, change it in all three.
CMAP = "Oranges"
GT4 = ["RR", "RA", "AA", "NC"]
# Rows drawn in the figure. The cohort's own no-calls are reported as a count in
# each panel's title, not as a fourth row: that row is ~all NC->NC, so it fills a
# quarter of the plate to restate one number, and it is not part of the question
# the matrix asks (given the cohort called it, did the replicate agree?). The log
# and the TSV keep the full 4x4.
GT_ROWS = ["RR", "RA", "AA"]
GW = -math.log10(5e-8)


def col(a):
    return C_ARM.get(a, "#c2452d")


def read_saige(path):
    out = {}
    with open(path) as fh:
        hdr = fh.readline().split()
        idx = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            f = line.split()
            if len(f) < len(hdr):
                continue
            try:
                out[f[idx["MarkerID"]]] = dict(
                    beta=float(f[idx["BETA"]]), se=float(f[idx["SE"]]),
                    p=float(f[idx["p.value.NA"]]), af=float(f[idx["AF_Allele2"]]),
                    miss=float(f[idx["MissingRate"]]),
                    a1=f[idx["Allele1"]], a2=f[idx["Allele2"]])
            except (ValueError, KeyError):
                continue
    return out


def read_rvtest(path):
    out = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(hdr) or "Pvalue" not in idx:
                continue
            try:
                out[f[0]] = dict(p=float(f[idx["Pvalue"]]),
                                 n=f[idx.get("N_INFORMATIVE", 0)],
                                 nvar=f[idx.get("NumVar", 0)])
            except ValueError:
                continue
    return out


def read_gt(path):
    """SAMPLE<TAB>VARIANT<TAB>DOSAGE -> {variant: {sample: 0|1|2|None}}."""
    out = defaultdict(dict)
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        out[f[1]][f[0]] = None if f[2] in ("NA", ".", "") else int(round(float(f[2])))
    return out


def mlog10(p):
    return -math.log10(p) if p and p > 0 else np.nan


def sci(x, sig=1):
    """LaTeX scientific notation for a figure: 4.0e-08 -> $4.0\\times10^{-8}$.

    Matplotlib's default and Python's '%E' both render the exponent as 'E-08',
    which is a programming convention, not a typographic one. Journals set the
    exponent as a superscript, and a figure that mixes '4.0E-08' with '5x10^-8'
    is using two notations for one idea on the same axis.
    """
    if x is None or x != x or x <= 0:
        return "n/a"
    e = int(math.floor(math.log10(abs(x))))
    m = x / (10 ** e)
    if abs(m - 1.0) < 10 ** (-sig) / 2:
        return f"$10^{{{e}}}$"
    return f"${m:.{sig}f}\\times10^{{{e}}}$"


def odds_ratio(beta, se):
    """exp(BETA) and its 95% CI — the number this study is quoted by.

    SAIGE reports BETA on the log-odds scale. Reporting BETA alone leaves the
    reader to convert before they can compare anything against the published
    result; an earlier version of this pipeline did exactly that and had its BETA
    confused with a crude allele-count OR reported elsewhere. Every OR in this
    pipeline is now exp(BETA) from the fitted model.
    """
    return math.exp(beta), math.exp(beta - 1.96 * se), math.exp(beta + 1.96 * se)


def confusion(base_gt, rep_gt, variants):
    """4x4 counts: baseline {RR,RA,AA,NC} x replicate {RR,RA,AA,NC}.

    The NC ROW is why this is 4x4 and not 3x4. The cohort's own genotypes are
    already missing at some samples — 3.3% at the lead variant — and a matrix that
    drops those rows hides how much of the grid was never data to begin with.
    Both axes are post-genotype-QC, so a no-call on either side is QC's verdict,
    not an accident of the file.
    """
    C = np.zeros((4, 4), dtype=int)
    for v in variants:
        b, r = base_gt.get(v, {}), rep_gt.get(v, {})
        for s in set(b) & set(r):
            i = 3 if b[s] is None else b[s]
            j = 3 if r[s] is None else r[s]
            C[i][j] += 1
    return C


def style(ax, title, xlab, ylab, grid_axis="both"):
    ax.set_facecolor(SURFACE)
    ax.set_title(title, color=INK, fontsize=10.5, loc="left", pad=10)
    ax.set_xlabel(xlab, color=MUTED, fontsize=9.5)
    ax.set_ylabel(ylab, color=MUTED, fontsize=9.5)
    ax.grid(axis=grid_axis, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED, labelsize=9)
    for k in ("top", "right"):
        ax.spines[k].set_visible(False)
    for k in ("left", "bottom"):
        ax.spines[k].set_color(AXIS)


def conf_block(L, tag, C, gt_label):
    tot = int(C.sum())
    lost = int(C[:3, 3].sum())
    gained = int(C[3, :3].sum())
    changed = int(tot - np.trace(C)) - lost - gained
    L += [f"   {tag}   {tot:,} {gt_label} compared",
          f"   {'':>11}{'rep RR':>10}{'rep RA':>10}{'rep AA':>10}{'rep NC':>10}"]
    for i, lab in enumerate(GT4):
        pre = f"base {lab}" + (" *" if lab == "NC" else "")
        L.append(f"   {pre:>11}" + "".join(f"{C[i][j]:>10,}" for j in range(4)))
    L += [f"   * the cohort was ALREADY missing here: {C[3].sum():,} "
          f"({100 * C[3].sum() / tot:.2f}% of the grid) — never data, not lost by down-sampling",
          f"     replicate lost a call the cohort had  : {lost:,}",
          f"     replicate made one the cohort lacked  : {gained:,}",
          f"     called differently on both sides      : {changed:,}",
          ""]
    return L


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snp", action="append", default=[], metavar="ARM:PATH")
    ap.add_argument("--gene", action="append", default=[], metavar="ARM:METHOD:PATH")
    ap.add_argument("--snp-gt", action="append", default=[], metavar="ARM:PATH")
    ap.add_argument("--gene-gt", action="append", default=[], metavar="ARM:PATH")
    ap.add_argument("--ac", action="append", default=[], metavar="PATH")
    ap.add_argument("--lead-variant", required=True)
    ap.add_argument("--lead-gene", required=True)
    ap.add_argument("--gene-sig", type=float, default=5.8e-6,
                    help="gene-based significance threshold, e.g. 0.05/n_genes_tested")
    ap.add_argument("--gene-sig-label", default="Bonferroni, 0.05/8,621 genes")
    ap.add_argument("--out-log", required=True)
    ap.add_argument("--out-fig-snp", required=True)
    ap.add_argument("--out-fig-gene", required=True)
    ap.add_argument("--out-tsv-snp", required=True)
    ap.add_argument("--out-tsv-gene", required=True)
    ap.add_argument("--out-tsv-confusion", required=True)
    args = ap.parse_args()

    snp, gene = {}, defaultdict(dict)
    for spec in args.snp:
        a, p = spec.split(":", 1)
        snp[a] = read_saige(p)
    for spec in args.gene:
        a, m, p = spec.split(":", 2)
        gene[m][a] = read_rvtest(p)
    snp_gt, gene_gt = {}, {}
    for spec in args.snp_gt:
        a, p = spec.split(":", 1)
        snp_gt[a] = read_gt(p)
    for spec in args.gene_gt:
        a, p = spec.split(":", 1)
        gene_gt[a] = read_gt(p)
    acs = defaultdict(dict)
    for path in args.ac:
        for line in open(path):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:
                acs[f[1]][f[0]] = (f[2], f[3])

    arms = [BASELINE] + sorted(a for a in snp if a != BASELINE)
    reps = [a for a in arms if a != BASELINE]
    lead = args.lead_variant
    b = snp.get(BASELINE, {}).get(lead)

    conf_snp, conf_gene = {}, {}
    gene_vars = sorted(gene_gt.get(BASELINE, {}))
    for a in reps:
        if BASELINE in snp_gt and a in snp_gt:
            conf_snp[a] = confusion(snp_gt[BASELINE], snp_gt[a], [lead])
        if BASELINE in gene_gt and a in gene_gt:
            conf_gene[a] = confusion(gene_gt[BASELINE], gene_gt[a], gene_vars)

    with open(args.out_tsv_confusion, "w") as fh:
        fh.write("ANALYSIS\tARM\tBASELINE_GT\tREP_GT\tN\tMEANING\n")
        for name, cm in (("snp_lead", conf_snp), ("gene_variants", conf_gene)):
            for a, C in sorted(cm.items()):
                for i, bi in enumerate(GT4):
                    for j, rj in enumerate(GT4):
                        m = ("cohort already missing" if bi == "NC" and rj == "NC" else
                             "cohort already missing; replicate called it" if bi == "NC" else
                             "replicate lost the call" if rj == "NC" else
                             "agree" if bi == rj else "called differently")
                        fh.write(f"{name}\t{a}\t{bi}\t{rj}\t{C[i][j]}\t{m}\n")

    with open(args.out_tsv_snp, "w") as fh:
        fh.write("ARM\tVARIANT\tALLELE1\tALLELE2\tAF_ALLELE2\tMISSING_RATE\tBETA\tSE\t"
                 "OR\tOR_L95\tOR_U95\tP_NONSPA\tNEGLOG10P\tIS_LEAD\n")
        for a in arms:
            for v in sorted(snp.get(a, {})):
                r = snp[a][v]
                o, lo, hi = odds_ratio(r["beta"], r["se"])
                fh.write(f"{a}\t{v}\t{r['a1']}\t{r['a2']}\t{r['af']:.6f}\t{r['miss']:.6f}\t"
                         f"{r['beta']:.6f}\t{r['se']:.6f}\t{o:.4f}\t{lo:.4f}\t{hi:.4f}\t"
                         f"{r['p']:.6E}\t{mlog10(r['p']):.4f}\t"
                         f"{'yes' if v == lead else 'no'}\n")
    with open(args.out_tsv_gene, "w") as fh:
        fh.write("ARM\tMETHOD\tGENE\tN_INFORMATIVE\tNUMVAR\tPVALUE\tNEGLOG10P\tIS_LEAD\n")
        for m in sorted(gene):
            for a in [BASELINE] + sorted(x for x in gene[m] if x != BASELINE):
                for g in sorted(gene[m].get(a, {})):
                    r = gene[m][a][g]
                    fh.write(f"{a}\t{m}\t{g}\t{r['n']}\t{r['nvar']}\t{r['p']:.6E}\t"
                             f"{mlog10(r['p']):.4f}\t{'yes' if g == args.lead_gene else 'no'}\n")

    # ── log ───────────────────────────────────────────────────────────────
    L = ["=" * W,
         " Does the association survive removing the depth difference?",
         "=" * W,
         f" generated : {datetime.now().isoformat(timespec='seconds')}",
         "",
         " Every control was sequenced on HiSeqX and no case was, so platform and phenotype",
         " are one axis: a depth-driven genotyping difference lands exactly where an",
         " association lands. Each replicate re-runs the association after every platform",
         " was brought down to the baseline depth and re-called.",
         "",
         " ALL GENOTYPES BELOW ARE POST-GENOTYPE-QC, on both sides. That is the point, not a",
         " detail: genotype QC is how depth reaches an odds ratio — a shallow genotype is the",
         " one QC turns into a no-call — so a check run before QC would test the wrong thing.",
         ""]
    shared = set()
    ors_r = []
    if b:
        o, lo, hi = odds_ratio(b["beta"], b["se"])
        L += ["-" * W,
              " 1. SNP-BASED  (SAIGE, non-SPA p.value.NA, post-GT-QC)",
              "-" * W,
              f"   Lead variant: {lead}",
              "",
              f"   {'arm':<10}{'A1/A2':<8}{'AF':>7}{'miss%':>7}{'BETA':>9}{'SE':>8}"
              f"{'OR (95% CI)':>22}{'p (non-SPA)':>14}{'vs base':>9}"]
        for a in arms:
            r = snp.get(a, {}).get(lead)
            if not r:
                L.append(f"   {a:<10}(not tested)")
                continue
            o_, lo_, hi_ = odds_ratio(r["beta"], r["se"])
            vs = "—" if a == BASELINE else f"{r['p'] / b['p']:.2f}x"
            L.append(f"   {a:<10}{r['a1'] + '/' + r['a2']:<8}{r['af']:>7.4f}"
                     f"{100 * r['miss']:>7.2f}{r['beta']:>9.4f}{r['se']:>8.4f}"
                     f"{f'{o_:.2f} ({lo_:.2f}-{hi_:.2f})':>22}{r['p']:>14.3E}{vs:>9}")
        ors_r = [odds_ratio(snp[a][lead]["beta"], snp[a][lead]["se"])[0]
                 for a in reps if lead in snp.get(a, {})]
        if ors_r:
            L += ["",
                  f"   OR   baseline {o:.2f} ({lo:.2f}-{hi:.2f})   replicates "
                  f"{min(ors_r):.2f} to {max(ors_r):.2f}",
                  "",
                  "   Read the OR, not the p. An effect that is really depth cannot keep its size",
                  "   once the depth is gone. A p-value moves on genotype noise alone, and at",
                  "   z~5.5 a few percent off the effect is most of an order of magnitude off p.",
                  ""]
        shared = set(snp.get(BASELINE, {}))
        for a in reps:
            shared &= set(snp.get(a, {}))
        if shared:
            sv = sorted(shared)
            bb = np.array([snp[BASELINE][v]["beta"] for v in sv])
            L += [f"   Across all {len(shared):,} variants in the region:",
                  f"   {'arm':<10}{'r(BETA)':>10}{'max |dBETA|':>14}"]
            for a in reps:
                rb = np.array([snp[a][v]["beta"] for v in sv])
                L.append(f"   {a:<10}{np.corrcoef(bb, rb)[0, 1]:>10.5f}"
                         f"{np.max(np.abs(rb - bb)):>14.4f}")
            L.append("")

    if conf_snp:
        L += ["-" * W,
              " 2. WHAT HAPPENED TO THE LEAD VARIANT'S GENOTYPES  (post-GT-QC both sides)",
              "-" * W,
              "   rows = the cohort's own call; cols = the replicate's.",
              ""]
        for a in reps:
            if a in conf_snp:
                L = conf_block(L, a, conf_snp[a], "samples")

    if gene:
        L += ["-" * W,
              " 3. GENE-BASED  (rvtest, MODERATE+HIGH stratum, post-GT-QC)",
              "-" * W,
              f"   Gene: {args.lead_gene}",
              "",
              f"   {'method':<10}{'arm':<10}{'N':>7}{'NumVar':>8}{'Pvalue':>14}"
              f"{'-log10p':>10}{'vs base':>10}"]
        for m in sorted(gene):
            gb = gene[m].get(BASELINE, {}).get(args.lead_gene)
            for a in [BASELINE] + sorted(x for x in gene[m] if x != BASELINE):
                r = gene[m].get(a, {}).get(args.lead_gene)
                if not r:
                    continue
                vs = "—" if a == BASELINE else (f"{r['p'] / gb['p']:.2f}x" if gb else "n/a")
                L.append(f"   {m:<10}{a:<10}{str(r['n']):>7}{str(r['nvar']):>8}"
                         f"{r['p']:>14.3E}{mlog10(r['p']):>10.2f}{vs:>10}")
            L.append("")

    if conf_gene:
        L += ["-" * W,
              " 4. WHAT HAPPENED TO THE GENE'S GENOTYPES  (post-GT-QC both sides)",
              "-" * W,
              f"   Pooled over the {len(gene_vars)} variants the gene test includes.",
              ""]
        for a in reps:
            if a in conf_gene:
                L = conf_block(L, a, conf_gene[a], "genotypes")

    if acs:
        L += ["-" * W,
              " 5. ALLELE COUNTS at the gene's variants",
              "-" * W,
              "   The variant set is held fixed at the reference's, so the arms differ by",
              "   genotypes alone — which makes AC where any change shows up. This gene rests",
              "   on so few carriers that losing one moves the p-value several fold.",
              "",
              f"   {'variant':<26}" + "".join(f"{a:>12}" for a in arms)]
        for v in sorted(acs):
            L.append(f"   {v:<26}" + "".join(f"{acs[v].get(a, ('—',))[0]:>12}" for a in arms))
        L.append("")

    L += ["-" * W,
          " 6. WHAT THIS CANNOT SAY",
          "-" * W,
          "   A replicate differs from the baseline for two reasons at once: the depth removed,",
          "   and the re-calling itself (a 426-sample joint call rather than the cohort-wide",
          "   one). This design bounds their SUM, not either alone — which is the useful",
          "   direction: if the sum is negligible, neither term can carry the signal.",
          "",
          "   Down-sampling levels the DEPTH, not the platform. A T7 read at baseline depth is",
          "   still a 100 bp DNBSEQ read. See 06_platform_report for that question.",
          "=" * W, ""]
    Path(args.out_log).write_text("\n".join(L))
    print("\n".join(L))

    # ─────────────────────────────────────────────────────────────────────
    # Figures. Two, because the SNP test rests on 970 variants across 2,193
    # samples and the gene test on 3 variants and about a dozen carriers: one
    # frame would lend the thinner evidence the confidence of the thicker.
    # ─────────────────────────────────────────────────────────────────────
    def panel_letter(ax, letter):
        """Publication convention: the letter belongs to the panel, outside the axes,
        so it survives being cropped into a multi-panel plate."""
        ax.text(-0.20, 1.06, letter, transform=ax.transAxes, fontsize=12,
                fontweight="bold", color=INK, va="bottom", ha="left")

    def block_head(fig, gs, letter, text, pad_in, FH):
        """The matrix grid is ONE panel, so it gets ONE letter — above the block, not
        on its first cell.

        Every cell of the grid is the same measurement (cohort vs replicate) repeated
        over replicates and variants; lettering them D, E, F… announces them as
        separate pieces of evidence and invites the reader to cite one alone. The
        cells keep their own titles to say which replicate and which variant they are.
        """
        y = gs.top + pad_in / FH
        fig.text(gs.left - 0.028, y, letter, fontsize=12, fontweight="bold",
                 color=INK, va="bottom", ha="left")
        fig.text(gs.left, y, text, fontsize=9.5, color=INK, va="bottom", ha="left")

    def conf_grid(fig, gs, entries):
        """One confusion matrix per cell + a single shared colourbar.

        Per variant, never pooled: pooling hides that these variants carry different
        numbers of carriers, and a matrix summed over them reads as one well-sampled
        variant rather than three thin ones.
        """
        im = None
        for k, (r0, c0, C, title, xlab, ylab) in enumerate(entries):
            ax = fig.add_subplot(gs[r0, c0])
            D = C[:3, :]                      # drop the cohort's own no-call row
            rowsum = D.sum(axis=1, keepdims=True)
            frac = np.divide(D, rowsum, out=np.zeros_like(D, dtype=float),
                             where=rowsum != 0)
            im = ax.imshow(frac, cmap=CMAP, vmin=0, vmax=1, aspect="auto")
            for i in range(len(GT_ROWS)):
                for j in range(4):
                    n = D[i][j]
                    ax.text(j, i, f"{n:,}" if n else "·", ha="center", va="center",
                            fontsize=8, color=INK if frac[i][j] < 0.6 else SURFACE,
                            fontweight="bold" if i == j and n else "normal")
            ax.set_xticks(range(4)); ax.set_xticklabels(GT4, fontsize=8)
            ax.set_yticks(range(len(GT_ROWS))); ax.set_yticklabels(GT_ROWS, fontsize=8)
            if xlab:
                ax.set_xlabel(xlab, color=MUTED, fontsize=8.5)
            if ylab:
                ax.set_ylabel(ylab, color=MUTED, fontsize=8.5)
            ax.set_title(title, color=INK, fontsize=9, loc="left", pad=5)
            # The line separates "called differently" from "not called at all".
            ax.axvline(2.5, color=INK, lw=1.2, ls="--")
            ax.tick_params(colors=MUTED, length=0)
            for sp in ax.spines.values():
                sp.set_edgecolor(AXIS); sp.set_linewidth(0.8)
        return im

    # ══════════════════════ FIGURE 1: SNP ══════════════════════
    # Top row taller than the matrices: A/B/C carry axes, ticks and a legend, the
    # matrices carry a 3x4 grid of numbers. Equal heights would starve the first
    # and pad the second.
    # The top row is placed explicitly and the gridspec covers only the matrices.
    # Sharing one grid between them puts the top row's tick labels, axis titles and
    # legend inside the cell below it — which is how B's axis label came to sit on
    # E's title. Two bands, sized for what each holds.
    FH = 7.8
    fig = plt.figure(figsize=(12.4, FH), facecolor=SURFACE)
    TOP_H = 2.1 / FH
    TOP_Y = 1.0 - (0.85 + 0.42 + 2.1) / FH
    gs = fig.add_gridspec(1, 3, wspace=0.30, left=0.075, right=0.885,
                          top=TOP_Y - 1.0 / FH, bottom=0.09)

    a = fig.add_axes([0.075, TOP_Y, 0.215, TOP_H])
    if b:
        ob, lob, hib = odds_ratio(b["beta"], b["se"])
        a.axhspan(lob, hib, color=C_ARM[BASELINE], alpha=0.07, zorder=0)
        a.axhline(ob, color=C_ARM[BASELINE], lw=1.1, ls="--", zorder=1)
    for i, arm in enumerate(arms):
        r = snp.get(arm, {}).get(lead)
        if not r:
            continue
        o_, lo_, hi_ = odds_ratio(r["beta"], r["se"])
        a.errorbar([i], [o_], yerr=[[o_ - lo_], [hi_ - o_]], fmt="o", ms=8,
                   capsize=5, lw=1.6, color=col(arm), ecolor=col(arm), zorder=3)
        a.annotate(f"{o_:.2f}", (i, o_), xytext=(10, 0), textcoords="offset points",
                   color=INK, fontsize=8.5, fontweight="bold", va="center")
    a.axhline(1.0, color=MUTED, lw=0.9, ls=":")
    a.set_xticks(range(len(arms))); a.set_xticklabels(arms, fontsize=8.5)
    a.set_xlim(-0.5, len(arms) - 0.35)
    style(a, f"Effect size at {lead}", "", "odds ratio (95% CI)")
    panel_letter(a, "A")

    bx = fig.add_axes([0.375, TOP_Y, 0.215, TOP_H])
    for i, arm in enumerate(arms):
        r = snp.get(arm, {}).get(lead)
        if not r:
            continue
        bx.bar([i], [mlog10(r["p"])], 0.55, color=col(arm))
        # The bars land close to the genome-wide line, so a label sits on it more
        # often than not; the pad lets the number win rather than the rule.
        bx.annotate(sci(r["p"]), (i, mlog10(r["p"])), xytext=(0, 3),
                    textcoords="offset points", ha="center", color=INK, fontsize=7.5,
                    zorder=4, bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.8))
    bx.axhline(GW, color=C_SIG, lw=1.2, ls="--")
    ys = [mlog10(snp[a2][lead]["p"]) for a2 in arms if lead in snp.get(a2, {})]
    if ys:
        bx.set_ylim(0, max(ys + [GW]) * 1.22)
    bx.set_xticks(range(len(arms))); bx.set_xticklabels(arms, fontsize=8.5)
    # The threshold is named in the title, not on the line: a bar stands at every x,
    # so an in-plot label lands on one of them or on the p-value above it.
    style(bx, f"Significance\ndashed: genome-wide {sci(5e-8)}", "",
          "$-\\log_{10}$ P (non-SPA)")
    panel_letter(bx, "B")

    # C belongs to the top band with A and B, not to the matrix gridspec below —
    # gs's only row is the matrices' row, so gs[0, 2] would put the scatter under F.
    c = fig.add_axes([0.675, TOP_Y, 0.215, TOP_H])
    if shared:
        sv = sorted(shared)
        bp = np.array([mlog10(snp[BASELINE][v]["p"]) for v in sv])
        for arm in reps:
            rp = np.array([mlog10(snp[arm][v]["p"]) for v in sv])
            ok = np.isfinite(bp) & np.isfinite(rp)
            c.scatter(bp, rp, s=9, alpha=0.4, linewidths=0, color=col(arm),
                      label=f"{arm}  $r$={np.corrcoef(bp[ok], rp[ok])[0, 1]:.4f}")
        lim = [0, float(np.nanmax(bp)) * 1.1]
        c.plot(lim, lim, color=INK, lw=1.1, ls="--", zorder=1)
        c.set_xlim(lim); c.set_ylim(lim)
        c.legend(fontsize=7.5, frameon=False, labelcolor=INK, loc="upper left")
    style(c, f"All {len(shared):,} variants in the region",
          "cohort  $-\\log_{10}$ P", "replicate  $-\\log_{10}$ P")
    panel_letter(c, "C")

    ent = []
    for k, arm in enumerate(reps[:3]):
        C = conf_snp.get(arm)
        if C is None:
            continue
        nc = int(C[3].sum()); tot = int(C.sum())
        # The cohort's own missingness stays as a number even though its row is not
        # drawn — it is the denominator the reader needs to judge the rest.
        ent.append((0, k, C,
                    f"{arm}   (cohort itself missing {nc:,}/{tot:,}, {100*nc/tot:.1f}%)",
                    f"{arm}, post-GT-QC", "cohort, post-GT-QC" if k == 0 else ""))
    im = conf_grid(fig, gs, ent)
    block_head(fig, gs, "D",
               f"Genotype concordance at {lead}, cohort vs each replicate",
               0.22, FH)
    if im is not None:
        cax = fig.add_axes([0.905, 0.13, 0.011, 0.22])
        cb = fig.colorbar(im, cax=cax)
        cb.set_label("fraction of the cohort's row", color=MUTED, fontsize=8)
        cb.ax.tick_params(colors=MUTED, labelsize=7.5)

    fig.suptitle(
        f"SNP-based association at {lead} after levelling sequencing depth\n"
        "All genotypes post-genotype-QC. Black, the cohort's own calls; coloured, re-called "
        "after every platform was brought to the baseline depth.\n"
        "D scores only the samples the cohort itself called; its own missing genotypes are "
        "given as a count per matrix, since they were never data to agree about.",
        color=INK, fontsize=10, x=0.008, ha="left", va="top", y=0.995)
    fig.savefig(args.out_fig_snp, dpi=200, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    # ══════════════════════ FIGURE 2: GENE ══════════════════════
    nv = len(gene_vars)
    # Same two-band layout as the SNP figure, for the same reason.
    # Budgeted in inches, then converted: a 3-line suptitle and a 2-line panel
    # title are a fixed physical height, while the figure grows with the number of
    # variants. Fractions alone would let them collide as soon as nv changes.
    nrow = max(nv, 1)
    FH = 3.4 + 1.75 * nrow
    fig = plt.figure(figsize=(12.4, FH), facecolor=SURFACE)
    TOP_H = 1.35 / FH                        # the A/B axes themselves
    TOP_Y = 1.0 - (1.10 + 0.42 + 1.35) / FH  # suptitle (4 lines) + panel title + axes
    # 1.15in of clearance: B's axis label (0.3), A's legend below its ticks (0.4),
    # and the matrices' own two-line titles (0.45). Each was measured by watching
    # which one landed on which.
    gs = fig.add_gridspec(nrow, 3, hspace=0.62, wspace=0.30,
                          left=0.075, right=0.885,
                          top=TOP_Y - 1.15 / FH, bottom=0.05)

    # Two panels across the top, not three. An earlier version drew "allele count"
    # and "carriers per variant" side by side before noticing they are the same
    # number: every carrier here is heterozygous, so AC == carrier count. Two
    # panels of one quantity is not two pieces of evidence.
    #
    # Both get explicit axes rather than gridspec cells. The grid below is three
    # columns wide because there are three replicates; the top row has two panels,
    # and forcing them into that grid leaves one cramped and the other short of the
    # left margin its full CHROM:POS:REF:ALT labels need.
    a = fig.add_axes([0.075, TOP_Y, 0.32, TOP_H])
    methods = sorted(gene)
    if methods:
        x = np.arange(len(methods))
        wd = 0.8 / max(len(arms), 1)
        for j, arm in enumerate(arms):
            ys = [mlog10(gene[m].get(arm, {}).get(args.lead_gene, {}).get("p", np.nan))
                  if gene[m].get(arm, {}).get(args.lead_gene) else np.nan for m in methods]
            a.bar(x + j * wd - 0.4 + wd / 2, ys, wd * 0.9, label=arm, color=col(arm))
        gsig = -math.log10(args.gene_sig)
        a.axhline(gsig, color=C_SIG, lw=1.2, ls="--", zorder=1)
        allp = [mlog10(gene[m][ar][args.lead_gene]["p"])
                for m in methods for ar in arms
                if gene[m].get(ar, {}).get(args.lead_gene)]
        a.set_xticks(x); a.set_xticklabels(methods, fontsize=8.5)
        if allp:
            a.set_ylim(0, max(allp + [gsig]) * 1.22)
        a.set_xlim(-0.5, len(methods) - 0.5)
        # Below the tick labels, not on them: -0.09 puts the legend row on the
        # method names.
        a.legend(fontsize=7.5, frameon=False, labelcolor=INK, ncol=4,
                 loc="upper center", bbox_to_anchor=(0.5, -0.22), columnspacing=1.2)
    # The threshold's value goes in the title, not next to its line: the bars span
    # the full width at every x, so an in-plot label lands on one of them.
    style(a, f"{args.lead_gene}, MODERATE+HIGH stratum\n"
             f"dashed: {args.gene_sig_label.split(',')[0]} "
             f"{sci(args.gene_sig)}, {args.gene_sig_label.split(',', 1)[1].strip()}",
          "", "$-\\log_{10}$ P")
    panel_letter(a, "A")

    # Horizontal bars: the labels are full CHROM:POS:REF:ALT IDs and only fit on a
    # vertical axis. A bare coordinate is not the variant — two alleles can share a
    # position — and every other table here is keyed on the ID, so an axis that
    # drops the alleles cannot be joined back to them.
    bx = fig.add_axes([0.60, TOP_Y, 0.285, TOP_H])
    if acs:
        vs = sorted(acs)
        y = np.arange(len(vs))
        ht = 0.8 / max(len(arms), 1)
        for j, arm in enumerate(arms):
            xsv = [float(acs[v][arm][0]) if acs[v].get(arm) else np.nan for v in vs]
            bx.barh(y + j * ht - 0.4 + ht / 2, xsv, ht * 0.9, label=arm, color=col(arm))
            for i, val in enumerate(xsv):
                if val == val:
                    bx.annotate(f"{val:.0f}", (val, y[i] + j * ht - 0.4 + ht / 2),
                                xytext=(3, 0), textcoords="offset points", va="center",
                                color=MUTED, fontsize=6.5)
        allv = [float(acs[v][ar][0]) for v in vs for ar in arms if acs[v].get(ar)]
        bx.set_xlim(0, max(allv) * 1.18)
        bx.set_yticks(y); bx.set_yticklabels(vs, fontsize=7.5)
        bx.set_ylim(len(vs) - 0.5, -0.5)
    style(bx, "Allele count per variant\nevery carrier is heterozygous, so this is also "
              "the carrier count", "allele count", "", grid_axis="x")
    bx.text(-0.40, 1.06, "B", transform=bx.transAxes, fontsize=12, fontweight="bold",
            color=INK, va="bottom", ha="left")

    ent = []
    for vi, v in enumerate(gene_vars):
        for k, arm in enumerate(reps[:3]):
            if arm not in gene_gt or BASELINE not in gene_gt:
                continue
            Cv = confusion(gene_gt[BASELINE], gene_gt[arm], [v])
            nc = int(Cv[3].sum()); tot = int(Cv.sum())
            ent.append((vi, k, Cv,
                        f"{v}\n{arm}   (cohort itself missing {nc:,}/{tot:,})",
                        f"{arm}, post-GT-QC" if vi == nv - 1 else "",
                        "cohort, post-GT-QC" if k == 0 else ""))
    im = conf_grid(fig, gs, ent)
    block_head(fig, gs, "C",
               "Genotype concordance per variant, cohort vs each replicate "
               "(rows: variants; columns: replicates)",
               0.38, FH)
    if im is not None:
        cax = fig.add_axes([0.905, 0.09, 0.011, 0.26])
        cb = fig.colorbar(im, cax=cax)
        cb.set_label("fraction of the cohort's row", color=MUTED, fontsize=8)
        cb.ax.tick_params(colors=MUTED, labelsize=7.5)

    # Hard-wrapped at ~130 characters. The suptitle is drawn from x=0.008 with no
    # wrapping of its own, so a longer line runs past the right edge and
    # bbox_inches="tight" then widens the saved canvas to hold it — which is what
    # left the matrices looking narrow inside a band of whitespace.
    fig.suptitle(
        f"Gene-based association at {args.lead_gene} after levelling sequencing depth\n"
        "All genotypes post-genotype-QC. One matrix per variant, never pooled: these variants "
        "differ in both carrier count and in the\n"
        "cohort's own missingness, and pooling would read as one well-sampled variant rather "
        "than three thin ones.\n"
        "This test rests on the carriers counted in B — a dozen of them. Read it beside the "
        "SNP figure, not as its equal.",
        color=INK, fontsize=10, x=0.008, ha="left", va="top", y=0.995)
    fig.savefig(args.out_fig_gene, dpi=200, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


if __name__ == "__main__":
    main()
