#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The external check, drawn: our control allele frequencies against the
#           published Japanese ones. TWO figures, from one pass over the tables.
#
#           typing_qc.png answers "was a call made, and how deeply did it
#           resolve". It cannot answer "is the call right", and neither can any
#           other output of this component. This is the closest substitute -- if
#           the typing were systematically mis-assigning alleles the frequency
#           spectrum would move off the identity line, and no call-rate metric
#           would notice.
#
#           WHY IT IS TWO FIGURES.
#             allele_frequency.png       the answer: pooled agreement, per-locus
#                                        concordance, per-locus unconfirmed share
#             allele_frequency_loci.png  the record: every locus, one panel each
#
#           NO CAPTION IS RENDERED INTO EITHER PNG. Every word of explanation
#           lives in report/hla_typing_report.qmd, written for a reader who has
#           never heard of HLA. The caption block was taking 42 % of the summary
#           canvas -- more than the plot it described.
#
#           THERE IS NO SPEARMAN HERE ANY MORE. See allele_freq_check.py: 70-96 %
#           of a locus's vector is an allele one source carries and the other does
#           not, so a rank correlation over it is dominated by ties at zero.
#           n_shared is drawn instead, which is the number that makes r readable.
# Component: hla.typing
# ---------------------------------------------------------------------------
import argparse
import string
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'analysis' / '_shared' / 'scripts'))
import plot_style as S                   # noqa: E402
import figure_doc                        # noqa: E402

ROW_IN = 0.190        # inches per locus row in the two bar panels
LOCI_NCOL = 5         # the record grid. 4 columns leaves THREE empty cells in the
                      # last row of a 13-locus figure -- 5.2 x 2.0 in of nothing,
                      # which reads as three missing panels. 5 columns leaves two
                      # at the end of a row, which reads as "that is all of them".
LOCI_PANEL_IN = 1.42  # per grid ROW of the record
FS_SMALL = 7.0        # every in-axes annotation in these figures
FS_AXIS = 10.0        # every axis label. The 14 pt 'slide' default is wider than
                      # the narrow panels here, so labels ran off the canvas.
LABEL_PAD_PT = 3.5
# Three labels, not six. An HLA frequency table is a handful of common alleles and a
# long tail at the origin, so naming more than the top few piles the text on top of
# itself exactly where the points are densest.
N_LABEL = 3
MIN_DIFF_TO_NAME = 0.02
N_OUTLIER_MAX = 2
MIN_CHR = 20          # control chromosomes below which nothing is drawn
R_GOOD = 0.98         # the concordance a locus has to clear to be read at face value


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--check', required=True, help='allele_frequency_check.tsv')
    p.add_argument('--summary', required=True, help='allele_frequency_summary.tsv')
    p.add_argument('--summary-secondary', default=None,
                   help='the SECOND panel\'s allele_frequency_summary.tsv. Panel (d) is '
                        'the reason both panels are run: it shows that they agree on '
                        'concordance, and on how much evidence each reaches that on.')
    p.add_argument('--reference-name-secondary', default='1000 Genomes JPT')
    p.add_argument('--reference-name', default='jMorp 61KJPN-HLA')
    p.add_argument('--flag-genes', default='DRB3',
                   help='loci to mark with a dagger on the figure. This is a LABEL; '
                        'whether a locus is comparable at all is decided from the two '
                        'denominators in allele_freq_check.py and read from the '
                        '`comparable` column')
    p.add_argument('--out-png', required=True, help='the summary figure')
    p.add_argument('--out-png-loci', required=True, help='the per-locus record')
    return p.parse_args()


def panel_head(ax, letter, title, *, pad=6.0, size=None):
    """Bold panel letter and a short title, as ONE left-aligned artist.

    plot_style.panel_tag() writes the letter at loc='left' and the title at
    loc='center', and its own docstring says not to pass a title below about
    2.2 in wide because the two collide. Most panels in this component are
    narrower than that. Written as one left-aligned string they cannot collide,
    and it reads as a label rather than as two competing ones.

    The title says WHAT the panel shows. It never says what to conclude -- that
    is the document's job, which is why no caption is rendered into these PNGs.
    """
    ax.set_title(rf'$\bf{{{letter}}}$  {title}', loc='left', pad=pad,
                 fontsize=size or (plt.rcParams['axes.titlesize'] - 2))


TITLE_H = 0.34        # inches reserved for the figure title


def lay_out(fig, axes, *, plot_h, top_pad, wspace, hspace, right=0.975, title=None):
    """Margins, canvas and figure title, with no caption block.

    caption_block() is the only thing in the shared toolkit that sets the canvas
    height and the measured left margin, so dropping it means doing that here.
    A transcription of what it does, minus the caption text -- but WITH a title:
    a figure lifted out of the report still has to say what it is, and each panel
    says what it shows. What none of them say is what to conclude; that is the
    document's job, and it is why the caption block is gone.
    """
    fig.set_layout_engine('none')
    xlab_in = (plt.rcParams['axes.labelsize']
               + plt.rcParams['xtick.labelsize'] + 8) * 1.5 / 72.0
    head = TITLE_H if title else 0.0
    fh = plot_h + head                       # the panels keep the size they were given
    fig.set_size_inches(S.COL_DOUBLE, fh, forward=True)
    kw = dict(right=right, top=1.0 - (top_pad + head) / fh,
              bottom=xlab_in / fh, wspace=wspace, hspace=hspace)
    left = S.fit_left_margin(fig, axes)
    fig.subplots_adjust(left=left, **kw)
    need = S.fit_left_margin(fig, axes)
    if need > left + 0.002:
        fig.subplots_adjust(left=need, **kw)
    if title:
        # WRAPPED. An unwrapped suptitle silently runs off both edges of the
        # canvas -- matplotlib never clips it and never warns. The character
        # estimate is deliberately conservative; a title that wraps to two lines
        # is fine, a title bleeding into the margin is not.
        fs = plt.rcParams['axes.titlesize'] + 1.5
        per_line = max(20, int(S.COL_DOUBLE * 72.0 * 0.88 / (fs * 0.56)))
        wrapped = textwrap.fill(title, per_line)
        n_lines = wrapped.count('\n') + 1
        if n_lines > 1:                       # give the extra line its own room
            fh += (n_lines - 1) * fs * 1.25 / 72.0
            fig.set_size_inches(S.COL_DOUBLE, fh, forward=True)
            fig.subplots_adjust(top=1.0 - (top_pad + head
                                           + (n_lines - 1) * fs * 1.25 / 72.0) / fh,
                                bottom=xlab_in / fh)
        fig.suptitle(wrapped, y=1.0 - 0.09 / fh, va='top', ha='center',
                     fontsize=fs, fontweight='bold', color=S.INK)


MARGIN_PX = 1.5  # slack between an estimated label box and the axes edge


def place_point_labels(ax, xs, ys, texts, *, marker_pt=3.0):
    """Label each point on whichever side is clearest of the OTHER points.

    Twenty-four candidates per point -- eight directions at three distances --
    scored against the BOXES of the other markers and of the labels already
    placed, with leaving the axes weighted far above any overlap and a final
    slide back inside so no label can be clipped. The first version offset every
    label to the upper right by a fixed amount and separated them in y only,
    which painted several across their own markers (`04:05P` rendered as
    `04:@5P`); the second scored eight directions at one distance against marker
    CENTRES, which still let `03:03P` land on `06:01P` and pushed `02:02P` out of
    its panel because a crowded origin left no good direction at that one radius.
    """
    fig = ax.figure
    fig.canvas.draw()
    pts = ax.transData.transform(np.column_stack([xs, ys]))
    r = (marker_pt + LABEL_PAD_PT) * fig.dpi / 72.0
    box = ax.get_window_extent()
    px2 = fig.dpi / 72.0

    def overlap(a, b):
        """Area two boxes share, 0 if they do not touch."""
        return (max(0.0, min(a[2], b[2]) - max(a[0], b[0]))
                * max(0.0, min(a[3], b[3]) - max(a[1], b[1])))

    # Eight directions at three distances. One radius is not enough: where the
    # points crowd the origin every direction at that radius is occupied, and the
    # best of eight bad candidates is still bad. Stepping out gives the crowded
    # panels somewhere to go without detaching the label from its point, which
    # the tie-break below keeps it from doing gratuitously.
    dirs = []
    for mult in (1.0, 1.7, 2.5):
        for k in range(8):
            a = k * np.pi / 4.0
            dx, dy = r * mult * np.cos(a), r * mult * np.sin(a)
            ha = 'left' if dx > 1e-9 else ('right' if dx < -1e-9 else 'center')
            va = 'bottom' if dy > 1e-9 else ('top' if dy < -1e-9 else 'center')
            dirs.append((dx, dy, ha, va, mult))

    mboxes = [(x - marker_pt * px2, y - marker_pt * px2,
               x + marker_pt * px2, y + marker_pt * px2) for x, y in pts]
    placed, out = [], []
    for i, (ppx, ppy) in enumerate(pts):
        # A margin, because the width is estimated from the character count while
        # the geometry check in verify.sh measures what was actually rendered.
        w = len(texts[i]) * FS_SMALL * 0.62 * px2 + 2 * MARGIN_PX
        h = FS_SMALL * 1.25 * px2 + 2 * MARGIN_PX
        best, best_score = None, -1e18
        for dx, dy, ha, va, mult in dirs:
            cx, cy = ppx + dx, ppy + dy
            x0 = cx if ha == 'left' else (cx - w if ha == 'right' else cx - w / 2)
            y0 = cy if va == 'bottom' else (cy - h if va == 'top' else cy - h / 2)
            cand = (x0, y0, x0 + w, y0 + h)
            outside = (max(0.0, box.x0 - x0) + max(0.0, cand[2] - box.x1)
                       + max(0.0, box.y0 - y0) + max(0.0, cand[3] - box.y1))
            hit = sum(overlap(cand, m) for j, m in enumerate(mboxes) if j != i)
            hit += sum(overlap(cand, q) for q in placed)
            # own marker last: a label may not sit on the point it names either
            hit += overlap(cand, mboxes[i])
            # Leaving the axes outweighs any overlap: a clipped label is drawn
            # into the neighbouring panel, which is worse than being crowded.
            # The last term is the tie-break that keeps a label near its point.
            score = (-(hit / (w * h)) * 1000.0 - outside * 5000.0
                     - (mult - 1.0) * 60.0)
            if score > best_score:
                best, best_score = (dx, dy, ha, va, cand), score
        dx, dy, ha, va, cand = best
        # Even the best candidate can cross an edge where the axes are narrower
        # than the label is wide. Slide it back in rather than let it be clipped.
        sx = max(0.0, box.x0 - cand[0]) - max(0.0, cand[2] - box.x1)
        sy = max(0.0, box.y0 - cand[1]) - max(0.0, cand[3] - box.y1)
        dx, dy = dx + sx, dy + sy
        cand = (cand[0] + sx, cand[1] + sy, cand[2] + sx, cand[3] + sy)
        placed.append(cand)
        out.append(ax.annotate(texts[i], (xs[i], ys[i]), textcoords='offset points',
                               xytext=(dx / px2, dy / px2),
                               ha=ha, va=va, fontsize=FS_SMALL, color=S.INK_SOFT))
    return out


def to_label(rows, n_label=N_LABEL, n_out=N_OUTLIER_MAX):
    """Which alleles to name: the commonest in the reference, plus real outliers."""
    if not len(rows):
        return set()
    named = set(rows.nlargest(min(n_label, len(rows)), 'freq_ref').allele)
    d = pd.to_numeric(rows.diff_ctrl_ref, errors='coerce').abs()
    out = rows.assign(_d=d)[d >= MIN_DIFF_TO_NAME].nlargest(n_out, '_d')
    return named | set(out.allele)


def unconfirmed_by_locus(chk):
    """Control chromosomes on an allele the reference does not carry, per locus.

    Derived from in_obs_only rather than read from a second table: pooled over
    loci this reproduces allele_frequency_sample.tsv exactly, so the two views
    cannot drift apart.
    """
    g = chk.groupby('gene').apply(
        lambda d: pd.Series({'n_chr': d.n_chr_ctrl.max(),
                             'unconf': d.loc[d.in_obs_only == 1, 'count_ctrl'].sum()}),
        include_groups=False)
    g['share'] = g.unconf / g.n_chr.replace(0, np.nan)
    return g


def scatter(ax, rows, *, s_pt=15, n_label=N_LABEL, n_out=N_OUTLIER_MAX, keep_gene=False):
    """Observed control frequency against the reference, with the identity line."""
    x, y = rows.freq_ref.to_numpy(float), rows.freq_ctrl.to_numpy(float)
    # 1.06, not 1.22: the old headroom left the upper 60 % of a square panel empty
    # because it was sized for labels that are now placed by measurement.
    hi = float(max(np.nanmax(x), np.nanmax(y))) * 1.06 + 0.005
    ax.plot([0, hi], [0, hi], color=S.REFERENCE, lw=1.0, ls='--', zorder=1)
    named = to_label(rows, n_label, n_out)
    is_out = rows.allele.isin(named).to_numpy()
    err = np.vstack([y - rows.ci_lo_ctrl.to_numpy(float),
                     rows.ci_hi_ctrl.to_numpy(float) - y])
    ax.errorbar(x, y, yerr=np.clip(err, 0, None), fmt='none',
                ecolor=S.NEUTRAL_D, elinewidth=0.7, zorder=2)
    ax.scatter(x[~is_out], y[~is_out], s=s_pt, color=S.DATA, alpha=0.75,
               edgecolors='none', zorder=3)
    ax.scatter(x[is_out], y[is_out], s=s_pt + 4, color=S.ACCENT, alpha=0.9,
               edgecolors='none', zorder=4)
    # Limits BEFORE placing labels: placement measures display coordinates.
    ax.set_xlim(0, hi); ax.set_ylim(0, hi)
    labels = [al if keep_gene else al.split('*', 1)[-1] for al in rows.allele[is_out]]
    if len(labels):
        place_point_labels(ax, x[is_out], y[is_out], labels,
                           marker_pt=np.sqrt((s_pt + 4) / np.pi))
    return hi


def locus_panel(ax, rows, gene, letter, *, xlab=False, ylab=False, flag=False):
    """One locus of the record figure."""
    n_chr = int(rows.n_chr_ctrl.max()) if len(rows) else 0
    if n_chr < MIN_CHR:
        ax.text(0.5, 0.5, 'too few control\nchromosomes', ha='center', va='center',
                transform=ax.transAxes, color=S.INK_SOFT, fontsize=FS_SMALL)
        ax.set_xticks([]); ax.set_yticks([])
        panel_head(ax, letter, f'HLA-{gene}' + (' †' if flag else ''),
                   size=FS_SMALL + 2)
        return
    scatter(ax, rows, s_pt=9, n_label=2, n_out=1)
    ax.tick_params(labelsize=FS_SMALL)
    ax.locator_params(nbins=3)
    if xlab:
        ax.set_xlabel('reference frequency', fontsize=FS_SMALL + 2)
    if ylab:
        ax.set_ylabel('observed frequency', fontsize=FS_SMALL + 2)
    panel_head(ax, letter, f'HLA-{gene}' + (' †' if flag else ''), size=FS_SMALL + 2)
    S.despine(ax)


def main():
    args = parse_args()
    S.setup_style('slide')

    chk = pd.read_csv(args.check, sep='\t')
    summ = pd.read_csv(args.summary, sep='\t')
    genes = list(summ.gene) if len(summ) else sorted(chk.gene.unique())
    if not genes:
        raise SystemExit(f'ABORT: {args.check} names no gene')
    flag_set = {g.strip() for g in args.flag_genes.split(',') if g.strip()} & set(genes)
    flagged = ' and '.join(f'HLA-{g}' for g in sorted(flag_set))

    n_ctrl = int(chk.n_chr_ctrl.max()) if len(chk) else 0
    n_ref = int(chk.n_chr_ref.max()) if len(chk) else 0
    unc = unconfirmed_by_locus(chk)
    have = (summ[pd.to_numeric(summ.pearson_r, errors='coerce').notna()
                 & (pd.to_numeric(summ.n_chr_ctrl, errors='coerce') >= MIN_CHR)]
            if len(summ) else summ)
    r_all = pd.to_numeric(summ.pearson_r, errors='coerce')
    # A locus whose two denominators mean different things is not counted toward
    # the headline, however good or bad its r looks -- see allele_freq_check.py.
    ok_cmp = (summ['comparable'] == 1 if 'comparable' in summ.columns
              else pd.Series(True, index=summ.index))
    n_cmp = int(ok_cmp.sum())
    n_good = int(((r_all >= R_GOOD) & ok_cmp).sum())
    not_cmp = sorted(summ.loc[~ok_cmp, 'gene'])

    # =====================================================================
    # FIGURE 1 -- the answer
    # =====================================================================
    sec = (pd.read_csv(args.summary_secondary, sep='\t')
           if args.summary_secondary and Path(args.summary_secondary).is_file() else None)
    shared_loci = ([g for g in have.gene if g in set(sec.gene)] if sec is not None else [])

    bar_h = max(1.5, len(have) * ROW_IN)
    cmp_h = max(1.10, len(shared_loci) * ROW_IN * 1.15) if shared_loci else 0.0
    plot_h = bar_h + 1.05 + (cmp_h + 0.75 if shared_loci else 0.0)
    fig = plt.figure(figsize=(S.COL_DOUBLE, plot_h))
    if shared_loci:
        gs = fig.add_gridspec(2, 3, width_ratios=[1.25, 1.0, 0.85],
                              height_ratios=[bar_h, cmp_h])
        ax_d = fig.add_subplot(gs[1, :])
    else:
        gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 0.85])
        ax_d = None
    ax_p, ax_r, ax_u = (fig.add_subplot(gs[0, i]) for i in range(3))

    # (a) every allele at every locus on one axis.
    # NO set_aspect('equal') HERE. It shrank this cell to a square inside a
    # gridspec column sized for something wider, which left panel (a) 0.73 in
    # shorter than its siblings and put the three bold panel letters on two
    # different lines. ROW_IN is chosen so the bar panels are close to this
    # panel's width, which keeps the identity line near 45 degrees anyway.
    pooled = chk[chk.n_chr_ctrl >= MIN_CHR].copy()
    scatter(ax_p, pooled, s_pt=9, n_label=2, n_out=3, keep_gene=True)
    ax_p.set_xlabel('reference frequency', fontsize=FS_AXIS)
    ax_p.set_ylabel('observed frequency', fontsize=FS_AXIS)
    ax_p.tick_params(labelsize=FS_SMALL + 2)
    panel_head(ax_p, 'a', 'all loci pooled')
    S.despine(ax_p)

    if len(have):
        idx = np.arange(len(have))
        hv_cmp = ((have['comparable'] == 1).to_numpy() if 'comparable' in have.columns
                  else np.ones(len(have), bool))
        # A non-comparable locus is drawn hollow: the bar is still there, because
        # hiding it would leave a reader wondering where the locus went, but it
        # cannot be read as a concordance.
        ax_r.barh(idx, have.pearson_r.astype(float), height=0.66,
                  color=np.where(hv_cmp, S.DATA_DARK, 'white'),
                  edgecolor=np.where(hv_cmp, S.DATA_DARK, S.NEUTRAL_D),
                  hatch=None, linewidth=0.9)
        ax_r.set_yticks(idx)
        ax_r.set_yticklabels([g + (' †' if not c else '')
                              for g, c in zip(have.gene, hv_cmp)],
                             fontstyle='italic', fontsize=FS_SMALL)
        ax_r.set_xlim(0, 1.02)
        ax_r.set_xticks([0, 0.5, R_GOOD])
        ax_r.axvline(R_GOOD, color=S.ACCENT, lw=0.8, ls=':', zorder=0)
        ax_r.invert_yaxis()
        ax_r.set_xlabel('Pearson $r$', fontsize=FS_AXIS)
        ax_r.tick_params(axis='x', labelsize=FS_SMALL + 1)

        u = unc.reindex(have.gene)
        ax_u.barh(idx, u['share'].to_numpy(float), height=0.66, color=S.NEUTRAL_D)
        ax_u.set_yticks(idx); ax_u.set_yticklabels([])
        ax_u.tick_params(axis='y', length=0)
        ax_u.set_ylim(*ax_r.get_ylim())
        ax_u.xaxis.set_major_formatter(lambda v, _: f'{v:.0%}')
        ax_u.locator_params(axis='x', nbins=4)
        ax_u.tick_params(axis='x', labelsize=FS_SMALL + 1)
        # Named in full: `unconfirmed` alone says nothing about what the share is
        # of, and the axis is a percentage of control chromosomes.
        ax_u.set_xlabel('% of chromosomes', fontsize=FS_AXIS)
    else:
        for ax in (ax_r, ax_u):
            ax.text(0.5, 0.5, 'too few control\nchromosomes', ha='center', va='center',
                    transform=ax.transAxes, color=S.INK_SOFT, fontsize=FS_SMALL)
            ax.set_xticks([]); ax.set_yticks([])
    panel_head(ax_r, 'b', 'concordance')
    panel_head(ax_u, 'c', 'unconfirmed')
    S.despine(ax_r); S.despine(ax_u)

    # ── (d) the same concordance under BOTH reference panels ────────────────
    # The tables for the second panel have existed since it was wired up and
    # nothing drew them: a reader was asked to accept that the two panels agree
    # without being shown it. The shared-allele count is printed on each bar
    # because that is what the agreement costs -- the smaller panel reaches the
    # same r on roughly a seventh of the evidence.
    axes_all = [ax_p, ax_r, ax_u]
    if ax_d is not None:
        sh = sec.set_index('gene')
        pr = have.set_index('gene')
        x = np.arange(len(shared_loci))
        wb = 0.36
        r1 = [float(pr.loc[g, 'pearson_r']) for g in shared_loci]
        r2 = [float(sh.loc[g, 'pearson_r']) for g in shared_loci]
        n1 = [int(pr.loc[g, 'n_shared']) for g in shared_loci]
        n2 = [int(sh.loc[g, 'n_shared']) for g in shared_loci]
        for xi, (v1, v2) in enumerate(zip(r1, r2)):
            ax_d.plot([xi, xi], [v1, v2], color=S.NEUTRAL_D, lw=1.0, zorder=1)
        ax_d.scatter(x, r1, s=110, color=S.DATA_DARK, zorder=3,
                     label=args.reference_name)
        ax_d.scatter(x, r2, s=42, color=S.DATA, zorder=4,
                     label=args.reference_name_secondary)
        for xi, (v1, v2, c1, c2) in enumerate(zip(r1, r2, n1, n2)):
            ax_d.annotate(f'{c1}', (xi, v1), textcoords='offset points',
                          xytext=(0, 7 if v1 >= v2 else -7), ha='center',
                          va='bottom' if v1 >= v2 else 'top',
                          fontsize=FS_SMALL, color=S.INK_SOFT)
            ax_d.annotate(f'{c2}', (xi, v2), textcoords='offset points',
                          xytext=(0, 7 if v2 > v1 else -7), ha='center',
                          va='bottom' if v2 > v1 else 'top',
                          fontsize=FS_SMALL, color=S.INK_SOFT)
        ax_d.set_xticks(x)
        ax_d.set_xticklabels(shared_loci, fontstyle='italic')
        lo_ = min(r1 + r2 + [R_GOOD])
        ax_d.set_ylim(lo_ - (1.0 - lo_) * 0.55, 1.0 + (1.0 - lo_) * 0.55)
        ax_d.set_xlim(-0.5, len(shared_loci) - 0.5)
        ax_d.axhline(R_GOOD, color=S.ACCENT, lw=0.8, ls=':', zorder=0)
        ax_d.set_ylabel('Pearson $r$', fontsize=FS_AXIS)
        ax_d.tick_params(labelsize=FS_SMALL + 1)
        S.legend_inside(ax_d, ax_d.get_legend_handles_labels()[0], loc='lower left',
                        fontsize=FS_SMALL + 1)
        panel_head(ax_d, 'd', 'the same loci under both panels')
        S.despine(ax_d)
        axes_all.append(ax_d)

    lay_out(fig, axes_all, plot_h=plot_h, top_pad=0.34,
            wspace=0.28, hspace=0.62,
            title=f'Do our control allele frequencies reproduce {args.reference_name}?')
    fig.savefig(args.out_png)
    plt.close(fig)

    worst = summ.loc[r_all.idxmin()] if len(summ) else None
    shared = dict(zip(summ.gene, summ.get('n_shared', pd.Series(dtype=int))))
    figure_doc.write_doc(
        args.out_png,
        title=f'HLA allele frequencies against {args.reference_name}',
        question='Do the alleles this component calls occur at the frequencies published '
                 'for a Japanese population — or has the typing moved the spectrum?',
        panels=[
            ('a', 'All loci pooled',
             f'Every 2-field allele seen in the controls or in the reference, at every '
             f'locus, on one axis. x is the reference frequency, y is ours, the dashed '
             f'line is equality. Error bars are 95 % Wilson intervals on the observed '
             f'frequency. Named in red: the two commonest and the three furthest from '
             f'the line. The per-locus version is {Path(args.out_png_loci).name}.'),
            ('b', 'Per-locus concordance',
             f'Pearson $r$ between the two frequency vectors, per locus, with a dotted '
             f'line at {R_GOOD:.2f}. There is deliberately no rank correlation beside it: '
             f'most of each vector is an allele only one source carries, so a rank '
             f'statistic over it measures panel coverage rather than agreement.'),
            ('c', 'Unconfirmed share by locus',
             'The fraction of control chromosomes carrying an allele the reference does '
             'not list, per locus, on panel (b)\'s rows. Pooled over loci it reproduces '
             'the cohort figure in typing_confound.png exactly.'),
        ],
        numbers={
            'loci compared': f'{len(summ)}',
            f'loci with r >= {R_GOOD:.2f}': f'{n_good}',
            'control samples': f'{n_ctrl // 2:,}',
            'reference individuals': f'{n_ref // 2:,}',
            'lowest r': (f'{float(worst.pearson_r):.4f} ({worst.gene})'
                         if worst is not None else 'n/a'),
            'highest r': f'{r_all.max():.4f}',
            'worst unconfirmed locus': f'{unc["share"].idxmax()} {unc["share"].max():.1%}',
            'best unconfirmed locus': f'{unc["share"].idxmin()} {unc["share"].min():.1%}',
            'alleles compared': f'{len(chk):,}',
            'alleles both sources carry, fewest': (f'{min(shared.values())}'
                                                   if shared else 'n/a'),
            'alleles both sources carry, most': (f'{max(shared.values())}'
                                                 if shared else 'n/a'),
        },
        tables=[('Per-locus concordance and unconfirmed share',
                 summ.merge(unc['share'].rename('unconfirmed_share').reset_index(),
                            on='gene', how='left'))],
        interpretation=(
            'Agreement here is evidence that the typing is not systematically '
            'mis-assigning alleles — the one thing every other QC output in this '
            'component is blind to, because call rate and field depth measure whether a '
            'call was produced, not whether it was right. Disagreement at a COMMON '
            'allele is the informative failure: it means a frequent haplotype is being '
            'read as something else. Disagreement in the tail is expected and mostly '
            'reflects which alleles the reference happens to carry. Panels (b) and (c) '
            'are deliberately on the same rows: a locus can have high concordance and '
            'still put a large share of chromosomes on alleles the reference never '
            'lists, and the two together say which of those is happening.'),
        limits=[
            'This is a comparison of POPULATION FREQUENCIES. It cannot say that any '
            'given sample was typed correctly, and a set of errors that happens to '
            'preserve the frequency spectrum is invisible to it. Only typing samples '
            'with known types measures accuracy — see docs/OPEN_QUESTIONS.md §2.',
            f'Only the {len(summ)} loci the reference carries can be checked at all. '
            f'Every other locus this component types has no external reference here, so '
            f'its calls are unverified rather than verified-and-passing.',
            'The reference is one Japanese population, not Japanese people in general. '
            'A locus where our cohort and the reference panel differ in ancestry will '
            'differ in frequency for reasons that have nothing to do with typing.',
            'Cases are excluded from the comparison and shown only in the table. They '
            'are a disease series, and the MHC is where a disease series is least '
            'expected to match a population reference.',
        ],
        methods_ref='../../docs/METHODS.md')

    # =====================================================================
    # FIGURE 2 -- the record
    # =====================================================================
    ncol = min(LOCI_NCOL, len(genes))
    nrow = -(-len(genes) // ncol)
    loci_h = LOCI_PANEL_IN * nrow + 0.55
    fig2, axes2 = plt.subplots(nrow, ncol, figsize=(S.COL_DOUBLE, loci_h), squeeze=False)
    flat = axes2.ravel()
    for i, (ax, g) in enumerate(zip(flat, genes)):
        # The x label goes on every panel in the LAST OCCUPIED ROW OF ITS COLUMN.
        # `i + ncol >= len(genes)` is the same thing, but it made panels (i) and
        # (j) carry a label while (f), (g) and (h) beside them did not, because
        # the final row is short. Reading it by column removes the asymmetry.
        last_in_col = (i + ncol >= len(genes))
        locus_panel(ax, chk[chk.gene == g].copy(), g, string.ascii_lowercase[i],
                    xlab=last_in_col or (i // ncol == nrow - 2 and i % ncol >= len(genes) % ncol
                                         and len(genes) % ncol),
                    ylab=(i % ncol == 0),            # leftmost column
                    flag=(g in flag_set))
    for ax in flat[len(genes):]:
        ax.axis('off')

    lay_out(fig2, [a for a in flat[:len(genes)]], plot_h=loci_h, top_pad=0.34,
            wspace=0.46, hspace=0.62,
            title='Observed against reference frequency, one panel per locus')
    fig2.savefig(args.out_png_loci)
    plt.close(fig2)

    figure_doc.write_doc(
        args.out_png_loci,
        title=f'HLA allele frequencies against {args.reference_name}, locus by locus',
        question='At which locus, and at which allele, do the observed frequencies and '
                 'the reference disagree?',
        panels=[(string.ascii_lowercase[i], f'HLA-{g}',
                 f'Every 2-field allele seen at HLA-{g} in the controls or in the '
                 f'reference. x is the reference frequency, y is ours, and the dashed '
                 f'line is equality. Error bars are 95 % Wilson intervals. Named in red: '
                 f'the two commonest and the one furthest from the line. The axis limits '
                 f'are set from this locus alone, so a panel cannot be compared with its '
                 f'neighbour by eye — that is what panel (b) of '
                 f'{Path(args.out_png).name} is for. '
                 + (f'The two sources share {shared.get(g)} alleles at this locus.'
                    if shared.get(g) is not None else ''))
                for i, g in enumerate(genes)],
        numbers={f'r, {r.gene}': f'{float(r.pearson_r):.4f}' for _, r in summ.iterrows()},
        interpretation=(
            'This figure is the evidence, not the argument; the argument is in '
            f'{Path(args.out_png).name}. Read a panel when a locus matters to a specific '
            'question — which allele is off the line, and by how much relative to its '
            'interval. A panel whose points sit on the line at every frequency is a '
            'locus whose common haplotypes are being read correctly.'
            + (f' † {flagged} is marked as not comparable: the reference assigns an '
               f'allele to every chromosome and encodes no absence for that gene, so its '
               f'denominator is not the same quantity as ours.' if flagged else '')),
        limits=[
            'Each panel is scaled to its own locus, so panel-to-panel comparison by eye '
            'is not meaningful. Use the concordance panel for that.',
            'A locus with few observed alleles has a short axis and looks tidier than a '
            'locus with many; tidiness here is allele diversity, not typing quality.',
        ],
        methods_ref='../../docs/METHODS.md')

    print(f'[plot_allele_freq] {len(genes)} loci, {n_ctrl // 2:,} controls vs '
          f'{args.reference_name} ({n_ref // 2:,} individuals); '
          f'{n_good}/{len(summ)} loci at r >= {R_GOOD:.2f}')
    print(f'[plot_allele_freq] -> {args.out_png} (summary), {args.out_png_loci} (record)')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:                       # noqa: BLE001
        print(f'Error in plot_allele_freq: {e}', file=sys.stderr)
        sys.exit(1)
