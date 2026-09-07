#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The one quality problem this component actually has, drawn.
#
#           typing_qc.png measures COMPLETENESS -- was a call made, how deeply
#           did it resolve -- and on this cohort both are near-saturated. Nothing
#           there can see a call that was made confidently and is wrong.
#
#           The frequency check can, and what it found is in
#           docs/OPEN_QUESTIONS.md section 3: the controls carry a materially
#           larger share of calls the reference panel never sees, and that
#           difference tracks phenotype. This figure is that finding, built to
#           answer the two questions a reader will ask about it -- is it depth,
#           and is it an artefact of the reference panel -- with data rather than
#           with prose.
#
#           IT IS NOT DEPTH: cases and controls sit at the same measured depth
#           and the gap survives matching on it. Depth does have an effect,
#           visible in panel (a) where phenotype is held fixed within the case
#           platforms, but there is no case-control depth difference for that
#           effect to act on.
#
#           IT IS NOT THE PANEL EITHER, and panel (d) is why that is now shown
#           rather than asserted. An allele the panel never carries is not
#           necessarily a wrong call: a panel samples a finite number of
#           chromosomes and a rare allele can be missing from it. So the ABSOLUTE
#           level is a property of the panel. Running the same calls against a
#           second, far smaller panel moves that level by a third and leaves the
#           RATIO between the groups where it was. The contrast is the readable
#           quantity; the level is not.
#
#           NO CAPTION IS RENDERED INTO THIS PNG. Every word of explanation lives
#           in report/hla_typing_report.qmd, where it can be written for a reader
#           who has never heard of HLA and can be translated. A caption block was
#           taking 31 % of this canvas; the figure is now only the figure.
#
#           NO NUMBER IS WRITTEN INTO THIS FILE either -- write_doc() still emits
#           the sidecar .md, and every value in it is an f-string field read from
#           the tables, because the last time a reference panel changed the
#           captions kept describing the old one for three months.
# Component: hla.typing
# ---------------------------------------------------------------------------
import argparse
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402
from scipy import stats                  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'analysis' / '_shared' / 'scripts'))
import plot_style as S                   # noqa: E402
import figure_doc                        # noqa: E402

PLOT_H = 5.15
# Point-label placement. A label is offset by ITS OWN marker radius plus this pad,
# then tried on four sides and scored against every OTHER marker -- the previous
# version offset only by its own radius to the right, which put the `NovaSeq 30x`
# label on top of the `G400RS 30x` marker and invited the reader to attach the
# wrong platform to the wrong point. That inverts panel (a)'s argument.
LABEL_PAD_PT = 4.0
FS_SMALL = 6.8        # one small size for every in-axes annotation in this figure
# The depth window panel (b) matches on. Chosen as the overlap of the two
# platforms that carry the comparison, not tuned to the answer.
WIN_LO, WIN_HI = 17.0, 21.0


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample', required=True,
                   help='allele_frequency_sample.tsv from the PRIMARY reference panel')
    p.add_argument('--reference-name', default='jMorp 61KJPN-HLA',
                   help='how the primary panel is named on the figure')
    p.add_argument('--sample-secondary', required=True,
                   help='allele_frequency_sample.tsv from the SECOND panel. Panel (d) is the '
                        'whole reason both panels are run: it separates what is a property of '
                        'the reference from what is a property of the typing.')
    p.add_argument('--reference-name-secondary', default='1000 Genomes JPT')
    p.add_argument('--control-group', default='AGP3K')
    p.add_argument('--out-png', required=True)
    p.add_argument('--out-table', required=True)
    return p.parse_args()


def share(d):
    """Chromosomes on an unconfirmed allele, over chromosomes called."""
    n = d.n_chr.sum()
    return (d.n_unconfirmed.sum() / n) if n else np.nan


def fisher(a, b):
    ta = [[a.n_unconfirmed.sum(), a.n_chr.sum() - a.n_unconfirmed.sum()],
          [b.n_unconfirmed.sum(), b.n_chr.sum() - b.n_unconfirmed.sum()]]
    if min(min(r) for r in ta) < 0 or not a.n_chr.sum() or not b.n_chr.sum():
        return np.nan
    return float(stats.fisher_exact(ta)[1])


def load(path, control_group):
    d = pd.read_csv(path, sep='\t', dtype={'sample_id': str})
    for col in ('n_chr', 'n_unconfirmed', 'group'):
        if col not in d.columns:
            raise SystemExit(f'ABORT: {path} has no `{col}` column')
    d['observed_depth'] = pd.to_numeric(d.get('observed_depth'), errors='coerce')
    ctrl = d[d.group == control_group]
    case = d[d.group != control_group]
    if not len(ctrl) or not len(case):
        raise SystemExit(f'ABORT: need both groups in {path}; found {sorted(set(d.group))}')
    return d, ctrl, case


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

def place_point_labels(ax, xs, ys, sizes, texts):
    """Label each point on whichever side is clearest of the OTHER points.

    Four candidates per point, each offset by that point's own marker radius so
    the text never lands on its own bubble. Each candidate is scored by its
    distance, in display space, to every other marker and to every label already
    placed; the best-scoring candidate wins. Deterministic, and with five points
    it is exact rather than heuristic.
    """
    fig = ax.figure
    fig.canvas.draw()
    trans = ax.transData
    pts = trans.transform(np.column_stack([xs, ys]))
    rad = np.sqrt(np.asarray(sizes, float) / np.pi) * fig.dpi / 72.0
    box = ax.get_window_extent()
    placed, out = [], []
    for i, (px, py) in enumerate(pts):
        r = rad[i] + LABEL_PAD_PT * fig.dpi / 72.0
        # Text width estimated from the string; the artist does not exist yet, and
        # a candidate that would push the label out of the axes has to be rejected
        # BEFORE it is drawn. Without this the widest label still ran into the
        # gutter, because staying clear of the other markers scored better than
        # staying inside the frame.
        w = len(texts[i]) * FS_SMALL * 0.60 * fig.dpi / 72.0
        h = FS_SMALL * 1.2 * fig.dpi / 72.0
        best, best_score = None, -1e18
        for dx, dy, ha, va in ((r, 0, 'left', 'center'), (-r, 0, 'right', 'center'),
                               (0, r, 'center', 'bottom'), (0, -r, 'center', 'top')):
            cx, cy = px + dx, py + dy
            x0 = cx if ha == 'left' else (cx - w if ha == 'right' else cx - w / 2)
            y0 = cy if va == 'bottom' else (cy - h if va == 'top' else cy - h / 2)
            outside = (max(0.0, box.x0 - x0) + max(0.0, (x0 + w) - box.x1)
                       + max(0.0, box.y0 - y0) + max(0.0, (y0 + h) - box.y1))
            d = [np.hypot(cx - pts[j][0], cy - pts[j][1]) - rad[j]
                 for j in range(len(pts)) if j != i]
            d += [np.hypot(cx - qx, cy - qy) for qx, qy in placed]
            score = (min(d) if d else 1e9) - 10.0 * outside
            if score > best_score:
                best, best_score = (dx, dy, ha, va, cx, cy), score
        dx, dy, ha, va, cx, cy = best
        placed.append((cx, cy))
        out.append(ax.annotate(texts[i], (xs[i], ys[i]), textcoords='offset points',
                               xytext=(dx / fig.dpi * 72.0, dy / fig.dpi * 72.0),
                               ha=ha, va=va, fontsize=FS_SMALL, color=S.INK_SOFT))
    return out


def stat_note(ax, text, *, y=0.985, va='top'):
    """A computed statistic INSIDE the axes.

    It used to be written with ax.set_title(), which panel_tag() then overwrote
    with the panel's short title -- both write loc='center' -- so every P value
    this script computed was silently absent from the figure it was computed for.
    """
    ax.text(0.5, y, text, transform=ax.transAxes, ha='center', va=va,
            fontsize=FS_SMALL, color=S.INK_SOFT, clip_on=False,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=1.2))


def main():
    args = parse_args()
    S.setup_style('slide')
    d, ctrl, case = load(args.sample, args.control_group)
    d2, ctrl2, case2 = load(args.sample_secondary, args.control_group)
    if 'platform' not in d.columns:
        raise SystemExit(f'ABORT: {args.sample} has no platform column')

    fig, axes = plt.subplots(2, 2, figsize=(S.COL_DOUBLE, PLOT_H))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()
    C_CTRL, C_CASE = S.DATA, S.ACCENT

    # -- (a) the whole story on one axis: share against MEASURED depth --------
    # One point per platform, so the two effects are separable by eye: within the
    # case platforms the points fall with depth, and the control platform sits far
    # above the case platform it shares a depth with.
    rows = []
    for pf, g in d.groupby('platform'):
        is_ctrl = (g.group == args.control_group).mean() > 0.5
        rows.append({'platform': pf, 'n': len(g), 'depth': g.observed_depth.median(),
                     'share': share(g), 'is_ctrl': is_ctrl})
    pl = pd.DataFrame(rows).sort_values('depth')
    for is_ctrl, col, lab in ((True, C_CTRL, 'controls'), (False, C_CASE, 'cases')):
        sel = pl[pl.is_ctrl == is_ctrl]
        ax_a.scatter(sel.depth, sel.share, s=np.sqrt(sel.n) * 9, color=col, alpha=0.85,
                     edgecolors='white', linewidths=0.8, zorder=3, label=lab)
    # The vendor prefix is dropped HERE ONLY, and only from the point labels: five
    # full names on a 2.4 in axis overlap into an unreadable band. The table and
    # every other panel keep the name in full.
    short = {q: str(q).split('-')[-1] if '-' in str(q) else str(q) for q in pl.platform}
    # Limits BEFORE placing labels: placement measures display coordinates, so the
    # axes have to be final first. The x range is widened symmetrically because a
    # label can now go on either side of its point.
    xpad = (pl.depth.max() - pl.depth.min()) * 0.16 + 2.0
    ax_a.set_xlim(pl.depth.min() - xpad, pl.depth.max() + xpad)
    ax_a.set_ylim(-pl.share.max() * 0.12, pl.share.max() * 1.30)
    ax_a.set_xlabel('median depth (×)')
    ax_a.set_ylabel('unconfirmed\nchromosomes')
    ax_a.yaxis.set_major_formatter(lambda v, _: f'{v:.0%}')
    place_point_labels(ax_a, pl.depth.to_numpy(), pl.share.to_numpy(),
                       (np.sqrt(pl.n) * 9).to_numpy(),
                       [short[q] for q in pl.platform])

    dd = d.dropna(subset=['observed_depth'])
    cc, ca = dd[dd.group == args.control_group], dd[dd.group != args.control_group]
    p_depth = float(stats.mannwhitneyu(ca.observed_depth, cc.observed_depth)[1]) \
        if len(cc) and len(ca) else np.nan
    panel_head(ax_a, 'a', 'share vs depth')
    # markerscale, because this is a SIZE-ENCODED scatter: without it the legend
    # reproduces the handles at their data size and the control swatch is as big
    # as the largest real datum, which a reader can legitimately read as a point.
    leg = S.legend_inside(ax_a, ax_a.get_legend_handles_labels()[0], loc='upper right',
                          fontsize=FS_SMALL + 1.5)
    if leg is not None:
        for h in leg.legend_handles:
            h.set_sizes([18])
    if p_depth == p_depth:
        stat_note(ax_a, f'{cc.observed_depth.median():.1f}× vs '
                        f'{ca.observed_depth.median():.1f}×, '
                        f'MW {S.p_tex(p_depth)}', y=0.02, va='bottom')
    S.despine(ax_a)

    # -- (b) the same comparison with depth held fixed -----------------------
    w = d[(d.observed_depth >= WIN_LO) & (d.observed_depth <= WIN_HI)]
    top = (w.groupby('platform').size().sort_values(ascending=False).index.tolist())[:2]
    bars, labs, cols = [], [], []
    for pf in top:
        g = w[w.platform == pf]
        bars.append(share(g))
        labs.append(f'{pf}\nn={len(g):,}  {g.observed_depth.median():.1f}×')
        cols.append(C_CTRL if (g.group == args.control_group).mean() > 0.5 else C_CASE)
    p_matched = fisher(w[w.platform == top[0]], w[w.platform == top[1]]) if len(top) == 2 else np.nan
    if len(bars) < 2:
        # An empty or one-sided window is a real outcome once the cohort is
        # restricted, and it has to be legible AS one rather than as a blank axis.
        ax_b.text(0.5, 0.5, f'no two platforms overlap\nat {WIN_LO:.0f}-{WIN_HI:.0f}x',
                  transform=ax_b.transAxes, ha='center', va='center',
                  fontsize=FS_SMALL, color=S.INK_SOFT)
        ax_b.set_xticks([])
    else:
        ax_b.bar(range(len(bars)), bars, color=cols, width=0.5)
        ax_b.set_xticks(range(len(bars)))
        ax_b.set_xticklabels(labs, fontsize=FS_SMALL)
        ax_b.set_ylim(0, max(bars) * 1.30)
        if p_matched == p_matched:
            stat_note(ax_b, f'Fisher {S.p_tex(p_matched)}')
    ax_b.yaxis.set_major_formatter(lambda v, _: f'{v:.0%}')
    # No y label here or on (d): all three share one quantity and one axis
    # formatter, and three copies of a two-line label ate 624 px of the gutter.
    panel_head(ax_b, 'b', f'matched at {WIN_LO:.0f}-{WIN_HI:.0f}×')
    S.despine(ax_b)

    # -- (c) per sample, so the reader sees a distribution and not a mean ----
    hi = int(d.n_unconfirmed.max())
    bins = np.arange(-0.5, hi + 1.5)
    for g, col, lab in ((ctrl, C_CTRL, f'controls (n={len(ctrl):,})'),
                        (case, C_CASE, f'cases (n={len(case):,})')):
        h, _ = np.histogram(g.n_unconfirmed, bins=bins)
        ax_c.step(np.arange(hi + 1), h / len(g), where='mid', color=col, lw=1.6, label=lab)
    n_lo, n_hi = int(d.n_chr.min()), int(d.n_chr.max())
    denom = f'{n_hi}' if n_lo == n_hi else f'{n_lo}-{n_hi}'
    ax_c.set_xlabel('unconfirmed chromosomes')
    ax_c.set_ylabel('fraction of samples')
    ax_c.set_xticks(range(0, hi + 1, 2))
    panel_head(ax_c, 'c', 'per sample')
    S.legend_inside(ax_c, ax_c.get_legend_handles_labels()[0], loc='upper right')
    S.despine(ax_c)

    # -- (d) the same contrast under both reference panels -------------------
    # The point of running two panels, made visible. The level is a property of
    # the panel; the ratio between the groups is not.
    panels = [(args.reference_name, ctrl, case),
              (args.reference_name_secondary, ctrl2, case2)]
    x = np.arange(len(panels))
    wbar = 0.34
    sh_c = [share(c) for _, c, _ in panels]
    sh_k = [share(k) for _, _, k in panels]
    ax_d.bar(x - wbar / 2, sh_c, wbar, color=C_CTRL)
    ax_d.bar(x + wbar / 2, sh_k, wbar, color=C_CASE)
    for i, (c, k) in enumerate(zip(sh_c, sh_k)):
        ax_d.annotate(f'{c / k:.2f}×', (i, max(c, k)), textcoords='offset points',
                      xytext=(0, 7), ha='center', fontsize=FS_SMALL, color=S.INK)
    ax_d.set_xticks(x)
    ax_d.set_xticklabels([nm for nm, _c, _k in panels], fontsize=FS_SMALL)
    ax_d.set_ylim(0, max(sh_c + sh_k) * 1.26)
    ax_d.yaxis.set_major_formatter(lambda v, _: f'{v:.0%}')
    panel_head(ax_d, 'd', 'both panels')
    S.despine(ax_d)

    # -- the table, which carries everything the four panels cannot ----------
    tab = pd.DataFrame([
        {'stratum': 'controls', 'n': len(ctrl), 'median_depth': round(ctrl.observed_depth.median(), 2),
         'n_chr': int(ctrl.n_chr.sum()), 'n_unconfirmed': int(ctrl.n_unconfirmed.sum()),
         'share_unconfirmed': round(share(ctrl), 5)},
        {'stratum': 'cases', 'n': len(case), 'median_depth': round(case.observed_depth.median(), 2),
         'n_chr': int(case.n_chr.sum()), 'n_unconfirmed': int(case.n_unconfirmed.sum()),
         'share_unconfirmed': round(share(case), 5)},
    ] + [
        {'stratum': f'platform:{r.platform}', 'n': int(r.n), 'median_depth': round(r.depth, 2),
         'n_chr': int(d[d.platform == r.platform].n_chr.sum()),
         'n_unconfirmed': int(d[d.platform == r.platform].n_unconfirmed.sum()),
         'share_unconfirmed': round(r.share, 5)} for _, r in pl.iterrows()
    ] + [
        {'stratum': f'matched {WIN_LO:.0f}-{WIN_HI:.0f}x:{pf}', 'n': int((w.platform == pf).sum()),
         'median_depth': round(w[w.platform == pf].observed_depth.median(), 2),
         'n_chr': int(w[w.platform == pf].n_chr.sum()),
         'n_unconfirmed': int(w[w.platform == pf].n_unconfirmed.sum()),
         'share_unconfirmed': round(share(w[w.platform == pf]), 5)} for pf in top
    ] + [
        {'stratum': f'{args.reference_name_secondary}:{lab}', 'n': len(g),
         'median_depth': round(g.observed_depth.median(), 2),
         'n_chr': int(g.n_chr.sum()), 'n_unconfirmed': int(g.n_unconfirmed.sum()),
         'share_unconfirmed': round(share(g), 5)}
        for lab, g in (('controls', ctrl2), ('cases', case2))
    ])
    p_2nd = fisher(ctrl2, case2)
    tab['fisher_p_matched'] = ''
    tab['mannwhitney_p_depth'] = ''
    tab['fisher_p_group'] = ''
    if len(top) == 2 and p_matched == p_matched:
        tab.loc[tab.index[-3], 'fisher_p_matched'] = f'{p_matched:.3g}'
    if p_depth == p_depth:
        tab.loc[tab.index[1], 'mannwhitney_p_depth'] = f'{p_depth:.3g}'
    p_1st = fisher(ctrl, case)
    if p_1st == p_1st:
        tab.loc[tab.index[1], 'fisher_p_group'] = f'{p_1st:.3g}'
    if p_2nd == p_2nd:
        tab.loc[tab.index[-1], 'fisher_p_group'] = f'{p_2nd:.3g}'
    tab.to_csv(args.out_table, sep='\t', index=False)

    n_loci = int(round(n_hi / 2))
    lay_out(fig, list(axes.ravel()), plot_h=PLOT_H, top_pad=0.34,
            wspace=0.34, hspace=0.52,
            title='Typing quality by group: depth, the reference panel, or neither?')
    fig.savefig(args.out_png)
    plt.close(fig)

    figure_doc.write_doc(
        args.out_png,
        title='Typing accuracy differs between the groups; it is not depth and not the panel',
        question='The controls are typed measurably worse than the cases. Is that sequencing '
                 'depth — which could be controlled for — or an artefact of the reference '
                 'panel, or something that can be neither corrected nor explained away?',
        panels=[
            ('a', 'Unconfirmed share against measured depth',
             'One point per sequencing platform, x its median measured depth, y the share of '
             'called chromosomes carrying an allele the reference panel never sees; area is '
             'proportional to √n. Read it twice: along the case platforms alone, phenotype is '
             'held fixed and the share falls as depth rises — depth is real. Then compare the '
             'control platform with the case platform beside it at the same depth. The note in '
             'the panel gives the two group medians and a Mann–Whitney test of their equality.'),
            ('b', 'The comparison at matched depth',
             f'The two platforms carrying the comparison, restricted to samples between '
             f'{WIN_LO:.0f}× and {WIN_HI:.0f}× so their median depths are equal. If depth were '
             f'the cause the bars would be level. Fisher exact on the 2×2 of unconfirmed against '
             f'confirmed chromosomes.'),
            ('c', 'Per sample',
             'How many of a sample\'s called chromosomes at the reference loci sit on an '
             'unconfirmed allele. Shown as a distribution because a mean hides that most '
             'samples have none and a tail has several. The denominator is printed on the axis '
             'and is read from the data, not assumed.'),
            ('d', 'Both reference panels',
             'The control and case shares computed twice from the same calls, once against each '
             'panel, with the control/case ratio printed above each pair. This is the panel that '
             'licenses the reading of every other one: if the ratio moved with the panel, the '
             'finding would be about the reference rather than about the typing.'),
        ],
        numbers={
            f'controls, {args.reference_name}': f'{share(ctrl):.2%}',
            f'cases, {args.reference_name}': f'{share(case):.2%}',
            'ratio, primary panel': f'{share(ctrl) / share(case):.2f}',
            f'controls, {args.reference_name_secondary}': f'{share(ctrl2):.2%}',
            f'cases, {args.reference_name_secondary}': f'{share(case2):.2%}',
            'ratio, secondary panel': f'{share(ctrl2) / share(case2):.2f}',
            'Fisher P, primary panel': f'{p_1st:.3g}',
            'Fisher P, secondary panel': f'{p_2nd:.3g}',
            'control median depth': f'{cc.observed_depth.median():.2f}×',
            'case median depth': f'{ca.observed_depth.median():.2f}×',
            'Mann-Whitney P, depth': f'{p_depth:.3g}',
            f'Fisher P, matched {WIN_LO:.0f}-{WIN_HI:.0f}x': (f'{p_matched:.3g}'
                                                              if p_matched == p_matched else 'n/a'),
        },
        tables=[('Every stratum behind the four panels, as written to typing_confound.tsv', tab)],
        interpretation=(
            'Three explanations are available and two of them are ruled out here. Depth '
            'genuinely affects typing — panel (a) shows it cleanly within the case platforms, '
            'where phenotype cannot be responsible — but cases and controls sit at the same '
            'measured depth and the gap survives matching on it, so that effect has nothing to '
            f'act on between the groups. The reference panel genuinely sets the level — panel '
            f'(d) shows the control level move from {share(ctrl):.2%} against '
            f'{args.reference_name} to {share(ctrl2):.2%} against '
            f'{args.reference_name_secondary} — but the ratio '
            f'between the groups is stable across both '
            f'({share(ctrl) / share(case):.2f} and {share(ctrl2) / share(case2):.2f}), so '
            f'the contrast is not an artefact of what the panel happens to carry. What is left '
            'is platform and cohort: every control is an AGP3K sample on HiSeqX, every case is '
            'a PH sample on something else. Those are perfectly confounded with each other and '
            'with phenotype, so which of read length, library chemistry, alignment reference or '
            'batch is responsible cannot be determined from these data. A common allele '
            'depleted more in controls than in cases reads as enrichment in cases — a false '
            'risk association from the typing alone.'),
        limits=[
            'The absolute unconfirmed share is not an error rate. An allele a panel never '
            'carries may still be real and rare; panel (d) measures how much that matters by '
            'changing the panel. Only the difference between groups, measured on the same '
            'panel, is interpretable.',
            'Only the loci a panel carries are visible. Loci absent from a panel are not '
            'scored against it at all, which is why the two panels contribute different '
            'numbers of chromosomes per sample and why each locus list travels with its panel.',
            'The cohort is ancestry-restricted, so an unconfirmed allele can no longer be '
            'explained by the sample not being mainland Japanese. That removes an explanation; '
            'it does not remove the platform confound, which is what remains.',
            'The figure identifies the confound; it does not correct it. No within-data control '
            'exists, which is why OPEN_QUESTIONS.md §2 — typing samples of known HLA type — is '
            'now the only remaining option.',
        ],
        methods_ref='../../docs/METHODS.md')

    print(f'[plot_typing_confound] {args.reference_name}: controls {share(ctrl):.2%} vs cases '
          f'{share(case):.2%} ({share(ctrl) / share(case):.2f}x, Fisher P = {p_1st:.3g})')
    print(f'[plot_typing_confound] {args.reference_name_secondary}: controls {share(ctrl2):.2%} '
          f'vs cases {share(case2):.2%} ({share(ctrl2) / share(case2):.2f}x, '
          f'Fisher P = {p_2nd:.3g})')
    print(f'[plot_typing_confound] matched-depth Fisher P = {p_matched:.3g}; '
          f'depth Mann-Whitney P = {p_depth:.3g} -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:                       # noqa: BLE001
        print(f'Error in plot_typing_confound: {e}', file=sys.stderr)
        sys.exit(1)
