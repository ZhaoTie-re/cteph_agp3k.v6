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
#           WHY IT IS TWO FIGURES. It used to be one: a 4x4 grid of thirteen
#           near-identical per-locus scatters with a concordance panel wedged into
#           the fourteenth cell, 10.8 inches of plot answering a question that
#           takes three panels. A reader asking "do the frequencies agree?" and a
#           reader auditing one locus want different things, and the single figure
#           served neither. So:
#             allele_frequency.png       the answer: pooled agreement, per-locus
#                                        concordance, per-locus unconfirmed share
#             allele_frequency_loci.png  the record: every locus, one panel each
#
#           It is a POPULATION comparison, not an accuracy measurement. That, and
#           the size of the reference, are ON the figure rather than only in the
#           document beside it.
#
#           NO NUMBER IS WRITTEN INTO THIS FILE. The reference is named and cited
#           by the caller and its size is read from the table, because the last
#           time the panel changed the captions kept describing the old one.
# Component: hla.typing
# ---------------------------------------------------------------------------
import argparse
import string
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'analysis' / '_shared' / 'scripts'))
import plot_style as S                   # noqa: E402
import figure_doc                        # noqa: E402

SUMMARY_H = 3.30      # the three-panel answer
PANEL_H = 1.95        # per grid ROW of the per-locus record
ROW_IN = 0.215        # inches per locus row in the two bar panels. Below ~0.20
                      # the italic locus labels touch at DPA1/DPB1.
# Three labels, not six. An HLA frequency table is a handful of common alleles and a
# long tail at the origin, so naming more than the top few piles the text on top of
# itself exactly where the points are densest. An allele off the line is named too,
# but only if it is off by enough to be worth reading -- below that the label is noise.
N_LABEL = 3
MIN_DIFF_TO_NAME = 0.02
N_OUTLIER_MAX = 2     # however many alleles clear MIN_DIFF_TO_NAME, name at most
                      # this many of them. Without the cap, a locus that genuinely
                      # disagrees -- or a run with too few controls to estimate
                      # anything -- names EVERY allele and the panel becomes a wall
                      # of text exactly where the points are.
MIN_CHR = 20          # control chromosomes below which nothing is drawn. At n=2 the
                      # only observable frequencies are 0, 0.5 and 1, every Wilson
                      # interval spans the axis, and a scatter of that is noise
                      # presented as a comparison.
R_GOOD = 0.98         # the concordance a locus has to clear to be read at face value


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--check', required=True, help='allele_frequency_check.tsv')
    p.add_argument('--summary', required=True, help='allele_frequency_summary.tsv')
    p.add_argument('--reference-name', default='jMorp 61KJPN-HLA',
                   help='what the reference is called on the figure')
    p.add_argument('--reference-cite', default='ToMMo jMorp; Tadaka et al. 2023',
                   help='the citation printed in the notes')
    p.add_argument('--flag-genes', default='DRB3',
                   help='loci whose comparison does not reconcile and must be marked on '
                        'the figure rather than left to a footnote. See '
                        'docs/OPEN_QUESTIONS.md.')
    p.add_argument('--out-png', required=True, help='the summary figure')
    p.add_argument('--out-png-loci', required=True, help='the per-locus record')
    return p.parse_args()


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
    hi = float(max(np.nanmax(x), np.nanmax(y))) * 1.22 + 0.01
    ax.plot([0, hi], [0, hi], color=S.REFERENCE, lw=1.0, ls='--', zorder=1)
    named = to_label(rows, n_label, n_out)
    is_out = rows.allele.isin(named).to_numpy()
    # Wilson intervals, so a point far from the line is visibly far relative to how
    # well it is estimated rather than in absolute terms.
    err = np.vstack([y - rows.ci_lo_ctrl.to_numpy(float),
                     rows.ci_hi_ctrl.to_numpy(float) - y])
    ax.errorbar(x, y, yerr=np.clip(err, 0, None), fmt='none',
                ecolor=S.NEUTRAL_D, elinewidth=0.7, zorder=2)
    ax.scatter(x[~is_out], y[~is_out], s=s_pt, color=S.DATA, alpha=0.75,
               edgecolors='none', zorder=3)
    ax.scatter(x[is_out], y[is_out], s=s_pt + 4, color=S.ACCENT, alpha=0.9,
               edgecolors='none', zorder=4)
    labels = [al if keep_gene else al.split('*', 1)[-1] for al in rows.allele[is_out]]
    annots = [ax.annotate(lab, (xi, yi), textcoords='offset points',
                          xytext=(-5, 3) if xi > 0.62 * hi else (5, 3),
                          ha='right' if xi > 0.62 * hi else 'left',
                          fontsize=6.5, color=S.INK_SOFT)
              for xi, yi, lab in zip(x[is_out], y[is_out], labels)]
    ax.set_xlim(0, hi); ax.set_ylim(0, hi)
    ax.set_aspect('equal', adjustable='box')
    if annots:
        S.spread_labels(ax, annots, axis='y')
    return hi


def locus_panel(ax, rows, gene, letter, *, xlab=False, ylab=False, flag=False):
    """One locus of the record figure."""
    n_chr = int(rows.n_chr_ctrl.max()) if len(rows) else 0
    if n_chr < MIN_CHR:
        msg = ('no control sample\nin this run' if not n_chr
               else f'only {n_chr} control chromosomes\n— too few to compare')
        ax.text(0.5, 0.5, msg, ha='center', va='center',
                transform=ax.transAxes, color=S.INK_SOFT, fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
        S.panel_tag(ax, letter, title=f'{gene}' + (' †' if flag else ''))
        return
    scatter(ax, rows)
    if xlab:
        # The reference is named in the caption; repeating it on four x axes only
        # makes them collide.
        ax.set_xlabel('reference')
    if ylab:
        ax.set_ylabel('observed')
    S.panel_tag(ax, letter, title=f'{gene}' + (' †' if flag else ''))
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
    n_good = int((r_all >= R_GOOD).sum())

    # =====================================================================
    # FIGURE 1 -- the answer
    # =====================================================================
    bar_h = max(1.5, len(have) * ROW_IN)
    plot_h = max(SUMMARY_H, bar_h + 1.05)
    fig = plt.figure(figsize=(S.COL_DOUBLE, plot_h))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 0.85])
    ax_p, ax_r, ax_u = (fig.add_subplot(gs[0, i]) for i in range(3))

    # (a) every allele at every locus on one axis. Thirteen scatters answering one
    # question is thirteen answers; this is the question.
    pooled = chk[chk.n_chr_ctrl >= MIN_CHR].copy()
    scatter(ax_p, pooled, s_pt=9, n_label=2, n_out=4, keep_gene=True)
    ax_p.set_xlabel('reference frequency')
    ax_p.set_ylabel('observed frequency')
    S.panel_tag(ax_p, 'a', title='all loci pooled')
    S.despine(ax_p)

    # (b) per-locus concordance, and (c) beside it on the same rows.
    if len(have):
        idx = np.arange(len(have))
        ax_r.barh(idx + 0.19, have.pearson_r.astype(float), height=0.36,
                  color=S.DATA_DARK, label='Pearson $r$')
        ax_r.barh(idx - 0.19, have.spearman_rho.astype(float), height=0.36,
                  color=S.DATA, label=r'Spearman $\rho$')
        ax_r.set_yticks(idx)
        ax_r.set_yticklabels([g + (' †' if g in flag_set else '') for g in have.gene],
                             fontstyle='italic')
        # rho can be negative here -- it weights every allele equally and an HLA
        # table is mostly a tail of ties -- so the axis has to admit that rather
        # than clip it to zero and make the bar look absent.
        lo = min(0.0, float(pd.to_numeric(have.spearman_rho, errors='coerce').min()))
        ax_r.set_xlim(min(-0.05, lo * 1.1), 1.02)
        ax_r.axvline(0, color=S.INK_SOFT, lw=0.6, zorder=0)
        ax_r.axvline(R_GOOD, color=S.ACCENT, lw=0.7, ls=':', zorder=0)
        ax_r.invert_yaxis()
        ax_r.tick_params(axis='y', labelsize=plt.rcParams['ytick.labelsize'] - 2)
        ax_r.legend(loc='upper center', bbox_to_anchor=(0.5, -0.055), ncol=2,
                    frameon=False, handlelength=1.1, columnspacing=1.0,
                    fontsize=plt.rcParams['legend.fontsize'] - 2)

        u = unc.reindex(have.gene)
        ax_u.barh(idx, u['share'].to_numpy(float), height=0.62, color=S.NEUTRAL_D)
        ax_u.set_yticks(idx); ax_u.set_yticklabels([])
        ax_u.set_ylim(*ax_r.get_ylim())
        ax_u.xaxis.set_major_formatter(lambda v, _: f'{v:.0%}')
        ax_u.set_xlabel('unconfirmed')
    else:
        for ax in (ax_r, ax_u):
            ax.text(0.5, 0.5, f'fewer than {MIN_CHR} control\nchromosomes',
                    ha='center', va='center', transform=ax.transAxes,
                    color=S.INK_SOFT, fontsize=8)
            ax.set_xticks([]); ax.set_yticks([])
    S.panel_tag(ax_r, 'b')
    S.panel_tag(ax_u, 'c')
    S.despine(ax_r); S.despine(ax_u)

    S.caption_block(
        fig, plot_h=plot_h,
        title=(f'Control frequencies reproduce {args.reference_name} at {n_good} of '
               f'{len(summ)} loci; the exceptions are named on the figure, not in a '
               f'footnote.'),
        panels=[
            'Every allele at every locus, observed against reference; dashed line is '
            'equality, bars are 95 % Wilson intervals.',
            r'Per-locus Pearson $r$ and Spearman $\rho$; dotted line at '
            + f'{R_GOOD:.2f}.',
            'Control chromosomes on an allele the reference does not list, on b\'s rows.',
        ],
        notes=(f'{n_ctrl // 2:,} controls, HLA-HD against IPD-IMGT/HLA 3.64.0, versus '
               f'{args.reference_name}, up to {n_ref // 2:,} individuals '
               f'({args.reference_cite}); alleles matched on IPD-IMGT P groups. '
               + (f'† {flagged}: does not reconcile, see OPEN_QUESTIONS. ' if flagged else '')
               + f'POPULATION FREQUENCIES, not per-sample accuracy.'),
        top_pad=0.30, wspace=0.30, hspace=0.40,
        left='auto', right=0.975, margin_axes=[ax_p, ax_r, ax_u])
    fig.savefig(args.out_png)
    plt.close(fig)

    worst = summ.loc[r_all.idxmin()] if len(summ) else None
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
             f'frequency — a normal-approximation interval goes below zero for the rare '
             f'alleles that make up most of an HLA table. Named in red: the two commonest '
             f'and the four furthest from the line. The per-locus version of this panel is '
             f'{Path(args.out_png_loci).name}.'),
            ('b', 'Per-locus concordance',
             f'Pearson $r$ and Spearman $\\rho$ between the two frequency vectors, per '
             f'locus, with a dotted line at {R_GOOD:.2f}. $r$ is dominated by the common '
             f'alleles and $\\rho$ weights every allele equally, so the two disagreeing '
             f'means the disagreement is in the tail — which is where a finite reference '
             f'is least informative.'),
            ('c', 'Unconfirmed share by locus',
             'The fraction of control chromosomes carrying an allele the reference does '
             'not list, per locus, drawn on panel (b)\'s rows so each locus\'s concordance '
             'and its unconfirmed share are read together. Derived from the same table as '
             '(a) and (b): pooled over loci it reproduces the cohort figure in '
             'typing_confound.png exactly.'),
        ],
        numbers={
            'loci compared': f'{len(summ)}',
            f'loci with r >= {R_GOOD:.2f}': f'{n_good}',
            'control samples': f'{n_ctrl // 2:,}',
            'reference individuals': f'{n_ref // 2:,}',
            'lowest r': (f'{float(worst.pearson_r):.4f} ({worst.gene})'
                         if worst is not None else 'n/a'),
            'highest r': f'{r_all.max():.4f}',
            'worst unconfirmed locus': (f'{unc["share"].idxmax()} '
                                        f'{unc["share"].max():.1%}'),
            'best unconfirmed locus': (f'{unc["share"].idxmin()} '
                                       f'{unc["share"].min():.1%}'),
            'alleles compared': f'{len(chk):,}',
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
            'reflects the reference\'s size. Panels (b) and (c) are deliberately on the '
            'same rows: a locus can have high concordance and still put a large share of '
            'chromosomes on alleles the reference never lists, and the two together say '
            'which of those is happening.'),
        limits=[
            'This is a comparison of POPULATION FREQUENCIES. It cannot say that any '
            'given sample was typed correctly, and a set of errors that happens to '
            'preserve the frequency spectrum is invisible to it. Only typing samples '
            'with known types measures accuracy — see docs/OPEN_QUESTIONS.md §2.',
            f'Only the {len(summ)} loci the reference carries can be checked at all. Every '
            f'other locus this component types has no external reference here, so its '
            f'calls are unverified rather than verified-and-passing.',
            f'The reference bounds error at the resolution its size allows: with up to '
            f'{n_ref // 2:,} individuals the standard error on a frequency near 0.4 is '
            f'about {0.5 / np.sqrt(max(n_ref, 1)):.4f}.',
            'Cases are excluded from the comparison and shown only in the table. They '
            'are a disease series, and the MHC is where a disease series is least '
            'expected to match a population reference.',
        ],
        methods_ref='../../docs/METHODS.md')

    # =====================================================================
    # FIGURE 2 -- the record
    # =====================================================================
    ncol = 3 if len(genes) <= 6 else 4
    nrow = -(-len(genes) // ncol)
    loci_h = PANEL_H * nrow
    fig2, axes2 = plt.subplots(nrow, ncol, figsize=(S.COL_DOUBLE, loci_h), squeeze=False)
    flat = axes2.ravel()
    for i, (ax, g) in enumerate(zip(flat, genes)):
        locus_panel(ax, chk[chk.gene == g].copy(), g, string.ascii_lowercase[i],
                    xlab=(i + ncol >= len(genes)),   # nothing plotted below it
                    ylab=(i % ncol == 0),            # leftmost column
                    flag=(g in flag_set))
    for ax in flat[len(genes):]:
        ax.axis('off')

    S.caption_block(
        fig2, plot_h=loci_h,
        # plot_style's default is 'abcdefgh' -- eight. Passed explicitly here rather
        # than widening the shared default, which every other component's figures use.
        letters=string.ascii_lowercase,
        title=(f'The per-locus record behind {Path(args.out_png).name}: observed control '
               f'frequency against {args.reference_name}, one panel per locus.'),
        notes=(f'Each locus is named on its own panel and scaled to its own range, so '
               f'panels are not comparable by eye; points on the dashed line agree. Bars '
               f'are 95 % Wilson intervals; alleles in red are the commonest plus any '
               f'differing by {MIN_DIFF_TO_NAME:.2f} or more. '
               + (f'† {flagged}: does not reconcile. ' if flagged else '')
               + f'{n_ctrl // 2:,} controls versus up to {n_ref // 2:,} reference '
                 f'individuals.'),
        top_pad=0.26, wspace=0.42, hspace=0.52,
        left='auto', right=0.975, margin_axes=list(flat))
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
                 f'line is equality. Error bars are 95 % Wilson intervals on the observed '
                 f'frequency. Named in red: the commonest alleles and the two furthest '
                 f'from the line. The axis limits are set from this locus alone, so a '
                 f'panel cannot be compared with its neighbour by eye — that is what the '
                 f'concordance panel of {Path(args.out_png).name} is for.')
                for i, g in enumerate(genes)],
        numbers={f'r, {r.gene}': f'{float(r.pearson_r):.4f}' for _, r in summ.iterrows()},
        interpretation=(
            'This figure is the evidence, not the argument; the argument is in '
            f'{Path(args.out_png).name}. Read a panel when a locus matters to a specific '
            'question — which allele is off the line, and by how much relative to its '
            'interval. A panel whose points sit on the line at every frequency is a '
            'locus whose common haplotypes are being read correctly.'),
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
