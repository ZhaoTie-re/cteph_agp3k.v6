#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The external check, drawn: our control allele frequencies against the
#           published Japanese ones.
#
#           typing_qc.png answers "was a call made, and how deeply did it resolve".
#           It cannot answer "is the call right", and neither can any other output of
#           this component. This figure is the closest available substitute — if the
#           typing were systematically mis-assigning alleles, the frequency spectrum
#           would move off the identity line, and no call-rate metric would notice.
#
#           It is a POPULATION comparison, not an accuracy measurement, and the
#           reference is 105 samples. Both facts are on the figure, not only in the
#           document beside it.
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

PANEL_H = 2.7         # per grid ROW at COL_DOUBLE; the grid follows the locus count
# Three labels, not six. An HLA frequency table is a handful of common alleles and a
# long tail at the origin, so naming more than the top few piles the text on top of
# itself exactly where the points are densest. An allele off the line is named too,
# but only if it is off by enough to be worth reading — below that the label is noise.
N_LABEL = 3
MIN_DIFF_TO_NAME = 0.02
N_OUTLIER_MAX = 2     # however many alleles clear MIN_DIFF_TO_NAME, name at most
                      # this many of them. Without the cap, a locus that genuinely
                      # disagrees — or a run with too few controls to estimate
                      # anything — names EVERY allele and the panel becomes a wall
                      # of text exactly where the points are.
MIN_CHR = 20          # control chromosomes below which nothing is drawn. At n=2 the
                      # only observable frequencies are 0, 0.5 and 1, every Wilson
                      # interval spans the axis, and a scatter of that is noise
                      # presented as a comparison.


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--check', required=True, help='allele_frequency_check.tsv')
    p.add_argument('--summary', required=True, help='allele_frequency_summary.tsv')
    p.add_argument('--population', default='JPT',
                   help='1kg_wide only; ignored when --reference-name is given')
    p.add_argument('--reference-name', default='jMorp 61KJPN-HLA',
                   help='what the reference is called on the figure')
    p.add_argument('--reference-cite', default='ToMMo jMorp; Tadaka et al. 2023',
                   help='the citation printed in the notes')
    p.add_argument('--flag-genes', default='DRB3',
                   help='loci whose comparison does not reconcile and must be marked on '
                        'the figure rather than left to a footnote. See '
                        'docs/OPEN_QUESTIONS.md.')
    p.add_argument('--out-png', required=True)
    return p.parse_args()


def to_label(rows):
    """Which alleles to name: the commonest in the reference, plus real outliers."""
    if not len(rows):
        return set()
    named = set(rows.nlargest(min(N_LABEL, len(rows)), 'freq_ref').allele)
    d = pd.to_numeric(rows.diff_ctrl_ref, errors='coerce').abs()
    out = rows.assign(_d=d)[d >= MIN_DIFF_TO_NAME].nlargest(N_OUTLIER_MAX, '_d')
    return named | set(out.allele)


def panel(ax, rows, gene, pop, letter, xlab=False, ylab=False, flag=False):
    """One locus: observed control frequency against the reference, with the line."""
    n_chr = int(rows.n_chr_ctrl.max()) if len(rows) else 0
    if n_chr < MIN_CHR:
        msg = ('no control sample\nin this run' if not n_chr
               else f'only {n_chr} control chromosomes\n— too few to compare')
        ax.text(0.5, 0.5, msg, ha='center', va='center',
                transform=ax.transAxes, color=S.INK_SOFT, fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
        S.panel_tag(ax, letter, title=f'{gene}' + (' †' if flag else ''))
        return

    x, y = rows.freq_ref.to_numpy(float), rows.freq_ctrl.to_numpy(float)
    hi = float(max(np.nanmax(x), np.nanmax(y))) * 1.22 + 0.01
    ax.plot([0, hi], [0, hi], color=S.REFERENCE, lw=1.0, ls='--', zorder=1)

    named = to_label(rows)
    is_out = rows.allele.isin(named).to_numpy()
    # Wilson intervals, so a point far from the line is visibly far relative to how
    # well it is estimated rather than in absolute terms.
    err = np.vstack([y - rows.ci_lo_ctrl.to_numpy(float),
                     rows.ci_hi_ctrl.to_numpy(float) - y])
    ax.errorbar(x, y, yerr=np.clip(err, 0, None), fmt='none',
                ecolor=S.NEUTRAL_D, elinewidth=0.7, zorder=2)
    ax.scatter(x[~is_out], y[~is_out], s=15, color=S.DATA, alpha=0.75,
               edgecolors='none', zorder=3)
    ax.scatter(x[is_out], y[is_out], s=19, color=S.ACCENT, alpha=0.9,
               edgecolors='none', zorder=4)

    annots = [ax.annotate(al.split('*', 1)[-1], (xi, yi), textcoords='offset points',
                          xytext=(5, 3), fontsize=6.5, color=S.INK_SOFT)
              for xi, yi, al in zip(x[is_out], y[is_out], rows.allele[is_out])]

    ax.set_xlim(0, hi); ax.set_ylim(0, hi)
    ax.set_aspect('equal', adjustable='box')
    if annots:
        S.spread_labels(ax, annots, axis='y')
    if xlab:
        # The reference is named in the title and again in the notes; repeating
        # it on four x axes only makes them collide.
        ax.set_xlabel('reference')
    if ylab:
        ax.set_ylabel('observed frequency')
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

    # The grid follows the locus count. It was 2x3 when the reference carried five
    # loci; jMorp carries thirteen, and a hard-coded grid silently drops the rest.
    n_cell = len(genes) + 1                     # + the concordance panel
    ncol = 3 if n_cell <= 6 else 4
    nrow = -(-n_cell // ncol)
    plot_h = PANEL_H * nrow
    fig, axes = plt.subplots(nrow, ncol, figsize=(S.COL_DOUBLE, plot_h), squeeze=False)
    flat = axes.ravel()
    # Loci whose comparison is known not to reconcile are marked ON the panel, so a
    # reader cannot take the scatter at face value and find the caveat only later.
    flag_set = {g.strip() for g in args.flag_genes.split(',') if g.strip()} & set(genes)
    flagged = ' and '.join(f'HLA-{g}' for g in sorted(flag_set))
    for i, (ax, g) in enumerate(zip(flat, genes)):
        panel(ax, chk[chk.gene == g].copy(), g, args.reference_name, chr(ord('a') + i),
              xlab=(i + ncol >= len(genes)),      # nothing plotted below it
              ylab=(i % ncol == 0),               # leftmost column
              flag=(g in flag_set))
    # Immediately after the last locus, not in the grid's final cell: with 13 loci in a
    # 4-wide grid the final cell is two columns away from anything, and the blank gap
    # reads as a missing panel.
    ax_s = flat[len(genes)]
    for ax in flat[len(genes) + 1:]:
        ax.axis('off')

    # The last cell is the per-locus concordance, which is the number a reader wants
    # after looking at five scatters.
    have = (summ[pd.to_numeric(summ.pearson_r, errors='coerce').notna()
                 & (pd.to_numeric(summ.n_chr_ctrl, errors='coerce') >= MIN_CHR)]
            if len(summ) else summ)
    if len(have):
        idx = np.arange(len(have))
        ax_s.barh(idx + 0.19, have.pearson_r.astype(float), height=0.36,
                  color=S.DATA_DARK, label='Pearson $r$')
        ax_s.barh(idx - 0.19, have.spearman_rho.astype(float), height=0.36,
                  color=S.DATA, label=r'Spearman $\rho$')
        ax_s.set_yticks(idx); ax_s.set_yticklabels(have.gene, fontstyle='italic')
        # The bars run to 1.0 when the typing agrees, so there is no free corner
        # inside the panel, and above it is the panel title. Below is the only space
        # left, and this panel has no x label to compete with.
        # rho can be negative here — it weights every allele equally and an HLA
        # table is mostly a tail of ties — so the axis has to admit that rather
        # than clip it to zero and make the bar look absent.
        lo = min(0.0, float(pd.to_numeric(have.spearman_rho, errors='coerce').min()))
        ax_s.set_xlim(min(-0.05, lo * 1.1), 1.02)
        ax_s.axvline(0, color=S.INK_SOFT, lw=0.6, zorder=0)
        ax_s.invert_yaxis()
        ax_s.tick_params(axis='y', labelsize=plt.rcParams['ytick.labelsize'] - 2)
        ax_s.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=1,
                    frameon=False, handlelength=1.1, labelspacing=0.3,
                    fontsize=plt.rcParams['legend.fontsize'] - 2)
    else:
        ax_s.text(0.5, 0.5, f'fewer than {MIN_CHR} control\nchromosomes — nothing\nto compare',
                  ha='center', va='center',
                  transform=ax_s.transAxes, color=S.INK_SOFT, fontsize=8)
        ax_s.set_xticks([]); ax_s.set_yticks([])
    # No centred title: at a quarter of COL_DOUBLE this panel is under the width
    # panel_tag documents as the limit, and 'concordance' runs into the letter.
    S.panel_tag(ax_s, chr(ord('a') + len(genes)))
    S.despine(ax_s)

    n_ctrl = int(chk.n_chr_ctrl.max()) if len(chk) else 0
    n_ref = int(chk.n_chr_ref.max()) if len(chk) else 0
    S.caption_block(
        fig, plot_h=plot_h,
        # plot_style's default is 'abcdefgh' — eight. Thirteen loci plus the
        # concordance panel is fourteen, and the default truncates rather than
        # warning. Passed explicitly here instead of widening the shared default,
        # which would invalidate every other component's figure cache.
        letters=string.ascii_lowercase,
        title=(f'Control allele frequencies against {args.reference_name} at '
               f'{len(genes)} loci: one point per allele, observed against reference, '
               f'and points on the dashed line agree.'),
        panels=[f'HLA-{g}.' for g in genes]
               + [r'Per-locus concordance: Pearson $r$ and Spearman $\rho$.'],
        notes=(f'Error bars are 95 % Wilson intervals on the observed frequency; '
               f'alleles named in red are the commonest plus any differing by '
               f'{MIN_DIFF_TO_NAME:.2f} or more. '
               f'Observed: {n_ctrl // 2} control samples typed by HLA-HD against '
               f'IPD-IMGT/HLA 3.64.0. Reference: {args.reference_name}, up to '
               f'{n_ref // 2:,} individuals ({args.reference_cite}); alleles are matched '
               f'on IPD-IMGT P groups, which is how the reference names them. '
               + (f'† {flagged}: marked because its comparison does not reconcile — the '
                  f'reference assigns every chromosome an allele there and we do not, and '
                  f'no mechanism tested explains the rest of the gap; see OPEN_QUESTIONS. '
                  if flagged else '')
               + f'This compares POPULATION FREQUENCIES and is not a measurement of '
                 f'per-sample accuracy: a set of errors that preserves the spectrum is '
                 f'invisible to it.'),
        top_pad=0.34, wspace=0.42, hspace=0.46,
        left='auto', right=0.965, margin_axes=list(flat))
    fig.savefig(args.out_png)
    plt.close(fig)

    figure_doc.write_doc(
        args.out_png,
        title=f'HLA allele frequencies against {args.reference_name}',
        question='Do the alleles this component calls occur at the frequencies published '
                 'for a Japanese population — or has the typing moved the spectrum?',
        panels=([(chr(ord('a') + i), f'HLA-{g}',
                  f'Every 2-field allele seen at HLA-{g} in the controls or in the '
                  f'reference. x is the {args.population} frequency, y is ours, and the '
                  f'dashed line is equality. Error bars are 95 % Wilson intervals on the '
                  f'observed frequency — a normal-approximation interval goes below zero '
                  f'for the rare alleles that make up most of an HLA table. Named in red: '
                  f'the commonest alleles and the two furthest from the line.')
                 for i, g in enumerate(genes)]
                + [(chr(ord('a') + len(genes)), 'Concordance',
                    'Pearson $r$ and Spearman $\\rho$ between the two frequency vectors '
                    'per locus. $r$ is dominated by the common alleles, $\\rho$ weights '
                    'every allele equally, so the two disagreeing means the disagreement '
                    'is in the tail.')]),
        interpretation=(
            'Agreement here is evidence that the typing is not systematically '
            'mis-assigning alleles — the one thing every other QC output in this '
            'component is blind to, because call rate and field depth measure whether a '
            'call was produced, not whether it was right. Disagreement at a COMMON '
            'allele is the informative failure: it means a frequent haplotype is being '
            'read as something else. Disagreement in the tail is expected and mostly '
            'reflects the reference\'s size.'),
        limits=[
            'This is a comparison of POPULATION FREQUENCIES. It cannot say that any '
            'given sample was typed correctly, and a set of errors that happens to '
            'preserve the frequency spectrum is invisible to it. Only typing samples '
            'with known types measures accuracy — see docs/OPEN_QUESTIONS.md §2.',
            f'The reference is {n_ref // 2} samples. The standard error on a frequency '
            f'near 0.4 is about 0.034, so this bounds gross error and nothing finer.',
            'Five loci only — A, B, C, DQB1, DRB1. The reference carries no DPB1, no '
            'DQA1 and no DRB3/4/5, so the loci with the least reliable typing in this '
            'component are exactly the ones it cannot check.',
            'Cases are excluded from the comparison and shown only in the table. They '
            'are a disease series, and the MHC is where a disease series is least '
            'expected to match a population reference.',
        ],
        methods_ref='../../docs/METHODS.md')

    print(f'[plot_allele_freq] {len(genes)} locus/loci, {n_ctrl // 2} controls '
          f'vs {args.population} -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in plot_allele_freq: {e}', file=sys.stderr)
        sys.exit(1)
