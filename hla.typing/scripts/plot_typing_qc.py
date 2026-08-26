#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : One figure answering the question this component cannot answer in
#           prose: does typing quality differ by sequencing platform?
#
#           It has to be read by platform, not by case/control. Cases and controls
#           were sequenced differently and with no overlap, so any difference
#           between the groups is confounded by construction. Platform is the axis
#           the technical gradient actually runs along, and if the metrics are flat
#           across it then read length and depth are not driving the typing.
# Component: hla.typing
# ---------------------------------------------------------------------------
import argparse
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

# 7.6, not 5.4: the two bar panels carry one row per gene, and 19 rows do not
# fit a 2.2 in axis — the labels overlapped into an unreadable band. The row
# count is data-dependent, so the height is set from it below.
PLOT_H_BASE = 3.2     # the two scatter panels, now 5 horizontal rows each
ROW_IN = 0.24         # inches per gene row; a 7 pt label plus its padding
                      # needs ~0.16 in, and crowding them was the defect this
                      # number fixes. 19 genes -> a 4.6 in bar block.


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample-qc', required=True)
    p.add_argument('--platform-qc', required=True)
    p.add_argument('--gene-qc', required=True)
    p.add_argument('--sites', required=True)
    p.add_argument('--out-png', required=True)
    return p.parse_args()


def order_platforms(df):
    """Platforms by sample count, largest first — the controls dominate and
    should anchor the eye."""
    return df.groupby('platform').size().sort_values(ascending=False).index.tolist()


def main():
    args = parse_args()
    S.setup_style('slide')

    sq = pd.read_csv(args.sample_qc, sep='\t')
    gq = pd.read_csv(args.gene_qc, sep='\t')
    sites = pd.read_csv(args.sites, sep='\t')
    if 'platform' not in sq.columns:
        raise SystemExit(f'ABORT: {args.sample_qc} has no platform column')

    plats = order_platforms(sq)
    colors = {p: S.SERIES[i % len(S.SERIES)] for i, p in enumerate(plats)}
    # Samples on one platform share few distinct values, so without jitter a
    # column of ten points draws as three and the spread is invisible. Seeded, so
    # the same data always renders identically.
    rng = np.random.default_rng(0)

    def jitter(n):
        return rng.uniform(-0.16, 0.16, n)

    n_rows = max(len(gq), int((sites.polymorphic > 0).sum()), 1)
    bar_h = max(1.6, n_rows * ROW_IN)
    plot_h = PLOT_H_BASE + bar_h
    fig, axes = plt.subplots(2, 2, figsize=(S.COL_DOUBLE, plot_h),
                             gridspec_kw={'height_ratios': [PLOT_H_BASE - 0.9, bar_h]})
    ax_call, ax_depth, ax_gene, ax_amb = axes.ravel()

    # (a) and (b) are HORIZONTAL, platform on y. The platform names are up to 17
    # characters ("DNBSeq-G400RS 30x") and five of them on a shared x-axis leave
    # 0.64 in each, which is not enough for a third of the name — they collided
    # into one unreadable smear. On y they get the left margin, which
    # fit_left_margin sizes to them, and the orientation matches (c) and (d).
    for ax, metric, label in [(ax_call, 'call_rate', 'call rate'),
                              (ax_depth, 'mean_field_depth', 'mean field depth')]:
        for i, pf in enumerate(plats):
            v = sq.loc[sq.platform == pf, metric].dropna()
            ax.scatter(v, i + jitter(len(v)), s=11, alpha=0.55, color=colors[pf],
                       edgecolors='none', rasterized=True)
            if len(v):
                ax.plot([v.median()] * 2, [i - 0.28, i + 0.28],
                        color=S.INK, lw=1.6, zorder=5)
        ax.set_yticks(range(len(plats)))
        ax.set_yticklabels(plats)
        ax.set_ylim(len(plats) - 0.5, -0.5)          # first platform at the top
        ax.set_xlabel(label)
        S.despine(ax)
    S.panel_tag(ax_call, 'a')
    S.panel_tag(ax_depth, 'b')
    ax_depth.set_yticklabels([])                     # (a) already names them

    # (c) call rate per gene
    gq = gq.sort_values('call_rate', ascending=True)
    ax_gene.barh(range(len(gq)), gq['call_rate'], color=S.DATA, height=0.72)
    ax_gene.set_yticks(range(len(gq)))
    ax_gene.set_yticklabels(gq['gene'], fontstyle='italic')
    ax_gene.set_xlabel('call rate')
    ax_gene.set_xlim(0, 1.02)
    S.panel_tag(ax_gene, 'c')
    S.despine(ax_gene)

    # (d) polymorphic residue positions per gene
    sp = sites[sites.polymorphic > 0].sort_values('polymorphic')
    ax_amb.barh(range(len(sp)), sp['polymorphic'], color=S.DATA_DARK, height=0.72)
    ax_amb.set_yticks(range(len(sp)))
    ax_amb.set_yticklabels(sp['gene'], fontstyle='italic')
    ax_amb.set_xlabel('polymorphic residue positions')
    S.panel_tag(ax_amb, 'd')
    S.despine(ax_amb)

    n_fail = int(sq['failed'].sum()) if 'failed' in sq.columns else 0
    S.caption_block(
        fig,
        title=(f'{len(sq):,} samples typed; call rate and field resolution are '
               f'compared across the {len(plats)} sequencing platforms, not across '
               f'case and control.'),
        panels=[
            'Per-sample call rate, one column per platform; the bar is the median.',
            'Mean field depth achieved — 2, 3 or 4 fields — by platform.',
            'Call rate per gene over all samples.',
            'Residue positions polymorphic in this cohort, per gene.',
        ],
        notes=(f'HLA-HD against a pinned IPD-IMGT/HLA 3.64.0 dictionary; residues '
               f'from IMGT protein alignments at IMGT numbering. '
               f'{n_fail} (sample, gene) typing failure(s). Platform is the axis to '
               f'read: it is fully confounded with case/control in this study, so a '
               f'case/control difference in these metrics cannot be interpreted.'),
        plot_h=plot_h, top_pad=0.30, wspace=0.46, hspace=0.52,
        left='auto', right=0.955, margin_axes=list(axes.ravel()))
    fig.savefig(args.out_png)
    plt.close(fig)

    figure_doc.write_doc(
        args.out_png,
        title='HLA typing quality, by sequencing platform',
        question='Does typing quality differ between the platforms, and therefore '
                 'between cases and controls, in a way that could imitate a finding?',
        panels=[
            ('a', 'Call rate by platform',
             'The fraction of (sample, gene) pairs that produced a genotype, per '
             'sample, grouped by platform. A platform whose distribution sits lower '
             'is producing fewer calls, and since platform is confounded with '
             'phenotype here that would be indistinguishable from a real difference.'),
            ('b', 'Field resolution by platform',
             'The mean number of IMGT fields resolved per call. Shorter reads carry '
             'less information for separating similar alleles, so a read-length '
             'effect appears here before it appears in the call rate.'),
            ('c', 'Call rate by gene',
             'Which loci are typed reliably. The pseudogenes and the DRB paralogues '
             'are expected to be low; the classical loci are not.'),
            ('d', 'Polymorphic residue positions by gene',
             'How many IMGT positions actually vary in this cohort — the number of '
             'testable sites each gene contributes.'),
        ],
        interpretation=(
            'Read this by platform. Cases were sequenced on DNBSeq and NovaSeq, '
            'controls entirely on HiSeqX, with no overlap, so any case/control '
            'difference in typing quality is confounded by construction and cannot '
            'be attributed. Platform is the axis on which a technical gradient is '
            'visible as itself. Flat metrics across platforms are evidence that read '
            'length and depth are not driving the typing; a gradient is the size of '
            'the problem, stated.'),
        limits=[
            'This figure says nothing about typing ACCURACY. It measures whether a '
            'call was produced and how deeply it resolved, not whether it is right. '
            'Accuracy needs an external truth set.',
            'A flat call rate does not rule out a systematic bias toward particular '
            'alleles, which would need the same truth set to detect.',
        ],
        methods_ref='../../docs/METHODS.md')

    print(f'[plot_typing_qc] {len(sq)} samples over {len(plats)} platform(s) -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in plot_typing_qc: {e}', file=sys.stderr)
        sys.exit(1)
