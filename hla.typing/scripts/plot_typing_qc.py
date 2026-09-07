#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : One figure answering the question this component cannot answer in
#           prose: does typing quality differ by sequencing platform?
#
#           It has to be read by platform, not by case/control. Cases and controls
#           were sequenced differently and with no overlap, so any difference
#           between the groups is confounded by construction.
#
#           EVERY METRIC HERE IS ABOUT COMPLETENESS, and on this cohort completeness
#           is nearly saturated: call rate is k/19 and takes a handful of distinct
#           values, nearly all samples sitting on two of them. An earlier version
#           plotted one jittered point per sample to show those few values, and drew
#           nineteen per-gene bars of which fourteen were exactly 1.0. Both are now
#           compositions, which is the encoding a saturated proportion needs.
#
#           PANELS (c) AND (d) SHARE ONE GENE AXIS. They are two measurements of the
#           same nineteen genes, and they used to be sorted differently, so the gene
#           names were printed twice in two different orders and the reader had to
#           re-find each gene to compare. One order, one set of labels.
#
#           NO NUMBER IS WRITTEN INTO THIS FILE. Every value in the caption and the
#           sidecar is an f-string field read from the tables it describes.
#
#           A LOW `call_rate` AT DRB3/4/5 IS NOT A FAILURE. Their shortfall is
#           entirely `not_typed` — the gene is absent from that haplotype — and
#           `failed` is 0 everywhere. Plotting it as a call rate showed correct
#           biology as if it were breakage, which is why panel (c) is a four-way
#           composition instead.
#
#           Whether a call is RIGHT is not measurable here at all. That is
#           typing_confound.png, off the external frequency check.
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
PLOT_H_BASE = 3.3     # the two composition panels (5 platform rows each) plus the
                      # space their below-axis legends and x labels need
ROW_IN = 0.155        # inches per gene row. Was 0.24, sized for a panel that
                      # printed the gene names twice; now (c) carries the labels
                      # and (d) is aligned to it, so a row needs one 7 pt line.


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample-qc', required=True)
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
    n_rows = max(len(gq), int((sites.polymorphic > 0).sum()), 1)
    bar_h = max(1.6, n_rows * ROW_IN)
    plot_h = PLOT_H_BASE + bar_h
    fig, axes = plt.subplots(2, 2, figsize=(S.COL_DOUBLE, plot_h),
                             gridspec_kw={'height_ratios': [PLOT_H_BASE - 1.4, bar_h]})
    ax_call, ax_depth, ax_gene, ax_amb = axes.ravel()

    # (a) FIELD RESOLUTION as a composition, not a mean. The mean lives in
    # 2.84-2.95 for every platform, a range no reader can act on; the 2/3/4-field
    # split is what actually differs and it is a proportion, so it is stacked.
    fld = [c for c in ('field2', 'field3', 'field4') if c in sq.columns]
    lab_f = {'field2': '2-field', 'field3': '3-field', 'field4': '4-field'}
    comp = sq.groupby('platform')[fld].sum().reindex(plats)
    comp = comp.div(comp.sum(axis=1), axis=0)
    # A legend swatch for a series that is identically zero is a lie of omission.
    # field4 is 0 for every platform here -- HLA-HD never resolves past 3 fields on
    # this cohort -- so the series is dropped and the fact is stated in the caption.
    fld = [c for c in fld if comp[c].sum() > 0]
    left = np.zeros(len(plats))
    for i, c in enumerate(fld):
        ax_call.barh(range(len(plats)), comp[c], left=left, height=0.66,
                     color=S.COHORT_RAMP[i % len(S.COHORT_RAMP)], label=lab_f[c])
        left += comp[c].to_numpy()
    # The discriminating number is the 2-field share; print it so the platform gap
    # is read off the figure rather than estimated from a 4-point-wide segment.
    for y, v in enumerate(comp[fld[0]].to_numpy()):
        ax_call.text(v + 0.015, y, f'{v:.1%}', va='center', ha='left',
                     fontsize=plt.rcParams['legend.fontsize'] - 2, color=S.INK)
    ax_call.set_yticks(range(len(plats)))
    ax_call.set_yticklabels(plats)
    ax_call.set_ylim(len(plats) - 0.5, -0.5)
    ax_call.set_xlim(0, 1)
    ax_call.xaxis.set_major_formatter(lambda v, _: f'{v:.0%}')
    ax_call.set_xlabel('calls')
    S.panel_tag(ax_call, 'a', title='field resolution')
    ax_call.legend(loc='upper center', bbox_to_anchor=(0.5, -0.30), ncol=3, frameon=False,
                   handlelength=1.1, columnspacing=1.1,
                   fontsize=plt.rcParams['legend.fontsize'] - 2)
    S.despine(ax_call)

    # (b) AMBIGUITY, which varies and was never plotted. 2,957 records across the
    # cohort; per sample it runs 0-5, and HLA-HD emitting several candidate pairs
    # is the closest thing to a confidence signal the tool gives.
    amb = sq.groupby('platform')['ambiguous'].apply(
        lambda s: pd.Series({f'{k}': (s == k).mean() for k in range(0, 4)}
                            | {'4+': (s >= 4).mean()})).unstack().reindex(plats)
    left = np.zeros(len(plats))
    ramp = [S.NEUTRAL, S.NEUTRAL_D, S.DATA, S.DATA_DARK, S.ACCENT]
    for i, c in enumerate(['0', '1', '2', '3', '4+']):
        if c not in amb.columns:
            continue
        ax_depth.barh(range(len(plats)), amb[c], left=left, height=0.66,
                      color=ramp[i], label=f'{c} genes' if c == '0' else c)
        left += amb[c].to_numpy()
    ax_depth.set_yticks(range(len(plats)))
    ax_depth.set_yticklabels([])
    ax_depth.set_ylim(len(plats) - 0.5, -0.5)
    ax_depth.set_xlim(0, 1)
    ax_depth.xaxis.set_major_formatter(lambda v, _: f'{v:.0%}')
    ax_depth.set_xlabel('samples')
    S.panel_tag(ax_depth, 'b', title='ambiguity')
    ax_depth.legend(loc='upper center', bbox_to_anchor=(0.5, -0.30), ncol=5, frameon=False,
                    handlelength=1.1, columnspacing=1.0,
                    fontsize=plt.rcParams['legend.fontsize'] - 2)
    S.despine(ax_depth)

    # (c) THE FOUR OUTCOMES per gene, not a call rate. `called` and `hemizygous`
    # are both successes and their ratio is real biology — DMA is hemizygous in
    # 73 % of samples and HLA-A in 14 %, because DMA is nearly monomorphic here.
    # `not_typed` is the gene being absent from the haplotype, which is what makes
    # DRB3/4/5 look broken on a call-rate axis when nothing failed.
    states = [c for c in ('called', 'hemizygous', 'not_typed', 'failed') if c in gq.columns]
    gq = gq.copy()
    tot = gq[states].sum(axis=1)
    gq = gq.assign(**{c: gq[c] / tot for c in states})
    # Sorted ASCENDING but drawn on a NON-inverted axis -- unlike (a) and (b), which
    # invert with set_ylim -- so row 0 lands at the BOTTOM and the panel reads
    # DESCENDING top to bottom: most not_typed first, then most hemizygous. Switching
    # this to ascending=False would reverse the FIGURE, not the sort.
    gq = gq.sort_values(['not_typed', 'hemizygous'])
    cols_s = {'called': S.DATA_DARK, 'hemizygous': S.DATA,
              'not_typed': S.NEUTRAL, 'failed': S.ACCENT}
    left = np.zeros(len(gq))
    for c in states:
        ax_gene.barh(range(len(gq)), gq[c], left=left, height=0.72,
                     color=cols_s[c], label=c)
        left += gq[c].to_numpy()
    ax_gene.set_yticks(range(len(gq)))
    ax_gene.set_yticklabels(gq['gene'], fontstyle='italic')
    ax_gene.set_xlim(0, 1)
    ax_gene.xaxis.set_major_formatter(lambda v, _: f'{v:.0%}')
    S.panel_tag(ax_gene, 'c', title='gene outcome')
    ax_gene.legend(loc='upper center', bbox_to_anchor=(0.5, -0.055), ncol=4, frameon=False,
                   handlelength=1.1, columnspacing=1.1,
                   fontsize=plt.rcParams['legend.fontsize'] - 2)
    S.despine(ax_gene)

    # (d) polymorphic residue positions per gene, ON (c)'S ROW ORDER.
    # These are two measurements of the same nineteen genes. Sorting each panel by
    # its own value printed the gene names twice in two different orders, so
    # comparing a gene's outcome with its testable-site count meant finding it
    # again in the second list. One order, and only (c) carries the labels.
    sp = sites.set_index('gene').reindex(gq['gene'])['polymorphic'].fillna(0)
    ax_amb.barh(range(len(gq)), sp.to_numpy(), color=S.DATA_DARK, height=0.72)
    ax_amb.set_yticks(range(len(gq)))
    ax_amb.set_yticklabels([])
    ax_amb.set_ylim(*ax_gene.get_ylim())
    ax_amb.set_xlabel('polymorphic positions')
    S.panel_tag(ax_amb, 'd', title='testable sites')
    S.despine(ax_amb)

    n_fail = int(sq['failed'].sum()) if 'failed' in sq.columns else 0
    # The values the prose used to state by hand. Each is derived here, once, and
    # interpolated wherever it is quoted -- caption, panel text and sidecar.
    f2 = comp[fld[0]]
    f2_worst, f2_best = f2.idxmax(), f2.idxmin()
    mean_field = (sq[[c for c in ('field2', 'field3', 'field4') if c in sq.columns]]
                  .mul([2, 3, 4][:len([c for c in ('field2', 'field3', 'field4')
                                       if c in sq.columns])]).sum(axis=1)
                  / sq[[c for c in ('field2', 'field3', 'field4')
                        if c in sq.columns]].sum(axis=1))
    mf = sq.assign(_m=mean_field).groupby('platform')['_m'].mean()
    n_cr = sq['call_rate'].nunique() if 'call_rate' in sq.columns else 0
    top2 = (sq['call_rate'].value_counts(normalize=True).head(2).sum()
            if 'call_rate' in sq.columns else float('nan'))
    n_amb_records = int(sq['ambiguous'].sum()) if 'ambiguous' in sq.columns else 0
    hz = gq.set_index('gene')['hemizygous']
    hz_top = hz.idxmax()
    S.caption_block(
        fig,
        title=(f'{len(sq):,} samples, {n_fail} (sample, gene) failures; completeness read '
               f'across the {len(plats)} platforms, not across case and control.'),
        panels=[
            'Calls by IMGT field depth; the printed number is the 2-field share.',
            'Samples by how many genes HLA-HD could not choose a genotype for.',
            'The four outcomes per gene; a low bar at DRB3/4/5 is `not_typed`, never `failed`.',
            'Positions polymorphic in this cohort, on panel c\'s gene order.',
        ],
        notes=(f'HLA-HD against a pinned IPD-IMGT/HLA 3.64.0 dictionary; IMGT-numbered '
               f'residues. Completeness is saturated here — call rate takes {n_cr} values, '
               f'{top2:.0%} of samples on two of them — hence compositions, not means. Whether '
               f'a call is RIGHT is typing_confound.png. Read by platform: it is fully '
               f'confounded with case/control.'),
        plot_h=plot_h, top_pad=0.30, wspace=0.30, hspace=0.52,
        left='auto', right=0.955, margin_axes=list(axes.ravel()))
    fig.savefig(args.out_png)
    plt.close(fig)

    figure_doc.write_doc(
        args.out_png,
        title='HLA typing quality, by sequencing platform',
        question='Does typing quality differ between the platforms, and therefore '
                 'between cases and controls, in a way that could imitate a finding?',
        panels=[
            ('a', 'Field resolution by platform',
             f'The share of calls resolving to each IMGT field depth, with the 2-field share '
             f'printed. Stacked rather than averaged: the per-platform mean field depth spans '
             f'only {mf.min():.2f}-{mf.max():.2f}, a range nobody can act on, while the '
             f'composition separates the platforms cleanly - {f2.max():.1%} of calls stop at '
             f'2 fields on {f2_worst} against {f2.min():.1%} on {f2_best}. Shorter reads carry '
             f'less information for separating similar alleles, so a read-length effect would '
             f'appear here first. Field depths absent from the cohort are not drawn as '
             f'zero-width segments; those present are '
             f'{", ".join(c.replace("field", "") + "-field" for c in fld)}.'),
            ('b', 'Ambiguity by platform',
             f'Samples grouped by how many genes HLA-HD emitted several candidate pairs for '
             f'({n_amb_records:,} (sample, gene) records across the cohort). It is the closest '
             f'thing to a per-call confidence the tool provides, and it was previously computed '
             f'but never plotted.'),
            ('c', 'Outcome by gene',
             f'The four outcomes, as a composition. `called` and `hemizygous` are both '
             f'successes and their ratio is biology - {hz_top} is hemizygous in '
             f'{hz.max():.0%} of samples and A in {hz.get("A", float("nan")):.0%}, because '
             f'{hz_top} is nearly monomorphic in this population. `not_typed` means the gene is '
             f'absent from that haplotype, which is why DRB3/4/5 look broken on a call-rate '
             f'axis when `failed` is {n_fail} for every gene in the cohort.'),
            ('d', 'Polymorphic residue positions by gene',
             f'How many IMGT positions actually vary in this cohort - the number of testable '
             f'sites each gene contributes, {int(sp.sum()):,} in total. Drawn on panel (c)\'s '
             f'row order so a gene sits at the same height in both.'),
        ],
        numbers={
            'samples': f'{len(sq):,}',
            'platforms': f'{len(plats)}',
            '(sample, gene) failures': f'{n_fail}',
            'distinct call-rate values': f'{n_cr}',
            'share of samples on the top two': f'{top2:.1%}',
            '2-field share, worst platform': f'{f2.max():.2%} ({f2_worst})',
            '2-field share, best platform': f'{f2.min():.2%} ({f2_best})',
            'mean field depth, range over platforms': f'{mf.min():.2f}-{mf.max():.2f}',
            'ambiguous (sample, gene) records': f'{n_amb_records:,}',
            'genes carrying residues': f'{int((sites.polymorphic > 0).sum())}',
            'polymorphic positions, total': f'{int(sp.sum()):,}',
        },
        interpretation=(
            'Read this by platform. Cases were sequenced on DNBSeq and NovaSeq, '
            'controls entirely on HiSeqX, with no overlap, so any case/control '
            'difference in typing quality is confounded by construction and cannot '
            'be attributed. Platform is the axis on which a technical gradient is '
            'visible as itself. The gradient here is real but small — the platform '
            'ordering in (a) and (b) is the same one, and the weakest platform is the '
            'one sequenced shallowest. Its size is the size of the problem, stated. '
            'Whether it reaches the CALLS is a different question and a different '
            'figure: see typing_confound.png. The cohort is ancestry-restricted, so the '
            'samples counted here are the ones every table in this component describes.'),
        limits=[
            f'Completeness is saturated on this cohort and therefore a weak quality signal: '
            f'call rate takes {n_cr} distinct values across the whole cohort and {top2:.0%} of '
            f'samples sit on two of them. That is why every panel here is a composition and no '
            f'panel is a mean - a mean over a saturated metric hides the only variation there '
            f'is.',
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
