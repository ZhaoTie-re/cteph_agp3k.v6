#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : What the GRM is actually accounting for — the relatedness among the
#           analysed samples, measured on the very markers the GRM is built from.
#             (a) the GRM off-diagonal distribution, and the diagonal beside it
#             (b) related pairs by degree
#             (c) what retaining them is worth: samples the fixed-effects design
#                 had to drop, which this one keeps
#
#           This figure exists because "the mixed model lets us keep relateds" is
#           the central design claim of the component, and a claim that is never
#           measured is only an assertion. Panel (c) is the claim as a number.
#
#           THE METRIC IS THE GRM ITSELF, computed on the same markers, so the
#           figure reports the values the model used rather than a proxy. The
#           diagonal gets its own panel because SAIGE overrides it
#           (--isDiagofKinSetAsOne=True) and what it was is otherwise lost.
#
#           Caveat carried in the sidecar: a GRM value conflates relatedness with
#           shared ancestry, so in ancestry-filtered cohorts the bulk shifts with
#           the filter and the degree bands are indicative, not diagnostic.
# Component: assoc_saige
# Used by : assoc_saige.nf  process COMPARE_RELATEDNESS
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
SHARED = Path(__file__).resolve().parent.parent.parent / '_shared' / 'scripts'
sys.path.insert(0, str(SHARED))
import plot_style as S
import figure_doc

PLOT_H = 5.2
# Expected GRM relationship per degree, weakest first so the bars read left to
# right in increasing relatedness. Twice the KING kinship boundaries.
DEG_ORDER = [('third', '3rd degree', 0.0884),
             ('second', '2nd degree', 0.177),
             ('first', '1st degree', 0.354),
             ('mz_or_dup', 'duplicate / MZ', 0.708)]


def parse_args():
    p = argparse.ArgumentParser(description='Relatedness among the analysed samples.')
    p.add_argument('--summary', action='append', required=True, metavar='COHORT=PATH')
    p.add_argument('--pairs', action='append', default=[], metavar='COHORT=PATH')
    p.add_argument('--hist', action='append', default=[], metavar='COHORT=PATH')
    p.add_argument('--diag', action='append', default=[], metavar='COHORT=PATH')
    p.add_argument('--cohort-order', default='')
    p.add_argument('--min-relationship', type=float, default=0.0884)
    p.add_argument('--out-png', required=True)
    p.add_argument('--out-table', required=True)
    return p.parse_args()


def spec(items):
    out = {}
    for s in items or []:
        c, _, path = s.partition('=')
        if path and Path(path).exists() and Path(path).stat().st_size:
            out[c] = path
    return out


def main():
    args = parse_args()
    S.setup_style('paper')
    summ, prs, hst = spec(args.summary), spec(args.pairs), spec(args.hist)
    dia = spec(args.diag)
    cohorts = [c for c in args.cohort_order.split(',') if c in summ] or list(summ)
    SHORT = S.shorten(cohorts)
    ccol = {c: S.COHORT_RAMP[i % len(S.COHORT_RAMP)] for i, c in enumerate(cohorts)}

    tab = pd.concat([pd.read_csv(summ[c], sep='\t') for c in cohorts], ignore_index=True)
    tab.to_csv(args.out_table, sep='\t', index=False)
    row = {r['cohort']: r for _, r in tab.iterrows()}

    fig = plt.figure(figsize=(S.COL_DOUBLE, PLOT_H))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0], width_ratios=[2.0, 1.0, 1.0])
    ax_d = fig.add_subplot(gs[0, :2])      # off-diagonal: relatedness between samples
    ax_g = fig.add_subplot(gs[0, 2])       # diagonal: relatedness to self
    ax_b = fig.add_subplot(gs[1, 0])
    ax_k = fig.add_subplot(gs[1, 1:])

    # ── (a) the whole pairwise distribution ────────────────────────────────
    for c in cohorts:
        if c not in hst:
            continue
        h = pd.read_csv(hst[c], sep='\t')
        mid = (h['lo'] + h['hi']) / 2
        ax_d.step(mid, h['n_pairs'].clip(lower=0.5), where='mid', lw=1.2,
                  color=ccol[c], label=SHORT[c], zorder=4)
    ax_d.set_yscale('log')
    ax_d.set_xlabel('GRM relationship (off-diagonal)')
    ax_d.set_ylabel('sample pairs')
    # The degree boundaries, so the tail is read against a standard rather than
    # by eye. Only the ones inside the plotted range are drawn.
    lo_x, hi_x = ax_d.get_xlim()
    for _key, label, thr in DEG_ORDER:
        if lo_x < thr < hi_x:
            ax_d.axvline(thr, color=S.REFERENCE, lw=0.8, ls=':', zorder=2)
            ax_d.annotate(label, xy=(thr, 1.0), xycoords=('data', 'axes fraction'),
                          xytext=(2, -3), textcoords='offset points', rotation=90,
                          ha='left', va='top', color=S.REFERENCE,
                          fontsize=plt.rcParams['legend.fontsize'] - 1.5)
    S.despine(ax_d, grid_axis='y')
    S.legend_above(ax_d, [Line2D([], [], color=ccol[c], lw=1.6, label=SHORT[c])
                          for c in cohorts], ncol=len(cohorts), y=1.04)

    # ── (a, right) the diagonal: each sample's relatedness to itself ───────
    # Reported because SAIGE replaces it with exactly 1.0 when it fits the null
    # (--isDiagofKinSetAsOne=True), so this is the only record of what the data
    # said. A diagonal far from 1 is inbreeding or a genotyping problem.
    any_diag = False
    for c in cohorts:
        if c not in dia:
            continue
        h = pd.read_csv(dia[c], sep='\t')
        if not len(h) or not h['n_samples'].sum():
            continue
        any_diag = True
        mid = (h['lo'] + h['hi']) / 2
        ax_g.step(mid, h['n_samples'].clip(lower=0.5), where='mid', lw=1.2,
                  color=ccol[c], zorder=4)
    ax_g.axvline(1.0, color=S.ACCENT, lw=1.0, ls='--', zorder=3)
    ax_g.set_yscale('log')
    ax_g.set_xlabel('GRM diagonal')
    ax_g.set_ylabel('samples')
    if not any_diag:
        ax_g.text(0.5, 0.5, 'no diagonal in this matrix', transform=ax_g.transAxes,
                  ha='center', va='center', color=S.REFERENCE,
                  fontsize=plt.rcParams['legend.fontsize'] - 1)
    S.despine(ax_g, grid_axis='y')

    # ── (b) related pairs by degree ────────────────────────────────────────
    keys = [k for k, _l, _t in DEG_ORDER]
    xs = np.arange(len(keys))
    w = 0.8 / max(len(cohorts), 1)
    for i, c in enumerate(cohorts):
        v = [int(row[c].get(f'n_pairs_{k}', 0)) for k in keys]
        ax_b.bar(xs + (i - (len(cohorts) - 1) / 2) * w, v, width=w * 0.9,
                 color=ccol[c], edgecolor='white', linewidth=0.6, zorder=3)
    ax_b.set_xticks(xs)
    ax_b.set_xticklabels([l for _k, l, _t in DEG_ORDER],
                         fontsize=plt.rcParams['xtick.labelsize'] - 1.5, rotation=20, ha='right')
    ax_b.set_ylabel('related pairs')
    S.despine(ax_b, grid_axis='y')

    # ── (c) what keeping them is worth ─────────────────────────────────────
    xs2 = np.arange(len(cohorts))
    have = [int(row[c]['n_samples_with_relative']) for c in cohorts]
    drop = [int(row[c].get('n_absent_from_fixed_model', 0)) for c in cohorts]
    ax_k.bar(xs2 - 0.19, have, width=0.36, color=[ccol[c] for c in cohorts],
             edgecolor='white', linewidth=0.6, zorder=3)
    ax_k.bar(xs2 + 0.19, drop, width=0.36, facecolor='white', linewidth=1.2,
             edgecolor=[ccol[c] for c in cohorts], zorder=3)
    S.value_labels(ax_k, xs2 - 0.19, have, [f'{v:,}' for v in have], offset=5)
    S.value_labels(ax_k, xs2 + 0.19, drop, [f'{v:,}' for v in drop], offset=5)
    ax_k.set_xticks(xs2)
    ax_k.set_xticklabels([SHORT[c] for c in cohorts],
                         fontsize=plt.rcParams['xtick.labelsize'] - 1)
    ax_k.set_ylabel('samples')
    ax_k.set_xlim(-0.6, len(cohorts) - 0.4)
    S.despine(ax_k, grid_axis='y')
    S.legend_above(ax_k, [
        Line2D([], [], marker='s', ls='none', markersize=7, markerfacecolor=S.NEUTRAL_D,
               markeredgecolor=S.NEUTRAL_D, label='has a relative in the analysis'),
        Line2D([], [], marker='s', ls='none', markersize=7, markerfacecolor='white',
               markeredgecolor=S.NEUTRAL_D, markeredgewidth=1.2,
               label='absent from the fixed-effects set')], ncol=1, y=1.04)

    tot_pairs = int(tab['n_related_pairs'].sum())
    tot_kept = int(tab.get('n_absent_from_fixed_model', pd.Series([0])).sum())
    fig.suptitle('Relatedness among the analysed samples',
                 fontsize=plt.rcParams['axes.titlesize'] + 1, fontweight='bold', y=0.988)
    S.caption_block(
        fig,
        title=(f'{tot_pairs} related pairs at GRM relationship $\\geq$ '
               f'{args.min_relationship:g} across {len(cohorts)} sample sets; {tot_kept} samples '
               f'are retained here that a relatedness-pruned design removes.'),
        panels=[
            ('GRM off-diagonal over every pair, log count. Dotted lines are the expected values '
             'for each degree of relationship.'),
            ('The GRM diagonal — each sample against itself. Red dashed line is 1.0, the value '
             'SAIGE substitutes when it fits the null.'),
            'Related pairs by degree.',
            ('Samples with at least one relative in the analysis, against those absent from the '
             'fixed-effects sample set.'),
        ],
        notes=('The matrix is the variance-standardized relationship matrix computed on exactly '
               'the LD-pruned markers the GRM is built from, so these are the values the null '
               'model used. A GRM relationship is twice a kinship coefficient. Full explanation: '
               + Path(args.out_png).stem + '.md'),
        plot_h=PLOT_H, top_pad=0.62, hspace=0.75, wspace=0.34, left='auto', right=0.985,
        margin_axes=[ax_d, ax_g, ax_b, ax_k])
    for ax, letter in ((ax_d, 'a'), (ax_g, 'b'), (ax_b, 'c'), (ax_k, 'd')):
        S.panel_tag(ax, letter, pad=S.strip_pad(ax))
    fig.savefig(args.out_png)
    plt.close(fig)

    allp = pd.concat([pd.read_csv(prs[c], sep='\t') for c in cohorts if c in prs],
                     ignore_index=True) if prs else pd.DataFrame()
    figure_doc.write_doc(
        args.out_png,
        title='Relatedness among the analysed samples',
        question=('What relatedness does the GRM encode among the samples this component analyses, '
                  'and what does retaining it buy over a design that prunes it away?'),
        interpretation=(
            'The mixed model exists so that related samples can be KEPT: the genetic relationship '
            'matrix accounts for the correlation between them instead of the design discarding one '
            'of every related pair. This figure measures that. Panel (a) shows the whole pairwise '
            'distribution, whose bulk sits near zero — the overwhelming majority of pairs are '
            'unrelated, and the informative content is the tail beyond the third-degree boundary. '
            'Panel (b) resolves that tail into degree classes. Panel (c) is the design claim as a '
            'number: how many samples carry a relative, and how many are present here but absent '
            'from the relatedness-pruned set the fixed-effects component analyses. Relatedness is '
            'accounted for, not removed, so these samples contribute to every estimate. Note that '
            'kinship is computed on the LD-pruned GRM markers, which is the matrix the null model '
            'actually fitted, not a separately ascertained marker set.'),
        panels=[
            ('a', 'GRM off-diagonal, whole distribution',
             'The variance-standardized relationship between every pair of analysed samples, '
             'binned, on a log count axis so the tail is visible against a bulk that is orders of '
             'magnitude larger. This is the matrix the null model fitted, computed on the same '
             'LD-pruned markers. Dotted vertical lines are the expected values per degree — '
             '0.0884, 0.177, 0.354, 0.708 — which are twice the corresponding KING kinship '
             'thresholds, because a GRM relationship is 2 x kinship. Negative values are expected: '
             'the estimator is not bounded below at zero and returns small negative values for '
             'pairs less related than the sample average.'),
            ('b', 'GRM diagonal',
             "Each sample's relationship with itself, which is 1 + inbreeding. It is reported "
             'because SAIGE replaces it with exactly 1.0 when fitting the null '
             '(--isDiagofKinSetAsOne=True), so this is the only record of what the data actually '
             'said. A diagonal far from 1 indicates inbreeding or a genotyping problem.'),
            ('c', 'Related pairs by degree',
             'Counts of pairs in each degree class, per cohort. A pair is assigned the strongest '
             'class whose threshold it clears.'),
            ('d', 'What retaining relatedness is worth',
             'Filled bars: samples with at least one relative in the analysis. Open bars: samples '
             'present in this component\'s sample set but absent from the fixed-effects one, which '
             'prunes relateds. The open bar is the direct cost the fixed-effects design pays and '
             'this one does not.'),
        ],
        numbers=[(f'{r["cohort"]}: {k}', r[k]) for _, r in tab.iterrows()
                 for k in ('n_samples', 'n_related_pairs', 'n_samples_with_relative',
                           'n_cases_with_relative', 'kinship_max', 'n_absent_from_fixed_model')
                 if k in tab.columns],
        tables=[('Every related pair', allp if len(allp) else None, 40)],
        reading=[
            'Read panel (a) on the log axis: the bulk near zero is every unrelated pair, and the '
            'question is only how far the right tail extends.',
            'Check panel (b) against 1.0. A systematic shift there is a property of the marker set '
            'or the samples, not of any one pair, and it is invisible in the fitted model because '
            'SAIGE overwrites the diagonal.',
            'A pair beyond the 1st-degree line is a parent-offspring or full-sibling pair. A '
            'duplicate or MZ twin would sit near 1.0 and would be worth investigating as a '
            'sample-handling question, not a biological one.',
            'In panel (d), compare the two bars within a cohort, not across cohorts.',
            'If cases carry relatives, the relatedness is not independent of phenotype and the GRM '
            'is doing more than nuisance correction — check n_cases_with_relative in the table.',
        ],
        limits=[
            'A GRM value conflates relatedness with SHARED ANCESTRY. These cohorts are '
            'ancestry-filtered subsets of one another, so the off-diagonal bulk shifts with the '
            'filter and the degree bands are indicative rather than diagnostic. A kinship '
            'estimator such as KING-robust is designed to be insensitive to structure and would '
            'classify degrees more reliably; what this figure gains in exchange is that it reports '
            'the values the model actually used.',
            'It says nothing about whether the GRM has absorbed the relatedness adequately. That '
            'question is answered by the calibration read-out, not here.',
            'The diagonal shown is the computed one. The fitted model does not use it.',
            'Degree classes are threshold calls on a noisy estimator. A pair near a boundary may '
            'belong to either class, and the counts should be read as approximate.',
        ],
        defs=['neff'],
        methods_ref='../../../docs/METHODS.md')
    print(f'[plot_relatedness] {len(cohorts)} cohorts, {tot_pairs} related pairs -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in plot_relatedness: {e}', file=sys.stderr)
        sys.exit(1)
