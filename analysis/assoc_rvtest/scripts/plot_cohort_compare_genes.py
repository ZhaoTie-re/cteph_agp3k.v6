#!/usr/bin/env python3
# ==============================================================================
# plot_cohort_compare_genes.py — the three cohorts side by side, at gene level.
#
# The three cohorts are NESTED: narrow subset of intermediate subset of full,
# sharing their cases entirely and most of their controls. Agreement between them
# is therefore NOT replication — a noise hit in the shared core reappears in all
# three by construction. What it does measure is ROBUSTNESS TO THE ANCESTRY
# FILTER, because what the wider cohorts add is precisely the samples the
# stricter filter rejected. Every panel and the caption say so, so the figure
# cannot be read as independent support.
#
#   (a) calibration: gene-level lambda per scan, and how many genes each tested
#   (b) every gene reaching a tier in any cohort, shown in all three
#   (c) in how many cohorts each such gene is a hit
# ==============================================================================
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_shared' / 'scripts'))
import plot_style as S          # noqa: E402
import figure_doc               # noqa: E402

CALLED = {'significant': dict(marker='D', size=34, filled=True, label='significant here'),
          'suggestive':  dict(marker='o', size=28, filled=True, label='suggestive here'),
          'not_a_hit':   dict(marker='o', size=20, filled=False, label='tested, no tier'),
          'not_tested':  dict(marker='x', size=24, filled=False, label='not tested here')}
ROW_IN = 0.170


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--genes', required=True, help='gene_crosscohort.tsv')
    p.add_argument('--scans', required=True, help='scan_qc_all.tsv')
    p.add_argument('--cohort-order', required=True)
    p.add_argument('--peak-stratum', default='moderate_high',
                   help='the stratum panels (b) and (c) show; others stay in the table')
    p.add_argument('--peak-method', default='skato')
    p.add_argument('--label-genes', type=int, default=18)
    p.add_argument('--out-png', required=True)
    return p.parse_args()


def calib_panel(ax, qc, order, colors):
    """Gene-level lambda per scan. One marker per (stratum, method), per cohort."""
    qc = qc.copy()
    # 'low_moderate_high / skato' is 26 characters on a rotated tick and lands on
    # the panel below; the initials are unambiguous and the caption spells them out.
    abbr = {'low_moderate_high': 'L+M+H', 'moderate_high': 'M+H', 'high': 'H'}
    qc['scan'] = qc['stratum'].map(lambda s: abbr.get(s, s)) + '\n' + qc['method']
    scans = sorted(qc['scan'].unique())
    xpos = {s: i for i, s in enumerate(scans)}
    for c in order:
        sel = qc[qc.cohort == c]
        if not len(sel):
            continue
        ax.scatter([xpos[s] for s in sel['scan']], sel['lambda_gc'],
                   s=34, color=colors[c], zorder=3, label=c.replace('_mainland', ''))
    ax.axhline(1.0, color=S.ACCENT, lw=1.0, zorder=2)
    ax.set_xticks(range(len(scans)))
    ax.set_xticklabels(scans, rotation=0, ha='center',
                       fontsize=plt.rcParams['xtick.labelsize'] - 1.5)
    ax.set_ylabel(f'gene-level {S.LAMBDA_GC}')
    S.despine(ax, grid_axis='y')


def forest_rows(g, order, cap):
    """Genes to draw, strongest first, capped. Returns (rows, n_total, n_shown)."""
    if not len(g):
        return [], 0, 0
    best = (g[g.called.isin(['significant', 'suggestive'])]
            .groupby('Gene')['Pvalue'].min().sort_values())
    genes = list(best.index)
    shown = genes[:max(0, cap)]
    rows = []
    for gene in shown:
        for c in order:
            sel = g[(g.Gene == gene) & (g.cohort == c)]
            rows.append((gene, c, sel.iloc[0] if len(sel) else None))
    return rows, len(genes), len(shown)


def forest_panel(ax, rows, order, colors):
    """-log10 P per cohort for each gene. A dot plot, not an effect forest: the
    gene-level test reports a P, and the burden odds ratio lives in the per-gene
    figures where its allele counts are visible beside it."""
    if not rows:
        ax.text(0.5, 0.5, 'no gene reached a tier', transform=ax.transAxes,
                ha='center', va='center', color=S.REFERENCE)
        ax.set_xticks([]); ax.set_yticks([])
        return
    labels, y = [], []
    for i, (gene, cohort, r) in enumerate(rows):
        yy = len(rows) - 1 - i
        y.append(yy)
        labels.append(f'{gene} · {cohort.replace("_mainland", "")}')
        if r is None or str(r['called']) == 'not_tested':
            st = CALLED['not_tested']
            ax.scatter([0.0], [yy], marker=st['marker'], s=st['size'],
                       color=S.NEUTRAL_D, zorder=3)
            continue
        st = CALLED.get(str(r['called']), CALLED['not_a_hit'])
        p = float(r['Pvalue']) if np.isfinite(r['Pvalue']) else np.nan
        if not np.isfinite(p) or p <= 0:
            continue
        ax.scatter([-np.log10(p)], [yy], marker=st['marker'], s=st['size'],
                   facecolors=colors[cohort] if st['filled'] else 'white',
                   edgecolors=colors[cohort], linewidths=1.3, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=plt.rcParams['ytick.labelsize'] - 1.5)
    ax.set_ylim(-0.8, len(rows) - 0.2)
    ax.set_xlabel(S.NEGLOG10P)
    S.despine(ax, grid_axis='x')


def robust_panel(ax, g, order, colors):
    """In how many cohorts each tiered gene is a hit — robustness, not replication."""
    hits = g[g.called.isin(['significant', 'suggestive'])]
    if not len(hits):
        ax.text(0.5, 0.5, 'no gene reached a tier', transform=ax.transAxes,
                ha='center', va='center', color=S.REFERENCE)
        ax.set_xticks([]); ax.set_yticks([])
        return {}
    n = hits.groupby('Gene')['cohort'].nunique().value_counts().sort_index()
    ks = list(range(1, len(order) + 1))
    vals = [int(n.get(k, 0)) for k in ks]
    ax.bar(ks, vals, 0.6, color=[S.COHORT_RAMP[min(k - 1, len(S.COHORT_RAMP) - 1)] for k in ks],
           zorder=3)
    for k, v in zip(ks, vals):
        if v:
            ax.annotate(str(v), xy=(k, v), xytext=(0, 3), textcoords='offset points',
                        ha='center', va='bottom',
                        fontsize=plt.rcParams['legend.fontsize'] - 1)
    ax.set_xticks(ks)
    ax.set_xticklabels([f'{k} of {len(order)}' for k in ks])
    ax.set_xlabel('cohorts in which it is a hit')
    ax.set_ylabel('genes')
    ax.set_ylim(0, max(vals + [1]) * 1.25)
    S.despine(ax, grid_axis='y')
    return dict(zip(ks, vals))


def main():
    args = parse_args()
    S.setup_style()
    order = [c.strip() for c in args.cohort_order.split(',') if c.strip()]
    colors = {c: S.COHORT_RAMP[i % len(S.COHORT_RAMP)] for i, c in enumerate(order)}
    g_all = pd.read_csv(args.genes, sep='\t')
    qc = pd.read_csv(args.scans, sep='\t')
    for c in ('Pvalue', 'FDR_BH'):
        if c in g_all:
            g_all[c] = pd.to_numeric(g_all[c], errors='coerce')

    g = g_all[(g_all.stratum == args.peak_stratum) & (g_all.method == args.peak_method)]
    rows, n_total, n_shown = forest_rows(g, order, args.label_genes)

    h_a, h_c = 1.85, 1.30
    h_b = max(0.9, len(rows) * ROW_IN) + 0.45
    plot_h = h_a + h_b + h_c + 1.30
    fig = plt.figure(figsize=(S.COL_DOUBLE, plot_h + 2.6))
    gs = fig.add_gridspec(3, 1, height_ratios=[h_a, h_b, h_c])
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[i]) for i in range(3))

    calib_panel(ax_a, qc, order, colors)
    S.panel_tag(ax_a, 'a')
    ax_a.legend(frameon=False, ncol=len(order), fontsize=plt.rcParams['legend.fontsize'] - 1,
                loc='upper left', bbox_to_anchor=(0.0, 1.02))
    forest_panel(ax_b, rows, order, colors)
    S.panel_tag(ax_b, 'b')
    counts = robust_panel(ax_c, g, order, colors)
    S.panel_tag(ax_c, 'c')

    lam_lo, lam_hi = qc['lambda_gc'].min(), qc['lambda_gc'].max()
    n_sig = g[g.called == 'significant']['Gene'].nunique()
    stratum_h = args.peak_stratum.replace('_', '+').upper()
    S.caption_block(
        fig, plot_h=plot_h,
        title=(f'Gene-level {S.LAMBDA_GC} runs {lam_lo:.2f}–{lam_hi:.2f} across the scans; '
               f'{n_total} gene(s) reach a tier in at least one cohort under {stratum_h} / '
               f'{args.peak_method.upper()}, {n_sig} of them at FDR < 0.05.'),
        panels=[
            f'Gene-level {S.LAMBDA_GC} for every impact stratum and test method, per cohort; '
            'the red line is 1. H, HIGH; M+H, MODERATE+HIGH; L+M+H, LOW+MODERATE+HIGH.',
            f'Every gene reaching a tier in any cohort, shown in all three — the '
            f'{n_shown} strongest of {n_total}. Marker shape gives that cohort’s own call.',
            'In how many cohorts each tiered gene is a hit.',
        ],
        notes=('The cohorts are NESTED and share their cases entirely, so (c) measures '
               'robustness to the ancestry filter, NOT replication: a chance hit in the '
               'shared core reappears in all three by construction. Panels (b) and (c) show '
               f'{stratum_h} / {args.peak_method.upper()}; every stratum and method is in '
               '`gene_crosscohort.tsv`. Full explanation: cohort_compare_genes.md'),
        top_pad=0.52, hspace=0.78, left='auto', right=0.985,
        margin_axes=[ax_a, ax_b, ax_c])
    fig.savefig(args.out_png)
    plt.close(fig)

    figure_doc.write_doc(
        args.out_png,
        title='Gene-based results across the three cohorts',
        question=('Does a gene-based hit survive relaxing the ancestry filter, and are the '
                  'scans themselves calibrated?'),
        panels=[
            ('a', 'Scan calibration',
             'Gene-level genomic-control inflation factor for each impact stratum and test '
             'method, per cohort. Computed over the gene P-values, not variant ones, so it '
             'measures whether the GENE-level null is calibrated.'),
            ('b', 'Every tiered gene, in every cohort',
             'The union of genes reaching a tier in any cohort, each shown in all three. '
             'A dot plot of -log10 P rather than an effect forest: the gene-level test '
             'reports a P, and the burden odds ratio belongs beside its allele counts, '
             'which is where the per-gene figures put it.'),
            ('c', 'Robustness to the ancestry filter',
             'How many of the three cohorts call each tiered gene a hit.'),
        ],
        reading=[
            'Read (c) as a sensitivity check, never as replication. The cohorts are nested '
            'and share every case, so a gene appearing in all three has not been confirmed '
            'by independent data — it has been shown not to depend on the ancestry cut.',
            'A gene called in the widest cohort ONLY is the one to distrust: what full adds '
            'over narrow is exactly the ancestry-outlier samples.',
            'Check (a) before (b). A stratum whose lambda is far from 1 has a mis-calibrated '
            'gene-level null, and its tiers are correspondingly optimistic.',
        ],
        limits=[
            'It cannot establish that any gene is real. No independent cohort exists for this '
            'phenotype, so nothing here is external validation.',
            'Tiers come from Benjamini-Hochberg within each scan, so a tier in the HIGH '
            'stratum (few genes tested) and one in MODERATE+HIGH are not the same evidence.',
        ],
        defs=['lambda_gc'],
        numbers=[('cohorts', ', '.join(order)),
                 ('stratum shown', args.peak_stratum), ('method shown', args.peak_method),
                 ('genes reaching a tier', n_total), ('genes shown', n_shown),
                 ('significant somewhere', n_sig),
                 ('gene-level lambda, min', round(float(lam_lo), 4)),
                 ('gene-level lambda, max', round(float(lam_hi), 4))]
        + [(f'genes called in {k} of {len(order)} cohorts', v) for k, v in counts.items()],
        interpretation=('Under a nested design the useful signal in a cross-cohort comparison '
                        'is instability, not agreement. A gene that appears only when the '
                        'ancestry filter is relaxed is a candidate artefact; one that holds '
                        'across all three has survived the only internal check available, '
                        'which is weaker than replication and should be described as such.'),
        methods_ref='../../docs/METHODS.md')
    print(f'[plot_cohort_compare_genes] {n_total} tiered gene(s), {n_shown} shown -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in plot_cohort_compare_genes: {e}', file=sys.stderr)
        sys.exit(1)
