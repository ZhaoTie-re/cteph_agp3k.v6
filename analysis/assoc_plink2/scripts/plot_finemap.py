#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : SuSiE fine-mapping figure for one locus.
#             (a) -log10 P against position, coloured by in-sample LD with the lead
#             (b) posterior inclusion probability, STACKED under (a) on the same
#                 genomic axis, credible-set members ringed
#             (c) the PIP-ordered credible set as a cumulative-mass curve, which is
#                 what says whether the signal resolves to a few variants or not
#           (a) over (b) rather than side by side: posterior mass is then read
#           directly under the association peak it belongs to, instead of the
#           reader re-locating each position by eye across two panels.
# Component: assoc_plink2
# Used by : assoc_plink2.nf  process PLOT_FINEMAP
# ---------------------------------------------------------------------------
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_style as S
import figure_doc

PLOT_H = 3.6          # inches of croppable plot block (titles + axes + x-labels)
# 3.1, not 2.4: three panels across 7.2 in were cramped enough that the CS
# bracket label and the PIP axis competed for the same rows, and at 2.4 in the
# figure was 42 % caption however short the text got.

# Credible sets are an unordered discrete variable -> the shared series.
CS_COLORS = S.SERIES


def parse_args():
    p = argparse.ArgumentParser(description='SuSiE fine-mapping figure for one locus.')
    p.add_argument('--pip', required=True)
    p.add_argument('--cs', required=True)
    p.add_argument('--json')
    p.add_argument('--ld-cohort', help='<locus>.ld_cohort.tsv for point colouring.')
    p.add_argument('--cohort', required=True)
    p.add_argument('--locus-id', required=True)
    p.add_argument('--lead-id', required=True)
    p.add_argument('--out-png', required=True)
    return p.parse_args()


def main():
    args = parse_args()
    S.setup_style('paper')
    pip = pd.read_csv(args.pip, sep='\t', dtype={'ID': str, 'CHROM': str})
    cs = pd.read_csv(args.cs, sep='\t') if Path(args.cs).stat().st_size else pd.DataFrame()
    meta = json.loads(Path(args.json).read_text()) if args.json and Path(args.json).exists() else {}

    ld = None
    if args.ld_cohort and Path(args.ld_cohort).exists():
        ld = pd.read_csv(args.ld_cohort, sep='\t', dtype={'ID': str}).set_index('ID')

    x = pip['POS'].to_numpy() / 1e6
    if ld is not None:
        al = ld.reindex(pip['ID'])
        r2 = al['r2'].to_numpy(dtype=float)
        known = al['state'].fillna('not_in_panel').to_numpy() != 'not_in_panel'
    else:
        r2, known = np.full(len(pip), np.nan), np.zeros(len(pip), dtype=bool)
    colors = S.ld_bin_colors(r2, known=known)
    order = np.argsort(np.nan_to_num(np.where(known, r2, -1.0)))

    # LEFT COLUMN: association over PIP, sharing the genomic axis, so posterior
    # mass is read directly under the association peak it belongs to. Side by side
    # the reader had to re-locate every position by eye. RIGHT: the cumulative
    # curve, spanning both rows.
    fig = plt.figure(figsize=(S.COL_DOUBLE, PLOT_H))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.4, 1.0], height_ratios=[1.0, 1.0])
    ax_p = fig.add_subplot(gs[0, 0])
    ax_pip = fig.add_subplot(gs[1, 0], sharex=ax_p)
    ax_cum = fig.add_subplot(gs[:, 1])

    # (a) association
    ax = ax_p
    y = -np.log10(pip['P'].to_numpy())
    ax.scatter(x[order], y[order], s=8, c=colors[order], linewidths=0.3,
               edgecolors='white', rasterized=True, zorder=3)
    lead = (pip['ID'] == args.lead_id).to_numpy()
    if lead.any():
        ax.scatter(x[lead], y[lead], s=34, marker='D', c=S.LEAD_COLOR,
                   edgecolors='white', linewidths=1.0, zorder=5)
    ax.axhline(-np.log10(S.GW_ALPHA), color=S.ACCENT, lw=1.1, zorder=2)
    ax.set_ylim(0, max(-np.log10(pip['P'].min()) * 1.16, -np.log10(S.GW_ALPHA) * 1.12))
    ax.set_ylabel(S.NEGLOG10P)
    plt.setp(ax.get_xticklabels(), visible=False)      # (b) carries the axis
    S.panel_tag(ax, 'a')
    S.despine(ax, grid_axis='y')

    # (b) PIP, credible-set members ringed in the set's colour
    ax = ax_pip
    pv = pip['pip'].to_numpy()
    ax.scatter(x[order], pv[order], s=8, c=colors[order],
               linewidths=0.3, edgecolors='white', rasterized=True, zorder=3)
    # Fitted to the DATA, not pinned at 1.16. At a top PIP of 0.43 a fixed ceiling
    # left three quarters of the panel empty and crushed every point into the
    # bottom strip; the label band above still needs room, hence 1.42.
    top = float(np.nanmax(pv)) if len(pv) else 1.0
    ax.set_ylim(-0.02 * top, max(top * 1.42, 0.12))
    span = (x.max() - x.min()) if len(x) else 1.0
    for i, (_, r) in (enumerate(cs.iterrows()) if len(cs) else []):
        members = (pip['cs'] == r['cs']).to_numpy()
        col = CS_COLORS[i % len(CS_COLORS)]
        ax.scatter(x[members], pv[members], s=46, facecolors='none', edgecolors=col,
                   linewidths=1.6, zorder=6)
        lo, hi = x[members].min(), x[members].max()
        # A credible set can span a few kb of a 500 kb window, at which point the
        # |--| bracket collapses to a single glyph and reads as a stray mark. Below
        # 2 % of the window, mark the cluster with one vertical line instead.
        if (hi - lo) < 0.02 * span:
            # A short segment over the set's OWN PIP range, not an axvline across
            # the whole panel — the marker should point at the cluster, not divide
            # the figure.
            ax.plot([(lo + hi) / 2] * 2, [0, pv[members].max() * 1.12],
                    color=col, lw=1.0, ls=':', zorder=4)
        else:
            ax.annotate('', xy=(lo, 0.88), xytext=(hi, 0.88),
                        xycoords=('data', 'axes fraction'),
                        arrowprops=dict(arrowstyle='|-|,widthA=0.35,widthB=0.35',
                                        color=col, lw=1.4))
        ax.text((lo + hi) / 2, 0.91, f"CS{int(r['cs'])} · {int(r['size'])} variant"
                                     f"{'' if r['size'] == 1 else 's'}",
                transform=ax.get_xaxis_transform(), ha='center', va='bottom',
                color=col, fontweight='bold', fontsize=plt.rcParams['legend.fontsize'] - 1)
    ax.set_ylabel('PIP')
    ax.set_xlabel(f'chr{pip.CHROM.iloc[0]} {S.POS_MB}'.replace('chrchr', 'chr'))
    S.panel_tag(ax, 'b')
    S.despine(ax, grid_axis='y')

    # (c) how fast posterior mass accumulates
    ax = ax_cum
    if len(cs):
        for i, (_, r) in enumerate(cs.iterrows()):
            members = pip[pip['cs'] == r['cs']].sort_values('pip', ascending=False)
            cum = np.cumsum(members['pip'].to_numpy())
            col = CS_COLORS[i % len(CS_COLORS)]
            ax.step(np.arange(1, len(cum) + 1), cum, where='post', color=col, lw=2.0)
            ax.scatter([1], [cum[0]], s=18, color=col, zorder=5)
            # Name the curve in place, anchored at its FIRST step. A legend for a
            # single credible set is a key to one thing, and place_legend put it in
            # the panel-title row; anchoring at the last step instead put it on the
            # 0.95 coverage line, which every curve ends beside.
            ax.annotate(f"CS{int(r['cs'])}", xy=(1, cum[0]), xytext=(8, 9),
                        textcoords='offset points', ha='left', va='bottom', color=col,
                        fontweight='bold', fontsize=plt.rcParams['legend.fontsize'])
        ax.axhline(0.95, color=S.REFERENCE, ls=':', lw=1.1)
        ax.text(0.02, 0.95, 'coverage 0.95', transform=ax.get_yaxis_transform(),
                ha='left', va='bottom', color=S.REFERENCE,
                fontsize=plt.rcParams['legend.fontsize'] - 1)
        ax.set_xlim(0.5, max(6, max(cs['size']) + 0.5))
        ax.set_ylim(0, 1.08)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    else:
        ax.text(0.5, 0.5, 'no credible set at\n95 % coverage', transform=ax.transAxes,
                ha='center', va='center', color=S.REFERENCE)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    ax.set_xlabel('variants included, PIP-ranked')
    ax.set_ylabel('cumulative posterior mass')
    S.panel_tag(ax, 'c')
    S.despine(ax)

    # Spans the whole figure width, not just panel (a); at 7 entries it will not
    # fit above a 2.2 in panel.
    fig.legend(handles=S.ld_legend_handles(include_unknown=True, lead=True),
               loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=7, frameon=False,
               title=f'in-sample {S.R2_LD} with the lead', handletextpad=0.4,
               columnspacing=1.1, fontsize=plt.rcParams['legend.fontsize'])

    n_cs = len(cs)
    top = cs.sort_values('top_pip', ascending=False).iloc[0] if n_cs else None
    S.caption_block(
        fig,
        title=(f'{args.cohort}, locus {args.locus_id}: SuSiE resolves the signal to '
               f'{n_cs} credible set'
               f"{'' if n_cs == 1 else 's'}"
               + (f", the largest-PIP set holding {int(top['size'])} variant"
                  f"{'' if top['size'] == 1 else 's'} (top PIP {top['top_pip']:.3f})." if n_cs else '.')),
        panels=[
            'Association over the locus, coloured by in-sample $r^{2}$ with the lead (purple diamond).',
            'Per-variant PIP on the same axis; rings mark credible-set membership.',
            'Cumulative posterior mass over the PIP-ranked variants of each set.',
        ],
        notes=(f"susie_rss, {meta.get('n_variants', len(pip)):,} variants, in-sample LD; "
               f"$L$ = {meta.get('L', 10)}, coverage {meta.get('coverage', 0.95)}, "
               f"converged = {meta.get('converged', 'NA')}. "
               f"Full explanation: {Path(args.out_png).stem}.md"),
        plot_h=PLOT_H, top_pad=0.62, wspace=0.40, hspace=0.22, left='auto', right=0.985,
        margin_axes=[ax_p, ax_pip, ax_cum])
    fig.savefig(args.out_png)
    plt.close(fig)

    figure_doc.write_doc(
        args.out_png,
        title=f'Fine-mapping — {args.cohort}, peak {args.locus_id}',
        question='Does the posterior concentrate on a few variants, or does LD spread it out?',
        panels=[
            ('a', 'Association',
             'The peak window, coloured by in-sample r2 with the lead (purple diamond). Shown so the '
             'PIP panel can be read against the evidence that produced it.'),
            ('b', 'Posterior inclusion probability',
             'Per-variant PIP from susie_rss on the summary statistics and the in-sample signed-r LD '
             'matrix, L = 10. Rings mark credible-set membership; the bracket spans a set\'s physical '
             'extent. A credible set is the smallest group of variants carrying 95% of the posterior '
             'mass for one signal.'),
            ('c', 'Resolution',
             'Cumulative posterior mass against the PIP-ranked members of each set. A curve reaching '
             '0.95 in a few steps means the signal is resolved to those variants; a slow curve means '
             'LD has spread the mass and the set cannot be narrowed at this sample size.'),
        ],
        numbers=[('cohort', args.cohort), ('peak', args.locus_id), ('lead variant', args.lead_id),
                 ('variants fine-mapped', int(meta.get('n_variants', len(pip)))),
                 ('GWAS n used', int(meta.get('n_gwas', 0))),
                 ('L (max signals)', meta.get('L')),
                 ('coverage', meta.get('coverage')),
                 ('purity filter min |r|', meta.get('min_abs_corr')),
                 ('converged', meta.get('converged')),
                 ('credible sets', n_cs),
                 ('largest-PIP set size', int(top['size']) if n_cs else None),
                 ('top PIP', float(top['top_pip']) if n_cs else None),
                 ('susieR version', meta.get('susieR_version'))],
        reading=[
            'A credible set of one or two variants is a resolved signal. A set of twenty means the '
            'data cannot distinguish among them, not that twenty variants are causal.',
            'Check that the set contains the lead. If it does not, the lead is a tag and the '
            'posterior prefers a neighbour.',
            'PIP is conditional on the LD matrix being the one the statistics came from, which is why '
            'in-sample LD is used here rather than a reference panel.',
        ],
        limits=[
            'SuSiE assumes the causal variant is present in the data. A causal variant not genotyped '
            'or filtered out cannot appear, and its posterior mass will be distributed over its tags.',
            'It cannot rank the biological plausibility of set members — only their statistical '
            'compatibility with the observed association pattern.',
        ],
        defs=['pip', 'ld_r2'],
        interpretation=('Fine-mapping is conditional on the LD matrix being the one the statistics came from, '
                   'which is why in-sample LD is used here. At this effective sample size the posterior is '
                   'driven by a small number of strongly associated variants, so a wide credible set should '
                   'be read as insufficient resolution rather than as evidence against a single causal variant.'),
        model=S.FORMULAS['susie'],
        methods_ref='../../../docs/METHODS.md')
    print(f'[plot_finemap] {args.cohort} {args.locus_id}: {n_cs} credible set(s) -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in plot_finemap: {e}', file=sys.stderr)
        sys.exit(1)
