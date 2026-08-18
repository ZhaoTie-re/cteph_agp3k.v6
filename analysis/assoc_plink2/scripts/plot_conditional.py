#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Conditional-analysis figure for one locus. The question the figure
#           answers is exactly one: does this locus carry one signal or more?
#             (a) -log10 P over the locus before and after each conditioning step
#             (b) the top statistic per round, which is the decision the stepwise
#                 procedure actually made, with the threshold it was compared to
#           A signal that collapses to the null once the lead is conditioned on
#           was tagging the lead; one that survives is independent.
# Component: assoc_plink2
# Used by : assoc_plink2.nf  process PLOT_CONDITIONAL
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
import plot_style as S
import figure_doc

PLOT_H = 2.9          # inches of croppable plot block (titles + axes + x-labels)
# 2.9, not 2.5: panel (b) carries two-line tick labels under each bar, which a
# 2.5 in block left no room for.

# Conditioning rounds are an ORDERED sequence, so round 0 takes the DATA colour and
# later rounds step through the shared discrete series.
ROUND_COLORS = S.SERIES + [S.REFERENCE]


def parse_args():
    p = argparse.ArgumentParser(description='Conditional analysis figure for one locus.')
    p.add_argument('--cond-dir', required=True)
    p.add_argument('--rounds', required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--locus-id', required=True)
    p.add_argument('--lead-id', required=True)
    p.add_argument('--alpha', type=float, default=S.GW_ALPHA)
    p.add_argument('--pc-label', default='',
                   help='Space the PCs were computed in, named in the model formula.')
    p.add_argument('--covar-label', default='SEX + PCs',
                   help='Human-readable covariate set, printed in the caption and the sidecar. '
                        'Configuration, not a property of the method — pass the run\'s own.')
    p.add_argument('--n-pcs', type=int, default=0,
                   help='Number of PCs in the covariate set; 0 leaves the formula\'s upper limit '
                        'as a generic K.')
    p.add_argument('--out-png', required=True)
    return p.parse_args()


def main():
    args = parse_args()
    S.setup_style('paper')
    rounds = pd.read_csv(args.rounds, sep='\t')
    files = sorted(Path(args.cond_dir).glob(f'{args.locus_id}.round*.tsv'))
    files = [f for f in files if not f.name.endswith('.rounds.tsv')]
    per_round = {}
    for f in files:
        rnd = int(f.name.rsplit('round', 1)[1].split('.')[0])
        per_round[rnd] = pd.read_csv(f, sep='\t', dtype={'CHROM': str, 'ID': str})

    fig, axes = plt.subplots(1, 2, figsize=(S.COL_DOUBLE, PLOT_H),
                             gridspec_kw={'width_ratios': [1.6, 1.0]})

    # (a) the locus, one trace per conditioning round
    ax = axes[0]
    thr = -np.log10(args.alpha)
    for rnd in sorted(per_round):
        d = per_round[rnd]
        d = d[d['P'].notna() & (d['P'] > 0)]
        col = ROUND_COLORS[rnd % len(ROUND_COLORS)]
        lbl = 'unconditioned' if rnd == 0 else f'conditioned on {rnd} variant{"" if rnd == 1 else "s"}'
        ax.scatter(d['POS'] / 1e6, -np.log10(d['P']), s=6 if rnd else 8,
                   c=col, alpha=0.95 if rnd == 0 else 0.75, linewidths=0,
                   rasterized=True, zorder=3 + rnd, label=lbl)
    if 0 in per_round:
        d0 = per_round[0]
        lead = d0[d0['ID'] == args.lead_id]
        if len(lead):
            ax.scatter(lead['POS'] / 1e6, -np.log10(lead['P']), s=42, marker='D',
                       facecolors='none', edgecolors=S.LEAD_COLOR, linewidths=1.8, zorder=9)
            # Inside the panel, not above the point: an offset label on the
            # tallest point lands squarely on the panel title.
            ax.annotate(args.lead_id, xy=(lead['POS'].iloc[0] / 1e6, -np.log10(lead['P'].iloc[0])),
                        xytext=(0.02, 0.95), textcoords='axes fraction', ha='left', va='top',
                        fontsize=plt.rcParams['legend.fontsize'] - 1, fontweight='bold',
                        color=S.LEAD_COLOR,
                        arrowprops=dict(arrowstyle='-', color=S.LEAD_COLOR, lw=0.6,
                                        shrinkA=2, shrinkB=4))
    ax.axhline(thr, color=S.ACCENT, lw=1.1, zorder=2)
    ax.text(0.995, thr, r'$5\times10^{-8}$', transform=ax.get_yaxis_transform(),
            ha='right', va='bottom', color=S.ACCENT,
            fontsize=plt.rcParams['legend.fontsize'] - 1)
    # Headroom so the lead label clears the threshold line it sits beside.
    ymax = max([-np.log10(d[d['P'] > 0]['P']).max() for d in per_round.values()] + [thr])
    ax.set_ylim(0, ymax * 1.22)
    chrom = str(per_round[0]['CHROM'].iloc[0]).replace('chr', '') if per_round else ''
    ax.set_xlabel(f'chr{chrom} {S.POS_MB}')
    ax.set_ylabel(S.NEGLOG10P)
    S.panel_tag(ax, 'a')
    S.despine(ax, grid_axis='y')

    # (b) the decision the procedure made
    ax = axes[1]
    r = rounds[rounds['top_p'].notna()].copy()
    if len(r):
        yv = -np.log10(r['top_p'].astype(float))
        cols = [ROUND_COLORS[int(i) % len(ROUND_COLORS)] for i in r['round']]
        # 0.42, not 0.6: with two rounds a 0.6-wide bar is ~1.5 in across and the
        # summary panel reads heavier than the locus panel it summarises.
        ax.bar(r['round'], yv, width=0.42, color=cols, edgecolor='white',
               linewidth=0.8, zorder=3)
        for _, rr in r.iterrows():
            ax.annotate(f"$P$ = {S.p_tex(rr['top_p'])}",
                        xy=(rr['round'], -np.log10(float(rr['top_p']))),
                        xytext=(0, 5), textcoords='offset points',
                        ha='center', va='bottom',
                        fontsize=plt.rcParams['legend.fontsize'] - 1.5,
                        color=S.INK, zorder=6)
        ax.set_ylim(0, max(yv.max() * 1.30, thr * 1.25))
        # The variant identified in each round goes on the SECOND LINE of that
        # round's tick label. Above the bars a 19-character ID printed across its
        # neighbour; rotated inside the bar it was clipped whenever the bar was
        # shorter than the text (round 1 here, at -log10 P = 3.5). Under the axis
        # it has the full column width and can collide with nothing.
        ids = dict(zip(r['round'].astype(int), r['top_id'].astype(str)))
        order = sorted(r['round'].astype(int))
        ax.set_xticks(order)
        ax.set_xticklabels([('unconditioned' if i == 0 else f'+{i}') + f'\n{ids.get(i, "")}'
                            for i in order],
                           fontsize=plt.rcParams['xtick.labelsize'] - 1.5)
    ax.axhline(thr, color=S.ACCENT, lw=1.1, zorder=4)
    ax.set_xlabel('conditioning round (variants held fixed)')
    ax.set_ylabel(f'top-variant {S.NEGLOG10P}')
    S.panel_tag(ax, 'b')
    S.despine(ax, grid_axis='y')

    handles = [Line2D([0], [0], marker='o', ls='none', color=ROUND_COLORS[i % len(ROUND_COLORS)],
                      markersize=7,
                      label='unconditioned' if i == 0 else f'conditioned on {i} variant{"" if i == 1 else "s"}')
               for i in sorted(per_round)]
    handles.append(Line2D([0], [0], marker='D', ls='none', markerfacecolor='none',
                          markeredgecolor=S.LEAD_COLOR, markeredgewidth=1.8, markersize=9,
                          label='lead variant'))
    S.legend_above(axes[0], handles, ncol=min(len(handles), 3), y=1.02)

    n_sig = int(rounds['passes_threshold'].astype(str).str.upper().isin(['TRUE', 'T']).sum())
    S.caption_block(
        fig,
        title=(f'{args.cohort}, locus {args.locus_id}: stepwise conditioning leaves '
               f'{n_sig} independent signal{"" if n_sig == 1 else "s"}.'),
        panels=[
            'Association across the locus at each conditioning round; round 0 is unconditioned.',
            'The statistic the procedure acted on: smallest $P$ remaining per round, against the '
            'same threshold.',
        ],
        notes=(f"{len(rounds)} round(s), max "
               f"{int(rounds['n_conditioned_on'].max()) if len(rounds) else 0} variant(s) conditioned "
               f"on; {args.covar_label} plus the conditioning genotypes. "
               f"Full explanation: {Path(args.out_png).stem}.md"),
        plot_h=PLOT_H, top_pad=0.50, wspace=0.32, left='auto', right=0.985, margin_axes=axes)
    fig.savefig(args.out_png)
    plt.close(fig)

    figure_doc.write_doc(
        args.out_png,
        title=f'Conditional analysis — {args.cohort}, peak {args.locus_id}',
        question='Does this peak carry one association signal, or more than one?',
        panels=[
            ('a', 'The peak, by conditioning round',
             'Round 0 is the unconditioned fit restricted to the peak window. Each later round adds '
             'the previous round\'s top variant to the covariate set and re-fits every variant in the '
             'window. Same samples, same model and the same covariate set as the genome-wide scan, so the '
             'rounds are comparable to it and to each other.'),
            ('b', 'The stepwise decision',
             'The statistic the procedure actually acted on: the smallest P remaining at each round, '
             'against the same genome-wide threshold. The loop stops at the first round in which '
             'nothing in the window clears it.'),
        ],
        numbers=[('cohort', args.cohort), ('peak', args.locus_id), ('lead variant', args.lead_id),
                 ('rounds run', len(rounds)),
                 ('independent signals', n_sig),
                 ('max variants conditioned on',
                  int(rounds['n_conditioned_on'].max()) if len(rounds) else 0)],
        reading=[
            'If the whole peak collapses after conditioning on the lead, every significant variant '
            'there was tagging one underlying signal — the expected result inside a single LD block.',
            'A residual peak that still clears the threshold is a second, independent signal, and the '
            'procedure will have added it to the conditioning set and continued.',
            'Panel (b) is the audit trail: it shows the number the stopping rule compared, round by '
            'round.',
        ],
        limits=[
            'It cannot separate two causal variants in near-perfect LD; conditioning on one removes '
            'both. That limit is a property of the sample, not of the method.',
            'Absence of a secondary signal at this N is weak evidence — the conditional test faces the '
            'same detection floor as the primary scan, applied to a residual effect.',
        ],
        defs=['conditional', 'gw_sig', 'model'],
        interpretation=('Collapse of the whole locus after conditioning on the lead means every significant variant '
                   'there was tagging one underlying signal — the expected outcome inside a single LD block. '
                   'A residual peak would mark a second, independent signal. Conditioning is performed on the '
                   'same samples, model and covariates as the genome-wide scan, so the rounds are comparable '
                   'to it and to each other.'),
        model=S.glm_formula(args.pc_label, args.n_pcs),
        methods_ref='../../../docs/METHODS.md')
    print(f'[plot_conditional] {args.cohort} {args.locus_id}: {n_sig} signal(s) -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in plot_conditional: {e}', file=sys.stderr)
        sys.exit(1)
