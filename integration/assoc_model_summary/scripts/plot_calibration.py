#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Calibration of two scan models side by side, per cohort:
#             a  lambda_GC per cohort
#             b  lambda_GC by minor-allele-frequency bin (is the inflation
#                genome-wide, or confined to some frequencies?)
#             c  variants at P < 1e-5 against the number expected by chance
#           Slide-ready: panel letters and a model legend, no caption block.
# Inputs  : tail_counts.tsv, lambda_by_maf.tsv (scan_tables.py)
# Outputs : PNG
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
import common as C  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cohorts', required=True)
    ap.add_argument('--cohort-labels', required=True)
    ap.add_argument('--cohort-colors', required=True)
    ap.add_argument('--models', required=True)
    ap.add_argument('--tail', required=True)
    ap.add_argument('--maf', required=True)
    ap.add_argument('--plot-style', required=True)
    ap.add_argument('--out-png', required=True)
    a = ap.parse_args()
    S = C.import_style(a.plot_style)
    psize = C.setup(S)
    cohorts, models = C.csv_list(a.cohorts), C.csv_list(a.models)
    coh = C.cohort_table(a.cohorts, a.cohort_labels, a.cohort_colors)
    tail = pd.read_csv(a.tail, sep='\t').set_index(['model', 'cohort'])
    maf = pd.read_csv(a.maf, sep='\t')
    xc = np.arange(len(cohorts))
    xlab = [coh[c]['label'].replace(' ', '\n', 1) for c in cohorts]

    W, H = 10.0, 3.8
    fig = plt.figure(figsize=(W, H))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.35, 1.1])
    fig.subplots_adjust(left=0.62 / W, right=1 - 0.12 / W, top=1 - 1.0 / H, bottom=0.62 / H, wspace=0.36)
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0, k]) for k in range(3))

    # a: lambda per cohort
    lams = []
    for k, m in enumerate(models):
        y = [float(tail.loc[(m, c), 'lambda_gc']) for c in cohorts]
        lams += y
        col = C.MODELS[m]['color']
        ax_a.plot(xc, y, color=col, lw=1.2, zorder=3)
        ax_a.scatter(xc, y, s=42, color=col, edgecolors='white', linewidths=0.8, zorder=4)
        for x, v in zip(xc, y):
            up = v >= max(float(tail.loc[(mm, cohorts[int(x)]), 'lambda_gc']) for mm in models)
            ax_a.annotate(f'{v:.3f}', (x, v), xytext=(0, 6 if up else -7), textcoords='offset points',
                          ha='center', va='bottom' if up else 'top', fontsize=plt.rcParams['legend.fontsize'],
                          color=col, fontweight='bold')
    ax_a.axhline(1.0, color=S.ACCENT, lw=1.0, zorder=2)
    lo, hi = min(lams + [1.0]), max(lams)
    ax_a.set_ylim(lo - 0.025, hi + 0.03)
    ax_a.set_xlim(-0.45, len(cohorts) - 0.55)
    ax_a.set_xticks(xc)
    ax_a.set_xticklabels(xlab)
    ax_a.set_ylabel(S.LAMBDA_GC)
    S.despine(ax_a, grid_axis='y')
    S.panel_tag(ax_a, 'a', 'Genomic inflation', size=psize)

    # b: lambda by MAF bin
    bins = [b for b in dict.fromkeys(maf.maf_bin) if b != '<0.01']     # <0.01: reported, not plotted
    xb = np.arange(len(bins))
    for m in models:
        for c in cohorts:
            d = maf[(maf.model == m) & (maf.cohort == c)].set_index('maf_bin').loc[bins]
            ax_b.plot(xb, d.lambda_gc, color=C.MODELS[m]['color'], lw=1.0, marker=coh[c]['marker'], ms=4.2,
                      mec='white', mew=0.6, zorder=3)
    ax_b.axhline(1.0, color=S.ACCENT, lw=1.0, zorder=2)
    ax_b.set_xticks(xb)
    ax_b.set_xticklabels(bins)
    ax_b.set_xlabel('minor-allele frequency')
    ax_b.set_ylabel(S.LAMBDA_GC)
    S.despine(ax_b, grid_axis='y')
    S.panel_tag(ax_b, 'b', f'{S.LAMBDA_GC} by allele frequency', size=psize)


    # c: tail counts against expectation
    wbar = 0.8 / len(models)
    exp = []
    for k, m in enumerate(models):
        y = [int(tail.loc[(m, c), 'n_p_lt_1e5']) for c in cohorts]
        exp += [float(tail.loc[(m, c), 'expected_p_lt_1e5']) for c in cohorts]
        x = xc + (k - (len(models) - 1) / 2) * wbar
        ax_c.bar(x, y, width=wbar * 0.92, color=C.MODELS[m]['color'], zorder=3)
        for xx, v in zip(x, y):
            ax_c.text(xx, v, f'{v}', ha='center', va='bottom', fontsize=plt.rcParams['legend.fontsize'],
                      color=S.INK, zorder=4)
    e = float(np.mean(exp))
    ax_c.axhline(e, color=S.INK_SOFT, lw=1.0, ls='--', zorder=4)
    ax_c.legend(handles=[Line2D([], [], color=S.INK_SOFT, lw=1.0, ls='--',
                                label=f'expected by chance ≈ {e:.0f}')],
                loc='upper left', frameon=False, fontsize=plt.rcParams['legend.fontsize'], borderaxespad=0.2)
    ax_c.set_xticks(xc)
    ax_c.set_xticklabels(xlab)
    ax_c.set_xlim(-0.5, len(cohorts) - 0.5)
    ax_c.set_ylim(0, max(tail.n_p_lt_1e5) * 1.15)
    ax_c.set_ylabel(r'variants with $P < 10^{-5}$')
    S.despine(ax_c, grid_axis='y')
    S.panel_tag(ax_c, 'c', r'Excess of $P < 10^{-5}$', size=psize)

    handles = [Line2D([], [], color=C.MODELS[m]['color'], marker='o', ms=6, lw=1.2,
                      label=f'{C.MODELS[m]["name"]} · {C.MODELS[m]["engine"]} · {C.MODELS[m]["covar"]}')
               for m in models]
    fig.legend(handles=handles, loc='upper center', ncol=len(models), frameon=False,
               bbox_to_anchor=(0.5, 1 - 0.04 / H), fontsize=plt.rcParams['legend.fontsize'] + 0.5,
               handlelength=1.8, columnspacing=2.2)
    fig.legend(handles=[Line2D([], [], color=S.INK_SOFT, marker=coh[c]['marker'], ls='', ms=5,
                               label=coh[c]['label']) for c in cohorts],
               loc='upper center', ncol=len(cohorts), frameon=False, bbox_to_anchor=(0.5, 1 - 0.30 / H),
               fontsize=plt.rcParams['legend.fontsize'], handletextpad=0.3, columnspacing=1.8,
               title='cohort (panel b)', title_fontsize=plt.rcParams['legend.fontsize'])
    fig.savefig(a.out_png)
    plt.close(fig)
    C.log(f'calibration {models} -> {a.out_png}')


if __name__ == '__main__':
    main()
