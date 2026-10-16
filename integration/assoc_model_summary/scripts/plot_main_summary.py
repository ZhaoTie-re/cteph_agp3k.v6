#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The main analysis at a glance, for one slide:
#             a  per cohort: cases and controls, lambda_GC, genome-wide and
#                suggestive loci
#             b  every locus genome-wide in any cohort or suggestive in all
#                three: OR (95% CI) of its representative variant in each
#                cohort, one row per locus, the cohorts side by side; a diamond
#                where that cohort called the variant genome-wide
# Inputs  : tail_counts.tsv (scan_tables.py, model grm_pc); assoc_saige
#           _comparison/tables/lead_crosscohort.tsv
# Outputs : PNG, loci table (TSV)
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

GROUPS = [('b', 'Genome-wide in at least one cohort'), ('c', 'Suggestive in all three cohorts')]


def tint(hex_color, f):
    """hex_color mixed with white; f = share of the colour."""
    h = hex_color.lstrip('#')
    rgb = np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)]) / 255.0
    return tuple(rgb * f + (1 - f))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cohorts', required=True)
    ap.add_argument('--cohort-labels', required=True)
    ap.add_argument('--cohort-colors', required=True)
    ap.add_argument('--model', default='grm_pc')
    ap.add_argument('--tail', required=True)
    ap.add_argument('--crosscohort', required=True, help='assoc_saige lead_crosscohort.tsv')
    ap.add_argument('--plot-style', required=True)
    ap.add_argument('--out-png', required=True)
    ap.add_argument('--out-loci', required=True)
    a = ap.parse_args()
    S = C.import_style(a.plot_style)
    psize = C.setup(S)
    cohorts = C.csv_list(a.cohorts)
    coh = C.cohort_table(a.cohorts, a.cohort_labels, a.cohort_colors)
    tail = pd.read_csv(a.tail, sep='\t')
    tail = tail[tail.model == a.model].set_index('cohort')

    x = pd.read_csv(a.crosscohort, sep='\t', dtype={'chrom': str, 'variant_id': str})
    x = x[x.panel.isin([g for g, _ in GROUPS]) & (x.locus_rep == 1)].copy()
    for c in ('OR', 'L95', 'U95', 'P'):
        x[c] = pd.to_numeric(x[c], errors='coerce')
    minp = x.groupby('locus_id').P.min()
    x['_minp'] = x.locus_id.map(minp)
    loci = (x.drop_duplicates('locus_id').assign(_g=lambda d: d.panel.map({g: k for k, (g, _) in enumerate(GROUPS)}))
            .sort_values(['_g', '_minp']))
    # rows: a header row per group, then one row per locus
    rows = []
    for g, title in GROUPS:
        sub = loci[loci.panel == g]
        if len(sub):
            rows.append(('header', title, None))
            rows += [('locus', r.locus_id, r) for _, r in sub.iterrows()]
    n = len(rows)

    W, H = 10.0, max(4.6, 0.92 + 0.27 * n)
    fig = plt.figure(figsize=(W, H))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 2.35])
    fig.subplots_adjust(left=0.62 / W, right=1 - 0.95 / W, top=1 - 0.78 / H, bottom=0.55 / H, wspace=0.62)
    ax_a, ax_b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])

    # a: composition, calibration and loci per cohort
    xc = np.arange(len(cohorts))
    for i, c in enumerate(cohorts):
        r = tail.loc[c]
        col = coh[c]['color']
        ax_a.bar(i, r.n_ctrl, bottom=r.n_case, width=0.62, color=tint(col, 0.35), zorder=3)
        ax_a.bar(i, r.n_case, width=0.62, color=col, zorder=3)
        ax_a.text(i, r.n_case / 2, f'{int(r.n_case)}', ha='center', va='center', color='white',
                  fontsize=plt.rcParams['legend.fontsize'], fontweight='bold', zorder=4)
        ax_a.text(i, r.n_case + r.n_ctrl / 2, f'{int(r.n_ctrl):,}', ha='center', va='center', color=S.INK,
                  fontsize=plt.rcParams['legend.fontsize'], zorder=4)
        ax_a.text(i, r.n_case + r.n_ctrl, f'{S.LAMBDA_GC} = {r.lambda_gc:.3f}\n'
                  f'{int(r.gw_loci)} GW $\\cdot$ {int(r.sg_loci)} SG', ha='center', va='bottom',
                  fontsize=plt.rcParams['legend.fontsize'] - 0.5, color=S.INK, zorder=4, linespacing=1.3)
    top = float((tail.n_case + tail.n_ctrl).max())
    ax_a.set_ylim(0, top * 1.22)
    ax_a.set_xlim(-0.6, len(cohorts) - 0.4)
    ax_a.set_xticks(xc)
    ax_a.set_xticklabels([coh[c]['label'].replace(' ', '\n', 1) for c in cohorts])
    ax_a.set_ylabel('samples')
    ax_a.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _p: f'{v:,.0f}'))
    S.despine(ax_a, grid_axis='y')
    S.panel_tag(ax_a, 'a', 'Cohorts', size=psize)
    ax_a.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=S.INK_SOFT, label='cases'),
                         plt.Rectangle((0, 0), 1, 1, color=tint('#454A50', 0.35), label='controls')],
                loc='upper left', frameon=False, fontsize=plt.rcParams['legend.fontsize'], ncol=2,
                handlelength=1.0, columnspacing=1.0, borderaxespad=0.1)

    # b: forest, one row per locus, the cohorts side by side
    off = np.linspace(0.24, -0.24, len(cohorts))
    out = []
    lo_all, hi_all = [], []
    for k, (kind, key, r) in enumerate(rows):
        y = -k
        if kind == 'header':
            ax_b.text(0.008, y, key, transform=ax_b.get_yaxis_transform(), ha='left', va='center',
                      fontsize=plt.rcParams['legend.fontsize'], fontstyle='italic', fontweight='bold',
                      color=S.INK_SOFT)
            if k:
                ax_b.axhline(y + 0.5, color=S.NEUTRAL_D, lw=0.7, zorder=1)
            continue
        if k % 2:
            ax_b.axhspan(y - 0.5, y + 0.5, color='#F5F6F7', zorder=0, lw=0)
        sub = x[x.locus_id == key].set_index('cohort')
        for j, c in enumerate(cohorts):
            if c not in sub.index or not np.isfinite(sub.loc[c, 'OR']):
                continue
            s = sub.loc[c]
            col = coh[c]['color']
            ax_b.plot([s.L95, s.U95], [y + off[j]] * 2, color=col, lw=1.4, solid_capstyle='butt', zorder=3)
            gw = s.called_peak == 'genome_wide'
            ax_b.scatter(s.OR, y + off[j], s=34 if gw else 20, marker='D' if gw else 'o', color=col,
                         edgecolors='white', linewidths=0.6, zorder=4)
            lo_all.append(s.L95)
            hi_all.append(s.U95)
            out.append(dict(group=r.panel, locus_id=key, gene=r.Gene, variant_id=r.variant_id, rsID=r.rsID,
                            cohort=c, OR=s.OR, L95=s.L95, U95=s.U95, P=s.P, called_peak=s.called_peak))
        gene = r.Gene if str(r.Gene) not in ('.', 'nan', '') else r.variant_id
        rs = r.rsID if str(r.rsID) not in ('.', 'nan', '') else r.variant_id
        ax_b.text(-0.012, y + 0.12, gene, transform=ax_b.get_yaxis_transform(), ha='right', va='center',
                  fontsize=plt.rcParams['ytick.labelsize'], fontstyle='italic', fontweight='bold', color=S.INK)
        ax_b.text(-0.012, y - 0.22, rs, transform=ax_b.get_yaxis_transform(), ha='right', va='center',
                  fontsize=plt.rcParams['ytick.labelsize'] - 1.5, color=S.INK_SOFT)
        best = sub.P.idxmin()
        ax_b.text(1.015, y, S.p_tex(float(sub.loc[best, 'P'])), transform=ax_b.get_yaxis_transform(),
                  ha='left', va='center', fontsize=plt.rcParams['legend.fontsize'], color=coh[best]['color'])
    ax_b.axvline(1.0, color=S.REFERENCE, lw=1.0, zorder=2)
    ax_b.set_xscale('log')
    ax_b.set_xlim(min(lo_all) / 1.15, max(hi_all) * 1.15)
    S.or_log_axis(ax_b, 'x')
    ax_b.set_ylim(-n + 0.5, 0.5)
    ax_b.set_yticks([])
    ax_b.set_xlabel(f'OR (95% CI), log scale')
    S.despine(ax_b, grid_axis='x')
    ax_b.spines['left'].set_visible(False)
    ax_b.text(1.015, 1.0, 'smallest $P$\n(its cohort)', transform=ax_b.transAxes, ha='left', va='bottom',
              fontsize=plt.rcParams['legend.fontsize'] - 0.5, color=S.INK_SOFT)
    S.panel_tag(ax_b, 'b', 'Loci in the main analysis', size=psize)
    handles = [Line2D([], [], color=coh[c]['color'], marker='o', lw=1.4, ms=4.5, label=coh[c]['label'])
               for c in cohorts]
    handles.append(Line2D([], [], color=S.INK_SOFT, marker='D', ls='', ms=5, label='genome-wide in that cohort'))
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 1 - 0.05 / H), ncol=len(handles),
               frameon=False, fontsize=plt.rcParams['legend.fontsize'] + 0.5, handlelength=1.6, columnspacing=1.6)
    fig.savefig(a.out_png)
    plt.close(fig)
    pd.DataFrame(out).to_csv(a.out_loci, sep='\t', index=False, float_format='%.6g')
    C.log(f'main summary: {len(loci)} loci -> {a.out_png}')


if __name__ == '__main__':
    main()
