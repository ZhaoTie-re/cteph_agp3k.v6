#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Headline minAC-sweep figure (2x2, publication grade). Reads the
#           merged per-minAC QC table (rate + OLS) and draws, vs the minAC
#           threshold kappa:
#             (a) beta_g  — Case-vs-Control burden difference, full range (V-shaped)
#             (b) beta_g  — zoom on the trough (the decision region)
#             (c) beta_d  — residual depth effect (DIAGNOSTIC; replaces the old
#                           Spearman panel; from the SAME OLS as beta_g)
#             (d) variants retained — the power cost of raising kappa (log y)
#           beta_g is the DECISION target (pick kappa at its non-significant
#           trough); beta_d is a diagnostic only (do NOT select on it). The
#           recommended kappa is marked on every panel. no-PC = primary model,
#           +ancestry-PC = sensitivity. All annotation lives in the caption below;
#           the shared legend sits above the panels — nothing overlaps the data.
# Project : cteph_agp3k.v6 / tuning.rv  (rare-variant depth-confounding QC)
# Used by : tuning.rv.nf  process PLOT_SWEEP
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import plot_style as ps


def parse_args():
    p = argparse.ArgumentParser(description='Headline minAC-sweep figure (beta_g decision, beta_d diagnostic).')
    p.add_argument('--qc-stats', required=True, help='Merged per-minAC QC table (all_qc_summary_stats.tsv).')
    p.add_argument('--out-png', required=True)
    p.add_argument('--recommended-minac', type=int, default=None, help='minAC to mark (from select_minac).')
    p.add_argument('--alpha', type=float, default=0.05, help='Significance level for the markers.')
    p.add_argument('--include-zero', action='store_true', help='Keep minAC=0 in the sweep.')
    return p.parse_args()


def _series(ax, x, beta, lo, hi, p, color, alpha):
    """CI ribbon + line; filled markers where significant (p<alpha), open where not.

    The marker half is delegated to plot_style.sig_markers so this figure and the
    cross-cohort figure cannot drift apart on what "filled" means.
    """
    ok = beta.notna()
    ax.fill_between(x[ok], lo[ok], hi[ok], color=color, alpha=0.13, zorder=1)
    ax.plot(x[ok], beta[ok], color=color, linewidth=2.0, alpha=0.85, zorder=2)
    ps.sig_markers(ax, x, beta, p, color, alpha=alpha, size=52, zorder=4, ok=ok)


def _mark(ax, xr):
    if xr is not None:
        ax.axvline(xr, color='#444444', linestyle=':', linewidth=1.6, zorder=5)


def main():
    args = parse_args()
    ps.setup_style('slide')

    df = pd.read_csv(args.qc_stats, sep='\t')
    need = {'MinAC_Threshold', 'Variant_Count', 'Beta_group_noPC', 'Beta_group_noPC_L',
            'Beta_group_noPC_U', 'P_group_noPC', 'Beta_depth_noPC', 'Beta_depth_noPC_L',
            'Beta_depth_noPC_U', 'P_depth_noPC'}
    miss = need - set(df.columns)
    if miss:
        raise ValueError(f'qc stats missing required columns: {sorted(miss)}')
    if not args.include_zero:
        df = df[df['MinAC_Threshold'] != 0]
    df = df.sort_values('MinAC_Threshold').reset_index(drop=True)
    if df.empty:
        raise ValueError('No rows to plot (check --include-zero / input).')

    x = df['MinAC_Threshold']
    xs = sorted(df['MinAC_Threshold'].unique())
    has_pc = ('Beta_group_pc' in df.columns) and df['Beta_group_pc'].notna().any()
    pc_lab = str(df['PC_source'].dropna().iloc[0]) if ('PC_source' in df.columns
                                                       and df['PC_source'].notna().any()) else 'PC'
    n_pc = int(df['N_PC'].dropna().iloc[0]) if ('N_PC' in df.columns and df['N_PC'].notna().any()) else None
    pc_txt = f'{pc_lab}' + (f', {n_pc} PCs' if n_pc else '')
    xr = args.recommended_minac
    kx = 'minAC threshold (minAC)'

    fig, axes = plt.subplots(2, 2, figsize=(14.0, 12.5))
    axA, axB, axC, axD = axes.flat

    # ── (a) beta_group, full range ───────────────────────────────────────────
    _series(axA, x, df['Beta_group_noPC'], df['Beta_group_noPC_L'], df['Beta_group_noPC_U'],
            df['P_group_noPC'], ps.BETA_COLOR, args.alpha)
    if has_pc:
        _series(axA, x, df['Beta_group_pc'], df['Beta_group_pc_L'], df['Beta_group_pc_U'],
                df['P_group_pc'], ps.GROUP_COLOR, args.alpha)
    axA.axhline(0, color='#888888', linewidth=1.1, zorder=0)
    axA.set_ylabel(r'OLS $\beta_{\mathrm{group}}$  (Case $-$ Control)')
    axA.set_xlabel(kx)
    ps.panel_tag(axA, 'a', 'Apparent phenotype effect (V-shaped)')
    _mark(axA, xr)
    ps.despine(axA)

    # ── (b) beta_group, trough zoom (decision region) ────────────────────────
    _series(axB, x, df['Beta_group_noPC'], df['Beta_group_noPC_L'], df['Beta_group_noPC_U'],
            df['P_group_noPC'], ps.BETA_COLOR, args.alpha)
    if has_pc:
        _series(axB, x, df['Beta_group_pc'], df['Beta_group_pc_L'], df['Beta_group_pc_U'],
                df['P_group_pc'], ps.GROUP_COLOR, args.alpha)
    axB.axhline(0, color='#888888', linewidth=1.1, zorder=0)
    zoom_max = max(6, (xr or 2) + 4)
    win = df[df['MinAC_Threshold'] <= zoom_max]
    axB.set_xlim(df['MinAC_Threshold'].min() - 0.4, zoom_max + 0.4)
    clo = ['Beta_group_noPC_L'] + (['Beta_group_pc_L'] if has_pc else [])
    chi = ['Beta_group_noPC_U'] + (['Beta_group_pc_U'] if has_pc else [])
    ylo, yhi = float(win[clo].min().min()), float(win[chi].max().max())
    pad = 0.18 * (yhi - ylo + 1e-9)
    axB.set_ylim(ylo - pad, yhi + pad)
    axB.set_ylabel(r'OLS $\beta_{\mathrm{group}}$  (Case $-$ Control)')
    axB.set_xlabel(kx)
    ps.panel_tag(axB, 'b', 'Trough — the decision region')
    _mark(axB, xr)
    if xr is not None:
        # Bottom-right: beta_group climbs left-to-right out of the trough, so the
        # lower-right corner is the one region of the zoom that is reliably empty —
        # the old top-left position sat right against the recommendation line.
        axB.text(0.97, 0.06, f'recommended minAC = {xr}', transform=axB.transAxes,
                 ha='right', va='bottom', fontsize=11, color='#333333', fontweight='bold')
    axB.set_xticks([v for v in xs if v <= zoom_max])
    ps.despine(axB)

    # ── (c) beta_depth, full range (DIAGNOSTIC — replaces Spearman) ───────────
    _series(axC, x, df['Beta_depth_noPC'], df['Beta_depth_noPC_L'], df['Beta_depth_noPC_U'],
            df['P_depth_noPC'], ps.BETA_COLOR, args.alpha)
    if has_pc and 'Beta_depth_pc' in df.columns:
        _series(axC, x, df['Beta_depth_pc'], df['Beta_depth_pc_L'], df['Beta_depth_pc_U'],
                df['P_depth_pc'], ps.GROUP_COLOR, args.alpha)
    axC.axhline(0, color='#888888', linewidth=1.1, zorder=0)
    axC.set_ylabel(r'OLS $\beta_{\mathrm{depth}}$  (per $+1$ SD depth)')
    axC.set_xlabel(kx)
    ps.panel_tag(axC, 'c', 'Residual depth effect (diagnostic)')
    _mark(axC, xr)
    ps.despine(axC)

    # ── (d) power cost ───────────────────────────────────────────────────────
    axD.plot(x, df['Variant_Count'], color=ps.COUNT_COLOR, linewidth=2.2, zorder=2)
    axD.scatter(x, df['Variant_Count'], color=ps.COUNT_COLOR, s=42, zorder=3,
                edgecolor='white', linewidth=0.7)
    axD.set_yscale('log')
    axD.set_ylabel('variants retained')
    axD.set_xlabel(kx)
    ps.panel_tag(axD, 'd', 'Power cost of raising minAC')
    _mark(axD, xr)
    ps.despine(axD)

    # 1..20 fully labelled crowds a half-width panel; label a readable subset and
    # leave the rest as minor ticks so no threshold is hidden. Panel b is the zoom,
    # where every integer fits and is worth reading.
    for ax in (axA, axC, axD):
        ps.sparse_int_ticks(ax, xs)

    # Panel tags are folded into each title (loc='left') to avoid colliding with the
    # tall y-axis labels; no separate letter needed.

    # Shared legend ABOVE the panels (off data): model colours + marker meaning.
    handles = [Line2D([], [], color=ps.BETA_COLOR, lw=2.6, label='no-PC (primary)')]
    if has_pc:
        handles.append(Line2D([], [], color=ps.GROUP_COLOR, lw=2.6, label=f'+PC ({pc_txt}, sensitivity)'))
    handles += [Line2D([], [], marker='o', color='#666', markerfacecolor='#666', linestyle='',
                       markersize=8, label=f'significant (p < {args.alpha:g})'),
                Line2D([], [], marker='o', color='#666', markerfacecolor='white', linestyle='',
                       markersize=8, markeredgewidth=1.6, label='not significant')]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.935),
               ncol=len(handles), frameon=False, handletextpad=0.5, columnspacing=1.6)

    fig.suptitle('minAC sweep — controlling depth/platform inflation of the burden test',
                 fontsize=16, y=0.99)
    ps.caption_block(
        fig, top=0.86, wspace=0.26, hspace=0.90,
        title=('The depth-adjusted Case–Control effect is V-shaped in minAC; its non-significant trough at the '
               f'recommended minAC = {xr} minimises the artifact while preserving power.'
               if xr is not None else
               'The depth-adjusted Case–Control effect is V-shaped in minAC; its non-significant trough minimises '
               'the artifact while preserving power.'),
        panels=[
            (r'OLS $\beta_{\mathrm{group}}$, the depth- and PC-adjusted Case$-$Control burden difference, across the '
             'full minAC range — inflated at low minAC, a trough, then a steep climb.'),
            r'the same $\beta_{\mathrm{group}}$ zoomed on its trough (the decision region); dotted line, recommended minAC.',
            r'OLS $\beta_{\mathrm{depth}}$, the residual depth effect from the same model — a diagnostic, not a selection axis.',
            'variants retained versus minAC (power cost; log scale).',
        ],
        interpret=(r'$\beta_{\mathrm{group}}$ is significantly inflated at very low minAC, reaches a non-significant '
                   'trough (the recommended minAC), then climbs again as high-count platform-differential variants '
                   r'dominate — so no threshold removes the confounding entirely; the trough merely minimises it at '
                   r'acceptable power. Crucially $\beta_{\mathrm{depth}}$ does NOT bottom out at the same minAC — '
                   r'its minimum sits well inside the range where $\beta_{\mathrm{group}}$ is already large — so it '
                   'must not drive selection (see decision_axis.png). The genome-wide inflation factor λ is the '
                   'final calibration check.'),
        notes=(f'{Path(args.qc_stats).name} (merged per-minAC OLS on the burden rate); ancestry PC source '
               f'{pc_txt}. Ribbons, 95% CI; filled markers, p < {args.alpha:g}; open markers, non-significant. '
               'no-PC is the primary model, +PC an ancestry sensitivity check. minAC = 0 is omitted '
               '(singleton-dominated) unless --include-zero.'),
        defs=['minac', 'rate', 'beta_group', 'beta_depth', 'ci_sig', 'models', 'variants_retained', 'lambda'],
        model=ps.FORMULAS['ols'])
    fig.savefig(args.out_png)
    print(f'Figure         : {args.out_png}')
    if xr is not None:
        print(f'recommended minAC marked: {xr}')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in plot_sweep: {e}', file=sys.stderr)
        sys.exit(1)
