#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Flag depth-confounded outlier samples for the rare-variant tuning
#           pre-check. Robustly regress the per-sample minor-allele burden on
#           mean depth (Theil-Sen), take the residual, and mark samples whose
#           residual robust-Z exceeds a threshold. Writes a PLINK2 --remove
#           list (FID<TAB>IID, no header) plus a diagnostic figure coloured by
#           case/control, target depth, and (if a platform file is given) by
#           sequencing platform.
# Project : cteph_agp3k.v6 / tuning.rv  (rare-variant depth-confounding QC)
# Used by : tuning.rv.nf  process DETECT_OUTLIER_SAMPLES
# ---------------------------------------------------------------------------
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import plot_style as ps


def parse_args():
    p = argparse.ArgumentParser(
        description='Flag depth-confounded outlier samples (SMinAC ~ MeanDP residual robust-Z).')
    p.add_argument('--sample-metrics', required=True,
                   help='sample_metrics.txt[.gz] from run_calc_metrics.sh (minAC 0, raw).')
    p.add_argument('--z-threshold', type=float, default=5.0,
                   help='|robust z| cutoff for residual outlier marking (default: 5.0).')
    p.add_argument('--out-remove', required=True,
                   help='Output PLINK2 --remove list (FID<TAB>IID, no header).')
    p.add_argument('--out-png', required=True,
                   help='Output diagnostic figure (PNG, 600 dpi).')
    p.add_argument('--platform-file', default=None,
                   help='Optional CSV/XLSX mapping sample ID -> platform (adds a by-platform panel).')
    p.add_argument('--platform-id-col', default='ID_JHRPv6',
                   help='ID column in --platform-file (default: ID_JHRPv6).')
    p.add_argument('--platform-col', default='WGS_Platform',
                   help='Platform column in --platform-file (default: WGS_Platform).')
    return p.parse_args()


def main():
    args = parse_args()
    ROBUST_Z_THRESHOLD = args.z_threshold

    # ── Load ────────────────────────────────────────────────────────────────
    df = pd.read_csv(args.sample_metrics, sep='\t')
    required_cols = ['MeanDP', 'SMinAC', 'Group', 'TargetDP',
                     'SNumHomRef', 'SNumHet', 'SNumHomAlt']
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        raise ValueError(f'Missing columns: {missing_cols}')

    # ── Burden RATE (per 1,000 callable sites); robust fit rate ~ MeanDP ──────
    # rate decouples burden from differential missingness (15x controls miss more).
    df_plot = df.copy()
    callable_n = (pd.to_numeric(df_plot['SNumHomRef'], errors='coerce')
                  + pd.to_numeric(df_plot['SNumHet'], errors='coerce')
                  + pd.to_numeric(df_plot['SNumHomAlt'], errors='coerce'))
    rate_num = pd.to_numeric(df_plot['SMinAC'], errors='coerce') / callable_n * 1000.0
    df_plot['rate_k'] = rate_num
    meandp_num = pd.to_numeric(df_plot['MeanDP'], errors='coerce')

    fit_mask = rate_num.notna() & meandp_num.notna()
    if fit_mask.sum() < 2:
        raise ValueError('Not enough finite (rate, MeanDP) pairs to fit a regression.')

    x_fit = meandp_num[fit_mask].to_numpy(dtype=float)
    y_fit = rate_num[fit_mask].to_numpy(dtype=float)
    try:
        from scipy.stats import theilslopes
        slope, intercept = theilslopes(y_fit, x_fit)[:2]   # robust (median-based) fit
        fit_kind = 'Theil-Sen'
    except Exception:
        slope, intercept = np.polyfit(x_fit, y_fit, 1)      # OLS fallback
        fit_kind = 'OLS'

    df_plot['rate_resid'] = rate_num - (slope * meandp_num + intercept)

    # ── Robust-z on the residual -> outlier flag ─────────────────────────────
    resid_num = df_plot['rate_resid']
    median_r = float(np.nanmedian(resid_num))
    mad_r = float(np.nanmedian(np.abs(resid_num - median_r)))
    if np.isnan(mad_r) or mad_r == 0:
        df_plot['resid_robust_z_abs'] = np.nan
    else:
        df_plot['resid_robust_z_abs'] = np.abs(0.6745 * (resid_num - median_r) / mad_r)

    df_plot['resid_is_outlier'] = df_plot['resid_robust_z_abs'] > ROBUST_Z_THRESHOLD
    outlier_df = df_plot[df_plot['resid_is_outlier']]
    print(f'Robust fit ({fit_kind}): rate = {slope:.4g} * MeanDP + {intercept:.4g}')
    print(f'Residual robust z outliers (|z| > {ROBUST_Z_THRESHOLD:g}): {len(outlier_df)}')

    # ── Export outliers as FID<TAB>IID (FID == IID, no header) ───────────────
    id_col = next((c for c in ('SampleID', 'IID', 'ID') if c in df_plot.columns),
                  df_plot.columns[0])
    outlier_ids = outlier_df[id_col].astype(str)
    with open(args.out_remove, 'w') as f:
        for sid in outlier_ids:
            f.write(f'{sid}\t{sid}\n')
    print(f'ID column used : {id_col}')
    print(f'Outlier file   : {args.out_remove}  ({len(outlier_ids)} samples)')

    # ── Optional platform annotation (adds the by-platform panel) ─────────────
    if args.platform_file:
        pf = args.platform_file
        pmap = (pd.read_excel(pf, dtype={args.platform_id_col: str}) if pf.lower().endswith(('.xlsx', '.xls'))
                else pd.read_csv(pf, dtype={args.platform_id_col: str}))
        if args.platform_id_col in pmap.columns and args.platform_col in pmap.columns:
            m = dict(zip(pmap[args.platform_id_col].astype(str).str.strip(),
                         pmap[args.platform_col].astype(str).str.strip()))
            df_plot['Platform'] = df_plot[id_col].astype(str).str.strip().map(m).fillna('NA')
            n_mapped = int((df_plot['Platform'] != 'NA').sum())
            print(f'Platform       : {args.platform_file}  (mapped {n_mapped}/{len(df_plot)})')
        else:
            print(f'[warn] platform file lacks {args.platform_id_col}/{args.platform_col}; skipping platform panel')

    # ── Color maps ───────────────────────────────────────────────────────────
    group_color_map = ps.group_palette(df_plot['Group'].fillna('NaN').astype(str).unique())

    targetdp_is_numeric = pd.api.types.is_numeric_dtype(df_plot['TargetDP'])
    if not targetdp_is_numeric:
        targetdp_color_map, _ = ps.target_palette(df_plot['TargetDP'].fillna('NaN').astype(str).unique())

    has_platform = 'Platform' in df_plot.columns
    if has_platform:
        platform_color_map = ps.categorical_palette(sorted(df_plot['Platform'].astype(str).unique()))

    # ── Publication-grade plot style (shared module) ─────────────────────────
    ps.setup_style('slide')
    plt.rcParams.update({'legend.title_fontsize': plt.rcParams['legend.fontsize'] + 0.5})

    # Transparency is the main lever against overplotting (many samples pile up
    # near residual 0): keep markers near their original size but lower alpha so
    # the dense band shows density instead of a solid blob.
    SCATTER_KW = dict(s=22, alpha=0.38, linewidths=0.0, zorder=3, rasterized=True)
    OUTLIER_KW = dict(s=95, facecolors='none', edgecolors='#1A1A1A', linewidths=1.4, zorder=6)
    x_lbl = 'mean depth'
    # Single-line label: the two-line form was tall enough that its top reached the
    # suptitle and overlapped it. Units and the "depth-adjusted" qualifier now live in
    # the caption's Definitions block, where the symbol is defined anyway.
    y_lbl = 'burden-rate residual'
    y_col = 'rate_resid'

    _x = meandp_num[np.isfinite(meandp_num)]
    _y = resid_num[np.isfinite(resid_num)]
    x_pad, y_pad = 0.04 * (_x.max() - _x.min()), 0.06 * (_y.max() - _y.min())
    XLIM = (_x.min() - x_pad, _x.max() + x_pad)
    YLIM = (_y.min() - y_pad, _y.max() + y_pad)

    band = (ROBUST_Z_THRESHOLD * mad_r / 0.6745) if (mad_r and not np.isnan(mad_r)) else np.nan
    lo, hi = median_r - band, median_r + band

    # ── Figure ────────────────────────────────────────────────────────────────
    # All panels show the SAME residual against the SAME depth axis, differing only in
    # colouring, so they share both scales: one y-label instead of three, and the width
    # freed by dropping the repeated tick labels goes back to the data.
    n_panels = 3 if has_platform else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(5.5 * n_panels, 7.0), sharey=True)

    def _frame(ax):
        ps.despine(ax)
        ax.set_xlim(*XLIM)
        ax.set_ylim(*YLIM)
        if np.isfinite(band):
            ax.axhspan(lo, hi, color='#9AA0A6', alpha=0.10, zorder=1)
            ax.axhline(lo, color='#9AA0A6', linestyle=(0, (5, 4)), linewidth=1.0, zorder=2)
            ax.axhline(hi, color='#9AA0A6', linestyle=(0, (5, 4)), linewidth=1.0, zorder=2)
        ax.axhline(0.0, color='#444444', linewidth=1.2, zorder=2)
        ax.set_xlabel(x_lbl)

    def _scatter_ordered(ax, col, cmap):
        # Draw the largest category first so smaller ones land on top (not buried).
        order = df_plot[col].fillna('NaN').astype(str).value_counts().index.tolist()
        for lbl in order:
            sub = df_plot[df_plot[col].fillna('NaN').astype(str) == lbl]
            ax.scatter(sub['MeanDP'], sub[y_col], color=cmap[lbl], label=lbl, **SCATTER_KW)

    # Colour key sits INSIDE each panel's empty upper-right corner (the residuals
    # hug y=0 and the outliers are on the left, so that corner is clear) — this
    # keeps the wide platform key off the suptitle. The outlier-ring and
    # inlier-band meaning goes once into the caption, not repeated per panel.
    # (a) coloured by Group
    ax = axes[0]
    _frame(ax)
    _scatter_ordered(ax, 'Group', group_color_map)
    if not outlier_df.empty:
        ax.scatter(outlier_df['MeanDP'], outlier_df[y_col], **OUTLIER_KW)
    ax.set_ylabel(y_lbl)
    _hg = [Line2D([0], [0], marker='o', linestyle='', color=group_color_map[g],
                  markersize=7, label=g) for g in group_color_map]
    ps.legend_inside(ax, _hg, loc='upper right', title='case / control group')

    # (b) coloured by TargetDP
    ax = axes[1]
    _frame(ax)
    if targetdp_is_numeric:
        sc = ax.scatter(df_plot['MeanDP'], df_plot[y_col], c=df_plot['TargetDP'],
                        cmap='viridis', **SCATTER_KW)
        cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
        cbar.set_label('Target depth', fontsize=11)
        cbar.ax.tick_params(labelsize=10, length=3, width=0.8)
        cbar.outline.set_linewidth(0.8)
        ax.set_title('target depth', loc='right', fontsize=plt.rcParams['legend.fontsize'])
    else:
        _scatter_ordered(ax, 'TargetDP', targetdp_color_map)
        _th = [Line2D([0], [0], marker='o', linestyle='', color=targetdp_color_map[t],
                      markersize=7, label=t) for t in targetdp_color_map]
        ps.legend_inside(ax, _th, loc='upper right', title='target depth')
    if not outlier_df.empty:
        ax.scatter(outlier_df['MeanDP'], outlier_df[y_col], **OUTLIER_KW)
    ax.set_ylabel(y_lbl)

    # (c) coloured by sequencing platform (optional)
    if has_platform:
        ax = axes[2]
        _frame(ax)
        _scatter_ordered(ax, 'Platform', platform_color_map)
        if not outlier_df.empty:
            ax.scatter(outlier_df['MeanDP'], outlier_df[y_col], **OUTLIER_KW)
        ax.set_ylabel(y_lbl)
        _ph = [Line2D([0], [0], marker='o', linestyle='', color=platform_color_map[p],
                      markersize=7, label=p) for p in platform_color_map]
        ps.legend_inside(ax, _ph, loc='upper right', title='sequencing platform')

    ps.share_y(axes, y_lbl)
    fig.suptitle('Sample QC — depth-adjusted burden-rate residual per sample', fontsize=16, y=0.97)
    for _ax, _L in zip(axes, 'abc'):
        ps.panel_tag(_ax, _L)

    _src = f'{Path(args.sample_metrics).name} (per-sample metrics, minAC 0)'
    if args.platform_file:
        _src += f'; platform from {Path(args.platform_file).name}'
    n_out = len(outlier_df)
    ps.caption_block(
        fig, top=0.905, wspace=0.10,
        title=('Sample-level outliers are flagged by a depth-adjusted burden-rate residual, which shows '
               'no case–control separation beyond the depth/platform trend.'),
        panels=[
            'per-sample depth-adjusted burden-rate residual, coloured by case/control group.',
            'the same residual, coloured by target sequencing depth (15× vs 30×).',
            'the same residual, coloured by sequencing platform.',
        ],
        interpret=(f'{n_out} sample(s) exceed |robust-$Z$| > {ROBUST_Z_THRESHOLD:g} (open rings) and are removed '
                   'before variant QC. The band is centred on zero with no phenotype separation; the residual '
                   'structure instead tracks target depth and platform (not case/control) — the fully confounded '
                   'design made visible, motivating the depth covariate and minAC sweep downstream.'),
        notes=(f'{_src}. $y$ = burden rate minus its Theil–Sen fit on mean depth (per $10^3$ callable sites); '
               f'open rings, |robust-$Z$| > {ROBUST_Z_THRESHOLD:g}; grey band, inlier region '
               f'($\\pm{ROBUST_Z_THRESHOLD:g}$ MAD). All three panels share one residual and depth scale.'),
        defs=['rate', 'robustz'],
        model=f'{ps.FORMULAS["rate"]}      {ps.FORMULAS["robustz"]}')
    fig.savefig(args.out_png)
    print(f'Figure         : {args.out_png}')


if __name__ == '__main__':
    main()
