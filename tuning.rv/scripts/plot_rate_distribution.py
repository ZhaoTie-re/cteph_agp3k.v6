#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Show what sample-QC does to the per-sample burden RATE and how that
#           rate is structured by the confounders. Four panels:
#             (a) rate distribution BEFORE vs AFTER outlier removal (density,
#                 with skewness / excess-kurtosis for each);
#             (b) rate by target depth (15x / 30x);
#             (c) rate by sequencing platform;
#             (d) rate by case/control group.
#           Because the design is fully confounded (cases 30x, controls 15x),
#           panels (b)-(d) should track each other — the apparent group effect
#           in (d) is the depth/platform effect in (b)/(c).
# Project : cteph_agp3k.v6 / tuning.rv  (rare-variant depth-confounding QC)
# Used by : tuning.rv.nf  process PLOT_RATE_DISTRIBUTION
# ---------------------------------------------------------------------------
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

import plot_style as ps


def parse_args():
    p = argparse.ArgumentParser(
        description='Burden-rate distribution before/after sample QC, and rate by depth/platform/group.')
    p.add_argument('--sample-metrics', required=True,
                   help='sample_metrics.txt[.gz] (minAC 0, raw) — the BEFORE-QC set.')
    p.add_argument('--removed-file', default=None,
                   help='PLINK2 --remove list (FID<TAB>IID) of QC-removed samples; AFTER = full minus these.')
    p.add_argument('--platform-file', default=None,
                   help='Optional CSV/XLSX mapping sample ID -> platform (adds the by-platform panel).')
    p.add_argument('--platform-id-col', default='ID_JHRPv6')
    p.add_argument('--platform-col', default='WGS_Platform')
    p.add_argument('--out-png', required=True, help='Output figure (PNG, 600 dpi).')
    return p.parse_args()


def _rate_k(df):
    """Per-sample burden rate = SMinAC / callable * 1000 (callable = non-missing genotypes)."""
    callable_n = (pd.to_numeric(df['SNumHomRef'], errors='coerce')
                  + pd.to_numeric(df['SNumHet'], errors='coerce')
                  + pd.to_numeric(df['SNumHomAlt'], errors='coerce'))
    return pd.to_numeric(df['SMinAC'], errors='coerce') / callable_n * 1000.0


def _load_removed(path):
    """First whitespace column of a headerless --remove list -> set of IDs (empty if unusable)."""
    try:
        r = pd.read_csv(path, sep=r'\s+', header=None, dtype=str)
        return set(r[0].astype(str).str.strip())
    except Exception:
        return set()


def _wrap_label(text, width=13):
    """Break a long category name onto two lines at a space or a hyphen.

    Platform names ('DNBSeq-G400RS 30x') are too wide to sit side by side upright,
    and rotating them made them collide instead; wrapping keeps them horizontal.
    """
    s = str(text)
    if len(s) <= width:
        return s
    cut = max(s.rfind(' ', 0, width + 1), s.rfind('-', 0, width + 1))
    if cut <= 0:
        return s
    return s[:cut].rstrip() + '\n' + s[cut:].lstrip('- ').strip()


def _violin(ax, groups, labels, colors, ylabel, rotate=0):
    """Violin + slim box + median dot for a list of 1-D arrays; annotate n and median.

    rotate is retained for call compatibility and now only selects a slightly smaller
    tick font for the many-level (platform) panel; labels are never rotated.
    """
    groups = [np.asarray(g, dtype=float) for g in groups]
    groups = [g[np.isfinite(g)] for g in groups]
    pos = np.arange(1, len(groups) + 1)
    keep = [i for i, g in enumerate(groups) if g.size]
    if keep:
        vp = ax.violinplot([groups[i] for i in keep], positions=pos[keep],
                           showextrema=False, widths=0.78)
        for body, i in zip(vp['bodies'], keep):
            body.set_facecolor(colors[i]); body.set_edgecolor(colors[i])
            body.set_alpha(0.35); body.set_linewidth(1.0)
        bp = ax.boxplot([groups[i] for i in keep], positions=pos[keep], widths=0.14,
                        showfliers=False, patch_artist=True, medianprops=dict(color='white', linewidth=1.4),
                        whiskerprops=dict(color='#444444', linewidth=1.0),
                        capprops=dict(color='#444444', linewidth=1.0),
                        boxprops=dict(linewidth=0.0))
        for patch, i in zip(bp['boxes'], keep):
            patch.set_facecolor(colors[i]); patch.set_alpha(0.95)
        for i in keep:
            med = float(np.median(groups[i]))
            # Value printed AT the median line (just right of the IQR box) so it
            # unambiguously labels the median — not the top of the violin. A light
            # halo keeps it readable over the body.
            # Printed AT the median line, just right of the IQR box, so it labels the
            # median unambiguously rather than the top of the violin. A light halo
            # keeps it readable where it crosses the body.
            ax.annotate(f'{med:.2f}', (pos[i] + 0.09, med), va='center', ha='left',
                        fontsize=9, color='#111111', zorder=8,
                        bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                                  edgecolor='none', alpha=0.65))
    ax.set_xticks(pos)
    # Always two lines, always upright. Rotated single-line platform labels
    # ('DNBSeq-G400RS 30x (n=39)' at 22 deg) overlapped each other by up to
    # 23,000 px^2 and ran off the left edge of the canvas; stacking the count under
    # a wrapped name keeps every label horizontal, readable and inside its slot.
    _w = 10 if rotate else 14          # many-level (platform) panel needs narrower labels
    ax.set_xticklabels([f'{_wrap_label(l, _w)}\n(n={g.size:,})' for l, g in zip(labels, groups)],
                       fontsize=plt.rcParams['xtick.labelsize'] - (2.0 if rotate else 0))
    ax.set_ylabel(ylabel, labelpad=6)
    ps.despine(ax, grid_axis='y')


def main():
    args = parse_args()
    ps.setup_style('slide')
    RATE_LBL = f'burden rate {ps.RATE}\n(per 1,000 callable sites)'

    df = pd.read_csv(args.sample_metrics, sep='\t')
    need = ['SampleID', 'Group', 'TargetDP', 'SNumHomRef', 'SNumHet', 'SNumHomAlt', 'SMinAC']
    miss = [c for c in need if c not in df.columns]
    if miss:
        raise ValueError(f'sample metrics missing columns: {miss}')
    df['SampleID'] = df['SampleID'].astype(str).str.strip()
    df['rate_k'] = _rate_k(df)

    # BEFORE vs AFTER outlier removal
    removed = _load_removed(args.removed_file) if args.removed_file else set()
    df['removed'] = df['SampleID'].isin(removed)
    before = df['rate_k'].dropna()
    after = df.loc[~df['removed'], 'rate_k'].dropna()
    n_removed = int(df['removed'].sum())

    # Optional platform annotation
    has_platform = False
    if args.platform_file:
        pf = args.platform_file
        pmap = (pd.read_excel(pf, dtype={args.platform_id_col: str}) if pf.lower().endswith(('.xlsx', '.xls'))
                else pd.read_csv(pf, dtype={args.platform_id_col: str}))
        if args.platform_id_col in pmap.columns and args.platform_col in pmap.columns:
            m = dict(zip(pmap[args.platform_id_col].astype(str).str.strip(),
                         pmap[args.platform_col].astype(str).str.strip()))
            df['Platform'] = df['SampleID'].map(m).fillna('NA')
            has_platform = True

    aft = df[~df['removed']].copy()          # analysis-ready set for the by-stratum panels

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 11.0))
    ax_a, ax_b, ax_c, ax_d = axes.flat

    # ── (a) rate distribution: before vs after ───────────────────────────────
    # Zoom to the bulk (removed outliers sit beyond the right edge; the tail they
    # form is quantified by the skewness / kurtosis annotation, not shown raw).
    lo = float(np.nanpercentile(before, 0.5))
    hi = float(np.nanpercentile(before, 99.7))
    span = hi - lo
    lo, hi = lo - 0.03 * span, hi + 0.03 * span
    bins = np.linspace(lo, hi, 46)
    ax_a.hist(before, bins=bins, density=True, color='#BBBBBB', alpha=0.85, zorder=2)
    ax_a.hist(after, bins=bins, density=True, histtype='step', color=ps.GROUP_COLOR,
              linewidth=2.2, zorder=4)
    for series, col in ((before, '#888888'), (after, ps.GROUP_COLOR)):
        ax_a.axvline(float(np.median(series)), color=col, linestyle='--', linewidth=1.3, zorder=5)
    sk_b, ku_b = float(stats.skew(before, nan_policy='omit')), float(stats.kurtosis(before, nan_policy='omit'))
    sk_a, ku_a = float(stats.skew(after, nan_policy='omit')), float(stats.kurtosis(after, nan_policy='omit'))
    # Colour key INSIDE the empty upper-right (the density is low in the tails);
    # skew/kurtosis reported in the caption. Frees the title slot for the panel letter.
    _ha = [Patch(facecolor='#BBBBBB', alpha=0.85, label=f'before QC (n = {before.size:,})'),
           Line2D([0], [0], color=ps.GROUP_COLOR, linewidth=2.4, label=f'after QC (n = {after.size:,})')]
    # Stays INSIDE the empty upper-right corner: the distribution's right tail is flat
    # there, and the title slot above is already taken by the panel letter + title.
    ps.legend_inside(ax_a, _ha, loc='upper right')
    ax_a.set_xlim(lo, hi)
    ax_a.set_xlabel(RATE_LBL.replace('\n', ' '), labelpad=6)
    ax_a.set_ylabel('density', labelpad=6)
    ps.despine(ax_a)

    # ── (b) rate by target depth ─────────────────────────────────────────────
    dpal, dorder = ps.target_palette(aft['TargetDP'].dropna().unique())
    dgroups = [aft.loc[aft['TargetDP'].astype(str) == lv, 'rate_k'].values for lv in dorder]
    _violin(ax_b, dgroups, [f'{lv}' for lv in dorder], [dpal[lv] for lv in dorder], RATE_LBL)
    ps.panel_tag(ax_b, 'b', 'Rate by target depth')   # x-title dropped — ticks name the levels

    # ── (c) rate by platform (or fallback note) ──────────────────────────────
    if has_platform:
        plats = (aft[aft['Platform'] != 'NA'].groupby('Platform')['rate_k'].median()
                 .sort_values(ascending=False).index.tolist())
        cpal = ps.categorical_palette(plats)
        cgroups = [aft.loc[aft['Platform'] == p, 'rate_k'].values for p in plats]
        _violin(ax_c, cgroups, plats, [cpal[p] for p in plats], RATE_LBL, rotate=22)
        ps.panel_tag(ax_c, 'c', 'Rate by platform')   # x-title dropped (tick labels are platform names)
    else:
        ax_c.axis('off')
        ax_c.text(0.5, 0.5, 'platform file not provided', transform=ax_c.transAxes,
                  ha='center', va='center', color='#888888', style='italic')

    # ── (d) rate by case/control group ───────────────────────────────────────
    gorder = [g for g in ('Control', 'Case') if g in set(aft['Group'])]
    gpal = ps.group_palette(gorder)
    ggroups = [aft.loc[aft['Group'] == g, 'rate_k'].values for g in gorder]
    _violin(ax_d, ggroups, gorder, [gpal[g] for g in gorder], RATE_LBL)
    ps.panel_tag(ax_d, 'd', 'Rate by group (apparent effect)')   # x-title dropped — ticks name the levels

    ps.panel_tag(ax_a, 'a', 'Rate before vs after QC')
    fig.suptitle('Burden-rate distribution and its confounding structure', fontsize=16, y=0.995)
    ps.caption_block(
        fig, top=0.90, wspace=0.24, hspace=0.62,
        title=('Sample QC restores a symmetric burden-rate distribution, and the apparent Case–Control shift '
               'mirrors the depth/platform split rather than phenotype.'),
        panels=[
            (f'burden-rate distribution before ($n$ = {before.size:,}) versus after ($n$ = {after.size:,}) sample '
             f'QC; QC removes {n_removed:,} sample(s), cutting skewness {sk_b:.2f}→{sk_a:.2f} and excess kurtosis '
             f'{ku_b:.1f}→{ku_a:.1f} (dashed lines, medians).'),
            'the same rate split by target sequencing depth (violin, density; box, IQR; printed value, median).',
            'the same rate split by sequencing platform.',
            'the same rate split by case/control group — the apparent phenotype contrast.',
        ],
        interpret=('the group medians in d track the depth/platform medians in b, c: because cases are 30× and '
                   'controls 15× with no platform overlap, the apparent Case–Control shift is a depth/platform '
                   'artifact, not a phenotype effect.'),
        notes=(f'{Path(args.sample_metrics).name} (per-sample rate, minAC 0)'
               + (f'; platform from {Path(args.platform_file).name}' if args.platform_file else '')
               + '. Skewness and excess (Fisher) kurtosis computed on the rate. Violin, kernel density; box, '
                 'IQR; white line and printed value, median; n per level on the axis.'),
        defs=['rate'],
        model=ps.FORMULAS['rate'])
    fig.savefig(args.out_png)
    print(f'Figure         : {args.out_png}')
    print(f'before/after n : {before.size:,} / {after.size:,}  (removed {n_removed:,})')
    print(f'skew  before/after : {sk_b:.3f} / {sk_a:.3f}')
    print(f'kurt  before/after : {ku_b:.3f} / {ku_a:.3f}')


if __name__ == '__main__':
    main()
