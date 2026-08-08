#!/usr/bin/env python3
# =============================================================================
# pc_source_viz.py
# -----------------------------------------------------------------------------
# Purpose : PC-source visualization for covariate selection (PC-only, per set).
#           Compares 3 PC sources (geno_pc / bbj_pc / bbj_mainland_pc) by the
#           McFadden pseudo-R^2 of PC-only logistic models vs an intercept-only
#           baseline, across k = 0..MAX_PC, for each mainland sample set.
# Project : cteph_agp3k.v6  (tuning.pc analysis; used by tuning.pc.rev1.ipynb)
# Note    : Visualization only — this module does NOT choose the number of PCs,
#           and age/sex are NOT part of the analysis.
# =============================================================================
"""Reusable data-loading, model-fitting and plotting helpers for the PC-source
comparison notebook. The notebook should only import and call these functions."""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import AutoMinorLocator
from sklearn.linear_model import LogisticRegression

# ── Configuration ─────────────────────────────────────────────────────────────
RESULTS_DEFAULT = Path(
    '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/12_model_inputs'
)
DEFAULT_SETS = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']
SOURCES = ['geno_pc', 'bbj_pc', 'bbj_mainland_pc']
SRC_FILES = {
    'geno_pc':         'geno_pc.sex.tsv',
    'bbj_pc':          'bbj_pc.sex.tsv',
    'bbj_mainland_pc': 'bbj_mainland_pc.sex.tsv',
}
SOURCE_LABELS = {
    'geno_pc':         'Genotype (in-sample) PCA',
    'bbj_pc':          'BBJ projection',
    'bbj_mainland_pc': 'BBJ-mainland projection',
}

ID_COLS   = ['#FID', 'IID']
PHENO_COL = 'PHENO1'
MAX_PC    = 20

# ── Plot styling ──────────────────────────────────────────────────────────────
COLORS = {'geno_pc': '#1B7F79', 'bbj_pc': '#222222', 'bbj_mainland_pc': '#B36A2E'}
MARKERS = {'geno_pc': 'D', 'bbj_pc': 'o', 'bbj_mainland_pc': 's'}
HEAT_LABELS = {
    'geno_pc':         r'$\sum_{j=1}^{k} PC^{\mathrm{GT}}_{j}$',
    'bbj_pc':          r'$\sum_{j=1}^{k} PC^{\mathrm{BBJ}}_{j}$',
    'bbj_mainland_pc': r'$\sum_{j=1}^{k} PC^{\mathrm{mBBJ}}_{j}$',
}
_DIVERGING = LinearSegmentedColormap.from_list('muted_diverging', ['#567A9A', '#F7F3ED', '#A8564A'], N=256)


# ── Data loading ──────────────────────────────────────────────────────────────
def load_sets(results_dir=RESULTS_DEFAULT, sets=None, src_files=None):
    """Return (cov, pheno) dicts: cov[set][source] -> DataFrame, pheno[set] -> DataFrame."""
    results_dir = Path(results_dir)
    sets = list(sets) if sets is not None else list(DEFAULT_SETS)
    src_files = dict(src_files) if src_files is not None else dict(SRC_FILES)

    cov, pheno = {}, {}
    for s in sets:
        base = results_dir / s
        cov[s] = {src: pd.read_csv(base / 'covariates' / fn, sep='\t') for src, fn in src_files.items()}
        pheno[s] = pd.read_csv(base / 'phenotype' / 'pheno.tsv', sep='\t')
    return cov, pheno


# ── Model fitting ─────────────────────────────────────────────────────────────
def merge_cov_pheno(cov_df, pheno_df):
    """Inner-join covariate and phenotype tables on ID columns."""
    return cov_df.merge(pheno_df[ID_COLS + [PHENO_COL]], on=ID_COLS, how='inner')


def detect_and_encode_pheno(df, pheno_col=PHENO_COL):
    """Return a binary 0/1 phenotype array (handles PLINK 1/2 coding and 0/1)."""
    y = df[pheno_col].values.copy()
    if set(y) <= {1, 2}:
        y = (y == 2).astype(int)
    return y


def get_pc_cols(df):
    """Return PC columns (start with 'PC', end with '_AVG'), sorted numerically."""
    cols = [c for c in df.columns if c.startswith('PC') and c.endswith('_AVG')]
    return sorted(cols, key=lambda c: int(c[2:].replace('_AVG', '')))


def _loglik_intercept_only(y):
    p = np.clip(np.mean(y), 1e-15, 1 - 1e-15)
    return np.sum(y * np.log(p) + (1 - y) * np.log(1 - p))


def fit_model_intercept_null(df, feature_cols):
    """McFadden pseudo-R^2 of a PC-only logistic model vs the intercept-only baseline."""
    y = detect_and_encode_pheno(df, PHENO_COL)
    ll_null = _loglik_intercept_only(y)
    if len(feature_cols) == 0:
        return 0.0
    X = df[feature_cols].values
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        clf = LogisticRegression(C=1e12, solver='lbfgs', max_iter=1000).fit(X, y)
        p = clf.predict_proba(X)
        ll_full = np.sum(np.log(np.clip(p[np.arange(len(y)), y], 1e-15, 1)))
    return 1.0 - ll_full / ll_null


def build_curve_intercept_null(df, max_pc=MAX_PC):
    """For k = 0..max_pc fit phenotype ~ first-k PCs (no other covariates)."""
    pc_cols = get_pc_cols(df)[:max_pc]
    rows, prev_r2 = [], None
    for k in range(0, max_pc + 1):
        r2 = fit_model_intercept_null(df, pc_cols[:k])
        step = (r2 - prev_r2) if prev_r2 is not None else float('nan')
        rows.append({'k_included': k, 'r2': r2, 'step_gain': step})
        prev_r2 = r2
    return pd.DataFrame(rows)


def build_all_curves(cov, pheno, sets=None, sources=None, max_pc=MAX_PC):
    """Return (merged, CURVES): merged[set][src] joined frame, CURVES[set][src] curve DataFrame."""
    sets = list(sets) if sets is not None else list(cov.keys())
    sources = list(sources) if sources is not None else list(SOURCES)
    merged = {s: {src: merge_cov_pheno(cov[s][src], pheno[s]) for src in sources} for s in sets}
    curves = {s: {src: build_curve_intercept_null(merged[s][src], max_pc) for src in sources} for s in sets}
    return merged, curves


def summary_table(curves, merged, sets=None, sources=None, k_grid=(1, 5, 10, 15, 20)):
    """Return (summary_df at selected k, curves_long) — descriptive, no PC-count decision."""
    sets = list(sets) if sets is not None else list(curves.keys())
    sources = list(sources) if sources is not None else list(SOURCES)

    rows = []
    for s in sets:
        for src in sources:
            r2_by_k = curves[s][src].set_index('k_included')['r2']
            row = {'set': s, 'source': src, 'n': merged[s][src].shape[0]}
            for k in k_grid:
                row[f'R2@{k}'] = round(float(r2_by_k.loc[k]), 4)
            rows.append(row)
    summary_df = pd.DataFrame(rows)

    long = []
    for s in sets:
        for src in sources:
            for _, r in curves[s][src].iterrows():
                long.append({
                    'set': s, 'source': src, 'k': int(r['k_included']),
                    'r2': float(r['r2']),
                    'step_gain': float(r['step_gain']) if pd.notna(r['step_gain']) else np.nan,
                })
    return summary_df, pd.DataFrame(long)


# ── Plotting ──────────────────────────────────────────────────────────────────
def _apply_rcparams(fs):
    plt.style.use('default')
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': fs['tick'],
        'axes.labelsize': fs['label'], 'xtick.labelsize': fs['tick'], 'ytick.labelsize': fs['tick'],
        'legend.fontsize': fs['legend'], 'axes.linewidth': 1.2,
        'xtick.major.width': 1.2, 'ytick.major.width': 1.2, 'xtick.minor.width': 0.7, 'ytick.minor.width': 0.7,
    })


def compute_shared_scales(curves, sets=None, sources=None, pad_frac=0.04):
    """Global y-range (R^2 curves) and colour scale (|ΔR^2|) across all sets/sources.

    Returns (ylim, vabs) to pass to plot_pc_source_comparison so every set uses one
    coordinate system, making the figures directly comparable.
    """
    sets = list(curves.keys()) if sets is None else list(sets)
    sources = list(SOURCES) if sources is None else list(sources)
    r2_all, step_all = [], []
    for s in sets:
        for src in sources:
            r2_all.append(curves[s][src]['r2'].to_numpy())
            step_all.append(curves[s][src]['step_gain'].to_numpy())
    r2_all = np.concatenate(r2_all); r2_all = r2_all[np.isfinite(r2_all)]
    step_all = np.concatenate(step_all); step_all = step_all[np.isfinite(step_all)]
    lo, hi = float(r2_all.min()), float(r2_all.max())
    span = (hi - lo) or 1.0
    ylim = (lo - pad_frac * span, hi + pad_frac * span)
    vabs = float(np.nanmax(np.abs(step_all))) or 1.0
    return ylim, vabs


def plot_pc_source_comparison(set_id, curves, n_samples, sources=None, max_pc=MAX_PC,
                              annotate='salient', annot_frac=0.20, ylim=None, vabs=None, show=True):
    """PPT-ready figure: cumulative R^2(k) curves + per-PC Delta-R^2 heatmap for one set.

    The legend is placed ABOVE the panels (outside the axes) so it never overlaps the curves.

    annotate : 'salient' (default) prints ΔR^2 only in cells with |ΔR^2| >= annot_frac * max|ΔR^2|
               (the rest is shown by colour + colorbar — clean for slides);
               'none' prints no numbers (colour only); 'all' prints every cell.
    ylim     : (lo, hi) shared y-range for the R^2 curves (None = auto per-set).
    vabs     : shared colour scale for the heatmap, i.e. vmin/vmax = -/+vabs (None = auto per-set).
               Use compute_shared_scales(CURVES) to get comparable (ylim, vabs) across all sets.
    """
    sources = list(sources) if sources is not None else list(SOURCES)
    fs = {'title': 24, 'label': 17.5, 'tick': 15, 'legend': 15, 'heat': 11, 'ytick_heat': 16}
    _apply_rcparams(fs)

    # Slightly wide, slide-friendly canvas; extra top margin reserved for title + legend.
    fig = plt.figure(figsize=(13.5, 9.4), facecolor='white')
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.65, 1.35], width_ratios=[1, 0.03],
                  hspace=0.07, wspace=0.02, top=0.80, left=0.135, right=0.965, bottom=0.105)
    ax1 = fig.add_subplot(gs[0, 0]); ax3 = fig.add_subplot(gs[1, 0], sharex=ax1); cax = fig.add_subplot(gs[1, 1])
    for ax in (ax1, ax3):
        ax.set_facecolor('#FCFCFA')

    # Upper: cumulative pseudo-R²
    for s in sources:
        g = curves[s].sort_values('k_included')
        ax1.plot(g['k_included'], g['r2'], color=COLORS[s], linewidth=2.8,
                 marker=MARKERS[s], markersize=7.5, markevery=2,
                 markeredgecolor='white', markeredgewidth=1.0, label=SOURCE_LABELS[s])
    ax1.set_ylabel("McFadden pseudo-$R^{2}$\n(vs intercept-only)", labelpad=8, linespacing=1.3)
    ax1.set_axisbelow(True)
    ax1.grid(axis='y', alpha=0.18, linewidth=0.7, linestyle=':', color='#9E9E9E')
    ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax1.tick_params(direction='in', length=5.0, width=1.1, top=True, right=True, which='major', colors='#333333')
    ax1.tick_params(direction='in', length=2.6, width=0.7, top=True, right=True, which='minor', colors='#666666')
    ax1.set_xlim(-0.5, max_pc + 0.5)
    if ylim is not None:
        ax1.set_ylim(*ylim)
    plt.setp(ax1.get_xticklabels(), visible=False)

    # Lower: marginal-gain heatmap (one row per source)
    mat, rl = [], []
    for s in sources:
        step = curves[s].set_index('k_included')['step_gain'].reindex(range(1, max_pc + 1)).fillna(0).values
        mat.append(step); rl.append(HEAT_LABELS[s])
    mat = np.array(mat)
    vloc = vabs if vabs is not None else (float(np.nanmax(np.abs(mat))) or 1.0)
    im = ax3.imshow(mat, aspect='auto', cmap=_DIVERGING, interpolation='nearest',
                    vmin=-vloc, vmax=vloc, extent=[0.5, max_pc + 0.5, len(sources) - 0.5, -0.5])
    if annotate != 'none':
        thresh = 0.0 if annotate == 'all' else annot_frac * vloc
        _salient = (annotate == 'salient')
        for ri in range(mat.shape[0]):
            for ci in range(mat.shape[1]):
                v = mat[ri, ci]
                if _salient and abs(v) < thresh:
                    continue
                col = '#F8F8F8' if abs(v) > 0.55 * vloc else '#2B2B2B'
                ax3.text(ci + 1, ri, f'{v + 0.0:.3f}', ha='center', va='center',
                         fontsize=fs['heat'] + (2 if _salient else 0),
                         fontweight='bold' if _salient else 'normal',
                         color=col, fontfamily='monospace', zorder=4)
    for sep in [i + 0.5 for i in range(len(sources))]:
        ax3.axhline(sep, color='#C7C7C7', linewidth=0.55, zorder=3, alpha=0.82)
    ax3.set_xlabel('Number of Principal Components Added ($k$)', labelpad=8)
    ax3.set_ylabel('PC source', labelpad=10)
    ax3.set_yticks(range(len(sources)))
    ytl = ax3.set_yticklabels(rl, fontsize=fs['ytick_heat'])
    for i, s in enumerate(sources):
        ytl[i].set_color(COLORS[s]); ytl[i].set_fontweight('bold')
    ax3.tick_params(axis='both', direction='in', length=5.0, width=1.1, top=True, right=True, colors='#333333')
    ax3.set_xlim(-0.5, max_pc + 0.5)
    ax3.set_xticks(np.arange(0, max_pc + 1, 1))

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(r'Incremental $\Delta R^{2}$ per added PC', rotation=90, labelpad=10, fontsize=fs['label'] - 1)
    cbar.ax.tick_params(labelsize=fs['tick'] - 1, colors='#4A4A4A')

    # Legend ABOVE the panels (outside the axes → never overlaps the curves); title above the legend.
    fig.canvas.draw()
    pos = ax1.get_position()
    cx = pos.x0 + pos.width / 2
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(cx, pos.y1 + 0.020),
               ncol=len(sources), frameon=True, framealpha=0.97, edgecolor='#D1D1D1',
               facecolor='white', borderpad=0.5, columnspacing=1.8, handlelength=2.4, handletextpad=0.6)
    fig.suptitle(f'PC-source comparison  ·  {set_id}  (n={n_samples:,})  ·  PC-only',
                 x=cx, y=0.965, fontsize=fs['title'], fontweight='bold')
    if show:
        plt.show()
    return fig
