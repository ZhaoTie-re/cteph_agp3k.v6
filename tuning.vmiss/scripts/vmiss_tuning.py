#!/usr/bin/env python3
# =============================================================================
# vmiss_tuning.py
# -----------------------------------------------------------------------------
# Purpose : Per-variant call-rate (VMISS) QC-threshold tuning by depth subgroup
#           (15X / 30X). Reads the pipeline's variant-QC summary, builds the VMISS
#           cumulative distribution per depth group, and uses the Kneedle elbow
#           algorithm to recommend a VMISS threshold. Ports cteph_agp3k.v5/
#           tuning.vmiss to v6, reusing the pipeline summary (no plink2 recompute).
# Project : cteph_agp3k.v6  (tuning.vmiss analysis; used by tuning.vmiss.rev1.ipynb)
# Note    : Recommends thresholds only — does not edit vqc_config_vmiss.json.
# =============================================================================
"""Reusable load / CDF / Kneedle / plotting helpers for the VMISS-tuning notebook.
The notebook should only import and call these functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ── Configuration ─────────────────────────────────────────────────────────────
# Per-variant VMISS summary produced by the v6 pipeline (RUN_VARIANT_QC) on the
# 07_sample_qc genotype: 37,209,295 variants; has VMISS/30X_VMISS/15X_VMISS/CASE/CTRL/MAF.
SUMMARY_DEFAULT = Path(
    '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/'
    '10_variant_qc/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.variant_qc_summary.tsv'
)

# Depth subgroups analysed (label -> summary column). v6 column names differ from v5's VMISS_15X.
DEPTH_COLS = {'15X': '15X_VMISS', '30X': '30X_VMISS'}
# Other VMISS columns available in the summary (for future overall / case-ctrl / MAF modes).
OTHER_COLS = {'overall': 'VMISS', 'case': 'CASE_VMISS', 'ctrl': 'CTRL_VMISS', 'maf': 'MAF'}

# QC-passed sample list + metrics table (07_sample_qc), used to size the depth subgroups so each
# subgroup's CDF bin can be derived from its own N (see subgroup_sizes / _auto_bin below).
KEEP_ID_DEFAULT = Path(
    '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/'
    '07_sample_qc/run_qc/cteph_agp3k_v6_wgs_merged.sample_qc.keep.id'
)
METRICS_DEFAULT = Path(
    '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/'
    '07_sample_qc/metrics/cteph_agp3k_v6_wgs_merged.sample_qc_metrics.tsv'
)

# Kneedle / CDF defaults.
CDF_STEP = 0.001
HIST_STEP = 0.005            # display histogram bin width (fine enough to show the elbow cleanly)
# Candidate CDF bin widths (all divide the 0.01 threshold grid, so a rounded knee is a bin edge and
# its marker lands on the CDF). _auto_bin picks the smallest that is >= 2/N for a subgroup.
CDF_BIN_CANDIDATES = [0.001, 0.0025, 0.005, 0.01]
# Fallback per-subgroup bins if subgroup sizes can't be read (auto-derivation is preferred).
CDF_STEP_BY_GROUP = {'15X': 0.001, '30X': 0.005}
KNEE_S = 1.0
# weight_x biases the knee toward smaller VMISS (stricter). Set to 8.0 for a conservative elbow on
# v6 data; the knee is stable for weight_x in [7,10]. Lower it to retain more variants.
KNEE_WEIGHT_X = 8.0
KNEE_WEIGHT_Y = 1.0
THRESHOLD_DECIMALS = 2
MIN_CUM_PCT = 5.0            # drop the steep <5% head before knee-finding (avoids false knees)

# Plot palette (tuning.pc house style)
GROUP_COLORS = {'15X': '#2C5A8A', '30X': '#B36A2E', 'overall': '#1B7F79',
                'case': '#C8102E', 'ctrl': '#005A9C', 'maf': '#6d6a61'}
HIST_COLOR = '#C7CDD4'
KNEE_COLOR = '#C44E52'


# ── Per-subgroup CDF bin derivation (from subgroup size N) ────────────────────
def _auto_bin(n):
    """Smallest threshold-aligned CDF bin spanning >= 2 achievable VMISS values (>= 2/N)."""
    if not n or n <= 0:
        return CDF_BIN_CANDIDATES[0]
    need = 2.0 / n
    for b in CDF_BIN_CANDIDATES:
        if b >= need:
            return b
    return CDF_BIN_CANDIDATES[-1]


def subgroup_sizes(keep_id=KEEP_ID_DEFAULT, metrics=METRICS_DEFAULT, cols=None):
    """QC-passed sample count per depth subgroup (label -> N), from keep.id ∩ metrics Target_DP.
    Returns {label: N} or None if the files can't be read / a subgroup is empty."""
    cols = dict(DEPTH_COLS) if cols is None else dict(cols)
    try:
        keep = pd.read_csv(keep_id, sep=r'\s+', header=None, engine='python')
        col = 1 if keep.shape[1] >= 2 else 0
        keep_iids = set(keep.iloc[:, col].astype(str))
        m = pd.read_csv(metrics, sep='\t')
        idc = 'IID' if 'IID' in m.columns else m.columns[1]
        dpc = 'Target_DP' if 'Target_DP' in m.columns else 'Target_Depth'
        norm = m.loc[m[idc].astype(str).isin(keep_iids), dpc].astype(str).str.strip().str.lower()
        out = {}
        for lab in cols:
            key = lab.lower().rstrip('x')                       # '15X' -> '15'
            out[lab] = int(norm.isin([lab.lower(), key, key + 'x']).sum())
        return out if all(v > 0 for v in out.values()) else None
    except Exception:  # noqa: BLE001
        return None


# ── Streaming histogram accumulation (memory-light: O(bins), not O(37M)) ──────
def stream_vmiss_histograms(summary_path=SUMMARY_DEFAULT, cols=None, steps=None,
                            cdf_step=CDF_STEP, chunksize=2_000_000):
    """Single streaming pass; per column accumulate a [0,1] histogram at that column's own CDF bin
    (steps[label], default cdf_step) plus running total/sum/min/max. Returns {label: {...}}."""
    cols = dict(DEPTH_COLS) if cols is None else dict(cols)
    steps = {} if steps is None else dict(steps)

    acc = {}
    for lab, c in cols.items():
        step = steps.get(lab, cdf_step)
        edges = np.round(np.arange(0.0, 1.0 + step, step), 6)
        acc[lab] = {'col': c, 'cdf_step': step, 'edges': edges,
                    'hist': np.zeros(len(edges) - 1, dtype=np.int64),
                    'total': 0, 'sum': 0.0, 'min': np.inf, 'max': -np.inf}

    usecols = list(dict.fromkeys(cols.values()))
    for chunk in pd.read_csv(summary_path, sep='\t', usecols=usecols, chunksize=chunksize):
        for lab, c in cols.items():
            v = pd.to_numeric(chunk[c], errors='coerce').to_numpy(dtype='float64')
            v = v[np.isfinite(v)]
            if v.size == 0:
                continue
            np.clip(v, 0.0, 1.0, out=v)
            acc[lab]['hist'] += np.histogram(v, bins=acc[lab]['edges'])[0]
            acc[lab]['total'] += int(v.size)
            acc[lab]['sum'] += float(v.sum())
            acc[lab]['min'] = min(acc[lab]['min'], float(v.min()))
            acc[lab]['max'] = max(acc[lab]['max'], float(v.max()))

    for lab in acc:
        tot = acc[lab]['total']
        acc[lab]['mean'] = (acc[lab]['sum'] / tot) if tot else float('nan')
        acc[lab]['median'] = _median_from_hist(acc[lab]['hist'], acc[lab]['edges']) if tot else float('nan')
        if tot == 0:
            acc[lab]['min'] = acc[lab]['max'] = float('nan')
    return acc


def _median_from_hist(hist, edges):
    cum = np.cumsum(hist)
    total = cum[-1]
    if total == 0:
        return float('nan')
    idx = int(np.searchsorted(cum, total / 2.0))
    idx = min(idx, len(edges) - 2)
    return float(edges[idx])


# ── Cumulative distribution + Kneedle knee ────────────────────────────────────
def cumulative(hist, edges):
    """Return (x=bin left edges, cum_counts, cum_pct) — proportion of variants <= x."""
    cum_counts = np.cumsum(hist)
    total = cum_counts[-1] if len(cum_counts) else 0
    x = edges[:-1]
    cum_pct = (cum_counts / total * 100.0) if total else np.zeros_like(cum_counts, dtype=float)
    return x, cum_counts, cum_pct


def _kneedle_knee(x, y, weight_x, weight_y):
    """Weighted Kneedle: knee = x maximizing (weight_y*y_norm - weight_x*x_norm) on a
    concave, increasing curve (x, y normalized to [0,1]). weight_x>1 biases toward a
    stricter/smaller knee (v5 used weight_x=4; lower it to retain more). Returns knee x or None."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xr = x.max() - x.min()
    yr = y.max() - y.min()
    if xr <= 0 or yr <= 0:
        return None
    xn = (x - x.min()) / xr
    yn = (y - y.min()) / yr
    d = weight_y * yn - weight_x * xn
    return float(x[int(np.argmax(d))])


def find_knee(x, cum_pct, hist, cdf_step, total,
              knee_s=KNEE_S, weight_x=KNEE_WEIGHT_X, weight_y=KNEE_WEIGHT_Y,
              decimals=THRESHOLD_DECIMALS, min_cum_pct=MIN_CUM_PCT):
    """Weighted-Kneedle knee on the >= min_cum_pct portion of the CDF. Returns a dict.
    (knee_s is kept for API compatibility but unused by this numpy implementation.)"""
    out = {'knee_found': False}
    if len(x) < 3:
        out['message'] = 'insufficient data points'
        return out

    thr = min_cum_pct
    mask = cum_pct >= thr
    if mask.sum() < 3:
        thr = 1.0
        mask = cum_pct >= thr
    if mask.sum() < 3:
        out['message'] = f'insufficient points after >= {thr}% filter'
        return out

    xf, cf = x[mask], cum_pct[mask]
    knee_raw = _kneedle_knee(xf, cf, weight_x, weight_y)
    if knee_raw is None:
        out['message'] = 'no knee detected'
        return out

    knee_vmiss = round(knee_raw, decimals)
    idx = int(round(knee_vmiss / cdf_step))          # bins [0, knee) -> retained (VMISS < knee)
    idx = max(0, min(idx, len(hist)))
    retained_n = int(hist[:idx].sum())
    retained_pct = (retained_n / total * 100.0) if total else 0.0
    out.update({'knee_found': True, 'knee_vmiss': knee_vmiss, 'knee_raw': knee_raw,
                'retained_n': retained_n, 'retained_pct': retained_pct, 'min_cum_pct': thr})
    return out


# ── Joint (dp-mode) retention ─────────────────────────────────────────────────
def compute_joint_retained(summary_path, thresholds, cols=None, chunksize=2_000_000):
    """Count variants passing ALL depth thresholds simultaneously (dp-mode PASS = AND of col<=thr).
    This equals the pipeline's PASS_VMISS set. thresholds: {label: thr}; cols: {label: column}.
    NA in any depth column fails (numpy nan<=thr is False). Returns (joint_n, total, joint_pct)
    with total = all variants."""
    cols = dict(DEPTH_COLS) if cols is None else dict(cols)
    labels = list(thresholds)
    usecols = [cols[l] for l in labels]
    joint = total = 0
    for chunk in pd.read_csv(summary_path, sep='\t', usecols=usecols, chunksize=chunksize):
        total += len(chunk)
        mask = np.ones(len(chunk), dtype=bool)
        for l in labels:
            v = pd.to_numeric(chunk[cols[l]], errors='coerce').to_numpy(dtype='float64')
            mask &= (v <= thresholds[l])
        joint += int(mask.sum())
    pct = (joint / total * 100.0) if total else 0.0
    return joint, total, pct


# ── Orchestration ─────────────────────────────────────────────────────────────
def analyze_depth(summary_path=SUMMARY_DEFAULT, cols=None, steps=None, cdf_step=CDF_STEP,
                  knee_s=KNEE_S, weight_x=KNEE_WEIGHT_X, weight_y=KNEE_WEIGHT_Y,
                  decimals=THRESHOLD_DECIMALS, min_cum_pct=MIN_CUM_PCT, chunksize=2_000_000,
                  verbose=True):
    """Stream the summary once, build per-depth CDF + Kneedle knee, and a recommendations table.
    `steps` sets each subgroup's CDF bin (default CDF_STEP_BY_GROUP) so 30X (few samples) is not sawtoothed."""
    cols = dict(DEPTH_COLS) if cols is None else dict(cols)
    if steps is None:
        sizes = subgroup_sizes(cols=cols)
        if sizes:
            steps = {lab: _auto_bin(sizes[lab]) for lab in cols}
            if verbose:
                print(f'[vmiss] subgroup N: {sizes}  -> CDF bins: {steps}')
        else:
            steps = {lab: CDF_STEP_BY_GROUP.get(lab, CDF_STEP) for lab in cols}
            if verbose:
                print(f'[vmiss] subgroup sizes unavailable; using fallback CDF bins: {steps}')
    else:
        steps = dict(steps)
    if verbose:
        print(f'[vmiss] reading {Path(summary_path).name} (streaming, chunksize={chunksize:,}) ...')
    acc = stream_vmiss_histograms(summary_path, cols=cols, steps=steps, cdf_step=cdf_step, chunksize=chunksize)

    results, rows = {}, []
    for lab in cols:
        a = acc[lab]
        x, cum_counts, cum_pct = cumulative(a['hist'], a['edges'])
        knee = find_knee(x, cum_pct, a['hist'], a['cdf_step'], a['total'],
                         knee_s=knee_s, weight_x=weight_x, weight_y=weight_y,
                         decimals=decimals, min_cum_pct=min_cum_pct)
        res = {'label': lab, **a, 'x': x, 'cum_counts': cum_counts, 'cum_pct': cum_pct, 'knee': knee}
        results[lab] = res
        rows.append({
            'group': lab,
            'vmiss_col': a['col'],
            'total_variants': a['total'],
            'median_vmiss': round(a['median'], 4),
            'mean_vmiss': round(a['mean'], 4),
            'recommended_vmiss': knee.get('knee_vmiss', float('nan')),
            'retained_pct': round(knee.get('retained_pct', float('nan')), 2),
            'retained_n': knee.get('retained_n', -1),
        })
        if verbose:
            k = knee.get('knee_vmiss', 'NA')
            print(f'[vmiss] {lab}: n={a["total"]:,}  median={a["median"]:.4f}  '
                  f'knee VMISS={k} (retained {knee.get("retained_pct", float("nan")):.1f}%)')

    # Joint retention (dp-mode PASS): variants passing ALL depth thresholds at once — equals the
    # pipeline's PASS_VMISS / vmiss_pass_variants set. Second streaming pass with the found knees.
    depth_knees = {lab: results[lab]['knee']['knee_vmiss'] for lab in cols
                   if results[lab].get('knee', {}).get('knee_found')}
    if len(depth_knees) >= 2:
        jn, jt, jp = compute_joint_retained(summary_path, depth_knees, cols=cols, chunksize=chunksize)
        combo = ' & '.join(f'{lab}<={depth_knees[lab]:g}' for lab in depth_knees)
        results['joint'] = {'thresholds': depth_knees, 'combo': combo,
                            'retained_n': jn, 'total': jt, 'retained_pct': jp}
        rows.append({'group': 'joint (dp)', 'vmiss_col': combo, 'total_variants': jt,
                     'median_vmiss': np.nan, 'mean_vmiss': np.nan, 'recommended_vmiss': np.nan,
                     'retained_pct': round(jp, 2), 'retained_n': jn})
        if verbose:
            print(f'[vmiss] joint (dp): {combo}  retained {jn:,} ({jp:.1f}% of {jt:,})')

    results['recommendations_df'] = pd.DataFrame(rows)
    return results


# ── Plotting (tuning.pc house style) ──────────────────────────────────────────
def _apply_rcparams(fs):
    plt.style.use('default')
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': fs['tick'],
        'axes.labelsize': fs['label'], 'xtick.labelsize': fs['tick'], 'ytick.labelsize': fs['tick'],
        'legend.fontsize': fs['legend'], 'axes.linewidth': 1.2,
        'xtick.major.width': 1.15, 'ytick.major.width': 1.15,
    })


def _coarsen(hist, cdf_step, hist_step):
    """Sum the fine histogram into coarser display bins; return (left_edges, counts, width)."""
    factor = max(1, int(round(hist_step / cdf_step)))
    n = (len(hist) // factor) * factor
    coarse = hist[:n].reshape(-1, factor).sum(axis=1)
    left = np.arange(len(coarse)) * (cdf_step * factor)
    return left, coarse, cdf_step * factor


def plot_vmiss_group(result, hist_step=HIST_STEP, xmax=None, show=True):
    """Publication/PPT-grade dual-axis figure for one depth group: histogram bars (left) +
    cumulative % line (right), with the Kneedle knee marked and the retained region shaded.
    Returns the Figure."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    lab = result['label']
    color = GROUP_COLORS.get(lab, '#1B7F79')
    fs = {'title': 22, 'label': 17, 'tick': 14.5, 'legend': 13.5, 'anno': 14}
    _apply_rcparams(fs)

    knee = result['knee']
    # Draw the cumulative on RIGHT bin edges so the curve at x means P(VMISS <= x) exactly; this makes
    # the knee marker (placed at the true retained fraction) sit on the line for any bin width.
    edges, hist = result['edges'], result['hist']
    cum = np.cumsum(hist)
    total_c = cum[-1] if len(cum) else 0
    x = edges[1:]
    cum_pct = (cum / total_c * 100.0) if total_c else np.zeros_like(cum, dtype=float)
    left, counts, width = _coarsen(result['hist'], result['cdf_step'], hist_step)
    kv = knee.get('knee_vmiss') if knee.get('knee_found') else None

    # x-view: focus on the informative low-VMISS region (up to ~99% cumulative), tight for slides.
    if xmax is None:
        hi_idx = int(np.searchsorted(cum_pct, 99.0))
        xv = float(x[min(hi_idx, len(x) - 1)]) if len(x) else 0.15
        xmax = min(0.20, max(0.10, xv, (kv or 0.05) * 4))

    fig, ax1 = plt.subplots(figsize=(11.0, 6.2), facecolor='white')
    ax1.set_facecolor('white')

    # Left axis: variant-count histogram (soft bars).
    m = left <= xmax
    ax1.bar(left[m], counts[m], width=width, align='edge', color=HIST_COLOR,
            edgecolor='white', linewidth=0.3, zorder=2)
    ax1.set_xlabel('Per-variant missing rate (VMISS)', labelpad=8)
    ax1.set_ylabel('Variant count', color='#4A5560', labelpad=8)
    ax1.tick_params(axis='y', colors='#4A5560')
    ax1.set_xlim(0, xmax)
    ax1.set_ylim(bottom=0)
    ax1.grid(axis='y', alpha=0.20, linestyle=':', linewidth=0.7, color='#B8B8B8', zorder=0)
    ax1.set_axisbelow(True)
    for sp in ('top',):
        ax1.spines[sp].set_visible(False)

    # Right axis: cumulative distribution + retained shading + knee.
    ax2 = ax1.twinx()
    if kv is not None:
        ax2.axvspan(0, kv, color=color, alpha=0.07, zorder=1)      # shade "retained" (VMISS <= knee)
    ax2.plot(x, cum_pct, color=color, linewidth=3.0, zorder=4, solid_capstyle='round')
    ax2.set_ylabel('Cumulative % of variants (VMISS ≤ x)', color=color, labelpad=10)
    ax2.tick_params(axis='y', colors=color)
    ax2.set_ylim(0, 101)
    ax2.spines['top'].set_visible(False)

    if kv is not None:
        ax2.axvline(kv, color=KNEE_COLOR, linestyle='--', linewidth=2.0, alpha=0.95, zorder=5)
        ax2.plot([kv], [knee['retained_pct']], marker='D', markersize=11,
                 markerfacecolor=KNEE_COLOR, markeredgecolor='white', markeredgewidth=1.4, zorder=6)
        ax2.annotate(f'threshold  VMISS = {kv:g}\n{knee["retained_pct"]:.1f}% retained '
                     f'({knee["retained_n"]:,})',
                     xy=(kv, knee['retained_pct']), xycoords='data',
                     xytext=(16, -46), textcoords='offset points',
                     fontsize=fs['anno'], color='#2B2B2B', ha='left', va='top',
                     bbox=dict(boxstyle='round,pad=0.45', facecolor='white',
                               edgecolor=KNEE_COLOR, linewidth=1.1, alpha=0.96),
                     zorder=7)

    # Combined legend (bars + cumulative line + threshold), clean and slide-friendly.
    handles = [
        Patch(facecolor=HIST_COLOR, edgecolor='white', label='Variant count'),
        Line2D([0], [0], color=color, linewidth=3.0, label='Cumulative %'),
        Line2D([0], [0], color=KNEE_COLOR, linewidth=2.0, linestyle='--', label=f'Threshold (knee) = {kv:g}' if kv is not None else 'Threshold'),
    ]
    ax2.legend(handles=handles, loc='center right', frameon=True, framealpha=0.96,
               edgecolor='#D8D8D8', borderpad=0.6, handlelength=1.9, fontsize=fs['legend'])

    ax1.set_title(f'{lab} depth subgroup  ·  VMISS call-rate threshold  '
                  f'(n = {result["total"]:,} variants)', fontsize=fs['title'], fontweight='bold', pad=12)
    fig.tight_layout()
    if show:
        plt.show()
    return fig


def plot_depth(results, hist_step=HIST_STEP, show=True):
    """Plot every depth-group result (skips summary keys like recommendations_df / joint)."""
    skip = {'recommendations_df', 'joint'}
    figs = []
    for lab in [k for k in results if k not in skip]:
        figs.append(plot_vmiss_group(results[lab], hist_step=hist_step, show=show))
    return figs
