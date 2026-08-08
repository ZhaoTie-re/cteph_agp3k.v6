#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Flag depth-differential variants for the rare-variant tuning
#           pre-check. Stream the per-variant depth-stratified missing rates
#           (VMISS_15X / VMISS_30X), build the |VMISS_30X - VMISS_15X| CDF,
#           locate the Kneedle knee, and exclude variants whose absolute
#           difference exceeds that knee. Writes a PLINK2 --exclude list
#           (VariantID, no header) plus a 3-panel diagnostic figure.
# Project : cteph_agp3k.v6 / tuning.rv  (rare-variant depth-confounding QC)
# Used by : tuning.rv.nf  process DETECT_DEPTHDIFF_VARIANTS
# ---------------------------------------------------------------------------
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, MaxNLocator
from kneed import KneeLocator   # modified build in env cteph_geno_pro (weight_x / weight_y)

import plot_style as ps

NEEDED = ['VMISS_15X', 'VMISS_30X']
ID_COL = 'VariantID'
CHUNK = 500_000


def parse_args():
    p = argparse.ArgumentParser(
        description='Flag depth-differential variants (Kneedle knee of |VMISS_30X - VMISS_15X| CDF).')
    p.add_argument('--variant-metrics', required=True,
                   help='variant_metrics.txt[.gz] from run_calc_metrics.sh (minAC 0, raw).')
    p.add_argument('--knee-s', type=float, default=1.0, help='Kneedle sensitivity S (default: 1.0).')
    p.add_argument('--knee-weight-x', type=float, default=1.0,
                   help='Kneedle weight_x (default: 1.0).')
    p.add_argument('--knee-weight-y', type=float, default=2.0,
                   help='Kneedle weight_y (default: 2.0).')
    p.add_argument('--out-exclude', required=True,
                   help='Output PLINK2 --exclude list (VariantID, no header).')
    p.add_argument('--out-png', required=True,
                   help='Output 3-panel diagnostic figure (PNG, 600 dpi).')
    p.add_argument('--out-stats', default=None,
                   help='Optional machine-readable sidecar TSV (knee, counts, retained_pct).')
    p.add_argument('--protect-file', default=None,
                   help='Optional file of VariantIDs (one per line) that must NEVER be excluded '
                        '(e.g. candidate signals) — they are dropped from the --exclude list even '
                        'if their |ΔVMISS| exceeds the knee.')
    return p.parse_args()


def _load_protect(path):
    """Read a protect list (one VariantID per line); return a set (empty if unusable)."""
    if not path:
        return set()
    try:
        with open(path) as fh:
            return {ln.strip() for ln in fh if ln.strip()}
    except Exception as e:
        print(f'[warn] could not read protect-file {path} ({e}); nothing protected.')
        return set()


def main():
    args = parse_args()
    KNEE_S, KNEE_WEIGHT_X, KNEE_WEIGHT_Y = args.knee_s, args.knee_weight_x, args.knee_weight_y
    _vsrc = args.variant_metrics

    # Confirm the columns exist before streaming the whole file.
    _head = pd.read_csv(_vsrc, sep='\t', nrows=0)
    _missing = [c for c in NEEDED if c not in _head.columns]
    if _missing:
        raise ValueError(f'Columns {_missing} not found in {Path(_vsrc).name}.')

    # Fixed bin edges (bounded ranges -> safe for streaming). Width 0.001.
    edges_signed = np.linspace(-1.0, 1.0, 2001)   # signed diff in [-1, 1]
    edges_abs = np.linspace(0.0, 1.0, 1001)        # |diff| in [0, 1]
    hist_signed = np.zeros(len(edges_signed) - 1, dtype=np.int64)
    hist_abs = np.zeros(len(edges_abs) - 1, dtype=np.int64)

    n_total = n_valid = 0
    sum_signed = sum_abs = 0.0
    min_signed, max_signed, max_abs = np.inf, -np.inf, 0.0

    for chunk in pd.read_csv(
        _vsrc, sep='\t', usecols=NEEDED, chunksize=CHUNK, compression='infer',
        dtype={'VMISS_15X': 'float32', 'VMISS_30X': 'float32'},
    ):
        signed = (chunk['VMISS_30X'] - chunk['VMISS_15X']).to_numpy(dtype=np.float64)
        n_total += signed.size
        signed = signed[np.isfinite(signed)]      # drop variants where a stratum had no samples
        n_valid += signed.size
        if not signed.size:
            continue
        absd = np.abs(signed)
        sum_signed += float(signed.sum())
        sum_abs += float(absd.sum())
        min_signed = min(min_signed, float(signed.min()))
        max_signed = max(max_signed, float(signed.max()))
        max_abs = max(max_abs, float(absd.max()))
        hist_signed += np.histogram(signed, bins=edges_signed)[0]
        hist_abs += np.histogram(absd, bins=edges_abs)[0]

    if n_valid == 0:
        raise ValueError('No variants have finite VMISS_15X and VMISS_30X simultaneously.')

    mean_signed = sum_signed / n_valid
    mean_abs = sum_abs / n_valid

    # CDF of |diff| from the accumulated histogram.
    cdf_abs = np.cumsum(hist_abs) / n_valid            # P(|diff| <= right edge)
    cdf_x = edges_abs[1:]
    quantiles = {q: float(cdf_x[np.searchsorted(cdf_abs, q, side='left')]) for q in (0.95, 0.99)}

    # X-range focused on the populated region (rounded to 0.01, with headroom).
    abs_hi = min(1.0, max(0.01, np.ceil(max_abs / 0.01) * 0.01))
    sign_hi = min(1.0, max(0.01, np.ceil(max(abs(min_signed), abs(max_signed)) / 0.01) * 0.01))

    # ── Kneedle knee of the |diff| CDF (concave, increasing) ──────────────────
    # Restrict to the populated range so the long flat tail (CDF==1) does not
    # dominate the [0,1] normalisation and bias the knee.
    _km = cdf_x <= abs_hi
    knee_thr = knee_cdf = None
    if int(_km.sum()) >= 2:
        try:
            knee_locator = KneeLocator(
                cdf_x[_km], cdf_abs[_km], S=KNEE_S, curve='concave', direction='increasing',
                weight_x=KNEE_WEIGHT_X, weight_y=KNEE_WEIGHT_Y,
            )
            knee_thr = knee_locator.knee                  # threshold on |diff|
            knee_cdf = knee_locator.knee_y                # CDF value at the knee
        except Exception as _e:
            print(f'[warn] KneeLocator failed ({_e}); no knee will be drawn.')

    if knee_thr is None:
        raise RuntimeError('Kneedle did not return a knee; cannot choose a threshold.')

    # Round the knee to 2 decimals (spec) and recompute the retained fraction at
    # the rounded threshold, so the figure/annotation match the actual cut used.
    knee_thr = round(float(knee_thr), 2)
    _ki = int(np.searchsorted(cdf_x, knee_thr, side='right')) - 1
    _ki = max(0, min(_ki, len(cdf_abs) - 1))
    knee_cdf = float(cdf_abs[_ki])

    retained_pct = 100.0 * float(knee_cdf)        # |diff| <= knee  -> kept
    removed_pct = 100.0 - retained_pct            # |diff| >  knee  -> removed

    print(f'variants total            : {n_total:,}')
    print(f'variants with valid diff  : {n_valid:,}  (NaN dropped: {n_total - n_valid:,})')
    print(f'signed diff  mean         : {mean_signed:.4g}   range [{min_signed:.4g}, {max_signed:.4g}]')
    print(f'abs    diff  mean / max   : {mean_abs:.4g} / {max_abs:.4g}')
    print(f'abs    diff  95th / 99th  : {quantiles[0.95]:.4g} / {quantiles[0.99]:.4g}')
    print(f'kneedle (S={KNEE_S}, wx={KNEE_WEIGHT_X}, wy={KNEE_WEIGHT_Y}) -> '
          f'knee |diff| = {knee_thr}  |  retained {retained_pct:.2f}% / removed {removed_pct:.2f}%')

    # ── Publication-grade figure (3 panels; shared style) ─────────────────────
    ps.setup_style('slide')

    BAR_COLOR, LINE_COLOR, KNEE_COLOR = ps.BAR_COLOR, ps.LINE_COLOR, ps.KNEE_COLOR
    cnt_fmt = ps.COUNT_FMT

    cnt_max = max(hist_signed.max(), hist_abs.max()) * 1.06

    ctr_s = 0.5 * (edges_signed[:-1] + edges_signed[1:])
    ctr_a = 0.5 * (edges_abs[:-1] + edges_abs[1:])

    # a and b are the same count axis (signed vs absolute difference), so they share
    # it; c is a cumulative fraction and keeps its own.
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 6.4))
    ax_s, ax_a, ax_c = axes
    ax_a.sharey(ax_s)

    def _style(ax, grid_axis='y'):
        ps.despine(ax, grid_axis=grid_axis)

    # (a) signed difference histogram
    ax_s.bar(ctr_s, hist_signed, width=np.diff(edges_signed), align='center',
             color=BAR_COLOR, edgecolor='white', linewidth=0.12, zorder=3)
    ax_s.axvline(0.0, color='#888888', linewidth=1.0, zorder=2)
    ax_s.set_xlim(-sign_hi, sign_hi)
    ax_s.set_ylim(0, cnt_max)
    ax_s.xaxis.set_major_locator(MaxNLocator(nbins=8, steps=[1, 2, 2.5, 5, 10]))
    ax_s.yaxis.set_major_formatter(cnt_fmt)
    ax_s.set_xlabel(ps.VMISS_SIGNED, labelpad=6)
    ax_s.set_ylabel('Number of variants', labelpad=6)
    ps.panel_tag(ax_s, 'a', 'Signed difference')
    _style(ax_s)

    # (b) absolute difference histogram (shared y-range with a)
    ax_a.bar(ctr_a, hist_abs, width=np.diff(edges_abs), align='center',
             color=BAR_COLOR, edgecolor='white', linewidth=0.12, zorder=3)
    ax_a.set_xlim(0, abs_hi)
    ax_a.set_ylim(0, cnt_max)
    ax_a.xaxis.set_major_locator(MaxNLocator(nbins=8, steps=[1, 2, 2.5, 5, 10]))
    ax_a.yaxis.set_major_formatter(cnt_fmt)
    # Echo the exclusion knee onto the |Δ| distribution so the cut is visible
    # against the histogram it acts on (variants to its right are excluded).
    ax_a.axvline(knee_thr, color=KNEE_COLOR, linestyle=(0, (4, 3)), linewidth=1.4, zorder=4)
    # Dropped clear of the title slot; at cnt_max*0.98 it sat against the panel title.
    ax_a.text(knee_thr, cnt_max * 0.88, f' knee = {knee_thr:.2f}', color=KNEE_COLOR,
              fontsize=10, fontweight='bold', ha='left', va='top')
    ax_a.set_xlabel(ps.VMISS_ABS, labelpad=6)
    ps.panel_tag(ax_a, 'b', 'Absolute difference')
    _style(ax_a)

    # (c) absolute difference CDF + kneedle knee
    ax_c.plot(cdf_x, cdf_abs, color=LINE_COLOR, linewidth=2.0, zorder=4)
    ax_c.fill_between(cdf_x, cdf_abs, color=LINE_COLOR, alpha=0.10, zorder=1)
    # Thin dotted guides drop from the knee to each axis (read the threshold on x,
    # the retained fraction on y); a single adjacent label — no crossing arrow.
    ax_c.vlines(knee_thr, 0, knee_cdf, color=KNEE_COLOR, linestyle=(0, (2, 2)), linewidth=1.2, zorder=5)
    ax_c.hlines(knee_cdf, 0, knee_thr, color=KNEE_COLOR, linestyle=(0, (2, 2)), linewidth=1.2, zorder=5)
    ax_c.plot([knee_thr], [knee_cdf], 'o', color=KNEE_COLOR, markersize=8,
              markeredgecolor='white', markeredgewidth=1.0, zorder=6)
    # Anchor the label just BELOW-right of the marker with a clamped top, so it stays
    # clear of the frame and the horizontal guide even when the knee sits near CDF=1.
    ax_c.annotate(f'knee = {knee_thr:.2f}\nkeeps {retained_pct:.1f}%',
                  xy=(knee_thr + 0.03 * abs_hi, min(knee_cdf - 0.03, 0.86)), xycoords='data',
                  ha='left', va='top', fontsize=10.5, color=KNEE_COLOR, fontweight='bold',
                  bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                            edgecolor=KNEE_COLOR, linewidth=0.8, alpha=0.92))
    ax_c.set_xlim(0, abs_hi)
    ax_c.set_ylim(0, 1.005)
    ax_c.xaxis.set_major_locator(MaxNLocator(nbins=8, steps=[1, 2, 2.5, 5, 10]))
    ax_c.yaxis.set_major_locator(MultipleLocator(0.2))
    ax_c.set_xlabel(ps.VMISS_ABS, labelpad=6)
    ax_c.set_ylabel('Cumulative fraction of variants', labelpad=6)
    ps.panel_tag(ax_c, 'c', 'Cumulative difference (CDF)')
    _style(ax_c)

    ps.share_y([ax_s, ax_a], 'Number of variants')
    fig.suptitle('Variant QC — per-variant depth-stratified missing-rate difference',
                 fontsize=16, y=0.965)
    ps.caption_block(
        fig, top=0.885, wspace=0.22,
        title=('Variant QC excludes only the heavy right tail of depth-differential missingness, set by the '
               f'Kneedle knee at Δ = {knee_thr:.2f} ({removed_pct:.1f}% of variants removed).'),
        panels=[
            r'histogram of the signed per-variant missing-rate difference, $V^{30\times}_v - V^{15\times}_v$.',
            'histogram of its absolute value Δ; the dashed line marks the exclusion knee.',
            (r'cumulative distribution $F(\Delta)$ with the Kneedle knee (marker); variants with Δ above it are '
             f'excluded — {retained_pct:.1f}% kept, {removed_pct:.1f}% removed.'),
        ],
        interpret=('the Δ distribution concentrates near zero with a heavy right tail; the data-driven knee isolates '
                   'that tail, dropping the variants whose missingness is depth-driven while retaining the bulk. '
                   'This suppresses depth-differential dropout before the burden test; protected candidate variants '
                   'are exempt from exclusion.'),
        notes=(f'per-variant depth-stratified missing rates recomputed on the sample-cleaned set '
               f'({Path(args.variant_metrics).name}); {n_valid:,} variants with both strata observed. '
               'a and b share one count axis; the dashed line in b and the marker in c are the same knee.'),
        defs=['delta_v', 'knee'],
        model=ps.FORMULAS['vmiss_diff'])
    fig.savefig(args.out_png)
    print(f'Figure         : {args.out_png}')

    # ── Export VariantIDs to REMOVE at the knee (streamed) ────────────────────
    _hcols = pd.read_csv(_vsrc, sep='\t', nrows=0).columns
    if ID_COL not in _hcols:
        raise ValueError(f"'{ID_COL}' not found in {Path(_vsrc).name}; columns: {list(_hcols)}")

    protect = _load_protect(args.protect_file)              # never exclude these VariantIDs
    n_scanned = n_finite = n_removed = 0
    n_prot_seen = n_prot_rescued = 0
    with open(args.out_exclude, 'w') as fout:
        for chunk in pd.read_csv(
            _vsrc, sep='\t', usecols=[ID_COL] + NEEDED, chunksize=CHUNK, compression='infer',
            dtype={'VMISS_15X': 'float32', 'VMISS_30X': 'float32'},
        ):
            ids_all = chunk[ID_COL].astype(str).to_numpy()
            diff = (chunk['VMISS_30X'] - chunk['VMISS_15X']).abs().to_numpy()
            finite = np.isfinite(diff)
            mask = finite & (diff > knee_thr)               # large depth-dependent missingness -> drop
            n_scanned += len(chunk)
            n_finite += int(finite.sum())
            if protect:
                is_prot = np.isin(ids_all, list(protect))
                n_prot_seen += int(is_prot.sum())
                n_prot_rescued += int((mask & is_prot).sum())   # would-be-excluded but protected
                mask = mask & ~is_prot                          # keep protected variants
            if mask.any():
                fout.write('\n'.join(ids_all[mask]) + '\n')
                n_removed += int(mask.sum())

    if n_finite == 0:
        raise ValueError('No variants with valid depth-stratified missingness were found.')

    pct_valid = 100.0 * n_removed / n_finite
    pct_all = 100.0 * n_removed / n_scanned if n_scanned else float('nan')
    print(f'threshold (knee)   : |VMISS_30X - VMISS_15X| > {knee_thr:.2f}')
    print(f'variants scanned   : {n_scanned:,}')
    print(f'variants w/ valid  : {n_finite:,}')
    print(f'variants to REMOVE : {n_removed:,}  ({pct_valid:.2f}% of valid, {pct_all:.2f}% of all)')
    if protect:
        print(f'protected variants : {len(protect)} requested; {n_prot_seen} present; '
              f'{n_prot_rescued} rescued from exclusion (kept despite |ΔVMISS| > knee)')
    print(f'Exclude file       : {args.out_exclude}')

    if args.out_stats:
        with open(args.out_stats, 'w') as fh:
            fh.write('metric\tvalue\n')
            fh.write(f'knee\t{knee_thr:.6g}\n')
            fh.write(f'n_total\t{n_scanned}\n')
            fh.write(f'n_valid\t{n_finite}\n')
            fh.write(f'n_removed\t{n_removed}\n')
            fh.write(f'removed_pct\t{pct_valid:.4f}\n')
            fh.write(f'retained_pct\t{100.0 - pct_valid:.4f}\n')
            fh.write(f'n_protected_requested\t{len(protect)}\n')
            fh.write(f'n_protected_present\t{n_prot_seen}\n')
            fh.write(f'n_protected_rescued\t{n_prot_rescued}\n')
        print(f'Stats sidecar      : {args.out_stats}')


if __name__ == '__main__':
    main()
