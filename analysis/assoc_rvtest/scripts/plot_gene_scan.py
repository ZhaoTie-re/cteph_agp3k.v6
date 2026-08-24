#!/usr/bin/env python3
# ==============================================================================
# plot_gene_scan.py — the scan figure for ONE gene-based association scan.
#
# The gene-level counterpart of _shared/scripts/plot_manhattan_qq.py, and it
# follows the same visual grammar (chromosome bands, both tier lines, QQ with a
# Beta band and lambda, caption below the panels).
#
# One thing is deliberately different. A variant-level Manhattan draws two FIXED
# lines, 5e-8 and 1e-6, because those thresholds are properties of the genome,
# not of the scan. A gene-based scan is tiered by Benjamini-Hochberg instead, and
# BH adapts to how many genes were tested — 164 genes in the HIGH stratum against
# ~8.6k in MODERATE+HIGH. The P that FDR 0.05 implies therefore differs between
# strata (measured here: 1.3e-4 vs 4.5e-5), so each scan draws ITS OWN
# BH-implied lines, read from scan_qc.tsv. Drawing a shared constant would say
# the strata were judged by one rule when they were not.
# ==============================================================================
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_shared' / 'scripts'))
import plot_style as S          # noqa: E402
import figure_doc               # noqa: E402

TIER_SIG, TIER_SUG = 'significant', 'suggestive'
PLOT_H = 4.3


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scan', required=True, help='gene_scan.tsv from gene_scan.py')
    p.add_argument('--qc', required=True, help='scan_qc.tsv from gene_scan.py')
    p.add_argument('--cohort', required=True)
    p.add_argument('--stratum', required=True)
    p.add_argument('--method', required=True)
    p.add_argument('--covar-label', default='SEX + PCs')
    p.add_argument('--label-genes', type=int, default=12,
                   help='How many hit genes to name, smallest P first.')
    p.add_argument('--out-png', required=True)
    return p.parse_args()


def chrom_key(c):
    c = str(c).replace('chr', '')
    return {'X': 23, 'Y': 24, 'MT': 25, 'M': 25}.get(c, int(c) if c.isdigit() else 99)


def manhattan(ax, df, qc):
    """Genes at their midpoints, chromosome-banded; returns the x lookup."""
    df = df.assign(_k=df['chrom'].map(chrom_key),
                   _mid=(df['start'] + df['end']) / 2).sort_values(['_k', '_mid'])
    offset, off, ticks, names, spans = 0, {}, [], [], []
    for i, (_k, g) in enumerate(df.groupby('_k', sort=True)):
        c = str(g['chrom'].iloc[0])
        off[c] = offset
        xs = offset + g['_mid'].to_numpy()
        ticks.append((xs[0] + xs[-1]) / 2)
        names.append(c)
        ax.scatter(xs, -np.log10(g['Pvalue'].to_numpy()), s=7.0,
                   c=S.CHROM_BANDS[i % 2], linewidths=0, rasterized=True, zorder=2)
        spans.append((xs[0], xs[-1]))
        offset = xs[-1] + 2.0e7
    x0, x1 = spans[0][0] - 1e7, spans[-1][1] + 1e7

    # The BH-implied lines for THIS scan. Either can be absent — a scan that
    # found nothing has no gene inside the tier and therefore no threshold; that
    # is a real outcome, not an error, so the line is simply not drawn.
    p_sig = _num(qc, 'p_at_significant')
    p_sug = _num(qc, 'p_at_suggestive')
    S.tier_lines(ax, p_sig, p_sug, label=False)
    box = dict(boxstyle='square,pad=0.15', facecolor='white', edgecolor='none', alpha=0.85)
    for p_val, col, lab in ((p_sig, S.ACCENT, 'FDR 0.05'), (p_sug, S.REFERENCE, 'FDR 0.10')):
        if p_val is None or not np.isfinite(p_val):
            continue
        ax.text(x1, -np.log10(p_val), f'{lab}  {S.p_tex(p_val)}',
                color=col, fontsize=plt.rcParams['legend.fontsize'] - 1.5,
                ha='right', va='bottom', zorder=7, bbox=box)

    ax.set_xlim(x0, x1)
    ax.set_ylabel(S.NEGLOG10P)
    top = max(-np.log10(df['Pvalue'].min()) * 1.12, 5.0)
    if np.isfinite(p_sig or np.nan):
        top = max(top, -np.log10(p_sig) * 1.15)
    ax.set_ylim(0, top)
    ax.set_xticks(ticks)
    ax.set_xticklabels(names)
    ax.set_xlabel('chromosome')
    S.thin_tick_labels(ax, 'x')
    S.despine(ax, grid_axis='y')
    return off


def _num(qc, col):
    if col not in qc.columns:
        return None
    v = pd.to_numeric(qc[col], errors='coerce').iloc[0]
    return float(v) if np.isfinite(v) else None


def label_hits(ax, df, off, n_label):
    """Name the strongest hit genes in a strip above the axes."""
    hits = df[df['tier'].isin([TIER_SIG, TIER_SUG])].nsmallest(max(0, n_label), 'Pvalue')
    if not len(hits):
        return 0
    xs = [off.get(str(r.chrom), 0) + (r.start + r.end) / 2 for _, r in hits.iterrows()]
    ys = [-np.log10(r.Pvalue) for _, r in hits.iterrows()]
    is_sig = (hits.tier == TIER_SIG).to_numpy()
    base = plt.rcParams['legend.fontsize']
    S.gene_labels(ax, xs, ys, list(hits['Gene'].astype(str)),
                  color=[S.ACCENT if g else S.REFERENCE for g in is_sig],
                  weight=['bold' if g else 'normal' for g in is_sig],
                  fontsize=[base if g else base - 1.0 for g in is_sig])
    return len(hits)


def qq(ax, p, lam):
    n = len(p)
    obs = -np.log10(np.sort(p))
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    i = np.arange(1, n + 1)
    lo = -np.log10(stats.beta.ppf(0.975, i, n - i + 1))
    hi = -np.log10(stats.beta.ppf(0.025, i, n - i + 1))
    ax.fill_between(exp, lo, hi, color=S.NEUTRAL, linewidth=0, zorder=1)
    m = max(exp.max(), obs.max())
    ax.plot([0, m], [0, m], color=S.ACCENT, lw=1.1, zorder=2)
    ax.scatter(exp, obs, s=6.0, c=S.DATA_DARK, linewidths=0, rasterized=True, zorder=3)
    ax.set_xlabel(f'expected {S.NEGLOG10P}')
    ax.set_ylabel(f'observed {S.NEGLOG10P}')
    lim = m * 1.05
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_box_aspect(1)
    ax.text(0.04, 0.95, f'{S.LAMBDA_GC} = {lam:.3f}', transform=ax.transAxes,
            ha='left', va='top', fontsize=plt.rcParams['legend.fontsize'] + 1,
            fontweight='bold')
    S.despine(ax)


def main():
    args = parse_args()
    S.setup_style()
    df = pd.read_csv(args.scan, sep='\t')
    qc = pd.read_csv(args.qc, sep='\t')
    for c in ('Pvalue', 'start', 'end'):
        df[c] = pd.to_numeric(df[c], errors='coerce')
    df = df[np.isfinite(df['Pvalue']) & np.isfinite(df['start'])].copy()
    if not len(df):
        raise SystemExit(f'{args.scan}: no gene with a usable P and position')

    fig, axes = plt.subplots(1, 2, figsize=(S.COL_DOUBLE, PLOT_H),
                             gridspec_kw={'width_ratios': [2.05, 1.0]})
    off = manhattan(axes[0], df, qc)
    n_lab = label_hits(axes[0], df, off, args.label_genes)
    S.panel_tag(axes[0], 'a')
    lam = float(pd.to_numeric(qc['lambda_gc'], errors='coerce').iloc[0])
    qq(axes[1], df['Pvalue'].to_numpy(), lam)
    S.panel_tag(axes[1], 'b')

    n_sig = int((df['tier'] == TIER_SIG).sum())
    n_sug = int((df['tier'] == TIER_SUG).sum())
    p_sig, p_sug = _num(qc, 'p_at_significant'), _num(qc, 'p_at_suggestive')
    top = df.nsmallest(1, 'Pvalue').iloc[0] if len(df) else None
    stratum_h = args.stratum.replace('_', '+').upper()

    S.caption_block(
        fig, plot_h=PLOT_H,
        title=(f'{args.cohort}, {stratum_h} variants, {args.method.upper()}: '
               f'{n_sig} gene{"" if n_sig == 1 else "s"} at FDR < 0.05'
               + (f', led by {top["Gene"]} ({S.p_tex(top["Pvalue"])})' if n_sig else '')
               + f'; {n_sug} further at FDR < 0.10.'),
        panels=[
            (f'Gene-based association across the genome; {len(df):,} genes tested. Solid line, '
             f'FDR 0.05; dashed, FDR 0.10 — both implied by Benjamini-Hochberg on this scan.'),
            (f'Quantile-quantile plot of the same gene {S.NEGLOG10P}, against a Beta(i, n-i+1) '
             f'95% null band.'),
        ],
        notes=(f'{args.method.upper()} on variants of impact {stratum_h}; {args.covar_label}. '
               f'Tiers are Benjamini-Hochberg over the {len(df):,} genes of THIS scan, so the '
               f'implied P differs between strata; here FDR 0.05 falls at '
               f'{S.p_tex(p_sig) if p_sig else "n/a"}. '
               f'{n_lab} gene(s) named. Full explanation: ../README.md'),
        top_pad=1.02, wspace=0.30, left='auto', right=0.985, margin_axes=axes)
    fig.savefig(args.out_png)
    plt.close(fig)

    figure_doc.write_stats(
        args.out_png, peak=f'{args.stratum}.{args.method}',
        values=[('cohort', args.cohort), ('impact stratum', args.stratum),
                ('test method', args.method),
                ('genes tested', len(df)),
                (S.LAMBDA_GC.replace('$', ''), round(lam, 4)),
                ('smallest P', float(df['Pvalue'].min())),
                ('top gene', str(top['Gene']) if top is not None else None),
                ('P at FDR 0.05', p_sig), ('P at FDR 0.10', p_sug),
                ('genes at FDR < 0.05', n_sig), ('genes at FDR < 0.10', n_sig + n_sug),
                ('genes named', n_lab)])
    print(f'[plot_gene_scan] {args.cohort} {args.stratum}/{args.method}: '
          f'{len(df):,} genes, lambda={lam:.3f}, {n_sig} significant -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in plot_gene_scan: {e}', file=sys.stderr)
        sys.exit(1)
