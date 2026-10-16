#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The cohort_manhattan.png figure of assoc_saige, re-drawn for slides
#           and for a model comparison: one row per cohort, one Manhattan per
#           model, and one QQ per row on which the models are overlaid.
#             --models A,B  -> [A Manhattan | B Manhattan | QQ (A and B)]
#             --models A    -> [A Manhattan | QQ]
#           Kept from cohort_manhattan.png, because they are what makes heights
#           comparable: one genomic axis for every panel, one -log10 P limit for
#           every Manhattan, one QQ limit for every QQ, the same threshold pair
#           (5e-8 red, 1e-5 grey dashed), alternating chromosome bands, the
#           pointwise 95% Beta band, dotted guides at every locus genome-wide in
#           any panel (named once above the top row), and a diamond at each
#           panel's own genome-wide leads. No caption block: the text is in the
#           family README.
# Inputs  : <cohort>.<model>.npz and .stats.json (reduce_scan.py); the upstream
#           model_peaks_annotation.tsv of every cohort
# Outputs : PNG, and a JSON of what the figure shows (for the README)
# ---------------------------------------------------------------------------
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.ticker import MultipleLocator

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

CHROM_GAP = 2.0e7
ROW0_EXTRA = 1.22
GW = 'genome_wide'


def load_scan(cohort, model):
    z = np.load(f'{cohort}.{model}.npz')
    return dict(chrom=z['chrom'], pos=z['pos'], nlp=z['nlp'], qq=z['qq'],
                stats=C.load_json(f'{cohort}.{model}.stats.json'))


def read_peaks(paths, cohorts, upstream):
    """{(model, cohort): genome-wide leads} from the upstream model_peaks_annotation.tsv files."""
    out = {}
    for model, files in paths.items():
        for c, f in zip(cohorts, C.csv_list(files)):
            d = pd.read_csv(f, sep='\t', dtype={'chrom': str, 'variant_id': str})
            d = d[(d.model == upstream[model]) & (d.tier == GW)].copy()
            d['P'] = pd.to_numeric(d.P)
            d['pos'] = pd.to_numeric(d.pos).astype(np.int64)
            out[(model, c)] = d
    return out


def label_for(r):
    for k in ('Gene', 'rsID'):
        v = str(r.get(k, '.'))
        if v not in ('.', '', 'nan'):
            return v
    return str(r['variant_id'])


def genome_axis(scans):
    """One cumulative offset per chromosome over every scan drawn, so a column is the same x in all panels."""
    lo, hi = {}, {}
    for s in scans:
        for ch in np.unique(s['chrom']):
            p = s['pos'][s['chrom'] == ch]
            lo[ch], hi[ch] = min(lo.get(ch, p.min()), p.min()), max(hi.get(ch, p.max()), p.max())
    off, ticks, names, cur = {}, [], [], 0.0
    for ch in sorted(lo):
        off[ch] = cur - lo[ch]
        ticks.append(cur + (hi[ch] - lo[ch]) / 2.0)
        names.append(str(ch))
        cur += hi[ch] - lo[ch] + CHROM_GAP
    return off, ticks, names, (-1e7, cur - CHROM_GAP + 1e7)


def overplot_filter(x, y, xlim, ylim, keep_above=3.0, nx=1200, ny=450):
    """Indices to draw: every point at -log10 P >= keep_above, and below it one point per cell of an
    nx x ny grid over the panel. A cell is ~2 px at 600 dpi against a marker ~15 px across, so the
    dropped points lie entirely under drawn ones and the image is unchanged, while the dense body
    (millions of points at P > 1e-3) shrinks to a few hundred thousand and renders in seconds."""
    hi = y >= keep_above
    gx = np.clip(((x - xlim[0]) / (xlim[1] - xlim[0]) * nx).astype(np.int64), 0, nx)
    gy = np.clip((y / ylim[1] * ny).astype(np.int64), 0, ny)
    _u, first = np.unique(gx[~hi] * (ny + 1) + gy[~hi], return_index=True)
    return np.concatenate([np.flatnonzero(~hi)[first], np.flatnonzero(hi)])


def draw_manhattan(S, ax, s, off, order, guides, gw, ylim, xlim, thresh_labels):
    x = np.zeros(len(s['pos']), dtype=float)
    band = np.zeros(len(s['pos']), dtype=np.int8)
    for ch in np.unique(s['chrom']):
        sel = s['chrom'] == ch
        x[sel] = off[ch] + s['pos'][sel]
        band[sel] = order[ch] % 2
    y = s['nlp'].astype(float)
    k = np.sort(overplot_filter(x, y, xlim, ylim))
    ax.scatter(x[k], y[k], s=3.2, c=band[k], cmap=ListedColormap(S.CHROM_BANDS), vmin=0, vmax=1,
               linewidths=0, rasterized=True, zorder=2)
    for g in guides:
        ax.axvline(g, color=S.NEUTRAL_D, lw=0.7, ls=':', zorder=1)
    ax.axhline(-np.log10(C.GW_ALPHA), color=S.ACCENT, lw=1.1, zorder=3)
    ax.axhline(-np.log10(C.SUGGEST_ALPHA), color=S.REFERENCE, lw=0.9, ls='--', zorder=3)
    if len(gw):
        st = S.TIER_STYLE[GW]
        ax.scatter([off[C.chrom_key(r.chrom)] + int(r.pos) for _, r in gw.iterrows()],
                   -np.log10(gw.P.astype(float)), s=st['size'], marker=st['marker'], c=st['color'],
                   edgecolors='white', linewidths=0.6, zorder=6)
    if thresh_labels:
        box = dict(boxstyle='square,pad=0.15', facecolor='white', edgecolor='none', alpha=0.85)
        for a, col in ((C.GW_ALPHA, S.ACCENT), (C.SUGGEST_ALPHA, S.REFERENCE)):
            ax.text(0.998, -np.log10(a), S.p_tex(a).replace('1.00', '1').replace('5.00', '5'),
                    transform=ax.get_yaxis_transform(), color=col,
                    fontsize=plt.rcParams['legend.fontsize'] - 1, ha='right', va='bottom', zorder=7, bbox=box)
    ax.set_ylim(*ylim)
    ax.yaxis.set_major_locator(MultipleLocator(2))
    S.despine(ax, grid_axis='y')


def draw_qq(S, ax, scans, colors, single):
    n0 = len(scans[0]['qq'])
    rank, exp, _ = C.qq_points(scans[0]['qq'])
    lo, hi = C.qq_band(rank, n0)
    ax.fill_between(exp, lo, hi, color=S.NEUTRAL, linewidth=0, zorder=1)
    top = exp.max()
    for k, (s, col) in enumerate(zip(scans, colors)):
        _r, e, o = C.qq_points(s['qq'])
        ax.scatter(e, o, s=4.0, c=col, linewidths=0, rasterized=True, zorder=3 + k)
        top = max(top, e.max(), o.max())
    ax.set_box_aspect(1)
    for k, (s, col) in enumerate(zip(scans, colors)):
        ax.text(0.04, 0.95 - 0.12 * k, f'{S.LAMBDA_GC} = {s["stats"]["lambda_gc"]:.3f}', transform=ax.transAxes,
                ha='left', va='top', fontsize=plt.rcParams['legend.fontsize'] + 1, fontweight='bold',
                color=S.INK if single else col, zorder=8)
    ax.set_ylabel(f'observed {S.NEGLOG10P}')
    S.despine(ax)
    return top


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cohorts', required=True)
    ap.add_argument('--cohort-labels', required=True)
    ap.add_argument('--cohort-colors', required=True)
    ap.add_argument('--models', required=True, help='one or two of fixed, grm_pc, grm_only')
    ap.add_argument('--peaks', action='append', required=True, metavar='MODEL=FILE,FILE,FILE',
                    help='upstream model_peaks_annotation.tsv per cohort, in cohort order')
    ap.add_argument('--upstream-models', required=True, help='fixed=additive,grm_pc=main,grm_only=detect')
    ap.add_argument('--column-titles', default='', help='MODEL=title overrides, "|"-separated')
    ap.add_argument('--plot-style', required=True)
    ap.add_argument('--width', type=float, default=0.0)
    ap.add_argument('--out-png', required=True)
    ap.add_argument('--out-json', required=True)
    a = ap.parse_args()
    S = C.import_style(a.plot_style)
    psize = C.setup(S)
    cohorts = C.csv_list(a.cohorts)
    coh = C.cohort_table(a.cohorts, a.cohort_labels, a.cohort_colors)
    models = C.csv_list(a.models)
    assert 1 <= len(models) <= 2 and all(m in C.MODELS for m in models), models
    single = len(models) == 1
    upstream = C.kv(a.upstream_models)
    titles = dict(t.split('=', 1) for t in a.column_titles.split('|') if '=' in t)
    peaks = read_peaks(C.kv(a.peaks), cohorts, upstream)

    data = {(m, c): load_scan(c, m) for c in cohorts for m in models}
    off, ticks, names, xlim = genome_axis(data.values())
    order = {ch: i for i, ch in enumerate(sorted(off))}
    best = {}
    for (m, c), d in peaks.items():
        if m not in models:
            continue
        for _, r in d.iterrows():
            lab = label_for(r)
            if lab not in best or r.P < best[lab]['P']:
                best[lab] = dict(label=lab, chrom=C.chrom_key(r.chrom), pos=int(r.pos), P=float(r.P),
                                 model=m, cohort=c)
    loci = sorted(best.values(), key=lambda d: (d['chrom'], d['pos']))
    guides = [off[d['chrom']] + d['pos'] for d in loci]
    top = max(float(s['nlp'].max()) for s in data.values())
    ylim = (0, max(8.4, top * 1.08))

    W = a.width or (10.0 if not single else 7.6)
    H = 6.5
    fig = plt.figure(figsize=(W, H))
    ncol = len(models) + 1
    gs = fig.add_gridspec(len(cohorts), ncol, height_ratios=[ROW0_EXTRA] + [1.0] * (len(cohorts) - 1),
                          width_ratios=[2.45] * len(models) + [1.0])
    fig.subplots_adjust(left=0.62 / W, right=1 - 0.12 / W, top=1 - 0.66 / H, bottom=0.48 / H,
                        wspace=0.22 if not single else 0.30, hspace=0.38)
    man = [[None] * len(models) for _ in cohorts]
    for i in range(len(cohorts)):
        for j in range(len(models)):
            man[i][j] = fig.add_subplot(gs[i, j], sharex=man[0][0] if (i or j) else None)
    qqs = [fig.add_subplot(gs[i, -1]) for i in range(len(cohorts))]

    qtop = 0.0
    for i, c in enumerate(cohorts):
        for j, m in enumerate(models):
            draw_manhattan(S, man[i][j], data[(m, c)], off, order, guides, peaks.get((m, c), pd.DataFrame()),
                           ylim, xlim, thresh_labels=(i == 0))
            if j == 0:
                man[i][j].set_ylabel(S.NEGLOG10P)
            if i < len(cohorts) - 1:
                man[i][j].tick_params(labelbottom=False)
        qtop = max(qtop, draw_qq(S, qqs[i], [data[(m, c)] for m in models],
                                 [C.MODELS[m]['color'] for m in models], single))
    man[0][0].set_xlim(*xlim)
    qlim = qtop * 1.05
    for i, q in enumerate(qqs):
        q.set_xlim(0, qlim)
        q.set_ylim(0, qlim)
        q.plot([0, qlim], [0, qlim], color=S.ACCENT, lw=1.1, zorder=2)
        if i == len(cohorts) - 1:
            q.set_xlabel(f'expected {S.NEGLOG10P}')
    for j in range(len(models)):
        ax = man[-1][j]
        ax.set_xticks(ticks)
        ax.set_xticklabels(names, fontsize=plt.rcParams['xtick.labelsize'] - 2)
        ax.set_xlabel('chromosome')

    # after the final layout: the label band takes its strip out of the top row only, then every
    # Manhattan is brought back to one data height (one shared y-limit needs one height)
    if loci:
        for j in range(len(models)):
            S.gene_labels(man[0][j], guides, [ylim[1]] * len(loci), [d['label'] for d in loci],
                          color=S.INK, weight='bold')
    S.equalise_row_heights(*[ax for row in man for ax in row])
    for i in range(len(cohorts)):
        p, pq = man[i][0].get_position(), qqs[i].get_position()
        qqs[i].set_position([pq.x0, p.y0, pq.width, p.height])
        qqs[i].set_anchor('N')
    for i, c in enumerate(cohorts):
        for j, m in enumerate(models):
            ax = man[i][j]
            pad = S.strip_pad(ax) if i == 0 else 7.0
            st = data[(m, c)]['stats']
            if j == 0:
                S.panel_tag(ax, f'{"abcdefgh"[i]}   {coh[c]["label"]}', size=psize - 1, pad=pad)
            n_txt = f"{st['n_case']:,} cases $\\cdot$ {st['n_ctrl']:,} controls"
            if single:
                n_txt += f" $\\cdot$ $N_{{\\mathrm{{eff}}}}$ = {st['n_eff']:,.0f}"
            ax.set_title(n_txt, loc='right', fontsize=plt.rcParams['legend.fontsize'] - 0.5, color=S.INK_SOFT,
                         pad=pad)
    S.thin_tick_labels(man[-1][0], 'x')
    if not single:
        S.thin_tick_labels(man[-1][1], 'x')

    # column titles: model (in its colour) and what it adjusts for
    y0 = 1 - 0.10 / H
    for j, m in enumerate(models):
        p = man[0][j].get_position()
        M = C.MODELS[m]
        head = titles.get(m, M['title'])
        fig.text(p.x0 + p.width / 2, y0, head, ha='center', va='top', fontweight='bold',
                 fontsize=plt.rcParams['axes.titlesize'] + 1, color=M['color'] if not single else S.INK)
        fig.text(p.x0 + p.width / 2, y0 - 0.19 / H, f'{M["covar"]} · {M["samples"]}', ha='center', va='top',
                 fontsize=plt.rcParams['legend.fontsize'] + 0.5, color=S.INK_SOFT)
    pq = qqs[0].get_position()
    fig.text(pq.x0 + pq.width / 2, y0, 'QQ', ha='center', va='top', fontweight='bold',
             fontsize=plt.rcParams['axes.titlesize'] + 1, color=S.INK)
    fig.savefig(a.out_png)
    plt.close(fig)

    gw_by_panel = {f'{m}:{c}': sorted({label_for(r) for _, r in peaks.get((m, c), pd.DataFrame()).iterrows()})
                   for c in cohorts for m in models}
    with open(a.out_json, 'w') as fh:
        json.dump(dict(models=models, cohorts=cohorts, ylim=ylim[1], qq_limit=qlim, loci=loci,
                       genome_wide_by_panel=gw_by_panel), fh, indent=2)
    C.log(f'{"/".join(models)}: {len(cohorts)} cohorts, loci {[d["label"] for d in loci]} -> {a.out_png}')


if __name__ == '__main__':
    main()
