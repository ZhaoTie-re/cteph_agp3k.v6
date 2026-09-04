#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The scan figure for ONE (cohort, stratum): both statistics drawn
#           against the SAME threshold, their calibration, and their agreement.
#
#           THE THRESHOLD IS READ, NEVER RECOMPUTED. denominator.<stratum>.tsv
#           carries the bonferroni_threshold MERGE_MAP fixed at the map barrier,
#           before any p-value existed, and SCAN stamped significant_bonferroni
#           from it. The value is taken from the file and cross-checked against
#           the one the scan was scored with; a disagreement is fatal.
#
#           THE TWO MANHATTANS SHARE ONE x MAPPING AND ONE y LIMIT, so a gene
#           sits at the same abscissa in both and a taller point is a stronger
#           point. Their label strips are the shared gwaslab idiom: a band above
#           the data, never inside it.
#
#           THE LABEL SET IS THE REJECTION SET, AND NOTHING ELSE. (a) and (b)
#           name every gene that clears at least one decision rule -- BH's
#           rejection set contains Bonferroni's, so `significant_bh` IS that set
#           -- and a panel that rejected nothing names nothing. Naming the k
#           strongest instead states a rule the reader cannot check and puts
#           unqualified genes in the same band as qualified ones.
#
#           THE NAMES ARE ALWAYS VERTICAL. gene_labels' `auto` rotation picks the
#           least intrusive angle that fits, which is 0 deg for one name and 90
#           for ten -- so two panels of the SAME figure were labelled differently.
#           Rotation is pinned so the family reads as one.
#
#           THE QQ CARRIES lambda_GC ONLY. The split by set size is a property of
#           a discrete statistic (anti-conservative near P = 0.5 for two-variant
#           genes) and lives in the family README's table, where it can be read
#           beside the number it explains.
#
#           THIS FIGURE WRITES NO SIDECAR. It is one member of a family of
#           cohort x stratum figures whose prose is identical; it drops its
#           numbers as <png>.stats.json and the CATALOGUE step writes the
#           family's README.md with one row per member.
# Component: assoc_gene
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt        # noqa: E402
import numpy as np                     # noqa: E402
import pandas as pd                    # noqa: E402
from matplotlib.lines import Line2D    # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_shared' / 'scripts'))
import plot_style as S                 # noqa: E402
import figure_doc                      # noqa: E402
import vocab as V                      # noqa: E402

METHOD_A, METHOD_B = V.METHODS
# The data box every Manhattan gets, in inches, in EVERY figure of the family.
# S.gene_labels takes its label band by shrinking the axes, so before this the
# box measured 0.46 in where ten long symbols were named and 0.89 in where three
# were -- the same variant set, half the panel. Here the band is reserved by
# growing that gridspec ROW instead, and the canvas grows with it.
MANHATTAN_H = 1.15
QQ_ROW_H = 2.008           # the square QQ / scatter row; it never grows
HSPACE, TOP_PAD = 0.50, 0.62
BOTTOM_IN = 0.479          # x tick + label slot under the last row; re-measured
PASS1_STRIP = 1.20         # generous, so gene_labels' 55 % band cap cannot bite
MAX_PASSES = 4
LABEL_ROTATION = 90        # pinned: one orientation for the whole family
LABEL_REPEL = 0.045        # at 90 deg the measured gap stops binding; this does
STRIP_TOL = 0.01           # inches; convergence of the achieved data box
CHROM_GAP = 2.0e7
EDGE_PAD = 1.0e7
REL_TOL = 1e-9


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gene-scan', required=True, help='gene_scan.tsv from SCAN')
    p.add_argument('--scan-qc', required=True, help='scan_qc.tsv from SCAN (lambda_gc)')
    p.add_argument('--denominator', required=True,
                   help='denominator.<stratum>.tsv; the Bonferroni line is READ from it')
    p.add_argument('--cohort', required=True)
    p.add_argument('--stratum', required=True)
    p.add_argument('--max-label-genes', type=int, default=24,
                   help='hard cap on the names in one Manhattan; the excess is counted in the '
                        'caption rather than crushed into the label strip')
    p.add_argument('--scan-scale', default=None,
                   help='scan_scale.tsv from SCAN_SCALE: the -log10 P maximum over the three '
                        'cohorts of each variant set, so the three figures of one set share a '
                        'y axis. Absent: the figure falls back to its own maximum.')
    p.add_argument('--out-png', required=True)
    p.add_argument('--out-stats', default=None,
                   help='<png>.stats.json for the family catalogue (default: beside the PNG)')
    return p.parse_args()


# ── Inputs ───────────────────────────────────────────────────────────────────
def read_scan(path, cohort, stratum):
    d = pd.read_csv(path, sep='\t', dtype={'cohort': str, 'stratum': str, 'method': str,
                                           'set_name': str, 'gene_symbol': str, 'chrom': str})
    d = d[(d['cohort'] == cohort) & (d['stratum'] == stratum)].copy()
    if not len(d):
        raise SystemExit(f'ABORT: {path} carries no row for {cohort}/{stratum}')
    bad = set(d['method']) - set(V.METHODS)
    if bad:
        raise SystemExit(f'ABORT: {path} carries method(s) {sorted(bad)}; expected {V.METHODS}')
    for c in ('pvalue', 'pos_min', 'pos_max', 'n_var_map', 'bonferroni_threshold'):
        d[c] = pd.to_numeric(d[c], errors='coerce')
    for c in ('significant_bonferroni', 'significant_bh'):
        d[c] = pd.to_numeric(d[c], errors='coerce').fillna(0).astype(int)
    return d


def read_threshold(path, cohort, stratum, scan):
    d = pd.read_csv(path, sep='\t', dtype={'cohort': str, 'stratum': str})
    if len(d) != 1:
        raise SystemExit(f'ABORT: {path} has {len(d)} data rows; the schema is exactly one')
    row = d.iloc[0]
    if str(row['cohort']) != cohort or str(row['stratum']) != stratum:
        raise SystemExit(f'ABORT: {path} is {row["cohort"]}/{row["stratum"]}, the figure is '
                         f'{cohort}/{stratum}')
    thr = float(row['bonferroni_threshold'])
    if not np.isfinite(thr) or thr <= 0:
        raise SystemExit(f'ABORT: {path} bonferroni_threshold is {row["bonferroni_threshold"]!r}')
    scored = pd.unique(scan['bonferroni_threshold'].dropna())
    bad = [float(v) for v in scored if abs(float(v) - thr) > REL_TOL * thr]
    if bad:
        raise SystemExit(f'ABORT: the scan was scored at {bad[0]:.12g} but the denominator says '
                         f'{thr:.12g}')
    return thr, row


def read_scan_scale(path, stratum):
    """The -log10 P maximum this variant set reached in ANY cohort, or None.

    The three cohorts of one variant set are the same statistic on nested data,
    so a reader compares point heights across them; that only works if the three
    figures share a y axis. The DATA maximum is what is shared -- each figure
    still applies its own Bonferroni floor, because the denominators differ and
    so do the threshold lines.
    """
    if not path:
        return None
    f = Path(path)
    if not f.is_file() or not f.stat().st_size:
        return None
    d = pd.read_csv(f, sep='\t', dtype={'stratum': str})
    d = d[d['stratum'] == stratum]
    if not len(d):
        return None
    v = pd.to_numeric(d['y_data_max'], errors='coerce').max()
    return float(v) if np.isfinite(v) else None


def read_qc(path, cohort, stratum):
    d = pd.read_csv(path, sep='\t', dtype={'cohort': str, 'stratum': str, 'method': str})
    d = d[(d['cohort'] == cohort) & (d['stratum'] == stratum)].copy()
    out = {}
    for _, r in d.iterrows():
        vals = {}
        for k in ('lambda_gc', 'lambda_gc_nvar2', 'lambda_gc_nvar_ge3', 'bh_threshold'):
            v = pd.to_numeric(r[k], errors='coerce')
            vals[k] = float(v) if np.isfinite(v) else None
        out[r['method']] = vals
    return out


def cell(scan, method):
    d = scan[scan['method'] == method].copy()
    d = d[np.isfinite(d['pvalue']) & (d['pvalue'] <= 1.0)]
    d['_y'] = V.nlp(d['pvalue'].to_numpy(dtype=float))
    return d


def chrom_key(c):
    c = str(c).replace('chr', '')
    return {'X': 23, 'Y': 24, 'MT': 25, 'M': 25}.get(c, int(c) if c.isdigit() else 99)


def genome_axis(scan):
    """One x position per gene, built from EVERY method's genes at once."""
    g = (scan.dropna(subset=['pos_min', 'pos_max'])
             .drop_duplicates('set_name')[['set_name', 'chrom', 'pos_min', 'pos_max']].copy())
    if not len(g):
        raise SystemExit('ABORT: no gene in the scan carries a position')
    g['_k'] = g['chrom'].map(chrom_key)
    g['_mid'] = (g['pos_min'] + g['pos_max']) / 2.0
    g = g.sort_values(['_k', '_mid'])
    xpos, band, ticks, names, offset, lo, hi = {}, {}, [], [], 0.0, None, None
    for i, (_k, sub) in enumerate(g.groupby('_k', sort=True)):
        xs = offset + sub['_mid'].to_numpy(dtype=float)
        xpos.update(dict(zip(sub['set_name'], xs)))
        band[str(sub['chrom'].iloc[0])] = i % 2
        ticks.append((xs[0] + xs[-1]) / 2.0)
        names.append(str(sub['chrom'].iloc[0]))
        lo = xs[0] if lo is None else lo
        hi = xs[-1]
        offset = xs[-1] + CHROM_GAP
    return xpos, band, ticks, names, (lo - EDGE_PAD, hi + EDGE_PAD)


# ── Panels ───────────────────────────────────────────────────────────────────
def threshold_line(ax, thr, x, text, colour, dashed, ha):
    """A decision rule, drawn once and named with its own value.

    The label's zorder is ABOVE the gene-label leader arms (zorder 9) and its
    backing is opaque: an arm crossing the threshold line is unavoidable -- the
    names live in a strip above the axes and their anchors are inside it -- so
    the label wins the pixel rather than being struck through.
    """
    y = -np.log10(thr)
    ax.axhline(y, color=colour, lw=1.0, zorder=4, ls=(0, (4, 2)) if dashed else '-')
    ax.text(x, y, f'{text} {S.p_tex(thr)}', color=colour, ha=ha, va='bottom',
            fontsize=V.MIN_FONT, zorder=11,
            bbox=dict(boxstyle='square,pad=0.15', facecolor='white', edgecolor='none',
                      alpha=1.0))


def manhattan(ax, d, geom, thr, ylim, xlabel=None, bh_thr=None):
    xpos, band, ticks, names, xlim = geom
    x = d['set_name'].map(xpos)
    ok = x.notna()
    d, x = d[ok], x[ok].to_numpy(dtype=float)
    colours = [S.CHROM_BANDS[band.get(str(c), 0)] for c in d['chrom']]
    ax.scatter(x, d['_y'].to_numpy(), s=6.0, c=colours, linewidths=0, rasterized=True, zorder=2)
    bonf = (d['significant_bonferroni'] == 1).to_numpy()
    bh_only = ((d['significant_bh'] == 1) & (d['significant_bonferroni'] != 1)).to_numpy()
    if bh_only.any():
        ax.scatter(x[bh_only], d.loc[bh_only, '_y'].to_numpy(), s=12.0, c=S.ACCENT_LT,
                   linewidths=0.4, edgecolors='white', zorder=5)
    if bonf.any():
        ax.scatter(x[bonf], d.loc[bonf, '_y'].to_numpy(), s=14.0, c=S.ACCENT,
                   linewidths=0.4, edgecolors='white', zorder=6)
    span = xlim[1] - xlim[0]
    # The BH cut can only lie at or below the Bonferroni line (BH's largest
    # rejected P is >= alpha / n), so labelling it at the LEFT edge and Bonferroni
    # at the RIGHT keeps the two from stacking however close the lines are.
    if bh_thr is not None and np.isfinite(bh_thr) and 0 < bh_thr <= 1:
        threshold_line(ax, bh_thr, xlim[0] + 0.004 * span, 'BH q < 0.05', S.ACCENT_LT,
                       True, 'left')
    else:
        # BH rejected nothing, so it HAS no cut to draw -- its largest rejected P
        # does not exist. Drawing neither line nor label left the panel looking as
        # though the second rule had been forgotten, so it is stated instead, on
        # the left where the BH label would have gone and below the Bonferroni
        # line so it cannot reach that label or the name band above the axes.
        ax.text(xlim[0] + 0.004 * span, -np.log10(thr), 'BH q < 0.05: no gene called',
                color=S.ACCENT_LT, ha='left', va='top', fontsize=V.MIN_FONT, zorder=11,
                bbox=dict(boxstyle='square,pad=0.15', facecolor='white', edgecolor='none',
                          alpha=1.0))
    threshold_line(ax, thr, xlim[1] - 0.004 * span, 'Bonferroni', S.ACCENT, False, 'right')
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks(ticks)
    ax.set_xticklabels(names)
    ax.set_ylabel(S.NEGLOG10P)
    if xlabel:
        ax.set_xlabel(xlabel)
    S.despine(ax, grid_axis='y')
    return d, x


def name_genes(ax, d, x, n_cap):
    """Name every called gene. A panel that called nothing names nothing.

    Returns (names, capped) so the caption can say whether the cap bit. Every
    name here clears at least one decision rule, so the two tiers are the only
    ones that exist: BH+Bonferroni in the full accent and bold, BH only in the
    lighter one.

    LABEL_ROTATION is pinned rather than left to gene_labels' `auto`, which
    measures and would give a one-name panel 0 deg and a ten-name panel 90 --
    two panels of one figure labelled differently. Note the side effect: at 90
    deg the measured per-label footprint is the text HEIGHT (~0.026 of the span),
    which falls below repel_force, so repel_force alone sets the spacing. It is
    raised here to hold the names apart, and gene_labels prints a WARNING to
    stderr if any adjacent pair still overlaps.
    """
    if not len(d):
        return [], False
    take = d[d['significant_bh'] == 1]
    if not len(take):
        return [], False
    take = take.assign(_x=x[(d['significant_bh'] == 1).to_numpy()])
    capped = len(take) > n_cap
    if capped:
        take = take.nsmallest(n_cap, 'pvalue')
    bonf = (take['significant_bonferroni'] == 1).to_numpy()
    S.gene_labels(ax, take['_x'].to_numpy(), take['_y'].to_numpy(),
                  [str(g) for g in take['gene_symbol']],
                  rotation=LABEL_ROTATION, repel_force=LABEL_REPEL,
                  color=[S.ACCENT if b else S.ACCENT_LT for b in bonf],
                  weight=['bold' if b else 'normal' for b in bonf])
    return list(take['gene_symbol'].astype(str)), capped


def qq(ax, cells, qc):
    handles, lim = [], 1.0
    n_max = max((len(d) for _m, d in cells), default=0)
    if n_max:
        exp, lo, hi = V.qq_band(n_max)
        ax.fill_between(exp, lo, hi, color=S.NEUTRAL, alpha=0.55, lw=0, zorder=1)
    for method, d in cells:
        if not len(d):
            continue
        obs = np.sort(d['_y'].to_numpy())[::-1]
        exp = V.qq_expected(len(obs))
        ax.plot(exp, obs, color=V.METHOD_COLOR[method], lw=1.1, zorder=3, rasterized=True)
        lim = max(lim, float(exp.max()), float(obs.max()))
        lam = qc.get(method, {}).get('lambda_gc')
        lam_s = '—' if lam is None else f'{lam:.2f}'
        handles.append(Line2D([], [], color=V.METHOD_COLOR[method], lw=1.6,
                              label=f'{V.METHOD_LABEL[method]}   {S.LAMBDA_GC} = {lam_s}'))
    handles.append(Line2D([], [], color=S.NEUTRAL, lw=6, alpha=0.8, label='null 95 % band'))
    lim *= 1.06
    ax.plot([0, lim], [0, lim], color=S.REFERENCE, lw=0.9, zorder=2)
    ax.set_box_aspect(1)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel(f'expected {S.NEGLOG10P}')
    ax.set_ylabel(f'observed {S.NEGLOG10P}')
    S.despine(ax)
    S.legend_inside(ax, handles, loc='lower right', ncol=1)


def method_scatter(ax, da, db, thr, bh_a, bh_b):
    """The same gene under both statistics, matched on set_name."""
    j = da[['set_name', 'gene_symbol', '_y', 'significant_bh']].merge(
        db[['set_name', '_y', 'significant_bh']], on='set_name', suffixes=('_a', '_b'))
    y_thr = -np.log10(thr)
    if not len(j):
        ax.text(0.5, 0.5, 'no gene returned by both statistics', transform=ax.transAxes,
                ha='center', va='center', color=S.INK_SOFT)
        S.despine(ax)
        return j, [], 1.0
    x, y = j['_y_a'].to_numpy(), j['_y_b'].to_numpy()
    lim = max(float(np.nanmax(x)), float(np.nanmax(y)), y_thr) * 1.12
    ax.plot([0, lim], [0, lim], color=S.REFERENCE, lw=0.9, ls=(0, (3, 2)), zorder=2)
    either = ((j['significant_bh_a'] == 1) | (j['significant_bh_b'] == 1)).to_numpy()
    ax.scatter(x[~either], y[~either], s=5.0, c=S.NEUTRAL_D, linewidths=0,
               rasterized=True, zorder=3)
    ax.scatter(x[either], y[either], s=14.0, c=S.ACCENT, linewidths=0.4, edgecolors='white',
               zorder=5)
    ax.axvline(y_thr, color=S.ACCENT, lw=1.0, zorder=4)
    ax.axhline(y_thr, color=S.ACCENT, lw=1.0, zorder=4)
    # Each statistic's own BH cut, on its own axis.
    if bh_a is not None and 0 < bh_a <= 1:
        ax.axvline(-np.log10(bh_a), color=S.ACCENT_LT, lw=0.9, ls=(0, (4, 2)), zorder=4)
    if bh_b is not None and 0 < bh_b <= 1:
        ax.axhline(-np.log10(bh_b), color=S.ACCENT_LT, lw=0.9, ls=(0, (4, 2)), zorder=4)
    ax.set_box_aspect(1)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel(f'{V.METHOD_LABEL[METHOD_A]} {S.NEGLOG10P}')
    ax.set_ylabel(f'{V.METHOD_LABEL[METHOD_B]} {S.NEGLOG10P}')
    S.despine(ax)
    # NO gene names here: they are already in (a) and (b), where they sit against
    # the genome and the thresholds. The accent is the statement this panel makes
    # -- a gene called by at least one rule under at least one statistic -- and it
    # is counted in the corner above the diagonal cloud, which is empty.
    n_either = int(either.sum())
    # The corner is CHOSEN by counting the points in each one, not assumed: the
    # cloud hugs the diagonal and the called genes sit above it, so which corner
    # is free depends on the scan.
    fx, fy = x / lim, y / lim
    corners = {(0.03, 0.97, 'left', 'top'): ((fx < 0.35) & (fy > 0.65)).sum(),
               (0.97, 0.97, 'right', 'top'): ((fx > 0.65) & (fy > 0.65)).sum(),
               (0.97, 0.03, 'right', 'bottom'): ((fx > 0.65) & (fy < 0.35)).sum(),
               (0.03, 0.03, 'left', 'bottom'): ((fx < 0.35) & (fy < 0.35)).sum()}
    cx, cy, ha, va = min(corners, key=corners.get)
    ax.text(cx, cy, f'{n_either} gene(s) called by\neither statistic',
            transform=ax.transAxes, ha=ha, va=va, fontsize=V.MIN_FONT,
            color=S.ACCENT, linespacing=1.3, zorder=9,
            bbox=dict(boxstyle='square,pad=0.2', facecolor='white', edgecolor='none',
                      alpha=0.85))
    return j, n_either, lim


def spearman(d):
    if len(d) < 3:
        return None
    r = pd.Series(d['n_var_map'].to_numpy(dtype=float)).corr(pd.Series(d['_y'].to_numpy()),
                                                              method='spearman')
    return float(r) if np.isfinite(r) else None


# ── Figure ────────────────────────────────────────────────────────────
def build_figure(D, *, strips, bottom_in):
    """One render of the four panels; `strips` reserves the gene-label bands.

    S.gene_labels takes its band by SHRINKING the axes, so the data box of a
    scan that names ten long symbols ends up half that of one that names three.
    Here the band is paid for by growing that gridspec ROW -- in absolute inches,
    so the two square panels below never move -- and the canvas grows with it.
    Every Manhattan in the family then has the same MANHATTAN_H data box.

    Returns the figure and what the render measured, for the next pass.
    """
    args, qc = D['args'], D['qc']
    rows = [MANHATTAN_H + strips[0], MANHATTAN_H + strips[1], QQ_ROW_H]
    # matplotlib puts (n-1) gaps of hspace x mean(row) between n rows.
    plot_h = TOP_PAD + sum(rows) * (1.0 + 2.0 * HSPACE / len(rows)) + bottom_in

    fig = plt.figure(figsize=(S.COL_DOUBLE, plot_h))
    gs = fig.add_gridspec(3, 2, height_ratios=rows, width_ratios=[1.0, 1.0])
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, :])
    ax_c = fig.add_subplot(gs[2, 0])
    ax_d = fig.add_subplot(gs[2, 1])

    bh_a, bh_b = D['bh_a'], D['bh_b']
    da_p, xa = manhattan(ax_a, D['da'], D['geom'], D['thr'], D['ylim'], bh_thr=bh_a)
    db_p, xb = manhattan(ax_b, D['db'], D['geom'], D['thr'], D['ylim'],
                         xlabel='chromosome', bh_thr=bh_b)
    qq(ax_c, ((METHOD_A, D['da']), (METHOD_B, D['db'])), qc)
    joint, n_either, _lim = method_scatter(ax_d, D['da'], D['db'], D['thr'], bh_a, bh_b)

    S.caption_block(
        fig, plot_h=plot_h,
        title=D['title'], panels=D['panels'], notes=D['notes'],
        letters='abcdefghijklmnop',
        top_pad=TOP_PAD, hspace=HSPACE, wspace=0.45, left='auto', right=0.93,
        margin_axes=[ax_a, ax_b, ax_c, ax_d])

    # ---- measured layout from here on.
    fig_h = fig.get_size_inches()[1]
    # The slack under the last row (its tick labels and axis label). Measured
    # rather than assumed, because it follows the type scale, and fed to the
    # next pass so the row heights it asks for are the ones it gets.
    bottom_measured = plot_h - (1.0 - ax_c.get_position().y0) * fig_h
    named_a, capped_a = name_genes(ax_a, da_p, xa, args.max_label_genes)
    named_b, capped_b = name_genes(ax_b, db_p, xb, args.max_label_genes)
    # What each band actually cost, and whether gene_labels' 55 % ceiling clipped
    # it -- a clipped band UNDERSTATES the inches the names need, so a pass that
    # hit the ceiling must not be trusted to size the next one.
    strips_req, hit_cap = [], []
    for ax in (ax_a, ax_b):
        frac = getattr(ax, '_label_strip_frac', 0.0)
        h_after = ax.get_position().height
        strips_req.append(frac * fig_h)
        hit_cap.append(frac >= 0.55 * (h_after + frac) - 1e-9)
    # SNAP, not S.equalise_row_heights. That helper brings a stack to its common
    # MINIMUM, which is right when the row height is whatever the layout produced;
    # here it is a family constant, and a shrink below it is the very defect this
    # pass removes. Each panel keeps its bottom edge, and the guarantee is then
    # measured rather than inferred from the gridspec arithmetic.
    want_frac = MANHATTAN_H / fig_h
    for ax in (ax_a, ax_b):
        pos = ax.get_position()
        if abs(pos.height - want_frac) > 1e-9:
            ax.set_position([pos.x0, pos.y0, pos.width, want_frac])
    boxes = [ax.get_position().height * fig_h for ax in (ax_a, ax_b)]
    for ax, letter, ttl in ((ax_a, 'a', V.METHOD_LABEL[METHOD_A]),
                            (ax_b, 'b', V.METHOD_LABEL[METHOD_B])):
        S.panel_tag(ax, letter, title=ttl, pad=S.strip_pad(ax))
    S.panel_tag(ax_c, 'c')
    S.panel_tag(ax_d, 'd')
    for ax in (ax_a, ax_b):
        S.thin_tick_labels(ax, 'x')
    return fig, dict(strips=tuple(strips_req), hit_cap=tuple(hit_cap), boxes=tuple(boxes),
                     bottom_in=bottom_measured, plot_h=plot_h,
                     named_a=named_a, named_b=named_b, capped=bool(capped_a or capped_b),
                     joint=joint, n_either=n_either)


def main():
    args = parse_args()
    S.setup_style('paper')
    scan = read_scan(args.gene_scan, args.cohort, args.stratum)
    thr, den = read_threshold(args.denominator, args.cohort, args.stratum, scan)
    qc = read_qc(args.scan_qc, args.cohort, args.stratum)
    da, db = cell(scan, METHOD_A), cell(scan, METHOD_B)
    for method, d in ((METHOD_A, da), (METHOD_B, db)):
        if not len(d):
            raise SystemExit(f'ABORT: {args.gene_scan} has no usable p-value for {method}')
    clabel = V.cohort_labels(V.COHORT_PREFERENCE).get(args.cohort, args.cohort)
    stratum_h = V.stratum_label(args.stratum)

    geom = genome_axis(scan)
    y_thr = -np.log10(thr)
    own_top = max(float(da['_y'].max()), float(db['_y'].max()))
    shared_top = read_scan_scale(args.scan_scale, args.stratum)
    data_top = max(own_top, shared_top) if shared_top is not None else own_top
    ytop = max(data_top, y_thr) * 1.14
    ylim = (0.0, max(ytop, y_thr + 0.6))

    bh_a = qc.get(METHOD_A, {}).get('bh_threshold')
    bh_b = qc.get(METHOD_B, {}).get('bh_threshold')
    n_sig_a = int((da['significant_bh'] == 1).sum())
    n_sig_b = int((db['significant_bh'] == 1).sum())
    n_bonf_a = int((da['significant_bonferroni'] == 1).sum())
    n_bonf_b = int((db['significant_bonferroni'] == 1).sum())
    top_a = da.nsmallest(1, 'pvalue').iloc[0]
    top_b = db.nsmallest(1, 'pvalue').iloc[0]

    def label_clause(method, n_drawn, n_sig):
        if not n_sig:
            return (f'{V.METHOD_LABEL[method]}, {n_drawn:,} genes; none called, so none '
                    f'named.')
        shown = min(n_sig, args.max_label_genes)
        tail = f' ({shown} of them named)' if n_sig > args.max_label_genes else ''
        return (f'{V.METHOD_LABEL[method]}, {n_drawn:,} genes; every called gene '
                f'named{tail}.')

    # A gene called by both statistics is in both rejection sets; the joint frame
    # method_scatter builds is an inner merge of the same two, so the intersection
    # over set names is that count and needs no render to compute.
    n_both = len(set(da.loc[da['significant_bh'] == 1, 'set_name'])
                 & set(db.loc[db['significant_bh'] == 1, 'set_name']))
    scale_note = (f' The y axis is this variant set\'s maximum over every cohort, so its '
                  f'figures are on one scale; the threshold lines are this cohort\'s own.'
                  if shared_top is not None else '')
    D = dict(args=args, qc=qc, da=da, db=db, geom=geom, thr=thr, ylim=ylim,
             bh_a=bh_a, bh_b=bh_b,
             title=(f'{clabel}, {stratum_h}: '
                    + (f'{n_sig_a} gene(s) called by CMC and {n_sig_b} by SKAT-O'
                       if n_sig_a or n_sig_b else 'no gene is called by either statistic')
                    + (f', {n_both} by both' if (n_sig_a or n_sig_b) else '')
                    + f'; strongest {top_a["gene_symbol"]} (CMC), '
                      f'{top_b["gene_symbol"]} (SKAT-O).'),
             panels=[
                 label_clause(METHOD_A, len(da), n_sig_a),
                 label_clause(METHOD_B, len(db), n_sig_b) + ' Same axes as (a).',
                 'Calibration of both statistics with the null 95 % band.',
                 'CMC against SKAT-O per gene; red = called by at least one rule under at '
                 'least one statistic.',
             ],
             notes=(f'{args.cohort}, {stratum_h}; {int(den["n_genes_mapped"]):,} genes with '
                    f'≥{int(den["min_num_var"])} variants; Bonferroni {S.p_tex(thr)} = '
                    f'{float(den["alpha"]):g}/{int(den["n_genes_mapped"]):,}.{scale_note} '
                    f'README.md'))

    # Render until every Manhattan really is MANHATTAN_H tall. The first pass
    # deliberately over-reserves: gene_labels clips a band at 55 % of its panel,
    # and a clipped band understates what the names need, so a pass that hit the
    # ceiling must not be trusted to size the next one.
    strips, bottom_in, fig, info = (PASS1_STRIP, PASS1_STRIP), BOTTOM_IN, None, None
    done = False
    for _pass in range(MAX_PASSES):
        if fig is not None:
            plt.close(fig)
        fig, info = build_figure(D, strips=strips, bottom_in=bottom_in)
        bottom_in = info['bottom_in']
        want = tuple(round(x, 4) for x in info['strips'])
        done = (not any(info['hit_cap'])
                and all(abs(b - MANHATTAN_H) <= STRIP_TOL for b in info['boxes'])
                and strips == want)
        # A band the ceiling clipped is a lower bound on what the names need, so
        # asking for the clipped value again would only clip it again.
        strips = tuple(max(w, cur * 1.6) if cap else w
                       for w, cur, cap in zip(want, strips, info['hit_cap']))
        if done:
            break
    if not done:
        print(f'[plot_gene_scan] WARNING: the gene-label bands did not settle in '
              f'{MAX_PASSES} passes; data boxes '
              f'{tuple(round(b, 4) for b in info["boxes"])} in, want {MANHATTAN_H} in',
              file=sys.stderr)
    named_a, named_b = info['named_a'], info['named_b']
    n_either, capped = info['n_either'], info['capped']
    fig.savefig(args.out_png)
    plt.close(fig)

    # ---- the numbers of this member, for the family catalogue.
    n_no_p = int(pd.to_numeric(scan['pvalue'], errors='coerce').isna().sum())
    values = [
        ('cohort', args.cohort), ('stratum', args.stratum),
        ('genes mapped', int(den['n_genes_mapped'])),
        ('MinNumVar', int(den['min_num_var'])),
        ('Bonferroni threshold', thr),
        ('BH cut, CMC', bh_a), ('BH cut, SKAT-O', bh_b),
        ('called (BH), CMC', n_sig_a), ('of which Bonferroni, CMC', n_bonf_a),
        ('called (BH), SKAT-O', n_sig_b), ('of which Bonferroni, SKAT-O', n_bonf_b),
        ('called by both', n_both),
        ('strongest, CMC', f'{top_a["gene_symbol"]} ({S.sci(float(top_a["pvalue"]))})'),
        ('strongest, SKAT-O', f'{top_b["gene_symbol"]} ({S.sci(float(top_b["pvalue"]))})'),
        ('lambda_GC, CMC', qc.get(METHOD_A, {}).get('lambda_gc')),
        ('lambda_GC n_var=2, CMC', qc.get(METHOD_A, {}).get('lambda_gc_nvar2')),
        ('lambda_GC n_var>=3, CMC', qc.get(METHOD_A, {}).get('lambda_gc_nvar_ge3')),
        ('lambda_GC, SKAT-O', qc.get(METHOD_B, {}).get('lambda_gc')),
        ('lambda_GC n_var=2, SKAT-O', qc.get(METHOD_B, {}).get('lambda_gc_nvar2')),
        ('lambda_GC n_var>=3, SKAT-O', qc.get(METHOD_B, {}).get('lambda_gc_nvar_ge3')),
        ('Spearman rho set size vs signal, CMC', spearman(da)),
        ('Spearman rho set size vs signal, SKAT-O', spearman(db)),
        ('genes named in (a)', ', '.join(named_a)),
        ('genes named in (b)', ', '.join(named_b)),
        ('label rule', 'every called gene, nothing else'),
        ('name cap reached', capped),
        ('genes called by either statistic (d)', n_either),
        ('gene tests with no P', n_no_p),
    ]
    stats_path = args.out_stats or str(Path(args.out_png).with_suffix('.stats.json'))
    figure_doc.write_stats(Path(stats_path).with_suffix('').with_suffix('.png')
                           if stats_path.endswith('.stats.json') else args.out_png,
                           peak=Path(args.out_png).stem, values=values)
    print(f'[plot_gene_scan] {args.cohort} {args.stratum}: CMC {len(da):,} genes '
          f'({n_sig_a} called), SKAT-O {len(db):,} ({n_sig_b} called), threshold {thr:.6g} '
          f'-> {args.out_png}, {stats_path}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:                       # noqa: BLE001
        print(f'Error in plot_gene_scan: {e}', file=sys.stderr)
        sys.exit(1)
