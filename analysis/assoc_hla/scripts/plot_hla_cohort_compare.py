#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : One HLA marker class, one model, the three NESTED sample sets.
#             (a) evidence  — -log10 P per cohort on one shared marker axis
#             (b) effect    — log(OR) between every pair of cohorts
#             (c) the hits  — OR (95% CI), the cohorts stacked per marker
#             (d) summary   — how many hits each cohort has, and which are its own
#
#           THE COHORTS ARE NESTED: narrow subset of intermediate subset of full.
#           They share every case and most controls, so this figure is NOT a
#           replication analysis and no panel of it may be read as one. What it
#           reports is how one estimate MOVES as the ancestry filter is relaxed —
#           an interval that narrows while the point estimate holds is a marker
#           gaining sample size; an estimate that drifts toward OR = 1 as the
#           filter opens is what a structure-driven signal does. The caption says
#           this in as many words, because a three-cohort forest is the exact
#           shape a reader has been trained to read as replication.
#
#           HLA MAKES THE NESTING WORSE, NOT BETTER. Class I and II allele
#           frequencies differ measurably between Japanese subpopulations, so the
#           ancestry filter that defines the three sets acts directly on the
#           markers being tested. A marker significant only in `narrow` is
#           therefore the case that most needs explaining, and (d) is what makes
#           those markers countable rather than something the reader must find by
#           eye in (c).
#
#           THE MARKER AXIS IS ORDERED BY GENE, THEN POSITION, and is the SAME
#           axis for all three cohorts, so a marker sits at one x in every series
#           and the panel can be read down as well as across. HLA markers have no
#           meaningful base-pair spacing — build_hla_markers.py mints a coordinate
#           per gene — so this is an ordered categorical axis, drawn as one, and
#           never as a genomic distance.
# Component: assoc_hla
# ---------------------------------------------------------------------------
import argparse
import itertools
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402
from matplotlib.lines import Line2D      # noqa: E402
from matplotlib.patches import Patch     # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
SHARED = Path(__file__).resolve().parent.parent.parent / '_shared' / 'scripts'
sys.path.insert(0, str(SHARED))
import plot_style as S                   # noqa: E402
import figure_doc                        # noqa: E402

NOMINAL = 0.05          # the dashed second line on the -log10 P axis
ROW_IN = 0.185          # inches per forest row
MAX_FOREST_MARKERS = 12  # a cap on the FIGURE only; --out-table holds them all
MAX_AXIS_LABELS = 8     # named on panel (a); more than this is a wall of text
MAX_NAMES_IN_D = 4      # cohort-specific markers named in (d) before "+N more"
CANON_NUM = ('POS', 'A1_FREQ', 'OBS_CT', 'OR', 'LOG(OR)_SE', 'L95', 'U95', 'Z_STAT', 'P')

# Three states, shape-encoded, so COLOUR stays free to carry the cohort — the
# same division of labour as plot_style.CALLED_STYLE, with this component's own
# wording: nothing here is "genome-wide", and saying so would import a threshold
# from an experiment that was not run.
STATE = {
    'significant':  dict(marker='D', size=34, filled=True,
                         label='significant in this cohort'),
    'not_significant': dict(marker='o', size=22, filled=False,
                            label='not significant there'),
    'no_fit':       dict(marker='x', size=26, filled=False,
                         label='no usable fit there'),
}


def parse_args():
    p = argparse.ArgumentParser(
        description='One HLA marker class and model across the three NESTED cohorts.')
    p.add_argument('--sumstat', action='append', required=True, metavar='COHORT=PATH',
                   help='repeatable; one sumstats.tsv per cohort, all of the same model '
                        'and marker class')
    p.add_argument('--cohort-order', default='',
                   help='comma-separated, in NESTING order (narrow first). The colour ramp '
                        'is ordered, so this decides which cohort is dark.')
    p.add_argument('--marker-class', required=True, choices=('allele', 'residue'))
    p.add_argument('--model', required=True, help='the single model these scans are from')
    p.add_argument('--model-label', default='',
                   help='how that model is named on the figure; defaults to --model')
    p.add_argument('--alpha', type=float, default=0.0,
                   help='significance threshold. 0 (the default) means Bonferroni 0.05 / M '
                        'over the M distinct markers with a usable fit in this class.')
    p.add_argument('--out-png', required=True)
    p.add_argument('--out-table', required=True, help='one row per (marker, cohort)')
    return p.parse_args()


def spec(items):
    """['cohort=path', …] -> {cohort: path}, existing non-empty files only."""
    out = {}
    for s in items or []:
        c, _, path = s.partition('=')
        if c and path and Path(path).exists() and Path(path).stat().st_size:
            out[c] = path
    return out


def read_sumstat(path, cohort):
    """One sumstats.tsv as a tidy frame; only rows with ERRCODE '.' and 0 < P <= 1."""
    d = pd.read_csv(path, sep='\t', dtype={'ID': str, 'ERRCODE': str, 'gene': str,
                                           'position': str, '#CHROM': str})
    for c in CANON_NUM:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors='coerce')
    for c in ('ID', 'gene', 'position'):
        if c not in d.columns:
            d[c] = ''
    err = d['ERRCODE'].fillna('.').astype(str).str.strip() if 'ERRCODE' in d.columns \
        else pd.Series('.', index=d.index)
    ok = err.eq('.') & d['P'].gt(0) & d['P'].le(1)
    out = pd.DataFrame({
        'marker': d['ID'].astype(str),
        'gene': d['gene'].fillna('').astype(str),
        'position': d['position'].fillna('').astype(str),
        'pos': pd.to_numeric(d['POS'], errors='coerce') if 'POS' in d.columns else np.nan,
        'cohort': cohort,
        'OR': d['OR'], 'L95': d.get('L95'), 'U95': d.get('U95'),
        'SE': d.get('LOG(OR)_SE'), 'A1_FREQ': d.get('A1_FREQ'),
        'OBS_CT': d.get('OBS_CT'), 'P': d['P']})
    out['log_or'] = np.log(out['OR'].where(out['OR'] > 0))
    return out[ok].reset_index(drop=True), int((~ok).sum())


def empty_panel(ax, msg):
    """An honest blank: say why there is nothing rather than draw an empty box."""
    ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha='center', va='center',
            color=S.INK_SOFT, fontsize=plt.rcParams['legend.fontsize'])
    ax.set_xticks([]); ax.set_yticks([])
    ax.grid(False)
    # No frame either: an empty axes box reads as a rendering failure, while a
    # bare sentence reads as the statement it is.
    for sp in ax.spines.values():
        sp.set_visible(False)


def marker_axis(long):
    """Marker -> x index, ordered by GENE then POSITION; plus the per-gene blocks.

    `POS` is the coordinate build_hla_markers.py mints per gene (its own sort
    order within the gene); `position` is the residue position id and is blank
    for an allele marker, so it cannot be the primary key. Ties fall back to the
    marker id, which makes the order deterministic run to run.
    """
    if not len(long):
        return {}, []
    m = (long.groupby('marker', as_index=False)
         .agg(gene=('gene', 'first'), pos=('pos', 'first'), position=('position', 'first')))
    m['pos'] = pd.to_numeric(m['pos'], errors='coerce')
    m = m.sort_values(['gene', 'pos', 'position', 'marker'],
                      na_position='last').reset_index(drop=True)
    xof = {mk: i for i, mk in enumerate(m['marker'])}
    blocks = [(g, int(idx.min()), int(idx.max()))
              for g, idx in m.groupby('gene', sort=False).groups.items()]
    return xof, blocks


def main():
    args = parse_args()
    S.setup_style('paper')
    ss = spec(args.sumstat)
    if not ss:
        raise SystemExit('ABORT: no --sumstat file exists and is non-empty')
    cohorts = [c for c in args.cohort_order.split(',') if c] or list(ss)
    cohorts = [c for c in cohorts if c in ss]
    if not cohorts:
        raise SystemExit('ABORT: --cohort-order selects no cohort that was given')
    SHORT = S.shorten(cohorts)
    ccol = {c: S.COHORT_RAMP[i % len(S.COHORT_RAMP)] for i, c in enumerate(cohorts)}
    MC = args.marker_class
    MODEL = args.model_label or args.model

    frames, n_dropped = [], 0
    for c in cohorts:
        d, bad = read_sumstat(ss[c], c)
        n_dropped += bad
        frames.append(d)
    long = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    n_markers = int(long['marker'].nunique()) if len(long) else 0
    alpha = args.alpha if args.alpha > 0 else (0.05 / max(n_markers, 1))
    long['significant'] = long['P'] < alpha if len(long) else pd.Series(dtype=bool)

    # How many cohorts called each marker, and in how many it was even fitted —
    # the two counts panel (d) is built from, and the honest denominator for
    # "shared": a marker cannot be shared with a cohort that never fitted it.
    n_sig_by_marker = (long[long.significant].groupby('marker')['cohort'].nunique()
                       if len(long) else pd.Series(dtype=int))
    n_fit_by_marker = (long.groupby('marker')['cohort'].nunique()
                       if len(long) else pd.Series(dtype=int))
    sig_markers = list(n_sig_by_marker.index)
    best = (long[long.marker.isin(sig_markers)].groupby('marker')['P'].min().sort_values()
            if sig_markers else pd.Series(dtype=float))
    drawn = list(best.index[:MAX_FOREST_MARKERS])
    capped = len(best) > len(drawn)

    # ── the table ───────────────────────────────────────────────────────────
    keep = list(best.index)
    if not keep and len(long):
        keep = list(long.groupby('marker')['P'].min().sort_values().index[:10])
    tab_cols = ['marker', 'gene', 'position', 'cohort', 'model', 'OR', 'L95', 'U95',
                'log_or', 'SE', 'P', 'A1_FREQ', 'OBS_CT', 'significant',
                'n_cohorts_significant', 'n_cohorts_fitted', 'shared_by_all_fitted']
    if len(long):
        tab = long[long.marker.isin(keep)].copy()
        tab['model'] = args.model
        tab['n_cohorts_significant'] = tab.marker.map(n_sig_by_marker).fillna(0).astype(int)
        tab['n_cohorts_fitted'] = tab.marker.map(n_fit_by_marker).fillna(0).astype(int)
        tab['shared_by_all_fitted'] = (tab.n_cohorts_significant == tab.n_cohorts_fitted)
        tab['_r'] = tab.marker.map({m: i for i, m in enumerate(keep)})
        tab['_c'] = tab.cohort.map({c: i for i, c in enumerate(cohorts)})
        tab = tab.sort_values(['_r', '_c'])[tab_cols]
    else:
        tab = pd.DataFrame(columns=tab_cols)
    tab.to_csv(args.out_table, sep='\t', index=False, na_rep='NA')

    # ── canvas: the forest is sized from its own row count ─────────────────
    rows = [(mk, c) for mk in drawn for c in cohorts]
    pairs = list(itertools.combinations(range(len(cohorts)), 2))
    n_named = min(len(best), MAX_AXIS_LABELS)
    # The label band is reserved by SHRINKING the axes (gene_labels never inflates
    # the data limits), so the row has to be given the height the band will take
    # or the axes is left half the size of its neighbours. Vertical marker ids at
    # this type scale need about an inch.
    h_a = 2.05 + (0.95 if n_named else 0.0)
    h_b = 2.15
    h_c = max(0.85, len(rows) * ROW_IN) + 0.70
    h_d = 1.30
    plot_h = h_a + h_b + h_c + h_d + 1.70

    fig = plt.figure(figsize=(S.COL_DOUBLE, plot_h))
    gs = fig.add_gridspec(4, 1, height_ratios=[h_a, h_b, h_c, h_d])
    ax_a = fig.add_subplot(gs[0])
    gsb = gs[1].subgridspec(1, max(len(pairs), 1), wspace=0.34)
    ax_b = [fig.add_subplot(gsb[0, i]) for i in range(max(len(pairs), 1))]
    ax_c = fig.add_subplot(gs[2])
    ax_d = fig.add_subplot(gs[3])

    # ── (a) one shared marker axis, one series per cohort ──────────────────
    xof, blocks = marker_axis(long)
    label_rows = []
    if not xof:
        empty_panel(ax_a, 'no marker with a usable fit')
    else:
        # Alternating gene bands: structure, so NEUTRAL, and behind everything.
        for k, (g, lo, hi) in enumerate(blocks):
            if k % 2:
                ax_a.axvspan(lo - 0.5, hi + 0.5, color=S.NEUTRAL, alpha=0.45,
                             linewidth=0, zorder=0)
        for c in cohorts:
            sub = long[long.cohort == c]
            if not len(sub):
                continue
            x = sub.marker.map(xof).to_numpy(float)
            y = -np.log10(sub['P'].to_numpy(float))
            hit = sub.significant.to_numpy()
            ax_a.scatter(x[~hit], y[~hit], s=11, color=ccol[c], alpha=0.85,
                         edgecolors='none', zorder=4, rasterized=True)
            # Shape and size carry "significant", colour carries the cohort, so
            # the two encodings never compete for the same channel.
            ax_a.scatter(x[hit], y[hit], s=STATE['significant']['size'], marker='D',
                         color=ccol[c], edgecolors='white', linewidths=0.5, zorder=6)
        S.tier_lines(ax_a, alpha, NOMINAL, label=True)
        ax_a.set_xlim(-0.5, len(xof) - 0.5)
        ax_a.set_ylabel(S.NEGLOG10P)
        ax_a.set_xticks([(lo + hi) / 2 for _g, lo, hi in blocks])
        ax_a.set_xticklabels([f'HLA-{g}' for g, _lo, _hi in blocks], fontstyle='italic')
        # Sixteen gene names on one axis at COL_DOUBLE do not fit: the genes with
        # few markers occupy a few pixels each and their labels overprint. Drop the
        # labels whose rendered boxes collide with a neighbour that was kept -- the
        # tick stays, so the block is still delimited, and the reader identifies the
        # gene from the block, not from a label stack.
        S.thin_tick_labels(ax_a, axis='x')
        ax_a.set_xlabel(f'{MC} markers, ordered by gene then position')
        S.despine(ax_a, grid_axis='y')
        # The markers to NAME above the data: the strongest, capped, because a
        # full HLA hit list does not fit a 7.2 in axis at a readable size.
        for mk in list(best.index[:MAX_AXIS_LABELS]):
            sub = long[long.marker == mk]
            label_rows.append((xof[mk], float(-np.log10(sub['P'].min())), mk))

    # ── (b) effect between every pair of cohorts ───────────────────────────
    rvals = {}
    if not pairs:
        empty_panel(ax_b[0], 'one cohort only —\nno pair to correlate')
    for k, (i, j) in enumerate(pairs):
        ax = ax_b[k]
        ci, cj = cohorts[i], cohorts[j]
        a = long[long.cohort == ci].set_index('marker')
        b = long[long.cohort == cj].set_index('marker')
        idx = a.index.intersection(b.index)
        x = a.loc[idx, 'log_or'].to_numpy(float) if len(idx) else np.array([])
        y = b.loc[idx, 'log_or'].to_numpy(float) if len(idx) else np.array([])
        good = np.isfinite(x) & np.isfinite(y) if len(x) else np.array([], dtype=bool)
        if int(good.sum()) < 2:
            empty_panel(ax, 'no marker fitted\nin both cohorts')
            ax.set_xlabel(SHORT[ci]); ax.set_ylabel(SHORT[cj])
            continue
        hit = ((a.loc[idx, 'significant'].to_numpy()) | (b.loc[idx, 'significant'].to_numpy()))
        # Coloured by the LARGER of the two cohorts, so the ramp still reads as
        # the nesting rather than as an arbitrary choice between the pair.
        ax.scatter(x[~hit], y[~hit], s=11, color=ccol[cj], alpha=0.7,
                   edgecolors='none', zorder=4, rasterized=True)
        ax.scatter(x[hit], y[hit], s=26, color=S.ACCENT, edgecolors='white',
                   linewidths=0.4, zorder=6)
        r = (float(np.corrcoef(x[good], y[good])[0, 1])
             if np.std(x[good]) > 0 and np.std(y[good]) > 0 else float('nan'))
        rvals[(ci, cj)] = r
        e = float(max(np.max(np.abs(x[good])), np.max(np.abs(y[good])), 0.1)) * 1.12
        ax.plot([-e, e], [-e, e], color=S.REFERENCE, lw=1.0, zorder=2)
        ax.axhline(0, color=S.NEUTRAL_D, lw=0.7, ls=':', zorder=1)
        ax.axvline(0, color=S.NEUTRAL_D, lw=0.7, ls=':', zorder=1)
        ax.set_xlim(-e, e); ax.set_ylim(-e, e)
        ax.set_box_aspect(1)
        # The cohort names go on the AXES, not in a centred title: at a third of
        # COL_DOUBLE a title like 'intermediate - full' collides with the panel
        # letter, which panel_tag documents as the limit.
        ax.set_xlabel(f'log OR, {SHORT[ci]}')
        ax.set_ylabel(f'log OR, {SHORT[cj]}')
        if np.isfinite(r):
            ax.annotate(f'$r$ = {r:.3f}', xy=(0.04, 0.94), xycoords='axes fraction',
                        ha='left', va='top', color=S.INK_SOFT,
                        fontsize=plt.rcParams['legend.fontsize'])
        S.despine(ax)

    # ── (c) forest: every marker significant in ANY cohort, in all of them ──
    states = []
    if not rows:
        empty_panel(ax_c, f'no {MC} marker clears $P$ < {S.mathsci(alpha)} in any cohort')
    else:
        yv = np.arange(len(rows))[::-1]
        iv = []
        for y, (mk, c) in zip(yv, rows):
            r = long[(long.marker == mk) & (long.cohort == c)]
            if not len(r):
                states.append('no_fit')
                ax_c.annotate('no usable fit in this cohort', xy=(0.02, y),
                              xycoords=('axes fraction', 'data'), ha='left', va='center',
                              color=S.INK_SOFT,
                              fontsize=plt.rcParams['legend.fontsize'] - 1)
                continue
            r = r.iloc[0]
            st = STATE['significant' if bool(r['significant']) else 'not_significant']
            states.append('significant' if bool(r['significant']) else 'not_significant')
            orv, l95, u95 = (float(r[k]) if pd.notna(r[k]) else np.nan
                             for k in ('OR', 'L95', 'U95'))
            if np.isfinite(l95) and np.isfinite(u95) and l95 > 0:
                ax_c.plot([l95, u95], [y, y], color=ccol[c], lw=2.1,
                          solid_capstyle='round', zorder=3)
                iv += [l95, u95]
            if np.isfinite(orv) and orv > 0:
                ax_c.scatter(orv, y, s=st['size'], marker=st['marker'],
                             facecolor=ccol[c] if st['filled'] else 'white',
                             edgecolors=ccol[c], linewidths=1.3, zorder=5)
                iv.append(orv)
        ax_c.axvline(1.0, color=S.REFERENCE, lw=1.0, zorder=2)
        base = plt.rcParams['ytick.labelsize']
        groups, seen = [], None
        for mk, c in rows:
            if seen != mk:
                g = long.loc[long.marker == mk, 'gene']
                gname = str(g.iloc[0]) if len(g) else ''
                groups.append(([
                    (mk, dict(fontstyle='normal', fontweight='bold', fontsize=base)),
                    (f'HLA-{gname}' if gname and gname != 'nan' else '',
                     dict(fontsize=base - 1.0, color=S.INK_SOFT, fontstyle='italic')),
                ], []))
                seen = mk
            groups[-1][1].append(SHORT[c])
        S.grouped_row_labels(ax_c, groups)
        ax_c.set_xscale('log')
        iv = [v for v in iv if np.isfinite(v) and v > 0]
        ax_c.set_xlim(min(iv) * 0.85 if iv else 0.5, (max(iv) if iv else 2.0) * 3.0)
        S.or_log_axis(ax_c)
        ax_c.set_xlabel(f'{S.OR_SYM} (95% CI) per copy of the {MC}, log scale — {MODEL}')
        for y, (mk, c) in zip(yv, rows):
            v = long.loc[(long.marker == mk) & (long.cohort == c), 'P']
            ax_c.annotate(f'$P$ = {S.p_tex(float(v.min())) if len(v) else "—"}',
                          xy=(0.995, y), xycoords=('axes fraction', 'data'),
                          ha='right', va='center', color=S.INK_SOFT,
                          fontsize=plt.rcParams['legend.fontsize'] - 0.5)
        keys = [k for k in STATE if k in set(states)]
        S.legend_above(ax_c, [Line2D([], [], ls='none', marker=STATE[k]['marker'],
                                     markersize=6.5, markeredgewidth=1.3,
                                     markerfacecolor=S.NEUTRAL_D if STATE[k]['filled']
                                     else 'white', markeredgecolor=S.NEUTRAL_D,
                                     label=STATE[k]['label']) for k in keys]
                       + [Line2D([], [], ls='-', lw=2.0, color=ccol[c], label=SHORT[c])
                          for c in cohorts],
                       ncol=min(4, len(keys) + len(cohorts)), pad_in=0.20)
        S.despine(ax_c, grid_axis='x')

    # ── (d) shared and cohort-specific hits ────────────────────────────────
    own = {}
    for c in cohorts:
        s = set(long.loc[(long.cohort == c) & long.significant, 'marker']) if len(long) else set()
        others = set()
        for c2 in cohorts:
            if c2 != c and len(long):
                others |= set(long.loc[(long.cohort == c2) & long.significant, 'marker'])
        own[c] = (s, s - others)
    n_all = len([mk for mk in sig_markers
                 if int(n_sig_by_marker.get(mk, 0)) == int(n_fit_by_marker.get(mk, 0))
                 and int(n_fit_by_marker.get(mk, 0)) == len(cohorts)])
    if not sig_markers:
        empty_panel(ax_d, 'no marker is significant in any cohort — '
                          'nothing to call shared or distinct')
    else:
        yv = np.arange(len(cohorts))[::-1]
        for y, c in zip(yv, cohorts):
            tot, only = own[c]
            n_shared = len(tot) - len(only)
            ax_d.barh(y, n_shared, height=0.56, color=ccol[c], zorder=3)
            ax_d.barh(y, len(only), left=n_shared, height=0.56, color='white',
                      edgecolor=ccol[c], linewidth=1.1, hatch='///', zorder=3)
            names = sorted(only)[:MAX_NAMES_IN_D]
            extra = len(only) - len(names)
            txt = ('only here: ' + ', '.join(names) + (f' (+{extra} more)' if extra else '')) \
                if names else 'none unique to this cohort'
            ax_d.annotate(f'{len(tot)}  \u00b7  {txt}', xy=(len(tot), y), xytext=(6, 0),
                          textcoords='offset points', ha='left', va='center',
                          color=S.INK_SOFT, fontsize=plt.rcParams['legend.fontsize'] - 0.5)
        ax_d.set_yticks(yv)
        ax_d.set_yticklabels([SHORT[c] for c in cohorts])
        ax_d.set_ylim(-0.7, len(cohorts) - 0.3)
        hi = max(len(t) for t, _o in own.values()) or 1
        # Room on the right for the names, which are the informative half of this
        # panel: the bar alone only says how many.
        ax_d.set_xlim(0, hi * 3.0)
        ax_d.set_xlabel(f'{MC} markers with $P$ < {S.mathsci(alpha)}')
        ax_d.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        S.legend_above(ax_d, [Patch(facecolor=S.NEUTRAL_D, edgecolor=S.NEUTRAL_D,
                                    label='also significant in another cohort'),
                              Patch(facecolor='white', edgecolor=S.NEUTRAL_D, hatch='///',
                                    label='significant in this cohort only')],
                       ncol=2, pad_in=0.16)
        S.despine(ax_d, grid_axis='x')

    # ── caption ─────────────────────────────────────────────────────────────
    per_cohort = ', '.join(f'{len(own[c][0])} ({SHORT[c]})' for c in cohorts)
    r_txt = ('; '.join(f'{SHORT[a]}–{SHORT[b]} $r$ = {v:.3f}'
                       for (a, b), v in rvals.items() if np.isfinite(v)))
    S.caption_block(
        fig, plot_h=plot_h,
        title=(f'HLA {MC} markers under the {MODEL} model in {len(cohorts)} NESTED sample '
               f'set(s): {len(sig_markers)} marker(s) clear $P$ < {S.mathsci(alpha)} in at '
               f'least one, {n_all} in all of them; per cohort {per_cohort}.'),
        panels=[
            (f'{S.NEGLOG10P} for every {MC} marker, one colour per cohort, on ONE marker axis '
             f'ordered by gene then position; bands separate the genes. The solid line is the '
             f'threshold, the dashed line nominal $P$ = {S.mathsci(NOMINAL)}. Diamonds are '
             f'markers clearing the threshold in that cohort.'),
            (f'log odds ratio for each pair of cohorts over the markers fitted in both; the '
             f'grey line is equality, and points significant in either member are red. '
             + (r_txt or 'One cohort only, so there is no pair to correlate.')),
            (f'{S.OR_SYM} and 95% CI for every marker significant in ANY cohort, drawn in ALL '
             f'of them; shape gives that cohort\'s own call. '
             + (f'The {len(drawn)} strongest of {len(sig_markers)} are drawn; all are in the '
                f'table. ' if capped else '')
             + 'The three rows of a group are the SAME marker in three overlapping samples.'),
            ('Per cohort, how many markers clear the threshold and how many of those are '
             'significant in that cohort and in no other; the cohort-specific ones are named.'),
        ],
        notes=(f'HLA {MC} markers, {MODEL} model, {n_markers} with a usable fit'
               + (f' ({n_dropped} excluded on ERRCODE)' if n_dropped else '')
               + f'. Threshold {S.sci(alpha)}'
               + ('' if args.alpha > 0 else f' = 0.05 / {n_markers}, a Bonferroni over this '
                                            f'marker class alone')
               + '. THE COHORTS ARE NESTED — each is a subset of the next and they share every '
                 'case — so a marker significant in all of them has survived the ancestry '
                 'filter, NOT been replicated in independent data. Agreement between these '
                 'cohorts is expected by construction and carries no independent evidence; '
                 'this figure reports how an estimate MOVES as the filter is relaxed. Full '
                 'explanation: ' + Path(args.out_png).stem + '.md'),
        top_pad=0.40, hspace=0.48, wspace=0.34, left='auto', right=0.955,
        margin_axes=[ax_a, ax_b[0], ax_c, ax_d])
    # AFTER caption_block: gene_labels measures the FINAL axes box to size its
    # strip and its leader arms, and caption_block resizes the figure.
    if label_rows:
        # anno_style='expand' — always vertical. An HLA marker id is 10-14
        # characters, and horizontal text anchored near the right-hand end of the
        # axis runs its own width off the canvas edge; rotated 90 degrees its
        # horizontal footprint is one line height and it cannot.
        S.gene_labels(ax_a, [x for x, _y, _t in label_rows],
                      [y for _x, y, _t in label_rows],
                      [t for _x, _y, t in label_rows],
                      anno_style='expand', color=S.ACCENT, weight='bold')
    if xof:
        # AFTER the strip is reserved, and clear of it: the legend is anchored to
        # the axes, so a legend placed before gene_labels moves down with the
        # shrinking axes and ends up inside the band of labels.
        S.legend_above(ax_a, [Line2D([], [], ls='none', marker='o', markersize=5.5,
                                     color=ccol[c], label=SHORT[c]) for c in cohorts]
                       + [Line2D([], [], ls='none', marker='D', markersize=6,
                                 markerfacecolor=S.NEUTRAL_D, markeredgecolor='white',
                                 color=S.NEUTRAL_D, label='significant')],
                       ncol=len(cohorts) + 1, full_width=True,
                       pad_in=0.16 + (S.strip_pad(ax_a) - 7.0) / 72.0)
    S.panel_tag(ax_a, 'a', pad=S.strip_pad(ax_a))
    S.panel_tag(ax_b[0], 'b')
    S.panel_tag(ax_c, 'c')
    S.panel_tag(ax_d, 'd')
    fig.savefig(args.out_png)
    plt.close(fig)

    figure_doc.write_doc(
        args.out_png,
        title=f'HLA {MC} scan across {len(cohorts)} nested sample sets',
        question=(f'How does the evidence for each HLA {MC} marker move as the ancestry filter '
                  f'is relaxed from the narrowest sample set to the widest?'),
        interpretation=(
            'The cohorts are NESTED: each is a subset of the next, and they share every case '
            'and most controls. Nothing in this figure is replication, and the three-row groups '
            'of panel (c) are the shape a reader is most likely to misread as such. What the '
            'figure does report is the behaviour of one estimate as the sample grows: an '
            'interval that narrows while the point estimate holds is a marker gaining sample '
            'size, and an estimate that drifts toward OR = 1 as the filter is relaxed is what a '
            'structure-driven signal does. The HLA region makes this sharper than elsewhere in '
            'the study, because class I and II allele frequencies differ between Japanese '
            'subpopulations, so the ancestry filter that defines the three sets acts directly on '
            'the markers being tested. A marker significant only in the narrowest set is '
            'therefore the case that most needs explaining, and panel (d) is what makes those '
            'countable instead of something the reader has to find by eye.'),
        panels=[
            ('a', 'Evidence on one marker axis',
             f'-log10 P for every {MC} marker with a usable fit, one series per cohort on ONE '
             f'shared marker axis. The axis is ordered by gene, then by the coordinate '
             f'build_hla_markers.py mints within the gene; alternating bands separate the '
             f'genes and the tick labels name them. HLA markers have no meaningful base-pair '
             f'spacing, so this is an ordered categorical axis and never a genomic distance. '
             f'The solid line is the significance threshold and the dashed line nominal '
             f'P = {NOMINAL}; both are drawn on every -log10 P axis in this study. Diamonds '
             f'mark markers clearing the threshold in that cohort — shape carries the call, '
             f'colour carries the cohort, so the two encodings never compete. At most '
             f'{MAX_AXIS_LABELS} markers are named above the data.'),
            ('b', 'Effect between cohort pairs',
             'log odds ratio in one cohort against another, one sub-panel per pair, over the '
             'markers fitted in both. The grey line is equality; dotted lines mark no effect. '
             'Points significant in either member of the pair are red. The annotated r is the '
             'Pearson correlation over every shared marker. These correlations are NOT '
             'replication statistics: the smaller cohort is a subset of the larger, so a high r '
             'is what the sample overlap alone would produce.'),
            ('c', 'The hits in every cohort',
             'A forest over (marker) x (cohort). Every marker significant in ANY cohort is '
             'drawn in ALL of them, because a nested design always has an estimate in the '
             'larger sets. Three states are distinguished, not two: filled diamond = '
             'significant in that cohort, open circle = fitted but not significant, and a row '
             'with no fit at all is annotated in words. A blank would be indistinguishable from '
             'a missing estimate. The group label names the marker and its gene; the row ticks '
             'name the cohort, so the marker id is written once rather than three times.'),
            ('d', 'Shared and cohort-specific',
             'Per cohort, the number of markers clearing the threshold, split into those also '
             'significant somewhere else (solid) and those significant in this cohort and no '
             'other (hatched), with the cohort-specific markers named to the right. Because the '
             'cohorts are nested, "shared" here means the same samples were counted twice — it '
             'is a bookkeeping statement about the three analyses, not evidence.'),
        ],
        numbers=([('marker class', MC), ('model', args.model),
                  ('cohorts (nesting order)', ', '.join(cohorts)),
                  ('markers with a usable fit', n_markers),
                  ('fits excluded on ERRCODE', n_dropped),
                  ('significance threshold', alpha),
                  ('threshold source',
                   'given with --alpha' if args.alpha > 0 else f'Bonferroni 0.05 / {n_markers}'),
                  ('markers significant in at least one cohort', len(sig_markers)),
                  ('markers significant in every cohort', n_all),
                  ('markers drawn in the forest', len(drawn))]
                 + [(f'significant, {c}', len(own[c][0])) for c in cohorts]
                 + [(f'significant only in {c}', len(own[c][1])) for c in cohorts]
                 + [(f'log-OR correlation {a} vs {b}', v) for (a, b), v in rvals.items()]),
        tables=[('Every significant marker, in every cohort', tab if len(tab) else None, 60)],
        reading=[
            'Read panel (a) down, not across: a marker is at the same x in all three series, so '
            'the vertical spread at one x is what the ancestry filter did to that marker.',
            'Do not read agreement between cohorts as replication anywhere in this figure. They '
            'are nested and share every case.',
            'In (c), read down a group. An interval that narrows while the point estimate holds '
            'is a marker gaining sample size; an estimate drifting toward 1 as the filter opens '
            'is what a structure-driven signal does.',
            'In (d), a hatched segment on the narrowest cohort is the finding that most needs '
            'explaining, because that set is the one most exposed to ancestry confounding.',
            'For every marker in every cohort, not only the significant ones, read the TSV '
            'beside this figure.',
        ],
        limits=[
            'It is not a replication analysis and cannot be one: the cohorts are nested and no '
            'independent HLA-typed cohort is available to this study.',
            'It cannot adjudicate between the cohorts. Which sample set is reported is a design '
            'decision, and this figure is one of its inputs, not its verdict.',
            f'It covers the {MC} class under the {args.model} model only. The other marker class '
            f'and the other model are separate analyses with their own figures and their own '
            f'multiple-testing burdens.',
            'The threshold treats the markers of this class as independent tests. HLA markers '
            'are in strong LD, so it is conservative in the number of tests and says nothing '
            'about the correlation between them.',
            'Panel (a) names at most the strongest few markers. An unnamed point is not a point '
            'the figure is claiming does not exist.',
        ],
        defs=['or', 'errcode', 'model'],
        model=S.FORMULAS['glm'])

    print(f'[plot_hla_cohort_compare] {MC} / {args.model}: {len(cohorts)} cohort(s), '
          f'{n_markers:,} marker(s), alpha {alpha:.3e}; {len(sig_markers)} significant in >=1 '
          f'cohort, {n_all} in all; {len(rows)} forest rows -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in plot_hla_cohort_compare: {e}', file=sys.stderr)
        sys.exit(1)
