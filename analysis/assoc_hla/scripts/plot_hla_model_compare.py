#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : One HLA marker class, both models, all three cohorts — does putting
#           the full GRM in the model change what the scan says?
#             (a) evidence   — -log10 P under each model, one panel per cohort
#             (b) effect     — log(OR) under each model, CI whiskers on the leads
#             (c) calibration— lambda_GC per cohort x model
#             (d) the hits   — OR (95% CI) under BOTH models, marker x cohort
#
#           WHY THIS FIGURE EXISTS. The HLA region is where population structure
#           is most likely to masquerade as association: the locus is the most
#           polymorphic in the genome, its allele frequencies differ sharply
#           between Japanese subpopulations, and the cases and controls of this
#           study were sequenced on different platforms. A fixed-effects logistic
#           model can only remove structure through the covariates it is given;
#           the random-effect model fits the SAME samples with a full genetic
#           relationship matrix, so it absorbs relatedness and fine-scale
#           structure the PCs cannot reach. The pair is therefore a read-out of
#           HOW MUCH OF EACH HIT IS STRUCTURE, and every panel is a different
#           projection of that one question.
#
#           A marker on the identity line in (a) and (b) is a marker the GRM does
#           not touch. A marker that leaves the line lost (or gained) its evidence
#           to relatedness, and it is (d) — the odds ratio and its interval under
#           both fits, side by side — that says whether the ESTIMATE moved or only
#           its standard error did.
#
#           NEITHER MODEL IS RANKED HERE. lambda_GC is a calibration read-out and
#           is never applied as a correction, exactly as in the common-variant
#           component; which model is reported is a design decision taken in
#           METHODS, not a contest won by the smaller lambda.
#
#           THE THRESHOLD IS A BONFERRONI OVER THE MARKER CLASS, not the
#           genome-wide constant: an HLA scan is a few hundred tests over one
#           region, and 5e-8 would be a threshold imported from an experiment
#           that was not run. Allele and residue markers are separate analyses
#           with separate burdens (hla_marker_qc.py writes one --extract list
#           each), so the denominator is this class's own marker count.
# Component: assoc_hla
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402
from matplotlib.lines import Line2D      # noqa: E402
from matplotlib.patches import Patch     # noqa: E402
from scipy.stats import chi2             # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
SHARED = Path(__file__).resolve().parent.parent.parent / '_shared' / 'scripts'
sys.path.insert(0, str(SHARED))
import plot_style as S                   # noqa: E402
import figure_doc                        # noqa: E402

# median of a 1-df chi-square: the denominator of lambda_GC, written out rather
# than called so the figure and the document quote the same constant.
CHI2_MEDIAN = 0.4549364231195729
NOMINAL = 0.05                # the second, dashed line on every -log10 P axis
ROW_IN = 0.185                # inches per forest row, to stay legible at 7.2 in
# A cap on the FIGURE only. Every significant marker is in --out-table; drawing
# 60 of them would produce a 20-inch canvas whose rows are unreadable, and the
# caption states the cap whenever it bites.
MAX_FOREST_MARKERS = 12
# The probe is drawn as the open member of every pair, so the reported model is
# never in competition with it for the eye. Same vocabulary as assoc_saige.
MODEL_STYLE = [dict(marker='o', filled=True), dict(marker='s', filled=False)]
CANON_NUM = ('POS', 'A1_FREQ', 'OBS_CT', 'OR', 'LOG(OR)_SE', 'L95', 'U95', 'Z_STAT', 'P')


def parse_args():
    p = argparse.ArgumentParser(
        description='Fixed vs random (full-GRM) HLA scan, one marker class, every cohort.')
    p.add_argument('--sumstat', action='append', required=True, metavar='COHORT/MODEL=PATH',
                   help='repeatable; KEY is <cohort>/<model>, PATH is that scan\'s sumstats.tsv')
    p.add_argument('--marker-class', required=True, choices=('allele', 'residue'),
                   help='the class these sumstats hold; it is a LABEL, not a filter — '
                        'hla_to_sumstats.py already wrote one file per class')
    p.add_argument('--cohort-order', default='',
                   help='comma-separated; empty means the order the --sumstat keys arrived in')
    p.add_argument('--model-order', default='fixed,random',
                   help='comma-separated; the FIRST is the reference model and is drawn filled')
    p.add_argument('--model-label', action='append', default=[], metavar='MODEL=LABEL',
                   help='repeatable; how a model is named on the figure')
    p.add_argument('--alpha', type=float, default=0.0,
                   help='significance threshold. 0 (the default) means Bonferroni 0.05 / M '
                        'over the M distinct markers with a usable fit in this class — an HLA '
                        'scan is a few hundred tests over one region, so the genome-wide '
                        'constant would be a threshold from an experiment that was not run.')
    p.add_argument('--out-png', required=True)
    p.add_argument('--out-table', required=True,
                   help='one row per (marker, cohort, model)')
    return p.parse_args()


def spec(items):
    """['cohort/model=path', …] -> {(cohort, model): path}, existing non-empty files only."""
    out = {}
    for s in items or []:
        key, _, path = s.partition('=')
        c, _, m = key.partition('/')
        if c and m and path and Path(path).exists() and Path(path).stat().st_size:
            out[(c, m)] = path
    return out


def read_sumstat(path, cohort, model):
    """One sumstats.tsv as a tidy frame, ERRCODE applied.

    The 16-column plink2 schema plus gene/position/cohort/model, written by
    hla_to_sumstats.py with na_rep='NA'. `position` is the RESIDUE position id
    (blank for an allele marker) and `POS` is the minted per-gene coordinate;
    both are kept, and ordering uses POS because it is the only numeric one.
    """
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
        'pos': d['POS'] if 'POS' in d.columns else np.nan,
        'cohort': cohort, 'model': model,
        'OR': d['OR'], 'L95': d.get('L95'), 'U95': d.get('U95'),
        'SE': d.get('LOG(OR)_SE'), 'A1_FREQ': d.get('A1_FREQ'),
        'OBS_CT': d.get('OBS_CT'), 'P': d['P'],
        'ERRCODE': err})
    out['log_or'] = np.log(out['OR'].where(out['OR'] > 0))
    return out[ok].reset_index(drop=True), int((~ok).sum())


def lambda_gc(p):
    """median(chi2) / chi2_{1,0.5}, from the P-values of a usable fit."""
    p = np.asarray(pd.to_numeric(pd.Series(p), errors='coerce').dropna(), dtype=float)
    p = p[(p > 0) & (p <= 1)]
    if len(p) < 2:
        return float('nan')
    return float(chi2.isf(np.median(p), 1) / CHI2_MEDIAN)


def empty_panel(ax, msg):
    """An honest blank: say why there is nothing, rather than draw an empty box."""
    ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha='center', va='center',
            color=S.INK_SOFT, fontsize=plt.rcParams['legend.fontsize'])
    ax.set_xticks([]); ax.set_yticks([])
    ax.grid(False)
    # No frame either: an empty axes box reads as a rendering failure, while a
    # bare sentence reads as the statement it is.
    for sp in ax.spines.values():
        sp.set_visible(False)


def logv(v):
    """log of a column of odds ratios, with non-positive values as NaN rather than -inf.

    An interval bound of 0 (complete separation) would otherwise render as an
    infinite whisker, which matplotlib draws across the whole panel.
    """
    a = np.asarray(v, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.log(np.where(a > 0, a, np.nan))


def ptxt(series):
    """A P-value for an annotation, or an em dash when there is no fit at all."""
    v = pd.to_numeric(pd.Series(series), errors='coerce').dropna()
    return S.p_tex(float(v.min())) if len(v) else '—'


def wide(long, cohort, models, value):
    """Marker x model table of one quantity, keeping only markers fitted under both."""
    sub = long[long.cohort == cohort]
    if not len(sub):
        return pd.DataFrame(columns=models)
    w = sub.pivot_table(index='marker', columns='model', values=value, aggfunc='first')
    for m in models:
        if m not in w.columns:
            w[m] = np.nan
    return w.dropna(subset=list(models))


def main():
    args = parse_args()
    S.setup_style('paper')
    ss = spec(args.sumstat)
    if not ss:
        raise SystemExit('ABORT: no --sumstat file exists and is non-empty')
    labels = dict(s.split('=', 1) for s in args.model_label if '=' in s)

    models = [m for m in args.model_order.split(',') if m] or sorted({m for _, m in ss})
    models = [m for m in models if any((c, m) in ss for c, _ in ss)]
    cohorts = [c for c in args.cohort_order.split(',') if c] or \
              list(dict.fromkeys(c for c, _ in ss))
    cohorts = [c for c in cohorts if any((c, m) in ss for m in models)]
    if not cohorts or not models:
        raise SystemExit('ABORT: --cohort-order / --model-order select nothing that was given')
    SHORT = S.shorten(cohorts)
    ccol = {c: S.COHORT_RAMP[i % len(S.COHORT_RAMP)] for i, c in enumerate(cohorts)}
    MC = args.marker_class
    ref = models[0]
    probe = models[1] if len(models) > 1 else None
    lab = {m: labels.get(m, m) for m in models}

    # ── read every scan once ────────────────────────────────────────────────
    frames, n_dropped = [], 0
    for c in cohorts:
        for m in models:
            path = ss.get((c, m))
            if not path:
                continue
            d, bad = read_sumstat(path, c, m)
            n_dropped += bad
            frames.append(d)
    long = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
        columns=['marker', 'gene', 'position', 'pos', 'cohort', 'model', 'OR', 'L95',
                 'U95', 'SE', 'A1_FREQ', 'OBS_CT', 'P', 'ERRCODE', 'log_or'])
    n_markers = int(long['marker'].nunique())
    # The multiple-testing burden is this marker class's own, counted over the
    # markers that actually produced a fit — not over what was submitted.
    alpha = args.alpha if args.alpha > 0 else (0.05 / max(n_markers, 1))
    long['significant'] = long['P'] < alpha

    lam = {(c, m): lambda_gc(long.loc[(long.cohort == c) & (long.model == m), 'P'])
           for c in cohorts for m in models}

    # ── which markers the forest draws ──────────────────────────────────────
    # Significant under EITHER model, in any cohort — the union, because a marker
    # that only one of the two fits calls is the interesting case, not a marker to
    # be dropped for disagreeing.
    sig = long[long.significant]
    best = (sig.groupby('marker')['P'].min().sort_values()
            if len(sig) else pd.Series(dtype=float))
    n_sig_markers = int(len(best))
    drawn = list(best.index[:MAX_FOREST_MARKERS])
    capped = n_sig_markers > len(drawn)

    # ── the table: one row per (marker, cohort, model) ──────────────────────
    # Every significant marker, in every cohort, under every model — including the
    # cohorts and models that did NOT call it, because "not significant here" is a
    # number, and a missing row is indistinguishable from a missing fit.
    keep = list(best.index)
    if not keep and len(long):
        # Nothing cleared the threshold. Report the ten strongest anyway, flagged
        # significant = False, so the file says what the scan's best evidence was
        # rather than being an empty file that could equally mean a crash.
        keep = list(long.groupby('marker')['P'].min().sort_values().index[:10])
    tab_cols = ['marker', 'gene', 'position', 'cohort', 'model', 'OR', 'L95', 'U95',
                'log_or', 'SE', 'P', 'A1_FREQ', 'OBS_CT', 'significant']
    tab = (long[long.marker.isin(keep)]
           .assign(_r=lambda d: d.marker.map({m: i for i, m in enumerate(keep)}),
                   _c=lambda d: d.cohort.map({c: i for i, c in enumerate(cohorts)}),
                   _m=lambda d: d.model.map({m: i for i, m in enumerate(models)}))
           .sort_values(['_r', '_c', '_m'])[tab_cols]
           if len(long) else pd.DataFrame(columns=tab_cols))
    tab.to_csv(args.out_table, sep='\t', index=False, na_rep='NA')

    # ── the forest's rows, counted before the canvas is sized ───────────────
    rows = [(mk, c) for mk in drawn for c in cohorts
            if len(long[(long.marker == mk) & (long.cohort == c)])]
    h_a = 1.95                                     # the per-cohort evidence row
    h_bc = 2.35                                    # effect concordance + lambda
    h_d = max(0.85, len(rows) * ROW_IN) + 0.70     # + the legend strip above it
    plot_h = h_a + h_bc + h_d + 1.85               # + inter-panel gaps and x labels

    fig = plt.figure(figsize=(S.COL_DOUBLE, plot_h))
    gs = fig.add_gridspec(3, 1, height_ratios=[h_a, h_bc, h_d])
    gsa = gs[0].subgridspec(1, len(cohorts), wspace=0.34)
    ax_a = [fig.add_subplot(gsa[0, i]) for i in range(len(cohorts))]
    gsb = gs[1].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.34)
    ax_b = fig.add_subplot(gsb[0, 0])
    ax_c = fig.add_subplot(gsb[0, 1])
    ax_d = fig.add_subplot(gs[2])

    thr = -np.log10(alpha)

    # ── (a) evidence, model against model, one panel per cohort ────────────
    lim, drawn_a = 0.0, []
    for i, c in enumerate(cohorts):
        ax = ax_a[i]
        w = wide(long, c, [ref, probe], 'P') if probe is not None else pd.DataFrame()
        if probe is None:
            empty_panel(ax, f'only one model\n({lab[ref]}) supplied')
        elif not len(w):
            empty_panel(ax, 'no marker fitted\nunder both models')
        else:
            x = -np.log10(w[ref].to_numpy(float))
            y = -np.log10(w[probe].to_numpy(float))
            hit = (w[ref].to_numpy(float) < alpha) | (w[probe].to_numpy(float) < alpha)
            ax.scatter(x[~hit], y[~hit], s=13, color=ccol[c], alpha=0.75,
                       edgecolors='none', zorder=4, rasterized=True)
            ax.scatter(x[hit], y[hit], s=26, color=S.ACCENT, alpha=0.95,
                       edgecolors='white', linewidths=0.4, zorder=6)
            lim = max(lim, float(np.nanmax(x)), float(np.nanmax(y)), thr)
            drawn_a.append(i)
        S.panel_tag(ax, 'a' if i == 0 else '', title=SHORT[c],
                    tsize=plt.rcParams['legend.fontsize'])
    # ONE limit over every cohort, applied after all of them are known: a
    # per-cohort limit would make a taller point mean a rescaled axis rather
    # than stronger evidence.
    lim = (lim if np.isfinite(lim) and lim > 0 else 1.0) * 1.10
    for i in drawn_a:
        ax = ax_a[i]
        ax.plot([0, lim], [0, lim], color=S.REFERENCE, lw=1.0, zorder=2)
        for f in (ax.axhline, ax.axvline):
            f(thr, color=S.ACCENT, lw=0.9, ls='--', zorder=3)
        ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        ax.set_box_aspect(1)
        S.despine(ax)
    if drawn_a and probe is not None:
        ax_a[drawn_a[0]].set_ylabel(f'{S.NEGLOG10P}, {lab[probe]}')
        ax_a[drawn_a[len(drawn_a) // 2]].set_xlabel(f'{S.NEGLOG10P}, {lab[ref]}')

    # ── (b) effect concordance, with intervals on the hits ─────────────────
    r_pooled, n_pairs = float('nan'), 0
    if probe is None:
        empty_panel(ax_b, f'only one model ({lab[ref]}) supplied —\nnothing to compare')
    else:
        xs_all, ys_all = [], []
        for c in cohorts:
            w = wide(long, c, [ref, probe], 'log_or')
            if not len(w):
                continue
            wp = wide(long, c, [ref, probe], 'P')
            hit = wp.reindex(w.index).lt(alpha).any(axis=1).fillna(False).to_numpy()
            x, y = w[ref].to_numpy(float), w[probe].to_numpy(float)
            xs_all.append(x); ys_all.append(y)
            n_pairs += len(x)
            ax_b.scatter(x[~hit], y[~hit], s=13, color=ccol[c], alpha=0.7,
                         edgecolors='none', zorder=4, rasterized=True)
            # Whiskers ON THE LEADS ONLY. Every marker would draw a few hundred
            # crossed intervals over a cloud whose shape is the message.
            if hit.any():
                lo_w = wide(long, c, [ref, probe], 'L95').reindex(w.index)
                hi_w = wide(long, c, [ref, probe], 'U95').reindex(w.index)
                lx, ux = logv(lo_w[ref]), logv(hi_w[ref])
                ly, uy = logv(lo_w[probe]), logv(hi_w[probe])
                ax_b.errorbar(x[hit], y[hit],
                              xerr=np.vstack([np.clip(x - lx, 0, None)[hit],
                                              np.clip(ux - x, 0, None)[hit]]),
                              yerr=np.vstack([np.clip(y - ly, 0, None)[hit],
                                              np.clip(uy - y, 0, None)[hit]]),
                              fmt='none', ecolor=S.NEUTRAL_D, elinewidth=0.7,
                              capsize=0, zorder=5)
                ax_b.scatter(x[hit], y[hit], s=30, color=S.ACCENT, edgecolors='white',
                             linewidths=0.4, zorder=7)
        X = np.concatenate(xs_all) if xs_all else np.array([])
        Y = np.concatenate(ys_all) if ys_all else np.array([])
        good = np.isfinite(X) & np.isfinite(Y) if len(X) else np.array([], dtype=bool)
        if int(good.sum()) >= 2:
            if good.sum() > 2 and np.std(X[good]) > 0 and np.std(Y[good]) > 0:
                r_pooled = float(np.corrcoef(X[good], Y[good])[0, 1])
            e = float(max(np.max(np.abs(X[good])), np.max(np.abs(Y[good])), 0.1)) * 1.12
            ax_b.plot([-e, e], [-e, e], color=S.REFERENCE, lw=1.0, zorder=2)
            ax_b.axhline(0, color=S.NEUTRAL_D, lw=0.7, ls=':', zorder=1)
            ax_b.axvline(0, color=S.NEUTRAL_D, lw=0.7, ls=':', zorder=1)
            ax_b.set_xlim(-e, e); ax_b.set_ylim(-e, e)
            ax_b.set_box_aspect(1)
            ax_b.set_xlabel(f'log OR, {lab[ref]}')
            ax_b.set_ylabel(f'log OR, {lab[probe]}')
            if np.isfinite(r_pooled):
                ax_b.annotate(f'$r$ = {r_pooled:.3f}', xy=(0.04, 0.94),
                              xycoords='axes fraction', ha='left', va='top',
                              fontsize=plt.rcParams['legend.fontsize'], color=S.INK_SOFT)
            S.despine(ax_b)
        else:
            empty_panel(ax_b, 'no marker fitted under both models')

    # ── (c) lambda_GC, cohort x model ──────────────────────────────────────
    xs = np.arange(len(cohorts), dtype=float)
    width = 0.72 / max(len(models), 1)
    finite = [v for v in lam.values() if np.isfinite(v)]
    if finite:
        for j, m in enumerate(models):
            st = MODEL_STYLE[j % len(MODEL_STYLE)]
            v = [lam.get((c, m), np.nan) for c in cohorts]
            off = (j - (len(models) - 1) / 2.0) * width
            ax_c.bar(xs + off, [x - 1.0 if np.isfinite(x) else np.nan for x in v],
                     bottom=1.0, width=width * 0.88,
                     color=[ccol[c] for c in cohorts] if st['filled'] else 'white',
                     edgecolor=[ccol[c] for c in cohorts],
                     linewidth=1.0, hatch=None if st['filled'] else '///', zorder=3)
        # 1.0 is the null, so the bars grow FROM it in both directions — a bar
        # from zero would put every lambda at visually the same height.
        ax_c.axhline(1.0, color=S.ACCENT, lw=1.1, zorder=4)
        lo, hi = min(finite + [1.0]), max(finite + [1.0])
        pad = max(0.03, (hi - lo) * 0.20)
        ax_c.set_ylim(lo - pad, hi + pad)
        # The value labels go on LAST, after the limit is set: they grow the axis
        # to make room for themselves, and a later set_ylim would undo that.
        for j, m in enumerate(models):
            off = (j - (len(models) - 1) / 2.0) * width
            for c, x in zip(cohorts, xs):
                val = lam.get((c, m), np.nan)
                if not np.isfinite(val):
                    continue
                # Outside the bar in the direction the bar points, so a deflated
                # lambda's label is not printed inside its own bar.
                S.value_labels(ax_c, [x + off], [val], [f'{val:.3f}'],
                               offset=4.0 if val >= 1.0 else -4.0,
                               va='bottom' if val >= 1.0 else 'top',
                               fontsize=plt.rcParams['legend.fontsize'] - 1.5)
        ax_c.set_xticks(xs)
        ax_c.set_xticklabels([SHORT[c] for c in cohorts])
        ax_c.set_xlim(-0.6, len(cohorts) - 0.4)
        ax_c.set_ylabel(S.LAMBDA_GC)
        S.despine(ax_c, grid_axis='y')
        S.legend_above(ax_c, [Patch(facecolor=S.NEUTRAL_D if MODEL_STYLE[j % 2]['filled']
                                    else 'white', edgecolor=S.NEUTRAL_D,
                                    hatch=None if MODEL_STYLE[j % 2]['filled'] else '///',
                                    label=lab[m]) for j, m in enumerate(models)],
                       ncol=len(models), pad_in=0.16)
    else:
        empty_panel(ax_c, 'no usable fit to compute\n' + S.LAMBDA_GC + ' from')

    # ── (d) the hits, both models, one row per marker x cohort ─────────────
    if not rows:
        empty_panel(ax_d, f'no {MC} marker clears $P$ < {S.mathsci(alpha)} '
                          f'under either model')
    else:
        yv = np.arange(len(rows))[::-1]
        iv = []
        for y, (mk, c) in zip(yv, rows):
            for j, m in enumerate(models):
                r = long[(long.marker == mk) & (long.cohort == c) & (long.model == m)]
                st = MODEL_STYLE[j % len(MODEL_STYLE)]
                dy = (j - (len(models) - 1) / 2.0) * 0.34
                if not len(r):
                    ax_d.annotate('not fitted', xy=(0.02, y + dy),
                                  xycoords=('axes fraction', 'data'), ha='left',
                                  va='center', color=S.INK_SOFT,
                                  fontsize=plt.rcParams['legend.fontsize'] - 1)
                    continue
                r = r.iloc[0]
                orv, l95, u95 = (float(r[k]) if pd.notna(r[k]) else np.nan
                                 for k in ('OR', 'L95', 'U95'))
                if np.isfinite(l95) and np.isfinite(u95) and l95 > 0:
                    ax_d.plot([l95, u95], [y + dy, y + dy], color=ccol[c], lw=1.9,
                              solid_capstyle='round', zorder=3)
                    iv += [l95, u95]
                if np.isfinite(orv) and orv > 0:
                    ax_d.scatter(orv, y + dy, s=30, marker=st['marker'],
                                 facecolor=ccol[c] if st['filled'] else 'white',
                                 edgecolors=ccol[c], linewidths=1.2, zorder=5)
                    iv.append(orv)
        ax_d.axvline(1.0, color=S.REFERENCE, lw=1.0, zorder=2)
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
        S.grouped_row_labels(ax_d, groups)
        ax_d.set_xscale('log')
        iv = [v for v in iv if np.isfinite(v) and v > 0]
        # A reserved right-hand column for the two P-values, so they cannot land
        # on the widest interval whatever it turns out to be.
        ax_d.set_xlim(min(iv) * 0.85 if iv else 0.5, (max(iv) if iv else 2.0) * 3.0)
        S.or_log_axis(ax_d)
        ax_d.set_xlabel(f'{S.OR_SYM} (95% CI) per copy of the {MC}, log scale')
        for y, (mk, c) in zip(yv, rows):
            txt = ' / '.join(
                ptxt(long.loc[(long.marker == mk) & (long.cohort == c)
                              & (long.model == m), 'P'])
                for m in models)
            ax_d.annotate(txt, xy=(0.995, y), xycoords=('axes fraction', 'data'),
                          ha='right', va='center', color=S.INK_SOFT,
                          fontsize=plt.rcParams['legend.fontsize'] - 1)
        S.legend_above(ax_d, [Line2D([], [], ls='none',
                                     marker=MODEL_STYLE[j % len(MODEL_STYLE)]['marker'],
                                     markersize=6.5, markeredgewidth=1.2,
                                     markeredgecolor=S.NEUTRAL_D,
                                     markerfacecolor=S.NEUTRAL_D
                                     if MODEL_STYLE[j % len(MODEL_STYLE)]['filled'] else 'white',
                                     label=lab[m]) for j, m in enumerate(models)]
                       + [Line2D([], [], ls='-', lw=2.0, color=ccol[c], label=SHORT[c])
                          for c in cohorts],
                       ncol=min(5, len(models) + len(cohorts)), pad_in=0.20)
        S.despine(ax_d, grid_axis='x')

    # ── caption ─────────────────────────────────────────────────────────────
    lam_txt = '; '.join(
        f"{SHORT[c]} " + ' vs '.join(f'{lam.get((c, m), float("nan")):.3f}' for m in models)
        for c in cohorts)
    conc = (f' Pooled log-{S.OR_SYM} correlation {r_pooled:.3f}.'
            if np.isfinite(r_pooled) else '')
    S.caption_block(
        fig, plot_h=plot_h,
        title=(f'HLA {MC} scan under a fixed-effects and a full-GRM model, {len(cohorts)} '
               f'nested cohort(s): {S.LAMBDA_GC} ' + lam_txt +
               f'; {n_sig_markers} marker(s) clear $P$ < {S.mathsci(alpha)} under either '
               f'model.' + conc),
        panels=[
            (f'{S.NEGLOG10P} under {lab[ref]} (x) against '
             f'{lab.get(probe, "the second model")} (y), one point per {MC} marker, one '
             f'sub-panel per cohort on one shared limit. The grey line is equality and the '
             f'dashed red lines are the threshold on each axis; markers clearing it under '
             f'either model are red.'),
            (f'log {S.OR_SYM} under the two models, all cohorts, with 95% CI whiskers on the '
             f'markers that clear the threshold. Distance from the grey line is what the GRM '
             f'changes about the ESTIMATE.'),
            (f'{S.LAMBDA_GC} for each cohort under each model, drawn as a deviation from 1 '
             f'(red line). Filled = {lab[ref]}, hatched = {lab.get(probe, "second model")}.'),
            (f'{S.OR_SYM} and 95% CI for every marker significant under EITHER model, in '
             f'every cohort, under BOTH models. Right-hand column: '
             + ' / '.join(lab[m] for m in models) + ' $P$.'
             + (f' The {len(drawn)} strongest of {n_sig_markers} significant markers are '
                f'drawn; all of them are in the table.' if capped else '')),
        ],
        notes=(f'HLA {MC} markers, {n_markers} with a usable fit'
               + (f' ({n_dropped} fit(s) excluded on ERRCODE)' if n_dropped else '')
               + f'. Threshold {S.sci(alpha)}'
               + ('' if args.alpha > 0 else f' = 0.05 / {n_markers}, a Bonferroni over this '
                                            f'marker class alone')
               + '. The two models fit the SAME samples and differ only in the random effect, '
                 'so agreement between them is not replication and carries no independent '
                 'evidence. ' + S.LAMBDA_GC + ' is reported, never applied. Full explanation: '
               + Path(args.out_png).stem + '.md'),
        top_pad=0.40, hspace=0.55, wspace=0.34, left='auto', right=0.955,
        margin_axes=[ax_a[0], ax_b, ax_c, ax_d])
    for ax, letter in ((ax_b, 'b'), (ax_c, 'c'), (ax_d, 'd')):
        S.panel_tag(ax, letter)
    fig.savefig(args.out_png)
    plt.close(fig)

    # ── the document beside it ──────────────────────────────────────────────
    figure_doc.write_doc(
        args.out_png,
        title=f'HLA {MC} scan: fixed-effects against full-GRM model',
        question=(f'Does putting a full genetic relationship matrix in the model change which '
                  f'HLA {MC} markers are associated, or by how much — and is either fit '
                  f'calibrated?'),
        interpretation=(
            'The two models fit the SAME samples with the SAME phenotype and covariates and '
            'differ only in whether a full genetic relationship matrix is included as a random '
            'effect. The difference between them is therefore what relatedness and fine-scale '
            'structure were contributing to the fixed-effects result, and nothing else. The HLA '
            'region is where that matters most in this study: it is the most polymorphic region '
            'of the genome, its allele frequencies differ between Japanese subpopulations, and '
            'the cases and the controls of this study were sequenced on different platforms. A '
            'marker sitting on the identity line in panels (a) and (b) is one the GRM does not '
            'touch. A marker that loses evidence under the GRM was, in part, being carried by '
            'relatedness; panel (d) separates the two ways that can happen, because an odds '
            'ratio that holds while its interval widens is a standard-error effect, while an '
            'odds ratio that moves toward 1 is the estimate itself deflating. Neither model is '
            'declared correct here and lambda_GC is never applied as a correction — it is a '
            'calibration read-out, as everywhere else in this study.'),
        panels=[
            ('a', 'Evidence, model against model',
             f'-log10 P under {lab[ref]} on x and {lab.get(probe, "the second model")} on y, '
             f'one point per {MC} marker, one sub-panel per cohort. All sub-panels share one '
             f'limit, so a marker can be compared between cohorts as well as between models. '
             f'The grey line is equality; the dashed red lines mark the significance threshold '
             f'on each axis, so whether a marker crosses under one model and not the other is '
             f'readable directly. Points clearing the threshold under either model are drawn in '
             f'red at a larger size.'),
            ('b', 'Effect concordance',
             'log odds ratio under the two models, pooled over cohorts and coloured by cohort, '
             'with 95% confidence whiskers drawn on the significant markers only — every marker '
             'would put a few hundred crossed intervals over a cloud whose shape is the whole '
             'message. Dotted lines mark log OR = 0 (no effect) on both axes. The annotated r is '
             'the Pearson correlation over every marker fitted under both models.'),
            ('c', 'Calibration',
             'lambda_GC = median(chi2)/0.4549 for each cohort under each model, computed from '
             'the P-values of markers with a usable fit (ERRCODE = "."). The bars are drawn as '
             'deviations from 1 rather than from 0, because a bar from zero renders every '
             'plausible lambda at visually the same height. An HLA scan is a few hundred tests '
             'over one region under strong LD, so this lambda is a much noisier and more '
             'structure-dependent quantity than a genome-wide one; it is read as a coarse '
             'read-out of whether one model is systematically shifted relative to the other, '
             'not as an inflation estimate.'),
            ('d', 'The hits under both models',
             'Every marker that clears the threshold under EITHER model, in every cohort, with '
             'its odds ratio and 95% CI under BOTH. One row group per marker; within a group '
             'one row per cohort, and within a row the two models are offset from each other, '
             'so the pair that must be compared is adjacent by construction. Colour is the '
             'cohort; the filled marker is ' + lab[ref] + '. The right-hand column carries both '
             'P-values in the order ' + ' / '.join(lab[m] for m in models) + '.'),
        ],
        numbers=([('marker class', MC),
                  ('cohorts', ', '.join(cohorts)),
                  ('models', ', '.join(f'{m} ({lab[m]})' for m in models)),
                  ('markers with a usable fit', n_markers),
                  ('fits excluded on ERRCODE', n_dropped),
                  ('significance threshold', alpha),
                  ('threshold source',
                   'given with --alpha' if args.alpha > 0 else f'Bonferroni 0.05 / {n_markers}'),
                  ('markers significant under either model', n_sig_markers),
                  ('markers drawn in the forest', len(drawn)),
                  ('forest rows', len(rows)),
                  ('pooled log-OR correlation', r_pooled)]
                 + [(f'lambda_GC {c} / {m}', lam.get((c, m), float('nan')))
                    for c in cohorts for m in models]
                 + [(f'markers significant, {c} / {m}',
                     int(((long.cohort == c) & (long.model == m) & long.significant).sum()))
                    for c in cohorts for m in models]),
        tables=[('Every significant marker, in every cohort, under both models',
                 tab if len(tab) else None, 60)],
        reading=[
            'Read panel (a) as one cloud per cohort against the identity line, not as three '
            'separate scans: the limits are shared.',
            'A red point off the identity line in (a) is a marker whose EVIDENCE depends on the '
            'model. Go to (d) and read its two intervals to see whether the estimate moved or '
            'only the standard error did.',
            'In (c), read the pair above each cohort, not the series: the distance between the '
            'two bars is what the GRM is worth there.',
            'In (d), read across a row for the model comparison and down a group for the cohort '
            'comparison. The cohorts are nested, so the three rows of a group share samples.',
            'Nothing here ranks the two models. Which one is reported is a design decision, '
            'stated in METHODS, not the one with the smaller lambda.',
        ],
        limits=[
            'It is not replication. The two models fit the same samples, so agreement between '
            'them says the result does not depend on the random effect — not that it is real.',
            'lambda_GC over a few hundred markers in one strongly-LD region is a coarse '
            'statistic. It cannot separate confounding from a genuine concentration of signal, '
            'and in the HLA the latter is expected.',
            f'Only the {MC} class is drawn. The other marker class is a separate analysis with '
            f'its own multiple-testing burden and its own figure.',
            'A marker absent from one model\'s output does not appear in that model\'s panels. '
            'Rows where one model produced no fit are annotated rather than dropped in (d).',
            'The threshold is a Bonferroni over the markers of this class and treats them as '
            'independent tests, which HLA markers are not — they are in strong LD, so this is '
            'conservative in the number of tests and says nothing about the correlation between '
            'them.',
        ],
        defs=['or', 'lambda_gc', 'errcode', 'model'],
        model=S.FORMULAS['lambda'])

    print(f'[plot_hla_model_compare] {MC}: {len(cohorts)} cohort(s) x {len(models)} model(s); '
          f'{n_markers:,} marker(s), alpha {alpha:.3e}; {n_sig_markers} significant '
          f'({len(drawn)} drawn, {len(rows)} forest rows) -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in plot_hla_model_compare: {e}', file=sys.stderr)
        sys.exit(1)
