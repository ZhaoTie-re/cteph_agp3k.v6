#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The reported model against its calibration probe, in every cohort.
#             (a) calibration — lambda_GC for each model in each cohort, with the
#                 QQ curves overlaid per cohort
#             (b) loci        — every peak lead's OR (95% CI) under BOTH models
#
#           The probe fits the SAME samples and the SAME GRM with the PCs
#           REMOVED. Two questions follow, and the two rows answer one each:
#
#           (a) HOW MUCH DO THE PCs STILL DO once a GRM is in the model? If the
#               probe's lambda is already near 1, the GRM has absorbed the
#               structure by itself and the PCs are close to redundant; if the
#               probe is inflated and the reported model is not, the PCs are
#               carrying real work the GRM cannot do.
#
#           (b) IS ANY SIGNAL HELD UP BY THE PCs? A lead whose effect is stable
#               across the two is not a structure artefact. One that only
#               survives with the PCs in the model, or only without them, is a
#               finding about the covariates rather than about the locus.
#
#           lambda is a read-out, never a correction — the same rule the rest of
#           the component follows.
# Component: assoc_saige
# Used by : assoc_saige.nf  process COMPARE_MODELS
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
SHARED = Path(__file__).resolve().parent.parent.parent / '_shared' / 'scripts'
sys.path.insert(0, str(SHARED))
import plot_style as S
import figure_doc

PLOT_H = 5.6
GWT_CAP = r'$5\times10^{-8}$'          # inches of croppable plot block
GW = 'genome_wide'
# The probe is a diagnostic, not a second result, so it is drawn as the open /
# lighter member of every pair and never competes with the reported model.
MODEL_STYLE = [dict(marker='o', filled=True), dict(marker='o', filled=False)]


def parse_args():
    p = argparse.ArgumentParser(description='Reported model vs calibration probe, per cohort.')
    p.add_argument('--sumstat', action='append', required=True, metavar='COHORT/MODEL=PATH')
    p.add_argument('--scan-qc', action='append', required=True, metavar='COHORT/MODEL=PATH')
    p.add_argument('--crosscohort', required=True,
                   help='lead_crosscohort.tsv. Panels (b) and (c) draw EXACTLY the variants '
                        'cohort_compare panels (b) and (c) do, read from its `panel` column, '
                        'so the two figures are about the same loci.')
    p.add_argument('--model-label', action='append', default=[], metavar='MODEL=LABEL')
    p.add_argument('--cohort-order', default='')
    p.add_argument('--model-order', default='')
    p.add_argument('--peak-model', default='main')
    p.add_argument('--alpha', type=float, default=S.GW_ALPHA)
    p.add_argument('--suggestive', type=float, default=S.SUGGEST_ALPHA)
    p.add_argument('--out-png', required=True)
    p.add_argument('--out-table', required=True)
    return p.parse_args()


def spec(items):
    """['cohort/model=path', …] -> {(cohort, model): path} for files that exist."""
    out = {}
    for s in items or []:
        key, _, path = s.partition('=')
        c, _, m = key.partition('/')
        if path and Path(path).exists() and Path(path).stat().st_size:
            out[(c, m)] = path
    return out


def qq_points(p, cap=60000):
    """Observed vs expected -log10 P, thinned for DRAWING only.

    Thinning a QQ is safe here and nowhere else in this component: the curve is
    read as a shape, the whole tail is kept, and only the dense lower arm — where
    thousands of points land on the same pixel — is sampled. lambda is computed
    from every P, not from these points.
    """
    p = np.asarray(p, dtype=float)
    p = p[(p > 0) & (p <= 1)]
    n = len(p)
    obs = -np.log10(np.sort(p))
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    if n > cap:
        keep = np.unique(np.concatenate([np.arange(min(5000, n)),
                                         np.linspace(0, n - 1, cap).astype(int)]))
        obs, exp = obs[keep], exp[keep]
    return exp, obs


def main():
    args = parse_args()
    S.setup_style('paper')
    ss, qcs = spec(args.sumstat), spec(args.scan_qc)
    labels = dict(s.split('=', 1) for s in args.model_label if '=' in s)

    cohorts = [c for c in args.cohort_order.split(',') if c] or \
              sorted({c for c, _ in ss})
    models = [m for m in args.model_order.split(',') if m] or \
             sorted({m for _, m in ss})
    cohorts = [c for c in cohorts if any((c, m) in ss for m in models)]
    SHORT = S.shorten(cohorts)
    ccol = {c: S.COHORT_RAMP[i % len(S.COHORT_RAMP)] for i, c in enumerate(cohorts)}

    # scan_qc.tsv is written once per cohort and carries a row per model, so the
    # same file arrives under every (cohort, model) key; read it once per cohort.
    lam = {}
    for c in cohorts:
        for m in models:
            path = qcs.get((c, m))
            if not path:
                continue
            d = pd.read_csv(path, sep='\t')
            row = d[d['model'] == m]
            if len(row):
                lam[(c, m)] = float(row.iloc[0]['lambda_gc'])

    fig = plt.figure(figsize=(S.COL_DOUBLE, PLOT_H))
    gs = fig.add_gridspec(2, len(cohorts) + 1, height_ratios=[1.0, 1.15],
                          width_ratios=[1.25] + [1.0] * len(cohorts))
    ax_lam = fig.add_subplot(gs[0, 0])
    ax_qq = [fig.add_subplot(gs[0, i + 1]) for i in range(len(cohorts))]
    ax_p = fig.add_subplot(gs[1, :max(1, len(cohorts) // 2 + 1)])
    ax_e = fig.add_subplot(gs[1, max(1, len(cohorts) // 2 + 1):])

    # ── (a) left: lambda for both models, per cohort ───────────────────────
    xs = np.arange(len(cohorts))
    ax_lam.axhline(1.0, color=S.ACCENT, lw=1.1, zorder=2)
    for j, m in enumerate(models):
        v = [lam.get((c, m), np.nan) for c in cohorts]
        st = MODEL_STYLE[j % len(MODEL_STYLE)]
        ax_lam.plot(xs, v, ls='-' if st['filled'] else '--', lw=0.9,
                    color=S.NEUTRAL_D, zorder=2)
        ax_lam.scatter(xs, v, s=46, marker=st['marker'],
                       facecolors=[ccol[c] for c in cohorts] if st['filled'] else 'white',
                       edgecolors=[ccol[c] for c in cohorts], linewidths=1.4, zorder=5)
    finite = [x for x in lam.values() if np.isfinite(x)]
    if finite:
        lo, hi = min(finite + [1.0]), max(finite + [1.0])
        pad = max(0.02, (hi - lo) * 0.35)
        ax_lam.set_ylim(lo - pad, hi + pad)
    ax_lam.set_xticks(xs)
    ax_lam.set_xticklabels([SHORT[c] for c in cohorts],
                           fontsize=plt.rcParams['xtick.labelsize'] - 1)
    ax_lam.set_ylabel(S.LAMBDA_GC)
    ax_lam.set_xlim(-0.55, len(cohorts) - 0.45)
    S.despine(ax_lam, grid_axis='y')

    # ── (a) right: one QQ per cohort, both models overlaid ─────────────────
    qmax = 0.0
    for i, c in enumerate(cohorts):
        ax = ax_qq[i]
        for j, m in enumerate(models):
            path = ss.get((c, m))
            if not path:
                continue
            d = pd.read_csv(path, sep='\t', usecols=['P', 'ERRCODE'],
                            dtype={'ERRCODE': str})
            d = d[d['ERRCODE'].fillna('.').eq('.')]
            exp, obs = qq_points(pd.to_numeric(d['P'], errors='coerce').dropna().to_numpy())
            qmax = max(qmax, float(max(exp.max(), obs.max())))
            st = MODEL_STYLE[j % len(MODEL_STYLE)]
            ax.plot(exp, obs, lw=1.3 if st['filled'] else 1.0,
                    ls='-' if st['filled'] else (0, (3, 2)),
                    color=ccol[c], alpha=1.0 if st['filled'] else 0.75, zorder=4 - j)
        ax.set_title(SHORT[c], fontsize=plt.rcParams['legend.fontsize'], pad=3)
        ax.set_box_aspect(1)
        if i == 0:
            ax.set_ylabel(f'observed {S.NEGLOG10P}')
        S.despine(ax)
    for ax in ax_qq:
        lim = qmax * 1.05 if qmax else 1.0
        ax.plot([0, lim], [0, lim], color=S.ACCENT, lw=1.0, zorder=2)
        ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    if ax_qq:
        ax_qq[len(ax_qq) // 2].set_xlabel(f'expected {S.NEGLOG10P}')

    # ── (b) every peak lead under both models, in every cohort ─────────────
    # BOTH TIERS. A forest could only hold a handful of loci; the scatter below
    # takes the whole peak list, which is the point — whether the PCs matter is a
    # question about the tier as a whole, not about the few strongest loci.
    cc = pd.read_csv(args.crosscohort, sep='\t', dtype={'variant_id': str, 'chrom': str})
    missing = [k for k in ('panel', 'locus_rep') if len(cc) and k not in cc.columns]
    if missing:
        raise SystemExit(f'lead_crosscohort.tsv has no {missing} column(s) — it was written by an '
                         'older cross_cohort.py; re-run CROSS_COHORT')
    # One row per LOCUS, at its representative variant — the same loci
    # cohort_compare draws. The table repeats each variant once per cohort, and
    # a locus can carry more than one lead, so both have to be collapsed here.
    panels = (cc[cc.panel.isin(['b', 'c']) & (cc.locus_rep.astype(int) == 1)]
              .drop_duplicates('variant_id')[['variant_id', 'chrom', 'pos', 'panel']]
              if len(cc) else pd.DataFrame(columns=['variant_id', 'chrom', 'pos', 'panel']))

    rows = []
    for c in cohorts:
        leads = panels
        # Read each model's scan ONCE per cohort, not once per lead: the inner
        # lookup used to re-read a 5 M-row file for every locus.
        stats = {}
        for m in models:
            path = ss.get((c, m))
            if path:
                d = pd.read_csv(path, sep='\t', usecols=['ID', 'OR', 'L95', 'U95', 'P'],
                                dtype={'ID': str})
                stats[m] = d.set_index('ID')
        for _, r in leads.iterrows():
            # No gene label: distance from the identity line is the whole message
            # and names would not survive the density. The per-locus figures and
            # cohort_compare carry them.
            row = {'cohort': c, 'panel': r['panel'],
                   'variant_id': r['variant_id'], 'chrom': r['chrom'], 'pos': r['pos']}
            ok = True
            for m in models:
                d = stats.get(m)
                if d is None or str(r['variant_id']) not in d.index:
                    ok = False
                    break
                h = d.loc[str(r['variant_id'])]
                row[f'OR_{m}'] = float(h['OR'])
                row[f'L95_{m}'] = float(h['L95'])
                row[f'U95_{m}'] = float(h['U95'])
                row[f'P_{m}'] = float(h['P'])
            if ok:
                rows.append(row)
    tab = pd.DataFrame(rows)
    tab.to_csv(args.out_table, sep='\t', index=False)

    ref, probe = (models + models)[0], (models + models)[1] if len(models) > 1 else models[0]
    if len(tab):
        # Two scatters, one per quantity. The identity line is "the covariates
        # changed nothing", so distance from it IS the answer, read over every
        # locus at once rather than over the handful a forest would fit.
        for ax, key, lab, logscale in (
                (ax_p, 'P', S.NEGLOG10P, False),
                (ax_e, 'OR', S.OR_SYM, True)):
            for c in cohorts:
                sub = tab[tab.cohort == c]
                for pan, mk, sz in (('b', 'D', 34), ('c', 'o', 18)):
                    s = sub[sub.panel == pan]
                    if not len(s):
                        continue
                    x = s[f'{key}_{ref}'].to_numpy(float)
                    y = s[f'{key}_{probe}'].to_numpy(float)
                    if key == 'P':
                        x, y = -np.log10(x), -np.log10(y)
                    ax.scatter(x, y, s=sz, marker=mk, c=ccol[c], edgecolors='white',
                               linewidths=0.5, alpha=0.9, zorder=5 if pan == 'b' else 4)
            lo = min(ax.get_xlim()[0], ax.get_ylim()[0])
            hi = max(ax.get_xlim()[1], ax.get_ylim()[1])
            if logscale:
                lo = max(lo, 0.05)
                ax.set_xscale('log'); ax.set_yscale('log')
                S.or_log_axis(ax, axis='x'); S.or_log_axis(ax, axis='y')
                ax.axhline(1.0, color=S.REFERENCE, lw=0.8, ls=':', zorder=2)
                ax.axvline(1.0, color=S.REFERENCE, lw=0.8, ls=':', zorder=2)
            ax.plot([lo, hi], [lo, hi], color=S.ACCENT, lw=1.0, zorder=3)
            ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
            ax.set_box_aspect(1)
            ax.set_xlabel(f'{lab}, {labels.get(ref, ref)}'
                          if len(labels.get(ref, ref)) < 22 else f'{lab}, {ref}')
            ax.set_ylabel(f'{lab}, {labels.get(probe, probe)}'
                          if len(labels.get(probe, probe)) < 22 else f'{lab}, {probe}')
            S.despine(ax)
        # The genome-wide line on the P panel, so "did it cross" stays readable.
        thr = -np.log10(args.alpha)
        for f in (ax_p.axhline, ax_p.axvline):
            f(thr, color=S.REFERENCE, lw=0.8, ls=':', zorder=2)
    else:
        for ax in (ax_p, ax_e):
            ax.text(0.5, 0.5, 'no peak lead to compare', transform=ax.transAxes,
                    ha='center', va='center', color=S.REFERENCE)
            ax.set_xticks([]); ax.set_yticks([])

    handles = [Line2D([], [], ls='none', marker=MODEL_STYLE[j % len(MODEL_STYLE)]['marker'],
                      markersize=7, markeredgewidth=1.3, markeredgecolor=S.NEUTRAL_D,
                      markerfacecolor=S.NEUTRAL_D if MODEL_STYLE[j % len(MODEL_STYLE)]['filled']
                      else 'white',
                      label=labels.get(m, m))
               for j, m in enumerate(models)]
    S.legend_above(ax_lam, handles, ncol=len(handles), y=1.06, full_width=True)

    fig.suptitle('Reported model against its calibration probe',
                 fontsize=plt.rcParams['axes.titlesize'] + 1, fontweight='bold', y=0.988)

    lam_txt = '; '.join(
        f"{SHORT[c]} " + ' vs '.join(f'{lam.get((c, m), float("nan")):.3f}' for m in models)
        for c in cohorts if any((c, m) in lam for m in models))
    S.caption_block(
        fig,
        title=(f'{S.LAMBDA_GC} with and without the principal components, on the same samples and '
               f'the same GRM: {lam_txt}.'),
        panels=[
            (f'{S.LAMBDA_GC} per cohort under each model; the dashed red line is 1. '
             f'Right, the matching quantile–quantile curves, one panel per cohort, both models '
             f'overlaid on one limit.'),
            (f'{S.NEGLOG10P} for every peak lead under each model; the red line is equality and '
             f'dotted lines mark {GWT_CAP}. Diamonds are cohort_compare panel (b) variants, '
             'circles panel (c).'),
            (f'{S.OR_SYM} for the same leads, log axes. Distance from the red line is what the '
             f'principal components change.'),
        ],
        notes=('The probe fits the same samples and the same full GRM with the PCs removed, so the '
               'difference between the pair is what the PCs contribute on top of the GRM. It defines '
               'no peaks and nothing downstream reads it. Full explanation: '
               + Path(args.out_png).stem + '.md'),
        plot_h=PLOT_H, top_pad=0.52, hspace=0.55, wspace=0.34, left='auto', right=0.985,
        margin_axes=[ax_lam, ax_p, ax_e])
    for ax, letter in ((ax_lam, 'a'), (ax_p, 'b'), (ax_e, 'c')):
        S.panel_tag(ax, letter)
    fig.savefig(args.out_png)
    plt.close(fig)

    figure_doc.write_doc(
        args.out_png,
        title='Reported model against its calibration probe',
        question=('How much of the population structure does the GRM absorb on its own, how much do '
                  'the principal components still remove on top of it, and does any signal depend on '
                  'them?'),
        interpretation=(
            'The two models fit the SAME samples and the SAME full genetic relationship matrix and '
            'differ only in whether the principal components are among the covariates. The pair is '
            'therefore a direct read-out of what the PCs contribute once relatedness and broad '
            'structure are already in the model. A probe whose lambda_GC is already near 1 says the '
            'GRM has done the work by itself; a probe that is inflated while the reported model is '
            'not says the PCs are removing structure the GRM cannot reach. In panel (b) a lead whose '
            'interval is stable across the pair is not an artefact of the covariate choice, while '
            'one that moves is a statement about the covariates rather than about the locus. The '
            'probe defines no peaks, receives no fine-mapping, conditional analysis or regional '
            'plot, and must not be reported as a second result.'),
        panels=[
            ('a', 'Calibration',
             'lambda_GC = median(chi2)/0.4549 for each model in each cohort, taken from that '
             "cohort's own scan QC table rather than recomputed. Right: the quantile-quantile curves "
             'behind those numbers, one panel per cohort with both models overlaid and every panel '
             'on one shared limit, so the curves are comparable between cohorts as well as between '
             'models. The curves are thinned for drawing only — the whole tail is kept and lambda '
             'is computed from every P-value.'),
            ('b', 'Evidence, model against probe',
             '-log10 P under the reported model on x and the probe on y, one point per locus per '
             'cohort. The loci are EXACTLY those cohort_compare draws: diamonds are its panel (b) '
             '(genome-wide in at least one cohort), circles its panel (c) (suggestive in all '
             'three); colour is the cohort. Both figures read the same `panel` column of '
             'lead_crosscohort.tsv, so a locus can be followed between them. The red line is '
             'equality: a point on it is a locus the principal components do not touch. Dotted '
             'lines mark the genome-wide threshold on each axis, so whether a locus crosses under '
             'one model and not the other is readable directly.'),
            ('c', 'Effect, model against probe',
             'The same leads by odds ratio, log axes, with dotted lines at OR = 1. Effects that sit '
             'on the identity line are unchanged by the covariates even where the evidence moves, '
             'which is the usual pattern when the difference is in the standard error rather than '
             'in the estimate.'),
        ],
        numbers=([('cohorts', ', '.join(cohorts)), ('models', ', '.join(models)),
                  ('peak model', args.peak_model)]
                 + [(f'lambda_GC {c} / {m}', v) for (c, m), v in sorted(lam.items())]
                 + [('leads compared', int(tab['variant_id'].nunique()) if len(tab) else 0),
                    ('rows in the forest', len(tab))]),
        tables=[('Every peak lead under both models', tab if len(tab) else None, 40)],
        reading=[
            'Read panel (a) as pairs, not as a series: the two markers above one cohort are the '
            'comparison, and the distance between them is what the PCs are worth there.',
            'A probe already near 1 means the GRM absorbed the structure by itself.',
            'In panels (b) and (c), read distance from the red identity line. A cloud lying on it '
            'means the covariates are doing nothing for these loci.',
            'A point far below the line in (b) but on it in (c) is a locus whose ESTIMATE is '
            'unchanged and whose EVIDENCE weakened — a standard-error effect, not a signal that '
            'depended on the covariates.',
            'lambda is reported and never applied, here as everywhere else in this component.',
        ],
        limits=[
            'It cannot say which model is correct. It reports the difference between them; the '
            'reported model is chosen on design grounds, not on which lambda is smaller.',
            'The probe is not an independent analysis. It shares every sample and the whole GRM '
            'with the reported model, so agreement between them carries no replication weight.',
            'Only the peak model defines loci, so panels (b) and (c) can only show leads that the '
            'reported model already found. A locus the probe alone would have flagged does not '
            'appear.',
            'Suggestive leads are shown to give the comparison a population to be read over. They '
            'remain a description of the scan, not findings.',
        ],
        defs=['lambda_gc', 'gw_sig', 'or'],
        model=S.FORMULAS['lambda'],
        methods_ref='../../../docs/METHODS.md')

    print(f'[plot_model_compare] {len(cohorts)} cohorts x {len(models)} models; '
          f'{tab["variant_id"].nunique() if len(tab) else 0} lead(s) compared -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in plot_model_compare: {e}', file=sys.stderr)
        sys.exit(1)
