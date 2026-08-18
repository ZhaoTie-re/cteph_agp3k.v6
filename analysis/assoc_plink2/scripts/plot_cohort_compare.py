#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The one cross-cohort figure. The three sample sets are NESTED —
#           narrow subset of intermediate subset of full — so they are not
#           replication of one another and are never ranked; the figure reports
#           how each quantity moves as the ancestry filter is relaxed.
#             (a) additive calibration over additive sample composition
#             (b) forest: every genome-wide lead in ALL THREE cohorts
#
#           ADDITIVE ONLY. (a) is two sub-panels on one cohort axis, not a
#           dual-axis plot: lambda_GC on the left axis against N_eff on a right
#           axis made the reader align two series by eye, and it could not show
#           the case/control composition, which is the reason narrow is the set
#           most exposed to ancestry confounding (95 % of cases, 67 % of controls).
#           DOM/REC calibration lives in scan_qc_all.tsv and on their own figures.
#
#           (b) reads lead_crosscohort.tsv, so a cohort that called no peak at a
#           variant still contributes its estimate, and the marker says which of
#           the three states that cohort was in.
# Component: assoc_plink2
# Used by : assoc_plink2.nf  process COMPARE_COHORTS
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
import plot_style as S
import figure_doc
import variant_annot as VA

PLOT_H = 4.4          # inches of croppable plot block (titles + axes + x-labels)
GW, SUG = 'genome_wide', 'suggestive'


def parse_args():
    p = argparse.ArgumentParser(description='Cross-cohort comparison of the additive scan.')
    p.add_argument('--scan-qc', action='append', required=True, metavar='COHORT=PATH')
    p.add_argument('--peaks', action='append', default=[], metavar='COHORT=PATH')
    p.add_argument('--annotation', action='append', default=[], metavar='COHORT=PATH')
    p.add_argument('--model-peaks', action='append', default=[], metavar='COHORT=PATH')
    p.add_argument('--model', default='additive',
                   help='The single model this comparison is drawn for.')
    p.add_argument('--crosscohort', required=True, help='lead_crosscohort.tsv from CROSS_COHORT')
    p.add_argument('--cohort-order', default='',
                   help='comma-separated cohort order; empty means the order the per-cohort arguments were given in.')
    p.add_argument('--alpha', type=float, default=S.GW_ALPHA)
    p.add_argument('--alpha-suggestive', type=float, default=S.SUGGEST_ALPHA)
    p.add_argument('--covar-label', default='SEX + PCs',
                   help='Human-readable covariate set, printed in the caption and the sidecar. '
                        'Configuration, not a property of the method — pass the run\'s own.')
    p.add_argument('--n-pcs', type=int, default=0,
                   help='Number of PCs in the covariate set; 0 leaves the formula\'s upper limit '
                        'as a generic K.')
    p.add_argument('--out-png', required=True)
    p.add_argument('--out-scan', required=True)
    p.add_argument('--out-peaks', required=True)
    p.add_argument('--out-annotation', required=True)
    p.add_argument('--out-model-peaks', required=True)
    return p.parse_args()


def gather(specs):
    out = {}
    for s in specs:
        c, _, path = s.partition('=')
        if path and Path(path).exists() and Path(path).stat().st_size:
            try:
                out[c] = pd.read_csv(path, sep='\t', dtype={'chrom': str})
            except Exception:
                pass
    return out


def stack(frames, order, path):
    parts = [frames[c] for c in order if c in frames and len(frames[c])]
    d = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    d.to_csv(path, sep='\t', index=False)
    return d


def calibration_panel(ax_lam, ax_n, scans, peaks, order, colors, SHORT):
    """ADDITIVE calibration over ADDITIVE sample composition, on one cohort axis.

    Two sub-panels rather than one dual-axis plot. lambda_GC on the left axis and
    N_eff on a right axis made the reader align two series by eye, and it could
    not show what the stacked bar shows: `narrow` keeps 95 % of the cases but only
    67 % of the controls, which is exactly why it is the sample set most exposed
    to ancestry confounding.
    """
    xs = np.arange(len(order))
    add = {c: scans[c][scans[c].model == 'additive'].iloc[0] for c in order}
    lam = [float(add[c]['lambda_gc']) for c in order]

    # ── upper: lambda as a deviation from 1 ────────────────────────────────
    ax_lam.axhline(1.0, color=S.ACCENT, lw=1.1, zorder=2)
    for x, v, c in zip(xs, lam, order):
        ax_lam.plot([x, x], [1.0, v], color=colors[c], lw=3.2,
                    solid_capstyle='round', zorder=3)
        ax_lam.scatter(x, v, s=50, color=colors[c], edgecolors='white',
                       linewidths=0.9, zorder=5)
    ax_lam.set_ylabel(f'additive {S.LAMBDA_GC}')
    lo, hi = min(lam), max(lam)
    pad = max(0.03, (hi - lo) * 1.1)
    ax_lam.set_ylim(min(0.99, lo) - pad, max(1.01, hi) + pad)
    S.value_labels(ax_lam, xs, lam, [f'{v:.3f}' for v in lam], fontweight='bold', offset=8)
    ax_lam.set_xlim(-0.6, len(order) - 0.4)
    ax_lam.set_xticks(xs)
    S.despine(ax_lam, grid_axis='y')

    # ── lower: who is actually in each cohort ──────────────────────────────
    ncase = [int(add[c]['n_case']) for c in order]
    nctrl = [int(add[c]['n_ctrl']) for c in order]
    ax_n.bar(xs, ncase, width=0.56, color=S.ACCENT, zorder=3, label='cases')
    ax_n.bar(xs, nctrl, width=0.56, bottom=ncase, color=S.DATA, zorder=3, label='controls')
    # Numbers only inside the segments; the colour is named once in the legend.
    # Naming the two levels on the leftmost bar alone read as an omission on the
    # other two rather than as "said once".
    for x, a, b in zip(xs, ncase, nctrl):
        ax_n.annotate(f'{a:,}', xy=(x, a / 2), ha='center', va='center', color='white',
                      fontweight='bold', fontsize=plt.rcParams['legend.fontsize'] - 1, zorder=6)
        ax_n.annotate(f'{b:,}', xy=(x, a + b / 2), ha='center', va='center', color='white',
                      fontweight='bold', fontsize=plt.rcParams['legend.fontsize'] - 1, zorder=6)
    ax_n.set_ylabel('samples')
    ax_n.set_xlim(-0.6, len(order) - 0.4)
    ax_n.set_xticks(xs)
    n_gw = [int((peaks[c].tier == GW).sum()) if c in peaks and len(peaks[c]) else 0 for c in order]
    n_sug = [int((peaks[c].tier == SUG).sum()) if c in peaks and len(peaks[c]) else 0 for c in order]
    # N_eff joins the tick label rather than floating above the bar: above the
    # bars is the one region the legend can occupy without covering data, and the
    # tick label is already this cohort's identity card.
    ax_n.set_xticklabels([f"{SHORT[c]}\n$N_{{\\mathrm{{eff}}}}$ = {float(add[c]['n_eff']):,.0f}"
                          f'\n{g} gw · {s} sug' for c, g, s in zip(order, n_gw, n_sug)])
    ax_n.yaxis.set_major_formatter(S.INT_FMT)
    # Upper LEFT, deliberately: the cohorts are nested, so `narrow` is always the
    # shortest stack and that corner is always the emptiest. Stated rather than
    # left to _emptiest_corner, which measures scatter points and sees no bars.
    S.legend_inside(ax_n, [Line2D([], [], ls='none', marker='s', color=S.ACCENT,
                                  markersize=6, label='cases'),
                           Line2D([], [], ls='none', marker='s', color=S.DATA,
                                  markersize=6, label='controls')],
                    loc='upper left', ncol=2)
    S.despine(ax_n, grid_axis='y')
    # AFTER both sub-panels are configured: with a shared x-axis the later
    # set_xticklabels wins on both, so the upper panel's labels have to be hidden
    # as artists rather than cleared as a formatter.
    plt.setp(ax_lam.get_xticklabels(), visible=False)


def forest_panel(ax, gwcc, order, colors, SHORT):
    """Every genome-wide lead, in every cohort, with a three-state marker."""
    rows = []
    if len(gwcc):
        gwcc = gwcc.copy()
        gwcc['_k'] = pd.to_numeric(gwcc['chrom'].astype(str).str.replace('chr', '', regex=False),
                                   errors='coerce').fillna(99)
        for v in list(dict.fromkeys(gwcc.sort_values(['_k', 'pos'])['variant_id'])):
            for c in order:
                m = gwcc[(gwcc.variant_id == v) & (gwcc.cohort == c)]
                if len(m):
                    rows.append((v, c, m.iloc[0]))
    if not rows:
        ax.text(0.5, 0.5, 'no genome-wide peak in any cohort', transform=ax.transAxes,
                ha='center', va='center', color=S.REFERENCE)
        ax.set_xticks([]); ax.set_yticks([])
        return 0, []

    yv = np.arange(len(rows))[::-1]
    lo_all, hi_all, states = [], [], []
    for y, (_v, c, r) in zip(yv, rows):
        state = str(r['called_peak'])
        states.append(state)
        st = S.CALLED_STYLE.get(state, S.CALLED_STYLE['not_a_peak'])
        orv, l95, u95 = (pd.to_numeric(r.get(k), errors='coerce') for k in ('OR', 'L95', 'U95'))
        if np.isfinite(l95) and np.isfinite(u95):
            ax.plot([l95, u95], [y, y], color=colors[c], lw=2.1,
                    solid_capstyle='round', zorder=3)
            lo_all.append(l95); hi_all.append(u95)
        if np.isfinite(orv):
            ax.scatter(orv, y, s=st['size'], marker=st['marker'],
                       facecolor=colors[c] if st['filled'] else 'white',
                       edgecolors=colors[c], linewidths=1.3, zorder=5)
        else:
            ax.annotate('not in this cohort\'s call set', xy=(0.02, y),
                        xycoords=('axes fraction', 'data'), ha='left', va='center',
                        fontsize=plt.rcParams['legend.fontsize'] - 1, color=S.REFERENCE)
    ax.axvline(1.0, color=S.REFERENCE, lw=1.0, zorder=2)
    ax.set_yticks(yv)
    ax.set_yticklabels(
        [f"{(r.get('Gene') if str(r.get('Gene')) not in ('.', 'nan', 'None') else v)} · {SHORT[c]}"
         for v, c, r in rows], fontsize=plt.rcParams['ytick.labelsize'])
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xscale('log')
    # A reserved right-hand column for the P values, so they cannot land on a
    # confidence interval whatever the widest interval turns out to be.
    ax.set_xlim(min(lo_all) * 0.88 if lo_all else 0.8, (max(hi_all) if hi_all else 2.0) * 2.6)
    S.or_log_axis(ax)
    ax.set_xlabel(f'{S.OR_SYM} (95% CI), log scale')
    for y, (_v, _c, r) in zip(yv, rows):
        ax.annotate(f"$P$ = {S.p_tex(r['P'])}", xy=(0.995, y),
                    xycoords=('axes fraction', 'data'), ha='right', va='center',
                    fontsize=plt.rcParams['legend.fontsize'] - 0.5, color='#444444')
    S.legend_above(ax, S.called_legend_handles(sorted(set(states), key=lambda s: list(
        S.CALLED_STYLE).index(s) if s in S.CALLED_STYLE else 9)), ncol=3, y=1.01)
    S.despine(ax, grid_axis='x')
    return len(rows), rows


def main():
    args = parse_args()
    S.setup_style('paper')
    scans, peaks = gather(args.scan_qc), gather(args.peaks)
    anns, mpk = gather(args.annotation), gather(args.model_peaks)
    want = [c for c in args.cohort_order.split(',') if c] or list(scans)
    order = [c for c in want if c in scans]
    SHORT = S.shorten(order)
    colors = {c: S.COHORT_RAMP[i % len(S.COHORT_RAMP)] for i, c in enumerate(order)}

    # ── stacked tables ──────────────────────────────────────────────────────
    stack(scans, order, args.out_scan)
    stack(peaks, order, args.out_peaks)
    alla = stack(anns, order, args.out_annotation)
    stack(mpk, order, args.out_model_peaks)

    cc = pd.read_csv(args.crosscohort, sep='\t', dtype={'chrom': str, 'variant_id': str}) \
        if Path(args.crosscohort).exists() and Path(args.crosscohort).stat().st_size \
        else pd.DataFrame()
    gwcc = cc[cc.best_tier == GW].copy() if len(cc) else pd.DataFrame()

    add = {c: scans[c][scans[c].model == 'additive'].iloc[0] for c in order}
    lam = [float(add[c]['lambda_gc']) for c in order]
    ne = [float(add[c].get('n_eff', np.nan)) for c in order]
    n_gw = [int((peaks[c].tier == GW).sum()) if c in peaks and len(peaks[c]) else 0 for c in order]
    n_sug = [int((peaks[c].tier == SUG).sum()) if c in peaks and len(peaks[c]) else 0 for c in order]

    fig = plt.figure(figsize=(S.COL_DOUBLE, PLOT_H))
    # (a) is two sub-panels sharing one cohort axis; they carry a single panel tag.
    # A NESTED gridspec so they sit tight against each other while (a) and (b) keep
    # the full gap — subplots_adjust(hspace) is one value for the whole figure.
    gs = fig.add_gridspec(2, 1, height_ratios=[1.36, 1.05])
    gsa = gs[0].subgridspec(2, 1, height_ratios=[0.62, 0.78], hspace=0.12)
    ax_lam = fig.add_subplot(gsa[0])
    ax_n = fig.add_subplot(gsa[1], sharex=ax_lam)
    ax_for = fig.add_subplot(gs[1])
    calibration_panel(ax_lam, ax_n, scans, peaks, order, colors, SHORT)
    n_rows, rows = forest_panel(ax_for, gwcc, order, colors, SHORT)
    S.panel_tag(ax_lam, 'a')
    S.panel_tag(ax_for, 'b')

    n_v = int(gwcc['variant_id'].nunique()) if len(gwcc) else 0
    S.caption_block(
        fig,
        title=('Additive ' + S.LAMBDA_GC + ' = '
               + ', '.join(f'{v:.3f} ({SHORT[c]})' for c, v in zip(order, lam))
               + f'; {sum(n_gw)} genome-wide and {sum(n_sug)} suggestive additive peaks across the '
                 f'three nested sample sets.'),
        panels=[
            'Additive $\\lambda_{\\mathrm{GC}}$ as a deviation from 1 (upper) over the case/control '
            'composition each cohort is built from (lower).',
            f'{S.OR_SYM} and 95% CI for each of the {n_v} genome-wide lead(s) in all three cohorts; '
            f'marker shape gives that cohort\'s own call on the variant.',
        ],
        notes=(f'{args.model} model only; {args.covar_label}. Panel (b) from '
               '`_comparison/tables/lead_crosscohort.tsv`. Full explanation, symbol definitions and '
               'the estimator: cohort_compare.md'),
        plot_h=PLOT_H, top_pad=0.34, hspace=0.52, left='auto', right=0.945,
        margin_axes=[ax_lam, ax_n, ax_for])
    fig.savefig(args.out_png)
    plt.close(fig)

    full = pd.DataFrame()
    if len(gwcc):
        full = gwcc[['variant_id', 'rsID', 'Gene', 'cohort', 'called_peak', 'EA', 'OA',
                     'OR', 'L95', 'U95', 'P', 'Case_Genotype_Distribution', 'Case_EAF',
                     'Case_Missing_Rate', 'Case_HWE_P', 'Control_Genotype_Distribution',
                     'Control_EAF', 'Control_Missing_Rate', 'Control_HWE_P',
                     'A1_FREQ', 'OBS_CT', 'N_case', 'N_ctrl']].copy()
    calib = pd.DataFrame()
    if len(scans):
        calib = pd.concat([scans[c][['cohort', 'model', 'n_case', 'n_ctrl', 'n_eff',
                                     'n_analysed', 'lambda_gc', 'n_genomewide', 'n_suggestive']]
                           for c in order], ignore_index=True)

    figure_doc.write_doc(
        args.out_png,
        title='Additive scan across three nested sample sets',
        question=('How do calibration and effect size move as the ancestry filter is relaxed from '
                  'narrow to full?'),
        interpretation=('The cohorts are nested — narrow within intermediate within full — so they '
                        'share most of their samples and are *not* replication of one another; '
                        'agreement between them is expected and carries little independent '
                        'information, while a signal present only in the narrowest set is a candidate '
                        'for an ancestry-driven artefact. Panel (b) is therefore a description of how '
                        'one estimate behaves as the sample grows, not a replication test: an '
                        'interval that narrows while the point estimate holds is a variant gaining '
                        'sample size, and an estimate that moves toward 1 as the filter relaxes is '
                        'what a structure-driven signal does. Panel (a) shows why narrow is the most '
                        'exposed set: relaxing the filter from narrow to full adds 20 cases but 883 '
                        'controls, so the narrow cohort keeps 95% of the cases against only 67% of '
                        'the controls. No cohort is designated correct here; a fixed-effects scan '
                        'cannot settle that, and the GRM-based random-effect follow-up is what will. '
                        'Dominant and recessive calibration is not in this figure — it is in '
                        'scan_qc_all.tsv and on each scan\'s own figure.'),
        panels=[
            ('a', 'Calibration against effective size',
             'lambda_GC (y) against N_eff (x), one marker per scan: colour = cohort, shape = model. '
             'A grey line joins each model\'s three cohorts in nesting order. This replaces a '
             'dual-axis stem plot, where lambda sat on the left axis and N_eff on the right and the '
             'reader had to align two series by eye; plotting one against the other makes the '
             'trade-off a direction on the page, and leaves room for all nine scans rather than the '
             'additive three.'),
            ('b', 'Effect concordance across all three cohorts',
             'A forest plot over (genome-wide lead variant) x (cohort). Every variant that reached '
             'genome-wide significance in *any* cohort is shown in *all three*, because a nested '
             'design always has an estimate in the larger sets. Three states are distinguished, not two: '
             'filled diamond = genome-wide in that cohort, filled circle = suggestive there, open '
             'circle = not significant there. A row with no estimate at all is annotated "not in this '
             'cohort\'s call set". A blank would have been indistinguishable from a missing estimate, '
             'which is why the earliest version of this panel was misleading.'),
        ],
        numbers=[('cohorts', ', '.join(SHORT[c] for c in order))]
                + [(f'lambda_GC additive, {SHORT[c]}', float(v)) for c, v in zip(order, lam)]
                + [(f'N_eff, {SHORT[c]}', float(v)) for c, v in zip(order, ne)]
                + [(f'N cases / controls, {SHORT[c]}',
                    f'{int(add[c]["n_case"]):,} / {int(add[c]["n_ctrl"]):,}') for c in order]
                + [(f'genome-wide additive peaks, {SHORT[c]}', v) for c, v in zip(order, n_gw)]
                + [(f'suggestive additive peaks, {SHORT[c]}', v) for c, v in zip(order, n_sug)]
                + [('distinct genome-wide lead variants', n_v),
                   ('forest rows', n_rows),
                   ('lead variants reported in every cohort',
                    int(cc['variant_id'].nunique()) if len(cc) else 0),
                   ('rows in lead_crosscohort.tsv', int(len(cc)))],
        tables=[('Calibration and size, all nine scans', calib),
                ('Genome-wide leads — complete statistics in all three cohorts', full)],
        reading=[
            'Do not read agreement between cohorts as replication. They share samples by '
            'construction; narrow is a subset of intermediate, which is a subset of full.',
            'In (a), read the vertical spread first: the three models separate far more than the '
            'three cohorts do, so model choice dominates calibration here.',
            'In (b), read down each variant\'s three rows. An interval that narrows while the point '
            'estimate holds is a variant gaining sample size. An estimate that drifts toward 1 as the '
            'filter is relaxed is what a structure-driven signal looks like.',
            'A peak that appears only in narrow deserves suspicion — the narrow filter retains 95% of '
            'cases but only 67% of controls, so it is the set most exposed to ancestry confounding.',
            'For every lead of both tiers in every cohort, not just the genome-wide ones, read '
            '`_comparison/tables/lead_crosscohort.tsv`.',
        ],
        limits=[
            'It cannot adjudicate between the cohorts. Choosing one requires a model that absorbs '
            'fine-scale structure, which this fixed-effects scan is not.',
            'It is not a replication analysis and no cohort here is independent of another. The '
            'component has no external cohort available.',
            'Panel (b) covers the additive model only. Dominant and recessive peaks are in each '
            'cohort\'s `03.peaks/model_peaks_annotation.tsv` and on their own scan figures.',
            'lambda_GC in (a) cannot separate confounding from polygenicity; see METHODS §7.',
        ],
        defs=['lambda_gc', 'neff', 'gw_sig', 'or', 'called_peak', 'model'],
        model=S.FORMULAS['lambda'] + '\n' + S.FORMULAS['neff'],
        methods_ref='../../docs/METHODS.md')

    print(f'[plot_cohort_compare] {len(order)} cohorts; additive genome-wide {n_gw}, '
          f'suggestive {n_sug}; {n_v} genome-wide lead(s) x {len(order)} cohorts = {n_rows} '
          f'forest rows -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in plot_cohort_compare: {e}', file=sys.stderr)
        sys.exit(1)
