#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The one cross-cohort figure. The three sample sets are NESTED —
#           narrow subset of intermediate subset of full — so they are not
#           replication of one another and are never ranked; the figure reports
#           how each quantity moves as the ancestry filter is relaxed.
#             (a) the peak model's calibration over sample composition
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

# The plot block is COMPUTED, not fixed: four blocks (calibration, the
# genome-wide forest, the capped suggestive forest, and the replication summary
# that keeps the cap honest) and two of them are sized by how many rows they
# actually hold. See main().
GW, SUG = 'genome_wide', 'suggestive'


def parse_args():
    p = argparse.ArgumentParser(description='Cross-cohort comparison of one model\'s scan.')
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


def model_row(scan_qc, model):
    """The one scan_qc row for `model`, with a legible error if it is absent.

    `.iloc[0]` on an empty selection raises an index error that says nothing
    about the cause, which is almost always a model name that does not exist in
    this run.
    """
    sel = scan_qc[scan_qc['model'] == model]
    if not len(sel):
        raise SystemExit(f"scan_qc has no row for model '{model}' "
                         f"(has: {sorted(scan_qc['model'].unique())})")
    return sel.iloc[0]


def stack(frames, order, path):
    parts = [frames[c] for c in order if c in frames and len(frames[c])]
    d = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    d.to_csv(path, sep='\t', index=False)
    return d


def calibration_panel(ax_lam, ax_n, scans, peaks, order, colors, SHORT, model):
    """ADDITIVE calibration over ADDITIVE sample composition, on one cohort axis.

    Two sub-panels rather than one dual-axis plot. lambda_GC on the left axis and
    N_eff on a right axis made the reader align two series by eye, and it could
    not show what the stacked bar shows: `narrow` keeps 95 % of the cases but only
    67 % of the controls, which is exactly why it is the sample set most exposed
    to ancestry confounding.
    """
    xs = np.arange(len(order))
    add = {c: model_row(scans[c], model) for c in order}
    lam = [float(add[c]['lambda_gc']) for c in order]

    # ── upper: lambda as a deviation from 1 ────────────────────────────────
    ax_lam.axhline(1.0, color=S.ACCENT, lw=1.1, zorder=2)
    for x, v, c in zip(xs, lam, order):
        ax_lam.plot([x, x], [1.0, v], color=colors[c], lw=3.2,
                    solid_capstyle='round', zorder=3)
        ax_lam.scatter(x, v, s=50, color=colors[c], edgecolors='white',
                       linewidths=0.9, zorder=5)
    ax_lam.set_ylabel(f'{model} {S.LAMBDA_GC}')
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


def forest_rows(cc, order):
    """(rows, n_loci) for a forest, WITHOUT drawing.

    Split from the drawing so the caller can count rows before the figure is
    created and size each forest's slot from its own content. A fixed height
    ratio is right only for the row count it was tuned against, and the two
    forests differ several-fold.

    No cap: (c) is bounded by its own definition — a lead has to be a peak in
    EVERY cohort to appear — so a top-N truncation on top of that would hide
    loci that are by construction few.
    """
    rows, n_loci = [], 0
    if len(cc):
        cc = cc.copy()
        cc['_k'] = pd.to_numeric(cc['chrom'].astype(str).str.replace('chr', '', regex=False),
                                 errors='coerce').fillna(99)
        cc['_p'] = pd.to_numeric(cc['P'], errors='coerce')
        # Strongest first, ranked on the SMALLEST P any of the three cohorts
        # reports — the number the reader is ranking on. `best` was already
        # sorted that way and the order was then thrown away by iterating in
        # genomic order instead.
        best = cc.groupby('variant_id')['_p'].min().sort_values()
        n_loci = len(best)
        sel = cc[cc.variant_id.isin(best.index)]
        for v in list(best.index):
            for c in order:
                m = sel[(sel.variant_id == v) & (sel.cohort == c)]
                if len(m):
                    rows.append((v, c, m.iloc[0]))
    return rows, n_loci


def forest_panel(ax, rows, order, colors, SHORT, empty_text=None, xlim=None,
                 alpha=None, suggestive=None, show_states=True):
    """Draw a forest from `forest_rows` output: one row per (lead x cohort)."""
    if not rows:
        ax.text(0.5, 0.5, empty_text or 'no peak in any cohort', transform=ax.transAxes,
                ha='center', va='center', color=S.REFERENCE)
        ax.set_xticks([]); ax.set_yticks([])
        return 0, []

    yv = np.arange(len(rows))[::-1]
    lo_all, hi_all, states = [], [], []
    for y, (_v, c, r) in zip(yv, rows):
        # State from the variant's OWN P in this cohort — the same quantity the
        # panel membership is decided on. It used to read `called_peak`, which
        # records whether the variant was that cohort's PEAK LEAD; a variant can
        # clear 1e-5 there and still not be the lead, and it was then drawn as
        # "not significant" inside a panel defined by clearing 1e-5 everywhere.
        # The figure contradicted itself.
        pv = pd.to_numeric(r.get('P'), errors='coerce')
        if not np.isfinite(pv):
            state = 'not_in_call_set'
        elif alpha is not None and pv < alpha:
            state = 'genome_wide'
        elif suggestive is not None and pv < suggestive:
            state = 'suggestive'
        else:
            state = 'not_a_peak'
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
                        fontsize=plt.rcParams['legend.fontsize'] - 1, color=S.INK_SOFT)
    ax.axvline(1.0, color=S.REFERENCE, lw=1.0, zorder=2)
    # A gene can carry more than one independent lead — two distinct loci then
    # render as two identically-labelled blocks and the reader cannot tell which
    # is which. Disambiguate by position, but only where the name actually
    # repeats, so the common case stays a bare gene symbol.
    def gene_of(v, r):
        g = str(r.get('Gene'))
        return g if g not in ('.', 'nan', 'None') else str(v)
    per_gene = {}
    for v, _c, r in rows:
        per_gene.setdefault(gene_of(v, r), set()).add(v)
    def group_of(v, r):
        g = gene_of(v, r)
        if len(per_gene.get(g, ())) > 1:
            g = f"{g} {r.get('chrom')}:{int(r.get('pos'))}"
        return g
    # The gene is a property of the BLOCK, not of each row: writing it once beside
    # its three cohorts states the grouping instead of repeating the name three
    # times and leaving the reader to infer which rows belong together.
    # Three stacked lines per locus: the gene symbol names it, the variant id
    # says WHICH variant all three rows are measuring — the locus is drawn at one
    # representative, so that has to be visible — and the rsID is the handle a
    # reader takes to an external resource.
    base = plt.rcParams['ytick.labelsize']
    groups, seen = [], None
    for v, c, r in rows:
        if seen != v:
            rs = S.rsid_text(r)
            groups.append(([
                (group_of(v, r), dict(fontstyle='italic', fontweight='bold', fontsize=base)),
                (str(v),         dict(fontsize=base - 1.0, color=S.INK_SOFT)),
                (rs,             dict(fontsize=base - 1.0, color=S.INK_SOFT)),
            ], []))
            seen = v
        groups[-1][1].append(SHORT[c])
    S.grouped_row_labels(ax, groups)
    ax.set_xscale('log')
    # A reserved right-hand column for the P values, so they cannot land on a
    # confidence interval whatever the widest interval turns out to be.
    # `xlim` is supplied so (b) and (c) share ONE range and the same OR sits at
    # the same x in both; each panel still draws its own ticks and axis title,
    # which sharex would have stripped from the upper one.
    ax.set_xlim(*(xlim if xlim else
                  (min(lo_all) * 0.88 if lo_all else 0.8,
                   (max(hi_all) if hi_all else 2.0) * 2.6)))
    S.or_log_axis(ax)
    ax.set_xlabel(f'{S.OR_SYM} (95% CI), log scale')
    for y, (_v, _c, r) in zip(yv, rows):
        ax.annotate(f"$P$ = {S.p_tex(r['P'])}", xy=(0.995, y),
                    xycoords=('axes fraction', 'data'), ha='right', va='center',
                    fontsize=plt.rcParams['legend.fontsize'] - 0.5, color=S.INK_SOFT)
    # y above the panel-tag band, not level with it: at y=1.01 the bold panel
    # letter rendered inside the first legend key.
    # (c)'s every point is suggestive by construction, so a state legend there
    # would be a constant — it is dropped and the caption says so instead.
    if show_states:
        S.legend_above(ax, S.called_legend_handles(sorted(set(states), key=lambda s: list(
            S.CALLED_STYLE).index(s) if s in S.CALLED_STYLE else 9)), ncol=3, pad_in=0.20)
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
    # Panel membership is assigned ONCE, in cross_cohort.py, and read here:
    #   b = genome-wide in at least one cohort
    #   c = suggestive in EVERY cohort, and not in b
    # Reading the column rather than re-deriving the rule is what keeps this
    # figure and model_compare on the same loci.
    missing = [k for k in ('panel', 'locus_id', 'locus_rep') if len(cc) and k not in cc.columns]
    if missing:
        raise SystemExit(f'lead_crosscohort.tsv has no {missing} column(s) — it was written by an '
                         'older cross_cohort.py; re-run CROSS_COHORT')
    # ONE row-group per locus, drawn at its representative. Without this a locus
    # whose lead differs between cohorts appears twice: chr6:32632724 led in two
    # cohorts and chr6:32632861, 137 bp away, in the third, and conditional
    # analysis says they are one signal.
    if len(cc):
        cc = cc[cc.locus_rep.astype(int) == 1].copy()
    gwcc = cc[cc.panel == 'b'].copy() if len(cc) else pd.DataFrame()
    sgcc = cc[cc.panel == 'c'].copy() if len(cc) else pd.DataFrame()

    add = {c: model_row(scans[c], args.model) for c in order}
    lam = [float(add[c]['lambda_gc']) for c in order]
    ne = [float(add[c].get('n_eff', np.nan)) for c in order]
    n_gw = [int((peaks[c].tier == GW).sum()) if c in peaks and len(peaks[c]) else 0 for c in order]
    n_sug = [int((peaks[c].tier == SUG).sum()) if c in peaks and len(peaks[c]) else 0 for c in order]

    # ── size the two forests from their own row counts ─────────────────────
    # A forest row needs a fixed amount of vertical space to stay readable, so
    # the panel heights are MEASURED from the content rather than tuned: the two
    # forests differ several-fold in length and one ratio cannot serve both. The
    # whole figure grows instead of the rows being crushed.
    rows_b, n_gw_loci = forest_rows(gwcc, order)
    rows_c, n_sug_loci = forest_rows(sgcc, order)
    ROW_IN = 0.175                    # inches per forest row, measured to stay legible
    h_a = 1.95                        # calibration block
    # +0.62 for the legend strip: it now sits ABOVE the panel-tag band, which it
    # used to collide with (the bold letter landed inside the first legend key).
    h_b = max(0.75, len(rows_b) * ROW_IN) + 0.62
    h_c = max(0.75, len(rows_c) * ROW_IN) + 0.62
    plot_h = h_a + h_b + h_c + 1.55                   # + inter-panel gaps and x labels

    fig = plt.figure(figsize=(S.COL_DOUBLE, plot_h))
    # (a) is two sub-panels sharing one cohort axis; they carry a single panel tag.
    # A NESTED gridspec so they sit tight against each other while (a) and (b) keep
    # the full gap — subplots_adjust(hspace) is one value for the whole figure.
    gs = fig.add_gridspec(3, 1, height_ratios=[h_a, h_b, h_c])
    gsa = gs[0].subgridspec(2, 1, height_ratios=[0.62, 0.78], hspace=0.12)
    ax_lam = fig.add_subplot(gsa[0])
    ax_n = fig.add_subplot(gsa[1], sharex=ax_lam)
    ax_for = fig.add_subplot(gs[1])
    ax_sug = fig.add_subplot(gs[2])
    calibration_panel(ax_lam, ax_n, scans, peaks, order, colors, SHORT, args.model)
    # ONE odds-ratio range for both forests, so a reader can carry a position
    # from (b) to (c). Computed over both panels' intervals before either draws.
    iv = [pd.to_numeric(r.get(k), errors='coerce')
          for rows in (rows_b, rows_c) for _v, _c, r in rows for k in ('L95', 'U95')]
    iv = [x for x in iv if np.isfinite(x) and x > 0]
    xlim = (min(iv) * 0.88, max(iv) * 2.6) if iv else None
    n_rows, rows = forest_panel(ax_for, rows_b, order, colors, SHORT,
                                empty_text='no genome-wide peak in any cohort', xlim=xlim,
                                alpha=args.alpha, suggestive=args.alpha_suggestive)
    forest_panel(ax_sug, rows_c, order, colors, SHORT,
                 empty_text='no lead is suggestive in every cohort', xlim=xlim,
                 alpha=args.alpha, suggestive=args.alpha_suggestive, show_states=False)
    for ax, letter in ((ax_lam, 'a'), (ax_for, 'b'), (ax_sug, 'c')):
        S.panel_tag(ax, letter)

    n_v = int(gwcc['variant_id'].nunique()) if len(gwcc) else 0
    S.caption_block(
        fig,
        title=('Additive ' + S.LAMBDA_GC + ' = '
               + ', '.join(f'{v:.3f} ({SHORT[c]})' for c, v in zip(order, lam))
               + f'; {sum(n_gw)} genome-wide and {sum(n_sug)} suggestive {args.model} peaks across the '
                 f'three nested sample sets.'),
        panels=[
            'Additive $\\lambda_{\\mathrm{GC}}$ as a deviation from 1 (upper) over the case/control '
            'composition each cohort is built from (lower).',
            f'{S.OR_SYM} and 95% CI for each of the {n_v} genome-wide lead(s) in every cohort; '
            f'marker shape gives that cohort\'s own call on the variant.',
            (f'The same, for the {n_sug_loci} locus/loci suggestive in ALL THREE cohorts and '
             f'genome-wide in none, so every point is suggestive by construction and carries no '
             f'state marker. Disjoint from (b).'),
        ],
        notes=(f'{args.model} model only; {args.covar_label}. Panels (b) and (c) from '
               '`_comparison/tables/lead_crosscohort.tsv`, partitioned by its `panel` column, so '
               'they cannot share a variant. One row-group per LOCUS, at its representative '
               'variant, so all three rows compare the SAME variant; the group label names the '
               'gene, that variant as CHROM:POS:REF:ALT, and its rsID, so which variant is being '
               'compared is never left to inference. Marker state is that '
               'variant\'s own P in that cohort. The cohorts are NESTED and share every case, so '
               'a locus in (c) has survived the ancestry filter, NOT been replicated. Full '
               'explanation, symbol definitions and the estimator: cohort_compare.md'),
        plot_h=plot_h, top_pad=0.34, hspace=0.42, left='auto', right=0.945,
        margin_axes=[ax_lam, ax_n, ax_for, ax_sug])
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
             f'one model\'s scans alone.'),
            ('b', 'Effect concordance across all three cohorts',
             'A forest plot over (genome-wide lead variant) x (cohort). Every variant that reached '
             'genome-wide significance in *any* cohort is shown in *all three*, because a nested '
             'design always has an estimate in the larger sets. Three states are distinguished, not two: '
             'filled diamond = genome-wide in that cohort, filled circle = suggestive there, open '
             'circle = not significant there. A row with no estimate at all is annotated "not in this '
             'cohort\'s call set". A blank would have been indistinguishable from a missing estimate, '
             'which is why the earliest version of this panel was misleading.'),
            ('c', 'Leads that hold across every cohort',
             'The same forest, for leads called SUGGESTIVE in all three cohorts and genome-wide in '
             'none. Membership of (b) and (c) is assigned once, in cross_cohort.py, and written to '
             'the `panel` column of lead_crosscohort.tsv, so the two panels cannot share a variant '
             'and model_compare can draw exactly the same loci. There is no top-N cap: the '
             'definition itself bounds the panel, because the variant has to clear the suggestive '
             'threshold in every cohort to appear. Each group is labelled with its gene symbol, '
             'the representative variant as CHROM:POS:REF:ALT, and that variant\'s rsID — the '
             'three rows are one variant measured three times, and the label says which one. '
             'What this panel does NOT show is replication. The cohorts are nested — narrow within '
             'intermediate within full — and share every case, so a lead present in all three has '
             'survived the ancestry filter, not been confirmed by independent data. An earlier '
             'version of this figure carried a fourth panel counting cohorts per lead; it is gone '
             'because this panel answers the same question directly and its axis label invited '
             'exactly the reading the design forbids.'),
        ],
        numbers=[('cohorts', ', '.join(SHORT[c] for c in order))]
                + [(f'lambda_GC {args.model}, {SHORT[c]}', float(v)) for c, v in zip(order, lam)]
                + [(f'N_eff, {SHORT[c]}', float(v)) for c, v in zip(order, ne)]
                + [(f'N cases / controls, {SHORT[c]}',
                    f'{int(add[c]["n_case"]):,} / {int(add[c]["n_ctrl"]):,}') for c in order]
                + [(f'genome-wide {args.model} peaks, {SHORT[c]}', v) for c, v in zip(order, n_gw)]
                + [(f'suggestive {args.model} peaks, {SHORT[c]}', v) for c, v in zip(order, n_sug)]
                + [('distinct genome-wide lead variants (b)', n_v),
                   ('leads suggestive in every cohort (c)', n_sug_loci),
                   ('forest rows in (b)', n_rows),
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
            f'Panel (b) covers the {args.model} model only. Any other model\'s peaks are in each '
            'cohort\'s `03.peaks/model_peaks_annotation.tsv` and on their own scan figures.',
            'lambda_GC in (a) cannot separate confounding from polygenicity; see METHODS §7.',
        ],
        defs=['lambda_gc', 'neff', 'gw_sig', 'or', 'called_peak', 'model'],
        model=S.FORMULAS['lambda'] + '\n' + S.FORMULAS['neff'],
        methods_ref='../../docs/METHODS.md')

    print(f'[plot_cohort_compare] {len(order)} cohorts; {args.model} genome-wide {n_gw}, '
          f'suggestive {n_sug}; {n_v} genome-wide lead(s) x {len(order)} cohorts = {n_rows} '
          f'forest rows -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in plot_cohort_compare: {e}', file=sys.stderr)
        sys.exit(1)
