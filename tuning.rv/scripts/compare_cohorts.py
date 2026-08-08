#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Cross-cohort comparison for the rare-variant depth-confounding
#           tuning. Overlays the per-minAC sweep of every cohort so SAMPLE-SELECTION
#           sensitivity is visible at a glance, in a 2-row grid (top = no-PC primary,
#           bottom = +ancestry-PC sensitivity) when PC columns are present:
#             (a) variants retained          vs minAC  (power cost; log y; shown once)
#             (b) OLS beta_g (no-PC)          vs minAC  (adjusted apparent effect; DECISION)
#             (c) OLS beta_d (no-PC)          vs minAC  (residual depth effect; DIAGNOSTIC)
#             (d) OLS beta_g (+PC)            vs minAC  (ancestry-adjusted apparent effect)
#             (e) OLS beta_d (+PC)            vs minAC  (ancestry-adjusted depth effect)
#           (If no PC columns exist it falls back to the 1x3 no-PC layout a/b/c.)
#           Primary model is no-PC (beta_g is the decision axis; beta_d a diagnostic).
#           Each cohort's recommended minAC (V-trough) is marked below its beta_g
#           panel — for the no-PC model under (b) and the +PC model under (d). A
#           per-cohort summary table is written alongside. This compares nested
#           SAMPLE sets — it does NOT break the platform confounding, identical in
#           every cohort.
# Project : cteph_agp3k.v6 / tuning.rv  (rare-variant depth-confounding QC)
# Used by : tuning.rv.nf  process COMPARE_COHORTS
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import plot_style as ps

# Distinct, colorblind-safe color per cohort (assigned in first-seen order).
COHORT_PALETTE = [ps.RHO_COLOR, ps.BETA_COLOR, ps.GROUP_COLOR, '#CC79A7', '#E69F00']


def _cohort_col(df):
    for c in ('Cohort', 'cohort', 'COHORT'):
        if c in df.columns:
            return c
    raise ValueError(f'no cohort column in {list(df.columns)}')


def _beta_cols(df):
    """PRIMARY = no-PC beta_group columns (isolates the technical artifact); fall back to +PC."""
    if 'Beta_group_noPC' in df.columns and df['Beta_group_noPC'].notna().any():
        return 'Beta_group_noPC', 'Beta_group_noPC_L', 'Beta_group_noPC_U', 'P_group_noPC', 'no-PC'
    return 'Beta_group_pc', 'Beta_group_pc_L', 'Beta_group_pc_U', 'P_group_pc', '+PC'


def _recommend(df, b, l, u, p, alpha=0.05, floor=1):
    """V-trough: smallest minAC (>= floor) where p>=alpha and 0 in CI; else argmin|beta|."""
    d = df[df['MinAC_Threshold'] >= floor]
    if d.empty or b not in d.columns:
        return None
    cal = (d[p].fillna(0.0) >= alpha) & (d[l] <= 0) & (d[u] >= 0)
    band = d.loc[cal, 'MinAC_Threshold']
    if len(band):
        return int(band.min())                       # left edge of the non-sig trough (max power)
    dd = d.dropna(subset=[b])
    return int(dd.loc[dd[b].abs().idxmin(), 'MinAC_Threshold']) if len(dd) else None  # least-confounded


def main():
    ap = argparse.ArgumentParser(description='Compare rare-variant depth-confounding across cohorts.')
    ap.add_argument('--stats-combined', required=True,
                    help='Long-format all-cohort per-minAC stats (with a Cohort column).')
    ap.add_argument('--summary-combined', default=None,
                    help='Long-format all-cohort summary (with a cohort column); optional.')
    ap.add_argument('--out-png', required=True)
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--alpha', type=float, default=0.05)
    args = ap.parse_args()

    try:
        stats = pd.read_csv(args.stats_combined, sep='\t')
        ccol = _cohort_col(stats)
        labels = list(dict.fromkeys(stats[ccol].tolist()))    # preserve first-seen order
        data = {lab: stats[stats[ccol] == lab].sort_values('MinAC_Threshold') for lab in labels}
        colors = {lab: COHORT_PALETTE[i % len(COHORT_PALETTE)] for i, lab in enumerate(labels)}
        b, l, u, p, beta_tag = _beta_cols(stats)          # primary (no-PC) cols — drive the table
        has_pc = ('Beta_group_pc' in stats.columns) and stats['Beta_group_pc'].notna().any()
        # Recommendation per model: no-PC drives the recommendation/table; +PC is a sensitivity.
        gnoPC = ('Beta_group_noPC', 'Beta_group_noPC_L', 'Beta_group_noPC_U', 'P_group_noPC')
        gpc = ('Beta_group_pc', 'Beta_group_pc_L', 'Beta_group_pc_U', 'P_group_pc')
        rec_noPC = {lab: _recommend(data[lab], *gnoPC, alpha=args.alpha) for lab in labels}
        rec_pc = {lab: _recommend(data[lab], *gpc, alpha=args.alpha) for lab in labels} if has_pc else {}
        rec = rec_noPC                                     # the recommendation (no-PC primary)

        ps.setup_style('slide')
        if has_pc:
            # Full 2x3 with a meaning on both axes: rows are the two models, columns are
            # beta_group / beta_depth / (model-independent or summary). That leaves no
            # vacant slot, and puts the decision summary bottom-right where a conclusion
            # belongs. The legend goes back above the panels.
            fig, axes = plt.subplots(2, 3, figsize=(17.0, 15.0))
            (axG0, axD0, axPow), (axG1, axD1, axFor) = axes
        else:                                              # fallback: 1x3 no-PC layout
            fig, (axG0, axD0, axPow) = plt.subplots(1, 3, figsize=(17.0, 6.9))
            axG1 = axD1 = axFor = None
        pc_lab = (str(stats['PC_source'].dropna().iloc[0])
                  if ('PC_source' in stats and stats['PC_source'].notna().any()) else 'none')
        n_pc = int(stats['N_PC'].dropna().iloc[0]) if ('N_PC' in stats and stats['N_PC'].notna().any()) else None
        pc_txt = f'{pc_lab}' + (f', {n_pc} PCs' if n_pc else '')
        kx = 'minAC threshold (minAC)'
        yg = r'OLS $\beta_{\mathrm{group}}$  (Case $-$ Control)'
        yd = r'OLS $\beta_{\mathrm{depth}}$  (per $+1$ SD depth)'

        def overlay(ax, col, ylabel, letter, title, logy=False, zero=False, pcol=None):
            """Overlay every cohort's trajectory; on the coefficient panels also encode
            whether each estimate is significant.

            Three cohorts x a CI ribbon would be unreadable, so significance is carried
            by the marker (filled = p < alpha, open = not) rather than by a band — the
            estimate, its uncertainty verdict and the cohort all stay legible at once.
            """
            for lab in labels:
                df = data[lab]
                if col not in df.columns:
                    continue
                d = df[df['MinAC_Threshold'] != 0]     # sweep panels drop the singleton-heavy minAC=0
                x, y = d['MinAC_Threshold'], d[col]
                ax.plot(x, y, color=colors[lab], linewidth=1.9, zorder=3)
                if pcol and pcol in d.columns:
                    ps.sig_markers(ax, x, y, d[pcol], colors[lab], alpha=args.alpha, size=34, zorder=4)
                else:
                    ax.plot(x, y, color=colors[lab], linestyle='', marker='o', markersize=4.5,
                            markeredgecolor='white', markeredgewidth=0.5, zorder=4)
            if zero:
                ax.axhline(0, color='#888888', linewidth=1.0, zorder=1)
            if logy:
                ax.set_yscale('log')
            ax.set_ylabel(ylabel)
            ax.set_xlabel(kx)
            ps.panel_tag(ax, letter, title)           # bold left letter + short centred title
            ps.despine(ax)

        overlay(axG0, 'Beta_group_noPC', yg, 'a', 'Apparent effect (no-PC)', zero=True,
                pcol='P_group_noPC')
        overlay(axD0, 'Beta_depth_noPC', yd, 'b', 'Depth effect (no-PC)', zero=True,
                pcol='P_depth_noPC')
        overlay(axPow, 'Variant_Count', 'variants retained', 'c', 'Power cost', logy=True)
        if has_pc:
            overlay(axG1, 'Beta_group_pc', yg, 'd', 'Apparent effect (+PC)', zero=True,
                    pcol='P_group_pc')
            overlay(axD1, 'Beta_depth_pc', yd, 'e', 'Depth effect (+PC)', zero=True,
                    pcol='P_depth_pc')
            # The primary model (top) and its sensitivity (bottom) are only comparable
            # by eye if the paired panels share a scale — otherwise a flat +PC curve can
            # look identical to a steep no-PC one.
            for a_top, a_bot in ((axG0, axG1), (axD0, axD1)):
                lo = min(a_top.get_ylim()[0], a_bot.get_ylim()[0])
                hi = max(a_top.get_ylim()[1], a_bot.get_ylim()[1])
                a_top.set_ylim(lo, hi)
                a_bot.set_ylim(lo, hi)

        # Each cohort's recommended minAC as a dodged marker near the BOTTOM of its
        # beta_g panel (dodged so triangles don't stack when cohorts share a minAC).
        from matplotlib.transforms import blended_transform_factory
        nlab = len(labels)

        def mark_recs(ax, recd):
            # Placed just BELOW the axis (clip_on=False) pointing up at the curve, not
            # at y=0.05 inside the panel where they landed on top of the trajectories
            # in exactly the low-minAC region the reader needs to see.
            tr = blended_transform_factory(ax.transData, ax.transAxes)
            for j, lab in enumerate(labels):
                xr = recd.get(lab)
                if xr is None:
                    continue
                dodge = (j - (nlab - 1) / 2.0) * 0.18
                ax.plot([xr + dodge], [-0.045], marker='^', markersize=10, color=colors[lab],
                        markeredgecolor='black', markeredgewidth=0.6, transform=tr,
                        clip_on=False, zorder=8)

        mark_recs(axG0, rec_noPC)
        if has_pc:
            mark_recs(axG1, rec_pc)

        # ── (f) decision summary: what survives all three filtering steps ────────
        # The sweep panels show trajectories; this reads off the single number the
        # whole pipeline exists to produce — beta_group at each cohort's OWN chosen
        # kappa, under both models, against zero.
        if has_pc:
            rows = []
            for lab in labels:
                for tag, (bc, lc, uc, pcc), recd in (
                        ('no-PC', gnoPC, rec_noPC), ('+PC', gpc, rec_pc)):
                    k = recd.get(lab)
                    if k is None:
                        continue
                    r = data[lab][data[lab]['MinAC_Threshold'] == k]
                    if not len(r):
                        continue
                    r = r.iloc[0]
                    # Two lines: a single-line label was wide enough to land inside the
                    # neighbouring panel's plotting area.
                    rows.append((f'{lab.split("_")[0]}\n{tag}, κ={k}', r[bc], r[lc], r[uc],
                                 r[pcc], colors[lab]))
            ypos = np.arange(len(rows))[::-1]                    # first row at the top
            # NB: do not name these b / l / u / p — those hold the primary COLUMN NAMES
            # from _beta_cols() and are read again when the summary table is built below.
            for y, (lbl, bval, lo_v, hi_v, pv, col) in zip(ypos, rows):
                sig = pd.notna(pv) and pv < args.alpha
                axFor.plot([lo_v, hi_v], [y, y], color=col, linewidth=2.0, zorder=3,
                           solid_capstyle='butt')
                axFor.plot([bval], [y], marker='o', markersize=8, zorder=4,
                           color=col if sig else 'white', markeredgecolor=col,
                           markeredgewidth=1.6 if not sig else 0.8)
            axFor.axvline(0, color='#888888', linewidth=1.1, zorder=1)
            axFor.set_yticks(ypos)
            axFor.set_yticklabels([r[0] for r in rows], fontsize=plt.rcParams['ytick.labelsize'] - 1)
            axFor.set_ylim(-0.7, len(rows) - 0.3)
            axFor.set_xlabel(yg)
            ps.panel_tag(axFor, 'f', 'Decision summary')
            ps.despine(axFor, grid_axis='x')

        # Dense minAC axis (1..MaxAC): label a readable subset, minor-tick the rest.
        xs = sorted({int(v) for df in data.values() for v in df['MinAC_Threshold'] if v != 0})
        for ax in [axPow, axG0, axD0] + ([axG1, axD1] if has_pc else []):
            ps.sparse_int_ticks(ax, xs)

        # Figure key. Cohort names are NOT abbreviated to their first token any more:
        # 'narrow' / 'intermediate' alone do not say what the set is, and the sample
        # size is the thing a reader wants next to the name.
        def _n_of(lab):
            s = data[lab]['Sample_Count'].dropna()
            return f'  (n = {int(s.iloc[0]):,})' if len(s) else ''

        handles = [Line2D([], [], color=colors[lab], lw=2.4, marker='o', markersize=6,
                          markeredgecolor='white', label=f'{lab}{_n_of(lab)}') for lab in labels]
        handles.append(Line2D([], [], color='#666666', marker='^', linestyle='', markersize=10,
                              markeredgecolor='black', label='recommended minAC'))
        handles += ps.sig_legend_handles(args.alpha)

        # Key above the panels: the grid is now full, so there is no spare slot for it.
        if has_pc:
            # One row: wrapped onto two, the key reached down into the top panels.
            fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.955),
                       ncol=len(handles), frameon=False, handletextpad=0.5, columnspacing=1.6)
            sup_y, cap_top = 0.988, 0.912
        else:
            fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.915),
                       ncol=len(handles), frameon=False, handletextpad=0.5, columnspacing=1.8)
            sup_y, cap_top = 0.985, 0.80

        fig.suptitle('Cross-cohort comparison — sample-selection robustness of the minAC choice',
                     fontsize=16, y=sup_y)
        rec_txt = ', '.join(f'{lab} {rec_noPC[lab]}' for lab in labels if rec_noPC.get(lab) is not None)
        # Grid: rows = model (no-PC primary / +PC sensitivity), columns = beta_group /
        # beta_depth / model-independent. Caption order must follow the grid, row-major.
        panels = [
            (r'OLS $\beta_{\mathrm{group}}$ (no-PC, the decision axis) versus minAC; coloured triangles below the '
             'axis mark each cohort’s recommended minAC — the non-significant V-trough.'),
            r'OLS $\beta_{\mathrm{depth}}$ (no-PC, a diagnostic — never a selection axis) versus minAC.',
            'variants retained versus the minAC threshold (power cost; log scale; model-independent).',
        ]
        if has_pc:
            panels += [
                (r'OLS $\beta_{\mathrm{group}}$ with ancestry PCs (+PC) versus minAC; triangles mark the +PC model’s '
                 'own recommended minAC. Paired with a on a common scale.'),
                r'OLS $\beta_{\mathrm{depth}}$ with ancestry PCs (+PC) versus minAC. Paired with b on a common scale.',
                (r'decision summary — $\beta_{\mathrm{group}}$ with its 95% CI at each cohort’s own recommended '
                 'minAC, under both models; the vertical line is zero.'),
            ]
        if has_pc:
            # Which cohorts the two models agree on is DERIVED, never asserted. The
            # sentence used to name refined/expanded/full_mainland literally, so any
            # --Cohorts subset produced a caption that contradicted its own figure.
            agree = [lab for lab in labels if rec_pc.get(lab) is not None and rec_pc[lab] == rec_noPC.get(lab)]
            disagree = [lab for lab in labels if rec_pc.get(lab) is not None and rec_pc[lab] != rec_noPC.get(lab)]
            if disagree:
                # Report the divergence as the measured perturbation it is, rather than
                # attributing a mechanism to the +PC model that nothing here computes.
                shifts = []
                for lab in labels:
                    k = rec_noPC.get(lab)
                    r = data[lab][data[lab]['MinAC_Threshold'] == k] if k is not None else None
                    if r is not None and len(r):
                        r = r.iloc[0]
                        b0, b1 = r['Beta_group_noPC'], r['Beta_group_pc']
                        if pd.notna(b0) and pd.notna(b1) and b0 != 0:
                            shifts.append(abs(b1 - b0) / abs(b0))
                stxt = ('; adding the PCs moves it by '
                        f'{min(shifts):.1f}-{max(shifts):.1f}x its own magnitude') if shifts else ''
                dtxt = '; '.join(f'{lab} ({rec_noPC.get(lab)} vs {rec_pc[lab]})' for lab in disagree)
                agreement = (f'the two models agree on the recommended minAC in {len(agree)}/{len(labels)} cohorts '
                             f'({", ".join(agree) if agree else "none"}) and differ in {dtxt}. At the recommended '
                             f'minAC the no-PC coefficient is the quantity being tuned{stxt}, in cohort-dependent '
                             'directions — so the PC-adjusted series is reported as a sensitivity check and does not '
                             'drive the threshold. ')
            else:
                agreement = ('the two models agree on the recommended minAC in every cohort, so the threshold does '
                             'not depend on whether ancestry is adjusted. ')
            interpret = (f'no-PC recommendations: {rec_txt}. ' + agreement +
                         'The near-identical trajectories across nested sample sets show the tuning is robust to '
                         'sample selection; this does NOT remove the platform confounding (cases 30×, controls 15×), '
                         'which only depth-matched down-sampling and the genome-wide inflation factor λ address.')
            notes = (f'{Path(args.stats_combined).name} — per-minAC OLS on the burden rate, one line per cohort. '
                     f'Top row = no-PC (primary); bottom row = +PC (ancestry sensitivity: {pc_txt}). Filled markers, '
                     f'p < {args.alpha:g}; open markers, non-significant. Triangles mark each model’s own '
                     'recommended minAC. minAC = 0 is omitted (singleton-dominated).')
        else:
            interpret = (f'the near-identical trajectories and shared recommended minAC ({rec_txt}) across nested '
                         'sample sets show the tuning is robust to sample selection. This does NOT remove the platform '
                         'confounding — cases are 30× and controls 15× in every cohort alike; only depth-matched '
                         'down-sampling and the genome-wide inflation factor λ address that residual axis.')
            notes = (f'{Path(args.stats_combined).name} — per-minAC OLS on the burden rate, one line per cohort '
                     f'(primary model: {beta_tag}; ancestry PCs available: {pc_txt}). Filled markers, '
                     f'p < {args.alpha:g}; open markers, non-significant.')
        ps.caption_block(
            fig, top=cap_top, wspace=0.34, hspace=0.52,
            title=('The minAC sweep and its recommended threshold are stable across the three nested cohorts under '
                   'both models — robust to sample selection, though not to the platform confounding.'),
            panels=panels, interpret=interpret, notes=notes,
            defs=['minac', 'rate', 'beta_group', 'beta_depth', 'ci_sig', 'models', 'variants_retained', 'lambda'],
            model=ps.FORMULAS['ols'])
        fig.savefig(args.out_png)
        print(f'Comparison figure : {args.out_png}')

        # ── Per-cohort summary table ─────────────────────────────────────────
        summ = None
        if args.summary_combined:
            summ = pd.read_csv(args.summary_combined, sep='\t')
            scol = _cohort_col(summ)

        rows = []
        for lab in labels:
            df = data[lab]
            row = {'cohort': lab, 'recommended_minac': rec[lab], 'beta_model': beta_tag}
            r0 = df[df['MinAC_Threshold'] == 0]
            if len(r0):
                r0 = r0.iloc[0]
                row['n_samples'] = int(r0['Sample_Count']) if pd.notna(r0.get('Sample_Count')) else None
                row['variant_count_minac0'] = int(r0['Variant_Count']) if pd.notna(r0.get('Variant_Count')) else None
                row['rate_mean_minac0'] = r0.get('Rate_Mean')
                row['beta_group_noPC_minac0'] = r0.get('Beta_group_noPC')
                row['beta_group_pc_minac0'] = r0.get('Beta_group_pc')
                row['beta_depth_noPC_minac0'] = r0.get('Beta_depth_noPC')
            rr = df[df['MinAC_Threshold'] == rec[lab]] if rec[lab] is not None else df.iloc[0:0]
            if len(rr):
                rr = rr.iloc[0]
                row['beta_group_at_rec'] = rr.get(b)
                row['p_group_at_rec'] = rr.get(p)
                row['variant_count_at_rec'] = int(rr['Variant_Count']) if pd.notna(rr.get('Variant_Count')) else None
            if summ is not None:
                srow = summ[summ[scol] == lab]
                if len(srow):
                    for k, v in srow.iloc[0].to_dict().items():
                        if k != scol and k not in row:
                            row[k] = v
            rows.append(row)

        out = pd.DataFrame(rows)
        out.to_csv(args.out_tsv, sep='\t', index=False)
        print(f'Comparison table  : {args.out_tsv}')
        with pd.option_context('display.width', 220, 'display.max_columns', None):
            print(out.to_string(index=False))

    except Exception as e:
        print(f'Error comparing cohorts: {e}', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
