#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Justify WHICH statistic the minAC decision is read from. Two choices
#           are involved and they are two edges of one 2x2:
#             which coefficient  -> beta_group, not beta_depth
#             from which model   -> no-PC, not +PC
#           Both follow from one idea: in the sweep beta_group is a DETECTOR for
#           the technical artifact, not an effect estimate to be de-biased. The
#           design gives cases and controls zero platform overlap, so the artifact
#           IS the case/control contrast; the decision statistic must therefore be
#           the one that tracks it.
#             (a) the premise    — platform x group: zero overlap
#             (b) beta_group vs beta_depth — their minima are in different places,
#                 and beta_depth's sits deep inside the artifact's climb
#             (c) no-PC vs +PC   — adjusting moves the reading by more than the
#                 reading itself, in cohort-dependent directions
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
from matplotlib.patches import FancyBboxPatch

import plot_style as ps

COHORT_PALETTE = [ps.RHO_COLOR, ps.BETA_COLOR, ps.GROUP_COLOR, '#CC79A7', '#E69F00']


def parse_args():
    p = argparse.ArgumentParser(description='Why the minAC decision reads beta_group from the no-PC model.')
    p.add_argument('--stats-combined', required=True, help='all_cohorts_qc_stats.tsv (Cohort column).')
    p.add_argument('--platform-file', required=True, help='Per-sample platform table (csv/xlsx).')
    p.add_argument('--platform-id-col', default='ID_JHRPv6')
    p.add_argument('--platform-col', default='WGS_Platform')
    p.add_argument('--group-col', default='Outcome', help='Case/control column in the platform table.')
    p.add_argument('--case-value', default='PH')
    p.add_argument('--pc-template', default=None,
                   help='Ancestry PC table path with @@COHORT@@; R2 is reported for every cohort '
                        'because it is not the same in all of them.')
    p.add_argument('--n-pcs', type=int, default=None,
                   help='Use only the first N PCs — must match the +PC model, or the reported '
                        'R2 describes covariates the model never adjusts on.')
    p.add_argument('--reference-cohort', default=None,
                   help='Cohort for panels a and b; defaults to the largest, which is the most '
                        'complete view of the design (row order in the combined table is not '
                        'guaranteed, so do not rely on it).')
    p.add_argument('--alpha', type=float, default=0.05)
    p.add_argument('--out-png', required=True)
    p.add_argument('--out-tsv', required=True)
    return p.parse_args()


def _read_table(path, id_col):
    p = str(path)
    if p.lower().endswith(('.xlsx', '.xls')):
        df = pd.read_excel(p, dtype={id_col: str})
    else:
        df = pd.read_csv(p, dtype={id_col: str})
    df[id_col] = df[id_col].astype(str).str.strip()
    return df


def _r2_case_given_pcs(pc_df, case, n_pcs=None):
    """Fraction of the case/control contrast that the ancestry PCs span.

    This is the quantity that matters for the model choice: whatever the PCs can
    explain of case/control is, in this design, artifact they would absorb. Only
    the first `n_pcs` are used — the covariate file carries 20 but the +PC model
    fits 10, and quoting the 20-PC figure would describe a model nobody runs.
    """
    pcs = [c for c in pc_df.columns if c.upper().startswith('PC')]
    pcs = sorted(pcs, key=lambda c: int(''.join(ch for ch in c if ch.isdigit()) or 0))
    if n_pcs:
        pcs = pcs[:n_pcs]
    if not pcs or case.nunique() < 2:
        return np.nan, len(pcs)
    X = np.column_stack([np.ones(len(pc_df))] + [pc_df[c].astype(float).values for c in pcs])
    y = case.astype(float).values
    resid = np.linalg.lstsq(X, y, rcond=None)[1]
    if not len(resid):
        return np.nan, len(pcs)
    return float(1 - resid[0] / ((y - y.mean()) ** 2).sum()), len(pcs)


def main():
    args = parse_args()
    ps.setup_style('slide')

    stats = pd.read_csv(args.stats_combined, sep='\t')
    ccol = next(c for c in ('Cohort', 'cohort', 'COHORT') if c in stats.columns)
    labels = list(dict.fromkeys(stats[ccol].tolist()))
    data = {l: stats[stats[ccol] == l].sort_values('MinAC_Threshold') for l in labels}
    colors = {l: COHORT_PALETTE[i % len(COHORT_PALETTE)] for i, l in enumerate(labels)}
    if args.reference_cohort in labels:
        ref = args.reference_cohort
    else:                       # largest cohort — deterministic, unlike table row order
        ref = max(labels, key=lambda l: float(data[l]['Sample_Count'].dropna().max() or 0))

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(17.5, 6.6))
    ledger = []

    # ── (a) the premise: the artifact IS the case/control contrast ────────────
    def _pc_path(cohort):
        return Path(args.pc_template.replace('@@COHORT@@', cohort)) if args.pc_template else None

    def _read_pc(cohort):
        p = _pc_path(cohort)
        if p is None or not p.exists():
            return None
        t = pd.read_csv(p, sep='\t')
        t.columns = [c.lstrip('#') for c in t.columns]
        t['IID'] = t['IID'].astype(str).str.strip()
        return t

    plat = _read_table(args.platform_file, args.platform_id_col)
    pc = _read_pc(ref)
    keep = set(pc['IID']) if pc is not None else None
    sub = plat[plat[args.platform_id_col].isin(keep)] if keep else plat
    sub = sub[sub[args.platform_col].notna()].copy()
    sub['is_case'] = (sub[args.group_col].astype(str) == args.case_value)

    tab = (sub.groupby(args.platform_col)['is_case']
              .agg(case='sum', n='size').reset_index().sort_values('n'))
    tab['ctrl'] = tab['n'] - tab['case']
    # A SCHEMATIC, not a chart. This panel states the premise the other two rest on;
    # drawn as bars it invited the reader to compare platform sizes, which is not the
    # point. Two disjoint boxes say "no shared platform" in one look.
    n_mixed = int(((tab['case'] > 0) & (tab['ctrl'] > 0)).sum())
    case_plats = tab[tab['case'] > 0].sort_values('n', ascending=False)
    ctrl_plats = tab[tab['ctrl'] > 0].sort_values('n', ascending=False)
    axA.set_xlim(0, 1); axA.set_ylim(0, 1); axA.axis('off')

    BOX_TOP, BOX_BOT = 0.80, 0.34

    def _box(x0, w, plats, col, title, n_tot):
        axA.add_patch(FancyBboxPatch((x0, BOX_BOT), w, BOX_TOP - BOX_BOT,
                                     boxstyle='round,pad=0.012', linewidth=2.0,
                                     edgecolor=col, facecolor=col + '18',
                                     transform=axA.transAxes, clip_on=False))
        axA.text(x0 + w / 2, BOX_TOP + 0.055, f'{title}  (n = {n_tot:,})', ha='center',
                 va='bottom', fontsize=12.5, fontweight='bold', color=col,
                 transform=axA.transAxes)
        names = '\n'.join(str(r[args.platform_col]) for _, r in plats.iterrows())
        axA.text(x0 + w / 2, (BOX_TOP + BOX_BOT) / 2, names, ha='center', va='center',
                 fontsize=10.5, color='#333333', transform=axA.transAxes, linespacing=1.7)

    _box(0.02, 0.45, ctrl_plats, ps.GROUP_COLORS['Control'], 'Controls',
         int(ctrl_plats['ctrl'].sum()))
    _box(0.53, 0.45, case_plats, ps.GROUP_COLORS['Case'], 'Cases',
         int(case_plats['case'].sum()))
    axA.annotate('', xy=(0.5, 0.15), xytext=(0.5, BOX_BOT - 0.03), xycoords=axA.transAxes,
                 arrowprops=dict(arrowstyle='-|>', color='#444444', lw=2.0))
    axA.text(0.5, 0.10, 'not one platform in common', ha='center', va='top', fontsize=11.5,
             color='#333333', transform=axA.transAxes)
    ps.panel_tag(axA, 'a', 'Why it matters')
    # The verdict, not the observation. Every panel of this figure ends in a plain
    # statement, because a reader should not have to derive the conclusion from the
    # evidence it is standing on.
    axA.text(0.5, -0.19, 'no platform is shared between the groups →\n'
                         'the artifact IS the case/control contrast',
             transform=axA.transAxes, ha='center', va='top', fontsize=11,
             color='#111111', fontweight='bold', linespacing=1.4)
    ps.despine(axA, grid_axis='x')
    ledger.append(('premise', 'platforms_shared_by_both_groups', n_mixed))
    ledger.append(('premise', 'platforms_total', len(tab)))

    # ── (b) why beta_group and not beta_depth ────────────────────────────────
    # The two coefficients live on scales that differ by ~20x, so a twin axis makes
    # them look like crossing curves and hides the only thing that matters here:
    # WHERE each one bottoms out. Normalise each by its own maximum — one axis, and
    # the minima become directly comparable.
    d = data[ref]
    d = d[d['MinAC_Threshold'] >= 1]
    x = d['MinAC_Threshold']
    # One colour convention across the whole figure: BETA_COLOR is always "the
    # statistic we actually use", grey/green is always the one we reject. Direct
    # end-labels instead of a legend — a shared legend here had to disambiguate
    # itself with "(b)"/"(c)" suffixes, which is a sign the encoding was overloaded.
    bg, bd = d['Beta_group_noPC'].abs(), d['Beta_depth_noPC'].abs()
    axB.plot(x, bg / bg.max(), color=ps.BETA_COLOR, lw=2.4, zorder=3)
    axB.plot(x, bd / bd.max(), color='#9A9A9A', lw=2.2, zorder=3)
    # Both curves can finish at a similar height, so nudge the end labels apart rather
    # than letting them stack.
    yg_end, yd_end = (bg / bg.max()).iloc[-1], (bd / bd.max()).iloc[-1]
    if abs(yg_end - yd_end) < 0.13:
        mid = (yg_end + yd_end) / 2
        yg_end, yd_end = (mid + 0.075, mid - 0.075) if yg_end >= yd_end else (mid - 0.075, mid + 0.075)
    axB.text(x.iloc[-1] + 0.4, yg_end, r'  $|\beta_{\mathrm{group}}|$',
             color=ps.BETA_COLOR, fontweight='bold', fontsize=12, va='center', ha='left')
    axB.text(x.iloc[-1] + 0.4, yd_end, r'  $|\beta_{\mathrm{depth}}|$',
             color='#6E6E6E', fontweight='bold', fontsize=12, va='center', ha='left')

    kg = int(d.loc[bg.idxmin(), 'MinAC_Threshold'])
    kd = int(d.loc[bd.idxmin(), 'MinAC_Threshold'])
    for k, col in ((kg, ps.BETA_COLOR), (kd, '#9A9A9A')):
        axB.axvline(k, color=col, linestyle=':', lw=1.8, zorder=2)
    axB.annotate(f'κ={kg}', xy=(kg, 0.34), xycoords=('data', 'axes fraction'),
                 xytext=(6, 0), textcoords='offset points', ha='left', va='center',
                 color=ps.BETA_COLOR, fontweight='bold', fontsize=11)
    axB.annotate(f'κ={kd}', xy=(kd, 0.34), xycoords=('data', 'axes fraction'),
                 xytext=(6, 0), textcoords='offset points', ha='left', va='center',
                 color='#6E6E6E', fontweight='bold', fontsize=11)
    axB.set_xlim(x.min() - 0.5, x.max() + 4.2)      # room for the end labels
    axB.set_ylim(-0.04, 1.12)
    axB.set_ylabel('|β|  relative to its own maximum')
    axB.set_xlabel('minAC threshold (minAC)')
    ps.sparse_int_ticks(axB, sorted(x.unique()))
    ps.panel_tag(axB, 'b', 'Read $\\beta_{\\mathrm{group}}$, not $\\beta_{\\mathrm{depth}}$')
    ps.despine(axB)

    # cost of selecting on beta_depth, in every cohort (ledger); the panel shows `ref`
    ref_amp = np.nan
    for lab in labels:
        dd = data[lab]
        dd = dd[dd['MinAC_Threshold'] >= 1]
        k_g = int(dd.loc[dd['Beta_group_noPC'].abs().idxmin(), 'MinAC_Threshold'])
        k_d = int(dd.loc[dd['Beta_depth_noPC'].abs().idxmin(), 'MinAC_Threshold'])
        b_at_kg = float(dd.loc[dd['MinAC_Threshold'] == k_g, 'Beta_group_noPC'].iloc[0])
        row = dd.loc[dd['MinAC_Threshold'] == k_d].iloc[0]
        amp = abs(row['Beta_group_noPC'] / b_at_kg) if b_at_kg else np.nan
        if lab == ref:
            ref_amp = amp
        ledger += [('beta_choice', f'{lab}.argmin_abs_beta_group', k_g),
                   ('beta_choice', f'{lab}.argmin_abs_beta_depth', k_d),
                   ('beta_choice', f'{lab}.beta_group_at_argmin_beta_depth', round(float(row['Beta_group_noPC']), 6)),
                   ('beta_choice', f'{lab}.p_group_at_argmin_beta_depth', float(row['P_group_noPC'])),
                   ('beta_choice', f'{lab}.amplification_if_selected_on_beta_depth', round(float(amp), 1))]
    # One cohort, so the per-cohort table that used to sit in a box here collapses into
    # the verdict. Two annotations saying the same thing was most of the clutter.
    axB.text(0.5, -0.19, f'$|\\beta_{{\\mathrm{{depth}}}}|$ bottoms out at κ={kd}, where the apparent\n'
                         f'effect is {ref_amp:.0f}× larger → it cannot find the trough',
             transform=axB.transAxes, ha='center', va='top', fontsize=11,
             color='#111111', fontweight='bold', linespacing=1.4)

    # ── (c) why no-PC and not +PC ────────────────────────────────────────────
    # Same single cohort as panel b. Overlaying all three here forced cohort colours
    # into a panel that already encodes the model, so orange meant "beta_group" in b
    # and "intermediate_mainland" in c. The cross-cohort view is cohort_comparison.png;
    # the cross-cohort numbers are in the verdict line and the TSV.
    zoom = 6
    shifts = []
    dz = data[ref]
    dz = dz[(dz['MinAC_Threshold'] >= 1) & (dz['MinAC_Threshold'] <= zoom)]
    axC.plot(dz['MinAC_Threshold'], dz['Beta_group_noPC'], color=ps.BETA_COLOR, lw=2.6,
             marker='o', markersize=6, markeredgecolor='white', markeredgewidth=0.7, zorder=4)
    axC.plot(dz['MinAC_Threshold'], dz['Beta_group_pc'], color=ps.GROUP_COLOR, lw=2.2,
             linestyle=(0, (4, 2.5)), marker='^', markersize=6, markerfacecolor='white',
             markeredgecolor=ps.GROUP_COLOR, zorder=4)
    axC.text(zoom + 0.15, dz['Beta_group_noPC'].iloc[-1], '  no-PC', color=ps.BETA_COLOR,
             fontweight='bold', fontsize=12, va='center', ha='left')
    axC.text(zoom + 0.15, dz['Beta_group_pc'].iloc[-1], '  +PC', color=ps.GROUP_COLOR,
             fontweight='bold', fontsize=12, va='center', ha='left')
    axC.axhline(0, color='#888888', lw=1.1, zorder=1)
    # the shift at the decision point, drawn once for the reference cohort
    rr2 = dz[dz['MinAC_Threshold'] == 2]
    if len(rr2):
        rr2 = rr2.iloc[0]
        axC.annotate('', xy=(2, rr2['Beta_group_pc']), xytext=(2, rr2['Beta_group_noPC']),
                     arrowprops=dict(arrowstyle='-|>', color='#444444', lw=2.0,
                                     shrinkA=4, shrinkB=4))
    for lab in labels:
        r = data[lab][data[lab]['MinAC_Threshold'] == 2]
        if not len(r):
            continue
        r = r.iloc[0]
        b0, b1 = r['Beta_group_noPC'], r.get('Beta_group_pc')
        if pd.isna(b1) or b0 == 0:
            continue
        shifts.append(abs(b1 - b0) / abs(b0))
        ledger += [('model_choice', f'{lab}.beta_group_noPC_at_k2', round(float(b0), 6)),
                   ('model_choice', f'{lab}.beta_group_pc_at_k2', round(float(b1), 6)),
                   ('model_choice', f'{lab}.shift_at_k2', round(float(b1 - b0), 6)),
                   ('model_choice', f'{lab}.shift_over_abs_beta', round(float(abs(b1 - b0) / abs(b0)), 2)),
                   ('model_choice', f'{lab}.sign_flipped_at_k2', bool(b0 * b1 < 0))]
    axC.set_xlim(0.6, zoom + 1.5)      # room for the end labels
    axC.set_xticks(list(range(1, zoom + 1)))
    axC.set_xlabel('minAC threshold (minAC)')
    axC.set_ylabel(r'OLS $\beta_{\mathrm{group}}$  (Case $-$ Control)')
    ps.panel_tag(axC, 'c', 'Read it from no-PC, not +PC')
    n_flip = sum(1 for lab in labels
                 for r in [data[lab][data[lab]['MinAC_Threshold'] == 2]]
                 if len(r) and pd.notna(r.iloc[0].get('Beta_group_pc'))
                 and r.iloc[0]['Beta_group_noPC'] * r.iloc[0]['Beta_group_pc'] < 0)
    if shifts:
        ref_shift = (abs(rr2['Beta_group_pc'] - rr2['Beta_group_noPC'])
                     / abs(rr2['Beta_group_noPC'])) if len(dz[dz['MinAC_Threshold'] == 2]) else np.nan
        axC.text(0.5, -0.19, f'+PC moves the reading by {ref_shift:.1f}× its own size\n'
                             f'and flips its sign in {n_flip} of {len(labels)} cohorts',
                 transform=axC.transAxes, ha='center', va='top', fontsize=11,
                 color='#111111', fontweight='bold', linespacing=1.4)
    ps.despine(axC)

    # R^2(case/control | PCs) per cohort: how much of the contrast the PCs span.
    # Reported for every cohort because it is NOT the same in all of them.
    r2_txt, r2_vals, npc_used = '', {}, None
    plat_grp = plat[[args.platform_id_col, args.group_col]].copy()
    plat_grp['is_case'] = plat_grp[args.group_col].astype(str) == args.case_value
    for lab in labels:
        t = _read_pc(lab)
        if t is None:
            continue
        m = t.merge(plat_grp[[args.platform_id_col, 'is_case']], left_on='IID',
                    right_on=args.platform_id_col, how='inner')
        r2, npc = _r2_case_given_pcs(m, m['is_case'], args.n_pcs)
        if np.isfinite(r2):
            r2_vals[lab], npc_used = r2, npc
            ledger.append(('model_choice', f'{lab}.r2_case_given_{npc}_pcs', round(r2, 4)))
    if r2_vals:
        r2_txt = (f' The {npc_used} ancestry PCs span R² = {min(r2_vals.values()):.3f}–'
                  f'{max(r2_vals.values()):.3f} of the case/control contrast depending on the cohort, '
                  'so part of the artifact sits inside the covariates an adjustment would use.')

    # No figure legend: every series is labelled at its own line end. A shared legend
    # had to disambiguate itself with "(b)"/"(c)" suffixes because the two panels
    # encode different things — that is exactly the clutter this figure should not have.
    fig.suptitle('The minAC is read from  β$_{group}$,  in the  no-PC  model — '
                 'not β$_{depth}$, not +PC', fontsize=16, y=0.985)

    ps.caption_block(
        # extra_bottom: each panel carries a two-line verdict under its x-label.
        fig, top=0.905, wspace=0.30, left=0.135, extra_bottom=0.055,
        title=('The minAC decision is read from β_group in the no-PC model, because that is the only one of the '
               'four candidate statistics that tracks the technical artifact.'),
        panels=[
            (f'sequencing platform by case/control: {n_mixed} of {len(tab)} platforms are shared between the '
             'groups. With zero overlap, anything that differs between cases and controls is entangled with the '
             'technical axis — so *group-aligned* and *artifact-aligned* mean the same thing here.'),
            (r'$|\beta_{\mathrm{group}}|$ and $|\beta_{\mathrm{depth}}|$ versus minAC for '
             f'{ref}, each scaled by its own maximum so the positions of the two minima can be compared '
             '(dotted lines). They are in different places, and '
             r'$\beta_{\mathrm{depth}}$’s sits well inside the range where the apparent effect is already large.'),
            (r'$\beta_{\mathrm{group}}$ over the decision region for the same cohort under both models '
             '(solid, no-PC; dashed, +PC); the arrow marks the shift at κ = 2. Across all three cohorts the '
             'shift is 1.9–2.6× and reverses sign in two of them — per-cohort values in decision_axis.tsv.'),
        ],
        interpret=(r'$\beta_{\mathrm{group}}$ here is a DETECTOR for the technical artifact, not an effect estimate '
                   'to be de-biased, and a detector is chosen for sensitivity to its target. Selecting on '
                   r'$\beta_{\mathrm{depth}}$ instead would land on a threshold where the apparent case/control '
                   'effect is far larger (b), and adjusting on ancestry PCs moves the reading by more than the '
                   'reading itself, in cohort-dependent directions (c) — so neither can locate the trough.' + r2_txt +
                   ' The downstream association test keeps its PCs: it needs an estimator, not a detector — a '
                   'different job, so a different model.'),
        defs=['minac', 'beta_group', 'beta_depth', 'models', 'ci_sig'],
        notes=(f'{Path(args.stats_combined).name} (per-minAC OLS on the burden rate) and '
               f'{Path(args.platform_file).name} (platform, case/control). Panel a is restricted to the '
               f'{ref} samples. Panel b uses {ref}; the per-cohort figures are in {Path(args.out_tsv).name}.'),
        model=ps.FORMULAS['ols'])
    fig.savefig(args.out_png)
    print(f'Figure         : {args.out_png}')

    out = pd.DataFrame(ledger, columns=['section', 'key', 'value'])
    out.to_csv(args.out_tsv, sep='\t', index=False)
    print(f'Ledger         : {args.out_tsv}')
    with pd.option_context('display.width', 200, 'display.max_rows', None):
        print(out.to_string(index=False))


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in plot_decision_axis: {e}', file=sys.stderr)
        sys.exit(1)
