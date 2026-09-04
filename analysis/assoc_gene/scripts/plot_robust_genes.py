#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The main figure. Every gene that reaches a robustness tier, one row
#           each, read across three row-aligned panels:
#
#             (a) the evidence -- its state in every cohort x statistic cell of
#                 both variant sets, and the tier each set awards it;
#             (b) the effect -- the CMC burden odds ratio with its 95 % CI in the
#                 report cohort, both variant sets;
#             (c) the people -- the fraction of controls and of cases carrying
#                 a qualifying allele, with the counts.
#
#           Rows are grouped by the gene's best tier, so the reader meets the
#           strongest claims first and can see, in the same row, what the claim
#           rests on: how many cells called it, how large the effect is, and how
#           many carriers produced it. A tier is a summary of (a); (b) and (c)
#           are what the summary table prints, drawn.
#
#           Four states are kept strictly apart in (a), because conflating them
#           overstates a result: Bonferroni + BH, BH only, tested but not called,
#           and NOT TESTED in that cell (fewer than MinNumVar variants in that
#           cohort's map). A blank usually means "no effect"; here it can mean
#           "never asked", and the two must not look alike.
#
#           THE HEIGHT IS COMPUTED, not fixed: a fixed inches-per-row keeps rows
#           readable whether three genes reach a tier or forty. The figure grows.
#
#           This is NOT replication. The cohorts are nested and differ almost
#           entirely in controls; agreement across them shows a signal survives
#           adding controls, and the README says so.
# Component: assoc_gene
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
from matplotlib.ticker import FixedLocator, ScalarFormatter   # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_shared' / 'scripts'))
import plot_style as S                   # noqa: E402
import figure_doc                        # noqa: E402
import vocab as V                        # noqa: E402

ROW_IN = 0.21              # inches per gene row
LEGEND_IN = 0.40           # (a)'s state key, and (b)+(c)'s two-row shared key
XLAB_IN = 0.48             # caption_block's x-label slot (rotated bottom ticks live there)
BLOCK_GAP = 1.4            # x units between the two stratum blocks in (a): wide enough
                           # that the two full stratum names cannot meet over them
# The tier is an ordered summary, so it takes a grey ramp -- the cell states own
# the red/blue, and two colour channels for two different things would compete.
TIER_FILL = {1: S.INK, 2: S.INK_SOFT, 3: S.NEUTRAL_D}
# The band that groups a tier's rows. grouped_row_labels' own banding is NEUTRAL
# at alpha 0.13 on alternate groups only, which composites to (249,250,250) on
# white -- six levels out of 255, invisible in print, and with three tiers it
# banded Tier 2 alone. Here every group gets one, over NEUTRAL_D and at an alpha
# that deepens with the tier, so the band joins the bar, the group label and the
# boxed number in ONE grey channel: (213,216,218) / (229,231,232) / (243,244,245).
TIER_BAND = {1: 0.42, 2: 0.26, 3: 0.12}
TIER_BAR_W = 2.6           # points: the solid tier bar at each panel's left edge
STRATUM_H = 0.15          # the stratum name over each block of (a)
# The row count is 3x the gene count now, so the canvas is capped and the
# overflow reported rather than letting it run past verify.sh's ceiling.
MAX_PLOT_H = 11.3
OR_RATIO_MAX = 1e4         # (b): a wider 95 % CI than this is not an estimate
DUMBBELL_FRAC = 0.50       # of (c)'s width; the counts take the rest
COUNT_X = 0.72             # the counts' centre, in (c)'s axes fraction
METHODS_REF = '../../../docs/METHODS.md'


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--evidence', required=True, help='evidence.tsv from ROBUSTNESS')
    p.add_argument('--tiers', required=True, help='tiers.tsv from ROBUSTNESS')
    p.add_argument('--cohort-samples', required=True, help='cohort_samples.tsv')
    p.add_argument('--denominators', required=True, help='denominators.tsv (README table)')
    p.add_argument('--report-cohort', required=True,
                   help='the cohort whose OR and carrier frequencies fill (b) and (c)')
    p.add_argument('--max-genes', type=int, default=40)
    p.add_argument('--alpha', type=float, default=0.05,
                   help='the alpha the scan called at; quoted in the caption\'s tier rule')
    p.add_argument('--out-png', default='robust_genes.png')
    p.add_argument('--out-md', default='README.md')
    return p.parse_args()


def fnum(v):
    x = pd.to_numeric(v, errors='coerce')
    return float(x) if np.isfinite(x) else None


def main():
    a = parse_args()
    S.setup_style('paper')
    ev = pd.read_csv(a.evidence, sep='\t', dtype=str)
    tiers = pd.read_csv(a.tiers, sep='\t', dtype=str)
    smp = pd.read_csv(a.cohort_samples, sep='\t', dtype=str)
    den = pd.read_csv(a.denominators, sep='\t', dtype=str)
    if ev.empty:
        raise SystemExit(f'ABORT: {a.evidence} is empty')
    bad = set(ev['state']) - set(V.STATES)
    if bad:
        raise SystemExit(f'ABORT: {a.evidence} carries state(s) {sorted(bad)}')
    n_of = {r['cohort']: (int(r['n_case']), int(r['n_control'])) for _, r in smp.iterrows()}
    if a.report_cohort not in n_of:
        raise SystemExit(f'ABORT: --report-cohort {a.report_cohort!r} is not in {a.cohort_samples}')
    tiers['tier'] = pd.to_numeric(tiers['tier'], errors='coerce').astype(int)

    cohorts = V.order_cohorts(set(ev['cohort']))
    clabel = V.cohort_labels(cohorts)
    # One cohort, two channels: (b) has its colour free and takes the ramp, (c)
    # has spent hue on case/control and takes the size. Both run narrow -> full.
    coh_col = V.cohort_colors(cohorts)
    strata = V.order_strata(ev)
    narrow = strata[-1]                      # the narrowest set: filled markers in (b)(c)
    rc = a.report_cohort

    # ---- rows: one per gene, best tier first, then the tiers file's order.
    best, first_seen = {}, {}
    for i, r in tiers.iterrows():
        sn, t = r['set_name'], int(r['tier'])
        first_seen.setdefault(sn, i)
        if sn not in best or t < best[sn][0]:
            best[sn] = (t, r['gene_symbol'], r['chrom'])
    order = sorted(best, key=lambda sn: (best[sn][0], first_seen[sn]))
    genes = order[:a.max_genes]
    dropped = order[a.max_genes:]
    tier_of = {(r['set_name'], r['stratum']): int(r['tier']) for _, r in tiers.iterrows()}
    ev_idx = ev.set_index(['set_name', 'stratum', 'cohort', 'method'])

    def cell_row(sn, s, c, m):
        try:
            r = ev_idx.loc[(sn, s, c, m)]
        except KeyError:
            return None
        return r.iloc[0] if isinstance(r, pd.DataFrame) else r

    # ---- the empty case keeps the DAG honest.
    if not genes:
        fig = plt.figure(figsize=(S.COL_DOUBLE, 2.0))
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, 'no gene reached a robustness tier', ha='center', va='center',
                transform=ax.transAxes, color=S.INK_SOFT)
        ax.set_axis_off()
        S.caption_block(fig, plot_h=2.0, title='No gene reaches a robustness tier.',
                        panels=['No tiered gene.'], notes='README.md', margin_axes=[ax])
        fig.savefig(a.out_png)
        plt.close(fig)
        figure_doc.write_doc(a.out_png, out_path=a.out_md,
                             subject=f'**Figure file:** `{Path(a.out_png).name}`',
                             title='Robust genes', question='Which genes reach a tier?',
                             panels=[('a', 'None', 'No gene was called in two or more cohorts.')],
                             methods_ref=METHODS_REF)
        print('[plot_robust_genes] no tiered gene')
        return

    # ---- geometry: ONE ROW PER (GENE x COHORT).
    # The cohort used to be a column of the matrix, which forced (b) and (c) to
    # pack three estimates into one gene row -- three overlapping intervals in a
    # 1.2 in panel, so only the report cohort could carry one. Making the cohort
    # a ROW gives every cohort its own interval and its own denominator, and it
    # is the idiom the project already uses for cross-cohort figures:
    # cohort_manhattan draws one cohort per row, and plot_cohort_compare's
    # forests repeat every locus once per cohort through grouped_row_labels.
    fig = plt.figure(figsize=(S.COL_DOUBLE, 4.0))
    fs_tick = plt.rcParams['ytick.labelsize']
    header_in = STRATUM_H + 0.10
    top_pad = 0.24 + LEGEND_IN + header_in
    # The row count is now 3x the gene count, so the drawn list is trimmed to
    # what fits rather than letting the canvas run past what verify.sh calls a
    # plausible figure. The genes trimmed here join the ones --max-genes already
    # dropped, and the caption says how many.
    per_gene_in = len(cohorts) * ROW_IN
    n_fit = max(1, int((MAX_PLOT_H - top_pad - XLAB_IN) // per_gene_in))
    if len(genes) > n_fit:
        dropped = list(dropped) + genes[n_fit:]
        genes = genes[:n_fit]
    n_rows = len(genes) * len(cohorts)
    plot_h = top_pad + n_rows * ROW_IN + XLAB_IN
    fig.set_size_inches(S.COL_DOUBLE, plot_h)
    # (a) needs only 2 methods x 2 variant sets + 2 tier columns now that the
    # cohort is a row, so it gives most of its old width to (b) and (c).
    gs = fig.add_gridspec(1, 3, width_ratios=[1.45, 1.0, 1.0])
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0, i]) for i in range(3))

    # ---- row labels: the gene names its block once, the ticks carry the cohort.
    groups = []
    for sn in genes:
        groups.append(([
            (best[sn][1], dict(fontstyle='italic', fontweight='bold', fontsize=fs_tick)),
            (V.TIER_LABEL[best[sn][0]], dict(fontsize=fs_tick - 0.5, color=S.INK_SOFT)),
        ], [clabel[c] for c in cohorts]))
    ticks = S.grouped_row_labels(ax_a, groups, band=False)
    y_of = {(sn, c): ticks[i * len(cohorts) + j]
            for i, sn in enumerate(genes) for j, c in enumerate(cohorts)}
    ylim = ax_a.get_ylim()

    # ---- the tier bands, one per gene block, on all three panels so a row reads
    # across them. Two marks, because a wash alone did not read: a solid bar in
    # (a)'s right gutter in the tier's own grey (the encoding) and a band behind
    # the block (the grouping), with a hairline between blocks.
    y = n_rows - 1
    for gi, sn in enumerate(genes):
        tier = best[sn][0]
        top, bottom = y, y - len(cohorts) + 1
        for ax in (ax_a, ax_b, ax_c):
            ax.axhspan(bottom - 0.5, top + 0.5, color=S.NEUTRAL_D, alpha=TIER_BAND[tier],
                       linewidth=0, zorder=0)
            if gi:
                ax.axhline(top + 0.5, color=S.NEUTRAL_D, lw=0.5, zorder=1)
        ax_a.plot([1.0, 1.0], [bottom - 0.30, top + 0.30], transform=ax_a.get_yaxis_transform(),
                  color=TIER_FILL[tier], lw=TIER_BAR_W, solid_capstyle='butt',
                  clip_on=False, zorder=6)
        y -= len(cohorts)

    # ---- (a) the evidence matrix: 2 statistics x 2 variant sets, + the tier.
    # Which of the four states the matrix actually draws. A key for a mark that
    # is not on the page sends the reader looking for it -- and `not_in_map` is
    # the one that usually is not there, because a gene absent from a cohort's
    # map cannot be called in it and so rarely reaches a tier at all.
    states_drawn = set()
    cols = []                                # (x, stratum, method | 'tier')
    x = 0.0
    for s in strata:
        for m in V.METHODS:
            cols.append((x, s, m)); x += 1.0
        cols.append((x, s, 'tier')); x += 1.0
        x += BLOCK_GAP
    x_last = cols[-1][0]
    for sn in genes:
        # The tier is a property of (gene, variant set), not of a cohort, so its
        # box is drawn ONCE per block, centred on the three cohort rows.
        y_mid = (y_of[(sn, cohorts[0])] + y_of[(sn, cohorts[-1])]) / 2.0
        for xx, s, m in cols:
            if m == 'tier':
                t = tier_of.get((sn, s))
                tested = any(cell_row(sn, s, c2, 'cmc') is not None
                             and cell_row(sn, s, c2, 'cmc')['state'] != 'not_in_map'
                             for c2 in cohorts)
                # Its own channel: a grey box that deepens with the tier, never
                # the red/blue the cell states own. 0 and — stay plain text,
                # because neither is a tier.
                if t:
                    ax_a.text(xx, y_mid, str(t), ha='center', va='center', fontsize=fs_tick,
                              fontweight='bold', color='white', zorder=5,
                              bbox=dict(boxstyle='square,pad=0.30', facecolor=TIER_FILL[t],
                                        edgecolor='none'))
                elif tested:
                    ax_a.text(xx, y_mid, '0', ha='center', va='center', fontsize=fs_tick,
                              color=S.INK_SOFT, zorder=4)
                else:
                    ax_a.text(xx, y_mid, '—', ha='center', va='center', fontsize=fs_tick,
                              color=S.REFERENCE, zorder=4)
                continue
            for c in cohorts:
                r = cell_row(sn, s, c, m)
                st = r['state'] if r is not None else 'not_in_map'
                states_drawn.add(st)
                V.draw_state(ax_a, xx, y_of[(sn, c)], st)
    # Structure: a rule between the two stratum blocks.
    per_block = len(V.METHODS) + 1
    for bi in range(1, len(strata)):
        xb = cols[bi * per_block][0] - (BLOCK_GAP + 1.0) / 2.0
        ax_a.axvline(xb, color=S.INK_SOFT, lw=0.8, zorder=1)
    ax_a.set_xlim(-0.6, x_last + 0.6)
    ax_a.set_xticks([xx for xx, *_ in cols])
    ax_a.set_xticklabels([V.METHOD_SHORT[m] if m != 'tier' else 'tier' for _x, _s, m in cols],
                         rotation=90, fontsize=fs_tick)
    ax_a.tick_params(axis='x', length=0)
    S.despine(ax_a, grid_axis='none')
    for sp in ('left', 'bottom'):
        ax_a.spines[sp].set_visible(False)

    # ---- (b) CMC OR with its 95 % CI, ONE COHORT PER ROW.
    # Every row now carries a full interval: with the cohort on the row there is
    # nothing to overlap, so the compromise of drawing only the report cohort's
    # interval is gone. Colour still follows the cohort (COHORT_RAMP, the house
    # grammar in plot_cohort_compare), which makes the ramp redundant with the
    # row label -- reinforcement, not a second variable.
    lo_all, hi_all, drawn_b = [], [], []
    for sn in genes:
        for c in cohorts:
            yy, col = y_of[(sn, c)], coh_col[c]
            for s in strata:
                dy = 0.17 if s == narrow else -0.17
                filled = s == narrow
                r = cell_row(sn, s, c, 'cmc')
                if r is None or r['state'] == 'not_in_map':
                    # THE SAME × AS (a). It used to be a '—' drawn at OR = 1, a
                    # data position that reads as an estimate of one, while the ×
                    # was spent on a different fact below. One rule now: the plot
                    # area holds estimates, the margin holds the reason there is
                    # none, and × means "not in that cohort's map" in both panels.
                    ax_b.plot([1.03], [yy + dy], transform=ax_b.get_yaxis_transform(),
                              marker='x', ls='none', markersize=3.2,
                              markeredgecolor=S.REFERENCE, markeredgewidth=1.2, zorder=5,
                              clip_on=False)
                    continue
                o, lo, hi = fnum(r['or']), fnum(r['or_l95']), fnum(r['or_u95'])
                nc, nk = fnum(r['n_carrier_case']), fnum(r['n_carrier_control'])
                # rvtest's Wald fit runs away when one group holds no minor
                # allele: an OR of order 1e5 with an interval spanning 1e-56 to
                # 1e66, which would collapse every real interval on a shared log
                # axis. Such a cell is marked in the margin, never drawn, and
                # never enters the limits.
                estimable = (o is not None and lo is not None and hi is not None and lo > 0
                             and nc and nk and hi / lo <= OR_RATIO_MAX)
                if not estimable:
                    # NOT an ×. The cell WAS tested -- it has a P, and one of the
                    # two here is even BH-called -- so borrowing the "not in map"
                    # glyph said the opposite of what (a) says on the same row.
                    # What is missing is only the effect estimate.
                    ax_b.text(1.03, yy + dy, 'n.e.', transform=ax_b.get_yaxis_transform(),
                              ha='left', va='center', fontsize=V.MIN_FONT,
                              color=S.REFERENCE, zorder=5, clip_on=False)
                    sk = cell_row(sn, s, c, 'skato')
                    # n_var_map has to be carried here too. Omitting it left the
                    # column empty for exactly these rows, and the renderer fills
                    # an empty cell with an em dash -- the glyph that means "not
                    # in map", on the two rows most likely to be misread that way.
                    drawn_b.append({'gene_symbol': best[sn][1], 'stratum': s, 'cohort': c,
                                    'or': 'not estimable', 'or_l95': 'NA', 'or_u95': 'NA',
                                    'cmc_state': r['state'], 'cmc_p': r['pvalue'],
                                    'n_var_map': r['n_var_map'],
                                    'skato_p': sk['pvalue'] if sk is not None else 'NA'})
                    continue
                ax_b.plot([lo, hi], [yy + dy, yy + dy], color=col, lw=1.0, zorder=2)
                ax_b.scatter([o], [yy + dy], s=20, marker='o',
                             facecolors=col if filled else 'white', edgecolors=col,
                             linewidths=0.9 if filled else 1.1, zorder=3)
                lo_all.append(lo); hi_all.append(hi)
                sk = cell_row(sn, s, c, 'skato')
                drawn_b.append({'gene_symbol': best[sn][1], 'stratum': s, 'cohort': c,
                                'or': o, 'or_l95': lo, 'or_u95': hi, 'cmc_state': r['state'],
                                'cmc_p': r['pvalue'],
                                'n_var_map': r['n_var_map'],
                                'skato_p': sk['pvalue'] if sk is not None else 'NA'})
    ax_b.axvline(1.0, color=S.REFERENCE, lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax_b.set_xscale('log')
    if lo_all:
        ax_b.set_xlim(min(lo_all) / 1.4, max(hi_all) * 1.4)
    else:
        ax_b.set_xlim(0.3, 5.0)
    S.or_log_axis(ax_b)
    tk = list(ax_b.get_xticks())
    lo_x, hi_x = ax_b.get_xlim()
    if not any(abs(t - 1.0) < 1e-9 for t in tk) and lo_x <= 1.0 <= hi_x:
        tk = sorted(set(tk) | {1.0})
        ax_b.xaxis.set_major_locator(FixedLocator(tk))
        ax_b.xaxis.set_major_formatter(ScalarFormatter())
    ax_b.set_ylim(*ylim)
    ax_b.set_yticks(ticks)
    ax_b.tick_params(axis='y', labelleft=False, length=0)
    ax_b.set_xlabel(f'CMC {S.OR_SYM} (95% CI)')
    S.despine(ax_b, grid_axis='x')

    # ---- (c) carriers: control % -> case %, EACH COHORT AGAINST ITS OWN N.
    # The cohort is the row, so the size ramp pass 4 needed to tell three
    # coincident markers apart is gone -- a second encoding of a variable the
    # row label already names would be noise.
    n_case, n_ctrl = n_of[rc]          # quoted in the caption's tier rule
    xmax, drawn_c = 0.0, []
    for sn in genes:
        s_use = next((s for s in (narrow, *strata)
                      if (r := cell_row(sn, s, rc, 'cmc')) is not None
                      and r['state'] != 'not_in_map'), None)
        if s_use is None:
            continue
        filled = s_use == narrow
        for c in cohorts:
            yy = y_of[(sn, c)]
            r = cell_row(sn, s_use, c, 'cmc')
            if r is None or r['state'] == 'not_in_map':
                ax_c.text(0.0, yy, '—', ha='left', va='center', fontsize=fs_tick,
                          color=S.REFERENCE, zorder=3)
                continue
            n_ca, n_ko = n_of[c]
            nc, nk = int(float(r['n_carrier_case'])), int(float(r['n_carrier_control']))
            pc, pk = 100.0 * nc / n_ca, 100.0 * nk / n_ko
            ax_c.plot([pk, pc], [yy, yy], color=S.NEUTRAL_D, lw=1.0, zorder=2)
            ax_c.scatter([pk], [yy], s=20, marker='o',
                         facecolors=V.CONTROL_COLOR if filled else 'white',
                         edgecolors=V.CONTROL_COLOR, linewidths=0.9, zorder=3)
            ax_c.scatter([pc], [yy], s=20, marker='o',
                         facecolors=V.CASE_COLOR if filled else 'white',
                         edgecolors=V.CASE_COLOR, linewidths=0.9, zorder=4)
            xmax = max(xmax, pc, pk)
            # Counts in a fixed right-hand column, this cohort's own: cases in
            # red, controls in blue, the group sizes in the caption.
            ax_c.text(COUNT_X - 0.04, yy, f'{nc}', transform=ax_c.get_yaxis_transform(),
                      ha='right', va='center', fontsize=V.MIN_FONT, color=V.CASE_COLOR,
                      zorder=5)
            ax_c.text(COUNT_X, yy, '·', transform=ax_c.get_yaxis_transform(), ha='center',
                      va='center', fontsize=V.MIN_FONT, color=S.INK_SOFT, zorder=5)
            ax_c.text(COUNT_X + 0.04, yy, f'{nk}', transform=ax_c.get_yaxis_transform(),
                      ha='left', va='center', fontsize=V.MIN_FONT, color=V.CONTROL_COLOR,
                      zorder=5)
            drawn_c.append({'gene_symbol': best[sn][1], 'stratum': s_use, 'cohort': c,
                            'carriers_case': nc, 'n_case': n_ca, 'pct_case': round(pc, 3),
                            'carriers_control': nk, 'n_control': n_ko,
                            'pct_control': round(pk, 3), 'n_var_map': r['n_var_map']})
    # The dumbbells occupy the left DUMBBELL_FRAC of the panel; the counts the rest.
    ax_c.set_xlim(0, max(xmax, 1.0) / DUMBBELL_FRAC)
    ax_c.set_ylim(*ylim)
    ax_c.set_yticks(ticks)
    ax_c.tick_params(axis='y', labelleft=False, length=0)
    ax_c.set_xlabel('carriers (%)')
    # A header for the count column, in the two group colours, so the pair of
    # numbers on every row needs no repeated label.
    y_head = ylim[1] + 0.02
    ax_c.text(COUNT_X - 0.04, y_head, 'case', transform=ax_c.get_yaxis_transform(),
              ha='right', va='bottom', fontsize=V.MIN_FONT, fontweight='bold',
              color=V.CASE_COLOR, clip_on=False)
    ax_c.text(COUNT_X, y_head, '·', transform=ax_c.get_yaxis_transform(), ha='center',
              va='bottom', fontsize=V.MIN_FONT, color=S.INK_SOFT, clip_on=False)
    ax_c.text(COUNT_X + 0.04, y_head, 'control', transform=ax_c.get_yaxis_transform(),
              ha='left', va='bottom', fontsize=V.MIN_FONT, fontweight='bold',
              color=V.CONTROL_COLOR, clip_on=False)
    S.despine(ax_c, grid_axis='x')

    # ---- caption, then the measured layout.
    n_t = {t: sum(1 for sn in order if best[sn][0] == t) for t in (1, 2, 3)}
    # How often the per-cohort maps disagree: the reason a carrier count can FALL
    # as controls are added (PES1 has 4 MODERATE+HIGH variants in the two smaller
    # cohorts and 3 in the largest), which would otherwise read as an error.
    n_not_est = sum(1 for d in drawn_b if d['or'] == 'not estimable')
    n_cells = len(drawn_b)
    # The caption names a mark only when the figure carries it, for the same
    # reason the legend does.
    unmapped = 'not_in_map' in states_drawn
    marks = ([f'× = not in that cohort\'s map, the same mark as (a)'] if unmapped else [])
    marks += ([f'n.e. = tested, but no {S.OR_SYM} estimate'] if n_not_est else [])
    margin_txt = f' In the margin: {"; ".join(marks)}.' if marks else ''
    tier_x = ' × = not in that cohort\'s map.' if unmapped else ''
    n_var_moves = sum(1 for sn in genes for st in strata
                      if len({(r['n_var_map'] if (r := cell_row(sn, st, c, 'cmc')) is not None
                               else None) for c in cohorts} - {None}) > 1)
    S.caption_block(
        fig, plot_h=plot_h,
        title=(f'{len(order)} gene(s) reach a robustness tier: {n_t[1]} Tier 1, {n_t[2]} Tier 2, '
               f'{n_t[3]} Tier 3.'),
        panels=[
            'Rows are gene × cohort. Cell state per statistic in both variant sets; the boxed '
            'number is that set\'s tier, written once per gene, and the band and bar carry the '
            'gene\'s best tier.',
            f'CMC burden {S.OR_SYM} with its 95% CI, one cohort per row, coloured darkest '
            f'{clabel[cohorts[0]]} to lightest {clabel[cohorts[-1]]}. Filled '
            f'{V.stratum_label(narrow)}, hollow {V.stratum_label(strata[0])}, dashed line '
            f'{S.OR_SYM} = 1.{margin_txt}',
            'Carriers in that row\'s own cohort, against that cohort\'s own N: blue control '
            'to red case, with both counts.',
        ],
        notes=(f'Called = BH q < {a.alpha:g} over that cohort × variant set\'s mapped genes. '
               f'Tier 1 = called in all {len(cohorts)} cohorts and by both statistics in ≥ 1; '
               f'Tier 2 = called in all {len(cohorts)}; Tier 3 = called in ≥ 2; 0 = tested, '
               f'called in < 2.{tier_x} Each cohort has its own minAC, '
               f'so a gene\'s mapped variant set is decided per cohort and its carrier counts '
               f'need not grow with the cohort ({n_var_moves} gene × variant set(s) here map a '
               f'different number of variants in different cohorts). {clabel[rc]}: '
               f'{n_case:,} cases, '
               f'{n_ctrl:,} controls. A cell with no minor allele in one group has an '
               f'unbounded Wald fit and so no {S.OR_SYM}: it keeps its state in (a), because it '
               f'was tested, and reads n.e. in (b) ({n_not_est} of {n_cells} cells here). '
               + (f'{len(dropped)} further gene(s) omitted for space. ' if dropped else '')
               + 'README.md'),
        letters='abcdefghijklmnop',
        top_pad=top_pad, wspace=0.30, left='auto', right=0.985, margin_axes=[ax_a])

    # Two keys, one variable each. The cell states belong to (a) and are keyed
    # above it in two columns so the strip stays inside that panel's width; the
    # variant-set fill is the only encoding (b) and (c) still share and is keyed
    # once across the pair. The cohort needs no key at all now that it is the row.
    # ---- the stratum names over their blocks of (a), MEASURED: the two full
    # names are ~1.14 in each and (a) is ~2.1 in, so they meet unless the block
    # gap holds them apart. Drawn here, after caption_block, because only now are
    # the axes positions final; the font shrinks if they would still touch.
    fs_str = fs_tick
    for _ in range(6):
        boxes = []
        for s_tag in strata:
            block = [c for c in cols if c[1] == s_tag]
            xc = (block[0][0] + block[-1][0]) / 2.0
            w, _h = V.text_extent_in(fig, V.stratum_label(s_tag), fs_str, weight='bold')
            boxes.append((xc, w))
        ax_w = ax_a.get_position().width * fig.get_size_inches()[0]
        per_unit = ax_w / (ax_a.get_xlim()[1] - ax_a.get_xlim()[0])
        gap = abs(boxes[1][0] - boxes[0][0]) * per_unit - (boxes[0][1] + boxes[1][1]) / 2.0
        if gap >= 0.05 or fs_str <= V.MIN_FONT:
            break
        fs_str -= 0.25
    for s_tag in strata:
        block = [c for c in cols if c[1] == s_tag]
        x0, x1 = block[0][0], block[-1][0]
        ax_a.annotate(V.stratum_label(s_tag), xy=((x0 + x1) / 2.0, 1.0),
                      xycoords=('data', 'axes fraction'), xytext=(0, 5.0),
                      textcoords='offset points', ha='center', va='bottom',
                      fontsize=fs_str, fontweight='bold', color=S.INK,
                      annotation_clip=False)
        ax_a.annotate('', xy=(x0 - 0.4, 1.0), xycoords=('data', 'axes fraction'),
                      xytext=(x1 + 0.4, 1.0), textcoords=('data', 'axes fraction'),
                      arrowprops=dict(arrowstyle='-', color=S.INK_SOFT, lw=0.7,
                                      shrinkA=0, shrinkB=0),
                      annotation_clip=False)
    keys = [st for st in V.STATES if st in states_drawn]
    S.legend_above(ax_a, V.state_handles(keys), ncol=2, pad_in=header_in + 0.06)
    # (b) and (c) get their keys back, but a key cannot go above either of them:
    # measured at this type scale the cohort row is 1.99 in and the variant-set
    # row 1.41 in, and each panel is 1.18 in wide. They share the cohort ramp
    # anyway, so ONE strip is anchored to (b) and stretched to (c)'s right edge.
    # (b) and (c) share the variant-set fill and nothing else now that the cohort
    # is a row, so ONE key serves them, centred across the pair: measured at this
    # type scale it is 1.41 in and neither panel alone is wide enough for it.
    span = ((ax_c.get_position().x1 - ax_b.get_position().x0)
            / ax_b.get_position().width)
    h_b = ax_b.get_position().height * fig.get_size_inches()[1]
    leg = ax_b.legend(handles=V.stratum_handles(strata), ncol=len(strata),
                      loc='lower center',
                      bbox_to_anchor=(0.0, 1.0 + (header_in + 0.06) / h_b, span, 0.0),
                      frameon=False, handletextpad=0.4, columnspacing=1.1,
                      borderpad=0.0, borderaxespad=0.0, fontsize=fs_tick)
    leg.set_in_layout(False)
    ax_b.add_artist(leg)
    pad = (header_in + LEGEND_IN) * 72.0 + 9.0
    for ax, letter in ((ax_a, 'a'), (ax_b, 'b'), (ax_c, 'c')):
        S.panel_tag(ax, letter, pad=pad)
    S.thin_tick_labels(ax_b, 'x')
    fig.savefig(a.out_png)
    plt.close(fig)

    # ---- README.
    den_tbl = den[den['cohort'].isin(cohorts) & den['stratum'].isin(strata)][
        ['cohort', 'stratum', 'n_genes_mapped', 'bonferroni_threshold']]
    tier_tbl = tiers[['gene_symbol', 'stratum', 'tier', 'n_cohorts_called', 'cohorts_called',
                      'n_cohorts_both_methods', 'called_bonferroni_anywhere',
                      'cohorts_not_in_map', 'min_p_cmc', 'min_p_skato']]
    figure_doc.write_doc(
        a.out_png, out_path=a.out_md,
        subject=f'**Figure file:** `{Path(a.out_png).name}`',
        title='Robust genes — tiers, effect sizes and carriers',
        question=('Which genes survive a change of statistic and the addition of controls, how '
                  'strong is their burden effect, and how many people carry it?'),
        panels=[
            ('a', 'Evidence in every cell, and the tier per variant set',
             'Rows are gene x cohort: every gene reaching Tier 1-3 is repeated once per cohort, '
             'so that (b) and (c) can give each cohort its own interval and its own '
             'denominator. Within each variant set (left LOW+MODERATE+HIGH, right '
             'MODERATE+HIGH) the two cells are the two statistics. Four states: a full red '
             'circle is above that cell\'s own Bonferroni threshold (and therefore also '
             'BH-called); a lighter red circle is BH-called only; a small blue dot is tested and '
             'called by neither; a grey cross is NOT tested -- the gene has fewer than '
             'MinNumVar variants in that cohort\'s map. The key shows only the states this run '
             'actually draws, so a mark a reader cannot find is never advertised; '
             + ('every cell here was tested, so no cross appears and the key has three states. '
                if not unmapped else 'all four occur here. ')
             + 'The boxed number after each block is the '
             'tier that variant set awards; it is a property of the gene, not of a cohort, so it '
             'is written once per block: 1 = called in every cohort and by both statistics in at '
             'least one; 2 = called in every cohort; 3 = called in at least two cohorts; '
             '0 = tested but called in fewer than two; -- = not tested in any cohort. The row '
             'band and the bar in the right gutter carry the gene\'s BEST tier over the variant '
             'sets, so a gene can sit in a Tier 1 band and carry a 0 in the set that did not '
             'call it. The tier channel is grey throughout: the red and blue belong to the cell '
             'states.'),
            ('b', 'CMC burden odds ratio',
             'exp(b) with exp(b +/- 1.96 SE) from rvtest --burden cmcWald, log axis, dashed line '
             'at OR = 1. One cohort per row, so every row carries its OWN full 95 % interval; '
             f'colour follows the cohort on the ordered ramp, darkest {clabel[cohorts[0]]} to '
             f'lightest {clabel[cohorts[-1]]}, which restates the row label rather than adding a '
             'variable. Two markers per row, one per variant set (filled '
             f'{V.stratum_label(narrow)}, hollow {V.stratum_label(strata[0])}). A cell whose '
             'Wald fit rvtest reset -- no minor allele in one group, which sends the interval '
             'to 1e-56 .. 1e66 -- has no usable estimate. It is marked n.e. in the margin, '
             'excluded from the axis limits, and reads "not estimable" in the table. n.e. is '
             'deliberately NOT the × of (a): × means the gene is not in that cohort\'s map and '
             'was never tested, whereas an n.e. cell WAS tested and keeps its state in (a). '
             'P2RY2 is the case here -- no control carries a MODERATE+HIGH allele in narrow or '
             'intermediate, which (c) shows as a carrier count of 0, yet the gene is BH-called '
             'in intermediate. Both facts are true, and the two marks now let the figure say '
             'both without contradicting itself.'),
            ('c', 'Carrier frequency',
             'The fraction of controls (blue) and of cases (red) carrying at least one '
             f'qualifying minor allele of the gene, for {V.stratum_label(narrow)} where the gene '
             'is mapped there and the wider set otherwise (hollow markers). Each row is scored '
             'against ITS OWN cohort\'s N, and the counts at the right are that cohort\'s. '
             'Carrier counts need not grow with the cohort: each cohort sets its own minAC, so a '
             'gene\'s mapped variant set is decided per cohort and a variant that passes in a '
             'smaller cohort can be dropped from a larger one. PES1 is the clearest case here -- '
             '4 MODERATE+HIGH variants in narrow and intermediate, 3 in full -- and the '
             'n_var_map column of the tables below carries the number for every row. A carrier '
             'holding several sites counts once here and once per site in the cumulative MAC of '
             'the summary table.'),
        ],
        interpretation=(
            f'{len(order)} gene(s) reach a tier. Tier 1 is the strongest reading this design '
            'offers; Tier 2 drops the both-statistics condition; Tier 3 asks only for two '
            'cohorts. Two caveats outrank the tiers. The cohorts are nested and differ almost '
            'entirely in controls, so agreement across them shows the signal survives adding '
            'controls — a sample-selection sensitivity analysis, not replication. And the '
            'variant sets are nested, so a gene tiered in both is one finding tested twice on '
            'overlapping sets. What the figure genuinely separates is the statistic: a gene '
            'called by SKAT-O and not by CMC (JAK2-like rows, with an OR near 1) carries a '
            'signal in a subset of its variants or in opposite directions, which a collapse '
            'cannot see. Cases and controls were sequenced on different platforms, so no gene '
            'here is established as a CTEPH gene.'),
        numbers=[('genes reaching a tier', len(order)), ('genes drawn', len(genes)),
                 ('Tier 1 genes', n_t[1]), ('Tier 2 genes', n_t[2]), ('Tier 3 genes', n_t[3]),
                 ('report cohort', rc), ('N case / N control', f'{n_case:,} / {n_ctrl:,}'),
                 ('cells not in a cohort\'s map', sum(1 for sn in genes for st in strata
                                                      for c in cohorts for m in V.METHODS
                                                      if (r := cell_row(sn, st, c, m)) is None
                                                      or r['state'] == 'not_in_map')),
                 ('cells with no OR estimate', f'{n_not_est} of {n_cells}'),
                 ('omitted for space', ', '.join(best[sn][1] for sn in dropped) or '—')],
        tables=[('Tiers (from tiers.tsv)', tier_tbl, 60),
                ('Odds ratios drawn in (b)', pd.DataFrame(drawn_b)),
                ('Carriers drawn in (c)', pd.DataFrame(drawn_c)),
                ('Denominator per cohort × variant set', den_tbl)],
        reading=[
            'Read a row across: the tier numbers summarise (a); (b) and (c) are what the '
            'summary table prints for that gene, drawn.',
            ('In (a), crosses are cells where the gene was never tested; they are not evidence '
             'of absence.' if unmapped else
             'No cell of (a) is a cross: every gene drawn here is in every cohort\'s map, which '
             'is expected -- a gene absent from a cohort\'s map cannot be called there and so '
             'rarely reaches a tier at all.'),
            'Compare the CMC and SKAT-O cells within one cohort block; that pair holds '
            'everything but the statistic fixed.',
            'A wide interval in (b) with few carriers in (c) is a statement about a handful '
            'of people.',
        ],
        limits=[
            'Not replication: the cohorts are nested and share nearly all cases.',
            'Not a CTEPH gene: the platform confound is in every cell.',
            'The OR is the collapsed-carrier burden; it says nothing about any single variant.',
        ],
        defs=[],
        model=('called = BH q < alpha over the cell\'s mapped genes; tiers from '
               'robust_genes.tier_of; OR = exp(beta) from rvtest --burden cmcWald'),
        methods_ref=METHODS_REF)
    print(f'[plot_robust_genes] {len(order)} tiered gene(s), {len(genes)} drawn, '
          f'plot block {plot_h:.2f} in -> {a.out_png}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:                       # noqa: BLE001
        print(f'Error in plot_robust_genes: {e}', file=sys.stderr)
        sys.exit(1)
