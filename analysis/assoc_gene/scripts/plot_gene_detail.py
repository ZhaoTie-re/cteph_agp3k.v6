#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : One figure per Tier 1 / Tier 2 gene: the variants that make up the
#           set and where the carriers are, over the effect size in every cell
#           of the design.
#
#           (a) answers "which variants carry this gene's signal", the one
#           thing a set-based p-value cannot say. Every variant of the widest
#           set the gene is mapped in is a stem: up for case carriers, down for
#           control carriers, marker by impact class. A gene whose signal sits in
#           one stem is a different finding from one whose carriers spread over
#           every stem, and the two have the same p-value.
#
#           (b) is the CMC odds ratio with its 95 % CI in every variant set x
#           cohort cell, with both p-values beside it. The OR is exp(beta) from
#           rvtest's own cmcWald on the same collapsed genotype.
#
#           THE GENE IS NAMED INSIDE THE CROP. The caption is below the crop
#           line; the header row names the gene, its span and its tiers above
#           everything else, so a cropped figure still says what it shows.
#
#           ONE DOCUMENT FOR THE FAMILY. The prose is identical for every gene;
#           the numbers are one table row each in README.md, written here after
#           the loop -- even when no gene is tiered, so the family directory
#           always says what it holds.
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
from matplotlib.ticker import (FixedFormatter, FixedLocator, FuncFormatter, MaxNLocator,
                               ScalarFormatter)  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_shared' / 'scripts'))
import plot_style as S                   # noqa: E402
import figure_doc                        # noqa: E402
import vocab as V                        # noqa: E402
# The exon/intron panel is the same one the regional figures draw, on the same
# Ensembl 86 tables: one drawing rule for gene structure in this project.
from plot_regional import draw_gene_track   # noqa: E402

DETAIL_TIERS = (1, 2)
PLOT_H = 7.0
STRIP_RATIO = 0.62         # gridspec height of ONE cohort strip in (a)
STRIP_HSPACE = 0.18        # between strips: they share the axis and stack tight
MAX_MODEL_ROWS = 3         # gene model rows; the target gene always keeps one
MODEL_ROW_RATIO = 0.16     # gridspec height per model row, plus the base below
MODEL_BASE_RATIO = 0.30
MODEL_HSPACE = 0.30        # (a) to gene model panel: a panel letter's worth
BLOCK_HSPACE = 0.28        # gene model panel to (c): axis label + column headers
METHODS_REF = '../../../docs/METHODS.md'
INDEX_COLUMNS = ('gene_symbol', 'set_name', 'chrom', 'pos_min', 'pos_max', 'tiers',
                 'stratum_drawn', 'n_variants_drawn', 'n_variants_per_cohort',
                 'n_genes_in_model_panel',
                 'n_multi_site_carriers_case', 'n_multi_site_carriers_control', 'png')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tiers', required=True)
    p.add_argument('--evidence', required=True)
    p.add_argument('--variant-detail', required=True)
    p.add_argument('--multi-detail', required=True)
    p.add_argument('--cohort-samples', required=True)
    p.add_argument('--report-cohort', required=True)
    p.add_argument('--exon-dir', default=None,
                   help='directory of <GENE>.exons.tsv from region_tracks.R; a missing or '
                        'empty file leaves the gene model panel with a stated blank rather '
                        'than failing the figure')
    p.add_argument('--out-dir', default='.')
    p.add_argument('--out-index', default='index.tsv')
    p.add_argument('--out-md', default='README.md')
    return p.parse_args()


def safe_name(s):
    return ''.join(ch if ch.isalnum() or ch in '-_.' else '_' for ch in s)


def fnum(v):
    x = pd.to_numeric(v, errors='coerce')
    return float(x) if np.isfinite(x) else None


def drop_overflowing_ticks(fig, ax, pad_in=0.02):
    """Blank an x tick label that would run past the figure's edge.

    The gene model panel's axis reaches the right margin, so whether its last
    tick label fits depends on the span drawn — a gene whose neighbours have
    short names leaves almost no pad. Measured after the layout is final; the
    tick itself stays, only its text goes.
    """
    fig.canvas.draw()
    w_px = fig.get_size_inches()[0] * fig.dpi
    pad_px = pad_in * fig.dpi
    labs = list(ax.get_xticklabels())
    texts = [t.get_text() for t in labs]
    changed = False
    for i, t in enumerate(labs):
        if not texts[i]:
            continue
        try:
            bb = t.get_window_extent(renderer=fig.canvas.get_renderer())
        except Exception:
            return
        if bb.x1 > w_px - pad_px or bb.x0 < pad_px:
            texts[i] = ''
            changed = True
    if changed:
        ax.xaxis.set_major_formatter(FixedFormatter(texts))


def read_exons(exon_dir, symbol, chrom):
    """The gene models to draw: the tested gene, plus genes overlapping its span.

    A window drawn around a 365 kb gene can hold a dozen named genes, and the
    figure is about ONE of them, so the panel keeps the tested gene and whatever
    overlaps it -- a read-through model, a nested gene -- and drops the
    neighbours. That also fixes the panel's height, which the gridspec must know
    before anything is drawn.

    Returns (frame or None, n_genes, note). `note` is the sentence to print in
    the panel when there is nothing to draw.
    """
    if not exon_dir:
        return None, 0, 'no gene model file was provided'
    path = Path(exon_dir) / f'{safe_name(symbol)}.exons.tsv'
    if not path.is_file() or path.stat().st_size == 0:
        return None, 0, f'no Ensembl 86 gene model for {symbol}'
    try:
        ex = pd.read_csv(path, sep='\t')
    except Exception:
        return None, 0, f'unreadable gene model file for {symbol}'
    need = {'gene', 'strand', 'gene_start', 'gene_end', 'exon_start', 'exon_end'}
    if not len(ex) or not need <= set(ex.columns):
        return None, 0, f'no Ensembl 86 gene model for {symbol}'
    for c in ('gene_start', 'gene_end', 'exon_start', 'exon_end'):
        ex[c] = pd.to_numeric(ex[c], errors='coerce')
    ex = ex.dropna(subset=['gene_start', 'gene_end', 'exon_start', 'exon_end'])
    tgt = ex[ex['gene'] == symbol]
    if not len(tgt):
        # The symbol is snpEff's; Ensembl 86 may not carry it (or carries it
        # under a synonym). Say so rather than drawing a neighbour as if it were
        # the tested gene.
        return None, 0, f'no Ensembl 86 gene model for {symbol}'
    t0, t1 = float(tgt['gene_start'].min()), float(tgt['gene_end'].max())
    spans = (ex.groupby('gene', sort=False)
               .agg(g0=('gene_start', 'min'), g1=('gene_end', 'max')).reset_index())
    over = spans[(spans['g0'] <= t1) & (spans['g1'] >= t0)]
    keep = [symbol] + [g for g in over['gene'] if g != symbol]
    if len(keep) > MAX_MODEL_ROWS:
        widest = (over[over['gene'] != symbol]
                  .assign(_w=lambda d: d['g1'] - d['g0'])
                  .nlargest(MAX_MODEL_ROWS - 1, '_w')['gene'].tolist())
        keep = [symbol] + widest
    return ex[ex['gene'].isin(keep)].copy(), len(keep), ''


def legend_height_pt(fig, leg):
    fig.canvas.draw()
    return leg.get_window_extent(fig.canvas.get_renderer()).height / fig.dpi * 72.0


def draw_gene(sn, symbol, chrom, tiers_of, ev, vd, md, strata, cohorts, clabel, n_of, rc,
              out_png, exon_dir=None):
    """One gene. Returns (stratum_drawn, variants df, cells table, n_multi_case, n_multi_ctrl)."""
    v_all = vd[(vd['set_name'] == sn) & (vd['cohort'] == rc)]
    rc_eff = rc
    if v_all.empty:
        alt = vd[vd['set_name'] == sn]
        if alt.empty:
            raise SystemExit(f'ABORT: {sn} is tiered but has no variant row anywhere')
        rc_eff = V.order_cohorts(set(alt['cohort']))[0]
        v_all = alt[alt['cohort'] == rc_eff]
    drawn_stratum = next(s for s in strata if (v_all['stratum'] == s).any())
    v = v_all[v_all['stratum'] == drawn_stratum].copy()
    for c in ('pos', 'n_carrier_case', 'n_carrier_control', 'mac_case', 'mac_control'):
        v[c] = pd.to_numeric(v[c], errors='coerce')
    v = v.sort_values('pos')
    # The SAME variants in the other two cohorts. The stratum is still chosen on
    # the report cohort -- it decides which set the figure is about, and letting
    # the union pick it would change that silently -- but the panel then draws
    # every cohort's carrier frequency for the variants that set contains. A
    # cohort's map is its own, so a variant can be missing from a narrower one;
    # that row is simply absent, which is not the same as a drawn 0 %.
    v_coh = {}
    for c in cohorts:
        w = vd[(vd['set_name'] == sn) & (vd['cohort'] == c)
               & (vd['stratum'] == drawn_stratum)].copy()
        for col in ('pos', 'n_carrier_case', 'n_carrier_control'):
            w[col] = pd.to_numeric(w[col], errors='coerce')
        v_coh[c] = w.set_index('variant_id')
    n_case, n_ctrl = n_of[rc_eff]
    m = md[(md['set_name'] == sn) & (md['cohort'] == rc_eff) & (md['stratum'] == drawn_stratum)]
    n_multi_case, n_multi_ctrl = int((m['group'] == 'case').sum()), int((m['group'] == 'control').sum())
    pos_min, pos_max = int(v['pos'].min()), int(v['pos'].max())

    exons, n_model_genes, model_note = read_exons(exon_dir, symbol, chrom)
    model_ratio = MODEL_BASE_RATIO + MODEL_ROW_RATIO * max(n_model_genes, 1)

    # Two blocks, not three rows: (a) and the gene model panel share one genomic
    # axis and belong together, so they sit in a nested grid with a gap wide
    # enough only for a panel letter, while the outer gap below them carries the
    # model panel's axis label and the forest's column headers. One flat 3-row
    # grid would give both gaps the same width and strand the short model panel
    # in white space.
    strips_total = STRIP_RATIO * len(cohorts)
    fig = plt.figure(figsize=(S.COL_DOUBLE, PLOT_H))
    gs = fig.add_gridspec(2, 2, height_ratios=[strips_total + model_ratio, 1.0],
                          width_ratios=[3.2, 1.3])
    gs_top = gs[0, :].subgridspec(2, 1, height_ratios=[strips_total, model_ratio],
                                  hspace=MODEL_HSPACE)
    # One strip per cohort. Three cohorts stacked on one strip as nested markers
    # put up to six coincident marks on a single variant and spent the size
    # channel telling them apart; a strip each says the same thing with position,
    # and a variant a narrower cohort does not carry reads directly as a gap.
    gs_str = gs_top[0].subgridspec(len(cohorts), 1, hspace=STRIP_HSPACE)
    ax_strip = {}
    for i, c in enumerate(cohorts):
        ax_strip[c] = fig.add_subplot(gs_str[i], sharex=ax_strip[cohorts[0]] if i else None)
    ax_a = ax_strip[cohorts[0]]
    ax_g = fig.add_subplot(gs_top[1], sharex=ax_a)     # the gene model panel
    ax_b = fig.add_subplot(gs[1, 0])
    ax_t = fig.add_subplot(gs[1, 1], sharey=ax_b)

    # ---- (a) carrier FREQUENCY per variant, in every cohort. x in Mb so the
    #      gene model panel below can share the axis to the pixel.
    # Counts put 438 cases and 2,630 controls on one axis, so a variant carried
    # by 0.9 % of controls out-stemmed one carried by 2.1 % of cases and the
    # panel read as control-dominated. The y axis is now the fraction of that
    # cohort's own group, which is the quantity (c) and 04.robust_genes also use.
    x = v['pos'].to_numpy(dtype=float) / 1e6
    pct = {}                       # cohort -> (case %, control %) aligned to v
    for c in cohorts:
        w = v_coh.get(c)
        n_ca, n_ko = n_of[c]
        cs, ks = [], []
        for vid in v['variant_id']:
            if w is None or vid not in w.index:
                cs.append(np.nan); ks.append(np.nan)      # not in this cohort's map
                continue
            r = w.loc[vid]
            cs.append(100.0 * float(r['n_carrier_case']) / n_ca)
            ks.append(-100.0 * float(r['n_carrier_control']) / n_ko)
        pct[c] = (np.array(cs), np.array(ks))
    # One strip per cohort, drawn against its OWN denominator. Shape is the
    # snpEff impact class in every strip now that nothing else needs the marker.
    for c in cohorts:
        ax = ax_strip[c]
        cs, ks = pct[c]
        ax.axhline(0, color=S.INK_SOFT, lw=0.8, zorder=1)
        ax.vlines(x, 0, np.nan_to_num(cs), color=V.CASE_COLOR, lw=1.0, alpha=0.55, zorder=2)
        ax.vlines(x, 0, np.nan_to_num(ks), color=V.CONTROL_COLOR, lw=1.0, alpha=0.55, zorder=2)
        for i, (_, r) in enumerate(v.iterrows()):
            if np.isnan(cs[i]):
                continue                      # not in this cohort's map: a gap
            mk = V.IMPACT_MARKER.get(r['impact'], 'o')
            ax.scatter([x[i]], [cs[i]], marker=mk, s=24, c=V.CASE_COLOR,
                       linewidths=0.5, edgecolors='white', zorder=4)
            ax.scatter([x[i]], [ks[i]], marker=mk, s=24, c=V.CONTROL_COLOR,
                       linewidths=0.5, edgecolors='white', zorder=4)
    # ONE y limit over all three strips, so the same variant is the same height
    # in every cohort and the strips can be read against each other. Top and
    # bottom still scale independently: case % runs several times control %, and
    # a symmetric axis would flatten the control half to nothing.
    all_up = np.concatenate([pct[c][0] for c in cohorts])
    all_dn = np.concatenate([pct[c][1] for c in cohorts])
    top = max(float(np.nanmax(all_up)) if np.isfinite(all_up).any() else 0.1, 0.1)
    bot = min(float(np.nanmin(all_dn)) if np.isfinite(all_dn).any() else -0.1, -0.1)
    for c in cohorts:
        ax_strip[c].set_ylim(bot * 1.3, top * 1.3)
    # The window shows the whole tested gene AND every variant: the variant span
    # only bounds the gene from the inside, and a panel cropped to it would cut
    # the exons the reader is here to see.
    lo_mb, hi_mb = float(x.min()), float(x.max())
    if exons is not None:
        tgt = exons[exons['gene'] == symbol]
        if len(tgt):
            lo_mb = min(lo_mb, float(tgt['gene_start'].min()) / 1e6)
            hi_mb = max(hi_mb, float(tgt['gene_end'].max()) / 1e6)
    span = max(hi_mb - lo_mb, 1e-4)
    # draw_gene_track writes a gene's name BESIDE its right end unless the gene
    # reaches the last 2 % of the window, so the right pad has to hold the widest
    # name it will write; its own row packer budgets 0.012 of the span per
    # character, and that is the number used here so the two cannot disagree.
    widest = max((len(str(g)) for g in exons['gene'].unique()), default=6) if exons is not None else 6
    xlim = (lo_mb - 0.04 * span, hi_mb + max(0.04, widest * 0.014) * span)
    ax_a.set_xlim(*xlim)
    dec = min(max(2, int(np.ceil(-np.log10(span))) + 1), 5)
    mid = cohorts[len(cohorts) // 2]
    for c in cohorts:
        ax = ax_strip[c]
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda val, _p: f'{val:.{dec}f}'))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda val, _p: f'{abs(val):g}'))
        ax.tick_params(axis='x', labelbottom=False)     # the gene panel carries them
        # The cohort names its own strip, inside it, where it cannot be confused
        # with a panel letter; only the middle strip carries the y label.
        ax.annotate(clabel[c], xy=(0.995, 0.94), xycoords='axes fraction',
                    ha='right', va='top', fontsize=V.MIN_FONT, fontweight='bold',
                    color=S.INK_SOFT, zorder=6)
        ax.set_ylabel('carriers, % of group   (case ↑ · control ↓)'
                      if c == mid else '')
        S.despine(ax, grid_axis='y')

    # ---- (b) OR per variant set x cohort, with the p-values in (t).
    # A table-style forest: each variant set is a bold header row over its
    # cohort rows, so the left margin carries one short label per row and no
    # group name beside the block (which would cost a third of the canvas).
    rows, labels, heads = [], [], []
    for s in strata:
        rows.append((s, None)); labels.append(V.stratum_label(s)); heads.append(True)
        for c in cohorts:
            rows.append((s, c)); labels.append(clabel[c]); heads.append(False)
    n_rows = len(rows)
    ys = list(range(n_rows - 1, -1, -1))
    ax_b.set_yticks(ys)
    ax_b.set_yticklabels(labels)
    for lab, h in zip(ax_b.get_yticklabels(), heads):
        lab.set_fontweight('bold' if h else 'normal')
    ax_b.set_ylim(-0.7, n_rows - 0.3)
    for bi, s in enumerate(strata):
        if bi % 2 == 1:
            band = [y for (ss, _c), y in zip(rows, ys) if ss == s]
            ax_b.axhspan(min(band) - 0.5, max(band) + 0.5, color=S.NEUTRAL, alpha=0.13,
                         lw=0, zorder=0)
    e = ev[ev['set_name'] == sn]
    lo_all, hi_all, cells_tbl = [], [], []
    for (s, c), yy in zip(rows, ys):
        if c is None:
            continue
        cmc = e[(e['stratum'] == s) & (e['cohort'] == c) & (e['method'] == 'cmc')]
        sk = e[(e['stratum'] == s) & (e['cohort'] == c) & (e['method'] == 'skato')]
        st_c = cmc['state'].iloc[0] if len(cmc) else 'not_in_map'
        st_s = sk['state'].iloc[0] if len(sk) else 'not_in_map'
        o = fnum(cmc['or'].iloc[0]) if len(cmc) else None
        lo = fnum(cmc['or_l95'].iloc[0]) if len(cmc) else None
        hi = fnum(cmc['or_u95'].iloc[0]) if len(cmc) else None
        pc = fnum(cmc['pvalue'].iloc[0]) if len(cmc) else None
        ps = fnum(sk['pvalue'].iloc[0]) if len(sk) else None
        col = V.STATE_COLOR[st_c]
        if st_c == 'not_in_map':
            V.draw_state(ax_b, 1.0, yy, 'not_in_map')
        elif o is not None and lo is not None and hi is not None and lo > 0:
            ax_b.plot([lo, hi], [yy, yy], color=col, lw=1.1, zorder=2)
            ax_b.scatter([o], [yy], s=26, color=col, zorder=3, edgecolors='white', linewidths=0.4)
            lo_all += [lo]; hi_all += [hi]
        else:
            ax_b.text(1.0, yy, 'OR not estimable', ha='center', va='center', fontsize=V.MIN_FONT,
                      color=S.INK_SOFT, zorder=3)
        # the two p-values, coloured by their own statistic's state
        for xf, p, st in ((0.02, pc, st_c), (0.55, ps, st_s)):
            txt = '—' if st == 'not_in_map' or p is None else S.p_tex(p)
            ax_t.text(xf, yy, txt, transform=ax_t.get_yaxis_transform(), ha='left', va='center',
                      fontsize=7, color=V.STATE_COLOR[st] if st != 'not_in_map' else S.REFERENCE)
        cells_tbl.append({'gene_symbol': symbol, 'stratum': s, 'cohort': c, 'cmc_state': st_c,
                          'skato_state': st_s, 'or': o, 'or_l95': lo, 'or_u95': hi,
                          'cmc_p': pc, 'skato_p': ps})
    ax_b.axvline(1.0, color=S.REFERENCE, lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax_b.set_xscale('log')
    # OR = 1 is always inside the range, so the reference line is always drawn.
    ax_b.set_xlim(min(lo_all + [1.0]) / 1.3, max(hi_all + [1.0]) * 1.3)
    S.or_log_axis(ax_b)
    tk = list(ax_b.get_xticks())
    lo_x, hi_x = ax_b.get_xlim()
    if not any(abs(t - 1.0) < 1e-9 for t in tk) and lo_x <= 1.0 <= hi_x:
        ax_b.xaxis.set_major_locator(FixedLocator(sorted(set(tk) | {1.0})))
        ax_b.xaxis.set_major_formatter(ScalarFormatter())
    ax_b.set_xlabel(f'CMC {S.OR_SYM} (95% CI)')
    S.despine(ax_b, grid_axis='x')
    ax_b.tick_params(axis='y', length=0)
    ax_t.set_axis_off()
    y_head = n_rows - 0.3 + 0.05
    for xf, head in ((0.02, 'CMC P'), (0.55, 'SKAT-O P')):
        ax_t.text(xf, y_head, head, transform=ax_t.get_yaxis_transform(), ha='left', va='bottom',
                  fontsize=7, fontweight='bold', color=S.INK, clip_on=False)

    # ---- caption with the header row, then the measured layout.
    # The header row is one line and the canvas clips it, so when both variant
    # sets award the same tier the stratum names are dropped rather than the
    # sentence being cut mid-word.
    tiers_seen = set(tiers_of.values())
    if len(tiers_seen) == 1 and len(tiers_of) > 1:
        tier_txt = f'Tier {tiers_seen.pop()} in both variant sets'
    else:
        tier_txt = ' · '.join(f'Tier {t} {V.stratum_label(s)}' for s, t in tiers_of.items())
    ors = [r['or'] for r in cells_tbl if r['or'] is not None]
    drawn_cmc = e[(e['stratum'] == drawn_stratum) & (e['cohort'] == rc_eff) & (e['method'] == 'cmc')]
    nc_d = int(float(drawn_cmc['n_carrier_case'].iloc[0])) if len(drawn_cmc) else int(up.sum())
    nk_d = int(float(drawn_cmc['n_carrier_control'].iloc[0])) if len(drawn_cmc) else int(-dn.sum())
    or_txt = (f'CMC OR {min(ors):.1f}–{max(ors):.1f} across cells' if ors
              else 'no estimable CMC OR')
    note_edge = ''
    if drawn_stratum != strata[-1] and not any((v_all['stratum'] == strata[-1])):
        note_edge = (f' The gene is not in the {V.stratum_label(strata[-1])} map of '
                     f'{clabel[rc_eff]} (fewer than MinNumVar HIGH/MODERATE variants).')
    multi_txt = (f' {n_multi_case} case / {n_multi_ctrl} control sample(s) carry ≥2 sites.'
                 if n_multi_case or n_multi_ctrl else '')
    S.caption_block(
        fig, plot_h=PLOT_H,
        header=(symbol, tier_txt),
        title=(f'{symbol}: {len(v)} variants, {nc_d} case / {nk_d} control carriers in '
               f'{clabel[rc_eff]}; {or_txt}.'),
        panels=[f'{len(v)} variant(s) of the {V.stratum_label(drawn_stratum)} set in '
                f'{clabel[rc_eff]}\'s map, one strip per cohort on one shared scale: carriers '
                f'as a percentage of that cohort\'s own group, cases up and controls down. '
                f'A variant a narrower cohort does not map is absent from its strip.',
                (f'Ensembl 86 gene models, {n_model_genes} gene(s) on the same axis.'
                 if n_model_genes else 'Gene models unavailable.'),
                f'CMC {S.OR_SYM} (95% CI) per variant set and cohort; marker and P colour '
                f'give that cell\'s state.'],
        notes=(f'Variants span chr{chrom}:{pos_min:,}–{pos_max:,}; (b) is drawn over that span and the gene\'s own extent. '
               f'Marker = snpEff impact; MODERATE+HIGH variants are the HIGH and MODERATE '
               f'markers. The two P columns are coloured by that statistic\'s own state in '
               f'the cell: dark red Bonferroni + BH, light red BH only, blue tested and not '
               f'called, — not in that cohort\'s map.{note_edge}{multi_txt} README.md'),
        letters='abcdefghijklmnop',
        top_pad=0.72, hspace=BLOCK_HSPACE, wspace=0.04, left='auto', right=0.985,
        margin_axes=[*ax_strip.values(), ax_b])
    handles_a = [Line2D([], [], ls='none', marker=V.IMPACT_MARKER[i], color=S.INK_SOFT,
                        markersize=5, label=i) for i in V.IMPACT_ORDER if (v['impact'] == i).any()]
    # The cohort is the strip now, so it needs no key, which leaves room for the
    # two group colours the y label alone had to carry.
    handles_a += [Line2D([], [], color=V.CASE_COLOR, lw=1.6, label='case'),
                  Line2D([], [], color=V.CONTROL_COLOR, lw=1.6, label='control')]
    leg_a = S.place_legend(ax_a, handles_a, y=1.02)
    pad_a = 7.0 + (legend_height_pt(fig, leg_a) + 3.0 if leg_a.get_bbox_to_anchor() else 0.0)
    S.panel_tag(ax_a, 'a', pad=pad_a)
    # The gene model panel measures its own labels (spread_labels), so it is
    # drawn after caption_block like every other measured pass.
    if exons is not None:
        draw_gene_track(ax_g, exons, xlim, chrom)
    else:
        draw_gene_track(ax_g, None, xlim, chrom)
        ax_g.texts[-1].set_text(model_note)
    S.panel_tag(ax_g, 'b')
    # (c) carries no legend strip: the three states colour BOTH its markers and
    # its two P columns, the caption states that rule once, and a strip above (c)
    # would sit on the gene panel's axis label.
    S.panel_tag(ax_b, 'c')
    # Pin the tick sets so thin_tick_labels' FixedFormatter cannot drift from them.
    for ax in (ax_g, ax_b):
        ax.xaxis.set_major_locator(FixedLocator(list(ax.get_xticks())))
    S.thin_tick_labels(ax_g, 'x')
    S.thin_tick_labels(ax_b, 'x')
    drop_overflowing_ticks(fig, ax_g)
    fig.savefig(out_png)
    plt.close(fig)

    v_tbl = v[['variant_id', 'impact', 'effect', 'n_carrier_case', 'n_carrier_control',
               'mac_case', 'mac_control', 'n_hom_case', 'n_hom_control']].assign(gene_symbol=symbol)
    n_var_coh = {c: int(v['variant_id'].isin(v_coh[c].index).sum()) for c in cohorts}
    return (drawn_stratum, rc_eff, pos_min, pos_max, v_tbl, pd.DataFrame(cells_tbl),
            n_multi_case, n_multi_ctrl, n_model_genes, n_var_coh)


def main():
    a = parse_args()
    S.setup_style('paper')
    out_dir = Path(a.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tiers = pd.read_csv(a.tiers, sep='\t', dtype=str)
    ev = pd.read_csv(a.evidence, sep='\t', dtype=str)
    vd = pd.read_csv(a.variant_detail, sep='\t', dtype=str)
    md = pd.read_csv(a.multi_detail, sep='\t', dtype=str)
    smp = pd.read_csv(a.cohort_samples, sep='\t', dtype=str)
    n_of = {r['cohort']: (int(r['n_case']), int(r['n_control'])) for _, r in smp.iterrows()}
    if a.report_cohort not in n_of:
        raise SystemExit(f'ABORT: --report-cohort {a.report_cohort!r} is not in {a.cohort_samples}')
    tiers['tier'] = pd.to_numeric(tiers['tier'], errors='coerce').astype(int)
    cohorts = V.order_cohorts(set(ev['cohort']))
    clabel = V.cohort_labels(cohorts)
    strata = V.order_strata(ev)

    pick = tiers[tiers['tier'].isin(DETAIL_TIERS)]
    genes = list(dict.fromkeys(pick['set_name']))
    index, v_tabs, c_tabs = [], [], []
    for sn in genes:
        sub = tiers[tiers['set_name'] == sn]
        symbol, chrom = sub['gene_symbol'].iloc[0], sub['chrom'].iloc[0]
        tiers_of = {s: int(t) for s, t in zip(sub['stratum'], sub['tier'])}
        tiers_of = {s: tiers_of[s] for s in strata if s in tiers_of}
        png = out_dir / f'{safe_name(symbol)}.png'
        drawn, rc_eff, pmin, pmax, v_tbl, c_tbl, n_mc, n_mk, n_gm, n_vc = draw_gene(
            sn, symbol, chrom, tiers_of, ev, vd, md, strata, cohorts, clabel, n_of,
            a.report_cohort, str(png), exon_dir=a.exon_dir)
        v_tabs.append(v_tbl); c_tabs.append(c_tbl)
        index.append({'gene_symbol': symbol, 'set_name': sn, 'chrom': chrom,
                      'pos_min': pmin, 'pos_max': pmax,
                      'tiers': ';'.join(f'{s}:{t}' for s, t in tiers_of.items()),
                      'stratum_drawn': drawn, 'n_variants_drawn': len(v_tbl),
                      # How many of the drawn variants each cohort's own map has;
                      # a narrower cohort can be missing some, which is why a
                      # marker is absent rather than drawn at 0 %.
                      'n_variants_per_cohort': ';'.join(f'{clabel[c]}:{n}'
                                                        for c, n in n_vc.items()),
                      'n_genes_in_model_panel': n_gm,
                      'n_multi_site_carriers_case': n_mc, 'n_multi_site_carriers_control': n_mk,
                      'png': png.name})
        print(f'[plot_gene_detail] {symbol}: {len(v_tbl)} variant(s) ({drawn}, {rc_eff}), '
              f'{n_mc}/{n_mk} multi-site carrier(s) -> {png}')

    with open(a.out_index, 'w') as fh:
        fh.write('\t'.join(INDEX_COLUMNS) + '\n')
        for r in index:
            fh.write('\t'.join(str(r[c]) for c in INDEX_COLUMNS) + '\n')
    index_df = pd.DataFrame(index, columns=INDEX_COLUMNS)
    v_df = pd.concat(v_tabs) if v_tabs else pd.DataFrame()
    c_df = pd.concat(c_tabs) if c_tabs else pd.DataFrame()
    n_case, n_ctrl = n_of[a.report_cohort]
    figure_doc.write_doc(
        None, out_path=a.out_md,
        subject=(f'**Figures:** `<GENE>.png` — {len(index)} gene figure(s), one per Tier 1 or '
                 f'Tier 2 gene ({", ".join(index_df["gene_symbol"]) or "none"}). One document '
                 'covers the family; each gene is one row of the tables below.'),
        title='Per-gene detail — Tier 1 and Tier 2 genes',
        question=('Which variants make up each robust gene\'s set, where are its carriers, and '
                  'does the burden effect hold in every cell of the design?'),
        panels=[
            ('a', 'Variants and carrier frequency, one strip per cohort',
             'Every variant of the widest set the gene is mapped in, in '
             f'{clabel[a.report_cohort]} (the report cohort, which decides WHICH set the figure '
             'is about), at its genomic position in Mb. The three strips are the three cohorts, '
             'drawn on ONE shared y scale so the same variant is the same height in each. The y '
             "axis is the percentage of that cohort's own group carrying at least one minor "
             'allele of the variant -- cases up, controls down. It is a percentage and not a '
             'count because the groups are of very different size: on a count axis a variant '
             'carried by 0.9 % of 2,630 controls out-stems one carried by 2.1 % of 438 cases, '
             'and the panel reads as control-dominated when the case excess is the finding. '
             'Marker shape is the snpEff impact class (diamond HIGH, circle MODERATE, square '
             "LOW). A variant absent from a narrower cohort's map draws nothing in that strip, "
             'which is not the same as a drawn 0 %: each cohort sets its own minAC, so its map '
             'is its own, and reading down the strips shows exactly which variants the larger '
             'sample adds. The x axis is shared with (b), so a marker sits exactly above the '
             'exon it falls in.'),
            ('b', 'Gene models',
             'The exon/intron structure of the tested gene and of any gene overlapping its '
             'span, from Ensembl 86 (GRCh38) — one representative transcript per gene, the '
             'one with the greatest total exonic length. Filled boxes are exons, the '
             'connecting line spans the introns and the chevrons give the transcribed strand. '
             'Genes that merely lie in the window without overlapping the tested gene are not '
             'drawn: the figure is about one gene, and its neighbours would cost rows without '
             'adding to it. Where Ensembl 86 has no model under the tested symbol the panel '
             'says so rather than drawing a neighbour in its place.'),
            ('c', 'CMC burden odds ratio in every cell, with both P-values',
             'exp(β) and exp(β ± 1.96 SE) from rvtest --burden cmcWald per variant set and '
             'cohort, log axis, dashed line at OR = 1. Marker colour is the CMC state of that '
             'cell (full red Bonferroni + BH, light red BH only, blue tested and not called); a '
             'grey cross is a cell in which the gene is not in the map; "OR not estimable" is a '
             'Wald fit rvtest reset. To the right, the CMC score P and the SKAT-O P of the same '
             'cell, EACH COLOURED BY ITS OWN STATISTIC\'S STATE IN THAT CELL: dark red = above '
             'that cell\'s Bonferroni threshold (and therefore BH-called), light red = '
             'BH-called only, blue = tested and called by neither, — = the gene is not in that '
             'cohort\'s map. So a row can carry a red SKAT-O P beside a blue CMC P: the two '
             'statistics disagree, and the marker shows what the burden estimate was.'),
        ],
        interpretation=(
            'Read (a) for how the signal is distributed: one tall stem is a single-variant '
            'finding wearing a gene-set P; many short stems is a burden. Read (b) for whether '
            'the effect is stable as controls are added (down a block) and as the set widens '
            '(between blocks). The cohorts are nested, so stability across them is not '
            'replication. Where a P is coloured but the OR sits near 1, the SKAT-O kernel is '
            'seeing structure the collapse cannot (a subset of variants, or opposite '
            'directions). Where samples carry two or more sites of the gene the caption says '
            'so; CMC counts such a sample once while the cumulative MAC counts every site.'),
        numbers=[('genes drawn', len(index)), ('tiers drawn', ', '.join(str(t) for t in DETAIL_TIERS)),
                 ('report cohort', a.report_cohort),
                 ('N case / N control', f'{n_case:,} / {n_ctrl:,}')],
        tables=[('Per-gene values', index_df),
                ('Cells drawn in (b) (from evidence.tsv)', c_df, 120),
                ('Variants drawn in (a) (from variant_detail.tsv)', v_df, 120)],
        reading=['(a) first: is this one variant or many?',
                 '(b) against (a): which exon does each stem fall in, and is the signal in '
                 'coding sequence at all?',
                 '(c) down a block: does adding controls move the OR?',
                 '(c) across blocks: does widening the set dilute it?',
                 'The P columns say which statistic called the cell, in that statistic\'s own '
                 'state colour; the marker says what the burden estimate is.'],
        limits=['Not replication: the cohorts are nested.',
                'Not a CTEPH gene: the platform confound is in every cell.',
                'The OR is the collapsed-carrier burden; it says nothing about any single '
                'variant.'],
        defs=[],
        model=('rvtest --burden cmc,cmcWald and --kernel skato, logistic, sex + ancestry PCs; '
               'gene models from Ensembl 86 (GRCh38) via the shared region_tracks.R'),
        methods_ref=METHODS_REF)
    print(f'[plot_gene_detail] {len(index)} gene figure(s) -> {a.out_index}, {a.out_md}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:                       # noqa: BLE001
        print(f'Error in plot_gene_detail: {e}', file=sys.stderr)
        sys.exit(1)
