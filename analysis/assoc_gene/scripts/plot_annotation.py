#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : What the annotation step did to this cohort's call set, on one
#           page: how every minAC-passing variant was classified, which
#           consequences carry the impact strata, and how many gene sets and
#           variants survive into the test.
#
#           THIS IS THE DENOMINATOR FIGURE, NOT A RESULT FIGURE. Every gene-based
#           p-value in this component is a statement about a set that was defined
#           here -- by snpEff's most-severe consequence, by the eligible-impact
#           list, and by the >= min_num_var rule. A reader who never sees this
#           page cannot tell whether a gene "was not significant" or was never
#           tested at all, and the Bonferroni threshold quoted everywhere else is
#           alpha / n_genes_mapped, a number set by panel (c) alone.
#
#           WHY MODIFIER AND UNANNOTATED ARE IN PANEL (a). They are excluded from
#           every test, and that is exactly why they must be drawn: the honest
#           statement is "x % of the call set was eligible for a gene test", and
#           a donut over the eligible impacts alone silently renormalises to
#           100 % and asserts the opposite.
#
#           ONE CONSEQUENCE PANEL, NOT THREE. The three impact classes are one
#           ordered variable, so their consequences are three blocks of one
#           panel, each block its own colour, with the impact written once beside
#           it. Bars are linear (length must stay proportional to count) and
#           every bar carries its count and its percent within the impact; the
#           labels, not the lengths, are how the rare consequences are read.
#
#           THIS FIGURE WRITES NO SIDECAR. It is one member of a per-cohort
#           family; it drops its numbers as <png>.stats.json and the CATALOGUE
#           step writes the family README with one row per cohort.
# Component: assoc_gene
# ---------------------------------------------------------------------------
import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402
from matplotlib.lines import Line2D      # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.ticker import FixedLocator, MaxNLocator   # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_shared' / 'scripts'))
import plot_style as S                   # noqa: E402
import figure_doc                        # noqa: E402
import vocab as V                        # noqa: E402

# The snpEff convention. Deliberately NOT the semantic palette: impact is a
# fixed five-level scale a reader knows from every snpEff summary, and it is
# used in this figure only. UNANNOTATED takes the neutral because it is an
# absence of information, not a fifth severity.
IMPACT_COLORS = {'HIGH': '#d62728', 'MODERATE': '#ff7f0e', 'LOW': '#2ca02c',
                 'MODIFIER': '#1f77b4', 'UNANNOTATED': S.NEUTRAL_D}
IMPACT_ORDER = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER', 'UNANNOTATED']
BAR_IMPACTS = ['HIGH', 'MODERATE', 'LOW']
ELIGIBLE = ('HIGH', 'MODERATE', 'LOW')

MAX_EFFECTS = 6           # per impact; the remainder folds into one 'Other' row
MAX_LABEL_CHARS = 30      # a consequence label wider than this is truncated
PLOT_H = 6.6
# The eligible fraction is one wedge on the ring; it is not an impact class, so
# it takes the accent layer rather than an impact colour.
ELIGIBLE_COLOR = S.ACCENT

ANNO_COLUMNS = ('Impact', 'Effect', 'Biotype', 'Count')
GENE_INDEX_COLUMNS = ('set_name', 'gene_symbol', 'chrom', 'n_var_map', 'pos_min', 'pos_max',
                      'span_bp', 'n_high', 'n_moderate', 'n_low')
DENOM_COLUMNS = ('cohort', 'stratum', 'min_num_var', 'n_genes_mapped',
                 'n_genes_dropped_lt_minvar', 'n_variants', 'n_variants_multigene',
                 'n_sets_chrom_split', 'n_positions_colliding', 'alpha',
                 'bonferroni_threshold')
DENOM_INTS = ('min_num_var', 'n_genes_mapped', 'n_genes_dropped_lt_minvar', 'n_variants',
              'n_variants_multigene', 'n_sets_chrom_split', 'n_positions_colliding')
MAP_QC_REQUIRED = ('chrom', 'n_variants_in', 'n_annotated', 'n_unannotated',
                   'n_modifier_only', 'n_map_rows', 'n_multigene_variants')
MAP_QC_INTS = ('n_variants_in', 'n_annotated', 'n_unannotated', 'n_modifier_only',
               'n_map_rows', 'n_multigene_variants')


def parse_args():
    p = argparse.ArgumentParser(
        description='One figure per cohort: how the call set was annotated, and how many '
                    'gene sets and variants that leaves for the gene-based test.')
    p.add_argument('--anno-counts', required=True)
    p.add_argument('--gene-index', required=True, nargs='+',
                   help='gene_index.<stratum>.tsv; the stratum is read from the FILE NAME')
    p.add_argument('--denominator', required=True, nargs='+',
                   help='denominator.<stratum>.tsv; one per stratum, matching --gene-index')
    p.add_argument('--map-qc', required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--out-png', required=True)
    p.add_argument('--out-stats', default=None,
                   help='<png>.stats.json for the family catalogue (default: beside the PNG)')
    p.add_argument('--out-dir', required=True,
                   help='annotation_summary.tsv and stratum_summary.tsv are written here')
    return p.parse_args()


# ── Reading, with every structural failure raised BEFORE anything is written ──
def read_anno_counts(path):
    d = pd.read_csv(path, sep='\t', dtype={'Impact': str, 'Effect': str, 'Biotype': str})
    if tuple(d.columns) != ANNO_COLUMNS:
        raise SystemExit(f'ABORT: {path} header is {list(d.columns)}, expected '
                         f'{list(ANNO_COLUMNS)}; it was not written by merge_gene_map.py')
    if not len(d):
        raise SystemExit(f'ABORT: {path} has no data row; there is nothing to draw')
    n = pd.to_numeric(d['Count'], errors='coerce')
    if n.isna().any():
        raise SystemExit(f'ABORT: {path} has a non-integer Count')
    d['Count'] = n.astype('int64')
    unknown = sorted(set(d['Impact'].astype(str)) - set(IMPACT_ORDER))
    if unknown:
        raise SystemExit(f'ABORT: {path} carries impact(s) {unknown}; expected only '
                         f'{IMPACT_ORDER}')
    return d


def read_map_qc(path):
    d = pd.read_csv(path, sep='\t', dtype=str)
    missing = [c for c in MAP_QC_REQUIRED if c not in d.columns]
    if missing:
        raise SystemExit(f'ABORT: {path} has no {missing} column(s)')
    chrom = d['chrom'].astype(str)
    tot = d[chrom == 'TOTAL']
    if len(tot) != 1:
        raise SystemExit(f'ABORT: {path} has {len(tot)} TOTAL row(s); exactly one is required')
    row = {k: v for k, v in tot.iloc[0].items()}
    for c in MAP_QC_INTS:
        row[c] = int(float(row[c]))
    return row, int((chrom != 'TOTAL').sum())


def stratum_of(path, kind):
    m = re.match(rf'^{kind}\.(?P<s>.+)\.tsv(\.gz)?$', Path(path).name)
    return m.group('s') if m else None


def read_denominators(paths):
    out = {}
    for path in paths:
        d = pd.read_csv(path, sep='\t', dtype=str)
        if tuple(d.columns) != DENOM_COLUMNS:
            raise SystemExit(f'ABORT: {path} header is {list(d.columns)}, expected '
                             f'{list(DENOM_COLUMNS)}')
        if len(d) != 1:
            raise SystemExit(f'ABORT: {path} has {len(d)} data rows; expected exactly one')
        row = {k: v for k, v in d.iloc[0].items()}
        for c in DENOM_INTS:
            row[c] = int(float(row[c]))
        name = str(row['stratum'])
        from_file = stratum_of(path, 'denominator')
        if from_file is not None and from_file != name:
            raise SystemExit(f'ABORT: {path} is named for stratum {from_file!r} but its '
                             f'stratum column says {name!r}')
        if name in out:
            raise SystemExit(f'ABORT: stratum {name!r} given twice with --denominator')
        out[name] = row
    return out


def read_gene_indexes(paths, dens):
    out = {}
    for path in paths:
        name = stratum_of(path, 'gene_index')
        if name is None:
            raise SystemExit(f'ABORT: {Path(path).name} is not named gene_index.<stratum>.tsv')
        if name in out:
            raise SystemExit(f'ABORT: stratum {name!r} given twice with --gene-index')
        d = pd.read_csv(path, sep='\t', dtype={'set_name': str, 'gene_symbol': str, 'chrom': str})
        if tuple(d.columns) != GENE_INDEX_COLUMNS:
            raise SystemExit(f'ABORT: {path} header is {list(d.columns)}, expected '
                             f'{list(GENE_INDEX_COLUMNS)}')
        for c in ('n_var_map', 'pos_min', 'pos_max', 'span_bp', 'n_high', 'n_moderate', 'n_low'):
            d[c] = pd.to_numeric(d[c], errors='coerce')
        if name not in dens:
            raise SystemExit(f'ABORT: gene_index for stratum {name!r} has no --denominator')
        if len(d) != dens[name]['n_genes_mapped']:
            raise SystemExit(f'ABORT: {path} has {len(d)} set(s) but its denominator says '
                             f'n_genes_mapped = {dens[name]["n_genes_mapped"]}')
        out[name] = d
    unpaired = sorted(set(dens) - set(out))
    if unpaired:
        raise SystemExit(f'ABORT: stratum/strata {unpaired} have a --denominator but no '
                         f'--gene-index')
    return out


# ── Measured helpers ─────────────────────────────────────────────────────────
def text_sizes_in(fig, texts, fontsize, weight='normal'):
    """[(width, height), ...] in inches for several strings, on ONE canvas draw."""
    try:
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        out = []
        for s in texts:
            t = fig.text(0, 0, s, fontsize=fontsize, fontweight=weight)
            bb = t.get_window_extent(renderer=r)
            t.remove()
            out.append((bb.width / fig.dpi, bb.height / fig.dpi))
        return out
    except Exception:
        return [(0.62 * max(len(ln) for ln in s.split('\n')) * fontsize / 72.0,
                 (1.0 + 1.2 * s.count('\n')) * 1.25 * fontsize / 72.0) for s in texts]


def box_inches(ax):
    fig = ax.figure
    try:
        fig.canvas.draw()
        bb = ax.get_window_extent(renderer=fig.canvas.get_renderer())
        return bb.width / fig.dpi, bb.height / fig.dpi
    except Exception:
        pos = ax.get_position()
        w, h = fig.get_size_inches()
        return pos.width * w, pos.height * h


def fit_label_room(ax, xs, widths_in, pad=0.02):
    """Grow the x limit until EVERY value label fits inside the panel, in one pass.

    A label is a fixed number of pixels wide, so the axis must be long enough
    that x sits at 1 - f of it, f being the label's fraction of the panel; the
    binding label is not necessarily the largest value, so every one is checked.
    """
    box_w, _h = box_inches(ax)
    lo, hi = ax.get_xlim()
    need = hi
    for x, w in zip(xs, widths_in):
        f = min(0.60, w / max(box_w, 0.1) + pad)
        need = max(need, lo + (x - lo) / (1.0 - f))
    if np.isfinite(need) and need > hi:
        ax.set_xlim(lo, need)


def widen_gap(fig, ax_left, ax_right, wspace, pad_in=0.12):
    """Grow the column gap until the right panel's y tick labels clear the left one.

    matplotlib's `wspace` is a fraction of the MEAN axes width, and that mean is
    not known until the layout is final: estimating it before caption_block put
    the consequence panel's ~1.2 in tick labels partly on top of the impact
    table. Measure the achieved gap against the rendered labels and correct once.
    """
    try:
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        labs = [t for t in ax_right.get_yticklabels() if t.get_text()]
        if not labs:
            return wspace
        need = max(t.get_window_extent(renderer=r).width for t in labs) / fig.dpi + pad_in
        a = ax_left.get_window_extent(renderer=r)
        b = ax_right.get_window_extent(renderer=r)
        have = (b.x0 - a.x1) / fig.dpi
        if have >= need:
            return wspace
        mean_w = (a.width + b.width) / 2.0 / fig.dpi
        new_ws = wspace + (need - have) / max(mean_w, 0.2)
        fig.subplots_adjust(wspace=new_ws)
        return new_ws
    except Exception:
        return wspace


def empty_panel(ax, msg):
    ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha='center', va='center',
            color=S.INK_SOFT, fontsize=plt.rcParams['legend.fontsize'])
    ax.set_xticks([]); ax.set_yticks([])
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


def effect_label(s):
    """A compound snpEff term as its first component with the extra count.

    The leading term is the most severe one; `splice_acceptor_variant&intron_variant`
    reads as `splice_acceptor_variant (+1)`. The full string is in
    annotation_summary.tsv.
    """
    s = str(s)
    if s.startswith('Other ('):
        return s
    parts = s.split('&')
    head = parts[0]
    if len(head) > MAX_LABEL_CHARS:
        head = head[:MAX_LABEL_CHARS - 1] + '…'
    return head if len(parts) == 1 else f'{head} (+{len(parts) - 1})'


def distinguished_label(s):
    """The same term when its short form collides with another row's.

    Two compound consequences can share their leading term
    (`splice_region_variant&intron_variant` and `…&synonymous_variant`), and
    then two rows of one impact block read alike. Here the second component —
    the first that differs — is spelled out, and the character budget is spent
    on it: the shared head is abbreviated first, the tail only if it still does
    not fit.
    """
    parts = str(s).split('&')
    if len(parts) == 1:
        return effect_label(s)
    head, tail = parts[0], parts[1]
    extra = f' (+{len(parts) - 2})' if len(parts) > 2 else ''
    budget = MAX_LABEL_CHARS - len(extra) - 3          # ' & '
    if len(tail) > budget - 4:
        tail = tail[:max(budget - 5, 3)] + '…'
    room = budget - len(tail)
    if len(head) > room:
        head = head[:max(room - 1, 3)] + '…'
    return f'{head} & {tail}{extra}'


def fold_effects(sub, max_effects=MAX_EFFECTS):
    g = (sub.groupby('Effect', as_index=False)['Count'].sum()
            .sort_values(['Count', 'Effect'], ascending=[False, True])
            .reset_index(drop=True))
    if len(g) <= max_effects:
        return g, 0
    top = g.head(max_effects).copy()
    rest = g.iloc[max_effects:]
    other = pd.DataFrame({'Effect': [f'Other ({len(rest)} types)'],
                          'Count': [int(rest['Count'].sum())]})
    return pd.concat([top, other], ignore_index=True), len(rest)


def survival(vals):
    v = np.sort(np.asarray(vals, dtype=float))
    v = v[np.isfinite(v)]
    if not len(v):
        return np.array([]), np.array([])
    x = np.unique(v)
    return x, 1.0 - np.searchsorted(v, x, side='left') / len(v)


# ── Panels ────────────────────────────────────────────────────────────────────
def panel_impact(ax, impact_totals, total):
    """(a) the whole call set by impact: a ring plus an aligned legend table.

    MODIFIER takes the great majority of a whole-genome call set, so the three
    ELIGIBLE classes are near-coincident slivers on the ring. Leader-line labels
    for them stack at the same angle and run into the neighbouring panel however
    they are spaced, so the ring carries only what it can show honestly --
    MODIFIER, UNANNOTATED and one accent wedge for everything eligible -- and the
    numbers move into a table beside or below it, the eligible classes indented
    under it in their own colours.

    THE LAYOUT IS CHOSEN BY MEASUREMENT. The table's three columns are sized from
    the rendered width of the strings that go in them; the type shrinks (never
    below the figure's floor) and, when the panel is too narrow to hold the table
    beside the ring, the table moves BELOW it. A fixed split overhangs the panel
    as soon as a count gains a digit or a neighbouring panel's tick labels widen
    the gutter -- and an overhang lands on those labels and reads as a fault.
    The total is the table's first row rather than text in the ring's hole, which
    cannot fit a seven-digit number at any readable size.

    Returns (n_rows_drawn, eligible_total).
    """
    fig = ax.figure
    fs0 = plt.rcParams['legend.fontsize']
    n_elig = int(sum(impact_totals.get(i, 0) for i in ELIGIBLE))

    # The ring: the excluded classes, then everything eligible as one wedge.
    ring = [(i, impact_totals.get(i, 0)) for i in ('MODIFIER', 'UNANNOTATED')
            if impact_totals.get(i, 0) > 0]
    ring.append(('eligible', n_elig))

    # The table: total, the ring wedges, then the eligible classes indented.
    rows = [(0, 'call set', total, S.INK, True, False)]
    for k, v in ring:
        col = ELIGIBLE_COLOR if k == 'eligible' else IMPACT_COLORS[k]
        rows.append((0, k, v, col, k == 'eligible', True))
        if k == 'eligible':
            for i in ELIGIBLE:
                if impact_totals.get(i, 0) > 0:
                    rows.append((1, i, impact_totals[i], IMPACT_COLORS[i], False, True))
    n = len(rows)
    labs = [('  ' if d else '') + lab for d, lab, *_ in rows]
    counts = [f'{v:,}' for _d, _l, v, *_ in rows]
    pcts = [f'{v / total:.1%}' for _d, _l, v, *_ in rows]

    box_w, box_h = box_inches(ax)
    sw_in, gap_in, indent_in, pad_in = 0.070, 0.075, 0.085, 0.10
    r_min = 0.15
    for step in range(11):
        fs = max(fs0 - 0.25 * step, V.MIN_FONT)
        sizes = text_sizes_in(fig, labs + counts + pcts, fs, 'bold')
        w_lab = max(w for w, _h in sizes[:n]) + indent_in
        w_cnt = max(w for w, _h in sizes[n:2 * n])
        w_pct = max(w for w, _h in sizes[2 * n:])
        line_h = max(h for _w, h in sizes) * 1.42
        tab_w = sw_in * 1.9 + w_lab + gap_in + w_cnt + gap_in + w_pct
        tab_h = n * line_h
        beside = tab_w + 2.0 * r_min + pad_in <= box_w
        stacked = tab_w <= box_w and tab_h + 2.0 * r_min + pad_in <= box_h
        if beside or stacked or fs <= V.MIN_FONT:
            break

    if beside:
        r_in = max(min(0.42 * box_h, (box_w - tab_w - pad_in) / 2.0), r_min)
        cx, cy = r_in, box_h / 2.0                       # ring centre, inches
        x_tab, y_top = 2.0 * r_in + pad_in, box_h / 2.0 + tab_h / 2.0 - line_h / 2.0
    else:
        r_in = max(min((box_h - tab_h - pad_in) / 2.0, 0.30 * box_w), r_min)
        cx, cy = box_w / 2.0, box_h - r_in - 0.02
        x_tab = max((box_w - tab_w) / 2.0, 0.0)
        y_top = box_h - 2.0 * r_in - pad_in - line_h / 2.0
    upi = 1.0 / r_in                                     # data units per inch; ring r = 1

    def dx(x_in):
        return (x_in - cx) * upi

    def dy(y_in):
        return (y_in - cy) * upi

    ax.pie([v for _k, v in ring], labels=None, startangle=90, counterclock=False,
           colors=[ELIGIBLE_COLOR if k == 'eligible' else IMPACT_COLORS[k] for k, _v in ring],
           wedgeprops=dict(width=0.42, edgecolor='white', linewidth=0.8),
           radius=1.0, center=(0.0, 0.0))
    # pie() imposes an equal aspect and its own limits; both are taken back so the
    # ring stays circular through the inch-based scale while the axes fills a
    # wide, short cell.
    ax.set_aspect('auto')
    ax.set_xlim(dx(0.0), dx(box_w))
    ax.set_ylim(dy(0.0), dy(box_h))

    x_cnt = x_tab + sw_in * 1.9 + w_lab + gap_in + w_cnt
    x_pct = x_cnt + gap_in + w_pct
    for i, (depth, label, value, colour, bold, swatch) in enumerate(rows):
        y = dy(y_top - i * line_h)
        x = dx(x_tab + depth * indent_in)
        if swatch:
            ax.add_patch(Rectangle((x, y - sw_in * upi / 2), sw_in * upi, sw_in * upi,
                                   facecolor=colour, edgecolor='none', zorder=6,
                                   clip_on=False))
        weight = 'bold' if bold else 'normal'
        ax.text(x + sw_in * 1.9 * upi, y, label, ha='left', va='center', fontsize=fs,
                fontweight=weight, color=colour if bold else S.INK, zorder=6, clip_on=False)
        ax.text(dx(x_cnt), y, f'{value:,}', ha='right', va='center', fontsize=fs,
                fontweight=weight, color=S.INK, zorder=6, clip_on=False)
        ax.text(dx(x_pct), y, f'{value / total:.1%}', ha='right', va='center',
                fontsize=fs, fontweight=weight, color=S.INK if bold else S.INK_SOFT,
                zorder=6, clip_on=False)
    return n, n_elig


def panel_effects(ax, anno):
    """(b) the consequences of the three eligible impacts, one block each."""
    groups, rows = [], []          # rows: (impact, effect, count, pct, is_other)
    meta = {}
    for imp in BAR_IMPACTS:
        sub = anno[anno['Impact'] == imp]
        if not len(sub):
            meta[imp] = (0, 0, 0)
            continue
        folded, n_other = fold_effects(sub)
        tot = int(sub['Count'].sum())
        members = []
        for _, r in folded.iterrows():
            is_other = str(r['Effect']).startswith('Other (')
            rows.append((imp, str(r['Effect']), int(r['Count']), r['Count'] / tot * 100,
                         is_other))
            members.append(effect_label(r['Effect']))
        # A short form shared by two rows of one block would make them read as
        # the same consequence; spell those out instead.
        dup = {m for m in members if members.count(m) > 1}
        if dup:
            members = [distinguished_label(e) if m in dup else m
                       for m, e in zip(members, folded['Effect'])]
        groups.append((imp, members))
        meta[imp] = (tot, n_other, len(folded))
    if not rows:
        empty_panel(ax, 'no eligible-impact variant in this call set')
        return [], [], [], meta
    # Table-style blocks: each impact is a bold, coloured header row over its
    # consequences, so the left margin carries one label per row and nothing
    # beside the block that could run into the neighbouring panels.
    labels, heads, ys_all, colours_lab = [], [], [], []
    y = sum(len(m) for _g, m in groups) + len(groups) - 1
    ys, block_bounds = [], []
    for gi, (imp, members) in enumerate(groups):
        top = y
        labels.append(imp); heads.append(True); colours_lab.append(IMPACT_COLORS[imp])
        ys_all.append(y); y -= 1
        for m in members:
            labels.append(m); heads.append(False); colours_lab.append(S.INK)
            ys_all.append(y); ys.append(y); y -= 1
        block_bounds.append((gi, top, y + 1))
    n_all = len(labels)
    ax.set_yticks(ys_all)
    ax.set_yticklabels(labels)
    for lab, h, col in zip(ax.get_yticklabels(), heads, colours_lab):
        lab.set_fontweight('bold' if h else 'normal')
        lab.set_color(col)
    ax.set_ylim(-0.7, n_all - 0.3)
    for gi, top, bottom in block_bounds:
        if gi % 2 == 1:
            ax.axhspan(bottom - 0.5, top + 0.5, color=S.NEUTRAL, alpha=0.13, lw=0, zorder=0)
    counts = np.array([r[2] for r in rows], dtype=float)
    for y, (imp, _e, c, _p, is_other) in zip(ys, rows):
        if is_other:
            ax.barh(y, c, height=0.62, color='white', edgecolor=IMPACT_COLORS[imp],
                    linewidth=0.9, hatch='///', zorder=3)
        else:
            ax.barh(y, c, height=0.62, color=IMPACT_COLORS[imp], zorder=3)
    ax.set_xlim(0, float(counts.max()) * 1.02)
    # Superscript count labels are wide; four ticks is what a 3 in axis holds.
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.xaxis.set_major_formatter(S.COUNT_FMT)
    ax.set_xlabel('variants (most severe annotation)')
    S.despine(ax, grid_axis='x')
    ax.tick_params(axis='y', length=0)
    texts = [f'{int(c):,} ({p:.1f}%)' for _i, _e, c, p, _o in rows]
    return counts, list(ys), texts, meta


def stacked_stratum(s):
    lab = V.stratum_label(s)
    if len(lab) > 12 and '+' in lab:
        head, tail = lab.rsplit('+', 1)
        return f'{head}\n+{tail}'
    return lab


def panel_strata(ax, order, dens, idx, colours):
    """(c) gene sets per stratum: produced by the map, and tested."""
    xs = np.arange(len(order))
    w = 0.36
    before = np.array([dens[s]['n_genes_mapped'] + dens[s]['n_genes_dropped_lt_minvar']
                       for s in order], dtype=float)
    tested = np.array([dens[s]['n_genes_mapped'] for s in order], dtype=float)
    ax.bar(xs - w / 2, before, width=w, color='white', edgecolor=S.NEUTRAL_D, linewidth=1.0,
           zorder=3)
    ax.bar(xs + w / 2, tested, width=w, color=[colours[s] for s in order], zorder=3)
    ax.set_xticks(xs)
    # Two lines for the stratum name so two groups fit a half-width panel.
    ax.set_xticklabels([f'{stacked_stratum(s)}\n{dens[s]["n_variants"]:,} variants'
                        for s in order])
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_ylim(0, float(before.max()) * 1.18 if len(before) else 1.0)
    ax.yaxis.set_major_formatter(S.INT_FMT)
    ax.set_ylabel('gene sets')
    S.despine(ax, grid_axis='y')
    ax.tick_params(axis='x', length=0)
    xs_lab = np.concatenate([xs - w / 2, xs + w / 2])
    ys_lab = np.concatenate([before, tested])
    texts = [f'{int(v):,}' for v in ys_lab]
    return xs_lab, ys_lab, texts


def panel_var_per_gene(ax, order, dens, idx, colours):
    """(d) the variants-per-set distribution, one survival curve per stratum."""
    drawn = 0
    for s in order:
        gi = idx[s]
        if not len(gi):
            continue
        x, y = survival(gi['n_var_map'].to_numpy(float))
        if not len(x):
            continue
        ax.plot(x, y, drawstyle='steps-pre', color=colours[s], lw=1.3, zorder=4,
                label=V.stratum_label(s))
        drawn += 1
    if not drawn:
        empty_panel(ax, 'no stratum retained a set')
        return 0
    mins = sorted({d['min_num_var'] for d in dens.values()})
    if len(mins) == 1:
        ax.axvline(mins[0], color=S.REFERENCE, ls='--', lw=0.9, zorder=2)
        ax.annotate(f'MinNumVar = {mins[0]}', xy=(mins[0], 0.02),
                    xycoords=('data', 'axes fraction'), xytext=(3, 0),
                    textcoords='offset points', rotation=90, ha='left', va='bottom',
                    color=S.REFERENCE, fontsize=V.MIN_FONT)
    ax.set_xscale('log')
    ax.set_ylim(0, 1.04)
    ax.set_xlabel('variants in the set (log scale)')
    ax.set_ylabel('fraction of tested sets\nwith at least this many')
    S.despine(ax)
    return drawn


# ── Tables ────────────────────────────────────────────────────────────────────
def write_annotation_summary(anno, out_dir):
    d = anno.copy()
    total = int(d['Count'].sum())
    within = d.groupby('Impact')['Count'].transform('sum')
    d['percent'] = (d['Count'] / total * 100).round(4)
    d['percent_within_impact'] = (d['Count'] / within * 100).round(4)
    rank = {k: i for i, k in enumerate(IMPACT_ORDER)}
    d = (d.assign(_r=d['Impact'].map(rank))
          .sort_values(['_r', 'Count', 'Effect', 'Biotype'], ascending=[True, False, True, True])
          .drop(columns='_r'))
    path = Path(out_dir) / 'annotation_summary.tsv'
    d.to_csv(path, sep='\t', index=False)
    return path, d


def write_stratum_summary(order, dens, idx, out_dir):
    rows = []
    for s in order:
        d, gi = dens[s], idx[s]
        v = gi['n_var_map'].to_numpy(float) if len(gi) else np.array([])
        q = np.percentile(v, [0, 25, 50, 75, 100]) if len(v) else [np.nan] * 5
        rows.append({
            'cohort': d['cohort'], 'stratum': s, 'min_num_var': d['min_num_var'],
            'n_genes_before_filter': d['n_genes_mapped'] + d['n_genes_dropped_lt_minvar'],
            'n_genes_mapped': d['n_genes_mapped'],
            'n_genes_dropped_lt_minvar': d['n_genes_dropped_lt_minvar'],
            'n_gene_variant_pairs': int(v.sum()) if len(v) else 0,
            'n_variants': d['n_variants'], 'n_variants_multigene': d['n_variants_multigene'],
            'n_sets_chrom_split': d['n_sets_chrom_split'],
            'n_positions_colliding': d['n_positions_colliding'],
            'n_chrom': int(gi['chrom'].nunique()) if len(gi) else 0,
            'var_per_set_min': q[0], 'var_per_set_p25': q[1], 'var_per_set_median': q[2],
            'var_per_set_p75': q[3], 'var_per_set_max': q[4],
            'alpha': d['alpha'], 'bonferroni_threshold': d['bonferroni_threshold'],
        })
    out = pd.DataFrame(rows)
    path = Path(out_dir) / 'stratum_summary.tsv'
    out.to_csv(path, sep='\t', index=False)
    return path, out


def main():
    args = parse_args()
    S.setup_style('paper')
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    anno = read_anno_counts(args.anno_counts)
    qc, n_chrom_qc = read_map_qc(args.map_qc)
    dens = read_denominators(args.denominator)
    idx = read_gene_indexes(args.gene_index, dens)
    # Widest stratum first everywhere in this component; the ramp is by nesting.
    order = sorted(dens, key=lambda s: (-dens[s]['n_genes_mapped'], s))
    colours = V.stratum_colors(order)
    clabel = V.cohort_labels(V.COHORT_PREFERENCE).get(args.cohort, args.cohort)

    impact_totals = {k: int(v) for k, v in anno.groupby('Impact')['Count'].sum().items()}
    total = int(anno['Count'].sum())
    if total <= 0:
        raise SystemExit(f'ABORT: {args.anno_counts} totals {total} variant(s)')
    n_eligible = int(sum(impact_totals.get(i, 0) for i in ELIGIBLE))
    pct_eligible = n_eligible / total * 100 if total else 0.0
    qc_gap = total - qc['n_variants_in']
    if qc_gap:
        print(f'[plot_annotation] WARNING: anno_counts totals {total:,} but map_qc TOTAL '
              f'n_variants_in = {qc["n_variants_in"]:,} (difference {qc_gap:+,})',
              file=sys.stderr)
    anno_path, anno_out = write_annotation_summary(anno, out_dir)
    strat_path, strat_out = write_stratum_summary(order, dens, idx, out_dir)

    # ── layout: donut, gene sets and set sizes down the left; consequences fill
    #    the right column, whose long labels live in the gutter between columns.
    fig = plt.figure(figsize=(S.COL_DOUBLE, PLOT_H))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.0, 1.05], height_ratios=[1.25, 1.0, 1.0])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[2, 0])
    ax_b = fig.add_subplot(gs[:, 1])

    n_legend_rows, _n_elig = panel_impact(ax_a, impact_totals, total)
    b_xs, b_ys, b_texts, bar_meta = panel_effects(ax_b, anno)
    c_xs, c_ys, c_texts = panel_strata(ax_c, order, dens, idx, colours)
    n_curves = panel_var_per_gene(ax_d, order, dens, idx, colours)

    # The gutter must hold (b)'s widest tick label plus its impact name; measured.
    tick_fs = plt.rcParams['ytick.labelsize']
    labels_b = [t.get_text() for t in ax_b.get_yticklabels()] or ['']
    w_tick = max(w for w, _h in text_sizes_in(fig, labels_b, tick_fs))
    w_name = max(w for w, _h in text_sizes_in(fig, BAR_IMPACTS, tick_fs, 'bold'))
    # A first estimate only: widen_gap() measures the achieved gap after the
    # layout is final and corrects it.
    ax_w_est = S.COL_DOUBLE * 0.34
    wspace = (w_tick + 0.16) / ax_w_est

    mins = sorted({d['min_num_var'] for d in dens.values()})
    minvar = str(mins[0]) if len(mins) == 1 else 'per stratum'
    S.caption_block(
        fig, plot_h=PLOT_H,
        title=(f'{clabel}: {pct_eligible:.1f} % of {total:,} minAC-passing variants can enter '
               f'a gene set.'),
        panels=['Call set by snpEff impact; eligible = HIGH + MODERATE + LOW.',
                f'Consequences by impact, top {MAX_EFFECTS} each.',
                f'Gene sets per stratum: produced by the map (hollow), tested at ≥{minvar} '
                'variants (filled).',
                'Set-size distribution of the tested sets.'],
        notes=(f'{args.cohort}; {total:,} variants over {n_chrom_qc} chromosome(s), each counted '
               f'once at its most severe annotation. Sets need ≥{minvar} variants. In (b) a '
               f'compound consequence is shown by its first (most severe) component with the '
               f'number of further components in brackets — splice_donor_variant (+1) is '
               f'splice_donor_variant&intron_variant. Where two would read alike the '
               f'differing component is spelled out; full strings in annotation_summary.tsv. '
               f'README.md'),
        letters='abcdefghijklmnop',
        top_pad=0.55, hspace=0.62, wspace=wspace, left='auto', right=0.985,
        margin_axes=[ax_a, ax_c, ax_d])

    # ---- measured layout from here on.
    wspace = widen_gap(fig, ax_a, ax_b, wspace)
    # (a)'s ring radius and table geometry are measured from the rendered axes
    # box, which caption_block and widen_gap have just changed; draw it once more
    # on a cleared axes rather than trusting the pre-caption geometry.
    ax_a.clear()
    ax_a.set_axis_off()
    n_legend_rows, _n_elig = panel_impact(ax_a, impact_totals, total)
    lab_fs = V.MIN_FONT
    if len(b_xs):
        fit_label_room(ax_b, b_xs, [w for w, _h in text_sizes_in(fig, b_texts, lab_fs)])
        S.value_labels(ax_b, b_xs, b_ys, b_texts, axis='x', offset=3.5, expand=False,
                       ha='left', va='center', fontsize=lab_fs, color=S.INK_SOFT)
    S.value_labels(ax_c, c_xs, c_ys, c_texts, axis='y', offset=2.5, expand=True,
                   fontsize=lab_fs, color=S.INK_SOFT)
    # (c) carries no legend: hollow = produced, filled = tested is said in the
    # caption, and a key would sit on the bars or on the value labels.
    if n_curves:
        S.legend_inside(ax_d, [Line2D([], [], color=colours[s], lw=1.6,
                                      label=V.stratum_label(s))
                               for s in order if len(idx[s])], loc='upper right')
    S.panel_tag(ax_a, 'a', 'call set by impact')
    S.panel_tag(ax_b, 'b', 'consequences by impact')
    S.panel_tag(ax_c, 'c', 'gene sets per stratum')
    S.panel_tag(ax_d, 'd', 'variants per tested set')
    ax_b.xaxis.set_major_locator(FixedLocator(list(ax_b.get_xticks())))
    S.thin_tick_labels(ax_b, 'x')
    fig.savefig(args.out_png)
    plt.close(fig)

    values = [
        ('cohort', args.cohort),
        ('variants in the call set', total),
        ('chromosomes', n_chrom_qc),
        ('annotated', qc['n_annotated']), ('unannotated', qc['n_unannotated']),
        ('MODIFIER-only', qc['n_modifier_only']),
        ('HIGH', impact_totals.get('HIGH', 0)), ('MODERATE', impact_totals.get('MODERATE', 0)),
        ('LOW', impact_totals.get('LOW', 0)), ('MODIFIER', impact_totals.get('MODIFIER', 0)),
        ('eligible (HIGH+MODERATE+LOW)', n_eligible),
        ('eligible percent', round(pct_eligible, 3)),
        ('distinct consequence terms', int(anno['Effect'].nunique())),
        ('multi-gene variants in the map', qc['n_multigene_variants']),
        ('MinNumVar', minvar),
        ('legend rows in (a)', n_legend_rows),
    ]
    for s in order:
        d = dens[s]
        values += [(f'{V.stratum_label(s)}: sets produced',
                    d['n_genes_mapped'] + d['n_genes_dropped_lt_minvar']),
                   (f'{V.stratum_label(s)}: sets tested', d['n_genes_mapped']),
                   (f'{V.stratum_label(s)}: variants tested', d['n_variants']),
                   (f'{V.stratum_label(s)}: Bonferroni', d['bonferroni_threshold']),
                   (f'{V.stratum_label(s)}: median variants per set',
                    float(np.median(idx[s]['n_var_map'])) if len(idx[s]) else None)]
    stats_path = args.out_stats or str(Path(args.out_png).with_suffix('.stats.json'))
    figure_doc.write_stats(args.out_png, peak=Path(args.out_png).stem, values=values)
    if args.out_stats and Path(args.out_stats).resolve() != \
            Path(args.out_png).with_suffix('.stats.json').resolve():
        Path(args.out_png).with_suffix('.stats.json').rename(args.out_stats)
    print(f'[plot_annotation] {args.cohort}: {total:,} variant(s), {n_eligible:,} eligible '
          f'({pct_eligible:.1f} %), '
          + ', '.join(f'{s}={dens[s]["n_genes_mapped"]:,} sets' for s in order)
          + f' -> {args.out_png}, {stats_path}')
    print(f'    {anno_path}  ({len(anno_out):,} row(s));  {strat_path}  ({len(strat_out):,} row(s))')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:                       # noqa: BLE001
        print(f'Error in plot_annotation: {e}', file=sys.stderr)
        sys.exit(1)
