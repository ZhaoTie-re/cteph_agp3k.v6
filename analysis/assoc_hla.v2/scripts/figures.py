#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The component's figures, one family per question, in paper order.
#           Main figures use the random model (the primary model); the fixed
#           model appears only in 07.models.
#
#   01.landscape    where in the MHC, and in how many cohorts. Per case set,
#                   two encodings of the same tests:
#                     bars.<case_set>.png      one vertical bar per test from
#                                              its smallest to its largest
#                                              -log10 P over the three cohorts,
#                                              coloured by n/3 at GEC
#                     manhattan.<case_set>.png the three cohorts as rows
#   02.signals      the independent signals: OR in every case-set comparison
#                   and cohort, and significance at each threshold
#   03.conditioning conditioning.png  each walk with a signal, before and after
#                   conditioning on its lead(s); carriers.png  the alleles of
#                   the gene that carry each residue signal
#   04.gwas_link    the genome-wide MHC lead against the HLA signal
#   05.dvt          log OR with DVT against without DVT, every test
#   06.thresholds   (supplementary) tests, Me and the family-wise error of each
#                   threshold
#   07.models       (supplementary) fixed against random model
#   08.qq           (supplementary) observed against expected -log10 P
#
#           ENCODINGS shared by every figure: n/3 at the LD-adjusted (GEC)
#           threshold as one hue light -> dark (vocab.CONS_*), the three
#           thresholds as solid / dashed / dotted lines, cohorts as an ordered
#           ramp, and one point per TEST (a residue mirror pair is drawn once).
#           The x axis of the region plots is genes in chromosomal order, then
#           markers within a gene; not a physical coordinate.
#
#           Every PNG gets a <name>.stats.json that catalogue.py turns into the
#           family README.
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402
from matplotlib.collections import LineCollection   # noqa: E402
from matplotlib.gridspec import GridSpec            # noqa: E402
from matplotlib.lines import Line2D                 # noqa: E402
from matplotlib.transforms import blended_transform_factory   # noqa: E402
from scipy import stats                  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_shared' / 'scripts'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_style as S                   # noqa: E402
import figure_doc                        # noqa: E402
import vocab as V                        # noqa: E402

LETTERS = 'abcdefghijklmnopqrstuvwxyz'
PM = V.PRIMARY_MODEL
MIN_BAND = {'allele': 9, 'residue': 20}          # x units reserved per gene, for its label
COMPARISON = {'all': 'all cases vs controls', 'dvt_yes': 'DVT yes vs controls',
              'dvt_no': 'DVT no vs controls', 'dvt_case_only': 'DVT yes vs DVT no'}
COH_ABBR = {'narrow_mainland': 'N', 'intermediate_mainland': 'I', 'full_mainland': 'F'}
FAMILIES = ('landscape', 'signals', 'conditioning', 'gwas_link', 'dvt', 'thresholds', 'models', 'qq')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--family', required=True, choices=FAMILIES)
    add_data_args(p)
    p.add_argument('--covar-label', default='sex + PC1-PC10')
    p.add_argument('--out-dir', default='.')
    return p.parse_args()


def add_data_args(p):
    p.add_argument('--assoc', nargs='+', required=True, help="every case set's assoc_all.tsv")
    p.add_argument('--summary', nargs='+', required=True, help="every case set's marker_summary.tsv")
    p.add_argument('--thresholds', nargs='+', required=True, help='every effective_tests.<class>.tsv')
    p.add_argument('--stepwise', nargs='+', required=True, help='every stepwise.<class>.tsv')
    p.add_argument('--signals', required=True)
    p.add_argument('--steps', required=True)
    p.add_argument('--lead-given', required=True)
    p.add_argument('--hla-given', required=True)
    p.add_argument('--gwas-leads', required=True)
    p.add_argument('--coords', required=True)


# ═══════════════════════════════════════════════════════════════ data
class Data:
    """Every table the figures and their READMEs read, loaded once."""

    def __init__(self, a):
        cat = lambda fs, **kw: pd.concat([pd.read_csv(f, sep='\t', **kw) for f in fs], ignore_index=True)  # noqa: E731
        self.assoc = cat(a.assoc, dtype={'marker': str, 'test_id': str})
        self.summ = cat(a.summary, dtype={'marker': str})
        self.thr = cat(a.thresholds)
        self.sw = cat(a.stepwise, dtype={'marker': str})
        self.sig = pd.read_csv(a.signals, sep='\t', dtype={'lead': str})
        self.steps = pd.read_csv(a.steps, sep='\t', dtype={'lead': str})
        self.lg = pd.read_csv(a.lead_given, sep='\t', dtype={'conditioned_on': str})
        self.hg = pd.read_csv(a.hla_given, sep='\t', dtype={'marker': str})
        gl = pd.read_csv(a.gwas_leads, sep='\t')
        self.gwas = gl[gl.status == 'kept'].reset_index(drop=True)
        self.coords = pd.read_csv(a.coords, sep='\t', dtype={'gene': str})
        self.gorder = dict(zip(self.coords.gene, self.coords.order))
        rank = {c: i for i, c in enumerate(V.CASE_SETS)}
        crank = {'residue': 0, 'allele': 1}
        s = self.sig[self.sig.model == PM].copy()
        s['_o'] = [(rank[c], crank[k], n) for c, k, n in zip(s.case_set, s['class'], s.signal)]
        self.signals = s.sort_values('_o').drop(columns='_o').reset_index(drop=True)

    def primary(self, cs, cls=None, cohort=None):
        d = self.assoc[(self.assoc.case_set == cs) & (self.assoc.model == PM)]
        if cls:
            d = d[d['class'] == cls]
        if cohort:
            d = d[d.cohort == cohort]
        return d

    def tests(self, cs, cls):
        """Markers that are their own test in at least one cohort, in plot order."""
        d = self.primary(cs, cls)
        names = set(d.loc[d.test_id == d.marker, 'marker'])
        return sorted(names, key=lambda n: (self.gorder[V.parse_marker(n, cls)[0]], V.sort_key(n, cls)))

    def rep(self, cs, cls, marker):
        """The plotted test standing for a marker (its mirror partner if that is the test)."""
        d = self.primary(cs, cls)
        t = d.loc[d.marker == marker, 'test_id'].dropna()
        return t.iloc[0] if len(t) else marker

    def threshold(self, cs, cls, kind, cohort=None):
        t = self.thr[(self.thr.case_set == cs) & (self.thr['class'] == cls)]
        if cohort:
            return float(t.loc[t.cohort == cohort, f'p_{kind}'].iloc[0])
        return float(t[f'p_{kind}'].median())

    def consistency(self, cs):
        return self.summ[self.summ.case_set == cs].set_index('marker')


def nlp(p):
    return -np.log10(np.clip(np.asarray(p, dtype=float), 1e-300, 1.0))


def majority(sig, orr):
    up = sum(1 for s, o in zip(sig, orr) if s and o > 1)
    down = sum(1 for s, o in zip(sig, orr) if s and o < 1)
    return max(up, down)


# ═══════════════════════════════════════════════════════════════ shared drawing
def xmap(names, cls, D):
    """x per marker (genes in chromosomal order, each at least MIN_BAND wide) and the gene bands."""
    by = {}
    for n in names:
        by.setdefault(V.parse_marker(n, cls)[0], []).append(n)
    xpos, bands, start = {}, [], 0.0
    for g in sorted(by, key=D.gorder.get):
        ms = by[g]
        width = max(len(ms), MIN_BAND[cls])
        step = width / len(ms)
        for i, n in enumerate(ms):
            xpos[n] = start + (i + 0.5) * step
        bands.append((g, start, start + width))
        start += width
    return xpos, bands


def gene_axis(ax, bands, labels=True):
    for i, (_g, lo, hi) in enumerate(bands):
        if i % 2 == 0:
            ax.axvspan(lo, hi, color=S.NEUTRAL, alpha=0.35, lw=0, zorder=0)
    ax.set_xlim(0, bands[-1][2])
    ax.set_xticks([(lo + hi) / 2 for _g, lo, hi in bands])
    ax.set_xticklabels([g for g, _lo, _hi in bands] if labels else [], rotation=90, fontstyle='italic',
                       fontsize=6.5)
    ax.tick_params(axis='x', length=0)
    ax.grid(False, axis='x')


def gwas_x(D, bands):
    """x of the GWAS lead among the gene bands: inside a gene, its band centre; else the boundary."""
    if D.gwas.empty:
        return None
    pos = float(D.gwas.pos.iloc[0])
    span = D.coords.set_index('gene')
    prev = 0.0
    for g, lo, hi in bands:
        s, e = span.at[g, 'gene_start'], span.at[g, 'gene_end']
        if s <= pos <= e:
            return (lo + hi) / 2
        if pos < s:
            return prev
        prev = hi
    return prev


def mark_gwas(ax, x):
    if x is None:
        return
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    ax.axvline(x, color=S.LEAD_COLOR, lw=0.7, ls=(0, (1, 2)), zorder=1)
    ax.plot([x], [1.0], marker='v', color=S.LEAD_COLOR, markersize=5, transform=tr, clip_on=False, zorder=7)


def thr_lines(ax, D, cs, cls, cohort=None, orient='h'):
    out = {}
    for t in V.THRESHOLDS:
        v = D.threshold(cs, cls, t, cohort)
        out[t] = v
        (ax.axhline if orient == 'h' else ax.axvline)(nlp(v), color=V.THRESHOLD_COLOR[t], lw=0.8,
                                                     ls=V.THRESHOLD_LS[t], zorder=2)
    return out


def range_bars(ax, xs, lo, hi, n):
    for k in (0, 1, 2, 3):
        sel = n == k
        if not sel.any():
            continue
        segs = [[(x, a), (x, b)] for x, a, b in zip(xs[sel], lo[sel], hi[sel])]
        ax.add_collection(LineCollection(segs, colors=V.CONS_COLOR[k], linewidths=V.CONS_WIDTH[k],
                                         zorder=3 + k, capstyle='round'))
        ax.scatter(xs[sel], hi[sel], s=V.CONS_SIZE[k], color=V.CONS_COLOR[k], lw=0, zorder=3 + k,
                   rasterized=k == 0)


def label(ax, x, y, text, dy=4):
    """Direct label beside a point; written to its left when the point is near the right edge."""
    lo, hi = ax.get_xlim()
    right = (x - lo) > 0.78 * (hi - lo)
    ax.annotate(text, (x, y), xytext=(-3 if right else 3, dy), textcoords='offset points',
                ha='right' if right else 'left', fontsize=6.5, color=S.INK, zorder=9,
                bbox=dict(fc='white', ec='none', alpha=0.85, pad=0.4))


def cons_handles():
    return [Line2D([], [], color=V.CONS_COLOR[k], lw=V.CONS_WIDTH[k] + 0.6, marker='o',
                   markersize=np.sqrt(V.CONS_SIZE[k]) * 0.9, label=V.CONS_LABEL[k]) for k in (3, 2, 1, 0)]


def thr_handles():
    return [Line2D([], [], color=V.THRESHOLD_COLOR[t], ls=V.THRESHOLD_LS[t], lw=0.9, label=V.THRESHOLD_LABEL[t])
            for t in V.THRESHOLDS]


def lead_handle():
    return Line2D([], [], color=S.LEAD_COLOR, marker='v', ls=(0, (1, 2)), lw=0.7, markersize=5,
                  label='GWAS MHC lead (position)')


def save(fig, out_dir, name, values):
    png = Path(out_dir) / f'{name}.png'
    fig.savefig(png, dpi=300)
    plt.close(fig)
    figure_doc.write_stats(png, peak=name, values=values)
    print(f'[figures] {name}.png')


def per_cohort(D, cs, cls, names, col='p'):
    """markers x cohorts table of one column, random model."""
    d = D.primary(cs, cls)
    return d[d.marker.isin(names)].pivot(index='marker', columns='cohort', values=col).reindex(
        index=names, columns=list(V.COHORTS))


def signal_rows(D, cs, cls):
    return D.signals[(D.signals.case_set == cs) & (D.signals['class'] == cls)]


# ═══════════════════════════════════════════════════════════════ 01.landscape
def fam_landscape(a, D):
    for cs in V.CASE_SETS:
        landscape_bars(a, D, cs)
        landscape_manhattan(a, D, cs)


def landscape_bars(a, D, cs):
    plot_h = 5.6
    fig, axes = plt.subplots(len(V.CLASSES), 1, figsize=(S.COL_DOUBLE, plot_h + 2), squeeze=False)
    axes = axes[:, 0]
    values, n_sig = [], {}
    for r, cls in enumerate(V.CLASSES):
        ax = axes[r]
        names = D.tests(cs, cls)
        xpos, bands = xmap(names, cls, D)
        P = per_cohort(D, cs, cls, names)
        lo, hi = nlp(P.max(axis=1)), nlp(P.min(axis=1))
        n = D.consistency(cs).loc[names, 'gec_n'].to_numpy(dtype=int)
        xs = np.array([xpos[m] for m in names])
        gene_axis(ax, bands)
        range_bars(ax, xs, lo, hi, n)
        thr = thr_lines(ax, D, cs, cls)
        mark_gwas(ax, gwas_x(D, bands))
        for s in signal_rows(D, cs, cls).itertuples():
            rep = D.rep(cs, cls, s.lead)
            i = names.index(rep)
            label(ax, xs[i], hi[i], f'{V.marker_label(s.lead, cls)}  {s.gec}')
        ax.set_ylim(0, max(np.nanmax(hi), nlp(thr['bonf'])) * 1.22)
        ax.set_ylabel(f'{V.CLASS_PLURAL[cls]}\n{S.NEGLOG10P}')
        S.despine(ax, grid_axis='y')
        ax.grid(False, axis='x')
        n_sig[cls] = {k: int((n == k).sum()) for k in (3, 2, 1)}
        values += [(f'{cls}: tests drawn', len(names))]
        values += [(f'{cls}: tests {k}/3 at GEC', n_sig[cls][k]) for k in (3, 2, 1)]
        values += [(f'{cls}: thresholds GEC / max-T / Bonferroni (median over cohorts)',
                    ' / '.join(f'{thr[t]:.2e}' for t in V.THRESHOLDS))]
    found = []
    for cls in V.CLASSES:
        s = signal_rows(D, cs, cls)
        found += [f'{V.marker_label(x.lead, cls)} ({x.gec})' for x in s.itertuples()]
    title = (f'{V.CASE_SET_LABEL[cs]}: ' +
             (f"independent signal{'s' if len(found) > 1 else ''} {', '.join(found)} at the LD-adjusted threshold."
              if found else 'no test reaches the LD-adjusted threshold in any cohort.'))
    S.caption_block(
        fig, plot_h=plot_h, title=title,
        panels=[f'{V.CLASS_PLURAL[c].capitalize()}: one bar per test, from its largest to its smallest P '
                'over the three cohorts (random model).' for c in V.CLASSES],
        notes=('Colour and width: cohorts significant at the LD-adjusted (GEC) threshold with the same direction. '
               'Lines: GEC (solid), exact max-T (dashed), conventional Bonferroni (dotted), medians over cohorts. '
               'Triangle: position of the genome-wide MHC lead. README.md'),
        top_pad=0.42, hspace=0.50, left='auto', right=0.985, margin_axes=list(axes))
    for r, cls in enumerate(V.CLASSES):
        S.panel_tag(axes[r], LETTERS[r])
    S.legend_above(axes[0], cons_handles() + [lead_handle()], ncol=5, full_width=True, fontsize=6.5)
    S.legend_inside(axes[1], thr_handles(), loc='upper left', fontsize=6)
    save(fig, a.out_dir, f'bars.{cs}', values)


def landscape_manhattan(a, D, cs):
    plot_h = 6.6
    fig = plt.figure(figsize=(S.COL_DOUBLE, plot_h + 2))
    gs = GridSpec(len(V.COHORTS), 2, figure=fig, width_ratios=[1, 1.75])
    axes = [[fig.add_subplot(gs[r, k]) for k in range(2)] for r in range(len(V.COHORTS))]
    values = []
    ymax = 0.0
    geo = {cls: xmap(D.tests(cs, cls), cls, D) for cls in V.CLASSES}
    for k, cls in enumerate(V.CLASSES):
        names = D.tests(cs, cls)
        xpos, bands = geo[cls]
        for r, coh in enumerate(V.COHORTS):
            ax = axes[r][k]
            e = D.primary(cs, cls, coh)
            e = e[e.marker.isin(names)]
            x = e.marker.map(xpos).to_numpy()
            y = nlp(e.p.fillna(1.0))
            sig = e.sig_gec.to_numpy(dtype=bool)
            gene_axis(ax, bands, labels=r == len(V.COHORTS) - 1)
            ax.scatter(x[~sig], y[~sig], s=4, color=S.DATA, lw=0, zorder=3, rasterized=True)
            ax.scatter(x[sig], y[sig], s=14, color=S.ACCENT, lw=0, zorder=4)
            thr = thr_lines(ax, D, cs, cls, coh)
            mark_gwas(ax, gwas_x(D, bands))
            for s in signal_rows(D, cs, cls).itertuples():
                rep = D.rep(cs, cls, s.lead)
                row = e[e.marker == rep]
                if len(row):
                    label(ax, xpos[rep], nlp(row.p.iloc[0]), V.marker_label(s.lead, cls))
            ymax = max(ymax, np.nanmax(y), nlp(thr['bonf']))
            S.despine(ax, grid_axis='y')
            ax.grid(False, axis='x')
            if k == 0:
                ax.set_ylabel(f'{V.COHORT_LABEL[coh]}\n{S.NEGLOG10P}')
            values.append((f'{cls} {V.COHORT_SHORT[coh]}: tests significant GEC / max-T / Bonferroni',
                           ' / '.join(str(int(e[f'sig_{t}'].sum())) for t in V.THRESHOLDS)))
    for row in axes:
        for ax in row:
            ax.set_ylim(0, ymax * 1.18)
    S.caption_block(
        fig, plot_h=plot_h,
        title=f'{V.CASE_SET_LABEL[cs]}: every test in each cohort (random model), alleles left, residues right.',
        panels=[f'{V.COHORT_LABEL[c]}, {V.CLASS_PLURAL[cl]}.' for c in V.COHORTS for cl in V.CLASSES],
        notes=('Red: significant at the cohort\'s LD-adjusted (GEC) threshold. Lines: GEC (solid), exact max-T '
               '(dashed), conventional Bonferroni (dotted), the cohort\'s own values. Triangle: genome-wide MHC '
               'lead. README.md'),
        letters=LETTERS, top_pad=0.40, hspace=0.30, wspace=0.10, left='auto', right=0.985,
        margin_axes=[r[0] for r in axes])
    for r, coh in enumerate(V.COHORTS):
        for k, cls in enumerate(V.CLASSES):
            S.panel_tag(axes[r][k], LETTERS[r * 2 + k], title=V.CLASS_PLURAL[cls] if r == 0 else None, tsize=8)
    S.legend_inside(axes[0][1], [Line2D([], [], ls='none', marker='o', color=S.ACCENT, markersize=4,
                                        label='significant (GEC)')] + thr_handles(), loc='upper left', fontsize=6)
    save(fig, a.out_dir, f'manhattan.{cs}', values)


# ═══════════════════════════════════════════════════════════════ 02.signals
def fam_signals(a, D):
    sigs = D.signals
    if sigs.empty:
        figure_doc.write_stats(Path(a.out_dir) / 'signals.png', peak='signals', values=[('signals', 0)])
        print('[figures] signals: no signal; no PNG')
        return
    groups, rows = [], []
    for s in sigs.itertuples():
        name = [(V.marker_label(s.lead, s[3]), dict(fontweight='bold')),
                (f'signal in {V.CASE_SET_SHORT[s.case_set]}', dict(fontsize=6, fontweight='normal',
                                                                     fontstyle='normal', color=S.INK_SOFT))]
        groups.append((name, [COMPARISON[c] for c in V.CASE_SETS]))
        rows += [(s.lead, s[3], s.case_set, c) for c in V.CASE_SETS]
    n = len(rows)
    plot_h = max(3.2, 0.24 * n + 1.0)
    fig = plt.figure(figsize=(S.COL_DOUBLE, plot_h + 2))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[2.2, 1])
    ax, axm = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    ys = S.grouped_row_labels(ax, groups, italic=False, fontsize=6.5)
    offs = {c: o for c, o in zip(V.COHORTS, (0.22, 0.0, -0.22))}
    xs_all, values = [], []
    for (lead, cls, own, comp), y in zip(rows, ys):
        e = D.primary(comp, cls)
        e = e[e.marker == lead].set_index('cohort')
        for coh in V.COHORTS:
            if coh not in e.index or pd.isna(e.at[coh, 'p']):
                continue
            r = e.loc[coh]
            yy = y + offs[coh]
            ax.plot([r['l95'], r['u95']], [yy, yy], color=V.COHORT_COLOR[coh], lw=1.0, zorder=2)
            ax.plot(r['or'], yy, ls='none', marker='o', markersize=3.8, color=V.COHORT_COLOR[coh],
                    markerfacecolor=V.COHORT_COLOR[coh] if r['sig_gec'] else 'white', markeredgewidth=0.9, zorder=3)
            xs_all += [r['l95'], r['u95']]
        # significance matrix: thresholds x cohorts, in this comparison's own thresholds
        for g, t in enumerate(V.THRESHOLDS):
            sig = [bool(e.at[c, f'sig_{t}']) if c in e.index else False for c in V.COHORTS]
            orr = [e.at[c, 'or'] if c in e.index else np.nan for c in V.COHORTS]
            for j, (c, sg) in enumerate(zip(V.COHORTS, sig)):
                axm.plot(g * 5 + j, y, ls='none', marker='o', markersize=4.2,
                         color=S.INK if sg else S.NEUTRAL_D, markerfacecolor=S.INK if sg else 'white',
                         markeredgewidth=0.8, zorder=3)
            k = majority(sig, orr)
            axm.text(g * 5 + 2.6, y, f'{k}/3', va='center', ha='left', fontsize=6.5, color=S.INK,
                     fontweight='bold' if k else 'normal')
        values.append((f'{lead} ({V.CASE_SET_SHORT[own]} signal), {COMPARISON[comp]}: OR N / I / F',
                       ' / '.join(f"{e.at[c, 'or']:.2f}" if c in e.index else '-' for c in V.COHORTS)))
    ax.axvline(1.0, color=S.REFERENCE, lw=0.8, zorder=1)
    ax.set_xscale('log')
    lo, hi = np.nanpercentile(xs_all, 1), np.nanpercentile(xs_all, 99)
    ax.set_xlim(max(lo * 0.9, 0.08), min(hi * 1.1, 12))
    S.or_log_axis(ax, 'x')
    ax.set_xlabel(f'{S.OR_SYM} per copy (95 % CI), random model')
    S.despine(ax, grid_axis='x')
    axm.set_ylim(ax.get_ylim())
    axm.set_xlim(-0.6, 14.4)
    axm.set_xticks([g * 5 + j for g in range(3) for j in range(3)])
    axm.set_xticklabels([COH_ABBR[c] for _ in range(3) for c in V.COHORTS], fontsize=6.5)
    axm.set_yticks([])
    for sp in ('left', 'top', 'right'):
        axm.spines[sp].set_visible(False)
    axm.grid(False)
    tr = blended_transform_factory(axm.transData, axm.transAxes)
    for g, t in enumerate(V.THRESHOLDS):
        axm.annotate({'gec': 'GEC', 'maxt': 'max-T', 'bonf': 'Bonferroni'}[t], xy=(g * 5 + 1.3, 0),
                     xycoords=tr, xytext=(0, -22), textcoords='offset points', ha='center', va='top',
                     fontsize=6.5, fontweight='bold' if t == 'gec' else 'normal', annotation_clip=False)
    S.caption_block(
        fig, plot_h=plot_h,
        title=(f"{len(sigs)} independent signal{'s' if len(sigs) > 1 else ''}, random model: " +
               '; '.join(f'{V.marker_label(s.lead, s[3])} in {V.CASE_SET_SHORT[s.case_set]} ({s.gec})'
                         for s in sigs.itertuples()) + '.'),
        panels=['Odds ratio per copy of each signal\'s lead in every comparison, one point per cohort '
                '(dark narrow, mid intermediate, light full); filled = significant at the LD-adjusted threshold.',
                'Significance in each cohort (N narrow, I intermediate, F full) at each threshold, filled = '
                'significant, and the count with the same direction (n/3).'],
        notes=('Thresholds are those of each comparison, cohort and class. The DVT-yes and DVT-no comparisons '
               'share their controls; DVT yes vs DVT no (cases only) is the test of a difference. README.md'),
        top_pad=0.40, wspace=0.05, left='auto', right=0.985, margin_axes=[ax])
    S.panel_tag(ax, 'a')
    S.panel_tag(axm, 'b')
    S.legend_above(ax, [Line2D([], [], color=V.COHORT_COLOR[c], marker='o', lw=1.0, markersize=3.8,
                               label=V.COHORT_SHORT[c]) for c in V.COHORTS], ncol=3, fontsize=6.5)
    save(fig, a.out_dir, 'signals', values)


# ═══════════════════════════════════════════════════════════════ 03.conditioning
def step_bars(D, cs, cls, step, names):
    w = D.sw[(D.sw.case_set == cs) & (D.sw.model == PM) & (D.sw['class'] == cls) & (D.sw.step == step)]
    P = w.pivot(index='marker', columns='cohort', values='p').reindex(index=names, columns=list(V.COHORTS))
    O = w.pivot(index='marker', columns='cohort', values='or').reindex(index=names, columns=list(V.COHORTS))
    G = w.pivot(index='marker', columns='cohort', values='sig_gec').reindex(index=names, columns=list(V.COHORTS))
    G = G.fillna(False).astype(bool)
    n = np.array([majority(G.loc[m].tolist(), O.loc[m].tolist()) for m in names])
    return P, n


def fam_conditioning(a, D):
    walks = []
    for (cs, cls), st in D.steps[D.steps.model == PM].groupby(['case_set', 'class'], sort=False):
        if st.entered.any():
            walks.append((cs, cls, st.sort_values('step')))
    rank = {c: i for i, c in enumerate(V.CASE_SETS)}
    walks.sort(key=lambda w: (rank[w[0]], 0 if w[1] == 'residue' else 1))
    if not walks:
        figure_doc.write_stats(Path(a.out_dir) / 'conditioning.png', peak='conditioning', values=[('walks', 0)])
    else:
        conditioning_png(a, D, walks)
    carriers_png(a, D)


def conditioning_png(a, D, walks):
    plot_h = 1.75 * len(walks) + 0.5
    fig, axes = plt.subplots(len(walks), 2, figsize=(S.COL_DOUBLE, plot_h + 2), squeeze=False)
    values, panels = [], []
    for r, (cs, cls, st) in enumerate(walks):
        names = D.tests(cs, cls)
        xpos, bands = xmap(names, cls, D)
        xs = np.array([xpos[m] for m in names])
        last = int(st.step.max())
        leads = st.loc[st.entered, 'lead'].tolist()
        ymax = 0.0
        for k, step in enumerate((0, last)):
            ax = axes[r][k]
            P, n = step_bars(D, cs, cls, step, names)
            lo, hi = nlp(P.max(axis=1).fillna(1.0)), nlp(P.min(axis=1).fillna(1.0))
            gene_axis(ax, bands)
            range_bars(ax, xs, lo, hi, n)
            thr = thr_lines(ax, D, cs, cls)
            for ld in leads:
                i = names.index(D.rep(cs, cls, ld))
                if step == 0:
                    label(ax, xs[i], hi[i], V.marker_label(ld, cls))
                else:
                    ax.plot([xs[i]], [1.0], marker='v', color=S.INK, markersize=4.5, clip_on=False, zorder=8,
                            transform=blended_transform_factory(ax.transData, ax.transAxes))
            ymax = max(ymax, np.nanmax(hi), nlp(thr['bonf']))
            S.despine(ax, grid_axis='y')
            ax.grid(False, axis='x')
            nn = {j: int((n == j).sum()) for j in (3, 2, 1)}
            values.append((f'{V.CASE_SET_SHORT[cs]} {cls} step {step}: drawn tests 3/3 / 2/3 / 1/3 at GEC',
                           f'{nn[3]} / {nn[2]} / {nn[1]}'))
        for ax in axes[r]:
            ax.set_ylim(0, ymax * 1.22)
        axes[r][0].set_ylabel(f'{V.CASE_SET_SHORT[cs]}\n{V.CLASS_PLURAL[cls]}\n{S.NEGLOG10P}')
        panels += [f'{V.CASE_SET_SHORT[cs]}, {V.CLASS_PLURAL[cls]}: marginal (step 0).',
                   f'{V.CASE_SET_SHORT[cs]}, {V.CLASS_PLURAL[cls]}: given {", ".join(V.marker_label(x, cls) for x in leads)} '
                   f'(step {last}).']
    left_over = {(cs, cls): int(st.loc[st.step == st.step.max(), 'n_sig_gec_any'].iloc[0]) for cs, cls, st in walks}
    n_sig = sum(int(st.entered.sum()) for _cs, _cls, st in walks)
    title = (f"{n_sig} signal{'s' if n_sig > 1 else ''} in {len(walks)} walk{'s' if len(walks) > 1 else ''}; "
             'given the signal lead(s), ' +
             ('no test remains significant at the LD-adjusted threshold in any cohort.'
              if not any(left_over.values()) else
              ', '.join(f'{V.CASE_SET_SHORT[c]} {k}: {n} test(s) remain significant' for (c, k), n in left_over.items()
                        if n) + '.'))
    S.caption_block(
        fig, plot_h=plot_h,
        title=title,
        panels=panels,
        notes=('Random model. Bars as in 01.landscape: from the largest to the smallest P over the three cohorts, '
               'coloured by n/3 at GEC. Triangle (right): the conditioned lead. README.md'),
        letters=LETTERS, top_pad=0.40, hspace=0.75, wspace=0.12, left='auto', right=0.985,
        margin_axes=[r[0] for r in axes])
    for r in range(len(walks)):
        for k in range(2):
            S.panel_tag(axes[r][k], LETTERS[r * 2 + k],
                        title=('marginal' if k == 0 else 'conditional on the lead') if r == 0 else None, tsize=8)
    S.legend_above(axes[0][0], cons_handles(), ncol=4, full_width=True, fontsize=6.5, y=1.18)
    save(fig, a.out_dir, 'conditioning', values)


def carriers_png(a, D):
    res = D.signals[D.signals['class'] == 'residue']
    if res.empty:
        figure_doc.write_stats(Path(a.out_dir) / 'carriers.png', peak='carriers', values=[('residue signals', 0)])
        return
    k = len(res)
    blocks = []
    for s in res.itertuples():
        gene = V.parse_marker(s.lead, 'residue')[0]
        carriers = [x.strip() for x in str(s.carrier_alleles).split(',') if x.strip() and x.strip() != 'nan']
        al = D.primary(s.case_set, 'allele')
        names = sorted({m for m in al.marker if m.startswith(f'{gene}*')}, key=lambda m: V.sort_key(m, 'allele'))
        blocks.append((s, carriers, [m for m in names if m in carriers], [m for m in names if m not in carriers]))
    nrow = max(len(c) + len(o) for _s, _cs, c, o in blocks)
    plot_h = max(3.0, 0.16 * nrow + 1.2)
    fig, axes = plt.subplots(1, k, figsize=(S.COL_DOUBLE, plot_h + 2), squeeze=False)
    axes = axes[0]
    values, panels = [], []
    offs = {c: o for c, o in zip(V.COHORTS, (0.22, 0.0, -0.22))}
    for i, (s, carriers, car, oth) in enumerate(blocks):
        ax = axes[i]
        lab = V.marker_label(s.lead, 'residue')
        groups = [('carriers', car), ('non-carriers', oth)]
        groups = [(g, m) for g, m in groups if m]
        ys = S.grouped_row_labels(ax, groups, italic=False, fontsize=6)
        order = [m for _g, ms in groups for m in ms]
        al = D.primary(s.case_set, 'allele')
        xs = []
        for m, y in zip(order, ys):
            e = al[al.marker == m].set_index('cohort')
            for coh in V.COHORTS:
                if coh not in e.index or pd.isna(e.at[coh, 'p']):
                    continue
                r = e.loc[coh]
                yy = y + offs[coh]
                ax.plot([r['l95'], r['u95']], [yy, yy], color=V.COHORT_COLOR[coh], lw=0.9, zorder=2)
                ax.plot(r['or'], yy, ls='none', marker='o', markersize=3.0, color=V.COHORT_COLOR[coh],
                        markerfacecolor=V.COHORT_COLOR[coh] if r['sig_gec'] else 'white',
                        markeredgewidth=0.8, zorder=3)
                xs.append(r['or'])
        ax.axvline(1.0, color=S.REFERENCE, lw=0.8, zorder=1)
        ax.set_xscale('log')
        ax.set_xlim(0.15, 6.5)
        S.or_log_axis(ax, 'x')
        ax.set_xlabel(f'{S.OR_SYM} per copy (95 % CI)')
        S.despine(ax, grid_axis='x')
        panels.append(f'{lab}, {V.CASE_SET_SHORT[s.case_set]}: the tested {V.parse_marker(s.lead, "residue")[0]} '
                      f'alleles, those carrying the residue above.')
        values += [(f'{lab} ({V.CASE_SET_SHORT[s.case_set]}): tested carrier alleles', ', '.join(car) or '-'),
                   (f'{lab} ({V.CASE_SET_SHORT[s.case_set]}): other tested alleles', len(oth))]
    S.caption_block(
        fig, plot_h=plot_h,
        title='Which alleles carry each residue signal: odds ratio of every tested allele of the gene, by cohort.',
        panels=panels,
        notes=('Random model, the signal\'s own comparison. Filled = significant at the LD-adjusted threshold. '
               'Carriers from hla.typing.v2 allele_match_map + residue_reference. CIs beyond the axis are cut. '
               'README.md'),
        top_pad=0.55, wspace=0.75, left='auto', right=0.985, margin_axes=[axes[0]])
    for i, (s, _c, _car, _oth) in enumerate(blocks):
        S.panel_tag(axes[i], LETTERS[i], title=f"{V.marker_label(s.lead, 'residue')}, {V.CASE_SET_SHORT[s.case_set]}",
                    tsize=8)
    S.legend_above(axes[0], [Line2D([], [], color=V.COHORT_COLOR[c], marker='o', lw=0.9, markersize=3.0,
                                    label=V.COHORT_SHORT[c]) for c in V.COHORTS], ncol=3, fontsize=6.5,
                   full_width=True, y=1.10)
    save(fig, a.out_dir, 'carriers', values)


# ═══════════════════════════════════════════════════════════════ 04.gwas_link
def fam_gwas_link(a, D):
    if D.gwas.empty:
        figure_doc.write_stats(Path(a.out_dir) / 'gwas_link.png', peak='gwas_link', values=[('GWAS leads', 0)])
        return
    ldid = D.gwas.lead_id.iloc[0]
    s1 = signal_rows(D, 'all', 'residue')
    sig1 = s1.lead.iloc[0] if len(s1) else None
    lg = D.lg[(D.lg.case_set == 'all') & (D.lg.model == PM) & (D.lg.lead == ldid)]
    single = lg[~lg.why.str.contains('jointly')]
    hg = D.hg[(D.hg.case_set == 'all') & (D.hg.model == PM) & (D.hg.lead == ldid)]
    plot_h = 4.6
    fig = plt.figure(figsize=(S.COL_DOUBLE, plot_h + 2))
    gs = GridSpec(2, 6, figure=fig, height_ratios=[1, 1.25])
    axa, axb = fig.add_subplot(gs[0, 0:3]), fig.add_subplot(gs[0, 3:6])
    axc = [fig.add_subplot(gs[1, 2 * i:2 * i + 2]) for i in range(3)]
    values = []
    yrow = {c: 2 - i for i, c in enumerate(V.COHORTS)}
    tag = {}
    for coh in V.COHORTS:
        d = single[single.cohort == coh]
        if d.empty:
            continue
        y = yrow[coh]
        m0 = nlp(d.lead_p_marginal.iloc[0])
        best = d.loc[d.r.abs().idxmax()]
        tag[coh] = best.conditioned_on
        pts = [(m0, 'marginal')]
        if sig1 is not None and (d.conditioned_on == sig1).any():
            pts.append((nlp(d.loc[d.conditioned_on == sig1, 'lead_p_given'].iloc[0]), 'sig1'))
        pts.append((nlp(best.lead_p_given), 'tag'))
        axa.plot([min(p for p, _ in pts), max(p for p, _ in pts)], [y, y], color=S.NEUTRAL_D, lw=0.8, zorder=1)
        for p, kind in pts:
            st = {'marginal': dict(marker='o', color=S.INK, mfc=S.INK),
                  'sig1': dict(marker='o', color=S.ACCENT, mfc=S.ACCENT),
                  'tag': dict(marker='s', color=S.REFERENCE, mfc='white')}[kind]
            axa.plot(p, y, ls='none', markersize=5, markeredgewidth=0.9, zorder=3, **st)
        values.append((f'{V.COHORT_SHORT[coh]}: GWAS lead P marginal / given {sig1} / given best LD tag '
                       f'({best.conditioned_on}, r2 {best.r ** 2:.2f})',
                       ' / '.join(f'{10 ** -p:.1e}' for p, _ in pts)))
        if sig1 is not None:
            h = hg[(hg.cohort == coh) & (hg.marker == sig1)]
            if len(h):
                p0, p1 = nlp(h.p_marginal.iloc[0]), nlp(h.p_given_lead.iloc[0])
                axb.plot([p1, p0], [y, y], color=S.NEUTRAL_D, lw=0.8, zorder=1)
                axb.plot(p0, y, ls='none', marker='o', color=S.ACCENT, markersize=5, zorder=3)
                axb.plot(p1, y, ls='none', marker='o', color=S.ACCENT, mfc='white', markersize=5,
                         markeredgewidth=0.9, zorder=3)
                axb.plot(nlp(D.threshold('all', 'residue', 'gec', coh)), y, ls='none', marker='|', color=S.ACCENT,
                         markersize=9, zorder=2)
                values.append((f'{V.COHORT_SHORT[coh]}: {sig1} P marginal / given the GWAS lead',
                               f'{10 ** -p0:.1e} / {10 ** -p1:.1e}'))
    fmt = lambda xs: '/'.join(f'{x:.1e}' for x in xs)                     # noqa: E731
    pm = [single[single.cohort == c].lead_p_marginal.iloc[0] for c in V.COHORTS if (single.cohort == c).any()]
    ps = [single[(single.cohort == c) & (single.conditioned_on == sig1)].lead_p_given.iloc[0]
          for c in V.COHORTS if ((single.cohort == c) & (single.conditioned_on == sig1)).any()]
    pt = [single[(single.cohort == c) & (single.conditioned_on == tag.get(c))].lead_p_given.iloc[0]
          for c in V.COHORTS if c in tag]
    gwas_title = (f'All cases, GWAS MHC lead {ldid}: P {fmt(pm)} marginally, {fmt(ps)} given '
                  f'{V.marker_label(sig1, "residue")}, {fmt(pt)} given its best LD tag (narrow / intermediate / full).'
                  if sig1 and ps else f'All cases, GWAS MHC lead {ldid}: P {fmt(pm)} marginally.')
    for ax in (axa, axb):                       # rows 2..0; the legend sits in the empty band below row 0
        ax.set_yticks([yrow[c] for c in V.COHORTS])
        ax.set_yticklabels([V.COHORT_SHORT[c] for c in V.COHORTS] if ax is axa else [])
        ax.set_ylim(-2.1, 2.5)
        ax.set_xlim(left=0)
        S.despine(ax, grid_axis='x')
    axa.axvline(nlp(1e-5), color=S.INK_SOFT, ls=':', lw=0.8, zorder=0)
    axa.set_xlabel(f'GWAS lead {S.NEGLOG10P}')
    axb.set_xlabel(f'{V.marker_label(sig1, "residue") if sig1 else "signal"} {S.NEGLOG10P}')
    axa.legend(handles=[Line2D([], [], ls='none', marker='o', color=S.INK, label='marginal'),
                        Line2D([], [], ls='none', marker='o', color=S.ACCENT,
                               label=f'given {V.marker_label(sig1, "residue")}' if sig1 else 'given signal'),
                        Line2D([], [], ls='none', marker='s', color=S.REFERENCE, mfc='white',
                               label='given its best LD tag'),
                        Line2D([], [], color=S.INK_SOFT, ls=':', label='GWAS suggestive (1e-5)')],
               loc='lower left', ncol=2, fontsize=6, frameon=True, framealpha=0.92, edgecolor='#C8CDD2')
    axb.legend(handles=[Line2D([], [], ls='none', marker='o', color=S.ACCENT, label='marginal'),
                        Line2D([], [], ls='none', marker='o', color=S.ACCENT, mfc='white', label='given the GWAS lead'),
                        Line2D([], [], ls='none', marker='|', color=S.ACCENT, markersize=8, label='GEC threshold')],
               loc='lower left', ncol=2, fontsize=6, frameon=True, framealpha=0.92, edgecolor='#C8CDD2')
    for i, coh in enumerate(V.COHORTS):
        ax = axc[i]
        d = single[single.cohort == coh]
        if d.empty:
            continue
        r2, y = d.r ** 2, nlp(d.lead_p_given)
        hit = (d.conditioned_on == sig1).to_numpy() if sig1 else np.zeros(len(d), bool)
        ax.scatter(r2[~hit], y[~hit], s=10, color=S.DATA, lw=0, zorder=3)
        ax.scatter(r2[hit], y[hit], s=22, color=S.ACCENT, lw=0, zorder=4)
        ax.axhline(nlp(d.lead_p_marginal.iloc[0]), color=S.INK, ls='--', lw=0.7, zorder=1)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, nlp(d.lead_p_given).max() * 1.25)
        if hit.any():
            label(ax, r2[hit].iloc[0], y[hit][0], V.marker_label(sig1, 'residue'), dy=-10)
        t = d[d.conditioned_on == tag.get(coh)]
        if len(t):
            label(ax, (t.r ** 2).iloc[0], nlp(t.lead_p_given.iloc[0]),
                  V.marker_label(t.conditioned_on.iloc[0], 'allele' if '*' in t.conditioned_on.iloc[0] else 'residue'))
        ax.set_xlabel(r'$r^{2}$ with the GWAS lead')
        if i == 0:
            ax.set_ylabel(f'GWAS lead {S.NEGLOG10P}\ngiven one HLA marker')
        S.despine(ax)
    S.caption_block(
        fig, plot_h=plot_h,
        title=gwas_title,
        panels=['GWAS lead per cohort: marginal, given the HLA signal lead, and given its best LD tag.',
                'HLA signal lead per cohort: marginal and given the GWAS lead.',
                'Narrow mainland: the lead given each single HLA marker tried, against that marker\'s r².',
                'Intermediate mainland, as c.', 'Full mainland, as c.'],
        notes=('All cases, random model. Markers tried: the 10 in highest |r| with the lead, the 10 with the '
               'smallest marginal P and the signal leads. Dashed: the lead\'s marginal −log10 P. The lead is '
               'mean-imputed here (METHODS §7.1). README.md'),
        top_pad=0.30, hspace=0.55, wspace=1.1, left='auto', right=0.985, margin_axes=[axa, axc[0]])
    for i, ax in enumerate([axa, axb] + axc):
        S.panel_tag(ax, LETTERS[i])
    save(fig, a.out_dir, 'gwas_link', values)


# ═══════════════════════════════════════════════════════════════ 05.dvt
def fam_dvt(a, D):
    plot_h = 4.8
    fig, axes = plt.subplots(len(V.CLASSES), len(V.COHORTS), figsize=(S.COL_DOUBLE, plot_h + 2), squeeze=False)
    values, panels = [], []
    lab_markers = {m for m in D.signals.lead}
    n_tot = {'dvt_yes': 0, 'dvt_no': 0, 'dvt_case_only': 0}
    for r, cls in enumerate(V.CLASSES):
        names = sorted(set(D.tests('dvt_yes', cls)) & set(D.tests('dvt_no', cls)))
        for k, coh in enumerate(V.COHORTS):
            ax = axes[r][k]
            g = {cs: D.primary(cs, cls, coh).set_index('marker').reindex(names)
                 for cs in ('dvt_yes', 'dvt_no', 'dvt_case_only')}
            ok = g['dvt_yes']['or'].notna() & g['dvt_no']['or'].notna()
            x, y = np.log(g['dvt_no'].loc[ok, 'or']), np.log(g['dvt_yes'].loc[ok, 'or'])
            sy = g['dvt_yes'].loc[ok, 'sig_gec'].fillna(False).astype(bool)
            sn = g['dvt_no'].loc[ok, 'sig_gec'].fillna(False).astype(bool)
            so = g['dvt_case_only'].loc[ok, 'sig_gec'].fillna(False).astype(bool)
            for cs_ in n_tot:                     # counted on each set's own tests, as in its tables
                e_ = D.primary(cs_, cls, coh)
                n_tot[cs_] += int(e_.loc[e_.sig_gec.astype(bool), 'test_id'].nunique())
            nn = ~sy & ~sn
            ax.scatter(x[nn], y[nn], s=5, color=S.NEUTRAL_D, lw=0, zorder=2, rasterized=True)
            ax.scatter(x[sn], y[sn], s=14, color=V.CASE_SET_COLOR['dvt_no'], lw=0, zorder=3)
            ax.scatter(x[sy], y[sy], s=14, color=V.CASE_SET_COLOR['dvt_yes'], lw=0, zorder=4)
            ax.scatter(x[so], y[so], s=46, facecolors='none', edgecolors=S.INK, linewidths=0.9, zorder=5)
            lim = max(np.abs(np.r_[x, y]).max(), 0.5) * 1.08
            ax.plot([-lim, lim], [-lim, lim], color=S.REFERENCE, lw=0.8, zorder=1)
            ax.axhline(0, color=S.NEUTRAL_D, lw=0.6, zorder=1)
            ax.axvline(0, color=S.NEUTRAL_D, lw=0.6, zorder=1)
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            for j, m in enumerate(sorted(lab_markers & set(x.index))):
                label(ax, x[m], y[m], V.marker_label(m, cls), dy=-10 if j % 2 == 0 else 5)
            ax.set_box_aspect(1)
            S.despine(ax)
            if r == len(V.CLASSES) - 1:
                ax.set_xlabel('log OR, DVT no vs controls')
            if k == 0:
                ax.set_ylabel(f'{V.CLASS_PLURAL[cls]}\nlog OR, DVT yes vs controls')
            rr = float(np.corrcoef(x, y)[0, 1])
            ax.text(0.04, 0.96, f'r = {rr:.2f}', transform=ax.transAxes, ha='left', va='top', fontsize=6.5,
                    color=S.INK_SOFT)
            values += [(f'{cls} {V.COHORT_SHORT[coh]}: tests, r(log OR)', f'{int(ok.sum())}, {rr:.3f}'),
                       (f'{cls} {V.COHORT_SHORT[coh]}: points drawn GEC-significant in DVT yes / DVT no / case-only',
                        f'{int(sy.sum())} / {int(sn.sum())} / {int(so.sum())}')]
            panels.append(f'{V.CLASS_PLURAL[cls]}, {V.COHORT_LABEL[coh]}.')
    S.caption_block(
        fig, plot_h=plot_h,
        title=('Cases with and without DVT, each against the same controls: '
               f"{n_tot['dvt_no']} and {n_tot['dvt_yes']} test-cohort pairs significant at the LD-adjusted threshold "
               f"without and with DVT; {n_tot['dvt_case_only']} in the case-only test."),
        panels=panels,
        notes=('Random model, one point per test. Blue / red: significant at the LD-adjusted threshold in DVT no / '
               'DVT yes. Ring: significant in the case-only test. Shared controls correlate the axes even '
               'under the null. README.md'),
        letters=LETTERS, top_pad=0.45, hspace=0.40, wspace=0.30, left='auto', right=0.985,
        margin_axes=[r[0] for r in axes])
    for r in range(len(V.CLASSES)):
        for k, coh in enumerate(V.COHORTS):
            S.panel_tag(axes[r][k], LETTERS[r * 3 + k], title=V.COHORT_LABEL[coh] if r == 0 else None, tsize=8)
    S.legend_above(axes[0][1], [
        Line2D([], [], ls='none', marker='o', color=V.CASE_SET_COLOR['dvt_no'], markersize=4,
               label='significant, DVT no'),
        Line2D([], [], ls='none', marker='o', color=V.CASE_SET_COLOR['dvt_yes'], markersize=4,
               label='significant, DVT yes'),
        Line2D([], [], ls='none', marker='o', color=S.NEUTRAL_D, markersize=3, label='neither'),
        Line2D([], [], ls='none', marker='o', mfc='none', mec=S.INK, markersize=6, label='case-only')],
        ncol=4, fontsize=6.5, full_width=True, y=1.16)
    save(fig, a.out_dir, 'dvt', values)


# ═══════════════════════════════════════════════════════════════ 06.thresholds
def fam_thresholds(a, D):
    t = D.thr.copy()
    rank = {c: i for i, c in enumerate(V.CASE_SETS)}
    crank = {c: i for i, c in enumerate(V.COHORTS)}
    t['_o'] = [(rank[c], crank[h]) for c, h in zip(t.case_set, t.cohort)]
    plot_h = 4.4
    fig = plt.figure(figsize=(S.COL_DOUBLE, plot_h + 2))
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1, 1, 0.85])
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    groups = [(V.CASE_SET_SHORT[c], [V.COHORT_SHORT[h] for h in V.COHORTS]) for c in V.CASE_SETS]
    style = {'tests': dict(marker='o', color=S.INK_SOFT, mfc='white'),
             'me_maxt': dict(marker='D', color=S.REFERENCE, mfc=S.REFERENCE),
             'me_gec': dict(marker='o', color=S.ACCENT, mfc=S.ACCENT)}
    lab = {'tests': 'tests (conventional Bonferroni)', 'me_maxt': 'Me implied by exact max-T',
           'me_gec': 'Me, GEC v0.2 (LD-adjusted)'}
    values = []
    for i, cls in enumerate(V.CLASSES):
        ax = axes[i]
        d = t[t['class'] == cls].sort_values('_o')
        ys = S.grouped_row_labels(ax, groups, italic=False, fontsize=6) if i == 0 else None
        if ys is None:
            ax.set_ylim(axes[0].get_ylim())
            ax.set_yticks(axes[0].get_yticks())
            ax.set_yticklabels([])
            ys = list(axes[0].get_yticks())
        for (_, row), y in zip(d.iterrows(), ys):
            ax.plot([row.me_gec, row.tests], [y, y], color=S.NEUTRAL_D, lw=0.8, zorder=1)
            for kcol in ('tests', 'me_maxt', 'me_gec'):
                ax.plot(row[kcol], y, ls='none', markersize=4, markeredgewidth=0.8, zorder=3, **style[kcol])
        ax.set_xlabel(f'number of tests, {V.CLASS_PLURAL[cls]}')
        S.despine(ax, grid_axis='x')
        values.append((f'{cls}: tests / Me max-T / Me GEC (range over families)',
                       f'{int(d.tests.min())}–{int(d.tests.max())} / {d.me_maxt.min():.1f}–{d.me_maxt.max():.1f} / '
                       f'{d.me_gec.min():.1f}–{d.me_gec.max():.1f}'))
    ax = axes[2]
    rng = np.random.default_rng(1)
    for j, cls in enumerate(V.CLASSES):
        d = t[t['class'] == cls]
        for kk, (col, st) in enumerate((('fwer_gec', dict(color=S.ACCENT, mfc=S.ACCENT, marker='o')),
                                        ('fwer_bonf', dict(color=S.INK_SOFT, mfc='white', marker='o')))):
            x = j + (kk - 0.5) * 0.36 + rng.uniform(-0.07, 0.07, len(d))
            ax.plot(x, d[col], ls='none', markersize=3.5, markeredgewidth=0.8, zorder=3, **st)
        values.append((f'{cls}: family-wise error under max-T, GEC / Bonferroni',
                       f'{d.fwer_gec.min():.3f}–{d.fwer_gec.max():.3f} / {d.fwer_bonf.min():.3f}–{d.fwer_bonf.max():.3f}'))
    ax.axhline(0.05, color=S.INK, lw=0.8, ls='--', zorder=1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['alleles', 'residues'], fontsize=6.5)
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(0, 0.07)
    ax.set_ylabel('FWER under max-T')
    S.despine(ax, grid_axis='y')
    S.caption_block(
        fig, plot_h=plot_h,
        title=('Family-wise error rate under max-T: LD-adjusted (GEC) '
               f'{t.fwer_gec.min():.3f}–{t.fwer_gec.max():.3f}; conventional Bonferroni '
               f"{t[t['class'] == 'allele'].fwer_bonf.min():.3f}–{t[t['class'] == 'allele'].fwer_bonf.max():.3f} "
               f"(alleles) and {t[t['class'] == 'residue'].fwer_bonf.min():.3f}–"
               f"{t[t['class'] == 'residue'].fwer_bonf.max():.3f} (residues)."),
        panels=['Alleles: tests, the Me implied by the exact max-T threshold, and GEC\'s Me, per case set and cohort.',
                'Amino-acid residues, as a.',
                'Family-wise error rate each threshold actually has under the max-T simulation, all 24 families.'],
        notes=('A mirror residue pair counts as one test. Max-T: 400,000 multivariate-normal draws on the '
               'covariate-adjusted genotype correlation. README.md'),
        top_pad=0.62, wspace=0.62, left='auto', right=0.985, margin_axes=[axes[0]])
    for i, ax in enumerate(axes):
        S.panel_tag(ax, LETTERS[i])
    S.legend_above(axes[0], [Line2D([], [], ls='none', markersize=4, label=lab[k], **style[k])
                             for k in ('tests', 'me_maxt', 'me_gec')], ncol=3, fontsize=6.5, full_width=True,
                   y=1.07)
    save(fig, a.out_dir, 'thresholds', values)


# ═══════════════════════════════════════════════════════════════ 07.models
def fam_models(a, D):
    plot_h = 4.8
    fig, axes = plt.subplots(2, len(V.COHORTS), figsize=(S.COL_DOUBLE, plot_h + 2), squeeze=False)
    values, panels = [], []
    disagree = []
    for k, coh in enumerate(V.COHORTS):
        for r, what in enumerate(('p', 'beta')):
            ax = axes[r][k]
            lo_all, hi_all = [], []
            for cs in V.CASE_SETS:
                d = D.assoc[(D.assoc.case_set == cs) & (D.assoc.cohort == coh) & (D.assoc.test_id == D.assoc.marker)]
                w = d.pivot_table(index=['class', 'marker'], columns='model', values=what).dropna()
                x, y = (nlp(w['fixed']), nlp(w['random'])) if what == 'p' else (w['fixed'], w['random'])
                ax.scatter(x, y, s=4, color=V.CASE_SET_COLOR[cs], lw=0, alpha=0.7, zorder=3, rasterized=True)
                lo_all += [np.min(x), np.min(y)]
                hi_all += [np.max(x), np.max(y)]
                if what == 'p':
                    rr = float(np.corrcoef(x, y)[0, 1])
                    values.append((f'{V.COHORT_SHORT[coh]} {V.CASE_SET_SHORT[cs]}: r(-log10 P)', round(rr, 4)))
                    g = d.pivot_table(index='marker', columns='model', values='sig_gec', aggfunc='first')
                    if {'fixed', 'random'} <= set(g.columns):
                        for m in g.index[g['fixed'].astype(bool) != g['random'].astype(bool)]:
                            disagree.append((coh, cs, m))
            lo, hi = min(lo_all), max(hi_all)
            pad = (hi - lo) * 0.06
            lo, hi = (0, hi + pad) if what == 'p' else (lo - pad, hi + pad)
            ax.plot([lo, hi], [lo, hi], color=S.REFERENCE, lw=0.8, zorder=1)
            ax.set_xlim(lo, hi)
            ax.set_ylim(lo, hi)
            ax.set_box_aspect(1)
            S.despine(ax)
            labt = S.NEGLOG10P if what == 'p' else r'$\beta$ (log OR)'
            ax.set_xlabel(f'fixed {labt}')
            if k == 0:
                ax.set_ylabel(f'random {labt}')
    for r, what in enumerate(('p', 'beta')):
        for coh in V.COHORTS:
            panels.append(f'{V.COHORT_LABEL[coh]}, ' + ('−log10 P' if what == 'p' else 'β') + ', fixed (x) vs random (y).')
    values.append(('tests significant at GEC in one model only (cohort, case set, marker)',
                   '; '.join(f'{V.COHORT_SHORT[c]} {V.CASE_SET_SHORT[s]} {m}' for c, s, m in disagree) or 'none'))
    S.caption_block(
        fig, plot_h=plot_h,
        title=(f"Fixed against random model: r(−log10 P) {min(v for k, v in values if 'r(-log10 P)' in k):.3f}–"
               f"{max(v for k, v in values if 'r(-log10 P)' in k):.3f} over the case sets and cohorts; "
               f"{len(disagree)} test-cohort pair{'s' if len(disagree) != 1 else ''} significant at GEC in one "
               'model only.'),
        panels=panels,
        notes=('One point per test, all four case sets (colour). Grey line: identity. The fixed model drops '
               'relatives and analyses complete cases; the random model mean-imputes. README.md'),
        letters=LETTERS, top_pad=0.45, hspace=0.45, wspace=0.28, left='auto', right=0.985,
        margin_axes=[r[0] for r in axes])
    for r in range(2):
        for k, coh in enumerate(V.COHORTS):
            S.panel_tag(axes[r][k], LETTERS[r * 3 + k], title=V.COHORT_LABEL[coh] if r == 0 else None, tsize=8)
    S.legend_above(axes[0][1], [Line2D([], [], ls='none', marker='o', color=V.CASE_SET_COLOR[c], markersize=4,
                                       label=V.CASE_SET_SHORT[c]) for c in V.CASE_SETS],
                   ncol=4, fontsize=6.5, full_width=True, y=1.16)
    save(fig, a.out_dir, 'models', values)


# ═══════════════════════════════════════════════════════════════ 08.qq
def qq_band(n):
    i = np.arange(1, n + 1)
    return (-np.log10((i - 0.5) / n), -np.log10(stats.beta.ppf(0.975, i, n - i + 1)),
            -np.log10(stats.beta.ppf(0.025, i, n - i + 1)))


def fam_qq(a, D):
    plot_h = 7.6
    fig, axes = plt.subplots(len(V.CASE_SETS), len(V.CLASSES), figsize=(S.COL_DOUBLE * 0.72, plot_h + 2),
                             squeeze=False)
    values, panels = [], []
    for r, cs in enumerate(V.CASE_SETS):
        for k, cls in enumerate(V.CLASSES):
            ax = axes[r][k]
            lim, n_max, txt = 1.0, 0, []
            for coh in V.COHORTS:
                d = D.primary(cs, cls, coh)
                p = d.loc[d.test_id == d.marker, 'p'].dropna()
                n = len(p)
                n_max = max(n_max, n)
                obs = np.sort(nlp(p))[::-1]
                exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
                ax.plot(exp, obs, lw=1.0, color=V.COHORT_COLOR[coh], zorder=3)
                lim = max(lim, exp.max(), obs.max())
                lam = float(np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1))
                txt.append(f'{COH_ABBR[coh]} {lam:.2f}')
                values.append((f'{V.CASE_SET_SHORT[cs]} {cls} {V.COHORT_SHORT[coh]}: tests, lambda_GC',
                               f'{n}, {lam:.3f}'))
            e, lo, hi = qq_band(n_max)
            ax.fill_between(e, lo, hi, color=S.NEUTRAL, alpha=0.55, lw=0, zorder=1)
            lim *= 1.06
            ax.plot([0, lim], [0, lim], lw=0.8, color=S.REFERENCE, zorder=2)
            ax.set_xlim(0, lim)
            ax.set_ylim(0, lim)
            ax.set_box_aspect(1)
            ax.text(0.04, 0.96, f'{S.LAMBDA_GC}\n' + '\n'.join(txt), transform=ax.transAxes, ha='left', va='top',
                    fontsize=6, color=S.INK_SOFT, linespacing=1.25)
            if r == len(V.CASE_SETS) - 1:
                ax.set_xlabel(f'expected {S.NEGLOG10P}')
            if k == 0:
                ax.set_ylabel(f'{V.CASE_SET_SHORT[cs]}\nobserved {S.NEGLOG10P}')
            S.despine(ax)
            panels.append(f'{V.CASE_SET_SHORT[cs]}, {V.CLASS_PLURAL[cls]}.')
    S.caption_block(
        fig, plot_h=plot_h,
        title='P values against the uniform null, per case set and class, three cohorts overlaid.',
        panels=panels,
        notes=('Random model, one curve per cohort, one point per test. Grey: pointwise 95 % band. λGC describes '
               'one strongly correlated region and is never applied. README.md'),
        letters=LETTERS, top_pad=0.45, hspace=0.40, wspace=0.30, left='auto', right=0.98,
        margin_axes=[r[0] for r in axes])
    for r in range(len(V.CASE_SETS)):
        for k, cls in enumerate(V.CLASSES):
            S.panel_tag(axes[r][k], LETTERS[r * 2 + k], title={'allele': 'alleles', 'residue': 'residues'}[cls]
                        if r == 0 else None, tsize=8)
    S.legend_above(axes[0][0], [Line2D([], [], color=V.COHORT_COLOR[c], lw=1.4, label=V.COHORT_SHORT[c])
                                for c in V.COHORTS], ncol=3, fontsize=6.5, full_width=True, y=1.2)
    save(fig, a.out_dir, 'qq', values)


def main():
    a = parse_args()
    S.setup_style('paper')
    D = Data(a)
    globals()[f'fam_{a.family}'](a, D)


if __name__ == '__main__':
    main()
