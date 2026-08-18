#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The three additive scans stacked on ONE genomic axis — the figure
#           that answers "what does the genome look like in each of the three
#           nested sample sets, side by side".
#             one row per cohort, in the configured order
#           Each row is one cohort: Manhattan left, its QQ right. The row header
#           carries the composition (cases, controls, N_eff) and the QQ carries
#           the calibration (lambda_GC) — the same division of labour, and the
#           same 2.5 : 1 grid, as the scan figure's top row.
#
#           THE MANHATTAN IS THE SCAN FIGURE'S PANEL (a), 4.08 x 1.77 in. Matching
#           that box by construction rather than by arithmetic is what makes a peak
#           the same SHAPE on both figures. A full-width row would be 6.50 in
#           across and squash the same data to aspect 4.6 against 2.3.
#
#           ADDITIVE ONLY, like cohort_compare.png. At double-column width a
#           3 x 3 cohort-by-model grid leaves each Manhattan ~2.2 in across, at
#           which 22 chromosomes are an unreadable smear; the per-model scans
#           are on their own figures.
#
#           FOUR THINGS ARE SHARED ACROSS THE ROWS, and each is a correctness
#           requirement rather than a preference:
#             - one cumulative offset map, so a column is the same x in all
#               three rows (per-cohort offsets would put one chromosome at three
#               different x, which is the comparison the figure exists to make);
#             - one y-limit, so a taller peak is a stronger peak and not a
#               rescaled axis;
#             - one data height per row, because with a shared y-limit unequal
#               heights re-introduce exactly the distortion the shared limit
#               removes;
#             - one threshold pair, drawn on every row and labelled once.
#
#           The genome-wide loci are named ONCE, above the top row, over the
#           UNION across cohorts, with a dotted guide dropping through every
#           panel. Naming them per row would repeat any locus that several
#           cohorts called and read as several unrelated bands. Each row still
#           marks its OWN genome-wide leads with a diamond, so a guide with a
#           diamond in one row and none in another is a locus that only some
#           cohorts called.
# Component: assoc_plink2
# Used by : assoc_plink2.nf  process COMPARE_MANHATTAN
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MultipleLocator

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_style as S
import figure_doc
# load(), chrom_key(), read_peaks() and label_for() come from the single-scan
# figure. Reusing load() in particular is not merely DRY: it applies the same
# ERRCODE / degenerate-fit filter as PLOT_SCAN, which is what guarantees this
# figure plots the IDENTICAL variant set as the per-cohort scan figures.
import plot_manhattan_qq as MQ

PLOT_H = 7.32         # inches of croppable plot block (titles + axes + x-labels)
# Sized so each Manhattan lands on 4.08 x 1.77 in — the exact box scan.<model>.png
# gives its panel (a), hence the exact same aspect and the same apparent peak
# shape on both figures.
# The TOP row is given a larger slot because it is the only one that gives up a
# label strip, and equalise_row_heights() then takes the strip out of every row.
# Sized so that row (a) minus its strip lands on the same height as (b) and (c);
# if a future label set needs a taller strip the equaliser still guarantees equal
# heights, it just leaves a little slack above (b) and (c).
ROW0_EXTRA = 1.15
CHROM_GAP = 2.0e7     # blank genome between chromosomes, as on the scan figure
GW = 'genome_wide'
NUMWORD = {1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five'}


def parse_args():
    p = argparse.ArgumentParser(description='Three cohorts\' additive scans on one genomic axis.')
    p.add_argument('--glm', action='append', required=True, metavar='COHORT=PATH')
    p.add_argument('--scan-qc', action='append', required=True, metavar='COHORT=PATH')
    p.add_argument('--model-peaks', action='append', default=[], metavar='COHORT=PATH')
    p.add_argument('--cohort-order', default='',
                   help='comma-separated cohort order; empty means the order the per-cohort arguments were given in.')
    p.add_argument('--model', default='additive')
    p.add_argument('--alpha', type=float, default=S.GW_ALPHA)
    p.add_argument('--suggestive', type=float, default=S.SUGGEST_ALPHA)
    p.add_argument('--anno-style', default='auto', choices=sorted(S.ANNO_STYLES))
    p.add_argument('--repel-force', type=float, default=S.REPEL_FORCE)
    p.add_argument('--pc-label', default='',
                   help='Space the PCs were computed in, named in the model formula.')
    p.add_argument('--covar-label', default='SEX + PCs',
                   help='Human-readable covariate set, printed in the caption and the sidecar. '
                        'Configuration, not a property of the method — pass the run\'s own.')
    p.add_argument('--n-pcs', type=int, default=0,
                   help='Number of PCs in the covariate set; 0 leaves the formula\'s upper limit '
                        'as a generic K.')
    p.add_argument('--out-png', required=True)
    return p.parse_args()




def spec(items):
    """['cohort=path', …] -> {cohort: path}, dropping empty or missing files."""
    out = {}
    for s in items or []:
        c, _, path = s.partition('=')
        if path and Path(path).exists() and Path(path).stat().st_size:
            out[c] = path
    return out


def load_reduced(path):
    """One cohort -> ({chrom: (POS, -log10 P)}, P), releasing the DataFrame.

    Reduced immediately, one cohort at a time. The per-chromosome float32 heights
    are what the Manhattan draws; the raw P vector is kept because the QQ panel is
    MQ.qq() unchanged, and feeding it P recovered from a float32 log would give
    this figure a visibly different tail from the scan figure it must agree with.
    """
    df, n_all, n_err = MQ.load(path)
    out = {str(c): (g['POS'].to_numpy(dtype=np.int64),
                    (-np.log10(g['P'].to_numpy())).astype(np.float32))
           for c, g in df.groupby('CHROM', sort=False)}
    return out, df['P'].to_numpy(), len(df), n_all, n_err


def genome_axis(per_cohort):
    """ONE cumulative offset map for every row: {chrom: offset}, ticks, names, xlim.

    The single-scan figure builds its offsets from its own dataframe, which is
    correct there because the axis serves one dataset. Doing that here would put
    a chromosome at a different x in each row and destroy the only comparison
    this figure exists to make, so the map is built from the UNION of per-chromosome
    extents — a column is then the same pixel in all three panels whatever each
    cohort's own call set happens to span.
    """
    chroms = sorted({c for d in per_cohort for c in d}, key=MQ.chrom_key)
    off, ticks, names, cur = {}, [], [], 0.0
    for c in chroms:
        lo = min(d[c][0].min() for d in per_cohort if c in d)
        hi = max(d[c][0].max() for d in per_cohort if c in d)
        off[c] = cur - lo
        ticks.append(cur + (hi - lo) / 2.0)
        names.append(c.replace('chr', ''))
        cur += (hi - lo) + CHROM_GAP
    return off, ticks, names, (-1e7, cur - CHROM_GAP + 1e7)


def union_loci(peaks, order):
    """The genome-wide leads of ANY cohort, once each, smallest P winning.

    Deduplicated by label rather than by variant: when two cohorts call the same
    locus their leads need not be the same variant, but it is one locus and
    belongs on the figure once.
    """
    best = {}
    for c in order:
        pk = peaks.get(c)
        if pk is None or not len(pk):
            continue
        for _, r in pk[pk.tier == GW].iterrows():
            lab = MQ.label_for(r)
            if lab not in best or float(r['P']) < best[lab]['P']:
                best[lab] = dict(label=lab, chrom=str(r['chrom']), pos=int(r['pos']),
                                 P=float(r['P']), cohort=c)
    return sorted(best.values(), key=lambda d: (MQ.chrom_key(d['chrom']), d['pos']))


def guide_example(loci, data, peaks, order):
    """A (label, cohort, P) case where a guide carries no diamond, FROM THE DATA.

    The sidecar used to name a locus and a P-value as its worked example, which
    made every run's documentation describe one particular dataset. The example is
    the same idea whatever the data, so it is looked up instead: the first union
    locus that some cohort did not call genome-wide, with that cohort's own P at
    the same position. Returns None when every cohort called everything.
    """
    for d in loci:
        for c in order:
            pk = peaks.get(c)
            if pk is not None and len(pk) and (
                    (pk.tier == GW)
                    & (pk.chrom.astype(str) == d['chrom'])
                    & (pk.pos.astype('int64') == d['pos'])).any():
                continue
            arr = data.get(c, {}).get(d['chrom'])
            if arr is None:
                continue
            pos, y = arr
            i = np.flatnonzero(pos == d['pos'])
            if len(i):
                return d['label'], c, float(10.0 ** -float(y[i[0]]))
    return None


def draw_row(ax, data, off, guides, gw, ylim, args, *, thresh_labels=False):
    """One cohort's Manhattan on the shared axis."""
    for i, c in enumerate(sorted(data, key=MQ.chrom_key)):
        pos, y = data[c]
        ax.scatter(off[c] + pos, y, s=3.2, c=S.CHROM_BANDS[i % 2],
                   linewidths=0, rasterized=True, zorder=2)
    # Guides first and underneath: they carry the eye from a called peak in one
    # row to the same column in the others, which is the vertical comparison.
    for x in guides:
        ax.axvline(x, color=S.NEUTRAL_D, lw=0.7, ls=':', zorder=1)
    ax.axhline(-np.log10(args.alpha), color=S.ACCENT, lw=1.1, zorder=3)
    ax.axhline(-np.log10(args.suggestive), color=S.REFERENCE, lw=0.9, ls='--', zorder=3)
    # This row's OWN genome-wide leads. A locus that is genome-wide elsewhere but
    # not here gets a guide and no diamond, which is the distinction the figure
    # is for.
    if len(gw):
        st = S.TIER_STYLE[GW]
        ax.scatter([off.get(str(r.chrom), 0) + int(r.pos) for _, r in gw.iterrows()],
                   -np.log10(gw['P'].astype(float)), s=st['size'], marker=st['marker'],
                   c=st['color'], edgecolors='white', linewidths=0.6, zorder=6)
    if thresh_labels:
        # Named once, on the top row. Three copies of the same two numbers is the
        # redundant text this layout drops. Placed in AXES fraction on x, not in
        # data coordinates: the shared x-limit is set after every row is drawn, so
        # a data-coordinate anchor lands wherever the autoscale happened to be and
        # is then clipped by the final limit.
        box = dict(boxstyle='square,pad=0.15', facecolor='white', edgecolor='none', alpha=0.85)
        for a, col in ((args.alpha, S.ACCENT), (args.suggestive, S.REFERENCE)):
            ax.text(0.998, -np.log10(a), S.p_tex(a).replace('1.00', '1').replace('5.00', '5'),
                    transform=ax.get_yaxis_transform(), color=col,
                    fontsize=plt.rcParams['legend.fontsize'] - 1, ha='right', va='bottom',
                    zorder=7, bbox=box)
    ax.set_ylim(*ylim)
    # Fixed locator, not the automatic one. The rows share a y-limit but their
    # heights can differ by hundredths of an inch, and that is enough to flip the
    # auto-locator from a step of 2 to a step of 1 on one row — three rows meant
    # to be read against each other must at least carry the same ticks.
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.set_ylabel(S.NEGLOG10P)
    S.despine(ax, grid_axis='y')


def header(ax, letter, cohort, q, pad):
    """Row header in the title slot: bold letter + cohort, then this row's numbers.

    Both titles live in the slot matplotlib already reserves above the axes, so
    neither can land on the data — the same rule that puts the gene names in a
    strip rather than inside the panel.
    """
    S.panel_tag(ax, f'{letter}   {cohort}', pad=pad)
    # Composition here, calibration in the QQ beside it — the same division of
    # labour as the scan figure, and it keeps lambda_GC from being printed twice
    # on one row.
    ax.set_title(f"{int(q['n_case']):,} cases $\\cdot$ {int(q['n_ctrl']):,} controls $\\cdot$ "
                 f"$N_{{\\mathrm{{eff}}}}$ = {float(q['n_eff']):,.0f}",
                 loc='right', fontsize=plt.rcParams['legend.fontsize'], color=S.INK, pad=pad)


def main():
    args = parse_args()
    S.setup_style('paper')
    glm_spec = spec(args.glm)
    want = [c for c in args.cohort_order.split(',') if c] or list(glm_spec)
    order = [c for c in want if c in glm_spec]
    glm, qcs, mpk = glm_spec, spec(args.scan_qc), spec(args.model_peaks)

    # One cohort at a time: load, reduce, release.
    data, pvals, stats = {}, {}, {}
    for c in order:
        data[c], pvals[c], n_kept, n_all, n_err = load_reduced(glm[c])
        stats[c] = dict(n_kept=n_kept, n_all=n_all, n_err=n_err)
    qc = {}
    for c, p in qcs.items():
        d = pd.read_csv(p, sep='\t')
        qc[c] = d[d['model'] == args.model].iloc[0]
    peaks = {c: MQ.read_peaks(p, args.model) for c, p in mpk.items()}

    SHORT = S.shorten(order)
    off, ticks, names, xlim = genome_axis([data[c] for c in order])
    loci = union_loci(peaks, order)
    guides = [off.get(d['chrom'], 0) + d['pos'] for d in loci]
    # ONE y-limit for every row: a taller peak has to mean a stronger peak.
    top = max(float(y.max()) for c in order for _, y in data[c].values())
    ylim = (0, max(8.4, top * 1.08))

    fig = plt.figure(figsize=(S.COL_DOUBLE, PLOT_H))
    # THE SAME GRID AS THE SCAN FIGURE'S TOP ROW, stacked once per cohort:
    # width_ratios [2.5, 1.0] with the same margins and wspace put each Manhattan
    # at 4.08 x 1.77 in — identical to scan.<model>.png panel (a). Matching the
    # aspect by construction rather than by arithmetic is what keeps a peak the
    # same SHAPE here as on the figure it is also drawn on; a full-width row is
    # 6.50 in across and squashes the same data to aspect 4.6 against 2.3.
    gs = fig.add_gridspec(len(order), 2,
                          height_ratios=[ROW0_EXTRA] + [1.0] * (len(order) - 1),
                          width_ratios=[2.5, 1.0])
    axes = [fig.add_subplot(gs[0, 0])]
    axes += [fig.add_subplot(gs[i, 0], sharex=axes[0]) for i in range(1, len(order))]
    qqs = [fig.add_subplot(gs[i, 1]) for i in range(len(order))]

    for i, (ax, c) in enumerate(zip(axes, order)):
        pk = peaks.get(c, pd.DataFrame())
        draw_row(ax, data[c], off, guides,
                 pk[pk.tier == GW] if len(pk) else pd.DataFrame(), ylim, args,
                 thresh_labels=(i == 0))
        # lambda_GC from scan_qc.tsv, not recomputed here, so the number beside
        # the QQ is the same one the header would have printed and the same one
        # scan.<model>.png and cohort_compare.png report.
        MQ.qq(qqs[i], pvals[c], float(qc[c]['lambda_gc']))
    axes[0].set_xlim(*xlim)
    # The three QQs are the same quantity for different cohorts, so they get ONE
    # limit too — otherwise the identity line sits at a different place in each
    # and the curves cannot be compared by eye any more than the Manhattans could.
    qlim = max(q.get_xlim()[1] for q in qqs)
    for i, q in enumerate(qqs):
        q.set_xlim(0, qlim)
        q.set_ylim(0, qlim)
        if i != len(order) - 1:
            q.set_xlabel('')
    # Chromosome numbers on the BOTTOM row only: the three rows are the same
    # quantity for different cohorts, so three identical axes is redundant text
    # and costs ~0.5 in of data height across the stack.
    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)
    axes[-1].set_xticks(ticks)
    axes[-1].set_xticklabels(names, fontsize=plt.rcParams['xtick.labelsize'] - 2)
    axes[-1].set_xlabel('chromosome')

    n_gw = {c: int((peaks[c].tier == GW).sum()) if c in peaks and len(peaks[c]) else 0
            for c in order}
    named = ', '.join(d['label'] for d in loci) or 'none'
    GWT = r'$5\times10^{-8}$'

    def panel_line(c):
        """What this row shows that its header does not already say.

        Not '<cohort> — N loci reach 5e-8': the cohort name and the counts are in
        the row header two lines above. The useful sentence names WHICH loci
        cross, and for a row where none does, which column comes closest.
        """
        pk = peaks.get(c, pd.DataFrame())
        gw = pk[pk.tier == GW] if len(pk) else pd.DataFrame()
        if len(gw):
            # Named left-to-right, matching the label band and the caption title.
            # Ordering by P instead would list them differently in three places.
            names_here = [d['label'] for d in loci
                          if d['label'] in {MQ.label_for(r) for _, r in gw.iterrows()}]
            return (', '.join(names_here or [MQ.label_for(r) for _, r in gw.iterrows()])
                    + f' cross{"" if len(gw) > 1 else "es"} {GWT}.')
        rest = pk[pk.tier != GW] if len(pk) else pd.DataFrame()
        if len(rest):
            r = rest.nsmallest(1, 'P').iloc[0]
            return (f'No locus crosses {GWT}; {MQ.label_for(r)} is the tallest column, at '
                    f'$P$ = {S.p_tex(float(r["P"]))}.')
        return f'No locus crosses {GWT}.'

    fig.suptitle(f'{args.model.capitalize()} scan across {NUMWORD.get(len(order), len(order))} nested sample sets',
                 fontsize=plt.rcParams['axes.titlesize'] + 1, fontweight='bold', y=0.988)

    S.caption_block(
        fig,
        title=(f'{args.model.capitalize()} scan in {NUMWORD.get(len(order), len(order))} nested sample sets on one genomic axis: '
               + (f'{named} reach {GWT} in at least one set — ' if loci
                  else f'no locus reaches {GWT} — ')
               + ', '.join(str(n_gw[c]) for c in order)
               + f' loci in {", ".join(SHORT[c] for c in order)}.'),
        panels=[panel_line(c) for c in order],
        notes=(f'Each row: Manhattan left, its QQ against the uniform null right, with '
               f'{S.LAMBDA_GC} and the pointwise 95% concentration band. Additive model only; '
               f'{args.covar_label}. One shared genomic axis, one shared {S.NEGLOG10P} limit '
               f'and one threshold pair, so heights are comparable between rows; dotted guides mark '
               f'every locus genome-wide in any cohort, diamonds mark the leads that cohort itself '
               f'called. Full explanation: {Path(args.out_png).stem}.md'),
        plot_h=PLOT_H, top_pad=0.46, hspace=0.20, wspace=0.30, left='auto', right=0.988,
        margin_axes=axes)

    # AFTER caption_block, which resizes the figure: gene_labels measures the
    # rendered axes box. ONE band, on the top row, over the union — and anchored
    # at the ceiling so the arm is a stub and the dotted guide is the leader.
    if loci:
        # INK, not ACCENT. In this palette ACCENT means "genome-wide", and these
        # names are the UNION across cohorts sitting above a row that may have
        # called none of them — row (a) calls neither. Red is left to the diamonds,
        # which do mean "this cohort called it"; the band only names the columns.
        S.gene_labels(axes[0], guides, [ylim[1]] * len(loci), [d['label'] for d in loci],
                      anno_style=args.anno_style, repel_force=args.repel_force,
                      color=S.INK, weight='bold')
    # Only the top row gave up a strip; bring the others to the same data height,
    # because a shared y-limit on unequal heights is not a shared scale.
    S.equalise_row_heights(*axes)
    # NOT align_panel_tops here. MQ.qq sets box_aspect(1), so a QQ's height is
    # fixed by its column width, and align_panel_tops would compare the Manhattan
    # to that squared box and shrink the MANHATTAN onto it — backwards, since the
    # Manhattans are the panels whose heights must stay equal and whose aspect
    # must match the scan figure. It cost rows (b) and (c) 0.01 in each, enough to
    # move them to a different y-tick locator than row (a). Give each QQ its row's
    # vertical extent instead and let it square itself against the top.
    for ax, q in zip(axes, qqs):
        p, pq = ax.get_position(), q.get_position()
        q.set_position([pq.x0, p.y0, pq.width, p.height])
        q.set_anchor('N')
    pad = S.strip_pad(axes[0])
    for i, (ax, c) in enumerate(zip(axes, order)):
        header(ax, 'abcdefgh'[i], c, qc[c], pad if i == 0 else 7.0)
    S.thin_tick_labels(axes[-1], 'x')
    fig.savefig(args.out_png)
    plt.close(fig)

    _ex = guide_example(loci, data, peaks, order)
    _lams = [float(qc[c]['lambda_gc']) for c in order]
    _n_sug = sum(int((peaks[c].tier != GW).sum()) if c in peaks and len(peaks[c]) else 0
                 for c in order)
    figure_doc.write_doc(
        args.out_png,
        title=(f'{args.model.capitalize()} scan across {NUMWORD.get(len(order), len(order))} nested sample sets — '
               + ', '.join(SHORT[c] for c in order)),
        question=('Which associations are shared across the three nested sample sets, and which '
                  'appear only as the ancestry filter is relaxed?'),
        interpretation=(
            'The three sample sets are NESTED — narrow within intermediate within full — so they '
            'share most of their samples and agreement between rows is expected rather than '
            'informative. This figure is not replication and must not be read as such. What it does '
            'show is where a locus sits in every set at once: a dotted guide marks each locus that '
            'is genome-wide in any cohort, and a diamond marks it only in the cohorts that actually '
            'called it, so a column that rises across the rows is a signal gaining sample size '
            'rather than a new finding. All three rows are drawn on one genomic axis, one '
            'lambda-independent y-limit and one threshold pair, and on equal data heights, so peak '
            'heights are directly comparable between rows.'),
        panels=[('abcdefgh'[i],
                 f'{c} — Manhattan and QQ',
                 f'LEFT — every analysed variant of the {c} additive scan at its genomic position, '
                 f'y = -log10 P, on the shared axis; drawn without thinning and rasterised. The box '
                 f'is 4.08 x 1.77 in, the same box scan.{args.model}.png gives its Manhattan, so a '
                 f'peak has the same shape on both. RIGHT — the same P-values ranked against a '
                 f'uniform null, with the pointwise 95% concentration band from the order-statistic '
                 f'Beta(i, n-i+1) distribution and this cohort\'s lambda_GC. All three QQs carry one '
                 f'limit, for the same reason the Manhattans do. The header counts and lambda_GC come '
                 f'from this cohort\'s own 02.scan/scan_qc.tsv rather than being recomputed, so they '
                 f'agree with scan.{args.model}.png and cohort_compare.png exactly. '
                 + ('Chromosome numbers are on this row, shared by all three.'
                    if i == len(order) - 1 else
                    'Chromosome numbers are on the bottom row, which this row shares.'))
                for i, c in enumerate(order)],
        numbers=([('model', args.model), ('cohorts', ', '.join(order))]
                 + [(f'{c}: cases', int(qc[c]['n_case'])) for c in order]
                 + [(f'{c}: controls', int(qc[c]['n_ctrl'])) for c in order]
                 + [(f'{c}: N_eff', float(qc[c]['n_eff'])) for c in order]
                 + [(f'{c}: lambda_GC', float(qc[c]['lambda_gc'])) for c in order]
                 + [(f'{c}: variants analysed', stats[c]['n_kept']) for c in order]
                 + [(f'{c}: excluded on ERRCODE / degenerate fit', stats[c]['n_err']) for c in order]
                 + [(f'{c}: genome-wide peaks', n_gw[c]) for c in order]
                 + [('shared -log10 P ceiling', round(ylim[1], 3)),
                    ('loci named (union, genome-wide in any cohort)', named),
                    ('annotation style', args.anno_style),
                    ('repel force', float(args.repel_force))]),
        tables=[('Loci genome-wide in at least one cohort',
                 pd.DataFrame(loci)[['label', 'chrom', 'pos', 'P', 'cohort']]
                 .rename(columns={'cohort': 'smallest_P_in'}) if loci else None)],
        reading=[
            'Read DOWN a guide, not across a row. The guide is the same genomic position in all '
            'three panels, so it shows what one locus does as the ancestry filter is relaxed.',
            ('A diamond means that cohort called the locus genome-wide. A guide with no diamond '
             'means the locus is there but below the threshold in that cohort'
             + (f' — {_ex[0]} in {_ex[1]} is the case here, at P = {_ex[2]:.3g}.' if _ex else '.')),
            'Compare peak heights between rows directly: the y-limit and the axes height are shared, '
            'so a taller column is a smaller P and not a rescaled panel.',
            (f'Read the header numbers together with the row, and the QQ beside it. lambda_GC spans '
             f'{min(_lams):.3f}-{max(_lams):.3f} across the {len(order)} sets; if that span is narrow '
             f'and the QQ curves have the same shape, the calibration does not change as the filter '
             f'is relaxed, so a column that grows down the figure is growing with N and not with '
             f'inflation.'),
        ],
        limits=[
            'The three sample sets are nested, so this is not replication and no independent cohort '
            'is available to this study. A locus present in all three rows has been seen once.',
            (f'Only the {args.model} model is drawn. Peaks called under the other genetic models are '
             f'on their own scan figures and are absent here by construction.'),
            (f'Suggestive peaks are not marked. The dashed line shows the tier; naming them across '
             f'{len(order)} rows would put {_n_sug} labels on the figure.'),
            'Variants plink2 could not fit cleanly are absent from every row; their counts and error '
            'codes are in each cohort\'s 02.scan/scan_qc.tsv.',
        ],
        defs=['model', 'lambda_gc', 'gw_sig', 'neff'],
        model=S.glm_formula(args.pc_label, args.n_pcs) + '\n' + S.FORMULAS['lambda'],
        methods_ref='../../../docs/METHODS.md')
    print(f'[plot_cohort_manhattan] {len(order)} cohorts, {args.model}: '
          + ', '.join(f'{SHORT[c]} {n_gw[c]}gw' for c in order)
          + f'; union={named} -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in plot_cohort_manhattan: {e}', file=sys.stderr)
        sys.exit(1)
