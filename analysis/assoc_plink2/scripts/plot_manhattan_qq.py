#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The scan figure for one cohort x model — everything that describes
#           that one scan, in one figure:
#             (a) Manhattan over all analysed variants, genome-wide loci named
#             (b) QQ with the 95 % concentration band and lambda_GC
#             (c) its peak leads on the SAME genomic axis as (a), gene-labelled
#             (d) lead |effect| against lead MAF
#
#           2x2: the two GENOMIC panels are the left column and share an x-axis,
#           so a peak in (c) sits directly under the Manhattan column it came
#           from; the two summary panels are the right column. Stacking all four
#           left a whole half-width of the figure empty.
#
#           Peaks come from THIS MODEL. Additive is the primary analysis and the
#           only model that defines the downstream fan-out, but a dominant or
#           recessive figure drawn with additive peaks says nothing about the
#           scan it is drawn from.
#
#           Every point is drawn — no thinning. Dropping variants below
#           -log10 P = 2 would change the QQ's lower arm, the very region lambda
#           is read from. The scatter is rasterised instead, so the file stays a
#           few MB while axes, ticks and text remain vector.
# Component: assoc_plink2
# Used by : assoc_plink2.nf  process PLOT_SCAN
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_style as S
import figure_doc

PLOT_H = 6.0          # inches of croppable plot block (titles + axes + x-labels)
MODEL_LABEL = {'additive': 'additive (ADD)', 'dominant': 'dominant (DOM)',
               'recessive': 'recessive (REC)'}
GW, SUG = 'genome_wide', 'suggestive'


def parse_args():
    p = argparse.ArgumentParser(description='Scan figure for one cohort x model.')
    p.add_argument('--glm', required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--model', required=True)
    p.add_argument('--model-peaks', help='model_peaks_annotation.tsv; rows for this model are drawn.')
    p.add_argument('--alpha', type=float, default=S.GW_ALPHA)
    p.add_argument('--suggestive', type=float, default=S.SUGGEST_ALPHA)
    p.add_argument('--label-suggestive', type=int, default=10,
                   help='How many suggestive leads to label, smallest P first.')
    p.add_argument('--anno-style', default='expand', choices=sorted(S.ANNO_STYLES),
                   help="gwaslab annotation style; 'expand' is this component's default.")
    p.add_argument('--repel-force', type=float, default=S.REPEL_FORCE,
                   help='Minimum label separation as a fraction of the plotted span.')
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


def load(path):
    """Tab-delimited -> C parser. Only the columns the figure needs are read."""
    head = pd.read_csv(path, sep='\t', nrows=0)
    cols = [c for c in ['#CHROM', 'POS', 'ID', 'OR', 'LOG(OR)_SE', 'P', 'ERRCODE']
            if c in head.columns]
    df = pd.read_csv(path, sep='\t', usecols=cols, dtype={'#CHROM': str, 'ID': str, 'ERRCODE': str})
    df = df.rename(columns={'#CHROM': 'CHROM', 'LOG(OR)_SE': 'SE'})
    n_all = len(df)
    for c in ('P', 'OR', 'SE'):
        if c in df:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    # Same usability rule as call_peaks.py: ERRCODE alone is not sufficient,
    # because plink2 leaves it '.' on degenerate fits that return OR = inf and
    # P underflowed to 0. Plotting those would put spurious points at the very
    # top of the Manhattan.
    ok = df['ERRCODE'].fillna('.').eq('.') if 'ERRCODE' in df else pd.Series(True, index=df.index)
    ok &= df['P'].notna() & (df['P'] > 0) & (df['P'] <= 1)
    if 'OR' in df:
        ok &= np.isfinite(df['OR']) & (df['OR'] > 0)
    if 'SE' in df:
        ok &= np.isfinite(df['SE']) & (df['SE'] > 0)
    df = df[ok]
    return df, n_all, n_all - len(df)


def chrom_key(c):
    c = str(c).replace('chr', '')
    return {'X': 23, 'Y': 24, 'MT': 25, 'M': 25}.get(c, int(c) if c.isdigit() else 99)


def read_peaks(path, model):
    """This model's peaks, both tiers, from model_peaks_annotation.tsv."""
    if not path or not Path(path).exists() or not Path(path).stat().st_size:
        return pd.DataFrame()
    d = pd.read_csv(path, sep='\t', dtype={'variant_id': str, 'model': str, 'chrom': str})
    return d[d['model'] == model].copy() if 'model' in d else d


def label_for(r):
    """Gene symbol, falling back to rsID and then the variant ID."""
    g, rs = str(r.get('Gene', '.')), str(r.get('rsID', '.'))
    if g not in ('.', '', 'nan'):
        return g
    if rs not in ('.', '', 'nan'):
        return rs
    return str(r['variant_id'])


def manhattan(ax, df, args):
    """Draw the Manhattan; return {chrom: cumulative offset}, tick positions, span."""
    df = df.assign(_k=df['CHROM'].map(chrom_key)).sort_values(['_k', 'POS'])
    offset, off, ticks, names, lo_hi = 0, {}, [], [], []
    for i, (_k, g) in enumerate(df.groupby('_k', sort=True)):
        c = str(g['CHROM'].iloc[0])
        off[c] = offset
        xs = offset + g['POS'].to_numpy()
        ticks.append((xs[0] + xs[-1]) / 2)
        names.append(c.replace('chr', ''))
        ax.scatter(xs, -np.log10(g['P'].to_numpy()), s=3.2,
                   c=S.CHROM_BANDS[i % 2], linewidths=0, rasterized=True, zorder=2)
        lo_hi.append((xs[0], xs[-1]))
        offset = xs[-1] + 2.0e7
    x0, x1 = lo_hi[0][0] - 1e7, lo_hi[-1][1] + 1e7

    ax.axhline(-np.log10(args.alpha), color=S.ACCENT, lw=1.2, zorder=3)
    ax.axhline(-np.log10(args.suggestive), color=S.REFERENCE, lw=0.9, ls='--', zorder=3)
    # Both labels sit at the right edge, over the last chromosome's column, so
    # they carry an opaque backing — a recessive scan's points reach that height.
    box = dict(boxstyle='square,pad=0.15', facecolor='white', edgecolor='none', alpha=0.85)
    for a, col in ((args.alpha, S.ACCENT), (args.suggestive, S.REFERENCE)):
        ax.text(x1, -np.log10(a), S.p_tex(a).replace('1.00', '1').replace('5.00', '5'),
                color=col, fontsize=plt.rcParams['legend.fontsize'] - 1,
                ha='right', va='bottom', zorder=7, bbox=box)
    ax.set_xlim(x0, x1)
    ax.set_ylabel(S.NEGLOG10P)
    # Natural range. The gene labels live in a strip ABOVE the axes (see
    # S.gene_labels), so nothing has to be reserved inside the data area.
    ax.set_ylim(0, max(8.4, -np.log10(df['P'].min()) * 1.12))
    S.despine(ax, grid_axis='y')
    return off, ticks, names


def peak_panel(ax, pk, off, ticks, names, args):
    """This scan's peak leads on the Manhattan's own x-axis, gene-labelled."""
    for tier in (SUG, GW):
        sel = pk[pk.tier == tier]
        if not len(sel):
            continue
        st = S.TIER_STYLE[tier]
        x = [off.get(str(r.chrom), 0) + int(r.pos) for _, r in sel.iterrows()]
        ax.scatter(x, -np.log10(sel['P'].astype(float)), s=st['size'], marker=st['marker'],
                   c=st['color'], edgecolors='white', linewidths=0.6,
                   zorder=6 if tier == GW else 5)
    ax.axhline(-np.log10(args.alpha), color=S.ACCENT, lw=1.1, zorder=3)
    ax.axhline(-np.log10(args.suggestive), color=S.REFERENCE, lw=0.9, ls='--', zorder=3)

    ax.set_xticks(ticks)
    ax.set_xticklabels(names, fontsize=plt.rcParams['xtick.labelsize'] - 2)
    ax.set_xlabel('chromosome')
    ax.set_ylabel(f'peak lead {S.NEGLOG10P}')
    # Natural range: the labels are in a strip above the axes, so the ceiling is
    # set by the DATA. It used to be 1.95x the largest -log10 P to hold two label
    # bands, which left the peaks in the bottom third of the panel.
    lo = -np.log10(args.suggestive) - 0.35
    hi = max(-np.log10(pk['P'].astype(float).min()) * 1.12, 8.0) if len(pk) else 8.0
    ax.set_ylim(lo, hi)
    S.despine(ax, grid_axis='y')


def _draw_band(ax, sel, off, args):
    """One label band for a set of peaks, tier carried by colour and weight.

    ONE band per panel, not one per tier. Two bands cost two strips — on the peak
    panel that was 0.98 in of text above a 1.00 in axes — and put two text
    orientations in the same panel. Merged, the repel also sees every label at
    once, so the genome-wide and suggestive names are spaced against each other
    instead of in two independent passes that happen to land in adjacent strips.
    """
    if not len(sel):
        return
    sel = sel.assign(_x=[off.get(str(r.chrom), 0) + int(r.pos)
                         for _, r in sel.iterrows()]).sort_values('_x')
    is_gw = (sel.tier == GW).to_numpy()
    base = plt.rcParams['legend.fontsize']
    S.gene_labels(
        ax, sel['_x'].to_numpy(), -np.log10(sel['P'].astype(float).to_numpy()),
        [label_for(r) for _, r in sel.iterrows()],
        anno_style=args.anno_style, repel_force=args.repel_force,
        color=[S.ACCENT if g else S.REFERENCE for g in is_gw],
        weight=['bold' if g else 'normal' for g in is_gw],
        fontsize=[base if g else base - 1.0 for g in is_gw])


def label_panels(ax_man, ax_pk, pk, off, args):
    """Gene labels on both genomic panels. MUST run AFTER caption_block.

    S.gene_labels measures the rendered axes box to size the strip and the leader
    arms, and caption_block resizes the figure and re-runs subplots_adjust, so
    labels placed before it land at the wrong height.

    The Manhattan carries the genome-wide loci only — it is the panel a reader
    looks at first and it used to be the one panel with no gene name on it. The
    peak panel carries those AND the strongest suggestive ones, in one band.
    """
    gw = pk[pk.tier == GW]
    sug = pk[pk.tier == SUG].nsmallest(max(0, args.label_suggestive), 'P')
    _draw_band(ax_man, gw, off, args)
    _draw_band(ax_pk, pd.concat([gw, sug]) if len(gw) or len(sug) else gw, off, args)


def qq(ax, p, lam):
    n = len(p)
    obs = -np.log10(np.sort(p))
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    # Beta(i, n-i+1) pointwise 95 % band — the concentration expected under the
    # null, so departures are judged against a stated reference, not by eye.
    i = np.arange(1, n + 1)
    lo = -np.log10(stats.beta.ppf(0.975, i, n - i + 1))
    hi = -np.log10(stats.beta.ppf(0.025, i, n - i + 1))
    ax.fill_between(exp, lo, hi, color=S.NEUTRAL, linewidth=0, zorder=1)
    m = max(exp.max(), obs.max())
    ax.plot([0, m], [0, m], color=S.ACCENT, lw=1.1, zorder=2)
    ax.scatter(exp, obs, s=4.0, c=S.DATA_DARK, linewidths=0, rasterized=True, zorder=3)
    ax.set_xlabel(f'expected {S.NEGLOG10P}')
    ax.set_ylabel(f'observed {S.NEGLOG10P}')
    # Square box and one limit for both axes. Observed and expected are the same
    # quantity in the same units, so anything else draws the identity line at an
    # angle that is not 45 degrees and misstates the departure.
    lim = m * 1.05
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_box_aspect(1)
    ax.text(0.04, 0.95, f'{S.LAMBDA_GC} = {lam:.3f}', transform=ax.transAxes,
            ha='left', va='top', fontsize=plt.rcParams['legend.fontsize'] + 1, fontweight='bold')
    S.despine(ax)


def effect_panel(ax, pk):
    """Lead |effect| against lead MAF — where this scan's peaks sit in the
    frequency/effect plane."""
    if not len(pk) or 'OR' not in pk or 'MAF' not in pk:
        ax.text(0.5, 0.5, 'no peaks in this scan', transform=ax.transAxes,
                ha='center', va='center', color=S.REFERENCE)
        ax.set_xticks([]); ax.set_yticks([])
        return
    d = pk.dropna(subset=['OR', 'MAF'])
    for tier in (SUG, GW):
        sel = d[d.tier == tier]
        if not len(sel):
            continue
        st = S.TIER_STYLE[tier]
        # |effect|, so a protective and a risk allele of equal strength plot at
        # the same height; direction is in the tables.
        eff = np.exp(np.abs(np.log(sel['OR'].astype(float))))
        ax.scatter(sel['MAF'].astype(float), eff, s=st['size'], marker=st['marker'],
                   c=st['color'], edgecolors='white', linewidths=0.6,
                   zorder=6 if tier == GW else 5)
    ax.set_xlabel(S.MAF_SYM)
    ax.set_ylabel(f'|effect| as {S.OR_SYM}')
    ax.set_xlim(0, 0.505)
    top = float(np.exp(np.abs(np.log(d['OR'].astype(float)))).max()) if len(d) else 2.0
    ax.set_ylim(1.0, min(max(top * 1.12, 2.0), 8.0))
    ax.set_yscale('log')
    S.or_log_axis(ax, axis='y')
    S.despine(ax)


def main():
    args = parse_args()
    S.setup_style('paper')
    df, n_all, n_err = load(args.glm)
    lam = float(stats.chi2.isf(np.median(df['P']), df=1) / 0.4549364231195729)
    pk = read_peaks(args.model_peaks, args.model)
    mlab = MODEL_LABEL.get(args.model, args.model)
    n_gw = int((pk.tier == GW).sum()) if len(pk) else 0
    n_sug = int((pk.tier == SUG).sum()) if len(pk) else 0

    fig = plt.figure(figsize=(S.COL_DOUBLE, PLOT_H))
    # 2x2. LEFT COLUMN = the two genomic panels, sharing an x-axis so a peak in
    # (c) sits under the Manhattan column it came from. RIGHT COLUMN = the two
    # summary panels. Stacking all four left half the figure width empty.
    # Row 2 gets the larger share because its label strip is taken out of it:
    # after the strip, (c)'s axes should match (a)'s rather than end up half its
    # height, which is what [1.0, 1.18] produced.
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.5], width_ratios=[2.5, 1.0])
    ax_man = fig.add_subplot(gs[0, 0])
    ax_qq = fig.add_subplot(gs[0, 1])
    ax_pk = fig.add_subplot(gs[1, 0], sharex=ax_man)
    ax_eff = fig.add_subplot(gs[1, 1])

    off, ticks, names = manhattan(ax_man, df, args)
    # sharex keeps (a) and (c) aligned to the pixel; both print their OWN
    # chromosome numbers. Only (c) carries the axis TITLE — two identical
    # 'chromosome' labels would be exactly the redundant text this layout drops.
    ax_man.set_xticks(ticks)
    ax_man.set_xticklabels(names, fontsize=plt.rcParams['xtick.labelsize'] - 2)
    peak_panel(ax_pk, pk, off, ticks, names, args)
    qq(ax_qq, df['P'].to_numpy(), lam)
    effect_panel(ax_eff, pk)

    fig.suptitle(f'{args.cohort} — {mlab}',
                 fontsize=plt.rcParams['axes.titlesize'] + 1, fontweight='bold', y=0.988)

    n_hits = int((df['P'] < args.alpha).sum())
    gw_names = ', '.join(label_for(r) for _, r in pk[pk.tier == GW].iterrows()) if n_gw else 'none'
    n_lab_sug = min(args.label_suggestive, n_sug)
    calib = ''
    if lam < 0.9:
        calib = (f' This scan is deflated ({S.LAMBDA_GC} = {lam:.3f}): the recessive coding leaves '
                 f'most variants with too few alt-homozygotes to fit, plink2 drops them on '
                 f'VIF_INFINITE, and what remains is tested against a null that is too narrow. '
                 f'Its peaks are not comparable with the additive ones.')
    elif lam > 1.05:
        calib = (f' {S.LAMBDA_GC} = {lam:.3f} is above 1 across the whole distribution, not just at '
                 f'the tail, which is residual structure rather than a handful of true signals; '
                 f'no genomic-control correction is applied to any reported statistic.')

    S.caption_block(
        fig,
        title=(f'{args.cohort} under the {mlab} model: {S.LAMBDA_GC} = {lam:.3f}, '
               f'{n_gw} genome-wide and {n_sug} suggestive peak'
               f'{"" if n_gw + n_sug == 1 else "s"}'
               + (f' ({gw_names}).' if n_gw else '.')),
        panels=[
            f'Manhattan of all {len(df):,} analysed variants, genome-wide loci named by gene; '
            f'red $5\\times10^{{-8}}$, dashed grey $1\\times10^{{-5}}$.',
            'Quantile\u2013quantile against the uniform null; grey band, pointwise 95% '
            'concentration.',
            f'This scan\'s peaks on the same axis as (a), merged at 250 kb; all {n_gw} '
            f'genome-wide and the {n_lab_sug} strongest of {n_sug} suggestive leads are named.',
            'Lead $|$effect$|$ against lead MAF, log axis.',
        ],
        notes=(f'{args.cohort}, {args.model} model, {args.covar_label}; '
               f'{len(df):,} variants analysed, {n_err:,} excluded. Full explanation, symbol '
               f'definitions and the estimator: {Path(args.out_png).stem}.md'),
        plot_h=PLOT_H, top_pad=0.42, hspace=0.30, wspace=0.30, left='auto', right=0.988,
        margin_axes=[ax_man, ax_pk])
    # After caption_block: both need the FINAL axes geometry to measure against.
    label_panels(ax_man, ax_pk, pk, off, args)
    # AFTER labelling: the strips have been taken out of (a) and (c), so bring
    # each row's second panel down to match, then set the panel letters with a pad
    # that clears the strip — the letter lives in the title slot, which follows the
    # axes, and would otherwise print below the gene names it introduces.
    # Bring each row's second panel down to the labelled panel's top. This was
    # what made (d) 1.00 in tall before — but the cause was row 2 being allocated
    # too little (height_ratio 1.18), not the alignment: at ratio 1.5 both panels
    # land near 1.9 in. Without it (d) stands 40 % taller than (c) and dominates a
    # row it is only the summary half of.
    S.align_panel_tops(ax_man, ax_qq)
    S.align_panel_tops(ax_pk, ax_eff)
    # One pad per ROW, the larger of the two, so the letters sit on one line even
    # though only the genomic panel gave up a strip.
    for row, letters in (((ax_man, ax_qq), 'ab'), ((ax_pk, ax_eff), 'cd')):
        pad = max(S.strip_pad(a) for a in row)
        for ax, letter in zip(row, letters):
            S.panel_tag(ax, letter, pad=pad)
    S.thin_tick_labels(ax_man, 'x')
    S.thin_tick_labels(ax_pk, 'x')
    fig.savefig(args.out_png)
    plt.close(fig)

    figure_doc.write_doc(
        args.out_png,
        title=f'Genome-wide scan — {args.cohort}, {mlab} model',
        question=('Where does this scan show association, what are its peaks, and is it calibrated '
                  'well enough for those signals to be believed?'),
        interpretation=(f'{n_err:,} of {n_all:,} variants carried a non-null ERRCODE or a degenerate '
                        f'fit and are excluded from every panel and from lambda_GC.{calib} The '
                        f'suggestive tier describes the shape of this scan and is not a list of '
                        f'findings; it receives no per-peak follow-up. Only the additive model '
                        f'defines the downstream fan-out — fine-mapping, conditional analysis and '
                        f'regional plots run on additive peaks alone — so dominant and recessive are '
                        f'reported here in full and go no further.'),
        panels=[
            ('a', 'Manhattan plot',
             'Every analysed variant at its genomic position, y = -log10 P. All points are drawn — no '
             'thinning — and the scatter is rasterised so the file stays small while axes and text '
             'remain vector. Thinning would distort the QQ plot in (c), which is read from the same '
             'P-values. Solid red line is genome-wide significance (5e-8); dashed grey is the '
             'suggestive threshold (1e-5). The x tick labels are on (b), which shares this axis.'),
            ('b', 'Peaks of this scan',
             'One marker per peak, at its lead variant, on the same x-axis as (a). Peaks are formed '
             'by merging variants that pass P < 1e-5 and lie within 250 kb; the lead is the '
             'smallest-P variant in the merged window. Distance merging is used rather than LD '
             'clumping so the peak definition needs no reference panel and cannot shift with the '
             'choice of LD sample. Labels: every genome-wide lead, plus the N smallest-P suggestive '
             'leads (N is --label-suggestive, default 10).'),
            ('c', 'Quantile-quantile plot',
             'The same P-values ranked against their expectation under a uniform null. The grey band '
             'is the pointwise 95% concentration from the order-statistic Beta(i, n-i+1) '
             'distribution, so departures are judged against a stated reference rather than by eye. '
             'lambda_GC = median(chi2)/0.4549 is printed in the panel.'),
            ('d', 'Effect against frequency',
             'MAF of each lead on x; |effect| as an odds ratio, exp|log OR|, on a log y axis. Taking '
             'the absolute value puts protective and risk alleles on one scale. The upward slope '
             'toward low MAF is a property of threshold selection, not of biology.'),
        ],
        numbers=[
            ('cohort', args.cohort), ('model', args.model),
            ('variants in file', n_all),
            ('excluded on ERRCODE / degenerate fit', n_err),
            ('variants analysed', len(df)),
            ('lambda_GC', float(lam)),
            ('variants P < 5e-8', n_hits),
            ('genome-wide peaks', n_gw),
            ('suggestive peaks', n_sug),
            ('suggestive peaks labelled', n_lab_sug),
            ('genome-wide loci', gw_names),
            ('annotation style', args.anno_style),
            ('repel force', float(args.repel_force)),
            ('smallest P', float(df['P'].min())),
        ],
        tables=[('Genome-wide peaks of this scan',
                 pk[pk.tier == GW][['peak_id', 'rsID', 'Gene', 'EA', 'OA', 'OR', 'L95', 'U95', 'P',
                                    'Case_Genotype_Distribution', 'Case_EAF',
                                    'Control_Genotype_Distribution', 'Control_EAF']]
                 if n_gw else None),
                ('Suggestive peaks (smallest P first)',
                 pk[pk.tier == SUG].nsmallest(15, 'P')[
                     ['peak_id', 'rsID', 'Gene', 'OR', 'L95', 'U95', 'P', 'MAF']]
                 if n_sug else None, 15)],
        reading=[
            'Check lambda_GC in (b) first. Near 1 means the bulk of the distribution is calibrated; '
            'the additive scans run ~1.12 and the recessive ~0.65, which is discussed in METHODS §7.',
            'In (b), a curve that lifts off the band only at the extreme right is what a true signal '
            'looks like. A curve that leaves the band along its whole length indicates residual '
            'structure rather than a handful of associations.',
            'Read (c) against (a): a peak marker sits directly below the Manhattan column it was '
            'called from, so the two can be compared without re-reading positions.',
            'Read the peaks as belonging to *this* model only. A locus shown on the recessive figure is '
            'a recessive result and is not carried into fine-mapping or conditional analysis, both of '
            'which run on the additive peaks alone.',
            'For a deflated scan, treat the peak count as an upper bound: the variants plink2 could '
            'fit are a biased subset of the call set.',
        ],
        limits=[
            'It cannot separate confounding from polygenicity. That needs LDSC, which is not usable '
            'at this effective sample size, or a random-effect model — see METHODS §7.',
            'Variants plink2 could not fit cleanly are absent from every panel, so the figure says '
            'nothing about them; their counts and error codes are in `02.scan/scan_qc.tsv`.',
            'A peak here is a threshold crossing in one scan of one nested sample set. It is not a '
            'replicated finding, and no independent cohort is available to this study.',
            'Panel (d) compares effects across frequencies but against no detection reference; this '
            'component makes no power claim.',
        ],
        defs=['model', 'lambda_gc', 'gw_sig', 'or', 'errcode'],
        model=S.glm_formula(args.pc_label, args.n_pcs) + '\n' + S.FORMULAS['lambda'],
        methods_ref='../../../docs/METHODS.md')
    print(f'[plot_manhattan_qq] {args.cohort}/{args.model}: lambda={lam:.4f} '
          f'n={len(df):,} err={n_err:,} hits={n_hits} peaks={n_gw}gw/{n_sug}sug '
          f'-> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in plot_manhattan_qq: {e}', file=sys.stderr)
        sys.exit(1)
