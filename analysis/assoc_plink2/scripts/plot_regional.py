#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : One regional plot per genome-wide locus, with the SAME association
#           statistics shown three times, coloured by LD from three independent
#           sources on one common binned scale:
#             (a) in-sample cohort LD   (b) population co-occurrence panel
#             (c) out-of-sample reference panel   — names come from --panel-label
#           Putting them in one figure with a shared x-axis makes LD concordance
#           readable at a glance; the reference produced three separate PNGs the
#           reader had to compare by eye.
#           Recombination rate is overlaid on each panel and a gene track sits
#           beneath, so the credible interval can be read against gene structure.
# Component: assoc_plink2
# Used by : assoc_plink2.nf  process PLOT_REGIONAL
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
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_style as S
import figure_doc

PLOT_H = 6.6          # inches of croppable plot block (titles + axes + x-labels)

# Panel titles are DEFAULTS for the three LD roles; the configured resource for
# each is named by --panel-label so no particular panel is written into the code.
SOURCE_TITLE = {
    'cohort':    'in-sample cohort LD',
    'tommo':     'population co-occurrence panel',
    '1000g_eas': 'out-of-sample reference panel',
}
# Background context, not a measurement: the recombination trace shares the panel
# with the association points and at full weight its spikes read as signal.
RECOMB_COLOR = '#9DBFD8'


def parse_args():
    p = argparse.ArgumentParser(description='Three-source regional plot for one locus.')
    p.add_argument('--sumstat', required=True)
    p.add_argument('--ld-dir', required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--locus-id', required=True)
    p.add_argument('--lead-id', required=True)
    p.add_argument('--exons')
    p.add_argument('--recomb')
    p.add_argument('--alpha', type=float, default=S.GW_ALPHA)
    p.add_argument('--panel-label', action='append', default=[], metavar='KEY=LABEL',
                   help='Display name per LD source key, e.g. tommo="gnomAD NFE".')
    p.add_argument('--out-png', required=True)
    return p.parse_args()


def read_ld(ld_dir, locus, source):
    f = Path(ld_dir) / f'{locus}.ld_{source}.tsv'
    if not f.exists():
        return None
    d = pd.read_csv(f, sep='\t', dtype={'ID': str})
    return d.set_index('ID')


def draw_panel(ax, ss, ld, source, args, ymax, label_lead=False):
    x = ss['POS'].to_numpy() / 1e6
    y = -np.log10(ss['P'].to_numpy())

    if ld is None:
        r2 = np.full(len(ss), np.nan)
        known = np.zeros(len(ss), dtype=bool)
        counts = {'measured': 0, 'below_threshold': 0, 'not_in_panel': len(ss)}
    else:
        al = ld.reindex(ss['ID'])
        r2 = al['r2'].to_numpy(dtype=float)
        state = al['state'].fillna('not_in_panel').to_numpy()
        known = state != 'not_in_panel'
        counts = pd.Series(state).value_counts().to_dict()

    colors = S.ld_bin_colors(r2, known=known)
    order = np.argsort(np.nan_to_num(np.where(known, r2, -1.0)))   # high LD drawn last
    ax.scatter(x[order], y[order], s=8, c=colors[order], linewidths=0.3,
               edgecolors='white', rasterized=True, zorder=3)

    lead = ss['ID'] == args.lead_id
    if lead.any():
        ax.scatter(x[lead], y[lead], s=34, marker='D', c=S.LEAD_COLOR,
                   edgecolors='white', linewidths=0.8, zorder=5)
        # Name the lead ONCE. It is the same variant in all three panels, so
        # repeating it three times only crowds the titles.
        if label_lead:
            ax.annotate(args.lead_id, xy=(x[lead][0], y[lead][0]), xytext=(0, 9),
                        textcoords='offset points', ha='center', fontweight='bold',
                        fontsize=plt.rcParams['legend.fontsize'] - 1,
                        color=S.LEAD_COLOR, zorder=6)

    ax.axhline(-np.log10(args.alpha), color=S.ACCENT, lw=1.1, zorder=2)
    ax.set_ylim(0, ymax)
    ax.set_ylabel(S.NEGLOG10P)
    S.despine(ax, grid_axis='y')

    # No in-panel coverage sentence. It is provenance, not a finding, and it sat
    # over the data at the top right of all three LD panels. The same counts go to
    # the sidecar's "Values in this rendering" table, which is where provenance
    # belongs; the `state` semantics are unchanged (below_threshold still draws in
    # the lowest bin, only not_in_panel draws grey — METHODS 10c).
    return counts


def draw_recomb(ax, recomb, xlim):
    """Recombination rate on a right-hand twin axis, drawn behind the points."""
    if recomb is None or not len(recomb):
        return None
    tw = ax.twinx()
    xm = (recomb['start'].to_numpy() + recomb['end'].to_numpy()) / 2e6
    tw.plot(xm, recomb['rate'].to_numpy(), color=RECOMB_COLOR, lw=0.8, alpha=0.55, zorder=1)
    tw.set_ylim(0, max(float(recomb['rate'].max()) * 1.15, 1.0))
    tw.set_xlim(*xlim)
    tw.set_ylabel('cM/Mb', color=RECOMB_COLOR, fontsize=plt.rcParams['axes.labelsize'] - 2)
    tw.tick_params(axis='y', colors=RECOMB_COLOR,
                   labelsize=plt.rcParams['ytick.labelsize'] - 2)
    tw.grid(False)
    tw.spines['top'].set_visible(False)
    tw.spines['right'].set_color(RECOMB_COLOR)
    tw.set_zorder(ax.get_zorder() - 1)
    ax.patch.set_visible(False)
    return tw


GENE_COLOR = '#33556E'
GENE_LABEL_COLOR = '#1B2B38'
EXON_H = 0.34          # exon box height, in row units


def draw_gene_track(ax, exons, xlim, chrom):
    """Exon/intron structure, one row per gene.

    Introns are a thin line carrying strand chevrons; exons are filled boxes on
    top of it. One representative transcript per gene (the longest coding one,
    chosen upstream in gene_utils.R), because drawing every transcript of a gene
    would need 22 rows to say the same thing.

    Only genes with an official symbol reach this function — clone-accession
    models such as AC007347.1 and RP11-357N13.3 are filtered in gene_utils.R,
    where the rule is shared with the lead-variant gene assignment.
    """
    ax.set_xlim(*xlim)
    ax.set_yticks([])
    ax.grid(False)
    for sp in ('top', 'right', 'left'):
        ax.spines[sp].set_visible(False)
    ax.set_xlabel(f'chr{str(chrom).replace("chr", "")} {S.POS_MB}')

    if exons is None or not len(exons):
        ax.text(0.5, 0.5, 'no named gene in window', transform=ax.transAxes,
                ha='center', va='center', color=S.REFERENCE,
                fontsize=plt.rcParams['legend.fontsize'] - 1)
        ax.set_ylim(0, 1)
        return 0

    span = xlim[1] - xlim[0]
    fs = plt.rcParams['legend.fontsize'] - 2.5
    genes = (exons.groupby('gene', sort=False)
             .agg(g0=('gene_start', 'min'), g1=('gene_end', 'max'),
                  strand=('strand', 'first'))
             .reset_index().sort_values('g0'))

    rows = []                                    # occupied right edge per row
    ypos = {}
    for _, g in genes.iterrows():
        x0, x1 = g.g0 / 1e6, g.g1 / 1e6
        need = max(x1, x0 + len(str(g.gene)) * span * 0.012) + span * 0.015
        for i, edge in enumerate(rows):
            if x0 > edge:
                rows[i] = need
                ypos[g.gene] = -i
                break
        else:
            rows.append(need)
            ypos[g.gene] = -(len(rows) - 1)

    # Intron line first, then chevrons placed ONLY inside introns wide enough to
    # hold one, then exons on top. Spacing chevrons by a fixed fraction of the
    # window instead would drop them on top of the exons of a compact gene —
    # RPGRIP1L packs 25 exons into 106 kb.
    min_exon_w = span * 0.0035        # so a 1 kb exon is still a visible tick
    chev_gap = span * 0.022
    for _, g in genes.iterrows():
        y = ypos[g.gene]
        x0, x1 = g.g0 / 1e6, g.g1 / 1e6
        ax.plot([x0, x1], [y, y], lw=1.1, color=GENE_COLOR, zorder=3,
                solid_capstyle='butt')

        ex = exons[exons['gene'] == g.gene].sort_values('exon_start')
        drawn = [(max(e.exon_start / 1e6 - min_exon_w / 2, x0),
                  max(e.exon_end / 1e6, e.exon_start / 1e6 + min_exon_w))
                 for _, e in ex.iterrows()]
        gaps = [(drawn[i][1], drawn[i + 1][0]) for i in range(len(drawn) - 1)]
        gaps += [(x0, drawn[0][0]), (drawn[-1][1], x1)]
        xs = [(a + b) / 2 for a, b in gaps
              if b - a > chev_gap and xlim[0] < (a + b) / 2 < xlim[1]]
        if xs:
            ax.plot(xs, np.full(len(xs), y), ls='none', zorder=4,
                    marker=('>' if g.strand == '+' else '<'),
                    ms=3.6, mfc=GENE_COLOR, mec=GENE_COLOR)

        for a, b in drawn:
            ax.add_patch(Rectangle((a, y - EXON_H / 2), b - a, EXON_H,
                                   facecolor=GENE_COLOR, edgecolor='none', zorder=5))

    above, beside = [], []
    for _, g in genes.iterrows():
        y = ypos[g.gene]
        x0, x1 = g.g0 / 1e6, g.g1 / 1e6
        vis0, vis1 = max(x0, xlim[0]), min(x1, xlim[1])
        # A gene wider than the plotted window has its right end
        # off-canvas, so a label anchored there disappears; put it above the
        # visible span instead.
        if x1 > xlim[1] - span * 0.02 or x0 < xlim[0] + span * 0.02:
            above.append(ax.annotate(
                g.gene, xy=((vis0 + vis1) / 2, y), xytext=(0, 8),
                textcoords='offset points', va='bottom', ha='center',
                fontsize=fs, fontstyle='italic', color=GENE_LABEL_COLOR, zorder=7))
        else:
            beside.append(ax.annotate(
                g.gene, xy=(x1, y), xytext=(4, 0), textcoords='offset points',
                va='center', ha='left', fontsize=fs,
                fontstyle='italic', color=GENE_LABEL_COLOR, zorder=7))
    # Row packing keeps two genes on one row apart, but a label placed ABOVE its
    # gene reaches into the row above, where another gene's own label may already
    # be. Separate those on measured boxes; the beside-labels sit inside their own
    # row and are already accounted for by the packer.
    ax.set_ylim(-len(rows) + 0.45, 0.9)
    S.spread_labels(ax, above, axis='x', pad_px=3.0)
    return len(genes)


def main():
    args = parse_args()
    S.setup_style('paper')
    ss = pd.read_csv(args.sumstat, sep='\t', dtype={'CHROM': str, 'ID': str})
    ss = ss[ss['P'].notna() & (ss['P'] > 0)]
    chrom = ss['CHROM'].iloc[0]
    xlim = (ss['POS'].min() / 1e6, ss['POS'].max() / 1e6)
    ymax = max(9.0, -np.log10(ss['P'].min()) * 1.22)

    exons = recomb = None
    if args.exons and Path(args.exons).exists():
        exons = pd.read_csv(args.exons, sep='\t')
        if not len(exons):
            exons = None
    if args.recomb and Path(args.recomb).exists():
        recomb = pd.read_csv(args.recomb, sep='\t')

    fig, axes = plt.subplots(4, 1, figsize=(S.COL_DOUBLE, PLOT_H), sharex=False,
                             gridspec_kw={'height_ratios': [1, 1, 1, 0.70]})
    all_counts = {}
    TITLE = dict(SOURCE_TITLE, **dict(s.split('=', 1) for s in args.panel_label if '=' in s))
    for ltr, (ax, src) in zip('abc', zip(axes[:3], ['cohort', 'tommo', '1000g_eas'])):
        ld = read_ld(args.ld_dir, args.locus_id, src)
        draw_recomb(ax, recomb, xlim)
        all_counts[src] = draw_panel(ax, ss, ld, src, args, ymax, label_lead=(ltr == 'a'))
        ax.set_xlim(*xlim)
        ax.tick_params(labelbottom=False)
        S.panel_tag(ax, ltr, title=TITLE[src])
    n_genes = draw_gene_track(axes[3], exons, xlim, chrom)
    S.panel_tag(axes[3], 'd')

    handles = S.ld_legend_handles(include_unknown=True, lead=True)
    if recomb is not None and len(recomb):
        handles.append(Line2D([0], [0], color=RECOMB_COLOR, lw=1.6, label='recombination rate'))
    # The 8-key strip needs room ABOVE it (canvas edge) and BELOW it (panel (a)'s
    # centred title). At y=1.22 with top_pad=0.42 in it rendered flush against the
    # top of the canvas and read as cut off; top_pad is raised to 0.85 in below.
    S.legend_above(axes[0], handles, title=f'{S.R2_LD} with the lead variant',
                   ncol=len(handles), y=1.30)

    tm = all_counts.get('tommo', {})
    S.caption_block(
        fig,
        title=(f'{args.cohort}, locus {args.locus_id}: identical association statistics coloured by LD '
               f'with {args.lead_id} from three independent sources.'),
        panels=[
            'In-sample LD — the same samples the statistics come from, and the LD SuSiE uses.',
            (f'{TITLE["tommo"]} co-occurrence; pairs below $r^{{2}}=0.2$ are unpublished, so those '
             f'{int(tm.get("below_threshold", 0)):,} variants are bounded into the lowest bin, not '
             f'missing. Grey is absent from the panel entirely.'),
            f'{TITLE["1000g_eas"]} — an out-of-sample population reference.',
            'Gene models (Ensembl 86): boxes are exons, chevrons give the strand, one representative '
            'transcript per gene.',
        ],
        notes=(f'Window chr{str(chrom).replace("chr", "")}:{int(ss.POS.min()):,}\u2013'
               f'{int(ss.POS.max()):,} ({len(ss):,} variants). Purple diamond, lead. Blue trace, '
               f'recombination rate (right axis). Full explanation: {Path(args.out_png).stem}.md'),
        plot_h=PLOT_H, top_pad=0.85, hspace=0.34, left='auto', right=0.918, margin_axes=axes)
    fig.savefig(args.out_png)
    plt.close(fig)

    cov = []
    for src, c in all_counts.items():
        cov += [(f'{src}: r2 measured', int(c.get('measured', 0))),
                (f'{src}: bounded r2 < 0.2', int(c.get('below_threshold', 0))),
                (f'{src}: not in panel', int(c.get('not_in_panel', 0)))]
    figure_doc.write_doc(
        args.out_png,
        title=f'Regional association — {args.cohort}, peak {args.locus_id}',
        question=('Is the LD structure that produces this peak a property of the East Asian '
                  'population, or an artefact of this particular sample?'),
        panels=[
            ('a', 'In-sample cohort LD',
             'The association statistics of the peak window, each variant coloured by r2 with the '
             'lead computed on the *same* samples the statistics came from. This is the LD matrix the '
             'SuSiE fine-mapping uses, so this panel shows exactly what fine-mapping saw.'),
            ('b', TITLE['tommo'],
             f'Identical statistics, recoloured by r2 from the {TITLE["tommo"]} co-occurrence tables. That '
             'resource publishes only pairs with r2 >= 0.2, so a variant present in the panel but '
             'returning no pair with the lead is *bounded* below 0.2 — it belongs in the lowest colour '
             'bin, not in grey. Grey is reserved for variants absent from the panel entirely. '
             'Treating the bound as missing data is what made this source look uninformative in the '
             'earlier implementation.'),
            ('c', TITLE['1000g_eas'],
             f'Identical statistics again, coloured by r2 in {TITLE["1000g_eas"]} — an out-of-sample '
             'population reference with no relationship to this cohort.'),
            ('d', 'Gene models',
             'Ensembl 86 / GRCh38. Filled boxes are exons, the connecting line spans the introns and '
             'chevrons give the transcribed strand. One representative transcript per gene (the '
             'longest protein-coding one), because drawing every transcript of a multi-transcript gene would '
             'need 22 rows to say the same thing. Only genes with an official symbol appear; '
             'clone-accession models such as AC007347.1 or RP11-357N13.3 name a sequencing clone '
             'rather than a gene and are suppressed.'),
        ],
        numbers=[('cohort', args.cohort), ('peak', args.locus_id), ('lead variant', args.lead_id),
                 ('window', f"chr{str(chrom).replace('chr','')}:{int(ss.POS.min()):,}-{int(ss.POS.max()):,}"),
                 ('variants in window', len(ss)),
                 ('genes drawn', int(n_genes))] + cov,
        reading=[
            'Compare the three colour patterns. If they agree, the LD block is a population property '
            'and the credible set can be trusted to reflect real correlation structure.',
            'If the in-sample panel alone shows tight LD, the correlation is being driven by sampling '
            'noise at this N and the credible set is correspondingly fragile.',
            'Read the lead against the gene track: whether it sits in an exon, an intron or between '
            'genes constrains which mechanisms are plausible.',
        ],
        limits=[
            'It does not identify a causal variant. LD colour is correlation with the lead, not '
            'evidence of function — that is what the fine-mapping figure addresses.',
            'The two external panels differ from each other and from the study samples in ancestry '
            'breadth and in genomic coverage. Disagreement between (b) and (c) may reflect that '
            'difference rather than an error in either.',
        ],
        defs=['ld_r2', 'gw_sig'],
        interpretation=('Agreement between the three panels means the LD structure is a property of the '
                   'population rather than of this sample; divergence localised to the in-sample panel would '
                   'indicate that the credible set is being shaped by sampling noise at this N. Each panel '
                   'prints its own coverage, so a sparse source cannot be mistaken for a low-LD region.'),
        methods_ref='../../../docs/METHODS.md')
    print(f'[plot_regional] {args.cohort} {args.locus_id} -> {args.out_png}')
    for k, v in all_counts.items():
        print(f'    {k:10s} {v}')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in plot_regional: {e}', file=sys.stderr)
        sys.exit(1)
