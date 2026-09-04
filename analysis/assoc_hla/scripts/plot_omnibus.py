#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The amino-acid OMNIBUS figure for one cohort x one model -- the
#           layer that makes this HLA fine mapping rather than a marker scan.
#           hla_omnibus.py asks, at every IMGT position, whether the residue
#           CONTENT of that position associates with disease (an m-1 df
#           likelihood-ratio test over the position's residues, Hirata et al.
#           Nat. Genet. 2019 / Nat. Protoc. 2023). This figure is that table.
#
#           WHY ONE SUB-AXES PER GENE AND NOT ONE LONG AXIS. The x quantity is
#           an IMGT AMINO-ACID POSITION, and IMGT numbering is per gene: DQB1:57
#           and A:57 are different residues of different proteins that happen to
#           share an integer. Laid end to end on one axis they would read as one
#           coordinate system, which is exactly the misreading the layout must
#           make impossible. Each gene therefore gets its own x-axis, and the
#           genes are ordered BY THEIR CHR6 POSITION (--gene-order) so the grid
#           still reads left-to-right, top-to-bottom along the MHC, matching
#           plot_hla_scan.py's shared axis.
#
#           WHY ONE SHARED Y LIMIT. -log10 P is the same quantity in every
#           panel, so a taller point has to be a stronger point. Per-gene limits
#           would rescale every panel to its own maximum and make the weakest
#           gene look exactly like the strongest.
#
#           THE SHADED BAND IS THE PEPTIDE-BINDING DOMAIN. Class I molecules
#           bind peptide through the alpha-1/alpha-2 domains, mature residues
#           1-180; class II through the single alpha-1 or beta-1 domain, mature
#           residues 1-90. That is where a FUNCTIONAL signal is expected, and
#           the band lets the reader see whether a hit is there or somewhere
#           with no binding-site interpretation -- without the figure asserting
#           anything about the hit, which is the caption's job. Positions are
#           mature-protein numbered, so the leader peptide is NEGATIVE and sits
#           left of the band by construction.
#
#           NOTHING IS SILENTLY DROPPED. A gene with no eligible position gets
#           an EMPTY AXES carrying a short note, never omission: an absent panel
#           is indistinguishable from a gene that was never in the analysis, and
#           the two mean opposite things. Positions that were tested but whose
#           fit was not usable (status != 'ok': separation, non-convergence, too
#           few complete cases) carry no P and cannot be plotted as points, so
#           they are drawn as a rug of ticks along the axis floor and counted in
#           the caption.
# Component: assoc_hla
# Used by : assoc_hla.nf  process PLOT_OMNIBUS
# ---------------------------------------------------------------------------
import argparse
import string
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt              # noqa: E402
import numpy as np                           # noqa: E402
import pandas as pd                          # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_shared' / 'scripts'))
import plot_style as S                       # noqa: E402
import figure_doc                            # noqa: E402

# The plot block GROWS WITH THE GRID rather than being fixed: the gene count is
# a property of the run (8 primary-assembly classical genes, more in a
# sensitivity arm that keeps DRB3/4/5), and a fixed height would crush the rows
# instead of adding one. Each row gets a stated number of inches and the block
# is the sum, so panels are the same size on every rendering.
ROW_H = 1.45            # inches per grid row (axes + its x tick labels)
BASE_H = 0.62           # inches above the first row: panel letters + gene titles
NCOLS_DEFAULT = 4

OMNIBUS_COLS = ['cohort', 'model', 'position', 'gene', 'n_residues', 'df', 'n',
                'reference_residue', 'residues', 'determined_rate', 'deviance', 'P', 'status']
QC_COLS = ['id', 'marker_class', 'gene', 'pos', 'mac', 'maf', 'in_reference', 'tested']

# Peptide-binding domain, in MATURE-PROTEIN residue numbers — the numbering the
# IMGT position field uses, which is why the leader peptide comes out negative.
# Class I binds through alpha-1 + alpha-2 (1-180); class II through the single
# alpha-1 or beta-1 domain (1-90).
PBD = {'I': (1, 180), 'II': (1, 90)}
CLASS_I = {'A', 'B', 'C', 'E', 'F', 'G'}


def parse_args():
    p = argparse.ArgumentParser(
        description='Amino-acid omnibus figure for one cohort x one model: one sub-axes '
                    'per HLA gene, omnibus -log10 P against IMGT residue position.')
    p.add_argument('--omnibus', required=True, help='hla_omnibus.py output')
    p.add_argument('--marker-qc', required=True,
                   help='hla_marker_qc.py output; gives the gene order fallback and the '
                        'tested-residue counts quoted in the sidecar')
    p.add_argument('--cohort', required=True)
    p.add_argument('--model', required=True, help='fixed | random')
    p.add_argument('--model-label', default='',
                   help='human-readable model name for the header; defaults to --model')
    p.add_argument('--alpha', type=float, required=True,
                   help='significance threshold for the omnibus tests, e.g. 0.05 / positions')
    p.add_argument('--gene-order', default='',
                   help='comma-separated gene symbols in chr6 order, e.g. A,C,B,DRB1,DQA1,'
                        'DQB1,DPA1,DPB1. Genes not named here are appended in marker-QC '
                        'position order; genes named but absent from the omnibus still get '
                        'an axes, so a gene can never disappear silently.')
    p.add_argument('--ncols', type=int, default=NCOLS_DEFAULT,
                   help='columns in the gene grid')
    p.add_argument('--label-top', type=int, default=4,
                   help='maximum positions named per gene, smallest P first')
    p.add_argument('--anno-style', default='auto', choices=sorted(S.ANNO_STYLES),
                   help='gwaslab label style passed to plot_style.gene_labels')
    p.add_argument('--out-png', required=True)
    return p.parse_args()


# ── Loading ──────────────────────────────────────────────────────────────────
def load_omnibus(path, args):
    """hla_omnibus.py's table, with the position string split into gene + IMGT number.

    `position` is 'GENE:N' and N can be NEGATIVE (the leader peptide, e.g.
    'A:-22'), so it is parsed on the LAST colon and read as a signed integer —
    stripping a '-' as if it were a separator would silently fold the leader
    onto the mature chain.
    """
    if not Path(path).exists():
        raise SystemExit(f'ABORT: {path} does not exist')
    d = pd.read_csv(path, sep='\t', dtype={'cohort': str, 'model': str, 'position': str,
                                           'gene': str, 'status': str,
                                           'reference_residue': str, 'residues': str})
    for c in ('position', 'gene', 'status'):
        if c not in d.columns:
            raise SystemExit(f'ABORT: {path} has no {c!r} column; it is not an hla_omnibus '
                             f'table')
    for c in ('n_residues', 'df', 'n', 'determined_rate', 'deviance', 'P'):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors='coerce')
    for col, want in (('cohort', args.cohort), ('model', args.model)):
        if col in d.columns and len(d) and set(d[col].dropna().unique()) - {want}:
            print(f'[plot_omnibus] WARNING: {path} carries {col} '
                  f'{sorted(set(d[col].dropna().unique()))}, figure is labelled {want!r}',
                  file=sys.stderr)
    d['gene'] = d['gene'].fillna('').astype(str)
    d['status'] = d['status'].fillna('').astype(str)
    d['aa_pos'] = [_imgt_pos(p) for p in d['position'].fillna('')]
    d = d[pd.notna(d['aa_pos'])].copy()
    d['aa_pos'] = d['aa_pos'].astype(int)
    # Only status == 'ok' rows carry a P. Everything else is a position that WAS
    # tested and whose fit could not be used; it has no y and is drawn as a rug.
    d['usable'] = (d['status'] == 'ok') & d['P'].notna() & (d['P'] > 0) & (d['P'] <= 1)
    return d


def _imgt_pos(pos):
    """'DQB1:57' -> 57, 'A:-22' -> -22, anything else -> None."""
    s = str(pos)
    if ':' not in s:
        return None
    tail = s.rsplit(':', 1)[1].strip()
    try:
        return int(tail)
    except ValueError:
        return None


def load_qc(path):
    d = pd.read_csv(path, sep='\t', dtype={'id': str, 'marker_class': str, 'gene': str})
    for c in ('id', 'marker_class', 'gene'):
        if c not in d.columns:
            raise SystemExit(f'ABORT: {path} has no {c!r} column')
    d = d[[c for c in QC_COLS if c in d.columns]].copy()
    for c in ('pos', 'mac', 'maf', 'in_reference'):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors='coerce')
    d['tested'] = (d['tested'].astype(str).str.lower().isin(('true', '1'))
                   if 'tested' in d.columns else True)
    return d


def gene_sequence(args, om, qc):
    """The genes to draw, in chr6 order, and where that order came from.

    --gene-order is authoritative because chr6 order is a property of the
    assembly and not of this cohort's marker set. What it does NOT get to do is
    shorten the figure: any gene present in the omnibus but unnamed is appended
    (in marker-QC position order), and any gene named but absent still gets its
    axes and its note.
    """
    named = [g.strip() for g in args.gene_order.split(',') if g.strip()]
    by_pos = []
    if 'pos' in qc.columns and len(qc):
        by_pos = list(qc.groupby('gene')['pos'].min().sort_values().index)
    extra = [g for g in by_pos if g and g not in named]
    extra += sorted(g for g in om['gene'].unique()
                    if g and g not in named and g not in extra)
    order = named + extra
    src = 'the --gene-order argument' if named else 'marker-QC positions'
    if not order:
        order = sorted(g for g in om['gene'].unique() if g)
        src = 'the omnibus table'
    return order, src, extra if named else []


def hla_class(gene):
    """Class I or class II, from the symbol. Anything else has no PBD to shade."""
    g = str(gene).upper()
    if g in CLASS_I:
        return 'I'
    if g.startswith('D'):
        return 'II'
    return None


# ── Panels ───────────────────────────────────────────────────────────────────
def gene_panel(ax, d, gene, args, ylim, thr):
    """One gene: omnibus -log10 P against IMGT amino-acid position.

    Returns the rows named in this panel, so the caption and the sidecar count
    the same labels the figure drew.
    """
    S.despine(ax, grid_axis='y')
    ok = d[d['usable']] if len(d) else d
    bad = d[~d['usable']] if len(d) else d

    if not len(ok) and not len(bad):
        # An EMPTY AXES WITH A NOTE, never a missing panel. A gene that vanished
        # from the grid is indistinguishable from a gene that was never in the
        # analysis, and the two mean opposite things. The peptide-binding band is
        # NOT drawn here: with no data there is no x-axis to read it against, and
        # a shaded region floating over blank axes is decoration that looks like
        # a measurement.
        ax.text(0.5, 0.5, 'no position\neligible', transform=ax.transAxes,
                ha='center', va='center', color=S.REFERENCE,
                fontsize=plt.rcParams['legend.fontsize'] - 0.5, linespacing=1.4)
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ('left', 'bottom'):
            ax.spines[sp].set_visible(False)
        ax.grid(False)
        return ok.head(0)

    cls = hla_class(gene)
    if cls:
        lo, hi = PBD[cls]
        # Behind the data and carrying no measurement, so NEUTRAL: it is ground,
        # not a datum.
        ax.axvspan(lo, hi, color=S.NEUTRAL, alpha=0.55, linewidth=0, zorder=0)

    ax.set_ylim(*ylim)
    ax.axhline(thr, color=S.ACCENT, lw=1.0, zorder=3)
    if len(ok):
        ax.scatter(ok['aa_pos'].to_numpy(), ok['y'].to_numpy(), s=11,
                   c=S.DATA_DARK, linewidths=0, zorder=4)
    if len(bad):
        # Tested, no usable fit, therefore NO P — so it cannot be a point. The
        # rug lives in a sliver BELOW zero, which is why ylim has a negative
        # floor: drawn at y = 0 it would sit inside the crowd of null positions
        # and read as a measurement of about P = 1, which is the one thing these
        # positions did not return.
        y0, y1 = ylim[0] * 0.92, ylim[0] * 0.28
        ax.vlines(bad['aa_pos'].to_numpy(), y0, y1, color=S.NEUTRAL_D, lw=0.7, zorder=2)

    lo_x = int(min([*d['aa_pos'], 1]))
    hi_x = int(max([*d['aa_pos'], 1]))
    if hi_x == lo_x:                 # a single position: give it an axis to sit on
        lo_x, hi_x = lo_x - 5, hi_x + 5
    pad = max(2.0, 0.03 * (hi_x - lo_x))
    ax.set_xlim(lo_x - pad, hi_x + pad)
    ax.xaxis.set_major_locator(plt.MaxNLocator(4, integer=True))
    return ok[ok['P'] < args.alpha].nsmallest(
        min(args.label_top, int((ok['P'] < args.alpha).sum())), 'P') if len(ok) else ok


def label_panel(ax, sel, args):
    """Name the positions clearing --alpha, in a strip ABOVE the data area.

    Only the IMGT number is written: the gene is already the panel's title, and
    repeating 'DQB1:' on every label in the DQB1 panel would cost two thirds of
    each label's width to say what the reader is already looking at.
    MUST run after caption_block, which resizes the figure.
    """
    if not len(sel):
        return
    sel = sel.sort_values('aa_pos')
    S.gene_labels(ax, sel['aa_pos'].to_numpy(dtype=float), sel['y'].to_numpy(),
                  [str(int(p)) for p in sel['aa_pos']],
                  anno_style=args.anno_style, color=S.ACCENT, weight='bold',
                  fontsize=plt.rcParams['legend.fontsize'] - 0.5)


# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    S.setup_style('paper')
    mlab = args.model_label or args.model
    thr = -np.log10(args.alpha) if 0 < args.alpha <= 1 else None
    if thr is None:
        raise SystemExit(f'ABORT: --alpha must be in (0, 1]; got {args.alpha}')

    qc = load_qc(args.marker_qc)
    om = load_omnibus(args.omnibus, args)
    om['y'] = np.where(om['usable'], -np.log10(om['P'].where(om['usable'], 1.0)), np.nan)
    genes, order_src, appended = gene_sequence(args, om, qc)
    if not genes:
        raise SystemExit('ABORT: no gene to draw — the omnibus table is empty and '
                         '--gene-order named nothing. An empty grid would publish silently.')

    # One shared y limit: -log10 P is the same quantity in every panel, so a
    # taller point must be a stronger point. The limit also has to CONTAIN the
    # threshold line, or a gene with no signal shows a line-less panel and its
    # points read as though nothing was being asked of them.
    ymax = float(om.loc[om['usable'], 'y'].max()) if om['usable'].any() else 0.0
    ytop = max(ymax * 1.15, thr + 0.8, 2.0)
    # The floor is negative on purpose: it is the sliver the "tested, no usable
    # fit" rug lives in, so a position with no P is never drawn where a P of
    # about 1 would be.
    ylim = (-0.055 * ytop, ytop)

    ncols = max(1, min(args.ncols, len(genes)))
    nrows = int(np.ceil(len(genes) / ncols))
    plot_h = BASE_H + nrows * ROW_H
    fig = plt.figure(figsize=(S.COL_DOUBLE, plot_h))
    # sharey, not just a shared limit: the grid is one scale, and a shared
    # locator keeps two rows differing by hundredths of an inch from flipping to
    # a different tick step and reading as two different axes.
    gs = fig.add_gridspec(nrows, ncols)
    axes, sels = [], {}
    for i, gene in enumerate(genes):
        r, c = divmod(i, ncols)
        ax = fig.add_subplot(gs[r, c], sharey=axes[0] if axes else None)
        axes.append(ax)
        d = om[om['gene'] == gene]
        sels[gene] = gene_panel(ax, d, gene, args, ylim, thr)
        if c == 0:
            ax.set_ylabel(S.NEGLOG10P)
        elif len(d):
            ax.tick_params(labelleft=False)
        # The x label goes on the LAST panel of each column, which is the bottom
        # one that exists — a ragged final row would otherwise leave its
        # columns' axes unnamed.
        if i + ncols >= len(genes) and len(d):
            ax.set_xlabel('IMGT residue position')

    # Ticks are chosen ONCE and inherited through sharey. Two rows differing by
    # hundredths of an inch otherwise flip to a different tick step and stop
    # reading as one scale; and the negative rug sliver must never carry a tick,
    # because there is no such thing as a negative -log10 P.
    ticks = plt.MaxNLocator(nbins=4, integer=True).tick_values(0, ytop)
    axes[0].set_yticks([t for t in ticks if 0 <= t <= ytop])

    n_ok = int(om['usable'].sum())
    n_bad = int((~om['usable']).sum())
    n_sig = int((om.loc[om['usable'], 'P'] < args.alpha).sum())
    sig = om[om['usable'] & (om['P'] < args.alpha)].sort_values('P')
    top = sig.iloc[0] if len(sig) else (
        om[om['usable']].nsmallest(1, 'P').iloc[0] if n_ok else None)
    in_band = 0
    if len(sig):
        for r in sig.itertuples():
            cls = hla_class(r.gene)
            if cls and PBD[cls][0] <= r.aa_pos <= PBD[cls][1]:
                in_band += 1

    # One caption line per gene panel. Terse on purpose: at 8-11 panels a
    # sentence each is a wall of text, and everything shared between them (what
    # the band is, what the rug is) belongs in the Data line, once.
    panel_lines = []
    for gene in genes:
        d = om[om['gene'] == gene]
        good = d[d['usable']]
        if not len(d):
            panel_lines.append(f'{gene}, no eligible position.')
            continue
        best = good.nsmallest(1, 'P').iloc[0] if len(good) else None
        panel_lines.append(
            f'{gene}, {len(good)} position' + ('' if len(good) == 1 else 's')
            + (f', strongest {int(best.aa_pos)} ({S.p_tex(float(best.P))})' if best is not None
               else ', none with a usable fit') + '.')

    hit_txt = (f'{n_sig} position' + ('' if n_sig == 1 else 's')
               + f' clear {S.mathsci(args.alpha)}'
               + (f', {in_band} of them inside a peptide-binding domain' if n_sig else '')
               if n_sig else f'no position clears {S.mathsci(args.alpha)}')
    top_txt = ('' if top is None else
               f'; strongest {top["position"]} (df = {int(top["df"]) if pd.notna(top["df"]) else "?"}, '
               f'P = {S.p_tex(float(top["P"]))})')

    S.caption_block(
        fig, plot_h=plot_h,
        title=f'{args.cohort}, HLA {mlab} model: {hit_txt}{top_txt}.',
        panels=panel_lines,
        # More than eight panels is the normal case here — eight classical genes
        # is already the whole default alphabet, and a sensitivity arm keeping
        # DRB3/4/5 needs eleven. The default 'abcdefgh' raises IndexError there.
        letters=string.ascii_lowercase,
        notes=(f'{args.cohort}, {mlab} model; omnibus likelihood-ratio test over the residues '
               f'at each IMGT position, deviance on m-1 df. {n_ok:,} position(s) with a usable '
               f'fit, {n_bad:,} tested without one (status != ok: separation, non-convergence '
               f'or too few complete cases) — those carry no P and are drawn as grey ticks '
               f'along the axis floor, not as points. Red line {S.sci(args.alpha)}. The shaded '
               f'band is the peptide-binding domain, where a functional signal is expected: '
               f'mature residues {PBD["I"][0]}-{PBD["I"][1]} for class I (A, B, C) and '
               f'{PBD["II"][0]}-{PBD["II"][1]} for class II (D-genes); positions are '
               f'mature-protein numbered, so the leader peptide is negative and falls left of '
               f'the band. Genes are in chr6 order, taken from {order_src}. Full explanation, '
               f'symbol definitions and the estimator: {Path(args.out_png).stem}.md'),
        header=(args.cohort, f'HLA amino-acid omnibus — {mlab} model', ''),
        top_pad=0.50, hspace=0.62, wspace=0.24, left='auto', right=0.985,
        margin_axes=axes)

    # After caption_block, which resizes the canvas: S.gene_labels measures the
    # rendered axes box to size its strip and its leader arms.
    for gene, ax in zip(genes, axes):
        label_panel(ax, sels[gene], args)
    for i, (gene, ax) in enumerate(zip(genes, axes)):
        S.panel_tag(ax, string.ascii_lowercase[i], title=gene, pad=S.strip_pad(ax))
    fig.savefig(args.out_png)
    plt.close(fig)

    # ── Documents ────────────────────────────────────────────────────────────
    rows = []
    for gene in genes:
        d = om[om['gene'] == gene]
        good = d[d['usable']]
        n_res = int(((qc['gene'] == gene) & (qc['marker_class'] == 'residue')
                     & qc['tested']).sum()) if len(qc) else 0
        rows.append({'gene': gene, 'class': hla_class(gene) or '—',
                     'positions_tested': int(len(d)), 'usable_fits': int(len(good)),
                     'residue_markers_tested': n_res,
                     'significant': int((good['P'] < args.alpha).sum()) if len(good) else 0,
                     'min_P': float(good['P'].min()) if len(good) else None,
                     'strongest_position': (str(good.nsmallest(1, 'P').iloc[0]['position'])
                                            if len(good) else '—'),
                     'labelled': int(len(sels[gene]))})
    per_gene = pd.DataFrame(rows)
    status_counts = (om.loc[~om['usable'], 'status'].value_counts().rename_axis('status')
                     .reset_index(name='positions') if n_bad else None)

    values = [('cohort', args.cohort), ('model', args.model),
              ('genes drawn', len(genes)), ('grid', f'{nrows} x {ncols}'),
              ('positions in the omnibus table', int(len(om))),
              ('positions with a usable fit', n_ok),
              ('positions tested without a usable fit', n_bad),
              ('alpha', float(args.alpha)),
              ('positions clearing alpha', n_sig),
              ('of those, inside a peptide-binding domain', in_band),
              ('strongest position', str(top['position']) if top is not None else 'none'),
              ('strongest P', float(top['P']) if top is not None else None),
              ('strongest deviance', float(top['deviance']) if top is not None else None),
              ('strongest df', int(top['df']) if top is not None and pd.notna(top['df'])
               else None),
              ('gene order source', order_src),
              ('genes appended to --gene-order', ', '.join(appended) or 'none'),
              ('y limit', round(float(ylim[1]), 3))]
    figure_doc.write_stats(args.out_png, peak=f'{args.cohort}/{args.model}', values=values)

    figure_doc.write_doc(
        args.out_png,
        title=f'Amino-acid omnibus — {args.cohort}, HLA {mlab} model',
        question=('At which amino-acid POSITIONS does the residue content of the HLA molecule '
                  'associate with disease — as opposed to which single residue does — and do '
                  'those positions fall in the peptide-binding domain, where a functional '
                  'explanation exists?'),
        panels=[(string.ascii_lowercase[i], f'HLA-{g}', _panel_doc(g, om, args))
                for i, g in enumerate(genes)],
        interpretation=(
            f'The omnibus is a JOINT test over a position\'s residues, so it answers a '
            f'different question from the residue markers in the scan figure: a position can '
            f'be significant here with no single residue significant there, when several '
            f'residues each carry part of the signal. {n_sig} position(s) clear '
            f'{S.sci(args.alpha)}'
            + (f', {in_band} of them inside a peptide-binding domain. '
               if n_sig else '. ')
            + 'The band is a PRIOR EXPECTATION drawn on the figure, not a filter and not a '
              'result: a hit inside it has a mechanism available to it, a hit outside it is '
              'not thereby refuted, and neither placement is evidence about the other. '
            + f'{n_bad} position(s) were tested without a usable fit and carry no P at all; '
              f'they are the rug ticks on the axis floor. Those are not null results — the m-1 '
              f'df test could not be evaluated there — so any count of "positions tested" that '
              f'includes them overstates what the figure actually measured. The positions '
              f'within a gene are also very far from independent: neighbouring residues sit on '
              f'the same haplotypes, so the number of hits in a gene is not a count of '
              f'independent signals.'),
        numbers=values,
        tables=[('Per gene', per_gene),
                ('Positions clearing alpha (smallest P first)',
                 sig[['position', 'gene', 'n_residues', 'df', 'n', 'reference_residue',
                      'determined_rate', 'deviance', 'P']] if len(sig) else None, 20),
                ('Positions tested without a usable fit', status_counts)],
        reading=[
            'Read the grid left to right, top to bottom: that is chr6 order, the same order '
            'the MHC scan figure lays out along its x-axis.',
            'Compare heights ACROSS panels freely — every panel is on one shared y limit, so a '
            'taller point is a stronger point wherever it sits.',
            'Read each x-axis as belonging to ITS OWN gene. IMGT numbering is per gene: '
            'position 57 of DQB1 and position 57 of A are different residues of different '
            'proteins, which is why there is no single shared x-axis here.',
            'Check whether a labelled position falls inside the shaded band. Inside, a '
            'peptide-binding mechanism is available; outside, the association still stands but '
            'has no binding-site explanation from this figure.',
            'Negative positions are the leader peptide, left of the band by construction, and '
            'are numbered that way by IMGT rather than being an error.',
            'A panel carrying only the note "no position eligible" means the gene had fewer '
            'than two residues surviving marker QC at any position, or a determination rate '
            'below the omnibus floor — not that the gene was left out.',
        ],
        limits=[
            'It does not say WHICH residue drives a significant position. The omnibus is a '
            'joint m-1 df test; the per-residue estimates are in the scan figure and the '
            'sumstats.',
            'It does not establish independence between positions, or between genes. Long-'
            'range LD across the MHC means a position can be significant purely by tagging a '
            'causal one elsewhere in the region; that is what the conditional rounds are for.',
            'The peptide-binding band is drawn from the canonical domain boundaries, not from '
            'a structure fitted to these data. It is context for the reader, never a test.',
            'A position with status != ok is not a null result. Its fit could not be '
            'evaluated, so the figure is silent about it beyond marking that it was tried.',
        ],
        defs=['model', 'gw_sig'],
        model=(r'$D = 2(\log L_{\mathrm{full}} - \log L_{\mathrm{null}}) \sim \chi^2_{m-1}$, '
               r'$m$ = residues at the position, one dropped as the reference level'),
        methods_ref='../docs/METHODS.md')

    print(f'[plot_omnibus] {args.cohort}/{args.model}: {len(genes)} gene(s) in a '
          f'{nrows}x{ncols} grid, {n_ok:,} usable fit(s), {n_bad:,} without, '
          f'{n_sig} clearing {S.sci(args.alpha)} -> {args.out_png}')


def _panel_doc(gene, om, args):
    """The sidecar's paragraph for one gene panel — the numbers are this rendering's."""
    d = om[om['gene'] == gene]
    good = d[d['usable']]
    cls = hla_class(gene)
    band = (f'The shaded band is the class {cls} peptide-binding domain, mature residues '
            f'{PBD[cls][0]}-{PBD[cls][1]}.' if cls else
            'No peptide-binding domain is shaded: the gene is not one of the classical class I '
            'or class II molecules this rule covers.')
    if not len(d):
        return (f'No position of HLA-{gene} was eligible for the omnibus, so the axes is empty '
                f'and says so. A position is eligible only when at least two of its residues '
                f'survive marker QC and its determination rate clears the floor. {band}')
    n_sig = int((good['P'] < args.alpha).sum()) if len(good) else 0
    best = good.nsmallest(1, 'P').iloc[0] if len(good) else None
    return (f'{len(d)} position(s) of HLA-{gene} entered the omnibus, {len(good)} with a usable '
            f'fit; x is the IMGT residue position, y the -log10 P of the m-1 df deviance test. '
            + (f'The strongest is {best["position"]} '
               f'(deviance {float(best["deviance"]):.2f} on {int(best["df"])} df, '
               f'P = {float(best["P"]):.3g}), and {n_sig} position(s) clear '
               f'{S.sci(args.alpha)}. ' if best is not None else
               'No position of this gene returned a usable fit. ')
            + band)


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:                       # noqa: BLE001
        print(f'Error in plot_omnibus: {e}', file=sys.stderr)
        sys.exit(1)
