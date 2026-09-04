#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The canonical MHC association figure for one cohort x one model --
#           the whole classical MHC on ONE chr6 axis, the two marker classes as
#           two tracks above the gene layout they sit in.
#
#           WHY THREE PANELS ON ONE SHARED AXIS. An HLA scan answers a
#           POSITIONAL question -- which GENE carries the signal -- and a
#           -log10 P track alone cannot answer it: 29.7-33.1 Mb is 3.4 Mb of
#           coordinates nobody reads HLA-DQB1 out of. Panel (c) is the key that
#           turns a peak position into a gene name, and it only works because it
#           sits on the SAME x-axis, to the pixel, as the two tracks above it.
#           Side by side, or on a separate figure, the reader has to re-locate
#           every peak by eye, which is the defect this layout exists to remove.
#
#           WHY THE TWO CLASSES ARE TWO PANELS AND NOT ONE PILE. An allele
#           marker and an amino-acid residue marker are different hypotheses
#           carrying different multiple-testing burdens -- hla_marker_qc.py
#           emits one --extract list per class precisely so the two scans stay
#           separate analyses. Drawn on one axis, a residue peak would appear
#           above an allele threshold it was never judged against. Each panel
#           therefore carries ITS OWN Bonferroni line, and the two shared
#           references (the effective-tests line and the genome-wide line the
#           rest of the paper uses) are drawn on both so the panels remain
#           comparable where comparison is legitimate.
#
#           COLOUR IS THE GENE. Points take an alternating two-tone by gene in
#           chr6 order (S.CHROM_BANDS), exactly as a Manhattan bands
#           chromosomes: adjacent genes separate without a legend, and every
#           point matches the colour of its box in (c). HLA genes are dense and
#           unequal in marker count, so a single colour would let HLA-B's
#           several hundred markers read as one continuous smear across HLA-C.
#
#           THE TYPING-ARTEFACT WARNING. A SIGNIFICANT ALLELE whose call the
#           reference allele panel never carries (marker QC in_reference == 0)
#           is ringed in S.ACCENT. hla.typing's confound figure showed the
#           controls carry roughly twice the share of such calls and that the
#           difference is differential by phenotype -- a common allele depleted
#           more in controls than in cases reads as case enrichment. A hit on
#           one of those markers has exactly the shape a typing artefact makes,
#           so the figure flags it rather than leaving the caveat to prose.
#           in_reference is an ALLELE-panel statement, so residue markers carry
#           no flag and are never ringed.
#
#           A MARKER POSITION IS A PLOTTING CONVENTION. build_hla_markers.py
#           places each marker at its gene's GRCh38 start plus its index within
#           the gene, so the x-axis is honest about WHICH GENE and says nothing
#           about which base. The caption states this; the figure must not be
#           read as fine mapping to a base pair.
# Component: assoc_hla
# Used by : assoc_hla.nf  process PLOT_HLA_SCAN
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt              # noqa: E402
import numpy as np                           # noqa: E402
import pandas as pd                          # noqa: E402
from matplotlib.patches import Rectangle     # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_shared' / 'scripts'))
import plot_style as S                       # noqa: E402
import figure_doc                            # noqa: E402

PLOT_H = 6.2            # inches of croppable plot block (titles + axes + x-label)
MB = 1.0e6

# The 16 canonical plink2 --glm columns plus the four hla_to_sumstats.py adds.
# Nothing outside this list is read, so a schema change is a loud KeyError here
# rather than a quietly mis-drawn panel.
SUMSTAT_COLS = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'A1', 'A1_FREQ', 'TEST', 'OBS_CT',
                'OR', 'LOG(OR)_SE', 'L95', 'U95', 'Z_STAT', 'P', 'ERRCODE',
                'gene', 'position', 'cohort', 'model']
QC_COLS = ['id', 'marker_class', 'gene', 'pos', 'mac', 'maf', 'in_reference', 'tested']


def parse_args():
    p = argparse.ArgumentParser(
        description='MHC association figure for one cohort x one model: allele track, '
                    'residue track and the gene layout, all on one chr6 axis.')
    p.add_argument('--sumstat-allele', required=True,
                   help='hla_to_sumstats.py output for --marker-class allele')
    p.add_argument('--sumstat-residue', required=True,
                   help='hla_to_sumstats.py output for --marker-class residue')
    p.add_argument('--marker-qc', required=True,
                   help='hla_marker_qc.py output; supplies gene extents and in_reference')
    p.add_argument('--cohort', required=True)
    p.add_argument('--model', required=True, help='fixed | random')
    p.add_argument('--model-label', default='',
                   help='human-readable model name for the title; defaults to --model')
    p.add_argument('--alpha-bonferroni-allele', type=float, required=True,
                   help='0.05 / (number of allele markers tested) — panel (a) decision line')
    p.add_argument('--alpha-bonferroni-residue', type=float, required=True,
                   help='0.05 / (number of residue markers tested) — panel (b) decision line')
    p.add_argument('--alpha-effective', type=float, default=None,
                   help='threshold from the EFFECTIVE number of independent tests; omitted '
                        'when not given rather than guessed')
    p.add_argument('--alpha-genomewide', type=float, default=S.GW_ALPHA,
                   help='the genome-wide line the rest of the study uses, for reference')
    p.add_argument('--label-top', type=int, default=8,
                   help='maximum markers named per track, smallest P first')
    p.add_argument('--anno-style', default='auto', choices=sorted(S.ANNO_STYLES),
                   help='gwaslab label style passed to plot_style.gene_labels')
    p.add_argument('--repel-force', type=float, default=S.REPEL_FORCE,
                   help='minimum label separation as a fraction of the plotted span')
    p.add_argument('--x-min-mb', type=float, default=29.7,
                   help='left edge of the shared chr6 axis, in Mb; widened if data fall outside')
    p.add_argument('--x-max-mb', type=float, default=33.1,
                   help='right edge of the shared chr6 axis, in Mb; widened if data fall outside')
    p.add_argument('--gene-box-min-kb', type=float, default=12.0,
                   help='minimum DRAWN width of a gene box in (c). A class I gene is 3-5 kb '
                        'and invisible on a 3.4 Mb axis; the box is centred on the true '
                        'extent and the exaggeration is stated in the caption.')
    p.add_argument('--out-png', required=True)
    return p.parse_args()


# ── Loading ──────────────────────────────────────────────────────────────────
def load_sumstat(path, cls, args):
    """One marker class's summary statistics, restricted to usable fits.

    The usability rule is hla_to_sumstats.py's own: it writes <ENGINE>_NO_FIT
    into ERRCODE for a non-finite OR/SE or an out-of-range P. Re-checking P here
    as well keeps the figure honest if a future engine leaves ERRCODE '.' on a
    degenerate fit, which is exactly what plink2 does genome-wide.
    """
    head = pd.read_csv(path, sep='\t', nrows=0)
    cols = [c for c in SUMSTAT_COLS if c in head.columns]
    missing = [c for c in ('POS', 'ID', 'P') if c not in cols]
    if missing:
        raise SystemExit(f'ABORT: {path} lacks {missing}; it is not a hla_to_sumstats.py file')
    dtype = {c: str for c in ('#CHROM', 'ID', 'ERRCODE', 'gene', 'position', 'cohort',
                              'model', 'REF', 'ALT', 'A1', 'TEST') if c in cols}
    d = pd.read_csv(path, sep='\t', usecols=cols, dtype=dtype)
    n_rows = len(d)
    for c in ('POS', 'P', 'OR', 'A1_FREQ', 'OBS_CT', 'L95', 'U95', 'Z_STAT', 'LOG(OR)_SE'):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors='coerce')
    for c in ('gene', 'position'):
        d[c] = d[c].fillna('') if c in d.columns else ''
    # The file is one cohort x one model by construction. Disagreeing labels mean
    # the wrong file was wired in, and a figure captioned with the arguments would
    # be a caption that lies — say so on stderr rather than silently relabelling.
    for col, want in (('cohort', args.cohort), ('model', args.model)):
        if col in d.columns and len(d) and set(d[col].dropna().unique()) - {want}:
            print(f'[plot_hla_scan] WARNING: {path} carries {col} '
                  f'{sorted(set(d[col].dropna().unique()))}, figure is labelled {want!r}',
                  file=sys.stderr)
    ok = d['ERRCODE'].fillna('.').eq('.') if 'ERRCODE' in d.columns \
        else pd.Series(True, index=d.index)
    ok &= d['P'].notna() & (d['P'] > 0) & (d['P'] <= 1) & d['POS'].notna()
    d = d[ok].copy()
    d['marker_class'] = cls
    d['x'] = d['POS'].to_numpy(dtype=float) / MB
    d['y'] = -np.log10(d['P'].to_numpy(dtype=float))
    return d, n_rows, n_rows - len(d)


def load_qc(path):
    """Marker QC, restricted to the columns this figure is allowed to read."""
    d = pd.read_csv(path, sep='\t', dtype={'id': str, 'marker_class': str, 'gene': str})
    for c in ('id', 'marker_class', 'gene', 'pos'):
        if c not in d.columns:
            raise SystemExit(f'ABORT: {path} has no {c!r} column')
    d = d[[c for c in QC_COLS if c in d.columns]].copy()
    for c in ('pos', 'mac', 'maf', 'in_reference'):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors='coerce')
    # pandas writes booleans as True/False text; read them back as booleans rather
    # than trusting the dtype inference of a column that may be all-True.
    d['tested'] = (d['tested'].astype(str).str.lower().isin(('true', '1'))
                   if 'tested' in d.columns else True)
    if 'in_reference' not in d.columns:
        d['in_reference'] = np.nan
    return d


def marker_label(mid, position):
    """What a marker is CALLED in a label band.

    'HLA_A*24:02' -> 'A*24:02'; 'AA_DQB1_57_D' -> 'DQB1:57D' (IMGT positions can
    be negative, e.g. 'AA_A_-22_I' -> 'A:-22I'). The 'HLA_'/'AA_' prefixes are
    id-space bookkeeping and cost a third of the label's width for nothing.
    """
    mid = str(mid)
    if mid.startswith('HLA_'):
        return mid[4:]
    if mid.startswith('AA_'):
        parts = mid.split('_')
        if len(parts) >= 4:
            return f'{parts[1]}:{parts[2]}{parts[-1]}'
    return str(position) if position else mid


def gene_table(qc, scans, args):
    """One row per gene: its marker extent on chr6, its two-tone colour, its counts.

    Extents come from the MARKER POSITIONS, not from a gene annotation: the
    component has no coordinate file at this stage, and the markers were placed
    at gene_start + index by build_hla_markers.py, so min(pos) is the gene start
    and the span is the gene's marker count in bases. That is the honest extent
    for THIS figure — it is where the plotted points actually are.
    """
    src = qc[['gene', 'pos']].dropna() if len(qc) else pd.DataFrame(columns=['gene', 'pos'])
    if not len(src):
        src = pd.concat([s[['gene', 'POS']].rename(columns={'POS': 'pos'}) for s in scans
                         if len(s)], ignore_index=True) if any(len(s) for s in scans) \
            else pd.DataFrame(columns=['gene', 'pos'])
    src = src[src['gene'].astype(str).str.len() > 0]
    if not len(src):
        return pd.DataFrame(columns=['gene', 'lo', 'hi', 'mid', 'colour'])
    g = (src.groupby('gene')['pos'].agg(['min', 'max']).rename(columns={'min': 'lo', 'max': 'hi'})
         .sort_values('lo').reset_index())
    g['lo'] /= MB
    g['hi'] /= MB
    g['mid'] = (g['lo'] + g['hi']) / 2.0
    # Alternating two-tone in chr6 order — a Manhattan's chromosome banding, at
    # gene scale. Two tones and not N hues: the reader needs adjacent genes to
    # SEPARATE, and (c) already names them, so a per-gene hue would be a second
    # encoding of information the figure already carries.
    g['colour'] = [S.CHROM_BANDS[i % len(S.CHROM_BANDS)] for i in range(len(g))]
    for cls, s in (('allele', scans[0]), ('residue', scans[1])):
        n = s.groupby('gene').size() if len(s) else pd.Series(dtype=int)
        g[f'n_{cls}'] = g['gene'].map(n).fillna(0).astype(int)
    return g


# ── Thresholds ───────────────────────────────────────────────────────────────
def threshold_set(args, alpha_main, main_name):
    """The lines one track is judged against, strongest rule last.

    Rule 7 of the shared figure contract: every -log10 P axis carries BOTH tiers.
    Here the tiers are this class's own Bonferroni line (the decision), the
    effective-tests line (a less conservative reference, drawn only when the
    caller supplies one) and the study-wide genome-wide line.
    """
    out = [(main_name, alpha_main, S.ACCENT, '-', 1.2)]
    if args.alpha_effective:
        out.append(('effective', args.alpha_effective, S.REFERENCE, '--', 0.9))
    if args.alpha_genomewide:
        out.append(('genome-wide', args.alpha_genomewide, S.NEUTRAL_D, ':', 1.0))
    return [t for t in out if t[1] and 0 < t[1] <= 1]


def draw_thresholds(ax, entries):
    for _name, alpha, colour, ls, lw in entries:
        ax.axhline(-np.log10(alpha), color=colour, lw=lw, ls=ls, zorder=3)


def label_thresholds(ax, entries):
    """Name each line at the right edge, separated BY MEASUREMENT.

    Bonferroni and effective-test thresholds routinely land within a line height
    of each other (0.05/2000 and 0.05/900 differ by 0.35 on this axis), so a
    label drawn at each line's own y would print one on top of the other. The
    positions are repelled in AXES FRACTION with plot_style's own 1-D repel, on a
    gap measured from the rendered axes height — not on a tuned offset.

    Must run AFTER caption_block: the axes box is resized by it, so a gap
    computed before it is a gap in the wrong units.
    """
    if not entries:
        return
    fig = ax.figure
    lo, hi = ax.get_ylim()
    if hi <= lo:
        return
    fs = plt.rcParams['legend.fontsize'] - 0.5
    try:
        h_in = ax.get_window_extent().height / fig.dpi
    except Exception:
        h_in = 1.5
    gap = min(0.5, (fs * 1.55 / 72.0) / max(h_in, 0.25))
    ys = np.array([(-np.log10(a) - lo) / (hi - lo) for _n, a, _c, _l, _w in entries])
    ys = S.repel_from_centre(ys, gap, gap / 2.0, 1.0 - gap / 2.0)
    box = dict(boxstyle='square,pad=0.14', facecolor='white', edgecolor='none', alpha=0.88)
    for (name, alpha, colour, _ls, _lw), y in zip(entries, ys):
        ax.text(0.999, y, f'{name} {S.mathsci(alpha)}', transform=ax.transAxes,
                ha='right', va='center', color=colour, fontsize=fs, zorder=8,
                bbox=box, clip_on=False)


# ── Panels ───────────────────────────────────────────────────────────────────
def scan_panel(ax, d, genes, entries, alpha_main, ylim, args, cls):
    """One marker class: -log10 P against chr6 position, coloured by gene.

    Returns the rows drawn with a typing-artefact ring, so the caption can count
    them without recomputing the rule.
    """
    S.despine(ax, grid_axis='y')
    ax.set_ylabel(S.NEGLOG10P)
    ax.set_ylim(*ylim)
    if not len(d):
        ax.text(0.5, 0.5, f'no {cls} marker survived QC in this cohort',
                transform=ax.transAxes, ha='center', va='center', color=S.REFERENCE,
                fontsize=plt.rcParams['legend.fontsize'])
        draw_thresholds(ax, entries)
        return d
    cmap = dict(zip(genes['gene'], genes['colour'])) if len(genes) else {}
    colours = [cmap.get(str(g), S.CHROM_BANDS[0]) for g in d['gene']]
    # Rasterised, as every dense scatter in this project is: the residue track is
    # several thousand points and the axes, ticks and text must stay vector.
    ax.scatter(d['x'].to_numpy(), d['y'].to_numpy(), s=9.0 if cls == 'allele' else 6.0,
               c=colours, linewidths=0, rasterized=True, zorder=4)
    draw_thresholds(ax, entries)

    art = d[(d['P'] < alpha_main) & (d['in_reference'] == 0)]
    if len(art):
        ax.scatter(art['x'].to_numpy(), art['y'].to_numpy(), s=46, facecolors='none',
                   edgecolors=S.ACCENT, linewidths=1.1, zorder=7)
    # NO in-panel legend for the ring, deliberately. Every free corner of a
    # -log10 P track is either data or the leader arm of a named marker, so a
    # legend box here lands on something whichever corner is chosen; the
    # Manhattan of the sibling component carries none for the same reason. The
    # ring is defined in the caption's (a) line and in the sidecar, which is
    # where this project defines its symbols.
    return art


def label_track(ax, d, alpha_main, args):
    """Name the top hits of one track, in a strip ABOVE the data area.

    Significant markers first, capped at --label-top; if nothing is significant
    the strongest three are still named, because a reader has to be able to see
    WHERE the maximum of a null-looking track sits without opening the TSV.
    MUST run after caption_block — S.gene_labels measures the rendered axes box.
    """
    if not len(d):
        return pd.DataFrame(columns=d.columns)
    sig = d[d['P'] < alpha_main]
    sel = (sig.nsmallest(max(0, args.label_top), 'P') if len(sig)
           else d.nsmallest(min(3, args.label_top, len(d)), 'P'))
    if not len(sel):
        return sel
    sel = sel.sort_values('x')
    is_sig = (sel['P'] < alpha_main).to_numpy()
    base = plt.rcParams['legend.fontsize']
    # One band, tier carried by ink alone (the shared contract's rule): bold
    # ACCENT for a marker that cleared its line, regular REFERENCE for one that
    # is only the local maximum.
    S.gene_labels(ax, sel['x'].to_numpy(), sel['y'].to_numpy(),
                  [marker_label(r.ID, r.position) for r in sel.itertuples()],
                  anno_style=args.anno_style, repel_force=args.repel_force,
                  color=[S.ACCENT if s else S.REFERENCE for s in is_sig],
                  weight=['bold' if s else 'normal' for s in is_sig],
                  fontsize=[base if s else base - 1.0 for s in is_sig])
    return sel


def gene_track(ax, genes, args):
    """The gene layout: one labelled box per HLA gene, at its real coordinates.

    This is the panel that makes the figure readable. It is deliberately thin —
    it carries no measurement, only the mapping from position to gene name — and
    it shares the x-axis with both tracks so a peak can be read straight down
    onto the gene it sits in.

    THE BOXES ARE WIDENED TO BE VISIBLE. A class I gene spans 3-5 kb, which is
    under a thousandth of this axis and would render as nothing at all. Each box
    is CENTRED ON ITS TRUE EXTENT and drawn at least --gene-box-min-kb wide; the
    exaggeration is stated in the caption, because an invisible box and a lying
    box are both worse than a stated convention.
    """
    S.despine(ax, grid_axis='none')
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_xlabel(f'chr6 {S.POS_MB}')
    if not len(genes):
        ax.text(0.5, 0.5, 'no gene extent available', transform=ax.transAxes,
                ha='center', va='center', color=S.REFERENCE,
                fontsize=plt.rcParams['legend.fontsize'])
        return
    min_w = max(args.gene_box_min_kb, 0.0) / 1000.0        # kb -> Mb
    for r in genes.itertuples():
        w = max(float(r.hi) - float(r.lo), min_w)
        ax.add_patch(Rectangle((float(r.mid) - w / 2.0, 0.58), w, 0.34,
                               facecolor=r.colour, edgecolor='white', linewidth=0.5,
                               zorder=4, clip_on=True))


def _cluster_spread(xs, gap, lo, hi):
    """Place gene names at `gap` apart with the SMALLEST total displacement.

    plot_style.repel_from_centre is the right tool for a label band, where every
    label is displaced anyway and only the order has to survive. It is the wrong
    tool here: it repels symmetrically, so spreading the five class II names at
    32.5-33.1 Mb drags HLA-A at 29.9 Mb — which has 1.3 Mb of empty axis around
    it — a full label width off its own box. On a key panel, a name that is not
    over the box it names is the one defect that cannot be tolerated.

    The standard block-merge instead: walk left to right, and whenever the next
    name would land inside its neighbour, merge the two into one run and place
    the run at the least-squares optimum of its members' wanted positions. A run
    absorbs a neighbour only when it actually has to, so an isolated gene never
    moves at all, and a run that grows is centred on its own members rather than
    pushed off one end.

    The clamp is then the same min/max sweep plot_style uses: pinning the last
    run to the right edge and pulling back with min() moves only the names the
    edge really binds, where translating the whole arrangement would drag every
    isolated name with it.
    """
    xs = np.asarray(xs, dtype=float)
    n = len(xs)
    if not n:
        return xs
    if gap * (n - 1) > (hi - lo):
        # Genuinely more names than the axis can hold at a readable size. Even
        # spacing is the least-bad answer; the caller should drop genes instead.
        out = np.linspace(lo, hi, n)
        res = np.empty_like(xs)
        res[np.argsort(xs)] = out
        return res
    order = np.argsort(xs)
    want = xs[order]
    blocks = []                       # [count, start] per run, left to right
    for p_i in want:
        blocks.append([1, float(p_i)])
        while len(blocks) >= 2:
            (n1, s1), (n2, s2) = blocks[-2], blocks[-1]
            if s2 >= s1 + gap * n1 - 1e-12:          # clears the run before it
                break
            # Least-squares start for the merged run: each member wants
            # start = its own position minus its offset within the run.
            s = (s1 * n1 + (s2 - gap * n1) * n2) / (n1 + n2)
            blocks[-2:] = [[n1 + n2, s]]
    q = np.concatenate([[s + gap * k for k in range(cnt)] for cnt, s in blocks])
    if q[-1] > hi:                    # pull back from the right, minimally
        q[-1] = hi
        for i in range(n - 2, -1, -1):
            q[i] = min(q[i], q[i + 1] - gap)
    if q[0] < lo:                     # and from the left
        q[0] = lo
        for i in range(1, n):
            q[i] = max(q[i], q[i - 1] + gap)
    res = np.empty_like(xs)
    res[order] = q
    return res


def gene_track_labels(ax, genes, fontsize=None):
    """Gene names UNDER their boxes, separated on measured text width.

    HLA-B and HLA-C are 80 kb apart and DQA1/DQB1 closer still; at 3.4 Mb across
    a 6.5 in axis their names overlap outright. The gap is measured from the
    rendered strings rather than estimated, and each displaced name keeps an
    L-shaped leader — vertical down from its own box, then horizontal to the
    text — so the displacement never costs the reader the anchor. A name that
    did not have to move gets no leader at all, because a leader that points
    straight down at the box it is already under is ink for nothing.

    Must run AFTER caption_block, which resizes the figure.
    """
    if not len(genes):
        return
    fig = ax.figure
    fs = fontsize or (plt.rcParams['legend.fontsize'] - 0.5)
    names = [str(g) for g in genes['gene']]
    try:
        w_in = ax.get_window_extent().width / fig.dpi
    except Exception:
        w_in = S.COL_DOUBLE * 0.9
    gap_in = max(_text_width_in(fig, n, fs) for n in names) + 0.07
    x0, x1 = ax.get_xlim()
    per_in = (x1 - x0) / w_in if w_in else 1.0
    gap = gap_in * per_in
    xs = genes['mid'].to_numpy(dtype=float)
    lab = _cluster_spread(xs, gap, x0 + gap / 2.0, x1 - gap / 2.0)
    for mid, lx, name, col in zip(xs, lab, names, genes['colour']):
        if abs(lx - mid) > 0.002 * (x1 - x0):
            ax.plot([mid, mid, lx], [0.56, 0.46, 0.46], color=S.NEUTRAL_D, lw=0.6,
                    solid_joinstyle='miter', zorder=3, clip_on=False)
        else:
            ax.plot([mid, mid], [0.56, 0.50], color=S.NEUTRAL_D, lw=0.6, zorder=3,
                    clip_on=False)
        ax.text(lx, 0.42, name, ha='center', va='top', fontsize=fs, fontstyle='italic',
                color=S.INK if col == S.CHROM_BANDS[0] else S.INK_SOFT,
                zorder=5, clip_on=False)


def _text_width_in(fig, s, fontsize):
    """Rendered width of a string in inches, without leaving it on the canvas.

    The gene names are repelled against each other, so the gap has to come from
    the real advance width; a chars-per-point estimate understates 'DPB1' and
    overstates 'A' and produces a band that is both crowded and sparse.
    """
    try:
        fig.canvas.draw()
        t = fig.text(0, 0, s, fontsize=fontsize, fontstyle='italic')
        w = t.get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.dpi
        t.remove()
        return w
    except Exception:
        return 0.62 * len(s) * fontsize / 72.0


# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    S.setup_style('paper')
    mlab = args.model_label or args.model

    qc = load_qc(args.marker_qc)
    al, n_al_rows, n_al_bad = load_sumstat(args.sumstat_allele, 'allele', args)
    re_, n_re_rows, n_re_bad = load_sumstat(args.sumstat_residue, 'residue', args)
    ref = dict(zip(qc['id'], qc['in_reference']))
    for d in (al, re_):
        d['in_reference'] = d['ID'].map(ref) if len(d) else pd.Series(dtype=float)

    genes = gene_table(qc, (al, re_), args)
    ent_a = threshold_set(args, args.alpha_bonferroni_allele, 'Bonferroni')
    ent_r = threshold_set(args, args.alpha_bonferroni_residue, 'Bonferroni')

    # ONE y limit for both tracks. They are the same quantity on the same axis, so
    # a taller peak has to be a stronger peak; per-panel limits would rescale the
    # residue track and invite exactly the comparison the reader should not make.
    alphas = [a for _n, a, _c, _l, _w in ent_a + ent_r]
    ymax = max([float(d['y'].max()) for d in (al, re_) if len(d)] or [0.0])
    ytop = max(ymax * 1.12, -np.log10(min(alphas)) + 0.7 if alphas else 0.0, 4.0)
    ylim = (0.0, ytop)

    xs = [d['x'] for d in (al, re_) if len(d)]
    xlo = min([args.x_min_mb] + [float(s.min()) - 0.05 for s in xs])
    xhi = max([args.x_max_mb] + [float(s.max()) + 0.05 for s in xs])

    fig = plt.figure(figsize=(S.COL_DOUBLE, PLOT_H))
    # Two equal tracks over a thin key. (a) and (b) get the same share because
    # both give up a label strip; (c) carries no data and needs only its boxes,
    # its names and the shared x tick labels.
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 0.34])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0], sharex=ax_a)
    ax_c = fig.add_subplot(gs[2, 0], sharex=ax_a)

    art_a = scan_panel(ax_a, al, genes, ent_a, args.alpha_bonferroni_allele, ylim, args, 'allele')
    art_b = scan_panel(ax_b, re_, genes, ent_r, args.alpha_bonferroni_residue, ylim, args,
                       'residue')
    gene_track(ax_c, genes, args)
    ax_a.set_xlim(xlo, xhi)
    # sharex keeps the three panels aligned to the pixel; only (c) prints the
    # numbers, so the axis is stated once instead of three times.
    for ax in (ax_a, ax_b):
        ax.tick_params(labelbottom=False)

    n_sig_a = int((al['P'] < args.alpha_bonferroni_allele).sum()) if len(al) else 0
    n_sig_r = int((re_['P'] < args.alpha_bonferroni_residue).sum()) if len(re_) else 0
    top = None
    for d, alpha in ((al, args.alpha_bonferroni_allele), (re_, args.alpha_bonferroni_residue)):
        if not len(d):
            continue
        r = d.loc[d['P'].idxmin()]
        if top is None or float(r['P']) < float(top['P']):
            top = r
    top_txt = ('none' if top is None else
               f'{marker_label(top["ID"], top["position"])} ({top["marker_class"]}, '
               f'{top["gene"] or "?"}, P = {S.p_tex(float(top["P"]))})')

    n_qc_a = int(((qc['marker_class'] == 'allele') & qc['tested']).sum()) if len(qc) else 0
    n_qc_r = int(((qc['marker_class'] == 'residue') & qc['tested']).sum()) if len(qc) else 0
    n_art = len(art_a) + len(art_b)
    thr_txt = (f'Bonferroni {S.sci(args.alpha_bonferroni_allele)} (allele) and '
               f'{S.sci(args.alpha_bonferroni_residue)} (residue)'
               + (f', effective-tests {S.sci(args.alpha_effective)}' if args.alpha_effective else '')
               + (f', genome-wide {S.sci(args.alpha_genomewide)}' if args.alpha_genomewide else ''))
    hit_txt = (f'{n_sig_a} allele and {n_sig_r} residue marker'
               f'{"" if n_sig_a + n_sig_r == 1 else "s"} clear their own Bonferroni threshold'
               if n_sig_a + n_sig_r else 'no marker clears its own Bonferroni threshold')

    S.caption_block(
        fig, plot_h=PLOT_H,
        title=(f'{args.cohort} under the HLA {mlab} model: {hit_txt}; strongest marker '
               f'{top_txt}.'),
        panels=[
            f'HLA allele markers, {S.NEGLOG10P} against chr6 position, coloured by gene in '
            f'alternating tones; an open red ring marks a significant allele the reference '
            f'allele panel never carries — a typing-artefact warning, not a result.',
            'Amino-acid residue markers on the same axis and the same '
            f'{S.NEGLOG10P} scale, so a taller point is a stronger point in both tracks.',
            'Gene layout: one box per HLA gene at its marker extent, which is what identifies '
            'the gene a peak above it sits in.',
        ],
        notes=(f'{args.cohort}, {mlab} model; {len(al):,} allele and {len(re_):,} residue '
               f'markers drawn of {n_qc_a:,} / {n_qc_r:,} passing marker QC '
               f'({n_al_bad:,} + {n_re_bad:,} excluded on ERRCODE or an unusable P). '
               f'Thresholds: {thr_txt}. Marker positions are a plotting convention — gene '
               f'start plus index within the gene — so the axis places a marker in the right '
               f'GENE and asserts nothing about the base. Gene boxes in (c) are centred on '
               f'their true extent and drawn at least {args.gene_box_min_kb:.0f} kb wide to be '
               f'visible at this scale. Full explanation, symbol definitions and the '
               f'estimator: {Path(args.out_png).stem}.md'),
        # The subject line is a caption_block HEADER, not a suptitle. It has to be
        # ABOVE the crop line — a cropped publication figure that no longer says
        # which cohort and model it is has lost the one thing the pixels cannot
        # carry — and the header row ADDS its height to the block, so it cannot
        # collide with the panel letters the way a suptitle at y≈1 does once a
        # label strip is reserved above panel (a).
        header=(args.cohort, f'HLA {mlab} model', ''),
        # top_pad holds the label strip AND the letter/title above (a); right
        # leaves canvas outside the axes, because gene_labels anchors its text
        # ha='left' and the rightmost name overhangs the last data position.
        top_pad=0.55, hspace=0.30, left='auto', right=0.955,
        margin_axes=[ax_a, ax_b, ax_c])

    # Everything below MEASURES the rendered figure, so all of it runs after
    # caption_block, which resizes the canvas and re-runs subplots_adjust.
    label_thresholds(ax_a, ent_a)
    label_thresholds(ax_b, ent_r)
    sel_a = label_track(ax_a, al, args.alpha_bonferroni_allele, args)
    sel_b = label_track(ax_b, re_, args.alpha_bonferroni_residue, args)
    # The two tracks share a y limit, and a shared limit on unequal axes heights
    # is not a shared scale — the label strips come out of (a) and (b)
    # independently, so bring the taller one down to match.
    S.equalise_row_heights(ax_a, ax_b)
    gene_track_labels(ax_c, genes)
    # Centred titles stay VERY short: they share the title slot with the
    # flush-left bold letter, and a long one walks straight into it.
    for ax, letter, ttl in ((ax_a, 'a', 'alleles'), (ax_b, 'b', 'residues'),
                            (ax_c, 'c', 'genes')):
        S.panel_tag(ax, letter, title=ttl, pad=S.strip_pad(ax))
    fig.savefig(args.out_png)
    plt.close(fig)

    # ── Documents ────────────────────────────────────────────────────────────
    def hit_rows(d, alpha, n=12):
        if not len(d):
            return None
        cols = [c for c in ('ID', 'gene', 'position', 'A1_FREQ', 'OBS_CT', 'OR', 'L95', 'U95',
                            'P', 'in_reference') if c in d.columns]
        return d.nsmallest(min(n, len(d)), 'P')[cols]

    values = [('cohort', args.cohort), ('model', args.model),
              ('allele markers drawn', len(al)), ('residue markers drawn', len(re_)),
              ('allele markers excluded', n_al_bad), ('residue markers excluded', n_re_bad),
              ('allele markers passing QC', n_qc_a), ('residue markers passing QC', n_qc_r),
              ('genes drawn', int(len(genes))),
              ('alpha Bonferroni allele', float(args.alpha_bonferroni_allele)),
              ('alpha Bonferroni residue', float(args.alpha_bonferroni_residue)),
              ('alpha effective', float(args.alpha_effective) if args.alpha_effective else None),
              ('alpha genome-wide', float(args.alpha_genomewide) if args.alpha_genomewide else None),
              ('allele markers significant', n_sig_a),
              ('residue markers significant', n_sig_r),
              ('significant alleles absent from the reference panel', n_art),
              ('markers named in (a)', int(len(sel_a))),
              ('markers named in (b)', int(len(sel_b))),
              ('strongest marker', str(top['ID']) if top is not None else 'none'),
              ('smallest P', float(top['P']) if top is not None else None),
              ('axis (Mb)', f'{xlo:.2f}–{xhi:.2f}')]
    figure_doc.write_stats(args.out_png, peak=f'{args.cohort}/{args.model}', values=values)

    figure_doc.write_doc(
        args.out_png,
        title=f'MHC association scan — {args.cohort}, HLA {mlab} model',
        question=('Where in the classical MHC does this cohort show association under this '
                  'model, in which GENE does it sit, and is any of it the shape a typing '
                  'artefact makes rather than a result?'),
        panels=[
            ('a', 'HLA allele markers',
             'One point per tested HLA allele marker, x its marker position on chr6, '
             'y = -log10 P from the association engine. Points are coloured by gene in '
             'alternating tones, in chr6 order, so adjacent genes separate without a legend '
             'and each point matches its box in (c). The solid red line is this class\'s own '
             f'Bonferroni threshold ({S.sci(args.alpha_bonferroni_allele)} = 0.05 over the '
             'allele markers tested), which is the decision rule for this panel. An open red '
             'ring marks a significant allele with in_reference == 0 in the marker QC — a '
             'two-field call the reference allele panel never carries.'),
            ('b', 'Amino-acid residue markers',
             'The same construction over the residue markers, on the SAME y limit as (a) so a '
             'taller point is a stronger point in either track, and against its OWN Bonferroni '
             f'line ({S.sci(args.alpha_bonferroni_residue)}). The two classes are separate '
             'analyses with separate multiple-testing burdens — that is why they are separate '
             'panels rather than one scatter. Residue markers carry no in_reference flag: the '
             'reference panel is an allele panel, so "not in it" is not a statement about a '
             'residue, and none is ringed.'),
            ('c', 'Gene layout',
             'One box per HLA gene, spanning the chr6 extent of that gene\'s markers, in the '
             'same two-tone colour the points above use. A name sits directly under its own '
             'box unless a neighbour is closer than one label width, in which case that run '
             'of names — and only that run — is spread at exactly one label width and each '
             'displaced name keeps an L-shaped leader back to its box. The boxes are centred '
             'on their true extent and drawn at least '
             f'{args.gene_box_min_kb:.0f} kb wide: a class I gene is 3-5 kb, under a '
             'thousandth of this axis, and would otherwise be invisible.'),
        ],
        interpretation=(
            f'{n_sig_a} allele and {n_sig_r} residue marker(s) clear their own Bonferroni '
            f'threshold. The two counts are NOT additive evidence: an allele marker and the '
            f'residue markers it carries are the same chromosomes counted twice, so a peak '
            f'appearing in both tracks at the same gene is one signal seen two ways, not two. '
            f'Reading the position off the x-axis alone is what (c) exists to prevent — the '
            f'marker coordinate is gene start plus index, so it identifies a GENE and no more. '
            + (f'{n_art} significant allele marker(s) are absent from the reference allele '
               f'panel and are ringed. hla.typing found that the controls carry roughly twice '
               f'the share of such calls as the cases, differentially by phenotype, so a hit '
               f'on one of them has the shape a typing artefact makes and cannot be reported '
               f'as a finding without that check.'
               if n_art else
               'No significant allele marker is absent from the reference allele panel, which '
               'is the outcome that leaves the typing-artefact route open but unexercised.')),
        numbers=values,
        tables=[('Strongest allele markers', hit_rows(al, args.alpha_bonferroni_allele), 12),
                ('Strongest residue markers', hit_rows(re_, args.alpha_bonferroni_residue), 12)],
        reading=[
            'Read (a) and (b) against their OWN red line. The two classes carry different '
            'numbers of tests, so the same P is significant in one panel and not in the other; '
            'that is the design, not an inconsistency.',
            'Take any peak and read straight down into (c). The gene box under it names the '
            'gene the signal is in. Nothing finer than the gene is claimed by this figure.',
            'Compare heights ACROSS the two panels freely — they share one y limit and, after '
            'the label strips are taken out, one axes height.',
            'Check whether a significant allele carries a red ring before quoting it. A ringed '
            'marker needs the typing check in hla.typing before it is a result.',
            'The dotted genome-wide line is a reference to the rest of the study, not this '
            'figure\'s decision rule. A marker between the Bonferroni line and it has cleared '
            'the burden that applies to it.',
        ],
        limits=[
            'It does not fine-map to a base. Marker positions are gene start plus index within '
            'the gene, so the axis resolves genes, not variants.',
            'It cannot separate one signal from another inside the MHC. Long-range LD across '
            'the region means a peak at one gene may be tagging a causal allele at a '
            'neighbouring one; that is what the conditional rounds and the omnibus figure are '
            'for.',
            'An allele absent from the reference panel is not thereby a wrong call — a small '
            'panel cannot sample a rare allele. The ring is a flag for follow-up, never a '
            'verdict.',
            'The two tracks are not independent evidence about each other, so no combined '
            'count of "significant markers" is meaningful.',
        ],
        defs=['model', 'or', 'gw_sig', 'errcode'],
        methods_ref='../docs/METHODS.md')

    print(f'[plot_hla_scan] {args.cohort}/{args.model}: allele {len(al):,} drawn '
          f'({n_sig_a} sig), residue {len(re_):,} drawn ({n_sig_r} sig), '
          f'{n_art} artefact-flagged, {len(genes)} gene(s) -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:                       # noqa: BLE001
        print(f'Error in plot_hla_scan: {e}', file=sys.stderr)
        sys.exit(1)
