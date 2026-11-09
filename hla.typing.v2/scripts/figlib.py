"""Figures of hla.typing.v2, drawn from the result tables.

Style: the project's shared plot_style (semantic palette: DATA = the measurement,
ACCENT = what must not be missed, REFERENCE = what it is judged against, NEUTRAL =
structure). Every figure is authored at a fixed width; `scale` is 'paper' for the
component's figures and 'slide' for the deck.
"""
import string
import sys
from pathlib import Path

import matplotlib
import matplotlib.font_manager

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402

sys.path.insert(0, str(C.PROJECT / 'analysis' / '_shared' / 'scripts'))
import plot_style as PS  # noqa: E402

LETTERS = string.ascii_lowercase
GROUPS = [('Class I', ['A', 'B', 'C', 'E', 'F', 'G']),
          ('Class II', ['DRA', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1', 'DMA', 'DMB', 'DOA', 'DOB']),
          ('DRB paralogues', ['DRB3', 'DRB4', 'DRB5'])]
ORDER = [g for _, gs in GROUPS for g in gs]

# categorical scales derived from the shared palette
STATE_COL = {'heterozygous': PS.DATA, 'homozygous': '#8FB4D4', 'present': PS.DATA_DARK,
             'absent': PS.NEUTRAL, 'ambiguous': PS.ACCENT_LT, 'no call': PS.ACCENT}
BASIS_COL = {'agreed by reads': PS.DATA, 'prior, decisive': '#8FB4D4',
             'prior, weak': PS.ACCENT_LT, 'P group': PS.SERIES[3], 'none': PS.ACCENT}
CONF_COL = {'high': PS.DATA, 'medium': '#8FB4D4', 'low': PS.ACCENT_LT, 'no call': PS.ACCENT}
TIER_COL = {1: PS.SERIES[2], 2: PS.SERIES[3], 3: PS.ACCENT}


FONT_DIR = C.ROOT / 'assets' / 'fonts'          # Liberation Sans: Arial metrics, SIL OFL
for _f in FONT_DIR.glob('LiberationSans-*.ttf'):
    matplotlib.font_manager.fontManager.addfont(str(_f))


def setup(scale):
    s = PS.setup_style(scale)
    plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'Arial', 'DejaVu Sans']
    return s


def gene_axis(ax, rows, fs):
    """Gene rows top to bottom with group bands; gene names in italics."""
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontstyle='italic', fontsize=fs)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    i = 0
    for k, (name, gs) in enumerate(GROUPS):
        gs = [g for g in gs if g in rows]
        if not gs:
            continue
        if k % 2 == 0:
            ax.axhspan(i - 0.5, i + len(gs) - 0.5, color='#F2F4F6', zorder=0, lw=0)
        i += len(gs)


def group_labels(fig, ax, rows, fs, x=-0.30):
    i = 0
    for name, gs in GROUPS:
        gs = [g for g in gs if g in rows]
        if not gs:
            continue
        mid = i + (len(gs) - 1) / 2
        ax.text(x, mid, name, transform=ax.get_yaxis_transform(), rotation=90, ha='center',
                va='center', fontsize=fs, color=PS.INK_SOFT)
        i += len(gs)


def stacked(ax, df, cols, colors, rows):
    left = np.zeros(len(rows))
    for c in cols:
        v = df.reindex(rows)[c].fillna(0).values * 100
        ax.barh(range(len(rows)), v, left=left, color=colors[c], height=0.72, lw=0, zorder=2)
        left += v
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.grid(axis='y', visible=False)


def letter(ax, k, s, dx=-0.02):
    ax.text(dx, 1.0, LETTERS[k], transform=ax.transAxes, fontsize=s['panel'], fontweight='bold',
            ha='right', va='bottom')


def callouts(ax, xs, ys, labels, x_text, ha, fs, color, gap):
    """Labels in a column at x_text, in the order of their points' y, at least `gap` apart
    (data units), each joined to its point by a thin leader line."""
    order = np.argsort(ys)
    yl = np.array(ys, float)[order]
    for i in range(1, len(yl)):
        yl[i] = max(yl[i], yl[i - 1] + gap)
    yl -= (yl.mean() - np.mean(np.array(ys)[order]))          # keep the column centred on its points
    for k, i in enumerate(order):
        ax.annotate(labels[i], (xs[i], ys[i]), xytext=(x_text, yl[k]), textcoords='data', fontsize=fs,
                    color=color, ha=ha, va='center',
                    arrowprops=dict(arrowstyle='-', color=PS.REFERENCE, lw=0.5, shrinkA=1, shrinkB=3))


def sqrt_ticks(lim):
    """Frequencies to label on a square-root axis whose upper end is sqrt(lim)."""
    cands = [0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0]
    out, last = [], -1
    for t in cands:
        v = np.sqrt(t)
        if v <= lim and v - last >= 0.16 * lim:
            out.append(t)
            last = v
    return out


def save(fig, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300)
    plt.close(fig)


# --------------------------------------------------------------------------- data
def load(results=C.RESULTS):
    d = {}
    d['ev'] = pd.read_csv(results / '02.evidence' / 'call_evidence.tsv', sep='\t', dtype=str).fillna('')
    d['gs'] = pd.read_csv(results / '02.evidence' / 'gene_summary.tsv', sep='\t').set_index('gene')
    d['sq'] = pd.read_csv(results / '02.evidence' / 'sample_qc.tsv', sep='\t')
    d['drb'] = pd.read_csv(results / '02.evidence' / 'drb_consistency.tsv', sep='\t', dtype=str)
    d['calls'] = pd.read_csv(results / '01.calls' / 'calls_long.tsv', sep='\t', dtype=str).fillna('')
    d['loc'] = pd.read_csv(results / '03.panels' / 'locus_concordance.tsv', sep='\t')
    d['al'] = pd.read_csv(results / '03.panels' / 'allele_concordance.tsv', sep='\t')
    d['hz'] = pd.read_csv(results / '03.panels' / 'homozygosity.tsv', sep='\t').set_index('gene')
    d['ceil'] = pd.read_csv(results / '03.panels' / 'panel_ceiling.tsv', sep='\t')
    d['cvp'] = pd.read_csv(results / '03.panels' / 'confidence_vs_panel.tsv', sep='\t')
    d['tiers'] = pd.read_csv(results / '04.confidence' / 'locus_tiers.tsv', sep='\t').set_index('gene')
    d['shift'] = pd.read_csv(results / '03.panels' / 'shift_direction.tsv', sep='\t')
    return d


def shares(d):
    """Per gene: state shares, 2-field basis shares, confidence shares."""
    ev = d['ev']
    st = ev.state.replace({'failed': 'no call', 'unresolved': 'no call', 'missing': 'no call'})
    state = pd.crosstab(ev.gene, st, normalize='index')
    al = pd.concat([ev[['gene', 'allele_1', 'basis_1', 'prior_ratio_1', 'copies']].set_axis(
                        ['gene', 'a', 'b', 'q', 'copies'], axis=1),
                    ev[['gene', 'allele_2', 'basis_2', 'prior_ratio_2', 'copies']].set_axis(
                        ['gene', 'a', 'b', 'q', 'copies'], axis=1).assign(second=1)])
    al = al[~((al.get('second') == 1) & (al.copies == '1'))]
    al = al[al.b != '']
    q = pd.to_numeric(al.q, errors='coerce')
    lab = np.select([al.b == 'agreed', (al.b == 'prior') & (q >= C.PRIOR_DOMINANCE),
                     al.b == 'prior', al.b == 'pgroup'],
                    ['agreed by reads', 'prior, decisive', 'prior, weak', 'P group'], 'none')
    basis = pd.crosstab(al.gene, lab, normalize='index')
    conf = pd.crosstab(ev.gene, ev.confidence.replace({'no_call': 'no call'}), normalize='index')
    return state, basis, conf


# --------------------------------------------------------------------------- figures
def fig_overview(d, out, scale='paper', width=PS.COL_DOUBLE):
    """Per gene: call state, how the 2-field allele was determined, call confidence."""
    s = setup(scale)
    state, basis, conf = shares(d)
    rows = [g for g in ORDER if g in state.index]
    fig, axes = plt.subplots(1, 3, figsize=(width, 0.21 * len(rows) * (1.35 if scale == 'slide' else 1) + 1.3),
                             sharey=True)
    plt.subplots_adjust(left=0.17, right=0.99, top=0.80, bottom=0.12, wspace=0.10)
    specs = [(state, [c for c in STATE_COL if c in state.columns], STATE_COL, 'Call state (% of samples)'),
             (basis, [c for c in BASIS_COL if c in basis.columns], BASIS_COL,
              'Source of the 2-field allele (% of alleles)'),
             (conf, [c for c in CONF_COL if c in conf.columns], CONF_COL, 'Call confidence (% of samples)')]
    for k, (ax, (df, cols, colors, xl)) in enumerate(zip(axes, specs)):
        stacked(ax, df, cols, colors, rows)
        gene_axis(ax, rows, s['tick'])
        ax.set_xlabel(xl, fontsize=s['label'])
        PS.despine(ax, grid_axis='x')
        handles = [Patch(color=colors[c], label=c) for c in cols]
        ax.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.0, 1.02), ncol=2,
                  fontsize=s['legend'], frameon=False, handlelength=1.0, columnspacing=0.8,
                  borderaxespad=0)
        letter(ax, k, s, dx=-0.03 if k else -0.30)
    group_labels(fig, axes[0], rows, s['tick'], x=-0.33)
    save(fig, out)


def fig_read_evidence(d, out, scale='paper', width=PS.COL_DOUBLE):
    """(a) unique-read balance of heterozygotes at obligate genes; (b) raw two-allele DRB3/4/5
    calls split by the copy number DRB1 implies; (c) DRB1-implied vs HLA-HD copy number."""
    s = setup(scale)
    ev = d['ev']
    fig = plt.figure(figsize=(width, width * (0.36 if scale == 'paper' else 0.40)))
    gsp = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.15, 1.10], wspace=0.40,
                           left=0.07, right=0.90, top=0.86, bottom=0.17)
    bins = np.linspace(0, 0.5, 26)
    # (a)
    ax = fig.add_subplot(gsp[0])
    h = ev[(ev.balance_assessed == 'True') & ev.gene.isin(C.OBLIGATE_GENES)]
    ms = pd.to_numeric(h.minor_share)
    ax.hist(ms, bins=bins, color=PS.DATA, lw=0, zorder=2)
    ax.axvline(C.MIN_MINOR_SHARE, color=PS.ACCENT, lw=1.1, ls='--', zorder=3)
    n_imb = int((ms < C.MIN_MINOR_SHARE).sum())
    ax.text(C.MIN_MINOR_SHARE + 0.01, 0.95, f'imbalanced\n{n_imb:,} of {len(ms):,}\n({n_imb / len(ms):.1%})',
            transform=ax.get_xaxis_transform(), va='top', ha='left', fontsize=s['legend'], color=PS.ACCENT)
    ax.set_xlabel('Minor allele share of unique reads', fontsize=s['label'])
    ax.set_ylabel('Heterozygous calls', fontsize=s['label'])
    ax.set_xlim(0, 0.5)
    PS.despine(ax, grid_axis='y')
    letter(ax, 0, s, dx=-0.12)
    ax.set_title('Obligate genes', fontsize=s['label'], loc='left', pad=4)
    # (b)
    ax = fig.add_subplot(gsp[1])
    c = d['calls']
    p = c[c.gene.isin(C.PARALOGUES) & (c.raw == 'two')].copy()
    u1, u2 = pd.to_numeric(p.unique_1), pd.to_numeric(p.unique_2)
    p['ms'] = np.minimum(u1, u2) / (u1 + u2)
    p = p[(u1 + u2) >= C.MIN_UNIQUE_READS]
    for e, col, lab in (('1', PS.ACCENT, '1 copy'), ('2', PS.DATA, '2 copies')):
        v = p[p.drb_expected == e].ms
        ax.hist(v, bins=bins, color=col, alpha=0.75 if e == '1' else 0.85, lw=0, zorder=2,
                label=f'{lab} (n = {len(v):,})')
    ax.axvline(C.MIN_MINOR_SHARE, color=PS.REFERENCE, lw=1.0, ls='--', zorder=3)
    ax.set_xlabel('Minor allele share of unique reads', fontsize=s['label'])
    ax.set_ylabel('Two-allele DRB3/4/5 calls', fontsize=s['label'])
    ax.set_xlim(0, 0.5)
    ax.legend(fontsize=s['legend'], frameon=False, loc='upper right', bbox_to_anchor=(1.0, 0.98),
              title='DRB1 implies', title_fontsize=s['legend'])
    PS.despine(ax, grid_axis='y')
    letter(ax, 1, s, dx=-0.12)
    ax.set_title('Two-allele DRB3/4/5 calls', fontsize=s['label'], loc='left', pad=4)
    # (c)
    ax = fig.add_subplot(gsp[2])
    drb = d['drb'].copy()
    drb['samples'] = drb.samples.astype(int)
    drb = drb[drb.expected != 'unknown']
    genes = C.PARALOGUES
    M = np.zeros((len(genes) * 3, 3))
    for i, g in enumerate(genes):
        for e in range(3):
            for o in range(3):
                v = drb[(drb.gene == g) & (drb.expected == str(e)) & (drb.observed == str(o))].samples.sum()
                M[i * 3 + e, o] = v
    rowsum = M.sum(axis=1, keepdims=True)
    F = np.divide(M, rowsum, out=np.zeros_like(M), where=rowsum > 0)
    ax.imshow(F, cmap=matplotlib.colors.LinearSegmentedColormap.from_list('b', ['#FFFFFF', PS.DATA]),
              vmin=0, vmax=1, aspect='auto')
    for i in range(M.shape[0]):
        for j in range(3):
            if rowsum[i, 0]:
                diag = (i % 3) == j
                ax.text(j, i, f'{int(M[i, j]):,}', ha='center', va='center', fontsize=s['legend'],
                        color='white' if F[i, j] > 0.55 else PS.INK, fontweight='bold' if diag else 'normal')
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['0', '1', '2'])
    ax.set_xlabel('Alleles reported by HLA-HD', fontsize=s['label'])
    ax.set_yticks(range(M.shape[0]))
    ax.set_yticklabels([str(e) for _ in genes for e in range(3)], fontsize=s['tick'])
    ax.set_ylabel('Copies implied by DRB1', fontsize=s['label'])
    for i, g in enumerate(genes):
        ax.axhline(i * 3 - 0.5, color=PS.INK, lw=0.8)
        tot = M[i * 3:(i + 1) * 3].sum()
        agree = sum(M[i * 3 + e, e] for e in range(3))
        ax.text(2.62, i * 3 + 0.55, g, ha='left', va='center', fontsize=s['tick'], fontstyle='italic',
                fontweight='bold', clip_on=False)
        ax.text(2.62, i * 3 + 1.45, f'{agree / tot:.0%} agree', ha='left', va='center',
                fontsize=s['legend'], color=PS.INK_SOFT, clip_on=False)
    ax.set_xlim(-0.5, 2.5)
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    letter(ax, 2, s, dx=-0.16)
    ax.set_title('Copy number: DRB1 rule vs HLA-HD', fontsize=s['label'], loc='left', pad=4)
    save(fig, out)


def fig_panel_frequencies(d, out, scale='paper', width=PS.COL_DOUBLE):
    """Our controls vs ToMMo, one panel per locus, square-root axes."""
    import nomenclature as N
    s = setup(scale)
    al = d['al'][d['al'].panel == 'ToMMo'].copy()
    loc = d['loc'][d['loc'].panel == 'ToMMo'].set_index('gene')
    vint = int(loc.panel_vintage.iloc[0])
    al['rel'] = [N.release_of(x) for x in al.allele]
    genes = [g for g in ORDER if g in loc.index]
    ncol = 5
    nrow = int(np.ceil((len(genes) + 1) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(width, width / ncol * nrow * 1.10))
    plt.subplots_adjust(left=0.07, right=0.99, top=0.96, bottom=0.08, wspace=0.34, hspace=0.46)
    axes = axes.ravel()
    ms = 9 if scale == 'paper' else 20
    for k, g in enumerate(genes):
        ax = axes[k]
        a = al[al.gene == g]
        m = loc.loc[g]
        lim = np.sqrt(max(a.freq_ours.max(), a.freq_panel.max()) * 1.08)
        ax.plot([0, lim], [0, lim], color=PS.REFERENCE, lw=0.8, zorder=1)
        absent = a.count_panel == 0
        late = absent & ~(a.rel <= vint)
        ax.scatter(np.sqrt(a.freq_panel[~absent]), np.sqrt(a.freq_ours[~absent]), s=ms, color=PS.DATA,
                   lw=0, zorder=3, alpha=0.85)
        ax.scatter(np.sqrt(a.freq_panel[absent & ~late]), np.sqrt(a.freq_ours[absent & ~late]), s=ms * 1.4,
                   marker='^', color=PS.ACCENT, lw=0, zorder=4)
        ax.scatter(np.sqrt(a.freq_panel[late]), np.sqrt(a.freq_ours[late]), s=ms * 1.4, marker='^',
                   facecolor='none', edgecolor=PS.REFERENCE, lw=0.7, zorder=4)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ticks = sqrt_ticks(lim)
        ax.set_xticks(np.sqrt(ticks))
        ax.set_yticks(np.sqrt(ticks))
        ax.set_xticklabels([f'{t:g}' for t in ticks], fontsize=s['tick'] - 1)
        ax.set_yticklabels([f'{t:g}' for t in ticks], fontsize=s['tick'] - 1)
        PS.despine(ax)
        ax.set_title(g, fontstyle='italic', fontsize=s['label'], loc='left', pad=2)
        bad = m.unconfirmed_assessable > C.TIER_MAX_UNCONFIRMED
        if m.frequency_comparable == 1:
            ax.text(0.97, 0.13, f"r = {m.r:.3f}", transform=ax.transAxes, ha='right', va='bottom',
                    fontsize=s['legend'] - 1, color=PS.INK_SOFT)
        else:
            ax.text(0.03, 0.97, 'frequencies\nnot comparable', transform=ax.transAxes, ha='left', va='top',
                    fontsize=s['legend'] - 1, color=PS.INK_SOFT)
        ax.text(0.97, 0.03, f"absent: {m.unconfirmed_assessable:.1%}", transform=ax.transAxes, ha='right',
                va='bottom', fontsize=s['legend'] - 1, color=PS.ACCENT if bad else PS.INK_SOFT,
                fontweight='bold' if bad else 'normal')
    for ax in axes[len(genes):]:
        ax.axis('off')
    axes[len(genes)].legend(handles=[
        Line2D([], [], marker='o', ls='', color=PS.DATA, label='allele in ToMMo'),
        Line2D([], [], marker='^', ls='', color=PS.ACCENT, label='absent from ToMMo'),
        Line2D([], [], marker='^', ls='', markerfacecolor='none', markeredgecolor=PS.REFERENCE,
               label=f'absent, named after\nToMMo\'s dictionary'),
        Line2D([], [], color=PS.REFERENCE, label='equal frequency')],
        loc='center left', frameon=False, fontsize=s['legend'], borderaxespad=0)
    fig.text(0.53, 0.012, 'ToMMo 61KJPN frequency (square-root scale)', ha='center', fontsize=s['label'])
    fig.text(0.008, 0.52, 'Frequency in our controls (square-root scale)', va='center', rotation=90,
             fontsize=s['label'])
    save(fig, out)


def fig_locus_summary(d, out, scale='paper', width=PS.COL_DOUBLE):
    """One row per locus: the three tier criteria and the tier."""
    s = setup(scale)
    t = d['tiers']
    rows = [g for g in ORDER if g in t.index]
    fig, axes = plt.subplots(1, 4, figsize=(width, 0.21 * len(rows) * (1.3 if scale == 'slide' else 1) + 1.15),
                             sharey=True, gridspec_kw={'width_ratios': [1.3, 1.3, 1.3, 0.42]})
    plt.subplots_adjust(left=0.165, right=0.99, top=0.88, bottom=0.15, wspace=0.16)
    msz = 4.2 if scale == 'paper' else 7
    na_kw = dict(ha='center', va='center', fontsize=s['legend'], color=PS.REFERENCE)

    def col(ok):
        return PS.DATA if ok else PS.ACCENT

    # (a) alleles absent from ToMMo although its dictionary had them
    ax = axes[0]
    for i, g in enumerate(rows):
        v = t.at[g, 'unconfirmed'] if 'unconfirmed' in t.columns else np.nan
        if pd.isna(v):
            ax.text(0.03, i, 'no panel', transform=ax.get_yaxis_transform(), **{**na_kw, 'ha': 'left'})
            continue
        ok = t.at[g, 'pass_unconfirmed'] == 1
        ax.plot(v, i, 'o', color=col(ok), ms=msz, zorder=4)
        vh = t.at[g, 'unconfirmed_high']
        if not pd.isna(vh):
            ax.plot(vh, i, 'o', mfc='white', mec=col(ok), ms=msz, mew=1.0, zorder=5)
    ax.axvline(C.TIER_MAX_UNCONFIRMED, color=PS.REFERENCE, ls='--', lw=0.9)
    ax.set_xlabel('Chromosomes on an allele\nabsent from ToMMo', fontsize=s['label'])
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1, decimals=0))
    ax.set_xlim(left=-0.002)
    # (b) dropout
    ax = axes[1]
    for i, g in enumerate(rows):
        v = t.at[g, 'dropout'] if 'dropout' in t.columns else np.nan
        if pd.isna(v):
            ax.text(0.03, i, 'not applicable', transform=ax.get_yaxis_transform(), **{**na_kw, 'ha': 'left'})
            continue
        ok = t.at[g, 'pass_dropout'] == 1
        lo, hi = max(t.at[g, 'dropout_lo'], -0.05), min(t.at[g, 'dropout_hi'], 0.15)
        ax.plot([lo, hi], [i, i], color=col(ok), lw=1.2, zorder=3)
        ax.plot(v, i, 'o', color=col(ok), ms=msz, zorder=4)
    ax.axvline(C.TIER_MAX_DROPOUT, color=PS.REFERENCE, ls='--', lw=0.9)
    ax.axvline(0, color=PS.INK, lw=0.6)
    ax.set_xlim(-0.05, 0.15)
    ax.set_xticks([0, 0.05, 0.10])
    ax.set_xlabel('Heterozygotes reported\nas homozygotes', fontsize=s['label'])
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1, decimals=0))
    # (c) high-confidence share
    ax = axes[2]
    for i, g in enumerate(rows):
        ok = t.at[g, 'pass_call_quality'] == 1
        ax.plot(t.at[g, 'share_high'], i, 'o', color=col(ok), ms=msz, zorder=4)
    ax.axvline(C.TIER_MIN_HIGH_CONF, color=PS.REFERENCE, ls='--', lw=0.9)
    ax.set_xlabel('High-confidence\ncalls', fontsize=s['label'])
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1, decimals=0))
    ax.set_xlim(0.75, 1.005)
    ax.set_xticks([0.8, 0.9, 1.0])
    # (d) tier
    ax = axes[3]
    for i, g in enumerate(rows):
        tr = int(t.at[g, 'tier'])
        ax.add_patch(matplotlib.patches.FancyBboxPatch((0.2, i - 0.34), 0.6, 0.68, boxstyle='round,pad=0.02',
                                                       color=TIER_COL[tr], lw=0))
        ax.text(0.5, i, str(tr), ha='center', va='center', color='white', fontweight='bold', fontsize=s['tick'])
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.set_xlabel('Tier', fontsize=s['label'])
    for k, ax in enumerate(axes):
        gene_axis(ax, rows, s['tick'])
        PS.despine(ax, grid_axis='x')
        if k == 3:
            ax.grid(False)
            ax.spines['bottom'].set_visible(False)
        letter(ax, k, s, dx=-0.03 if k else -0.30)
    group_labels(fig, axes[0], rows, s['tick'], x=-0.34)
    fig.legend(handles=[Line2D([], [], marker='o', ls='', color=PS.DATA, label='criterion met'),
                        Line2D([], [], marker='o', ls='', color=PS.ACCENT, label='criterion not met'),
                        Line2D([], [], marker='o', ls='', mfc='white', mec=PS.DATA,
                               label='high-confidence calls only'),
                        Line2D([], [], color=PS.REFERENCE, ls='--', label='threshold')],
               loc='upper center', ncol=4, frameon=False, fontsize=s['legend'], bbox_to_anchor=(0.55, 0.995))
    save(fig, out)


def fig_shift_direction(d, out, scale='paper', width=PS.COL_ONEHALF):
    """Alleles >= 1 point apart between our controls and ToMMo: where does 1KG JPT lie?"""
    s = setup(scale)
    sh = d['shift'].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(width, width * 0.82))
    plt.subplots_adjust(left=0.14, right=0.97, top=0.88, bottom=0.13)
    lim = max(sh.diff_ours.abs().max(), sh.diff_kg.abs().max()) * 1.6
    for x0 in (0, -lim):
        ax.add_patch(matplotlib.patches.Rectangle((x0, x0), lim, lim, color='#EAF1F7', zorder=0, lw=0))
    ax.axhline(0, color=PS.INK, lw=0.6)
    ax.axvline(0, color=PS.INK, lw=0.6)
    same = sh.kg_side == 'ours'
    sz = 16 if scale == 'paper' else 34
    ax.scatter(sh.diff_ours[same], sh.diff_kg[same], s=sz, color=PS.DATA, zorder=3,
               label=f'1KG JPT on our side ({int(same.sum())})')
    ax.scatter(sh.diff_ours[~same], sh.diff_kg[~same], s=sz, color=PS.ACCENT, zorder=3,
               label=f"1KG JPT on ToMMo's side ({int((~same).sum())})")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    lab = sh[sh.diff_ours.abs() >= 0.019]
    fs = s['legend'] - 1
    for side, x_text, ha in (('neg', -lim * 0.97, 'left'), ('pos', lim * 0.97, 'right')):
        q = lab[lab.diff_ours < 0] if side == 'neg' else lab[lab.diff_ours > 0]
        if len(q):
            callouts(ax, q.diff_ours.values, q.diff_kg.values, q.allele.values, x_text, ha, fs,
                     PS.INK_SOFT, gap=lim * 0.085)
    ax.set_xlabel('Our controls − ToMMo (frequency)', fontsize=s['label'])
    ax.set_ylabel('1KG JPT − ToMMo (frequency)', fontsize=s['label'])
    ax.legend(loc='lower left', bbox_to_anchor=(0, 1.0), ncol=2, frameon=False, fontsize=s['legend'],
              handletextpad=0.3, columnspacing=1.2, borderaxespad=0.2)
    PS.despine(ax)
    save(fig, out)


def fig_panel_ceiling(d, out, scale='paper', width=PS.COL_ONEHALF):
    """D against 1KG JPT: ours vs ToMMo's, with the sampling expectation."""
    s = setup(scale)
    c = d['ceil'].set_index('gene')
    kg = d['loc'][d['loc'].panel == '1KG JPT'].set_index('gene')
    genes = [g for g in ORDER if g in c.index]
    fig, ax = plt.subplots(figsize=(width, width * 0.62))
    plt.subplots_adjust(left=0.13, right=0.98, top=0.86, bottom=0.14)
    x = np.arange(len(genes))
    w = 0.34
    ax.bar(x - w / 2, c.loc[genes, 'D_ours_vs_kg'], w, color=PS.DATA, label='our controls vs 1KG JPT', zorder=2)
    ax.bar(x + w / 2, c.loc[genes, 'D_tommo_vs_kg'], w, color='#8FB4D4', label='ToMMo vs 1KG JPT', zorder=2)
    ax.scatter(x - w / 2, kg.loc[genes, 'D_null_mean'], marker='_', s=120 if scale == 'paper' else 260,
               color=PS.INK, zorder=3, label='expected from sampling alone')
    ax.set_xticks(x)
    ax.set_xticklabels(genes, fontstyle='italic')
    ax.set_ylabel('Discordance D', fontsize=s['label'])
    PS.despine(ax, grid_axis='y')
    ax.legend(frameon=False, fontsize=s['legend'], loc='lower left', bbox_to_anchor=(0, 1.0), ncol=2,
              handlelength=1.2, columnspacing=1.0)
    save(fig, out)


def fig_confidence_vs_panel(d, out, scale='paper', width=PS.COL_ONEHALF):
    """Share of chromosomes on alleles absent from ToMMo, by call confidence."""
    s = setup(scale)
    cv = d['cvp']
    pooled = cv[cv.gene == 'all ToMMo loci'].set_index('confidence').reindex(['high', 'medium', 'low']).dropna()
    fig, ax = plt.subplots(figsize=(width, width * 0.58))
    plt.subplots_adjust(left=0.14, right=0.97, top=0.92, bottom=0.18)
    x = np.arange(len(pooled))
    k, n = pooled.n_unconfirmed.values, pooled.n_chr.values
    p = k / n
    lo, hi = zip(*[wilson(a, b) for a, b in zip(k, n)])
    ax.bar(x, p, color=[CONF_COL[c] for c in pooled.index], width=0.6, zorder=2)
    ax.errorbar(x, p, yerr=[p - np.array(lo), np.array(hi) - p], fmt='none', ecolor=PS.INK, lw=1, capsize=3,
                zorder=3)
    for i in range(len(x)):
        ax.text(x[i], hi[i] + 0.004, f'{p[i]:.1%}\n({int(k[i]):,} / {int(n[i]):,})', ha='center', va='bottom',
                fontsize=s['legend'])
    ax.set_xticks(x)
    ax.set_xticklabels([f'{c} confidence' for c in pooled.index])
    ax.set_ylabel('Chromosomes on an allele\nabsent from ToMMo', fontsize=s['label'])
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1, decimals=0))
    ax.set_ylim(0, max(hi) * 1.35)
    PS.despine(ax, grid_axis='y')
    save(fig, out)


def fig_sample_qc(d, out, scale='paper', width=PS.COL_DOUBLE):
    """(a) read pairs per sample; (b) reported genes without a confident call per sample."""
    s = setup(scale)
    sq = d['sq']
    fig, axes = plt.subplots(1, 2, figsize=(width, width * 0.32))
    plt.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=0.2, wspace=0.28)
    ax = axes[0]
    v = sq.read_pairs / 1e3
    ax.hist(v, bins=np.logspace(np.log10(v.min() * 0.95), np.log10(v.max() * 1.05), 60), color=PS.DATA,
            lw=0, zorder=2)
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: f'{x:,.0f}'))
    cut = np.percentile(sq.read_pairs, C.OUTLIER_READ_PAIRS_PCTL) / 1e3
    ax.axvline(cut, color=PS.ACCENT, ls='--', lw=1)
    ax.text(0.98, 0.95, f'dashed line: 1st percentile\n({cut:,.0f}k read pairs)', transform=ax.transAxes,
            va='top', ha='right', fontsize=s['legend'], color=PS.ACCENT)
    ax.set_xlabel('Read pairs given to HLA-HD (thousands)', fontsize=s['label'])
    ax.set_ylabel('Samples', fontsize=s['label'])
    PS.despine(ax, grid_axis='y')
    letter(ax, 0, s, dx=-0.1)
    ax = axes[1]
    cnt = sq.low_conf.value_counts().sort_index()
    ax.bar(cnt.index, cnt.values, color=PS.DATA, width=0.8, zorder=2)
    for xx, yy in cnt.items():
        ax.text(xx, yy, f'{yy:,}', ha='center', va='bottom', fontsize=s['legend'] - 1)
    ax.set_xlabel('Reported genes with a low-confidence or no call', fontsize=s['label'])
    ax.set_ylabel('Samples', fontsize=s['label'])
    ax.set_xticks(cnt.index)
    PS.despine(ax, grid_axis='y')
    letter(ax, 1, s, dx=-0.1)
    save(fig, out)


def wilson(k, n, z=1.959963984540054):
    if not n:
        return np.nan, np.nan
    p, dd = k / n, 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return max(0, (c - h) / dd), min(1, (c + h) / dd)
