#!/usr/bin/env python3
"""Deck figures: concept schematics and result figures, drawn at slide size (12.6 in wide = the slide's
content width, so 11-13 pt figure type reads as 11-13 pt on the slide).

One visual system throughout: class I in blue, class II in green, peptides and anything that must not be
missed in red, reference data in gray. Structures are drawn from PDB files with molview.py.

Usage: python3 make_figures.py            (reads ../results and ../assets; writes figures/)
"""
import csv
import json
import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.patches import (Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Polygon,  # noqa: E402
                                Rectangle)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'scripts'))
sys.path.insert(0, str(HERE))
import config as C  # noqa: E402
import figlib as F  # noqa: E402
import molview as MV  # noqa: E402
import plot_style as PS  # noqa: E402

OUT = HERE / 'figures'
OUT.mkdir(exist_ok=True)
ASSETS = C.ROOT / 'assets'
S = json.loads((HERE / 'stats.json').read_text())

INK, SOFT, GREY, NEUTRAL = PS.INK, PS.INK_SOFT, PS.REFERENCE, PS.NEUTRAL
BLUE, BLUE_D, BLUE_L, BLUE_XL = PS.DATA, PS.DATA_DARK, '#8FB4D4', '#E3ECF5'
GREEN, GREEN_D, GREEN_L, GREEN_XL = PS.SERIES[2], '#2A5E4B', '#9CCBB5', '#E2F0E8'
RED, RED_L = PS.ACCENT, PS.ACCENT_LT
OCHRE, OCHRE_L = PS.SERIES[3], '#EFE3C8'
PALE = '#F4F6F8'
CLASS_I = ['A', 'B', 'C', 'E', 'F', 'G']
CLASSICAL = {'A', 'B', 'C', 'DRA', 'DRB1', 'DRB3', 'DRB4', 'DRB5', 'DQA1', 'DQB1', 'DPA1', 'DPB1'}
DPI = 250


def setup():
    F.setup('slide')
    plt.rcParams.update({'axes.titleweight': 'bold'})


def canvas(w, h):
    setup()
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, w)
    ax.set_ylim(0, h)
    ax.axis('off')
    return fig, ax


def rbox(ax, x, y, w, h, fc, ec=None, lw=1.0, r=0.06, z=2, ls='-'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f'round,pad=0,rounding_size={r}', fc=fc,
                                ec=ec or fc, lw=lw, zorder=z, ls=ls))


def arrow(ax, x0, y0, x1, y1, col=INK, lw=1.4, rad=0.0, ms=14, z=5, ls='-'):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=ms, color=col, lw=lw,
                                 zorder=z, connectionstyle=f'arc3,rad={rad}', ls=ls))


def tag(ax, k, x, y, s=16):
    ax.text(x, y, 'abcdefgh'[k], fontsize=s, fontweight='bold', color=INK, va='top', ha='left')


def save(fig, name):
    fig.savefig(OUT / name, dpi=DPI)
    plt.close(fig)


def gi(ax, x, y, gene, **kw):
    """Gene symbol in italics."""
    return ax.text(x, y, gene, fontstyle='italic', **kw)


# =========================================================================== 1 function
def _membrane(ax, x0, x1, y, col='#C9A66B'):
    xs = np.linspace(x0, x1, 200)
    ax.fill_between(xs, y - 0.07, y + 0.07, color='#E9D9B8', lw=0, zorder=1)
    ax.plot(xs, np.full_like(xs, y + 0.07), color=col, lw=1.2, zorder=1)
    ax.plot(xs, np.full_like(xs, y - 0.07), color=col, lw=1.2, zorder=1)


def _peptide(ax, x, y, n=5, r=0.055, gap=0.12):
    for k in range(n):
        ax.add_patch(Circle((x + (k - (n - 1) / 2) * gap, y), r, fc=RED, ec='white', lw=0.6, zorder=7))


def _class_i(ax, x, y):
    """Class I on a membrane at height y: alpha3 stalk + beta2m, alpha1/alpha2 platform, peptide."""
    ax.plot([x, x], [y + 0.07, y + 0.35], color=BLUE, lw=3, zorder=3)
    ax.add_patch(Ellipse((x, y + 0.55), 0.36, 0.42, fc=BLUE_L, ec=BLUE_D, lw=0.8, zorder=4))
    ax.add_patch(Circle((x + 0.32, y + 0.42), 0.15, fc='#C9CED4', ec=GREY, lw=0.8, zorder=4))
    for dx in (-0.19, 0.19):
        ax.add_patch(Ellipse((x + dx, y + 0.95), 0.40, 0.36, fc=BLUE, ec=BLUE_D, lw=0.8, zorder=5))
    _peptide(ax, x, y + 1.13)


def _class_ii(ax, x, y):
    """Class II: alpha and beta chains, alpha1/beta1 platform, longer peptide."""
    for dx, c_low, c_top in ((-0.19, GREEN_L, '#6DB592'), (0.19, '#8FC1A9', GREEN)):
        ax.plot([x + dx, x + dx], [y + 0.07, y + 0.35], color=GREEN, lw=3, zorder=3)
        ax.add_patch(Ellipse((x + dx, y + 0.55), 0.34, 0.40, fc=c_low, ec=GREEN_D, lw=0.8, zorder=4))
        ax.add_patch(Ellipse((x + dx, y + 0.95), 0.40, 0.36, fc=c_top, ec=GREEN_D, lw=0.8, zorder=5))
    _peptide(ax, x, y + 1.13, n=7, gap=0.105)


def _tcell(ax, x, y, col, coreceptor, reach):
    """T cell above the HLA: membrane arc, T-cell receptor touching the peptide, and the co-receptor reaching
    the lower domain of the HLA molecule at `reach` (CD8 to alpha3, CD4 to beta2)."""
    cell = Ellipse((x, y + 1.25), 4.2, 2.5, fc=PALE, ec=col, lw=1.6, zorder=1)
    ax.add_patch(cell)
    cell.set_clip_path(Rectangle((x - 2.2, y - 0.1), 4.4, 0.95, transform=ax.transData))
    for dx in (-0.1, 0.1):
        ax.plot([x + dx, x + dx], [y, y - 0.45], color=INK, lw=2.6, zorder=6, solid_capstyle='round')
    ax.add_patch(FancyBboxPatch((x - 0.21, y - 0.62), 0.42, 0.2, boxstyle='round,pad=0,rounding_size=0.08',
                                fc='#5B6168', ec=INK, lw=0.6, zorder=6))
    ax.plot([x + 0.85, reach[0]], [y + 0.04, reach[1]], color=OCHRE, lw=3.2, zorder=3, solid_capstyle='round')
    ax.text(x + 0.82, y - 0.3, coreceptor, fontsize=10.5, color=OCHRE, fontweight='bold', va='center')
    ax.text(x - 0.3, y - 0.3, 'T-cell\nreceptor', fontsize=10, color=INK, ha='right', va='center')


def hla_function():
    """What HLA molecules do: class I shows peptides from inside the cell to CD8+ T cells, class II
    peptides from proteins taken up from outside to CD4+ T cells (Medhasi & Chantratita 2022)."""
    fig, ax = canvas(12.6, 4.75)
    MEM = 1.62
    for k, (ttl, col, tname, cor, cell, cellsub) in enumerate((
            ('HLA class I  (HLA-A, -B, -C)', BLUE, 'CD8+ cytotoxic T cell', 'CD8',
             'Presenting cell: nearly any nucleated cell', 'peptides from proteins made inside the cell'),
            ('HLA class II  (HLA-DR, -DQ, -DP)', GREEN, 'CD4+ helper T cell', 'CD4',
             'Presenting cell: macrophage, B cell, dendritic cell', 'peptides from proteins taken up from outside'))):
        ox = 0.25 + k * 6.35
        cx = ox + 2.75
        ax.text(ox, 4.55, ttl, fontsize=14, fontweight='bold', color=col, va='center')
        _membrane(ax, ox + 0.2, ox + 5.9, MEM)
        ax.fill_between([ox + 0.2, ox + 5.9], 0.05, MEM - 0.07, color='#FBF6EC', lw=0, zorder=0)
        ax.text(ox + 0.35, 0.5, cell, fontsize=11, color=INK, fontweight='bold', va='center')
        ax.text(ox + 0.35, 0.2, cellsub, fontsize=10.5, color=SOFT, va='center', fontstyle='italic')
        (_class_i if k == 0 else _class_ii)(ax, cx, MEM)
        reach = (cx + 0.24, MEM + 0.66) if k == 0 else (cx + 0.38, MEM + 0.62)
        _tcell(ax, cx, MEM + 1.88, col, cor, reach)
        ax.text(cx, MEM + 2.27, tname, ha='center', fontsize=12, fontweight='bold', color=col, zorder=3)
        ax.text(cx - 0.62, MEM + 0.85, 'HLA', fontsize=11.5, fontweight='bold', color=col, ha='right', va='center')
        ax.annotate('peptide', (cx - 0.3, MEM + 1.13), (cx - 1.25, MEM + 1.38), fontsize=10.5, color=RED,
                    va='center', ha='right', arrowprops=dict(arrowstyle='-', color=RED, lw=0.8))
        if k == 0:
            px, py = ox + 0.45, 1.12
            t = np.linspace(0, 1, 60)
            ax.plot(px + 1.0 * t, py + 0.08 * np.sin(t * 6 * np.pi), color=OCHRE, lw=2.2, zorder=3)
            ax.text(px, py + 0.22, 'own or viral protein', fontsize=9.5, color=SOFT)
            arrow(ax, px + 1.08, py, px + 1.42, py, col=SOFT, lw=1.0, ms=10)
            _peptide(ax, px + 1.72, py, n=3, r=0.045, gap=0.1)
            arrow(ax, px + 1.9, py + 0.06, cx - 0.12, MEM + 0.32, col=SOFT, lw=1.0, rad=-0.25, ms=10)
        else:
            bx, by = ox + 0.85, 2.55
            ax.add_patch(FancyBboxPatch((bx - 0.32, by - 0.13), 0.64, 0.26, boxstyle='round,pad=0,rounding_size=0.13',
                                        fc=OCHRE_L, ec=OCHRE, lw=1.2, zorder=3))
            ax.text(bx, by + 0.25, 'microbe', ha='center', fontsize=9.5, color=SOFT)
            vx, vy = ox + 1.55, 1.05
            ax.add_patch(Circle((vx, vy), 0.26, fc='white', ec=OCHRE, lw=1.2, zorder=3))
            _peptide(ax, vx, vy, n=3, r=0.045, gap=0.1)
            arrow(ax, bx + 0.1, by - 0.15, vx - 0.12, vy + 0.28, col=SOFT, lw=1.0, rad=0.2, ms=10)
            ax.text(vx + 0.36, vy, 'taken up and\ncut into peptides', fontsize=9.5, color=SOFT, va='center')
            arrow(ax, vx + 0.2, vy + 0.2, cx - 0.36, MEM + 0.3, col=SOFT, lw=1.0, rad=-0.3, ms=10)
    save(fig, 'concept_hla_function.png')


# =========================================================================== 2 structure
def al_(x):
    """Allele name for a figure label."""
    return x


def _entropy(gene):
    """Shannon entropy (bits) at each IMGT position across the named alleles of `gene` (IPD-IMGT/HLA 3.64.0),
    over alleles with a determined residue there."""
    ref = pd.read_csv(C.D_EXPORTS / 'residue_reference' / f'{gene}.tsv', sep='\t', index_col=0, dtype=str)
    H = {}
    for col in ref.columns:
        v = ref[col]
        v = v[~v.isin(['*', '.'])]
        p = v.value_counts(normalize=True).values
        H[int(col)] = float(-(p * np.log2(p)).sum()) if len(p) else 0.0
    return H, len(ref)


def _label(ax, xy, txt, xyt, col=INK, fs=11, bold=True):
    ax.annotate(txt, xy, xyt, fontsize=fs, color=col, fontweight='bold' if bold else 'normal', ha='center',
                va='center', arrowprops=dict(arrowstyle='-', color=GREY, lw=0.7, shrinkA=2, shrinkB=2), zorder=20)


def hla_structure():
    """(a) class I (HLA-A*02:01 with a viral peptide, PDB 1HHK) and (b) class II (HLA-DR1 with an influenza
    peptide, PDB 1DLH), side views by domain; (c) the class I groove from above, colored by how much the
    HLA-A alleles differ at each position."""
    setup()
    fig = plt.figure(figsize=(12.6, 4.6))
    axs = [fig.add_axes(r) for r in ([0.0, 0.0, 0.29, 0.91], [0.29, 0.0, 0.29, 0.91], [0.585, 0.15, 0.40, 0.76])]
    fig.text(0.012, 0.965, '(a)  Class I  (HLA-A*02:01)', fontsize=12.5, fontweight='bold', color=BLUE, va='center')
    fig.text(0.302, 0.965, '(b)  Class II  (HLA-DR1)', fontsize=12.5, fontweight='bold', color=GREEN, va='center')
    fig.text(0.59, 0.965, '(c)  Groove viewed from above, as a T cell sees it', fontsize=12.5, fontweight='bold',
             color=INK, va='center')
    for ax in axs:
        ax.set_aspect('equal')
        ax.axis('off')
    pdb = ASSETS / 'pdb'
    c1 = MV.read_pdb(pdb / '1HHK.pdb', 'ABC')
    c2 = MV.read_pdb(pdb / '1DLH.pdb', 'ABC')
    # (a) class I
    up = MV.centroid(c1, [('A', 1, 180)]) - MV.centroid(c1, [('A', 182, 275), ('B', 0, 100)])
    u, x, z = MV.frame(up, MV.direction(c1, 'A', 57, 84))
    v1 = MV.tilt((x, u, z), 25)
    dom1 = {'a1': BLUE, 'a2': '#7FAAD3', 'a3': '#BACFE4', 'b2m': '#CDD1D6'}

    def col1(c, r):
        return dom1['b2m'] if c == 'B' else dom1['a1'] if r <= 90 else dom1['a2'] if r <= 180 else dom1['a3']
    ax = axs[0]
    MV.draw(ax, c1, v1, col1, chains='AB', width=5.0)
    MV.peptide(ax, c1, 'C', v1, RED, size=34)
    x0, x1, y0, y1 = MV.extent(c1, v1)
    ax.set_xlim(x0 - 14, x1 + 14)
    ax.set_ylim(y0 - 13, y1 + 6)
    P = lambda sel: tuple(MV.centroid(c1, sel) @ np.array([v1[0], v1[1]]).T)   # noqa: E731
    _label(ax, P([('A', 1, 85)]), 'α1', (x0 - 7, y1 - 6), BLUE)
    _label(ax, P([('A', 95, 178)]), 'α2', (x1 + 7, y1 - 4), BLUE)
    _label(ax, P([('A', 185, 270)]), 'α3', (x0 - 7, y0 + 20), BLUE)
    _label(ax, P([('B', 1, 99)]), 'β2m', (x1 + 8, y0 + 18), GREY)
    pep = c1['C']['xyz'].mean(axis=0) @ np.array([v1[0], v1[1]]).T
    _label(ax, tuple(pep), 'peptide', (x1 + 6, y1 + 2), RED)
    ax.plot([x0 - 6, x1 + 6], [y0 - 5, y0 - 5], color='#C9A66B', lw=3, solid_capstyle='butt')
    ax.text((x0 + x1) / 2, y0 - 7.5, 'toward the cell membrane', ha='center', va='top', fontsize=9.5, color=SOFT)
    # (b) class II
    up2 = MV.centroid(c2, [('A', 1, 80), ('B', 1, 92)]) - MV.centroid(c2, [('A', 82, 182), ('B', 94, 190)])
    u2, x2, z2 = MV.frame(up2, MV.direction(c2, 'B', 56, 90))
    v2 = MV.tilt((x2, u2, z2), 25)
    dom2 = {'a1': '#6DB592', 'a2': '#BCDFCC', 'b1': GREEN_D, 'b2': '#93C4AC'}

    def col2(c, r):
        if c == 'A':
            return dom2['a1'] if r <= 81 else dom2['a2']
        return dom2['b1'] if r <= 94 else dom2['b2']
    ax = axs[1]
    MV.draw(ax, c2, v2, col2, chains='AB', width=5.0)
    MV.peptide(ax, c2, 'C', v2, RED, size=34)
    x0, x1, y0, y1 = MV.extent(c2, v2)
    ax.set_xlim(x0 - 14, x1 + 14)
    ax.set_ylim(y0 - 13, y1 + 6)
    P2 = lambda sel: tuple(MV.centroid(c2, sel) @ np.array([v2[0], v2[1]]).T)   # noqa: E731
    _label(ax, P2([('A', 3, 80)]), 'α1', (x0 - 7, y1 - 8), GREEN)
    _label(ax, P2([('B', 3, 92)]), 'β1', (x1 + 7, y1 - 4), GREEN_D)
    _label(ax, P2([('A', 85, 180)]), 'α2', (x0 - 7, y0 + 18), GREEN)
    _label(ax, P2([('B', 95, 190)]), 'β2', (x1 + 7, y0 + 18), GREEN)
    pep = c2['C']['xyz'].mean(axis=0) @ np.array([v2[0], v2[1]]).T
    _label(ax, tuple(pep), 'peptide', (x0 - 4, y1 + 3), RED)
    ax.plot([x0 - 6, x1 + 6], [y0 - 5, y0 - 5], color='#C9A66B', lw=3, solid_capstyle='butt')
    ax.text((x0 + x1) / 2, y0 - 7.5, 'toward the cell membrane', ha='center', va='top', fontsize=9.5, color=SOFT)
    # (c) groove from above, colored by variability among HLA-A alleles
    H, n_al = _entropy('A')
    cm = LinearSegmentedColormap.from_list('var', ['#E6E8EB', '#EFC36A', RED, '#5A1715'])
    HMAX = 2.0
    ax = axs[2]
    vt = (x, z, u)
    pep2 = c1['C']['xyz'] @ np.array([vt[0], vt[1]]).T
    ang = np.degrees(np.arctan2(*(pep2[-1] - pep2[0])[::-1]))
    a = np.radians(-ang)                                       # rotate in the screen plane: peptide horizontal
    vt = (np.cos(a) * vt[0] - np.sin(a) * vt[1], np.sin(a) * vt[0] + np.cos(a) * vt[1], vt[2])
    MV.draw(ax, c1, vt, lambda c, r: cm(min(H.get(r, 0.0) / HMAX, 1.0)), chains='A', keep=lambda c, r: r <= 180,
            fade=0.22, width=6.5)
    px, py = MV.peptide(ax, c1, 'C', vt, '#202020', size=30)
    sel = c1['A']['xyz'][c1['A']['res'] <= 180]
    sx, sy = sel @ vt[0], sel @ vt[1]
    ax.set_xlim(sx.min() - 4, sx.max() + 16)
    ax.set_ylim(sy.min() - 3, sy.max() + 6)
    pos = {r: (p @ vt[0], p @ vt[1]) for r, p in zip(c1['A']['res'], c1['A']['xyz'])}
    hel1 = np.mean([pos[r] for r in range(57, 85)], axis=0)
    hel2 = np.mean([pos[r] for r in range(140, 176)], axis=0)
    ax.text(sx.max() + 1.5, hel1[1], 'α1 helix', fontsize=11, fontweight='bold', color=INK, va='center')
    ax.text(sx.max() + 1.5, hel2[1], 'α2 helix', fontsize=11, fontweight='bold', color=INK, va='center')
    ax.text(sx.max() + 1.5, py.mean(), 'peptide', fontsize=11, fontweight='bold', color='#202020', va='center')
    top = []
    for p in sorted((p for p in H if 1 <= p <= 180 and p in pos), key=lambda p: -H[p]):
        if all(np.hypot(pos[p][0] - pos[q][0], pos[p][1] - pos[q][1]) > 5.5 for q in top):
            top.append(p)
        if len(top) == 6:
            break
    for p in top:
        ax.text(pos[p][0], pos[p][1], str(p), fontsize=8.5, color='white', fontweight='bold', ha='center',
                va='center', zorder=40, bbox=dict(boxstyle='circle,pad=0.18', fc='#3A2A2A', ec='none', alpha=0.85))
    cax = fig.add_axes([0.61, 0.075, 0.17, 0.028])
    cb = plt.colorbar(plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(0, HMAX)), cax=cax, orientation='horizontal')
    cb.set_ticks([0, 1, 2])
    cb.ax.tick_params(labelsize=9.5)
    cb.outline.set_visible(False)
    fig.text(0.61, 0.125, f'Variability among the {n_al:,} HLA-A alleles (bits)', fontsize=10, color=INK, va='bottom')
    fig.text(0.795, 0.088, 'numbered: the most variable\npositions (IMGT numbering)', fontsize=9.5, color=SOFT,
             va='center')
    save(fig, 'concept_hla_structure.png')
    return {'n_alleles_A': n_al, 'top_positions_A': top}


# =========================================================================== 3 map
def mhc_map():
    """Chromosome 6 ideogram (UCSC hg38 bands) and the HLA genes at their GRCh38 positions (Ensembl), with
    the region whose reads were given to HLA-HD."""
    co = pd.read_csv(C.PROJECT / 'analysis' / 'assoc_hla' / 'results' / '00.prep' / 'hla_gene_coords.tsv', sep='\t')
    co['mid'] = (co.gene_start + co.gene_end) / 2e6
    bands = pd.read_csv(ASSETS / 'ucsc' / 'cytoBand.hg38.chr6.txt', sep='\t', header=None,
                        names=['chrom', 'start', 'end', 'name', 'stain'])
    fig, ax = canvas(12.6, 4.75)
    # ideogram
    L = bands.end.max() / 1e6
    X0, X1, YI, HI = 0.5, 12.1, 4.2, 0.26
    sx = lambda mb: X0 + (X1 - X0) * mb / L   # noqa: E731
    shade = {'gneg': '#FFFFFF', 'gpos25': '#D9D9D9', 'gpos50': '#A6A6A6', 'gpos75': '#737373', 'gpos100': '#404040',
             'acen': '#C98B8B', 'gvar': '#BFBFBF', 'stalk': '#E0E0E0'}
    for b in bands.itertuples():
        ax.add_patch(Rectangle((sx(b.start / 1e6), YI), sx(b.end / 1e6) - sx(b.start / 1e6), HI,
                               fc=shade.get(b.stain, '#FFFFFF'), ec='none', zorder=2))
    rbox(ax, X0, YI, X1 - X0, HI, 'none', ec=GREY, lw=1.0, r=0.12, z=3)
    ax.text(X0, YI + HI + 0.1, 'chromosome 6', fontsize=11, color=SOFT)
    R0, R1 = 29.54, 33.42                                  # region extracted for typing (config: HLA region)
    ax.add_patch(Rectangle((sx(R0), YI - 0.05), sx(R1) - sx(R0), HI + 0.1, fc='none', ec=RED, lw=2.0, zorder=4))
    ax.text(sx(R1) + 0.1, YI + HI + 0.1, '6p22.1–p21.32: the HLA region', fontsize=11, color=RED)
    # zoom
    Z0, Z1, ZX0, ZX1, YG = 29.55, 33.25, 0.5, 12.1, 1.95
    zx = lambda mb: ZX0 + (ZX1 - ZX0) * (mb - Z0) / (Z1 - Z0)   # noqa: E731
    ax.plot([sx(R0), zx(Z0)], [YI - 0.05, 3.45], color=RED, lw=0.8, ls='--')
    ax.plot([sx(R1), zx(Z1)], [YI - 0.05, 3.45], color=RED, lw=0.8, ls='--')
    b_end = co.loc[co.gene == 'B', 'gene_end'].item() / 1e6
    dra = co.loc[co.gene == 'DRA', 'gene_start'].item() / 1e6
    for lab, a_, b_, col in (('class I region', Z0, b_end + 0.05, BLUE), ('class III region', b_end + 0.05, dra - 0.05, GREY),
                             ('class II region', dra - 0.05, Z1, GREEN)):
        ax.add_patch(Rectangle((zx(a_), 1.55), zx(b_) - zx(a_), 1.9, fc=col, alpha=0.08, ec='none', zorder=0))
        ax.text((zx(a_) + zx(b_)) / 2, 3.18, lab, ha='center', fontsize=12, fontweight='bold', color=col)
    ax.text((zx(b_end) + zx(dra)) / 2, 3.08, 'complement, TNF, …; no HLA molecules', ha='center', va='top',
            fontsize=9.5, color=SOFT)
    ax.plot([ZX0, ZX1], [YG, YG], color=INK, lw=1.2, zorder=1)
    for t in np.arange(30.0, 33.01, 0.5):
        ax.plot([zx(t)] * 2, [1.55, 1.62], color=GREY, lw=0.8)
        ax.text(zx(t), 1.45, f'{t:.1f}', ha='center', va='top', fontsize=9.5, color=SOFT)
    ax.text(ZX1, 1.45, 'Mb', ha='left', va='top', fontsize=9.5, color=SOFT)
    # genes: a tick at the true position; class I labels above their gene, class II labels fanned out over the
    # gene-free class III interval in two rows, each joined to its gene by a leader line
    rows = co.sort_values('mid').reset_index(drop=True)
    c1 = rows[rows.gene.isin(CLASS_I)]
    c2 = rows[~rows.gene.isin(CLASS_I)].reset_index(drop=True)
    lab_x = {}
    last = -9
    for g in c1.itertuples():
        xm = zx(g.mid)
        lane = 1 if xm - last < 0.45 else 0
        last = xm if lane == 0 else last
        lab_x[g.gene] = (xm, YG + 0.42 + lane * 0.36)
    xa, xb = zx(31.75), ZX1 - 0.25
    for i, g in enumerate(c2.itertuples()):
        lab_x[g.gene] = (xa + (xb - xa) * i / (len(c2) - 1), YG + 0.42 + (i % 2) * 0.36)
    for g in rows.itertuples():
        col = BLUE if g.gene in CLASS_I else GREEN
        cls = g.gene in CLASSICAL
        xm = zx(g.mid)
        lx, ly = lab_x[g.gene]
        ax.add_patch(Rectangle((zx(g.gene_start / 1e6), YG - 0.09), max(zx(g.gene_end / 1e6) - zx(g.gene_start / 1e6), 0.03),
                               0.18, fc=col, ec=col, zorder=3))
        ax.plot([xm, xm, lx], [YG + 0.09, YG + 0.2, ly - 0.02], color=col, lw=0.8, zorder=2)
        gi(ax, lx, ly, g.gene, ha='center', va='bottom', fontsize=10.5, color='white' if cls else col, fontweight='bold',
           zorder=4, bbox=dict(boxstyle='round,pad=0.18', fc=col if cls else 'white', ec=col, lw=1.0))
    ax.text(0.5, 0.9, 'Filled labels: genes of the classical molecules HLA-A, -B, -C, -DR, -DQ, and -DP, which present '
            'peptides to T cells. Open labels: other HLA genes.', fontsize=10.5, color=INK)
    ax.text(0.5, 0.55, 'DRB3 and DRB4 are not on the GRCh38 primary chromosome 6, which carries a DR haplotype with '
            'DRB5. Red box: reads extracted for typing (29.54–33.42 Mb).', fontsize=10.5, color=INK)
    ax.text(0.5, 0.2, 'Positions: Ensembl, GRCh38; bands: UCSC hg38.', fontsize=9.5, color=SOFT)
    save(fig, 'concept_mhc_map.png')


# =========================================================================== 4 inheritance
HAP = {'a': ('#3C6E9F', ['A*24:02', 'B*52:01', 'DRB1*15:02']), 'b': ('#9DBEDC', ['A*02:01', 'B*35:01', 'DRB1*04:05']),
       'c': ('#B08A3C', ['A*11:01', 'B*15:01', 'DRB1*09:01']), 'd': ('#E5CF9C', ['A*26:01', 'B*40:02', 'DRB1*08:03'])}


def _chrom_pair(ax, x, y, haps, w=0.34, h=1.25, labels=True):
    """Two chromosome-6 homologues as vertical bars, each carrying one haplotype; alleles printed beside."""
    for k, hp in enumerate(haps):
        col, als = HAP[hp]
        xx = x + k * (w + 0.12)
        rbox(ax, xx, y, w, h, col, ec='white', lw=1, r=0.12, z=3)
        ax.text(xx + w / 2, y + h + 0.1, hp, ha='center', fontsize=11, fontweight='bold', color=col)
    if labels:
        for j, a in enumerate(HAP[haps[0]][1]):
            ax.text(x - 0.08, y + h - 0.25 - j * 0.38, a, ha='right', va='center', fontsize=9, color=INK)
        for j, a in enumerate(HAP[haps[1]][1]):
            ax.text(x + 2 * w + 0.2, y + h - 0.25 - j * 0.38, a, ha='left', va='center', fontsize=9, color=INK)


def inheritance():
    """Pedigree: each child receives one HLA haplotype from each parent as a block (alleles illustrative);
    siblings share 2, 1 or 0 haplotypes with probabilities 25%, 50%, 25% (Choo 2007)."""
    fig, ax = canvas(12.6, 4.7)
    # parents
    for x, sym, haps, lab in ((2.2, 'sq', 'ab', 'father'), (6.2, 'ci', 'cd', 'mother')):
        if sym == 'sq':
            ax.add_patch(Rectangle((x - 0.2, 3.85), 0.4, 0.4, fc='white', ec=INK, lw=1.4))
        else:
            ax.add_patch(Circle((x, 4.05), 0.2, fc='white', ec=INK, lw=1.4))
        ax.text(x, 4.38, lab, fontsize=11.5, ha='center', va='bottom', fontweight='bold')
        _chrom_pair(ax, x - 0.4, 2.35, haps, h=1.2)
    ax.plot([2.4, 6.0], [4.05, 4.05], color=INK, lw=1.2)
    ax.plot([4.2, 4.2], [4.05, 2.1], color=INK, lw=1.2)
    ax.plot([0.85, 7.55], [2.1, 2.1], color=INK, lw=1.2)
    for k, haps in enumerate(('ac', 'ad', 'bc', 'bd')):
        x = 0.85 + k * 2.23
        ax.plot([x, x], [2.1, 1.95], color=INK, lw=1.2)
        ax.add_patch(Polygon([(x, 1.95), (x + 0.15, 1.8), (x, 1.65), (x - 0.15, 1.8)], fc='white', ec=INK, lw=1.2))
        ax.text(x + 0.22, 1.8, f'{haps[0]} + {haps[1]}', fontsize=11, va='center', fontweight='bold')
        _chrom_pair(ax, x - 0.4, 0.35, haps, h=1.0, w=0.3, labels=False)
        ax.text(x + 0.42, 0.6, '25%', fontsize=11, color=SOFT, va='center')
    ax.text(0.3, 0.08, 'Alleles are examples. Each bar is one chromosome 6; its color marks the haplotype it carries.',
            fontsize=9.5, color=SOFT)
    # sibling sharing
    ax2 = fig.add_axes([0.715, 0.2, 0.26, 0.6])
    setup()
    v = [25, 50, 25]
    ax2.bar([0, 1, 2], v, color=[BLUE, BLUE_L, NEUTRAL], width=0.62)
    for i, t in enumerate(v):
        ax2.text(i, t + 1.5, f'{t}%', ha='center', fontsize=12, fontweight='bold')
    ax2.set_xticks([0, 1, 2])
    ax2.set_xticklabels(['both', 'one', 'none'], fontsize=11.5)
    ax2.set_ylim(0, 62)
    ax2.set_yticks([])
    ax2.set_xlabel('HLA haplotypes shared by two siblings', fontsize=11.5)
    ax2.grid(False)
    for sp in ('top', 'right', 'left'):
        ax2.spines[sp].set_visible(False)
    save(fig, 'concept_inheritance.png')


# =========================================================================== 5 diversity
def _history():
    """Named alleles per gene in every IPD-IMGT/HLA release of the pinned allele history."""
    with open(C.IMGT_HISTORY) as fh:
        rd = csv.reader(ln for ln in fh if not ln.startswith('#'))
        head = next(rd)
        rels = [int(r) for r in head[1:]]
        cnt = {r: Counter() for r in rels}
        for row in rd:
            for r, nm in zip(rels, row[1:]):
                nm = nm.strip()
                if '*' in nm and nm != 'NA':
                    cnt[r][nm.split('*')[0]] += 1
    return cnt


def diversity():
    """(a) alleles named per reported gene, release 3.64.0; (b) growth of the named alleles across releases."""
    cnt = _history()
    rels = sorted(cnt)
    last = cnt[max(rels)]
    setup()
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.6, 4.35), gridspec_kw={'width_ratios': [1, 1.15]})
    plt.subplots_adjust(left=0.07, right=0.985, top=0.9, bottom=0.14, wspace=0.28)
    genes = [g for g in F.ORDER]
    y = np.arange(len(genes))
    v = [last.get(g, 0) for g in genes]
    cols = [(BLUE if g in CLASS_I else GREEN) if g in CLASSICAL else (BLUE_L if g in CLASS_I else GREEN_L) for g in genes]
    a.barh(y, v, color=cols, height=0.72)
    for i, t in enumerate(v):
        a.text(t + 150, i, f'{t:,}', va='center', fontsize=10, color=INK)
    a.set_yticks(y)
    a.set_yticklabels(genes, fontstyle='italic', fontsize=10.5)
    a.set_ylim(len(genes) - 0.5, -0.5)
    a.set_xlim(0, max(v) * 1.2)
    a.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda t, _: f'{t:,.0f}'))
    a.set_xlabel(f'Alleles named (IPD-IMGT/HLA {rel_label(max(rels))})', fontsize=12)
    a.legend(handles=[Line2D([], [], marker='s', ls='', ms=10, color=BLUE, label='class I'),
                      Line2D([], [], marker='s', ls='', ms=10, color=GREEN, label='class II'),
                      Line2D([], [], marker='s', ls='', ms=10, color=GREY, alpha=0.5, label='light: non-classical')],
             frameon=False, fontsize=10.5, loc='lower right')
    PS.despine(a, grid_axis='x')
    a.set_title('(a) Alleles per gene', loc='left', fontsize=13, fontweight='bold')
    # growth since release 3.0.0 (January 2010 nomenclature; HLA-C was "Cw" before it)
    rels3 = [r for r in rels if r >= 3000]
    show = [('B', BLUE_D), ('C', BLUE_L), ('A', BLUE), ('DRB1', GREEN_D), ('DPB1', GREEN), ('DQB1', GREEN_L)]
    xr = np.arange(len(rels3))
    ends = {}
    for g, col in show:
        yy = [cnt[r].get(g, 0) for r in rels3]
        b.plot(xr, yy, color=col, lw=2.2)
        ends[g] = yy[-1]
    ymax = max(ends.values()) * 1.08
    lab_y = dict(zip(sorted(ends, key=ends.get), np.array(sorted(ends.values()), float)))
    keys = sorted(ends, key=ends.get)
    for i in range(1, len(keys)):
        lab_y[keys[i]] = max(lab_y[keys[i]], lab_y[keys[i - 1]] + 0.055 * ymax)
    for g, col in show:
        b.annotate(g, (xr[-1], ends[g]), (xr[-1] + 2.2, lab_y[g]), fontstyle='italic', fontsize=11, color=col,
                   va='center', fontweight='bold', arrowprops=dict(arrowstyle='-', color=col, lw=0.6))
    ticks = [i for i, r in enumerate(rels3) if (r // 10) % 10 == 0]
    b.set_xticks(ticks)
    b.set_xticklabels([rel_label(rels3[i]) for i in ticks], fontsize=10.5)
    b.set_xlim(-0.5, len(rels3) + 6)
    b.set_ylim(0, ymax)
    b.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda t, _: f'{t / 1000:.0f}k' if t else '0'))
    b.set_xlabel('IPD-IMGT/HLA release', fontsize=12)
    b.set_ylabel('Alleles named per gene', fontsize=12)
    PS.despine(b, grid_axis='y')
    t0, t1 = sum(cnt[rels3[0]].values()), sum(cnt[rels3[-1]].values())
    b.text(0.02, 0.97, f'All HLA genes: {t0:,} alleles in {rel_label(rels3[0])}, {t1:,} in {rel_label(rels3[-1])}',
           transform=b.transAxes, fontsize=10.5, color=SOFT, va='top')
    b.set_title('(b) Growth across database releases', loc='left', fontsize=13, fontweight='bold')
    fig.savefig(OUT / 'fig_diversity.png', dpi=DPI)
    plt.close(fig)
    first = 3000
    return {'first_release': rel_label(first), 'last_release': rel_label(max(rels)),
            'total_first': sum(cnt[first].values()), 'total_last': sum(cnt[max(rels)].values()),
            'B_first': cnt[first]['B'], 'B_last': cnt[max(rels)]['B']}


def rel_label(r):
    return f'{r // 1000}.{(r // 10) % 100}.{r % 10}'


# =========================================================================== 6 prior reports
def prior_cteph():
    """HLA alleles reported in Japanese CTEPH, both studies on the same two scales, from the full texts in ../ref:
    carrier frequencies (Table 2 of each) and odds ratios with 95% CI exactly as the authors report them (none
    recomputed). Tanabe 2005, all CTEPH, 80 patients and 678 controls: OR in Table 2, 95% CI in the Results text.
    Kominami 2009, CTEPH without DVT, 99 patients and 380 controls: OR and 95% CI in Table 2."""
    setup()
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=True)
    plt.subplots_adjust(left=0.27, right=0.97, top=0.88, bottom=0.2, wspace=0.1)
    # allele, study, carriers % (patients, controls), OR (95% CI)
    rows = [('B*52:01', 'Tanabe 2005', 40, 24, 2.14, 1.29, 3.55),
            ('B*52:01', 'Kominami 2009, no DVT', 45.5, 25.3, 2.47, 1.56, 3.90),
            ('DPB1*02:02', 'Tanabe 2005', 19, 6, 3.41, 1.71, 6.74),
            ('DPB1*02:02', 'Kominami 2009, no DVT', 19.2, 4.5, 5.07, 2.52, 10.19)]
    ys = [0.75, 1.75, 3.55, 4.55]
    for y, (_al, _st, p, c, o, lo, hi) in zip(ys, rows):
        a.plot([c, p], [y, y], color=NEUTRAL, lw=4, solid_capstyle='round', zorder=1)
        a.plot(c, y, 'o', ms=9, color=GREY, zorder=2)
        a.plot(p, y, 'o', ms=9, color=RED, zorder=3)
        a.text(p + 2.5, y, f'{p:g}%', va='center', fontsize=10, color=RED, fontweight='bold')
        a.text(c, y + 0.42, f'{c:g}%', va='center', ha='center', fontsize=9.5, color=SOFT)
        b.plot([lo, hi], [y, y], color=RED, lw=2)
        b.plot(o, y, 's', ms=8, color=RED)
        b.text(hi * 1.1, y, f'{o:.2f}', va='center', fontsize=10, color=RED, fontweight='bold')
    a.set_yticks(ys)
    a.set_yticklabels([r[1] for r in rows], fontsize=10)
    for name, y in (('B*52:01', 0.0), ('DPB1*02:02', 2.8)):
        a.text(-0.03, y, name, transform=a.get_yaxis_transform(), ha='right', va='center', fontsize=11,
               fontweight='bold', color=INK)
    a.set_ylim(5.2, -0.45)
    a.set_xlim(0, 62)
    a.set_xlabel('Carriers (%)', fontsize=11)
    a.legend(handles=[Line2D([], [], marker='o', ls='', color=RED, ms=7, label='patients'),
                      Line2D([], [], marker='o', ls='', color=GREY, ms=7, label='controls')],
             frameon=False, fontsize=10, loc='lower left', bbox_to_anchor=(-0.02, 0.98), ncol=2, handletextpad=0.2,
             columnspacing=0.8)
    PS.despine(a, grid_axis='x')
    b.axvline(1, color=GREY, lw=1, ls='--')
    b.set_xscale('log')
    b.set_xticks([1, 2, 5, 10])
    b.set_xticklabels(['1', '2', '5', '10'])
    b.set_xlim(0.8, 20)
    b.set_xlabel('Odds ratio (95% CI)', fontsize=11)
    PS.despine(b, grid_axis='x')
    b.tick_params(axis='y', length=0)
    save(fig, 'fig_prior_cteph.png')


# =========================================================================== 7 nomenclature
def nomenclature():
    """Anatomy of an allele name; P group (A*02:09 and A*02:66 differ from A*02:01:01:01 only in exon 4, at
    residues 236 and 225, and all three are in A*02:01P — IPD-IMGT 3.64.0)."""
    fig, ax = canvas(12.6, 4.35)
    parts = [('HLA-A', INK, 'gene'), ('*', SOFT, None), ('24', BLUE, '1'), (':', SOFT, None),
             ('02', BLUE, '2'), (':', SOFT, None), ('01', GREY, '3'), (':', SOFT, None), ('01', GREY, '4'),
             ('N', RED, 'N')]
    x = 0.35
    boxes = {}
    for t, c, tg in parts:
        txt = ax.text(x, 3.4, t, fontsize=36, color=c, fontweight='bold', va='center', ha='left')
        fig.canvas.draw()
        bb = txt.get_window_extent().transformed(ax.transData.inverted())
        boxes[tg] = bb
        if tg:
            ax.plot([bb.x0 + 0.03, bb.x1 - 0.03], [2.85, 2.85], color=c if c != SOFT else INK, lw=2.2)
            ax.text((bb.x0 + bb.x1) / 2, 2.6, tg, ha='center', va='center', fontsize=13, color=c, fontweight='bold')
        x = bb.x1 + 0.02
    x0, x1 = boxes['1'].x0 - 0.06, boxes['2'].x1 + 0.06
    rbox(ax, x0, 2.98, x1 - x0, 0.85, 'none', ec=BLUE, lw=1.4, r=0.08, ls='--')
    ax.text((x0 + x1) / 2, 4.0, '2-field resolution = the protein', ha='center', fontsize=12, color=BLUE,
            fontweight='bold')
    rows = [('gene', INK, 'HLA-A, HLA-B, …'), ('1', BLUE, 'allele group (often the serological antigen)'),
            ('2', BLUE, 'a distinct protein (amino-acid sequence)'), ('3', GREY, 'synonymous change in the coding sequence'),
            ('4', GREY, 'change outside the coding sequence'), ('N', RED, 'suffix: not expressed (null allele)')]
    for k, (tg, c, txt) in enumerate(rows):
        y = 2.15 - k * 0.34
        ax.text(0.35, y, tg, fontsize=12, color=c, fontweight='bold', va='center')
        ax.text(1.1, y, txt, fontsize=12, color=INK, va='center')
    ox = 6.85
    ax.text(ox, 4.0, 'P group: identical protein where the peptide is held', fontsize=13, fontweight='bold', color=INK)
    ex = [('exon 1', 0.45), ('exon 2', 0.9), ('exon 3', 0.9), ('exon 4', 0.9), ('exon 5', 0.45)]
    alle = [('A*02:01:01:01', None, ''), ('A*02:09', 3, 'residue 236'), ('A*02:66', 3, 'residue 225')]
    for k, (al, diff, lab) in enumerate(alle):
        y = 3.15 - k * 0.6
        ax.text(ox, y + 0.12, al, fontsize=11.5, va='center', color=INK)
        xx = ox + 1.55
        for j, (en, w) in enumerate(ex):
            rbox(ax, xx, y, w, 0.24, '#CFE0EF' if j in (1, 2) else '#EEF0F2', ec=GREY, lw=0.6, r=0.02)
            if diff is not None and j == diff:
                xm = xx + w * (0.62 if k == 1 else 0.45)
                ax.plot([xm, xm], [y - 0.02, y + 0.26], color=RED, lw=2.6, zorder=4)
                ax.text(xm, y - 0.07, lab, fontsize=9.5, color=RED, va='top', ha='center')
            if k == 0:
                ax.text(xx + w / 2, y + 0.42, en, ha='center', fontsize=10, color=SOFT)
            xx += w + 0.05
    ax.add_patch(Rectangle((ox + 2.02, 1.72), 1.88, 1.8, fill=False, ec=BLUE, lw=1.6, ls='--', zorder=6))
    ax.text(ox + 2.96, 1.6, 'peptide-binding domain (exons 2 and 3)', ha='center', va='top', fontsize=10.5, color=BLUE)
    ax.text(ox, 1.0, 'The three proteins differ only outside exons 2 and 3, so they form one P group,\n'
            'A*02:01P (IPD-IMGT/HLA 3.64.0). ToMMo reports alleles at this resolution.', fontsize=11, color=SOFT, va='top')
    save(fig, 'concept_nomenclature.png')


# =========================================================================== 8 DR haplotypes
PAR_COL = {'DRB3': OCHRE, 'DRB4': PS.SERIES[4], 'DRB5': BLUE}


def _drb_track(ax, x, y, drb1, par, w1=1.45, w2=0.85, h=0.36, fs=10.5):
    """One chromosome's DR region: DRB1 box, then the second DRB gene (or an empty slot)."""
    ax.plot([x - 0.12, x + w1 + w2 + 0.42], [y + h / 2] * 2, color=GREY, lw=1.2, zorder=1)
    rbox(ax, x, y, w1, h, GREEN, r=0.05)
    ax.text(x + w1 / 2, y + h / 2, drb1, ha='center', va='center', fontsize=fs, color='white', fontweight='bold',
            fontstyle='normal' if '*' in drb1 else 'italic')
    x2 = x + w1 + 0.3
    if par:
        rbox(ax, x2, y, w2, h, PAR_COL[par], r=0.05)
        gi(ax, x2 + w2 / 2, y + h / 2, par, ha='center', va='center', fontsize=fs, color='white', fontweight='bold')
    else:
        rbox(ax, x2, y, w2, h, 'white', ec=GREY, lw=1, r=0.05, ls='--')
        ax.text(x2 + w2 / 2, y + h / 2, 'empty', ha='center', va='center', fontsize=fs - 0.5, color=SOFT)


def dr_haplotypes():
    """Left: what one chromosome 6 carries, by DRB1 group (hla.alleles.org DRB haplotypes). Right: a person has
    two chromosomes, so DRB3, DRB4, and DRB5 are each present in 0, 1 or 2 copies (examples illustrative)."""
    fig, ax = canvas(12.6, 4.55)
    # left: the rule
    ax.text(0.2, 4.3, 'One chromosome 6', fontsize=13.5, fontweight='bold', color=INK, va='center')
    ax.text(0.2, 3.95, 'slot 1 always holds DRB1; its group decides slot 2', fontsize=10.5, color=SOFT, va='center')
    rows = [('DRB1*15, *16', 'DR51', 'DRB5'), ('DRB1*03, *11, *12, *13, *14', 'DR52', 'DRB3'),
            ('DRB1*04, *07, *09', 'DR53', 'DRB4'), ('DRB1*01, *08, *10', 'DR1, DR8, DR10', None)]
    ax.text(0.2, 3.5, 'DRB1 group', fontsize=10, color=SOFT, va='center')
    ax.text(2.75, 3.5, 'slot 1', fontsize=10, color=SOFT, va='center', ha='center')
    ax.text(3.95, 3.5, 'slot 2', fontsize=10, color=SOFT, va='center', ha='center')
    ax.text(4.75, 3.5, 'haplotype', fontsize=10, color=SOFT, va='center')
    for k, (grp, name, par) in enumerate(rows):
        y = 2.9 - k * 0.72
        ax.text(0.2, y + 0.18, grp, fontsize=10.5, color=INK, va='center')
        _drb_track(ax, 2.2, y, 'DRB1', par, w1=1.1, w2=0.85)
        ax.text(4.75, y + 0.18, name, fontsize=10.5, color=INK, va='center', fontweight='bold')
    ax.text(0.2, 0.25, 'DRA, on every chromosome, is not drawn. Exceptions to the rule occur.', fontsize=9.5,
            color=SOFT)
    ax.plot([6.05, 6.05], [0.2, 4.4], color=NEUTRAL, lw=1)
    # right: three people
    ox = 6.35
    ax.text(ox, 4.3, 'One person: two chromosomes, add up', fontsize=13.5, fontweight='bold', color=INK, va='center')
    people = [(('DRB1*04:05', 'DRB4'), ('DRB1*15:02', 'DRB5')),
              (('DRB1*04:05', 'DRB4'), ('DRB1*09:01', 'DRB4')),
              (('DRB1*01:01', None), ('DRB1*08:03', None))]
    ax.text(ox + 4.47, 4.0, 'copies per person', fontsize=10, color=SOFT, va='center', ha='center')
    for j, g in enumerate(('DRB3', 'DRB4', 'DRB5')):
        gi(ax, ox + 3.75 + j * 0.72, 3.68, g, fontsize=10, color=PAR_COL[g], ha='center', va='center', fontweight='bold')
    for k, person in enumerate(people):
        y0 = 3.08 - k * 1.05
        ax.text(ox - 0.05, y0 - 0.05, 'ABC'[k], fontsize=12, fontweight='bold', color=INK, va='center', ha='right')
        for i, (d1, par) in enumerate(person):
            _drb_track(ax, ox + 0.15, y0 - i * 0.42, d1, par, w1=1.45, w2=0.75, h=0.32, fs=10)
        cnt = {g: sum(1 for _, p in person if p == g) for g in ('DRB3', 'DRB4', 'DRB5')}
        for j, g in enumerate(('DRB3', 'DRB4', 'DRB5')):
            cx = ox + 3.75 + j * 0.72
            rbox(ax, cx - 0.25, y0 - 0.3, 0.5, 0.5, PAR_COL[g] if cnt[g] else 'white', ec=PAR_COL[g], lw=1.2, r=0.06)
            ax.text(cx, y0 - 0.05, str(cnt[g]), ha='center', va='center', fontsize=13, fontweight='bold',
                    color='white' if cnt[g] else INK)
    ax.text(ox, 0.25, 'Person C carries none of the three: 0 copies is a correct result, not a failed typing.',
            fontsize=10, color=INK)
    save(fig, 'concept_dr_haplotypes.png')


# =========================================================================== 10 reads
def _exon2(names):
    """Exon 2 coding sequence of each allele from the pinned IPD-IMGT/HLA 3.64.0 nucleotide alignment."""
    import re
    from collections import defaultdict
    seqs, order = defaultdict(str), []
    for ln in open(ASSETS / 'imgt' / 'A_nuc.txt'):
        m = re.match(r'^ (A\*\S+)\s+(.*)$', ln.rstrip('\n'))
        if m:
            if m.group(1) not in seqs:
                order.append(m.group(1))
            seqs[m.group(1)] += m.group(2).replace(' ', '')
    ref = seqs[order[0]]
    ex = {n: ''.join(r if c == '-' else c for c, r in zip(seqs[n], ref)).split('|')[1] for n in names}
    keep = [i for i in range(len(ex[names[0]])) if any(ex[n][i] != '.' for n in names)]
    return {n: ''.join(ex[n][i] for i in keep) for n in names}


def _letters(ax, x0, y, seq, w, col=INK, hi=(), hicol=RED, fs=10.5, bold_hi=True):
    for k, ch in enumerate(seq):
        on = k in hi
        ax.text(x0 + (k + 0.5) * w, y, ch, ha='center', va='center', fontsize=fs, family='monospace',
                color=hicol if on else col, fontweight='bold' if (on and bold_hi) else 'normal', zorder=5)


def short_reads():
    """Reads over a real stretch of HLA-A exon 2: A*24:02 and A*11:01 differ at two of these 48 bases
    (IPD-IMGT/HLA 3.64.0); reads shortened for display. Right: a look-alike gene (schematic)."""
    A1, A2 = 'A*24:02:01:01', 'A*11:01:01:01'
    ex = _exon2([A1, A2])
    diff = [i for i, (p, q) in enumerate(zip(ex[A1], ex[A2])) if p != q]
    st = max(0, diff[0] - 20)
    s1, s2 = ex[A1][st:st + 48], ex[A2][st:st + 48]
    d = [i - st for i in diff if st <= i < st + 48]
    fig, ax = canvas(12.6, 3.95)
    ax.set_ylim(0.45, 4.4)                                      # drop the empty band under the caption
    X0, W = 2.05, 0.125
    ax.text(0.25, 4.15, "One person's two HLA-A alleles", fontsize=13.5, fontweight='bold', color=INK)
    ax.text(0.25, 3.82, f'48 bases of exon 2, real sequence; the alleles differ at {len(d)} of them (red)',
            fontsize=10.5, color=SOFT)
    for k in d:
        ax.add_patch(Rectangle((X0 + k * W, 1.2), W, 2.35, fc='#F6DCDA', ec='none', zorder=0))
    for lab, sq, y, col in ((A1[:7], s1, 3.32, BLUE), (A2[:7], s2, 2.95, OCHRE)):
        ax.text(X0 - 0.15, y, lab, ha='right', va='center', fontsize=11.5, color=col, fontweight='bold')
        _letters(ax, X0, y, sq, W, hi=d)
    ax.plot([0.25, X0 + 48 * W], [2.68, 2.68], color=NEUTRAL, lw=0.8)
    rows = [('fits both', GREY, 2.42), (f'fits {A1[:7]} only', BLUE_D, 1.92), (f'fits {A2[:7]} only', '#7A5A1E', 1.42)]
    for lab, col, y in rows:
        ax.text(X0 - 0.15, y, lab, ha='right', va='center', fontsize=10.5, color=col, fontweight='bold')
    reads = [(0, 15, None, 2.42), (33, 15, None, 2.42), (13, 15, 1, 1.92), (18, 15, 2, 1.42)]
    for start, ln, who, y in reads:
        seq = (s1 if who != 2 else s2)[start:start + ln]
        fc, ec = ('#EEF0F2', GREY) if who is None else (BLUE_XL, BLUE) if who == 1 else ('#F7EEDC', OCHRE)
        rbox(ax, X0 + start * W - 0.02, y - 0.17, ln * W + 0.04, 0.34, fc, ec=ec, lw=1.2, r=0.08, z=3)
        _letters(ax, X0 + start * W, y, seq, W, hi=[k - start for k in d if start <= k < start + ln],
                 hicol=BLUE_D if who == 1 else '#7A5A1E', fs=10)
    ax.text(X0 + 48 * W + 0.1, 2.42, 'no\ninformation', va='center', fontsize=10, color=SOFT)
    ax.text(X0 + 48 * W + 0.1, 1.67, 'tell the\nalleles apart', va='center', fontsize=10, color=INK)
    ax.text(0.25, 0.75, 'Each box is a read, aligned under the alleles; reads are shortened for display (real reads are '
            'about 150 bases).', fontsize=10.5, color=SOFT)
    # right: look-alike gene (schematic letters)
    ox = 9.55
    ax.plot([ox - 0.3, ox - 0.3], [1.0, 4.2], color=NEUTRAL, lw=1)
    ax.text(ox, 4.15, 'A look-alike gene', fontsize=13.5, fontweight='bold', color=INK)
    ax.text(ox, 3.82, 'schematic', fontsize=10.5, color=SOFT)
    g1, g2 = 'GGTATTTCTCCACATCC', 'GGTATTTCTCCACATGC'
    w2 = 0.125
    ax.text(ox, 3.32, 'HLA-A', fontsize=11.5, color=BLUE, fontweight='bold', va='center', fontstyle='italic')
    ax.text(ox, 2.95, 'HLA-Y', fontsize=11.5, color=GREY, fontweight='bold', va='center', fontstyle='italic')
    _letters(ax, ox + 0.75, 3.32, g1, w2, hi=[15])
    _letters(ax, ox + 0.75, 2.95, g2, w2, hi=[15])
    rbox(ax, ox + 0.75 + 1 * w2 - 0.02, 2.25, 11 * w2 + 0.04, 0.34, '#EEF0F2', ec=GREY, lw=1.2, r=0.08, z=3)
    _letters(ax, ox + 0.75 + w2, 2.42, g2[1:12], w2, fs=10)
    ax.text(ox, 1.92, 'A read from the pseudogene\nHLA-Y fits HLA-A as well: it\ncan be counted for the wrong\n'
            'gene and create a false allele.', fontsize=10.5, color=INK, va='top')
    ax.set_xlim(0, 12.6)
    save(fig, 'concept_short_reads.png')
    return {'diff_in_window': len(d), 'window_start_in_exon2': st + 1}


# =========================================================================== 11 HLA-HD
def residue_concept():
    """Residue columns, on real HLA-A sequences of the IPD-IMGT/HLA 3.64.0 protein alignment (the rows the exports
    match each 2-field allele to): (a) why a residue is tested as well as alleles — alleles of three groups share
    Y at position 9; (b) how one person's cells are filled; (c) when a cell is blank — a person in the cohort whose
    second allele is sequenced only over the groove exons — and why a 0 there would be wrong."""
    ref = pd.read_csv(C.D_EXPORTS / 'residue_reference' / 'A.tsv', sep='\t', index_col=0, dtype=str)
    mm = pd.read_csv(C.D_EXPORTS / 'allele_match_map.tsv', sep='\t').set_index('allele')
    aa = lambda a, ps: ''.join(ref.loc[mm.at[a, 'matched_to'], str(p)] for p in ps)
    rcols = pd.read_csv(C.D_EXPORTS / 'residue_dosage.tsv', sep='\t', nrows=0).columns
    seen = lambda p: sorted(c.split(':')[2] for c in rcols if c.startswith(f'A:{p}:'))
    fig, ax = canvas(12.6, 4.55)
    W, BAND = 0.27, '#FBE9E7'
    for x, h in ((0.12, '(a)  Why residues as well as alleles'), (4.42, '(b)  How a residue column is filled'),
                 (8.62, '(c)  When a cell is blank')):
        ax.text(x, 4.36, h, fontsize=12.5, fontweight='bold', color=INK, va='center')
    for x in (4.27, 8.47):
        ax.plot([x, x], [0.25, 4.5], color=NEUTRAL, lw=0.8)

    def ruler(x0, ps, mark, y):
        for k, q in enumerate(ps):
            ax.text(x0 + (k + 0.5) * W, y, str(q), ha='center', va='center', fontsize=8, rotation=0 if q < 100 else 90,
                    color=RED if q == mark else SOFT, fontweight='bold' if q == mark else 'normal')

    def cell(x, y, w, h, v, head):
        fc = {0: 'white', 1: GREEN_L, 2: GREEN_D, None: '#EEEEEE'}[v]
        ax.add_patch(Rectangle((x, y), w, h, fc=fc, ec=GREEN if v else NEUTRAL, lw=0.8, zorder=3,
                               hatch='////' if v is None else None))
        ax.text(x + w / 2, y + h / 2, '' if v is None else str(v), ha='center', va='center', fontsize=11,
                color=INK, zorder=4)
        ax.text(x + w / 2, y + h + 0.1, head, ha='center', va='bottom', fontsize=9, color=GREEN_D)

    # (a) three allele groups, one amino acid at position 9
    P = list(range(5, 14)); k9 = P.index(9)
    alleles = ['A*02:01', 'A*02:06', 'A*11:01', 'A*26:01', 'A*24:02']
    x0, y0, dy = 1.05, 3.45, 0.4
    ruler(x0, P, 9, 3.88)
    ax.add_patch(Rectangle((x0 + k9 * W, y0 - (len(alleles) - 1) * dy - 0.2), W, (len(alleles) - 1) * dy + 0.4,
                           fc=BAND, ec='none', zorder=1))
    ys = []
    for i, a in enumerate(alleles):
        y = y0 - i * dy
        sq = aa(a, P)
        ax.text(x0 - 0.1, y, a, ha='right', va='center', fontsize=10, color=INK)
        _letters(ax, x0, y, sq, W, hi=[k9], hicol=GREEN_D if sq[k9] == 'Y' else INK, fs=10.5)
        if sq[k9] == 'Y':
            ys.append(y)
    bx = x0 + len(P) * W + 0.1
    ax.plot([bx, bx + 0.1, bx + 0.1, bx], [ys[0] + 0.15, ys[0] + 0.15, ys[-1] - 0.15, ys[-1] - 0.15], color=GREEN_D,
            lw=1.2)
    ym = (ys[0] + ys[-1]) / 2
    ax.text(bx + 0.2, ym + 0.1, '9 Y', fontsize=12, fontweight='bold', color=GREEN_D, va='center')
    ax.text(bx + 0.2, ym - 0.2, 'one\ncolumn', fontsize=8.5, color=SOFT, va='center', linespacing=1.1)
    ax.text(0.12, 1.12, 'Three allele groups (A*02, A*11, A*26), one amino acid:\nthe residue column 9 Y tests these alleles together,\n'
            'so a signal can be traced to one position.', fontsize=10, color=INK, va='center', linespacing=1.3)
    ax.text(0.12, 0.42, 'A*02:01 and A*02:06 share the first field (A*02), not the amino acid.', fontsize=9.5,
            color=SOFT, va='center')

    # (b) one person: each allele's amino acid, counted over the two copies
    one = ['A*24:02', 'A*11:01']
    xb = 5.35
    ax.text(4.42, 3.98, f'One person: {one[0]} / {one[1]}', fontsize=10.5, color=INK, va='center')
    ax.text(4.42, 3.66, 'amino acids from the IPD-IMGT/HLA protein alignment, not the reads', fontsize=9,
            color=SOFT, va='center')
    ruler(xb, P, 9, 3.36)
    ax.add_patch(Rectangle((xb + k9 * W, 2.56), W, 0.62, fc=BAND, ec='none', zorder=1))
    for i, a in enumerate(one):
        y = 3.03 - i * 0.32
        ax.text(xb - 0.1, y, a, ha='right', va='center', fontsize=10, color=INK)
        _letters(ax, xb, y, aa(a, P), W, hi=[k9], hicol=GREEN_D, fs=10.5)
    arrow(ax, xb + (k9 + 0.5) * W, 2.5, xb + (k9 + 0.5) * W, 2.02)
    res9 = seen(9)
    mine = [aa(a, [9]) for a in one]
    cw = 0.62
    cx = xb + (k9 + 0.5) * W - len(res9) * cw / 2
    for j, r in enumerate(res9):
        cell(cx + j * cw, 1.38, cw - 0.06, 0.42, mine.count(r), f'9 {r}')
    ax.text(cx + len(res9) * cw + 0.08, 1.59, f'= {len(mine)}', fontsize=11, color=INK, va='center')
    ax.text(4.42, 0.72, 'One column per amino acid seen at the position;\neach cell counts the copies that carry it: '
            '0, 1, or 2.', fontsize=10, color=INK, va='center', linespacing=1.3)

    # (c) a person in the cohort whose second allele is not sequenced at 184
    two = ['A*24:02', 'A*24:82']
    Q = list(range(180, 189)); k4 = Q.index(184)
    xc = 9.55
    ax.text(8.62, 3.98, f'One person in the cohort: {two[0]} / {two[1]}', fontsize=10.5, color=INK, va='center')
    ruler(xc, Q, 184, 3.55)
    ax.add_patch(Rectangle((xc + k4 * W, 2.72), W, 0.62, fc=BAND, ec='none', zorder=1))
    got = []
    for i, a in enumerate(two):
        y = 3.19 - i * 0.32
        sq = aa(a, Q)
        got.append(sq[k4])
        ax.text(xc - 0.1, y, a, ha='right', va='center', fontsize=10, color=INK)
        _letters(ax, xc, y, sq, W, hi=[k4], hicol=RED if sq[k4] == '*' else GREEN_D, fs=10.5)
    ax.text(8.62, 2.5, f'* = not sequenced ({two[1]} covers only the groove exons)', fontsize=9,
            color=SOFT, va='center')
    r184 = seen(184)                                   # the amino acids seen at 184 in the cohort
    known = [g for g in got if g != '*']               # the copy whose amino acid the alignment names
    lo = {r: known.count(r) for r in r184}             # the unsequenced copy adds 0 or 1 to one of them
    cw, cx = 0.6, 10.95
    for j, r in enumerate(r184):
        ax.text(cx + j * cw + (cw - 0.06) / 2, 2.2, f'184 {r}', ha='center', va='center', fontsize=9.5, color=GREEN_D)
    rows = [('True count', INK, 'the unsequenced copy is ' + ' or '.join(r184),
             [f'{lo[r]} or {lo[r] + 1}' for r in r184], 'truth'),
            ('Known copy only', RED, 'counts the unsequenced copy as 0', [lo[r] for r in r184], 'wrong'),
            ('Written: blank', GREEN_D, 'unknown, so left out of the test', [None] * len(r184), 'right')]
    for k, (lab, col, why, vals, kind) in enumerate(rows):
        yb = 1.62 - k * 0.6
        ax.text(8.62, yb + 0.29, lab, fontsize=10.5, fontweight='bold', color=col, va='center')
        ax.text(8.62, yb + 0.03, why, fontsize=9, color=col, va='center')
        for j, v in enumerate(vals):
            x = cx + j * cw
            if kind == 'truth':
                ax.add_patch(Rectangle((x, yb), cw - 0.06, 0.38, fc='white', ec=INK, lw=0.8, ls='--', zorder=3))
                ax.text(x + (cw - 0.06) / 2, yb + 0.19, v, ha='center', va='center', fontsize=9, color=INK, zorder=4)
            elif kind == 'wrong':
                ax.add_patch(Rectangle((x, yb), cw - 0.06, 0.38, fc=GREEN_L if v else 'white', ec=RED, lw=1.4,
                                       zorder=3))
                ax.text(x + (cw - 0.06) / 2, yb + 0.19, str(v), ha='center', va='center', fontsize=11, color=INK,
                        zorder=4)
            else:
                ax.add_patch(Rectangle((x, yb), cw - 0.06, 0.38, fc='#EEEEEE', ec=NEUTRAL, lw=0.8, zorder=3,
                                       hatch='////'))
    fig.savefig(OUT / 'concept_residues.png', dpi=DPI)
    plt.close(fig)


def hlahd():
    """How HLA-HD chooses an allele pair, simplified: each read is compared with every named allele; each pair
    is scored by the reads it explains (a read is explained if either allele fits it); the best pair wins."""
    alleles = ['A*24:02', 'A*24:20', 'A*26:01', 'A*11:01']
    fits = [{0, 1, 2, 3}, {0, 1}, {0}, {2}, {2, 3}, {0, 1, 2}, {2}, {0, 3}]      # illustrative
    pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    score = {p: sum(1 for f in fits if p[0] in f or p[1] in f) for p in pairs}
    best = sorted(pairs, key=lambda p: -score[p])
    fig, ax = canvas(12.6, 4.5)
    for k, (x, h) in enumerate(zip((0.2, 4.95, 9.05), ('Compare every read with every allele',
                                                        'Score each pair: reads it explains',
                                                        'Report the best pair'))):
        ax.add_patch(Circle((x + 0.17, 4.2), 0.17, fc=BLUE_D, ec='none'))
        ax.text(x + 0.17, 4.2, str(k + 1), ha='center', va='center', fontsize=11, color='white', fontweight='bold')
        ax.text(x + 0.45, 4.2, h, fontsize=12.5, fontweight='bold', color=INK, va='center')
    # 1: read x allele matrix
    cx0, cw, ry0, rh = 1.35, 0.78, 3.45, 0.34
    for j, a in enumerate(alleles):
        ax.text(cx0 + (j + 0.5) * cw, 3.78, a, ha='center', fontsize=10, color=INK, fontweight='bold')
    for i, f in enumerate(fits):
        y = ry0 - i * rh
        ax.plot([0.3, 0.95], [y, y], color=BLUE, lw=3.2, solid_capstyle='round')
        ax.text(1.05, y, f'read {i + 1}', va='center', fontsize=9.5, color=SOFT)
        for j in range(4):
            on = j in f
            ax.add_patch(Circle((cx0 + (j + 0.5) * cw, y), 0.1, fc=BLUE if on else 'white',
                                ec=BLUE if on else NEUTRAL, lw=1.0, zorder=3))
    ax.text(0.3, 0.42, 'filled dot: the read fits that allele', fontsize=10, color=SOFT)
    arrow(ax, 4.55, 2.3, 4.9, 2.3)
    # 2: pair scores
    for k, p in enumerate(best[:4]):
        y = 3.45 - k * 0.62
        lab = f'{alleles[p[0]]} + {alleles[p[1]]}'
        ax.text(5.05, y, lab, va='center', fontsize=10.5, color=INK, fontweight='bold' if k == 0 else 'normal')
        ax.add_patch(Rectangle((6.95, y - 0.14), 1.5 * score[p] / len(fits), 0.28, color=BLUE if k == 0 else BLUE_L,
                               lw=0))
        ax.text(6.95 + 1.5 * score[p] / len(fits) + 0.08, y, f'{score[p]}/{len(fits)}', va='center', fontsize=10.5,
                color=INK)
    ax.text(5.05, 0.75, 'a pair explains a read if either\nof its alleles fits the read', fontsize=10, color=SOFT,
            va='center')
    arrow(ax, 8.7, 2.3, 9.0, 2.3)
    # 3: result
    p = best[0]
    rbox(ax, 9.1, 2.75, 3.35, 0.9, BLUE_XL, ec=BLUE, lw=1.2)
    ax.text(10.78, 3.2, f'{alleles[p[0]]} / {alleles[p[1]]}', ha='center', va='center', fontsize=15, fontweight='bold',
            color=BLUE_D)
    rbox(ax, 9.1, 1.4, 3.35, 1.05, '#FBEFE3', ec=OCHRE, lw=1.0)
    ax.text(9.25, 2.33, 'If pairs tie, a 2-field name is\nwritten only if the tied candidates\nagree at 2 fields, '
            'HLA-HD\'s\nfrequency table (8 genes) picks\none, or they share one P group.', fontsize=10, color=INK,
            va='top', linespacing=1.25)
    ax.text(9.1, 0.72, 'Simplified: HLA-HD scores each\nexon separately, over all named alleles.', fontsize=10,
            color=SOFT, va='center')
    save(fig, 'concept_hlahd.png')


# =========================================================================== 12 the check
def frequency_check():
    """Why agreement with a population reference checks typing (illustration, invented frequencies)."""
    setup()
    fig, axs = plt.subplots(1, 2, figsize=(12.6, 3.75))
    plt.subplots_adjust(left=0.07, right=0.98, top=0.86, bottom=0.17, wspace=0.28)
    names = ['A*24:02', 'A*02:01', 'A*26:01', 'A*11:01', 'A*31:01', 'A*33:03', 'A*02:06', 'A*24:20']
    ref = np.array([0.36, 0.11, 0.08, 0.09, 0.08, 0.07, 0.09, 0.001])
    good = np.array([0.35, 0.115, 0.083, 0.088, 0.078, 0.071, 0.087, 0.0015])
    bad = good.copy()
    bad[0], bad[-1] = 0.27, 0.08
    for k, (ax, us, ttl, col) in enumerate(((axs[0], good, 'Typing consistent with the population', BLUE),
                                            (axs[1], bad, 'Typing with a systematic error', RED))):
        lim = np.sqrt(0.42)
        ax.plot([0, lim], [0, lim], color=GREY, lw=1, zorder=1)
        ax.scatter(np.sqrt(ref), np.sqrt(us), s=55, color=col, zorder=3, edgecolor='white', lw=0.6)
        ticks = [0.01, 0.05, 0.1, 0.2, 0.4]
        ax.set_xticks(np.sqrt(ticks))
        ax.set_yticks(np.sqrt(ticks))
        ax.set_xticklabels([f'{t:g}' for t in ticks])
        ax.set_yticklabels([f'{t:g}' for t in ticks])
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_xlabel('Allele frequency in the Japanese reference', fontsize=12)
        ax.set_ylabel('Allele frequency in the controls', fontsize=12)
        ax.set_title(ttl, loc='left', fontsize=13, color=col)
        PS.despine(ax)
        ax.set_aspect('equal')
        if k == 0:
            ax.text(np.sqrt(0.12), np.sqrt(0.045), 'on the diagonal', fontsize=11, color=BLUE_D)
        else:
            ax.annotate('too low', (np.sqrt(0.36), np.sqrt(0.27)), (np.sqrt(0.2), np.sqrt(0.13)), fontsize=11, color=RED,
                        ha='center', arrowprops=dict(arrowstyle='-|>', color=RED, lw=1))
            ax.annotate('too high', (np.sqrt(0.001), np.sqrt(0.08)), (np.sqrt(0.02), np.sqrt(0.2)), fontsize=11,
                        color=RED, arrowprops=dict(arrowstyle='-|>', color=RED, lw=1))
    fig.text(0.985, 0.95, 'illustration', ha='right', fontsize=9.5, color=SOFT)
    save(fig, 'concept_frequency_check.png')


# =========================================================================== 13 completeness
def fig_completeness():
    """Per gene: share of the 3,101 people with a 2-field genotype (absent DRB3/4/5 counts as a genotype),
    one field only, no call."""
    parts = S['completeness']['parts_by_gene']
    genes = list(F.ORDER)
    setup()
    fig, ax = plt.subplots(figsize=(6.5, 5.45))
    plt.subplots_adjust(left=0.15, right=0.86, top=0.88, bottom=0.12)
    y = np.arange(len(genes))
    n = {g: sum(parts[g].values()) for g in genes}
    # a 2-field genotype and an absent DRB3/4/5 are both complete results: one colour, so that the truncated
    # axis cannot hide part of a result
    segs = [(('two_field', 'absent'), 'genotype (2-field, or gene absent)', BLUE),
            (('one_field',), 'first field only', RED_L), (('no_call',), 'no call', RED)]
    left = np.zeros(len(genes))
    for keys, lab, col in segs:
        v = np.array([100 * sum(parts[g][k] for k in keys) / n[g] for g in genes])
        ax.barh(y, v, left=left, color=col, height=0.7, label=lab, lw=0)
        left += v
    for i, g in enumerate(genes):
        two = 100 * (parts[g]['two_field'] + parts[g]['absent']) / n[g]
        ax.text(101.2, i, f'{two:.1f}%', va='center', fontsize=10.5, color=INK if two >= 95 else RED,
                fontweight='normal' if two >= 95 else 'bold', clip_on=False)
    ax.set_yticks(y)
    ax.set_yticklabels(genes, fontstyle='italic', fontsize=10.5)
    ax.set_ylim(len(genes) - 0.5, -0.5)
    ax.set_xlim(80, 100)
    ax.set_xticks([80, 85, 90, 95, 100])
    ax.set_xticklabels(['80%', '85%', '90%', '95%', '100%'])
    ax.set_xlabel('People (n = 3,101); axis starts at 80%', fontsize=11.5)
    ax.legend(loc='lower left', bbox_to_anchor=(0, 1.01), ncol=2, frameon=False, fontsize=10, handlelength=1,
              columnspacing=1.0, borderaxespad=0)
    ax.text(1.03, 1.0, 'with a\ngenotype', transform=ax.transAxes, ha='left', va='bottom', fontsize=10, color=SOFT)
    for k, (name, gs) in enumerate(F.GROUPS):
        i0 = genes.index(gs[0])
        if k % 2 == 0:
            ax.axhspan(i0 - 0.5, i0 + len(gs) - 0.5, color='#F2F4F6', zorder=0, lw=0)
    PS.despine(ax, grid_axis='x')
    fig.savefig(OUT / 'fig_completeness.png', dpi=DPI)
    plt.close(fig)


# =========================================================================== 14 frequencies
def fig_panel_coverage():
    """Is each allele called in the controls also carried by ToMMo? Counted two ways, side by side: (a) every
    distinct allele (P group) once, (b) every copy once. Same genes, same three classes; the last row pools the
    compared genes. The rare tail fills much of (a) and almost none of (b)."""
    cov = S['coverage']
    genes = [g for g in F.ORDER if g in cov]
    segs = [('in_panel', 'carried by ToMMo', BLUE),
            ('absent', f'not carried, though named by release {rel_label(S["tommo_vintage"])}', RED),
            ('too_new', f'named after release {rel_label(S["tommo_vintage"])}: cannot be in the panel', GREY)]
    rows = genes + ['All']
    pooled = {key: {k: sum(cov[g][key][k] for g in genes) for k, _l, _c in segs} for key in ('alleles', 'copies')}
    setup()
    fig = plt.figure(figsize=(12.6, 4.6))
    panels = (('alleles', '(a) Distinct alleles', 'Share of the distinct alleles (P-group resolution)', 'alleles'),
              ('copies', '(b) Copies', 'Share of the copies', 'copies'))
    H = 0.66
    for j, (key, ttl, xl, unit) in enumerate(panels):
        ax = fig.add_axes([0.065 + j * 0.355, 0.13, 0.255, 0.70])
        for i, g in enumerate(rows):
            d = pooled[key] if g == 'All' else cov[g][key]
            tot, left = sum(d.values()), 0.0
            for k, lab, col in segs:
                v = 100 * d[k] / tot
                ax.barh(i, v, left=left, color=col, height=H, lw=0, label=lab if (i == 0 and j == 0) else None)
                left += v
            bold = 'bold' if g == 'All' else 'normal'
            ax.text(102, i, f'{tot:,}', va='center', ha='left', fontsize=9.5, color=INK, fontweight=bold,
                    clip_on=False)
            if g == 'All':                                    # the pooled share the conclusion quotes
                ax.text(2, i, f'{100 * d["in_panel"] / tot:.1f}%', va='center', ha='left', fontsize=9.5,
                        color='white', fontweight='bold')
        ax.text(102, -0.95, unit, va='center', ha='left', fontsize=9, color=SOFT, clip_on=False)
        ax.axhline(len(genes) - 0.5, color=INK, lw=0.8)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels(rows if j == 0 else [], fontsize=11)
        for lab in ax.get_yticklabels():
            lab.set_fontstyle('normal' if lab.get_text() == 'All' else 'italic')
            if lab.get_text() == 'All':
                lab.set_fontweight('bold')
        ax.set_ylim(len(rows) - 0.5, -0.5)
        ax.set_xlim(0, 100)
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_xticklabels(['0%', '25%', '50%', '75%', '100%'])
        ax.set_xlabel(xl, fontsize=11)
        ax.set_title(ttl, loc='left', fontsize=12.5, fontweight='bold', color=INK, pad=14)
        for k, (_name, gs) in enumerate(F.GROUPS):
            inside = [g for g in gs if g in genes]
            if inside and k % 2 == 0:
                i0 = genes.index(inside[0])
                ax.axhspan(i0 - 0.5, i0 + len(inside) - 0.5, color='#F2F4F6', zorder=0, lw=0)
        PS.despine(ax, grid_axis='x')
        ax.tick_params(axis='y', length=0)                   # after despine, which restores the tick marks
        if j == 0:
            ax.legend(loc='lower left', bbox_to_anchor=(0, 1.11), ncol=3, frameon=False, fontsize=10,
                      handlelength=1.1, columnspacing=1.4, borderaxespad=0)
    fig.text(0.79, 0.80, 'Not compared', fontsize=11.5, fontweight='bold', color=INK)
    fig.text(0.79, 0.755, 'DRA, DMA, DMB, DOA, DOB, and\nDRB5 have no entry in ToMMo.\n\n'
             'DRB3 and DRB4 are excluded:\nToMMo assigns a DRB3 allele to\nevery haplotype although only\n'
             'a quarter carry the gene, and its\nDRB4 entries do not match its\nown DRB1. Their presence here\n'
             'follows DRB1, so they are also\nexcluded from the association.', fontsize=10, color=SOFT,
             va='top', linespacing=1.3)
    fig.savefig(OUT / 'fig_panel_coverage.png', dpi=DPI)
    plt.close(fig)


def fig_panel_frequencies():
    """One panel per comparable gene: every allele (P group) of that gene, controls vs ToMMo 61KJPN, square-root
    axes, with the gene's Pearson r. Class I on the top row, class II on the bottom row."""
    al = pd.read_csv(C.RESULTS / '03.panels' / 'allele_concordance.tsv', sep='\t')
    loc = pd.read_csv(C.RESULTS / '03.panels' / 'locus_concordance.tsv', sep='\t')
    al = al[al.panel == 'ToMMo']
    loc = loc[loc.panel == 'ToMMo'].set_index('gene')
    comp = [g for g in F.ORDER if g in loc.index and loc.at[g, 'frequency_comparable'] == 1]
    rows = [[g for g in comp if g in CLASS_I], [g for g in comp if g not in CLASS_I]]
    ncol = max(len(r) for r in rows)
    setup()
    fig, axes = plt.subplots(2, ncol, figsize=(12.6, 5.3))
    plt.subplots_adjust(left=0.065, right=0.99, top=0.9, bottom=0.12, wspace=0.42, hspace=0.55)
    for i, genes in enumerate(rows):
        for j in range(ncol):
            ax = axes[i, j]
            if j >= len(genes):
                ax.axis('off')
                if i == 1 and j == len(genes):              # the free slot states what is not compared
                    ax.text(-0.05, 1.0, 'Not compared', transform=ax.transAxes, fontsize=11.5, fontweight='bold',
                            va='top', color=INK)
                    ax.text(-0.05, 0.84, 'DRA, DMA, DMB,\nDOA, DOB, DRB5:\nnot in ToMMo\n\nDRB3, DRB4:\n'
                            'presence follows\nDRB1; ToMMo\nentries unusable', transform=ax.transAxes, fontsize=10,
                            va='top', color=SOFT, linespacing=1.25)
                continue
            g = genes[j]
            d = al[al.gene == g]
            col = BLUE if g in CLASS_I else GREEN
            lim = np.sqrt(max(d.freq_ours.max(), d.freq_panel.max()) * 1.1)
            ax.plot([0, lim], [0, lim], color=GREY, lw=0.9, zorder=1)
            ax.scatter(np.sqrt(d.freq_panel), np.sqrt(d.freq_ours), s=22, color=col, alpha=0.8, lw=0, zorder=3)
            ticks, last = [], -1.0                     # labeled frequencies, far enough apart on the sqrt axis
            for t in (0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8):
                if np.sqrt(t) <= lim and np.sqrt(t) - last >= 0.24 * lim:
                    ticks.append(t)
                    last = np.sqrt(t)
            ax.set_xticks(np.sqrt(ticks))
            ax.set_yticks(np.sqrt(ticks))
            ax.set_xticklabels([f'{t:g}' for t in ticks], fontsize=9.5)
            ax.set_yticklabels([f'{t:g}' for t in ticks], fontsize=9.5)
            ax.set_xlim(0, lim)
            ax.set_ylim(0, lim)
            ax.set_aspect('equal')
            PS.despine(ax)
            ax.set_title(g, loc='left', fontsize=13, fontstyle='italic', color=col, pad=3)
            ax.text(0.97, 0.05, f"r = {loc.at[g, 'r']:.3f}", transform=ax.transAxes, ha='right', va='bottom',
                    fontsize=10.5, color=INK)
            O = S.get('freq_outlier')                     # the one allele far off the diagonal
            if O and O['gene'] == g:
                ax.annotate(O['allele'], (np.sqrt(O['freq_panel']), np.sqrt(O['freq_ours'])),
                            (lim * 0.47, lim * 0.31), fontsize=9.5, color=RED, ha='left', va='center',
                            arrowprops=dict(arrowstyle='-|>', color=RED, lw=1.0, shrinkA=4, shrinkB=3))
    fig.text(0.065, 0.965, 'Class I', fontsize=12, fontweight='bold', color=BLUE)
    fig.text(0.065, 0.505, 'Class II', fontsize=12, fontweight='bold', color=GREEN)
    fig.text(0.53, 0.02, 'Allele frequency in ToMMo 61KJPN (61,424 individuals; square-root scale)', ha='center',
             fontsize=12)
    fig.text(0.008, 0.5, 'Allele frequency in the controls typed at the gene', va='center',
             rotation=90, fontsize=12)
    fig.savefig(OUT / 'fig_panel_frequencies.png', dpi=DPI)
    plt.close(fig)


# =========================================================================== 15 association inputs
def _dosage_box(ax, x, y, w, h, label, value, col, fs=12):
    """One column of the matrix: its name above, its 0/1/2 value inside."""
    ax.text(x + w / 2, y + h + 0.1, label, ha='center', va='bottom', fontsize=9.5, color=SOFT)
    on = isinstance(value, int) and value > 0
    rbox(ax, x, y, w, h, col if on else 'white', ec=col, lw=1.1, r=0.05)
    ax.text(x + w / 2, y + h / 2, str(value), ha='center', va='center', fontsize=fs,
            color='white' if on else INK, fontweight='bold')


def fig_assoc_inputs():
    """A real slice of the matrix the association reads: people × markers, over several genes, with allele
    columns, an 'other' count (a person's remaining copies at the gene; not a column of the matrix) and the
    residue columns of whole positions. Shading is the 0/1/2 dosage; per gene and person the alleles shown plus
    'other' add up to 2, as do the residue columns of each position."""
    M = S['matrix_excerpt']
    G, ppl = M['genes'], M['people']
    setup()
    fig, ax = canvas(12.6, 5.78)
    ax.set_ylim(-1.0, 4.9)
    ax.text(0.1, 4.68, 'The matrix the association reads: full mainland cohort', fontsize=14, fontweight='bold',
            color=INK, va='center')
    ax.text(0.1, 4.36, f'{len(ppl)} of its {M["n_samples"]:,} people × {len(G)} of its {M["n_genes_total"]} genes. '
            'Each column is one marker and each value the number of copies — the encoding the genome-wide scan '
            'uses for a SNP.', fontsize=10.5, color=SOFT, va='center')
    CW, RH, Y1 = 0.2, 0.2, 2.92                        # cell size; Y1 = bottom edge of the top row
    X0, GAP, SUBGAP, PGAP = 1.02, 0.26, 0.1, 0.08      # PGAP separates column groups within allele / residue blocks
    GREY_L, GREY_D = '#C3C8CD', '#6B7279'
    SHADE = {0: ('white', NEUTRAL), 1: (BLUE_L, BLUE), 2: (BLUE_D, BLUE_D)}
    SHADE_R = {0: ('white', NEUTRAL), 1: (GREEN_L, GREEN), 2: (GREEN_D, GREEN_D)}
    SHADE_O = {0: ('white', NEUTRAL), 1: (GREY_L, GREY), 2: (GREY_D, GREY_D)}
    yl = Y1 - len(ppl) * RH + 0.05

    def column(x, name, vals, shade, style='normal'):
        ax.text(x + CW / 2, Y1 + RH + 0.06, name, fontsize=7, color=SOFT, rotation=90, ha='center', va='bottom',
                fontstyle=style)
        for i, v in enumerate(vals):
            fc, ec = shade[v]
            ax.add_patch(Rectangle((x + 0.012, Y1 - i * RH + 0.012), CW - 0.024, RH - 0.024, fc=fc, ec=ec, lw=0.5))

    x = X0
    spans = []
    for gi, g in enumerate(G):
        gx0 = x
        # alleles, then 'other' after a small gap
        for k, nm in enumerate(g['alleles']):
            column(x, nm, [p['alleles'][gi][k] for p in ppl], SHADE)
            x += CW
        ax.plot([gx0, x - 0.02], [yl] * 2, color=BLUE, lw=1.4)
        ax.text((gx0 + x) / 2, Y1 - len(ppl) * RH - 0.12, 'alleles', fontsize=8.5, color=BLUE, ha='center',
                va='top')
        x += PGAP
        column(x, 'other', [p['other'][gi] for p in ppl], SHADE_O, 'italic')
        ax.plot([x, x + CW - 0.02], [yl] * 2, color=GREY, lw=1.4)
        x += CW + SUBGAP
        # residues, one group per position
        rx0, k = x, 0
        for n in g['residue_groups']:
            for _ in range(n):
                column(x, g['residues'][k], [p['residues'][gi][k] for p in ppl], SHADE_R)
                k += 1
                x += CW
            ax.plot([x - n * CW, x - 0.02], [yl] * 2, color=GREEN, lw=1.4)
            x += PGAP
        x -= PGAP
        ax.text((rx0 + x) / 2, Y1 - len(ppl) * RH - 0.12, 'residues', fontsize=8.5, color=GREEN, ha='center',
                va='top')
        x += GAP
        spans.append((gx0, x - GAP, g))
    for gx0, gx1, g in spans:
        ax.text((gx0 + gx1) / 2, Y1 + 1.0, f'HLA-{g["gene"]}', fontsize=12, fontweight='bold', color=INK,
                ha='center', fontstyle='italic')
        ax.text((gx0 + gx1) / 2, Y1 + 0.85, f'{g["n_alleles"]} allele and {g["n_residues"]} residue columns tested',
                fontsize=8.5, color=SOFT, ha='center')
        ax.plot([gx0, gx1 - 0.02], [Y1 + 0.79] * 2, color=NEUTRAL, lw=1.2)
    for i in range(len(ppl)):
        ax.text(X0 - 0.12, Y1 - i * RH + RH / 2, f'Person {i + 1}', fontsize=8.5, color=SOFT, ha='right',
                va='center')
    # legend and caption
    ly = Y1 - len(ppl) * RH - 0.52
    ax.text(X0 - 0.12, ly, 'Copies:', fontsize=10, color=INK, ha='right', va='center')
    bx = X0
    for shade, name in ((SHADE, 'allele column'), (SHADE_R, 'residue column'), (SHADE_O, '"other" (not a column)')):
        for v in (0, 1, 2):
            fc, ec = shade[v]
            ax.add_patch(Rectangle((bx, ly - 0.08), 0.18, 0.18, fc=fc, ec=ec, lw=0.8))
            ax.text(bx + 0.22, ly, str(v), fontsize=10, color=INK, va='center')
            bx += 0.42
        ax.text(bx - 0.04, ly, name, fontsize=10, color=SOFT, va='center')
        bx += 1.75
    ax.text(X0 - 0.12, ly - 0.4, 'Real data. Column names: an allele column is the 2-field allele; a residue column '
            'is the position and the amino acid, so "9 F" is phenylalanine at position 9.', fontsize=9.5, color=SOFT,
            va='center')
    ax.text(X0 - 0.12, ly - 0.68, f'Shown per gene: the {len(G[0]["alleles"])} most frequent of the allele columns '
            'tested, and two positions with every amino acid seen there as a column (a gap separates them).',
            fontsize=9.5, color=SOFT, va='center')
    ax.text(X0 - 0.12, ly - 0.96, '"other" is not a column of the matrix: it counts the person\'s remaining copies, of '
            'alleles not shown here or of alleles below the MAF floor, which have no column.', fontsize=9.5,
            color=SOFT, va='center')
    ax.text(X0 - 0.12, ly - 1.24, 'Per gene, each person therefore adds up to 2 copies: over the alleles shown plus '
            '"other", and over the residue columns at each position (one amino acid per copy).', fontsize=9.5,
            color=SOFT, va='center')
    save(fig, 'fig_assoc_inputs.png')


def fig_outlier_panel():
    """The DRB1 panel of the frequency slide, enlarged: every DRB1 allele of either table (square-root axes),
    with the dot off the diagonal and its two DRB1*15 neighbours named."""
    O = S['freq_outlier']
    al = pd.read_csv(C.RESULTS / '03.panels' / 'allele_concordance.tsv', sep='\t')
    d = al[(al.panel == 'ToMMo') & (al.gene == O['gene'])]
    named = {O['allele']: RED, O['neighbours'][0]['allele']: GREEN_D, O['neighbours'][1]['allele']: GREEN_D}
    setup()
    fig, ax = plt.subplots(figsize=(4.6, 4.1))
    plt.subplots_adjust(left=0.17, right=0.96, top=0.9, bottom=0.15)
    lim = np.sqrt(max(d.freq_ours.max(), d.freq_panel.max()) * 1.1)
    ax.plot([0, lim], [0, lim], color=GREY, lw=0.8, zorder=1)
    rest = d[~d.allele.isin(named)]
    ax.scatter(np.sqrt(rest.freq_panel), np.sqrt(rest.freq_ours), s=22, color=GREEN, alpha=0.75, lw=0, zorder=2)
    # label places in empty parts of the panel (fractions of the axis): the outlier above it, the 15:02 dot to
    # the upper left, the 15:01 dot to the lower right
    spots = {O['allele']: (0.62, 0.22, 'center'), O['neighbours'][0]['allele']: (0.42, 0.95, 'right'),
             O['neighbours'][1]['allele']: (0.99, 0.5, 'right')}
    for a, c in named.items():
        r = d[d.allele == a].iloc[0]
        x, y = np.sqrt(r.freq_panel), np.sqrt(r.freq_ours)
        ax.scatter([x], [y], s=46, color=c, lw=0, zorder=3)
        fx, fy, ha = spots[a]
        ax.annotate(f'{a}\nToMMo {100 * r.freq_panel:.1f}%, ours {100 * r.freq_ours:.1f}%', (x, y),
                    (fx * lim, fy * lim), ha=ha, va='center', fontsize=8.5, color=c,
                    arrowprops=dict(arrowstyle='-', color=c, lw=0.7))
    ticks = [0, 0.01, 0.05, 0.1, 0.2]
    ax.set_xticks(np.sqrt(ticks)); ax.set_yticks(np.sqrt(ticks))
    ax.set_xticklabels([f'{x:g}' for x in ticks], fontsize=9); ax.set_yticklabels([f'{x:g}' for x in ticks], fontsize=9)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel('Frequency in ToMMo 61KJPN', fontsize=10)
    ax.set_ylabel('Frequency in our controls', fontsize=10)
    r_ = S['tommo'][O['gene']]['r'] if 'r' in S['tommo'][O['gene']] else None
    ax.set_title(f'{O["gene"]} panel of the frequency figure, enlarged' + (f' (r = {r_:.3f})' if r_ else ''), loc='left',
                 fontsize=10.5, fontweight='bold')
    PS.despine(ax, grid_axis='both')
    save(fig, 'fig_outlier_panel.png')


def fig_outlier_group():
    """The DRB1*15 group on both sides (stats 'freq_outlier'): where ToMMo puts 2.3% of its DRB1 copies on
    DRB1*15:13, our controls carry none, yet the group totals are close; the split inside the group differs.
    Shares of each side's DRB1 copies."""
    O = S['freq_outlier']
    nb = {n['allele']: n for n in O['neighbours']}
    a1, a2 = O['neighbours'][0]['allele'], O['neighbours'][1]['allele']          # DRB1*15:02, DRB1*15:01
    sides = [('ToMMo 61KJPN', 'panel', O['freq_panel']), ('Our controls', 'ours', O['freq_ours'])]
    segs = [(a1, GREEN_D), (a2, GREEN_L), (O['allele'], RED), ('other', GREY)]
    setup()
    fig, ax = plt.subplots(figsize=(5.0, 1.9))
    plt.subplots_adjust(left=0.24, right=0.86, top=0.70, bottom=0.27)
    for i, (lab, key, f_out) in enumerate(sides):
        vals = [nb[a1][key], nb[a2][key], f_out]
        vals.append(max(O['group'][key] - sum(vals), 0))
        left = 0
        for (name, col), v in zip(segs, vals):
            ax.barh(i, 100 * v, left=left, height=0.56, color=col, lw=0)
            if 100 * v >= 2:
                ax.text(left + 50 * v, i, f'{100 * v:.1f}', ha='center', va='center', fontsize=9, color='white',
                        fontweight='bold')
            left += 100 * v
        ax.text(left + 0.4, i, f'{100 * O["group"][key]:.1f}%', va='center', fontsize=9.5, color=INK)
    ax.set_yticks([0, 1]); ax.set_yticklabels([s[0] for s in sides], fontsize=9.5)
    ax.set_ylim(1.6, -0.6)
    ax.set_xlim(0, 22)
    ax.set_xlabel('% of DRB1 copies', fontsize=9.5)
    ax.tick_params(axis='x', labelsize=8.5)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c in segs]
    names = [n.split('*')[1] if '*' in n else 'other DRB1*15' for n, _ in segs]
    ax.legend(handles, names, loc='lower left', bbox_to_anchor=(-0.02, 1.02), ncol=4, frameon=False, fontsize=8.5,
              handlelength=1, columnspacing=0.9, borderaxespad=0)
    ax.set_title('The DRB1*15 group', loc='left', fontsize=10.5, fontweight='bold', pad=22)
    PS.despine(ax, grid_axis='x')
    ax.tick_params(axis='y', length=0)
    save(fig, 'fig_outlier_group.png')


def fig_worked_example():
    """One gene through the ToMMo comparison (stats 'worked_example'): every allele at P-group resolution found in
    either table, drawn as on the frequency slide (square-root axes), colored by where it is found. Alleles one table
    lacks sit at 0 on that axis; an inset magnifies the corner where they lie, with how many share each point."""
    W = S['worked_example']
    rows = W['rows']
    col = {'both': GREEN, 'ours': OCHRE, 'tommo': GREY}
    lab = {'both': f'in both tables ({sum(r["where"] == "both" for r in rows)})',
           'ours': f'only in our controls ({sum(r["where"] == "ours" for r in rows)})',
           'tommo': f'only in ToMMo ({sum(r["where"] == "tommo" for r in rows)})'}
    setup()
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    plt.subplots_adjust(left=0.15, right=0.97, top=0.86, bottom=0.14)
    lim = np.sqrt(max(max(r['freq_ours'], r['freq_tommo']) for r in rows) * 1.12)
    ax.plot([0, lim], [0, lim], color=GREY, lw=0.8, zorder=1)
    for w in ('both', 'ours', 'tommo'):
        d = [r for r in rows if r['where'] == w]
        ax.scatter([np.sqrt(r['freq_tommo']) for r in d], [np.sqrt(r['freq_ours']) for r in d], s=34, color=col[w],
                   alpha=0.85, lw=0, zorder=3, label=lab[w], clip_on=False)
    ticks = [0, 0.01, 0.05, 0.1, 0.2, 0.4]
    ax.set_xticks(np.sqrt(ticks)); ax.set_yticks(np.sqrt(ticks))
    ax.set_xticklabels([f'{x:g}' for x in ticks], fontsize=9.5); ax.set_yticklabels([f'{x:g}' for x in ticks], fontsize=9.5)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel('Frequency in ToMMo 61KJPN', fontsize=10.5)
    ax.set_ylabel('Frequency in our controls', fontsize=10.5)
    ax.set_title(f'{W["gene"]}: {W["either"]} dots, r = {W["r"]:.3f}', loc='left', fontsize=12, fontweight='bold',
                 fontstyle='normal', pad=26)
    ax.legend(loc='lower left', bbox_to_anchor=(0, 1.0), ncol=3, frameon=False, fontsize=8.5, handletextpad=0.1,
              columnspacing=0.8, borderaxespad=0.2)
    PS.despine(ax, grid_axis='both')
    # the corner where the one-sided alleles lie, magnified (linear, in % of copies)
    ins = ax.inset_axes([0.08, 0.56, 0.32, 0.32])           # the empty upper left, clear of every dot
    for w in ('ours', 'tommo'):
        d = [r for r in rows if r['where'] == w]
        xs = [100 * r['freq_tommo'] for r in d]
        ys = [100 * r['freq_ours'] for r in d]
        ins.scatter(xs, ys, s=26, color=col[w], lw=0, zorder=3, clip_on=False)
        for (x, y), n in pd.Series(list(zip(xs, ys))).value_counts().items():
            if n > 1:
                ins.text(x + 0.00015, y, f'×{n}', fontsize=8.5, color=col[w], va='center')
    ins.set_xlim(0, 0.004); ins.set_ylim(0, 0.08)
    ins.set_xticks([0, 0.002, 0.004]); ins.set_yticks([0, 0.04, 0.08])
    ins.tick_params(labelsize=8)
    ins.set_title('near 0, in % of copies', fontsize=8.5, color=SOFT, pad=3)
    for sp in ins.spines.values():
        sp.set_color(NEUTRAL)
    save(fig, 'fig_worked_example.png')


def broken_panel(hi, lo, x, series, hi_lim, lo_lim, label=None):
    """One column of grouped bars drawn across a broken y-axis: `hi` carries the one tall class, `lo` the range
    a floor acts in, which would otherwise be flattened against it. `series` is [(legend label, color, values)];
    only the cohort named in `label` is given printed values, so the classes stay readable."""
    for ax in (hi, lo):
        for k, (lab, col, v) in enumerate(series):
            ax.bar(x + (k - 1) * 0.27, v, width=0.26, color=col, label=lab)
        ax.grid(axis='x', visible=False)
    hi.set_ylim(*hi_lim)
    lo.set_ylim(*lo_lim)
    hi.tick_params(bottom=False, labelbottom=False)
    # the break marks: drawn in figure inches, not in each axes' own fraction, so that the two halves — which
    # have different heights — get marks of the same length and the same slope
    fig = hi.figure
    fw, fh = fig.get_size_inches()
    dx, dy = 0.075, 0.055
    for ax, edge in ((hi, 'y0'), (lo, 'y1')):
        bb = ax.get_position()
        bx, by = bb.x0, getattr(bb, edge)
        fig.add_artist(plt.Line2D([bx - dx / fw, bx + dx / fw], [by - dy / fh, by + dy / fh],
                                  transform=fig.transFigure, color=INK, lw=1))
    for k, (lab, col, v) in enumerate(series):
        if label is not None and lab != label:
            continue
        for xi, vi in zip(x, v):
            ax = hi if vi >= hi_lim[0] else lo
            span = (hi_lim[1] - hi_lim[0]) if ax is hi else (lo_lim[1] - lo_lim[0])
            ax.text(xi + (k - 1) * 0.27, vi + span * 0.035, f'{vi:,}', ha='center', va='bottom', fontsize=9,
                    color=col)
    for ax in (hi, lo):
        PS.despine(ax, grid_axis='y')
    hi.spines['bottom'].set_visible(False)
    lo.spines['top'].set_visible(False)


def limits(vals, tall):
    """y limits for the two halves of a broken axis: `tall` is the class drawn on the upper half."""
    top = [v[tall] for v in vals]
    rest = [v[i] for v in vals for i in range(len(v)) if i != tall]
    return (min(top) * 0.93, max(top) * 1.10), (0, max(rest) * 1.32 or 1)


def fig_mac_distribution():
    """MAF of every marker the call-rate floor leaves: markers per MAF class, alleles and residues, three
    cohorts (colors of the cohort slides), on a broken y-axis so the classes the floor keeps stay readable
    beside the <0.1% class. Markers absent from a cohort are not drawn."""
    B = S['assoc_basis']['cohorts']
    co = [('narrow_mainland', 'narrow mainland', '#0571B0'),
          ('intermediate_mainland', 'intermediate mainland', '#008837'),
          ('full_mainland', 'full mainland', '#B2182B')]
    labels = list(B['full_mainland']['distribution']['allele']['bins'])
    tested = [i for i, t in enumerate(labels) if t not in ('<0.1%', '0.1-0.5%')]      # kept at MAF >= 0.5 %
    setup()
    fig, axes = plt.subplots(2, 2, figsize=(12.6, 4.15), sharex='col',
                             gridspec_kw={'height_ratios': [1, 2.15], 'hspace': 0.09, 'wspace': 0.17})
    plt.subplots_adjust(left=0.068, right=0.99, top=0.775, bottom=0.155)
    fig.text(0.068, 0.985, 'Cohort (genomes):', fontsize=10.5, color=SOFT, va='top')
    x = np.arange(len(labels))
    for j, (cl, ttl) in enumerate((('allele', 'Allele columns'), ('residue', 'Residue columns'))):
        hi, lo = axes[0, j], axes[1, j]
        series = [(f'{lab} ({B[key]["n"]:,})', col, [B[key]['distribution'][cl]['bins'][t] for t in labels])
                  for key, lab, col in co]
        hl, ll = limits([v for _, _, v in series], labels.index('<0.1%'))
        broken_panel(hi, lo, x, series, hl, ll, label=f'full mainland ({B["full_mainland"]["n"]:,})')
        lo.axvspan(min(tested) - 0.5, len(labels) - 0.4, color='#F2F4F6', zorder=0, lw=0)
        lo.text(len(labels) - 0.45, ll[1] * 0.97, 'tested: MAF ≥ 0.5%', fontsize=10, color=SOFT, va='top',
                ha='right')
        lo.set_xticks(x)
        lo.set_xticklabels([t.replace('>=', '≥') for t in labels], fontsize=10.5)
        lo.set_xlabel('MAF: copies of the less common form, over the copies typed', fontsize=11.5)
        n_cl = B['full_mainland']['distribution'][cl]['live_after_cr']
        hi.set_title(f'{ttl} ({n_cl:,} after call rate, full mainland)', loc='left', fontsize=13, pad=8)
    fig.text(0.068, 0.885, 'The 16 genes of the association, after call rate ≥ 95% in cases and in controls; DRB3, '
             'DRB4, and DRB5 are left out', fontsize=10, color=SOFT, va='top')
    axes[1, 0].set_ylabel('Markers', fontsize=12)
    axes[1, 0].yaxis.set_label_coords(-0.085, 0.72)
    h, l_ = axes[1, 0].get_legend_handles_labels()
    fig.legend(h, l_, loc='upper left', ncol=3, frameon=False, fontsize=10.5, bbox_to_anchor=(0.158, 1.025))
    fig.savefig(OUT / 'fig_mac_distribution.png', dpi=DPI)
    plt.close(fig)


def fig_call_rate():
    """Call rate of every marker, alleles and residues apart, three cohorts. Broken y-axis: the ">= 0.99" class
    holds most markers, so the low classes the floor acts on would otherwise be invisible. Split by kind because
    the floor acts on one kind only: no allele of any cohort falls below 0.95."""
    B = S['assoc_basis']['cohorts']
    co = [('narrow_mainland', 'narrow mainland', '#0571B0'),
          ('intermediate_mainland', 'intermediate mainland', '#008837'),
          ('full_mainland', 'full mainland', '#B2182B')]
    labels = ['<0.90', '0.90-0.95', '0.95-0.98', '0.98-0.99', '>=0.99']
    setup()
    fig, axes = plt.subplots(2, 2, figsize=(12.6, 4.15), sharex='col',
                             gridspec_kw={'height_ratios': [1, 2.15], 'hspace': 0.09, 'wspace': 0.17})
    plt.subplots_adjust(left=0.068, right=0.99, top=0.775, bottom=0.155)
    fig.text(0.068, 0.985, 'Cohort (genomes):', fontsize=10.5, color=SOFT, va='top')
    x = np.arange(len(labels))
    for j, (cl, ttl) in enumerate((('allele', 'Allele columns'), ('residue', 'Residue columns'))):
        hi, lo = axes[0, j], axes[1, j]
        series = [(f'{lab} ({B[key]["n"]:,})', col,
                   [B[key]['distribution'][cl]['call_rate_bins'][t] for t in labels]) for key, lab, col in co]
        hl, ll = limits([v for _, _, v in series], labels.index('>=0.99'))
        broken_panel(hi, lo, x, series, hl, ll, label=f'full mainland ({B["full_mainland"]["n"]:,})')
        lo.axvspan(-0.6, 1.5, color='#FBE9E7', zorder=0, lw=0)
        lost = B['full_mainland']['distribution'][cl]['cr_dropped']['0.95']
        lo.text(-0.52, ll[1] * 0.96, f'removed at call rate ≥ 95%: {lost:,}', fontsize=10.5, color=RED, va='top')
        lo.set_xticks(x)
        lo.set_xticklabels([t.replace('>=', '≥') for t in labels], fontsize=11)
        lo.set_xlabel('Call rate in the less complete group (cases or controls)', fontsize=11.5)
        n_cl = B['full_mainland']['distribution'][cl]['live']
        hi.set_title(f'{ttl} ({n_cl:,} varying, full mainland)', loc='left', fontsize=13, pad=8)
    fig.text(0.068, 0.885, 'The 16 genes of the association; DRB3, DRB4, and DRB5 are left out. Printed numbers '
             'are the full mainland cohort.', fontsize=10, color=SOFT, va='top')
    axes[1, 0].set_ylabel('Markers', fontsize=12)
    axes[1, 0].yaxis.set_label_coords(-0.085, 0.72)
    h, l_ = axes[1, 0].get_legend_handles_labels()
    fig.legend(h, l_, loc='upper left', ncol=3, frameon=False, fontsize=10.5, bbox_to_anchor=(0.158, 1.025))
    fig.savefig(OUT / 'fig_call_rate.png', dpi=DPI)
    plt.close(fig)


if __name__ == '__main__':
    out = {}
    hla_function()
    out['structure'] = hla_structure()
    mhc_map()
    inheritance()
    out['diversity'] = diversity()
    prior_cteph()
    nomenclature()
    dr_haplotypes()
    out['reads'] = short_reads()
    hlahd()
    residue_concept()
    frequency_check()
    fig_completeness()
    fig_panel_coverage()
    fig_panel_frequencies()
    fig_outlier_group()
    fig_outlier_panel()
    fig_worked_example()
    fig_assoc_inputs()
    fig_mac_distribution()
    fig_call_rate()
    (HERE / 'figure_facts.json').write_text(json.dumps(out, indent=1, default=int) + '\n')
    print('figures written to', OUT)
