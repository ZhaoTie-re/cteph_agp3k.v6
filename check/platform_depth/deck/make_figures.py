#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The figures of the summary deck that the component's figure families cannot supply at slide size. Each
#           is drawn at the size of its slot from the same result tables, with the component's encodings and the
#           shared project style (analysis/_shared/scripts/plot_style.py). Why each is redrawn:
#             forest_platform.full_mainland, forest_sources.full_mainland  results/figures/03.platform/forest.<cohort> and
#                 05.controls/sources.<cohort>: their legend overlaps the panel title, and at slide size their text is
#                 about 6 pt. Same estimates and colors; the Q P printed is the one the evidence rule reads
#                 (platforms: empirical P against the null SNPs; sources: Q P). The common markers (MAF >= 5 %, where the
#                 rule applies) and the rare ones form two blocks in the matrix order; the rare rows are pale except the
#                 points the slides name (INPP5D; ADARB2 against ACC); each group has a shape as well as a color.
#             freq_panel (not drawn since the flag slides were restructured, user 2026-11-06)  06.evidence/
#                 lead_frequencies.<cohort>: one panel there is 1.6 in wide; cases by platform only, not by registration
#                 study (user 2026-11-05).
#             cram_ADARB2, cram_LINC01146  the two flagged leads (deck/cram_check, stats.json cram, linc2): ADARB2, how a
#                 carrier is read (the minisatellite drawn to scale in bp), carriers called against in the reads, the
#                 signature reads per genome; LINC01146, each group's A frequency in the calls and in the reads, and
#                 each case platform's odds ratio with the strata combined. Full mainland.
#             carriers_<gene>  06.evidence/gene_carriers.<cohort>: a crop loses the group labels (they sit on the first
#                 column only); the driver panels are new (followup.gene_carriers.tsv, site rows); full mainland. JAK2:
#                 called carriers for the gene set and one read-based count for V617F, no "(+n)" (the fixed decision).
#             evidence_matrix  06.evidence/evidence_matrix: its text is about 6 pt at slide size; same cells and colors;
#                 a flag: bold navy border and "!", its row on a SIG band; the gene sets' arrays: one cell with its count; N I F over each check.
#             lead_frequencies_cohorts.{1,2}, gene_carriers_cohorts  06.evidence/{lead_frequencies,
#                 gene_carriers}.<cohort> side by side: the three cohorts in one panel per lead or gene set (appendix;
#                 the eight gene sets on one slide, M/H on top, L/M/H below).
#             loo_regions  the regions of the leave-out plane (flag, below the bar, not testable, consistent) with the
#                 spread an unchanged effect would have, at one real reading (ADARB2 without DNBSeq-T7, narrow).
#             loo_scatter.{platform,sources}  06.evidence/z_observed_expected (one check per slide, Q2 and Q4): its labels are 4.5 pt and overlap in panel b. Same points
#                 (evidence.loo.tsv) and the same flag rule (scripts/checks.py), drawn at slot size; flagged points
#                 labeled with offsets that do not collide; the readings below the testability bar with |z| < 1.96 are
#                 ringed and named on a gray band (stats.json loo_below_bar).
#             array_panels.full_mainland  04.array/array_evidence.<cohort>: its text is 5-5.5 pt at slide size and point labels
#                 overlap. Same tables (array.evidence.tsv); panel a: one named row per signal, WGS open and array filled.
#             v57_trace  06.evidence/v57_replicates: its cell text is 4.2 pt (under 4 pt at slide width). Same table
#                 (followup.hla_v57.tsv, full mainland), cells and colors, the arm plan's arms only; drawn at its slot
#                 in the DQB1:57:V dossier (under the verdict), the cells' text naming a missing dosage, so no legend;
#                 not drawn when no down-sampled case changes its 57 V dosage. A case is named by its platform, never
#                 by its sample ID.
#             Drawn for the deck's grammar (tmp/ppt/hla_assoc), in its palette (cohorts narrow / intermediate / full;
#             cases dark red, controls gray, reference panels black diamonds; SIG marks a flag):
#             fig_q1..fig_q4  the four questions with hypothetical values (a real signal against an artifact), at the
#                 size of their cell on the aims slide, so every text is 9 pt.
#             background_depth  depth (median, IQR) and read length of the five sequencing groups, full mainland
#                 (01.depth/platform_depth.tsv).
#             emulation_dumbbell  Q1's emulation as before -> after against the controls, per case group (full
#                 mainland): a, effective depth (01.depth/depth_per_sample.tsv; after = fraction kept x depth from
#                 fractions.all.tsv, achieved depth checked by gate G8); b, genotypes failing a QC rule (stats.json qc_fail).
#             ds_zscale  Q1's signals on the |z| scale, one panel per cohort: the original |z| and the median of the
#                 replicates at the controls' depth, every signal (gene sets: the test that called the set); no thresholds,
#                 no ranges; the drop must equal the matrix reading.
#             v57_pg  02.downsampling/hla_dqb1_57v with the three cohorts in their colors and text at slide size
#                 (ds.hla_dqb1_57v.tsv).
#             depth_platforms  results/figures/01.depth/depth_platforms in the deck's grammar (the family figure colors
#                 two platforms in the cohort colors): per-genome depth, fraction kept, read length (01.depth/).
# Inputs  : the results root and the report tables as collect_stats.py reads them (PD_RESULTS, PD_TABLES)
# Outputs : figures/*.png
# Usage   : /home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3 make_figures.py
# ---------------------------------------------------------------------------
import json
import math
import os
import re
import shutil
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402,F401
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from scipy.stats import norm  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'figures')
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(HERE, '..', 'scripts'))
import design_figures as DF  # noqa: E402  (PD_RESULTS / PD_TABLES)
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'analysis', '_shared', 'scripts'))
import plot_style as S  # noqa: E402

R = str(DF.RES)
CK = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']
CNAME = {'narrow_mainland': 'Narrow mainland', 'intermediate_mainland': 'Intermediate mainland',
         'full_mainland': 'Full mainland'}
CSHORT = {'narrow_mainland': 'N', 'intermediate_mainland': 'I', 'full_mainland': 'F'}
STRATA = [('only_t7', 'DNBSeq-T7'), ('only_novaseq', 'NovaSeq'), ('only_g400rs30', 'DNBSeq-G400RS 30x'),
          ('only_g400rs15', 'DNBSeq-G400RS 15x')]
SOURCES = [('nagahama', 'Nagahama'), ('bbj', 'BBJ'), ('htlv1', 'HTLV-1'), ('acc', 'ACC'), ('lungca', 'Lung cancer')]
# the deck palette of tmp/ppt/hla_assoc: cohorts narrow / intermediate / full; groups (cases dark red, controls grey,
# reference panels black diamonds); platforms and sources (hla_assoc's platform palette); SIG marks a flag
NAVY, INK_, GREY_, SIG = '#1F3864', '#212121', '#595959', '#FCEFD2'
INK_SOFT_ = '#454A50'
PALE_ = '#F4F7FB'
COH_COL = {'narrow_mainland': '#0571B0', 'intermediate_mainland': '#008837', 'full_mainland': '#B2182B'}
CASE_C, CTRL_C, REF_C = '#8E1B22', '#8C8C8C', S.INK
PLAT_COL = ['#1F3864', '#4F81BD', '#8DB0DC', '#B07C1E']               # DNBSeq-T7, NovaSeq, G400RS 30x, G400RS 15x
SRC_COL = ['#1F3864', '#4F81BD', '#8DB0DC', '#B07C1E', '#7F7F7F']    # Nagahama, BBJ, HTLV-1, ACC, lung cancer
OCHRE = '#B07C1E'


def ptitle(ax, letter, text, size=10.5):
    # a panel letter (bold) and a short centered title, as in tmp/ppt/hla_assoc
    ax.set_title(rf'$\bf{{{letter}}}$  {text}', fontsize=size, loc='center')


# review round 9 (D-06): gene symbols in italics, digits too. Mathtext's \it and \mathit keep digits upright
# ("STBD" + "1"), so a gene label is drawn as two artists: the gene in italic text, the rest upright.
def gene_label(ax, x, y, gene, rest, **kw):
    """the gene in italics, the rest upright; right-aligned at x (data coordinates); call once the axes are placed"""
    t = ax.text(x, y, rest, ha='right', **kw)
    bb = t.get_window_extent(ax.figure.canvas.get_renderer())
    gap = 0.28 * kw.get('fontsize', 10) * ax.figure.dpi / 72                          # one space
    xg = ax.transData.inverted().transform((bb.x0 - gap, 0))[0]
    return ax.text(xg, y, gene, ha='right', style='italic', **kw)


def gene_ticklabels(ax, labels, fontsize, color=None, pad_pt=4.0):
    """y tick labels whose first word is a gene symbol (their tick text left empty): labels = [(y, gene, rest)]; the
    gene in italics, the rest upright, packed and right-aligned at the axis"""
    from matplotlib.offsetbox import TextArea, HPacker, AnnotationBbox
    color = color or INK_
    for y, gene, rest in labels:
        kids = [TextArea(gene, textprops=dict(fontsize=fontsize, style='italic', color=color))]
        if rest:
            kids.append(TextArea(rest, textprops=dict(fontsize=fontsize, color=color)))
        box = HPacker(children=kids, align='baseline', pad=0, sep=0.28 * fontsize)
        ax.add_artist(AnnotationBbox(box, (0, y), xycoords=ax.get_yaxis_transform(), xybox=(-pad_pt, 0),
                                     boxcoords='offset points', box_alignment=(1, 0.5), frameon=False, pad=0,
                                     annotation_clip=False))
COHORT_MK = {'narrow_mainland': 'o', 'intermediate_mainland': 's', 'full_mainland': '^'}
COHORT_OFF = {'narrow_mainland': 0.24, 'intermediate_mainland': 0.0, 'full_mainland': -0.24}


def tsv(*p):
    return pd.read_csv(os.path.join(R, *p), sep='\t')


def rtsv(name):
    return pd.read_csv(DF.report_file(name), sep='\t')


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, f'{name}.png'), dpi=300, bbox_inches='tight', pad_inches=0.04)
    plt.close(fig)


# ── forests: the case-platform strata and the control sources ──────────────
def deck_signals(df):
    """the evidence rows the deck shows: the HLA residue position (DQB1 position 57) is left out (user 2026-11-05:
    redundant with DQB1:57:V); the pipeline still tests and reports it"""
    return df[df.kind != 'hla_position'] if 'kind' in df.columns else df[~df.signal.str.contains('position')]


def _matrix_order():
    # the signals in the order of the evidence matrix (evidence.long.tsv), used by every per-signal figure
    return list(dict.fromkeys(deck_signals(rtsv('evidence.long.tsv')).signal))


# rows the slide text names are drawn at full strength even when rare (platform: INPP5D; sources: ADARB2 vs ACC)
# review round 12 (D-10): slide 24 names no rare lead, so ADARB2 / ACC is faded like the other rare leads
FOREST_NAMED = {'platform': {'INPP5D': None}, 'sources': {}}
PLAT_MK = ['o', 's', '^', 'v']
SRC_MK = ['o', 's', '^', 'v', 'h']


# review round 13 (F-01): drawn at the slot's size (placed at 0.985 of its drawn size, was 0.92), intervals 1.4 pt,
# every estimate drawn above the pooled diamond (white edge)
def forest(kind, cohort, w=5.95, h=5.15):
    psig = tsv('10.platform', '_summary', 'platform.signals.tsv')
    rows = psig[(psig.cohort == 'full_mainland') & psig.role.isin(['Lead', 'HLA residue'])]
    mk = dict(zip(rows.signal, rows.marker))
    het_p = tsv('10.platform', cohort, 'platform.heterogeneity.tsv').set_index('variant_id')
    if kind == 'platform':
        st = tsv('10.platform', cohort, 'platform.strata.tsv')
        groups = [(k, v.replace('DNBSeq-G400RS', 'G400RS')) for k, v in STRATA]
        col_q, het, series, shapes = 'p_Q_empirical', het_p, PLAT_COL, PLAT_MK
        xlabel, pooled = 'OR (95% CI), one platform’s cases vs. all controls', 'pooled'   # review round 12 (D-11)
    else:
        st = rtsv('controls.strata.tsv')
        st = st[st.cohort == cohort]
        groups = [(f'ctrlonly_{k}', v) for k, v in SOURCES]
        het = rtsv('controls.heterogeneity.tsv')
        het = het[het.cohort == cohort].set_index('variant_id')
        col_q, series, shapes = 'p_Q', SRC_COL, SRC_MK
        xlabel, pooled = 'OR (95% CI), all cases vs. one control source', 'all controls'
    named = FOREST_NAMED[kind]
    order = [s for s in _matrix_order() if s in mk]
    maf = {s: min(het_p.loc[mk[s]].af_pooled, 1 - het_p.loc[mk[s]].af_pooled) for s in order}
    common = [s for s in order if maf[s] >= 0.05]          # the rule applies at MAF >= 5% (scripts/checks.py)
    rare = [s for s in order if maf[s] < 0.05]
    S.setup_style('slide')
    plt.rcParams.update({'xtick.labelsize': 10.5, 'ytick.labelsize': 10.5})
    fig, ax = plt.subplots(figsize=(w, h))
    gap, ypos, y = 1.45, {}, len(order) + 1.45
    for s in common:
        ypos[s] = y
        y -= 1
    y -= gap
    for s in rare:
        ypos[s] = y
        y -= 1
    off = np.linspace(-0.27, 0.27, len(groups))
    for s in order:
        vid, y = mk[s], ypos[s]
        for k, (sub, lab) in enumerate(groups):
            s_ = st[(st.variant_id == vid) & (st.subset == sub)]
            if not len(s_):
                continue
            s_ = s_.iloc[0]
            full = s in common or (s in named and (named[s] is None or lab in named[s]))
            a = 1.0 if full else 0.6
            lo, hi = max(s_.or_l95, 0.051), min(s_.or_u95, 99)
            ax.plot([lo, hi], [y + off[k]] * 2, color=series[k], lw=1.4, alpha=a, solid_capstyle='butt')
            # review round 7: an interval cut by the axis ends in an arrowhead; an estimate beyond the axis is
            # drawn at its edge as '>' (no value label)
            if s_.or_u95 > 99:
                ax.annotate('', xy=(99, y + off[k]), xytext=(60, y + off[k]),
                            arrowprops=dict(arrowstyle='-|>', color=series[k], lw=1.4, alpha=a, mutation_scale=7))
            if s_.or_l95 < 0.051:
                ax.annotate('', xy=(0.051, y + off[k]), xytext=(0.085, y + off[k]),
                            arrowprops=dict(arrowstyle='-|>', color=series[k], lw=1.4, alpha=a, mutation_scale=7))
            if s_['or'] > 99:
                ax.plot(92, y + off[k], '>', ms=4.6, color=series[k], alpha=a, mew=0, zorder=3)
            ax.plot(s_['or'] if s_['or'] <= 99 else np.nan, y + off[k], shapes[k], ms=4.6, color=series[k],
                    alpha=a, mew=0, zorder=5, label=lab if s == order[0] else None)
        hp = het_p.loc[vid]
        a = 1.0 if (s in common or (s in named and named[s] is None)) else 0.6
        # review round 13 (F-01): a white edge parts the diamond from the estimates drawn above it; ms 5.6 keeps its outer
        # size as before (5.2 + a 1-pt edge of its own color), so the reference keeps most of its weight
        ax.plot(np.exp(hp.beta_pooled), y + 0.42, 'D', ms=5.6, color=S.INK, mec='white', mew=0.6, alpha=a, zorder=4,
                label=pooled if s == order[0] else None)
        h_ = het.loc[vid]
        c_ = S.INK if s in common else S.NEUTRAL_D
        if s in common:                           # review round 8 (D-05): Q P and I² for the markers the rule reads
            ax.text(1.02, y, f'{h_[col_q]:.2g}', transform=ax.get_yaxis_transform(), fontsize=10, va='center', color=c_)
            ax.text(1.2, y, f'{100 * h_.I2:.0f}%', transform=ax.get_yaxis_transform(), fontsize=10, va='center', color=c_)
    top = len(order) + gap
    for xx, t in ((1.02, 'Q P'), (1.2, 'I²')):
        ax.text(xx, top + 0.98, t, transform=ax.get_yaxis_transform(), fontsize=10, va='center', fontweight='bold')
    # block labels in a white box above the rows, so neither the OR = 1 line nor the frame crosses them
    box = dict(facecolor='white', edgecolor='none', pad=1.2)
    ax.text(0.064, ypos[common[0]] + 0.98, f'MAF ≥ 5%: the rule applies ({len(common)})', fontsize=10, color=NAVY,
            fontweight='bold', va='center', bbox=box, zorder=6)
    ax.text(0.064, ypos[rare[0]] + 0.98, f'MAF < 5%: outside the rule, for orientation ({len(rare)})', fontsize=10,
            color=S.INK_SOFT, fontstyle='italic', va='center', bbox=box, zorder=6)
    ax.axhline(ypos[rare[0]] + 1.0 + gap / 2, color=S.NEUTRAL_D, lw=0.7)
    ax.axvline(1, color=S.REFERENCE, lw=0.9)
    ax.set_xscale('log')
    ax.set_xlim(0.05, 100)
    ax.set_xticks([0.1, 0.3, 1, 3, 10, 30, 100], ['0.1', '0.3', '1', '3', '10', '30', '100'])
    ax.minorticks_off()
    ax.set_ylim(min(ypos.values()) - 0.6, top + 1.35)
    ax.set_yticks([ypos[s] for s in order], order)
    for t, s in zip(ax.get_yticklabels(), order):
        t.set_fontstyle('normal' if s.startswith('DQB1') else 'italic')
        if s in rare and s not in named:
            t.set_color(S.INK_SOFT)
    ax.set_xlabel(xlabel, fontsize=11)
    S.despine(ax, 'x')
    ax.grid(axis='y', visible=False)
    ax.legend(frameon=False, ncol=3, loc='lower left', bbox_to_anchor=(-0.02, 1.0), handletextpad=0.2,
              columnspacing=1.0, borderaxespad=0.2, fontsize=10)
    fig.subplots_adjust(right=0.78)
    save(fig, f'forest_{kind}.{cohort}')


# ── lead frequencies by group (ADARB2, LINC01146), full mainland ───────────
GROUP_ORDER = [
    ('case | DNBSeq-T7 30x', 'DNBSeq-T7 cases'), ('case | NovaSeq 30x', 'NovaSeq cases'),
    ('case | DNBSeq-G400RS 30x', 'DNBSeq-G400RS 30x cases'), ('case | DNBseq-G400RS 15x', 'DNBSeq-G400RS 15x cases'),
    ('case | all', 'all cases'),
    ('control | nagahama', 'Nagahama'), ('control | bbj', 'BBJ'), ('control | htlv1', 'HTLV-1'),
    ('control | acc', 'ACC'), ('control | lungca', 'Lung cancer'), ('control | all', 'all controls'),
    ('reference | ToMMo 61KJPN', 'ToMMo 61KJPN'), ('reference | 1KG JPT', '1KG JPT')]


def freq_panel(sig, w=4.55, h=4.95, cohort='full_mainland'):
    # drawn at the size of its slot (the tight box adds the count column): placed at about 1.0, every text >= 10 pt
    af = rtsv('followup.gwas_af.tsv')
    a = af[(af.signal == sig) & (af.cohort == cohort)].set_index('group')
    S.setup_style('slide')
    plt.rcParams.update({'xtick.labelsize': 10.5, 'ytick.labelsize': 10.5, 'legend.fontsize': 10})
    fig, ax = plt.subplots(figsize=(w, h))
    n = len(GROUP_ORDER)
    for i, (g, lab) in enumerate(GROUP_ORDER):
        y = n - i
        r = a.loc[g]
        col = CASE_C if g.startswith('case') else CTRL_C if g.startswith('control') else REF_C
        mk = 'D' if g.startswith('reference') else 'o'
        big = g.endswith('| all')
        ax.plot(r.alt_freq, y, mk, ms=8 if big else 6, color=col, zorder=3)
        if not pd.isna(r.n):
            ax.text(1.03, y, f'{2 * int(r.n):,}', transform=ax.get_yaxis_transform(), fontsize=10, va='center',
                    ha='left', color=S.INK_SOFT)
    ax.axvline(a.loc['control | all'].alt_freq, color=CTRL_C, lw=0.8, ls=':')
    ax.axvline(a.loc['case | all'].alt_freq, color=CASE_C, lw=0.8, ls=':')
    for yb in (n - 4.5, n - 10.5):
        ax.axhline(yb, color=S.NEUTRAL_D, lw=0.6)
    ax.set_yticks([n - i for i in range(n)], [g[1].replace('DNBSeq-G400RS', 'G400RS') for g in GROUP_ORDER])
    ax.set_ylim(0.6, n + 1.0)
    ax.set_xlabel('ALT allele frequency', fontsize=10.5)
    ax.text(1.03, n + 0.68, 'alleles', transform=ax.get_yaxis_transform(), fontsize=10, va='center', fontweight='bold')
    ax.locator_params(axis='x', nbins=5)
    ax.tick_params(axis='x', labelsize=10)
    S.despine(ax, 'x')
    ax.grid(axis='y', visible=False)
    fig.subplots_adjust(right=0.84)
    save(fig, f'freq_{sig}')


# ── carriers by group: gene sets and a driver variant, full mainland ───────
CAR_ORDER = [
    ('case | DNBSeq-T7 30x', 'DNBSeq-T7'), ('case | NovaSeq 30x', 'NovaSeq'),
    ('case | DNBSeq-G400RS 30x', 'G400RS 30x'), ('case | DNBseq-G400RS 15x', 'G400RS 15x'), ('case | all', 'all cases'),
    ('control | nagahama', 'Nagahama'), ('control | bbj', 'BBJ'), ('control | htlv1', 'HTLV-1'),
    ('control | acc', 'ACC'), ('control | lungca', 'Lung cancer'), ('control | all', 'all controls')]


def carriers(name, panels, w, h, xmax, cohort='full_mainland'):
    # panels: (gene, stratum, unit, title, mode); mode 'plus' prints the possible carriers set missing as (+n) (the
    # germline sets), 'called' the called carriers only, 'reads' one count, called or set missing (JAK2 V617F: the
    # read-based count of the fixed decision, no composite)
    gc = rtsv('followup.gene_carriers.tsv')
    gc = gc[gc.cohort == cohort]
    S.setup_style('slide')
    plt.rcParams.update({'xtick.labelsize': 10.5, 'ytick.labelsize': 10.5, 'axes.titlesize': 12})
    fig, axes = plt.subplots(1, len(panels), figsize=(w, h), sharey=True)
    n = len(CAR_ORDER) + 1
    for ax, (gene, stratum, unit, title, mode) in zip(np.atleast_1d(axes), panels):
        d = gc[(gc.gene == gene) & (gc.stratum == stratum) & (gc.unit == unit)].set_index('group')
        for i, (g, lab) in enumerate(CAR_ORDER):
            y = n - i
            r = d.loc[g]
            k = int(r.carriers + (r.hidden_as_missing if mode == 'reads' else 0))
            col = CASE_C if g.startswith('case') else CTRL_C
            ax.plot(100 * k / r.n, y, 'o', ms=7 if g.endswith('| all') else 5, color=col, zorder=3)
            t = f'{k}/{int(r.n):,}' + (f' (+{int(r.hidden_as_missing)})' if mode == 'plus' and r.hidden_as_missing else '')
            ax.annotate(t, (100 * k / r.n, y), xytext=(6, 0), textcoords='offset points', fontsize=9.5,
                        va='center', color=S.INK_SOFT)
        tom = d.loc['control | all']
        ax.plot(100 * tom.tommo_expected_carriers / tom.n, 1, 'D', ms=6, color=REF_C, zorder=3)
        for yb in (n - 4.5, n - 10.5):
            ax.axhline(yb, color=S.NEUTRAL_D, lw=0.6)
        ax.set_xlim(-0.3, xmax)
        ax.set_title(title, fontsize=12)
        ax.set_xlabel('Carriers (%)', fontsize=10.5)
        ax.grid(axis='y', visible=False)
        S.despine(ax, 'x')
    axes = np.atleast_1d(axes)
    axes[0].set_yticks([n - i for i in range(len(CAR_ORDER))] + [1], [c[1] for c in CAR_ORDER] + ['ToMMo-predicted'])
    axes[0].set_ylim(0.4, n + 0.6)
    fig.tight_layout(w_pad=0.6)
    save(fig, f'carriers_{name}')


# ── the evidence matrix at slide size ─────────────────────────────────────────────────────────────────────────────
# user 2026-11-06 ("现在的证据矩阵太让人误解了 ... 重新审视整个证据链条"): drawn from stats.json matrix2, the registry's
# seven checks shown as seven readouts by question, with the same rules and flags (collect_stats.py checks that they
# reproduce the registry's readings):
#   - leaving out a case platform is split: without the DNBSeq-T7 cases (most of the cases; threat T3) and without one
#     small platform (T2). One rule over four groups had read as if the DNBSeq-T7 cases were tested;
#   - Q and the random-effects estimate share one column (both compare the platforms' own estimates);
#   - a below-the-bar reading is a dot in its cell; the threats each readout tests stand under its name;
#   - each flag carries its follow-up fact; a bottom row gives, per readout, the signals the rule could read.
M2_NAMES = {'depth': 'Read\ndown-\nsampling', 'platform_est': 'Platform\nestimates\n(Q, RE)',
            'lopo_small': 'Without\none small\nplatform', 'lopo_t7': 'Without\nDNBSeq-T7\ncases', 'array': 'SNP\narrays',
            'sources_Q': 'Source\nestimates\n(Q)', 'sources_loo': 'Without\none\nsource'}
M2_FLAGNAME = {'depth': 'down-sampling', 'platform_est': 'platform estimates', 'lopo_small': 'without one small platform',
               'lopo_t7': 'without DNBSeq-T7', 'array': 'SNP arrays', 'sources_Q': 'source estimates',
               'sources_loo': 'without one source'}
M2_QNAME = {'Q1': 'Q1 · Depth', 'Q2': 'Q2 · Platform', 'Q3': 'Q3 · Technology', 'Q4': 'Q4 · Control source'}


def matrix(w=12.7, h=5.35):
    """The evidence matrix as one systematic table: columns = the seven readouts grouped by question (Q1-Q4), each
    with the threats it tests (T1-T5) and three cells (N / I / F; the arrays: N_arr / I_arr / F_arr; the gene sets'
    arrays: one cell); rows = the signals grouped by type. Right: each flag with its readout, cohorts and follow-up
    fact. Bottom: flag readings and the signals each readout could read; the legend with every reading's count."""
    st = json.load(open(os.path.join(HERE, 'stats.json')))
    M2 = st['matrix2']
    cols = [c['key'] for c in M2['columns']]
    qof = {c['key']: c['question'] for c in M2['columns']}
    thr = {c['key']: c['threats'] for c in M2['columns']}
    kind = {r['signal']: r['kind'] for r in st['signals']}
    rank = {'gwas': 0, 'hla_residue': 1, 'gene': 2}
    sigs = sorted([s for s in _matrix_order() if s in M2['rows']], key=lambda s: rank[kind[s]])
    if set(sigs) != set(M2['rows']) or set(cols) != set(M2_NAMES):
        sys.exit('matrix: matrix2 and the drawing disagree on the signals or the readouts')
    CN3 = {'narrow_mainland': 'N', 'intermediate_mainland': 'I', 'full_mainland': 'F'}
    # geometry, data units (the axes fill the figure)
    CELL, GAP, GAPS = 0.85, 0.1, {'array': 0.45}           # the arrays' cells wider apart: subscripted letters

    def gap(c):
        return GAPS.get(c, GAP)

    def cw(c):
        return 3 * CELL + 2 * gap(c)
    IN_G, BETWEEN = 0.55, 2.3                               # between readouts of a question, between questions
    xs, x = {}, 0.0
    for i, c in enumerate(cols):
        if i:
            x += IN_G if qof[c] == qof[cols[i - 1]] else BETWEEN
        xs[c] = x
        x += cw(c)
    FLAG_X, FW, LAB_X = x + 1.0, 12.6, -0.5

    def cells(s):
        return [v for cell in M2['rows'][s].values() for v in cell.values()]
    flagged = [s for s in sigs if any(v['r'] == 'flag' for v in cells(s))]

    def flag_line(s):
        parts = []
        for c in cols:
            cs = [CN3[k] for k in CK if M2['rows'][s][c].get(k, {}).get('r') == 'flag']
            if cs:
                parts.append(f"{M2_FLAGNAME[c][0].upper()}{M2_FLAGNAME[c][1:]}: {', '.join(cs)}")
        return '; '.join(parts)

    def follow_up(s):
        """the fact the follow-up adds to a flag (the dossier slides give the rest)"""
        if s == 'ADARB2':
            e = st['rare_platform_exact']['ADARB2|t7|full_mainland']
            if (e['alt'] != 0 or M2['rows'][s]['array']['full_mainland']['r'] != 'not checkable'
                    or st['cram']['adarb2']['full_mainland']['T7']['read'] == 0):
                sys.exit('matrix: the ADARB2 follow-up no longer holds')
            return 'T7 reads carry it; T7 calls miss it (CRAM)'
        if s == 'LINC01146':
            a = [r for r in st['array_rows'] if r['signal'] == s]
            if not a or any(r['verdict'] not in ('consistent', 'consistent on arrays') or r['conc_cases'] < 0.995 for r in a):
                sys.exit('matrix: the LINC01146 follow-up no longer holds')
            # review round 8 (A-09): the flag fires in N and I; 'OR < 1' holds in every cohort
            if any(st['linc2'][c]['strata']['g400rs15']['or'] >= 1 for c in CK):
                sys.exit('matrix: the LINC01146 follow-up no longer holds (G400RS 15x)')
            return 'from G400RS 15x (OR < 1); T7 calls confirmed'
        return ''
    # rows, top down: each type's head and its signals; a flagged signal's row is taller (two lines on the right)
    TALL = 1.2                                   # review round 8 (D-01): the flag's two lines apart
    kinds = [('gwas', 'GWAS leads'), ('hla_residue', 'HLA residue'), ('gene', 'Gene sets, M/H')]
    rows, y = [], 0.0
    for k, lab in kinds:
        ks = [s for s in sigs if kind[s] == k]
        if not ks:
            continue
        rows.append(('head', f'{lab} ({len(ks)})', y))
        y -= 0.95
        for s in ks:
            ex = TALL if s in flagged else 0.0
            y -= ex / 2
            rows.append(('sig', s, y))
            y -= 1.0 + ex / 2
    y_last = y + 1.0
    S.setup_style('slide')
    fig, ax = plt.subplots(figsize=(w, h))
    ax.axis('off')
    ax.set_xlim(LAB_X - 6.4, FLAG_X + FW)
    fig.subplots_adjust(left=0.003, right=0.997, bottom=0.003, top=0.997)

    def box(x0, yy, wd, fc, ec='none', lw=0.0, z=2):
        ax.add_patch(plt.Rectangle((x0, yy - 0.38), wd, 0.76, facecolor=fc, edgecolor=ec, lw=lw, zorder=z))

    def dot(xc, yy, on_dark):
        ax.plot([xc], [yy], marker='o', ms=3.8, color='white' if on_dark else NAVY, mec='none', zorder=5)
    R_END = FLAG_X + FW - 0.2
    for kind_, item, yy in rows:
        if kind_ == 'head':
            ax.text(LAB_X - 6.2, yy - 0.05, item, ha='left', va='center', fontsize=9.5, fontweight='bold', color=NAVY)
            ax.plot([LAB_X - 6.2, R_END], [yy - 0.48, yy - 0.48], color=S.NEUTRAL_D, lw=0.6)
            continue
        s = item
        fl = s in flagged
        if fl:                                    # a flagged signal's row on a SIG band
            ax.add_patch(plt.Rectangle((LAB_X - 6.3, yy - 0.5 - TALL / 2), R_END + 0.3 - (LAB_X - 6.3), 1.0 + TALL,
                                       facecolor=SIG, lw=0, zorder=0))
        if kind[s] == 'gene':                     # the gene in italics, the set upright
            g, _ = s.split(' (')
            gene_label(ax, LAB_X, yy, g, 'M/H', va='center', fontsize=9.5, color=INK_)
        else:                                     # plain italic text: mathtext would turn a hyphen into a minus
            ax.text(LAB_X, yy, s, ha='right', va='center', fontsize=9.5, color=NAVY if fl else INK_,
                    fontweight='bold' if fl else 'normal', style='normal' if s.startswith('DQB1') else 'italic')
        for c in cols:
            cell, x0 = M2['rows'][s][c], xs[c]
            if 'merged' in cell:                  # the gene sets' arrays: one cell for the three nested subsets
                v = cell['merged']
                if v['r'] == 'reported':          # described, no rule: an outlined cell with the count
                    box(x0, yy, cw(c), 'white', ec=S.NEUTRAL_D, lw=0.8)
                    ax.text(x0 + cw(c) / 2, yy, f"{v['count'][0]} of {v['count'][1]}", ha='center', va='center',
                            fontsize=9.5, color=NAVY, zorder=3)
                else:
                    box(x0, yy, cw(c), S.NEUTRAL)
                continue
            for m, co in enumerate(CK):
                v, cx = cell[co], x0 + m * (CELL + gap(c))
                r = v['r']
                if r == 'not applicable':         # a dash: the readout does not apply
                    ax.plot([cx + 0.35 * CELL, cx + 0.65 * CELL], [yy, yy], color=S.NEUTRAL_D, lw=0.9, zorder=2)
                    continue
                fc = {'consistent': S.DATA, 'flag': SIG}.get(r, S.NEUTRAL)
                box(cx, yy, CELL, fc, ec=NAVY if r == 'flag' else 'none', lw=1.8 if r == 'flag' else 0.0,
                    z=3 if r == 'flag' else 2)
                if r == 'flag':
                    if v['below']:
                        sys.exit(f'matrix: {s} {c} {co}: a flag cell with a below-the-bar mark')
                    ax.text(cx + CELL / 2, yy, '!', ha='center', va='center', fontsize=10, fontweight='bold', color=NAVY,
                            zorder=4)
                elif v['below']:
                    dot(cx + CELL / 2, yy, r == 'consistent')
        if fl:                                    # the flag: readout and cohorts, then the follow-up fact
            ax.text(FLAG_X, yy + 0.55, flag_line(s), ha='left', va='center', fontsize=9.5, color=NAVY, fontweight='bold')
            ax.text(FLAG_X, yy - 0.55, follow_up(s), ha='left', va='center', fontsize=9.5, color=INK_)
        # review round 8 (D-14): no mark for an unflagged signal ('–' also meant 'does not apply' in the cells)
    # headers: question, readout (three lines), the threats it tests, cohort letters (the arrays: the array subsets)
    y_top = rows[0][2] + 0.75
    for q in dict.fromkeys(qof[c] for c in cols):
        cs = [c for c in cols if qof[c] == q]
        xa, xb = xs[cs[0]], xs[cs[-1]] + cw(cs[-1])
        ax.text((xa + xb) / 2, y_top + 5.95, M2_QNAME[q], ha='center', va='center', fontsize=10, fontweight='bold',
                color=NAVY)
        ax.plot([xa, xb], [y_top + 5.4, y_top + 5.4], color=NAVY, lw=1.4)
    for c in cols:
        ax.text(xs[c] + cw(c) / 2, y_top + 3.65, M2_NAMES[c], ha='center', va='center', fontsize=9.5, color=INK_,
                linespacing=1.0)
        ax.text(xs[c] + cw(c) / 2, y_top + 1.4, ', '.join(thr[c]), ha='center', va='center', fontsize=9.5,
                color=GREY_, style='italic')
        for m, co in enumerate(CK):
            lab_ = CN3[co] if c != 'array' else rf"$\mathbf{{{CN3[co]}}}_{{\rm arr}}$"
            ax.text(xs[c] + m * (CELL + gap(c)) + CELL / 2, y_top + 0.5, lab_, ha='center', va='center',
                    fontsize=9.5, fontweight='bold', color=COH_COL[co])
    ax.text(FLAG_X, y_top + 3.65, 'Flags and follow-up', ha='left', va='center', fontsize=9.5, fontweight='bold',
            color=NAVY)
    ax.text(LAB_X, y_top + 1.4, 'Threats tested', ha='right', va='center', fontsize=9.5, color=GREY_, style='italic')
    # bottom: flag readings and the signals the rule could read, per readout
    yb = y_last - 1.15
    ax.plot([LAB_X - 6.2, R_END], [yb + 0.62, yb + 0.62], color=NAVY, lw=0.8)
    ax.text(LAB_X, yb, 'Flag readings', ha='right', va='center', fontsize=9.5, fontweight='bold', color=NAVY)
    yr = yb - 1.0
    ax.text(LAB_X, yr, 'Signals read', ha='right', va='center', fontsize=9.5, fontweight='bold', color=NAVY)
    for c in cols:
        n_ = M2['flags'][c]
        ax.text(xs[c] + cw(c) / 2, yb, str(n_), ha='center', va='center', fontsize=9.5, fontweight='bold',
                color=NAVY if n_ else S.NEUTRAL_D)
        rd3 = [M2['read'][c][k] for k in CK]
        ax.text(xs[c] + cw(c) / 2, yr, str(rd3[0]) if len(set(rd3)) == 1 else ' / '.join(map(str, rd3)), ha='center',
                va='center', fontsize=9.5, fontweight='bold', color=NAVY)
    ax.text(FLAG_X, yb, f"{sum(M2['flags'].values())} in {len(flagged)} signals", ha='left', va='center', fontsize=9.5,
            fontweight='bold', color=NAVY)
    ax.text(FLAG_X, yr, f'of {len(sigs)} (N / I / F where they differ)', ha='left', va='center', fontsize=9.5,
            color=INK_)
    # legend with the count of every reading, wrapped into rows within the table's width
    n = M2['cells']
    bbv = [v for s_ in sigs for cell in M2['rows'][s_].values() for v in cell.values() if v.get('below')]
    n_bb_ok = sum(v['r'] == 'consistent' for v in bbv)
    n_bb_nt = sum(v['r'] in ('not testable', 'not checkable') for v in bbv)
    if n_bb_ok + n_bb_nt != M2['below']:
        sys.exit('matrix: a below-the-bar dot on a cell that is neither consistent nor not testable')
    rend = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    x_lo, x_hi = LAB_X - 6.2, R_END
    yl = yr - 1.3
    items = [('box', S.DATA, None, False, f"consistent ({n.get('consistent', 0)})"),
             ('box', SIG, NAVY, False, f"flag ({n.get('flag', 0)})"),
             ('box', S.NEUTRAL, None, False, f"not testable ({n.get('not testable', 0)}) or not checkable "
                                             f"({n.get('not checkable', 0)})"),
             # a below-the-bar dot on either kind of cell (user 2026-11-06: the white dot on blue had no legend):
             # on a not-testable cell the reading itself is below the bar; on a consistent cell the other groups
             # were read and hold, and one left-out group is below the bar
             ('box', S.NEUTRAL, None, True, f"below the bar ({n_bb_nt}): nt, |z$_{{exp}}$| ≥ 1.96, |z| < 1.96"),
             ('box', S.DATA, None, True, f"consistent, one left-out group below the bar ({n_bb_ok})"),
             ('dash', None, None, False, f"does not apply ({n.get('not applicable', 0)})"),
             ('box', 'white', S.NEUTRAL_D, False, f"gene carriers confirmed on the imputed arrays, one cell "
                                                  f"({n.get('reported', 0)})")]
    xl = x_lo
    for kind_, fc, ec, below, lab in items:
        t_ = ax.text(0, yl, lab, ha='left', va='center', fontsize=9.5, color=INK_)
        bb = t_.get_window_extent(rend).transformed(inv)
        if xl > x_lo and xl + 1.2 + bb.width > x_hi:        # the next row
            xl, yl = x_lo, yl - 1.0
        t_.set_position((xl + 1.2, yl))
        if kind_ == 'box':
            ax.add_patch(plt.Rectangle((xl, yl - 0.32), 0.9, 0.64, facecolor=fc, edgecolor=ec or 'none',
                                       lw=1.4 if ec == NAVY else 0.8))
            if ec == NAVY:
                ax.text(xl + 0.45, yl, '!', ha='center', va='center', fontsize=9.5, fontweight='bold', color=NAVY)
            if below:
                dot(xl + 0.45, yl, fc == S.DATA)
        else:
            ax.plot([xl + 0.25, xl + 0.65], [yl, yl], color=S.NEUTRAL_D, lw=0.9)
        if xl + 1.2 + bb.width > x_hi:
            sys.exit(f'matrix: the legend item "{lab}" is wider than the table')
        xl += 1.2 + bb.width + 1.0
    ax.set_ylim(yl - 0.7, y_top + 6.5)
    # nothing may reach past the table's frame (else the saved figure is cropped around it, off center on the slide)
    fig.canvas.draw()
    fr = ax.get_window_extent(rend)
    out_ = [t.get_text() for t in ax.texts if not fr.contains(*t.get_window_extent(rend).p0)
            or not fr.contains(*t.get_window_extent(rend).p1)]
    if out_:
        sys.exit(f'matrix: text outside the frame: {out_}')
    save(fig, 'evidence_matrix')


# ── the follow-up readings and their rules, three cohorts ───────────────────
def _gene_label(sig):
    g, st = sig.split(' (')
    return rf'$\it{{{g}}}$ ' + ('M/H' if st.startswith('moderate') else 'L/M/H')


def _cohort_handles(color=S.INK):
    return [Line2D([], [], marker=COHORT_MK[c], ls='', color=color, mfc='white' if c == 'full_mainland' else color,
                   ms=6, label=CNAME[c].split()[0]) for c in CK]




# ── the three cohorts side by side: lead frequencies and gene carriers by group (appendix) ─────────────────────
def lead_freq_cohorts(part, w=12.7, h=5.45):
    af = rtsv('followup.gwas_af.tsv')
    leads = list(dict.fromkeys(af.signal))
    half = leads[:8] if part == 1 else leads[8:]
    S.setup_style('slide')
    plt.rcParams.update({'xtick.labelsize': 9, 'ytick.labelsize': 9.5, 'axes.titlesize': 11})
    fig, axes = plt.subplots(2, 4, figsize=(w, h), sharey=True)
    n = len(GROUP_ORDER)
    for ax, sig in zip(axes.flat, half):
        for i, (g, lab) in enumerate(GROUP_ORDER):
            y = n - i
            col = CASE_C if g.startswith('case') else CTRL_C if g.startswith('control') else REF_C
            for c in CK:
                r = af[(af.signal == sig) & (af.cohort == c) & (af.group == g)]
                if not len(r):
                    continue
                if g.startswith('reference'):
                    if c == 'full_mainland':
                        ax.plot(r.alt_freq.iloc[0], y, 'D', ms=4.5, color=REF_C, zorder=3)
                    continue
                ax.plot(r.alt_freq.iloc[0], y + COHORT_OFF[c], COHORT_MK[c], ms=3.6, color=col,
                        mfc='white' if c == 'full_mainland' else col, mew=0.8, zorder=3)
        for yb in (n - 4.5, n - 10.5):
            ax.axhline(yb, color=S.NEUTRAL_D, lw=0.6)
        ax.set_title(sig, fontstyle='italic')
        ax.grid(axis='y', visible=False)
        ax.locator_params(axis='x', nbins=4)
        S.despine(ax, 'x')
    axes[0, 0].set_yticks([n - i for i in range(n)], [g[1] for g in GROUP_ORDER])
    axes[0, 0].set_ylim(0.5, n + 0.5)
    for ax in axes[1]:
        ax.set_xlabel('ALT allele frequency', fontsize=9.5)
    fig.legend(handles=_cohort_handles(S.INK), loc='upper right', ncol=3, frameon=False, fontsize=9.5,
               bbox_to_anchor=(1.0, 1.03))
    fig.tight_layout(h_pad=0.8, w_pad=0.5)
    save(fig, f'lead_frequencies_cohorts.{part}')


def gene_carriers_cohorts(w=12.7, h=5.3):
    """Both variant sets of the four Tier 1-2 genes on one slide: the M/H sets on top, the L/M/H sets below."""
    gc = rtsv('followup.gene_carriers.tsv')
    genes = ['JAK2', 'OSBP', 'RPH3A', 'STBD1']
    sets = [(g, 'moderate_high') for g in genes] + [(g, 'low_moderate_high') for g in genes]
    S.setup_style('slide')
    plt.rcParams.update({'xtick.labelsize': 9, 'ytick.labelsize': 9, 'axes.titlesize': 10.5})
    fig, axes = plt.subplots(2, 4, figsize=(w, h), sharey=True)
    axes = axes.ravel()
    n = len(CAR_ORDER) + 1
    for ax, (gene, st) in zip(axes, sets):
        d = gc[(gc.gene == gene) & (gc.stratum == st) & (gc.unit == f'{gene} (gene)')]
        xmax = 0
        for i, (g, lab) in enumerate(CAR_ORDER):
            y = n - i
            col = CASE_C if g.startswith('case') else CTRL_C
            for c in CK:
                r = d[(d.cohort == c) & (d.group == g)]
                if not len(r):
                    continue
                r = r.iloc[0]
                xmax = max(xmax, 100 * r.carrier_share)
                ax.plot(100 * r.carrier_share, y + COHORT_OFF[c], COHORT_MK[c], ms=4, color=col,
                        mfc='white' if c == 'full_mainland' else col, mew=0.8, zorder=3)
                if r.hidden_as_missing and gene != 'JAK2':
                    xmax = max(xmax, 100 * (r.carriers + r.hidden_as_missing) / r.n)
                    ax.plot(100 * (r.carriers + r.hidden_as_missing) / r.n, y + COHORT_OFF[c], '|', ms=6, color=col,
                            mew=1.2, zorder=2)
        for c in CK:
            t = d[(d.cohort == c) & (d.group == 'control | all')]
            if len(t):
                ax.plot(100 * t.tommo_expected_carriers.iloc[0] / t.n.iloc[0], 1 + COHORT_OFF[c], 'D', ms=3.6,
                        color=REF_C, mfc='white' if c == 'full_mainland' else REF_C, zorder=3)
        for yb in (n - 4.5, n - 10.5):
            ax.axhline(yb, color=S.NEUTRAL_D, lw=0.6)
        ax.set_xlim(-0.3, max(4, xmax * 1.15))
        ax.set_title(f'$\\it{{{gene}}}$ ' + ('M/H' if st == 'moderate_high' else 'L/M/H'), pad=3)
        if st == 'low_moderate_high':
            ax.set_xlabel('Carriers (%)', fontsize=9.5)
        ax.grid(axis='y', visible=False)
        S.despine(ax, 'x')
    for ax in axes[[0, 4]]:
        ax.set_yticks([n - i for i in range(len(CAR_ORDER))] + [1], [c[1] for c in CAR_ORDER] + ['ToMMo-predicted'])
        ax.set_ylim(0.5, n + 0.5)
    hs = _cohort_handles(S.INK) + [Line2D([], [], marker='|', ls='', color=S.INK, ms=7, mew=1.2,
                                          label='with possible carriers set missing')]
    fig.legend(handles=hs, loc='upper right', ncol=4, frameon=False, fontsize=9.5, bbox_to_anchor=(1.0, 1.035))
    fig.tight_layout(h_pad=0.9, w_pad=0.6)
    save(fig, 'gene_carriers_cohorts')


# ── leave one group out: z observed against z expected, every cohort ──────────
def _short_set(sig):
    if ' (' in sig:
        g, st = sig[:-1].split(' (')
        return f'{g} ' + ('M/H' if st.startswith('moderate') else 'L/M/H')
    return sig


def loo_scatter(chk, w=7.15, h=3.2):
    # |z_exp| below void_z: the left-out group holds almost all of a reading's information (stats.json const, so the
    # slide's caption names the same readings)
    VOID_Z = json.load(open(os.path.join(HERE, 'stats.json')))['const']['void_z']
    """one check (lopo: platform left out; sources_loo: source left out), one panel per cohort"""
    lo = rtsv('evidence.loo.tsv').dropna(subset=['z_observed', 'z_expected'])
    lo['x'] = lo.z_expected.abs()
    lo['y'] = lo.z_observed * np.sign(lo.z_expected)          # signed to the pooled effect
    zx, zo = DF.CK.LOO_Z_EXPECTED, DF.CK.LOO_Z_OBSERVED
    S.setup_style('slide')
    fig, axes = plt.subplots(1, 3, figsize=(w, h), sharex=True, sharey=True)
    lim = (0, max(7.5, float(lo.x.max()) + 0.3))
    panels = []                                   # per panel: markers, labels and boxes, for the labels below the bar
    nflag = 0                                     # review round 8 (D-04): "flag" in the legend only if a flag is drawn
    if True:
        for col, c in enumerate(CK):
            ax = axes[col]
            x = lo[(lo.check == chk) & (lo.cohort == c)]
            ax.add_patch(plt.Rectangle((-0.3, -3), zx + 0.3, lim[1] + 3, facecolor='#F5F5F5', edgecolor='none', zorder=0))
            ax.add_patch(plt.Rectangle((zx, -3), lim[1], 3 + zo, color=SIG, alpha=0.9, lw=0, zorder=0))
            # below the bar: an unchanged effect would stay significant (|z_exp| >= 1.96), the rule does not read it
            ax.add_patch(plt.Rectangle((zo, -3), zx - zo, 3 + zo, facecolor='#DDE1E6', edgecolor='none', zorder=0))
            cc = COH_COL[c]
            ax.plot(lim, lim, color=S.REFERENCE, lw=0.7, zorder=1)
            # review round 13 (F-02): the bar lines as on loo_regions; gene sets filled (every point is a check), so the
            # only open marks are the readings below the bar
            ax.axhline(zo, color=GREY_, lw=0.8, ls=':')
            ax.axvline(zx, color=GREY_, lw=0.8, ls=':')
            flags, nts, marks_ = [], [], []
            for kind, mk in [('gwas', 'o'), ('hla_residue', 's'), ('gene', '^')]:
                y = x[x.kind == kind]
                bad = (y.x >= zx) & (y.y < zo)                    # the evidence-matrix rule (scripts/checks.py)
                bb = (y.x < zx) & (y.x >= zo) & (y.y < zo)          # not testable, yet |z| < 1.96
                ok = ~bad & ~bb
                ax.scatter(y.x[ok], y.y[ok], s=14, marker=mk, facecolor=cc,
                           edgecolor=cc, lw=0.7, zorder=3)
                ax.scatter(y.x[bad], y.y[bad], s=34, marker=mk, facecolor=NAVY, edgecolor=NAVY, lw=1.3, zorder=4)  # round 12, D-07
                ax.scatter(y.x[bb], y.y[bb], s=30, marker=mk, facecolor='white', edgecolor=GREY_, lw=1.1, zorder=4)
                marks_ += ([(a, b, 14) for a, b in zip(y.x[ok], y.y[ok])] + [(a, b, 34) for a, b in zip(y.x[bad], y.y[bad])]
                           + [(a, b, 30) for a, b in zip(y.x[bb], y.y[bb])])
                flags += [(p.signal, p.x, p.y) for p in y[bad].itertuples()]
                nts += [(p.signal, p.x, p.y) for p in y[bb].itertuples()]
            jk = x[x.signal.str.startswith('JAK2')]
            obst, jbox = [], None
            if len(jk):                                      # JAK2's burden is weak: every reading not testable
                from matplotlib.patches import FancyBboxPatch
                x0b, x1b = max(0.05, float(jk.x.min()) - 0.25), float(jk.x.max()) + 0.25
                y0b, y1b = float(jk.y.min()) - 0.35, float(jk.y.max()) + 0.35
                ax.add_patch(FancyBboxPatch((x0b, y0b), x1b - x0b, y1b - y0b, boxstyle='round,pad=0,rounding_size=0.25',
                                            fill=False, ec=GREY_, lw=0.8, ls='--', zorder=2))
                # review round 13: the bar line at 1.96, darker since F-02, must not run through the label; lift the
                # label clear of it when it would sit on the line (the source panels)
                ly = y1b + 0.15 if not zo - 0.6 < y1b + 0.15 < zo + 0.2 else zo + 0.2
                obst.append(ax.text(x0b, ly, 'JAK2\n(nt)', fontsize=9, color=GREY_, style='italic',
                                    ha='left', va='bottom', linespacing=1.0, zorder=5))
                jbox = (x0b, y0b, x1b, y1b)
            # a reading with almost no information left (|z_exp| < 0.5: the left-out group holds nearly all of it, e.g.
            # OSBP M/H without Nagahama, whose control carriers are all there) is named under its point (user
            # 2026-11-06: an unexplained point at 0); the labels below the bar keep clear of it
            void = x[(x.x < VOID_Z) & ~x.signal.str.startswith('JAK2')]
            if len(void) > 1:
                sys.exit(f'loo_scatter: more than one reading without information in {c}: {list(void.signal)}')
            for p in void.itertuples():
                obst.append(ax.text(p.x + 0.1, p.y - 0.3, _short_set(p.signal), fontsize=9, color=GREY_,
                                    style='italic', ha='left', va='top', zorder=5))
            for k, (sig, px, py) in enumerate(sorted(flags, key=lambda v: -v[2])):
                obst.append(ax.annotate(sig, (px, py), xytext=(px + 0.35, 1.2 - 1.15 * k), fontsize=9, color=NAVY,
                            style='italic', fontweight='bold', va='center',
                                        arrowprops=dict(arrowstyle='-', color=NAVY, lw=0.5)))
            nflag += len(flags)
            panels.append({'ax': ax, 'marks': marks_, 'nts': nts, 'obst': obst, 'box': jbox})
            ax.set_xlim(-0.3, lim[1])
            ax.set_ylim(-3, lim[1])
            ax.tick_params(labelsize=9)
            ptitle(ax, 'abc'[col], CNAME[c], size=10)
            ax.title.set_color(cc)
            ax.set_xlabel('|z| expected', fontsize=9.5)
        axes[0].set_ylabel('z observed, signed\nto the pooled effect', fontsize=9.5)
    hs = [Line2D([], [], marker='o', ls='', color=GREY_, label='GWAS lead'),
          Line2D([], [], marker='s', ls='', color=GREY_, label='DQB1:57:V'),
          Line2D([], [], marker='^', ls='', color=GREY_, label='gene set'),
          Line2D([], [], marker='o', ls='', mfc=NAVY, mec=NAVY, mew=1.3, label='flag'),
          Line2D([], [], marker='o', ls='', mfc='white', mec=GREY_, mew=1.1, label='below the bar')]
    hs = [h_ for h_ in hs if h_.get_label() != 'flag' or nflag]
    fig.legend(handles=hs, fontsize=9, frameon=False, loc='lower center', ncol=len(hs), bbox_to_anchor=(0.5, 0.995),
               handletextpad=0.2, columnspacing=0.9)
    fig.tight_layout(h_pad=0.6, w_pad=0.4)
    _label_below_bar(fig, panels)
    save(fig, f"loo_scatter.{'platform' if chk == 'lopo' else 'sources'}")


def _label_below_bar(fig, panels):
    """the labels of the readings below the bar (review round 7: the leaders crossed labels and other leaders). The
    square marker names DQB1:57:V, so only gene sets and GWAS leads get a label: one row each in the empty lower left,
    the leader from the label's right end or its top edge; the first arrangement (row order, alignment, leader start)
    that keeps clear of every marker, label, leader and the JAK2 box is drawn; else a label left or right of the point"""
    import itertools
    from matplotlib.transforms import Bbox
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    ppt = fig.dpi / 72

    def seg_dist(p, a, b):
        ab, ap = b - a, p - a
        t = np.clip(np.dot(ap, ab) / max(np.dot(ab, ab), 1e-9), 0, 1)
        return float(np.hypot(*(a + t * ab - p)))

    def ccw(a, b, c):
        return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])

    def seg_cross(a, b, c, d):
        return ccw(a, c, d) != ccw(b, c, d) and ccw(a, b, c) != ccw(a, b, d)

    def seg_box(a, b, bb):
        if bb.contains(*a) or bb.contains(*b):
            return True
        cs = [np.array(c) for c in ((bb.x0, bb.y0), (bb.x1, bb.y0), (bb.x1, bb.y1), (bb.x0, bb.y1))]
        return any(seg_cross(a, b, cs[i], cs[(i + 1) % 4]) for i in range(4))

    def disc_box(q, r, bb):
        return bb.x0 - r < q[0] < bb.x1 + r and bb.y0 - r < q[1] < bb.y1 + r
    rows = [-0.9, -1.6, -2.3]
    for P in panels:
        ax, T = P['ax'], P['ax'].transData
        inv = T.inverted()
        base = [a_.get_window_extent(rend).expanded(1.0, 1.0) for a_ in P['obst']]
        if P['box'] is not None:
            x0, y0, x1, y1 = P['box']
            base.append(Bbox([T.transform((x0, y0)), T.transform((x1, y1))]))
        marks = [(T.transform((u, v)), np.sqrt(sz) / 2 * ppt) for u, v, sz in P['marks']]
        frame = ax.get_window_extent(rend)
        items = [n for n in P['nts'] if not n[0].startswith('DQB1')]
        r_t = np.sqrt(30) / 2 * ppt

        def fits(bb, seg, tgt, boxes, leaders):
            others = [(q, r) for q, r in marks if np.hypot(*(q - tgt)) > 1.0]
            if not (frame.contains(bb.x0, bb.y0) and frame.contains(bb.x1, bb.y1)):
                return False
            if any(bb.overlaps(b2) for b2 in boxes) or any(disc_box(q, r + 1.5, bb) for q, r in others + [(tgt, r_t)]):
                return False
            if seg is None:
                return True
            return (all(seg_dist(q, *seg) > r + 2.0 for q, r in others) and not any(seg_box(*seg, b2) for b2 in boxes)
                    and not any(seg_cross(*seg, *l2) for l2 in leaders))

        def place(sig, px, py, k, aligns, boxes, leaders):
            """one label in row k, else in a later row, then an earlier one (user 2026-11-06: a label named under its
            point, OSBP M/H at 0, can hold row k), else beside the point: the first candidate that fits, as (text,
            bbox, leader)"""
            tgt = T.transform((px, py))
            kw = dict(fontsize=9, color=GREY_, va='center', zorder=5, style='italic')
            for kr in list(range(k, len(rows))) + list(range(k)):
                for al in aligns:
                    t = ax.text(px - 0.15 if al == 'right' else -0.25, rows[kr], _short_set(sig), ha=al, **kw)
                    bb = t.get_window_extent(rend)
                    if inv.transform((bb.x0, bb.y0))[0] < -0.28:          # too long to end left of its point
                        t.remove()
                        continue
                    for st in ('rm', 0.75, 0.5, 0.25):
                        a_ = (np.array([bb.x1 + 2.0, (bb.y0 + bb.y1) / 2]) if st == 'rm'
                              else np.array([bb.x0 + st * (bb.x1 - bb.x0), bb.y1 + 1.0]))
                        u_ = (tgt - a_) / np.hypot(*(tgt - a_))
                        seg = (a_, tgt - u_ * (r_t + 1.0))
                        if fits(bb, seg, tgt, boxes, leaders):
                            return t, bb, seg
                    t.remove()
            for side, dx, dy, ha in (('left', -0.2, 0.0, 'right'), ('right', 0.2, 0.18, 'left')):
                t = ax.text(px + dx, py + dy, _short_set(sig), ha=ha, **kw)
                bb = t.get_window_extent(rend)
                if fits(bb, None, tgt, boxes, leaders):
                    return t, bb, None
                t.remove()
            return None
        placed = None
        for pref in (('right', 'left'), ('left', 'right')):
            for order in itertools.permutations(items):
                boxes, leaders, arts = list(base), [], []
                for k, (sig, px, py) in enumerate(order):
                    got = place(sig, px, py, k, pref, boxes, leaders)
                    if got is None:
                        break
                    t, bb, seg = got
                    arts.append(t)
                    boxes.append(bb.expanded(1.0, 1.0))
                    if seg is not None:
                        leaders.append(seg)
                else:
                    placed = (arts, leaders)
                    break
                for t in arts:
                    t.remove()
            if placed:
                break
        if placed is None:
            sys.exit(f'loo_scatter: no arrangement for the labels {[n[0] for n in items]}')
        for seg in placed[1]:
            (xa, ya), (xb, yb) = inv.transform(seg[0]), inv.transform(seg[1])
            ax.plot([xa, xb], [ya, yb], color=GREY_, lw=0.7, zorder=3)


def loo_regions(w=5.6, h=4.5):
    st = json.load(open(os.path.join(HERE, 'stats.json')))
    zx, zo = DF.CK.LOO_Z_EXPECTED, DF.CK.LOO_Z_OBSERVED
    ex = next(r for r in st['loo_flag_basis'] if r['signal'] == 'ADARB2' and r['cohort'] == 'narrow_mainland')
    ex_x, ex_y = abs(ex['z_exp']), ex['z'] * (1 if ex['z_exp'] >= 0 else -1)
    lo = rtsv('evidence.loo.tsv')
    sd = float(lo[(lo.signal == 'ADARB2') & (lo.check == 'lopo') & (lo.cohort == 'narrow_mainland') &
                  (lo.left_out == 't7')].sd_conditional.iloc[0])
    S.setup_style('slide')
    fig, ax = plt.subplots(figsize=(w, h))
    X1, Y0, Y1 = 6.0, -1.0, 6.0
    ax.add_patch(plt.Rectangle((zx, Y0), X1 - zx, zo - Y0, color=SIG, lw=0, zorder=0))
    ax.add_patch(plt.Rectangle((0, Y0), zx, Y1 - Y0, facecolor='#F5F5F5', edgecolor='none', zorder=0))
    ax.add_patch(plt.Rectangle((zo, Y0), zx - zo, zo - Y0, facecolor='#DDE1E6', edgecolor='none', zorder=0))
    ax.plot([0, X1], [0, X1], color=S.REFERENCE, lw=1.0, zorder=1)
    ax.axhline(zo, color=GREY_, lw=0.9, ls=':', zorder=1)
    for x in (zo, zx):
        ax.axvline(x, color=GREY_, lw=0.9, ls=':', zorder=1)
    kw = dict(fontsize=10, va='center', ha='center', zorder=5)
    ax.text(5.5, 1.55, 'flag', color=NAVY, fontweight='bold', **kw)
    ax.text((zo + zx) / 2, -0.55, 'below\nthe bar', color=INK_, linespacing=1.0, **{**kw, 'fontsize': 9})
    ax.annotate('', xy=(0.05, 5.35), xytext=(zx - 0.05, 5.35),
                arrowprops=dict(arrowstyle='<->', color=GREY_, lw=0.8), zorder=5)
    ax.text(zx / 2, 5.55, r'not testable (nt): |z$_{\rm exp}$| < ' + f'{zx}', color=GREY_, fontsize=9, ha='center',
            va='bottom', zorder=5)
    ax.text(zo / 2, 3.2, 'even an\nunchanged effect\nis not significant', color=GREY_, fontsize=9,
            ha='center', va='center', linespacing=1.0, zorder=5)
    ax.text(5.0, 2.35, 'consistent', color=INK_, **kw)
    ax.text(5.55, 4.95, 'unchanged\neffect', color=GREY_, fontsize=9, ha='center', va='top', linespacing=1.0)
    # the spread of z if the effect were unchanged, drawn sideways at the example's z_exp; its tail below 1.96 shaded
    yy = np.linspace(ex_x - 3.2 * sd, ex_x + 3.2 * sd, 200)
    dens = norm.pdf(yy, ex_x, sd)
    k = 0.75 / dens.max()
    ax.plot([ex_x, ex_x], [yy[0], yy[-1]], color=NAVY, lw=0.6, zorder=3)          # the baseline of the sideways bell
    ax.plot(ex_x + k * dens, yy, color=NAVY, lw=1.1, zorder=3)
    ax.text(ex_x + 0.85, ex_x, 'z if the effect\nwere unchanged', fontsize=9, color=NAVY, ha='left',   # round 12, D-09
            va='center', linespacing=1.0, zorder=5)
    tail = yy <= ex_y
    ax.fill_betweenx(yy[tail], ex_x, ex_x + k * dens[tail], color=NAVY, alpha=0.8, lw=0, zorder=3)
    ax.plot(ex_x, ex_x, 'o', ms=6, mfc='white', mec=NAVY, mew=1.3, zorder=4)
    ax.plot(ex_x, ex_y, 'o', ms=8, mfc=NAVY, mec=NAVY, mew=1.5, zorder=4)
    ax.annotate(f'z {ex_y:.2f}, z$_{{\\rm exp}}$ {ex_x:.2f}; '
                f'chance {ex["p_unchanged"]:.3f}', (ex_x, ex_y), xytext=(3.45, -0.45), fontsize=9,
                color=NAVY, va='center', ha='left', arrowprops=dict(arrowstyle='-', color=NAVY, lw=0.6), zorder=6)
    ax.set_xlim(0, X1)
    ax.set_ylim(Y0, Y1)
    ax.set_xticks([0, zo, zx, 4, 6], ['0', f'{zo}', f'{zx}', '4', '6'])
    ax.set_yticks([0, zo, 4, 6], ['0', f'{zo}', '4', '6'])
    ax.tick_params(labelsize=10)
    ax.set_xlabel(r'|z| expected without the platform (z$_{\rm exp}$)', fontsize=10.5)
    ax.set_ylabel('z observed without the platform', fontsize=10.5)
    ax.grid(False)
    S.despine(ax, 'xy')
    fig.tight_layout()
    save(fig, 'loo_regions')


# ── SNP arrays: the three kinds of evidence of one cohort ──────────────────────
def array_z(w=7.7, h=4.75, fs=10):
    """Q3, GWAS leads and DQB1:57:V (user 2026-11-05: the rule is the z difference only, flag if the opposite
    direction or |z WGS - z array| > 3; visualize the two z): a, full mainland, z on WGS (open) and on the array data
    (filled), one marker per signal, signed so that the WGS direction is positive (the opposite direction shaded: flag);
    b, z WGS - z array in the three cohorts, same sign, the flag zone |dz| > 3 shaded"""
    ev = rtsv('array.evidence.tsv')
    ev = ev[ev.kind != 'gene']
    sigs = [s_ for s_ in _matrix_order() if s_ in set(ev.signal)]
    marks = {'typed': 'o', 'imputed': 'D', 'LD proxy': 'P', 'imputed (HIBAG)': 's'}
    S.setup_style('slide')
    fig, axes = plt.subplots(1, 2, figsize=(w, h), sharey=True, gridspec_kw={'width_ratios': [1.15, 1.0]})
    ax = axes[0]
    cc = COH_COL['full_mainland']
    ax.axvspan(-2.2, 0, color=SIG, lw=0, zorder=0)
    ax.text(-1.1, len(sigs) - 0.45, 'opposite:\nflag', ha='center', va='top', fontsize=fs - 1, color=NAVY)
    f = ev[ev.cohort == 'full_mainland']
    for i, s_ in enumerate(sigs):
        y = len(sigs) - 1 - i
        r = f[f.signal == s_].iloc[0]
        if r.evidence not in marks or not np.isfinite(r.z_array):
            ax.text(0.15, y, 'not on the arrays', va='center', fontsize=fs - 1, color=GREY_, style='italic')
            continue
        sg = np.sign(r.z_wgs)
        zw, za = abs(r.z_wgs), sg * r.z_array
        ax.plot([zw, za], [y, y], color=cc, lw=1.0, alpha=0.75, zorder=2)
        ax.scatter(zw, y, s=38, marker=marks[r.evidence], facecolor='white', edgecolor=cc, lw=1.2, zorder=3)
        ax.scatter(za, y, s=34, marker=marks[r.evidence], color=cc, lw=0, zorder=4)
    ax.axvline(0, color=S.REFERENCE, lw=0.8)
    ax.set_xlim(-2.2, 6.6)
    ax.set_yticks(range(len(sigs))[::-1], sigs, fontsize=fs - 0.5)
    for tl, s_ in zip(ax.get_yticklabels(), sigs):
        tl.set_fontstyle('normal' if s_.startswith('DQB1') else 'italic')
    ax.set_ylim(-0.7, len(sigs) - 0.3)
    ax.tick_params(axis='x', labelsize=fs)
    ax.tick_params(axis='y', length=0)
    ax.grid(axis='y', visible=False)
    ax.set_xlabel('z: WGS (open), array data (filled)', fontsize=fs)
    ptitle(ax, 'a', r'F$_{\rm arr}$', size=fs + 0.5)
    S.despine(ax, 'x')
    ax = axes[1]
    lim = CONST_ARRAY_DZ = 3.0
    C3_ = {'narrow_mainland': r'N$_{\rm arr}$', 'intermediate_mainland': r'I$_{\rm arr}$', 'full_mainland': r'F$_{\rm arr}$'}
    ax.axvspan(-4.6, -lim, color=SIG, lw=0, zorder=0)
    ax.axvspan(lim, 4.6, color=SIG, lw=0, zorder=0)
    for xv in (-lim, lim):
        ax.axvline(xv, color=NAVY, lw=0.9, ls='--', zorder=1)
    ax.text(3.8, len(sigs) - 0.45, 'flag', ha='center', va='top', fontsize=fs - 1, color=NAVY)
    ax.text(-3.8, len(sigs) - 0.45, 'flag', ha='center', va='top', fontsize=fs - 1, color=NAVY)
    ax.axvline(0, color=S.REFERENCE, lw=0.8)
    for c in CK:
        x = ev[ev.cohort == c]
        for i, s_ in enumerate(sigs):
            r = x[x.signal == s_]
            if not len(r) or not np.isfinite(r.iloc[0].z_array):
                continue
            r = r.iloc[0]
            d = np.sign(r.z_wgs) * (r.z_wgs - r.z_array)
            y = len(sigs) - 1 - i + COHORT_OFF[c]
            ax.scatter(d, y, s=22, marker='o', color=COH_COL[c],
                       facecolor=COH_COL[c], lw=1.0, zorder=3)   # review round 12 (D-08): every marker is a check, filled
    ax.set_xlim(-4.6, 4.6)
    # review round 13 (F-07): plain tick values (the FixedLocator printed '−3.0', '0.0'); upright subscripts, as on the slides
    ax.set_xticks([-3, -1.5, 0, 1.5, 3], ['−3', '−1.5', '0', '1.5', '3'])
    ax.tick_params(axis='x', labelsize=fs)
    ax.tick_params(axis='y', length=0)
    ax.grid(axis='y', visible=False)
    ax.set_xlabel(r'z$_{\rm WGS}$ − z$_{\rm array}$ (WGS direction > 0)', fontsize=fs)
    ax.legend(handles=[Line2D([], [], marker='o', ls='', color=COH_COL[c],
                              mfc=COH_COL[c], ms=5.5, label=C3_[c])
                       for c in CK], loc='lower left', bbox_to_anchor=(0.17, 0.0), ncol=1, fontsize=fs - 1,
              frameon=False, handletextpad=0.15, columnspacing=0.6, borderaxespad=0.1)
    ptitle(ax, 'b', 'Difference', size=fs + 0.5)
    S.despine(ax, 'x')
    fig.tight_layout(w_pad=0.6)
    save(fig, 'array_z')


def array_genes(w=5.5, h=3.4, fs=10, cohort='full_mainland'):
    """Q3, gene drivers (user 2026-11-05: full mainland, the imputed arrays, a description: no verdict, no expected
    number): the WGS case carriers of each driver and those the imputed arrays confirm (dosage >= 0.5)"""
    ev = rtsv('array.evidence.tsv')
    c = ev[(ev.cohort == cohort) & (ev.kind == 'gene') & (ev.verdict == 'counted')].copy()
    S.setup_style('slide')
    fig, ax = plt.subplots(figsize=(w, h))
    d = c.iloc[::-1]
    yy = np.arange(len(d))
    ax.barh(yy, d.carriers_wgs, color='#EAD3D4', edgecolor=CASE_C, lw=0.8, height=0.62, label='WGS case carriers')  # legend order as drawn
    ax.barh(yy, d.carriers_found, color=CASE_C, height=0.62, label='confirmed on the imputed arrays')
    for y, (f, n) in enumerate(zip(d.carriers_found, d.carriers_wgs)):
        ax.text(n + 0.25, y, f'{int(f)} of {int(n)}', va='center', fontsize=fs - 0.5, color=INK_)
    ax.set_yticks(yy, ['' for _ in d.marker], fontsize=fs - 0.5)
    gene_ticklabels(ax, [(y, m.split(' ')[0], m.split(' ', 1)[1].replace('p.', '')) for y, m in zip(yy, d.marker)],
                    fs - 0.5, pad_pt=3.5)
    ax.set_xlim(0, max(d.carriers_wgs) + 3.2)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=4, integer=True))
    ax.tick_params(axis='x', labelsize=fs)
    ax.tick_params(axis='y', length=0)
    ax.grid(axis='y', visible=False)
    ax.set_xlabel('Case carriers', fontsize=fs)
    ax.legend(loc='lower right', fontsize=fs - 1, frameon=False, ncol=1, handlelength=1.2, labelspacing=0.3)
    ax.set_title(r'Gene drivers, F$_{\rm arr}$ cases', fontsize=fs + 0.5, loc='center')
    S.despine(ax, 'x')
    fig.tight_layout()
    save(fig, 'array_genes')


def hibag_algorithm(w=12.7, h=4.25):
    """The HIBAG typing algorithm as a worked example, drawn like the HLA-HD slide of hla.typing.v2 (user
    2026-11-05): 1 learn from reference haplotypes (each classifier a bootstrap sample and a random set of SNPs),
    2 every classifier scores the allele pairs of one person from that person's SNPs, 3 the average gives the best pair
    and its posterior; residue 57 read from the two alleles. Illustrative values; the model's real numbers are on the
    slide."""
    from matplotlib.patches import FancyBboxPatch
    S.setup_style('slide')
    fig, ax = plt.subplots(figsize=(w, h))
    ax.set_xlim(0, w)
    ax.set_ylim(0, h)
    ax.axis('off')
    BLUE, LBLUE, TAN, TANE = '#3A6EA5', '#8FB4D4', '#FBEEE1', '#C8A15A'

    def head(x, n, t):
        ax.add_patch(plt.Circle((x + 0.17, h - 0.22), 0.16, color=NAVY, zorder=3))
        ax.text(x + 0.17, h - 0.22, str(n), ha='center', va='center', color='white', fontsize=11, fontweight='bold', zorder=4)
        ax.text(x + 0.42, h - 0.22, t, ha='left', va='center', fontsize=11.5, fontweight='bold', color=INK_)

    def dot(x, y, on, r=0.085):
        ax.add_patch(plt.Circle((x, y), r, facecolor=BLUE if on else 'white', edgecolor=BLUE if on else '#BFBFBF',
                                lw=1.0, zorder=3))

    # 1 · learn from the reference people
    head(0.0, 1, 'Learn from typed reference people')
    snps = 'abcdef'
    haps = [('●○●●○●', '05:01'), ('●○●○○●', '05:01'), ('○●○●●○', '06:01'), ('○●○●●●', '06:01'),
            ('○●●●○○', '06:04'), ('●●○○●○', '04:01')]
    xs = [1.3 + 0.36 * i for i in range(6)]
    y0, dy = 3.35, 0.36
    ax.text(xs[0] - 0.3, y0 + 0.32, 'SNP', ha='right', va='bottom', fontsize=9.2, fontweight='bold', color=INK_)
    for i, c in enumerate(snps):
        ax.text(xs[i], y0 + 0.32, c, ha='center', va='bottom', fontsize=9.4, fontweight='bold', color=INK_)
    sub = {0, 2, 5}                                       # the SNPs one classifier draws
    for i in sub:
        ax.add_patch(plt.Rectangle((xs[i] - 0.16, y0 - dy * 5 - 0.2), 0.32, dy * 5 + 0.42, facecolor='#E1EBF5',
                                   edgecolor='none', zorder=1))
    ax.text(xs[-1] + 0.6, y0 + 0.32, 'DQB1', ha='center', va='bottom', fontsize=9.4, fontweight='bold', color=INK_)
    for r, (pat, al) in enumerate(haps):
        y = y0 - dy * r
        ax.plot([0.12, 0.55], [y, y], color=BLUE, lw=2.6, solid_capstyle='round')
        ax.text(0.62, y, f'hap {r + 1}', ha='left', va='center', fontsize=9.2, color=GREY_)
        for i, ch in enumerate(pat):
            dot(xs[i], y, ch == '●')
        ax.text(xs[-1] + 0.6, y, al, ha='center', va='center', fontsize=9.6, color=INK_, fontweight='bold')
    ax.text(0.12, y0 - dy * 6 + 0.05, 'filled: the minor allele; shaded: the SNPs one classifier\n'
            'draws, from a bootstrap sample of the people; every\nclassifier learns which SNP haplotypes go with which '
            'allele', ha='left', va='top', fontsize=9.2, color=GREY_, linespacing=1.3)
    ax.annotate('', xy=(4.55, 2.45), xytext=(4.22, 2.45),
                arrowprops=dict(arrowstyle='-|>', color=INK_, lw=1.6, mutation_scale=14))

    # 2 · every classifier scores the pairs of one person
    x2 = 4.7
    head(x2, 2, 'Score every allele pair for one person')
    ax.text(x2 + 0.05, 3.45, 'one person, SNP a–f:', ha='left', va='center', fontsize=9.4, color=INK_)
    for i, g in enumerate([1, 1, 1, 2, 1, 1]):
        ax.text(x2 + 1.85 + 0.27 * i, 3.45, str(g), ha='center', va='center', fontsize=10, fontweight='bold', color=BLUE)
    ax.text(x2 + 0.05, 3.18, '(copies of the minor allele, from the array)', ha='left', va='center', fontsize=9.2,
            color=GREY_)
    cls = [('classifier 1', [('05:01 + 06:01', 0.80), ('05:01 + 06:04', 0.15), ('04:01 + 06:01', 0.05)]),
           ('classifier 2', [('05:01 + 06:01', 0.70), ('05:01 + 04:01', 0.25), ('05:01 + 06:04', 0.05)])]
    yb = 2.75
    for name, pairs in cls:
        ax.text(x2 + 0.05, yb, name, ha='left', va='center', fontsize=9.4, fontweight='bold', color=INK_)
        for k, (pr, v) in enumerate(pairs):
            y = yb - 0.27 * (k + 1)
            ax.text(x2 + 0.2, y, pr, ha='left', va='center', fontsize=9.2, color=INK_)
            ax.add_patch(plt.Rectangle((x2 + 1.55, y - 0.09), 1.6 * v, 0.18, facecolor=BLUE if k == 0 else LBLUE,
                                       edgecolor='none'))
            ax.text(x2 + 1.6 + 1.6 * v, y, f'{v:.2f}', ha='left', va='center', fontsize=9.2, color=INK_)
        yb -= 1.15
    ax.text(x2 + 0.05, 0.42, '… and so for every classifier', ha='left', va='center', fontsize=9.2, color=GREY_,
            style='italic')
    ax.annotate('', xy=(8.6, 2.45), xytext=(8.27, 2.45),
                arrowprops=dict(arrowstyle='-|>', color=INK_, lw=1.6, mutation_scale=14))

    # 3 · average, report the best pair, read residue 57
    x3 = 8.75
    head(x3, 3, 'Average; report the best pair')
    avg = [('05:01 + 06:01', 0.75), ('05:01 + 04:01', 0.13), ('05:01 + 06:04', 0.10)]
    ax.text(x3 + 0.05, 3.45, 'average over the classifiers', ha='left', va='center', fontsize=9.4, fontweight='bold',
            color=INK_)
    for k, (pr, v) in enumerate(avg):
        y = 3.18 - 0.27 * k
        ax.text(x3 + 0.2, y, pr, ha='left', va='center', fontsize=9.2, color=INK_, fontweight='bold' if k == 0 else 'normal')
        ax.add_patch(plt.Rectangle((x3 + 1.55, y - 0.09), 1.6 * v, 0.18, facecolor=BLUE if k == 0 else LBLUE,
                                   edgecolor='none'))
        ax.text(x3 + 1.6 + 1.6 * v, y, f'{v:.2f}', ha='left', va='center', fontsize=9.2, color=INK_)
    ax.add_patch(FancyBboxPatch((x3 + 0.05, 1.62), 3.85, 0.62, boxstyle='round,pad=0.02,rounding_size=0.06',
                                facecolor='#E3ECF6', edgecolor=BLUE, lw=1.4))
    ax.text(x3 + 1.97, 2.0, 'DQB1*05:01 / 06:01', ha='center', va='center', fontsize=13, fontweight='bold', color=NAVY)
    ax.text(x3 + 1.97, 1.75, 'posterior 0.75', ha='center', va='center', fontsize=9.4, color=INK_)
    ax.text(x3 + 0.05, 1.38, 'residue 57: V (05:01) and D (06:01) → one copy of 57 V', ha='left', va='center',
            fontsize=9.4, color=INK_)
    ax.add_patch(FancyBboxPatch((x3 + 0.05, 0.2), 3.85, 0.92, boxstyle='round,pad=0.02,rounding_size=0.06',
                                facecolor=TAN, edgecolor=TANE, lw=1.3))
    ax.text(x3 + 0.17, 0.66, 'Only alleles of the reference people can be\ncalled; a low posterior marks an uncertain '
            'call;\nno sequencing read is used.', ha='left', va='center', fontsize=9.2, color=INK_, linespacing=1.3)
    fig.subplots_adjust(0, 0, 1, 1)
    save(fig, 'hibag_algorithm')


def hibag_scheme(w=6.3, h=4.6):
    """HIBAG in one picture (user 2026-11-05): a model trained once on people with known HLA types relates SNP
    haplotypes to DQB1 alleles (an ensemble of classifiers on bootstrap samples and random SNP subsets); each of our
    people is called from their array SNPs, and residue 57 is read from the two alleles"""
    st = json.load(open(os.path.join(HERE, 'stats.json')))
    M, HB = st['hibag_model'], st['array_hibag']['full_mainland']['run']
    S.setup_style('slide')
    fig, ax = plt.subplots(figsize=(w, h))
    ax.set_xlim(0, w)
    ax.set_ylim(0, h)
    ax.axis('off')

    def box(x, y, bw, bh, head, lines, hc=NAVY, fill='white', fs=9.4):
        ax.add_patch(plt.Rectangle((x, y), bw, bh, facecolor=fill, edgecolor=hc, lw=1.1, zorder=2))
        ax.add_patch(plt.Rectangle((x, y + bh - 0.32), bw, 0.32, facecolor=hc, edgecolor=hc, lw=1.1, zorder=3))
        ax.text(x + 0.08, y + bh - 0.16, head, ha='left', va='center', fontsize=fs + 0.4, color='white',
                fontweight='bold', zorder=4)
        ax.text(x + 0.08, y + bh - 0.42, '\n'.join(lines), ha='left', va='top', fontsize=fs, color=INK_,
                linespacing=1.3, zorder=4)

    def arrow(x0, y0, x1, y1):
        ax.annotate('', xy=(x1, y1), xytext=(x0, y0), zorder=1,
                    arrowprops=dict(arrowstyle='-|>', color=NAVY, lw=1.3, mutation_scale=13, shrinkA=0, shrinkB=0))

    ax.text(0, h - 0.05, 'Trained once, by the model\'s authors', ha='left', va='top', fontsize=10, fontweight='bold',
            color=GREY_)
    box(0, 2.55, 2.75, 1.55, 'Reference people', [f"{M['training_samples']} people with DQB1 typed", 'by classical methods',
                                                   'and genotyped SNPs', '(East and South Asian;', 'HapMap CHB and JPT)'],
        hc='#6E6E6E')
    box(3.35, 2.55, 2.95, 1.55, f"{M['classifiers']} classifiers", [
        'each: a bootstrap sample of', f"the people and ~{M['snps_per_classifier']:.0f} of the {M['snps']} SNPs;",
        'SNP haplotypes and their DQB1', f"alleles ({M['alleles']} alleles in all)"], hc='#6E6E6E')
    arrow(2.75, 3.3, 3.35, 3.3)
    ax.text(0, 2.3, 'Applied to each of our people', ha='left', va='top', fontsize=10, fontweight='bold', color=GREY_)
    box(0, 0.05, 1.95, 1.95, 'Array SNPs', ['one person,', 'xMHC SNPs typed', 'on the arrays',
                                            f"({HB['model_snps_on_array']} of the {M['snps']}", 'model SNPs)'])
    box(2.2, 0.05, 2.0, 1.95, 'DQB1 call', ['each classifier gives', 'the chance of every', 'allele pair; averaged:',
                                            'best pair and its', 'posterior'])
    box(4.45, 0.05, 1.85, 1.95, 'Residue 57', ['IMGT sequence of', 'each allele: V or', 'not; 57 V copies', '0, 1 or 2'],
        fill=PALE_)
    arrow(1.95, 1.0, 2.2, 1.0)
    arrow(4.2, 1.0, 4.45, 1.0)
    arrow(4.85, 2.55, 3.3, 2.0)
    fig.subplots_adjust(0, 0, 1, 1)
    save(fig, 'hibag_scheme')


def q3_flow(w=12.7, h=2.7):
    """Q3 in one picture (user 2026-11-05; compact for the design slide): the QC'd arrays, the Minimac4-imputed arrays
    (autosomes), the one marker each tested signal is read from, and what is compared in the same people"""
    st = json.load(open(os.path.join(HERE, 'stats.json')))
    AS, MC, DR, HB = st['array_sets'], st['array_match_counts'], st['array_drivers'], st['array_hibag']['full_mainland']
    S.setup_style('slide')
    fig, ax = plt.subplots(figsize=(w, h))
    ax.set_xlim(0, w)
    ax.set_ylim(0, h)
    ax.axis('off')
    FS = 9.6

    def box(x, y, bw, bh, head, lines, fill='white', hc=NAVY, fs=FS):
        ax.add_patch(plt.Rectangle((x, y), bw, bh, facecolor=fill, edgecolor=hc, lw=1.1, zorder=2))
        ax.add_patch(plt.Rectangle((x, y + bh - 0.3), bw, 0.3, facecolor=hc, edgecolor=hc, lw=1.1, zorder=3))
        ax.text(x + 0.08, y + bh - 0.15, head, ha='left', va='center', fontsize=fs + 0.3, color='white',
                fontweight='bold', zorder=4)
        ax.text(x + 0.08, y + bh - 0.38, '\n'.join(lines), ha='left', va='top', fontsize=fs, color=INK_,
                linespacing=1.3, zorder=4)

    def arrow(x0, y0, x1, y1, ls='-'):
        ax.annotate('', xy=(x1, y1), xytext=(x0, y0), zorder=1,
                    arrowprops=dict(arrowstyle='-|>', color=NAVY, lw=1.3, mutation_scale=13, shrinkA=0, shrinkB=0, ls=ls))

    def step(x, txt):
        ax.text(x, h - 0.04, txt, ha='left', va='top', fontsize=10.5, fontweight='bold', color=NAVY)

    fmt = lambda v: f'{v:,}'  # noqa: E731
    y1, bh1 = 0.89, 1.2
    x1, w1 = 0.0, 2.35
    step(x1, '1 · QC\'d arrays')
    box(x1, y1, w1, bh1, 'Genotyped SNPs', ['ASA24v1: the cases', 'JSA24v1: the controls',
                                           f"{fmt(AS['qc']['samples'])} people",
                                           f"{fmt(AS['qc']['variants'])} variants"])
    x2, w2 = 2.95, 2.35
    step(x2, '2 · Imputation')
    box(x2, y1, w2, bh1, 'Imputed arrays', ['from the QC\'d SNPs', 'panel: AGP3K whole genomes',   # round 12, C-32
                                           f"chromosomes {AS['imputed']['chromosomes']}",
                                           f"{fmt(AS['imputed']['variants'])} variants"])
    arrow(x1 + w1, y1 + 0.6, x2, y1 + 0.6)
    # 3 · one marker per signal
    x3, w3 = 5.95, 3.55
    step(x3, '3 · One marker per signal')
    ys, bh3 = [1.73, 0.89, 0.05], 0.72
    rows = [(f"{MC['typed'] + MC['imputed'] + MC['LD proxy'] + MC['none']} GWAS leads",
             ['typed, else imputed (R² ≥ 0.8),', 'else an LD proxy (r² ≥ 0.8)']),
            ('DQB1:57:V', ['HIBAG on the typed xMHC SNPs',
                           f"({HB['run']['model_snps_on_array']} of {HB['run']['model_snps']} model SNPs)"]),
            (f"{DR['genes']} genes, full mainland", [f"{DR['n']} driver variants: their WGS case", 'carriers on the imputed arrays'])]
    for y, (head, lines) in zip(ys, rows):
        box(x3, y, w3, bh3, head, lines, fill=PALE_)          # review round 13 (F-05): one type size (was fs=9.3)
    arrow(x2 + w2, y1 + 1.05, x3, ys[0] + bh3 / 2)
    arrow(x2 + w2, y1 + 0.35, x3, ys[2] + bh3 / 2)
    xs, yb = x1 + w1 / 2, 0.45
    ax.plot([xs, xs, x3 - 0.3, x3 - 0.3], [y1, yb, yb, ys[1] + bh3 / 2], color=NAVY, lw=1.3, ls='--', zorder=1)
    arrow(x3 - 0.3, ys[1] + bh3 / 2, x3, ys[1] + bh3 / 2, ls='--')
    ax.text((xs + x3 - 0.3) / 2, yb - 0.05, 'typed xMHC SNPs, not imputed', ha='center', va='top', fontsize=9.5,
            color=GREY_, style='italic')
    # 4 · what is compared
    x4 = 10.1
    w4 = w - x4 - 0.02
    step(x4, '4 · Same signal?')
    box(x4, 0.89, w4, 1.56, 'Association', ['re-run in the same people,', 'same model, vs. WGS:',
                                           'z compared: contradicted', 'or consistent'])
    box(x4, 0.05, w4, 0.72, 'Carriers', ['WGS case carriers the', 'imputed arrays confirm'])
    arrow(x3 + w3, ys[0] + bh3 / 2, x4, ys[0] + bh3 / 2)
    arrow(x3 + w3, ys[1] + bh3 / 2, x4, ys[1] + bh3 / 2)
    arrow(x3 + w3, ys[2] + bh3 / 2, x4, ys[2] + bh3 / 2)
    fig.subplots_adjust(0, 0, 1, 1)
    save(fig, 'q3_flow')


# ── the four questions: hypothetical values, what an artifact would look like against a real signal ────────────
def fig_questions(w=2.62, h=1.45, fs=9.3):
    # drawn at the size of its cell on the aims slide (w x h inches), so that every text is 9 pt
    real, art = NAVY, OCHRE

    def canvas():
        S.setup_style('slide')
        fig, ax = plt.subplots(figsize=(w, h))
        ax.tick_params(labelsize=fs, length=2.5, pad=1.5)
        return fig, ax

    def legend(ax, loc='upper right'):
        ax.legend(handles=[Line2D([], [], color=real, marker='o', ms=3.5, lw=1.2, label='real'),
                           Line2D([], [], color=art, marker='o', ms=3.5, lw=1.2, label='artifact')],
                  loc=loc, fontsize=fs, frameon=False, handlelength=1.0, handletextpad=0.3, borderaxespad=0.1,
                  labelspacing=0.2)

    def done(fig, ax, name):
        S.despine(ax, 'y')
        fig.tight_layout(pad=0.1)
        save(fig, name)
    # Q1 depth: |z| at the original depth and at the controls' depth
    fig, ax = canvas()
    ax.plot([0, 1], [5.0, 4.85], '-o', color=real, ms=3.5, lw=1.3)
    ax.plot([0, 1], [5.0, 3.1], '-o', color=art, ms=3.5, lw=1.3)
    ax.annotate('', xy=(1.1, 3.1), xytext=(1.1, 5.0), arrowprops=dict(arrowstyle='<->', color=GREY_, lw=0.7))
    ax.text(1.16, 4.05, 'drop\n> 1', fontsize=fs, color=GREY_, va='center')
    ax.set_xticks([0, 1], ['original', "controls' depth"], fontsize=fs)
    ax.set_xlim(-0.25, 1.55)
    ax.set_ylim(2.4, 5.6)
    ax.set_yticks([3, 4, 5])
    ax.set_ylabel('|z|', fontsize=fs, labelpad=4)
    legend(ax, 'lower left')
    done(fig, ax, 'fig_q1')
    # Q2 platform: OR in each case platform against all controls
    fig, ax = canvas()
    x = np.arange(4)
    ax.errorbar(x - 0.13, [1.8, 1.7, 1.9, 1.75], yerr=0.35, fmt='o', ms=3.5, color=real, lw=1.0, capsize=0)
    ax.errorbar(x + 0.13, [2.7, 1.05, 0.95, 1.0], yerr=0.35, fmt='o', ms=3.5, color=art, lw=1.0, capsize=0)
    ax.axhline(1, color=GREY_, lw=0.7, ls=':')
    ax.set_xticks(x, ['T7', 'NovaSeq', '30x', '15x'], fontsize=fs)
    ax.set_ylim(0.4, 3.3)
    ax.set_yticks([1, 2, 3])
    ax.set_ylabel('OR', fontsize=fs, labelpad=1)
    legend(ax)
    done(fig, ax, 'fig_q2')
    # Q3 technology: OR on arrays against OR on WGS, same people
    fig, ax = canvas()
    ax.plot([0.8, 3.0], [0.8, 3.0], color=GREY_, lw=0.7, ls=':')
    ax.scatter([1.4, 2.0, 2.6], [1.35, 2.05, 2.5], s=14, color=real, zorder=3)
    ax.scatter([2.3], [1.05], s=14, color=art, zorder=3)
    ax.annotate('', xy=(2.3, 1.12), xytext=(2.3, 2.22), arrowprops=dict(arrowstyle='->', color=art, lw=0.8))
    ax.set_xlim(0.8, 3.0)
    ax.set_ylim(0.8, 3.0)
    ax.set_xticks([1, 2, 3])
    ax.set_yticks([1, 2, 3])
    ax.set_xlabel('OR, WGS', fontsize=fs, labelpad=1)
    ax.set_ylabel('OR, arrays', fontsize=fs, labelpad=1)
    legend(ax, 'upper left')
    done(fig, ax, 'fig_q3')
    # Q4 control source: OR of all cases against one control source
    fig, ax = canvas()
    x = np.arange(5)
    ax.errorbar(x - 0.13, [1.7, 1.8, 1.65, 1.75, 1.7], yerr=0.3, fmt='o', ms=3.5, color=real, lw=1.0, capsize=0)
    ax.errorbar(x + 0.13, [2.7, 1.05, 1.0, 0.95, 1.05], yerr=0.3, fmt='o', ms=3.5, color=art, lw=1.0, capsize=0)
    ax.axhline(1, color=GREY_, lw=0.7, ls=':')
    ax.set_xticks(x, ['1', '2', '3', '4', '5'], fontsize=fs)
    ax.set_xlabel('Control source', fontsize=fs, labelpad=1)
    ax.set_ylim(0.4, 3.3)
    ax.set_yticks([1, 2, 3])
    ax.set_ylabel('OR', fontsize=fs, labelpad=1)
    legend(ax)
    done(fig, ax, 'fig_q4')


# ── the background: depth and read length of each sequencing group, full mainland ─────────────────────────────
def background_depth(w=4.75, h=3.15):
    # the median, read length and n in three aligned columns right of the axes, so no line runs through a label
    d = tsv('01.depth', 'platform_depth.tsv')
    d = d[d.cohort == 'full_mainland'].set_index('platform')
    order = [('HiSeqX 15x', 'Controls · HiSeq X'), ('DNBSeq-T7 30x', 'Cases · DNBSeq-T7'), ('NovaSeq 30x', 'Cases · NovaSeq'),
             ('DNBSeq-G400RS 30x', 'Cases · G400RS 30x'), ('DNBseq-G400RS 15x', 'Cases · G400RS 15x')]
    S.setup_style('slide')
    fig, ax = plt.subplots(figsize=(w, h))
    n = len(order)
    tx = ax.get_yaxis_transform()
    for i, (p_, lab) in enumerate(order):
        r = d.loc[p_]
        y = n - 1 - i
        col = CTRL_C if p_.startswith('HiSeq') else CASE_C
        ax.plot([r.depth_q25, r.depth_q75], [y, y], color=col, lw=2.4, solid_capstyle='butt')
        ax.plot(r.depth_median, y, 'o', ms=8, color=col, zorder=3)
        ax.text(1.04, y, f'{r.depth_median:.1f}x', transform=tx, fontsize=9.5, va='center', ha='left',
                color=GREY_ if p_.startswith('HiSeq') else col, fontweight='bold')   # review round 13 (F-06): legible grey
        ax.text(1.27, y, f'{r.readlen_median:.0f} bp', transform=tx, fontsize=9.5, va='center', ha='left', color=INK_)
        ax.text(1.52, y, f'{int(r.n):,}', transform=tx, fontsize=9.5, va='center', ha='left', color=INK_)
    for xx, t in ((1.04, 'depth'), (1.27, 'reads'), (1.52, 'n')):
        ax.text(xx, n - 0.45, t, transform=tx, fontsize=9.5, va='bottom', ha='left', fontweight='bold', color=GREY_)
    ax.axvline(d.loc['HiSeqX 15x'].depth_median, color=NAVY, lw=0.9, ls='--', zorder=1)
    ax.set_yticks(range(n)[::-1], [o[1] for o in order], fontsize=10)
    ax.set_ylim(-0.5, n - 0.45)
    ax.set_xlim(10, 38)
    ax.set_xticks([10, 15, 20, 25, 30, 35])
    ax.set_xlabel('Effective depth (x): median and IQR', fontsize=10)
    ax.tick_params(axis='x', labelsize=9.5)
    S.despine(ax, 'x')
    ax.grid(axis='y', visible=False)
    fig.subplots_adjust(right=0.64)
    save(fig, 'background_depth')


# ── Q1 emulation as "before → after, against the controls" (user 2026-11-04: slides 10–11 more direct): a, the depth
#    of each case group; b, the share of genotypes failing a QC rule; one grammar for both panels ───────────────────
def emulation_dumbbell(w=8.0, h=3.6):
    """Q1's positive control, full mainland, medians only: each down-sampled case group from all reads (open) to the
    controls' depth (filled), against the controls' median; no value labels (the numbers are in the slide's table)
    and no G400RS 15x row (not down-sampled; the caption says so): user 2026-11-04, 'less clutter'"""
    st = json.load(open(os.path.join(HERE, 'stats.json')))
    dps = pd.read_csv(os.path.join(R, '01.depth', 'depth_per_sample.tsv'), sep='\t')
    dps = dps[dps.in_full_mainland == 1]
    fr = pd.read_csv(os.path.join(R, '01.depth', 'fractions.all.tsv'), sep='\t')
    qf = st['qc_fail']['full_mainland']
    groups = [('DNBSeq-T7 30x', 'DNBSeq-T7'), ('NovaSeq 30x', 'NovaSeq'), ('DNBSeq-G400RS 30x', 'G400RS 30x')]
    n = len(groups)
    ctrl = dps[dps.platform == 'HiSeqX 15x'].probe_depth
    S.setup_style('slide')
    fig, axes = plt.subplots(1, 2, figsize=(w, h), sharey=True, gridspec_kw={'wspace': 0.12})
    MS = 10

    def pair(ax, y, x0, x1):
        ax.plot(x0, y, 'o', ms=MS, mfc='white', mec=CASE_C, mew=1.6, zorder=4)
        ax.plot(x1, y, 'o', ms=MS, color=CASE_C, zorder=5)
        if abs(x1 - x0) > 1e-9:
            ax.annotate('', xy=(x1, y), xytext=(x0, y), zorder=2,
                        arrowprops=dict(arrowstyle='-|>', color=CASE_C, lw=1.4, shrinkA=7, shrinkB=7, mutation_scale=12))

    def ref(ax, x, txt):
        ax.axvline(x, color=GREY_, ls='--', lw=1.2, zorder=1)
        ax.text(x, n + 0.55, txt, fontsize=10.5, color=GREY_, ha='center', va='center',
                bbox=dict(fc='white', ec='none', pad=1.0))

    ax = axes[0]
    med = ctrl.median()
    ref(ax, med, f'controls {med:.1f}x')
    for i, (p, lab) in enumerate(groups):
        pair(ax, n - i, dps[dps.platform == p].probe_depth.median(), fr[fr.platform == p].expected_depth.median())
    ax.set_xlim(12, 37)
    ax.set_xlabel('Effective depth (x)', fontsize=11)
    ptitle(ax, 'a', 'Depth', size=11.5)

    ax = axes[1]
    cq = 100 * qf['controls']
    ref(ax, cq, f'controls {cq:.1f}%')
    for i, (p, lab) in enumerate(groups):
        pair(ax, n - i, 100 * qf[f'baseline|{p}'], float(np.median([100 * x for x in qf[f'reps|{p}']])))
    ax.set_xlim(-0.2, 3.2)
    ax.set_xlabel('Genotypes failing a QC rule (%)', fontsize=11)
    ptitle(ax, 'b', 'Missing genotypes', size=11.5)
    axes[0].set_yticks([n - i for i in range(n)], [lab for _, lab in groups], fontsize=11)
    axes[0].set_ylim(0.5, n + 0.85)
    for a in axes:
        a.tick_params(axis='x', labelsize=10.5)
        a.grid(axis='y', visible=False)
        S.despine(a, 'x')
        a.spines['left'].set_visible(False)
        a.tick_params(axis='y', length=0, pad=4)
    hs = [Line2D([], [], marker='o', ls='', ms=MS, mfc='white', mec=CASE_C, mew=1.6, label='all reads'),
          Line2D([], [], marker='o', ls='', ms=MS, color=CASE_C, label="at the controls' depth")]
    axes[1].legend(handles=hs, loc='lower right', bbox_to_anchor=(1.0, 1.08), ncol=2, frameon=False, fontsize=10.5,
                   handletextpad=0.3, columnspacing=1.4)
    save(fig, 'emulation_dumbbell')

# ── Q1's signals on the |z| scale, medians only, decluttered: where each signal was and where it is at the controls'
#    depth; one panel per cohort; no thresholds, no ranges (user, 2026-11-04) ──────────────────────────────────────
def ds_zscale(w=8.6, h=5.0):
    """every signal on the |z| scale, three cohorts: the original genotypes and the replicates at the controls'
    depth; GWAS leads, the HLA residue DQB1:57:V, and the gene sets (the test that called the set: the smaller
    full-mainland P of CMC and SKAT-O)"""
    st = json.load(open(os.path.join(HERE, 'stats.json')))
    reps = replicate_names()
    dsg = rtsv('ds.gwas.tsv')
    dsg = dsg[(dsg.imputation == 'best_guess') & (dsg.role == 'Lead')]
    dsh = rtsv('ds.hla_dqb1_57v.tsv')
    dsn = rtsv('ds.genes.tsv')

    def zp(p):
        return float(norm.isf(float(p) / 2))
    order = _matrix_order()
    leads = [s_ for s_ in order if not (s_.startswith('DQB1') or ' (' in s_)]
    hla = [s_ for s_ in order if s_.startswith('DQB1')]
    genes = [s_ for s_ in order if ' (' in s_]
    test = {}
    for g in genes:
        gene, sr = g[:-1].split(' (')
        x = dsn[(dsn.gene == gene) & (dsn.stratum == sr) & (dsn.cohort == 'full_mainland')]
        test[g] = x.sort_values('pvalue_original').iloc[0].method        # the test that called the set

    def row(sig, c):
        """original |z| and the replicates' |z| of one signal in one cohort"""
        if sig in leads:
            d = dsg[(dsg.unit == sig) & (dsg.cohort == c)].iloc[0]
            po, pr = d.p_original, [d[f'p_{r}'] for r in reps]
        elif sig == 'DQB1:57:V':
            d = dsh[dsh.cohort == c].iloc[0]
            po, pr = d.p_original, [d[f'p_{r}'] for r in reps]
        else:
            gene, sr = sig[:-1].split(' (')
            d = dsn[(dsn.gene == gene) & (dsn.stratum == sr) & (dsn.cohort == c) & (dsn.method == test[sig])].iloc[0]
            po, pr = d.pvalue_original, [d[f'pvalue_{r}'] for r in reps]
        return zp(po), [zp(p) for p in pr]
    # rows top to bottom: a header row per block, then its signals
    blocks = [('GWAS leads', leads), ('HLA', hla), ('Gene sets', genes)]
    ypos, labels, heads, y = {}, [], [], 0
    for name, sigs in blocks:
        heads.append((y, name)); y -= 1
        for s_ in sigs:
            ypos[s_] = y; y -= 1
    vals = {(s_, c): row(s_, c) for s_ in ypos for c in CK}
    for s_ in leads + hla:                                   # the figure must read as the matrix does
        for c in CK:
            z0, zr = vals[s_, c]
            if abs(st['long'][s_]['depth'][c]['value'] - (z0 - float(np.median(zr)))) > 1e-6:
                sys.exit(f'ds_zscale: {s_} {c}: the drop differs from the matrix reading')
    allz = [v for (z0, zr) in vals.values() for v in [z0] + zr]
    lo, hi = min(allz) - 0.15, max(allz) + 0.15
    S.setup_style('slide')
    fig, axes = plt.subplots(1, 3, figsize=(w, h), sharey=True, gridspec_kw={'wspace': 0.07})
    for ax, c in zip(axes, CK):
        for yy, _ in heads[1:]:
            ax.axhline(yy + 0.5, color=S.NEUTRAL_D, lw=0.7, zorder=0)
        for s_, yy in ypos.items():
            z0, zr = vals[s_, c]
            ax.plot(z0, yy, 'o', ms=8, mfc='white', mec=COH_COL[c], mew=1.4, zorder=4)
            ax.plot(float(np.median(zr)), yy, 'o', ms=5, color=COH_COL[c], zorder=5)   # the median only
        ax.set_xlim(lo, hi)
        ax.set_xticks([t for t in (3, 4, 5) if lo <= t <= hi])
        ax.tick_params(axis='x', labelsize=10.5)
        ptitle(ax, 'abc'[CK.index(c)], CNAME[c].split()[0], size=11)   # review round 12 (C-15): lettered panels
        ax.title.set_color(COH_COL[c])
        ax.grid(axis='y', visible=False)
        S.despine(ax, 'x')
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis='y', length=0, pad=3)
    ticks = [(yy, nm, 'head') for yy, nm in heads] + [(yy, s_, 'sig') for s_, yy in ypos.items()]
    ticks.sort(key=lambda t: -t[0])

    def lab(s_):
        if s_ in genes:
            return _gene_label(s_)                           # the test is named in the caption
        return s_
    axes[0].set_ylim(y + 0.5, 0.6)
    axes[0].set_yticks([t[0] for t in ticks], [t[1] if t[2] == 'head' else '' if t[1] in genes else lab(t[1])
                                               for t in ticks], fontsize=10.5)
    gene_ticklabels(axes[0], [(t[0], t[1].split(' (')[0], 'M/H') for t in ticks if t[2] == 'sig' and t[1] in genes],
                    10.5, pad_pt=3.0)
    for t, (yy, nm, kind) in zip(axes[0].get_yticklabels(), ticks):
        if kind == 'head':
            t.set_fontweight('bold'); t.set_color(NAVY); t.set_fontsize(11)
        elif nm in leads:
            t.set_fontstyle('italic')
    hs = [Line2D([], [], marker='o', ls='', ms=8, mfc='white', mec=INK_, mew=1.4, label='original'),
          Line2D([], [], marker='o', ls='', ms=5, color=INK_, label="controls' depth (median of 3 replicates)")]
    fig.supxlabel('|z|', fontsize=11, y=0.02)
    axes[1].legend(handles=hs, loc='lower center', bbox_to_anchor=(0.5, 1.07), ncol=2, frameon=False, fontsize=10.5,
                   handletextpad=0.4, columnspacing=1.6)
    save(fig, 'ds_zscale')

def v57_pg(w=5.2, h=4.45):
    v = rtsv('ds.hla_dqb1_57v.tsv')
    reps = replicate_names()
    S.setup_style('slide')
    fig, ax = plt.subplots(figsize=(w, h))
    for j, c in enumerate(CK):
        r = v[v.cohort == c].iloc[0]
        pg = [r[f'p_over_gec_{a}'] for a in reps]
        ax.plot([j + 0.2] * 2, [min(pg), max(pg)], color=COH_COL[c], lw=2.4)
        ax.plot(j - 0.2, r.p_over_gec_original, 'o', ms=9, mfc='white', mec=COH_COL[c], mew=1.6)
        ax.plot(j, r.p_over_gec_baseline, 's', ms=8, color=COH_COL[c])
    ax.axhline(1, color=NAVY, ls='--', lw=1.0)
    ax.text(2.45, 1.04, 'GEC threshold', ha='right', va='bottom', fontsize=9.5, color=NAVY)
    ax.set_yscale('log')
    ax.set_yticks([0.1, 0.2, 0.4, 0.6, 1.0], ['0.1', '0.2', '0.4', '0.6', '1'])
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_ylim(0.1, 1.3)
    ax.set_xticks(range(3), [CNAME[c].split()[0] for c in CK], fontsize=10.5)
    for tl, c in zip(ax.get_xticklabels(), CK):
        tl.set_color(COH_COL[c])
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylabel('P / GEC threshold (below 1: significant)', fontsize=10)
    ax.tick_params(axis='y', labelsize=9.5)
    hs = [Line2D([], [], marker='o', ls='', mfc='white', mec=GREY_, label='original'),
          Line2D([], [], marker='s', ls='', color=GREY_, label='f = 1 baseline'),
          Line2D([], [], color=GREY_, lw=2.2, label=f'{len(reps)} replicates')]
    ax.legend(handles=hs, fontsize=9.5, frameon=False, ncol=3, loc='lower left', bbox_to_anchor=(0, 1.0),
              handletextpad=0.3, columnspacing=1.0)
    S.despine(ax, 'y')
    fig.tight_layout()
    save(fig, 'v57_pg')


# ── DQB1:57:V: the down-sampled cases whose 57 V dosage differs between arms ──
def replicate_names():
    ap = tsv('03.calls', 'downsampling', 'arm_plan.tsv')
    return sorted(a for a in ap[ap.role == 'deep'].arm.unique() if re.fullmatch(r'r\d+', a))


PSHORT = {'DNBSeq-T7 30x': 'DNBSeq-T7', 'NovaSeq 30x': 'NovaSeq', 'DNBSeq-G400RS 30x': 'G400RS 30x',
          'DNBseq-G400RS 15x': 'G400RS 15x'}


def v57_trace(w=6.6, row=0.42, lab_w=1.5, top=0.27, bottom=0.04, legend=False):
    arms = ['original', 'baseline'] + replicate_names()
    h = rtsv('followup.hla_v57.tsv')
    h = h[(h.cohort == 'full_mainland') & h.arm.isin(arms)]
    # the follow-up's own criterion (signal_followup.hla_v57), within the arm plan: the dosage differs between two
    # arms, a missing dosage counted as a value; collect_stats.py applies the same one (v57_traced)
    samples = [s for s, x in h.groupby('sample_id', sort=True) if x.v57_dosage.astype(str).nunique() > 1]
    out = os.path.join(OUT, 'v57_trace.png')
    if not samples:                                         # nothing to draw: the slide says so in words
        if os.path.exists(out):
            os.remove(out)
        return
    # rows are named by the case's platform, never by its sample ID (the deck shows no IDs)
    plat = pd.read_csv(os.path.join(R, '00.samples', 'sample_table.tsv'), sep='\t', dtype=str).set_index('sample_id').platform
    names = [f'{PSHORT.get(plat.get(s, ""), "down-sampled")} case' for s in samples]
    names = [f'{nm} {names[:i].count(nm) + 1}' if names.count(nm) > 1 else nm for i, nm in enumerate(names)]
    S.setup_style('slide')
    n = len(samples)
    H = top + row * n + bottom
    fig = plt.figure(figsize=(w + lab_w - 0.95, H))
    W = w + lab_w - 0.95
    ax = fig.add_axes([lab_w / W, bottom / H, (W - lab_w - 0.02) / W, row * n / H])
    for i, s in enumerate(samples):
        x = h[h.sample_id == s].set_index('arm')
        y = n - 1 - i
        for j, ar in enumerate(arms):
            r = x.loc[ar] if ar in x.index else None
            ch = r.change_vs_baseline if r is not None else np.nan
            col = OCHRE if ch > 0 else NAVY if ch < 0 else S.NEUTRAL if np.isfinite(ch) else 'white'
            ax.add_patch(plt.Rectangle((j + 0.02, y + 0.05), 0.96, 0.9, facecolor=col, edgecolor=S.NEUTRAL_D, lw=0.5))
            tc = 'white' if np.isfinite(ch) and ch != 0 else S.INK
            call = (str(r.dqb1_call).replace('DQB1*', '').strip('/') if r is not None and isinstance(r.dqb1_call, str)
                    else 'not typed')
            v = 'no dosage' if r is None or not np.isfinite(r.v57_dosage) else f'{r.v57_dosage:.0f} V'
            conf = f' · {r.confidence}' if r is not None and isinstance(r.confidence, str) else ''
            ax.text(j + 0.5, y + 0.68, call, fontsize=9, ha='center', va='center', color=tc)
            ax.text(j + 0.5, y + 0.3, v + conf, fontsize=9, ha='center', va='center', color=tc)
    ax.set_xlim(0, len(arms))
    ax.set_ylim(0, n)
    ax.set_xticks(np.arange(len(arms)) + 0.5, arms, fontsize=10)
    ax.set_yticks(np.arange(n)[::-1] + 0.5, names, fontsize=10)
    ax.xaxis.tick_top()
    ax.tick_params(length=0, pad=3)
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    used = {('as the baseline', S.NEUTRAL), ('no dosage', 'white')}
    hs = [(lb, c) for lb, c in (('as the baseline', S.NEUTRAL), ('one V copy more', OCHRE), ('one V copy fewer', NAVY),
                                ('no dosage', 'white'))]
    chs = h[h.sample_id.isin(samples)].change_vs_baseline
    if (chs > 0).any():
        used.add(('one V copy more', OCHRE))
    if (chs < 0).any():
        used.add(('one V copy fewer', NAVY))
    hs = [(lb, c) for lb, c in hs if (lb, c) in used]       # the legend lists only the cell kinds drawn
    if legend:
        ax.legend([plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor=S.NEUTRAL_D, lw=0.5) for _, c in hs],
                  [lb for lb, _ in hs], loc='upper center', bbox_to_anchor=(0.5, -0.04), ncol=len(hs), frameon=False,
                  fontsize=9, handlelength=1.2, columnspacing=1.2)
    save(fig, 'v57_trace')


# ── the appendix depth figure in the deck's grammar (results/figures/01.depth/depth_platforms colors two platforms
#    in the cohort colors and the cases blue): per-genome depth, fraction kept, read length ─────────────────────
def depth_platforms(w=12.7, h=4.4):
    d = tsv('01.depth', 'depth_per_sample.tsv')
    d = d[d.in_full_mainland == 1]
    fr = tsv('01.depth', 'fractions.all.tsv')
    target = float(json.load(open(os.path.join(R, '01.depth', 'depth_summary.json')))['target_depth'])
    order = [('HiSeqX 15x', 'HiSeq X\n(controls)'), ('DNBSeq-T7 30x', 'DNBSeq-T7'), ('NovaSeq 30x', 'NovaSeq'),
             ('DNBSeq-G400RS 30x', 'G400RS\n30x'), ('DNBseq-G400RS 15x', 'G400RS\n15x')]
    pcol = {'DNBSeq-T7 30x': PLAT_COL[0], 'NovaSeq 30x': PLAT_COL[1], 'DNBSeq-G400RS 30x': PLAT_COL[2]}
    S.setup_style('slide')
    fig, axes = plt.subplots(1, 3, figsize=(w, h), gridspec_kw={'width_ratios': [1.3, 0.95, 1.2]})
    ax = axes[0]
    rng = np.random.default_rng(1)
    for i, (p_, lab) in enumerate(order):
        x = d[d.platform == p_].probe_depth.values
        col = CTRL_C if p_.startswith('HiSeq') else CASE_C
        ax.scatter(i + rng.uniform(-0.28, 0.28, len(x)), x, s=4, color=col, alpha=0.45 if len(x) > 500 else 0.7, lw=0)
        ax.plot([i - 0.32, i + 0.32], [np.median(x)] * 2, color=S.INK, lw=2.0, zorder=3)
    ax.axhline(target, color=NAVY, ls='--', lw=1.0)
    ax.text(2.5, target + 0.6, f'D* = {target:.1f}x', color=NAVY, fontsize=10, ha='center', va='bottom')
    ax.set_xticks(range(len(order)), [o[1] for o in order], fontsize=10)
    ax.set_ylabel('Effective depth (x)', fontsize=10.5)
    ax.tick_params(axis='y', labelsize=10)
    ax.grid(axis='x', visible=False)
    S.despine(ax, 'y')
    ptitle(ax, 'a', 'Effective depth per genome', size=11)
    ax = axes[1]
    bins = np.linspace(0.35, 1.0, 27)
    for p_, lab in (('DNBSeq-G400RS 30x', 'G400RS 30x'), ('NovaSeq 30x', 'NovaSeq'), ('DNBSeq-T7 30x', 'DNBSeq-T7')):
        f = fr[fr.platform == p_].fraction.values
        ax.hist(f, bins=bins, color=pcol[p_], alpha=0.95, label=f'{lab} ({len(f)})')
    ax.set_xlabel('Fraction of read pairs kept', fontsize=10.5)
    ax.set_ylabel('Down-sampled cases', fontsize=10.5)
    ax.tick_params(labelsize=10)
    ax.legend(fontsize=10, frameon=False, loc='upper left')
    S.despine(ax, 'y')
    ptitle(ax, 'b', 'Fraction kept', size=11)
    ax = axes[2]
    xs = np.arange(len(order))
    ax.bar(xs - 0.19, [d[d.platform == p_].readlen_median.median() for p_, _ in order], width=0.38, color='#595959',
           label='measured')
    ax.bar(xs + 0.19, [d[d.platform == p_].readlen_implied.median() for p_, _ in order], width=0.38, color='#BFBFBF',
           label="implied by the sheet's depth")
    ax.set_xticks(xs, ['HiSeq X\n(controls)', 'DNBSeq-\nT7', 'NovaSeq', 'G400RS\n30x', 'G400RS\n15x'], fontsize=10)
    ax.set_ylabel('Read length (bp)', fontsize=10.5)
    ax.tick_params(axis='y', labelsize=10)
    ax.set_ylim(0, 190)
    ax.legend(fontsize=10, frameon=False, loc='upper center')
    ax.grid(axis='x', visible=False)
    S.despine(ax, 'y')
    ptitle(ax, 'c', 'Read-length audit', size=11)
    fig.tight_layout(w_pad=1.2)
    save(fig, 'depth_platforms')



# ── the two flagged leads (user 2026-11-06: "在27和28，我需要你查看cram原文件来判明问题"; restructured the same day:
# "26～29 重构 … 只要能直观，系统地说明我们的问题就行") ──────────────────────────────────────────────────────────────
def _cram_table():
    st = json.load(open(os.path.join(HERE, 'stats.json')))
    d = pd.read_csv(os.path.join(HERE, 'cram_check', 'cram_scan.tsv'), sep='\t', dtype={'sample_id': str})
    return st, d[d.in_full_mainland == 1].copy()


def _wilson(k, n, z=1.959964):
    """95% Wilson interval of a proportion k / n"""
    p, d = k / n, 1 + z * z / n
    c, h = (p + z * z / (2 * n)) / d, z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def _wilson_af(p, n):
    """95% Wilson interval of an allele frequency p over n alleles (user 2026-11-07: both flag figures use Wilson's
    interval; it replaced a Wald interval here); p * n must be a whole allele count"""
    k = round(p * n)
    if abs(p * n - k) > 1e-6:
        raise ValueError(f'allele frequency {p} over {n} alleles is not a whole count')
    return _wilson(k, n)


def _calls_reads(ax, rows, xlim, xticks, xlabel, ref, fmt, seps=()):
    """the two flag figures share one design (user 2026-11-07: "26～29页重构我们不需要这么复杂，就看看这两个信号，说明
    其被flag的来源"): one row per group, its value in the calls (open, 95% CI) above its value in the reads (filled,
    95% CI): open = the original analysis, filled = the check, as on the Q1, Q2 and Q3 figures (review round 9, D-01);
    dashed: the controls in the reads; bold: the cases without T7, the reading that raised the flag; shaded (in regular
    weight, review round 12: one emphasis each, as the captions key them): the group the flag comes from, banded across
    the whole row by _band_rows once the layout is final; values at the right"""
    n = len(rows)
    for i, r in enumerate(rows):
        y = n - 1 - i
        col = r['color']
        ax.plot([r['c_lo'], r['c_hi']], [y + 0.14] * 2, color=col, lw=1.0, alpha=0.8, zorder=2, solid_capstyle='butt')
        ax.scatter([r['c']], [y + 0.14], s=36, facecolor='white', edgecolor=col, lw=1.3, zorder=3)
        ax.plot([r['r_lo'], r['r_hi']], [y - 0.16] * 2, color=col, lw=1.3, zorder=2, solid_capstyle='butt')
        ax.scatter([r['r']], [y - 0.16], s=36, color=col, zorder=3, lw=0)
        ax.text(1.02, y, f"{fmt(r['c'])} | {fmt(r['r'])}", transform=ax.get_yaxis_transform(), va='center', ha='left',
                fontsize=9.5, color=INK_, fontweight='bold' if r.get('bold') else 'normal')
    ax.text(1.02, n - 0.45, 'calls ○ | reads ●', transform=ax.get_yaxis_transform(), va='bottom', ha='left',
            fontsize=9.5, color=GREY_)
    for k in seps:
        ax.axhline(n - 1 - k - 0.5, color=S.NEUTRAL_D, lw=0.6, zorder=1)
    ax.axvline(ref, color=GREY_, lw=1.0, ls='--', zorder=1)
    ax.set_yticks(range(n)[::-1], [r['label'] for r in rows], fontsize=9.5)
    for t, r in zip(ax.get_yticklabels(), rows):
        if r.get('bold'):
            t.set_fontweight('bold')
    ax.set_xlim(*xlim)
    ax.set_xticks(xticks)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlabel(xlabel, fontsize=9.5)
    ax.tick_params(axis='x', labelsize=9.5)
    ax.tick_params(axis='y', length=0)
    ax.grid(axis='y', visible=False)


def _band_rows(fig, rows):
    """review round 9 (D-04): the SIG band across the whole row (tick label, plot and values), as the key rows of the
    matrix and the tables; drawn after the layout is final. Review round 12 (D-06): one x-extent for every panel of a
    figure (the union of their labels and values), so the bands of panels a and b, and of both flag figures, line up"""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    axs = [a for a, _ in rows]
    x0 = min(t.get_window_extent(r).x0 for a in axs for t in a.get_yticklabels() if t.get_text())
    x1 = max([t.get_window_extent(r).x1 for a in axs for t in a.texts] + [a.get_window_extent(r).x1 for a in axs])
    for a, y in rows:
        inv = a.transAxes.inverted()
        f0, f1 = inv.transform((x0 - 4, 0))[0], inv.transform((x1 + 4, 0))[0]
        a.add_patch(plt.Rectangle((f0, y - 0.5), f1 - f0, 1.0, transform=a.get_yaxis_transform(), facecolor=SIG, lw=0,
                                  zorder=0, clip_on=False))


def cram_adarb2(w=5.9, h=4.6):
    """ADARB2, full mainland: a, carriers (% of samples) in the calls (before QC) and in the reads (>= 4 signature
    reads), by group, with the cases without T7 (the reading that raised the flag); b, the HaplotypeCaller census
    (cram_check/hc_census): the share of each group whose candidate haplotypes carry the deletion, 95% Wilson CI (the
    14-genome pilot of cram_check/hc_debug drew one square per genome). Drawn at the flag slides' common width (review
    round 9, D-03)"""
    st, d = _cram_table()
    k = st['cram']['min_reads']
    called, read = d.ad_call.isin(['0/1', '1/1']), d.ad_alt >= k
    t7 = d.platform == 'DNBSeq-T7 30x'
    grp = [('Controls, HiSeq X', d.status == 'control', CTRL_C, {}),
           ('All cases', d.status == 'case', CASE_C, {}),
           ('T7 cases', t7, CASE_C, {'band': True}),
           ('Cases without T7', (d.status == 'case') & ~t7, CASE_C, {'bold': True})]
    rows = []
    for lab, m, col, o in grp:
        nn = int(m.sum())
        kc, kr = int((m & called).sum()), int((m & read).sum())
        lc, hc = _wilson(kc, nn)
        lr, hr = _wilson(kr, nn)
        rows.append({'label': f'{lab} ({nn:,})', 'c': 100 * kc / nn, 'c_lo': 100 * lc, 'c_hi': 100 * hc,
                     'r': 100 * kr / nn, 'r_lo': 100 * lr, 'r_hi': 100 * hr, 'color': col, **o})
    S.setup_style('slide')
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(w, h), gridspec_kw={'height_ratios': [4.0, 2.6]})
    _calls_reads(ax, rows, (-0.45, 16), [0, 5, 10, 15], 'Calls and reads, 95% Wilson CI', rows[0]['r'],
                 lambda x: f'{x:.1f}', seps=(1,))
    ptitle(ax, 'a', 'Carriers (% of samples)', size=10)   # user 2026-11-07
    # b: the HaplotypeCaller census (cram_check/hc_census, user 2026-11-07: every carrier in the reads or the calls, and
    # 20 random non-carriers per platform): the share whose candidate haplotypes carry the deletion, 95% Wilson CI
    # review round 12: the first two rows are the carriers in the reads (>= min_reads signature reads) on 150-bp reads,
    # split by the calls (190 + 69 = 259), so the panel shows the like-for-like comparison of step 4 (190 of 259 against
    # 0 of 32 on T7; A-01, A-02, B-01, D-05); the census's one genome called with fewer signature reads than the bar is
    # counted in the notes. The rows that are not about T7 are drawn in ink: they mix cases and controls (D-01, A-15)
    H = st['cram']['hc_census']['groups']
    p150 = st['cram']['adarb2']['full_mainland']['p150']
    if H['called']['n'] != p150['called_read'] or H['called']['n'] + H['missed']['n'] != p150['read']:
        raise SystemExit('cram_adarb2: panel b rows are no longer the read carriers on 150-bp reads')
    brow = [(f"Carriers called, 150-bp reads ({H['called']['n']:,})", H['called'], INK_, False),
            (f"Carriers missed, 150-bp reads ({H['missed']['n']:,})", H['missed'], INK_, False),
            (f"T7 carriers, 100-bp reads ({H['t7']['n']:,})", H['t7'], CASE_C, True),
            (f"Non-carriers, random ({H['non']['n']:,})", H['non'], INK_, False)]
    for i, (lab, g, col, key) in enumerate(brow):
        y = len(brow) - 1 - i
        lo, hi = _wilson(g['del'], g['n'])
        bx.plot([100 * lo, 100 * hi], [y, y], color=col, lw=1.3, zorder=2, solid_capstyle='butt')
        bx.scatter([100 * g['del'] / g['n']], [y], s=36, color=col, zorder=3, lw=0)
        bx.text(1.02, y, f"{g['del']} of {g['n']}", transform=bx.get_yaxis_transform(), va='center', ha='left',
                fontsize=9.5, color=INK_)
    bx.set_xlim(-3, 110)                     # review round 13 (F-03): the 100% point clear of its value text (was 103)
    bx.set_xticks([0, 25, 50, 75, 100])
    bx.set_ylim(-0.6, len(brow) - 0.4)
    bx.set_yticks(range(len(brow))[::-1], [b[0] for b in brow], fontsize=9.5)
    bx.tick_params(axis='y', length=0)
    bx.tick_params(axis='x', labelsize=9.5)
    bx.grid(axis='y', visible=False)
    bx.set_xlabel('HaplotypeCaller re-run, 95% Wilson CI', fontsize=9.5)
    ptitle(bx, 'b', 'Deletion reconstructed (% of samples)', size=10)
    fig.tight_layout(h_pad=0.9)
    _band_rows(fig, [(ax, 1), (bx, 1)])   # the T7 rows (the census adds the non-carriers below panel b's)
    save(fig, 'cram_ADARB2')


def cram_linc(w=6.0, h=4.55):
    """LINC01146, full mainland: the A frequency of each group in the calls (before QC) and in the reads (genomes with
    >= 8 reads, mapping and base quality >= 20), with 95% Wilson CI: all cases, DNBSeq-T7, the cases without T7 (the reading
    that raised the flag) and their three platforms; the genotypes hold in every group, and the G400RS 15x cases sit
    below the controls"""
    st, d = _cram_table()
    CL = st['cram']['linc']['full_mainland']
    grp = [('Controls, HiSeq X', 'controls', CTRL_C, {}),
           ('All cases', 'cases', CASE_C, {}),
           ('T7 cases', 'T7', CASE_C, {}),
           ('Cases without T7', 'other cases', CASE_C, {'bold': True}),
           ('NovaSeq cases', 'NovaSeq', CASE_C, {}),
           ('G400RS 30x cases', 'G400RS 30x', CASE_C, {}),
           ('G400RS 15x cases', 'G400RS 15x', CASE_C, {'band': True})]
    rows = []
    for lab, key, col, o in grp:
        g = CL[key]
        lc, hc = _wilson_af(g['af_called'], 2 * g['n'])
        lr, hr = _wilson_af(g['af_reads'], 2 * g['depth_ok'])
        rows.append({'label': f"{lab} ({g['n']:,})", 'c': g['af_called'], 'c_lo': lc, 'c_hi': hc, 'r': g['af_reads'],
                     'r_lo': lr, 'r_hi': hr, 'color': col, **o})
    S.setup_style('slide')
    fig, ax = plt.subplots(figsize=(w, h))
    _calls_reads(ax, rows, (0.36, 0.84), [0.4, 0.5, 0.6, 0.7, 0.8], 'Calls and reads, 95% Wilson CI',
                 CL['controls']['af_reads'], lambda x: f'{x:.2f}', seps=(1, 3))
    ax.set_title('Frequency of allele A', fontsize=10, loc='center')
    fig.tight_layout()
    _band_rows(fig, [(ax, 0)])
    save(fig, 'cram_LINC01146')


if __name__ == '__main__':
    forest('platform', 'full_mainland')                     # narrow and intermediate: the Q table of the appendix
    forest('sources', 'full_mainland')
    matrix()
    loo_scatter('lopo')                                     # Q2: platform left out
    loo_regions()                                           # Q2: how a leave-out reading is judged
    loo_scatter('sources_loo')                              # Q4: source left out
    fig_questions()
    background_depth()
    emulation_dumbbell()
    ds_zscale()
    array_z(w=7.7, h=4.75, fs=10)          # z on WGS and arrays (F) and their difference (N / I / F)
    array_genes(w=5.5, h=3.4, fs=10)
    q3_flow()
    hibag_algorithm()
    # HLA-HD as hla.typing.v2 presents it (its slide 'The HLA-HD Typing Algorithm'): the same figure, copied unchanged
    shutil.copyfile(os.path.join(HERE, '..', '..', '..', 'hla.typing.v2', 'deck', 'figures', 'concept_hlahd.png'),
                    os.path.join(OUT, 'concept_hlahd.png'))
    # freq_ADARB2, freq_LINC01146: retired with the restructured flag slides (user 2026-11-06, "26～29 重构"); the
    # frequencies by group stay in the notes, and the reads figures below carry each flag
    cram_adarb2()                                           # the two flagged leads: calls against reads
    cram_linc()
    # depth_platforms(): retired with its appendix slide (user 2026-11-05)
    # figures the deck no longer shows: the former crops, the P-scale lead view, the narrow and intermediate forests and
    # array panels, the two-part carrier figures
    stale = [f for f in os.listdir(OUT) if re.match(r'(lead_frequencies|gene_carriers)\.(top|bottom)\.png$', f) or
             re.match(r'(hibag_scheme|loo_scatter|ds_leads|emulation_depth|ds_change|ds_delta|gene_carriers_cohorts\.[12]|(forest_platform|forest_sources|array_panels)\.'
                      r'(narrow|intermediate)_mainland)\.png$', f) or
             # the slides removed on 2026-11-05 (dossiers of unflagged signals, frequency and carrier appendix)
             re.match(r'(v57_trace|v57_pg|carriers_\w+|lead_frequencies_cohorts\.[12]|gene_carriers_cohorts|depth_platforms|array_panels\.\w+)\.png$', f) or
             # the frequency panels of the flag slides before the restructure (user 2026-11-06)
             re.match(r'freq_(ADARB2|LINC01146)\.png$', f)]
    for f in stale:
        os.remove(os.path.join(OUT, f))
    print('figures written to', OUT)
