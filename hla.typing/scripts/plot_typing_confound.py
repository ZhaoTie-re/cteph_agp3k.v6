#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The one quality problem this component actually has, drawn.
#
#           typing_qc.png measures COMPLETENESS -- was a call made, how deeply
#           did it resolve -- and on this cohort both are near-saturated. Nothing
#           there can see a call that was made confidently and is wrong.
#
#           The frequency check can, and what it found is in
#           docs/OPEN_QUESTIONS.md section 3: the controls carry a materially
#           larger share of calls the reference panel never sees, and that
#           difference tracks phenotype. This figure is that finding, built to
#           answer the two questions a reader will ask about it -- is it depth,
#           and is it an artefact of the reference panel -- with data rather than
#           with prose.
#
#           IT IS NOT DEPTH: cases and controls sit at the same measured depth
#           and the gap survives matching on it. Depth does have an effect,
#           visible in panel (a) where phenotype is held fixed within the case
#           platforms, but there is no case-control depth difference for that
#           effect to act on.
#
#           IT IS NOT THE PANEL EITHER, and panel (d) is why that is now shown
#           rather than asserted. An allele the panel never carries is not
#           necessarily a wrong call: a panel samples a finite number of
#           chromosomes and a rare allele can be missing from it. So the ABSOLUTE
#           level is a property of the panel. Running the same calls against a
#           second, far smaller panel moves that level by a third and leaves the
#           RATIO between the groups where it was. The contrast is the readable
#           quantity; the level is not.
#
#           NO NUMBER IS WRITTEN INTO THIS FILE. Every value in the caption is an
#           f-string field read from the tables, and each panel is named by the
#           caller, because the last time a reference panel changed the captions
#           kept describing the old one for three months.
# Component: hla.typing
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402
from scipy import stats                  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'analysis' / '_shared' / 'scripts'))
import plot_style as S                   # noqa: E402
import figure_doc                        # noqa: E402

PLOT_H = 5.6
# The depth window panel (b) matches on. Chosen as the overlap of the two
# platforms that carry the comparison, not tuned to the answer.
WIN_LO, WIN_HI = 17.0, 21.0


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample', required=True,
                   help='allele_frequency_sample.tsv from the PRIMARY reference panel')
    p.add_argument('--reference-name', default='jMorp 61KJPN-HLA',
                   help='how the primary panel is named on the figure')
    p.add_argument('--sample-secondary', required=True,
                   help='allele_frequency_sample.tsv from the SECOND panel. Panel (d) is the '
                        'whole reason both panels are run: it separates what is a property of '
                        'the reference from what is a property of the typing.')
    p.add_argument('--reference-name-secondary', default='1000 Genomes JPT')
    p.add_argument('--control-group', default='AGP3K')
    p.add_argument('--out-png', required=True)
    p.add_argument('--out-table', required=True)
    return p.parse_args()


def share(d):
    """Chromosomes on an unconfirmed allele, over chromosomes called."""
    n = d.n_chr.sum()
    return (d.n_unconfirmed.sum() / n) if n else np.nan


def fisher(a, b):
    ta = [[a.n_unconfirmed.sum(), a.n_chr.sum() - a.n_unconfirmed.sum()],
          [b.n_unconfirmed.sum(), b.n_chr.sum() - b.n_unconfirmed.sum()]]
    if min(min(r) for r in ta) < 0 or not a.n_chr.sum() or not b.n_chr.sum():
        return np.nan
    return float(stats.fisher_exact(ta)[1])


def load(path, control_group):
    d = pd.read_csv(path, sep='\t', dtype={'sample_id': str})
    for col in ('n_chr', 'n_unconfirmed', 'group'):
        if col not in d.columns:
            raise SystemExit(f'ABORT: {path} has no `{col}` column')
    d['observed_depth'] = pd.to_numeric(d.get('observed_depth'), errors='coerce')
    ctrl = d[d.group == control_group]
    case = d[d.group != control_group]
    if not len(ctrl) or not len(case):
        raise SystemExit(f'ABORT: need both groups in {path}; found {sorted(set(d.group))}')
    return d, ctrl, case


def stat_note(ax, text, *, y=0.985, va='top'):
    """A computed statistic INSIDE the axes.

    It used to be written with ax.set_title(), which panel_tag() then overwrote
    with the panel's short title -- both write loc='center' -- so every P value
    this script computed was silently absent from the figure it was computed for.
    """
    ax.text(0.5, y, text, transform=ax.transAxes, ha='center', va=va,
            fontsize=6.8, color=S.INK_SOFT)


def main():
    args = parse_args()
    S.setup_style('slide')
    d, ctrl, case = load(args.sample, args.control_group)
    d2, ctrl2, case2 = load(args.sample_secondary, args.control_group)
    if 'platform' not in d.columns:
        raise SystemExit(f'ABORT: {args.sample} has no platform column')

    fig, axes = plt.subplots(2, 2, figsize=(S.COL_DOUBLE, PLOT_H))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()
    C_CTRL, C_CASE = S.DATA, S.ACCENT

    # -- (a) the whole story on one axis: share against MEASURED depth --------
    # One point per platform, so the two effects are separable by eye: within the
    # case platforms the points fall with depth, and the control platform sits far
    # above the case platform it shares a depth with.
    rows = []
    for pf, g in d.groupby('platform'):
        is_ctrl = (g.group == args.control_group).mean() > 0.5
        rows.append({'platform': pf, 'n': len(g), 'depth': g.observed_depth.median(),
                     'share': share(g), 'is_ctrl': is_ctrl})
    pl = pd.DataFrame(rows).sort_values('depth')
    for is_ctrl, col, lab in ((True, C_CTRL, 'controls'), (False, C_CASE, 'cases')):
        s = pl[pl.is_ctrl == is_ctrl]
        ax_a.scatter(s.depth, s.share, s=np.sqrt(s.n) * 9, color=col, alpha=0.85,
                     edgecolors='white', linewidths=0.8, zorder=3, label=lab)
    # The vendor prefix is dropped HERE ONLY, and only from the point labels: five
    # full names at 6 pt on a 3.5 in axis overlap into an unreadable band. The
    # table and every other panel keep the name in full.
    short = {p: str(p).split('-')[-1] if '-' in str(p) else str(p) for p in pl.platform}
    ann = [ax_a.annotate(short[r.platform], (r.depth, r.share),
                         textcoords='offset points',
                         xytext=(np.sqrt(np.sqrt(r.n) * 9 / np.pi) + 4.0, 3), fontsize=6.2,
                         color=S.INK_SOFT) for _, r in pl.iterrows()]
    ax_a.set_xlim(pl.depth.min() - 3, pl.depth.max() + 7)
    ax_a.set_ylim(-pl.share.max() * 0.10, pl.share.max() * 1.42)
    if ann:
        S.spread_labels(ax_a, ann, axis='y')
    ax_a.set_xlabel('measured depth, platform median (×)')
    ax_a.set_ylabel('unconfirmed\nchromosomes')
    ax_a.yaxis.set_major_formatter(lambda v, _: f'{v:.0%}')

    dd = d.dropna(subset=['observed_depth'])
    cc, ca = dd[dd.group == args.control_group], dd[dd.group != args.control_group]
    p_depth = float(stats.mannwhitneyu(ca.observed_depth, cc.observed_depth)[1]) \
        if len(cc) and len(ca) else np.nan
    S.panel_tag(ax_a, 'a', title='platform, not depth')
    S.legend_inside(ax_a, ax_a.get_legend_handles_labels()[0], loc='upper right')
    if p_depth == p_depth:
        stat_note(ax_a, f'group depths {cc.observed_depth.median():.1f}× vs '
                        f'{ca.observed_depth.median():.1f}×, '
                        f'Mann–Whitney {S.p_tex(p_depth)}', y=0.015, va='bottom')
    S.despine(ax_a)

    # -- (b) the same comparison with depth held fixed -----------------------
    w = d[(d.observed_depth >= WIN_LO) & (d.observed_depth <= WIN_HI)]
    top = (w.groupby('platform').size().sort_values(ascending=False).index.tolist())[:2]
    bars, labs, cols = [], [], []
    for pf in top:
        g = w[w.platform == pf]
        bars.append(share(g))
        labs.append(f'{pf}\nn={len(g):,}  {g.observed_depth.median():.1f}×')
        cols.append(C_CTRL if (g.group == args.control_group).mean() > 0.5 else C_CASE)
    p_matched = fisher(w[w.platform == top[0]], w[w.platform == top[1]]) if len(top) == 2 else np.nan
    if len(bars) < 2:
        # An empty or one-sided window is a real outcome once the cohort is
        # restricted, and it has to be legible AS one rather than as a blank axis.
        ax_b.text(0.5, 0.5, f'no two platforms overlap\nat {WIN_LO:.0f}–{WIN_HI:.0f}×',
                  transform=ax_b.transAxes, ha='center', va='center',
                  fontsize=7.5, color=S.INK_SOFT)
        ax_b.set_xticks([])
    else:
        ax_b.bar(range(len(bars)), bars, color=cols, width=0.55)
        ax_b.set_xticks(range(len(bars)))
        ax_b.set_xticklabels(labs, fontsize=6.5)
        ax_b.set_ylim(0, max(bars) * 1.34)
        if p_matched == p_matched:
            stat_note(ax_b, f'Fisher {S.p_tex(p_matched)}')
    ax_b.yaxis.set_major_formatter(lambda v, _: f'{v:.0%}')
    ax_b.set_ylabel('unconfirmed\nchromosomes')
    S.panel_tag(ax_b, 'b', title=f'matched {WIN_LO:.0f}–{WIN_HI:.0f}×')
    S.despine(ax_b)

    # -- (c) per sample, so the reader sees a distribution and not a mean ----
    hi = int(d.n_unconfirmed.max())
    bins = np.arange(-0.5, hi + 1.5)
    for g, col, lab in ((ctrl, C_CTRL, f'controls (n={len(ctrl):,})'),
                        (case, C_CASE, f'cases (n={len(case):,})')):
        h, _ = np.histogram(g.n_unconfirmed, bins=bins)
        ax_c.step(np.arange(hi + 1), h / len(g), where='mid', color=col, lw=1.6, label=lab)
    n_lo, n_hi = int(d.n_chr.min()), int(d.n_chr.max())
    denom = f'{n_hi}' if n_lo == n_hi else f'{n_lo}–{n_hi}'
    ax_c.set_xlabel(f'unconfirmed chromosomes (of {denom})')
    ax_c.set_ylabel('fraction of samples')
    ax_c.set_xticks(range(0, hi + 1, 2))
    S.panel_tag(ax_c, 'c', title='per sample')
    S.legend_inside(ax_c, ax_c.get_legend_handles_labels()[0], loc='upper right')
    S.despine(ax_c)

    # -- (d) the same contrast under both reference panels -------------------
    # The point of running two panels, made visible. The level is a property of
    # the panel; the ratio between the groups is not.
    panels = [(args.reference_name, ctrl, case),
              (args.reference_name_secondary, ctrl2, case2)]
    x = np.arange(len(panels))
    wbar = 0.34
    sh_c = [share(c) for _, c, _ in panels]
    sh_k = [share(k) for _, _, k in panels]
    ax_d.bar(x - wbar / 2, sh_c, wbar, color=C_CTRL)
    ax_d.bar(x + wbar / 2, sh_k, wbar, color=C_CASE)
    for i, (c, k) in enumerate(zip(sh_c, sh_k)):
        ax_d.annotate(f'{c / k:.2f}×', (i, max(c, k)), textcoords='offset points',
                      xytext=(0, 7), ha='center', fontsize=7.0, color=S.INK)
    ax_d.set_xticks(x)
    ax_d.set_xticklabels([f'{nm}\n{int(dd_.n_chr.sum()):,} chr'
                          for (nm, c_, k_), dd_ in zip(panels, (d, d2))], fontsize=6.5)
    ax_d.set_ylim(0, max(sh_c + sh_k) * 1.26)
    ax_d.yaxis.set_major_formatter(lambda v, _: f'{v:.0%}')
    ax_d.set_ylabel('unconfirmed\nchromosomes')
    S.panel_tag(ax_d, 'd', title='both panels')
    S.despine(ax_d)

    # -- the table, which carries everything the four panels cannot ----------
    tab = pd.DataFrame([
        {'stratum': 'controls', 'n': len(ctrl), 'median_depth': round(ctrl.observed_depth.median(), 2),
         'n_chr': int(ctrl.n_chr.sum()), 'n_unconfirmed': int(ctrl.n_unconfirmed.sum()),
         'share_unconfirmed': round(share(ctrl), 5)},
        {'stratum': 'cases', 'n': len(case), 'median_depth': round(case.observed_depth.median(), 2),
         'n_chr': int(case.n_chr.sum()), 'n_unconfirmed': int(case.n_unconfirmed.sum()),
         'share_unconfirmed': round(share(case), 5)},
    ] + [
        {'stratum': f'platform:{r.platform}', 'n': int(r.n), 'median_depth': round(r.depth, 2),
         'n_chr': int(d[d.platform == r.platform].n_chr.sum()),
         'n_unconfirmed': int(d[d.platform == r.platform].n_unconfirmed.sum()),
         'share_unconfirmed': round(r.share, 5)} for _, r in pl.iterrows()
    ] + [
        {'stratum': f'matched {WIN_LO:.0f}-{WIN_HI:.0f}x:{pf}', 'n': int((w.platform == pf).sum()),
         'median_depth': round(w[w.platform == pf].observed_depth.median(), 2),
         'n_chr': int(w[w.platform == pf].n_chr.sum()),
         'n_unconfirmed': int(w[w.platform == pf].n_unconfirmed.sum()),
         'share_unconfirmed': round(share(w[w.platform == pf]), 5)} for pf in top
    ] + [
        {'stratum': f'{args.reference_name_secondary}:{lab}', 'n': len(g),
         'median_depth': round(g.observed_depth.median(), 2),
         'n_chr': int(g.n_chr.sum()), 'n_unconfirmed': int(g.n_unconfirmed.sum()),
         'share_unconfirmed': round(share(g), 5)}
        for lab, g in (('controls', ctrl2), ('cases', case2))
    ])
    p_2nd = fisher(ctrl2, case2)
    tab['fisher_p_matched'] = ''
    tab['mannwhitney_p_depth'] = ''
    tab['fisher_p_group'] = ''
    if len(top) == 2 and p_matched == p_matched:
        tab.loc[tab.index[-3], 'fisher_p_matched'] = f'{p_matched:.3g}'
    if p_depth == p_depth:
        tab.loc[tab.index[1], 'mannwhitney_p_depth'] = f'{p_depth:.3g}'
    p_1st = fisher(ctrl, case)
    if p_1st == p_1st:
        tab.loc[tab.index[1], 'fisher_p_group'] = f'{p_1st:.3g}'
    if p_2nd == p_2nd:
        tab.loc[tab.index[-1], 'fisher_p_group'] = f'{p_2nd:.3g}'
    tab.to_csv(args.out_table, sep='\t', index=False)

    n_loci = int(round(n_hi / 2))
    S.caption_block(
        fig, plot_h=PLOT_H,
        title=(f'Controls carry {share(ctrl):.1%} of called chromosomes on an allele the '
               f'reference never sees, cases {share(case):.1%} — not depth, not the panel.'),
        panels=[
            'Share against measured depth, one point per platform, area ∝ √n.',
            f'The two platforms of the comparison at {WIN_LO:.0f}–{WIN_HI:.0f}×, equal medians.',
            'Per sample, so the spread shows rather than a mean.',
            'Both panels, control/case ratio above each pair.',
        ],
        notes=(f'"Unconfirmed" = a call the named panel does not carry. {args.reference_name}, '
               f'{n_loci} loci, {int(d.n_chr.sum()):,} chr; {args.reference_name_secondary}, '
               f'{int(round(int(d2.n_chr.max()) / 2))} loci, {int(d2.n_chr.sum()):,} chr. Only '
               f'the group contrast is readable, never the level. docs/OPEN_QUESTIONS.md §3.'),
        top_pad=0.32, wspace=0.46, hspace=0.60,
        left='auto', right=0.965, margin_axes=list(axes.ravel()))
    fig.savefig(args.out_png)
    plt.close(fig)

    figure_doc.write_doc(
        args.out_png,
        title='Typing accuracy differs between the groups; it is not depth and not the panel',
        question='The controls are typed measurably worse than the cases. Is that sequencing '
                 'depth — which could be controlled for — or an artefact of the reference '
                 'panel, or something that can be neither corrected nor explained away?',
        panels=[
            ('a', 'Unconfirmed share against measured depth',
             'One point per sequencing platform, x its median measured depth, y the share of '
             'called chromosomes carrying an allele the reference panel never sees; area is '
             'proportional to √n. Read it twice: along the case platforms alone, phenotype is '
             'held fixed and the share falls as depth rises — depth is real. Then compare the '
             'control platform with the case platform beside it at the same depth. The note in '
             'the panel gives the two group medians and a Mann–Whitney test of their equality.'),
            ('b', 'The comparison at matched depth',
             f'The two platforms carrying the comparison, restricted to samples between '
             f'{WIN_LO:.0f}× and {WIN_HI:.0f}× so their median depths are equal. If depth were '
             f'the cause the bars would be level. Fisher exact on the 2×2 of unconfirmed against '
             f'confirmed chromosomes.'),
            ('c', 'Per sample',
             'How many of a sample\'s called chromosomes at the reference loci sit on an '
             'unconfirmed allele. Shown as a distribution because a mean hides that most '
             'samples have none and a tail has several. The denominator is printed on the axis '
             'and is read from the data, not assumed.'),
            ('d', 'Both reference panels',
             'The control and case shares computed twice from the same calls, once against each '
             'panel, with the control/case ratio printed above each pair. This is the panel that '
             'licenses the reading of every other one: if the ratio moved with the panel, the '
             'finding would be about the reference rather than about the typing.'),
        ],
        numbers={
            f'controls, {args.reference_name}': f'{share(ctrl):.2%}',
            f'cases, {args.reference_name}': f'{share(case):.2%}',
            'ratio, primary panel': f'{share(ctrl) / share(case):.2f}',
            f'controls, {args.reference_name_secondary}': f'{share(ctrl2):.2%}',
            f'cases, {args.reference_name_secondary}': f'{share(case2):.2%}',
            'ratio, secondary panel': f'{share(ctrl2) / share(case2):.2f}',
            'Fisher P, primary panel': f'{p_1st:.3g}',
            'Fisher P, secondary panel': f'{p_2nd:.3g}',
            'control median depth': f'{cc.observed_depth.median():.2f}×',
            'case median depth': f'{ca.observed_depth.median():.2f}×',
            'Mann-Whitney P, depth': f'{p_depth:.3g}',
            f'Fisher P, matched {WIN_LO:.0f}-{WIN_HI:.0f}x': (f'{p_matched:.3g}'
                                                              if p_matched == p_matched else 'n/a'),
        },
        tables=[('Every stratum behind the four panels, as written to typing_confound.tsv', tab)],
        interpretation=(
            'Three explanations are available and two of them are ruled out here. Depth '
            'genuinely affects typing — panel (a) shows it cleanly within the case platforms, '
            'where phenotype cannot be responsible — but cases and controls sit at the same '
            'measured depth and the gap survives matching on it, so that effect has nothing to '
            f'act on between the groups. The reference panel genuinely sets the level — panel '
            f'(d) shows the control level move from {share(ctrl):.2%} against '
            f'{args.reference_name} to {share(ctrl2):.2%} against '
            f'{args.reference_name_secondary} — but the ratio '
            f'between the groups is stable across both '
            f'({share(ctrl) / share(case):.2f} and {share(ctrl2) / share(case2):.2f}), so '
            f'the contrast is not an artefact of what the panel happens to carry. What is left '
            'is platform and cohort: every control is an AGP3K sample on HiSeqX, every case is '
            'a PH sample on something else. Those are perfectly confounded with each other and '
            'with phenotype, so which of read length, library chemistry, alignment reference or '
            'batch is responsible cannot be determined from these data. A common allele '
            'depleted more in controls than in cases reads as enrichment in cases — a false '
            'risk association from the typing alone.'),
        limits=[
            'The absolute unconfirmed share is not an error rate. An allele a panel never '
            'carries may still be real and rare; panel (d) measures how much that matters by '
            'changing the panel. Only the difference between groups, measured on the same '
            'panel, is interpretable.',
            'Only the loci a panel carries are visible. Loci absent from a panel are not '
            'scored against it at all, which is why the two panels contribute different '
            'numbers of chromosomes per sample and why each locus list travels with its panel.',
            'The cohort is ancestry-restricted, so an unconfirmed allele can no longer be '
            'explained by the sample not being mainland Japanese. That removes an explanation; '
            'it does not remove the platform confound, which is what remains.',
            'The figure identifies the confound; it does not correct it. No within-data control '
            'exists, which is why OPEN_QUESTIONS.md §2 — typing samples of known HLA type — is '
            'now the only remaining option.',
        ],
        methods_ref='../../docs/METHODS.md')

    print(f'[plot_typing_confound] {args.reference_name}: controls {share(ctrl):.2%} vs cases '
          f'{share(case):.2%} ({share(ctrl) / share(case):.2f}x, Fisher P = {p_1st:.3g})')
    print(f'[plot_typing_confound] {args.reference_name_secondary}: controls {share(ctrl2):.2%} '
          f'vs cases {share(case2):.2%} ({share(ctrl2) / share(case2):.2f}x, '
          f'Fisher P = {p_2nd:.3g})')
    print(f'[plot_typing_confound] matched-depth Fisher P = {p_matched:.3g}; '
          f'depth Mann-Whitney P = {p_depth:.3g} -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:                       # noqa: BLE001
        print(f'Error in plot_typing_confound: {e}', file=sys.stderr)
        sys.exit(1)
