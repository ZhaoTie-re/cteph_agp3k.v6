#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The one quality problem this component actually has, drawn.
#
#           typing_qc.png measures COMPLETENESS — was a call made, how deeply did
#           it resolve — and on this cohort both are near-saturated: call rate
#           takes five distinct values across 3,569 samples and 97 % sit on two of
#           them. Nothing there can see a call that was made confidently and is
#           wrong.
#
#           The frequency check can, and what it found is in
#           docs/OPEN_QUESTIONS.md §3: the controls carry roughly twice the share
#           of calls the 1000 Genomes JPT panel never sees, and that difference is
#           differential by phenotype. This figure is that finding, and it is
#           built to answer the ONE question a reader will ask about it — is this
#           depth? — with the data rather than with prose.
#
#           IT IS NOT DEPTH, and panels (a) and (d) are the argument. Cases and
#           controls sit at the same measured depth; the gap survives matching on
#           it. Depth does have an effect, visible in panel (a) where phenotype is
#           held fixed within the case platforms, but there is no case–control
#           depth difference for that effect to act on.
#
#           An allele the panel never carries is not necessarily a wrong call —
#           105 samples cannot sample a 0.3 % allele. The SHARE is the signal, and
#           only the comparison between groups is interpretable.
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

PLOT_H = 6.0
# The depth window panel (b) matches on. Chosen to be the overlap of the two
# platforms that carry the comparison, not tuned to the answer.
WIN_LO, WIN_HI = 17.0, 21.0


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample', required=True, help='allele_frequency_sample.tsv')
    p.add_argument('--summary', required=True, help='allele_frequency_summary.tsv')
    p.add_argument('--control-group', default='AGP3K')
    p.add_argument('--population', default='JPT')
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


def main():
    args = parse_args()
    S.setup_style('slide')
    d = pd.read_csv(args.sample, sep='\t', dtype={'sample_id': str})
    d['observed_depth'] = pd.to_numeric(d.get('observed_depth'), errors='coerce')
    if 'platform' not in d.columns:
        raise SystemExit(f'ABORT: {args.sample} has no platform column')
    ctrl = d[d.group == args.control_group]
    case = d[d.group != args.control_group]
    if not len(ctrl) or not len(case):
        raise SystemExit(f'ABORT: need both groups; found {sorted(set(d.group))}')

    fig, axes = plt.subplots(2, 2, figsize=(S.COL_DOUBLE, PLOT_H))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()
    C_CTRL, C_CASE = S.DATA, S.ACCENT

    # ── (a) the whole story on one axis: share against MEASURED depth ────────
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
                         textcoords='offset points', xytext=(8, 4), fontsize=6.2,
                         color=S.INK_SOFT) for _, r in pl.iterrows()]
    ax_a.set_xlim(pl.depth.min() - 3, pl.depth.max() + 7)
    ax_a.set_ylim(0, pl.share.max() * 1.30)
    if ann:
        S.spread_labels(ax_a, ann, axis='y')
    ax_a.set_xlabel('measured depth, platform median (×)')
    ax_a.set_ylabel('unconfirmed\nchromosomes')
    ax_a.yaxis.set_major_formatter(lambda v, _: f'{v:.0%}')
    S.panel_tag(ax_a, 'a', title='platform, not depth')
    S.legend_inside(ax_a, ax_a.get_legend_handles_labels()[0], loc='upper right')
    S.despine(ax_a)

    # ── (b) the same comparison with depth held fixed ────────────────────────
    w = d[(d.observed_depth >= WIN_LO) & (d.observed_depth <= WIN_HI)]
    top = (w.groupby('platform').size().sort_values(ascending=False).index.tolist())[:2]
    bars, labs, cols = [], [], []
    for pf in top:
        g = w[w.platform == pf]
        bars.append(share(g))
        labs.append(f'{pf}\nn={len(g):,}  {g.observed_depth.median():.1f}×')
        cols.append(C_CTRL if (g.group == args.control_group).mean() > 0.5 else C_CASE)
    p_matched = fisher(w[w.platform == top[0]], w[w.platform == top[1]]) if len(top) == 2 else np.nan
    ax_b.bar(range(len(bars)), bars, color=cols, width=0.55)
    ax_b.set_xticks(range(len(bars)))
    ax_b.set_xticklabels(labs, fontsize=6.5)
    ax_b.yaxis.set_major_formatter(lambda v, _: f'{v:.0%}')
    if p_matched == p_matched:
        ax_b.set_title(f'Fisher {S.p_tex(p_matched)}', fontsize=7.5, color=S.INK_SOFT, pad=2)
    S.panel_tag(ax_b, 'b', title=f'matched {WIN_LO:.0f}–{WIN_HI:.0f}×')
    S.despine(ax_b)

    # ── (c) per sample, so the reader sees a distribution and not a mean ─────
    hi = int(d.n_unconfirmed.max())
    bins = np.arange(-0.5, hi + 1.5)
    for g, col, lab in ((ctrl, C_CTRL, f'controls (n={len(ctrl):,})'),
                        (case, C_CASE, f'cases (n={len(case):,})')):
        h, _ = np.histogram(g.n_unconfirmed, bins=bins)
        ax_c.step(np.arange(hi + 1), h / len(g), where='mid', color=col, lw=1.6, label=lab)
    ax_c.set_xlabel('unconfirmed chromosomes (of 10)')
    ax_c.set_ylabel('fraction of samples')
    ax_c.set_xticks(range(0, hi + 1, 2))
    S.panel_tag(ax_c, 'c', title='per sample')
    S.legend_inside(ax_c, ax_c.get_legend_handles_labels()[0], loc='upper right')
    S.despine(ax_c)

    # ── (d) the depth distributions the downsampling plan assumed differed ───
    dd = d.dropna(subset=['observed_depth'])
    cc, ca = dd[dd.group == args.control_group], dd[dd.group != args.control_group]
    p_depth = float(stats.mannwhitneyu(ca.observed_depth, cc.observed_depth)[1]) \
        if len(cc) and len(ca) else np.nan
    lo, hiq = np.nanpercentile(dd.observed_depth, [0.5, 99.5])
    dbins = np.linspace(lo, hiq, 46)
    for g, col, lab in ((cc, C_CTRL, 'controls'), (ca, C_CASE, 'cases')):
        ax_d.hist(g.observed_depth.clip(lo, hiq), bins=dbins, density=True,
                  color=col, alpha=0.55, label=f'{lab}  med {g.observed_depth.median():.1f}×')
    ax_d.set_xlabel('measured depth (×)')
    ax_d.set_ylabel('density')
    if p_depth == p_depth:
        ax_d.set_title(f'Mann–Whitney {S.p_tex(p_depth)}', fontsize=7.5, color=S.INK_SOFT, pad=2)
    S.panel_tag(ax_d, 'd', title='depth by group')
    S.legend_inside(ax_d, ax_d.get_legend_handles_labels()[0], loc='upper right')
    S.despine(ax_d)

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
    ])
    tab['fisher_p_matched'] = ''
    tab.loc[tab.index[-1], 'fisher_p_matched'] = f'{p_matched:.3g}' if p_matched == p_matched else ''
    tab['mannwhitney_p_depth'] = ''
    tab.loc[tab.index[1], 'mannwhitney_p_depth'] = f'{p_depth:.3g}' if p_depth == p_depth else ''
    tab.to_csv(args.out_table, sep='\t', index=False)

    S.caption_block(
        fig, plot_h=PLOT_H,
        title=(f'Controls carry {share(ctrl):.1%} of called chromosomes on an allele the '
               f'{args.population} panel never sees, against {share(case):.1%} in cases — '
               f'and the difference is not sequencing depth.'),
        panels=[
            'Unconfirmed share against measured depth, one point per platform, area ∝ √n. '
            'Within the case platforms the share falls with depth; the control platform sits '
            'far above the case platform of the same depth.',
            f'The same two platforms restricted to {WIN_LO:.0f}–{WIN_HI:.0f}×, where their '
            f'median depths are equal.',
            'Per-sample counts, so the spread is visible rather than a mean.',
            'The depth distributions the retired downsampling control assumed differed.',
        ],
        notes=(f'"Unconfirmed" = a 2-field call absent from the 1000 Genomes {args.population} '
               f'panel (105 samples). That is not proof a call is wrong — 105 samples cannot '
               f'sample a rare allele — so only the COMPARISON between groups is interpretable, '
               f'never the absolute level. Depth is Observed_Depth from the workbook; the '
               f'30x/15x platform labels are targets and disagree with it. Five loci: '
               f'A, B, C, DQB1, DRB1. See docs/OPEN_QUESTIONS.md §3.'),
        top_pad=0.34, wspace=0.46, hspace=0.72,
        left='auto', right=0.965, margin_axes=list(axes.ravel()))
    fig.savefig(args.out_png)
    plt.close(fig)

    figure_doc.write_doc(
        args.out_png,
        title='Typing accuracy differs between the groups, and it is not depth',
        question='The controls are typed measurably worse than the cases. Is that sequencing '
                 'depth — which could be controlled for — or something that cannot be?',
        panels=[
            ('a', 'Unconfirmed share against measured depth',
             'One point per sequencing platform, x its median measured depth, y the share of '
             'called chromosomes carrying an allele the reference panel never sees; area is '
             'proportional to √n. Read it twice: along the case platforms alone, phenotype is '
             'held fixed and the share falls as depth rises — depth is real. Then compare the '
             'control platform with the case platform beside it at the same depth.'),
            ('b', 'The comparison at matched depth',
             f'The two platforms carrying the comparison, restricted to samples between '
             f'{WIN_LO:.0f}× and {WIN_HI:.0f}× so their median depths are equal. If depth were '
             f'the cause the bars would be level.'),
            ('c', 'Per sample',
             'How many of a sample\'s ten called chromosomes at the five reference loci sit on '
             'an unconfirmed allele. Shown as a distribution because a mean hides that most '
             'samples have none and a tail has four or more.'),
            ('d', 'Depth by group',
             'The measured-depth distributions of cases and controls with a Mann–Whitney test. '
             'This panel is why the downsampling control was retired: there is no depth '
             'difference to remove.'),
        ],
        interpretation=(
            'Two effects are present and only one of them is a confounder. Depth genuinely '
            'affects typing — panel (a) shows it cleanly within the case platforms, where '
            'phenotype cannot be responsible. But cases and controls are at the same measured '
            'depth (panel d), so that effect has nothing to act on between the groups, and the '
            'gap survives matching on depth (panel b). What is left is platform and cohort: '
            'every control is an AGP3K sample on HiSeqX, every case is a PH sample on something '
            'else. Those are perfectly confounded with each other and with phenotype, so which '
            'of read length, library chemistry, alignment reference or batch is responsible '
            'cannot be determined from these data. A common allele depleted more in controls '
            'than in cases reads as enrichment in cases — a false risk association from the '
            'typing alone.'),
        limits=[
            'The reference panel is 105 samples. An allele it never carries may still be real '
            'and rare, so the absolute level of the unconfirmed share is not an error rate. '
            'Only the difference between groups, measured on the same panel, is interpretable.',
            'Five loci only — A, B, C, DQB1, DRB1. DPA1, DPB1, DQA1 and DRB3/4/5 have no '
            'external reference and are invisible here.',
            'The figure identifies the confound; it does not correct it. No within-data control '
            'exists, which is why OPEN_QUESTIONS.md §2 — typing samples of known HLA type — is '
            'now the only remaining option.',
        ],
        methods_ref='../../docs/METHODS.md')

    print(f'[plot_typing_confound] controls {share(ctrl):.2%} vs cases {share(case):.2%}; '
          f'matched-depth Fisher P = {p_matched:.3g}; depth Mann-Whitney P = {p_depth:.3g} '
          f'-> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in plot_typing_confound: {e}', file=sys.stderr)
        sys.exit(1)
