#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : One figure answering the question this component cannot answer in
#           prose: does typing quality differ by sequencing platform?
#
#           It has to be read by platform, not by case/control. Cases and controls
#           were sequenced differently and with no overlap, so any difference
#           between the groups is confounded by construction.
#
#           EVERY METRIC HERE IS ABOUT COMPLETENESS, and on this cohort completeness
#           is nearly saturated: call rate is k/19 and takes a handful of distinct
#           values, nearly all samples sitting on two of them. An earlier version
#           plotted one jittered point per sample to show those few values, and drew
#           nineteen per-gene bars of which fourteen were exactly 1.0. Both are now
#           compositions, which is the encoding a saturated proportion needs.
#
#           PANELS (c) AND (d) SHARE ONE GENE AXIS. They are two measurements of the
#           same nineteen genes, and they used to be sorted differently, so the gene
#           names were printed twice in two different orders and the reader had to
#           re-find each gene to compare. One order, one set of labels.
#
#           NO CAPTION IS RENDERED INTO THIS PNG. Every word of explanation lives
#           in report/hla_typing_report.qmd, written for a reader who has never
#           heard of HLA. The caption block was taking 34 % of this canvas.
#
#           NO NUMBER IS WRITTEN INTO THIS FILE either -- write_doc() still emits
#           the sidecar, and every value in it is read from the tables.
#
#           A LOW `call_rate` AT DRB3/4/5 IS NOT A FAILURE. Their shortfall is
#           entirely `not_typed` — the gene is absent from that haplotype — and
#           `failed` is 0 everywhere. Plotting it as a call rate showed correct
#           biology as if it were breakage, which is why panel (c) is a four-way
#           composition instead.
#
#           Whether a call is RIGHT is not measurable here at all. That is
#           typing_confound.png, off the external frequency check.
# Component: hla.typing
# ---------------------------------------------------------------------------
import argparse
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'analysis' / '_shared' / 'scripts'))
import plot_style as S                   # noqa: E402
import figure_doc                        # noqa: E402

# 7.6, not 5.4: the two bar panels carry one row per gene, and 19 rows do not
# fit a 2.2 in axis — the labels overlapped into an unreadable band. The row
# count is data-dependent, so the height is set from it below.
PLOT_H_BASE = 2.75    # the two composition panels (5 platform rows each) plus the
                      # space their below-axis legends and x labels need
ROW_IN = 0.200        # inches per gene row = 14.4 pt of pitch. It was 0.155, and
                      # the tick font was left at the 12 pt 'slide' default, so 19
                      # labels in an 8.8 pt pitch fused into three solid runs of
                      # ink. A row's font is now set explicitly below; the two
                      # numbers have to be chosen together or this recurs.
FS_ROW = 7.5          # gene row labels
FS_SMALL = 7.0        # legends, in-bar values, platform names


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample-qc', required=True)
    p.add_argument('--gene-qc', required=True)
    p.add_argument('--sites', required=True)
    p.add_argument('--out-png', required=True)
    return p.parse_args()


def panel_head(ax, letter, title, *, pad=6.0, size=None):
    """Bold panel letter and a short title, as ONE left-aligned artist.

    plot_style.panel_tag() writes the letter at loc='left' and the title at
    loc='center', and its own docstring says not to pass a title below about
    2.2 in wide because the two collide. Most panels in this component are
    narrower than that. Written as one left-aligned string they cannot collide,
    and it reads as a label rather than as two competing ones.

    The title says WHAT the panel shows. It never says what to conclude -- that
    is the document's job, which is why no caption is rendered into these PNGs.
    """
    ax.set_title(rf'$\bf{{{letter}}}$  {title}', loc='left', pad=pad,
                 fontsize=size or (plt.rcParams['axes.titlesize'] - 2))


TITLE_H = 0.34        # inches reserved for the figure title


def lay_out(fig, axes, *, plot_h, top_pad, wspace, hspace, right=0.975, title=None):
    """Margins, canvas and figure title, with no caption block.

    caption_block() is the only thing in the shared toolkit that sets the canvas
    height and the measured left margin, so dropping it means doing that here.
    A transcription of what it does, minus the caption text -- but WITH a title:
    a figure lifted out of the report still has to say what it is, and each panel
    says what it shows. What none of them say is what to conclude; that is the
    document's job, and it is why the caption block is gone.
    """
    fig.set_layout_engine('none')
    xlab_in = (plt.rcParams['axes.labelsize']
               + plt.rcParams['xtick.labelsize'] + 8) * 1.5 / 72.0
    head = TITLE_H if title else 0.0
    fh = plot_h + head                       # the panels keep the size they were given
    fig.set_size_inches(S.COL_DOUBLE, fh, forward=True)
    kw = dict(right=right, top=1.0 - (top_pad + head) / fh,
              bottom=xlab_in / fh, wspace=wspace, hspace=hspace)
    left = S.fit_left_margin(fig, axes)
    fig.subplots_adjust(left=left, **kw)
    need = S.fit_left_margin(fig, axes)
    if need > left + 0.002:
        fig.subplots_adjust(left=need, **kw)
    if title:
        # WRAPPED. An unwrapped suptitle silently runs off both edges of the
        # canvas -- matplotlib never clips it and never warns. The character
        # estimate is deliberately conservative; a title that wraps to two lines
        # is fine, a title bleeding into the margin is not.
        fs = plt.rcParams['axes.titlesize'] + 1.5
        per_line = max(20, int(S.COL_DOUBLE * 72.0 * 0.88 / (fs * 0.56)))
        wrapped = textwrap.fill(title, per_line)
        n_lines = wrapped.count('\n') + 1
        if n_lines > 1:                       # give the extra line its own room
            fh += (n_lines - 1) * fs * 1.25 / 72.0
            fig.set_size_inches(S.COL_DOUBLE, fh, forward=True)
            fig.subplots_adjust(top=1.0 - (top_pad + head
                                           + (n_lines - 1) * fs * 1.25 / 72.0) / fh,
                                bottom=xlab_in / fh)
        fig.suptitle(wrapped, y=1.0 - 0.09 / fh, va='top', ha='center',
                     fontsize=fs, fontweight='bold', color=S.INK)


def order_platforms(df):
    """Platforms by sample count, largest first — the controls dominate and
    should anchor the eye."""
    return df.groupby('platform').size().sort_values(ascending=False).index.tolist()


def short_platform(name):
    """Drop the vendor prefix, HERE ONLY, and only from the row labels.

    `DNBSeq-G400RS 30x` at 12 pt is 1.86 in wide, and fit_left_margin sizes ONE
    left margin for the whole figure from the widest label in it -- which left a
    1.4 x 3.3 in blank column beside the gene panels below. The full names are in
    the sidecar and in every table.
    """
    t = str(name)
    return t.split('-', 1)[1] if '-' in t else t


def main():
    args = parse_args()
    S.setup_style('slide')

    sq = pd.read_csv(args.sample_qc, sep='\t')
    gq = pd.read_csv(args.gene_qc, sep='\t')
    sites = pd.read_csv(args.sites, sep='\t')
    if 'platform' not in sq.columns:
        raise SystemExit(f'ABORT: {args.sample_qc} has no platform column')

    plats = order_platforms(sq)
    n_rows = max(len(gq), int((sites.polymorphic > 0).sum()), 1)
    bar_h = max(1.6, n_rows * ROW_IN)
    plot_h = PLOT_H_BASE + bar_h
    fig, axes = plt.subplots(2, 2, figsize=(S.COL_DOUBLE, plot_h),
                             gridspec_kw={'height_ratios': [PLOT_H_BASE - 1.4, bar_h]})
    ax_call, ax_depth, ax_gene, ax_amb = axes.ravel()

    # (a) FIELD RESOLUTION as a composition, not a mean. The mean lives in
    # 2.84-2.95 for every platform, a range no reader can act on; the 2/3/4-field
    # split is what actually differs and it is a proportion, so it is stacked.
    fld = [c for c in ('field2', 'field3', 'field4') if c in sq.columns]
    lab_f = {'field2': '2-field', 'field3': '3-field', 'field4': '4-field'}
    short = [short_platform(q) for q in plats]
    comp = sq.groupby('platform')[fld].sum().reindex(plats)
    comp = comp.div(comp.sum(axis=1), axis=0)
    # A legend swatch for a series that is identically zero is a lie of omission.
    # field4 is 0 for every platform here -- HLA-HD never resolves past 3 fields on
    # this cohort -- so the series is dropped and the fact is stated in the caption.
    fld = [c for c in fld if comp[c].sum() > 0]
    left = np.zeros(len(plats))
    for i, c in enumerate(fld):
        ax_call.barh(range(len(plats)), comp[c], left=left, height=0.66,
                     color=S.COHORT_RAMP[i % len(S.COHORT_RAMP)], label=lab_f[c])
        left += comp[c].to_numpy()
    # The discriminating number is the 2-field share; print it so the platform gap
    # is read off the figure rather than estimated from a 4-point-wide segment.
    # WHITE, and inside the next segment. Drawn in ink at v + 0.015 it landed on
    # the mid-blue 3-field bar with 58-70 % of each glyph box over the fill, which
    # is unreadable on a 100 %-stacked bar -- there is no white space to land in.
    for y, v in enumerate(comp[fld[0]].to_numpy()):
        ax_call.text(v + 0.012, y, f'{v:.1%}', va='center', ha='left',
                     fontsize=FS_SMALL, color='white', fontweight='bold')
    ax_call.set_yticks(range(len(plats)))
    ax_call.set_yticklabels(short, fontsize=FS_ROW + 1.5)
    ax_call.set_ylim(len(plats) - 0.5, -0.5)
    ax_call.set_xlim(0, 1)
    ax_call.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax_call.tick_params(axis='x', labelsize=FS_ROW + 1.5)
    ax_call.xaxis.set_major_formatter(lambda v, _: f'{v:.0%}')
    ax_call.set_xlabel('calls')
    panel_head(ax_call, 'a', 'field resolution', pad=18.0)
    ax_call.legend(loc='lower center', bbox_to_anchor=(0.5, 1.005), ncol=len(fld),
                   frameon=False, handlelength=1.0, columnspacing=0.9,
                   fontsize=FS_SMALL)
    S.despine(ax_call)

    # (b) AMBIGUITY, which varies and was never plotted. 2,957 records across the
    # cohort; per sample it runs 0-5, and HLA-HD emitting several candidate pairs
    # is the closest thing to a confidence signal the tool gives.
    amb = sq.groupby('platform')['ambiguous'].apply(
        lambda s: pd.Series({f'{k}': (s == k).mean() for k in range(0, 4)}
                            | {'4+': (s >= 4).mean()})).unstack().reindex(plats)
    left = np.zeros(len(plats))
    ramp = [S.NEUTRAL, S.NEUTRAL_D, S.DATA, S.DATA_DARK, S.ACCENT]
    for i, c in enumerate(['0', '1', '2', '3', '4+']):
        # A swatch for a category that is zero on EVERY platform is the same lie
        # of omission panel (a) guards against. A category that is zero on some
        # platforms and not others stays: that is data, not absence.
        if c not in amb.columns or not amb[c].sum():
            continue
        ax_depth.barh(range(len(plats)), amb[c], left=left, height=0.66,
                      color=ramp[i], label=c)
        left += amb[c].to_numpy()
    ax_depth.set_yticks(range(len(plats)))
    ax_depth.set_yticklabels([])
    ax_depth.tick_params(axis='y', length=0)
    ax_depth.set_ylim(len(plats) - 0.5, -0.5)
    ax_depth.set_xlim(0, 1)
    ax_depth.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax_depth.tick_params(axis='x', labelsize=FS_ROW + 1.5)
    ax_depth.xaxis.set_major_formatter(lambda v, _: f'{v:.0%}')
    ax_depth.set_xlabel('samples')
    panel_head(ax_depth, 'b', 'ambiguity', pad=18.0)
    ax_depth.legend(loc='lower center', bbox_to_anchor=(0.5, 1.005), ncol=5,
                    frameon=False, handlelength=0.9, columnspacing=0.8,
                    fontsize=FS_SMALL)
    S.despine(ax_depth)

    # (c) THE FOUR OUTCOMES per gene, not a call rate. `called` and `hemizygous`
    # are both successes and their ratio is real biology — DMA is hemizygous in
    # 73 % of samples and HLA-A in 14 %, because DMA is nearly monomorphic here.
    # `not_typed` is the gene being absent from the haplotype, which is what makes
    # DRB3/4/5 look broken on a call-rate axis when nothing failed.
    # `failed` is dropped when it is zero for every gene -- the caption says there
    # were none, and a red swatch with no red pixels in the panel invites a reader
    # to look for something that is not there. Panel (a) has applied this guard
    # since it was written; this panel did not.
    states = [c for c in ('called', 'hemizygous', 'not_typed', 'failed')
              if c in gq.columns and gq[c].sum()]
    gq = gq.copy()
    tot = gq[states].sum(axis=1)
    gq = gq.assign(**{c: gq[c] / tot for c in states})
    # Sorted ASCENDING but drawn on a NON-inverted axis -- unlike (a) and (b), which
    # invert with set_ylim -- so row 0 lands at the BOTTOM and the panel reads
    # DESCENDING top to bottom: most not_typed first, then most hemizygous. Switching
    # this to ascending=False would reverse the FIGURE, not the sort.
    gq = gq.sort_values(['not_typed', 'hemizygous'])
    cols_s = {'called': S.DATA_DARK, 'hemizygous': S.DATA,
              'not_typed': S.NEUTRAL, 'failed': S.ACCENT}
    left = np.zeros(len(gq))
    for c in states:
        ax_gene.barh(range(len(gq)), gq[c], left=left, height=0.72,
                     color=cols_s[c], label=c)
        left += gq[c].to_numpy()
    ax_gene.set_yticks(range(len(gq)))
    ax_gene.set_yticklabels(gq['gene'], fontstyle='italic', fontsize=FS_ROW)
    ax_gene.set_xlim(0, 1)
    ax_gene.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax_gene.tick_params(axis='x', labelsize=FS_ROW + 1.5)
    ax_gene.xaxis.set_major_formatter(lambda v, _: f'{v:.0%}')
    ax_gene.set_xlabel('samples')
    panel_head(ax_gene, 'c', 'outcome per gene', pad=18.0)
    ax_gene.legend(loc='lower center', bbox_to_anchor=(0.5, 1.005), ncol=len(states),
                   frameon=False, handlelength=0.9, columnspacing=0.8,
                   fontsize=FS_SMALL)
    S.despine(ax_gene)

    # (d) polymorphic residue positions per gene, ON (c)'S ROW ORDER.
    # These are two measurements of the same nineteen genes. Sorting each panel by
    # its own value printed the gene names twice in two different orders, so
    # comparing a gene's outcome with its testable-site count meant finding it
    # again in the second list. One order, and only (c) carries the labels.
    sp = sites.set_index('gene').reindex(gq['gene'])['polymorphic'].fillna(0)
    ax_amb.barh(range(len(gq)), sp.to_numpy(), color=S.DATA_DARK, height=0.72)
    ax_amb.set_yticks(range(len(gq)))
    ax_amb.set_yticklabels([])
    ax_amb.set_ylim(*ax_gene.get_ylim())
    ax_amb.tick_params(axis='x', labelsize=FS_ROW + 1.5)
    ax_amb.set_xlabel('positions')
    panel_head(ax_amb, 'd', 'testable positions', pad=18.0)
    S.despine(ax_amb)

    n_fail = int(sq['failed'].sum()) if 'failed' in sq.columns else 0
    # The values the prose used to state by hand. Each is derived here, once, and
    # interpolated wherever it is quoted -- caption, panel text and sidecar.
    f2 = comp[fld[0]]
    f2_worst, f2_best = f2.idxmax(), f2.idxmin()
    mean_field = (sq[[c for c in ('field2', 'field3', 'field4') if c in sq.columns]]
                  .mul([2, 3, 4][:len([c for c in ('field2', 'field3', 'field4')
                                       if c in sq.columns])]).sum(axis=1)
                  / sq[[c for c in ('field2', 'field3', 'field4')
                        if c in sq.columns]].sum(axis=1))
    mf = sq.assign(_m=mean_field).groupby('platform')['_m'].mean()
    n_cr = sq['call_rate'].nunique() if 'call_rate' in sq.columns else 0
    top2 = (sq['call_rate'].value_counts(normalize=True).head(2).sum()
            if 'call_rate' in sq.columns else float('nan'))
    n_amb_records = int(sq['ambiguous'].sum()) if 'ambiguous' in sq.columns else 0
    hz = gq.set_index('gene')['hemizygous']
    hz_top = hz.idxmax()
    # top_pad has to clear the legends now parked in the title slot.
    lay_out(fig, list(axes.ravel()), plot_h=plot_h, top_pad=0.60,
            wspace=0.26, hspace=0.30, right=0.965,
            title='Typing completeness, by sequencing platform and by gene')
    fig.savefig(args.out_png)
    plt.close(fig)

    figure_doc.write_doc(
        args.out_png,
        title='HLA typing quality, by sequencing platform',
        question='Does typing quality differ between the platforms, and therefore '
                 'between cases and controls, in a way that could imitate a finding?',
        panels=[
            ('a', 'Field resolution by platform',
             f'The share of calls resolving to each IMGT field depth, with the 2-field share '
             f'printed. Stacked rather than averaged: the per-platform mean field depth spans '
             f'only {mf.min():.2f}-{mf.max():.2f}, a range nobody can act on, while the '
             f'composition separates the platforms cleanly - {f2.max():.1%} of calls stop at '
             f'2 fields on {f2_worst} against {f2.min():.1%} on {f2_best}. Shorter reads carry '
             f'less information for separating similar alleles, so a read-length effect would '
             f'appear here first. Field depths absent from the cohort are not drawn as '
             f'zero-width segments; those present are '
             f'{", ".join(c.replace("field", "") + "-field" for c in fld)}.'),
            ('b', 'Ambiguity by platform',
             f'Samples grouped by how many genes HLA-HD emitted several candidate pairs for '
             f'({n_amb_records:,} (sample, gene) records across the cohort). It is the closest '
             f'thing to a per-call confidence the tool provides, and it was previously computed '
             f'but never plotted.'),
            ('c', 'Outcome by gene',
             f'The four outcomes, as a composition. `called` and `hemizygous` are both '
             f'successes and their ratio is biology - {hz_top} is hemizygous in '
             f'{hz.max():.0%} of samples and A in {hz.get("A", float("nan")):.0%}, because '
             f'{hz_top} is nearly monomorphic in this population. `not_typed` means the gene is '
             f'absent from that haplotype, which is why DRB3/4/5 look broken on a call-rate '
             f'axis when `failed` is {n_fail} for every gene in the cohort.'),
            ('d', 'Polymorphic residue positions by gene',
             f'How many IMGT positions actually vary in this cohort - the number of testable '
             f'sites each gene contributes, {int(sp.sum()):,} in total. Drawn on panel (c)\'s '
             f'row order so a gene sits at the same height in both.'),
        ],
        numbers={
            'samples': f'{len(sq):,}',
            'platforms': f'{len(plats)}',
            '(sample, gene) failures': f'{n_fail}',
            'distinct call-rate values': f'{n_cr}',
            'share of samples on the top two': f'{top2:.1%}',
            '2-field share, worst platform': f'{f2.max():.2%} ({f2_worst})',
            '2-field share, best platform': f'{f2.min():.2%} ({f2_best})',
            'mean field depth, range over platforms': f'{mf.min():.2f}-{mf.max():.2f}',
            'ambiguous (sample, gene) records': f'{n_amb_records:,}',
            'genes carrying residues': f'{int((sites.polymorphic > 0).sum())}',
            'polymorphic positions, total': f'{int(sp.sum()):,}',
        },
        interpretation=(
            'Read this by platform. Cases were sequenced on DNBSeq and NovaSeq, '
            'controls entirely on HiSeqX, with no overlap, so any case/control '
            'difference in typing quality is confounded by construction and cannot '
            'be attributed. Platform is the axis on which a technical gradient is '
            'visible as itself. The gradient here is real but small — the platform '
            'ordering in (a) and (b) is the same one, and the weakest platform is the '
            'one sequenced shallowest. Its size is the size of the problem, stated. '
            'Whether it reaches the CALLS is a different question and a different '
            'figure: see typing_confound.png. The cohort is ancestry-restricted, so the '
            'samples counted here are the ones every table in this component describes.'),
        limits=[
            f'Completeness is saturated on this cohort and therefore a weak quality signal: '
            f'call rate takes {n_cr} distinct values across the whole cohort and {top2:.0%} of '
            f'samples sit on two of them. That is why every panel here is a composition and no '
            f'panel is a mean - a mean over a saturated metric hides the only variation there '
            f'is.',
            'This figure says nothing about typing ACCURACY. It measures whether a '
            'call was produced and how deeply it resolved, not whether it is right. '
            'Accuracy needs an external truth set.',
            'A flat call rate does not rule out a systematic bias toward particular '
            'alleles, which would need the same truth set to detect.',
        ],
        methods_ref='../../docs/METHODS.md')

    print(f'[plot_typing_qc] {len(sq)} samples over {len(plats)} platform(s) -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in plot_typing_qc: {e}', file=sys.stderr)
        sys.exit(1)
