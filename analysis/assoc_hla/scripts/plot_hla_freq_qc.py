#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The typing artefact check, applied to the HITS THEMSELVES.
#             (a) our allele frequency against the reference panel's
#             (b) -log10 P against minor-allele COUNT
#             (c) the significant alleles the reference has NEVER observed
#
#           EVERY OTHER FIGURE IN THIS COMPONENT ASSUMES THE TYPING IS RIGHT.
#           An association test cannot tell a real allele from a mis-called one:
#           if HLA-HD systematically assigns a handful of case chromosomes to an
#           allele they do not carry, the result is a perfectly clean, highly
#           significant odds ratio for an allele that does not exist in these
#           people. Nothing in a P-value, an interval or a lambda detects that.
#           The only external handle available is the population frequency: an
#           allele we call at 1 % that a large Japanese reference panel calls at
#           1 % is a plausible allele, and one the reference has never seen at
#           all is a candidate artefact — and the association scan will happily
#           report the second as a finding.
#
#           WHY (b) IS ON THE SAME FIGURE. The other way a hit can be an artefact
#           is arithmetic rather than biological: at MAC in the tens, a handful of
#           chromosomes moving between arms produces a large odds ratio and a
#           small P. Plotting evidence against count says whether the hit list is
#           tracking the data or tracking the rarity floor. The two panels are
#           the same question asked of the same alleles, which is why they share
#           a figure rather than sitting in two.
#
#           (c) IS THE ONE PANEL THAT MUST NEVER BE OMITTED. `in_reference = 0`
#           means the reference panel assigned NO chromosome to that allele's P
#           group. An allele a panel of tens of thousands of Japanese has never
#           observed, appearing at a frequency high enough to be tested here and
#           significant, is more likely a typing artefact than a discovery, and
#           the figure says so in its own caption rather than leaving it to the
#           document. When no significant allele is flagged, that is a positive
#           result and the panel prints it — a blank panel would read as a
#           missing analysis.
#
#           ALLELE MARKERS ONLY. `freq_reference` is blank for an allele the
#           reference has never observed AND for every residue marker, because
#           the reference is an ALLELE panel and "not in it" is not a statement
#           about a residue. Those two blanks mean different things and are never
#           merged here: residues are excluded from the figure, absent alleles
#           are its subject.
# Component: assoc_hla
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402
from matplotlib.lines import Line2D      # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
SHARED = Path(__file__).resolve().parent.parent.parent / '_shared' / 'scripts'
sys.path.insert(0, str(SHARED))
import plot_style as S                   # noqa: E402
import figure_doc                        # noqa: E402

NOMINAL = 0.05        # the dashed second line on the -log10 P axis
ROW_IN = 0.30         # inches per row of panel (c)
MAX_FLAGGED_ROWS = 14  # a cap on the FIGURE only; the caption states it if it bites
QC_COLS = ('id', 'marker_class', 'gene', 'mac', 'maf', 'af', 'in_reference',
           'freq_reference', 'tested', 'call_rate')
CANON_NUM = ('POS', 'A1_FREQ', 'OBS_CT', 'OR', 'LOG(OR)_SE', 'L95', 'U95', 'Z_STAT', 'P')


def parse_args():
    p = argparse.ArgumentParser(
        description='Frequency and count QC of the HLA allele hits themselves.')
    p.add_argument('--sumstat', required=True,
                   help='one cohort/model sumstats.tsv, ALLELE class only')
    p.add_argument('--marker-qc', required=True,
                   help='marker_qc.tsv: id, marker_class, gene, mac, maf, af, in_reference, '
                        'freq_reference, tested, call_rate')
    p.add_argument('--cohort', required=True)
    p.add_argument('--model', required=True)
    p.add_argument('--alpha', type=float, default=0.0,
                   help='significance threshold. 0 (the default) means Bonferroni 0.05 / M '
                        'over the M alleles with a usable fit.')
    p.add_argument('--reference-name', default='jMorp 61KJPN-HLA',
                   help='what the frequency reference is called on the figure')
    p.add_argument('--reference-n', type=int, default=61424,
                   help='individuals in that reference, quoted in the caption. The claim '
                        '"never observed" is only as strong as this number, so it is printed '
                        'from the parameter rather than written into the prose.')
    p.add_argument('--out-png', required=True)
    return p.parse_args()


def empty_panel(ax, msg):
    """An honest blank: say why there is nothing rather than draw an empty box."""
    ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha='center', va='center',
            color=S.INK_SOFT, fontsize=plt.rcParams['legend.fontsize'])
    ax.set_xticks([]); ax.set_yticks([])
    ax.grid(False)
    # No frame either: an empty axes box reads as a rendering failure, while a
    # bare sentence reads as the statement it is.
    for sp in ax.spines.values():
        sp.set_visible(False)


def read_inputs(args):
    """sumstats x marker_qc, allele markers only, ERRCODE applied.

    Joined on the marker id, which build_hla_markers.py minted and
    hla_to_sumstats.py asserts against the .bim, so the key is the same string on
    both sides by construction.
    """
    d = pd.read_csv(args.sumstat, sep='\t', dtype={'ID': str, 'ERRCODE': str, 'gene': str,
                                                   'position': str, '#CHROM': str})
    for c in CANON_NUM:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors='coerce')
    err = d['ERRCODE'].fillna('.').astype(str).str.strip() if 'ERRCODE' in d.columns \
        else pd.Series('.', index=d.index)
    ok = err.eq('.') & d['P'].gt(0) & d['P'].le(1)
    n_dropped = int((~ok).sum())
    d = d[ok].copy()
    d['marker'] = d['ID'].astype(str)

    q = pd.read_csv(args.marker_qc, sep='\t', dtype={'id': str, 'marker_class': str,
                                                     'gene': str})
    missing = [c for c in QC_COLS if c not in q.columns]
    if missing:
        raise SystemExit(f'ABORT: {args.marker_qc} has no {missing} column(s); it was not '
                         f'written by hla_marker_qc.py / build_hla_markers.py')
    if 'marker_class' in q.columns:
        q = q[q['marker_class'].astype(str) == 'allele'].copy()
    for c in ('mac', 'maf', 'af', 'call_rate', 'freq_reference'):
        # freq_reference is BLANK both for an allele the reference never saw and
        # for every residue; coerced to NaN, and the residues are already gone.
        q[c] = pd.to_numeric(q[c], errors='coerce')
    q['in_reference'] = pd.to_numeric(q['in_reference'], errors='coerce')
    m = d.merge(q[list(QC_COLS)], left_on='marker', right_on='id', how='inner',
                suffixes=('', '_qc'))
    return m, n_dropped, len(d), len(q)


def main():
    args = parse_args()
    S.setup_style('paper')
    m, n_dropped, n_fit, n_qc = read_inputs(args)
    REF = args.reference_name
    n_alleles = int(m['marker'].nunique())
    alpha = args.alpha if args.alpha > 0 else (0.05 / max(n_alleles, 1))
    m['significant'] = m['P'] < alpha
    m['absent'] = m['in_reference'].fillna(-1) == 0        # never observed by the reference
    m['neglogp'] = -np.log10(m['P'].to_numpy(float))

    has_ref = m['freq_reference'].notna()
    flagged = (m[m.significant & m.absent].sort_values('P')
               if len(m) else m.iloc[0:0])
    n_flag_all = len(flagged)
    flagged = flagged.head(MAX_FLAGGED_ROWS)
    capped = n_flag_all > len(flagged)

    h_ab = 2.75
    h_c = max(0.80, len(flagged) * ROW_IN) + 0.35
    plot_h = h_ab + h_c + 1.60

    fig = plt.figure(figsize=(S.COL_DOUBLE, plot_h))
    gs = fig.add_gridspec(2, 1, height_ratios=[h_ab, h_c])
    # A narrow strip for the alleles that HAVE no reference frequency, beside the
    # scatter that needs one. They cannot be drawn at x = 0 — the reference does
    # not say the frequency is zero, it says the allele was never seen — and
    # dropping them would hide the whole subject of panel (c).
    gsa = gs[0].subgridspec(1, 3, width_ratios=[0.20, 1.0, 1.0], wspace=0.30)
    ax_ax = fig.add_subplot(gsa[0, 0])
    ax_a = fig.add_subplot(gsa[0, 1], sharey=ax_ax)
    ax_b = fig.add_subplot(gsa[0, 2])
    ax_c = fig.add_subplot(gs[1])

    thr = -np.log10(alpha)

    # ── (a) our frequency against the reference's ──────────────────────────
    if not len(m):
        empty_panel(ax_ax, '')
        empty_panel(ax_a, 'no allele with a usable fit\nand a QC row')
    elif not has_ref.any():
        # Every allele is absent from the reference, or the panel was never
        # attached. Either way the identity plot has no x to draw against, and
        # saying so is the honest rendering.
        empty_panel(ax_a, f'no allele has a frequency in\n{REF} — nothing to compare '
                          f'against')
        empty_panel(ax_ax, '')
    else:
        sub = m[has_ref]
        hit = sub.significant.to_numpy()
        x = sub['freq_reference'].to_numpy(float)
        y = sub['af'].to_numpy(float)
        hi = float(np.nanmax(np.concatenate([x, y]))) * 1.12 + 0.005
        ax_a.plot([0, hi], [0, hi], color=S.REFERENCE, lw=1.0, ls='--', zorder=2)
        ax_a.scatter(x[~hit], y[~hit], s=13, color=S.DATA, alpha=0.75,
                     edgecolors='none', zorder=4, rasterized=True)
        ax_a.scatter(x[hit], y[hit], s=30, color=S.ACCENT, edgecolors='white',
                     linewidths=0.4, zorder=6)
        # A name on a point in the right-hand half is written to the LEFT of its
        # marker: an HLA allele id is wide enough to run out of the panel and
        # into the next one, and an annotation is not clipped by default.
        ann = [ax_a.annotate(t, (xi, yi), textcoords='offset points',
                             xytext=(-5 if xi > 0.55 * hi else 5, 3),
                             ha='right' if xi > 0.55 * hi else 'left',
                             fontsize=plt.rcParams['legend.fontsize'] - 1.5,
                             color=S.INK_SOFT)
               for xi, yi, t in zip(x[hit], y[hit], sub.marker[hit])]
        if ann:
            S.spread_labels(ax_a, ann, axis='y')
        ax_a.set_xlim(0, hi); ax_a.set_ylim(0, hi)
        ax_a.set_xlabel(f'frequency in {REF}')
        plt.setp(ax_a.get_yticklabels(), visible=False)
        S.despine(ax_a)

        # the strip: alleles with no reference frequency at all
        nof = m[~has_ref]
        if len(nof):
            # Deterministic rank offsets, never a random jitter (FIGURES.md rule 6).
            xs = np.linspace(0.28, 0.72, len(nof)) if len(nof) > 1 else np.array([0.5])
            hit2 = nof.significant.to_numpy()
            yv = nof['af'].to_numpy(float)
            ax_ax.scatter(xs[~hit2], yv[~hit2], s=13, marker='s', facecolors='white',
                          edgecolors=S.DATA, linewidths=0.8, zorder=4)
            ax_ax.scatter(xs[hit2], yv[hit2], s=30, marker='s', facecolors='white',
                          edgecolors=S.ACCENT, linewidths=1.2, zorder=6)
        else:
            # Saying "none" is a result; an empty strip is ambiguous between
            # "none" and "this was not computed".
            ax_ax.annotate('none', xy=(0.5, 0.5), xycoords='axes fraction', ha='center',
                           va='center', color=S.INK_SOFT,
                           fontsize=plt.rcParams['legend.fontsize'] - 1)
        ax_ax.set_xlim(0, 1)
        ax_ax.set_xticks([])
        ax_ax.set_xlabel('absent')
        ax_ax.set_ylabel(f'{args.cohort} allele frequency')
        S.despine(ax_ax, grid_axis='y')

    # ── (b) evidence against count ─────────────────────────────────────────
    if not len(m):
        empty_panel(ax_b, 'no allele with a usable fit')
    else:
        mac = pd.to_numeric(m['mac'], errors='coerce').to_numpy(float)
        y = m['neglogp'].to_numpy(float)
        good = np.isfinite(mac) & (mac > 0) & np.isfinite(y)
        hit = m.significant.to_numpy() & good
        ab = m.absent.to_numpy() & good
        for sel, kw in ((good & ~hit & ~ab, dict(s=13, color=S.DATA, alpha=0.75,
                                                 edgecolors='none', zorder=4)),
                        (good & ~hit & ab, dict(s=18, marker='s', facecolors='white',
                                                edgecolors=S.DATA, linewidths=0.8, zorder=5)),
                        (hit & ~ab, dict(s=30, color=S.ACCENT, edgecolors='white',
                                         linewidths=0.4, zorder=6)),
                        (hit & ab, dict(s=34, marker='s', facecolors='white',
                                        edgecolors=S.ACCENT, linewidths=1.3, zorder=7))):
            if sel.any():
                ax_b.scatter(mac[sel], y[sel], **kw)
        if good.any():
            ax_b.set_xscale('log')
            ax_b.set_xlim(max(1.0, float(np.nanmin(mac[good])) * 0.7),
                          float(np.nanmax(mac[good])) * 1.5)
        S.tier_lines(ax_b, alpha, NOMINAL, label=True)
        ax_b.set_xlabel('minor-allele count (MAC), log scale')
        ax_b.set_ylabel(S.NEGLOG10P)
        S.despine(ax_b)
        # place_legend, not legend_above: three keys naming a reference panel do
        # not fit one row above a half-width panel, and legend_above deliberately
        # does not wrap — it would be clipped. This switches PLACEMENT instead,
        # into whichever corner holds the fewest points.
        S.place_legend(ax_b, [
            Line2D([], [], ls='none', marker='o', markersize=5.5, color=S.DATA,
                   label='in the reference'),
            Line2D([], [], ls='none', marker='s', markersize=5.5, markerfacecolor='white',
                   markeredgecolor=S.DATA, label='never observed there'),
            Line2D([], [], ls='none', marker='o', markersize=6, color=S.ACCENT,
                   label='significant')],
            data=(mac[good], y[good]) if good.any() else None)

    # ── (c) the significant alleles the reference has never observed ───────
    if not len(flagged):
        # A positive result, and it must be PRINTED. A blank panel here reads as
        # an analysis that was not done.
        msg = ('No allele clearing the threshold is absent from '
               f'{REF}: every significant allele is one the reference panel has observed.'
               if len(m) else 'no allele with a usable fit and a QC row')
        empty_panel(ax_c, msg)
    else:
        yv = np.arange(len(flagged))[::-1]
        af = flagged['af'].to_numpy(float)
        ax_c.barh(yv, af, height=0.5, color='white', edgecolor=S.ACCENT, linewidth=1.2,
                  hatch='///', zorder=3)
        ax_c.set_yticks(yv)
        ax_c.set_yticklabels(flagged.marker.tolist(), fontweight='bold')
        ax_c.set_ylim(-0.7, len(flagged) - 0.3)
        ax_c.set_xlim(0, float(np.nanmax(af)) * 3.6 if np.isfinite(np.nanmax(af)) else 1.0)
        ax_c.set_xlabel(f'{args.cohort} allele frequency (the reference has none)')
        for y, (_i, r) in zip(yv, flagged.iterrows()):
            orv = float(r['OR']) if pd.notna(r['OR']) else float('nan')
            l95 = float(r['L95']) if pd.notna(r.get('L95')) else float('nan')
            u95 = float(r['U95']) if pd.notna(r.get('U95')) else float('nan')
            ci = (f' ({l95:.2f}–{u95:.2f})' if np.isfinite(l95) and np.isfinite(u95) else '')
            ax_c.annotate(f"MAC {int(r['mac']):,}  ·  call rate {float(r['call_rate']):.3f}"
                          f"  ·  {S.OR_SYM} {orv:.2f}{ci}  ·  $P$ = {S.p_tex(r['P'])}",
                          xy=(float(r['af']), y), xytext=(6, 0), textcoords='offset points',
                          ha='left', va='center', color=S.INK_SOFT,
                          fontsize=plt.rcParams['legend.fontsize'] - 0.5)
        S.despine(ax_c, grid_axis='x')

    # ── caption ─────────────────────────────────────────────────────────────
    n_sig = int(m.significant.sum()) if len(m) else 0
    n_absent = int(m.absent.sum()) if len(m) else 0
    verdict = (f'{n_flag_all} of the {n_sig} significant allele(s) are absent from a panel of '
               f'{args.reference_n:,} individuals and are candidate typing artefacts, not '
               f'findings.' if n_flag_all else
               f'None of the {n_sig} significant allele(s) is absent from the reference panel.')
    S.caption_block(
        fig, plot_h=plot_h,
        title=(f'HLA allele hits against {REF}, {args.cohort} / {args.model}: ' + verdict),
        panels=[
            (f'cohort allele frequency against the frequency in {REF}, one point per tested '
             f'allele; the dashed line is equality and alleles with $P$ < {S.mathsci(alpha)} '
             f'are red and named. The narrow left-hand strip holds the {n_absent} allele(s) '
             f'the reference has NEVER observed, which have no x to be drawn at.'),
            (f'{S.NEGLOG10P} against minor-allele count, log axis. Open squares are the '
             f'alleles absent from the reference. The solid line is the threshold, the dashed '
             f'line nominal $P$ = {S.mathsci(NOMINAL)}.'),
            ('Every significant allele the reference has never observed, with its frequency '
             'here, its count, its call rate and its effect. '
             + (f'The {len(flagged)} strongest of {n_flag_all} are drawn. ' if capped else '')
             + 'An allele that a reference panel of '
             + f'{args.reference_n:,} Japanese individuals has never seen is more likely a '
               'typing artefact than a discovery, whatever its $P$-value.'),
        ],
        notes=(f'{args.cohort}, {args.model} model, HLA ALLELE markers only: {n_alleles} with '
               f'a usable fit and a QC row'
               + (f' ({n_dropped} excluded on ERRCODE)' if n_dropped else '')
               + f'. Threshold {S.sci(alpha)}'
               + ('' if args.alpha > 0 else f' = 0.05 / {n_alleles}, a Bonferroni over the '
                                            f'allele class alone')
               + f'. Reference: {REF}, {args.reference_n:,} individuals, matched on IPD-IMGT P '
                 f'groups; `freq_reference` is blank both for an allele the reference never '
                 f'observed and for every residue marker, and the two are never merged — '
                 f'residues are not in this figure. This is a POPULATION frequency comparison '
                 f'and not a measurement of per-sample typing accuracy. Full explanation: '
               + Path(args.out_png).stem + '.md'),
        top_pad=0.42, hspace=0.55, wspace=0.30, left='auto', right=0.965,
        margin_axes=[ax_ax, ax_b, ax_c])
    for ax, letter in ((ax_ax, 'a'), (ax_b, 'b'), (ax_c, 'c')):
        S.panel_tag(ax, letter)
    fig.savefig(args.out_png)
    plt.close(fig)

    cols = ['marker', 'gene', 'af', 'freq_reference', 'in_reference', 'mac', 'maf',
            'call_rate', 'OR', 'L95', 'U95', 'P']
    hits = (m[m.significant].sort_values('P')[cols] if len(m) else pd.DataFrame(columns=cols))
    figure_doc.write_doc(
        args.out_png,
        title=f'HLA allele hits against {REF} — {args.cohort} / {args.model}',
        question=('Are the alleles this scan reports as associated alleles that a large '
                  'reference population actually carries, at the frequency we call them — or '
                  'are they typing artefacts and rarity artefacts that an association test '
                  'cannot distinguish from findings?'),
        interpretation=(
            'An association test cannot tell a real allele from a mis-called one. If the typing '
            'systematically assigns a handful of case chromosomes to an allele they do not '
            'carry, the result is a clean odds ratio with a small P-value for an allele that is '
            'not there, and no P-value, interval or lambda_GC will reveal it. The only external '
            'handle available is population frequency. Panel (a) is that check drawn over every '
            'tested allele: agreement with the reference means the typing has not moved the '
            'frequency spectrum, and a significant allele sitting on the identity line is one '
            'the reference independently agrees exists at that frequency. Panel (b) covers the '
            'other artefact class, which is arithmetic rather than biological: at counts in the '
            'tens a few chromosomes moving between arms produce a large odds ratio and a small '
            'P. Panel (c) is the strongest single caution the component can issue — an allele '
            f'that {args.reference_n:,} Japanese individuals have never been observed to carry, '
            'called here often enough to be tested and significant, is more likely a typing '
            'artefact than a discovery. None of this proves any individual call wrong; it '
            'identifies which hits must not be reported without further evidence.'),
        panels=[
            ('a', 'Frequency against the reference',
             f'One point per tested allele: the frequency in {args.cohort} on y against the '
             f'frequency in {REF} on x, with the dashed identity line. Alleles clearing the '
             f'significance threshold are red and named. The narrow left-hand strip holds the '
             f'alleles the reference has NEVER observed: they have no reference frequency, and '
             f'drawing them at x = 0 would assert that the reference measured them at zero when '
             f'it did not measure them at all. Their x positions inside the strip are spread '
             f'deterministically by rank and carry no meaning.'),
            ('b', 'Evidence against count',
             '-log10 P against minor-allele count on a log axis. A hit list that tracks the '
             'count floor rather than the data shows up here as evidence concentrating at the '
             'left-hand end. Open squares mark alleles absent from the reference, so the two '
             'artefact classes can be seen together: a significant allele that is both rare and '
             'unknown to the reference is the weakest kind of hit this component can produce.'),
            ('c', 'Significant and never observed in the reference',
             'Every allele that clears the threshold and has in_reference = 0, with its '
             'frequency in this cohort, its minor-allele count, its call rate, its odds ratio '
             'and its P-value. in_reference = 0 means the reference assigned NO chromosome to '
             "that allele's P group. When the panel is empty it says so in words: that is a "
             'positive result, and a blank panel would read as an analysis that was not run.'),
        ],
        numbers=[('cohort', args.cohort), ('model', args.model),
                 ('marker class', 'allele'),
                 ('reference panel', f'{REF} ({args.reference_n:,} individuals)'),
                 ('alleles with a usable fit and a QC row', n_alleles),
                 ('fits excluded on ERRCODE', n_dropped),
                 ('alleles in the sumstats file', n_fit),
                 ('allele rows in marker_qc.tsv', n_qc),
                 ('alleles with a reference frequency', int(has_ref.sum()) if len(m) else 0),
                 ('alleles the reference has never observed', n_absent),
                 ('significance threshold', alpha),
                 ('threshold source',
                  'given with --alpha' if args.alpha > 0 else f'Bonferroni 0.05 / {n_alleles}'),
                 ('significant alleles', n_sig),
                 ('significant AND absent from the reference', n_flag_all)],
        tables=[('Every significant allele, with its reference frequency and count',
                 hits if len(hits) else None, 40),
                ('Significant alleles the reference has never observed',
                 flagged[cols] if len(flagged) else None, MAX_FLAGGED_ROWS)],
        reading=[
            'Start with panel (c). If it names an allele, that allele does not go into a result '
            'sentence without independent evidence, whatever panels (a) and (b) show.',
            'In (a), read the red points only. A red point on the dashed line is a significant '
            'allele whose frequency the reference agrees with; a red point far above it is an '
            'allele we call more often than the reference does, which is what a typing bias '
            'toward a common allele looks like.',
            'In (b), read left to right. Evidence that appears only at the smallest counts is '
            'evidence about the count floor.',
            'The left-hand strip in (a) is not "frequency zero in the reference". It is "the '
            'reference never observed this allele", which is a different statement and the '
            'reason the strip is separated from the axis.',
            'Nothing here measures per-sample typing accuracy. See the typing component for '
            'what the frequency comparison can and cannot establish.',
        ],
        limits=[
            'It compares POPULATION frequencies. A set of typing errors that happens to '
            'preserve the frequency spectrum is invisible to it, and it can never say that a '
            'particular sample was typed correctly.',
            'Agreement with the reference is not evidence that an association is real. It only '
            'removes one specific way for it to be false.',
            f'The reference is {args.reference_n:,} individuals from a Japanese panel matched on '
            f'IPD-IMGT P groups. A P group is coarser than a 2-field allele, so an allele that '
            f'differs from a reference allele only below P-group resolution is compared at the '
            f'group level, not at its own.',
            'in_reference = 0 is a flag, not a verdict. A genuinely rare allele absent from the '
            'reference by sampling alone is possible; the flag says the hit needs external '
            'evidence, not that it is wrong.',
            'Residue markers are excluded entirely. The reference is an allele panel, so a '
            'blank freq_reference on a residue means "not applicable", not "never observed", '
            'and merging the two would manufacture flags.',
        ],
        defs=['or', 'eaf', 'errcode'],
        model=S.FORMULAS['glm'])

    print(f'[plot_hla_freq_qc] {args.cohort}/{args.model}: {n_alleles:,} allele(s), '
          f'alpha {alpha:.3e}; {n_sig} significant, {n_absent} absent from {REF}, '
          f'{n_flag_all} both -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in plot_hla_freq_qc: {e}', file=sys.stderr)
        sys.exit(1)
