#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Calibration of both statistics in every cell of the design, on one
#           canvas: variant sets down, cohorts across, a QQ per cell.
#
#           A QQ plot answers one question -- does the null hold -- and the only
#           reference it needs is the identity line and the null's own
#           pointwise band. No threshold line is drawn: a threshold is a decision
#           about individual genes and belongs on the scan figures; on a QQ it
#           reads as a claim about the distribution, which it is not.
#
#           lambda_GC is the one number each cell states, per statistic. Its split
#           by set size (n_var = 2 against >= 3) is a diagnostic of the discrete
#           statistic, not of the scan, and it is tabulated in the README rather
#           than drawn: two lambdas per curve on a small panel say less than one.
#
#           Every cell carries its OWN denominator, so heights are never compared
#           between cells; the README table gives n_genes_mapped and the threshold
#           per cell for the reader who needs them.
# Component: assoc_gene
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

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_shared' / 'scripts'))
import plot_style as S                   # noqa: E402
import figure_doc                        # noqa: E402
import vocab as V                        # noqa: E402

PLOT_H = 5.0
LETTERS = 'abcdefghijklmnop'
METHODS_REF = '../../../docs/METHODS.md'


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gene-scan-all', required=True)
    p.add_argument('--scan-qc-all', required=True)
    p.add_argument('--out-png', default='calibration.png')
    p.add_argument('--out-md', default='README.md')
    return p.parse_args()


def qq_curve(ax, nlp_values, color):
    obs = np.sort(np.asarray(nlp_values, dtype=float))[::-1]
    exp = V.qq_expected(len(obs))
    ax.plot(exp, obs, lw=1.0, color=color, zorder=3, rasterized=True)
    return float(max(exp.max(), obs.max()))


def main():
    a = parse_args()
    S.setup_style('paper')
    scan = pd.read_csv(a.gene_scan_all, sep='\t', dtype=str)
    qc = pd.read_csv(a.scan_qc_all, sep='\t', dtype=str)
    if scan.empty:
        raise SystemExit(f'ABORT: {a.gene_scan_all} is empty')
    bad = set(scan['method']) - set(V.METHODS)
    if bad:
        raise SystemExit(f'ABORT: {a.gene_scan_all} carries method(s) {sorted(bad)}')
    scan['nlp'] = V.nlp(pd.to_numeric(scan['pvalue'], errors='coerce').fillna(1.0))
    scan = scan[np.isfinite(pd.to_numeric(scan['pvalue'], errors='coerce'))]
    for c in ('lambda_gc', 'lambda_gc_nvar2', 'lambda_gc_nvar_ge3', 'n_genes_mapped',
              'bonferroni_threshold'):
        qc[c] = pd.to_numeric(qc[c], errors='coerce')

    cohorts = V.order_cohorts(set(scan['cohort']))
    clabel = V.cohort_labels(cohorts)
    strata = V.order_strata(scan)

    fig, axes = plt.subplots(len(strata), len(cohorts),
                             figsize=(S.COL_DOUBLE, PLOT_H + 2.0), squeeze=False)
    cells, corner = [], {}
    for r, stratum in enumerate(strata):
        for k, cohort in enumerate(cohorts):
            ax = axes[r][k]
            lim, lam, n_max = 1.0, {}, 0
            for m in V.METHODS:
                d = scan[(scan['cohort'] == cohort) & (scan['stratum'] == stratum) &
                         (scan['method'] == m)]
                if d.empty:
                    continue
                lim = max(lim, qq_curve(ax, d['nlp'], V.METHOD_COLOR[m]))
                n_max = max(n_max, len(d))
                q = qc[(qc['cohort'] == cohort) & (qc['stratum'] == stratum) & (qc['method'] == m)]
                if len(q):
                    lam[m] = q.iloc[0]
            if n_max:
                # The null's pointwise 95 % band, under the curves.
                exp, lo, hi = V.qq_band(n_max)
                ax.fill_between(exp, lo, hi, color=S.NEUTRAL, alpha=0.55, lw=0, zorder=1)
            lim *= 1.06
            ax.plot([0, lim], [0, lim], lw=0.8, color=S.REFERENCE, zorder=2)
            ax.set_box_aspect(1)
            ax.set_xlim(0, lim)
            ax.set_ylim(0, lim)
            if r == len(strata) - 1:
                ax.set_xlabel(f'expected {S.NEGLOG10P}')
            if k == 0:
                ax.set_ylabel(f'{V.stratum_label(stratum)}\nobserved {S.NEGLOG10P}')
            S.despine(ax)
            corner[(r, k)] = lam
            row = {'cohort': cohort, 'stratum': stratum}
            any_q = next(iter(lam.values()), None)
            row['n_genes_mapped'] = int(any_q['n_genes_mapped']) if any_q is not None else None
            row['bonferroni_threshold'] = (float(any_q['bonferroni_threshold'])
                                           if any_q is not None else None)
            for m in V.METHODS:
                for col in ('lambda_gc', 'lambda_gc_nvar2', 'lambda_gc_nvar_ge3'):
                    row[f'{col}_{m}'] = float(lam[m][col]) if m in lam else None
            cells.append(row)

    lam_all = [c[f'lambda_gc_{m}'] for c in cells for m in V.METHODS
               if c[f'lambda_gc_{m}'] is not None]
    lam_rng = f'{min(lam_all):.2f}–{max(lam_all):.2f}' if lam_all else 'n/a'

    S.caption_block(
        fig, plot_h=PLOT_H,
        title=(f'{len(cohorts)} cohorts × {len(strata)} variant sets: '
               f'{S.LAMBDA_GC} spans {lam_rng} over the {len(cells)} scans.'),
        panels=[f'{V.stratum_label(s)}, {clabel[c]}.' for s in strata for c in cohorts],
        notes=('Grey: pointwise 95 % band of the null; each cell has its own gene family and '
               'threshold, so heights are not comparable between cells. README.md'),
        letters=LETTERS, top_pad=0.50, hspace=0.40, wspace=0.28, left='auto', right=0.985,
        margin_axes=[ax for row in axes for ax in row])

    # ---- everything below measures the rendered figure.
    for r, stratum in enumerate(strata):
        for k, cohort in enumerate(cohorts):
            ax = axes[r][k]
            lam = corner[(r, k)]
            parts = [f'{V.METHOD_SHORT[m]} {float(lam[m]["lambda_gc"]):.2f}'
                     for m in V.METHODS if m in lam]
            one = f'{S.LAMBDA_GC}  ' + ' · '.join(parts)
            w_in, _h = V.text_extent_in(fig, one, V.MIN_FONT)
            ax_w = ax.get_window_extent().width / fig.dpi
            txt = one if w_in <= 0.9 * ax_w else f'{S.LAMBDA_GC}\n' + '\n'.join(parts)
            # Upper-left corner: the curves run along the diagonal and the band
            # hugs it, so this corner is empty by construction.
            ax.text(0.04, 0.96, txt, transform=ax.transAxes, ha='left', va='top',
                    fontsize=V.MIN_FONT, color=S.INK_SOFT, zorder=6, linespacing=1.3)
            S.panel_tag(ax, LETTERS[r * len(cohorts) + k],
                        title=clabel[cohort] if r == 0 else None, tsize=8)
    handles = [Line2D([], [], color=V.METHOD_COLOR[m], lw=1.6, label=V.METHOD_LABEL[m])
               for m in V.METHODS]
    handles.append(Line2D([], [], color=S.NEUTRAL, lw=6, alpha=0.8, label='null 95 % band'))
    # Lower-right corner of the first cell: below the identity line, always empty.
    S.legend_inside(axes[0][0], handles, loc='lower right', ncol=1)
    fig.savefig(a.out_png)
    plt.close(fig)

    cell_tbl = pd.DataFrame(cells)
    figure_doc.write_doc(
        a.out_png, out_path=a.out_md,
        subject=f'**Figure file:** `{Path(a.out_png).name}`',
        title='Calibration across cohorts and variant sets',
        question=('Is the null calibrated for both statistics in every cohort × variant-set cell '
                  'of the design, and does calibration change as controls are added or the '
                  'impact filter is loosened?'),
        panels=[(LETTERS[r * len(cohorts) + k],
                 f'{V.stratum_label(strata[r])}, {clabel[cohorts[k]]}',
                 'Observed against expected −log10 P for CMC burden and SKAT-O over the genes '
                 'each returned in this cell; expected quantiles are (i − 0.5)/n. The grey band '
                 'is the pointwise 95 % interval of the uniform null, Beta(i, n − i + 1). The '
                 'corner text gives each statistic\'s genomic-control inflation factor, read '
                 'from scan_qc.tsv. No threshold is drawn on a QQ; each cell\'s own Bonferroni '
                 'threshold and gene count are in the table below.')
                for r in range(len(strata)) for k in range(len(cohorts))],
        interpretation=(
            f'λ_GC spans {lam_rng} over the {len(cells)} cells. Read across a row: the cohorts '
            'are nested and differ almost entirely in controls, so a trend as the sample set '
            'widens is a statement about the controls added, not a replication. Read down a '
            'column: the variant sets are nested, and the wider set\'s extra genes carry only '
            'LOW-impact variants. Within a cell the two curves share genes, samples and '
            'denominator, so their difference is the statistic. λ is reported, not targeted: '
            'the statistic is a count over a handful of minor alleles, and its median-based '
            'λ is inflated by the two-variant genes (table column lambda_gc_nvar2 against '
            'lambda_gc_nvar_ge3); the residual 1.1–1.2 over larger sets is consistent with '
            'the platform confound and cannot be separated from it here (METHODS §4.4).'),
        numbers=[('cohorts', len(cohorts)), ('variant sets', len(strata)),
                 ('cells drawn', len(cells)), ('lambda_GC range', lam_rng)],
        tables=[('Denominator and calibration per cell', cell_tbl)],
        reading=[
            'Compare the two curves within a cell first; that comparison holds everything '
            'but the statistic fixed.',
            'Then read across a row (adding controls) and down a column (loosening the '
            'impact filter).',
            'Never compare heights between cells: each has its own gene family and threshold.',
        ],
        limits=[
            'It is not replication; the cohorts are nested and share nearly all cases.',
            'It does not establish any association; the platform confound is in every curve.',
            'λ_GC over non-independent, discrete gene tests is a coarse instrument.',
        ],
        defs=[],
        model=('per cell: rvtest CMC burden and SKAT-O, logistic, sex + ancestry PCs, on the '
               'cell\'s own snpEff map; band = Beta(i, n − i + 1) pointwise 95 % interval'),
        methods_ref=METHODS_REF)
    print(f'[plot_grid] {len(strata)} × {len(cohorts)} = {len(cells)} cell(s), '
          f'lambda_GC {lam_rng} -> {a.out_png}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:                       # noqa: BLE001
        print(f'Error in plot_grid: {e}', file=sys.stderr)
        sys.exit(1)
