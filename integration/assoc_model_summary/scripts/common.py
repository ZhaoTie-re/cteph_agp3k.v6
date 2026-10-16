#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Shared vocabulary and helpers of assoc_model_summary.
#             - the three scans compared (fixed, grm_pc, grm_only): labels,
#               covariate text and colours, so no figure or README can drift
#             - cohort order, display names and colours
#             - lambda_GC, MAF bins, the QQ thinning and its Beta band
#             - the shared plot_style of analysis/_shared, imported from the
#               staged file so the look matches cohort_manhattan.png
# Inputs  : -
# Outputs : -
# ---------------------------------------------------------------------------
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

CHI2_MEDIAN_1DF = 0.4549364231195729
GW_ALPHA, SUGGEST_ALPHA = 5e-8, 1e-5
# Study-wide MAF. The upstream filter is MAF >= 0.01 in the fixed-model CONTROLS, so a few percent of
# variants are below 0.01 overall; they are kept in their own bin (reported, not plotted).
MAF_EDGES = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
MAF_LABELS = ['<0.01', '0.01–0.02', '0.02–0.05', '0.05–0.10', '0.10–0.20', '0.20–0.50']

# The scans this component compares. 'upstream' is the model name inside the
# upstream component's tables (model_peaks_annotation.tsv, scan_qc.tsv).
MODELS = {
    'fixed': dict(name='Fixed model', engine='PLINK2 logistic regression', title='Fixed model (PLINK2 logistic)',
                  covar='SEX + 10 PCs', samples='relatives removed (PI_HAT > 0.2)',
                  short='Fixed (PLINK2)', color='#B08A3C', component='assoc_plink2'),
    'grm_pc': dict(name='Random model', engine='SAIGE logistic mixed model', title='Random model (SAIGE mixed model)',
                   covar='SEX + 10 PCs + full GRM', samples='relatives retained',
                   short='GRM + PCs (SAIGE)', color='#1B3B5F', component='assoc_saige'),
    'grm_only': dict(name='Random model, no PCs', engine='SAIGE logistic mixed model',
                     title='Random model without PCs (SAIGE)',
                     covar='SEX + full GRM', samples='relatives retained',
                     short='GRM only (SAIGE)', color='#3F8A6E', component='assoc_saige'),
}

COHORT_MARKERS = ['o', 's', '^', 'D', 'v']


def import_style(plot_style_path):
    """The shared plot_style module, from its staged copy (a symlink into analysis/_shared/scripts)."""
    d = str(Path(plot_style_path).resolve().parent)
    if d not in sys.path:
        sys.path.insert(0, d)
    import plot_style
    return plot_style


def setup(S, bump=1):
    """cohort_manhattan's rcParams ('paper' scale), one point larger: these figures are projected on a
    ~10-inch slide rather than printed at 7.2 inches."""
    import matplotlib.pyplot as plt
    S.setup_style('paper')
    for k in ('font.size', 'axes.titlesize', 'axes.labelsize', 'xtick.labelsize', 'ytick.labelsize',
              'legend.fontsize'):
        plt.rcParams[k] = plt.rcParams[k] + bump
    return 9 + bump                                    # panel-letter size


def csv_list(s):
    return [x.strip() for x in str(s).split(',') if x.strip()]


def kv(items):
    """['a=b', 'c=d'] or 'a=b,c=d' -> {a: b, c: d}."""
    if isinstance(items, str):
        items = csv_list(items)
    out = {}
    for it in items or []:
        k, _, v = it.partition('=')
        out[k.strip()] = v.strip()
    return out


def cohort_table(cohorts, labels, colors):
    cohorts, labels, colors = csv_list(cohorts), csv_list(labels), csv_list(colors)
    assert len(cohorts) == len(labels) == len(colors), (cohorts, labels, colors)
    return {c: dict(label=l, color='#' + col.lstrip('#'), marker=COHORT_MARKERS[i % len(COHORT_MARKERS)])
            for i, (c, l, col) in enumerate(zip(cohorts, labels, colors))}


def lambda_gc(p):
    p = np.asarray(p, dtype=float)
    return float(stats.chi2.isf(np.median(p), 1) / CHI2_MEDIAN_1DF)


def maf_bin(maf):
    """Bin index: (lo, hi] except the first, which is [0, 0.01)."""
    return np.digitize(maf, MAF_EDGES[1:-1], right=False)


def qq_points(nlp_desc, n_tail=20000, n_body=30000):
    """Thinned QQ points from the observed -log10 P sorted in DESCENDING order.

    Every one of the n_tail smallest P is kept (the tail is the part that is read); the body is a
    log-spaced sample of ranks, which draws the same curve. The ranks are returned so the band is
    computed at exactly the plotted points."""
    n = len(nlp_desc)
    head = np.arange(1, min(n, n_tail) + 1)
    body = np.unique(np.geomspace(min(n, n_tail) + 1, n, num=n_body).astype(np.int64)) if n > n_tail else []
    rank = np.unique(np.concatenate([head, np.asarray(body, dtype=np.int64)]))
    exp = -np.log10((rank - 0.5) / n)
    obs = nlp_desc[rank - 1].astype(float)
    return rank, exp, obs


def qq_band(rank, n):
    """Pointwise 95% band of the k-th smallest of n uniform P: Beta(k, n - k + 1)."""
    lo = -np.log10(stats.beta.ppf(0.975, rank, n - rank + 1))
    hi = -np.log10(stats.beta.ppf(0.025, rank, n - rank + 1))
    return lo, hi


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def chrom_key(c):
    c = str(c).replace('chr', '')
    return {'X': 23, 'Y': 24, 'MT': 25, 'M': 25}.get(c, int(c) if c.isdigit() else 99)


def fmt_p(p):
    """P for plain text: 3.8 × 10^-8 written as 3.8e-8."""
    if p >= 0.001:
        return f'{p:.3g}'
    m, e = f'{p:.1e}'.split('e')
    return f'{m} × 10^{int(e)}'


def p_md(p):
    """P for Markdown with a superscript exponent."""
    if p >= 0.001:
        return f'{p:.3g}'
    m, e = f'{p:.1e}'.split('e')
    return f'{m} × 10<sup>{int(e)}</sup>'


def log(msg):
    print(f'[assoc_model_summary] {msg}', file=sys.stderr, flush=True)
