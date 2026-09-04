#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The one vocabulary every assoc_gene figure draws with.
#
#           Five figure scripts used to carry five copies of the method names,
#           the cohort short forms, the stratum labels and the state styles,
#           and the copies drifted ('interm.' in one figure, 'intermediate' in
#           the next; 'LMH' here, 'LOW+MODERATE+HIGH' there). A reader moving
#           between figures should meet one word for one thing, so every name,
#           colour and marker that means something in this component is defined
#           here once and imported.
#
#           TIER MEANS THE ROBUSTNESS TIER. Tier 1 / 2 / 3 are the integers
#           robust_genes.tier_of assigns to a (gene, stratum) from the BH calls
#           across cohorts and statistics. plot_style.TIER_STYLE and
#           plot_style.tier_lines belong to the variant components (genome-wide
#           / suggestive peaks on a -log10 P axis) and are never imported here:
#           a figure that drew them would be making a claim this component does
#           not make.
#
#           NO ABBREVIATIONS. Strata are written LOW+MODERATE+HIGH and
#           MODERATE+HIGH; cohorts as the suffix-stripped tags plot_style.shorten
#           computes from the run (narrow / intermediate / full). Where a label
#           does not fit it is rotated or the layout is widened, never shortened.
#
#           COLOUR BY MEANING, from the shared palette: the two statistics are
#           SERIES[0] / SERIES[1]; the four cell states are ACCENT (Bonferroni +
#           BH), ACCENT_LT (BH only), DATA (tested, not called) and REFERENCE
#           (not in map); cases are ACCENT and controls DATA wherever the two
#           groups are contrasted; strata, an ordered nesting, take COHORT_RAMP.
# Component: assoc_gene
# ---------------------------------------------------------------------------
import math
import sys
from pathlib import Path

import numpy as np
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_shared' / 'scripts'))
import plot_style as S                 # noqa: E402

# ── statistics ───────────────────────────────────────────────────────────────
METHODS = ('cmc', 'skato')
METHOD_LABEL = {'cmc': 'CMC burden', 'skato': 'SKAT-O'}
METHOD_SHORT = {'cmc': 'CMC', 'skato': 'SKAT-O'}
METHOD_COLOR = {'cmc': S.SERIES[0], 'skato': S.SERIES[1]}

# ── the four states of one (gene, cohort, stratum, statistic) cell ──────────
# The tuple is repeated LITERALLY in robust_genes.py; verify.sh §8 asserts the two
# are identical, so a state added on one side cannot go unnoticed on the other.
STATES = ('called_bonferroni', 'called_bh_only', 'not_significant', 'not_in_map')
CALLED_STATES = ('called_bonferroni', 'called_bh_only')
STATE_LABEL = {'called_bonferroni': 'Bonferroni + BH',
               'called_bh_only': 'BH only',
               'not_significant': 'tested, not called',
               'not_in_map': 'not in map'}
STATE_COLOR = {'called_bonferroni': S.ACCENT, 'called_bh_only': S.ACCENT_LT,
               'not_significant': S.DATA, 'not_in_map': S.REFERENCE}
# marker, size and fill; colour comes from STATE_COLOR so text and marks agree.
STATE_STYLE = {'called_bonferroni': dict(marker='o', s=30, filled=True),
               'called_bh_only':    dict(marker='o', s=22, filled=True),
               'not_significant':   dict(marker='o', s=9,  filled=True),
               'not_in_map':        dict(marker='x', s=16, filled=False)}


def state_handles(states=STATES):
    """Legend keys for the cell states, in the canonical order."""
    out = []
    for st in states:
        sty = STATE_STYLE[st]
        out.append(Line2D([], [], ls='none', marker=sty['marker'], color=STATE_COLOR[st],
                          markerfacecolor=STATE_COLOR[st] if sty['filled'] else 'none',
                          markeredgewidth=1.2 if st == 'not_in_map' else 0.6,
                          markersize=math.sqrt(sty['s']) * 0.95, label=STATE_LABEL[st]))
    return out


def draw_state(ax, x, y, state, zorder=3):
    """One cell mark, the same everywhere."""
    sty = STATE_STYLE[state]
    col = STATE_COLOR[state]
    if sty['filled']:
        ax.scatter([x], [y], marker=sty['marker'], s=sty['s'], c=col, linewidths=0,
                   zorder=zorder)
    else:
        ax.scatter([x], [y], marker=sty['marker'], s=sty['s'], c=col, linewidths=1.2,
                   zorder=zorder)


# ── cohorts and strata ───────────────────────────────────────────────────────
# The nesting order is a preference, not a requirement: unknown tags follow.
COHORT_PREFERENCE = ('narrow_mainland', 'intermediate_mainland', 'full_mainland')


def order_cohorts(seen):
    pref = [c for c in COHORT_PREFERENCE if c in seen]
    return pref + sorted(set(seen) - set(pref))


def cohort_labels(cohorts):
    """cohort tag -> the short form shared by every figure (suffix-stripped)."""
    return S.shorten(list(cohorts))


def order_strata(frame, stratum_col='stratum', set_col='set_name'):
    """Widest stratum first: the one with the most distinct tested sets."""
    n = frame.drop_duplicates([stratum_col, set_col]).groupby(stratum_col).size()
    return list(n.sort_values(ascending=False).index)


def stratum_label(tag):
    return tag.replace('_', '+').upper()


def stratum_colors(strata):
    """Ordered strata -> COHORT_RAMP by nesting (narrowest set darkest)."""
    ramp = S.COHORT_RAMP
    return {s: ramp[min(i, len(ramp) - 1)] for i, s in enumerate(reversed(list(strata)))}


# Where colour is free the cohort takes the COHORT_RAMP, narrowest darkest
# (plot_robust_genes (b), which also puts the cohort on the row, so the ramp
# restates the label rather than carrying it alone).
def cohort_colors(cohorts):
    """Ordered cohorts -> COHORT_RAMP, narrowest darkest."""
    ramp = S.COHORT_RAMP
    return {c: ramp[min(i, len(ramp) - 1)] for i, c in enumerate(cohorts)}


def stratum_handles(strata, labels=None, color=None):
    """Filled = the narrowest set, hollow = the wider one; the house rule."""
    col = color or S.INK_SOFT
    lab = labels or {s: stratum_label(s) for s in strata}
    out = []
    for i, s in enumerate(reversed(list(strata))):          # narrowest first
        out.append(Line2D([], [], ls='none', marker='o', color=col,
                          markerfacecolor=col if i == 0 else 'none',
                          markeredgewidth=0.6 if i == 0 else 1.1,
                          markersize=4.0, label=lab[s]))
    return out


# ── groups, impacts, tiers ───────────────────────────────────────────────────
CASE_COLOR, CONTROL_COLOR = S.ACCENT, S.DATA
IMPACT_ORDER = ('HIGH', 'MODERATE', 'LOW')
IMPACT_MARKER = {'HIGH': 'D', 'MODERATE': 'o', 'LOW': 's'}
TIER_LABEL = {1: 'Tier 1', 2: 'Tier 2', 3: 'Tier 3'}
TIER_DEF = {1: 'called in every cohort, and by both statistics in at least one',
            2: 'called in every cohort by at least one statistic',
            3: 'called in at least two cohorts'}

# ── numbers ──────────────────────────────────────────────────────────────────
MIN_FONT = 6.5          # no text in any figure below this size
P_FLOOR = 1e-300


def nlp(p):
    """-log10 P with the underflow floor, elementwise."""
    p = np.asarray(p, dtype=float)
    return -np.log10(np.clip(p, P_FLOOR, 1.0))


def qq_expected(n):
    """Expected -log10 P for n ranked tests, (i - 0.5) / n."""
    return -np.log10((np.arange(1, n + 1) - 0.5) / n)


def qq_band(n, alpha=0.05):
    """(expected, lower, upper) -log10 P of the pointwise Beta(i, n-i+1) band."""
    from scipy.stats import beta
    i = np.arange(1, n + 1)
    lo = -np.log10(beta.ppf(1 - alpha / 2, i, n - i + 1))
    hi = -np.log10(beta.ppf(alpha / 2, i, n - i + 1))
    return qq_expected(n), lo, hi


def text_extent_in(fig, text, fontsize, rotation=0.0, **kw):
    """(width, height) in inches of `text` as it would render on `fig`."""
    t = fig.text(0, 0, text, fontsize=fontsize, rotation=rotation, **kw)
    fig.canvas.draw()
    bb = t.get_window_extent(fig.canvas.get_renderer())
    t.remove()
    return bb.width / fig.dpi, bb.height / fig.dpi


def fmt_mb(span_bp):
    """A tick formatter for genomic positions in Mb with decimals suited to the span."""
    from matplotlib.ticker import FuncFormatter
    d = max(2, int(math.ceil(-math.log10(max(span_bp, 1) / 1e6))) + 1)
    d = min(d, 5)
    return FuncFormatter(lambda v, _pos: f'{v / 1e6:.{d}f}')
