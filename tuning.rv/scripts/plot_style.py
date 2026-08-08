#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Single source of truth for the tuning.rv figure aesthetic. One
#           canonical publication / slide style (Arial, 600 dpi, despined,
#           Wong 2011 colorblind-safe palette, italic math symbols), plus the
#           shared helpers every plotter uses (despine, panel letters, group /
#           target palettes, compact KPI metric strip). Imported as a sibling
#           module by the other scripts (ScriptDir is on sys.path[0] when
#           Nextflow runs `python3 ${ScriptDir}/x.py`), so no packaging needed.
# Project : cteph_agp3k.v6 / tuning.rv  (rare-variant depth-confounding QC)
# Used by : detect_outlier_samples.py, plot_rate_distribution.py,
#           detect_depthdiff_variants.py, plot_sweep.py, compare_cohorts.py
# ---------------------------------------------------------------------------
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D

# ── Palette (Wong 2011, colorblind-safe) ─────────────────────────────────────
GROUP_COLORS = {'Case': '#D55E00', 'PH': '#D55E00', 'Control': '#0072B2', 'NaN': '#7F7F7F'}
GROUP_FALLBACK = ['#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#F0E442']
# Target-depth is a 2-level variable (15x / 30x); distinct HUES read far better
# than two shades of blue (and stay legible at low scatter alpha). Wong-safe:
# blue (15x), amber (30x), then teal / pink as fallbacks for any extra levels.
TARGETDP_COLORS = ['#0072B2', '#E69F00', '#009E73', '#CC79A7']

# Accent colors used by the single-series / trend figures.
COUNT_COLOR = '#333333'
RHO_COLOR = '#0072B2'     # blue      — depth-burden correlation
BETA_COLOR = '#D55E00'    # vermilion — no-PC primary (depth-adjusted) series
GROUP_COLOR = '#009E73'   # green     — +ancestry-PC sensitivity series
BAR_COLOR = '#3B7DB3'
LINE_COLOR = '#1F4E79'
KNEE_COLOR = '#2E8B57'

# ── Academic math symbols (single source; identical in figures + docs) ───────
# Per-sample burden rate, its regression terms, the sweep threshold, and the
# per-variant depth-stratified missingness gap. Keep these consistent with the
# "Symbols & notation" table in docs/METHODS.md.
S_MINAC   = r'$\mathit{S}_{minAC}$'
D_MEAN    = r'mean depth'                          # readable axis phrase
RATE      = r'$\mathit{r}$'                        # per-sample minor-allele burden rate (per 1,000 callable)
BETA_G    = r'$\beta_{\mathrm{group}}$'            # OLS Case-vs-Control coefficient (word subscript = readable)
BETA_D    = r'$\beta_{\mathrm{depth}}$'            # OLS depth coefficient (per +1 SD depth)
VMISS_SIGNED = r'missing-rate difference  (30$\times$ $-$ 15$\times$)'
VMISS_ABS    = r'|missing-rate difference|   ($\Delta$)'

# ── What each symbol MEANS IN THIS STUDY (not a generic glossary) ────────────
# Rendered into figure captions by caption_block(defs=[...]). Every entry states
# the quantity, its units where applicable, and its role in the tuning decision,
# because the design is fully confounded (cases 30x, controls 15x, no platform
# overlap) and a symbol read generically would be read wrongly here.
# Keep in lock-step with the "Symbols & notation" table in docs/METHODS.md.
SYMBOL_DEFS = {
    'rate': (r'$r_i$, per-sample minor-allele burden rate '
             r'$(n_{\mathrm{het}}+2n_{\mathrm{hom}})/n_{\mathrm{call}}\times10^{3}$, i.e. variants per 1,000 '
             r'callable sites — the response variable throughout.'),
    'minac': (r'$\kappa$ (minAC), the minimum minor-allele count a variant must reach to be kept before $r$ is '
              r'recomputed — the sweep axis, and the quantity being tuned.'),
    'beta_group': (r'$\beta_{\mathrm{group}}$, OLS coefficient on the Case indicator = the APPARENT '
                   r'Case$-$Control difference in $r$. Cases are 30$\times$ and controls 15$\times$ with no '
                   r'platform overlap, so at low $\kappa$ it measures the depth/platform artifact, not phenotype; '
                   r'it is the DECISION axis and $\kappa$ is taken at its non-significant trough.'),
    'beta_depth': (r'$\beta_{\mathrm{depth}}$, OLS coefficient on standardized mean depth $z(D)$ = residual depth '
                   r'sensitivity of $r$ given group (and PCs). DIAGNOSTIC only: its minimum sits at a different '
                   r'$\kappa$ from that of $\beta_{\mathrm{group}}$ — well inside the range where the apparent '
                   r'case/control effect is already large — so it must never drive $\kappa$.'),
    'ci_sig': (r'significance, non-significant requires BOTH $p\geq\alpha$ and $0\in[L,U]$ (95% CI) — the same '
               r'rule that selects the recommended $\kappa$.'),
    'models': (r'models, no-PC ($r\sim\mathrm{group}+z(D)$) is PRIMARY: here $\beta_{\mathrm{group}}$ is a '
               r'DETECTOR for the artifact, and with zero platform overlap the artifact is the case/control '
               r'contrast, so a covariate carrying that contrast hides part of what the detector is for. The '
               r'ancestry PCs carry some of it, and adjusting on them moves the reading by more than the reading '
               r'itself, in cohort-dependent directions; +PC is therefore a SENSITIVITY check only. The '
               r'downstream association test keeps its PCs — it needs an estimator, not a detector.'),
    'delta_v': (r'$\Delta_v=|V^{30\times}_v-V^{15\times}_v|$, per-variant gap between missing rates computed '
                r'separately in the 30$\times$ and 15$\times$ samples — how differentially a variant drops out '
                r'with depth, the artifact the burden test would otherwise absorb.'),
    'knee': (r'knee, the Kneedle elbow of the CDF $F(\Delta)$ above which $\Delta$ becomes a heavy tail; variants '
             r'beyond it are excluded (protected candidates exempt).'),
    'robustz': (r'robust-$Z_i=0.6745(\varepsilon_i-\mathrm{median}\,\varepsilon)/\mathrm{MAD}(\varepsilon)$ on the '
                r'Theil$-$Sen depth residual $\varepsilon_i=r_i-\hat{r}(D_i)$ — median/MAD based, so extreme '
                r'samples cannot inflate their own threshold.'),
    'lambda': (r'$\lambda$, genome-wide inflation factor of the downstream rare-variant test — the final '
               r'calibration check: this tuning bounds the artifact, $\lambda$ confirms it.'),
    'variants_retained': (r'variants retained, call-set size surviving $\kappa$ — the power cost of raising it.'),
}

# ── Formatters ────────────────────────────────────────────────────────────────
INT_FMT = FuncFormatter(lambda x, _: format(int(x), ','))


def _sci_fmt(v, _):
    """Scientific axis labels in superscript form (e.g. 1.6×10⁷), publication-style.

    Small magnitudes (|v| < 1000) are shown as plain integers; 0 stays '0'.
    """
    import math
    if v == 0:
        return '0'
    if abs(v) < 1000:
        return f'{v:g}'
    e = int(math.floor(math.log10(abs(v))))
    m = v / 10 ** e
    ms = f'{m:.1f}'.rstrip('0').rstrip('.')      # 1-decimal mantissa, drop trailing .0
    return rf'${ms}\times10^{{{e}}}$'


# Superscript scientific formatter for count / large-magnitude axes.
# (Name kept as COUNT_FMT so existing consumers pick it up unchanged.)
COUNT_FMT = FuncFormatter(_sci_fmt)

# ── Type scales: one definition, two sizes ───────────────────────────────────
#   slide : larger type so a figure stays legible when projected (PPT / talk).
#   paper : tighter type for a dense multi-panel manuscript figure.
_SCALES = {
    'slide': dict(base=13, title=15, label=14, tick=12, legend=11, panel=16),
    'paper': dict(base=11, title=12.5, label=12, tick=10.5, legend=10, panel=14),
}


def setup_style(scale='slide'):
    """Apply the canonical rcParams at the requested type scale; return the scale dict."""
    s = _SCALES[scale]
    plt.rcParams.update({
        'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
        'savefig.dpi': 600, 'savefig.bbox': 'tight', 'figure.dpi': 150,
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'mathtext.fontset': 'dejavusans',
        'font.size': s['base'], 'axes.titlesize': s['title'], 'axes.labelsize': s['label'],
        'xtick.labelsize': s['tick'], 'ytick.labelsize': s['tick'], 'legend.fontsize': s['legend'],
        'axes.linewidth': 1.0, 'axes.edgecolor': 'black',
        'text.color': 'black', 'axes.labelcolor': 'black',
        'xtick.color': 'black', 'ytick.color': 'black',
        'xtick.direction': 'out', 'ytick.direction': 'out',
        'axes.grid': True, 'grid.color': '#DDDDDD', 'grid.linewidth': 0.8, 'grid.alpha': 1.0,
        'axes.spines.top': False, 'axes.spines.right': False,
    })
    return s


# ── Axis / annotation helpers ─────────────────────────────────────────────────
def despine(ax, grid_axis='both'):
    """Recessive grid, top/right spines off, outward ticks — the shared axis look."""
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if grid_axis in ('both', 'x', 'y'):
        ax.grid(True, which='major', axis=grid_axis)
    else:
        ax.grid(False)
    ax.tick_params(direction='out', length=4.5, width=1.0)


def panel_letters(axes, letters='abcdefgh', size=None, dx=-0.07, dy=1.04, ha='right', inside=False):
    """Bold (a)/(b)/… panel tags in one consistent position/size across all figures.

    inside=True places the tag in the panel's empty upper-left corner (use when an
    outside tag would collide with a wide/2-line y-axis label or the suptitle).
    """
    size = size or _SCALES['slide']['panel']
    for ax, letter in zip(axes, letters):
        if inside:
            ax.text(0.02, 0.97, f'({letter})', transform=ax.transAxes, fontsize=size,
                    fontweight='bold', ha='left', va='top', zorder=25)
        else:
            ax.text(dx, dy, f'({letter})', transform=ax.transAxes,
                    fontsize=size, fontweight='bold', ha=ha, va='bottom')


def panel_tag(ax, letter, title=None, *, size=None, tsize=None, pad=7):
    """Nature/Science panel label: a **bold letter** flush-left in the title slot,
    with an optional short centred title. Both live ABOVE the axes box (matplotlib
    reserves the space), so neither can ever overlap the data, ticks, or y-label.

    The bold left letter and the centred title are two independent title artists
    (loc='left' / loc='center'); keep the centred title short so its left edge
    stays clear of the letter on narrow panels.
    """
    size = size or (_SCALES['slide']['panel'] - 1)
    ax.set_title(f'{letter}', loc='left', fontweight='bold', fontsize=size, pad=pad)
    if title:
        ax.set_title(title, loc='center', fontsize=tsize or plt.rcParams['axes.titlesize'], pad=pad)


# General categorical palette (Wong 8-safe) for arbitrary discrete levels (e.g. platform).
CATEGORICAL = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#F0E442', '#999999']


def categorical_palette(levels):
    """Map discrete levels (first-seen order) to distinct colorblind-safe colors."""
    lv = list(dict.fromkeys(str(x) for x in levels))
    return {l: CATEGORICAL[i % len(CATEGORICAL)] for i, l in enumerate(lv)}


def group_palette(levels):
    """Map phenotype-group levels to fixed Wong colors (unknowns cycle the fallback)."""
    pal, fb = {}, 0
    for g in levels:
        g = str(g)
        if g in GROUP_COLORS:
            pal[g] = GROUP_COLORS[g]
        else:
            pal[g] = GROUP_FALLBACK[fb % len(GROUP_FALLBACK)]
            fb += 1
    return pal


def target_palette(levels):
    """Numeric-aware sequential palette for TargetDP levels (NaN -> grey); returns (pal, ordered)."""
    s = pd.Series([str(x) for x in levels])
    num = s.str.extract(r'(\d+(?:\.\d+)?)', expand=False).astype(float)
    order = pd.DataFrame({'lvl': [str(x) for x in levels], 'num': num})
    order['_na'] = order['num'].isna()
    ordered = order.sort_values(['_na', 'num', 'lvl'])['lvl'].tolist()
    pal = {}
    for i, lvl in enumerate(ordered):
        pal[lvl] = '#7F7F7F' if lvl.lower() == 'nan' else TARGETDP_COLORS[i % len(TARGETDP_COLORS)]
    return pal, ordered


def metric_strip(ax, items, title=None, value_color='#111111', label_color='#666666'):
    """Render a horizontal KPI row on a dedicated (wide, short) axis.

    items : list of (label, value) — value is drawn bold above a muted label.
    """
    ax.axis('off')
    n = max(len(items), 1)
    y_val, y_lbl = (0.60, 0.20) if not title else (0.52, 0.14)
    if title:
        ax.text(0.0, 0.94, title, transform=ax.transAxes, ha='left', va='top',
                fontsize=plt.rcParams['legend.fontsize'], fontweight='bold', color=label_color)
    for i, (label, value) in enumerate(items):
        xc = (i + 0.5) / n
        ax.text(xc, y_val, str(value), transform=ax.transAxes, ha='center', va='center',
                fontsize=plt.rcParams['axes.titlesize'], fontweight='bold', color=value_color)
        ax.text(xc, y_lbl, str(label), transform=ax.transAxes, ha='center', va='center',
                fontsize=plt.rcParams['legend.fontsize'], color=label_color)
        if i > 0:
            ax.axvline(i / n, ymin=0.18, ymax=0.82, color='#E2E2E2', linewidth=1.0)


# ── Shared formula bank (academic notation; single source for figures + docs) ─
# Keep in lock-step with the "Symbols & notation" table in docs/METHODS.md.
FORMULAS = {
    'rate': r'$r_i = (n_{\mathrm{het}} + 2\,n_{\mathrm{hom}})\,/\,n_{\mathrm{call}} \times 10^{3}$      '
            r'($n_{\mathrm{hom}}$ = homozygous-alt; $n_{\mathrm{call}}$ = non-missing genotypes)',
    'ols':  r'$r_i = \beta_0 + \beta_{\mathrm{group}}\,\mathrm{Case}_i + \beta_{\mathrm{depth}}\,z(D_i)'
            r'\;[\,+\,\sum_k \beta_k\,\mathrm{PC}_{k,i}\,] + \varepsilon_i$      '
            r'($\mathrm{Case}_i\!\in\!\{0,1\}$; $z(D)$ = standardized mean depth; $\varepsilon_i\!\sim\!N(0,\sigma^2)$)',
    'betadp': r'$\beta_{\mathrm{depth}}$ = burden-rate change per $+1$ SD mean depth '
              r'(group- [and PC-] adjusted; same OLS as $\beta_{\mathrm{group}}$)',
    'robustz': r'$\varepsilon_i = r_i - \hat{r}(D_i)$ (Theil–Sen fit of $r$ on depth $D$);      '
               r'robust-$Z_i = 0.6745\,(\varepsilon_i - \mathrm{median}\,\varepsilon)\,/\,\mathrm{MAD}(\varepsilon)$',
    'vmiss_diff': r'$\Delta_v = |\,V^{30\times}_v - V^{15\times}_v\,|$ ($V^{d}$ = variant missing rate in stratum $d$);      '
                  r'exclusion threshold = Kneedle knee of the CDF $F(\Delta)$',
}


def legend_above(ax, handles, *, title=None, ncol=None, y=1.02, fontsize=None):
    """Place a legend ABOVE the axes (off the data region), as a horizontal strip.

    Use this instead of an in-panel corner legend so keys never sit over points/
    curves. `title` doubles as the panel caption when the axes has no set_title.
    """
    return ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, y),
                     ncol=ncol or len(handles), frameon=False,
                     fontsize=fontsize or plt.rcParams['legend.fontsize'],
                     handletextpad=0.5, columnspacing=1.4, borderaxespad=0.0,
                     title=title)


def legend_inside(ax, handles, *, loc='upper right', title=None, ncol=1, fontsize=None):
    """Compact, lightly-framed legend INSIDE the axes, placed over an empty corner.

    Preferred when an above-panel strip would be too wide (many discrete levels,
    e.g. sequencing platform) or would collide with the suptitle. The white,
    semi-opaque frame keeps keys legible over the grid; position it over a region
    with no data (e.g. the upper corner of a residual plot).
    """
    fs = fontsize or (plt.rcParams['legend.fontsize'] - 0.5)
    leg = ax.legend(handles=handles, loc=loc, ncol=ncol, fontsize=fs, title=title,
                    frameon=True, framealpha=0.92, borderpad=0.55, labelspacing=0.35,
                    handletextpad=0.5, borderaxespad=0.6)
    fr = leg.get_frame()
    fr.set_edgecolor('#C8C8C8'); fr.set_linewidth(0.8); fr.set_facecolor('white')
    if leg.get_title() is not None and leg.get_title().get_text():
        leg.get_title().set_fontweight('bold')
        leg.get_title().set_fontsize(fs + 0.5)
    leg.set_zorder(20)
    return leg


def sig_markers(ax, x, y, p, color, *, alpha=0.05, size=52, zorder=4, ok=None):
    """Filled marker where p < alpha, open (white-filled) marker where it is not.

    THE single definition of "filled = significant" for this project — every sweep
    panel must encode significance the same way, whether or not it also draws a CI
    ribbon. Returns nothing; draws two scatter layers on `ax`.

    ok : optional boolean mask of finite/plottable points (defaults to y.notna()).
    """
    ok = y.notna() if ok is None else ok
    sig = (p.fillna(1.0) < alpha) & ok
    ax.scatter(x[sig], y[sig], color=color, s=size, zorder=zorder,
               edgecolor='white', linewidth=0.8)
    ax.scatter(x[~sig & ok], y[~sig & ok], facecolors='white', edgecolors=color,
               s=size, linewidth=1.6, zorder=zorder)


def sig_legend_handles(alpha=0.05, color='#666666'):
    """The two-key legend entry explaining the filled/open significance convention."""
    return [Line2D([], [], marker='o', color=color, markerfacecolor=color, linestyle='',
                   markersize=8, label=f'significant (p < {alpha:g})'),
            Line2D([], [], marker='o', color=color, markerfacecolor='white', linestyle='',
                   markersize=8, markeredgewidth=1.6, label='not significant')]


def share_y(axes, label=None, *, keep=0):
    """Label the y-axis once for a row of axes that share a scale.

    Repeating an identical (often two-line) y-label on every panel wastes width and
    is what pushes tall labels up into the suptitle. Keeps the label on `axes[keep]`,
    clears it elsewhere, and hides the inner tick labels.
    """
    axes = list(axes)
    for i, ax in enumerate(axes):
        if i == keep:
            if label is not None:
                ax.set_ylabel(label)
        else:
            ax.set_ylabel('')
            ax.tick_params(axis='y', labelleft=False)
    return axes


def sparse_int_ticks(ax, values, *, axis='x', keep=(1, 2, 3, 4, 5), step=5):
    """Label a readable subset of dense integer ticks; unlabelled minor ticks mark the rest.

    Labels every value in `keep` plus multiples of `step` (…5, 10, 15, 20); the
    omitted integers become minor ticks, so a 1..20 axis stays uncrowded without
    hiding any level. Use on the minAC sweep axes when panels are narrow.
    """
    vals = sorted({int(v) for v in values})
    major = [v for v in vals if v in keep or v % step == 0]
    minor = [v for v in vals if v not in major]
    if axis == 'x':
        ax.set_xticks(major); ax.set_xticklabels([str(v) for v in major])
        ax.set_xticks(minor, minor=True)
    else:
        ax.set_yticks(major); ax.set_yticklabels([str(v) for v in major])
        ax.set_yticks(minor, minor=True)
    ax.tick_params(axis=axis, which='minor', length=2.5, width=0.8, color='black')


def _wrap_chars(fig, fs, left, right):
    """Measure how many caption characters actually fit on one figure-wide line.

    A fixed chars-per-inch guess is wrong whenever the type scale or font changes,
    and guessing low silently inflates the caption by extra wrapped lines, which
    then squeezes the panels. Measure the real advance width with the renderer and
    keep a margin for the mathtext runs (β, subscripts) that are wider than plain
    text.
    """
    probe = ('OLS beta and the depth-adjusted Case-Control burden difference across '
             'the full minAC range, then a steep climb.')
    try:
        r = fig.canvas.get_renderer()
        t = fig.text(0, 0, probe, fontsize=fs)
        per_char = (t.get_window_extent(renderer=r).width / fig.dpi) / len(probe)
        t.remove()
        usable = fig.get_size_inches()[0] * (right - 0.012)   # canvas x=0.012 → right margin
        return max(80, int(0.93 * usable / per_char))         # 0.93: slack for mathtext runs
    except Exception:
        return max(80, int(fig.get_size_inches()[0] * 12.5))


def caption_block(fig, *, title=None, panels=None, interpret=None, defs=None, notes=None, model=None,
                  letters='abcdefgh', wrap=None, top=0.90, left=0.07, right=0.985,
                  wspace=0.26, hspace=0.42, extra_bottom=0.0):
    """Journal (Nature/Science/Elsevier)-style caption rendered BELOW the axes.

    A thin rule separates the panels from the caption. Reading order (top → bottom):
        <title>          bold declarative sentence — the figure's finding
        a, <panels[0]>   bold run-in panel letter + what that panel shows / how computed
        b, <panels[1]>   …one entry per panel, in order
        <interpret>      plain concluding sentence(s) — the objective discussion (探讨)
        Definitions: …   what each symbol MEANS IN THIS STUDY (keys into SYMBOL_DEFS)
        Data: <notes>    data source, marker definitions, sample sizes, statistics
        Model: <model>   the estimator / equation (professional notation, mathtext)

    `panels` is a list of strings (one per panel); the leading bold letter is added
    automatically from `letters`. `defs` is a list of SYMBOL_DEFS keys — a figure
    must define every symbol it plots, so the reader never has to infer whether
    β is an effect size or an artifact measure. Deterministic margins (layout engine
    off + subplots_adjust) keep the caption clear of the x-labels and the suptitle
    clear of the panel row. Call once, right before savefig.
    """
    import re
    import textwrap
    fw, fh = fig.get_size_inches()
    fs = max(8, plt.rcParams['legend.fontsize'] - 1)
    if wrap is None:
        wrap = _wrap_chars(fig, fs, left, right)
    def _wrap(txt, width=None):
        """Wrap WITHOUT ever breaking inside a $…$ mathtext span.

        textwrap breaks on whitespace, and mathtext such as $p \\geq \\alpha$ contains
        spaces — a break there renders as literal '$p' / '\\geq \\alpha$' garbage. Mask
        the spaces inside math spans, wrap, then restore them.
        """
        masked = re.sub(r'\$[^$]*\$', lambda m: m.group(0).replace(' ', '\x00'), txt)
        return [ln.replace('\x00', ' ') for ln in textwrap.wrap(masked, width or wrap)]

    def _runin(label, txt, sep, width=None):
        """One bold run-in entry ('a, …' or 'Data: …'), wrapped to figure width."""
        lines = _wrap(f'{label}{sep}{txt}', width) or [f'{label}{sep}'.rstrip()]
        lines[0] = f'$\\mathbf{{{label}}}${sep}' + lines[0][len(label) + len(sep):]
        return lines

    def _build_body(width):
        """Caption lines, top → bottom, at a given wrap width."""
        b = []
        for i, p in enumerate(panels or []):
            b.extend(_runin(letters[i], p, ', ', width))   # Nature panel run-in: "a, …"
        if interpret:
            b.extend(_wrap(interpret, width))              # plain concluding discussion
        if defs:
            # Symbols are defined for THIS analysis; keys resolve against SYMBOL_DEFS so
            # a symbol cannot drift between two figures, or between figure and METHODS.
            b.extend(_runin('Definitions', '  '.join(SYMBOL_DEFS.get(k, k) for k in defs), ': ', width))
        if notes:
            b.extend(_runin('Data', notes, ': ', width))
        if model:
            b.append(f'$\\mathbf{{Model}}$: {model}')       # not wrapped (mathtext)
        return b

    def _fs(key, dflt):
        v = plt.rcParams.get(key, dflt)
        return v if isinstance(v, (int, float)) else dflt

    def _height_in(artist):
        """Rendered height of a text artist, in inches — a font metric, so it does
        not change when the figure is resized. That is what makes the two-pass
        layout below stable."""
        try:
            bb = artist.get_window_extent(renderer=fig.canvas.get_renderer())
            return bb.height / fig.dpi
        except Exception:
            return 0.0

    try:
        fig.set_layout_engine('none')           # deterministic margins from here
    except Exception:
        pass

    # ── Pass 1: draw the caption, measure what it actually needs (in inches) ──
    # Estimating multi-line text height from fontsize x linespacing drifts badly as
    # the caption grows, so every gap below is measured, never arithmetic.
    max_w_in = fw * (right - 0.012)          # usable text width, canvas x=0.012 → right

    def _fit(render, fontsize, weight, color, seed):
        """Draw text, then shrink the wrap width until it really fits the canvas.

        A chars-per-inch estimate cannot be trusted across weights and mathtext —
        bold runs ~15 % wider, which is what pushed the bold caption title off the
        right edge. `render(width)` returns the lines at a candidate width; measure,
        re-wrap, repeat. Converges in one or two passes.
        """
        w, t = seed, None
        for _ in range(6):
            if t is not None:
                t.remove()
            t = fig.text(0.012, 0.0, '\n'.join(render(w)), ha='left', va='bottom',
                         fontsize=fontsize, fontweight=weight, color=color, linespacing=1.5)
            try:
                wid = t.get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.dpi
            except Exception:
                return t
            if wid <= max_w_in or w <= 40:
                return t
            w = max(40, int(w * max_w_in / wid) - 1)
        return t

    t_body = t_title = None
    if _build_body(wrap):
        t_body = _fit(_build_body, fs, 'normal', '#333333', wrap)
    if title:
        # The title is set larger AND bold, so it needs its own narrower wrap width —
        # left unwrapped it silently ran off the right edge of the canvas.
        t_title = _fit(lambda w: _wrap(title, w), fs + 1.5, 'bold', '#111111',
                       max(40, int(wrap * fs / (fs + 1.5))))

    pad_in = 0.10                                        # canvas edge → caption text
    gap_in = 0.55 * fs / 72.0                            # body ↔ title
    sep_in = 0.016 * fh                                  # caption ↔ rule ↔ x-labels
    xlab_in = (_fs('axes.labelsize', 14) + _fs('xtick.labelsize', 11) + 8) * 1.5 / 72.0
    h_body = _height_in(t_body) if t_body is not None else 0.0
    h_title = _height_in(t_title) if t_title is not None else 0.0

    cap_in = pad_in + h_body + ((gap_in + h_title) if t_title is not None else 0.0)
    need_in = cap_in + 2 * sep_in + xlab_in + extra_bottom * fh

    # ── Pass 2: if the caption cannot fit in its share, GROW the figure ──────
    # The previous clamp (bottom = min(0.60, …)) silently let the axes overrun a tall
    # caption instead of making room — the caption is content, not decoration.
    MAX_FRAC = 0.50
    if need_in / fh > MAX_FRAC:
        fh = need_in / MAX_FRAC
        fig.set_size_inches(fw, fh, forward=True)
        sep_in = 0.016 * fh

    if t_body is not None:
        t_body.set_position((0.012, pad_in / fh))
    if t_title is not None:
        t_title.set_position((0.012, (pad_in + h_body + gap_in) / fh))

    sep_y = (cap_in + sep_in) / fh                      # thin rule under the panels
    bottom = sep_y + (sep_in + xlab_in + extra_bottom * fh) / fh
    try:
        fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom,
                            wspace=wspace, hspace=hspace)
    except Exception:
        pass
    # thin separator rule delineating the figure from its caption (print-figure style)
    fig.add_artist(Line2D([left, right], [sep_y, sep_y], transform=fig.transFigure,
                          color='#CCCCCC', linewidth=0.8, zorder=0))
