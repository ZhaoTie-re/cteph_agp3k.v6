#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Write the companion .md that sits beside every figure this component
#           produces. The caption inside a figure has to stay short enough not to
#           crowd the panels; the sidecar is where the full explanation lives —
#           what each panel plots, how every statistic is defined, what the
#           numbers actually came out as for THIS cohort/peak, how to read the
#           figure, and what it cannot answer.
#
#           One writer for all figures, and ALWAYS the same eight sections in the
#           same order, so a reader who has read one sidecar knows where to look
#           in every other one:
#
#             1  The question this figure answers
#             2  Panels
#             3  Interpretation             (moved out of the figure caption)
#             4  Values in this rendering
#             5  Full statistics            (optional; rendered tables)
#             6  How to read it
#             7  What this figure does NOT establish
#             8  Symbols + Model            (moved out of the figure caption)
#
#           Sections 3 and 8 used to be printed inside the PNG. The caption grew
#           to 63-70 % of the canvas, which is the wrong shape for a figure being
#           laid out as a slide, so they live here instead. Nothing was dropped.
# Component: assoc_plink2
# ---------------------------------------------------------------------------
from pathlib import Path

import plot_style as S


def _fmt(v):
    if v is None:
        return '—'
    if isinstance(v, float):
        if v != v:                      # NaN
            return '—'
        if v and (abs(v) < 1e-3 or abs(v) >= 1e5):
            return f'{v:.3e}'
        return f'{v:,.4g}'
    if isinstance(v, int):
        return f'{v:,}'
    return str(v)


def _plain(t):
    """Figure mathtext -> readable prose. The sidecar is markdown, not a figure."""
    return (t.replace('$\\lambda_{\\mathrm{GC}}$', 'lambda_GC')
             .replace('$\\mathrm{OR}$', 'OR')
             .replace('$N_{\\mathrm{eff}}$', 'N_eff')
             .replace('$\\mathrm{EAF}$', 'EAF')
             .replace('$r^{2}$', 'r^2')
             .replace('\\times', 'x').replace('\\geq', '>=').replace('\\leq', '<=')
             .replace('\\mathrm', '').replace('\\lambda', 'lambda').replace('\\chi', 'chi')
             .replace('$', '').replace('{', '').replace('}', '').replace('\\', ''))


def _table(df, max_rows=None):
    """A DataFrame as a GitHub markdown table, values formatted like the rest."""
    d = df if max_rows is None else df.head(max_rows)
    cols = list(d.columns)
    out = ['| ' + ' | '.join(str(c) for c in cols) + ' |',
           '|' + '|'.join('---' for _ in cols) + '|']
    for _, r in d.iterrows():
        out.append('| ' + ' | '.join(_fmt(r[c]) for c in cols) + ' |')
    if max_rows is not None and len(df) > max_rows:
        out.append(f'| … | {len(df) - max_rows} further rows in the TSV | ' +
                   ' | '.join('' for _ in cols[2:]) + ' |')
    return out


def write_doc(png_path, *, title, question, panels, interpretation=None, numbers=None,
              tables=None, reading=None, limits=None, defs=None, model=None,
              methods_ref='../../../docs/METHODS.md'):
    """Write `<png>.md` next to `<png>`.

    Parameters
    ----------
    title    : one line naming the figure and its subject
    question : the single question this figure exists to answer
    panels   : list of (panel_letter, heading, body) — body explains what is
               plotted AND how the statistic is computed
    numbers  : dict or list of (label, value) — the concrete values behind this
               particular rendering, so the file is checkable against the figure
    tables   : list of (caption, DataFrame[, max_rows]) — full statistics that a
               figure cannot hold. This is where a per-cohort statistics block
               goes, rather than being crushed into a panel.
    interpretation : the objective discussion that used to sit in the figure
               caption — what the panels mean together, and what they do not
    reading  : list of strings — how to read it, in order
    limits   : list of strings — what the figure cannot establish
    defs     : list of SYMBOL_DEFS keys to expand into a glossary
    model    : the estimator / equation, in mathtext as the figure used to print it
    """
    png = Path(png_path)
    out = png.with_suffix('.md')
    L = [f'# {title}', '', f'**Figure file:** `{png.name}`', '',
         '## The question this figure answers', '', question, '']

    L += ['## Panels', '']
    for letter, heading, body in panels:
        L += [f'**({letter}) {heading}**', '', body, '']

    if interpretation:
        L += ['## Interpretation', '', _plain(interpretation), '']

    if numbers:
        items = numbers.items() if isinstance(numbers, dict) else numbers
        L += ['## Values in this rendering', '',
              '| quantity | value |', '|---|---|']
        L += [f'| {k} | {_fmt(v)} |' for k, v in items]
        L += ['']

    if tables:
        L += ['## Full statistics', '']
        for entry in tables:
            cap, df = entry[0], entry[1]
            max_rows = entry[2] if len(entry) > 2 else None
            if df is None or not len(df):
                continue
            L += [f'**{cap}**', '']
            L += _table(df, max_rows)
            L += ['']

    if reading:
        L += ['## How to read it', '']
        L += [f'{i}. {s}' for i, s in enumerate(reading, 1)]
        L += ['']

    if limits:
        L += ['## What this figure does *not* establish', '']
        L += [f'- {s}' for s in limits]
        L += ['']

    if defs:
        L += ['## Symbols', '']
        for k in defs:
            d = S.SYMBOL_DEFS.get(k)
            if d:
                lab, txt = d
                L += [f'- **{_plain(lab)}** — {_plain(txt)}', '']

    if model:
        L += ['## Model', '', '```', _plain(model), '```', '']

    L += ['---', '', f'Methods and rationale: [`METHODS.md`]({methods_ref})', '']
    out.write_text('\n'.join(L))
    return out
