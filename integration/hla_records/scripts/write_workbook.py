#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The workbook hla_records.xlsx, written with XlsxWriter in the house
#           style of the integration workbooks. Presentation only: every number is
#           taken as is from the tables of typing_tables.py / assoc_tables.py.
#           Sheets in the order of the two decks (Typing part = hla_typing.pptx,
#           Association part = hla_assoc.pptx), Contents first (with the map of
#           every slide to its sheet, read from the decks at build time) and
#           Sources and checks last (source files and the result of every check;
#           the only sheet that names files). P values are printed as the slides
#           print them (deckfmt); detail sheets hold numeric P for sorting.
# Inputs  : --typing (typing tables) --assoc (association tables) --decks
#           (decks.json) --slide-map --checks (checks.json) --inputs (inputs.tsv)
# Outputs : --out (xlsx), --meta (workbook.json: per sheet, its blocks and cells)
# Component: hla_records
# ---------------------------------------------------------------------------
import argparse
import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xlsxwriter

sys.path.insert(0, str(Path(__file__).resolve().parent))
import deckfmt as F  # noqa: E402
from hr_common import COHORTS, COH_CODE, COH_LABEL, COMP_LABEL, COMPARISONS, LEAD57  # noqa: E402
from hr_xlsx import (CALL_FILL, FAINT, GREY, HEADER_FILL, INK, INK2, MEDIUM, MISSING, NONE_, NOTE, THIN,  # noqa: E402
                     Table, abbrev, blank, md, missing_block, new_sheet, p_parts, print_setup, sup_parts,
                     title, write_notes)

COLORS = {'narrow_mainland': '#0571B0', 'intermediate_mainland': '#008837', 'full_mainland': '#B2182B'}
TINTS = {'narrow_mainland': '#EAF3F9', 'intermediate_mainland': '#E9F5EE', 'full_mainland': '#FBEDEE'}
NAVY = '#1F3864'
NOTE_W = 150
ABBR = {'BBJ': 'BioBank Japan', 'CI': 'confidence interval', 'CTEPH': 'chronic thromboembolic pulmonary hypertension',
        'DVT': 'deep vein thrombosis', 'FWER': 'family-wise error rate', 'GEC': 'Genetic Type I Error Calculator',
        'GRM': 'genetic relationship matrix', 'GWAS': 'genome-wide association study', 'HLA': 'human leukocyte antigen',
        'IPD-IMGT/HLA': 'the HLA sequence database', 'LD': 'linkage disequilibrium', 'LOCO': 'leave one chromosome out',
        'LRT': 'likelihood-ratio test', 'MAC': 'minor allele count', 'MAF': 'minor allele frequency',
        'Me': 'effective number of independent tests', 'MHC': 'major histocompatibility complex',
        'N / I / F': 'narrow / intermediate / full mainland cohort', 'OR': 'odds ratio', 'PC': 'principal component',
        'SAIGE': 'scalable and accurate implementation of generalized mixed model', 'SPA': 'saddlepoint approximation',
        'ToMMo': 'Tohoku Medical Megabank Organization (61KJPN panel)', 'df': 'degrees of freedom',
        'λGC': 'genomic inflation factor', '1KG JPT': '1000 Genomes Project, Japanese in Tokyo',
        'max-T': 'maximum test statistic over all tests', 'SNP': 'single-nucleotide polymorphism'}


SLN = {}                                       # (deck, slide key) -> slide number, from decks.json


def sl(deck, *keys):
    """'slide 20', 'slides 20–22' or 'slides 2, 36' for slide keys of a deck, so captions follow the decks."""
    ns = sorted(SLN[(deck, k)] for k in keys)
    runs_, a = [], ns[0]
    for x, y in zip(ns, ns[1:] + [None]):
        if y != x + 1:
            runs_.append(f'{a}–{x}' if x > a else f'{a}')
            a = y
    return ('slide ' if len(ns) == 1 else 'slides ') + ', '.join(runs_)


def sa(*keys):
    return sl('association', *keys)


def sty(*keys):
    return sl('typing', *keys)


def ab(*keys):
    return abbrev({k: ABBR[k] for k in keys})


# ── cell writers ────────────────────────────────────────────────────────────
def marker_parts(m):
    m = str(m)
    if '*' in m or ':' not in m:
        return [m]
    p = m.split(':')
    p[1] = p[1].replace('-', '−')
    return [('i', p[0]), ' ' + ' '.join(p[1:])]


def mlab(m):
    """A marker in light markup: 'DQB1:57:V' -> '*DQB1* 57 V' (typographic minus); alleles and SNPs unchanged."""
    m = str(m)
    if '*' in m or ':' not in m or m.startswith('chr'):
        return m
    p = m.split(':')
    return f"*{p[0]}* " + ' '.join(p[1:]).replace('-', '−')


def mlist(ms, sep=', '):
    return sep.join(mlab(m) for m in str(ms).replace(';', ',').split(',') if m.strip())


def ci_text(o, l_, u):
    if blank(o):
        return None
    return f'{F.sig3(o)} ({F.sig3(l_)}–{F.sig3(u)})'


NA = '\x00na'                       # not applicable: written as an en dash, never as missing
NUMFMT = {'int': '#,##0', 'f1': '0.0;"−"0.0', 'f2': '0.00;"−"0.00', 'f3': '0.000;"−"0.000', 'f4': '0.0000',
          'pct0': '0%', 'pct1': '0.0%;"−"0.0%', 'pct2': '0.00%', 'pnum': '0.00E+00', 'or': '0.00'}


def C(key, leaf, kind='text', w=10, blk=None, sub=None, note=None, coh=None, align=None, val=None, hi=None, vline=None,
      span=1, wrap=False):
    return dict(key=key, leaf=leaf, kind=kind, w=w, blk=blk, sub=sub, note=note, coh=coh, align=align, val=val, hi=hi,
                vline=vline, span=span, wrap=wrap)


def lines_needed(text, width):
    """Lines a header text takes when wrapped at a column width (Arial 9 bold: about 1.1 characters per unit)."""
    cap, n, cur = max(int(width * 1.1) - 1, 4), 1, 0
    for w in str(text).split():
        if cur and cur + 1 + len(w) > cap:
            n, cur = n + 1, len(w)
        else:
            cur += (1 if cur else 0) + len(w)
    return n


def plain_text(x):
    return x if isinstance(x, str) else ''.join(p if isinstance(p, str) else p[1] for p in x)


def put(sh, caption, cols, rows, span=None, group=None, bold=None, freeze=None, autofilter=False, key=None):
    """One booktabs block: a caption line, three header rows, the body. rows: list of dicts. A column with span k
    takes k sheet columns (merged), for long text, so that it does not widen the columns of the other blocks."""
    r = sh.row
    if caption:
        sh.rich(r, 0, md(caption), {'align': 'left'}, b=True, color=NAVY, sz=10)
        r += 1
    phys, start = [], []
    for c in cols:
        start.append(len(phys))
        k = c['span']
        if k > 1 and c['blk'] is not None:
            raise SystemExit(f'{sh.name} / {caption}: a spanned column cannot sit under a spanner')
        phys += [c] + [dict(c, cont=True)] * (k - 1)
    spec = [dict(blk=c['blk'], sub=c['sub'], leaf=c['leaf'] if isinstance(c['leaf'], list) else md(c['leaf'])
                 if '*' in c['leaf'] else c['leaf'], note=c['note'], w=c['w'], coh=c['coh'], vline=c['vline'],
                 cont=c.get('cont', False)) for c in phys]
    # header heights: the leaf row grows until every leaf fits; spanner rows grow for long spanner texts
    need = 12.75
    for c in cols:
        lines = lines_needed(plain_text(c['leaf']) + (' x' if c['note'] else ''), c['w'] * c['span'])
        room = 20 + 12.75 if c['blk'] is None else (12.75 if c['sub'] is None else 0)
        need = max(need, lines * 11.5 + 3 - room)
    T = Table(sh, spec, r, span or {}, COLORS, TINTS, header_h=need, set_width=False)
    for lvl, field, base in ((0, 'blk', 20), (1, 'sub', 12.75)):
        h = base
        for g in {c[field] for c in cols if c[field]}:
            w = sum(c['w'] for c in cols if c[field] == g)
            h = max(h, lines_needed(plain_text(span[g][0]) if (lvl == 0 and span and g in span) else g, w) * 12 + 4)
        if h > base:
            sh.ws.set_row(r + lvl, h)
    for j, c in enumerate(phys):
        if c['span'] == 1:
            sh.widths[j] = max(sh.widths.get(j, 0), c['w'])
    body0 = r + 3
    for c in cols:   # a column without a single value is a wrong key, never data
        if rows and all(blank(c['val'](row) if c['val'] else row.get(c['key'])) for row in rows):
            raise SystemExit(f"{sh.name} / {caption}: column '{c['key']}' is empty in every row")
    cells = []
    for i, row in enumerate(rows):
        rr = body0 + i
        gtop = bool(group and i > 0 and row.get(group) != rows[i - 1].get(group))
        last = i == len(rows) - 1
        isbold = bool(bold and bold(row))
        nl = 1
        seen = []
        for jl, c in enumerate(cols):
            j = start[jl]
            v = c['val'](row) if c['val'] else row.get(c['key'])
            d = T.rules(j, {'align': c['align'] or ('left' if c['kind'] in ('text', 'gene', 'marker', 'md', 'comp',
                                                                                'cohort') else 'center')},
                        group_top=gtop, last=last)
            if c['hi'] and c['hi'](row):
                d.update({'bg_color': CALL_FILL, 'pattern': 1})
            if c['wrap']:
                d.update({'text_wrap': True, 'valign': 'top'})
                if isinstance(v, str) and v != NA:
                    nl = max(nl, lines_needed(v, c['w'] * c['span'] * 1.15))
            if c['span'] > 1:
                sh.ws.merge_range(rr, j, rr, j + c['span'] - 1, '', sh.fmt(**d))
            seen.append(write_cell(sh, rr, j, c['kind'], v, d, isbold))
        cells.append(seen)
        if nl > 1:
            sh.ws.set_row(rr, 12 * nl + 1)
    sh.blocks.append({'caption': caption, 'key': key or caption, 'header_row': r, 'first_row': body0,
                      'n_rows': len(rows), 'cols': [c['key'] for c in cols], 'kinds': [c['kind'] for c in cols],
                      'pcols': start, 'row_keys': [row.get('_key') for row in rows], 'cells': cells,
                      'raw': {c['key']: [row.get(c['key']) for row in rows] for c in cols
                              if c['kind'] in NUMFMT or c['kind'] == 'p'} if sh.name.endswith('(detail)') else {}})
    if freeze:
        sh.ws.freeze_panes(body0, freeze)
    if autofilter and rows:
        sh.ws.autofilter(r + 2, 0, body0 + len(rows) - 1, len(phys) - 1)
    sh.ncol = max(sh.ncol, len(phys))
    sh.row = body0 + len(rows) + 2
    return body0


SUP = str.maketrans('−0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')


def shown(parts):
    """The text a cell shows, superscript runs as Unicode superscripts (what the audit reads back)."""
    return ''.join(p if isinstance(p, str) else (p[1].translate(SUP) if p[0] == 'sup' else p[1]) for p in parts)


def write_cell(sh, r, c, kind, v, d, bold=False):
    """Write one body cell; return what it shows (a number, or its text) for workbook_cells.json.gz."""
    if isinstance(v, float) and not np.isfinite(v):
        v = None
    if isinstance(v, str) and v == NA:
        sh.rich(r, c, [NONE_], dict(d, align='center'), color=GREY)
        return NONE_
    if v is None or (isinstance(v, str) and v == ''):
        sh.rich(r, c, [MISSING], dict(d, align='center'), color=FAINT, b=bold)
        return MISSING
    if kind in NUMFMT:
        f = sh.fmt(**dict(d, num_format=NUMFMT[kind], bold=True if bold else None))
        sh.ws.write_number(r, c, float(v), f)
        return float(v)
    if kind == 'p':
        parts = p_parts(float(v))
    elif kind == 'gene':
        parts = [('i', str(v))]
    elif kind == 'marker':
        parts = marker_parts(v)
    elif kind == 'md':
        parts = md(str(v))
    elif kind == 'bool':
        parts = ['Yes' if v in (True, 'True', 'true', 1, '1') else 'No']
    elif kind == 'comp':
        parts = [COMP_LABEL.get(v, str(v))]
    elif kind == 'cohort':
        parts = [(COLORS[v], COH_LABEL[v])] if v in COLORS else [str(v)]
    else:
        parts = sup_parts(str(v)) if any(ch in str(v) for ch in '⁻⁰¹²³⁴⁵⁶⁷⁸⁹') else [str(v)]
    sh.rich(r, c, parts, d, b=bold)
    return shown(parts)


def coh_span(extra=None):
    s = {COH_LABEL[c]: ([(COLORS[c], COH_LABEL[c])], COLORS[c], c) for c in COHORTS}
    s.update(extra or {})
    return s


def finish(sh, notes, footer):
    row = sh.row
    row, lens = write_notes(sh, row, notes, NOTE_W)
    for j in range(sh.ncol):
        sh.ws.set_column(j, j, sh.widths.get(j, 12))
    print_setup(sh.ws, row, max(sh.ncol, 1), footer)
    sh.last_row = row
    return lens


def cohort_cols(fields, sep=True):
    """[(key, leaf, kind, w)] -> the same columns for each cohort, under its coloured spanner."""
    cols = []
    for c in COHORTS:
        for k, leaf, kind, w in fields:
            cols.append(C(f'{k}_{c}', leaf, kind, w, blk=COH_LABEL[c], coh=c))
    return cols


def kv(v):
    """A value of a key-value block: integers with thousands separators, also inside text (not in allele names
    or release numbers)."""
    t = str(v)
    if re.fullmatch(r'\d+', t):
        return f'{int(t):,}'
    return re.sub(r'(?<![\d.*:])(\d{4,})(?![\d.:])', lambda m: f'{int(m.group(1)):,}', t)


def records(df):
    return [dict((k, (None if (isinstance(v, float) and not np.isfinite(v)) else v)) for k, v in r.items())
            for r in df.to_dict('records')]


# ── the book ────────────────────────────────────────────────────────────────
class Book:
    def __init__(self, a):
        self.a = a
        self.wb = xlsxwriter.Workbook(a.out, {'nan_inf_to_errors': False})
        self.wb.set_properties({'title': 'HLA typing and HLA association in CTEPH: record of the analyses',
                                'author': 'CTEPH WGS project'})
        self.sheets = []
        self.T = Path(a.typing)
        self.A = Path(a.assoc)
        self.tf = json.loads((self.T / 'typing_facts.json').read_text())
        self.af = json.loads((self.A / 'assoc_facts.json').read_text())
        self.decks = json.loads(Path(a.decks).read_text())
        SLN.update({(d, x['key']): x['number'] for d, v in self.decks.items() for x in v['slides']})
        self.smap = pd.read_csv(a.slide_map, sep='\t', dtype=str).fillna('')
        parts = [json.loads(Path(f).read_text()) for f in a.checks]
        self.checks = {'all_passed': all(x['all_passed'] for x in parts),
                       'checks': [c for x in parts for c in x['checks']],
                       'deck_files': [d for x in parts for d in x.get('deck_files', [])]}

    def t(self, name):
        return pd.read_csv(self.T / name, sep='\t', dtype={'position': str, 'residue': str})

    def s(self, name):
        return pd.read_csv(self.A / name, sep='\t', dtype={'position': str, 'residue': str, 'marker': str})

    def sheet(self, name, ttl, desc, part):
        sh = new_sheet(self.wb, name)
        title(sh, ttl, desc)
        sh.row, sh.widths, sh.ncol, sh.part, sh.ttl = 3, {}, 1, part, ttl
        self.sheets.append(sh)
        return sh

    def slides_of(self, sheet):
        out = {'typing': [], 'association': []}
        for r in self.smap.itertuples():
            if sheet in r.sheets.split(';'):
                n = next(s['number'] for s in self.decks[r.deck]['slides'] if s['key'] == r.key)
                out[r.deck].append(n)
        return out


# ── Contents ────────────────────────────────────────────────────────────────
SHEET_INFO = {}


def contents(B, order):
    sh = next(s for s in B.sheets if s.name == 'Contents')
    rows = []
    for name in order:
        if name in ('Contents', 'Sources and checks'):
            continue
        sl = B.slides_of(name)
        rows.append({'sheet': name, 'part': 'Typing' if name.startswith('T') else 'Association',
                     'typing': ', '.join(map(str, sorted(sl['typing']))) or NA,
                     'association': ', '.join(map(str, sorted(sl['association']))) or NA,
                     'per': SHEET_INFO[name][0], 'what': SHEET_INFO[name][1]})
    cols = [C('sheet', 'Sheet', 'text', 26), C('part', 'Part', 'text', 11),
            C('typing', 'Typing deck', 'text', 20, blk='Slides', note='a'),
            C('association', 'Association deck', 'text', 24, blk='Slides', note='a'),
            C('per', 'One row per', 'text', 28), C('what', 'What it records', 'text', 12, span=8)]
    body0 = put(sh, 'The sheets', cols, rows, span={'Slides': ([('b', 'Slides')], INK, None)}, group='part')
    for i, r in enumerate(rows):
        sh.ws.write_url(body0 + i, 0, f"internal:'{r['sheet']}'!A1", string=r['sheet'],
                        cell_format=sh.fmt(align='left', font_color='#1F4E79', underline=1))
    # every slide, with its sheet or the reason it is not tabulated
    rows = []
    for deck in ('typing', 'association'):
        for s in B.decks[deck]['slides']:
            m = B.smap[(B.smap.deck == deck) & (B.smap.key == s['key'])]
            if m.empty:
                raise SystemExit(f'slide map: no entry for {deck} slide {s["number"]} ({s["key"]})')
            m = m.iloc[0]
            rows.append({'deck': B.decks[deck]['name'], 'number': s['number'], 'section': s['section'],
                         'title': s['title'] or 'Cover', 'sheets': m.sheets.replace(';', '; ') if m.sheets != '-' else NA,
                         'reason': m.reason or NA, '_key': f'{deck}|{s["key"]}'})
    cols = [C('deck', 'Deck', 'text', 26), C('number', 'Slide', 'int', 11), C('section', 'Section', 'text', 20),
            C('title', 'Title', 'text', 26, span=3), C('sheets', 'Recorded in', 'text', 12, span=3),
            C('reason', 'Note', 'text', 12, span=6, wrap=True)]
    put(sh, 'Every slide of the two decks', cols, rows, group='deck')
    tf, af = B.tf, B.af
    notes = [('a', ['Slides, the slide numbers of the two decks as they stand at this build: ', ('i', 'HLA typing'),
                    f" ({B.decks['typing']['n']} slides) and ", ('i', 'HLA association'),
                    f" ({B.decks['association']['n']} slides); the second table lists every slide with its sheet, or why "
                    'it has none (–, no slide or no sheet).']),
             ('', [('b', 'Conventions. '), 'Cohorts: narrow, intermediate and full mainland (N / I / F), nested ancestry '
                   'selections of one study; agreement across them is robustness to the ancestry definition, not '
                   'replication. Comparisons: all cases vs. controls (main), with DVT vs. controls, without DVT vs. '
                   'controls, with vs. without DVT among cases (supplementary). Single markers (allele and residue '
                   'tests): random model, LD-adjusted threshold (GEC, primary) and conventional Bonferroni '
                   '(sensitivity). Whole positions: fixed-model omnibus test, Bonferroni over the positions tested. '
                   'The exact max-T (parametric bootstrap) is the benchmark of the thresholds, never a threshold.']),
             ('', [('b', 'Numbers. '), 'P values are printed as on the slides (two significant digits from 0.001, '
                   'three below); counts are counts of tests where a mirror pair (two residue columns that are '
                   'complements in every sample) counts once; every number shared with a slide equals it (sheet '
                   'Sources and checks).']),
             ('', [('b', 'Not included. '), 'Nothing per sample (no sample identifiers anywhere); figures; the '
                   'per-replicate bootstrap values; analyses the decks dropped (Benjamini–Hochberg, gene-level '
                   'omnibus tests, the multivariate-normal max-T).']),
             ab('CTEPH', 'DVT', 'GEC', 'HLA', 'LD', 'N / I / F', 'max-T')]
    return finish(sh, notes, 'Contents')


# ── Typing sheets ───────────────────────────────────────────────────────────
def t1(B):
    sh = B.sheet('T1 Typing run', 'HLA typing: the run and the genomes typed',
                 'HLA-HD on the HLA reads of every genome · counts per sequencing platform and group · '
                 'no sample was removed', 'Typing')
    run = records(B.t('t1_run.tsv'))
    put(sh, 'The typing run', [C('item', 'Item', 'text', 14, span=3), C('value', 'Value', 'text', 8, span=2,
                                                                           align='right')],
        [{**r, 'value': kv(r['value']), '_key': r['deck_key']} for r in run])
    pl = records(B.t('t1_platforms.tsv'))
    for r in pl:
        r['platform_label'] = r['platform'].replace('DNBseq', 'DNBSeq').replace('HiSeqX', 'HiSeq X').replace(' 30x', ' 30×') \
            .replace(' 15x', ' 15×')
    cols = [C('group', 'Group', 'text', 12), C('platform_label', 'Platform', 'text', 20),
            C('samples', 'Genomes', 'int', 8), C('male', 'Male', 'int', 7, blk='Sex'), C('female', 'Female', 'int', 7, blk='Sex'),
            C('read_pairs_median', 'Median', 'int', 9, blk='HLA read pairs', note='a'),
            C('read_pairs_min', 'Minimum', 'int', 9, blk='HLA read pairs', note='a'),
            C('depth_median', 'Median', 'f1', 8, blk='Observed depth', note='b'),
            C('depth_min', 'Minimum', 'f1', 8, blk='Observed depth', note='b'),
            C('depth_max', 'Maximum', 'f1', 8, blk='Observed depth', note='b'),
            C('outliers', 'All', 'int', 7, blk='Sample outliers', note='c'),
            C('out_read_pairs', 'Few read pairs', 'int', 9, blk='Sample outliers', note='c'),
            C('out_no2field', '≥ 3 genes without a genotype', 'int', 12, blk='Sample outliers', note='c'),
            C('out_failed_obligate', 'Gene failed', 'int', 8, blk='Sample outliers', note='c'),
            C('calls', 'Calls', 'int', 8, blk='Calls (genomes × 19 genes)', note='d'),
            C('genotype_share', 'Genotype', 'pct1', 9, blk='Calls (genomes × 19 genes)', note='d',
              hi=lambda r: r['genotype_share'] < 0.97),
            C('high_share', 'High confidence', 'pct1', 10, blk='Calls (genomes × 19 genes)', note='d',
              hi=lambda r: r['high_share'] < 0.92)]
    sp = {k: ([('b', k)], INK, None) for k in ('Sex', 'HLA read pairs', 'Observed depth', 'Sample outliers',
                                               'Calls (genomes × 19 genes)')}
    put(sh, 'Genomes by group and sequencing platform', cols, pl, span=sp, group='group',
        bold=lambda r: r['platform'] in ('All cases', 'All genomes'))
    reasons = B.tf['outlier_reasons']
    g4 = next(r for r in pl if r['group'] == 'Cases' and r['platform'].lower().startswith('dnbseq-g400rs 15'))
    allc = next(r for r in pl if r['platform'] == 'All cases')
    notes = [('a', ['HLA read pairs, read pairs extracted for typing from the HLA region, its alternate contigs, the '
                    'HLA contigs and unmapped reads.']),
             ('b', ['Observed depth, mean genome-wide depth of the genome as sequenced; the platforms are as recorded, '
                    'never merged.']),
             ('c', ['Sample outliers are listed, never removed: HLA read pairs below the 1st percentile of all genomes '
                    f"({reasons.get('read_pairs<p1', 0)}), three or more reported genes without a two-field genotype "
                    f"({reasons.get('no_2field>=3', 0)}), or a gene that every haplotype carries failed "
                    f"({reasons.get('failed_obligate', 0)}); a genome can meet more than one rule, so the reasons add up "
                    'to more than the outliers.']),
             ('d', ['Calls, one per genome and reported gene. Genotype, a two-field genotype (for DRB3, DRB4, DRB5 also '
                    '"absent"); High confidence, no weakness in the read evidence. Shaded: genotype share below 97% or '
                    'high-confidence share below 92%.']),
             ('', [('b', 'Platform and case status. '), 'Every control was sequenced on HiSeq X at 15× and every case on '
                   f"DNBSeq or NovaSeq; the {g4['samples']} DNBSeq-G400RS 15× cases hold {g4['outliers']} of the "
                   f"{allc['outliers']} case outliers and the lowest call quality (planned platform checks; association "
                   'sheet A10).']),
             ab('HLA', 'IPD-IMGT/HLA')]
    return finish(sh, notes, 'T1 Typing run')


def t2(B):
    sh = B.sheet('T2 DRB3, DRB4, DRB5', 'HLA typing: DRB3, DRB4 and DRB5 copy number against DRB1',
                 'Expected copies follow from the two DRB1 alleles (DR51: DRB5; DR52: DRB3; DR53: DRB4) · '
                 'one row per gene', 'Typing')
    d = records(B.t('t2_drb.tsv'))
    cols = [C('gene', 'Gene', 'gene', 8), C('with_drb1', 'Genomes with DRB1 known', 'int', 12, note='a'),
            C('agree', 'Copies as DRB1 predicts', 'int', 12, blk='Before reconciliation', note='a'),
            C('agree_share', 'Share', 'pct1', 8, blk='Before reconciliation'),
            C('two_where_one', 'Two alleles where one expected', 'int', 12, blk='Before reconciliation'),
            C('present_where_none', 'Present where none expected', 'int', 12, blk='Before reconciliation'),
            C('no_flag', 'No flag', 'int', 9, blk='Reconciliation flags', note='b'),
            C('reconciled', 'Second allele dropped', 'int', 10, blk='Reconciliation flags', note='b'),
            C('unexpected', 'Set absent', 'int', 9, blk='Reconciliation flags', note='b'),
            C('missing', 'Expected, not typed', 'int', 9, blk='Reconciliation flags', note='b'),
            C('drb1_unknown', 'DRB1 unknown', 'int', 9, blk='Reconciliation flags', note='b'),
            C('with_result', 'With a result', 'int', 10, blk='Exported copies', note='c'),
            C('exported_match', 'Equal to DRB1', 'int', 10, blk='Exported copies', note='c'),
            C('exported_blank', 'Blank', 'int', 8, blk='Exported copies', note='c')]
    sp = {k: ([('b', k)], INK, None) for k in ('Before reconciliation', 'Reconciliation flags', 'Exported copies')}
    put(sh, 'Copies against DRB1', cols, d, span=sp)
    two = records(B.t('t2_two_allele.tsv'))
    put(sh, 'Two-allele calls with ≥ 10 unique reads', [
        C('expected_copies', 'Copies DRB1 predicts', 'int', 14), C('two_allele_calls', 'Calls of two alleles', 'int', 12),
        C('minor_below_15pct', 'Minor allele < 15% of unique reads', 'int', 18, note='d')], two)
    notes = [('a', ['Genomes with DRB1 known, both DRB1 alleles typed; Copies as DRB1 predicts, HLA-HD\'s raw call has '
                    'the number of copies the DR group implies.']),
             ('b', ['Reconciliation (typing step 01): two alleles where one copy is expected keep the one with more '
                    'unique reads; a gene typed where none is expected is set absent; an expected gene not typed is '
                    'left missing.']),
             ('c', ['Exported copies, the sum of the gene\'s allele dosages in the typing export; every genome with a '
                    'result carries the copies DRB1 implies.']),
             ('d', ['A raw second allele carried by fewer than 15% of the unique reads is the usual sign of a paralogue '
                    'read misplaced, not of a second copy.']),
             ab('HLA')]
    return finish(sh, notes, 'T2 DRB3, DRB4, DRB5')


def t3(B):
    sh = B.sheet('T3 Completeness', 'HLA typing: completeness, read evidence and locus tier per gene',
                 'All 3,101 genomes · 19 reported genes · calls by result, by how the allele was set, by read support '
                 'and by confidence', 'Typing')
    g = B.t('t3_genes.tsv')
    tot = B.tf['totals']
    rows = records(g)
    for r in rows:
        r['criteria'] = (r['criteria_failed'] or '').replace('call_quality', 'call quality').replace(
            'unconfirmed', 'absent from ToMMo').replace(';', ', ') or 'none'
    total = {'gene': 'All', 'calls': tot['calls'], 'genotypes': tot['genotypes'],
             'genotype_share': tot['genotypes'] / tot['calls'], 'result_share': tot['result'] / tot['calls'],
             'imgt_alleles': tot['imgt_alleles'], 'assoc': NA, 'median_min_reads': NA, 'tier': NA, 'criteria': NA}
    for k in ('heterozygous', 'homozygous', 'present', 'absent', 'ambiguous', 'no_result', 'basis_agreed', 'basis_prior',
              'basis_pgroup', 'weak_prior', 'low_support', 'het_assessed', 'imbalanced', 'conf_high', 'conf_medium',
              'conf_low', 'conf_no_call'):
        total[k] = int(g[k].sum())
    total['share_high'] = total['conf_high'] / tot['calls']
    rows.append(total)
    cols = [C('gene', 'Gene', 'gene', 7), C('assoc', 'Tested', 'bool', 7, note='a'),
            C('imgt_alleles', 'Alleles named', 'int', 9, note='b'),
            C('calls', 'Calls', 'int', 7),
            C('heterozygous', 'Hetero­zygous', 'int', 9, blk='Calls by result'), C('homozygous', 'Homo­zygous', 'int', 9, blk='Calls by result'),
            C('present', 'Present', 'int', 8, blk='Calls by result', note='c'), C('absent', 'Absent', 'int', 8, blk='Calls by result', note='c'),
            C('ambiguous', 'First field only', 'int', 8, blk='Calls by result'), C('no_result', 'No result', 'int', 7, blk='Calls by result'),
            C('genotype_share', 'Genotype', 'pct1', 9, blk='Calls by result', note='d'),
            C('basis_agreed', 'Reads agree', 'int', 9, blk='Allele set by', note='e'), C('basis_prior', 'Frequency table', 'int', 9, blk='Allele set by', note='e'),
            C('basis_pgroup', 'P group', 'int', 8, blk='Allele set by', note='e'),
            C('median_min_reads', 'Median reads, weaker allele', 'f1', 10, blk='Read evidence'),
            C('het_assessed', 'Hetero­zygotes assessed', 'int', 10, blk='Read evidence', note='f'),
            C('imbalanced', 'Imbalanced', 'int', 9, blk='Read evidence', note='f'),
            C('conf_high', 'High', 'int', 8, blk='Confidence', note='g'), C('conf_medium', 'Medium', 'int', 8, blk='Confidence'),
            C('conf_low', 'Low', 'int', 7, blk='Confidence'), C('conf_no_call', 'No call', 'int', 7, blk='Confidence'),
            C('share_high', 'High, %', 'pct1', 8, blk='Confidence'),
            C('tier', 'Tier', 'int', 6, blk='Locus tier', note='h'), C('criteria', 'Criteria failed', 'text', 30, blk='Locus tier')]
    sp = {k: ([('b', k)], INK, None) for k in ('Calls by result', 'Allele set by', 'Read evidence', 'Confidence', 'Locus tier')}
    put(sh, 'Per gene', cols, rows, span=sp, bold=lambda r: r['gene'] == 'All', freeze=1)
    notes = [('a', ['Tested, the gene enters the association; DRB3, DRB4 and DRB5 do not, as their presence follows DRB1.']),
             ('b', ['Alleles named, alleles of the gene in IPD-IMGT/HLA (release history 3.64.0).']),
             ('c', ['Present / Absent, DRB3, DRB4 and DRB5 only: one copy typed, or no copy (a result).']),
             ('d', ['Genotype, share of calls with a two-field genotype (DRB3/4/5 absent counts).']),
             ('e', [f"Allele set by, over allele slots ({tot['slots']:,} in all): the candidates agree to two fields "
                    f"({tot['agreed']:,}); HLA-HD's frequency table chose ({tot['prior']:,}, of which {tot['prior_decisive']:,} "
                    f"with a ratio ≥ 100); all candidates in one P group ({tot['pgroup']:,}); none ({tot['none']:,})."]),
             ('f', [f"Imbalanced, heterozygote whose minor allele holds < 15% of ≥ 10 unique reads; genes every haplotype "
                    f"carries: {tot['imbalanced_obligate']:,} of {tot['het_assessed_obligate']:,}."]),
             ('g', ['Confidence: high, no weakness; medium, one; low, two or more, or first field only; no call, no result.']),
             ('h', [f"Locus tier, 1 + criteria failed (absent share in ToMMo ≤ 2%, dropout lower bound ≤ 3%, high-confidence "
                    f"share ≥ 95%), at most 3; genes without a panel start at 2. Tiers: "
                    f"{', '.join(f'{k}: {v}' for k, v in tot['tier_counts'].items())} genes."]),
             *missing_block([('–', 'All', 'not applicable (a per-gene property)')]),
             ab('HLA', 'IPD-IMGT/HLA', 'ToMMo')]
    return finish(sh, notes, 'T3 Completeness')


def t4(B):
    sh = B.sheet('T4 ToMMo comparison', 'HLA typing: our controls against ToMMo 61KJPN and 1KG JPT',
                 'Controls only (2,662), at P-group resolution · null alleles dropped on both sides · the frequency '
                 'comparison of the typing deck', 'Typing')
    h = records(B.t('t4_harmonization.tsv'))
    put(sh, 'Making the two tables comparable', [C('item', 'Step', 'text', 9, span=7),
                                                 C('value', 'Value', 'text', 9, span=7)],
        [{**r, 'value': kv(r['value']), '_key': r['deck_key']} for r in h])
    p = B.t('t4_panels.tsv')
    rows = records(p)
    for r in rows:
        r['comparable'] = 'Yes' if r['frequency_comparable'] == 1 else 'No'
        for k in ('two_field', 'ours_keys', 'tommo_keys', 'alleles_in_panel', 'alleles_absent', 'alleles_too_new',
                  'copies_absent', 'r_shared'):
            if blank(r.get(k)) and (r['panel'] != 'ToMMo 61KJPN' or r['comparable'] == 'No'):
                r[k] = NA
    cols = [C('panel', 'Panel', 'text', 14), C('gene', 'Gene', 'gene', 7), C('comparable', 'Compared', 'text', 9, note='a'),
            C('n_samples', 'Controls', 'int', 8), C('n_chr_ours', 'Ours', 'int', 8, blk='Copies'),
            C('n_chr_panel', 'Panel', 'int', 9, blk='Copies'),
            C('two_field', 'Two-field alleles', 'int', 9, blk='Our alleles'), C('ours_keys', 'P-group keys', 'int', 9, blk='Our alleles', note='b'),
            C('tommo_keys', 'ToMMo keys', 'int', 9, blk='Keys'), C('n_shared', 'In both', 'int', 8, blk='Keys'),
            C('n_alleles', 'In either', 'int', 8, blk='Keys'),
            C('alleles_in_panel', 'In ToMMo', 'int', 8, blk='Presence of our keys', note='c'),
            C('alleles_absent', 'Absent', 'int', 8, blk='Presence of our keys', note='c'),
            C('alleles_too_new', 'Named later', 'int', 8, blk='Presence of our keys', note='c'),
            C('copies_absent', 'Copies absent', 'int', 9, blk='Presence of our keys'),
            C('r', 'r', 'f3', 7, blk='Frequency agreement', note='d'), C('r_shared', 'r, keys in both', 'f3', 9, blk='Frequency agreement'),
            C('D', 'D', 'f3', 7, blk='Frequency agreement', note='d'), C('excess_D', 'Excess D', 'f3', 8, blk='Frequency agreement', note='d')]
    sp = {k: ([('b', k)], INK, None) for k in ('Copies', 'Our alleles', 'Keys', 'Presence of our keys', 'Frequency agreement')}
    put(sh, 'Per gene and panel', cols, rows, span=sp, group='panel')
    o = records(B.t('t4_outlier.tsv'))
    for r in o:
        pct_row = r['item'].endswith((', frequency', 'group, our controls', 'group, ToMMo'))
        r['value'] = F.pct(float(r['value']), 2) if pct_row else kv(r['value'])
        r['_key'] = r['item']
    put(sh, f'The one allele off the diagonal ({sty("outlier")})', [C('item', 'Item', 'text', 9, span=7),
                                                          C('value', 'Value', 'text', 9, span=7)], o)
    hz = records(B.t('t4_homozygosity.tsv'))
    put(sh, 'Homozygosity and dropout (genes every haplotype carries; not on the slides)', [
        C('gene', 'Gene', 'gene', 7), C('n_samples', 'Genomes', 'int', 8), C('obs_hom', 'Observed', 'pct1', 9, blk='Homozygotes'),
        C('exp_hom', 'Expected', 'pct1', 9, blk='Homozygotes', note='e'), C('ratio', 'Ratio', 'f2', 7, blk='Homozygotes'),
        C('dropout', 'Estimate', 'pct1', 9, blk='Allele dropout', note='e'), C('dropout_lo', 'Lower 95%', 'pct1', 9, blk='Allele dropout'),
        C('dropout_hi', 'Upper 95%', 'pct1', 9, blk='Allele dropout')], hz,
        span={k: ([('b', k)], INK, None) for k in ('Homozygotes', 'Allele dropout')})
    cv = B.t('t4_confidence_vs_panel.tsv')
    cv = records(cv[cv.gene == 'all ToMMo loci'])
    put(sh, 'Absent from ToMMo, by call confidence (11 loci pooled; not on the slides)', [
        C('confidence', 'Confidence', 'text', 12), C('n_chr', 'Copies', 'int', 9), C('n_unconfirmed', 'Absent from ToMMo', 'int', 10),
        C('unconfirmed_share', 'Share', 'pct2', 8)], cv)
    sh_ = B.tf['shift']
    notes = [('a', ['Compared: DRB3 and DRB4 keep their rows but their frequencies are not comparable (ToMMo gives every '
                    'chromosome a DRB3 allele; its non-null DRB4 exceeds the DR53 haplotypes of its own DRB1); DRA, DMA, '
                    'DMB, DOA, DOB and DRB5 have no panel.']),
             ('b', ['P-group keys, our two-field alleles at P-group resolution (alleles identical in protein over the '
                    'peptide-binding domain share one key), as ToMMo names its alleles.']),
             ('c', ['Presence of our keys: in ToMMo; absent though ToMMo could have carried them; or named after ToMMo\'s '
                    'newest allele (release 3.21.0), so absent whatever their frequency.']),
             ('d', ['r, Pearson correlation of allele frequencies over the keys in either table; D, total variation '
                    'distance; Excess D, D beyond the sampling expectation (1,000 draws).']),
             ('e', ['Expected homozygosity under Hardy–Weinberg from our allele frequencies; dropout, the share of '
                    'heterozygotes called homozygous that the excess implies.']),
             ('', [('b', 'Shift direction. '), f"Of the {sh_['n']} alleles at least one point from ToMMo at the 1KG loci, "
                   f"1KG JPT lies on our side for {sh_['kg_on_our_side']} and on ToMMo's for {sh_['kg_on_tommo_side']}."]),
             *missing_block([('–', 'Per gene and panel', 'not assessed: the presence of our keys is assessed against '
                              'ToMMo only, and not for DRB3 and DRB4, whose frequencies are not comparable'),
                             ('—', 'r, keys in both', 'fewer than three keys in both tables')]),
             ab('1KG JPT', 'HLA', 'ToMMo')]
    return finish(sh, notes, 'T4 ToMMo comparison')


def t5(B):
    sh = B.sheet('T5 Marker QC', 'HLA typing: the allele and residue columns the association reads',
                 'Encoding, the floors and the stepwise record per cohort · call rate ≥ 95% in cases and in controls, '
                 'then MAF ≥ 0.5% · DRB3, DRB4 and DRB5 not tested', 'Typing')
    rec = B.t('t5_record.tsv')
    stages = list(dict.fromkeys(rec.stage))
    names = dict(zip(rec.stage, rec.stage_name))
    rows = []
    for st in stages:
        row = {'stage': names[st], '_key': st}
        for c in COHORTS:
            for k in ('allele', 'residue'):
                row[f'{k}_{c}'] = int(rec[(rec.stage == st) & (rec.cohort == c) & (rec['class'] == k)]['columns'].iloc[0])
        rows.append(row)
    cols = [C('stage', 'Step', 'text', 40)] + cohort_cols([('allele', 'Allele columns', 'int', 10),
                                                           ('residue', 'Residue columns', 'int', 10)])
    put(sh, f'Stepwise record ({sty("record")})', cols, rows, span=coh_span(), bold=lambda r: r['_key'] == 'maf')
    g = B.t('t5_genes.tsv')
    rows = []
    for gene in list(dict.fromkeys(g.gene)):
        row = {'gene': gene}
        for c in COHORTS:
            d = g[(g.gene == gene) & (g.cohort == c)].iloc[0]
            row[f'a_typed_{c}'] = d.allele_typed
            row[f'a_maf_{c}'] = d.allele_maf
            row[f'cov_{c}'] = d.allele_copies_covered
            row[f'r_typed_{c}'] = d.residue_typed
            row[f'r_cr_{c}'] = d.residue_callrate
            row[f'r_maf_{c}'] = d.residue_maf
            if gene in ('DRB3', 'DRB4', 'DRB5'):
                for k in ('a_maf', 'cov', 'r_cr', 'r_maf'):
                    row[f'{k}_{c}'] = NA
        rows.append(row)
    fields = [('a_typed', 'Typed', 'int', 7), ('a_maf', 'Tested', 'int', 7), ('cov', 'Copies covered', 'pct1', 8),
              ('r_typed', 'Typed', 'int', 7), ('r_cr', 'Call rate', 'int', 7), ('r_maf', 'Tested', 'int', 7)]
    cols = [C('gene', 'Gene', 'gene', 7)]
    for c in COHORTS:
        for k, leaf, kind, w in fields:
            cols.append(C(f'{k}_{c}', leaf, kind, w, blk=COH_LABEL[c], sub='Allele columns' if k in ('a_typed', 'a_maf', 'cov')
                          else 'Residue columns', coh=c, note='a' if k == 'cov' else None))
    put(sh, f'Per gene ({sty("alleles", "residues")})', cols, rows, span=coh_span(), bold=lambda r: r['gene'] == 'All', freeze=1)
    fl = B.t('t5_floors.tsv')
    rows = records(fl)
    cols = [C('cohort', 'Cohort', 'cohort', 20), C('class', 'Columns', 'text', 9),
            C('columns', 'Exported', 'int', 8), C('not_varying', 'Not varying', 'int', 8), C('varying', 'Varying', 'int', 8),
            C('call_rate_min', 'Lowest', 'pct1', 8, blk='Call rate of the varying columns', note='b'),
            C('call_rate_median', 'Median', 'pct1', 8, blk='Call rate of the varying columns'),
            C('lost_at_95', 'Lost at 95%', 'int', 8, blk='Call rate of the varying columns'),
            C('lost_at_98', 'Lost at 98%', 'int', 8, blk='Call rate of the varying columns'),
            C('lost_at_99', 'Lost at 99%', 'int', 8, blk='Call rate of the varying columns'),
            C('pass_call_rate', 'Pass', 'int', 8, blk='Call rate of the varying columns'),
            C('maf_<0.1%', '< 0.1%', 'int', 7, blk='MAF of the columns passing the call rate'),
            C('maf_0.1-0.5%', '0.1–0.5%', 'int', 7, blk='MAF of the columns passing the call rate'),
            C('maf_0.5-1%', '0.5–1%', 'int', 7, blk='MAF of the columns passing the call rate'),
            C('maf_1-5%', '1–5%', 'int', 7, blk='MAF of the columns passing the call rate'),
            C('maf_5-10%', '5–10%', 'int', 7, blk='MAF of the columns passing the call rate'),
            C('maf_10-25%', '10–25%', 'int', 7, blk='MAF of the columns passing the call rate'),
            C('maf_>=25%', '≥ 25%', 'int', 7, blk='MAF of the columns passing the call rate'),
            C('kept_maf05', 'MAF ≥ 0.5% (applied)', 'int', 10, blk='Kept under each floor', note='c'),
            C('floor_mac_maf05', 'Smallest MAC kept', 'int', 9, blk='Kept under each floor'),
            C('kept_mac20', 'MAC ≥ 20', 'int', 8, blk='Kept under each floor'),
            C('kept_maf1', 'MAF ≥ 1%', 'int', 8, blk='Kept under each floor')]
    sp = {k: ([('b', k)], INK, None) for k in ('Call rate of the varying columns', 'MAF of the columns passing the call rate',
                                               'Kept under each floor')}
    put(sh, f'The two floors ({sty("callrate", "macmaf")})', cols, rows, span=sp, group='cohort')
    comp = records(B.t('t5_completeness.tsv'))
    for r in comp:
        r['worst_position'] = str(r['worst_position']).replace('-', '−')
    cols = [C('gene', 'Gene', 'gene', 7), C('all_genomes', 'All genomes', 'pct2', 9, blk='Genotype share'),
            C('cases_full', 'Cases', 'pct2', 8, blk='Genotype share', note='d'), C('controls_full', 'Controls', 'pct2', 8, blk='Genotype share', note='d'),
            C('narrow_mainland', 'N', 'pct2', 8, blk='Genotype share'), C('intermediate_mainland', 'I', 'pct2', 8, blk='Genotype share'),
            C('full_mainland', 'F', 'pct2', 8, blk='Genotype share'),
            C('worst_position', 'Position', 'text', 8, blk='Least complete residue position', align='center'),
            C('worst_position_call_rate', 'Share determined', 'pct2', 10, blk='Least complete residue position')]
    put(sh, 'Gene completeness: the ceiling of every column\'s call rate', cols, comp,
        span={k: ([('b', k)], INK, None) for k in ('Genotype share', 'Least complete residue position')})
    ex = B.tf['markers']
    notes = [('a', ['Copies covered, share of the gene\'s typed allele copies that belong to a tested allele column; the '
                    'rest are alleles below the MAF floor, which enter the association through the residue columns only. '
                    'All row: also the residue copies covered at the positions kept.']),
             ('b', ['Call rate of a column, the lower of its share of genomes with a genotype in cases and in controls '
                    '(the rule is applied in each group, because platform and depth differ between them).']),
             ('c', ['Kept under each floor, columns passing the call rate and the floor; MAF ≥ 0.5% is the floor '
                    'applied, the others are shown for comparison.']),
             ('d', ['Cases and controls of the full mainland cohort.']),
             ('', [('b', 'Encoding. '), f"Exports: {ex['allele_columns_exported']:,} allele columns and "
                   f"{ex['residue_columns_exported']:,} residue columns at {ex['polymorphic_positions_assoc_genes']:,} "
                   'polymorphic positions of the 16 tested genes (each column counts copies, 0, 1 or 2; a residue left '
                   'undetermined by the alignment leaves the whole position blank for that genome). '
                   f"{ex['positions_losing_over_1pct']} of {ex['residue_positions']:,} residue positions lose more than 1% "
                   f"of the gene's completeness to the alignment, {ex['positions_losing_over_5pct']} more than 5%."]),
             *missing_block([('–', 'Per gene', 'DRB3, DRB4 and DRB5 are not tested (their presence follows DRB1)')]),
             ab('MAC', 'MAF', 'HLA', 'N / I / F')]
    return finish(sh, notes, 'T5 Marker QC')


def t6(B):
    sh = B.sheet('T6 P groups', 'HLA typing: P groups, the DQA1 example and tested columns that share a P group',
                 'ToMMo reports P groups; the association tests two-field alleles · full mainland unless stated', 'Typing')
    sc = records(B.t('t6_split_counts.tsv'))
    put(sh, f'Tested allele columns that share a P group, per cohort ({sty("pgsplit")})', [
        C('cohort', 'Cohort', 'cohort', 20), C('tested', 'Allele columns tested', 'int', 12),
        C('groups', 'P groups', 'int', 9, blk='Shared by ≥ 2 tested columns'), C('columns', 'Columns', 'int', 9, blk='Shared by ≥ 2 tested columns'),
        C('groups_compared', 'P groups', 'int', 9, blk='Of them in genes compared with ToMMo'),
        C('columns_compared', 'Columns', 'int', 9, blk='Of them in genes compared with ToMMo')], sc,
        span={k: ([('b', k)], INK, None) for k in ('Shared by ≥ 2 tested columns', 'Of them in genes compared with ToMMo')})
    rows = records(B.t('t6_split_groups.tsv'))
    for r in rows:
        if r['compared'] in (False, 'False', 0):
            r['tommo_freq'] = NA
        tf = [float(x) for x in str(r['tested_freq_controls']).split('; ')]
        r['tested_list'] = '; '.join(f'{a} ({F.pct(f, 2)})' for a, f in zip(str(r['tested']).split('; '), tf))
    put(sh, f'The P groups with more than one tested column, full mainland ({sty("pgsplit")})', [
        C('gene', 'Gene', 'gene', 7), C('p_group', 'P group', 'text', 14), C('compared', 'Compared with ToMMo', 'bool', 10),
        C('tested_list', 'Tested columns (frequency in controls)', 'text', 11, span=5, wrap=True),
        C('others_n', 'Other alleles', 'int', 8, note='a'), C('others_freq', 'Their frequency', 'pct2', 9),
        C('group_freq_controls', 'Our controls', 'pct2', 9, blk='P group frequency'),
        C('tommo_freq', 'ToMMo', 'pct2', 9, blk='P group frequency')], rows,
        span={'P group frequency': ([('b', 'P group frequency')], INK, None)})
    ex = records(B.t('t6_dqa1.tsv'))
    for r in ex:
        r['where_label'] = {'both': 'Both', 'ours': 'Ours only', 'tommo': 'ToMMo only'}[r['where']]
        if r['where'] == 'tommo':
            r['members'] = NA
    put(sh, f'One gene through the comparison: DQA1 ({sty("pgexample")})', [
        C('p_group', 'P group', 'text', 14), C('where_label', 'Found in', 'text', 11),
        C('copies_ours', 'Ours', 'int', 8, blk='Copies'), C('copies_tommo', 'ToMMo', 'int', 9, blk='Copies'),
        C('freq_ours', 'Ours', 'pct2', 9, blk='Frequency'), C('freq_tommo', 'ToMMo', 'pct2', 9, blk='Frequency'),
        C('members', 'Our two-field alleles in it (copies)', 'text', 11, span=7)], ex,
        span={k: ([('b', k)], INK, None) for k in ('Copies', 'Frequency')}, group='where_label')
    notes = [('a', ['Other alleles, our alleles of the same P group without a column of their own (below the MAF '
                    'floor); the P group\'s frequency counts them.']),
             ('', ['A P group holds alleles identical in protein over the peptide-binding domain; ToMMo names alleles by '
                   'P group, so the comparison with ToMMo cannot tell apart tested columns of one group.']),
             *missing_block([('–', 'ToMMo', 'the gene is not compared with ToMMo'),
                             ('–', 'Our two-field alleles', 'a P group found in ToMMo only')]),
             ab('HLA', 'MAF', 'ToMMo')]
    return finish(sh, notes, 'T6 P groups')


GENE_ORDER = ['A', 'B', 'C', 'E', 'F', 'G', 'DRA', 'DRB1', 'DRB3', 'DRB4', 'DRB5', 'DQA1', 'DQB1', 'DPA1', 'DPB1', 'DMA',
              'DMB', 'DOA', 'DOB']


def allele_key(a):
    """Natural order of allele names: gene in the order of the sheets, then the fields as numbers."""
    g, rest = a.split('*')
    f = re.findall(r'\d+|\D+', rest)
    return (GENE_ORDER.index(g), [int(x) if x.isdigit() else x for x in f if x != ':'])


STATUS_ORDER = ['Tested', 'Below MAF 0.5%', 'Call rate below 95%', 'Gene not tested', 'Not varying']


def t7(B):
    sh = B.sheet('T7 Alleles (detail)', 'HLA typing: every exported allele column',
                 '1,995 two-field alleles · per cohort: copies in cases and controls, call rate, MAF and whether the '
                 'association tests it · sortable (filter on the header row)', 'Typing')
    d = B.t('t7_alleles.tsv')
    rows = sorted(records(d), key=lambda r: allele_key(r['marker']))
    for r in rows:
        fr = r.get('first_release')
        r['release'] = f'{str(int(fr))[0]}.{str(int(fr))[1:3]}.{str(int(fr))[3]}' if fr else None
        r['_key'] = r['marker']
        null = r['marker'].endswith('N')
        for k, ok in (('tommo_freq', 'tommo_compared'), ('kg_freq', 'kg_compared')):
            if not r[ok] or null:
                r[k] = NA
    cols = [C('gene', 'Gene', 'gene', 7), C('marker', 'Allele', 'text', 14), C('p_group', 'P group', 'text', 13, note='a'),
            C('release', 'First named', 'text', 9, note='b'), C('match', 'Sequence', 'text', 11, note='c')]
    for c, t in zip(COHORTS, 'nif'):
        cols += [C(f'{t}_status', 'Status', 'text', 15, blk=COH_LABEL[c], coh=c, note='d'),
                 C(f'{t}_copies_cases', 'Copies, cases', 'int', 8, blk=COH_LABEL[c], coh=c),
                 C(f'{t}_copies_controls', 'Copies, controls', 'int', 9, blk=COH_LABEL[c], coh=c),
                 C(f'{t}_cr_cases', 'Call rate, cases', 'pct1', 8, blk=COH_LABEL[c], coh=c),
                 C(f'{t}_cr_controls', 'Call rate, controls', 'pct1', 9, blk=COH_LABEL[c], coh=c),
                 C(f'{t}_maf', 'MAF', 'pct2', 8, blk=COH_LABEL[c], coh=c)]
    cols += [C('tommo_freq', 'ToMMo', 'pct2', 8, blk='P group in the panels', note='e'),
             C('kg_freq', '1KG JPT', 'pct2', 8, blk='P group in the panels', note='e')]
    put(sh, None, cols, rows, span=coh_span({'P group in the panels': ([('b', 'P group in the panels')], INK, None)}),
        group='gene', freeze=2, autofilter=True)
    notes = [('a', ['P group as ToMMo writes it (alleles identical in protein over the peptide-binding domain, IPD-IMGT/HLA '
                    '3.64.0); an allele in no P group keeps its own name.']),
             ('b', ['First named, the IPD-IMGT/HLA release in which an allele with this two-field name first appeared.']),
             ('c', ['Sequence, how the allele was matched to the protein alignment for its residues: exact, or its first '
                    'deeper allele.']),
             ('d', ['Status, the first rule the column fails, in this order: not varying in the cohort, gene not tested '
                    '(DRB3/4/5), call rate below 95% in cases or in controls, MAF below 0.5%; otherwise tested.']),
             ('e', ['Frequency of the allele\'s P group in ToMMo 61KJPN and 1KG JPT (genes they report).']),
             *missing_block([('–', 'P group in the panels', 'not compared: the gene is not in the panel or its frequencies '
                              'are not comparable (DRB3, DRB4), or a null allele (dropped from the comparison)'),
                             ('—', 'P group in the panels', 'the P group is in neither our controls nor the panel (the '
                              'comparison covers the alleles of our controls)'),
                             ('—', 'First named', 'not in the release history')]),
             ab('1KG JPT', 'IPD-IMGT/HLA', 'MAF', 'ToMMo')]
    return finish(sh, notes, 'T7 Alleles (detail)')


def t8(B):
    sh = B.sheet('T8 Residues (detail)', 'HLA typing: every exported residue column',
                 '2,872 amino-acid residues at their IPD-IMGT/HLA positions · per cohort: call rate, MAF and whether the '
                 'association tests it · sortable (filter on the header row)', 'Typing')
    d = B.t('t8_residues.tsv')
    rows = sorted(records(d), key=lambda r: (GENE_ORDER.index(r['gene']), int(r['position']), str(r['residue'])))
    for r in rows:
        r['pos_label'] = str(r['position']).replace('-', '−')
        r['_key'] = f"{r['gene']}:{r['position']}:{r['residue']}"
    cols = [C('gene', 'Gene', 'gene', 7), C('pos_label', 'Position', 'text', 8, align='center', note='a'),
            C('residue', 'Residue', 'text', 7, align='center')]
    for c, t in zip(COHORTS, 'nif'):
        cols += [C(f'{t}_status', 'Status', 'text', 15, blk=COH_LABEL[c], coh=c, note='b'),
                 C(f'{t}_copies_cases', 'Copies, cases', 'int', 8, blk=COH_LABEL[c], coh=c),
                 C(f'{t}_copies_controls', 'Copies, controls', 'int', 9, blk=COH_LABEL[c], coh=c),
                 C(f'{t}_cr_cases', 'Call rate, cases', 'pct1', 8, blk=COH_LABEL[c], coh=c),
                 C(f'{t}_cr_controls', 'Call rate, controls', 'pct1', 9, blk=COH_LABEL[c], coh=c),
                 C(f'{t}_maf', 'MAF', 'pct2', 8, blk=COH_LABEL[c], coh=c)]
    put(sh, None, cols, rows, span=coh_span(), group='gene', freeze=3, autofilter=True)
    notes = [('a', ['Position, IPD-IMGT/HLA numbering of the mature protein (negative: leader peptide).']),
             ('b', ['Status as in sheet T7 Alleles (detail). A residue undetermined in the alignment for one of a genome\'s '
                    'alleles leaves the whole position blank for that genome, which lowers the call rate.']),
             ab('IPD-IMGT/HLA', 'MAF')]
    return finish(sh, notes, 'T8 Residues (detail)')


# ── Association sheets ──────────────────────────────────────────────────────
TT_LABEL = {'allele': 'Alleles', 'residue': 'Residues', 'position': 'Positions'}
TT_MODEL = {'allele': 'Alleles (random model, GEC)', 'residue': 'Residues (random model, GEC)',
            'position': 'Positions (fixed model, Bonferroni)'}
B_ = lambda k: ([('b', k)], INK, None)  # noqa: E731


def bspan(*keys):
    return {k: B_(k) for k in keys}


def ci(r, o='or', l_='l95', u='u95'):
    return ci_text(r.get(o), r.get(l_), r.get(u))


def a1(B):
    sh = B.sheet('A1 Design', 'HLA association: cohorts, comparisons, tests and models',
                 'Three nested cohorts × four comparisons · alleles and residues one at a time (random model), whole '
                 'positions by an omnibus test (fixed model)', 'Association')
    s = records(B.s('a1_samples.tsv'))
    cols = [C('comparison', 'Comparison', 'comp', 30), C('cohort', 'Cohort', 'cohort', 20),
            C('analysis', 'Analysis', 'text', 17),
            C('random_samples', 'Samples', 'int', 8, blk='Random model', sub='Single markers', note='a'),
            C('random_cases', 'Case side', 'int', 8, blk='Random model', sub='Single markers'),
            C('random_controls', 'Control side', 'int', 8, blk='Random model', sub='Single markers'),
            C('fixed_samples', 'Samples', 'int', 8, blk='Fixed model', sub='Positions', note='b'),
            C('fixed_cases', 'Case side', 'int', 8, blk='Fixed model', sub='Positions'),
            C('fixed_controls', 'Control side', 'int', 8, blk='Fixed model', sub='Positions'),
            C('relatives_removed', 'Relatives removed', 'int', 9, blk='Fixed model', sub='Positions'),
            C('null_model', 'Null model', 'text', 24, note='c'), C('variance_ratio', 'Variance ratio', 'f3', 9, note='c')]
    put(sh, f'Samples per comparison and cohort ({sa("cohorts")})', cols, s, span=bspan('Random model', 'Fixed model'),
        group='comparison', bold=lambda r: r['comparison'] == 'all')
    dv = [r for r in s if r['comparison'] == 'all']
    put(sh, f'History of DVT among the cases ({sa("cohorts")})', [
        C('cohort', 'Cohort', 'cohort', 20), C('random_cases', 'Cases', 'int', 8),
        C('dvt_yes', 'With DVT', 'int', 8, blk='History of DVT'), C('dvt_no', 'Without DVT', 'int', 9, blk='History of DVT'),
        C('dvt_unknown', 'Not recorded', 'int', 9, blk='History of DVT', note='d')], dv, span=bspan('History of DVT'))
    g = records(B.s('a1_tests_by_gene.tsv'))
    cols = [C('gene', 'Gene', 'gene', 7)] + cohort_cols([('allele', 'Alleles', 'int', 8), ('residue', 'Residues', 'int', 8),
                                                         ('position', 'Positions', 'int', 8)])
    put(sh, f'Tests per gene, all cases vs. controls ({sa("markers")})', cols, g, span=coh_span(),
        bold=lambda r: r['gene'] == 'All')
    m = records(B.s('a1_methods.tsv'))
    put(sh, f'Models and settings ({sa("design", "models", "testCount")})', [C('scope', 'Applies to', 'text', 30), C('item', 'Item', 'text', 20),
                                                     C('value', 'Setting', 'text', 10, span=10, wrap=True)], m,
        group='scope')
    notes = [('a', ['Random model: every sample of the comparison, relatives included (the genetic relationship matrix '
                    'accounts for them).']),
             ('b', ['Fixed model: unrelated samples only (one of each related pair removed); the fixed-model single-marker '
                    'tests use the same samples.']),
             ('c', ['Null model, the model without the marker; the variance ratio calibrates the score test of each marker '
                    'to the full genetic relationship matrix.']),
             ('d', ['Not recorded: the case enters all cases vs. controls only.']),
             ('', ['The case sides of the four comparisons: all cases; cases with a history of DVT; cases without one; the '
                   'case-only comparison sets cases with DVT against cases without (no controls). Tests per gene: a mirror '
                   'pair counts once; the counts differ slightly between cohorts because the marker floors are applied in '
                   'each cohort (sheet T5).']),
             ab('BBJ', 'DVT', 'GEC', 'GRM', 'LD', 'LOCO', 'LRT', 'N / I / F', 'OR', 'PC', 'SAIGE', 'SPA', 'df')]
    return finish(sh, notes, 'A1 Design')


def ival(lo, hi, f):
    return None if blank(lo) or blank(hi) else f'{f(lo)}–{f(hi)}'


def a2(B):
    sh = B.sheet('A2 Thresholds', 'HLA association: thresholds per scan and their benchmark',
                 'Single markers: GEC = 0.05 / Me (primary), Bonferroni = 0.05 / tests (sensitivity) · positions: '
                 'Bonferroni over the positions tested · the exact max-T benchmarks the single-marker thresholds',
                 'Association')
    d = B.s('a2_scans.tsv')
    rows = records(d)
    for r in rows:
        r['cls'] = TT_LABEL[r['class']]
        r['_key'] = f"{r['comparison']}|{r['cohort']}|{r['class']}"
    cols = [C('comparison', 'Comparison', 'comp', 30), C('cohort', 'Cohort', 'cohort', 20), C('cls', 'Test type', 'text', 9),
            C('markers', 'Columns', 'int', 8, blk='Tests', note='a'), C('mirror_pairs', 'Mirror pairs', 'int', 8, blk='Tests', note='a'),
            C('tests', 'Tests', 'int', 7, blk='Tests'),
            C('me_gec', 'Me', 'f1', 8, blk='GEC (primary)', note='b'), C('me_share', 'Me / tests', 'pct1', 8, blk='GEC (primary)'),
            C('p_gec', 'Threshold', 'p', 12, blk='GEC (primary)'),
            C('p_bonf', 'Threshold', 'p', 12, blk='Bonferroni (sensitivity)'),
            C('bonf_stricter', 'Stricter than GEC, ×', 'f2', 9, blk='Bonferroni (sensitivity)')]
    put(sh, f'Single-marker thresholds per scan ({sa("testCount", "gecMe", "meTable")})', cols, rows,
        span=bspan('Tests', 'GEC (primary)', 'Bonferroni (sensitivity)'), group='comparison',
        bold=lambda r: r['comparison'] == 'all')
    for r in rows:
        for t in ('score', 'wald'):
            r[f'thr_{t}_iv'] = F.pspan([r[f'boot_p_maxt_{t}_lo'], r[f'boot_p_maxt_{t}_hi']])
            r[f'gec_{t}_iv'] = ival(r[f'boot_fwer_gec_{t}_lo'], r[f'boot_fwer_gec_{t}_hi'], lambda x: F.pct(x, 1))
            r[f'bonf_{t}_iv'] = ival(r[f'boot_fwer_bonf_{t}_lo'], r[f'boot_fwer_bonf_{t}_hi'], lambda x: F.pct(x, 1))
    cols = [C('comparison', 'Comparison', 'comp', 30), C('cohort', 'Cohort', 'cohort', 20), C('cls', 'Test type', 'text', 9),
            C('boot_draws', 'Null draws', 'int', 8, note='c')]
    for t, lab in (('score', 'Score test'), ('wald', 'Wald test')):
        cols += [C(f'boot_p_maxt_{t}', 'Threshold', 'p', 11, blk=lab, sub='Exact max-T', note='d'),
                 C(f'thr_{t}_iv', '95% interval', 'text', 19, blk=lab, sub='Exact max-T', align='center'),
                 C(f'boot_fwer_gec_{t}', 'FWER', 'pct1', 7, blk=lab, sub='At GEC', note='e'),
                 C(f'gec_{t}_iv', '95% interval', 'text', 11, blk=lab, sub='At GEC', align='center'),
                 C(f'boot_fwer_bonf_{t}', 'FWER', 'pct1', 7, blk=lab, sub='At Bonferroni'),
                 C(f'bonf_{t}_iv', '95% interval', 'text', 11, blk=lab, sub='At Bonferroni', align='center'),
                 C(f'boot_replicates_with_na_{t}', 'Draws with a test undefined', 'int', 10, blk=lab, note='f')]
    put(sh, f'Exact max-T benchmark per scan ({sa("benchMethod", "maxtWhy", "benchmark")})', cols, rows,
        span=bspan('Score test', 'Wald test'), group='comparison', bold=lambda r: r['comparison'] == 'all')
    pz = records(B.s('a2_positions.tsv'))
    put(sh, f'Position thresholds per scan ({sa("testCount")})', [
        C('comparison', 'Comparison', 'comp', 30), C('cohort', 'Cohort', 'cohort', 20),
        C('positions', 'Polymorphic positions', 'int', 11), C('tested', 'Tested', 'int', 8, note='g'),
        C('not_tested', 'Not tested', 'int', 8, note='g'), C('separated', 'Separated', 'int', 9, note='g'),
        C('threshold', 'Bonferroni threshold', 'p', 13)], pz, group='comparison', bold=lambda r: r['comparison'] == 'all')
    gc = records(B.s('a2_gec_check.tsv'))
    for r in gc:
        r['cls'] = TT_LABEL[r['class']]
    put(sh, f'Me recomputed from the genotypes, all cases vs. controls, full mainland ({sa("gec")})', [
        C('cls', 'Test type', 'text', 9), C('tests', 'Tests', 'int', 7), C('blocks', 'Blocks', 'int', 7, note='h'),
        C('duplicates', 'Duplicates set aside', 'int', 10, note='h'), C('me_gec', 'GEC', 'f1', 8, blk='Me'),
        C('me_recomputed', 'Recomputed', 'f1', 10, blk='Me'), C('max_block_diff', 'Largest difference in a block', 'f3', 12, blk='Me')],
        gc, span=bspan('Me'))
    notes = [('a', ['Columns, allele or residue columns tested; a mirror pair is two residue columns that are complements in '
                    'every sample (a two-residue position) and counts as one test.']),
             ('b', ['Me, the effective number of independent tests from the correlation of the scan\'s markers (GEC v0.2, '
                    'blocks of at most 100 markers); the GEC threshold is 0.05 / Me.']),
             ('c', ['Null draws: case status drawn (Bernoulli) from the null logistic model of sex and PC1–PC10 on the '
                    'comparison\'s unrelated samples, so the markers have no effect; every test of the scan is rerun on '
                    'each draw. Score test, the score test without SPA (the P of the random model\'s calls); Wald test, the '
                    'fixed model\'s test.']),
             ('d', ['Exact max-T, the benchmark: the 0.05 quantile of the smallest P per draw, the threshold that holds the '
                    'family-wise error rate at exactly 0.05 for this scan; it calls no marker. 95% interval from 1,000 '
                    'resamples of the draws.']),
             ('e', ['FWER, share of draws whose smallest P is below the threshold; Wilson 95% interval.']),
             ('f', ['Draws with a test undefined: draws in which a test returned no P (non-convergence or separation, '
                    'rare markers in the Wald test); the smallest P of such a draw is taken over the tests that returned '
                    'one.']),
             ('g', ['Tested, positions with at least two residue groups after residues below 1% are pooled; separated, '
                    'tested but with an estimate at the boundary (its P is kept, not trusted, and left out of λGC).']),
             ('h', ['Blocks as GEC set them; duplicates, markers with |r| ≥ 0.9999 to one already in the block (GEC '
                    'counts them once). Recomputed with the published polynomial of GEC on the same unrelated samples.']),
             ab('FWER', 'GEC', 'Me', 'PC', 'SPA', 'max-T', 'λGC')]
    return finish(sh, notes, 'A2 Thresholds')


SCAN_COLS = lambda: [  # noqa: E731
    C('tt', 'Test type', 'text', 30), C('tests', 'Tests', 'int', 7),
    C('threshold_primary', 'Threshold', 'p', 12, blk='Primary', note='a'),
    C('sig_primary', 'Significant', 'int', 9, blk='Primary', hi=lambda r: (r.get('sig_primary') or 0) > 0),
    C('threshold_bonf', 'Threshold', 'p', 12, blk='Bonferroni (sensitivity)'),
    C('sig_bonf', 'Significant', 'int', 9, blk='Bonferroni (sensitivity)'),
    C('sig_spa', 'Primary, SPA P', 'int', 9, note='b'),
    C('min_p', 'P', 'p', 12, blk='Smallest P'), C('top', 'Marker', 'marker', 14, blk='Smallest P'),
    C('top_or', 'OR', 'f2', 7, blk='Smallest P'), C('lambda_gc', 'λGC', 'f2', 7, note='c')]


def scan_rows(df):
    rows = records(df)
    for r in rows:
        r['tt'] = TT_MODEL[r['test_type']]
        r['_key'] = f"{r['comparison']}|{r['cohort']}|{r['test_type']}"
        if r['test_type'] == 'position':
            r['threshold_bonf'] = r['sig_bonf'] = r['sig_spa'] = r['top_or'] = NA
    return rows


def a3(B):
    sh = B.sheet('A3 Q1 Association', 'HLA association, Q1: which alleles, residues or positions are associated?',
                 'All cases vs. controls in three cohorts · the main analysis', 'Association')
    sc = B.s('a3_scans.tsv')
    for c in COHORTS:
        d = sc[(sc.comparison == 'all') & (sc.cohort == c)]
        cols = [C('cohort', 'Cohort', 'cohort', 20)] + SCAN_COLS()
        put(sh, None if c != COHORTS[0] else f'Every scan ({sa("resAllele", "resResidue", "resPos")})', cols, scan_rows(d),
            span=bspan('Primary', 'Bonferroni (sensitivity)', 'Smallest P'))
        sh.row -= 1
    sh.row += 1
    top = B.s('a3_top.tsv')
    for k in ('allele', 'residue', 'position'):
        d = records(top[top.test_type == k])
        for r in d:
            r['_key'] = f"{k}|{r['marker']}"
            if k == 'position':
                r['res_df'] = f"{r['residues']} ({int(r['df'])})"
        cols = [C('rank', 'Rank', 'int', 5), C('marker', {'allele': 'Allele', 'residue': 'Residue',
                                                              'position': 'Position'}[k], 'marker', 14)]
        if k == 'position':
            cols.append(C('res_df', 'Residues (df)', 'text', 16, align='center'))
        else:
            cols += [C('freq_case', 'Cases', 'pct1', 8, blk='Frequency', note='d'),
                     C('freq_control', 'Controls', 'pct1', 8, blk='Frequency', note='d')]
        cols += cohort_cols([('p', 'P', 'p', 12)] + ([] if k == 'position' else [('or', 'OR', 'f2', 7)]))
        cols += [C('n3_primary', 'Primary', 'text', 8, blk='Significant in', val=lambda r: f"{int(r['n3_primary'])}/3",
                   hi=lambda r: r['n3_primary'] > 0)]
        if k != 'position':
            cols.append(C('n3_bonf', 'Bonferroni', 'text', 9, blk='Significant in', val=lambda r: f"{int(r['n3_bonf'])}/3"))
        cap = {'allele': f'The five strongest alleles ({sa("resAllele")})', 'residue': 'The five strongest residues, one per '
               f'position ({sa("resResidue")})', 'position': f'The five strongest positions ({sa("resPos")})'}[k]
        put(sh, cap, cols, d, span=coh_span(bspan('Frequency', 'Significant in')))
    notes = [('a', ['Primary threshold: GEC for alleles and residues (random model); Bonferroni over the positions tested '
                    'for positions (fixed-model omnibus test). The conventional Bonferroni of single markers is a '
                    'sensitivity analysis.']),
             ('b', ['Primary, SPA P: markers passing the primary threshold when the saddlepoint-approximated P is used '
                    'instead of the score-test P.']),
             ('c', ['λGC, median χ² of the scan\'s tests (each P read as one degree of freedom) over its expectation; the '
                    'HLA markers are chosen around a known region and correlated, so λGC above 1 here is not a measure of '
                    'confounding (sheet A6).']),
             ('d', ['Frequency, copy frequency in the full mainland cohort.']),
             ('', ['Significant in n/3: number of cohorts in which the marker passes the threshold; the cohorts are nested '
                   'and share most cases, so agreement shows robustness to the ancestry definition, not replication.']),
             *missing_block([('–', None, 'not applicable (positions have one threshold; the omnibus test gives no single '
                              'OR)')]),
             ab('GEC', 'OR', 'SPA', 'df', 'λGC')]
    return finish(sh, notes, 'A3 Q1 Association')


def three_cols(fields, sig=None):
    cols = []
    for c in COHORTS:
        for k, leaf, kind, w in fields:
            cols.append(C(f'{k}_{c}', leaf, kind, w, blk=COH_LABEL[c], coh=c,
                          hi=(lambda r, c=c: bool(r.get(f'sig_{c}'))) if (sig and k == 'p') else None))
    return cols


def step_rows(B, df):
    thr = scan_thresholds(B)
    order = ['allele', 'residue', 'position']
    rows = sorted(records(df), key=lambda r: (COMPARISONS.index(r['comparison']), order.index(r['test_type']), r['step']))
    for r in rows:
        r['tt'] = TT_MODEL[r['test_type']]
        r['cond'] = mlist(r['conditioned_on'], ' + ') if r['conditioned_on'] else 'nothing'
        r['outcome'] = 'Entered' if r['entered'] in (True, 'True') else \
            'Stop: ' + str(r['reason']).replace('no unit significant', 'no position significant')
        for c in COHORTS:
            r[f'sig_{c}'] = not blank(r[f'p_{c}']) and r[f'p_{c}'] < thr[(r['comparison'], c, r['test_type'])]
        if r['test_type'] == 'position':
            for c in COHORTS:
                r[f'or_{c}'] = NA
        r['_key'] = f"{r['comparison']}|{r['test_type']}|{r['step']}"
    return rows


def scan_thresholds(B):
    """(comparison, cohort, test type) -> the scan's primary threshold."""
    sc = B.s('a3_scans.tsv')
    return {(r.comparison, r.cohort, r.test_type): r.threshold_primary for r in sc.itertuples()}


STEP_COLS = lambda: ([C('tt', 'Test type', 'text', 30), C('step', 'Step', 'int', 5), C('cond', 'Given', 'md', 16),  # noqa: E731
                      C('lead', 'Strongest marker', 'marker', 14)]
                     + three_cols([('p', 'P', 'p', 12), ('or', 'OR', 'f2', 7)], sig=True)
                     + [C('n3', 'Significant in', 'text', 9), C('outcome', 'Outcome', 'text', 46)])



def a4(B):
    sh = B.sheet('A4 Q2 Conditioning', 'HLA association, Q2: how many independent signals are there?',
                 'Stepwise conditioning, all cases vs. controls · alleles and residues: random model, GEC · positions: '
                 'omnibus test given the whole of the positions entered, Bonferroni', 'Association')
    st = B.s('a4_steps.tsv')
    put(sh, f'The steps ({sa("resCond")})', STEP_COLS(), step_rows(B, st[st.comparison == 'all']), span=coh_span(), group='tt')
    af = B.s('a4_after.tsv')
    af = af[af.comparison == 'all']
    rows = []
    for k, lab in (('allele', 'Alleles, given DQB1 57 V'), ('residue', 'Residues, given DQB1 57 V (step 1)'),
                   ('position', 'Positions, given position 57 (step 1)')):
        r = {'tt': lab, '_key': k}
        for c in COHORTS:
            d = af[(af.test_type == k) & (af.cohort == c)].iloc[0]
            r[f'p_{c}'], r[f'top_{c}'], r[f'sig_{c}'] = d.min_p, d.top, int(d.sig)
        rows.append(r)
    put(sh, f'What is left once the signal is entered ({sa("resCond", "sumAnswers")})', [C('tt', 'Scan', 'text', 34)] + three_cols(
        [('p', 'Smallest P', 'p', 12), ('top', 'Marker', 'marker', 14), ('sig', 'Significant', 'int', 9)]), rows,
        span=coh_span())
    b = B.af['b114']
    put(sh, f'HLA-B position 114, the second strongest position, given DQB1 57 V ({sa("resCond")})', [
        C('cohort', 'Cohort', 'cohort', 20), C('before', 'Alone', 'p', 12, blk='P of position 114'),
        C('after', 'Given DQB1 57 V', 'p', 14, blk='P of position 114')],
        [{'cohort': c, **b[c]} for c in COHORTS], span=bspan('P of position 114'))
    notes = [('', ['Each step enters the marker with the smallest P in any cohort if it passes the threshold in at least '
                   'one; the next step tests every other marker given the markers entered. Shaded, P passes the threshold '
                   'in that cohort. No allele passes at step 0, so the allele walk stops there; the allele scan given DQB1 57 V '
                   'is the conditional scan of the summary slide. Positions are tested given the whole of position 57.']),
             ('', ['Given: for positions, the residues of the position entered enter as one block (the omnibus test given '
                   'the position).']),
             *missing_block([('–', None, 'not applicable (the omnibus test gives no single OR)')]),
             ab('GEC', 'OR')]
    return finish(sh, notes, 'A4 Q2 Conditioning')


def a5(B):
    sh = B.sheet('A5 Q3 GWAS lead', 'HLA association, Q3: does the HLA signal explain the genome-wide MHC peak?',
                 'The lead SNP of the genome-wide scan in the MHC, conditioned on HLA markers and the reverse · all '
                 'cases vs. controls unless stated', 'Association')
    ld = records(B.s('a5_leads.tsv'))
    for r in ld:
        r['status_l'] = 'Kept' if r['status'] == 'kept' else 'Not kept'
        r['reason_l'] = r['reason'].replace('P >= 1e-05', 'P ≥ 1 × 10⁻⁵').replace('P < 1e-05', 'P < 1 × 10⁻⁵') \
            .replace('narrow_mainland', 'narrow').replace('intermediate_mainland', 'intermediate') \
            .replace('full_mainland', 'full')
    put(sh, f'Candidate leads: MHC SNPs of the genome-wide scan ({sa("resGwas")})', [
        C('lead', 'SNP', 'text', 20), C('status_l', 'Lead', 'text', 8), C('reason_l', 'Rule', 'text', 16, note='a', span=3)]
        + three_cols([('p', 'P', 'p', 12), ('or', 'OR', 'f2', 7), ('freq', 'Frequency', 'pct1', 9)]),
        ld, span=coh_span())
    lk = records(B.s('a5_link.tsv'))
    cols = [C('cohort', 'Cohort', 'cohort', 20),
            C('lead_p', 'Alone', 'p', 12, blk='P of the GWAS lead', note='b'),
            C('lead_given_57v', 'Given DQB1 57 V', 'p', 12, blk='P of the GWAS lead'),
            C('lead_given_pos57', 'Given position 57', 'p', 12, blk='P of the GWAS lead'),
            C('tag', 'Best LD tag', 'marker', 13, blk='P of the GWAS lead', note='c'),
            C('tag_r2', 'r²', 'f2', 6, blk='P of the GWAS lead'),
            C('lead_given_tag', 'Given the tag', 'p', 12, blk='P of the GWAS lead'),
            C('r_57v', 'r, lead and 57 V', 'f2', 8),
            C('v57_p', 'Alone', 'p', 12, blk='P of DQB1 57 V'), C('v57_given_lead', 'Given the lead', 'p', 12, blk='P of DQB1 57 V'),
            C('pos57_p', 'Alone', 'p', 12, blk='P of position 57'),
            C('pos57_given_lead', 'Given the lead', 'p', 12, blk='P of position 57'),
            C('hla_min_given_lead', 'Smallest P', 'p', 12, blk='Every HLA marker given the lead', note='d'),
            C('hla_top_given_lead', 'Marker', 'marker', 13, blk='Every HLA marker given the lead'),
            C('most_explaining', 'Explains the lead most', 'marker', 12, note='e')]
    put(sh, f'The lead and DQB1 57 V remove each other\'s evidence ({sa("resGwas")})', cols, lk,
        span=bspan('P of the GWAS lead', 'P of DQB1 57 V', 'P of position 57', 'Every HLA marker given the lead'))
    gl = B.s('a5_given_lead.tsv')
    rows = records(gl[gl.comparison == 'all'])          # Q3 is asked of the main comparison, as on the GWAS-lead slide
    for r in rows:
        r['tt'] = TT_MODEL[r['test_type']]
        r['_key'] = f"{r['comparison']}|{r['cohort']}|{r['test_type']}"
    cols = [C('cohort', 'Cohort', 'cohort', 20), C('tt', 'Test type', 'text', 30),
            C('tests', 'Tests', 'int', 7), C('threshold', 'Threshold', 'p', 12),
            C('sig_before', 'Alone', 'int', 7, blk='Significant'),
            C('sig_given', 'Given the lead', 'int', 9, blk='Significant'),
            C('min_p_given', 'P', 'p', 12, blk='Smallest P given the lead'),
            C('top_given', 'Marker', 'marker', 13, blk='Smallest P given the lead'),
            C('lambda_before', 'Alone', 'f2', 7, blk='λGC'), C('lambda_given', 'Given the lead', 'f2', 9, blk='λGC')]
    put(sh, f'Every scan given the GWAS lead, all cases vs. controls ({sa("resGwas")})', cols, rows,
        span=bspan('Significant', 'Smallest P given the lead', 'λGC'), group='cohort')
    lg = records(B.s('a5_lead_given.tsv'))
    WHY = {'ld_tag': 'LD tag', 'top_p': 'strongest marker', 'signal_lead': 'signal lead',
           'signal_leads_jointly': 'all signal leads together', 'residue_lead': 'residue signal',
           'omnibus_lead': 'position signal'}
    for r in lg:
        w = str(r['why'])
        if w.startswith('position_jointly:'):
            r['why_l'] = 'whole position of the ' + ' and '.join(WHY[x] for x in w.split(':', 1)[1].split(','))
            r['cond_l'] = mlab(r['position']) + ' (whole position)'
            r['r'] = r['r2'] = NA
        else:
            r['why_l'] = '; '.join(WHY[x] for x in w.split(','))
            r['cond_l'] = mlist(r['conditioned_on'], ' + ')
        if w == 'signal_leads_jointly' or w.startswith('position_jointly:'):
            r['r'] = r['r2'] = r['marker_p_marginal'] = r['marker_p_given_lead'] = NA
        r['_key'] = f"{r['comparison']}|{r['cohort']}|{r['conditioned_on']}"
    cols = [C('comparison', 'Comparison', 'comp', 30), C('cohort', 'Cohort', 'cohort', 20),
            C('cond_l', 'HLA marker', 'md', 30), C('why_l', 'Chosen as', 'text', 16, note='f', span=2),
            C('r', 'r with the lead', 'f3', 8), C('r2', 'r²', 'f3', 7),
            C('lead_p_marginal', 'Alone', 'p', 12, blk='P of the GWAS lead'),
            C('lead_p_given', 'Given the marker', 'p', 12, blk='P of the GWAS lead'),
            C('lead_or_given', 'OR given the marker', 'f2', 9, blk='P of the GWAS lead'),
            C('marker_p_marginal', 'Alone', 'p', 12, blk='P of the HLA marker'),
            C('marker_p_given_lead', 'Given the lead', 'p', 12, blk='P of the HLA marker')]
    put(sh, 'The GWAS lead given each of its LD tags and strongest HLA markers (random model; not tabulated on the '
            'slides)', cols, lg, span=bspan('P of the GWAS lead', 'P of the HLA marker'), group='cohort')
    notes = [('a', ['Rule: a lead is an MHC SNP of the genome-wide scan with P < 1 × 10⁻⁵ in all three cohorts; one SNP '
                    'qualifies.']),
             ('b', ['The lead retested in the HLA run (random model, same samples and covariates); given a marker, with the '
                    'marker in the model; given position 57, with every residue column of the position but the most '
                    'frequent. Position 57 alone and given the lead: the fixed-model omnibus test.']),
             ('c', ['Best LD tag, the HLA marker with the largest r² with the lead in that cohort.']),
             ('d', ['Every HLA marker given the lead: alleles, residues and positions together, smallest P of all.']),
             ('e', ['Explains the lead most, the HLA marker that raises the lead\'s P the most when conditioned on.']),
             ('f', ['Chosen as: the ten HLA markers with the largest |r| with the lead (LD tag), the ten with the smallest '
                    'P alone (strongest marker), and the signal leads of the comparison, one at a time and together; whole '
                    'positions enter with every residue column but the most frequent.']),
             *missing_block([('–', 'r, r², P of the HLA marker', 'not applicable (several columns conditioned on '
                              'together)')]),
             ab('DVT', 'GWAS', 'LD', 'MHC', 'OR', 'SAIGE', 'SNP', 'λGC')]
    return finish(sh, notes, 'A5 Q3 GWAS lead')


LAM_ORDER = ['allele', 'residue', 'position', 'position, 1 df', 'position, ≥ 2 df', 'GWAS SNPs in the MHC',
             'GWAS SNPs outside the MHC']
LAM_LABEL = {'allele': 'Alleles (random model)', 'residue': 'Residues (random model)',
             'position': 'Positions (fixed model)', 'position, 1 df': 'Positions with 1 df',
             'position, ≥ 2 df': 'Positions with ≥ 2 df', 'GWAS SNPs in the MHC': 'Genome-wide scan, SNPs in the MHC',
             'GWAS SNPs outside the MHC': 'Genome-wide scan, SNPs outside the MHC'}
COND_ORDER = ['none', 'GWAS lead', 'DQB1 57 V', 'position 57']
COND_LABEL = {'none': 'Nothing', 'GWAS lead': 'The GWAS lead', 'DQB1 57 V': 'DQB1 57 V', 'position 57': 'Position 57'}


def a6(B):
    sh = B.sheet('A6 Inflation', 'HLA association: λGC of the HLA scans before and after conditioning',
                 'All cases vs. controls · λGC = median χ² over its expectation, over tests, each P read as 1 df',
                 'Association')
    d = B.s('a6_lambda.tsv')
    rows = []
    for t in LAM_ORDER:
        for k in COND_ORDER:
            x = d[(d.test_type == t) & (d.conditioning == k)]
            if x.empty:
                continue
            r = {'tt': LAM_LABEL[t], 'cond': COND_LABEL[k], '_key': f'{t}|{k}'}
            for c in COHORTS:
                y = x[x.cohort == c].iloc[0]
                r[f'tests_{c}'], r[f'before_{c}'] = y.tests, y.before
                r[f'after_{c}'] = NA if k == 'none' else y.after
            rows.append(r)
    cols = [C('tt', 'Scan', 'text', 34), C('cond', 'Given', 'text', 14)] + three_cols(
        [('tests', 'Tests', 'int', 9), ('before', 'Alone', 'f2', 7), ('after', 'Given', 'f2', 7)])
    put(sh, f'λGC per scan and conditioning ({sa("resQQ", "resLambda")})', cols, rows, span=coh_span(), group='tt')
    m = B.af['positions_multi_df']
    notes = [('', ['Alone is computed over the tests that remain after conditioning (the conditioned marker and its mirror '
                   'leave the scan), so it can differ slightly from the scan\'s own λGC. HLA markers are a correlated '
                   'region around a known signal, so λGC above 1 before conditioning is expected. Positions with ≥ 2 df: '
                   f"{m['n']} in full mainland, {m['in_dr_dq']} of them in DRB1, DQA1 or DQB1; reading their P as 1 df "
                   'inflates λGC.']),
             ('', ['Genome-wide scan: the common-variant scan of the same samples (MAF ≥ 0.5%); the MHC is 29.6–33.4 Mb '
                   'on chromosome 6. Outside the MHC, λGC near 1 shows the scan itself is not inflated.']),
             *missing_block([('–', None, 'not applicable (nothing conditioned on)')]),
             ab('GWAS', 'MAF', 'MHC', 'SNP', 'df', 'λGC')]
    return finish(sh, notes, 'A6 Inflation')


def a7(B):
    sh = B.sheet('A7 Summary', 'HLA association: the main finding, DQB1 57 V',
                 'All cases vs. controls · one protective signal at DQB1 position 57, in every cohort and in both models',
                 'Association')
    v = records(B.s('a7_v57.tsv'))
    for r in v:
        r['ci_r'], r['ci_f'] = ci(r), ci(r, 'or_fixed', 'l95_fixed', 'u95_fixed')
        r['chi_df'] = f"{F.fix(r['chisq_pos57'], 1)} ({int(r['df_pos57'])})"
    cols = [C('cohort', 'Cohort', 'cohort', 20),
            C('ci_r', 'OR (95% CI)', 'text', 18, blk='Random model', sub='DQB1 57 V', align='center'),
            C('p', 'P', 'p', 12, blk='Random model', sub='DQB1 57 V', hi=lambda r: r['p'] < r['gec']),
            C('gec', 'GEC', 'p', 12, blk='Random model', sub='Threshold'),
            C('bonf', 'Bonferroni', 'p', 12, blk='Random model', sub='Threshold'),
            C('p_spa', 'SPA P', 'p', 12, blk='Random model', note='a'),
            C('ci_f', 'OR (95% CI)', 'text', 18, blk='Fixed model', sub='DQB1 57 V', align='center'),
            C('p_fixed', 'P', 'p', 12, blk='Fixed model', sub='DQB1 57 V'),
            C('chi_df', 'χ² (df)', 'text', 10, blk='Fixed model', sub='Position 57', align='center'),
            C('p_pos57', 'P', 'p', 12, blk='Fixed model', sub='Position 57', hi=lambda r: r['p_pos57'] < r['threshold_pos57']),
            C('threshold_pos57', 'Bonferroni', 'p', 12, blk='Fixed model', sub='Position 57')]
    put(sh, f'DQB1 57 V and position 57 in both models ({sa("sumFind")})', cols, v, span=bspan('Random model', 'Fixed model'))
    g = records(B.s('a7_genotypes.tsv'))
    for r in g:
        r['grp'] = r['group'].capitalize()
    put(sh, f'Copies of valine at DQB1 57 per person ({sa("sumFind")})', [
        C('cohort', 'Cohort', 'cohort', 20), C('grp', 'Group', 'text', 9), C('analyzed', 'People', 'int', 8),
        C('undetermined', 'Undetermined', 'int', 10, note='b'), C('n0', '0', 'int', 7, blk='Copies'),
        C('n1', '1', 'int', 7, blk='Copies'), C('n2', '2', 'int', 7, blk='Copies'),
        C('carriers', 'Carriers', 'int', 8), C('copy_freq', 'Copy frequency', 'pct1', 9)], g,
        span=bspan('Copies'), group='cohort')
    rs = records(B.s('a7_residues.tsv'))
    for r in rs:
        r['ci'] = ci(r)
        r['name_l'] = r['name']
        pooled = 'rare allele' in r['name']
        r['name_l'] = re.sub(r'^1 rare alleles$', '1 rare allele', r['name'])
        r['what'] = 'Residue' if r['level'] == 'residue' else ('Alleles' if pooled else 'Allele')
        if pooled:
            r['ci'] = r['p'] = NA
        r['_key'] = f"{r['cohort']}|{r['name']}"
    cols = [C('cohort', 'Cohort', 'cohort', 20), C('residue', 'Residue at 57', 'text', 8, align='center'),
            C('what', 'Level', 'text', 8), C('name_l', 'Marker', 'marker', 13),
            C('n_alleles', 'Alleles', 'int', 7, note='c'), C('copies_case', 'Copies', 'int', 7, blk='Cases'),
            C('freq_case', 'Frequency', 'pct1', 9, blk='Cases'), C('copies_control', 'Copies', 'int', 7, blk='Controls'),
            C('freq_control', 'Frequency', 'pct1', 9, blk='Controls'),
            C('ci', 'OR (95% CI)', 'text', 18, align='center', note='d'), C('p', 'P', 'p', 12, note='d')]
    put(sh, f'Each residue at position 57 and its alleles ({sa("sumAlleles")})', cols, rs, span=bspan('Cases', 'Controls'),
        group='residue', bold=lambda r: r['level'] == 'residue')
    pv_ = B.af['panels_v57']
    put(sh, f'Valine at DQB1 57 in our groups and in the population panels ({sa("sumAlleles", "panels57groups")})', [
        C('src', 'Source', 'text', 34), C('chr', 'Copies typed', 'int', 10), C('f', 'Frequency of valine', 'pct1', 11)],
        [{'src': 'Our cases, full mainland', 'chr': pv_['cases_chr'], 'f': pv_['cases']},
         {'src': 'Our controls, full mainland', 'chr': pv_['controls_chr'], 'f': pv_['controls']},
         {'src': 'ToMMo 61KJPN (Tohoku)', 'chr': pv_['tommo_chr'], 'f': pv_['tommo']},
         {'src': '1KG JPT (Tokyo)', 'chr': pv_['kg_chr'], 'f': pv_['kg']}])
    # the three answers as on the summary slide, with the numbers behind them
    v3 = {r['cohort']: r for r in v}
    lk = {r['cohort']: r for r in records(B.s('a5_link.tsv'))}
    af = B.s('a4_after.tsv')
    af = af[af.comparison == 'all']
    left = {c: float(af[af.cohort == c].min_p.min()) for c in COHORTS}
    nif = lambda f: ' / '.join(f(c) for c in COHORTS)  # noqa: E731
    car = ' and '.join(B.af['carrier_alleles_57v'])
    ans = [{'q': 'Q1  Association', 'ask': 'Which alleles, residues, or positions are associated?',
            'a': 'DQB1 57 V, and position 57 as a whole, in all three cohorts — and no allele',
            'n': f"*DQB1* 57 V: P {nif(lambda c: F.pv(v3[c]['p']))} (GEC {nif(lambda c: F.pv(v3[c]['gec']))}); passes "
                 f"conventional Bonferroni in narrow mainland only. Position 57: P "
                 f"{nif(lambda c: F.pv(v3[c]['p_pos57']))}. No allele passes; the tested alleles with valine are {car}."},
           {'q': 'Q2  Independence', 'ask': 'How many independent signals are there?',
            'a': 'One. Conditioning on DQB1 57 V leaves no allele, residue, or position significant',
            'n': f"Smallest P left (alleles, residues, positions together): {nif(lambda c: F.pv(left[c]))}."},
           {'q': 'Q3  The GWAS signal', 'ask': 'Does the HLA signal explain the genome-wide MHC peak?',
            'a': 'Yes. The GWAS lead and DQB1 57 V remove each other\'s evidence; the lead\'s best LD tag does not',
            'n': f"The lead: P {nif(lambda c: F.pv(lk[c]['lead_p']))}, given *DQB1* 57 V "
                 f"{nif(lambda c: F.pv(lk[c]['lead_given_57v']))}; *DQB1* 57 V given the lead "
                 f"{nif(lambda c: F.pv(lk[c]['v57_given_lead']))}."}]
    put(sh, f'The three main questions ({sa("sumAnswers")})', [
        C('q', 'Question', 'text', 20, wrap=True), C('ask', '', 'text', 14, span=3, wrap=True),
        C('a', 'Answer', 'text', 14, span=3, wrap=True), C('n', 'Numbers (N / I / F)', 'md', 12, span=4, wrap=True)],
        ans)
    notes = [('a', ['SPA P, the saddlepoint-approximated P of the same score test (sensitivity analysis).']),
             ('b', ['Undetermined, the residue at 57 not determined for one of the person\'s alleles; those people enter '
                    'the random model with the mean copy number and leave the fixed model.']),
             ('c', ['Alleles, two-field alleles with this residue at 57; a residue row pools them all, an allele row is '
                    'one tested allele column.']),
             ('d', ['Random model, per copy, each marker alone; alleles below the MAF floor have no test of their own '
                    '(pooled rows, and named alleles without a value).']),
             ('', ['Caveats, as on the summary slide: the cohorts are nested and share most cases, so three of three shows '
                   'robustness to the ancestry definition, not replication; sequencing platform and depth coincide with '
                   'case status (sheet A10). ToMMo and 1KG JPT sample other regions of Japan (Tohoku, Tokyo).']),
             *missing_block([('–', 'Each residue at position 57', 'pooled alleles have no single test'),
                             ('—', 'Each residue at position 57', 'the allele is not a tested column (below the MAF '
                              'floor)')]),
             ab('1KG JPT', 'CI', 'GEC', 'GWAS', 'LD', 'MAF', 'MHC', 'N / I / F', 'OR', 'SPA', 'ToMMo', 'df')]
    return finish(sh, notes, 'A7 Summary')


G57_MARKERS = [('DQB1:87:Y', 'DQB1 87 Y (residue)', 'residue'), ('DQB1:125', 'DQB1 125 (whole position)', 'position'),
               ('DQB1:87', 'DQB1 87 (whole position)', 'position'), ('DQB1:57', 'DQB1 57 (whole position)', 'position'),
               ('DPB1*05:01', 'DPB1*05:01 (allele)', 'allele'), ('DQB1*05:03', 'DQB1*05:03 (allele)', 'allele'),
               ('DQB1*05:03|all', 'DQB1*05:03 (allele), all cases vs. controls', 'allele')]


def a8(B):
    sh = B.sheet('A8 DVT', 'HLA association, supplementary: a history of DVT',
                 'With DVT vs. controls, without DVT vs. controls, and with vs. without DVT among cases · the same '
                 'tests and thresholds as the main analysis · a pre-specified supplementary analysis', 'Association')
    s = records(B.s('a8_strata.tsv'))
    for r in s:
        r['ci'] = ci(r)
    put(sh, f'DQB1 57 V in the four comparisons ({sa("suppDvt")})', [
        C('comparison', 'Comparison', 'comp', 30), C('cohort', 'Cohort', 'cohort', 20),
        C('n_case', 'Case side', 'int', 8, blk='People'), C('n_control', 'Control side', 'int', 8, blk='People'),
        C('freq_case', 'Case side', 'pct1', 8, blk='Frequency of valine'),
        C('freq_control', 'Control side', 'pct1', 8, blk='Frequency of valine'),
        C('ci', 'OR (95% CI)', 'text', 18, align='center'), C('p', 'P', 'p', 12, hi=lambda r: r['p'] < r['gec']),
        C('gec', 'GEC', 'p', 12)], s, span=bspan('People', 'Frequency of valine'), group='comparison',
        bold=lambda r: r['comparison'] == 'all')
    sc = B.s('a3_scans.tsv')
    sc = sc[sc.comparison != 'all']
    cols = [C('comparison', 'Comparison', 'comp', 30), C('cohort', 'Cohort', 'cohort', 20)] + SCAN_COLS()
    put(sh, f'Every scan of the three DVT comparisons ({sa("suppNoDvtScan", "suppDvtYes", "suppDvtOnly")})', cols, scan_rows(sc),
        span=bspan('Primary', 'Bonferroni (sensitivity)', 'Smallest P'), group='comparison')
    st = B.s('a4_steps.tsv')
    put(sh, f'Stepwise conditioning in the DVT comparisons ({sa("suppNoDvtWhat", "suppDvtYes", "suppDvtOnly")})',
        [C('comparison', 'Comparison', 'comp', 30)] + STEP_COLS(), step_rows(B, st[st.comparison != 'all']),
        span=coh_span(), group='comparison')
    # without DVT: the HLA-DQ signal
    npos = records(B.s('a8_nodvt_positions.tsv'))
    for r in npos:
        r['unit_l'] = r['unit']
    put(sh, f'Without DVT vs. controls: the positions that pass Bonferroni in all three cohorts ({sa("suppNoDvtScan")})', [
        C('unit_l', 'Position', 'marker', 12)] + cohort_cols([('p', 'P', 'p', 12)]), npos, span=coh_span())
    lr = records(B.s('a8_lrt_vs_score.tsv'))
    put(sh, 'Without DVT vs. controls: DQB1 −27 A, −9 I, −5 S (one pattern), fixed-model LRT against the random-model '
            f'score test ({sa("suppNoDvtWhat")})', [
        C('cohort', 'Cohort', 'cohort', 20), C('p_lrt_fixed', 'Fixed model, LRT', 'p', 13, blk='P'),
        C('p_score_random', 'Random model, score test', 'p', 14, blk='P'), C('ratio', 'Ratio', 'f1', 7, note='d')], lr,
        span=bspan('P'))
    nd = records(B.s('a8_nodvt_dqb1.tsv'))
    RNAME = {'v_other': 'Other alleles with V at 57', 'y_other': 'Other alleles with Y at 87',
             'F': 'Alleles with F at 87', 'L': 'Alleles with L at 87', 'unresolved': 'Alleles without 57 or 87'}
    for r in nd:
        named = '*' in r['kind']
        r['kind_l'] = r['kind'] if named else RNAME[r['kind']]
        if not named:
            r['or'] = r['p'] = NA
        r['_key'] = r['kind']
    no125 = {k: sum(int(r[f'no125_{k}']) for r in nd) for k in ('alleles', 'case', 'control')}
    put(sh, f'Without DVT vs. controls, full mainland: DQB1 alleles by their residues at 57, 87 and 125 ({sa("suppNoDvtWhat")})', [
        C('kind_l', 'Allele or group', 'text', 30), C('n_alleles', 'Alleles', 'int', 7),
        C('res_57', '57', 'text', 6, blk='Residue at', align='center'), C('res_87', '87', 'text', 6, blk='Residue at', align='center'),
        C('res_125', '125', 'text', 6, blk='Residue at', align='center'),
        C('copies_case', 'Cases', 'int', 7, blk='Copies'), C('copies_control', 'Controls', 'int', 8, blk='Copies'),
        C('freq_case', 'Cases', 'pct1', 8, blk='Frequency'), C('freq_control', 'Controls', 'pct1', 8, blk='Frequency'),
        C('or', 'OR', 'f2', 7, note='e'), C('p', 'P', 'p', 12, note='e')], nd,
        span=bspan('Residue at', 'Copies', 'Frequency'))
    gy = records(B.s('a8_57v_given_87y.tsv'))
    put(sh, 'Without DVT vs. controls: DQB1 57 V given DQB1 87 Y (fixed model; not on the slide face)', [
        C('cohort', 'Cohort', 'cohort', 20), C('p', 'Alone', 'p', 12, blk='P of DQB1 57 V'),
        C('p_given_87y', 'Given 87 Y', 'p', 12, blk='P of DQB1 57 V'), C('or_given_87y', 'OR given 87 Y', 'f2', 9)], gy,
        span=bspan('P of DQB1 57 V'))
    # without DVT given DQB1 57 V (supplementary run)
    g = B.s('a8_given57v.tsv').set_index('cohort')
    rows = []
    for k in ('allele', 'residue', 'position'):
        r = {'tt': TT_MODEL[k], '_key': k}
        for c in COHORTS:
            x = g.loc[c]
            r[f'before_{c}'], r[f'after_{c}'] = int(x[f'{k}_before']), int(x[f'{k}_after'])
            r[f'minp_{c}'], r[f'top_{c}'] = x[f'{k}_min_p'], x[f'{k}_top']
            r[f'thr_{c}'] = x[f'{k}_gec'] if k != 'position' else x['position_threshold']
        rows.append(r)
    put(sh, f'Without DVT vs. controls, given DQB1 57 V: passing before and after ({sa("suppNoDvtNew")})',
        [C('tt', 'Test type', 'text', 30)] + three_cols([('before', 'Alone', 'int', 6), ('after', 'Given', 'int', 6),
                                                         ('thr', 'Threshold', 'p', 12), ('minp', 'Smallest P given', 'p', 12),
                                                         ('top', 'Marker', 'marker', 12)]), rows, span=coh_span())
    rows = []
    for m, lab, k in G57_MARKERS:
        r = {'m': mlab(m.split('|')[0]) + ' (' + lab.split(' (', 1)[1], '_key': m}
        for c in COHORTS:
            x = g.loc[c]
            r[f'p_{c}'], r[f'pg_{c}'] = x[f'{m}|p'], x[f'{m}|p_given']
            r[f'or_{c}'] = x[f'{m}|or'] if k != 'position' else NA
            r[f'org_{c}'] = x[f'{m}|or_given'] if k != 'position' else NA
        rows.append(r)
    put(sh, f'Without DVT vs. controls, given DQB1 57 V: the leads of {sa("suppNoDvtWhat")} ({sa("suppNoDvtNew")})', [C('m', 'Marker', 'md', 40)]
        + three_cols([('p', 'P alone', 'p', 12), ('pg', 'P given', 'p', 12), ('or', 'OR alone', 'f2', 7),
                      ('org', 'OR given', 'f2', 7)]), rows, span=coh_span())
    sl = records(B.s('a8_spa_lost.tsv'))
    put(sh, 'Calls that do not hold with the SPA P (all four comparisons)', [
        C('comparison', 'Comparison', 'comp', 30), C('cohort', 'Cohort', 'cohort', 20), C('marker', 'Marker', 'marker', 13),
        C('p', 'Score test', 'p', 12, blk='P'), C('p_spa', 'SPA', 'p', 12, blk='P'), C('threshold', 'GEC', 'p', 12)], sl,
        span=bspan('P'))
    ld = B.af['ld_57_87_without_dvt']
    spa = B.af['spa_calls']
    notes = [('a', ['Primary threshold: GEC for alleles and residues (random model); Bonferroni over the positions tested '
                    'for positions (fixed-model omnibus test).']),
             ('b', ['Primary, SPA P: markers passing the primary threshold with the saddlepoint-approximated P.']),
             ('c', ['λGC, median χ² of the scan\'s tests (each P read as one degree of freedom) over its expectation.']),
             ('d', ['Ratio, random-model P over fixed-model P; the fixed-model LRT is the more powerful test here, and the '
                    'random model is the one that calls single markers.']),
             ('e', ['OR and P of the allele column (random model, each alone); groups have no single test. Grouped rows '
                    f"list their members' residues; {no125['alleles']} rare members are not aligned at 125 "
                    f"({no125['control']} control and {no125['case']} case copies)."]),
             ('', [f"LD of 57 V and 87 Y without DVT (full mainland): r² {F.fix(ld['r2'], 3)} (its ceiling given the two "
                   f"frequencies {F.fix(ld['r2_ceiling'], 3)}), D′ {F.fix(ld['dprime'], 3)}; valine at 57 without tyrosine "
                   f"at 87 only on {', '.join(ld['v_not_y'])} ({ld['v_not_y_control']} control copies)."]),
             ('', ['Single-marker calls at GEC, all three cohorts together, with the score-test P / the SPA P: '
                   + '; '.join(f"{COMP_LABEL[k]} {v['calls']} / {v['calls_spa']}" for k, v in spa.items()) + '.']),
             ('', ['Given DQB1 57 V: a run of the association pipeline outside the component, with the pipeline\'s scripts '
                   'unchanged; every P alone equals the scan\'s P '
                   f"({B.af['given57v_marginal_rows_checked']:,} rows checked)."]),
             *missing_block([('–', None, 'not applicable (whole positions and allele groups give no single OR or P; '
                              'positions have one threshold)'),
                             ('—', 'DQB1 alleles by their residues', 'the allele is not a tested column in this comparison '
                              '(below the MAF floor)')]),
             ab('CI', 'DVT', 'GEC', 'LD', 'LRT', 'MAF', 'OR', 'SPA', 'λGC')]
    return finish(sh, notes, 'A8 DVT')


def a9(B):
    sh = B.sheet('A9 Reported alleles', 'HLA association, supplementary: the alleles reported in Japanese CTEPH',
                 'HLA-B*52:01 and DPB1*02:02 · the published results as printed, ours in the same comparisons',
                 'Association')
    pb = records(B.s('a9_published.tsv'))
    for r in pb:
        r['ci'] = ci(r, 'or_')
        r['car_c'] = f"{int(r['case_car'])} / {int(r['case_n'])}"
        r['car_k'] = f"{int(r['ctrl_car'])} / {int(r['ctrl_n'])}"
    put(sh, f'The published results ({sa("prior", "suppPrior")})', [
        C('study', 'Study', 'text', 14), C('set', 'Cases', 'text', 26), C('allele', 'Allele', 'text', 11),
        C('car_c', 'Cases', 'text', 9, blk='Carriers', align='center'), C('car_k', 'Controls', 'text', 9, blk='Carriers', align='center'),
        C('ci', 'OR (95% CI)', 'text', 17, align='center'), C('p', 'P', 'p', 12), C('pc', 'Corrected P', 'f4', 9, note='a'),
        C('case_freq_hwe', 'Cases', 'pct1', 8, blk='Copy frequency', note='b'),
        C('ctrl_freq_hwe', 'Controls', 'pct1', 8, blk='Copy frequency', note='b')], pb,
        span=bspan('Carriers', 'Copy frequency'))
    ou = records(B.s('a9_ours.tsv'))
    for r in ou:
        r['ci'] = ci(r)
    put(sh, f'Ours, per copy ({sa("suppPrior")})', [
        C('comparison', 'Comparison', 'comp', 30), C('cohort', 'Cohort', 'cohort', 20), C('allele', 'Allele', 'text', 11),
        C('freq_case', 'Case side', 'pct1', 8, blk='Frequency'), C('freq_control', 'Control side', 'pct1', 8, blk='Frequency'),
        C('ci', 'OR (95% CI)', 'text', 18, blk='Random model', align='center'), C('p', 'P', 'p', 12, blk='Random model'),
        C('or_fixed', 'OR', 'f2', 7, blk='Fixed model'), C('p_fixed', 'P', 'p', 12, blk='Fixed model')], ou,
        span=bspan('Frequency', 'Random model', 'Fixed model'), group='comparison')
    ca = records(B.s('a9_carriers.tsv'))
    for r in ca:
        r['ci'] = ci(r)
        r['car_c'] = f"{int(r['case_car'])} / {int(r['case_n'])}"
        r['car_k'] = f"{int(r['ctrl_car'])} / {int(r['ctrl_n'])}"
    put(sh, f'Ours, as carriers, the published measure ({sa("suppPrior")})', [
        C('comparison', 'Comparison', 'comp', 30), C('cohort', 'Cohort', 'cohort', 20), C('allele', 'Allele', 'text', 11),
        C('car_c', 'Case side', 'text', 10, blk='Carriers', align='center'),
        C('car_k', 'Control side', 'text', 11, blk='Carriers', align='center'),
        C('ci', 'OR (95% CI)', 'text', 18, align='center', note='c'), C('p', 'P', 'p', 12, note='c'),
        C('compared_with', 'Compared with', 'text', 14),
        C('power_published', 'Power for the published OR', 'pct1', 11, note='d'),
        C('detectable_or80', 'OR detectable at 80%', 'f2', 10, note='d')], ca, span=bspan('Carriers'), group='comparison')
    pr = B.af['pop_reported']
    rows = [{'allele': a, **{k: (NA if v is None else v) for k, v in d.items()}} for a, d in pr.items()]
    put(sh, f'The two alleles in our controls and in the population panels (P groups; {sa("suppPrior")})', [
        C('allele', 'Allele', 'text', 11), C('controls_pgroup', 'Our controls', 'pct1', 10, blk='Frequency'),
        C('tommo', 'ToMMo', 'pct1', 8, blk='Frequency'), C('kg', '1KG JPT', 'pct1', 8, blk='Frequency'),
        C('tommo_vs_controls_p', 'ToMMo vs. ours', 'p', 12, blk='P of the difference', note='e'),
        C('kg_vs_controls_p', '1KG JPT vs. ours', 'p', 12, blk='P of the difference', note='e')], rows,
        span=bspan('Frequency', 'P of the difference'))
    notes = [('a', ['Corrected P as published.']),
             ('b', ['Copy frequency implied by the published carrier frequency under Hardy–Weinberg, 1 − √(1 − carrier '
                    'frequency), for comparison with ours.']),
             ('c', ['Carriers vs. non-carriers on the unrelated samples, unadjusted: OR with Woolf\'s 95% CI, Fisher\'s '
                    'exact P.']),
             ('d', ['Power at α = 0.05 (two-sided) to detect the published carrier OR with our numbers, and the carrier OR '
                    'we could detect with 80% power.']),
             ('e', ['Two-sided test of equal frequency of the P group.']),
             ('', ['Tanabe 2005 compared CTEPH with or without DVT; Kominami 2009, CTEPH without DVT; ours is compared with '
                   'each in the matching comparison. Not replicated in any cohort.']),
             *missing_block([('–', None, '1KG JPT does not type DPB1')]),
             ab('1KG JPT', 'CI', 'CTEPH', 'DVT', 'OR', 'ToMMo')]
    return finish(sh, notes, 'A9 Reported alleles')


def a10(B):
    sh = B.sheet('A10 Appendix', 'HLA association, appendix: test counting, model agreement, platform, panels',
                 'Mirror pairs · fixed against random model · sequencing platform and case status · position 57 by P '
                 'group in ToMMo and 1KG JPT', 'Association')
    mi = records(B.s('a10_mirror.tsv'))
    put(sh, f'Mirror pairs of residue columns ({sa("tests")})', [
        C('comparison', 'Comparison', 'comp', 30), C('cohort', 'Cohort', 'cohort', 20),
        C('markers', 'Residue columns', 'int', 9), C('mirror_pairs', 'Mirror pairs', 'int', 8), C('tests', 'Tests', 'int', 7)],
        mi, group='comparison')
    co = records(B.s('a10_concordance.tsv'))
    for r in co:
        r['tt'] = TT_LABEL[r['test_type']]
    put(sh, f'Fixed against random model, single markers ({sa("fixed")})', [
        C('comparison', 'Comparison', 'comp', 30), C('cohort', 'Cohort', 'cohort', 20), C('tt', 'Test type', 'text', 9),
        C('tests', 'Tests', 'int', 7), C('r_nlp', '−log₁₀ P', 'f3', 8, blk='Correlation', note='a'),
        C('r_beta', 'log OR', 'f3', 8, blk='Correlation', note='a'),
        C('median_ratio', 'Median ratio of −log₁₀ P', 'f2', 10, note='a')], co, span=bspan('Correlation'), group='comparison')
    oo = records(B.s('a10_one_model_only.tsv'))
    for r in oo:
        r['tt'] = TT_LABEL[r['test_type']]
        r['in_'] = r['significant_in'].capitalize() + ' model only'
    put(sh, f'Markers significant in one model only ({sa("fixed")})', [
        C('comparison', 'Comparison', 'comp', 30), C('cohort', 'Cohort', 'cohort', 20), C('tt', 'Test type', 'text', 9),
        C('marker', 'Marker', 'marker', 13), C('in_', 'Passes GEC in', 'text', 16), C('p_fixed', 'Fixed', 'p', 12, blk='P'),
        C('p_random', 'Random', 'p', 12, blk='P'), C('gec', 'GEC', 'p', 12)], oo, span=bspan('P'), group='comparison')
    pl = B.s('a10_platform.tsv')
    rows = []
    for p in list(dict.fromkeys(pl.platform)):
        r = {'platform': p.replace('DNBseq', 'DNBSeq').replace('HiSeqX', 'HiSeq X').replace(' 30x', ' 30×')
             .replace(' 15x', ' 15×'), '_key': p}
        for c in COHORTS:
            x = pl[(pl.platform == p) & (pl.cohort == c)].iloc[0]
            for k in ('controls', 'cases_dvt_yes', 'cases_dvt_no', 'cases_dvt_unknown'):
                r[f'{k}_{c}'] = int(x[k])
        rows.append(r)
    put(sh, f'Sequencing platform by group ({sa("platform")})', [C('platform', 'Platform', 'text', 20)] + cohort_cols(
        [('controls', 'Controls', 'int', 8), ('cases_dvt_yes', 'Cases with DVT', 'int', 8),
         ('cases_dvt_no', 'Cases without DVT', 'int', 9), ('cases_dvt_unknown', 'Cases, DVT not recorded', 'int', 9)]),
        rows, span=coh_span())
    ch = records(B.s('a10_platform_chi2.tsv'))
    put(sh, f'Platform against history of DVT among cases with it recorded ({sa("platform")})', [
        C('cohort', 'Cohort', 'cohort', 20), C('chi2', 'χ²', 'f1', 8), C('df', 'df', 'int', 5), C('p', 'P', 'p', 12)], ch)
    pg = records(B.s('a10_pos57_pgroups.tsv'))
    for r in pg:
        r['ours'] = r['our_alleles'] or NA
        r['later'] = 'Yes' if r['named_later_than_tommo'] in (True, 'True') else 'No'
        r['_key'] = r['p_group']
    put(sh, f'Position 57 by P group: our cases and controls, ToMMo and 1KG JPT ({sa("panels57", "panels57groups")})', [
        C('residue_57', 'Residue at 57', 'text', 8, align='center'), C('p_group', 'P group', 'text', 14),
        C('ours', 'Our two-field alleles (copies)', 'text', 12, span=5), C('copies_case', 'Cases', 'int', 7, blk='Copies'),
        C('copies_control', 'Controls', 'int', 8, blk='Copies'), C('tommo', 'ToMMo', 'int', 8, blk='Copies'),
        C('kg', '1KG JPT', 'int', 8, blk='Copies'), C('later', 'Named after ToMMo', 'text', 9, note='b')], pg,
        span=bspan('Copies'), group='residue_57')
    mf, pt = B.af['mirror_full'], B.af['pos57_totals']
    notes = [('a', ['Correlation over the tests of the scan of −log₁₀ P and of log OR between the two models; median ratio, '
                    'the median over the tests of −log₁₀ P of the random model over that of the fixed model.']),
             ('b', ['Named after ToMMo: first named in IPD-IMGT/HLA after ToMMo\'s release '
                    f"({pt['tommo_vintage']}), so ToMMo could not report it."]),
             ('', [f"Mirror pairs, full mainland, all cases: {mf['pairs']} of the {mf['two_residue_positions']} "
                   f"two-residue positions; the other {mf['near_not_merged']} are complements in all but "
                   f"{mf['break_min']}–{mf['break_max']} people (an undetermined residue) and are tested separately. The two "
                   'columns of a pair give the same P and opposite log OR (largest difference in −log₁₀ P '
                   f"{F.fix(mf['max_dlog10p'], 1)})."]),
             ('', [f"Copies typed at position 57: cases {pt['cases_chr']:,}, controls {pt['controls_chr']:,}, ToMMo "
                   f"{pt['tommo_chr']:,}, 1KG JPT {pt['kg_chr']} ({pt['kg_people']} people; {pt['kg_missing']} copies "
                   'without a DQB1 call). Platform: every control was sequenced on HiSeq X at 15×, every case on DNBSeq or '
                   'NovaSeq; two checks are planned (platform-sensitivity analysis; down-sampling of depth).']),
             *missing_block([('–', 'Our two-field alleles', 'none of ours in the P group')]),
             ab('1KG JPT', 'DVT', 'GEC', 'IPD-IMGT/HLA', 'OR', 'ToMMo', 'df')]
    return finish(sh, notes, 'A10 Appendix')


def a11(B):
    sh = B.sheet('A11 Single markers (detail)', 'HLA association: every allele and residue test',
                 'Four comparisons × three cohorts · random model (primary) and fixed model side by side · P as numbers '
                 'for sorting · sortable (filter on the header row)', 'Association')
    d = B.s('a11_single_markers.tsv')
    rows = records(d)
    for r in rows:
        r['cls'] = TT_LABEL[r['class']]
        r['pos_l'] = NA if r['class'] == 'allele' else str(r['position']).replace('-', '−')
        r['res_l'] = NA if r['class'] == 'allele' else r['residue']
        r['mirror'] = r['mirror_of'] if r['mirror_of'] else NA
        r['f_note'] = {'.': NA, 'SEPARATION,ALT1': 'Separation', 'SEPARATION,REF': 'Separation'}.get(r['f_errcode'],
                                                                                                    r['f_errcode'])
        r['_key'] = f"{r['comparison']}|{r['cohort']}|{r['marker']}"
    cols = [C('comparison', 'Comparison', 'comp', 28), C('cohort', 'Cohort', 'cohort', 19), C('cls', 'Test type', 'text', 8),
            C('marker', 'Marker', 'marker', 13), C('gene', 'Gene', 'gene', 6), C('pos_l', 'Position', 'text', 7, align='center'),
            C('res_l', 'Residue', 'text', 6, align='center'), C('mirror', 'Mirror of', 'marker', 12, note='a'),
            C('n_case', 'Case side', 'int', 7, blk='People'), C('n_control', 'Control side', 'int', 7, blk='People'),
            C('copies_case', 'Case side', 'int', 7, blk='Copies'), C('copies_control', 'Control side', 'int', 7, blk='Copies'),
            C('freq_case', 'Case side', 'pct2', 8, blk='Frequency'), C('freq_control', 'Control side', 'pct2', 8, blk='Frequency'),
            C('mac', 'MAC', 'int', 6),
            C('missing_case', 'Case side', 'pct1', 7, blk='Undetermined'), C('missing_control', 'Control side', 'pct1', 7, blk='Undetermined'),
            C('r_or', 'OR', 'or', 7, blk='Random model (primary)'), C('r_l95', 'Lower 95%', 'or', 7, blk='Random model (primary)'),
            C('r_u95', 'Upper 95%', 'or', 7, blk='Random model (primary)'),
            C('r_p', 'P', 'pnum', 9, blk='Random model (primary)', hi=lambda r: r['r_sig_gec'] in (True, 'True')),
            C('r_p_spa', 'SPA P', 'pnum', 9, blk='Random model (primary)'),
            C('f_or', 'OR', 'or', 7, blk='Fixed model'), C('f_l95', 'Lower 95%', 'or', 7, blk='Fixed model'),
            C('f_u95', 'Upper 95%', 'or', 7, blk='Fixed model'), C('f_p', 'P', 'pnum', 9, blk='Fixed model'),
            C('f_note', 'Note', 'text', 9, blk='Fixed model', note='b'),
            C('gec', 'GEC', 'pnum', 9, blk='Threshold'), C('bonf', 'Bonferroni', 'pnum', 9, blk='Threshold'),
            C('r_sig_gec', 'GEC', 'bool', 6, blk='Significant, random model'),
            C('r_sig_bonf', 'Bonferroni', 'bool', 8, blk='Significant, random model'),
            C('f_sig_gec', 'Fixed model, GEC', 'bool', 8, note='c'),
            C('n3_gec', 'Cohorts significant, GEC', 'int', 9, note='d')]
    put(sh, None, cols, rows, span=bspan('People', 'Copies', 'Frequency', 'Undetermined', 'Random model (primary)',
                                         'Fixed model', 'Threshold', 'Significant, random model'),
        freeze=4, autofilter=True)
    notes = [('a', ['Mirror of: the column is the complement of the named one in every sample; the pair is one test (both '
                    'rows are listed, with the same P).']),
             ('b', ['Note: separation, the fixed-model estimate is at the boundary (a group without any copy); its P is not '
                    'trusted.']),
             ('c', ['Fixed model, GEC: the fixed-model P passes the same GEC threshold (sensitivity analysis).']),
             ('d', ['Cohorts significant, GEC: in how many of the three cohorts the marker passes GEC in this comparison.']),
             ('', ['People, those in the random model; the fixed model uses unrelated samples (sheet A1). Copies and '
                   'frequency are over called copies; undetermined, share of people without a call. Shaded P: passes GEC.']),
             *missing_block([('–', None, 'not applicable (alleles have no position or residue; no mirror; no note)'),
                             ('—', 'Fixed model', 'no estimate (the marker does not vary among the unrelated samples)')]),
             ab('GEC', 'MAC', 'OR', 'SPA')]
    return finish(sh, notes, 'A11 Single markers (detail)')


def a12(B):
    sh = B.sheet('A12 Positions (detail)', 'HLA association: every whole-position test',
                 'Fixed-model omnibus LRT, all residues of a position at once · four comparisons × three cohorts · P as '
                 'numbers for sorting · sortable (filter on the header row)', 'Association')
    d = B.s('a12_positions.tsv')
    rows = records(d)
    for r in rows:
        r['pos_l'] = str(r['position']).replace('-', '−')
        e = r['errcode']
        r['status'] = {'.': 'Tested', 'NOT_TESTED': 'Not tested', 'SEPARATION': 'Separation'}.get(e, e)
        r['others'] = (r['categories'] or '').replace('rare', 'pooled rare').replace(',', ', ') or None
        r['ors_l'] = (r['ors'] or '').replace(',', '; ').replace('rare', 'pooled rare').replace(':', ' ') or None
        r['_key'] = f"{r['comparison']}|{r['cohort']}|{r['unit']}"
        if e == 'NOT_TESTED':
            for k in ('n', 'cases', 'controls', 'df', 'chisq', 'p', 'others', 'ors_l', 'reference'):
                r[k] = NA
    cols = [C('comparison', 'Comparison', 'comp', 28), C('cohort', 'Cohort', 'cohort', 19), C('gene', 'Gene', 'gene', 6),
            C('pos_l', 'Position', 'text', 7, align='center'), C('status', 'Status', 'text', 10, note='a'),
            C('reference', 'Reference residue', 'text', 9, align='center', note='b'),
            C('others', 'Other residue groups', 'text', 16, note='b'),
            C('n', 'People', 'int', 7, blk='Unrelated samples'), C('cases', 'Case side', 'int', 7, blk='Unrelated samples'),
            C('controls', 'Control side', 'int', 7, blk='Unrelated samples'),
            C('df', 'df', 'int', 4, blk='Omnibus LRT'), C('chisq', 'χ²', 'f2', 7, blk='Omnibus LRT'),
            C('p', 'P', 'pnum', 9, blk='Omnibus LRT', hi=lambda r: r['sig_bonf'] in (True, 'True')),
            C('threshold', 'Bonferroni', 'pnum', 9), C('sig_bonf', 'Significant', 'bool', 8),
            C('n3', 'Cohorts significant', 'int', 9, note='c'), C('ors_l', 'OR per residue group', 'text', 30, note='d')]
    put(sh, None, cols, rows, span=bspan('Unrelated samples', 'Omnibus LRT'), freeze=4, autofilter=True)
    notes = [('a', ['Status: not tested, fewer than two residue groups once residues below 1% are pooled; separation, an '
                    'estimate at the boundary (P kept, not trusted, outside λGC).']),
             ('b', ['Reference, the most frequent residue; the other groups each enter the model with one copy column '
                    '("pooled rare", the residues below 1% together).']),
             ('c', ['Cohorts significant: in how many of the three cohorts the position passes Bonferroni in this '
                    'comparison.']),
             ('d', ['OR of each residue group against the reference, from the same model.']),
             ('', ['Shaded P: passes Bonferroni.']),
             *missing_block([('–', None, 'not applicable (position not tested)')]),
             ab('LRT', 'OR', 'df', 'λGC')]
    return finish(sh, notes, 'A12 Positions (detail)')


# ── Sources and checks ──────────────────────────────────────────────────────
def sources(B):
    sh = B.sheet('Sources and checks', 'Sources of every sheet and the checks of every record',
                 'For the analysts: the upstream runs, the files read, and each check with the records it compared · '
                 'the only sheet that names files', 'Sources')
    up = B.up
    rows = [{'what': r['component'], 'item': r['item'], 'value': r['value']} for r in up['runs']]
    put(sh, 'Upstream runs', [C('what', 'Component', 'text', 24), C('item', 'Item', 'text', 24),
                              C('value', 'Value', 'text', 12, span=10, wrap=True)], rows, group='what')
    rows = [{**r, 'rows_': r.get('rows')} for r in up['files']]
    put(sh, 'Files read', [C('sheets', 'Sheets', 'text', 24), C('role', 'Content', 'text', 24, wrap=True),
                           C('path', 'File', 'text', 12, span=7, wrap=True), C('modified', 'Modified', 'text', 12, span=2),
                           C('md5', 'MD5', 'text', 12, span=3)], rows, group='sheets')
    ck = B.checks['checks']
    rows = []
    for c in ck:
        rows.append({'group': c['group'], 'sheet': c['sheet'], 'what': c['what'], 'compared': c['compared'],
                     'mismatches': c['mismatches'], 'result': 'Passed' if c['passed'] else 'FAILED',
                     'route': c.get('route') or NA, '_key': c['id']})
    put(sh, 'Checks of the records (before the workbook is written)', [
        C('group', 'Check', 'text', 24), C('sheet', 'Sheets', 'text', 24, wrap=True),
        C('what', 'What is compared', 'text', 12, span=4, wrap=True), C('route', 'Independent route', 'text', 12, span=3,
                                                                        wrap=True),
        C('compared', 'Values compared', 'int', 12), C('mismatches', 'Mismatches', 'int', 12,
                                                         hi=lambda r: r['mismatches'] > 0),
        C('result', 'Result', 'text', 12)], rows, group='group')
    dk = B.checks.get('deck_files', [])
    if dk:
        put(sh, 'Deck data the agreement check read', [C('deck', 'Deck', 'text', 24), C('path', 'File', 'text', 12,
                                                                                         span=8, wrap=True),
                                                        C('md5', 'MD5', 'text', 12, span=3)], dk, group='deck')
    notes = [('', [('b', 'After writing. '), 'The workbook itself is re-read cell by cell and compared with the tables it '
                   'was written from (every number, every P string against the decks\' formatter), scanned for internal '
                   'names on the content sheets, for any sample identifier on every sheet, and for empty or undefined '
                   'cells; the result is in audit.json next to the workbook.']),
             ('', [('b', 'Not presented. '), 'The pipeline\'s multivariate-normal max-T, its Benjamini–Hochberg calls and '
                   'its gene-level omnibus tests are computed upstream but were dropped from the decks; they are not read.'])]
    return finish(sh, notes, 'Sources and checks')


SHEET_INFO.update({
    'T1 Typing run': ('setting; group × platform', 'the typing run; genomes, reads, depth, outliers, call quality by platform'),
    'T2 DRB3, DRB4, DRB5': ('gene', 'copy number of DRB3, DRB4, DRB5 against the DRB1 alleles; reconciliation'),
    'T3 Completeness': ('gene (19)', 'calls by result, basis of the allele, read evidence, confidence, locus tier'),
    'T4 ToMMo comparison': ('step; gene × panel', 'harmonization, presence in ToMMo, frequency agreement, DRB1*15:13, homozygosity'),
    'T5 Marker QC': ('step; gene; cohort × class', 'encoding, call-rate and MAF floors, stepwise record, columns per gene'),
    'T6 P groups': ('cohort; P group', 'tested columns that share a P group; DQA1 through the comparison'),
    'T7 Alleles (detail)': ('allele (1,995)', 'every allele column per cohort: copies, call rate, MAF, status, panels'),
    'T8 Residues (detail)': ('residue column (2,872)', 'every residue column per cohort: copies, call rate, MAF, status'),
    'A1 Design': ('comparison × cohort; gene', 'samples of both models, history of DVT, tests per gene, models and settings'),
    'A2 Thresholds': ('scan', 'Me, GEC and Bonferroni thresholds, exact max-T benchmark and FWER, Me recomputed'),
    'A3 Q1 Association': ('scan; marker', 'significant per scan in the main comparison; the five strongest of each test type'),
    'A4 Q2 Conditioning': ('step', 'stepwise conditioning; what is left given DQB1 57 V'),
    'A5 Q3 GWAS lead': ('cohort; scan; HLA marker', 'the GWAS lead and DQB1 57 V conditioned on each other; every scan given the lead'),
    'A6 Inflation': ('scan × conditioning', 'λGC before and after conditioning; genome-wide SNPs in and outside the MHC'),
    'A7 Summary': ('cohort; residue', 'DQB1 57 V in both models, genotypes, residues at 57, panels, answers to Q1–Q3'),
    'A8 DVT': ('comparison × cohort; marker', 'the three DVT comparisons; without DVT given DQB1 57 V'),
    'A9 Reported alleles': ('study; comparison × cohort', 'B*52:01 and DPB1*02:02: published, ours per copy and as carriers'),
    'A10 Appendix': ('scan; platform; P group', 'mirror pairs, model agreement, platform by group, position 57 by P group'),
    'A11 Single markers (detail)': ('test (9,844)', 'every allele and residue test, both models side by side'),
    'A12 Positions (detail)': ('position (3,552)', 'every whole-position test'),
})

ORDER = ['Contents', 'T1 Typing run', 'T2 DRB3, DRB4, DRB5', 'T3 Completeness', 'T4 ToMMo comparison', 'T5 Marker QC',
         'T6 P groups', 'T7 Alleles (detail)', 'T8 Residues (detail)', 'A1 Design', 'A2 Thresholds', 'A3 Q1 Association',
         'A4 Q2 Conditioning', 'A5 Q3 GWAS lead', 'A6 Inflation', 'A7 Summary', 'A8 DVT', 'A9 Reported alleles',
         'A10 Appendix', 'A11 Single markers (detail)', 'A12 Positions (detail)', 'Sources and checks']
BUILD = {'T1 Typing run': t1, 'T2 DRB3, DRB4, DRB5': t2, 'T3 Completeness': t3, 'T4 ToMMo comparison': t4,
         'T5 Marker QC': t5, 'T6 P groups': t6, 'T7 Alleles (detail)': t7, 'T8 Residues (detail)': t8, 'A1 Design': a1,
         'A2 Thresholds': a2, 'A3 Q1 Association': a3, 'A4 Q2 Conditioning': a4, 'A5 Q3 GWAS lead': a5,
         'A6 Inflation': a6, 'A7 Summary': a7, 'A8 DVT': a8, 'A9 Reported alleles': a9, 'A10 Appendix': a10,
         'A11 Single markers (detail)': a11, 'A12 Positions (detail)': a12, 'Sources and checks': sources}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for k in ('typing', 'assoc', 'decks', 'slide_map', 'upstream', 'out', 'meta', 'cells'):
        ap.add_argument('--' + k.replace('_', '-'), required=True)
    ap.add_argument('--checks', required=True, nargs='+', help='deck_check.json and recheck.json')
    a = ap.parse_args()
    B = Book(a)
    B.up = json.loads(Path(a.upstream).read_text())
    if not B.checks.get('all_passed'):
        raise SystemExit('checks.json: not all checks passed; the workbook is not written')
    B.sheet('Contents', 'HLA typing and HLA association in CTEPH: record of the analyses',
            'One workbook for both components, in the order of the two decks · every number shared with a slide equals '
            'it', 'Contents')
    for name in ORDER[1:]:
        BUILD[name](B)
    contents(B, ORDER)
    B.wb.close()
    # workbook.json: the layout (sheets, blocks, rows, columns); workbook_cells.json.gz: what every body cell shows
    meta = {'sheets': [{'name': s.name, 'part': s.part, 'title': s.ttl, 'last_row': s.last_row,
                        'blocks': [{k: v for k, v in b.items() if k not in ('cells', 'raw')} for b in s.blocks]}
                       for s in B.sheets]}
    Path(a.meta).write_text(json.dumps(meta, indent=1, ensure_ascii=False, default=str) + '\n')
    cells = {s.name: [{'first_row': b['first_row'], 'pcols': b['pcols'], 'cells': b['cells'], 'cols': b['cols'],
                       'row_keys': b['row_keys'], 'raw': b['raw']} for b in s.blocks] for s in B.sheets}
    with gzip.open(a.cells, 'wt') as f:
        json.dump(cells, f, ensure_ascii=False, default=lambda x: None if isinstance(x, float) else str(x))
    print(f'[write_workbook] {len(B.sheets)} sheets -> {a.out}')


if __name__ == '__main__':
    main()
