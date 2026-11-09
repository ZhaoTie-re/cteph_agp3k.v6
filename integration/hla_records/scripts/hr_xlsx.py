#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : House-style building blocks of the workbook, copied (not imported)
#           from integration/gene_signals/scripts/write_table.py so that a later
#           change there cannot change this workbook: an XlsxWriter sheet with a
#           format cache and rich text, the booktabs Table (spanner / sub-spanner /
#           leaf header rows, cohort colours, rules), notes as short lines below
#           the table, the 'Missing values' block and the Abbreviations line.
#           Added here: a generic block writer driven by column specs, and P values
#           in the decks' format (deckfmt.pv) with superscript exponents.
# Component: hla_records
# ---------------------------------------------------------------------------
import re

import numpy as np

import deckfmt as F

FONT = 'Arial'
INK, INK2, GREY, FAINT, RULE_GREY, NOTE = '#1A1A1A', '#595959', '#7F7F7F', '#A6A6A6', '#BFBFBF', '#404040'
HEADER_FILL, CALL_FILL = '#F2F2F2', '#FCEFD2'
THIN, MEDIUM, HAIR = 1, 2, 7
MISSING, NONE_ = '—', '–'
SUPMAP = str.maketrans('⁻⁰¹²³⁴⁵⁶⁷⁸⁹', '−0123456789')


def blank(x):
    return x is None or (isinstance(x, float) and not np.isfinite(x)) or (isinstance(x, str) and x == '')


def sup_parts(text):
    """'4.57 × 10⁻⁵' -> ['4.57 × 10', ('sup', '−5')]: Unicode superscripts become superscript runs."""
    out, buf, sup = [], '', ''
    for ch in str(text):
        if ch in '⁻⁰¹²³⁴⁵⁶⁷⁸⁹':
            if buf:
                out.append(buf)
                buf = ''
            sup += ch.translate(SUPMAP)
        else:
            if sup:
                out.append(('sup', sup))
                sup = ''
            buf += ch
    if buf:
        out.append(buf)
    if sup:
        out.append(('sup', sup))
    return out


def p_parts(p):
    return None if blank(p) else sup_parts(F.pv(float(p)))


def text_len(parts):
    return sum(len(p if isinstance(p, str) else p[1]) for p in parts)


def plain(v, sup=False):
    """Text of a cell read back with rich_text=True (sup=True: only its superscript runs)."""
    from openpyxl.cell.rich_text import CellRichText
    if isinstance(v, CellRichText):
        return ''.join((b if isinstance(b, str) else b.text) for b in v
                       if (not isinstance(b, str) and b.font is not None and b.font.vertAlign == 'superscript') == sup)
    return '' if (v is None or sup) else str(v)


def wrap_parts(parts, width):
    tokens = []
    for p in parts:
        st, tx = ('', p) if isinstance(p, str) else p
        words = tx.split(' ')
        for k, wd in enumerate(words):
            tokens.append((st, wd + (' ' if k < len(words) - 1 else '')))
    lines, cur, n = [], [], 0
    for st, tx in tokens:
        if tx == '':
            continue
        if n + len(tx.rstrip()) > width and cur and st not in ('sup', 'sub') and cur[-1][1].endswith(' '):
            lines.append(cur)
            cur, n = [], 0
        cur.append((st, tx))
        n += len(tx)
    if cur:
        lines.append(cur)
    out = []
    for ln in lines:
        ln[-1] = (ln[-1][0], ln[-1][1].rstrip())
        out.append([tx if st == '' else (st, tx) for st, tx in ln])
    return out


def md(text):
    """Light markup for notes and cells: *x* italic, **x** bold, ^x^ superscript."""
    parts = []
    for tok in re.split(r'(\*\*[^*]+\*\*|\*[^*]+\*|\^[^^]+\^)', str(text)):
        if not tok:
            continue
        if tok.startswith('**'):
            parts.append(('b', tok[2:-2]))
        elif tok.startswith('*'):
            parts.append(('i', tok[1:-1]))
        elif tok.startswith('^'):
            parts.append(('sup', tok[1:-1]))
        else:
            parts.append(tok)
    return parts


class Sheet:
    """An XlsxWriter worksheet with a format cache and a rich-text helper (gene_signals)."""

    def __init__(self, wb, ws, name):
        self.wb, self.ws, self.name, self._fmt = wb, ws, name, {}
        self.blocks = []

    def fmt(self, **kw):
        props = {'font_name': FONT, 'font_size': 9, 'font_color': INK, 'valign': 'vcenter'}
        props.update({k: v for k, v in kw.items() if v is not None})
        key = tuple(sorted(props.items()))
        if key not in self._fmt:
            self._fmt[key] = self.wb.add_format(props)
        return self._fmt[key]

    def font(self, style, b, i, color, sz):
        props = {'font_name': FONT, 'font_size': sz, 'font_color': color, 'bold': bool(b), 'italic': bool(i)}
        if style == 'i':
            props['italic'] = True
        elif style == 'b':
            props['bold'] = True
        elif style == 'bi':
            props['bold'], props['italic'] = True, True
        elif style == 'sup':
            props['font_script'] = 1
        elif style == 'sub':
            props['font_script'] = 2
        elif style.startswith('#'):
            props['font_color'], props['bold'] = style, True
        key = ('font',) + tuple(sorted(props.items()))
        if key not in self._fmt:
            self._fmt[key] = self.wb.add_format(props)
        return self._fmt[key]

    def rich(self, r, c, parts, cell, b=False, i=False, color=INK, sz=9):
        frags = []
        for p in parts:
            st, tx = ('', p) if isinstance(p, str) else p
            if tx == '':
                continue
            if frags and frags[-1][0] == st:
                frags[-1][1] += tx
            else:
                frags.append([st, tx])
        if not frags:
            self.ws.write_blank(r, c, None, self.fmt(**cell))
        elif len(frags) == 1:
            st, tx = frags[0]
            f = dict(cell, font_size=sz, font_color=st if st.startswith('#') else color,
                     bold=True if (b or st in ('b', 'bi') or st.startswith('#')) else None,
                     italic=True if (i or st in ('i', 'bi')) else None,
                     font_script=1 if st == 'sup' else 2 if st == 'sub' else None)
            self.ws.write_string(r, c, tx, self.fmt(**f))
        else:
            args = []
            for st, tx in frags:
                args += [self.font(st, b, i, color, sz), tx]
            self.ws.write_rich_string(r, c, *args, self.fmt(**cell))


class Table:
    """A booktabs table (gene_signals): spanner / sub-spanner / leaf header rows, then body rows with rules.
    spec: dicts {blk, sub, leaf, note, w, coh, vline}; span: {blk: (parts, colour, cohort or None)}."""

    def __init__(self, sh, spec, top, span, colors, tints, header_h=None, set_width=True):
        self.sh, self.spec, self.top, self.colors, self.tints = sh, spec, top, colors, tints
        self.ncol = len(spec)
        self.first_coh = {}
        for j, s in enumerate(spec):
            c = s.get('coh')
            if c and (j == 0 or spec[j - 1].get('coh') != c):
                self.first_coh[j] = c
        if set_width:
            for j, s in enumerate(spec):
                sh.ws.set_column(j, j, s['w'])
        self._header(span, header_h)
        self.row0 = top + 3

    def _hcell(self, j, row, **extra):
        s = self.spec[j]
        d = {'align': 'center', 'valign': 'bottom', 'text_wrap': True, 'bold': True}
        c = s.get('coh')
        if c:
            d.update({'bg_color': self.tints[c], 'pattern': 1})
            if j in self.first_coh:
                d.update({'left': MEDIUM, 'left_color': self.colors[c]})
        if s.get('vline'):
            d.update({'left': s['vline'], 'left_color': RULE_GREY})
        if row == self.top:
            d.update({'top': MEDIUM, 'top_color': INK, 'valign': 'vcenter'})
        if row == self.top + 2:
            d.update({'bottom': MEDIUM, 'bottom_color': INK})
        d.update(extra)
        return d

    def _span_rule(self, coh):
        return {'bottom': MEDIUM, 'bottom_color': self.colors[coh]} if coh else {'bottom': THIN, 'bottom_color': INK}

    def _leaf(self, j):
        s = self.spec[j]
        leaf = s['leaf'] if isinstance(s['leaf'], list) else [('b', s['leaf'])]
        return leaf + [('sup', s.get('note') or '')]

    def _header(self, span, header_h):
        sh, ws, spec, R0 = self.sh, self.sh.ws, self.spec, self.top
        j = 0
        while j < self.ncol:
            blk = spec[j]['blk']
            k = j
            while k + 1 < self.ncol and spec[k + 1]['blk'] == blk and blk is not None:
                k += 1
            if blk is None:
                # added here: a column may span several sheet columns ('cont' entries follow its first one)
                k = j
                while k + 1 < self.ncol and spec[k + 1].get('cont'):
                    k += 1
                f = self._hcell(j, R0, bottom=MEDIUM, bottom_color=INK, valign='bottom')
                ws.merge_range(R0, j, R0 + 2, k, '', sh.fmt(**f))
                sh.rich(R0, j, self._leaf(j), f, b=True)
                j = k + 1
                continue
            parts, color, coh_blk = span[blk]
            for jj in range(j, k + 1):
                ws.write_blank(R0, jj, None, sh.fmt(**self._hcell(jj, R0, **self._span_rule(coh_blk))))
            f = self._hcell(j, R0, **self._span_rule(coh_blk))
            if k > j:
                ws.merge_range(R0, j, R0, k, '', sh.fmt(**f))
            sh.rich(R0, j, parts, f, b=True, color=color, sz=10)
            m = j
            while m <= k:
                sub = spec[m]['sub']
                if sub is None:
                    f = self._hcell(m, R0 + 1)
                    ws.merge_range(R0 + 1, m, R0 + 2, m, '', sh.fmt(**f))
                    sh.rich(R0 + 1, m, self._leaf(m), f, b=True)
                    m += 1
                    continue
                n = m
                while n + 1 <= k and spec[n + 1]['sub'] == sub:
                    n += 1
                note = next((spec[q]['note'] for q in range(m, n + 1) if spec[q].get('note')), '')
                f = self._hcell(m, R0 + 1, bottom=THIN, bottom_color=INK)
                for q in range(m, n + 1):
                    ws.write_blank(R0 + 1, q, None, sh.fmt(**self._hcell(q, R0 + 1, bottom=THIN, bottom_color=INK)))
                if n > m:
                    ws.merge_range(R0 + 1, m, R0 + 1, n, '', sh.fmt(**f))
                sub_parts = sub if isinstance(sub, list) else [('b', sub)]
                sh.rich(R0 + 1, m, sub_parts + [('sup', note if n > m else '')], f, b=True)
                for q in range(m, n + 1):
                    leaf = spec[q]['leaf']
                    lp = leaf if isinstance(leaf, list) else [('b', leaf)]
                    sh.rich(R0 + 2, q, lp + ([('sup', spec[q].get('note') or '')] if n == m else []),
                            self._hcell(q, R0 + 2), b=True)
                m = n + 1
            j = k + 1
        ws.set_row(R0, 20)
        if header_h:
            ws.set_row(R0 + 2, header_h)

    def rules(self, j, d, group_top=False, soft_top=False, last=False):
        c = self.first_coh.get(j)
        if c:
            d.update({'left': MEDIUM, 'left_color': self.colors[c]})
        elif self.spec[j].get('vline'):
            d.update({'left': self.spec[j]['vline'], 'left_color': RULE_GREY})
        if group_top:
            d.update({'top': THIN, 'top_color': RULE_GREY})
        elif soft_top:
            d.update({'top': HAIR, 'top_color': RULE_GREY})
        if last:
            d.update({'bottom': MEDIUM, 'bottom_color': INK})
        return d


def write_notes(sh, row, notes, width):
    lens = []
    for letter, parts in notes:
        for k, ln in enumerate(wrap_parts(parts, width - 3)):
            lead = [('sup', letter), ' '] if (letter and k == 0) else []
            sh.rich(row, 0, lead + ln, {'align': 'left', 'indent': 1 if k else None}, color=NOTE)
            lens.append(text_len(lead + ln))
            row += 1
    return row, lens


def title(sh, text, sub):
    sh.rich(0, 0, md(text) if isinstance(text, str) else text, {'align': 'left'}, b=True, sz=11)
    sh.rich(1, 0, md(sub) if isinstance(sub, str) else sub, {'align': 'left'}, i=True, color=INK2)
    sh.ws.set_row(0, 18)


def new_sheet(wb, name, zoom=90):
    sh = Sheet(wb, wb.add_worksheet(name), name)
    sh.ws.hide_gridlines(2)
    sh.ws.set_zoom(zoom)
    sh.ws.set_default_row(12.75)
    return sh


def print_setup(ws, last_row, ncol, footer, repeat=None, landscape=True, paper=9, scale=None):
    if landscape:
        ws.set_landscape()
    ws.set_paper(paper)
    if scale:
        ws.set_print_scale(scale)
    else:
        ws.fit_to_pages(1, 0)
    if repeat:
        ws.repeat_rows(*repeat)
    ws.print_area(0, 0, last_row, ncol - 1)
    ws.set_margins(left=0.4, right=0.4)
    ws.set_footer(f'&L{footer}&RPage &P of &N')


def missing_block(items):
    out = [('', [('b', 'Missing values')])]
    for sym, where, text in items:
        out.append(('', ['   ', ('b', sym)] + ([f', in {where}: '] if where else [': ']) + md(text)))
    return out


def abbrev(d):
    return ('', ['Abbreviations: ' + '; '.join(f'{k}, {v}' for k, v in sorted(d.items(), key=lambda kv:
                                                                               kv[0].lower())) + '.'])
