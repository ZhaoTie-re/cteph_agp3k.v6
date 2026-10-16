#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The deliverable workbook, written with XlsxWriter in the house style
#           of the common-variant table (saige_signals.xlsx): booktabs rules,
#           three header levels with footnote letters, cohorts in their figure
#           colours, italic gene symbols, P as m x 10^-e, footnotes as short lines
#           below each table. Presentation only: every number is taken as is from
#           the assembled tables. Nine sheets, each after the sheets it refers to:
#             Rare-variant QC          the rare-variant QC table (layout of the
#                                      project's QC summary)
#             Variant sets             per cohort: the analysed rare variants by
#                                      impact, their mapping to genes, and the genes
#                                      and variants each variant set tests
#             Consequence terms        per cohort: the snpEff terms within each
#                                      impact (with totals), analysed and tested
#             Gene-based scans         per cohort x set: genes tested, thresholds,
#                                      significant genes and lambda_GC of both
#                                      statistics side by side (as observed)
#             Gene-based signals       gene x cohort: variant burden and SKAT-O /
#                                      CMC results of both variant sets
#             Gene-based variants      gene x variant, in the columns of the common-
#                                      variant table, plus per cohort whether the
#                                      gene test used the variant, and why not
#             Leave-one-variant-out    gene x variant (the order of Gene-based
#                                      variants): -log10 P with all variants and
#                                      without each, both sets x statistics x
#                                      cohorts, tiered sets marked (descriptive)
#             Gene-based Open Targets  every Open Targets record behind the summary
#                                      columns (as the common-variant table's sheet)
#             Multi-site carriers      the carriers of minor alleles at >= 2 of a
#                                      gene's tested variants, by sample ID
#           The workbook names no file, script, parameter or pipeline. Sample IDs and
#           project names appear only in sheet Multi-site carriers (asserted against
#           the real sample IDs on every other sheet).
# Inputs  : the assembled tables, Open Targets detail / meta, the rare-variant QC
#           json per cohort, the variant-set tables, the upstream manifest
# Outputs : gene_signals.xlsx, workbook.json
# ---------------------------------------------------------------------------
import argparse
import json
import re
import sys
import warnings
import zipfile
from math import floor, log10
from pathlib import Path

import numpy as np
import pandas as pd
import xlsxwriter

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gs_common import MISSING, csv_list, load_json, log, read_table, read_tsv  # noqa: E402

FORBIDDEN_TOKENS = ['_mainland', 'moderate_high', 'low_moderate_high', 'set_name', 'assoc_gene', 'tuning.rv',
                    'gene_scan', 'variant_stats', 'set_stats', 'Control_', 'Lung_cancer', 'JHRP', 'RADDAR', 'Old PH',
                    'bbj_mainland', 'pc1_avg', '.tsv', '.py', '.nf', '/LARGE', '/home', 'nextflow', 'rvtest ',
                    'not_in_callset', 'no_carrier', 'gap_filter', 'p.value']
# project names that sheet Multi-site carriers shows as recorded in the sample sheet
SAMPLE_SHEET_TOKENS = ('Control_', 'Lung_cancer', 'RADDAR', 'Old PH')
FONT = 'Arial'
INK, INK2, GREY, FAINT, RULE_GREY, NOTE = '#1A1A1A', '#595959', '#7F7F7F', '#A6A6A6', '#BFBFBF', '#404040'
IMPACT_STYLE = {'HIGH': ('#9C0006', True), 'MODERATE': ('#9C5700', True), 'LOW': (INK, False), 'MODIFIER': (GREY, False)}
HEADER_FILL = '#F2F2F2'
THIN, MEDIUM, HAIR = 1, 2, 7
NUM_AF = '[=0]0;[<0.0001]"<0.0001";0.0000'
METHOD_LABEL = {'skato': 'SKAT-O', 'cmc': 'CMC burden'}
METHOD_SHORT = {'skato': 'SKAT-O', 'cmc': 'CMC'}
SH_SIG, SH_VAR, SH_OT, SH_SETS, SH_TERMS, SH_SCANS, SH_QC, SH_MS, SH_LOVO = (
    'Gene-based signals', 'Gene-based variants', 'Gene-based Open Targets', 'Variant sets', 'Consequence terms',
    'Gene-based scans', 'Rare-variant QC', 'Multi-site carriers', 'Leave-one-variant-out')
# every sheet after the sheets it refers to: input QC -> what the tests were given -> the genome-wide scans ->
# the genes -> their variants -> which variants carry the signal -> external evidence on them -> the individuals
SHEETS = [SH_QC, SH_SETS, SH_TERMS, SH_SCANS, SH_SIG, SH_VAR, SH_LOVO, SH_OT, SH_MS]
REASON_TEXT = {'mac': 'No: MAC {mac}', 'no_carrier': 'No: no carrier', 'common': 'No: common (MAF {maf})',
               'gap_filter': 'No: depth-gap filter'}
# sheet Gene-based Open Targets (as the common-variant table's Open Targets sheet)
QTL_COL = 'cis / trans (to target TSS)'
POP_COL = 'Population (% of sample)'
OT_H = 2
GNOMAD_POP = {'AFR': 'African / African American', 'AMI': 'Amish', 'AMR': 'admixed American',
              'ASJ': 'Ashkenazi Jewish', 'CSA': 'Central / South Asian', 'EAS': 'East Asian', 'FIN': 'Finnish',
              'MID': 'Middle Eastern', 'NFE': 'non-Finnish European', 'SAS': 'South Asian'}
OT_SOURCES = ('ClinVar', 'GWAS', 'molQTL', 'E2G')
VARIANT_FILL = '#F2F2F2'
LOCUS_FILL = '#E4E4E4'


# ── text helpers ─────────────────────────────────────────────────────────────
def p_parts(p):
    """P for print: 3.84 × 10⁻⁸ below 0.001, else 2-3 decimals."""
    if p is None or not np.isfinite(p):
        return None
    if p < 1e-3:
        e = floor(log10(p))
        m = p / 10 ** e
        if round(m, 2) >= 10:
            m, e = m / 10, e + 1
        return [f'{m:.2f} × 10', ('sup', f'−{-e}')]
    return [f'{p:.2f}' if p >= 0.01 else f'{p:.3f}']


def sci_parts(x):
    e = floor(log10(x))
    m = x / 10 ** e
    head = '10' if abs(m - 1) < 1e-9 else f'{m:.2f} × 10'
    return [head, ('sup', f'−{-e}')]


def sig3(x):
    return f'{x:.0f}' if x >= 100 else f'{x:.1f}' if x >= 10 else f'{x:.2f}'


def so(term):
    return str(term).replace('&', ' & ').replace('_', ' ')


def text_len(parts):
    return sum(len(p if isinstance(p, str) else p[1]) for p in parts)


def blank(x):
    return x is None or (isinstance(x, float) and not np.isfinite(x)) or (isinstance(x, str) and x == '')


def s_(x):
    return '' if blank(x) else str(x)


def is_true(x):
    return str(x) in ('True', 'true', '1', '1.0')


def plain(v, sup=False):
    """Text of a cell read back with rich_text=True, without its superscript runs (sup=True: only those)."""
    from openpyxl.cell.rich_text import CellRichText
    if isinstance(v, CellRichText):
        return ''.join((b if isinstance(b, str) else b.text) for b in v
                       if (not isinstance(b, str) and b.font is not None and b.font.vertAlign == 'superscript') == sup)
    return '' if (v is None or sup) else str(v)


def wrap_parts(parts, width):
    """Split rich parts (str / (style, str)) into lines of <= width characters, at spaces."""
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


def platform_label(p):
    """Sequencer labels as the manufacturers write them; the platform groups are never merged."""
    p = re.sub(r'DNB[Ss]eq', 'DNBSEQ', str(p)).replace('HiSeqX', 'HiSeq X')
    return re.sub(r'(\d+)x$', r'\1×', p)


def source_label(s):
    return str(s).replace('Control_', '').replace('_', ' ')


def stratum_label(tag, short=False):
    parts = tag.split('_')
    if short:
        return '/'.join(p[0].upper() for p in parts)
    t = '/'.join(parts) + '-impact variants'
    return t[0].upper() + t[1:]


def dist_text(bp):
    return f'{bp:,.0f} bp' if bp < 1e3 else f'{bp / 1e3:,.0f} kb' if bp < 1e6 else f'{bp / 1e6:.1f} Mb'


def qtl_text(x):
    """cis / trans of a molQTL record with its evidence (as the common-variant table)."""
    if x.qtl_class in ('', 'unknown'):
        t = 'unknown (no gene location)'
    elif pd.isna(x.qtl_distance):
        t = f'{x.qtl_class} (gene on {x.target_chrom})'
    else:
        t = f'{x.qtl_class} ({dist_text(x.qtl_distance)})'
    if (x.qtl_class, x.ot_is_trans) in (('cis', 'true'), ('trans', 'false')):
        t += f' · Open Targets: {"trans" if x.ot_is_trans == "true" else "not trans"}'
    return t


def source_text(p_src, b_src):
    t = {'variant': 'this variant', 'lead': 'set lead'}
    p, b = t.get(str(p_src), '—'), t.get(str(b_src), '—')
    return p if p == b else f'P: {p}; β: {b}'


# ── the writer ───────────────────────────────────────────────────────────────
class Sheet:
    """An XlsxWriter worksheet with a format cache and a rich-text helper."""

    def __init__(self, wb, ws):
        self.wb, self.ws, self._fmt = wb, ws, {}

    def fmt(self, **kw):
        props = {'font_name': FONT, 'font_size': 9, 'font_color': INK, 'valign': 'vcenter'}
        props.update({k: v for k, v in kw.items() if v is not None})
        key = tuple(sorted(props.items()))
        if key not in self._fmt:
            self._fmt[key] = self.wb.add_format(props)
        return self._fmt[key]

    def font(self, style, b, i, color, sz):
        """A font-only format for one rich-text fragment; style: '', 'i', 'b', 'bi' (bold italic), 'sup', 'sub',
        '#RRGGBB' (bold colour)."""
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
    """A booktabs table: spanner / sub-spanner / leaf header rows, then body rows with rules.
    spec: list of dicts {blk, sub, leaf (str or parts), note, w, coh}; blk None = a column without spanner;
    span: {blk: (parts, colour, cohort or None)}."""

    def __init__(self, sh, spec, top, span, colors, tints, header_h=None):
        self.sh, self.spec, self.top, self.colors, self.tints = sh, spec, top, colors, tints
        self.ncol = len(spec)
        self.first_coh = {}
        for j, s in enumerate(spec):
            c = s.get('coh')
            if c and (j == 0 or spec[j - 1].get('coh') != c):
                self.first_coh[j] = c
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
                f = self._hcell(j, R0, bottom=MEDIUM, bottom_color=INK, valign='bottom')
                ws.merge_range(R0, j, R0 + 2, j, '', sh.fmt(**f))
                sh.rich(R0, j, self._leaf(j), f, b=True)
                j += 1
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
                sh.rich(R0 + 1, m, sub_parts + [('sup', note)], f, b=True)
                for q in range(m, n + 1):
                    leaf = spec[q]['leaf']
                    sh.rich(R0 + 2, q, leaf if isinstance(leaf, list) else [('b', leaf)], self._hcell(q, R0 + 2),
                            b=True)
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
    sh.rich(0, 0, [('b', text)], {'align': 'left'}, sz=11)
    sh.rich(1, 0, sub if isinstance(sub, list) else [sub], {'align': 'left'}, i=True, color=INK2)
    sh.ws.set_row(0, 18)


def new_sheet(wb, name, zoom=90):
    sh = Sheet(wb, wb.add_worksheet(name))
    sh.ws.hide_gridlines(2)
    sh.ws.set_zoom(zoom)
    sh.ws.set_default_row(12.75)
    return sh


def print_setup(ws, last_row, ncol, footer, repeat=None, landscape=True, paper=8, scale=None):
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
    """The meaning of every placeholder, one line per symbol and column: [(symbol, where, meaning), ...]."""
    out = [('', [('b', 'Missing values')])]
    for sym, where, text in items:
        out.append(('', ['   ', ('b', sym)] + ([f', in {where}: '] if where else [': ']) + text))
    return out


def abbrev(d):
    return ('', ['Abbreviations: ' + '; '.join(f'{k}, {v}' for k, v in sorted(d.items(), key=lambda kv:
                                                                               kv[0].lower())) + '.'])


ABBR = {'BH': 'Benjamini–Hochberg', 'CI': 'confidence interval', 'CMC': 'combined multivariate and collapsing '
        '(burden) test', 'MAC': 'minor allele count', 'MAF': 'minor allele frequency', 'OR': 'odds ratio',
        'SKAT-O': 'optimal sequence kernel association test', 'PC': 'principal component',
        'AGP3K': 'control cohort', 'CTEPH': 'chronic thromboembolic pulmonary hypertension',
        'EA/OA': 'effect/other allele', 'EAF': 'effect-allele frequency', 'EAS': 'East Asian',
        'HWE': 'Hardy–Weinberg equilibrium', 'JPT': 'Japanese in Tokyo', 'L2G': 'locus-to-gene',
        'NFE': 'non-Finnish European', 'QTL': 'quantitative trait locus', 'E2G': 'enhancer-to-gene',
        'M/H': 'moderate/high-impact variant set', 'L/M/H': 'low/moderate/high-impact variant set'}


# ── the workbook ─────────────────────────────────────────────────────────────
def build(path, D, a):
    cohorts, strata = csv_list(a.cohorts), csv_list(a.strata)
    codes = dict(zip(cohorts, csv_list(a.cohort_codes)))
    labels = dict(zip(cohorts, csv_list(a.cohort_labels)))
    colors = {c: '#' + x.upper() for c, x in zip(cohorts, csv_list(a.cohort_colors))}
    tints = {c: '#' + x.upper() for c, x in zip(cohorts, csv_list(a.cohort_tints))}
    code_color = {codes[c]: colors[c] for c in cohorts}
    G, VL, SC, genes, up = D['G'], D['VL'], D['SC'], D['genes'], D['up']
    otd, otm, QC, TERMS, FLOW, MSD = D['otd'], D['otm'], D['qc'], D['terms'], D['flow'], D['ms']
    ot_on = bool(otm.get('queried'))
    methods = list(METHOD_LABEL)
    widest = max(strata, key=lambda s: int(VL[f'in_set_{s}'].sum()))
    narrow_s = [s for s in strata if s != widest][0]
    assert not (VL[f'in_set_{narrow_s}'] & ~VL[f'in_set_{widest}']).any(), 'the variant sets must be nested'
    strata = [narrow_s, widest]
    minac = {c: int(QC[c]['minac']['recommended_minac']) for c in cohorts}
    alpha = float(up['alpha'])
    n_pcs = len([c for c in up['covariates'].split(',') if c.lower().startswith('pc')])
    min_var = int(up['min_num_var'])
    g_order = genes.sort_values('gene_order')
    by_gsc = G.set_index(['set_name', 'stratum', 'cohort'])
    n_case = {c: int(G[G.cohort == c].n_case.iloc[0]) for c in cohorts}
    n_ctrl = {c: int(G[G.cohort == c].n_ctrl.iloc[0]) for c in cohorts}
    cls = {r.set_name: {'B': 'Single cohort', 'S': 'Sensitivity only'}.get(r['class'], r['class'])
           for _, r in genes.iterrows()}
    sens_min = a.sens_min_num_var
    sens_sub = f'Sensitivity: ≥ {sens_min} variants'
    cls_n = pd.Series(list(cls.values())).value_counts().to_dict()
    code_parts = []
    for k, c in enumerate(cohorts):
        code_parts += ([' · '] if k else []) + [(colors[c], codes[c])]
    cohort_text = '; '.join(f'{codes[c]}, {labels[c].lower()}' for c in cohorts)
    nested = ' ⊂ '.join(labels[c].lower() for c in cohorts)
    mac_txt = ' / '.join(sorted({str(v) for v in minac.values()}))

    def cohort_parts(lead):
        """The three cohorts: codes, sizes, nesting, and why agreement across them is not replication."""
        return [('b', lead), ': ' + '; '.join(f'{codes[c]}, {labels[c].lower()} ({n_case[c]:,} cases / {n_ctrl[c]:,} '
                                               'controls)' for c in cohorts)
                + f'. Nested ({" ⊂ ".join(codes[c] for c in cohorts)}) ancestry-filtered mainland Japanese sample sets, '
                'unrelated within each; they share samples, so agreement across them is not independent replication.']
    wb = xlsxwriter.Workbook(path, {'nan_inf_to_errors': False})
    label = f'{a.table_label}. ' if a.table_label.strip() else ''
    wb.set_properties({'title': f'{label}CTEPH gene-based rare-variant association signals',
                       'author': 'CTEPH WGS project'})
    span_s = {s: ([('b', stratum_label(s))], INK, None) for s in strata}
    stats = {}
    # every sheet created up front, in the workbook's order; each section fills its own
    sheets = {}
    for name in SHEETS:
        if name == SH_OT and not ot_on:
            continue
        sheets[name] = wb.add_worksheet(name) if name == SH_QC else new_sheet(wb, name)

    # ════════════════════════════════════════════════════════════════════════
    # 1. Gene-based signals: one row per gene x variant set x cohort, the two sets stacked under each gene
    # ════════════════════════════════════════════════════════════════════════
    sh = sheets[SH_SIG]
    ws = sh.ws
    title(sh, f'{label}CTEPH gene-based rare-variant association signals across three nested ancestry-filtered '
              'cohorts', f'{len(genes)} genes (' + ', '.join(f'{cls_n[k]} {k}' for k in ('Tier 1', 'Tier 2', 'Tier 3',
                                                                                         'Single cohort',
                                                                                         'Sensitivity only')
                                                         if cls_n.get(k)) + ') · per gene, each variant set with its '
          'three cohorts · shaded and bold, BH-significant (q < 0.05); †, below the Bonferroni threshold')
    CALL_FILL = '#FCEFD2'
    spec = [dict(blk=None, sub=None, leaf='Gene', note='a', w=10), dict(blk=None, sub=None, leaf='Chr', w=4.5),
            dict(blk=None, sub=None, leaf='Class', note='a', w=10.5),
            dict(blk=None, sub=None, leaf='Variant set', note='b', w=8),
            dict(blk=None, sub=None, leaf='Cohort', note='c', w=18),
            dict(blk='bur', sub=None, leaf='Variants', note='d', w=7.5),
            dict(blk='bur', sub='Carriers, n (%)', leaf='Case', note='d', w=10.5),
            dict(blk='bur', sub='Carriers, n (%)', leaf='Control', note='d', w=10.5),
            dict(blk='bur', sub=None, leaf='Cumulative MAC case / control', note='d', w=10),
            dict(blk='bur', sub=None, leaf='Cumulative MAF case / control', note='d', w=13.5),
            dict(blk='bur', sub=None, leaf='Multi-site carriers case / control', note='e', w=9.5),
            dict(blk='test', sub='SKAT-O', leaf='ρ', note='f', w=5),
            dict(blk='test', sub='SKAT-O', leaf='P', note='f', w=10.5),
            dict(blk='test', sub='SKAT-O', leaf='q', note='f', w=8.5),
            dict(blk='test', sub='CMC burden', leaf='OR (95% CI)', note='g', w=16),
            dict(blk='test', sub='CMC burden', leaf='P', note='g', w=10.5),
            dict(blk='test', sub='CMC burden', leaf='q', note='g', w=8.5),
            dict(blk='test', sub=None, leaf='Tier', note='a', w=5.5),
            dict(blk='sens', sub=None, leaf='SKAT-O q', note='h', w=10),
            dict(blk='sens', sub=None, leaf='CMC q', note='h', w=10),
            dict(blk='sens', sub=None, leaf='Tier', note='h', w=5.5)]
    span = {'bur': ([('b', 'Variant burden')], INK, None),
            'test': ([('b', f'Gene-based tests (primary: genes with ≥ {min_var} variants)')], INK, None),
            'sens': ([('b', f'Sensitivity (≥ {sens_min} variants)')], INK, None)}
    T = Table(sh, spec, 2, span, colors, tints, header_h=39)
    nc = len(cohorts)
    per_gene = nc * len(strata)
    r = T.row0
    last = r + len(g_order) * per_gene - 1
    n_sep = 0
    for _, gr in g_order.iterrows():
        sn = gr.set_name
        rg = r
        g_last = rg + per_gene - 1 == last
        for si, s in enumerate(strata):
            rs = r
            s_last = rs + nc - 1 == last
            for ci, c in enumerate(cohorts):
                gtop = si == 0 and ci == 0 and rg > T.row0
                stop = si > 0 and ci == 0
                is_last = r == last
                x = by_gsc.loc[(sn, s, c)]

                def put(j, value, align='center', num=None, **kw):
                    if blank(value):
                        ws.write_string(r, j, MISSING, sh.fmt(**T.rules(j, {'align': 'center', 'font_color': FAINT},
                                                                      group_top=gtop, soft_top=stop, last=is_last)))
                        return
                    d = T.rules(j, {'align': align, 'num_format': num}, group_top=gtop, soft_top=stop, last=is_last)
                    d.update(kw)
                    if isinstance(value, str):
                        ws.write_string(r, j, value, sh.fmt(**d))
                    else:
                        ws.write_number(r, j, float(value), sh.fmt(**d))

                def rich(j, parts, align='left', fill=None, **kw):
                    d = T.rules(j, {'align': align}, group_top=gtop, soft_top=stop, last=is_last)
                    if fill:
                        d.update({'bg_color': fill, 'pattern': 1})
                    sh.rich(r, j, parts, d, **kw)

                def pcell(j, p, called, bonf):
                    rich(j, p_parts(float(p)) + ([('sup', '†')] if bonf else []), align='right',
                         fill=CALL_FILL if called else None, b=bool(called))
                rich(4, [(colors[c], labels[c])])
                if not x.tested:
                    for j in range(5, 17):
                        put(j, None)
                    put(18, None)
                    put(19, None)
                else:
                    put(5, int(x.n_var), num='0')
                    put(6, f'{int(x.carriers_case)} ({100 * x.carriers_case / x.n_case:.1f}%)', align='right')
                    put(7, f'{int(x.carriers_ctrl)} ({100 * x.carriers_ctrl / x.n_ctrl:.1f}%)', align='right')
                    put(8, f'{int(x.mac_case)} / {int(x.mac_ctrl)}')
                    put(9, f'{x.maf_case:.4f} / {x.maf_ctrl:.4f}')
                    mm = int(x.multi_case) + int(x.multi_ctrl)
                    put(10, f'{int(x.multi_case)} / {int(x.multi_ctrl)}', font_color=None if mm else GREY)
                    put(11, float(x.skato_rho), num='0.0#')
                    pcell(12, x.skato_p, x.skato_bh, x.skato_bonferroni)
                    rich(13, p_parts(float(x.skato_q)), align='right')
                    if x.carriers_case == 0 or x.carriers_ctrl == 0:
                        n_sep += 1
                        put(14, None)
                    else:
                        put(14, f'{sig3(x.cmc_or)} ({sig3(x.cmc_l95)}–{sig3(x.cmc_u95)})', align='right')
                    pcell(15, x.cmc_p, x.cmc_bh, x.cmc_bonferroni)
                    rich(16, p_parts(float(x.cmc_q)), align='right')
                    for j, m in ((18, 'skato'), (19, 'cmc')):
                        if bool(x.in_sens):
                            rich(j, p_parts(float(x[f'{m}_q_sens'])) + ([('sup', '†')] if bool(x[f'{m}_bonferroni_sens'])
                                                                         else []),
                                 align='right', fill=CALL_FILL if bool(x[f'{m}_bh_sens']) else None,
                                 b=bool(x[f'{m}_bh_sens']))
                        else:
                            put(j, None)
                r += 1
            # merged per set: the variant set, its tier and its sensitivity tier
            x0 = by_gsc.loc[(sn, s, cohorts[0])]
            for j, val, kw in ((3, stratum_label(s, short=True), {'bold': True}),
                               (17, (str(int(x0.tier)) if int(x0.tier) else '–'),
                                {'bold': True if int(x0.tier) else None, 'font_color': INK if int(x0.tier) else GREY}),
                               (20, (MISSING if pd.isna(x0.tier_sens) else (str(int(x0.tier_sens)) if x0.tier_sens
                                                                            else '–')),
                                {'bold': True if (not pd.isna(x0.tier_sens) and x0.tier_sens) else None,
                                 'font_color': FAINT if pd.isna(x0.tier_sens) else (INK if x0.tier_sens else GREY)})):
                f = T.rules(j, dict({'align': 'center', 'valign': 'vcenter'}, **kw), group_top=si == 0 and rg > T.row0,
                            soft_top=si > 0, last=s_last)
                ws.merge_range(rs, j, rs + nc - 1, j, val, sh.fmt(**f))
        # merged per gene: gene, chromosome, class
        for j, parts, kw in ((0, [('i', gr.gene_symbol)], {'b': True}), (1, [str(gr.chrom)], {}),
                             (2, [cls[sn]], {'b': cls[sn].startswith('Tier')})):
            f = T.rules(j, {'align': 'left' if j != 1 else 'center', 'valign': 'top'}, group_top=rg > T.row0,
                        last=g_last)
            ws.merge_range(rg, j, rg + per_gene - 1, j, '', sh.fmt(**f))
            sh.rich(rg, j, parts, f, **kw)
    assert r - 1 == last
    ws.freeze_panes(T.row0, 5)
    thr = [float(x) for x in SC.bonferroni_threshold]
    thr_s = [float(x) for x in SC.bonferroni_threshold_sens]
    n_untested = int((~G.tested).sum())
    notes = [
        ('a', [('b', 'Gene'), ', the snpEff gene symbol the variants were mapped to. ', ('b', 'Class'), ', the better '
               'tier of the two variant sets; Single cohort, not tiered but Bonferroni-significant in one cohort; '
               'Sensitivity only, selected by the sensitivity analysis alone (note h). ', ('b', 'Tier'),
               f', per variant set, from the BH calls (q < {alpha:g}): 1, called in all three cohorts and by both '
               'statistics in at least one; 2, called in all three cohorts; 3, called in two cohorts. Genes shown: a '
               'tier in either set, or Bonferroni-significant in any cohort, set and statistic.']),
        ('b', [('b', 'Variant set'), ': M/H, moderate/high-impact variants; L/M/H, low/moderate/high-impact variants '
               '(sheet ' + SH_SETS + ').']),
        ('c', cohort_parts('Cohort')),
        ('d', [('b', 'Variant burden'), ' of the set as the cohort tested it. Variants, its qualifying variants (MAC ≥ '
               + mac_txt + ' in the cohort); Carriers, samples with ≥ 1 minor allele at any of them, n (% of cases or '
               'of controls); Cumulative MAC, minor alleles summed over the variants; Cumulative MAF, cumulative MAC / '
               '(2 × group size). The minor allele is defined on the pooled samples.']),
        ('e', [('b', 'Multi-site carriers'), ', carriers with minor alleles at two or more of the gene\'s variants, '
               'counted once among the carriers; listed by sample in sheet ' + SH_MS + '.']),
        ('f', [('b', 'SKAT-O'), ', Beta(1, 25) weights on MAF. ρ, the fitted mix of kernel (0) and burden (1) test; '
               'P, score-test P; q, BH q-value over the genes tested in that cohort and set.']),
        ('g', [('b', 'CMC burden'), ', carrier status (≥ 1 minor allele at any of the gene\'s variants) as the '
               'predictor. OR (95% CI), from the Wald fit of the same model; P, score test; q, as in note f.']),
        ('h', [('b', 'Sensitivity analysis'), f', the same P judged within the genes with {sens_min} or more '
               'qualifying variants: BH q over that family; Tier as in note a, from these calls.']),
        ('', [('b', 'Model'), f': logistic regression adjusted for sex and {n_pcs} principal components; missing '
              'genotypes imputed to the mean. No correction across statistics, variant sets or cohorts.']),
        ('', [('b', 'Marks'), f': shaded and bold, BH-significant (q < {alpha:g}); †, below the Bonferroni '
              f'threshold, {alpha:g} / genes tested (', *sci_parts(min(thr)), ' to ', *sci_parts(max(thr)),
              '; sensitivity, ', *sci_parts(min(thr_s)), ' to ', *sci_parts(max(thr_s)), '; sheet ' + SH_SCANS + ').']),
        *missing_block([
            ('—', 'OR (95% CI)', ['not estimable: no carrier in one group.']),
            ('—', 'the sensitivity columns', [f'the cohort tested the gene on fewer than {sens_min} variants, so the '
                                              'gene is outside that family; in its Tier, no cohort did.'])]
            + ([('—', 'every column of a cohort', [f'the cohort did not test the gene (fewer than {min_var} '
                                                   'qualifying variants).'])] if n_untested else [])
            + [('–', 'Tier', ['not tiered (called in fewer than two cohorts).'])]),
        abbrev({k: ABBR[k] for k in ('BH', 'CI', 'CMC', 'MAC', 'MAF', 'OR', 'SKAT-O')}),
    ]
    rr, lens = write_notes(sh, last + 2, notes, 150)
    stats['signals'] = {'rows': last - T.row0 + 1, 'or_not_estimable': n_sep, 'max_note': max(lens)}
    print_setup(ws, rr - 1, T.ncol, f'{label}Gene-based signals', repeat=(2, 4))

    # ════════════════════════════════════════════════════════════════════════
    # Multi-site carriers: one row per gene and sample, with the sample's sequencing and clinical information
    # ════════════════════════════════════════════════════════════════════════
    sh = sheets[SH_MS]
    ws = sh.ws
    n_samples = MSD.IID.nunique()
    title(sh, f'{label}Multi-site carriers of the gene-based signals',
          f'{n_samples} samples carrying minor alleles at two or more of a gene\'s tested variants · one row per gene '
          'and sample, with the sample\'s sequencing and clinical information · cases first')
    spec = [dict(blk=None, sub=None, leaf='Sample ID', note='a', w=13),
            dict(blk=None, sub=None, leaf='Group', w=8),
            dict(blk=None, sub=None, leaf='Variants carried (genotype)', note='b', w=58),
            dict(blk='ms', sub=None, leaf=stratum_label(narrow_s, short=True), note='c', w=8),
            dict(blk='ms', sub=None, leaf=stratum_label(widest, short=True), note='c', w=8),
            dict(blk='seq', sub=None, leaf='Project', w=16), dict(blk='seq', sub=None, leaf='Platform', note='d', w=17),
            dict(blk='seq', sub=None, leaf='Target depth', w=7.5),
            dict(blk='seq', sub=None, leaf='Observed depth', note='d', w=8.5),
            dict(blk='clin', sub=None, leaf='Sex', w=5),
            dict(blk='clin', sub=None, leaf='Age at DNA collection', note='e', w=10),
            dict(blk='clin', sub=None, leaf='Age at diagnosis', note='f', w=9),
            dict(blk='clin', sub=None, leaf='DVT', note='f', w=6),
            dict(blk='clin', sub=None, leaf='Max mPAP (mmHg)', note='f', w=10)]
    span = {'ms': ([('b', 'Multi-site carrier in')], INK, None), 'seq': ([('b', 'Sequencing')], INK, None),
            'clin': ([('b', 'Clinical')], INK, None)}
    T = Table(sh, spec, 2, span, colors, tints, header_h=26)
    hg = VL.drop_duplicates(['set_name', 'variant_id']).set_index(['set_name', 'variant_id'])
    genes_ms = [g for g in g_order.set_name if g in set(MSD.set_name)]
    r = T.row0
    last = r + len(MSD) + len(genes_ms) - 1
    for sn in genes_ms:
        xg = MSD[MSD.set_name == sn]
        n_ca, n_co = int((xg.group == 'case').sum()), int((xg.group == 'ctrl').sum())
        hfmt = dict(T.rules(0, {}, group_top=r > T.row0), bg_color=HEADER_FILL, pattern=1, align='left')
        sh.rich(r, 0, [('i', xg.gene_symbol.iloc[0]), f'   {len(xg)} sample{"s" if len(xg) != 1 else ""}: {n_ca} '
                       f'case{"s" if n_ca != 1 else ""}, {n_co} control{"s" if n_co != 1 else ""}'], hfmt, b=True)
        for j in range(1, T.ncol):
            ws.write_blank(r, j, None, sh.fmt(**dict(T.rules(j, {}, group_top=r > T.row0), bg_color=HEADER_FILL,
                                                    pattern=1)))
        r += 1
        for _, x in xg.iterrows():
            is_last = r == last
            case = x.group == 'case'

            def put(j, value, align='left', num=None, **kw):
                if blank(value):
                    ws.write_string(r, j, MISSING, sh.fmt(**T.rules(j, {'align': 'center', 'font_color': FAINT},
                                                                  last=is_last)))
                    return
                d = T.rules(j, {'align': align, 'num_format': num}, last=is_last)
                d.update(kw)
                if isinstance(value, str):
                    ws.write_string(r, j, value, sh.fmt(**d))
                else:
                    ws.write_number(r, j, float(value), sh.fmt(**d))

            def rich(j, parts, align='left', **kw):
                sh.rich(r, j, parts, T.rules(j, {'align': align}, last=is_last), **kw)
            put(0, x.IID, bold=True if case else None)
            put(1, 'Case' if case else 'Control', bold=True if case else None)
            parts = []
            for k, (vid, dose, imp) in enumerate(zip(x.variant_ids.split(';'), x.dosages.split(';'),
                                                     x.impacts.split(';'))):
                h = hg.loc[(x.set_name, vid)]
                name = s_(h.hgvs_p) or s_(h.hgvs_c) or vid
                c_, b_ = IMPACT_STYLE.get(imp, (INK, False))
                parts += ([' · '] if k else []) + [(c_ if b_ else '', name), f' ({"hom" if dose == "2" else "het"})']
            rich(2, parts)
            for q, st in enumerate(strata):
                cp = []
                for c in cohorts:
                    if codes[c] in s_(x[f'cohorts_{st}']):
                        cp += ([' '] if cp else []) + [(colors[c], codes[c])]
                if cp:
                    rich(3 + q, cp, align='center')
                else:
                    put(3 + q, '–', align='center', font_color=GREY)
            put(5, s_(x.project))
            put(6, s_(x.platform))
            put(7, s_(x.target_depth), align='center')
            put(8, float(x.observed_depth), num='0.0', align='right')
            put(9, s_(x.sex), align='center')
            put(10, None if blank(x.age_at_dna) else float(x.age_at_dna), num='0.#', align='right')
            if case:
                put(11, None if blank(x.age_at_diagnosis) else float(x.age_at_diagnosis), num='0', align='right')
                put(12, s_(x.dvt) or None, align='center')
                put(13, None if blank(x.max_mpap) else float(x.max_mpap), num='0', align='right')
            else:
                for j in (11, 12, 13):
                    put(j, 'n/a', align='center', font_color=GREY)
            r += 1
    assert r - 1 == last
    ws.freeze_panes(T.row0, 1)
    none_ms = [g for g in g_order.gene_symbol if g not in set(MSD.gene_symbol)]
    notes = [
        ('a', [('b', 'Sample ID'), ': samples with minor alleles at two or more of a gene\'s tested variants in at '
               'least one cohort and variant set, i.e. the multi-site carriers counted in sheet ' + SH_SIG + '. A grey '
               'row per gene, cases first; genes without any: ' + (', '.join(none_ms) if none_ms else 'none') + '.']),
        ('b', [('b', 'Variants carried'), ': the gene\'s tested variants at which the sample carries a minor allele, '
               'in genomic order, as HGVS on the transcript of sheet ' + SH_VAR + '; colour, impact (High dark red, '
               'Moderate dark orange, Low black); het, one minor allele; hom, two.']),
        ('c', [('b', 'Multi-site carrier in'), ': the cohorts whose test of the gene in that variant set counts the '
               'sample as a multi-site carrier.']),
        ('d', [('b', 'Project, Platform'), ': as recorded in the sample sheet (platform names verbatim). ',
               ('b', 'Observed depth'), ': mean sequencing depth of the sample.']),
        ('e', [('b', 'Age at DNA collection'), ' (years): cases, from the project\'s sample table; controls, from the '
               'AGP3K age and sex release.']),
        ('f', [('b', 'Clinical'), ' (cases only, from the project\'s sample table): age at diagnosis; DVT, deep-vein '
               'thrombosis; Max mPAP, highest mean pulmonary arterial pressure (mmHg).']),
        ('', cohort_parts('Cohorts')),
        *missing_block([('—', '', ['not recorded.']), ('n/a', '', ['not applicable to controls.']),
                        ('–', 'Multi-site carrier in', ['no cohort (e.g. the second site is a low-impact variant, '
                                                         'tested only in L/M/H).'])]),
        abbrev({'DVT': 'deep-vein thrombosis', 'mPAP': 'mean pulmonary arterial pressure', **ABBR_SET}),
    ]
    rr, lens = write_notes(sh, last + 2, notes, 140)
    stats['multisite'] = {'rows': len(MSD), 'samples': int(n_samples), 'genes_without': none_ms}
    print_setup(ws, rr - 1, T.ncol, f'{label}Multi-site carriers', repeat=(2, 4))

    # ════════════════════════════════════════════════════════════════════════
    # Leave-one-variant-out: every gene in both variant sets side by side (tier marked); -log10 P with all variants
    # (gene row) and with each variant left out (variant rows)
    # ════════════════════════════════════════════════════════════════════════
    sh = sheets[SH_LOVO]
    ws = sh.ws
    LV = D['lovo']
    # shading by the fall in -log10 P when the variant is left out (delta); ColorBrewer Oranges, Blues for a rise
    # a diverging key (ColorBrewer RdBu): red, -log10 P falls without the variant; blue, it rises; three steps each
    D_FALL = ((3.0, '#D6604D', '#FFFFFF'), (2.0, '#F4A582', INK), (1.0, '#FDDBC7', INK))
    D_RISE = ((3.0, '#4393C3', '#FFFFFF'), (2.0, '#92C5DE', INK), (1.0, '#D1E5F0', INK))
    GENE_FILL = '#EDEDED'
    TIER_FILL, SENS_FILL = '#A9D08E', '#E2EFDA'
    lv = LV.set_index(['set_name', 'stratum', 'cohort', 'method', 'k'])
    lv_v = LV[LV.k > 0].set_index(['set_name', 'stratum', 'cohort', 'method', 'dropped'])
    hg = VL.drop_duplicates(['set_name', 'variant_id']).set_index(['set_name', 'variant_id'])
    by_vc = VL.set_index(['set_name', 'variant_id', 'cohort'])
    meths = ('skato', 'cmc')
    title(sh, f'{label}Leave-one-variant-out analysis of the gene-based signals',
          f'{len(g_order)} genes × 2 variant sets × 2 statistics × 3 cohorts · −log10 P with all variants (grey row) '
          'and with each variant left out · red, a fall; blue, a rise · Tier marks the tiered sets · descriptive')

    # Consequence and Impact exactly as sheet Gene-based variants writes them
    spec = [dict(blk=None, sub=None, leaf='Gene / HGVS', note='a', w=24),
            dict(blk=None, sub=None, leaf='Variant', note='b', w=19),
            dict(blk=None, sub=None, leaf='rsID', note='b', w=12.5),
            dict(blk=None, sub=None, leaf='Consequence', note='c', w=35),
            dict(blk=None, sub=None, leaf='Impact', note='c', w=9)]
    J_VID, J_RS, J_CQ, J_IMP = 1, 2, 3, 4
    j_rs = J_IMP                                             # the last identity column
    spec += [dict(blk='car', sub=None, leaf=[(colors[c], codes[c])], w=7.5, vline=THIN if i == 0 else None)
             for i, c in enumerate(cohorts)]
    j_tier, j_lp = {}, {}
    for s in strata:
        j_tier[s] = len(spec)
        spec.append(dict(blk=s, sub=None, leaf='Tier', note='e', w=5, vline=THIN))
        for mi, m in enumerate(meths):
            for i, c in enumerate(cohorts):
                j_lp[(s, m, c)] = len(spec)
                spec.append(dict(blk=s, sub=f'{METHOD_SHORT[m]}, −log10 P', leaf=[(colors[c], codes[c])], w=6,
                                 vline=HAIR if i == 0 else None))
    span = {'car': ([('b', 'Carriers, case / control'), ('sup', 'd')], INK, None)}
    for s in strata:
        span[s] = ([('b', f'{stratum_label(s).replace(" variants", " set")} ({stratum_label(s, short=True)})'),
                    ('sup', 'f')], INK, None)
    T = Table(sh, spec, 2, span, colors, tints, header_h=None)
    j_car = {c: j_rs + 1 + i for i, c in enumerate(cohorts)}

    r = T.row0
    head_rows, n_two, n_car_missing, n_bins, n_rows_v, lovo_order, lovo_annot = [], 0, {}, {}, 0, [], {}
    n_tier = {'tier': 0, 'sens': 0, 'none': 0}
    for _, gr in g_order.iterrows():
        sn, gene = gr.set_name, gr.gene_symbol
        g3 = {s: {c: by_gsc.loc[(sn, s, c)] for c in cohorts} for s in strata}
        used_any = VL[VL.set_name == sn].groupby('variant_id')[[f'in_set_{s}' for s in strata]].any()
        vids = list(used_any[used_any.any(axis=1)].index)
        # the variants in genomic order, as in sheet Gene-based variants (asserted below)
        vids = sorted(vids, key=lambda v: (int(hg.loc[(sn, v)].pos), v))
        lovo_order += [(sn, v) for v in vids]
        # the gene row: per set, its tier and the published tests (all variants)
        lab = [('bi', gene), f'   {cls[sn]} · variants ']
        two_any = False
        for q, s in enumerate(strata):
            nv = [int(g3[s][c].n_var) if g3[s][c].tested else None for c in cohorts]
            two = any(n == 2 for n in nv if n is not None)
            two_any |= two
            lab += ([' · '] if q else []) + [f'{stratum_label(s, short=True)} '
                                             + (str(nv[0]) if len(set(nv)) == 1 else
                                                '/'.join(MISSING if n is None else str(n) for n in nv))] \
                + ([' ', ('b', '‡')] if two else [])
        n_two += two_any
        gfmt = lambda j, **kw: {**T.rules(j, {'align': 'center'}, group_top=r > T.row0), 'bg_color': GENE_FILL,
                                'pattern': 1, **kw}
        sh.rich(r, 0, lab, gfmt(0, align='left'))
        for j in (*range(1, j_rs + 1), *j_car.values()):
            ws.write_blank(r, j, None, sh.fmt(**gfmt(j)))
        for s in strata:
            x0 = g3[s][cohorts[0]]
            t, ts = int(x0.tier), x0.tier_sens
            tiered = bool(t) or (not pd.isna(ts) and bool(ts))
            if t:
                ws.write_string(r, j_tier[s], str(t), sh.fmt(**gfmt(j_tier[s], bold=True, bg_color=TIER_FILL)))
                n_tier['tier'] += 1
            elif not pd.isna(ts) and ts:
                ws.write_string(r, j_tier[s], f'S{int(ts)}', sh.fmt(**gfmt(j_tier[s], bold=True, bg_color=SENS_FILL)))
                n_tier['sens'] += 1
            else:
                ws.write_string(r, j_tier[s], '–', sh.fmt(**gfmt(j_tier[s], font_color=GREY)))
                n_tier['none'] += 1
            for c in cohorts:
                x = g3[s][c]
                for m in meths:
                    j = j_lp[(s, m, c)]
                    if not x.tested:
                        ws.write_string(r, j, MISSING, sh.fmt(**gfmt(j, font_color=FAINT)))
                        continue
                    p0 = float(lv.loc[(sn, s, c, m, 0)].p)
                    assert p0 == float(x[f'{m}_p']), (sn, s, c, m)
                    ws.write_number(r, j, -np.log10(p0), sh.fmt(**gfmt(j, bold=True if tiered else None,
                                                                        num_format='0.0', align='right')))
        head_rows.append((r + 1, gene))
        r += 1
        for vid in vids:
            h = hg.loc[(sn, vid)]

            def cell(j, **kw):
                d_ = T.rules(j, {'align': 'center'})
                d_.update({k_: v_ for k_, v_ in kw.items() if v_ is not None})
                return d_

            def put(j, value, **kw):
                if blank(value):
                    ws.write_string(r, j, MISSING, sh.fmt(**cell(j, font_color=FAINT)))
                elif isinstance(value, str):
                    ws.write_string(r, j, value, sh.fmt(**cell(j, **kw)))
                else:
                    ws.write_number(r, j, float(value), sh.fmt(**cell(j, **kw)))
            put(0, s_(h.hgvs_p) or s_(h.hgvs_c) or vid, align='left', indent=2)
            lovo_annot[(sn, vid)] = (so(h.effect), str(h.impact).title())
            put(J_VID, vid, align='left')
            put(J_RS, s_(h.rsID) or None, align='left')
            put(J_CQ, so(h.effect), align='left')
            c_, b_ = IMPACT_STYLE.get(h.impact, (INK, False))
            put(J_IMP, str(h.impact).title(), align='left', font_color=c_, bold=b_ or None)
            used_c = {c: any(bool(by_vc.loc[(sn, vid, c)][f'in_set_{s}']) for s in strata) for c in cohorts}
            for c in cohorts:
                v = by_vc.loc[(sn, vid, c)]
                if v.callset_reason == 'no_carrier':              # absent from the call set: no carrier in the cohort
                    put(j_car[c], '0 / 0', font_color=GREY)
                elif blank(v.case_carriers) or blank(v.ctrl_carriers):
                    n_car_missing[v.callset_reason] = n_car_missing.get(v.callset_reason, 0) + 1
                    put(j_car[c], None)
                else:
                    put(j_car[c], f'{int(v.case_carriers)} / {int(v.ctrl_carriers)}',
                        font_color=None if used_c[c] else GREY)
            for s in strata:
                in_set = bool(used_any.loc[vid, f'in_set_{s}'])
                ws.write_blank(r, j_tier[s], None, sh.fmt(**cell(j_tier[s])))
                for c in cohorts:
                    used = bool(by_vc.loc[(sn, vid, c)][f'in_set_{s}'])
                    for m in meths:
                        j = j_lp[(s, m, c)]
                        if not in_set:                            # outside this set (a low-impact variant in M/H)
                            ws.write_blank(r, j, None, sh.fmt(**cell(j)))
                            continue
                        if not used:
                            put(j, None)
                            continue
                        x = lv_v.loc[(sn, s, c, m, vid)]
                        drop = float(x.delta_log10)              # = -log10 P(all) - (-log10 P(without))
                        fall = next(((f, fc) for thr, f, fc in D_FALL if drop >= thr), None)
                        rise = next(((f, fc) for thr, f, fc in D_RISE if -drop >= thr), None)
                        if fall:
                            kw = {'bg_color': fall[0], 'pattern': 1, 'font_color': fall[1], 'bold': True}
                        elif rise:
                            kw = {'bg_color': rise[0], 'pattern': 1, 'font_color': rise[1]}
                        else:
                            kw = {}
                        n_bins[kw.get('bg_color', '')] = n_bins.get(kw.get('bg_color', ''), 0) + 1
                        put(j, -np.log10(float(x.p)), num_format='0.0', align='right', **kw)
            r += 1
            n_rows_v += 1
    last = r - 1
    for j in range(T.ncol):                                  # the closing rule
        ws.write_blank(r, j, None, sh.fmt(**{'top': MEDIUM, 'top_color': INK}))
    assert set(n_car_missing) <= {'common', 'gap_filter'}, n_car_missing
    assert sum(n_bins.values()) == int((LV.k > 0).sum()), 'every re-test in a cell'
    # the colour key, top right of the title row (above the frozen panes, so always in view):
    # strongest fall -> no change -> strongest rise
    kr = 0
    key = ([('falls ≥ 3', *D_FALL[0][1:], True), ('falls 2–3', *D_FALL[1][1:], True),
            ('falls 1–2', *D_FALL[2][1:], True), ('change < 1', None, INK, False),
            ('rises 1–2', *D_RISE[2][1:], False), ('rises 2–3', *D_RISE[1][1:], False),
            ('rises ≥ 3', *D_RISE[0][1:], False)])
    assert [t for t, *_ in D_FALL] == [t for t, *_ in D_RISE] == [3.0, 2.0, 1.0]
    for q, (txt, f, fc, b) in enumerate(key):
        j0 = T.ncol - 2 * len(key) + 2 * q
        kf = {'align': 'center', 'font_size': 8, 'font_color': fc, 'bold': b or None,
              'border': 1, 'border_color': RULE_GREY, **({'bg_color': f, 'pattern': 1} if f else {})}
        ws.merge_range(kr, j0, kr, j0 + 1, txt, sh.fmt(**kf))
    ws.freeze_panes(T.row0, 1)
    n_tests = int((LV.k > 0).sum())
    notes = [
        ('a', [('b', 'Grey row'), ': the gene, its class (sheet ' + SH_SIG + ') and the variants each set tested per '
               'cohort (one number when equal); ‡, a cohort tested the gene on two variants in that set, so each test '
               'without one is a test of the single remaining variant. ', ('b', 'Variant rows'), ': every variant either '
               'set tests in any cohort, in the order of sheet ' + SH_VAR + ', as HGVS protein change (coding change '
               'where there is none).']),
        ('b', [('b', 'Variant'), ': GRCh38, CHROM:POS:REF:ALT (ALT, the effect allele). ', ('b', 'rsID'), ': as in '
               'sheet ' + SH_VAR + '.']),
        ('c', [('b', 'Consequence, Impact'), ': the variant\'s most severe snpEff term in the tested gene and its '
               'impact, as in sheet ' + SH_VAR + ' (same colours).']),
        ('d', [('b', 'Carriers'), ': samples with the minor allele among the cohort\'s cases / controls (minor allele '
               'defined on the pooled samples, as in the tests). Grey, a variant the cohort\'s tests did not use.']),
        ('e', [('b', 'Tier'), ', of the gene in that set, as in sheet ' + SH_SIG + ': 1–3; S3, tier 3 in the '
               f'sensitivity analysis only (genes with ≥ {sens_min} variants); –, not tiered. A tiered set has its '
               '−log10 P in bold on the grey row.']),
        ('f', [('b', '−log10 P'), ' of the gene test. Grey row: with all the gene\'s variants in that set, the '
               'published P (sheet ' + SH_SIG + '), reproduced exactly. Variant rows: the same test with that variant '
               f'left out: the same samples, sex and {n_pcs} principal components as covariates, missing genotypes '
               'imputed to the mean; SKAT-O with Beta(1, 25) weights; CMC, score test.']),
        ('', [('b', 'Shading'), ' (key above the table): the change from the grey row, in steps of 1 (a tenfold change '
              'in P). Red, −log10 P falls without the variant: the variant contributes to the association (bold). '
              'Blue, it rises: the association is stronger without the variant. Darker, a larger change (1 to < 2, 2 to '
              '< 3, ≥ 3); none, a change below 1.']),
        ('', cohort_parts('Cohorts')),
        *missing_block([
            ('—', '−log10 P', ['the cohort\'s test did not use the variant in that set. The Carriers cell of the cohort '
                               f'shows why: grey count (e.g. 0 / 1), minor allele count below {mac_txt} there; grey '
                               '0 / 0, no carrier there; —, not among the cohort\'s rare variants. Reason per cohort: '
                               'sheet ' + SH_VAR + ', column Used in test.']),
            ('—', 'Carriers', [f'not among the cohort\'s rare variants: common there (control MAF ≥ {a.rare_maf:g}) or '
                               'removed by the depth-gap filter.']),
            ('—', 'rsID', ['none in ToMMo or dbSNP 151.']),
            ('Empty', 'the M/H columns', ['a low-impact variant, which the M/H set does not test.'])]),
        ('', [('b', 'Descriptive'), f': {n_tests:,} tests, with no multiple-testing correction; their P are not '
              'compared with the significance thresholds of the scans.']),
        abbrev({**{k: ABBR[k] for k in ('CMC', 'MAC', 'MAF', 'SKAT-O')}, **ABBR_SET}),
    ]
    rr, lens = write_notes(sh, r + 1, notes, 150)
    stats['lovo'] = {'rows': last - T.row0 + 1, 'genes': len(g_order), 'variant_rows': n_rows_v, 'tests': n_tests,
                     'gene_sets_tiered': n_tier, 'genes_two_variants': n_two,
                     'shading': {k or 'none': v for k, v in n_bins.items()}, 'head_rows': head_rows,
                     'max_note': max(lens)}
    print_setup(ws, rr - 1, T.ncol, f'{label}Leave-one-variant-out', repeat=(2, 4))

    # ════════════════════════════════════════════════════════════════════════
    # 2. Gene-based variants (gene x variant), in the columns of the common-variant table
    # ════════════════════════════════════════════════════════════════════════
    sh = sheets[SH_VAR]
    ws = sh.ws
    W = VL.drop_duplicates(['set_name', 'variant_id']).copy()
    W = W.merge(genes[['set_name', 'gene_order']], on='set_name', suffixes=('', '_g'))
    W = W.sort_values(['gene_order', 'pos', 'variant_id']).reset_index(drop=True)
    assert lovo_order == list(zip(W.set_name, W.variant_id)), 'Leave-one-variant-out lists the variants of this sheet'
    by_vc = VL.set_index(['set_name', 'variant_id', 'cohort'])
    n_rows = len(W)
    title(sh, f'{label}Variants of the gene-based signals in three nested ancestry-filtered cohorts',
          f'{len(genes)} genes · {n_rows} variants (' + ', '.join(f'{int((W.impact == k).sum())} {k.lower()}'
                                                                 for k in ('HIGH', 'MODERATE', 'LOW')
                                                                 if (W.impact == k).any())
          + ' impact) · GRCh38 · genes in the order of sheet Gene-based signals, variants in genomic order')
    spec = [dict(blk='var', sub=None, leaf='Tested gene', note='a', w=11),
            dict(blk='var', sub=None, leaf='Class', note='a', w=10),
            dict(blk='var', sub=None, leaf='Set', note='a', w=6.5),
            dict(blk='var', sub=None, leaf='Variant', note='b', w=19),
            dict(blk='var', sub=None, leaf='rsID', note='b', w=13),
            dict(blk='var', sub=None, leaf='Chr', w=4.5), dict(blk='var', sub=None, leaf='Position', w=11),
            dict(blk='var', sub=None, leaf='EA/OA', w=9), dict(blk='var', sub=None, leaf='Gene', note='b', w=12),
            dict(blk='var', sub=None, leaf='Used in', note='a', w=8),
            dict(blk='fun', sub=None, leaf='Consequence', w=30), dict(blk='fun', sub=None, leaf='Impact', w=9),
            dict(blk='fun', sub=None, leaf='Transcript', w=30), dict(blk='fun', sub=None, leaf='HGVS', w=27),
            dict(blk='fun', sub=None, leaf='Other consequences', w=40)]
    if ot_on:
        spec += [dict(blk='ot', sub=None, leaf='ClinVar', w=30), dict(blk='ot', sub=None, leaf='GWAS credible sets', w=60),
                 dict(blk='ot', sub=None, leaf='L2G', w=13), dict(blk='ot', sub=None, leaf='molQTL credible sets', w=56),
                 dict(blk='ot', sub=None, leaf='Enhancer-to-gene', w=44)]
    AFW = 12.5
    af_spec = [('This study', 'CTEPH', 'af_cteph'), ('This study', 'AGP3K', 'af_agp3k'), ('Japanese', 'ToMMo', 'tommo_af'),
               ('Japanese', '1KG JPT', 'kg_jpt_af'), ('East Asian', '1KG EAS', 'kg_eas_af'),
               ('East Asian', 'gnomAD EAS', 'gnomad_eas_af'), ('European', 'gnomAD NFE', 'gnomad_nfe_af')]
    spec += [dict(blk='af', sub=sb, leaf=lf, w=AFW, key=k) for sb, lf, k in af_spec]
    cohort_cols = [(None, 'Used in test', 'g', 17), ('Genotypes', 'CTEPH', 'h', 14), ('Genotypes', 'AGP3K', 'h', 16),
                   ('EAF', 'CTEPH', None, 7), ('EAF', 'AGP3K', None, 7), ('Missing', 'CTEPH', 'h', 7),
                   ('Missing', 'AGP3K', 'h', 7), ('HWE P', 'CTEPH', 'h', 8), ('HWE P', 'AGP3K', 'h', 8)]
    for c in cohorts:
        spec += [dict(blk=c, sub=sb, leaf=lf, note=nt, w=wd, coh=c) for sb, lf, nt, wd in cohort_cols]
    span = {'var': ([('b', 'Variant'), ('sup', 'a')], INK, None),
            'fun': ([('b', 'Functional annotation'), ('sup', 'c')], INK, None),
            'ot': ([('b', 'Open Targets'), ('sup', 'd')], INK, None),
            'af': ([('b', 'Effect-allele frequency'), ('sup', 'e')], INK, None)}
    for c in cohorts:
        span[c] = ([(colors[c], labels[c]), ('sup', 'f'), f'   {n_case[c]:,} cases / {n_ctrl[c]:,} controls'],
                   colors[c], c)
    T = Table(sh, spec, 2, span, colors, tints)
    J = {s['leaf']: j for j, s in enumerate(spec) if isinstance(s['leaf'], str) and s['blk'] in ('var', 'fun', 'ot')}
    j_af = [j for j, s in enumerate(spec) if s['blk'] == 'af']
    j_coh = {c: next(j for j, s in enumerate(spec) if s.get('coh') == c) for c in cohorts}
    FROZEN = J['Used in'] + 1
    r = T.row0
    last = r + n_rows - 1
    gene_seen, n_used = set(), {c: 0 for c in cohorts}
    for _, w in W.iterrows():
        sn, vid = w.set_name, w.variant_id
        first = sn not in gene_seen
        gene_seen.add(sn)
        is_last = r == last

        def put(j, value, align='center', num=None, **kw):
            if blank(value):
                ws.write_string(r, j, MISSING, sh.fmt(**T.rules(j, {'align': 'center', 'font_color': FAINT},
                                                              group_top=first and r > T.row0, last=is_last)))
                return
            d = T.rules(j, {'align': align, 'num_format': num}, group_top=first and r > T.row0, last=is_last)
            d.update(kw)
            if isinstance(value, str):
                ws.write_string(r, j, value, sh.fmt(**d))
            else:
                ws.write_number(r, j, float(value), sh.fmt(**d))

        def rich(j, parts, align='left', **kw):
            sh.rich(r, j, parts, T.rules(j, {'align': align}, group_top=first and r > T.row0, last=is_last), **kw)

        written = set()
        if first:
            rich(0, [('i', w.gene_symbol)], b=True)
            put(2 - 1, cls[sn], align='left', bold=True if cls[sn].startswith('Tier') else None)
        else:
            for j in (0, 1):
                ws.write_blank(r, j, None, sh.fmt(**T.rules(j, {}, last=is_last)))
        put(J['Set'], stratum_label(narrow_s if w.impact != 'LOW' else widest, short=True))
        put(J['Variant'], vid, align='left')
        if s_(w.rsID):
            put(J['rsID'], w.rsID, align='left')
        elif s_(w.rsID_unresolved):
            rich(J['rsID'], [w.rsID_unresolved, ('sup', '‡')])
        else:
            put(J['rsID'], None)
        put(J['Chr'], int(w.chrom) if str(w.chrom).isdigit() else str(w.chrom), num='0')
        put(J['Position'], int(w.pos), num='#,##0', align='right')
        put(J['EA/OA'], f'{w.EA}/{w.OA}')
        rich(J['Gene'], [('i', s_(w.Gene_EnsDb))]) if s_(w.Gene_EnsDb) else put(J['Gene'], None)
        used = [c for c in cohorts if bool(by_vc.loc[(sn, vid, c)][f'in_set_{widest}'])]
        if used:
            parts = []
            for k, c in enumerate(used):
                parts += ([' '] if k else []) + [(colors[c], codes[c])]
            rich(J['Used in'], parts, align='center')
        else:
            put(J['Used in'], '–', font_color=GREY)
        # functional annotation (the tested gene's consequence; other terms over all transcripts)
        put(J['Consequence'], so(w.effect), align='left')
        c_, b_ = IMPACT_STYLE.get(w.impact, (INK, False))
        put(J['Impact'], w.impact.title(), font_color=c_, bold=b_ or None)
        assert lovo_annot[(sn, vid)] == (so(w.effect), w.impact.title()), (sn, vid, 'as in Leave-one-variant-out')
        fid, bt = s_(w.transcript), s_(w.transcript_biotype)
        if fid:
            tail = f' ({so(bt)})' if bt and bt != 'protein_coding' else ''
            rich(J['Transcript'], [('i', s_(w.tx_gene) or w.gene_symbol), f' · {fid}{tail}'])
        else:
            put(J['Transcript'], None)
        hc, hp = s_(w.hgvs_c), s_(w.hgvs_p)
        put(J['HGVS'], (hc + (f' ({hp})' if hp else '')) or None, align='left')
        oc = s_(w.other_consequences)
        if oc:
            put(J['Other consequences'], ', '.join(so(t) for t in oc.split(',')), align='left', font_color=INK2)
        else:
            ws.write_blank(r, J['Other consequences'], None,
                           sh.fmt(**T.rules(J['Other consequences'], {}, group_top=first and r > T.row0, last=is_last)))
        # Open Targets: the number of records, then up to OtTopN of them (as the common-variant table)
        if ot_on:
            jo = J['ClinVar']
            if not is_true(w.ot_in_index):
                for jj in range(jo, jo + 5):
                    put(jj, None)
            else:
                n_f, n_t = int(w.ot_gwas_focus_n), int(w.ot_qtl_trans_n)
                cells = {jo: (int(w.ot_clinvar_n), '', s_(w.ot_clinvar)),
                         jo + 1: (int(w.ot_gwas_n), f' ({n_f} thrombosis)' if n_f else '', s_(w.ot_gwas_top)),
                         jo + 3: (int(w.ot_qtl_n), f' ({n_t} trans)' if n_t else '', s_(w.ot_qtl_top)),
                         jo + 4: (int(w.ot_e2g_n), '', s_(w.ot_e2g_top))}
                for jj, (n, extra, top_) in cells.items():
                    put(jj, f'{n}{extra} · {top_}' if n else None, align='left')
                if s_(w.ot_l2g):
                    g_, sc_ = w.ot_l2g.split(' ', 1)
                    rich(jo + 2, [('i', g_), ' ' + sc_])
                else:
                    put(jo + 2, None)
        # frequencies
        for jj in j_af:
            key = spec[jj]['key']
            val = w[key]
            kw = {'italic': True, 'font_color': GREY} \
                if (key == 'tommo_af' and pd.notna(val) and s_(w.tommo_filter) != 'PASS') else {}
            put(jj, None if pd.isna(val) else float(val), num=NUM_AF, align='right', **kw)
        if a.af_bars:
            ws.conditional_format(r, j_af[0], r, j_af[-1],
                                  {'type': 'data_bar', 'bar_color': '#C9C9C9', 'bar_solid': True,
                                   'bar_no_border': True, 'min_type': 'num', 'min_value': 0, 'max_type': 'max'})
        # cohort blocks
        for c in cohorts:
            x = by_vc.loc[(sn, vid, c)]
            j0 = j_coh[c]
            why = s_(x[f'set_reason_{widest}'])
            if x[f'in_set_{widest}']:
                n_used[c] += 1
                put(j0, 'Yes')
            else:
                mac = '' if pd.isna(x.mac) else int(x.mac)
                maf = '' if pd.isna(x.ctrl_maf_split) else f'{x.ctrl_maf_split:.4f}'
                put(j0, REASON_TEXT[why].format(mac=mac, maf=maf), italic=True, font_color=GREY)
            if not x.in_callset:
                for q in range(1, 9):
                    put(j0 + q, None)
                continue
            gf = '{}/{}/{}/{}'
            put(j0 + 1, gf.format(int(x.case_homRef), int(x.case_het), int(x.case_homAlt), int(x.case_missing_n)))
            put(j0 + 2, gf.format(int(x.ctrl_homRef), int(x.ctrl_het), int(x.ctrl_homAlt), int(x.ctrl_missing_n)))
            put(j0 + 3, float(x.case_alt_af), num=NUM_AF, align='right')
            put(j0 + 4, float(x.ctrl_alt_af), num=NUM_AF, align='right')
            put(j0 + 5, float(x.case_missing_rate), num='0.0%', align='right')
            put(j0 + 6, float(x.ctrl_missing_rate), num='0.0%', align='right')
            for jj, hv in ((j0 + 7, x.case_hwe_p), (j0 + 8, x.ctrl_hwe_p)):
                hp_ = p_parts(float(hv)) if not blank(hv) else None
                if hp_:
                    rich(jj, hp_, align='right')
                else:
                    put(jj, None)
        r += 1
    assert r - 1 == last
    ws.freeze_panes(T.row0, FROZEN)
    tom_lbl = a.tommo_label
    notes = [
        ('a', [('b', 'Tested gene, Class'), ': as in sheet ' + SH_SIG + ', on the gene\'s first row. ', ('b', 'Set'),
               ': M/H, a moderate- or high-impact variant, tested in both sets; L/M/H, a low-impact variant, tested in '
               'the L/M/H set only. ', ('b', 'Used in'), ': the cohorts whose gene test used the variant. Rows: every '
               'variant any cohort tested; genes in the order of sheet ' + SH_SIG + ', variants in genomic order.']),
        ('b', [('b', 'Variant'), ': GRCh38, CHROM:POS:REF:ALT; EA/OA, effect allele (ALT) / other allele (REF). ',
               ('b', 'rsID'), f': the ToMMo {tom_lbl} record\'s; where it lists several (a multi-allelic position), '
               'the one dbSNP 151 assigns to this allele; for a variant absent from ToMMo, dbSNP 151\'s'
               + ('; ‡, ToMMo lists several and dbSNP 151 assigns none of them to this allele (all shown)'
                  if (W.rsID_unresolved != '').any() else '') + '. ', ('b', 'Gene'), ': the Ensembl 86 gene '
               'overlapping the variant, else the nearest within 1 Mb.']),
        ('c', [('b', 'Functional annotation'), ' (snpEff, in the tested gene): Consequence and Impact, the gene\'s '
               'most severe term, as in the gene map; Transcript and HGVS, its most severe transcript (protein-coding '
               'first); Other consequences, every other term over all transcripts and features at the variant. '
               '"sequence feature" marks a protein-feature annotation over the position, not a coding change; '
               '"structural interaction variant" marks a residue in contact within a protein structure, which snpEff '
               'ranks HIGH (the coding change is then under Other consequences).']),
    ]
    if ot_on:
        top = {2: 'two', 3: 'three', 4: 'four', 5: 'five'}.get(otm.get('top_n'), str(otm.get('top_n')))
        cw, td = otm.get('cis_window', 1_000_000) / 1e6, otm.get('trans_distance', 5_000_000) / 1e6
        notes.append(('d', [('b', 'Open Targets'), f' Platform (release {otm["data_version"]}), per variant. Each '
                            f'cell: the number of records · up to {top} of them. ', ('b', 'ClinVar'), ': clinical '
                            'significance: condition, alphabetical. ', ('b', 'GWAS credible sets'), ': sets containing '
                            'the variant (in brackets, those with a thrombosis-related trait: thrombosis, embolism, '
                            'venous, pulmonary hypertension or artery, coagulation, fibrin, haemostasis) · the traits of '
                            'the sets with the highest posterior probability (PP) of the variant, with that PP. ',
                            ('b', 'L2G'), ': the gene most often ranked first by locus-to-gene across those sets (median '
                            'score). ', ('b', 'molQTL credible sets'), ': sets (in brackets, those in trans) · the '
                            'records with the highest PP (QTL type with a cis- or trans- prefix, gene, tissue / cell '
                            f'type, PP); cis, within {cw:g} Mb of the gene\'s transcription start site; trans, another '
                            f'chromosome or more than {td:g} Mb away. ', ('b', 'Enhancer-to-gene'), ': target genes · '
                            'those with the highest score (score, best tissue / cell type, number of tissues / cell '
                            f'types predicting the link). Every record: sheet {SH_OT}, same variants in the same '
                            'order.']))
    notes += [
        ('e', [('b', 'Effect-allele frequency'), f' (ALT). This study: {labels[a.af_cohort].lower()} cases (CTEPH) '
               f'and controls (AGP3K). ToMMo, {tom_lbl}; 1KG, 1000 Genomes {a.kg_label}, JPT (n = {a.kg_n_jpt}) and '
               f'EAS (n = {a.kg_n_eas}); gnomAD, {a.gnomad_label}, EAS and NFE. Grey italic, a ToMMo record failing '
               'its filters.' + (' Bars scaled to the row maximum.' if a.af_bars else '')]),
        ('f', cohort_parts('Cohort blocks')),
        ('g', [('b', 'Used in test'), ': Yes, the cohort\'s gene test used the variant; otherwise the reason. MAC k, '
               f'minor allele count k below the cohort\'s threshold ({mac_txt}); no carrier, no minor allele among '
               f'the cohort\'s samples; common (MAF x), control MAF x ≥ {a.rare_maf:g} in that cohort; depth-gap '
               'filter, removed by the rare-variant QC (sheet ' + SH_QC + ').']),
        ('h', [('b', 'Genotypes'), ': OA/OA / OA/EA / EA/EA / missing. ', ('b', 'Missing'), ': share of the group '
               'not called. ', ('b', 'HWE P'), ': exact test within the group.']),
        *missing_block([
            ('—', 'rsID', ['none in ToMMo or dbSNP 151.']),
            ('—', 'Gene', ['no Ensembl 86 gene within 1 Mb.']),
            ('—', 'Transcript and HGVS', ['no transcript annotation.'])]
            + ([('—', 'the Open Targets columns', ['no record of that kind; in all five, the variant is not in the '
                                                   'Open Targets variant index.'])] if ot_on else [])
            + [('—', 'the frequencies', ['absent from the reference panel; for This study, not among the '
                                         f'{labels[a.af_cohort].lower()} cohort\'s rare variants.']),
               ('—', 'a cohort block (Genotypes to HWE P)', ['the variant is not in the cohort\'s rare-variant call set '
                                                              '(Used in test: no carrier, common or depth-gap filter).']),
               ('–', 'Used in', ['no cohort used the variant.'])]),
        abbrev({k: ABBR[k] for k in ('AGP3K', 'CTEPH', 'EA/OA', 'EAF', 'EAS', 'HWE', 'JPT', 'L2G', 'MAC', 'MAF',
                                     'NFE', 'QTL', 'M/H', 'L/M/H')}),
    ]
    width = int((sum(spec[k]['w'] for k in range(FROZEN)) * 7 + 5) / 5.8)
    rr, lens = write_notes(sh, last + 2, notes, width)
    stats['variants'] = {'rows': n_rows, 'used': n_used, 'max_note': max(lens), 'note_width': width}
    ws.set_landscape()
    ws.set_paper(8)
    ws.set_print_scale(55)
    ws.repeat_rows(2, 4)
    ws.repeat_columns(0, FROZEN - 1)
    ws.set_v_pagebreaks([j_coh[c] for c in cohorts])
    ws.print_area(0, 0, rr - 1, T.ncol - 1)
    ws.set_margins(left=0.4, right=0.4)
    ws.set_footer(f'&L{label}Gene-based variants&RPage &P of &N')

    # ════════════════════════════════════════════════════════════════════════
    # 3. Gene-based Open Targets (every record, as the common-variant table's sheet)
    # ════════════════════════════════════════════════════════════════════════
    n_det, variant_rows, gene_rows = 0, [], []
    if ot_on:
        s2 = sheets[SH_OT]
        w2 = s2.ws
        cols2 = [('Gene', 12, None), ('Set', 7, None), ('Variant', 22, None), ('rsID', 13, None),
                 ('Source', 9, 'a'), ('Type', 9, 'a'), ('Trait, condition or target gene', 46, None),
                 (QTL_COL, 18, 'd'), ('QTL tissue / cell type', 26, 'e'), ('Study', 15, 'b'),
                 ('Publication', 22, 'b'), ('Study N', 9, 'b'), ('Cohorts', 20, 'b'), (POP_COL, 20, 'b'),
                 ('P', 10, 'c'), ('Beta', 8, 'c'), ('P / Beta from', 22, 'c'), ('Set lead variant', 32, 'c'),
                 ('PP', 7, 'e'), ('Fine-mapping (set)', 11, 'e'), ('L2G gene (set)', 11, 'e'),
                 ('L2G score (set)', 9, 'e'), ('E2G score', 8, 'f'), ('E2G best tissue / cell type', 26, 'f'),
                 ('E2G tissues / cell types (n)', 11, 'f'), ('E2G interval', 9, 'f'),
                 ('Clinical significance', 22, None), ('Review status', 30, None)]
        J2 = {h: k for k, (h, _wd, _n) in enumerate(cols2)}
        for k, (_h, wd, _n) in enumerate(cols2):
            w2.set_column(k, k, wd)
        s2.rich(0, 0, [('b', f'{label}Open Targets Platform records for the variants in sheet {SH_VAR} '
                              f'(data release {otm["data_version"]})')], {'align': 'left'}, sz=11)
        s2.rich(1, 0, [f'One row per variant of sheet {SH_VAR}, in the same order, grouped under a header row per '
                       'gene, each with its Open Targets records grouped below it. Outline buttons: 1, genes; 2, all '
                       'variants (as opened); 3, every record.'], {'align': 'left'}, i=True, color=INK2)
        H = OT_H
        for k, (h, _wd, note) in enumerate(cols2):
            s2.rich(H, k, [('b', h), ('sup', note or '')], {'align': 'center', 'valign': 'bottom', 'text_wrap': True,
                                                           'top': MEDIUM, 'top_color': INK, 'bottom': MEDIUM,
                                                           'bottom_color': INK}, b=True)
        w2.set_row(H, 39)
        d = otd.copy()
        d['_c'] = d.category.map({k: n for n, k in enumerate(OT_SOURCES)})
        d['_k'] = np.where(d.category.isin(['GWAS', 'molQTL']), -d.pp.fillna(-1),
                           np.where(d.category == 'E2G', -d.e2g_score.fillna(0), 0.0))
        cv = d.category == 'ClinVar'
        pair = list(zip(d.loc[cv, 'clinical_significance'].replace('', 'not provided'),
                        d.loc[cv, 'trait_or_target'].replace('', 'not specified')))
        rank = {p: n for n, p in enumerate(sorted(set(pair)))}
        d.loc[cv, '_k'] = [rank[p] for p in pair]
        d = d.sort_values(['_c', '_k'], kind='stable')
        by_var = {vid: g for vid, g in d.groupby('variant_id', sort=False)}
        rr2, gene_prev = H + 1, None
        for _, wv in W.iterrows():
            vid = wv.variant_id
            recs = by_var.get(vid, d.iloc[0:0])
            if wv.set_name != gene_prev:
                lv = W[W.set_name == wv.set_name]
                n_idx = int(lv.ot_in_index.map(is_true).sum())
                hfmt = dict({'top': THIN, 'top_color': RULE_GREY} if gene_prev is not None else {},
                            bg_color=LOCUS_FILL, pattern=1, align='left')
                s2.rich(rr2, 0, [('i', wv.gene_symbol)], hfmt, b=True)
                w2.write_blank(rr2, 1, None, s2.fmt(**hfmt))
                s2.rich(rr2, 2, [f'{len(lv)} variant{"s" if len(lv) != 1 else ""} · {n_idx} in the Open Targets '
                                 'variant index'], hfmt, color=INK2)
                for k in range(3, len(cols2)):
                    w2.write_blank(rr2, k, None, s2.fmt(**hfmt))
                w2.set_row(rr2, None, None, {'level': 0})
                gene_rows.append((rr2, wv.set_name))
                rr2 += 1
            vfmt = dict(bg_color=VARIANT_FILL, pattern=1, align='left')
            w2.write_blank(rr2, 0, None, s2.fmt(**vfmt))
            s2.rich(rr2, 1, [stratum_label(narrow_s if wv.impact != 'LOW' else widest, short=True)], vfmt)
            s2.rich(rr2, 2, [vid], vfmt)
            s2.rich(rr2, 3, [s_(wv.rsID)], vfmt)
            if not is_true(wv.ot_in_index):
                s2.rich(rr2, 4, [('i', 'Not in the Open Targets variant index')], vfmt, color=GREY)
            elif not len(recs):
                s2.rich(rr2, 4, [('i', 'In the Open Targets variant index; no record')], vfmt, color=GREY)
            else:
                n = recs.category.value_counts()
                n_tr = int((recs.qtl_class == 'trans').sum())
                s2.rich(rr2, 4, [' · '.join(f'{k} {int(n.get(k, 0))}' + (f' ({n_tr} trans)' if k == 'molQTL' and n_tr
                                                                          else '') for k in OT_SOURCES)], vfmt)
            for k in range(5, len(cols2)):
                w2.write_blank(rr2, k, None, s2.fmt(**vfmt))
            w2.set_row(rr2, None, None, {'level': 1, 'collapsed': bool(len(recs))})
            variant_rows.append((rr2, vid))
            gene_prev = wv.set_name
            rr2 += 1
            for _, x in recs.iterrows():
                def p2(k, val, num=None, align='left', **kw):
                    f = dict(align=align, num_format=num, **kw)
                    if val is None or (isinstance(val, float) and not np.isfinite(val)) or val == '':
                        w2.write_blank(rr2, k, None, s2.fmt(**f))
                    elif isinstance(val, str):
                        w2.write_string(rr2, k, val, s2.fmt(**f))
                    else:
                        w2.write_number(rr2, k, float(val), s2.fmt(**f))
                cs_rec = x.category in ('GWAS', 'molQTL')
                p2(J2['Gene'], '')
                p2(J2['Set'], '')
                p2(J2['Variant'], vid, font_color=FAINT)
                p2(J2['rsID'], '')
                p2(J2['Source'], x.category)
                p2(J2['Type'], x.type if x.type != x.category and x.category != 'ClinVar' else '')
                jt = J2['Trait, condition or target gene']
                if x.category in ('molQTL', 'E2G'):
                    s2.rich(rr2, jt, [('i', ' '.join(x.trait_or_target.split()))], dict(align='left'))
                else:
                    p2(jt, ' '.join(x.trait_or_target.split()), bold=True if x.focus else None)
                p2(J2[QTL_COL], qtl_text(x) if x.category == 'molQTL' else '')
                p2(J2['QTL tissue / cell type'], x.biosample if x.category == 'molQTL' else '')
                p2(J2['Study'], x.study_id)
                p2(J2['Publication'], f'{x.first_author} {x.year}'.strip() if x.first_author else '')
                p2(J2['Study N'], x.n_samples if pd.notna(x.n_samples) else None, num='#,##0', align='right')
                if not cs_rec:
                    p2(J2['Cohorts'], '')
                elif x.cohorts:
                    p2(J2['Cohorts'], x.cohorts)
                elif x.category == 'molQTL' and x.project:
                    p2(J2['Cohorts'], x.project, italic=True, font_color=GREY)
                else:
                    p2(J2['Cohorts'], '—')
                p2(J2[POP_COL], (x.ld_population or '—') if cs_rec else '')
                lead_p, lead_b = x.p_source == 'lead', x.beta_source == 'lead'
                pp_ = p_parts(x.p_value) if pd.notna(x.p_value) else None
                if pp_:
                    s2.rich(rr2, J2['P'], pp_, dict(align='right'), i=lead_p, color=GREY if lead_p else INK)
                else:
                    p2(J2['P'], None)
                p2(J2['Beta'], x.beta if pd.notna(x.beta) else None, num='0.000', align='right',
                   **({'italic': True, 'font_color': GREY} if lead_b else {}))
                p2(J2['P / Beta from'], source_text(x.p_source, x.beta_source) if cs_rec else '')
                if cs_rec:
                    p2(J2['Set lead variant'], 'this variant' if x.lead_is_variant else
                       (f'{x.lead_variant} ({x.lead_rsid.replace(";", "/")})' if x.lead_rsid else x.lead_variant))
                else:
                    p2(J2['Set lead variant'], '')
                p2(J2['PP'], x.pp if pd.notna(x.pp) else None, num='0.000', align='right')
                p2(J2['Fine-mapping (set)'], x.method)
                if x.l2g_gene:
                    s2.rich(rr2, J2['L2G gene (set)'], [('i', x.l2g_gene)], dict(align='left'))
                else:
                    p2(J2['L2G gene (set)'], '')
                p2(J2['L2G score (set)'], x.l2g_score if pd.notna(x.l2g_score) else None, num='0.00', align='right')
                p2(J2['E2G score'], x.e2g_score if pd.notna(x.e2g_score) else None, num='0.00', align='right')
                p2(J2['E2G best tissue / cell type'], x.biosample if x.category == 'E2G' else '')
                p2(J2['E2G tissues / cell types (n)'], x.n_biosamples if pd.notna(x.n_biosamples) else None,
                   num='#,##0', align='right')
                p2(J2['E2G interval'], x.interval_type)
                p2(J2['Clinical significance'], x.clinical_significance)
                p2(J2['Review status'], x.review_status)
                w2.set_row(rr2, None, None, {'level': 2, 'hidden': True})
                rr2 += 1
                n_det += 1
        n_ot = sum(len(by_var.get(v, [])) for v in W.variant_id)
        assert len(variant_rows) == n_rows and n_det == n_ot, (len(variant_rows), n_rows, n_det, n_ot)
        assert set(otd.variant_id) <= set(W.variant_id), 'Open Targets records of variants not in the table'
        stats['ot_records'] = n_det
        w2.outline_settings(True, False, True, False)
        for k in range(len(cols2)):
            w2.write_blank(rr2, k, None, s2.fmt(top=MEDIUM, top_color=INK))
        ot_notes = [
            ('a', [('b', 'Records'), ': ClinVar, evidence on the variant; GWAS and molQTL, credible sets of GWAS and '
                   'molecular-QTL studies containing the variant (Type, the QTL type); E2G, enhancer-to-gene '
                   'predictions, one row per target gene. Under each variant: ClinVar alphabetical, GWAS and molQTL by '
                   f'PP, E2G by score; the first ones are those named in sheet {SH_VAR}. Thrombosis-related GWAS traits '
                   'in bold. ', ('b', 'Set'), ': as in sheet ' + SH_VAR + '.']),
            ('b', [('b', 'Study'), ', the Open Targets study identifier; ', ('b', 'Publication'), ', first author and '
                   'year; ', ('b', 'Study N'), ', sample size; ', ('b', 'Cohorts'), ', as listed by Open Targets (grey '
                   'italic: the study\'s project, for studies listing none, i.e. molQTL studies); ', ('b', 'Population'),
                   ', share of the study sample by gnomAD reference population, as Open Targets assigns them for LD.']),
            ('c', [('b', 'P, Beta'), ': the variant\'s own in that study where reported, otherwise those of the '
                   'credible set\'s lead variant (grey italic; source in P / Beta from, lead in Set lead variant).']),
            ('d', [('b', 'cis / trans'), f': the variant against the molQTL target gene. cis, within {cw:g} Mb of the '
                   'gene\'s TSS (start of its Ensembl canonical transcript, end on the minus strand); trans, another '
                   f'chromosome or more than {td:g} Mb away; neither, {cw:g}–{td:g} Mb. Open Targets\' own flag (a '
                   'credible set is trans when every variant in it is) is added where it disagrees.']),
            ('e', [('b', 'PP'), ', posterior probability of the variant in the credible set. ', ('b', 'Fine-mapping, '
                   'L2G gene, L2G score'), ': properties of the credible set (method; its top locus-to-gene '
                   'prediction), not of the variant. ', ('b', 'QTL tissue / cell type'), ': that of the molQTL study.']),
            ('f', [('b', 'E2G'), ': per target gene, the best score, the tissue / cell type giving it, the number of '
                   'tissues / cell types predicting the link, and the type of the predicted interval.']),
            *missing_block([('—', 'Cohorts', ['none listed by Open Targets.']),
                            ('—', 'Population', ['not reported.']),
                            ('—', 'P / Beta from', ['neither the variant\'s nor the lead\'s P and Beta reported (P and '
                                                    'Beta then empty).'])])]
        pops = {p.split(' ')[0] for s in otd.ld_population if s for p in s.split(', ')}
        abbr2 = {'E2G': 'enhancer-to-gene', 'L2G': 'locus-to-gene', 'LD': 'linkage disequilibrium',
                 'PP': 'posterior probability', 'TSS': 'transcription start site',
                 'QTL': 'quantitative trait locus (eQTL, pQTL, sQTL, tuQTL: expression, protein, splicing, '
                        'transcript-usage; sc-, single-cell)', **ABBR_SET,
                 **{p: GNOMAD_POP.get(p, 'gnomAD population') for p in pops}}
        ot_width = int((sum(wd for _h, wd, _n in cols2[:4]) * 7 + 5) / 5.8)
        rn, _l = write_notes(s2, rr2 + 2, ot_notes + [abbrev(abbr2)], ot_width)
        if n_det:
            w2.autofilter(H, 0, rr2 - 1, len(cols2) - 1)
        w2.freeze_panes(H + 1, 4)
        w2.set_landscape()
        w2.set_paper(8)
        w2.repeat_rows(H, H)
        w2.print_area(0, 0, rn - 1, len(cols2) - 1)

    # ════════════════════════════════════════════════════════════════════════
    # Variant sets: each step on one row, its variants and its genes side by side
    # ════════════════════════════════════════════════════════════════════════
    sh = sheets[SH_SETS]
    ws = sh.ws
    assert min_var == 2, 'the gene rows (1 / 2 / >= 3 variants) are laid out for a minimum of 2 variants per gene'
    title(sh, f'{label}From the analysed rare variants to the genes tested',
          'per cohort · each step on one row, with its variants (left) and the genes they fall in (right) · a '
          'variant enters the gene-based tests through the genes it is mapped to')
    F = FLOW.set_index(['cohort', 'section', 'stratum', 'item']).value

    def fv(c, sec, item, s=''):
        return int(F.loc[(c, sec, s, item)])
    all_c = lambda f: {c: f(c) for c in cohorts}                                          # noqa: E731
    imps_of = {narrow_s: ('HIGH', 'MODERATE'), widest: ('HIGH', 'MODERATE', 'LOW')}
    assert set(imps_of[widest]) == {'HIGH', 'MODERATE', 'LOW'}
    qual_of = {s: ' or '.join([', '.join(imps_of[s][:-1]), imps_of[s][-1]]) for s in strata}
    spec = [dict(blk=None, sub=None, leaf='Step', w=46), dict(blk=None, sub=None, leaf='Definition', w=62)]
    spec += [dict(blk='var', sub=None, leaf=[(colors[c], labels[c])], w=12.5) for c in cohorts]
    spec += [dict(blk='gene', sub=None, leaf=[(colors[c], labels[c])], w=12.5) for c in cohorts]
    span = {'var': ([('b', 'Variants'), ('sup', 'f')], INK, None), 'gene': ([('b', 'Genes'), ('sup', 'f')], INK, None)}
    T = Table(sh, spec, 2, span, colors, tints, header_h=26)
    examples = {'HIGH': 'e.g. stop gained, frameshift, splice donor or acceptor',
                'MODERATE': 'e.g. missense, in-frame insertion or deletion',
                'LOW': 'e.g. synonymous, splice region',
                'MODIFIER': 'e.g. intronic, intergenic, UTR'}
    V_ = lambda sec, item, s='': all_c(lambda c: fv(c, sec, item, s))                        # noqa: E731
    rows = [('h', [('b', 'Analysed rare variants'), ('sup', 'a')], None, None, None),
            ('b', 'Analysed variants', 'Passed the rare-variant QC (sheet ' + SH_QC + ')',
             V_('analysed', 'variants'), None),
            ('s', [('i', 'By the highest impact among the variant\'s snpEff annotations'), ('sup', 'b')], None, None,
             None)]
    for imp in ('HIGH', 'MODERATE', 'LOW', 'MODIFIER'):
        rows.append(('r', f'   {imp}', examples[imp], V_('most_severe', imp), None))
    if any(fv(c, 'most_severe', 'none') for c in cohorts):
        rows.append(('r', '   No annotation', 'No snpEff record for the allele', V_('most_severe', 'none'), None))
    not_mapped = all_c(lambda c: fv(c, 'analysed', 'variants') - fv(c, 'assigned', 'variants'))
    mod_only = all_c(lambda c: fv(c, 'most_severe', 'MODIFIER') + fv(c, 'most_severe', 'none'))
    rows += [('h', [('b', 'Mapping to genes'), ('sup', 'c')], None, None, None),
             ('b', 'Not mapped to a gene', '', not_mapped, None),
             ('r', '   MODIFIER annotations only', 'No HIGH, MODERATE or LOW annotation', mod_only, None),
             ('r', '   HIGH, MODERATE or LOW annotation without a gene', 'Only on transcription-factor binding-site '
              'motifs', V_('assigned', 'eligible_impact_no_gene'), None),
             ('b', 'Mapped to genes', 'A HIGH, MODERATE or LOW consequence in a named gene', V_('assigned', 'variants'),
              V_('set', 'genes_with_variant', widest)),
             ('r', '   in more than one gene', 'Overlapping genes: readthrough transcripts, UGT1A, protocadherins',
              V_('assigned', 'multigene_variants'), None)]
    for s in strata:
        q = qual_of[s]
        rows += [('h', [('b', f'{stratum_label(s).replace(" variants", "")} set ({q} consequences)'), ('sup', 'd')],
                  None, None, None),
                 ('b', 'Qualifying', f'A {q} consequence in the gene', V_('set', 'qualifying_variants', s),
                  V_('set', 'genes_with_variant', s)),
                 ('r', '   → not tested', 'The gene has only 1 qualifying variant',
                  V_('set', 'untested_only_variants', s), V_('set', 'genes_one_variant', s)),
                 ('b', '   → tested', 'The gene has 2 or more qualifying variants', V_('set', 'variants_tested', s),
                  V_('set', 'genes_tested', s)),
                 ('r', '         genes with 2 variants', '', None, V_('set', 'genes_2', s)),
                 ('r', '         genes with 3 or more variants', '', None, V_('set', 'genes_ge3', s))]
        for imp in imps_of[s]:
            rows.append(('r', f'         {imp} variants', 'By the most severe consequence among the tested genes',
                         V_('set_impact', imp, s), None))
        rows.append(('r', [('', '         variants in more than one tested gene'), ('sup', 'e')],
                     'Enter the test of each of these genes', V_('set', 'multigene_variants', s), None))
    r = T.row0
    last = r + len(rows) - 1
    j_gene = 2 + len(cohorts)
    for kind, step, rule, vv, gg in rows:
        is_last = r == last
        parts = step if isinstance(step, list) else [step]
        if kind in ('h', 's'):
            for j in range(T.ncol):
                d = T.rules(j, {'align': 'left', 'bg_color': HEADER_FILL if kind == 'h' else None,
                                'pattern': 1 if kind == 'h' else None}, group_top=kind == 'h' and r > T.row0,
                            last=is_last)
                if j == 0:
                    sh.rich(r, 0, parts, d, b=kind == 'h', color=INK if kind == 'h' else INK2, i=kind == 's')
                else:
                    ws.write_blank(r, j, None, sh.fmt(**d))
            r += 1
            continue
        b = kind == 'b'
        sh.rich(r, 0, parts, T.rules(0, {'align': 'left'}, last=is_last), b=b)
        if rule:
            ws.write_string(r, 1, rule, sh.fmt(**T.rules(1, {'align': 'left', 'font_color': INK2}, last=is_last)))
        else:
            ws.write_blank(r, 1, None, sh.fmt(**T.rules(1, {}, last=is_last)))
        for j0, vals in ((2, vv), (j_gene, gg)):
            for k, c in enumerate(cohorts):
                left = {'left': THIN, 'left_color': RULE_GREY} if (j0 == j_gene and k == 0) else {}
                if vals is None:
                    ws.write_blank(r, j0 + k, None, sh.fmt(**T.rules(j0 + k, dict(left), last=is_last)))
                else:
                    ws.write_number(r, j0 + k, vals[c], sh.fmt(**T.rules(j0 + k, dict(left, align='right',
                                                                                       num_format='#,##0',
                                                                                       bold=b or None),
                                                                         last=is_last)))
        r += 1
    # each step adds up, in variants and in genes
    for c in cohorts:
        an = fv(c, 'analysed', 'variants')
        assert sum(fv(c, 'most_severe', i) for i in ('HIGH', 'MODERATE', 'LOW', 'MODIFIER', 'none')) == an
        assert fv(c, 'assigned', 'variants') + mod_only[c] + fv(c, 'assigned', 'eligible_impact_no_gene') == an
        assert fv(c, 'set', 'qualifying_variants', widest) == fv(c, 'assigned', 'variants')
        for s in strata:
            assert fv(c, 'set', 'qualifying_variants', s) == fv(c, 'set', 'untested_only_variants', s) + \
                fv(c, 'set', 'variants_tested', s)
            assert fv(c, 'set', 'genes_with_variant', s) == fv(c, 'set', 'genes_one_variant', s) + \
                fv(c, 'set', 'genes_tested', s)
            assert fv(c, 'set', 'genes_tested', s) == fv(c, 'set', 'genes_2', s) + fv(c, 'set', 'genes_ge3', s)
            assert sum(fv(c, 'set_impact', i, s) for i in imps_of[s]) == fv(c, 'set', 'variants_tested', s)
    ws.freeze_panes(T.row0, 0)
    notes = [
        ('a', [('b', 'Analysed variants'), ': the cohort\'s rare variants after the rare-variant QC (control MAF < '
               f'{a.rare_maf:g}, MAC ≥ {mac_txt}; sheet ' + SH_QC + ').']),
        ('b', [('b', 'Highest impact'), ': snpEff annotates a variant on every transcript and feature it overlaps, each '
               'annotation with a consequence term and an impact (HIGH > MODERATE > LOW > MODIFIER); the variant is '
               'counted under its highest. The terms within each impact: sheet ' + SH_TERMS + '.']),
        ('c', [('b', 'Mapping to genes'), ', as in the gene-based tests: per gene the variant is annotated in, the '
               'gene\'s most severe consequence is kept, and the variant is mapped when it is HIGH, MODERATE or LOW. '
               '"Without a gene": annotations on transcription-factor binding-site motifs name no gene, so a variant '
               'whose only HIGH, MODERATE or LOW annotations are of that kind is not mapped. "More than one gene": '
               'overlapping genes, e.g. readthrough transcripts sharing exons with their parent genes (ARHGAP8 and '
               'PRR5-ARHGAP8) and clusters with shared exons (UGT1A, protocadherins).']),
        ('d', [('b', 'Per variant set'), ': a variant qualifies with a consequence of the set\'s impacts in a gene; a '
               f'gene is tested when it holds {min_var} or more qualifying variants. A variant is tested when at '
               'least one of its genes is, and not tested when all its genes hold a single qualifying variant. Tested '
               f'genes with {min_var} and with ≥ {sens_min} variants: the two groups of λ', ('sub', 'GC'),
               ' in sheet ' + SH_SCANS + '.']),
        ('e', [('b', 'More than one tested gene'), ': such a variant enters the test of each gene; it is counted once '
               'in the Variants columns.']),
        ('f', [('b', 'Variants, Genes'), ': distinct variants; distinct genes (snpEff gene symbols).']),
        ('', cohort_parts('Cohorts')),
        *missing_block([('Empty', 'Variants or Genes', ['the step concerns only genes, or only variants.'])]),
        abbrev({k: ABBR[k] for k in ('MAC', 'MAF')}),
    ]
    rr, lens = write_notes(sh, last + 2, notes, 150)
    stats['variant_sets'] = {'rows': len(rows)}
    print_setup(ws, rr - 1, T.ncol, f'{label}Variant sets', paper=9)
    ws.fit_to_pages(1, 1)

    # ════════════════════════════════════════════════════════════════════════
    # Consequence terms: what each impact holds, per cohort (with totals)
    # ════════════════════════════════════════════════════════════════════════
    sh = sheets[SH_TERMS]
    ws = sh.ws
    title(sh, f'{label}snpEff consequence terms of the analysed rare variants and of the tested genes',
          'per cohort · each impact with its total and its terms; counts among the analysed variants and among the '
          'variants of the tested genes in each variant set')
    spec2 = [dict(blk=None, sub=None, leaf='Impact', w=11), dict(blk=None, sub=None, leaf='Consequence term', note='a',
                                                                    w=80)]
    for c in cohorts:
        spec2 += [dict(blk=c, sub=None, leaf='Analysed', note='b', w=11, coh=c),
                  dict(blk=c, sub='Tested genes', leaf=stratum_label(narrow_s, short=True), note='c', w=9, coh=c),
                  dict(blk=c, sub='Tested genes', leaf=stratum_label(widest, short=True), note='c', w=9, coh=c)]
    span2 = {c: ([(colors[c], labels[c])], colors[c], c) for c in cohorts}
    T2 = Table(sh, spec2, 2, span2, colors, tints, header_h=14)
    terms = TERMS.copy()
    terms['rank'] = terms.impact.map({'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}).fillna(9)
    piv = {c: terms[terms.cohort == c].set_index(['impact', 'effect']) for c in cohorts}
    allt = terms.groupby(['rank', 'impact', 'effect'], as_index=False).analysed.sum() \
        .sort_values(['rank', 'analysed'], ascending=[True, False])
    impacts = list(dict.fromkeys(allt.impact))
    in_set = lambda imp, s: imp in (('HIGH', 'MODERATE') if s == narrow_s else ('HIGH', 'MODERATE', 'LOW'))  # noqa

    def cell_vals(c, imp, effect=None):
        t = terms[(terms.cohort == c) & (terms.impact == imp)]
        if effect is not None:
            t = t[t.effect == effect]
        out = [int(t.analysed.sum())]
        for s in strata:
            out.append(int(t[s].sum()) if in_set(imp, s) else None)
        return out
    r = T2.row0
    n_t = len(allt)
    last2 = r + n_t + len(impacts)          # term rows + one total row per impact + the grand total
    for imp in impacts:
        # the impact's total row, then its terms
        is_last = False
        d0 = dict(bg_color=HEADER_FILL, pattern=1)
        sh.rich(r, 0, [('b', imp)], T2.rules(0, dict(d0, align='left'), group_top=r > T2.row0))
        sh.rich(r, 1, [('b', f'All {imp} terms')], T2.rules(1, dict(d0, align='left'), group_top=r > T2.row0))
        for k, c in enumerate(cohorts):
            for q, v in enumerate(cell_vals(c, imp)):
                d = T2.rules(2 + 3 * k + q, dict(d0, align='right' if v is not None else 'center', num_format='#,##0',
                                                 bold=True if v is not None else None,
                                                 font_color=GREY if v is None else None), group_top=r > T2.row0)
                if v is None:
                    ws.write_string(r, 2 + 3 * k + q, '–', sh.fmt(**d))
                else:
                    ws.write_number(r, 2 + 3 * k + q, v, sh.fmt(**d))
        r += 1
        for _, t in allt[allt.impact == imp].iterrows():
            ws.write_blank(r, 0, None, sh.fmt(**T2.rules(0, {})))
            ws.write_string(r, 1, '   ' + so(t.effect), sh.fmt(**T2.rules(1, {'align': 'left'})))
            for k, c in enumerate(cohorts):
                for q, v in enumerate(cell_vals(c, imp, t.effect)):
                    d = T2.rules(2 + 3 * k + q, {'align': 'right' if v is not None else 'center',
                                                 'num_format': '#,##0', 'font_color': GREY if v is None else None})
                    if v is None:
                        ws.write_string(r, 2 + 3 * k + q, '–', sh.fmt(**d))
                    else:
                        ws.write_number(r, 2 + 3 * k + q, v, sh.fmt(**d))
            r += 1
    # grand total
    d0 = dict(bold=True, top=THIN, top_color=INK)
    sh.rich(r, 0, [('b', 'All')], T2.rules(0, dict(d0, align='left'), last=True))
    sh.rich(r, 1, [('b', 'All analysed variants / all variants of the tested genes')],
            T2.rules(1, dict(d0, align='left'), last=True))
    for k, c in enumerate(cohorts):
        tot = [fv(c, 'analysed', 'variants')] + [fv(c, 'set', 'variants_tested', s) for s in strata]
        for q, v in enumerate(tot):
            ws.write_number(r, 2 + 3 * k + q, v, sh.fmt(**T2.rules(2 + 3 * k + q, dict(d0, align='right',
                                                                                        num_format='#,##0'),
                                                                   last=True)))
    assert r == last2, (r, last2)
    for c in cohorts:                     # the totals are those of sheet Variant sets
        tc = terms[terms.cohort == c]
        assert int(tc.analysed.sum()) == fv(c, 'analysed', 'variants')
        for imp in ('HIGH', 'MODERATE', 'LOW', 'MODIFIER'):
            assert cell_vals(c, imp)[0] == fv(c, 'most_severe', imp), (c, imp)
        for s in strata:
            assert int(tc[s].sum()) == fv(c, 'set', 'variants_tested', s)
            for imp in ('HIGH', 'MODERATE', 'LOW'):
                if in_set(imp, s):
                    assert cell_vals(c, imp)[1 + strata.index(s)] == fv(c, 'set_impact', imp, s), (c, s, imp)
    ws.freeze_panes(T2.row0, 2)
    notes = [
        ('a', [('b', 'Consequence term'), ': snpEff terms as written; several consequences in one annotation form a '
               'combined term (joined by "&"). "sequence feature", a protein-feature annotation over the position, '
               'not a coding change; "structural interaction variant" and "protein protein contact", snpEff\'s '
               'HIGH-impact terms for residues in contact within a protein structure or complex; "gene fusion" and '
               f'"bidirectional gene fusion" name two genes and are never tested (fewer than {min_var} variants under '
               'that name).']),
        ('b', [('b', 'Analysed'), ': each analysed variant once, under the most severe of all its annotations (as in '
               'sheet ' + SH_SETS + '); the impact totals are that sheet\'s impact rows.']),
        ('c', [('b', 'Tested genes'), ': each variant of a tested gene once, under its most severe consequence in the '
               'tested genes of that set; it can sit under another term than in Analysed (e.g. MODERATE in an '
               'untested gene, LOW in a tested one).']),
        ('', cohort_parts('Cohorts')),
        *missing_block([('–', 'Tested genes', ['the impact is not in the set (M/H: HIGH and MODERATE; L/M/H: HIGH, '
                                                'MODERATE and LOW).'])]),
        abbrev({k: ABBR[k] for k in ('M/H', 'L/M/H')}),
    ]
    rr, lens = write_notes(sh, last2 + 2, notes, 150)
    stats['terms'] = {'rows': n_t, 'impacts': len(impacts)}
    print_setup(ws, rr - 1, T2.ncol, f'{label}Consequence terms', paper=9)

    # ════════════════════════════════════════════════════════════════════════
    # Gene-based scans: one row per cohort and variant set; the primary family (genes with >= MinNumVar
    # variants) and the sensitivity family (>= SensMinNumVar) side by side, each with both statistics
    # ════════════════════════════════════════════════════════════════════════
    sh = sheets[SH_SCANS]
    ws = sh.ws
    title(sh, f'{label}Gene-based scans: genes tested, significance thresholds, significant genes and λGC',
          [f'one row per cohort and variant set · primary analysis (genes with ≥ {min_var} variants) and sensitivity '
           f'analysis (genes with ≥ {sens_min} variants) side by side · λ', ('sub', 'GC'), ' as observed; no '
           'genomic-control correction is applied'])
    lam = [('b', 'λ'), ('sub', 'GC')]
    fams = [('prim', f'Primary analysis: genes with ≥ {min_var} variants'),
            ('sens', f'Sensitivity analysis: genes with ≥ {sens_min} variants')]
    spec = [dict(blk=None, sub=None, leaf='Cohort', w=18), dict(blk=None, sub=None, leaf='Variant set', w=26)]
    for fk, _ft in fams:
        spec += [dict(blk=fk, sub=None, leaf='Genes tested', note='a', w=9),
                 dict(blk=fk, sub=None, leaf='Bonferroni threshold', note='b', w=11)]
        for m in methods:
            ml = METHOD_SHORT[m]
            spec += [dict(blk=fk, sub=ml, leaf='Smallest P', w=10),
                     dict(blk=fk, sub=ml, leaf='Bonferroni', note='c', w=8.5),
                     dict(blk=fk, sub=ml, leaf='BH', note='c', w=6),
                     dict(blk=fk, sub=ml, leaf='Largest significant P (BH)', note='d', w=10.5),
                     dict(blk=fk, sub=ml, leaf=lam, note='e', w=7)]
    span = {fk: ([('b', ft)], INK, None) for fk, ft in fams}
    T = Table(sh, spec, 2, span, colors, tints, header_h=39)
    SCi = SC.set_index(['cohort', 'stratum', 'method'])
    j_f = {fk: 2 + 12 * i for i, (fk, _ft) in enumerate(fams)}
    r = T.row0
    last = r + len(cohorts) * len(strata) - 1
    for c in cohorts:
        for si, s in enumerate(strata):
            gtop = si == 0 and r > T.row0
            stop = si > 0
            is_last = r == last
            x0 = SCi.loc[(c, s, methods[0])]

            def put(j, value, align='center', num=None, **kw):
                d = T.rules(j, {'align': align, 'num_format': num}, group_top=gtop, soft_top=stop, last=is_last)
                d.update(kw)
                if isinstance(value, str):
                    ws.write_string(r, j, value, sh.fmt(**d))
                else:
                    ws.write_number(r, j, float(value), sh.fmt(**d))

            def rich(j, parts, align='right', **kw):
                sh.rich(r, j, parts, T.rules(j, {'align': align}, group_top=gtop, soft_top=stop, last=is_last), **kw)
            if si == 0:
                rich(0, [(colors[c], labels[c])], align='left')
            else:
                ws.write_blank(r, 0, None, sh.fmt(**T.rules(0, {}, last=is_last)))
            put(1, stratum_label(s), align='left')
            for fk, _ft in fams:
                j0 = j_f[fk]
                if fk == 'prim':
                    n_g, thr = int(x0.n_genes_mapped), float(x0.bonferroni_threshold)
                else:
                    n_g, thr = int(x0.n_genes_sens), float(x0.bonferroni_threshold_sens)
                put(j0, n_g, num='#,##0', align='right')
                rich(j0 + 1, sci_parts(thr))
                for mi, m in enumerate(methods):
                    x = SCi.loc[(c, s, m)]
                    jm = j0 + 2 + 5 * mi
                    if fk == 'prim':
                        assert int(x.n_genes_lambda) == int(x.n_genes_mapped) == n_g, (c, s, m)
                        mp, nbf, nbh, bht, lm = (float(x.min_p), int(x.n_significant_bonferroni),
                                                 int(x.n_significant_bh), x.bh_threshold, float(x.lambda_gc))
                    else:
                        assert int(x.n_genes_sens) == int(x.n_genes_lambda_ge3) == n_g, (c, s, m)
                        mp, nbf, nbh, bht, lm = (float(x.min_p_sens), int(x.n_significant_bonferroni_sens),
                                                 int(x.n_significant_bh_sens), x.bh_threshold_sens,
                                                 float(x.lambda_gc_nvar_ge3))
                    rich(jm, p_parts(mp))
                    put(jm + 1, nbf, num='0', bold=True if nbf else None, font_color=None if nbf else GREY)
                    put(jm + 2, nbh, num='0', bold=True if nbh else None, font_color=None if nbh else GREY)
                    if nbh:
                        rich(jm + 3, p_parts(float(bht)))
                    else:
                        put(jm + 3, '–', font_color=GREY)
                    put(jm + 4, lm, num='0.000', align='right', bold=True)
            r += 1
    ws.freeze_panes(T.row0, 2)
    lam_rng = {(fk, m): (min(SC[SC.method == m][col].astype(float)), max(SC[SC.method == m][col].astype(float)))
               for fk, col in (('prim', 'lambda_gc'), ('sens', 'lambda_gc_nvar_ge3')) for m in methods}
    notes = [
        ('a', [('b', 'Genes tested'), f': primary analysis, the genes with {min_var} or more qualifying variants '
               f'(MAC ≥ {mac_txt} in the cohort; sheet ' + SH_SETS + '); sensitivity analysis, those with '
               f'{sens_min} or more. The tests and P are the same in both; only the family in which significance is '
               'judged differs.']),
        ('b', [('b', 'Bonferroni threshold'), f': {alpha:g} / genes tested in the family, for both statistics.']),
        ('c', [('b', 'Significant genes'), ' per cohort, variant set, statistic and family: Bonferroni, P below the '
               f'threshold; BH, Benjamini–Hochberg q < {alpha:g} over the family. Every Bonferroni-significant gene is '
               'also BH-significant. No correction across statistics, variant sets or cohorts.']),
        ('d', [('b', 'Largest significant P (BH)'), ': the largest P among the BH-significant genes, i.e. the P '
               'cut-off BH applied in that scan.']),
        ('e', [('b', 'λ'), ('sub', 'GC'), ': median of the 1-df χ² equivalents of the family\'s gene P / 0.455 (its '
               'expected value under no association). Ranges over cohorts and variant sets: '
               + '; '.join(f'{METHOD_SHORT[m]}, {lam_rng[("prim", m)][0]:.3f}–{lam_rng[("prim", m)][1]:.3f} '
                           f'(≥ {min_var} variants) and {lam_rng[("sens", m)][0]:.3f}–{lam_rng[("sens", m)][1]:.3f} '
                           f'(≥ {sens_min} variants)' for m in methods)
               + '. Reported as observed; no P is corrected.']),
        ('', cohort_parts('Cohorts')),
        *missing_block([('–', 'Largest significant P (BH)', [f'no BH-significant gene (every q ≥ {alpha:g}).'])]),
        abbrev({k: ABBR[k] for k in ('BH', 'CMC', 'MAC', 'SKAT-O')}),
    ]
    rr, lens = write_notes(sh, last + 2, notes, 160)
    print_setup(ws, rr - 1, T.ncol, f'{label}Gene-based scans', repeat=(2, 4))

    # ════════════════════════════════════════════════════════════════════════
    # 7. Rare-variant QC (the layout of the project's rare-variant QC summary)
    # ════════════════════════════════════════════════════════════════════════
    write_qc_sheet(wb, sheets[SH_QC], QC, cohorts, labels, colors, tints, n_pcs)
    wb.close()
    stats['ot_on'] = ot_on
    stats['n_det'] = n_det
    stats['variant_rows'] = variant_rows
    stats['gene_rows'] = gene_rows
    return stats


ABBR_SET = {'M/H': ABBR['M/H'], 'L/M/H': ABBR['L/M/H']}


def write_qc_sheet(wb, ws, QC, cohorts, labels, colors, tints, n_pcs):
    """Sheet Rare-variant QC: title, the QC flow and the basis of the minor allele count threshold, footnotes and
    references; ruled cells, navy group labels; title and tables in Arial, footnotes and references in DengXian."""
    C = [QC[c] for c in cohorts]

    def same(xs):
        assert len(set(map(str, xs))) == 1, f'differs between cohorts: {xs}'
        return xs[0]
    KNEE = same([c['variant_qc']['knee'] for c in C])
    KAPPA = int(same([c['minac']['recommended_minac'] for c in C]))
    ZTHR = same([c['sample_qc']['z_threshold'] for c in C])
    fmt_n = lambda x: f'{int(x):,}'                                                        # noqa: E731
    fmt_s = lambda d: f"{d['all']:,} ({d['case']:,}/{d['control']:,})"                    # noqa: E731
    SUP = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')

    def band(s):
        v = [int(x) for x in str(s).split(',')]
        run_ = all(v[i] == v[i - 1] + 1 for i in range(1, len(v)))
        return str(v[0]) if len(v) == 1 else (f'{v[0]}–{v[-1]}' if run_ else ', '.join(map(str, v)))

    def pval(p):
        if p >= 0.01:
            return f'{p:.2f}'
        m, e = f'{p:.0e}'.split('e')
        return f'{m} × 10{str(int(e)).translate(SUP)}'
    gap = lambda g: f"{g['pct']:+.2f}% (P {pval(g['p'])})".replace('-', '−')               # noqa: E731
    NAVY, MID, PALE, RULE, GREY_ = '#1F3864', '#BDD7EE', '#DDEBF7', '#808080', '#404040'
    CC = [(colors[c], tints[c]) for c in cohorts]
    COH = [(c, labels[c].title()) for c in cohorts]
    ws.hide_gridlines(2)
    base = {'font_name': 'Arial', 'font_size': 12, 'valign': 'vcenter', 'align': 'center', 'text_wrap': True,
            'border': 1, 'border_color': RULE}
    F = lambda **k: wb.add_format({**base, **k})                                            # noqa: E731
    title_ = wb.add_format({'font_name': 'Arial', 'font_size': 16, 'bold': True, 'font_color': NAVY})
    subtitle = wb.add_format({'font_name': 'Arial', 'font_size': 13, 'bold': True, 'font_color': NAVY})
    hdr = F(bold=True, font_color='#FFFFFF', bg_color=NAVY)
    group = F(bold=True, font_color='#FFFFFF', bg_color=NAVY)
    label_ = F(bold=True, bg_color=PALE)
    label_mid = F(bold=True, bg_color=MID)
    body = F()
    body_l = F(align='left')
    bold = F(bold=True)
    c_hdr = [F(bold=True, font_color='#FFFFFF', bg_color=c) for c, _ in CC]
    c_sub = [F(bold=True, font_color=c, bg_color=t) for c, t in CC]
    c_out = [F(bold=True, font_color=c, bg_color=t) for c, t in CC]
    note = wb.add_format({'font_name': 'DengXian', 'font_size': 10, 'font_color': GREY_, 'valign': 'vcenter'})
    for col, w in zip('ABCDEFGHIJK', [1.5, 15, 38, 16, 42, 22, 13.5, 22, 13.5, 22, 13.5]):
        ws.set_column(f'{col}:{col}', w)
    ws.set_row(0, 24)
    ws.write(0, 1, 'Rare-Variant QC in v6', title_)
    r = 2
    ws.set_row(r, 20)
    ws.set_row(r + 1, 20)
    for col, t in [(1, 'Step'), (2, 'Exclude'), (3, 'Threshold'), (4, 'Threshold set by')]:
        ws.merge_range(r, col, r + 1, col, t, hdr)
    for i, (_, name) in enumerate(COH):
        c0 = 5 + 2 * i
        ws.merge_range(r, c0, r, c0 + 1, name, c_hdr[i])
        ws.write(r + 1, c0, 'Samples (case/ctrl)', c_sub[i])
        ws.write(r + 1, c0 + 1, 'Variants', c_sub[i])
    r += 2
    kappa1 = [c['sweep']['1'] for c in C]
    for c, k1 in zip(C, kappa1):
        vq = c['variant_qc']
        assert vq['n_in'] - vq['n_removed'] - (vq['n_mono_after_sample_qc'] - vq['n_mono_removed_by_delta']) == k1
    steps = [
        ('Input', 'Fixed-model samples; variants with control MAF < 0.01', None, None,
         [fmt_s(c['n_in']) for c in C], [fmt_n(c['variant_qc']['n_in']) for c in C], False),
        ('Sample QC', 'Samples whose burden rateᵃ departs from that expected at their depth', f'|robust Zᵇ| > {ZTHR:g}',
         'Theil–Sen line of the burden rate on mean depth [1]; robust Z of the residual',
         [fmt_s(c['n_kept']) for c in C],
         [fmt_n(c['variant_qc']['n_in'] - c['variant_qc']['n_mono_after_sample_qc']) + 'ᶜ' for c in C], False),
        ('Variant QC', 'Variants missed far more often in one depth group than in the other', f'Δᵈ > {KNEE:g}',
         'Knee of the Δ distribution (Kneedle [2]), the same in every cohort',
         [fmt_s(c['n_kept']) for c in C], [fmt_n(v) for v in kappa1], False),
        ('Allele count', 'Variants with a minor allele count below κ (singletons)', f'κ = {KAPPA}',
         'Smallest κ with no significant case–control burden gap (table below)',
         [fmt_s(c['n_kept']) for c in C], [fmt_n(c['minac']['variant_count_at_rec']) for c in C], True),
    ]
    for name, excl, thr, how, samples, variants, final in steps:
        ws.set_row(r, 44)
        ws.write(r, 1, name, group)
        if thr is None:
            ws.merge_range(r, 2, r, 4, excl, body)
        else:
            ws.write(r, 2, excl, body_l)
            ws.write(r, 3, thr, bold)
            ws.write(r, 4, how, body_l)
        for i in range(len(C)):
            fmt = c_out[i] if final else body
            ws.write(r, 5 + 2 * i, samples[i], fmt)
            ws.write(r, 6 + 2 * i, variants[i], fmt)
        r += 1
    ws.set_row(r, 18)
    ws.merge_range(r, 1, r, 10, 'Bold, last row: the samples and variants passed to the gene-based tests.', note)
    r += 2
    ws.set_row(r, 22)
    ws.write(r, 1, 'Basis of the Minor Allele Count Threshold κ', subtitle)
    r += 1
    ws.set_row(r, 20)
    ws.merge_range(r, 1, r, 4, 'Case–control burden gapᵉ', hdr)
    for i, (_, name) in enumerate(COH):
        ws.merge_range(r, 5 + 2 * i, r, 6 + 2 * i, name, c_hdr[i])
    r += 1
    chk = lambda c, k, m: c['model_check'][str(k)][m]                                      # noqa: E731
    rows = [
        ('Gap at κ = 1 (singletons kept)', [gap(chk(c, 1, 'decision')) for c in C], False, label_),
        (f'Gap at κ = {KAPPA}', [gap(chk(c, KAPPA, 'decision')) for c in C], False, label_),
        ('κ with no significant gap (P ≥ 0.05)', [band(c['minac']['nonsig_band']) for c in C], False, label_),
        ('κ chosen: the smallest of them', [c['minac']['recommended_minac'] for c in C], True, label_),
        (f'Gap at κ = {KAPPA}, association covariates', [gap(chk(c, KAPPA, 'association')) for c in C], False,
         label_mid),
    ]
    g0 = r
    for lab, vals, is_bold, lf in rows:
        ws.set_row(r, 22)
        ws.merge_range(r, 2, r, 4, lab, lf)
        for i, v in enumerate(vals):
            ws.merge_range(r, 5 + 2 * i, r, 6 + 2 * i, v, c_out[i] if is_bold else body)
        r += 1
    ws.merge_range(g0, 1, g0 + 3, 1, 'Decision model', group)
    ws.write(g0 + 4, 1, 'Check', group)
    r += 1
    foot = [
        'a: Burden rate r = 1,000 × minor alleles / called genotypes of a sample; mean depth, the sample\'s observed '
        'mean coverage.',
        'b: Robust Z = 0.6745 × (residual − median) / MAD; MAD, median absolute deviation.',
        'c: Variants left monomorphic by the sample removal are removed.',
        'd: Δ = |VMISS(30x) − VMISS(15x)|; VMISS, share of samples with no genotype at a variant, within the 15x or the '
        '30x samples (target depth).',
        'e: Case–control difference in r, % of the mean (OLS, P in brackets). Decision model: r ~ case + z(mean '
        f'depth) [3]; association covariates: r ~ case + sex + {n_pcs} BBJ-mainland PCs, as in the gene-based tests.',
        '',
        'References',
        '[1] Sen, P. K. (1968). Estimates of the regression coefficient based on Kendall\'s tau. Journal of the '
        'American Statistical Association, 63(324), 1379–1389.',
        '[2] Satopaa, V., Albrecht, J., Irwin, D., & Raghavan, B. (2011). Finding a "Kneedle" in a haystack: '
        'Detecting knee points in system behavior. In 2011 31st International Conference on Distributed Computing '
        'Systems Workshops (pp. 166–171). IEEE.',
        '[3] Garner, C. (2011). Confounded by sequencing depth in association studies of rare alleles. Genetic '
        'Epidemiology, 35(4), 261–268.',
    ]
    for t in foot:
        ws.set_row(r, 15)
        if t:
            ws.write(r, 1, t, note)
        r += 1
    ws.print_area(0, 0, r - 1, 10)
    ws.set_landscape()
    ws.set_paper(9)
    ws.fit_to_pages(1, 1)
    ws.set_margins(left=0.2, right=0.2, top=0.2, bottom=0.2)


def verify(path, D, a, stats, sample_ids):
    """Re-open the workbook: sheets, rows against the tables, header letters against footnotes, text scan (no
    missing-value artefacts, internal names or sample IDs), and the XML Excel is strict about."""
    from openpyxl import load_workbook
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)          # openpyxl drops the data-bar extension on reading
        wb = load_workbook(path, rich_text=True)
    want = SHEETS if stats['ot_on'] else [s for s in SHEETS if s != SH_OT]
    assert wb.sheetnames == want, wb.sheetnames
    genes = D['genes'].sort_values('gene_order')
    cohorts = csv_list(a.cohorts)
    ws = wb[SH_SIG]
    got = [plain(ws.cell(6 + k * len(cohorts) * 2, 1).value) for k in range(len(genes))]   # 2 variant sets per gene
    assert got == genes.gene_symbol.tolist(), got
    ws = wb[SH_LOVO]
    for r, gene in stats['lovo']['head_rows']:
        assert plain(ws.cell(r, 1).value).startswith(gene + ' '), (r, gene, plain(ws.cell(r, 1).value)[:40])
    assert not any(d.hidden for d in ws.row_dimensions.values()), 'no hidden rows'
    ws = wb[SH_VAR]
    n_v = stats['variants']['rows']
    labels_ = [plain(ws.cell(r, 1).value) for r in range(6, 6 + n_v)]
    assert [x for x in labels_ if x] == genes.gene_symbol.tolist(), 'one gene label per gene, in order'
    if stats['ot_on']:
        wo = wb[SH_OT]
        sig = [str(ws.cell(r, 4 + 0).value) for r in range(6, 6 + n_v)]       # rsID column of the variants sheet
        ids = [wo.cell(r + 1, 3).value for r, _v in stats['variant_rows']]
        vids = [v for _r, v in stats['variant_rows']]
        assert ids == vids and len(ids) == n_v, 'Open Targets variant rows differ from the variants sheet'
        rec = [r for r in range(OT_H + 2, wo.max_row + 1) if wo.cell(r, 5).value in OT_SOURCES]
        assert len(rec) == stats['n_det'], (len(rec), stats['n_det'])
        assert all(wo.row_dimensions[r].outlineLevel == 2 and wo.row_dimensions[r].hidden for r in rec)
        del sig
    for sheet in wb.worksheets:
        if sheet.title != SH_QC:
            head = set()
            for r in (3, 4, 5):
                for c in range(1, sheet.max_column + 1):
                    head |= set(plain(sheet.cell(r, c).value, sup=True))
            notes, notes_text, has_dash = set(), '', False
            for r in range(6, sheet.max_row + 1):
                v = sheet.cell(r, 1).value
                if v is None:
                    continue
                sup, txt = plain(v, sup=True), plain(v)
                if sup and txt.startswith(' ') and len(sup) == 1:
                    notes.add(sup)
                notes_text += txt
            for row in sheet.iter_rows(min_row=6):
                for c in row:
                    if c.value is not None and plain(c.value) == MISSING:
                        has_dash = True
            head -= {'†', '‡', ''}
            assert head <= notes, (sheet.title, sorted(head), sorted(notes))
            if has_dash:
                assert '—,' in notes_text or '—:' in notes_text, (sheet.title, 'a — without its reason')
        for row in sheet.iter_rows():
            for c in row:
                if c.value is None:
                    continue
                s = plain(c.value) + plain(c.value, sup=True)
                assert s not in ('nan', 'None') and 'nan' not in s.split() and 'NaN' not in s, \
                    (sheet.title, c.coordinate, s)
                assert '…' not in s, (sheet.title, c.coordinate, s[:80])
                allowed = SAMPLE_SHEET_TOKENS if sheet.title == SH_MS else ()
                bad = [t for t in FORBIDDEN_TOKENS if t in s and t not in allowed]
                assert not bad, (sheet.title, c.coordinate, bad, s[:80])
                if sheet.title != SH_MS:                  # sample IDs belong to sheet Multi-site carriers only
                    assert not (set(re.findall(r'[A-Za-z0-9_-]+', s)) & sample_ids), (sheet.title, c.coordinate,
                                                                                      'sample ID')
    z = zipfile.ZipFile(path)
    for name in z.namelist():
        if name.startswith('xl/worksheets/sheet'):
            x = z.read(name).decode()
            assert 'inlineStr' not in x, f'inline strings in {name}'
            # a frozen pane needs a split (Excel offers to repair the file otherwise)
            for pane in re.findall(r'<pane [^>]*>', x):
                assert 'xSplit=' in pane or 'ySplit=' in pane, (name, pane)
            # merged ranges must not overlap (Excel repairs those too)
            cells = set()
            for ref in re.findall(r'<mergeCell ref="([A-Z]+\d+:[A-Z]+\d+)"', x):
                from openpyxl.utils import range_boundaries
                c0, r0_, c1, r1 = range_boundaries(ref)
                block = {(rr_, cc_) for rr_ in range(r0_, r1 + 1) for cc_ in range(c0, c1 + 1)}
                assert not (block & cells), (name, 'overlapping merged ranges', ref)
                cells |= block
    ss = z.read('xl/sharedStrings.xml').decode()
    bare = [t for p, t in re.findall(r'<t( xml:space="preserve")?>([^<]*)</t>', ss) if not p and t != t.strip()]
    assert not bare, bare[:3]


def read_multisite(path):
    """gene_signals.multisite.tsv: IDs, depths, codes and lists as text (numeric IDs keep their leading zeros)."""
    m = read_tsv(path)
    for c in ('observed_depth', 'age_at_dna', 'age_at_diagnosis', 'max_mpap', 'gene_order', 'n_variants_carried'):
        m[c] = pd.to_numeric(m[c].replace('', np.nan))
    return m


def read_lovo(path):
    """gene_signals.lovo.tsv: flags as booleans, the left-out variant as text ('' for the test on all variants)."""
    L = read_table(path)
    for c in ('primary', 'top'):
        L[c] = L[c].astype(str) == 'True'
    L['dropped'] = L.dropped.fillna('').astype(str)
    L['k'] = L.k.astype(int)
    return L


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--genes-long', required=True)
    ap.add_argument('--variants-long', required=True)
    ap.add_argument('--scans', required=True)
    ap.add_argument('--genes', required=True)
    ap.add_argument('--checks', required=True)
    ap.add_argument('--upstream-manifest', required=True)
    ap.add_argument('--opentargets-detail', required=True)
    ap.add_argument('--opentargets-meta', required=True)
    ap.add_argument('--rv-qc', nargs='+', required=True, help='rv_qc.<cohort>.json')
    ap.add_argument('--variant-terms', required=True)
    ap.add_argument('--variant-flow', required=True)
    ap.add_argument('--multisite', required=True, help='gene_signals.multisite.tsv (by sample ID)')
    ap.add_argument('--samples', nargs='+', required=True, help='samples.<cohort>.tsv: IDs the workbook must not name')
    ap.add_argument('--lovo', required=True, help='gene_signals.lovo.tsv (leave-one-variant-out)')
    ap.add_argument('--cohorts', required=True)
    ap.add_argument('--cohort-codes', required=True)
    ap.add_argument('--cohort-labels', required=True)
    ap.add_argument('--cohort-colors', required=True)
    ap.add_argument('--cohort-tints', required=True)
    ap.add_argument('--strata', required=True)
    ap.add_argument('--af-cohort', required=True)
    ap.add_argument('--af-bars', default='true')
    ap.add_argument('--rare-maf', type=float, required=True)
    ap.add_argument('--sens-min-num-var', type=int, required=True)
    ap.add_argument('--table-label', default='')
    ap.add_argument('--tommo-label', required=True)
    ap.add_argument('--kg-label', required=True)
    ap.add_argument('--kg-n-eas', type=int, required=True)
    ap.add_argument('--kg-n-jpt', type=int, required=True)
    ap.add_argument('--gnomad-label', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--out-stats', required=True)
    a = ap.parse_args()
    a.af_bars = str(a.af_bars).lower() in ('true', '1', 'yes')
    for opt in ('cohort_codes', 'cohort_labels', 'cohort_colors', 'cohort_tints'):
        assert len(csv_list(getattr(a, opt))) == len(csv_list(a.cohorts)), opt
    checks = load_json(a.checks)
    assert checks['all_passed']
    VL = read_table(a.variants_long)
    for c in ('rsID', 'rsID_unresolved', 'Gene_EnsDb', 'tx_gene', 'transcript', 'transcript_biotype', 'hgvs_c',
              'hgvs_p', 'other_consequences', 'tommo_filter', 'ot_clinvar', 'ot_gwas_top', 'ot_l2g', 'ot_qtl_top',
              'ot_e2g_top', 'callset_reason'):
        if c in VL:
            VL[c] = VL[c].fillna('').astype(str)
    for c in [c for c in VL.columns if c.startswith('in_set_')] + ['in_callset']:
        VL[c] = VL[c].astype(str) == 'True'
    for c in ('ot_gwas_n', 'ot_gwas_focus_n', 'ot_qtl_n', 'ot_qtl_trans_n', 'ot_e2g_n', 'ot_clinvar_n'):
        if c in VL:
            VL[c] = VL[c].fillna(0)
    VL['chrom'] = VL.chrom.astype(str)
    G = read_table(a.genes_long)
    for c in ('tested', 'skato_bh', 'skato_bonferroni', 'cmc_bh', 'cmc_bonferroni', 'in_sens', 'skato_bh_sens',
              'skato_bonferroni_sens', 'cmc_bh_sens', 'cmc_bonferroni_sens'):
        G[c] = G[c].astype(str) == 'True'
    otd = read_table(a.opentargets_detail)
    for c in ('category', 'type', 'trait_or_target', 'biosample', 'study_id', 'project', 'first_author', 'method',
              'l2g_gene', 'interval_type', 'clinical_significance', 'review_status', 'datasource', 'p_source',
              'beta_source', 'lead_variant', 'lead_rsid', 'target_chrom', 'qtl_class', 'cohorts', 'ld_population'):
        otd[c] = otd[c].fillna('').astype(str)
    otd['lead_is_variant'] = otd.lead_is_variant.astype(str).isin(['True', 'true', '1', '1.0'])
    otd['ot_is_trans'] = otd.ot_is_trans.map(lambda x: '' if pd.isna(x) else str(x).lower())
    otd['year'] = otd.year.map(lambda y: '' if pd.isna(y) else str(int(y)) if isinstance(y, float) else str(y))
    otd['focus'] = otd.focus.astype(str).isin(['True', 'true', '1', '1.0']) if len(otd) else pd.Series(dtype=bool)
    qc = {}
    for f in a.rv_qc:
        j = load_json(f)
        qc[j['cohort']] = j
    flow = read_table(a.variant_flow)
    flow['stratum'] = flow.stratum.fillna('').astype(str)
    D = {'G': G, 'VL': VL, 'SC': read_table(a.scans), 'genes': read_table(a.genes), 'up': load_json(a.upstream_manifest),
         'otd': otd, 'otm': load_json(a.opentargets_meta), 'qc': qc, 'terms': read_table(a.variant_terms),
         'flow': flow, 'ms': read_multisite(a.multisite), 'lovo': read_lovo(a.lovo)}
    stats = build(a.out, D, a)
    sample_ids = set().union(*(set(read_tsv(f).IID) for f in a.samples))      # as text: IDs keep leading zeros
    verify(a.out, D, a, stats, sample_ids)
    out = {k: v for k, v in stats.items() if k not in ('variant_rows', 'gene_rows')}
    out['lovo'] = {k: v for k, v in out['lovo'].items() if k != 'head_rows'}
    json.dump(out, open(a.out_stats, 'w'), indent=1, default=str)
    log(f'workbook: {len(SHEETS)} sheets; ' + '; '.join(f'{k} {v}' for k, v in out.items()) + f' -> {a.out}')


if __name__ == '__main__':
    main()
