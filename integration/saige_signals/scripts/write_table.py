#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The deliverable workbook, written with XlsxWriter (Excel-conformant
#           output: shared strings, rich strings, print titles, page breaks).
#           Presentation only — every number is taken as is from
#           saige_signals.long.tsv and the Open Targets detail table.
#             Sheet SHEET_MAIN      the academic table: booktabs rules, three
#                                   header levels with footnote letters, cohort
#                                   blocks in the ancestry-refinement figure's
#                                   colours, italic gene symbols, P as m x 10^-e,
#                                   footnotes as short lines below the table
#             Sheet "Open Targets"  every Open Targets record behind the summary
#                                   columns, one row per record
#           The workbook names no file, script, parameter or pipeline (asserted).
# Inputs  : saige_signals.long.tsv, loci.tsv, checks.json, the upstream manifest,
#           opentargets_detail.tsv, opentargets_meta.json
# Outputs : saige_signals.xlsx
# ---------------------------------------------------------------------------
import argparse
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
from sig_common import MISSING, as_bool, csv_list, load_json, log, read_table, read_tsv, sig4  # noqa: E402

FORBIDDEN_TOKENS = ['cohort_compare', 'panel (b)', 'panel (c)', 'lead_crosscohort', 'locus_rep', 'assoc_saige',
                    'random_model', 'bbj_mainland', '_mainland', 'p.value', 'no_spa', 'JHRP', 'assoc_gene',
                    '.tsv', '.py', '/LARGE', '/home', 'susie_finemap', 'variant_annot', 'nextflow']
FONT = 'Arial'
INK, INK2, GREY, FAINT, RULE_GREY, NOTE = '#1A1A1A', '#595959', '#7F7F7F', '#A6A6A6', '#BFBFBF', '#404040'
IMPACT_STYLE = {'HIGH': ('#9C0006', True), 'MODERATE': ('#9C5700', True), 'LOW': (INK, False), 'MODIFIER': (GREY, False)}
NUM_AF = '[=0]0;[<0.001]"<0.001";0.000'
NUM_PIP = '[=0]0.000;[<0.001]"<0.001";0.000'
THIN, MEDIUM = 1, 2
ROLE_TEXT = {'Lead': 'Lead', 'Alt lead': 'Alt. lead', 'CS member': 'CS variant'}
SHEET_MAIN = 'Common-variant signals'
SHEET_OT = 'Common-variant Open Targets'
QTL_COL = 'cis / trans (to target TSS)'                 # sheet SHEET_OT, molQTL records
POP_COL = 'Population (% of sample)'                    # sheet SHEET_OT, GWAS / molQTL records
OT_H = 2                                                  # header row of sheet SHEET_OT (0-based)
GNOMAD_POP = {'AFR': 'African / African American', 'AMI': 'Amish', 'AMR': 'admixed American',
              'ASJ': 'Ashkenazi Jewish', 'CSA': 'Central / South Asian', 'EAS': 'East Asian', 'FIN': 'Finnish',
              'MID': 'Middle Eastern', 'NFE': 'non-Finnish European', 'SAS': 'South Asian'}
OT_SOURCES = ('ClinVar', 'GWAS', 'molQTL', 'E2G')       # record order under a variant (sheet SHEET_OT)
VARIANT_FILL = '#F2F2F2'                                  # variant rows of sheet Open Targets
LOCUS_FILL = '#E4E4E4'                                    # locus header rows of sheet Open Targets


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
    head = '10' if abs(m - 1) < 1e-9 else f'{m:g} × 10'
    return [head, ('sup', f'−{-e}')]


def so(term):
    return str(term).replace('&', ' & ').replace('_', ' ')


def text_len(parts):
    return sum(len(p if isinstance(p, str) else p[1]) for p in parts)


def source_text(p_src, b_src):
    """Where a credible-set record's P and beta come from: the variant itself or the set's lead."""
    t = {'variant': 'this variant', 'lead': 'set lead'}
    p, b = t.get(str(p_src), '—'), t.get(str(b_src), '—')
    return p if p == b else f'P: {p}; β: {b}'


def plain(v, sup=False):
    """Text of a cell read back with rich_text=True, without its superscript runs (sup=True: only those)."""
    from openpyxl.cell.rich_text import CellRichText
    if isinstance(v, CellRichText):
        return ''.join((b if isinstance(b, str) else b.text) for b in v
                       if (not isinstance(b, str) and b.font is not None and b.font.vertAlign == 'superscript') == sup)
    return '' if (v is None or sup) else str(v)


def dist_text(bp):
    return f'{bp:,.0f} bp' if bp < 1e3 else f'{bp / 1e3:,.0f} kb' if bp < 1e6 else f'{bp / 1e6:.1f} Mb'


def qtl_text(x):
    """cis / trans of a molQTL record with its evidence (distance to the target gene's TSS, or the gene's
    chromosome); Open Targets' set-level flag added where it contradicts the class."""
    if x.qtl_class in ('', 'unknown'):
        t = 'unknown (no gene location)'
    elif pd.isna(x.qtl_distance):
        t = f'{x.qtl_class} (gene on {x.target_chrom})'
    else:
        t = f'{x.qtl_class} ({dist_text(x.qtl_distance)})'
    if (x.qtl_class, x.ot_is_trans) in (('cis', 'true'), ('trans', 'false')):
        t += f' · Open Targets: {"trans" if x.ot_is_trans == "true" else "not trans"}'
    return t


def wrap_parts(parts, width):
    """Split rich parts (str / (style, str)) into lines of <= width characters, at spaces."""
    tokens = []                                       # (style, word + following space)
    for p in parts:
        st, tx = ('', p) if isinstance(p, str) else p
        words = tx.split(' ')
        for k, wd in enumerate(words):
            tokens.append((st, wd + (' ' if k < len(words) - 1 else '')))
    lines, cur, n = [], [], 0
    for st, tx in tokens:
        if tx == '':
            continue
        # a superscript sticks to the word before it; a line breaks only before a new word
        if n + len(tx.rstrip()) > width and cur and st != 'sup' and cur[-1][1].endswith(' '):
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
        """A font-only format for one rich-text fragment."""
        props = {'font_name': FONT, 'font_size': sz, 'font_color': color, 'bold': bool(b), 'italic': bool(i)}
        if style == 'i':
            props['italic'] = True
        elif style == 'b':
            props['bold'] = True
        elif style == 'sup':
            props['font_script'] = 1
        elif style.startswith('#'):
            props['font_color'], props['bold'] = style, True
        key = ('font',) + tuple(sorted(props.items()))
        if key not in self._fmt:
            self._fmt[key] = self.wb.add_format(props)
        return self._fmt[key]

    def rich(self, r, c, parts, cell, b=False, i=False, color=INK, sz=9):
        """parts: str (plain) or (style, str), style in {'i', 'b', 'sup', '#RRGGBB'}; cell: format dict."""
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
                     bold=True if (b or st == 'b' or st.startswith('#')) else None,
                     italic=True if (i or st == 'i') else None, font_script=1 if st == 'sup' else None)
            self.ws.write_string(r, c, tx, self.fmt(**f))
        else:
            args = []
            for st, tx in frags:
                args += [self.font(st, b, i, color, sz), tx]
            self.ws.write_rich_string(r, c, *args, self.fmt(**cell))


def build(path, long, loci, up, otd, otm, a):
    cohorts = csv_list(a.cohorts)
    codes = dict(zip(cohorts, csv_list(a.cohort_codes)))
    labels = dict(zip(cohorts, csv_list(a.cohort_labels)))
    colors = {c: '#' + x.upper() for c, x in zip(cohorts, csv_list(a.cohort_colors))}
    tints = {c: '#' + x.upper() for c, x in zip(cohorts, csv_list(a.cohort_tints))}
    code_color = {codes[c]: colors[c] for c in cohorts}
    ot_on = bool(otm.get('queried'))

    # columns: (block, sub-spanner, leaf, footnote letter, width, frequency key)
    spec = [('var', None, 'Locus', 'a', 12, ''), ('var', None, 'Class', 'a', 6, ''),
            ('var', None, 'Role', 'a', 9, ''), ('var', None, 'rsID', 'h', 24, ''),
            ('var', None, 'Chr', None, 4.5, ''), ('var', None, 'Position', None, 11, ''),
            ('var', None, 'EA/OA', None, 12, ''), ('var', None, 'Gene', 'h', 18, ''),
            ('var', None, 'CS in', 'a', 7, ''),
            ('fun', None, 'Consequence', None, 30, ''), ('fun', None, 'Impact', None, 9, ''),
            ('fun', None, 'Transcript', None, 38, ''), ('fun', None, 'HGVS', None, 27, ''),
            ('fun', None, 'Other consequences', None, 46, '')]
    if ot_on:
        spec += [('ot', None, 'ClinVar', None, 30, ''), ('ot', None, 'GWAS credible sets', None, 72, ''),
                 ('ot', None, 'L2G', None, 13, ''), ('ot', None, 'molQTL credible sets', None, 64, ''),
                 ('ot', None, 'Enhancer-to-gene', None, 50, '')]
    # one width for every frequency column, so the bars of a row are comparable; wide enough for the widest
    # header ('gnomAD NFE', bold) on one line
    AFW = 12.5
    spec += [('af', 'This study', 'CTEPH', None, AFW, 'af_cteph'), ('af', 'This study', 'AGP3K', None, AFW, 'af_agp3k'),
             ('af', 'Japanese', 'ToMMo', None, AFW, 'tommo_af'), ('af', 'Japanese', '1KG JPT', None, AFW, 'kg_jpt_af'),
             ('af', 'East Asian', '1KG EAS', None, AFW, 'kg_eas_af'),
             ('af', 'East Asian', 'gnomAD EAS', None, AFW, 'gnomad_eas_af'),
             ('af', 'European', 'gnomAD NFE', None, AFW, 'gnomad_nfe_af')]
    cohort_cols = [(None, 'Peak', 'c', 5.5), (None, 'OR (95% CI)', 'c', 16), (None, 'P', 'c', 10),
                   (None, 'PIP', 'd', 6.5), (None, 'CS', 'e', 34),
                   ('Genotypes', 'CTEPH', 'f', 14), ('Genotypes', 'AGP3K', 'f', 16),
                   ('EAF', 'CTEPH', None, 6.5), ('EAF', 'AGP3K', None, 6.5),
                   ('Missing', 'CTEPH', 'f', 7), ('Missing', 'AGP3K', 'f', 7),
                   ('HWE P', 'CTEPH', 'f', 8), ('HWE P', 'AGP3K', 'f', 8)]
    for c in cohorts:
        spec += [(c, sub, leaf, note, w, '') for sub, leaf, note, w in cohort_cols]
    ncol = len(spec)
    blocks = []
    for j, s in enumerate(spec):
        if not blocks or blocks[-1][0] != s[0]:
            blocks.append([s[0], j, j])
        else:
            blocks[-1][2] = j
    block_first = {b[1]: b[0] for b in blocks}
    j_ot = next((b[1] for b in blocks if b[0] == 'ot'), None)
    j_af = [j for j, s in enumerate(spec) if s[0] == 'af']
    FROZEN = 7                                            # columns A:G stay on screen

    wb = xlsxwriter.Workbook(path, {'nan_inf_to_errors': False})
    label = f'{a.table_label}. ' if a.table_label.strip() else ''     # optional, e.g. 'Table 2'
    wb.set_properties({'title': f'{label}CTEPH association signals and credible sets',
                       'author': 'CTEPH GWAS project'})
    sh = Sheet(wb, wb.add_worksheet(SHEET_MAIN))
    ws = sh.ws
    ws.hide_gridlines(2)
    ws.set_zoom(90)
    ws.set_default_row(12.75)
    for j, s in enumerate(spec):
        ws.set_column(j, j, s[4])

    ROW_T, ROW_D, ROW_S, ROW_SS, ROW_H, ROW0 = 0, 1, 2, 3, 4, 5
    wide = long.drop_duplicates('variant_id').set_index('variant_id')
    wide = wide.sort_values(['order', 'role_rank', 'max_pip', 'pos'], ascending=[True, True, False, True])
    by_c = {c: long[long.cohort == c].set_index('variant_id') for c in cohorts}
    n_case = {c: int(by_c[c].n_case.iloc[0]) for c in cohorts}
    n_ctrl = {c: int(by_c[c].n_ctrl.iloc[0]) for c in cohorts}
    n_rows = len(wide)
    last = ROW0 + n_rows - 1
    n_lead = int((wide.role != 'CS member').sum())
    cls_n = loci.tier.value_counts().to_dict()

    # ── title ──────────────────────────────────────────────────────────────
    sh.rich(ROW_T, 0, [('b', f'{label}CTEPH association signals across three nested ancestry-filtered '
                             'cohorts and their 95% credible sets')], {'align': 'left'}, sz=11)
    sh.rich(ROW_D, 0, [f'{len(loci)} loci ({cls_n.get("GW", 0)} GW, {cls_n.get("SG", 0)} SG) · {n_lead} lead and '
                       f'{n_rows - n_lead} credible-set variants · GRCh38 · credible-set variants are grouped under '
                       'their lead (outline, left)'], {'align': 'left'}, i=True, color=INK2)
    ws.set_row(ROW_T, 18)

    # ── header ─────────────────────────────────────────────────────────────
    span_text = {'var': ('Variant', 'a'), 'fun': ('Functional annotation', 'h'),
                 'ot': ('Open Targets', 'i'), 'af': ('Effect-allele frequency', 'g')}
    first_of = {b[0]: b[1] for b in blocks}

    def hcell(blk, j, row, **extra):
        d = {'align': 'center', 'valign': 'bottom', 'text_wrap': True, 'bold': True}
        if blk in cohorts:
            d.update({'bg_color': tints[blk], 'pattern': 1})
            if j == first_of[blk]:
                d.update({'left': MEDIUM, 'left_color': colors[blk]})
        if row == ROW_S:
            coh = blk in cohorts
            d.update({'top': MEDIUM, 'top_color': INK, 'valign': 'vcenter',
                      'bottom': MEDIUM if coh else THIN, 'bottom_color': colors[blk] if coh else INK})
        if row == ROW_H:
            d.update({'bottom': MEDIUM, 'bottom_color': INK})
        d.update(extra)
        return d

    for blk, j0, j1 in blocks:
        if blk in cohorts:
            parts = [(colors[blk], labels[blk]), ('sup', 'b'),
                     f'   {n_case[blk]:,} cases / {n_ctrl[blk]:,} controls']
            color = colors[blk]
        else:
            t, note = span_text[blk]
            parts, color = [('b', t), ('sup', note)], INK
        for jj in range(j0, j1 + 1):                         # every cell of the spanner carries its rules
            ws.write_blank(ROW_S, jj, None, sh.fmt(**hcell(blk, jj, ROW_S)))
        if j1 > j0:
            ws.merge_range(ROW_S, j0, ROW_S, j1, '', sh.fmt(**hcell(blk, j0, ROW_S)))
        sh.rich(ROW_S, j0, parts, hcell(blk, j0, ROW_S), b=True, color=color, sz=10)
    j = 0
    while j < ncol:
        blk, sub = spec[j][0], spec[j][1]
        if sub is None:
            leaf, note = spec[j][2], spec[j][3]
            f = hcell(blk, j, ROW_H)
            ws.merge_range(ROW_SS, j, ROW_H, j, '', sh.fmt(**f))
            sh.rich(ROW_SS, j, [('b', leaf), ('sup', note or '')], f, b=True)
            j += 1
            continue
        k = j
        while k + 1 < ncol and spec[k + 1][0] == blk and spec[k + 1][1] == sub:
            k += 1
        note = next((s[3] for s in spec[j:k + 1] if s[3]), '')
        f = hcell(blk, j, ROW_SS, bottom=THIN, bottom_color=INK)
        for jj in range(j, k + 1):
            ws.write_blank(ROW_SS, jj, None, sh.fmt(**hcell(blk, jj, ROW_SS, bottom=THIN, bottom_color=INK)))
        if k > j:
            ws.merge_range(ROW_SS, j, ROW_SS, k, '', sh.fmt(**f))
        sh.rich(ROW_SS, j, [('b', sub), ('sup', note)], f, b=True)
        for jj in range(j, k + 1):
            sh.rich(ROW_H, jj, [('b', spec[jj][2])], hcell(blk, jj, ROW_H), b=True)
        j = k + 1
    ws.set_row(ROW_S, 20)
    ws.set_row(ROW_H, 14)

    # ── body ───────────────────────────────────────────────────────────────
    def cs_text(x, role):
        st = x.cs_status
        muted = {'italic': True, 'font_color': GREY}
        if st == 'in_cs':
            return f'CS{int(x.cs_index)} ({int(x.cs_size)})', {}
        if st == 'below_cut':                         # the PIP is in the column to the left
            if int(x.n_cs) == 1:
                return f'Not in CS1 (PIP < {sig4(x.cs_cut)})', {'font_color': INK2}
            cuts = [(k, sig4(float(v))) for k, v in (t.split(':') for t in x.cs_cuts_all.split(';'))]
            return (f'Not in {", ".join(f"CS{k}" for k, _v in cuts)} (PIP < {", ".join(v for _k, v in cuts)})',
                    {'font_color': INK2})
        if st == 'no_cs_null':
            return 'No CS: no signal', muted
        if st == 'no_cs_impure':                      # lead rows: the (strongest effect's) dropped 95% set
            if role == 'CS member':
                return 'No CS: low LD', muted
            return (f'No CS: low LD ({int(x.dropped_set_size):,} variants; min |r| {x.dropped_set_min_r:.2f})',
                    muted)
        if st == 'not_tested':
            return 'Not tested', muted
        return 'Not fine-mapped', muted

    locus_seen, r = set(), ROW0
    for vid, w in wide.iterrows():
        role = w.role
        bold, ital = role == 'Lead', role == 'Alt lead'
        first = w.locus_id not in locus_seen
        locus_seen.add(w.locus_id)
        written = set()

        def rules(j, d):
            """Row rules: grey rule above a locus, bottom rule after the last row, cohort-coloured left rule."""
            blk = block_first.get(j)
            if blk in cohorts:
                d.update({'left': MEDIUM, 'left_color': colors[blk]})
            if first and r > ROW0:
                d.update({'top': THIN, 'top_color': RULE_GREY})
            if r == last:
                d.update({'bottom': MEDIUM, 'bottom_color': INK})
            return d

        def put(j, value, num=None, align='center', **kw):
            written.add(j)
            if value is None or (isinstance(value, float) and not np.isfinite(value)):
                ws.write_string(r, j, MISSING, sh.fmt(**rules(j, {'align': 'center', 'font_color': FAINT})))
                return
            d = rules(j, {'align': align, 'bold': bold or None, 'italic': ital or None, 'num_format': num})
            d.update(kw)
            if isinstance(value, str):
                ws.write_string(r, j, value, sh.fmt(**d))
            else:
                ws.write_number(r, j, float(value), sh.fmt(**d))

        def rich(j, parts, align='left', b=None, i=None, color=INK):
            written.add(j)
            sh.rich(r, j, parts, rules(j, {'align': align}), b=bold if b is None else b,
                    i=ital if i is None else i, color=color)

        if first:
            rich(0, [('i', w.label)], b=True)
            put(1, w.tier, bold=True)
        put(2, ROLE_TEXT[role], align='left')
        if w.rsID not in ('.', ''):
            put(3, w.rsID, align='left')
        elif w.rsID_unresolved not in ('.', ''):
            rich(3, [w.rsID_unresolved.replace(';', '/'), ('sup', '‡')])
        else:
            put(3, None)
        put(4, int(w.chrom) if str(w.chrom).isdigit() else w.chrom, num='0')
        put(5, int(w.pos), num='#,##0', align='right')
        put(6, f'{w.EA}/{w.OA}')
        rich(7, [('i', w.Gene_EnsDb)])
        csin = w.cs_cohorts if isinstance(w.cs_cohorts, str) else ''
        if csin:
            parts = []
            for k, ch in enumerate(csin):
                parts += ([' '] if k else []) + [(code_color[ch], ch)]
            rich(8, parts, align='center')
        else:
            put(8, None)
        # functional annotation
        eff = w.snpEff_Effect if w.snpEff_Effect not in ('.', '') else None
        put(9, so(eff) if eff else None, align='left')
        imp = w.snpEff_Impact if w.snpEff_Impact not in ('.', '') else None
        if imp:
            c_, b_ = IMPACT_STYLE.get(imp, (INK, False))
            put(10, imp.title(), font_color=c_, bold=(b_ or bold) or None)
        else:
            put(10, None)
        gsym = w.snpEff_Gene if w.snpEff_Gene not in ('.', '') else ''
        fid = w.snpEff_FeatureID if w.snpEff_FeatureID not in ('.', '') else ''
        bt = w.snpEff_Biotype if w.snpEff_Biotype not in ('.', '') else ''
        if eff and 'intergenic' in eff:
            rich(11, [('i', gsym.replace('-', '–', 1) if gsym.count('-') == 1 else gsym), ' (intergenic)'])
        elif gsym or fid:
            tail = f' ({so(bt)})' if bt and bt != 'protein_coding' else ''
            rich(11, [('i', gsym), ' · ' if gsym and fid and gsym != fid else '', fid if fid != gsym else '', tail])
        else:
            put(11, None)
        hc = w.snpEff_HGVS_c if w.snpEff_HGVS_c not in ('.', '') else ''
        hp = w.snpEff_HGVS_p if w.snpEff_HGVS_p not in ('.', '') else ''
        put(12, (hc + (f' ({hp})' if hp else '')) or None, align='left')
        shown = set(str(eff).split('&')) if eff else set()
        allt = [t for x in str(w.snpEff_all_effects).split(',') for t in x.split('&') if t not in ('', '.')]
        others = [so(t) for t in dict.fromkeys(allt) if t not in shown]
        if others:
            put(13, ', '.join(others), align='left', font_color=INK2, bold=None)
        # Open Targets: one line per cell, the number of records then up to OtTopN of them, unwrapped (text
        # beyond the column is clipped on screen and read in full in the formula bar)
        if ot_on:
            jo = j_ot
            if str(w.ot_in_index) not in ('True', 'true', '1', '1.0'):
                for jj in range(jo, jo + 5):
                    put(jj, None)
            else:
                n_f, n_t = int(w.ot_gwas_focus_n), int(w.ot_qtl_trans_n)
                cells = {jo: (int(w.ot_clinvar_n), '', w.ot_clinvar),
                         jo + 1: (int(w.ot_gwas_n), f' ({n_f} thrombosis)' if n_f else '', w.ot_gwas_top),
                         jo + 3: (int(w.ot_qtl_n), f' ({n_t} trans)' if n_t else '', w.ot_qtl_top),
                         jo + 4: (int(w.ot_e2g_n), '', w.ot_e2g_top)}
                for jj, (n, extra, top) in cells.items():
                    put(jj, f'{n}{extra} · {top}' if n else None, align='left', bold=None, italic=None)
                if w.ot_l2g:
                    g_, s_ = w.ot_l2g.split(' ', 1)
                    rich(jo + 2, [('i', g_), ' ' + s_], b=False, i=False)
                else:
                    put(jo + 2, None)
        # frequencies
        for jj in j_af:
            key = spec[jj][5]
            val = w[key]
            kw = {'italic': True, 'font_color': GREY, 'bold': None} \
                if (key == 'tommo_af' and w.tommo_match != 'absent' and w.tommo_filter != 'PASS') else {}
            put(jj, None if pd.isna(val) else float(val), num=NUM_AF, align='right', **kw)
        if a.af_bars:
            ws.conditional_format(r, j_af[0], r, j_af[-1],
                                  {'type': 'data_bar', 'bar_color': '#C9C9C9', 'bar_solid': True,
                                   'bar_no_border': True, 'min_type': 'num', 'min_value': 0, 'max_type': 'max'})
        # cohort blocks
        for blk, j0, _j1 in blocks:
            if blk not in cohorts:
                continue
            x = by_c[blk].loc[vid]
            if role != 'CS member':
                if x.called_peak in ('genome_wide', 'suggestive'):
                    put(j0, 'GW' if x.called_peak == 'genome_wide' else 'SG',
                        bold=True if (x.called_peak == 'genome_wide' or bold) else None)
                else:
                    put(j0, '–', font_color=GREY, bold=None)
            put(j0 + 1, None if pd.isna(x.OR) else f'{x.OR:.2f} ({x.L95:.2f}–{x.U95:.2f})', align='right')
            pp = p_parts(x.P)
            if pp:
                rich(j0 + 2, pp, align='right')
            else:
                put(j0 + 2, None)
            put(j0 + 3, None if pd.isna(x.pip) else float(x.pip), num=NUM_PIP, align='right',
                **({'font_color': GREY} if x.cs_status == 'no_cs_null' else {}))
            txt, kw = cs_text(x, role)
            if 'italic' in kw:
                kw = dict(kw, bold=None)
            put(j0 + 4, txt, align='left', **kw)
            gf = '{}/{}/{}/{}'
            put(j0 + 5, gf.format(int(x.case_homOA), int(x.case_het), int(x.case_homEA), int(x.case_missing_n)))
            put(j0 + 6, gf.format(int(x.ctrl_homOA), int(x.ctrl_het), int(x.ctrl_homEA), int(x.ctrl_missing_n)))
            put(j0 + 7, float(x.case_eaf), num='0.000', align='right')
            put(j0 + 8, float(x.ctrl_eaf), num='0.000', align='right')
            put(j0 + 9, float(x.case_missing_rate), num='0.0%', align='right')
            put(j0 + 10, float(x.ctrl_missing_rate), num='0.0%', align='right')
            for jj, hv in ((j0 + 11, x.case_hwe_p), (j0 + 12, x.ctrl_hwe_p)):
                hp_ = p_parts(hv)
                if hp_:
                    rich(jj, hp_, align='right')
                else:
                    put(jj, None)
        # empty cells still carry the row's rules
        for jj in range(ncol):
            if jj not in written:
                ws.write_blank(r, jj, None, sh.fmt(**rules(jj, {})))
        if role == 'CS member':
            ws.set_row(r, None, None, {'level': 1})
        r += 1
    assert r - 1 == last
    ws.outline_settings(True, False, True, False)        # summary (lead) rows above their groups
    ws.freeze_panes(ROW0, FROZEN)

    # ── footnotes: short lines in column A, inside the frozen columns ───────
    s = up['susie']
    flank = int(up['peak_flank_bp']) // 1000
    p_gw, p_sg = float(up['p_genomewide']), float(up['p_suggestive'])
    approx = {'no_spa': 'normal approximation', 'spa': 'saddlepoint approximation'}[up['p_value_source']]
    grm = 'sparse GRM' if up.get('sparse_grm') else 'full GRM'
    loco = ', LOCO' if up.get('loco') else ''
    af_lab = labels[a.af_cohort].lower()
    notes = [
        ('a', ['Loci: GW, lead P < ', *sci_parts(p_gw), ' in ≥ 1 cohort; SG, lead P < ', *sci_parts(p_sg),
               ' in all three cohorts. Peaks: variants with P < ', *sci_parts(p_sg), f' within ±{flank} kb '
               f'(lead, minimum P); cohort leads within {flank} kb form one locus unless conditionally independent. '
               'Lead (bold), minimum P across cohorts; Alt. lead (italic), another cohort\'s peak lead at the locus; '
               'CS variant, member of a 95% credible set in ≥ 1 cohort (CS in: '
               + '; '.join(f'{codes[c]}, {labels[c].lower()}' for c in cohorts) + ').']),
        ('b', ['Nested ancestry-filtered mainland Japanese sample sets ('
               + ' ⊂ '.join(labels[c].lower() for c in cohorts) + '); related individuals retained; overlapping '
               'samples, hence not independent replication.']),
        ('c', [f'SAIGE logistic mixed model ({grm}{loco}), additive, adjusted for sex and {up["n_pcs"]} PCs; OR per '
               f'effect allele; P, score test ({approx}). Peak: the cohort\'s call where the variant is its peak '
               'lead (GW, SG); –, not the cohort\'s lead.']),
        ('d', [f'SuSiE (susie_rss; L = {s["L"]}; in-sample LD; ±{flank} kb). {s["coverage"] * 100:g}% credible sets '
               f'retained if min |r| ≥ {s["min_abs_corr"]} (smallest absolute pairwise LD in the set). PIP, '
               'posterior inclusion probability.']),
        ('e', ['CS, the variant\'s status in that cohort\'s SuSiE fit of the locus. CSk (n): in 95% credible set k '
               'of n variants. Not in CSk (PIP < y): the cohort has a credible set, but the variant\'s PIP is below '
               'y, the smallest PIP in the set. No CS: no signal: SuSiE found no effect at the locus (all PIP = 0). '
               'No CS: low LD: SuSiE found a signal, but its 95% set failed the purity rule (min |r| ≥ '
               f'{s["min_abs_corr"]}, note d) and was discarded; lead rows give that set\'s size and min |r| (with '
               'several effects, those of the strongest). Not fine-mapped: no peak of that cohort at the locus, or '
               'the variant outside its window. Not tested: absent from that cohort\'s scan.']),
        ('f', ['Genotypes, OA/OA / OA/EA / EA/EA / missing; Missing, proportion of the group; HWE P, exact test '
               'within the group.']),
        ('g', [f'Effect-allele frequency. CTEPH, AGP3K: {af_lab}. ToMMo, {a.tommo_label}; 1KG, {a.kg_label}, '
               f'JPT (n = {a.kg_n_jpt}) and EAS (n = {a.kg_n_eas}); gnomAD, {a.gnomad_label}, EAS and NFE. —, '
               'absent from the reference; grey italic, record failing ToMMo filters.'
               + (' Bars scaled to the row maximum.' if a.af_bars else '')]),
        ('h', ['rsID: ToMMo 60KJPN; where it lists several at a multi-allelic position, the one dbSNP 151 '
               'assigns to this allele, or the only one it does not assign to another allele. Gene: Ensembl 86 '
               'gene overlapping the variant, else the nearest within 1 Mb. Consequence: most severe snpEff '
               'annotation; other consequences over all transcripts.']
         + ([' ', ('sup', '‡'), ' rsIDs of the position that could not be attributed to this allele.']
            if (wide.rsID_unresolved.fillna('.') != '.').any() else [])),
    ]
    if ot_on:
        top = {2: 'two', 3: 'three', 4: 'four', 5: 'five'}.get(otm.get('top_n'), str(otm.get('top_n')))
        cw, td = otm.get('cis_window', 1_000_000) / 1e6, otm.get('trans_distance', 5_000_000) / 1e6
        notes.append(('i', [f'Open Targets Platform (data release {otm["data_version"]}), queried per variant. Each '
                            f'cell: the number of records · up to {top} of them, chosen as follows. ClinVar: '
                            'clinical significance: condition, in alphabetical order. GWAS credible sets: sets '
                            'containing the variant (in brackets, those with a thrombosis-related trait: '
                            'thrombosis, embolism, venous, pulmonary hypertension or artery, coagulation, fibrin, '
                            'haemostasis) · the distinct traits of the sets with the highest posterior probability '
                            '(PP) of the variant, with that PP. L2G: the gene most often ranked first by '
                            'locus-to-gene across those sets (median score). molQTL credible sets: sets (in brackets, '
                            'those in trans) · the records with the highest PP (QTL type with a cis- or trans- prefix, '
                            f'gene, tissue / cell type, PP); cis, the variant within {cw:g} Mb of the gene\'s '
                            'transcription start site (that of its Ensembl canonical transcript); trans, on another '
                            f'chromosome or more than {td:g} Mb from it; no prefix, {cw:g}–{td:g} Mb (neither) or gene '
                            'location unknown. Enhancer-to-gene: '
                            'target genes · those with the highest score (score, best-scoring tissue / cell type, '
                            'number of tissues / cell types predicting the link). Cells are not wrapped: select one to '
                            f'read text beyond the column. Every record: sheet {SHEET_OT}, same variants in the same '
                            'order, records ranked as here. —, no record; across all five columns, variant absent '
                            'from the Open Targets variant index.']))
    abbrev = ['Abbreviations: AGP3K, control cohort; CI, confidence interval; CS, credible set; CTEPH, chronic '
              'thromboembolic pulmonary hypertension; EA/OA, effect/other allele; EAF, effect-allele frequency; EAS, '
              'East Asian; GRM, genetic relationship matrix; GW, genome-wide; HWE, Hardy–Weinberg equilibrium; JPT, '
              'Japanese in Tokyo; L2G, locus-to-gene; LOCO, leave-one-chromosome-out; NFE, non-Finnish European; OR, '
              'odds ratio; PC, principal component; PIP, posterior inclusion probability; QTL, quantitative trait '
              'locus; SG, suggestive.']
    frozen_px = sum(spec[k][4] for k in range(FROZEN)) * 7 + 5
    width = int(frozen_px / 5.8)                          # Arial 9 averages ~5.3 px per character
    rr = last + 2
    note_lines = []
    for letter, parts in notes + [('', abbrev)]:
        for k, ln in enumerate(wrap_parts(parts, width - 3)):
            lead = [('sup', letter), ' '] if (letter and k == 0) else []
            sh.rich(rr, 0, lead + ln, {'align': 'left', 'indent': 1 if k else None}, color=NOTE)
            note_lines.append(text_len(lead + ln))
            rr += 1
    last_row = rr - 1

    # ── print ──────────────────────────────────────────────────────────────
    ws.set_landscape()
    ws.set_paper(8)                                       # A3
    ws.set_print_scale(60)
    ws.repeat_rows(ROW_S, ROW_H)
    ws.repeat_columns(0, FROZEN - 1)
    ws.set_v_pagebreaks([b[1] for b in blocks if b[0] in cohorts])
    ws.print_area(0, 0, last_row, ncol - 1)
    ws.set_margins(left=0.4, right=0.4)
    ws.set_footer(f'&L{label}CTEPH association signals&RPage &P of &N')

    # ── sheet 2: every Open Targets record ──────────────────────────────────
    n_det, variant_rows, locus_rows, ot_note_lines, ot_width = 0, [], [], [], 0
    if ot_on:
        s2 = Sheet(wb, wb.add_worksheet(SHEET_OT))
        w2 = s2.ws
        w2.hide_gridlines(2)
        w2.set_zoom(90)
        w2.set_default_row(12.75)
        # one meaning per column: the molQTL study's tissue / cell type, and the E2G columns grouped
        # under one prefix (best-scoring tissue / cell type, number predicting the link); third element: the
        # column's note letter (notes below the table, as in SHEET_MAIN)
        cols2 = [('Locus', 12, None), ('Role', 9, None), ('Variant', 22, None), ('rsID', 13, None),
                 ('Source', 9, 'a'), ('Type', 9, 'a'), ('Trait, condition or target gene', 46, None),
                 (QTL_COL, 18, 'd'), ('QTL tissue / cell type', 26, 'e'), ('Study', 15, 'b'),
                 ('Publication', 22, 'b'), ('Study N', 9, 'b'), ('Cohorts', 20, 'b'), (POP_COL, 20, 'b'),
                 ('P', 10, 'c'), ('Beta', 8, 'c'), ('P / Beta from', 22, 'c'), ('Set lead variant', 32, 'c'),
                 ('PP', 7, 'e'), ('Fine-mapping (set)', 11, 'e'), ('L2G gene (set)', 11, 'e'),
                 ('L2G score (set)', 9, 'e'), ('E2G score', 8, 'f'), ('E2G best tissue / cell type', 26, 'f'),
                 ('E2G tissues / cell types (n)', 11, 'f'), ('E2G interval', 9, 'f'),
                 ('Clinical significance', 22, None), ('Review status', 30, None)]
        J = {h: k for k, (h, _wd, _n) in enumerate(cols2)}        # record cells are written by column name
        for k, (_h, wd, _n) in enumerate(cols2):
            w2.set_column(k, k, wd)
        s2.rich(0, 0, [('b', f'Open Targets Platform records for the variants in sheet {SHEET_MAIN} '
                              f'(data release {otm["data_version"]})')], {'align': 'left'}, sz=11)
        s2.rich(1, 0, [f'One row per variant of sheet {SHEET_MAIN}, in the same order (Lead bold, Alt. lead italic), '
                       'grouped under a header row per locus, each with its Open Targets records grouped below it. '
                       'Outline buttons: 1, loci; 2, all variants (as opened); 3, every record.'],
                {'align': 'left'}, i=True, color=INK2)
        H = OT_H
        for k, (h, _wd, note) in enumerate(cols2):
            s2.rich(H, k, [('b', h), ('sup', note or '')], {'align': 'center', 'valign': 'bottom', 'text_wrap': True,
                                                           'top': MEDIUM, 'top_color': INK, 'bottom': MEDIUM,
                                                           'bottom_color': INK}, b=True)
        w2.set_row(H, 39)                                 # three lines: the narrow E2G count header
        # records ranked exactly as the Signals summaries pick them (stable sorts on the fetch order):
        # ClinVar alphabetical, GWAS and molQTL by PP (missing last), E2G by score
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
        # Outline (summary rows above): locus header rows level 0, variant rows level 1, records level 2. Each
        # group has its own summary row, so a lead's records fold exactly like a CS variant's.
        rr2, locus_prev = H + 1, None
        for vid, wv in wide.iterrows():
            recs = by_var.get(vid, d.iloc[0:0])
            level = 1
            if wv.locus_id != locus_prev:
                lv = wide[wide.locus_id == wv.locus_id]
                n_role = lv.role.value_counts()
                roles = [f'{n_role[k]} {t}{"s" if t == "CS variant" and n_role[k] != 1 else ""}'
                         for k, t in (('Lead', 'lead'), ('Alt lead', 'alt. lead'), ('CS member', 'CS variant'))
                         if n_role.get(k, 0)]
                n_idx = int(lv.ot_in_index.astype(str).isin(['True', 'true', '1', '1.0']).sum())
                hfmt = dict({'top': THIN, 'top_color': RULE_GREY} if locus_prev is not None else {},
                            bg_color=LOCUS_FILL, pattern=1, align='left')
                s2.rich(rr2, 0, [('i', wv.label)], hfmt, b=True)
                w2.write_blank(rr2, 1, None, s2.fmt(**hfmt))
                s2.rich(rr2, 2, [f'{len(lv)} variant{"s" if len(lv) != 1 else ""} ({", ".join(roles)}) · '
                                 f'{n_idx} in the Open Targets variant index'], hfmt, color=INK2)
                for k in range(3, len(cols2)):            # blank, so the text runs across
                    w2.write_blank(rr2, k, None, s2.fmt(**hfmt))
                w2.set_row(rr2, None, None, {'level': 0})
                locus_rows.append((rr2, wv.locus_id))
                rr2 += 1
            vb, vi = wv.role == 'Lead', wv.role == 'Alt lead'
            vfmt = dict(bg_color=VARIANT_FILL, pattern=1, align='left')
            w2.write_blank(rr2, 0, None, s2.fmt(**vfmt))   # the locus label is on the header row
            s2.rich(rr2, 1, [ROLE_TEXT[wv.role]], vfmt, b=vb, i=vi)
            s2.rich(rr2, 2, [vid], vfmt, b=vb, i=vi)
            s2.rich(rr2, 3, [wv.rsID if wv.rsID not in ('.', '') else ''], vfmt, b=vb, i=vi)
            if str(wv.ot_in_index) not in ('True', 'true', '1', '1.0'):
                s2.rich(rr2, 4, [('i', 'Not in the Open Targets variant index')], vfmt, color=GREY)
            elif not len(recs):
                s2.rich(rr2, 4, [('i', 'In the Open Targets variant index; no record')], vfmt, color=GREY)
            else:
                n = recs.category.value_counts()
                n_tr = int((recs.qtl_class == 'trans').sum())
                s2.rich(rr2, 4, [' · '.join(f'{k} {int(n.get(k, 0))}' + (f' ({n_tr} trans)' if k == 'molQTL' and n_tr
                                                                          else '') for k in OT_SOURCES)], vfmt)
            for k in range(5, len(cols2)):                # blank, so the summary text runs across
                w2.write_blank(rr2, k, None, s2.fmt(**vfmt))
            w2.set_row(rr2, None, None, {'level': level, 'collapsed': bool(len(recs))})
            variant_rows.append((rr2, vid, level))
            locus_prev = wv.locus_id
            rr2 += 1
            for _, x in recs.iterrows():
                base = {}

                def p2(k, val, num=None, align='left', **kw):
                    f = dict(base, align=align, num_format=num, **kw)
                    if val is None or (isinstance(val, float) and not np.isfinite(val)) or val == '':
                        w2.write_blank(rr2, k, None, s2.fmt(**f))
                    elif isinstance(val, str):
                        w2.write_string(rr2, k, val, s2.fmt(**f))
                    else:
                        w2.write_number(rr2, k, float(val), s2.fmt(**f))
                cs_rec = x.category in ('GWAS', 'molQTL')
                p2(J['Locus'], '')
                p2(J['Role'], '')
                p2(J['Variant'], vid, font_color=FAINT)     # the variant again, for filtered views
                p2(J['rsID'], '')
                p2(J['Source'], x.category)
                p2(J['Type'], x.type if x.type != x.category and x.category != 'ClinVar' else '')
                jt = J['Trait, condition or target gene']
                if x.category in ('molQTL', 'E2G'):
                    s2.rich(rr2, jt, [('i', ' '.join(x.trait_or_target.split()))], dict(base, align='left'))
                else:
                    p2(jt, ' '.join(x.trait_or_target.split()), bold=True if x.focus else None)   # as in Signals
                p2(J[QTL_COL], qtl_text(x) if x.category == 'molQTL' else '')
                p2(J['QTL tissue / cell type'], x.biosample if x.category == 'molQTL' else '')
                p2(J['Study'], x.study_id)
                p2(J['Publication'], f'{x.first_author} {x.year}'.strip() if x.first_author else '')
                p2(J['Study N'], x.n_samples if pd.notna(x.n_samples) else None, num='#,##0', align='right')
                # cohorts as Open Targets lists them, else (molQTL) the study's project in grey italic, else —
                if not cs_rec:
                    p2(J['Cohorts'], '')
                elif x.cohorts:
                    p2(J['Cohorts'], x.cohorts)
                elif x.category == 'molQTL' and x.project:
                    p2(J['Cohorts'], x.project, italic=True, font_color=GREY)
                else:
                    p2(J['Cohorts'], '—')
                p2(J[POP_COL], (x.ld_population or '—') if cs_rec else '')
                # P and beta: the variant's own, else the credible set lead's (grey italic), else blank
                lead_p, lead_b = x.p_source == 'lead', x.beta_source == 'lead'
                pp_ = p_parts(x.p_value) if pd.notna(x.p_value) else None
                if pp_:
                    s2.rich(rr2, J['P'], pp_, dict(base, align='right'), i=lead_p, color=GREY if lead_p else INK)
                else:
                    p2(J['P'], None)
                p2(J['Beta'], x.beta if pd.notna(x.beta) else None, num='0.000', align='right',
                   **({'italic': True, 'font_color': GREY} if lead_b else {}))
                p2(J['P / Beta from'], source_text(x.p_source, x.beta_source) if cs_rec else '')
                if cs_rec:
                    p2(J['Set lead variant'], 'this variant' if x.lead_is_variant else
                       (f'{x.lead_variant} ({x.lead_rsid.replace(";", "/")})' if x.lead_rsid else x.lead_variant))
                else:
                    p2(J['Set lead variant'], '')
                p2(J['PP'], x.pp if pd.notna(x.pp) else None, num='0.000', align='right')
                p2(J['Fine-mapping (set)'], x.method)
                if x.l2g_gene:
                    s2.rich(rr2, J['L2G gene (set)'], [('i', x.l2g_gene)], dict(base, align='left'))
                else:
                    p2(J['L2G gene (set)'], '')
                p2(J['L2G score (set)'], x.l2g_score if pd.notna(x.l2g_score) else None, num='0.00', align='right')
                p2(J['E2G score'], x.e2g_score if pd.notna(x.e2g_score) else None, num='0.00', align='right')
                p2(J['E2G best tissue / cell type'], x.biosample if x.category == 'E2G' else '')
                p2(J['E2G tissues / cell types (n)'], x.n_biosamples if pd.notna(x.n_biosamples) else None,
                   num='#,##0', align='right')
                p2(J['E2G interval'], x.interval_type)
                p2(J['Clinical significance'], x.clinical_significance)
                p2(J['Review status'], x.review_status)
                w2.set_row(rr2, None, None, {'level': 2, 'hidden': True})   # collapsed under the variant
                rr2 += 1
                n_det += 1
        assert len(variant_rows) == n_rows and n_det == len(otd), (len(variant_rows), n_rows, n_det, len(otd))
        w2.outline_settings(True, False, True, False)    # variant (summary) rows above their records
        for k in range(len(cols2)):
            w2.write_blank(rr2, k, None, s2.fmt(top=MEDIUM, top_color=INK))
        # notes: lettered short lines in column A below the table, inside the frozen columns (as in SHEET_MAIN)
        ot_notes = [
            ('a', ['Records: ClinVar, ClinVar evidence on the variant; GWAS and molQTL, credible sets of GWAS and '
                   'molecular QTL studies containing the variant (Type, the QTL type); E2G, enhancer-to-gene '
                   'predictions, one row per target gene. Under each variant: ClinVar in alphabetical order, GWAS '
                   f'and molQTL by PP, E2G by score; the first ones are those named in sheet {SHEET_MAIN}. '
                   'Thrombosis-related GWAS traits in bold.']),
            ('b', ['Study, Open Targets study identifier; Publication, first author and year; Study N, the study\'s '
                   'sample size. Cohorts, as listed by Open Targets; grey italic, the study\'s project where it lists '
                   'none (molQTL studies); —, none listed. Population, share of the study sample by gnomAD reference '
                   'population, as Open Targets assigns them for LD; —, not reported (eQTL Catalogue studies).']),
            ('c', ['P and Beta: the variant\'s own in that study when Open Targets reports them, otherwise those of '
                   'the credible set\'s lead variant (grey italic; source in P / Beta from, lead in Set lead '
                   'variant); blank when neither is reported (— in P / Beta from).']),
            ('d', [f'cis / trans: the variant against the molQTL\'s target gene. cis, within {cw:g} Mb of the '
                   'gene\'s TSS (start of its Ensembl canonical transcript, end on the minus strand); trans, on '
                   f'another chromosome or more than {td:g} Mb away; neither, {cw:g}–{td:g} Mb. Open Targets\' own '
                   'flag (a credible set is trans when every variant in it is) is added where it disagrees.']),
            ('e', ['PP, posterior probability of the variant in the credible set. Fine-mapping, L2G gene and L2G '
                   'score belong to the credible set (fine-mapping method; its top locus-to-gene prediction), not to '
                   'the variant. QTL tissue / cell type, that of the molQTL study.']),
            ('f', ['E2G: per target gene, the best score, the tissue / cell type giving it, the number of tissues / '
                   'cell types predicting the link, and the type of the predicted interval.'])]
        pops = {p.split(' ')[0] for s in otd.ld_population if s for p in s.split(', ')}
        abbr = {'E2G': 'enhancer-to-gene', 'L2G': 'locus-to-gene', 'LD': 'linkage disequilibrium',
                'PP': 'posterior probability', 'TSS': 'transcription start site',
                'QTL': 'quantitative trait locus (eQTL, pQTL, sQTL, tuQTL: expression, protein, splicing, '
                       'transcript-usage; sc-, single-cell)',
                **{p: GNOMAD_POP.get(p, 'gnomAD population') for p in pops}}
        ot_abbrev = ['Abbreviations: ' + '; '.join(f'{k}, {v}' for k, v in sorted(abbr.items(), key=lambda kv: kv[0].lower()))
                     + '.']
        ot_width = int((sum(wd for _h, wd, _n in cols2[:4]) * 7 + 5) / 5.8)     # the frozen columns A:D
        rn = rr2 + 2
        for letter, parts in ot_notes + [('', ot_abbrev)]:
            for k, ln in enumerate(wrap_parts(parts, ot_width - 3)):
                lead = [('sup', letter), ' '] if (letter and k == 0) else []
                s2.rich(rn, 0, lead + ln, {'align': 'left', 'indent': 1 if k else None}, color=NOTE)
                ot_note_lines.append(text_len(lead + ln))
                rn += 1
        if n_det:
            w2.autofilter(H, 0, rr2 - 1, len(cols2) - 1)
        w2.freeze_panes(H + 1, 4)
        w2.set_landscape()
        w2.set_paper(8)
        w2.repeat_rows(H, H)
    wb.close()
    return ncol, n_rows, note_lines, width, n_det, variant_rows, locus_rows, ot_note_lines, ot_width


def verify(path, ncol, n_rows, colors, af_bars, ot_on, n_det, variant_rows, locus_rows):
    """Re-open the workbook: structure, colours, text scan (no truncated or wrapped text), the Open Targets
    sheet against Signals (same variants, same order, outline), and the XML Excel is strict about."""
    from openpyxl import load_workbook
    with warnings.catch_warnings():                   # openpyxl drops the data-bar extension on reading
        warnings.simplefilter('ignore', UserWarning)
        wb = load_workbook(path, rich_text=True)
    assert wb.sheetnames == ([SHEET_MAIN, SHEET_OT] if ot_on else [SHEET_MAIN]), wb.sheetnames
    ws = wb[SHEET_MAIN]
    assert ws.max_column == ncol, (ws.max_column, ncol)
    n_rules = sum(len(cf.rules) for cf in ws.conditional_formatting)
    assert n_rules == (n_rows if af_bars else 0), n_rules
    rules = {str(ws.cell(3, j).border.bottom.color.rgb)[-6:].upper() for j in range(1, ncol + 1)
             if ws.cell(3, j).border.bottom is not None and ws.cell(3, j).border.bottom.color is not None}
    assert {c.upper() for c in colors} <= rules, (colors, rules)
    body = range(6, 6 + n_rows)                             # Signals: one line per variant
    assert all((ws.row_dimensions[r].height or 12.75) <= 12.75 for r in body), 'Signals rows: default height only'
    assert not any('\n' in str(ws.cell(r, c).value) for r in body for c in range(1, ncol + 1)), 'line break'
    # every CS cell in one of the formats of footnote e; lead rows give the dropped set of a low-LD fit
    cs_cols = [k for k in range(1, ncol + 1) if plain(ws.cell(4, k).value) == 'CS']
    cs_fmt = re.compile(r'CS\d+ \(\d+\)|Not in CS\d+(, CS\d+)* \(PIP < [\d.]+(, [\d.]+)*\)|No CS: no signal|'
                        r'No CS: low LD( \([\d,]+ variants; min \|r\| \d\.\d\d\))?|Not fine-mapped|Not tested')
    assert len(cs_cols) == len(colors), cs_cols
    bad = [(r, k, ws.cell(r, k).value) for r in body for k in cs_cols
           if not cs_fmt.fullmatch(str(ws.cell(r, k).value))
           or (ws.cell(r, k).value == 'No CS: low LD' and ws.cell(r, 3).value != 'CS variant')]
    assert not bad, bad[:3]
    if ot_on:
        wo = wb[SHEET_OT]
        # the variant rows are the Signals rows, in order (chr:pos:REF:ALT rebuilt from Signals' Chr, Position,
        # EA/OA with EA = ALT)
        sig = [f'chr{ws.cell(r, 5).value}:{int(ws.cell(r, 6).value)}:{str(ws.cell(r, 7).value).split("/")[1]}:'
               f'{str(ws.cell(r, 7).value).split("/")[0]}' for r in body]
        assert [wo.cell(r + 1, 3).value for r, _v, _l in variant_rows] == sig, 'variant rows differ from Signals'
        for r, _v, lvl in variant_rows:
            d = wo.row_dimensions[r + 1]
            assert (d.outlineLevel or 0) == lvl == 1 and not d.hidden, (r, lvl, d.outlineLevel, d.hidden)
        for r, _l in locus_rows:
            d = wo.row_dimensions[r + 1]
            assert not d.outlineLevel and not d.hidden and wo.cell(r + 1, 2).value is None, r
        h1, b1 = OT_H + 1, OT_H + 2                         # 1-based header row and first body row
        rec = [r for r in range(b1, wo.max_row + 1) if wo.cell(r, 5).value in OT_SOURCES]
        assert len(rec) == n_det, (len(rec), n_det)
        assert all(wo.row_dimensions[r].outlineLevel == 2 and wo.row_dimensions[r].hidden for r in rec)
        head = [plain(wo.cell(h1, k).value) for k in range(1, wo.max_column + 1)]
        assert head[0] == 'Locus' and 'Records' not in plain(wo.cell(OT_H, 1).value), 'header not on its row'
        jsrc = head.index('P / Beta from') + 1
        cs = [r for r in rec if wo.cell(r, 5).value in ('GWAS', 'molQTL')]
        assert all(wo.cell(r, jsrc).value for r in cs), 'a credible-set record without its P / Beta source'
        jcls = head.index(QTL_COL) + 1                     # cis / trans on every molQTL record, and only there
        assert all(bool(wo.cell(r, jcls).value) == (wo.cell(r, 5).value == 'molQTL') for r in rec), 'cis / trans'
        for col in ('Cohorts', POP_COL):                   # on every credible-set record (— if none), only there
            jc = head.index(col) + 1
            assert all(bool(wo.cell(r, jc).value) == (wo.cell(r, 5).value in ('GWAS', 'molQTL')) for r in rec), col
        # every group has its own summary row: each run of level >= k starts right under a row of level k - 1
        lvl = {r: wo.row_dimensions[r].outlineLevel or 0 for r in range(b1, wo.max_row + 1)}
        for r in range(b1 + 1, wo.max_row + 1):
            if lvl[r] > lvl[r - 1]:
                assert lvl[r] == lvl[r - 1] + 1, ('outline jumps a level', r, lvl[r - 1], lvl[r])
        # every header letter has its note below the table, and every note a header
        table_end = max(rec + [r + 1 for r, _v, _l in variant_rows]) + 1      # the closing rule
        sup_head = {plain(wo.cell(h1, k).value, sup=True) for k in range(1, wo.max_column + 1)} - {''}
        sup_note = {plain(wo.cell(r, 1).value, sup=True) for r in range(table_end + 1, wo.max_row + 1)} - {''}
        assert sup_head and sup_head == sup_note, (sorted(sup_head), sorted(sup_note))
    for sheet in wb.worksheets:
        for row in sheet.iter_rows():
            for c in row:
                if c.value is None:
                    continue
                s = str(c.value)
                assert s not in ('nan', 'None') and 'nan' not in s.split(), (sheet.title, c.coordinate, s)
                assert '…' not in s, (sheet.title, c.coordinate, s[:80])
                bad = [t for t in FORBIDDEN_TOKENS if t in s]
                assert not bad, (sheet.title, c.coordinate, bad, s[:80])
    z = zipfile.ZipFile(path)
    # the seven frequency columns share one width (comparable bars within a row)
    widths = {}
    for lo, hi, wd in re.findall(r'<col min="(\d+)" max="(\d+)" width="([\d.]+)"',
                                 z.read('xl/worksheets/sheet1.xml').decode()):
        widths.update({k: float(wd) for k in range(int(lo), int(hi) + 1)})
    af = [k for k in range(1, ncol + 1) if ws.cell(5, k).value in ('ToMMo', '1KG JPT', '1KG EAS', 'gnomAD EAS',
                                                                   'gnomAD NFE')]
    af = list(range(af[0] - 2, af[-1] + 1))
    assert len(af) == 7 and len({widths[k] for k in af}) == 1, [(k, widths.get(k)) for k in af]
    wbx = z.read('xl/workbook.xml').decode()
    t = re.search(r'name="_xlnm.Print_Titles" localSheetId="0">([^<]+)<', wbx)
    assert t and t.group(1).split(',')[0].split('!')[1].startswith('$A'), wbx     # columns first, as Excel writes
    for name in z.namelist():
        if name.startswith('xl/worksheets/sheet'):
            x = z.read(name).decode()
            assert 'inlineStr' not in x, f'inline strings in {name}'
            assert all(m == '1048575' for m in re.findall(r'<brk [^>]*max="(\d+)"', x.split('<rowBreaks')[0])), name
    ss = z.read('xl/sharedStrings.xml').decode()           # edge whitespace kept only with xml:space="preserve"
    bare = [t for p, t in re.findall(r'<t( xml:space="preserve")?>([^<]*)</t>', ss) if not p and t != t.strip()]
    assert not bare, bare[:3]
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--long', required=True)
    ap.add_argument('--loci', required=True)
    ap.add_argument('--checks', required=True)
    ap.add_argument('--upstream-manifest', required=True)
    ap.add_argument('--opentargets-detail', required=True)
    ap.add_argument('--opentargets-meta', required=True)
    ap.add_argument('--cohorts', required=True)
    ap.add_argument('--cohort-codes', required=True)
    ap.add_argument('--cohort-labels', required=True)
    ap.add_argument('--cohort-colors', required=True, help='hex without #, one per cohort')
    ap.add_argument('--cohort-tints', required=True)
    ap.add_argument('--af-cohort', required=True)
    ap.add_argument('--table-label', default='Table')
    ap.add_argument('--af-bars', default='true')
    ap.add_argument('--tommo-label', required=True)
    ap.add_argument('--kg-label', required=True)
    ap.add_argument('--kg-n-eas', type=int, required=True)
    ap.add_argument('--kg-n-jpt', type=int, required=True)
    ap.add_argument('--gnomad-label', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    a.af_bars = str(a.af_bars).lower() in ('true', '1', 'yes')
    for opt in ('cohort_codes', 'cohort_labels', 'cohort_colors', 'cohort_tints'):
        assert len(csv_list(getattr(a, opt))) == len(csv_list(a.cohorts)), opt
    long = read_table(a.long)
    for c in ('peak_id', 'fm_reason', 'cs_cuts_all', 'cs_status', 'called_peak', 'rsID', 'rsID_unresolved',
              'Gene_EnsDb', 'snpEff_Effect', 'snpEff_Impact', 'snpEff_Gene', 'snpEff_FeatureID', 'snpEff_Biotype',
              'snpEff_HGVS_c', 'snpEff_HGVS_p', 'snpEff_all_effects', 'tommo_filter', 'tommo_match', 'cs_cohorts',
              'label', 'role', 'tier', 'locus_id', 'ot_in_index', 'ot_clinvar', 'ot_gwas_top', 'ot_l2g',
              'ot_qtl_top', 'ot_e2g_top', 'ot_rsids'):
        long[c] = long[c].fillna('').astype(str)
    for c in ('in_cs', 'fine_mapped'):
        long[c] = as_bool(long[c])
    for c in ('ot_gwas_n', 'ot_gwas_focus_n', 'ot_qtl_n', 'ot_qtl_trans_n', 'ot_e2g_n', 'ot_clinvar_n'):
        long[c] = long[c].fillna(0)
    long['chrom'] = long.chrom.astype(str)
    loci = read_tsv(a.loci)
    checks = load_json(a.checks)
    up = load_json(a.upstream_manifest)
    otm = load_json(a.opentargets_meta)
    otd = read_table(a.opentargets_detail)
    for c in ('category', 'type', 'trait_or_target', 'biosample', 'study_id', 'project', 'first_author', 'method',
              'l2g_gene', 'interval_type', 'clinical_significance', 'review_status', 'datasource', 'p_source',
              'beta_source', 'lead_variant', 'lead_rsid', 'target_chrom', 'qtl_class', 'cohorts', 'ld_population'):
        otd[c] = otd[c].fillna('').astype(str)
    otd['lead_is_variant'] = otd.lead_is_variant.astype(str).isin(['True', 'true', '1', '1.0'])
    otd['ot_is_trans'] = otd.ot_is_trans.map(lambda x: '' if pd.isna(x) else str(x).lower())   # Open Targets' flag
    otd['year'] = otd.year.map(lambda y: '' if pd.isna(y) else str(int(y)) if isinstance(y, float) else str(y))
    otd['focus'] = as_bool(otd.focus.fillna(False)) if len(otd) else pd.Series(dtype=bool)
    assert checks['all_passed'] and checks['p_value_source'] == up['p_value_source']
    ncol, n_rows, note_lines, width, n_det, variant_rows, locus_rows, ot_note_lines, ot_width = \
        build(a.out, long, loci, up, otd, otm, a)
    assert max(note_lines) <= width + 2, (max(note_lines), width)
    assert not ot_note_lines or max(ot_note_lines) <= ot_width + 2, (max(ot_note_lines), ot_width)
    verify(a.out, ncol, n_rows, csv_list(a.cohort_colors), a.af_bars, bool(otm.get('queried')), n_det, variant_rows,
           locus_rows)
    log(f'table: {n_rows} variant rows x {ncol} columns, {len(note_lines)} footnote lines (<= {max(note_lines)} of '
        f'{width} characters); Open Targets sheet: {len(locus_rows)} loci, {len(variant_rows)} variant rows, '
        f'{n_det} records, {len(ot_note_lines)} note lines (<= {max(ot_note_lines, default=0)} of {ot_width}) '
        f'-> {a.out}')


if __name__ == '__main__':
    main()
