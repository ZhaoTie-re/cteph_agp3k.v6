#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The audit of the written workbook. Re-opens hla_records.xlsx (openpyxl,
#           read only) and checks, sheet by sheet:
#             cells     every body cell shows what the writer meant to write
#                       (workbook_cells.json.gz), at the row and column of its block;
#                       on the detail sheets every number and P also equals the
#                       table it came from, row by row through the row keys
#             P         every displayed P has the decks' format
#             notes     every header letter has its note; every '—' (missing) and
#                       '–' (not applicable) on a sheet has its line under 'Missing
#                       values'
#             text      no nan / None / True / False / '…'; no internal names on the
#                       content sheets; no sample identifier on any sheet
#             XML       the parts Excel is strict about (inline strings, frozen panes,
#                       overlapping merges, untrimmed shared strings)
# Inputs  : --xlsx --meta (workbook.json) --cells (workbook_cells.json.gz) --typing
#           --assoc (the table folders) --sample-ids (a table with column sample_id)
# Outputs : --out (audit.json); exits non-zero if any check fails
# Component: hla_records
# ---------------------------------------------------------------------------
import argparse
import gzip
import json
import re
import sys
import warnings
import zipfile
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook
from openpyxl.cell.rich_text import CellRichText
from openpyxl.utils import range_boundaries

sys.path.insert(0, str(Path(__file__).resolve().parent))
import deckfmt as F  # noqa: E402

SUP = str.maketrans('−0123456789-', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹⁻')
FORBIDDEN = ['_mainland', 'dvt_yes', 'dvt_no', 'dvt_case_only', 'case_set', 'assoc_hla', 'hla.typing', 'hla_records',
             'typing_tables', 'assoc_tables', 'p_gec', 'p_bonf', 'NOT_TESTED', 'SEPARATION', 'errcode', 'v_other',
             'y_other', 'deck_key', 'test_id', 'maxt_boot', 'effective_tests', 'omnibus_', 'p.value', 'pheno_covar',
             'JHRP', 'bbj_mainland', 'PC1_AVG', '.tsv', '.json', '.py', '.nf', '/LARGE', '/home', 'nextflow', 'results/',
             'run17', 'run16', 'hla_given', 'lead_given']
P_RE = re.compile(r'^(0\.\d+|1\.0|[1-9]\.\d+ × 10⁻[⁰¹²³⁴⁵⁶⁷⁸⁹]+|0\.0\d+)$')
DETAIL = {'T7 Alleles (detail)': ('01', 't7_alleles.tsv', lambda d: d.marker),
          'T8 Residues (detail)': ('01', 't8_residues.tsv', lambda d: d.gene + ':' + d.position + ':' + d.residue),
          'A11 Single markers (detail)': ('02', 'a11_single_markers.tsv', lambda d: d.comparison + '|' + d.cohort + '|' + d.marker),
          'A12 Positions (detail)': ('02', 'a12_positions.tsv', lambda d: d.comparison + '|' + d.cohort + '|' + d.unit)}


def text(v):
    """A cell as shown: rich runs joined, superscript runs as Unicode superscripts."""
    if isinstance(v, CellRichText):
        out = ''
        for b in v:
            if isinstance(b, str):
                out += b
            elif b.font is not None and b.font.vertAlign == 'superscript':
                out += b.text.translate(SUP)
            else:
                out += b.text
        return out
    return '' if v is None else str(v)


def sups(v):
    if isinstance(v, CellRichText):
        return ''.join(b.text for b in v if not isinstance(b, str) and b.font is not None and b.font.vertAlign == 'superscript')
    return ''


class Audit:
    def __init__(self):
        self.checks = []

    def add(self, cid, what, n, bad):
        self.checks.append({'id': cid, 'what': what, 'compared': n, 'failures': len(bad), 'passed': not bad,
                            'detail': bad[:30]})
        print(f'[audit] {cid}: {n} checked, {len(bad)} failures')
        for b in bad[:5]:
            print('   ', b)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for k in ('xlsx', 'meta', 'cells', 'typing', 'assoc', 'sample_ids', 'out'):
        ap.add_argument('--' + k.replace('_', '-'), required=True)
    a = ap.parse_args()
    meta = json.loads(Path(a.meta).read_text())
    with gzip.open(a.cells, 'rt') as f:
        cells = json.load(f)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        wb = load_workbook(a.xlsx, rich_text=True)
    A = Audit()
    names = [s['name'] for s in meta['sheets']]
    A.add('sheets', 'sheet names and order', 1, [] if wb.sheetnames == names else [f'{wb.sheetnames} != {names}'])

    # cells as meant, and the P format
    n, bad, n_p, bad_p = 0, [], 0, []
    for sname, blocks in cells.items():
        ws = wb[sname]
        for b in blocks:
            for i, row in enumerate(b['cells']):
                for j, want in enumerate(row):
                    c = ws.cell(b['first_row'] + i + 1, b['pcols'][j] + 1)
                    n += 1
                    if isinstance(want, (int, float)) and not isinstance(want, bool):
                        ok = isinstance(c.value, (int, float)) and abs(c.value - want) <= 1e-12 * max(1, abs(want))
                    else:
                        ok = text(c.value) == want
                    if not ok:
                        bad.append(f'{sname}!{c.coordinate}: shows {text(c.value)[:40]!r}, meant {str(want)[:40]!r}')
            kinds = next(x['kinds'] for s in meta['sheets'] if s['name'] == sname for x in s['blocks']
                         if x['first_row'] == b['first_row'])
            for j, k in enumerate(kinds):
                if k != 'p':
                    continue
                for i in range(len(b['cells'])):
                    t = text(ws.cell(b['first_row'] + i + 1, b['pcols'][j] + 1).value)
                    if t in ('—', '–'):
                        continue
                    n_p += 1
                    if not P_RE.match(t):
                        bad_p.append(f'{sname} row {b["first_row"] + i + 1}: {t!r}')
    A.add('cells', 'every body cell shows the value the writer meant, at its row and column', n, bad)
    A.add('p-format', 'every P printed in the decks\' format (two significant digits from 0.001, m.mm × 10⁻ⁿ below)',
          n_p, bad_p)

    # detail sheets: every number and P against its source table, through the row keys
    n, bad = 0, []
    for sname, (part, tbl, keyf) in DETAIL.items():
        folder = a.typing if part == '01' else a.assoc
        d = pd.read_csv(Path(folder) / tbl, sep='\t', dtype={'marker': str, 'position': str, 'residue': str, 'unit': str})
        d.index = keyf(d)
        for b in cells[sname]:
            for col, vals in b['raw'].items():
                if col not in d.columns:
                    continue
                src = d.loc[b['row_keys'], col].tolist()
                for k, v, s in zip(b['row_keys'], vals, src):
                    n += 1
                    if (pd.isna(s) and v is None) or v == '\x00na':     # not applicable: the notes give the reason
                        continue
                    try:
                        ok = abs(float(v) - float(s)) <= 1e-9 * max(1.0, abs(float(s)))
                    except (TypeError, ValueError):
                        ok = False
                    if not ok:
                        bad.append(f'{sname} {k} {col}: workbook {v}, table {s}')
            n += 1
            if len(b['row_keys']) != len(d) or set(b['row_keys']) != set(d.index):
                bad.append(f'{sname}: rows {len(b["row_keys"])} vs table {len(d)}')
    A.add('detail', 'detail sheets: every number and P equals its table (row by row)', n, bad)

    # notes: header letters, and the reasons for '—' and '–'
    bad, n = [], 0
    for s in meta['sheets']:
        ws = wb[s['name']]
        head, body_rows = set(), set()
        for b in s['blocks']:
            for r in range(b['header_row'] + 1, b['header_row'] + 4):
                for c in range(1, ws.max_column + 1):
                    head |= set(sups(ws.cell(r, c).value))
            body_rows |= set(range(b['first_row'] + 1, b['first_row'] + b['n_rows'] + 1))
        notes, notes_text, dash, endash = set(), '', False, False
        for r in range(s['blocks'][-1]['first_row'] + s['blocks'][-1]['n_rows'] + 1 if s['blocks'] else 1, ws.max_row + 1):
            v = ws.cell(r, 1).value
            if v is None:
                continue
            sp, tx = sups(v), text(v)
            if sp and len(sp) == 1 and tx.startswith(sp.translate(SUP)):
                notes.add(sp)
            notes_text += tx + '\n'
        for r in body_rows:
            for c in range(1, ws.max_column + 1):
                t = text(ws.cell(r, c).value)
                dash |= t == '—'
                endash |= t == '–'
        n += 1
        missing = sorted(head - notes)
        if missing:
            bad.append(f"{s['name']}: header letters without a note: {missing}")
        if dash and '—,' not in notes_text and '—:' not in notes_text:
            bad.append(f"{s['name']}: a — without its reason")
        if endash and '–,' not in notes_text and '–:' not in notes_text and s['name'] != 'Contents':
            bad.append(f"{s['name']}: a – without its reason")
    A.add('notes', 'every header letter has its note; every — and – has its reason under Missing values', n, bad)

    # text: artefacts, internal names (content sheets), sample identifiers (all sheets)
    ids = set(pd.read_csv(a.sample_ids, sep='\t', dtype=str, usecols=['sample_id']).sample_id)
    n, bad = 0, []
    for ws in wb.worksheets:
        content = ws.title != 'Sources and checks'
        for row in ws.iter_rows():
            for c in row:
                if c.value is None:
                    continue
                s = text(c.value)
                n += 1
                words = set(re.findall(r'[A-Za-z0-9_.-]+', s))
                if words & {'nan', 'NaN', 'None', 'True', 'False', 'inf'} or '…' in s or '\x00' in s:
                    bad.append(f'{ws.title}!{c.coordinate}: artefact in {s[:60]!r}')
                if content:
                    hit = [t for t in FORBIDDEN if t in s]
                    if hit:
                        bad.append(f'{ws.title}!{c.coordinate}: internal name {hit} in {s[:60]!r}')
                if words & ids:
                    bad.append(f'{ws.title}!{c.coordinate}: sample identifier')
    A.add('text', 'no artefacts; no internal names on the content sheets; no sample identifier on any sheet '
          f'({len(ids):,} identifiers)', n, bad)

    # XML Excel is strict about
    bad, n = [], 0
    z = zipfile.ZipFile(a.xlsx)
    for name in z.namelist():
        if not name.startswith('xl/worksheets/sheet'):
            continue
        n += 1
        x = z.read(name).decode()
        if 'inlineStr' in x:
            bad.append(f'{name}: inline strings')
        for pane in re.findall(r'<pane [^>]*>', x):
            if 'xSplit=' not in pane and 'ySplit=' not in pane:
                bad.append(f'{name}: frozen pane without a split')
        seen = set()
        for ref in re.findall(r'<mergeCell ref="([A-Z]+\d+:[A-Z]+\d+)"', x):
            c0, r0, c1, r1 = range_boundaries(ref)
            blk = {(r, c) for r in range(r0, r1 + 1) for c in range(c0, c1 + 1)}
            if blk & seen:
                bad.append(f'{name}: overlapping merged range {ref}')
            seen |= blk
    ss = z.read('xl/sharedStrings.xml').decode()
    bare = [t for p, t in re.findall(r'<t( xml:space="preserve")?>([^<]*)</t>', ss) if not p and t != t.strip()]
    bad += [f'untrimmed shared string {t[:30]!r}' for t in bare[:5]]
    A.add('xml', 'no inline strings, no split-less frozen pane, no overlapping merges, trimmed shared strings', n, bad)

    import hashlib
    out = {'all_passed': all(c['passed'] for c in A.checks), 'xlsx_md5': hashlib.md5(Path(a.xlsx).read_bytes()).hexdigest(),
           'checks': A.checks}
    Path(a.out).write_text(json.dumps(out, indent=1, ensure_ascii=False) + '\n')
    if not out['all_passed']:
        raise SystemExit('audit: FAILED')
    print('[audit] all checks passed')


if __name__ == '__main__':
    main()
