#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The two decks as they stand: per slide, its number, section and title
#           (with the italic runs of gene symbols), read from the .pptx files,
#           and the deck's slide registry (the ORDER array of its build script)
#           to tie each slide to its key. The workbook's Contents sheet maps every
#           sheet to these slides, so the map follows the decks when they change.
# Inputs  : --typing-pptx --typing-js --assoc-pptx --assoc-js
# Outputs : --out (decks.json)
# Component: hla_records
# ---------------------------------------------------------------------------
import argparse
import hashlib
import json
import re

from lxml import etree
from pptx import Presentation

NS = {'p': 'http://schemas.openxmlformats.org/presentationml/2006/main',
      'p14': 'http://schemas.microsoft.com/office/powerpoint/2010/main',
      'a': 'http://schemas.openxmlformats.org/drawingml/2006/main'}


def order_keys(js):
    m = re.search(r'const ORDER = \[(.*?)\];', open(js).read(), flags=re.S)
    if not m:
        raise SystemExit(f'{js}: no ORDER array')
    return re.findall(r'"([^"]+)"', m.group(1))


def title_runs(slide):
    """The slide's title: the placeholder title if present, else the one text box whose runs are 24 pt underlined
    (the decks' title style). Runs as [text, italic]."""
    shapes = [s for s in slide.shapes if s.has_text_frame]
    for s in shapes:
        if s.is_placeholder and s.placeholder_format.type is not None and 'TITLE' in str(s.placeholder_format.type):
            return [[r.text, bool(r.font.italic)] for p in s.text_frame.paragraphs for r in p.runs]
    for s in shapes:
        runs = [r for p in s.text_frame.paragraphs for r in p.runs]
        if runs and all(r.font.size is not None and r.font.size.pt == 24 for r in runs if r.text.strip()):
            return [[r.text, bool(r.font.italic)] for r in runs]
    return []


def sections(prs):
    root = prs.part._element
    out = {}
    ids = [s.slide_id for s in prs.slides]
    for sec in root.iter('{%s}section' % NS['p14']):
        for sid in sec.iter('{%s}sldId' % NS['p14']):
            out[int(sid.get('id'))] = sec.get('name')
    return [out.get(i, '') for i in ids]


def deck(pptx, js, name):
    prs = Presentation(pptx)
    keys = order_keys(js)
    secs = sections(prs)
    slides = []
    for n, (s, sec) in enumerate(zip(prs.slides, secs), 1):
        runs = title_runs(s)
        slides.append({'number': n, 'key': keys[n - 1] if n <= len(keys) else None, 'section': sec,
                       'title': ''.join(t for t, _ in runs).strip(), 'runs': runs})
    if len(slides) != len(keys):
        raise SystemExit(f'{name}: {len(slides)} slides but {len(keys)} keys in ORDER')
    md5 = hashlib.md5(open(pptx, 'rb').read()).hexdigest()
    return {'name': name, 'slides': slides, 'n': len(slides), 'md5': md5}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--typing-pptx', required=True)
    ap.add_argument('--typing-js', required=True)
    ap.add_argument('--assoc-pptx', required=True)
    ap.add_argument('--assoc-js', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    out = {'typing': deck(a.typing_pptx, a.typing_js, 'HLA typing'),
           'association': deck(a.assoc_pptx, a.assoc_js, 'HLA association')}
    for k, d in out.items():
        missing = [s['number'] for s in d['slides'] if not s['title'] and s['number'] != 1]
        if missing:
            raise SystemExit(f'{k}: slides without a title: {missing}')
    json.dump(out, open(a.out, 'w'), indent=1, ensure_ascii=False)
    print(f"[deck_index] typing {out['typing']['n']} slides, association {out['association']['n']} slides")


if __name__ == '__main__':
    main()
