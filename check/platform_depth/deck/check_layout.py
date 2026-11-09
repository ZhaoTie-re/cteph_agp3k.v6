#!/usr/bin/env python3
"""Layout checks for the deck; both must report 0 before a preview is saved.

    python3 check_layout.py hla_typing.pptx hla_typing.preview.pdf

1. Rendered text (the PDF): two text lines drawn on top of each other.
2. Shapes (the .pptx): one element reaching into another, or past the slide's lower edge. A table is measured at
   the height its wrapped text gives it, and a text box at the height of its text, because both grow beyond their
   declared size. A label inside its own rectangle is the callout idiom and is not reported. The estimate is
   deliberately generous, so a reported overlap of a few hundredths of an inch should be looked at in the render
   before it is acted on."""
import math
import re
import subprocess
import sys
import zipfile
from xml.etree import ElementTree as ET

A = 'http://schemas.openxmlformats.org/drawingml/2006/main'
P = 'http://schemas.openxmlformats.org/presentationml/2006/main'
EMU = 914400.0
SLIDE_H = 7.5

def text_of(el):
    return ' '.join(t.text or '' for t in el.iter(f'{{{A}}}t')).strip()

def table_height(gf):
    tbl = gf.find(f'.//{{{A}}}tbl')
    cols = [int(c.get('w')) / EMU for c in tbl.iter(f'{{{A}}}gridCol')]
    total = 0.0
    for tr in tbl.iter(f'{{{A}}}tr'):
        declared = int(tr.get('h')) / EMU
        need = declared
        ci = 0
        for tc in tr.findall(f'{{{A}}}tc'):
            if tc.get('hMerge') == '1':          # covered by the gridSpan cell before it
                continue
            if tc.get('vMerge') == '1':          # covered by the rowSpan cell above it
                ci += int(tc.get('gridSpan', '1'))
                continue
            span = int(tc.get('gridSpan', '1'))
            w = sum(cols[ci:ci + span]); ci += span
            pr = tc.find(f'{{{A}}}tcPr')
            ml = int(pr.get('marL', 91440)) / EMU if pr is not None else 0.1
            mr = int(pr.get('marR', 91440)) / EMU if pr is not None else 0.1
            mt = int(pr.get('marT', 45720)) / EMU if pr is not None else 0.05
            mb = int(pr.get('marB', 45720)) / EMU if pr is not None else 0.05
            sz = [int(r.get('sz')) / 100 for r in tc.iter(f'{{{A}}}rPr') if r.get('sz')]
            pt = max(sz) if sz else 12
            txt = text_of(tc)
            if not txt:                          # an empty cell keeps its row's stored height (spacer rows)
                continue
            per_line = max(1, int((w - ml - mr) / (0.6 * pt / 72)))
            n = max(1, math.ceil(len(txt) / per_line)) if txt else 1
            # a cell merged over several rows shares its height among them
            need = max(need, (n * 1.25 * pt / 72 + mt + mb) / int(tc.get('rowSpan', '1')))
        total += need
    return total

def boxes(xml):
    root = ET.fromstring(xml)
    tree = root.find(f'.//{{{P}}}cSld/{{{P}}}spTree')
    out = []
    for el in list(tree):
        tag = el.tag.split('}')[1]
        if tag not in ('sp', 'pic', 'graphicFrame'):
            continue
        xf = el.find(f'.//{{{P}}}xfrm') if tag == 'graphicFrame' else el.find(f'.//{{{A}}}xfrm')
        if xf is None:
            continue
        off = xf.find(f'{{{A}}}off'); ext = xf.find(f'{{{A}}}ext')
        x, y = int(off.get('x')) / EMU, int(off.get('y')) / EMU
        w, h = int(ext.get('cx')) / EMU, int(ext.get('cy')) / EMU
        nv = el.find(f'.//{{{P}}}cNvPr')
        name = nv.get('name', '') if nv is not None else ''
        # a copy of a .doc.rev6 table: made in PowerPoint (S, Q, P decks), it is drawn at its stored size
        if name.startswith('SRC:') and name.split(':')[1] in 'SQP':
            out.append(('table' if tag == 'graphicFrame' else tag, x, y, w, h, text_of(el)[:40], name))
            continue
        if tag == 'graphicFrame' and el.find(f'.//{{{A}}}tbl') is not None:
            h = max(h, table_height(el)); tag = 'table'
        elif tag == 'sp' and el.find(f'.//{{{A}}}t') is not None:
            body = el.find(f'.//{{{P}}}txBody')
            anchor = body.find(f'{{{A}}}bodyPr').get('anchor', 't') if body is not None else 't'
            need = 0.1
            for par in el.iter(f'{{{A}}}p'):
                sz = [int(r.get('sz')) / 100 for r in par.iter(f'{{{A}}}rPr') if r.get('sz')]
                pt = max(sz) if sz else 18
                txt = ' '.join(t.text or '' for t in par.iter(f'{{{A}}}t'))
                per_line = max(1, int((w - 0.1) / (0.55 * pt / 72)))
                need += max(1, math.ceil(len(txt) / per_line)) * 1.2 * pt / 72
            if need > h and anchor == 't':
                h = need
        out.append((tag, x, y, w, h, text_of(el)[:40], name))
    return out

def contains(a, b, e=.02):
    return a[1] - e <= b[1] and a[2] - e <= b[2] and a[1] + a[3] + e >= b[1] + b[3] and a[2] + a[4] + e >= b[2] + b[4]



def check_shapes(pptx, tol=0.05):
    z = zipfile.ZipFile(pptx)
    names = sorted((n for n in z.namelist() if re.match(r'ppt/slides/slide\d+\.xml$', n)),
                   key=lambda n: int(re.search(r'(\d+)', n).group(1)))
    bad = 0
    for i, n in enumerate(names, 1):
        bs = boxes(z.read(n))
        for b in bs:
            if b[5] and b[2] + b[4] > SLIDE_H + 0.1:
                print(f'slide {i:2d}  PAST EDGE by {b[2] + b[4] - SLIDE_H:.2f} in   [{b[0]} "{b[5][:44]}"]'); bad += 1
        for j in range(len(bs)):
            for k in range(j + 1, len(bs)):
                a, b = bs[j], bs[k]
                if contains(a, b) or contains(b, a) or not (a[5] or b[5]):
                    continue
                if a[6].startswith('SRC:') and a[6] == b[6]:   # parts of one copied drawing overlap by design
                    continue
                # a callout label is drawn inside its rectangle; its text estimate can poke out of it
                if a[0] == 'sp' and not a[5] and b[0] == 'sp' and a[1] <= b[1] and a[1] + a[3] >= b[1] + b[3]:
                    continue
                ox = min(a[1] + a[3], b[1] + b[3]) - max(a[1], b[1])
                oy = min(a[2] + a[4], b[2] + b[4]) - max(a[2], b[2])
                if ox > tol and oy > tol:
                    print(f'slide {i:2d}  overlap {oy:.2f} in  [{a[0]} y{a[2]:.2f} "{a[5][:30]}"]'
                          f' <> [{b[0]} y{b[2]:.2f} "{b[5][:30]}"]'); bad += 1
    return bad


def check_render(pdf):
    xml = subprocess.run(['pdftotext', '-bbox-layout', pdf, '-'], capture_output=True, text=True).stdout
    xml = re.sub(r'<!DOCTYPE[^>]*>', '', xml)
    xml = re.sub(r'\sxmlns="[^"]+"', '', xml, count=1)
    root = ET.fromstring(xml)
    bad = 0
    for pno, page in enumerate(root.iter('page'), 1):
        H = float(page.get('height'))
        lines = []
        for ln in page.iter('line'):
            b = [float(ln.get(k)) for k in ('xMin', 'yMin', 'xMax', 'yMax')]
            words = ' '.join(w.text or '' for w in ln.iter('word'))
            lines.append((b, words))
        for b, w in lines:                                      # text past the bottom edge
            if b[3] > H - 1:
                print(f'slide {pno:2d}  OFF-PAGE  "{w[:60]}"'); bad += 1
        for i in range(len(lines)):
            for j in range(i + 1, len(lines)):
                (a, wa), (b, wb) = lines[i], lines[j]
                oy = min(a[3], b[3]) - max(a[1], b[1])
                ox = min(a[2], b[2]) - max(a[0], b[0])
                ha, hb = a[3] - a[1], b[3] - b[1]
                if oy > 0.35 * min(ha, hb) and ox > 2:
                    print(f'slide {pno:2d}  COLLIDE  "{wa[:45]}"  <>  "{wb[:45]}"'); bad += 1
    return bad


if __name__ == '__main__':
    n1 = check_render(sys.argv[2])
    n2 = check_shapes(sys.argv[1])
    print(f'rendered text: {n1} collision(s); shapes: {n2} problem(s)')
    sys.exit(1 if n1 or n2 else 0)
