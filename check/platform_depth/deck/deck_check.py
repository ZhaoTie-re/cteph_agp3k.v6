#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Every number printed on the summary deck (slide text, footers and speaker notes) must come from
#           stats.json, the file collect_stats.py fills from results/. A number matches when stats.json holds a value
#           that prints as it (any decimals, percent, thousands separators, significant digits, or m × 10^e), or when
#           the number occurs in a string of stats.json (metric texts, variant IDs, references). Single digits are not
#           checked; the notes' "Timing" and "Backup" sections (minutes, slide numbers) are skipped. Also text rules:
#           American spelling; no DiseaseLens; no composite numbers ("142 + 679"); no internal path or the software
#           name SAIGE on a slide, except on the verbatim rules slide (apRules) and the models slide (models), which
#           names each signal type's published model.
# Inputs  : platform_depth_summary.pptx; slide_keys.json; stats.json
# Outputs : a report on stdout; exit 1 if any number is unmatched or any text rule fails
# Usage   : python3 deck_check.py [platform_depth_summary.pptx]
# ---------------------------------------------------------------------------
import json
import math
import os
import re
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
SUP = str.maketrans('⁻⁰¹²³⁴⁵⁶⁷⁸⁹', '-0123456789')
NUM = re.compile(r'\d+(?:\.\d+)? × 10\^-?\d+|(?<![\d.] × )(?<![\d.])10\^-?\d+|\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d+(?:\.\d+)?')
HEADS = {'Timing', 'Say', 'Numbers', 'If asked', 'Backup', 'Source'}
VERBATIM = {'apRules'}
SOFTWARE_OK = {'apRules', 'models'}
# Sample IDs (user 2026-11-07): only the appendix slide of the ADARB2 census lists carriers by sample. Its IDs and its
# per-sample numbers come from the internal census table (git-ignored; never stats.json, which git may publish): the IDs
# are removed before numbers are read, and its numbers must be values of that table. An ID anywhere else fails.
ID_SLIDES = {'apHcIds'}
CENSUS = os.path.join(HERE, 'cram_check', 'hc_census', 'hc_census.tsv')
CENSUS_NUMS = ('signature_reads', 'region_reads', 'lead_gq', 'ga_1747168_gq', 'ga_1747199_gq')


def census():
    """the census table's sample IDs (a regex, longest first) and its per-sample numbers"""
    if not os.path.exists(CENSUS):
        return None, set()
    import csv
    ids, nums = set(), set()
    with open(CENSUS) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            ids.add(r['sample_id'])
            for k in CENSUS_NUMS:
                if r.get(k, '') not in ('', 'nan'):
                    nums.add(str(int(float(r[k]))))
    rx = re.compile(r'(?<![\w])(?:' + '|'.join(re.escape(x) for x in sorted(ids, key=len, reverse=True)) + r')(?![\w])')
    return rx, nums


def norm(s):
    s = s.replace(' ', ' ').replace(' ', ' ').replace(' ', ' ')
    s = re.sub(r'([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)', lambda m: '^' + m.group(1).translate(SUP), s)
    return s.replace('−', '-').replace('–', '-')


def variants(v):
    out = set()
    a = abs(float(v))
    for d in range(0, 7):
        out.add(f'{a:.{d}f}')
    for d in range(0, 4):
        out.add(f'{100 * a:.{d}f}')
    if a >= 1000:
        for d in range(0, 3):
            out.add(f'{a:,.{d}f}')
    if a > 0 and math.isfinite(a):
        for sig in range(1, 5):
            g = f'{a:.{sig}g}'
            if 'e' not in g:
                out.add(g)
        e = math.floor(math.log10(a))
        for ee in (e, e + 1):
            for dp in range(0, 4):
                out.add(f'{a / 10 ** ee:.{dp}f} × 10^{ee}')
    return out


def pool_of(stats):
    pool = set()

    def walk(x):
        if isinstance(x, bool) or x is None:
            return
        if isinstance(x, (int, float)):
            pool.update(variants(x))
        elif isinstance(x, str):
            for t in re.findall(r'\d+(?:\.\d+)?(?:e-?\d+)?', x):
                pool.add(t)
                try:
                    pool.update(variants(float(t)))
                except ValueError:
                    pass
            for t in re.findall(r'\d{1,3}(?:,\d{3})+', x):
                pool.add(t)
        elif isinstance(x, dict):
            for k, v in x.items():
                walk(k)
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
    walk(stats)
    return pool


def runs_text(xml):
    """Paragraph texts of one part; superscript runs get a '^' prefix."""
    out = []
    for p in re.findall(r'<a:p>(.*?)</a:p>', xml, flags=re.S):
        buf = []
        for attrs, txt in re.findall(r'<a:r>\s*<a:rPr([^>]*)/?>(?:.*?</a:rPr>)?\s*<a:t>(.*?)</a:t>\s*</a:r>', p, flags=re.S):
            t = txt.replace('&amp;', '&').replace('&lt;', '<').replace('&gt;', '>').replace('&apos;', "'").replace('&quot;', '"')
            buf.append(('^' + t) if 'baseline="30000"' in attrs else t)
        if buf:
            out.append(''.join(buf))
    return out


def deck_parts(path):
    z = zipfile.ZipFile(path)
    pres = z.read('ppt/presentation.xml').decode('utf8')
    rels = z.read('ppt/_rels/presentation.xml.rels').decode('utf8')
    rid = {}
    for m in re.finditer(r'<Relationship ([^>]+)/>', rels):
        a = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
        rid[a['Id']] = a['Target']
    slides = []
    for r in re.findall(r'<p:sldId [^>]*r:id="(rId\d+)"', pres):
        name = 'ppt/' + rid[r].split('ppt/')[-1].lstrip('/')
        xml = z.read(name).decode('utf8')
        srels = z.read(name.replace('slides/', 'slides/_rels/') + '.rels').decode('utf8')
        nm = re.search(r'Target="\.\./notesSlides/(notesSlide\d+\.xml)"', srels)
        notes = z.read(f'ppt/notesSlides/{nm.group(1)}').decode('utf8') if nm else ''
        slides.append((runs_text(xml), runs_text(notes)))
    return slides


def tokens(text):
    t = norm(text)
    out = set()
    for tok in NUM.findall(t):
        if re.fullmatch(r'10\^-?\d+', tok):
            tok = '1 × ' + tok
        out.add(tok)
    return out


BRITISH = re.compile(r'\b(analys(?:e|ed|ing)|colours?|labell(?:ed|ing)|haemat\w*|haemo\w*|behaviours?|centres?|'
                     r'favour\w*|modell(?:ed|ing)|programmes?|tumours?|signall(?:ed|ing)|catalogue|grey|'
                     r'organis(?:e|ed|ation)\w*|randomis\w*|recognis\w*|summaris\w*|characteris\w*|normalis\w*|'
                     r'optimis\w*|utilis\w*|minimis\w*|maximis\w*|emphasis(?:e|ed))\b', re.I)


def check_titles(issues):
    """review round 8 (D-02): every title on one line in PowerPoint. Width in Liberation Sans Bold 24 pt (the PDF's
    stand-in for DengXian; italic runs in Bold Italic); the title box is 13.0 in wide, and the limit 11.9 in leaves room
    for a DengXian up to 9% wider. Reads slide_titles.json, which the build writes with each title's markup."""
    p = os.path.join(HERE, 'slide_titles.json')
    fdir = '/home/b/b37974/opt/lo/root/opt/libreoffice25.8/share/fonts/truetype'
    if not (os.path.exists(p) and os.path.isdir(fdir)):
        issues.append('title widths not checked: slide_titles.json or the Liberation fonts are missing')
        return
    from PIL import ImageFont
    B = ImageFont.truetype(os.path.join(fdir, 'LiberationSans-Bold.ttf'), 1000)
    BI = ImageFont.truetype(os.path.join(fdir, 'LiberationSans-BoldItalic.ttf'), 1000)
    for i, r in enumerate(json.load(open(p)), start=1):
        if 'title' not in r:                      # the cover: a text box of its own, no title placeholder
            continue
        parts = re.split(r'(?<!\\)\*', re.sub(r'[~^]', '', r['title']))
        w = sum((BI if j % 2 else B).getlength(x) for j, x in enumerate(parts)) / 1000 * 24 / 72
        if w > 11.9:
            issues.append(f'slide {i} {r["key"]}: title {w:.2f} in wide (limit 11.9): {r["title"]}')


def main():
    deck = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'platform_depth_summary.pptx')
    stats = json.load(open(os.path.join(HERE, 'stats.json')))
    pool = {norm(x) for x in pool_of(stats)}
    keys = json.load(open(os.path.join(HERE, 'slide_keys.json')))
    slides = deck_parts(deck)
    if len(keys) != len(slides):
        sys.exit(f'slide_keys.json ({len(keys)}) does not match the deck ({len(slides)} slides)')
    bad, issues, n = [], [], 0
    id_rx, census_nums = census()
    if any(k in ID_SLIDES for k in keys) and id_rx is None:
        sys.exit(f'{CENSUS} is missing: the sample IDs of the appendix cannot be checked')
    for i, ((body, notes), key) in enumerate(zip(slides, keys), start=1):
        if id_rx is not None:
            hit = sorted({m.group(0) for t in body + notes for m in id_rx.finditer(t)})
            if hit and key not in ID_SLIDES:
                issues.append(f'slide {i} {key}: {len(hit)} sample ID(s) outside the ID appendix')
            body, notes = [id_rx.sub(' ', t) for t in body], [id_rx.sub(' ', t) for t in notes]
        slide_pool = pool | census_nums if key in ID_SLIDES else pool
        kept, skip = [], False
        for line in notes:
            if line in HEADS:
                skip = line in ('Timing', 'Backup')
            if not skip:
                kept.append(line)
        # file paths in the notes' Source lines name tables, not numbers
        kept = [re.sub(r'\S*/\S*|\S+\.(?:tsv|json|py|png|md)\b', ' ', x) for x in kept]
        for where, texts in (('slide', body), ('notes', kept)):
            for tok in sorted(tokens(' \n '.join(texts))):
                if re.fullmatch(r'\d', tok):
                    continue
                n += 1
                if tok not in slide_pool:
                    bad.append((i, key, where, tok))
        sl, nt = ' \n '.join(body), ' \n '.join(notes)
        for where, t in (('slide', sl), ('notes', nt)):
            for m in BRITISH.finditer(t):
                issues.append(f'slide {i} {where}: British spelling "{m.group(0)}"')
            if re.search(r'diseaselens', t, re.I):
                issues.append(f'slide {i} {where}: DiseaseLens mentioned')
        for m in re.finditer(r'\b\d[\d,]* \+ \d[\d,]*\b', sl):
            issues.append(f'slide {i}: composite number "{m.group(0)}"')
        if key not in VERBATIM:
            if re.search(r'/LARGE|\btmp/|\.py\b|\.tsv\b|\.json\b', sl):
                issues.append(f'slide {i}: internal path on the slide')
        if key not in SOFTWARE_OK:
            if 'SAIGE' in sl:
                issues.append(f'slide {i}: software name SAIGE on the slide')
    check_titles(issues)
    print(f'{n} numbers checked on {len(slides)} slides; {len(bad)} not found in stats.json')
    for i, key, where, tok in bad:
        print(f'  slide {i:2d} {key:12s} {where:5s} {tok}')
    print(f'text rules: {len(issues)} issue(s)')
    for x in issues:
        print('  ' + x)
    sys.exit(1 if bad or issues else 0)


if __name__ == '__main__':
    main()
