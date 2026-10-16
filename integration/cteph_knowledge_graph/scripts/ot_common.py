#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Helpers shared by the cteph_knowledge_graph step scripts: the Open
#           Targets GraphQL client (retries with backoff, as in
#           integration/gene_signals/scripts/opentargets.py), gzipped JSON
#           I/O, the TSV conventions of the other integration components, and
#           the gene themes and population labels of the CTEPH literature
#           (shared by the knowledge base and the workbook).
#           No paths, no project constants: every script receives its inputs
#           on the command line from cteph_knowledge_graph.nf.
# ---------------------------------------------------------------------------
import datetime as dt
import gzip
import json
import re
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd
import yaml

MISSING = '—'   # em dash: "not available / not applicable"
META = '{ meta { apiVersion { x y z } dataVersion { year month iteration } } }'
# Open Targets datasource ids -> the labels printed in the workbook and the knowledge base
DATASOURCES = {'europepmc': 'Literature', 'impc': 'Animal model', 'clinical_precedence': 'Clinical',
               'gwas_credible_sets': 'GWAS'}


# Open Targets / gnomAD population codes of a study's LD reference
ANCESTRY = {'nfe': 'European', 'fin': 'Finnish', 'eas': 'East Asian', 'afr': 'African', 'amr': 'admixed American',
            'sas': 'South Asian', 'csa': 'Central / South Asian', 'mid': 'Middle Eastern', 'asj': 'Ashkenazi Jewish'}
STAGE = {'APPROVAL': 'Approval', 'PHASE_4': 'Phase 4', 'PHASE_3': 'Phase 3', 'PHASE_2_3': 'Phase 2/3',
         'PHASE_2': 'Phase 2', 'PHASE_1_2': 'Phase 1/2', 'PHASE_1': 'Phase 1', 'EARLY_PHASE_1': 'Early phase 1',
         'PRECLINICAL': 'Preclinical', 'UNKNOWN': 'Unknown'}
CAVEATS = [
    'Open Targets associations are scores that aggregate evidence; they are not causal claims.',
    'Literature evidence is text mining (Europe PMC): a gene symbol can be mapped to the wrong gene (for example, '
    '"TF" for tissue factor is mapped to transferrin, TF). The words a target was matched as are listed beside it.',
    'Locus-to-gene (L2G) scores are predictions of the causal gene at a GWAS locus, not proof.',
    'GWAS P values and effects are those of each credible set\'s lead variant.',
    'Mouse-model evidence links a knockout\'s phenotypes to the disease\'s clinical features; it is not a model of '
    'the disease itself.',
    'Related-disease studies overlap in their samples (for example, FinnGen and UK Biobank endpoints); genes are '
    'therefore counted by independent publications, not by studies.',
]
SUP = str.maketrans('0123456789-', '⁰¹²³⁴⁵⁶⁷⁸⁹⁻')


def stage_label(s):
    s = str(s or '')
    return STAGE.get(s, s.replace('_', ' ').capitalize())


def phase_label(phase, stage=''):
    """A trial's phase as registered ('PHASE2/PHASE3', 'EARLY_PHASE1', 'NA') -> 'Phase 2/3' ...; the harmonised
    stage when the registry gives none."""
    s = str(phase or '').upper()
    digits = re.findall(r'PHASE_?(\d)', s)
    if digits:
        return ('Early phase ' if s.startswith('EARLY') else 'Phase ') + '/'.join(digits)
    if s in ('NA', 'N/A', 'NOT_APPLICABLE'):
        return 'phase not applicable'
    return stage_label(stage) if stage and str(stage).upper() != 'UNKNOWN' else 'phase not given'


def ancestry_label(s):
    """'nfe 70%; eas 30%' -> 'European 70%, East Asian 30%'."""
    out = []
    for part in [x.strip() for x in str(s or '').split(';') if x.strip()]:
        code, _, pct = part.partition(' ')
        out.append(f'{ANCESTRY.get(code, code)} {pct}'.strip())
    return ', '.join(out)


def p_unicode(mant, exp):
    """P from Open Targets' mantissa / exponent as 7.0 × 10⁻²⁷ (no float underflow)."""
    if mant is None or exp is None or pd.isna(mant) or pd.isna(exp):
        return ''
    m, e = float(mant), int(exp)
    if m >= 9.95:
        m, e = m / 10, e + 1
    return f'{m:.1f} × 10{str(e).translate(SUP)}'


def pub_label(pid):
    """A literature ID as printed: PubMed IDs as 'PMID n'; Europe PMC preprints (PPR...) and full texts (PMC...)
    by their Europe PMC ID."""
    pid = str(pid or '')
    return f'PMID {pid}' if pid.isdigit() else f'Europe PMC {pid}' if pid else ''


def pub_url(pid):
    pid = str(pid or '')
    if pid.isdigit():
        return f'https://pubmed.ncbi.nlm.nih.gov/{pid}/'
    m = re.match(r'(PPR|PMC)\d+$', pid)
    return f'https://europepmc.org/article/{m.group(1)}/{pid}' if m else ''


def variant_label(vid):
    """Open Targets 9_133261703_A_G -> chr9:133261703:A:G."""
    return 'chr' + str(vid).replace('_', ':') if vid else ''


def log(msg):
    print(f'[{dt.datetime.now():%H:%M:%S}] {msg}', flush=True)


def gql(url, query, variables=None, tries=5):
    """One GraphQL request. Retries server errors, 429 and network failures with exponential backoff;
    any other HTTP error or a GraphQL error raises at once."""
    body = json.dumps({'query': query, 'variables': variables or {}}).encode()
    err = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=body, headers={'Content-Type': 'application/json'})
            d = json.load(urllib.request.urlopen(req, timeout=180))
            if d.get('errors'):
                raise RuntimeError(f'GraphQL errors: {d["errors"]}')
            return d['data']
        except urllib.error.HTTPError as e:
            if e.code < 500 and e.code != 429:
                raise RuntimeError(f'HTTP {e.code}: {e.read().decode()[:800]}')
            err = e
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            err = e
        time.sleep(2 ** k * 3)
    raise RuntimeError(f'Open Targets API unreachable after {tries} tries: {err}')


def release(url):
    m = gql(url, META)['meta']
    dv, av = m['dataVersion'], m['apiVersion']
    data = f'{dv["year"]}.{dv["month"]}' + (f'.{dv["iteration"]}' if dv.get('iteration') else '')
    return {'data_version': data, 'api_version': f'{av["x"]}.{av["y"]}.{av["z"]}'}


def write_json(obj, path):
    path = Path(path)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'wt') as fh:
        json.dump(obj, fh, indent=None if path.suffix == '.gz' else 1, ensure_ascii=False)


def load_json(path):
    path = Path(path)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as fh:
        return json.load(fh)


def read_tsv(path):
    """Every column as text; '' stays '' (never NaN). Callers convert with num()."""
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, na_values=[])


def write_tsv(df, path):
    """P columns in scientific notation; floats with 10 significant digits."""
    out = df.copy()
    for c in out.columns:
        if (c == 'p' or c.endswith('_p')) and pd.api.types.is_numeric_dtype(out[c]):
            out[c] = out[c].map(lambda x: '' if pd.isna(x) else f'{float(x):.6e}')
    out.to_csv(path, sep='\t', index=False, float_format='%.10g')


def num(s):
    return pd.to_numeric(s, errors='coerce')


def csv_list(s):
    """'a,b,c' -> ['a', 'b', 'c'] (the form in which the .nf passes lists)."""
    return [x.strip() for x in str(s).split(',') if x.strip()]


# ── the CTEPH literature: populations and gene themes ───────────────────────
EAS_COUNTRIES = ('Japan', 'China', 'Korea', 'Taiwan', 'Hong Kong', 'Mongolia')
OTHER_GENES, NO_GENE = 'Other genes', 'No gene named'
POSITIVE = ('associated', 'enriched')
CURATED = 'Curated'   # display name of the curated / Mendelian gene evidence ("Mendelian" in the tables)


def study_population(country, ancestry=''):
    """The population a study names: its country, else its ancestry, else 'not stated'."""
    for v in (country, ancestry):
        v = str(v or '').strip()
        if v and v != 'not stated':
            return v
    return 'not stated'


def is_japanese(s):
    """A literature study of a Japanese population: flagged by the extraction, or done in Japan."""
    return s.japanese == 'True' or 'Japan' in str(s.country)


def is_east_asian(s):
    return is_japanese(s) or s.east_asian == 'True' or any(c in str(s.country) for c in EAS_COUNTRIES)


def load_themes(path):
    """The gene themes (config/gene_themes.yaml); a gene may be listed in one theme only."""
    themes = yaml.safe_load(Path(path).read_text())['themes']
    seen = {}
    for t in themes:
        assert {'key', 'name', 'description', 'genes'} <= set(t), t
        for g in t['genes']:
            assert g not in seen, f'{g} is in two themes: {seen[g]}, {t["key"]}'
            seen[g] = t['key']
    return themes


def theme_order(themes):
    return [t['name'] for t in themes] + [OTHER_GENES, NO_GENE]


def assign_themes(findings, themes):
    """The findings with `theme` and `label` (the gene symbol; for a finding without one, its gene text).
    A symbol joins the theme that lists it; a finding without a symbol joins the first theme whose genes or
    words occur in its gene text."""
    by_gene = {g: t['name'] for t in themes for g in t['genes']}

    def hit(word, text, case):
        return re.search(rf'(?<![A-Za-z0-9]){re.escape(word)}(?![A-Za-z0-9])', text, 0 if case else re.I)

    def one(r):
        if r.symbol:
            return by_gene.get(r.symbol, OTHER_GENES)
        for t in themes:
            if any(hit(g, r.gene, True) for g in t['genes']) or any(hit(w, r.gene, False)
                                                                    for w in t.get('words', [])):
                return t['name']
        return NO_GENE
    f = findings.copy()
    f['theme'] = [one(r) for r in f.itertuples()]
    f['label'] = [r.symbol or (r.gene if t != NO_GENE else 'no gene named') for r, t in zip(f.itertuples(), f.theme)]
    return f


def gene_tallies(f):
    """One row per theme x gene label of assign_themes() output: studies by result."""
    rows = []
    for (th, lab), g in f.groupby(['theme', 'label'], sort=False):
        n = lambda res: int(g[g.result.isin(res)].pmid.nunique())  # noqa: E731
        rows.append({'theme': th, 'label': lab, 'symbol': g.symbol.iloc[0], 'n_studies': int(g.pmid.nunique()),
                     'n_findings': len(g), 'associated': n(POSITIVE), 'not_associated': n(('not associated',)),
                     'identified': n(('identified',)), 'other': n(('other',))})
    return pd.DataFrame(rows, columns=['theme', 'label', 'symbol', 'n_studies', 'n_findings', 'associated',
                                       'not_associated', 'identified', 'other'])


def tally_text(r):
    """'associated in 3 studies, not associated in 1, variants identified in 2'."""
    parts = [(r.associated, 'associated in'), (r.not_associated, 'not associated in'),
             (r.identified, 'variants identified in'), (r.other, 'other result in')]
    parts = [f'{t} {n}' for n, t in parts if n]
    if not parts:
        return ''
    word = 'study' if [n for n in (r.associated, r.not_associated, r.identified, r.other) if n][0] == 1 else 'studies'
    return f'{parts[0]} {word}' + ''.join(f', {p}' for p in parts[1:])
