#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Open Targets Platform annotation of every table variant:
#             ClinVar            evidence of the EVA/ClinVar datasource
#             GWAS credible sets credible sets of GWAS studies containing the
#                                variant (posterior probability of the variant,
#                                top locus-to-gene prediction of each set)
#             molQTL credible sets eQTL / pQTL / sQTL / tuQTL (bulk and single-cell), each
#                                classed cis / trans by the variant's distance to the
#                                target gene's TSS (Open Targets' isTransQtl kept as a check)
#             Enhancer-to-gene   E2G predictions linking the variant's interval to genes,
#                                collapsed to one row per target gene (best biosample)
#           One GraphQL query per variant (chr:pos:REF:ALT -> chrom_pos_ref_alt).
#           A variant absent from the Open Targets variant index has no record of
#           any kind there; it is reported as such, not as an error.
# Inputs  : variants.tsv (gene x variant rows; each variant is queried once)
# Origin  : integration/saige_signals/scripts/opentargets.py, unchanged but for the helper
#           module and the de-duplication of variants that map to two genes.
# Outputs : opentargets.tsv         one row per variant (counts + top records, full
#                                   names, records separated by '; ')
#           opentargets_detail.tsv  one row per Open Targets record (E2G: per target gene).
#                                   Credible-set P / beta are the variant's own when
#                                   Open Targets reports them, else the set lead's
#                                   (p_source / beta_source = variant | lead | ''); the
#                                   study's cohorts and population structure as listed
#                                   by Open Targets
#           opentargets_meta.json   API / data release, query settings
#           raw/<variant>.json      the API responses, for audit
# ---------------------------------------------------------------------------
import argparse
import datetime as dt
import json
import re
import sys
import time
import urllib.error
import urllib.request
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gs_common import load_variants, log, write_tsv  # noqa: E402

QTL_TYPES = ['eqtl', 'pqtl', 'sqtl', 'tuqtl', 'sceqtl', 'scpqtl', 'scsqtl', 'sctuqtl']
QTL_LABEL = {'eqtl': 'eQTL', 'pqtl': 'pQTL', 'sqtl': 'sQTL', 'tuqtl': 'tuQTL', 'sceqtl': 'sc-eQTL',
             'scpqtl': 'sc-pQTL', 'scsqtl': 'sc-sQTL', 'sctuqtl': 'sc-tuQTL'}

# pValueMantissa / pValueExponent / beta of a credible set are those of its LEAD variant (`variant`); the
# annotated variant's own statistics, when Open Targets has them, are in its `locus` row. isTransQtl is Open
# Targets' set-level flag (every variant of the set on another chromosome or > 5 Mb from the target gene's
# TSS; null for GWAS); the target's canonical transcript gives the TSS of our own, per-variant class.
CS_FIELDS = '''count rows { studyLocusId studyType pValueMantissa pValueExponent beta finemappingMethod isTransQtl
      variant { id rsIds }
      study { id traitFromSource projectId publicationFirstAuthor publicationDate nSamples condition
              cohorts ldPopulationStructure { ldPopulation relativeSampleSize }
              target { id approvedSymbol canonicalTranscript { chromosome start end strand }
                       genomicLocation { chromosome start end strand } }
              biosample { biosampleName } }
      locus(variantIds: [$id]) { rows { posteriorProbability is95CredibleSet beta pValueMantissa pValueExponent } }
      l2GPredictions(page: {index: 0, size: 1}) { rows { score target { approvedSymbol } } } }'''
QUERY = '''query($id: String!, $index: Int!, $size: Int!) { variant(variantId: $id) { id rsIds
  gwas: credibleSets(studyTypes: [gwas], page: {index: $index, size: $size}) { %s }
  qtl: credibleSets(studyTypes: [%s], page: {index: $index, size: $size}) { %s }
  enhancerToGenes(page: {index: $index, size: $size}) { count rows { target { approvedSymbol } score biosampleName
      datasourceId intervalType distanceToTss } }
} }''' % (CS_FIELDS, ', '.join(QTL_TYPES), CS_FIELDS)
CLINVAR = '''query($id: String!, $cursor: String) { variant(variantId: $id) {
  evidences(datasourceIds: ["eva"], size: 500, cursor: $cursor) { count cursor rows {
      disease { name } diseaseFromSource clinicalSignificances confidence studyId } } } }'''
META = '{ meta { apiVersion { x y z } dataVersion { year month iteration } } }'


def gql(url, query, variables=None, tries=4):
    body = json.dumps({'query': query, 'variables': variables or {}}).encode()
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=body, headers={'Content-Type': 'application/json'})
            d = json.load(urllib.request.urlopen(req, timeout=120))
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


def ot_id(vid):
    c, p, ref, alt = vid.split(':')
    return f'{c.replace("chr", "")}_{p}_{ref}_{alt}'


def fetch_variant(url, oid, size):
    """All pages of the three paginated lists, then the ClinVar cursor."""
    out = {'gwas': [], 'qtl': [], 'enhancerToGenes': [], 'clinvar': [], 'rsIds': []}
    idx, counts = 0, None
    while True:
        d = gql(url, QUERY, {'id': oid, 'index': idx, 'size': size})['variant']
        if d is None:
            return None
        out['rsIds'] = d.get('rsIds') or []
        if counts is None:          # a page past the end of one list reports count 0: keep page 1's
            counts = {k: d[k]['count'] for k in ('gwas', 'qtl', 'enhancerToGenes')}
        for k in counts:
            out[k].extend(d[k]['rows'])
        if all(len(out[k]) >= counts[k] or not d[k]['rows'] for k in counts):
            break
        idx += 1
    cursor = None
    while True:
        e = gql(url, CLINVAR, {'id': oid, 'cursor': cursor})['variant']['evidences']
        out['clinvar'].extend(e['rows'])
        cursor = e.get('cursor')
        if not cursor or len(out['clinvar']) >= e['count'] or not e['rows']:
            break
    out['counts'] = {**counts, 'clinvar': len(out['clinvar'])}
    return out


def pp_of(row):
    rows = (row.get('locus') or {}).get('rows') or []
    return (rows[0]['posteriorProbability'], rows[0]['is95CredibleSet']) if rows else (np.nan, None)


def p_of(row):
    m, e = row.get('pValueMantissa'), row.get('pValueExponent')
    return m * 10.0 ** e if m is not None and e is not None else np.nan


def own_or_lead(variant_value, lead_value, is_lead):
    """(value, source): the variant's own statistic, else the set lead's, else missing."""
    if is_lead and pd.notna(lead_value):
        return lead_value, 'variant'
    if pd.notna(variant_value):
        return variant_value, 'variant'
    if pd.notna(lead_value):
        return lead_value, 'lead'
    return np.nan, ''


def ot_to_vid(oid):
    c, p, ref, alt = oid.split('_')
    return f'chr{c}:{p}:{ref}:{alt}'


QTL_PREFIX = {'cis': 'cis-', 'trans': 'trans-'}       # 'neither' (between the two distances), 'unknown': none


def target_tss(target):
    """('chrN', TSS candidates) of a QTL's target gene: the start of its canonical transcript on the + strand,
    the end on the - strand (the gene's own span without a canonical transcript); both ends when the strand
    is unknown (the nearer one is used)."""
    loc = (target or {}).get('canonicalTranscript') or (target or {}).get('genomicLocation')
    if not loc:
        return None, ()
    s, e = int(loc['start']), int(loc['end'])
    return f'chr{loc["chromosome"]}', {'POSITIVE': (s,), 'NEGATIVE': (e,)}.get(loc.get('strand'), (s, e))


def qtl_class(vid, chrom, tss, cis_window, trans_distance):
    """(class, |distance| bp, TSS used) of the variant against the target gene: cis, same chromosome within
    cis_window of the TSS; trans, another chromosome or beyond trans_distance; neither, in between; unknown,
    no gene location."""
    if chrom is None:
        return 'unknown', np.nan, np.nan
    vc, vp = vid.split(':')[:2]
    if vc != chrom:
        return 'trans', np.nan, tss[0] if len(tss) == 1 else np.nan
    t = min(tss, key=lambda x: abs(int(vp) - x))
    d = abs(int(vp) - t)
    return ('cis' if d <= cis_window else 'trans' if d > trans_distance else 'neither'), d, t


def pp_text(pp):
    return f'{pp:.2f}' if pp >= 0.005 else '<0.01'


def full(s):
    """Names are never truncated; ';' becomes ',' so that '; ' stays a safe record separator."""
    return re.sub(r'\s+', ' ', str(s)).strip().replace(';', ',')


def ld_text(structure):
    """'NFE 96%, EAS 2%, AFR 2%': the study sample's share by gnomAD reference population, as Open Targets
    assigns them for LD, largest first ('<1%' below half a percent; the code alone when no share is given)."""
    def share(x):
        f = x.get('relativeSampleSize')
        return '' if f is None else ' <1%' if f < 0.005 else f' {f:.0%}'
    rows = sorted(structure or [], key=lambda x: -(x.get('relativeSampleSize') or 0))
    return ', '.join(f'{x["ldPopulation"].upper()}{share(x)}' for x in rows)


def detail_rows(vid, rec, focus, cis_window, trans_distance):
    rows = []
    for cat, key in (('GWAS', 'gwas'), ('molQTL', 'qtl')):
        for r in rec[key]:
            st = r.get('study') or {}
            if cat == 'molQTL':
                t_chrom, t_tss = target_tss(st.get('target'))
                cls, dist, tss = qtl_class(vid, t_chrom, t_tss, cis_window, trans_distance)
            else:
                t_chrom, cls, dist, tss = None, '', np.nan, np.nan
            pp, is95 = pp_of(r)
            l2g = ((r.get('l2GPredictions') or {}).get('rows') or [{}])[0]
            trait = st.get('traitFromSource') or ''
            lead = r.get('variant') or {}
            lead_vid = ot_to_vid(lead['id']) if lead.get('id') else ''
            is_lead = lead_vid == vid
            own = ((r.get('locus') or {}).get('rows') or [{}])[0]
            var_p = p_of(own)
            var_beta = np.nan if own.get('beta') is None else own.get('beta')
            lead_p, lead_beta = p_of(r), (np.nan if r.get('beta') is None else r.get('beta'))
            p_val, p_src = own_or_lead(var_p, lead_p, is_lead)
            b_val, b_src = own_or_lead(var_beta, lead_beta, is_lead)
            rows.append({
                'variant_id': vid, 'category': cat,
                'type': 'GWAS' if cat == 'GWAS' else QTL_LABEL.get(r.get('studyType'), r.get('studyType')),
                'trait_or_target': trait if cat == 'GWAS' else ((st.get('target') or {}).get('approvedSymbol') or trait),
                'biosample': (st.get('biosample') or {}).get('biosampleName') or '',
                'condition': st.get('condition') or '',
                'study_id': st.get('id') or '', 'project': st.get('projectId') or '',
                'first_author': st.get('publicationFirstAuthor') or '',
                'year': (st.get('publicationDate') or '')[:4], 'n_samples': st.get('nSamples'),
                'cohorts': ', '.join(full(c) for c in st.get('cohorts') or []),
                'ld_population': ld_text(st.get('ldPopulationStructure')),
                'p_value': p_val, 'beta': b_val, 'p_source': p_src, 'beta_source': b_src,
                'lead_variant': lead_vid, 'lead_rsid': ';'.join(lead.get('rsIds') or []), 'lead_is_variant': is_lead,
                'variant_p': var_p, 'variant_beta': var_beta, 'lead_p': lead_p, 'lead_beta': lead_beta,
                'pp': pp, 'in_95cs': is95, 'method': r.get('finemappingMethod') or '',
                'l2g_gene': (l2g.get('target') or {}).get('approvedSymbol') or '', 'l2g_score': l2g.get('score'),
                'target_chrom': t_chrom or '', 'target_tss': tss, 'qtl_distance': dist, 'qtl_class': cls,
                'ot_is_trans': r.get('isTransQtl'),
                'e2g_score': np.nan, 'interval_type': '', 'distance_to_tss': np.nan,
                'clinical_significance': '', 'review_status': '', 'datasource': r.get('studyLocusId') or '',
                'focus': bool(cat == 'GWAS' and focus.search(trait))})
    # E2G predictions come per gene x biosample (thousands for a promoter variant): one row per
    # target gene, the best-scoring biosample, and the number of biosamples predicting the link
    best = {}
    for r in rec['enhancerToGenes']:
        gene = (r.get('target') or {}).get('approvedSymbol') or ''
        b = best.setdefault(gene, {'n': 0, 'bios': set(), 'r': r})
        b['n'] += 1
        b['bios'].add(r.get('biosampleName') or '')
        if (r.get('score') or 0) > (b['r'].get('score') or 0):
            b['r'] = r
    for gene, b in best.items():
        r = b['r']
        rows.append({'variant_id': vid, 'category': 'E2G', 'type': 'E2G', 'trait_or_target': gene,
                     'biosample': r.get('biosampleName') or '', 'condition': '', 'study_id': '', 'project': '',
                     'first_author': '', 'year': '', 'n_samples': np.nan, 'p_value': np.nan, 'beta': np.nan,
                     'pp': np.nan, 'in_95cs': None, 'method': '', 'l2g_gene': '', 'l2g_score': np.nan,
                     'e2g_score': r.get('score'), 'n_biosamples': len(b['bios']), 'n_predictions': b['n'],
                     'interval_type': r.get('intervalType') or '', 'distance_to_tss': r.get('distanceToTss'),
                     'clinical_significance': '', 'review_status': '', 'datasource': r.get('datasourceId') or '',
                     'focus': False})
    for r in rec['clinvar']:
        rows.append({'variant_id': vid, 'category': 'ClinVar', 'type': 'ClinVar',
                     'trait_or_target': (r.get('disease') or {}).get('name') or r.get('diseaseFromSource') or '',
                     'biosample': '', 'condition': r.get('diseaseFromSource') or '', 'study_id': r.get('studyId') or '',
                     'project': '', 'first_author': '', 'year': '', 'n_samples': np.nan, 'p_value': np.nan,
                     'beta': np.nan, 'pp': np.nan, 'in_95cs': None, 'method': '', 'l2g_gene': '', 'l2g_score': np.nan,
                     'e2g_score': np.nan, 'interval_type': '', 'distance_to_tss': np.nan,
                     'clinical_significance': '; '.join(r.get('clinicalSignificances') or []),
                     'review_status': r.get('confidence') or '', 'datasource': 'eva', 'focus': False})
    return rows


def summarise(vid, rec, det, top_n):
    base = {'variant_id': vid, 'ot_in_index': rec is not None}
    empty = {'ot_rsids': '', 'ot_clinvar_n': 0, 'ot_clinvar': '', 'ot_gwas_n': 0, 'ot_gwas_focus_n': 0,
             'ot_gwas_top': '', 'ot_l2g': '', 'ot_qtl_n': 0, 'ot_qtl_trans_n': 0, 'ot_qtl_top': '', 'ot_e2g_n': 0,
             'ot_e2g_top': ''}
    if rec is None:
        return {**base, **empty}
    d = pd.DataFrame(det) if det else pd.DataFrame(columns=['category'])
    out = {**base, **empty, 'ot_rsids': ';'.join(rec['rsIds'])}
    g = d[d.category == 'GWAS'].copy() if len(d) else d
    if len(g):
        out['ot_gwas_n'] = int(g.datasource.nunique())
        out['ot_gwas_focus_n'] = int(g[g.focus].datasource.nunique())
        # by posterior probability alone: thrombosis-related traits are counted, not ranked first
        g = g.assign(_pp=g.pp.fillna(-1)).sort_values('_pp', ascending=False, kind='stable')
        seen, parts = set(), []
        for _, r in g.iterrows():
            k = r.trait_or_target.lower()
            if k in seen:
                continue
            seen.add(k)
            parts.append(f'{full(r.trait_or_target)} ({pp_text(r.pp)})' if np.isfinite(r._pp) and r._pp >= 0
                         else full(r.trait_or_target))
            if len(parts) == top_n:
                break
        out['ot_gwas_top'] = '; '.join(parts)
        l2g = g[g.l2g_gene != '']
        if len(l2g):
            # most frequent top gene; ties to the higher median score, then the name (order-independent)
            med = {k: float(np.nanmedian(x.l2g_score.astype(float))) for k, x in l2g.groupby('l2g_gene')}
            gene = min(Counter(l2g.l2g_gene).items(), key=lambda kv: (-kv[1], -med[kv[0]], kv[0]))[0]
            out['ot_l2g'] = f'{gene} ({med[gene]:.2f})'
    q = d[d.category == 'molQTL'].copy() if len(d) else d
    if len(q):
        out['ot_qtl_n'] = int(q.datasource.nunique())
        out['ot_qtl_trans_n'] = int(q[q.qtl_class == 'trans'].datasource.nunique())   # counted, not ranked
        q = q.assign(_pp=q.pp.fillna(-1)).sort_values('_pp', ascending=False, kind='stable')
        seen, parts = set(), []
        for _, r in q.iterrows():
            k = (r.type, r.trait_or_target, r.biosample)
            if k in seen:
                continue
            seen.add(k)
            parts.append(f'{QTL_PREFIX.get(r.qtl_class, "")}{r.type} {r.trait_or_target} '
                         f'({full(r.biosample)}, {pp_text(r.pp)})')
            if len(parts) == top_n:
                break
        out['ot_qtl_top'] = '; '.join(parts)
    e = d[d.category == 'E2G'].copy() if len(d) else d
    if len(e):
        out['ot_e2g_n'] = int(e.trait_or_target.nunique())
        e = e.sort_values('e2g_score', ascending=False, kind='stable')
        out['ot_e2g_top'] = '; '.join(
            f'{r.trait_or_target} {r.e2g_score:.2f} ({full(r.biosample)}'
            + (f', {int(r.n_biosamples)} tissues / cell types)' if r.n_biosamples > 1 else ')')
            for _, r in e.head(top_n).iterrows())
    c = d[d.category == 'ClinVar'] if len(d) else d
    if len(c):
        out['ot_clinvar_n'] = len(c)
        pairs = sorted({(r.clinical_significance or 'not provided', r.trait_or_target or 'not specified')
                        for _, r in c.iterrows()})
        out['ot_clinvar'] = '; '.join(f'{s}: {full(t)}' for s, t in pairs[:top_n])    # the count says how many
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--variants', required=True)
    ap.add_argument('--api', required=True)
    ap.add_argument('--enabled', default='true')
    ap.add_argument('--focus-traits', required=True,
                    help='regex; GWAS sets whose trait matches are counted as thrombosis-related')
    ap.add_argument('--top-n', type=int, default=3, help='records per summary, by posterior probability')
    ap.add_argument('--page-size', type=int, default=100, help='larger pages are refused as too expensive')
    ap.add_argument('--cis-window', type=int, default=1_000_000,
                    help='molQTL cis: the variant within this distance (bp) of the target gene TSS')
    ap.add_argument('--trans-distance', type=int, default=5_000_000,
                    help='molQTL trans: another chromosome, or farther than this (bp); in between: neither')
    ap.add_argument('--out', required=True)
    ap.add_argument('--out-detail', required=True)
    ap.add_argument('--out-meta', required=True)
    a = ap.parse_args()
    v = load_variants(a.variants).drop_duplicates('variant_id').reset_index(drop=True)   # one query per variant
    focus = re.compile(a.focus_traits, re.I)
    enabled = str(a.enabled).lower() in ('true', '1', 'yes')
    raw = Path('raw')
    raw.mkdir(exist_ok=True)
    assert 0 < a.cis_window <= a.trans_distance, (a.cis_window, a.trans_distance)
    meta = {'queried': enabled, 'api': a.api, 'focus_traits': a.focus_traits, 'top_n': a.top_n,
            'cis_window': a.cis_window, 'trans_distance': a.trans_distance,
            'queried_at': dt.datetime.now().isoformat(timespec='seconds')}
    summ, det = [], []
    if enabled:
        m = gql(a.api, META)['meta']
        meta['api_version'] = '.'.join(str(m['apiVersion'][k]) for k in ('x', 'y', 'z'))
        dv = m['dataVersion']
        meta['data_version'] = f'{dv["year"]}.{dv["month"]}' + (f'.{dv["iteration"]}' if dv.get('iteration') else '')
        for i, vid in enumerate(v.variant_id, 1):
            rec = fetch_variant(a.api, ot_id(vid), a.page_size)
            (raw / f'{vid.replace(":", "_")}.json').write_text(json.dumps(rec))
            rows = detail_rows(vid, rec, focus, a.cis_window, a.trans_distance) if rec else []
            if rec:
                for cat, key in (('GWAS', 'gwas'), ('molQTL', 'qtl'), ('E2G', 'enhancerToGenes'), ('ClinVar', 'clinvar')):
                    got = sum(r.get('n_predictions', 1) for r in rows if r['category'] == cat)
                    assert got == rec['counts'][key], (vid, cat, got, rec['counts'][key])
            det.extend(rows)
            summ.append(summarise(vid, rec, rows, a.top_n))
            if i % 25 == 0:
                log(f'  {i}/{len(v)} variants queried')
    else:
        summ = [{'variant_id': vid, 'ot_in_index': None} for vid in v.variant_id]
    s = pd.DataFrame(summ)
    d = pd.DataFrame(det, columns=['variant_id', 'category', 'type', 'trait_or_target', 'biosample', 'condition',
                                   'study_id', 'project', 'first_author', 'year', 'n_samples', 'cohorts',
                                   'ld_population', 'p_value', 'beta',
                                   'p_source', 'beta_source', 'lead_variant', 'lead_rsid', 'lead_is_variant',
                                   'variant_p', 'variant_beta', 'lead_p', 'lead_beta', 'pp', 'in_95cs', 'method', 'l2g_gene', 'l2g_score',
                                   'target_chrom', 'target_tss', 'qtl_distance', 'qtl_class', 'ot_is_trans', 'e2g_score', 'n_biosamples',
                                   'n_predictions', 'interval_type', 'distance_to_tss', 'clinical_significance',
                                   'review_status', 'datasource', 'focus'])
    assert len(s) == len(v) and s.variant_id.is_unique
    if enabled:
        absent = s[~s.ot_in_index.astype(bool)]
        assert (absent[['ot_gwas_n', 'ot_qtl_n', 'ot_e2g_n', 'ot_clinvar_n']] == 0).all().all()
        meta.update({'n_variants': len(s), 'n_in_index': int(s.ot_in_index.sum()),
                     'n_with_gwas': int((s.ot_gwas_n > 0).sum()), 'n_with_focus_gwas': int((s.ot_gwas_focus_n > 0).sum()),
                     'n_with_qtl': int((s.ot_qtl_n > 0).sum()), 'n_with_e2g': int((s.ot_e2g_n > 0).sum()),
                     'n_with_clinvar': int((s.ot_clinvar_n > 0).sum()), 'n_records': len(d)})
        cs = d[d.category.isin(['GWAS', 'molQTL'])]
        meta['stat_sources'] = {k: {(src or 'missing'): int(n) for src, n in cs[k].fillna('').value_counts().items()}
                                for k in ('p_source', 'beta_source')}
        log(f'credible-set records {len(cs)}: P from {meta["stat_sources"]["p_source"]}; '
            f'beta from {meta["stat_sources"]["beta_source"]}')
        # the studies' cohorts and population structure, as far as Open Targets reports them
        meta['study_annotation'] = {
            cat: {'studies': int(g.study_id.nunique()),
                  'with_cohorts': int(g[g.cohorts.fillna('') != ''].study_id.nunique()),
                  'with_population': int(g[g.ld_population.fillna('') != ''].study_id.nunique()),
                  'records_without_cohorts': int((g.cohorts.fillna('') == '').sum()),
                  'records_without_population': int((g.ld_population.fillna('') == '').sum())}
            for cat, g in cs.groupby('category')}
        log(f'studies: {meta["study_annotation"]}')
        # cis / trans of the molQTL records against Open Targets' set-level flag (shown where they contradict)
        q = d[d.category == 'molQTL']
        flag = q.ot_is_trans.map({True: 'true', False: 'false'}).fillna('null')
        meta['qtl_class'] = {k: int(n) for k, n in q.qtl_class.value_counts().items()}
        meta['qtl_class_vs_isTransQtl'] = {f'{c} / {f}': int(n) for (c, f), n in sorted(Counter(zip(q.qtl_class, flag)).items())}
        meta['qtl_contradictions'] = int(((q.qtl_class == 'cis') & (flag == 'true')).sum()
                                         + ((q.qtl_class == 'trans') & (flag == 'false')).sum())
        log(f'molQTL records {len(q)}: {meta["qtl_class"]} (cis <= {a.cis_window:,} bp from the target TSS; trans, '
            f'another chromosome or > {a.trans_distance:,} bp); against isTransQtl {meta["qtl_class_vs_isTransQtl"]}, '
            f'{meta["qtl_contradictions"]} contradicting')
        log(f'Open Targets {meta["api_version"]} (data {meta["data_version"]}): {meta["n_in_index"]}/{len(s)} in the '
            f'variant index; GWAS {meta["n_with_gwas"]} (thrombosis-related {meta["n_with_focus_gwas"]}), molQTL '
            f'{meta["n_with_qtl"]}, E2G {meta["n_with_e2g"]}, ClinVar {meta["n_with_clinvar"]}; {len(d)} records')
    write_tsv(s, a.out)
    write_tsv(d, a.out_detail)
    with open(a.out_meta, 'w') as fh:
        json.dump(meta, fh, indent=2)


if __name__ == '__main__':
    main()
