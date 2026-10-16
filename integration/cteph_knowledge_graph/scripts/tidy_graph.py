#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The long tables and the graph of one snapshot, every count asserted
#           against its source:
#             conditions.tsv         the conditions, why they relate to CTEPH, their
#                                    genetic data (GWAS studies, credible sets, East
#                                    Asian studies, Mendelian genes)
#             relations.tsv          one row per piece of relation evidence
#             colocalisation.tsv     CTEPH locus x colocalising GWAS credible set
#             gwas_studies.tsv       every GWAS study of CTEPH and the conditions:
#                                    design (n, cases, controls, cohorts), ancestry,
#                                    East Asian flag, summary statistics
#             gwas_loci.tsv          every credible set: lead variant, its own P / beta,
#                                    fine-mapping, top three locus-to-gene genes
#             mendelian.tsv          condition x gene x source: score, classifications
#                                    (ClinGen / Gene2Phenotype class, PanelApp rating, ...),
#                                    ClinVar review status and clinical significance,
#                                    specific diseases, PMIDs
#             mendelian_records.tsv  every curated / Mendelian / burden evidence row
#             cteph_targets.tsv      CTEPH's associated targets, their literature score
#                                    and Reactome top-level themes
#             literature.tsv         text-mined target x publication with sentences
#             cteph_studies / findings / methods / screening .tsv
#                                    the genetic studies of CTEPH (extract_genetics.py)
#             gene_condition.tsv     gene x condition: the kinds of genetic evidence
#             nodes.tsv, edges.tsv, graph.json   the knowledge graph
#             tidy.json              counts and release
# Inputs  : the fetch, resolve and extraction outputs
# Outputs : as above, in --out-dir
# ---------------------------------------------------------------------------
import argparse
import json
import math
import re
import shutil
import sys
from collections import OrderedDict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import load_json, log, num, read_tsv, write_tsv  # noqa: E402

CORE = 'EFO_0022954'
SECTION_RANK = {'abstract': 0, 'results': 1, 'concl': 2, 'discuss': 3, 'intro': 4, 'other': 5, 'table': 6}
SOURCE_LABEL = {'clingen': 'ClinGen', 'gene2phenotype': 'Gene2Phenotype', 'genomics_england': 'PanelApp',
                'orphanet': 'Orphanet', 'eva': 'ClinVar', 'uniprot_literature': 'UniProt literature',
                'uniprot_variants': 'UniProt variants', 'gene_burden': 'Gene burden'}


def sentences(ev, k):
    """(title, words the target was matched as, n distinct sentences, the k best sentences)."""
    title, seen, matched, rows = '', set(), OrderedDict(), []
    for i, s in enumerate(ev.get('textMiningSentences') or []):
        text = ' '.join(re.sub(r'<[^>]+>', ' ', s['text']).split())
        if s.get('tStart') is not None and s.get('tEnd') is not None:
            term = s['text'][s['tStart']:s['tEnd']].strip()
            if term:
                matched[term] = None
        if s['section'] == 'title':
            title = title or text
            continue
        if text in seen:
            continue
        seen.add(text)
        rows.append((SECTION_RANK.get(s['section'], 7), i, text))
    return title, list(matched), len(rows), [t for _r, _i, t in sorted(rows)[:k]]


def ancestry(pops):
    return '; '.join(f'{p["ldPopulation"]} {100 * p["relativeSampleSize"]:.0f}%'
                     for p in sorted(pops or [], key=lambda p: -p['relativeSampleSize']))


def eas_share(pops):
    return sum(p['relativeSampleSize'] for p in pops or [] if p['ldPopulation'] == 'eas')


def log10p(m, e):
    return math.log10(m) + e if m and m > 0 and e is not None else None


def study_rows(node, name, studies, eas_min):
    srows, lrows = [], []
    for s in studies:
        pops = s.get('ldPopulationStructure') or []
        base = {'node_id': node, 'condition': name, 'study_id': s['id'], 'project': s.get('projectId') or '',
                'trait': s.get('traitFromSource') or '', 'first_author': s.get('publicationFirstAuthor') or '',
                'year': (s.get('publicationDate') or '')[:4], 'journal': s.get('publicationJournal') or '',
                'pmid': s.get('pubmedId') or '', 'ancestry': ancestry(pops),
                'eas': eas_share(pops) >= eas_min}
        srows.append({**base, 'n_samples': s.get('nSamples'), 'n_cases': s.get('nCases'),
                      'n_controls': s.get('nControls'), 'cohorts': '; '.join(s.get('cohorts') or []),
                      'summary_statistics': bool(s.get('hasSumstats')),
                      'n_credible_sets': len(s['credible_sets'])})
        for c in s['credible_sets']:
            v = c.get('variant') or {}
            l2g = (c.get('l2GPredictions') or {}).get('rows') or []
            row = {**base, 'study_locus_id': c['studyLocusId'], 'variant_id': v.get('id', ''),
                   'rsid': ' '.join(v.get('rsIds') or []), 'p_mantissa': c.get('pValueMantissa'),
                   'p_exponent': c.get('pValueExponent'),
                   'log10_p': log10p(c.get('pValueMantissa'), c.get('pValueExponent')), 'beta': c.get('beta'),
                   'finemapping': c.get('finemappingMethod') or '', 'confidence': c.get('confidence') or ''}
            for k in range(3):
                row[f'l2g{k + 1}_symbol'] = l2g[k]['target']['approvedSymbol'] if k < len(l2g) else ''
                row[f'l2g{k + 1}_ensembl'] = l2g[k]['target']['id'] if k < len(l2g) else ''
                row[f'l2g{k + 1}_score'] = l2g[k]['score'] if k < len(l2g) else None
            lrows.append(row)
    return srows, lrows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cteph-raw', required=True)
    ap.add_argument('--conditions', required=True, help='resolve.py outputs')
    ap.add_argument('--diseases', nargs='*', default=[], help='disease.<id>.json.gz')
    ap.add_argument('--extract', required=True, help='extract_genetics.py outputs')
    ap.add_argument('--sentences', type=int, default=3)
    ap.add_argument('--eas-min', type=float, default=0.5, help='East Asian share of a study to flag it')
    ap.add_argument('--l2g-min', type=float, default=0.5, help='top locus-to-gene score for a gene-condition link')
    ap.add_argument('--out-dir', default='.')
    a = ap.parse_args()
    raw, R, X, out = Path(a.cteph_raw), Path(a.conditions), Path(a.extract), Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    meta, fc = load_json(raw / 'meta.json'), load_json(raw / 'fetch_cteph.json')
    assoc, ev = load_json(raw / 'associations.json.gz'), load_json(raw / 'evidence.json.gz')
    tg, cg = load_json(raw / 'targets.json.gz'), load_json(raw / 'gwas.json.gz')
    cond, rel, co = read_tsv(R / 'conditions.tsv'), read_tsv(R / 'relations.tsv'), read_tsv(R / 'colocalisation.tsv')
    assert len(assoc) == fc['associations'] and len(ev) == fc['evidence']
    cteph_name = load_json(raw / 'disease.json')['name']

    # ── CTEPH targets, themes, text-mined literature ─────────────────────────
    themes = {t['id']: sorted({p['topLevelTerm'] for p in t.get('pathways') or [] if p.get('topLevelTerm')})
              for t in tg}
    lit = []
    for e in (x for x in ev if x['datasourceId'] == 'europepmc'):
        title, matched, n, best = sentences(e, a.sentences)
        for pmid in e['literature'] or ['']:
            lit.append({'ensembl_id': e['target']['id'], 'symbol': e['target']['approvedSymbol'], 'pmid': pmid,
                        'year': e.get('publicationYear'), 'title': title, 'matched_as': ' | '.join(matched),
                        'n_sentences': n, 'sentences': ' || '.join(best), 'evidence_score': e['score']})
    lit = pd.DataFrame(lit).sort_values(['symbol', 'year', 'pmid'], ascending=[True, False, True])
    assert not lit.duplicated(['ensembl_id', 'pmid']).any()
    ln = lit.groupby('ensembl_id').agg(n_pmids=('pmid', 'nunique'), first_year=('year', 'min'),
                                       last_year=('year', 'max'))
    trows = []
    for k, r in enumerate(sorted(assoc, key=lambda r: (-r['score'], r['target']['approvedSymbol'])), 1):
        t, ds = r['target'], {x['id']: x['score'] for x in r['datasourceScores']}
        trows.append({'rank': k, 'ensembl_id': t['id'], 'symbol': t['approvedSymbol'], 'name': t.get('approvedName')
                      or '', 'score': r['score'], 'ds_europepmc': ds.get('europepmc'),
                      'ds_gwas_credible_sets': ds.get('gwas_credible_sets'),
                      'n_pmids': int(ln.n_pmids.get(t['id'], 0)), 'first_year': ln.first_year.get(t['id']),
                      'last_year': ln.last_year.get(t['id']), 'themes': '; '.join(themes.get(t['id'], []))})
    targets = pd.DataFrame(trows)

    # ── GWAS of CTEPH and of the conditions ─────────────────────────────────
    sr, lr = study_rows(CORE, cteph_name, cg['studies'], a.eas_min)
    assert len(lr) == fc['credible_sets']
    mrows, mrec = [], []
    names = dict(zip(cond.node_id, cond.name))
    fetched = set()
    for f in sorted(a.diseases):
        d = load_json(f)
        node = d['disease']['id']
        assert d['release']['data_version'] == meta['data_version'], f'{f}: another release'
        fetched.add(node)
        s1, l1 = study_rows(node, names.get(node, d['disease']['name']), d['studies'], a.eas_min)
        assert len(l1) == d['n_credible_sets']
        sr += s1
        lr += l1
        for e in d['mendelian_evidence']:
            src = SOURCE_LABEL.get(e['datasourceId'], e['datasourceId'])
            conf = e.get('confidence') or ''   # ClinVar: the review status, not a classification
            mrec.append({'node_id': node, 'condition': names.get(node, d['disease']['name']), 'source': src,
                         'ensembl_id': e['target']['id'], 'symbol': e['target']['approvedSymbol'],
                         'specific_disease': (e.get('disease') or {}).get('name', ''),
                         'disease_from_source': e.get('diseaseFromSource') or '', 'score': e['score'],
                         'confidence': '' if src == 'ClinVar' else conf,
                         'review_status': conf if src == 'ClinVar' else '',
                         'clinical_significance': '; '.join(e.get('clinicalSignificances') or []),
                         'allele_origin': '; '.join(e.get('alleleOrigins') or []),
                         'variant': ((e.get('variant') or {}).get('rsIds') or [''])[0] or e.get('variantRsId') or
                         (e.get('variant') or {}).get('id', ''),
                         'pmids': ' '.join(e.get('literature') or []), 'cohort': e.get('cohortShortName') or '',
                         'statistical_method': e.get('statisticalMethod') or '',
                         'p_mantissa': e.get('pValueMantissa'), 'p_exponent': e.get('pValueExponent'),
                         'odds_ratio': e.get('oddsRatio'), 'ancestry': e.get('ancestry') or '',
                         'sample_size': e.get('studySampleSize')})
        for ds, ts in d['mendelian_targets'].items():
            for t in ts:
                mrows.append({'node_id': node, 'condition': names.get(node, d['disease']['name']),
                              'source': SOURCE_LABEL.get(ds, ds), 'ensembl_id': t['id'], 'symbol': t['symbol'],
                              'source_score': t['score']})
    want = set(cond[cond.genetics.eq('full') & cond.ot_ids.ne('')].node_id)
    assert fetched == want, (sorted(want - fetched), sorted(fetched - want))
    studies, loci = pd.DataFrame(sr), pd.DataFrame(lr)
    mrec = pd.DataFrame(mrec)
    M = pd.DataFrame(mrows)
    if len(M):
        agg = mrec.groupby(['node_id', 'source', 'ensembl_id']).agg(
            n_records=('score', 'size'),
            classifications=('confidence', lambda x: '; '.join(sorted({v for v in x if v}))),
            review_status=('review_status', lambda x: '; '.join(sorted({v for v in x if v}))),
            clinical_significance=('clinical_significance',
                                   lambda x: '; '.join(sorted({s for v in x for s in v.split('; ') if s}))),
            specific_diseases=('specific_disease', lambda x: '; '.join(sorted({v for v in x if v}))),
            pmids=('pmids', lambda x: ' '.join(sorted({p for v in x for p in v.split() if p}))),
            best_log10_p=('p_mantissa', lambda x: None)).reset_index()
        M = M.merge(agg, on=['node_id', 'source', 'ensembl_id'], how='left')
        M['n_records'] = M.n_records.fillna(0).astype(int)
        M = M.sort_values(['condition', 'symbol', 'source']).drop(columns='best_log10_p')

    # ── colocalisation: East Asian flag of the colocalising study ────────────
    co['eas'] = co.ancestry.map(lambda s: sum(float(x.split()[1].rstrip('%')) / 100 for x in s.split('; ')
                                              if x.startswith('eas ')) >= a.eas_min if s else False)

    # ── genetic studies of CTEPH (extraction) ────────────────────────────────
    cst, cfi = read_tsv(X / 'studies.tsv'), read_tsv(X / 'findings.tsv')
    cmf, csc = read_tsv(X / 'methods.tsv'), read_tsv(X / 'screening.tsv')

    # ── gene x condition evidence ────────────────────────────────────────────
    gc = []
    for r in loci.itertuples():
        if r.l2g1_symbol and num(pd.Series([r.l2g1_score]))[0] >= a.l2g_min:
            gc.append({'ensembl_id': r.l2g1_ensembl, 'symbol': r.l2g1_symbol, 'node_id': r.node_id,
                       'condition': r.condition, 'evidence': 'GWAS', 'detail': f'{r.rsid or r.variant_id} '
                       f'({r.study_id}, L2G {float(r.l2g1_score):.2f})', 'eas': r.eas})
    for r in M.itertuples() if len(M) else []:
        gc.append({'ensembl_id': r.ensembl_id, 'symbol': r.symbol, 'node_id': r.node_id, 'condition': r.condition,
                   'evidence': 'Mendelian', 'detail': f'{r.source} ({r.classifications or r.clinical_significance or "score " + format(float(r.source_score), ".2f")})',
                   'eas': False})
    jp = dict(zip(cst.pmid, cst.japanese)) if len(cst) else {}
    for r in cfi.itertuples():
        if r.ensembl_id and r.result in ('associated', 'identified', 'enriched'):
            gc.append({'ensembl_id': r.ensembl_id, 'symbol': r.symbol, 'node_id': CORE, 'condition': cteph_name,
                       'evidence': 'Literature', 'detail': f'PMID {r.pmid}: {r.result} {r.variant}'.strip(),
                       'eas': jp.get(r.pmid) == 'True'})
    GC = pd.DataFrame(gc)
    GCs = (GC.groupby(['ensembl_id', 'symbol', 'node_id', 'condition', 'evidence'])
           .agg(n=('detail', 'size'), details=('detail', lambda x: '; '.join(list(dict.fromkeys(x))[:5])),
                eas=('eas', 'any')).reset_index())

    # ── the graph ────────────────────────────────────────────────────────────
    nodes, edges = {}, []

    def node(i, label, kind, **kw):
        if i not in nodes:
            nodes[i] = {'id': i, 'label': label, 'type': kind, **kw}

    node(CORE, cteph_name, 'condition', category='CTEPH')
    for r in cond.itertuples():
        node(r.node_id, r.name, 'condition', category=r.category)
    for r in rel.itertuples():
        src = f'PMID:{r.pmid}' if r.pmid else r.citation
        if r.pmid:
            node(src, r.citation, 'publication')
        edges.append({'source': r.node_id, 'target': CORE, 'relation': f'related ({r.evidence_type})',
                      'detail': r.detail, 'evidence': src, 'strength': r.strength})
    for r in GCs.itertuples():
        node(r.ensembl_id, r.symbol, 'gene')
        edges.append({'source': r.ensembl_id, 'target': r.node_id, 'relation': f'genetic evidence ({r.evidence})',
                      'detail': r.details, 'evidence': r.evidence, 'strength': r.n})
    for r in co.itertuples():
        a_, b_ = f'variant:{r.cteph_locus}', f'variant:{r.other_lead}'
        node(a_, r.cteph_locus, 'variant')
        node(b_, r.other_lead, 'variant')
        edges.append({'source': a_, 'target': b_, 'relation': 'colocalises with',
                      'detail': f'{r.condition or r.trait} ({r.study_id})', 'evidence': r.method,
                      'strength': r.clpp})
    for r in cst.itertuples():
        node(f'PMID:{r.pmid}', f'{r.first_author} {r.year}', 'publication')
        edges.append({'source': f'PMID:{r.pmid}', 'target': CORE, 'relation': 'genetic study of',
                      'detail': r.design, 'evidence': 'PubMed', 'strength': r.n_cases})
    N, E = pd.DataFrame(nodes.values()), pd.DataFrame(edges)

    # ── write ────────────────────────────────────────────────────────────────
    cond['n_gwas_studies'] = cond.node_id.map(studies.groupby('node_id').size()).fillna(0).astype(int)
    cond['n_credible_sets'] = cond.node_id.map(loci.groupby('node_id').size()).fillna(0).astype(int)
    cond['n_eas_studies'] = cond.node_id.map(studies[studies.eas].groupby('node_id').size()).fillna(0).astype(int)
    cond['n_mendelian_genes'] = cond.node_id.map(M.groupby('node_id').ensembl_id.nunique() if len(M)
                                                 else pd.Series(dtype=int)).fillna(0).astype(int)
    for name, df in [('conditions', cond), ('relations', rel), ('colocalisation', co), ('gwas_studies', studies),
                     ('gwas_loci', loci), ('mendelian', M), ('mendelian_records', mrec), ('cteph_targets', targets),
                     ('literature', lit), ('gene_condition', GCs), ('nodes', N), ('edges', E)]:
        write_tsv(df, out / f'{name}.tsv')
    for f in ('studies', 'findings', 'methods', 'screening'):
        shutil.copy(X / f'{f}.tsv', out / f'cteph_{f}.tsv')
    (out / 'graph.json').write_text(json.dumps({'nodes': N.to_dict('records'), 'edges': E.to_dict('records')},
                                               ensure_ascii=False, default=str))
    stats = {'release': meta['data_version'], 'api_version': meta['api_version'], 'fetch_date': meta['fetch_date'],
             'fetched_at': meta['fetched_at'], 'conditions': len(cond),
             'conditions_by_source': {s: int(cond.sources.str.contains(s).sum())
                                      for s in ('clinical', 'ontology', 'colocalisation')},
             'relations': len(rel), 'colocalisations': len(co), 'gwas_studies': len(studies),
             'credible_sets': len(loci), 'eas_studies': int(studies.eas.sum()), 'mendelian_links': len(M),
             'mendelian_records': len(mrec), 'cteph_targets': len(targets), 'literature_rows': len(lit),
             'literature_pmids': int(lit.pmid.nunique()), 'cteph_genetic_records': len(csc),
             'cteph_genetic_studies': len(cst), 'cteph_findings': len(cfi), 'methods_fields': len(cmf),
             'gene_condition': len(GCs), 'graph_nodes': len(N), 'graph_edges': len(E),
             'extract': json.loads((X / 'extract.json').read_text())}
    (out / 'tidy.json').write_text(json.dumps(stats, indent=1))
    log(f'{len(cond)} conditions, {len(studies)} GWAS studies, {len(loci)} credible sets, {len(M)} Mendelian links, '
        f'{len(cst)} CTEPH genetic studies; graph {len(N)} nodes / {len(E)} edges')


if __name__ == '__main__':
    main()
