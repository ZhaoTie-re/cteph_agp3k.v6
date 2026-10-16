#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The conditions of the knowledge graph and the evidence that relates
#           each to CTEPH. A condition enters by at least one of:
#             clinical       the curated seed (config/related_conditions.yaml): a
#                            guideline, review or study whose PubMed abstract states
#                            the relation; every quote is checked against the
#                            abstract fetched for this snapshot
#             ontology       the disease ontology of Open Targets: the parent of CTEPH,
#                            its siblings, and the children of pulmonary arterial
#                            hypertension
#             colocalisation a GWAS credible set of the condition colocalises with a
#                            CTEPH GWAS locus (eCAVIAR CLPP >= --clpp-min, or coloc
#                            H4 >= --h4-min); graded by the CTEPH loci it shares
#           Every seed identifier is checked against Open Targets. Each condition
#           gets a genetics mode: 'full' (its GWAS and curated / Mendelian
#           evidence are fetched), 'coloc_only' (only the colocalising credible
#           sets, from the CTEPH fetch) or 'none' (no Open Targets record: e.g.
#           procedures such as splenectomy).
# Inputs  : the seed, the seed's PubMed records, the CTEPH GWAS with colocalisations
# Outputs : conditions.tsv (one row per condition), relations.tsv (one row per
#           piece of evidence), colocalisation.tsv (one row per colocalising pair),
#           fetch_list.tsv (the conditions whose genetics are fetched), seed_check.tsv
# ---------------------------------------------------------------------------
import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from llm import in_text  # noqa: E402
from ot_common import gql, load_json, log, write_tsv  # noqa: E402

DISEASES = '''query($ids: [String!]!) { diseases(efoIds: $ids) { id name description parents { id name }
  children { id name } } }'''
CORE = 'EFO_0022954'


def records(api, ids):
    out = {}
    ids = sorted(set(ids))
    for i in range(0, len(ids), 100):
        for d in gql(api, DISEASES, {'ids': ids[i:i + 100]})['diseases']:
            out[d['id']] = d
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--seed', required=True, help='config/related_conditions.yaml')
    ap.add_argument('--seed-pubmed', required=True, help="the seed PMIDs' PubMed records (pubmed.py --mode ids)")
    ap.add_argument('--cteph-gwas', required=True, help='gwas.json.gz of the CTEPH fetch (with colocalisations)')
    ap.add_argument('--api', required=True)
    ap.add_argument('--clpp-min', type=float, default=0.01)
    ap.add_argument('--h4-min', type=float, default=0.8)
    ap.add_argument('--out-dir', default='.')
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    seed = yaml.safe_load(Path(a.seed).read_text())
    pub = {r['pmid']: r for r in load_json(a.seed_pubmed)['records']}

    # ── clinical evidence: every quote must be in its abstract ───────────────
    conds, rel, check, alias = {}, [], [], {}
    for c in seed['conditions']:
        ids = c.get('ot_ids') or []
        node = ids[0] if ids else f'custom:{c["key"]}'
        alias.update({i: node for i in ids})
        conds[node] = {'node_id': node, 'name': c['name'], 'category': c['category'], 'ot_ids': ' '.join(ids),
                       'genetics': c.get('genetics', 'full' if ids else 'none'),
                       'mendelian': bool(c.get('mendelian', bool(ids))), 'genes': ' '.join(c.get('genes') or []),
                       'seed_key': c['key'], 'sources': {'clinical'}}
        for e in (c.get('evidence') or []) + [dict(x, east_asian=True) for x in (c.get('east_asian') or [])]:
            r = pub.get(str(e['pmid']))
            ok = r is not None and in_text(e['quote'], r['title'] + ' ' + r['abstract'])
            check.append({'node_id': node, 'pmid': str(e['pmid']), 'found': r is not None, 'quote_in_abstract': ok})
            if not ok:
                continue
            rel.append({'node_id': node, 'evidence_type': 'clinical', 'relation': e.get('direction', ''),
                        'detail': e['statement'], 'quote': e['quote'], 'design': e.get('design', ''),
                        'population': e.get('population', ''), 'east_asian': bool(e.get('east_asian')),
                        'pmid': str(e['pmid']), 'citation': f'{r["first_author"]} {r["year"]}, {r["journal"]}',
                        'strength': '', 'loci': ''})
    chk = pd.DataFrame(check)
    bad = chk[~chk.quote_in_abstract]
    assert not len(bad), f'seed quotes not found in their PubMed abstracts:\n{bad.to_string()}'
    no_ev = [n for n in conds if not any(r['node_id'] == n for r in rel)]
    assert not no_ev, f'seed conditions without verified evidence: {no_ev}'

    # ── ontology ─────────────────────────────────────────────────────────────
    recs = records(a.api, [CORE] + [i for c in seed['conditions'] for i in (c.get('ot_ids') or [])])
    missing = [i for c in seed['conditions'] for i in (c.get('ot_ids') or []) if i not in recs]
    assert not missing, f'seed IDs not in Open Targets: {missing}'
    core = recs[CORE]
    onto = []
    for p in core['parents']:
        onto.append((p['id'], 'parent of CTEPH'))
        for s in records(a.api, [p['id']])[p['id']]['children']:
            if s['id'] != CORE:
                onto.append((s['id'], f'sibling of CTEPH (child of {p["name"]})'))
    pah = 'MONDO_0015924'
    for ch in records(a.api, [pah])[pah]['children']:
        onto.append((ch['id'], 'child of pulmonary arterial hypertension'))
    orec = records(a.api, [i for i, _ in onto])
    for i0, why in onto:
        i = alias.get(i0, i0)
        if i not in conds:
            conds[i] = {'node_id': i, 'name': orec[i0]['name'], 'category': 'PH group', 'ot_ids': i0,
                        'genetics': 'full', 'mendelian': True, 'genes': '', 'seed_key': '', 'sources': set()}
        conds[i]['sources'].add('ontology')
        rel.append({'node_id': i, 'evidence_type': 'ontology', 'relation': why, 'detail': f'Open Targets disease '
                    f'ontology: {orec[i0]["name"]} is the {why}', 'quote': '', 'design': '', 'population': '',
                    'east_asian': False, 'pmid': '', 'citation': 'EFO / MONDO via Open Targets', 'strength': '', 'loci': ''})

    # ── colocalisation with the CTEPH loci ───────────────────────────────────
    g = load_json(a.cteph_gwas)
    co = []
    for s in g['studies']:
        for c in s['credible_sets']:
            locus = (c['variant'] or {}).get('rsIds', [''])[0] or c['variant']['id']
            gene = ((c.get('l2GPredictions') or {}).get('rows') or [{}])[0].get('target', {}).get('approvedSymbol', '')
            for x in c.get('colocalisation') or []:
                if not ((x.get('clpp') or 0) >= a.clpp_min or (x.get('h4') or 0) >= a.h4_min):
                    continue
                o = x['otherStudyLocus']
                st = o['study']
                for d in st.get('diseases') or [{'id': '', 'name': st.get('traitFromSource') or ''}]:
                    co.append({'cteph_locus': locus, 'cteph_gene': gene, 'cteph_study_locus_id': c['studyLocusId'],
                               'node_id': d['id'], 'condition': d['name'], 'study_id': st['id'],
                               'study_type': x['rightStudyType'], 'trait': st.get('traitFromSource') or '',
                               'first_author': st.get('publicationFirstAuthor') or '',
                               'pub_date': st.get('publicationDate') or '', 'pmid': st.get('pubmedId') or '',
                               'project': st.get('projectId') or '', 'n_samples': st.get('nSamples'),
                               'ancestry': '; '.join(f'{p["ldPopulation"]} {100 * p["relativeSampleSize"]:.0f}%'
                                                     for p in st.get('ldPopulationStructure') or []),
                               'other_study_locus_id': o['studyLocusId'],
                               'other_lead': ((o.get('variant') or {}).get('rsIds') or [''])[0]
                               or (o.get('variant') or {}).get('id', ''),
                               'p_mantissa': o.get('pValueMantissa'), 'p_exponent': o.get('pValueExponent'),
                               'method': x['colocalisationMethod'], 'clpp': x.get('clpp'), 'h4': x.get('h4'),
                               'n_coloc_variants': x.get('numberColocalisingVariants'),
                               'beta_ratio_sign': x.get('betaRatioSignAverage')})
    co = pd.DataFrame(co)
    agg = defaultdict(lambda: {'loci': set(), 'studies': set(), 'clpp': 0.0, 'h4': 0.0, 'name': ''})
    for r in co[co.node_id.ne('')].itertuples():
        a_ = agg[alias.get(r.node_id, r.node_id)]
        a_['loci'].add(r.cteph_locus)
        a_['studies'].add(r.study_id)
        a_['clpp'] = max(a_['clpp'], r.clpp or 0)
        a_['h4'] = max(a_['h4'], r.h4 if pd.notna(r.h4) else 0)
        a_['name'] = r.condition
    for i, v in agg.items():
        if i == CORE:
            continue
        if i not in conds:
            conds[i] = {'node_id': i, 'name': v['name'], 'category': 'Shared genetic signal', 'ot_ids': i,
                        'genetics': 'coloc_only', 'mendelian': False, 'genes': '', 'seed_key': '', 'sources': set()}
        conds[i]['sources'].add('colocalisation')
        grade = 'multi-locus' if len(v['loci']) >= 2 else 'single-locus'
        rel.append({'node_id': i, 'evidence_type': 'colocalisation', 'relation': f'shared genetic signal ({grade})',
                    'detail': f'{len(v["studies"])} GWAS colocalise with CTEPH loci {", ".join(sorted(v["loci"]))}',
                    'quote': '', 'design': '', 'population': '', 'east_asian': False, 'pmid': '',
                    'citation': 'Open Targets colocalisation', 'strength': f'{v["clpp"]:.3f}',
                    'loci': ' '.join(sorted(v['loci'])), 'n_loci': len(v['loci']), 'n_studies': len(v['studies'])})

    C = pd.DataFrame(conds.values())
    C['sources'] = C.sources.map(lambda s: ' '.join(sorted(s)))
    R = pd.DataFrame(rel)
    for c in ('n_loci', 'n_studies'):
        if c not in R:
            R[c] = None
    C['n_clinical'] = C.node_id.map(R[R.evidence_type == 'clinical'].groupby('node_id').size()).fillna(0).astype(int)
    C['coloc_loci'] = C.node_id.map(R[R.evidence_type == 'colocalisation'].set_index('node_id').n_loci).fillna(0)
    C['coloc_loci'] = C.coloc_loci.astype(int)
    C = C.sort_values(['category', 'name'])
    fetch = C[C.genetics.eq('full') & C.ot_ids.ne('')][['node_id', 'mendelian']]
    write_tsv(C, out / 'conditions.tsv')
    write_tsv(R, out / 'relations.tsv')
    write_tsv(co, out / 'colocalisation.tsv')
    write_tsv(chk, out / 'seed_check.tsv')
    fetch.assign(mendelian=fetch.mendelian.map(lambda x: 'true' if x else 'false')).to_csv(
        out / 'fetch_list.tsv', sep='\t', index=False, header=False)
    log(f'{len(C)} conditions ({C.sources.str.contains("clinical").sum()} clinical, '
        f'{C.sources.str.contains("ontology").sum()} ontology, {C.sources.str.contains("colocalisation").sum()} '
        f'colocalisation); {len(R)} evidence rows; {len(fetch)} genetics fetches')


if __name__ == '__main__':
    main()
