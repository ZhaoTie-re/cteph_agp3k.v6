#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The knowledge base of one snapshot (Markdown + JSON), in the order
#           of the workbook's questions. Every section is built from the tables;
#           the generated summaries (summarise.py) are shown only when grounded,
#           always marked "Generated summary" and followed by the source records
#           they rest on.
#             1  Scope (and the questions)
#             2  Q1 Conditions related to CTEPH, and why (by category, with the
#                evidence)
#             3  Q2 Genes shared by CTEPH and the related conditions
#             4  Q3 Genetics of CTEPH (GWAS loci and their colocalisations; the
#                literature by gene theme; curated facts)
#             5  Q4 Methods of the genetic studies
#             6  Q5 Genetics of the related conditions (GWAS, curated genes)
#             7  Q6 Mechanisms (text-mined CTEPH targets by Reactome theme)
#             8  Q7 East Asian evidence
#             9  Changes since the previous snapshot
#             10 Caveats
# Inputs  : the tidy tables, summaries, changes, the seed and its PubMed
#           records, the gene themes
# Outputs : cteph_knowledge_graph.md, cteph_knowledge_graph.json
# ---------------------------------------------------------------------------
import argparse
import json
import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import (NO_GENE, ancestry_label, assign_themes, gene_tallies, is_east_asian, is_japanese,  # noqa: E402
                       load_json, load_themes, log, num, p_unicode, read_tsv, study_population, tally_text,
                       theme_order, variant_label)

CORE = 'EFO_0022954'
CATEGORY_ORDER = ['Precursor / spectrum', 'Risk factor', 'Comorbidity', 'Intermediate trait', 'PH group',
                  'Shared genetic signal']
QUESTIONS = ['Which conditions are related to CTEPH, and on what evidence?',
             'Which genes do CTEPH and its related conditions share?',
             'What is known about the genetics of CTEPH?',
             'How were the genetic studies of CTEPH done?',
             'What are the genetics of each related condition?',
             'Which mechanisms does the CTEPH literature describe?',
             'What evidence comes from East Asian and Japanese populations?']
CAVEATS = [
    'Relations to CTEPH are evidence of association (clinical studies, ontology, shared genetic signal), not of '
    'causation. Many classic risk factors were found against non-thromboembolic PH controls; against PE survivors '
    'most were not enriched (Dodson 2022, PMID 35000431).',
    'Colocalisation here is eCAVIAR on PICS fine-mapping (CLPP); a trait sharing only the ABO locus mostly reflects '
    'the pleiotropy of ABO, so single-locus relations are weaker than multi-locus ones.',
    'Locus-to-gene (L2G) scores are predictions of the causal gene, not proof.',
    'Curated gene evidence comes from curated resources (ClinGen, PanelApp, Orphanet, ClinVar, UniProt, '
    'Gene2Phenotype, gene burden) for the condition or its descendants in the ontology; it is mostly rare-variant, '
    'Mendelian evidence.',
    'The literature studies were screened and extracted by a language model from PubMed abstracts and open-access '
    'Methods sections; every gene, variant, number and quote was checked against the source text, and items that '
    'failed were dropped. Studies behind paywalls contribute their abstract only.',
    'Gene themes of the CTEPH literature are an editorial grouping; "associated in two or more studies" counts '
    'studies, not independent samples.',
    'Generated summaries are syntheses of the listed source records; they cite only those records, but they are not '
    'a systematic review.',
    'Text-mined target links can map an abbreviation to the wrong gene; the matched words are listed in the workbook.',
]


def md_table(cols, rows):
    esc = lambda x: str(x).replace('|', '/').replace('\n', ' ')  # noqa: E731
    return ['| ' + ' | '.join(cols) + ' |', '|' + '|'.join('---' for _ in cols) + '|'] + \
        ['| ' + ' | '.join(esc(x) for x in r) + ' |' for r in rows]


def counts(series):
    return ', '.join(f'{k} {v}' for k, v in series.value_counts().items())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tables', required=True)
    ap.add_argument('--summaries', required=True)
    ap.add_argument('--changes', required=True)
    ap.add_argument('--changes-json', required=True)
    ap.add_argument('--seed', required=True)
    ap.add_argument('--seed-pubmed', required=True)
    ap.add_argument('--themes', required=True, help='gene themes of the CTEPH literature (YAML)')
    ap.add_argument('--l2g-min', type=float, default=0.5)
    ap.add_argument('--top-genes', type=int, default=8, help='GWAS genes listed per related condition')
    ap.add_argument('--max-changes', type=int, default=40)
    ap.add_argument('--out-md', default='cteph_knowledge_graph.md')
    ap.add_argument('--out-json', default='cteph_knowledge_graph.json')
    a = ap.parse_args()
    T = Path(a.tables)
    st = json.loads((T / 'tidy.json').read_text())
    cond, rel, co = read_tsv(T / 'conditions.tsv'), read_tsv(T / 'relations.tsv'), read_tsv(T / 'colocalisation.tsv')
    studies, loci, men = read_tsv(T / 'gwas_studies.tsv'), read_tsv(T / 'gwas_loci.tsv'), read_tsv(T / 'mendelian.tsv')
    gc = read_tsv(T / 'gene_condition.tsv')
    cst, cfi = read_tsv(T / 'cteph_studies.tsv'), read_tsv(T / 'cteph_findings.tsv')
    cmf = read_tsv(T / 'cteph_methods.tsv')
    summ = read_tsv(a.summaries)
    changes, chj = read_tsv(a.changes), json.loads(Path(a.changes_json).read_text())
    seed = yaml.safe_load(Path(a.seed).read_text())
    themes = load_themes(a.themes)
    pub = {r['pmid']: r for r in load_json(a.seed_pubmed)['records']}
    cite = lambda p: f'{pub[p]["first_author"].split()[0]} {pub[p]["year"]}, PMID {p}' if p in pub else f'PMID {p}'  # noqa: E731
    names = dict(zip(cond.node_id, cond.name))
    cst['population'] = [study_population(r.country, r.ancestry) for r in cst.itertuples()]
    cst['jp'] = [is_japanese(r) for r in cst.itertuples()]
    cst['eas'] = [is_east_asian(r) for r in cst.itertuples()]
    K = {'release': st['release'], 'fetch_date': st['fetch_date'], 'questions': QUESTIONS,
         'counts': {k: v for k, v in st.items() if isinstance(v, (int, dict))}}
    md = [f'# CTEPH knowledge graph (Open Targets {st["release"]}, PubMed; fetched {st["fetch_date"]})', '',
          'Generated by `integration/cteph_knowledge_graph`; do not edit by hand. Sections are built from the '
          'snapshot tables; paragraphs headed **Generated summary** were written by a language model from the '
          'records listed with them and checked against those records.', '']

    # 1 scope ───────────────────────────────────────────────────────────────
    bs = st['conditions_by_source']
    md += ['## 1. Scope', '',
           f'- {st["conditions"]} conditions related to CTEPH (EFO_0022954): {bs["clinical"]} by clinical evidence '
           f'(curated, PubMed-verified), {bs["ontology"]} by the disease ontology, {bs["colocalisation"]} by a GWAS '
           'signal shared with a CTEPH locus (a condition can have several).',
           f'- Genetics: {st["gwas_studies"]} GWAS studies ({st["eas_studies"]} East Asian) with {st["credible_sets"]:,} '
           f'credible sets; {st["mendelian_links"]} condition-gene links from curated resources.',
           f'- CTEPH literature: {st["cteph_genetic_records"]} PubMed records screened, {st["cteph_genetic_studies"]} '
           f'genetic studies of CTEPH extracted ({st["cteph_findings"]} findings, {st["methods_fields"]} methods '
           'fields).',
           f'- Graph: {st["graph_nodes"]} nodes, {st["graph_edges"]} edges (`graph.json`).', '',
           'Sections 2 to 8 answer the questions of the workbook, in its order (workbook sheets 1 to 7):', '']
    md += [f'{k}. {q}' for k, q in enumerate(QUESTIONS, 1)] + ['']

    # 2 related conditions ──────────────────────────────────────────────────
    md += ['## 2. Conditions related to CTEPH, and why', '', f'*{QUESTIONS[0]}*', '']
    K['conditions'] = []
    for cat in CATEGORY_ORDER:
        sub = cond[cond.category == cat]
        if not len(sub):
            continue
        md += [f'### {cat}', '']
        if cat == 'Shared genetic signal':
            cr = rel[(rel.evidence_type == 'colocalisation') & rel.node_id.isin(sub.node_id)].copy()
            cr['n_loci'] = num(cr.n_loci)
            multi, single = cr[cr.n_loci >= 2], cr[cr.n_loci < 2]
            md += [f'Traits whose GWAS signal colocalises with CTEPH loci and that have no clinical or ontology '
                   f'relation: {len(sub)}. Shared with two or more CTEPH loci ({len(multi)}):', '']
            md += md_table(['Trait', 'CTEPH loci', 'Colocalising GWAS', 'Max CLPP'],
                           [[names[r.node_id], r.loci, r.n_studies, r.strength]
                            for r in multi.sort_values('n_loci', ascending=False).itertuples()])
            md += ['', f'Shared with one CTEPH locus only ({len(single)}; mostly the pleiotropic ABO locus): ' +
                   ', '.join(f'{names[r.node_id]} ({r.loci})' for r in single.sort_values('loci').itertuples())
                   + '.', '']
            continue
        rows = []
        for c in sub.itertuples():
            ev = rel[rel.node_id == c.node_id]
            cl = ev[ev.evidence_type == 'clinical']
            other = ev[ev.evidence_type != 'clinical']
            items = [f'{r.detail} ({cite(r.pmid)}{"; East Asian" if r.east_asian == "True" else ""})'
                     for r in cl.itertuples()]
            items += [r.detail for r in other.itertuples()]
            K['conditions'].append({'condition': c.name, 'id': c.node_id, 'category': cat, 'sources': c.sources,
                                    'evidence': items})
            rows.append([c.name, c.sources.replace(' ', ', '), '<br>'.join(items)])
        md += md_table(['Condition', 'Evidence types', 'Evidence'], rows) + ['']

    # 3 shared genes ────────────────────────────────────────────────────────
    gc['kind'] = gc.evidence.replace({'Mendelian': 'curated', 'Literature': 'literature'})
    gcs = gc.groupby(['ensembl_id', 'symbol']).agg(
        n=('node_id', lambda x: len(set(x) - {CORE})), cteph=('node_id', lambda x: CORE in set(x)),
        conds=('condition', lambda x: '; '.join(sorted(set(x) - {'chronic thromboembolic pulmonary hypertension'}))),
        kinds=('kind', lambda x: '/'.join(sorted(set(x))))).reset_index()
    ct = gcs[gcs.cteph & (gcs.n >= 1)].sort_values(['n', 'symbol'], ascending=[False, True])
    other = gcs[~gcs.cteph & (gcs.n >= 3)].sort_values(['n', 'symbol'], ascending=[False, True])
    K['shared_genes'] = ct.to_dict('records')
    md += ['## 3. Genes shared by CTEPH and the related conditions', '', f'*{QUESTIONS[1]}*', '',
           f'{int(gcs.cteph.sum())} genes have genetic evidence in CTEPH (its GWAS, top L2G ≥ {a.l2g_min}, or a '
           f'genetic study in the literature); {len(ct)} of them also in at least one related condition (GWAS or '
           'curated evidence):', '']
    md += md_table(['Gene', 'Related conditions', 'Evidence kinds', 'Which'],
                   [[r.symbol, r.n, r.kinds, r.conds[:200]] for r in ct.itertuples()]) + ['']
    md += [f'Without CTEPH evidence, {len(other)} genes have genetic evidence in three or more related conditions; '
           'the most widely shared: ' + ', '.join(f'{r.symbol} ({r.n})' for r in other.head(15).itertuples()) + '.',
           '']

    # 4 genetics of CTEPH ───────────────────────────────────────────────────
    md += ['## 4. Genetics of CTEPH', '', f'*{QUESTIONS[2]}*', '', '### GWAS (Open Targets)', '']
    cl = loci[loci.node_id == CORE]
    for sid, s in cl.groupby('study_id'):
        r0 = studies[(studies.node_id == CORE) & (studies.study_id == sid)].iloc[0]
        nc, nn = num(pd.Series([r0.n_cases, r0.n_controls]))
        md += [f'{r0.first_author} {r0.year}, {r0.journal} (PMID {r0.pmid}; {sid}; {nc:,.0f} cases / '
               f'{nn:,.0f} controls; {ancestry_label(r0.ancestry)}):', '']
        rows = []
        for r in s.itertuples():
            cc = co[co.cteph_locus == r.rsid]
            traits = cc.condition.replace('', pd.NA).dropna().value_counts()
            rows.append([r.rsid, variant_label(r.variant_id), p_unicode(r.p_mantissa, r.p_exponent), r.finemapping,
                         ', '.join(f'{getattr(r, f"l2g{k}_symbol")} {float(getattr(r, f"l2g{k}_score")):.2f}'
                                   for k in (1, 2, 3) if getattr(r, f'l2g{k}_symbol')),
                         ', '.join(f'{t} ({n})' for t, n in traits.head(8).items())])
        md += md_table(['Lead variant', 'Position', 'P', 'Fine-mapping', 'Locus-to-gene', 'Colocalising traits '
                        '(GWAS)'], rows) + ['']
    fx = assign_themes(cfi, themes).merge(cst[['pmid', 'year', 'first_author', 'eas']], on='pmid', how='left')
    gt = gene_tallies(fx)
    rep = gt[(gt.associated >= 2) & (gt.theme != NO_GENE)].sort_values('associated', ascending=False)
    md += ['### Genetic studies in the literature (PubMed), by gene theme', '',
           f'{len(cst)} studies, {len(cfi)} findings on {gt.label.nunique()} genes or gene sets. Associated in two or '
           'more studies: ' + ', '.join(f'{r.label} ({r.associated})' for r in rep.itertuples()) + '. Gene themes '
           'are an editorial grouping (config/gene_themes.yaml).', '']
    sf = summ[(summ.kind == 'findings') & (summ.grounded == 'True')]
    if len(sf):
        md += ['**Generated summary (findings).**', '', sf.iloc[0].summary, '']
    desc = {t['name']: t['description'] for t in themes}
    K['themes'] = []
    for th in theme_order(themes):
        g = gt[gt.theme == th].sort_values(['associated', 'n_studies', 'label'], ascending=[False, False, True])
        if not len(g):
            continue
        f = fx[fx.theme == th]
        md += [f'#### {th}', '', (desc.get(th, '') + ' ').lstrip() + f'{len(g)} genes, {f.pmid.nunique()} studies, '
               f'{len(f)} findings.', '']
        rows = []
        for r in g.itertuples():
            ff = f[f.label == r.label].assign(_y=lambda d: num(d.year)).sort_values('_y')
            refs = '; '.join(dict.fromkeys(f'{x.first_author} {x.year} ({x.variant + ", " if x.variant else ""}'
                                           f'{x.result}; PMID {x.pmid})' for x in ff.itertuples()))
            rows.append([r.label, tally_text(r), refs])
        md += md_table(['Gene', 'Studies', 'Findings (variant, result)'], rows) + ['']
        K['themes'].append({'theme': th, 'genes': g.to_dict('records')})
    md += ['### Curated facts', ''] + [f'- {f["topic"]}: {f["statement"]} ({cite(str(f["pmid"]))})'
                                       for f in seed.get('cteph_facts') or []] + ['']
    K['cteph_genetics'] = {'gwas_loci': len(cl), 'literature_studies': len(cst), 'findings': len(cfi),
                           'facts': seed.get('cteph_facts') or []}

    # 5 methods ─────────────────────────────────────────────────────────────
    sm = summ[(summ.kind == 'methods') & (summ.grounded == 'True')]
    gw = studies.copy()
    gw['n_samples'] = num(gw.n_samples)
    md += ['## 5. Methods of the genetic studies', '', f'*{QUESTIONS[3]}*', '', '### CTEPH literature', '']
    if len(cst):
        md += ['- Designs (a study can have several): ' + counts(cst.design.str.split('; ').explode()) + '.',
               '- Populations: ' + counts(cst.population) + '.',
               f'- Cases per study: median {num(cst.n_cases).median():.0f}, range {num(cst.n_cases).min():.0f}–'
               f'{num(cst.n_cases).max():.0f} (where stated).',
               f'- Methods fields ({len(cmf)}: platform, sequencing, quality control, statistics, covariates, '
               f'multiple testing): from the open-access Methods section for {cst.methods_source.eq("full text").sum()} '
               f'studies, from the abstract for {cst.methods_source.eq("abstract").sum()} (no full text).', '']
        rows = [[s.year, s.first_author, s.design, s.population, s.n_cases or '—', s.n_controls or '—',
                 s.methods_source or '—', s.pmid] for s in cst.assign(_y=num(cst.year)).sort_values('_y').itertuples()]
        md += md_table(['Year', 'First author', 'Design', 'Population', 'Cases', 'Controls', 'Methods from', 'PMID'],
                       rows) + ['']
    if len(sm):
        md += ['**Generated summary (methods).**', '', sm.iloc[0].summary, '']
    md += ['### GWAS in Open Targets (CTEPH and the related conditions)', '',
           f'- {len(gw)} studies; sample size median {gw.n_samples.median():,.0f} (range {gw.n_samples.min():,.0f}–'
           f'{gw.n_samples.max():,.0f}); with summary statistics {gw.summary_statistics.eq("True").sum()}.',
           '- Fine-mapping of the credible sets: ' + counts(loci.finemapping) + '.',
           '- Largest ancestry group of each study: ' + counts(gw.ancestry.str.split('; ').str[0].str.split().str[0])
           + '.', '']

    # 6 genetics of related conditions ──────────────────────────────────────
    md += ['## 6. Genetics of the related conditions', '', f'*{QUESTIONS[4]}*', '']
    rows = []
    loci['l2g1_score_n'] = num(loci.l2g1_score)
    for c in cond[cond.genetics == 'full'].sort_values(['category', 'name']).itertuples():
        L = loci[(loci.node_id == c.node_id) & (loci.l2g1_score_n >= a.l2g_min)]
        top = L.groupby('l2g1_symbol').study_id.nunique().sort_values(ascending=False).head(a.top_genes)
        M = men[men.node_id == c.node_id] if len(men) else men
        mg = M.groupby('symbol').source.apply(lambda x: '/'.join(sorted(set(x)))) if len(M) else pd.Series(dtype=str)
        rows.append([c.name, c.category, f'{c.n_gwas_studies} ({c.n_eas_studies} EAS)', c.n_credible_sets,
                     ', '.join(f'{g} ({n})' for g, n in top.items()) or '—',
                     ', '.join(f'{g} [{s}]' for g, s in list(mg.items())[:12]) or '—'])
    md += md_table(['Condition', 'Category', 'GWAS studies', 'Credible sets', 'Top GWAS genes (studies)',
                    'Curated genes [source]'], rows) + ['']

    # 7 mechanisms ──────────────────────────────────────────────────────────
    md += ['## 7. Mechanisms (text-mined CTEPH literature, by Reactome theme)', '', f'*{QUESTIONS[5]}*', '']
    mech = summ[summ.kind == 'mechanism'].copy()
    mech['n'] = num(mech.n_targets)
    for r in mech.sort_values('n', ascending=False).itertuples():
        md += [f'### {r.title} ({int(r.n)} targets)', '']
        if r.grounded == 'True':
            md += ['**Generated summary.** ' + r.summary, '']
        else:
            md += [f'(Summary withheld: {r.problems}.)', '']
        md += [f'Source targets: {r.genes}. Source PMIDs: {r.pmids}.', '']
    K['mechanisms'] = summ[summ.kind == 'mechanism'].to_dict('records')

    # 8 East Asian evidence ─────────────────────────────────────────────────
    eas = studies[studies.eas == 'True']
    ea_rel = rel[rel.east_asian == 'True']
    ea_lit = cst[cst.eas]
    md += ['## 8. East Asian evidence', '', f'*{QUESTIONS[6]}*', '',
           f'- Clinical evidence from East Asian populations ({len(ea_rel)}): ' +
           '; '.join(f'{names[r.node_id]}: {r.detail} ({cite(r.pmid)})' for r in ea_rel.itertuples()) + '.',
           f'- Genetic studies of CTEPH in East Asian populations ({len(ea_lit)}, Japanese {int(ea_lit.jp.sum())}): ' +
           '; '.join(f'{r.first_author} {r.year} ({r.population}; {r.design}; {r.gene_symbols or r.genes or "no gene"}'
                     f'; PMID {r.pmid})' for r in ea_lit.sort_values('year').itertuples()) + '.',
           f'- GWAS of the related conditions with a mostly East Asian sample ({len(eas)}): ' +
           '; '.join(f'{r.condition}: {r.study_id} ({r.first_author or r.project} {r.year}, n = {r.n_samples})'
                     for r in eas.itertuples()) + '.', '']

    # 9 changes ─────────────────────────────────────────────────────────────
    md += ['## 9. Changes since the previous snapshot', '']
    if chj['previous'] is None:
        md += ['First snapshot: nothing to compare with.', '']
    else:
        md += [f'Compared with {chj["previous"]} (release {chj["previous_release"]}): ' +
               ('; '.join(f'{k.replace("|", " ")} {v}' for k, v in chj['counts'].items()) or 'no change') + '.', '']
        if len(changes):
            md += md_table(['Category', 'Change', 'Item', 'Detail'],
                           [[r.category, r.change, r.item, r.detail] for r in changes.head(a.max_changes).itertuples()])
            md += ['']
    K['changes'] = chj

    md += ['## 10. Caveats', ''] + [f'- {c}' for c in CAVEATS] + ['']
    K['caveats'] = CAVEATS
    Path(a.out_md).write_text('\n'.join(md))
    Path(a.out_json).write_text(json.dumps(K, indent=1, ensure_ascii=False, default=str))
    log(f'knowledge: {len(md)} lines')


if __name__ == '__main__':
    main()
