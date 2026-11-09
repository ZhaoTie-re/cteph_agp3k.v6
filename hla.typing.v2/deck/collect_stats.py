#!/usr/bin/env python3
"""Every number the deck states, read from the component's result tables -> stats.json."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'scripts'))
import config as C  # noqa: E402
import nomenclature as N  # noqa: E402


def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else float(x)
    return x


# cohorts and floors come from the pipeline's config (step 06 applies them); the deck only reads the result
COHORTS, MODEL_INPUTS = C.COHORTS, C.MODEL_INPUTS
MIN_MAF, MIN_CALL_RATE = C.MIN_MAF, C.MIN_CALL_RATE
ASSOC = C.D_ASSOC


def freq_outlier(R, gene):
    """The point that stands out in one gene's panel, and why it stands out. Distance is measured the way the
    figure draws it (square-root axes), so the allele named here is the one the eye picks out. The explanation
    is assembled from the dictionary itself: the name the release gives the allele, the release that changed
    that name, whether it has a P group, whether HLA-HD's frequency table carries it, the allele it cannot be
    told apart from in the protein alignment, and the 2-field group totals on both sides."""
    al = pd.read_csv(R / '03.panels' / 'allele_concordance.tsv', sep='\t')
    al = al[(al.panel == 'ToMMo') & (al.gene == gene)].copy()
    al['vis'] = (np.sqrt(al.freq_ours) - np.sqrt(al.freq_panel)).abs()
    o = al.sort_values('vis', ascending=False).iloc[0]
    out = {'gene': gene, 'allele': o.allele, 'freq_ours': round(float(o.freq_ours), 5),
           'freq_panel': round(float(o.freq_panel), 5), 'count_ours': int(o.count_ours),
           'count_panel': int(o.count_panel), 'runner_up': round(float(al.vis.sort_values().iloc[-2]), 4)}

    # the name the pinned release gives it, and the release in which that name last changed
    hist = pd.read_csv(C.IMGT_HISTORY, comment='#', dtype=str).set_index('HLA_ID')
    stem = o.allele.split('*')[1]
    row = next((hist.loc[i] for i in hist.index
                if str(hist.loc[i].iloc[0]).startswith(f'{gene}*{stem}')
                and str(hist.loc[i].iloc[0]).rstrip('QNLSCA') == f'{gene}*{stem}'), None)
    if row is not None:
        out['name_now'] = row.iloc[0]
        changed = next((rel for rel, v in row.items() if v != row.iloc[0]), None)
        out['name_before'] = row[changed] if changed else None
        out['renamed_in'] = int(hist.columns[list(hist.columns).index(changed) - 1]) if changed else None

    name = out.get('name_now', o.allele)
    # P group: an empty third field means IPD puts the allele in no group, so the P-group step cannot merge it
    pg = [ln.split(';') for ln in C.P_GROUPS.read_text().splitlines() if ln.startswith(f'{gene}*')]
    out['p_group'] = next((f[2] for f in pg if name.split('*')[1] in f[1].split('/') and f[2]), None)
    # HLA-HD's frequency table decides which allele is reported when the reads cannot separate a pair
    freq = (C.HLAHD_FREQ / f'{gene}_count.txt')
    out['in_hlahd_freq'] = bool(freq.exists() and f'{gene}*{stem}' in
                               [w.split('\t')[0] for w in freq.read_text().split('\n')])
    # the allele it cannot be told apart from: identical wherever both protein sequences are determined
    ref = pd.read_csv(R / '05.exports' / 'residue_reference' / f'{gene}.tsv', sep='\t', index_col=0)
    if name in ref.index:
        a = ref.loc[name]
        bad = {'*', '.', '-', 'X'}
        best = None
        for other in ref.index:
            two = ':'.join(str(other).split(':')[:2])
            if two == name or two in (out.get('name_before'), name.rstrip('QNLSCA')):
                continue
            b = ref.loc[other]
            both = a.notna() & b.notna() & ~a.isin(bad) & ~b.isin(bad)
            n = int(both.sum())
            d = int((a[both] != b[both]).sum())
            if n and (best is None or (d, -n) < (best[1], -best[2])):
                best = (two, d, n)
        if best and best[1] == 0:
            out['identical_to'] = {'allele': best[0], 'positions': best[2]}
    # the 2-field group both sides agree on, even though they split it differently
    grp = al[al.allele.str.startswith(f'{gene}*{stem.split(":")[0]}:')]
    out['group'] = {'name': f'{gene}*{stem.split(":")[0]}', 'ours': round(float(grp.freq_ours.sum()), 4),
                    'panel': round(float(grp.freq_panel.sum()), 4)}
    out['neighbours'] = [{'allele': r.allele, 'ours': round(float(r.freq_ours), 4),
                          'panel': round(float(r.freq_panel), 4)}
                         for _, r in grp.sort_values('freq_ours', ascending=False).head(2).iterrows()]
    return out


def panel_coverage(R, S):
    """Per compared gene: the distinct 2-field alleles we call, split by what ToMMo does with them, and the same
    split counted in copies. An allele first named after the panel's newest allele cannot be in the panel whatever
    its true frequency, so it is kept apart from the alleles the panel could have carried and does not."""
    vint = S['tommo_vintage']
    a = pd.read_csv(R / '03.panels' / 'allele_concordance.tsv', sep='\t')
    a = a[(a.panel == 'ToMMo') & (a.count_ours > 0)].copy()
    rel = a.allele.map(N.release_of)
    a['state'] = np.where(a.count_panel > 0, 'in_panel',
                          np.where(rel.isna() | (rel > vint), 'too_new', 'absent'))
    out = {}
    for g in [k for k, v in S['tommo'].items() if v['frequency_comparable'] == 1]:
        d = a[a.gene == g]
        alleles = d.groupby('state').allele.count()
        copies = d.groupby('state').count_ours.sum()
        out[g] = {'alleles': {k: int(alleles.get(k, 0)) for k in ('in_panel', 'absent', 'too_new')},
                  'copies': {k: int(copies.get(k, 0)) for k in ('in_panel', 'absent', 'too_new')}}
    # the 2-field names (null alleles left out, as in the comparison) that the controls carry, which the
    # comparison merges into the P groups above; the association matrices keep these 2-field names
    ex = pd.read_csv(R / '05.exports' / 'allele_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    ph = pd.read_csv(MODEL_INPUTS / 'full_mainland' / 'phenotype' / 'pheno.tsv', sep='\t', dtype={'IID': str})
    seen = ex.loc[ph.loc[ph.PHENO1 == 0, 'IID']].sum()
    for g in out:
        out[g]['two_field'] = int(((seen > 0) & (seen.index.str.split('*').str[0] == g)
                                   & ~seen.index.str.endswith('N')).sum())
    return out


def assoc_basis(R):
    """What the association reads, per cohort — taken from step 06's own record (06.assoc_inputs/markers.tsv,
    summary.json), not recomputed, so every marker count on the slides is the pipeline's. Distributions are over
    the columns of the 16 genes that vary in the cohort; the MAF classes over those that also pass the call-rate
    floor, the order the pipeline applies them in."""
    mk = pd.read_csv(ASSOC / 'markers.tsv', sep='\t', dtype={'position': str, 'residue': str})
    summ = json.loads((ASSOC / 'summary.json').read_text())
    ad_cols = pd.read_csv(R / '05.exports' / 'allele_dosage.tsv', sep='\t', nrows=0).columns[1:]
    rd_cols = pd.read_csv(R / '05.exports' / 'residue_dosage.tsv', sep='\t', nrows=0).columns[1:]
    sites = pd.read_csv(R / '05.exports' / 'residue_sites.tsv', sep='\t')
    genes = list(C.ASSOC_GENES)
    out = {'genes': genes, 'polymorphic_positions': int(sites[sites.gene.isin(genes)].polymorphic.sum()),
           'allele_columns': len(ad_cols), 'residue_columns': len(rd_cols), 'cohorts': {},
           'exported_per_gene': {g: {'alleles': int((ad_cols.str.split('*').str[0] == g).sum()),
                                     'residues': int((rd_cols.str.split(':').str[0] == g).sum())} for g in genes}}
    edges = [0, 0.001, 0.005, 0.01, 0.05, 0.10, 0.25, 1.0]
    labels = ['<0.1%', '0.1-0.5%', '0.5-1%', '1-5%', '5-10%', '10-25%', '>=25%']
    for co in COHORTS:
        sc = summ['cohorts'][co]
        m = mk[mk.cohort == co]
        tested = m[m.status == 'tested']
        dist = {}
        for kind in ('allele', 'residue'):
            k = m[(m['class'] == kind) & m.gene.isin(genes)]
            live = k[k.status != 'not_varying']                       # the columns a floor can act on
            # the call rate the rule uses: the lower of the two groups' (cases, controls)
            cr = live[['call_rate_cases', 'call_rate_controls']].min(axis=1)
            ok = live[cr >= MIN_CALL_RATE]
            floors = {'MAC>=20': ok.mac >= 20, 'MAF>=0.5%': ok.maf >= MIN_MAF, 'MAF>=1%': ok.maf >= 0.01}
            dist[kind] = {
                'columns': len(k), 'mac0': int((k.status == 'not_varying').sum()),
                'live': len(live), 'live_after_cr': len(ok),
                'call_rate_min': round(float(cr.min()), 4), 'call_rate_median': round(float(cr.median()), 4),
                'call_rate_bins': {'<0.90': int((cr < 0.90).sum()), '0.90-0.95': int(((cr >= 0.90) & (cr < 0.95)).sum()),
                                   '0.95-0.98': int(((cr >= 0.95) & (cr < 0.98)).sum()),
                                   '0.98-0.99': int(((cr >= 0.98) & (cr < 0.99)).sum()), '>=0.99': int((cr >= 0.99).sum())},
                'cr_dropped': {t_: int((cr < t_).sum()) for t_ in (0.95, 0.98, 0.99)},
                'bins': {str(b_): int(v) for b_, v in pd.cut(ok.maf, edges, right=False, labels=labels)
                         .value_counts().reindex(labels).items()},
                'retained': {f: int(v.sum()) for f, v in floors.items()},
                'floor_mac': {f: int(ok.mac[v].min()) for f, v in floors.items()}}
        if (dist['allele']['retained']['MAF>=0.5%'], dist['residue']['retained']['MAF>=0.5%']) != \
                (sc['alleles_tested'], sc['residues_tested']):
            raise SystemExit(f'{co}: the record and the matrices of step 06 disagree')
        out['cohorts'][co] = {
            'n': sc['samples'], 'cases': sc['cases'], 'controls': sc['controls'], 'not_typed': 0,
            'alleles': {g: int(((tested['class'] == 'allele') & (tested.gene == g)).sum()) for g in genes},
            'residues': {g: int(((tested['class'] == 'residue') & (tested.gene == g)).sum()) for g in genes},
            'alleles_total': sc['alleles_tested'], 'residues_total': sc['residues_tested'], 'distribution': dist}
    return out


def drb_matrix(R, rec, ev):
    """Per paralogue: in how many samples with a result the exported copies (sum of the allele dosages)
    equal the copies DRB1 implies, and how many samples are blank."""
    ad = pd.read_csv(R / '05.exports' / 'allele_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    out = {}
    for g in C.PARALOGUES:
        tot = ad[[c for c in ad.columns if c.split('*')[0] == g]].sum(axis=1, min_count=1)
        r = rec[rec.gene == g].set_index('sample_id')
        e = ev[ev.gene == g].set_index('sample_id')
        exp = pd.to_numeric(r.expected, errors='coerce')
        ids = e.index[e.state.isin(['present', 'absent']) & exp.reindex(e.index).notna()]
        out[g] = {'with_result': len(ids), 'match': int((tot[ids] == exp[ids]).sum()), 'blank': int(tot.isna().sum()),
                  'with_genotype': int(tot.notna().sum()), 'with_drb1': int(exp.notna().sum())}
        # people with a genotype but no DRB1 prediction keep HLA-HD's own call
        out[g]['no_drb1'] = out[g]['with_genotype'] - out[g]['with_result']
        assert out[g]['match'] == out[g]['with_result'], g
    return out


def read_support(R):
    """Reads behind a heterozygous call, per gene (medians): all reads assigned to the pair, and the share that
    fits both alleles and so cannot tell them apart. The classical genes are the eight that present peptides."""
    e = pd.read_csv(R / '02.evidence' / 'call_evidence.tsv', sep='\t', low_memory=False)
    h = e[e.state == 'heterozygous'].assign(tot=lambda d: d.unique_1 + d.unique_2 + d.shared)
    h = h.assign(sh=h.shared / h.tot)
    by = h.groupby('gene').agg(reads=('tot', 'median'), shared=('sh', 'median'))
    classical = ['A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
    other = [g for g in by.index if g not in classical]
    return {'reads_min': float(by.reads.min()), 'reads_max': float(by.reads.max()),
            'shared_classical': [float(by.shared[classical].min()), float(by.shared[classical].max())],
            'shared_other': [float(by.shared[other].min()), float(by.shared[other].max())],
            'other_genes': other}


def drb_examples(rec, calls):
    """One real case per outcome of the DRB1 check (first in sample order; no sample ids kept): DRB1 alleles,
    HLA-HD's call (2-field), reads assigned and reads unique to each allele, and what was recorded."""
    cl = calls[calls.gene.isin(C.PARALOGUES)].set_index(['sample_id', 'gene'])
    flags = rec.pivot(index='sample_id', columns='gene', values='flag')
    d1 = rec.drop_duplicates('sample_id').set_index('sample_id')[['drb1_1', 'drb1_2']]
    num = lambda x: int(float(x)) if x not in ('', None) else 0          # noqa: E731

    def case(sid, gene):
        c = cl.loc[(sid, gene)]
        return {'drb1': d1.loc[sid].tolist(), 'gene': gene, 'raw': [a for a in (c.raw_allele_1, c.raw_allele_2) if a],
                'reads': [num(c.reads_1), num(c.reads_2)], 'unique': [num(c.unique_1), num(c.unique_2)],
                'final': c.allele_1, 'copies': num(c.copies), 'state': c.state}
    out = {}
    pair = d1.drb1_1 + d1.drb1_2
    for sid in d1.index[pair.str.contains('15:02') & pair.str.contains('04:05')]:
        if (flags.loc[sid] == '').all() and cl.loc[(sid, 'DRB4')].allele_1 and cl.loc[(sid, 'DRB5')].allele_1:
            out['agree'] = [case(sid, 'DRB4'), case(sid, 'DRB5')]
            break
    for sid in rec[(rec.gene == 'DRB5') & (rec.flag == 'drb_reconciled')].sample_id:
        c = cl.loc[(sid, 'DRB5')]
        if c.raw_allele_1 and c.raw_allele_2 and c.allele_1:
            out['extra'] = case(sid, 'DRB5')
            break
    for sid in rec[(rec.gene == 'DRB3') & (rec.flag == 'drb_unexpected')].sample_id:
        if cl.loc[(sid, 'DRB3')].raw_allele_1:
            out['unexpected'] = case(sid, 'DRB3')
            break
    sid = rec[(rec.gene == 'DRB4') & (rec.flag == 'drb_missing')].sample_id.iloc[0]
    out['missing'] = case(sid, 'DRB4')
    return out


def encoding_example(R, gene='A', a1='A*24:02', a2='A*11:01', positions=('9', '62')):
    """The worked example of slide 19: how many columns one gene contributes, and the residues the two example
    alleles carry at the shown positions (IPD-IMGT/HLA 3.64.0 protein alignment)."""
    head = pd.read_csv(R / '05.exports' / 'allele_dosage.tsv', sep='\t', nrows=0).columns
    rhead = pd.read_csv(R / '05.exports' / 'residue_dosage.tsv', sep='\t', nrows=0).columns
    sites = pd.read_csv(R / '05.exports' / 'residue_sites.tsv', sep='\t').set_index('gene')
    ref = pd.read_csv(R / '05.exports' / 'residue_reference' / f'{gene}.tsv', sep='\t', index_col=0, dtype=str)
    full = {a: next(x for x in ref.index if x.startswith(a + ':')) for a in (a1, a2)}
    return {'gene': gene, 'alleles': [a1, a2],
            'allele_columns': int(sum(1 for c in head if c.startswith(gene + '*'))),
            'residue_columns': int(sum(1 for c in rhead if c.startswith(gene + ':'))),
            'polymorphic_positions': int(sites.at[gene, 'polymorphic']),
            'residues': {p: {a: ref.at[full[a], p] for a in (a1, a2)} for p in positions},
            'residues_seen': {p: sorted(set(ref[p]) - {'*', '.'}) for p in positions}}


MATRIX_GENES = ('A', 'B', 'DRB1', 'DQB1')


def matrix_excerpt(R, genes=MATRIX_GENES, n_people=10, n_alleles=5, n_res_cols=6, anchor=(('A', '9'),), seed=7):
    """A real slice of the matrices the association reads (06.assoc_inputs/full_mainland): a few people × several
    genes, with allele columns and residue columns. Within a gene the most frequent alleles are taken; residues are
    taken as whole positions, every column of the position, so that a person's columns at a position add up to 2
    as on the residue slide. Only complete positions are eligible (every amino acid seen there kept a column; at
    the others a carrier of a residue below the MAF floor adds up to less than 2); per gene the two most variable
    that fit in n_res_cols columns, HLA-A always with position 9 of the residue slide. People are drawn at random
    from the full cohort; no sample ids are kept."""
    F = ASSOC / 'full_mainland'
    ad = pd.read_csv(F / 'allele_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    rd = pd.read_csv(F / 'residue_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    pos = {}                                                          # position -> its residue columns
    for c in rd.columns:
        pos.setdefault(c.rsplit(':', 1)[0], []).append(c)
    complete, short_people = {}, []
    for p, cols in pos.items():
        s = rd[cols].sum(axis=1, min_count=len(cols)).dropna()
        complete[p] = bool((s == 2).all())
        if not complete[p]:
            short_people.append(int((s < 2).sum()))
    anchor = dict(anchor)
    out, used = [], []
    for g in genes:
        a = ad[[c for c in ad.columns if c.split('*')[0] == g]]
        acols = a.sum().sort_values(ascending=False).head(n_alleles).index.tolist()
        r = rd[[c for c in rd.columns if c.split(':')[0] == g]]
        het = {p: 1 - float(((rd[cols].mean() / 2) ** 2).sum()) for p, cols in pos.items()
               if p.split(':')[0] == g and complete[p]}
        pairs = [(p1, p2) for p1 in het for p2 in het if p1 < p2 and len(pos[p1]) + len(pos[p2]) <= n_res_cols
                 and (g not in anchor or f'{g}:{anchor[g]}' in (p1, p2))]
        best = sorted(max(pairs, key=lambda t: het[t[0]] + het[t[1]]), key=lambda p: int(p.split(':')[1]))
        rcols = [c for p in best for c in pos[p]]
        out.append({'gene': g, 'alleles': [c.split('*')[1] for c in acols],
                    'residues': [c.split(':', 1)[1].replace(':', ' ') for c in rcols],
                    'residue_groups': [len(pos[p]) for p in best],
                    'n_alleles': int(a.shape[1]), 'n_residues': int(r.shape[1])})
        used.append((acols, rcols))
    keep = ad.index[~ad[[c for a_, _ in used for c in a_]].isna().any(axis=1)
                    & ~rd[[c for _, r_ in used for c in r_]].isna().any(axis=1)]
    ids = list(pd.Series(keep).sample(n_people, random_state=seed))
    # 'other': a person's remaining copies at the gene — alleles with a column not shown here, and alleles below
    # the MAF floor, which have no column. Not a column of the matrix; read from the person's genotype (export).
    ex = pd.read_csv(R / '05.exports' / 'allele_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    ex = ex.loc[ad.index]
    for g in out:
        gcols = [c for c in ex.columns if c.split('*')[0] == g['gene']]
        none = [c for c in gcols if c not in ad.columns]
        g['no_column_copies'] = float(ex[none].sum().sum() / ex[gcols].sum().sum())
    rows = [{'alleles': [[int(ad.at[i, c]) for c in a_] for a_, _ in used],
             'other': [int(ex.loc[i, [c for c in ex.columns if c.split('*')[0] == g['gene']]].sum())
                       - int(ad.loc[i, a_].sum()) for g, (a_, _) in zip(out, used)],
             'residues': [[int(rd.at[i, c]) for c in r_] for _, r_ in used]} for i in ids]
    for row in rows:                                                  # the properties the slice is chosen to show
        for g, av, o, rv in zip(out, row['alleles'], row['other'], row['residues']):
            assert sum(av) + o == 2, (g['gene'], av, o)
            k = 0
            for n in g['residue_groups']:
                assert sum(rv[k:k + n]) == 2, (g['gene'], rv)
                k += n
    return {'genes': out, 'people': rows, 'n_genes_total': len(C.ASSOC_GENES), 'n_samples': len(ad),
            'positions': len(pos), 'positions_complete': sum(complete.values()),
            'short_people_median': float(np.median(short_people)), 'short_people_max': max(short_people)}


STAGES = [('typed', 'exported', 'Typed by HLA-HD', 'columns that vary in the cohort'),
          ('genes', 'genes', 'DRB3/4/5 excluded', 'their presence follows DRB1'),
          ('callrate', 'call_rate', 'Call rate ≥ 95%', 'in cases and in controls'),
          ('maf', 'maf', 'MAF ≥ 0.5%', 'the frequency floor')]


def allele_coverage(R):
    """Share of the allele copies typed in a cohort that belong to an allele with a column in the final matrix
    (06.assoc_inputs), per gene of the association and over all 16. The rest are copies of alleles below the MAF
    floor: they have no allele column and enter the association through the residue columns only."""
    ex = pd.read_csv(R / '05.exports' / 'allele_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    out = {}
    for co in COHORTS:
        ad = pd.read_csv(ASSOC / co / 'allele_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
        e = ex.loc[ad.index]
        assert (e[ad.columns].fillna(-1).astype(int) == ad.fillna(-1).astype(int)).all().all()
        typed, kept = {}, {}
        for gene in C.ASSOC_GENES:
            gcols = [c for c in e.columns if c.split('*')[0] == gene]
            typed[gene] = int(e[gcols].sum().sum())
            kept[gene] = int(e[[c for c in gcols if c in ad.columns]].sum().sum())
        out[co] = {'by_gene': {gene: kept[gene] / typed[gene] for gene in C.ASSOC_GENES},
                   'all': sum(kept.values()) / sum(typed.values())}
    return out


def residue_coverage(R):
    """At the positions that keep at least one residue column in a cohort's final matrix, the share of the copies
    (two per person whose residue there is determined) that fall in a kept column; the rest are amino acids below
    the MAF floor. A residue left undetermined blanks the whole position, so the copies are the non-blank cells."""
    ex = pd.read_csv(R / '05.exports' / 'residue_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    out = {}
    for co in COHORTS:
        rd = pd.read_csv(ASSOC / co / 'residue_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
        e = ex.loc[rd.index]
        assert (e[rd.columns].fillna(-1).astype(int) == rd.fillna(-1).astype(int)).all().all()
        kept_pos = {c.rsplit(':', 1)[0] for c in rd.columns}
        at_kept = e[[c for c in e.columns if c.rsplit(':', 1)[0] in kept_pos]]
        out[co] = {'positions': len(kept_pos), 'all': float(rd.sum().sum() / at_kept.sum().sum())}
    return out


def pgroup_splits(R, S):
    """The P groups in which a cohort's final matrix tests more than one allele column. The alleles of a P group
    share their protein over the peptide-binding exons and differ only outside them, so the ToMMo comparison,
    made on P groups, cannot check these splits. Counts per cohort, over the 16 genes and over the genes compared
    with ToMMo; the rows from the full mainland cohort, each column with its frequency (copies over copies typed)."""
    pm = pd.read_csv(R / '05.exports' / 'allele_pgroup_map.tsv', sep='\t').set_index('allele').p_group
    mk = pd.read_csv(ASSOC / 'markers.tsv', sep='\t')
    compared = set(S['coverage'])
    out = {'per_cohort': {}, 'rows': []}
    for co in COHORTS:
        m = mk[(mk.cohort == co) & (mk['class'] == 'allele') & (mk.status == 'tested')].copy()
        cols = pd.read_csv(ASSOC / co / 'allele_dosage.tsv', sep='\t', nrows=0).columns[1:]
        assert set(m.marker) == set(cols), co                       # the record and the matrix agree
        m['pg'] = m.marker.map(pm)
        assert m.pg.notna().all(), co
        size = m.groupby('pg').marker.transform('size')
        sel = m[size > 1]
        cmp_sel = sel[sel.gene.isin(compared)]
        out['per_cohort'][co] = {'tested': len(m), 'groups': int(sel.pg.nunique()), 'columns': len(sel),
                                 'groups_compared': int(cmp_sel.pg.nunique()), 'columns_compared': len(cmp_sel)}
        if co != 'full_mainland':
            continue
        sel = sel.assign(freq=sel.copies / (2 * sel.called), rank=sel.gene.map(list(C.ASSOC_GENES).index))
        for pg, d in sel.sort_values(['rank', 'pg']).groupby('pg', sort=False):
            rep = pg[:-1] if pg.endswith('P') else pg               # the group's representative name
            d = d.assign(first=d.marker != rep).sort_values(['first', 'freq'], ascending=[True, False])
            out['rows'].append({'gene': d.gene.iloc[0], 'p_group': pg, 'compared': d.gene.iloc[0] in compared,
                                'columns': [{'allele': a, 'freq': float(f)} for a, f in zip(d.marker, d.freq)]})
    return out


def harmonization(R, S):
    """How the controls' calls and the ToMMo table are made like for like before any comparison (step 03 and
    panels.py): the same genes, P groups read from one IPD-IMGT 3.64.0 file on both sides, null alleles dropped
    on both sides, a control counted at a gene only with two 2-field alleles, and the panel's dictionary age.
    Every count is over the genes whose frequencies are compared; examples are the largest cases in the data."""
    genes = list(S['coverage'])
    loc = pd.read_csv(R / '03.panels' / 'locus_concordance.tsv', sep='\t')
    loc = loc[(loc.panel == 'ToMMo') & loc.gene.isin(genes)]
    tm = pd.read_csv(C.TOMMO)
    tm = tm[tm.gene.isin(genes)]
    null = tm.allele.str.rstrip('P').str.endswith('N')
    named_p = tm.allele.str.endswith('P') & ~null
    ex = pd.read_csv(R / '05.exports' / 'allele_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    ph = pd.read_csv(MODEL_INPUTS / 'full_mainland' / 'phenotype' / 'pheno.tsv', sep='\t', dtype={'IID': str})
    seen = ex.loc[ph.loc[ph.PHENO1 == 0, 'IID']].sum()
    seen = seen[(seen > 0) & seen.index.str.split('*').str[0].isin(genes) & ~seen.index.str.endswith('N')]
    merged = {a: N.pg_of(a) for a in seen.index if N.pg_of(a) != a}
    ex_merge = max(merged, key=lambda a: seen[a])                   # our commonest allele named by another group
    ac = pd.read_csv(R / '03.panels' / 'allele_concordance.tsv', sep='\t')
    ac = ac[(ac.panel == 'ToMMo') & ac.gene.isin(genes) & (ac.count_ours > 0) & (ac.count_panel > 0)]
    ac = ac[~ac.allele.map(N.pg_label).str.endswith('P')]            # alone in its group, carried by both sides
    lone = ac.sort_values('count_panel', ascending=False).iloc[0]
    tnull = tm[null].sort_values('count', ascending=False).iloc[0]
    cov = S['coverage']
    return {
        'genes': genes,
        'controls_min': int(loc.n_samples.min()), 'controls_max': int(loc.n_samples.max()),
        'ours_two_field': sum(cov[g]['two_field'] for g in genes),
        'ours_groups': sum(sum(cov[g]['alleles'].values()) for g in genes),
        'tommo_entries': int(len(tm)), 'tommo_named_p': int(named_p.sum()),
        'tommo_lone': int((~named_p & ~null).sum()), 'tommo_null_entries': int(null.sum()),
        'null_copies_ours': int(loc.null_chr_ours.sum()), 'null_copies_tommo': int(loc.null_chr_panel.sum()),
        'example_merge': {'allele': ex_merge, 'group': N.pg_label(ex_merge), 'copies': int(seen[ex_merge])},
        'example_lone': {'allele': lone.allele, 'copies_ours': int(lone.count_ours), 'copies_tommo': int(lone.count_panel)},
        'example_null': {'allele': tnull.allele, 'copies': int(tnull['count'])},
        'too_new_groups': sum(cov[g]['alleles']['too_new'] for g in genes),
        # the rescaling worked once: ToMMo's commonest HLA-A entry over all copies, then over protein-coding ones
        'renorm': (lambda d: {'allele': d.sort_values('count').allele.iloc[-1], 'count': int(d['count'].max()),
                              'total': int(d['count'].sum()),
                              'nonnull': int(d.loc[~d.allele.str.rstrip('P').str.endswith('N'), 'count'].sum())})(
            pd.read_csv(C.TOMMO).query("gene == 'A'")),
    }


def pgroup_mapping(R, S):
    """How the controls' 2-field alleles and the ToMMo table meet at P-group resolution, per compared gene, and
    the frequencies behind every P group in which the association tests more than one allele column. Frequencies
    are shares of the controls' non-null copies at the gene, from the controls with a 2-field genotype there
    (the denominator of the ToMMo comparison, checked against locus_concordance.tsv); ToMMo's are its own, nulls
    dropped (allele_concordance.tsv)."""
    genes = list(S['coverage'])
    tm = pd.read_csv(C.TOMMO)
    ph = pd.read_csv(MODEL_INPUTS / 'full_mainland' / 'phenotype' / 'pheno.tsv', sep='\t', dtype={'IID': str})
    ex = pd.read_csv(R / '05.exports' / 'allele_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    X = ex.loc[ph.loc[ph.PHENO1 == 0, 'IID']]
    loc = pd.read_csv(R / '03.panels' / 'locus_concordance.tsv', sep='\t')
    loc = loc[loc.panel == 'ToMMo'].set_index('gene')
    ac = pd.read_csv(R / '03.panels' / 'allele_concordance.tsv', sep='\t')
    ac = ac[ac.panel == 'ToMMo'].set_index(['gene', 'allele'])

    def copies(g):
        """The controls' non-null copies per 2-field allele at gene g, and their total."""
        cols = [c for c in X.columns if c.split('*')[0] == g and not c.endswith('N')]
        s = X[cols].sum()
        return s[s > 0], float(s.sum())

    per_gene = []
    for g in genes:
        t_ = tm[tm.gene == g]
        nonnull = ~t_.allele.str.rstrip('P').str.endswith('N')
        s, tot = copies(g)
        assert int(tot) == int(loc.loc[g, 'n_chr_ours']), g             # the comparison's own denominator
        keys = pd.Series(s.index.map(N.pg_of), index=s.index)
        size = keys.map(keys.value_counts())
        cmp_ = ac.loc[g]                                                  # every allele of either table, one row each
        both = cmp_[(cmp_.count_ours > 0) & (cmp_.count_panel > 0)]
        per_gene.append({'gene': g,
                         'tommo_p': int((t_.allele.str.endswith('P') & nonnull).sum()),
                         'tommo_plain': int((~t_.allele.str.endswith('P') & nonnull).sum()),
                         'ours': int(len(s)), 'ours_shared': int((size > 1).sum()),
                         'shared_groups': int(keys[size > 1].nunique()), 'ours_keys': int(keys.nunique()),
                         'in_tommo': int(sum(ac.loc[(g, k), 'count_panel'] > 0 for k in keys.unique())),
                         'tommo_keys': int((cmp_.count_panel > 0).sum()), 'either': int(len(cmp_)),
                         # r over the alleles both tables carry: a sensitivity check, undefined with 2 points
                         'r_shared': float(np.corrcoef(both.freq_ours, both.freq_panel)[0, 1]) if len(both) > 2 else None})
        pg_ = per_gene[-1]
        assert pg_['ours_keys'] == sum(S['coverage'][g]['alleles'].values()), g
        assert pg_['either'] == int(loc.loc[g, 'n_alleles']) and pg_['in_tommo'] == int(loc.loc[g, 'n_shared']), g
        assert pg_['either'] == pg_['ours_keys'] + pg_['tommo_keys'] - pg_['in_tommo'], g

    rows = []
    for r in S['pgroup_splits']['rows']:
        g, key = r['gene'], N.pg_of(r['columns'][0]['allele'])
        s, tot = copies(g)
        members = s[[N.pg_of(a) == key for a in s.index]]
        tested = [c['allele'] for c in r['columns']]
        others = members.drop([a for a in tested if a in members.index])
        row = {'gene': g, 'p_group': r['p_group'], 'compared': r['compared'],
               'tested': [{'allele': a, 'freq': float(s.get(a, 0) / tot)} for a in tested],
               'others_n': int(len(others)), 'others_freq': float(others.sum() / tot),
               'ours_freq': float(members.sum() / tot)}
        if r['compared']:
            row['tommo_freq'] = float(ac.loc[(g, key), 'freq_panel'])
            assert int(ac.loc[(g, key), 'count_ours']) == int(members.sum()), (g, key)   # same copies as the comparison
        rows.append(row)
    return {'per_gene': per_gene, 'rows': rows}


def worked_example(R, S, gene='DQA1'):
    """One gene followed through the ToMMo comparison: every allele at P-group resolution found in either table,
    whether it is in both, only in our controls, or only in ToMMo, its frequency on each side, and the 2-field
    alleles of the controls merged into it. DQA1 is small enough to show whole and has P groups that merge."""
    ac = pd.read_csv(R / '03.panels' / 'allele_concordance.tsv', sep='\t')
    a = ac[(ac.panel == 'ToMMo') & (ac.gene == gene)].copy()
    ph = pd.read_csv(MODEL_INPUTS / 'full_mainland' / 'phenotype' / 'pheno.tsv', sep='\t', dtype={'IID': str})
    ex = pd.read_csv(R / '05.exports' / 'allele_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    X = ex.loc[ph.loc[ph.PHENO1 == 0, 'IID']]
    s = X[[c for c in X.columns if c.split('*')[0] == gene and not c.endswith('N')]].sum()
    s = s[s > 0]
    members = {}
    for al_, n in s.items():
        members.setdefault(N.pg_of(al_), []).append((al_, int(n)))
    rows = []
    for _, r in a.iterrows():
        where = 'both' if r.count_ours > 0 and r.count_panel > 0 else ('ours' if r.count_ours > 0 else 'tommo')
        m = sorted(members.get(r.allele, []), key=lambda x: -x[1])
        assert sum(n for _, n in m) == int(r.count_ours), r.allele              # the comparison's own counts
        rows.append({'allele': N.pg_label(r.allele), 'where': where, 'freq_ours': float(r.freq_ours),
                     'freq_tommo': float(r.freq_panel), 'copies_ours': int(r.count_ours),
                     'copies_tommo': int(r.count_panel), 'members': [{'allele': x, 'copies': n} for x, n in m]})
    rows.sort(key=lambda d: ({'both': 0, 'ours': 1, 'tommo': 2}[d['where']], -d['freq_tommo'], -d['freq_ours']))
    loc = pd.read_csv(R / '03.panels' / 'locus_concordance.tsv', sep='\t')
    loc = loc[(loc.panel == 'ToMMo') & (loc.gene == gene)].iloc[0]
    pm = next(x for x in S['pgroup_mapping']['per_gene'] if x['gene'] == gene)
    assert len(rows) == pm['either'] and sum(r['where'] == 'both' for r in rows) == pm['in_tommo']
    return {'gene': gene, 'rows': rows, 'r': float(loc.r), 'tommo_copies': int(loc.n_chr_panel), 'two_field': pm['ours'], 'ours_keys': pm['ours_keys'],
            'tommo_keys': pm['tommo_keys'], 'both': pm['in_tommo'], 'either': pm['either'],
            'ours_shared': pm['ours_shared'], 'shared_groups': pm['shared_groups']}


def residue_loss(R):
    """How much call rate the residue positions of the 16 association genes lose to the alignment, full mainland
    cohort: per exported position, the gene's 2-field completeness minus the share of people whose residue there
    is determined. Counts of positions losing more than 1% and more than 5%."""
    ph = pd.read_csv(MODEL_INPUTS / 'full_mainland' / 'phenotype' / 'pheno.tsv', sep='\t', dtype={'IID': str})
    ad = pd.read_csv(R / '05.exports' / 'allele_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    rd = pd.read_csv(R / '05.exports' / 'residue_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    ad, rd = ad.loc[ph.IID], rd.loc[ph.IID]
    loss = []
    for g in C.ASSOC_GENES:
        comp = ad[[c for c in ad.columns if c.split('*')[0] == g]].notna().all(axis=1).mean()
        cols = [c for c in rd.columns if c.split(':')[0] == g]
        for pos_, cc in pd.Series(cols).groupby(pd.Series(cols).map(lambda c: c.rsplit(':', 1)[0])):
            loss.append(comp - rd[list(cc)].notna().all(axis=1).mean())
    loss = pd.Series(loss)
    return {'positions': int(len(loss)), 'over_1pct': int((loss > 0.01).sum()), 'over_5pct': int((loss > 0.05).sum())}


def qc_per_gene(R):
    """The record of the marker QC, per cohort, step and gene, as step 06 wrote it (qc_record.tsv). A column
    counts at a step while it varies in that cohort (at least one copy of its less common form). The first step
    covers the 19 reported genes, the later ones the 16 of the association; a gene with nothing left counts 0."""
    rec = pd.read_csv(ASSOC / 'qc_record.tsv', sep='\t')
    out = {}
    for co in COHORTS:
        stages = {}
        for key, step, _, _ in STAGES:
            r = rec[(rec.cohort == co) & (rec.step == step)]
            genes = C.REPORTED_GENES if key == 'typed' else C.ASSOC_GENES
            stages[key] = {plural: {g: int(r.loc[(r['class'] == kind) & (r.gene == g), 'columns'].sum()) for g in genes}
                           for kind, plural in (('allele', 'alleles'), ('residue', 'residues'))}
        out[co] = stages
    return {'stages': [{'key': k, 'name': n, 'rule': r} for k, _, n, r in STAGES], 'per_cohort': out,
            'genes_all': list(C.REPORTED_GENES), 'genes_kept': list(C.ASSOC_GENES)}


def floor_rationale(R, cohorts):
    """The evidence the deck shows for the two floors: per-gene completeness (the ceiling every marker of a gene
    inherits as its call rate), and the patient chromosomes expected at each candidate MAF floor."""
    ad = pd.read_csv(R / '05.exports' / 'allele_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    genes = [g for g in C.REPORTED_GENES if g not in C.PARALOGUES]
    comp, per_cohort = {}, {}
    for g in genes:
        cols = [c for c in ad.columns if c.split('*')[0] == g]
        comp[g] = round(float(ad[cols].sum(axis=1, min_count=1).notna().mean()), 4)
    rd = pd.read_csv(R / '05.exports' / 'residue_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    rd = rd.loc[:, [c.split(':')[0] in genes for c in rd.columns]]
    pos = rd.columns.str.rsplit(':', n=1).str[0]
    worst = {}
    for pid in pd.unique(pos):                       # lowest per-position call rate within each gene
        cols = rd.columns[pos == pid]
        sm = rd[cols].sum(axis=1, min_count=len(cols))
        cr = float(((sm == 2).sum()) / len(rd))
        g = pid.split(':')[0]
        if g not in worst or cr < worst[g]:
            worst[g] = round(cr, 4)
    for co in COHORTS:
        ids = pd.read_csv(MODEL_INPUTS / co / 'phenotype' / 'pheno.tsv', sep='\t', dtype={'IID': str}).IID
        d = ad.loc[[i for i in ids if i in ad.index]]
        per_cohort[co] = {g: round(float(d[[c for c in d.columns if c.split('*')[0] == g]].sum(axis=1, min_count=1)
                                        .notna().mean()), 4) for g in genes}
    exp = {}
    for co, v in cohorts.items():
        share = v['cases'] / v['n']
        exp[co] = {f'{f}': round(f * 2 * v['n'] * share, 1) for f in (0.005, 0.01)}
        exp[co]['mac20'] = round(20 * share, 1)
    ph = pd.read_csv(MODEL_INPUTS / 'full_mainland' / 'phenotype' / 'pheno.tsv', sep='\t', dtype={'IID': str}).set_index('IID')
    d = ad.loc[ph.index]
    by_group = {}
    for g in genes:
        has = d[[c for c in d.columns if c.split('*')[0] == g]].sum(axis=1, min_count=1).notna()
        by_group[g] = {'cases': round(float(has[ph.PHENO1 == 1].mean()), 4),
                       'controls': round(float(has[ph.PHENO1 == 0].mean()), 4)}
    low_g = min(genes, key=lambda g: min(by_group[g].values()))
    low_grp = min(by_group[low_g], key=by_group[low_g].get)
    return {'gene_completeness': comp, 'gene_completeness_by_group': by_group,
            'lowest_group_gene': low_g, 'lowest_group': low_grp, 'lowest_group_value': by_group[low_g][low_grp], 'gene_completeness_by_cohort': per_cohort, 'worst_residue_position': worst,
            'lowest_gene': min(comp, key=comp.get), 'lowest': min(comp.values()),
            'expected_patient_chromosomes': exp}


def main():
    R = C.RESULTS
    prov = json.loads((R / '00.inputs' / 'provenance.json').read_text())
    si = pd.read_csv(R / '00.inputs' / 'sample_inputs.tsv', sep='\t')
    gs = pd.read_csv(R / '02.evidence' / 'gene_summary.tsv', sep='\t')
    ev = pd.read_csv(R / '02.evidence' / 'call_evidence.tsv', sep='\t', dtype=str).fillna('')
    sq = pd.read_csv(R / '02.evidence' / 'sample_qc.tsv', sep='\t')
    drb = pd.read_csv(R / '02.evidence' / 'drb_consistency.tsv', sep='\t', dtype=str)
    rec = pd.read_csv(R / '01.calls' / 'drb_reconciliation.tsv', sep='\t', dtype=str).fillna('')
    calls = pd.read_csv(R / '01.calls' / 'calls_long.tsv', sep='\t', dtype=str).fillna('')
    loc = pd.read_csv(R / '03.panels' / 'locus_concordance.tsv', sep='\t')
    hz = pd.read_csv(R / '03.panels' / 'homozygosity.tsv', sep='\t')
    cv = pd.read_csv(R / '03.panels' / 'confidence_vs_panel.tsv', sep='\t')
    sh = pd.read_csv(R / '03.panels' / 'shift_direction.tsv', sep='\t')
    unc = pd.read_csv(R / '03.panels' / 'unconfirmed_alleles.tsv', sep='\t')
    ceil = pd.read_csv(R / '03.panels' / 'panel_ceiling.tsv', sep='\t')
    tiers = pd.read_csv(R / '04.confidence' / 'locus_tiers.tsv', sep='\t')
    ad = pd.read_csv(R / '05.exports' / 'allele_dosage.tsv', sep='\t', nrows=1)
    rd = pd.read_csv(R / '05.exports' / 'residue_dosage.tsv', sep='\t', nrows=1)

    S = {}
    S['cohort'] = {'n': len(si), 'groups': si.group.value_counts().to_dict(),
                   'read_pairs_median': int(si.read_pairs.median()),
                   'read_pairs_p1': int(np.percentile(si.read_pairs, C.OUTLIER_READ_PAIRS_PCTL)),
                   'read_pairs_min': int(si.read_pairs.min()), 'read_pairs_max': int(si.read_pairs.max()),
                   'validation_ok': int((si.validation == 'ok').sum())}
    S['run'] = {'hlahd_version': prov['hlahd_version'], 'regions': prov['regions'],
                'genes_typed': len(set(calls.gene)), 'genes_reported': len(C.REPORTED_GENES),
                'imgt': '3.64.0', 'gene_split': '3.50.0'}
    al = pd.concat([ev[['gene', 'basis_1', 'allele_1', 'prior_ratio_1', 'copies']].set_axis(
                        ['gene', 'basis', 'allele', 'q', 'copies'], axis=1),
                    ev[['gene', 'basis_2', 'allele_2', 'prior_ratio_2', 'copies']].set_axis(
                        ['gene', 'basis', 'allele', 'q', 'copies'], axis=1).assign(second=1)])
    al = al[~((al.second == 1) & (al.copies == '1'))]
    al = al[al.basis != '']
    q = pd.to_numeric(al.q, errors='coerce')
    S['basis'] = {'total': len(al), 'agreed': int((al.basis == 'agreed').sum()),
                  'prior': int((al.basis == 'prior').sum()),
                  'prior_decisive': int(((al.basis == 'prior') & (q >= C.PRIOR_DOMINANCE)).sum()),
                  'pgroup': int((al.basis == 'pgroup').sum()), 'none': int((al.basis == 'none').sum()),
                  'prior_by_gene': al[al.basis == 'prior'].gene.value_counts().to_dict()}
    S['genes'] = {r.gene: {k: (None if pd.isna(v) else v) for k, v in r._asdict().items() if k != 'Index'}
                  for r in gs.itertuples()}
    S['states'] = ev.state.value_counts().to_dict()
    S['confidence'] = ev.confidence.value_counts().to_dict()
    S['calls_reported'] = len(ev)
    h = ev[(ev.balance_assessed == 'True') & ev.gene.isin(C.OBLIGATE_GENES)]
    S['imbalance'] = {'assessed': len(h), 'imbalanced': int((h.imbalanced == 'True').sum())}
    p = calls[calls.gene.isin(C.PARALOGUES) & (calls.raw == 'two')].copy()
    u1, u2 = pd.to_numeric(p.unique_1), pd.to_numeric(p.unique_2)
    p = p[(u1 + u2) >= C.MIN_UNIQUE_READS]
    ms = np.minimum(u1, u2) / (u1 + u2)
    S['paralogue_two'] = {e: {'n': int((p.drb_expected == e).sum()),
                              'below': int(((p.drb_expected == e) & (ms[p.index] < C.MIN_MINOR_SHARE)).sum())}
                          for e in ('1', '2')}
    drb['samples'] = drb.samples.astype(int)
    S['drb'] = {}
    for g in C.PARALOGUES:
        d = drb[(drb.gene == g) & (drb.expected != 'unknown')]
        agree = int(d[d.expected == d.observed].samples.sum())
        S['drb'][g] = {'agree': agree, 'total': int(d.samples.sum()),
                       'flags': rec[rec.gene == g].flag.value_counts().to_dict(),
                       'two_where_one': int(d[(d.expected == '1') & (d.observed == '2')].samples.sum()),
                       'present_where_none': int(d[(d.expected == '0') & (d.observed != '0')].samples.sum())}
    S['samples'] = {'outliers': int(sq.outlier.sum()),
                    'reasons': sq.outlier_reason[sq.outlier].str.split(';').explode().value_counts().to_dict(),
                    'low_conf_hist': sq.low_conf.value_counts().sort_index().to_dict()}
    tom = loc[loc.panel == 'ToMMo'].set_index('gene')
    kg = loc[loc.panel == '1KG JPT'].set_index('gene')
    S['tommo'] = {g: {k: tom.at[g, k] for k in ('n_samples', 'n_chr_ours', 'r', 'D', 'D_null_mean', 'excess_D',
                                                 'D_shared', 'D_shared_high', 'unconfirmed',
                                                 'unconfirmed_assessable', 'unconfirmed_after_vintage',
                                                 'unconfirmed_assessable_high', 'frequency_comparable')}
                  for g in tom.index}
    S['tommo_vintage'] = int(tom.panel_vintage.iloc[0])
    S['freq_outlier'] = freq_outlier(R, 'DRB1')
    S['coverage'] = panel_coverage(R, S)
    S['kg_vintage'] = int(kg.panel_vintage.iloc[0])
    S['kg'] = {g: {k: kg.at[g, k] for k in ('n_chr_panel', 'r', 'D', 'D_null_mean', 'excess_D')} for g in kg.index}
    S['ceiling'] = ceil.set_index('gene').to_dict(orient='index')
    S['hom'] = hz.set_index('gene')[['obs_hom', 'exp_hom', 'ratio', 'dropout', 'dropout_lo', 'dropout_hi']] \
        .to_dict(orient='index')
    pooled = cv[cv.gene == 'all ToMMo loci'].set_index('confidence')
    S['conf_vs_panel'] = pooled[['n_chr', 'n_unconfirmed', 'unconfirmed_share']].to_dict(orient='index')
    S['shift'] = {'n': len(sh), 'ours': int((sh.kg_side == 'ours').sum()),
                  'tommo': int((sh.kg_side == 'ToMMo').sum())}
    S['unconfirmed_alleles'] = {'n_alleles': len(unc), 'n_chr': int(unc.count_ours.sum()),
                                'assessable_alleles': int(unc.assessable.sum()),
                                'assessable_chr': int(unc[unc.assessable == 1].count_ours.sum())}
    S['tiers'] = tiers.set_index('gene').to_dict(orient='index')
    S['tier_counts'] = tiers.tier.value_counts().sort_index().to_dict()
    S['exports'] = {'alleles': ad.shape[1] - 1, 'residues': rd.shape[1] - 1}
    S['thresholds'] = {k: getattr(C, k) for k in (
        'MIN_ALLELE_READS', 'MIN_MINOR_SHARE', 'MIN_UNIQUE_READS', 'PRIOR_DOMINANCE',
        'OUTLIER_READ_PAIRS_PCTL', 'OUTLIER_MIN_UNRESOLVED', 'TIER_MAX_UNCONFIRMED', 'TIER_MIN_HIGH_CONF',
        'TIER_MAX_DROPOUT', 'TIER_MAX_EXCESS_D')}
    import csv
    from collections import Counter
    cnt = Counter()
    with open(C.IMGT_HISTORY) as fh:
        rd = csv.reader(ln for ln in fh if not ln.startswith('#'))
        next(rd)
        for r in rd:
            nm = r[1].strip()
            if '*' in nm and nm != 'NA':
                cnt[nm.split('*')[0]] += 1
    S['imgt_alleles'] = {'total': sum(cnt.values()), **{g: cnt[g] for g in C.REPORTED_GENES}}
    # completeness: a result (any state but no call) and a 2-field genotype (absent DRB3/4/5 counts)
    two = ev.state.isin(['heterozygous', 'homozygous', 'present', 'absent'])
    res = ~ev.state.isin(C.NO_CALL_STATES)
    S['completeness'] = {'calls': len(ev), 'result': int(res.sum()), 'two_field': int(two.sum()),
                         'by_gene': {g: float(two[ev.gene == g].mean()) for g in C.REPORTED_GENES},
                         'parts_by_gene': {g: {'two_field': int((two & (ev.gene == g) & (ev.state != 'absent')).sum()),
                                               'absent': int(((ev.gene == g) & (ev.state == 'absent')).sum()),
                                               'one_field': int(((ev.gene == g) & (ev.state == 'ambiguous')).sum()),
                                               'no_call': int(((ev.gene == g) & ~res).sum())}
                                           for g in C.REPORTED_GENES}}
    S['tommo_n'] = 61424
    S['kg_n'] = 105
    S['assoc_basis'] = assoc_basis(R)
    S['drb_matrix'] = drb_matrix(R, rec, ev)
    S['drb_examples'] = drb_examples(rec, calls)
    S['floor_rationale'] = floor_rationale(R, S['assoc_basis']['cohorts'])
    S['encoding_example'] = encoding_example(R)
    S['matrix_excerpt'] = matrix_excerpt(R)
    S['qc_per_gene'] = qc_per_gene(R)
    S['allele_coverage'] = allele_coverage(R)
    S['pgroup_splits'] = pgroup_splits(R, S)
    S['pgroup_mapping'] = pgroup_mapping(R, S)
    S['worked_example'] = worked_example(R, S)
    S['residue_loss'] = residue_loss(R)
    S['harmonization'] = harmonization(R, S)
    S['read_support'] = read_support(R)
    S['residue_coverage'] = residue_coverage(R)
    (HERE / 'stats.json').write_text(json.dumps(clean(S), indent=1) + '\n')
    print('stats.json written;', {k: S[k] for k in ('basis', 'imbalance', 'paralogue_two', 'shift',
                                                     'tier_counts', 'exports', 'unconfirmed_alleles')})


if __name__ == '__main__':
    main()
