#!/usr/bin/env bash
# hla.typing.v2 — check every output against the rules in docs/METHODS.md.
# Prints one line per check; exits non-zero if any fails.
set -euo pipefail
cd "$(dirname "$0")"
exec /home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3 - <<'EOF'
import hashlib, json, re, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, 'scripts')
import config as C
import nomenclature as N

fails = []
def check(name, ok, detail=''):
    print(f"[{'ok' if ok else 'FAIL'}] {name}{(' — ' + detail) if detail else ''}")
    if not ok:
        fails.append(name)

R = C.RESULTS
ids = [l.split()[1] for l in C.COHORT_KEEP.read_text().splitlines() if l.strip()]
n = len(ids)

# 00 inputs
man = pd.read_csv(R / '00.inputs/archive_manifest.tsv', sep='\t', dtype=str)
check('archive: one per cohort sample', sorted(man.sample_id) == sorted(ids), f'{len(man)} vs {n}')
check('archive: 75 members each', set(man.n_members) == {'75'})
bad = [s for s, m in zip(man.sample_id, man.md5)
       if hashlib.md5((R / f'00.inputs/typing/{s}.tar.gz').read_bytes()).hexdigest() != m]
check('archive: md5 matches the manifest', not bad, f'{len(bad)} differ')
si = pd.read_csv(R / '00.inputs/sample_inputs.tsv', sep='\t', dtype=str)
check('provenance: validation ok for all', (si.validation == 'ok').all())
check('provenance: one HLA-HD version', si.hlahd_version.nunique() == 1, si.hlahd_version.iloc[0])
prov = json.loads((R / '00.inputs/provenance.json').read_text())
for k, v in prov['pinned'].items():
    if 'md5' in v and Path(v['path']).is_file():
        ok = hashlib.md5(Path(v['path']).read_bytes()).hexdigest() == v['md5']
        check(f'pinned file unchanged: {k}', ok)

# 01 calls
cl = pd.read_csv(R / '01.calls/calls_long.tsv', sep='\t', dtype=str).fillna('')
check('calls: samples x 33 genes', len(cl) == n * 33, f'{len(cl)} rows')
check('calls: 19 reported genes', sorted(cl[cl.reported == '1'].gene.unique()) == sorted(C.REPORTED_GENES))
rep = cl[cl.reported == '1']
two_field = re.compile(r'^[A-Z0-9]+\*\d+:\d+[NLSQ]?$')
al = pd.concat([rep.allele_1, rep.allele_2])
al = al[al != '']
check('calls: every allele is a 2-field name', al.map(lambda a: bool(two_field.match(a))).all(),
      f'{(~al.map(lambda a: bool(two_field.match(a)))).sum()} not')
nocall = rep.state.isin(C.NO_CALL_STATES + ['absent'])
check('calls: no allele on a no-call or absent state', ((rep[nocall].allele_1 == '') & (rep[nocall].allele_2 == '')).all())
called = rep.state.isin(['heterozygous', 'homozygous'])
check('calls: heterozygous/homozygous have both alleles', ((rep[called].allele_1 != '') & (rep[called].allele_2 != '')).all())
check('calls: heterozygous alleles differ, homozygous equal',
      (rep[rep.state == 'heterozygous'].allele_1 != rep[rep.state == 'heterozygous'].allele_2).all()
      and (rep[rep.state == 'homozygous'].allele_1 == rep[rep.state == 'homozygous'].allele_2).all())
par = rep[rep.gene.isin(C.PARALOGUES) & (rep.drb_expected != '') & (rep.state == 'present')]
check('calls: DRB3/4/5 copies follow DRB1', (par.copies.astype(int) == par.drb_expected.astype(int)).all())
check('calls: DRB3/4/5 absent where DRB1 implies none',
      (rep[rep.gene.isin(C.PARALOGUES) & (rep.drb_expected == '0')].state == 'absent').all())
check('calls: no 1-field allele written for an ambiguous call',
      (rep[rep.state == 'ambiguous'][['allele_1', 'allele_2']].apply(lambda r: all(not x or two_field.match(x) for x in r), axis=1)).all())

# 02 evidence
ev = pd.read_csv(R / '02.evidence/call_evidence.tsv', sep='\t', dtype=str).fillna('')
check('evidence: one row per sample x reported gene', len(ev) == n * 19)
check('evidence: confidence classes', set(ev.confidence) <= {'high', 'medium', 'low', 'no_call'})
nc = ev.state.isin(C.NO_CALL_STATES)
check('evidence: no call <=> no_call', ((ev.confidence == 'no_call') == nc).all())
hi = ev[ev.confidence == 'high']
check('evidence: high-confidence calls carry no weakness', (hi.weaknesses == '').all())
imb = ev[ev.imbalanced == 'True']
check('evidence: imbalanced <=> minor share < threshold',
      (pd.to_numeric(imb.minor_share) < C.MIN_MINOR_SHARE).all()
      and (pd.to_numeric(imb.unique_total) >= C.MIN_UNIQUE_READS).all())
sq = pd.read_csv(R / '02.evidence/sample_qc.tsv', sep='\t')
check('sample QC: one row per sample', len(sq) == n)

# 03 panels
loc = pd.read_csv(R / '03.panels/locus_concordance.tsv', sep='\t')
ctrl = int((si.group == C.CONTROL_GROUP).sum())
check('panels: controls only', (loc.n_samples <= ctrl).all(), f'max {loc.n_samples.max()} of {ctrl}')
check('panels: ToMMo and 1KG loci', sorted(loc[loc.panel == 'ToMMo'].gene) == sorted(C.TOMMO_LOCI)
      and sorted(loc[loc.panel == '1KG JPT'].gene) == sorted(C.KG_LOCI))
check('panels: assessable + after-vintage = all absent',
      np.allclose(loc.unconfirmed, loc.unconfirmed_assessable + loc.unconfirmed_after_vintage))
unc = pd.read_csv(R / '03.panels/unconfirmed_alleles.tsv', sep='\t')
tl = loc[loc.panel == 'ToMMo'].set_index('gene')
agg = unc[unc.assessable == 1].groupby('gene').count_ours.sum()
ok = all(abs(agg.get(g, 0) / tl.at[g, 'n_chr_ours'] - tl.at[g, 'unconfirmed_assessable']) < 1e-5 for g in tl.index)
check('panels: absent-allele list adds up to the locus shares', ok)

# 04 tiers
t = pd.read_csv(R / '04.confidence/locus_tiers.tsv', sep='\t').set_index('gene')
check('tiers: one per reported gene', sorted(t.index) == sorted(C.REPORTED_GENES))
re_t = {}
for g, r in t.iterrows():
    f = 0
    if not pd.isna(r.get('unconfirmed')):
        f += r.unconfirmed > C.TIER_MAX_UNCONFIRMED
    if not pd.isna(r.get('dropout_lo')):
        f += r.dropout_lo > C.TIER_MAX_DROPOUT
    f += r.share_high < C.TIER_MIN_HIGH_CONF
    re_t[g] = min((1 if r.external_panel != 'none' else 2) + f, 3)
check('tiers: recomputed from the thresholds', all(re_t[g] == t.at[g, 'tier'] for g in t.index))

# 05 exports
ad = pd.read_csv(R / '05.exports/allele_dosage.tsv', sep='\t', index_col=0, dtype={'sample_id': str})
rd = pd.read_csv(R / '05.exports/residue_dosage.tsv', sep='\t', index_col=0, dtype={'sample_id': str})
check('exports: sample order = sample_inputs', list(ad.index) == list(si.sample_id) == list(rd.index))
vals = pd.unique(ad.values.ravel())
check('exports: allele dosage in {0,1,2,blank}', set(v for v in vals if not pd.isna(v)) <= {0, 1, 2})
vals = pd.unique(rd.values.ravel())
check('exports: residue dosage in {0,1,2,blank}', set(v for v in vals if not pd.isna(v)) <= {0, 1, 2})
ok = True
for g in C.OBLIGATE_GENES:
    cols = [c for c in ad.columns if c.startswith(g + '*')]
    s = ad[cols].sum(axis=1, min_count=1)
    ok &= set(s.dropna().astype(int)) <= {2}
check('exports: obligate genes sum to 2 copies where called', ok)
for g in C.PARALOGUES:
    cols = [c for c in ad.columns if c.startswith(g + '*')]
    s = ad[cols].sum(axis=1, min_count=1).dropna().astype(int)
    check(f'exports: {g} dosage sums are 0, 1 or 2', set(s) <= {0, 1, 2}, str(sorted(set(s))))
um = pd.read_csv(R / '05.exports/allele_unmatched.tsv', sep='\t')
mm = pd.read_csv(R / '05.exports/allele_match_map.tsv', sep='\t')
check('exports: alleles matched to their own or a deeper alignment row, never a shorter one',
      set(mm.how.str.replace(r'\(.*\)', '', regex=True)) <= {'exact', 'prefix', 'no_match'})
ok = True
for a, g in zip(um.allele, um.gene):
    carriers = ad.index[ad[a].fillna(0) > 0]
    gcols = [c for c in rd.columns if c.split(':')[0] == g]
    ok &= bool(rd.loc[carriers, gcols].isna().all().all())
check('exports: carriers of an allele absent from the alignment are blank at every residue of the gene', ok,
      f'{len(um)} such allele(s): {", ".join(um.allele)}')
rob = rd.loc[:, [c.split(':')[0] in C.OBLIGATE_GENES for c in rd.columns]]
rpos = rob.columns.str.rsplit(':', n=1).str[0]
psum = rob.T.groupby(rpos.values).sum(min_count=1).T
pany = rob.notna().T.groupby(rpos.values).any().T
pall = rob.notna().T.groupby(rpos.values).all().T
check('exports: obligate residues sum to 2 at a called position, else the whole position is blank',
      bool(((psum == 2) | psum.isna()).all().all()) and bool((pany == pall).all().all()))

# 06 association inputs (written here; no test is run)
mk = pd.read_csv(R / '06.assoc_inputs/markers.tsv', sep='\t', dtype={'position': str, 'residue': str})
summ = json.loads((R / '06.assoc_inputs/summary.json').read_text())
check('assoc inputs: thresholds as configured', summ['thresholds'] == {'min_call_rate': C.MIN_CALL_RATE,
                                                                        'min_maf': C.MIN_MAF})
for co in C.COHORTS:
    ph = pd.read_csv(C.MODEL_INPUTS / co / 'phenotype/pheno.tsv', sep='\t', dtype={'IID': str})
    m = mk[mk.cohort == co].set_index('marker')
    check(f'assoc inputs {co}: every exported column assessed once',
          len(m) == ad.shape[1] + rd.shape[1] and m.index.is_unique)
    cr_low = (m.call_rate_cases < C.MIN_CALL_RATE) | (m.call_rate_controls < C.MIN_CALL_RATE)
    re_status = np.select([m.mac < 1, ~m.gene.isin(C.ASSOC_GENES), cr_low, m.maf < C.MIN_MAF],
                          ['not_varying', 'excluded_gene', 'call_rate', 'maf'], 'tested')
    check(f'assoc inputs {co}: status = first rule failed', (re_status == m.status.values).all())
    for kind, name in (('allele', 'allele_dosage.tsv'), ('residue', 'residue_dosage.tsv')):
        f = R / '06.assoc_inputs' / co / name
        d = pd.read_csv(f, sep='\t', index_col=0, dtype={'sample_id': str})
        want = list(m[(m['class'] == kind) & (m.status == 'tested')].index)
        check(f'assoc inputs {co} {kind}: samples = phenotype file, in its order', list(d.index) == ph.IID.tolist())
        check(f'assoc inputs {co} {kind}: columns = the tested markers', list(d.columns) == want, f'{d.shape[1]}')
        an = 2 * d.notna().sum(); ac = d.sum(); mac = np.minimum(ac, an - ac)
        cs = (ph.set_index('IID').PHENO1.reindex(d.index) == 1).values
        check(f'assoc inputs {co} {kind}: call rate (in cases and in controls) and MAF recomputed from the matrix pass',
              bool(((d[cs].notna().mean() >= C.MIN_CALL_RATE) & (d[~cs].notna().mean() >= C.MIN_CALL_RATE)
                    & (mac / an >= C.MIN_MAF)).all()))
        check(f'assoc inputs {co} {kind}: md5 as recorded',
              hashlib.md5(f.read_bytes()).hexdigest() == summ['cohorts'][co]['md5'][name])

# figures
figs = sorted(p.name for p in (R / 'figures').glob('*.png'))
check('figures: 8 present', len(figs) == 8, ', '.join(figs))

print(f'\n{len(fails)} check(s) failed' if fails else '\nall checks passed')
sys.exit(1 if fails else 0)
EOF
