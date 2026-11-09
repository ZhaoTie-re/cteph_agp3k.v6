#!/usr/bin/env python3
"""Read-only checks of check/platform_depth (run by verify.sh). Exit 1 on any failure.

Nothing is written except to stdout. Each check prints PASS / FAIL with its numbers. Checks 13-16 re-derive
evidence-matrix readings from the source tables with code of their own (not evidence_matrix.py's or
report_lib.py's), and hold the registry (scripts/checks.py), the workflow, the docs and the published rules to
each other.
"""
import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / 'results'
FAIL = []
sys.path.insert(0, str(ROOT / 'scripts'))
import checks as CK  # noqa: E402  the registry (data only)
import frozen  # noqa: E402  the frozen pieces of the cached stages (read-only check)
COHORTS = ('narrow_mainland', 'intermediate_mainland', 'full_mainland')


def check(name, ok, detail=''):
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")
    if not ok:
        FAIL.append(name)


def md5(p):
    h = hashlib.md5()
    with open(p, 'rb') as fh:
        for c in iter(lambda: fh.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--label', default='downsampling')
    a = ap.parse_args()

    # 1. configuration: resources only through --rsc, memory = cores x 4571M
    cfg = (ROOT / 'nextflow.config').read_text()
    check('no cpus/memory directives', not re.search(r'^\s*(cpus|memory)\s*=', cfg, re.M))
    bad = [m for m in re.findall(r"--rsc p=1:t=(\d+):c=(\d+):m=(\d+)M", cfg) if int(m[2]) != int(m[1]) * 4571]
    check('--rsc memory = cores x 4571M', not bad, str(bad))

    # 2. reused inputs unchanged since the manifest was taken
    man = ROOT / 'docs' / 'INPUTS.md5'
    if man.exists():
        diffs = []
        for ln in man.read_text().splitlines():
            h, p = ln.split('  ', 1)
            if not Path(p).exists() or md5(p) != h:
                diffs.append(p)
        check('reused inputs unchanged (docs/INPUTS.md5)', not diffs, f'{len(diffs)} changed: {diffs[:3]}')

    # 3. samples, depth, fractions
    st = pd.read_csv(RES / '00.samples' / 'sample_table.tsv', sep='\t', dtype=str)
    full = st[st.in_full_mainland == '1']
    check('full mainland 439 cases / 2,662 controls', (full.status == 'case').sum() == 439 and (full.status == 'control').sum() == 2662)
    fr = pd.read_csv(RES / '01.depth' / 'fractions.tsv', sep='\t')
    check('89 deep cases, 0 < f < 1', len(fr) == 89 and fr.fraction.between(0, 1, inclusive='neither').all(),
          f'f {fr.fraction.min():.3f}-{fr.fraction.max():.3f}')
    fa = pd.read_csv(RES / '01.depth' / 'fractions.all.tsv', sep='\t', dtype={'sample_id': str})
    t7 = full[(full.status == 'case') & (full.platform == 'DNBSeq-T7 30x')]
    dep_ = pd.read_csv(RES / '01.depth' / 'depth_per_sample.tsv', sep='\t', dtype={'sample_id': str})
    dstar = dep_[(dep_.in_full_mainland == 1) & (dep_.status == 'control')].probe_depth.median()
    t7_med = dep_[dep_.sample_id.isin(t7.sample_id)].probe_depth.median()
    ft7 = fa[fa.platform == 'DNBSeq-T7 30x'].fraction
    check('down-sampled cases: the 89 deep cases as fractions.tsv + every full-mainland DNBSeq-T7 case at one '
          'fraction D* / their median',
          len(fa) == 89 + len(t7) and set(fa.sample_id) == set(fr.sample_id.astype(str)) | set(t7.sample_id) and
          ft7.nunique() == 1 and abs(ft7.iloc[0] - round(dstar / t7_med, 5)) < 1e-9,
          f'{len(fa)} cases; DNBSeq-T7 f = {ft7.iloc[0]} (D* {dstar:.4f} / {t7_med:.4f})')
    dep = pd.read_csv(RES / '01.depth' / 'depth_per_sample.tsv', sep='\t', dtype={'sample_id': str})
    check('depth measured for every genome', set(dep.sample_id) == set(st.sample_id), f'{len(dep)} genomes')

    # 4. gates
    gates = sorted((RES / '_gates').glob('G[1-5]*.json')) + sorted((RES / '_gates' / a.label).glob('G*.json'))
    for g in gates:
        d = json.loads(g.read_text())
        check(f'gate {g.relative_to(RES)}', bool(d.get('pass')))
    need = {'G1': 3, 'G2': 6, 'G3': 1, 'G4': 1, 'G5': 3}
    for k, n in need.items():
        got = len(list((RES / '_gates').glob(f'{k}*.json')))
        check(f'{k} present for every cohort/set', got == n, f'{got} of {n}')

    # 5. platform subsets: strata add up, controls shared
    for c in ('narrow_mainland', 'intermediate_mainland', 'full_mainland'):
        f = RES / '08.subsets' / c / 'subsets' / 'subset_counts.tsv'
        if not f.exists():
            check(f'{c}: subset counts present', False, str(f.relative_to(RES)))
            continue
        s = pd.read_csv(f, sep='\t').set_index('subset')
        strata = s[s.index.str.startswith('only_')]
        lopo = s[s.index.str.startswith('without_')]
        ok = (strata.cases_random.sum() == s.loc['pooled', 'cases_random'] and s.controls_random.nunique() == 1 and
              all(lopo.loc[f'without_{x[5:]}', 'cases_random'] == s.loc['pooled', 'cases_random'] - strata.loc[x, 'cases_random']
                  for x in strata.index))
        check(f'{c}: strata sum to pooled, controls shared, LOPO = pooled - stratum', ok)

    # 6. heterogeneity recomputed independently for the GWAS leads and DQB1:57:V
    for c in ('narrow_mainland', 'intermediate_mainland', 'full_mainland'):
        hf = RES / '10.platform' / c / 'platform.heterogeneity.tsv'
        if not hf.exists():
            check(f'{c}: platform heterogeneity present', False, str(hf.relative_to(RES)))
            continue
        het = pd.read_csv(hf, sep='\t', dtype={'variant_id': str}).set_index('variant_id')
        strt = pd.read_csv(RES / '10.platform' / c / 'platform.strata.tsv', sep='\t', dtype={'variant_id': str})
        cnt = pd.read_csv(RES / '08.subsets' / c / 'subsets' / 'subset_counts.tsv', sep='\t').set_index('subset')
        n0 = cnt.loc['pooled', 'controls_random']
        worst = 0.0
        for vid in list(het.index[:25]) + (['DQB1:57:V'] if 'DQB1:57:V' in het.index else []):
            x = strt[(strt.variant_id == vid) & strt.subset.str.startswith('only_')].set_index('subset')
            x = x[np.isfinite(x.beta) & (x.se > 0)]
            if len(x) < 2:
                continue
            nk = np.array([cnt.loc[k, 'cases_random'] for k in x.index], float)
            R = np.sqrt(np.outer(nk, nk) / np.outer(nk + n0, nk + n0))
            np.fill_diagonal(R, 1.0)
            S = R * np.outer(x.se, x.se)
            Si = np.linalg.inv(S)
            one = np.ones(len(x))
            mu = one @ Si @ x.beta.values / (one @ Si @ one)
            q = (x.beta.values - mu) @ Si @ (x.beta.values - mu)
            worst = max(worst, abs(q - het.loc[vid, 'Q']) / max(abs(het.loc[vid, 'Q']), 1e-9),
                        abs(stats.chi2.sf(q, len(x) - 1) - het.loc[vid, 'p_Q']))
        check(f'{c}: Q and its P recomputed independently', worst < 1e-6, f'max deviation {worst:.2e}')

    # 7. array evidence: gate G10 passed; one row per signal in every cohort (the verdicts are results, printed)
    gf = RES / '11.array' / 'array.gate.json'
    g10 = json.loads(gf.read_text()) if gf.exists() else {'pass': False, 'tables': {}}
    check('gate G10: the array re-run in full mainland reproduces check/array_crosscheck',
          g10['pass'] and len(g10['tables']) == 7 and all(t['pass'] for t in g10['tables'].values()),
          f'{len(g10["tables"])} tables')
    af = RES / '11.array' / 'array.signals.tsv'
    if not af.exists():
        check('array evidence present', False, str(af.relative_to(RES)))
    else:
        s = pd.read_csv(af, sep='\t')
        sig = pd.read_csv(RES / '02.targets' / 'signal_variants.tsv', sep='\t', dtype=str)
        want = set(sig[(sig.kind == 'gwas') & (sig.role == 'Lead')].unit) | {'DQB1:57:V'}
        genes = {'STBD1', 'OSBP', 'RPH3A', 'JAK2'}                       # gene carriers: full mainland only
        ok = all(set(s[s.cohort == c].signal) == (want | genes if c == 'full_mainland' else want) and
                 (s.cohort == c).sum() == len(want) + (len(genes) if c == 'full_mainland' else 0) for c in COHORTS)
        check('array evidence: one row per signal in every cohort (genes in full mainland)', ok and set(s.cohort) == set(COHORTS),
              f'{len(s)} rows; ' + ', '.join(f'{k} {v}' for k, v in s.array_verdict.value_counts().items()))

    # 8. control sources: the subsets add up; one heterogeneity row per headline marker
    for c in ('narrow_mainland', 'intermediate_mainland', 'full_mainland'):
        f = RES / '08.subsets' / c / 'control_sources' / 'subsets' / 'subset_counts.tsv'
        if not f.exists():
            check(f'{c}: control-source subsets present', False, str(f.relative_to(RES)))
            continue
        s = pd.read_csv(f, sep='\t').set_index('subset')
        pooled = pd.read_csv(RES / '08.subsets' / c / 'subsets' / 'subset_counts.tsv', sep='\t').set_index('subset')
        only = s[s.index.str.startswith('ctrlonly_')]
        wo = s[s.index.str.startswith('ctrlwithout_')]
        n0 = pooled.loc['pooled', 'controls_random']
        ok = (only.controls_random.sum() == n0 and s.cases_random.eq(pooled.loc['pooled', 'cases_random']).all() and
              all(wo.loc[f'ctrlwithout_{x[9:]}', 'controls_random'] == n0 - only.loc[x, 'controls_random'] for x in only.index))
        check(f'{c}: control sources sum to the controls, cases in every subset', ok, f'{len(only)} sources')
    cs = RES / '12.controls' / 'controls.signals.tsv'
    if not cs.exists():
        check('control-source summary present', False, str(cs.relative_to(RES)))
    else:
        s = pd.read_csv(cs, sep='\t')
        sig = pd.read_csv(RES / '02.targets' / 'signal_variants.tsv', sep='\t', dtype=str)
        n_lead = sig[(sig.kind == 'gwas') & (sig.role == 'Lead')].unit.nunique() + 1
        check('control sources: every lead and DQB1:57:V in every cohort', len(s) == 3 * n_lead, f'{len(s)} rows')

    # 9. follow-up: the re-applied genotype filters reproduce the published genotypes (else the audit of
    #    carriers set missing does not stand); every lead and every tiered gene set present, in every cohort
    ff = RES / '13.followup' / 'followup.gene_filters.tsv'
    if not ff.exists():
        check('follow-up present', False, str(ff.relative_to(RES)))
    else:
        f = pd.read_csv(ff, sep='\t')
        check('follow-up: re-applied filters = published genotypes at every gene variant, every cohort',
              len(f) > 0 and set(f.cohort) == set(COHORTS) and (f.genotype_agreement == 1).all() and
              (f.missing_published == f.missing_reapplied).all(), f'{len(f)} cohort x variant rows')
        fa = pd.read_csv(RES / '13.followup' / 'followup.gwas_af.tsv', sep='\t')
        sig = pd.read_csv(RES / '02.targets' / 'signal_variants.tsv', sep='\t', dtype=str)
        leads = set(sig[(sig.kind == 'gwas') & (sig.role == 'Lead')].unit)
        check('follow-up: every GWAS lead has group frequencies in every cohort',
              all(set(fa[fa.cohort == c].signal) == leads for c in COHORTS), f'{fa.signal.nunique()} leads')

    # 10. evidence matrix: one row per signal, every check read, readings from the fixed vocabulary
    vocab = {'consistent', 'flag', 'not testable', 'not applicable', 'not checkable', 'reported', 'genome-wide', 'suggestive', 'below 1e-5',
             'significant', 'not significant'}
    ef = RES / '14.evidence' / 'evidence.signals.tsv'
    if not ef.exists():
        check('evidence matrix present', False, str(ef.relative_to(RES)))
    else:
        e = pd.read_csv(ef, sep='\t')
        lg = pd.read_csv(RES / '14.evidence' / 'evidence.long.tsv', sep='\t')
        rules = pd.read_csv(RES / '14.evidence' / 'evidence.rules.tsv', sep='\t')
        sig = pd.read_csv(RES / '02.targets' / 'signal_variants.tsv', sep='\t', dtype=str)
        n_want = (sig[(sig.kind == 'gwas') & (sig.role == 'Lead')].unit.nunique() + 2 +
                  len(sig[(sig.kind == 'gene') & sig.stratum.isin(CK.GENE_SETS_READ)][['unit', 'stratum']].drop_duplicates()))
        check('evidence matrix: one row per signal, every rule used, readings in the vocabulary',
              len(e) == n_want and e.signal.is_unique and set(rules.check) == set(lg.check) and set(lg.reading) <= vocab,
              f'{len(e)} signals; flagged: ' + '; '.join(f'{r.signal} ({r.flags})' for r in e.itertuples() if r.flags != 'none'))

    # 11. figures: one README per family
    fam = [d for d in (RES / 'figures').glob('0*') if d.is_dir()] if (RES / 'figures').exists() else []
    check('one README per figure family', fam and all((d / 'README.md').exists() for d in fam), f'{len(fam)} families')

    # 12. gate G7: every swap left the untouched genomes byte-identical and swapped every down-sampled case the
    #     cohort's fileset holds
    deep = set(pd.read_csv(RES / '03.calls' / a.label / 'deep_cases.tsv', sep='\t', dtype=str).sample_id)
    sw = sorted((RES / '04.genotypes' / a.label).glob('*/*/swap.*.json'))
    arms = {f.parts[-3] for f in sw}
    ok = len(sw) == 2 * len(COHORTS) * len(arms)
    for f in sw:
        r = json.loads(f.read_text())
        kind, cohort = f.name.split('.')[1], f.parts[-2]
        fam = (RES / '04.genotypes' / 'original' / cohort / f'{kind}.fam').read_text().split('\n')
        held = sum(1 for ln in fam if ln.strip() and ln.split()[1] in deep)
        ok &= r['untouched_identical'] and r['samples_swapped'] == held
    check('gate G7: swaps leave untouched genomes identical, every down-sampled case of the fileset swapped', ok,
          f'{len(sw)} swap reports, {len(arms)} arms, {len(deep)} down-sampled cases')

    # 13. evidence-matrix readings re-derived from the source tables (own code)
    ev = pd.read_csv(RES / '14.evidence' / 'evidence.long.tsv', sep='\t')
    evs = pd.read_csv(RES / '14.evidence' / 'evidence.signals.tsv', sep='\t', keep_default_na=False)
    sig = pd.read_csv(RES / '02.targets' / 'signal_variants.tsv', sep='\t', dtype=str)
    lead = sig[(sig.kind == 'gwas') & (sig.role == 'Lead')].drop_duplicates('unit').set_index('unit').variant_id
    hla57 = 'DQB1:57:V'
    n_saige = len(lead) + 1

    def reading(s_, chk, c):
        x = ev[(ev.signal == s_) & (ev.check == chk) & (ev.cohort == c)]
        return (x.reading.iloc[0], x.value.iloc[0]) if len(x) == 1 else (None, np.nan)

    def loo(z_obs, z_exp):
        z_obs, z_exp = np.asarray(z_obs, float), np.asarray(z_exp, float)
        ok = np.isfinite(z_obs) & np.isfinite(z_exp) & (np.abs(z_exp) >= 2.8)
        if not ok.any():
            return 'not testable'
        weak = (np.abs(z_obs) < 1.96) | (z_obs * z_exp <= 0)
        return 'flag' if (ok & weak).any() else 'consistent'

    dsg = pd.read_csv(RES / '07.downsampling' / a.label / 'summary' / 'ds.gwas.tsv', sep='\t')
    dsg = dsg[(dsg.role == 'Lead') & (dsg.imputation == 'best_guess')]
    reps = [c for c in dsg.columns if re.fullmatch(r'p_r\d+', c)]
    cl = pd.read_csv(RES / '12.controls' / 'controls.loo.tsv', sep='\t', dtype={'variant_id': str})
    n = {'depth': 0, 'platform_Q': 0, 'lopo': 0, 'sources_loo': 0}
    bad = []
    for c in COHORTS:
        het = pd.read_csv(RES / '10.platform' / c / 'platform.heterogeneity.tsv', sep='\t',
                          dtype={'variant_id': str}).set_index('variant_id')
        lo = pd.read_csv(RES / '10.platform' / c / 'platform.lopo.tsv', sep='\t', dtype={'variant_id': str})
        for unit, vid in list(lead.items()) + [(hla57, hla57)]:
            if unit != hla57:
                d = dsg[(dsg.unit == unit) & (dsg.cohort == c)].iloc[0]
                dz = stats.norm.isf(d.p_original / 2) - stats.norm.isf(np.median(d[reps].astype(float)) / 2)
                got, val = reading(unit, 'depth', c)
                n['depth'] += 1
                if got != ('flag' if dz > 1 else 'consistent') or not abs(val - dz) < 1e-9:
                    bad.append(f'depth {unit} {c}')
            h = het.loc[vid]
            maf = min(h.af_pooled, 1 - h.af_pooled)
            want = 'not applicable' if maf < 0.05 else 'flag' if h.p_Q_empirical < 0.05 / n_saige else 'consistent'
            n['platform_Q'] += 1
            if reading(unit, 'platform_Q', c)[0] != want:
                bad.append(f'platform_Q {unit} {c}')
            for chk, t in (('lopo', lo), ('sources_loo', cl[cl.cohort == c])):
                x = t[t.variant_id == vid]
                n[chk] += 1
                if reading(unit, chk, c)[0] != loo(x.z_obs, x.z_expected):
                    bad.append(f'{chk} {unit} {c}')
    check('evidence matrix: depth, platform_Q, lopo, sources_loo of the SAIGE signals re-derived', not bad,
          ', '.join(f'{k} {v}' for k, v in n.items()) + (f'; differ: {bad[:4]}' if bad else ''))

    arr = pd.read_csv(RES / '11.array' / 'array.signals.tsv', sep='\t')
    to_reading = {'consistent on arrays': 'consistent', 'carriers counted': 'reported', 'not checkable': 'not checkable'}
    bad, k = [], 0
    for r_ in arr.itertuples():
        rows = ev[(ev.check == 'array') & (ev.cohort == r_.cohort) &
                  ((ev.signal == r_.signal) | ev.signal.str.startswith(f'{r_.signal} ('))]
        k += 1
        want = 'flag' if 'contradict' in r_.array_verdict else to_reading.get(r_.array_verdict, r_.array_verdict)
        if not len(rows) or (rows.reading != want).any():
            bad.append(f'array {r_.signal} {r_.cohort}')
    ga = ev[(ev.check == 'array') & (ev.kind == 'gene')]
    if not (len(ga) and set(ga.cohort) == {'full_mainland'} and ga.signal.is_unique and
            set(ga.reading) <= {'reported', 'not checkable'}):
        bad.append('array: the gene sets need one cell each (full mainland, reported or not checkable)')
    check('evidence matrix: array readings re-derived in every cohort', not bad,
          f'{k} readings' + (f'; differ: {bad[:4]}' if bad else ''))
    one = (ev.check == 'array') & (ev.kind == 'gene')            # the gene sets' arrays: one cell (nested cohorts)
    check('evidence matrix: every check read in all three cohorts (the gene sets\' arrays: one cell)',
          all(set(x.cohort) == set(COHORTS) for _, x in ev[(ev.check != 'discovery') & ~one].groupby(['signal', 'check'])),
          f'{ev[["signal", "check"]].drop_duplicates().shape[0]} signal x check pairs')

    bad = []
    for r in evs.itertuples():
        for chk in [c for c in CK.RULES if c != 'discovery']:
            by = ev[(ev.signal == r.signal) & (ev.check == chk)].set_index('cohort').reading
            by = by.reindex([c for c in COHORTS if c in by.index])
            if not len(by):
                want = ''
            else:
                nf = int((by == 'flag').sum())
                want = 'flag' if nf >= CK.min_cohorts(chk) else 'flag in 1/3' if nf == 1 else by.iloc[0] if by.nunique() == 1 else \
                    'consistent' if (by == 'consistent').any() else ' / '.join(by)
            if getattr(r, chk) != want:
                bad.append(f'{r.signal}:{chk}')
        fl = [chk for chk in CK.RULES if getattr(r, chk) == 'flag']
        if r.flags != (', '.join(fl) if fl else 'none'):
            bad.append(f'{r.signal}:flags')
    check('evidence matrix: per-signal readings = the 2-of-3 rule over the cohorts', not bad,
          f'{len(evs)} signals' + (f'; differ: {bad[:4]}' if bad else ''))

    # 14. the registry against the published rules, the matrix columns, the workflow and docs/PIPELINE.md
    rules = pd.read_csv(RES / '14.evidence' / 'evidence.rules.tsv', sep='\t')
    check('registry: published rules = checks.RULES, in column order',
          list(rules.check) == list(CK.RULES) and list(rules.rule) == list(CK.RULES.values()), f'{len(rules)} rules')
    check('registry: evidence.signals.tsv columns = the registry\'s checks',
          list(evs.columns) == ['signal', 'kind'] + list(CK.RULES) + ['flags'])
    nf = (ROOT / 'platform_depth.nf').read_text()
    procs = re.findall(r'^process (\w+) \{', nf, re.M)
    uncached = re.findall(r'^process (\w+) \{\n\s+cache false', nf, re.M)
    check('registry: every workflow process listed once; report processes = the `cache false` ones (but ARRAY_COHORTS)',
          sorted(procs) == sorted(p.name for p in CK.PROCESSES) and len(procs) == len(set(procs)) and
          sorted(uncached) == sorted(set(CK.REPORT_PROCESSES) - set(CK.REPORT_CACHED)) and
          sorted(CK.REPORT_PROCESSES) == sorted(frozen.REPORT), f'{len(procs)} processes')
    check('docs/PIPELINE.md written from the registry (python3 scripts/checks.py markdown)',
          (ROOT / 'docs' / 'PIPELINE.md').exists() and (ROOT / 'docs' / 'PIPELINE.md').read_text() == CK.markdown())
    check('registry: reading vocabulary = verify\'s', set(CK.READINGS) | set(CK.DISCOVERY_READINGS) == vocab)

    # 15. table conventions of the report tables: a header of unique names, missing = empty cell, known cohorts
    tabs = sorted((RES / '07.downsampling' / a.label / 'summary').glob('*.tsv'))
    tabs += sorted(f for d in ('10.platform', '11.array', '12.controls', '13.followup', '14.evidence')
                   for f in (RES / d).rglob('*.tsv'))
    bad = []
    for f in tabs:
        t = pd.read_csv(f, sep='\t', dtype=str, keep_default_na=False)
        if (t.columns.str.startswith('Unnamed').any() or not t.columns.is_unique or
                t.isin(['nan', 'NaN', 'inf', '-inf', 'None', 'NA']).any().any() or
                ('cohort' in t.columns and not set(t.cohort) <= set(COHORTS))):
            bad.append(str(f.relative_to(RES)))
    check('report tables: header, missing values empty, cohort names known', not bad, f'{len(tabs)} tables; {bad[:3]}')

    # 16. the frozen pieces of the cached stages unchanged (docs/FROZEN.manifest.tsv)
    want = {ln.split('\t')[1]: ln.split('\t')[2:] for ln in frozen.MANIFEST.read_text().splitlines()[1:]}
    have = {r[1]: list(r[2:]) for r in frozen.current()}
    diff = sorted(k for k in set(want) | set(have) if want.get(k) != have.get(k))
    check('frozen pieces unchanged (docs/FROZEN.manifest.tsv)', not diff, f'{len(want)} items; changed: {diff[:3]}')

    # 17. the names kept twice agree: the workflow's params, common.py, control_subsets.py, the signal table
    import ast
    import common as K
    from control_subsets import CONTROL_SOURCES
    par = {m.group(1): ast.literal_eval(m.group(2)) for m in
           re.finditer(r"^params\.(Cohorts|Strata|CasePlatformSlugs|ControlSources)\s*=\s*(\[[^\]]*\])", nf, re.M)}
    want = {'Cohorts': list(K.COHORTS), 'Strata': sorted(set(sig.loc[sig.kind == 'gene', 'stratum']), reverse=True),
            'CasePlatformSlugs': [K.PLATFORM_SLUG[p] for p in K.CASE_PLATFORMS],
            'ControlSources': list(CONTROL_SOURCES.values())}
    import report_lib as RL
    check('names agree: workflow params = common.py, control_subsets.py, report_lib.py, the signal table',
          par == want and RL.VARIANT_SETS == want['Strata'] and list(RL.SOURCE_LABEL) == want['ControlSources'] and
          set(RL.PLATFORM_SHORT) == set(K.PLATFORM_SLUG) and list(RL.COHORT_SHORT) == want['Cohorts'],
          '; '.join(f'{k} {v}' for k, v in par.items()))

    print(f'\n{len(FAIL)} failed' if FAIL else '\nALL PASSED')
    sys.exit(1 if FAIL else 0)


if __name__ == '__main__':
    main()
