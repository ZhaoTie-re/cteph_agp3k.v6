#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The independent re-check of the workbook's tables. Each record is
#           recomputed by another route from the per-sample or per-run files, with
#           code of its own (this script imports nothing from the table scripts):
#             typing   call states and confidence per gene from the genotype and
#                      confidence tables; every allele and residue column's copies,
#                      call rates, MAF and status from the dosage matrices and the
#                      case status; the stepwise record from those statuses
#             assoc.   samples from the phenotype files; every single-marker row
#                      against the per-run association tables (both models) and
#                      every position row against the per-comparison omnibus
#                      tables; thresholds as 0.05 / Me and 0.05 / tests; significance
#                      counts and the smallest P per scan from the detail rows; λGC
#                      with its own implementation; mirror pairs from the residue
#                      matrices; the DQB1 57 V genotypes from the matrices
#             cross    a number on two sheets is the same number
# Inputs  : --typing --assoc (the table folders) --typing-results --assoc-results
# Outputs : --out (recheck.json)
# Component: hla_records
# ---------------------------------------------------------------------------
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

COH = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']
CS = ['all', 'dvt_yes', 'dvt_no', 'dvt_case_only']
NOT_TESTED_GENES = {'DRB3', 'DRB4', 'DRB5'}


class Book:
    def __init__(self):
        self.checks = []

    def add(self, cid, group, sheet, what, route, pairs, rel=1e-6):
        """pairs: list of (item, workbook value, recomputed value)."""
        bad = [{'item': str(i), 'workbook': str(a), 'deck': str(b)} for i, a, b in pairs if not equal(a, b, rel)]
        self.checks.append({'id': cid, 'group': group, 'sheet': sheet, 'what': what, 'route': route,
                            'compared': len(pairs), 'mismatches': len(bad), 'passed': not bad and len(pairs) > 0,
                            'detail': bad[:50]})
        print(f'[recheck] {cid}: {len(pairs)} compared, {len(bad)} mismatches')


def isna(x):
    return x is None or (isinstance(x, float) and np.isnan(x)) or (isinstance(x, str) and x == '')


def equal(a, b, rel):
    if isna(a) or isna(b):
        return isna(a) and isna(b)
    if isinstance(a, (bool, np.bool_)) or isinstance(b, (bool, np.bool_)):
        return str(a).lower() == str(b).lower()
    try:
        x, y = float(a), float(b)
        return abs(x - y) <= rel * max(abs(x), abs(y)) + 1e-12
    except (TypeError, ValueError):
        return str(a) == str(b)


def read(p, **kw):
    return pd.read_csv(p, sep='\t', **kw)


def lam(p):
    p = np.asarray(p, dtype=float)
    p = p[np.isfinite(p)]
    return float(np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1))


# ── typing ──────────────────────────────────────────────────────────────────
def typing(B, T, TR):
    g = read(TR / '01.calls/genotypes.tsv', dtype=str)
    conf = read(TR / '05.exports/call_confidence.tsv', dtype=str)
    t3 = read(Path(T) / 't3_genes.tsv').set_index('gene')
    pairs = []
    for gene in t3.index:
        st = g[f'{gene}_state'].fillna('').value_counts()
        for k in ('heterozygous', 'homozygous', 'present', 'absent', 'ambiguous'):
            pairs.append((f'{gene} {k}', t3.at[gene, k], int(st.get(k, 0))))
        pairs.append((f'{gene} calls', t3.at[gene, 'calls'], len(g)))
        cc = conf[gene].fillna('').value_counts()
        for k in ('high', 'medium', 'low'):
            pairs.append((f'{gene} confidence {k}', t3.at[gene, f'conf_{k}'], int(cc.get(k, 0))))
    B.add('recheck|T3', 'Independent recheck', 'T3 Completeness', 'calls by result and by confidence, per gene',
          'counted from the per-sample genotype and confidence tables', pairs)

    # every exported column, per cohort: copies, call rates, MAF and status from the matrices and case status
    status_name = {'tested': 'Tested', 'maf': 'Below MAF 0.5%', 'call_rate': 'Call rate below 95%',
                   'excluded_gene': 'Gene not tested', 'not_varying': 'Not varying'}
    record = {}
    for kind, sheet, tbl in (('allele', 'T7 Alleles (detail)', 't7_alleles.tsv'),
                             ('residue', 'T8 Residues (detail)', 't8_residues.tsv')):
        X = read(TR / f'05.exports/{kind}_dosage.tsv', dtype={'sample_id': str}).set_index('sample_id')
        W = read(Path(T) / tbl, dtype={'marker': str, 'position': str, 'residue': str})
        W['col'] = W.marker if kind == 'allele' else W.gene + ':' + W.position + ':' + W.residue
        W = W.set_index('col')
        pairs = []
        for c, tag in zip(COH, 'nif'):
            ph = read(Path(B.assoc_results) / 'all' / c / '00.prep/pheno_covar.random.tsv', dtype={'IID': str})
            M = X.loc[ph.IID.values]
            case = ph.PHENO1.values == 1
            called = M.notna()
            n_case, n_ctrl = case.sum(), (~case).sum()
            cr_case = called[case].sum() / n_case
            cr_ctrl = called[~case].sum() / n_ctrl
            cp_case = M[case].sum()
            cp_ctrl = M[~case].sum()
            ncall = called.sum()
            copies = cp_case + cp_ctrl
            mac = np.minimum(copies, 2 * ncall - copies)
            maf = mac / (2 * ncall)
            nuniq = M.nunique()
            gene = pd.Series([x.split('*')[0] if kind == 'allele' else x.split(':')[0] for x in M.columns], index=M.columns)
            st = pd.Series('tested', index=M.columns)
            st[maf < 0.005] = 'maf'
            st[(cr_case < 0.95) | (cr_ctrl < 0.95)] = 'call_rate'
            st[gene.isin(NOT_TESTED_GENES)] = 'excluded_gene'
            st[nuniq <= 1] = 'not_varying'
            for col in M.columns:
                pairs += [(f'{c} {col} status', W.at[col, f'{tag}_status'], status_name[st[col]]),
                          (f'{c} {col} copies cases', W.at[col, f'{tag}_copies_cases'], cp_case[col]),
                          (f'{c} {col} copies controls', W.at[col, f'{tag}_copies_controls'], cp_ctrl[col]),
                          (f'{c} {col} call rate cases', W.at[col, f'{tag}_cr_cases'], cr_case[col]),
                          (f'{c} {col} call rate controls', W.at[col, f'{tag}_cr_controls'], cr_ctrl[col]),
                          (f'{c} {col} MAF', W.at[col, f'{tag}_maf'], round(float(maf[col]), 6))]
            # the stepwise record: columns left after each rule, applied in order
            v = st != 'not_varying'
            g2 = v & (st != 'excluded_gene')
            g3 = g2 & (st != 'call_rate')
            g4 = st == 'tested'
            record[(c, kind)] = {'typed': int(v.sum()), 'genes': int(g2.sum()), 'callrate': int(g3.sum()),
                                 'maf': int(g4.sum())}
        pairs = [(i, a, b) for i, a, b in pairs]
        B.add(f'recheck|{sheet[:2]}', 'Independent recheck', sheet,
              'every column in every cohort: status, copies in cases and controls, call rates, MAF',
              'recomputed from the dosage matrices of all genomes and the case status of the cohort (MAF to six '
              'decimals, as the typing component prints it)', pairs,
              rel=1e-5)
    rec = read(Path(T) / 't5_record.tsv')
    pairs = [(f'{r.cohort} {r.stage} {r._5}', r.columns, record[(r.cohort, r._5)][r.stage]) for r in rec.itertuples()]
    B.add('recheck|T5', 'Independent recheck', 'T5 Marker QC', 'the stepwise record per cohort and class',
          'counted from the statuses recomputed above', pairs)
    return record


# ── association ─────────────────────────────────────────────────────────────
def association(B, T, A, AR, TR, record):
    # samples
    s = read(Path(A) / 'a1_samples.tsv')
    pairs = []
    for r in s.itertuples():
        for m in ('random', 'fixed'):
            ph = read(AR / r.comparison / r.cohort / f'00.prep/pheno_covar.{m}.tsv')
            pairs += [(f'{r.comparison} {r.cohort} {m} samples', getattr(r, f'{m}_samples'), len(ph)),
                      (f'{r.comparison} {r.cohort} {m} cases', getattr(r, f'{m}_cases'), int((ph.PHENO1 == 1).sum())),
                      (f'{r.comparison} {r.cohort} {m} controls', getattr(r, f'{m}_controls'), int((ph.PHENO1 == 0).sum()))]
    B.add('recheck|A1', 'Independent recheck', 'A1 Design', 'samples of both models, per comparison and cohort',
          'counted from the phenotype files of each run', pairs)

    # single markers, cell by cell against the per-run tables
    D = read(Path(A) / 'a11_single_markers.tsv', dtype={'marker': str, 'position': str, 'residue': str})
    D['k'] = D.comparison + '|' + D.cohort + '|' + D.marker
    D = D.set_index('k')
    pairs = []
    n_rows = 0
    for cs in CS:
        for c in COH:
            for cls in ('allele', 'residue'):
                for m, pre in (('random', 'r_'), ('fixed', 'f_')):
                    u = read(AR / cs / c / f'03.tables/assoc.{m}.{cls}.tsv', dtype={'marker': str, 'position': str,
                                                                                      'residue': str})
                    for r in u.to_dict('records'):
                        r = pd.Series(r)
                        k = f'{cs}|{c}|{r.marker}'
                        if k not in D.index:
                            pairs.append((k, 'absent', 'present'))
                            continue
                        w = D.loc[k]
                        n_rows += m == 'random'
                        for f in ('or', 'l95', 'u95', 'p'):
                            pairs.append((f'{k} {pre}{f}', w[f'{pre}{f}'], r[f]))
                        if m == 'random':
                            pairs += [(f'{k} p_spa', w.r_p_spa, r.p_spa), (f'{k} copies case', w.copies_case, r.copies_case),
                                      (f'{k} copies control', w.copies_control, r.copies_control),
                                      (f'{k} n case', w.n_case, r.n_case), (f'{k} n control', w.n_control, r.n_control),
                                      (f'{k} freq case', w.freq_case, r.freq_case), (f'{k} mac', w.mac, r.mac),
                                      (f'{k} test', w.test_id, r.test_id)]
    pairs.append(('rows', len(D), n_rows))
    B.add('recheck|A11', 'Independent recheck', 'A11 Single markers (detail)',
          'every allele and residue row of both models (OR, CI, P, SPA P, counts, frequencies)',
          'compared cell by cell with the association tables of each of the 48 runs', pairs)

    # positions, cell by cell against the per-comparison omnibus tables
    P = read(Path(A) / 'a12_positions.tsv', dtype={'unit': str, 'position': str})
    P['k'] = P.comparison + '|' + P.cohort + '|' + P.unit
    P = P.set_index('k')
    pairs = []
    for cs in CS:
        u = read(AR / cs / '_omnibus/omnibus.residue.tsv', dtype={'unit': str, 'position': str})
        for r in u.itertuples():
            k = f'{cs}|{r.cohort}|{r.unit}'
            w = P.loc[k]
            for f in ('n', 'df', 'chisq', 'p', 'threshold', 'errcode'):
                pairs.append((f'{k} {f}', w[f], getattr(r, f)))
    pairs.append(('rows', len(P), sum(len(read(AR / cs / '_omnibus/omnibus.residue.tsv')) for cs in CS)))
    B.add('recheck|A12', 'Independent recheck', 'A12 Positions (detail)', 'every position row (n, df, χ², P, threshold)',
          'compared cell by cell with the omnibus table of each comparison', pairs)

    # thresholds: 0.05 / Me and 0.05 / tests from each run's effective-test table
    sc = read(Path(A) / 'a2_scans.tsv')
    pairs = []
    for r in sc.itertuples():
        e = read(AR / r.comparison / r.cohort / f'01.markers/{r._3}/effective_tests.{r._3}.tsv').iloc[0]
        pairs += [(f'{r.comparison} {r.cohort} {r._3} GEC', r.p_gec, 0.05 / e.me_gec),
                  (f'{r.comparison} {r.cohort} {r._3} Bonferroni', r.p_bonf, 0.05 / e.tests),
                  (f'{r.comparison} {r.cohort} {r._3} tests', r.tests, e.tests)]
    B.add('recheck|A2', 'Independent recheck', 'A2 Thresholds', 'GEC = 0.05 / Me and Bonferroni = 0.05 / tests per scan',
          'recomputed from Me and the test count of each run (Me is printed by GEC to six digits)', pairs, rel=1e-5)

    # significance counts, smallest P and λGC per scan, from the detail rows
    a3 = read(Path(A) / 'a3_scans.tsv')
    Dt = D.reset_index()
    Dt = Dt[Dt.test_id == Dt.marker]
    pairs = []
    for r in a3.itertuples():
        if r.test_type == 'position':
            d = P.reset_index()
            d = d[(d.comparison == r.comparison) & (d.cohort == r.cohort) & (d.errcode != 'NOT_TESTED')]
            ok = d[d.errcode == '.']
            pairs += [(f'{r.comparison} {r.cohort} position tests', r.tests, len(d)),
                      (f'{r.comparison} {r.cohort} position sig', r.sig_primary, int((d.p < d.threshold).sum())),
                      (f'{r.comparison} {r.cohort} position min p', r.min_p, ok.p.min()),
                      (f'{r.comparison} {r.cohort} position λ', r.lambda_gc, lam(ok.p))]
            continue
        d = Dt[(Dt.comparison == r.comparison) & (Dt.cohort == r.cohort) & (Dt['class'] == r.test_type)]
        pairs += [(f'{r.comparison} {r.cohort} {r.test_type} tests', r.tests, len(d)),
                  (f'{r.comparison} {r.cohort} {r.test_type} sig GEC', r.sig_primary, int((d.r_p < d.gec).sum())),
                  (f'{r.comparison} {r.cohort} {r.test_type} sig Bonferroni', r.sig_bonf, int((d.r_p < d.bonf).sum())),
                  (f'{r.comparison} {r.cohort} {r.test_type} sig SPA', r.sig_spa, int((d.r_p_spa < d.gec).sum())),
                  (f'{r.comparison} {r.cohort} {r.test_type} min p', r.min_p, d.r_p.min()),
                  (f'{r.comparison} {r.cohort} {r.test_type} top', r.top, d.loc[d.r_p.idxmin(), 'marker']),
                  (f'{r.comparison} {r.cohort} {r.test_type} λ', r.lambda_gc, lam(d.r_p))]
    B.add('recheck|A3', 'Independent recheck', 'A3 Q1 Association; A8 DVT',
          'per scan: tests, significant (primary, Bonferroni, SPA), smallest P and its marker, λGC',
          'counted from the detail rows (one row per test) with λGC from its own χ² implementation', pairs)

    # mirror pairs from the residue matrices (all cases vs. controls: the cohort's genomes)
    pairs = []
    for c in COH:
        R_ = read(TR / f'06.assoc_inputs/{c}/residue_dosage.tsv', dtype={'sample_id': str}).set_index('sample_id')
        ph = read(AR / 'all' / c / '00.prep/pheno_covar.random.tsv', dtype={'IID': str})
        R_ = R_.loc[ph.IID.values]
        pos = pd.Series([x.rsplit(':', 1)[0] for x in R_.columns], index=R_.columns)
        n = 0
        for p_, cols in pos.groupby(pos).groups.items():
            cols = list(cols)
            for i in range(len(cols)):
                for j in range(i + 1, len(cols)):
                    a, b = R_[cols[i]], R_[cols[j]]
                    both = a.notna() & b.notna()
                    if (a.isna() == b.isna()).all() and ((a[both] + b[both]) == 2).all():
                        n += 1
        w = sc[(sc.comparison == 'all') & (sc.cohort == c) & (sc['class'] == 'residue')].iloc[0]
        pairs.append((f'{c} mirror pairs', w.mirror_pairs, n))
    B.add('recheck|A10', 'Independent recheck', 'A2 Thresholds; A10 Appendix',
          'mirror pairs of residue columns, all cases vs. controls',
          'columns of one position that are complements (both undetermined, or summing to 2) in every genome', pairs)

    # DQB1 57 V copies per person
    gt = read(Path(A) / 'a7_genotypes.tsv')
    pairs = []
    for c in COH:
        R_ = read(TR / f'06.assoc_inputs/{c}/residue_dosage.tsv', dtype={'sample_id': str},
                  usecols=['sample_id', 'DQB1:57:V']).set_index('sample_id')
        ph = read(AR / 'all' / c / '00.prep/pheno_covar.random.tsv', dtype={'IID': str})
        v = R_.loc[ph.IID.values, 'DQB1:57:V']
        for grp, msk in (('cases', ph.PHENO1.values == 1), ('controls', ph.PHENO1.values == 0)):
            x = v[msk]
            w = gt[(gt.cohort == c) & (gt.group == grp)].iloc[0]
            pairs += [(f'{c} {grp} n', w.analyzed, len(x)), (f'{c} {grp} undetermined', w.undetermined, int(x.isna().sum())),
                      (f'{c} {grp} 0', w.n0, int((x == 0).sum())), (f'{c} {grp} 1', w.n1, int((x == 1).sum())),
                      (f'{c} {grp} 2', w.n2, int((x == 2).sum())), (f'{c} {grp} carriers', w.carriers, int((x > 0).sum())),
                      (f'{c} {grp} frequency', w.copy_freq, x.sum() / (2 * x.notna().sum()))]
    B.add('recheck|A7', 'Independent recheck', 'A7 Summary', 'copies of valine at DQB1 57 per person, cases and controls',
          'counted from the residue matrix of each cohort and the case status', pairs)

    # cross-sheet: one number, one value
    pairs = []
    t5 = read(Path(T) / 't5_record.tsv')
    for r in sc[sc.comparison == 'all'].itertuples():
        x = t5[(t5.cohort == r.cohort) & (t5['class'] == r._3) & (t5.stage == 'maf')].iloc[0]
        pairs.append((f'T5 kept = A2 columns {r.cohort} {r._3}', x.columns, r.markers))
    a3i = a3.set_index(['comparison', 'cohort', 'test_type'])
    for r in sc.itertuples():
        pairs.append((f'A2 = A3 tests {r.comparison} {r.cohort} {r._3}', r.tests, a3i.loc[(r.comparison, r.cohort, r._3), 'tests']))
    pz = read(Path(A) / 'a2_positions.tsv')
    for r in pz.itertuples():
        pairs.append((f'A2 = A3 positions {r.comparison} {r.cohort}', r.tested, a3i.loc[(r.comparison, r.cohort, 'position'), 'tests']))
    tg = read(Path(A) / 'a1_tests_by_gene.tsv')
    for c in COH:
        for cls in ('allele', 'residue'):
            pairs.append((f'A1 sum = A2 tests {c} {cls}', int(tg[tg.gene != 'All'][f'{cls}_{c}'].sum()),
                          a3i.loc[('all', c, cls), 'tests']))
        pairs.append((f'A1 sum = A3 positions {c}', int(tg[tg.gene != 'All'][f'position_{c}'].sum()),
                      a3i.loc[('all', c, 'position'), 'tests']))
    v = read(Path(A) / 'a7_v57.tsv').set_index('cohort')
    st = read(Path(A) / 'a8_strata.tsv')
    for r in st.itertuples():
        w = D.loc[f'{r.comparison}|{r.cohort}|DQB1:57:V']
        pairs += [(f'A8 = A11 {r.comparison} {r.cohort} p', r.p, w.r_p), (f'A8 = A11 {r.comparison} {r.cohort} OR', r._7, w.r_or)]
        if r.comparison == 'all':
            pairs += [(f'A7 = A11 {r.cohort} p', v.at[r.cohort, 'p'], w.r_p),
                      (f'A7 = A11 {r.cohort} fixed p', v.at[r.cohort, 'p_fixed'], w.f_p)]
    for r in read(Path(A) / 'a7_v57.tsv').itertuples():
        w = P.loc[f'all|{r.cohort}|DQB1:57']
        pairs.append((f'A7 = A12 {r.cohort} position 57', r.p_pos57, w.p))
    B.add('cross|all', 'Cross-sheet', 'T5, A1–A3, A7, A8, A11, A12',
          'a number shown on two sheets: columns kept, tests per scan, tests per gene, DQB1 57 V and position 57',
          'the two sheets\' tables read side by side', pairs)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for k in ('typing', 'assoc', 'typing_results', 'assoc_results', 'out'):
        ap.add_argument('--' + k.replace('_', '-'), required=True)
    a = ap.parse_args()
    B = Book()
    B.assoc_results = Path(a.assoc_results)
    record = typing(B, a.typing, Path(a.typing_results))
    association(B, a.typing, a.assoc, Path(a.assoc_results), Path(a.typing_results), record)
    out = {'all_passed': all(x['passed'] for x in B.checks), 'checks': B.checks}
    Path(a.out).write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str) + '\n')
    print(f"[recheck] {sum(x['compared'] for x in B.checks)} values compared, "
          f"{sum(x['mismatches'] for x in B.checks)} mismatches")
    for x in B.checks:
        for d in x['detail'][:6]:
            print(f"  {x['sheet']}: {d['item']}: workbook {d['workbook']} | recomputed {d['deck']}")


if __name__ == '__main__':
    main()
