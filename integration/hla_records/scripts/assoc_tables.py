#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The Association part of the workbook (sheets A1-A12), assembled from
#           the finished analysis/assoc_hla.v2 run, with the association deck's
#           definitions (counts over tests, a mirror pair once; lambda_GC over
#           tests, each P read as a 1-df chi-square; the random model for alleles
#           and residues, the fixed-model omnibus for positions; GEC and Bonferroni
#           for single markers, Bonferroni for positions; the exact max-T as the
#           benchmark only). The pipeline's multivariate-normal max-T is not used.
#             A1 design            samples per comparison, tests per gene, methods
#             A2 thresholds        per scan: tests, Me, thresholds; exact max-T;
#                                  GEC recomputed from the genotypes
#             A3 Q1 association    per scan: what passes; the strongest markers
#             A4 Q2 conditioning   stepwise walks; what is left
#             A5 Q3 GWAS lead      the lead given HLA markers and the reverse
#             A6 inflation         lambda_GC before / after conditioning
#             A7 summary           DQB1 57 V and position 57 in full
#             A8 DVT               the three DVT comparisons; without DVT given 57 V
#             A9 reported alleles  B*52:01, DPB1*02:02
#             A10 appendix         mirror pairs, fixed vs. random, platform,
#                                  position 57 by P group against the panels
#             A11 / A12            every single-marker / position test (detail)
# Inputs  : --assoc --typing --model-inputs --cram --saige --given57v --plink2
#           --gwas-lambda (gwas_lambda.py output)
# Outputs : --out/a*.tsv, --out/assoc_facts.json
# Component: hla_records
# ---------------------------------------------------------------------------
import argparse
import glob
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hr_common import (COHORTS, COMPARISONS, COMP_LABEL, COMP_SIDES, DPB, LEAD57, LEAD87, POS125,  # noqa: E402
                       die, log, write_json, write_tsv)

POS57_SET = 'DQB1:57:A;DQB1:57:S;DQB1:57:V'
V57_ALLELES = ['DQB1*05:01', 'DQB1*06:04', 'DQB1*06:09']
A0503 = 'DQB1*05:03'
NODVT_NAMED = ['DQB1*05:01', 'DQB1*06:04', 'DQB1*06:09', 'DQB1*05:03', 'DQB1*05:02']
K_COMPONENTS = {'narrow_mainland': 9, 'intermediate_mainland': 12, 'full_mainland': 17}
# as printed in the two papers (the deck's transcription from the full texts): Tanabe 2005 Table 2 and Results;
# Kominami 2009 Table 2 (carrier counts recovered from its percentages; they reproduce its ORs)
PUBLISHED = [
    dict(study='Tanabe 2005', set='CTEPH with or without DVT', allele='B*52:01', case_car=32, case_n=80, ctrl_car=161,
         ctrl_n=678, or_=2.14, l95=1.29, u95=3.55, p=0.0016, pc=0.037),
    dict(study='Tanabe 2005', set='CTEPH with or without DVT', allele='DPB1*02:02', case_car=15, case_n=80, ctrl_car=43,
         ctrl_n=678, or_=3.41, l95=1.71, u95=6.74, p=0.00008, pc=0.0014),
    dict(study='Kominami 2009', set='CTEPH without DVT', allele='B*52:01', case_car=45, case_n=99, ctrl_car=96,
         ctrl_n=380, or_=2.47, l95=1.56, u95=3.90, p=8.6e-5, pc=0.016),
    dict(study='Kominami 2009', set='CTEPH without DVT', allele='DPB1*02:02', case_car=19, case_n=99, ctrl_car=17,
         ctrl_n=380, or_=5.07, l95=2.52, u95=10.19, p=7.5e-7, pc=0.00014)]
REP_ALLELES = ['B*52:01', 'DPB1*02:02']


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assoc', required=True, help='analysis/assoc_hla.v2 directory')
    p.add_argument('--typing', required=True, help='hla.typing.v2 directory')
    p.add_argument('--model-inputs', required=True)
    p.add_argument('--cram', required=True, help='cram.v6.summary.csv (WGS_Platform)')
    p.add_argument('--saige', required=True, help='analysis/assoc_saige directory')
    p.add_argument('--given57v', required=True, help='the supplementary without-DVT run given DQB1:57:V')
    p.add_argument('--gwas-lambda', required=True, nargs='+', help='gwas_lambda.<cohort>.json of this component')
    p.add_argument('--plink2', required=True)
    p.add_argument('--out', required=True)
    return p.parse_args()


def lam_gc(p):
    """lambda_GC from P values: median chi-square (each P read as 1 df) over its expectation (the decks')."""
    p = np.asarray(p, dtype=float)
    p = p[np.isfinite(p)]
    return float(np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1)) if len(p) else np.nan


class Assoc:
    def __init__(self, a):
        self.a = a
        self.R = Path(a.assoc) / 'results'
        self.TY = Path(a.typing)
        sys.path.insert(0, str(self.TY / 'scripts'))
        import config as TC            # noqa: E402
        import nomenclature as TN      # noqa: E402
        import panels as TP            # noqa: E402
        self.TC, self.TN, self.TP = TC, TN, TP
        R = self.R
        self.A = pd.concat([pd.read_csv(R / cs / '_comparison/assoc_all.tsv', sep='\t', dtype={'marker': str})
                            for cs in COMPARISONS], ignore_index=True)
        self.O = pd.read_csv(R / '_omnibus/omnibus_all.tsv', sep='\t', dtype={'unit': str, 'position': str})
        self.thr = pd.concat([pd.read_csv(R / cs / c / '01.markers' / k / f'effective_tests.{k}.tsv', sep='\t')
                              for cs in COMPARISONS for c in COHORTS for k in ('allele', 'residue')], ignore_index=True)
        self.boot = pd.concat([pd.read_csv(R / cs / c / '01.markers' / k / f'maxt_boot.{k}.tsv', sep='\t')
                               for cs in COMPARISONS for c in COHORTS for k in ('allele', 'residue')], ignore_index=True)
        self.manifest = json.loads((R / '_run_info/run_manifest.json').read_text())
        if self.manifest.get('component') != 'assoc_hla.v2':
            die('the association run manifest is not that of assoc_hla.v2')

    def t(self, cs, c, k, f):
        return float(self.thr.loc[(self.thr.case_set == cs) & (self.thr.cohort == c) & (self.thr['class'] == k), f].iloc[0])

    def tests(self, model='random'):
        """One row per test (a mirror pair once)."""
        A = self.A
        return A[(A.model == model) & (A.test_id == A.marker)]

    def pheno(self, cs, c, s='random'):
        return pd.read_csv(self.R / cs / c / '00.prep' / f'pheno_covar.{s}.tsv', sep='\t', dtype={'IID': str})


# ── A1 ───────────────────────────────────────────────────────────────────────
def a1(X):
    R, a = X.R, X.a
    case_sets = pd.read_csv(R / '_reference/case_sets.tsv', sep='\t', dtype=str)
    rows = []
    for cs in COMPARISONS:
        for c in COHORTS:
            pc = json.loads((R / cs / c / '00.prep/pheno_coding.json').read_text())
            nm = json.loads((R / cs / c / '00.prep/null_model.json').read_text())
            ra, fx = X.pheno(cs, c, 'random'), X.pheno(cs, c, 'fixed')
            row = {'comparison': cs, 'cohort': c, 'case_side': COMP_SIDES[cs][0], 'control_side': COMP_SIDES[cs][1],
                   'analysis': 'Main (Q1–Q3)' if cs == 'all' else 'Supplementary (Q4)', 'k': K_COMPONENTS[c],
                   'random_samples': len(ra), 'random_cases': int((ra.PHENO1 == 1).sum()),
                   'random_controls': int((ra.PHENO1 == 0).sum()), 'fixed_samples': len(fx),
                   'fixed_cases': int((fx.PHENO1 == 1).sum()), 'fixed_controls': int((fx.PHENO1 == 0).sum()),
                   'relatives_removed': len(ra) - len(fx),
                   'null_model': 'genome-wide scan, reused' if cs == 'all' else 'fitted for this comparison',
                   'variance_ratio': float(nm['variance_ratio'][0])}
            for k in ('fixed', 'random'):
                if pc[k]['samples'] != row[f'{k}_samples'] or pc[k]['cases'] != row[f'{k}_cases']:
                    die(f'A1 {cs} {c}: sample counts differ from the pipeline record')
            if cs == 'all':
                ph = pd.read_csv(Path(a.model_inputs) / c / 'phenotype/pheno.tsv', sep='\t', dtype=str)
                cases = case_sets[case_sets.IID.isin(ph.loc[ph.PHENO1 == '1', 'IID'])]
                for k in ('yes', 'no', 'unknown'):
                    row[f'dvt_{k}'] = int((cases.dvt == k).sum())
                if row['dvt_yes'] + row['dvt_no'] + row['dvt_unknown'] != row['random_cases']:
                    die(f'A1 {c}: DVT strata do not add up to the cases')
            rows.append(row)
    samples = pd.DataFrame(rows)
    # GRM markers of the random model (the genome-wide scan's LD-pruned set)
    grm = {c: sum(1 for _ in open(Path(a.saige) / f'results/{c}/00.prep/grm.prune.in')) for c in COHORTS}
    # tests per gene and cohort, all cases: alleles and residues (one per mirror pair), positions tested
    T = X.tests()
    T = T[T.case_set == 'all']
    O = X.O[(X.O.case_set == 'all') & (X.O.errcode != 'NOT_TESTED')]
    genes = list(X.TC.ASSOC_GENES)
    tg = []
    for g in genes:
        row = {'gene': g}
        for c in COHORTS:
            d = T[(T.cohort == c) & (T.gene == g)]
            row[f'allele_{c}'] = int((d['class'] == 'allele').sum())
            row[f'residue_{c}'] = int((d['class'] == 'residue').sum())
            row[f'position_{c}'] = int(((O.cohort == c) & (O.gene == g)).sum())
        tg.append(row)
    tg = pd.DataFrame(tg)
    tot = {'gene': 'All'}
    for col in tg.columns[1:]:
        tot[col] = int(tg[col].sum())
    tg = pd.concat([tg, pd.DataFrame([tot])], ignore_index=True)
    for c in COHORTS:
        for k in ('allele', 'residue'):
            if tot[f'{k}_{c}'] != int(X.t('all', c, k, 'tests')):
                die(f'A1 {c} {k}: tests per gene do not add up to the scan')
    # the methods record (slide 7), each statement checked against the run manifest
    M = X.manifest
    checks = [('random_model', ['--LOCO=TRUE', '--impute_method=mean']), ('random_model_p', ['p.value.NA']),
              ('fixed_model', ['no-firth']), ('covariates', ['PC1_AVG-PC10_AVG', 'SEX']),
              ('omnibus', ['Bonferroni'.lower()]), ('tests', ['mirror pair'])]
    for key, toks in checks:
        txt = json.dumps(M[key]).lower()
        for t_ in toks:
            if t_.lower() not in txt:
                die(f'A1: the run manifest no longer says "{t_}" under {key}; revise the methods record')
    methods = pd.DataFrame([
        ('Both models', 'Covariates', 'sex and PC1–PC10 of the projection onto BBJ mainland'),
        ('Both models', 'Outcome', 'case side (1) vs. control side (0) of each comparison'),
        ('Single markers (random model)', 'Model', 'logistic mixed model, SAIGE 1.1.9; full GRM from genome-wide '
         'markers, chromosome 6 left out (LOCO); relatives kept'),
        ('Single markers (random model)', 'GRM markers (N / I / F)', ' / '.join(f'{grm[c]:,}' for c in COHORTS)),
        ('Single markers (random model)', 'Test and P', 'score test of one allele or residue column; P without SPA '
         '(the saddlepoint P is kept as a sensitivity analysis)'),
        ('Single markers (random model)', 'Effect', 'OR = exp(β) per copy (SAIGE\'s one-step estimate)'),
        ('Single markers (random model)', 'Undetermined copies', 'replaced by the mean copy number'),
        ('Single markers (random model)', 'Null model', 'all cases vs. controls: the genome-wide scan\'s; '
         'DVT comparisons: fitted here with the same GRM'),
        ('Whole positions (fixed model)', 'Model', 'fixed-effect logistic regression on unrelated samples; one copy '
         'column per residue (residues below 1% pooled), the most frequent residue as reference'),
        ('Whole positions (fixed model)', 'Test', 'likelihood-ratio test of all residue terms; df = residues − 1'),
        ('Single markers, sensitivity', 'Fixed model', 'the same allele and residue tests as fixed-effect logistic '
         'regression (no Firth correction), complete cases per marker'),
        ('Thresholds', 'Single markers', 'LD-adjusted (GEC v0.2): 0.05 / Me per scan (primary); conventional '
         'Bonferroni 0.05 / tests (sensitivity)'),
        ('Thresholds', 'Whole positions', 'Bonferroni: 0.05 / positions tested per scan'),
        ('Thresholds', 'Benchmark', 'exact max-T by parametric bootstrap, 10,000 null phenotypes; no call depends on it'),
        ('Thresholds', 'Counting', 'a mirror pair (two residue columns that are complements in every sample) counts '
         'as one test'),
        ('Consistency', 'Across cohorts', 'significant in n/3 nested cohorts, same direction; robustness to the '
         'ancestry definition, not replication')], columns=['scope', 'item', 'value'])
    return samples, tg, methods, grm


# ── A2 ───────────────────────────────────────────────────────────────────────
def a2(X):
    thr, boot = X.thr.copy(), X.boot.copy()
    keep = ['case_set', 'cohort', 'class', 'markers', 'polymorphic', 'mirror_pairs', 'tests', 'gec_samples', 'me_gec',
            'p_gec', 'p_bonf', 'gec_blocks']
    sc = thr[keep].rename(columns={'case_set': 'comparison'})
    sc['me_share'] = sc.me_gec / sc.tests
    sc['bonf_stricter'] = sc.p_gec / sc.p_bonf
    if not np.allclose(sc.p_gec, 0.05 / sc.me_gec, rtol=1e-6) or not np.allclose(sc.p_bonf, 0.05 / sc.tests, rtol=1e-6):
        die('A2: thresholds are not 0.05/Me and 0.05/tests')
    bcols = ['samples', 'cases', 'controls', 'draws', 'cases_mean', 'cases_sd']
    for s in ('score', 'wald'):
        bcols += [f'p_maxt_{s}', f'p_maxt_{s}_lo', f'p_maxt_{s}_hi', f'me_maxt_{s}', f'fwer_gec_{s}', f'fwer_gec_{s}_lo',
                  f'fwer_gec_{s}_hi', f'fwer_bonf_{s}', f'fwer_bonf_{s}_lo', f'fwer_bonf_{s}_hi', f'replicates_with_na_{s}']
    b = boot[['case_set', 'cohort', 'class'] + bcols].rename(columns={'case_set': 'comparison'})
    b = b.rename(columns={c: f'boot_{c}' for c in bcols})
    sc = sc.merge(b, on=['comparison', 'cohort', 'class'], how='left', validate='1:1')
    if sc.boot_draws.isna().any():
        die('A2: a scan without its bootstrap benchmark')
    # draws beyond each threshold (FWER x draws), as slide 15 counts them
    for s in ('score', 'wald'):
        sc[f'boot_beyond_gec_{s}'] = (sc[f'boot_fwer_gec_{s}'] * sc.boot_draws).round().astype(int)
        sc[f'boot_beyond_bonf_{s}'] = (sc[f'boot_fwer_bonf_{s}'] * sc.boot_draws).round().astype(int)
    # positions: Bonferroni over the positions tested (separated fits counted, never significant)
    O = X.O
    pos = []
    for cs in COMPARISONS:
        for c in COHORTS:
            d = O[(O.case_set == cs) & (O.cohort == c)]
            tested = d[d.errcode != 'NOT_TESTED']
            thr_ = float(tested.threshold.iloc[0])
            if not np.isclose(thr_, 0.05 / len(tested), rtol=1e-6):
                die(f'A2 {cs} {c}: the position threshold is not 0.05 / positions tested')
            pos.append({'comparison': cs, 'cohort': c, 'positions': len(d), 'tested': len(tested),
                        'not_tested': int((d.errcode == 'NOT_TESTED').sum()),
                        'separated': int((d.errcode == 'SEPARATION').sum()), 'threshold': thr_})
    pos = pd.DataFrame(pos)
    return sc, pos


def gec_check(X):
    """GEC's Me recomputed from the genotypes of the all-cases full-mainland scans (the deck's slide-11 check)."""
    R = X.R
    poly = lambda r: (0.2982 * r ** 6 - 0.0127 * r ** 5 + 0.0588 * r ** 4 + 0.0099 * r ** 3   # noqa: E731
                      + 0.6281 * r ** 2 - 0.0009 * r)
    out = []
    for cls in ('allele', 'residue'):
        d = R / 'all' / 'full_mainland' / '01.markers' / cls
        ts = pd.read_csv(d / f'tests.{cls}.tsv', sep='\t', dtype=str, keep_default_na=False)
        with tempfile.TemporaryDirectory(dir='.') as tmp:
            ext = Path(tmp) / 'gec.extract'
            ext.write_text('\n'.join(ts.loc[ts.in_gec == '1', 'marker']) + '\n')
            subprocess.run([X.a.plink2, '--bfile', str(d / cls), '--keep', str(R / 'all/full_mainland/00.prep/fixed.keep'),
                            '--extract', str(ext), '--mac', '1', '--export', 'A', '--out', str(Path(tmp) / 'g')],
                           check=True, stdout=subprocess.DEVNULL)
            raw = pd.read_csv(Path(tmp) / 'g.raw', sep=r'\s+', na_values=['NA'])
        x = raw.iloc[:, 6:]
        x.columns = [c.rsplit('_', 1)[0] for c in x.columns]
        bim = pd.read_csv(d / f'{cls}.bim', sep='\t', header=None, names=['c', 'id', 'cm', 'pos', 'a1', 'a2'])
        pos = bim.set_index('id').pos
        blk = pd.read_csv(d / f'gec.{cls}.block.txt', sep='\t')
        mine, dup = [], 0
        for b in blk.itertuples(index=False):
            ids = [i for i in x.columns if b.StartPos <= pos[i] <= b.EndPos]
            if len(ids) == 1:
                mine.append(1.0)
                continue
            r = np.nan_to_num(x[ids].corr().to_numpy())
            np.fill_diagonal(r, 1.0)
            keep = []
            for j in range(len(ids)):
                if all(abs(r[j, k]) < 0.9999 for k in keep):
                    keep.append(j)
            dup += len(ids) - len(keep)
            rho = poly(r[np.ix_(keep, keep)])
            np.fill_diagonal(rho, 1.0)
            lam = np.clip(np.linalg.eigvalsh(rho), 0, None)
            mine.append(float(len(keep) - np.sum(np.where(lam > 1, lam - 1, 0))))
        gec = blk['#Effective'].to_numpy()
        out.append({'class': cls, 'tests': int(sum(blk['#Observed'])), 'blocks': int(len(blk)), 'duplicates': int(dup),
                    'me_recomputed': float(sum(mine)), 'me_gec': float(gec.sum()),
                    'max_block_diff': float(np.max(np.abs(np.array(mine) - gec)))})
    return pd.DataFrame(out)


# ── A3 ───────────────────────────────────────────────────────────────────────
def a3(X):
    T = X.tests()
    rows = []
    for cs in COMPARISONS:
        for c in COHORTS:
            for k in ('allele', 'residue'):
                d = T[(T.case_set == cs) & (T.cohort == c) & (T['class'] == k)]
                g, b = X.t(cs, c, k, 'p_gec'), X.t(cs, c, k, 'p_bonf')
                top = d.loc[d.p.idxmin()]
                if int((d.p < g).sum()) != int(d.sig_gec.sum()) or int((d.p < b).sum()) != int(d.sig_bonf.sum()):
                    die(f'A3 {cs} {c} {k}: significance flags differ from P < threshold')
                rows.append({'comparison': cs, 'cohort': c, 'test_type': k, 'tests': len(d), 'threshold_primary': g,
                             'threshold_bonf': b, 'sig_primary': int((d.p < g).sum()), 'sig_bonf': int((d.p < b).sum()),
                             'sig_spa': int((d.p_spa < g).sum()), 'min_p': float(top.p), 'top': top.marker,
                             'top_or': float(top['or']), 'lambda_gc': lam_gc(d.p)})
            d = X.O[(X.O.case_set == cs) & (X.O.cohort == c) & (X.O.errcode != 'NOT_TESTED')]
            ok = d[d.errcode == '.']
            top = ok.loc[ok.p.idxmin()]
            rows.append({'comparison': cs, 'cohort': c, 'test_type': 'position', 'tests': len(d),
                         'threshold_primary': float(d.threshold.iloc[0]), 'threshold_bonf': np.nan,
                         'sig_primary': int(d.sig_bonf.sum()), 'sig_bonf': np.nan, 'sig_spa': np.nan,
                         'min_p': float(top.p), 'top': top.unit, 'top_or': np.nan, 'lambda_gc': lam_gc(ok.p)})
    scans = pd.DataFrame(rows)
    # the five strongest of each test type, all cases (slides 18-20): alleles; residues one per position (its
    # strongest residue); positions; ranked by the smallest P in any cohort
    A = T[T.case_set == 'all']
    top = []
    for k in ('allele', 'residue'):
        D = A[A['class'] == k]
        w = D.pivot_table(index='marker', columns='cohort', values='p', aggfunc='first')
        order = w.min(axis=1).sort_values().index
        if k == 'residue':
            order = pd.Series(order, index=order.str.rsplit(':', n=1).str[0]).groupby(level=0, sort=False).first()
            order = order.loc[order.index.unique()].tolist()
        for rank, m in enumerate(list(order)[:5], 1):
            d = D[D.marker == m].set_index('cohort')
            f = d.loc['full_mainland']
            row = {'test_type': k, 'rank': rank, 'marker': m, 'gene': f.gene, 'freq_case': float(f.freq_case),
                   'freq_control': float(f.freq_control),
                   'n3_primary': int(sum(d.loc[c].p < X.t('all', c, k, 'p_gec') for c in COHORTS)),
                   'n3_bonf': int(sum(d.loc[c].p < X.t('all', c, k, 'p_bonf') for c in COHORTS))}
            for c in COHORTS:
                row[f'p_{c}'], row[f'or_{c}'] = float(d.loc[c].p), float(d.loc[c]['or'])
            top.append(row)
    O = X.O[(X.O.case_set == 'all') & (X.O.errcode == '.')]
    w = O.pivot_table(index='unit', columns='cohort', values='p', aggfunc='first')
    for rank, u in enumerate(w.min(axis=1).sort_values().index[:5], 1):
        d = O[O.unit == u].set_index('cohort')
        f = d.loc['full_mainland']
        row = {'test_type': 'position', 'rank': rank, 'marker': u, 'gene': f.gene,
               'residues': ', '.join(sorted([f.reference] + str(f.categories).split(','))), 'df': int(f.df),
               'n3_primary': int(d.sig_bonf.astype(bool).sum())}
        for c in COHORTS:
            row[f'p_{c}'] = float(d.loc[c].p)
        top.append(row)
    return scans, pd.DataFrame(top)


# ── A4 ───────────────────────────────────────────────────────────────────────
def a4(X):
    R = X.R
    st = pd.read_csv(R / '_conditional/steps_all.tsv', sep='\t', dtype={'lead': str})
    st = st[st.model == 'random']
    sg = pd.read_csv(R / '_conditional/signals_all.tsv', sep='\t', dtype={'lead': str})
    sg = sg[sg.model == 'random'].set_index(['case_set', 'class', 'lead'])
    rows = []
    for r in st.itertuples():
        row = {'comparison': r.case_set, 'test_type': r._3, 'step': int(r.step),
               'conditioned_on': '' if r.conditioned_on == '-' else r.conditioned_on, 'lead': r.lead,
               'entered': bool(r.entered), 'reason': r.reason}
        for c, s in zip(COHORTS, ('narrow', 'intermediate', 'full')):
            row[f'p_{c}'] = getattr(r, f'p_{s}')
            row[f'or_{c}'] = getattr(r, f'or_{s}')
            row[f'sig_{c}'] = getattr(r, f'n_sig_gec_{s}')
        row['n3'] = r.gec
        key = (r.case_set, r._3, r.lead)
        if r.entered and key in sg.index:
            s = sg.loc[key]
            row['members'] = int(s.n_members)
            row['carrier_alleles'] = '' if pd.isna(s.carrier_alleles) else str(s.carrier_alleles)
        rows.append(row)
    ps = pd.read_csv(R / '_omnibus/omnibus_steps_all.tsv', sep='\t', dtype={'lead': str})
    psg = pd.read_csv(R / '_omnibus/omnibus_signals_all.tsv', sep='\t', dtype={'lead': str}).set_index(['case_set', 'lead'])
    for r in ps.itertuples():
        row = {'comparison': r.case_set, 'test_type': 'position', 'step': int(r.step),
               'conditioned_on': '' if r.conditioned_on == '-' else r.conditioned_on,
               'lead': '' if pd.isna(r.lead) else r.lead, 'entered': bool(r.entered), 'reason': r.reason}
        for c in COHORTS:
            row[f'p_{c}'] = getattr(r, f'p_{c}')
            row[f'sig_{c}'] = getattr(r, f'n_sig_{c}')
        row['n3'] = r.consistency
        if r.entered and (r.case_set, r.lead) in psg.index:
            row['members'] = int(psg.loc[(r.case_set, r.lead)].n_members)
        if not row['lead']:
            # the omnibus walk records no unit when none passes; as the residue walks do, report the step's
            # strongest position (smallest P in any cohort) and its P in each cohort, from the step's own scan
            if r.step == 0:
                d = X.O[(X.O.case_set == r.case_set) & (X.O.errcode == '.')]
            else:
                d = pd.read_csv(R / r.case_set / '_omnibus/omnibus_stepwise.residue.tsv', sep='\t', dtype={'unit': str})
                d = d[(d.step == r.step) & (d.errcode == '.')]
            w = d.pivot_table(index='unit', columns='cohort', values='p', aggfunc='first')
            u = w.min(axis=1).idxmin()
            ds = d[d.unit == u].set_index('cohort')
            if ds.sig_bonf.astype(bool).any():
                die(f'A4 {r.case_set} position step {r.step}: {u} passes but the walk recorded no lead')
            row['lead'] = u
            for c in COHORTS:
                row[f'p_{c}'], row[f'sig_{c}'] = float(ds.loc[c].p), 0
            row['n3'] = '0/3'
        rows.append(row)
    steps = pd.DataFrame(rows)
    # what is left after the first step, per comparison x test type x cohort (from the walks' own scans)
    after = []
    for cs in COMPARISONS:
        for k in ('allele', 'residue'):
            f = R / cs / '_conditional/random' / f'stepwise.{k}.tsv'
            s = pd.read_csv(f, sep='\t', dtype={'marker': str})
            for c in COHORTS:
                d = s[(s.step == 1) & (s.cohort == c) & (s.errcode == '.')]
                if len(d):
                    top = d.loc[d.p.idxmin()]
                    after.append({'comparison': cs, 'test_type': k, 'cohort': c, 'min_p': float(top.p),
                                  'top': top.marker, 'sig': int(d.sig_gec.sum())})
        s = pd.read_csv(R / cs / '_omnibus/omnibus_stepwise.residue.tsv', sep='\t', dtype={'unit': str})
        for c in COHORTS:
            d = s[(s.step == 1) & (s.cohort == c) & (s.errcode == '.')]
            if len(d):
                top = d.loc[d.p.idxmin()]
                after.append({'comparison': cs, 'test_type': 'position', 'cohort': c, 'min_p': float(top.p),
                              'top': top.unit, 'sig': int(d.sig_bonf.astype(bool).sum())})
    # no allele enters the all-cases walk, so the allele scan given DQB1 57 V comes from HLA_GIVEN_SIGNAL (the deck's
    # summary slide reads it there; one test per mirror pair as everywhere)
    T = X.tests()
    for c in COHORTS:
        reps = set(T[(T.case_set == 'all') & (T.cohort == c) & (T['class'] == 'allele')].marker)
        g = pd.read_csv(R / f'all/{c}/04.conditional/random/hla_given_signal.tsv', sep='\t', dtype={'marker': str})
        g = g[(g.conditioned_on == LEAD57) & (g['class'] == 'allele') & (g.errcode == '.') & g.marker.isin(reps)]
        if len(g) != len(reps):
            die(f'A4 alleles given {LEAD57}, {c}: {len(g)} of {len(reps)} tests')
        top = g.loc[g.p_given.idxmin()]
        after.append({'comparison': 'all', 'test_type': 'allele', 'cohort': c, 'min_p': float(top.p_given),
                      'top': top.marker, 'sig': int((g.p_given < X.t('all', c, 'allele', 'p_gec')).sum())})
    # B 114 before and after position 57 (all cases)
    s = pd.read_csv(R / 'all/_omnibus/omnibus_stepwise.residue.tsv', sep='\t', dtype={'unit': str})
    b114 = {c: {'before': float(s[(s.step == 0) & (s.cohort == c) & (s.unit == 'B:114')].p.iloc[0]),
                'after': float(s[(s.step == 1) & (s.cohort == c) & (s.unit == 'B:114')].p.iloc[0])} for c in COHORTS}
    return steps, pd.DataFrame(after), b114


# ── A5 ───────────────────────────────────────────────────────────────────────
def a5(X):
    R = X.R
    gl = pd.read_csv(R / '_conditional/gwas_leads.tsv', sep='\t')
    leads = []
    for r in gl.itertuples():
        row = {'lead': r.lead_id, 'position': int(r.pos), 'status': r.status, 'reason': r.reason}
        for c in COHORTS:
            row[f'p_{c}'] = getattr(r, f'p_{c}')
            row[f'or_{c}'] = getattr(r, f'or_{c}')
            row[f'freq_{c}'] = getattr(r, f'a1_freq_{c}')
        leads.append(row)
    leads = pd.DataFrame(leads)
    kept = gl[gl.status == 'kept'].lead_id.tolist() if 'kept' in set(gl.status) else gl.lead_id[:1].tolist()
    lg = pd.read_csv(R / '_conditional/lead_given_hla_all.tsv', sep='\t',
                     dtype={'conditioned_on': str, 'position': str, 'position_reference': str})
    lg = lg[lg.model == 'random'].copy()
    lg['r2'] = lg.r ** 2
    lg = lg.rename(columns={'case_set': 'comparison'})
    cols = ['comparison', 'cohort', 'lead', 'conditioned_on', 'why', 'position', 'r', 'r2', 'lead_p_marginal',
            'lead_or_marginal', 'lead_p_given', 'lead_or_given', 'marker_p_marginal', 'marker_p_given_lead']
    lg = lg[cols]
    # what passes given the lead, per comparison x cohort x test type (single markers: random model, GEC; positions:
    # the omnibus with the lead as a covariate, Bonferroni)
    hg = pd.read_csv(R / '_conditional/hla_given_lead_all.tsv', sep='\t', dtype={'marker': str})
    hg = hg[hg.model == 'random']
    reps = X.tests()
    reps = {(r.case_set, r.cohort, r.marker) for r in reps.itertuples()}
    given = []
    for cs in COMPARISONS:
        og = pd.read_csv(R / cs / '_omnibus/omnibus_given.residue.tsv', sep='\t', dtype={'unit': str})
        og = og[og.why == 'gwas_lead']
        for c in COHORTS:
            for k in ('allele', 'residue'):
                d = hg[(hg.case_set == cs) & (hg.cohort == c) & (hg['class'] == k) & (hg.errcode == '.')]
                d = d[[(cs, c, m) in reps for m in d.marker]]
                g = X.t(cs, c, k, 'p_gec')
                top = d.loc[d.p_given_lead.idxmin()]
                given.append({'comparison': cs, 'cohort': c, 'test_type': k, 'tests': len(d), 'threshold': g,
                              'sig_before': int((d.p_marginal < g).sum()), 'sig_given': int((d.p_given_lead < g).sum()),
                              'passing_given': '; '.join(d[d.p_given_lead < g].marker), 'min_p_given': float(top.p_given_lead),
                              'top_given': top.marker, 'lambda_before': lam_gc(d.p_marginal),
                              'lambda_given': lam_gc(d.p_given_lead)})
            d = og[(og.cohort == c) & (og.errcode == '.')]
            top = d.loc[d.p.idxmin()]
            given.append({'comparison': cs, 'cohort': c, 'test_type': 'position', 'tests': len(d),
                          'threshold': float(d.threshold.iloc[0]), 'sig_before': int((d.p_marginal < d.threshold).sum()),
                          'sig_given': int(d.sig_bonf.sum()), 'passing_given': '; '.join(d[d.sig_bonf].unit),
                          'min_p_given': float(top.p), 'top_given': top.unit, 'lambda_before': lam_gc(d.p_marginal),
                          'lambda_given': lam_gc(d.p)})
    given = pd.DataFrame(given)
    # the link to DQB1 57 V, all cases (summary.tsv): the lead alone, given 57 V, given position 57, given its best
    # LD tag; 57 V alone and given the lead; position 57 alone and given the lead
    S0 = pd.read_csv(R / '_conditional/summary.tsv', sep='\t')
    S0 = S0[(S0.case_set == 'all') & (S0.model == 'random')].set_index('cohort')
    L = [c.split('|')[0] for c in S0.columns if c.endswith('|marginal_p')]
    if len(L) != 1 or (S0.residue_signal1 != LEAD57).any():
        die(f'A5: the GWAS link changed: leads {L}')
    L = L[0]
    og = pd.read_csv(R / 'all/_omnibus/omnibus_given.residue.tsv', sep='\t', dtype={'unit': str})
    og = og[(og.unit == 'DQB1:57') & (og.conditioned_on == L)].set_index('cohort')
    om = X.O[(X.O.case_set == 'all') & (X.O.unit == 'DQB1:57')].set_index('cohort')
    ld = pd.read_csv(R / '_conditional/ld_lead_hla.tsv', sep='\t')
    link = []
    for c in COHORTS:
        g = lambda k: S0.loc[c, f'{L}|{k}']  # noqa: E731
        v = X.A[(X.A.case_set == 'all') & (X.A.model == 'random') & (X.A.cohort == c) & (X.A.marker == LEAD57)].iloc[0]
        link.append({'cohort': c, 'lead': L, 'lead_p': float(g('marginal_p')), 'lead_given_57v': float(g('p_given_residue_signal1')),
                     'lead_given_pos57': float(g('p_given_omnibus_position1')), 'tag': g('best_ld_tag'),
                     'tag_r2': float(g('best_ld_tag_r2')), 'lead_given_tag': float(g('p_given_best_ld_tag')),
                     'v57_p': float(v.p), 'v57_given_lead': float(g('residue_signal1_p_given_lead')),
                     'pos57_p': float(om.loc[c, 'p']), 'pos57_given_lead': float(og.loc[c, 'p']),
                     'pos57_threshold': float(om.loc[c, 'threshold']), 'hla_min_given_lead': float(g('hla_min_p_given_lead')),
                     'hla_top_given_lead': g('hla_top_given_lead'), 'most_explaining': g('most_explaining'),
                     'r_57v': float(ld[(ld.cohort == c) & (ld.marker == LEAD57)].r.iloc[0])})
    return leads, lg, given, pd.DataFrame(link), kept


# ── A6 ───────────────────────────────────────────────────────────────────────
def a6(X, gw):
    R = X.R
    T = X.tests()
    A = T[T.case_set == 'all']
    reps = {k: set(d.marker) for k, d in A.groupby(['cohort', 'class'])}
    lead = pd.read_csv(R / '_conditional/gwas_leads.tsv', sep='\t')
    H = pd.read_csv(R / '_conditional/hla_given_lead_all.tsv', sep='\t')
    H = H[(H.case_set == 'all') & (H.model == 'random') & (H.errcode == '.')]
    Sg = pd.concat([pd.read_csv(R / f'all/{c}/04.conditional/random/hla_given_signal.tsv', sep='\t') for c in COHORTS])
    Sg = Sg[Sg.errcode == '.']
    Og = pd.read_csv(R / 'all/_omnibus/omnibus_given.residue.tsv', sep='\t', dtype={'unit': str})
    Og = Og[Og.errcode == '.']
    P = pd.read_csv(R / 'all/_omnibus/omnibus_stepwise.residue.tsv', sep='\t', dtype={'unit': str})
    P = P[P.errcode == '.']
    L = [c for c in set(Og.conditioned_on) if c.startswith('chr')]
    if len(L) != 1:
        die(f'A6: one GWAS lead expected among the omnibus conditioning sets, found {L}')
    L = L[0]
    rows = []
    for c in COHORTS:
        for cls in ('allele', 'residue'):
            d = A[(A.cohort == c) & (A['class'] == cls)]
            rows.append({'test_type': cls, 'cohort': c, 'conditioning': 'none', 'tests': len(d), 'before': lam_gc(d.p),
                         'after': np.nan})
            d = H[(H.cohort == c) & (H['class'] == cls) & H.marker.isin(reps[(c, cls)])]
            rows.append({'test_type': cls, 'cohort': c, 'conditioning': 'GWAS lead', 'tests': len(d),
                         'before': lam_gc(d.p_marginal), 'after': lam_gc(d.p_given_lead)})
            for k, cond in (('DQB1 57 V', LEAD57), ('position 57', POS57_SET)):
                d = Sg[(Sg.cohort == c) & (Sg['class'] == cls) & (Sg.conditioned_on == cond) & Sg.marker.isin(reps[(c, cls)])]
                rows.append({'test_type': cls, 'cohort': c, 'conditioning': k, 'tests': len(d),
                             'before': lam_gc(d.p_marginal), 'after': lam_gc(d.p_given)})
        o = X.O[(X.O.case_set == 'all') & (X.O.cohort == c) & (X.O.errcode == '.')]
        rows.append({'test_type': 'position', 'cohort': c, 'conditioning': 'none', 'tests': len(o), 'before': lam_gc(o.p),
                     'after': np.nan})
        for k, cond in (('GWAS lead', L), ('DQB1 57 V', LEAD57)):
            d = Og[(Og.cohort == c) & (Og.conditioned_on == cond)]
            rows.append({'test_type': 'position', 'cohort': c, 'conditioning': k, 'tests': len(d),
                         'before': lam_gc(d.p_marginal), 'after': lam_gc(d.p)})
        d0, d1 = P[(P.cohort == c) & (P.step == 0)], P[(P.cohort == c) & (P.step == 1)]
        d0 = d0[d0.unit.isin(set(d1.unit))]
        rows.append({'test_type': 'position', 'cohort': c, 'conditioning': 'position 57', 'tests': len(d1),
                     'before': lam_gc(d0.p), 'after': lam_gc(d1.p)})
        # positions by df, before and given 57 V
        G = Og[(Og.cohort == c) & (Og.conditioned_on == LEAD57)]
        for k, m in (('position, 1 df', o.df == 1), ('position, ≥ 2 df', o.df >= 2)):
            u = set(o[m].unit)
            rows.append({'test_type': k, 'cohort': c, 'conditioning': 'DQB1 57 V', 'tests': int(m.sum()),
                         'before': lam_gc(o[m].p), 'after': lam_gc(G[G.unit.isin(u)].p)})
        for k in ('mhc', 'rest'):
            rows.append({'test_type': 'GWAS SNPs ' + ('in the MHC' if k == 'mhc' else 'outside the MHC'), 'cohort': c,
                         'conditioning': 'none', 'tests': gw[c][f'{k}_snps'], 'before': gw[c][k], 'after': np.nan})
    lam = pd.DataFrame(rows)
    f = X.O[(X.O.case_set == 'all') & (X.O.cohort == 'full_mainland') & (X.O.errcode == '.') & (X.O.df >= 2)]
    multi = {'n': int(len(f)), 'in_dr_dq': int(f.unit.str.split(':').str[0].isin(['DRB1', 'DQB1', 'DQA1']).sum())}
    return lam, multi


# ── A7 ───────────────────────────────────────────────────────────────────────
def res_at(X, pos):
    """{2-field allele: residue at the position} through the reference sequence the typing matched it to."""
    ref = pd.read_csv(X.TY / 'results/05.exports/residue_reference/DQB1.tsv', sep='\t', dtype=str).set_index('allele')
    mm = pd.read_csv(X.TY / 'results/05.exports/allele_match_map.tsv', sep='\t', dtype=str)
    mm = dict(zip(mm.allele, mm.matched_to))

    def f(a):
        v = ref.at[mm[a], pos] if a in mm and mm[a] in ref.index else None
        return None if v in (None, '.', '*') else v
    return f


def dqb1_res57(X):
    """Residue at DQB1 57 of an allele or P-group name (any member's residue: a P group shares the peptide-binding
    domain), as the deck's summary slides read it."""
    ref = pd.read_csv(X.TY / 'results/05.exports/residue_reference/DQB1.tsv', sep='\t', dtype=str).set_index('allele')['57']

    def r(name):
        base = name[:-1] if name.endswith('P') else name
        vals = {ref[a] for a in ref.index if (a == base or a.startswith(base + ':')) and ref[a] not in ('.', '*')}
        return vals.pop() if len(vals) == 1 else None
    return r


def a7(X):
    R, TY = X.R, X.TY
    r57 = dqb1_res57(X)
    A = X.A[X.A.case_set == 'all']
    AR = A[A.model == 'random'].set_index(['cohort', 'marker'])
    AF = A[A.model == 'fixed'].set_index(['cohort', 'marker'])
    al_all = pd.read_csv(TY / 'results/05.exports/allele_dosage.tsv', sep='\t', dtype={'sample_id': str},
                         usecols=lambda k: k == 'sample_id' or k.startswith('DQB1*')).set_index('sample_id')
    O = X.O[(X.O.case_set == 'all') & (X.O.unit == 'DQB1:57')].set_index('cohort')
    geno, res_rows, al_rows = [], [], []
    RES = ['V', 'D', 'A', 'S']
    for c in COHORTS:
        ph = X.pheno('all', c).set_index('IID')
        rs = pd.read_csv(TY / f'results/06.assoc_inputs/{c}/residue_dosage.tsv', sep='\t', dtype={'sample_id': str},
                         usecols=['sample_id'] + [f'DQB1:57:{r}' for r in RES]).set_index('sample_id').loc[ph.index]
        al = al_all.loc[ph.index]
        rr, rf, o = AR.loc[(c, LEAD57)], AF.loc[(c, LEAD57)], O.loc[c]
        for grp, k in ((1, 'cases'), (0, 'controls')):
            v = rs.loc[ph.index[ph.PHENO1 == grp], LEAD57]
            called = v.dropna()
            row = {'cohort': c, 'group': k, 'analyzed': len(v), 'called': len(called), 'undetermined': int(v.isna().sum()),
                   'n0': int((called == 0).sum()), 'n1': int((called == 1).sum()), 'n2': int((called == 2).sum()),
                   'carriers': int((called >= 1).sum()), 'copy_freq': float(called.sum() / (2 * len(called)))}
            kk = 'case' if grp == 1 else 'control'
            if not np.isclose(row['copy_freq'], rr[f'freq_{kk}'], atol=1e-6):
                die(f'A7 {c} {k}: the matrix and the association disagree on the copy frequency of DQB1 57 V')
            geno.append(row)
        res_rows.append({'cohort': c, 'or': float(rr['or']), 'l95': float(rr.l95), 'u95': float(rr.u95), 'p': float(rr.p),
                         'gec': X.t('all', c, 'residue', 'p_gec'), 'bonf': X.t('all', c, 'residue', 'p_bonf'),
                         'or_fixed': float(rf['or']), 'l95_fixed': float(rf.l95), 'u95_fixed': float(rf.u95),
                         'p_fixed': float(rf.p), 'p_spa': float(rr.p_spa), 'chisq_pos57': float(o.chisq),
                         'df_pos57': int(o.df), 'p_pos57': float(o.p), 'threshold_pos57': float(o.threshold)})
        # the four residues at 57 and their alleles (slide 26): alleles >= 0.5% of control copies each, the rest pooled
        typed = al.loc[al.notna().any(axis=1)]
        ctrl = typed.loc[[i for i in typed.index if ph.loc[i, 'PHENO1'] == 0]]
        for r in RES:
            cols = [a_ for a_ in al.columns if r57(a_) == r]
            ar = AR.loc[(c, f'DQB1:57:{r}')]
            row = {'cohort': c, 'residue': r, 'level': 'residue', 'name': f'DQB1:57:{r}',
                   'n_alleles': int((al[cols].sum() > 0).sum()), 'or': float(ar['or']), 'l95': float(ar.l95),
                   'u95': float(ar.u95), 'p': float(ar.p)}
            for grp, k in ((1, 'case'), (0, 'control')):
                d = rs.loc[ph.index[ph.PHENO1 == grp], f'DQB1:57:{r}'].dropna()
                row[f'copies_{k}'], row[f'freq_{k}'] = int(d.sum()), float(d.sum() / (2 * len(d)))
            al_rows.append(row)
            named = [a_ for a_ in cols if ctrl[a_].sum() / (2 * len(ctrl)) >= 0.005]
            rare = [a_ for a_ in cols if a_ not in named and typed[a_].sum() > 0]
            for name, members in [(a_, [a_]) for a_ in sorted(named)] + ([(f'{len(rare)} rare alleles', rare)] if rare else []):
                e = {'cohort': c, 'residue': r, 'level': 'allele' if len(members) == 1 and name in named else 'rare',
                     'name': name, 'n_alleles': len(members)}
                for grp, k in ((1, 'case'), (0, 'control')):
                    idx = [i for i in typed.index if ph.loc[i, 'PHENO1'] == grp]
                    cop = float(typed.loc[idx, members].sum().sum())
                    e[f'copies_{k}'], e[f'freq_{k}'] = int(cop), cop / (2 * len(idx))
                if name in named and (c, name) in AR.index:
                    a_ = AR.loc[(c, name)]
                    e.update({'or': float(a_['or']), 'l95': float(a_.l95), 'u95': float(a_.u95), 'p': float(a_.p)})
                al_rows.append(e)
    # valine 57 against the two Japanese panels (slide 26 lower table), full mainland
    T = pd.read_csv(X.TC.TOMMO)
    T = T[T.gene == 'DQB1'].copy()
    r57s = lambda g: r57(g) or next((r57(g + x) for x in 'QLSCAN' if r57(g + x)), None)  # noqa: E731
    T['r57'] = T.allele.map(r57s)
    K = pd.read_csv(TY / 'results/03.panels/allele_concordance.tsv', sep='\t')
    K = K[(K.panel == '1KG JPT') & (K.gene == 'DQB1')].copy()
    n_kg = int(round((K.count_panel / K.freq_panel).replace([np.inf], np.nan).dropna().median()))
    K['r57'] = K.allele.map(r57s)
    g = [r for r in geno if r['cohort'] == 'full_mainland']
    panels = {'cases': next(r for r in g if r['group'] == 'cases')['copy_freq'],
              'controls': next(r for r in g if r['group'] == 'controls')['copy_freq'],
              'tommo': float(T[T.r57 == 'V'].frequency.sum()), 'tommo_chr': int(round(T['count'].sum())),
              'kg': float(K[K.r57 == 'V'].count_panel.sum() / n_kg), 'kg_chr': n_kg,
              'cases_chr': 2 * next(r for r in g if r['group'] == 'cases')['called'],
              'controls_chr': 2 * next(r for r in g if r['group'] == 'controls')['called']}
    sg = pd.read_csv(R / '_conditional/signals_all.tsv', sep='\t')
    sg = sg[(sg.case_set == 'all') & (sg.model == 'random') & (sg.lead == LEAD57)]
    carriers = [x.strip() for x in str(sg.carrier_alleles.iloc[0]).split(',')]
    return pd.DataFrame(geno), pd.DataFrame(res_rows), pd.DataFrame(al_rows), panels, carriers


# ── A8 ───────────────────────────────────────────────────────────────────────
def a8(X):
    R, a = X.R, X.a
    T = X.tests()
    strata = []
    for cs in COMPARISONS:
        for c in COHORTS:
            d = X.A[(X.A.case_set == cs) & (X.A.cohort == c) & (X.A.model == 'random')]
            own = d[d.test_id == d.marker]
            m = d[d.marker == LEAD57].iloc[0]
            strata.append({'comparison': cs, 'cohort': c, 'n_case': int(own.n_case.max()), 'n_control': int(own.n_control.max()),
                           'freq_case': float(m.freq_case), 'freq_control': float(m.freq_control), 'or': float(m['or']),
                           'l95': float(m.l95), 'u95': float(m.u95), 'p': float(m.p), 'gec': X.t(cs, c, 'residue', 'p_gec')})
    strata = pd.DataFrame(strata)
    # without DVT: what passes (slide 29): positions passing in all three cohorts; the one-df position -27, LRT vs score
    O = X.O[(X.O.case_set == 'dvt_no') & (X.O.errcode == '.')]
    sets = {c: set(O[(O.cohort == c) & O.sig_bonf].unit) for c in COHORTS}
    common = sorted(set.intersection(*sets.values()), key=lambda u: (u.split(':')[0], int(u.split(':')[1])))
    pos3 = pd.DataFrame([{'unit': u, **{f'p_{c}': float(O[(O.cohort == c) & (O.unit == u)].p.iloc[0]) for c in COHORTS}}
                         for u in common])
    A27 = X.A[(X.A.case_set == 'dvt_no') & (X.A.model == 'random') & X.A.marker.astype(str).str.startswith('DQB1:-27:')]
    lrt = []
    for c in COHORTS:
        o = O[(O.cohort == c) & (O.unit == 'DQB1:-27')].iloc[0]
        if int(o.df) != 1:
            die('A8: DQB1 -27 is no longer a one-df position')
        lrt.append({'cohort': c, 'p_lrt_fixed': float(o.p), 'p_score_random': float(A27[A27.cohort == c].p.iloc[0])})
    lrt = pd.DataFrame(lrt)
    lrt['ratio'] = lrt.p_score_random / lrt.p_lrt_fixed
    # the DQB1 alleles behind the without-DVT leads (slide 30): residues at 57, 87, 125 via the allele map
    ph = X.pheno('dvt_no', 'full_mainland').set_index('IID')
    al = pd.read_csv(X.TY / 'results/05.exports/allele_dosage.tsv', sep='\t', dtype={'sample_id': str},
                     usecols=lambda k: k == 'sample_id' or k.startswith('DQB1*')).set_index('sample_id')
    al = al.loc[ph.index].dropna(how='all')
    cp = {g: al.loc[al.index.intersection(ph.index[ph.PHENO1 == v])].sum() for g, v in (('case', 1), ('control', 0))}
    n2 = {g: int(2 * len(al.index.intersection(ph.index[ph.PHENO1 == v]))) for g, v in (('case', 1), ('control', 0))}
    seen = [x for x in al.columns if cp['case'][x] + cp['control'][x] > 0]
    f57, f87, f125 = res_at(X, '57'), res_at(X, '87'), res_at(X, '125')
    R3 = {x: {'57': f57(x), '87': f87(x), '125': f125(x)} for x in seen}
    AD = X.A[(X.A.case_set == 'dvt_no') & (X.A.model == 'random') & (X.A.cohort == 'full_mainland')].set_index('marker')

    def kind(x):
        if x in NODVT_NAMED:
            return x
        r = R3[x]
        if r['57'] is None or r['87'] is None:
            return 'unresolved'
        return {'Y': 'v_other' if r['57'] == 'V' else 'y_other', 'F': 'F', 'L': 'L'}[r['87']]
    rows = []
    for k in NODVT_NAMED + ['v_other', 'y_other', 'F', 'L', 'unresolved']:
        mem = [x for x in seen if kind(x) == k]
        if not mem:
            continue
        miss = [x for x in mem if R3[x]['125'] is None]
        e = {'kind': k, 'n_alleles': len(mem),
             **{f'res_{p}': ', '.join(sorted({R3[x][p] for x in mem if R3[x][p]})) for p in ('57', '87', '125')},
             'copies_case': int(sum(cp['case'][x] for x in mem)), 'copies_control': int(sum(cp['control'][x] for x in mem)),
             'no125_alleles': len(miss), 'no125_case': int(sum(cp['case'][x] for x in miss)),
             'no125_control': int(sum(cp['control'][x] for x in miss))}
        e['freq_case'], e['freq_control'] = e['copies_case'] / n2['case'], e['copies_control'] / n2['control']
        if k in NODVT_NAMED and k in AD.index:
            e['or'], e['p'] = float(AD.at[k, 'or']), float(AD.at[k, 'p'])
        rows.append(e)
    nodvt = pd.DataFrame(rows)
    tot = cp['case'] + cp['control']
    N = n2['case'] + n2['control']
    pv_ = sum(tot[x] for x in seen if R3[x]['57'] == 'V') / N
    py_ = sum(tot[x] for x in seen if R3[x]['87'] == 'Y') / N
    pvy = sum(tot[x] for x in seen if R3[x]['57'] == 'V' and R3[x]['87'] == 'Y') / N
    dd = pvy - pv_ * py_
    ld = {'r2': dd * dd / (pv_ * (1 - pv_) * py_ * (1 - py_)), 'dprime': dd / min(pv_ * (1 - py_), (1 - pv_) * py_),
          'r2_ceiling': pv_ * (1 - py_) / (py_ * (1 - pv_)), 'freq_v': pv_, 'freq_y': py_, 'n_cases': n2['case'] // 2,
          'n_controls': n2['control'] // 2,
          'v_not_y': sorted(x for x in seen if R3[x]['57'] == 'V' and R3[x]['87'] not in ('Y', None)),
          'v_not_y_control': int(sum(cp['control'][x] for x in seen if R3[x]['57'] == 'V' and R3[x]['87'] not in ('Y', None)))}
    # without DVT given DQB1 57 V: the supplementary run (outside the pipeline), every marginal P equal to the scan
    G = Path(a.given57v)
    Og = pd.read_csv(G / 'omnibus/omnibus_given.residue.tsv', sep='\t', dtype={'unit': str})
    Og = Og[(Og.conditioned_on == LEAD57) & (Og.errcode == '.')]
    Oa = X.O[(X.O.case_set == 'dvt_no') & (X.O.errcode == '.')]
    g57 = []
    checked = 0
    for c in COHORTS:
        g = pd.read_csv(G / c / 'hla_given_signal.tsv', sep='\t', dtype={'marker': str})
        g = g[(g.conditioned_on == LEAD57) & (g.why == 'signal_lead')].copy()
        ad_ = X.A[(X.A.case_set == 'dvt_no') & (X.A.model == 'random') & (X.A.cohort == c)].set_index('marker')
        chk = g[g.marker.isin(ad_.index)]
        if len(chk) != len(ad_) or not np.allclose(chk.p_marginal, ad_.loc[chk.marker, 'p'], rtol=1e-4):
            die(f'A8 {c}: the supplementary run does not reproduce the scan')
        checked += len(chk)
        g['test_id'] = g.marker.map(ad_.test_id)
        g = g[g.test_id == g.marker]
        o, m = Oa[Oa.cohort == c], Og[Og.cohort == c]
        if len(m) != len(o) or not np.allclose(m.set_index('unit').p_marginal.loc[o.unit].values, o.p.values, rtol=1e-4):
            die(f'A8 {c}: the supplementary omnibus run does not reproduce the scan')
        row = {'cohort': c}
        for k in ('allele', 'residue'):
            t_ = X.t('dvt_no', c, k, 'p_gec')
            d = g[(g['class'] == k) & (g.errcode == '.')]
            top = d.loc[d.p_given.idxmin()]
            row.update({f'{k}_before': int(((g['class'] == k) & (g.p_marginal < t_)).sum()),
                        f'{k}_after': int((d.p_given < t_).sum()), f'{k}_min_p': float(top.p_given), f'{k}_top': top.marker,
                        f'{k}_gec': t_})
        top = m.loc[m.p.idxmin()]
        row.update({'position_before': int(o.sig_bonf.sum()), 'position_after': int(m.sig_bonf.sum()),
                    'position_min_p': float(top.p), 'position_top': top.unit, 'position_threshold': float(m.threshold.iloc[0])})
        for k in (LEAD87, DPB, A0503):
            r = g[g.marker == k].iloc[0]
            row.update({f'{k}|p': float(r.p_marginal), f'{k}|p_given': float(r.p_given), f'{k}|or': float(r.or_marginal),
                        f'{k}|or_given': float(r.or_given)})
        ga = pd.read_csv(X.R / f'all/{c}/04.conditional/random/hla_given_signal.tsv', sep='\t', dtype={'marker': str})
        r = ga[(ga.conditioned_on == LEAD57) & (ga.why == 'signal_lead') & (ga.marker == A0503) & (ga.errcode == '.')].iloc[0]
        row.update({f'{A0503}|all|p': float(r.p_marginal), f'{A0503}|all|p_given': float(r.p_given),
                    f'{A0503}|all|or': float(r.or_marginal), f'{A0503}|all|or_given': float(r.or_given)})
        for u in (POS125, 'DQB1:87', 'DQB1:57'):
            r = m[m.unit == u].iloc[0]
            row.update({f'{u}|p': float(r.p_marginal), f'{u}|p_given': float(r.p), f'{u}|df': int(r.df)})
        g57.append(row)
    g57 = pd.DataFrame(g57)
    # the reverse direction within the comparison: 57 V given 87 Y (pipeline)
    rev = []
    for c in COHORTS:
        h = pd.read_csv(X.R / f'dvt_no/{c}/04.conditional/random/hla_given_signal.tsv', sep='\t', dtype={'marker': str})
        r = h[(h.marker == LEAD57) & (h.conditioned_on == LEAD87) & (h.errcode == '.')].iloc[0]
        rev.append({'cohort': c, 'p': float(r.p_marginal), 'p_given_87y': float(r.p_given), 'or_given_87y': float(r.or_given)})
    # SPA sensitivity of the calls (GEC), every comparison
    spa = []
    for cs in COMPARISONS:
        d = T[T.case_set == cs].copy()
        d['gec'] = [X.t(cs, c, k, 'p_gec') for c, k in zip(d.cohort, d['class'])]
        na, sp = d.p < d.gec, d.p_spa < d.gec
        for r in d[na & ~sp].itertuples():
            spa.append({'comparison': cs, 'cohort': r.cohort, 'marker': r.marker, 'p': float(r.p), 'p_spa': float(r.p_spa),
                        'threshold': float(r.gec)})
    spa_counts = {}
    for cs in COMPARISONS:
        d = T[T.case_set == cs]
        g_ = np.array([X.t(cs, c, k, 'p_gec') for c, k in zip(d.cohort, d['class'])])
        spa_counts[cs] = {'calls': int((d.p.values < g_).sum()), 'calls_spa': int((d.p_spa.values < g_).sum()),
                          'both': int(((d.p.values < g_) & (d.p_spa.values < g_)).sum())}
    return strata, pos3, lrt, nodvt, ld, g57, pd.DataFrame(rev), pd.DataFrame(spa), spa_counts, checked


# ── A9 ───────────────────────────────────────────────────────────────────────
def carrier_power(or_, n1, n0, p0, alpha=0.05):
    p1 = or_ * p0 / (1 - p0 + or_ * p0)
    se = np.sqrt(1 / (n1 * p1 * (1 - p1)) + 1 / (n0 * p0 * (1 - p0)))
    z, zc = abs(np.log(or_)) / se, stats.norm.ppf(1 - alpha / 2)
    return float(stats.norm.cdf(z - zc) + stats.norm.cdf(-z - zc))


def a9(X):
    TY, R = X.TY, X.R
    pan = pd.read_csv(TY / 'results/03.panels/allele_concordance.tsv', sep='\t')
    pop = {}
    for a_ in REP_ALLELES:
        r = pan[(pan.panel == 'ToMMo') & (pan.allele == a_)].iloc[0]
        k = pan[(pan.panel == '1KG JPT') & (pan.allele == a_)]
        pop[a_] = {'tommo': float(r.freq_panel), 'tommo_chr': int(round(r.count_panel / r.freq_panel)),
                   'controls_pgroup': float(r.freq_ours), 'kg': float(k.freq_panel.iloc[0]) if len(k) else np.nan,
                   'kg_chr': int(round(k.count_panel.iloc[0] / k.freq_panel.iloc[0])) if len(k) else np.nan,
                   'tommo_vs_controls_p': float(r.fisher_p),
                   'kg_vs_controls_p': float(k.fisher_p.iloc[0]) if len(k) else np.nan}
    pub = pd.DataFrame(PUBLISHED)
    for g in ('case', 'ctrl'):
        pub[f'{g}_freq_hwe'] = 1 - np.sqrt(1 - pub[f'{g}_car'] / pub[f'{g}_n'])
    # our per-copy results (random and fixed model) for all cases and without DVT, and copy frequencies
    rows = []
    for cs in ('all', 'dvt_no'):
        for c in COHORTS:
            for a_ in REP_ALLELES:
                r = X.A[(X.A.case_set == cs) & (X.A.cohort == c) & (X.A.marker == a_)].set_index('model')
                rows.append({'comparison': cs, 'cohort': c, 'allele': a_, 'freq_case': float(r.loc['random', 'freq_case']),
                             'freq_control': float(r.loc['random', 'freq_control']), 'or': float(r.loc['random', 'or']),
                             'l95': float(r.loc['random', 'l95']), 'u95': float(r.loc['random', 'u95']),
                             'p': float(r.loc['random', 'p']), 'or_fixed': float(r.loc['fixed', 'or']),
                             'p_fixed': float(r.loc['fixed', 'p'])})
    ours = pd.DataFrame(rows)
    # carriers on the papers' scale (2x2 on the fixed-model samples, unadjusted), with power at the published OR
    car = []
    for c in COHORTS:
        d = pd.read_csv(TY / f'results/06.assoc_inputs/{c}/allele_dosage.tsv', sep='\t', dtype=str).rename(columns={'sample_id': 'IID'})
        for cs, study in (('all', 'Tanabe 2005'), ('dvt_no', 'Kominami 2009')):
            ph = pd.read_csv(R / cs / c / '00.prep/pheno_covar.fixed.tsv', sep='\t', dtype=str)
            m = ph[['IID', 'PHENO1']].merge(d[['IID'] + REP_ALLELES], on='IID')
            y = m.PHENO1 == '1'
            for a_ in REP_ALLELES:
                g = pd.to_numeric(m[a_].replace('', np.nan), errors='coerce')
                ok, carry = g.notna(), g >= 1
                a1, n1 = int((carry & y & ok).sum()), int((y & ok).sum())
                c1, n0 = int((carry & ~y & ok).sum()), int((~y & ok).sum())
                b1, d1 = n1 - a1, n0 - c1
                o, se = a1 * d1 / (b1 * c1), np.sqrt(1 / a1 + 1 / b1 + 1 / c1 + 1 / d1)
                pb = next(x for x in PUBLISHED if x['study'] == study and x['allele'] == a_)
                p0 = c1 / n0
                pw = lambda o_: carrier_power(o_, n1, n0, p0)  # noqa: E731
                car.append({'cohort': c, 'comparison': cs, 'allele': a_, 'case_car': a1, 'case_n': n1, 'ctrl_car': c1,
                            'ctrl_n': n0, 'or': float(o), 'l95': float(np.exp(np.log(o) - 1.96 * se)),
                            'u95': float(np.exp(np.log(o) + 1.96 * se)),
                            'p': float(stats.fisher_exact([[a1, b1], [c1, d1]])[1]), 'compared_with': study,
                            'power_published': pw(pb['or_']),
                            'detectable_or80': float(optimize.brentq(lambda o_: pw(o_) - 0.8, 1.0001, 50))})
    return pub, ours, pd.DataFrame(car), pop


# ── A10 ──────────────────────────────────────────────────────────────────────
def a10(X):
    R, a = X.R, X.a
    # mirror pairs per scan, and the full-mainland check (slide 35)
    mp = X.thr[['case_set', 'cohort', 'class', 'markers', 'mirror_pairs', 'tests']].rename(columns={'case_set': 'comparison'})
    mp = mp[mp['class'] == 'residue']
    d = R / 'all' / 'full_mainland' / '01.markers' / 'residue'
    ts = pd.read_csv(d / 'tests.residue.tsv', sep='\t', dtype=str, keep_default_na=False)
    ts['pos'] = ts.marker.str.rsplit(':', n=1).str[0]
    A = X.A[(X.A.case_set == 'all') & (X.A.cohort == 'full_mainland') & (X.A.model == 'random') &
            (X.A['class'] == 'residue')].set_index('marker')
    m = ts[ts.mirror_of != '']
    pairs = m[m.marker < m.mirror_of]
    dp = max(abs(np.log10(A.at[r.marker, 'p']) - np.log10(A.at[r.mirror_of, 'p'])) for r in pairs.itertuples())
    dor = max(abs(np.log(A.at[r.marker, 'or']) + np.log(A.at[r.mirror_of, 'or'])) for r in pairs.itertuples())
    two = ts.groupby('pos').marker.apply(list)
    two = two[two.str.len() == 2]
    near = [v for k, v in two.items() if k not in set(m.pos)]
    with tempfile.TemporaryDirectory(dir='.') as tmp:
        subprocess.run([a.plink2, '--bfile', str(d / 'residue'), '--export', 'A', '--out', str(Path(tmp) / 'g')],
                       check=True, stdout=subprocess.DEVNULL)
        raw = pd.read_csv(Path(tmp) / 'g.raw', sep=r'\s+', na_values=['NA'])
    x = raw.iloc[:, 6:]
    x.columns = [c.rsplit('_', 1)[0] for c in x.columns]
    breaks = [int(((x[u] + x[v])[x[u].notna() & x[v].notna()] != 2).sum()) for u, v in near]
    mirror = {'columns': int(len(ts)), 'tests': int((ts.test_id == ts.marker).sum()), 'pairs': int(len(pairs)),
              'two_residue_positions': int(len(two)), 'near_not_merged': int(len(near)),
              'break_min': int(min(breaks)) if breaks else 0, 'break_max': int(max(breaks)) if breaks else 0,
              'samples': int(len(x)), 'max_dlog10p': float(dp), 'max_logor_sum': float(dor)}
    # fixed against random (slide 36), over tests
    conc, only = [], []
    for cs in COMPARISONS:
        for c in COHORTS:
            for k in ('allele', 'residue'):
                d_ = X.A[(X.A.case_set == cs) & (X.A.cohort == c) & (X.A['class'] == k) & (X.A.test_id == X.A.marker)]
                p = d_.pivot_table(index='marker', columns='model', values='p', aggfunc='first').dropna()
                b = d_.pivot_table(index='marker', columns='model', values='beta', aggfunc='first').reindex(p.index)
                nl = -np.log10(p)
                conc.append({'comparison': cs, 'cohort': c, 'test_type': k, 'tests': len(p),
                             'r_nlp': float(np.corrcoef(nl.fixed, nl.random)[0, 1]),
                             'r_beta': float(np.corrcoef(b.fixed, b.random)[0, 1]),
                             'median_ratio': float(np.median(nl.random / nl.fixed.where(nl.fixed > 0)))})
                gec = X.t(cs, c, k, 'p_gec')
                g = (p < gec)
                for mk in g.index[g.fixed != g.random]:
                    only.append({'comparison': cs, 'cohort': c, 'test_type': k, 'marker': mk,
                                 'significant_in': 'fixed' if g.at[mk, 'fixed'] else 'random',
                                 'p_fixed': float(p.at[mk, 'fixed']), 'p_random': float(p.at[mk, 'random']), 'gec': gec})
    # platform by DVT status within cases (slide 38), and by group
    case_sets = pd.read_csv(R / '_reference/case_sets.tsv', sep='\t', dtype=str)
    cram = pd.read_csv(a.cram, dtype=str)
    plat, chi = [], []
    for c in COHORTS:
        ph = pd.read_csv(Path(a.model_inputs) / c / 'phenotype/pheno.tsv', sep='\t', dtype=str)
        mm = ph.merge(cram[['ID_JHRPv6', 'WGS_Platform']], left_on='IID', right_on='ID_JHRPv6', how='left')
        if mm.WGS_Platform.isna().any():
            die(f'A10 {c}: samples without a platform')
        cases = case_sets[case_sets.IID.isin(ph.loc[ph.PHENO1 == '1', 'IID'])].merge(
            cram[['ID_JHRPv6', 'WGS_Platform']], left_on='IID', right_on='ID_JHRPv6', how='left')
        full = pd.crosstab(cases.WGS_Platform, cases.dvt).reindex(columns=['yes', 'no', 'unknown'], fill_value=0)
        ctrl = mm[mm.PHENO1 == '0'].WGS_Platform.value_counts()
        for pf in sorted(set(full.index) | set(ctrl.index)):
            plat.append({'cohort': c, 'platform': pf, 'controls': int(ctrl.get(pf, 0)),
                         'cases_dvt_yes': int(full.at[pf, 'yes']) if pf in full.index else 0,
                         'cases_dvt_no': int(full.at[pf, 'no']) if pf in full.index else 0,
                         'cases_dvt_unknown': int(full.at[pf, 'unknown']) if pf in full.index else 0})
        r = stats.chi2_contingency(full[['yes', 'no']].values)
        chi.append({'cohort': c, 'chi2': float(r[0]), 'df': int(r[2]), 'p': float(r[1])})
    plat = pd.DataFrame(plat)
    plat['cases_all'] = plat.cases_dvt_yes + plat.cases_dvt_no + plat.cases_dvt_unknown
    return mirror, mp, pd.DataFrame(conc), pd.DataFrame(only), plat, pd.DataFrame(chi)


def pos57_panels(X):
    """DQB1 position 57 by P group: our full-mainland cases and controls against ToMMo 61KJPN and 1KG JPT (slides
    39-40); the deck's pos57_panels, ported."""
    TC, TN, TP = X.TC, X.TN, X.TP
    r57_ = dqb1_res57(X)
    r57 = lambda g: r57_(g) or next((r57_(g + x) for x in 'QLSCAN' if r57_(g + x)), None)  # noqa: E731
    tom, tnull = TP.tommo_counts()
    kg, _ = TP.kg_counts()
    T, K = tom['DQB1'], kg['DQB1']
    td = pd.read_csv(TC.TOMMO)
    tnames = {}
    for a_ in td[td.gene == 'DQB1'].allele:
        if not a_.rstrip('P').endswith('N'):
            tnames.setdefault(TN.panel_pg(a_), []).append(a_)
    kd = pd.read_csv(TC.KG, sep='\t', dtype=str)
    kd = kd[kd.Population == TC.KG_POPULATION]
    kcalls, kmiss, kdrop = {}, 0, 0
    for a_ in pd.concat([kd['HLA-DQB1 1'], kd['HLA-DQB1 2']]):
        s_ = '' if pd.isna(a_) else str(a_).strip()
        if s_ in ('', '-', 'NA'):
            kmiss += 1
            continue
        m_ = TP.KG_CLEAN.match(s_)
        if not m_ or TN.is_null(f'DQB1*{m_.group(1)}'):
            kdrop += 1
            continue
        nm_ = f'DQB1*{m_.group(1)}'
        kcalls.setdefault(TN.pg_of(nm_), {}).setdefault(nm_, 0)
        kcalls[TN.pg_of(nm_)][nm_] += 1
    if {g: sum(v.values()) for g, v in kcalls.items()} != K:
        die('A10: the 1KG JPT calls do not reproduce the typing component\'s counts')
    ph = X.pheno('all', 'full_mainland').set_index('IID')
    al = pd.read_csv(X.TY / 'results/05.exports/allele_dosage.tsv', sep='\t', dtype={'sample_id': str},
                     usecols=lambda k: k == 'sample_id' or k.startswith('DQB1*')).set_index('sample_id').loc[ph.index]
    typed = al.loc[al.notna().any(axis=1)]
    ours, n2 = {}, {}
    for g, k in ((1, 'case'), (0, 'control')):
        d = typed.loc[[i for i in typed.index if ph.loc[i, 'PHENO1'] == g]]
        n2[k] = int(2 * len(d))
        ours[k] = {}
        for a_ in d.columns:
            pg = TN.pg_of(a_)
            ours[k][pg] = ours[k].get(pg, 0) + float(d[a_].sum())
    members = {}
    for a_ in typed.columns:
        if typed[a_].sum() > 0:
            members.setdefault(TN.pg_of(a_), []).append(a_)
    groups = set(T) | set(K) | set(ours['case']) | set(ours['control'])
    if [g for g in groups if r57(g) is None]:
        die('A10: P groups without a residue at 57')
    num = lambda g: [int(x) for x in g.split('*')[1].split(':')]  # noqa: E731
    vint = TP.tommo_vintage()
    rel = lambda v: f'{v // 1000}.{(v % 1000) // 10:02d}.{v % 10}'  # noqa: E731
    rows = []
    for g in sorted(groups, key=lambda g: ('VDAS'.index(r57(g)), num(g))):
        mem = sorted(members.get(g, []), key=num)
        rels = [TN.release_of(a_) for a_ in mem]
        rows.append({'p_group': TN.pg_label(g), 'residue_57': r57(g),
                     'our_alleles': ', '.join(f'{a_.split("*")[1]} ({int(typed[a_].sum())})' for a_ in mem),
                     'n_our_alleles': len(mem), 'copies_case': int(ours['case'].get(g, 0)),
                     'copies_control': int(ours['control'].get(g, 0)), 'tommo': int(T.get(g, 0)), 'kg': int(K.get(g, 0)),
                     'tommo_names': ', '.join(sorted(tnames.get(g, []))),
                     'kg_calls': ', '.join(f'{n_} ({c_})' for n_, c_ in sorted(kcalls.get(g, {}).items())),
                     'named_later_than_tommo': bool(mem and all((r_ or 0) > vint for r_ in rels)),
                     'first_release': ', '.join(rel(r_) for r_ in rels if r_)})
    pg = pd.DataFrame(rows)
    nT, nK = int(sum(T.values())), int(sum(K.values()))
    totals = {'cases_chr': n2['case'], 'controls_chr': n2['control'], 'tommo_chr': nT, 'kg_chr': nK,
              'tommo_vintage': rel(vint), 'kg_vintage': rel(TP.kg_vintage()), 'tommo_null_copies': int(tnull['DQB1']),
              'kg_people': int(len(kd)), 'kg_missing': kmiss, 'kg_dropped': kdrop,
              'dictionary_ours': re.search(r'\d+\.\d+\.\d+', TC.HLAHD_DICT.name).group(0)}
    for k, col in (('case', 'copies_case'), ('control', 'copies_control'), ('tommo', 'tommo'), ('kg', 'kg')):
        if int(pg[col].sum()) != totals[{'case': 'cases_chr', 'control': 'controls_chr', 'tommo': 'tommo_chr',
                                         'kg': 'kg_chr'}[k]]:
            die(f'A10: the P groups of {k} do not add up to its chromosomes')
    return pg, totals


# ── A11 / A12 ────────────────────────────────────────────────────────────────
def a11(X):
    A = X.A
    keys = ['case_set', 'cohort', 'class', 'marker', 'gene', 'position', 'residue', 'test_id', 'mirror_of']
    r = A[A.model == 'random'].set_index(keys)
    f = A[A.model == 'fixed'].set_index(keys)
    out = r[['n_case', 'n_control', 'called_case', 'called_control', 'copies_case', 'copies_control', 'freq_case',
             'freq_control', 'mac', 'missing_case', 'missing_control', 'or', 'l95', 'u95', 'p', 'p_spa', 'errcode']].copy()
    out = out.rename(columns={c: f'r_{c}' for c in ('or', 'l95', 'u95', 'p', 'p_spa', 'errcode')})
    ff = f[['or', 'l95', 'u95', 'p', 'errcode', 'n_case', 'n_control']].rename(
        columns={c: f'f_{c}' for c in ('or', 'l95', 'u95', 'p', 'errcode', 'n_case', 'n_control')})
    out = out.join(ff, how='outer').reset_index()
    if len(out) != len(A) // 2:
        die('A11: the random and fixed rows do not pair up')
    out['gec'] = [X.t(cs, c, k, 'p_gec') for cs, c, k in zip(out.case_set, out.cohort, out['class'])]
    out['bonf'] = [X.t(cs, c, k, 'p_bonf') for cs, c, k in zip(out.case_set, out.cohort, out['class'])]
    out['r_sig_gec'] = out.r_p < out.gec
    out['r_sig_bonf'] = out.r_p < out.bonf
    out['f_sig_gec'] = out.f_p < out.gec
    # n/3: cohorts significant at GEC (random model) with the direction of the full mainland
    n3 = out[out.r_sig_gec].groupby(['case_set', 'class', 'marker']).size()
    out['n3_gec'] = [int(n3.get((cs, k, m), 0)) for cs, k, m in zip(out.case_set, out['class'], out.marker)]
    out = out.rename(columns={'case_set': 'comparison'})
    order = {cs: i for i, cs in enumerate(COMPARISONS)}
    out['_o'] = out.comparison.map(order)
    out['_c'] = out.cohort.map({c: i for i, c in enumerate(COHORTS)})
    out = out.sort_values(['_o', '_c', 'class', 'gene', 'marker']).drop(columns=['_o', '_c'])
    return out


def a12(X):
    O = X.O.rename(columns={'case_set': 'comparison'}).copy()
    n3 = O[O.sig_bonf.astype(bool)].groupby(['comparison', 'unit']).size()
    O['n3'] = [int(n3.get((cs, u), 0)) for cs, u in zip(O.comparison, O.unit)]
    O['_o'] = O.comparison.map({cs: i for i, cs in enumerate(COMPARISONS)})
    O['_c'] = O.cohort.map({c: i for i, c in enumerate(COHORTS)})
    return O.sort_values(['_o', '_c', 'gene', 'unit']).drop(columns=['_o', '_c', 'class', 'pooled'], errors='ignore')


def main():
    a = parse_args()
    o = Path(a.out)
    X = Assoc(a)
    gw = {k: v for f in a.gwas_lambda for k, v in json.loads(Path(f).read_text()).items()}
    if sorted(gw) != sorted(COHORTS):
        die(f'gwas lambda: cohorts {sorted(gw)}')
    samples, tg, methods, grm = a1(X)
    write_tsv(samples, o / 'a1_samples.tsv')
    write_tsv(tg, o / 'a1_tests_by_gene.tsv')
    write_tsv(methods, o / 'a1_methods.tsv')
    sc, pos = a2(X)
    write_tsv(sc, o / 'a2_scans.tsv')
    write_tsv(pos, o / 'a2_positions.tsv')
    gc = gec_check(X)
    write_tsv(gc, o / 'a2_gec_check.tsv')
    scans, top = a3(X)
    write_tsv(scans, o / 'a3_scans.tsv')
    write_tsv(top, o / 'a3_top.tsv')
    steps, after, b114 = a4(X)
    write_tsv(steps, o / 'a4_steps.tsv')
    write_tsv(after, o / 'a4_after.tsv')
    leads, lg, given, link, kept = a5(X)
    write_tsv(leads, o / 'a5_leads.tsv')
    write_tsv(lg, o / 'a5_lead_given.tsv')
    write_tsv(given, o / 'a5_given_lead.tsv')
    write_tsv(link, o / 'a5_link.tsv')
    lam, multi = a6(X, gw)
    write_tsv(lam, o / 'a6_lambda.tsv')
    geno, v57, res57, panels, carriers = a7(X)
    write_tsv(geno, o / 'a7_genotypes.tsv')
    write_tsv(v57, o / 'a7_v57.tsv')
    write_tsv(res57, o / 'a7_residues.tsv')
    strata, pos3, lrt, nodvt, ld, g57, rev, spa, spa_counts, checked = a8(X)
    write_tsv(strata, o / 'a8_strata.tsv')
    write_tsv(pos3, o / 'a8_nodvt_positions.tsv')
    write_tsv(lrt, o / 'a8_lrt_vs_score.tsv')
    write_tsv(nodvt, o / 'a8_nodvt_dqb1.tsv')
    write_tsv(g57, o / 'a8_given57v.tsv')
    write_tsv(rev, o / 'a8_57v_given_87y.tsv')
    write_tsv(spa, o / 'a8_spa_lost.tsv')
    pub, ours, car, pop = a9(X)
    write_tsv(pub, o / 'a9_published.tsv')
    write_tsv(ours, o / 'a9_ours.tsv')
    write_tsv(car, o / 'a9_carriers.tsv')
    mirror, mp, conc, only, plat, chi = a10(X)
    write_tsv(mp, o / 'a10_mirror.tsv')
    write_tsv(conc, o / 'a10_concordance.tsv')
    write_tsv(only, o / 'a10_one_model_only.tsv')
    write_tsv(plat, o / 'a10_platform.tsv')
    write_tsv(chi, o / 'a10_platform_chi2.tsv')
    pg57, pg_tot = pos57_panels(X)
    write_tsv(pg57, o / 'a10_pos57_pgroups.tsv')
    d11 = a11(X)
    write_tsv(d11, o / 'a11_single_markers.tsv')
    d12 = a12(X)
    write_tsv(d12, o / 'a12_positions.tsv')
    write_json({'grm_markers': grm, 'b114': b114, 'gwas_leads_kept': kept, 'positions_multi_df': multi,
                'panels_v57': panels, 'carrier_alleles_57v': carriers, 'ld_57_87_without_dvt': ld,
                'spa_calls': spa_counts, 'given57v_marginal_rows_checked': checked, 'pop_reported': pop,
                'mirror_full': mirror, 'pos57_totals': pg_tot, 'gwas_lambda': gw,
                'run': {'manifest_component': X.manifest['component']}}, o / 'assoc_facts.json')
    log(f'association tables written: A11 {len(d11)} rows, A12 {len(d12)} rows')


if __name__ == '__main__':
    main()
