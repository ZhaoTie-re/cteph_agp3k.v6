#!/usr/bin/env python3
"""DQB1:57:V on the SNP arrays, imputed with HIBAG (report stage, run by ARRAY_EVIDENCE after array_cohorts.py).

The array check of DQB1:57:V uses HIBAG's published Asian two-field model (user 2026-11-05), not a model trained on
our own typing: the model never saw the WGS HLA-HD calls, so the comparison is independent in cases and controls
alike. It stays a within-person comparison (array imputation against WGS typing, same people); no external allele
frequency is read.

Input: the xMHC (chr6:28-34 Mb, hg38) of the QC'd array set for the full-mainland people with arrays (DNBSeq-T7
cases on ASA24v1, HiSeq X controls on JSA24v1; array_crosscheck's keep list), restricted to SNPs called in >= 98% of
the cases AND of the controls, so both groups are imputed from the same SNPs. HIBAG imputes each person on its own,
so one run serves the three nested cohorts; each cohort then reads its own members (fixed-model set).

Per cohort (<cohorts>/<cohort>/07.hibag/summary.json; counts and estimates only):
  two-field DQB1 allele concordance with the WGS typing (per allele copy), cases and controls;
  57 V: residue 57 of each two-field allele from hla.typing.v2's IMGT reference; best-guess agreement, r2 of the
  expected dosage (posterior-weighted) with the typed dosage, frequencies, typed-minus-imputed residual cases vs
  controls (Welch t); the association (case ~ dosage + sex + 10 PCs, standardized), typed and HIBAG dosage.
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as K  # noqa: E402

FIXED = 'cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.fixed_ready'
ARRAY_SET = K.PROJECT.parent / 'cteph_agp3k.v5/array/results/10.variant_qc_filtered/cteph_agp3k.array.sqc.vqc'
ASSETS = K.ROOT / 'assets' / 'hibag'
RSCRIPT = '/home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/Rscript'
HLA = K.PROJECT / 'hla.typing.v2' / 'results'
COVARS = ['SEX'] + [f'PC{i}_AVG' for i in range(1, 11)]
V57 = 'DQB1:57:V'
XMHC = ('6', 28_000_000, 34_000_000)
MAX_MISS = 0.02


def ids(p):
    return set(pd.read_csv(p, sep='\t', header=None, dtype=str)[1])


def run(cmd):
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if r.returncode:
        sys.exit(f'ABORT: {cmd[0]} failed\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}')
    return r


def residue57():
    ref = pd.read_csv(HLA / '05.exports' / 'residue_reference' / 'DQB1.tsv', sep='\t', dtype=str, usecols=['allele', '57'])
    ref['two'] = ref.allele.str.replace('DQB1*', '', regex=False).str.split(':').str[:2].str.join(':')
    ok = ref[ref['57'].str.fullmatch(r'[A-Z]')]
    if (ok.groupby('two')['57'].nunique() > 1).any():
        sys.exit('ABORT: a two-field DQB1 allele with more than one residue at 57')
    return ok.groupby('two')['57'].first()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--array', required=True, help='check/array_crosscheck/results (or the component folder)')
    ap.add_argument('--plink2', required=True)
    ap.add_argument('--cohorts', required=True, help='array_cohorts.py output folder')
    ap.add_argument('--work', default='hibag_work')
    ap.add_argument('--threads', type=int, default=4)
    a = ap.parse_args()
    w, out = Path(a.work), Path(a.cohorts)
    w.mkdir(parents=True, exist_ok=True)
    # --array is check/array_crosscheck/results in the pipeline (as for array_cohorts.py); the component folder works too
    xc = Path(a.array)
    inp = (xc / '01.inputs') if (xc / '01.inputs').is_dir() else (xc / 'results' / '01.inputs')
    cases, controls = ids(inp / 'keep.cases.txt'), ids(inp / 'keep.controls.txt')
    md5 = (ASSETS / 'MD5SUMS').read_text().split()[0]
    got = run(['md5sum', ASSETS / 'Asian-HLA4-hg38.RData']).stdout.split()[0]
    if got != md5:
        sys.exit('ABORT: the HIBAG model does not match assets/hibag/MD5SUMS')

    # the xMHC SNPs called in both groups, for every person with arrays
    (w / 'xmhc.range').write_text(f'{XMHC[0]}\t{XMHC[1]}\t{XMHC[2]}\txmhc\n')
    for g in ('cases', 'controls'):
        run([a.plink2, '--bfile', ARRAY_SET, '--keep', inp / f'keep.{g}.txt', '--extract', 'range', w / 'xmhc.range',
             '--missing', 'variant-only', '--out', w / f'miss.{g}', '--threads', a.threads])
    mc = pd.read_csv(w / 'miss.cases.vmiss', sep='\t')
    mk = pd.read_csv(w / 'miss.controls.vmiss', sep='\t')
    both = set(mc[mc.F_MISS <= MAX_MISS].ID) & set(mk[mk.F_MISS <= MAX_MISS].ID)
    (w / 'both.ids').write_text('\n'.join(sorted(both)) + '\n')
    run([a.plink2, '--bfile', ARRAY_SET, '--keep', inp / 'keep.txt', '--extract', w / 'both.ids', '--snps-only',
         'just-acgt', '--max-alleles', '2', '--make-bed', '--out', w / 'arr.xmhc', '--threads', a.threads])
    run([RSCRIPT, Path(__file__).resolve().parent / 'array_hibag.R', w / 'arr.xmhc', ASSETS / 'Asian-HLA4-hg38.RData',
         ASSETS / 'Rlib', w / 'hibag_DQB1', a.threads])
    runinfo = dict(l.split('\t') for l in (w / 'hibag_DQB1.run.tsv').read_text().splitlines())
    runinfo['xmhc_snps_both_groups'] = str(len(both))

    # HIBAG calls -> two-field alleles and the 57 V dosage (best pair; expected over every pair)
    r57 = residue57()

    def is_v(al):
        al = ':'.join(str(al).replace('DQB1*', '').split(':')[:2])
        return np.nan if al not in r57.index else float(r57[al] == 'V')
    best = pd.read_csv(w / 'hibag_DQB1.best.tsv', sep='\t', dtype={'sample.id': str, 'allele1': str, 'allele2': str})
    best = best.set_index('sample.id')
    pp = pd.read_csv(w / 'hibag_DQB1.postprob.tsv', sep='\t', dtype={'sample_id': str})
    pp[['a1', 'a2']] = pp.pair.str.split('/', expand=True)
    pp['v'] = pp.a1.map(is_v) + pp.a2.map(is_v)
    unmapped = sorted((set(pp.a1) | set(pp.a2)) - set(r57.index))   # alleles without a residue at 57 (none expected)
    ok = pp.dropna(subset=['v'])
    exp_v = (ok.v * ok.prob).groupby(ok.sample_id).sum() / ok.prob.groupby(ok.sample_id).sum()
    H = pd.DataFrame({'h1': best.allele1, 'h2': best.allele2, 'h_prob': best.prob})
    H['h_v_best'] = H.h1.map(is_v) + H.h2.map(is_v)
    H['h_v_exp'] = exp_v.reindex(H.index)
    g = pd.read_csv(HLA / '01.calls' / 'genotypes.tsv', sep='\t', dtype=str, usecols=['sample_id', 'DQB1_1', 'DQB1_2'])
    g = g.set_index('sample_id')
    two = lambda x: ':'.join(str(x).replace('DQB1*', '').split(':')[:2]) if isinstance(x, str) else np.nan  # noqa: E731
    H = H.join(pd.DataFrame({'w1': g.DQB1_1.map(two), 'w2': g.DQB1_2.map(two)}), how='left')

    def matched(r):
        wa, ha = [r.w1, r.w2], [r.h1, r.h2]
        if any(pd.isna(x) for x in wa + ha):
            return np.nan
        m = 0
        for x in wa:
            if x in ha:
                ha.remove(x)
                m += 1
        return m
    H['matched'] = H.apply(matched, axis=1)

    for c in K.COHORTS:
        mi = K.MODEL_INPUTS / c
        fixed = set(pd.read_csv(mi / 'genotype' / 'fixed_model' / f'{FIXED}.fam', sep=r'\s+', header=None, dtype=str)[1])
        typed = pd.read_csv(HLA / '06.assoc_inputs' / c / 'residue_dosage.tsv', sep='\t', usecols=['sample_id', V57],
                            dtype={'sample_id': str}).set_index('sample_id')[V57]
        d = H[H.index.isin(fixed) & H.index.isin(cases | controls)].copy()
        d['w_v'] = typed.reindex(d.index)
        d = d.dropna(subset=['w_v', 'h_v_exp'])
        is_case = d.index.isin(sorted(cases))
        rk, rc = (d.w_v - d.h_v_exp)[is_case], (d.w_v - d.h_v_exp)[~is_case]
        tt = stats.ttest_ind(rk, rc, equal_var=False)
        S = {'tool': f"HIBAG {runinfo['hibag_version']}, published Asian two-field model (hg38, matched by position)",
             'run': runinfo, 'model_md5': md5, 'n_cases': int(is_case.sum()), 'n_controls': int((~is_case).sum())}
        for lab, sel in (('cases', is_case), ('controls', ~is_case)):
            x = d[sel]
            S[lab] = {'allele_concordance': float(x.matched.sum() / (2 * x.matched.notna().sum())),
                      'v57_best_agreement': float((x.h_v_best == x.w_v).mean()),
                      'v57_r2': float(np.corrcoef(x.w_v, x.h_v_exp)[0, 1] ** 2),
                      'v57_freq_typed': float(x.w_v.mean() / 2), 'v57_freq_hibag': float(x.h_v_exp.mean() / 2),
                      'posterior_median': float(x.h_prob.median())}
        S['residual'] = {'mean_cases': float(rk.mean()), 'mean_controls': float(rc.mean()),
                         'welch_t': float(tt.statistic), 'p': float(tt.pvalue)}
        cov = pd.read_csv(mi / 'covariates' / 'bbj_mainland_pc.sex.tsv', sep='\t', dtype={'IID': str}).set_index('IID')
        base = cov.reindex(d.index)[COVARS]
        base = (base - base.mean()) / base.std()
        y = is_case.astype(int)

        def logit(dose):
            X = sm.add_constant(pd.concat([pd.Series(dose.values, index=d.index, name='d'), base], axis=1))
            f = sm.Logit(y, X).fit(disp=0)
            b, se = f.params['d'], f.bse['d']
            return {'or': float(np.exp(b)), 'l95': float(np.exp(b - 1.96 * se)), 'u95': float(np.exp(b + 1.96 * se)),
                    'p': float(f.pvalues['d'])}
        S['assoc_typed'] = logit(d.w_v)
        S['assoc_hibag'] = logit(d.h_v_exp)
        S['unmapped_alleles'] = unmapped
        o = out / c / '07.hibag'
        o.mkdir(parents=True, exist_ok=True)
        (o / 'summary.json').write_text(json.dumps(S, indent=2) + '\n')
        print(c, json.dumps({k: S[k] for k in ('n_cases', 'n_controls', 'assoc_typed', 'assoc_hibag')}))


if __name__ == '__main__':
    main()
