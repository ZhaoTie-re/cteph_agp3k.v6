#!/usr/bin/env python3
"""The SNP-array cross-check in each of the three cohorts (report stage, run by ARRAY_EVIDENCE before
array_evidence.py).

check/array_crosscheck compared WGS with SNP arrays in the full-mainland people who have both (unrelated set: 299
DNBSeq-T7 cases typed on ASA24v1, 1,384 HiSeq X controls typed on JSA24v1; NovaSeq and DNBSeq-G400RS genomes have
no array). The cohorts are nested, so narrow and intermediate mainland are the same people's data with fewer
people: this script takes array_crosscheck's inputs (read only), keeps each cohort's unrelated members (its
fixed-model set) and re-runs every step of array_crosscheck (scripts 01, 03, 04, 05, 06) with that cohort's
phenotype and covariates. The steps are the same code path in every cohort: the shared sites, the logistic model
(sex + 10 BBJ-mainland PCs, Firth fallback), the proxy chosen by r2 in the cohort's controls, the DQB1:57:V dosage
predicted from array SNPs by ridge regression trained on the cohort's controls (out of fold), the rare drivers typed
on ASA24v1, and the imputed arrays (Minimac4 v4.1.4 against the AGP3K whole-genome reference panel, which holds the
controls' own genomes, so only the cases' imputed genotypes are independent of WGS).

Gate G10: in full mainland the re-run must reproduce array_crosscheck's published tables (array.gate.json); the
script exits 1 otherwise, before any cohort is compared.

Writes, per cohort, the tables array_evidence.py reads, in array_crosscheck's layout (counts and estimates only):
  <out>/<cohort>/04.compare/{signals.tsv, signals.direct_on_array.tsv, summary.json}
  <out>/<cohort>/05.rare_direct/rare_direct.tsv
  <out>/<cohort>/06.imputed/{leads.tsv, drivers.tsv, recovery_by_group.tsv, concordance.tsv}
  <out>/array.gate.json
Work files (plink sets and exports keyed by sample ID) stay in <work>, which the process does not publish.
"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as K  # noqa: E402

FIXED = 'cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.fixed_ready'
ARRAY_SET = K.PROJECT.parent / 'cteph_agp3k.v5/array/results/10.variant_qc_filtered/cteph_agp3k.array.sqc.vqc'
IMPUTED = K.PROJECT.parent / 'cteph_agp3k/array_imputed_rev1/results/04.sample_qc/cteph_agp3k.imputed_array.sqc'
EUROPEAN_LOCI = K.PROJECT / 'tmp/ppt/overview/european_loci.json'
COVARS = ['SEX'] + [f'PC{i}_AVG' for i in range(1, 11)]
GLM = ['--pheno-name', 'PHENO1', '--1', '--covar-name', 'SEX', 'PC1_AVG-PC10_AVG', '--covar-variance-standardize',
       '--glm', 'firth-fallback', 'hide-covar',
       'cols=chrom,pos,ref,alt,a1freq,a1freqcc,firth,test,nobs,orbeta,se,ci,tz,p']
DRIVERS = {'chr4:76306990:G:C': 'STBD1 c.220+1G>C', 'chr4:76309559:G:T': 'STBD1 p.Trp212Cys',
           'chr4:76309487:C:G': 'STBD1 p.Ser188Arg', 'chr12:112890858:C:T': 'RPH3A p.Arg544Cys',
           'chr11:59576952:G:A': 'OSBP p.Pro712Ser', 'chr11:59594182:C:T': 'OSBP p.Arg462Gln',
           'chr9:5073770:G:T': 'JAK2 p.Val617Phe', 'chr9:5073770:G:A': 'JAK2 c.1849G>A (other allele)'}
DIRECT_LABEL = {'chr12:112890858': 'RPH3A p.Arg544Cys', 'chr11:59576952': 'OSBP p.Pro712Ser',
                'chr9:5073770': 'JAK2 p.Val617Phe'}
BINS = [(0, 0, '0'), (1, 1, '1'), (2, 3, '2-3'), (4, 10, '4-10'), (11, 10**6, '>=11')]
V57 = 'DQB1:57:V'


class Tools:
    def __init__(self, plink2, plink1, threads, memory):
        self.p2, self.p1, self.t, self.m = plink2, plink1, str(threads), str(memory)

    def run(self, args, log_ok=None, allow_fail=False):
        r = subprocess.run(args, capture_output=True, text=True)
        if r.returncode and not allow_fail:
            sys.exit(f'ABORT: {" ".join(map(str, args[:6]))} ... exit {r.returncode}\n{r.stdout[-2000:]}{r.stderr[-2000:]}')
        if log_ok:
            log, pat = log_ok
            if not re.search(pat, Path(log).read_text()):
                sys.exit(f'ABORT: "{pat}" not in {log} (plink does not stop when a named covariate is absent)')
        return r

    def plink2(self, *a, **kw):
        return self.run([self.p2, *map(str, a), '--threads', self.t, '--memory', self.m], **kw)

    def plink1(self, *a, **kw):
        return self.run([self.p1, *map(str, a), '--threads', self.t, '--memory', self.m], **kw)


def glm(path):
    d = pd.read_csv(path, sep='\t').rename(columns={'#CHROM': 'CHROM'})
    return d[d.TEST == 'ADD'].set_index('ID')


def raw(path):
    """a plink .raw file: counted alleles per sample (index IID), columns = variant IDs, and each one's counted allele"""
    r = pd.read_csv(path, sep=r'\s+').drop(columns=['FID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE']).set_index('IID')
    r.index = r.index.astype(str)
    allele = {c.rsplit('_', 1)[0]: c.rsplit('_', 1)[1] for c in r.columns}
    r.columns = [c.rsplit('_', 1)[0] for c in r.columns]
    return r, allele


def alt_counts(path):
    """a .raw file as ALT-allele counts or dosages, columns = variant IDs (chr:pos:ref:alt)"""
    r = pd.read_csv(path, sep='\t').drop(columns=['FID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'])
    r['IID'] = r.IID.astype(str)
    r = r.set_index('IID')
    out = {}
    for c in r.columns:
        vid, counted = c.rsplit('_', 1)
        out[vid] = r[c] if counted == vid.split(':')[3] else 2 - r[c]
    return pd.DataFrame(out)


def lam(p):
    return float(np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1))


def binlab(n):
    return next(lab for lo, hi, lab in BINS if lo <= n <= hi)


# ── steps 01 and 03: people, association, missingness, HWE, concordance ─────────────────────────────────────────
def prepare(T, xc, cohort, w):
    I = xc / '01.inputs'
    mi = K.MODEL_INPUTS / cohort
    pheno, cov = mi / 'phenotype' / 'pheno.tsv', mi / 'covariates' / 'bbj_mainland_pc.sex.tsv'
    fixed = set(pd.read_csv(mi / 'genotype' / 'fixed_model' / f'{FIXED}.fam', sep=r'\s+', header=None, dtype=str)[1])
    keep = pd.read_csv(I / 'keep.txt', sep='\t', header=None, dtype=str)
    keep = keep[keep[1].isin(fixed)]
    ph = pd.read_csv(pheno, sep='\t', dtype=str).set_index('IID').PHENO1
    w.mkdir(parents=True, exist_ok=True)
    keep.to_csv(w / 'keep.txt', sep='\t', header=False, index=False)
    for g, v in (('cases', '1'), ('controls', '0')):
        keep[keep[1].map(ph) == v].to_csv(w / f'keep.{g}.txt', sep='\t', header=False, index=False)
    glm_args = ['--pheno', pheno, *GLM[:3], '--covar', cov, *GLM[3:]]
    ok = (r'11 covariates loaded',)
    for s in ('wgs', 'arr'):
        T.plink2('--bfile', I / f'{s}.shared', '--keep', w / 'keep.txt', *glm_args, '--out', w / s,
                 log_ok=(w / f'{s}.log', ok[0]))
        for g in ('cases', 'controls'):
            T.plink2('--bfile', I / f'{s}.shared', '--keep', w / f'keep.{g}.txt', '--missing', 'variant-only',
                     '--out', w / f'{s}.{g}')
        T.plink2('--bfile', I / f'{s}.shared', '--keep', w / 'keep.controls.txt', '--hardy', '--out', w / f'{s}.controls')
    wgs_full = K.MODEL_INPUTS / 'full_mainland' / 'genotype' / 'fixed_model' / FIXED
    T.plink2('--bfile', wgs_full, '--keep', w / 'keep.txt', '--extract', I / 'wgs_window.ids', *glm_args,
             '--out', w / 'wgs.window', log_ok=(w / 'wgs.window.log', ok[0]))
    T.plink2('--bfile', ARRAY_SET, '--keep', w / 'keep.txt', '--extract', 'range', I / 'lead_windows.range', *glm_args,
             '--out', w / 'arr.window', log_ok=(w / 'arr.window.log', ok[0]))
    for g in ('cases', 'controls'):
        for s in ('wgs', 'arr'):
            T.plink2('--bfile', I / f'{s}.shared', '--keep', w / f'keep.{g}.txt', '--make-bed', '--out', w / f'tmp.{s}.{g}')
        T.plink1('--bfile', w / f'tmp.wgs.{g}', '--bmerge', w / f'tmp.arr.{g}', '--merge-mode', '7',
                 '--out', w / f'wgs_vs_arr.{g}', allow_fail=True)
        for f in w.glob(f'tmp.*.{g}.*'):
            f.unlink()
    return pheno, cov


# ── step 04: genome-wide, signals, DQB1:57:V ─────────────────────────────────────────────────────────────────────
def eur(e, cohort):
    """a European lead's OR and P in this cohort, as the overview reports them, for the ALT allele of its ID"""
    if not e:
        return np.nan, np.nan
    f = e.get(cohort) or (e.get('recovered') or {}).get('direct', {}).get(cohort)
    if not f:
        return np.nan, np.nan
    return (1 / f['or_'] if e['minor'] == e['ref'] else f['or_']), f['p']


def compare(xc, cohort, w, cov_path, out):
    I = xc / '01.inputs'
    S = {}
    cases = set(pd.read_csv(w / 'keep.cases.txt', sep='\t', header=None, dtype=str)[1])
    controls = set(pd.read_csv(w / 'keep.controls.txt', sep='\t', header=None, dtype=str)[1])
    S['samples'] = {'cases': len(cases), 'controls': len(controls)}

    gw, ga = glm(w / 'wgs.PHENO1.glm.logistic.hybrid'), glm(w / 'arr.PHENO1.glm.logistic.hybrid')
    D = gw.join(ga, lsuffix='_w', rsuffix='_a', how='inner')
    flip = D.A1_w != D.A1_a
    D['z_a'] = np.where(flip, -D.Z_STAT_a, D.Z_STAT_a)
    D['or_a'] = np.where(flip, 1 / D.OR_a, D.OR_a)
    D['fca_a'] = np.where(flip, 1 - D.A1_CASE_FREQ_a, D.A1_CASE_FREQ_a)
    D['fco_a'] = np.where(flip, 1 - D.A1_CTRL_FREQ_a, D.A1_CTRL_FREQ_a)
    for t in ('wgs', 'arr'):
        for g in ('cases', 'controls'):
            m = pd.read_csv(w / f'{t}.{g}.vmiss', sep='\t').set_index('ID')
            D[f'miss_{t}_{g}'] = m.F_MISS.reindex(D.index)
        h = pd.read_csv(w / f'{t}.controls.hardy', sep='\t').set_index('ID')
        D[f'hwe_{t}'] = h.P.reindex(D.index)
    maf = np.minimum(D.A1_FREQ_w, 1 - D.A1_FREQ_w)
    qc = ((maf >= 0.01) & D.Z_STAT_w.notna() & D.Z_STAT_a.notna() & (D[['hwe_wgs', 'hwe_arr']].min(axis=1) >= 1e-6))
    S['sites'] = {'shared': int(len(D)), 'qc_pass': int(qc.sum()), 'allele_flipped': int(flip.sum())}
    for g, n in (('cases', len(cases)), ('controls', len(controls))):
        p = w / f'wgs_vs_arr.{g}.diff'
        dd = pd.read_csv(p, sep=r'\s+', usecols=[0]) if p.exists() and p.stat().st_size else pd.DataFrame({'SNP': []})
        D[f'disc_{g}'] = dd.SNP.value_counts().reindex(D.index).fillna(0)
        nonmiss = n * (1 - D[f'miss_wgs_{g}']) * (1 - D[f'miss_arr_{g}'])
        D[f'disc_rate_{g}'] = D[f'disc_{g}'] / nonmiss
        log = (w / f'wgs_vs_arr.{g}.log').read_text()
        m = re.search(r'(\d+) overlapping calls, (\d+) nonmissing in both filesets\.\s+(\d+) concordant, for a '
                      r'concordance rate of ([0-9]*\.[0-9]+)', log.replace('\n', ' '))
        S[f'concordance_all_shared_{g}'] = ({'nonmissing_both': int(m.group(2)), 'concordant': int(m.group(3)),
                                             'rate': float(m.group(4))} if m else None)
    Q = D[qc].copy()
    dz = Q.z_a - Q.Z_STAT_w
    S['genome_wide'] = {
        'lambda_wgs': lam(Q.P_w), 'lambda_array': lam(Q.P_a),
        'r_z': float(np.corrcoef(Q.Z_STAT_w, Q.z_a)[0, 1]),
        'slope_z_array_on_wgs': float(np.polyfit(Q.Z_STAT_w, Q.z_a, 1)[0]),
        'abs_dz_gt3': int((dz.abs() > 3).sum()), 'abs_dz_gt4': int((dz.abs() > 4).sum()),
        'wgs_p_lt_1e-4': int((Q.P_w < 1e-4).sum()),
        'wgs_p_lt_1e-4_array_same_direction': int(((Q.P_w < 1e-4) & (np.sign(Q.z_a) == np.sign(Q.Z_STAT_w))).sum()),
        'wgs_p_lt_1e-4_array_p_lt_1e-2': int(((Q.P_w < 1e-4) & (Q.P_a < 1e-2)).sum()),
        'median_disc_rate_cases': float(Q.disc_rate_cases.median()),
        'median_disc_rate_controls': float(Q.disc_rate_controls.median()),
        'mean_disc_rate_cases': float(Q.disc_rate_cases.mean()),
        'mean_disc_rate_controls': float(Q.disc_rate_controls.mean()),
        'sites_disc_rate_gt_1pct_cases': int((Q.disc_rate_cases > 0.01).sum()),
        'sites_disc_rate_gt_1pct_controls': int((Q.disc_rate_controls > 0.01).sum()),
        'r_absdz_vs_disc_rate_gap': float(stats.spearmanr(dz.abs(),
                                                          (Q.disc_rate_cases - Q.disc_rate_controls).abs())[0]),
    }

    # signals: each lead, its best array proxy (r2 in the cohort's controls), the association at the proxy
    WR, wal = raw(I / 'wgs.window.raw')
    gww, gaw = glm(w / 'wgs.window.PHENO1.glm.logistic.hybrid'), glm(w / 'arr.window.PHENO1.glm.logistic.hybrid')
    EU = {e['region']: e for e in json.loads(EUROPEAN_LOCI.read_text())}
    T = pd.read_csv(I / 'targets.tsv', sep='\t', dtype={'chrom': str})
    long = pd.read_csv(K.SAIGE_SIGNALS / 'saige_signals.long.tsv', sep='\t', low_memory=False)
    pub = long[long.cohort == cohort].drop_duplicates('variant_id').set_index('variant_id')
    ctrl_idx = WR.index.intersection(sorted(controls))
    rows = []
    for t in T[T.role == 'Lead'].itertuples():
        e = EU.get(t.signal)
        if t.variant_id not in WR:
            a_ = gaw.loc[t.variant_id] if t.variant_id in gaw.index else None
            fo = eur(e, cohort)[0]
            if a_ is not None and a_.A1 != t.variant_id.split(':')[3]:
                fo = 1 / fo
            rows.append(dict(signal=t.signal, kind=t.kind, lead=t.variant_id, rsid=t.rsid, lead_source=t.lead_source,
                             full_or=fo, full_p=eur(e, cohort)[1], lead_a1=a_.A1 if a_ is not None else None,
                             sub_arr_or_lead=a_.OR if a_ is not None else np.nan,
                             sub_arr_p_lead=a_.P if a_ is not None else np.nan,
                             sub_arr_freq_cases_lead=a_.A1_CASE_FREQ if a_ is not None else np.nan,
                             sub_arr_freq_controls_lead=a_.A1_CTRL_FREQ if a_ is not None else np.nan,
                             note='lead not in the WGS analysis set; typed on the array'))
            continue
        x = WR.loc[ctrl_idx, t.variant_id]
        win = Q[(Q.CHROM_w.astype(str) == t.chrom) & (Q.POS_w.between(t.pos - 500_000, t.pos + 500_000))]
        cand = [v for v in win.index if v in WR]
        if not cand:
            rows.append(dict(signal=t.signal, lead=t.variant_id, note="no QC'd array site within 500 kb"))
            continue
        rr = WR.loc[ctrl_idx, cand].corrwith(x)
        best = (rr ** 2).idxmax()
        r = rr[best]
        pr = Q.loc[best]
        ld = gww.loc[t.variant_id] if t.variant_id in gww.index else None
        lead_sign = 1 if (ld is not None and ld.A1 == wal[t.variant_id]) else -1
        proxy_sign = 1 if pr.A1_w == wal[best] else -1
        o = lambda v: v if np.sign(r) * lead_sign * proxy_sign > 0 else 1 / v  # noqa: E731,B023
        lead_a1 = ld.A1 if ld is not None else None
        fo = pub.OR.get(t.variant_id, eur(e, cohort)[0])
        fea = pub.EA.get(t.variant_id, t.variant_id.split(':')[3])
        fo = fo if (lead_a1 is None or fea == lead_a1) else 1 / fo
        near = win.assign(r2=(WR.loc[ctrl_idx, cand].corrwith(x) ** 2).reindex(win.index))
        near = near[near.r2 >= 0.2]
        rows.append(dict(
            signal=t.signal, kind=t.kind, lead=t.variant_id, rsid=t.rsid, lead_source=t.lead_source,
            full_or=fo, full_p=pub.P.get(t.variant_id, eur(e, cohort)[1]), lead_a1=lead_a1,
            sub_wgs_or_lead=ld.OR if ld is not None else np.nan, sub_wgs_p_lead=ld.P if ld is not None else np.nan,
            proxy=best, proxy_r2=float(r ** 2), proxy_dist_kb=round(abs(int(pr.POS_w) - t.pos) / 1000, 1),
            sub_wgs_or_proxy=o(pr.OR_w), sub_wgs_p_proxy=pr.P_w, sub_arr_or_proxy=o(pr.or_a), sub_arr_p_proxy=pr.P_a,
            disc_cases=pr.disc_rate_cases, disc_controls=pr.disc_rate_controls,
            miss_wgs_cases=pr.miss_wgs_cases, miss_wgs_controls=pr.miss_wgs_controls,
            miss_arr_cases=pr.miss_arr_cases, miss_arr_controls=pr.miss_arr_controls,
            region_sites_r2ge02=int(len(near)),
            region_r_z=float(np.corrcoef(near.Z_STAT_w, near.z_a)[0, 1]) if len(near) >= 3 else np.nan,
            region_max_abs_dz=float((near.z_a - near.Z_STAT_w).abs().max()) if len(near) else np.nan))
    SG = pd.DataFrame(rows)
    out.mkdir(parents=True, exist_ok=True)
    SG.to_csv(out / 'signals.tsv', sep='\t', index=False)
    direct = T[T.variant_id.isin(Q.index)].merge(
        Q[['OR_w', 'or_a', 'P_w', 'P_a', 'disc_rate_cases', 'disc_rate_controls']], left_on='variant_id', right_index=True)
    direct.to_csv(out / 'signals.direct_on_array.tsv', sep='\t', index=False)
    S['signals'] = {'leads': int(len(SG)), 'with_proxy_r2_ge_0.8': int((SG.proxy_r2 >= 0.8).sum()),
                    'with_proxy_r2_ge_0.5': int((SG.proxy_r2 >= 0.5).sum()),
                    'credible_set_members_on_array': int(direct.variant_id.nunique())}

    # DQB1:57:V: dosage predicted from array SNPs (ridge, trained on the cohort's controls; out of fold for them)
    res = pd.read_csv(K.PROJECT / 'hla.typing.v2' / 'results' / '06.assoc_inputs' / cohort / 'residue_dosage.tsv',
                      sep='\t', usecols=['sample_id', V57], dtype={'sample_id': str}).set_index('sample_id')[V57]
    MX, _ = raw(I / 'arr.mhc.raw')
    MX = MX.loc[MX.index.isin(cases | controls)]
    feat = [c for c in MX.columns if c in Q.index]
    MX = MX[feat]
    y = res.reindex(MX.index)
    ok = y.notna()
    MX, y = MX[ok], y[ok]
    is_case = MX.index.isin(sorted(cases))
    Xc, Xk = MX[~is_case], MX[is_case]
    mu, sd = Xc.mean(), Xc.std().replace(0, 1)
    z = lambda X: ((X.fillna(mu) - mu) / sd).values  # noqa: E731
    alphas = np.logspace(-1, 4, 30)
    oof = np.zeros((~is_case).sum())
    for tr, te in KFold(10, shuffle=True, random_state=1).split(Xc):
        m = RidgeCV(alphas=alphas).fit(z(Xc.iloc[tr]), y[~is_case].iloc[tr])
        oof[te] = m.predict(z(Xc.iloc[te]))
    model = RidgeCV(alphas=alphas).fit(z(Xc), y[~is_case])
    pred_case = model.predict(z(Xk))
    P = pd.Series(np.r_[oof, pred_case], index=list(Xc.index) + list(Xk.index)).clip(0, 2).reindex(MX.index)
    yc, yk, pc, pk = y[~is_case], y[is_case], P[~is_case], P[is_case]
    resid_c, resid_k = yc - pc, yk - pk
    tt = stats.ttest_ind(resid_k, resid_c, equal_var=False)
    H = {'features': len(feat), 'alpha': float(model.alpha_),
         'n_cases': int(is_case.sum()), 'n_controls': int((~is_case).sum()),
         'r2_controls_out_of_fold': float(np.corrcoef(yc, pc)[0, 1] ** 2),
         'r2_cases_trained_on_controls': float(np.corrcoef(yk, pk)[0, 1] ** 2),
         'freq_typed_cases': float(yk.mean() / 2), 'freq_typed_controls': float(yc.mean() / 2),
         'freq_array_pred_cases': float(pk.mean() / 2), 'freq_array_pred_controls': float(pc.mean() / 2),
         'resid_mean_cases': float(resid_k.mean()), 'resid_mean_controls': float(resid_c.mean()),
         'resid_diff_welch_t': float(tt.statistic), 'resid_diff_p': float(tt.pvalue),
         'best_call_agreement_cases': float((np.rint(pk) == yk).mean()),
         'best_call_agreement_controls': float((np.rint(pc) == yc).mean())}
    cov = pd.read_csv(cov_path, sep='\t', dtype={'IID': str}).set_index('IID')
    base = cov.reindex(MX.index)[COVARS]
    base = (base - base.mean()) / base.std()
    yy = is_case.astype(int)

    def logit(dose):
        X = sm.add_constant(pd.concat([pd.Series(dose.values, index=MX.index, name='d'), base], axis=1))
        f = sm.Logit(yy, X).fit(disp=0)
        b, se = f.params['d'], f.bse['d']
        return {'or': float(np.exp(b)), 'l95': float(np.exp(b - 1.96 * se)), 'u95': float(np.exp(b + 1.96 * se)),
                'p': float(f.pvalues['d'])}
    H['assoc_typed'] = logit(y)
    H['assoc_array_predicted'] = logit(P)
    S['hla_dqb1_57V'] = H
    (out / 'summary.json').write_text(json.dumps(S, indent=2) + '\n')
    return SG


# ── step 05: rare drivers typed directly on one array ───────────────────────────────────────────────────────────
def rare_direct(xc, cohort, w, out):
    O = xc / '05.rare_direct'
    keep = set(pd.read_csv(w / 'keep.txt', sep='\t', header=None, dtype=str)[1])
    fixed = set(pd.read_csv(K.MODEL_INPUTS / cohort / 'genotype' / 'fixed_model' / f'{FIXED}.fam', sep=r'\s+',
                            header=None, dtype=str)[1])
    ph = pd.read_csv(K.MODEL_INPUTS / cohort / 'phenotype' / 'pheno.tsv', sep='\t', dtype={'IID': str}).set_index('IID').PHENO1

    def rd(path):
        r = pd.read_csv(path, sep=r'\s+').drop(columns=['FID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'])
        r['IID'] = r.IID.astype(str)
        return r.set_index('IID')
    W = rd(O / 'wgs.raw')
    W = W[W.index.isin(fixed)]                 # the cohort's unrelated set
    rows = []
    for arr in ('illuminaASA24v1', 'illuminaJSA24v1'):
        p = O / f'{arr}.raw'
        if not p.exists():
            continue
        R = rd(p)
        common = R.index.intersection(W.index).intersection(ph.index)
        for ca in R.columns:
            aid, a_counted = ca.rsplit('_', 1)
            pos = ':'.join(aid.split(':')[:2])
            ref, alt = aid.split(':')[2:4]
            for cw in [c for c in W.columns if c.startswith(pos + ':')]:
                wid, w_counted = cw.rsplit('_', 1)
                if set(wid.split(':')[2:4]) != {ref, alt}:
                    continue
                a = R.loc[common, ca] if a_counted == alt else 2 - R.loc[common, ca]
                g_ = W.loc[common, cw] if w_counted == alt else 2 - W.loc[common, cw]
                for gv, lab in ((1, 'cases'), (0, 'controls')):
                    idx = ph.loc[common][ph.loc[common] == gv].index
                    tab = pd.crosstab(g_.loc[idx].fillna(-1).astype(int), a.loc[idx].fillna(-1).astype(int))
                    rows.append(dict(variant=DIRECT_LABEL.get(pos, pos), array_id=aid, array=arr[8:], group=lab,
                                     n=len(idx), wgs_carriers=int((g_.loc[idx] >= 1).sum()),
                                     array_carriers=int((a.loc[idx] >= 1).sum()),
                                     both_carriers=int(((g_.loc[idx] >= 1) & (a.loc[idx] >= 1)).sum()),
                                     wgs_only=int(((g_.loc[idx] >= 1) & (a.loc[idx] == 0)).sum()),
                                     array_only=int(((g_.loc[idx] == 0) & (a.loc[idx] >= 1)).sum()),
                                     missing_wgs=int(g_.loc[idx].isna().sum()), missing_array=int(a.loc[idx].isna().sum()),
                                     table='; '.join(f'WGS {i}/array {j}: {tab.loc[i, j]}' for i in tab.index
                                                     for j in tab.columns if tab.loc[i, j])))
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / 'rare_direct.tsv', sep='\t', index=False)
    return keep


# ── step 06: the imputed arrays ──────────────────────────────────────────────────────────────────────────────────
def imputed(T, xc, cohort, w, pheno, cov, SG, out):
    O = xc / '06.imputed'
    T.plink2('--pfile', IMPUTED, '--no-psam-pheno', '--keep', w / 'keep.txt', '--extract', O / 'leads.ids',
             '--pheno', pheno, *GLM[:3], '--covar', cov, *GLM[3:], '--out', w / 'imp.leads',
             log_ok=(w / 'imp.leads.log', r'11 covariates loaded'))
    ph = pd.read_csv(pheno, sep='\t', dtype={'IID': str}).set_index('IID').PHENO1
    keep = set(pd.read_csv(w / 'keep.txt', sep='\t', header=None, dtype=str)[1])
    W = alt_counts(O / 'wgs.raw')
    D = alt_counts(O / 'imp.raw')
    ids = W.index.intersection(D.index).intersection(ph.index)
    ids = ids[ids.isin(keep)]
    W, D, ph = W.loc[ids], D.loc[ids], ph.loc[ids]
    case, ctrl = (ph == 1).values, (ph == 0).values

    def info(path):
        t = pd.read_csv(path, sep='\t', header=None, names=['id', 'info']).set_index('id')['info']
        get = lambda s, k: next((float(x.split('=')[1]) for x in s.split(';') if x.startswith(k + '=')), np.nan)  # noqa: E731
        return pd.DataFrame({'r2': t.map(lambda s: get(s, 'R2')), 'typed': t.map(lambda s: 'TYPED' in s.split(';'))})
    inf = pd.concat([info(O / 'imp.drivers_window.info'), info(O / 'imp.leads.info')])
    inf = inf[~inf.index.duplicated()]
    shared = [v for v in W.columns if v in D.columns]

    rare = [v for v in (O / 'rare_shared.ids').read_text().split() if v in shared and not inf.typed.get(v, True)]
    rows = []
    for v in rare:
        g, d = W[v].values, D[v].values
        if np.nanmean(g) / 2 > 0.5:
            g, d = 2 - g, 2 - d
        cc = int(np.nansum(g[ctrl] >= 1))
        for grp, m in (('cases', case), ('controls', ctrl)):
            car, non = m & (g >= 1), m & (g == 0)
            rows.append(dict(variant=v, group=grp, control_carriers_wgs=cc, bin=binlab(cc), wgs_carriers=int(car.sum()),
                             recovered=int((car & (d >= 0.5)).sum()), wgs_noncarriers=int(non.sum()),
                             imputed_only=int((non & (d >= 0.5)).sum())))
    R = pd.DataFrame(rows)
    rec = (R.groupby(['group', 'bin'], sort=False)
           .agg(variants=('variant', lambda s: int((R.loc[s.index, 'wgs_carriers'] > 0).sum())),
                wgs_carriers=('wgs_carriers', 'sum'), recovered=('recovered', 'sum'),
                wgs_noncarriers=('wgs_noncarriers', 'sum'), imputed_only=('imputed_only', 'sum')).reset_index())
    rec['recovery'] = rec.recovered / rec.wgs_carriers.where(rec.wgs_carriers > 0)
    rec['false_carrier_rate'] = rec.imputed_only / rec.wgs_noncarriers
    rec['bin'] = pd.Categorical(rec.bin, [b[2] for b in BINS], ordered=True)
    rec = rec.sort_values(['group', 'bin'])
    out.mkdir(parents=True, exist_ok=True)
    rec.to_csv(out / 'recovery_by_group.tsv', sep='\t', index=False)
    expected = rec[rec.group == 'cases'].astype({'bin': str}).set_index('bin').recovery

    drv = []
    for v, lab in DRIVERS.items():
        if v not in W.columns:
            continue
        g = W[v]
        d = D[v] if v in D.columns else None
        cc = int((g[ctrl] >= 1).sum())
        for grp, m in (('cases', case), ('controls', ctrl)):
            car = m & (g.values >= 1)
            row = dict(variant=lab, id=v, in_imputed_set=d is not None, typed=bool(inf.typed.get(v, False)),
                       r2=inf.r2.get(v, np.nan), group=grp, n=int(m.sum()), wgs_carriers=int(car.sum()),
                       wgs_missing=int(np.isnan(g.values[m]).sum()))
            if d is not None:
                dv = d.values
                row.update(imputed_carriers=int((m & (dv >= 0.5)).sum()), both=int((car & (dv >= 0.5)).sum()),
                           wgs_only=int((car & (dv < 0.5)).sum()),
                           imputed_only=int((m & (g.values == 0) & (dv >= 0.5)).sum()),
                           mean_dosage_wgs_carriers=round(float(dv[car].mean()), 3) if car.any() else np.nan,
                           imputed_dosage_sum=round(float(dv[m].sum()), 2))
            if grp == 'cases':
                row['expected_case_recovery'] = expected.get(binlab(cc), np.nan)
            drv.append(row)
    pd.DataFrame(drv).to_csv(out / 'drivers.tsv', sep='\t', index=False)

    sig = pd.read_csv(SG, sep='\t')                  # as array_crosscheck: the written table, read back
    sig = sig[sig.lead_source == 'wgs']
    gl = pd.read_csv(w / 'imp.leads.PHENO1.glm.logistic.hybrid', sep='\t').set_index('ID')
    lead_rows = []
    for r in sig.itertuples():
        v = r.lead
        if v not in gl.index:
            lead_rows.append(dict(signal=r.signal, lead=v, rsid=r.rsid, in_imputed_set=False))
            continue
        e = gl.loc[v]
        a1 = str(r.lead_a1)
        or_imp = e.OR if e.A1 == a1 else 1 / e.OR
        alt = v.split(':')[3]
        g, d = W[v], D[v]
        sgn = 1 if a1 == alt else -1
        fq = lambda x, m: float((x[m] if sgn == 1 else 2 - x[m]).mean() / 2)  # noqa: E731,B023
        r2 = lambda m: float(np.corrcoef(g.values[m][~np.isnan(g.values[m])],  # noqa: E731,B023
                                         d.values[m][~np.isnan(g.values[m])])[0, 1] ** 2)
        lead_rows.append(dict(signal=r.signal, lead=v, rsid=r.rsid, a1=a1, in_imputed_set=True,
                              typed=bool(inf.typed.get(v, False)), info_r2=inf.r2.get(v, np.nan),
                              or_wgs=r.sub_wgs_or_lead, p_wgs=r.sub_wgs_p_lead, or_imputed=or_imp, p_imputed=e.P,
                              freq_cases_wgs=fq(g, case), freq_cases_imputed=fq(d, case),
                              freq_controls_wgs=fq(g, ctrl), freq_controls_imputed=fq(d, ctrl),
                              r2_cases=r2(case), r2_controls=r2(ctrl), proxy_r2_on_array=r.proxy_r2))
    pd.DataFrame(lead_rows).to_csv(out / 'leads.tsv', sep='\t', index=False)

    # genotype concordance with WGS (user 2026-11-05: one Minimac R2 per variant, and the concordance in cases and in
    # controls): the best-guess imputed genotype (dosage < 0.5 -> 0, < 1.5 -> 1, else 2) against the WGS genotype, in
    # the people with both; also among those with a non-reference genotype on either side. A table of its own, not
    # one of gate G10's (array_crosscheck has none)
    con = []
    for kind, v, lab in ([('lead', r['lead'], r['signal']) for r in lead_rows if r.get('in_imputed_set')] +
                         [('driver', v, lab) for v, lab in DRIVERS.items() if v in W.columns and v in D.columns]):
        g, d = W[v].values.astype(float), D[v].values.astype(float)
        bg = np.where(d < 0.5, 0, np.where(d < 1.5, 1, 2))
        for grp, m in (('cases', case), ('controls', ctrl)):
            k = m & ~np.isnan(g) & ~np.isnan(d)
            nr = k & ((g >= 1) | (bg >= 1))
            con.append(dict(kind=kind, label=lab, variant=v, group=grp, n=int(k.sum()),
                            concordant=int((bg[k] == g[k]).sum()),
                            concordance=float((bg[k] == g[k]).mean()) if k.any() else np.nan,
                            nonref=int(nr.sum()),
                            nonref_concordance=float((bg[nr] == g[nr]).mean()) if nr.any() else np.nan))
    pd.DataFrame(con).to_csv(out / 'concordance.tsv', sep='\t', index=False)


# ── gate G10: full mainland reproduces array_crosscheck ─────────────────────────────────────────────────────────
def same_table(a, b):
    """(identical?, largest relative difference, what differs) of two tables read as published"""
    x, y = pd.read_csv(a, sep='\t', dtype=str, keep_default_na=False), pd.read_csv(b, sep='\t', dtype=str, keep_default_na=False)
    if list(x.columns) != list(y.columns) or len(x) != len(y):
        return False, np.inf, f'shape {x.shape} vs {y.shape}'
    worst, what = 0.0, ''
    for c in x.columns:
        xn, yn = pd.to_numeric(x[c], errors='coerce'), pd.to_numeric(y[c], errors='coerce')
        num = xn.notna() & yn.notna()
        if (x[c][~num] != y[c][~num]).any():
            return False, np.inf, f'column {c} differs as text'
        if num.any():
            d = (np.abs(xn[num] - yn[num]) / np.maximum(np.abs(yn[num]), 1e-12)).max()
            if d > worst:
                worst, what = float(d), c
    return worst <= 1e-6, worst, what


def same_json(a, b, path=''):
    """the largest relative difference over the numbers of b also present in a"""
    worst = 0.0
    if isinstance(b, dict):
        for k, v in b.items():
            if k not in a:
                return np.inf
            worst = max(worst, same_json(a[k], v, f'{path}.{k}'))
    elif isinstance(b, (int, float)) and not isinstance(b, bool):
        worst = abs(a - b) / max(abs(b), 1e-12) if np.isfinite(b) else (0.0 if a == b or (a != a and b != b) else np.inf)
    elif a != b:
        return np.inf
    return worst


def gate(xc, full_out):
    rec = {'gate': 'G10', 'what': "the array re-run in full mainland reproduces check/array_crosscheck's tables",
           'rule': 'every table: same rows, columns and text; numbers within 1e-6 relative', 'tables': {}}
    pairs = [('04.compare/signals.tsv',) * 2, ('04.compare/signals.direct_on_array.tsv',) * 2,
             ('05.rare_direct/rare_direct.tsv',) * 2, ('06.imputed/leads.tsv',) * 2, ('06.imputed/drivers.tsv',) * 2,
             ('06.imputed/recovery_by_group.tsv',) * 2]
    ok = True
    for mine, pub in pairs:
        same, d, what = same_table(full_out / mine, xc / pub)
        rec['tables'][pub] = {'pass': bool(same), 'max_rel_diff': d, 'column': what}
        ok &= same
    js = json.loads((xc / '04.compare' / 'summary.json').read_text())
    js.get('hla_dqb1_57V', {}).pop('tags', None)
    d = same_json(json.loads((full_out / '04.compare' / 'summary.json').read_text()), js)
    rec['tables']['04.compare/summary.json'] = {'pass': bool(d <= 1e-6), 'max_rel_diff': d}
    ok &= d <= 1e-6
    rec['pass'] = bool(ok)
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--array', type=Path, default=K.PROJECT / 'check' / 'array_crosscheck' / 'results',
                    help='check/array_crosscheck/results (read only)')
    ap.add_argument('--plink2', default='/home/b/b37974/plink2')
    ap.add_argument('--plink1', default='/home/b/b37974/plink')
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--memory', type=int, default=3500, help='MB for plink (below the task\'s --rsc memory)')
    ap.add_argument('--work', type=Path, default=Path('array_work'))
    ap.add_argument('--out', type=Path, default=Path('cohorts'))
    a = ap.parse_args()
    T = Tools(a.plink2, a.plink1, a.threads, a.memory)
    xc = a.array
    gate_rec = None
    for cohort in [K.COHORTS[-1]] + K.COHORTS[:-1]:            # full mainland first: gate G10 before the others
        w, o = a.work / cohort, a.out / cohort
        pheno, cov = prepare(T, xc, cohort, w)
        compare(xc, cohort, w, cov, o / '04.compare')
        rare_direct(xc, cohort, w, o / '05.rare_direct')
        imputed(T, xc, cohort, w, pheno, cov, o / '04.compare' / 'signals.tsv', o / '06.imputed')
        s = json.loads((o / '04.compare' / 'summary.json').read_text())['samples']
        print(f'[array_cohorts] {cohort}: {s["cases"]} cases, {s["controls"]} controls')
        if cohort == 'full_mainland':
            gate_rec = gate(xc, o)
            a.out.mkdir(parents=True, exist_ok=True)
            (a.out / 'array.gate.json').write_text(json.dumps(gate_rec, indent=2, default=float) + '\n')
            print(f'[array_cohorts] gate G10: {"pass" if gate_rec["pass"] else "FAIL"}: ' +
                  ', '.join(f'{k} {v["max_rel_diff"]:.1e}' for k, v in gate_rec['tables'].items()))
            if not gate_rec['pass']:
                sys.exit('ABORT: gate G10 failed (array re-run differs from check/array_crosscheck)')


if __name__ == '__main__':
    main()
