#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The data checks verify.sh runs, recomputed from the inputs rather
#           than read back from the pipeline's own summaries wherever possible.
#           Read-only: nothing under results/ is written.
#
#   §1  input validation passed for every cohort (validation.json)
#   §12 case sets: recomputed from the sample sheet's DVT status and the
#       phenotype file (all controls + the set's cases; for dvt_case_only the
#       cases of known status only, PHENO1 recoded 1 = DVT yes / 0 = DVT no)
#   §2  PLINK round-trip: every .bed decodes to the typing matrix, cell by cell
#   §3  engine logs: samples, cases, controls, covariates, --1, LOCO, impute
#   §4  every tested marker has exactly one result per scan
#   §5  marker_info counts recomputed from the matrices
#   §6  tests recomputed from the matrix (polymorphic markers, mirror pairs);
#       GEC v0.2 re-run on an input rebuilt here reproduces Me; the three
#       thresholds follow from tests / Me / max-T; the GEC threshold's
#       family-wise error under max-T is close to 0.05; every sig_* flag and
#       count recomputed from P
#   §7  EVERY fixed-model marker refitted exactly (statsmodels Logit, Newton,
#       tol 1e-10) on the same samples and covariates. plink2's logistic fit
#       stops short of the MLE for a few markers (METHODS §3): the INFO line
#       counts markers off by > 1e-3 in log10 P. PASS requires that the ORs
#       agree (1e-3 relative) and that no call at any threshold would change
#       if the exact P replaced plink2's
#   §8  SAIGE rows: BETA = Tstat/var, p.value.NA from Tstat/sqrt(var), OR = exp(BETA)
#   §9  every null model (reused for 'all', fitted for the DVT sets) passed
#       CHECK_NULL, and its samples equal the case set's random-model samples
#   §10 _comparison tables equal the per-cohort tables stacked, per case set;
#       the n/3 consistency of every marker and threshold recomputed from the
#       random-model rows; _case_sets/scan_qc_case_sets.tsv holds all 48 scans
#   §11 trace.txt: one completed task per expected process instance
#   §14 benchmark (MAXT_BOOT): replicates, samples, cases and tests as the
#       scan's; the observed phenotype's Wald rerun equals ASSOC_FIXED's P on
#       every test; its score P recomputed for 10 tests by weighted least
#       squares; every exact threshold and FWER recomputed from the replicates
#   §15 omnibus tests of residue positions (no gene-level omnibus output may exist): for 5 units per scan the categories are
#       rebuilt from the matrix by the documented rule and the likelihood-ratio
#       test refitted (statsmodels); the Bonferroni threshold and every call
#       recomputed; the stepwise rule (candidates, lead, members) recomputed
#       step by step; _omnibus/omnibus_all.tsv = the per-case-set tables
#   §16 HLA | signal and positions | marker (HLA_GIVEN_SIGNAL, OMNIBUS_GIVEN): their marginal values equal
#       CONDITIONAL's and OMNIBUS's; the lead given each set equals lead_given_hla.tsv; residues given a
#       residue signal lead equal STEPWISE's step 1 for that lead
#   §13 conditional analyses: the pooled PLINK set decodes to the matrices and
#       the lead's random_model genotypes, SENTINEL first; the stepwise step-0
#       P of every class marker equals the class scan's; the joint selection
#       rule (candidates, lead, entry, members) recomputed step by step; every
#       signal lead tested in a cohort is among that cohort's lead | HLA
#       conditioning sets; every lead | whole-position set holds the position's
#       residues but its most frequent, for the residue signal leads' and the
#       omnibus leads' positions; and, for case set 'all', SAIGE with the genome-wide
#       scan's default imputation (best_guess) on the pooled set reproduces
#       the scan's own P for each lead (same genotypes, same null model)
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import warnings
from scipy import stats

CASE_SETS = ('all', 'dvt_yes', 'dvt_no', 'dvt_case_only')
COHORTS = ('narrow_mainland', 'intermediate_mainland', 'full_mainland')
CLASSES = ('allele', 'residue')
COVARS = ['SEX'] + [f'PC{i}_AVG' for i in range(1, 11)]
THRESHOLDS = ('gec', 'maxt', 'bonf')
BOOT_DRAWS = 10000                               # params.MaxTBootDraws
MAX_SIGNALS = 5                                  # params.CondMaxSignals
EXPECTED_TASKS = {'GENE_COORDS': 1, 'CASE_SETS': 1, 'NULL_INFO_ALL': 3, 'VALIDATE_INPUTS': 3,
                  'PREP': 12, 'FIT_NULL': 9, 'NULL_INFO_FIT': 9, 'CHECK_NULL': 12,
                  'BUILD_MARKERS': 24, 'EFFECTIVE_TESTS': 24, 'MAXT_BOOT': 24, 'ASSOC_FIXED': 24, 'ASSOC_RANDOM': 24,
                  'COLLECT': 48,
                  'COMPARE': 4, 'FIGURES': 8, 'CASE_SET_COMPARE': 1,
                  'GWAS_LEADS': 1, 'EXTRACT_LEAD': 3, 'BUILD_POOLED': 12, 'STEPWISE': 16,
                  'CONDITIONAL': 24, 'HLA_GIVEN_SIGNAL': 24, 'OMNIBUS_GIVEN': 4,
                  'CONDITIONAL_COMPARE': 1, 'OMNIBUS': 4, 'OMNIBUS_SUMMARY': 1,
                  'WRITE_RUN_MANIFEST': 1}

N_PASS = N_FAIL = 0


def note(msg):
    print(f'  INFO  {msg}')


def report(ok, msg):
    global N_PASS, N_FAIL
    if ok:
        N_PASS += 1
        print(f'  PASS  {msg}')
    else:
        N_FAIL += 1
        print(f'  FAIL  {msg}')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results', required=True)
    p.add_argument('--typing', required=True, help='hla.typing.v2/results/06.assoc_inputs')
    p.add_argument('--model-inputs', required=True)
    p.add_argument('--refit-per-scan', type=int, default=10)
    p.add_argument('--plink2', default='/home/b/b37974/plink2')
    p.add_argument('--saige-rscript', default='/home/b/b37974/anaconda3/envs/saige/bin/Rscript')
    p.add_argument('--saige-step2', default='/home/b/b37974/anaconda3/envs/saige/bin/step2_SPAtests.R')
    p.add_argument('--gwas-dir', required=True, help='analysis/assoc_saige/results')
    p.add_argument('--tmp', required=True, help='scratch directory for the re-runs')
    p.add_argument('--gec-jar', required=True)
    p.add_argument('--java', default='/usr/bin/java')
    return p.parse_args()


def read_bed(prefix):
    """PLINK1 SNP-major .bed -> (samples x variants) count of the bim column-5 allele, NaN = missing."""
    fam = pd.read_csv(f'{prefix}.fam', sep=r'\s+', header=None, dtype=str)
    bim = pd.read_csv(f'{prefix}.bim', sep='\t', header=None, dtype=str)
    n, m = len(fam), len(bim)
    raw = np.fromfile(f'{prefix}.bed', dtype=np.uint8)
    if raw[:3].tolist() != [0x6C, 0x1B, 0x01]:
        raise SystemExit(f'{prefix}.bed is not SNP-major PLINK1')
    bpv = (n + 3) // 4
    codes = np.unpackbits(raw[3:].reshape(m, bpv)[:, :, None], axis=2, bitorder='little')
    codes = codes.reshape(m, bpv * 4, 2)[:, :n, :]
    two = codes[:, :, 0] + 2 * codes[:, :, 1]        # 0 hom A1, 1 missing, 2 het, 3 hom A2
    lut = np.array([2.0, np.nan, 1.0, 0.0])
    return fam[1].tolist(), bim, lut[two].T


def recompute_tests(names, dose, cls):
    """(polymorphic, mirror partner, test_id) per marker, from the matrix rows of the case set."""
    called = dose.notna()
    copies = dose.sum()
    mac = np.minimum(copies, 2 * called.sum() - copies)
    poly = (mac > 0).to_dict()
    mirror = {n: '' for n in names}
    tid = {n: (n if poly[n] else '') for n in names}
    if cls == 'residue':
        by = {}
        for n in names:
            if poly[n]:
                g, pos, _ = n.split(':')
                by.setdefault((g, pos), []).append(n)
        for ms in by.values():
            if len(ms) == 2:
                x, y = dose[ms[0]], dose[ms[1]]
                if x.isna().equals(y.isna()) and np.allclose((x + y).dropna(), 2.0):
                    mirror[ms[0]], mirror[ms[1]] = ms[1], ms[0]
                    tid[ms[1]] = ms[0]
    return poly, mirror, tid


def gec_rerun(a, prefix, keep, reps, out):
    """GEC v0.2 on an input rebuilt here; (Me, markers)."""
    import subprocess
    Path(f'{out}.extract').write_text('\n'.join(reps) + '\n')
    subprocess.run([a.plink2, '--bfile', prefix, '--keep', keep, '--extract', f'{out}.extract', '--mac', '1',
                    '--make-bed', '--out', out], check=True, stdout=subprocess.DEVNULL)
    bim = pd.read_csv(f'{out}.bim', sep='\t', header=None, dtype=str)
    bim[4] = bim[4].replace('P', 'T')
    bim[5] = bim[5].replace('P', 'T')
    bim.to_csv(f'{out}.bim', sep='\t', header=False, index=False)
    r = subprocess.run([a.java, '-Xmx4g', '-jar', a.gec_jar, '--effect-number', '--plink-binary', out, '--genome',
                        '--maf', '0', '--no-web', '--out', out], check=True, capture_output=True, text=True)
    m = re.search(r'total estimated effective number of indpendent markers is ([0-9.Ee+-]+) for the (\d+) markers',
                  r.stdout + r.stderr)
    return (float(m.group(1)), int(m.group(2))) if m else (np.nan, 0)


def majority(sig, orr):
    up = sum(1 for x, o in zip(sig, orr) if x and o > 1)
    down = sum(1 for x, o in zip(sig, orr) if x and o < 1)
    return max(up, down)


def main():
    a = parse_args()
    R = Path(a.results)
    T = Path(a.typing)
    summ = json.load(open(T / 'summary.json'))['cohorts']

    print('§1 input validation')
    for c in COHORTS:
        v = json.load(open(R / 'all' / c / '00.prep' / 'validation.json'))
        report(v['n_failed'] == 0 and v['n_checks'] >= 30, f"{c}: {v['n_checks']} checks, {v['n_failed']} failed")

    dvt = pd.read_csv(R / '_reference' / 'case_sets.tsv', sep='\t', dtype=str).set_index('IID')['dvt']
    n_scans = 0
    for cs in CASE_SETS:
        all_assoc, all_qc = [], []
        for c in COHORTS:
            n_scans += check_case_set(a, R / cs, cs, c, T, summ, dvt, all_assoc, all_qc)

        print(f'§10 {cs}: _comparison tables')
        cmp = pd.read_csv(R / cs / '_comparison' / 'assoc_all.tsv', sep='\t', dtype={'marker': str})
        stacked = pd.concat(all_assoc, ignore_index=True)
        key = ['cohort', 'class', 'model', 'marker']
        j = cmp.merge(stacked, on=key, suffixes=('', '_s'))
        report(len(cmp) == len(stacked) == len(j) and np.allclose(j.p, j.p_s, rtol=1e-5, equal_nan=True)
               and set(cmp.case_set) == {cs},
               f'assoc_all.tsv = 12 per-cohort tables ({len(cmp)} rows)')
        qa = pd.read_csv(R / cs / '_comparison' / 'scan_qc_all.tsv', sep='\t')
        report(len(qa) == 12, f'scan_qc_all.tsv has 12 scans ({len(qa)})')
        ms = pd.read_csv(R / cs / '_comparison' / 'marker_summary.tsv', sep='\t', dtype={'marker': str})
        rnd = stacked[stacked.model == 'random'].set_index(['class', 'marker', 'cohort'])
        bad = 0
        for r in ms.itertuples():
            for t in THRESHOLDS:
                sig, orr = [], []
                for c in COHORTS:
                    k = (r[ms.columns.get_loc('class') + 1], r.marker, c)
                    sig.append(bool(rnd.at[k, f'sig_{t}']) if k in rnd.index else False)
                    orr.append(rnd.at[k, 'or'] if k in rnd.index else np.nan)
                n = majority(sig, orr)
                bad += int(getattr(r, f'{t}_n') != n or getattr(r, t) != f'{n}/3')
        hits = pd.read_csv(R / cs / '_comparison' / 'hits.tsv', sep='\t', dtype={'marker': str})
        want_hits = ms[(ms[[f'{t}_n' for t in THRESHOLDS]].sum(axis=1) > 0) | (ms.gec_n_fixed > 0)]
        report(bad == 0 and set(hits.marker + hits['class']) == set(want_hits.marker + want_hits['class']),
               f'n/3 consistency recomputed for {len(ms)} markers x 3 thresholds (random model, same direction); '
               f'hits.tsv = {len(hits)} markers significant somewhere')

    print('§10 _case_sets')
    q36 = pd.read_csv(R / '_case_sets' / 'scan_qc_case_sets.tsv', sep='\t')
    n_all = 12 * len(CASE_SETS)
    report(len(q36) == n_all and len(q36.drop_duplicates(['case_set', 'cohort', 'class', 'model'])) == n_all,
           f'scan_qc_case_sets.tsv has {n_all} distinct scans ({len(q36)})')

    check_conditional(a, R, T, dvt)
    check_omnibus(R, T)
    check_given(R)

    print('§11 trace.txt')
    tr = pd.read_csv(R / '_run_info' / 'trace.txt', sep='\t')
    tr['proc'] = tr['name'].str.extract(r'^([A-Z_]+)')
    good = tr[tr.status.isin(['COMPLETED', 'CACHED'])].groupby('proc').size().to_dict()
    for p, n in EXPECTED_TASKS.items():
        report(good.get(p, 0) == n, f'{p}: {good.get(p, 0)} of {n}')

    print(f'\n{N_PASS} passed, {N_FAIL} failed')
    sys.exit(1 if N_FAIL else 0)


def check_conditional(a, R, T, dvt):
    import subprocess
    tmp = Path(a.tmp)
    tmp.mkdir(parents=True, exist_ok=True)
    gl = pd.read_csv(R / '_conditional' / 'gwas_leads.tsv', sep='\t')
    leads = gl.loc[gl.status == 'kept', 'lead_id'].tolist()
    report(len(leads) >= 1, f"§13 GWAS MHC leads kept: {', '.join(leads)}; excluded: "
           f"{', '.join(gl.loc[gl.status != 'kept', 'lead_id']) or 'none'}")
    for c in COHORTS:
        # the lead's genotypes straight from the random_model set, once per cohort
        geno = sorted((Path(a.model_inputs) / c / 'genotype' / 'random_model').glob('*.random_model.bed'))[0]
        (tmp / 'ids').write_text('\n'.join(leads) + '\n')
        subprocess.run([a.plink2, '--bfile', str(geno)[:-4], '--extract', str(tmp / 'ids'), '--make-bed',
                        '--out', str(tmp / f'lead_{c}')], check=True, stdout=subprocess.DEVNULL)
        lids, lbim, lg = read_bed(str(tmp / f'lead_{c}'))
        for cs in CASE_SETS:
            base = R / cs / c / '04.conditional'
            pr = pd.read_csv(R / cs / c / '00.prep' / 'pheno_covar.random.tsv', sep='\t', dtype={'IID': str})
            print(f'§13 {cs} {c}: pooled set')
            iids, bim, g = read_bed(str(base / 'pooled' / 'pooled'))
            report(bim.iloc[0, 1] == 'SENTINEL' and iids == pr.IID.tolist(),
                   f'SENTINEL is record 1; samples = case set ({len(iids)})')
            ok = True
            for cls in CLASSES:
                mat = pd.read_csv(T / c / f'{cls}_dosage.tsv', sep='\t', dtype=str, keep_default_na=False)
                mat = mat.set_index('sample_id').loc[pr.IID]
                dose = mat.replace('', np.nan).astype(float)
                gd = pd.DataFrame(g, columns=bim[1].tolist())[list(mat.columns)]
                ok &= np.array_equal(np.nan_to_num(gd.values, nan=-1), np.nan_to_num(dose.values, nan=-1))
            # lead: pooled counts the ALT (column 5 of pooled.bim); the random_model .bed counts its column 5
            lsrc = pd.DataFrame(lg, index=lids, columns=lbim[1].tolist()).loc[pr.IID]
            for ld in leads:
                col = bim[1].tolist().index(ld)
                same_a1 = bim.iloc[col, 4] == lbim.set_index(1).at[ld, 4]
                pooled_ld = g[:, col] if same_a1 else 2 - g[:, col]
                ok &= np.array_equal(np.nan_to_num(pooled_ld, nan=-1), np.nan_to_num(lsrc[ld].values, nan=-1))
            report(ok, 'HLA columns = the matrices; lead = the random_model genotypes, cell by cell')

            for model in ('fixed', 'random'):
                cj = json.load(open(base / model / 'conditional.json'))
                lgv = pd.read_csv(base / model / 'lead_given_hla.tsv', sep='\t',
                                  dtype={'conditioned_on': str, 'position': str, 'position_reference': str})
                tested = set(bim[1])
                want_sl = [m for m in sum((open(R / cs / '_conditional' / model / f'leads.{k}.txt').read().split()
                                           for k in CLASSES), []) if m in tested]
                got = set(lgv.loc[lgv.why.str.contains('signal_lead'), 'conditioned_on'])
                report(sorted(cj['signal_leads']) == sorted(want_sl) and set(want_sl) <= got,
                       f'§13 {cs} {c} {model}: lead | HLA includes the {len(want_sl)} signal lead(s) tested here')
                # whole positions: every residue column of the position but its most frequent in the pooled set,
                # for the residue signal leads' positions and the omnibus position leads
                names = bim[1].tolist()
                om = pd.read_csv(R / cs / '_omnibus' / 'omnibus_signals.residue.tsv', sep='\t', dtype={'lead': str})
                want_pos = {m.rsplit(':', 1)[0] for m in want_sl if '*' not in m} | set(om.lead.dropna())
                ok_pos, n_pos = True, 0
                for u in sorted(want_pos):
                    cols = [m for m in names if '*' not in m and m != 'SENTINEL' and m not in leads
                            and m.rsplit(':', 1)[0] == u]
                    if not cols:
                        continue
                    fr = {m: np.nanmean(g[:, names.index(m)]) for m in cols}
                    ref = max(cols, key=fr.get) if len(cols) > 1 else ''
                    exp = sorted(m for m in cols if m != ref)
                    rows = lgv[lgv.why.str.startswith('position_jointly') & (lgv.position.astype(str) == u)]
                    ok_pos &= len(rows) == len(leads) and all(sorted(str(x).split(';')) == exp for x in rows.conditioned_on)
                    ok_pos &= all(str(x) == (ref or 'nan') for x in rows.position_reference)
                    n_pos += 1
                report(ok_pos and n_pos == len(lgv[lgv.why.str.startswith('position_jointly')]) // max(len(leads), 1),
                       f'§13 {cs} {c} {model}: lead | whole position for {n_pos} position(s) '
                       f'({", ".join(sorted(want_pos)) or "none"}): every residue but the most frequent')

            if cs == 'all':
                nm = a.gwas_dir + f'/{c}/01.assoc/main/null'
                out = tmp / f'bg_{c}.txt'
                subprocess.run([a.saige_rscript, a.saige_step2, f'--bedFile={base}/pooled/pooled.bed',
                                f'--bimFile={base}/pooled/pooled.bim', f'--famFile={base}/pooled/pooled.fam',
                                '--AlleleOrder=alt-first', '--chrom=6', f'--GMMATmodelFile={nm}.rda',
                                f'--varianceRatioFile={nm}.varianceRatio.txt', f'--SAIGEOutputFile={out}',
                                '--is_Firth_beta=FALSE', '--LOCO=TRUE', '--is_output_moreDetails=TRUE'],
                               check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                bg = pd.read_csv(out, sep=r'\s+').set_index('MarkerID')
                for ld in leads:
                    gp = float(gl.loc[gl.lead_id == ld, f'p_{c}'].iloc[0])
                    report(abs(bg.at[ld, 'p.value.NA'] / gp - 1) < 1e-4,
                           f'§13 {c} {ld}: pooled set + best_guess reproduces the genome-wide scan '
                           f"(P {bg.at[ld, 'p.value.NA']:.4e} vs {gp:.4e})")

    for cs in CASE_SETS:
        for model in ('fixed', 'random'):
            for cls in CLASSES:
                check_stepwise(R, cs, model, cls)


def check_stepwise(R, cs, model, cls):
    """§13: the joint selection rule recomputed from stepwise.<class>.tsv."""
    d = R / cs / '_conditional' / model
    st = pd.read_csv(d / f'steps.{cls}.tsv', sep='\t', dtype={'lead': str})
    sw = pd.read_csv(d / f'stepwise.{cls}.tsv', sep='\t', dtype={'marker': str})
    sg = pd.read_csv(d / f'signals.{cls}.tsv', sep='\t', dtype={'lead': str})
    thr = {c: pd.read_csv(R / cs / c / '01.markers' / cls / f'effective_tests.{cls}.tsv', sep='\t').iloc[0]
           for c in COHORTS}
    # step 0 = the class scans
    s0 = sw[sw.step == 0].set_index(['cohort', 'marker'])
    ref = pd.concat([pd.read_csv(R / cs / c / '03.tables' / f'assoc.{model}.{cls}.tsv', sep='\t',
                                 dtype={'marker': str}).assign(cohort=c) for c in COHORTS]).set_index(['cohort', 'marker'])
    j = s0[['p']].join(ref[['p']], rsuffix='_scan', how='outer')
    if model == 'random':
        ok0 = len(j) == len(ref) and np.allclose(j['p'], j['p_scan'], rtol=1e-4, equal_nan=True)
    else:
        dd = np.abs(np.log10(j['p']) - np.log10(j['p_scan']))
        note(f'§13 {cs} {model} {cls}: pooled vs class-scan plink2 P, {int((dd > 1e-3).sum())} of {len(j)} '
             'off by > 1e-3 in log10 P')
        ok0 = len(j) == len(ref)
    # every step: candidates, lead, entry, members
    ok, cond, prev, members = True, [], set(), {}
    for k, row in st.sort_values('step').iterrows():
        e = sw[sw.step == row.step]
        P = e.pivot(index='marker', columns='cohort', values='p').reindex(columns=list(COHORTS))
        sig = pd.DataFrame({c: P[c].lt(float(thr[c].p_gec)) for c in COHORTS})
        sig_any = set(sig.index[sig.any(axis=1)])
        cand = [m for m in sig_any if m not in cond]
        maxp = P.fillna(1.0).max(axis=1)
        entered = bool(cand) and row.step < MAX_SIGNALS
        ok &= bool(row.entered) == entered
        if cand:
            ok &= row.lead in cand and np.isclose(maxp[row.lead], maxp[cand].min(), rtol=1e-9)
        if row.step > 0 and cond:
            members[cond[-1]] = sorted(prev - sig_any)
        if entered:
            cond.append(row.lead)
        prev = sig_any
    okm = list(sg.lead) == cond and all(
        sorted(str(r.members).split(';')) == sorted(members.get(r.lead, []))
        if r.n_members else not members.get(r.lead) for r in sg.itertuples())
    report(ok0 and ok and okm,
           f'§13 {cs} {model} {cls}: step 0 = class scans; selection rule recomputed over {len(st)} step(s); '
           f"signals {', '.join(f'{r.lead} ({r.gec})' for r in sg.itertuples()) or 'none'}")


def check_maxt_boot(a, R, c, cls, tag, th, coding, names, dose, mat):
    """§14 for one scan: the bootstrap benchmark of MAXT_BOOT."""
    print(f'§14 {tag}: benchmark (exact max-T by parametric bootstrap)')
    md = R / c / '01.markers' / cls
    bt = pd.read_csv(md / f'maxt_boot.{cls}.tsv', sep='\t').iloc[0]
    mp = pd.read_csv(md / f'maxt_boot.{cls}.minp.tsv.gz', sep='\t')
    ob = pd.read_csv(md / f'maxt_boot.{cls}.observed.tsv', sep='\t', dtype={'marker': str}).set_index('marker')
    fx = pd.read_csv(R / c / '00.prep' / 'pheno_covar.fixed.tsv', sep='\t', dtype={'IID': str})
    rel = lambda x, y: abs(x / y - 1) < 1e-4                                       # noqa: E731
    report(int(bt.draws) == BOOT_DRAWS == len(mp) and int(bt.samples) == len(fx)
           and int(bt.cases) == coding['fixed']['cases'] and int(bt.controls) == coding['fixed']['controls']
           and int(bt.tests) == int(th.tests) == len(ob) and rel(bt.p_gec, th.p_gec) and rel(bt.p_bonf, th.p_bonf),
           f"{int(bt.draws)} replicates on {int(bt.samples)} fixed-model samples ({int(bt.cases)} cases), "
           f"{int(bt.tests)} tests, thresholds as the scan's")
    fixed = pd.read_csv(md.parent.parent / '02.assoc' / 'fixed' / cls / f'fixed.{cls}.PHENO1.glm.logistic', sep='\t',
                        dtype={'ID': str}, na_values=['NA']).set_index('ID').P.reindex(ob.index)
    same = (ob.p_wald.isna() == fixed.isna()).all() and np.allclose(ob.p_wald.dropna(), fixed.dropna(), rtol=1e-4)
    report(same, f'observed phenotype rerun: Wald P = ASSOC_FIXED P on all {len(ob)} tests')
    y = fx['PHENO1'].astype(float).values
    X = sm.add_constant(fx[COVARS].astype(float).values)
    mu = sm.GLM(y, X, family=sm.families.Binomial()).fit().fittedvalues
    w = mu * (1 - mu)
    idx = mat.set_index('sample_id').index.get_indexer(fx['IID'])
    worst = 0.0
    for mk in ob.index[:: max(1, len(ob) // 10)][:10]:
        g = dose[mk].values[idx]
        g = np.where(np.isnan(g), np.nanmean(g), g)
        gt = g - sm.WLS(g, X, weights=w).fit().fittedvalues
        z = np.sum(gt * (y - mu)) / np.sqrt(np.sum(w * gt * gt))
        worst = max(worst, abs(np.log10(2 * stats.norm.sf(abs(z))) - np.log10(ob.at[mk, 'p_score'])))
    report(worst < 1e-3, f'observed phenotype: score P of 10 tests recomputed by weighted least squares '
                         f'(max |delta log10 P| {worst:.1e})')
    ok = True
    for s in ('wald', 'score'):
        m = mp[f'min_p_{s}'].values
        ok &= rel(bt[f'p_maxt_{s}'], np.quantile(m, 0.05)) and bt[f'p_maxt_{s}_lo'] <= bt[f'p_maxt_{s}'] <= bt[f'p_maxt_{s}_hi']
        for t in ('gec', 'bonf'):
            ok &= abs(bt[f'fwer_{t}_{s}'] - (m < th[f'p_{t}']).mean()) < 1e-6
    report(ok, f"thresholds and FWER recomputed from the replicates: FWER at GEC Wald {bt.fwer_gec_wald:.3f}, "
               f"score {bt.fwer_gec_score:.3f}; at Bonferroni Wald {bt.fwer_bonf_wald:.3f}, score {bt.fwer_bonf_score:.3f}")


def omnibus_categories(xu):
    """The documented category rule, re-implemented: dict label -> copies (0 off the called rows), reference."""
    called = ~xu.isna().any(axis=1).values
    src = xu.fillna(0.0)
    cat = {c: np.where(called, src[c].values, 0.0) for c in src.columns}
    unl = np.where(called, 2.0 - src.sum(axis=1).values, 0.0)
    if unl.max() > 1e-9:
        cat['unlisted'] = np.clip(unl, 0, 2)
    n2 = 2.0 * called.sum()
    fr = {k: v.sum() / n2 for k, v in cat.items()}
    own = {k: v for k, v in cat.items() if fr[k] >= 0.01}
    pool = [k for k in cat if 0 < fr[k] < 0.01]
    if pool and sum(fr[k] for k in pool) >= 0.005:
        own['rare'] = sum(cat[k] for k in pool)
    if len(own) < 2:
        return None, None, called
    ref = max(own, key=lambda k: own[k].sum())
    return {k: v for k, v in own.items() if k != ref}, ref, called


def close(x, y, rtol=1e-5):
    """Equal as written: both missing, or relative difference below rtol."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    both = np.isnan(x) & np.isnan(y)
    return bool(np.all(both | (np.abs(x - y) <= rtol * np.maximum(np.abs(x), np.abs(y)))))


def check_given(R):
    """§16: HLA_GIVEN_SIGNAL and OMNIBUS_GIVEN against the tables they extend."""
    for cs in CASE_SETS:
        for c in COHORTS:
            for model in ('fixed', 'random'):
                base = R / cs / c / '04.conditional' / model
                tag = f'§16 {cs} {c} {model}'
                if not (base / 'hla_given_signal.tsv').exists():
                    report(False, f'{tag}: hla_given_signal.tsv missing')
                    continue
                g = pd.read_csv(base / 'hla_given_signal.tsv', sep='\t', dtype={'marker': str})
                if g.empty:
                    note(f'{tag}: no signal lead, nothing to condition on')
                    continue
                hl = pd.read_csv(base / 'hla_given_lead.tsv', sep='\t', dtype={'marker': str})
                m0 = hl.drop_duplicates('marker').set_index('marker').p_marginal
                gm = g[g['class'] != 'gwas_lead'].drop_duplicates('marker').set_index('marker').p_marginal
                ok_m = close(gm.reindex(m0.index), m0)
                lg = pd.read_csv(base / 'lead_given_hla.tsv', sep='\t')
                lead_ok, n_sets = True, 0
                for (cond, ld), d in g[g['class'] == 'gwas_lead'].groupby(['conditioned_on', 'marker']):
                    ref = lg[(lg.conditioned_on == cond) & (lg.lead == ld)].lead_p_given
                    if len(ref):
                        n_sets += 1
                        lead_ok &= close(d.p_given.iloc[0], ref.iloc[0])
                st_ok, n_st = True, 0
                stw = R / cs / '_conditional' / model / 'stepwise.residue.tsv'
                if stw.exists():
                    W = pd.read_csv(stw, sep='\t', dtype={'marker': str})
                    for lead_ in g.loc[g.why == 'signal_lead', 'conditioned_on'].unique():
                        w = W[(W.cohort == c) & (W.step == 1) & (W.errcode == '.')]
                        steps = pd.read_csv(R / cs / '_conditional' / model / 'steps.residue.tsv', sep='\t')
                        if not ((steps.step == 1) & (steps.conditioned_on == lead_)).any():
                            continue
                        x = g[(g.conditioned_on == lead_) & (g['class'] == 'residue')].set_index('marker').p_given
                        w = w.set_index('marker').p
                        n_st += 1
                        st_ok &= close(x.reindex(w.index), w)
                report(ok_m and lead_ok and st_ok,
                       f'{tag}: marginal = CONDITIONAL ({len(m0)} markers); lead given {n_sets} set(s) = '
                       f'lead_given_hla; residues given {n_st} residue lead(s) = STEPWISE step 1')
        f = R / cs / '_omnibus' / 'omnibus_given.residue.tsv'
        if not f.exists():
            report(False, f'§16 {cs}: omnibus_given.residue.tsv missing')
            continue
        og = pd.read_csv(f, sep='\t', dtype={'unit': str})
        om = pd.read_csv(R / cs / '_omnibus' / 'omnibus.residue.tsv', sep='\t', dtype={'unit': str})
        j = og.merge(om[['cohort', 'unit', 'p']].rename(columns={'p': 'p_omnibus'}), on=['cohort', 'unit'], how='left')
        # every unit OMNIBUS tested is there, for every conditioning marker (OMNIBUS also lists, as NOT_TESTED, the
        # units of other cohorts; the cohort's own set is what is retested)
        sets = og.groupby(['cohort', 'conditioned_on']).unit.apply(frozenset)
        tested = om[om.p.notna()].groupby('cohort').unit.apply(frozenset)
        complete = all(tested[c_] <= u for (c_, _m), u in sets.items()) and \
            all(len(set(sets.xs(c_, level='cohort'))) == 1 for c_ in og.cohort.unique())
        report(close(j.p_marginal, j.p_omnibus) and complete,
               f"§16 {cs}: omnibus given {', '.join(og.conditioned_on.unique())}: marginal = OMNIBUS for every unit; "
               "every tested unit present for every marker")


def check_omnibus(R, T):
    """§15: the omnibus tests of residue positions (a main analysis; fixed model, Bonferroni)."""
    stacked = []
    stale = [str(f.relative_to(R)) for cs in CASE_SETS for f in (R / cs / '_omnibus').glob('omnibus*.allele.tsv')]
    report(not stale, 'no gene-level (allele) omnibus output' + (f': {stale}' if stale else ''))
    for cs in CASE_SETS:
        for cls in ('residue',):
            tag = f'{cs} {cls}'
            print(f'§15 {tag}: omnibus tests')
            od = R / cs / '_omnibus'
            om = pd.read_csv(od / f'omnibus.{cls}.tsv', sep='\t', dtype={'unit': str, 'position': str})
            lg = pd.read_csv(od / f'omnibus_stepwise.{cls}.tsv', sep='\t', dtype={'unit': str})
            st = pd.read_csv(od / f'omnibus_steps.{cls}.tsv', sep='\t', dtype={'lead': str})
            sg = pd.read_csv(od / f'omnibus_signals.{cls}.tsv', sep='\t', dtype={'lead': str, 'members': str})
            stacked.append(om)
            ok_thr, worst, n_refit, df_ok = True, 0.0, 0, True
            for c in COHORTS:
                o = om[om.cohort == c].set_index('unit')
                tested = o.errcode != 'NOT_TESTED'
                thr = 0.05 / int(tested.sum())
                ok_thr &= np.allclose(o.threshold, thr, rtol=1e-4)
                want = o.p.fillna(1.0).lt(thr) & o.errcode.eq('.')
                ok_thr &= bool((o.sig_bonf.astype(bool) == want).all())
                fx = pd.read_csv(R / cs / c / '00.prep' / 'pheno_covar.fixed.tsv', sep='\t', dtype={'IID': str})
                mat = pd.read_csv(T / c / f'{cls}_dosage.tsv', sep='\t', dtype=str, keep_default_na=False)
                mat = mat.set_index('sample_id').loc[fx.IID]
                y = fx.PHENO1.astype(float).values
                Cv = sm.add_constant(fx[COVARS].astype(float).values)
                pick = o.index[tested & o.p.notna() & o.errcode.eq('.')]
                for u in pick[np.linspace(0, len(pick) - 1, min(5, len(pick))).astype(int)]:
                    cols = [m for m in mat.columns if (m.rsplit(':', 1)[0] if cls == 'residue' else m.split('*')[0]) == u]
                    xu = mat[cols].replace('', np.nan).astype(float)
                    xu.columns = [m.rsplit(':', 1)[1] if cls == 'residue' else m for m in cols]
                    cats, ref, called = omnibus_categories(xu)
                    Xa = np.column_stack(list(cats.values()))[called]
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore')
                        f0 = sm.Logit(y[called], Cv[called]).fit(disp=0, method='newton', maxiter=200, tol=1e-10)
                        f1 = sm.Logit(y[called], np.column_stack([Cv[called], Xa])).fit(disp=0, method='newton',
                                                                                         maxiter=200, tol=1e-10)
                    chisq = 2 * (f1.llf - f0.llf)
                    worst = max(worst, abs(chisq - o.at[u, 'chisq']) / max(1.0, chisq))
                    df_ok &= int(o.at[u, 'df']) == Xa.shape[1] and o.at[u, 'reference'] == ref
                    n_refit += 1
            report(ok_thr, 'Bonferroni threshold = 0.05 / units tested, per cohort; every call recomputed '
                           '(P < threshold and a clean fit)')
            report(worst < 1e-4 and df_ok, f'{n_refit} units refitted from the matrix (categories, reference, df; '
                                           f'max relative chi-square difference {worst:.1e})')
            ok, entered = True, []
            for k in range(len(st)):
                cur = lg[lg.step == k]
                if cur.empty:
                    break
                sig = cur[cur.sig_bonf.astype(bool)]
                cand = [u for u in dict.fromkeys(sig.unit) if u not in entered]
                row = st[st.step == k].iloc[0]
                if not cand:
                    ok &= not bool(row.entered)
                    break
                pv = cur.assign(pp=cur.p.fillna(1.0)).groupby('unit')
                mx, mn = pv.pp.max(), cur.groupby('unit').p.min()
                lead = min(cand, key=lambda u: (mx[u], mn[u]))
                ok &= bool(row.entered) and row.lead == lead
                nxt = lg[lg.step == k + 1]
                before = set(sig.unit)
                after = set(nxt[nxt.sig_bonf.astype(bool)].unit) - set(entered + [lead])
                members = [u for u in dict.fromkeys(cur.unit) if u in before and u not in after]
                srow = sg[sg.lead == lead]
                ok &= len(srow) == 1 and set(str(srow.members.iloc[0]).split(';')) == set(members)
                entered.append(lead)
            report(ok and len(sg) == len(entered), f'stepwise rule recomputed: {len(entered)} signal(s) '
                                                   f'({", ".join(entered) or "none"})')
    print('§15 _omnibus summary')
    oa = pd.read_csv(R / '_omnibus' / 'omnibus_all.tsv', sep='\t', dtype={'unit': str, 'position': str})
    st_all = pd.concat(stacked, ignore_index=True)
    report(len(oa) == len(st_all) and set(oa.case_set) == set(CASE_SETS) and set(oa['class']) == {'residue'},
           f'omnibus_all.tsv = the 4 per-case-set tables, residues only ({len(oa)} rows)')


def check_case_set(a, R, cs, c, T, summ, dvt, all_assoc, all_qc):
    """§2-§9 and §12 for one (case set, cohort); returns the number of scans checked."""
    n_scans = 0
    prep = R / c / '00.prep'
    coding = json.load(open(prep / 'pheno_coding.json'))
    pr = pd.read_csv(prep / 'pheno_covar.random.tsv', sep='\t', dtype={'IID': str, '#FID': str})
    keep = set(pd.read_csv(prep / 'fixed.keep', sep='\t', header=None, dtype=str)[1])
    nm = json.load(open(prep / 'null_model.json'))
    ns = pd.read_csv(prep / 'null_model.samples.tsv', sep='\t', dtype=str)
    pheno = pd.read_csv(Path(a.model_inputs) / c / 'phenotype' / 'pheno.tsv', sep='\t', dtype=str)

    print(f'§12 {cs} {c}: case set')
    status = pheno['IID'].map(dvt).fillna('unknown')
    is_case = pheno['PHENO1'] == '1'
    if cs == 'all':
        want = pheno
    elif cs == 'dvt_case_only':
        want = pheno[is_case & status.isin(['yes', 'no'])].copy()
        want['PHENO1'] = (status[want.index] == 'yes').astype(int).astype(str)
    else:
        want = pheno[~is_case | (status == cs.split('_')[1])]
    report(pr['IID'].tolist() == want['IID'].tolist()
           and pr['PHENO1'].astype(str).tolist() == want['PHENO1'].tolist(),
           f'random-model samples and PHENO1 = {cs} recomputed from the sheet, phenotype-file order ({len(pr)})')
    report(int((pr.PHENO1 == 1).sum()) == coding['random']['cases'],
           f"cases {coding['random']['cases']}, controls {coding['random']['controls']}")

    print(f'§9 {cs} {c}: null model')
    nc = json.load(open(prep / 'null_check.json'))
    report(nc['n_failed'] == 0, f"CHECK_NULL: {nc['n_checks']} checks, {nc['n_failed']} failed")
    report(set(ns.sample_id) == set(pr.IID) and len(ns) == len(pr),
           f'null model samples = random-model samples ({len(ns)})')

    for cls in CLASSES:
        tag = f'{cs} {c} {cls}'
        mat = pd.read_csv(T / c / f'{cls}_dosage.tsv', sep='\t', dtype=str, keep_default_na=False)
        mat = mat[mat.sample_id.isin(set(pr.IID))].reset_index(drop=True)
        names = list(mat.columns[1:])
        dose = mat[names].replace('', np.nan).astype(float)

        print(f'§2 {tag}: PLINK round-trip')
        iids, bim, g = read_bed(str(R / c / '01.markers' / cls / cls))
        report(iids == mat['sample_id'].tolist(), f'.fam samples = matrix rows ({len(iids)})')
        report(set(bim[4]) == {'P'} and set(bim[5]) == {'A'}, 'bim alleles: column 5 P, column 6 A')
        gdf = pd.DataFrame(g, columns=bim[1].tolist())[names]
        same = np.array_equal(np.isnan(gdf.values), np.isnan(dose.values)) and \
            np.array_equal(np.nan_to_num(gdf.values, nan=-1), np.nan_to_num(dose.values, nan=-1))
        report(same and len(bim) == len(names) == int(summ[c][f'{cls}s_tested']),
               f'{len(names)} markers x {len(iids)} samples equal cell by cell')

        print(f'§5 {tag}: marker_info recomputed')
        info = pd.read_csv(R / c / '01.markers' / cls / 'marker_info.tsv', sep='\t', dtype={'marker': str})
        y = pr['PHENO1'].astype(int).values
        ok_all = True
        for sset, mask in (('random', np.ones(len(pr), bool)), ('fixed', pr['IID'].isin(keep).values)):
            inf = info[info.sample_set == sset].set_index('marker').loc[names]
            for grp, gy in (('case', 1), ('control', 0)):
                sel = mask & (y == gy)
                called = dose[sel].notna().sum().values
                copies = dose[sel].sum().values
                ok_all &= np.array_equal(inf[f'called_{grp}'].values, called)
                ok_all &= np.array_equal(inf[f'copies_{grp}'].values, copies.astype(int))
                ok_all &= bool((inf[f'n_{grp}'] == sel.sum()).all())
        report(ok_all, 'called / copies / n per group, both sample sets')

        print(f'§6 {tag}: tests and thresholds')
        md = R / c / '01.markers' / cls
        th = pd.read_csv(md / f'effective_tests.{cls}.tsv', sep='\t').iloc[0]
        tt = pd.read_csv(md / f'tests.{cls}.tsv', sep='\t', dtype=str, keep_default_na=False).set_index('marker')
        order = bim[1].tolist()                     # file order: a pair's test_id is its first member there
        poly, mirror, tid = recompute_tests(order, dose[order], cls)
        reps = [n for n in order if tid[n] == n]
        report(list(tt.index) == order and (tt.polymorphic.astype(int).astype(bool) == pd.Series(poly)[order]).all()
               and (tt.mirror_of == pd.Series(mirror)[order]).all() and (tt.test_id == pd.Series(tid)[order]).all()
               and int(th.tests) == len(reps),
               f"tests recomputed: {len(reps)} ({int(th.mirror_pairs)} mirror pairs, {int(th.monomorphic)} monomorphic)")
        me, n_in = gec_rerun(a, str(md / cls), str(prep / 'fixed.keep'), reps, str(Path(a.tmp) / f'gec_{cs}_{c}_{cls}'))
        glog = open(md / f'gec.{cls}.log').read()
        report('v0.2' in glog and abs(me - th.me_gec_raw) < 1e-3 and n_in == int(th.gec_markers)
               and abs(th.me_gec - (th.me_gec_raw + th.gec_added)) < 1e-6 and n_in + int(th.gec_added) == len(reps),
               f'GEC v0.2 re-run reproduces Me {me:.2f} on {n_in} tests ({int(th.gec_added)} added as independent)')
        rel = lambda x, y: abs(x / y - 1) < 1e-4                                   # noqa: E731
        report(rel(th.p_bonf, 0.05 / len(reps)) and rel(th.p_gec, 0.05 / th.me_gec) and 1 <= th.me_gec <= len(reps)
               and th.p_maxt_lo <= th.p_maxt <= th.p_maxt_hi and th.p_bonf <= th.p_maxt and rel(th.me_maxt, 0.05 / th.p_maxt),
               f'thresholds: GEC {th.p_gec:.3g} (Me {th.me_gec:.1f}), max-T {th.p_maxt:.3g} (Me {th.me_maxt:.1f}), '
               f'Bonferroni {th.p_bonf:.3g}')
        report(0.035 <= th.fwer_gec <= 0.065,
               f'family-wise error of the GEC threshold under max-T: {th.fwer_gec:.3f} (Bonferroni {th.fwer_bonf:.3f})')

        check_maxt_boot(a, R, c, cls, tag, th, coding, names, dose, mat)

        for model in ('fixed', 'random'):
            t = f'{tag} {model}'
            assoc = pd.read_csv(R / c / '03.tables' / f'assoc.{model}.{cls}.tsv', sep='\t', dtype={'marker': str})
            qc = pd.read_csv(R / c / '03.tables' / f'scan_qc.{model}.{cls}.tsv', sep='\t')
            all_assoc.append(assoc)
            all_qc.append(qc)
            n_scans += 1
            ad = R / c / '02.assoc' / model / cls

            print(f'§3 {t}: engine log')
            want = coding[model]
            if model == 'fixed':
                log = open(ad / f'fixed.{cls}.log').read()
                mm = re.search(r'(\d+) cases and (\d+) controls remaining after main filters', log)
                report(mm and int(mm.group(1)) == want['cases'] and int(mm.group(2)) == want['controls'],
                       f"plink2 analysed {mm.group(1) if mm else '?'} cases / {mm.group(2) if mm else '?'} controls "
                       f"= {want['cases']} / {want['controls']}")
                report(bool(re.search(r'\b11 covariates loaded', log)), 'plink2 loaded 11 covariates')
                report(bool(re.search(r'^\s+--1\s*$', log, re.M)), 'plink2 ran with --1 (0/1 phenotype)')
                report(bool(re.search(r'--glm .*no-firth', log)), 'plink2 ran no-firth')
            else:
                log = open(ad / f'random.{cls}.log').read()
                raw = pd.read_csv(ad / f'random.{cls}.saige.txt', sep=r'\s+', dtype={'MarkerID': str})
                report(set(raw.N_case) == {want['cases']} and set(raw.N_ctrl) == {want['controls']},
                       f"SAIGE N_case / N_ctrl = {want['cases']} / {want['controls']} on every row")
                report('[1] "mean"' in log and bool(re.search(r'Leave chromosome\s+6\s+out', log)),
                       'SAIGE: impute_method mean, LOCO chromosome 6')
                vr = re.search(r'variance Ratio null is\s+([0-9.]+)', log)
                report(vr and abs(float(vr.group(1)) - nm['variance_ratio'][0]) < 1e-6,
                       f"SAIGE variance ratio {vr.group(1) if vr else '?'} = null model's")

                print(f'§8 {t}: SAIGE rows')
                beta_ok = np.allclose(raw.BETA, raw.Tstat / raw['var'], rtol=1e-4, atol=1e-8)
                pna = 2 * stats.norm.sf(np.abs(raw.Tstat) / np.sqrt(raw['var']))
                p_ok = np.allclose(-np.log10(raw['p.value.NA']), -np.log10(pna), rtol=1e-3, atol=1e-4)
                j = assoc.set_index('marker').loc[raw.MarkerID]
                or_ok = np.allclose(j['or'].values, np.exp(raw.BETA.values), rtol=1e-4)
                pp_ok = np.allclose(j['p'].values, raw['p.value.NA'].values, rtol=1e-4)
                report(beta_ok and p_ok and or_ok and pp_ok,
                       'BETA = Tstat/var; p.value.NA = 2·Φ(−|Tstat|/√var); table or = exp(BETA), p = p.value.NA')

            print(f'§4 {t}: one result per tested marker')
            report(len(assoc) == len(names) and assoc.marker.is_unique and set(assoc.marker) == set(names),
                   f'{len(assoc)} rows, {len(names)} tested markers')

            print(f'§6 {t}: calls')
            q0 = qc.iloc[0]
            pf = assoc['p'].fillna(1.0).values
            ok6 = all((assoc[f'sig_{k}'].values == (pf < th[f'p_{k}'])).all()
                      and int(q0[f'n_{k}']) == int(assoc[f'sig_{k}'].sum())
                      and int(q0[f'tests_{k}']) == assoc.loc[assoc[f'sig_{k}'], 'test_id'].nunique()
                      and rel(q0[f'p_{k}'], th[f'p_{k}']) for k in THRESHOLDS)
            report(ok6 and (assoc.set_index('marker').test_id == tt.test_id).all(),
                   'sig_gec / sig_maxt / sig_bonf and counts recomputed (' +
                   ', '.join(f"{k} {int(assoc[f'sig_{k}'].sum())}" for k in THRESHOLDS) + ')')

            if model == 'fixed':
                print(f'§7 {t}: every marker refitted exactly')
                fx = pd.read_csv(R / c / '00.prep' / 'pheno_covar.fixed.tsv', sep='\t', dtype={'IID': str})
                idx = mat.set_index('sample_id').index.get_indexer(fx['IID'])
                yv = fx['PHENO1'].astype(int).values
                Cv = fx[COVARS].astype(float).values
                ai = assoc.set_index('marker')
                ex_p = pd.Series(np.nan, index=ai.index)
                ex_or = pd.Series(np.nan, index=ai.index)
                for mk in ai.index[ai.p.notna()]:
                    x = dose[mk].values[idx]
                    ok = ~np.isnan(x)
                    X = sm.add_constant(np.column_stack([x[ok], Cv[ok]]))
                    try:
                        with warnings.catch_warnings():
                            warnings.simplefilter('ignore')
                            fit = sm.Logit(yv[ok], X).fit(disp=0, method='newton', maxiter=200, tol=1e-10)
                        ex_p[mk], ex_or[mk] = fit.pvalues[1], np.exp(fit.params[1])
                    except Exception:
                        pass
                dd = np.abs(np.log10(ex_p) - np.log10(ai.p))
                note(f'plink2 vs exact fit: {int((dd > 1e-3).sum())} of {int(ai.p.notna().sum())} marker(s) off by '
                     f'> 1e-3 in log10 P (max {np.nanmax(dd) if dd.notna().any() else 0:.2g}: '
                     f'{dd.idxmax() if dd.notna().any() else "-"})')
                p_ex = ex_p.where(ex_p.notna(), ai.p).fillna(1.0)
                n_change = int(sum(((p_ex < th[f'p_{k}']) != ai[f'sig_{k}']).sum() for k in THRESHOLDS))
                or_ok = np.nanmax(np.abs(ex_or / ai['or'] - 1)) < 1e-3 if ex_or.notna().any() else True
                report(n_change == 0 and or_ok and ex_p[ai.p.notna()].notna().all(),
                       f'{int(ai.p.notna().sum())} markers refitted; calls at all three thresholds unchanged with '
                       f'exact P ({n_change} would change); OR within 1e-3')
    return n_scans


if __name__ == '__main__':
    main()
