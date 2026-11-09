#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The significance thresholds of one (case set, cohort, class), fixed
#           from the genotypes alone, BEFORE any association P exists.
#
#           TESTS. A marker with no minor allele in the case set's samples is
#           not tested. Two residues at one IMGT position whose dosages are
#           exact complements in every sample (x2 = 2 - x1, blanks in the same
#           samples: a position with exactly two residues) are ONE test reported
#           twice, with reciprocal ORs and the same P; the first in file order
#           is the pair's test_id. Alleles have no such pairs.
#
#           THREE THRESHOLDS, each at family-wise alpha:
#             gec   PRIMARY. alpha / Me, Me = the effective number of independent
#                   tests from GEC v0.2 (Li et al. 2012), the original standalone
#                   release, archived in tools/gec-0.2 (sha256 checked). Input: the
#                   tests (one per mirror pair), the case set's fixed-model samples
#                   (unrelated, as GEC requires), --maf 0 (v0.2 drops MAF < 0.01 by
#                   default; every tested marker must count). A test monomorphic in
#                   those samples cannot enter GEC and is added to Me as one
#                   independent test.
#             maxt  VALIDATION. The exact family-wise threshold: the alpha
#                   quantile of min P under the joint null, simulated as max |Z|
#                   over the tests from a multivariate normal whose correlation is
#                   that of the covariate-adjusted dosages (sex + 10 PCs regressed
#                   out; blanks mean-imputed, as in the random model) over the
#                   case set's random-model samples. Also reported: the family-wise
#                   error rate the gec and bonf thresholds actually give.
#             bonf  SENSITIVITY. alpha / tests (conventional Bonferroni).
#
#           BLAS threads follow OMP_NUM_THREADS / OPENBLAS_NUM_THREADS, which the
#           calling process sets before Python starts.
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import hashlib
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import vocab as V                      # noqa: E402

CHUNK = 20000


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case-set', required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--class', dest='cls', required=True, choices=V.CLASSES)
    p.add_argument('--bfile', required=True, help="the class's PLINK prefix (BUILD_MARKERS)")
    p.add_argument('--pheno-random', required=True, help='pheno_covar.random.tsv (covariates)')
    p.add_argument('--keep', required=True, help='fixed.keep: unrelated samples for GEC')
    p.add_argument('--covar-cols', required=True)
    p.add_argument('--plink2', required=True)
    p.add_argument('--java', required=True)
    p.add_argument('--gec-jar', required=True)
    p.add_argument('--gec-sha256', required=True)
    p.add_argument('--alpha', type=float, default=0.05)
    p.add_argument('--draws', type=int, default=400000)
    p.add_argument('--seed', type=int, default=20261023)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--out-thresholds', required=True)
    p.add_argument('--out-tests', required=True)
    return p.parse_args()


def run(cmd, log=None):
    if log:
        with open(log, 'w') as fh:
            subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.STDOUT)
    else:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def read_dosage(a):
    """P dosage (samples x markers, NaN = blank) from the PLINK set, in file order."""
    run([a.plink2, '--bfile', a.bfile, '--export', 'A', '--threads', str(a.threads), '--out', 'dosage'])
    raw = pd.read_csv('dosage.raw', sep='\t', dtype={'FID': str, 'IID': str}, na_values=['NA'])
    bim = pd.read_csv(f'{a.bfile}.bim', sep='\t', header=None, dtype={1: str})
    cols = raw.columns[6:]
    names, x = [], []
    for c in cols:
        name, allele = c.rsplit('_', 1)
        v = raw[c].to_numpy(dtype=float)
        if allele == 'A':
            v = 2.0 - v
        elif allele != 'P':
            sys.exit(f'ERROR: plink2 counted allele {allele!r} for {name}; expected A or P')
        names.append(name)
        x.append(v)
    if names != bim[1].tolist():
        sys.exit('ERROR: exported markers differ from the .bim, or are in another order')
    return raw['IID'].tolist(), names, np.column_stack(x)


def tests_of(names, x, cls):
    """polymorphic flag, test_id and mirror partner per marker."""
    called = ~np.isnan(x)
    copies = np.nansum(x, axis=0)
    mac = np.minimum(copies, 2 * called.sum(axis=0) - copies)
    poly = mac > 0
    test_id = {n: (n if p else '') for n, p in zip(names, poly)}
    mirror = {n: '' for n in names}
    if cls == 'residue':
        idx = {}
        for i, n in enumerate(names):
            if poly[i]:
                g, pos, _ = V.parse_marker(n, cls)
                idx.setdefault((g, pos), []).append(i)
        for ii in idx.values():
            if len(ii) != 2:
                continue
            i, j = ii
            if (np.array_equal(called[:, i], called[:, j])
                    and np.allclose(x[called[:, i], i] + x[called[:, j], j], 2.0)):
                mirror[names[i]], mirror[names[j]] = names[j], names[i]
                test_id[names[j]] = names[i]
    return poly, test_id, mirror


def gec(a, reps):
    """Me from GEC v0.2 on the unrelated samples; (me_used, details)."""
    h = hashlib.sha256(Path(a.gec_jar).read_bytes()).hexdigest()
    if h != a.gec_sha256:
        sys.exit(f'ERROR: {a.gec_jar} sha256 {h} differs from the archived {a.gec_sha256}')
    Path('gec.extract').write_text('\n'.join(reps) + '\n')
    run([a.plink2, '--bfile', a.bfile, '--keep', a.keep, '--extract', 'gec.extract', '--mac', '1',
         '--make-bed', '--threads', str(a.threads), '--out', 'gec_in'])
    bim = pd.read_csv('gec_in.bim', sep='\t', header=None, dtype={1: str})
    if not set(bim[4]) | set(bim[5]) <= {'A', 'P'}:
        sys.exit('ERROR: GEC input alleles are not A/P')
    bim[4] = bim[4].replace('P', 'T')        # GEC reads nucleotide alleles only
    bim[5] = bim[5].replace('P', 'T')
    bim.to_csv('gec_in.bim', sep='\t', header=False, index=False)
    n_in = len(bim)
    n_samples = sum(1 for _ in open('gec_in.fam'))
    run([a.java, '-Xmx4g', '-jar', a.gec_jar, '--effect-number', '--plink-binary', 'gec_in', '--genome',
         '--maf', '0', '--no-web', '--out', 'gec'], log='gec.log')
    txt = Path('gec.log').read_text()
    m_read = re.search(r'number of SNPs on chromosome 6\s+in map file \S+ is (\d+)', txt)
    m_pass = re.search(r'(\d+) SNPs \(out of (\d+)\) pass the fliter', txt)    # printed only if a filter acts
    m_me = re.search(r'total estimated effective number of indpendent markers is ([0-9.Ee+-]+) for the (\d+) markers',
                     txt)
    if not m_read or not m_me:
        sys.exit('ERROR: GEC log does not report the markers read and the effective number (see gec.log)')
    used = {int(m_read.group(1)), int(m_me.group(2))} | ({int(m_pass.group(1))} if m_pass else set())
    if used != {n_in}:
        sys.exit(f'ERROR: GEC read {m_read.group(1)} and used {m_me.group(2)} of {n_in} input markers; '
                 'every test must enter')
    me_raw = float(m_me.group(1))
    blocks = pd.read_csv('gec.block.txt', sep='\t')
    added = len(reps) - n_in                  # tests monomorphic in the unrelated samples
    return me_raw + added, {'gec_samples': n_samples, 'gec_markers': n_in, 'gec_added': added,
                            'me_gec_raw': me_raw, 'gec_blocks': len(blocks), 'gec_sha256': h}


def max_t(a, iids, x, reps_idx, thr_gec, thr_bonf):
    """Exact family-wise threshold by simulation of max |Z| (see header)."""
    cov = pd.read_csv(a.pheno_random, sep='\t', dtype={'#FID': str, 'IID': str}).set_index('IID')
    cols = a.covar_cols.split(',')
    missing = [i for i in iids if i not in cov.index]
    if missing or len(cov) != len(iids):
        sys.exit(f'ERROR: covariates cover {len(cov)} samples, the PLINK set {len(iids)} ({len(missing)} absent)')
    xs = x[:, reps_idx].copy()
    mu = np.nanmean(xs, axis=0)
    xs = np.where(np.isnan(xs), mu, xs)
    z = np.column_stack([np.ones(len(iids)), cov.loc[iids, cols].astype(float).to_numpy()])
    resid = xs - z @ np.linalg.lstsq(z, xs, rcond=None)[0]
    r = np.corrcoef(resid.T)
    w, vec = np.linalg.eigh(r)
    load = vec * np.sqrt(np.clip(w, 0, None))
    rng = np.random.default_rng(a.seed)
    mx = np.empty(a.draws)
    for s in range(0, a.draws, CHUNK):
        k = min(CHUNK, a.draws - s)
        mx[s:s + k] = np.abs(rng.standard_normal((k, len(reps_idx))) @ load.T).max(axis=1)
    q = 1 - a.alpha
    thr = 2 * stats.norm.sf(np.quantile(mx, q))
    boot = [2 * stats.norm.sf(np.quantile(rng.choice(mx, a.draws), q)) for _ in range(50)]
    lo, hi = np.percentile(boot, [2.5, 97.5])
    fwer = lambda t: float((mx > stats.norm.isf(t / 2)).mean())       # noqa: E731
    return {'maxt_samples': len(iids), 'maxt_draws': a.draws, 'maxt_seed': a.seed,
            'p_maxt': thr, 'p_maxt_lo': lo, 'p_maxt_hi': hi, 'me_maxt': a.alpha / thr,
            'fwer_gec': fwer(thr_gec), 'fwer_bonf': fwer(thr_bonf)}


def main():
    a = parse_args()
    iids, names, x = read_dosage(a)
    poly, test_id, mirror = tests_of(names, x, a.cls)
    reps = [n for n in names if test_id[n] == n]
    reps_idx = [names.index(n) for n in reps]

    me, gd = gec(a, reps)
    gec_in = set(pd.read_csv('gec_in.bim', sep='\t', header=None, dtype={1: str})[1])
    thr_gec, thr_bonf = a.alpha / me, a.alpha / len(reps)
    mt = max_t(a, iids, x, reps_idx, thr_gec, thr_bonf)

    pd.DataFrame({'marker': names, 'polymorphic': poly.astype(int),
                  'test_id': [test_id[n] for n in names], 'mirror_of': [mirror[n] for n in names],
                  'in_gec': [int(n in gec_in) for n in names]}).to_csv(a.out_tests, sep='\t', index=False)
    row = {'case_set': a.case_set, 'cohort': a.cohort, 'class': a.cls, 'alpha': a.alpha,
           'markers': len(names), 'polymorphic': int(poly.sum()), 'monomorphic': int((~poly).sum()),
           'mirror_pairs': int(sum(1 for n in names if mirror[n] and test_id[n] != n)),
           'tests': len(reps), **{k: gd[k] for k in ('gec_samples', 'gec_markers', 'gec_added', 'me_gec_raw')},
           'me_gec': me, 'p_gec': thr_gec, 'gec_blocks': gd['gec_blocks'], **mt, 'p_bonf': thr_bonf,
           'gec_version': 'GEC v0.2 (standalone)', 'gec_sha256': gd['gec_sha256']}
    pd.DataFrame([row]).to_csv(a.out_thresholds, sep='\t', index=False, float_format='%.6g')
    print(f"[effective_tests] {a.case_set} {a.cohort} {a.cls}: {len(names)} markers, {len(reps)} tests "
          f"({row['mirror_pairs']} mirror pairs, {row['monomorphic']} monomorphic); GEC Me {me:.2f} -> "
          f"{thr_gec:.3g}; max-T {mt['p_maxt']:.3g} (Me {mt['me_maxt']:.1f}); Bonferroni {thr_bonf:.3g}; "
          f"FWER at GEC {mt['fwer_gec']:.3f}, at Bonferroni {mt['fwer_bonf']:.3f}")


if __name__ == '__main__':
    main()
