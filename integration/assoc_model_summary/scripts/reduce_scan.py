#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : One genome-wide scan (cohort x model) reduced to what the figures
#           need, and its calibration statistics:
#             - npz: chromosome, position and -log10 P of every usable variant
#               (Manhattan), and -log10 P sorted descending (QQ)
#             - stats.json: lambda_GC overall and by MAF bin, the number of
#               variants at P < 1e-5 and < 5e-8 against their expectation, min P,
#               sample sizes
#             - top.tsv: every variant at P < 1e-3 (for joins between models)
#           Usable variant = ERRCODE '.', 0 < P <= 1, finite OR > 0 and SE > 0 —
#           the rule of the upstream call_peaks.py / plot_manhattan_qq.py, so
#           lambda_GC and the counts must equal the upstream scan_qc.tsv (checked).
# Inputs  : PLINK2-style sumstats (#CHROM POS ID A1_FREQ OR LOG(OR)_SE P ERRCODE),
#           the upstream scan_qc.tsv
# Outputs : <cohort>.<model>.npz, <cohort>.<model>.stats.json, <cohort>.<model>.top.tsv
# ---------------------------------------------------------------------------
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import GW_ALPHA, MAF_EDGES, MAF_LABELS, SUGGEST_ALPHA, chrom_key, lambda_gc, log, maf_bin  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--sumstats', required=True)
    ap.add_argument('--scan-qc', required=True, help='upstream 02.scan/scan_qc.tsv of this cohort')
    ap.add_argument('--cohort', required=True)
    ap.add_argument('--model', required=True, help='this component\'s name: fixed / grm_pc / grm_only')
    ap.add_argument('--upstream-model', required=True, help='the model name in scan_qc.tsv')
    ap.add_argument('--top-p', type=float, default=1e-3)
    ap.add_argument('--out-prefix', required=True)
    a = ap.parse_args()

    cols = ['#CHROM', 'POS', 'ID', 'A1_FREQ', 'OR', 'LOG(OR)_SE', 'P', 'ERRCODE']
    head = pd.read_csv(a.sumstats, sep='\t', nrows=0).columns
    missing = [c for c in cols if c not in head]
    assert not missing, f'{a.sumstats}: missing columns {missing}'
    d = pd.read_csv(a.sumstats, sep='\t', usecols=cols,
                    dtype={'#CHROM': str, 'ID': str, 'ERRCODE': str, 'POS': np.int64})
    n_all = len(d)
    for c in ('A1_FREQ', 'OR', 'LOG(OR)_SE', 'P'):
        d[c] = pd.to_numeric(d[c], errors='coerce')
    ok = (d.ERRCODE.fillna('.') == '.') & d.P.notna() & (d.P > 0) & (d.P <= 1) \
        & np.isfinite(d.OR) & (d.OR > 0) & np.isfinite(d['LOG(OR)_SE']) & (d['LOG(OR)_SE'] > 0)
    d = d[ok].reset_index(drop=True)
    n = len(d)
    p = d.P.to_numpy()
    lam = lambda_gc(p)
    maf = np.minimum(d.A1_FREQ.to_numpy(), 1 - d.A1_FREQ.to_numpy())
    b = maf_bin(maf)
    by_maf = []
    for k, lab in enumerate(MAF_LABELS):
        sel = b == k
        by_maf.append(dict(bin=lab, lo=MAF_EDGES[k], hi=MAF_EDGES[k + 1], n=int(sel.sum()),
                           lambda_gc=lambda_gc(p[sel]) if sel.any() else None))
    n_sug, n_gw = int((p < SUGGEST_ALPHA).sum()), int((p < GW_ALPHA).sum())

    # the upstream scan_qc.tsv must describe the same scan
    q = pd.read_csv(a.scan_qc, sep='\t')
    q = q[q.model == a.upstream_model]
    assert len(q) == 1, f'{a.scan_qc}: {len(q)} rows for model {a.upstream_model}'
    q = q.iloc[0]
    assert q.cohort == a.cohort, (q.cohort, a.cohort)
    assert abs(float(q.lambda_gc) - lam) < 1e-4, f'lambda {lam:.5f} vs upstream {q.lambda_gc}'
    assert int(q.n_suggestive) == n_sug and int(q.n_genomewide) == n_gw, \
        (n_sug, q.n_suggestive, n_gw, q.n_genomewide)
    assert int(q.n_analysed) == n, (n, q.n_analysed)

    nlp = (-np.log10(p)).astype(np.float32)
    chrom = d['#CHROM'].map(chrom_key).to_numpy().astype(np.int16)
    np.savez(f'{a.out_prefix}.npz', chrom=chrom, pos=d.POS.to_numpy().astype(np.int64), nlp=nlp,
             qq=np.sort(nlp)[::-1].copy())
    top = d[d.P < a.top_p][['ID', '#CHROM', 'POS', 'A1_FREQ', 'OR', 'P']].rename(columns={'#CHROM': 'CHROM'})
    top.to_csv(f'{a.out_prefix}.top.tsv', sep='\t', index=False, float_format='%.6g')
    out = dict(cohort=a.cohort, model=a.model, upstream_model=a.upstream_model, n_all=n_all, n_excluded=n_all - n,
               n=n, lambda_gc=lam, lambda_by_maf=by_maf, n_p_lt_suggestive=n_sug, n_p_lt_genomewide=n_gw,
               expected_p_lt_suggestive=n * SUGGEST_ALPHA, expected_p_lt_genomewide=n * GW_ALPHA,
               min_p=float(p.min()), n_case=int(q.n_case), n_ctrl=int(q.n_ctrl), n_eff=float(q.n_eff),
               upstream_lambda_gc=float(q.lambda_gc))
    with open(f'{a.out_prefix}.stats.json', 'w') as fh:
        json.dump(out, fh, indent=2)
    log(f'{a.cohort} {a.model}: {n:,} variants ({n_all - n} excluded), lambda {lam:.4f}, '
        f'P<1e-5 {n_sug} (exp {n * SUGGEST_ALPHA:.0f}), P<5e-8 {n_gw}')


if __name__ == '__main__':
    main()
