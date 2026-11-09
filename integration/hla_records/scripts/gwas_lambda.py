#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : lambda_GC of the genome-wide random-model scan (the scan the GWAS MHC
#           lead comes from), for one cohort: SNPs with MAF >= 0.5% inside the MHC
#           window chr6:29.6-33.4 Mb and all others, each P read as a 1-df
#           chi-square (the association deck's definition, slides 23-24).
# Inputs  : --sumstats (genome-wide SAIGE sumstats.tsv of the cohort), --cohort
# Outputs : --out (json: mhc, rest, mhc_snps, rest_snps)
# Component: hla_records
# ---------------------------------------------------------------------------
import argparse
import json

import numpy as np
import pandas as pd
from scipy import stats

MHC_SPAN = (29.6e6, 33.4e6)
MIN_MAF = 0.005


def lam_gc(p):
    p = np.asarray(p, dtype=float)
    p = p[np.isfinite(p)]
    return float(np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--sumstats', required=True)
    ap.add_argument('--cohort', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    g = pd.read_csv(a.sumstats, sep='\t', usecols=['#CHROM', 'POS', 'A1_FREQ', 'P'], dtype={'#CHROM': str})
    g = g.dropna(subset=['P'])
    g = g[np.minimum(g.A1_FREQ, 1 - g.A1_FREQ) >= MIN_MAF]
    inm = (g['#CHROM'] == '6') & g.POS.between(*MHC_SPAN)
    out = {a.cohort: {'mhc': lam_gc(g.P[inm]), 'rest': lam_gc(g.P[~inm]), 'mhc_snps': int(inm.sum()),
                      'rest_snps': int((~inm).sum()), 'mhc_span_mb': [MHC_SPAN[0] / 1e6, MHC_SPAN[1] / 1e6],
                      'min_maf': MIN_MAF}}
    with open(a.out, 'w') as fh:
        json.dump(out, fh, indent=1)
    print(f'[gwas_lambda] {a.cohort}: MHC {out[a.cohort]["mhc"]:.3f} ({out[a.cohort]["mhc_snps"]} SNPs), '
          f'rest {out[a.cohort]["rest"]:.3f}')


if __name__ == '__main__':
    main()
