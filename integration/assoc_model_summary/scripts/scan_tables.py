#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The calibration numbers of a set of scans, as tables:
#             tail_counts.tsv    per cohort x model: samples, lambda_GC, variants
#                                at P < 1e-5 and < 5e-8 against expectation,
#                                min P, genome-wide / suggestive loci (upstream
#                                peak calls) and the genome-wide lead labels
#             lambda_by_maf.tsv  per cohort x model x MAF bin: lambda_GC
#             gw_leads.tsv       every genome-wide lead (model, cohort, gene, variant, P)
# Inputs  : <cohort>.<model>.stats.json (reduce_scan.py); upstream
#           model_peaks_annotation.tsv per cohort
# Outputs : tail_counts.tsv, lambda_by_maf.tsv, gw_leads.tsv
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402


def peak_counts(files, cohorts, upstream_model, model, leads):
    out = {}
    for c, f in zip(cohorts, C.csv_list(files)):
        d = pd.read_csv(f, sep='\t', dtype={'chrom': str})
        d = d[d.model == upstream_model]
        gw = d[d.tier == 'genome_wide']
        leads += [dict(model=model, cohort=c, Gene=r.Gene, rsID=r.rsID, variant_id=r.variant_id, P=r.P, OR=r.OR)
                  for _, r in gw.iterrows()]
        labels = [g if str(g) not in ('.', '', 'nan') else r for g, r in zip(gw.Gene.astype(str), gw.rsID.astype(str))]
        out[c] = dict(gw_loci=len(gw), sg_loci=int((d.tier == 'suggestive').sum()),
                      gw_leads=', '.join(labels))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cohorts', required=True)
    ap.add_argument('--models', required=True)
    ap.add_argument('--peaks', action='append', required=True, metavar='MODEL=FILE,FILE,FILE')
    ap.add_argument('--upstream-models', required=True)
    ap.add_argument('--out-tail', required=True)
    ap.add_argument('--out-maf', required=True)
    ap.add_argument('--out-leads', required=True)
    a = ap.parse_args()
    cohorts, models = C.csv_list(a.cohorts), C.csv_list(a.models)
    upstream, pk = C.kv(a.upstream_models), C.kv(a.peaks)
    rows, mrows, leads = [], [], []
    for m in models:
        counts = peak_counts(pk[m], cohorts, upstream[m], m, leads)
        for c in cohorts:
            s = C.load_json(f'{c}.{m}.stats.json')
            rows.append(dict(cohort=c, model=m, model_label=C.MODELS[m]['short'], n_case=s['n_case'],
                             n_ctrl=s['n_ctrl'], n_eff=s['n_eff'], n_variants=s['n'], lambda_gc=s['lambda_gc'],
                             n_p_lt_1e5=s['n_p_lt_suggestive'], expected_p_lt_1e5=s['expected_p_lt_suggestive'],
                             ratio_p_lt_1e5=s['n_p_lt_suggestive'] / s['expected_p_lt_suggestive'],
                             n_p_lt_5e8=s['n_p_lt_genomewide'], min_p=s['min_p'], **counts[c]))
            for b in s['lambda_by_maf']:
                mrows.append(dict(cohort=c, model=m, maf_bin=b['bin'], maf_lo=b['lo'], maf_hi=b['hi'], n=b['n'],
                                  lambda_gc=b['lambda_gc']))
    pd.DataFrame(rows).to_csv(a.out_tail, sep='\t', index=False, float_format='%.6g')
    pd.DataFrame(mrows).to_csv(a.out_maf, sep='\t', index=False, float_format='%.6g')
    pd.DataFrame(leads, columns=['model', 'cohort', 'Gene', 'rsID', 'variant_id', 'P', 'OR']) \
        .to_csv(a.out_leads, sep='\t', index=False, float_format='%.6g')
    C.log(f'tables for {models}: {len(rows)} scans')


if __name__ == '__main__':
    main()
