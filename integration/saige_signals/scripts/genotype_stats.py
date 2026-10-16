#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Genotype counts, effect-allele frequency, missingness and HWE P of
#           every table variant in one cohort, for cases and controls apart.
#           Uses the component-shared variant_annot.group_stats (plink2
#           --geno-counts --hardy), the function behind assoc_saige's own lead
#           tables, on the genotype set the scan ran on.
# Inputs  : variants.tsv; the cohort's PLINK files (.fam column 6: 1 control, 2 case)
# Outputs : geno.<cohort>.tsv  one row per variant
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sig_common import import_annot, load_variants, log, write_tsv  # noqa: E402


def fetch_geno(v, cohort, bfile, plink2, threads, va):
    fam = pd.read_csv(f'{bfile}.fam', sep=r'\s+', header=None,
                      names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'],
                      dtype={'FID': str, 'IID': str})
    n_case, n_ctrl = int((fam.PHENO == 2).sum()), int((fam.PHENO == 1).sum())
    assert set(fam.PHENO) == {1, 2} and n_case > 0 and n_ctrl > 0, (cohort, set(fam.PHENO))
    bim_ids = set(pd.read_csv(f'{bfile}.bim', sep='\t', header=None, usecols=[1], dtype=str)[1])
    missing = set(v.variant_id) - bim_ids
    assert not missing, (cohort, missing)
    ids = sorted(v.variant_id)
    wd = Path('geno') / cohort
    case, nc = va.group_stats(plink2, str(bfile), fam, wd, 'case', 2, ids, threads=threads)
    ctrl, nk = va.group_stats(plink2, str(bfile), fam, wd, 'control', 1, ids, threads=threads)
    # group_stats returns {} silently when plink2 fails
    assert (nc, nk) == (n_case, n_ctrl) and len(case) == len(v) and len(ctrl) == len(v), (len(case), len(ctrl))
    out = []
    for _, r in v.iterrows():
        rec = {'variant_id': r.variant_id, 'cohort': cohort, 'n_case': n_case, 'n_ctrl': n_ctrl}
        for grp, d, n in (('case', case, n_case), ('ctrl', ctrl, n_ctrl)):
            g = d[r.variant_id]
            assert g['alt'] == r.EA and g['ref'] == r.OA, (r.variant_id, g['ref'], g['alt'])
            hom_oa, het, hom_ea, called = int(g['hom_ref']), int(g['het']), int(g['hom_alt']), int(g['called'])
            miss = n - called
            assert miss >= 0 and abs(g['missing_rate'] - miss / n) < 1e-9, (r.variant_id, cohort, grp)
            eaf = (2 * hom_ea + het) / (2 * called) if called else np.nan
            rec.update({f'{grp}_homOA': hom_oa, f'{grp}_het': het, f'{grp}_homEA': hom_ea,
                        f'{grp}_missing_n': miss, f'{grp}_called': called,
                        f'{grp}_geno': f'{hom_oa} / {het} / {hom_ea} / {miss}',
                        f'{grp}_eaf': eaf, f'{grp}_missing_rate': miss / n,
                        f'{grp}_hwe_p': float(g['hwe_p'])})
        out.append(rec)
    log(f'genotype statistics: {cohort} ({n_case} cases / {n_ctrl} controls, {len(v)} variants)')
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--variants', required=True)
    ap.add_argument('--cohort', required=True)
    ap.add_argument('--bfile', required=True, help='PLINK prefix')
    ap.add_argument('--plink2', required=True)
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--annot-lib', required=True, help='the staged _shared/scripts/variant_annot.py')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    va = import_annot(a.annot_lib)
    out = fetch_geno(load_variants(a.variants), a.cohort, a.bfile, a.plink2, a.threads, va)
    write_tsv(out, a.out)


if __name__ == '__main__':
    main()
