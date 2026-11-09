#!/usr/bin/env python3
"""An arm's HLA association matrices on the ORIGINAL tested columns.

hla.typing.v2 step 06 re-applies its floors (call rate >= 95 % in cases and in controls, MAF >= 0.5 %)
to whatever it is given, so re-typed samples could change which columns are tested and with them the
GEC thresholds. The re-test keeps the main analysis' tests: the columns of the published
06.assoc_inputs/<cohort>/<class>_dosage.tsv, in its order, taken from the arm's 05 exports, rows in
the cohort's phenotype order, written as step 06 writes them (Int64, blank = undetermined). On the
ORIGINAL exports this reproduces the published files byte for byte (checked in --check mode).
"""
import argparse
import hashlib
import sys
from pathlib import Path

import pandas as pd

V2 = Path('/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/hla.typing.v2')
MODEL_INPUTS = Path('/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/12_model_inputs')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--exports', required=True, help="the arm's 05.exports directory")
    ap.add_argument('--cohort', required=True)
    ap.add_argument('--class', dest='cls', required=True, choices=['allele', 'residue'])
    ap.add_argument('--out', required=True)
    ap.add_argument('--check', action='store_true', help='abort unless identical to the published matrix')
    a = ap.parse_args()

    pub = V2 / 'results' / '06.assoc_inputs' / a.cohort / f'{a.cls}_dosage.tsv'
    cols = pub.open().readline().rstrip('\n').split('\t')[1:]
    src = pd.read_csv(Path(a.exports) / f'{a.cls}_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    missing = [c for c in cols if c not in src.columns]
    if missing:
        sys.exit(f'ABORT: {len(missing)} tested columns are not in the exports, e.g. {missing[:3]}')
    ph = pd.read_csv(MODEL_INPUTS / a.cohort / 'phenotype' / 'pheno.tsv', sep='\t', dtype={'IID': str})
    d = src.loc[ph.IID.tolist(), cols].astype('Int64')
    d.index.name = 'sample_id'
    out = Path(a.out)
    d.to_csv(out, sep='\t', na_rep='')
    m = hashlib.md5(out.read_bytes()).hexdigest()
    ref = hashlib.md5(pub.read_bytes()).hexdigest()
    print(f'[hla_matrices] {a.cohort} {a.cls}: {d.shape[0]} samples x {d.shape[1]} columns; md5 {m} '
          f'({"identical to" if m == ref else "differs from"} the published matrix)')
    if a.check and m != ref:
        sys.exit('ABORT: the ORIGINAL exports do not reproduce the published matrix')


if __name__ == '__main__':
    main()
