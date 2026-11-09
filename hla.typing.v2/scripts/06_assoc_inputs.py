#!/usr/bin/env python3
"""Step 06 — the matrices the association analysis reads, one set per cohort. No test is run here.

Input: the step-05 exports (sample x 2-field allele, sample x residue; 0/1/2 copies; blank = no
genotype, including a residue IPD-IMGT/HLA leaves undetermined). For each cohort of the association
(config.COHORTS, phenotype files of wgs.auto.par), every exported column is measured within the
cohort and kept or dropped by these rules, in this order; a column's `status` is the first it fails:

  not_varying     no copy of its less common form in the cohort (it cannot inform a test)
  excluded_gene   DRB3, DRB4 or DRB5 (presence follows DRB1)
  call_rate       share with a genotype at the marker < MIN_CALL_RATE in cases or in controls (measured in
                  each group, because the sequencing platform and depth differ between them)
  maf             copies of the less common form over the copies typed < MIN_MAF
  tested          passes all four

Samples are the cohort's phenotype file, in its order; every one must be typed. No sample is
dropped: the rules act on columns only.

Outputs (results/06.assoc_inputs/):
  markers.tsv                   one row per cohort and exported column: class, gene, position,
                                residue, samples, called, call rate (all, cases, controls), copies,
                                MAC, MAF, status
  qc_record.tsv                 columns left after each step, per cohort, gene and class
  summary.json                  thresholds, and per cohort the samples, cases, controls, columns
                                tested, and the md5 of each matrix
  <cohort>/allele_dosage.tsv    cohort samples x tested allele columns, 0/1/2, blank = no genotype
  <cohort>/residue_dosage.tsv   cohort samples x tested residue columns, the same
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402

STEPS = [('exported', 'varies in the cohort'),
         ('genes', 'DRB3, DRB4 and DRB5 left out'),
         ('call_rate', f'call rate >= {C.MIN_CALL_RATE:g} in cases and in controls'),
         ('maf', f'MAF >= {C.MIN_MAF:g}')]
# a column is still in after step k when its status is not one of the first k failures
FAILS = ['not_varying', 'excluded_gene', 'call_rate', 'maf']


def describe(cols, kind):
    """gene, position and residue of each column name."""
    if kind == 'allele':
        return pd.DataFrame({'gene': cols.str.split('*').str[0], 'position': '', 'residue': ''}, index=cols)
    parts = cols.str.split(':')
    return pd.DataFrame({'gene': parts.str[0], 'position': parts.str[1], 'residue': parts.str[2]}, index=cols)


def measure(d, kind, cohort, case):
    """Per column of one class, within one cohort: the counts, the floors' measures, and the status."""
    n = len(d)
    called = d.notna().sum()
    copies = d.sum()
    an = 2 * called
    mac = np.minimum(copies, an - copies)
    m = describe(d.columns, kind)
    m.insert(0, 'class', kind)
    m.insert(0, 'cohort', cohort)
    m['samples'] = n
    m['called'] = called.astype(int)
    m['call_rate'] = (called / n).round(6)
    m['call_rate_cases'] = d[case].notna().mean().round(6)
    m['call_rate_controls'] = d[~case].notna().mean().round(6)
    m['copies'] = copies.astype(int)
    m['mac'] = mac.astype(int)
    m['maf'] = (mac / an.where(an > 0)).round(6)
    status = pd.Series('tested', index=d.columns)
    status[m.maf < C.MIN_MAF] = 'maf'
    status[(m.call_rate_cases < C.MIN_CALL_RATE) | (m.call_rate_controls < C.MIN_CALL_RATE)] = 'call_rate'
    status[~m.gene.isin(C.ASSOC_GENES)] = 'excluded_gene'
    status[m.mac < 1] = 'not_varying'
    m['status'] = status
    m.index.name = 'marker'
    return m


def write_matrix(d, path):
    d = d.astype('Int64')
    d.index.name = 'sample_id'
    d.to_csv(path, sep='\t', na_rep='')
    return hashlib.md5(path.read_bytes()).hexdigest()


def main():
    out = C.D_ASSOC
    out.mkdir(parents=True, exist_ok=True)
    ad = pd.read_csv(C.D_EXPORTS / 'allele_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    rd = pd.read_csv(C.D_EXPORTS / 'residue_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    markers, record, summary = [], [], {
        'thresholds': {'min_call_rate': C.MIN_CALL_RATE, 'min_maf': C.MIN_MAF},
        'genes': C.ASSOC_GENES, 'rules': [r for _, r in STEPS], 'cohorts': {}}

    for co in C.COHORTS:
        ph = pd.read_csv(C.MODEL_INPUTS / co / 'phenotype' / 'pheno.tsv', sep='\t', dtype={'IID': str})
        missing = sorted(set(ph.IID) - set(ad.index))
        if missing:
            sys.exit(f'{co}: {len(missing)} samples of the phenotype file are not typed')
        ids = ph.IID.tolist()
        case = (ph.PHENO1 == 1).values
        m = pd.concat([measure(ad.loc[ids], 'allele', co, case), measure(rd.loc[ids], 'residue', co, case)])
        markers.append(m)

        for k, (step, _) in enumerate(STEPS):
            left = m[~m.status.isin(FAILS[:k + 1])]
            for (g, kind), v in left.groupby(['gene', 'class']).size().items():
                record.append({'cohort': co, 'step': step, 'gene': g, 'class': kind, 'columns': int(v)})

        (out / co).mkdir(exist_ok=True)
        keep = m[m.status == 'tested']
        md5 = {}
        for kind, src, name in (('allele', ad, 'allele_dosage.tsv'), ('residue', rd, 'residue_dosage.tsv')):
            cols = keep[keep['class'] == kind].index
            md5[name] = write_matrix(src.loc[ids, cols], out / co / name)
        summary['cohorts'][co] = {
            'samples': len(ids), 'cases': int((ph.PHENO1 == 1).sum()), 'controls': int((ph.PHENO1 == 0).sum()),
            'alleles_tested': int((keep['class'] == 'allele').sum()),
            'residues_tested': int((keep['class'] == 'residue').sum()), 'md5': md5}
        print(f'{co}: {len(ids)} samples; tested {summary["cohorts"][co]["alleles_tested"]} alleles, '
              f'{summary["cohorts"][co]["residues_tested"]} residues')

    pd.concat(markers).to_csv(out / 'markers.tsv', sep='\t')
    pd.DataFrame(record).to_csv(out / 'qc_record.tsv', sep='\t', index=False)
    (out / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')


if __name__ == '__main__':
    main()
