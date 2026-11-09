#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : One PLINK-ready VCF per (case set, cohort) holding every marker the
#           conditional analyses need together: the cohort's alleles AND
#           residues, the genome-wide scan's MHC lead variant(s), and a sentinel.
#
#           HLA markers are coded as in build_markers.py (REF A = absent, ALT
#           P = present; blank -> ./.) at the same positions (gene start + rank
#           within gene and class). The GWAS lead keeps its own REF/ALT and its
#           real GRCh38 position, read from the random_model genotypes; its
#           effect allele is the ALT, as in the genome-wide scan.
#
#           SENTINEL. SAIGE 1.1.9 step 2 run with --condition overwrites the
#           output row of the FIRST marker in the PLINK file with another
#           marker's statistics (seen here: the first residue took a DQB1
#           residue's allele count). The first record is therefore a sentinel --
#           a copy of the first HLA marker's genotypes under the ID SENTINEL at
#           position 1 -- whose row every reader discards. conditional.py also
#           checks, in every conditional run, that every marker's marginal
#           statistics equal the unconditional run's.
#
#           LD. pooled_markers.tsv carries, for every HLA marker, its Pearson r
#           with each lead's ALT dosage over the case set's samples (pairwise
#           complete) and over the PHENO1 = 0 group alone.
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import vocab as V                      # noqa: E402

GT = {'0': '0/0', '1': '0/1', '2': '1/1', '': './.'}
DOSE = {'0/0': 0.0, '0/1': 1.0, '1/0': 1.0, '1/1': 2.0, '0|0': 0.0, '0|1': 1.0, '1|0': 1.0, '1|1': 2.0}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case-set', required=True, choices=V.CASE_SETS)
    p.add_argument('--cohort', required=True)
    p.add_argument('--allele', required=True)
    p.add_argument('--residue', required=True)
    p.add_argument('--lead-vcf', required=True, help='the lead variant(s), all cohort samples, IID sample names')
    p.add_argument('--coords', required=True)
    p.add_argument('--pheno', required=True, help="pheno_covar.random.tsv of the case set")
    p.add_argument('--pheno-col', default='PHENO1')
    p.add_argument('--out-vcf', required=True)
    p.add_argument('--out-markers', required=True)
    return p.parse_args()


def read_lead_vcf(path):
    samples, recs = None, []
    for ln in open(path):
        if ln.startswith('##'):
            continue
        f = ln.rstrip('\n').split('\t')
        if ln.startswith('#'):
            samples = f[9:]
            continue
        gts = [g.split(':')[0] for g in f[9:]]
        recs.append({'chrom': f[0], 'pos': int(f[1]), 'id': f[2], 'ref': f[3], 'alt': f[4],
                     'gt': dict(zip(samples, gts))})
    if not recs:
        sys.exit(f'ERROR: no lead variant in {path}')
    return recs


def main():
    a = parse_args()
    pheno = pd.read_csv(a.pheno, sep='\t', dtype=str)
    ids = pheno['IID'].tolist()
    coords = pd.read_csv(a.coords, sep='\t', dtype={'gene': str}).set_index('gene')

    rows, dose = [], {}
    for cls, path in (('allele', a.allele), ('residue', a.residue)):
        m = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False).set_index('sample_id')
        missing = set(ids) - set(m.index)
        if missing:
            sys.exit(f'ERROR: {len(missing)} case-set samples are not {cls} matrix rows')
        m = m.loc[ids]
        info = []
        for name in m.columns:
            gene, pos, aa = V.parse_marker(name, cls)
            info.append((name, gene, pos, aa, int(coords.at[gene, 'order']), V.sort_key(name, cls)))
        info.sort(key=lambda t: (t[4], t[5]))
        rank = {}
        for name, gene, pos, aa, _, _ in info:
            rank[gene] = rank.get(gene, 0) + 1
            rows.append({'marker': name, 'class': cls, 'gene': gene, 'position': pos, 'residue': aa,
                         'chrom': '6', 'bp': int(coords.at[gene, 'gene_start']) + rank[gene],
                         'ref': 'A', 'alt': 'P', 'gts': [GT[v] for v in m[name].tolist()]})
            dose[name] = m[name].replace('', np.nan).astype(float).to_numpy()

    leads = read_lead_vcf(a.lead_vcf)
    for r in leads:
        absent = [s for s in ids if s not in r['gt']]
        if absent:
            sys.exit(f"ERROR: {len(absent)} case-set samples missing from the lead VCF for {r['id']}")
        gts = [r['gt'][s] if r['gt'][s] in DOSE else './.' for s in ids]
        rows.append({'marker': r['id'], 'class': 'gwas_lead', 'gene': '', 'position': None, 'residue': '',
                     'chrom': r['chrom'].replace('chr', ''), 'bp': r['pos'], 'ref': r['ref'], 'alt': r['alt'],
                     'gts': gts})
        dose[r['id']] = np.array([DOSE.get(g, np.nan) for g in gts])

    tab = pd.DataFrame(rows)
    if not tab['marker'].is_unique:
        sys.exit('ERROR: a marker ID occurs twice in the pooled set')
    tab = tab.sort_values(['bp', 'class', 'marker'], kind='stable').reset_index(drop=True)
    first = tab.iloc[0]

    with open(a.out_vcf, 'w') as fh:
        fh.write('##fileformat=VCFv4.2\n##source=assoc_hla.v2/build_pooled.py\n')
        fh.write(f'##assoc_hla_v2=case set {a.case_set}; cohort {a.cohort}; pooled alleles + residues + '
                 'GWAS lead(s); HLA REF A = absent, ALT P = present; SENTINEL = copy of the first HLA '
                 'marker, discarded by every reader\n')
        fh.write('##contig=<ID=6>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(ids) + '\n')
        fh.write(f"6\t1\tSENTINEL\t{first.ref}\t{first.alt}\t.\t.\t.\tGT\t" + '\t'.join(first.gts) + '\n')
        for r in tab.itertuples():
            fh.write(f'{r.chrom}\t{r.bp}\t{r.marker}\t{r.ref}\t{r.alt}\t.\t.\t.\tGT\t' + '\t'.join(r.gts) + '\n')

    # ── LD of every HLA marker with every lead ──────────────────────────────
    y = pheno[a.pheno_col].astype(int).to_numpy()
    out = tab.drop(columns='gts').copy()
    out.insert(0, 'case_set', a.case_set)
    out.insert(1, 'cohort', a.cohort)
    out['file_order'] = np.arange(1, len(out) + 1) + 1          # the sentinel is record 1

    def corr(x, z, mask):
        ok = mask & ~np.isnan(x) & ~np.isnan(z)
        if ok.sum() < 3 or np.std(x[ok]) == 0 or np.std(z[ok]) == 0:
            return np.nan
        return float(np.corrcoef(x[ok], z[ok])[0, 1])

    every = np.ones(len(y), dtype=bool)
    for r in leads:
        g = dose[r['id']]
        out[f"r_{r['id']}"] = [corr(dose[mk], g, every) for mk in out['marker']]
        out[f"r_group0_{r['id']}"] = [corr(dose[mk], g, y == 0) for mk in out['marker']]
    out['position'] = out['position'].astype('Int64')
    out.to_csv(a.out_markers, sep='\t', index=False, float_format='%.6g')
    n_hla = int((out['class'] != 'gwas_lead').sum())
    best = {r['id']: out.loc[out['class'] != 'gwas_lead', f"r_{r['id']}"].abs().max() for r in leads}
    print(f'[build_pooled] {a.case_set} {a.cohort}: {n_hla} HLA markers + {len(leads)} lead(s), '
          f'{len(ids)} samples; max |r| with lead: ' + ', '.join(f'{k} {v:.3f}' for k, v in best.items()))


if __name__ == '__main__':
    main()
