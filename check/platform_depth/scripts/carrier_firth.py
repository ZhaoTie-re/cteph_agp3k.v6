#!/usr/bin/env python3
"""Per-subset carrier odds ratio of each Tier 1-2 gene (descriptive; CMC's collapsed indicator).

A sample carries the gene when it holds the ORIGINAL minor allele at one or more of the gene's tested
variants in that variant set (missing genotypes count as non-carrying, as the collapsed indicator of
CMC does after mean imputation rounds to zero for rare variants). plink2 --glm firth with sex + 10 PCs
on the subset's rare-variant samples gives a finite OR even when a platform stratum has no carrier;
with 0-2 carriers per stratum these estimates are descriptive only.
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bedio as B  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bfile', required=True, help='the ORIGINAL gene fileset of the cohort')
    ap.add_argument('--sets', nargs='+', required=True, help='<stratum>=<set.txt from gene_retest.py prepare>')
    ap.add_argument('--pheno', required=True, help='pheno_covar.fixed.tsv of the pooled cohort (#FID IID PHENO1 SEX PC..)')
    ap.add_argument('--include', nargs='+', required=True, help='<subset>=<rv_include.txt>')
    ap.add_argument('--plink2', required=True)
    ap.add_argument('--covar-cols', required=True)
    ap.add_argument('--cohort', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    fam = B.read_fam(a.bfile)
    bim = B.read_bim(a.bfile)
    codes = B.unpack(B.read_rows(a.bfile, list(range(len(bim))), len(fam)), len(fam))
    ph = pd.read_csv(a.pheno, sep='\t', dtype=str)
    in_ph = np.array([s in set(ph.IID) for s in fam])
    pos = {f[1]: (f[1].split(':')[0], int(f[1].split(':')[1])) for f in bim}

    ind = {}
    for spec in a.sets:
        stratum, path = spec.split('=', 1)
        for ln in Path(path).read_text().splitlines():
            if not ln.strip():
                continue
            gene, rng = ln.split('\t')
            want = {(f'chr{r.split(":")[0]}', int(r.split(':')[1].split('-')[0])) for r in rng.split(',')}
            car = np.zeros(len(fam), bool)
            for i, f in enumerate(bim):
                if pos[f[1]] not in want:
                    continue
                c = codes[i]
                n5 = 2 * (c[in_ph] == B.HOM5).sum() + (c[in_ph] == B.HET).sum()
                n6 = 2 * (c[in_ph] == B.HOM6).sum() + (c[in_ph] == B.HET).sum()
                minor_hom = B.HOM5 if n5 <= n6 else B.HOM6
                car |= (c == B.HET) | (c == minor_hom)
            ind[f'{gene}|{stratum}'] = car
    names = list(ind)
    # a PLINK fileset of the indicators: carrier = heterozygous for allele 'C' (column 5)
    codes_ind = np.vstack([np.where(ind[n], B.HET, B.HOM6).astype(np.uint8) for n in names])
    with open('carriers.bim', 'w') as fh:
        for k, n in enumerate(names):
            fh.write(f'1\t{n.replace("|", "__")}\t0\t{k + 1}\tC\tA\n')
    with open('carriers.bed', 'wb') as fh:
        fh.write(B.MAGIC)
        fh.write(B.pack(codes_ind).tobytes())
    Path('carriers.fam').write_text(Path(f'{a.bfile}.fam').read_text())

    rows = []
    for spec in a.include:
        subset, path = spec.split('=', 1)
        inc = set(Path(path).read_text().split())
        keep = ph[ph.IID.isin(inc)][['#FID', 'IID']]
        keep.to_csv('keep.txt', sep='\t', index=False, header=False)
        subprocess.run([a.plink2, '--bfile', 'carriers', '--keep', 'keep.txt', '--pheno', a.pheno, '--pheno-name',
                        'PHENO1', '--1', '--covar', a.pheno, '--covar-name', a.covar_cols, '--glm', 'firth',
                        'hide-covar', 'cols=+a1countcc,+totallelecc', '--ci', '0.95', '--out', f'firth.{subset}',
                        '--threads', '1'], check=True, capture_output=True)
        g = pd.read_csv(f'firth.{subset}.PHENO1.glm.firth', sep='\t', dtype={'ID': str})
        case_ids = set(ph[(ph.IID.isin(inc)) & (ph.PHENO1 == '1')].IID)
        ctrl_ids = set(ph[(ph.IID.isin(inc)) & (ph.PHENO1 == '0')].IID)
        if set(g.A1) != {'C'}:
            sys.exit(f'ABORT: plink2 tested another allele than the carrier allele: {set(g.A1)}')
        for _, r in g.iterrows():
            name = r.ID.replace('__', '|')
            gene, stratum = name.split('|')
            car = ind[name]
            rows.append({'cohort': a.cohort, 'subset': subset, 'gene': gene, 'stratum': stratum,
                         'carriers_case': int(sum(car[i] for i, s in enumerate(fam) if s in case_ids)),
                         'carriers_control': int(sum(car[i] for i, s in enumerate(fam) if s in ctrl_ids)),
                         'n_case': len(case_ids), 'n_control': len(ctrl_ids),
                         'or_firth': r['OR'], 'l95': r['L95'], 'u95': r['U95'], 'se_log_or': r['LOG(OR)_SE'],
                         'p_firth': r['P'], 'errcode': r.get('ERRCODE', '.')})
    out = pd.DataFrame(rows)
    out.to_csv(a.out, sep='\t', index=False)
    print(out.to_string(index=False))


if __name__ == '__main__':
    main()
