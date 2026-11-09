#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : One cohort x class matrix -> a VCF of hard genotype calls that both
#           engines read (via one PLINK set), plus the per-marker counts every
#           result row is reported with.
#
#           CODING. Each matrix cell is the number of a sample's two copies of
#           the gene that carry the marker (an allele, or one amino acid at one
#           IMGT position): 0 -> 0/0, 1 -> 0/1, 2 -> 1/1, blank -> ./. (the
#           copy was undetermined at typing; never imputed here). REF is 'A'
#           (absent) and ALT is 'P' (present), so every effect is per copy
#           carrying the marker, in both engines.
#
#           POSITIONS. Markers are placed on chromosome 6 at gene start + rank
#           of the marker within its gene (residues by IMGT position, then
#           residue letter; alleles by their numeric fields). The positions are
#           unique labels in chromosomal gene order, not physical positions of
#           the residues; LOCO in the random model needs only the chromosome.
#
#           SAMPLES. The PLINK set holds exactly the case set's random-model
#           samples (all controls + the set's cases), in phenotype-file order.
#
#           COUNTS. marker_info.tsv is written once per sample set -- 'random'
#           (all samples of the set) and 'fixed' (its fixed_model keep list) -- so the
#           counts reported beside a P value are those of the samples that
#           model analysed. A marker's minor allele count (mac) is over the
#           called copies of that set.
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


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort', required=True)
    p.add_argument('--case-set', required=True, choices=V.CASE_SETS)
    p.add_argument('--class', dest='cls', required=True, choices=V.CLASSES)
    p.add_argument('--matrix', required=True)
    p.add_argument('--coords', required=True, help='hla_gene_coords.tsv')
    p.add_argument('--pheno', required=True, help="pheno_covar.random.tsv (the case set's samples)")
    p.add_argument('--pheno-col', default='PHENO1')
    p.add_argument('--fixed-keep', required=True)
    p.add_argument('--out-vcf', required=True)
    p.add_argument('--out-info', required=True)
    return p.parse_args()


def counts(dose, y, rows):
    """Per-marker counts over the samples in `rows` (boolean mask)."""
    out = {}
    for grp, sel in (('case', rows & (y == 1)), ('control', rows & (y == 0))):
        d = dose[sel]
        called = (~np.isnan(d)).sum(axis=0)
        copies = np.nansum(d, axis=0)
        out[f'n_{grp}'] = np.full(dose.shape[1], int(sel.sum()))
        out[f'called_{grp}'] = called
        out[f'copies_{grp}'] = copies.astype(int)
        with np.errstate(invalid='ignore', divide='ignore'):
            out[f'freq_{grp}'] = copies / (2 * called)
            out[f'missing_{grp}'] = 1 - called / sel.sum()
    tot_called = out['called_case'] + out['called_control']
    tot_copies = out['copies_case'] + out['copies_control']
    out['freq'] = tot_copies / (2 * tot_called)
    out['mac'] = np.minimum(tot_copies, 2 * tot_called - tot_copies).astype(int)
    return out


def main():
    a = parse_args()
    m = pd.read_csv(a.matrix, sep='\t', dtype=str, keep_default_na=False)
    pheno = pd.read_csv(a.pheno, sep='\t', dtype=str)
    # The case set's samples, in the matrix's (= the cohort phenotype file's) order.
    missing = set(pheno['IID']) - set(m['sample_id'])
    if missing:
        sys.exit(f'ERROR: {len(missing)} case-set samples are not matrix rows')
    m = m[m['sample_id'].isin(set(pheno['IID']))].reset_index(drop=True)
    if m['sample_id'].tolist() != pheno['IID'].tolist():
        sys.exit('ERROR: matrix rows are not the case-set samples in order')
    keep = set(pd.read_csv(a.fixed_keep, sep='\t', header=None, dtype=str)[1])
    coords = pd.read_csv(a.coords, sep='\t', dtype={'gene': str}).set_index('gene')

    # ── marker order and positions ──────────────────────────────────────────
    rows = []
    for name in m.columns[1:]:
        gene, pos, aa = V.parse_marker(name, a.cls)
        if gene not in coords.index:
            sys.exit(f'ERROR: gene {gene} of {name} has no coordinates')
        rows.append({'marker': name, 'gene': gene, 'position': pos, 'residue': aa,
                     'gene_order': int(coords.at[gene, 'order']),
                     'key': V.sort_key(name, a.cls)})
    info = pd.DataFrame(rows).sort_values(['gene_order', 'key'], kind='stable').reset_index(drop=True)
    info['rank'] = info.groupby('gene').cumcount() + 1
    info['chrom'] = '6'
    info['bp'] = [int(coords.at[g, 'gene_start']) + r for g, r in zip(info.gene, info['rank'])]
    if not info['bp'].is_unique:
        sys.exit('ERROR: two markers share a position; gene starts are too close for the rank offsets')
    if not info['bp'].is_monotonic_increasing:
        sys.exit('ERROR: positions are not increasing in the VCF order')

    # ── VCF of hard calls ───────────────────────────────────────────────────
    samples = m['sample_id'].tolist()
    with open(a.out_vcf, 'w') as fh:
        fh.write('##fileformat=VCFv4.2\n')
        fh.write('##source=assoc_hla.v2/build_markers.py\n')
        fh.write(f'##assoc_hla_v2=case set {a.case_set}; cohort {a.cohort}; class {a.cls}; REF A = copy does not carry '
                 'the marker, ALT P = copy carries it; ./. = undetermined at typing\n')
        fh.write('##contig=<ID=6>\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples) + '\n')
        for r in info.itertuples():
            gts = '\t'.join(GT[v] for v in m[r.marker].tolist())
            fh.write(f'6\t{r.bp}\t{r.marker}\tA\tP\t.\t.\t.\tGT\t{gts}\n')

    # ── per-marker counts, per sample set ───────────────────────────────────
    dose = m[info['marker']].replace('', np.nan).astype(float).to_numpy()
    y = pheno[a.pheno_col].astype(int).to_numpy()
    sets = {'random': np.ones(len(m), dtype=bool),
            'fixed': m['sample_id'].isin(keep).to_numpy()}
    frames = []
    for sname, mask in sets.items():
        c = counts(dose, y, mask)
        f = info[['marker', 'gene', 'position', 'residue', 'chrom', 'bp', 'rank']].copy()
        f.insert(0, 'class', a.cls)
        f.insert(0, 'cohort', a.cohort)
        f.insert(0, 'case_set', a.case_set)
        f['sample_set'] = sname
        for k in ('n_case', 'n_control', 'called_case', 'called_control', 'copies_case',
                  'copies_control', 'freq_case', 'freq_control', 'freq', 'mac',
                  'missing_case', 'missing_control'):
            f[k] = c[k]
        frames.append(f)
    out = pd.concat(frames, ignore_index=True)
    out['position'] = out['position'].astype('Int64')
    out.to_csv(a.out_info, sep='\t', index=False, float_format='%.6g')
    print(f'[build_markers] {a.case_set} {a.cohort} {a.cls}: {len(info)} markers, {len(samples)} samples, '
          f'{int(sets["fixed"].sum())} in the fixed set; blank cells {int(np.isnan(dose).sum())}')


if __name__ == '__main__':
    main()
