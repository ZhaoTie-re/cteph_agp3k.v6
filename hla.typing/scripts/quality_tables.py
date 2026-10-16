#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : One quality row per column of the two matrices an association reads.
#
#           allele_quality.tsv  : one row per column of allele_dosage.tsv
#           residue_quality.tsv : one row per column of residue_dosage.tsv
#
#           Everything needed to decide whether a column is worth testing was
#           already computable from the delivered tables, and was therefore
#           computed over and over by whoever opened them -- keyed by P group in
#           allele_frequency_check.tsv rather than by the column name that is
#           actually regressed on, and, for residues, not published at all: the
#           only per-position measurement in the component lived inside
#           build_facts.py and was thrown away as four summary numbers.
#
#           These tables carry NO verdict column. They report the chromosome
#           counts and the flags; where the line between usable and not lies is
#           the analyst's decision and depends on the model being fitted, and a
#           `tier` column here would be this component quietly making it.
#
#           Nothing is recomputed that another script already computes:
#           wilson() comes from allele_freq_check.py and the reference index from
#           residue_matrix.py. A second copy of either is how the number in the
#           prose and the number in the pipeline come to disagree.
# Component: hla.typing
# ---------------------------------------------------------------------------
import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SAMPLE_ID = {'sample_id': str}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--allele-dosage', required=True, help='allele_dosage.tsv')
    p.add_argument('--residue-dosage', required=True, help='residue_dosage.tsv')
    p.add_argument('--pgroup-map', required=True, help='allele_pgroup_map.tsv')
    p.add_argument('--ambiguity', required=True, help='allele_ambiguity.tsv')
    p.add_argument('--gene-qc', required=True, help='typing_qc_gene.tsv')
    p.add_argument('--summary', required=True,
                   help='allele_frequency_summary.tsv, for the comparable flag')
    p.add_argument('--manifest', required=True, help='sample_manifest.tsv')
    p.add_argument('--reference-dir', required=True, help='residue_reference/')
    p.add_argument('--control-group', default='AGP3K',
                   help='the manifest group the frequencies are reported for')
    p.add_argument('--no-dosage-genes', default='DRB3,DRB4,DRB5',
                   help='genes whose dosage columns are wrong by construction, '
                        'because a hemizygous call is written as two copies; see '
                        'docs/OPEN_QUESTIONS.md section 1')
    p.add_argument('--field-depth', type=int, default=2,
                   help='field depth of the dosage columns, for matching the '
                        'ambiguity table back onto them')
    p.add_argument('--out-allele', default='allele_quality.tsv')
    p.add_argument('--out-residue', default='residue_quality.tsv')
    return p.parse_args()


def sibling(name):
    """Import another script from this directory rather than copying it."""
    path = Path(__file__).resolve().parent / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    argv, sys.argv = sys.argv, [str(path)]
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.argv = argv
    return mod


def numeric(df):
    """The matrices write '' where a gene was not typed. That is not a zero."""
    return df.apply(pd.to_numeric, errors='coerce')


def allele_table(a):
    afc = sibling('allele_freq_check.py')
    dos = pd.read_csv(a.allele_dosage, sep='\t', dtype=SAMPLE_ID).set_index('sample_id')
    dos = numeric(dos)
    man = pd.read_csv(a.manifest, sep='\t', dtype=str).set_index('sample_id')
    ctrl = dos.index.isin(man.index[man['group'] == a.control_group])
    case = ~ctrl

    pg = pd.read_csv(a.pgroup_map, sep='\t', dtype=str).set_index('allele')
    gq = pd.read_csv(a.gene_qc, sep='\t').set_index('gene')
    sm = pd.read_csv(a.summary, sep='\t').set_index('gene')
    no_dosage = {g.strip() for g in a.no_dosage_genes.split(',') if g.strip()}

    # A chromosome carrying this allele came from a call the program could not
    # resolve on its own. Counted per allele, from the per-sample table that is
    # never published: docs/OPEN_QUESTIONS.md section 5 is why it is worth knowing.
    amb = pd.read_csv(a.ambiguity, sep='\t', dtype=str)
    amb_n = {}
    for col in ('reported_1', 'reported_2'):
        if col not in amb.columns:
            continue
        for v in amb[col].dropna():
            k = afc.trim(str(v), a.field_depth)
            amb_n[k] = amb_n.get(k, 0) + 1

    # The denominator is the gene's own called chromosomes, taken from the matrix
    # itself: an allele's columns sum to 2 per typed sample, so the gene's columns
    # summed over a group ARE that group's called chromosomes at that gene.
    gene_of = {c: c.split('*', 1)[0] for c in dos.columns}
    by_gene = {}
    for c, g in gene_of.items():
        by_gene.setdefault(g, []).append(c)
    denom = {g: (float(dos.loc[ctrl, cols].sum().sum()),
                 float(dos.loc[case, cols].sum().sum()))
             for g, cols in by_gene.items()}

    rows = []
    for col in dos.columns:
        g = gene_of[col]
        nc = float(dos.loc[ctrl, col].sum())
        na = float(dos.loc[case, col].sum())
        d_ctrl = denom[g][0]
        freq = nc / d_ctrl if d_ctrl else float('nan')
        lo, hi = afc.wilson(nc, d_ctrl) if d_ctrl else (float('nan'), float('nan'))
        in_ref = pg['in_reference'].get(col, '')
        freq_ref = pg['freq_reference'].get(col, '')
        n_amb = amb_n.get(col, 0)
        total = nc + na
        comparable = ('' if g not in sm.index
                      else int(sm.at[g, 'comparable']) if 'comparable' in sm.columns else 1)
        valid = 0 if g in no_dosage else 1
        note = []
        if str(in_ref) == '0':
            note.append('no published panel lists this allele')
        elif in_ref == '':
            note.append('locus not covered by the reference panel')
        if comparable == 0:
            note.append('locus not comparable against the panel')
        if not valid:
            note.append('hemizygous written as two copies; do not use this dosage')
        if n_amb:
            note.append(f'{n_amb} chromosome(s) from a resolved ambiguity')
        rows.append({
            'gene': g,
            'allele': col,
            'p_group': pg['p_group'].get(col, ''),
            'n_chr_ctrl': int(nc),
            'n_chr_case': int(na),
            'n_chr_total': int(total),
            'n_chr_ctrl_gene': int(d_ctrl),
            'freq_ctrl': round(freq, 6) if d_ctrl else '',
            'ci_lo_ctrl': round(lo, 6) if d_ctrl else '',
            'ci_hi_ctrl': round(hi, 6) if d_ctrl else '',
            'in_reference': in_ref,
            'freq_reference': freq_ref,
            'n_chr_ambiguous': n_amb,
            'share_ambiguous': round(n_amb / total, 4) if total else '',
            'gene_call_rate': gq['call_rate'].get(g, ''),
            'gene_comparable': comparable,
            'gene_dosage_valid': valid,
            'note': '; '.join(note),
        })
    out = pd.DataFrame(rows).sort_values(['gene', 'allele'])
    out.to_csv(a.out_allele, sep='\t', index=False)
    print(f'[quality_tables] {len(out):,} allele row(s) -> {a.out_allele}')
    return out


def residue_table(a):
    rm = sibling('residue_matrix.py')
    dos = pd.read_csv(a.residue_dosage, sep='\t', dtype=SAMPLE_ID).set_index('sample_id')
    dos = numeric(dos)
    man = pd.read_csv(a.manifest, sep='\t', dtype=str).set_index('sample_id')
    ctrl = dos.index.isin(man.index[man['group'] == a.control_group])
    case = ~ctrl
    no_dosage = {g.strip() for g in a.no_dosage_genes.split(',') if g.strip()}

    # Which of THIS cohort's alleles carry each amino acid. This is the quantity
    # the report's claim about residue tests rests on -- that a mis-call to a
    # neighbour usually does not move the residue -- made checkable per residue
    # instead of asserted in general.
    ad = pd.read_csv(a.allele_dosage, sep='\t', nrows=1)
    carried = {}
    for g in sorted({c.split(':', 1)[0] for c in dos.columns}):
        ref_path = Path(a.reference_dir) / f'{g}.tsv'
        if not ref_path.exists():
            continue
        ref = pd.read_csv(ref_path, sep='\t', index_col=0, dtype=str)
        exact, by_prefix = rm.build_index(ref)
        for al in [c for c in ad.columns if c.split('*', 1)[0] == g]:
            m, _how = rm.resolve(al, exact, by_prefix)
            if m is None:
                continue
            row = ref.loc[m]
            for pos, aa in row.items():
                carried[f'{g}:{pos}:{aa}'] = carried.get(f'{g}:{pos}:{aa}', 0) + 1

    # A position whose dosages do not sum to 2 has a chromosome the reference
    # could not place. build_facts.py used to compute this and publish only the
    # summary; it now reads it back from here.
    by_pos = {}
    for col in dos.columns:
        gene, position, _aa = col.split(':', 2)
        by_pos.setdefault(f'{gene}:{position}', []).append(col)
    undet = {k: float(1 - (dos[v].sum(axis=1) == 2).mean()) for k, v in by_pos.items()}
    pos_ctrl = {k: float(dos.loc[ctrl, v].sum().sum()) for k, v in by_pos.items()}

    rows = []
    for col in dos.columns:
        gene, position, aa = col.split(':', 2)
        key = f'{gene}:{position}'
        nc = float(dos.loc[ctrl, col].sum())
        na = float(dos.loc[case, col].sum())
        d = pos_ctrl[key]
        n_al = carried.get(col, 0)
        valid = 0 if gene in no_dosage else 1
        note = []
        if not valid:
            note.append('hemizygous written as two copies; do not use this dosage')
        if n_al == 1:
            note.append('carried by a single allele in this cohort')
        if undet[key] > 0.05:
            note.append('position undetermined on more than 5 % of samples')
        rows.append({
            'gene': gene,
            'position': int(position),
            'aa': aa,
            'residue': col,
            'n_chr_ctrl': int(nc),
            'n_chr_case': int(na),
            'n_chr_total': int(nc + na),
            'n_chr_ctrl_position': int(d),
            'freq_ctrl': round(nc / d, 6) if d else '',
            'n_alleles_carrying': n_al,
            'position_undetermined_rate': round(undet[key], 6),
            'gene_dosage_valid': valid,
            'note': '; '.join(note),
        })
    out = pd.DataFrame(rows).sort_values(['gene', 'position', 'aa'])
    out.to_csv(a.out_residue, sep='\t', index=False)
    print(f'[quality_tables] {len(out):,} residue row(s) -> {a.out_residue}')
    return out


def main():
    a = parse_args()
    al = allele_table(a)
    rs = residue_table(a)
    print(f'[quality_tables] {int((al.in_reference == "0").sum()):,} allele(s) that no '
          f'panel lists; {int((rs.n_alleles_carrying == 1).sum()):,} residue(s) carried '
          f'by one allele')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:                       # noqa: BLE001
        print(f'Error in quality_tables: {e}', file=sys.stderr)
        sys.exit(1)
