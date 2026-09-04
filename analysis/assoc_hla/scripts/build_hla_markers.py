#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : hla.typing's dosage matrices -> a VCF that plink2 turns into the
#           .bed/.bim/.fam both association engines read.
#
#           WHY A VCF AND NOT A HAND-WRITTEN .bed. The allele order decides
#           which allele the effect is reported for, and getting it wrong
#           inverts every odds ratio silently. A VCF states REF and ALT
#           explicitly and plink2 carries them through, so the effect allele is
#           declared rather than inferred from a bit-packing convention. The
#           .bed is then plink2's own, not ours.
#
#           ENCODING. One biallelic marker per HLA allele and per amino-acid
#           residue, REF = 'A' (absent), ALT = 'P' (present) -- the SNP2HLA
#           convention. Dosage is an exact integer count of that allele or
#           residue on the two chromosomes, so a hard call loses nothing:
#             2 -> 1/1   1 -> 0/1   0 -> 0/0   no call -> ./.
#           ALT = P is therefore the effect allele, and a positive log-odds
#           means carrying the HLA allele raises risk.
#
#           POSITIONS ARE A PLOTTING CONVENTION, NOT A GENOMIC CLAIM. Each
#           marker sits at its gene's GRCh38 start plus its index within the
#           gene. That puts every marker inside the right gene, in a stable
#           order, on the real MHC axis -- which is what the regional figures
#           need. It does NOT assert that the allele is caused by a variant at
#           that base. Positions being genuinely on chr6 does buy one real
#           thing: SAIGE's --LOCO=TRUE holds chr6 out of the null these markers
#           are tested against, so the GRM cannot contain the signal.
# Component: assoc_hla
# ---------------------------------------------------------------------------
import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

REF_ABSENT, ALT_PRESENT = 'A', 'P'


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--allele-dosage', required=True, help='04.residues/allele_dosage.tsv')
    p.add_argument('--residue-dosage', required=True, help='04.residues/residue_dosage.tsv')
    p.add_argument('--gene-coords', required=True, help='hla_gene_coords.R output')
    p.add_argument('--fam', required=True,
                   help="the cohort's .fam — defines WHICH samples and IN WHAT ORDER")
    p.add_argument('--pgroup-map', default=None,
                   help='05.qc/allele_pgroup_map.tsv, for the in_reference flag')
    p.add_argument('--residue-reference', default=None,
                   help='04.residues/residue_reference/ — needed to carry the allele-level '
                        'in_reference flag down to the residue markers the alleles carry')
    p.add_argument('--drop-genes', default='DRB3,DRB4,DRB5',
                   help='genes held out of the primary analysis. See docs/METHODS.md: their '
                        'dosages are known-wrong (hemizygotes encoded as homozygotes), DRB3 '
                        'fails the jMorp frequency check, and DRB3/DRB4 have no GRCh38 '
                        'primary-assembly locus.')
    p.add_argument('--keep-only-genes', default=None,
                   help='inverse of --drop-genes, for the sensitivity arm')
    p.add_argument('--cohort', required=True)
    p.add_argument('--out-vcf', default='hla_markers.vcf')
    p.add_argument('--out-map', default='marker_map.tsv')
    return p.parse_args()


def read_dosage(path):
    """sample_id MUST be a string: 107 of the 3,569 ids carry leading zeros and
    pandas silently turns 0000063134 into 63134, losing them on the join."""
    d = pd.read_csv(path, sep='\t', dtype={'sample_id': str})
    if d.columns[0] != 'sample_id':
        raise SystemExit(f'ABORT: {path} first column is {d.columns[0]!r}, not sample_id')
    return d


def gene_of(col, kind):
    return col.split('*', 1)[0] if kind == 'allele' else col.split(':', 1)[0]


def marker_id(col, kind):
    if kind == 'allele':
        return f'HLA_{col}'                      # HLA_A*24:02
    gene, pos, aa = col.split(':', 2)
    return f'AA_{gene}_{pos}_{aa}'               # AA_A_9_F, AA_A_-22_I


def sort_key(col, kind):
    """Stable within-gene order, so a marker's position never moves between runs."""
    if kind == 'allele':
        return (0, col)
    gene, pos, aa = col.split(':', 2)
    return (1, int(pos), aa)


def main():
    args = parse_args()

    fam = pd.read_csv(args.fam, sep=r'\s+', header=None, dtype=str,
                      names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'])
    order = list(fam.IID)
    coords = pd.read_csv(args.gene_coords, sep='\t', dtype={'chrom': str})
    cmap = dict(zip(coords.gene, coords.gene_start))

    drop = {g.strip() for g in args.drop_genes.split(',') if g.strip()}
    keep = ({g.strip() for g in args.keep_only_genes.split(',') if g.strip()}
            if args.keep_only_genes else None)

    frames, meta = {}, []
    for path, kind in ((args.allele_dosage, 'allele'), (args.residue_dosage, 'residue')):
        d = read_dosage(path).set_index('sample_id')
        missing = [s for s in order if s not in d.index]
        if missing:
            raise SystemExit(f'ABORT: {len(missing)} cohort sample(s) absent from {path}, '
                             f'first: {missing[:3]}. The cohorts are 100 % HLA-typed, so '
                             f'this means an id-type mismatch, not a real gap.')
        d = d.reindex(order)
        cols = []
        for c in d.columns:
            g = gene_of(c, kind)
            if g in drop or (keep is not None and g not in keep):
                continue
            if g not in cmap:                    # DRB3/DRB4: no primary-assembly locus
                continue
            cols.append(c)
        frames[kind] = d[cols].apply(pd.to_numeric, errors='coerce')
        meta += [(c, kind, gene_of(c, kind)) for c in cols]

    # AN UNDETERMINED RESIDUE IS NOT AN ABSENT ONE.
    # IMGT writes '*' for an unsequenced residue and '.' for a gap, so a
    # chromosome can carry no determined residue at a position. Its dosage is 0
    # at every residue there -- indistinguishable, in the matrix, from carrying a
    # different residue. That is wrong, and it is wrong in the worst possible
    # direction: typing quality in this cohort is DIFFERENTIAL between cases and
    # controls (hla.typing OPEN_QUESTIONS §3), so a residue whose complement is
    # partly "unknown" picks up the typing artefact as an association.
    #
    # Measured on the first run, before this: HLA-A position 175 is 99.67 % G, and
    # of the 21 samples not homozygous G, 5 -- 24 % -- were undetermined rather
    # than carrying A or R. The marker came out as the strongest residue signal in
    # the study, at P = 2.5e-8, on twenty chromosomes.
    #
    # So a sample whose residue dosages at a position do not sum to 2 has at least
    # one chromosome of unknown residue, and is set MISSING at every residue of
    # that position. It is not that we know it lacks them; we do not know.
    rframes = frames['residue']
    determined = {}
    if len(rframes.columns):
        by_pos = {}
        for c in rframes.columns:
            by_pos.setdefault(c.rsplit(':', 1)[0], []).append(c)
        n_masked = 0
        for pos_id, cols in by_pos.items():
            s = rframes[cols].sum(axis=1, skipna=False)
            seen = s[s.notna()]
            # Recorded BEFORE the mask. Afterwards every retained sample sums to 2
            # by construction, so measuring it then would report 1.0 everywhere and
            # hide exactly the data-quality fact this column exists to carry.
            determined[pos_id] = round(float((seen == 2).mean()), 6) if len(seen) else np.nan
            bad = s.notna() & (s != 2)
            if bad.any():
                rframes.loc[bad, cols] = np.nan
                n_masked += int(bad.sum())
        frames['residue'] = rframes
        print(f'[build_hla_markers] masked {n_masked:,} (sample, position) pair(s) whose '
              f'residues do not sum to 2 -- at least one chromosome undetermined')

    m = pd.DataFrame(meta, columns=['key', 'marker_class', 'gene'])
    m['_sort'] = [sort_key(k, c) for k, c in zip(m.key, m.marker_class)]
    m = m.sort_values(['gene', '_sort']).drop(columns='_sort').reset_index(drop=True)
    m['pos'] = m.groupby('gene').cumcount() + m.gene.map(cmap)
    m['id'] = [marker_id(k, c) for k, c in zip(m.key, m.marker_class)]
    m['chrom'] = '6'
    m = m.sort_values(['pos']).reset_index(drop=True)
    if m.id.duplicated().any():
        dup = m.id[m.id.duplicated()].tolist()[:5]
        raise SystemExit(f'ABORT: duplicate marker id(s): {dup}')

    # AC/AN over the cohort, for the QC stage and for the map.
    D = pd.concat([frames['allele'], frames['residue']], axis=1)[list(m.key)]
    arr = D.to_numpy(dtype=float)
    ac = np.nansum(arr, axis=0)
    an = 2.0 * np.sum(~np.isnan(arr), axis=0)
    m['ac'] = ac.astype(int)
    m['an'] = an.astype(int)
    with np.errstate(invalid='ignore', divide='ignore'):
        af = np.where(an > 0, ac / an, np.nan)
    m['af'] = np.round(af, 6)
    m['mac'] = np.minimum(ac, an - ac).astype(int)
    m['maf'] = np.round(np.minimum(af, 1 - af), 6)
    m['call_rate'] = np.round(an / (2.0 * len(order)), 6)

    # PER-POSITION RESIDUE DETERMINATION. IMGT's alignment writes '*' for an
    # unsequenced residue and '.' for a gap, so an allele carrying either
    # contributes NO residue at that position and the position's dosages sum to 1
    # or 0 rather than 2. A residue dosage is therefore a count among chromosomes
    # with a DETERMINED residue, not among all typed chromosomes. Median
    # undetermined rate is 0.03 %, but 157 of 1,115 positions exceed 5 % and
    # DQA1:56 reaches 24.8 % -- there the m-1 df omnibus is testing a position
    # a quarter of whose chromosomes have no residue to contribute, so
    # hla_marker_qc.py holds those out of it.
    m['position'] = np.where(m.marker_class == 'residue',
                             m.key.str.rsplit(':', n=1).str[0], '')
    m['determined_rate'] = m.position.map(determined)
    # RESIDUE MARKERS INHERIT THE ALLELE ARTEFACT WITHOUT INHERITING THE FLAG.
    # in_reference marks an allele the 61,424-person Japanese reference panel has
    # never observed -- the strongest typing-artefact warning available. A residue
    # marker carries no such name, but it is CARRIED BY alleles, and if the
    # chromosomes contributing it are mostly unconfirmed alleles then the residue
    # signal is the allele artefact wearing different clothes.
    #
    # Measured, and this is why the column exists: before it, the strongest
    # residue signal in the study was AA_A_178_T at P = 8e-9, and 15 of the 20
    # chromosomes driving it carried A*02:783 -- one allele, in_reference = 0.
    # The allele-level flag would have caught it; the residue layer had no flag
    # to catch it with.
    #
    # frac_unconfirmed_allele is the share of the marker's MINOR-state
    # chromosomes whose allele is absent from the reference panel -- the minor
    # state because that is what an association at low count is made of. It is a
    # WARNING, not a filter: a real Japanese allele the panel never sampled lands
    # here too.
    if args.pgroup_map and args.residue_reference and len(frames['residue'].columns):
        pgm = pd.read_csv(args.pgroup_map, sep='\t')
        unconf = set(pgm.loc[pgm.in_reference == 0, 'allele'])
        ad_raw = read_dosage(args.allele_dosage).set_index('sample_id').reindex(order)
        frac = {}
        for gene in sorted({c.split(':', 1)[0] for c in frames['residue'].columns}):
            ref_path = Path(args.residue_reference) / f'{gene}.tsv'
            if not ref_path.exists():
                continue
            ref = pd.read_csv(ref_path, sep='\t', dtype=str)
            ref['f2'] = [f"{gene}*" + ':'.join(a.split('*')[1].split(':')[:2])
                         for a in ref.allele]
            ref = ref.drop_duplicates('f2').set_index('f2')
            acols = [c for c in ad_raw.columns
                     if c.split('*', 1)[0] == gene and c in ref.index]
            if not acols:
                continue
            # Total chromosome count per allele across the cohort, not the
            # per-sample matrix: what is wanted is how many CHROMOSOMES carrying
            # this residue came from an unconfirmed allele.
            ac_per_allele = (ad_raw[acols].apply(pd.to_numeric, errors='coerce')
                             .fillna(0).to_numpy().sum(axis=0))
            bad = np.array([c in unconf for c in acols], dtype=float)
            for col in [c for c in frames['residue'].columns if c.split(':', 1)[0] == gene]:
                _g, pos, aa = col.split(':', 2)
                if pos not in ref.columns:
                    continue
                col_ref = ref.loc[acols, pos].to_numpy()
                carries = (col_ref == aa).astype(float)
                # An association at low count is driven by the MINOR state, so
                # that is the side to interrogate. AA_A_178_T is 99.6 % T: asking
                # what carries T says nothing, while asking what the twenty
                # non-T chromosomes are is the whole question -- fifteen of them
                # were A*02:783, one allele the reference panel has never seen.
                # 'determined' excludes '*' and '.', which are neither state.
                other = ((col_ref != aa) & (col_ref != '*') & (col_ref != '.')).astype(float)
                side = carries if float(ac_per_allele @ carries) <= float(ac_per_allele @ other) \
                    else other
                tot = float(ac_per_allele @ side)
                frac[col] = (round(float(ac_per_allele @ (side * bad)) / tot, 6)
                             if tot > 0 else np.nan)
        m_frac = frac
    else:
        m_frac = {}
    m['frac_unconfirmed_allele'] = m.key.map(m_frac)


    if args.pgroup_map:
        pg = pd.read_csv(args.pgroup_map, sep='\t')
        pg['key'] = pg.allele
        m = m.merge(pg[['key', 'p_group', 'in_reference', 'freq_reference']],
                    on='key', how='left')
    for c in ('p_group', 'in_reference', 'freq_reference'):
        if c not in m.columns:
            m[c] = ''
    # Blank, not 0, for a residue marker: the reference panel is an ALLELE panel,
    # so "not in it" is not a statement about a residue. Int64 keeps the allele
    # flag printing as 0/1 rather than 0.0/1.0.
    m['in_reference'] = pd.to_numeric(m.in_reference, errors='coerce').astype('Int64')
    m.loc[m.marker_class != 'allele', 'in_reference'] = pd.NA
    m['cohort'] = args.cohort

    # ---- VCF -------------------------------------------------------------
    gt = np.full(arr.shape, './.', dtype=object)
    gt[arr == 0] = '0/0'
    gt[arr == 1] = '0/1'
    gt[arr == 2] = '1/1'
    other = (~np.isnan(arr)) & ~np.isin(arr, [0, 1, 2])
    if other.any():
        raise SystemExit(f'ABORT: {int(other.sum())} dosage value(s) are not 0, 1, 2 or '
                         f'missing. These matrices are exact chromosome counts; anything '
                         f'else means the wrong file was passed.')
    gt = pd.DataFrame(gt, columns=list(D.columns))[list(m.key)].to_numpy()

    with open(args.out_vcf, 'w') as fh:
        fh.write('##fileformat=VCFv4.2\n')
        fh.write('##source=assoc_hla/build_hla_markers.py\n')
        fh.write('##contig=<ID=6>\n')
        fh.write(f'##ALT_meaning=P: the HLA allele or residue is present on that chromosome\n')
        fh.write(f'##REF_meaning=A: it is absent\n')
        fh.write('##INFO=<ID=HLACLASS,Number=1,Type=String,Description="allele or residue">\n')
        fh.write('##INFO=<ID=HLAGENE,Number=1,Type=String,Description="HLA gene">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
                 + '\t'.join(order) + '\n')
        for j, r in enumerate(m.itertuples()):
            fh.write(f'6\t{r.pos}\t{r.id}\t{REF_ABSENT}\t{ALT_PRESENT}\t.\t.\t'
                     f'HLACLASS={r.marker_class};HLAGENE={r.gene}\tGT\t'
                     + '\t'.join(gt[:, j]) + '\n')

    m.drop(columns=['key']).assign(key=m.key).to_csv(args.out_map, sep='\t', index=False)

    n_a = int((m.marker_class == 'allele').sum())
    n_r = int((m.marker_class == 'residue').sum())
    print(f'[build_hla_markers] {args.cohort}: {len(order):,} samples x {len(m):,} markers '
          f'({n_a:,} allele, {n_r:,} residue) over {m.gene.nunique()} gene(s) '
          f'-> {args.out_vcf}')
    print(f'    genes dropped by request: {", ".join(sorted(drop)) or "none"}')
    absent = sorted(set(coords.gene) - set(m.gene))
    if absent:
        print(f'    genes with coordinates but no surviving marker: {", ".join(absent)}')
    print(f'    call rate min {m.call_rate.min():.4f}, median {m.call_rate.median():.4f}')
    dr = m.loc[m.marker_class == 'residue', ['position', 'determined_rate']].drop_duplicates()
    if len(dr):
        print(f'    residue positions {len(dr):,}; determination rate median '
              f'{dr.determined_rate.median():.4f}, {(dr.determined_rate < 0.95).sum()} below 0.95')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in build_hla_markers: {e}', file=sys.stderr)
        sys.exit(1)
