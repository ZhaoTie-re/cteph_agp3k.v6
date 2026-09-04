#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Decide which HLA markers are tested, and record why every one that
#           is not was dropped. Emits one --extract list per marker class, so
#           the allele scan and the residue scan are separate analyses with
#           separate multiple-testing burdens rather than one mixed pile.
#
#           THE FREQUENCY FLOOR IS A COUNT, NOT A FREQUENCY. Nat. Protoc. 2023
#           recommends MAF >= 1 %, but its reason is imputation noise and these
#           calls are typed directly, so that reason does not apply here. What
#           does apply is the case count: 419-439. MAF 1 % is ~44 chromosomes in
#           the narrow cohort and about 8 on the case side, which is too thin to
#           estimate an odds ratio from. MAC >= 20 is the primary floor and
#           MAF >= 1 % is kept as a sensitivity arm, both reported.
#
#           RESIDUE DETERMINATION. IMGT writes '*' for an unsequenced residue and
#           '.' for a gap, so at some positions a chromosome contributes no
#           residue at all. Positions below --min-determined are still tested
#           marker by marker -- a residue either is or is not on a chromosome --
#           but are held OUT OF THE OMNIBUS, whose m-1 df assumes the position is
#           observed. See docs/METHODS.md.
# Component: assoc_hla
# ---------------------------------------------------------------------------
import argparse
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--marker-map', required=True, help='build_hla_markers.py output')
    p.add_argument('--cohort', required=True)
    p.add_argument('--min-mac', type=int, default=20,
                   help='primary floor: minor-allele COUNT')
    p.add_argument('--min-maf', type=float, default=0.0,
                   help='additional frequency floor; 0 disables. The MAF >= 0.01 '
                        'sensitivity arm sets this and leaves --min-mac at 0.')
    p.add_argument('--min-call-rate', type=float, default=0.95)
    p.add_argument('--min-determined', type=float, default=0.95,
                   help='per-position residue determination rate below which the position '
                        'is excluded from the OMNIBUS (not from single-marker tests)')
    p.add_argument('--arm', default='primary', help='label recorded in the output')
    p.add_argument('--out-qc', default='marker_qc.tsv')
    p.add_argument('--out-prefix', default='tested',
                   help='writes <prefix>.allele.txt and <prefix>.residue.txt for --extract')
    return p.parse_args()


def main():
    args = parse_args()
    m = pd.read_csv(args.marker_map, sep='\t')
    for c in ('id', 'marker_class', 'gene', 'mac', 'maf', 'call_rate'):
        if c not in m.columns:
            raise SystemExit(f'ABORT: {args.marker_map} has no {c!r} column')

    m['fail_mac'] = m.mac < args.min_mac
    m['fail_maf'] = (m.maf < args.min_maf) if args.min_maf > 0 else False
    m['fail_call_rate'] = m.call_rate < args.min_call_rate
    m['tested'] = ~(m.fail_mac | m.fail_maf | m.fail_call_rate)

    # Omnibus eligibility is a PROPERTY OF THE POSITION, not of one residue: a
    # position is testable only if at least two of its residues survive, because
    # the m-1 df test needs m >= 2, and only if the position is actually observed.
    det_ok = m.determined_rate.fillna(1.0) >= args.min_determined
    n_tested = (m[m.tested & (m.marker_class == 'residue')]
                .groupby('position').size().rename('n_residue_tested'))
    m = m.merge(n_tested, left_on='position', right_index=True, how='left')
    m['n_residue_tested'] = m.n_residue_tested.fillna(0).astype(int)
    m['omnibus_eligible'] = ((m.marker_class == 'residue') & m.tested
                             & det_ok & (m.n_residue_tested >= 2))

    m['drop_reason'] = ''
    m.loc[m.fail_call_rate, 'drop_reason'] = 'call_rate'
    m.loc[m.fail_maf & (m.drop_reason == ''), 'drop_reason'] = 'maf'
    m.loc[m.fail_mac & (m.drop_reason == ''), 'drop_reason'] = 'mac'
    m['arm'] = args.arm
    m['cohort'] = args.cohort
    m.to_csv(args.out_qc, sep='\t', index=False)

    for cls in ('allele', 'residue'):
        sel = m[(m.marker_class == cls) & m.tested]
        with open(f'{args.out_prefix}.{cls}.txt', 'w') as fh:
            fh.write('\n'.join(sel.id) + ('\n' if len(sel) else ''))

    n_a = int(((m.marker_class == 'allele') & m.tested).sum())
    n_r = int(((m.marker_class == 'residue') & m.tested).sum())
    n_p = int(m.loc[m.omnibus_eligible, 'position'].nunique())
    print(f'[hla_marker_qc] {args.cohort} / {args.arm}: MAC >= {args.min_mac}'
          + (f', MAF >= {args.min_maf}' if args.min_maf > 0 else '')
          + f', call rate >= {args.min_call_rate}')
    print(f'    allele   {n_a:5,} of {int((m.marker_class == "allele").sum()):5,} tested')
    print(f'    residue  {n_r:5,} of {int((m.marker_class == "residue").sum()):5,} tested')
    print(f'    omnibus  {n_p:5,} position(s) with >= 2 tested residues and determination '
          f'>= {args.min_determined}')
    drops = m.loc[~m.tested, 'drop_reason'].value_counts()
    if len(drops):
        print('    dropped: ' + ', '.join(f'{k} {v:,}' for k, v in drops.items()))
    if n_a == 0 or n_r == 0:
        raise SystemExit('ABORT: a marker class has no tested marker; the filters or the '
                         'marker map are wrong, and an empty scan would publish silently.')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in hla_marker_qc: {e}', file=sys.stderr)
        sys.exit(1)
