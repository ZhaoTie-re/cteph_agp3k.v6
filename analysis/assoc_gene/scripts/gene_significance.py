#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The decision rules for ONE (cohort, stratum, method), stated as a
#           three-row table beside the scan they judged.
#
#           Two rules are applied and both are named here so a reader never has
#           to infer them from the scan:
#
#             bonferroni  p < alpha / n_genes_mapped. n is the number of genes
#                         with >= MinNumVar variants in THIS cohort x stratum
#                         map, FIXED AT THE MAP BARRIER before any association
#                         ran. The number is the same for cmc and skato, so it is
#                         the one threshold the two methods can be compared on.
#             bh          Benjamini-Hochberg q < alpha over the SAME family; the
#                         genes rvtest did not return enter at p = 1. Its cut
#                         (the largest p rejected) depends on the method's own
#                         p-value distribution, so it is a different number per
#                         method and must not be compared across them.
#             nominal     alpha itself. Reference only; never a decision rule,
#                         and no figure in this component draws a line at it.
#
#           WHY THE THRESHOLD IS RECOMPUTED AND THEN CHECKED, NOT COPIED. The
#           value is alpha / n_genes_mapped, and both files this step reads
#           already carry it: the denominator wrote it and the scan judged every
#           gene by it. A corrupted denominator must not propagate, so the number
#           is derived again here -- and then REQUIRED to match both files, so a
#           table can never state a rule the scan did not apply.
# Component: assoc_gene
# ---------------------------------------------------------------------------
import argparse
import sys

METHODS = ('cmc', 'skato')

# The scan_qc.tsv schema as scan.py writes it -- the only p-value summary this
# step sees. Read by name and asserted, so a moved column is an abort here and
# not a mis-labelled count in a published table.
QC_COLUMNS = ('cohort', 'stratum', 'method', 'n_genes_mapped', 'n_genes_returned',
              'n_genes_missing', 'lambda_gc', 'lambda_gc_nvar2', 'lambda_gc_nvar_ge3',
              'min_p', 'n_significant_bonferroni', 'n_significant_bh', 'bh_threshold',
              'bonferroni_threshold', 'min_num_var')

DEN_COLUMNS = ('cohort', 'stratum', 'min_num_var', 'n_genes_mapped',
               'n_genes_dropped_lt_minvar', 'n_variants', 'n_variants_multigene',
               'n_sets_chrom_split', 'n_positions_colliding', 'alpha', 'bonferroni_threshold')

# significance.<method>.<stratum>.tsv, pinned. `primary` says whether a hit can
# come from this row; `comparable_across_methods` says whether this row's value
# is the same number for cmc and skato. Bonferroni is both; BH is primary but
# not comparable; nominal is neither.
OUT_COLUMNS = ('threshold', 'value', 'basis', 'primary', 'comparable_across_methods',
               'cohort', 'stratum', 'method', 'n_genes_mapped', 'n_genes_returned',
               'n_genes_missing', 'min_num_var', 'lambda_gc', 'min_p', 'n_significant',
               'alpha')

NA = 'NA'
NA_TOKENS = {'', 'NA', 'na', 'NaN', 'nan', 'None', 'null'}
SIGDIG = '.12g'
# The tolerance is on MERGE_MAP's 12-significant-digit round trip, not on any
# physical quantity; anything looser would let a genuinely different denominator
# through.
REL_TOL = 1e-9


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scan-qc', required=True, help='SCAN scan_qc.tsv; one row per method')
    p.add_argument('--denominator', required=True,
                   help='denominator.<stratum>.tsv from MERGE_MAP; the map barrier, one row')
    p.add_argument('--cohort', required=True)
    p.add_argument('--stratum', required=True)
    p.add_argument('--method', required=True, choices=list(METHODS))
    p.add_argument('--alpha', type=float, default=0.05,
                   help='family-wise error rate; must agree with the denominator')
    p.add_argument('--out', default='significance.tsv')
    return p.parse_args()


def clean(v):
    v = v.strip()
    return NA if v in NA_TOKENS else v


def read_tsv(path, required, label):
    """[{col: value}], header asserted by name. Both inputs are a handful of rows."""
    with open(path) as fh:
        first = fh.readline()
        if not first.strip():
            raise SystemExit(f'ABORT: {path} is empty; {label} must carry a header and at '
                             f'least one data row')
        head = [c.lstrip('#') for c in first.rstrip('\n').split('\t')]
        missing = [c for c in required if c not in head]
        if missing:
            raise SystemExit(f'ABORT: {path} has no {missing} column(s); found {head}. '
                             f'This is not {label}.')
        rows = []
        for ln, raw in enumerate(fh, 2):
            line = raw.rstrip('\n')
            if not line:
                continue
            f = line.split('\t')
            if len(f) != len(head):
                raise SystemExit(f'ABORT: {path}:{ln} has {len(f)} field(s), expected '
                                 f'{len(head)}')
            rows.append(dict(zip(head, f)))
    if not rows:
        raise SystemExit(f'ABORT: {path} carries a header and no data row; {label} is empty')
    return rows


def as_float(v, path, col):
    v = clean(v)
    if v == NA:
        return None
    try:
        return float(v)
    except ValueError:
        raise SystemExit(f'ABORT: {path} {col} is {v!r}, not a number')


def as_int(v, path, col):
    v = clean(v)
    if v == NA:
        return None
    try:
        return int(v)
    except ValueError:
        raise SystemExit(f'ABORT: {path} {col} is {v!r}, not an integer')


def fmt(x):
    return NA if x is None else format(x, SIGDIG)


def agree(a, b):
    return a is not None and b is not None and abs(a - b) <= REL_TOL * max(abs(a), abs(b), 1e-300)


def pick_qc_row(rows, args, path):
    """The one scan_qc row for this cell. Anything else is a mis-wired DAG, not a
    result to be worked around: SIGNIFICANCE fans out over the two methods, so a
    silent miss would publish a file describing a scan that never ran."""
    same_cell = [r for r in rows if r['cohort'] == args.cohort and r['stratum'] == args.stratum]
    if not same_cell:
        seen = sorted({(r['cohort'], r['stratum']) for r in rows})
        raise SystemExit(f'ABORT: {path} carries no row for {args.cohort}/{args.stratum}; it '
                         f'holds {seen}. The wrong scan_qc.tsv was staged.')
    hits = [r for r in same_cell if r['method'] == args.method]
    if not hits:
        have = sorted({r['method'] for r in same_cell})
        raise SystemExit(f'ABORT: {path} has no row for method={args.method}. It reports '
                         f'{have}. SCAN returned nothing for this cell.')
    if len(hits) > 1:
        raise SystemExit(f'ABORT: {path} has {len(hits)} rows for method={args.method}; '
                         f'scan_qc is one row per method, so two scans have been concatenated '
                         f'and it is ambiguous which one this table would describe.')
    return hits[0]


def check_denominator(den, qc, args, den_path, qc_path):
    """The map barrier as MERGE_MAP fixed it, cross-checked against what SCAN was
    handed. These must be the same MERGE_MAP task or the "one denominator, both
    methods" claim in the basis string is simply false."""
    if den['cohort'] != args.cohort or den['stratum'] != args.stratum:
        raise SystemExit(f'ABORT: {den_path} is {den["cohort"]}/{den["stratum"]}, not '
                         f'{args.cohort}/{args.stratum}')
    n_mapped = as_int(den['n_genes_mapped'], den_path, 'n_genes_mapped')
    if n_mapped is None or n_mapped < 0:
        raise SystemExit(f'ABORT: {den_path} n_genes_mapped is {den["n_genes_mapped"]!r}')
    qc_mapped = as_int(qc['n_genes_mapped'], qc_path, 'n_genes_mapped')
    if qc_mapped != n_mapped:
        raise SystemExit(f'ABORT: n_genes_mapped is {n_mapped} in {den_path} but {qc_mapped} '
                         f'in {qc_path}. SCAN and this step were handed denominators from '
                         f'different MERGE_MAP tasks, so the two methods were NOT judged '
                         f'against one denominator.')
    min_num_var = as_int(den['min_num_var'], den_path, 'min_num_var')
    qc_minvar = as_int(qc['min_num_var'], qc_path, 'min_num_var')
    if min_num_var is None or qc_minvar != min_num_var:
        raise SystemExit(f'ABORT: min_num_var is {den["min_num_var"]!r} in {den_path} but '
                         f'{qc["min_num_var"]!r} in {qc_path}; the two sides counted a '
                         f'different gene set.')
    den_alpha = as_float(den['alpha'], den_path, 'alpha')
    if not agree(den_alpha, args.alpha):
        raise SystemExit(f'ABORT: --alpha {args.alpha:g} but {den_path} fixed the threshold at '
                         f'alpha={den["alpha"]!r}. The basis string and the threshold the scan '
                         f'was actually judged at would disagree.')
    return n_mapped, min_num_var


def main():
    args = parse_args()
    if args.alpha <= 0 or args.alpha >= 1:
        raise SystemExit(f'ABORT: --alpha must lie in (0, 1); got {args.alpha}')

    qc_rows = read_tsv(args.scan_qc, QC_COLUMNS, 'scan_qc.tsv')
    den_rows = read_tsv(args.denominator, DEN_COLUMNS, 'denominator.<stratum>.tsv')
    if len(den_rows) != 1:
        raise SystemExit(f'ABORT: {args.denominator} has {len(den_rows)} data rows; the '
                         f'denominator is exactly one row per stratum, so more than one means '
                         f'two strata were concatenated and the threshold is ambiguous.')
    den = den_rows[0]
    qc = pick_qc_row(qc_rows, args, args.scan_qc)
    n_mapped, min_num_var = check_denominator(den, qc, args, args.denominator, args.scan_qc)

    value = args.alpha / n_mapped if n_mapped else None
    if value is None:
        print(f'[gene_significance] WARNING: {args.stratum} retained no gene at '
              f'--min-num-var {min_num_var}; the bonferroni row is NA and nothing in this '
              f'cell can be called significant', file=sys.stderr)

    # Derived again, then required to match both files that already used it --
    # SCAN set gene_scan.significant_bonferroni from exactly this number.
    for path, src in ((args.denominator, den), (args.scan_qc, qc)):
        got = as_float(src['bonferroni_threshold'], path, 'bonferroni_threshold')
        if value is None:
            if got is not None:
                raise SystemExit(f'ABORT: {path} carries bonferroni_threshold={got:g} while '
                                 f'n_genes_mapped is 0; alpha / 0 is not a threshold.')
        elif not agree(got, value):
            raise SystemExit(f'ABORT: alpha / n_genes_mapped = {fmt(value)} but {path} carries '
                             f'bonferroni_threshold={src["bonferroni_threshold"]!r}. The scan '
                             f'was judged at that number; publishing this one would state a '
                             f'decision rule that was never applied.')

    alpha_s = format(args.alpha, 'g')
    n_missing_s = clean(qc['n_genes_missing'])
    rows = [
        {'threshold': 'bonferroni', 'value': fmt(value), 'primary': 1,
         'comparable_across_methods': 1,
         'n_significant': clean(qc['n_significant_bonferroni']),
         'basis': (f'{alpha_s} / {n_mapped} genes with >= {min_num_var} variants in this '
                   f'cohort x stratum map; fixed before any p-value; identical for cmc and '
                   f'skato')},
        # `value` is the largest p BH rejected, NA when it rejected nothing; it is
        # SCAN's number, copied verbatim, because no p-vector reaches this step.
        {'threshold': 'bh', 'value': clean(qc['bh_threshold']), 'primary': 1,
         'comparable_across_methods': 0,
         'n_significant': clean(qc['n_significant_bh']),
         'basis': (f'Benjamini-Hochberg q < {alpha_s} over the same {n_mapped} gene(s) '
                   f'(>= {min_num_var} variants); the {n_missing_s} gene(s) rvtest did not '
                   f'return enter at p = 1; the cut is method-specific, do NOT compare it '
                   f'between cmc and skato')},
        {'threshold': 'nominal', 'value': fmt(args.alpha), 'primary': 0,
         'comparable_across_methods': 1, 'n_significant': NA,
         'basis': 'uncorrected, for reference only; never a decision rule'},
    ]
    # lambda_gc and min_p describe the scan, not the row, so they are SCAN's
    # carried verbatim on all three. n_significant is per row: a `bh` row
    # silently carrying the Bonferroni count would state a rule that was never
    # applied, and `nominal` stays NA because a count under a non-rule would be
    # read as one.
    common = {
        'cohort': args.cohort, 'stratum': args.stratum, 'method': args.method,
        'n_genes_mapped': n_mapped,
        'n_genes_returned': clean(qc['n_genes_returned']),
        'n_genes_missing': clean(qc['n_genes_missing']),
        'min_num_var': min_num_var,
        'lambda_gc': clean(qc['lambda_gc']),
        'min_p': clean(qc['min_p']),
        'alpha': alpha_s,
    }
    with open(args.out, 'w') as fh:
        fh.write('\t'.join(OUT_COLUMNS) + '\n')
        for r in rows:
            assert not (set(common) & set(r)), sorted(set(common) & set(r))
            r.update(common)
            fh.write('\t'.join(str(r[c]) for c in OUT_COLUMNS) + '\n')

    cell = f'{args.cohort}/{args.stratum}/{args.method}'
    print(f'[gene_significance] {cell} -> {args.out}')
    print(f'    bonferroni  p < {fmt(value)}  = {alpha_s} / {n_mapped:,} gene(s) with >= '
          f'{min_num_var} mapped variant(s)  <- decision rule, same for both methods')
    print(f'    bh          q < {alpha_s} over the same {n_mapped:,} gene(s), cut at p = '
          f'{clean(qc["bh_threshold"])}  <- decision rule, THIS method only')
    print(f'    nominal     p < {fmt(args.alpha)}  reference only, never a decision rule')
    missing = as_int(qc['n_genes_missing'], args.scan_qc, 'n_genes_missing')
    returned = as_int(qc['n_genes_returned'], args.scan_qc, 'n_genes_returned')
    if missing:
        print(f'    {returned:,} of {n_mapped:,} mapped gene(s) returned; the {missing:,} '
              f'missing still cost a test')
    print(f'    lambda_GC = {clean(qc["lambda_gc"])} (n_var=2: {clean(qc["lambda_gc_nvar2"])}, '
          f'>=3: {clean(qc["lambda_gc_nvar_ge3"])}), min p = {clean(qc["min_p"])}, '
          f'{clean(qc["n_significant_bonferroni"])} gene(s) below the bonferroni threshold, '
          f'{clean(qc["n_significant_bh"])} at BH q < {alpha_s} '
          f'(the first set is contained in the second)')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in gene_significance: {e}', file=sys.stderr)
        sys.exit(1)
