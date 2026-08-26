#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Allele calls + the IMGT residue reference -> the genotype matrices an
#           association reads, at both levels it is normally run at.
#
#           allele_dosage.tsv     : one row per sample, one column per 2-field
#               allele, value 0/1/2. The classical allele test, which is the first
#               thing every HLA paper reports. It needs no reference — it is a
#               re-encoding of allele_calls.tsv — but it is produced HERE so that it
#               shares the residue matrices' sample order exactly.
#
#           residue_diplotype.tsv : one row per sample, one column per
#               (gene, IMGT position), value the two residues sorted. Human
#               readable, and SORTED so a heterozygote is one value and not two:
#               concatenating in call order makes 'LX' and 'XL' different
#               genotypes, which splits every heterozygote across two categories.
#
#           residue_dosage.tsv    : one row per sample, one column per
#               (gene, position, residue), value 0/1/2. This is what a regression
#               takes; it is order-free by construction.
#
#           How each allele was matched is recorded ONCE PER ALLELE, in
#           allele_match_map.tsv, because that is what it is a property of. Recording
#           it per (sample, allele) — as this did — makes a table of ~120,000 rows for
#           this cohort, 99 % of them successful prefix matches, in which the handful
#           of alleles with no match at all is invisible. allele_unmatched.tsv now
#           holds only those.
#
#           An allele is matched to the reference at the deepest field depth that
#           resolves it, then by prefix. HLA-HD reports 2-, 3- and 4-field names
#           and IMGT indexes all of them, so truncating the REFERENCE (as one
#           obvious implementation does) merges null alleles — `A*03:01:01:02N`,
#           a non-expressed protein — into the functional allele beside them.
#           Matching by prefix keeps the distinction and reports what was used.
# Component: hla.typing
# ---------------------------------------------------------------------------
import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

MISSING = '.'


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--calls', required=True, help='allele_calls.tsv')
    p.add_argument('--reference-dir', required=True, help='residue_reference/')
    p.add_argument('--genes', required=True, help='comma-separated')
    p.add_argument('--out-diplotype', default='residue_diplotype.tsv')
    p.add_argument('--out-dosage', default='residue_dosage.tsv')
    p.add_argument('--out-sites', default='residue_sites.tsv')
    p.add_argument('--out-unmatched', default='allele_unmatched.tsv')
    p.add_argument('--out-match-map', default='allele_match_map.tsv')
    p.add_argument('--out-allele-dosage', default='allele_dosage.tsv')
    p.add_argument('--allele-field-depth', type=int, default=2,
                   help='field depth for the allele dosage columns. 2 is the '
                        'resolution HLA association is conventionally reported at, '
                        'and the depth every call reaches.')
    p.add_argument('--min-count', type=int, default=1,
                   help='drop a residue seen on fewer than this many chromosomes '
                        '(0 keeps everything; the default keeps anything observed)')
    return p.parse_args()


def build_index(ref):
    """{allele: row} plus a prefix index for shallower calls.

    IMGT names nest: A*24:02 is a prefix of A*24:02:01:01. A call at 2 fields is
    matched to the alphabetically first 4-field allele under it, which is IMGT's
    own representative convention and is deterministic.
    """
    exact = {a: i for i, a in enumerate(ref.index)}
    by_prefix = defaultdict(list)
    for a in ref.index:
        f = a.split(':')
        for d in range(1, len(f) + 1):
            by_prefix[':'.join(f[:d])].append(a)
    return exact, by_prefix


def resolve(allele, exact, by_prefix):
    """(matched_allele, how) or (None, reason)."""
    if allele in exact:
        return allele, 'exact'
    cands = by_prefix.get(allele)
    if cands:
        return sorted(cands)[0], f'prefix({len(cands)})'
    # A deeper call than the reference carries: walk back a field at a time.
    f = allele.split(':')
    for d in range(len(f) - 1, 0, -1):
        key = ':'.join(f[:d])
        if key in exact:
            return key, f'truncated_to_{d}'
        if key in by_prefix:
            return sorted(by_prefix[key])[0], f'truncated_to_{d}'
    return None, 'no_match'


def allele_dosage(calls, genes, depth, min_count):
    """{allele: [0/1/2 per sample]} at `depth` fields, over every gene in `genes`.

    A hemizygous call is two copies of one allele, which is what collect_alleles.py
    writes and is right for the classical loci and wrong for DRB3/4/5 — see
    docs/OPEN_QUESTIONS.md §1. The defect is the same one the residue matrix carries;
    it is not introduced here and it is not corrected here either.
    """
    out = {}
    for g in genes:
        c1, c2, cs = f'{g}_1', f'{g}_2', f'{g}_state'
        if c1 not in calls.columns:
            continue
        per_sample = []
        for _, r in calls.iterrows():
            if r.get(cs, '') in ('failed', 'not_typed') or not r[c1] or not r[c2]:
                per_sample.append(None)
            else:
                per_sample.append((trim(r[c1], depth), trim(r[c2], depth)))
        seen = sorted({a for v in per_sample if v for a in v})
        for al in seen:
            col = [(v.count(al) if v else '') for v in per_sample]
            if sum(x for x in col if x != '') >= min_count:
                out[al] = col
    return out


def trim(allele, depth):
    """`A*24:02:01:02L` at depth 2 -> `A*24:02`. Shallower calls are left alone."""
    f = allele.split(':')
    return ':'.join(f[:depth]) if len(f) > depth else allele


def main():
    args = parse_args()
    calls = pd.read_csv(args.calls, sep='\t', dtype=str).fillna('')
    if 'sample_id' not in calls.columns:
        raise SystemExit(f'ABORT: {args.calls} has no sample_id column')
    genes = [g.strip() for g in args.genes.split(',') if g.strip()]
    samples = calls['sample_id'].tolist()

    diplo, dosage, sites, unmatched = {}, {}, [], []
    # (gene, allele) -> [matched_to, how, n_chromosomes]. One row per ALLELE, not
    # per (sample, allele): how an allele matches the reference is a property of the
    # allele.
    match_map = {}

    for g in genes:
        ref_path = Path(args.reference_dir) / f'{g}.tsv'
        if not ref_path.exists():
            # No IMGT protein alignment for this gene. Recorded here and in the
            # sites table; never silently skipped.
            sites.append({'gene': g, 'positions': 0, 'polymorphic': 0,
                          'note': 'no IMGT protein alignment'})
            continue
        ref = pd.read_csv(ref_path, sep='\t', index_col=0, dtype=str)
        exact, by_prefix = build_index(ref)
        c1, c2, cs = f'{g}_1', f'{g}_2', f'{g}_state'
        if c1 not in calls.columns:
            continue

        # sample -> (row1, row2) of residues, or None where the call is unusable
        per_sample = {}
        for _, r in calls.iterrows():
            sid, a1, a2, st = r['sample_id'], r[c1], r[c2], r.get(cs, '')
            if st in ('failed', 'not_typed') or not a1 or not a2:
                per_sample[sid] = None
                continue
            rows = []
            for a in (a1, a2):
                key = (g, a)
                if key in match_map:
                    match_map[key][2] += 1
                    m, how = match_map[key][0], match_map[key][1]
                else:
                    m, how = resolve(a, exact, by_prefix)
                    match_map[key] = [m or '', how, 1]
                if m is None:
                    # No match at ANY depth. This is the only thing worth a per-sample
                    # row, because it is the only one that costs a genotype.
                    unmatched.append({'sample_id': sid, 'gene': g, 'allele': a,
                                      'reason': how})
                    rows.append(None)
                else:
                    rows.append(ref.loc[m])
            per_sample[sid] = rows if all(x is not None for x in rows) else None

        n_poly = 0
        for pos in ref.columns:
            vals = {}
            for sid in samples:
                rows = per_sample.get(sid)
                if rows is None:
                    vals[sid] = None
                    continue
                vals[sid] = (rows[0][pos], rows[1][pos])
            observed = {a for v in vals.values() if v for a in v} - {MISSING, '*'}
            if len(observed) < 2:
                continue                     # monomorphic here; nothing to test
            n_poly += 1
            key = f'{g}:{pos}'
            diplo[key] = [(''.join(sorted(vals[s])) if vals[s] else '') for s in samples]
            for aa in sorted(observed):
                col = [(vals[s].count(aa) if vals[s] else '') for s in samples]
                n_chr = sum(x for x in col if x != '')
                if n_chr >= args.min_count:
                    dosage[f'{g}:{pos}:{aa}'] = col
        sites.append({'gene': g, 'positions': len(ref.columns),
                      'polymorphic': n_poly, 'note': ''})
        print(f'[residue_matrix] {g:6s} {len(ref.columns):4d} positions, '
              f'{n_poly:4d} polymorphic in this cohort')

    # Every gene with calls, not only the 19 that carry residues: the classical
    # allele test is run on all of them.
    called_genes = [c[:-6] for c in calls.columns if c.endswith('_state')]
    ad = allele_dosage(calls, called_genes, args.allele_field_depth, args.min_count)

    def write(path, cols):
        df = pd.DataFrame(cols, index=samples) if cols else pd.DataFrame(index=samples)
        df.index.name = 'sample_id'
        df.to_csv(path, sep='\t')
        return df.shape

    s1 = write(args.out_diplotype, diplo)
    s2 = write(args.out_dosage, dosage)
    s3 = write(args.out_allele_dosage, ad)
    pd.DataFrame(sites).to_csv(args.out_sites, sep='\t', index=False)

    # Both tables get their header even when they have no rows: an empty file is
    # not readable by pandas, and 'no allele failed to match' is a result.
    pd.DataFrame(unmatched, columns=['sample_id', 'gene', 'allele', 'reason']) \
        .to_csv(args.out_unmatched, sep='\t', index=False)
    pd.DataFrame([{'gene': g, 'allele': a, 'matched_to': v[0], 'how': v[1],
                   'n_chromosomes': v[2]}
                  for (g, a), v in sorted(match_map.items())],
                 columns=['gene', 'allele', 'matched_to', 'how', 'n_chromosomes']) \
        .to_csv(args.out_match_map, sep='\t', index=False)

    print(f'[residue_matrix] allele dosage {s3[0]} x {s3[1]} -> {args.out_allele_dosage}'
          f'  ({args.allele_field_depth}-field, {len(called_genes)} genes)')
    print(f'[residue_matrix] diplotype    {s1[0]} x {s1[1]} -> {args.out_diplotype}')
    print(f'[residue_matrix] dosage       {s2[0]} x {s2[1]} -> {args.out_dosage}')
    n_exact = sum(1 for v in match_map.values() if v[1] == 'exact')
    print(f'    {len(match_map)} distinct (gene, allele) matched to the reference: '
          f'{n_exact} exact, {len(match_map) - n_exact} by prefix or truncation '
          f'-> {args.out_match_map}')
    if unmatched:
        print(f'    {len(unmatched)} call(s) with NO match at any depth '
              f'-> {args.out_unmatched}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in residue_matrix: {e}', file=sys.stderr)
        sys.exit(1)
