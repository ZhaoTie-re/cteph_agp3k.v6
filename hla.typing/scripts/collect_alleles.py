#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Every sample's HLA-HD result -> one allele table, one ambiguity table,
#           and one status table.
#
#           Samples are identified BY FILENAME. HLA-HD names its output
#           `<sample>_final.result.txt`, which is already unique, so the results are
#           read straight out of the directory they were staged into. The earlier
#           design passed `SAMPLE=res_N` pairs alongside files staged as `res_*`,
#           which made two independent orderings have to agree — and Nextflow does
#           not number `res_*` at all when the collection holds one file.
#
#           Parsed BY GENE NAME, never by row position. When pick_up_allele cannot
#           read an estimate it writes `Couldn't read result file.` on a line of its
#           own, carrying no gene name, so a positional parser silently attributes
#           every later gene to the wrong locus.
#
#           Four outcomes are kept apart, because collapsing them loses the one
#           distinction a reader needs:
#             called      two alleles
#             hemizygous  one allele and a '-' (HLA-HD's own convention)
#             not_typed   'Not typed' — the locus is absent from this sample
#             failed      the gene has no row at all, or its row is the unreadable
#                         token. This is a PIPELINE failure, not biology, and it is
#                         a differential-missingness confounder if it is not recorded.
#
#           Ambiguity is recorded, not silently resolved against the cohort. HLA-HD
#           emits every candidate when it cannot choose; resolving those by their
#           frequency IN THE CURRENT BATCH makes a sample's genotype depend on who
#           else was run with it, so the same sample genotypes differently in a
#           pilot and in the full cohort. Here the candidates are written out, and
#           the reported call is trimmed to the field depth they agree on.
# Component: hla.typing
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

FAIL_TOKEN = "Couldn't read result file."
NOT_TYPED = 'Not typed'


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--result-dir', default='.',
                   help='directory holding one <sample>_final.result.txt per sample')
    p.add_argument('--result-suffix', default='_final.result.txt',
                   help="HLA-HD's own naming; the sample id is the filename without it")
    p.add_argument('--expect-n', type=int, default=0,
                   help='abort unless exactly this many results are found. 0 disables. '
                        'Without it a partially staged directory produces a short table '
                        'and says nothing.')
    p.add_argument('--gene-split', required=True)
    p.add_argument('--out-calls', default='allele_calls.tsv')
    p.add_argument('--out-ambiguity', default='allele_ambiguity.tsv')
    p.add_argument('--out-status', default='typing_status.tsv')
    return p.parse_args()


def strip_prefix(a):
    return a[4:] if a.startswith('HLA-') else a


def trim(allele, depth):
    """`A*24:02:01:02L` at depth 2 -> `A*24:02`."""
    f = allele.split(':')
    return ':'.join(f[:depth]) if len(f) > depth else allele


def resolve_ambiguous(cands):
    """The deepest field depth at which every candidate PAIR is the same genotype.

    HLA-HD lists candidates as consecutive pairs, and it does not keep the two
    slots in a consistent order between them: G*01:04:03/G*01:01:08 and
    G*01:01:01/G*01:04:01 are the same genotype with the slots swapped. Comparing
    slot-wise therefore collapses to the genus (`G*01`) when the real answer is
    `G*01:01 / G*01:04`. Pairs are compared as UNORDERED sets instead, and the
    deepest depth at which they all agree is reported.

    Returns (allele_1, allele_2, depth) or (…, 0) when even one field disagrees.
    """
    pairs = [(cands[i], cands[i + 1]) for i in range(0, len(cands) - 1, 2)]
    if not pairs:
        return cands[0], cands[0], 0
    max_depth = max(len(a.split(':')) for a in cands)
    for depth in range(max_depth, 0, -1):
        trimmed = {frozenset((trim(a, depth), trim(b, depth))) for a, b in pairs}
        if len(trimmed) == 1:
            a, b = sorted(next(iter(trimmed)) if len(next(iter(trimmed))) == 2
                          else list(next(iter(trimmed))) * 2)
            return a, b, depth
    return trim(cands[0], 1), trim(cands[1], 1), 0


def parse_result(path):
    """{gene: [alleles]} plus the count of unreadable-token lines."""
    genes, n_fail = {}, 0
    for ln in Path(path).read_text(errors='replace').splitlines():
        ln = ln.rstrip('\n')
        if not ln.strip():
            continue
        if ln.strip() == FAIL_TOKEN:
            n_fail += 1                      # no gene name to attribute it to
            continue
        f = ln.split('\t')
        genes[f[0]] = [x for x in f[1:] if x.strip()]
    return genes, n_fail


def main():
    args = parse_args()
    want = [ln.split()[0] for ln in Path(args.gene_split).read_text().splitlines()
            if ln.strip() and not ln.startswith('#')]
    if not want:
        raise SystemExit(f'ABORT: no gene named in {args.gene_split}')

    rd = Path(args.result_dir)
    if not rd.is_dir():
        raise SystemExit(f'ABORT: --result-dir {rd} is not a directory')
    spec = {p.name[:-len(args.result_suffix)]: p
            for p in sorted(rd.glob(f'*{args.result_suffix}'))}
    if not spec:
        raise SystemExit(f'ABORT: no *{args.result_suffix} in {rd.resolve()}')
    if args.expect_n and len(spec) != args.expect_n:
        raise SystemExit(f'ABORT: {len(spec)} result file(s) in {rd.resolve()}, '
                         f'expected {args.expect_n}. A short table here would look '
                         f'exactly like a cohort that is genuinely smaller.')

    calls, ambig, status = [], [], []
    for sid in sorted(spec):
        genes, n_fail = parse_result(spec[sid])

        row, n = {'sample_id': sid}, {'called': 0, 'hemizygous': 0,
                                      'not_typed': 0, 'failed': 0, 'ambiguous': 0}
        for g in want:
            a = [strip_prefix(x) for x in genes.get(g, [])]
            if g not in genes:
                state, a1, a2 = 'failed', '', ''
            elif a and a[0] == NOT_TYPED:
                state, a1, a2 = 'not_typed', '', ''
            elif len(a) >= 2 and a[1] == '-':
                # KNOWN-WRONG for DRB3/4/5, deliberately kept — see
                # docs/OPEN_QUESTIONS.md §1. HLA-HD's '-' means homozygous at the
                # classical loci and hemizygous at the DRB paralogues, and filling
                # it as below is right for the first and wrong for the second. It
                # matches the reference implementation, and changing one without
                # the other would make the two incomparable.
                state, a1, a2 = 'hemizygous', a[0], a[0]
            elif len(a) > 2:
                # More than one candidate pair. Report the depth they agree on and
                # keep every candidate in the ambiguity table.
                state = 'called'
                n['ambiguous'] += 1
                a1, a2, depth = resolve_ambiguous(a)
                ambig.append({'sample_id': sid, 'gene': g,
                              'n_candidates': len(a), 'candidates': ';'.join(a),
                              'reported_1': a1, 'reported_2': a2,
                              'agreed_field_depth': depth})
            elif len(a) == 2:
                state, a1, a2 = 'called', a[0], a[1]
            else:
                state, a1, a2 = 'failed', '', ''
            n[state] += 1
            row[f'{g}_1'], row[f'{g}_2'], row[f'{g}_state'] = a1, a2, state
        calls.append(row)
        status.append({'sample_id': sid, 'genes': len(want), **n,
                       'unreadable_result_lines': n_fail})

    # The header is passed in, not taken from rows[0]. A table with no rows used to
    # be written as an empty file, and pandas cannot read one: an entirely
    # unambiguous cohort made typing_qc.py die with EmptyDataError.
    def write(path, cols, rows):
        with Path(path).open('w') as fh:
            fh.write('\t'.join(cols) + '\n')
            for r in rows:
                fh.write('\t'.join(str(r.get(k, '')) for k in cols) + '\n')

    call_cols = ['sample_id'] + [f'{g}_{s}' for g in want for s in ('1', '2', 'state')]
    write(args.out_calls, call_cols, calls)
    write(args.out_ambiguity,
          ['sample_id', 'gene', 'n_candidates', 'candidates',
           'reported_1', 'reported_2', 'agreed_field_depth'], ambig)
    write(args.out_status,
          ['sample_id', 'genes', 'called', 'hemizygous', 'not_typed', 'failed',
           'ambiguous', 'unreadable_result_lines'], status)

    tot = {k: sum(s[k] for s in status) for k in
           ('called', 'hemizygous', 'not_typed', 'failed', 'ambiguous')}
    print(f'[collect_alleles] {len(calls)} sample(s) x {len(want)} gene(s) -> {args.out_calls}')
    print(f'    called {tot["called"]}  hemizygous {tot["hemizygous"]}  '
          f'not_typed {tot["not_typed"]}  FAILED {tot["failed"]}')
    print(f'    ambiguous calls recorded: {tot["ambiguous"]} -> {args.out_ambiguity}')
    if tot['failed']:
        worst = sorted(status, key=lambda s: -s['failed'])[:5]
        print('    samples with failed genes (a pipeline failure, not biology):')
        for s in worst:
            if s['failed']:
                print(f'        {s["sample_id"]}  {s["failed"]} gene(s)')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in collect_alleles: {e}', file=sys.stderr)
        sys.exit(1)
