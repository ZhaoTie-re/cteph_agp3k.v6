#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Decide whether an HLA-HD run actually succeeded, and exit non-zero
#           when it did not.
#
#           This stage exists because hlahd.sh cannot fail. Its estimation step is
#           a generated shell script with no `set -e`, and it checks none of its
#           children. When hla_est is killed — an out-of-memory kill on HLA-A and
#           HLA-B is the common case — hlahd.sh still exits 0, the workflow sees
#           success, and a truncated result is published. In the reference
#           implementation that happened to 6 of 348 samples, reproducibly, with
#           every exit code 0 and an empty stderr.
#
#           The truncation is not obvious downstream either: pick_up_allele writes
#           `Couldn't read result file.` on a line of its OWN, with no gene name,
#           so a parser reading by row position silently shifts every later gene.
#
#           Four checks. Any one failing is a non-zero exit, which is what lets the
#           workflow retry the sample with more memory instead of publishing it.
# Component: hla.typing
# ---------------------------------------------------------------------------
import argparse
import re
import sys
from pathlib import Path

FAIL_TOKEN = "Couldn't read result file."
KILL_RE = re.compile(r'\bKilled\b|oom-kill|Out of memory|std::bad_alloc', re.I)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample', required=True)
    p.add_argument('--result-dir', required=True, help='<sample>/result')
    p.add_argument('--log-dir', help='<sample>/log')
    p.add_argument('--run-log', help='the hlahd.sh stdout/stderr capture')
    p.add_argument('--gene-split', required=True,
                   help='HLA_gene.split.*.txt — defines which genes must be present')
    p.add_argument('--out', default=None, help='per-sample status TSV')
    return p.parse_args()


def expected_genes(path):
    genes = []
    for ln in Path(path).read_text().splitlines():
        ln = ln.strip()
        if ln and not ln.startswith('#'):
            genes.append(ln.split()[0])
    if not genes:
        raise SystemExit(f'ABORT: no gene named in {path}')
    return genes


def main():
    args = parse_args()
    rd = Path(args.result_dir)
    problems, status = [], {}

    final = rd / f'{args.sample}_final.result.txt'
    if not final.exists():
        problems.append(f'no {final.name}')
        genes_seen, n_fail_token = set(), 0
    else:
        lines = [l.rstrip('\n') for l in final.read_text().splitlines() if l.strip()]
        # Count the token BEFORE reading genes: it occupies a whole line and carries
        # no gene name, which is exactly why row-position parsing is unsafe.
        n_fail_token = sum(1 for l in lines if l.strip() == FAIL_TOKEN)
        genes_seen = {l.split('\t', 1)[0] for l in lines if l.strip() != FAIL_TOKEN}
        if n_fail_token:
            problems.append(f'{n_fail_token} line(s) of "{FAIL_TOKEN}" — '
                            f'pick_up_allele could not read an .est.txt')

    want = expected_genes(args.gene_split)
    missing = [g for g in want if g not in genes_seen]
    if missing:
        problems.append(f'{len(missing)} gene(s) absent from the result: '
                        f'{", ".join(missing[:8])}{" …" if len(missing) > 8 else ""}')

    # An .est.txt per gene is what pick_up_allele consumes. A missing one is the
    # upstream half of the same failure, and names the gene the token does not.
    missing_est = [g for g in want if not (rd / f'{args.sample}_{g}.est.txt').exists()]
    if missing_est:
        problems.append(f'{len(missing_est)} missing .est.txt: '
                        f'{", ".join(missing_est[:8])}{" …" if len(missing_est) > 8 else ""}')

    killed = []
    for p in filter(None, [args.run_log]):
        if Path(p).exists():
            killed += [ln.strip()[:120] for ln in Path(p).read_text(errors='replace').splitlines()
                       if KILL_RE.search(ln)]
    if killed:
        problems.append(f'{len(killed)} kill/OOM line(s) in the run log: {killed[0]}')

    status.update(sample_id=args.sample,
                  genes_expected=len(want), genes_present=len(want) - len(missing),
                  missing_genes=';'.join(missing) or '.',
                  missing_est=';'.join(missing_est) or '.',
                  unreadable_result_lines=n_fail_token,
                  kill_lines=len(killed),
                  status='ok' if not problems else 'FAILED')

    if args.out:
        Path(args.out).write_text(
            '\t'.join(status) + '\n' + '\t'.join(str(v) for v in status.values()) + '\n')

    if problems:
        print(f'[validate_typing] {args.sample}: FAILED', file=sys.stderr)
        for p in problems:
            print(f'    {p}', file=sys.stderr)
        print('  HLA-HD reports success regardless; this exit code is the only signal.',
              file=sys.stderr)
        sys.exit(1)

    print(f'[validate_typing] {args.sample}: ok — {len(want)} genes, '
          f'{len(want)} .est.txt, no kill lines')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in validate_typing: {e}', file=sys.stderr)
        sys.exit(1)
