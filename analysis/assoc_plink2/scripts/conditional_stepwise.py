#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Stepwise conditional analysis inside one genome-wide locus.
#           Round 0 is the unconditioned fit restricted to the locus window;
#           each subsequent round adds the previous round's top variant to
#           --condition-list and re-fits. The loop stops when nothing inside the
#           window still passes the genome-wide threshold. The number of rounds
#           that yielded a signal is the number of independent signals.
# Component: assoc_plink2
# Used by : assoc_plink2.nf  process CONDITIONAL
# ---------------------------------------------------------------------------
import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

MODEL_FLAG = {'additive': '', 'dominant': 'dominant', 'recessive': 'recessive'}


def parse_args():
    p = argparse.ArgumentParser(description='Stepwise conditional analysis for one locus.')
    p.add_argument('--plink2', required=True)
    p.add_argument('--bfile', required=True)
    p.add_argument('--pheno', required=True)
    p.add_argument('--pheno-name', default='PHENO1')
    p.add_argument('--covar', required=True)
    p.add_argument('--covar-name', required=True)
    p.add_argument('--model', default='additive', choices=list(MODEL_FLAG))
    p.add_argument('--firth-mode', default='no-firth')
    p.add_argument('--cohort', required=True)
    p.add_argument('--locus-id', required=True)
    p.add_argument('--chrom', required=True)
    p.add_argument('--start', type=int, required=True)
    p.add_argument('--end', type=int, required=True)
    p.add_argument('--lead-id', required=True)
    p.add_argument('--p-threshold', type=float, default=5e-8)
    p.add_argument('--max-rounds', type=int, default=5)
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--out-dir', required=True)
    return p.parse_args()


def run_glm(args, out_prefix, condition_ids, work):
    """One plink2 GLM restricted to the locus window, optionally conditioned."""
    cmd = [args.plink2, '--bfile', args.bfile,
           '--chr', str(args.chrom).replace('chr', ''),
           '--from-bp', str(args.start), '--to-bp', str(args.end),
           '--pheno', args.pheno, '--pheno-name', args.pheno_name,
           '--covar', args.covar, '--covar-name', args.covar_name,
           '--glm', 'omit-ref', 'hide-covar']
    if MODEL_FLAG[args.model]:
        cmd.append(MODEL_FLAG[args.model])
    if args.firth_mode:
        cmd.append(args.firth_mode)
    cmd += ['--ci', '0.95', '--threads', str(args.threads), '--out', str(out_prefix)]

    if condition_ids:
        cfile = work / f'{out_prefix.name}.condlist'
        cfile.write_text('\n'.join(condition_ids) + '\n')
        cmd += ['--condition-list', str(cfile)]

    res = subprocess.run(cmd, capture_output=True, text=True)
    hits = sorted(out_prefix.parent.glob(f'{out_prefix.name}.*.glm.logistic*'))
    if not hits:
        # plink2 exits 0 with no association file when every variant in the
        # window is collinear with the conditioning set — a real, informative
        # outcome, so it must not be treated as a crash.
        print(f'[conditional] no output for {out_prefix.name}\n{res.stderr[-2000:]}', file=sys.stderr)
        return None
    df = pd.read_csv(hits[0], sep='\t', dtype={'#CHROM': str, 'ID': str, 'ERRCODE': str})
    df = df.rename(columns={'#CHROM': 'CHROM', 'LOG(OR)_SE': 'SE'})
    df['P'] = pd.to_numeric(df['P'], errors='coerce')
    return df[df['ERRCODE'].fillna('.').eq('.') & df['P'].notna()]


def main():
    args = parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    work = out / '_work'
    work.mkdir(exist_ok=True)

    signals, rounds = [], []
    condition = []
    for rnd in range(args.max_rounds + 1):
        pfx = work / f'{args.locus_id}.round{rnd}'
        df = run_glm(args, pfx, condition, work)
        if df is None or not len(df):
            rounds.append({'round': rnd, 'n_conditioned_on': len(condition),
                           'conditioned_on': ','.join(condition) or '.',
                           'n_variants': 0, 'top_id': '.', 'top_p': None,
                           'top_or': None, 'passes_threshold': False})
            break
        top = df.loc[df['P'].idxmin()]
        passes = bool(top['P'] < args.p_threshold)
        rounds.append({'round': rnd, 'n_conditioned_on': len(condition),
                       'conditioned_on': ','.join(condition) or '.',
                       'n_variants': len(df), 'top_id': top['ID'],
                       'top_p': top['P'], 'top_or': top.get('OR'),
                       'passes_threshold': passes})
        # Keep every round's statistics: the conditional plot needs them all,
        # not just the round that produced a signal.
        df.assign(round=rnd, conditioned_on=','.join(condition) or '.') \
          .to_csv(out / f'{args.locus_id}.round{rnd}.tsv', sep='\t', index=False)
        if not passes:
            break
        signals.append({'cohort': args.cohort, 'locus_id': args.locus_id, 'signal': len(signals) + 1,
                        'variant_id': top['ID'], 'chrom': top['CHROM'], 'pos': int(top['POS']),
                        'a1': top.get('A1', ''), 'a1_freq': top.get('A1_FREQ'),
                        'p': top['P'], 'or': top.get('OR'), 'se': top.get('SE'),
                        'conditioned_on': ','.join(condition) or '.'})
        condition.append(top['ID'])

    pd.DataFrame(rounds).assign(cohort=args.cohort, locus_id=args.locus_id) \
      .to_csv(out / f'{args.locus_id}.rounds.tsv', sep='\t', index=False)
    pd.DataFrame(signals, columns=['cohort', 'locus_id', 'signal', 'variant_id', 'chrom', 'pos',
                                   'a1', 'a1_freq', 'p', 'or', 'se', 'conditioned_on']) \
      .to_csv(out / f'{args.locus_id}.signals.tsv', sep='\t', index=False)

    # The scratch GLM files are ~1 file per round over a 500 kb window; small,
    # but there is no reason to publish them and every reason not to stage them.
    for f in work.glob('*.glm.logistic*'):
        f.unlink()

    print(f'[conditional] {args.cohort} {args.locus_id}: {len(signals)} independent signal(s) '
          f'over {len(rounds)} round(s); lead {args.lead_id}')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in conditional_stepwise: {e}', file=sys.stderr)
        sys.exit(1)
